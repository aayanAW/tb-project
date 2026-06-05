"""
phase6_environmental/posterior_propagation.py — Posterior-propagation robustness analysis.

Propagates the MCMC posterior from Phase 5 (omega, dG_spacer) through the
Gillespie simulation to verify that CV(asymmetric) > CV(symmetric) holds
robustly across the posterior distribution.

For each posterior sample (omega_i, dG_spacer_i):
  1. Use the partition function to compute equilibrium state probabilities
     at the estimated Mce3R concentration.
  2. Derive effective Kd values that incorporate cooperativity:
       - The doubly-bound state weight is omega * exp(-beta*dG_spacer) times
         what it would be without cooperativity.
       - This modifies the effective dissociation when the other site is
         already occupied: Kd_eff = Kd_base / (omega * exp(-dG_spacer/kT))
       - We approximate this in the Gillespie engine by adjusting k_off_strong
         and k_off_weak by the cooperativity factor.
  3. Run short Gillespie simulations for asymmetric and symmetric architectures.
  4. Compute CV for each, record delta_CV = CV_asym - CV_sym.

Output:  results/phase6/posterior_propagation.csv
"""

import sys
import os
import time
import numpy as np
import csv

sys.path.insert(0, '/Users/aayanalwani/tb project/mce3r_stochastic')
from config.parameters import PARAMS, kT_KCAL
from phase2_simulation.operator_model import OperatorModel
from phase2_simulation.gillespie_engine import run_population


# ---- Configuration ----
N_POSTERIOR_SAMPLES = 50
N_CELLS = 2000
MASTER_SEED = 77777
OUTPUT_DIR = '/Users/aayanalwani/tb project/mce3r_stochastic/results/phase6'
POSTERIOR_PATH = '/Users/aayanalwani/tb project/mce3r_stochastic/results/phase5/mcmc_posteriors.npz'


def load_posterior_samples(path, n_samples, seed=42):
    """
    Load and thin MCMC posterior samples for omega and dG_spacer.

    Parameters
    ----------
    path : str
        Path to mcmc_posteriors.npz.
    n_samples : int
        Number of thinned samples to return.
    seed : int
        Random seed for reproducible thinning.

    Returns
    -------
    omega_samples : np.ndarray shape (n_samples,)
    dG_spacer_samples : np.ndarray shape (n_samples,)
    """
    data = np.load(path)
    omega_all = data['omega_samples']
    dG_spacer_all = data['dG_spacer_samples']
    n_total = len(omega_all)

    # Thin by selecting evenly-spaced indices with jitter for decorrelation
    rng = np.random.default_rng(seed)
    step = max(1, n_total // (n_samples * 10))  # large thinning factor
    thinned_omega = omega_all[::step]
    thinned_dG = dG_spacer_all[::step]

    # Randomly select n_samples from the thinned set
    n_avail = len(thinned_omega)
    if n_avail < n_samples:
        # If not enough after thinning, just subsample from the full chain
        idx = rng.choice(n_total, size=n_samples, replace=False)
        return omega_all[idx], dG_spacer_all[idx]

    idx = rng.choice(n_avail, size=n_samples, replace=False)
    return thinned_omega[idx], thinned_dG[idx]


def compute_cv(proteins):
    """Compute coefficient of variation from protein counts."""
    p = np.asarray(proteins, dtype=np.float64)
    mean = np.mean(p)
    std = np.std(p)
    if mean <= 0:
        return np.nan
    return std / mean


def run_single_condition(Kd_strong, Kd_weak, block_strong, block_weak,
                         n_cells, seed):
    """
    Run Gillespie for one architecture with given parameters.

    Returns CV of final protein distribution.
    """
    model = OperatorModel(
        Kd_strong=Kd_strong,
        Kd_weak=Kd_weak,
        block_strong=block_strong,
        block_weak=block_weak,
    )
    arrays = model.get_numba_arrays()

    result = run_population(
        arrays, n_cells=n_cells, master_seed=seed,
        record_traces=0, n_trace_cells=0
    )
    return compute_cv(result['proteins'])


def cooperativity_factor(omega, dG_spacer, kT=kT_KCAL):
    """
    Compute the cooperativity enhancement factor for doubly-bound state.

    When one site is already occupied, the effective Kd for the second site
    is modified by: Kd_eff = Kd_base / coop_factor

    This means k_off_eff = k_off_base / coop_factor (lower off-rate = tighter binding).

    Parameters
    ----------
    omega : float
        Cooperativity factor (>1 = positive cooperativity).
    dG_spacer : float
        Spacer interaction energy (kcal/mol, negative = favorable).
    kT : float
        Thermal energy (kcal/mol).

    Returns
    -------
    coop_factor : float
        Multiplicative factor for the second binding event.
        coop_factor > 1 means cooperative (tighter binding).
    """
    return omega * np.exp(-dG_spacer / kT)


def effective_Kd(Kd_base, omega, dG_spacer, kT=kT_KCAL):
    """
    Compute effective Kd when the other site is already occupied.

    The cooperativity factor reduces the effective Kd:
      Kd_eff = Kd_base / (omega * exp(-dG_spacer/kT))

    Clip to a minimum of 0.01 nM to avoid numerical issues.
    """
    cf = cooperativity_factor(omega, dG_spacer, kT)
    kd_eff = Kd_base / cf
    return max(kd_eff, 0.01)


def run_posterior_propagation():
    """
    Main function: propagate MCMC posterior through Gillespie simulations.

    For each posterior sample:
      - Asymmetric architecture: Kd_strong=2.4, Kd_weak=49.0
        (cooperativity modifies effective binding in doubly-occupied state)
      - Symmetric architecture: Kd = geometric mean = 10.84 for both sites

    The cooperativity factor from (omega, dG_spacer) modifies the effective
    dissociation constant for the second binding event. We approximate this
    by adjusting the Kd values passed to the OperatorModel.

    Specifically, for a given (omega, dG_spacer):
      - The effective Kd for the second site to bind (when one is already
        occupied) is reduced by the cooperativity factor.
      - We model this as a weighted-average effective Kd:
          Kd_eff_strong = (Kd_strong + effective_Kd(Kd_strong, omega, dG_spacer)) / 2
          Kd_eff_weak   = (Kd_weak   + effective_Kd(Kd_weak,   omega, dG_spacer)) / 2
        This averages the "first binding" (no cooperativity) and "second binding"
        (with cooperativity) Kd values to get an overall effective Kd.

    Note: Since both architectures use the same (omega, dG_spacer), the
    cooperativity affects both equally. The noise difference arises from
    the asymmetry in Kd_strong vs Kd_weak.
    """
    print("=" * 70)
    print("POSTERIOR PROPAGATION ROBUSTNESS ANALYSIS")
    print("=" * 70)

    # Load posterior samples
    print(f"\nLoading {N_POSTERIOR_SAMPLES} posterior samples from {POSTERIOR_PATH}...")
    omega_samples, dG_spacer_samples = load_posterior_samples(
        POSTERIOR_PATH, N_POSTERIOR_SAMPLES, seed=42
    )
    print(f"  omega:     median={np.median(omega_samples):.3f}, "
          f"range=[{omega_samples.min():.3f}, {omega_samples.max():.3f}]")
    print(f"  dG_spacer: median={np.median(dG_spacer_samples):.3f}, "
          f"range=[{dG_spacer_samples.min():.3f}, {dG_spacer_samples.max():.3f}]")

    # Base parameters
    Kd_strong_base = PARAMS['Kd_strong']       # 2.4 nM
    Kd_weak_base = PARAMS['Kd_weak']           # 49.0 nM
    Kd_sym_base = PARAMS['Kd_symmetric']       # ~10.84 nM
    block_strong = PARAMS['block_strong']       # 0.85
    block_weak = PARAMS['block_weak']           # 0.50
    block_sym = PARAMS['block_symmetric']       # ~0.652

    # Results storage
    results = []
    t_start = time.time()

    print(f"\nRunning {N_POSTERIOR_SAMPLES} posterior samples "
          f"x 2 architectures x {N_CELLS} cells each...")
    print(f"{'Sample':>6}  {'omega':>8}  {'dG_sp':>8}  {'coop_f':>8}  "
          f"{'CV_asym':>8}  {'CV_sym':>8}  {'dCV':>8}  {'time':>6}")
    print("-" * 74)

    for i in range(N_POSTERIOR_SAMPLES):
        t_sample = time.time()
        omega_i = float(omega_samples[i])
        dG_spacer_i = float(dG_spacer_samples[i])

        # Cooperativity factor
        cf = cooperativity_factor(omega_i, dG_spacer_i)

        # Effective Kd values incorporating cooperativity
        # For the asymmetric architecture:
        #   When strong site is occupied, effective Kd for weak site binding
        #   is reduced by cooperativity. And vice versa.
        #   We use the average of "first binding Kd" and "cooperative Kd"
        #   to represent the overall effective binding strength.
        Kd_strong_eff = (Kd_strong_base + effective_Kd(Kd_strong_base, omega_i, dG_spacer_i)) / 2.0
        Kd_weak_eff = (Kd_weak_base + effective_Kd(Kd_weak_base, omega_i, dG_spacer_i)) / 2.0

        # For the symmetric architecture (same logic, but both sites identical):
        Kd_sym_eff = (Kd_sym_base + effective_Kd(Kd_sym_base, omega_i, dG_spacer_i)) / 2.0

        # Seed for this sample (different for asym vs sym to avoid correlation)
        seed_asym = MASTER_SEED + i * 2
        seed_sym = MASTER_SEED + i * 2 + 1

        # Run asymmetric
        cv_asym = run_single_condition(
            Kd_strong=Kd_strong_eff,
            Kd_weak=Kd_weak_eff,
            block_strong=block_strong,
            block_weak=block_weak,
            n_cells=N_CELLS,
            seed=seed_asym,
        )

        # Run symmetric
        cv_sym = run_single_condition(
            Kd_strong=Kd_sym_eff,
            Kd_weak=Kd_sym_eff,
            block_strong=block_sym,
            block_weak=block_sym,
            n_cells=N_CELLS,
            seed=seed_sym,
        )

        delta_cv = cv_asym - cv_sym
        elapsed = time.time() - t_sample

        results.append({
            'sample_idx': i,
            'omega': omega_i,
            'dG_spacer': dG_spacer_i,
            'coop_factor': cf,
            'Kd_strong_eff': Kd_strong_eff,
            'Kd_weak_eff': Kd_weak_eff,
            'Kd_sym_eff': Kd_sym_eff,
            'cv_asym': cv_asym,
            'cv_sym': cv_sym,
            'delta_cv': delta_cv,
        })

        print(f"{i+1:>6}  {omega_i:>8.3f}  {dG_spacer_i:>8.3f}  {cf:>8.3f}  "
              f"{cv_asym:>8.4f}  {cv_sym:>8.4f}  {delta_cv:>+8.4f}  {elapsed:>5.1f}s")

    total_time = time.time() - t_start
    print(f"\nTotal runtime: {total_time:.1f}s ({total_time/60:.1f} min)")

    # ---- Save results ----
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    csv_path = os.path.join(OUTPUT_DIR, 'posterior_propagation.csv')
    fieldnames = ['sample_idx', 'omega', 'dG_spacer', 'coop_factor',
                  'Kd_strong_eff', 'Kd_weak_eff', 'Kd_sym_eff',
                  'cv_asym', 'cv_sym', 'delta_cv']
    with open(csv_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(results)
    print(f"\nResults saved to {csv_path}")

    # ---- Summary statistics ----
    cv_asym_arr = np.array([r['cv_asym'] for r in results])
    cv_sym_arr = np.array([r['cv_sym'] for r in results])
    delta_cv_arr = np.array([r['delta_cv'] for r in results])

    n_valid = np.sum(np.isfinite(delta_cv_arr))
    n_asym_wins = np.sum(delta_cv_arr > 0)
    frac_asym_wins = n_asym_wins / n_valid if n_valid > 0 else 0.0

    mean_delta = np.nanmean(delta_cv_arr)
    ci_lo = np.nanpercentile(delta_cv_arr, 2.5)
    ci_hi = np.nanpercentile(delta_cv_arr, 97.5)

    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"  Posterior samples analyzed:        {N_POSTERIOR_SAMPLES}")
    print(f"  Valid (non-NaN) results:           {n_valid}")
    print(f"  Cells per condition:               {N_CELLS}")
    print(f"")
    print(f"  CV(asymmetric):  mean={np.nanmean(cv_asym_arr):.4f}, "
          f"std={np.nanstd(cv_asym_arr):.4f}")
    print(f"  CV(symmetric):   mean={np.nanmean(cv_sym_arr):.4f}, "
          f"std={np.nanstd(cv_sym_arr):.4f}")
    print(f"")
    print(f"  delta_CV = CV(asym) - CV(sym):")
    print(f"    Mean:           {mean_delta:+.4f}")
    print(f"    95% CI:         [{ci_lo:+.4f}, {ci_hi:+.4f}]")
    print(f"    Median:         {np.nanmedian(delta_cv_arr):+.4f}")
    print(f"")
    print(f"  Fraction where CV(asym) > CV(sym): {frac_asym_wins:.1%} "
          f"({n_asym_wins}/{n_valid})")
    print(f"")

    robust = frac_asym_wins > 0.95
    if robust:
        print(f"  RESULT: ROBUST — >{95}% of posterior samples show CV(asym) > CV(sym)")
    else:
        print(f"  RESULT: NOT ROBUST at 95% threshold "
              f"(only {frac_asym_wins:.1%} of samples)")

    return results, {
        'n_samples': N_POSTERIOR_SAMPLES,
        'n_valid': int(n_valid),
        'n_cells': N_CELLS,
        'frac_asym_wins': float(frac_asym_wins),
        'mean_delta_cv': float(mean_delta),
        'ci_lo': float(ci_lo),
        'ci_hi': float(ci_hi),
        'robust': robust,
    }


# ---------------------------------------------------------------------------
# Self-tests
# ---------------------------------------------------------------------------
if __name__ == '__main__':
    results, summary = run_posterior_propagation()

    print("\n" + "=" * 70)
    print("SELF-TESTS")
    print("=" * 70)

    n_pass = 0
    n_fail = 0

    def sanity(name, cond, msg=""):
        global n_pass, n_fail
        if cond:
            print(f"  SANITY PASS: {name} {msg}")
            n_pass += 1
        else:
            print(f"  SANITY FAIL: {name} {msg}")
            n_fail += 1

    # Test 1: At least 50 rows in output
    sanity("At least 50 rows",
           len(results) >= 50,
           f"-- got {len(results)} rows")

    # Test 2: >80% of samples show CV_asym > CV_sym
    sanity(">80% samples show CV(asym) > CV(sym)",
           summary['frac_asym_wins'] > 0.80,
           f"-- got {summary['frac_asym_wins']:.1%}")

    # Test 3: No NaN values in delta_cv
    n_nan = sum(1 for r in results if np.isnan(r['delta_cv']))
    sanity("No NaN values in delta_cv",
           n_nan == 0,
           f"-- found {n_nan} NaN values")

    # Test 4: All omega values positive
    all_omega_pos = all(r['omega'] > 0 for r in results)
    sanity("All omega values positive",
           all_omega_pos)

    # Test 5: CSV file exists and has correct number of rows
    csv_path = os.path.join(OUTPUT_DIR, 'posterior_propagation.csv')
    if os.path.exists(csv_path):
        with open(csv_path) as f:
            reader = csv.reader(f)
            n_rows = sum(1 for _ in reader) - 1  # minus header
        sanity("CSV has correct row count",
               n_rows == len(results),
               f"-- CSV has {n_rows} rows, expected {len(results)}")
    else:
        sanity("CSV file exists", False, "-- file not found")

    # Test 6: CV values are in reasonable range (0-2)
    all_reasonable = all(
        0 < r['cv_asym'] < 2 and 0 < r['cv_sym'] < 2
        for r in results
    )
    sanity("All CV values in reasonable range (0, 2)",
           all_reasonable)

    print(f"\n{'=' * 50}")
    print(f"Sanity: {n_pass} pass, {n_fail} fail")
    if n_fail > 0:
        print("STOPPING: Sanity check failed.")
        sys.exit(1)
    else:
        print("All checks passed.")
        sys.exit(0)
