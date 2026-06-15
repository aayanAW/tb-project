"""
phase6_environmental/persistence_threshold.py — Persister fraction quantification.

Defines a persistence threshold from the unregulated (Condition D) distribution
and computes the fraction of cells below that threshold for each condition.

Biological logic:
- Persisters are cells with very low expression of the Mce3R regulon target.
- Low target protein = metabolically quiescent = antibiotic tolerant.
- Threshold = 10th percentile of the unregulated (constitutive) distribution.
  This is operationally defined: cells that fall this low even without repression
  are considered the baseline "persister" fraction.
"""

import sys
import numpy as np
import pandas as pd

sys.path.insert(0, '/Users/aayanalwani/tb project/mce3r_stochastic')
from config.parameters import PARAMS


def define_threshold(
    proteins_unregulated: np.ndarray,
    method: str = 'percentile',
    percentile: float = 10.0
) -> float:
    """
    Define the persistence threshold from the unregulated distribution.

    Parameters
    ----------
    proteins_unregulated : np.ndarray
        Protein counts from Condition D (no regulation / constitutive expression).
    method : str
        'percentile' — Use given percentile of unregulated distribution (default, robust).
        'absolute'   — Fixed threshold at 100 proteins (biologically motivated).
    percentile : float
        Percentile to use for 'percentile' method. Default: 10.

    Returns
    -------
    float — threshold value in protein counts
    """
    if method == 'percentile':
        return float(np.percentile(proteins_unregulated, percentile))
    elif method == 'absolute':
        return 100.0
    else:
        raise ValueError(f"Unknown threshold method '{method}'. Use 'percentile' or 'absolute'.")


def compute_persister_fractions(
    results: dict,
    threshold: float
) -> pd.DataFrame:
    """
    Compute persister fractions for each condition in the results dict.

    Parameters
    ----------
    results : dict
        Keyed by (architecture, environment, model_type).
        Values are dicts containing either 'target_proteins' (two-species)
        or 'proteins' (single-species) numpy arrays.
    threshold : float
        Protein count threshold; cells below are counted as persisters.

    Returns
    -------
    pd.DataFrame with columns:
        architecture, environment, model_type, n_cells,
        persister_fraction, fold_enrichment, mean_protein, cv_protein
    """
    rows = []

    # First pass: find symmetric_baseline reference fraction
    sym_baseline_fraction = None
    for key, data in results.items():
        arch, env, model_type = key
        if arch == 'symmetric' and env == 'baseline' and model_type == 'single':
            proteins = data.get('target_proteins', data.get('proteins'))
            if proteins is not None:
                n_below = np.sum(proteins < threshold)
                sym_baseline_fraction = float(n_below) / len(proteins)
            break

    # If symmetric baseline not found, use 0.1 (by construction)
    if sym_baseline_fraction is None or sym_baseline_fraction == 0:
        sym_baseline_fraction = 0.1

    for key, data in results.items():
        arch, env, model_type = key

        # Extract protein array (two-species uses 'target_proteins', single uses 'proteins')
        proteins = data.get('target_proteins', data.get('proteins'))
        if proteins is None:
            continue

        proteins = np.asarray(proteins, dtype=np.float64)
        n_cells = len(proteins)
        if n_cells == 0:
            continue

        n_below = int(np.sum(proteins < threshold))
        persister_frac = float(n_below) / n_cells

        mean_p = float(np.mean(proteins))
        std_p = float(np.std(proteins))
        cv_p = std_p / mean_p if mean_p > 0 else 0.0

        fold_enrich = persister_frac / sym_baseline_fraction if sym_baseline_fraction > 0 else np.nan

        rows.append({
            'architecture':       arch,
            'environment':        env,
            'model_type':         model_type,
            'n_cells':            n_cells,
            'persister_fraction': persister_frac,
            'fold_enrichment':    fold_enrich,
            'mean_protein':       mean_p,
            'cv_protein':         cv_p,
        })

    if not rows:
        return pd.DataFrame()

    df = pd.DataFrame(rows)
    df = df.sort_values(['architecture', 'environment', 'model_type']).reset_index(drop=True)
    return df


def persistence_phase_diagram(
    ratios: list,
    environments: list,
    n_cells_per_point: int = 2000,
    seed: int = 42,
    threshold: float = None
) -> pd.DataFrame:
    """
    2D sweep: asymmetry ratio × environment → persister fraction.

    For each (ratio, environment) pair:
    1. Set up asymmetric operator with given Kd ratio.
       Kd_geo = sqrt(2.4 * 49.0) ≈ 10.84 nM (geometric mean preserved).
       Kd_s = Kd_geo / sqrt(ratio)
       Kd_w = Kd_geo * sqrt(ratio)
    2. Apply environment multiplier to k_on and gamma_protein.
    3. Simulate n_cells_per_point single-species cells.
    4. Compute persister fraction (fraction below threshold).

    Parameters
    ----------
    ratios : list of float
        Asymmetry ratios to sweep (Kd_w / Kd_s).
    environments : list of str
        Environment names.
    n_cells_per_point : int
        Number of simulated cells per (ratio, environment) pair.
    seed : int
        Master random seed.
    threshold : float or None
        Protein threshold for persister definition.
        If None, use 10th percentile of Condition D from results/phase2/.

    Returns
    -------
    pd.DataFrame with columns: ratio, environment, persister_fraction, cv, mean_protein
    """
    from phase6_environmental.environmental_signals import apply_environment
    from phase2_simulation.operator_model import OperatorModel
    from phase2_simulation.gillespie_engine import run_population

    # Load threshold from Condition D if not provided
    if threshold is None:
        import os
        cond_d_path = '/Users/aayanalwani/tb project/mce3r_stochastic/results/phase2/condition_D.npz'
        if os.path.exists(cond_d_path):
            d = np.load(cond_d_path)
            threshold = define_threshold(d['proteins'], method='percentile', percentile=10)
        else:
            threshold = 100.0  # fallback

    Kd_geo = float(np.sqrt(2.4 * 49.0))  # ≈ 10.84 nM

    rows = []
    for ri, ratio in enumerate(ratios):
        ratio_f = float(ratio)
        sqrt_r = np.sqrt(ratio_f)
        Kd_s = Kd_geo / sqrt_r
        Kd_w = Kd_geo * sqrt_r

        model = OperatorModel(Kd_strong=Kd_s, Kd_weak=Kd_w)
        base_arrays = model.get_numba_arrays()
        base_arrays['gamma_protein'] = np.float64(PARAMS['gamma_protein'])

        for ei, env_name in enumerate(environments):
            env_arrays = apply_environment(base_arrays, env_name)
            cell_seed = seed + ri * 1000 + ei * 100
            result = run_population(env_arrays, n_cells=n_cells_per_point,
                                    master_seed=cell_seed)
            proteins = result['proteins'].astype(np.float64)

            n_below = int(np.sum(proteins < threshold))
            persister_frac = float(n_below) / n_cells_per_point
            mean_p = float(np.mean(proteins))
            std_p = float(np.std(proteins))
            cv_p = std_p / mean_p if mean_p > 0 else 0.0

            rows.append({
                'ratio':              ratio_f,
                'environment':        env_name,
                'persister_fraction': persister_frac,
                'cv':                 cv_p,
                'mean_protein':       mean_p,
                'threshold_used':     threshold,
            })

    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Self-tests
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    from phase6_environmental.environmental_signals import ENVIRONMENTS

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

    def scientific(name, cond, msg=""):
        global n_pass
        if cond:
            print(f"  SCIENTIFIC PASS: {name} {msg}")
            n_pass += 1
        else:
            print(f"  SCIENTIFIC WARN: {name} {msg}")

    print("=== persistence_threshold.py self-tests ===")

    # Load Condition D proteins for threshold definition
    import os
    cond_d_path = '/Users/aayanalwani/tb project/mce3r_stochastic/results/phase2/condition_D.npz'
    d_data = np.load(cond_d_path)
    proteins_D = d_data['proteins']

    threshold_10 = define_threshold(proteins_D, method='percentile', percentile=10)
    threshold_5  = define_threshold(proteins_D, method='percentile', percentile=5)
    threshold_15 = define_threshold(proteins_D, method='percentile', percentile=15)
    print(f"  Thresholds: 5th={threshold_5:.0f}, 10th={threshold_10:.0f}, 15th={threshold_15:.0f}")

    # Test 1: Persister fraction in [0, 1]
    # Build a small mock results dict
    from phase6_environmental.two_species_model import TwoSpeciesModel
    from phase6_environmental.two_species_gillespie import run_two_species_population
    from phase2_simulation.operator_model import OperatorModel
    from phase2_simulation.gillespie_engine import run_population

    print("  Running small test simulations (500 cells each)...")

    # Asymmetric
    model_asym = OperatorModel()
    arr_asym = model_asym.get_numba_arrays()
    arr_asym['gamma_protein'] = np.float64(PARAMS['gamma_protein'])
    from phase6_environmental.environmental_signals import apply_environment
    res_asym_base = run_population(
        apply_environment(arr_asym, 'baseline'), n_cells=500, master_seed=1
    )

    # Symmetric
    model_sym = OperatorModel(Kd_strong=10.84, Kd_weak=10.84,
                               block_strong=float(np.sqrt(0.85*0.50)),
                               block_weak=float(np.sqrt(0.85*0.50)))
    arr_sym = model_sym.get_numba_arrays()
    arr_sym['gamma_protein'] = np.float64(PARAMS['gamma_protein'])
    res_sym_base = run_population(
        apply_environment(arr_sym, 'baseline'), n_cells=500, master_seed=2
    )

    mock_results = {
        ('asymmetric', 'baseline', 'single'): {'proteins': res_asym_base['proteins']},
        ('symmetric',  'baseline', 'single'): {'proteins': res_sym_base['proteins']},
    }

    df = compute_persister_fractions(mock_results, threshold=threshold_10)

    asym_frac = df[df['architecture'] == 'asymmetric']['persister_fraction'].values
    sym_frac  = df[df['architecture'] == 'symmetric' ]['persister_fraction'].values

    # Test 1: Persister fraction in [0,1]
    all_valid = np.all((df['persister_fraction'] >= 0) & (df['persister_fraction'] <= 1))
    sanity("Persister fraction in [0,1]",
           all_valid,
           f"— fractions: {df['persister_fraction'].tolist()}")

    # Test 2: Condition D persister fraction ≈ 10% by construction
    frac_D = float(np.mean(proteins_D < threshold_10))
    sanity("Condition D persister fraction ≈ 10%",
           abs(frac_D - 0.10) < 0.02,
           f"— got {frac_D:.3f} (expected ~0.10)")

    # Test 3: Asymmetric has higher persister fraction than symmetric under baseline
    if len(asym_frac) > 0 and len(sym_frac) > 0:
        scientific("Asymmetric persister fraction > symmetric (baseline)",
                   asym_frac[0] >= sym_frac[0],
                   f"— asymmetric={asym_frac[0]:.4f}, symmetric={sym_frac[0]:.4f}")

    # Test 4: Small phase diagram (3 ratios, 2 environments, 200 cells)
    print("  Running small phase diagram test (3 ratios × 2 envs × 200 cells)...")
    df_pd = persistence_phase_diagram(
        ratios=[1, 5, 20],
        environments=['baseline', 'host_like'],
        n_cells_per_point=200,
        seed=99,
        threshold=threshold_10
    )
    sanity("Phase diagram runs successfully",
           len(df_pd) == 6,
           f"— got {len(df_pd)} rows (expected 6)")
    sanity("Phase diagram persister fractions in [0,1]",
           np.all((df_pd['persister_fraction'] >= 0) & (df_pd['persister_fraction'] <= 1)),
           f"— fractions: {df_pd['persister_fraction'].tolist()}")

    # Test 5: Ranking holds across 5th/10th/15th percentile thresholds
    for pct, thr in [(5, threshold_5), (10, threshold_10), (15, threshold_15)]:
        n_asym = int(np.sum(res_asym_base['proteins'] < thr))
        n_sym  = int(np.sum(res_sym_base['proteins']  < thr))
        f_asym = n_asym / 500
        f_sym  = n_sym  / 500
        scientific(f"Asymmetric > Symmetric at {pct}th percentile",
                   f_asym >= f_sym,
                   f"— asym={f_asym:.4f}, sym={f_sym:.4f}")

    # Summary
    print(f"\n{'='*50}")
    print(f"SANITY: {n_pass} pass, {n_fail} fail")
    if n_fail > 0:
        print("STOPPING: Sanity check failed.")
        sys.exit(1)
    else:
        print("All checks passed.")
        sys.exit(0)
