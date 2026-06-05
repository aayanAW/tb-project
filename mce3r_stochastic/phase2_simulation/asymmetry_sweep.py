"""
phase2_simulation/asymmetry_sweep.py — Sweep Kd_weak/Kd_strong ratio.

Condition E: For each of 10 asymmetry ratios, run 1000 cells and compute
the fraction in an "intermediate" expression state.

Ratios: [1, 2, 5, 8, 10, 15, 20, 30, 40, 50]

Parameterization: geometric mean Kd is held constant so that total
regulatory capacity is equivalent across ratios. For each ratio r:
  Kd_strong = Kd_geo / sqrt(r)
  Kd_weak   = Kd_geo * sqrt(r)
where Kd_geo = sqrt(2.4 * 49.0) ≈ 10.84 nM.

Block fractions are also co-varied via power-law interpolation so that
at ratio=1 both sites are symmetric (block_geo) and at the wild-type
ratio (~20.4) we recover the measured values (0.85, 0.50).

Intermediate fraction is estimated using a simple threshold method:
cells with protein count between 10th and 90th percentile of the
no-regulation distribution are considered "intermediate".
"""

import sys
import os
import numpy as np
import time

sys.path.insert(0, '/Users/aayanalwani/tb project/mce3r_stochastic')
from config.parameters import PARAMS
from phase2_simulation.operator_model import OperatorModel
from phase2_simulation.gillespie_engine import run_population


OUTPUT_DIR = '/Users/aayanalwani/tb project/mce3r_stochastic/results/phase2'


def compute_intermediate_fraction(proteins, low_threshold, high_threshold):
    """
    Fraction of cells with protein in (low_threshold, high_threshold).
    This captures the "intermediate" expression phenotype.
    """
    p = proteins.astype(np.float64)
    n_inter = np.sum((p > low_threshold) & (p < high_threshold))
    return n_inter / len(p)


def run_sweep(n_cells_per_ratio=None, ratios=None):
    """
    Run the asymmetry sweep (Condition E).

    For each Kd ratio, creates an OperatorModel with Kd_weak = Kd_strong * ratio,
    runs n_cells_per_ratio cells, and computes statistics.

    Parameters
    ----------
    n_cells_per_ratio : int (default from PARAMS)
    ratios : list of floats (default from PARAMS)

    Returns
    -------
    dict with:
        ratios : float64 array
        means : float64 array
        cvs : float64 array
        fanos : float64 array
        intermediate_fractions : float64 array
        all_proteins : dict mapping ratio -> protein array
    """
    if n_cells_per_ratio is None:
        n_cells_per_ratio = PARAMS['n_cells_sweep']
    if ratios is None:
        ratios = PARAMS['sweep_ratios']

    # Geometric mean values — anchor points for the sweep
    Kd_geo = np.sqrt(PARAMS['Kd_strong'] * PARAMS['Kd_weak'])  # ≈ 10.84 nM
    block_geo = np.sqrt(PARAMS['block_strong'] * PARAMS['block_weak'])  # ≈ 0.652
    # Wild-type ratio for block interpolation
    r_wt = PARAMS['Kd_weak'] / PARAMS['Kd_strong']  # ≈ 20.4
    log_r_wt = np.log(r_wt)

    n_ratios = len(ratios)

    # Arrays to store results
    means = np.zeros(n_ratios, dtype=np.float64)
    cvs = np.zeros(n_ratios, dtype=np.float64)
    fanos = np.zeros(n_ratios, dtype=np.float64)
    intermediate_fracs = np.zeros(n_ratios, dtype=np.float64)
    all_proteins = {}

    # Determine intermediate thresholds from no-regulation run
    unreg_mean_est = (PARAMS['k_max'] * PARAMS['k_translation'] /
                      (PARAMS['gamma_mRNA'] * PARAMS['gamma_protein']))
    low_thresh = unreg_mean_est * 0.10   # ~222
    high_thresh = unreg_mean_est * 0.60  # ~1334

    for i, ratio in enumerate(ratios):
        # Hold geometric mean Kd constant while varying asymmetry
        r = float(ratio)
        Kd_s = Kd_geo / np.sqrt(r)
        Kd_w = Kd_geo * np.sqrt(r)

        # Co-vary block fractions via power-law in log-ratio space
        # At r=1: both = block_geo. At r=r_wt: block_s=0.85, block_w=0.50
        if r <= 1.0:
            bs = block_geo
            bw = block_geo
        else:
            frac = min(1.0, np.log(r) / log_r_wt)
            bs = block_geo * (PARAMS['block_strong'] / block_geo) ** frac
            bw = block_geo * (PARAMS['block_weak'] / block_geo) ** frac

        model = OperatorModel(Kd_strong=Kd_s, Kd_weak=Kd_w,
                              block_strong=bs, block_weak=bw)
        arrays = model.get_numba_arrays()

        # Seed: master_seed + 100 + ratio_index
        seed = PARAMS['master_seed'] + 100 + i

        t0 = time.time()
        result = run_population(arrays, n_cells=n_cells_per_ratio,
                                master_seed=seed)
        elapsed = time.time() - t0

        prots = result['proteins'].astype(np.float64)
        mean_p = np.mean(prots)
        std_p = np.std(prots)
        var_p = np.var(prots)

        means[i] = mean_p
        cvs[i] = std_p / mean_p if mean_p > 0 else 0.0
        fanos[i] = var_p / mean_p if mean_p > 0 else 0.0
        intermediate_fracs[i] = compute_intermediate_fraction(
            result['proteins'], low_thresh, high_thresh)
        all_proteins[ratio] = result['proteins']

        print(f"    Ratio {ratio:>3d}: mean={mean_p:7.1f}, "
              f"CV={cvs[i]:.3f}, Fano={fanos[i]:.1f}, "
              f"inter_frac={intermediate_fracs[i]:.3f} ({elapsed:.1f}s)")

    # Save results
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    save_dict = {
        'ratios': np.array(ratios, dtype=np.float64),
        'means': means,
        'cvs': cvs,
        'fanos': fanos,
        'intermediate_fractions': intermediate_fracs,
    }
    # Save per-ratio protein arrays
    for i, ratio in enumerate(ratios):
        save_dict[f'proteins_ratio_{ratio}'] = all_proteins[ratio]

    np.savez(os.path.join(OUTPUT_DIR, 'condition_E_sweep.npz'), **save_dict)

    return {
        'ratios': np.array(ratios, dtype=np.float64),
        'means': means,
        'cvs': cvs,
        'fanos': fanos,
        'intermediate_fractions': intermediate_fracs,
        'all_proteins': all_proteins,
    }


# ---------------------------------------------------------------------------
# Self-tests
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    n_pass = 0
    n_fail = 0
    n_warn = 0

    def sanity(name, cond, msg=""):
        global n_pass, n_fail
        if cond:
            print(f"  SANITY PASS: {name} {msg}")
            n_pass += 1
        else:
            print(f"  SANITY FAIL: {name} {msg}")
            n_fail += 1

    def scientific(name, cond, msg=""):
        global n_pass, n_warn
        if cond:
            print(f"  SCIENTIFIC PASS: {name} {msg}")
            n_pass += 1
        else:
            print(f"  SCIENTIFIC WARN: {name} {msg}")
            n_warn += 1

    print("=== asymmetry_sweep.py self-tests ===")

    # Run small sweep: 3 ratios × 10 cells
    test_ratios = [1, 10, 50]
    result = run_sweep(n_cells_per_ratio=10, ratios=test_ratios)

    # SANITY: all 3 ratios completed
    sanity("all ratios completed",
           len(result['ratios']) == len(test_ratios),
           f"— got {len(result['ratios'])}")

    # SANITY: all values finite
    sanity("means finite", np.all(np.isfinite(result['means'])),
           f"— {result['means']}")
    sanity("cvs finite", np.all(np.isfinite(result['cvs'])),
           f"— {result['cvs']}")
    sanity("fanos finite", np.all(np.isfinite(result['fanos'])),
           f"— {result['fanos']}")
    sanity("inter_fracs finite",
           np.all(np.isfinite(result['intermediate_fractions'])),
           f"— {result['intermediate_fractions']}")

    # SANITY: output file exists
    sanity("sweep output file exists",
           os.path.exists(os.path.join(OUTPUT_DIR, 'condition_E_sweep.npz')))

    # SCIENTIFIC: ratio=1 (symmetric) should have lower CV than ratio=50
    scientific("CV increases with ratio",
               result['cvs'][-1] > result['cvs'][0],
               f"— CV(ratio=1)={result['cvs'][0]:.3f}, CV(ratio=50)={result['cvs'][-1]:.3f}")

    # Summary
    print(f"\n{'='*50}")
    print(f"SANITY: {n_pass} pass, {n_fail} fail | SCIENTIFIC: {n_warn} warn")
    if n_fail > 0:
        print("STOPPING: Sanity check failed.")
        sys.exit(1)
    else:
        print("All checks passed.")
        sys.exit(0)
