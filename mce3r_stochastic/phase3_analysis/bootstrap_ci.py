"""
phase3_analysis/bootstrap_ci.py — Bootstrap confidence intervals for noise metrics.

Computes bootstrap CIs for: CV, Fano, intermediate fraction, persister fraction.
n_bootstrap=1000, ci=0.95
Output: results/phase3/bootstrap_results.csv
"""

import sys
import os
import numpy as np
import pandas as pd

sys.path.insert(0, '/Users/aayanalwani/tb project/mce3r_stochastic')
from config.parameters import PARAMS

PROJECT_ROOT = '/Users/aayanalwani/tb project/mce3r_stochastic'
PHASE2_DIR = os.path.join(PROJECT_ROOT, 'results', 'phase2')
PHASE3_DIR = os.path.join(PROJECT_ROOT, 'results', 'phase3')


def load_proteins(condition):
    """Load protein array from Phase 2 .npz file."""
    path = os.path.join(PHASE2_DIR, f'condition_{condition}.npz')
    data = np.load(path)
    return data['proteins'].astype(np.float64)


def compute_cv(x):
    """Coefficient of variation."""
    m = np.mean(x)
    if m == 0:
        return np.nan
    return np.std(x, ddof=1) / m


def compute_fano(x):
    """Fano factor."""
    m = np.mean(x)
    if m == 0:
        return np.nan
    return np.var(x, ddof=1) / m


def compute_intermediate_fraction(x, low_q=0.25, high_q=0.75):
    """Fraction of cells in intermediate expression range (between 25th and 75th percentile of full range)."""
    if len(x) == 0:
        return np.nan
    xmin, xmax = np.min(x), np.max(x)
    rng = xmax - xmin
    if rng == 0:
        return 1.0
    low_thresh = xmin + low_q * rng
    high_thresh = xmin + high_q * rng
    return np.mean((x >= low_thresh) & (x <= high_thresh))


def compute_persister_fraction(x, proteins_B=None):
    """
    Persister fraction: cells below mean_B + 2*std_B threshold (secondary method).
    If proteins_B is None, use 10th percentile as fallback.
    """
    if proteins_B is not None and len(proteins_B) > 1:
        threshold = np.mean(proteins_B) + 2.0 * np.std(proteins_B, ddof=1)
    else:
        threshold = np.percentile(x, 10)
    return np.mean(x < threshold)


def bootstrap_statistic(x, stat_func, n_bootstrap=1000, ci=0.95, seed=None,
                         extra_arg=None):
    """
    Compute bootstrap CI for a statistic.

    Returns: (point_estimate, ci_lower, ci_upper)
    """
    rng = np.random.default_rng(seed)
    n = len(x)
    point = stat_func(x) if extra_arg is None else stat_func(x, extra_arg)

    boot_vals = np.zeros(n_bootstrap)
    for i in range(n_bootstrap):
        idx = rng.integers(0, n, size=n)
        sample = x[idx]
        if extra_arg is not None:
            boot_vals[i] = stat_func(sample, extra_arg)
        else:
            boot_vals[i] = stat_func(sample)

    alpha = 1.0 - ci
    lo = np.nanpercentile(boot_vals, 100 * alpha / 2)
    hi = np.nanpercentile(boot_vals, 100 * (1 - alpha / 2))
    return point, lo, hi


def run_bootstrap_ci(n_bootstrap=None, ci=None):
    """Compute bootstrap CIs for all conditions and metrics."""
    if n_bootstrap is None:
        n_bootstrap = PARAMS.get('n_bootstrap', 1000)
    if ci is None:
        ci = PARAMS.get('ci_level', 0.95)
    seed = PARAMS.get('bootstrap_seed', 123)

    os.makedirs(PHASE3_DIR, exist_ok=True)

    # Load condition B for persister threshold
    proteins_B = load_proteins('B')

    rows = []
    for cond in ['A', 'B', 'C', 'D']:
        proteins = load_proteins(cond)

        # CV
        pt, lo, hi = bootstrap_statistic(proteins, compute_cv,
                                          n_bootstrap=n_bootstrap, ci=ci, seed=seed)
        rows.append({'condition': cond, 'metric': 'CV',
                     'point': pt, 'CI_lower': lo, 'CI_upper': hi})

        # Fano
        pt, lo, hi = bootstrap_statistic(proteins, compute_fano,
                                          n_bootstrap=n_bootstrap, ci=ci, seed=seed + 1)
        rows.append({'condition': cond, 'metric': 'Fano',
                     'point': pt, 'CI_lower': lo, 'CI_upper': hi})

        # Intermediate fraction
        pt, lo, hi = bootstrap_statistic(proteins, compute_intermediate_fraction,
                                          n_bootstrap=n_bootstrap, ci=ci, seed=seed + 2)
        rows.append({'condition': cond, 'metric': 'intermediate_fraction',
                     'point': pt, 'CI_lower': lo, 'CI_upper': hi})

        # Persister fraction (using condition B threshold)
        def persister_func(x, ref=proteins_B):
            return compute_persister_fraction(x, ref)

        pt, lo, hi = bootstrap_statistic(proteins, persister_func,
                                          n_bootstrap=n_bootstrap, ci=ci, seed=seed + 3)
        rows.append({'condition': cond, 'metric': 'persister_fraction',
                     'point': pt, 'CI_lower': lo, 'CI_upper': hi})

    df = pd.DataFrame(rows)
    out_path = os.path.join(PHASE3_DIR, 'bootstrap_results.csv')
    df.to_csv(out_path, index=False)
    print(f"  Wrote {out_path}")
    return df


# ---------------------------------------------------------------------------
# Self-tests
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    n_sanity_pass = 0
    n_sanity_fail = 0
    n_science_pass = 0
    n_science_warn = 0

    def sanity(name, cond, msg=""):
        global n_sanity_pass, n_sanity_fail
        if cond:
            print(f"  SANITY PASS: {name} {msg}")
            n_sanity_pass += 1
        else:
            print(f"  SANITY FAIL: {name} {msg}")
            n_sanity_fail += 1

    def scientific(name, cond, msg=""):
        global n_science_pass, n_science_warn
        if cond:
            print(f"  SCIENTIFIC PASS: {name} {msg}")
            n_science_pass += 1
        else:
            print(f"  SCIENTIFIC WARN: {name} {msg}")
            n_science_warn += 1

    print("=== bootstrap_ci.py self-tests ===")

    df = run_bootstrap_ci()

    # SANITY: CI_lower <= point <= CI_upper for all rows
    all_ci_valid = True
    for _, row in df.iterrows():
        lo, pt, hi = row['CI_lower'], row['point'], row['CI_upper']
        if not (lo <= pt <= hi):
            # Bootstrap CIs can rarely not bracket point with small N
            # but check for gross errors
            if not np.isnan(lo) and not np.isnan(hi):
                all_ci_valid = False
                print(f"    CI violation: {row['condition']} {row['metric']}: "
                      f"{lo:.4f} <= {pt:.4f} <= {hi:.4f}")
    sanity("CI_lower <= point <= CI_upper", all_ci_valid)

    # SANITY: all values finite
    for col in ['point', 'CI_lower', 'CI_upper']:
        sanity(f"{col} all finite", np.all(np.isfinite(df[col].values)),
               f"— {col}")

    # SANITY: correct number of rows (4 conditions x 4 metrics = 16)
    sanity("16 rows", len(df) == 16, f"— got {len(df)}")

    # SANITY: output file exists
    out_path = os.path.join(PHASE3_DIR, 'bootstrap_results.csv')
    sanity("output file exists", os.path.exists(out_path))

    # SANITY: CI width > 0 (non-degenerate)
    ci_widths = df['CI_upper'].values - df['CI_lower'].values
    sanity("CI widths >= 0", np.all(ci_widths >= 0),
           f"— min_width={np.min(ci_widths):.4f}")

    print(f"\n{'='*50}")
    print(f"SANITY: {n_sanity_pass} pass, {n_sanity_fail} fail")
    print(f"SCIENTIFIC: {n_science_pass} pass, {n_science_warn} warn")
    if n_sanity_fail > 0:
        print("STOPPING: Sanity check failed.")
        sys.exit(1)
    else:
        print("All checks passed.")
        sys.exit(0)
