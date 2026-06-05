"""
phase3_analysis/statistical_tests.py — Statistical tests on Phase 2 simulation data.

KS test: A vs B, A vs C, A vs D, B vs D, C vs D
GMM LR test: LL_total = model.score(X) * N, LR = 2*(LL3_total - LL2_total)
Cohen's d with correct pooled SD formula
Output: results/phase3/statistical_tests.csv
"""

import sys
import os
import warnings
import numpy as np
import pandas as pd
from scipy import stats

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


def cohens_d(x1, x2):
    """
    Cohen's d with pooled SD.
    pooled_sd = sqrt(((n1-1)*var1 + (n2-1)*var2) / (n1+n2-2))
    Guard divide-by-zero.
    """
    n1, n2 = len(x1), len(x2)
    var1 = np.var(x1, ddof=1) if n1 > 1 else 0.0
    var2 = np.var(x2, ddof=1) if n2 > 1 else 0.0

    denom_df = n1 + n2 - 2
    if denom_df <= 0:
        return np.nan

    pooled_var = ((n1 - 1) * var1 + (n2 - 1) * var2) / denom_df
    pooled_sd = np.sqrt(pooled_var)

    if pooled_sd == 0:
        return np.nan

    return (np.mean(x1) - np.mean(x2)) / pooled_sd


def gmm_lr_test(x, k_null=2, k_alt=3):
    """
    GMM likelihood ratio test.
    LL_total = model.score(X) * N (TOTAL log-likelihood, not average).
    LR = 2*(LL_alt_total - LL_null_total)
    Returns (LR_stat, p_value_approx)
    """
    from sklearn.mixture import GaussianMixture

    X = x.reshape(-1, 1)
    N = len(x)

    # Add dither
    std_x = np.std(X)
    if std_x > 0:
        X_dither = X + np.random.default_rng(42).normal(0, std_x * 1e-4, size=X.shape)
    else:
        X_dither = X.copy()

    def fit_and_score(n_comp):
        try:
            gmm = GaussianMixture(
                n_components=n_comp,
                reg_covar=PARAMS.get('gmm_reg_covar', 1e-3),
                n_init=PARAMS.get('gmm_n_init', 5),
                random_state=PARAMS.get('gmm_random_state', 42),
                max_iter=300,
            )
            gmm.fit(X_dither)
            # score() returns average log-likelihood; multiply by N for total
            ll_total = gmm.score(X_dither) * N
            return ll_total
        except Exception:
            return -np.inf

    ll_null = fit_and_score(k_null)
    ll_alt = fit_and_score(k_alt)

    lr_stat = 2.0 * (ll_alt - ll_null)
    if lr_stat < 0:
        lr_stat = 0.0  # alt should fit at least as well

    # Approximate p-value using chi2 (heuristic — not exact for GMM boundary)
    # df = difference in free parameters: 3-comp has 3 more params than 2-comp
    # (1 weight + 1 mean + 1 variance = 3 params per component)
    df_diff = 3  # approximate
    p_value = 1.0 - stats.chi2.cdf(lr_stat, df=df_diff) if lr_stat > 0 else 1.0

    return lr_stat, p_value


def run_statistical_tests():
    """Run all statistical tests and save CSV."""
    os.makedirs(PHASE3_DIR, exist_ok=True)

    # Load all conditions
    data = {}
    for cond in ['A', 'B', 'C', 'D']:
        data[cond] = load_proteins(cond)

    rows = []

    # --- KS tests ---
    pairs = [('A', 'B'), ('A', 'C'), ('A', 'D'), ('B', 'D'), ('C', 'D')]
    for c1, c2 in pairs:
        ks_stat, ks_p = stats.ks_2samp(data[c1], data[c2])
        d = cohens_d(data[c1], data[c2])
        rows.append({
            'test': 'KS',
            'comparison': f'{c1}_vs_{c2}',
            'statistic': ks_stat,
            'p_value': ks_p,
            'cohens_d': d,
        })
        print(f"  KS {c1} vs {c2}: stat={ks_stat:.4f}, p={ks_p:.4g}, d={d:.4f}")

    # --- GMM LR tests (for each condition individually) ---
    for cond in ['A', 'B', 'C', 'D']:
        lr_stat, lr_p = gmm_lr_test(data[cond], k_null=2, k_alt=3)
        rows.append({
            'test': 'GMM_LR_2v3',
            'comparison': f'{cond}',
            'statistic': lr_stat,
            'p_value': lr_p,
            'cohens_d': np.nan,
        })
        print(f"  GMM LR (2v3) {cond}: LR={lr_stat:.4f}, p={lr_p:.4g}")

    # --- GMM LR 1 vs 2 ---
    for cond in ['A', 'B', 'C', 'D']:
        lr_stat, lr_p = gmm_lr_test(data[cond], k_null=1, k_alt=2)
        rows.append({
            'test': 'GMM_LR_1v2',
            'comparison': f'{cond}',
            'statistic': lr_stat,
            'p_value': lr_p,
            'cohens_d': np.nan,
        })
        print(f"  GMM LR (1v2) {cond}: LR={lr_stat:.4f}, p={lr_p:.4g}")

    df = pd.DataFrame(rows)
    out_path = os.path.join(PHASE3_DIR, 'statistical_tests.csv')
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

    print("=== statistical_tests.py self-tests ===")

    df = run_statistical_tests()

    # SANITY: all p-values in [0, 1]
    p_vals = df['p_value'].values
    sanity("p-values in [0,1]",
           np.all((p_vals >= 0) & (p_vals <= 1)),
           f"— range: [{np.min(p_vals):.4g}, {np.max(p_vals):.4g}]")

    # SANITY: all statistics >= 0
    stats_vals = df['statistic'].values
    sanity("statistics >= 0", np.all(stats_vals >= 0),
           f"— min: {np.min(stats_vals):.4g}")

    # SANITY: KS statistics in [0, 1]
    ks_rows = df[df['test'] == 'KS']
    ks_stats = ks_rows['statistic'].values
    sanity("KS stats in [0,1]",
           np.all((ks_stats >= 0) & (ks_stats <= 1)),
           f"— range: [{np.min(ks_stats):.4f}, {np.max(ks_stats):.4f}]")

    # SANITY: Cohen's d is finite for KS comparisons
    ks_d = ks_rows['cohens_d'].values
    sanity("Cohen's d finite", np.all(np.isfinite(ks_d)),
           f"— range: [{np.min(ks_d):.4f}, {np.max(ks_d):.4f}]")

    # SANITY: correct number of rows (5 KS + 4 LR_2v3 + 4 LR_1v2 = 13)
    sanity("correct row count", len(df) == 13, f"— got {len(df)}")

    # SANITY: output file exists
    out_path = os.path.join(PHASE3_DIR, 'statistical_tests.csv')
    sanity("output file exists", os.path.exists(out_path))

    # SCIENTIFIC: A vs D should show significant difference (large d)
    d_ad = ks_rows[ks_rows['comparison'] == 'A_vs_D']['cohens_d'].values[0]
    scientific("A vs D large effect", abs(d_ad) > 0.5,
               f"— Cohen's d={d_ad:.4f}")

    print(f"\n{'='*50}")
    print(f"SANITY: {n_sanity_pass} pass, {n_sanity_fail} fail")
    print(f"SCIENTIFIC: {n_science_pass} pass, {n_science_warn} warn")
    if n_sanity_fail > 0:
        print("STOPPING: Sanity check failed.")
        sys.exit(1)
    else:
        print("All checks passed.")
        sys.exit(0)
