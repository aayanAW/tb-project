"""
phase3_analysis/noise_metrics.py — Compute noise metrics and GMM fits for Phase 2 data.

For each condition (A, B, C, D): mean, variance, CV, Fano factor, bimodality coefficient.
Fit 1,2,3-component GMMs with robustness safeguards.
Select best model by BIC.
Output: results/phase3/noise_metrics.csv
"""

import sys
import os
import warnings
import numpy as np
import pandas as pd
from scipy.stats import skew, kurtosis

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


def bimodality_coefficient(x):
    """
    Bimodality coefficient: BC = (skewness^2 + 1) / (kurtosis_excess + 3).
    BC > 5/9 (~0.555) suggests bimodality.
    """
    n = len(x)
    if n < 4:
        return np.nan
    s = skew(x)
    k = kurtosis(x, fisher=True)  # excess kurtosis
    denom = k + 3.0
    if denom == 0:
        return np.nan
    return (s ** 2 + 1.0) / denom


def fit_gmm_robust(X, n_components, random_state=None):
    """
    Fit a GMM with robustness safeguards:
    - Add dither noise
    - reg_covar=1e-3
    - n_init=5
    - try/except with BIC=inf fallback
    """
    from sklearn.mixture import GaussianMixture

    # Dither noise to prevent singular covariance
    X_dither = X.copy()
    std_x = np.std(X_dither)
    if std_x > 0:
        X_dither = X_dither + np.random.default_rng(random_state or 42).normal(
            0, std_x * 1e-4, size=X_dither.shape
        )

    try:
        gmm = GaussianMixture(
            n_components=n_components,
            reg_covar=PARAMS.get('gmm_reg_covar', 1e-3),
            n_init=PARAMS.get('gmm_n_init', 5),
            random_state=random_state or PARAMS.get('gmm_random_state', 42),
            max_iter=300,
        )
        gmm.fit(X_dither)
        bic = gmm.bic(X_dither)
        return gmm, bic
    except Exception as e:
        warnings.warn(f"GMM fit failed for n_components={n_components}: {e}")
        return None, np.inf


def compute_noise_metrics(proteins, label):
    """Compute noise metrics for a single condition."""
    n = len(proteins)
    mean_val = np.mean(proteins)
    var_val = np.var(proteins, ddof=1) if n > 1 else 0.0
    cv = np.sqrt(var_val) / mean_val if mean_val > 0 else np.nan
    fano = var_val / mean_val if mean_val > 0 else np.nan
    bc = bimodality_coefficient(proteins)

    # GMM fitting
    X = proteins.reshape(-1, 1)
    best_bic = np.inf
    best_k = 1
    best_gmm = None
    for k in PARAMS.get('n_gmm_components', [1, 2, 3]):
        gmm, bic = fit_gmm_robust(X, k, random_state=PARAMS.get('gmm_random_state', 42))
        if bic < best_bic:
            best_bic = bic
            best_k = k
            best_gmm = gmm

    return {
        'condition': label,
        'n_cells': n,
        'mean': mean_val,
        'variance': var_val,
        'CV': cv,
        'Fano': fano,
        'bimodality_coeff': bc,
        'gmm_best_k': best_k,
        'gmm_best_bic': best_bic,
    }


def run_noise_metrics():
    """Compute noise metrics for all conditions and save CSV."""
    os.makedirs(PHASE3_DIR, exist_ok=True)
    rows = []
    for cond in ['A', 'B', 'C', 'D']:
        proteins = load_proteins(cond)
        row = compute_noise_metrics(proteins, cond)
        rows.append(row)
    df = pd.DataFrame(rows)
    out_path = os.path.join(PHASE3_DIR, 'noise_metrics.csv')
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

    print("=== noise_metrics.py self-tests ===")

    df = run_noise_metrics()

    # SANITY: all values finite
    for col in ['mean', 'variance', 'CV', 'Fano']:
        vals = df[col].values
        sanity(f"{col} all finite", np.all(np.isfinite(vals)),
               f"— {col}: {vals}")

    # SANITY: correct number of conditions
    sanity("4 conditions", len(df) == 4, f"— got {len(df)}")

    # SANITY: means are positive
    sanity("means positive", np.all(df['mean'].values > 0),
           f"— means: {df['mean'].values}")

    # SANITY: gmm_best_k in [1,2,3]
    sanity("gmm_best_k valid", all(k in [1, 2, 3] for k in df['gmm_best_k']),
           f"— best_k: {df['gmm_best_k'].values}")

    # SANITY: output file exists
    out_path = os.path.join(PHASE3_DIR, 'noise_metrics.csv')
    sanity("output file exists", os.path.exists(out_path))

    # SCIENTIFIC: CV(A) > CV(B) — asymmetric operator should produce more noise
    cv_a = df.loc[df['condition'] == 'A', 'CV'].values[0]
    cv_b = df.loc[df['condition'] == 'B', 'CV'].values[0]
    scientific("CV(A) > CV(B)", cv_a > cv_b,
               f"— CV_A={cv_a:.4f}, CV_B={cv_b:.4f}")

    # SCIENTIFIC: Fano(D) close to expected unregulated
    fano_d = df.loc[df['condition'] == 'D', 'Fano'].values[0]
    expected_fano = PARAMS['expected_fano_unregulated']
    scientific("Fano(D) ~ expected_unregulated",
               expected_fano * 0.3 < fano_d < expected_fano * 3.0,
               f"— Fano_D={fano_d:.2f}, expected~{expected_fano:.2f}")

    print(f"\n{'='*50}")
    print(f"SANITY: {n_sanity_pass} pass, {n_sanity_fail} fail")
    print(f"SCIENTIFIC: {n_science_pass} pass, {n_science_warn} warn")
    if n_sanity_fail > 0:
        print("STOPPING: Sanity check failed.")
        sys.exit(1)
    else:
        print("All checks passed.")
        sys.exit(0)
