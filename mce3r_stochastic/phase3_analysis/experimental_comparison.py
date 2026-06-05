"""
phase3_analysis/experimental_comparison.py — Compare simulation to published data.

Persister threshold: GMM intersection (primary) and 2-sigma (secondary).
Predicted fold-change: mean_protein_D / mean_protein_A.
Compare to published values from Santangelo 2009.
Output: appended columns to noise_metrics.csv and logged.
"""

import sys
import os
import warnings
import numpy as np
import pandas as pd

sys.path.insert(0, '/Users/aayanalwani/tb project/mce3r_stochastic')
from config.parameters import PARAMS

PROJECT_ROOT = '/Users/aayanalwani/tb project/mce3r_stochastic'
PHASE2_DIR = os.path.join(PROJECT_ROOT, 'results', 'phase2')
PHASE3_DIR = os.path.join(PROJECT_ROOT, 'results', 'phase3')
EXTERNAL_DIR = os.path.join(PROJECT_ROOT, 'data', 'external')


def load_proteins(condition):
    """Load protein array from Phase 2 .npz file."""
    path = os.path.join(PHASE2_DIR, f'condition_{condition}.npz')
    data = np.load(path)
    return data['proteins'].astype(np.float64)


def gmm_intersection_threshold(proteins, n_components=2):
    """
    Find intersection point between two GMM components (primary persister threshold).
    Returns threshold value, or None if GMM fails or has < 2 components.
    """
    from sklearn.mixture import GaussianMixture

    X = proteins.reshape(-1, 1)
    std_x = np.std(X)
    if std_x > 0:
        X_dither = X + np.random.default_rng(42).normal(0, std_x * 1e-4, size=X.shape)
    else:
        X_dither = X.copy()

    try:
        gmm = GaussianMixture(
            n_components=n_components,
            reg_covar=PARAMS.get('gmm_reg_covar', 1e-3),
            n_init=PARAMS.get('gmm_n_init', 5),
            random_state=PARAMS.get('gmm_random_state', 42),
            max_iter=300,
        )
        gmm.fit(X_dither)
    except Exception as e:
        warnings.warn(f"GMM intersection: fit failed: {e}")
        return None

    if n_components < 2:
        return None

    # Sort components by mean
    means = gmm.means_.flatten()
    sorted_idx = np.argsort(means)
    mu1 = means[sorted_idx[0]]
    mu2 = means[sorted_idx[1]]

    # Find intersection by evaluating log-probabilities on a grid
    x_grid = np.linspace(mu1, mu2, 1000).reshape(-1, 1)
    log_probs = gmm.predict_proba(x_grid)

    # Find where probability of lower component crosses higher component
    prob_low = log_probs[:, sorted_idx[0]]
    prob_high = log_probs[:, sorted_idx[1]]
    diff = prob_low - prob_high

    # Find zero crossing
    sign_changes = np.where(np.diff(np.sign(diff)))[0]
    if len(sign_changes) == 0:
        # No crossing — use midpoint
        return (mu1 + mu2) / 2.0

    # Use first crossing
    idx = sign_changes[0]
    threshold = x_grid[idx, 0]
    return threshold


def two_sigma_threshold(proteins_B):
    """
    Secondary persister threshold: mean_B + 2*std_B.
    """
    return np.mean(proteins_B) + 2.0 * np.std(proteins_B, ddof=1)


def run_experimental_comparison():
    """Compare simulation predictions to published data."""
    os.makedirs(PHASE3_DIR, exist_ok=True)

    proteins_A = load_proteins('A')
    proteins_B = load_proteins('B')
    proteins_D = load_proteins('D')

    results = {}

    # --- Persister fraction ---
    # PRIMARY METHOD: GMM component weights.
    # Fit 2-component GMM to Condition A. The HIGH-expression component represents
    # cells with de-repressed mce3 operon (persisters). Persister fraction = weight
    # of the high component. This avoids threshold artifacts when components overlap.
    from sklearn.mixture import GaussianMixture

    X_A = proteins_A.reshape(-1, 1)
    std_A = np.std(X_A)
    X_dither = X_A + np.random.default_rng(42).normal(0, std_A * 1e-4, size=X_A.shape) if std_A > 0 else X_A.copy()

    try:
        gmm = GaussianMixture(
            n_components=2,
            reg_covar=PARAMS.get('gmm_reg_covar', 1e-3),
            n_init=PARAMS.get('gmm_n_init', 5),
            random_state=PARAMS.get('gmm_random_state', 42),
            max_iter=300,
        )
        gmm.fit(X_dither)
        means = gmm.means_.flatten()
        weights = gmm.weights_.flatten()
        sorted_idx = np.argsort(means)  # low component first
        high_idx = sorted_idx[-1]       # highest mean = de-repressed = persister
        persister_frac_gmm = weights[high_idx]
        persister_mean_high = means[high_idx]
        persister_mean_low = means[sorted_idx[0]]
        gmm_thresh = gmm_intersection_threshold(proteins_A, n_components=2)
    except Exception as e:
        warnings.warn(f"GMM persister fraction: {e}")
        persister_frac_gmm = np.nan
        persister_mean_high = np.nan
        persister_mean_low = np.nan
        gmm_thresh = None

    # SECONDARY METHOD: 2-sigma threshold from Condition B
    sigma2_thresh = two_sigma_threshold(proteins_B)
    persister_frac_2sigma = np.mean(proteins_A > sigma2_thresh)

    # Use GMM weight as primary
    persister_frac_A = persister_frac_gmm
    threshold = gmm_thresh  # still report intersection for reference

    results['gmm_threshold'] = gmm_thresh
    results['sigma2_threshold'] = sigma2_thresh
    results['persister_threshold'] = threshold
    results['persister_fraction_A'] = persister_frac_A
    results['persister_fraction_2sigma'] = persister_frac_2sigma
    results['gmm_component_low_mean'] = persister_mean_low
    results['gmm_component_high_mean'] = persister_mean_high

    print(f"  Persister fraction (GMM weight of high component): {persister_frac_A:.4f}")
    print(f"  GMM low mean: {persister_mean_low:.1f}, high mean: {persister_mean_high:.1f}")
    print(f"  Persister fraction (2-sigma, cells > {sigma2_thresh:.1f}): {persister_frac_2sigma:.4f}")
    print(f"  GMM intersection threshold (reference): {gmm_thresh}")

    # --- Predicted fold-change ---
    mean_A = np.mean(proteins_A)
    mean_D = np.mean(proteins_D)
    fold_change = mean_D / mean_A if mean_A > 0 else np.nan
    results['mean_A'] = mean_A
    results['mean_D'] = mean_D
    results['fold_change_D_over_A'] = fold_change
    print(f"  Predicted fold-change (D/A): {fold_change:.2f}")

    # --- Compare to published ---
    ext_path = os.path.join(EXTERNAL_DIR, 'santangelo2009_foldchange.csv')
    if os.path.exists(ext_path):
        pub_df = pd.read_csv(ext_path)
        yrbE3A_row = pub_df[pub_df['gene'] == 'yrbE3A']
        if len(yrbE3A_row) > 0:
            pub_fc = yrbE3A_row['fold_change'].values[0]
            results['published_fold_change'] = pub_fc
            results['fold_change_ratio'] = fold_change / pub_fc if pub_fc > 0 else np.nan
            print(f"  Published fold-change (yrbE3A): {pub_fc}")
            print(f"  Prediction/Published ratio: {results['fold_change_ratio']:.2f}")
    else:
        results['published_fold_change'] = np.nan
        results['fold_change_ratio'] = np.nan

    # --- Append to noise_metrics.csv if it exists ---
    nm_path = os.path.join(PHASE3_DIR, 'noise_metrics.csv')
    if os.path.exists(nm_path):
        nm_df = pd.read_csv(nm_path)
        nm_df['persister_threshold'] = np.nan
        nm_df['persister_fraction'] = np.nan
        nm_df['fold_change_vs_D'] = np.nan

        # For condition A: use GMM weight. For others: use 2-sigma threshold.
        for _, row in nm_df.iterrows():
            cond = row['condition']
            prots = load_proteins(cond)
            if cond == 'A':
                nm_df.loc[nm_df['condition'] == cond, 'persister_fraction'] = persister_frac_A
            elif sigma2_thresh is not None:
                nm_df.loc[nm_df['condition'] == cond, 'persister_fraction'] = np.mean(prots > sigma2_thresh)
            nm_df.loc[nm_df['condition'] == cond, 'fold_change_vs_D'] = mean_D / np.mean(prots) if np.mean(prots) > 0 else np.nan

        nm_df.loc[nm_df['condition'] == 'A', 'persister_threshold'] = threshold
        nm_df.to_csv(nm_path, index=False)
        print(f"  Updated {nm_path} with persister/fold-change columns")

    # --- Save comparison log ---
    log_df = pd.DataFrame([results])
    log_path = os.path.join(PHASE3_DIR, 'experimental_comparison.csv')
    log_df.to_csv(log_path, index=False)
    print(f"  Wrote {log_path}")

    return results


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

    print("=== experimental_comparison.py self-tests ===")

    results = run_experimental_comparison()

    # SANITY: fold_change is finite and positive
    fc = results.get('fold_change_D_over_A', np.nan)
    sanity("fold_change finite and positive",
           np.isfinite(fc) and fc > 0,
           f"— fold_change={fc:.4f}")

    # SANITY: means are positive
    sanity("mean_A > 0", results.get('mean_A', 0) > 0,
           f"— mean_A={results.get('mean_A', 0):.2f}")
    sanity("mean_D > 0", results.get('mean_D', 0) > 0,
           f"— mean_D={results.get('mean_D', 0):.2f}")

    # SANITY: threshold is finite
    thresh = results.get('persister_threshold')
    sanity("threshold finite", thresh is not None and np.isfinite(thresh),
           f"— threshold={thresh}")

    # SANITY: output files exist
    sanity("comparison csv exists",
           os.path.exists(os.path.join(PHASE3_DIR, 'experimental_comparison.csv')))

    # SCIENTIFIC: persister fraction in [0.01%, 10%] = [0.0001, 0.10]
    pf = results.get('persister_fraction_A', np.nan)
    scientific("persister fraction in [0.01%, 10%]",
               0.0001 <= pf <= 0.10,
               f"— persister_fraction={pf:.4f} ({pf*100:.2f}%)")

    # SCIENTIFIC: fold-change within order of magnitude of published
    ratio = results.get('fold_change_ratio', np.nan)
    if np.isfinite(ratio):
        scientific("fold-change within 10x of published",
                   0.1 < ratio < 10.0,
                   f"— prediction/published={ratio:.2f}")
    else:
        scientific("fold-change comparison", False, "— no published data")

    print(f"\n{'='*50}")
    print(f"SANITY: {n_sanity_pass} pass, {n_sanity_fail} fail")
    print(f"SCIENTIFIC: {n_science_pass} pass, {n_science_warn} warn")
    if n_sanity_fail > 0:
        print("STOPPING: Sanity check failed.")
        sys.exit(1)
    else:
        print("All checks passed.")
        sys.exit(0)
