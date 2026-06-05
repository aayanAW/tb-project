"""
phase3_analysis/sensitivity_analysis.py — Parameter sensitivity analysis.

Vary 7 parameters +/-50% in 10 steps, 1000 cells each (production).
Self-test uses 10 cells, 3 steps for speed.
Computes CV and n_modes at each step.
Output: results/phase3/sensitivity_data.csv
"""

import sys
import os
import warnings
import numpy as np
import pandas as pd

sys.path.insert(0, '/Users/aayanalwani/tb project/mce3r_stochastic')
from config.parameters import PARAMS

PROJECT_ROOT = '/Users/aayanalwani/tb project/mce3r_stochastic'
PHASE3_DIR = os.path.join(PROJECT_ROOT, 'results', 'phase3')

# Parameters to vary
SENSITIVITY_PARAMS = ['k_max', 'k_translation', 't_half_mRNA', 't_half_protein',
                       'block_strong', 'block_weak', 'k_on']


def _run_single_point(param_name, param_value, n_cells, seed):
    """
    Run Gillespie simulation with one parameter varied.
    Returns protein array.
    """
    from phase2_simulation.operator_model import OperatorModel
    from phase2_simulation.gillespie_engine import run_population

    # Build model overrides based on which parameter is varied
    model_kwargs = {}
    extra_params = {}

    if param_name == 'k_max':
        model_kwargs['k_max'] = param_value
    elif param_name == 'k_translation':
        extra_params['k_translation'] = param_value
    elif param_name == 't_half_mRNA':
        extra_params['gamma_mRNA'] = np.log(2) / param_value
    elif param_name == 't_half_protein':
        extra_params['gamma_protein'] = np.log(2) / param_value
    elif param_name == 'block_strong':
        model_kwargs['block_strong'] = param_value
    elif param_name == 'block_weak':
        model_kwargs['block_weak'] = param_value
    elif param_name == 'k_on':
        model_kwargs['k_on'] = param_value

    model = OperatorModel(**model_kwargs)
    arrays = model.get_numba_arrays()

    # For parameters not handled by OperatorModel, we need to monkey-patch
    # the run_population call. The engine reads k_translation, gamma_mRNA,
    # gamma_protein from PARAMS. We temporarily override them.
    old_vals = {}
    for k, v in extra_params.items():
        old_vals[k] = PARAMS[k]
        PARAMS[k] = v

    try:
        result = run_population(arrays, n_cells=n_cells, master_seed=seed,
                                record_traces=0, n_trace_cells=0)
    finally:
        # Restore original PARAMS
        for k, v in old_vals.items():
            PARAMS[k] = v

    return result['proteins'].astype(np.float64)


def count_gmm_modes(x, max_k=3):
    """Count number of GMM modes by BIC selection."""
    from sklearn.mixture import GaussianMixture

    X = x.reshape(-1, 1)
    std_x = np.std(X)
    if std_x > 0:
        X_dither = X + np.random.default_rng(42).normal(0, std_x * 1e-4, size=X.shape)
    else:
        X_dither = X.copy()

    best_bic = np.inf
    best_k = 1
    for k in range(1, max_k + 1):
        try:
            gmm = GaussianMixture(
                n_components=k,
                reg_covar=1e-3,
                n_init=5,
                random_state=42,
                max_iter=300,
            )
            gmm.fit(X_dither)
            bic = gmm.bic(X_dither)
            if bic < best_bic:
                best_bic = bic
                best_k = k
        except Exception:
            pass
    return best_k


def run_sensitivity(n_steps=None, n_cells=None):
    """
    Run sensitivity analysis for all 7 parameters.
    Returns DataFrame.
    """
    if n_steps is None:
        n_steps = PARAMS.get('sensitivity_steps', 10)
    if n_cells is None:
        n_cells = PARAMS.get('sensitivity_cells', 1000)

    sensitivity_range = PARAMS.get('sensitivity_range', 0.50)
    os.makedirs(PHASE3_DIR, exist_ok=True)

    rows = []
    seed_counter = 5000

    for param_name in SENSITIVITY_PARAMS:
        base_val = PARAMS[param_name]
        lo = base_val * (1.0 - sensitivity_range)
        hi = base_val * (1.0 + sensitivity_range)

        # For block fractions, clamp to [0.01, 0.99]
        if param_name.startswith('block_'):
            lo = max(0.01, lo)
            hi = min(0.99, hi)
        # For k_on, keep positive
        if param_name == 'k_on':
            lo = max(1e-6, lo)

        steps = np.linspace(lo, hi, n_steps)

        for step_val in steps:
            proteins = _run_single_point(param_name, step_val, n_cells, seed_counter)
            seed_counter += n_cells + 1

            mean_val = np.mean(proteins)
            var_val = np.var(proteins, ddof=1) if len(proteins) > 1 else 0.0
            cv = np.sqrt(var_val) / mean_val if mean_val > 0 else np.nan
            n_modes = count_gmm_modes(proteins)

            rows.append({
                'parameter': param_name,
                'base_value': base_val,
                'test_value': step_val,
                'fold_change': step_val / base_val if base_val != 0 else np.nan,
                'n_cells': n_cells,
                'mean': mean_val,
                'CV': cv,
                'n_modes': n_modes,
            })
            print(f"    {param_name}={step_val:.4g}: CV={cv:.4f}, modes={n_modes}")

    df = pd.DataFrame(rows)
    out_path = os.path.join(PHASE3_DIR, 'sensitivity_data.csv')
    df.to_csv(out_path, index=False)
    print(f"  Wrote {out_path}")
    return df


# ---------------------------------------------------------------------------
# Self-tests (small: 10 cells, 3 steps for speed)
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

    print("=== sensitivity_analysis.py self-tests ===")
    print("  Using small test: n_cells=10, n_steps=3")

    df = run_sensitivity(n_steps=3, n_cells=10)

    # SANITY: no NaN in CV
    sanity("no NaN in CV", not df['CV'].isna().any(),
           f"— NaN count: {df['CV'].isna().sum()}")

    # SANITY: correct dimensions (7 params x 3 steps = 21 rows)
    expected_rows = len(SENSITIVITY_PARAMS) * 3
    sanity("correct row count", len(df) == expected_rows,
           f"— expected {expected_rows}, got {len(df)}")

    # SANITY: all n_modes in [1,2,3]
    sanity("n_modes valid", all(m in [1, 2, 3] for m in df['n_modes']),
           f"— modes: {df['n_modes'].unique()}")

    # SANITY: all means > 0
    sanity("means positive", np.all(df['mean'].values > 0),
           f"— min mean: {df['mean'].min():.2f}")

    # SANITY: output file exists
    out_path = os.path.join(PHASE3_DIR, 'sensitivity_data.csv')
    sanity("output file exists", os.path.exists(out_path))

    # SCIENTIFIC: CV varies across steps (not all identical)
    for param in SENSITIVITY_PARAMS[:2]:  # Just check first 2 for speed
        sub = df[df['parameter'] == param]
        cv_range = sub['CV'].max() - sub['CV'].min()
        scientific(f"CV varies for {param}", cv_range > 0,
                   f"— range={cv_range:.4f}")

    print(f"\n{'='*50}")
    print(f"SANITY: {n_sanity_pass} pass, {n_sanity_fail} fail")
    print(f"SCIENTIFIC: {n_science_pass} pass, {n_science_warn} warn")
    if n_sanity_fail > 0:
        print("STOPPING: Sanity check failed.")
        sys.exit(1)
    else:
        print("All checks passed.")
        sys.exit(0)
