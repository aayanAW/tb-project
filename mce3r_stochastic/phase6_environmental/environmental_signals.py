"""
phase6_environmental/environmental_signals.py — Environmental signal definitions.

Defines how cholesterol and acidic pH affect Mce3R effective concentration
and protein noise via multipliers on k_on and gamma_protein.

Biological basis:
- Cholesterol reduces Mce3R DNA-binding affinity (3-fold, mce3r_multiplier=0.3)
- Acidic pH (phagosomal) mildly reduces Mce3R activity + increases protein turnover
- Combined host-like stress = strong derepression + high noise
"""

import sys
import numpy as np

sys.path.insert(0, '/Users/aayanalwani/tb project/mce3r_stochastic')

ENVIRONMENTS = {
    'baseline': {
        'description': 'Glycerol, neutral pH (in vitro)',
        'mce3r_multiplier': 1.0,
        'noise_scale': 1.0,
    },
    'cholesterol': {
        'description': 'Cholesterol as sole carbon source',
        'mce3r_multiplier': 0.3,
        'noise_scale': 1.0,
    },
    'acidic_pH': {
        'description': 'pH 5.5 (phagosomal)',
        'mce3r_multiplier': 0.7,
        'noise_scale': 1.2,
    },
    'host_like': {
        'description': 'Cholesterol + acidic pH (macrophage phagosome)',
        'mce3r_multiplier': 0.2,
        'noise_scale': 1.3,
    },
}


def apply_environment(model_arrays: dict, environment: str) -> dict:
    """
    Return a COPY of model_arrays with environmental modifications applied.

    Modifications:
    - k_on scaled by mce3r_multiplier  (reduces effective Mce3R binding)
    - gamma_protein scaled by noise_scale  (increases protein turnover under stress)

    Does NOT mutate the input dict.

    Parameters
    ----------
    model_arrays : dict
        Dict from OperatorModel.get_numba_arrays() or TwoSpeciesModel.get_numba_arrays().
        Must contain 'k_on' (float64) and 'gamma_protein' (float64).
    environment : str
        One of 'baseline', 'cholesterol', 'acidic_pH', 'host_like'.

    Returns
    -------
    dict — modified copy of model_arrays
    """
    if environment not in ENVIRONMENTS:
        raise ValueError(f"Unknown environment '{environment}'. "
                         f"Must be one of: {list(ENVIRONMENTS.keys())}")

    env = ENVIRONMENTS[environment]
    mce3r_mult = env['mce3r_multiplier']
    noise_scale = env['noise_scale']

    # Deep copy: copy all values; arrays are copied, scalars are immutable
    modified = {}
    for key, val in model_arrays.items():
        if isinstance(val, np.ndarray):
            modified[key] = val.copy()
        else:
            modified[key] = val

    # Apply multipliers
    modified['k_on'] = np.float64(model_arrays['k_on'] * mce3r_mult)
    modified['gamma_protein'] = np.float64(model_arrays['gamma_protein'] * noise_scale)

    # Also update derived k_off values if present (they depend on k_on in binding engine)
    # Note: k_off_strong and k_off_weak are UNBINDING rates (Kd*k_on).
    # They should NOT be scaled by mce3r_multiplier — only k_on (binding rate) changes.
    # The multiplier reduces the forward binding rate, not the unbinding rate.
    # This is the correct biophysical interpretation: ligand reduces on-rate affinity.

    return modified


# ---------------------------------------------------------------------------
# Self-tests
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    import sys
    from config.parameters import PARAMS
    from phase2_simulation.operator_model import OperatorModel

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

    print("=== environmental_signals.py self-tests ===")

    model = OperatorModel()
    base_arrays = model.get_numba_arrays()
    # Add gamma_protein (not in OperatorModel by default, but needed for apply_environment)
    base_arrays['gamma_protein'] = np.float64(PARAMS['gamma_protein'])

    # Test 1: All environments produce non-negative rates
    all_non_neg = True
    for env_name in ENVIRONMENTS:
        modified = apply_environment(base_arrays, env_name)
        if modified['k_on'] < 0 or modified['gamma_protein'] < 0:
            all_non_neg = False
            print(f"    FAIL: {env_name} produced negative rate: "
                  f"k_on={modified['k_on']}, gamma_protein={modified['gamma_protein']}")
    sanity("All environments produce non-negative rates", all_non_neg)

    # Test 2: host_like k_on < baseline k_on
    baseline_mod = apply_environment(base_arrays, 'baseline')
    hostlike_mod = apply_environment(base_arrays, 'host_like')
    sanity("host_like k_on < baseline k_on",
           hostlike_mod['k_on'] < baseline_mod['k_on'],
           f"— host_like k_on={hostlike_mod['k_on']:.6f}, "
           f"baseline k_on={baseline_mod['k_on']:.6f}")

    # Test 3: Input not mutated
    original_k_on = base_arrays['k_on']
    _ = apply_environment(base_arrays, 'cholesterol')
    sanity("Input not mutated",
           base_arrays['k_on'] == original_k_on,
           f"— original k_on={original_k_on:.6f}")

    # Summary
    print(f"\n{'='*50}")
    print(f"SANITY: {n_pass} pass, {n_fail} fail")
    if n_fail > 0:
        print("STOPPING: Sanity check failed.")
        sys.exit(1)
    else:
        print("All checks passed.")
        sys.exit(0)
