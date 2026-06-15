"""
phase6_environmental/mutual_information.py — Mutual information estimation.

Computes I(Environment; Target_protein) — how much information about the
environmental condition is encoded in the target gene expression level.

Methods:
  - 'histogram': Classic binned estimator with equal-width bins.
  - 'ksg': Same as histogram but with adaptive bin count (Sturges rule).
           For a full KSG implementation, sklearn would be needed; here
           we use the adaptive histogram approach which is robust at our
           sample sizes.

Both methods add a small epsilon before log to avoid log(0).
"""

import sys
import numpy as np
import pandas as pd

sys.path.insert(0, '/Users/aayanalwani/tb project/mce3r_stochastic')
from config.parameters import PARAMS


def estimate_mutual_information(
    distributions: dict,
    n_bins: int = 50,
    method: str = 'histogram'
) -> float:
    """
    Estimate I(Environment; Target_protein) in bits.

    Parameters
    ----------
    distributions : dict
        Mapping env_name -> np.ndarray of protein counts (integers or floats).
        All arrays need not have the same length.
    n_bins : int
        Number of histogram bins (used for 'histogram' method).
        For 'ksg', bins are set adaptively via Sturges rule.
    method : str
        'histogram' — fixed n_bins.
        'ksg'       — adaptive bins (Sturges rule: 1 + log2(n_samples)).

    Returns
    -------
    float — MI in bits (>= 0)
    """
    env_names = list(distributions.keys())
    n_envs = len(env_names)
    if n_envs == 0:
        return 0.0

    # Pool all data to find global range for shared bin edges
    all_data = np.concatenate([distributions[e] for e in env_names])
    n_total = len(all_data)
    if n_total == 0:
        return 0.0

    # Adaptive bin count for 'ksg' method
    if method == 'ksg':
        n_bins = max(2, int(1 + np.log2(n_total)))

    data_min = float(np.min(all_data))
    data_max = float(np.max(all_data))

    # If all values identical, MI = 0
    if data_max == data_min:
        return 0.0

    # Build bin edges (shared across all environments)
    bin_edges = np.linspace(data_min, data_max, n_bins + 1)

    # P(env) = fraction of total samples from each environment
    n_per_env = np.array([len(distributions[e]) for e in env_names], dtype=np.float64)
    p_env = n_per_env / n_total  # shape: (n_envs,)

    # Histogram each environment with shared bin edges
    hist_matrix = np.zeros((n_envs, n_bins), dtype=np.float64)
    for i, env in enumerate(env_names):
        counts, _ = np.histogram(distributions[env], bins=bin_edges)
        hist_matrix[i, :] = counts.astype(np.float64)

    # P(bin | env) = hist_matrix[i, j] / n_per_env[i]
    # P(env, bin) = P(bin | env) * P(env) = hist_matrix[i, j] / n_total
    p_joint = hist_matrix / n_total  # shape: (n_envs, n_bins)

    # P(bin) = sum over envs of P(env, bin)
    p_bin = np.sum(p_joint, axis=0)  # shape: (n_bins,)

    # MI = sum_{env, bin} P(env, bin) * log2(P(env, bin) / (P(env) * P(bin)))
    eps = 1e-10
    mi = 0.0
    for i in range(n_envs):
        for j in range(n_bins):
            pij = p_joint[i, j]
            pi  = p_env[i]
            pj  = p_bin[j]
            if pij > 0 and pj > 0:
                mi += pij * np.log2((pij + eps) / (pi * pj + eps))

    # Clip to 0 to avoid tiny numerical negatives
    return float(max(0.0, mi))


def mutual_information_vs_concentration(
    architecture: str,
    concentrations_nM: np.ndarray,
    environments: list,
    n_cells: int = 5000,
    seed: int = 42
) -> pd.DataFrame:
    """
    Compute MI(Environment; Target_protein) across a range of Mce3R concentrations.

    For each concentration, the operator k_on is set to concentration-dependent value
    by treating effective_k_on = k_on * (concentration / reference_conc).
    We simulate the single-species model (target protein in each environment) and
    compute MI between the environment label and the protein count distribution.

    Parameters
    ----------
    architecture : str
        One of 'asymmetric', 'symmetric', 'single_site'.
    concentrations_nM : np.ndarray
        Array of Mce3R concentrations to sweep (nM).
    environments : list of str
        Environment names to include in MI computation.
    n_cells : int
        Number of cells per (concentration, environment) combination.
    seed : int
        Master random seed.

    Returns
    -------
    pd.DataFrame with columns: concentration, MI_bits, architecture
    """
    from phase6_environmental.environmental_signals import apply_environment, ENVIRONMENTS as ENV_DICT
    from phase2_simulation.operator_model import OperatorModel
    from phase2_simulation.gillespie_engine import run_population

    # Architecture Kd parameters
    arch_params = {
        'asymmetric':  {'Kd_strong': 2.4,   'Kd_weak': 49.0},
        'symmetric':   {'Kd_strong': 10.84,  'Kd_weak': 10.84},
        'single_site': {'Kd_strong': 2.4,   'Kd_weak': 1e12},
    }
    if architecture not in arch_params:
        raise ValueError(f"Unknown architecture '{architecture}'")

    params = arch_params[architecture]
    rows = []

    for ci, conc_nM in enumerate(concentrations_nM):
        # Scale k_on proportionally to concentration
        # Reference: PARAMS['k_on'] is already in nM^-1 min^-1 units.
        # We adjust so the effective binding rate at this concentration matches.
        # k_on is a second-order rate constant; we keep it fixed and instead
        # think of concentration as a multiplier on the effective on-rate.
        # Simplification: k_on_eff = k_on * conc_nM / reference_conc
        # (This is a physiological interpretation for population-level MI.)
        reference_conc = 332.0  # nM (from Phase 5 mcmc_summary.json)
        k_on_scaled = PARAMS['k_on'] * (conc_nM / reference_conc)

        model = OperatorModel(
            Kd_strong=params['Kd_strong'],
            Kd_weak=params['Kd_weak'],
            k_on=k_on_scaled
        )
        base_arrays = model.get_numba_arrays()
        base_arrays['gamma_protein'] = np.float64(PARAMS['gamma_protein'])

        distributions = {}
        for ei, env_name in enumerate(environments):
            env_arrays = apply_environment(base_arrays, env_name)
            cell_seed = seed + ci * 1000 + ei * 100
            result = run_population(env_arrays, n_cells=n_cells, master_seed=cell_seed)
            distributions[env_name] = result['proteins'].astype(np.float64)

        mi = estimate_mutual_information(distributions, method='histogram')
        rows.append({
            'concentration': float(conc_nM),
            'MI_bits': mi,
            'architecture': architecture,
        })

    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Self-tests
# ---------------------------------------------------------------------------
if __name__ == "__main__":
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

    print("=== mutual_information.py self-tests ===")

    # Test 1: MI >= 0
    np.random.seed(0)
    dist1 = {'A': np.random.normal(100, 10, 500), 'B': np.random.normal(200, 10, 500)}
    mi1 = estimate_mutual_information(dist1, n_bins=50, method='histogram')
    sanity("MI >= 0 for different distributions",
           mi1 >= 0,
           f"— MI = {mi1:.4f} bits")

    # Test 2: MI for identical distributions ≈ 0 (< 0.01 bits)
    identical_data = np.random.normal(100, 20, 1000)
    dist_identical = {'A': identical_data.copy(), 'B': identical_data.copy()}
    mi2 = estimate_mutual_information(dist_identical, n_bins=50, method='histogram')
    sanity("MI for identical distributions ≈ 0",
           mi2 < 0.01,
           f"— MI = {mi2:.6f} bits (< 0.01)")

    # Test 3: MI for very different distributions > 0
    dist_different = {
        'low':  np.ones(500) * 10,
        'high': np.ones(500) * 500,
    }
    mi3 = estimate_mutual_information(dist_different, n_bins=50, method='histogram')
    sanity("MI for very different distributions > 0",
           mi3 > 0,
           f"— MI = {mi3:.4f} bits")

    # Test 4: KSG method also returns non-negative MI
    mi_ksg = estimate_mutual_information(dist1, method='ksg')
    sanity("KSG method MI >= 0",
           mi_ksg >= 0,
           f"— MI_ksg = {mi_ksg:.4f} bits")

    # Summary
    print(f"\n{'='*50}")
    print(f"SANITY: {n_pass} pass, {n_fail} fail")
    if n_fail > 0:
        print("STOPPING: Sanity check failed.")
        sys.exit(1)
    else:
        print("All checks passed.")
        sys.exit(0)
