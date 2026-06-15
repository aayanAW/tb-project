"""
phase6_environmental/run_environmental_conditions.py — Run all 3 arch × 4 env = 12 conditions.

For each of the 3 architectures × 4 environments, runs both:
  - Single-species simulations (existing gillespie_engine)
  - Two-species simulations (new two_species_gillespie)

Total: 24 conditions.

Also runs 5 trace cells (asymmetric, host_like, two-species) for fig14.

Output: results/phase6/env_condition_{arch}_{env}_{model}.npz  (24 files)
        results/phase6/two_species_results.npz
"""

import sys
import os
import time
import numpy as np

sys.path.insert(0, '/Users/aayanalwani/tb project/mce3r_stochastic')
from config.parameters import PARAMS
from phase2_simulation.operator_model import OperatorModel
from phase2_simulation.gillespie_engine import run_population
from phase6_environmental.two_species_model import TwoSpeciesModel
from phase6_environmental.two_species_gillespie import (
    run_two_species_population,
    simulate_two_species_cell_with_trace
)
from phase6_environmental.environmental_signals import apply_environment


OUTPUT_DIR = '/Users/aayanalwani/tb project/mce3r_stochastic/results/phase6'

ARCHITECTURES = {
    'asymmetric':  {'Kd_strong': 2.4,   'Kd_weak': 49.0},
    'symmetric':   {'Kd_strong': 10.84, 'Kd_weak': 10.84},
    'single_site': {'Kd_strong': 2.4,   'Kd_weak': 1e12},
}

ENVIRONMENTS = ['baseline', 'cholesterol', 'acidic_pH', 'host_like']


def ensure_output_dir():
    """Create output directory if it doesn't exist."""
    os.makedirs(OUTPUT_DIR, exist_ok=True)


def compute_stats(proteins):
    """Compute mean, std, CV, Fano for a protein array."""
    p = np.asarray(proteins, dtype=np.float64)
    mean = float(np.mean(p))
    std  = float(np.std(p))
    var  = float(np.var(p))
    cv   = std / mean if mean > 0 else 0.0
    fano = var / mean if mean > 0 else 0.0
    return {'mean': mean, 'std': std, 'cv': cv, 'fano': fano}


def run_all_environmental(
    n_cells: int = 10000,
    master_seed: int = 12345
) -> dict:
    """
    Run 3 architectures × 4 environments × 2 models = 24 conditions.

    Single-species: use existing gillespie_engine.run_population() with
                    model arrays modified by apply_environment().
    Two-species:    use two_species_gillespie.run_two_species_population().

    Also saves 5 trace cells for asymmetric, host_like, two-species (fig14).

    Parameters
    ----------
    n_cells : int
        Number of cells per condition.
    master_seed : int
        Master random seed; each condition gets a derived seed.

    Returns
    -------
    dict keyed by (architecture, environment, model_type) -> data_dict
    """
    ensure_output_dir()
    results = {}

    arch_list = list(ARCHITECTURES.keys())
    env_list = ENVIRONMENTS

    condition_index = 0
    for ai, arch in enumerate(arch_list):
        arch_params = ARCHITECTURES[arch]

        # Build single-species model arrays
        model_single = OperatorModel(
            Kd_strong=arch_params['Kd_strong'],
            Kd_weak=arch_params['Kd_weak']
        )
        base_single = model_single.get_numba_arrays()
        base_single['gamma_protein'] = np.float64(PARAMS['gamma_protein'])

        # Build two-species model
        model_two = TwoSpeciesModel(
            Kd_strong=arch_params['Kd_strong'],
            Kd_weak=arch_params['Kd_weak']
        )

        for ei, env in enumerate(env_list):
            condition_index += 1
            cell_seed_single = master_seed + ai * 100 + ei * 10
            cell_seed_two    = master_seed + ai * 100 + ei * 10 + 500

            # --- Single-species run ---
            print(f"  [{condition_index}/24] {arch} × {env} (single-species)...")
            t0 = time.time()
            env_single = apply_environment(base_single, env)
            res_single = run_population(
                env_single, n_cells=n_cells, master_seed=cell_seed_single
            )
            elapsed_s = time.time() - t0
            stats_s = compute_stats(res_single['proteins'])
            print(f"    Done in {elapsed_s:.1f}s — mean={stats_s['mean']:.1f}, "
                  f"CV={stats_s['cv']:.3f}")

            key_single = (arch, env, 'single')
            results[key_single] = {
                'proteins':    res_single['proteins'],
                'mRNAs':       res_single['mRNAs'],
                'op_states':   res_single['op_states'],
                'target_proteins': res_single['proteins'],  # alias for persistence_threshold
                'stats':       stats_s,
            }

            # Save NPZ
            fname_single = os.path.join(OUTPUT_DIR,
                                         f'env_condition_{arch}_{env}_single.npz')
            np.savez(fname_single,
                     proteins=res_single['proteins'],
                     mRNAs=res_single['mRNAs'],
                     op_states=res_single['op_states'])

            # --- Two-species run ---
            print(f"  [{condition_index}/24] {arch} × {env} (two-species)...")
            t0 = time.time()
            res_two = run_two_species_population(
                model_two, n_cells=n_cells,
                master_seed=cell_seed_two, environment=env
            )
            elapsed_t = time.time() - t0
            stats_t = compute_stats(res_two['target_proteins'])
            print(f"    Done in {elapsed_t:.1f}s — mean target={stats_t['mean']:.1f}, "
                  f"CV={stats_t['cv']:.3f}")

            key_two = (arch, env, 'two_species')
            results[key_two] = {**res_two, 'stats': stats_t}

            # Save NPZ
            fname_two = os.path.join(OUTPUT_DIR,
                                      f'env_condition_{arch}_{env}_two_species.npz')
            np.savez(fname_two,
                     mce3r_proteins=res_two['mce3r_proteins'],
                     target_proteins=res_two['target_proteins'],
                     mce3r_mRNAs=res_two['mce3r_mRNAs'],
                     target_mRNAs=res_two['target_mRNAs'],
                     op_states=res_two['op_states'])

    # --- Trace cells for fig14: asymmetric, host_like, two-species ---
    print("\n  Recording 5 trace cells (asymmetric, host_like, two-species)...")
    trace_arch = 'asymmetric'
    trace_env  = 'host_like'
    trace_model = TwoSpeciesModel(
        Kd_strong=ARCHITECTURES[trace_arch]['Kd_strong'],
        Kd_weak=ARCHITECTURES[trace_arch]['Kd_weak']
    )
    trace_arrays = trace_model.get_numba_arrays()
    trace_arrays_env = apply_environment(trace_arrays, trace_env)

    t_max     = np.float64(PARAMS['t_max'])
    t_burn_in = np.float64(PARAMS['t_burn_in'])
    trace_interval  = np.float64(PARAMS['trace_interval'])
    trace_max_pts   = np.int64(PARAMS['trace_max_points'])

    save_traces = {}
    for tc in range(5):
        tc_seed = master_seed + 9000 + tc
        (mp, tp, mm, tm, os_,
         tt, tmce3r, ttgt, n_tr) = simulate_two_species_cell_with_trace(
            trace_arrays_env['k_txn'],
            trace_arrays_env['n_bound'],
            trace_arrays_env['k_on'],
            trace_arrays_env['k_off_strong'],
            trace_arrays_env['k_off_weak'],
            trace_arrays_env['k_translation'],
            trace_arrays_env['gamma_mRNA'],
            trace_arrays_env['gamma_protein'],
            trace_arrays_env['nM_per_molecule'],
            trace_arrays_env['alpha_mce3r'],
            trace_arrays_env['alpha_target'],
            t_max, t_burn_in, trace_interval, trace_max_pts, np.int64(tc_seed)
        )
        save_traces[f'trace_time_{tc}']   = tt[:n_tr]
        save_traces[f'trace_mce3r_{tc}']  = tmce3r[:n_tr]
        save_traces[f'trace_target_{tc}'] = ttgt[:n_tr]
        save_traces[f'n_trace_{tc}']      = np.array([n_tr])
        print(f"    Trace cell {tc}: n_points={n_tr}, "
              f"final mce3r={mp}, target={tp}")

    save_traces['n_trace_cells'] = np.array([5])
    save_traces['architecture']  = np.array([trace_arch])
    save_traces['environment']   = np.array([trace_env])

    two_species_path = os.path.join(OUTPUT_DIR, 'two_species_results.npz')
    np.savez(two_species_path, **save_traces)
    print(f"  Saved two_species_results.npz")

    return results


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

    def scientific(name, cond, msg=""):
        global n_pass
        if cond:
            print(f"  SCIENTIFIC PASS: {name} {msg}")
            n_pass += 1
        else:
            print(f"  SCIENTIFIC WARN: {name} {msg}")

    print("=== run_environmental_conditions.py self-tests ===")
    print("  Running 24 conditions (n_cells=500 for fast testing)...")

    results = run_all_environmental(n_cells=500, master_seed=12345)

    # Test 1: All 24 conditions complete
    sanity("24 conditions in results dict",
           len(results) == 24,
           f"— got {len(results)} keys")

    # Test 2: All NPZ files exist and loadable
    all_files_ok = True
    for arch in ARCHITECTURES:
        for env in ENVIRONMENTS:
            for model_type in ['single', 'two_species']:
                fname = os.path.join(OUTPUT_DIR,
                                     f'env_condition_{arch}_{env}_{model_type}.npz')
                if not os.path.exists(fname):
                    all_files_ok = False
                    print(f"    MISSING: {fname}")
                else:
                    try:
                        d = np.load(fname)
                        _ = d.files
                    except Exception as e:
                        all_files_ok = False
                        print(f"    UNLOADABLE: {fname} — {e}")
    sanity("All 24 NPZ files exist and loadable", all_files_ok)

    two_species_path = os.path.join(OUTPUT_DIR, 'two_species_results.npz')
    sanity("two_species_results.npz exists",
           os.path.exists(two_species_path))

    # Test 3: CV(asymmetric) > CV(symmetric) under every environment
    # (using single-species results)
    cv_asym_gt_sym = True
    for env in ENVIRONMENTS:
        stats_asym = results[('asymmetric', env, 'single')]['stats']
        stats_sym  = results[('symmetric',  env, 'single')]['stats']
        if stats_asym['cv'] <= stats_sym['cv']:
            cv_asym_gt_sym = False
            print(f"    WARN: CV asymmetric ({stats_asym['cv']:.4f}) NOT > "
                  f"CV symmetric ({stats_sym['cv']:.4f}) under {env}")
        else:
            print(f"    {env}: CV_asym={stats_asym['cv']:.4f} > "
                  f"CV_sym={stats_sym['cv']:.4f} OK")
    scientific("CV(asymmetric) > CV(symmetric) under every environment",
               cv_asym_gt_sym)

    # Test 4: host_like mean < baseline mean for asymmetric architecture
    mean_asym_base = results[('asymmetric', 'baseline',  'single')]['stats']['mean']
    mean_asym_host = results[('asymmetric', 'host_like', 'single')]['stats']['mean']
    scientific("host_like mean > baseline mean for asymmetric (derepression)",
               mean_asym_host > mean_asym_base,
               f"— baseline={mean_asym_base:.1f}, host_like={mean_asym_host:.1f}")

    # Summary
    print(f"\n{'='*50}")
    print(f"SANITY: {n_pass} pass, {n_fail} fail")
    if n_fail > 0:
        print("STOPPING: Sanity check failed.")
        sys.exit(1)
    else:
        print("All checks passed.")
        sys.exit(0)
