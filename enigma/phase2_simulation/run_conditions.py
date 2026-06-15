"""
phase2_simulation/run_conditions.py — Run all 5 simulation conditions.

Conditions:
  A — Asymmetric (wild-type Mce3R operator)
  B — Symmetric (both sites use strong Kd)
  C — Single-site (weak site disabled via Kd=1e12)
  D — No regulation (k_on=0, k_txn=k_max for all states)
  E — Asymmetry sweep (10 Kd ratios, 1000 cells each)

Saves .npz files to results/phase2/ directory.
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


def ensure_output_dir():
    """Create output directory if it doesn't exist."""
    os.makedirs(OUTPUT_DIR, exist_ok=True)


def compute_stats(proteins):
    """Compute mean, std, CV, Fano for a protein array."""
    p = proteins.astype(np.float64)
    mean = np.mean(p)
    std = np.std(p)
    var = np.var(p)
    cv = std / mean if mean > 0 else 0.0
    fano = var / mean if mean > 0 else 0.0
    return {'mean': mean, 'std': std, 'cv': cv, 'fano': fano}


def run_condition_D(n_cells=None):
    """
    Condition D: No regulation (constitutive expression).
    k_on=0, k_txn=[k_max]*4 for all states.
    Seed: master_seed + 3.
    """
    if n_cells is None:
        n_cells = PARAMS['n_cells_main']

    print(f"  Condition D (no regulation): {n_cells} cells...")
    model = OperatorModel(k_on=PARAMS['k_on_disabled'])
    arrays = model.get_numba_arrays()
    # Override: all transcription rates = k_max (no blocking)
    arrays['k_txn'] = np.array([PARAMS['k_max']] * 4, dtype=np.float64)

    seed = PARAMS['master_seed'] + 3
    t0 = time.time()
    result = run_population(arrays, n_cells=n_cells, master_seed=seed)
    elapsed = time.time() - t0

    stats = compute_stats(result['proteins'])
    print(f"    Done in {elapsed:.1f}s — mean={stats['mean']:.1f}, "
          f"CV={stats['cv']:.3f}, Fano={stats['fano']:.2f}")

    # Save
    ensure_output_dir()
    np.savez(os.path.join(OUTPUT_DIR, 'condition_D.npz'),
             proteins=result['proteins'],
             mRNAs=result['mRNAs'],
             op_states=result['op_states'])

    return {**result, 'stats': stats, 'condition': 'D', 'elapsed': elapsed}


def run_condition_B(n_cells=None):
    """
    Condition B: Symmetric operator (both sites strong Kd).
    Seed: master_seed + 1.
    """
    if n_cells is None:
        n_cells = PARAMS['n_cells_main']

    print(f"  Condition B (symmetric): {n_cells} cells...")
    model = OperatorModel(
        Kd_strong=PARAMS['Kd_symmetric'],
        Kd_weak=PARAMS['Kd_symmetric'],
        block_strong=PARAMS['block_symmetric'],
        block_weak=PARAMS['block_symmetric']
    )
    arrays = model.get_numba_arrays()

    seed = PARAMS['master_seed'] + 1
    n_trace = PARAMS['n_trace_cells']
    t0 = time.time()
    result = run_population(arrays, n_cells=n_cells, master_seed=seed,
                            record_traces=1, n_trace_cells=n_trace)
    elapsed = time.time() - t0

    stats = compute_stats(result['proteins'])
    print(f"    Done in {elapsed:.1f}s — mean={stats['mean']:.1f}, "
          f"CV={stats['cv']:.3f}, Fano={stats['fano']:.2f}")

    # Save
    ensure_output_dir()
    save_dict = {
        'proteins': result['proteins'],
        'mRNAs': result['mRNAs'],
        'op_states': result['op_states'],
    }
    # Save traces
    for i, (tt, tp) in enumerate(result['traces']):
        save_dict[f'trace_time_{i}'] = tt
        save_dict[f'trace_prot_{i}'] = tp
    save_dict['n_traces'] = np.array([len(result['traces'])])

    np.savez(os.path.join(OUTPUT_DIR, 'condition_B.npz'), **save_dict)

    return {**result, 'stats': stats, 'condition': 'B', 'elapsed': elapsed}


def run_condition_C(n_cells=None):
    """
    Condition C: Single-site (weak site disabled via Kd=1e12).
    Seed: master_seed + 2.
    """
    if n_cells is None:
        n_cells = PARAMS['n_cells_main']

    print(f"  Condition C (single-site): {n_cells} cells...")
    model = OperatorModel(Kd_weak=PARAMS['Kd_weak_disabled'])
    arrays = model.get_numba_arrays()

    seed = PARAMS['master_seed'] + 2
    n_trace = PARAMS['n_trace_cells']
    t0 = time.time()
    result = run_population(arrays, n_cells=n_cells, master_seed=seed,
                            record_traces=1, n_trace_cells=n_trace)
    elapsed = time.time() - t0

    stats = compute_stats(result['proteins'])
    print(f"    Done in {elapsed:.1f}s — mean={stats['mean']:.1f}, "
          f"CV={stats['cv']:.3f}, Fano={stats['fano']:.2f}")

    # Save
    ensure_output_dir()
    save_dict = {
        'proteins': result['proteins'],
        'mRNAs': result['mRNAs'],
        'op_states': result['op_states'],
    }
    for i, (tt, tp) in enumerate(result['traces']):
        save_dict[f'trace_time_{i}'] = tt
        save_dict[f'trace_prot_{i}'] = tp
    save_dict['n_traces'] = np.array([len(result['traces'])])

    np.savez(os.path.join(OUTPUT_DIR, 'condition_C.npz'), **save_dict)

    return {**result, 'stats': stats, 'condition': 'C', 'elapsed': elapsed}


def run_condition_A(n_cells=None):
    """
    Condition A: Asymmetric (wild-type Mce3R operator).
    Seed: master_seed + 0.
    """
    if n_cells is None:
        n_cells = PARAMS['n_cells_main']

    print(f"  Condition A (asymmetric): {n_cells} cells...")
    model = OperatorModel()  # default = wild-type
    arrays = model.get_numba_arrays()

    seed = PARAMS['master_seed'] + 0
    n_trace = PARAMS['n_trace_cells']
    t0 = time.time()
    result = run_population(arrays, n_cells=n_cells, master_seed=seed,
                            record_traces=1, n_trace_cells=n_trace)
    elapsed = time.time() - t0

    stats = compute_stats(result['proteins'])
    print(f"    Done in {elapsed:.1f}s — mean={stats['mean']:.1f}, "
          f"CV={stats['cv']:.3f}, Fano={stats['fano']:.2f}")

    # Save
    ensure_output_dir()
    save_dict = {
        'proteins': result['proteins'],
        'mRNAs': result['mRNAs'],
        'op_states': result['op_states'],
    }
    for i, (tt, tp) in enumerate(result['traces']):
        save_dict[f'trace_time_{i}'] = tt
        save_dict[f'trace_prot_{i}'] = tp
    save_dict['n_traces'] = np.array([len(result['traces'])])

    np.savez(os.path.join(OUTPUT_DIR, 'condition_A.npz'), **save_dict)

    return {**result, 'stats': stats, 'condition': 'A', 'elapsed': elapsed}


def run_condition_E(n_cells_per_ratio=None):
    """
    Condition E: Asymmetry sweep across 10 Kd ratios.
    Delegates to asymmetry_sweep module.
    """
    if n_cells_per_ratio is None:
        n_cells_per_ratio = PARAMS['n_cells_sweep']

    from phase2_simulation.asymmetry_sweep import run_sweep

    print(f"  Condition E (asymmetry sweep): {len(PARAMS['sweep_ratios'])} ratios × {n_cells_per_ratio} cells...")
    t0 = time.time()
    sweep_result = run_sweep(n_cells_per_ratio=n_cells_per_ratio)
    elapsed = time.time() - t0

    print(f"    Done in {elapsed:.1f}s")
    return {**sweep_result, 'condition': 'E', 'elapsed': elapsed}


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

    print("=== run_conditions.py self-tests ===")

    # Use small cell counts for testing
    N_TEST = 10

    # Run D first (validates engine)
    res_d = run_condition_D(n_cells=N_TEST)
    sanity("D: proteins shape", res_d['proteins'].shape == (N_TEST,))
    sanity("D: output file exists",
           os.path.exists(os.path.join(OUTPUT_DIR, 'condition_D.npz')))

    # Run B
    res_b = run_condition_B(n_cells=N_TEST)
    sanity("B: proteins shape", res_b['proteins'].shape == (N_TEST,))
    sanity("B: traces recorded", len(res_b['traces']) > 0,
           f"— {len(res_b['traces'])} traces")
    sanity("B: output file exists",
           os.path.exists(os.path.join(OUTPUT_DIR, 'condition_B.npz')))

    # Run C
    res_c = run_condition_C(n_cells=N_TEST)
    sanity("C: proteins shape", res_c['proteins'].shape == (N_TEST,))
    sanity("C: output file exists",
           os.path.exists(os.path.join(OUTPUT_DIR, 'condition_C.npz')))

    # Run A
    res_a = run_condition_A(n_cells=N_TEST)
    sanity("A: proteins shape", res_a['proteins'].shape == (N_TEST,))
    sanity("A: output file exists",
           os.path.exists(os.path.join(OUTPUT_DIR, 'condition_A.npz')))

    # Scientific: CV(A) > CV(B) — this is the key biological prediction
    # With 10 cells, this may not hold, so it's a warning
    scientific("CV(A) > CV(B) (key prediction)",
               res_a['stats']['cv'] > res_b['stats']['cv'],
               f"— CV_A={res_a['stats']['cv']:.3f}, CV_B={res_b['stats']['cv']:.3f}")

    # Scientific: unregulated mean > regulated means
    scientific("D mean > A mean",
               res_d['stats']['mean'] > res_a['stats']['mean'],
               f"— D={res_d['stats']['mean']:.1f}, A={res_a['stats']['mean']:.1f}")

    # Summary
    print(f"\n{'='*50}")
    print(f"SANITY: {n_pass} pass, {n_fail} fail | SCIENTIFIC: {n_warn} warn")
    if n_fail > 0:
        print("STOPPING: Sanity check failed.")
        sys.exit(1)
    else:
        print("All checks passed.")
        sys.exit(0)
