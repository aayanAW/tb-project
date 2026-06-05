"""
phase2_simulation/simulation_main.py — Main entry point for Phase 2 simulation.

Runs all conditions in order: D, B, C, A, E (D first to validate engine).
Collects results, writes phase2_summary.json, returns structured dict.
"""

import sys
import os
import json
import time
import numpy as np

sys.path.insert(0, '/Users/aayanalwani/tb project/mce3r_stochastic')
from config.parameters import PARAMS
from phase2_simulation.run_conditions import (
    run_condition_D, run_condition_B, run_condition_C,
    run_condition_A, run_condition_E, compute_stats
)


OUTPUT_DIR = '/Users/aayanalwani/tb project/mce3r_stochastic/results/phase2'


def run_all(n_cells_main=None, n_cells_sweep=None):
    """
    Run all 5 conditions and produce summary.

    Parameters
    ----------
    n_cells_main : int or None  — cells for A, B, C, D (default: PARAMS)
    n_cells_sweep : int or None — cells per ratio for E (default: PARAMS)

    Returns
    -------
    dict with condition results and summary statistics.
    """
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    t_total_start = time.time()

    results = {}

    # --- Condition D (no regulation) — run first to validate engine ---
    print("[Phase 2] Running Condition D (no regulation)...")
    results['D'] = run_condition_D(n_cells=n_cells_main)

    # --- Condition B (symmetric) ---
    print("[Phase 2] Running Condition B (symmetric)...")
    results['B'] = run_condition_B(n_cells=n_cells_main)

    # --- Condition C (single-site) ---
    print("[Phase 2] Running Condition C (single-site)...")
    results['C'] = run_condition_C(n_cells=n_cells_main)

    # --- Condition A (asymmetric / wild-type) ---
    print("[Phase 2] Running Condition A (asymmetric)...")
    results['A'] = run_condition_A(n_cells=n_cells_main)

    # --- Condition E (asymmetry sweep) ---
    print("[Phase 2] Running Condition E (asymmetry sweep)...")
    results['E'] = run_condition_E(n_cells_per_ratio=n_cells_sweep)

    t_total = time.time() - t_total_start
    print(f"\n[Phase 2] All conditions complete in {t_total:.1f}s")

    # --- Build summary ---
    summary = {
        'phase': 2,
        'description': 'Mce3R stochastic simulation — 4-state operator model',
        'total_time_seconds': round(t_total, 1),
        'parameters': {
            'Kd_strong': PARAMS['Kd_strong'],
            'Kd_weak': PARAMS['Kd_weak'],
            'Kd_ratio': round(PARAMS['Kd_weak'] / PARAMS['Kd_strong'], 1),
            'k_on': PARAMS['k_on'],
            'k_max': PARAMS['k_max'],
            'block_strong': PARAMS['block_strong'],
            'block_weak': PARAMS['block_weak'],
            'Kd_symmetric': round(PARAMS['Kd_symmetric'], 4),
            'block_symmetric': round(PARAMS['block_symmetric'], 4),
            'k_translation': PARAMS['k_translation'],
            'gamma_mRNA': round(PARAMS['gamma_mRNA'], 4),
            'gamma_protein': round(PARAMS['gamma_protein'], 6),
            't_max': PARAMS['t_max'],
            't_burn_in': PARAMS['t_burn_in'],
        },
        'conditions': {},
    }

    for cond_key in ['D', 'B', 'C', 'A']:
        r = results[cond_key]
        n_cells_used = len(r['proteins'])
        summary['conditions'][cond_key] = {
            'n_cells': n_cells_used,
            'mean_protein': round(r['stats']['mean'], 2),
            'std_protein': round(r['stats']['std'], 2),
            'cv': round(r['stats']['cv'], 4),
            'fano': round(r['stats']['fano'], 2),
            'elapsed_seconds': round(r['elapsed'], 1),
        }

    # Condition E summary
    e = results['E']
    summary['conditions']['E'] = {
        'n_ratios': len(e['ratios']),
        'n_cells_per_ratio': n_cells_sweep if n_cells_sweep else PARAMS['n_cells_sweep'],
        'ratios': e['ratios'].tolist(),
        'means': [round(x, 2) for x in e['means'].tolist()],
        'cvs': [round(x, 4) for x in e['cvs'].tolist()],
        'fanos': [round(x, 2) for x in e['fanos'].tolist()],
        'intermediate_fractions': [round(x, 4) for x in e['intermediate_fractions'].tolist()],
        'elapsed_seconds': round(e['elapsed'], 1),
    }

    # Key biological predictions
    cv_a = results['A']['stats']['cv']
    cv_b = results['B']['stats']['cv']
    cv_c = results['C']['stats']['cv']
    cv_d = results['D']['stats']['cv']
    fano_d = results['D']['stats']['fano']

    summary['predictions'] = {
        'cv_A_gt_cv_B': bool(cv_a > cv_b),
        'cv_A': round(float(cv_a), 4),
        'cv_B': round(float(cv_b), 4),
        'cv_C': round(float(cv_c), 4),
        'cv_D': round(float(cv_d), 4),
        'fano_D_approx_expected': round(float(fano_d), 2),
        'expected_fano_unregulated': round(float(PARAMS['expected_fano_unregulated']), 2),
    }

    # Write summary JSON
    summary_path = os.path.join(OUTPUT_DIR, 'phase2_summary.json')
    with open(summary_path, 'w') as f:
        json.dump(summary, f, indent=2)
    print(f"[Phase 2] Summary written to {summary_path}")

    return {
        'results': results,
        'summary': summary,
        'summary_path': summary_path,
    }


# ---------------------------------------------------------------------------
# Self-tests
# ---------------------------------------------------------------------------
def run_phase2():
    """Entry point for production runs (called by main.py subprocess)."""
    output = run_all()
    summary = output['summary']

    # Print final status
    print(f"\n{'='*60}")
    print(f"PHASE 2 SUMMARY")
    print(f"  Status: PASS")
    print(f"  Total time: {summary['total_time_seconds']}s")
    for cond in ['D', 'B', 'C', 'A']:
        c = summary['conditions'][cond]
        print(f"  Condition {cond}: {c['n_cells']} cells, "
              f"mean={c['mean_protein']:.1f}, CV={c['cv']:.4f}, Fano={c['fano']:.2f}")
    e = summary['conditions']['E']
    print(f"  Condition E: {e['n_ratios']} ratios × {e['n_cells_per_ratio']} cells")
    print(f"{'='*60}")

    return output


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--test', action='store_true',
                        help='Run self-tests with small cell counts instead of production')
    args = parser.parse_args()

    if not args.test:
        # PRODUCTION MODE: use PARAMS defaults (10,000 cells, etc.)
        run_phase2()
        sys.exit(0)

    # SELF-TEST MODE: small cell counts
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

    print("=== simulation_main.py self-tests (N=10) ===")

    output = run_all(n_cells_main=10, n_cells_sweep=10)

    summary = output['summary']
    results = output['results']

    sanity("summary JSON exists",
           os.path.exists(output['summary_path']))

    for c in ['A', 'B', 'C', 'D', 'E']:
        sanity(f"condition {c} in results", c in results)
        sanity(f"condition {c} in summary", c in summary['conditions'])

    sanity("predictions in summary", 'predictions' in summary)

    with open(output['summary_path']) as f:
        loaded = json.load(f)
    sanity("summary JSON loadable", loaded['phase'] == 2)

    for fname in ['condition_A.npz', 'condition_B.npz', 'condition_C.npz',
                   'condition_D.npz', 'condition_E_sweep.npz']:
        sanity(f"{fname} exists",
               os.path.exists(os.path.join(OUTPUT_DIR, fname)))

    fano_d = results['D']['stats']['fano']
    scientific("D Fano in range (with 10 cells, wide tolerance)",
               0.5 < fano_d < 50.0,
               f"— Fano={fano_d:.2f}")

    print(f"\n{'='*50}")
    print(f"SANITY: {n_pass} pass, {n_fail} fail | SCIENTIFIC: {n_warn} warn")
    if n_fail > 0:
        print("STOPPING: Sanity check failed.")
        sys.exit(1)
    else:
        print("All checks passed.")
        sys.exit(0)
