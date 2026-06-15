"""
Evolutionary simulation orchestrator.

Runs genetic algorithm under two selection regimes:
  1. growth_only — favors high expression (should converge toward symmetric/minimal regulation)
  2. persistence — favors growth + persister fraction (should converge toward asymmetry)

Compares final evolved architectures to Mce3R wild-type (Kd_ratio = 20.4).
"""

import sys
import os
import json
import time
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS_DIR = os.path.join(PROJECT_ROOT, 'results', 'phase_evo')


def run_phase_evo(pop_size=50, n_generations=100, n_cells=200, seed=42):
    """Run evolutionary simulation under both selection regimes."""
    os.makedirs(RESULTS_DIR, exist_ok=True)
    start = time.time()

    from phase_evo.genetic_algorithm import run_ga
    from config.parameters import PARAMS

    mce3r_ratio = PARAMS['Kd_weak'] / PARAMS['Kd_strong']  # 20.4
    persister_threshold = 112.0

    summary = {
        'description': 'Evolutionary simulation of operator architecture under selection',
        'pop_size': pop_size,
        'n_generations': n_generations,
        'n_cells_per_eval': n_cells,
        'mce3r_wt_ratio': mce3r_ratio,
        'persister_threshold': persister_threshold,
        'n_sanity_pass': 0,
        'n_sanity_fail': 0,
        'n_science_expected': 0,
        'n_science_unexpected': 0,
        'errors': [],
    }

    # --- Regime 1: Growth Only ---
    print("\n[Evo] Regime 1/2: Growth-only selection...")
    try:
        t0 = time.time()
        traj_g, pop_g = run_ga(
            regime='growth_only', pop_size=pop_size, n_generations=n_generations,
            n_cells=n_cells, seed=seed, persister_threshold=persister_threshold,
        )
        elapsed_g = time.time() - t0

        traj_g.to_csv(os.path.join(RESULTS_DIR, 'trajectory_growth.csv'), index=False)
        pop_g.to_csv(os.path.join(RESULTS_DIR, 'final_pop_growth.csv'), index=False)

        summary['growth_elapsed_sec'] = round(elapsed_g, 1)
        summary['growth_final_mean_ratio'] = float(pop_g['Kd_ratio'].mean())
        summary['growth_final_mean_cv'] = float(pop_g['cv'].mean())
        summary['growth_final_best_fitness'] = float(pop_g['fitness'].max())

        # Sanity: fitness values finite
        if all(np.isfinite(traj_g['mean_fitness'])):
            summary['n_sanity_pass'] += 1
        else:
            summary['n_sanity_fail'] += 1
    except Exception as e:
        summary['errors'].append(f'Growth regime failed: {e}')
        import traceback; traceback.print_exc()

    # --- Regime 2: Persistence ---
    print("\n[Evo] Regime 2/2: Persistence selection...")
    try:
        t0 = time.time()
        traj_p, pop_p = run_ga(
            regime='persistence', pop_size=pop_size, n_generations=n_generations,
            n_cells=n_cells, seed=seed + 100000, persister_threshold=persister_threshold,
        )
        elapsed_p = time.time() - t0

        traj_p.to_csv(os.path.join(RESULTS_DIR, 'trajectory_persistence.csv'), index=False)
        pop_p.to_csv(os.path.join(RESULTS_DIR, 'final_pop_persistence.csv'), index=False)

        summary['persist_elapsed_sec'] = round(elapsed_p, 1)
        summary['persist_final_mean_ratio'] = float(pop_p['Kd_ratio'].mean())
        summary['persist_final_mean_cv'] = float(pop_p['cv'].mean())
        summary['persist_final_best_fitness'] = float(pop_p['fitness'].max())

        if all(np.isfinite(traj_p['mean_fitness'])):
            summary['n_sanity_pass'] += 1
        else:
            summary['n_sanity_fail'] += 1
    except Exception as e:
        summary['errors'].append(f'Persistence regime failed: {e}')
        import traceback; traceback.print_exc()

    # --- Compare Regimes ---
    try:
        g_ratio = summary.get('growth_final_mean_ratio', 1.0)
        p_ratio = summary.get('persist_final_mean_ratio', 1.0)

        # Science: persistence regime should evolve higher Kd ratio (more asymmetric)
        if p_ratio > g_ratio:
            summary['n_science_expected'] += 1
            summary['persist_more_asymmetric'] = True
        else:
            summary['n_science_unexpected'] += 1
            summary['persist_more_asymmetric'] = False

        # Science: persistence-evolved ratio should be closer to Mce3R WT
        g_dist = abs(np.log(g_ratio / mce3r_ratio))
        p_dist = abs(np.log(p_ratio / mce3r_ratio))
        if p_dist < g_dist:
            summary['n_science_expected'] += 1
            summary['persist_closer_to_wt'] = True
        else:
            summary['n_science_unexpected'] += 1
            summary['persist_closer_to_wt'] = False

        summary['ratio_comparison'] = {
            'growth_mean_ratio': float(g_ratio),
            'persist_mean_ratio': float(p_ratio),
            'mce3r_wt_ratio': mce3r_ratio,
            'growth_log_dist_to_wt': float(g_dist),
            'persist_log_dist_to_wt': float(p_dist),
        }

        print(f"\n  Growth final mean Kd_ratio:      {g_ratio:.1f}")
        print(f"  Persistence final mean Kd_ratio: {p_ratio:.1f}")
        print(f"  Mce3R wild-type Kd_ratio:        {mce3r_ratio:.1f}")
    except Exception as e:
        summary['errors'].append(f'Comparison failed: {e}')

    # --- Finalize ---
    elapsed = time.time() - start
    summary['wall_time_sec'] = round(elapsed, 1)
    summary['status'] = 'complete' if not summary['errors'] else 'complete_with_errors'

    summary_path = os.path.join(RESULTS_DIR, 'evolution_summary.json')
    with open(summary_path, 'w') as f:
        json.dump(summary, f, indent=2, default=str)

    print(f"\n[Evo] Complete in {elapsed:.0f}s")
    print(f"  Sanity: {summary['n_sanity_pass']} pass, {summary['n_sanity_fail']} fail")
    print(f"  Science: {summary['n_science_expected']} expected, {summary['n_science_unexpected']} unexpected")

    return summary


if __name__ == '__main__':
    run_phase_evo()
