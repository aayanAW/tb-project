"""
Unified figure generation module for all 19 publication figures.

Replaces the old phase4_figures, phase7_extended_figures, and phase10_new_figures
modules with a single orchestrator.
"""

import sys
import os
import time
import json
import traceback

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT_DIR = os.path.join(BASE, 'results', 'figures')

FIGURES = [
    'fig1_binding_sites',
    'fig2_distributions',
    'fig3_asymmetry_sweep',
    'fig4_sensitivity',
    'fig5_validation',
    'fig6_single_cell_traces',
    'fig7_bic_comparison',
    'fig8_methods_architecture',
    'fig9_repression_curves',
    'fig10_cooperativity',
    'fig11_environmental_distributions',
    'fig12_mutual_information',
    'fig13_persistence_phase',
    'fig14_two_species_traces',
    'fig15_autocorrelation',
    'fig16_dwell_times',
    'fig17_power_spectrum',
    'fig18_dose_response',
    'fig19_summary_dashboard',
]


def run_figures(figures=None):
    """
    Generate all (or selected) figures.

    Parameters
    ----------
    figures : list of str, optional
        Figure names to generate (e.g. ['fig1_binding_sites', 'fig2_distributions']).
        If None, generates all 19.

    Returns
    -------
    dict : summary with n_success, n_error, per-figure status.
    """
    os.makedirs(OUT_DIR, exist_ok=True)

    from figures.figure_style import apply_style
    apply_style()

    targets = figures or FIGURES
    results = {}
    n_success = 0
    n_error = 0
    start = time.time()

    for name in targets:
        print(f"  [{name}] generating...")
        t0 = time.time()
        try:
            mod = __import__(f'figures.{name}', fromlist=['generate'])
            result = mod.generate()
            elapsed = time.time() - t0

            # Check output exists
            expected_png = os.path.join(OUT_DIR, f'{name}.png')
            if isinstance(result, str) and os.path.exists(result):
                expected_png = result
            elif isinstance(result, dict) and result.get('output_file'):
                expected_png = result['output_file']

            if os.path.exists(expected_png) and os.path.getsize(expected_png) > 1000:
                results[name] = {'status': 'success', 'elapsed': round(elapsed, 2)}
                n_success += 1
            else:
                results[name] = {'status': 'no_output', 'elapsed': round(elapsed, 2)}
                n_error += 1
        except Exception as e:
            elapsed = time.time() - t0
            results[name] = {'status': 'error', 'error': str(e), 'elapsed': round(elapsed, 2)}
            n_error += 1
            traceback.print_exc()

    total_elapsed = time.time() - start

    summary = {
        'n_figures': len(targets),
        'n_success': n_success,
        'n_error': n_error,
        'wall_time_sec': round(total_elapsed, 2),
        'figures': results,
    }

    summary_path = os.path.join(OUT_DIR, 'figures_summary.json')
    with open(summary_path, 'w') as f:
        json.dump(summary, f, indent=2, default=str)

    print(f"\n  Figures complete: {n_success}/{len(targets)} success, "
          f"{n_error} error, {total_elapsed:.1f}s total")
    return summary


# Backward compatibility alias
run_phase4 = run_figures


if __name__ == '__main__':
    run_figures()
