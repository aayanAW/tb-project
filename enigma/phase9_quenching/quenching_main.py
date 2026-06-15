"""
Phase 9: Noise Quenching — Master Orchestrator

Coordinates the symmetrization sweep and IC50 computation.
Predicts the therapeutic dose-response for noise-reducing interventions.
"""

import sys
import os
import json
import time
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PHASE9_DIR = os.path.join(PROJECT_ROOT, 'results', 'phase9')


def run_phase9(n_doses=20, n_cells=5000):
    """Run all Phase 9 noise quenching analyses."""
    os.makedirs(PHASE9_DIR, exist_ok=True)
    start = time.time()

    summary = {
        'phase': 9,
        'description': 'Noise quenching: symmetrization sweep and IC50',
        'status': 'running',
        'n_sanity_pass': 0,
        'n_sanity_fail': 0,
        'n_science_expected': 0,
        'n_science_unexpected': 0,
        'output_files': [],
        'errors': [],
    }

    # --- Step 1: Symmetrization Sweep ---
    print("\n[Phase 9] Step 1/1: Symmetrization sweep...")
    try:
        from phase9_quenching.symmetrization_sweep import run_symmetrization_sweep
        df, ic50, info = run_symmetrization_sweep(n_doses=n_doses, n_cells=n_cells)

        summary['output_files'].append('symmetrization_sweep.csv')

        # Sanity: CV values all positive and finite
        if all(df['cv'] > 0) and all(np.isfinite(df['cv'])):
            summary['n_sanity_pass'] += 1
        else:
            summary['n_sanity_fail'] += 1

        # Sanity: CV decreases monotonically (approximately)
        cvs = df['cv'].values
        n_decreasing = sum(cvs[i] >= cvs[i+1] for i in range(len(cvs)-1))
        if n_decreasing >= len(cvs) * 0.6:  # at least 60% monotonic
            summary['n_sanity_pass'] += 1
        else:
            summary['n_sanity_fail'] += 1
            summary['errors'].append(f'CV not approximately monotonic: {n_decreasing}/{len(cvs)-1} decreasing')

        # Science: CV(dose=0) > CV(dose=1)
        if info['cv_asym'] > info['cv_sym']:
            summary['n_science_expected'] += 1
        else:
            summary['n_science_unexpected'] += 1

        # Science: IC50 should be in (0, 1)
        import numpy as np
        if np.isfinite(ic50) and 0 < ic50 < 1:
            summary['n_science_expected'] += 1
            summary['ic50_valid'] = True
        else:
            summary['n_science_unexpected'] += 1
            summary['ic50_valid'] = False

        # Record key results
        summary['cv_asymmetric'] = float(info['cv_asym'])
        summary['cv_symmetric'] = float(info['cv_sym'])
        summary['delta_cv'] = float(info['delta_cv_max'])
        summary['ic50_of_noise'] = float(ic50) if np.isfinite(ic50) else None
        summary['n_doses'] = n_doses
        summary['n_cells_per_dose'] = n_cells

        # Persister fraction at endpoints
        pf_start = float(df.iloc[0]['persister_fraction'])
        pf_end = float(df.iloc[-1]['persister_fraction'])
        summary['persister_fraction_asymmetric'] = pf_start
        summary['persister_fraction_symmetric'] = pf_end
        summary['persister_fold_reduction'] = pf_start / pf_end if pf_end > 0 else None

    except Exception as e:
        summary['errors'].append(f'Symmetrization sweep failed: {str(e)}')
        import traceback
        traceback.print_exc()

    # --- Finalize ---
    elapsed = time.time() - start
    summary['wall_time_sec'] = round(elapsed, 2)
    summary['status'] = 'complete' if not summary['errors'] else 'complete_with_errors'

    summary_path = os.path.join(PHASE9_DIR, 'phase9_summary.json')
    with open(summary_path, 'w') as f:
        json.dump(summary, f, indent=2, default=str)
    summary['output_files'].append('phase9_summary.json')

    print(f"\n[Phase 9] Complete in {elapsed:.1f}s")
    print(f"  Sanity: {summary['n_sanity_pass']} pass, {summary['n_sanity_fail']} fail")
    print(f"  Science: {summary['n_science_expected']} expected, {summary['n_science_unexpected']} unexpected")
    if summary['errors']:
        print(f"  Errors: {summary['errors']}")

    return summary


if __name__ == '__main__':
    run_phase9()
