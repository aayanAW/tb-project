"""
Phase 8: Temporal Noise Dynamics — Master Orchestrator

Coordinates three temporal analyses:
1. Autocorrelation function and tau_c extraction
2. Operator state dwell times (theoretical)
3. Power spectral density analysis

Key scientific question: Does operator asymmetry create longer noise memory?
"""

import sys
import os
import json
import time

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PHASE8_DIR = os.path.join(PROJECT_ROOT, 'results', 'phase8')


def run_phase8():
    """Run all Phase 8 temporal analyses."""
    os.makedirs(PHASE8_DIR, exist_ok=True)
    start = time.time()

    summary = {
        'phase': 8,
        'description': 'Temporal noise dynamics: autocorrelation, dwell times, PSD',
        'status': 'running',
        'n_sanity_pass': 0,
        'n_sanity_fail': 0,
        'n_science_expected': 0,
        'n_science_unexpected': 0,
        'output_files': [],
        'errors': [],
    }

    # --- Step 1: Autocorrelation ---
    print("\n[Phase 8] Step 1/3: Autocorrelation analysis...")
    try:
        from phase8_temporal.autocorrelation import run_autocorrelation
        acf_results, acf_df = run_autocorrelation()

        summary['output_files'].append('autocorrelation.csv')
        summary['output_files'].append('acf_data.npz')

        # Sanity: tau_c values are finite and positive
        for cond in ['A', 'B', 'C']:
            if cond in acf_results:
                tc = acf_results[cond]['tau_c_mean_acf']
                if tc > 0 and tc < 100000:
                    summary['n_sanity_pass'] += 1
                else:
                    summary['n_sanity_fail'] += 1
                    summary['errors'].append(f'tau_c({cond}) out of range: {tc}')

        # Science: tau_c(A) comparison with tau_c(B)
        if 'A' in acf_results and 'B' in acf_results:
            tc_a = acf_results['A']['tau_c_mean_acf']
            tc_b = acf_results['B']['tau_c_mean_acf']
            summary['tau_c_A_min'] = float(tc_a)
            summary['tau_c_B_min'] = float(tc_b)
            summary['tau_c_A_hr'] = float(tc_a / 60)
            summary['tau_c_B_hr'] = float(tc_b / 60)
            summary['tau_c_ratio_A_over_B'] = float(tc_a / tc_b) if tc_b > 0 else None
            # Record whether asymmetric has longer memory
            if tc_a > tc_b:
                summary['n_science_expected'] += 1
                summary['tau_c_A_gt_B'] = True
            else:
                summary['n_science_unexpected'] += 1
                summary['tau_c_A_gt_B'] = False

    except Exception as e:
        summary['errors'].append(f'Autocorrelation failed: {str(e)}')
        import traceback
        traceback.print_exc()

    # --- Step 2: Dwell Times ---
    print("\n[Phase 8] Step 2/3: Dwell time analysis...")
    try:
        from phase8_temporal.dwell_times import run_dwell_analysis
        dwell_df = run_dwell_analysis()

        summary['output_files'].append('dwell_times.csv')

        # Sanity: all dwell times positive
        if all(dwell_df['dwell_state0_min'] > 0):
            summary['n_sanity_pass'] += 1
        else:
            summary['n_sanity_fail'] += 1

        # Science: asymmetric has different dwell distribution than symmetric
        row_a = dwell_df[dwell_df['condition'] == 'A'].iloc[0]
        row_b = dwell_df[dwell_df['condition'] == 'B'].iloc[0]

        # In asymmetric: State 1 (strong) dwell != State 2 (weak) dwell
        asym_ratio = row_a['dwell_state1_min'] / row_a['dwell_state2_min']
        sym_ratio = row_b['dwell_state1_min'] / row_b['dwell_state2_min']
        summary['dwell_asymmetry_ratio_A'] = float(asym_ratio)
        summary['dwell_asymmetry_ratio_B'] = float(sym_ratio)

        # For asymmetric: ratio should differ significantly from 1
        if abs(asym_ratio - 1.0) > abs(sym_ratio - 1.0):
            summary['n_science_expected'] += 1
            summary['dwell_asym_more_different'] = True
        else:
            summary['n_science_unexpected'] += 1
            summary['dwell_asym_more_different'] = False

    except Exception as e:
        summary['errors'].append(f'Dwell times failed: {str(e)}')
        import traceback
        traceback.print_exc()

    # --- Step 3: Power Spectral Density ---
    print("\n[Phase 8] Step 3/3: Power spectral density analysis...")
    try:
        from phase8_temporal.power_spectrum import run_psd_analysis
        psd_df, psd_data = run_psd_analysis()

        summary['output_files'].append('power_spectrum.csv')
        summary['output_files'].append('psd_data.npz')

        # Sanity: total power positive
        if all(psd_df['total_power'] > 0):
            summary['n_sanity_pass'] += 1
        else:
            summary['n_sanity_fail'] += 1

        # Science: asymmetric should have more low-frequency power
        row_a = psd_df[psd_df['condition'] == 'A'].iloc[0]
        row_b = psd_df[psd_df['condition'] == 'B'].iloc[0]
        summary['low_freq_power_A'] = float(row_a['low_freq_fraction_500min'])
        summary['low_freq_power_B'] = float(row_b['low_freq_fraction_500min'])

    except Exception as e:
        summary['errors'].append(f'PSD failed: {str(e)}')
        import traceback
        traceback.print_exc()

    # --- Finalize ---
    elapsed = time.time() - start
    summary['wall_time_sec'] = round(elapsed, 2)
    summary['status'] = 'complete' if not summary['errors'] else 'complete_with_errors'

    summary_path = os.path.join(PHASE8_DIR, 'phase8_summary.json')
    with open(summary_path, 'w') as f:
        json.dump(summary, f, indent=2, default=str)
    summary['output_files'].append('phase8_summary.json')

    print(f"\n[Phase 8] Complete in {elapsed:.1f}s")
    print(f"  Sanity: {summary['n_sanity_pass']} pass, {summary['n_sanity_fail']} fail")
    print(f"  Science: {summary['n_science_expected']} expected, {summary['n_science_unexpected']} unexpected")
    if summary['errors']:
        print(f"  Errors: {summary['errors']}")

    return summary


if __name__ == '__main__':
    run_phase8()
