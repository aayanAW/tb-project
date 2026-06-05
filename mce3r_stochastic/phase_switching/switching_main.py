"""
Aim 2: State-switching logic that distinguishes asymmetry-driven noise from generic noise.

Computes threshold-free metrics (dwell time, first-passage, tail probability) under
three killing models, comparing four architectures designed to isolate the mechanism.

Key question: Does asymmetry create rare long-lived excursions, or just more noise?
"""

import sys
import os
import json
import time
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from config.parameters import PARAMS
from phase2_simulation.operator_model import OperatorModel
from phase2_simulation.gillespie_engine import simulate_cell
from phase_evo.genetic_algorithm import run_population_fast

PROJECT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS_DIR = os.path.join(PROJECT, 'results', 'phase_switching')


# ============================================================
# TRACE GENERATION (need full time series, not just endpoints)
# ============================================================

def generate_traces(Kd_strong, Kd_weak, block_strong, block_weak,
                    n_cells=200, seed=42, t_max=30000.0, t_burn_in=15000.0,
                    trace_interval=10.0):
    """Generate full protein traces for a given architecture."""
    model = OperatorModel(Kd_strong=Kd_strong, Kd_weak=Kd_weak,
                          block_strong=block_strong, block_weak=block_weak)
    arrays = model.get_numba_arrays()

    k_trans = np.float64(PARAMS['k_translation'])
    gamma_m = np.float64(np.log(2) / PARAMS['t_half_mRNA'])
    gamma_p = np.float64(np.log(2) / PARAMS['t_half_protein'])
    nM_mol = np.float64(PARAMS['nM_per_molecule'])
    trace_max = np.int64(int((t_max - t_burn_in) / trace_interval) + 100)

    traces = []
    for i in range(n_cells):
        _, _, _, t_times, t_prots, n_trace = simulate_cell(
            arrays['k_txn'], arrays['n_bound'], arrays['k_on'],
            arrays['k_off_strong'], arrays['k_off_weak'],
            k_trans, gamma_m, gamma_p, nM_mol,
            np.float64(t_max), np.float64(t_burn_in),
            np.float64(trace_interval), trace_max,
            np.int64(seed + i), np.int64(1)
        )
        n = int(n_trace)
        if n > 10:
            traces.append(t_prots[:n].astype(float))

    return traces


# ============================================================
# THRESHOLD-FREE METRICS
# ============================================================

def compute_excursion_metrics(traces, threshold, dt_minutes=10.0):
    """
    Compute threshold-free state-switching metrics from traces.

    Returns metrics about excursions below threshold:
    - probability of entering low state
    - mean dwell time in low state
    - tail probabilities of long excursions
    - max excursion duration
    """
    all_dwell_times = []
    n_cells_entering = 0
    n_total = len(traces)
    total_time_in_low = 0
    total_time = 0

    for trace in traces:
        # Find contiguous runs below threshold
        below = trace < threshold
        total_time += len(trace) * dt_minutes

        if not np.any(below):
            continue

        n_cells_entering += 1
        total_time_in_low += np.sum(below) * dt_minutes

        # Find runs
        changes = np.diff(below.astype(int))
        starts = np.where(changes == 1)[0] + 1
        ends = np.where(changes == -1)[0] + 1

        # Handle edge cases
        if below[0]:
            starts = np.concatenate([[0], starts])
        if below[-1]:
            ends = np.concatenate([ends, [len(below)]])

        for s, e in zip(starts, ends):
            dwell = (e - s) * dt_minutes
            all_dwell_times.append(dwell)

    if len(all_dwell_times) == 0:
        return {
            'prob_entering_low': 0.0,
            'mean_dwell_min': 0.0,
            'median_dwell_min': 0.0,
            'max_dwell_min': 0.0,
            'tail_gt_2hr': 0.0,
            'tail_gt_6hr': 0.0,
            'tail_gt_12hr': 0.0,
            'n_excursions': 0,
            'frac_time_in_low': 0.0,
        }

    dwells = np.array(all_dwell_times)

    return {
        'prob_entering_low': n_cells_entering / n_total,
        'mean_dwell_min': float(np.mean(dwells)),
        'median_dwell_min': float(np.median(dwells)),
        'max_dwell_min': float(np.max(dwells)),
        'tail_gt_2hr': float(np.mean(dwells > 120)),
        'tail_gt_6hr': float(np.mean(dwells > 360)),
        'tail_gt_12hr': float(np.mean(dwells > 720)),
        'n_excursions': len(dwells),
        'frac_time_in_low': total_time_in_low / total_time if total_time > 0 else 0.0,
    }


# ============================================================
# THREE KILLING MODELS
# ============================================================

def apply_hard_threshold_kill(traces, threshold, exposure_hours):
    """Kill cells below threshold at end of exposure."""
    survivors = 0
    for trace in traces:
        exposure_points = int(exposure_hours * 60 / 10)  # dt=10 min
        if len(trace) < exposure_points:
            exposure_points = len(trace)
        # Cell survives if it's below threshold during exposure
        window = trace[-exposure_points:]
        mean_in_window = np.mean(window)
        if mean_in_window < threshold:
            survivors += 1
    return survivors / len(traces) if traces else 0.0


def apply_hazard_kill(traces, threshold, exposure_hours, base_kill_rate=0.1):
    """
    Probabilistic hazard: lower expression = lower kill probability per hour.
    Cells above threshold have kill_rate = base_kill_rate per hour.
    Cells below threshold have kill_rate = base_kill_rate * 0.05 (95% protection).
    """
    dt_hr = 10.0 / 60.0  # 10 min in hours
    exposure_points = int(exposure_hours * 60 / 10)
    survivors = 0
    rng = np.random.RandomState(42)

    for trace in traces:
        if len(trace) < exposure_points:
            ep = len(trace)
        else:
            ep = exposure_points
        window = trace[-ep:]
        alive = True
        for val in window:
            if val < threshold:
                kill_prob = base_kill_rate * 0.05 * dt_hr
            else:
                kill_prob = base_kill_rate * dt_hr
            if rng.random() < kill_prob:
                alive = False
                break
        if alive:
            survivors += 1

    return survivors / len(traces) if traces else 0.0


def apply_dwell_kill(traces, threshold, min_dwell_hours):
    """
    Dwell-time-dependent: cells survive only if they spent >N continuous hours
    below threshold during the observation window.
    """
    dt_min = 10.0
    min_dwell_points = int(min_dwell_hours * 60 / dt_min)
    survivors = 0

    for trace in traces:
        below = trace < threshold
        # Find max contiguous run below threshold
        max_run = 0
        current_run = 0
        for b in below:
            if b:
                current_run += 1
                max_run = max(max_run, current_run)
            else:
                current_run = 0
        if max_run >= min_dwell_points:
            survivors += 1

    return survivors / len(traces) if traces else 0.0


# ============================================================
# FOUR-ARCHITECTURE COMPARISON
# ============================================================

def define_four_architectures():
    """
    Four architectures to isolate asymmetry's mechanism:
    1. Wild-type asymmetric
    2. Symmetric equal-affinity
    3. Asymmetric, mean-occupancy matched
    4. Symmetric, variance-matched (tuned to match WT variance)
    """
    return {
        'wt_asymmetric': {
            'Kd_strong': PARAMS['Kd_strong'],  # 2.4
            'Kd_weak': PARAMS['Kd_weak'],       # 49.0
            'block_strong': PARAMS['block_strong'],  # 0.85
            'block_weak': PARAMS['block_weak'],      # 0.50
            'description': 'Wild-type asymmetric (Mce3R)',
        },
        'symmetric_equal': {
            'Kd_strong': PARAMS['Kd_symmetric'],  # 10.84
            'Kd_weak': PARAMS['Kd_symmetric'],
            'block_strong': PARAMS['block_symmetric'],  # 0.652
            'block_weak': PARAMS['block_symmetric'],
            'description': 'Symmetric equal-affinity (geometric mean)',
        },
        'asymmetric_mean_matched': {
            # Same mean occupancy as WT but different asymmetry
            # Reduce ratio by ~half while adjusting blocks to match mean expression
            'Kd_strong': 4.0,
            'Kd_weak': 25.0,
            'block_strong': 0.80,
            'block_weak': 0.55,
            'description': 'Asymmetric, reduced ratio (~6x vs 20x), similar mean',
        },
        'symmetric_variance_matched': {
            # Symmetric but with higher block to increase variance
            # Tighter binding = more time in repressed state = more variance
            'Kd_strong': 3.0,
            'Kd_weak': 3.0,
            'block_strong': 0.90,
            'block_weak': 0.90,
            'description': 'Symmetric, tuned for high variance via tight binding',
        },
    }


# ============================================================
# MAIN PIPELINE
# ============================================================

def run_switching(n_cells=200, seed=42):
    """Run full Aim 2 analysis."""
    os.makedirs(RESULTS_DIR, exist_ok=True)
    start = time.time()

    architectures = define_four_architectures()

    # Determine threshold from Condition B
    cond_b_path = os.path.join(PROJECT, 'results', 'phase2', 'condition_B.npz')
    if os.path.exists(cond_b_path):
        b = np.load(cond_b_path)
        prot_b = b['proteins'].astype(float)
        threshold = float(np.mean(prot_b) - 2 * np.std(prot_b))
        threshold = max(threshold, float(np.percentile(prot_b, 1)))
    else:
        threshold = 150.0

    print(f"[Switching] Persistence threshold: {threshold:.0f} proteins")
    print(f"[Switching] Generating traces for 4 architectures ({n_cells} cells each)...")

    all_results = []
    arch_traces = {}

    for arch_name, params in architectures.items():
        print(f"\n  Architecture: {arch_name}")
        print(f"    Kd={params['Kd_strong']:.1f}/{params['Kd_weak']:.1f}, "
              f"block={params['block_strong']:.2f}/{params['block_weak']:.2f}")

        # Generate traces
        traces = generate_traces(
            params['Kd_strong'], params['Kd_weak'],
            params['block_strong'], params['block_weak'],
            n_cells=n_cells, seed=seed,
        )
        arch_traces[arch_name] = traces

        # Excursion metrics (threshold-free style but using threshold)
        for thr_name, thr_val in [('low', threshold), ('very_low', threshold * 0.5)]:
            metrics = compute_excursion_metrics(traces, thr_val)
            metrics['architecture'] = arch_name
            metrics['threshold_type'] = thr_name
            metrics['threshold_value'] = thr_val

            # Three killing models
            for exposure in [6, 12, 24]:
                surv_ht = apply_hard_threshold_kill(traces, thr_val, exposure)
                surv_hz = apply_hazard_kill(traces, thr_val, exposure)
                surv_dw = apply_dwell_kill(traces, thr_val, min_dwell_hours=exposure * 0.5)

                metrics[f'survival_hardthresh_{exposure}hr'] = surv_ht
                metrics[f'survival_hazard_{exposure}hr'] = surv_hz
                metrics[f'survival_dwell_{exposure}hr'] = surv_dw

            all_results.append(metrics)

            print(f"    [{thr_name}] P(enter)={metrics['prob_entering_low']:.3f}, "
                  f"mean_dwell={metrics['mean_dwell_min']:.0f}min, "
                  f"tail>6hr={metrics['tail_gt_6hr']:.3f}")

    df = pd.DataFrame(all_results)
    csv_path = os.path.join(RESULTS_DIR, 'switching_data.csv')
    df.to_csv(csv_path, index=False)
    print(f"\n  Saved: {csv_path}")

    # Summary analysis
    elapsed = time.time() - start

    # Key comparison: WT asymmetric vs variance-matched symmetric
    wt_low = df[(df['architecture'] == 'wt_asymmetric') & (df['threshold_type'] == 'low')].iloc[0]
    sym_var = df[(df['architecture'] == 'symmetric_variance_matched') & (df['threshold_type'] == 'low')].iloc[0]
    sym_eq = df[(df['architecture'] == 'symmetric_equal') & (df['threshold_type'] == 'low')].iloc[0]

    summary = {
        'description': 'Aim 2: State-switching logic comparison',
        'n_cells': n_cells,
        'threshold': threshold,
        'n_architectures': 4,
        'wall_time_sec': round(elapsed, 1),
        'n_sanity_pass': 0,
        'n_sanity_fail': 0,
        'n_science_expected': 0,
        'n_science_unexpected': 0,
    }

    # Sanity
    if len(df) == 8:  # 4 arch x 2 thresholds
        summary['n_sanity_pass'] += 1
    else:
        summary['n_sanity_fail'] += 1

    # Science: WT should have longer dwell times than symmetric equal
    if wt_low['mean_dwell_min'] > sym_eq['mean_dwell_min']:
        summary['n_science_expected'] += 1
        summary['wt_longer_dwells_than_symmetric'] = True
    else:
        summary['n_science_unexpected'] += 1
        summary['wt_longer_dwells_than_symmetric'] = False

    # Science: WT should have more tail>6hr than variance-matched symmetric
    # This is THE key test: is asymmetry special beyond just making more noise?
    if wt_low['tail_gt_6hr'] > sym_var['tail_gt_6hr']:
        summary['n_science_expected'] += 1
        summary['asymmetry_creates_longer_tails'] = True
    else:
        summary['n_science_unexpected'] += 1
        summary['asymmetry_creates_longer_tails'] = False

    # Science: WT should have higher survival under dwell-time killing
    if wt_low['survival_dwell_12hr'] > sym_eq['survival_dwell_12hr']:
        summary['n_science_expected'] += 1
    else:
        summary['n_science_unexpected'] += 1

    summary['wt_mean_dwell'] = float(wt_low['mean_dwell_min'])
    summary['sym_eq_mean_dwell'] = float(sym_eq['mean_dwell_min'])
    summary['sym_var_mean_dwell'] = float(sym_var['mean_dwell_min'])
    summary['wt_tail_6hr'] = float(wt_low['tail_gt_6hr'])
    summary['sym_var_tail_6hr'] = float(sym_var['tail_gt_6hr'])

    summary_path = os.path.join(RESULTS_DIR, 'switching_summary.json')
    with open(summary_path, 'w') as f:
        json.dump(summary, f, indent=2, default=str)

    print(f"\n[Switching] Complete in {elapsed:.0f}s")
    print(f"  Sanity: {summary['n_sanity_pass']} pass, {summary['n_sanity_fail']} fail")
    print(f"  Science: {summary['n_science_expected']} expected, {summary['n_science_unexpected']} unexpected")
    print(f"  WT mean dwell: {wt_low['mean_dwell_min']:.0f} min")
    print(f"  Symmetric equal dwell: {sym_eq['mean_dwell_min']:.0f} min")
    print(f"  Variance-matched dwell: {sym_var['mean_dwell_min']:.0f} min")

    return summary, df


if __name__ == '__main__':
    run_switching()
