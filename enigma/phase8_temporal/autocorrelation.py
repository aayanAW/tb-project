"""
Phase 8: Autocorrelation analysis of protein expression traces.

Computes the normalized autocorrelation function C(tau) for each condition
and extracts the autocorrelation time tau_c (time for C to decay to 1/e).

Key question: Does operator asymmetry create longer noise memory?
"""

import sys
import os
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from config.parameters import PARAMS


PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PHASE2_DIR = os.path.join(PROJECT_ROOT, 'results', 'phase2')
PHASE8_DIR = os.path.join(PROJECT_ROOT, 'results', 'phase8')


def load_traces(condition, max_traces=20):
    """Load protein traces from Phase 2 NPZ files."""
    path = os.path.join(PHASE2_DIR, f'condition_{condition}.npz')
    if not os.path.exists(path):
        return None, None
    data = np.load(path, allow_pickle=True)

    n_traces_arr = data.get('n_traces', None)
    if n_traces_arr is not None:
        n_traces = int(n_traces_arr[0]) if hasattr(n_traces_arr, '__len__') else int(n_traces_arr)
    else:
        n_traces = 0
        while f'trace_time_{n_traces}' in data:
            n_traces += 1

    n_traces = min(n_traces, max_traces)
    times_list = []
    prots_list = []

    for i in range(n_traces):
        t_key = f'trace_time_{i}'
        p_key = f'trace_prot_{i}'
        if t_key not in data or p_key not in data:
            continue
        t = data[t_key]
        p = data[p_key]
        # filter to valid (non-zero time) entries
        valid = t > 0
        if valid.sum() < 10:
            continue
        times_list.append(t[valid])
        prots_list.append(p[valid].astype(float))

    return times_list, prots_list


def normalized_autocorrelation(x):
    """
    Compute normalized autocorrelation C(tau) = <delta_x(t) delta_x(t+tau)> / <delta_x^2>.

    Uses FFT-based computation for efficiency.
    Returns C(tau) for tau = 0, dt, 2*dt, ... up to N//2.
    """
    x = np.asarray(x, dtype=float)
    x = x - np.mean(x)
    var = np.var(x)
    if var < 1e-12:
        return np.zeros(len(x) // 2)

    n = len(x)
    # FFT-based autocorrelation
    fft_x = np.fft.fft(x, n=2 * n)
    acf = np.fft.ifft(fft_x * np.conj(fft_x)).real[:n]
    # Normalize by number of overlapping pairs and variance
    acf = acf / (np.arange(n, 0, -1) * var)

    return acf[:n // 2]


def extract_tau_c(acf, dt):
    """
    Extract autocorrelation time tau_c: time for C(tau) to first drop below 1/e.

    Parameters
    ----------
    acf : ndarray
        Normalized autocorrelation function (C[0] should be ~1.0).
    dt : float
        Time step between trace points (minutes).

    Returns
    -------
    tau_c : float
        Autocorrelation time in minutes. NaN if never crosses 1/e.
    """
    threshold = 1.0 / np.e
    below = np.where(acf < threshold)[0]
    if len(below) == 0:
        return np.nan
    idx = below[0]
    # Linear interpolation for sub-sample accuracy
    if idx == 0:
        return 0.0
    c_prev = acf[idx - 1]
    c_curr = acf[idx]
    frac = (c_prev - threshold) / (c_prev - c_curr) if c_prev != c_curr else 0.5
    tau_c = (idx - 1 + frac) * dt
    return tau_c


def integrated_autocorrelation_time(acf, dt):
    """
    Compute integrated autocorrelation time: tau_int = dt * sum(C(tau)) for tau >= 0.

    Uses truncation at first negative crossing to avoid noise accumulation.
    """
    # Truncate at first negative value
    neg = np.where(acf < 0)[0]
    if len(neg) > 0:
        acf_trunc = acf[:neg[0]]
    else:
        acf_trunc = acf

    if len(acf_trunc) == 0:
        return 0.0

    tau_int = dt * np.sum(acf_trunc)
    return tau_int


def run_autocorrelation():
    """
    Compute autocorrelation functions and tau_c for conditions A, B, C.

    Returns
    -------
    results : dict
        Per-condition autocorrelation results.
    df : pd.DataFrame
        Summary table saved to CSV.
    """
    os.makedirs(PHASE8_DIR, exist_ok=True)

    conditions = ['A', 'B', 'C']
    labels = {
        'A': 'Asymmetric (Mce3R)',
        'B': 'Symmetric control',
        'C': 'Single-site control',
    }

    results = {}
    rows = []
    all_acfs = {}

    for cond in conditions:
        times_list, prots_list = load_traces(cond)
        if times_list is None or len(times_list) == 0:
            print(f"  [WARN] No traces for Condition {cond}")
            continue

        dt = times_list[0][1] - times_list[0][0]  # time step (minutes)
        acfs = []
        tau_cs = []
        tau_ints = []

        for i, prot in enumerate(prots_list):
            acf = normalized_autocorrelation(prot)
            acfs.append(acf)
            tc = extract_tau_c(acf, dt)
            ti = integrated_autocorrelation_time(acf, dt)
            tau_cs.append(tc)
            tau_ints.append(ti)

        # Average ACF across cells
        min_len = min(len(a) for a in acfs)
        acf_matrix = np.array([a[:min_len] for a in acfs])
        mean_acf = np.mean(acf_matrix, axis=0)
        std_acf = np.std(acf_matrix, axis=0)

        # tau_c from mean ACF
        tau_c_mean = extract_tau_c(mean_acf, dt)
        tau_int_mean = integrated_autocorrelation_time(mean_acf, dt)

        # Per-cell statistics
        tau_cs_valid = [t for t in tau_cs if np.isfinite(t)]
        tau_ints_valid = [t for t in tau_ints if np.isfinite(t)]

        result = {
            'condition': cond,
            'label': labels[cond],
            'n_traces': len(prots_list),
            'dt_min': dt,
            'tau_c_mean_acf': tau_c_mean,
            'tau_c_median_cells': np.median(tau_cs_valid) if tau_cs_valid else np.nan,
            'tau_c_std_cells': np.std(tau_cs_valid) if tau_cs_valid else np.nan,
            'tau_int_mean_acf': tau_int_mean,
            'tau_int_median_cells': np.median(tau_ints_valid) if tau_ints_valid else np.nan,
            'mean_acf': mean_acf,
            'std_acf': std_acf,
            'lags_min': np.arange(min_len) * dt,
        }
        results[cond] = result
        all_acfs[cond] = mean_acf

        rows.append({
            'condition': cond,
            'label': labels[cond],
            'n_traces': len(prots_list),
            'dt_min': dt,
            'tau_c_mean_acf_min': tau_c_mean,
            'tau_c_mean_acf_hr': tau_c_mean / 60 if np.isfinite(tau_c_mean) else np.nan,
            'tau_c_median_cells_min': result['tau_c_median_cells'],
            'tau_c_std_cells_min': result['tau_c_std_cells'],
            'tau_int_mean_acf_min': tau_int_mean,
            'tau_int_median_cells_min': result['tau_int_median_cells'],
        })

        print(f"  Condition {cond} ({labels[cond]}):")
        print(f"    tau_c (mean ACF)    = {tau_c_mean:.1f} min ({tau_c_mean/60:.1f} hr)")
        print(f"    tau_c (median cell) = {result['tau_c_median_cells']:.1f} min")
        print(f"    tau_int (mean ACF)  = {tau_int_mean:.1f} min ({tau_int_mean/60:.1f} hr)")

    df = pd.DataFrame(rows)
    csv_path = os.path.join(PHASE8_DIR, 'autocorrelation.csv')
    df.to_csv(csv_path, index=False)
    print(f"  Saved: {csv_path}")

    # Save ACF arrays for figure generation
    acf_path = os.path.join(PHASE8_DIR, 'acf_data.npz')
    save_dict = {}
    for cond, res in results.items():
        save_dict[f'mean_acf_{cond}'] = res['mean_acf']
        save_dict[f'std_acf_{cond}'] = res['std_acf']
        save_dict[f'lags_{cond}'] = res['lags_min']
    np.savez(acf_path, **save_dict)

    return results, df


if __name__ == '__main__':
    print("=== Phase 8: Autocorrelation Analysis ===")
    results, df = run_autocorrelation()
    print("\nSummary:")
    print(df.to_string(index=False))
