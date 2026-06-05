"""
Phase 8: Power Spectral Density (PSD) analysis of protein expression traces.

Computes PSD to identify characteristic frequencies of expression fluctuations.
Low-frequency power indicates slow switching (long-lived states), while
high-frequency power indicates fast fluctuations (Poisson noise).

The asymmetric operator is predicted to have more low-frequency power
(slower switching between partially-bound states).
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


def compute_psd(x, dt):
    """
    Compute one-sided Power Spectral Density via FFT.

    Parameters
    ----------
    x : ndarray
        Time series (protein counts).
    dt : float
        Sampling interval (minutes).

    Returns
    -------
    freqs : ndarray
        Frequency array (cycles per minute).
    psd : ndarray
        Power spectral density (protein^2 / (cycles/min)).
    """
    x = np.asarray(x, dtype=float)
    x = x - np.mean(x)
    n = len(x)

    fft_vals = np.fft.rfft(x)
    psd = (2.0 * dt / n) * np.abs(fft_vals) ** 2
    freqs = np.fft.rfftfreq(n, d=dt)

    # Exclude DC component
    return freqs[1:], psd[1:]


def corner_frequency(freqs, psd):
    """
    Estimate corner frequency where PSD transitions from flat to decay.

    Uses the frequency where cumulative power reaches 50% of total.
    """
    cumulative = np.cumsum(psd)
    total = cumulative[-1]
    if total < 1e-12:
        return np.nan
    idx = np.searchsorted(cumulative, 0.5 * total)
    if idx >= len(freqs):
        return freqs[-1]
    return freqs[idx]


def low_freq_power_fraction(freqs, psd, cutoff_period_min=500):
    """
    Fraction of total power below cutoff frequency.

    cutoff_period_min=500 means frequencies < 1/500 cycles/min.
    Higher fraction = more slow switching.
    """
    cutoff_freq = 1.0 / cutoff_period_min
    low_mask = freqs <= cutoff_freq
    if psd.sum() < 1e-12:
        return np.nan
    return psd[low_mask].sum() / psd.sum()


def run_psd_analysis():
    """
    Compute PSD for conditions A, B, C and compare spectral properties.
    """
    os.makedirs(PHASE8_DIR, exist_ok=True)

    from phase8_temporal.autocorrelation import load_traces

    conditions = ['A', 'B', 'C']
    labels = {'A': 'Asymmetric', 'B': 'Symmetric', 'C': 'Single-site'}

    rows = []
    psd_data = {}

    for cond in conditions:
        times_list, prots_list = load_traces(cond)
        if times_list is None or len(times_list) == 0:
            print(f"  [WARN] No traces for Condition {cond}")
            continue

        dt = times_list[0][1] - times_list[0][0]
        all_psds = []
        freqs_ref = None

        for prot in prots_list:
            freqs, psd = compute_psd(prot, dt)
            all_psds.append(psd)
            if freqs_ref is None:
                freqs_ref = freqs

        # Average PSD across cells
        min_len = min(len(p) for p in all_psds)
        psd_matrix = np.array([p[:min_len] for p in all_psds])
        mean_psd = np.mean(psd_matrix, axis=0)
        freqs_ref = freqs_ref[:min_len]

        fc = corner_frequency(freqs_ref, mean_psd)
        lf_frac = low_freq_power_fraction(freqs_ref, mean_psd, cutoff_period_min=500)
        lf_frac_1000 = low_freq_power_fraction(freqs_ref, mean_psd, cutoff_period_min=1000)

        psd_data[cond] = {
            'freqs': freqs_ref,
            'mean_psd': mean_psd,
        }

        rows.append({
            'condition': cond,
            'label': labels[cond],
            'n_traces': len(prots_list),
            'corner_freq_per_min': fc,
            'corner_period_min': 1.0 / fc if fc > 0 and np.isfinite(fc) else np.nan,
            'low_freq_fraction_500min': lf_frac,
            'low_freq_fraction_1000min': lf_frac_1000,
            'total_power': mean_psd.sum(),
        })

        print(f"  Condition {cond} ({labels[cond]}):")
        print(f"    Corner frequency: {fc:.6f} /min (period: {1/fc:.0f} min)" if np.isfinite(fc) and fc > 0 else f"    Corner frequency: N/A")
        print(f"    Low-freq power (<1/500 min): {lf_frac:.3f}")
        print(f"    Low-freq power (<1/1000 min): {lf_frac_1000:.3f}")

    df = pd.DataFrame(rows)
    csv_path = os.path.join(PHASE8_DIR, 'power_spectrum.csv')
    df.to_csv(csv_path, index=False)
    print(f"  Saved: {csv_path}")

    # Save PSD arrays for figure generation
    save_dict = {}
    for cond, data in psd_data.items():
        save_dict[f'freqs_{cond}'] = data['freqs']
        save_dict[f'psd_{cond}'] = data['mean_psd']
    np.savez(os.path.join(PHASE8_DIR, 'psd_data.npz'), **save_dict)

    return df, psd_data


if __name__ == '__main__':
    print("=== Phase 8: Power Spectral Density Analysis ===")
    df, _ = run_psd_analysis()
    print("\nSummary:")
    print(df.to_string(index=False))
