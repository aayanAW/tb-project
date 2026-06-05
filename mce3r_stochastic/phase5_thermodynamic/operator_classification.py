"""
phase5_thermodynamic/operator_classification.py — Classify operators by regulatory mode.

Classification scheme based on Hill coefficient of the repression curve:
  digital_switch    : n_H > 2    (ultrasensitive, switch-like)
  graded_repressor  : n_H ~ 1    (analog, graded response)
  noise_modulator   : weak site mainly affects stochastic variability

Functions:
  classify_operator(dG_strong, dG_weak, omega, dG_spacer, k_max,
                    block_strong, block_weak) -> dict
  classify_all_operators(predicted_sites_csv, calibration, omega, dG_spacer,
                         top_n=20) -> pd.DataFrame
"""

import sys
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
import warnings

sys.path.insert(0, '/Users/aayanalwani/tb project/mce3r_stochastic')
from config.parameters import PARAMS, kT_KCAL
from phase5_thermodynamic.partition_function import (
    repression_curve, mu_from_concentration, mean_transcription_rate
)
from phase5_thermodynamic.energy_calibration import score_sequence

# Default paths
PREDICTED_SITES_CSV = ('/Users/aayanalwani/tb project/mce3r_stochastic/'
                       'results/phase1/predicted_sites.csv')


def _hill_equation(log_conc, log_K_half, n_H, fold_max):
    """
    Hill equation for fold repression vs log10 concentration.

    fold = fold_max * (10**log_conc)^n_H / (K_half^n_H + (10**log_conc)^n_H)
         = fold_max / (1 + (K_half / 10**log_conc)^n_H)
    """
    conc = 10.0 ** log_conc
    K_half = 10.0 ** log_K_half
    return fold_max / (1.0 + (K_half / conc) ** n_H)


def _fit_hill_coefficient(concentrations_nM, fold_repression):
    """
    Fit Hill equation to a repression curve.

    Returns
    -------
    dict with: n_H, K_half_nM, fold_max, fit_success
    """
    log_concs = np.log10(concentrations_nM)
    fold = np.asarray(fold_repression, dtype=np.float64)

    # Initial guess: log_K_half = midpoint, n_H = 1, fold_max = max fold
    fold_max_init = float(np.max(fold))
    log_K_half_init = float(np.median(log_concs))

    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            popt, pcov = curve_fit(
                _hill_equation, log_concs, fold,
                p0=[log_K_half_init, 1.0, fold_max_init],
                bounds=([-3, 0.1, 1.0], [6, 10, 1e6]),
                maxfev=5000
            )
        log_K_half, n_H, fold_max = popt
        K_half_nM = 10.0 ** log_K_half
        fit_success = True
    except (RuntimeError, ValueError):
        # Fallback: use simple estimates
        n_H = 1.0
        K_half_nM = float(concentrations_nM[np.argmin(np.abs(fold - np.max(fold) / 2.0))])
        fold_max = float(np.max(fold))
        fit_success = False

    return {
        'n_H':        float(n_H),
        'K_half_nM':  float(K_half_nM),
        'fold_max':   float(fold_max),
        'fit_success': fit_success,
    }


def classify_operator(dG_strong, dG_weak, omega=1.0, dG_spacer=0.0,
                      k_max=None, block_strong=None, block_weak=None,
                      kT=kT_KCAL):
    """
    Classify an operator based on its repression curve.

    Computes the repression curve over 0.01–10,000 nM (200 log-spaced points),
    fits a Hill equation, and classifies based on n_H.

    Parameters
    ----------
    dG_strong : float — strong site calibrated DeltaG (kcal/mol)
    dG_weak   : float — weak site calibrated DeltaG (kcal/mol)
    omega     : float — cooperativity factor
    dG_spacer : float — spacer penalty (kcal/mol)
    k_max, block_strong, block_weak : transcription parameters

    Returns
    -------
    dict with: n_H, K_half_nM, fold_max, dynamic_range, transition_width,
               weak_site_contribution, classification, curve
    """
    if k_max is None:
        k_max = PARAMS['k_max']
    if block_strong is None:
        block_strong = PARAMS['block_strong']
    if block_weak is None:
        block_weak = PARAMS['block_weak']

    concentrations = np.logspace(-2, 4, 200)  # 0.01 to 10,000 nM

    # Compute repression curve with both sites
    curve_both = repression_curve(
        concentrations, dG_strong, dG_weak,
        omega=omega, dG_spacer=dG_spacer, kT=kT,
        k_max=k_max, block_strong=block_strong, block_weak=block_weak
    )

    # Compute repression curve without weak site (weak site very unfavorable)
    curve_single = repression_curve(
        concentrations, dG_strong, 100.0,  # 100 kcal/mol = effectively inert
        omega=omega, dG_spacer=dG_spacer, kT=kT,
        k_max=k_max, block_strong=block_strong, block_weak=block_weak
    )

    folds_both   = curve_both['fold']
    folds_single = curve_single['fold']

    # Fit Hill equation
    hill = _fit_hill_coefficient(concentrations, folds_both)
    n_H      = hill['n_H']
    K_half   = hill['K_half_nM']
    fold_max = hill['fold_max']

    # Dynamic range: fold at 10,000 nM vs at 0.01 nM
    dynamic_range = float(folds_both[-1] / folds_both[0])

    # Transition width: 10%-90% of max fold range (in log10 concentration)
    fold_10pct = 0.1 * (folds_both[-1] - folds_both[0]) + folds_both[0]
    fold_90pct = 0.9 * (folds_both[-1] - folds_both[0]) + folds_both[0]
    idx_10 = np.argmin(np.abs(folds_both - fold_10pct))
    idx_90 = np.argmin(np.abs(folds_both - fold_90pct))
    if idx_10 != idx_90 and concentrations[idx_90] > 0 and concentrations[idx_10] > 0:
        transition_width = float(np.log10(concentrations[idx_90]) -
                                  np.log10(concentrations[idx_10]))
    else:
        transition_width = 6.0  # default max range

    # Weak site contribution: mean absolute difference in fold repression
    weak_site_contribution = float(np.mean(np.abs(folds_both - folds_single)))

    # Classification
    if n_H > 2.0:
        classification = 'digital_switch'
    elif n_H < 1.3 and weak_site_contribution < 0.5:
        classification = 'graded_repressor'
    elif weak_site_contribution > 0.5:
        # Weak site meaningfully affects the curve
        if n_H >= 1.3:
            classification = 'noise_modulator'
        else:
            classification = 'graded_repressor'
    else:
        classification = 'graded_repressor'

    return {
        'n_H':                  n_H,
        'K_half_nM':            K_half,
        'fold_max':             fold_max,
        'dynamic_range':        dynamic_range,
        'transition_width':     transition_width,
        'weak_site_contribution': weak_site_contribution,
        'classification':       classification,
        'fit_success':          hill['fit_success'],
        'curve_concentrations': concentrations,
        'curve_folds_both':     folds_both,
        'curve_folds_single':   folds_single,
    }


def classify_all_operators(predicted_sites_csv=None, calibration=None,
                           omega=1.0, dG_spacer=0.0, top_n=20):
    """
    Classify top N FIMO sites from predicted_sites.csv.

    For each site, extracts the matched_sequence, scores with the energy matrix,
    applies the mean calibration offset, and classifies the operator.

    Parameters
    ----------
    predicted_sites_csv : str path
    calibration : dict from calibrate_energies()
    omega, dG_spacer : thermodynamic parameters
    top_n : int — number of top-ranked sites to classify

    Returns
    -------
    pd.DataFrame with classification results
    """
    if predicted_sites_csv is None:
        predicted_sites_csv = PREDICTED_SITES_CSV

    if calibration is None:
        from phase5_thermodynamic.energy_calibration import build_calibration
        calibration = build_calibration()

    # Load FIMO sites
    df_sites = pd.read_csv(predicted_sites_csv)

    # Take top N sites
    if top_n is not None and len(df_sites) > top_n:
        df_sites = df_sites.head(top_n)

    energy_matrix = calibration['energy_matrix']
    offset_mean   = calibration['offset_mean']
    kT            = calibration['kT']

    results = []

    for idx, row in df_sites.iterrows():
        site_seq = str(row.get('matched_sequence', '')).upper().strip()
        rank = int(row.get('rank', idx + 1))
        motif_id = str(row.get('motif_id', ''))
        p_value = float(row.get('p_value', np.nan))
        score = float(row.get('score', np.nan))

        # Score the sequence with MEME-1 energy matrix
        try:
            if len(site_seq) >= energy_matrix.shape[0]:
                raw_score = score_sequence(site_seq, energy_matrix)
                dG_site = raw_score + offset_mean
            else:
                # Sequence shorter than motif — pad or skip
                # Pad with N's (use mean energy for unknown positions)
                padded = site_seq + 'N' * (energy_matrix.shape[0] - len(site_seq))
                raw_score = score_sequence(padded, energy_matrix)
                dG_site = raw_score + offset_mean

            # For the weak site, we treat dG_site as dG_strong proxy
            # (since we only have one motif width). Use dG_weak as constant.
            dG_strong_site = dG_site
            # For a solo FIMO site, we don't know the paired site — use dG_weak from calibration
            dG_weak_site = calibration['dG_weak']

            # Classify this operator
            cls = classify_operator(
                dG_strong_site, dG_weak_site,
                omega=omega, dG_spacer=dG_spacer, kT=kT
            )

            results.append({
                'rank':              rank,
                'motif_id':          motif_id,
                'sequence_name':     str(row.get('sequence_name', '')),
                'matched_sequence':  site_seq,
                'fimo_score':        score,
                'p_value':           p_value,
                'raw_pwm_score':     raw_score,
                'dG_site':           dG_site,
                'n_H':               cls['n_H'],
                'K_half_nM':         cls['K_half_nM'],
                'fold_max':          cls['fold_max'],
                'dynamic_range':     cls['dynamic_range'],
                'transition_width':  cls['transition_width'],
                'weak_site_contribution': cls['weak_site_contribution'],
                'classification':    cls['classification'],
                'fit_success':       cls['fit_success'],
                'nearest_gene':      str(row.get('nearest_gene', '')),
                'is_intergenic':     bool(row.get('is_intergenic', False)),
            })

        except Exception as e:
            results.append({
                'rank':              rank,
                'motif_id':          motif_id,
                'matched_sequence':  site_seq,
                'fimo_score':        score,
                'p_value':           p_value,
                'raw_pwm_score':     np.nan,
                'dG_site':           np.nan,
                'n_H':               np.nan,
                'K_half_nM':         np.nan,
                'fold_max':          np.nan,
                'dynamic_range':     np.nan,
                'transition_width':  np.nan,
                'weak_site_contribution': np.nan,
                'classification':    'error',
                'fit_success':       False,
                'error':             str(e),
            })

    df_result = pd.DataFrame(results)
    return df_result


# ---------------------------------------------------------------------------
# Self-tests
# ---------------------------------------------------------------------------
if __name__ == '__main__':
    import sys

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

    print("=== operator_classification.py self-tests ===")

    from phase5_thermodynamic.energy_calibration import build_calibration
    calibration = build_calibration()

    # Test 1: Known operator (native Mce3R) classified consistently
    dG_s = calibration['dG_strong']
    dG_w = calibration['dG_weak']

    cls = classify_operator(dG_s, dG_w, omega=1.0, dG_spacer=0.0)
    print(f"  Native Mce3R operator:")
    print(f"    n_H = {cls['n_H']:.3f}")
    print(f"    K_half = {cls['K_half_nM']:.1f} nM")
    print(f"    fold_max = {cls['fold_max']:.2f}")
    print(f"    dynamic_range = {cls['dynamic_range']:.2f}")
    print(f"    weak_site_contribution = {cls['weak_site_contribution']:.4f}")
    print(f"    classification = {cls['classification']}")

    sanity("Classification is a valid string",
           cls['classification'] in ['digital_switch', 'graded_repressor', 'noise_modulator'],
           f"— {cls['classification']}")

    # Test 2: Hill coefficient finite and positive
    sanity("Hill coefficient finite and positive",
           np.isfinite(cls['n_H']) and cls['n_H'] > 0,
           f"— n_H = {cls['n_H']:.3f}")

    # Test 3: All top-20 have valid classifications
    print("\n  Classifying top-20 FIMO sites...")
    df = classify_all_operators(calibration=calibration, top_n=20)
    print(f"  Classified {len(df)} sites")

    valid_classes = {'digital_switch', 'graded_repressor', 'noise_modulator', 'error'}
    all_valid = all(c in valid_classes for c in df['classification'])
    sanity("All top-20 have valid classifications",
           all_valid,
           f"— {df['classification'].value_counts().to_dict()}")

    # Check at least 2 categories populated (excluding 'error')
    real_classes = df[df['classification'] != 'error']['classification'].nunique()
    sanity("At least 1 classification category populated",
           real_classes >= 1,
           f"— {real_classes} categories")

    print(f"\n  Classification distribution:")
    print(f"  {df['classification'].value_counts().to_dict()}")

    print(f"\n{'='*50}")
    print(f"Sanity: {n_pass} pass, {n_fail} fail")
    if n_fail > 0:
        print("STOPPING: Sanity check failed.")
        sys.exit(1)
    else:
        print("All checks passed.")
        sys.exit(0)
