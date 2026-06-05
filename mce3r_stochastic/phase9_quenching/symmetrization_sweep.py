"""
Phase 9: Symmetrization sweep — dose-response for noise quenching.

Simulates a gradual transition from asymmetric (wild-type) to symmetric operator
by interpolating Kd and block values. At each "dose" of symmetrization, measures
CV and persister fraction to build a therapeutic dose-response curve.

The key output is an IC50 of noise: the degree of symmetrization required to
reduce the excess noise (delta-CV) by 50%.
"""

import sys
import os
import time
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from config.parameters import PARAMS
from phase2_simulation.operator_model import OperatorModel
from phase2_simulation.gillespie_engine import run_population


PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PHASE9_DIR = os.path.join(PROJECT_ROOT, 'results', 'phase9')


def interpolate_parameters(frac_symmetric):
    """
    Interpolate between wild-type asymmetric and symmetric operator.

    frac_symmetric = 0.0: fully asymmetric (wild-type)
    frac_symmetric = 1.0: fully symmetric (geometric mean)

    Interpolation is done in log-space for Kd values (geometric interpolation)
    and linear for block fractions.

    Parameters
    ----------
    frac_symmetric : float
        Degree of symmetrization [0, 1].

    Returns
    -------
    dict with Kd_strong, Kd_weak, block_strong, block_weak
    """
    # Log-space interpolation for Kd
    log_Kd_s_wt = np.log(PARAMS['Kd_strong'])
    log_Kd_w_wt = np.log(PARAMS['Kd_weak'])
    log_Kd_sym = np.log(PARAMS['Kd_symmetric'])

    log_Kd_s = log_Kd_s_wt + frac_symmetric * (log_Kd_sym - log_Kd_s_wt)
    log_Kd_w = log_Kd_w_wt + frac_symmetric * (log_Kd_sym - log_Kd_w_wt)

    Kd_s = np.exp(log_Kd_s)
    Kd_w = np.exp(log_Kd_w)

    # Linear interpolation for block fractions
    block_s = PARAMS['block_strong'] + frac_symmetric * (PARAMS['block_symmetric'] - PARAMS['block_strong'])
    block_w = PARAMS['block_weak'] + frac_symmetric * (PARAMS['block_symmetric'] - PARAMS['block_weak'])

    return {
        'Kd_strong': Kd_s,
        'Kd_weak': Kd_w,
        'block_strong': block_s,
        'block_weak': block_w,
    }


def compute_stats(proteins):
    """Compute noise metrics from protein distribution."""
    proteins = proteins.astype(float)
    mean = np.mean(proteins)
    std = np.std(proteins, ddof=1)
    cv = std / mean if mean > 0 else np.nan
    fano = np.var(proteins, ddof=1) / mean if mean > 0 else np.nan
    return {'mean': mean, 'std': std, 'cv': cv, 'fano': fano}


def persister_fraction(proteins, threshold):
    """Fraction of cells below persistence threshold."""
    return float(np.sum(proteins < threshold) / len(proteins))


def run_symmetrization_sweep(n_doses=20, n_cells=5000, seed=55555):
    """
    Sweep from fully asymmetric (dose=0) to fully symmetric (dose=1).

    Parameters
    ----------
    n_doses : int
        Number of symmetrization steps.
    n_cells : int
        Cells per dose point.
    seed : int
        Base random seed.

    Returns
    -------
    df : pd.DataFrame
        Dose-response table.
    ic50 : float
        IC50 of noise (frac_symmetric where delta-CV drops to 50%).
    """
    os.makedirs(PHASE9_DIR, exist_ok=True)

    # Define threshold: use calibrated value from Phase 6 v2 (absolute_calibrated)
    # or fall back to a biologically meaningful value based on Condition B baseline
    cond_b_path = os.path.join(PROJECT_ROOT, 'results', 'phase2', 'condition_B.npz')
    if os.path.exists(cond_b_path):
        b = np.load(cond_b_path)
        prot_b = b['proteins'].astype(float)
        # Threshold: mean_B - 2*std_B (low-expression tail of symmetric baseline)
        threshold = float(np.mean(prot_b) - 2 * np.std(prot_b))
        threshold = max(threshold, float(np.percentile(prot_b, 0.1)))  # floor at 0.1th percentile
    else:
        threshold = 112.0  # calibrated value from Phase 6

    doses = np.linspace(0.0, 1.0, n_doses)
    rows = []

    # Get reference CV values
    cv_asym = None  # will be dose=0
    cv_sym = None   # will be dose=1

    print(f"  Running {n_doses} dose points x {n_cells} cells...")
    start = time.time()

    for i, frac in enumerate(doses):
        params = interpolate_parameters(frac)

        model = OperatorModel(
            Kd_strong=params['Kd_strong'],
            Kd_weak=params['Kd_weak'],
            block_strong=params['block_strong'],
            block_weak=params['block_weak'],
        )
        model_arrays = model.get_numba_arrays()

        result = run_population(model_arrays, n_cells, seed + i * 1000)
        proteins = result['proteins']

        stats = compute_stats(proteins)
        pf = persister_fraction(proteins, threshold)

        if i == 0:
            cv_asym = stats['cv']
        if i == n_doses - 1:
            cv_sym = stats['cv']

        # Kd ratio at this dose
        kd_ratio = params['Kd_weak'] / params['Kd_strong']

        rows.append({
            'dose_index': i,
            'frac_symmetric': frac,
            'Kd_strong_nM': params['Kd_strong'],
            'Kd_weak_nM': params['Kd_weak'],
            'Kd_ratio': kd_ratio,
            'block_strong': params['block_strong'],
            'block_weak': params['block_weak'],
            'mean_protein': stats['mean'],
            'std_protein': stats['std'],
            'cv': stats['cv'],
            'fano': stats['fano'],
            'persister_fraction': pf,
            'n_cells': n_cells,
            'threshold': threshold,
        })

        if (i + 1) % 5 == 0 or i == 0:
            print(f"    Dose {i+1}/{n_doses}: frac={frac:.2f}, "
                  f"ratio={kd_ratio:.1f}, CV={stats['cv']:.4f}, "
                  f"persist={pf:.4f}")

    elapsed = time.time() - start
    print(f"  Sweep complete in {elapsed:.1f}s")

    df = pd.DataFrame(rows)

    # Compute delta-CV relative to symmetric
    delta_cv_max = cv_asym - cv_sym if cv_asym and cv_sym else 0
    df['delta_cv'] = df['cv'] - cv_sym if cv_sym else df['cv']
    df['noise_reduction_frac'] = 1.0 - (df['delta_cv'] / delta_cv_max) if delta_cv_max > 0 else 0.0

    # Find IC50: frac_symmetric where noise_reduction_frac = 0.5
    ic50 = np.nan
    if delta_cv_max > 0:
        nrf = df['noise_reduction_frac'].values
        doses_arr = df['frac_symmetric'].values
        crossings = np.where(np.diff(np.sign(nrf - 0.5)))[0]
        if len(crossings) > 0:
            idx = crossings[0]
            # Linear interpolation
            f1, f2 = doses_arr[idx], doses_arr[idx + 1]
            n1, n2 = nrf[idx], nrf[idx + 1]
            ic50 = f1 + (0.5 - n1) * (f2 - f1) / (n2 - n1) if n2 != n1 else f1

    csv_path = os.path.join(PHASE9_DIR, 'symmetrization_sweep.csv')
    df.to_csv(csv_path, index=False)
    print(f"  Saved: {csv_path}")

    print(f"\n  Results:")
    print(f"    CV(asymmetric, dose=0): {cv_asym:.4f}")
    print(f"    CV(symmetric, dose=1):  {cv_sym:.4f}")
    print(f"    Delta-CV:               {delta_cv_max:.4f}")
    print(f"    IC50 of noise:          {ic50:.3f}" if np.isfinite(ic50) else "    IC50: not found")

    return df, ic50, {'cv_asym': cv_asym, 'cv_sym': cv_sym, 'delta_cv_max': delta_cv_max}


if __name__ == '__main__':
    print("=== Phase 9: Symmetrization Sweep ===")
    df, ic50, info = run_symmetrization_sweep(n_doses=20, n_cells=5000)
