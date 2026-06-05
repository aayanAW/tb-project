"""
phase5_thermodynamic/thermodynamic_main.py — Orchestrate all Phase 5 steps.

Sequence:
  1. Create results/phase5/ directory
  2. Load MEME PWM, calibrate energies → energy_parameters.json
  3. Compute repression curves for 4 architectures → repression_curves.csv
  4. Run MCMC → mcmc_posteriors.npz, mcmc_summary.json
  5. Classify top-20 operators → operator_classifications.csv
  6. (Optional) DNA shape
  7. Write phase5_summary.json

Usage:
  python -m phase5_thermodynamic.thermodynamic_main
"""

import sys
import os
import json
import time
import numpy as np
import pandas as pd
import warnings

sys.path.insert(0, '/Users/aayanalwani/tb project/mce3r_stochastic')

RESULTS_DIR = '/Users/aayanalwani/tb project/mce3r_stochastic/results/phase5'
PHASE2_SUMMARY = '/Users/aayanalwani/tb project/mce3r_stochastic/results/phase2/phase2_summary.json'


def run_phase5(
    results_dir=RESULTS_DIR,
    n_walkers=None,
    n_steps=None,
    n_burn=None,
    seed=None,
    top_n=20,
    run_mcmc_flag=True,
    verbose=True,
):
    """
    Run all Phase 5 steps.

    Parameters
    ----------
    results_dir : str — output directory
    n_walkers, n_steps, n_burn, seed : MCMC parameters
    top_n : int — number of FIMO sites to classify
    run_mcmc_flag : bool — set False to skip MCMC (for quick testing)
    verbose : bool

    Returns
    -------
    dict — standard summary contract
    """
    from config.parameters import PARAMS, MCMC_N_WALKERS, MCMC_N_STEPS, MCMC_N_BURN, MCMC_SEED

    if n_walkers is None:
        n_walkers = MCMC_N_WALKERS
    if n_steps is None:
        n_steps = MCMC_N_STEPS
    if n_burn is None:
        n_burn = MCMC_N_BURN
    if seed is None:
        seed = MCMC_SEED

    os.makedirs(results_dir, exist_ok=True)
    t_start = time.time()

    n_sanity_pass = 0
    n_sanity_fail = 0
    n_science_expected = 0
    n_science_unexpected = 0
    output_files = []

    def sanity(name, cond, val=""):
        nonlocal n_sanity_pass, n_sanity_fail
        if cond:
            if verbose:
                print(f"  SANITY PASS: {name} {val}")
            n_sanity_pass += 1
        else:
            if verbose:
                print(f"  SANITY FAIL: {name} {val}")
            n_sanity_fail += 1

    def scientific(name, cond, expected, val=""):
        nonlocal n_science_expected, n_science_unexpected
        status = "EXPECTED" if cond else "UNEXPECTED"
        if verbose:
            print(f"  SCIENCE {status}: {name} (expect: {expected}) {val}")
        if cond:
            n_science_expected += 1
        else:
            n_science_unexpected += 1

    # ================================================================
    # Step 1: Load MEME PWM and calibrate energies
    # ================================================================
    if verbose:
        print("\n[Phase 5] Step 1: Load MEME PWM and calibrate energies...")

    from phase5_thermodynamic.energy_calibration import (
        load_pwm_from_meme, pwm_to_energy_matrix, calibrate_energies, build_calibration
    )

    calibration = build_calibration()
    motif_width = calibration['motif_width']

    sanity("PWM loaded",
           calibration['energy_matrix'] is not None, "")
    sanity("Energy matrix no NaN/Inf",
           not np.any(np.isnan(calibration['energy_matrix'])) and
           not np.any(np.isinf(calibration['energy_matrix'])),
           f"— shape {calibration['energy_matrix'].shape}")
    sanity("Kd_strong reproduced",
           abs(calibration['Kd_strong_pred'] - 2.4) / 2.4 < 0.01,
           f"— {calibration['Kd_strong_pred']:.4f} nM")
    sanity("Kd_weak reproduced",
           abs(calibration['Kd_weak_pred'] - 49.0) / 49.0 < 0.01,
           f"— {calibration['Kd_weak_pred']:.4f} nM")

    # Save energy parameters
    energy_params = {
        'motif_width':      motif_width,
        'dG_strong':        calibration['dG_strong'],
        'dG_weak':          calibration['dG_weak'],
        'dG_strong_exp':    calibration['dG_strong_exp'],
        'dG_weak_exp':      calibration['dG_weak_exp'],
        'offset_strong':    calibration['offset_strong'],
        'offset_weak':      calibration['offset_weak'],
        'offset_mean':      calibration['offset_mean'],
        'Kd_strong_exp':    calibration['Kd_strong_exp'],
        'Kd_weak_exp':      calibration['Kd_weak_exp'],
        'Kd_strong_pred':   calibration['Kd_strong_pred'],
        'Kd_weak_pred':     calibration['Kd_weak_pred'],
        'kT':               calibration['kT'],
        'strong_seq':       calibration['strong_seq'],
        'weak_seq':         calibration['weak_seq'],
    }
    energy_params_path = os.path.join(results_dir, 'energy_parameters.json')
    with open(energy_params_path, 'w') as f:
        json.dump(energy_params, f, indent=2)
    output_files.append(energy_params_path)
    if verbose:
        print(f"  Saved: energy_parameters.json")
        print(f"    dG_strong = {calibration['dG_strong']:.3f} kcal/mol (Kd={calibration['Kd_strong_exp']} nM)")
        print(f"    dG_weak   = {calibration['dG_weak']:.3f} kcal/mol (Kd={calibration['Kd_weak_exp']} nM)")

    # ================================================================
    # Step 2: Compute repression curves for 4 architectures
    # ================================================================
    if verbose:
        print("\n[Phase 5] Step 2: Computing repression curves...")

    from phase5_thermodynamic.partition_function import (
        compare_architectures, repression_curve as pf_repression_curve
    )

    concentrations = np.logspace(-2, 4, 200)  # 0.01 to 10,000 nM
    arch_results = compare_architectures(concentrations, calibration, omega=1.0)

    # Build DataFrame
    rows = []
    for arch_name, curve in arch_results.items():
        for i in range(len(concentrations)):
            rows.append({
                'architecture':    arch_name,
                'concentration_nM': concentrations[i],
                'k_mean':           curve['k_mean'][i],
                'fold_repression':  curve['fold'][i],
                'prob_U':           curve['prob_U'][i],
                'prob_S':           curve['prob_S'][i],
                'prob_W':           curve['prob_W'][i],
                'prob_D':           curve['prob_D'][i],
            })
    df_curves = pd.DataFrame(rows)

    curves_path = os.path.join(results_dir, 'repression_curves.csv')
    df_curves.to_csv(curves_path, index=False)
    output_files.append(curves_path)
    if verbose:
        print(f"  Saved: repression_curves.csv ({len(df_curves)} rows)")

    # Science check: asymmetric should differ from symmetric
    asym_folds = arch_results['native_asymmetric']['fold']
    sym_folds  = arch_results['symmetric']['fold']
    scientific("Asymmetric differs from symmetric architecture",
               not np.allclose(asym_folds, sym_folds, atol=0.1),
               expected="True",
               val=f"— max diff={np.max(np.abs(asym_folds - sym_folds)):.3f}")

    # Science check: prob_U + prob_S + prob_W + prob_D = 1
    prob_sums = (arch_results['native_asymmetric']['prob_U'] +
                 arch_results['native_asymmetric']['prob_S'] +
                 arch_results['native_asymmetric']['prob_W'] +
                 arch_results['native_asymmetric']['prob_D'])
    sanity("Probabilities sum to 1 in repression curve",
           np.allclose(prob_sums, 1.0, atol=1e-10),
           f"— max deviation = {np.max(np.abs(prob_sums - 1.0)):.2e}")

    # ================================================================
    # Step 3: Run MCMC (Bayesian inference of omega and dG_spacer)
    # ================================================================
    if verbose:
        print("\n[Phase 5] Step 3: Running MCMC cooperativity inference...")

    from phase5_thermodynamic.cooperativity_inference import (
        run_mcmc, build_data_dict
    )

    mcmc_data = build_data_dict(
        phase2_summary_path=PHASE2_SUMMARY,
        calibration=calibration
    )

    if run_mcmc_flag:
        mcmc_result = run_mcmc(
            mcmc_data,
            n_walkers=n_walkers,
            n_steps=n_steps,
            n_burn=n_burn,
            seed=seed
        )

        omega_fitted   = mcmc_result['omega_median']
        dG_spacer_fitted = mcmc_result['dG_spacer_median']

        # Save posteriors
        posteriors_path = os.path.join(results_dir, 'mcmc_posteriors.npz')
        np.savez(
            posteriors_path,
            flat_chain=mcmc_result['flat_chain'],
            chain=mcmc_result['chain'],
            omega_samples=np.exp(mcmc_result['flat_chain'][:, 0]),
            dG_spacer_samples=mcmc_result['flat_chain'][:, 1],
        )
        output_files.append(posteriors_path)

        # Save MCMC summary
        mcmc_summary = {
            'omega_median':        float(omega_fitted),
            'omega_ci_16':         float(mcmc_result['omega_ci'][0]),
            'omega_ci_84':         float(mcmc_result['omega_ci'][1]),
            'ln_omega_median':     float(mcmc_result['ln_omega_median']),
            'ln_omega_ci_16':      float(mcmc_result['ln_omega_ci'][0]),
            'ln_omega_ci_84':      float(mcmc_result['ln_omega_ci'][1]),
            'dG_spacer_median':    float(dG_spacer_fitted),
            'dG_spacer_ci_16':     float(mcmc_result['dG_spacer_ci'][0]),
            'dG_spacer_ci_84':     float(mcmc_result['dG_spacer_ci'][1]),
            'acceptance_fraction': float(mcmc_result['acceptance_fraction']),
            'r_hat_ln_omega':      float(mcmc_result['r_hat'].get('ln_omega', np.nan)),
            'r_hat_dG_spacer':     float(mcmc_result['r_hat'].get('dG_spacer', np.nan)),
            'n_walkers':           n_walkers,
            'n_steps':             n_steps,
            'n_burn':              n_burn,
            'mce3r_conc_nM':       mcmc_data['mce3r_conc_nM'],
        }
        mcmc_summary_path = os.path.join(results_dir, 'mcmc_summary.json')
        with open(mcmc_summary_path, 'w') as f:
            json.dump(mcmc_summary, f, indent=2)
        output_files.append(mcmc_summary_path)

        # Sanity checks on MCMC
        sanity("Acceptance fraction in [0.05, 0.95]",
               0.05 <= mcmc_result['acceptance_fraction'] <= 0.95,
               f"— {mcmc_result['acceptance_fraction']:.3f}")
        sanity("omega_median > 0",
               omega_fitted > 0,
               f"— {omega_fitted:.3f}")
        rhat_ok = (mcmc_result['r_hat'].get('ln_omega', np.inf) < 1.2 and
                   mcmc_result['r_hat'].get('dG_spacer', np.inf) < 1.2)
        scientific("MCMC converged (R-hat < 1.2)",
                   rhat_ok,
                   expected=True,
                   val=f"— R-hat: ln_omega={mcmc_result['r_hat'].get('ln_omega',np.nan):.3f}, "
                       f"dG_spacer={mcmc_result['r_hat'].get('dG_spacer',np.nan):.3f}")

        if verbose:
            print(f"  Saved: mcmc_posteriors.npz, mcmc_summary.json")
    else:
        # Skip MCMC, use defaults
        omega_fitted     = 1.0
        dG_spacer_fitted = 0.0
        mcmc_result      = None
        if verbose:
            print("  MCMC skipped (run_mcmc_flag=False); using omega=1.0, dG_spacer=0.0")

    # ================================================================
    # Step 4: Classify top-20 operators
    # ================================================================
    if verbose:
        print("\n[Phase 5] Step 4: Classifying top-20 FIMO operators...")

    from phase5_thermodynamic.operator_classification import classify_all_operators

    df_cls = classify_all_operators(
        calibration=calibration,
        omega=omega_fitted,
        dG_spacer=dG_spacer_fitted,
        top_n=top_n
    )

    cls_path = os.path.join(results_dir, 'operator_classifications.csv')
    df_cls.to_csv(cls_path, index=False)
    output_files.append(cls_path)

    n_classified = len(df_cls)
    n_error = (df_cls['classification'] == 'error').sum()
    n_categories = df_cls[df_cls['classification'] != 'error']['classification'].nunique()

    sanity("All operators classified",
           n_classified == top_n,
           f"— {n_classified}/{top_n}")
    sanity("No classification errors",
           n_error == 0,
           f"— {n_error} errors")

    scientific("Multiple classification categories populated",
               n_categories >= 2,
               expected=True,
               val=f"— {n_categories} categories: {df_cls['classification'].value_counts().to_dict()}")

    if verbose:
        print(f"  Saved: operator_classifications.csv ({n_classified} sites)")
        print(f"  Classification distribution: {df_cls['classification'].value_counts().to_dict()}")

    # ================================================================
    # Step 5: (Optional) DNA shape
    # ================================================================
    if verbose:
        print("\n[Phase 5] Step 5: DNA shape features (optional)...")

    from phase5_thermodynamic.dna_shape import compute_dna_shape

    dna_shape_result = compute_dna_shape()
    dna_shape_status = dna_shape_result['status']
    if verbose:
        if dna_shape_status == 'ok':
            print(f"  DNA shape computed successfully")
        else:
            print(f"  DNA shape skipped: {dna_shape_result.get('reason', 'unknown')}")

    # ================================================================
    # Step 6: Write phase5_summary.json
    # ================================================================
    if verbose:
        print("\n[Phase 5] Step 6: Writing phase5_summary.json...")

    wall_time = time.time() - t_start

    # Scientific predictions validation
    # Prediction 1: P(S) + P(W) at K_half > 0.3 for asymmetric
    from phase5_thermodynamic.partition_function import state_probabilities, mu_from_concentration
    kT = calibration['kT']
    dG_s = calibration['dG_strong']
    dG_w = calibration['dG_weak']
    K_half_nM = df_cls.loc[df_cls['rank'] == 1, 'K_half_nM'].values
    K_half = float(K_half_nM[0]) if len(K_half_nM) > 0 else 10.0
    mu_Khalf = mu_from_concentration(K_half, kT)
    probs_Khalf = state_probabilities(dG_s, dG_w, mu_Khalf, omega=omega_fitted,
                                       dG_spacer=dG_spacer_fitted, kT=kT)
    P_intermediate = float(probs_Khalf[1] + probs_Khalf[2])

    scientific("Prediction 1: P(S)+P(W) > 0.3 at K_half for asymmetric",
               P_intermediate > 0.3,
               expected=True,
               val=f"— P_intermediate = {P_intermediate:.3f}")

    # Prediction 7: Top-20 operators have >=2 categories
    scientific("Prediction 7: >=2 classification categories",
               n_categories >= 2,
               expected=True,
               val=f"— {n_categories} categories")

    summary = {
        'phase': 5,
        'description': 'Mce3R thermodynamic model — PWM energy calibration + MCMC cooperativity inference',
        'status': 'complete' if n_sanity_fail == 0 else 'partial',
        'wall_time_sec': round(wall_time, 1),
        'n_sanity_pass': n_sanity_pass,
        'n_sanity_fail': n_sanity_fail,
        'n_science_expected': n_science_expected,
        'n_science_unexpected': n_science_unexpected,
        'output_files': [os.path.basename(f) for f in output_files],
        'energy_calibration': {
            'motif_width':   motif_width,
            'dG_strong':     calibration['dG_strong'],
            'dG_weak':       calibration['dG_weak'],
            'Kd_strong_exp': calibration['Kd_strong_exp'],
            'Kd_weak_exp':   calibration['Kd_weak_exp'],
            'kT':            calibration['kT'],
        },
        'mcmc': {
            'omega_median':     float(omega_fitted),
            'dG_spacer_median': float(dG_spacer_fitted),
            'skipped':          not run_mcmc_flag,
        } if run_mcmc_flag and mcmc_result is not None else {
            'omega_median': 1.0, 'dG_spacer_median': 0.0, 'skipped': True
        },
        'operator_classifications': {
            'n_classified':  int(n_classified),
            'n_categories':  int(n_categories),
            'distribution':  df_cls['classification'].value_counts().to_dict(),
        },
        'dna_shape': {
            'status': dna_shape_status,
            'reason': dna_shape_result.get('reason', ''),
        },
        'scientific_predictions': {
            'P_intermediate_at_Khalf': float(P_intermediate),
            'n_classification_categories': int(n_categories),
        },
    }

    summary_path = os.path.join(results_dir, 'phase5_summary.json')
    with open(summary_path, 'w') as f:
        json.dump(summary, f, indent=2)
    output_files.append(summary_path)

    if verbose:
        print(f"\n{'='*60}")
        print(f"Phase 5 complete in {wall_time:.1f}s")
        print(f"  Sanity: {n_sanity_pass} pass, {n_sanity_fail} fail")
        print(f"  Science: {n_science_expected} expected, {n_science_unexpected} unexpected")
        print(f"  Output files: {len(output_files)}")

    return summary


# ---------------------------------------------------------------------------
# Main entrypoint
# ---------------------------------------------------------------------------
if __name__ == '__main__':
    print("=" * 60)
    print("Phase 5: Thermodynamic Model")
    print("=" * 60)

    summary = run_phase5()

    print("\n--- Phase 5 Summary ---")
    print(f"Status:          {summary['status']}")
    print(f"Wall time:       {summary['wall_time_sec']}s")
    print(f"Sanity pass/fail: {summary['n_sanity_pass']}/{summary['n_sanity_fail']}")
    print(f"Science expected/unexpected: "
          f"{summary['n_science_expected']}/{summary['n_science_unexpected']}")
    print(f"omega (MCMC):    {summary['mcmc']['omega_median']:.3f}")
    print(f"dG_spacer (MCMC):{summary['mcmc']['dG_spacer_median']:.3f} kcal/mol")
    print(f"N classifications: {summary['operator_classifications']['n_classified']}")
    print(f"N categories:    {summary['operator_classifications']['n_categories']}")

    if summary['n_sanity_fail'] > 0:
        print(f"\nWARNING: {summary['n_sanity_fail']} sanity checks failed!")
        sys.exit(1)
    else:
        print("\nAll sanity checks passed.")
        sys.exit(0)
