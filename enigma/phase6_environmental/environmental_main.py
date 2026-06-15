"""
phase6_environmental/environmental_main.py — Phase 6 orchestrator.

Runs the full Phase 6 pipeline:
  1. Create results/phase6/
  2. Load Phase 5 omega from results/phase5/mcmc_summary.json
  3. Load Condition D proteins for threshold definition
  4. Run all 3 arch × 4 env × 2 model = 24 environmental conditions
  5. Compute mutual information (MI) by architecture × concentration
  6. Define persistence threshold + compute fractions for all conditions
  7. Generate persistence phase diagram (ratios × environments)
  8. Write all outputs to results/phase6/
  9. Write phase6_summary.json

Scientific predictions verified:
  1. CV(asymmetric) > CV(symmetric) under EVERY environment
  2. host_like produces higher mean expression than baseline (derepression)
  3. Pearson correlation between mce3r_protein and target_protein is NEGATIVE (autoregulation)
  4. Asymmetric has higher persister fraction than symmetric under host_like
  5. MI(asymmetric) >= MI(symmetric) at physiological Mce3R concentrations
"""

import sys
import os
import json
import time
import numpy as np
import pandas as pd

sys.path.insert(0, '/Users/aayanalwani/tb project/mce3r_stochastic')
from config.parameters import PARAMS
from phase6_environmental.run_environmental_conditions import (
    run_all_environmental, ARCHITECTURES, ENVIRONMENTS, OUTPUT_DIR
)
from phase6_environmental.mutual_information import (
    estimate_mutual_information, mutual_information_vs_concentration
)
from phase6_environmental.persistence_threshold import (
    define_threshold, compute_persister_fractions, persistence_phase_diagram
)


PHASE5_SUMMARY = '/Users/aayanalwani/tb project/mce3r_stochastic/results/phase5/mcmc_summary.json'
PHASE2_COND_D  = '/Users/aayanalwani/tb project/mce3r_stochastic/results/phase2/condition_D.npz'
PHASE2_COND_A  = '/Users/aayanalwani/tb project/mce3r_stochastic/results/phase2/condition_A.npz'


def run_phase6(
    n_cells: int = 10000,
    master_seed: int = 12345,
    n_cells_mi: int = 5000,
    n_cells_phase_diagram: int = 2000
) -> dict:
    """
    Orchestrate Phase 6 pipeline.

    Parameters
    ----------
    n_cells : int
        Cells per condition for main 24-condition run.
    master_seed : int
        Master random seed.
    n_cells_mi : int
        Cells per (concentration, environment) for MI sweep.
    n_cells_phase_diagram : int
        Cells per (ratio, environment) for phase diagram.

    Returns
    -------
    dict — standard summary dict:
        status, n_sanity_pass, n_sanity_fail, n_science_expected,
        n_science_unexpected, n_warn, mock_used, wall_time_sec, output_files
    """
    wall_start = time.time()

    n_sanity_pass    = 0
    n_sanity_fail    = 0
    n_science_expected   = 0
    n_science_unexpected = 0
    n_warn = 0
    output_files = []

    os.makedirs(OUTPUT_DIR, exist_ok=True)

    # =========================================================================
    # Step 1: Load Phase 5 outputs
    # =========================================================================
    print("\n=== Phase 6: Environmental Stochastic Extensions ===")
    print("\nStep 1: Loading Phase 5 outputs...")

    with open(PHASE5_SUMMARY, 'r') as f:
        phase5 = json.load(f)
    omega = phase5.get('omega_median', 1.083)
    dG_spacer = phase5.get('dG_spacer_median', -0.148)
    mce3r_conc_nM = phase5.get('mce3r_conc_nM', 332.0)
    print(f"  omega_median = {omega:.4f}")
    print(f"  dG_spacer_median = {dG_spacer:.4f} kcal/mol")
    print(f"  mce3r_conc_nM = {mce3r_conc_nM:.1f} nM")

    # =========================================================================
    # Step 2: Load Condition D proteins for threshold definition
    # =========================================================================
    print("\nStep 2: Loading Condition D distribution for persistence threshold...")
    d_data = np.load(PHASE2_COND_D)
    proteins_D = d_data['proteins']

    threshold_10 = define_threshold(proteins_D, method='percentile', percentile=10)
    threshold_5  = define_threshold(proteins_D, method='percentile', percentile=5)
    threshold_15 = define_threshold(proteins_D, method='percentile', percentile=15)
    print(f"  Threshold (5th pct) = {threshold_5:.0f}")
    print(f"  Threshold (10th pct) = {threshold_10:.0f}")
    print(f"  Threshold (15th pct) = {threshold_15:.0f}")

    # Sanity: Condition D persister fraction ≈ 10%
    frac_D_10 = float(np.mean(proteins_D < threshold_10))
    ok = abs(frac_D_10 - 0.10) < 0.02
    if ok:
        n_sanity_pass += 1
        print(f"  SANITY PASS: Condition D persister fraction ≈ 10% "
              f"(got {frac_D_10:.3f})")
    else:
        n_sanity_fail += 1
        print(f"  SANITY FAIL: Condition D persister fraction = {frac_D_10:.3f}, "
              f"expected ~0.10")

    # =========================================================================
    # Step 3: Run all 24 environmental conditions
    # =========================================================================
    print(f"\nStep 3: Running 24 conditions ({n_cells} cells each)...")
    results = run_all_environmental(n_cells=n_cells, master_seed=master_seed)

    # Verify all 24 files written
    all_files_ok = True
    for arch in ARCHITECTURES:
        for env in ENVIRONMENTS:
            for mt in ['single', 'two_species']:
                fpath = os.path.join(OUTPUT_DIR,
                                     f'env_condition_{arch}_{env}_{mt}.npz')
                if os.path.exists(fpath):
                    output_files.append(fpath)
                else:
                    all_files_ok = False
                    print(f"  MISSING: {fpath}")

    two_sp_path = os.path.join(OUTPUT_DIR, 'two_species_results.npz')
    if os.path.exists(two_sp_path):
        output_files.append(two_sp_path)

    if all_files_ok:
        n_sanity_pass += 1
        print("  SANITY PASS: All 24 NPZ files exist")
    else:
        n_sanity_fail += 1
        print("  SANITY FAIL: Some NPZ files missing")

    # =========================================================================
    # Step 4: Scientific prediction 1 — CV(asymmetric) > CV(symmetric)
    # =========================================================================
    print("\nStep 4: Checking CV predictions...")
    cv_fails = []
    for env in ENVIRONMENTS:
        stats_asym = results[('asymmetric', env, 'single')]['stats']
        stats_sym  = results[('symmetric',  env, 'single')]['stats']
        if stats_asym['cv'] > stats_sym['cv']:
            n_science_expected += 1
            print(f"  SCIENCE PASS: CV(asym) > CV(sym) under {env} "
                  f"({stats_asym['cv']:.4f} > {stats_sym['cv']:.4f})")
        else:
            n_science_unexpected += 1
            cv_fails.append(env)
            print(f"  SCIENCE UNEXPECTED: CV(asym)={stats_asym['cv']:.4f} NOT > "
                  f"CV(sym)={stats_sym['cv']:.4f} under {env}")

    # =========================================================================
    # Step 5: Scientific prediction 2 — host_like mean > baseline (derepression)
    # =========================================================================
    print("\nStep 5: Checking derepression predictions...")
    for arch in ARCHITECTURES:
        mean_base = results[(arch, 'baseline',  'single')]['stats']['mean']
        mean_host = results[(arch, 'host_like', 'single')]['stats']['mean']
        if mean_host > mean_base:
            n_science_expected += 1
            print(f"  SCIENCE PASS: host_like mean > baseline for {arch} "
                  f"({mean_host:.1f} > {mean_base:.1f})")
        else:
            n_science_unexpected += 1
            print(f"  SCIENCE UNEXPECTED: host_like mean ({mean_host:.1f}) NOT > "
                  f"baseline ({mean_base:.1f}) for {arch}")

    # =========================================================================
    # Step 6: Scientific prediction 3 — Pearson correlation (mce3r vs target) < 0
    # =========================================================================
    print("\nStep 6: Checking mce3r–target correlation (autoregulation)...")
    pearson_results = {}
    for arch in ['asymmetric', 'symmetric']:
        for env in ['baseline', 'host_like']:
            res_two = results[(arch, env, 'two_species')]
            mce3r_p = res_two['mce3r_proteins'].astype(np.float64)
            target_p = res_two['target_proteins'].astype(np.float64)
            if np.std(mce3r_p) > 0 and np.std(target_p) > 0:
                r = float(np.corrcoef(mce3r_p, target_p)[0, 1])
            else:
                r = 0.0
            pearson_results[(arch, env)] = r
            if r < 0:
                n_science_expected += 1
                print(f"  SCIENCE PASS: Pearson(mce3r, target) = {r:.4f} < 0 "
                      f"[{arch}, {env}] (autoregulation)")
            else:
                n_science_unexpected += 1
                print(f"  SCIENCE UNEXPECTED: Pearson(mce3r, target) = {r:.4f} >= 0 "
                      f"[{arch}, {env}] — autoregulation not detected")

    # =========================================================================
    # Step 7: Compute persister fractions for all conditions
    # =========================================================================
    print("\nStep 7: Computing persister fractions...")

    # Build results dict in the format compute_persister_fractions expects
    pf_results = {}
    for key, data in results.items():
        pf_results[key] = data

    df_pf = compute_persister_fractions(pf_results, threshold=threshold_10)

    # Scientific prediction 4: asymmetric > symmetric persister fraction under host_like
    asym_hl = df_pf[
        (df_pf['architecture'] == 'asymmetric') &
        (df_pf['environment']  == 'host_like') &
        (df_pf['model_type']   == 'single')
    ]
    sym_hl = df_pf[
        (df_pf['architecture'] == 'symmetric') &
        (df_pf['environment']  == 'host_like') &
        (df_pf['model_type']   == 'single')
    ]
    if len(asym_hl) > 0 and len(sym_hl) > 0:
        f_asym = asym_hl['persister_fraction'].values[0]
        f_sym  = sym_hl['persister_fraction'].values[0]
        if f_asym >= f_sym:
            n_science_expected += 1
            print(f"  SCIENCE PASS: Asymmetric persister fraction >= symmetric under host_like "
                  f"({f_asym:.4f} >= {f_sym:.4f})")
        else:
            n_science_unexpected += 1
            print(f"  SCIENCE UNEXPECTED: Asymmetric persister fraction ({f_asym:.4f}) "
                  f"< symmetric ({f_sym:.4f}) under host_like")

    # Save persister fractions CSV
    pf_path = os.path.join(OUTPUT_DIR, 'persistence_fractions.csv')
    df_pf.to_csv(pf_path, index=False)
    output_files.append(pf_path)
    print(f"  Saved {pf_path} ({len(df_pf)} rows)")

    # =========================================================================
    # Step 8: Compute mutual information vs concentration
    # =========================================================================
    print("\nStep 8: Computing mutual information vs concentration...")

    # Use 10 concentrations spanning physiological range
    concentrations_nM = np.array([10, 30, 60, 100, 150, 200, 332, 500, 800, 1200],
                                  dtype=np.float64)

    mi_dfs = []
    for arch in ['asymmetric', 'symmetric', 'single_site']:
        print(f"  Computing MI for {arch} ({len(concentrations_nM)} concentrations)...")
        df_mi = mutual_information_vs_concentration(
            architecture=arch,
            concentrations_nM=concentrations_nM,
            environments=ENVIRONMENTS,
            n_cells=n_cells_mi,
            seed=master_seed + 200
        )
        mi_dfs.append(df_mi)

    df_mi_all = pd.concat(mi_dfs, ignore_index=True)

    # Scientific prediction 5: MI(asymmetric) >= MI(symmetric) at ~332 nM
    mi_asym_physio = df_mi_all[
        (df_mi_all['architecture'] == 'asymmetric') &
        (df_mi_all['concentration'] >= 300) & (df_mi_all['concentration'] <= 400)
    ]['MI_bits'].values
    mi_sym_physio = df_mi_all[
        (df_mi_all['architecture'] == 'symmetric') &
        (df_mi_all['concentration'] >= 300) & (df_mi_all['concentration'] <= 400)
    ]['MI_bits'].values

    if len(mi_asym_physio) > 0 and len(mi_sym_physio) > 0:
        ma = float(np.mean(mi_asym_physio))
        ms = float(np.mean(mi_sym_physio))
        if ma >= ms:
            n_science_expected += 1
            print(f"  SCIENCE PASS: MI(asymmetric)={ma:.4f} >= MI(symmetric)={ms:.4f} "
                  f"at physiological concentration")
        else:
            n_science_unexpected += 1
            print(f"  SCIENCE UNEXPECTED: MI(asymmetric)={ma:.4f} < MI(symmetric)={ms:.4f} "
                  f"at physiological concentration")

    mi_path = os.path.join(OUTPUT_DIR, 'mutual_information.csv')
    df_mi_all.to_csv(mi_path, index=False)
    output_files.append(mi_path)
    print(f"  Saved {mi_path} ({len(df_mi_all)} rows)")

    # =========================================================================
    # Step 9: Persistence phase diagram
    # =========================================================================
    print("\nStep 9: Generating persistence phase diagram...")

    ratios = PARAMS['sweep_ratios']  # [1, 2, 5, 8, 10, 15, 20, 30, 40, 50]

    df_pd = persistence_phase_diagram(
        ratios=ratios,
        environments=ENVIRONMENTS,
        n_cells_per_point=n_cells_phase_diagram,
        seed=master_seed + 500,
        threshold=threshold_10
    )

    pd_path = os.path.join(OUTPUT_DIR, 'phase_diagram.csv')
    df_pd.to_csv(pd_path, index=False)
    output_files.append(pd_path)
    print(f"  Saved {pd_path} ({len(df_pd)} rows)")

    # =========================================================================
    # Step 10: Sensitivity analysis — threshold at 5th/10th/15th percentile
    # =========================================================================
    print("\nStep 10: Sensitivity analysis for threshold percentiles...")
    sens_rows = []
    for pct_label, thr in [('p05', threshold_5), ('p10', threshold_10), ('p15', threshold_15)]:
        df_s = compute_persister_fractions(pf_results, threshold=thr)
        df_s['threshold_pct'] = pct_label
        sens_rows.append(df_s)
    df_sens = pd.concat(sens_rows, ignore_index=True)
    sens_path = os.path.join(OUTPUT_DIR, 'persistence_sensitivity.csv')
    df_sens.to_csv(sens_path, index=False)
    output_files.append(sens_path)
    print(f"  Saved {sens_path} ({len(df_sens)} rows)")

    # =========================================================================
    # Step 11: Write phase6_summary.json
    # =========================================================================
    wall_time = time.time() - wall_start

    # Check sanity: all output files loadable
    n_missing = sum(1 for f in output_files if not os.path.exists(f))
    if n_missing == 0:
        n_sanity_pass += 1
    else:
        n_sanity_fail += 1
        print(f"  SANITY FAIL: {n_missing} output files missing")

    summary = {
        'status':               'complete',
        'n_sanity_pass':        n_sanity_pass,
        'n_sanity_fail':        n_sanity_fail,
        'n_science_expected':   n_science_expected,
        'n_science_unexpected': n_science_unexpected,
        'n_warn':               n_warn,
        'mock_used':            False,
        'wall_time_sec':        round(wall_time, 2),
        'output_files':         [os.path.basename(f) for f in output_files],
        'n_conditions':         24,
        'n_cells_per_condition': n_cells,
        'phase5_omega':         omega,
        'phase5_dG_spacer':     dG_spacer,
        'mce3r_conc_nM':        mce3r_conc_nM,
        'threshold_10th_pct':   float(threshold_10),
        'threshold_5th_pct':    float(threshold_5),
        'threshold_15th_pct':   float(threshold_15),
        'cv_prediction_failures': cv_fails,
        'pearson_mce3r_target': {
            f'{arch}_{env}': r
            for (arch, env), r in pearson_results.items()
        },
    }

    summary_path = os.path.join(OUTPUT_DIR, 'phase6_summary.json')
    with open(summary_path, 'w') as f:
        json.dump(summary, f, indent=2)
    output_files.append(summary_path)

    print(f"\n{'='*60}")
    print(f"Phase 6 Complete in {wall_time:.1f}s")
    print(f"  Sanity:    {n_sanity_pass} pass, {n_sanity_fail} fail")
    print(f"  Science:   {n_science_expected} expected, "
          f"{n_science_unexpected} unexpected")
    print(f"  Outputs:   {len(output_files)} files in {OUTPUT_DIR}")
    print(f"  Summary:   {summary_path}")

    return summary


# ---------------------------------------------------------------------------
# Self-tests / main execution
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Run Phase 6 Environmental Analysis")
    parser.add_argument('--n-cells', type=int, default=10000,
                        help='Cells per condition (default: 10000)')
    parser.add_argument('--n-cells-mi', type=int, default=5000,
                        help='Cells per MI concentration point (default: 5000)')
    parser.add_argument('--n-cells-pd', type=int, default=2000,
                        help='Cells per phase diagram point (default: 2000)')
    parser.add_argument('--seed', type=int, default=12345)
    args = parser.parse_args()

    summary = run_phase6(
        n_cells=args.n_cells,
        master_seed=args.seed,
        n_cells_mi=args.n_cells_mi,
        n_cells_phase_diagram=args.n_cells_pd
    )

    # Exit with error code if sanity checks failed
    if summary['n_sanity_fail'] > 0:
        print(f"\nSTOPPING: {summary['n_sanity_fail']} sanity check(s) failed.")
        sys.exit(1)
    else:
        print("\nAll sanity checks passed.")
        sys.exit(0)
