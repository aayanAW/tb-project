"""
config/parameters.py — All tunable parameters for the Mce3R project.

Every parameter includes its source reference.
This is the single source of truth for all configurable values.
"""

import numpy as np

PARAMS = {
    # === Genome Information ===
    'H37Rv_accession': 'NC_000962.3',
    'H37Rv_genome_size': 4_411_532,
    'H37Rv_gc_content': 0.656,
    'M_bovis_accession': 'NC_002945.4',
    'M_marinum_accession': 'NC_010612.1',

    # === Ortholog Identifiers (for cross-species MEME input) ===
    # yrbE3A orthologs in other species — used by extract_upstream.py
    'yrbE3A_bovis_locus': 'Mb1997',
    'yrbE3A_marinum_locus': 'MMAR_2522',

    # === Mce3R Binding Parameters (Panagoda et al. 2024) ===
    'Kd_strong': 2.4,              # nM, Table 1, Probe A
    'Kd_weak': 49.0,               # nM, Table 1, Probe C
    'k_on': 0.0167,                # nM^-1 min^-1, Stormo & Zhao 2010
                                   # (10^6 M^-1 s^-1 × 60 s/min × 10^-9 M/nM)

    # === Cell Volume & Concentration Conversion ===
    'cell_volume_fL': 1.0,         # femtoliters, standard Mtb cell volume
    'nM_per_molecule': 1.66,       # 1 molecule in 1 fL ≈ 1.66 nM
                                   # (from 1/(Avogadro × 1e-15 L) × 1e9)

    # === Transcription Model ===
    'k_max': 0.15,                 # mRNA/min, estimated from Mtb transcriptomics
    'block_strong': 0.85,          # fraction blocked when strong site occupied
    'block_weak': 0.50,            # fraction blocked when weak site occupied

    # === Translation & Degradation ===
    'k_translation': 0.5,          # protein per mRNA per min, Taniguchi 2010 adjusted
    't_half_mRNA': 9.5,            # min, Rustad et al. 2013 NAR
    't_half_protein': 1500.0,      # min (~25 hr), dominated by Mtb dilution rate

    # === Simulation Settings ===
    'n_cells_main': 50_000,        # cells per condition (A, B, C, D)
    'n_cells_sweep': 5_000,        # cells per sweep point (Condition E)
    't_max': 30_000,               # min (500 hr ≈ 25 doubling times, 20× protein half-life)
    't_burn_in': 15_000,           # min, discard first half (10× protein half-life for steady state)

    # === Trace Recording (for fig6) ===
    'n_trace_cells': 20,           # number of cells to save full trajectories for
    'trace_interval': 10.0,        # minutes between recorded time points
    'trace_max_points': 1_600,     # pre-allocated: (t_max-t_burn_in)/trace_interval = 1500 + margin

    # === Reproducibility Seeds ===
    'master_seed': 42,
    'gmm_random_state': 42,
    'bootstrap_seed': 123,
    # Condition-level seeds: master_seed + condition_index (A=0, B=1, C=2, D=3)
    # Sweep-level seeds: master_seed + 100 + ratio_index

    # === Asymmetry Sweep (Condition E) ===
    'sweep_ratios': [1, 2, 5, 8, 10, 15, 20, 30, 40, 50],
    'mce3r_actual_ratio': 20.4,    # Kd_weak/Kd_strong = 49.0/2.4

    # === Symmetric Control (Condition B) ===
    # Geometric mean of strong/weak values so total regulatory capacity is equivalent
    # to the asymmetric case. This ensures a fair noise comparison.
    'Kd_symmetric': float(np.sqrt(2.4 * 49.0)),    # ≈ 10.84 nM, geometric mean
    'block_symmetric': float(np.sqrt(0.85 * 0.50)), # ≈ 0.652, geometric mean

    # === Single-Site Control (Condition C) ===
    'Kd_weak_disabled': 1e12,      # nM, effectively infinite — weak site never occupied

    # === No-Regulation Control (Condition D) ===
    'k_on_disabled': 0.0,          # nM^-1 min^-1, no binding at all

    # === Sensitivity Analysis ===
    'sensitivity_range': 0.50,     # ±50% variation
    'sensitivity_steps': 10,       # number of steps per parameter
    'sensitivity_cells': 1_000,    # cells per step

    # === Statistical Parameters ===
    'n_bootstrap': 1_000,
    'ci_level': 0.95,
    'n_gmm_components': [1, 2, 3],
    'gmm_reg_covar': 1e-3,
    'gmm_n_init': 5,
    'shuffle_n': 1_000,
    'shuffle_alpha': 0.01,

    # === MEME/FIMO Parameters ===
    'meme_mod': 'zoops',
    'meme_minw': 6,
    'meme_maxw': 110,
    'meme_nmotifs': 5,
    'fimo_thresh': 1e-4,
    'upstream_length': 200,
    'known_spacer_bp': 53,              # Paper: two 25bp sites separated by 53bp
    'weak_site_start_in_operator': 10,  # Approximate start of weak site in 123bp operator
    'strong_site_start_in_operator': 88, # Approximate start of strong site in 123bp operator
    'mast_ev_threshold': 100,           # MAST E-value threshold
    'mast_mt_threshold': 1e-4,          # MAST motif p-value threshold

    # === Operator Sequence (PDB 9B7Y, Panagoda 2024) ===
    'operator_sequence': (
        'GCCCCGCGCTATAGGATACTAGCAAGATACATCATAGCCAATATATGCCAGT'
        'TTGCATTGCTATTTACCGATCAGTTGTCCAAGCAATCGCGTATTGGCTATG'
        'GACATCAGCGGTTCTGCCGC'
    ),
    'operator_length': 123,
    'operator_region_h37rv_start': 2_207_477,  # approximate intergenic region start
    'operator_region_h37rv_end': 2_207_699,    # approximate intergenic region end
    # NOTE: span is ~222 bp (larger than 123-bp operator). This is the broader
    # intergenic window containing the operator, not exact 123-bp boundaries.

    # === Gene Identifiers ===
    'mce3R_gene': 'Rv1963c',
    'yrbE3A_gene': 'Rv1964',
    'mce3_operon_genes': [f'Rv{i}' for i in range(1964, 1978)],
    'regulon_genes_1': ['Rv1933c', 'Rv1934c', 'Rv1935c'],
    'regulon_genes_2': [f'Rv{i}' for i in range(1936, 1942)],
}

# === Derived parameters (calculated, not configurable) ===
PARAMS['gamma_mRNA'] = np.log(2) / PARAMS['t_half_mRNA']           # 0.0730 min^-1
PARAMS['gamma_protein'] = np.log(2) / PARAMS['t_half_protein']      # 0.000462 min^-1
PARAMS['burst_size'] = PARAMS['k_translation'] / PARAMS['gamma_mRNA']  # ~6.85
PARAMS['expected_fano_unregulated'] = 1 + PARAMS['burst_size']       # ~7.85

# Derived binding rate constants
PARAMS['k_off_strong'] = PARAMS['Kd_strong'] * PARAMS['k_on']  # 2.4 × 0.0167 = 0.0401 min^-1
PARAMS['k_off_weak'] = PARAMS['Kd_weak'] * PARAMS['k_on']      # 49.0 × 0.0167 = 0.818 min^-1

# Transcription rates per operator state (length-4 float64 array)
PARAMS['k_txn'] = np.array([
    PARAMS['k_max'],                                                                     # State 0: both empty
    PARAMS['k_max'] * (1 - PARAMS['block_strong']),                                      # State 1: strong occupied
    PARAMS['k_max'] * (1 - PARAMS['block_weak']),                                        # State 2: weak occupied
    PARAMS['k_max'] * (1 - PARAMS['block_strong']) * (1 - PARAMS['block_weak']),         # State 3: both occupied
])


# ============================================================
# PHASE 5-7 PARAMETERS (Thermodynamic + Environmental Extensions)
# ============================================================
TEMPERATURE_K = 310.0
kT_KCAL = 0.001987 * TEMPERATURE_K  # ≈ 0.616 kcal/mol
GC_CONTENT = PARAMS['H37Rv_gc_content']
BACKGROUND_FREQ = {
    'A': (1 - GC_CONTENT) / 2,
    'C': GC_CONTENT / 2,
    'G': GC_CONTENT / 2,
    'T': (1 - GC_CONTENT) / 2,
}
MCMC_N_WALKERS = 32
MCMC_N_STEPS = 5000
MCMC_N_BURN = 1000
MCMC_SEED = 42
STRONG_SITE_START = 98
STRONG_SITE_END = 123
WEAK_SITE_START = 0
WEAK_SITE_END = 25
ENV_COLORS = {
    'baseline': '#6B7280',
    'cholesterol': '#F59E0B',
    'acidic_pH': '#EF4444',
    'host_like': '#7C3AED',
}


if __name__ == "__main__":
    """Self-test block for parameters.py."""
    import sys

    n_sanity_pass = 0
    n_sanity_fail = 0

    # Sanity check 1: All Kd values positive
    try:
        assert PARAMS['Kd_strong'] > 0 and PARAMS['Kd_weak'] > 0
        print(f"SANITY PASS: Kd values positive — Kd_strong={PARAMS['Kd_strong']}, Kd_weak={PARAMS['Kd_weak']}")
        n_sanity_pass += 1
    except AssertionError:
        print(f"SANITY FAIL: Kd values must be positive")
        n_sanity_fail += 1

    # Sanity check 2: k_off values match calculation
    try:
        expected_k_off_s = PARAMS['Kd_strong'] * PARAMS['k_on']
        expected_k_off_w = PARAMS['Kd_weak'] * PARAMS['k_on']
        assert abs(PARAMS['k_off_strong'] - expected_k_off_s) < 1e-10
        assert abs(PARAMS['k_off_weak'] - expected_k_off_w) < 1e-10
        print(f"SANITY PASS: k_off values match — k_off_strong={PARAMS['k_off_strong']:.4f}, k_off_weak={PARAMS['k_off_weak']:.4f}")
        n_sanity_pass += 1
    except AssertionError:
        print(f"SANITY FAIL: k_off values do not match Kd × k_on calculation")
        n_sanity_fail += 1

    # Sanity check 3: Transcription rates array has 4 elements, all positive
    try:
        assert len(PARAMS['k_txn']) == 4
        assert all(PARAMS['k_txn'] > 0)
        print(f"SANITY PASS: k_txn array length=4, all positive — {PARAMS['k_txn']}")
        n_sanity_pass += 1
    except (AssertionError, Exception) as e:
        print(f"SANITY FAIL: k_txn array issue — {e}")
        n_sanity_fail += 1

    # Sanity check 4: Transcription rates monotonic (state 0 > state 2 > state 1 > state 3)
    try:
        k = PARAMS['k_txn']
        assert k[0] > k[2] > k[1] > k[3]
        print(f"SANITY PASS: k_txn monotonic — {k[0]:.4f} > {k[2]:.4f} > {k[1]:.4f} > {k[3]:.4f}")
        n_sanity_pass += 1
    except AssertionError:
        print(f"SANITY FAIL: k_txn not monotonic as expected")
        n_sanity_fail += 1

    # Sanity check 5: Derived rates are reasonable
    try:
        assert 0.05 < PARAMS['gamma_mRNA'] < 0.15
        assert 0.0001 < PARAMS['gamma_protein'] < 0.001
        assert 5 < PARAMS['burst_size'] < 10
        print(f"SANITY PASS: Derived rates — gamma_mRNA={PARAMS['gamma_mRNA']:.4f}, gamma_protein={PARAMS['gamma_protein']:.6f}, burst_size={PARAMS['burst_size']:.2f}")
        n_sanity_pass += 1
    except AssertionError:
        print(f"SANITY FAIL: Derived rates out of expected range")
        n_sanity_fail += 1

    # Sanity check 6: Simulation time is consistent
    try:
        assert PARAMS['t_max'] == 30_000
        assert PARAMS['t_burn_in'] == 15_000
        assert PARAMS['t_burn_in'] < PARAMS['t_max']
        max_trace = int((PARAMS['t_max'] - PARAMS['t_burn_in']) / PARAMS['trace_interval'])
        assert max_trace <= PARAMS['trace_max_points']
        print(f"SANITY PASS: Simulation time — t_max={PARAMS['t_max']}, t_burn_in={PARAMS['t_burn_in']}, max_trace_points={max_trace}")
        n_sanity_pass += 1
    except AssertionError:
        print(f"SANITY FAIL: Simulation time parameters inconsistent")
        n_sanity_fail += 1

    # Sanity check 7: Operator sequence length
    try:
        assert len(PARAMS['operator_sequence']) == PARAMS['operator_length']
        print(f"SANITY PASS: Operator sequence length = {len(PARAMS['operator_sequence'])}")
        n_sanity_pass += 1
    except AssertionError:
        print(f"SANITY FAIL: Operator sequence length mismatch — got {len(PARAMS['operator_sequence'])}, expected {PARAMS['operator_length']}")
        n_sanity_fail += 1

    # Sanity check 8: nM_per_molecule calculation
    try:
        # 1 molecule in 1 fL: concentration = 1 / (6.022e23 * 1e-15) * 1e9 nM ≈ 1.66 nM
        expected = 1.0 / (6.022e23 * 1e-15) * 1e9
        assert abs(PARAMS['nM_per_molecule'] - expected) < 0.01
        print(f"SANITY PASS: nM_per_molecule = {PARAMS['nM_per_molecule']:.2f} (calculated {expected:.2f})")
        n_sanity_pass += 1
    except AssertionError:
        print(f"SANITY FAIL: nM_per_molecule incorrect")
        n_sanity_fail += 1

    # Summary
    print(f"\n{'='*50}")
    print(f"Sanity checks: {n_sanity_pass} PASS, {n_sanity_fail} FAIL")
    if n_sanity_fail > 0:
        print("STOPPING: Software sanity check failed.")
        sys.exit(1)
    else:
        print("All checks passed.")
        sys.exit(0)
