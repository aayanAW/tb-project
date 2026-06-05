"""
Aim 1: Local sequence-phenotype tradeoff landscape around the Mce3R operator.

Enumerates all single-base mutations across both 25bp half-sites, predicts
affinity changes via three model tiers, simulates expression, and maps the
Pareto front between growth and persistence.

Key question: Does wild-type sit at a tradeoff optimum?
"""

import sys
import os
import json
import time
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from config.parameters import PARAMS
from phase5_thermodynamic.energy_calibration import build_calibration, score_sequence
from phase2_simulation.operator_model import OperatorModel
from phase_evo.genetic_algorithm import run_population_fast

PROJECT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS_DIR = os.path.join(PROJECT, 'results', 'phase_landscape')
BASES = ['A', 'C', 'G', 'T']


def enumerate_single_mutants(sequence):
    """Generate all single-base mutations of a sequence."""
    mutants = []
    seq = sequence.upper()
    for i in range(len(seq)):
        wt_base = seq[i]
        for alt in BASES:
            if alt == wt_base:
                continue
            mut_seq = seq[:i] + alt + seq[i+1:]
            mutants.append({
                'position': i,
                'wt_base': wt_base,
                'mut_base': alt,
                'mutation': f'{wt_base}{i+1}{alt}',
                'sequence': mut_seq,
            })
    return mutants


def classify_mutation(position, strong_start, strong_end, weak_start, weak_end):
    """Classify mutation by which half-site it affects."""
    if weak_start <= position < weak_end:
        return 'weak_site'
    elif strong_start <= position < strong_end:
        return 'strong_site'
    else:
        return 'spacer'


def predict_kd_three_tiers(mut_seq, wt_seq, calibration, position, site_type):
    """
    Predict Kd change for a mutation using three model tiers.

    Tier 1: Additive energy matrix (simple PWM scoring)
    Tier 2: Structure-informed (position-weighted: contact positions 2x penalty)
    Tier 3: Uncertainty envelope (sample ddG +/- 0.5 kcal/mol)
    """
    kT = calibration['kT']
    energy_matrix = calibration['energy_matrix']
    motif_width = calibration['motif_width']
    offset_mean = calibration['offset_mean']

    # Score wild-type and mutant
    wt_score = score_sequence(wt_seq, energy_matrix)
    mut_score = score_sequence(mut_seq, energy_matrix)
    ddG = mut_score - wt_score  # positive = weaker binding

    # Tier 1: simple additive
    if site_type == 'strong_site':
        wt_kd = PARAMS['Kd_strong']
    elif site_type == 'weak_site':
        wt_kd = PARAMS['Kd_weak']
    else:
        wt_kd = np.sqrt(PARAMS['Kd_strong'] * PARAMS['Kd_weak'])

    tier1_kd = wt_kd * np.exp(ddG / kT)

    # Tier 2: structure-informed (DNA-protein contact positions get 2x weight)
    # Contact positions: positions 2-5, 8-12, 16-19 in a 21bp motif (major groove contacts)
    contact_positions = set(range(2, 6)) | set(range(8, 13)) | set(range(16, 20))
    pos_in_motif = position % motif_width
    weight = 2.0 if pos_in_motif in contact_positions else 0.5
    tier2_ddG = ddG * weight
    tier2_kd = wt_kd * np.exp(tier2_ddG / kT)

    # Tier 3: uncertainty envelope (sample 5 values around ddG)
    tier3_kds = []
    for delta in [-0.5, -0.25, 0.0, 0.25, 0.5]:
        t3_kd = wt_kd * np.exp((ddG + delta) / kT)
        tier3_kds.append(t3_kd)

    return {
        'ddG': ddG,
        'tier1_kd': np.clip(tier1_kd, 0.01, 1e6),
        'tier2_kd': np.clip(tier2_kd, 0.01, 1e6),
        'tier3_kd_low': np.clip(min(tier3_kds), 0.01, 1e6),
        'tier3_kd_mid': np.clip(tier3_kds[2], 0.01, 1e6),
        'tier3_kd_high': np.clip(max(tier3_kds), 0.01, 1e6),
    }


def simulate_variant(Kd_strong, Kd_weak, n_cells=500, seed=42, persister_threshold=112.0):
    """Simulate one operator variant and return growth + persistence metrics."""
    Kd_strong = np.clip(Kd_strong, 0.01, 1e6)
    Kd_weak = np.clip(Kd_weak, 0.01, 1e6)

    model = OperatorModel(Kd_strong=Kd_strong, Kd_weak=Kd_weak)
    arrays = model.get_numba_arrays()
    proteins = run_population_fast(arrays, n_cells, seed, t_max=10000.0, t_burn_in=5000.0)

    pf = proteins.astype(float)
    mean_prot = np.mean(pf)
    cv = np.std(pf, ddof=1) / mean_prot if mean_prot > 0 else 0.0
    persister_frac = float(np.sum(pf < persister_threshold) / len(pf))
    growth_score = mean_prot / 500.0

    return {
        'mean_protein': mean_prot,
        'cv': cv,
        'persister_fraction': persister_frac,
        'growth_score': growth_score,
    }


def is_pareto_efficient(costs):
    """Find Pareto-efficient points (maximize both dimensions)."""
    n = len(costs)
    is_efficient = np.ones(n, dtype=bool)
    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            if costs[j][0] >= costs[i][0] and costs[j][1] >= costs[i][1] and \
               (costs[j][0] > costs[i][0] or costs[j][1] > costs[i][1]):
                is_efficient[i] = False
                break
    return is_efficient


def run_landscape(n_cells=500, seed=42):
    """
    Full Aim 1 pipeline: enumerate mutants, predict affinities, simulate, map landscape.
    """
    os.makedirs(RESULTS_DIR, exist_ok=True)
    start = time.time()

    print("[Landscape] Loading calibration...")
    calibration = build_calibration()
    operator = PARAMS['operator_sequence'].upper()

    # Extract half-site sequences
    strong_seq = operator[88:123]  # last 35bp (contains 25bp strong site)
    weak_seq = operator[0:25]      # first 25bp

    # Use motif-width windows for scoring
    mw = calibration['motif_width']
    strong_window = operator[123-mw:123]  # last motif_width bp
    weak_window = operator[0:mw]          # first motif_width bp

    # Enumerate mutations in both half-sites
    print("[Landscape] Enumerating single-base mutants...")
    strong_mutants = enumerate_single_mutants(strong_window)
    weak_mutants = enumerate_single_mutants(weak_window)

    # Classify and tag
    for m in strong_mutants:
        m['site'] = 'strong_site'
        m['operator_position'] = 123 - mw + m['position']
    for m in weak_mutants:
        m['site'] = 'weak_site'
        m['operator_position'] = m['position']

    all_mutants = strong_mutants + weak_mutants
    print(f"  {len(strong_mutants)} strong-site mutants, {len(weak_mutants)} weak-site mutants")
    print(f"  {len(all_mutants)} total mutants + 1 wild-type")

    # Predict affinities for each mutant
    print("[Landscape] Predicting affinities (3 tiers)...")
    for m in all_mutants:
        wt_window = strong_window if m['site'] == 'strong_site' else weak_window
        kd_pred = predict_kd_three_tiers(
            m['sequence'], wt_window, calibration, m['position'], m['site']
        )
        m.update(kd_pred)

    # Simulate wild-type
    print("[Landscape] Simulating wild-type...")
    wt_result = simulate_variant(PARAMS['Kd_strong'], PARAMS['Kd_weak'], n_cells, seed)
    wt_result['mutation'] = 'WT'
    wt_result['site'] = 'wild_type'
    wt_result['tier1_kd'] = PARAMS['Kd_strong'] if True else PARAMS['Kd_weak']
    wt_result['Kd_strong_used'] = PARAMS['Kd_strong']
    wt_result['Kd_weak_used'] = PARAMS['Kd_weak']

    # Simulate all mutants (Tier 1 only for speed; Tier 2/3 for key metrics)
    print(f"[Landscape] Simulating {len(all_mutants)} mutants ({n_cells} cells each)...")
    rows = []

    for i, m in enumerate(all_mutants):
        if m['site'] == 'strong_site':
            Kd_s = m['tier1_kd']
            Kd_w = PARAMS['Kd_weak']
        else:
            Kd_s = PARAMS['Kd_strong']
            Kd_w = m['tier1_kd']

        sim = simulate_variant(Kd_s, Kd_w, n_cells, seed + i * 10)

        row = {
            'mutation': m['mutation'],
            'site': m['site'],
            'operator_position': m['operator_position'],
            'wt_base': m['wt_base'],
            'mut_base': m['mut_base'],
            'ddG': m['ddG'],
            'tier1_kd': m['tier1_kd'],
            'tier2_kd': m['tier2_kd'],
            'tier3_kd_low': m['tier3_kd_low'],
            'tier3_kd_high': m['tier3_kd_high'],
            'Kd_strong_used': Kd_s,
            'Kd_weak_used': Kd_w,
            'Kd_ratio': Kd_w / Kd_s if Kd_s > 0 else 999,
            'mean_protein': sim['mean_protein'],
            'cv': sim['cv'],
            'persister_fraction': sim['persister_fraction'],
            'growth_score': sim['growth_score'],
        }
        rows.append(row)

        if (i + 1) % 30 == 0:
            print(f"  {i+1}/{len(all_mutants)} done...")

    # Add wild-type row
    wt_row = {
        'mutation': 'WT', 'site': 'wild_type', 'operator_position': -1,
        'wt_base': '-', 'mut_base': '-', 'ddG': 0.0,
        'tier1_kd': PARAMS['Kd_strong'], 'tier2_kd': PARAMS['Kd_strong'],
        'tier3_kd_low': PARAMS['Kd_strong'], 'tier3_kd_high': PARAMS['Kd_strong'],
        'Kd_strong_used': PARAMS['Kd_strong'], 'Kd_weak_used': PARAMS['Kd_weak'],
        'Kd_ratio': PARAMS['Kd_weak'] / PARAMS['Kd_strong'],
        'mean_protein': wt_result['mean_protein'], 'cv': wt_result['cv'],
        'persister_fraction': wt_result['persister_fraction'],
        'growth_score': wt_result['growth_score'],
    }
    rows.append(wt_row)

    df = pd.DataFrame(rows)

    # Pareto analysis
    print("[Landscape] Computing Pareto front...")
    costs = list(zip(df['growth_score'].values, df['persister_fraction'].values))
    pareto_mask = is_pareto_efficient(costs)
    df['is_pareto'] = pareto_mask

    # Check if WT is near Pareto front
    wt_idx = df[df['mutation'] == 'WT'].index[0]
    wt_is_pareto = pareto_mask[wt_idx]

    # Distance to nearest Pareto point if not on front
    if not wt_is_pareto:
        pareto_pts = df[df['is_pareto']][['growth_score', 'persister_fraction']].values
        wt_pt = np.array([wt_row['growth_score'], wt_row['persister_fraction']])
        # Normalize
        g_range = df['growth_score'].max() - df['growth_score'].min()
        p_range = df['persister_fraction'].max() - df['persister_fraction'].min()
        if g_range > 0 and p_range > 0:
            dists = np.sqrt(((pareto_pts[:, 0] - wt_pt[0]) / g_range) ** 2 +
                           ((pareto_pts[:, 1] - wt_pt[1]) / p_range) ** 2)
            wt_pareto_distance = float(np.min(dists))
        else:
            wt_pareto_distance = 0.0
    else:
        wt_pareto_distance = 0.0

    # Mutation class analysis
    strong_muts = df[df['site'] == 'strong_site']
    weak_muts = df[df['site'] == 'weak_site']

    # Mutations that erase asymmetry (bring ratio toward 1)
    df['erases_asymmetry'] = df['Kd_ratio'].apply(lambda r: r < 5.0)
    n_erasing = df['erases_asymmetry'].sum()

    # Save
    csv_path = os.path.join(RESULTS_DIR, 'landscape_data.csv')
    df.to_csv(csv_path, index=False)
    print(f"  Saved: {csv_path}")

    # Ortholog conservation analysis
    print("[Landscape] Checking ortholog conservation...")
    conservation = check_ortholog_conservation()

    # Summary
    elapsed = time.time() - start
    summary = {
        'description': 'Aim 1: Local sequence-phenotype tradeoff landscape',
        'n_mutants': len(all_mutants),
        'n_strong_site': len(strong_mutants),
        'n_weak_site': len(weak_mutants),
        'n_cells_per_variant': n_cells,
        'n_pareto_points': int(pareto_mask.sum()),
        'wt_is_pareto': bool(wt_is_pareto),
        'wt_pareto_distance': wt_pareto_distance,
        'wt_growth_score': wt_row['growth_score'],
        'wt_persister_fraction': wt_row['persister_fraction'],
        'wt_cv': wt_row['cv'],
        'wt_Kd_ratio': wt_row['Kd_ratio'],
        'n_mutations_erasing_asymmetry': int(n_erasing),
        'mean_cv_strong_muts': float(strong_muts['cv'].mean()),
        'mean_cv_weak_muts': float(weak_muts['cv'].mean()),
        'conservation': conservation,
        'wall_time_sec': round(elapsed, 1),
        'n_sanity_pass': 0,
        'n_sanity_fail': 0,
        'n_science_expected': 0,
        'n_science_unexpected': 0,
    }

    # Sanity checks
    if len(df) == len(all_mutants) + 1:
        summary['n_sanity_pass'] += 1
    else:
        summary['n_sanity_fail'] += 1

    if all(df['cv'] > 0) and all(np.isfinite(df['cv'])):
        summary['n_sanity_pass'] += 1
    else:
        summary['n_sanity_fail'] += 1

    # Science: WT near Pareto front
    if wt_is_pareto or wt_pareto_distance < 0.15:
        summary['n_science_expected'] += 1
        summary['wt_near_pareto'] = True
    else:
        summary['n_science_unexpected'] += 1
        summary['wt_near_pareto'] = False

    summary_path = os.path.join(RESULTS_DIR, 'landscape_summary.json')
    with open(summary_path, 'w') as f:
        json.dump(summary, f, indent=2, default=str)

    print(f"\n[Landscape] Complete in {elapsed:.0f}s")
    print(f"  WT on Pareto front: {wt_is_pareto} (distance: {wt_pareto_distance:.3f})")
    print(f"  Pareto points: {pareto_mask.sum()}/{len(df)}")
    print(f"  WT growth={wt_row['growth_score']:.3f}, persist={wt_row['persister_fraction']:.3f}")

    return summary, df


def check_ortholog_conservation():
    """Compare operator across M. tb, M. bovis, M. marinum."""
    conservation_path = os.path.join(PROJECT, 'results', 'phase1', 'conservation_status.csv')
    if not os.path.exists(conservation_path):
        return {'status': 'no_data'}

    df = pd.read_csv(conservation_path)
    n_both = int(df['both_conserved'].sum()) if 'both_conserved' in df.columns else 0
    n_bovis = int(df.get('bovis_conserved', pd.Series([0])).sum())

    return {
        'status': 'available',
        'n_sites_checked': len(df),
        'n_conserved_both_species': n_both,
        'n_conserved_bovis': n_bovis,
        'interpretation': 'High conservation in M. bovis (closely related) supports functional constraint; '
                         'conservation in M. marinum (distant) supports ancient selection.'
    }


if __name__ == '__main__':
    run_landscape()
