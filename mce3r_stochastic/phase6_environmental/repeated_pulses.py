"""
Aim 3 extension: Repeated antibiotic pulse survival model.

Simulates 3-5 rounds of antibiotic exposure with regrowth between rounds.
Tests whether asymmetry helps survive repeated treatment, not just one snapshot.

Also expands the environment set to include oxidative stress and pulsed schedules.
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
RESULTS_DIR = os.path.join(PROJECT, 'results', 'phase6')


ENVIRONMENTS = {
    'baseline': {'k_on_mult': 1.0, 'gamma_mult': 1.0},
    'cholesterol': {'k_on_mult': 0.3, 'gamma_mult': 1.0},
    'acidic_pH': {'k_on_mult': 0.7, 'gamma_mult': 1.2},
    'oxidative': {'k_on_mult': 0.5, 'gamma_mult': 1.4},
    'combined_macrophage': {'k_on_mult': 0.2, 'gamma_mult': 1.3},
}

PULSE_DURATIONS_HR = [2, 6, 12, 24]
N_ROUNDS = 4
KD_RATIOS = [1.0, 2.0, 5.0, 10.0, 20.4, 30.0, 50.0]


def simulate_population_with_env(Kd_strong, Kd_weak, k_on_mult, gamma_mult,
                                  n_cells=500, seed=42):
    """Simulate a population under environmental stress."""
    model = OperatorModel(
        Kd_strong=Kd_strong, Kd_weak=Kd_weak,
        k_on=PARAMS['k_on'] * k_on_mult,
    )
    arrays = model.get_numba_arrays()

    k_trans = np.float64(PARAMS['k_translation'])
    gamma_m = np.float64(np.log(2) / PARAMS['t_half_mRNA'])
    gamma_p = np.float64(np.log(2) / PARAMS['t_half_protein'] / gamma_mult)
    nM_mol = np.float64(PARAMS['nM_per_molecule'])

    proteins = np.zeros(n_cells, dtype=np.int64)
    for i in range(n_cells):
        prot, _, _, _, _, _ = simulate_cell(
            arrays['k_txn'], arrays['n_bound'], arrays['k_on'],
            arrays['k_off_strong'], arrays['k_off_weak'],
            k_trans, gamma_m, gamma_p, nM_mol,
            np.float64(10000), np.float64(5000),
            np.float64(100), np.int64(200),
            np.int64(seed + i), np.int64(0)
        )
        proteins[i] = prot

    return proteins.astype(float)


def apply_antibiotic_pulse(proteins, threshold, kill_prob_above=0.8):
    """
    One round of antibiotic exposure.
    Cells below threshold survive with high probability (persisters).
    Cells above threshold die with kill_prob_above probability.
    """
    rng = np.random.RandomState(hash(tuple(proteins[:10].astype(int))) % (2**31))
    survivors = []
    for p in proteins:
        if p < threshold:
            # Persister — survives with 95% probability
            if rng.random() > 0.05:
                survivors.append(p)
        else:
            # Susceptible — dies with kill_prob_above probability
            if rng.random() > kill_prob_above:
                survivors.append(p)
    return np.array(survivors) if survivors else np.array([0.0])


def regrow_population(survivors, target_size, mean_protein, cv):
    """
    Simulate regrowth from survivors back to target population size.
    New cells drawn from distribution centered on survivor characteristics.
    """
    if len(survivors) == 0 or (len(survivors) == 1 and survivors[0] == 0):
        return np.array([])

    n_new = target_size - len(survivors)
    if n_new <= 0:
        return survivors[:target_size]

    # Survivors seed the new population — resample with noise
    rng = np.random.RandomState(42)
    new_cells = rng.choice(survivors, size=n_new, replace=True)
    noise = rng.normal(0, cv * np.mean(survivors), size=n_new)
    new_cells = np.clip(new_cells + noise, 1, None)

    return np.concatenate([survivors, new_cells])


def run_repeated_pulses(Kd_strong, Kd_weak, env_name, env_params,
                        n_rounds=4, n_cells=500, seed=42):
    """
    Run multi-round antibiotic treatment simulation.

    Returns survival fraction after each round.
    """
    # Initial population
    proteins = simulate_population_with_env(
        Kd_strong, Kd_weak,
        env_params['k_on_mult'], env_params['gamma_mult'],
        n_cells=n_cells, seed=seed
    )

    mean_prot = np.mean(proteins)
    cv = np.std(proteins) / mean_prot if mean_prot > 0 else 0.1

    # Threshold from Condition B baseline
    threshold = mean_prot * 0.4  # cells below 40% of mean are "persisters"

    round_data = []
    current_pop = proteins.copy()
    initial_size = len(current_pop)

    for r in range(n_rounds):
        pre_size = len(current_pop)
        if pre_size == 0:
            round_data.append({
                'round': r, 'pre_size': 0, 'survivors': 0,
                'survival_frac': 0.0, 'cumulative_survival': 0.0,
            })
            continue

        survivors = apply_antibiotic_pulse(current_pop, threshold)
        n_surv = len(survivors)
        survival_frac = n_surv / pre_size if pre_size > 0 else 0.0
        cumulative = n_surv / initial_size if initial_size > 0 else 0.0

        round_data.append({
            'round': r,
            'pre_size': pre_size,
            'survivors': n_surv,
            'survival_frac': survival_frac,
            'cumulative_survival': cumulative,
        })

        # Regrow for next round
        current_pop = regrow_population(survivors, n_cells, mean_prot, cv)

    final_survival = round_data[-1]['cumulative_survival'] if round_data else 0.0
    return round_data, final_survival


def run_phase_diagram_extended(n_cells=500, seed=42):
    """
    Build extended phase diagram: Kd ratio x environment x pulse duration.
    """
    os.makedirs(RESULTS_DIR, exist_ok=True)
    start = time.time()

    Kd_geo = PARAMS['Kd_symmetric']  # 10.84 nM
    rows = []

    total = len(KD_RATIOS) * len(ENVIRONMENTS)
    done = 0

    print(f"[Pulses] Running {len(KD_RATIOS)} ratios x {len(ENVIRONMENTS)} environments...")

    for ratio in KD_RATIOS:
        Kd_s = Kd_geo / np.sqrt(ratio)
        Kd_w = Kd_geo * np.sqrt(ratio)

        for env_name, env_params in ENVIRONMENTS.items():
            done += 1
            round_data, final_survival = run_repeated_pulses(
                Kd_s, Kd_w, env_name, env_params,
                n_rounds=N_ROUNDS, n_cells=n_cells, seed=seed + done * 100
            )

            # Also compute single-snapshot metrics
            proteins = simulate_population_with_env(
                Kd_s, Kd_w,
                env_params['k_on_mult'], env_params['gamma_mult'],
                n_cells=n_cells, seed=seed + done * 100 + 50
            )
            mean_prot = float(np.mean(proteins))
            cv = float(np.std(proteins) / mean_prot) if mean_prot > 0 else 0.0

            row = {
                'Kd_ratio': ratio,
                'Kd_strong': Kd_s,
                'Kd_weak': Kd_w,
                'environment': env_name,
                'k_on_mult': env_params['k_on_mult'],
                'gamma_mult': env_params['gamma_mult'],
                'mean_protein': mean_prot,
                'cv': cv,
                'growth_score': mean_prot / 500.0,
                'multi_round_survival': final_survival,
                'round1_survival': round_data[0]['survival_frac'] if round_data else 0,
                'round4_survival': round_data[-1]['survival_frac'] if len(round_data) >= 4 else 0,
                'n_rounds': N_ROUNDS,
            }
            rows.append(row)

            if done % 10 == 0:
                print(f"  {done}/{total}: ratio={ratio:.1f}, env={env_name}, "
                      f"multi_surv={final_survival:.4f}")

    df = pd.DataFrame(rows)
    csv_path = os.path.join(RESULTS_DIR, 'extended_phase_diagram.csv')
    df.to_csv(csv_path, index=False)

    elapsed = time.time() - start

    # Summary
    wt_ratio = 20.4
    wt_rows = df[np.abs(df['Kd_ratio'] - wt_ratio) < 1.0]
    sym_rows = df[np.abs(df['Kd_ratio'] - 1.0) < 0.5]

    summary = {
        'description': 'Aim 3: Extended phase diagram with repeated pulses',
        'n_ratios': len(KD_RATIOS),
        'n_environments': len(ENVIRONMENTS),
        'n_rounds': N_ROUNDS,
        'n_cells': n_cells,
        'wall_time_sec': round(elapsed, 1),
        'n_sanity_pass': 0,
        'n_sanity_fail': 0,
        'n_science_expected': 0,
        'n_science_unexpected': 0,
    }

    # Sanity
    if len(df) == len(KD_RATIOS) * len(ENVIRONMENTS):
        summary['n_sanity_pass'] += 1
    else:
        summary['n_sanity_fail'] += 1

    if all(df['cv'] > 0) and all(np.isfinite(df['cv'])):
        summary['n_sanity_pass'] += 1
    else:
        summary['n_sanity_fail'] += 1

    # Science: asymmetry should help survival under stress more than baseline
    for env in ['cholesterol', 'acidic_pH', 'combined_macrophage']:
        env_data = df[df['environment'] == env]
        if len(env_data) > 0:
            asym = env_data[np.abs(env_data['Kd_ratio'] - 20.4) < 1.0]
            sym = env_data[np.abs(env_data['Kd_ratio'] - 1.0) < 0.5]
            if len(asym) > 0 and len(sym) > 0:
                if asym.iloc[0]['multi_round_survival'] > sym.iloc[0]['multi_round_survival']:
                    summary['n_science_expected'] += 1
                else:
                    summary['n_science_unexpected'] += 1

    # Where is asymmetry beneficial vs harmful?
    benefit_envs = []
    for env in ENVIRONMENTS:
        env_data = df[df['environment'] == env]
        asym = env_data[np.abs(env_data['Kd_ratio'] - 20.4) < 1.0]
        sym = env_data[np.abs(env_data['Kd_ratio'] - 1.0) < 0.5]
        if len(asym) > 0 and len(sym) > 0:
            a_surv = asym.iloc[0]['multi_round_survival']
            s_surv = sym.iloc[0]['multi_round_survival']
            benefit_envs.append({
                'environment': env,
                'asymmetric_survival': float(a_surv),
                'symmetric_survival': float(s_surv),
                'asymmetry_beneficial': a_surv > s_surv,
                'fold_advantage': float(a_surv / s_surv) if s_surv > 0 else float('inf'),
            })

    summary['environment_benefits'] = benefit_envs

    summary_path = os.path.join(RESULTS_DIR, 'extended_phase_summary.json')
    with open(summary_path, 'w') as f:
        json.dump(summary, f, indent=2, default=str)

    print(f"\n[Pulses] Complete in {elapsed:.0f}s")
    print(f"  Sanity: {summary['n_sanity_pass']} pass, {summary['n_sanity_fail']} fail")
    print(f"  Science: {summary['n_science_expected']} expected, {summary['n_science_unexpected']} unexpected")

    if benefit_envs:
        print("\n  Asymmetry benefit by environment:")
        for b in benefit_envs:
            marker = '+' if b['asymmetry_beneficial'] else '-'
            print(f"    [{marker}] {b['environment']}: asym={b['asymmetric_survival']:.4f} "
                  f"sym={b['symmetric_survival']:.4f} ({b['fold_advantage']:.1f}x)")

    return summary, df


if __name__ == '__main__':
    run_phase_diagram_extended()
