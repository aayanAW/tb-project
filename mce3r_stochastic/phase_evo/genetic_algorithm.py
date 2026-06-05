"""
Genetic algorithm for evolving operator architectures under selection pressure.

Each individual is an operator parameterized by (Kd_strong, Kd_weak, block_strong, block_weak).
Two selection regimes:
  - growth_only: fitness = mean protein (favors minimal regulation)
  - persistence: fitness = growth + persister fraction (favors noise-generating asymmetry)

Key prediction: persistence selection evolves asymmetric operators; growth-only does not.
"""

import sys
import os
import time
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from config.parameters import PARAMS
from phase2_simulation.operator_model import OperatorModel
from phase2_simulation.gillespie_engine import simulate_cell


# ============================================================
# FAST GILLESPIE WRAPPER (shared with phase_genome)
# ============================================================

def run_population_fast(model_arrays, n_cells, master_seed, t_max=10000.0, t_burn_in=5000.0):
    """
    Fast Gillespie population simulation with custom timing.

    Calls simulate_cell directly (bypasses run_population's hardcoded PARAMS timing).
    No trace recording. Returns only final protein counts.
    """
    k_txn = model_arrays['k_txn']
    n_bound = model_arrays['n_bound']
    k_on = model_arrays['k_on']
    k_off_s = model_arrays['k_off_strong']
    k_off_w = model_arrays['k_off_weak']
    k_trans = np.float64(PARAMS['k_translation'])
    gamma_m = np.float64(np.log(2) / PARAMS['t_half_mRNA'])
    gamma_p = np.float64(np.log(2) / PARAMS['t_half_protein'])
    nM_mol = np.float64(PARAMS['nM_per_molecule'])

    proteins = np.zeros(n_cells, dtype=np.int64)
    for i in range(n_cells):
        prot, _, _, _, _, _ = simulate_cell(
            k_txn, n_bound, k_on, k_off_s, k_off_w,
            k_trans, gamma_m, gamma_p, nM_mol,
            np.float64(t_max), np.float64(t_burn_in),
            np.float64(100.0), np.int64(200),
            np.int64(master_seed + i), np.int64(0)
        )
        proteins[i] = prot

    return proteins


# ============================================================
# GA PRIMITIVES
# ============================================================

KD_MIN, KD_MAX = 0.1, 1000.0
BLOCK_MIN, BLOCK_MAX = 0.05, 0.99


def initialize_population(pop_size, rng):
    """Create random initial population. Kd log-uniform, block uniform."""
    pop = []
    for _ in range(pop_size):
        ind = {
            'Kd_strong': np.exp(rng.uniform(np.log(0.5), np.log(200.0))),
            'Kd_weak': np.exp(rng.uniform(np.log(0.5), np.log(200.0))),
            'block_strong': rng.uniform(0.1, 0.95),
            'block_weak': rng.uniform(0.1, 0.95),
        }
        # Enforce Kd_strong <= Kd_weak (strong = tighter binding = lower Kd)
        if ind['Kd_strong'] > ind['Kd_weak']:
            ind['Kd_strong'], ind['Kd_weak'] = ind['Kd_weak'], ind['Kd_strong']
            ind['block_strong'], ind['block_weak'] = ind['block_weak'], ind['block_strong']
        pop.append(ind)
    return pop


def mutate(individual, mutation_rate, rng):
    """Mutate an individual. Kd in log-space, block in linear."""
    child = dict(individual)
    for key in ['Kd_strong', 'Kd_weak']:
        if rng.random() < mutation_rate:
            child[key] *= np.exp(rng.normal(0, 0.3))
            child[key] = np.clip(child[key], KD_MIN, KD_MAX)
    for key in ['block_strong', 'block_weak']:
        if rng.random() < mutation_rate:
            child[key] += rng.normal(0, 0.05)
            child[key] = np.clip(child[key], BLOCK_MIN, BLOCK_MAX)
    # Re-enforce ordering
    if child['Kd_strong'] > child['Kd_weak']:
        child['Kd_strong'], child['Kd_weak'] = child['Kd_weak'], child['Kd_strong']
        child['block_strong'], child['block_weak'] = child['block_weak'], child['block_strong']
    return child


def crossover(parent1, parent2, rng):
    """Uniform crossover: each parameter from random parent."""
    child = {}
    for key in ['Kd_strong', 'Kd_weak', 'block_strong', 'block_weak']:
        child[key] = parent1[key] if rng.random() < 0.5 else parent2[key]
    if child['Kd_strong'] > child['Kd_weak']:
        child['Kd_strong'], child['Kd_weak'] = child['Kd_weak'], child['Kd_strong']
        child['block_strong'], child['block_weak'] = child['block_weak'], child['block_strong']
    return child


def tournament_select(population, fitnesses, tournament_size, rng):
    """Select one individual via tournament."""
    indices = rng.choice(len(population), size=tournament_size, replace=False)
    best_idx = indices[np.argmax([fitnesses[i] for i in indices])]
    return dict(population[best_idx])


# ============================================================
# FITNESS EVALUATION
# ============================================================

def evaluate_fitness(individual, regime, n_cells, seed, persister_threshold=112.0):
    """
    Evaluate fitness of one operator architecture.

    regime='growth_only': fitness = mean_protein / 500
    regime='persistence': fitness = 0.7 * growth + 0.3 * persister_fraction
    """
    model = OperatorModel(
        Kd_strong=individual['Kd_strong'],
        Kd_weak=individual['Kd_weak'],
        block_strong=individual['block_strong'],
        block_weak=individual['block_weak'],
    )
    model_arrays = model.get_numba_arrays()
    proteins = run_population_fast(model_arrays, n_cells, seed)

    proteins_f = proteins.astype(float)
    mean_prot = np.mean(proteins_f)
    std_prot = np.std(proteins_f, ddof=1) if len(proteins_f) > 1 else 0.0
    cv = std_prot / mean_prot if mean_prot > 0 else 0.0
    persister_frac = float(np.sum(proteins_f < persister_threshold)) / len(proteins_f)

    growth_score = mean_prot / 500.0

    if regime == 'growth_only':
        fitness = growth_score
    elif regime == 'persistence':
        fitness = 0.7 * growth_score + 0.3 * persister_frac
    else:
        fitness = growth_score

    return {
        'fitness': fitness,
        'mean_protein': mean_prot,
        'cv': cv,
        'persister_fraction': persister_frac,
        'Kd_ratio': individual['Kd_weak'] / individual['Kd_strong'] if individual['Kd_strong'] > 0 else 1.0,
    }


# ============================================================
# MAIN GA LOOP
# ============================================================

def run_ga(regime, pop_size=50, n_generations=100, n_cells=200,
           mutation_rate=0.8, tournament_size=3, elite_count=2,
           seed=42, persister_threshold=112.0):
    """
    Run genetic algorithm under specified selection regime.

    Returns
    -------
    trajectory : pd.DataFrame
        Per-generation statistics.
    final_pop : pd.DataFrame
        Final generation individuals with parameters and fitness.
    """
    rng = np.random.RandomState(seed)
    population = initialize_population(pop_size, rng)

    trajectory_rows = []

    print(f"  GA [{regime}]: {pop_size} pop x {n_generations} gen x {n_cells} cells")

    for gen in range(n_generations):
        gen_seed = seed + gen * pop_size * 10

        # Evaluate all individuals
        evals = []
        for i, ind in enumerate(population):
            ev = evaluate_fitness(ind, regime, n_cells, gen_seed + i * 10, persister_threshold)
            evals.append(ev)

        fitnesses = [e['fitness'] for e in evals]
        kd_ratios = [e['Kd_ratio'] for e in evals]
        cvs = [e['cv'] for e in evals]

        best_idx = np.argmax(fitnesses)

        row = {
            'generation': gen,
            'mean_fitness': np.mean(fitnesses),
            'best_fitness': fitnesses[best_idx],
            'mean_Kd_ratio': np.mean(kd_ratios),
            'std_Kd_ratio': np.std(kd_ratios),
            'best_Kd_ratio': kd_ratios[best_idx],
            'mean_cv': np.mean(cvs),
            'best_cv': cvs[best_idx],
            'mean_persister': np.mean([e['persister_fraction'] for e in evals]),
            'diversity': np.std([np.log(r + 1e-6) for r in kd_ratios]),
        }
        trajectory_rows.append(row)

        if gen % 20 == 0 or gen == n_generations - 1:
            print(f"    Gen {gen:3d}: fitness={row['mean_fitness']:.3f} "
                  f"best={row['best_fitness']:.3f} "
                  f"Kd_ratio={row['mean_Kd_ratio']:.1f} "
                  f"CV={row['mean_cv']:.4f}")

        # Create next generation
        # Sort by fitness for elitism
        sorted_indices = np.argsort(fitnesses)[::-1]
        new_population = []

        # Elitism
        for i in range(elite_count):
            new_population.append(dict(population[sorted_indices[i]]))

        # Fill rest via tournament + crossover + mutation
        while len(new_population) < pop_size:
            p1 = tournament_select(population, fitnesses, tournament_size, rng)
            p2 = tournament_select(population, fitnesses, tournament_size, rng)
            child = crossover(p1, p2, rng)
            child = mutate(child, mutation_rate, rng)
            new_population.append(child)

        population = new_population

    # Final evaluation for output
    final_rows = []
    for i, ind in enumerate(population):
        ev = evaluate_fitness(ind, regime, n_cells, seed + 999999 + i * 10, persister_threshold)
        final_rows.append({
            'individual': i,
            'Kd_strong': ind['Kd_strong'],
            'Kd_weak': ind['Kd_weak'],
            'block_strong': ind['block_strong'],
            'block_weak': ind['block_weak'],
            'Kd_ratio': ev['Kd_ratio'],
            'fitness': ev['fitness'],
            'mean_protein': ev['mean_protein'],
            'cv': ev['cv'],
            'persister_fraction': ev['persister_fraction'],
        })

    trajectory = pd.DataFrame(trajectory_rows)
    final_pop = pd.DataFrame(final_rows)

    return trajectory, final_pop
