"""
phase6_environmental/two_species_gillespie.py — Numba-accelerated two-species Gillespie engine.

18 reactions for the Mce3R + Target two-species circuit:

Operator transitions (8):
  0: U→S  (state 0→1): propensity = k_on * free_mce3r_nM  [if op_state==0]
  1: S→U  (state 1→0): propensity = k_off_strong           [if op_state==1]
  2: U→W  (state 0→2): propensity = k_on * free_mce3r_nM  [if op_state==0]
  3: W→U  (state 2→0): propensity = k_off_weak             [if op_state==2]
  4: W→D  (state 2→3): propensity = k_on * free_mce3r_nM  [if op_state==2]
  5: D→W  (state 3→2): propensity = k_off_strong           [if op_state==3]
  6: S→D  (state 1→3): propensity = k_on * free_mce3r_nM  [if op_state==1]
  7: D→S  (state 3→1): propensity = k_off_weak             [if op_state==3]

Transcription (2):
  8:  Mce3R mRNA production: k_txn[op_state] * alpha_mce3r
  9:  Target mRNA production: k_txn[op_state] * alpha_target

Translation (2):
  10: Mce3R translation: k_translation * mce3r_mRNA
  11: Target translation: k_translation * target_mRNA

mRNA decay (2):
  12: Mce3R mRNA decay: gamma_mRNA * mce3r_mRNA
  13: Target mRNA decay: gamma_mRNA * target_mRNA

Protein decay (2):
  14: Mce3R protein decay: gamma_protein * mce3r_protein
  15: Target protein decay: gamma_protein * target_protein

Burst placeholders (2):
  16, 17: propensity = 0.0

where free_mce3r_nM = max(0, mce3r_protein - n_bound[op_state]) * nM_per_molecule

Critical: Mce3R protein count dynamically determines binding propensity (TRUE autoregulation).
"""

import sys
import numpy as np
import numba

sys.path.insert(0, '/Users/aayanalwani/tb project/mce3r_stochastic')
from config.parameters import PARAMS


@numba.njit
def simulate_two_species_cell(
    k_txn,
    n_bound,
    k_on,
    k_off_strong,
    k_off_weak,
    k_translation,
    gamma_mRNA,
    gamma_protein,
    nM_per_molecule,
    alpha_mce3r,
    alpha_target,
    t_max,
    t_burn_in,
    seed
):
    """
    Run one Gillespie SSA trajectory for a two-species cell.

    All inputs are scalars or numpy arrays (Numba-safe — no dicts/classes).

    Parameters
    ----------
    k_txn          : float64[4]  — transcription rate per operator state
    n_bound        : int64[4]    — Mce3R molecules bound per state
    k_on           : float64     — binding rate constant (nM^-1 min^-1)
    k_off_strong   : float64     — unbinding rate, strong site (min^-1)
    k_off_weak     : float64     — unbinding rate, weak site (min^-1)
    k_translation  : float64     — translation rate (protein/mRNA/min)
    gamma_mRNA     : float64     — mRNA degradation rate (min^-1)
    gamma_protein  : float64     — protein degradation rate (min^-1)
    nM_per_molecule: float64     — nM per single molecule in 1 fL cell
    alpha_mce3r    : float64     — transcription scale for Mce3R direction
    alpha_target   : float64     — transcription scale for Target direction
    t_max          : float64     — total simulation time (min)
    t_burn_in      : float64     — burn-in time; record state only after this (min)
    seed           : int64       — random seed

    Returns
    -------
    (final_mce3r_protein, final_target_protein,
     final_mce3r_mRNA, final_target_mRNA, final_op_state) : tuple of int
    """
    np.random.seed(seed)

    # Initial conditions: pre-loaded to speed approach to steady state
    op_state       = 0
    mce3r_mRNA     = 1
    mce3r_protein  = 50
    target_mRNA    = 1
    target_protein = 50

    t = 0.0

    # Propensity array (18 reactions)
    a = np.zeros(18, dtype=np.float64)

    while t < t_max:
        # Compute free Mce3R protein available for binding
        nb = n_bound[op_state]
        free_mce3r = mce3r_protein - nb
        if free_mce3r < 0:
            free_mce3r = 0
        free_mce3r_nM = free_mce3r * nM_per_molecule

        # --- Operator transitions ---
        # Reaction 0: U→S (0→1): bind strong site from unbound state
        a[0] = k_on * free_mce3r_nM if op_state == 0 else 0.0
        # Reaction 1: S→U (1→0): unbind strong site
        a[1] = k_off_strong if op_state == 1 else 0.0
        # Reaction 2: U→W (0→2): bind weak site from unbound state
        a[2] = k_on * free_mce3r_nM if op_state == 0 else 0.0
        # Reaction 3: W→U (2→0): unbind weak site
        a[3] = k_off_weak if op_state == 2 else 0.0
        # Reaction 4: W→D (2→3): bind strong when weak occupied
        a[4] = k_on * free_mce3r_nM if op_state == 2 else 0.0
        # Reaction 5: D→W (3→2): unbind strong when both occupied
        a[5] = k_off_strong if op_state == 3 else 0.0
        # Reaction 6: S→D (1→3): bind weak when strong occupied
        a[6] = k_on * free_mce3r_nM if op_state == 1 else 0.0
        # Reaction 7: D→S (3→1): unbind weak when both occupied
        a[7] = k_off_weak if op_state == 3 else 0.0

        # --- Transcription (both genes, same operator state) ---
        # Reaction 8: Mce3R mRNA production
        a[8] = k_txn[op_state] * alpha_mce3r
        # Reaction 9: Target mRNA production
        a[9] = k_txn[op_state] * alpha_target

        # --- Translation ---
        # Reaction 10: Mce3R translation
        a[10] = k_translation * mce3r_mRNA
        # Reaction 11: Target translation
        a[11] = k_translation * target_mRNA

        # --- mRNA decay ---
        # Reaction 12: Mce3R mRNA decay
        a[12] = gamma_mRNA * mce3r_mRNA
        # Reaction 13: Target mRNA decay
        a[13] = gamma_mRNA * target_mRNA

        # --- Protein decay ---
        # Reaction 14: Mce3R protein decay
        a[14] = gamma_protein * mce3r_protein
        # Reaction 15: Target protein decay
        a[15] = gamma_protein * target_protein

        # --- Burst placeholders (future use) ---
        a[16] = 0.0
        a[17] = 0.0

        # Total propensity
        a_total = 0.0
        for i in range(18):
            a_total += a[i]

        if a_total <= 0.0:
            t = t_max
            break

        # Draw time to next reaction (exponential waiting time)
        r1 = np.random.random()
        dt = -np.log(r1) / a_total
        t += dt

        if t >= t_max:
            break

        # Select reaction via linear search
        r2 = np.random.random() * a_total
        cumsum = 0.0
        reaction = -1
        for i in range(18):
            cumsum += a[i]
            if cumsum >= r2:
                reaction = i
                break
        if reaction == -1:
            # Fallback: pick reaction with largest propensity
            max_a = -1.0
            for i in range(18):
                if a[i] > max_a:
                    max_a = a[i]
                    reaction = i

        # Execute reaction (update state)
        if reaction == 0:
            op_state = 1   # U→S
        elif reaction == 1:
            op_state = 0   # S→U
        elif reaction == 2:
            op_state = 2   # U→W
        elif reaction == 3:
            op_state = 0   # W→U
        elif reaction == 4:
            op_state = 3   # W→D
        elif reaction == 5:
            op_state = 2   # D→W
        elif reaction == 6:
            op_state = 3   # S→D
        elif reaction == 7:
            op_state = 1   # D→S
        elif reaction == 8:
            mce3r_mRNA += 1
        elif reaction == 9:
            target_mRNA += 1
        elif reaction == 10:
            mce3r_protein += 1
        elif reaction == 11:
            target_protein += 1
        elif reaction == 12:
            if mce3r_mRNA > 0:
                mce3r_mRNA -= 1
        elif reaction == 13:
            if target_mRNA > 0:
                target_mRNA -= 1
        elif reaction == 14:
            if mce3r_protein > 0:
                mce3r_protein -= 1
        elif reaction == 15:
            if target_protein > 0:
                target_protein -= 1
        # reactions 16, 17: no-op (placeholders)

    return (mce3r_protein, target_protein, mce3r_mRNA, target_mRNA, op_state)


@numba.njit
def simulate_two_species_cell_with_trace(
    k_txn,
    n_bound,
    k_on,
    k_off_strong,
    k_off_weak,
    k_translation,
    gamma_mRNA,
    gamma_protein,
    nM_per_molecule,
    alpha_mce3r,
    alpha_target,
    t_max,
    t_burn_in,
    trace_interval,
    trace_max_points,
    seed
):
    """
    Run one Gillespie SSA trajectory for a two-species cell, recording traces.

    Same as simulate_two_species_cell but also records time-series at trace_interval.

    Returns
    -------
    (final_mce3r_protein, final_target_protein,
     final_mce3r_mRNA, final_target_mRNA, final_op_state,
     trace_times, trace_mce3r, trace_target, n_trace)
    """
    np.random.seed(seed)

    op_state       = 0
    mce3r_mRNA     = 1
    mce3r_protein  = 50
    target_mRNA    = 1
    target_protein = 50

    t = 0.0

    # Trace arrays (pre-allocated)
    trace_times  = np.zeros(trace_max_points, dtype=np.float64)
    trace_mce3r  = np.zeros(trace_max_points, dtype=np.float64)
    trace_target = np.zeros(trace_max_points, dtype=np.float64)
    n_trace = 0
    next_trace_time = t_burn_in

    a = np.zeros(18, dtype=np.float64)

    while t < t_max:
        nb = n_bound[op_state]
        free_mce3r = mce3r_protein - nb
        if free_mce3r < 0:
            free_mce3r = 0
        free_mce3r_nM = free_mce3r * nM_per_molecule

        a[0] = k_on * free_mce3r_nM if op_state == 0 else 0.0
        a[1] = k_off_strong if op_state == 1 else 0.0
        a[2] = k_on * free_mce3r_nM if op_state == 0 else 0.0
        a[3] = k_off_weak if op_state == 2 else 0.0
        a[4] = k_on * free_mce3r_nM if op_state == 2 else 0.0
        a[5] = k_off_strong if op_state == 3 else 0.0
        a[6] = k_on * free_mce3r_nM if op_state == 1 else 0.0
        a[7] = k_off_weak if op_state == 3 else 0.0
        a[8] = k_txn[op_state] * alpha_mce3r
        a[9] = k_txn[op_state] * alpha_target
        a[10] = k_translation * mce3r_mRNA
        a[11] = k_translation * target_mRNA
        a[12] = gamma_mRNA * mce3r_mRNA
        a[13] = gamma_mRNA * target_mRNA
        a[14] = gamma_protein * mce3r_protein
        a[15] = gamma_protein * target_protein
        a[16] = 0.0
        a[17] = 0.0

        a_total = 0.0
        for i in range(18):
            a_total += a[i]

        if a_total <= 0.0:
            t = t_max
            break

        r1 = np.random.random()
        dt = -np.log(r1) / a_total
        t += dt

        if t >= t_max:
            break

        # Record trace after burn-in
        if t >= t_burn_in:
            while next_trace_time <= t and n_trace < trace_max_points:
                trace_times[n_trace] = next_trace_time
                trace_mce3r[n_trace] = float(mce3r_protein)
                trace_target[n_trace] = float(target_protein)
                n_trace += 1
                next_trace_time += trace_interval

        r2 = np.random.random() * a_total
        cumsum = 0.0
        reaction = -1
        for i in range(18):
            cumsum += a[i]
            if cumsum >= r2:
                reaction = i
                break
        if reaction == -1:
            max_a = -1.0
            for i in range(18):
                if a[i] > max_a:
                    max_a = a[i]
                    reaction = i

        if reaction == 0:
            op_state = 1
        elif reaction == 1:
            op_state = 0
        elif reaction == 2:
            op_state = 2
        elif reaction == 3:
            op_state = 0
        elif reaction == 4:
            op_state = 3
        elif reaction == 5:
            op_state = 2
        elif reaction == 6:
            op_state = 3
        elif reaction == 7:
            op_state = 1
        elif reaction == 8:
            mce3r_mRNA += 1
        elif reaction == 9:
            target_mRNA += 1
        elif reaction == 10:
            mce3r_protein += 1
        elif reaction == 11:
            target_protein += 1
        elif reaction == 12:
            if mce3r_mRNA > 0:
                mce3r_mRNA -= 1
        elif reaction == 13:
            if target_mRNA > 0:
                target_mRNA -= 1
        elif reaction == 14:
            if mce3r_protein > 0:
                mce3r_protein -= 1
        elif reaction == 15:
            if target_protein > 0:
                target_protein -= 1

    # Flush remaining trace points up to t_max
    while next_trace_time <= t_max and n_trace < trace_max_points:
        trace_times[n_trace] = next_trace_time
        trace_mce3r[n_trace] = float(mce3r_protein)
        trace_target[n_trace] = float(target_protein)
        n_trace += 1
        next_trace_time += trace_interval

    return (
        mce3r_protein, target_protein, mce3r_mRNA, target_mRNA, op_state,
        trace_times, trace_mce3r, trace_target, n_trace
    )


def run_two_species_population(
    model,
    n_cells: int,
    master_seed: int,
    environment: str = 'baseline'
) -> dict:
    """
    Run n_cells two-species simulations for a given model and environment.

    Parameters
    ----------
    model : TwoSpeciesModel
        Pre-built model with all parameters.
    n_cells : int
        Number of cells to simulate.
    master_seed : int
        Master random seed; cell i uses seed = master_seed + i.
    environment : str
        One of 'baseline', 'cholesterol', 'acidic_pH', 'host_like'.
        Environmental modifications are applied to model arrays before simulation.

    Returns
    -------
    dict with:
        mce3r_proteins : int64[n_cells]
        target_proteins : int64[n_cells]
        mce3r_mRNAs : int64[n_cells]
        target_mRNAs : int64[n_cells]
        op_states : int64[n_cells]
    """
    from phase6_environmental.environmental_signals import apply_environment

    arrays = model.get_numba_arrays()
    arrays = apply_environment(arrays, environment)

    k_txn           = arrays['k_txn']
    n_bound         = arrays['n_bound']
    k_on_val        = arrays['k_on']
    k_off_s         = arrays['k_off_strong']
    k_off_w         = arrays['k_off_weak']
    k_translation   = arrays['k_translation']
    gamma_mRNA      = arrays['gamma_mRNA']
    gamma_protein   = arrays['gamma_protein']
    nM_per_molecule = arrays['nM_per_molecule']
    alpha_mce3r     = arrays['alpha_mce3r']
    alpha_target    = arrays['alpha_target']

    t_max     = np.float64(PARAMS['t_max'])
    t_burn_in = np.float64(PARAMS['t_burn_in'])

    mce3r_proteins  = np.zeros(n_cells, dtype=np.int64)
    target_proteins = np.zeros(n_cells, dtype=np.int64)
    mce3r_mRNAs     = np.zeros(n_cells, dtype=np.int64)
    target_mRNAs    = np.zeros(n_cells, dtype=np.int64)
    op_states       = np.zeros(n_cells, dtype=np.int64)

    for i in range(n_cells):
        cell_seed = master_seed + i
        mp, tp, mm, tm, os_ = simulate_two_species_cell(
            k_txn, n_bound, k_on_val, k_off_s, k_off_w,
            k_translation, gamma_mRNA, gamma_protein,
            nM_per_molecule, alpha_mce3r, alpha_target,
            t_max, t_burn_in, cell_seed
        )
        mce3r_proteins[i]  = mp
        target_proteins[i] = tp
        mce3r_mRNAs[i]     = mm
        target_mRNAs[i]    = tm
        op_states[i]       = os_

    return {
        'mce3r_proteins':  mce3r_proteins,
        'target_proteins': target_proteins,
        'mce3r_mRNAs':     mce3r_mRNAs,
        'target_mRNAs':    target_mRNAs,
        'op_states':       op_states,
    }


# ---------------------------------------------------------------------------
# Self-tests
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    from phase6_environmental.two_species_model import TwoSpeciesModel

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

    def scientific(name, cond, msg=""):
        global n_pass
        if cond:
            print(f"  SCIENTIFIC PASS: {name} {msg}")
            n_pass += 1
        else:
            print(f"  SCIENTIFIC WARN: {name} {msg}")

    print("=== two_species_gillespie.py self-tests ===")

    model = TwoSpeciesModel()
    arrays = model.get_numba_arrays()

    # Test 1: Numba compiles (small run: 10 cells, short time)
    print("  Compiling Numba JIT (first call — may take ~30 sec)...")
    result = run_two_species_population(model, n_cells=10, master_seed=42,
                                        environment='baseline')
    sanity("Numba compiles and runs 10 cells", True,
           "— JIT compilation succeeded")

    # Test 2: All species counts >= 0
    sanity("mce3r_proteins >= 0",
           np.all(result['mce3r_proteins'] >= 0),
           f"— min={result['mce3r_proteins'].min()}")
    sanity("target_proteins >= 0",
           np.all(result['target_proteins'] >= 0),
           f"— min={result['target_proteins'].min()}")
    sanity("mce3r_mRNAs >= 0",
           np.all(result['mce3r_mRNAs'] >= 0),
           f"— min={result['mce3r_mRNAs'].min()}")
    sanity("target_mRNAs >= 0",
           np.all(result['target_mRNAs'] >= 0),
           f"— min={result['target_mRNAs'].min()}")

    # Test 3: op_state in {0,1,2,3}
    sanity("op_states in {0,1,2,3}",
           np.all((result['op_states'] >= 0) & (result['op_states'] <= 3)),
           f"— unique states: {np.unique(result['op_states'])}")

    # Test 4: Mean mce3r_protein > 0 (autoregulation doesn't collapse to zero)
    mean_mce3r = np.mean(result['mce3r_proteins'])
    sanity("Mean mce3r_protein > 0",
           mean_mce3r > 0,
           f"— mean={mean_mce3r:.1f}")

    # Test 5: Trace-recording variant compiles and works
    print("  Testing trace-recording variant...")
    arr = model.get_numba_arrays()
    (mp, tp, mm, tm, os_, tt, tmce3r, ttgt, n_tr) = simulate_two_species_cell_with_trace(
        arr['k_txn'], arr['n_bound'], arr['k_on'], arr['k_off_strong'], arr['k_off_weak'],
        arr['k_translation'], arr['gamma_mRNA'], arr['gamma_protein'],
        arr['nM_per_molecule'], arr['alpha_mce3r'], arr['alpha_target'],
        np.float64(1000.0), np.float64(500.0),
        np.float64(10.0), np.int64(100), np.int64(99)
    )
    sanity("Trace variant compiles and runs", True)
    sanity("Trace has points", n_tr > 0, f"— n_trace={n_tr}")

    # Summary
    print(f"\n{'='*50}")
    print(f"SANITY: {n_pass} pass, {n_fail} fail")
    if n_fail > 0:
        print("STOPPING: Sanity check failed.")
        sys.exit(1)
    else:
        print("All checks passed.")
        sys.exit(0)
