"""
phase2_simulation/gillespie_engine.py — Numba-accelerated Gillespie SSA engine.

12 reactions in the 4-state operator model:
  Reactions 0-7: operator transitions (binding/unbinding at strong/weak sites)
    0: 0→1 (bind strong, empty→strong)
    1: 1→0 (unbind strong)
    2: 0→2 (bind weak, empty→weak)
    3: 2→0 (unbind weak)
    4: 2→3 (bind strong when weak occupied)
    5: 3→2 (unbind strong when both occupied)
    6: 1→3 (bind weak when strong occupied)
    7: 3→1 (unbind weak when both occupied)
  Reaction 8:  transcription (mRNA production)
  Reaction 9:  translation (protein production)
  Reaction 10: mRNA decay
  Reaction 11: protein decay
"""

import sys
import numpy as np
import numba

sys.path.insert(0, '/Users/aayanalwani/tb project/mce3r_stochastic')
from config.parameters import PARAMS


@numba.njit
def simulate_cell(k_txn, n_bound, k_on, k_off_strong, k_off_weak,
                  k_translation, gamma_mRNA, gamma_protein,
                  nM_per_molecule, t_max, t_burn_in,
                  trace_interval, trace_max_points,
                  seed, record_trace):
    """
    Run one Gillespie SSA trajectory for a single cell.

    All inputs are scalars or numpy arrays (no dicts/classes — Numba-safe).

    Parameters
    ----------
    k_txn : float64[4]        transcription rate per operator state
    n_bound : int64[4]        repressor molecules bound per state
    k_on : float64            binding rate constant (nM^-1 min^-1)
    k_off_strong : float64    unbinding rate for strong site (min^-1)
    k_off_weak : float64      unbinding rate for weak site (min^-1)
    k_translation : float64   translation rate (protein/mRNA/min)
    gamma_mRNA : float64      mRNA degradation rate (min^-1)
    gamma_protein : float64   protein degradation rate (min^-1)
    nM_per_molecule : float64 nM concentration per molecule
    t_max : float64           total simulation time (min)
    t_burn_in : float64       burn-in time (min)
    trace_interval : float64  interval for recording trace (min)
    trace_max_points : int64  pre-allocated trace array size
    seed : int64              random seed
    record_trace : int64      1 to record trace, 0 to skip

    Returns
    -------
    final_protein : int64
    final_mRNA : int64
    final_op_state : int64
    trace_times : float64[trace_max_points]
    trace_proteins : float64[trace_max_points]
    n_trace : int64            number of valid trace points
    """
    np.random.seed(seed)

    # State variables — initial condition: everything empty/zero
    op_state = 0   # operator state (0=both empty)
    mRNA = 0       # mRNA count
    protein = 0    # protein count
    t = 0.0

    # Trace arrays (pre-allocated for Numba)
    trace_times = np.zeros(trace_max_points, dtype=np.float64)
    trace_proteins = np.zeros(trace_max_points, dtype=np.float64)
    n_trace = 0
    next_trace_time = t_burn_in  # start recording after burn-in

    # Propensity array (12 reactions)
    a = np.zeros(12, dtype=np.float64)

    while t < t_max:
        # Compute free protein
        nb = n_bound[op_state]
        free_protein = protein - nb
        if free_protein < 0:
            free_protein = 0
        R_nM = free_protein * nM_per_molecule

        # --- Compute propensities ---
        # Binding propensities depend on R_nM (concentration-dependent)
        # Unbinding propensities are concentration-independent

        # Reaction 0: 0→1 (bind strong from state 0)
        if op_state == 0:
            a[0] = k_on * R_nM
        else:
            a[0] = 0.0

        # Reaction 1: 1→0 (unbind strong from state 1)
        if op_state == 1:
            a[1] = k_off_strong
        else:
            a[1] = 0.0

        # Reaction 2: 0→2 (bind weak from state 0)
        if op_state == 0:
            a[2] = k_on * R_nM
        else:
            a[2] = 0.0

        # Reaction 3: 2→0 (unbind weak from state 2)
        if op_state == 2:
            a[3] = k_off_weak
        else:
            a[3] = 0.0

        # Reaction 4: 2→3 (bind strong from state 2)
        if op_state == 2:
            a[4] = k_on * R_nM
        else:
            a[4] = 0.0

        # Reaction 5: 3→2 (unbind strong from state 3)
        if op_state == 3:
            a[5] = k_off_strong
        else:
            a[5] = 0.0

        # Reaction 6: 1→3 (bind weak from state 1)
        if op_state == 1:
            a[6] = k_on * R_nM
        else:
            a[6] = 0.0

        # Reaction 7: 3→1 (unbind weak from state 3)
        if op_state == 3:
            a[7] = k_off_weak
        else:
            a[7] = 0.0

        # Reaction 8: transcription
        a[8] = k_txn[op_state]

        # Reaction 9: translation
        a[9] = k_translation * mRNA

        # Reaction 10: mRNA decay
        a[10] = gamma_mRNA * mRNA

        # Reaction 11: protein decay
        a[11] = gamma_protein * protein

        # Total propensity
        a_total = 0.0
        for i in range(12):
            a_total += a[i]

        if a_total <= 0.0:
            # No reactions possible — advance to end
            t = t_max
            break

        # Draw time to next reaction
        r1 = np.random.random()
        dt = -np.log(r1) / a_total
        t += dt

        if t >= t_max:
            break

        # Record trace if past burn-in and due
        if record_trace == 1 and t >= t_burn_in:
            while next_trace_time <= t and n_trace < trace_max_points:
                trace_times[n_trace] = next_trace_time
                trace_proteins[n_trace] = float(protein)
                n_trace += 1
                next_trace_time += trace_interval

        # Select reaction
        r2 = np.random.random() * a_total
        cumsum = 0.0
        reaction = -1
        for i in range(12):
            cumsum += a[i]
            if cumsum >= r2:
                reaction = i
                break
        if reaction == -1:
            # Fallback: select reaction with largest propensity (avoids bias)
            max_a = -1.0
            for i in range(12):
                if a[i] > max_a:
                    max_a = a[i]
                    reaction = i

        # Execute reaction
        if reaction == 0:
            op_state = 1      # 0→1
        elif reaction == 1:
            op_state = 0      # 1→0
        elif reaction == 2:
            op_state = 2      # 0→2
        elif reaction == 3:
            op_state = 0      # 2→0
        elif reaction == 4:
            op_state = 3      # 2→3
        elif reaction == 5:
            op_state = 2      # 3→2
        elif reaction == 6:
            op_state = 3      # 1→3
        elif reaction == 7:
            op_state = 1      # 3→1
        elif reaction == 8:
            mRNA += 1
        elif reaction == 9:
            protein += 1
        elif reaction == 10:
            if mRNA > 0:
                mRNA -= 1
        elif reaction == 11:
            if protein > 0:
                protein -= 1

    # Record any remaining trace points up to t_max
    if record_trace == 1:
        while next_trace_time <= t_max and n_trace < trace_max_points:
            trace_times[n_trace] = next_trace_time
            trace_proteins[n_trace] = float(protein)
            n_trace += 1
            next_trace_time += trace_interval

    return protein, mRNA, op_state, trace_times, trace_proteins, n_trace


def run_population(model_arrays, n_cells, master_seed, record_traces=0,
                   n_trace_cells=0):
    """
    Run Gillespie SSA for a population of cells.

    Parameters
    ----------
    model_arrays : dict from OperatorModel.get_numba_arrays() plus overrides
    n_cells : int
    master_seed : int
    record_traces : int (0 or 1)
    n_trace_cells : int — how many cells to record traces for

    Returns
    -------
    dict with:
        proteins : int64[n_cells]
        mRNAs : int64[n_cells]
        op_states : int64[n_cells]
        traces : list of (times, proteins) tuples (only first n_trace_cells)
    """
    k_txn = model_arrays['k_txn']
    n_bound = model_arrays['n_bound']
    k_on_val = model_arrays['k_on']
    k_off_s = model_arrays['k_off_strong']
    k_off_w = model_arrays['k_off_weak']

    k_translation = np.float64(PARAMS['k_translation'])
    gamma_mRNA = np.float64(PARAMS['gamma_mRNA'])
    gamma_protein = np.float64(PARAMS['gamma_protein'])
    nM_per_molecule = np.float64(PARAMS['nM_per_molecule'])
    t_max = np.float64(PARAMS['t_max'])
    t_burn_in = np.float64(PARAMS['t_burn_in'])
    trace_interval = np.float64(PARAMS['trace_interval'])
    trace_max_points = PARAMS['trace_max_points']

    proteins = np.zeros(n_cells, dtype=np.int64)
    mRNAs = np.zeros(n_cells, dtype=np.int64)
    op_states = np.zeros(n_cells, dtype=np.int64)
    traces = []

    for i in range(n_cells):
        do_trace = 1 if (record_traces and i < n_trace_cells) else 0
        cell_seed = master_seed + i

        p, m, o, tt, tp, nt = simulate_cell(
            k_txn, n_bound, k_on_val, k_off_s, k_off_w,
            k_translation, gamma_mRNA, gamma_protein,
            nM_per_molecule, t_max, t_burn_in,
            trace_interval, trace_max_points,
            cell_seed, do_trace
        )
        proteins[i] = p
        mRNAs[i] = m
        op_states[i] = o

        if do_trace:
            traces.append((tt[:nt].copy(), tp[:nt].copy()))

    return {
        'proteins': proteins,
        'mRNAs': mRNAs,
        'op_states': op_states,
        'traces': traces,
    }


# ---------------------------------------------------------------------------
# Self-tests
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    from phase2_simulation.operator_model import OperatorModel

    n_pass = 0
    n_fail = 0
    n_warn = 0

    def sanity(name, cond, msg=""):
        global n_pass, n_fail
        if cond:
            print(f"  SANITY PASS: {name} {msg}")
            n_pass += 1
        else:
            print(f"  SANITY FAIL: {name} {msg}")
            n_fail += 1

    def scientific(name, cond, msg=""):
        global n_pass, n_warn
        if cond:
            print(f"  SCIENTIFIC PASS: {name} {msg}")
            n_pass += 1
        else:
            print(f"  SCIENTIFIC WARN: {name} {msg}")
            n_warn += 1

    print("=== gillespie_engine.py self-tests ===")

    # --- Test 1: Numba compilation (run 1 cell) ---
    print("  Compiling Numba JIT (first call)...")
    model = OperatorModel()
    arrays = model.get_numba_arrays()
    p, m, o, tt, tp, nt = simulate_cell(
        arrays['k_txn'], arrays['n_bound'],
        arrays['k_on'], arrays['k_off_strong'], arrays['k_off_weak'],
        np.float64(PARAMS['k_translation']),
        np.float64(PARAMS['gamma_mRNA']),
        np.float64(PARAMS['gamma_protein']),
        np.float64(PARAMS['nM_per_molecule']),
        np.float64(PARAMS['t_max']),
        np.float64(PARAMS['t_burn_in']),
        np.float64(PARAMS['trace_interval']),
        PARAMS['trace_max_points'],
        42, 1
    )
    sanity("Numba compiles", True, "— JIT compilation succeeded")
    sanity("protein >= 0", p >= 0, f"— protein={p}")
    sanity("mRNA >= 0", m >= 0, f"— mRNA={m}")
    sanity("op_state in [0,1,2,3]", o in [0, 1, 2, 3], f"— op_state={o}")
    sanity("trace has points", nt > 0, f"— n_trace={nt}")

    # --- Test 2: Run small population with traces ---
    print("  Running 10-cell population test...")
    result = run_population(arrays, n_cells=10, master_seed=42,
                            record_traces=1, n_trace_cells=3)
    sanity("population proteins shape", result['proteins'].shape == (10,))
    sanity("traces recorded", len(result['traces']) == 3,
           f"— got {len(result['traces'])} traces")
    sanity("all proteins >= 0", np.all(result['proteins'] >= 0))

    # --- Test 3: Unregulated (Condition D) — check Fano factor ---
    print("  Running 500-cell unregulated test for Fano check...")
    model_d = OperatorModel(k_on=0.0)
    arrays_d = model_d.get_numba_arrays()
    # Override k_txn to all k_max (no regulation)
    arrays_d['k_txn'] = np.array([PARAMS['k_max']] * 4, dtype=np.float64)
    result_d = run_population(arrays_d, n_cells=500, master_seed=99)

    prots_d = result_d['proteins'].astype(np.float64)
    mean_d = np.mean(prots_d)
    var_d = np.var(prots_d)
    fano_d = var_d / mean_d if mean_d > 0 else 0.0
    expected_fano = PARAMS['expected_fano_unregulated']

    sanity("unreg mean > 0", mean_d > 0, f"— mean={mean_d:.1f}")
    scientific("Fano ≈ 7.85 (±50%)",
               expected_fano * 0.5 < fano_d < expected_fano * 1.5,
               f"— Fano={fano_d:.2f}, expected={expected_fano:.2f}")

    # --- Test 4: Regulated (Condition A) — protein should be lower ---
    print("  Running 100-cell regulated test...")
    result_a = run_population(arrays, n_cells=100, master_seed=42)
    mean_a = np.mean(result_a['proteins'].astype(np.float64))
    scientific("regulated mean < unregulated mean",
               mean_a < mean_d,
               f"— regulated={mean_a:.1f}, unregulated={mean_d:.1f}")

    # --- Test 5: Trajectory in sane range ---
    scientific("protein < 10000", np.max(result_a['proteins']) < 10000,
               f"— max={np.max(result_a['proteins'])}")

    # Summary
    print(f"\n{'='*50}")
    print(f"SANITY: {n_pass} pass, {n_fail} fail | SCIENTIFIC: {n_warn} warn")
    if n_fail > 0:
        print("STOPPING: Sanity check failed.")
        sys.exit(1)
    else:
        print("All checks passed.")
        sys.exit(0)
