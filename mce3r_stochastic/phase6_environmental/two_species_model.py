"""
phase6_environmental/two_species_model.py — Two-species Mce3R circuit model.

Models the divergent transcription unit:
    Mce3R (autorepressor) <-- [Operator] --> Target (downstream gene)

Both genes are controlled by the SAME operator state (divergent transcription).
Mce3R protein feeds back into its own binding propensity (TRUE autoregulation).

Species:
    1. mce3r_mRNA
    2. mce3r_protein  (the repressor — binds operator)
    3. target_mRNA
    4. target_protein
    5. op_state        (0=both empty, 1=strong, 2=weak, 3=both)
"""

import sys
import numpy as np

sys.path.insert(0, '/Users/aayanalwani/tb project/mce3r_stochastic')
from config.parameters import PARAMS


class TwoSpeciesModel:
    """
    Two-species Mce3R circuit model extending OperatorModel logic.

    Computes all rate arrays for the two-species Gillespie engine.

    Parameters
    ----------
    Kd_strong : float, optional
        Dissociation constant (nM) for strong operator site. Default: PARAMS['Kd_strong'].
    Kd_weak : float, optional
        Dissociation constant (nM) for weak operator site. Default: PARAMS['Kd_weak'].
    k_on : float, optional
        Binding rate constant (nM^-1 min^-1). Default: PARAMS['k_on'].
    block_strong : float, optional
        Transcription block fraction when strong site occupied. Default: PARAMS['block_strong'].
    block_weak : float, optional
        Transcription block fraction when weak site occupied. Default: PARAMS['block_weak'].
    k_max : float, optional
        Maximum transcription rate (mRNA/min). Default: PARAMS['k_max'].
    alpha_mce3r : float
        Relative transcription strength for Mce3R direction. Default: 1.0.
    alpha_target : float
        Relative transcription strength for Target direction. Default: 1.0.
    """

    def __init__(
        self,
        Kd_strong=None,
        Kd_weak=None,
        k_on=None,
        block_strong=None,
        block_weak=None,
        k_max=None,
        alpha_mce3r=1.0,
        alpha_target=1.0,
    ):
        self.Kd_strong   = Kd_strong   if Kd_strong   is not None else PARAMS['Kd_strong']
        self.Kd_weak     = Kd_weak     if Kd_weak     is not None else PARAMS['Kd_weak']
        self.k_on        = k_on        if k_on        is not None else PARAMS['k_on']
        self.block_strong = block_strong if block_strong is not None else PARAMS['block_strong']
        self.block_weak  = block_weak  if block_weak  is not None else PARAMS['block_weak']
        self.k_max       = k_max       if k_max       is not None else PARAMS['k_max']
        self.alpha_mce3r = float(alpha_mce3r)
        self.alpha_target = float(alpha_target)

        # Derived unbinding rates
        self.k_off_strong = self.Kd_strong * self.k_on
        self.k_off_weak   = self.Kd_weak   * self.k_on

        # Transcription rate per operator state (length-4 float64)
        k = self.k_max
        bs = self.block_strong
        bw = self.block_weak
        self.k_txn = np.array([
            k,                          # State 0: both empty
            k * (1.0 - bs),             # State 1: strong occupied
            k * (1.0 - bw),             # State 2: weak occupied
            k * (1.0 - bs) * (1.0 - bw),  # State 3: both occupied
        ], dtype=np.float64)

        # Number of Mce3R molecules bound per state
        self.n_bound = np.array([0, 1, 1, 2], dtype=np.int64)

    def get_numba_arrays(self) -> dict:
        """
        Return all arrays needed by two_species_gillespie.simulate_two_species_cell().

        Returns
        -------
        dict with keys:
            k_txn          : float64[4]  — base transcription rate per state
            n_bound        : int64[4]    — Mce3R molecules bound per state
            k_on           : float64     — binding rate (nM^-1 min^-1)
            k_off_strong   : float64     — unbinding from strong site (min^-1)
            k_off_weak     : float64     — unbinding from weak site (min^-1)
            k_translation  : float64     — translation rate (protein/mRNA/min)
            gamma_mRNA     : float64     — mRNA degradation rate (min^-1)
            gamma_protein  : float64     — protein degradation rate (min^-1)
            nM_per_molecule: float64     — nM per single molecule in cell
            alpha_mce3r    : float64     — transcription strength for Mce3R direction
            alpha_target   : float64     — transcription strength for Target direction
        """
        return {
            'k_txn':           self.k_txn.copy(),
            'n_bound':         self.n_bound.copy(),
            'k_on':            np.float64(self.k_on),
            'k_off_strong':    np.float64(self.k_off_strong),
            'k_off_weak':      np.float64(self.k_off_weak),
            'k_translation':   np.float64(PARAMS['k_translation']),
            'gamma_mRNA':      np.float64(PARAMS['gamma_mRNA']),
            'gamma_protein':   np.float64(PARAMS['gamma_protein']),
            'nM_per_molecule': np.float64(PARAMS['nM_per_molecule']),
            'alpha_mce3r':     np.float64(self.alpha_mce3r),
            'alpha_target':    np.float64(self.alpha_target),
        }


# ---------------------------------------------------------------------------
# Self-tests
# ---------------------------------------------------------------------------
if __name__ == "__main__":
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

    print("=== two_species_model.py self-tests ===")

    # Default model (asymmetric, wild-type)
    model = TwoSpeciesModel()
    arrays = model.get_numba_arrays()

    # Test 1: Arrays have correct numpy dtypes
    sanity("k_txn dtype float64",
           arrays['k_txn'].dtype == np.float64,
           f"— got {arrays['k_txn'].dtype}")
    sanity("n_bound dtype int64",
           arrays['n_bound'].dtype == np.int64,
           f"— got {arrays['n_bound'].dtype}")
    sanity("k_on is float64 scalar",
           isinstance(arrays['k_on'], (np.float64, float)),
           f"— got type {type(arrays['k_on'])}")
    sanity("gamma_protein is float64 scalar",
           isinstance(arrays['gamma_protein'], (np.float64, float)),
           f"— got type {type(arrays['gamma_protein'])}")
    sanity("alpha_mce3r is float64 scalar",
           isinstance(arrays['alpha_mce3r'], (np.float64, float)),
           f"— got {arrays['alpha_mce3r']}")
    sanity("alpha_target is float64 scalar",
           isinstance(arrays['alpha_target'], (np.float64, float)),
           f"— got {arrays['alpha_target']}")

    # Test 2: alpha_mce3r = alpha_target = 1.0 → both use same k_txn array (equal scaling)
    model_sym_alpha = TwoSpeciesModel(alpha_mce3r=1.0, alpha_target=1.0)
    arr_sym = model_sym_alpha.get_numba_arrays()
    # Both directions scale k_txn by the same alpha (1.0), so scaled rates are identical
    mce3r_scaled = arr_sym['k_txn'] * arr_sym['alpha_mce3r']
    target_scaled = arr_sym['k_txn'] * arr_sym['alpha_target']
    sanity("alpha_mce3r=alpha_target=1.0 → same effective k_txn for both",
           np.allclose(mce3r_scaled, target_scaled),
           f"— mce3r={mce3r_scaled}, target={target_scaled}")

    # Test 3: n_bound = [0, 1, 1, 2]
    sanity("n_bound = [0, 1, 1, 2]",
           np.array_equal(arrays['n_bound'], np.array([0, 1, 1, 2])),
           f"— got {arrays['n_bound']}")

    # Test 4: k_off derived correctly
    expected_k_off_s = model.Kd_strong * model.k_on
    expected_k_off_w = model.Kd_weak * model.k_on
    sanity("k_off_strong = Kd_strong * k_on",
           abs(arrays['k_off_strong'] - expected_k_off_s) < 1e-10,
           f"— got {arrays['k_off_strong']:.6f}, expected {expected_k_off_s:.6f}")
    sanity("k_off_weak = Kd_weak * k_on",
           abs(arrays['k_off_weak'] - expected_k_off_w) < 1e-10,
           f"— got {arrays['k_off_weak']:.6f}, expected {expected_k_off_w:.6f}")

    # Test 5: k_txn monotonic (state0 > state2 > state1 > state3)
    k = arrays['k_txn']
    sanity("k_txn monotonic",
           k[0] > k[2] > k[1] > k[3],
           f"— {k[0]:.4f} > {k[2]:.4f} > {k[1]:.4f} > {k[3]:.4f}")

    # Test 6: All rates non-negative
    sanity("All rates non-negative",
           all([
               np.all(arrays['k_txn'] >= 0),
               arrays['k_on'] >= 0,
               arrays['k_off_strong'] >= 0,
               arrays['k_off_weak'] >= 0,
               arrays['k_translation'] >= 0,
               arrays['gamma_mRNA'] >= 0,
               arrays['gamma_protein'] >= 0,
           ]),
           "— checked k_txn, k_on, k_offs, k_translation, gammas")

    # Summary
    print(f"\n{'='*50}")
    print(f"SANITY: {n_pass} pass, {n_fail} fail")
    if n_fail > 0:
        print("STOPPING: Sanity check failed.")
        sys.exit(1)
    else:
        print("All checks passed.")
        sys.exit(0)
