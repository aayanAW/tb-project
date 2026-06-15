"""
phase2_simulation/operator_model.py — OperatorModel class for Mce3R operator states.

4-state operator model:
  State 0: both sites empty        (n_bound=0)
  State 1: strong site occupied    (n_bound=1)
  State 2: weak site occupied      (n_bound=1)
  State 3: both sites occupied     (n_bound=2)

Pre-computes arrays suitable for passing into Numba @njit functions.
"""

import sys
import numpy as np

sys.path.insert(0, '/Users/aayanalwani/tb project/mce3r_stochastic')
from config.parameters import PARAMS


class OperatorModel:
    """Pre-compute operator transition arrays for the Gillespie engine."""

    def __init__(self, Kd_strong=None, Kd_weak=None, k_on=None,
                 block_strong=None, block_weak=None, k_max=None):
        """
        Build operator model. All parameters default to PARAMS if not given.
        This allows condition-specific overrides (e.g. Condition C sets Kd_weak=1e12).
        """
        self.Kd_strong = Kd_strong if Kd_strong is not None else PARAMS['Kd_strong']
        self.Kd_weak = Kd_weak if Kd_weak is not None else PARAMS['Kd_weak']
        self.k_on = k_on if k_on is not None else PARAMS['k_on']
        self.block_strong = block_strong if block_strong is not None else PARAMS['block_strong']
        self.block_weak = block_weak if block_weak is not None else PARAMS['block_weak']
        self.k_max = k_max if k_max is not None else PARAMS['k_max']

        # Derived k_off values
        self.k_off_strong = self.Kd_strong * self.k_on
        self.k_off_weak = self.Kd_weak * self.k_on

        # n_bound per state: state 0->0, 1->1, 2->1, 3->2
        self.n_bound_lookup = np.array([0, 1, 1, 2], dtype=np.int64)

    def get_transcription_rates(self):
        """Return length-4 float64 array of transcription rates per operator state."""
        k = self.k_max
        bs = self.block_strong
        bw = self.block_weak
        return np.array([
            k,                          # State 0: both empty
            k * (1.0 - bs),             # State 1: strong occupied
            k * (1.0 - bw),             # State 2: weak occupied
            k * (1.0 - bs) * (1.0 - bw) # State 3: both occupied
        ], dtype=np.float64)

    def get_binding_rate_constants(self):
        """Return dict with k_on, k_off_strong, k_off_weak."""
        return {
            'k_on': self.k_on,
            'k_off_strong': self.k_off_strong,
            'k_off_weak': self.k_off_weak,
        }

    def get_numba_arrays(self):
        """
        Return all arrays needed by the Numba Gillespie engine.

        Returns:
            dict with keys:
                k_txn: float64[4] — transcription rate per state
                n_bound: int64[4] — number of repressor molecules bound per state
                k_on: float64 scalar
                k_off_strong: float64 scalar
                k_off_weak: float64 scalar
        """
        return {
            'k_txn': self.get_transcription_rates(),
            'n_bound': self.n_bound_lookup.copy(),
            'k_on': np.float64(self.k_on),
            'k_off_strong': np.float64(self.k_off_strong),
            'k_off_weak': np.float64(self.k_off_weak),
        }


# ---------------------------------------------------------------------------
# Self-tests
# ---------------------------------------------------------------------------
if __name__ == "__main__":
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

    print("=== operator_model.py self-tests ===")

    # --- Default model (Condition A: asymmetric) ---
    model_a = OperatorModel()
    rates_a = model_a.get_transcription_rates()
    binding_a = model_a.get_binding_rate_constants()
    numba_a = model_a.get_numba_arrays()

    # SANITY: rates array shape
    sanity("k_txn length", len(rates_a) == 4, f"— got {len(rates_a)}")

    # SANITY: all rates positive
    sanity("k_txn all positive", np.all(rates_a > 0),
           f"— min={rates_a.min():.6f}")

    # SANITY: rates monotonic (state0 > state2 > state1 > state3)
    sanity("k_txn monotonic", rates_a[0] > rates_a[2] > rates_a[1] > rates_a[3],
           f"— {rates_a[0]:.4f} > {rates_a[2]:.4f} > {rates_a[1]:.4f} > {rates_a[3]:.4f}")

    # SANITY: k_off matches Kd * k_on
    sanity("k_off_strong matches",
           abs(binding_a['k_off_strong'] - PARAMS['Kd_strong'] * PARAMS['k_on']) < 1e-10,
           f"— {binding_a['k_off_strong']:.6f}")
    sanity("k_off_weak matches",
           abs(binding_a['k_off_weak'] - PARAMS['Kd_weak'] * PARAMS['k_on']) < 1e-10,
           f"— {binding_a['k_off_weak']:.6f}")

    # SANITY: n_bound_lookup
    sanity("n_bound_lookup",
           np.array_equal(model_a.n_bound_lookup, np.array([0, 1, 1, 2])),
           f"— {model_a.n_bound_lookup}")

    # SANITY: numba arrays have correct types
    sanity("numba k_txn dtype", numba_a['k_txn'].dtype == np.float64)
    sanity("numba n_bound dtype", numba_a['n_bound'].dtype == np.int64)

    # --- Condition C: single-site (weak disabled) ---
    model_c = OperatorModel(Kd_weak=PARAMS['Kd_weak_disabled'])
    rates_c = model_c.get_transcription_rates()
    sanity("Condition C: weak k_off very large",
           model_c.k_off_weak > 1e9,
           f"— k_off_weak={model_c.k_off_weak:.2e}")

    # --- Condition D: no regulation ---
    model_d = OperatorModel(k_on=PARAMS['k_on_disabled'])
    rates_d = model_d.get_transcription_rates()
    sanity("Condition D: k_on=0", model_d.k_on == 0.0)
    # All txn rates should equal k_max for no-reg, but the model still
    # computes from block fractions. Caller should set k_txn=[k_max]*4.
    # This is handled in run_conditions.py.

    # --- Condition B: symmetric (both sites use geometric mean Kd and block) ---
    model_b = OperatorModel(Kd_strong=PARAMS['Kd_symmetric'],
                            Kd_weak=PARAMS['Kd_symmetric'],
                            block_strong=PARAMS['block_symmetric'],
                            block_weak=PARAMS['block_symmetric'])
    rates_b = model_b.get_transcription_rates()
    scientific("Condition B: symmetric rates",
               abs(rates_b[1] - rates_b[2]) < 1e-10,
               f"— state1={rates_b[1]:.4f}, state2={rates_b[2]:.4f}")

    # Summary
    print(f"\n{'='*50}")
    print(f"SANITY: {n_pass} pass, {n_fail} fail | SCIENTIFIC: {n_warn} warn")
    if n_fail > 0:
        print("STOPPING: Sanity check failed.")
        sys.exit(1)
    else:
        print("All checks passed.")
        sys.exit(0)
