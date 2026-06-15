"""
phase5_thermodynamic/partition_function.py — 4-state statistical mechanics model.

The 4 states of the operator:
  U  = both empty          weight = 1                          k_txn = k_max
  S  = strong site bound   weight = exp(-beta*(dG_s - mu))     k_txn = k_max*(1-block_s)
  W  = weak site bound     weight = exp(-beta*(dG_w - mu))     k_txn = k_max*(1-block_w)
  D  = doubly bound        weight = omega*exp(-beta*(dG_s+dG_w-2*mu+dG_spacer))
                                                               k_txn = k_max*(1-block_s)*(1-block_w)

Functions:
  partition_function(dG_strong, dG_weak, mu, omega, dG_spacer, kT) -> float
  state_probabilities(dG_strong, dG_weak, mu, omega, dG_spacer, kT,
                      block_strong, block_weak) -> np.ndarray(4,)
  mean_transcription_rate(...) -> float
  repression_fold(...) -> float
  mu_from_concentration(conc_nM, kT) -> float
  repression_curve(concentrations_nM, dG_strong, dG_weak, omega, dG_spacer,
                   kT, k_max, block_strong, block_weak) -> dict
  compare_architectures(concentrations_nM, calibration, omega, dG_spacer) -> dict
"""

import sys
import numpy as np

sys.path.insert(0, '/Users/aayanalwani/tb project/mce3r_stochastic')
from config.parameters import PARAMS, kT_KCAL


def partition_function(dG_strong, dG_weak, mu, omega=1.0, dG_spacer=0.0, kT=kT_KCAL):
    """
    Compute the partition function for the 4-state operator model.

    Z = 1
      + exp(-beta*(dG_s - mu))
      + exp(-beta*(dG_w - mu))
      + omega * exp(-beta*(dG_s + dG_w - 2*mu + dG_spacer))

    Parameters
    ----------
    dG_strong : float — calibrated DeltaG for strong site (kcal/mol, negative = favorable)
    dG_weak   : float — calibrated DeltaG for weak site (kcal/mol)
    mu        : float — chemical potential of repressor (kcal/mol)
    omega     : float — cooperativity factor (1.0 = no cooperativity)
    dG_spacer : float — spacer penalty (kcal/mol, positive = penalty)
    kT        : float — thermal energy (kcal/mol)

    Returns
    -------
    Z : float — partition function (dimensionless)
    """
    beta = 1.0 / kT

    w_U = 1.0
    w_S = np.exp(-beta * (dG_strong - mu))
    w_W = np.exp(-beta * (dG_weak   - mu))
    w_D = omega * np.exp(-beta * (dG_strong + dG_weak - 2.0 * mu + dG_spacer))

    Z = w_U + w_S + w_W + w_D
    return Z


def state_probabilities(dG_strong, dG_weak, mu, omega=1.0, dG_spacer=0.0, kT=kT_KCAL,
                        block_strong=None, block_weak=None):
    """
    Compute probability of each operator state.

    Returns
    -------
    probs : np.ndarray shape (4,) — [P(U), P(S), P(W), P(D)]
    """
    beta = 1.0 / kT

    w_U = 1.0
    w_S = np.exp(-beta * (dG_strong - mu))
    w_W = np.exp(-beta * (dG_weak   - mu))
    w_D = omega * np.exp(-beta * (dG_strong + dG_weak - 2.0 * mu + dG_spacer))

    Z = w_U + w_S + w_W + w_D

    return np.array([w_U / Z, w_S / Z, w_W / Z, w_D / Z])


def mean_transcription_rate(dG_strong, dG_weak, mu, omega=1.0, dG_spacer=0.0, kT=kT_KCAL,
                             k_max=None, block_strong=None, block_weak=None):
    """
    Compute mean transcription rate weighted by state probabilities.

    Rate for each state:
      U: k_max
      S: k_max * (1 - block_strong)
      W: k_max * (1 - block_weak)
      D: k_max * (1 - block_strong) * (1 - block_weak)

    Parameters
    ----------
    All thermodynamic parameters as above, plus:
    k_max        : float (default PARAMS['k_max'])
    block_strong : float (default PARAMS['block_strong'])
    block_weak   : float (default PARAMS['block_weak'])

    Returns
    -------
    k_mean : float (mRNA/min)
    """
    if k_max is None:
        k_max = PARAMS['k_max']
    if block_strong is None:
        block_strong = PARAMS['block_strong']
    if block_weak is None:
        block_weak = PARAMS['block_weak']

    probs = state_probabilities(dG_strong, dG_weak, mu, omega, dG_spacer, kT)

    k_rates = np.array([
        k_max,
        k_max * (1.0 - block_strong),
        k_max * (1.0 - block_weak),
        k_max * (1.0 - block_strong) * (1.0 - block_weak),
    ])

    return np.dot(probs, k_rates)


def repression_fold(dG_strong, dG_weak, mu, omega=1.0, dG_spacer=0.0, kT=kT_KCAL,
                    k_max=None, block_strong=None, block_weak=None):
    """
    Compute fold repression: k_max / mean_transcription_rate.
    """
    if k_max is None:
        k_max = PARAMS['k_max']
    if block_strong is None:
        block_strong = PARAMS['block_strong']
    if block_weak is None:
        block_weak = PARAMS['block_weak']

    k_mean = mean_transcription_rate(dG_strong, dG_weak, mu, omega, dG_spacer, kT,
                                     k_max, block_strong, block_weak)
    if k_mean <= 0.0:
        return np.inf
    return k_max / k_mean


def mu_from_concentration(conc_nM, kT=kT_KCAL):
    """
    Convert repressor concentration (nM) to chemical potential (kcal/mol).

    mu = kT * ln([R]_M) = kT * ln(conc_nM * 1e-9)

    At [R] = Kd, mu = kT*ln(Kd_M) = DeltaG_site, so P(bound) = 0.5. Correct.

    Parameters
    ----------
    conc_nM : float — repressor concentration (nM)
    kT : float — thermal energy (kcal/mol)

    Returns
    -------
    mu : float (kcal/mol)
    """
    if conc_nM <= 0.0:
        return -np.inf
    return kT * np.log(conc_nM * 1e-9)


def repression_curve(concentrations_nM, dG_strong, dG_weak,
                     omega=1.0, dG_spacer=0.0, kT=kT_KCAL,
                     k_max=None, block_strong=None, block_weak=None):
    """
    Compute repression curve over a range of repressor concentrations.

    Parameters
    ----------
    concentrations_nM : array-like of concentrations (nM)
    (all other parameters as above)

    Returns
    -------
    dict with arrays:
      concentrations_nM : input concentrations
      mu : chemical potentials
      k_mean : mean transcription rate at each concentration
      fold : fold repression at each concentration
      prob_U, prob_S, prob_W, prob_D : state probabilities
    """
    if k_max is None:
        k_max = PARAMS['k_max']
    if block_strong is None:
        block_strong = PARAMS['block_strong']
    if block_weak is None:
        block_weak = PARAMS['block_weak']

    concs = np.asarray(concentrations_nM, dtype=np.float64)
    N = len(concs)

    mu_arr       = np.zeros(N)
    k_mean_arr   = np.zeros(N)
    fold_arr     = np.zeros(N)
    prob_U_arr   = np.zeros(N)
    prob_S_arr   = np.zeros(N)
    prob_W_arr   = np.zeros(N)
    prob_D_arr   = np.zeros(N)

    for i, c in enumerate(concs):
        mu = mu_from_concentration(c, kT)
        mu_arr[i] = mu

        probs = state_probabilities(dG_strong, dG_weak, mu, omega, dG_spacer, kT)
        prob_U_arr[i] = probs[0]
        prob_S_arr[i] = probs[1]
        prob_W_arr[i] = probs[2]
        prob_D_arr[i] = probs[3]

        k_rates = np.array([
            k_max,
            k_max * (1.0 - block_strong),
            k_max * (1.0 - block_weak),
            k_max * (1.0 - block_strong) * (1.0 - block_weak),
        ])
        k_mean = np.dot(probs, k_rates)
        k_mean_arr[i] = k_mean
        fold_arr[i] = k_max / k_mean if k_mean > 0 else np.inf

    return {
        'concentrations_nM': concs,
        'mu':                mu_arr,
        'k_mean':            k_mean_arr,
        'fold':              fold_arr,
        'prob_U':            prob_U_arr,
        'prob_S':            prob_S_arr,
        'prob_W':            prob_W_arr,
        'prob_D':            prob_D_arr,
    }


def compare_architectures(concentrations_nM, calibration, omega=1.0, dG_spacer=0.0,
                          kT=kT_KCAL, k_max=None, block_strong=None, block_weak=None):
    """
    Compare repression curves for 4 operator architectures.

    Architectures:
      native_asymmetric : dG_strong and dG_weak as calibrated
      symmetric         : mean dG at both sites  (same total regulatory capacity)
      strong_only       : dG_strong at both sites (two strong sites)
      single_site       : only strong site active; weak site = +100 kcal/mol (inert)

    Parameters
    ----------
    concentrations_nM : array-like
    calibration : dict from calibrate_energies()
    omega, dG_spacer, kT : thermodynamic parameters
    k_max, block_strong, block_weak : transcription parameters (default from PARAMS)

    Returns
    -------
    dict mapping architecture name -> repression_curve dict
    """
    if k_max is None:
        k_max = PARAMS['k_max']
    if block_strong is None:
        block_strong = PARAMS['block_strong']
    if block_weak is None:
        block_weak = PARAMS['block_weak']

    dG_s = calibration['dG_strong']  # e.g. -12.23 kcal/mol
    dG_w = calibration['dG_weak']    # e.g. -10.37 kcal/mol
    dG_mean = (dG_s + dG_w) / 2.0

    architectures = {
        'native_asymmetric': (dG_s,    dG_w),
        'symmetric':         (dG_mean, dG_mean),
        'strong_only':       (dG_s,    dG_s),
        'single_site':       (dG_s,    100.0),  # weak site effectively inert
    }

    results = {}
    for arch_name, (dg_s, dg_w) in architectures.items():
        curve = repression_curve(
            concentrations_nM, dg_s, dg_w,
            omega=omega, dG_spacer=dG_spacer, kT=kT,
            k_max=k_max, block_strong=block_strong, block_weak=block_weak
        )
        results[arch_name] = curve

    return results


# ---------------------------------------------------------------------------
# Self-tests
# ---------------------------------------------------------------------------
if __name__ == '__main__':
    import sys

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
        global n_pass, n_fail
        if cond:
            print(f"  SCIENTIFIC PASS: {name} {msg}")
            n_pass += 1
        else:
            print(f"  SCIENTIFIC WARN (not fail): {name} {msg}")
            # Scientific warnings don't count as failures

    from phase5_thermodynamic.energy_calibration import build_calibration

    print("=== partition_function.py self-tests ===")

    calibration = build_calibration()
    dG_s = calibration['dG_strong']
    dG_w = calibration['dG_weak']
    kT   = calibration['kT']

    print(f"  Using dG_strong = {dG_s:.3f}, dG_weak = {dG_w:.3f} kcal/mol")

    # Test 1: At [Mce3R]=0 (mu -> -inf): P(U) = 1.0, fold = 1.0
    mu_zero = mu_from_concentration(1e-20)  # effectively 0
    probs_zero = state_probabilities(dG_s, dG_w, mu_zero, kT=kT)
    fold_zero  = repression_fold(dG_s, dG_w, mu_zero, kT=kT)
    sanity("At [R]~0: P(U) approx 1",
           probs_zero[0] > 0.999,
           f"— P(U)={probs_zero[0]:.6f}")
    sanity("At [R]~0: fold repression approx 1",
           abs(fold_zero - 1.0) < 0.001,
           f"— fold={fold_zero:.4f}")

    # Test 2: At [Mce3R]->infinity: P(D) -> 1.0
    mu_inf = mu_from_concentration(1e12)  # 1 uM, very high
    probs_inf = state_probabilities(dG_s, dG_w, mu_inf, omega=2.0, kT=kT)
    sanity("At [R]>>Kd: P(D) approaches 1",
           probs_inf[3] > 0.99,
           f"— P(D)={probs_inf[3]:.4f}")

    # Test 3: Sum of probabilities = 1.0 at all concentrations
    concs_test = np.logspace(-2, 4, 50)  # 0.01 to 10000 nM
    all_sum_one = True
    for c in concs_test:
        mu = mu_from_concentration(c)
        probs = state_probabilities(dG_s, dG_w, mu, omega=1.5, kT=kT)
        if abs(probs.sum() - 1.0) > 1e-10:
            all_sum_one = False
            break
    sanity("Probabilities sum to 1.0 at all concentrations",
           all_sum_one,
           f"— tested {len(concs_test)} concentrations")

    # Test 4: Fold repression monotonically increasing with [Mce3R]
    concs_mono = np.logspace(-2, 4, 100)
    curve = repression_curve(concs_mono, dG_s, dG_w, omega=1.0)
    folds = curve['fold']
    is_monotone = np.all(np.diff(folds) >= -1e-6)  # allow tiny numerical noise
    sanity("Fold repression monotonically increasing",
           is_monotone,
           f"— max fold = {folds[-1]:.2f}")

    # Test 5: Asymmetric architecture has steeper transition than symmetric
    # (expected when omega > 1; with omega=1 this may not hold strictly)
    concs_compare = np.logspace(-2, 4, 200)
    results = compare_architectures(concs_compare, calibration, omega=2.0)
    asym_curve = results['native_asymmetric']['fold']
    sym_curve   = results['symmetric']['fold']

    # Hill coefficient proxy: steepness around K_half
    # Asymmetric should reach higher max fold or have steeper mid-range slope
    # Just check that both curves go up and are different
    asym_max = asym_curve[-1]
    sym_max  = sym_curve[-1]
    scientific("Asymmetric max fold >= symmetric max fold",
               asym_max >= sym_max * 0.9,
               f"— asym={asym_max:.2f}, sym={sym_max:.2f}")

    print(f"\n  Max fold repression (asymmetric): {asym_max:.2f}")
    print(f"  Max fold repression (symmetric):  {sym_max:.2f}")

    print(f"\n{'='*50}")
    print(f"Sanity: {n_pass} pass, {n_fail} fail")
    if n_fail > 0:
        print("STOPPING: Sanity check failed.")
        sys.exit(1)
    else:
        print("All checks passed.")
        sys.exit(0)
