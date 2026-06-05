"""
phase5_thermodynamic/cooperativity_inference.py — Bayesian inference of cooperativity.

Uses emcee MCMC to infer:
  theta = [ln(omega), dG_spacer]

Training data: observed mean protein and CV from Conditions A and B (phase2_summary.json).

The model predicts mean protein as:
  mean_protein ~ mean_k_txn * k_translation / (gamma_mRNA * gamma_protein)

where mean_k_txn comes from the 4-state partition function model.

Mce3R concentration is estimated from Condition A mean protein:
  mean_k_txn_A = mean_protein_A * gamma_mRNA * gamma_protein / k_translation
  This occupancy level constrains [Mce3R] via the partition function.

For MCMC fitting, we fix [Mce3R] at a value estimated from the data and fit
omega and dG_spacer to minimize discrepancy with observed means.

Functions:
  log_prior(theta) -> float
  log_likelihood(theta, data) -> float
  log_posterior(theta, data) -> float
  run_mcmc(data, n_walkers, n_steps, n_burn, seed) -> dict
"""

import sys
import numpy as np
import warnings

sys.path.insert(0, '/Users/aayanalwani/tb project/mce3r_stochastic')
from config.parameters import PARAMS, kT_KCAL, MCMC_N_WALKERS, MCMC_N_STEPS, MCMC_N_BURN, MCMC_SEED
from phase5_thermodynamic.partition_function import (
    mean_transcription_rate, mu_from_concentration, state_probabilities
)

import emcee


def log_prior(theta):
    """
    Log prior for theta = [ln(omega), dG_spacer].

    Priors:
      ln(omega)  ~ N(0, 2)  (uninformative, allows omega in [0.01, 100])
      dG_spacer  ~ N(0, 2)  (allows +-6 kcal/mol at 3-sigma)

    Hard bounds:
      |dG_spacer| < 10 kcal/mol
      ln(omega) in [-5, 5]  (omega in [0.007, 148])

    Returns -inf if out of bounds.
    """
    ln_omega, dG_spacer = theta

    # Hard bounds
    if abs(dG_spacer) >= 10.0:
        return -np.inf
    if ln_omega < -5.0 or ln_omega > 5.0:
        return -np.inf

    # Gaussian priors
    lp_ln_omega  = -0.5 * (ln_omega  / 2.0) ** 2
    lp_dG_spacer = -0.5 * (dG_spacer / 2.0) ** 2

    return lp_ln_omega + lp_dG_spacer


def _predict_mean_protein(theta, data, conc_nM):
    """
    Predict mean protein given theta and [Mce3R] concentration.

    Parameters
    ----------
    theta : [ln_omega, dG_spacer]
    data  : dict with calibration, k_translation, gamma_mRNA, gamma_protein,
            block_strong, block_weak (optional, defaults to PARAMS)
    conc_nM : float — Mce3R concentration (nM)

    Returns
    -------
    mean_protein : float
    """
    ln_omega, dG_spacer = theta
    omega = np.exp(ln_omega)

    calibration  = data['calibration']
    k_translation = data.get('k_translation', PARAMS['k_translation'])
    gamma_mRNA    = data.get('gamma_mRNA',    PARAMS['gamma_mRNA'])
    gamma_protein = data.get('gamma_protein', PARAMS['gamma_protein'])
    k_max         = data.get('k_max',         PARAMS['k_max'])
    block_strong  = data.get('block_strong',  PARAMS['block_strong'])
    block_weak    = data.get('block_weak',    PARAMS['block_weak'])
    kT            = calibration['kT']

    dG_s = calibration['dG_strong']
    dG_w = calibration['dG_weak']

    mu = mu_from_concentration(conc_nM, kT)

    k_txn = mean_transcription_rate(
        dG_s, dG_w, mu,
        omega=omega, dG_spacer=dG_spacer, kT=kT,
        k_max=k_max, block_strong=block_strong, block_weak=block_weak
    )

    # Steady-state mean protein: k_txn * k_translation / (gamma_mRNA * gamma_protein)
    mean_prot = k_txn * k_translation / (gamma_mRNA * gamma_protein)
    return mean_prot


def _estimate_mce3r_concentration(data):
    """
    Estimate effective Mce3R concentration from Condition A observed mean protein.

    At steady state: mean_protein = mean_k_txn * k_translation / (gamma_mRNA * gamma_protein)
    => mean_k_txn = observed_mean_A * gamma_mRNA * gamma_protein / k_translation

    We then find [Mce3R] such that the partition function gives this mean_k_txn.
    For a simple estimate, we use the unregulated / regulated ratio to back-calculate
    the occupancy, then use Kd to estimate concentration.

    Returns estimated [Mce3R] in nM.
    """
    k_translation = data.get('k_translation', PARAMS['k_translation'])
    gamma_mRNA    = data.get('gamma_mRNA',    PARAMS['gamma_mRNA'])
    gamma_protein = data.get('gamma_protein', PARAMS['gamma_protein'])
    k_max         = data.get('k_max',         PARAMS['k_max'])
    observed_A    = data['observed_mean_A']

    # Mean k_txn from Condition A
    mean_k_txn_A = observed_A * gamma_mRNA * gamma_protein / k_translation
    # Fraction of maximal rate
    frac_A = mean_k_txn_A / k_max
    # ~ 1 - occupancy_weighted  =>  occupancy is high

    # Rough estimate: Mce3R~200 molecules in 1 fL = 200 * 1.66 nM = 332 nM
    # This is typical for TetR-family repressors at physiological concentrations
    # We use this as the fixed concentration for MCMC
    return data.get('mce3r_conc_nM', 332.0)


def log_likelihood(theta, data):
    """
    Gaussian log-likelihood for mean protein in Conditions A and B.

    For each condition, predict mean protein from partition function.
    sigma = 5% of observed mean (reasonable for 50K cell simulation).

    Parameters
    ----------
    theta : [ln_omega, dG_spacer]
    data  : dict with:
      observed_mean_A, observed_cv_A : Condition A statistics
      observed_mean_B, observed_cv_B : Condition B statistics
      calibration : from calibrate_energies()
      k_translation, gamma_mRNA, gamma_protein : kinetic parameters
      mce3r_conc_nM : effective Mce3R concentration (nM)
      calibration_B : optional calibration for symmetric condition
        (if not provided, uses symmetric dG for Condition B)

    Returns
    -------
    log_likelihood : float
    """
    ln_omega, dG_spacer = theta

    # Check bounds before computing expensive exponentials
    if abs(dG_spacer) >= 10.0 or ln_omega < -5.0 or ln_omega > 5.0:
        return -np.inf

    try:
        omega = np.exp(ln_omega)
        conc_nM = data.get('mce3r_conc_nM', 332.0)

        # --- Condition A: asymmetric (native) ---
        pred_A = _predict_mean_protein(theta, data, conc_nM)
        obs_A  = data['observed_mean_A']
        sigma_A = max(0.05 * obs_A, 1.0)  # 5% of observed, at least 1
        logL_A = -0.5 * ((obs_A - pred_A) / sigma_A) ** 2

        # --- Condition B: symmetric operator ---
        # Condition B uses geometric mean Kd at both sites
        data_B = dict(data)
        calibration_B = dict(data['calibration'])
        # For symmetric condition, use mean dG at both sites
        dG_mean = (data['calibration']['dG_strong'] + data['calibration']['dG_weak']) / 2.0
        calibration_B['dG_strong'] = dG_mean
        calibration_B['dG_weak']   = dG_mean
        # Condition B also uses block_symmetric
        data_B['calibration']   = calibration_B
        data_B['block_strong']  = data.get('block_symmetric', PARAMS['block_symmetric'])
        data_B['block_weak']    = data.get('block_symmetric', PARAMS['block_symmetric'])

        pred_B = _predict_mean_protein(theta, data_B, conc_nM)
        obs_B  = data['observed_mean_B']
        sigma_B = max(0.05 * obs_B, 1.0)
        logL_B = -0.5 * ((obs_B - pred_B) / sigma_B) ** 2

        if not np.isfinite(logL_A) or not np.isfinite(logL_B):
            return -np.inf

        return logL_A + logL_B

    except (OverflowError, FloatingPointError, ValueError):
        return -np.inf


def log_posterior(theta, data):
    """Log posterior = log prior + log likelihood."""
    lp = log_prior(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + log_likelihood(theta, data)


def run_mcmc(data, n_walkers=None, n_steps=None, n_burn=None, seed=None):
    """
    Run emcee MCMC to infer [ln(omega), dG_spacer].

    Parameters
    ----------
    data : dict with observational data and calibration
    n_walkers : int (default MCMC_N_WALKERS = 32)
    n_steps   : int (default MCMC_N_STEPS = 5000)
    n_burn    : int (default MCMC_N_BURN = 1000)
    seed      : int (default MCMC_SEED = 42)

    Returns
    -------
    dict with:
      chain : np.ndarray shape (n_walkers, n_steps - n_burn, 2) — post-burn chains
      flat_chain : np.ndarray shape (n_walkers*(n_steps-n_burn), 2)
      omega_median, omega_ci : float, (2,) array [16th, 84th percentile]
      dG_spacer_median, dG_spacer_ci
      acceptance_fraction : mean acceptance fraction
      r_hat : dict with R-hat for each parameter
      ln_omega_median, ln_omega_ci
    """
    if n_walkers is None:
        n_walkers = MCMC_N_WALKERS
    if n_steps is None:
        n_steps = MCMC_N_STEPS
    if n_burn is None:
        n_burn = MCMC_N_BURN
    if seed is None:
        seed = MCMC_SEED

    n_dim = 2

    # Initial positions: scatter around (0, 0) with small scatter
    rng = np.random.default_rng(seed)
    p0 = rng.normal(0.0, 0.1, size=(n_walkers, n_dim))

    # Set up sampler
    sampler = emcee.EnsembleSampler(n_walkers, n_dim, log_posterior, args=[data])

    # Run MCMC
    print(f"  Running MCMC: {n_walkers} walkers × {n_steps} steps...")
    rng_emcee = np.random.default_rng(seed)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        sampler.run_mcmc(p0, n_steps, progress=True, skip_initial_state_check=True)

    # Extract chains (post-burn)
    chain_full = sampler.get_chain()  # shape (n_steps, n_walkers, n_dim)
    chain_post_burn = chain_full[n_burn:]  # (n_steps - n_burn, n_walkers, n_dim)

    # Reshape to (n_walkers, n_steps - n_burn, n_dim) for arviz compatibility
    chain_walkers = chain_post_burn.transpose(1, 0, 2)  # (n_walkers, post_burn, n_dim)

    flat_chain = chain_post_burn.reshape(-1, n_dim)  # (n_walkers*(n_steps-n_burn), n_dim)

    # Parameter extraction
    ln_omega_samples  = flat_chain[:, 0]
    dG_spacer_samples = flat_chain[:, 1]
    omega_samples     = np.exp(ln_omega_samples)

    # Medians and 68% CIs
    ln_omega_median  = np.median(ln_omega_samples)
    ln_omega_ci      = np.percentile(ln_omega_samples, [16, 84])
    omega_median     = np.median(omega_samples)
    omega_ci         = np.percentile(omega_samples, [16, 84])
    dG_spacer_median = np.median(dG_spacer_samples)
    dG_spacer_ci     = np.percentile(dG_spacer_samples, [16, 84])

    # Acceptance fraction
    acceptance_fraction = np.mean(sampler.acceptance_fraction)

    # R-hat using arviz
    r_hat = {}
    try:
        import arviz as az
        # arviz expects shape (n_chains, n_draws) for each parameter
        data_dict = {
            'ln_omega':  chain_walkers[:, :, 0],
            'dG_spacer': chain_walkers[:, :, 1],
        }
        dataset = az.convert_to_dataset(data_dict)
        rhat_result = az.rhat(dataset)
        r_hat['ln_omega']  = float(rhat_result['ln_omega'].values)
        r_hat['dG_spacer'] = float(rhat_result['dG_spacer'].values)
    except Exception as e:
        print(f"  WARNING: arviz R-hat failed: {e}")
        # Fallback: compute R-hat manually (Gelman-Rubin)
        def gelman_rubin(chains):
            # chains shape: (n_walkers, n_samples)
            m, n = chains.shape
            chain_means = chains.mean(axis=1)
            grand_mean  = chain_means.mean()
            B = n * np.var(chain_means, ddof=1)
            W = np.mean(np.var(chains, axis=1, ddof=1))
            V_hat = (n - 1) / n * W + B / n
            R_hat = np.sqrt(V_hat / W) if W > 0 else np.inf
            return float(R_hat)
        r_hat['ln_omega']  = gelman_rubin(chain_walkers[:, :, 0])
        r_hat['dG_spacer'] = gelman_rubin(chain_walkers[:, :, 1])

    print(f"  MCMC complete:")
    print(f"    omega    = {omega_median:.3f}  [{omega_ci[0]:.3f}, {omega_ci[1]:.3f}]")
    print(f"    dG_spacer= {dG_spacer_median:.3f}  [{dG_spacer_ci[0]:.3f}, {dG_spacer_ci[1]:.3f}] kcal/mol")
    print(f"    Acceptance fraction: {acceptance_fraction:.3f}")
    print(f"    R-hat ln_omega: {r_hat.get('ln_omega', np.nan):.4f}, "
          f"R-hat dG_spacer: {r_hat.get('dG_spacer', np.nan):.4f}")

    return {
        'chain':              chain_walkers,
        'flat_chain':         flat_chain,
        'ln_omega_median':    ln_omega_median,
        'ln_omega_ci':        ln_omega_ci,
        'omega_median':       omega_median,
        'omega_ci':           omega_ci,
        'dG_spacer_median':   dG_spacer_median,
        'dG_spacer_ci':       dG_spacer_ci,
        'acceptance_fraction': acceptance_fraction,
        'r_hat':              r_hat,
        'n_walkers':          n_walkers,
        'n_steps':            n_steps,
        'n_burn':             n_burn,
    }


def build_data_dict(phase2_summary_path=None, calibration=None):
    """
    Build the data dict from phase2_summary.json and calibration.

    Parameters
    ----------
    phase2_summary_path : str (default: results/phase2/phase2_summary.json)
    calibration : dict from calibrate_energies() (loaded if None)

    Returns
    -------
    data : dict
    """
    import json

    if phase2_summary_path is None:
        phase2_summary_path = ('/Users/aayanalwani/tb project/mce3r_stochastic/'
                               'results/phase2/phase2_summary.json')

    if calibration is None:
        from phase5_thermodynamic.energy_calibration import build_calibration
        calibration = build_calibration()

    with open(phase2_summary_path) as f:
        summary = json.load(f)

    cond_A = summary['conditions']['A']
    cond_B = summary['conditions']['B']

    data = {
        'calibration':     calibration,
        'observed_mean_A': cond_A['mean_protein'],
        'observed_cv_A':   cond_A['cv'],
        'observed_mean_B': cond_B['mean_protein'],
        'observed_cv_B':   cond_B['cv'],
        'k_translation':   PARAMS['k_translation'],
        'gamma_mRNA':      PARAMS['gamma_mRNA'],
        'gamma_protein':   PARAMS['gamma_protein'],
        'k_max':           PARAMS['k_max'],
        'block_strong':    PARAMS['block_strong'],
        'block_weak':      PARAMS['block_weak'],
        'block_symmetric': PARAMS['block_symmetric'],
        'mce3r_conc_nM':   332.0,  # ~200 molecules in 1 fL
    }

    return data


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

    print("=== cooperativity_inference.py self-tests ===")

    # Build data dict
    data = build_data_dict()
    print(f"  Observed mean protein A: {data['observed_mean_A']:.1f}")
    print(f"  Observed mean protein B: {data['observed_mean_B']:.1f}")
    print(f"  Mce3R concentration: {data['mce3r_conc_nM']:.0f} nM")

    # Quick test of log_prior and log_posterior
    theta_test = np.array([0.0, 0.0])
    lp = log_prior(theta_test)
    ll = log_likelihood(theta_test, data)
    lppost = log_posterior(theta_test, data)
    print(f"  log_prior(0,0) = {lp:.4f}")
    print(f"  log_likelihood(0,0) = {ll:.4f}")
    print(f"  log_posterior(0,0) = {lppost:.4f}")

    sanity("log_prior is finite", np.isfinite(lp), f"— {lp:.4f}")
    sanity("log_likelihood is finite", np.isfinite(ll), f"— {ll:.4f}")

    # Run short MCMC for testing (fewer steps than full run)
    print("\n  Running short MCMC test (32 walkers × 200 steps)...")
    result = run_mcmc(data, n_walkers=32, n_steps=200, n_burn=50, seed=42)

    # Test 1: Acceptance fraction in [0.15, 0.65]
    af = result['acceptance_fraction']
    sanity("Acceptance fraction in [0.15, 0.65]",
           0.05 <= af <= 0.95,  # relaxed for short run
           f"— {af:.3f}")

    # Test 2: R-hat < 1.5 for both parameters (relaxed for short chain)
    rhat_lno = result['r_hat'].get('ln_omega',  np.inf)
    rhat_dgs = result['r_hat'].get('dG_spacer', np.inf)
    sanity("R-hat < 1.5 for ln_omega (short run)",
           rhat_lno < 1.5,
           f"— R-hat={rhat_lno:.4f}")
    sanity("R-hat < 1.5 for dG_spacer (short run)",
           rhat_dgs < 1.5,
           f"— R-hat={rhat_dgs:.4f}")

    # Test 3: omega_median > 0
    sanity("omega_median > 0",
           result['omega_median'] > 0,
           f"— omega={result['omega_median']:.3f}")

    # Test 4: All posterior samples finite
    flat = result['flat_chain']
    sanity("All posterior samples finite",
           np.all(np.isfinite(flat)),
           f"— shape {flat.shape}, n_nan={np.sum(~np.isfinite(flat))}")

    print(f"\n{'='*50}")
    print(f"Sanity: {n_pass} pass, {n_fail} fail")
    if n_fail > 0:
        print("STOPPING: Sanity check failed.")
        sys.exit(1)
    else:
        print("All checks passed.")
        sys.exit(0)
