# PHASE 2 ARCHITECTURE: Thermodynamic Model + Environmental Stochastic Extensions

## Extending the Mce3R Asymmetric Operator Pipeline

---

## CLEAN EXTENSION REQUIREMENT

**This document specifies NEW code to add to an EXISTING, COMPLETE project.**

- Project directory: `/Users/aayanalwani/tb project/mce3r_stochastic/`
- The existing pipeline (Phases 1-4) is COMPLETE and PASSING. Do NOT modify existing files unless explicitly stated.
- Read `HANDOFF.md` first for full context on what exists.
- Read `PROJECT_PROPOSAL.md` for the scientific motivation.
- All new code integrates with the existing `config/parameters.py` and reuses existing modules where specified.

---

## TABLE OF CONTENTS

1. [What Already Exists](#1-what-already-exists)
2. [What This Build Adds](#2-what-this-build-adds)
3. [New Directory Structure](#3-new-directory-structure)
4. [Environment & New Dependencies](#4-environment--new-dependencies)
5. [Phase 5: Thermodynamic Model (Aim 2)](#5-phase-5-thermodynamic-model)
6. [Phase 6: Environmental Stochastic Extensions (Aim 4)](#6-phase-6-environmental-stochastic-extensions)
7. [Phase 7: Extended Figures](#7-phase-7-extended-figures)
8. [Updated Orchestrator](#8-updated-orchestrator)
9. [Parameter Additions](#9-parameter-additions)
10. [Self-Testing Requirements](#10-self-testing-requirements)
11. [Critical Design Decisions](#11-critical-design-decisions)
12. [Scientific Predictions to Validate](#12-scientific-predictions-to-validate)

---

## 1. WHAT ALREADY EXISTS

The following is COMPLETE and WORKING. Do not rebuild it.

### Existing Phases

| Phase | Purpose | Status | Key Output |
|-------|---------|--------|------------|
| Phase 1 | MEME/FIMO motif discovery + conservation | COMPLETE | `predicted_sites.csv`, `conservation_status.csv` |
| Phase 2 | Gillespie SSA: 4 conditions + sweep | COMPLETE | `condition_A/B/C/D.npz`, `condition_E_sweep.npz` |
| Phase 3 | Statistical validation | COMPLETE | `noise_metrics.csv`, `statistical_tests.csv`, `bootstrap_results.csv` |
| Phase 4 | 7 publication figures | COMPLETE | `results/figures/fig1-8_*.png` |

### Existing Key Results (use these, do not recompute)

| Condition | Architecture | Kd (nM) | Mean Protein | CV | Fano |
|-----------|-------------|---------|-------------|------|------|
| A | Asymmetric (wild-type) | 2.4 / 49.0 | 197.3 | 0.1866 | 6.87 |
| B | Symmetric (geometric mean) | 10.84 / 10.84 | 292.6 | 0.1567 | 7.18 |
| C | Single-site (weak disabled) | 2.4 / 1e12 | 341.8 | 0.1495 | 7.63 |
| D | No regulation | 0 / 0 | 2224.0 | 0.0591 | 7.77 |

### Existing Code You MUST Reuse

| Module | What to reuse |
|--------|--------------|
| `config/parameters.py` | All existing parameters. ADD new ones, don't modify existing. |
| `phase2_simulation/operator_model.py` | `OperatorModel` class for building 4-state models |
| `phase2_simulation/gillespie_engine.py` | `simulate_cell()` (numba) and `run_population()` for stochastic simulations |
| `phase2_simulation/run_conditions.py` | `compute_stats()` helper function |
| `phase1_pipeline/run_fimo.py` | FIMO results for new operator candidates |
| `results/phase2/condition_A.npz` | Existing simulation data (load, don't rerun) |

### Existing Design Decisions (RESPECT THESE)

1. **Symmetric control = geometric mean Kd.** `Kd_symmetric = sqrt(2.4 * 49.0) = 10.84 nM`
2. **Asymmetry sweep holds geometric mean constant.** `Kd_s = Kd_geo / sqrt(r)`, `Kd_w = Kd_geo * sqrt(r)`
3. **The simulated "protein" is an abstract downstream output**, not literal Mce3R abundance (see FULL_PROJECT_ARCHITECTURE.md Section 2.5).
4. **Numba JIT for Gillespie.** All inner-loop simulation code must use `@numba.njit`.

---

## 2. WHAT THIS BUILD ADDS

### Phase 5 — Thermodynamic Model (Aim 2 from PROJECT_PROPOSAL.md)

Build a statistical-mechanical model of the 4-state operator that:
- Converts PWM scores to binding free energies, calibrated to experimental Kd values
- Computes partition functions and occupancy probabilities as a function of Mce3R concentration
- Fits cooperativity (omega) and spacer penalty via MCMC
- Predicts repression curves for native vs hypothetical architectures
- Classifies operators as digital switches, graded repressors, or noise modulators

### Phase 6 — Environmental Stochastic Extensions (Aim 4 from PROJECT_PROPOSAL.md)

Extend the existing Gillespie engine to:
- Add environmental signals (cholesterol, acidic pH) that modulate effective Mce3R concentration
- Build a two-species circuit (Mce3R repressor + downstream target) with explicit autoregulation
- Calculate mutual information between environment and gene output
- Define and quantify persistence threshold crossing under different architectures and environments
- Compare native asymmetric vs symmetric vs single-site operators under host-like stress

### Phase 7 — Extended Figures

Generate 5-6 new publication figures for the extended results.

---

## 3. NEW DIRECTORY STRUCTURE

```
mce3r_stochastic/                          # EXISTING root
├── [all existing files unchanged]
│
├── phase5_thermodynamic/                  # NEW
│   ├── __init__.py
│   ├── thermodynamic_main.py             # Phase 5 orchestrator
│   ├── energy_calibration.py             # PWM → ΔG conversion
│   ├── partition_function.py             # 4-state Z, occupancy, repression curves
│   ├── cooperativity_inference.py        # MCMC fitting of omega and spacer penalty
│   ├── operator_classification.py        # Classify operators by regulatory mode
│   └── dna_shape.py                      # DNA structural features (optional, see notes)
│
├── phase6_environmental/                  # NEW
│   ├── __init__.py
│   ├── environmental_main.py             # Phase 6 orchestrator
│   ├── environmental_signals.py          # Cholesterol/pH signal models
│   ├── two_species_model.py              # Mce3R + downstream target circuit
│   ├── two_species_gillespie.py          # Extended Gillespie engine (numba)
│   ├── mutual_information.py             # I(input; output) estimation
│   ├── persistence_threshold.py          # Threshold definition + counting
│   └── run_environmental_conditions.py   # Run all environmental × architecture combos
│
├── phase7_extended_figures/               # NEW
│   ├── __init__.py
│   ├── extended_figures_main.py          # Phase 7 orchestrator
│   ├── fig9_repression_curves.py         # Thermodynamic repression predictions
│   ├── fig10_cooperativity.py            # MCMC posteriors + operator heatmap
│   ├── fig11_environmental_distributions.py  # Expression under host-like stress
│   ├── fig12_mutual_information.py       # MI vs Mce3R concentration
│   ├── fig13_persistence_phase.py        # Phase diagram: persistence vs asymmetry × environment
│   └── fig14_two_species_traces.py       # Two-species circuit time traces
│
├── results/
│   ├── [existing phase1-3 results unchanged]
│   ├── phase5/                           # NEW
│   │   ├── energy_parameters.json
│   │   ├── repression_curves.csv
│   │   ├── mcmc_posteriors.npz
│   │   ├── operator_classifications.csv
│   │   └── phase5_summary.json
│   ├── phase6/                           # NEW
│   │   ├── env_condition_*.npz           # Per-environment × architecture results
│   │   ├── two_species_results.npz
│   │   ├── mutual_information.csv
│   │   ├── persistence_fractions.csv
│   │   └── phase6_summary.json
│   └── extended_figures/                 # NEW
│       ├── fig9-14_*.png
│       └── phase7_summary.json
│
└── logs/
    ├── [existing logs unchanged]
    ├── phase5_thermodynamic.log          # NEW
    ├── phase6_environmental.log          # NEW
    └── phase7_extended_figures.log       # NEW
```

---

## 4. ENVIRONMENT & NEW DEPENDENCIES

Add to the existing `mce3r` conda environment:

```bash
conda activate mce3r
pip install emcee corner arviz
```

| Package | Version | Purpose |
|---------|---------|---------|
| emcee | >=3.1 | MCMC ensemble sampler for cooperativity inference |
| corner | >=2.2 | Corner plots for MCMC posteriors |
| arviz | >=0.15 | MCMC diagnostics (R-hat, ESS) |

**Already installed (reuse):** numpy, scipy, pandas, matplotlib, seaborn, scikit-learn, numba, biopython

**NOT required:** DNAshapeR (R package). DNA shape features are optional — see Section 11.

---

## 5. PHASE 5: THERMODYNAMIC MODEL

### 5.1 Dependency

Phase 5 depends on:
- `config/parameters.py` (Kd values, operator sequence, k_max, block fractions)
- `results/phase1/predicted_sites.csv` (FIMO-predicted binding sites for new operators)
- `results/phase1/meme_output/meme.txt` (PWM from MEME)

Phase 5 does NOT depend on Phase 2 simulation data. It can run in parallel with Phases 1-4 in principle, but for simplicity run it sequentially after Phase 4.

### 5.2 `energy_calibration.py`

**Purpose:** Convert PWM log-odds scores to binding free energies (ΔG) and calibrate to experimental Kd values.

**Functions:**

```python
def load_pwm_from_meme(meme_txt_path: str) -> list[np.ndarray]:
    """
    Parse meme.txt to extract position-weight matrices.
    Returns list of PWM arrays, each shape (motif_length, 4) with columns [A, C, G, T].
    Use Bio.motifs.meme module or manual parsing.
    """

def pwm_to_energy_matrix(pwm: np.ndarray, background: dict = None) -> np.ndarray:
    """
    Convert PWM frequencies to position-specific energy contributions.

    For each position i and base b:
        epsilon(i, b) = -kT * ln(pwm[i,b] / background[b])

    where kT = 0.616 kcal/mol at 37°C (310 K), R = 1.987e-3 kcal/(mol·K).

    Args:
        pwm: shape (L, 4) frequency matrix
        background: dict with keys 'A','C','G','T', default = H37Rv GC content (0.656)
                    so background = {'A': 0.172, 'C': 0.328, 'G': 0.328, 'T': 0.172}

    Returns:
        energy_matrix: shape (L, 4) in kcal/mol. Lower = more favorable binding.
    """

def score_sequence(seq: str, energy_matrix: np.ndarray) -> float:
    """
    Sum position-specific energies for a given sequence.
    Returns total ΔG_seq in kcal/mol.
    """

def calibrate_energies(
    known_strong_seq: str,
    known_weak_seq: str,
    energy_matrix: np.ndarray,
    Kd_strong: float = 2.4,   # nM
    Kd_weak: float = 49.0,    # nM
    T: float = 310.0           # Kelvin (37°C)
) -> dict:
    """
    Calibrate the energy scale so that:
        ΔG_strong = RT * ln(Kd_strong / 1e9)  [convert nM to M]
        ΔG_weak = RT * ln(Kd_weak / 1e9)

    The raw PWM energies have an arbitrary offset. We compute:
        offset = ΔG_strong_experimental - ΔG_strong_raw

    Then for any new sequence:
        ΔG_calibrated = ΔG_raw + offset

    And the predicted Kd:
        Kd_predicted = exp(ΔG_calibrated / RT) * 1e9  (in nM)

    Returns dict with:
        - offset: float (kcal/mol)
        - dG_strong: float (calibrated, kcal/mol)
        - dG_weak: float (calibrated, kcal/mol)
        - Kd_strong_predicted: float (nM, should match 2.4)
        - Kd_weak_predicted: float (nM, should match 49.0)
        - energy_matrix_calibrated: np.ndarray
        - T: float
        - RT: float
    """

def predict_Kd_for_sequence(seq: str, calibration: dict) -> float:
    """
    Given a new binding site sequence, predict its Kd in nM
    using the calibrated energy matrix.
    """
```

**Known sequences for calibration:**

The 123-bp operator is in `config/parameters.py` as `OPERATOR_SEQUENCE`. The strong site is the downstream 25 bp and the weak site is the upstream 25 bp. Extract them as:

```python
# From Panagoda 2024, the operator is:
# 5'-[upstream_weak_25bp]---[53bp_spacer]---[downstream_strong_25bp]-3'
# Positions within the 123-bp operator:
strong_site = OPERATOR_SEQUENCE[-25:]    # last 25 bp (downstream, high affinity)
weak_site = OPERATOR_SEQUENCE[:25]       # first 25 bp (upstream, low affinity)
```

**IMPORTANT:** Verify this assignment by checking that the PWM raw score for the strong site is LOWER (more favorable) than the weak site. If reversed, swap the indexing.

**Self-tests (4):**
1. Calibrated Kd_strong within 10% of 2.4 nM
2. Calibrated Kd_weak within 10% of 49.0 nM
3. Energy matrix has no NaN or Inf values
4. Background frequencies sum to 1.0

---

### 5.3 `partition_function.py`

**Purpose:** Compute the 4-state partition function, occupancy probabilities, and mean transcription rate as functions of Mce3R concentration.

**The Model:**

The promoter exists in 4 states:

| State | Strong site | Weak site | Weight in Z | Transcription rate |
|-------|------------|-----------|-------------|-------------------|
| 0 (U) | Empty | Empty | 1 | k_max |
| 1 (S) | Bound | Empty | exp(-β(ΔG_s - μ)) | k_max × (1 - block_s) |
| 2 (W) | Empty | Bound | exp(-β(ΔG_w - μ)) | k_max × (1 - block_w) |
| 3 (D) | Bound | Bound | ω × exp(-β(ΔG_s + ΔG_w - 2μ + ΔG_spacer)) | k_max × (1 - block_s)(1 - block_w) |

where:
- β = 1/(kT), with kT = 0.616 kcal/mol at 37°C
- μ = μ₀ + kT × ln([Mce3R] / [Mce3R]_ref) is the chemical potential of Mce3R
- ω is the cooperativity factor (ω > 1 = positive cooperativity, ω < 1 = anti-cooperative)
- ΔG_spacer penalizes non-optimal spacing between sites

**Functions:**

```python
def partition_function(
    dG_strong: float,       # kcal/mol (negative = favorable)
    dG_weak: float,         # kcal/mol
    mu: float,              # chemical potential (kcal/mol)
    omega: float = 1.0,     # cooperativity
    dG_spacer: float = 0.0, # spacer penalty (kcal/mol, positive = unfavorable)
    kT: float = 0.616       # kcal/mol at 37°C
) -> float:
    """Compute Z = 1 + exp(-β(ΔG_s - μ)) + exp(-β(ΔG_w - μ)) + ω·exp(-β(ΔG_s + ΔG_w - 2μ + ΔG_spacer))"""

def state_probabilities(
    dG_strong: float, dG_weak: float, mu: float,
    omega: float = 1.0, dG_spacer: float = 0.0, kT: float = 0.616
) -> np.ndarray:
    """
    Returns array of shape (4,) with P(U), P(S), P(W), P(D).
    P(state_i) = weight_i / Z
    """

def mean_transcription_rate(
    dG_strong: float, dG_weak: float, mu: float,
    k_max: float, block_strong: float, block_weak: float,
    omega: float = 1.0, dG_spacer: float = 0.0, kT: float = 0.616
) -> float:
    """
    <k_txn> = sum_i P(state_i) * k_txn_i
    """

def repression_fold(
    dG_strong: float, dG_weak: float, mu: float,
    k_max: float, block_strong: float, block_weak: float,
    omega: float = 1.0, dG_spacer: float = 0.0, kT: float = 0.616
) -> float:
    """
    Fold repression = k_max / <k_txn>
    When fully unbound, fold = 1. When fully repressed, fold → 1/((1-bs)(1-bw)).
    """

def mu_from_concentration(
    concentration_nM: float,
    mu_ref: float = 0.0,
    conc_ref_nM: float = 1.0,
    kT: float = 0.616
) -> float:
    """
    Convert Mce3R concentration (nM) to chemical potential.
    μ = μ_ref + kT × ln(concentration / conc_ref)
    """

def repression_curve(
    concentrations_nM: np.ndarray,
    dG_strong: float, dG_weak: float,
    k_max: float, block_strong: float, block_weak: float,
    omega: float = 1.0, dG_spacer: float = 0.0, kT: float = 0.616
) -> dict:
    """
    Compute mean transcription rate and fold repression at each Mce3R concentration.

    Args:
        concentrations_nM: array of Mce3R concentrations to evaluate

    Returns dict:
        - concentrations: input array
        - mean_k_txn: array of mean transcription rates
        - fold_repression: array of fold repression values
        - prob_unbound: array of P(U) at each concentration
        - prob_strong: array of P(S)
        - prob_weak: array of P(W)
        - prob_double: array of P(D)
    """

def compare_architectures(
    concentrations_nM: np.ndarray,
    calibration: dict,
    omega: float = 1.0,
    dG_spacer: float = 0.0
) -> dict:
    """
    Compute repression curves for 4 architectures:
    1. Native asymmetric (ΔG_s, ΔG_w from calibration)
    2. Symmetric (ΔG = mean of ΔG_s and ΔG_w at both sites)
    3. Strong-only (ΔG_s at both sites — NOT the correct symmetric control, but informative)
    4. Single-site (only strong site active, weak site ΔG → +∞)

    Returns dict keyed by architecture name, each containing repression_curve output.
    """
```

**Self-tests (5):**
1. At [Mce3R] = 0: P(U) = 1.0, fold_repression = 1.0
2. At [Mce3R] → ∞: P(D) → 1.0 (if ω ≥ 1)
3. Sum of state probabilities = 1.0 at all concentrations
4. Fold repression is monotonically increasing with [Mce3R]
5. Asymmetric operator has steeper transition (larger Hill-like coefficient) than symmetric at same geometric mean affinity

---

### 5.4 `cooperativity_inference.py`

**Purpose:** Fit cooperativity (ω) and spacer penalty (ΔG_spacer) using MCMC.

**Training data:** We have limited experimental data. Use these constraints:
1. The known operator reproduces Kd_strong = 2.4 nM and Kd_weak = 49 nM (from calibration)
2. The Gillespie simulation with the existing parameters produces the observed CV and mean protein values
3. Any newly discovered operators from Phase 1 FIMO hits provide additional data points

**Strategy:** Bayesian inference with informative priors.

```python
def log_prior(theta: np.ndarray) -> float:
    """
    theta = [ln(omega), dG_spacer]

    Priors:
    - ln(omega) ~ Normal(0, 2)  [centered on no cooperativity, wide]
      This means omega ~ LogNormal(0, 2), ranging from ~0.02 to ~50
    - dG_spacer ~ Normal(0, 2) kcal/mol [centered on no penalty, wide]

    Returns -inf if outside bounds:
    - omega must be > 0 (ln(omega) unconstrained)
    - |dG_spacer| < 10 kcal/mol
    """

def log_likelihood(theta: np.ndarray, data: dict) -> float:
    """
    Compute log-likelihood of observed data given parameters.

    data dict contains:
    - 'observed_mean_A': float (197.3, from existing simulation)
    - 'observed_cv_A': float (0.1866)
    - 'observed_mean_B': float (292.6)
    - 'observed_cv_B': float (0.1567)
    - 'concentrations_nM': np.ndarray (sweep of Mce3R concentrations)
    - 'calibration': dict (from energy_calibration)

    For each observed condition, compute predicted mean transcription rate
    from the thermodynamic model, convert to predicted mean protein
    (mean_protein ≈ k_txn × k_translation / (gamma_mRNA × gamma_protein)),
    and compute Gaussian log-likelihood:
        logL += -0.5 * ((observed - predicted) / sigma)^2

    sigma estimated from bootstrap CIs of existing results.
    """

def log_posterior(theta: np.ndarray, data: dict) -> float:
    """log_prior(theta) + log_likelihood(theta, data)"""

def run_mcmc(
    data: dict,
    n_walkers: int = 32,
    n_steps: int = 5000,
    n_burn: int = 1000,
    seed: int = 42
) -> dict:
    """
    Run emcee ensemble MCMC sampler.

    Returns dict:
    - chain: np.ndarray shape (n_walkers, n_steps, 2) [ln(omega), dG_spacer]
    - flat_chain: np.ndarray shape (n_walkers * (n_steps - n_burn), 2)
    - omega_median: float
    - omega_ci: (float, float)  [16th, 84th percentile]
    - dG_spacer_median: float
    - dG_spacer_ci: (float, float)
    - acceptance_fraction: float (should be 0.2-0.5)
    - r_hat: np.ndarray (Gelman-Rubin convergence diagnostic, should be < 1.1)
    """
```

**MCMC settings:**
- 32 walkers, 5000 steps, 1000 burn-in
- Initial positions: scatter around (ln(omega)=0, dG_spacer=0) with scale 0.1
- Target acceptance fraction: 0.2–0.5 (adjust proposal scale via emcee's stretch move)
- Convergence: R-hat < 1.1, ESS > 100

**Self-tests (4):**
1. Acceptance fraction in [0.15, 0.65]
2. R-hat < 1.2 for both parameters
3. omega_median > 0
4. Posterior samples are finite (no NaN/Inf)

---

### 5.5 `operator_classification.py`

**Purpose:** Use the fitted thermodynamic model to classify operators by regulatory behavior.

```python
def classify_operator(
    dG_strong: float,
    dG_weak: float,
    omega: float,
    dG_spacer: float,
    k_max: float,
    block_strong: float,
    block_weak: float
) -> dict:
    """
    Classify an operator based on its repression curve properties.

    Compute the repression curve over [Mce3R] = 0.01 to 10,000 nM (200 log-spaced points).

    Metrics:
    1. Hill coefficient (n_H): fit fold_repression vs [Mce3R] to Hill equation
       fold = fold_max * [R]^n / (K_half^n + [R]^n)
       n_H > 2 → ultrasensitive/digital switch
       n_H ≈ 1 → graded
       n_H < 1 → subsensitive

    2. Dynamic range: fold_repression(max) / fold_repression(min)

    3. Transition width: [Mce3R] range over which fold goes from 10% to 90% of max
       Narrow → switch-like, Wide → graded

    4. Weak-site contribution:
       Compare fold_repression with both sites vs strong-site-only.
       If difference < 20% → weak site is marginal
       If difference > 50% → weak site is a cooperative enhancer
       If difference is primarily in variance not mean → noise modulator

    Returns dict:
        - classification: str ('digital_switch' | 'graded_repressor' | 'noise_modulator')
        - hill_coefficient: float
        - dynamic_range: float
        - transition_width_nM: float
        - weak_site_contribution: float (fractional increase in repression)
        - K_half: float (nM, concentration at half-maximal repression)
    """

def classify_all_operators(
    predicted_sites_csv: str,
    calibration: dict,
    omega: float,
    dG_spacer: float,
    top_n: int = 20
) -> pd.DataFrame:
    """
    Load top N FIMO-predicted operators from Phase 1.
    For each, extract strong/weak site sequences, compute ΔG, classify.

    Returns DataFrame with columns:
    rank, start, stop, dG_strong, dG_weak, Kd_strong_predicted, Kd_weak_predicted,
    classification, hill_coefficient, dynamic_range, K_half
    """
```

**Self-tests (3):**
1. Known operator classified consistently with its noise behavior
2. Hill coefficient is finite and positive
3. All top-20 operators have valid classifications

---

### 5.6 `dna_shape.py` (OPTIONAL)

**Purpose:** Compute DNA structural features that may influence binding.

**IMPORTANT:** This module is OPTIONAL. DNAshapeR is an R/Bioconductor package that requires R installation. If R is not available, skip this module entirely. The thermodynamic model works without it — ΔG_spacer from MCMC already captures the net effect.

If implemented:
- Use `subprocess` to call R with a script that runs DNAshapeR on the operator sequences
- Extract: minor groove width (MGW), propeller twist (ProT), roll, helix twist (HelT)
- Include mean values as additional features in the energy model

If NOT implemented:
- Report in the summary that DNA shape features were not included
- Note this as a future direction

**Self-tests (1):**
1. If R is available, shape features are finite. If not, module reports "skipped" gracefully.

---

### 5.7 `thermodynamic_main.py`

```python
def run_phase5() -> dict:
    """
    Orchestrate Phase 5:
    1. Load MEME PWM and calibrate energies (energy_calibration.py)
    2. Compute repression curves for 4 architectures (partition_function.py)
    3. Run MCMC for cooperativity inference (cooperativity_inference.py)
    4. Classify top-20 operators (operator_classification.py)
    5. (Optional) Compute DNA shape features (dna_shape.py)
    6. Write results to results/phase5/
    7. Write phase5_summary.json

    Returns dict with standard contract:
    {
        'status': 'COMPLETE' | 'FAILED',
        'n_sanity_pass': int,
        'n_sanity_fail': int,
        'n_science_expected': int,
        'n_science_unexpected': int,
        'n_warn': int,
        'mock_used': bool,
        'wall_time_sec': float,
        'output_files': list[str],
    }
    """
```

**Output files:**
- `results/phase5/energy_parameters.json` — Calibrated energies, ΔG values, predicted Kds
- `results/phase5/repression_curves.csv` — Repression vs [Mce3R] for each architecture
- `results/phase5/mcmc_posteriors.npz` — Full MCMC chain + flat chain
- `results/phase5/mcmc_summary.json` — omega, dG_spacer medians + CIs + diagnostics
- `results/phase5/operator_classifications.csv` — Top-20 operator properties
- `results/phase5/phase5_summary.json` — Standard summary

---

## 6. PHASE 6: ENVIRONMENTAL STOCHASTIC EXTENSIONS

### 6.1 Dependency

Phase 6 depends on:
- `config/parameters.py` (all parameters)
- `phase2_simulation/gillespie_engine.py` (reuse `simulate_cell` as template)
- `phase2_simulation/operator_model.py` (reuse `OperatorModel`)
- `results/phase5/mcmc_summary.json` (fitted omega and dG_spacer)
- `results/phase2/condition_A.npz` (for comparison, load existing)

Phase 6 MUST run after Phase 5 (needs cooperativity parameter).

### 6.2 `environmental_signals.py`

**Purpose:** Define how cholesterol and acidic pH affect Mce3R effective concentration.

**Biological basis:**
- Mce3R is a repressor. When cholesterol is present, Mce3R ligand (likely a cholesterol derivative) binds Mce3R and reduces its DNA-binding affinity → derepression.
- Acidic pH (pH 5.5 in phagosome) may affect Mce3R protein stability or folding.
- Combined stress (cholesterol + acid) represents the host macrophage environment.

```python
# Environmental conditions to simulate
ENVIRONMENTS = {
    'baseline': {
        'description': 'Glycerol, neutral pH (in vitro)',
        'mce3r_multiplier': 1.0,        # no change to effective [Mce3R]
        'noise_scale': 1.0,             # no additional noise
    },
    'cholesterol': {
        'description': 'Cholesterol as sole carbon source',
        'mce3r_multiplier': 0.3,        # cholesterol derivative reduces Mce3R binding ~3-fold
        'noise_scale': 1.0,
    },
    'acidic_pH': {
        'description': 'pH 5.5 (phagosomal)',
        'mce3r_multiplier': 0.7,        # mild reduction in Mce3R activity
        'noise_scale': 1.2,             # acid stress adds stochastic protein degradation
    },
    'host_like': {
        'description': 'Cholesterol + acidic pH (macrophage phagosome)',
        'mce3r_multiplier': 0.2,        # combined strong derepression
        'noise_scale': 1.3,             # combined stress increases noise
    },
}

def apply_environment(
    model_arrays: dict,
    environment: str
) -> dict:
    """
    Modify model arrays to reflect environmental condition.

    The mce3r_multiplier scales the effective k_on (binding rate):
        k_on_effective = k_on * mce3r_multiplier

    The noise_scale multiplies the protein degradation rate:
        gamma_protein_effective = gamma_protein * noise_scale

    This is a phenomenological approximation. The multiplier values
    are order-of-magnitude estimates based on:
    - Cholesterol reduces TetR-family repressor binding ~3-10 fold (various refs)
    - Acid stress increases protein turnover ~20-30% (Vandal 2008)

    Returns modified copy of model_arrays (does NOT mutate input).
    """
```

**IMPORTANT CAVEAT:** The exact multiplier values are educated guesses. The key scientific question is NOT "what is the exact fold-change?" but "does asymmetry amplify noise MORE under stress than under baseline?" The qualitative prediction is robust to the exact parameter values — verify this via sensitivity analysis.

**Self-tests (2):**
1. All environments produce valid model_arrays (no negative rates)
2. `host_like` has lower k_on_effective than `baseline`

---

### 6.3 `two_species_model.py`

**Purpose:** Define a two-species circuit: Mce3R (repressor) + Target (downstream gene).

**The circuit:**

```
Mce3R ──┐
        │ represses (via operator)
        ▼
    [Operator] → Target mRNA → Target protein
        ▲
        │ also represses (autoregulation)
        │
    Mce3R ──┘
```

Both Mce3R and Target are controlled by the SAME operator (they are divergently transcribed from the mce3R-yrbE3A intergenic region).

**Species:**
1. `mce3r_mRNA` — Mce3R mRNA count
2. `mce3r_protein` — Mce3R protein count (this IS the repressor)
3. `target_mRNA` — Downstream target mRNA count
4. `target_protein` — Downstream target protein count
5. `op_state` — Operator state (0, 1, 2, 3)

**Reactions (18 total):**

| # | Reaction | Propensity |
|---|----------|-----------|
| 0-3 | Operator: U→S, S→U, U→W, W→U | k_on × [Mce3R_free]_nM, k_off_s, k_on × [Mce3R_free]_nM, k_off_w |
| 4-5 | Operator: S→D, D→S | k_on × [Mce3R_free]_nM, k_off_w |
| 6-7 | Operator: W→D, D→W | k_on × [Mce3R_free]_nM, k_off_s |
| 8 | Mce3R transcription | k_txn[op_state] × alpha_mce3r |
| 9 | Target transcription | k_txn[op_state] × alpha_target |
| 10 | Mce3R translation | k_translation × mce3r_mRNA |
| 11 | Target translation | k_translation × target_mRNA |
| 12 | Mce3R mRNA decay | gamma_mRNA × mce3r_mRNA |
| 13 | Target mRNA decay | gamma_mRNA × target_mRNA |
| 14 | Mce3R protein decay | gamma_protein × mce3r_protein |
| 15 | Target protein decay | gamma_protein × target_protein |
| 16 | Mce3R mRNA burst (optional) | 0 (placeholder for future) |
| 17 | Target mRNA burst (optional) | 0 (placeholder for future) |

where:
- `[Mce3R_free]_nM = max(0, mce3r_protein - n_bound[op_state]) * nM_per_molecule`
- `alpha_mce3r` and `alpha_target` are relative transcription strengths for each direction (both default to 1.0, but they can differ since the genes are divergent)

```python
class TwoSpeciesModel:
    """
    Two-species Mce3R circuit model.

    Extends OperatorModel with:
    - Separate mRNA/protein species for repressor and target
    - Explicit autoregulation (Mce3R protein feeds back to binding)
    - Separate transcription rates for each direction (alpha_mce3r, alpha_target)
    """

    def __init__(
        self,
        Kd_strong: float = None,
        Kd_weak: float = None,
        k_on: float = None,
        block_strong: float = None,
        block_weak: float = None,
        k_max: float = None,
        alpha_mce3r: float = 1.0,
        alpha_target: float = 1.0,
    ):
        ...

    def get_numba_arrays(self) -> dict:
        """
        Returns all arrays needed by two_species_gillespie.simulate_two_species_cell().
        Same format as OperatorModel.get_numba_arrays() plus:
        - alpha_mce3r: float
        - alpha_target: float
        """
```

**Self-tests (3):**
1. Model arrays have correct dtypes for numba
2. When alpha_mce3r = alpha_target = 1.0, the two species have similar mean expression
3. Mce3R protein count feeds back into binding propensity (not fixed)

---

### 6.4 `two_species_gillespie.py`

**Purpose:** Numba-accelerated Gillespie engine for the two-species model.

**CRITICAL:** This must be a NEW `@numba.njit` function. Do NOT modify the existing `gillespie_engine.py`. The existing engine simulates one species; this simulates four species (2 mRNA + 2 protein) + operator state.

```python
@numba.njit
def simulate_two_species_cell(
    k_txn: np.ndarray,          # float64[4] transcription rates per state
    n_bound: np.ndarray,        # int64[4] number of Mce3R bound per state
    k_on: float,
    k_off_strong: float,
    k_off_weak: float,
    k_translation: float,
    gamma_mRNA: float,
    gamma_protein: float,
    nM_per_molecule: float,
    alpha_mce3r: float,
    alpha_target: float,
    t_max: float,
    t_burn_in: float,
    seed: int
) -> tuple:
    """
    Returns: (
        final_mce3r_protein: int,
        final_target_protein: int,
        final_mce3r_mRNA: int,
        final_target_mRNA: int,
        final_op_state: int
    )

    The key difference from simulate_cell:
    - State vector: (op_state, mce3r_mRNA, mce3r_protein, target_mRNA, target_protein)
    - FREE Mce3R for binding: mce3r_protein - n_bound[op_state]
    - Both transcription reactions use the SAME operator state
    - Mce3R protein count directly determines binding propensity (TRUE autoregulation)
    """
```

```python
def run_two_species_population(
    model: TwoSpeciesModel,
    n_cells: int,
    master_seed: int,
    environment: str = 'baseline'
) -> dict:
    """
    Run n_cells simulations of the two-species model.

    Returns dict:
    - mce3r_proteins: np.ndarray[n_cells]
    - target_proteins: np.ndarray[n_cells]
    - mce3r_mRNAs: np.ndarray[n_cells]
    - target_mRNAs: np.ndarray[n_cells]
    - op_states: np.ndarray[n_cells]
    """
```

**Self-tests (4):**
1. Numba compiles without error
2. All species counts >= 0
3. Operator state in {0, 1, 2, 3}
4. Mean mce3r_protein > 0 (autoregulation doesn't collapse to zero — negative feedback stabilizes)

---

### 6.5 `mutual_information.py`

**Purpose:** Estimate mutual information between environmental input and gene expression output.

```python
def estimate_mutual_information(
    distributions: dict,
    n_bins: int = 50,
    method: str = 'ksg'
) -> float:
    """
    Estimate I(Environment; Target_protein).

    Args:
        distributions: dict mapping environment_name → np.ndarray of target_protein counts
        n_bins: number of bins for histogram method
        method: 'histogram' or 'ksg' (Kraskov-Stogbauer-Grassberger, more accurate)

    For 'histogram' method:
    1. Pool all distributions with environment labels
    2. Compute joint distribution P(env, protein_bin)
    3. Compute marginals P(env) and P(protein_bin)
    4. MI = sum P(env, bin) * log2(P(env, bin) / (P(env) * P(bin)))

    For 'ksg' method:
    1. Use scipy.special.digamma-based KSG estimator
    2. More robust for continuous variables with limited samples

    Returns MI in bits.
    """

def mutual_information_vs_concentration(
    architecture: str,
    concentrations_nM: np.ndarray,
    environments: list[str],
    n_cells: int = 5000,
    seed: int = 42
) -> pd.DataFrame:
    """
    For each Mce3R concentration, simulate cells under each environment,
    compute MI between environment label and target protein count.

    Returns DataFrame: concentration, MI_bits, architecture
    """
```

**Self-tests (2):**
1. MI >= 0 (mutual information is non-negative)
2. MI for identical distributions = 0

---

### 6.6 `persistence_threshold.py`

**Purpose:** Define persistence threshold and quantify persister fractions.

```python
def define_threshold(
    proteins_unregulated: np.ndarray,
    method: str = 'percentile'
) -> float:
    """
    Define the persistence threshold.

    Methods:
    - 'percentile': 10th percentile of unregulated (Condition D) distribution.
      Cells below this are in a low-expression state consistent with slow growth.
    - 'gmm': Fit 2-component GMM to Condition A, use intersection point.
    - 'absolute': Fixed threshold at 100 proteins (biologically motivated).

    Default: 'percentile' — most robust, doesn't assume bimodality.
    """

def compute_persister_fractions(
    results: dict,
    threshold: float
) -> pd.DataFrame:
    """
    For each architecture × environment combination, compute:
    - n_below_threshold: cells with target_protein < threshold
    - persister_fraction: n_below / n_total
    - fold_enrichment: persister_fraction / persister_fraction_symmetric_baseline

    Returns DataFrame: architecture, environment, n_cells, persister_fraction,
                       fold_enrichment, mean_target, cv_target
    """

def persistence_phase_diagram(
    ratios: np.ndarray,
    environments: list[str],
    n_cells_per_point: int = 2000,
    seed: int = 42
) -> pd.DataFrame:
    """
    2D sweep: asymmetry ratio × environment → persister fraction.

    For each (ratio, environment) pair:
    1. Set up asymmetric operator with given ratio (same geometric mean Kd logic as Condition E)
    2. Apply environment multiplier
    3. Simulate n_cells
    4. Compute persister fraction

    Returns DataFrame: ratio, environment, persister_fraction, cv, mean_target
    """
```

**Self-tests (3):**
1. Persister fraction in [0, 1]
2. No-regulation condition has persister fraction consistent with threshold definition
3. Asymmetric operator has higher persister fraction than symmetric under same environment

---

### 6.7 `run_environmental_conditions.py`

**Purpose:** Run all architecture × environment combinations.

```python
ARCHITECTURES = ['asymmetric', 'symmetric', 'single_site']
ENVIRONMENTS = ['baseline', 'cholesterol', 'acidic_pH', 'host_like']

def run_all_environmental(
    n_cells: int = 10000,
    master_seed: int = 12345
) -> dict:
    """
    Run 3 architectures × 4 environments = 12 conditions.
    For each, run both single-species (existing engine) and two-species (new engine).

    Total simulations: 12 × 2 = 24 conditions.

    For single-species: use existing gillespie_engine.run_population() with modified k_on.
    For two-species: use two_species_gillespie.run_two_species_population().

    Save each condition's results as env_condition_{arch}_{env}_{model}.npz

    Returns dict of dicts keyed by (architecture, environment, model_type).
    """
```

**Cell counts:**
- Single-species: 10,000 cells per condition (fast, ~2 sec each)
- Two-species: 10,000 cells per condition (slower, ~10 sec each with numba)
- Phase diagram: 2,000 cells per point × 10 ratios × 4 environments = 80,000 cells

**Estimated wall time:** ~10 minutes total for Phase 6.

**Self-tests (4):**
1. All 24 conditions complete without error
2. All NPZ files written and loadable
3. CV(asymmetric) > CV(symmetric) under every environment
4. Host-like environment produces lower mean expression than baseline (derepression)

---

### 6.8 `environmental_main.py`

```python
def run_phase6() -> dict:
    """
    Orchestrate Phase 6:
    1. Run all environmental conditions (run_environmental_conditions.py)
    2. Compute mutual information (mutual_information.py)
    3. Define persistence threshold and compute fractions (persistence_threshold.py)
    4. Generate persistence phase diagram
    5. Write results to results/phase6/
    6. Write phase6_summary.json

    Returns standard summary dict.
    """
```

**Output files:**
- `results/phase6/env_condition_*.npz` — Per-condition simulation data
- `results/phase6/two_species_results.npz` — Two-species circuit data
- `results/phase6/mutual_information.csv` — MI values for each architecture × concentration
- `results/phase6/persistence_fractions.csv` — Persister fractions for all conditions
- `results/phase6/phase_diagram.csv` — Ratio × environment → persister fraction
- `results/phase6/phase6_summary.json`

---

## 7. PHASE 7: EXTENDED FIGURES

### 7.1 Figure List

| Figure | File | Content |
|--------|------|---------|
| Fig 9 | `fig9_repression_curves.py` | Repression fold vs [Mce3R] for 4 architectures. X-axis: log10([Mce3R] nM), Y-axis: fold repression. 4 colored curves. Highlight transition region. |
| Fig 10 | `fig10_cooperativity.py` | (A) Corner plot of MCMC posteriors for ω and ΔG_spacer. (B) Heatmap: repression strength vs (ω, ΔG_spacer). |
| Fig 11 | `fig11_environmental_distributions.py` | 3×4 panel grid: rows = architectures, columns = environments. Each panel: target protein histogram/density. Color: same as existing fig_style.py. |
| Fig 12 | `fig12_mutual_information.py` | MI (bits) vs [Mce3R] for asymmetric vs symmetric vs single-site. Show that asymmetric transmits more information. |
| Fig 13 | `fig13_persistence_phase.py` | 2D heatmap: X = asymmetry ratio, Y = environment (categorical), color = persister fraction. Annotate the native operator position. |
| Fig 14 | `fig14_two_species_traces.py` | (A) Single-cell time traces of Mce3R protein AND target protein for 5 cells under host-like conditions. (B) Scatter plot: Mce3R vs target protein (10,000 cells), colored by operator state. |

### 7.2 Style

Reuse `phase4_figures/figure_style.py` for colors and style. Import it as:
```python
from phase4_figures.figure_style import COLORS, STYLE, CONDITION_LABELS, apply_style
```

Add new colors for environments:
```python
ENV_COLORS = {
    'baseline': '#6B7280',      # Gray
    'cholesterol': '#F59E0B',   # Amber
    'acidic_pH': '#EF4444',     # Red
    'host_like': '#7C3AED',     # Violet
}
```

All figures: 300 dpi PNG, saved to `results/extended_figures/`.

**Self-tests per figure (1 each):**
1. File exists, size > 10 KB, is valid PNG

---

## 8. UPDATED ORCHESTRATOR

Modify `main.py` to support the new phases.

**IMPORTANT:** Do NOT change the existing Phase 1-4 logic. ADD new logic after Phase 4.

Add to `main.py`:

```python
# After Phase 4 completes:
if args.phase is None or args.phase == 5:
    # Phase 5: Thermodynamic model
    p5 = subprocess.Popen([sys.executable, '-m', 'phase5_thermodynamic.thermodynamic_main'], ...)
    p5.wait()

if args.phase is None or args.phase == 6:
    # Phase 6: Environmental extensions (depends on Phase 5)
    p6 = subprocess.Popen([sys.executable, '-m', 'phase6_environmental.environmental_main'], ...)
    p6.wait()

if args.phase is None or args.phase == 7:
    # Phase 7: Extended figures (depends on Phase 6)
    p7 = subprocess.Popen([sys.executable, '-m', 'phase7_extended_figures.extended_figures_main'], ...)
    p7.wait()
```

**Updated dependency graph:**
```
Phase 1 (bioinformatics) ──┐
                           ├── Phase 3 (analysis) ── Phase 4 (figures)
Phase 2 (simulation) ──────┘
                                                          │
                                                          ▼
                                                     Phase 5 (thermodynamic)
                                                          │
                                                          ▼
                                                     Phase 6 (environmental)
                                                          │
                                                          ▼
                                                     Phase 7 (extended figures)
```

Phases 5, 6, 7 are sequential. Phase 5 needs Phase 1 outputs (MEME PWM). Phase 6 needs Phase 5 outputs (omega). Phase 7 needs Phase 6 outputs.

---

## 9. PARAMETER ADDITIONS

Add these to `config/parameters.py` in a new section at the bottom. Do NOT modify any existing parameters.

```python
# ============================================================
# PHASE 5-7 PARAMETERS (Thermodynamic + Environmental Extensions)
# ============================================================

# Thermodynamic model
TEMPERATURE_K = 310.0                    # 37°C in Kelvin
kT = 0.001987 * TEMPERATURE_K            # kcal/mol (R = 1.987e-3 kcal/(mol·K))
BACKGROUND_FREQ = {                       # H37Rv nucleotide frequencies
    'A': (1 - GC_CONTENT) / 2,           # 0.172
    'C': GC_CONTENT / 2,                 # 0.328
    'G': GC_CONTENT / 2,                 # 0.328
    'T': (1 - GC_CONTENT) / 2,           # 0.172
}

# MCMC settings
MCMC_N_WALKERS = 32
MCMC_N_STEPS = 5000
MCMC_N_BURN = 1000
MCMC_SEED = 42

# Operator site positions within 123-bp operator
STRONG_SITE_START = 98                    # 0-indexed, last 25 bp
STRONG_SITE_END = 123
WEAK_SITE_START = 0
WEAK_SITE_END = 25

# Environmental model
ENV_CHOLESTEROL_MULTIPLIER = 0.3          # fold reduction in effective Mce3R binding
ENV_ACID_PH_MULTIPLIER = 0.7
ENV_HOST_LIKE_MULTIPLIER = 0.2
ENV_ACID_NOISE_SCALE = 1.2               # fold increase in protein degradation
ENV_HOST_NOISE_SCALE = 1.3

# Two-species model
ALPHA_MCE3R = 1.0                        # relative Mce3R transcription strength
ALPHA_TARGET = 1.0                       # relative target transcription strength

# Phase 6 simulation settings
N_CELLS_ENVIRONMENTAL = 10000            # per architecture × environment
N_CELLS_PHASE_DIAGRAM = 2000             # per ratio × environment point
N_CELLS_MI = 5000                        # per MI calculation
PHASE6_MASTER_SEED = 54321

# Persistence threshold
PERSISTENCE_THRESHOLD_METHOD = 'percentile'  # 'percentile', 'gmm', or 'absolute'
PERSISTENCE_PERCENTILE = 10                   # 10th percentile of unregulated

# Mutual information
MI_N_BINS = 50
MI_METHOD = 'histogram'                  # 'histogram' or 'ksg'
MI_CONCENTRATIONS = np.logspace(-1, 4, 30)  # 0.1 to 10,000 nM, 30 points

# Figure settings (extend existing)
ENV_COLORS = {
    'baseline': '#6B7280',
    'cholesterol': '#F59E0B',
    'acidic_pH': '#EF4444',
    'host_like': '#7C3AED',
}
```

---

## 10. SELF-TESTING REQUIREMENTS

Every module must include a `if __name__ == '__main__':` block that runs self-tests.

Format:
```python
if __name__ == '__main__':
    print("[Phase 5] Running self-tests for module_name...")
    n_pass, n_fail = 0, 0

    # Test 1: description
    try:
        assert condition, "explanation"
        print("  [PASS] Test 1: description")
        n_pass += 1
    except Exception as e:
        print(f"  [FAIL] Test 1: description — {e}")
        n_fail += 1

    print(f"\n[Phase 5] module_name: {n_pass} passed, {n_fail} failed")
    sys.exit(1 if n_fail > 0 else 0)
```

### Test Counts

| Module | Expected tests |
|--------|---------------|
| energy_calibration.py | 4 |
| partition_function.py | 5 |
| cooperativity_inference.py | 4 |
| operator_classification.py | 3 |
| dna_shape.py | 1 |
| environmental_signals.py | 2 |
| two_species_model.py | 3 |
| two_species_gillespie.py | 4 |
| mutual_information.py | 2 |
| persistence_threshold.py | 3 |
| run_environmental_conditions.py | 4 |
| Each figure file | 1 |
| **Total** | ~41 |

---

## 11. CRITICAL DESIGN DECISIONS

### Decision 1: Cooperativity is fitted, not assumed

Do NOT hardcode ω = 1. The whole point of Aim 2 is to infer cooperativity from the data. Start with a wide prior and let MCMC determine the posterior. If the posterior is centered near ω = 1, that's a result (no cooperativity). If ω >> 1, the sites are cooperative.

### Decision 2: Environmental multipliers are phenomenological

The exact values of `ENV_CHOLESTEROL_MULTIPLIER` etc. are not precisely known. The key prediction is QUALITATIVE: asymmetric operators amplify noise MORE under stress. Run sensitivity analysis on the multiplier values to show robustness.

### Decision 3: Two-species model uses same operator for both genes

Mce3R and yrbE3A are divergently transcribed from the same intergenic operator. Both genes see the SAME operator state. This means operator state switching simultaneously affects both the repressor and the target — creating correlated noise.

### Decision 4: DNAshapeR is optional

Do not let the pipeline fail because R is not installed. Check for R at the start of Phase 5 and skip `dna_shape.py` if unavailable. All other Phase 5 modules work without it.

### Decision 5: Reuse existing Gillespie for single-species comparisons

For the environmental × architecture grid with single-species models, reuse the existing `gillespie_engine.simulate_cell()` with modified k_on (via environmental multipliers). This ensures consistency with Phases 2-4.

### Decision 6: Mutual information uses histogram method by default

KSG is more accurate but slower and harder to implement correctly. Start with the histogram method (simpler, more transparent). If the estimates are noisy, switch to KSG.

---

## 12. SCIENTIFIC PREDICTIONS TO VALIDATE

These are the predictions from `PROJECT_PROPOSAL.md` that Phases 5-7 must test:

| # | Prediction | How to test | Expected result |
|---|-----------|-------------|-----------------|
| 1 | Weak site creates intermediate occupancy regime | Phase 5: Plot P(states) vs [Mce3R]. At intermediate concentrations, the asymmetric operator should have more time in partially-bound states than symmetric. | P(S) + P(W) > 0.3 at K_half for asymmetric; < 0.1 for symmetric |
| 2 | Native operator has steeper transition | Phase 5: Compare effective Hill coefficients. | n_H(asymmetric) > n_H(symmetric) if ω > 1 |
| 3 | Asymmetry amplifies noise MORE under stress | Phase 6: Compare CV(asymmetric)/CV(symmetric) ratio across environments. | Ratio increases from baseline to host_like |
| 4 | Two-species model shows correlated noise | Phase 6: Compute Pearson correlation between Mce3R and target protein across cells. | Negative correlation (when Mce3R high → target low, and vice versa) |
| 5 | Asymmetric operator transmits more information | Phase 6: MI(asymmetric) > MI(symmetric) at physiologically relevant [Mce3R]. | MI difference > 0.1 bits |
| 6 | Persistence fraction increases with asymmetry under stress | Phase 6: Phase diagram shows persister_fraction increases with ratio, especially under host_like. | Monotonic increase, steeper under stress |
| 7 | Operators classify into distinct functional categories | Phase 5: At least 2 of 3 classification categories populated among top-20 operators. | ≥2 categories |

**If a prediction FAILS:**
- Report it honestly as an unexpected result.
- Do NOT tune parameters to force the prediction to work.
- Investigate WHY — a failed prediction is scientifically interesting.
- Flag it in the summary with `n_science_unexpected += 1`.

---

## 13. BUILD ORDER

Write files in this exact order:

1. Add new parameters to `config/parameters.py` (append, don't modify existing)
2. `phase5_thermodynamic/__init__.py`
3. `phase5_thermodynamic/energy_calibration.py`
4. `phase5_thermodynamic/partition_function.py`
5. `phase5_thermodynamic/cooperativity_inference.py`
6. `phase5_thermodynamic/operator_classification.py`
7. `phase5_thermodynamic/dna_shape.py`
8. `phase5_thermodynamic/thermodynamic_main.py`
9. `phase6_environmental/__init__.py`
10. `phase6_environmental/environmental_signals.py`
11. `phase6_environmental/two_species_model.py`
12. `phase6_environmental/two_species_gillespie.py`
13. `phase6_environmental/mutual_information.py`
14. `phase6_environmental/persistence_threshold.py`
15. `phase6_environmental/run_environmental_conditions.py`
16. `phase6_environmental/environmental_main.py`
17. `phase7_extended_figures/__init__.py`
18. `phase7_extended_figures/fig9_repression_curves.py`
19. `phase7_extended_figures/fig10_cooperativity.py`
20. `phase7_extended_figures/fig11_environmental_distributions.py`
21. `phase7_extended_figures/fig12_mutual_information.py`
22. `phase7_extended_figures/fig13_persistence_phase.py`
23. `phase7_extended_figures/fig14_two_species_traces.py`
24. `phase7_extended_figures/extended_figures_main.py`
25. Update `main.py` (append Phase 5-7 logic)
26. Run full pipeline: `python main.py`

Write ONE file at a time. Run self-tests after each file. Do NOT proceed to the next file until the current one passes.

---

## 14. RUNNING THE EXTENDED PIPELINE

```bash
conda activate mce3r
cd "/Users/aayanalwani/tb project/mce3r_stochastic"

# Install new dependencies
pip install emcee corner arviz

# Run everything (Phases 1-7)
python main.py

# Or run only new phases
python main.py --phase 5
python main.py --phase 6
python main.py --phase 7
```

**Expected total wall time:** ~15-20 minutes (existing 11 min + new ~5-10 min)

---

## 15. SUMMARY

This architecture extends the existing COMPLETE Mce3R pipeline with three new phases:

| Phase | What it does | Key output |
|-------|-------------|------------|
| 5 | Thermodynamic model: PWM→energy, partition function, MCMC cooperativity, operator classification | repression_curves.csv, mcmc_posteriors.npz, operator_classifications.csv |
| 6 | Environmental extensions: cholesterol/pH signals, two-species circuit, mutual information, persistence thresholds | persistence_fractions.csv, mutual_information.csv, phase_diagram.csv |
| 7 | 6 new publication figures | fig9-14_*.png |

The central question shifts from "does asymmetry increase noise?" (answered: yes, 20%) to "WHY does the architecture produce noise, and does that noise create persisters under host conditions?"
