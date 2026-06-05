# MCE3R Stochastic Gene Expression: Architecture Document

**Project:** ENIGMA — Expression Noise In Gene-regulatory Mechanisms and Architecture
**System:** Mce3R repressor, *Mycobacterium tuberculosis*
**Version:** 4.0 (March 2026)

---

## Table of Contents

1. [What This Project Actually Does](#1-what-this-project-actually-does)
2. [Aims (Honest Version)](#2-aims-honest-version)
3. [What's Missing and Why](#3-whats-missing-and-why)
4. [Architecture Overview](#4-architecture-overview)
5. [Configuration](#5-configuration)
6. [Phase 1: Motif Scanning for Mce3R](#6-phase-1-motif-scanning-for-mce3r)
7. [Phase 2: Gillespie Stochastic Simulation](#7-phase-2-gillespie-stochastic-simulation)
8. [Phase 3: Statistical Analysis](#8-phase-3-statistical-analysis)
9. [Phase 5: Thermodynamic Model & MCMC](#9-phase-5-thermodynamic-model--mcmc)
10. [Phase 6: Environmental Perturbations & Two-Species](#10-phase-6-environmental-perturbations--two-species)
11. [Phase 8: Temporal Dynamics](#11-phase-8-temporal-dynamics)
12. [Phase 9: Noise Quenching](#12-phase-9-noise-quenching)
13. [Figures](#13-figures)
14. [Data Flow & Dependencies](#14-data-flow--dependencies)
15. [Known Limitations](#15-known-limitations)
16. [Parameters](#16-parameters)

---

## 1. What This Project Actually Does

One thing: quantifies how the **20.4-fold binding-site affinity asymmetry** of the Mce3R operator affects gene expression noise in *M. tuberculosis*.

The operator has two binding sites (Kd = 2.4 nM and 49.0 nM, Panagoda et al. 2024). We compare this wild-type asymmetric architecture against a symmetric control (geometric mean Kd at both sites) and a single-site ablation. The comparison is done via:

- **Gillespie SSA** (50,000 cells per condition, Numba-accelerated)
- **4-state operator model** (unbound, strong-only, weak-only, both-bound)
- **Thermodynamic partition function** with MCMC-inferred cooperativity
- **Environmental stress modeling** (4 conditions x 3 architectures)
- **Temporal analysis** (autocorrelation, dwell times)
- **Noise quenching dose-response** (symmetrization sweep with IC50)

The central result: asymmetry generates ~19% more noise (CV 0.187 vs 0.157), the autocorrelation time is 12% longer (18.7 hr vs 16.7 hr), and the IC50 of noise quenching is 30% symmetrization.

This is a **single-TF computational study**, not a genome-wide analysis.

---

## 2. Aims (Honest Version)

| Aim | What It Actually Does | Novel? |
|-----|----------------------|--------|
| **Aim 1** | Scans H37Rv for Mce3R binding sites using MEME/FIMO on 3 yrbE3A orthologs. This is **motif scanning for one TF**, not genome-wide operator classification. | No. Standard bioinformatics. |
| **Aim 2** | Gillespie simulation comparing 4 operator architectures (asymmetric, symmetric, single-site, unregulated). Thermodynamic partition function with MCMC cooperativity inference. | **Partially.** First quantitative noise characterization of a structurally resolved asymmetric bacterial operator with experimentally measured Kd values. |
| **Aim 3a** | Statistical characterization: CV, Fano, bootstrap CIs, sensitivity analysis, GMM model selection. Environmental perturbations (cholesterol, acidic pH, host-like). Two-species autoregulatory circuit. | No. Standard stochastic gene expression analysis. The environmental and two-species extensions add scope but not novelty. |
| **Aim 3b** | Temporal dynamics: autocorrelation time, dwell times, PSD. | **Partially.** Connecting operator architecture to noise temporal memory is not well-explored for asymmetric operators. |
| **Aim 3c** | Noise quenching: symmetrization sweep, IC50 of noise. | **Yes, if framed correctly.** First dose-response curve for noise quenching in a TB system. |

**What's genuinely novel:** Items 2, 3b, and 3c. Everything else is scaffolding.

---

## 3. What's Missing and Why

### 3.1 Evolutionary Simulation (Former Aim 4)

The earlier proposal included a genetic algorithm simulating persistence selection pressure on operator architecture. This was the strongest novelty angle: showing that asymmetry could be evolutionarily selected for. It was dropped in favor of more analysis phases and figures. **This was a bad trade.** The evolutionary simulation would have answered "why asymmetry?" rather than just "what does asymmetry do?"

**Status:** Not implemented. Should be the next priority if this project continues.

### 3.2 Genome-Wide Operator Classification

The original pitch was to classify noncanonical operator architectures across the entire M. tuberculosis genome. Phase 1 only scans for one TF's motif. A true genome-wide analysis would require: (1) scanning all ~200 TetR-family regulators in Mtb, (2) classifying operators as palindromic vs asymmetric vs tandem, (3) correlating architecture with expression noise. None of this is done.

**Status:** Not implemented. Phase 1 title has been corrected to "Motif Scanning for Mce3R."

### 3.3 Experimental Validation

All results are computational predictions. No wet-lab validation exists. The fold-change comparison to Santangelo 2009 is the closest thing to experimental grounding, and it's indirect.

---

## 4. Architecture Overview

### 4.1 Directory Structure

```
mce3r_stochastic/
|-- main.py                      # Orchestrator (Phases 1-4, parallel)
|-- config/parameters.py         # All parameters, single source of truth
|
|-- phase1_pipeline/             # Motif scanning (MEME/FIMO)
|-- phase2_simulation/           # Gillespie SSA engine
|-- phase3_analysis/             # Statistical analysis
|-- phase5_thermodynamic/        # Partition function, MCMC
|-- phase6_environmental/        # Environmental + two-species
|-- phase8_temporal/             # Autocorrelation, dwell times, PSD
|-- phase9_quenching/            # Symmetrization sweep, IC50
|-- figures/                     # All 19 publication figures (unified)
|
|-- data/                        # Genomes, sequences, external data
|-- results/                     # All outputs organized by phase
+-- logs/                        # Execution logs
```

### 4.2 File Counts

| Module | Files | Notes |
|--------|:---:|-------|
| config | 2 | Parameters + init |
| phase1_pipeline | 8 | Could be 3-4 |
| phase2_simulation | 6 | Reasonable |
| phase3_analysis | 8 | Could be 4-5 |
| phase5_thermodynamic | 7 | Could be 4 |
| phase6_environmental | 13 | persistence_threshold has 3 versions that should be 1 |
| phase8_temporal | 4 | Fine |
| phase9_quenching | 3 | Fine |
| figures | 22 | 19 figure scripts + style + orchestrator + init |
| **Total** | **~73** | |

### 4.3 Phases

7 computation phases + 1 unified figure module:

| Phase | Module | What It Computes |
|-------|--------|-----------------|
| 1 | phase1_pipeline | Mce3R motif scanning (MEME/FIMO) |
| 2 | phase2_simulation | Gillespie SSA, 4 conditions + sweep |
| 3 | phase3_analysis | Noise metrics, bootstrap, sensitivity, validation |
| 5 | phase5_thermodynamic | Partition function, MCMC cooperativity |
| 6 | phase6_environmental | 24 environmental conditions, two-species, MI, persistence |
| 8 | phase8_temporal | Autocorrelation, dwell times, PSD |
| 9 | phase9_quenching | Symmetrization sweep, IC50 |
| -- | figures | All 19 publication figures (reads from computation phases) |

Figures are not a "phase." They're visualization of existing results.

---

## 5. Configuration

All parameters in `config/parameters.py`. Key values:

### Binding Kinetics (Panagoda et al. 2024)
| Parameter | Value | Units |
|-----------|-------|-------|
| Kd_strong | 2.4 | nM |
| Kd_weak | 49.0 | nM |
| k_on | 0.0167 | nM^-1 min^-1 |
| block_strong | 0.85 | fraction |
| block_weak | 0.50 | fraction |

### Gene Expression
| Parameter | Value | Units | Source |
|-----------|-------|-------|--------|
| k_max | 0.15 | mRNA/min | Mtb transcriptomics |
| k_translation | 0.5 | protein/mRNA/min | Taniguchi 2010 |
| mRNA t_half | 9.5 | min | Rustad 2013 |
| protein t_half | 1500 | min (~25 hr) | Dilution-dominated |

### Controls
| Control | Design | Rationale |
|---------|--------|-----------|
| Symmetric (B) | Kd = sqrt(2.4 x 49) = 10.84 nM both sites | Geometric mean preserves total regulatory capacity |
| Single-site (C) | Kd_weak = 1e12 (disabled) | Tests weak-site contribution |
| Unregulated (D) | k_on = 0 | Constitutive baseline |

### Simulation
- 50,000 cells per condition (A/B/C/D)
- 5,000 cells per sweep point (Condition E, 10 ratios)
- t_max = 30,000 min, burn-in = 15,000 min
- 20 trace cells per condition (dt = 10 min)
- Master seed: 42

### MCMC
- 32 walkers, 5,000 steps, 1,000 burn-in (emcee)
- Inferred: ln(omega), dG_spacer

### Transcription Rate Array (4 states)
```
State 0 (empty):       0.1500 mRNA/min (100%)
State 1 (strong):      0.0225 mRNA/min (15%)
State 2 (weak):        0.0750 mRNA/min (50%)
State 3 (both):        0.01125 mRNA/min (7.5%)
```

Self-validated at import time. Exit(1) on any inconsistency.

---

## 6. Phase 1: Motif Scanning for Mce3R

**What it is:** MEME/FIMO pipeline to find Mce3R binding sites in H37Rv.
**What it is NOT:** Genome-wide operator classification.

### Pipeline

```
download_genomes -> extract_upstream (200bp x ~4000 CDS) -> MEME (3 orthologs)
-> FIMO (genome scan, p < 1e-4) -> conservation_check (M. bovis, M. marinum)
```

### Key Design Decisions
- MEME input: 3 yrbE3A orthologs only (H37Rv, M. bovis, M. marinum)
- No `-pal` flag (allows asymmetric motifs)
- Mock fallback if MEME Suite not installed
- Conservation threshold: 80% identity over full motif

### Outputs
- `predicted_sites.csv` — ~1,400 FIMO hits
- `conservation_status.csv` — top 50 sites with cross-species identity
- `meme.txt` — PWM (feeds Phase 5 energy calibration)

---

## 7. Phase 2: Gillespie Stochastic Simulation

The computational core. 12-reaction SSA with Numba JIT compilation.

### Operator Model (`operator_model.py`)

Class `OperatorModel` pre-computes state-dependent transcription rates and binding parameters. Condition-specific overrides (symmetric Kd, disabled weak site, etc.) are passed as constructor arguments.

### Gillespie Engine (`gillespie_engine.py`)

**12 Reactions:**
- 0-7: Operator binding/unbinding (4-state transitions, concentration-dependent)
- 8: Transcription (state-dependent rate)
- 9: Translation
- 10-11: mRNA/protein decay

Free repressor concentration computed dynamically: R_nM = (protein - n_bound[state]) x nM_per_molecule

### Conditions

| Condition | Description | Key Parameter Override |
|-----------|-------------|----------------------|
| A | Wild-type asymmetric | Kd = 2.4 / 49.0 nM |
| B | Symmetric control | Kd = 10.84 / 10.84 nM |
| C | Single-site | Kd_weak = 1e12 |
| D | Unregulated | k_on = 0 |
| E | Sweep | 10 Kd ratios, constant geometric mean |

### Asymmetry Sweep (Condition E)
Holds geometric mean Kd constant while varying ratio. Isolates asymmetry from total regulatory capacity. Block fractions co-vary via power-law interpolation in log-ratio space.

### Outputs
- `condition_A/B/C/D.npz` — 50,000 cells each (proteins, mRNAs, op_states, traces)
- `condition_E_sweep.npz` — 10 ratios x 5,000 cells
- `phase2_summary.json`

---

## 8. Phase 3: Statistical Analysis

Six sub-analyses, all reading Phase 2 outputs:

| Module | What It Does |
|--------|-------------|
| `noise_metrics.py` | CV, Fano, bimodality coefficient, GMM fit (k=1,2,3) per condition |
| `bootstrap_ci.py` | 1,000-sample percentile bootstrap for CV, Fano, persister fraction |
| `statistical_tests.py` | KS test, Cohen's d, GMM likelihood-ratio test (pairwise) |
| `sensitivity_analysis.py` | 7 parameters varied +/-50%, 10 steps, 1,000 cells/step |
| `experimental_comparison.py` | Fold-change vs Santangelo 2009 |
| `negative_controls.py` | Shuffle test, symmetric control, Poisson baseline |

### Key Outputs
- `noise_metrics.csv` — 4 rows: mean, CV, Fano, GMM_k per condition
- `statistical_tests.csv` — 13 pairwise tests
- `sensitivity_data.csv` — 70 rows (7 params x 10 steps)

---

## 9. Phase 5: Thermodynamic Model & MCMC

### Energy Calibration (`energy_calibration.py`)

MEME PWM -> position-specific energy matrix -> calibrate to experimental Kd values using **two independent offsets** (one per binding site). Mean offset for novel sequences.

### Partition Function (`partition_function.py`)

4-state statistical mechanics:

```
Z = 1 + exp(-B(dG_s - mu)) + exp(-B(dG_w - mu)) + omega * exp(-B(dG_s + dG_w - 2mu + dG_spacer))
```

Computes state probabilities, mean transcription rate, repression fold, and dose-response curves for 4 architectures.

### MCMC Cooperativity (`cooperativity_inference.py`)

Fits omega (cooperativity) and dG_spacer to observed mean protein levels in Conditions A and B using emcee. Gaussian priors, R-hat convergence check.

### Operator Classification (`operator_classification.py`)

Classifies top 20 FIMO sites by Hill coefficient: digital_switch (n_H > 2), graded_repressor (n_H ~ 1), noise_modulator (significant weak-site contribution).

### Key Outputs
- `mcmc_posteriors.npz` — full chain (32 x 5000 x 2)
- `mcmc_summary.json` — omega median + 68% CI
- `repression_curves.csv` — 4 architectures x 200 concentrations
- `operator_classifications.csv` — top 20 sites

---

## 10. Phase 6: Environmental Perturbations & Two-Species

### Environments

| Environment | k_on multiplier | gamma_protein multiplier | Biological basis |
|-------------|:---:|:---:|------------------|
| baseline | 1.0 | 1.0 | In vitro |
| cholesterol | 0.3 | 1.0 | Reduced Mce3R binding |
| acidic_pH | 0.7 | 1.2 | Phagosomal stress |
| host_like | 0.2 | 1.3 | Combined macrophage stress |

### Condition Matrix
3 architectures x 4 environments x 2 models (single-species, two-species) = 24 simulations, 10,000 cells each.

### Two-Species Model
Divergent transcription: Mce3R and Target share one operator. Mce3R autoregulates (negative feedback). 18-reaction Numba Gillespie engine with dynamic free-repressor concentration.

### Additional Analyses
- **Mutual information:** I(Environment; Target_protein) via histogram entropy
- **Persistence threshold:** 3 methods (v2): tail percentile, absolute calibrated, fold-change
- **Posterior propagation:** 50 MCMC samples through Gillespie — 100% confirm CV(asym) > CV(sym)
- **TetR benchmark:** Biological control using real symmetric operator (Kd = 2.0 nM)

### Key Outputs
- 24 NPZ files (one per condition)
- `mutual_information.csv`, `persistence_fractions_v2.csv`, `phase_diagram.csv`
- `posterior_propagation.csv`, `tetR_benchmark.csv`

---

## 11. Phase 8: Temporal Dynamics

### Autocorrelation
FFT-based normalized ACF on 20 traces per condition. Extracts tau_c (1/e crossing) and integrated autocorrelation time.

**Result:** tau_c(A) = 18.7 hr, tau_c(B) = 16.7 hr. Asymmetric has 12% longer noise memory.

### Dwell Times
Theoretical mean dwell times from exit rates. Asymmetric operator has distinct State 1 vs State 2 dwells (0.181 vs 0.159 min); symmetric does not (0.121 = 0.121).

### Power Spectral Density
FFT-based PSD. All conditions dominated by low-frequency power (>96% below 1/500 cycles/min), consistent with slow protein dynamics (t_half = 25 hr).

### Key Outputs
- `autocorrelation.csv`, `acf_data.npz`
- `dwell_times.csv`
- `power_spectrum.csv`, `psd_data.npz`

---

## 12. Phase 9: Noise Quenching

### Symmetrization Sweep
20 dose points from fully asymmetric (frac=0) to fully symmetric (frac=1). Log-space Kd interpolation, linear block interpolation. 5,000 cells per dose.

### IC50 of Noise
- delta_cv = CV(asym) - CV(sym) = 0.0273
- IC50 = frac where delta_cv drops to 50% = **0.301** (30.1% symmetrization)
- Persister fraction drops from 55.7% to 1.5% across sweep (**36.6x reduction**)

### Key Outputs
- `symmetrization_sweep.csv` — 20 rows with CV, persister fraction, Kd ratios
- `phase9_summary.json`

---

## 13. Figures

19 figures across 3 generation modules. **This is too many for one paper.** A realistic publication would use 4-6 main figures + supplementary. The current set is comprehensive but should be triaged for any manuscript.

### Main Figures (Phase 4)
| Fig | Content | Would keep for paper? |
|-----|---------|:---:|
| 1 | Circular genome + FIMO sites | Yes (methods) |
| 2 | Protein distributions + GMM | **Yes (central result)** |
| 3 | CV vs asymmetry ratio sweep | **Yes (central result)** |
| 4 | Sensitivity tornado + heatmap | Supplementary |
| 5 | Experimental validation | Yes (validation) |
| 6 | Single-cell traces | Supplementary |
| 7 | BIC model comparison | Supplementary |
| 8 | Pipeline schematic | Maybe (methods) |

### Extended Figures (Phase 7)
| Fig | Content | Would keep? |
|-----|---------|:---:|
| 9 | Repression curves | Supplementary |
| 10 | MCMC posterior | Supplementary |
| 11 | Environmental distributions (3x4 grid) | Yes if environmental story told |
| 12 | Mutual information | Supplementary |
| 13 | Persistence phase diagram | **Yes (if persistence is a main result)** |
| 14 | Two-species traces | Supplementary |

### New Figures (Phase 10)
| Fig | Content | Would keep? |
|-----|---------|:---:|
| 15 | Autocorrelation + tau_c | **Yes (temporal novelty)** |
| 16 | Dwell times + occupancy | Supplementary |
| 17 | Power spectral density | Supplementary |
| 18 | Noise quenching dose-response | **Yes (quenching novelty)** |
| 19 | Summary dashboard | No (presentation only) |

**Recommended paper figures:** 2, 3, 5, 13, 15, 18 = 6 figures. Everything else is supplementary or cut.

---

## 14. Data Flow & Dependencies

```
Phase 1 (motif scan) ---------> Phase 5 (thermo, needs PWM)
                                    |
Phase 2 (Gillespie) --+            v
                       +--> Phase 3 (stats) --> Phase 6 (environmental, needs MCMC + Cond D)
                       |                            |
                       +--> Phase 8 (temporal,      +--> Phase 7 (figs 9-14)
                       |    needs traces)
                       |
                       +--> Phase 9 (quenching,
                            needs Gillespie engine)

Figure phases (4, 7, 10) read from their respective computation phases.
```

**Execution:**
1. Phase 1 + Phase 2 run in **parallel**
2. Phase 3 requires Phase 2
3. Phase 5 requires Phase 1 + Phase 2
4. Phase 6 requires Phase 5
5. Phase 8 requires Phase 2 traces
6. Phase 9 requires Phase 2 engine
7. Figures at any point after their data phases

---

## 15. Known Limitations

1. **All computational, no wet-lab validation.** Predictions only.
2. **Single-TF study.** Phase 1 is motif scanning, not genome-wide classification.
3. **Evolutionary simulation not implemented.** The strongest novelty angle (Aim 4: genetic algorithm under persistence selection) was dropped.
4. **Wide MCMC posterior.** omega = 1.08 [0.20-5.89]. Limited calibration data (2 conditions).
5. **91 Python files is overengineered.** Could be ~40-50 with cleanup.
6. **19 figures is too many.** 6 for a paper, rest supplementary.
7. **Environmental conditions are scalar multipliers**, not mechanistic models.
8. **tau_c difference (12%) may be dominated by protein half-life** rather than operator architecture.
9. **CV difference (0.031) is moderate.** Other noise sources (TA modules, metabolic fluctuations) are likely larger.
10. **Persister threshold is model-dependent.** Different methods give different fractions.

---

## 16. Parameters

Full parameter reference in `config/parameters.py` (72+ parameters with literature citations). Key categories:

- **Binding kinetics** (5 params): Kd_strong, Kd_weak, k_on, k_off_s, k_off_w
- **Transcription** (3): k_max, block_strong, block_weak
- **Translation & decay** (4): k_translation, gamma_mRNA, gamma_protein, burst_size
- **Simulation** (6): n_cells, t_max, t_burn_in, seeds, traces
- **Controls** (4): Kd_symmetric, block_symmetric, Kd_disabled, k_on_disabled
- **MCMC** (4): walkers, steps, burn-in, seed
- **Thermodynamic** (4): temperature, kT, background freq
- **Environmental** (4 environments x 2 multipliers)
- **MEME/FIMO** (11): modes, widths, thresholds
- **Operator** (5): sequence, length, coordinates

---

*Document version: 4.0, March 30, 2026*
*This is a single-TF computational study, not a genome-wide analysis.*
*Evolutionary simulation (strongest novelty angle) is not yet implemented.*
