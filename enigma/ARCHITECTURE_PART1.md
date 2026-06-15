# ENIGMA Full Project Architecture v2

**Title:** Temporal Noise Dynamics and Therapeutic Noise Quenching in the Asymmetric Mce3R Operator of *Mycobacterium tuberculosis*

**Acronym:** ENIGMA (Expression Noise In Gene-regulatory Mechanisms and Architecture)

**Author:** Aayan Alwani

**Document scope:** This is the FIRST HALF of the definitive technical specification. Every algorithm, every equation, every parameter, every file, every function signature. Written for a developer who needs to reimplement the entire project from scratch.

---

## Table of Contents

- [1. Project Overview](#1-project-overview)
  - [1.1 Description](#11-description)
  - [1.2 Central Hypothesis](#12-central-hypothesis)
  - [1.3 Aims and Stretch Goals](#13-aims-and-stretch-goals)
  - [1.4 The 10 Phases and Dependency Graph](#14-the-10-phases-and-dependency-graph)
  - [1.5 What Is Novel vs. Established](#15-what-is-novel-vs-established)
- [2. Configuration -- Every Parameter](#2-configuration----every-parameter)
  - [2.1 Primary Parameters Table](#21-primary-parameters-table)
  - [2.2 Derived Parameters](#22-derived-parameters)
  - [2.3 k_txn Array Computation](#23-k_txn-array-computation)
  - [2.4 Phase 5-7 Parameters](#24-phase-5-7-parameters)
- [3. Phase 1: Motif Discovery and Genome-Wide Scanning](#3-phase-1-motif-discovery-and-genome-wide-scanning)
  - [3.1 Genome Acquisition](#31-genome-acquisition-download_genomespy)
  - [3.2 Upstream Sequence Extraction](#32-upstream-sequence-extraction-extract_upstreampy)
  - [3.3 Background Model](#33-background-model)
  - [3.4 MEME Motif Discovery](#34-meme-motif-discovery-run_memepy)
  - [3.5 FIMO Genome-Wide Scanning](#35-fimo-genome-wide-scanning-run_fimopy)
  - [3.6 Conservation Analysis](#36-conservation-analysis-conservation_checkpy)
  - [3.7 Pipeline Orchestrator](#37-pipeline-orchestrator-pipeline_mainpy)
  - [3.8 Phase 1 Self-Tests](#38-phase-1-self-tests)
- [4. Phase 2: Gillespie Stochastic Simulation](#4-phase-2-gillespie-stochastic-simulation)
  - [4.1 Operator Model](#41-operator-model-operator_modelpy)
  - [4.2 Gillespie Engine](#42-gillespie-engine-gillespie_enginepy)
  - [4.3 Run Conditions](#43-run-conditions-run_conditionspy)
  - [4.4 Asymmetry Sweep](#44-asymmetry-sweep-asymmetry_sweeppy)
  - [4.5 Simulation Main](#45-simulation-main-simulation_mainpy)
  - [4.6 Phase 2 Self-Tests](#46-phase-2-self-tests)
- [5. Phase 3: Statistical Analysis](#5-phase-3-statistical-analysis)
  - [5.1 Noise Metrics](#51-noise-metrics-noise_metricspy)
  - [5.2 Bootstrap Confidence Intervals](#52-bootstrap-confidence-intervals-bootstrap_cipy)
  - [5.3 Statistical Tests](#53-statistical-tests-statistical_testspy)
  - [5.4 Sensitivity Analysis](#54-sensitivity-analysis-sensitivity_analysispy)
  - [5.5 Negative Controls](#55-negative-controls-negative_controlspy)
  - [5.6 Experimental Comparison](#56-experimental-comparison-experimental_comparisonpy)
  - [5.7 GMM / BIC Analysis](#57-gmm--bic-analysis)
  - [5.8 Analysis Orchestrator](#58-analysis-orchestrator-analysis_mainpy)
- [6. Phase 5: Thermodynamic Calibration](#6-phase-5-thermodynamic-calibration)
  - [6.1 Energy Calibration](#61-energy-calibration-energy_calibrationpy)
  - [6.2 Partition Function](#62-partition-function-partition_functionpy)
  - [6.3 MCMC Cooperativity Inference](#63-mcmc-cooperativity-inference-cooperativity_inferencepy)
  - [6.4 Operator Classification](#64-operator-classification-operator_classificationpy)
  - [6.5 DNA Shape](#65-dna-shape-dna_shapepy)
  - [6.6 Thermodynamic Main](#66-thermodynamic-main-thermodynamic_mainpy)
- [7. Phase 6: Environmental Extensions](#7-phase-6-environmental-extensions)
  - [7.1 Environmental Signals](#71-environmental-signals-environmental_signalspy)
  - [7.2 Two-Species Model](#72-two-species-model-two_species_modelpy)
  - [7.3 Two-Species Gillespie](#73-two-species-gillespie-two_species_gillespiepy)
  - [7.4 Mutual Information](#74-mutual-information-mutual_informationpy)
  - [7.5 Persistence Threshold](#75-persistence-threshold-persistence_thresholdpy)
  - [7.6 Run Environmental Conditions](#76-run-environmental-conditions-run_environmental_conditionspy)
  - [7.7 Environmental Main](#77-environmental-main-environmental_mainpy)
  - [7.8 Phase 6 Self-Tests](#78-phase-6-self-tests)

---

## 1. Project Overview

### 1.1 Description

This project quantifies the gene expression noise produced by the asymmetric Mce3R operator in *Mycobacterium tuberculosis* and characterizes its temporal dynamics and therapeutic quenching response. The Mce3R operator has two binding sites with a 20.4-fold affinity difference (Kd = 2.4 nM for the strong site vs. Kd = 49 nM for the weak site, measured by Panagoda et al. 2024 via electrophoretic mobility shift assays). The project uses stochastic simulation (Gillespie SSA), thermodynamic modeling (partition functions with MCMC-inferred cooperativity), and environmental perturbation analysis to predict that asymmetry generates approximately 19% more expression noise (CV) than an equivalent symmetric architecture, and that this noise may contribute to antibiotic persistence through phenotypic heterogeneity.

### 1.2 Central Hypothesis

The 20.4-fold binding site asymmetry of the Mce3R operator generates quantifiably more gene expression noise than symmetric alternatives, and this noise exhibits architecture-dependent temporal properties that may contribute to antibiotic persistence.

### 1.3 Aims and Stretch Goals

| Aim | Description | Phases |
|-----|-------------|--------|
| **Aim 1** | Discover Mce3R binding sites genome-wide via cross-species motif discovery (MEME) and genome scanning (FIMO) | Phase 1 |
| **Aim 2** | Model repression quantitatively using a 4-state operator model with Gillespie SSA and thermodynamic partition functions | Phases 2, 5 |
| **Aim 3a** | Characterize noise amplitude: does asymmetry generate more noise than symmetry? | Phases 3, 6 |
| **Aim 3b** (stretch) | Characterize temporal dynamics: autocorrelation time, dwell times, power spectral density | Phase 8 |
| **Aim 3c** (stretch) | Predict noise quenching: symmetrization sweep, IC50 of noise, dose-response | Phase 9 |
| **Visualization** | Publication-quality figures 1-19 | Phases 4, 7, 10 |

### 1.4 The 10 Phases and Dependency Graph

```
Phase 1 (motif discovery) ──────────────────────────┐
                                                     │
Phase 2 (Gillespie SSA) ───────┐                    │
                                ├── Phase 3 (stats)  │
Phase 2 (Gillespie SSA) ───────┘        │           │
                                         │           │
                                Phase 4 (figs 1-8)   │
                                                     │
                                    Phase 5 (thermo) ←┘
                                         │
                                Phase 6 (environmental)
                                         │
                                Phase 7 (figs 9-14)
                                         │
                                Phase 8 (temporal dynamics)
                                         │
                                Phase 9 (noise quenching)
                                         │
                                Phase 10 (figs 15-19)
```

**Execution model:**
- Phases 1 and 2 run **in parallel** as subprocesses from `main.py`
- Phase 3 runs **sequentially** after Phase 2 completes (requires Phase 2 NPZ outputs)
- Phase 4 runs **sequentially** after Phase 3 (requires `noise_metrics.csv`, `bootstrap_results.csv`)
- Phase 5 loads Phase 2 `phase2_summary.json` (for MCMC fitting targets) and Phase 1 `meme.txt` (for PWM)
- Phase 6 loads Phase 5 `mcmc_summary.json` (omega, dG_spacer medians) and Phase 2 `condition_D.npz` (for threshold)
- Phase 8 loads Phase 6 trace data
- Phase 9 uses Phase 2 Gillespie engine with Phase 5 calibrated parameters

### 1.5 What Is Novel vs. Established

**Already established -- do NOT claim as novel:**

| Claim | Citation |
|-------|----------|
| Cis-regulatory architecture shapes noise | Chowdhury et al. 2021, *Frontiers in Genetics* |
| Multiple binding sites modulate noise | Lengyel & Morelli 2017, *Physical Review E* |
| Noise drives persistence in M. tuberculosis | Quigley & Lewis 2022, *Microbiology Spectrum* |
| Mce3R deletion increases persister frequency | Pandey et al. 2023, *Research in Microbiology* |
| Mce3R has an asymmetric operator | Panagoda et al. 2024, *ACS Chemical Biology* |

**Novel contributions of this project:**

1. First quantitative noise characterization of a structurally resolved asymmetric bacterial operator with real Kd values
2. First prediction of noise temporal memory (autocorrelation time, dwell time) from operator architecture
3. First therapeutic dose-response curve for noise quenching in a TB regulatory system

---

## 2. Configuration -- Every Parameter

All parameters live in a single file: `config/parameters.py`. This is the single source of truth.

### 2.1 Primary Parameters Table

| Parameter | Value | Unit | Source | Phases Used |
|-----------|-------|------|--------|-------------|
| `H37Rv_accession` | `'NC_000962.3'` | -- | NCBI RefSeq | 1 |
| `H37Rv_genome_size` | `4_411_532` | bp | NCBI | 1 |
| `H37Rv_gc_content` | `0.656` | fraction | NCBI | 1, 5 |
| `M_bovis_accession` | `'NC_002945.4'` | -- | NCBI RefSeq | 1 |
| `M_marinum_accession` | `'NC_010612.1'` | -- | NCBI RefSeq | 1 |
| `yrbE3A_bovis_locus` | `'Mb1997'` | -- | NCBI annotation | 1 |
| `yrbE3A_marinum_locus` | `'MMAR_2522'` | -- | NCBI annotation | 1 |
| `Kd_strong` | `2.4` | nM | Panagoda et al. 2024, Table 1, Probe A | 2, 3, 5, 6 |
| `Kd_weak` | `49.0` | nM | Panagoda et al. 2024, Table 1, Probe C | 2, 3, 5, 6 |
| `k_on` | `0.0167` | nM^-1 min^-1 | Stormo & Zhao 2010 (10^6 M^-1 s^-1 x 60 s/min x 10^-9 M/nM) | 2, 3, 5, 6 |
| `cell_volume_fL` | `1.0` | femtoliters | Standard Mtb cell volume | 2, 6 |
| `nM_per_molecule` | `1.66` | nM/molecule | 1/(Avogadro x 1e-15 L) x 1e9 | 2, 6 |
| `k_max` | `0.15` | mRNA/min | Estimated from Mtb transcriptomics | 2, 3, 5, 6 |
| `block_strong` | `0.85` | fraction | Estimated | 2, 3, 5, 6 |
| `block_weak` | `0.50` | fraction | Estimated | 2, 3, 5, 6 |
| `k_translation` | `0.5` | protein/mRNA/min | Taniguchi et al. 2010, adjusted | 2, 3, 5, 6 |
| `t_half_mRNA` | `9.5` | min | Rustad et al. 2013, *NAR* | 2, 3 |
| `t_half_protein` | `1500.0` | min (~25 hr) | Dilution-dominated | 2, 3 |
| `n_cells_main` | `50_000` | cells | -- | 2 |
| `n_cells_sweep` | `5_000` | cells | -- | 2 |
| `t_max` | `30_000` | min (500 hr) | ~25 doubling times, 20x protein half-life | 2, 6 |
| `t_burn_in` | `15_000` | min | Discard first half; 10x protein half-life | 2, 6 |
| `n_trace_cells` | `20` | cells | -- | 2 |
| `trace_interval` | `10.0` | min | -- | 2, 6 |
| `trace_max_points` | `1_600` | points | (t_max - t_burn_in)/trace_interval = 1500 + margin | 2, 6 |
| `master_seed` | `42` | -- | -- | 2, 3, 6 |
| `gmm_random_state` | `42` | -- | -- | 3 |
| `bootstrap_seed` | `123` | -- | -- | 3 |
| `sweep_ratios` | `[1, 2, 5, 8, 10, 15, 20, 30, 40, 50]` | -- | -- | 2, 6 |
| `mce3r_actual_ratio` | `20.4` | -- | Kd_weak/Kd_strong = 49.0/2.4 | 2 |
| `Kd_symmetric` | `sqrt(2.4 * 49.0) = 10.8395...` | nM | Geometric mean | 2, 5 |
| `block_symmetric` | `sqrt(0.85 * 0.50) = 0.6519...` | fraction | Geometric mean | 2, 5 |
| `Kd_weak_disabled` | `1e12` | nM | Effectively infinite -- weak site never occupied | 2 |
| `k_on_disabled` | `0.0` | nM^-1 min^-1 | No binding at all | 2 |
| `sensitivity_range` | `0.50` | fraction | +/- 50% variation | 3 |
| `sensitivity_steps` | `10` | -- | Steps per parameter | 3 |
| `sensitivity_cells` | `1_000` | cells | Cells per step | 3 |
| `n_bootstrap` | `1_000` | -- | Bootstrap resamples | 3 |
| `ci_level` | `0.95` | -- | 95% confidence interval | 3 |
| `n_gmm_components` | `[1, 2, 3]` | -- | Components to fit | 3 |
| `gmm_reg_covar` | `1e-3` | -- | Regularization covariance | 3, 5 |
| `gmm_n_init` | `5` | -- | Number of GMM initializations | 3, 5 |
| `shuffle_n` | `1_000` | -- | Shuffle test permutations | 3 |
| `shuffle_alpha` | `0.01` | -- | Significance threshold for shuffle test | 3 |
| `meme_mod` | `'zoops'` | -- | MEME motif occurrence model | 1 |
| `meme_minw` | `6` | bp | MEME minimum motif width | 1 |
| `meme_maxw` | `110` | bp | MEME maximum motif width | 1 |
| `meme_nmotifs` | `5` | -- | MEME number of motifs to find | 1 |
| `fimo_thresh` | `1e-4` | -- | FIMO p-value threshold | 1 |
| `upstream_length` | `200` | bp | Upstream region extraction length | 1 |
| `known_spacer_bp` | `53` | bp | Two 25bp sites separated by 53bp | 1 |
| `weak_site_start_in_operator` | `10` | bp | Approximate start of weak site in 123bp operator | 1 |
| `strong_site_start_in_operator` | `88` | bp | Approximate start of strong site in 123bp operator | 1 |
| `mast_ev_threshold` | `100` | -- | MAST E-value threshold | 1 |
| `mast_mt_threshold` | `1e-4` | -- | MAST motif p-value threshold | 1 |
| `operator_length` | `123` | bp | Full operator length | 1, 5 |
| `operator_region_h37rv_start` | `2_207_477` | bp | Approximate intergenic region start | 1 |
| `operator_region_h37rv_end` | `2_207_699` | bp | Approximate intergenic region end | 1 |
| `mce3R_gene` | `'Rv1963c'` | -- | Mce3R gene identifier | 1 |
| `yrbE3A_gene` | `'Rv1964'` | -- | yrbE3A gene identifier | 1 |
| `mce3_operon_genes` | `['Rv1964', 'Rv1965', ..., 'Rv1977']` | -- | 14 genes Rv1964-Rv1977 | 1 |
| `regulon_genes_1` | `['Rv1933c', 'Rv1934c', 'Rv1935c']` | -- | First regulon group | 1 |
| `regulon_genes_2` | `['Rv1936', 'Rv1937', ..., 'Rv1941']` | -- | Second regulon group (Rv1936-Rv1941) | 1 |

**Operator sequence** (123 bp, from PDB 9B7Y, Panagoda 2024):

```
GCCCCGCGCTATAGGATACTAGCAAGATACATCATAGCCAATATATGCCAGT
TTGCATTGCTATTTACCGATCAGTTGTCCAAGCAATCGCGTATTGGCTATG
GACATCAGCGGTTCTGCCGC
```

### 2.2 Derived Parameters

All derived parameters are computed at module load time in `config/parameters.py`:

```python
gamma_mRNA    = ln(2) / t_half_mRNA      = 0.6931 / 9.5    = 0.07296 min^-1
gamma_protein = ln(2) / t_half_protein    = 0.6931 / 1500.0 = 0.0004621 min^-1
burst_size    = k_translation / gamma_mRNA = 0.5 / 0.07296  = 6.852
expected_fano_unregulated = 1 + burst_size = 1 + 6.852      = 7.852

k_off_strong  = Kd_strong * k_on = 2.4 * 0.0167  = 0.04008 min^-1
k_off_weak    = Kd_weak * k_on   = 49.0 * 0.0167 = 0.8183 min^-1
```

**Condition-level seeds:**
- Condition A: `master_seed + 0` = 42
- Condition B: `master_seed + 1` = 43
- Condition C: `master_seed + 2` = 44
- Condition D: `master_seed + 3` = 45
- Sweep ratio index `i`: `master_seed + 100 + i`

### 2.3 k_txn Array Computation

The transcription rate array has 4 elements, one per operator state. It is computed as a `float64[4]` numpy array:

```python
k_txn = np.array([
    k_max,                                                 # State 0: both empty
    k_max * (1 - block_strong),                            # State 1: strong occupied
    k_max * (1 - block_weak),                              # State 2: weak occupied
    k_max * (1 - block_strong) * (1 - block_weak),         # State 3: both occupied
])
```

With default values:

```
State 0 (unbound):     k_txn = 0.150     mRNA/min
State 1 (strong only): k_txn = 0.150 * (1 - 0.85) = 0.150 * 0.15 = 0.0225 mRNA/min
State 2 (weak only):   k_txn = 0.150 * (1 - 0.50) = 0.150 * 0.50 = 0.0750 mRNA/min
State 3 (both bound):  k_txn = 0.150 * 0.15 * 0.50 = 0.01125 mRNA/min
```

**Ordering invariant:** `k_txn[0] > k_txn[2] > k_txn[1] > k_txn[3]` (verified by self-test).

### 2.4 Phase 5-7 Parameters

These are defined as module-level constants in `config/parameters.py`, outside the PARAMS dict:

| Parameter | Value | Unit | Purpose |
|-----------|-------|------|---------|
| `TEMPERATURE_K` | `310.0` | K | Physiological temperature (37 C) |
| `kT_KCAL` | `0.001987 * 310.0 = 0.61597` | kcal/mol | Thermal energy |
| `GC_CONTENT` | `0.656` | fraction | From `PARAMS['H37Rv_gc_content']` |
| `BACKGROUND_FREQ` | `{'A': 0.172, 'C': 0.328, 'G': 0.328, 'T': 0.172}` | -- | `(1-GC)/2` and `GC/2` |
| `MCMC_N_WALKERS` | `32` | -- | emcee walkers |
| `MCMC_N_STEPS` | `5000` | -- | MCMC total steps |
| `MCMC_N_BURN` | `1000` | -- | MCMC burn-in steps |
| `MCMC_SEED` | `42` | -- | MCMC random seed |
| `STRONG_SITE_START` | `98` | bp (0-indexed) | Strong site start in operator |
| `STRONG_SITE_END` | `123` | bp (0-indexed, exclusive) | Strong site end in operator |
| `WEAK_SITE_START` | `0` | bp (0-indexed) | Weak site start in operator |
| `WEAK_SITE_END` | `25` | bp (0-indexed, exclusive) | Weak site end in operator |
| `ENV_COLORS` | `{'baseline': '#6B7280', 'cholesterol': '#F59E0B', 'acidic_pH': '#EF4444', 'host_like': '#7C3AED'}` | hex | Figure colors |

**Note:** The energy_calibration.py module uses different site positions that were refined from FIMO alignment:
- Strong site: `op[74:95]` (GTTGTCCAAGCAATCGCGTAT) -- corresponds to MEME-1
- Weak site: `op[51:72]` (TTTGCATTGCTATTTACCGAT) -- corresponds to MEME-2

---

## 3. Phase 1: Motif Discovery and Genome-Wide Scanning

**Directory:** `phase1_pipeline/`
**Orchestrator:** `pipeline_main.py`
**Output directory:** `results/phase1/`

### 3.1 Genome Acquisition (`download_genomes.py`)

**Purpose:** Download reference genomes and annotations from NCBI for three mycobacterial species.

**Module-level constants:**
```python
Entrez.email = "mce3r_project@example.com"
GENOME_DIR = os.path.join(PROJECT_ROOT, 'data', 'genomes')
```

**Downloads (5 files):**

| File | Accession | Format | API Call | Output Path |
|------|-----------|--------|----------|-------------|
| H37Rv genome | NC_000962.3 | FASTA | `Entrez.efetch(db="nucleotide", rettype="fasta")` | `data/genomes/H37Rv.fasta` |
| H37Rv annotation | NC_000962.3 | GenBank | `Entrez.efetch(rettype="gbwithparts")` | `data/genomes/H37Rv.gb` |
| M. bovis genome | NC_002945.4 | FASTA | `Entrez.efetch(rettype="fasta")` | `data/genomes/M_bovis.fasta` |
| M. marinum genome | NC_010612.1 | FASTA | `Entrez.efetch(rettype="fasta")` | `data/genomes/M_marinum.fasta` |
| H37Rv GFF3 | NC_000962.3 | GFF3 | `Entrez.efetch(rettype="gff3")` | `data/genomes/H37Rv.gff` |

**Functions:**

```python
def download_genome_fasta(accession: str, outpath: str, retries: int = 3) -> None
def download_genome_genbank(accession: str, outpath: str, retries: int = 3) -> None
def download_gff(accession: str, outpath: str, retries: int = 3) -> None
def run_download() -> dict  # Returns {'h37rv_fasta': path, 'h37rv_gb': path, ...}
```

**Caching logic:** If file exists and `os.path.getsize(outpath) > 1000`, skip download.

**Retry logic:** 3 attempts with 3-second sleep between retries.

**Self-tests (when run as `__main__`):**
- Sanity: All 5 output files exist and > 100 bytes
- Sanity: H37Rv FASTA parseable by `SeqIO.read(path, 'fasta')`
- Sanity: GFF contains > 100 CDS lines
- Science: H37Rv genome size approximately 4,411,532 bp
- Science: H37Rv GC content approximately 0.656

### 3.2 Upstream Sequence Extraction (`extract_upstream.py`)

**Purpose:** Extract 200 bp upstream of every CDS in H37Rv, plus ortholog upstream regions from M. bovis and M. marinum for cross-species MEME input.

**Key functions:**

```python
def load_genome(fasta_path: str) -> SeqRecord
def extract_circular(seq_str: str, start: int, length: int, genome_len: int) -> str
def reverse_complement(seq_str: str) -> str
def parse_gff_cds(gff_path: str) -> list[dict]
def deduplicate_cds(features: list) -> list
def extract_upstream_regions(features, genome_seq, upstream_len, genome_len) -> list[tuple]
def write_fasta(sequences: list[tuple], outpath: str) -> None
def find_ortholog_upstream(genome_seq, genome_len, genbank_path, target_locus_tags, upstream_len) -> tuple | None
def run_extract_upstream() -> dict
```

**How orthologs are identified:**

Locus tags used for substring matching against GenBank feature qualifiers (`locus_tag`, `old_locus_tag`, `gene`):
- H37Rv: `Rv1964` (yrbE3A gene)
- M. bovis: `Mb1997`
- M. marinum: `MMAR_2522`

The matching is case-insensitive substring matching: `target_lower in lt.lower()`.

**Extraction logic (200 bp upstream of start codon):**

For plus-strand genes:
```python
up_start = feat['start0'] - upstream_len  # 0-based
if up_start < 0:
    up_start = up_start % genome_len  # circular wrapping
seq = extract_circular(genome_seq, up_start, upstream_len, genome_len)
```

For minus-strand genes:
```python
up_start = feat['end0']  # upstream is AFTER the end on the minus strand
seq = extract_circular(genome_seq, up_start, upstream_len, genome_len)
seq = reverse_complement(seq)
```

**GFF coordinate conversion:** GFF3 uses 1-based inclusive coordinates. Converted to 0-based immediately:
```python
start0 = start_1based - 1
end0 = end_1based       # 0-based exclusive end (equals 1-based end)
```

**Circular genome handling:**
```python
def extract_circular(seq_str, start, length, genome_len):
    start = start % genome_len
    end = start + length
    if end <= genome_len:
        return seq_str[start:end]
    else:
        return seq_str[start:] + seq_str[:end - genome_len]
```

**Additional GenBank downloads:** If `M_bovis.gb` or `M_marinum.gb` do not exist, the module downloads them on-the-fly using `Entrez.efetch(rettype="gbwithparts")`.

**Outputs:**

| File | Contents |
|------|----------|
| `data/sequences/all_upstream_200bp.fasta` | ~3500-4500 sequences, each exactly 200 bp |
| `data/sequences/meme_input_orthologs.fasta` | 2-3 sequences (H37Rv + M. bovis + M. marinum yrbE3A orthologs) |
| `data/sequences/known_operator.fasta` | Single 123 bp sequence: the Mce3R operator from Panagoda 2024 |

### 3.3 Background Model

Generated by `run_meme.py` via `generate_background_model()`.

**Primary method:** If `fasta-get-markov` (MEME Suite) is available:
```bash
fasta-get-markov -dna <input_fasta> <bg_path>
```

**Fallback method:** Compute 0-order Markov model from H37Rv GC content:
```python
gc = PARAMS['H37Rv_gc_content']  # 0.656
at = 1.0 - gc                    # 0.344

# Written to background.model:
A  0.172000
C  0.328000
G  0.328000
T  0.172000
```

**Output file:** `results/phase1/background.model`

### 3.4 MEME Motif Discovery (`run_meme.py`)

**Purpose:** Discover conserved motifs in the yrbE3A ortholog upstream regions from 3 species.

**Critical design decision:** MEME input is the 3 ortholog sequences (`meme_input_orthologs.fasta`), NOT the full ~4000 upstream sequences.

**MEME availability check:** `shutil.which('meme')`

**Conda PATH fix:** The module prepends `sys.prefix + '/bin'` to `PATH` so conda-installed MEME is found.

**Real MEME command:**
```bash
meme <input_fasta> \
  -dna \
  -mod zoops \
  -revcomp \
  -minw 6 \
  -maxw 110 \
  -nmotifs 5 \
  -bfile <bg_path> \
  -oc <output_dir>
```

**Note:** The actual flags used in code are `-minw 6 -maxw 110` (from PARAMS), but the docstring mentions `-minw 20 -maxw 30 -nmotifs 3`. The code values from PARAMS take precedence. The `-pal` flag is explicitly NOT used (as stated in the module docstring).

**Return code handling:** MEME may return non-zero due to `meme_xml_to_html` Ghostscript issues. The code checks if `meme.txt` exists and is > 100 bytes even if the return code is non-zero.

**Mock output (when MEME not installed):**

Generates a minimal `meme.txt` with two motifs:
1. Motif 1: 25-bp PWM from the first 25 bp of the operator sequence (strong site). Consensus probability = 0.85 at each position.
2. Motif 2: 25-bp PWM from operator[49:74] (weak site region). Consensus probability = 0.80 at each position.

Background frequencies in mock output use the same GC-content-derived values.

**Output directory structure:**
```
results/phase1/meme_output/
  meme.txt      # Main output (parsed by Phase 5)
  meme.xml      # XML format (used by FIMO if real MEME ran)
  meme.html     # HTML summary (or mock placeholder)
  logo*.eps     # Sequence logos (real MEME only)
```

**Key function signatures:**
```python
def generate_background_model(input_fasta: str, bg_path: str) -> str
def run_real_meme(input_fasta: str, output_dir: str, bg_path: str) -> str
def generate_mock_meme_output() -> str
def run_meme_step() -> dict  # Returns {'meme_output_dir', 'meme_txt', 'background_model', 'mock_used'}
def parse_meme_motif_widths(meme_txt_path: str) -> list[int]
```

### 3.5 FIMO Genome-Wide Scanning (`run_fimo.py`)

**Purpose:** Scan the entire H37Rv genome for all occurrences of the discovered motifs.

**FIMO availability check:** `shutil.which('fimo')`

**Real FIMO command:**
```bash
fimo \
  --thresh 1e-4 \
  --oc <output_dir> \
  [--bfile <bg_path>]    # if background.model exists
  <meme_txt> <H37Rv.fasta>
```

**FIMO TSV parsing (`parse_fimo_tsv`):** Reads `fimo.tsv`, skipping lines starting with `#`. Extracts columns: `motif_id`, `sequence_name`, `start`, `stop`, `strand`, `score`, `p-value`, `q-value`, `matched_sequence`.

**Mock output (when FIMO not installed):** Generates 10 plausible 25-bp hits with known positions including:
- Hit at operator position 2,207,477 with p-value 1.2e-10
- Hit at operator position 2,207,500 with p-value 8.5e-11
- Additional hits near mce3 regulon genes

**Gene annotation (`annotate_and_filter_hits`):**

Parses H37Rv GFF for gene/CDS features. For each FIMO hit:
1. Finds nearest gene by minimum distance from hit position to gene boundaries
2. Computes distance to nearest gene
3. Determines if hit is intergenic (not overlapping any gene/CDS feature)

```python
def find_nearest_gene(pos: int, genes: list) -> tuple[tuple, int]
def is_intergenic(pos: int, genes: list) -> bool
```

**Output:** `results/phase1/predicted_sites.csv`

Columns: `rank`, `motif_id`, `sequence_name`, `start`, `stop`, `strand`, `score`, `p_value`, `q_value`, `matched_sequence`, `nearest_gene`, `nearest_gene_name`, `distance_to_gene`, `is_intergenic`

Rows are sorted by p-value (ascending = best first).

### 3.6 Conservation Analysis (`conservation_check.py`)

**Purpose:** Check conservation of the top 50 predicted binding sites in M. bovis and M. marinum genomes.

**Algorithm:** Sliding-window percent-identity search using numpy vectorization.

```python
def sliding_window_search(
    query: str,
    target_genome: str,
    min_identity: float = 0.80,
    min_coverage: float = 0.80
) -> tuple[float, int, str]
```

**Implementation details:**
1. Convert sequences to numpy byte arrays using `np.frombuffer(seq.encode('ascii'), dtype=np.uint8)`
2. For each strand (+ and reverse complement):
   - Build match count array by iterating over motif positions: `match_counts += (target_arr[j:j+n_positions] == q_arr[j]).astype(np.int32)`
   - Compute identities: `identities = match_counts / query_len`
   - Find best match: `idx = np.argmax(identities)`
3. Return best identity across both strands

**Conservation threshold:** Identity >= 0.80 (80%) = conserved.

**Input:** Top 50 sites from `predicted_sites.csv` (already sorted by p-value).

**Output:** `results/phase1/conservation_status.csv`

Columns: `rank`, `start`, `stop`, `strand`, `p_value`, `matched_sequence`, `nearest_gene`, `bovis_identity`, `bovis_position`, `bovis_strand`, `bovis_conserved`, `marinum_identity`, `marinum_position`, `marinum_strand`, `marinum_conserved`, `both_conserved`

### 3.7 Pipeline Orchestrator (`pipeline_main.py`)

**Step execution order:**
1. `download_genomes.run_download()`
2. `extract_upstream.run_extract_upstream()`
3. `run_meme.run_meme_step()`
4. `run_fimo.run_fimo_step()`
5. `conservation_check.run_conservation_check()`

Each step is wrapped in a try/except. Failures increment `n_sanity_fail` but do not halt subsequent steps.

**Mock tracking:** If any step's result dict contains `mock_used: True`, the pipeline sets `mock_used = True` and increments `n_warn`.

**Post-pipeline checks:**
1. Verify 5 key output files exist (H37Rv genome, all upstream, MEME output, predicted sites, conservation)
2. Check predicted_sites.csv contains a hit within 500 bp of the known operator region (2,207,477 - 2,207,699)
3. Check conservation_status.csv shows some sites conserved in M. bovis

**Output:** `results/phase1/phase1_summary.json` with keys: `status`, `n_sanity_pass`, `n_sanity_fail`, `n_science_expected`, `n_science_unexpected`, `n_warn`, `mock_used`, `wall_time_sec`, `output_files`

**Entry point function:**
```python
def run_pipeline() -> dict
```

Also aliased as `run_phase1 = run_pipeline` for import by `main.py`.

### 3.8 Phase 1 Self-Tests

**Pipeline-level sanity checks (5):**
1. Summary dict has all required keys
2. `phase1_summary.json` matches returned dict
3. All pipeline steps completed without sanity failures
4. Wall time < 30 minutes
5. Output files list >= 3 files

**Pipeline-level science checks (1):**
1. Some scientific expectations were met (n_science_expected > 0)

**Per-module self-tests (run via `__main__`):**
- download_genomes: 3 sanity + 2 science
- extract_upstream: 6 sanity + 2 science
- run_meme: 5 sanity + 2 science
- run_fimo: 5 sanity + 2 science
- conservation_check: 5 sanity + 2 science

---

## 4. Phase 2: Gillespie Stochastic Simulation

**Directory:** `phase2_simulation/`
**Orchestrator:** `simulation_main.py`
**Output directory:** `results/phase2/`

### 4.1 Operator Model (`operator_model.py`)

**Class:** `OperatorModel`

**Constructor:**
```python
class OperatorModel:
    def __init__(
        self,
        Kd_strong: float = None,   # defaults to PARAMS['Kd_strong']
        Kd_weak: float = None,     # defaults to PARAMS['Kd_weak']
        k_on: float = None,        # defaults to PARAMS['k_on']
        block_strong: float = None, # defaults to PARAMS['block_strong']
        block_weak: float = None,   # defaults to PARAMS['block_weak']
        k_max: float = None,       # defaults to PARAMS['k_max']
    )
```

**Computed attributes:**
```python
self.k_off_strong = self.Kd_strong * self.k_on  # 2.4 * 0.0167 = 0.04008 min^-1
self.k_off_weak   = self.Kd_weak   * self.k_on  # 49.0 * 0.0167 = 0.8183 min^-1
self.n_bound_lookup = np.array([0, 1, 1, 2], dtype=np.int64)
```

**4-state operator model:**

| State | Description | n_bound | k_txn formula |
|-------|-------------|---------|---------------|
| 0 | Both sites empty | 0 | `k_max` = 0.150 |
| 1 | Strong site occupied | 1 | `k_max * (1 - block_strong)` = 0.0225 |
| 2 | Weak site occupied | 1 | `k_max * (1 - block_weak)` = 0.075 |
| 3 | Both sites occupied | 2 | `k_max * (1 - block_strong) * (1 - block_weak)` = 0.01125 |

**Methods:**

```python
def get_transcription_rates(self) -> np.ndarray  # float64[4]
def get_binding_rate_constants(self) -> dict      # {'k_on', 'k_off_strong', 'k_off_weak'}
def get_numba_arrays(self) -> dict                # All arrays for Numba engine
```

`get_numba_arrays()` returns:
```python
{
    'k_txn': float64[4],
    'n_bound': int64[4],
    'k_on': float64,
    'k_off_strong': float64,
    'k_off_weak': float64,
}
```

### 4.2 Gillespie Engine (`gillespie_engine.py`)

**THE COMPLETE REACTION SYSTEM -- 12 reactions:**

| Reaction | Transition | Propensity | Condition |
|----------|------------|------------|-----------|
| 0 | 0 -> 1 (bind strong from empty) | `k_on * R_nM` | `op_state == 0` |
| 1 | 1 -> 0 (unbind strong) | `k_off_strong` | `op_state == 1` |
| 2 | 0 -> 2 (bind weak from empty) | `k_on * R_nM` | `op_state == 0` |
| 3 | 2 -> 0 (unbind weak) | `k_off_weak` | `op_state == 2` |
| 4 | 2 -> 3 (bind strong when weak occupied) | `k_on * R_nM` | `op_state == 2` |
| 5 | 3 -> 2 (unbind strong when both occupied) | `k_off_strong` | `op_state == 3` |
| 6 | 1 -> 3 (bind weak when strong occupied) | `k_on * R_nM` | `op_state == 1` |
| 7 | 3 -> 1 (unbind weak when both occupied) | `k_off_weak` | `op_state == 3` |
| 8 | mRNA synthesis | `k_txn[op_state]` | always |
| 9 | protein synthesis | `k_translation * mRNA` | always |
| 10 | mRNA decay | `gamma_mRNA * mRNA` | always |
| 11 | protein decay | `gamma_protein * protein` | always |

**How R_nM is computed (single-species model):**
```python
nb = n_bound[op_state]             # 0, 1, 1, or 2
free_protein = protein - nb
if free_protein < 0:
    free_protein = 0
R_nM = free_protein * nM_per_molecule   # 1.66 nM per molecule
```

**Note:** In the single-species model, the protein being simulated IS the Mce3R repressor. There is no separate "Mce3R concentration" parameter. The protein count itself, converted to nM, determines binding propensity. This means single-species is NOT truly constitutive -- it models a gene product that also acts as its own repressor.

**The SSA loop:**

```python
# 1. Compute all 12 propensities
a_total = sum(a[0:12])

if a_total <= 0.0:
    break

# 2. Draw time to next reaction (exponential waiting time)
r1 = np.random.random()
dt = -np.log(r1) / a_total
t += dt

# 3. Record trace if past burn-in and due
if record_trace == 1 and t >= t_burn_in:
    while next_trace_time <= t and n_trace < trace_max_points:
        trace_times[n_trace] = next_trace_time
        trace_proteins[n_trace] = float(protein)
        n_trace += 1
        next_trace_time += trace_interval

# 4. Select reaction via linear search
r2 = np.random.random() * a_total
cumsum = 0.0
reaction = -1
for i in range(12):
    cumsum += a[i]
    if cumsum >= r2:
        reaction = i
        break

# 5. Execute reaction (update state variables)
# ... switch on reaction number ...
```

**Fallback for reaction selection:** If linear search fails (numerical edge case), select the reaction with the largest propensity.

**Numba @njit decoration:**

```python
@numba.njit
def simulate_cell(
    k_txn, n_bound, k_on, k_off_strong, k_off_weak,
    k_translation, gamma_mRNA, gamma_protein,
    nM_per_molecule, t_max, t_burn_in,
    trace_interval, trace_max_points,
    seed, record_trace
) -> tuple[int, int, int, np.ndarray, np.ndarray, int]
```

All inputs are scalars or numpy arrays (no dicts, classes, or Python objects). This is required for Numba compatibility. First call triggers JIT compilation (~30 seconds on M4 MacBook Air).

**Initial conditions:** `op_state = 0`, `mRNA = 0`, `protein = 0`, `t = 0.0`

**Seed handling:** Each cell uses `np.random.seed(cell_seed)` where `cell_seed = master_seed + i` (cell index).

**Population runner:**

```python
def run_population(
    model_arrays: dict,
    n_cells: int,
    master_seed: int,
    record_traces: int = 0,
    n_trace_cells: int = 0
) -> dict
```

Returns: `{'proteins': int64[n_cells], 'mRNAs': int64[n_cells], 'op_states': int64[n_cells], 'traces': list}`

Traces are recorded only for the first `n_trace_cells` cells when `record_traces=1`.

### 4.3 Run Conditions (`run_conditions.py`)

**Condition A -- Asymmetric (wild-type):**
```python
OperatorModel()  # all defaults
seed = master_seed + 0  # = 42
n_cells = 50,000
# Records traces for n_trace_cells = 20 cells
```

**Condition B -- Symmetric:**
```python
OperatorModel(
    Kd_strong=10.84,        # Kd_symmetric = sqrt(2.4 * 49.0)
    Kd_weak=10.84,
    block_strong=0.652,     # block_symmetric = sqrt(0.85 * 0.50)
    block_weak=0.652
)
seed = master_seed + 1  # = 43
```

**Condition C -- Single-site (weak disabled):**
```python
OperatorModel(Kd_weak=1e12)  # effectively infinite Kd
seed = master_seed + 2  # = 44
```

**Condition D -- No regulation:**
```python
model = OperatorModel(k_on=0.0)
arrays = model.get_numba_arrays()
arrays['k_txn'] = np.array([k_max] * 4, dtype=np.float64)  # Override: all states = k_max
seed = master_seed + 3  # = 45
```

**Critical:** For Condition D, `k_txn` is explicitly overridden to `[k_max, k_max, k_max, k_max]`. The `OperatorModel` still computes block-reduced rates, but these are replaced. Additionally, `k_on = 0.0` means no binding reactions fire.

**Statistics computed per condition:**
```python
def compute_stats(proteins):
    p = proteins.astype(np.float64)
    mean = np.mean(p)
    std = np.std(p)
    var = np.var(p)
    cv = std / mean if mean > 0 else 0.0
    fano = var / mean if mean > 0 else 0.0
    return {'mean': mean, 'std': std, 'cv': cv, 'fano': fano}
```

**NPZ output per condition:**

Each condition saves to `results/phase2/condition_{A,B,C,D}.npz` containing:
- `proteins`: int64 array
- `mRNAs`: int64 array
- `op_states`: int64 array
- `trace_time_{i}`, `trace_prot_{i}`: float64 arrays for each trace cell (A, B, C only)
- `n_traces`: int array

### 4.4 Asymmetry Sweep (`asymmetry_sweep.py`)

**Ratios swept:** `[1, 2, 5, 8, 10, 15, 20, 30, 40, 50]` (10 values)

**How Kd values are computed for each ratio r:**

The geometric mean Kd is held constant:
```python
Kd_geo = sqrt(2.4 * 49.0) = 10.84 nM

Kd_strong = Kd_geo / sqrt(r)
Kd_weak   = Kd_geo * sqrt(r)
```

Verification: `Kd_strong * Kd_weak = Kd_geo^2 = 2.4 * 49.0 = 117.6` for all ratios.

**Block fraction co-variation via power-law interpolation:**

At ratio=1: both blocks = `block_geo = sqrt(0.85 * 0.50) = 0.652`
At ratio=r_wt (~20.4): recovers measured values (0.85, 0.50)

```python
r_wt = Kd_weak / Kd_strong  # 49.0 / 2.4 = 20.4167
log_r_wt = np.log(r_wt)

if r <= 1.0:
    bs = block_geo
    bw = block_geo
else:
    frac = min(1.0, np.log(r) / log_r_wt)
    bs = block_geo * (block_strong / block_geo) ** frac
    bw = block_geo * (block_weak / block_geo) ** frac
```

**Intermediate fraction thresholds:**
```python
unreg_mean_est = k_max * k_translation / (gamma_mRNA * gamma_protein)
# = 0.15 * 0.5 / (0.07296 * 0.0004621) = 2224.6
low_thresh  = unreg_mean_est * 0.10  # ~222
high_thresh = unreg_mean_est * 0.60  # ~1335
```

**Cells per ratio:** 5,000 (from `PARAMS['n_cells_sweep']`)

**Seed per ratio:** `master_seed + 100 + i` where i is the ratio index (0-9)

**Output:** `results/phase2/condition_E_sweep.npz` containing:
- `ratios`: float64 array
- `means`, `cvs`, `fanos`, `intermediate_fractions`: float64 arrays
- `proteins_ratio_{r}`: int64 array per ratio value

### 4.5 Simulation Main (`simulation_main.py`)

**Orchestration order:** D -> B -> C -> A -> E

Condition D runs first to validate the engine (unregulated = simplest case).

**Summary JSON (`results/phase2/phase2_summary.json`):**

```json
{
    "phase": 2,
    "description": "Mce3R stochastic simulation -- 4-state operator model",
    "total_time_seconds": ...,
    "parameters": { ... },  // All key PARAMS values
    "conditions": {
        "D": {"n_cells": 50000, "mean_protein": ..., "std_protein": ..., "cv": ..., "fano": ..., "elapsed_seconds": ...},
        "B": { ... },
        "C": { ... },
        "A": { ... },
        "E": {"n_ratios": 10, "n_cells_per_ratio": 5000, "ratios": [...], "means": [...], ...}
    },
    "predictions": {
        "cv_A_gt_cv_B": true,  // THE key biological prediction
        "cv_A": ...,
        "cv_B": ...,
        "cv_C": ...,
        "cv_D": ...,
        "fano_D_approx_expected": ...,
        "expected_fano_unregulated": 7.85
    }
}
```

**Entry points:**
```python
def run_all(n_cells_main=None, n_cells_sweep=None) -> dict
def run_phase2() -> dict  # Production wrapper called by main.py subprocess
```

### 4.6 Phase 2 Self-Tests

**Per-module self-tests:**

- operator_model.py: 10 sanity + 1 scientific
- gillespie_engine.py: 8 sanity + 3 scientific
- run_conditions.py: 10 sanity + 2 scientific
- asymmetry_sweep.py: 6 sanity + 1 scientific
- simulation_main.py: 12 sanity + 1 scientific

**Key scientific validations:**
1. CV(A) > CV(B) -- asymmetric operator produces more noise
2. Fano(D) near expected_fano_unregulated (7.85, within 50%)
3. Mean(D) > Mean(A) -- unregulated mean exceeds regulated mean
4. CV increases with Kd ratio in asymmetry sweep

---

## 5. Phase 3: Statistical Analysis

**Directory:** `phase3_analysis/`
**Orchestrator:** `analysis_main.py`
**Output directory:** `results/phase3/`

### 5.1 Noise Metrics (`noise_metrics.py`)

**Purpose:** Compute summary noise statistics and GMM fits for each condition.

**Function:**
```python
def compute_noise_metrics(proteins: np.ndarray, label: str) -> dict
def run_noise_metrics() -> pd.DataFrame
```

**Metrics computed per condition (A, B, C, D):**

| Metric | Formula |
|--------|---------|
| mean | `np.mean(proteins)` |
| variance | `np.var(proteins, ddof=1)` (sample variance) |
| CV | `sqrt(variance) / mean` |
| Fano | `variance / mean` |
| bimodality_coeff | `(skewness^2 + 1) / (kurtosis_excess + 3)` |

**Bimodality coefficient:** BC > 5/9 (~0.555) suggests bimodality. Uses `scipy.stats.skew()` and `scipy.stats.kurtosis(fisher=True)`.

**GMM fitting (per condition):**

For k in [1, 2, 3]:
```python
gmm = GaussianMixture(
    n_components=k,
    reg_covar=1e-3,       # PARAMS['gmm_reg_covar']
    n_init=5,             # PARAMS['gmm_n_init']
    random_state=42,      # PARAMS['gmm_random_state']
    max_iter=300,
)
```

Robustness: Adds dither noise (`std * 1e-4`) before fitting to prevent singular covariance. If fit fails, returns BIC = infinity.

Best model selected by lowest BIC.

**Output:** `results/phase3/noise_metrics.csv` with columns: `condition`, `n_cells`, `mean`, `variance`, `CV`, `Fano`, `bimodality_coeff`, `gmm_best_k`, `gmm_best_bic`

### 5.2 Bootstrap Confidence Intervals (`bootstrap_ci.py`)

**Algorithm:**
```python
def bootstrap_statistic(x, stat_func, n_bootstrap=1000, ci=0.95, seed=None):
    rng = np.random.default_rng(seed)
    point = stat_func(x)
    boot_vals = np.zeros(n_bootstrap)
    for i in range(n_bootstrap):
        idx = rng.integers(0, len(x), size=len(x))  # resample WITH replacement
        boot_vals[i] = stat_func(x[idx])
    lo = np.nanpercentile(boot_vals, 100 * (1-ci)/2)     # 2.5th percentile
    hi = np.nanpercentile(boot_vals, 100 * (1-(1-ci)/2)) # 97.5th percentile
    return point, lo, hi
```

**Metrics bootstrapped:** CV, Fano, intermediate_fraction, persister_fraction

**Seeds:**
- CV: `bootstrap_seed` = 123
- Fano: `bootstrap_seed + 1` = 124
- intermediate_fraction: `bootstrap_seed + 2` = 125
- persister_fraction: `bootstrap_seed + 3` = 126

**Persister threshold for bootstrap:** `mean(proteins_B) + 2 * std(proteins_B)`

**Output:** `results/phase3/bootstrap_results.csv` with columns: `condition`, `metric`, `point`, `CI_lower`, `CI_upper`

Total rows: 4 conditions x 4 metrics = 16 rows.

### 5.3 Statistical Tests (`statistical_tests.py`)

**Kolmogorov-Smirnov test:**

Pairs tested: (A,B), (A,C), (A,D), (B,D), (C,D)

```python
ks_stat, ks_p = scipy.stats.ks_2samp(data[c1], data[c2])
```

**Cohen's d:**

```python
def cohens_d(x1, x2):
    pooled_var = ((n1-1)*var1 + (n2-1)*var2) / (n1+n2-2)
    pooled_sd = sqrt(pooled_var)
    return (mean(x1) - mean(x2)) / pooled_sd
```

Uses `ddof=1` for sample variance.

**GMM Likelihood Ratio test:**

For each condition, compares 2-component vs 3-component and 1-component vs 2-component GMMs.

```python
def gmm_lr_test(x, k_null=2, k_alt=3):
    # LL_total = model.score(X) * N  (TOTAL log-likelihood, not average)
    lr_stat = 2 * (LL_alt_total - LL_null_total)
    # Approximate p-value via chi2 with df_diff=3
    p_value = 1.0 - scipy.stats.chi2.cdf(lr_stat, df=3)
```

**Output:** `results/phase3/statistical_tests.csv` with columns: `test`, `comparison`, `statistic`, `p_value`, `cohens_d`

Total rows: 5 KS + 4 GMM_LR_2v3 + 4 GMM_LR_1v2 = 13 rows.

### 5.4 Sensitivity Analysis (`sensitivity_analysis.py`)

**Parameters swept (7):**

```python
SENSITIVITY_PARAMS = ['k_max', 'k_translation', 't_half_mRNA', 't_half_protein',
                      'block_strong', 'block_weak', 'k_on']
```

**Note:** This is 7 parameters, not 8 as mentioned in some documentation. `Kd_strong` and `Kd_weak` are NOT directly swept here.

**Range:** +/- 50% (`sensitivity_range = 0.50`). Clamped for block fractions to [0.01, 0.99] and k_on to [1e-6, inf].

**Steps:** 10 per parameter (`sensitivity_steps = 10`)

**Cells:** 1,000 per step (`sensitivity_cells = 1_000`)

**For `t_half_mRNA` and `t_half_protein`:** Converted to gamma rates before simulation:
```python
gamma_mRNA = np.log(2) / t_half_mRNA
gamma_protein = np.log(2) / t_half_protein
```

These are temporarily monkey-patched into the PARAMS dict, then restored after each simulation.

**Output:** `results/phase3/sensitivity_data.csv` with columns: `parameter`, `base_value`, `test_value`, `fold_change`, `n_cells`, `mean`, `CV`, `n_modes`

Total rows: 7 params x 10 steps = 70 rows (production).

### 5.5 Negative Controls (`negative_controls.py`)

**Three control tests:**

**1. Shuffle test:** Mononucleotide-preserving shuffle of predicted binding sites.
- Takes up to 10 sites from the operator sequence (25 bp each, with 12-13 bp overlap)
- For each site: shuffle 1,000 times (`shuffle_n`), compute GC content as scoring metric
- p-value: fraction of shuffles with GC content >= real site
- Uses mock sites from operator sequence (not actual FIMO output)

```python
def mononucleotide_shuffle(seq, rng=None):
    chars = list(seq)
    rng.shuffle(chars)
    return ''.join(chars)
```

**2. Symmetric TetR control:** Run simulation with Kd_strong = Kd_weak = 5 nM, block_strong = block_weak = 0.85. Expect unimodal distribution (GMM best_k = 1).

```python
OperatorModel(Kd_strong=5.0, Kd_weak=5.0, block_strong=0.85, block_weak=0.85)
# seed = master_seed + 9000
# n_cells = 100 (default)
```

**3. Poisson baseline:** Verify Condition D Fano factor in [5, 12].

**Output:** `results/phase3/negative_controls.csv` with columns: `test`, `detail`, `value`, `p_value`, `pass`, `note`

### 5.6 Experimental Comparison (`experimental_comparison.py`)

**Fold-change prediction:**
```python
fold_change = mean_D / mean_A
# Published (Santangelo et al. 2009): 8.5x for yrbE3A
# Predicted: mean_D / mean_A (typically ~11.3x)
# Ratio: fold_change / 8.5 (typically ~1.33x)
```

**Persister fraction (primary method -- GMM component weights):**

Fit a 2-component GMM to Condition A protein distribution. The HIGH-expression component represents de-repressed cells (persisters). Persister fraction = weight of the high component.

```python
gmm = GaussianMixture(n_components=2, reg_covar=1e-3, n_init=5, random_state=42, max_iter=300)
gmm.fit(X_dither)
# Sort components by mean; high-mean component weight = persister fraction
```

**Persister fraction (secondary method -- 2-sigma threshold):**
```python
sigma2_thresh = mean(proteins_B) + 2.0 * std(proteins_B, ddof=1)
persister_frac_2sigma = mean(proteins_A > sigma2_thresh)
```

**GMM intersection threshold (reference):** Find the concentration where the two GMM components have equal probability by evaluating `predict_proba()` on a grid between the two component means.

**Output:** `results/phase3/experimental_comparison.csv` and updates `noise_metrics.csv` with persister columns.

### 5.7 GMM / BIC Analysis

Performed within `noise_metrics.py` as part of `compute_noise_metrics()`:

For each condition, fit 1, 2, and 3-component Gaussian mixtures:
```python
for k in [1, 2, 3]:
    gmm = GaussianMixture(
        n_components=k,
        reg_covar=1e-3,
        n_init=5,
        random_state=42,
        max_iter=300,
    )
    gmm.fit(X_dither)
    bic = gmm.bic(X_dither)
```

Select best by lowest BIC. Result saved as `gmm_best_k` and `gmm_best_bic` in `noise_metrics.csv`.

### 5.8 Analysis Orchestrator (`analysis_main.py`)

**Step execution order:**
1. Noise Metrics
2. Bootstrap CIs
3. Sensitivity Analysis
4. Statistical Tests
5. Experimental Comparison
6. Negative Controls

```python
def run_phase3(self_test_mode: bool = False) -> dict
```

When `self_test_mode=True`, sensitivity analysis uses `n_steps=3, n_cells=10` and negative controls use `n_shuffle=100, n_cells_symmetric=20`.

**Output:** `results/phase3/phase3_summary.json` with standard summary dict keys.

---

## 6. Phase 5: Thermodynamic Calibration

**Directory:** `phase5_thermodynamic/`
**Orchestrator:** `thermodynamic_main.py`
**Output directory:** `results/phase5/`

### 6.1 Energy Calibration (`energy_calibration.py`)

**Purpose:** Convert MEME PWM to per-position binding free energies and calibrate against experimental Kd values.

**Berg-von Hippel model:**

For each position i and base b in the PWM:
```
epsilon(i, b) = -kT * ln((pwm[i,b] + pseudocount) / bg[b])
```

Where:
- `pwm[i,b]` = frequency at position i, base b (ACGT order)
- `bg[b]` = background frequency for base b
- `kT` = 0.616 kcal/mol at T = 310 K
- `pseudocount` = 0.001 (to avoid log(0))

After adding pseudocount, rows are renormalized: `pwm_norm = (pwm + pseudocount) / row_sums`

**Scoring a sequence against the energy matrix:**
```python
def score_sequence(seq, energy_matrix) -> float:
    # If len(seq) == motif_len: sum energies at each position
    # If len(seq) > motif_len: sliding window, return lowest (most favorable) score
    # If len(seq) < motif_len: raise ValueError
```

Unknown bases ('N') use the mean energy at that position.

**Calibration procedure:**

The operator has two binding sites corresponding to **different** MEME motifs, so a single-offset approach is physically incorrect. The code uses **two independent offsets**:

```python
# Experimental DeltaG from Kd:
dG_strong_exp = kT * ln(Kd_strong * 1e-9)  # e.g., 0.616 * ln(2.4e-9) = -12.23 kcal/mol
dG_weak_exp   = kT * ln(Kd_weak * 1e-9)    # e.g., 0.616 * ln(49e-9) = -10.37 kcal/mol

# Raw PWM scores:
dG_strong_raw = score_sequence(strong_seq, energy_matrix)
dG_weak_raw   = score_sequence(weak_seq, energy_matrix)

# Independent offsets:
offset_strong = dG_strong_exp - dG_strong_raw
offset_weak   = dG_weak_exp - dG_weak_raw
offset_mean   = (offset_strong + offset_weak) / 2.0

# Calibrated DeltaG values (exact by construction):
dG_strong = dG_strong_exp
dG_weak   = dG_weak_exp
```

**Binding site positions in operator (0-indexed):**
- Strong site: `op[74:95]` = GTTGTCCAAGCAATCGCGTAT (matches MEME-1, Kd = 2.4 nM)
- Weak site: `op[51:72]` = TTTGCATTGCTATTTACCGAT (matches MEME-2, Kd = 49 nM)

**Predicting Kd for unknown sequences:**
```python
def predict_Kd_for_sequence(seq, calibration) -> float:
    raw_score = score_sequence(seq, energy_matrix)
    dG_calib = raw_score + offset_mean     # uses mean offset (conservative)
    Kd_pred = exp(dG_calib / kT) * 1e9     # convert back to nM
    return Kd_pred
```

**Convenience builder:**
```python
def build_calibration(meme_path=MEME_TXT) -> dict
```

Returns a dict with all calibration info including `energy_matrix`, both dG values, both offsets, predicted Kd values, etc.

### 6.2 Partition Function (`partition_function.py`)

**4-state statistical mechanics model:**

```
Z = w_U + w_S + w_W + w_D

w_U = 1                                           (unbound)
w_S = exp(-beta * (dG_strong - mu))                (strong only)
w_W = exp(-beta * (dG_weak - mu))                  (weak only)
w_D = omega * exp(-beta * (dG_strong + dG_weak - 2*mu + dG_spacer))  (both bound)

where beta = 1/kT
```

**State probabilities:**
```
P_i = w_i / Z
```

**Chemical potential from concentration:**
```python
def mu_from_concentration(conc_nM, kT=kT_KCAL):
    return kT * ln(conc_nM * 1e-9)
```

At `[R] = Kd`, `mu = kT * ln(Kd_M) = dG_site`, so `P(bound) = 0.5`. Correct by construction.

**Mean transcription rate:**
```python
k_mean = sum_i P_i * k_txn_i
# where k_txn_i are the same 4-state rates as in the Gillespie model
```

**Fold repression:**
```
fold = k_max / k_mean
```

**Architecture comparison (`compare_architectures`):**

4 architectures compared:

| Architecture | dG_strong | dG_weak |
|-------------|-----------|---------|
| native_asymmetric | dG_s (calibrated) | dG_w (calibrated) |
| symmetric | (dG_s + dG_w) / 2 | (dG_s + dG_w) / 2 |
| strong_only | dG_s | dG_s |
| single_site | dG_s | +100 kcal/mol (inert) |

**Key function signatures:**
```python
def partition_function(dG_strong, dG_weak, mu, omega=1.0, dG_spacer=0.0, kT=kT_KCAL) -> float
def state_probabilities(dG_strong, dG_weak, mu, omega=1.0, dG_spacer=0.0, kT=kT_KCAL) -> np.ndarray  # shape (4,)
def mean_transcription_rate(...) -> float
def repression_fold(...) -> float
def mu_from_concentration(conc_nM, kT=kT_KCAL) -> float
def repression_curve(concentrations_nM, dG_strong, dG_weak, omega=1.0, dG_spacer=0.0, ...) -> dict
def compare_architectures(concentrations_nM, calibration, omega=1.0, dG_spacer=0.0, ...) -> dict
```

### 6.3 MCMC Cooperativity Inference (`cooperativity_inference.py`)

**Purpose:** Use Bayesian inference (emcee MCMC) to estimate cooperativity factor omega and spacer penalty dG_spacer from Phase 2 simulation outputs.

**Free parameters:** `theta = [ln(omega), dG_spacer]`

**emcee EnsembleSampler settings:**
- Walkers: 32 (`MCMC_N_WALKERS`)
- Total steps: 5,000 (`MCMC_N_STEPS`)
- Burn-in: 1,000 (`MCMC_N_BURN`)
- Seed: 42 (`MCMC_SEED`)
- Dimensionality: 2

**Prior (`log_prior`):**
```python
# Gaussian priors:
lp_ln_omega  = -0.5 * (ln_omega / 2.0)^2   # N(0, 2)
lp_dG_spacer = -0.5 * (dG_spacer / 2.0)^2  # N(0, 2)

# Hard bounds:
|dG_spacer| < 10 kcal/mol
ln(omega) in [-5, 5]  # omega in [0.007, 148]
```

**Likelihood (`log_likelihood`):**

Compares predicted mean protein against observed values from Phase 2 summary for Conditions A and B:

```python
# Condition A (asymmetric):
pred_A = mean_k_txn_A * k_translation / (gamma_mRNA * gamma_protein)
sigma_A = max(0.05 * obs_A, 1.0)  # 5% of observed mean
logL_A = -0.5 * ((obs_A - pred_A) / sigma_A)^2

# Condition B (symmetric): uses mean dG at both sites + block_symmetric
pred_B = mean_k_txn_B * k_translation / (gamma_mRNA * gamma_protein)
sigma_B = max(0.05 * obs_B, 1.0)
logL_B = -0.5 * ((obs_B - pred_B) / sigma_B)^2

return logL_A + logL_B
```

**Mce3R concentration:** Fixed at 332 nM (~200 molecules in 1 fL cell). This is typical for TetR-family repressors.

**Initial positions:**
```python
p0 = rng.normal(0.0, 0.1, size=(n_walkers, 2))
```

**Convergence diagnostics:**
- R-hat via `arviz.rhat()` (falls back to manual Gelman-Rubin if arviz fails)
- Acceptance fraction: `np.mean(sampler.acceptance_fraction)`

**Output files:**
- `results/phase5/mcmc_posteriors.npz`: `flat_chain`, `chain`, `omega_samples`, `dG_spacer_samples`
- `results/phase5/mcmc_summary.json`: medians, 16th/84th percentile CIs, R-hat values, acceptance fraction

**Data dict builder:**
```python
def build_data_dict(phase2_summary_path=None, calibration=None) -> dict
```

Loads Phase 2 summary JSON and combines with calibration data and kinetic parameters.

### 6.4 Operator Classification (`operator_classification.py`)

**Purpose:** For each FIMO-predicted binding site, compute a repression curve and classify the operator type.

**Classification procedure per site:**

1. Score the matched_sequence against the PWM energy matrix using `score_sequence()`
2. Apply mean calibration offset to get `dG_site`
3. Use `dG_site` as `dG_strong` and the calibrated `dG_weak` as the paired site
4. Compute repression curve over 0.01-10,000 nM (200 log-spaced points)
5. Also compute repression curve with only the strong site (weak = +100 kcal/mol)
6. Fit Hill equation to the dual-site curve

**Hill equation:**
```python
def _hill_equation(log_conc, log_K_half, n_H, fold_max):
    conc = 10^log_conc
    K_half = 10^log_K_half
    return fold_max / (1 + (K_half / conc)^n_H)
```

Fitted using `scipy.optimize.curve_fit` with bounds: `log_K_half in [-3, 6]`, `n_H in [0.1, 10]`, `fold_max in [1, 1e6]`.

**Classification rules:**

| Condition | Classification |
|-----------|---------------|
| n_H > 2.0 | `digital_switch` |
| n_H < 1.3 AND weak_site_contribution < 0.5 | `graded_repressor` |
| weak_site_contribution > 0.5 AND n_H >= 1.3 | `noise_modulator` |
| else | `graded_repressor` |

Where `weak_site_contribution = mean(abs(folds_both - folds_single))`.

**Output:** `results/phase5/operator_classifications.csv`

Columns: `rank`, `motif_id`, `sequence_name`, `matched_sequence`, `fimo_score`, `p_value`, `raw_pwm_score`, `dG_site`, `n_H`, `K_half_nM`, `fold_max`, `dynamic_range`, `transition_width`, `weak_site_contribution`, `classification`, `fit_success`, `nearest_gene`, `is_intergenic`

### 6.5 DNA Shape (`dna_shape.py`)

**Status:** Skipped if DNAshapeR R package is not installed.

**Would compute:** Minor groove width, propeller twist, roll from the operator sequence using the DNAshapeR R package.

The `compute_dna_shape()` function returns `{'status': 'skipped', 'reason': '...'}` when the R package is unavailable.

### 6.6 Thermodynamic Main (`thermodynamic_main.py`)

**Step execution order:**

1. Load MEME PWM and calibrate energies -> `energy_parameters.json`
2. Compute repression curves for 4 architectures -> `repression_curves.csv`
3. Run MCMC cooperativity inference -> `mcmc_posteriors.npz`, `mcmc_summary.json`
4. Classify top-20 FIMO operators -> `operator_classifications.csv`
5. (Optional) DNA shape features
6. Write `phase5_summary.json`

**Entry point:**
```python
def run_phase5(
    results_dir=RESULTS_DIR,
    n_walkers=None, n_steps=None, n_burn=None, seed=None,
    top_n=20,
    run_mcmc_flag=True,
    verbose=True,
) -> dict
```

**Repression curves output:** `results/phase5/repression_curves.csv` with columns: `architecture`, `concentration_nM`, `k_mean`, `fold_repression`, `prob_U`, `prob_S`, `prob_W`, `prob_D`

Total rows: 4 architectures x 200 concentrations = 800 rows.

---

## 7. Phase 6: Environmental Extensions

**Directory:** `phase6_environmental/`
**Orchestrator:** `environmental_main.py`
**Output directory:** `results/phase6/`

### 7.1 Environmental Signals (`environmental_signals.py`)

**4 environments defined:**

| Environment | mce3r_multiplier | noise_scale | Biological Basis |
|-------------|-----------------|-------------|------------------|
| `baseline` | 1.0 | 1.0 | Glycerol, neutral pH (in vitro) |
| `cholesterol` | 0.3 | 1.0 | Cholesterol as sole carbon source; reduces Mce3R DNA-binding affinity 3-fold |
| `acidic_pH` | 0.7 | 1.2 | pH 5.5 (phagosomal); mildly reduces Mce3R activity + increases protein turnover |
| `host_like` | 0.2 | 1.3 | Cholesterol + acidic pH (macrophage phagosome); strong derepression + high noise |

**How multipliers are applied (`apply_environment`):**

```python
def apply_environment(model_arrays: dict, environment: str) -> dict:
    modified = deep_copy(model_arrays)
    modified['k_on'] = model_arrays['k_on'] * mce3r_multiplier
    modified['gamma_protein'] = model_arrays['gamma_protein'] * noise_scale
    return modified
```

**Critical biophysics note:** `k_off_strong` and `k_off_weak` are NOT scaled. The multiplier reduces only the forward binding rate (`k_on`), not the unbinding rate. This correctly models reduced on-rate affinity, not faster dissociation.

The function does NOT mutate the input dict -- it returns a deep copy.

### 7.2 Two-Species Model (`two_species_model.py`)

**Class:** `TwoSpeciesModel`

Models the divergent transcription unit where Mce3R (autorepressor) and Target gene are controlled by the SAME operator.

```python
class TwoSpeciesModel:
    def __init__(
        self,
        Kd_strong=None, Kd_weak=None, k_on=None,
        block_strong=None, block_weak=None, k_max=None,
        alpha_mce3r=1.0,     # transcription strength for Mce3R direction
        alpha_target=1.0,    # transcription strength for Target direction
    )
```

**Species tracked:**
1. `mce3r_mRNA` -- Mce3R mRNA count
2. `mce3r_protein` -- Mce3R protein count (the repressor itself)
3. `target_mRNA` -- target gene mRNA count
4. `target_protein` -- target gene protein count
5. `op_state` -- operator state (0-3)

**`n_bound` array:** `[0, 1, 1, 2]` -- molecules of Mce3R sequestered per operator state.

**`get_numba_arrays()` returns:**
```python
{
    'k_txn':           float64[4],
    'n_bound':         int64[4],
    'k_on':            float64,
    'k_off_strong':    float64,
    'k_off_weak':      float64,
    'k_translation':   float64,
    'gamma_mRNA':      float64,
    'gamma_protein':   float64,
    'nM_per_molecule': float64,
    'alpha_mce3r':     float64,
    'alpha_target':    float64,
}
```

### 7.3 Two-Species Gillespie (`two_species_gillespie.py`)

**18 reactions (NOT 12):**

| Reaction | Description | Propensity | Condition |
|----------|-------------|------------|-----------|
| 0 | U -> S (bind strong from empty) | `k_on * free_mce3r_nM` | `op_state == 0` |
| 1 | S -> U (unbind strong) | `k_off_strong` | `op_state == 1` |
| 2 | U -> W (bind weak from empty) | `k_on * free_mce3r_nM` | `op_state == 0` |
| 3 | W -> U (unbind weak) | `k_off_weak` | `op_state == 2` |
| 4 | W -> D (bind strong when weak occupied) | `k_on * free_mce3r_nM` | `op_state == 2` |
| 5 | D -> W (unbind strong when both bound) | `k_off_strong` | `op_state == 3` |
| 6 | S -> D (bind weak when strong occupied) | `k_on * free_mce3r_nM` | `op_state == 1` |
| 7 | D -> S (unbind weak when both bound) | `k_off_weak` | `op_state == 3` |
| 8 | Mce3R mRNA production | `k_txn[op_state] * alpha_mce3r` | always |
| 9 | Target mRNA production | `k_txn[op_state] * alpha_target` | always |
| 10 | Mce3R translation | `k_translation * mce3r_mRNA` | always |
| 11 | Target translation | `k_translation * target_mRNA` | always |
| 12 | Mce3R mRNA decay | `gamma_mRNA * mce3r_mRNA` | always |
| 13 | Target mRNA decay | `gamma_mRNA * target_mRNA` | always |
| 14 | Mce3R protein decay | `gamma_protein * mce3r_protein` | always |
| 15 | Target protein decay | `gamma_protein * target_protein` | always |
| 16 | Placeholder | 0.0 | -- |
| 17 | Placeholder | 0.0 | -- |

**FREE Mce3R computation (TRUE autoregulation):**
```python
free_mce3r = max(0, mce3r_protein - n_bound[op_state])
free_mce3r_nM = free_mce3r * nM_per_molecule
```

This creates a negative feedback loop: more Mce3R protein -> more operator binding -> more repression -> less Mce3R transcription -> less Mce3R protein.

**Initial conditions (pre-loaded to speed equilibration):**
```python
op_state       = 0
mce3r_mRNA     = 1
mce3r_protein  = 50
target_mRNA    = 1
target_protein = 50
```

**Two Numba-decorated functions:**

```python
@numba.njit
def simulate_two_species_cell(...) -> tuple[int, int, int, int, int]
# Returns: (mce3r_protein, target_protein, mce3r_mRNA, target_mRNA, op_state)

@numba.njit
def simulate_two_species_cell_with_trace(...) -> tuple[int, int, int, int, int, ndarray, ndarray, ndarray, int]
# Same but also returns trace_times, trace_mce3r, trace_target, n_trace
```

**Population runner:**
```python
def run_two_species_population(
    model: TwoSpeciesModel,
    n_cells: int,
    master_seed: int,
    environment: str = 'baseline'
) -> dict
```

Applies environmental modifications via `apply_environment()` before simulation.

### 7.4 Mutual Information (`mutual_information.py`)

**Purpose:** Compute I(Environment; Target_protein) -- how much information about the environment is encoded in expression level.

**Histogram-based estimator:**

```python
def estimate_mutual_information(
    distributions: dict,  # env_name -> protein count array
    n_bins: int = 50,
    method: str = 'histogram'  # or 'ksg'
) -> float  # MI in bits
```

**Algorithm:**

1. Pool all data to find global range
2. Build shared bin edges: `np.linspace(data_min, data_max, n_bins + 1)`
3. `P(env) = n_per_env / n_total` (prior probability of each environment)
4. `P(env, bin) = hist_matrix[i, j] / n_total` (joint probability)
5. `P(bin) = sum over envs of P(env, bin)` (marginal)
6. `MI = sum_{env, bin} P(env, bin) * log2(P(env, bin) / (P(env) * P(bin)))` with epsilon = 1e-10 to avoid log(0)
7. Clip result to >= 0

**KSG method:** Uses adaptive bins via Sturges rule: `n_bins = max(2, int(1 + log2(n_total)))`. Same histogram algorithm otherwise.

**MI vs concentration sweep:**

```python
def mutual_information_vs_concentration(
    architecture: str,     # 'asymmetric', 'symmetric', 'single_site'
    concentrations_nM: np.ndarray,
    environments: list,
    n_cells: int = 5000,
    seed: int = 42
) -> pd.DataFrame
```

For each concentration, scales `k_on` proportionally:
```python
reference_conc = 332.0  # nM
k_on_scaled = PARAMS['k_on'] * (conc_nM / reference_conc)
```

Simulates using single-species model with `apply_environment()` for each environment, then computes MI.

### 7.5 Persistence Threshold (`persistence_threshold.py`)

**Threshold definition:**

```python
def define_threshold(
    proteins_unregulated: np.ndarray,
    method: str = 'percentile',  # or 'absolute'
    percentile: float = 10.0
) -> float
```

Methods:
- `'percentile'`: `np.percentile(proteins_unregulated, 10)` from Condition D
- `'absolute'`: Fixed at 100 proteins

**Persister fraction computation:**

```python
def compute_persister_fractions(results: dict, threshold: float) -> pd.DataFrame
```

For each `(architecture, environment, model_type)` combination:
- `persister_fraction = n_below_threshold / n_cells`
- `fold_enrichment = persister_fraction / symmetric_baseline_fraction`

Reference: symmetric baseline fraction (architecture='symmetric', env='baseline', model='single').

**Output columns:** `architecture`, `environment`, `model_type`, `n_cells`, `persister_fraction`, `fold_enrichment`, `mean_protein`, `cv_protein`

**Phase diagram:**

```python
def persistence_phase_diagram(
    ratios: list,           # asymmetry ratios
    environments: list,     # environment names
    n_cells_per_point: int = 2000,
    seed: int = 42,
    threshold: float = None  # auto from Condition D if None
) -> pd.DataFrame
```

For each (ratio, environment) pair: constructs operator with `Kd_s = Kd_geo / sqrt(r)`, `Kd_w = Kd_geo * sqrt(r)`, applies environment, simulates, computes persister fraction.

**Sensitivity:** The orchestrator sweeps thresholds at 5th, 10th, and 15th percentiles to verify ranking stability.

### 7.6 Run Environmental Conditions (`run_environmental_conditions.py`)

**24 conditions:** 3 architectures x 4 environments x 2 model types

**Architectures:**
```python
ARCHITECTURES = {
    'asymmetric':  {'Kd_strong': 2.4,   'Kd_weak': 49.0},
    'symmetric':   {'Kd_strong': 10.84, 'Kd_weak': 10.84},
    'single_site': {'Kd_strong': 2.4,   'Kd_weak': 1e12},
}
```

**Environments:** `['baseline', 'cholesterol', 'acidic_pH', 'host_like']`

**Model types:** single-species (existing Gillespie engine) and two-species

**Cells per condition:** 10,000 (default)

**Seed scheme:**
- Single-species: `master_seed + arch_index*100 + env_index*10`
- Two-species: `master_seed + arch_index*100 + env_index*10 + 500`

**Trace recording:** 5 trace cells for (asymmetric, host_like, two-species) using `simulate_two_species_cell_with_trace()`. Seeds: `master_seed + 9000 + trace_cell_index`.

**Output files:**
- 24 NPZ files: `results/phase6/env_condition_{arch}_{env}_{model}.npz`
- 1 trace file: `results/phase6/two_species_results.npz`

Single-species NPZ contains: `proteins`, `mRNAs`, `op_states`

Two-species NPZ contains: `mce3r_proteins`, `target_proteins`, `mce3r_mRNAs`, `target_mRNAs`, `op_states`

### 7.7 Environmental Main (`environmental_main.py`)

**Orchestration steps (11):**

1. Load Phase 5 outputs (`mcmc_summary.json`: omega_median, dG_spacer_median, mce3r_conc_nM)
2. Load Condition D proteins for threshold definition; compute thresholds at 5th/10th/15th percentile
3. Run all 24 environmental conditions (`run_all_environmental()`)
4. Scientific check: CV(asymmetric) > CV(symmetric) under every environment
5. Scientific check: host_like mean > baseline mean for all architectures (derepression)
6. Scientific check: Pearson correlation between mce3r_protein and target_protein is negative (autoregulation)
7. Compute persister fractions for all 24 conditions; check asymmetric > symmetric under host_like
8. Compute mutual information vs concentration for 3 architectures x 10 concentrations; check MI(asymmetric) >= MI(symmetric) at physiological concentration
9. Generate persistence phase diagram (10 ratios x 4 environments x 2000 cells/point)
10. Sensitivity analysis: recompute persister fractions at 5th/10th/15th percentile thresholds
11. Write `phase6_summary.json`

**MI concentrations swept:** `[10, 30, 60, 100, 150, 200, 332, 500, 800, 1200]` nM

**Entry point:**
```python
def run_phase6(
    n_cells: int = 10000,
    master_seed: int = 12345,
    n_cells_mi: int = 5000,
    n_cells_phase_diagram: int = 2000
) -> dict
```

**Output files:**
- `persistence_fractions.csv`
- `mutual_information.csv`
- `phase_diagram.csv`
- `persistence_sensitivity.csv`
- `phase6_summary.json`

### 7.8 Phase 6 Self-Tests

**Sanity checks (3):**
1. Condition D persister fraction approximately 10% (within 2%)
2. All 24 NPZ files exist and are loadable
3. All output files exist

**Science checks (12 total, 1 typically unexpected):**

| # | Prediction | Expected Result |
|---|------------|----------------|
| 1-4 | CV(asymmetric) > CV(symmetric) under each of 4 environments | PASS (all 4) |
| 5-7 | host_like mean > baseline mean for each of 3 architectures | PASS (all 3) |
| 8-11 | Pearson(mce3r, target) < 0 for (asymmetric, baseline), (asymmetric, host_like), (symmetric, baseline), (symmetric, host_like) | PASS (all 4) |
| 12 | MI(asymmetric) >= MI(symmetric) at 332 nM | **Typically UNEXPECTED** -- MI(asymmetric) is often LESS than MI(symmetric) because higher noise obscures the environmental signal |

The MI result being unexpected is a genuine scientific finding: asymmetric operators trade information fidelity for noise-driven phenotypic bet-hedging.

---

*End of ARCHITECTURE_PART1.md*
