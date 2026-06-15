# FULL_PROJECT_ARCHITECTURE.md

## Asymmetric Operator Architecture Drives Gene Expression Noise and Antibiotic Persistence in *Mycobacterium tuberculosis*

### A Computational Study of Mce3R Regulation

---

## CLEAN BUILD REQUIREMENT

**This project MUST be built entirely from scratch.**

- Project directory: `/Users/aayanalwani/tb project/mce3r_stochastic/`
- Do NOT read, reference, copy, or reuse any existing code from `mce3r_project/` or any other directory
- Do NOT import from, read, or reference any scripts in `mce3r_project/`
- Write every function, every parameter definition, every analysis script, and every visualization from scratch
- If a concept overlaps with something that might exist elsewhere, write a new implementation — do not copy
- The ONLY shared resource is this `FULL_PROJECT_ARCHITECTURE.md` file

---

## TABLE OF CONTENTS

1. [Project Overview](#1-project-overview)
2. [Scientific Background & Hypothesis](#2-scientific-background--hypothesis)
3. [Build Strategy](#3-build-strategy)
4. [Directory Structure](#4-directory-structure)
5. [Environment & Dependencies](#5-environment--dependencies)
6. [Phase 1: Bioinformatics Pipeline](#6-phase-1-bioinformatics-pipeline)
7. [Phase 2: Gillespie Stochastic Simulation](#7-phase-2-gillespie-stochastic-simulation)
8. [Phase 3: Statistical Validation](#8-phase-3-statistical-validation)
9. [Phase 4: Publication Figures](#9-phase-4-publication-figures)
10. [Main Orchestrator (main.py)](#10-main-orchestrator-mainpy)
11. [Self-Testing Requirements](#11-self-testing-requirements)
12. [Parameter Reference Table](#12-parameter-reference-table)
13. [Literature References](#13-literature-references)
14. [Verification Checklist](#14-verification-checklist)

---

## 1. PROJECT OVERVIEW

**Research Question:** Does the physical shape (asymmetry) of the Mce3R DNA binding site control how many *M. tuberculosis* bacteria randomly become antibiotic-tolerant persister cells?

**Core Logic Chain:**
```
Asymmetric operator (Panagoda 2024)
    → 4 distinct DNA occupancy states with different transcription rates
    → Increased gene expression noise (stochastic switching)
    → Bimodal/trimodal protein distribution
    → Subpopulation crosses persister threshold
    → Predicted increase in persister frequency
    → Consistent with experimental observation (Pandey 2023)
```

**Two Computational Phases:**
- **Phase 1 — Bioinformatics:** Discover Mce3R binding sites genome-wide using MEME/FIMO
- **Phase 2 — Simulation:** Gillespie SSA modeling asymmetric operator → noise → persistence

**Deliverables:**
- Publication-quality figures (300 dpi PNG) — 7 figures total
- CSV data files for all results
- BUILD_REPORT.md auto-generated after each run
- All scripts self-testing with PASS/FAIL checks

---

## 2. SCIENTIFIC BACKGROUND & HYPOTHESIS

### 2.1 The Mce3R System

**Mce3R (Rv1963c)** is a TetR-family transcriptional repressor in *M. tuberculosis* H37Rv. It represses the mce3 operon (Rv1964-Rv1977), which encodes an ABC transporter-like lipid/cholesterol import system. Mce3R also regulates two additional transcriptional units: Rv1933c-Rv1935c and Rv1936-Rv1941 (lipid metabolism enzymes including FadE17-FadE18).

**CRITICAL: Mce3R is autoregulatory.** The mce3R gene (Rv1963c) is divergently transcribed from its own operator region, meaning Mce3R protein represses its own transcription. This negative autoregulation is a fundamental feature of the model — the repressor concentration is NOT a fixed parameter but a dynamic variable that depends on the current protein count in the Gillespie simulation.

**Unprecedented structural feature:** Mce3R is a double-TFR (TetR-fold repeat) protein — each monomer contains two fused TetR-like domains (M1: residues 1-205, M2: residues 206-406), each with its own HTH DNA-binding motif. The homodimer binds a **nonpalindromic** operator containing two 25-bp binding sites separated by 53 bp, with dramatically different affinities:

| Site | Kd | Description |
|------|----|-------------|
| Downstream (strong) | 2.4 ± 0.7 nM | High-affinity site |
| Upstream (weak) | ~49 nM (Probe C) to >100 nM (Probe B) | Low-affinity site |

**Source:** Panagoda, Balázsi & Sampson (2024) ACS Chem. Biol. 19:2580-2592. PDB: 9B7Y (2.51 Å cryo-EM).

### 2.2 The 123-bp Operator Sequence (from PDB 9B7Y)

```
5'-GCCCCGCGCTATAGGATACTAGCAAGATACATCATAGCCAATATATGCCAGTTTGCATTGCTATTTACCGATCAGTTGTCCAAGCAATCGCGTATTGGCTATGGACATCAGCGGTTCTGCCGC-3'
```

Located in the mce3R-yrbE3A intergenic region at approximately H37Rv coordinates 2,207,477–2,207,699.

### 2.3 Connection to Persistence

Pandey et al. (2023, Res. Microbiol. 174:104082) showed that Δmce3R increases antibiotic persister frequency. Flentie et al. (2019, ACS Infect. Dis.) found 6-azasteroid compounds targeting the Mce3R pathway enhanced isoniazid activity 16-fold and bedaquiline ~50-fold.

### 2.4 Noise–Persistence Theory

Balázsi's work (Cell 2011, Nat. Commun. 2019) established that:
- **Positive feedback** amplifies noise → bistability → drug resistance
- **Negative feedback** suppresses noise → narrows distributions
- In CHO cells, high-noise circuits facilitated survival at high drug concentrations
- In mycobacteria, rel expression shows bimodal distribution driving persistence (Sureka/Balázsi 2008 PLoS ONE)

**Our hypothesis:** The asymmetric Mce3R operator creates an intermediate regulatory architecture — not fully ON or OFF — that generates excess noise compared to a symmetric operator, increasing the fraction of cells that stochastically cross a persistence threshold.

### 2.5 Model Interpretation: What "Protein" Represents

**IMPORTANT CLARIFICATION:** In the minimal Gillespie model (Phase 2), the simulated "protein" species is an **abstract persistence-linked downstream output** controlled by the Mce3R operator logic. It is NOT literal intracellular Mce3R abundance.

The autoregulatory negative feedback loop (Mce3R represses its own operator) motivates the dynamic, protein-count-dependent binding propensities in the simulation. However, the protein distribution used for threshold analysis, fold-change comparisons, and multimodality testing should be interpreted as downstream functional expression of the mce3 regulon (e.g., cholesterol import machinery, lipid metabolism enzymes) — the output that is biologically linked to persistence.

In this minimal abstraction, the single simulated protein plays a dual modeling role: it acts as the effective regulatory species controlling operator occupancy and also serves as the proxy readout for downstream persistence-linked output. This is a reduced phenomenological model of regulatory architecture, not a literal two-species biochemical representation. The key prediction — that asymmetric operators produce more noise than symmetric ones — is architecture-dependent and holds regardless of whether the readout is the repressor itself or a downstream target.

---

## 3. BUILD STRATEGY

### 3.1 Sequential Code Generation, Parallel Runtime Execution

**BUILD ORDER (Sequential file generation):**
You MUST write all code sequentially in this order:
1. `config/parameters.py`
2. All `phase1_pipeline/` files
3. All `phase2_simulation/` files
4. All `phase3_analysis/` files
5. All `phase4_figures/` files
6. `main.py` (orchestrator)

Do NOT attempt to spawn background agents, parallel build processes, or concurrent file generation. Write one file at a time.

**RUNTIME PARALLELISM (in main.py):**
When the user runs `python main.py`, Phase 1 and Phase 2 have no dependencies on each other and MUST execute simultaneously using Python's `subprocess` module:

```python
import subprocess
import sys

# Launch Phase 1 and Phase 2 in parallel
p1 = subprocess.Popen([sys.executable, '-m', 'phase1_pipeline.pipeline_main'],
                      stdout=open('logs/phase1_pipeline.log', 'w'),
                      stderr=subprocess.STDOUT)
p2 = subprocess.Popen([sys.executable, '-m', 'phase2_simulation.simulation_main'],
                      stdout=open('logs/phase2_gillespie.log', 'w'),
                      stderr=subprocess.STDOUT)

# Wait for both to finish
p1.wait()
p2.wait()

# Check return codes before proceeding to Phase 3
if p1.returncode != 0 or p2.returncode != 0:
    write_build_report(failed=True, ...)
    sys.exit(1)
```

Phase 3 runs only after Phase 2 succeeds. Phase 4 runs only after Phase 3 succeeds. These are sequential.

### 3.2.1 Phase Orchestrator Return Contract

Each phase orchestrator function (`run_phase1()`, `run_phase2()`, `run_phase3()`, `run_phase4()`) MUST return a structured dictionary with at least these keys:

```python
{
    'status': 'COMPLETE' | 'FAILED',         # overall phase status
    'n_sanity_pass': int,                     # software sanity checks passed
    'n_sanity_fail': int,                     # software sanity checks failed
    'n_science_expected': int,                # scientific expectations met
    'n_science_unexpected': int,              # scientific expectations not met
    'n_warn': int,                            # warnings (e.g., mock output used)
    'mock_used': bool,                        # True if MEME/FIMO mock was used
    'wall_time_sec': float,                   # elapsed time for this phase
    'output_files': list[str],                # paths to all generated output files
}
```

`main.py` uses these returned dictionaries to generate BUILD_REPORT.md. For Phase 1 and Phase 2, which run as subprocesses, each phase must also write a small machine-readable summary file (e.g., `results/phase1/phase1_summary.json` and `results/phase2/phase2_summary.json`) containing the same structured fields. `main.py` reads these summary files after subprocess completion to populate BUILD_REPORT.md.

### 3.2 Dependency Graph

```
INDEPENDENT (run in parallel via subprocess):
  ├── Phase 1: bioinformatics pipeline
  └── Phase 2: Gillespie simulation engine

SEQUENTIAL (must wait):
  └── Phase 3: statistical validation (depends on Phase 2)
      └── Phase 4: publication figures (depends on Phase 3)
```

### 3.3 Agent Communication & Logging

Each phase writes to its own log file:
```
logs/phase1_pipeline.log
logs/phase2_gillespie.log
logs/phase3_analysis.log
logs/phase4_figures.log
```

**Log format:**
```
[TIMESTAMP] [PHASE] [LEVEL] Message
[2024-01-15 14:23:01] [Phase 2] [INFO] Starting Gillespie SSA for 5 conditions...
[2024-01-15 14:23:45] [Phase 2] [PASS] Fano factor = 1.12 (expected 0.8-1.5)
[2024-01-15 14:24:00] [Phase 2] [FAIL] Mean protein = 52 (expected 100-500)
```

**Console prints (high-level only):**
```
[Phase 2] Starting Gillespie SSA for 5 conditions...
[Phase 2] Condition A (asymmetric): 10000/10000 cells complete
[Phase 2] PASS: Fano factor = 1.12 (expected 0.8-1.5)
```

**BUILD_REPORT.md** (auto-generated by orchestrator after completion):
```markdown
# Build Report
- Timestamp: YYYY-MM-DD HH:MM:SS
- Status: COMPLETE | FAILED_AT_PHASE_X
- Phase 1: X/Y checks passed, wall time: Xs [REAL | MOCK]
- Phase 2: X/Y checks passed, wall time: Xs
- Phase 3: X/Y checks passed, wall time: Xs
- Phase 4: X/Y checks passed, wall time: Xs
- Output files: [list of all generated files with sizes]
- Next steps: [if failed, what to fix]
```

---

## 4. DIRECTORY STRUCTURE

```
mce3r_stochastic/
├── main.py                          # Orchestrator: runs everything
├── environment.yml                  # Conda environment specification
├── FULL_PROJECT_ARCHITECTURE.md     # This file
├── BUILD_REPORT.md                  # Auto-generated after each run
├── working.md                       # Current status / handoff notes
│
├── config/
│   ├── __init__.py
│   └── parameters.py                # ALL tunable parameters in one place
│
├── data/
│   ├── genomes/                     # Downloaded genome FASTA files
│   │   ├── H37Rv.fasta              # NC_000962.3
│   │   ├── H37Rv.gff                # Gene annotations
│   │   ├── M_bovis.fasta            # NC_002945.4
│   │   └── M_marinum.fasta          # NC_010612.1
│   ├── sequences/                   # Extracted sequences
│   │   ├── upstream_regions.fasta   # 200bp upstream of all H37Rv genes
│   │   ├── meme_input.fasta         # Orthologous yrbE3A upstream from 3 species
│   │   └── known_operator.fasta     # The 123bp operator from PDB 9B7Y
│   └── external/                    # Published data for validation
│       └── santangelo2009_foldchange.csv  # Fold-changes from Santangelo 2009
│
├── phase1_pipeline/
│   ├── __init__.py
│   ├── download_genomes.py          # Fetch genomes from NCBI via Bio.Entrez
│   ├── extract_upstream.py          # Extract promoter regions (circular genome aware)
│   ├── run_meme.py                  # MEME motif discovery wrapper (with mock fallback)
│   ├── run_fimo.py                  # FIMO genome-wide scan wrapper (with mock fallback)
│   ├── conservation_check.py        # Cross-species local sequence conservation validation
│   └── pipeline_main.py            # Phase 1 orchestrator
│
├── phase2_simulation/
│   ├── __init__.py
│   ├── gillespie_engine.py          # Core SSA implementation (Numba JIT compiled)
│   ├── operator_model.py            # 4-state operator with dynamic protein-dependent binding
│   ├── run_conditions.py            # Run all 5 conditions (A, B, C, D, E)
│   ├── asymmetry_sweep.py           # Sweep asymmetry ratio 1-50 (Condition E)
│   └── simulation_main.py          # Phase 2 orchestrator
│
├── phase3_analysis/
│   ├── __init__.py
│   ├── noise_metrics.py             # CV, Fano factor, bimodality coefficient, GMM
│   ├── bootstrap_ci.py              # Bootstrap confidence intervals
│   ├── sensitivity_analysis.py      # Parameter sensitivity (tornado plot data)
│   ├── statistical_tests.py         # KS test, likelihood ratio test, GMM
│   ├── experimental_comparison.py   # Compare to Pandey 2023, Santangelo 2009
│   ├── negative_controls.py         # Shuffle test, symmetric control, Poisson check
│   └── analysis_main.py            # Phase 3 orchestrator
│
├── phase4_figures/
│   ├── __init__.py
│   ├── figure_style.py              # Shared matplotlib style settings
│   ├── fig1_binding_sites.py        # Genome-wide binding site map
│   ├── fig2_distributions.py        # Protein distributions (all conditions)
│   ├── fig3_asymmetry_sweep.py      # Sweep: asymmetry ratio vs intermediate fraction
│   ├── fig4_sensitivity.py          # Tornado plot + heatmap
│   ├── fig5_validation.py           # Experimental comparison + negative controls
│   ├── fig6_single_cell_traces.py   # Time-series traces for representative cells
│   ├── fig7_bic_comparison.py       # BIC comparison across conditions
│   └── figures_main.py             # Phase 4 orchestrator
│
├── logs/
│   ├── phase1_pipeline.log
│   ├── phase2_gillespie.log
│   ├── phase3_analysis.log
│   └── phase4_figures.log
│
├── results/
│   ├── phase1/
│   │   ├── meme_output/             # Raw MEME output directory
│   │   ├── fimo_output/             # Raw FIMO output directory
│   │   ├── predicted_sites.csv      # Ranked binding site predictions
│   │   └── conservation_status.csv  # Cross-species conservation
│   ├── phase2/
│   │   ├── condition_A_asymmetric.npz    # Raw simulation data
│   │   ├── condition_B_symmetric.npz
│   │   ├── condition_C_single_site.npz
│   │   ├── condition_D_unregulated.npz
│   │   ├── condition_E_sweep.npz
│   │   ├── traces/                       # Subsampled trajectories for fig6
│   │   │   ├── traces_A.npz
│   │   │   ├── traces_B.npz
│   │   │   └── traces_C.npz
│   │   └── steady_state_distributions.csv
│   ├── phase3/
│   │   ├── noise_metrics.csv
│   │   ├── bootstrap_results.csv
│   │   ├── sensitivity_data.csv
│   │   ├── statistical_tests.csv
│   │   └── negative_controls.csv
│   └── figures/
│       ├── fig1_binding_sites.png
│       ├── fig2_distributions.png
│       ├── fig3_asymmetry_sweep.png
│       ├── fig4_sensitivity.png
│       ├── fig5_validation.png
│       ├── fig6_single_cell_traces.png
│       └── fig7_bic_comparison.png
│
└── tests/                           # Integration tests (optional)
    └── test_integration.py
```

---

## 5. ENVIRONMENT & DEPENDENCIES

### 5.1 Conda Environment (environment.yml)

```yaml
name: mce3r_project
channels:
  - conda-forge
  - bioconda
  - defaults
dependencies:
  - python=3.11
  - numpy>=1.24
  - scipy>=1.11
  - pandas>=2.0
  - matplotlib>=3.7
  - seaborn>=0.12
  - biopython>=1.81
  - scikit-learn>=1.3      # For GMM fitting
  - statsmodels>=0.14      # For statistical tests
  - numba>=0.58            # REQUIRED: JIT compilation for Gillespie loop
  - tqdm>=4.65             # Progress bars
```

**MEME Suite** is assumed pre-installed and available on PATH (`meme`, `fimo` commands). If not available, scripts fall back to mock output (see Section 6.3).

### 5.2 Setup Command

```bash
conda env create -f environment.yml
conda activate mce3r_project
python main.py
```

---

## 6. PHASE 1: BIOINFORMATICS PIPELINE

### 6.1 download_genomes.py

**Purpose:** Download reference genomes from NCBI.

**Inputs:** None (accession numbers hardcoded in parameters.py)

**Outputs:**
- `data/genomes/H37Rv.fasta` — M. tuberculosis H37Rv (NC_000962.3)
- `data/genomes/H37Rv.gff` — Gene annotations
- `data/genomes/M_bovis.fasta` — M. bovis AF2122/97 (NC_002945.4)
- `data/genomes/M_marinum.fasta` — M. marinum M (NC_010612.1)

**Algorithm:**

Do NOT use raw `requests` to download from NCBI. Use `Bio.Entrez` which handles rate limiting and retries properly:

```python
from Bio import Entrez, SeqIO

Entrez.email = "mce3r_project@example.com"  # Required by NCBI policy

# Download genome
handle = Entrez.efetch(db="nucleotide", id=accession, rettype="fasta", retmode="text")
record = SeqIO.read(handle, "fasta")
handle.close()

# Download GFF
handle = Entrez.efetch(db="nucleotide", id=accession, rettype="gff3", retmode="text")
```

1. Download all 3 genomes via Bio.Entrez
2. Validate: check file sizes > 4 MB for H37Rv, confirm sequence starts with expected characters
3. Calculate and log GC content — must be ~65.6% for H37Rv

**Self-checks:**
- H37Rv FASTA contains exactly 1 sequence of ~4,411,532 bp
- GC content is between 65.0% and 66.0%
- All 3 genome files exist and are non-empty

### 6.2 extract_upstream.py

**Purpose:** Extract 200 bp upstream of annotated genes for motif scanning.

**Inputs:**
- `data/genomes/H37Rv.fasta`
- `data/genomes/H37Rv.gff`
- `data/genomes/M_bovis.fasta`, `M_marinum.fasta` (for ortholog extraction only)

**Outputs:**
- `data/sequences/upstream_regions.fasta` — 200 bp upstream of ALL H37Rv annotated gene starts (~4,000 sequences)
- `data/sequences/meme_input.fasta` — 200 bp upstream of yrbE3A (Rv1964) ortholog from all 3 species (3 sequences)
- `data/sequences/known_operator.fasta` — The 123 bp operator sequence from PDB 9B7Y

**Algorithm:**

**GFF COORDINATE INDEXING (CRITICAL):**
GFF coordinates are 1-based inclusive. Python string slicing is 0-based end-exclusive. You MUST convert before slicing:
```python
start_0 = start_1_based - 1       # Convert to 0-based
end_0_exclusive = end_1_based      # Already correct for Python slicing
```
Never mix raw GFF coordinates with Python slice indices.

1. Parse GFF to get gene start/end coordinates and strand for all CDS features. Convert to 0-based immediately upon parsing.
2. For each gene, extract 200 bp upstream (accounting for strand):
   - Plus strand: upstream region is BEFORE gene start
   - Minus strand: upstream region is AFTER gene end (reverse complement)
3. **CIRCULAR GENOME HANDLING (BOTH ENDS):** M. tuberculosis has a circular chromosome. Handle wrapping at BOTH the origin AND terminus:
   ```python
   genome_length = len(genome_seq)
   L = upstream_length  # 200

   # Plus strand: upstream is before start_0
   if strand == '+':
       if start_0 - L < 0:
           # Wraps around origin (gene near coordinate 0)
           upstream = genome_seq[genome_length + (start_0 - L):] + genome_seq[:start_0]
       else:
           upstream = genome_seq[start_0 - L : start_0]

   # Minus strand: upstream is after end_0_exclusive (reverse complement)
   elif strand == '-':
       if end_0_exclusive + L > genome_length:
           # Wraps around terminus (gene near end of genome)
           upstream = genome_seq[end_0_exclusive:] + genome_seq[:L - (genome_length - end_0_exclusive)]
       else:
           upstream = genome_seq[end_0_exclusive : end_0_exclusive + L]
       upstream = reverse_complement(upstream)
   ```
4. Handle edge cases: overlapping genes (truncate if upstream overlaps another CDS)
5. For MEME input: extract yrbE3A ortholog upstream from M. bovis and M. marinum using known locus tags or manually specified ortholog identifiers defined in parameters.py (do not use remote BLAST)
6. Write the known 123-bp operator as a separate FASTA

**Self-checks:**
- upstream_regions.fasta contains between 3,500 and 4,500 sequences
- All sequences are exactly 200 bp (no sequence shorter than 200bp due to missed wrapping)
- The upstream region of Rv1964 (yrbE3A) contains the known operator sequence (substring match or local sequence alignment)
- meme_input.fasta contains exactly 3 sequences
- No sequences contain characters other than A, C, G, T, N
- **Wrapping test (plus strand):** Verify correct extraction for a synthetic plus-strand gene at coordinate 50 (should wrap around origin to grab from genome end)
- **Wrapping test (minus strand):** Verify correct extraction for a synthetic minus-strand gene ending at genome_length - 50 (should wrap around terminus)
- **Length assertion:** Assert every extracted sequence has length == upstream_length (200 bp exactly)

### 6.3 run_meme.py

**Purpose:** De novo motif discovery from orthologous upstream regions.

**Inputs:** `data/sequences/meme_input.fasta` (3 sequences from 3 species)

**Outputs:** `results/phase1/meme_output/` (standard MEME output directory containing meme.txt, meme.html, etc.)

**MEME SUITE AVAILABILITY CHECK:**
At the start, check if the required command exists:

```python
import shutil

def check_meme_installed():
    if shutil.which('meme') is None:
        print("WARNING: 'meme' command not found on PATH.")
        print("Generating mock output for downstream testing.")
        generate_mock_meme_output()
        return False
    return True
```

**Mock output requirements:**
- `generate_mock_meme_output()` creates `results/phase1/meme_output/meme.txt` with a single valid motif entry — a 25-bp PWM derived from the downstream (strong) binding site within the 123-bp operator. The mock motif must be 20-30 bp wide (matching MEME's configured search range), NOT the full 123-bp operator.
- Mock outputs MUST be clearly labeled in the log: `[Phase 1] [WARNING] Using MOCK output — MEME Suite not installed`
- The BUILD_REPORT.md must flag `Phase 1: MOCK (MEME Suite not available)` so it's never mistaken for real results

**Command constructed and executed (if MEME is available):**
```bash
meme data/sequences/meme_input.fasta \
  -dna \
  -mod zoops \
  -revcomp \
  -minw 20 \
  -maxw 30 \
  -nmotifs 3 \
  -oc results/phase1/meme_output/ \
  -bfile <generated_background_file>
```

**CRITICAL: Do NOT use -pal flag.** Mce3R binds a nonpalindromic operator. Using -pal would force palindrome symmetry and miss the asymmetry.

**Background model:**
- Generate a 2nd-order Markov background model from the H37Rv genome using `fasta-get-markov` (MEME Suite tool):
  ```bash
  fasta-get-markov -m 2 data/genomes/H37Rv.fasta data/genomes/H37Rv.bg
  ```
- This accounts for the 65.6% GC content of H37Rv

**Self-checks:**
- If MEME Suite is missing: print WARNING (not FAIL) and test only downstream parsing logic against mock output
- If MEME Suite IS present:
  - MEME output directory exists and contains meme.txt
  - At least 1 motif found (check meme.txt for "MOTIF" lines)
  - Top motif E-value < 0.05
  - Top motif width is between 20 and 30 bp
  - Top motif aligns to one of the two ~25-bp binding sites within the known 123-bp operator region (not the full 123-bp sequence — MEME finds single-site motifs of 20-30 bp)

### 6.4 run_fimo.py

**Purpose:** Scan the entire H37Rv genome for matches to the discovered motif.

**Inputs:**
- `results/phase1/meme_output/meme.txt` (motif definition)
- `data/genomes/H37Rv.fasta` (target genome)
- `data/genomes/H37Rv.bg` (background model)

**Outputs:**
- `results/phase1/fimo_output/` (standard FIMO output directory)
- `results/phase1/predicted_sites.csv` (parsed and ranked predictions)

**FIMO AVAILABILITY CHECK:**
Same pattern as run_meme.py — check `shutil.which('fimo')`. If missing, `generate_mock_fimo_output()` creates `results/phase1/fimo_output/fimo.tsv` with 10 plausible 25-bp hits at realistic H37Rv coordinates, including hits near the known mce3 operator region and known regulon genes. Log WARNING.

**Command (if FIMO available):**
```bash
fimo --thresh 1e-4 \
     --bfile data/genomes/H37Rv.bg \
     --oc results/phase1/fimo_output/ \
     results/phase1/meme_output/meme.txt \
     data/genomes/H37Rv.fasta
```

**Post-processing (Python):**
1. Parse `fimo_output/fimo.tsv`
2. For each hit, determine:
   - Nearest downstream gene (from GFF annotations)
   - Distance to gene start
   - Whether it falls in an intergenic region
   - Strand relative to the nearest gene
3. Filter: keep only hits in intergenic regions within 300 bp upstream of a gene start
4. Rank by FIMO p-value (ascending)
5. Add columns: gene_name, gene_product, distance_to_start, strand_match
6. Flag known Mce3R-regulated genes (Rv1964-Rv1977, Rv1933c-Rv1935c, Rv1936-Rv1941)
7. Write to `predicted_sites.csv`

**Self-checks:**
- If FIMO missing: WARNING, test parsing logic against mock output only
- If FIMO present:
  - FIMO output contains at least 5 hits at p < 1e-4
  - The known operator site (upstream of Rv1964/yrbE3A) is in the top 5 hits
  - Known Mce3R-regulated genes (Rv1933c-Rv1935c, Rv1936-Rv1941) are recovered
  - Total hits at p < 1e-4 is between 5 and 500 (sanity range)
  - No hits have p-value exactly 0 (would indicate error)

### 6.5 conservation_check.py

**Purpose:** Validate top FIMO predictions by checking cross-species conservation.

**Inputs:**
- `results/phase1/predicted_sites.csv`
- `data/genomes/M_bovis.fasta`
- `data/genomes/M_marinum.fasta`

**Outputs:**
- `results/phase1/conservation_status.csv` — predicted_sites.csv with added columns: conserved_M_bovis (bool), conserved_M_marinum (bool), conservation_score (0-2)

**Algorithm:**

**MUST use local sequence search only. Do NOT use Bio.Blast.NCBIWWW or any remote BLAST service.**

1. For each of the top 50 predicted sites (by p-value):
   - Extract the binding site sequence (~25 bp) from H37Rv
   - Perform a **local sliding-window percent-identity search** against the full target genome sequence:
     - Slide the motif across the target genome one base at a time
     - At each position, compute percent identity (matches / motif_length × 100)
     - Record the best-matching position and its identity score
   - A site is **"conserved"** if the best local match has **≥80% identity over ≥80% of motif length**
2. Conservation score = number of species (0, 1, or 2) where the site is conserved
3. Merge with predicted_sites.csv

**Self-checks:**
- The known mce3 operator site has conservation_score = 2 (conserved in both species)
- At least 3 of the top 10 predictions are conserved in at least one species
- Output CSV has same number of rows as input (top 50)

### 6.6 pipeline_main.py

**Purpose:** Phase 1 orchestrator — runs all pipeline steps in sequence.

```python
def run_phase1():
    """Execute the complete bioinformatics pipeline."""
    # Step 1: Download genomes
    # Step 2: Extract upstream regions
    # Step 3: Generate background model
    # Step 4: Run MEME (or mock)
    # Step 5: Run FIMO (or mock)
    # Step 6: Conservation check
    # Step 7: Log summary
```

**Logging:** All output to `logs/phase1_pipeline.log`

---

## 7. PHASE 2: GILLESPIE STOCHASTIC SIMULATION

### 7.1 operator_model.py

**Purpose:** Define the 4-state operator model and transition rates. **CRITICAL: Binding rates depend on current free protein count, not a fixed concentration.**

**The 4 Operator States:**

| State | Description | Transcription Rate |
|-------|-------------|-------------------|
| 0 | Both sites empty | k_txn[0] = k_max = 0.15 mRNA/min |
| 1 | Strong site occupied | k_txn[1] = k_max × (1 - block_strong) = 0.15 × 0.15 = 0.0225 mRNA/min |
| 2 | Weak site occupied | k_txn[2] = k_max × (1 - block_weak) = 0.15 × 0.50 = 0.075 mRNA/min |
| 3 | Both sites occupied | k_txn[3] = k_max × (1 - block_strong) × (1 - block_weak) = 0.15 × 0.15 × 0.50 = 0.01125 mRNA/min |

**Binding/Unbinding Kinetics:**

The operator transitions between states based on repressor binding and unbinding:

```
State 0 ⇌ State 1    (strong site: k_on × [R_free], k_off_s)
State 0 ⇌ State 2    (weak site: k_on × [R_free], k_off_w)
State 1 ⇌ State 3    (weak site binds when strong already bound)
State 2 ⇌ State 3    (strong site binds when weak already bound)
```

**Where [R_free] is the current free protein concentration in nM, computed dynamically at each Gillespie step.**

**Dissociation constants:**
- Strong site: Kd_strong = 2.4 nM (Panagoda 2024)
- Weak site: Kd_weak = 49 nM (Panagoda 2024, Probe C measurement)

**Converting Kd to k_on/k_off:**
- k_on = 0.0167 nM⁻¹ min⁻¹ (Stormo & Zhao 2010; 10⁶ M⁻¹ s⁻¹ converted: 10⁶ M⁻¹ s⁻¹ × 60 s/min × 10⁻⁹ M/nM = 0.0167 nM⁻¹ min⁻¹)
- k_off_strong = Kd_strong × k_on = 2.4 × 0.0167 = 0.0401 min⁻¹
- k_off_weak = Kd_weak × k_on = 49.0 × 0.0167 = 0.818 min⁻¹

**Dynamic repressor concentration [R_free]:**
- NOT a fixed parameter. Computed at each Gillespie step from the current protein count.
- Convert molecule count to nM: `R_nM = free_protein_count × nM_per_molecule`
- `nM_per_molecule = 1.66` (1 molecule in 1 fL ≈ 1.66 nM, from 1/(Avogadro × 1e-15 L) × 1e9)
- `cell_volume_fL = 1.0` (standard Mtb cell volume in femtoliters)

**Free protein calculation:**
In this minimal model with a SINGLE operator copy per cell, the number of bound Mce3R molecules depends on operator state:
- State 0: 0 bound
- State 1: 1 bound (on strong site)
- State 2: 1 bound (on weak site)
- State 3: 2 bound (on both sites)
- `free_protein = max(0, protein - n_bound)`

**Implementation:**

```python
class OperatorModel:
    """4-state operator model for Mce3R binding."""

    def __init__(self, params):
        self.Kd_strong = params['Kd_strong']  # nM
        self.Kd_weak = params['Kd_weak']      # nM
        self.k_on = params['k_on']            # nM^-1 min^-1
        self.k_max = params['k_max']           # mRNA/min
        self.block_strong = params['block_strong']
        self.block_weak = params['block_weak']
        self.nM_per_molecule = params['nM_per_molecule']

    def get_transcription_rates(self):
        """Return length-4 array of transcription rates indexed by state."""
        rates = np.zeros(4)
        rates[0] = self.k_max
        rates[1] = self.k_max * (1 - self.block_strong)
        rates[2] = self.k_max * (1 - self.block_weak)
        rates[3] = self.k_max * (1 - self.block_strong) * (1 - self.block_weak)
        return rates

    def get_binding_rate_constants(self):
        """Return rate constants (NOT propensities).
        Gillespie engine multiplies by current protein concentration."""
        k_off_strong = self.Kd_strong * self.k_on   # min^-1
        k_off_weak = self.Kd_weak * self.k_on       # min^-1
        return {
            'k_on': self.k_on,              # nM^-1 min^-1
            'k_off_strong': k_off_strong,   # min^-1
            'k_off_weak': k_off_weak,       # min^-1
        }
```

**Self-checks:**
- Transcription rates array is monotonically: rates[0] > rates[2] > rates[1] > rates[3]
- All rates are positive
- Detailed balance: k_on × k_on (0→1→3) path = k_on × k_on (0→2→3) path ✓
- k_off_strong = 0.0401 min⁻¹ (verify calculation)
- k_off_weak = 0.818 min⁻¹ (verify calculation)

### 7.2 gillespie_engine.py

**Purpose:** Core Gillespie Stochastic Simulation Algorithm (SSA) implementation.

**PERFORMANCE REQUIREMENT:**
The `simulate_cell()` function and all functions it calls (calculate_propensities, execute_reaction) MUST be decorated with `@numba.njit` to compile to machine code.

**Constraints for Numba compatibility:**
- No Python objects (dicts, classes) inside the JIT loop. Use only NumPy arrays and scalar primitives (int64, float64).
- Store operator transition rates as a pre-computed NumPy array, not a dict.
- Store transcription rates as a length-4 float64 array indexed by operator state.
- The OperatorModel class is used OUTSIDE the loop to pre-compute arrays; those arrays are passed INTO the Numba-compiled function.
- Use `np.random.random()` (Numba-compatible) not `np.random.Generator`.

**Species tracked per cell:**
- `operator_state` — integer 0-3 (which sites are occupied)
- `mRNA` — integer count of mRNA molecules
- `protein` — integer count of protein molecules
- `free_protein` — `max(0, protein - n_bound)` where n_bound depends on operator_state:
  - State 0: 0 bound; State 1: 1 bound; State 2: 1 bound; State 3: 2 bound

**Reactions:**

| # | Reaction | Propensity | Description |
|---|----------|-----------|-------------|
| 1-8 | State i → State j | See below (protein-dependent) | Operator binding/unbinding |
| 9 | ∅ → mRNA | k_txn[operator_state] | Transcription (rate depends on operator state) |
| 10 | mRNA → mRNA + Protein | k_translation × mRNA | Translation |
| 11 | mRNA → ∅ | gamma_mRNA × mRNA | mRNA degradation |
| 12 | Protein → ∅ | gamma_protein × protein | Protein degradation + dilution |

**Binding propensity calculation (DYNAMIC, protein-dependent):**
```python
# Convert free protein count to concentration in nM
R_nM = free_protein_count * nM_per_molecule

# Binding propensities (concentration-dependent)
# bind strong: rate = k_on * R_nM (if strong site is unoccupied)
# bind weak: rate = k_on * R_nM (if weak site is unoccupied)
# Unbinding propensities are concentration-independent (just k_off)
```

**Reconciled Parameters (from literature):**

| Parameter | Symbol | Value | Source |
|-----------|--------|-------|--------|
| Max transcription rate | k_max | 0.15 mRNA/min | Estimated from Mtb transcriptomics |
| Translation rate | k_translation | 0.5 protein/(mRNA·min) | Taniguchi 2010 adjusted for Mtb |
| mRNA half-life | t_half_mRNA | 9.5 min | Rustad et al. 2013 NAR |
| mRNA degradation rate | gamma_mRNA | ln(2)/9.5 = 0.0730 min⁻¹ | Calculated |
| Protein half-life (effective) | t_half_protein | 1500 min (~25 hr) | Mtb doubling ~20h + degradation |
| Protein degradation rate | gamma_protein | ln(2)/1500 = 0.000462 min⁻¹ | Calculated |
| Binding on-rate | k_on | 0.0167 nM⁻¹ min⁻¹ | Stormo & Zhao 2010 (10⁶ M⁻¹ s⁻¹) |
| Strong site Kd | Kd_strong | 2.4 nM | Panagoda 2024 |
| Weak site Kd | Kd_weak | 49.0 nM | Panagoda 2024 (Probe C) |
| Strong site k_off | k_off_strong | 0.0401 min⁻¹ | Kd_strong × k_on |
| Weak site k_off | k_off_weak | 0.818 min⁻¹ | Kd_weak × k_on |
| Strong site blocking | block_strong | 0.85 | Model parameter |
| Weak site blocking | block_weak | 0.50 | Model parameter |
| Cell volume | cell_volume_fL | 1.0 fL | Standard Mtb |
| nM per molecule | nM_per_molecule | 1.66 nM | 1/(Avogadro × 1e-15 L) × 1e9 |

**Algorithm (standard Gillespie direct method, Numba-compiled):**

```python
from numba import njit

@njit
def simulate_cell(k_txn, k_translation, gamma_mRNA, gamma_protein,
                  k_on, k_off_strong, k_off_weak, nM_per_molecule,
                  t_max, t_burn_in, save_trace, trace_interval):
    """Core Gillespie loop — compiled to machine code by Numba.

    Args:
        k_txn: length-4 float64 array of transcription rates per operator state
        save_trace: bool — if True, record trajectory at trace_interval
        trace_interval: float — minutes between recorded trace points

    NUMBA NOTE: Trace arrays must be pre-allocated to max size before the loop.
    Max trace points = (t_max - t_burn_in) / trace_interval = (4200 - 2100) / 10 = 210.
    Pre-allocate: trace_times = np.empty(250, dtype=np.float64)   # 250 for safety margin
                  trace_proteins = np.empty(250, dtype=np.int64)
    Use a counter to track actual points recorded, then truncate after loop returns.
    You CANNOT dynamically grow arrays inside @njit.

    Returns:
        final_protein: int — protein level at end of simulation
        final_mRNA: int — mRNA level at end of simulation
        final_op_state: int — operator state at end
        trace_times: float64 array (empty if save_trace=False)
        trace_proteins: int64 array (empty if save_trace=False)
    """
    # Initialize
    operator_state = 0  # start with both sites empty (consistent with protein=0)
    mRNA = 0
    protein = 0
    t = 0.0

    # n_bound lookup: state -> number of bound Mce3R molecules
    n_bound_lookup = np.array([0, 1, 1, 2], dtype=np.int64)

    while t < t_max:
        # Compute free protein
        n_bound = n_bound_lookup[operator_state]
        free_protein = max(0, protein - n_bound)
        R_nM = free_protein * nM_per_molecule

        # 1. Calculate all propensities
        # ... (pure numeric, no Python objects)
        # Binding: k_on * R_nM (for unoccupied sites)
        # Unbinding: k_off_strong or k_off_weak (for occupied sites)

        a_total = sum(propensities)
        if a_total == 0:
            break

        # 2. Draw waiting time
        tau = -np.log(np.random.random()) / a_total

        # 3. Select and execute reaction
        # ... standard Gillespie selection

        # 4. Advance time
        t += tau

        # 5. Record trace if saving
        if save_trace and t >= t_burn_in:
            # Record at trace_interval spacing
            ...

    return final_protein, final_mRNA, final_op_state, trace_times, trace_proteins
```

**Simulation parameters:**
- Mtb doubling time ≈ 20 hours = 1,200 min
- `t_max` = 4,200 min (70 hours ≈ 3.5 doubling times) — ensures steady-state is reached
- `t_burn_in` = 2,100 min (first half discarded as transient)
- Record protein level at the END of each trajectory (steady-state snapshot)
- For distribution: collect final protein value from N independent cells
- For traces (first `n_trace_cells` cells only): record protein every `trace_interval` minutes after burn-in
- Max trace points per cell = (4,200 - 2,100) / 10 = 210. Pre-allocate arrays to 250 for safety margin.

**Self-checks:**
- Single cell trajectory: mRNA count stays in range [0, 50] (sanity)
- Single cell trajectory: protein count stays in range [0, 5000] (sanity)
- Operator state transitions occur (not stuck in one state)
- With no regulation (all states have k_max): Fano factor of protein ≈ 1 + k_translation/gamma_mRNA ≈ 1 + 0.5/0.073 ≈ 7.85 (two-stage burst model)
- Mean mRNA at steady state (no regulation) ≈ k_max/gamma_mRNA = 0.15/0.073 ≈ 2.05

### 7.3 run_conditions.py

**Purpose:** Run all 5 simulation conditions.

**Condition A — Asymmetric Operator (the real Mce3R system):**
- Kd_strong = 2.4 nM, Kd_weak = 49 nM
- block_strong = 0.85, block_weak = 0.50
- N = 10,000 cells (first 20 save traces)
- Scientific expectation: broad or multimodal protein distribution

**Condition B — Symmetric Operator (strong-strong control):**
- Kd_strong = 2.4 nM, Kd_weak = 2.4 nM (both sites use strong affinity)
- block_strong = 0.85, block_weak = 0.85 (both block equally)
- N = 10,000 cells (first 20 save traces)
- Scientific expectation: tighter, unimodal distribution (less noise)

**Condition C — Single Site (minimal control):**
- Simulates a single-binding-site operator (strong site only, Kd = 2.4 nM, block = 0.85)
- **MUST use the same 4-state Numba engine as all other conditions.** Do NOT build a separate 2-state engine.
- Disable the weak site numerically: set Kd_weak = 1e12 nM (so k_off_weak = 1e12 × k_on = enormous, weak site is effectively never occupied)
- The 4-state engine still has states 0-3, but states 2 and 3 (weak site occupied) will have near-zero probability
- Effectively behaves as a 2-state system (state 0 ↔ state 1) while using the shared codebase
- N = 10,000 cells (first 20 save traces)
- Scientific expectation: bimodal (ON/OFF) but NO intermediate state
- **Shows that two-site architecture, not just repression, generates the intermediate subpopulation**

**Condition D — No Regulation (Poisson baseline):**
- **MUST use the same 4-state Numba engine as all other conditions.** Do NOT build a separate engine.
- Implement no-regulation numerically: set all 4 transcription rates to k_max (i.e., k_txn = [k_max, k_max, k_max, k_max]) AND set all binding propensities to zero (k_on = 0, so no operator transitions ever occur — operator stays in initial state 0 and all states produce the same output regardless).
- N = 10,000 cells
- Scientific expectation: Fano factor ≈ 1 + burst_size (Poisson-like, unimodal)
- This is the negative control confirming that the operator model is what creates excess noise

**Condition E — Asymmetry Sweep:**
- Kd_strong = 2.4 nM (fixed)
- Kd_weak swept: ratio = Kd_weak/Kd_strong ∈ {1, 2, 5, 8, 10, 15, 20, 30, 40, 50}
- So Kd_weak ∈ {2.4, 4.8, 12, 19.2, 24, 36, 48, 72, 96, 120} nM
- **Decision: keep block_strong=0.85, block_weak=0.50 fixed. Only sweep Kd_weak.**
- **NOTE:** At ratio=1, the sweep still has asymmetric blocking (0.85 vs 0.50) even though Kd values are equal. This is intentionally different from Condition B (which uses symmetric blocking 0.85/0.85). The sweep isolates the effect of binding affinity asymmetry alone, while Condition B represents full symmetry in both affinity and blocking.
- N = 1,000 cells per ratio (10 ratios × 1,000 = 10,000 cells total)
- Output: intermediate_fraction vs asymmetry_ratio

**Output files:**
- `results/phase2/condition_A_asymmetric.npz` — arrays: protein_final (10000,), mRNA_final (10000,), operator_states_final (10000,)
- `results/phase2/condition_B_symmetric.npz` — same structure
- `results/phase2/condition_C_single_site.npz` — same structure
- `results/phase2/condition_D_unregulated.npz` — same structure
- `results/phase2/condition_E_sweep.npz` — arrays: ratios (10,), protein_final (10, 1000), intermediate_fractions (10,)
- `results/phase2/traces/traces_A.npz` — arrays: times (T,), proteins (20, T) for 20 cells
- `results/phase2/traces/traces_B.npz` — same
- `results/phase2/traces/traces_C.npz` — same
- `results/phase2/steady_state_distributions.csv` — summary table with mean, variance, CV, Fano for each condition

**Self-checks (per condition):**
- All 10,000 (or 1,000) cells completed without error
- Mean protein > 0 (cells are expressing)
- Variance > 0 (not deterministic)
- Condition D: Fano factor between 5 and 12 (expected ~7.85)
- Condition A: CV > Condition B CV (asymmetric is noisier)
- Condition A: CV > Condition C CV (asymmetric noisier than single-site)
- Condition C: GMM selects 2 components (bimodal, not trimodal)
- Condition E: intermediate fraction at ratio=1 < intermediate fraction at ratio=50 (monotonic trend)
- No NaN or Inf values in any output

### 7.4 asymmetry_sweep.py

**Purpose:** Dedicated sweep of asymmetry ratio — extracted from run_conditions.py for clarity. This is Condition E.

**Algorithm:**
```python
RATIOS = [1, 2, 5, 8, 10, 15, 20, 30, 40, 50]

for ratio in RATIOS:
    params['Kd_weak'] = params['Kd_strong'] * ratio
    # Run 1000 cells
    proteins = [simulate_cell(params) for _ in range(1000)]
    # Calculate intermediate fraction using GMM
    intermediate_frac = calculate_intermediate_fraction(proteins)
```

**Persister threshold definition (GMM intersection method):**

**Primary method (GMM intersection):**
- Fit 2-component or 3-component GMM to Condition A protein distribution (whichever BIC selects)
- Find intersection point(s) between the Gaussian components
- Persister threshold = intersection point between the lowest and middle components (for 3-component) or between the two components (for 2-component)
- Cells above this threshold are in the intermediate or high state = predicted persisters

**Secondary method (report alongside for robustness):**
- Threshold = mean_protein_condition_B + 2 × std_protein_condition_B
- Report: "Persister fraction is X% (GMM intersection) or Y% (2σ threshold)"

**Plot annotation:** Mark the real Mce3R asymmetry ratio (~20.4, since 49/2.4 ≈ 20.4) with a vertical dashed line labeled "Mce3R (Panagoda 2024)".

**Self-checks:**
- 10 ratio values × 1,000 cells = 10,000 simulations completed
- Intermediate fraction increases with ratio (positive correlation)
- Intermediate fraction at ratio=1 is < 5% (symmetric case)
- Intermediate fraction at ratio=50 is > intermediate fraction at ratio=1

### 7.5 simulation_main.py

**Purpose:** Phase 2 orchestrator.

```python
def run_phase2():
    """Execute all simulation conditions."""
    # Step 1: Initialize operator model, pre-compute rate arrays for Numba
    # Step 2: Run Condition D (unregulated) — fastest, validates engine
    # Step 3: Run Condition B (symmetric) — establishes baseline
    # Step 4: Run Condition C (single site) — minimal control
    # Step 5: Run Condition A (asymmetric) — main result
    # Step 6: Run Condition E (sweep) — dose-response
    # Step 7: Save all outputs including traces
    # Step 8: Log summary with PASS/FAIL checks
```

---

## 8. PHASE 3: STATISTICAL VALIDATION

<!-- DEPENDS_ON: Phase 2 -->

### 8.1 noise_metrics.py

**Purpose:** Calculate noise metrics for each condition (A, B, C, D).

**Metrics computed:**
- **Mean (μ):** average protein level
- **Variance (σ²):** variance of protein levels
- **Coefficient of Variation (CV):** σ/μ — normalized noise
- **Fano Factor:** σ²/μ — excess noise over Poisson
- **Bimodality Coefficient:** (skewness² + 1) / (kurtosis + 3 × (n-1)²/((n-2)(n-3))) — values > 0.555 suggest bimodality
- **Number of modes:** from Gaussian Mixture Model (GMM) fitting with BIC selection

**GMM ROBUSTNESS FOR DISCRETE DATA:**
Protein counts from the Gillespie simulation are discrete integers. Feeding highly skewed discrete data into sklearn.mixture.GaussianMixture often causes singular covariance matrices or convergence failures.

**REQUIRED:** Before fitting any GMM, apply BOTH of these safeguards:
1. Add continuous dither noise to break ties:
    ```python
    data_continuous = data.astype(float) + np.random.uniform(-0.5, 0.5, size=len(data))
    ```
2. Initialize GaussianMixture with reg_covar=1e-3:
    ```python
    gmm = GaussianMixture(n_components=k, reg_covar=1e-3, n_init=5, random_state=42)
    ```

Also wrap GMM fitting in a try/except — if fitting fails for a given n_components, assign that model BIC = +inf so it is never selected:
```python
try:
    gmm.fit(data_continuous.reshape(-1, 1))
    bic = gmm.bic(data_continuous.reshape(-1, 1))
except Exception:
    bic = np.inf
```

**GMM fitting:**
- Fit 1, 2, and 3-component GMMs using `sklearn.mixture.GaussianMixture`
- Select best model by BIC (Bayesian Information Criterion)
- Report: number of components, means, variances, weights, BIC values

**Output:** `results/phase3/noise_metrics.csv`

| Condition | Mean | Variance | CV | Fano | Bimodality_Coeff | N_modes | GMM_BIC_1 | GMM_BIC_2 | GMM_BIC_3 |
|-----------|------|----------|----|----|------------------|---------|-----------|-----------|-----------|

**Self-checks:**
- CV(A) > CV(B) — asymmetric noisier than symmetric
- CV(A) > CV(C) — asymmetric noisier than single-site
- Fano(D) ≈ 1 + burst_size — unregulated matches theory
- Fano(A) > Fano(B) — asymmetric has more excess noise
- Bimodality coefficient: A > B
- All values are finite and positive

### 8.2 bootstrap_ci.py

**Purpose:** Bootstrap confidence intervals for all noise metrics.

**Algorithm:**
```python
def bootstrap_metric(data, metric_func, n_bootstrap=1000, ci=0.95):
    """
    Bootstrap confidence interval for a metric.

    Args:
        data: array of protein levels (N cells)
        metric_func: function that computes the metric (e.g., np.std/np.mean for CV)
        n_bootstrap: number of bootstrap resamples
        ci: confidence level

    Returns:
        (point_estimate, ci_lower, ci_upper)
    """
    n = len(data)
    boot_values = []
    for _ in range(n_bootstrap):
        resample = np.random.choice(data, size=n, replace=True)
        boot_values.append(metric_func(resample))

    point = metric_func(data)
    alpha = (1 - ci) / 2
    ci_lower = np.percentile(boot_values, alpha * 100)
    ci_upper = np.percentile(boot_values, (1 - alpha) * 100)
    return point, ci_lower, ci_upper
```

- Apply to: CV, Fano factor, intermediate fraction, persister fraction
- For each condition (A, B, C, D)
- n_bootstrap = 1,000
- Report 95% CIs

**Output:** `results/phase3/bootstrap_results.csv`

**Self-checks:**
- CI_lower < point_estimate < CI_upper for all metrics
- CI width is reasonable (not zero, not wider than the mean)
- CIs for CV(A) and CV(B) do not overlap (significant difference)

### 8.3 sensitivity_analysis.py

**Purpose:** Assess robustness of results to parameter choices.

**Parameters to vary (one at a time, ±50% in 10 steps):**

| Parameter | Baseline | Range |
|-----------|----------|-------|
| k_max | 0.15 | 0.075 – 0.225 |
| k_translation | 0.5 | 0.25 – 0.75 |
| t_half_mRNA | 9.5 min | 4.75 – 14.25 |
| t_half_protein | 1500 min | 750 – 2250 |
| block_strong | 0.85 | 0.70 – 0.95 |
| block_weak | 0.50 | 0.30 – 0.70 |
| k_on | 0.0167 | 0.00835 – 0.02505 |

**Note:** R_conc is no longer a parameter (it's dynamic). 7 parameters total.

**Algorithm:**
1. For each parameter:
   - Create 10 evenly-spaced values across the range
   - For each value, run 1,000 cells with the asymmetric model
   - Compute CV and number of GMM modes
2. Store results for tornado plot and heatmap

**Outputs:**
- `results/phase3/sensitivity_data.csv` — columns: parameter, value, CV, Fano, n_modes, intermediate_fraction
- Data for tornado plot: for each parameter, the change in CV when varied from -50% to +50%

**Self-checks:**
- 7 parameters × 10 steps × 1,000 cells = 70,000 simulations completed
- No NaN values in output
- At least 4 of 7 parameters: multimodality persists across >50% of the range (robustness)
- Tornado plot: at least one parameter has >20% impact on CV (to identify key drivers)

### 8.4 statistical_tests.py

**Purpose:** Formal statistical comparison of conditions.

**Tests:**

1. **Kolmogorov-Smirnov test (2-sample):**
   - Compare protein distributions: A vs B, A vs C, A vs D, B vs D, C vs D
   - `scipy.stats.ks_2samp(proteins_A, proteins_B)`
   - Report: KS statistic, p-value

2. **GMM model selection and likelihood ratio (heuristic):**
   - **PRIMARY criterion: BIC (Bayesian Information Criterion)** — select the number of GMM components (1, 2, or 3) that minimizes BIC. This is the authoritative model selection method.
   - **Secondary/descriptive: Likelihood ratio statistic** — reported alongside BIC for additional context, but NOT used as the primary decision criterion.
   - **CRITICAL implementation note:** `sklearn.mixture.GaussianMixture.score(X)` returns the **average** log-likelihood per sample, NOT total.
   - Must compute total log-likelihood: `LL_total = model.score(X) * N` where N = len(X)
   - LR statistic = 2 × (LL3_total - LL2_total)
   - Do NOT use raw `.score(X)` values directly in the LR formula without multiplying by N.
   - Degrees of freedom = 3 (one extra mean + variance + weight)
   - p-value from chi-squared distribution: `scipy.stats.chi2.sf(LR_stat, df=3)`
   - **Limitation note:** This is an approximate heuristic comparison. The regularity conditions for a standard chi-square LR test are not strictly satisfied at GMM component boundaries (parameters on boundary of parameter space). Report the p-value with this caveat.

3. **Effect size (Cohen's d):**
   - Between mean protein of condition A vs B, A vs C
   - **CRITICAL: Use correct pooled standard deviation, NOT averaged SDs:**
     ```python
     pooled_sd = np.sqrt(((n1 - 1) * var1 + (n2 - 1) * var2) / (n1 + n2 - 2))
     d = (mean1 - mean2) / pooled_sd
     ```
   - Guard against divide-by-zero if pooled_sd == 0 (return d = 0.0 or np.inf with warning)

**Output:** `results/phase3/statistical_tests.csv`

| Test | Comparison | Statistic | p_value | Significant |
|------|-----------|-----------|---------|------------|

**Self-checks:**
- KS test A vs B: p < 0.05 (distributions are significantly different)
- All p-values are between 0 and 1
- Report honestly if any test yields p > 0.05

### 8.5 experimental_comparison.py

**Purpose:** Compare simulation predictions to published experimental data.

**Comparisons:**

1. **Persister frequency:**
   - Our prediction: fraction of cells above persister threshold in condition A
   - **Use GMM intersection threshold (primary)** and 2σ threshold (secondary)
   - Experimental: Pandey et al. 2023 reported increased persister frequency in Δmce3R
   - Note: we predict the WILD-TYPE persister fraction (repressor present); Δmce3R would have no repression (closer to condition D)
   - Report: "Wild-type predicted persister fraction: X% (GMM) / Y% (2σ). Unregulated (Δmce3R proxy) fraction: Z%. Ratio: W."

2. **Gene expression fold-change:**
   - Santangelo et al. 2009 reported fold-changes for Mce3R regulon genes in Δmce3R vs WT
   - Our model predicts mean expression in regulated (condition A) vs unregulated (condition D)
   - Predicted fold-change = mean_protein_D / mean_protein_A
   - Compare to published values (stored in `data/external/santangelo2009_foldchange.csv`)

**Output section in** `results/phase3/noise_metrics.csv` (appended) and logged.

**Self-checks:**
- Predicted persister fraction is between 0.01% and 10% (biologically reasonable range)
- Predicted fold-change is between 2 and 50 (reasonable derepression)
- Explicit caveat logged: "Simplified model; quantitative match not expected"

### 8.6 negative_controls.py

**Purpose:** Confirm results are not artifacts.

**Control 1 — Sequence shuffle test (for Phase 1):**
- For each top 10 predicted binding site:
  - Shuffle the DNA sequence 1,000 times (randomly permuting the sequence while preserving mononucleotide composition, i.e., use `numpy.random.permutation` on the list of characters)
  - Re-score each shuffle against the MEME motif using the PWM
  - Real sequence score must exceed 99th percentile of shuffled scores (p < 0.01)

**Control 2 — Symmetric TetR control:**
- Run the Gillespie simulation with parameters from a generic symmetric TetR system:
  - Kd_strong = Kd_weak = 5 nM (typical TetR palindrome)
  - block_strong = block_weak = 0.85
  - All other parameters same
- Confirm: does NOT produce multimodality (GMM selects 1 component)

**Control 3 — Poisson baseline verification:**
- Condition D (no regulation) should have Fano factor ≈ 1 + burst_size
- burst_size = k_translation / gamma_mRNA = 0.5 / 0.073 ≈ 6.85
- Expected Fano ≈ 7.85
- Accept if Fano is within [5, 12]

**Output:** `results/phase3/negative_controls.csv`

**Self-checks:**
- Shuffle test: ≥8 of top 10 sites pass at p < 0.01
- Symmetric TetR: GMM selects 1 component (not bimodal)
- Poisson baseline: Fano within [5, 12]

### 8.7 analysis_main.py

**Purpose:** Phase 3 orchestrator.

```python
def run_phase3():
    """Execute all statistical analyses."""
    # Step 1: Compute noise metrics (with GMM robustness safeguards)
    # Step 2: Bootstrap CIs
    # Step 3: Sensitivity analysis (longest step — 70K simulations)
    # Step 4: Statistical tests
    # Step 5: Experimental comparison (using GMM intersection threshold)
    # Step 6: Negative controls
    # Step 7: Log summary
```

---

## 9. PHASE 4: PUBLICATION FIGURES

<!-- DEPENDS_ON: Phase 3 -->

### 9.1 figure_style.py

**Shared style settings for all figures:**

```python
import matplotlib.pyplot as plt
import matplotlib as mpl

STYLE = {
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'font.size': 10,
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'Helvetica'],
    'axes.labelsize': 12,
    'axes.titlesize': 13,
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
    'legend.fontsize': 9,
    'figure.figsize': (7, 5),
    'axes.linewidth': 1.0,
    'lines.linewidth': 1.5,
    'axes.spines.top': False,
    'axes.spines.right': False,
}

COLORS = {
    'asymmetric': '#0D9488',     # Teal
    'symmetric': '#F97316',      # Coral/Orange
    'single_site': '#8B5CF6',    # Purple
    'unregulated': '#6B7280',    # Gray
    'sweep_cmap': 'viridis',
    'mce3r_marker': '#E9C46A',   # Gold for Mce3R point on sweep
}
```

### 9.2 Figure Descriptions

**Figure 1 — Genome-wide Binding Site Map** (`fig1_binding_sites.py`)
- Circular genome plot of H37Rv with predicted Mce3R binding sites marked
- Color-code by FIMO p-value (darker = more significant)
- Label known Mce3R regulon genes
- Inset: bar chart of top 15 predicted sites ranked by p-value
- Data source: `results/phase1/predicted_sites.csv`

**Figure 2 — Protein Distributions: All Conditions** (`fig2_distributions.py`)
- 2×3 or 3×2 panel figure:
  - (A) Histogram of protein levels, Condition A (asymmetric) with GMM overlay
  - (B) Histogram of protein levels, Condition B (symmetric) with GMM overlay
  - (C) Histogram of protein levels, Condition C (single site) with GMM overlay
  - (D) Histogram of protein levels, Condition D (unregulated)
  - (E) Box plot comparing CV, Fano, bimodality coefficient across A, B, C, D with bootstrap error bars
  - (F) Overlay of all conditions on single axis (normalized density)
- Persister threshold (GMM intersection) marked as vertical dashed line on histograms
- Data source: `results/phase2/condition_*.npz`, `results/phase3/noise_metrics.csv`

**Figure 3 — Asymmetry Sweep** (`fig3_asymmetry_sweep.py`)
- X-axis: asymmetry ratio (Kd_weak/Kd_strong), log scale
- Y-axis: intermediate state fraction (or persister fraction)
- Each point = 1,000 cells; error bars from bootstrap CIs
- Vertical dashed line at ratio ≈ 20.4 labeled "Mce3R (Panagoda 2024)"
- Show monotonic increase: more asymmetry → more intermediates
- Data source: `results/phase2/condition_E_sweep.npz`

**Figure 4 — Parameter Sensitivity** (`fig4_sensitivity.py`)
- (A) Tornado plot: horizontal bars showing impact of each parameter (±50%) on CV
  - Parameters sorted by impact magnitude
  - Teal bars = increase CV, orange bars = decrease CV
- (B) Heatmap: parameter value (x) vs parameter name (y), color = number of GMM modes
  - Shows which parameter regions produce multimodality
- Data source: `results/phase3/sensitivity_data.csv`

**Figure 5 — Validation & Controls** (`fig5_validation.py`)
- (A) Scatter plot: predicted fold-change vs published fold-change (Santangelo 2009)
  - Diagonal line = perfect agreement; annotate R² value
- (B) Bar chart: predicted persister fraction (WT) vs unregulated (Δmce3R proxy)
  - With bootstrap error bars; show both GMM and 2σ thresholds
- (C) Shuffle test results: histogram of shuffled scores with real score marked
  - For the top binding site
- (D) Symmetric TetR control distribution vs asymmetric Mce3R distribution
- Data source: `results/phase3/negative_controls.csv`, `results/phase3/noise_metrics.csv`

**Figure 6 — Single-Cell Traces** (`fig6_single_cell_traces.py`)
- 3-panel figure: one panel each for Conditions A, B, C
- Each panel: protein count (y-axis) vs time in minutes (x-axis) for 3 representative cells
- Show ~2,100 minutes of post-burn-in trajectory (the full recording window)
- For Condition A, select cells that show switching between states (visually demonstrates stochastic transitions)
- Thin lines, transparency for overlap
- Color: match condition colors from COLORS dict
- Data source: `results/phase2/traces/traces_*.npz`

**Figure 7 — BIC Comparison** (`fig7_bic_comparison.py`)
- Grouped bar chart
- X-axis: condition (A, B, C, D)
- Y-axis: BIC value
- 3 bars per condition: 1-component GMM, 2-component GMM, 3-component GMM
- Lowest BIC = best model (annotate with star)
- This directly supports the multimodality claim
- Data source: `results/phase3/noise_metrics.csv` (GMM_BIC columns)

### 9.3 figures_main.py

```python
def run_phase4():
    """Generate all 7 publication figures."""
    apply_style()
    # Generate each figure (fig1 through fig7)
    # Save to results/figures/ as 300 dpi PNG
    # Log which figures were generated successfully
```

---

## 10. MAIN ORCHESTRATOR (main.py)

**Purpose:** Top-level script that runs the entire project end-to-end with runtime parallelism for Phase 1 + Phase 2.

```python
#!/usr/bin/env python3
"""
main.py — Orchestrator for the Mce3R Asymmetric Operator Project

Usage:
    python main.py              # Run everything (Phase 1+2 parallel, then 3, then 4)
    python main.py --phase 1    # Run only Phase 1
    python main.py --phase 2    # Run only Phase 2
    python main.py --phase 3    # Run only Phase 3 (requires Phase 2 output)
    python main.py --phase 4    # Run only Phase 4 (requires Phase 3 output)
"""

import argparse
import os
import subprocess
import sys
import time
from datetime import datetime

def main():
    parser = argparse.ArgumentParser(description='Mce3R Asymmetric Operator Project')
    parser.add_argument('--phase', type=int, choices=[1, 2, 3, 4],
                        help='Run specific phase only')
    args = parser.parse_args()

    # Create directory structure
    create_directories()

    # Initialize logging
    setup_logging()

    build_status = {}
    start_time = time.time()

    try:
        if args.phase is None:
            # PARALLEL EXECUTION: Phase 1 and Phase 2 simultaneously
            print("[Orchestrator] Launching Phase 1 + Phase 2 in parallel...")
            p1 = subprocess.Popen(
                [sys.executable, '-m', 'phase1_pipeline.pipeline_main'],
                stdout=open('logs/phase1_pipeline.log', 'w'),
                stderr=subprocess.STDOUT
            )
            p2 = subprocess.Popen(
                [sys.executable, '-m', 'phase2_simulation.simulation_main'],
                stdout=open('logs/phase2_gillespie.log', 'w'),
                stderr=subprocess.STDOUT
            )
            p1.wait()
            p2.wait()
            build_status['phase1'] = {'returncode': p1.returncode}
            build_status['phase2'] = {'returncode': p2.returncode}

            if p1.returncode != 0 or p2.returncode != 0:
                raise RuntimeError(f"Parallel phases failed: P1={p1.returncode}, P2={p2.returncode}")

            # SEQUENTIAL: Phase 3 depends on Phase 2
            print("[Orchestrator] Starting Phase 3: Statistical Validation")
            from phase3_analysis.analysis_main import run_phase3
            build_status['phase3'] = run_phase3()

            # SEQUENTIAL: Phase 4 depends on Phase 3
            print("[Orchestrator] Starting Phase 4: Publication Figures")
            from phase4_figures.figures_main import run_phase4
            build_status['phase4'] = run_phase4()

        elif args.phase == 1:
            from phase1_pipeline.pipeline_main import run_phase1
            build_status['phase1'] = run_phase1()
        elif args.phase == 2:
            from phase2_simulation.simulation_main import run_phase2
            build_status['phase2'] = run_phase2()
        elif args.phase == 3:
            verify_phase2_outputs_exist()
            from phase3_analysis.analysis_main import run_phase3
            build_status['phase3'] = run_phase3()
        elif args.phase == 4:
            verify_phase3_outputs_exist()
            from phase4_figures.figures_main import run_phase4
            build_status['phase4'] = run_phase4()

    except Exception as e:
        build_status['error'] = str(e)
        print(f"[Orchestrator] FATAL ERROR: {e}")

    # Write BUILD_REPORT.md
    elapsed = time.time() - start_time
    write_build_report(build_status, elapsed)
    print(f"[Orchestrator] Complete. Total time: {elapsed/60:.1f} minutes")
    print(f"[Orchestrator] See BUILD_REPORT.md for details")

if __name__ == '__main__':
    main()
```

---

## 11. SELF-TESTING REQUIREMENTS

### 11.1 Two Classes of Checks

Every script must include TWO distinct classes of checks in its `if __name__ == "__main__"` block:

**CLASS 1 — SOFTWARE / SANITY CHECKS (hard-fail the build if any fails):**
- No exceptions during execution
- Output files exist and are non-empty
- Array lengths/dimensions correct
- No NaN or Inf values
- Data types correct
- Values within physically possible ranges
- Numba compilation succeeds

**CLASS 2 — SCIENTIFIC EXPECTATION CHECKS (logged but do NOT hard-fail):**
- CV(A) > CV(B), bimodality/trimodality outcomes
- Monotonic asymmetry trend
- Non-overlapping CIs
- Fano factor near expected theoretical value
- Any outcome that depends on the scientific hypothesis being correct

Scientific expectation checks are logged as `EXPECTED: [description] | OBSERVED: [value]` and included in BUILD_REPORT.md, but they only produce a WARNING, not a build-stopping FAIL. If a scientific check fails, it may indicate the hypothesis is wrong or parameters need adjustment — that is a scientific finding, not a software bug.

### 11.2 Self-Test Template

```python
if __name__ == "__main__":
    """Self-test block."""
    import sys

    n_sanity_pass = 0
    n_sanity_fail = 0
    n_science_expected = 0
    n_science_unexpected = 0

    # === SOFTWARE / SANITY CHECKS (hard-fail) ===
    try:
        result = function_under_test(known_input)
        if expected_lower <= result <= expected_upper:
            print(f"SANITY PASS: [test name] — {result}")
            n_sanity_pass += 1
        else:
            print(f"SANITY FAIL: [test name] — {result} (expected {expected_lower}-{expected_upper})")
            n_sanity_fail += 1
    except Exception as e:
        print(f"SANITY FAIL: [test name] — Exception: {e}")
        n_sanity_fail += 1

    # === SCIENTIFIC EXPECTATION CHECKS (warn only) ===
    if cv_a > cv_b:
        print(f"EXPECTED: CV(A) > CV(B) | OBSERVED: {cv_a:.3f} > {cv_b:.3f}")
        n_science_expected += 1
    else:
        print(f"UNEXPECTED: CV(A) <= CV(B) | OBSERVED: {cv_a:.3f} <= {cv_b:.3f}")
        n_science_unexpected += 1

    # Summary
    print(f"\n{'='*50}")
    print(f"Sanity checks: {n_sanity_pass} PASS, {n_sanity_fail} FAIL")
    print(f"Scientific checks: {n_science_expected} EXPECTED, {n_science_unexpected} UNEXPECTED")
    if n_sanity_fail > 0:
        print("STOPPING: Software sanity check failed.")
        sys.exit(1)
    elif n_science_unexpected > 0:
        print("WARNING: Some scientific expectations not met. Review results.")
        sys.exit(0)  # Do NOT fail the build
    else:
        print("All checks passed.")
        sys.exit(0)
```

### 11.3 Required Checks Per Module

| Module | Sanity Checks (hard-fail) | Scientific Expectations (warn) |
|--------|--------------------------|-------------------------------|
| operator_model.py | All rates positive; k_off values match calculation; Numba array shapes correct | Transcription rates monotonic; detailed balance holds |
| gillespie_engine.py | No exceptions; Numba compilation succeeds; mRNA in [0,50]; protein in [0,5000] | Mean mRNA ≈ 2 (no regulation); Fano ≈ 7.85 (no regulation) |
| run_conditions.py | All 5 conditions complete; no NaN/Inf; output files exist | CV(A) > CV(B); CV(A) > CV(C); Fano(D) in [5,12] |
| asymmetry_sweep.py | 10 ratios × 1000 cells complete; no NaN | Monotonic intermediate fraction trend |
| noise_metrics.py | All metrics finite; GMM fitting doesn't crash | CV(A) > CV(B); bimodality coefficient A > B |
| bootstrap_ci.py | CI_lower < point < CI_upper (mathematical requirement) | CIs for CV(A) and CV(B) do not overlap |
| sensitivity_analysis.py | No NaN values; correct array dimensions | Multimodality persists in >50% of parameter range |
| statistical_tests.py | All p-values in [0,1]; no exceptions | KS p-value A vs B < 0.05 |
| negative_controls.py | No exceptions; outputs parseable | Symmetric TetR unimodal; Poisson Fano in [5,12] |
| download_genomes.py | Files exist; H37Rv length ≈ 4,411,532 | GC ≈ 65.6% |
| extract_upstream.py | 3500-4500 sequences; all exactly 200bp; wrapping tests pass | Known operator found in Rv1964 upstream |
| run_meme.py | Output dir exists (or mock flagged) | If real: E-value < 0.05, width 20-30bp |
| run_fimo.py | Output parseable (or mock flagged) | If real: known operator in top 5 |
| conservation_check.py | Output CSV has correct columns and row count | Known operator conserved in both species |

---

## 12. PARAMETER REFERENCE TABLE

All parameters are defined in `config/parameters.py` — the single source of truth.

```python
"""
config/parameters.py — All tunable parameters for the Mce3R project.

Every parameter includes its source reference.
"""

PARAMS = {
    # === Genome Information ===
    'H37Rv_accession': 'NC_000962.3',         # NCBI RefSeq
    'H37Rv_genome_size': 4_411_532,            # bp
    'H37Rv_gc_content': 0.656,                 # fraction
    'M_bovis_accession': 'NC_002945.4',        # NCBI RefSeq
    'M_marinum_accession': 'NC_010612.1',      # NCBI RefSeq

    # === Mce3R Binding Parameters (Panagoda et al. 2024) ===
    'Kd_strong': 2.4,          # nM, Table 1, Probe A
    'Kd_weak': 49.0,           # nM, Table 1, Probe C — direct measurement, most defensible
    'k_on': 0.0167,            # nM^-1 min^-1, Stormo & Zhao 2010 (10^6 M^-1 s^-1 converted)
    # NOTE: R_conc is NOT a parameter. Repressor concentration is dynamic,
    # computed from current free protein count at each Gillespie step.

    # === Cell Volume & Concentration Conversion ===
    'cell_volume_fL': 1.0,           # femtoliters, standard Mtb cell volume
    'nM_per_molecule': 1.66,         # 1 molecule in 1 fL ≈ 1.66 nM
                                     # (from 1/(Avogadro × 1e-15 L) × 1e9)

    # === Transcription Model ===
    'k_max': 0.15,             # mRNA/min, estimated from Mtb transcriptomics
    'block_strong': 0.85,      # fraction blocked when strong site occupied
    'block_weak': 0.50,        # fraction blocked when weak site occupied

    # === Translation & Degradation (Rustad et al. 2013; Mtb literature) ===
    'k_translation': 0.5,      # protein per mRNA per min, Taniguchi 2010 adjusted
    't_half_mRNA': 9.5,        # min, Rustad et al. 2013 NAR
    't_half_protein': 1500.0,  # min (~25 hr), dominated by Mtb dilution rate

    # === Simulation Settings ===
    'n_cells_main': 10_000,     # cells per condition (A, B, C, D)
    'n_cells_sweep': 1_000,     # cells per sweep point (Condition E)
    't_max': 4_200,             # min (70 hours ≈ 3.5 Mtb doubling times at ~20hr/doubling)
    't_burn_in': 2_100,         # min, discard first half as transient

    # === Trace Recording (for fig6) ===
    'n_trace_cells': 20,        # number of cells to save full trajectories for
    'trace_interval': 10.0,     # minutes between recorded time points

    # === Asymmetry Sweep (Condition E) ===
    'sweep_ratios': [1, 2, 5, 8, 10, 15, 20, 30, 40, 50],
    'mce3r_actual_ratio': 20.4,  # Kd_weak/Kd_strong = 49.0/2.4 ≈ 20.4

    # === Symmetric Control (Condition B) ===
    'Kd_symmetric': 2.4,        # nM, both sites use strong affinity
    'block_symmetric': 0.85,    # both sites block equally

    # === Sensitivity Analysis ===
    'sensitivity_range': 0.50,  # ±50% variation
    'sensitivity_steps': 10,    # number of steps per parameter
    'sensitivity_cells': 1_000, # cells per step

    # === Statistical Parameters ===
    'n_bootstrap': 1_000,       # bootstrap resamples
    'ci_level': 0.95,           # confidence interval level
    'n_gmm_components': [1, 2, 3],  # GMM components to test
    'gmm_reg_covar': 1e-3,     # regularization for GMM on discrete data
    'gmm_n_init': 5,           # number of GMM initializations
    'shuffle_n': 1_000,         # shuffles per binding site
    'shuffle_alpha': 0.01,      # significance threshold for shuffle test

    # === Reproducibility Seeds ===
    'master_seed': 42,             # Master seed for all stochastic operations
    'gmm_random_state': 42,        # Passed to GaussianMixture(random_state=...)
    'bootstrap_seed': 123,         # Seed for bootstrap resampling
    # Condition-level seeds are derived: seed_condition_X = master_seed + condition_index
    # Sweep-level seeds: seed_sweep_ratio_i = master_seed + 100 + i
    # This ensures reproducibility: identical seeds → identical results.
    # To generate a different realization, change master_seed only.

    # === MEME/FIMO Parameters ===
    'meme_mod': 'zoops',        # zero or one occurrence per sequence
    'meme_minw': 20,            # minimum motif width
    'meme_maxw': 30,            # maximum motif width
    'meme_nmotifs': 3,          # number of motifs to find
    'fimo_thresh': 1e-4,        # FIMO p-value threshold
    'upstream_length': 200,     # bp upstream of gene start

    # === Operator Sequence (PDB 9B7Y, Panagoda 2024) ===
    'operator_sequence': (
        'GCCCCGCGCTATAGGATACTAGCAAGATACATCATAGCCAATATATGCCAGT'
        'TTGCATTGCTATTTACCGATCAGTTGTCCAAGCAATCGCGTATTGGCTATG'
        'GACATCAGCGGTTCTGCCGC'
    ),
    'operator_length': 123,     # bp
    'operator_region_h37rv_start': 2_207_477,  # approximate intergenic region start (NOT exact 123-bp boundary)
    'operator_region_h37rv_end': 2_207_699,    # approximate intergenic region end (NOTE: span is ~222 bp, larger than 123-bp operator; this is the broader intergenic window containing the operator)

    # === Gene Identifiers ===
    'mce3R_gene': 'Rv1963c',
    'yrbE3A_gene': 'Rv1964',
    'mce3_operon_genes': [f'Rv{i}' for i in range(1964, 1978)],
    'regulon_genes_1': ['Rv1933c', 'Rv1934c', 'Rv1935c'],
    'regulon_genes_2': [f'Rv{i}' for i in range(1936, 1942)],
}

# Derived parameters (calculated, not configurable)
import numpy as np
PARAMS['gamma_mRNA'] = np.log(2) / PARAMS['t_half_mRNA']       # 0.0730 min^-1
PARAMS['gamma_protein'] = np.log(2) / PARAMS['t_half_protein']  # 0.000462 min^-1
PARAMS['burst_size'] = PARAMS['k_translation'] / PARAMS['gamma_mRNA']  # ~6.85
PARAMS['expected_fano_unregulated'] = 1 + PARAMS['burst_size']  # ~7.85

# Derived binding rate constants
PARAMS['k_off_strong'] = PARAMS['Kd_strong'] * PARAMS['k_on']  # 2.4 × 0.0167 = 0.0401 min^-1
PARAMS['k_off_weak'] = PARAMS['Kd_weak'] * PARAMS['k_on']      # 49.0 × 0.0167 = 0.818 min^-1

# Transcription rates per operator state
PARAMS['k_txn'] = np.array([
    PARAMS['k_max'],                                                                    # State 0: both empty
    PARAMS['k_max'] * (1 - PARAMS['block_strong']),                                     # State 1: strong occupied
    PARAMS['k_max'] * (1 - PARAMS['block_weak']),                                       # State 2: weak occupied
    PARAMS['k_max'] * (1 - PARAMS['block_strong']) * (1 - PARAMS['block_weak']),        # State 3: both occupied
])
```

---

## 13. LITERATURE REFERENCES

**CITATION HANDLING RULE:** Do not modify, "fix," or re-verify literature citations during code generation. Use citations exactly as written in this document. Treat citation text as frozen project metadata, not something to reinterpret.

<!-- CITATION VERIFICATION NEEDED (for manual human review only — NOT for code to act on):
1. Santangelo 2009: Handoff says "Microbiology 155:2245-2255", architecture says "Microbiology 155:882-891". One is wrong. Verify correct page numbers before submission.
2. Persister phenotype paper: Handoff cites "Srivastava et al. 2023, Microbiol. Spectrum". Architecture cites "Pandey et al. 2023, Res. Microbiol. 174:104082". Determine if these are the same finding or two separate papers. If two papers, cite both.
3. Santangelo 2002 vs 2008: Architecture Section 2.1 cites Santangelo 2002 (Microbiology 148:3044-3055) for original characterization, but handoff cites Santangelo 2008 (BMC Microbiol. 8:38) for the finding that mce3R represses mce3 but NOT mce1/2/4. These may be different papers — keep both if so.
-->

| # | Reference | Used For |
|---|-----------|----------|
| 1 | Panagoda NT, Balázsi G, Sampson NS (2024) ACS Chem. Biol. 19:2580-2592 | Asymmetric operator, Kd values, operator sequence, PDB 9B7Y |
| 2 | Pandey M et al. (2023) Res. Microbiol. 174:104082 | Δmce3R increases persister frequency |
| 3 | Flentie K et al. (2019) ACS Infect. Dis. 5(7):1239-1254 | Mce3R pathway enhances antibiotic activity |
| 4 | Santangelo MP et al. (2002) Microbiology 148:3044-3055 | Original Mce3R characterization |
| 5 | Santangelo MP et al. (2009) Microbiology 155:882-891 | Mce3R regulon, fold-change data |
| 6 | Santangelo MP et al. (2008) BMC Microbiol. 8:38 | mce3R represses mce3 but not mce1/2/4 |
| 7 | Balázsi G, van Oudenaarden A, Collins JJ (2011) Cell 144:910-925 | Noise → cell fate decisions review |
| 8 | Farquhar KS et al. (2019) Nat. Commun. 10:2766 | Noise circuits drive drug resistance |
| 9 | Sureka K et al. (2008) PLoS ONE 3(3):e1771 | Bimodal rel expression in mycobacteria |
| 10 | Rustad TR et al. (2013) Nucleic Acids Res. 41:509-517 | Mtb mRNA half-life = 9.5 min |
| 11 | Taniguchi Y et al. (2010) Science 329:533-538 | E. coli gene expression noise parameters |
| 12 | Ramos JL et al. (2005) Microbiol. Mol. Biol. Rev. 69:326-356 | TetR family review |
| 13 | Gillespie DT (1977) J. Phys. Chem. 81:2340-2361 | SSA algorithm |
| 14 | Stormo GD, Zhao Y (2010) Nat. Rev. Genet. 11:751-760 | DNA-protein binding kinetics, k_on value |

---

## 14. VERIFICATION CHECKLIST

After the full build, verify these end-to-end criteria:

**Note:** Outcome-based checks listed below are scientific expectation checks unless explicitly labeled as software/sanity checks. Unexpected scientific outcomes should generate warnings, not build-stopping failures.

### Phase 1 Checks
- [ ] H37Rv genome downloaded, ~4.4 Mb, GC ≈ 65.6%
- [ ] 3,500-4,500 upstream regions extracted (no sequences <200bp)
- [ ] MEME finds ≥1 motif with E-value < 0.05 (or MOCK flagged)
- [ ] FIMO recovers known operator in top 5 hits (or MOCK flagged)
- [ ] Known regulon genes (Rv1933c-Rv1935c, Rv1936-Rv1941) recovered
- [ ] Conservation check: known operator conserved in M. bovis + M. marinum
- [ ] predicted_sites.csv exists with ranked predictions

### Phase 2 Checks
- [ ] All 5 conditions complete (A: 10K, B: 10K, C: 10K, D: 10K, E: 10×1K cells)
- [ ] CV(Asymmetric) > CV(Symmetric)
- [ ] CV(Asymmetric) > CV(Single-site)
- [ ] Condition C (single-site): GMM selects 2 components (bimodal, not trimodal)
- [ ] Fano(Unregulated) ≈ 7.85 (within [5, 12])
- [ ] Asymmetry sweep shows monotonic increase in intermediate fraction
- [ ] Mce3R ratio ≈ 20.4 marked on sweep plot
- [ ] Trace data saved for 20 cells per condition (A, B, C)
- [ ] No NaN/Inf values in any output
- [ ] Numba JIT compilation succeeds without error

### Phase 3 Checks
- [ ] Bootstrap CIs: lower < point < upper for all metrics
- [ ] CIs for CV(A) and CV(B) do NOT overlap
- [ ] KS test A vs B: p < 0.05
- [ ] KS test A vs C: p < 0.05
- [ ] Sensitivity: multimodality persists across >50% of parameter range for most parameters
- [ ] Tornado plot identifies top parameter drivers
- [ ] Shuffle test: ≥8/10 top sites pass
- [ ] Symmetric TetR control: unimodal
- [ ] Poisson baseline: Fano in [5, 12]
- [ ] Predicted persister fraction: 0.01%-10% (both GMM and 2σ methods)
- [ ] All p-values reported honestly (including non-significant ones)
- [ ] GMM fitting handles discrete data without crashes

### Phase 4 Checks
- [ ] 7 figures generated as 300 dpi PNG
- [ ] All figures have axis labels, legends, consistent style
- [ ] Colors: teal (asymmetric), orange (symmetric), purple (single-site), gray (unregulated)
- [ ] Figure 3 has Mce3R dashed line annotation at ratio ≈ 20.4
- [ ] Figure 6 shows representative single-cell traces with visible switching
- [ ] Figure 7 shows BIC comparison with star on best model
- [ ] No blank/empty figures

### Integration Checks
- [ ] `python main.py` runs end-to-end without error
- [ ] Phase 1 and Phase 2 run in parallel (subprocess)
- [ ] BUILD_REPORT.md generated with all phases PASS
- [ ] BUILD_REPORT.md flags MOCK if MEME Suite was unavailable
- [ ] All output files exist in results/

---

*Architecture document for sequential build + parallel runtime execution.*
*Project: Mce3R Asymmetric Operator → Noise → Persistence*
*PI: Dr. Gábor Balázsi, Stony Brook University*
*Student Researcher: Aayan Alwani*
