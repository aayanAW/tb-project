# Mce3R Binding Motif Discovery Pipeline

A Python bioinformatics pipeline for identifying DNA binding motifs of **Mce3R**, a TetR-family transcription factor in *Mycobacterium tuberculosis*.

## Biological Background

Mce3R represses the **mce3 operon** (Rv1970–Rv1978), which encodes mammalian cell entry (Mce) proteins critical for *M. tuberculosis* virulence and host-pathogen interactions. Mce3R belongs to the **TetR transcription factor family**, which bind DNA as homodimers and typically recognize **palindromic inverted repeat sequences** (~15–20 bp) in target promoters.

This pipeline identifies and characterizes the Mce3R operator sequence by:
- Discovering conserved motifs in promoter/intergenic regions (MEME)
- Scanning for motif occurrences genome-wide (FIMO)
- Analyzing palindromic symmetry, inter-site spacing, and motif architecture

## Project Structure

```
mce3r_project/
├── data/
│   ├── raw/                    # Input FASTA sequences
│   │   └── sequences.fasta
│   └── processed/              # Analysis outputs
│       ├── sequence_metadata.csv
│       ├── annotated_hits.csv
│       └── spacing_analysis.csv
├── scripts/
│   ├── utils.py                # Shared helpers
│   ├── generate_sequences.py   # Synthetic sequence generator
│   ├── run_meme.py             # MEME wrapper + simulation fallback
│   ├── run_fimo.py             # FIMO wrapper + Python scanner fallback
│   ├── analyze_results.py      # Statistical analysis
│   └── visualize.py            # Figure generation
├── results/
│   ├── motifs/                 # MEME output (meme.txt, etc.)
│   ├── scans/                  # FIMO output (fimo.tsv, etc.)
│   └── figures/                # PNG and SVG plots
├── main.py                     # Pipeline orchestrator
├── requirements.txt
└── README.md
```

## Installation

### 1. Python dependencies

```bash
pip install -r requirements.txt
```

Required: `biopython`, `pandas`, `numpy`, `matplotlib`, `seaborn`

### 2. MEME Suite (optional)

The pipeline works **without MEME Suite** using built-in Python simulation fallbacks. For production use with real data, install MEME Suite:

```bash
# Recommended (conda)
conda install -c bioconda meme

# Or download from:
# https://meme-suite.org/meme/doc/install.html
```

## Quick Start

Run the entire pipeline (simulation mode if MEME Suite not installed):

```bash
cd mce3r_project
python main.py
```

## Usage

### Run all steps
```bash
python main.py
```

### Run specific steps only
```bash
python main.py --steps generate meme fimo
python main.py --steps analyze visualize
```

### Force simulation mode (skip real MEME/FIMO)
```bash
python main.py --simulate
```

### Skip steps with existing output
```bash
python main.py --skip-existing
```

### Custom parameters
```bash
python main.py --n-sequences 50 --n-with-motif 30 --seed 123 --pvalue-threshold 1e-4
```

### All options
```
--steps          Pipeline steps to run (default: all)
--seed           Random seed for reproducibility (default: 42)
--n-sequences    Number of synthetic sequences (default: 25)
--n-with-motif   Sequences with embedded motif (default: 14)
--threads        CPU threads for MEME (default: 4)
--pvalue-threshold  FIMO p-value cutoff (default: 1e-3)
--simulate       Force simulation mode
--skip-existing  Skip steps with existing output
```

## Pipeline Steps

### 1. Generate Sequences (`generate_sequences.py`)
Creates synthetic *M. tuberculosis*-like promoter sequences:
- **25 sequences**, 200–300 bp each
- **65% GC content** (matching M. tb H37Rv genome)
- **First-order Markov chain** with dinucleotide bias for realism
- **14 sequences** contain an embedded Mce3R motif (`TTGACANNNNNTGTCAA`)
- FASTA headers encode ground-truth metadata for benchmarking

### 2. MEME Motif Discovery (`run_meme.py`)
```
meme sequences.fasta -dna -mod zoops -nmotifs 5 -minw 6 -maxw 20 -oc results/motifs/
```
- `-mod zoops`: Zero or one occurrence per sequence (appropriate for promoter analysis)
- **Simulation fallback**: writes synthetic `meme.txt` with position weight matrix

### 3. FIMO Motif Scanning (`run_fimo.py`)
```
fimo --oc results/scans/ results/motifs/meme.txt data/raw/sequences.fasta
```
- Scans both strands
- **Simulation fallback**: Python log-odds scorer with empirical p-values
- Writes canonical `fimo.tsv` format regardless of mode

### 4. Analysis (`analyze_results.py`)
- **Inter-site spacing**: gap between adjacent binding sites in same sequence
- **Palindrome detection**: scores each site for half-site symmetry (TetR homodimer model)
- **Direct repeat detection**: distinguishes palindromes from direct repeats
- **Performance metrics**: precision/recall/F1 against synthetic ground truth

### 5. Visualization (`visualize.py`)
| Figure | Description |
|--------|-------------|
| `score_distribution` | FIMO log-odds scores colored by strand |
| `position_distribution` | Motif positions along each sequence |
| `spacing_histogram` | Gap between adjacent sites (highlights 30–50 bp range) |
| `palindrome_scores` | Symmetry scores for all predicted sites |
| `gc_vs_hits` | Sequence GC content vs. hit count (bias check) |
| `sequence_logo` | Stacked bar chart approximating sequence logo |

All figures saved as `.png` (300 DPI) and `.svg`.

## Output Files

| File | Description |
|------|-------------|
| `data/raw/sequences.fasta` | Input sequences (generated or user-provided) |
| `data/processed/sequence_metadata.csv` | Per-sequence metadata (GC, length, ground truth) |
| `results/motifs/meme.txt` | MEME output: position weight matrix |
| `results/scans/fimo.tsv` | FIMO output: binding site hit table |
| `data/processed/annotated_hits.csv` | FIMO hits joined with sequence metadata |
| `data/processed/spacing_analysis.csv` | Inter-site spacing for multi-hit sequences |
| `results/scans/summary_statistics.json` | Full summary statistics + performance metrics |
| `results/figures/*.png` | Publication-quality figures |

## Simulation Mode

When MEME Suite is not installed, the pipeline automatically uses simulation mode:

- **MEME simulation**: Writes `meme.txt` with a synthetic Position Weight Matrix derived from the consensus `TTGACANNNNNTGTCAA` with realistic Dirichlet noise
- **FIMO simulation**: Python log-odds scanner that:
  1. Parses the PWM from `meme.txt`
  2. Converts to log-odds scores using M. tb background frequencies
  3. Scans both strands with a sliding window
  4. Computes empirical p-values from 10,000 background windows
  5. Writes output in canonical FIMO TSV format

Simulation mode produces **scientifically valid output** — not a stub. The scanner will find hits preferentially in sequences where the motif was embedded.

## Using Real Sequence Data

To replace synthetic sequences with real *M. tuberculosis* promoter sequences:

1. Download promoter regions from NCBI (genome accession NC_000962.3 for H37Rv)
2. Format as FASTA with one sequence per promoter/intergenic region
3. Place in `data/raw/sequences.fasta`
4. Skip the generate step:

```bash
python main.py --steps meme fimo analyze visualize
```

### Recommended sequence sources for Mce3R:
- Upstream regions (−300 to +1) of Mce3R-regulated genes
- ChIP-seq peak regions (if available)
- Known TetR operator sequences from related *Mycobacterium* species

## Tiered Prioritization of Predicted Binding Sites

Every predicted binding site receives a `priority_tier` label that reflects the biological confidence of that prediction. Sites are sorted tier-first (Tier 1 → Tier 2 → Tier 3), then by motif score within each tier.

### Tier 1 — Known Mce3R Regulatory Region (Validation Set)

These are the genes whose promoters are **experimentally established** as Mce3R targets:

| Locus | Gene | Role |
|-------|------|------|
| Rv1963c | mce3R | The repressor itself (autoregulation) |
| Rv1964–Rv1965 | yrbE3A, yrbE3B | Inner membrane permeases of mce3 transporter |
| Rv1966–Rv1972 | mce3A–mce3F | Core mce3 operon structural genes |
| Rv1973–Rv1978 | Accessory genes | Additional mce3 operon members |

A hit in a Tier 1 promoter is **self-validating** — it constitutes direct computational evidence for an Mce3R operator.

### Tier 2 — Biologically Plausible Targets

Genes that share metabolic context with Mce3R regulation or belong to related Mce-family systems:

- **Mce1 operon** (Rv0169–Rv0178): Fatty acid uptake, functionally analogous to mce3
- **Mce2 operon** (Rv0586–Rv0594): Related lipid import system
- **Mce4 operon** (Rv3499c–Rv3513c): Cholesterol import, critical for intracellular survival
- **Cholesterol catabolism** (kshA/Rv3526, kshB/Rv3577, kstD/Rv3537, hsa genes): Downstream of Mce4 cholesterol import
- **igr operon** (Rv3543c–Rv3548c): Required for intracellular growth on cholesterol
- **Oxidative stress response** (katG/Rv1908c, ahpC/Rv2428, sodA/Rv3846, sodC/Rv0432, tpx/Rv1932)
- Any gene with lipid/fatty acid/cholesterol-related product annotation

Tier 2 hits are **strong candidates for experimental validation**.

### Tier 3 — Exploratory (Genome-Wide)

All other predicted sites across the M. tuberculosis genome. These may represent:
- Secondary regulatory targets not yet characterized for Mce3R
- False positives from the motif scan
- Genes regulated by Mce3R in specific conditions only

Tier 3 hits require **independent experimental validation** before biological conclusions can be drawn.

### Output Column: `priority_tier`

The `predicted_mce3r_sites.csv` table includes a `priority_tier` column with one of three values:

| Value | Meaning |
|-------|---------|
| `Tier1_known_operator` | mce3 operon — highest confidence |
| `Tier2_plausible_target` | Lipid/cholesterol/stress genes — strong candidates |
| `Tier3_exploratory` | All other genome-wide predictions |

Example:
```
sequence_name,start,stop,strand,matched_sequence,score,priority_tier
Rv1963c,45,61,+,TTGACAGTGTCTATCAA,14.25,Tier1_known_operator
Rv3526,88,104,-,TTGACATCGGGTGTCAA,13.71,Tier2_plausible_target
Rv0118c,98,114,-,TTTACGCGGTCTATCAA,13.20,Tier3_exploratory
```

## Interpreting Results

### FIMO TSV (`fimo.tsv`)
Each row is a predicted Mce3R binding site:
- `score`: Log-odds score (higher = better motif match)
- `p-value`: Probability of seeing this score by chance
- `strand`: Which strand the site is on (+/−)
- `matched_sequence`: The actual sequence at this site

### Palindrome Score
- Range: 0.0 (fully asymmetric) to 1.0 (perfect palindrome)
- **Threshold ≥ 0.8**: Consistent with TetR homodimer binding
- TetR operators: left half-site is reverse complement of right half-site

### Spacing Histogram
- Peak at 30–50 bp: Suggests tandem operator pairs
- Mce3R may use cooperative binding (two dimers at tandem operators)

### Performance Metrics (synthetic data)
- **Precision**: Fraction of predicted positive sequences that truly contain the motif
- **Recall**: Fraction of motif-containing sequences that FIMO detected
- **F1**: Harmonic mean of precision and recall

## Requirements

- Python ≥ 3.8
- biopython ≥ 1.81
- pandas ≥ 2.0
- numpy ≥ 1.24
- matplotlib ≥ 3.7
- seaborn ≥ 0.12
- MEME Suite ≥ 5.0 (optional)
