# Mce3R Genome-Wide Operator Scan (corrected, validated)

A bioinformatics pipeline that predicts and validates genome-wide binding sites of
**Mce3R**, the TetR-family transcriptional repressor of _Mycobacterium tuberculosis_
(Rv1963c). The operator model is built from **real operator-region sequence** and the
genome-wide scan is validated against the **known regulon** with pre-registered gates.

> **History.** An earlier version of this project produced circular results: it ran MEME on
> 9 positive-only promoters, scanned for a hand-typed operator consensus via a
> self-referential Python "FIMO simulation" (p-value floor == reporting threshold, q == p,
> no FDR), and flagged ~30% of the genome ranked by a hard-coded known-gene prior. That
> code is removed; its outputs are archived under `results/_deprecated_circular/`. This
> README documents the corrected pipeline.

## Biology (single source of truth: `scripts/mce3r_biology.py`)

- **Mce3R (Rv1963c)** is a TetR-family / "double-TFR" repressor. Its operator is
  **nonpalindromic**: two ~25 bp asymmetric sites ~53 bp apart in the mce3R–yrbE3A
  promoter (Panagoda, Balázsi & Sampson 2024, _ACS Chem Biol_ 19:2580; PDB 9B7Y; Kd 2.4 nM).
- **Regulon:** mce3 structural operon (Rv1964–Rv1971) + the divergent lipid/redox cluster
  Rv1933c–Rv1935c and Rv1936–Rv1941 + autoregulation (Santangelo 2009, _Microbiology_
  155:2245). Mce3R does **not** regulate mce1/mce2/mce4 (Santangelo 2008, _BMC Microbiol_
  8:38) — these are used as negative controls.

## Requirements

- Python ≥ 3.10 with `biopython`, `pandas`, `numpy` (`scikit-learn` optional; a local
  AUPRC/AUROC fallback is built in).
- **MEME Suite** (`meme`, `fimo`, `tomtom`, `fasta-get-markov`) — required, no simulation
  fallback. On Apple Silicon:
  ```bash
  CONDA_SUBDIR=osx-64 conda create -y -n meme_x64 -c conda-forge -c bioconda meme
  ```
  The wrappers locate this env automatically (or set `MEME_BIN`).

## Run

```bash
python main.py                     # full pipeline + gates
python main.py --steps scan validate rank
python scripts/fetch_orthologs.py  # refresh mycobacterial yrbE3A orthologs from NCBI
```

## Pipeline steps

| Step             | What it does                                                                                                                       |
| ---------------- | ---------------------------------------------------------------------------------------------------------------------------------- |
| `ensure-genome`  | Extract the real H37Rv FASTA (NC_000962.3) from the GenBank file                                                                   |
| `extract`        | Genome-wide promoters + divergent IGRs; always includes both mapped operator IGRs                                                  |
| `background`     | Order-2 genomic Markov background (`fasta-get-markov`) for MEME/FIMO                                                               |
| `operator-model` | **Knowledge-based** operator PWM from real operator-region sequence (M. tb IGR + orthologs) via anchored MEME — no typed consensus |
| `denovo`         | **De novo** MEME on the bound IGRs vs a genomic control set, then Tomtom vs the knowledge PWM (Gate G1)                            |
| `scan`           | Real FIMO genome-wide with genomic background; q-values (BH FDR); score > 0 enforced                                               |
| `conservation`   | FIMO the operator across mycobacterial orthologs (Gate G5)                                                                         |
| `validate`       | Gates G2–G4 against the known regulon                                                                                              |
| `rank`           | Ranked predicted-site table by q-value/score (tier is a descriptive label only)                                                    |

## Validation gates (pre-registered — see `PREREGISTRATION.md`)

| Gate                                | Criterion                                                          |
| ----------------------------------- | ------------------------------------------------------------------ |
| G1 operator recovered, not injected | de novo motif matches knowledge PWM (Tomtom q < 0.05)              |
| G2 known sites recovered            | both operator IGRs hit at q < 0.05 & score > 0, top 5%             |
| G3 specificity                      | genome-wide promoter hit fraction ≤ 5%                             |
| G4 regulon enrichment               | AUPRC > 0.5 and empirical p < 0.05; negative controls ≈ background |
| G5 conservation                     | operator detected in ≥ 8 mycobacterial orthologs                   |

Ranking uses FIMO q-value/score **only**; the regulon tier is never a ranking input.
Negative-log-odds (score ≤ 0) windows are never reported as sites. Any gate failure is
reported honestly (see `PREREGISTRATION.md` → Deviations) rather than tuned away.

## Key outputs

- `results/motifs/operator_knowledge/operator_knowledge.meme` — operator PWM (real sequence)
- `results/motifs/denovo/` — de novo motif + Tomtom cross-check
- `results/scans/fimo.tsv` — genome-wide hits (genomic background, q-values)
- `results/scans/validation_report.txt`, `results/scans/gates.json` — gate results
- `data/processed/predicted_mce3r_sites.csv` — ranked sites (score > 0), tier as a label

Predicted sites outside the known regulon are **hypotheses for experimental validation**
(EMSA/ChIP), not conclusions.
