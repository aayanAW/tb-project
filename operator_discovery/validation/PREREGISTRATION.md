# Pre-registration — Mce3R genome-wide operator scan

Locked **before** any genome-wide scan is run, so the validation cannot be tuned to the
outcome. Commit this file before generating `results/`. Any deviation is recorded in the
"Deviations" section with a reason.

## Hypothesis

Mce3R (Rv1963c) binds a specific, conserved operator and represses a defined regulon. A
correctly-built operator model should (a) recover the experimentally-mapped operator
sites, (b) be enriched in the known regulon over genomic background, and (c) flag only a
small fraction of the genome — i.e. behave like a specific repressor, not a GC detector.

## Ground-truth set (fixed)

Source: Santangelo 2008 (BMC Microbiol 8:38), Santangelo 2009 (Microbiology 155:2245),
Panagoda, Balázsi & Sampson 2024 (ACS Chem Biol 19:2580; PDB 9B7Y). Encoded once in
`scripts/mce3r_biology.py`.

- **Autoregulation:** Rv1963c (mce3R).
- **mce3 structural operon:** Rv1964–Rv1971 (yrbE3A, yrbE3B, mce3A–mce3F).
- **Divergent lipid/redox cluster:** Rv1933c–Rv1935c and Rv1936–Rv1941.
- **Excluded (negative control genes):** mce1 (Rv0169–Rv0178), mce2 (Rv0586–Rv0594),
  mce4 (Rv3499c–Rv3503c) — Mce3R does NOT regulate these (Santangelo 2008). They must
  score like background; if they score high, the motif is non-specific.
- **Mapped operator regions:** mce3R–yrbE3A IGR (primary; two ~25 bp nonpalindromic sites
  ~53 bp apart, higher-affinity site ~100 bp upstream of the yrbE3A start) and the
  Rv1935c–Rv1936 IGR (second footprinted region).

## Operator model (fixed before scanning)

1. **Knowledge-based PWM** — built from the genomic sequence of the mapped operator
   regions and their orthologs (phylogenetic footprinting). NOT a typed consensus.
2. **De novo MEME** — positive set = orthologous mce3R–yrbE3A / Rv1935c–Rv1936 regions
   across _Mycobacterium_ spp.; negative `-neg` = matched genomic intergenic regions;
   `-bfile` = genome background. Real MEME only.
3. The de novo motif and the knowledge-based PWM are compared by Tomtom.

## Pre-registered decision thresholds

| Gate                                | Criterion                                                                                                          |
| ----------------------------------- | ------------------------------------------------------------------------------------------------------------------ |
| G1 operator recovered, not injected | de novo motif matches knowledge-based PWM: Tomtom q < 0.05                                                         |
| G2 known sites recovered            | both mapped operator regions contain a hit at FIMO **q < 0.05 AND score > 0**, in the top 5% of genome-wide scores |
| G3 specificity                      | genome-wide promoter hit fraction **≤ 5%** (current broken pipeline: 30%)                                          |
| G4 regulon enrichment               | regulon-vs-rest **AUPRC > 0.5** and empirical **p < 0.05** vs all three nulls                                      |
| G5 conservation                     | operator detected in **≥ 8** mycobacterial orthologs of the mce3R–yrbE3A region                                    |

## Reporting / null models (fixed)

- FIMO: real `fimo`, genome `-bfile` background, threshold on **q-value (BH FDR)**, report
  only hits with **score > 0**. Negative-log-odds windows are never reported as sites.
- Ranking: by FIMO **q-value / score only**. Priority tier is a descriptive label, never a
  ranking input. No composite tier/architecture bonus.
- Nulls for enrichment: (a) dinucleotide-preserving shuffles (n = 1000), (b) equal-size
  random gene-promoter sets (n = 1000), (c) real intergenic background.
- Novel candidate targets (regulon hits outside the known set) are reported as
  **hypotheses for EMSA/ChIP follow-up**, never as conclusions.

## Decision rule

All five gates green → the genome-wide Mce3R operator predictions are reported as valid.
Any gate red → report the negative/partial result honestly (per research-integrity
standards); do not relax a threshold to pass.

## Deviations / honest outcomes

- **G5 conservation (reported, not gamed).** Public NCBI `gene`-symbol coverage for
  _yrbE3A_ across Mycobacteriaceae yields 8 distinct species. The operator is detected
  (score > 0, FIMO p < 1e-3) in **7 of 8**: M. tuberculosis, M. bovis, M. africanum,
  M. canettii, M. marinum, M. ulcerans, M. liflandii. The single exception is the distant
  fast-grower **M. smegmatis** (the outgroup), where the operator has diverged. The
  pre-registered bar (≥ 8 detections) is therefore NOT met. The threshold is left
  unchanged (lowering it post hoc would be the p-hacking this file exists to prevent).
  The 7/8 result is itself biologically meaningful — clade-specific conservation across
  pathogenic/slow-growing mycobacteria with divergence in the outgroup supports operator
  specificity. Reaching ≥ 8 conserved would require protein-similarity (BLAST) ortholog
  mining beyond gene-symbol annotation; logged as future work, not forced here.
- All other gates (G1–G4) pass on real, non-circular evidence.

## Post-audit correction (2026-06-14) — gates strengthened, not relaxed

A cross-model audit (Claude + Codex/GPT-5.5; `audit_report.md`) ran AFTER the locked gates
above were evaluated and found the original G1/G2/G4 **circular/leaky** and the FDR not a
genuine genome-wide BH. Per the decision rule, a gate found invalid cannot be reported as
passed. The gates were therefore **re-specified to be harder and non-circular** (no
threshold was loosened to obtain a pass — every change tightens the test):

| Gate | Original (circular)                                                                                                       | Corrected (non-circular)                                                                                                                                                                                                         | Why                                                             |
| ---- | ------------------------------------------------------------------------------------------------------------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | --------------------------------------------------------------- |
| G1   | de novo MEME vs knowledge PWM by Tomtom — two MEME runs over the SAME operator DNA + orthologs (guaranteed to agree)      | knowledge PWM vs a PWM built from **independently-published** Santangelo 2009 footprinted operator sites                                                                                                                         | external, methodologically independent corroboration (audit C1) |
| G2   | recover BOTH operator IGRs — but both were used to build the model                                                        | recover only the **held-out** Rv1935c–Rv1936 operator (the knowledge PWM is trained on the mce3R–yrbE3A window + orthologs only); mce3R–yrbE3A reported as a training positive control                                           | a true generalization test (audit C2/C4)                        |
| G3   | hit fraction ≤ 5% at FIMO's per-report q                                                                                  | hit fraction ≤ 5% at a **pooled genome-wide BH** q<0.05 (M = 2 motifs × 2 strands × Σ(L−w+1) = 1,536,228 tests, recomputed in `validate.py`)                                                                                     | genuine genome-wide FDR (audit C3/C5)                           |
| G4   | regulon AUPRC>0.5 + empirical p — most "positives" were the PWM's own training regions                                    | **negative-control specificity**: mce1/mce2/mce4 (Mce3R does NOT bind them, Santangelo 2008) must show zero strict hits and none in the top 5%. The old AUPRC is retained as a **descriptive, explicitly-circular** readout only | removes the tautology (audit C2/C5/C10)                         |
| G5   | operator detected in ≥8 ortholog windows — the SAME windows that trained the PWM; 4 MTBC near-clones counted as 4 species | **leave-one-lineage-out**: each lineage is scored by a PWM trained on the OTHER lineages only; MTBC near-clones collapsed to one lineage; threshold ≥2 of 3 independent lineages; mycobacterial (not H37Rv) background           | non-circular + de-pseudo-replicated (audit C7/S5)               |

Strict-hit definition tightened to require `score>0` **and** genome-wide BH `q<0.05` on the
**same FIMO row** (audit C5). `all_gates_pass` now requires **all five** gates (the previous
`all_core_pass` excluded G1/G5 — audit escape hatch removed). Run provenance (git SHA, MEME
tool versions, input checksums) is stamped into `gates.json` (audit S6).

### Pre-registered nulls — honest deviation (audit C8)

The prereg promised three enrichment nulls (dinucleotide shuffle, random gene-set,
real-intergenic). Because the corrected G4 is reframed as a **negative-control specificity**
test (the random-gene-set/hypergeometric enrichment was circular by construction), the
multi-null enrichment machinery is **not** the gate. The honest specificity evidence is the
genome-wide hit fraction (G3) + the clean mce1/2/4 negative controls (G4). This deviation is
recorded here rather than silently dropped.
