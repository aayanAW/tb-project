# Slide Deck Content: Transcription Factor Binding Site Discovery

Student research project, Laughney Lab, Weill Cornell Medicine.
Standard bioinformatics pipeline for identifying DNA binding motifs.

All data below is from:
- Panagoda, Balázsi & Sampson (2024) ACS Chemical Biology (peer-reviewed, PDB 9B7Y)
- Santangelo et al. (2009) Microbiology (peer-reviewed)
- Eram Kabir (2021) lab rotation analysis

---

## Slide 1 — Title

**Title:** Computational Identification of Transcription Factor Binding Sites

**Subtitle:** A standard bioinformatics pipeline applied to a mycobacterial reference genome

**Author:** Aayan Alwani, Laughney Lab

**Stats:**
- 1,442 candidate sites identified
- 5 upstream regions analyzed
- 3 reference genomes compared
- 50/50 top sites cross-species conserved

---

## Slide 2 — Background

A transcriptional repressor called Mce3R in *Mycobacterium tuberculosis* H37Rv has an unusual DNA binding architecture, recently characterized by cryo-EM (Panagoda 2024, PDB 9B7Y).

The binding site has two halves with different affinities:
- Half-site A: Kd = 2.4 nM
- Half-site B: Kd = 49.0 nM
- 20.4-fold difference in binding strength

To study how this architecture affects gene regulation, we needed to identify all potential binding sites across the genome using standard computational tools.

---

## Slide 3 — Pipeline

5 Python scripts in sequence:

| Step | Script | Task | Output |
|------|--------|------|--------|
| 1 | download_genomes.py | Retrieve reference genomes from NCBI | FASTA + GenBank |
| 2 | extract_upstream.py | Extract 200bp upstream of coding sequences | Upstream FASTA |
| 3 | run_meme.py | De novo motif discovery | Position Weight Matrix |
| 4 | run_fimo.py | Scan genome for motif matches | Predicted sites |
| 5 | conservation_check.py | Cross-species conservation check | Conservation table |

---

## Slide 4 — Step 1: Genome Download

**Script:** download_genomes.py

Downloads three publicly available reference genomes from the NCBI Entrez database:
- *M. tuberculosis* H37Rv (NC_000962.3)
- *M. bovis* AF2122/97 (NC_002945.4)
- *M. marinum* M (NC_010612.1)

Using three species allows cross-validation: if a predicted site is conserved across species that diverged millions of years ago, it is likely to be functional rather than a false positive.

---

## Slide 5 — Step 2: Upstream Region Extraction

**Script:** extract_upstream.py

Parses the GFF3 annotation of H37Rv to identify all protein-coding genes. For each gene, extracts the 200 bp region upstream of the start codon, where regulatory elements are typically located.

Handles:
- Circular genome wrapping at the origin
- Strand-aware extraction (reverse-complements minus-strand genes)

Also extracts five specific intergenic regions from a prior lab analysis (Kabir 2021):
- Region 1: mce3R–yrbE3A (897 bp)
- Region 2: echA13–Rv1936 (224 bp)
- Region 3: Rv1941–Rv1942c (209 bp)
- Region 4: Mce3F–Rv1972 (1,335 bp)
- Region 5: Rv1944c–Rv1945 (54 bp)

---

## Slide 6 — Step 3: MEME Motif Discovery

**Script:** run_meme.py

Runs MEME (Bailey et al. 2009, *Nucleic Acids Research*) on the three orthologous upstream sequences to discover a consensus DNA pattern.

Settings:
- ZOOPS mode (zero or one occurrence per sequence)
- Motif width: 20–30 bp
- Both DNA strands searched
- Palindromic constraint NOT used (standard for asymmetric motifs)

Output: a Position Weight Matrix (PWM) — a standard probability representation of the motif.

**Visual:** Use logo_mce3r_fimo.png to show the resulting sequence logo.

---

## Slide 7 — Step 4: Genome Scan

**Script:** run_fimo.py

FIMO (Grant et al. 2011, *Bioinformatics*) scans the entire H37Rv genome against the PWM from Step 3. Reports every position with p < 1e-4.

**Results:**
- 1,442 candidate sites total
- 222 in intergenic regions
- Score range: −2.55 to 34.86
- p-value range: 1.1e-11 to 9.9e-05

**Top 10 candidate sites:**

| Rank | Nearest gene | p-value | Score | Intergenic |
|------|--------------|---------|-------|:---:|
| 1 | Rv1964 (yrbE3A) | 1.1e-11 | 34.9 | Yes |
| 2 | Rv1964 (yrbE3A) | 3.1e-11 | 33.2 | Yes |
| 3 | Rv1962A (vapB35) | 4.8e-11 | 32.6 | No |
| 4 | Rv0304c (PPE5) | 1.5e-07 | 15.8 | No |
| 5 | Rv0251c (hsp) | 3.0e-07 | 14.8 | No |
| 6 | Rv2378c (mbtG) | 3.4e-07 | 14.8 | No |
| 7 | Rv0412c | 3.5e-07 | 11.9 | Yes |
| 8 | Rv1393c | 3.5e-07 | 11.9 | No |
| 9 | Rv0603 | 8.5e-07 | 13.3 | Yes |
| 10 | Rv2476c (gdh) | 9.2e-07 | 13.3 | No |

The top two candidate sites match the two binding sites already characterized by cryo-EM in Panagoda 2024 — independent validation that the pipeline works.

---

## Slide 8 — Step 5: Cross-species Conservation

**Script:** conservation_check.py

Takes the top 50 candidate sites and checks whether each is conserved in *M. bovis* and *M. marinum*. Uses sliding-window sequence alignment with ≥80% identity threshold.

**Results:**
- 50 / 50 conserved in *M. bovis* (100%)
- 44 / 50 conserved in both species (88%)

Sites conserved across distantly related species are more likely to be functionally important.

---

## Slide 9 — Cross-validation with published analyses

Three independent analyses of the same genomic region converge on the same binding site predictions:

1. Our pipeline (2026)
2. Kabir (2021) lab rotation — MEME on five intergenic regions
3. Santangelo et al. (2009) — published regulon characterization
4. Panagoda et al. (2024) — cryo-EM structure (PDB 9B7Y)

In Region 1 (mce3R–yrbE3A, 897 bp), all three sources identify motifs in the same ~100 bp window:

| Source | Position | Length | Identity |
|--------|:---:|:---:|:---:|
| Panagoda structural binding region | 674 | 123 bp | 100% |
| Kabir MEME motif | 691 | 49 bp | 100% |
| Santangelo motifs 1–4 | 244, 320, 690, 766 | 16 bp each | 100% |

**Visual:** Display the three sequence logos (logo_mce3r_fimo.png, logo_eram_motif2.png, logo_santangelo.png) as three stacked rows, then dna_motif_map.png below them.

---

## Slide 10 — Region 2 Validation

A second intergenic region previously characterized by Santangelo 2009 (Region 2: echA13–Rv1936, 224 bp) was also examined.

Santangelo motif 4 (TATTGGCTATGGACAT) is present at position 61 of Region 2 with 100% identity — confirming the pipeline can detect this second regulatory region.

---

## Slide 11 — Candidate sites beyond the characterized region

The pipeline also identifies candidate sites at genes not included in the originally characterized regulon:

- **Rank 3 — Rv1962A (vapB35):** adjacent toxin-antitoxin locus
- **Rank 4 — Rv0304c (PPE5):** PE/PPE family protein
- **Rank 6 — Rv2378c (mbtG):** mycobactin biosynthesis
- **Rank 38 — Rv1963c:** the mce3R gene itself, suggesting autoregulation

**Next steps:**
- Rerun MEME with extended motif widths (up to 120 bp)
- Improve coverage of Region 2 with a region-specific PWM
- Expand the DNA annotation map to show all five regions

---

## Key numbers (do not change)

- 1,442 candidate sites
- 222 intergenic (15%)
- 50/50 conserved in M. bovis
- 44/50 conserved in both species
- Panagoda Kd: 2.4 nM and 49.0 nM
- 20.4-fold asymmetry
- Region 1 = 897 bp
- Region 2 = 224 bp

---

## Citations

- Panagoda, Balázsi & Sampson (2024) *ACS Chemical Biology* 19:2580–2592. PDB 9B7Y.
- Santangelo et al. (2009) *Microbiology* 155(7):2245–2255.
- Bailey et al. (2009) *Nucleic Acids Research* 37:W202–W208. (MEME Suite)
- Grant et al. (2011) *Bioinformatics* 27:1017–1018. (FIMO)
- Kabir (2021) Laughney Lab internal analysis.
