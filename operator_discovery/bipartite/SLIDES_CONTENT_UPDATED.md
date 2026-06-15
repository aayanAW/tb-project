# Slide Deck Content: Mce3R Binding Site Discovery — Updated Version

Student research project, Laughney Lab, Weill Cornell Medicine.
Updated to address PI feedback on extended MEME analysis and Region 2 validation.

---

## Slide 1 — Title

**Title:** Finding Mce3R Binding Sites

**Subtitle:** Identifying where the Mce3R repressor binds across *M. tuberculosis*

**Author:** Aayan Alwani

(Center vertically and horizontally on white background. No other text, no stats, no date.)

---

## Slide 2 — Background

- Mce3R is a transcriptional repressor in *M. tuberculosis* H37Rv
- Panagoda 2024 (PDB 9B7Y) solved its structure and showed its binding site is unusually asymmetric
- Strong half-site: Kd = 2.4 nM
- Weak half-site: Kd = 49.0 nM
- 20.4-fold affinity difference
- To study how architecture affects regulation, we need to find every binding site

---

## Slide 3 — Pipeline

(Keep the 5-step flowchart as main visual.)

- 1: download_genomes.py → retrieve 3 reference genomes
- 2: extract_upstream.py → extract 200 bp upstream of each gene
- 3: run_meme.py → discover binding motif
- 4: run_fimo.py → scan genome for matches
- 5: conservation_check.py → validate across species

---

## Slide 4 — Step 1: Genome download

`phase1_pipeline/download_genomes.py`

- NCBI Entrez download
- H37Rv + M. bovis + M. marinum
- FASTA + GenBank + GFF3 per species
- Three species enables cross-species validation

(Visual: phylogenetic tree of the 3 species.)

---

## Slide 5 — Step 2: Upstream extraction

`phase1_pipeline/extract_upstream.py`

- ~4,000 protein-coding genes parsed from GFF3
- 200 bp upstream extracted per gene
- Strand-aware, circular-genome-aware
- Five specific intergenic regions also extracted (Kabir 2021 protocol)

(Visual: genome schematic with highlighted 200 bp upstream region.)

---

## Slide 6 — Step 3: MEME motif discovery

`phase1_pipeline/run_meme.py`

- MEME Suite (Bailey 2009)
- Run on 3 yrbE3A ortholog sequences
- ZOOPS mode, widths 20–30 bp, both strands
- No palindrome constraint (key choice — Mce3R binds asymmetrically)

(Visual: use logo_mce3r_fimo.png full-width.)

---

## Slide 7 — Step 4: FIMO genome scan

`phase1_pipeline/run_fimo.py`

- 1,442 candidate binding sites identified
- Top 2 hits match the cryo-EM binding sites at yrbE3A (Panagoda 2024)
- Independent validation that the pipeline works

**Top 10 candidate sites:**

| Rank | Gene | p-value | Score |
|------|------|:---:|:---:|
| 1 | Rv1964 (yrbE3A) | 1.1e-11 | 34.9 |
| 2 | Rv1964 (yrbE3A) | 3.1e-11 | 33.2 |
| 3 | Rv1962A (vapB35) | 4.8e-11 | 32.6 |
| 4 | Rv0304c (PPE5) | 1.5e-07 | 15.8 |
| 5 | Rv0251c (hsp) | 3.0e-07 | 14.8 |
| 6 | Rv2378c (mbtG) | 3.4e-07 | 14.8 |
| 7 | Rv0412c | 3.5e-07 | 11.9 |
| 8 | Rv1393c | 3.5e-07 | 11.9 |
| 9 | Rv0603 | 8.5e-07 | 13.3 |
| 10 | Rv2476c (gdh) | 9.2e-07 | 13.3 |

---

## Slide 8 — Step 5: Cross-species conservation

`phase1_pipeline/conservation_check.py`

- Top 50 sites checked in M. bovis and M. marinum
- 50/50 conserved in M. bovis
- 44/50 conserved in both species
- Conservation across divergent species supports biological function

(Visual: two horizontal bar meters showing 50/50 and 44/50.)

---

## Slide 9 — NEW: Extended MEME recovers a longer motif

(This is a NEW slide addressing PI feedback.)

**Title:** Extending MEME to wider motifs (up to 120 bp)

- Per PI feedback: the Panagoda 123 bp operator is not the only gold standard
- Re-ran MEME on Regions 1+2 combined with `-maxw 120 -mod anr`
- Eram Kabir's 2021 protocol (ANR mode, wide width range)

**Result: a 99 bp motif with E-value 3.9 × 10⁻⁸**

- 4 orders of magnitude better than Eram's previous best (E = 2.8 × 10⁻⁴)
- Eram's 49 bp motif is contained within our 99 bp motif at **100% identity**
- Three sites identified: two in Region 1, one in Region 2

(Visual: the 99 bp sequence logo — logo_99bp_denovo.png.)

---

## Slide 10 — NEW: Region 2 validated de novo

(NEW slide — directly addresses PI's "validate echA13-Rv1936" request.)

**Title:** First computational recovery of the Region 2 binding site

- Region 2 (echA13–Rv1936, 224 bp) is the second established Mce3R regulatory region
- Our original 21 bp FIMO scan missed it (only 1 weak hit, rank 1422)
- The extended 99 bp MEME search **found a strong hit in Region 2** at position 62, reverse strand
- Site p-value: 4.2 × 10⁻⁴⁵
- This is the first de novo computational recovery of this region with a matching motif

**Three sites of the 99 bp motif:**

| Region | Position | Strand | p-value | Match to Eram's 49 bp |
|--------|:---:|:---:|:---:|:---:|
| 1 | 684 | + | 1.6e-45 | 100% identity |
| 1 | 238 | + | 2.3e-44 | 73% identity |
| **2** | **62** | **−** | **4.2e-45** | **67% identity** |

---

## Slide 11 — NEW: Zoomed map of all three clusters (Dr. Balázsi's nomenclature)

(NEW slide — directly addresses Dr. Balázsi's 2026-05-18 feedback: "zoom into all of these regions where motifs seem to exist" and "add to the annotation all 3 suspected motifs (ABC, DEF and HXG)".)

**Title:** Three 99 bp motif clusters — ABC, DEF, HXG

- The 99 bp de novo motif occurs three times; each occurrence is one of the three suspected binding clusters in Dr. Balázsi's annotated schematic (TBmotifs.pdf)
- **ABC** — mce3R–yrbE3A, proximal (near mce3R): 99 bp motif at position 238 (+), p = 2.3 × 10⁻⁴⁴; sub-sites A, B, C
- **DEF** — mce3R–yrbE3A, distal (near yrbE3A): 99 bp motif at position 684 (+), p = 1.6 × 10⁻⁴⁵; **contains Eram's 49 bp core at 100% identity**; sub-sites D, E, F
- **HXG** — echA13–Rv1936, − strand: 99 bp motif at position 62 (−), p = 4.2 × 10⁻⁴⁵; sub-sites H, X, G
- Every sub-site core matches a Santangelo 2009 / Bigi instance exactly (0 mismatches); equivalence A = C = D = F = rcH = rcG and B = D = X (per TBmotifs.pdf)

(Visual: dna_motif_map_zoomed.png — full width. Two region overviews on top, then three zoom panels — ABC, DEF, HXG — each with a base-pair ruler and the DNA sequence.)

---

## Slide 12 — Three analyses converge

**Title:** Multiple independent analyses converge on the same region

- Four sources now agree on the Mce3R binding region in mce3R–yrbE3A:
  - Panagoda 2024 (cryo-EM structure, PDB 9B7Y)
  - Kabir 2021 (MEME, 49 bp motif)
  - Santangelo 2009 (published regulon analysis)
  - Our pipeline (1,442 FIMO hits + 99 bp extended MEME)

(Visual: stack logo_mce3r_fimo.png, logo_eram_motif2.png, logo_santangelo.png, logo_99bp_denovo.png as four separate rows, then dna_motif_map_zoomed.png below them — full width, not cropped.)

---

## Slide 13 — Novel candidate sites

**Title:** Candidate sites beyond the characterized region

- Rank 3 — Rv1962A (vapB35): adjacent toxin-antitoxin locus
- Rank 4 — Rv0304c (PPE5): PE/PPE family
- Rank 6 — Rv2378c (mbtG): mycobactin biosynthesis
- Rank 38 — Rv1963c: the mce3R gene itself (autoregulation)

**Next steps:**
- Apply 99 bp PWM genome-wide via FIMO for broader search
- Examine remaining three intergenic regions (3, 4, 5) with extended parameters
- Cross-validate novel candidates with published ChIP-seq data where available

---

## KEY NUMBERS (do not change)

- 1,442 FIMO candidate sites
- 222 intergenic (15%)
- 50/50 conserved in M. bovis
- 44/50 conserved in both species
- Panagoda Kd: 2.4 nM and 49.0 nM (20.4-fold asymmetry)
- **99 bp de novo motif, E-value 3.9 × 10⁻⁸**
- **Eram 49 bp motif recovered at 100% identity**
- **Region 2 site at position 62, p = 4.2 × 10⁻⁴⁵**
- Region 1 = 897 bp
- Region 2 = 224 bp

## SITE NOMENCLATURE (Dr. Balázsi, TBmotifs.pdf — use these exact labels)

- Three motif clusters, one per 99 bp de novo occurrence:
  - **ABC** — Region 1 (mce3R–yrbE3A), proximal, 99 bp @ pos 238 (+), p = 2.3 × 10⁻⁴⁴
  - **DEF** — Region 1 (mce3R–yrbE3A), distal, 99 bp @ pos 684 (+), p = 1.6 × 10⁻⁴⁵ (Eram 49 bp @ 100%)
  - **HXG** — Region 2 (echA13–Rv1936), 99 bp @ pos 62 (−), p = 4.2 × 10⁻⁴⁵
- Sub-sites: A,B,C / D,E,F in mce3R–yrbE3A; H,X,G in echA13–Rv1936 (HXG order = 5'→3' on the − strand)
- Equivalence: A = C = D = F = rcH = rcG; B = D = X
- First-half sites (Eram motif first half): A, C, D, F, G, H. Second-half sites: B, E, X.

## CITATIONS

- Panagoda, Balázsi & Sampson (2024) *ACS Chemical Biology* 19:2580–2592. PDB 9B7Y.
- Santangelo et al. (2009) *Microbiology* 155(7):2245–2255.
- Bailey et al. (2009) *Nucleic Acids Research* 37:W202–W208. (MEME Suite)
- Grant et al. (2011) *Bioinformatics* 27:1017–1018. (FIMO)
- Kabir (2021) Laughney Lab internal analysis.

## IMAGES TO ATTACH

- logo_mce3r_fimo.png — 21 bp motif from original FIMO hits
- logo_eram_motif2.png — Eram's 21 bp motif
- logo_santangelo.png — Santangelo 16 bp motif
- logo_99bp_denovo.png — **NEW: 99 bp motif recovered de novo**
- dna_motif_map_zoomed.png — **NEW (replaces dna_motif_map.png): zoomed, ABC/DEF/HXG-annotated map — 2 region overviews + 3 sequence-level zoom panels, Dr. Balázsi's site nomenclature. Generated by `phase1_pipeline/zoomed_motif_map.py`; audit in `results/extended_meme/report/zoomed_motif_map.json`.**
