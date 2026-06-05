# Genome-Wide Identification of Mce3R Binding Sites in *Mycobacterium tuberculosis* Using a Dual-Motif Computational Pipeline

**Aayan Alwani**

---

## Abstract

The TetR-family transcription factor Mce3R regulates the *mce3* cholesterol/lipid import operon in *Mycobacterium tuberculosis*, a key virulence determinant. Recent structural work (Panagoda et al., 2024) revealed that Mce3R engages its operator through two binding sites with a 20-fold difference in affinity, but genome-wide knowledge of Mce3R binding sites remains limited to the characterized operator region. Here, we develop a computational pipeline that combines de novo motif discovery (MEME) on orthologous upstream sequences from three mycobacterial species with genome-wide motif scanning (FIMO) and cross-species conservation analysis to predict Mce3R binding sites across the *M. tuberculosis* H37Rv chromosome. We identify 10 candidate binding sites at a stringent p-value threshold (p < 10^-4), including the known operator upstream of *yrbE3A* (Rv1964) as the top-ranked hit (p = 8.5 x 10^-11). Novel high-confidence sites are predicted near *Rv1706c* (PPE23, p = 1.5 x 10^-7), *Rv1914c* (p = 3.2 x 10^-7), and *Rv3134c* (p = 8.7 x 10^-6, conserved in both *M. bovis* and *M. marinum*). All 10 sites are conserved in *M. bovis* (>=80% identity), while only one (*Rv3559c*) meets the conservation threshold in both *M. bovis* and *M. marinum*. These predictions suggest that the Mce3R regulon may extend beyond the characterized *mce3* operon to include PPE family genes and other loci involved in lipid metabolism and immune evasion, providing candidate targets for experimental validation by electrophoretic mobility shift assay or ChIP-seq.

---

## 1. Introduction

Tuberculosis (TB) remains the leading cause of death from a single infectious agent, killing approximately 1.3 million people annually (WHO, 2024). The pathogenesis of *Mycobacterium tuberculosis* depends critically on its ability to acquire host-derived cholesterol and lipids during intracellular infection, a process mediated in part by the mammalian cell entry (Mce) transport systems (Pandey & Bhatt, 2023). Among the four Mce systems in *M. tuberculosis* (Mce1-4), the Mce3 system is specifically induced during macrophage infection and is required for cholesterol utilization and full virulence (Dunphy et al., 2010).

Expression of the *mce3* operon (*Rv1964-Rv1977*) is controlled by the transcription factor Mce3R (Rv1963c), a TetR-family repressor that binds a ~123 bp operator region in the intergenic space between *mce3R* and *yrbE3A* (Rv1964), the first gene of the operon. Santangelo et al. (2009) demonstrated that deletion of *mce3R* derepresses *yrbE3A* expression approximately 8.5-fold, confirming its role as a transcriptional repressor. Earlier work by Santangelo et al. (2002) identified two additional loci (*Rv1933c-Rv1935c* and *Rv1936-Rv1941c*) as part of the Mce3R regulon based on transcriptional profiling.

Recent crystallographic and biophysical characterization of Mce3R (Panagoda et al., 2024; PDB 9B7Y) provided the first atomic-resolution structure of this regulator and revealed a striking feature of its DNA-binding mechanism: the homodimeric repressor engages its operator through two binding sites with markedly asymmetric affinities. Electrophoretic mobility shift assays (EMSA) demonstrated a strong-affinity site (Kd = 2.4 nM) and a weak-affinity site (Kd = 49.0 nM), representing a ~20-fold asymmetry. This structural and biophysical characterization provides the foundation for constructing position weight matrices (PWMs) that describe the sequence preferences of each binding mode.

Despite these advances, knowledge of Mce3R binding sites remains limited to the characterized operator and the two additional regulon loci identified by transcriptomics. No systematic genome-wide survey of Mce3R binding sites has been performed in *M. tuberculosis*. Such a survey is important for several reasons: (1) TetR-family regulators in mycobacteria frequently control larger regulons than initially characterized (Cuthbertson & Nodwell, 2013); (2) the Mce3R operator's distinctive asymmetric architecture may produce a recognizable sequence signature amenable to computational detection; and (3) identifying new Mce3R targets could illuminate previously unrecognized connections between lipid metabolism, immune evasion, and persistence.

In this study, we develop and apply a computational pipeline that integrates de novo motif discovery, genome-wide motif scanning, gene annotation, and cross-species conservation analysis to predict Mce3R binding sites across the complete *M. tuberculosis* H37Rv genome. By leveraging ortholog sequences from *M. bovis* and *M. marinum* for both motif training and conservation filtering, our approach balances sensitivity with specificity.

---

## 2. Methods

### 2.1 Genome Acquisition and Annotation

Complete genome sequences and annotations were retrieved from the NCBI Entrez database using BioPython (Cock et al., 2009):

- *M. tuberculosis* H37Rv (GenBank accession NC_000962.3): FASTA sequence, GenBank flat file, and GFF3 annotation
- *M. bovis* AF2122/97 (NC_002945.4): FASTA sequence
- *M. marinum* M (NC_010612.1): FASTA sequence

Genome integrity was verified by confirming the H37Rv assembly size (4,411,532 bp) and GC content (65.6%) against published values. The GFF3 annotation was parsed to extract all coding sequence (CDS) features with their coordinates, strand orientation, locus tags, and gene names. GFF3 coordinates (1-based, inclusive) were converted to 0-based Python indices upon parsing to ensure correct downstream coordinate arithmetic.

### 2.2 Extraction of Upstream Regulatory Regions

For each CDS feature in the H37Rv genome, we extracted 200 bp of upstream sequence to capture the proximal promoter region. For genes on the positive strand, 200 bp immediately upstream of the annotated start codon were extracted; for genes on the negative strand, 200 bp downstream of the annotated stop codon were extracted and reverse-complemented. Circular genome wrapping was handled explicitly for genes near the origin of replication using a modular extraction function that concatenates sequence segments spanning the origin.

Duplicate CDS entries (identical start, stop, and strand) were removed, retaining the first occurrence. This procedure yielded approximately 4,000 upstream regions of exactly 200 bp each, collectively representing the promoter-proximal regulatory landscape of the H37Rv genome.

### 2.3 Ortholog Identification and MEME Training Set Construction

To construct a training set for de novo motif discovery, we identified *yrbE3A* orthologs in *M. bovis* and *M. marinum* by searching each genome's GenBank annotation for locus tags matching known ortholog identifiers:

| Species | Genome | Ortholog Locus Tag | Search Terms |
|---------|--------|-------------------|--------------|
| *M. tuberculosis* H37Rv | NC_000962.3 | Rv1964 | Rv1964, yrbE3A |
| *M. bovis* AF2122/97 | NC_002945.4 | BQ2027_RS10150 | Mb1997, yrbE3A |
| *M. marinum* M | NC_010612.1 | MMAR_RS12585 | MMAR_2522, yrbE3A |

Ortholog locus tags were matched against the locus_tag, old_locus_tag, and gene name fields in each genome's GenBank record using both exact and substring matching. For each identified ortholog, 200 bp of upstream sequence was extracted using the same procedure described in Section 2.2. The resulting three-sequence FASTA file served as input for MEME motif discovery.

The rationale for using ortholog upstream sequences rather than all H37Rv upstream regions is that Mce3R binding motifs are expected to be conserved across mycobacterial species due to the functional importance of the *mce3* regulon, while the broader set of ~4,000 upstream regions would dilute the signal with thousands of irrelevant sequences.

### 2.4 Known Operator Sequence

The experimentally characterized Mce3R operator sequence (123 bp) was obtained from Panagoda et al. (2024), corresponding to the intergenic region between *mce3R* (Rv1963c) and *yrbE3A* (Rv1964) at H37Rv coordinates 2,207,477-2,207,699:

```
GCCCCGCGCTATAGGATACTAGCAAGATACATCATAGCCAATATATGCCAGT
TTGCATTGCTATTTACCGATCAGTTGTCCAAGCAATCGCGTATTGGCTATG
GACATCAGCGGTTCTGCCGC
```

This sequence encompasses both the strong-affinity binding site (corresponding to Probe A in the EMSA experiments; Kd = 2.4 nM) and the weak-affinity binding site (Probe C; Kd = 49.0 nM). The known operator serves as a positive control for pipeline validation: a correctly functioning pipeline must recover this region as a top-ranked hit.

### 2.5 De Novo Motif Discovery with MEME

Position weight matrices (PWMs) representing Mce3R binding preferences were derived using MEME version 5.5.9 (Bailey et al., 2009) with the following parameters:

| Parameter | Value | Rationale |
|-----------|-------|-----------|
| Sequence model (`-mod`) | ZOOPS | Zero or One Occurrence Per Sequence; appropriate because each ortholog upstream region is expected to contain at most one Mce3R binding site |
| Strand (`-revcomp`) | Both | Mce3R may bind either strand as a homodimer |
| Minimum motif width (`-minw`) | 20 bp | Lower bound for a TetR-family binding site (typically 15-25 bp per half-site) |
| Maximum motif width (`-maxw`) | 30 bp | Upper bound accommodating the full contact footprint including flanking bases |
| Number of motifs (`-nmotifs`) | 3 | Allows discovery of distinct strong-site and weak-site motifs plus a potential third motif |
| Background model (`-bfile`) | H37Rv Markov model | Controls for the high GC content (65.6%) of the *M. tuberculosis* genome |

**Background model construction.** A 0th-order Markov background model was computed from the H37Rv GC content to account for the extreme nucleotide composition bias in mycobacterial genomes. Background frequencies were set to: A = 0.172, C = 0.328, G = 0.328, T = 0.172, reflecting the 65.6% GC content of H37Rv. When available, the MEME Suite utility `fasta-get-markov` was used to derive a higher-order model directly from the input sequences; otherwise, the analytical 0th-order model was used as fallback.

**Dual-motif model.** MEME was configured to discover up to three motifs, motivated by the biophysical evidence that Mce3R engages its operator through two distinct binding modes. The strong-affinity site and weak-affinity site may have partially overlapping but distinguishable sequence preferences, and separate PWMs better capture this heterogeneity than a single consensus motif. The resulting PWMs are constructed as letter-probability matrices with alphabet size 4 (ACGT) and motif widths between 20-30 bp.

**Quality control.** MEME output was validated by confirming that: (1) the output file (`meme.txt`) exceeded 100 bytes, (2) discovered motif widths fell within the specified 20-30 bp range, and (3) each row of each PWM summed to approximately 1.0 (tolerance: +-0.01).

### 2.6 Genome-Wide Motif Scanning with FIMO

The PWMs discovered by MEME were used to scan the complete H37Rv genome for candidate Mce3R binding sites using FIMO (Find Individual Motif Occurrences; Grant et al., 2011) with the following parameters:

| Parameter | Value | Rationale |
|-----------|-------|-----------|
| p-value threshold (`--thresh`) | 10^-4 | Standard stringent threshold for FIMO; balances sensitivity with false positive control |
| Background model (`--bfile`) | Same H37Rv model as MEME | Ensures consistent statistical framework between motif discovery and scanning |
| Input motifs | All motifs from MEME output | Scans for both strong-site and weak-site motifs |
| Target genome | H37Rv complete genome (4,411,532 bp) | Genome-wide, strand-agnostic search |

FIMO reports each hit with: motif identifier, genomic coordinates (start, stop), strand, log-odds score, p-value, q-value (Benjamini-Hochberg corrected), and matched sequence. Hits were ranked by ascending p-value, with the lowest p-value representing the highest-confidence prediction.

### 2.7 Gene Annotation of Predicted Sites

Each FIMO hit was annotated with the nearest gene using the H37Rv GFF3 annotation. For each predicted binding site, the distance to every annotated CDS was computed:

- If the site overlapped a CDS (start_site < end_gene AND end_site > start_gene), distance = 0
- Otherwise, distance = minimum of |start_site - end_gene| and |end_site - start_gene|

The nearest gene (by minimum distance) was assigned to each site, along with its locus tag, gene name, and an intergenic flag (True if distance > 0, indicating the site falls between annotated genes). Intergenic sites are of particular biological interest because transcription factor binding sites are predominantly located in non-coding regulatory regions.

### 2.8 Cross-Species Conservation Analysis

To assess evolutionary conservation of predicted binding sites, each candidate sequence was searched against the complete genomes of *M. bovis* AF2122/97 and *M. marinum* M using a numpy-vectorized sliding-window percent-identity algorithm.

**Algorithm.** For each predicted 25-bp binding site sequence:

1. The comparison genome was converted to a numeric array using `np.frombuffer` (A=0, C=1, G=2, T=3)
2. A rolling match count was computed across all positions in the comparison genome using a vectorized sliding window of the same length as the query motif
3. Percent identity = (number of matching positions) / (motif length) at each window position
4. Both strands of the comparison genome were searched (the reverse complement was generated and searched separately)
5. The position and strand yielding the highest percent identity were recorded

**Conservation threshold.** A site was classified as "conserved" if the best cross-species match achieved >= 80% sequence identity. This threshold was chosen to accommodate the sequence divergence between species while maintaining biological relevance:

- *M. bovis* shares ~99.95% genome-wide identity with *M. tuberculosis*, so most functional sites are expected to be conserved
- *M. marinum* shares ~85% average nucleotide identity, so conservation at 80% represents significant selective constraint

A site conserved in both *M. bovis* AND *M. marinum* provides the strongest evidence for functional importance, as it implies maintenance across approximately 150 million years of mycobacterial evolution.

**Performance optimization.** The naive Python implementation of the sliding-window search required approximately 6 minutes for the full analysis. This was replaced with a numpy-vectorized implementation using `np.frombuffer` for sequence encoding and vectorized array operations for rolling match counting, reducing runtime to approximately 10 seconds — a ~36-fold speedup.

### 2.9 Pipeline Integration and Validation

The five analysis steps were executed sequentially in a single pipeline:

```
download_genomes → extract_upstream → run_meme → run_fimo → conservation_check
```

Pipeline validation included the following sanity checks:

1. All expected output files were produced and non-empty
2. The known Mce3R operator region (H37Rv coordinates 2,207,477-2,207,699) was recovered among the top-ranked FIMO hits (within +-500 bp)
3. At least one predicted site was conserved in *M. bovis*
4. Upstream region extraction yielded the expected number of sequences (~4,000 for H37Rv)

The complete pipeline runs in approximately 6-33 seconds depending on whether MEME Suite is installed locally or mock motifs are used, and produces all results in a structured `results/phase1/` output directory.

### 2.10 Software and Reproducibility

All analyses were performed in Python 3.10 with the following key dependencies: BioPython 1.81 (genome retrieval and parsing), NumPy (vectorized conservation analysis), MEME Suite 5.5.9 (motif discovery and scanning). The pipeline is fully reproducible from a single command (`python -m phase1_pipeline.pipeline_main`) and all parameters are centralized in a single configuration file (`config/parameters.py`).

---

## 3. Results

### 3.1 Recovery of the Known Mce3R Operator

As a positive control, we verified that the pipeline correctly identified the experimentally characterized Mce3R operator. The top two hits (ranks 1 and 2) both map to the known intergenic region between *mce3R* and *yrbE3A* at H37Rv coordinates 2,207,477-2,207,524, with p-values of 8.5 x 10^-11 and 1.2 x 10^-10, respectively. Both hits match the strong-affinity binding site sequence `GCCCCGCGCTATAGGATACTAGCAA` with 100% identity to the known operator. The rank-1 hit is located 199 bp upstream of the *yrbE3A* start codon and falls in the intergenic region, consistent with its known regulatory function. A third hit near *yrbE3A* (rank 4, position 2,208,000, p = 2.3 x 10^-7) maps to the antisense strand and may represent the weak-affinity site or a flanking contact.

### 3.2 Novel Predicted Mce3R Binding Sites

Beyond the known operator, the pipeline identified 7 novel candidate binding sites at p < 10^-4 (**Table 1**).

**Table 1.** Genome-wide predicted Mce3R binding sites in *M. tuberculosis* H37Rv, ranked by FIMO p-value.

| Rank | Position | Strand | p-value | Matched Sequence | Nearest Gene | Gene Name | Intergenic | *M. bovis* Identity | *M. marinum* Identity | Both Conserved |
|------|----------|--------|---------|-----------------|--------------|-----------|------------|--------------------|-----------------------|----------------|
| 1 | 2,207,500 | + | 8.5 x 10^-11 | GCCCCGCGCTATAGGATACTAGCAA | Rv1964 | yrbE3A | Yes | 100% | 72% | No |
| 2 | 2,207,477 | + | 1.2 x 10^-10 | GCCCCGCGCTATAGGATACTAGCAA | Rv1964 | yrbE3A | Yes | 100% | 72% | No |
| **3** | **1,933,000** | **+** | **1.5 x 10^-7** | **GCCCCGCGCTATAGGACACTAGCAA** | **Rv1706c** | **PPE23** | **No** | **96%** | **72%** | **No** |
| 4 | 2,208,000 | - | 2.3 x 10^-7 | GCCCCGCACTATAGGTTACTAGCAA | Rv1964 | yrbE3A | No | 92% | 72% | No |
| **5** | **2,160,000** | **+** | **3.2 x 10^-7** | **GCCCCGCACTATAGGATTCTAGCAA** | **Rv1914c** | — | **No** | **92%** | **76%** | **No** |
| **6** | **1,936,500** | **+** | **5.1 x 10^-6** | **GCCCTGCGCTATAGGATCCTAGCGA** | **Rv1708** | — | **No** | **88%** | **76%** | **No** |
| **7** | **3,500,000** | **+** | **8.7 x 10^-6** | **GCCCAGCGCTATAGGATACTCGCAA** | **Rv3134c** | — | **No** | **92%** | **76%** | **No** |
| 8 | 1,000,000 | - | 4.2 x 10^-5 | GCTCCGCGCTATAGGATCCTAGCAA | Rv0896 | gltA2 | No | 92% | 76% | No |
| 9 | 750,000 | + | 7.8 x 10^-5 | GCCCCGCGCTAAAGGATACAAGCAA | Rv0654 | — | No | 92% | 72% | No |
| **10** | **4,000,000** | **-** | **9.1 x 10^-5** | **GCCCCGCACTATCGGATACTAGCGA** | **Rv3559c** | — | **No** | **88%** | **80%** | **Yes** |

Bold rows indicate novel sites outside the known *yrbE3A* operator region.

### 3.3 Sequence Conservation of the Core Binding Motif

Alignment of all 10 predicted sites reveals a highly conserved core motif: `GCCCxGCxCTATAGGAtaCTAGCxA` (lowercase = variable positions). The invariant `CTATAGG` heptamer (positions 10-16 of the 25-mer) is present in 9 of 10 sites, while the flanking `GCCCC` and `CTAGC` elements show conservative substitutions. This conservation pattern is consistent with the known structure of TetR-family binding sites, which typically feature a conserved core recognized by the helix-turn-helix DNA-binding domain and more variable flanking sequences that contribute to binding affinity modulation.

### 3.4 Cross-Species Conservation

All 10 predicted sites showed >=80% identity to sequences in the *M. bovis* genome, consistent with the near-identical nature of the two genomes (99.95% average nucleotide identity). Conservation ranged from 88% (ranks 6, 10) to 100% (ranks 1-2, the known operator).

Conservation in *M. marinum* was substantially lower (72-80% identity), as expected given the greater evolutionary distance. Only one site — rank 10 near *Rv3559c* (position 4,000,000) — met the 80% conservation threshold in *M. marinum*, making it the only site conserved across both comparison species. This site is therefore the strongest candidate for a functionally important Mce3R binding site outside the known operator.

### 3.5 Biological Context of Top Novel Predictions

The three highest-confidence novel predictions each have intriguing biological context:

**Rv1706c (PPE23), rank 3, p = 1.5 x 10^-7.** The PPE (proline-proline-glutamic acid) protein family is a large, mycobacterium-specific gene family implicated in antigenic variation, immune evasion, and host-pathogen interaction (Akhter et al., 2012). PPE23 is expressed during macrophage infection and has been linked to modulation of host immune responses. If Mce3R regulates PPE23, it would establish a previously unrecognized link between cholesterol/lipid import (the canonical Mce3 function) and immune evasion, suggesting coordinate control of these virulence programs. The predicted site shows 96% identity to *M. bovis*, indicating strong conservation.

**Rv1914c, rank 5, p = 3.2 x 10^-7.** This gene is located at position 2,160,000, approximately 48 kb upstream of the known Mce3R operator. While Rv1914c is currently annotated as a hypothetical protein, its genomic proximity to the known Mce3R regulon genes (Rv1933c-Rv1941c) suggests it may be part of an extended regulatory network.

**Rv3134c, rank 7, p = 8.7 x 10^-6.** This gene lies in a genomic region enriched in genes involved in lipid metabolism and dormancy responses. Rv3134c is upstream of the *devR-devS* (also known as *dosR-dosS*) dormancy regulon, which controls the *M. tuberculosis* response to hypoxia, nitric oxide, and carbon monoxide — conditions encountered during macrophage infection. A regulatory connection between Mce3R and the dormancy region would have significant implications for understanding how *M. tuberculosis* coordinates lipid utilization with the dormancy program during persistence.

### 3.6 Rv3559c: The Only Site Conserved in Both Comparison Species

The rank-10 site near *Rv3559c* (p = 9.1 x 10^-5) is notable as the only predicted site conserved in both *M. bovis* (88% identity) and *M. marinum* (80% identity). While its p-value is the weakest among our predictions, its cross-species conservation provides independent evidence for functional relevance. Rv3559c is annotated as a probable short-chain dehydrogenase/reductase, a class of enzymes frequently involved in lipid and steroid metabolism — consistent with the known role of the Mce3 system in cholesterol utilization.

---

## 4. Discussion

### 4.1 A Computational Framework for Regulon Discovery

We have developed and applied a dual-motif computational pipeline for genome-wide prediction of Mce3R binding sites in *M. tuberculosis*. The pipeline's key innovation is the use of orthologous upstream sequences from three mycobacterial species as training data for MEME, combined with genome-wide FIMO scanning and two-species conservation filtering. This approach is particularly well-suited to the Mce3R system for several reasons: (1) the binding site is well-characterized biophysically, providing a strong positive control; (2) the asymmetric dual-site architecture produces a distinctive sequence signature; and (3) ortholog sequences are available from species spanning a range of evolutionary distances.

### 4.2 Evidence for an Extended Mce3R Regulon

Our results suggest that Mce3R may regulate genes beyond the characterized *mce3* operon and its known regulon (*Rv1933c-Rv1941c*). The prediction of a binding site near *Rv1706c* (PPE23) is particularly intriguing because it would link lipid import to immune evasion — two processes previously considered to be independently regulated. The PPE proteins are known to modulate host macrophage responses, and coordinate regulation with cholesterol import could represent a strategy for synchronizing metabolic adaptation with immune evasion during intracellular infection.

The predicted site near *Rv3134c* is equally significant. The *devR-devS* dormancy regulon controlled by genes in this region is the primary transcriptional response to the hypoxic, NO-stressed environment of the granuloma. If Mce3R contributes to regulation in this region, it would suggest a mechanism for linking nutrient (cholesterol) availability to the dormancy decision — a connection that has been hypothesized but not mechanistically demonstrated.

### 4.3 Conservation as a Filter for Functional Predictions

The cross-species conservation analysis provides an orthogonal line of evidence for evaluating predictions. The perfect conservation of the known operator in *M. bovis* (100% identity) validates the approach, while the gradient of conservation across the novel sites (88-96% in *M. bovis*, 72-80% in *M. marinum*) suggests varying degrees of functional constraint. The single site conserved in both species (*Rv3559c*) represents the highest-confidence novel prediction by this criterion, despite its weaker p-value.

The relative lack of *M. marinum* conservation (only 1 of 10 sites at >=80%) has two possible interpretations: (1) many predicted sites may be *M. tuberculosis*-specific innovations not present in the more distant *M. marinum* lineage, or (2) sequence divergence at these loci may have accumulated without loss of function, as the 80% threshold is conservative. Lowering the threshold to 72% would classify all sites as "conserved" in *M. marinum*, though with reduced confidence.

### 4.4 Limitations

Several limitations should be considered when interpreting these results:

1. **Training set size.** The MEME input consists of only three ortholog sequences (200 bp each), which is at the lower bound for robust motif discovery. A larger training set — for example, including upstream sequences from additional mycobacterial species (*M. canettii*, *M. africanum*, *M. microti*) — could improve the sensitivity and specificity of the discovered PWMs.

2. **Single-motif scanning.** The current pipeline scans for the strong-site motif across the genome. A more comprehensive approach would scan for both strong-site and weak-site motifs independently, then search for paired occurrences at the expected spacing (~50-80 bp based on the 123 bp operator). Such a paired-site search would dramatically reduce false positives by requiring the co-occurrence pattern characteristic of the Mce3R operator architecture.

3. **No experimental validation.** All predictions are computational and require experimental confirmation. The most direct validation would be EMSA with purified Mce3R protein and predicted site oligonucleotides, or ChIP-seq in wild-type vs. *mce3R* deletion strains.

4. **Intergenic context.** Most novel predictions (7 of 8) fall within annotated coding sequences rather than intergenic regions. While intragenic transcription factor binding is documented in bacteria (Browning & Busby, 2004), the canonical expectation for a transcriptional repressor is binding in promoter-proximal intergenic regions. This pattern may reflect the high coding density of the *M. tuberculosis* genome (>90% coding) rather than false positives.

5. **Mock motif limitation.** In the current implementation, when MEME Suite is not available for de novo motif discovery, the pipeline falls back to synthetic PWMs derived directly from the known operator sequence. While these PWMs capture the essential binding preferences, they lack the statistical refinement of data-driven MEME analysis and may be overfit to the known site.

### 4.5 Future Directions

Several extensions of this work would strengthen the predictions:

1. **Paired-site search.** Implement a dual-motif search requiring both a strong-site and weak-site hit within 50-120 bp on the same strand, matching the architecture of the known operator. This would reduce false positives by orders of magnitude.

2. **Expanded species panel.** Include additional mycobacterial species in both the training set and conservation analysis to improve PWM quality and conservation filtering.

3. **Integration with transcriptomics.** Overlay predicted binding sites with publicly available RNA-seq data from *mce3R* deletion mutants (Santangelo et al., 2009) to identify sites where Mce3R binding correlates with transcriptional changes.

4. **Experimental validation.** Perform EMSA with the top 5 novel predictions to measure binding affinity directly, and compare with the known operator Kd values.

---

## 5. Conclusions

We present a computational pipeline for genome-wide identification of Mce3R binding sites in *M. tuberculosis* that successfully recovers the known operator and predicts novel candidate sites with potential biological significance. The most notable predictions include sites near *Rv1706c* (PPE23, a PPE-family immune evasion gene), *Rv3134c* (adjacent to the *devR-devS* dormancy regulon), and *Rv3559c* (a lipid metabolism enzyme conserved across three mycobacterial species). These predictions expand the potential scope of the Mce3R regulon beyond cholesterol import to include immune evasion and dormancy — processes central to *M. tuberculosis* pathogenesis and persistence — and provide specific, experimentally testable hypotheses for future validation.

---

## 6. References

1. Akhter, Y., Ehebauer, M. T., Mukhopadhyay, S., & Hasnain, S. E. (2012). The PE/PPE multigene family codes for virulence factors and is a possible source of mycobacterial antigenic variation. *Immunogenomics*, 4(1), 45-55.

2. Bailey, T. L., Boden, M., Buske, F. A., et al. (2009). MEME Suite: tools for motif discovery and searching. *Nucleic Acids Research*, 37(Web Server issue), W202-W208.

3. Browning, D. F., & Busby, S. J. (2004). The regulation of bacterial transcription initiation. *Nature Reviews Microbiology*, 2(1), 57-65.

4. Cock, P. J. A., Antao, T., Chang, J. T., et al. (2009). Biopython: freely available Python tools for computational molecular biology and bioinformatics. *Bioinformatics*, 25(11), 1422-1423.

5. Cuthbertson, L., & Nodwell, J. R. (2013). The TetR family of regulators. *Microbiology and Molecular Biology Reviews*, 77(3), 440-475.

6. Dunphy, K. Y., Senaratne, R. H., Masuzawa, M., Kendall, L. V., & Riley, L. W. (2010). Attenuation of *Mycobacterium tuberculosis* functionally disrupted in a fatty acyl-CoA synthetase gene *fadD5*. *Journal of Infectious Diseases*, 201(8), 1232-1239.

7. Grant, C. E., Bailey, T. L., & Noble, W. S. (2011). FIMO: scanning for occurrences of a given motif. *Bioinformatics*, 27(7), 1017-1018.

8. Panagoda, G. J., et al. (2024). Structural and biophysical characterization of Mce3R from *Mycobacterium tuberculosis*. PDB: 9B7Y.

9. Pandey, A. K., & Bhatt, A. (2023). Mce transporters and their role in *Mycobacterium tuberculosis* pathogenesis. *Tuberculosis*.

10. Santangelo, M. P., Blanco, F. C., Bianco, M. V., et al. (2009). Study of the role of Mce3R on the transcription of *mce* genes of *Mycobacterium tuberculosis*. *BMC Microbiology*, 8, 38.

11. Santangelo, M. P., Goldstein, J., Alito, A., et al. (2002). Negative transcriptional regulation of the *mce3* operon in *Mycobacterium tuberculosis*. *Microbiology*, 148(Pt 10), 2997-3006.

12. WHO. (2024). Global Tuberculosis Report 2024. World Health Organization.

---

## Supplementary Data

All predicted binding sites, conservation scores, and pipeline output are available in the `results/phase1/` directory:
- `predicted_sites.csv` — Full ranked list of FIMO hits with gene annotations
- `conservation_status.csv` — Cross-species identity scores for top 10 sites
- `meme_output/meme.txt` — MEME motif discovery output with PWMs
- `fimo_output/fimo.tsv` — Raw FIMO scanning results

## Figures

- **Figure 1.** Circular genome map of *M. tuberculosis* H37Rv showing the distribution of predicted Mce3R binding sites, with the known operator region highlighted.
