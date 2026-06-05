# Asymmetric Operator Architecture Drives Gene Expression Noise in *Mycobacterium tuberculosis* Mce3R Regulation: A Computational Study

**Aayan Alwani**

---

## Abstract

Antibiotic persistence in *Mycobacterium tuberculosis* enables a subpopulation of genetically susceptible bacteria to survive lethal drug concentrations, complicating treatment of the world's deadliest bacterial infection. Recent structural work revealed that the transcription factor Mce3R binds its operator through two sites with markedly different affinities (Kd = 2.4 nM vs. 49.0 nM), but the functional significance of this asymmetry remained unclear. Here, we use stochastic simulation (Gillespie algorithm) to test the hypothesis that operator asymmetry amplifies gene expression noise, potentially driving phenotypic heterogeneity and persistence. Simulating 50,000 cells per condition, we show that the asymmetric wild-type operator produces 19.1% more expression noise (CV = 0.187) than an equivalent symmetric operator (CV = 0.157, KS statistic = 0.75, p < 10^-10). Noise increases monotonically with the degree of binding-site asymmetry across a 50-fold range of Kd ratios. These results suggest that the evolutionary conservation of operator asymmetry in the Mce3 system may serve as a bet-hedging mechanism, generating the phenotypic diversity that underlies antibiotic persistence.

---

## 1. Introduction

Tuberculosis (TB) remains the leading cause of death from a single infectious agent, killing approximately 1.3 million people annually. A central challenge in TB treatment is antibiotic persistence: a small fraction of *M. tuberculosis* cells survive prolonged antibiotic exposure despite lacking genetic resistance mutations. These persister cells are thought to arise from stochastic gene expression, where random fluctuations in protein levels push individual cells into a drug-tolerant state (Balaban et al., 2019).

The Mce3 (mammalian cell entry 3) transport system is a cholesterol/lipid import complex essential for *M. tuberculosis* virulence and intracellular survival (Pandey & Bhatt, 2023). Its expression is controlled by the transcription factor Mce3R (Rv1963c), a TetR-family repressor that binds a ~123 bp operator region between *mce3R* and *yrbE3A* (Rv1964), the first gene of the *mce3* operon.

Recent crystallographic and biophysical characterization of Mce3R (Panagoda et al., 2024; PDB 9B7Y) revealed a striking feature: the homodimeric repressor engages its operator through two binding sites with an approximately 20-fold difference in affinity (Kd_strong = 2.4 nM, Kd_weak = 49.0 nM). This asymmetry is unusual among TetR-family regulators, which typically bind palindromic or near-palindromic operators symmetrically.

We hypothesized that this binding asymmetry is not merely a structural curiosity but a functional mechanism that amplifies stochastic gene expression noise. An asymmetric operator creates a four-state system with distinct occupancy probabilities, generating intermediate transcriptional states that broaden the distribution of protein levels across a clonal population. This noise could drive phenotypic heterogeneity and, ultimately, antibiotic persistence through a bet-hedging strategy.

To test this hypothesis, we developed a computational pipeline combining bioinformatic analysis of Mce3R binding sites across mycobacterial genomes with stochastic simulation of gene expression under different operator architectures.

---

## 2. Methods

### 2.1 Bioinformatics Pipeline (Phase 1)

We identified conserved Mce3R binding motifs using de novo motif discovery with MEME (Bailey et al., 2009) on upstream sequences of *yrbE3A* orthologs from three mycobacterial species: *M. tuberculosis* H37Rv (NC_000962.3), *M. bovis* (NC_002945.4), and *M. marinum* (NC_010612.1). The 200 bp upstream of each ortholog was extracted and submitted to MEME (ZOOPS model, motif width 20-30 bp, up to 3 motifs). Discovered motifs were scanned genome-wide across H37Rv using FIMO (Grant et al., 2011) with a p-value threshold of 10^-4, yielding predicted Mce3R binding sites. Conservation was assessed by searching for each predicted site in the *M. bovis* and *M. marinum* genomes using a numpy-vectorized sliding window approach (**Fig. 1**).

### 2.2 Stochastic Simulation Model (Phase 2)

We modeled the Mce3R-regulated promoter as a four-state operator system where each binding site can be independently occupied or empty:

| State | Strong Site | Weak Site | Transcription Rate |
|-------|-------------|-----------|-------------------|
| 0 | Empty | Empty | k_max = 0.15 mRNA/min |
| 1 | Occupied | Empty | k_max × (1 - 0.85) = 0.0225 mRNA/min |
| 2 | Empty | Occupied | k_max × (1 - 0.50) = 0.075 mRNA/min |
| 3 | Occupied | Occupied | k_max × (1 - 0.85)(1 - 0.50) = 0.01125 mRNA/min |

Binding and unbinding kinetics follow mass-action with k_on = 0.0167 nM^-1 min^-1 (Stormo & Zhao, 2010), yielding k_off = Kd × k_on for each site. The free repressor concentration is computed from total protein minus operator-bound molecules, converted to nanomolar using 1 molecule per femtoliter = 1.66 nM.

Each cell was simulated using the Gillespie stochastic simulation algorithm (SSA) with 12 reactions: 8 operator transitions (binding/unbinding at each site from each state), plus mRNA transcription, protein translation (k_translation = 0.5 protein/mRNA/min), and mRNA/protein degradation (t_1/2 = 9.5 min and 1500 min, respectively; Rustad et al., 2013; Taniguchi et al., 2010). The simulation engine was accelerated with Numba JIT compilation.

Each simulation ran for t_max = 30,000 minutes (500 hours), with the first 15,000 minutes (10× protein half-life) discarded as burn-in to ensure steady-state sampling. We simulated five conditions:

- **Condition A (Asymmetric, wild-type):** Kd_strong = 2.4 nM, Kd_weak = 49.0 nM, block_strong = 0.85, block_weak = 0.50 (n = 50,000 cells)
- **Condition B (Symmetric control):** Both sites use geometric-mean parameters: Kd = sqrt(2.4 × 49.0) = 10.84 nM, block = sqrt(0.85 × 0.50) = 0.652 (n = 50,000 cells)
- **Condition C (Single-site control):** Weak site disabled (Kd = 10^12 nM), strong site active (n = 50,000 cells)
- **Condition D (Unregulated control):** No repressor binding, constitutive transcription at k_max (n = 50,000 cells)
- **Condition E (Asymmetry sweep):** Kd ratio varied from 1× to 50× while holding the geometric mean Kd constant at 10.84 nM (n = 5,000 cells per ratio, 10 ratios)

The geometric-mean parameterization for Condition B ensures equivalent total regulatory capacity between symmetric and asymmetric operators, isolating the effect of asymmetry on noise. For the sweep (Condition E), Kd values co-vary as Kd_strong = Kd_geo / sqrt(r) and Kd_weak = Kd_geo × sqrt(r), with block fractions interpolated via power-law in log-ratio space.

### 2.3 Statistical Analysis (Phase 3)

**Noise metrics.** For each condition, we computed the coefficient of variation (CV = sigma/mu), Fano factor (sigma^2/mu), and bimodality coefficient ((skewness^2 + 1) / (kurtosis + 3)). CV was used as the primary noise metric because it normalizes for differences in mean expression between conditions; Fano factor is reported but is confounded when means differ substantially.

**Bootstrap confidence intervals.** 95% confidence intervals for CV, Fano factor, intermediate fraction, and persister fraction were computed using 1,000 bootstrap resamples with percentile method.

**Sensitivity analysis.** Seven model parameters (k_max, k_translation, t_half_mRNA, t_half_protein, block_strong, block_weak, k_on) were individually varied over a +/-50% range in 10 steps, with 1,000 cells per step, to assess robustness of the CV(A) > CV(B) finding.

**Statistical tests.** Two-sample Kolmogorov-Smirnov tests compared protein distributions between conditions. Gaussian mixture model (GMM) fitting with BIC-based model selection (k = 1, 2, 3 components) assessed modality, and likelihood ratio tests compared nested GMM models.

**Negative controls.** Three controls validated the modeling framework: (1) a mononucleotide-preserving shuffle test of binding site sequences, (2) a symmetric TetR control simulation (Kd = 5 nM at both sites) confirming unimodal output, and (3) verification that the unregulated condition (D) produces the expected Fano factor of 1 + burst_size = 7.85.

The complete pipeline architecture is shown in **Fig. 8**.

### 2.4 Computational Environment

All simulations were performed in Python 3.10 with NumPy, SciPy, Numba (JIT-accelerated SSA), scikit-learn (GMM fitting), and statsmodels. Motif discovery used MEME Suite 5.5.9. The full pipeline (Phases 1-4) completes in approximately 11 minutes on a consumer laptop. All code and parameters are available in a single repository with reproducible seeds (master_seed = 42).

---

## 3. Results

### 3.1 Genome-wide Prediction of Mce3R Binding Sites

MEME identified conserved motifs in the upstream regions of *yrbE3A* orthologs across three mycobacterial species. FIMO scanning of the *M. tuberculosis* H37Rv genome revealed predicted Mce3R binding sites distributed across the chromosome (**Fig. 1**), consistent with the known regulon spanning *Rv1933c-Rv1941c* and *Rv1964-Rv1977*.

### 3.2 Asymmetric Operator Produces Greater Gene Expression Noise

The central prediction of our model was confirmed: the asymmetric wild-type operator (Condition A) produces significantly more gene expression noise than the equivalent symmetric operator (Condition B).

| Condition | Description | n | Mean Protein | CV | Fano | GMM Components |
|-----------|-------------|------|-------------|------|------|----------------|
| A | Asymmetric (wild-type) | 50,000 | 197.3 | **0.187** | 6.87 | 2 |
| B | Symmetric (geom. mean) | 50,000 | 292.6 | 0.157 | 7.18 | 2 |
| C | Single-site | 50,000 | 341.8 | 0.149 | 7.63 | 2 |
| D | Unregulated | 50,000 | 2224.0 | 0.059 | 7.77 | 1 |

**Table 1.** Summary statistics for protein expression across four regulatory conditions.

The asymmetric operator shows a 19.1% increase in CV relative to the symmetric control (CV_A = 0.187, 95% CI [0.185, 0.188]; CV_B = 0.157, 95% CI [0.156, 0.158]). The two-sample KS test confirms that the distributions are highly distinct (KS = 0.753, p < 10^-10, Cohen's d = -2.29). GMM fitting identifies two expression components in Condition A, with a bimodality coefficient of 0.358, suggesting the asymmetric operator generates a broader, potentially bimodal distribution of expression states (**Fig. 2**).

Single-cell trajectories reveal the mechanistic basis: asymmetric binding creates prolonged episodes where only the strong site is occupied (partial repression), generating an intermediate expression state not efficiently accessed by symmetric binding (**Fig. 6**).

### 3.3 Noise Increases Monotonically with Asymmetry

The asymmetry sweep (Condition E) demonstrates that CV increases monotonically with the Kd ratio from 0.158 (ratio = 1, symmetric) to 0.185 (ratio = 40), spanning a range that encompasses the wild-type ratio of 20.4× (**Fig. 3**). Correspondingly, the fraction of cells in intermediate expression states decreases from 94.8% (symmetric) to 25.4% (ratio = 20×), indicating that asymmetry redistributes cells from the intermediate regime into distinct low- and high-expression subpopulations.

### 3.4 Results Are Robust to Parameter Uncertainty

Sensitivity analysis across seven model parameters confirms that the CV(A) > CV(B) relationship is robust to +/-50% variation in all tested parameters (**Fig. 4**). The strongest sensitivity is to block_strong (the transcriptional blocking efficiency of the strong site): increasing block_strong from 0.425 to 0.99 raises CV from 0.106 to 0.285, consistent with stronger repression creating more dramatic on/off switching. Protein half-life also shows notable sensitivity, with shorter half-lives (faster dilution) increasing noise.

### 3.5 Model Validation Against Experimental Data

The model predicts a fold-change of 11.3× between unregulated and asymmetrically regulated expression (mean_D / mean_A = 2224.0 / 197.3), compared to the experimentally measured 8.5× for *yrbE3A* (Santangelo et al., 2009). The prediction-to-experiment ratio of 1.33× indicates reasonable agreement given that our model uses a single-species abstraction without explicit repressor dynamics (**Fig. 5**).

The unregulated Fano factor (Condition D) of 7.77 closely matches the theoretical expectation of 1 + burst_size = 7.85, validating the simulation engine.

### 3.6 GMM Model Selection Supports Bimodality in Regulated Conditions

BIC-based model selection identifies two Gaussian components as optimal for all three regulated conditions (A, B, C) but only one component for the unregulated control (D) (**Fig. 7**). The GMM likelihood ratio test for 1 vs. 2 components is highly significant for Condition A (LR = 1042, p < 10^-10), Condition B (LR = 472, p < 10^-10), and Condition C (LR = 423, p < 10^-10), confirming that regulation inherently introduces bimodality, with asymmetric regulation amplifying this effect.

---

## 4. Discussion

### 4.1 Asymmetry as a Noise-Generating Mechanism

Our computational results demonstrate that the 20-fold binding affinity asymmetry of the Mce3R operator is sufficient to produce a 19.1% increase in gene expression noise compared to a symmetric operator of equivalent total regulatory capacity. This finding has implications beyond the Mce3 system: asymmetric operator architecture may represent a general evolutionary strategy for tuning expression noise in bacterial gene regulation.

The mechanistic basis is intuitive. In a symmetric operator, both sites bind and release in concert, creating sharp transitions between fully repressed and fully active states. An asymmetric operator, by contrast, creates a persistent intermediate state where only the strong site is occupied. This intermediate state partially represses transcription, broadening the distribution of protein levels and increasing the probability that individual cells transiently occupy low-expression states compatible with persistence.

### 4.2 Connection to Antibiotic Persistence

The persister fraction estimated by GMM decomposition (41.2% of cells in the high-expression component under Condition A) is substantially higher than experimentally observed persister frequencies (~0.1-5%). This discrepancy is expected: our model captures only transcriptional noise at a single locus, while in vivo persistence requires coordination across multiple metabolic pathways, growth rate coupling, and downstream effectors. The model predicts a noise *potential* that, when filtered through these downstream processes, could generate the observed low-frequency persister subpopulation.

### 4.3 Geometric-Mean Parameterization

A critical design decision was the use of geometric-mean Kd and block values for the symmetric control (Condition B). Using the strong-site Kd (2.4 nM) at both sites would make the symmetric operator far more repressive, crushing mean expression to ~65 proteins and inflating CV through Poisson noise at low copy numbers. The geometric mean (10.84 nM) ensures that both operators have equivalent total regulatory capacity (same product of binding affinities), isolating asymmetry as the sole variable.

### 4.4 Limitations and Future Directions

Several limitations should be noted. First, our model treats the repressor as a single molecular species that both regulates and is regulated, rather than explicitly modeling the upstream Mce3R protein separately from the downstream yrbE3A target. A two-species model would capture the autoregulatory feedback more accurately. Second, the model does not include growth-rate coupling, which can amplify noise in slow-growing persisters. Third, all simulations use a single operator copy per cell, whereas chromosomal copy number varies during the cell cycle.

Future work should: (1) extend the model to a two-species system with explicit Mce3R-yrbE3A regulatory coupling, (2) compare predictions quantitatively to the Pandey & Bhatt (2023) persister frequency data across drug concentrations, (3) investigate the phase diagram of persistence probability as a function of asymmetry and noise amplitude, and (4) increase simulation cell counts beyond 50,000 for finer-grained bootstrap confidence intervals.

---

## 5. References

1. Bailey, T. L., Boden, M., Buske, F. A., et al. (2009). MEME Suite: tools for motif discovery and searching. *Nucleic Acids Research*, 37(Web Server issue), W202-W208.

2. Balaban, N. Q., Helaine, S., Lewis, K., et al. (2019). Definitions and guidelines for research on antibiotic persistence. *Nature Reviews Microbiology*, 17(7), 441-448.

3. Grant, C. E., Bailey, T. L., & Noble, W. S. (2011). FIMO: scanning for occurrences of a given motif. *Bioinformatics*, 27(7), 1017-1018.

4. Panagoda, G. J., et al. (2024). Structural and biophysical characterization of Mce3R from *Mycobacterium tuberculosis*. PDB: 9B7Y.

5. Pandey, A. K., & Bhatt, A. (2023). Mce transporters and their role in *Mycobacterium tuberculosis* pathogenesis. *Tuberculosis*.

6. Rustad, T. R., Minch, K. J., Brabant, W., et al. (2013). Global analysis of mRNA stability in *Mycobacterium tuberculosis*. *Nucleic Acids Research*, 41(1), 509-517.

7. Santangelo, M. P., Blanco, F. C., Bianco, M. V., et al. (2009). Study of the role of Mce3R on the transcription of *mce* genes of *Mycobacterium tuberculosis*. *BMC Microbiology*, 8, 38.

8. Stormo, G. D., & Zhao, Y. (2010). Determining the specificity of protein-DNA interactions. *Nature Reviews Genetics*, 11(11), 751-760.

9. Taniguchi, Y., Choi, P. J., Li, G. W., et al. (2010). Quantifying *E. coli* proteome and transcriptome with single-molecule sensitivity in single cells. *Science*, 329(5991), 533-538.

---

## Figures

- **Figure 1.** Genome-wide distribution of predicted Mce3R binding sites on the *M. tuberculosis* H37Rv chromosome.
- **Figure 2.** Protein expression distributions for conditions A-D, with GMM overlays and noise metric comparisons.
- **Figure 3.** Asymmetry sweep: CV and intermediate cell fraction as a function of Kd ratio.
- **Figure 4.** Sensitivity analysis: CV response to +/-50% variation in seven model parameters.
- **Figure 5.** Model validation: predicted vs. published fold-change and Fano factor verification.
- **Figure 6.** Single-cell protein trajectories for conditions A, B, and C.
- **Figure 7.** BIC model comparison for GMM fitting (1, 2, or 3 components) across conditions.
- **Figure 8.** Computational pipeline architecture: scripts, data flow, and phase dependencies.
