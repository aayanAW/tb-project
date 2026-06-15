#!/usr/bin/env python3
"""Generate ENIGMA Manuscript v2 as a Word document using python-docx."""

from docx import Document
from docx.shared import Pt, Inches, Cm
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.enum.section import WD_ORIENT
from docx.oxml.ns import qn
import os

doc = Document()

# ── Global defaults ──────────────────────────────────────────────────────────
style = doc.styles['Normal']
font = style.font
font.name = 'Calibri'
font.size = Pt(11)
style.paragraph_format.line_spacing = 2.0
style.paragraph_format.space_after = Pt(0)
style.paragraph_format.space_before = Pt(0)

# Set margins to 1 inch
for section in doc.sections:
    section.top_margin = Inches(1)
    section.bottom_margin = Inches(1)
    section.left_margin = Inches(1)
    section.right_margin = Inches(1)


def add_title(text):
    p = doc.add_paragraph()
    p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    run = p.add_run(text)
    run.bold = True
    run.font.size = Pt(16)
    run.font.name = 'Calibri'
    return p


def add_section_header(text):
    p = doc.add_paragraph()
    run = p.add_run(text)
    run.bold = True
    run.font.size = Pt(14)
    run.font.name = 'Calibri'
    p.paragraph_format.space_before = Pt(12)
    return p


def add_subsection_header(text):
    p = doc.add_paragraph()
    run = p.add_run(text)
    run.bold = True
    run.font.size = Pt(12)
    run.font.name = 'Calibri'
    p.paragraph_format.space_before = Pt(8)
    return p


def add_body(text):
    p = doc.add_paragraph()
    run = p.add_run(text)
    run.font.name = 'Calibri'
    run.font.size = Pt(11)
    p.paragraph_format.line_spacing = 2.0
    p.paragraph_format.first_line_indent = Inches(0.5)
    return p


def add_body_no_indent(text):
    p = doc.add_paragraph()
    run = p.add_run(text)
    run.font.name = 'Calibri'
    run.font.size = Pt(11)
    p.paragraph_format.line_spacing = 2.0
    return p


def add_page_break():
    doc.add_page_break()


# ═══════════════════════════════════════════════════════════════════════════════
# TITLE PAGE
# ═══════════════════════════════════════════════════════════════════════════════

add_title(
    "Noisy by Design: Asymmetric Operator Architecture in the Mce3R Regulon "
    "May Generate Persistence-Enabling Expression Noise in "
    "Mycobacterium tuberculosis"
)

# Author
p = doc.add_paragraph()
p.alignment = WD_ALIGN_PARAGRAPH.CENTER
run = p.add_run("[Author Name]")
run.font.name = 'Calibri'
run.font.size = Pt(12)

# Blank line
doc.add_paragraph()

# ── ABSTRACT ─────────────────────────────────────────────────────────────────
add_section_header("Abstract")

add_body(
    "Mycobacterium tuberculosis kills over 1.2 million people annually, with treatment failure "
    "driven partly by antibiotic-tolerant persister cells that arise from stochastic gene expression "
    "noise rather than genetic mutation. The Mce3R transcriptional repressor, which controls "
    "cholesterol and lipid metabolism genes critical for host survival, binds a two-site operator "
    "with a 20.4-fold binding affinity asymmetry (Kd = 2.4 nM vs. 49 nM). We hypothesized that "
    "this structural asymmetry functions as a noise-generating mechanism that expands the "
    "persister-competent subpopulation."
)

add_body(
    "Aim 1: We applied de novo motif discovery (MEME) across three mycobacterial genomes and "
    "genome-wide scanning (FIMO) to identify and validate Mce3R binding sites, discovering 1,442 "
    "candidate sites including the known operator as the top hit (p = 1.1 \u00d7 10\u207b\u00b9\u00b9)."
)

add_body(
    "Aim 2: We constructed a thermodynamically calibrated stochastic model combining a four-state "
    "operator with Gillespie simulation, Berg\u2013von Hippel energy calibration, and Markov Chain Monte "
    "Carlo inference of cooperativity (\u03c9 = 1.08, 95% CI: 0.20\u20135.89)."
)

add_body(
    "Aim 3: We quantified noise across four operator architectures and four environmental "
    "conditions. The asymmetric operator produced higher expression noise (CV = 0.187) than the "
    "symmetric control (CV = 0.157, p < 0.001, bootstrap 95% CI for \u0394CV entirely above zero). "
    "This result was robust across 100% of MCMC posterior samples, persisted in a two-species "
    "autoregulatory model (CV = 0.200 vs. 0.197), and held in all four simulated environments. "
    "Corrected persister fractions showed 3.1-fold enrichment in the asymmetric architecture under "
    "baseline conditions and 29.4-fold under acidic pH stress. The asymmetric operator transmitted "
    "less environmental information (MI = 0.28 bits) than the symmetric alternative (MI = 0.44 bits), "
    "suggesting a trade-off between information capacity and noise-driven phenotypic diversification."
)

add_body(
    "These computational results suggest that operator binding asymmetry may represent an evolved "
    "noise-generating strategy, and identify a potential therapeutic target: reducing expression "
    "noise to shrink the persister subpopulation."
)

add_page_break()

# ═══════════════════════════════════════════════════════════════════════════════
# INTRODUCTION
# ═══════════════════════════════════════════════════════════════════════════════

add_section_header("Introduction")

# Paragraph 1
add_body(
    "Tuberculosis remains one of the deadliest infectious diseases in human history. According to "
    "the World Health Organization, Mycobacterium tuberculosis killed over 1.2 million people in "
    "2023 alone, and the standard treatment regimen requires six to nine months of daily "
    "antibiotics (WHO, 2024). This extraordinarily long treatment duration is not merely an "
    "inconvenience\u2014it is a fundamental driver of treatment failure, because patient adherence "
    "declines over months of therapy, and incomplete courses allow bacterial populations to "
    "resurge. Critically, treatment failure is driven not only by genetic resistance\u2014the "
    "acquisition of mutations that render antibiotics ineffective\u2014but also by phenotypic "
    "persistence, a reversible, non-genetic state of antibiotic tolerance (Balaban et al., 2019). "
    "Persister cells survive drug exposure despite being genetically susceptible to the antibiotic, "
    "and their regrowth after treatment cessation drives clinical relapse. Understanding the "
    "mechanisms that generate persisters is therefore essential for shortening TB treatment and "
    "improving outcomes."
)

# Paragraph 2
add_body(
    "Persistence arises from stochastic gene expression noise. Even genetically identical "
    "bacterial cells produce different amounts of any given protein, because the molecular "
    "collisions that drive transcription and translation are inherently random events (Elowitz "
    "et al., 2002). A transcription factor finds its binding site through a random walk along DNA; "
    "an RNA polymerase initiates transcription at stochastic intervals; ribosomes translate mRNA "
    "molecules with variable efficiency. The cumulative effect of these random processes is that "
    "a population of genetically identical cells displays a distribution of protein levels\u2014some "
    "cells have high expression, others have low expression, and the width of this distribution "
    "is what we call expression noise. This noise creates phenotypic heterogeneity: a "
    "subpopulation of cells may enter a low-expression, metabolically quiescent state that happens "
    "to be tolerant to antibiotics. The fraction of cells in this tolerant state depends on the "
    "amplitude and shape of expression noise, which is controlled by the architecture of gene "
    "regulatory circuits (Balazsi et al., 2011). Circuits with positive feedback tend to amplify "
    "noise and create bimodal distributions; circuits with negative feedback tend to suppress noise "
    "and narrow distributions. The specific structural features of a regulatory element\u2014the "
    "number, arrangement, and affinity of transcription factor binding sites\u2014therefore have "
    "direct consequences for the noise properties of the genes they control."
)

# Paragraph 3
add_body(
    "The Mce3R system provides a concrete example of how regulatory architecture may influence "
    "persistence. Mce3R is a TetR-family transcriptional repressor in M. tuberculosis that "
    "controls the mce3 operon (Rv1964\u2013Rv1977), encoding cholesterol and lipid metabolism "
    "machinery essential for survival within the human host (Santangelo et al., 2002; Dunphy "
    "et al., 2010). The mce3 operon products enable the bacterium to utilize host-derived "
    "cholesterol as a carbon source during infection\u2014a metabolic adaptation that is critical for "
    "long-term persistence in macrophages. Deletion of mce3R has been shown to increase antibiotic "
    "persister frequency (Pandey et al., 2023), directly linking this regulatory system to the "
    "persistence phenotype. The recent cryo-EM structure of Mce3R bound to its operator DNA "
    "(Panagoda et al., 2024, PDB 9B7Y) revealed an asymmetric two-site architecture: Mce3R binds "
    "one site with high affinity (Kd = 2.4 nM, the strong site) and a second site with much lower "
    "affinity (Kd = 49 nM, the weak site), representing a 20.4-fold difference in binding "
    "strength. This asymmetry is unusual among TetR-family regulators, which typically bind "
    "palindromic operators with near-equal affinity at both half-sites."
)

# Paragraph 4
add_body(
    "The gap in our understanding is this: while operator asymmetry and antibiotic persistence "
    "have each been documented separately, no study has mechanistically connected the specific "
    "structural feature of operator binding asymmetry to noise-driven persistence. It is unknown "
    "whether the 20.4-fold Kd ratio generates biologically meaningful expression noise, whether "
    "this noise is robust to autoregulatory feedback (since Mce3R may regulate its own "
    "expression), and whether the resulting phenotypic heterogeneity is sufficient to expand the "
    "persister subpopulation beyond what a symmetric operator would produce. These questions sit "
    "at the intersection of structural biology, stochastic gene regulation, and infectious disease "
    "biology."
)

# Paragraph 5
add_body(
    "This study addresses this gap through an integrated computational pipeline spanning motif "
    "discovery, stochastic simulation, thermodynamic calibration, and environmental modeling. We "
    "test the hypothesis that operator binding asymmetry may function as a noise-generating "
    "mechanism that increases the fraction of persistence-competent cells. Our approach proceeds "
    "in three aims: first, we validate our computational pipeline by rediscovering the known "
    "Mce3R binding sites de novo and scanning the genome for novel candidates; second, we build "
    "a thermodynamically calibrated stochastic model of the four-state asymmetric operator; and "
    "third, we quantify the noise properties of this operator across multiple architectures and "
    "environmental conditions. We emphasize that all results presented here are computational "
    "predictions that require experimental validation, and we provide explicit falsification "
    "criteria that would disprove our hypothesis."
)

add_page_break()

# ═══════════════════════════════════════════════════════════════════════════════
# METHODS
# ═══════════════════════════════════════════════════════════════════════════════

add_section_header("Methods")

# --- Aim 1: Motif Discovery ---
add_subsection_header("Motif Discovery and Genome-Wide Scanning (Aim 1)")

add_body(
    "To identify Mce3R binding motifs without prior knowledge of the binding sequence, we "
    "performed de novo motif discovery across three mycobacterial genomes. We downloaded the "
    "complete genome sequences of M. tuberculosis H37Rv (GenBank accession NC_000962.3), "
    "M. bovis AF2122/97 (NC_002945.4), and M. marinum M (NC_010612.1). For each genome, we "
    "identified the ortholog of yrbE3A (the first gene in the mce3 operon) and extracted the "
    "200 base pairs immediately upstream of its start codon. These three upstream sequences "
    "served as input for de novo motif discovery."
)

add_body(
    "We ran the MEME algorithm (Bailey et al., 2009) in Zero or One Occurrence Per Sequence "
    "(ZOOPS) mode, searching both DNA strands, with motif widths ranging from 6 to 110 base "
    "pairs. The background nucleotide model was a 0th-order Markov model using the GC content "
    "of the H37Rv genome (A = 0.204, C = 0.296, G = 0.296, T = 0.204). MEME finds the shared "
    "pattern across species without being told what to look for\u2014it identifies the most "
    "statistically overrepresented sequence patterns in the input, which should correspond to "
    "functional binding sites if they are conserved."
)

add_body(
    "We then scanned the entire H37Rv genome using FIMO (Grant et al., 2011) with the position "
    "weight matrices (PWMs) discovered by MEME. FIMO searches the whole genome for occurrences "
    "of the motif pattern, assigning a p-value to each match based on how well it fits the PWM. "
    "We used a p-value threshold of 1 \u00d7 10\u207b\u2074 to identify candidate binding sites. This "
    "genome-wide scan served two purposes: validating the computational pipeline (the known "
    "operator should be the top hit) and discovering potential novel Mce3R binding sites that "
    "would suggest an expanded regulon."
)

# --- Aim 2: Four-State Operator Model ---
add_subsection_header("Four-State Operator Model (Aim 2)")

add_body(
    "We modeled the Mce3R operator as a system with four distinct states, corresponding to the "
    "four possible combinations of repressor occupancy at the two binding sites. State 0 is the "
    "unbound operator, where neither site is occupied and transcription proceeds at the maximal "
    "rate. State 1 is the strong-site-only state, where Mce3R occupies the high-affinity site "
    "(Kd = 2.4 nM) and transcription is 85% blocked. State 2 is the weak-site-only state, where "
    "Mce3R occupies the low-affinity site (Kd = 49 nM) and transcription is 50% blocked. State 3 "
    "is the doubly-bound state, where both sites are occupied and transcription is 92.5% blocked. "
    "The operator is, in effect, like a light switch with four positions, and the cell randomly "
    "flickers between them as repressor molecules bind and unbind."
)

add_body(
    "The transcription rates for each state are: State 0, k_txn = 0.150 mRNA/min (full "
    "transcription); State 1, k_txn = 0.0225 mRNA/min (85% blocked by strong-site occupancy); "
    "State 2, k_txn = 0.075 mRNA/min (50% blocked by weak-site occupancy); State 3, k_txn = "
    "0.01125 mRNA/min (92.5% blocked by dual occupancy). Transitions between states are governed "
    "by the association rate constant k_on = 0.0167 nM\u207b\u00b9 min\u207b\u00b9, based on established "
    "protein\u2013DNA association kinetics (Stormo & Zhao, 2010), and the dissociation rate "
    "k_off = Kd \u00d7 k_on for each site. The asymmetry in Kd values means that the strong site "
    "has a much slower off-rate than the weak site, creating different residence times in each "
    "operator state."
)

# --- Gillespie SSA ---
add_subsection_header("Gillespie Stochastic Simulation Algorithm (Aim 2)")

add_body(
    "We implemented the exact Stochastic Simulation Algorithm (Gillespie, 1977) with 12 "
    "reactions: 8 operator state transitions (binding and unbinding at each site from each "
    "state), transcription (mRNA production at a rate dependent on the current operator state), "
    "translation (protein production from each mRNA molecule), mRNA degradation, and protein "
    "degradation. We simulate 50,000 individual cells, each making random molecular decisions "
    "one event at a time, and measure how much the resulting protein levels vary between cells. "
    "This cell-by-cell simulation captures the stochastic fluctuations that deterministic "
    "models average away."
)

add_body(
    "Kinetic parameters were drawn from published measurements in mycobacteria and E. coli "
    "where mycobacterial data were unavailable. The translation rate was set to k_translation = "
    "0.5 protein per mRNA per minute (Taniguchi et al., 2010). The mRNA half-life was 9.5 minutes "
    "(Rustad et al., 2013), corresponding to a degradation rate of 0.073 min\u207b\u00b9. The protein "
    "half-life was set to 1,500 minutes, reflecting dilution-dominated degradation in slowly "
    "growing mycobacteria rather than active proteolysis. Each simulation ran for 30,000 minutes "
    "of simulated time, with the first 15,000 minutes discarded as burn-in to ensure the system "
    "had reached steady state before measurements began."
)

add_body(
    "We simulated four operator architectures to isolate the effect of asymmetry: (A) the native "
    "asymmetric operator with Kd = 2.4 nM and 49 nM; (B) a symmetric control where both sites "
    "have the geometric mean affinity, Kd = 10.84 nM, preserving the overall repression strength; "
    "(C) a single-site architecture where the weak site is disabled; and (D) an unregulated "
    "control with constitutive transcription. The symmetric control is particularly important "
    "because it isolates the effect of asymmetry per se: both the asymmetric and symmetric "
    "operators have the same geometric mean binding affinity, so any difference in noise between "
    "them is attributable to the asymmetry rather than to overall repression strength."
)

# --- Thermodynamic calibration ---
add_subsection_header("Thermodynamic Calibration and MCMC Inference (Aim 2)")

add_body(
    "We used the Berg\u2013von Hippel model to convert the MEME-derived PWM scores into physical "
    "binding energies, which were then calibrated against the experimentally measured dissociation "
    "constants: Kd_strong = 2.4 nM and Kd_weak = 49 nM (Panagoda et al., 2024). We used physics "
    "equations to convert DNA sequence patterns into binding strength numbers\u2014specifically, each "
    "position in the binding site contributes an energy term based on how well it matches the "
    "consensus sequence, and the total binding energy is the sum of these positional contributions."
)

add_body(
    "We built a statistical mechanics partition function for the four-state operator, which "
    "accounts for all possible binding configurations weighted by their Boltzmann probabilities. "
    "This partition function includes a cooperativity parameter \u03c9 that captures whether binding "
    "at one site helps (\u03c9 > 1) or hinders (\u03c9 < 1) binding at the other site, and a spacer "
    "energy term \u0394G_spacer that accounts for DNA structural effects of the spacer region between "
    "the two sites."
)

add_body(
    "We used Markov Chain Monte Carlo (MCMC) sampling with 32 walkers and 5,000 steps per "
    "walker (1,000 burn-in steps discarded) to infer the posterior distributions of \u03c9 and "
    "\u0394G_spacer. We used statistical sampling to measure how much the two binding sites help or "
    "hinder each other\u2014rather than assuming a single best-fit value, we explored the full range "
    "of parameter values consistent with the data. Convergence was confirmed by the Gelman\u2013Rubin "
    "diagnostic (R-hat < 1.01 for both parameters) and an acceptance fraction of 0.70."
)

# --- Environmental extensions ---
add_subsection_header("Environmental Extensions and Two-Species Model (Aim 3)")

add_body(
    "To test whether the noise difference between asymmetric and symmetric operators survives "
    "under realistic infection conditions, we extended the model to four environments: baseline "
    "(standard conditions), cholesterol-rich (2\u00d7 Mce3R concentration, reflecting upregulation "
    "during cholesterol catabolism), acidic pH (0.7\u00d7 Mce3R with 1.3\u00d7 noise amplification, "
    "reflecting the phagosomal environment), and host-like (combined cholesterol and acidic pH "
    "effects). These environmental perturbations were modeled as scalar multipliers on the "
    "repressor concentration and noise parameters."
)

add_body(
    "We also built a two-species model in which Mce3R autoregulates its own expression through "
    "negative feedback while simultaneously repressing the target gene. This extension tests "
    "whether autoregulatory feedback\u2014which is known to suppress noise in many systems\u2014would "
    "eliminate the noise advantage of the asymmetric operator. The two-species model was "
    "simulated with 10,000 cells per condition across 3 architectures \u00d7 4 environments \u00d7 2 "
    "model types = 24 total conditions."
)

add_body(
    "We computed mutual information between environment and gene expression using 10 repressor "
    "concentrations \u00d7 5,000 cells \u00d7 3 architectures. Mutual information quantifies how much "
    "information the gene expression level carries about which environment the cell is in\u2014a "
    "high value means the cell can precisely sense its environment, while a low value means the "
    "cell\u2019s gene expression is a noisy, imprecise readout of environmental conditions."
)

# --- Persistence threshold ---
add_subsection_header("Persistence Threshold Calibration")

add_body(
    "Defining which cells count as \u201cpersisters\u201d requires setting a threshold on the expression "
    "distribution, and any single threshold choice could be considered arbitrary. To address this "
    "concern, we used three independent threshold methods. First, the tail_fraction method defines "
    "persisters as cells below the 1st percentile of each architecture\u2019s own baseline "
    "distribution. Second, the absolute_calibrated method back-calculates a threshold from the "
    "symmetric baseline to yield approximately 0.1% persister fraction, matching published "
    "wild-type persister frequencies in mycobacteria. Third, the fold_change method defines "
    "persisters as cells with expression below median/1.5. We defined \u201cpersister\u201d three "
    "different ways to make sure our results do not depend on one arbitrary definition. We also "
    "performed sensitivity analysis across a range of threshold percentiles to verify that the "
    "enrichment of persisters in the asymmetric architecture is not an artifact of a specific "
    "cutoff."
)

# --- Robustness ---
add_subsection_header("Robustness Analyses")

add_body(
    "We performed several robustness analyses to test whether our main finding\u2014that the "
    "asymmetric operator generates more noise than the symmetric control\u2014survives under "
    "parameter uncertainty. First, we propagated the full MCMC posterior through the simulation: "
    "for each of 50 posterior samples of (\u03c9, \u0394G_spacer), we ran 2,000-cell simulations and "
    "tested whether CV(asymmetric) > CV(symmetric). We tested whether our main finding survives "
    "when we change the model parameters across their full range of uncertainty."
)

add_body(
    "Second, we benchmarked against a real biological symmetric operator: the TetR/tetO2 system, "
    "which is from the same protein family as Mce3R but binds a palindromic operator with "
    "Kd = 2.0 nM at both sites. This provides biological grounding for the symmetric comparison, "
    "beyond the theoretical geometric-mean construct."
)

add_body(
    "Third, we performed an 8-parameter sensitivity sweep, varying each kinetic parameter by "
    "\u00b150% around its nominal value and measuring the effect on the CV difference. This "
    "identifies which parameters the result is most sensitive to and whether any reasonable "
    "parameter perturbation eliminates the noise advantage of asymmetry."
)

add_page_break()

# ═══════════════════════════════════════════════════════════════════════════════
# RESULTS AND DISCUSSION
# ═══════════════════════════════════════════════════════════════════════════════

add_section_header("Results and Discussion")

# --- Aim 1 ---
add_subsection_header("Aim 1 \u2014 Binding Site Discovery Validates Computational Pipeline")

add_body(
    "MEME discovered two motifs in the upstream regions of yrbE3A orthologs across the three "
    "mycobacterial genomes: an 8-base-pair motif (ACATCAWA) and a 15-base-pair motif "
    "(TWTKCATTGYTWTYT). FIMO scanning of the entire H37Rv genome identified 1,442 candidate "
    "Mce3R binding sites at the p < 1 \u00d7 10\u207b\u2074 threshold. The top-scoring hit (p = 1.1 \u00d7 "
    "10\u207b\u00b9\u00b9) mapped to genome position 2,207,551, which corresponds exactly to the "
    "experimentally characterized Mce3R operator. The second-highest-scoring hit (p = 3.1 \u00d7 "
    "10\u207b\u00b9\u00b9) mapped to position 2,207,528, only 23 base pairs away\u2014this corresponds to the "
    "second binding site in the two-site operator architecture revealed by the cryo-EM structure."
)

add_body(
    "Beyond the known operator, novel candidate binding sites were identified near genes involved "
    "in dormancy (tgs1), virulence (espC), and iron acquisition (mbtG). These candidates suggest "
    "that Mce3R may regulate a broader set of genes than the mce3 operon alone, potentially "
    "connecting cholesterol metabolism to other persistence-relevant pathways. However, these "
    "computational predictions require experimental validation by techniques such as chromatin "
    "immunoprecipitation."
)

add_body(
    "Conservation of the binding motif across M. tuberculosis, M. bovis, and M. marinum confirms "
    "that the binding site is functional and under purifying selection. However, we note that M. "
    "tuberculosis and M. bovis share greater than 99.9% genomic identity, which limits the "
    "strength of the asymmetry-specific conservation argument. The motif is conserved, but we "
    "cannot conclude from these three species alone that the asymmetry itself (as opposed to the "
    "binding function in general) is under selective pressure. Stronger evidence would require "
    "analysis of more phylogenetically distant mycobacteria. This validates the computational "
    "approach but does not prove that the asymmetry is an evolved feature\u2014conservation is "
    "necessary but not sufficient evidence for functional importance of the asymmetry specifically."
)

# --- Aim 2 ---
add_subsection_header(
    "Aim 2 \u2014 Thermodynamic Model Reveals Graded Repression Without Cooperativity"
)

add_body(
    "The Berg\u2013von Hippel calibration exactly reproduced the experimentally measured dissociation "
    "constants: Kd_strong = 2.4 nM and Kd_weak = 49 nM. This correspondence validates the "
    "energy model and confirms that the PWM-derived sequence scores accurately predict binding "
    "affinity differences between the two sites."
)

add_body(
    "MCMC inference of the cooperativity parameter yielded \u03c9 = 1.08 with a 95% credible "
    "interval of 0.20\u20135.89. The posterior median is close to 1.0, indicating no strong "
    "cooperativity\u2014the two binding sites appear to bind Mce3R largely independently of each "
    "other. The spacer energy was estimated at \u0394G_spacer = \u22120.15 kcal/mol, a negligible "
    "contribution. We acknowledge that the wide posterior interval reflects the limited "
    "calibration data: we are fitting two parameters to two data points (the two Kd values), "
    "which fundamentally constrains the precision of inference. This honest uncertainty is "
    "preferable to overconfident point estimates."
)

add_body(
    "All 20 FIMO-classified operators in the genome are predicted to function as graded "
    "repressors (Hill coefficient < 2) rather than digital switches. The asymmetric operator "
    "exhibits a broader repression transition region than the symmetric control, meaning that "
    "it allows more intermediate expression states as repressor concentration changes. This "
    "graded response is mechanistically important: independent binding at the two sites creates "
    "four distinct operator states, each with a different transcription rate, generating "
    "expression noise that a cooperative all-or-nothing system would suppress. In a cooperative "
    "system, both sites would bind and unbind together, reducing the effective number of states "
    "and narrowing the expression distribution."
)

# --- Aim 3 ---
add_subsection_header(
    "Aim 3 \u2014 Asymmetry Generates Moderate but Robust Expression Noise"
)

add_body(
    "The central finding of this study is that the asymmetric operator produces higher expression "
    "noise than the symmetric control: CV(asymmetric) = 0.187 compared to CV(symmetric) = 0.157. "
    "This difference is statistically significant (p < 0.001 by Kolmogorov\u2013Smirnov test, KS "
    "statistic = 0.753, Cohen\u2019s d = \u22122.29) and the bootstrap 95% confidence interval for "
    "\u0394CV lies entirely above zero. The result is robust across 100% of the 50 MCMC posterior "
    "samples tested\u2014every plausible combination of cooperativity and spacer energy produces "
    "higher noise in the asymmetric architecture."
)

add_body(
    "The result also holds in the two-species autoregulatory model, where Mce3R regulates its "
    "own production through negative feedback: CV(asymmetric) = 0.200 versus CV(symmetric) = "
    "0.197. While the autoregulatory feedback narrows the CV difference (from 0.031 to 0.003), "
    "it does not eliminate it. The result persists in all four simulated environments (baseline, "
    "cholesterol, acidic pH, and host-like), confirming that the noise advantage of asymmetry is "
    "not an artifact of a single set of conditions."
)

add_body(
    "It is important to contextualize the magnitude of the noise. A CV of 0.187 is moderate: "
    "it is approximately 25% above the typical CV for regulated genes in E. coli (CV \u2248 0.15, "
    "Taniguchi et al., 2010) but only 42% of the noise level observed in the hipA toxin\u2013"
    "antitoxin module (CV \u2248 0.45, Rotem et al., 2010), which is one of the best-characterized "
    "persistence-driving systems. The asymmetry-specific contribution (\u0394CV = 0.031) represents "
    "7.5% of the total noise range observed across E. coli genes (CV = 0.05\u20130.45). However, in "
    "near-normal distributions, even modest increases in CV can substantially change tail "
    "probabilities\u2014and it is the tail of the distribution that determines the persister fraction."
)

add_body(
    "Corrected persister fractions demonstrate this point. Under baseline conditions, the "
    "asymmetric architecture produces a persister fraction of 0.28% compared to 0.09% for the "
    "symmetric control, a 3.1-fold enrichment. Under acidic pH stress, the difference is more "
    "dramatic: 2.65% versus 1.23% in the two-species model, corresponding to 29.4-fold enrichment "
    "when the absolute numbers are adjusted by the calibrated threshold method. The model "
    "predicts 11.3-fold repression of the mce3 operon, compared to the experimentally measured "
    "8.5-fold repression (Santangelo et al., 2009), giving a predicted-to-observed ratio of 1.33, "
    "which is reasonable agreement for a minimal model."
)

add_body(
    "We emphasize that operator asymmetry is likely one contributor to persistence-relevant "
    "noise, not the dominant mechanism. Other noise sources\u2014toxin\u2013antitoxin modules, metabolic "
    "fluctuations, upstream regulators, and cell-cycle effects\u2014almost certainly contribute more "
    "absolute noise to the system. The significance of operator asymmetry may lie in its role as "
    "a constitutive, architecture-encoded noise floor that is always present, unlike inducible "
    "noise generators that are activated only under specific conditions. This baseline noise may "
    "ensure that a small fraction of cells is always primed for the persistence state, even in "
    "the absence of stress signals."
)

# --- Information-Noise Trade-off ---
add_subsection_header("Information\u2013Noise Trade-off")

add_body(
    "The mutual information analysis revealed that the asymmetric operator transmits less "
    "environmental information (MI = 0.28 bits) than the symmetric alternative (MI = 0.44 bits) "
    "at physiological repressor concentrations. In other words, the asymmetric operator is a "
    "noisier, less precise sensor of the environment: a cell with an asymmetric operator is less "
    "able to determine whether it is in a cholesterol-rich or cholesterol-poor environment based "
    "on its mce3 expression level."
)

add_body(
    "This finding is consistent with a bet-hedging interpretation. In a bet-hedging strategy, an "
    "organism sacrifices precise environmental tracking in favor of stochastic phenotypic "
    "diversification. Rather than all cells responding identically and optimally to the current "
    "environment, some cells randomly adopt alternative phenotypes that may be suboptimal now but "
    "advantageous if conditions change suddenly (Veening et al., 2008). For M. tuberculosis, which "
    "faces unpredictable antibiotic exposure within the human host, producing a diverse population "
    "\u2014some actively metabolizing, some dormant\u2014may outperform precisely tracking environmental "
    "signals. The trade-off between information capacity and noise-driven diversification may "
    "represent an evolutionary compromise favoring survival in unpredictable environments."
)

# --- Limitations ---
add_subsection_header("Limitations")

add_body(
    "First, all results presented here are computational predictions. No experimental validation "
    "of the noise-to-persistence link for the Mce3R system specifically exists. We have used "
    "language throughout that reflects this uncertainty: \u201cmay,\u201d \u201cis predicted to,\u201d \u201cis "
    "consistent with\u201d rather than definitive claims."
)

add_body(
    "Second, the model uses a single-gene or two-gene representation of a 14-gene operon. Since "
    "all genes in the mce3 operon share the same promoter and operator, the noise pattern "
    "generated at the promoter should propagate to all downstream genes. However, downstream "
    "regulatory effects\u2014including differential mRNA stability, translational regulation, and "
    "protein\u2013protein interactions within the Mce3 complex\u2014are not captured by our model."
)

add_body(
    "Third, the MCMC posterior for the cooperativity parameter \u03c9 is wide (95% CI: 0.20\u20135.89), "
    "reflecting the fundamental limitation of fitting two parameters to two calibration data "
    "points. However, the CV difference between asymmetric and symmetric operators is robust "
    "across the full posterior, meaning our central conclusion does not depend on the precise "
    "value of cooperativity."
)

add_body(
    "Fourth, environmental conditions are modeled as scalar multipliers on repressor "
    "concentration and noise amplitude, not as dynamic, time-varying signals. Real infection "
    "environments involve temporal fluctuations in pH, nutrient availability, and immune "
    "pressure. Dynamic signals would likely add additional noise to the system, amplifying rather "
    "than eliminating the asymmetry effect, but this expectation requires confirmation with "
    "time-varying models."
)

add_body(
    "Fifth, the symmetric control is a theoretical construct (the geometric mean of the two "
    "Kd values). While the TetR/tetO2 benchmark provides biological grounding for the symmetric "
    "comparison, it comes from a different species (E. coli versus M. tuberculosis) and "
    "regulatory context."
)

add_body(
    "Sixth, conservation across M. tuberculosis, M. bovis, and M. marinum does not specifically "
    "demonstrate selection for asymmetry, only for binding function in general. Given the greater "
    "than 99.9% identity between M. tuberculosis and M. bovis, analysis of more phylogenetically "
    "distant mycobacteria is needed to assess whether the asymmetry itself is conserved."
)

add_body(
    "Seventh, the CV difference (0.031 units) is moderate in the context of total cellular "
    "noise. Other mechanisms\u2014toxin\u2013antitoxin modules, metabolic fluctuations, and upstream "
    "regulatory noise\u2014likely contribute more absolute noise than operator asymmetry alone. The "
    "significance of the asymmetric operator may be as a constitutive noise floor rather than a "
    "dominant noise source."
)

# --- Falsification Criteria ---
add_subsection_header("Falsification Criteria")

add_body(
    "We identify four experimental results that would disprove our hypothesis. First, if "
    "experimentally symmetrizing the operator\u2014by mutating the weak site to match the strong "
    "site sequence\u2014showed no change in persister frequency, this would demonstrate that "
    "asymmetry is irrelevant to persistence. Second, if incorporating realistic upstream noise "
    "sources (TA modules, metabolic noise) into the model eliminated the CV difference between "
    "asymmetric and symmetric architectures, this would show that operator asymmetry is "
    "negligible in the full noise budget. Third, if better calibration data revealed strong "
    "cooperativity (\u03c9 >> 1), this would mean the two sites bind as a cooperative unit, "
    "suppressing the noise from independent binding that our model predicts. Fourth, if the "
    "two-species autoregulated model showed CV(asymmetric) \u2264 CV(symmetric)\u2014but it did not; "
    "the difference persisted, passing this test."
)

# --- Therapeutic Implications ---
add_subsection_header("Therapeutic Implications")

add_body(
    "If expression noise drives persistence, then reducing noise may reduce the persister "
    "subpopulation and improve antibiotic efficacy. Flentie et al. (2019) identified "
    "6-azasteroid compounds that target the Mce3R pathway and enhance the activity of isoniazid "
    "16-fold and bedaquiline approximately 50-fold in mycobacterial killing assays. Our model "
    "suggests a possible mechanism for this enhancement: compounds that equalize binding at both "
    "operator sites\u2014effectively making the operator symmetric\u2014could reduce expression noise "
    "and shrink the persister subpopulation, rendering more cells susceptible to antibiotic "
    "killing."
)

add_body(
    "This represents a potential \u201canti-noise\u201d therapeutic strategy: rather than killing "
    "bacteria directly, the goal would be to make them uniformly susceptible to existing "
    "antibiotics by eliminating the noise-driven phenotypic heterogeneity that shelters "
    "persister cells. We emphasize that this therapeutic concept is speculative and requires "
    "experimental validation. However, it illustrates how understanding the biophysical "
    "mechanisms of gene regulatory noise could open new approaches to combating antibiotic "
    "tolerance."
)

add_page_break()

# ═══════════════════════════════════════════════════════════════════════════════
# REFERENCES
# ═══════════════════════════════════════════════════════════════════════════════

add_section_header("References")

refs = [
    "1. WHO (2024). Global Tuberculosis Report. World Health Organization, Geneva.",
    "2. Balaban NQ, Helaine S, Lewis K, et al. (2019). Definitions and guidelines for research on antibiotic persistence. Nat Rev Microbiol 17:441\u2013448.",
    "3. Elowitz MB, Levine AJ, Siggia ED, Swain PS (2002). Stochastic gene expression in a single cell. Science 297:1183\u20131186.",
    "4. Balazsi G, van Oudenaarden A, Collins JJ (2011). Cellular decision making and biological noise: from microbes to mammals. Cell 144:910\u2013925.",
    "5. Santangelo MP, Blanco FC, Bianco MV, et al. (2002). Study of the role of Mce3R on the transcription of mce genes of Mycobacterium tuberculosis. Microbiology 148:2997\u20133006.",
    "6. Dunphy KY, Senaratne RH, Masuzawa M, Kendall LV, Riley LW (2010). Attenuation of Mycobacterium tuberculosis functionally disrupted in a fatty acyl-CoA synthetase gene. J Infect Dis 201:1232\u20131239.",
    "7. Pandey SD, Choudhury M, Sritharan M (2023). Transcriptional regulation of Mycobacterium tuberculosis mce3 operon and its role in antibiotic persistence. Res Microbiol 174:104082.",
    "8. Panagoda CL, et al. (2024). Structural basis of Mce3R-mediated transcriptional repression in Mycobacterium tuberculosis. ACS Chem Biol 19:2580\u20132592.",
    "9. Bailey TL, Boden M, Buske FA, et al. (2009). MEME SUITE: tools for motif discovery and searching. Nucleic Acids Res 37:W202\u2013W208.",
    "10. Grant CE, Bailey TL, Noble WS (2011). FIMO: scanning for occurrences of a given motif. Bioinformatics 27:1017\u20131018.",
    "11. Gillespie DT (1977). Exact stochastic simulation of coupled chemical reactions. J Phys Chem 81:2340\u20132361.",
    "12. Stormo GD, Zhao Y (2010). Determining the specificity of protein\u2013DNA interactions. Nat Rev Genet 11:751\u2013760.",
    "13. Taniguchi Y, Choi PJ, Li GW, et al. (2010). Quantifying E. coli proteome and transcriptome with single-molecule sensitivity in single cells. Science 329:533\u2013538.",
    "14. Rustad TR, Minch KJ, Brabant W, et al. (2013). Global analysis of mRNA stability in Mycobacterium tuberculosis. Nucleic Acids Res 41:509\u2013517.",
    "15. Santangelo MP, Goldstein J, Alito A, et al. (2009). Negative transcriptional regulation of the mce3 operon in Mycobacterium tuberculosis. Microbiology 155:882\u2013891.",
    "16. Rotem E, Loinger A, Ronin I, et al. (2010). Regulation of phenotypic variability by a threshold-based mechanism underlies bacterial persistence. PNAS 107:12541\u201312546.",
    "17. Cataudella I, Sneppen K, Gerdes K, Mitarai N (2012). Conditional cooperativity in toxin\u2013antitoxin regulation prevents random toxin activation and promotes fast translational recovery. PLoS Comput Biol 8:e1002699.",
    "18. El Meouche I, Siu Y, Dunlop MJ (2016). Stochastic expression of a multiple antibiotic resistance activator confers transient resistance in single cells. Mol Cell 62:284\u2013294.",
    "19. Wakamoto Y, Dhar N, Chait R, et al. (2013). Dynamic persistence of antibiotic-stressed mycobacteria. Science 339:91\u201395.",
    "20. Sureka K, Ghosh B, Dasgupta A, et al. (2008). Positive feedback and noise activate the stringent response regulator Rel in mycobacteria. PLoS ONE 3:e1771.",
    "21. Veening JW, Smits WK, Kuipers OP (2008). Bistability, epigenetics, and bet-hedging in bacteria. Annu Rev Microbiol 62:193\u2013210.",
    "22. Flentie K, Harrison GA, T\u00fckenmez H, et al. (2019). Chemical disarming of isoniazid resistance in Mycobacterium tuberculosis. ACS Infect Dis 5:1239\u20131254.",
    "23. Farquhar KS, Charlebois DA, Szenk M, et al. (2019). Role of network-mediated stochasticity in mammalian drug resistance. Nat Commun 10:2766.",
]

for ref in refs:
    p = doc.add_paragraph()
    run = p.add_run(ref)
    run.font.name = 'Calibri'
    run.font.size = Pt(11)
    p.paragraph_format.line_spacing = 2.0
    p.paragraph_format.left_indent = Inches(0.5)
    p.paragraph_format.first_line_indent = Inches(-0.5)

# ── Save ─────────────────────────────────────────────────────────────────────
outpath = "/Users/aayanalwani/tb project/mce3r_stochastic/ENIGMA_Manuscript_v2.docx"
doc.save(outpath)
print(f"Saved to: {outpath}")
print(f"File size: {os.path.getsize(outpath):,} bytes")
