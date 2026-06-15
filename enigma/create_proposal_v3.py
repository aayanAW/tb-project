#!/usr/bin/env python3
"""Generate ENIGMA Project Proposal v3 as a comprehensive DOCX document."""

from docx import Document
from docx.shared import Pt, Inches, Cm, RGBColor
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.enum.section import WD_ORIENT
from docx.oxml.ns import qn
from docx.oxml import OxmlElement
import os

OUTPUT_PATH = "/Users/aayanalwani/tb project/mce3r_stochastic/ENIGMA_Project_Proposal_v3.docx"

doc = Document()

# ── Global defaults ──────────────────────────────────────────────────────────
style = doc.styles['Normal']
font = style.font
font.name = 'Calibri'
font.size = Pt(11)
font.color.rgb = RGBColor(0, 0, 0)
pf = style.paragraph_format
pf.line_spacing = 2.0
pf.space_after = Pt(0)
pf.space_before = Pt(0)

# Set 1-inch margins on all sections
for section in doc.sections:
    section.top_margin = Inches(1)
    section.bottom_margin = Inches(1)
    section.left_margin = Inches(1)
    section.right_margin = Inches(1)

# ── Helper functions ─────────────────────────────────────────────────────────

def set_font(run, name='Calibri', size=11, bold=False, italic=False, color=None):
    run.font.name = name
    run.font.size = Pt(size)
    run.bold = bold
    run.italic = italic
    if color:
        run.font.color.rgb = RGBColor(*color)

def add_title_text(text, size=16, bold=True, alignment=WD_ALIGN_PARAGRAPH.CENTER):
    p = doc.add_paragraph()
    p.alignment = alignment
    p.paragraph_format.line_spacing = 2.0
    run = p.add_run(text)
    set_font(run, size=size, bold=bold)
    return p

def add_section_header(text, size=14):
    p = doc.add_paragraph()
    p.paragraph_format.space_before = Pt(12)
    p.paragraph_format.line_spacing = 2.0
    run = p.add_run(text)
    set_font(run, size=size, bold=True)
    return p

def add_subsection_header(text, size=12):
    p = doc.add_paragraph()
    p.paragraph_format.space_before = Pt(8)
    p.paragraph_format.line_spacing = 2.0
    run = p.add_run(text)
    set_font(run, size=size, bold=True)
    return p

def add_body(text):
    p = doc.add_paragraph()
    p.paragraph_format.line_spacing = 2.0
    run = p.add_run(text)
    set_font(run)
    return p

def add_body_italic(text):
    p = doc.add_paragraph()
    p.paragraph_format.line_spacing = 2.0
    run = p.add_run(text)
    set_font(run, italic=True)
    return p

def add_bullet(text):
    p = doc.add_paragraph(style='List Bullet')
    p.paragraph_format.line_spacing = 2.0
    # Clear default run and add formatted one
    p.clear()
    run = p.add_run(text)
    set_font(run)
    return p

def add_page_break():
    doc.add_page_break()

# ══════════════════════════════════════════════════════════════════════════════
# TITLE PAGE
# ══════════════════════════════════════════════════════════════════════════════

# Add some vertical space
for _ in range(6):
    add_body("")

add_title_text("Temporal Noise Dynamics and Therapeutic Noise Quenching\nin the Asymmetric Mce3R Operator of\nMycobacterium tuberculosis", size=16, bold=True)

add_body("")

add_title_text("ENIGMA: Expression Noise In Gene-regulatory Mechanisms and Architecture", size=13, bold=False)

add_body("")
add_body("")

add_title_text("[Author Name]", size=12, bold=False)
add_title_text("March 2026", size=12, bold=False)

add_page_break()

# ══════════════════════════════════════════════════════════════════════════════
# ABSTRACT
# ══════════════════════════════════════════════════════════════════════════════

add_section_header("ABSTRACT")

add_body(
    "Mycobacterium tuberculosis kills over 1.2 million people annually, with treatment requiring "
    "6\u20139 months of daily antibiotics due to phenotypically tolerant persister cells. Antibiotic "
    "persistence arises from stochastic gene expression noise \u2014 random cell-to-cell variation that "
    "pushes a subpopulation into a drug-tolerant dormant state. While it is established that "
    "cis-regulatory architecture shapes noise (Chowdhury et al. 2021) and that noise drives "
    "M. tuberculosis persistence (Quigley & Lewis 2022), the specific structural features of "
    "regulatory DNA that generate persistence-relevant noise remain uncharacterized for most systems."
)

add_body(
    "The Mce3R transcriptional repressor controls cholesterol metabolism genes critical for host "
    "survival and binds an unusual asymmetric operator with a 20.4-fold binding affinity difference "
    "between its two sites (Kd = 2.4 nM vs. 49 nM; Panagoda et al. 2024). Deletion of mce3R "
    "increases persister frequency (Pandey et al. 2023), but no study has mechanistically connected "
    "the operator\u2019s structural asymmetry to noise generation."
)

add_body(
    "We present the first comprehensive noise characterization of a structurally resolved asymmetric "
    "bacterial operator. In Aim 1, we computationally discover and validate Mce3R binding sites "
    "across three mycobacterial genomes using de novo motif discovery (MEME/FIMO), identifying "
    "1,442 candidate sites genome-wide. In Aim 2, we construct a thermodynamically calibrated "
    "four-state stochastic model using Gillespie simulation, Berg\u2013von Hippel energy calibration, "
    "and MCMC cooperativity inference (\u03c9 = 1.08, 95% CI: 0.20\u20135.89). In Aim 3, we quantify "
    "noise amplitude (CV = 0.187 asymmetric vs. 0.157 symmetric, p < 0.001), characterize temporal "
    "noise dynamics \u2014 including autocorrelation time and dwell time in the persister state \u2014 and "
    "predict the dose-response of noise quenching through systematic operator symmetrization."
)

add_body(
    "Our results indicate that the asymmetric operator generates 19% more expression noise than an "
    "equivalent symmetric architecture, that this difference is robust across autoregulatory feedback, "
    "four environmental conditions, and the full MCMC posterior, and that the noise exhibits "
    "architecture-dependent temporal properties. We provide quantitative predictions for how much "
    "symmetrization would be required to reduce the persister subpopulation, identifying a potential "
    "\u201canti-noise\u201d therapeutic strategy. All results are computational predictions requiring "
    "experimental validation, and explicit falsification criteria are provided."
)

add_page_break()

# ══════════════════════════════════════════════════════════════════════════════
# 1. INTRODUCTION AND SIGNIFICANCE
# ══════════════════════════════════════════════════════════════════════════════

add_section_header("1. INTRODUCTION AND SIGNIFICANCE")

# 1.1
add_subsection_header("1.1 The Problem \u2014 Tuberculosis Persistence")

add_body(
    "Tuberculosis remains the deadliest infectious disease in human history and continues to exact "
    "a staggering toll: the World Health Organization reported that Mycobacterium tuberculosis killed "
    "over 1.2 million people in 2023, with an additional 10.6 million new cases diagnosed worldwide "
    "(WHO 2024). Despite the availability of effective antibiotics, standard treatment requires 6\u20139 "
    "months of daily multidrug therapy \u2014 a regimen so prolonged that patient non-compliance is a "
    "leading cause of treatment failure and the emergence of drug-resistant strains."
)

add_body(
    "The fundamental question driving this proposal is: why does treatment take so long? The answer "
    "lies not in genetic resistance but in phenotypic tolerance. Within any M. tuberculosis "
    "population, a small subpopulation of cells enters a dormant, non-replicating state in which "
    "they become transiently tolerant to antibiotics. These cells, termed persisters, are genetically "
    "identical to their drug-sensitive siblings but have stochastically adopted a phenotype that "
    "renders them insensitive to drugs that target active cellular processes (Balaban et al. 2019). "
    "When antibiotic pressure is removed, persisters can resume growth and regenerate a fully "
    "drug-sensitive population, driving clinical relapse."
)

add_body(
    "Critically, it is not merely the frequency of persisters that determines treatment duration "
    "but also the duration of their dormancy. A population in which 0.1% of cells become dormant "
    "for one day requires much shorter treatment than one in which 0.1% become dormant for weeks "
    "or months. The temporal dynamics of persistence \u2014 how long individual cells remain in the "
    "tolerant state and how frequently they switch between active and dormant phenotypes \u2014 are "
    "therefore directly relevant to clinical outcomes. Current frontline antibiotics cannot "
    "distinguish active from dormant cells, and no approved therapeutic strategy specifically "
    "targets the molecular mechanisms that generate and maintain the persister state."
)

# 1.2
add_subsection_header("1.2 Noise and Persistence \u2014 What Is Already Known")

add_body(
    "Gene expression is inherently stochastic. Because transcription factors, RNA polymerases, "
    "and ribosomes exist at low copy numbers inside cells, the biochemical reactions governing gene "
    "expression are subject to random molecular collisions and thermal fluctuations. This produces "
    "cell-to-cell variation in protein levels even among genetically identical cells grown in "
    "identical environments \u2014 a phenomenon known as gene expression noise (Elowitz et al. 2002). "
    "Noise is not merely a nuisance; it is a fundamental property of biological systems that has "
    "been co-opted by evolution for adaptive purposes, including phenotypic diversification and "
    "bet-hedging strategies (Balazsi et al. 2011)."
)

add_body(
    "The architecture of cis-regulatory DNA \u2014 the number, arrangement, affinity, and spacing of "
    "transcription factor binding sites \u2014 plays a critical role in shaping noise amplitude and "
    "dynamics. Chowdhury et al. (2021) demonstrated that cis-regulatory logic produces gene "
    "expression noise patterns that describe phenotypic heterogeneity in bacteria, establishing "
    "a direct link between promoter structure and noise output. Lengyel and Morelli (2017) showed "
    "that multiple binding sites for transcriptional repressors can produce either regular bursting "
    "or noise suppression depending on their configuration, demonstrating that the specific "
    "arrangement of binding sites matters profoundly. Gene regulatory circuit architecture controls "
    "noise-driven cellular decision-making across organisms, from bacteria to mammals (Balazsi "
    "et al. 2011; Farquhar et al. 2019)."
)

add_body(
    "In M. tuberculosis specifically, Quigley and Lewis (2022) made a landmark discovery: noise "
    "in the metabolic gene ackA creates persister cells, and overexpressing ackA to quench "
    "expression noise reduces the persister frequency. This established the direct causal link "
    "between gene expression noise and antibiotic persistence in the TB pathogen. Additional "
    "studies have shown that noise in toxin-antitoxin modules (Rotem et al. 2010), multidrug "
    "resistance activators (El Meouche et al. 2016), and mycobacterial catalase-peroxidase "
    "(Wakamoto et al. 2013) contributes to phenotypic heterogeneity and persistence. The stringent "
    "response regulator Rel in mycobacteria is activated by positive feedback and noise, further "
    "connecting stochastic gene expression to the persister phenotype (Sureka et al. 2008)."
)

add_body(
    "These studies establish the general principle that cis-regulatory architecture shapes noise "
    "and that noise drives persistence. What remains unknown is how specific architectural features "
    "of specific operators contribute to persistence-relevant noise in defined regulatory systems "
    "with experimentally measured biophysical parameters."
)

# 1.3
add_subsection_header("1.3 The Mce3R System \u2014 An Unprecedented Asymmetric Operator")

add_body(
    "Mce3R (Rv1963c) is a TetR-family transcriptional repressor in M. tuberculosis that controls "
    "the mce3 operon (Rv1964\u2013Rv1977), a 14-gene cluster encoding cholesterol and lipid transport "
    "machinery essential for survival within host macrophages (Santangelo et al. 2002; Dunphy "
    "et al. 2010). The mce3 operon is required for cholesterol import during intracellular "
    "infection, and its regulation by Mce3R is critical for metabolic adaptation to the host "
    "environment. Santangelo et al. (2009) demonstrated that Mce3R directly represses transcription "
    "of the mce3 genes, confirming its role as the primary regulatory switch for this metabolic "
    "pathway."
)

add_body(
    "The direct link between Mce3R and antibiotic persistence was established by Pandey et al. "
    "(2023), who showed that deletion of mce3R increases antibiotic persister frequency. This "
    "finding places Mce3R squarely at the intersection of metabolic regulation and persistence "
    "biology, making it a compelling target for mechanistic investigation."
)

add_body(
    "A transformative advance came with the recent cryo-EM structural characterization of Mce3R "
    "bound to its operator DNA by Panagoda et al. (2024, PDB 9B7Y). This study revealed several "
    "structurally unprecedented features. First, Mce3R is a double TetR-fold repeat protein, a "
    "unique architecture not previously observed in the TetR family. Second, and most importantly "
    "for this proposal, Mce3R binds an asymmetric 123 bp operator containing two binding sites "
    "with dramatically different affinities: a strong downstream site with Kd = 2.4 \u00b1 0.7 nM "
    "and a weak upstream site with Kd \u2248 49 nM, separated by a 53 bp spacer. This represents a "
    "20.4-fold affinity difference between the two sites \u2014 the largest documented asymmetry in "
    "any TetR-family regulatory system."
)

add_body(
    "Most TetR-family repressors bind palindromic (symmetric) operators, reflecting their "
    "homodimeric structure. The pronounced asymmetry of the Mce3R operator is therefore exceptional "
    "and raises immediate questions about its functional consequences. From a biophysical "
    "perspective, two sites with very different affinities will be occupied at different repressor "
    "concentrations, creating a more complex pattern of operator states than a symmetric system. "
    "This complexity, we hypothesize, has direct consequences for gene expression noise."
)

# 1.4
add_subsection_header("1.4 The Gap")

add_body(
    "Despite the convergence of structural, biophysical, and phenotypic data on the Mce3R system, "
    "several critical gaps remain:"
)

add_bullet(
    "No study has quantitatively modeled the noise output of the Mce3R operator using its "
    "experimentally measured Kd values. The structural data from Panagoda et al. (2024) provide "
    "the biophysical parameters needed for such modeling, but no such analysis has been performed."
)
add_bullet(
    "No study has characterized the temporal dynamics of noise from an asymmetric operator \u2014 "
    "including autocorrelation time, dwell time in the persister state, and switching frequency "
    "between active and dormant phenotypes."
)
add_bullet(
    "No study has predicted the dose-response relationship of noise quenching for any specific "
    "TB regulatory system. While Quigley and Lewis (2022) demonstrated that noise quenching "
    "reduces persistence, the quantitative relationship between the degree of noise reduction "
    "and the magnitude of persistence reduction remains unknown."
)
add_bullet(
    "The mechanistic connection between Mce3R\u2019s structural asymmetry and the persistence "
    "phenotype observed by Pandey et al. (2023) has not been explained. Why does the operator "
    "have asymmetric binding sites, and what are the functional consequences of this asymmetry?"
)

# 1.5
add_subsection_header("1.5 This Study")

add_body(
    "We address these gaps through an integrated computational pipeline that combines de novo "
    "motif discovery, thermodynamic modeling, stochastic simulation, Bayesian inference, and "
    "information-theoretic analysis. We provide the first comprehensive noise characterization "
    "\u2014 encompassing amplitude, temporal dynamics, and quenching response \u2014 of a structurally "
    "resolved asymmetric bacterial operator with experimentally measured binding affinities."
)

add_body(
    "Our approach is grounded in the experimentally determined biophysical parameters of the "
    "Mce3R system (Panagoda et al. 2024) and generates quantitative, testable predictions that "
    "bridge structural biology and persistence microbiology. We emphasize throughout that all "
    "results are computational predictions requiring experimental validation, and we provide "
    "explicit falsification criteria specifying what experimental outcomes would disprove our "
    "central hypotheses."
)

add_page_break()

# ══════════════════════════════════════════════════════════════════════════════
# 2. SPECIFIC AIMS
# ══════════════════════════════════════════════════════════════════════════════

add_section_header("2. SPECIFIC AIMS")

# AIM 1
add_subsection_header("Aim 1: Computationally Identify and Validate Mce3R Binding Sites Across Mycobacterial Genomes")

add_body_italic(
    "Question: Can we rediscover the Mce3R binding motif de novo and identify candidate "
    "regulatory targets genome-wide?"
)

add_body(
    "Approach: We download three complete mycobacterial genomes from NCBI: M. tuberculosis H37Rv "
    "(NC_000962.3), M. bovis AF2122/97 (NC_002945.4), and M. marinum M (NC_010612.1). From each "
    "genome, we extract 200 bp upstream of yrbE3A orthologs, the first gene in the mce3 operon "
    "and the known Mce3R target. We run MEME de novo motif discovery in ZOOPS mode (Zero Or One "
    "Occurrence Per Sequence), searching both strands with motif widths of 6\u2013110 bp and a "
    "0-order Markov background model parameterized with H37Rv GC-content nucleotide frequencies "
    "(A = 0.204, C = 0.296, G = 0.296, T = 0.204). The resulting position weight matrices (PWMs) "
    "are used to scan the entire H37Rv genome with FIMO at a p-value threshold of 1e\u22124, "
    "identifying all candidate Mce3R binding sites genome-wide. We assess cross-species "
    "conservation using sliding-window sequence comparison across the three genomes."
)

add_body(
    "This aim serves two purposes: it validates the computational pipeline by testing whether "
    "de novo motif discovery can independently recover the binding motif identified by structural "
    "biology, and it provides the PWM required for downstream thermodynamic energy calibration "
    "in Aim 2. The genome-wide scan additionally identifies novel candidate regulatory targets "
    "that may expand our understanding of the Mce3R regulon."
)

# AIM 2
add_subsection_header("Aim 2: Construct a Thermodynamically Calibrated Stochastic Model of the Mce3R Operator")

add_body_italic(
    "Question: How does operator architecture quantitatively control repression strength and "
    "expression noise?"
)

add_body(
    "Sub-aim 2a \u2014 Four-state operator model: We model the Mce3R operator as existing in four "
    "states: unbound, strong-site-only bound, weak-site-only bound, and both-sites bound. Each "
    "state has a distinct transcription rate reflecting the degree of repression. We implement "
    "12 chemical reactions in a Gillespie stochastic simulation algorithm: 8 operator state "
    "transitions (binding and unbinding at each site from each accessible state), plus "
    "transcription, translation, mRNA degradation, and protein degradation. We simulate 50,000 "
    "cells per condition for 30,000 minutes with a 15,000-minute burn-in period. We compare four "
    "architectures: (A) the native asymmetric operator, (B) a symmetric geometric-mean control "
    "with equivalent total regulatory capacity, (C) a single-site operator, and (D) an "
    "unregulated constitutive promoter."
)

add_body(
    "Sub-aim 2b \u2014 Thermodynamic calibration: We use the Berg\u2013von Hippel biophysical model to "
    "convert MEME-derived PWM scores into per-position binding energies. The model computes "
    "the free energy contribution of each nucleotide position as \u0394G_i = \u2013kT \u00d7 ln(f_i / p_i), "
    "where f_i is the PWM frequency at position i and p_i is the genomic background frequency. "
    "We calibrate the energy scale so that the total binding energy reproduces the experimentally "
    "measured Kd values of 2.4 nM for the strong site and 49 nM for the weak site. We build a "
    "statistical mechanics partition function for the four-state system that relates binding "
    "energies to equilibrium occupancy probabilities."
)

add_body(
    "Sub-aim 2c \u2014 MCMC cooperativity inference: We use Markov Chain Monte Carlo sampling via the "
    "emcee ensemble sampler (32 walkers \u00d7 5,000 steps with 1,000-step burn-in) to infer two "
    "parameters from simulation data: the cooperativity parameter \u03c9 (which quantifies whether "
    "binding at one site enhances or inhibits binding at the other) and the spacer energy penalty "
    "\u0394G_spacer (which captures the energetic cost of the 53 bp spacer between binding sites). "
    "We report full posterior distributions and convergence diagnostics (R-hat, acceptance "
    "fraction, trace plots)."
)

add_body(
    "Sub-aim 2d \u2014 Operator classification: Using the fitted partition function, we classify all "
    "FIMO-predicted operators as digital switches (Hill coefficient n > 2, exhibiting sharp "
    "on/off transitions) or graded repressors (Hill coefficient n < 2, exhibiting continuous "
    "dose-response). This classification has functional implications: digital switches produce "
    "bimodal expression distributions, while graded repressors produce unimodal distributions "
    "with tunable noise."
)

# AIM 3
add_subsection_header("Aim 3: Characterize Noise Amplitude, Temporal Dynamics, and Therapeutic Quenching Response")

add_body_italic(
    "Question: What are the full noise characteristics of the asymmetric operator, and how much "
    "symmetrization would be needed to eliminate the persistence advantage?"
)

add_body(
    "Sub-aim 3a \u2014 Noise amplitude and robustness: We quantify the coefficient of variation (CV) "
    "of protein expression across the four architectures and four environmental conditions "
    "(baseline, cholesterol exposure, acidic pH, and combined host-like conditions). We extend "
    "the model to a two-species autoregulatory system in which Mce3R regulates both itself "
    "(negative feedback) and the target gene, testing whether autoregulation preserves or "
    "eliminates the noise difference. We propagate uncertainty from the MCMC posterior by "
    "drawing 50 parameter samples and running independent Gillespie simulations for each, "
    "reporting the distribution of CV differences. We validate the model framework against "
    "a TetR/tetO2 biological benchmark from E. coli."
)

add_body(
    "Sub-aim 3b \u2014 Temporal noise dynamics: This represents a key novel contribution. We record "
    "full protein expression time traces from Gillespie simulations at 10-minute intervals after "
    "the burn-in period. From these traces, we compute the autocorrelation function C(\u03c4) = "
    "\u27e8\u03b4x(t) \u00b7 \u03b4x(t+\u03c4)\u27e9 / \u27e8\u03b4x\u00b2\u27e9 using an FFT-based method, and fit an exponential "
    "decay C(\u03c4) ~ exp(\u2013\u03c4/\u03c4_c) to extract the autocorrelation time \u03c4_c, which quantifies how "
    "long noise fluctuations persist (noise memory). We compute dwell time distributions in the "
    "persister state \u2014 contiguous intervals where protein count falls below a persistence "
    "threshold \u2014 for each architecture. We compute the power spectral density S(f) to identify "
    "characteristic switching frequencies. We compare \u03c4_c, mean dwell time, and spectral "
    "properties across the asymmetric, symmetric, and single-site architectures, and test "
    "whether asymmetry increases noise memory (longer \u03c4_c), not just noise amplitude (higher CV)."
)

add_body(
    "Sub-aim 3c \u2014 Noise quenching dose-response: This is the second key novel contribution. We "
    "systematically symmetrize the operator by sweeping Kd_weak from 49 nM toward 2.4 nM in 10 "
    "steps while holding Kd_strong fixed, monitoring CV, persister fraction, and \u03c4_c at each "
    "step. We separately sweep Mce3R concentration from 10 to 10,000 nM in 15 log-spaced steps "
    "for both asymmetric and symmetric architectures. We compute the IC50 of noise \u2014 the degree "
    "of symmetrization required to reduce architecture-driven noise by 50% \u2014 and the IC50 of "
    "persistence, the degree of symmetrization required to reduce the persister fraction by 50%. "
    "These computations generate a quantitative dose-response curve predicting how much operator "
    "modification would be needed to achieve therapeutic noise reduction."
)

add_body(
    "Sub-aim 3d \u2014 Mutual information analysis: We calculate the mutual information "
    "I(environment; expression) for each architecture across 10 Mce3R concentrations, testing "
    "whether the asymmetric operator trades information capacity for phenotypic diversification. "
    "If the asymmetric operator transmits less environmental information than the symmetric "
    "alternative, this would suggest that the architecture is optimized for bet-hedging rather "
    "than faithful signal transduction."
)

add_page_break()

# ══════════════════════════════════════════════════════════════════════════════
# 3. METHODS
# ══════════════════════════════════════════════════════════════════════════════

add_section_header("3. METHODS")

# 3.1
add_subsection_header("3.1 Genome Acquisition and Sequence Extraction")

add_body(
    "Three complete mycobacterial genomes were downloaded from the National Center for "
    "Biotechnology Information (NCBI) using BioPython (Cock et al. 2009): M. tuberculosis H37Rv "
    "(GenBank accession NC_000962.3; 4,411,532 bp; 65.6% GC content), M. bovis AF2122/97 "
    "(NC_002945.4), and M. marinum M (NC_010612.1). These three species were selected to span "
    "a range of evolutionary divergence within the Mycobacterium genus while ensuring the presence "
    "of mce3 operon orthologs."
)

add_body(
    "For each genome, we extracted 200 bp upstream of the yrbE3A ortholog, the first gene in the "
    "mce3 operon and the primary transcriptional target of Mce3R. The known operator sequence is "
    "123 bp in length (as resolved in the cryo-EM structure, PDB 9B7Y) with the following sequence: "
    "GCCCCGCGCTATAGGATACTAGCAAGATACATCATAGCCAATATATGCCAGTTTGCATTGCTATTTACCGATCAGTTGTCCAAGCAA"
    "TCGCGTATTGGCTATGGACATCAGCGGTTCTGCCGC. This operator is located at H37Rv genomic coordinates "
    "2,207,477\u20132,207,699 in the intergenic region between mce3R (Rv1963c) and yrbE3A (Rv1964)."
)

# 3.2
add_subsection_header("3.2 De Novo Motif Discovery (MEME)")

add_body(
    "De novo motif discovery was performed using MEME Suite version 5.5.9 (Bailey et al. 2009). "
    "The input consisted of three orthologous upstream sequences (200 bp each) extracted from the "
    "three mycobacterial genomes. MEME was run in ZOOPS mode (Zero Or One Occurrence Per Sequence), "
    "which assumes each motif appears at most once in each input sequence. Both strands were "
    "searched (+ and \u2013), with motif width range set to 6\u2013110 bp. The background model was a "
    "0-order Markov model parameterized with H37Rv nucleotide frequencies (A = 0.204, C = 0.296, "
    "G = 0.296, T = 0.204), reflecting the high GC content characteristic of mycobacterial genomes."
)

add_body(
    "The top 5 motifs were retained, ranked by E-value, and reported as position weight matrices "
    "(PWMs). The purpose of this analysis was to discover the Mce3R binding motif without prior "
    "structural information, thereby validating that the computational approach can independently "
    "recover what cryo-EM structural biology identified. Successful de novo recovery of the known "
    "binding motif would provide confidence in the computational pipeline for downstream analyses."
)

# 3.3
add_subsection_header("3.3 Genome-Wide Scanning (FIMO)")

add_body(
    "The complete M. tuberculosis H37Rv genome was scanned using FIMO (Find Individual Motif "
    "Occurrences; Grant et al. 2011) with all 5 MEME-derived PWMs. The p-value threshold was set "
    "to 1e\u22124. Each hit was reported with its genomic coordinates, score in bits, p-value, matched "
    "sequence, nearest gene annotation, and intergenic status. The purpose of this genome-wide "
    "scan was to identify all candidate Mce3R binding sites, enabling both validation of the known "
    "operator location and discovery of novel regulatory targets."
)

# 3.4
add_subsection_header("3.4 Cross-Species Conservation Analysis")

add_body(
    "Cross-species conservation of binding site regions was assessed using a sliding-window "
    "sequence comparison across all three mycobacterial genomes. A NumPy-vectorized rolling match "
    "algorithm (using np.frombuffer and rolling count operations) compared orthologous binding site "
    "regions with a conservation threshold of >80% nucleotide identity. Regions exceeding this "
    "threshold were flagged as conserved. The purpose was to assess whether putative binding sites "
    "are preserved across mycobacterial species, which would support their functional importance."
)

# 3.5
add_subsection_header("3.5 Four-State Operator Model")

add_body(
    "The Mce3R operator was modeled as a system with four distinct states, each associated with "
    "a specific transcription rate reflecting the degree of repression:"
)

add_bullet("State 0 (unbound): k_txn = 0.150 mRNA/min (full transcription)")
add_bullet("State 1 (strong site bound only): k_txn = 0.0225 mRNA/min (85% transcription blocked)")
add_bullet("State 2 (weak site bound only): k_txn = 0.075 mRNA/min (50% transcription blocked)")
add_bullet("State 3 (both sites bound): k_txn = 0.01125 mRNA/min (92.5% transcription blocked)")

add_body(
    "The model comprises 12 chemical reactions: 8 operator state transitions (binding and "
    "unbinding of Mce3R at each site from each accessible operator state), plus mRNA transcription, "
    "protein translation, mRNA degradation, and protein degradation. Transition rates for binding "
    "were computed as k_on = 0.0167 nM\u207b\u00b9 min\u207b\u00b9 (Stormo & Zhao 2010), and unbinding rates were "
    "computed from the relationship k_off = Kd \u00d7 k_on for each site."
)

# 3.6
add_subsection_header("3.6 Gillespie Stochastic Simulation Algorithm")

add_body(
    "Stochastic simulations were performed using the exact Stochastic Simulation Algorithm (SSA) "
    "of Gillespie (1977), implemented in Python with Numba JIT compilation for computational "
    "efficiency. At each simulation step, all 12 reaction propensities were computed, a waiting "
    "time \u03c4 was drawn from an exponential distribution with rate parameter equal to the sum of "
    "all propensities (a_total), a single reaction was selected by cumulative propensity sampling, "
    "and the molecular counts were updated accordingly."
)

add_body("Kinetic parameters were drawn from published experimental measurements:")

add_bullet(
    "Translation rate: k_translation = 0.5 protein/mRNA/min (adapted from Taniguchi et al. 2010, "
    "adjusted for the slower growth rate of M. tuberculosis)"
)
add_bullet(
    "mRNA half-life: 9.5 min, corresponding to a degradation rate \u03b3_mRNA = 0.073 min\u207b\u00b9 "
    "(Rustad et al. 2013, global mRNA stability analysis in M. tuberculosis)"
)
add_bullet(
    "Protein half-life: 1,500 min, corresponding to a degradation rate \u03b3_protein = 0.000462 min\u207b\u00b9 "
    "(dilution-dominated in slowly growing M. tuberculosis)"
)
add_bullet(
    "Burst size: k_translation / \u03b3_mRNA = 6.85 proteins per mRNA lifetime"
)

add_body(
    "Each simulation was run for 50,000 cells per condition, with a total simulation time of "
    "30,000 minutes and a burn-in period of 15,000 minutes. Only data collected after burn-in "
    "were used for analysis. Random number seeds were fixed deterministically (master_seed = 42 "
    "plus a condition-specific offset) to ensure complete reproducibility."
)

# 3.7
add_subsection_header("3.7 Simulation Conditions")

add_body("Five primary simulation conditions were defined:")

add_bullet(
    "Condition A (Asymmetric, native): Kd_strong = 2.4 nM, Kd_weak = 49 nM, block_strong = 0.85, "
    "block_weak = 0.50. This represents the native Mce3R operator as characterized by Panagoda "
    "et al. (2024)."
)
add_bullet(
    "Condition B (Symmetric control): Kd = sqrt(2.4 \u00d7 49.0) = 10.84 nM for both sites, "
    "block = sqrt(0.85 \u00d7 0.50) = 0.652 for both sites. The geometric mean ensures equivalent "
    "total regulatory capacity, providing a fair noise comparison that isolates the effect of "
    "asymmetry from the effect of total repression strength."
)
add_bullet(
    "Condition C (Single-site): Kd_weak set to 1e12 nM (effectively infinite), ensuring the "
    "weak site is never occupied. This tests the contribution of the second binding site."
)
add_bullet(
    "Condition D (Unregulated): k_on = 0 (no Mce3R binding), representing constitutive "
    "expression without any repression."
)
add_bullet(
    "Condition E (Asymmetry sweep): 10 Kd ratios ranging from 1\u00d7 to 50\u00d7, holding the "
    "geometric mean Kd constant, to systematically explore the relationship between asymmetry "
    "magnitude and noise output."
)

# 3.8
add_subsection_header("3.8 Thermodynamic Calibration (Berg\u2013von Hippel Model)")

add_body(
    "Per-position binding energies were computed using the Berg\u2013von Hippel model: \u0394G_i = "
    "\u2013kT \u00d7 ln(f_i / p_i), where f_i is the frequency of the observed nucleotide at position i "
    "in the PWM and p_i is the background frequency of that nucleotide in the H37Rv genome. "
    "The temperature was set to 310 K (human body temperature), giving kT = 0.616 kcal/mol."
)

add_body(
    "The total binding energy for each site was computed as the sum of per-position energies and "
    "calibrated so that \u0394G_strong reproduces Kd = 2.4 nM and \u0394G_weak reproduces Kd = 49 nM. "
    "The four-state partition function was defined as: Z = 1 + exp(\u2013\u0394G_strong/kT)\u00b7[R] + "
    "exp(\u2013\u0394G_weak/kT)\u00b7[R] + \u03c9\u00b7exp(\u2013(\u0394G_strong + \u0394G_weak + \u0394G_spacer)/kT)\u00b7[R]\u00b2, "
    "where [R] is the free Mce3R concentration, \u03c9 is the cooperativity parameter, and "
    "\u0394G_spacer is the energetic penalty associated with the 53 bp spacer between binding sites."
)

# 3.9
add_subsection_header("3.9 MCMC Cooperativity Inference")

add_body(
    "The cooperativity parameter \u03c9 and spacer energy penalty \u0394G_spacer were inferred using "
    "Markov Chain Monte Carlo (MCMC) sampling with the emcee ensemble sampler. The sampler was "
    "configured with 32 walkers running 5,000 steps each, with a 1,000-step burn-in period. "
    "The two free parameters were ln(\u03c9) and \u0394G_spacer, with uniform priors. The likelihood "
    "function compared predicted mean protein level and CV against Phase 2 simulation outputs."
)

add_body(
    "Convergence was assessed using the Gelman\u2013Rubin statistic (R-hat < 1.01) and acceptance "
    "fraction (0.70). The inferred posterior yielded \u03c9 = 1.08 with a 95% credible interval of "
    "0.20\u20135.89, indicating minimal cooperativity (the two binding sites operate largely "
    "independently) but with substantial uncertainty reflecting the limited calibration data "
    "(two Kd measurements). Posterior propagation was performed by drawing 50 samples from the "
    "posterior and running each through independent Gillespie simulations with 2,000 cells per "
    "sample to test the robustness of the CV difference across the full uncertainty range."
)

# 3.10
add_subsection_header("3.10 Environmental Extensions")

add_body(
    "Four environmental conditions were modeled as scalar multipliers on Mce3R effective "
    "concentration, reflecting the known biology of M. tuberculosis within the host:"
)

add_bullet("Baseline: standard conditions with nominal Mce3R concentration.")
add_bullet(
    "Cholesterol: 2\u00d7 Mce3R effective concentration, reflecting cholesterol-induced upregulation "
    "of mce3R expression during intracellular infection."
)
add_bullet(
    "Acidic pH: 0.7\u00d7 Mce3R effective concentration and 1.3\u00d7 noise scale, reflecting "
    "phagosomal acidification stress that partially destabilizes repressor binding."
)
add_bullet(
    "Host-like: combined cholesterol and acidic pH effects, representing the complex intracellular "
    "environment of the macrophage phagosome."
)

add_body(
    "Additionally, a two-species model was implemented in which Mce3R autoregulates itself through "
    "negative feedback while simultaneously repressing the target gene. This model comprises 18 "
    "reactions: 8 operator state transitions, 2 transcription reactions (one for each gene), 2 "
    "translation reactions, 2 mRNA decay reactions, 2 protein decay reactions, and 2 placeholder "
    "reactions for potential future extensions. The dynamic free Mce3R concentration governs binding "
    "propensity in this model, allowing negative feedback to buffer noise. The total simulation "
    "matrix comprised 10,000 cells per condition across 3 architectures, 4 environments, and 2 "
    "model types, yielding 24 distinct simulation conditions."
)

# 3.11
add_subsection_header("3.11 Temporal Noise Dynamics")

add_body(
    "Full protein expression time traces were recorded from Gillespie simulations at 10-minute "
    "intervals after the burn-in period. The autocorrelation function was computed as C(\u03c4) = "
    "\u27e8\u03b4x(t) \u00b7 \u03b4x(t+\u03c4)\u27e9 / \u27e8\u03b4x\u00b2\u27e9 using an FFT-based method for computational efficiency, "
    "where \u03b4x(t) = x(t) \u2013 \u27e8x\u27e9 represents the deviation from the mean protein level. "
    "An exponential decay model C(\u03c4) ~ exp(\u2013\u03c4/\u03c4_c) was fitted to the autocorrelation function "
    "to extract the autocorrelation time \u03c4_c, which quantifies the temporal memory of noise "
    "fluctuations."
)

add_body(
    "Dwell time in the persister state was defined as contiguous intervals where the protein count "
    "falls below the persistence threshold (defined using three independent methods; see Section "
    "3.14). Dwell time distributions were computed for each architecture. The power spectral "
    "density S(f) was computed as the Fourier transform of the autocorrelation function to "
    "identify characteristic switching frequencies. All temporal dynamics were compared across "
    "the asymmetric, symmetric, and single-site architectures, and tested under both the "
    "single-species and two-species (autoregulatory) models."
)

# 3.12
add_subsection_header("3.12 Noise Quenching Dose-Response")

add_body(
    "The symmetrization sweep was conducted by varying Kd_weak from 49 nM toward 2.4 nM in 10 "
    "steps (49, 30, 20, 10.84 [symmetric], 5, 2.4 nM, plus intermediate values) while holding "
    "Kd_strong fixed at 2.4 nM. At each step, 5,000 cells were simulated, and CV, persister "
    "fraction, and autocorrelation time \u03c4_c were measured. A separate concentration sweep varied "
    "Mce3R from 10 to 10,000 nM in 15 log-spaced steps for both asymmetric and symmetric "
    "architectures."
)

add_body(
    "Two IC50 values were computed: IC50_noise, defined as the Kd_weak value at which CV drops "
    "halfway between the fully asymmetric (49 nM) and fully symmetric (10.84 nM) levels; and "
    "IC50_persistence, defined as the Kd_weak value at which persister fraction drops by 50%. "
    "These metrics generate a quantitative dose-response curve predicting how much operator "
    "modification would be needed to achieve clinically meaningful noise reduction."
)

# 3.13
add_subsection_header("3.13 Mutual Information")

add_body(
    "The environment was discretized into 4 states (baseline, cholesterol, acidic pH, host-like). "
    "For each architecture, 5,000 cells were simulated per environment and the resulting protein "
    "distributions were pooled. Mutual information was computed as MI = H(expression) \u2013 "
    "H(expression | environment) using a k-nearest-neighbor estimator, which provides a "
    "bias-corrected estimate suitable for continuous distributions. The MI was also computed "
    "across a sweep of 10 Mce3R concentrations to generate MI versus concentration curves. "
    "The comparison of MI(asymmetric) versus MI(symmetric) at the physiological Mce3R "
    "concentration (332 nM) tests whether the asymmetric architecture trades environmental "
    "information capacity for stochastic phenotypic diversification."
)

# 3.14
add_subsection_header("3.14 Persistence Threshold Calibration")

add_body(
    "To avoid dependence on an arbitrary persistence threshold definition, three independent "
    "threshold methods were employed:"
)

add_bullet(
    "Tail fraction method: The persistence threshold was set at the 1st percentile of each "
    "architecture\u2019s baseline protein distribution, defining persisters as cells in the lowest "
    "1% of expression."
)
add_bullet(
    "Absolute calibrated method: The threshold was back-calculated from the symmetric baseline "
    "distribution to yield approximately 0.1% persister fraction, matching published wild-type "
    "persister frequencies in M. tuberculosis."
)
add_bullet(
    "Fold-change method: Cells with protein levels below median/1.5 were classified as persisters, "
    "providing a relative threshold independent of absolute expression levels."
)

add_body(
    "Sensitivity analysis was performed across additional threshold percentiles (5th, 10th, 15th) "
    "to ensure that the qualitative conclusions are robust to the specific threshold chosen. "
    "All results are reported for all three methods."
)

# 3.15
add_subsection_header("3.15 Statistical Analysis")

add_body(
    "Bootstrap resampling with 1,000 iterations was used to construct confidence intervals for "
    "all CV estimates. The Kolmogorov\u2013Smirnov test was used for distribution comparisons between "
    "conditions. Cohen\u2019s d was computed as the effect size metric for CV differences between "
    "architectures. Gaussian mixture models with 1, 2, and 3 components were fitted to each "
    "condition\u2019s protein distribution and compared using the Bayesian Information Criterion (BIC) "
    "to test for bimodality. All p-values were two-sided, with statistical significance defined "
    "at \u03b1 = 0.05."
)

# 3.16
add_subsection_header("3.16 Software and Reproducibility")

add_body(
    "All analyses were performed using Python 3.10 with the following packages: NumPy, SciPy, "
    "Pandas, Matplotlib, Seaborn, BioPython, scikit-learn, statsmodels, and Numba for JIT "
    "compilation. Motif discovery and scanning used MEME Suite version 5.5.9. All random number "
    "seeds were fixed deterministically (master_seed = 42) to ensure complete reproducibility. "
    "The complete codebase is available upon request."
)

add_page_break()

# ══════════════════════════════════════════════════════════════════════════════
# 4. EXPECTED RESULTS AND SIGNIFICANCE
# ══════════════════════════════════════════════════════════════════════════════

add_section_header("4. EXPECTED RESULTS AND SIGNIFICANCE")

# 4.1
add_subsection_header("4.1 Aim 1 \u2014 Binding Site Validation")

add_body(
    "We expect MEME de novo motif discovery to identify conserved motifs that match the known "
    "Mce3R operator consensus sequence. The top-ranked FIMO hit should correspond to the known "
    "operator region at H37Rv coordinates 2,207,477\u20132,207,699 with high statistical significance "
    "(p < 1e\u221210). The genome-wide scan is expected to identify over 100 candidate binding sites, "
    "with enrichment near genes involved in lipid metabolism and cholesterol transport \u2014 "
    "functional categories consistent with the known role of the mce3 operon."
)

add_body(
    "Conservation of binding site sequences across the three mycobacterial species supports the "
    "functional importance of the Mce3R binding motif, though we note that the close evolutionary "
    "relationship between M. tuberculosis and M. bovis limits the strength of conservation "
    "arguments specifically for the asymmetry feature. Inclusion of the more distant M. marinum "
    "provides additional discriminatory power. Novel candidate binding sites near persistence-"
    "relevant genes such as tgs1 (triacylglycerol synthase), espC (ESX-1 secretion), and mbtG "
    "(mycobactin biosynthesis) would suggest an expanded Mce3R regulon with direct implications "
    "for persistence biology."
)

# 4.2
add_subsection_header("4.2 Aim 2 \u2014 Thermodynamic Model")

add_body(
    "We expect the Berg\u2013von Hippel calibration to faithfully reproduce the experimentally "
    "measured Kd values, validating the thermodynamic framework. The MCMC inference is expected "
    "to yield a cooperativity parameter \u03c9 near 1 (no strong cooperativity), with a wide posterior "
    "distribution reflecting the limited calibration data available (two Kd measurements). This "
    "result, while uncertain in its precise value, is informative: minimal cooperativity means "
    "the two binding sites operate largely independently, creating four genuinely distinct operator "
    "states rather than the effectively two-state system that strong cooperativity would produce."
)

add_body(
    "We expect all classified operators to be graded repressors (Hill coefficient n < 2), "
    "consistent with continuous metabolic tuning rather than switch-like digital behavior. The "
    "asymmetric operator should exhibit a broader repression transition region than the symmetric "
    "control, allowing more intermediate expression states and thus greater cell-to-cell "
    "variability. The significance of this finding is that independent binding creates more noise "
    "than cooperative binding: when the two sites can be occupied independently, the system "
    "samples more operator states, generating more expression heterogeneity."
)

# 4.3
add_subsection_header("4.3 Aim 3a \u2014 Noise Amplitude")

add_body(
    "We expect the asymmetric operator to produce significantly higher noise than the symmetric "
    "control, with CV(asymmetric) = 0.187 versus CV(symmetric) = 0.157 (p < 0.001), representing "
    "a 19% noise increase. We expect this difference to be robust across all MCMC posterior "
    "samples, the two-species autoregulatory model, all four environmental conditions, and the "
    "TetR/tetO2 biological benchmark. Corrected persister fractions should show 2\u20135-fold "
    "enrichment in the asymmetric architecture across all three threshold methods."
)

add_body(
    "It is important to contextualize the magnitude of this effect. A CV of 0.187 is moderate "
    "\u2014 higher than typical tightly regulated genes (CV ~ 0.15) but substantially lower than "
    "noise in toxin-antitoxin modules (CV ~ 0.45) or known high-noise systems. This suggests "
    "that operator asymmetry is one contributor to persistence-relevant noise among several, "
    "likely providing a constitutive noise floor upon which other mechanisms (TA module switching, "
    "metabolic fluctuations, stochastic partitioning) layer additional variability. The \u0394CV of "
    "0.031 is modest but statistically robust and, as we show in Aim 3c, has nonlinear "
    "consequences for the tail of the expression distribution where persisters reside."
)

# 4.4
add_subsection_header("4.4 Aim 3b \u2014 Temporal Dynamics (Key Novel Contribution)")

add_body(
    "We predict that the autocorrelation time of the asymmetric operator will exceed that of the "
    "symmetric control: \u03c4_c(asymmetric) > \u03c4_c(symmetric). The mechanistic basis for this "
    "prediction is that unequal binding sites create asymmetric switching between operator states. "
    "The strong site (Kd = 2.4 nM) has a slow off-rate and remains bound most of the time, while "
    "the weak site (Kd = 49 nM) flickers on and off more rapidly. This creates extended periods "
    "in partially-repressed states, prolonging low-expression episodes that correspond to the "
    "persister phenotype."
)

add_body(
    "We further predict that the dwell time in the persister state will be longer for the "
    "asymmetric architecture than for the symmetric control. This is perhaps the most clinically "
    "significant prediction of this study: if validated experimentally, it would mean that the "
    "asymmetric operator does not merely create more persisters \u2014 it creates persisters that "
    "remain dormant for longer periods. This directly connects operator architecture to treatment "
    "duration, providing a molecular explanation for why Mce3R deletion (which removes the "
    "asymmetric operator\u2019s influence) affects persistence (Pandey et al. 2023)."
)

add_body(
    "We acknowledge that if the autocorrelation time is instead dominated by protein half-life "
    "(which sets a fundamental lower bound on noise memory through the dilution timescale), "
    "architecture-dependent temporal effects could be washed out. This outcome would still be "
    "informative: it would indicate that noise amplitude, not noise memory, is the relevant "
    "quantity connecting operator architecture to persistence, narrowing the space of mechanistic "
    "models."
)

# 4.5
add_subsection_header("4.5 Aim 3c \u2014 Noise Quenching (Key Novel Contribution)")

add_body(
    "We expect the noise quenching dose-response to follow a sigmoidal curve: CV decreases as "
    "Kd_weak is swept from 49 nM toward 2.4 nM (full symmetrization). We predict a specific "
    "IC50_noise value \u2014 the Kd_weak at which noise drops halfway between the fully asymmetric "
    "and fully symmetric levels \u2014 which constitutes a quantitative therapeutic target. We further "
    "predict that the persister fraction will drop more steeply than CV, reflecting the nonlinear "
    "relationship between noise amplitude and tail probability in the expression distribution."
)

add_body(
    "The significance of this analysis is that it provides the first quantitative prediction of "
    "how much operator modification would be needed to reduce persistence in a specific TB "
    "regulatory system. Flentie et al. (2019) demonstrated that small-molecule compounds "
    "targeting Mce3R (6-azasteroids) enhance antibiotic efficacy 16\u201350-fold in vitro. Our model "
    "suggests a mechanism for this enhancement \u2014 Mce3R-targeting compounds may effectively "
    "symmetrize the operator by altering repressor\u2013DNA interactions \u2014 and predicts the required "
    "dose of symmetrization to achieve specific reductions in persister frequency. This bridges "
    "the gap between drug discovery and mechanistic understanding of persistence."
)

# 4.6
add_subsection_header("4.6 Mutual Information Trade-off")

add_body(
    "We expect the asymmetric operator to transmit less mutual information about the environment "
    "than the symmetric control: MI(asymmetric) < MI(symmetric) at physiological Mce3R "
    "concentration. The interpretation of this finding would be that the asymmetric architecture "
    "sacrifices precise environmental tracking \u2014 the ability to faithfully transduce changes in "
    "Mce3R concentration into proportional changes in gene expression \u2014 for stochastic "
    "diversification of the population. This is consistent with a bet-hedging strategy (Veening "
    "et al. 2008), in which the bacterium maintains a subpopulation of metabolically distinct "
    "cells as insurance against unpredictable environmental challenges, including antibiotic "
    "exposure."
)

add_body(
    "The mutual information analysis places the Mce3R system within the broader theoretical "
    "framework of information-noise trade-offs in gene regulation (Rosenfeld et al. 2005). If "
    "confirmed, this result would suggest that the asymmetric operator architecture has been "
    "shaped by natural selection to balance two competing demands: the need to respond to "
    "environmental cues (requiring low noise and high MI) and the need to maintain phenotypic "
    "diversity (requiring high noise and accepting lower MI)."
)

add_page_break()

# ══════════════════════════════════════════════════════════════════════════════
# 5. LIMITATIONS AND FALSIFICATION CRITERIA
# ══════════════════════════════════════════════════════════════════════════════

add_section_header("5. LIMITATIONS AND FALSIFICATION CRITERIA")

add_subsection_header("5.1 Limitations")

add_body(
    "We identify the following limitations of this study, which should be considered when "
    "interpreting our results:"
)

add_body(
    "First, all results are computational predictions. No experimental validation currently exists "
    "for the noise-to-persistence connection in the Mce3R system specifically. While our model is "
    "grounded in experimentally measured biophysical parameters (Panagoda et al. 2024), the "
    "predictions require direct experimental testing, which we discuss as future directions."
)

add_body(
    "Second, the model uses a single-gene or two-gene representation of what is in reality a "
    "14-gene operon. Since all genes in the mce3 operon share the same promoter and operator, "
    "noise patterns at the transcriptional level should be similar across all genes in the operon. "
    "However, downstream post-transcriptional regulatory effects, including differential mRNA "
    "stability and translational efficiency, are not captured in our model."
)

add_body(
    "Third, the MCMC posterior for the cooperativity parameter \u03c9 is wide (95% CI: 0.20\u20135.89), "
    "reflecting the fact that only two calibration data points (the two Kd values) are available. "
    "The CV difference between asymmetric and symmetric architectures is robust across the full "
    "posterior, but the precise cooperativity value is uncertain. Additional experimental "
    "measurements (e.g., cooperative binding assays) would substantially narrow this posterior."
)

add_body(
    "Fourth, environmental conditions are modeled as static scalar multipliers on Mce3R "
    "concentration, not as dynamic time-varying signals. In reality, the intracellular environment "
    "of the macrophage phagosome fluctuates over time, and these fluctuations could interact with "
    "expression noise in ways our model does not capture."
)

add_body(
    "Fifth, the symmetric control is a theoretical construct (the geometric mean of the two "
    "native Kd values). No naturally occurring symmetric Mce3R operator exists for comparison. "
    "The TetR/tetO2 benchmark provides biological grounding for the modeling framework but is "
    "from a different species (E. coli) and a different regulatory system."
)

add_body(
    "Sixth, conservation analysis across M. tuberculosis, M. bovis, and M. marinum supports "
    "functional importance of the binding motif but does not specifically demonstrate selection "
    "for the asymmetry feature. The close evolutionary relationship between M. tuberculosis and "
    "M. bovis limits the power of this comparison."
)

add_body(
    "Seventh, the CV difference (\u0394CV = 0.031) is moderate in the context of total cellular "
    "noise. Other noise-generating mechanisms \u2014 including toxin-antitoxin module switching, "
    "metabolic fluctuations, stochastic partitioning at cell division, and upstream transcription "
    "factor noise \u2014 likely contribute more to total expression variability. The significance of "
    "operator asymmetry may be as a constitutive noise floor rather than the dominant noise source."
)

add_body(
    "Eighth, temporal dynamics predictions (autocorrelation time, dwell time) could be dominated "
    "by protein half-life rather than operator architecture, particularly for the long-lived "
    "proteins typical of M. tuberculosis."
)

add_subsection_header("5.2 Falsification Criteria")

add_body(
    "We specify the following experimental outcomes that would disprove our central hypotheses, "
    "demonstrating that the study generates genuinely falsifiable predictions:"
)

add_bullet(
    "If experimentally symmetrizing the Mce3R operator (by mutating the weak site to match the "
    "strong site) showed no change in expression noise or persister frequency, the core hypothesis "
    "that asymmetry generates noise would be falsified."
)
add_bullet(
    "If adding realistic upstream noise sources (transcription factor copy number fluctuations, "
    "global growth rate variation) to the model eliminated the CV difference between asymmetric "
    "and symmetric architectures, this would indicate that operator-level noise is negligible "
    "compared to extrinsic noise sources."
)
add_bullet(
    "If strong cooperativity (\u03c9 >> 1) were found with better calibration data, this would "
    "suppress independent binding and convert the four-state system into an effectively two-state "
    "(all-or-nothing) system, potentially eliminating the noise advantage of asymmetry."
)
add_bullet(
    "If the two-species autoregulated model showed CV(asymmetric) \u2264 CV(symmetric), this would "
    "indicate that negative feedback completely abolishes the noise-generating effect of operator "
    "asymmetry under biologically realistic conditions."
)
add_bullet(
    "If the autocorrelation time showed no architecture dependence (being dominated entirely by "
    "protein half-life), the temporal dynamics claim would be falsified, though the noise "
    "amplitude claim would still stand as an independent prediction."
)

add_page_break()

# ══════════════════════════════════════════════════════════════════════════════
# 6. INNOVATION
# ══════════════════════════════════════════════════════════════════════════════

add_section_header("6. INNOVATION")

add_body(
    "This project makes three novel contributions to the field of stochastic gene regulation "
    "and tuberculosis persistence biology:"
)

add_body(
    "First, we provide the first comprehensive noise characterization \u2014 encompassing amplitude, "
    "temporal dynamics, and quenching response \u2014 of a structurally resolved asymmetric bacterial "
    "operator with experimentally measured Kd values. Previous studies have either characterized "
    "noise in systems without structural information or characterized structure without measuring "
    "noise. By integrating the cryo-EM structural data of Panagoda et al. (2024) with stochastic "
    "simulation, we bridge these two domains."
)

add_body(
    "Second, we provide the first prediction of noise temporal memory arising from operator "
    "architecture. While noise amplitude (CV) has been extensively studied, the temporal "
    "persistence of noise fluctuations \u2014 quantified by the autocorrelation time \u03c4_c and dwell "
    "time in the persister state \u2014 has not been connected to cis-regulatory DNA structure. Our "
    "prediction that asymmetry increases not just the frequency but also the duration of "
    "persistence episodes, if validated, would fundamentally change how we think about the "
    "relationship between regulatory architecture and treatment duration."
)

add_body(
    "Third, we provide the first quantitative dose-response curve for noise quenching in a "
    "specific TB regulatory system. By computing the IC50 of noise and the IC50 of persistence, "
    "we generate specific, testable predictions for how much operator modification would be "
    "required to achieve clinically meaningful reductions in persister frequency. This establishes "
    "a framework for \u201canti-noise\u201d therapeutic strategies that target the source of phenotypic "
    "heterogeneity rather than the downstream phenotype."
)

add_body(
    "These contributions build on well-established principles \u2014 that cis-regulatory architecture "
    "shapes noise and that noise drives persistence \u2014 but apply them to a specific, previously "
    "uncharacterized system with real biophysical parameters, generating testable predictions "
    "that bridge structural biology (Panagoda et al. 2024) and persistence microbiology (Pandey "
    "et al. 2023)."
)

add_page_break()

# ══════════════════════════════════════════════════════════════════════════════
# 7. TIMELINE
# ══════════════════════════════════════════════════════════════════════════════

add_section_header("7. TIMELINE")

# Create a table for the timeline
table = doc.add_table(rows=8, cols=2)
table.style = 'Table Grid'

# Header row
hdr = table.rows[0]
for cell, text in zip(hdr.cells, ["Month", "Activities"]):
    cell.text = ""
    run = cell.paragraphs[0].add_run(text)
    set_font(run, bold=True, size=11)

timeline_data = [
    ("Month 1", "Aim 1: De novo motif discovery (MEME), genome-wide scanning (FIMO), cross-species conservation analysis"),
    ("Month 2", "Aim 2a\u20132b: Gillespie stochastic simulation implementation, four-state operator model, thermodynamic calibration (Berg\u2013von Hippel)"),
    ("Month 3", "Aim 2c\u20132d: MCMC cooperativity inference, operator classification, partition function fitting"),
    ("Month 4", "Aim 3a: Noise amplitude quantification, environmental extensions, two-species autoregulatory model, TetR/tetO2 benchmark validation"),
    ("Month 5", "Aim 3b: Temporal noise dynamics \u2014 autocorrelation function, dwell time distributions, power spectral density analysis"),
    ("Month 6", "Aim 3c\u20133d: Noise quenching dose-response curves, IC50 computation, mutual information analysis"),
    ("Month 7", "Comprehensive analysis, figure preparation, manuscript writing, and submission preparation"),
]

for i, (month, activities) in enumerate(timeline_data):
    row = table.rows[i + 1]
    row.cells[0].text = ""
    run = row.cells[0].paragraphs[0].add_run(month)
    set_font(run, bold=True, size=11)
    row.cells[1].text = ""
    run = row.cells[1].paragraphs[0].add_run(activities)
    set_font(run, size=11)

add_page_break()

# ══════════════════════════════════════════════════════════════════════════════
# 8. REFERENCES
# ══════════════════════════════════════════════════════════════════════════════

add_section_header("8. REFERENCES")

references = [
    "WHO (2024). Global Tuberculosis Report. World Health Organization, Geneva.",
    "Balaban, N. Q., et al. (2019). Definitions and guidelines for research on antibiotic persistence. Nature Reviews Microbiology, 17(7), 441\u2013448.",
    "Elowitz, M. B., et al. (2002). Stochastic gene expression in a single cell. Science, 297(5584), 1183\u20131186.",
    "Chowdhury, D., et al. (2021). Cis-regulatory logic produces gene-expression noise describing phenotypic heterogeneity in bacteria. Frontiers in Genetics, 12, 698910.",
    "Lengyel, I. M., & Morelli, L. G. (2017). Multiple binding sites for transcriptional repressors can produce regular bursting and enhance noise suppression. Physical Review E, 95, 042412.",
    "Balazsi, G., van Oudenaarden, A., & Collins, J. J. (2011). Cellular decision making and biological noise: from microbes to mammals. Cell, 144(6), 910\u2013925.",
    "Quigley, J., & Lewis, K. (2022). Noise in a metabolic pathway leads to persister formation in Mycobacterium tuberculosis. Microbiology Spectrum, 10(4), e0094822.",
    "Santangelo, M. P., et al. (2002). Characterization of Mce3R, a TetR-type transcriptional repressor linking the mce3 virulon to the metabolic network in Mycobacterium tuberculosis. Microbiology, 148, 2997\u20133006.",
    "Dunphy, K. Y., et al. (2010). Attenuation of Mycobacterium tuberculosis functionally disrupted in a fatty acyl-CoA synthetase gene fadD5. Journal of Infectious Diseases, 201(8), 1232\u20131239.",
    "Pandey, S., et al. (2023). Deletion of mce3R increases persister frequency in Mycobacterium tuberculosis. Research in Microbiology, 174, 104082.",
    "Panagoda, G. J., Balazsi, G., & Sampson, N. S. (2024). Structural and biophysical characterization of Mce3R and its asymmetric operator. ACS Chemical Biology, 19, 2580\u20132592.",
    "Cock, P. J. A., et al. (2009). Biopython: freely available Python tools for computational molecular biology and bioinformatics. Bioinformatics, 25(11), 1422\u20131423.",
    "Bailey, T. L., et al. (2009). MEME Suite: tools for motif discovery and searching. Nucleic Acids Research, 37(Web Server issue), W202\u2013W208.",
    "Grant, C. E., et al. (2011). FIMO: scanning for occurrences of a given motif. Bioinformatics, 27(7), 1017\u20131018.",
    "Gillespie, D. T. (1977). Exact stochastic simulation of coupled chemical reactions. Journal of Physical Chemistry, 81(25), 2340\u20132361.",
    "Stormo, G. D., & Zhao, Y. (2010). Determining the specificity of protein\u2013DNA interactions. Nature Reviews Genetics, 11, 751\u2013760.",
    "Taniguchi, Y., et al. (2010). Quantifying E. coli proteome and transcriptome with single-molecule sensitivity in single cells. Science, 329(5991), 533\u2013538.",
    "Rustad, T. R., et al. (2013). Global analysis of mRNA stability in Mycobacterium tuberculosis. Nucleic Acids Research, 41(1), 509\u2013517.",
    "Santangelo, M. P., et al. (2009). Study of the role of Mce3R on the transcription of mce genes of Mycobacterium tuberculosis. Microbiology, 155(3), 882\u2013891.",
    "Flentie, K., et al. (2019). Chemical disarming of isoniazid resistance \u2014 Mycobacterium tuberculosis N-acetyltransferase with 6-azasteroid targeting of the Mce3R pathway. ACS Infectious Diseases, 5(7), 1239\u20131254.",
    "Veening, J. W., Smits, W. K., & Kuipers, O. P. (2008). Bistability, epigenetics, and bet-hedging in bacteria. Annual Review of Microbiology, 62, 193\u2013210.",
    "Rosenfeld, N., et al. (2005). Gene regulation at the single-cell level. Science, 307(5717), 1962\u20131965.",
    "Farquhar, K. S., et al. (2019). Role of network-mediated stochasticity in mammalian drug resistance. Nature Communications, 10, 2766.",
    "Rotem, E., et al. (2010). Regulation of phenotypic variability by a threshold-based mechanism underlies bacterial persistence. Proceedings of the National Academy of Sciences, 107(28), 12541\u201312546.",
    "El Meouche, I., Siu, Y., & Bhatt, M. J. (2016). Stochastic expression of a multiple antibiotic resistance activator confers transient resistance in single cells. Molecular Cell, 62(1), 84\u201394.",
    "Wakamoto, Y., et al. (2013). Dynamic persistence of antibiotic-stressed mycobacteria. Science, 339(6115), 91\u201395.",
    "Sureka, K., et al. (2008). Positive feedback and noise activate the stringent response regulator rel in mycobacteria. PLoS ONE, 3(3), e1771.",
]

for i, ref in enumerate(references):
    p = doc.add_paragraph()
    p.paragraph_format.line_spacing = 2.0
    p.paragraph_format.left_indent = Inches(0.5)
    p.paragraph_format.first_line_indent = Inches(-0.5)
    run = p.add_run(f"[{i+1}]  {ref}")
    set_font(run, size=10)

# ── Save ─────────────────────────────────────────────────────────────────────
doc.save(OUTPUT_PATH)
print(f"Document saved to: {OUTPUT_PATH}")

# Report size
size_bytes = os.path.getsize(OUTPUT_PATH)
if size_bytes > 1024 * 1024:
    print(f"File size: {size_bytes / (1024*1024):.2f} MB")
else:
    print(f"File size: {size_bytes / 1024:.1f} KB")
