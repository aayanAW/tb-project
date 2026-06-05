#!/usr/bin/env python3
"""
Generate ENIGMA Project Proposal v4 — The Persistence Regulatory Code.
Definitive project proposal incorporating all revisions.
"""

from docx import Document
from docx.shared import Pt, Inches, Cm, RGBColor
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.enum.section import WD_ORIENT
from docx.oxml.ns import qn, nsdecls
from docx.oxml import parse_xml
import os

OUTPUT_PATH = "/Users/aayanalwani/tb project/mce3r_stochastic/ENIGMA_Project_Proposal_v4.docx"

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

for section in doc.sections:
    section.top_margin = Inches(1)
    section.bottom_margin = Inches(1)
    section.left_margin = Inches(1)
    section.right_margin = Inches(1)


# ── Helper functions ─────────────────────────────────────────────────────────
def add_title(text, size=16):
    p = doc.add_paragraph()
    p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    p.paragraph_format.line_spacing = 2.0
    run = p.add_run(text)
    run.bold = True
    run.font.size = Pt(size)
    run.font.name = 'Calibri'
    return p


def add_subtitle(text, size=14):
    p = doc.add_paragraph()
    p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    p.paragraph_format.line_spacing = 2.0
    run = p.add_run(text)
    run.italic = True
    run.font.size = Pt(size)
    run.font.name = 'Calibri'
    return p


def add_section_header(text, size=14):
    p = doc.add_paragraph()
    p.paragraph_format.line_spacing = 2.0
    p.paragraph_format.space_before = Pt(12)
    run = p.add_run(text)
    run.bold = True
    run.font.size = Pt(size)
    run.font.name = 'Calibri'
    return p


def add_subsection_header(text, size=12):
    p = doc.add_paragraph()
    p.paragraph_format.line_spacing = 2.0
    p.paragraph_format.space_before = Pt(8)
    run = p.add_run(text)
    run.bold = True
    run.font.size = Pt(size)
    run.font.name = 'Calibri'
    return p


def add_body(text):
    p = doc.add_paragraph()
    p.paragraph_format.line_spacing = 2.0
    run = p.add_run(text)
    run.font.size = Pt(11)
    run.font.name = 'Calibri'
    return p


def add_body_runs(runs_spec):
    """Add a paragraph with mixed formatting. runs_spec is list of (text, bold, italic)."""
    p = doc.add_paragraph()
    p.paragraph_format.line_spacing = 2.0
    for text, bold, italic in runs_spec:
        run = p.add_run(text)
        run.font.size = Pt(11)
        run.font.name = 'Calibri'
        run.bold = bold
        run.italic = italic
    return p


def add_bullet(text, level=0):
    p = doc.add_paragraph(style='List Bullet')
    p.paragraph_format.line_spacing = 2.0
    if level > 0:
        p.paragraph_format.left_indent = Inches(0.5 * (level + 1))
    # Clear default and add formatted run
    p.clear()
    run = p.add_run(text)
    run.font.size = Pt(11)
    run.font.name = 'Calibri'
    return p


def page_break():
    doc.add_page_break()


# =============================================================================
# TITLE PAGE
# =============================================================================
# Add vertical spacing
for _ in range(6):
    p = doc.add_paragraph()
    p.paragraph_format.line_spacing = 2.0

add_title("The Persistence Regulatory Code:", 16)
add_title("Genome-Wide Operator Asymmetry, Temporal Noise Memory,", 16)
add_title("and Therapeutic Noise Quenching in", 16)
add_title("Mycobacterium tuberculosis", 16)

p = doc.add_paragraph()
p.paragraph_format.line_spacing = 2.0

add_subtitle("ENIGMA: Expression Noise In Gene-regulatory Mechanisms and Architecture", 14)

for _ in range(4):
    p = doc.add_paragraph()
    p.paragraph_format.line_spacing = 2.0

p = doc.add_paragraph()
p.alignment = WD_ALIGN_PARAGRAPH.CENTER
p.paragraph_format.line_spacing = 2.0
run = p.add_run("Project Proposal v4 — Final Reframed Version")
run.font.size = Pt(12)
run.font.name = 'Calibri'
run.italic = True

page_break()

# =============================================================================
# ABSTRACT
# =============================================================================
add_section_header("ABSTRACT")

add_body(
    "Mycobacterium tuberculosis kills over 1.2 million people annually. Treatment requires "
    "6-9 months of daily antibiotics because phenotypically tolerant persister cells, generated "
    "by stochastic gene expression noise, survive drug exposure and drive relapse. While it is "
    "established that cis-regulatory architecture shapes noise (Chowdhury et al. 2021) and that "
    "noise drives M. tuberculosis persistence (Quigley & Lewis 2022), no study has systematically "
    "characterized how operator-level DNA architecture contributes to persistence-relevant noise "
    "across the M. tuberculosis genome."
)

add_body(
    "We address this through a four-aim computational framework anchored by the Mce3R system, "
    "a TetR-family repressor with an unprecedented 20.4-fold binding site asymmetry (Kd = 2.4 nM "
    "vs. 49 nM; Panagoda et al. 2024) whose deletion increases persister frequency (Pandey et al. "
    "2023)."
)

add_body(
    "Aim 1: We discover and validate Mce3R binding sites using de novo motif discovery (MEME/FIMO) "
    "across three mycobacterial genomes, then extend genome-wide by scoring all 1,442 predicted "
    "operator sites for binding asymmetry and cross-referencing against published Tn-seq persistence "
    "gene screens. We test whether persistence-critical genes are enriched for high-asymmetry operators."
)

add_body(
    "Aim 2: We construct a thermodynamically calibrated four-state stochastic model using Gillespie "
    "simulation, Berg-von Hippel energy calibration, and MCMC cooperativity inference (omega = 1.08, "
    "95% CI: 0.20-5.89), establishing that the Mce3R operator functions as a graded repressor "
    "without strong cooperativity."
)

add_body(
    "Aim 3: We characterize noise across three dimensions: (a) amplitude, with CV = 0.187 asymmetric "
    "vs. 0.157 symmetric, robust across autoregulation, four environments, and the full MCMC posterior; "
    "(b) temporal dynamics, including autocorrelation time, dwell time in the persister state, and "
    "multi-generational lineage memory through cell division; (c) therapeutic quenching, with "
    "dose-response curves predicting how much operator symmetrization is needed to reduce the persister "
    "subpopulation, generating quantitative targets for anti-noise interventions."
)

add_body(
    "Aim 4 (Stretch): We wrap the Gillespie engine in an evolutionary simulation using a fast "
    "surrogate fitness model, testing whether asymmetric operator architecture is preferentially "
    "selected under antibiotic pulse pressure in persistence-linked genes versus housekeeping genes."
)

add_body(
    "Results indicate that operator asymmetry generates moderate but robust expression noise "
    "(delta-CV = 0.031, p < 0.001), that persistence-associated regulons may be enriched for "
    "high-asymmetry architectures, and that the asymmetric operator trades environmental information "
    "capacity (MI = 0.28 bits) for phenotypic diversification. All results are computational "
    "predictions requiring experimental validation, with explicit falsification criteria provided."
)

page_break()

# =============================================================================
# 1. INTRODUCTION AND SIGNIFICANCE
# =============================================================================
add_section_header("1. INTRODUCTION AND SIGNIFICANCE")

# 1.1
add_subsection_header("1.1 The Problem: Tuberculosis Persistence")

add_body(
    "Tuberculosis remains one of humanity's most devastating infectious diseases. According to the "
    "World Health Organization's 2024 Global Tuberculosis Report, Mycobacterium tuberculosis kills "
    "over 1.2 million people every year, making it the leading bacterial cause of death worldwide. "
    "Despite the availability of effective antibiotics, tuberculosis treatment requires an extraordinary "
    "6-9 months of daily multi-drug therapy. This prolonged regimen is not because the drugs fail to "
    "kill M. tuberculosis cells; rather, it exists because a small subpopulation of phenotypically "
    "tolerant persister cells survives drug exposure and drives disease relapse."
)

add_body(
    "Persister cells are fundamentally different from genetically resistant mutants. They carry no "
    "resistance mutations and are genetically identical to drug-sensitive cells in the same population. "
    "Instead, persisters are cells that happened to be in a metabolically dormant or slow-growing state "
    "at the moment antibiotics arrived (Balaban et al. 2019). Because most antibiotics target active "
    "cellular processes such as cell wall synthesis, DNA replication, or protein translation, cells that "
    "have transiently shut down these processes are effectively invisible to the drugs. Once antibiotic "
    "pressure is removed, persisters can resume growth and re-establish the infection."
)

add_body(
    "The clinical consequences are severe. Treatment duration is determined not by how quickly drugs "
    "kill actively growing cells, but by how many persisters form and how long they remain dormant. "
    "The extended treatment creates enormous burdens: poor patient compliance, increased side effects, "
    "greater healthcare costs, and paradoxically, the extended exposure to subtherapeutic drug levels "
    "during periods of non-compliance promotes the evolution of true genetic resistance. If we could "
    "reduce the persister subpopulation or shorten the duration of persister dormancy, we could "
    "dramatically shorten treatment and improve outcomes for millions of patients."
)

add_body(
    "Current drugs cannot distinguish active from dormant cells. What is needed is a fundamentally "
    "new approach: understanding the molecular mechanisms that generate persisters in the first place, "
    "and developing strategies to suppress persister formation at its source. This is the central "
    "motivation of our project."
)

# 1.2
add_subsection_header("1.2 Noise and Persistence: Established Principles")

add_body(
    "The following principles are well established in the literature and represent the scientific "
    "foundation upon which our project builds. We state them explicitly to distinguish established "
    "knowledge from our novel contributions."
)

add_body(
    "Gene expression is inherently stochastic. Even in genetically identical cells grown under "
    "identical conditions, protein and mRNA levels vary from cell to cell due to the random nature "
    "of molecular events such as transcription factor binding, RNA polymerase initiation, and "
    "ribosome loading. This was definitively demonstrated by Elowitz et al. (2002) using dual-reporter "
    "systems in E. coli, and has since been confirmed across all domains of life. The magnitude of this "
    "stochastic variation, termed gene expression noise, is not random with respect to gene regulatory "
    "architecture: it is shaped by the cis-regulatory elements that control gene expression."
)

add_body(
    "Cis-regulatory architecture, encompassing promoter structure, operator design, and binding site "
    "arrangement, shapes noise amplitude and dynamics. Chowdhury et al. (2021) demonstrated that the "
    "number, affinity, and arrangement of transcription factor binding sites within a regulatory region "
    "directly influence the coefficient of variation of downstream gene expression. Lengyel and Morelli "
    "(2017) showed theoretically that multiple binding sites modulate noise in configuration-dependent "
    "ways, with the spatial arrangement of sites mattering as much as their number. Balazsi et al. "
    "(2011) established that regulatory circuit architecture, including feedback loops, autoregulation, "
    "and multi-layer cascades, controls noise-driven cellular decisions, determining whether noise is "
    "amplified, filtered, or converted into discrete phenotypic states."
)

add_body(
    "In M. tuberculosis specifically, noise in metabolic genes creates persisters. Quigley and Lewis "
    "(2022) demonstrated this directly by showing that stochastic fluctuations in the expression of "
    "the metabolic gene ackA generate a subpopulation of cells with reduced metabolic activity that "
    "are tolerant to antibiotics. Critically, they showed that overexpressing ackA to reduce noise "
    "decreased persister frequency, providing direct causal evidence that noise drives persistence in "
    "M. tuberculosis. Wakamoto et al. (2013) provided complementary evidence by showing that "
    "stochastic variation in catalase-peroxidase (KatG) expression creates subpopulations with "
    "differential isoniazid susceptibility."
)

add_body(
    "The general principle is therefore clear: regulatory architecture shapes noise, and noise drives "
    "persistence. What remains unknown is which specific architectural features matter most, how much "
    "each feature contributes quantitatively to noise, and whether these features are systematically "
    "enriched in persistence-critical genes across the M. tuberculosis genome. Our project addresses "
    "these open questions."
)

# 1.3
add_subsection_header("1.3 The Mce3R System: A Case Study")

add_body(
    "The Mce3R system provides an ideal model for investigating how operator-level DNA architecture "
    "contributes to persistence-relevant noise. Mce3R (Rv1963c) is a TetR-family transcriptional "
    "repressor that controls the mce3 operon (Rv1964-Rv1977), which encodes cholesterol and lipid "
    "transport machinery essential for survival within macrophages (Santangelo et al. 2002; Dunphy "
    "et al. 2010). The mce3 system is one of four mammalian cell entry (mce) operons in M. tuberculosis, "
    "and its products are required for the bacterium to access host-derived lipids as carbon sources "
    "during intracellular infection."
)

add_body(
    "Two recent findings make the Mce3R system directly relevant to persistence. First, Pandey et al. "
    "(2023) demonstrated that deletion of the mce3R gene increases persister frequency, establishing "
    "a direct genetic link between this regulatory system and the persister phenotype. Second, Panagoda "
    "et al. (2024) solved the cryo-EM structure of Mce3R bound to its operator DNA (PDB 9B7Y), "
    "revealing a remarkable architectural feature: the Mce3R operator is a 123 base pair bipartite "
    "element containing two binding sites with dramatically different affinities."
)

add_body(
    "The strong binding site has a dissociation constant (Kd) of 2.4 plus or minus 0.7 nM, while the "
    "weak binding site has a Kd of approximately 49 nM. This represents a 20.4-fold difference in "
    "binding affinity, the largest documented asymmetry in any TetR-family system. The two sites are "
    "separated by a 53 base pair spacer. Mce3R itself is a double TetR-fold repeat protein, meaning "
    "each monomer contains two tandem copies of the canonical TetR helix-turn-helix domain, creating "
    "an unusually large DNA-binding surface."
)

add_body(
    "This system is uniquely suited for computational study because it provides what most theoretical "
    "noise models lack: real, structurally characterized dissociation constants for both binding sites, "
    "derived from biophysical measurements on the actual protein-DNA complex. Rather than exploring "
    "noise as a function of hypothetical parameters, we can model noise using the true binding "
    "energetics of a biologically important regulatory system."
)

# 1.4
add_subsection_header("1.4 The Gap: Three Open Questions")

add_body(
    "Our project addresses three open questions that emerge from the intersection of the established "
    "principles above and the specific biology of the Mce3R system."
)

add_body(
    "Question 1: Does the specific 20.4-fold Kd asymmetry produce quantifiably more noise than "
    "symmetric alternatives? This is a single-system question that has been partially explored in the "
    "theoretical literature for generic operators, but has never been addressed for this specific "
    "operator with its real, cryo-EM-derived binding parameters. The answer will establish whether "
    "the Mce3R operator's architecture is a noise-generating feature or a noise-neutral structural "
    "consequence of other functional constraints."
)

add_body(
    "Question 2: Is operator asymmetry enriched in persistence-critical genes across the M. tuberculosis "
    "genome? This genome-wide question is completely unexplored. If high-asymmetry operators are "
    "systematically deployed at genes involved in persistence, this would suggest that M. tuberculosis "
    "has evolved a regulatory code in which operator architecture encodes the propensity for "
    "phenotypic diversification."
)

add_body(
    "Question 3: What are the temporal properties of asymmetry-driven noise, and how much "
    "symmetrization would be needed to reduce persistence? This therapeutic question is also completely "
    "unexplored. Understanding not just how much noise the asymmetric operator generates, but how long "
    "that noise persists through time and across cell divisions, is essential for predicting the "
    "duration of persister dormancy and the efficacy of potential anti-noise interventions."
)

# 1.5
add_subsection_header("1.5 This Study")

add_body(
    "We address all three questions through a four-aim computational framework. Aim 1 performs "
    "genome-wide binding site discovery, asymmetry scoring, and persistence enrichment analysis. "
    "Aim 2 constructs a thermodynamically calibrated stochastic model of the Mce3R operator using "
    "Gillespie simulation, Berg-von Hippel energy calibration, and MCMC cooperativity inference. "
    "Aim 3 characterizes noise across amplitude, temporal dynamics, multi-generational lineage memory, "
    "and therapeutic quenching dose-response. Aim 4, designated as a stretch goal, builds an "
    "evolutionary simulation testing whether operator asymmetry is positively selected under "
    "persistence pressure."
)

add_body(
    "All results are computational predictions with explicit falsification criteria. The project "
    "bridges structural biology (Panagoda et al. 2024), persistence microbiology (Pandey et al. 2023), "
    "noise theory (Balazsi et al. 2011, Chowdhury et al. 2021), and the direct causal link between "
    "noise and persistence (Quigley & Lewis 2022) into an integrated framework that generates testable "
    "predictions for experimental validation."
)

page_break()

# =============================================================================
# 2. SPECIFIC AIMS
# =============================================================================
add_section_header("2. SPECIFIC AIMS")

# AIM 1
add_subsection_header("AIM 1: Discover Binding Sites, Score Genome-Wide Asymmetry, and Test Persistence Enrichment")

add_body_runs([
    ("Sub-aim 1a: Motif discovery and validation. ", True, False),
    ("We will perform de novo motif discovery using MEME (Multiple EM for Motif Elicitation) across "
     "the upstream regulatory regions of Mce3R-controlled genes in three mycobacterial genomes: "
     "M. tuberculosis H37Rv, M. bovis, and M. marinum. MEME will be run in ZOOPS (Zero or One "
     "Occurrence Per Sequence) mode, searching both strands, with motif widths ranging from 6 to "
     "110 base pairs. The resulting position weight matrices (PWMs) will be used for genome-wide "
     "scanning with FIMO (Find Individual Motif Occurrences) across the complete H37Rv genome at a "
     "p-value threshold of 1e-4. Cross-species conservation analysis will assess whether predicted "
     "sites are maintained across the three genomes, providing evidence for functional binding. "
     "Validation will be performed by confirming that the top FIMO hit recovers the known 123 bp "
     "Mce3R operator at its established genomic coordinates.", False, False),
])

add_body_runs([
    ("Sub-aim 1b: Genome-wide asymmetry scoring (NEW). ", True, False),
    ("For all 1,442 FIMO-predicted binding sites, we will compute a comprehensive set of architectural "
     "features designed to capture operator asymmetry. The primary metric is an asymmetry score "
     "defined as the absolute value of log2(PWM_score_site1 / PWM_score_site2) for paired binding "
     "sites occurring within a 150 base pair window, which identifies potential bipartite operators. "
     "Additional features include a palindromicity index (defined as 1 minus the reverse complement "
     "similarity score), spacer length between paired sites, total information content of the binding "
     "motif (sum of per-position information in bits), and conservation score (fraction identity "
     "across three species at each locus).", False, False),
])

add_body(
    "We will then cross-reference each site's nearest gene against published persistence gene lists, "
    "primarily the comprehensive Tn-seq screen of DeJesus et al. (2017), supplemented by published "
    "CRISPRi data. The central statistical test is a Mann-Whitney U test comparing asymmetry scores "
    "of persistence genes versus non-persistence genes. If persistence genes show significantly higher "
    "median asymmetry scores, this constitutes evidence that M. tuberculosis has evolved a persistence "
    "regulatory code in which high-asymmetry operators are systematically deployed at genes that "
    "benefit from noise-driven phenotypic diversification."
)

add_body_runs([
    ("Sub-aim 1c: Noise prediction for top candidates. ", True, False),
    ("Using the calibrated Gillespie model from Aim 2, we will simulate noise for the top 20-50 "
     "highest-asymmetry operators identified in Sub-aim 1b. For each operator, we will compute the "
     "predicted coefficient of variation (CV) and persister fraction, testing whether predicted noise "
     "correlates with persistence gene classification. This analysis connects the genome-wide "
     "architectural survey directly to functional noise predictions.", False, False),
])

# AIM 2
add_subsection_header("AIM 2: Construct a Thermodynamically Calibrated Stochastic Model")

add_body_runs([
    ("Sub-aim 2a: Four-state operator model. ", True, False),
    ("We will construct a stochastic model of the Mce3R operator based on four distinct binding "
     "states: state 0 (unbound, both sites empty), state 1 (strong site occupied, weak site empty), "
     "state 2 (weak site occupied, strong site empty), and state 3 (both sites occupied). The model "
     "comprises 12 reactions in the Gillespie stochastic simulation algorithm (SSA): 8 operator "
     "state transitions (binding and unbinding at each site from each state), plus transcription, "
     "translation, mRNA decay, and protein decay. We will simulate 50,000 cells under four primary "
     "conditions: asymmetric (wild-type Kd values), symmetric (both sites set to the geometric mean "
     "Kd of 10.84 nM), single-site (strong site only), and unregulated (constitutive expression). "
     "State-specific transcription rates reflect the degree of repression imposed by each occupancy "
     "state: state 0 produces 0.150 mRNA per minute (full expression), state 1 produces 0.0225 mRNA "
     "per minute (85 percent blocked by the strong site), state 2 produces 0.075 mRNA per minute "
     "(50 percent blocked by the weak site), and state 3 produces 0.01125 mRNA per minute "
     "(92.5 percent blocked by both sites).", False, False),
])

add_body_runs([
    ("Sub-aim 2b: Berg-von Hippel energy calibration. ", True, False),
    ("We will convert MEME-derived position weight matrix scores to physical binding energies using "
     "the Berg-von Hippel framework. Per-position binding energies are computed as dG_i = -kT times "
     "ln(f_i / p_i), where f_i is the observed frequency at position i in the PWM and p_i is the "
     "background frequency, evaluated at T = 310 K (physiological temperature). The energy scale is "
     "calibrated against the experimentally measured dissociation constants: dG_strong = -12.23 "
     "kcal/mol corresponding to Kd = 2.4 nM, and dG_weak = -10.37 kcal/mol corresponding to "
     "Kd = 49 nM. This calibration anchors all subsequent thermodynamic calculations to real "
     "biophysical measurements.", False, False),
])

add_body_runs([
    ("Sub-aim 2c: MCMC cooperativity inference. ", True, False),
    ("We will infer the cooperativity parameter omega and the spacer energy contribution dG_spacer "
     "using Markov Chain Monte Carlo (MCMC) sampling implemented with the emcee package. The "
     "ensemble sampler will use 32 walkers for 5,000 steps with 1,000 steps of burn-in. Convergence "
     "will be assessed using the Gelman-Rubin R-hat statistic (requiring R-hat < 1.01) and "
     "acceptance fraction. The posterior distributions of omega and dG_spacer will be reported, "
     "allowing us to propagate parameter uncertainty through all downstream noise predictions.", False, False),
])

add_body_runs([
    ("Sub-aim 2d: Operator classification. ", True, False),
    ("All FIMO-predicted operator sites from Aim 1 will be classified as either graded repressors "
     "or digital switches based on their effective Hill coefficient, computed from the thermodynamic "
     "model. A Hill coefficient below 2 indicates graded repression (analog-like response), while "
     "a Hill coefficient above 2 indicates switch-like (digital) behavior. This classification "
     "determines the functional character of each operator and its potential contribution to noise.", False, False),
])

# AIM 3
add_subsection_header("AIM 3: Characterize Noise Amplitude, Temporal Dynamics, Lineage Memory, and Quenching")

add_body_runs([
    ("Sub-aim 3a: Noise amplitude and robustness. ", True, False),
    ("We will compute the coefficient of variation (CV) across a comprehensive matrix of conditions: "
     "4 operator architectures (asymmetric, symmetric, single-site, unregulated) times 4 environments "
     "(baseline, cholesterol-rich, acidic pH, host-like combined) times 2 model types (single-gene "
     "and two-species autoregulatory). To ensure robustness to parameter uncertainty, we will "
     "propagate noise predictions through 50 MCMC posterior samples, testing whether CV_asymmetric "
     "exceeds CV_symmetric across the full posterior distribution. A TetR/tetO2 biological benchmark "
     "will provide an independent validation point. Mutual information analysis will quantify the "
     "information capacity of each architecture across 3 architectures and 10 Mce3R concentrations, "
     "testing the noise-information trade-off hypothesis. Persister fractions will be computed using "
     "a calibrated threshold matched to published wild-type persistence frequencies.", False, False),
])

add_body_runs([
    ("Sub-aim 3b: Temporal noise dynamics (NEW). ", True, False),
    ("We will characterize the temporal structure of noise by recording full expression traces from "
     "Gillespie simulations at 10-minute intervals for 15,000 minutes post-burn-in, yielding 1,500 "
     "time points per cell for 1,000 cells per architecture per environment. From these traces, we "
     "will compute the autocorrelation function C(tau) = <delta_x(t) times delta_x(t+tau)> / "
     "<delta_x squared>, computed via fast Fourier transform for efficiency. We will fit exponential "
     "decay models C(tau) approximately equal to exp(-tau/tau_c) to extract the autocorrelation time "
     "tau_c for each architecture, which quantifies how long noise fluctuations persist in time. "
     "We will compute dwell time distributions, defined as the durations of contiguous intervals "
     "during which protein levels remain below the persistence threshold. Power spectral density "
     "S(f) will be computed as the Fourier transform of C(tau), identifying characteristic switching "
     "frequencies. All metrics will be compared across asymmetric, symmetric, and single-site "
     "architectures.", False, False),
])

add_body(
    "The central prediction is that tau_c(asymmetric) will exceed tau_c(symmetric), meaning that "
    "noise fluctuations persist longer in the asymmetric architecture. The mechanism is that the "
    "strong site stays bound for extended periods (slow k_off = 0.04 per minute), while the weak "
    "site flickers rapidly (fast k_off = 0.82 per minute), creating asymmetric switching dynamics "
    "with extended periods in partially repressed states. If confirmed, this means asymmetry affects "
    "not just the probability of entering the persister state, but the duration of persister dormancy, "
    "which has direct implications for treatment duration."
)

add_body_runs([
    ("Sub-aim 3c: Multi-generational lineage memory (NEW). ", True, False),
    ("We will simulate cell division events at 1,500-minute intervals, reflecting the approximately "
     "25-hour doubling time of M. tuberculosis. At each division event, the daughter cell inherits "
     "half of the parent's protein count (modeled as binomial partitioning with probability 0.5 per "
     "molecule) plus the current operator binding state. We will track full lineage trees of 64-128 "
     "cells spanning 6-7 generations from single ancestor cells. Starting from 100 lineages per "
     "architecture, each initiated from a random ancestor state, we will measure the persistence "
     "lineage duration, defined as the number of generations until all descendants of a persister "
     "ancestor have exited the persister state.", False, False),
])

add_body(
    "We will also compute the fraction of terminal cells at generation 7 that remain in the persister "
    "state given a persister ancestor. The prediction is that the asymmetric operator produces longer "
    "persistence lineages because its wider protein distribution means that some daughters are born "
    "deep in the persister zone after binomial partitioning and take longer to escape. This analysis "
    "connects operator architecture to trans-generational persistence memory, a previously unmodeled "
    "dimension of tuberculosis latency."
)

add_body_runs([
    ("Sub-aim 3d: Noise quenching dose-response (NEW). ", True, False),
    ("We will generate therapeutic dose-response curves by systematically varying two parameters. "
     "First, a symmetrization sweep: the weak site Kd will be varied from 49 nM to 2.4 nM in 15 "
     "steps, while the strong site Kd remains fixed at 2.4 nM, progressively symmetrizing the "
     "operator. Second, a concentration sweep: Mce3R concentration will be varied from 10 to "
     "10,000 nM in 15 log-spaced steps for both architectures. For each combination, we will compute "
     "CV, persister fraction, and mean protein level from 5,000 simulated cells. From the "
     "symmetrization sweep, we will extract IC50_noise (the Kd_weak value at which CV drops to "
     "halfway between its maximum and minimum values) and IC50_persistence (the Kd_weak value at "
     "which persister fraction drops by 50 percent). Sigmoidal dose-response curves will be fitted "
     "to generate quantitative therapeutic targets.", False, False),
])

add_body(
    "Flentie et al. (2019) demonstrated that Mce3R-targeting compounds enhance antibiotic efficacy "
    "by 16 to 50-fold. Our model predicts the mechanism underlying this enhancement and provides "
    "quantitative predictions for the required dose of symmetrization to achieve clinically "
    "meaningful reductions in persister frequency. This represents the first quantitative prediction "
    "of noise-quenching dose-response in a tuberculosis regulatory system."
)

# AIM 4
add_subsection_header("AIM 4 (STRETCH): Evolutionary Simulation of Asymmetry Under Persistence Selection")

add_body(
    "As a stretch goal, we will build a fast surrogate fitness model by fitting a polynomial function "
    "CV(Kd_ratio, spacer, block) from existing sweep data, then wrap this surrogate in a genetic "
    "algorithm to simulate the evolution of operator architecture under different selection pressures. "
    "The population will consist of 200 virtual operators, each defined by mutable parameters: Kd "
    "ratio (range 1-100), spacer length (10-100 bp), and blocking fraction (0.3-0.95). Mutation "
    "will be modeled as Gaussian perturbation with sigma equal to 5 percent of the parameter range. "
    "Selection will retain the top 50 percent of the population by fitness, with the remainder "
    "regenerated through crossover and mutation."
)

add_body(
    "Two fitness functions will be applied. Fitness function A (persistence selection): fitness equals "
    "the fraction of cells below the expression threshold after a simulated antibiotic pulse, modeling "
    "an environment where persistence confers survival advantage. Fitness function B (growth-rate "
    "selection): fitness equals mean expression level, modeling an environment where faster growth "
    "is advantageous and persistence is irrelevant. We will run both fitness functions for 1,000 "
    "generations and compare the evolved Kd ratios."
)

add_body(
    "The prediction is that persistence selection will drive evolution toward higher operator "
    "asymmetry, while growth-rate selection will not. This will be tested separately for "
    "persistence-pathway genes and housekeeping genes. The novel claim is that asymmetric operator "
    "architecture is positively selected as a noise generator specifically under persistence pressure. "
    "This aim is labeled as exploratory: a positive result strengthens the overall project narrative "
    "significantly, while an inconclusive result does not weaken the core findings from Aims 1-3."
)

page_break()

# =============================================================================
# 3. METHODS
# =============================================================================
add_section_header("3. METHODS")

# 3.1
add_subsection_header("3.1 Genome Acquisition")

add_body(
    "Three mycobacterial genomes will be acquired from NCBI using BioPython (Cock et al. 2009): "
    "M. tuberculosis H37Rv (accession NC_000962.3, 4,411,532 bp, 65.6% GC content), M. bovis "
    "(accession NC_002945.4), and M. marinum (accession NC_010612.1). For motif discovery, 200 bp "
    "upstream of yrbE3A orthologs will be extracted from each genome using annotated gene coordinates. "
    "The known Mce3R operator is located at H37Rv coordinates 2,207,477-2,207,699, spanning 123 bp "
    "within the yrbE3A-mce3R intergenic region. Orthologous regions in M. bovis and M. marinum will "
    "be identified by reciprocal best BLAST hit of the yrbE3A coding sequence."
)

# 3.2
add_subsection_header("3.2 MEME Motif Discovery")

add_body(
    "De novo motif discovery will be performed using MEME Suite version 5.5.9 (Bailey et al. 2009). "
    "MEME will be run in ZOOPS (Zero or One Occurrence Per Sequence) mode, searching both DNA strands, "
    "with motif widths ranging from 6 to 110 base pairs. A zero-order Markov background model will be "
    "used with nucleotide frequencies reflecting the GC-rich M. tuberculosis genome: A = 0.204, "
    "C = 0.296, G = 0.296, T = 0.204. The top 5 motifs ranked by E-value will be retained for "
    "genome-wide scanning. Position weight matrices (PWMs) will be extracted from each motif for "
    "subsequent FIMO analysis."
)

# 3.3
add_subsection_header("3.3 FIMO Genome-Wide Scan")

add_body(
    "All 5 PWMs from MEME will be scanned against the complete M. tuberculosis H37Rv genome using "
    "FIMO (Find Individual Motif Occurrences; Grant et al. 2011). The p-value threshold will be set "
    "at 1e-4, balancing sensitivity against false positive rate. Output for each hit includes genomic "
    "coordinates, strand, PWM match score, p-value, and nearest gene annotation. A total of 1,442 "
    "candidate binding sites are expected based on preliminary analysis at this threshold. Each site "
    "will be annotated with its distance to the nearest gene start codon and the functional category "
    "of that gene."
)

# 3.4
add_subsection_header("3.4 Genome-Wide Asymmetry Scoring (NEW)")

add_body(
    "For each FIMO hit, we will determine whether paired binding sites exist within a 150 bp window, "
    "indicating a potential bipartite operator analogous to the Mce3R system. For paired sites, we "
    "will compute the following per-site features:"
)

add_body(
    "(a) Asymmetry score: defined as |log2(PWM_score_strong / PWM_score_weak)| for the two sites in "
    "each pair. A score of 0 indicates perfect symmetry; higher values indicate greater asymmetry. "
    "The Mce3R operator serves as a calibration point with a known asymmetry corresponding to its "
    "20.4-fold Kd difference."
)

add_body(
    "(b) Palindromicity index: defined as 1 minus (reverse complement alignment score / maximum "
    "possible score). TetR-family proteins typically bind palindromic sequences; deviations from "
    "palindromicity may correlate with binding asymmetry."
)

add_body(
    "(c) Spacer length: the distance in base pairs between the centers of paired binding sites. The "
    "Mce3R operator has a 53 bp spacer, which is unusually long for TetR-family systems."
)

add_body(
    "(d) Total information content: the sum of per-position information in bits across the binding "
    "motif, reflecting the specificity of the protein-DNA interaction."
)

add_body(
    "(e) Conservation score: the fraction of nucleotide identity across the three mycobacterial "
    "species at each locus, providing evidence for functional constraint."
)

add_body(
    "For unpaired sites (single FIMO hit with no nearby partner within 150 bp), the asymmetry score "
    "is defined as 0, reflecting a single-site architecture. Cross-referencing will match each site's "
    "nearest gene against persistence gene lists from the comprehensive Tn-seq screen of DeJesus "
    "et al. (2017) and published CRISPRi data. Statistical testing will use the Mann-Whitney U test "
    "for comparing asymmetry score distributions between persistence and non-persistence genes, and "
    "Fisher's exact test for enrichment of high-asymmetry sites (defined as the top quartile) near "
    "persistence genes versus expectation by chance."
)

# 3.5
add_subsection_header("3.5 Four-State Operator Model")

add_body(
    "The stochastic model represents the Mce3R operator as occupying one of four binding states: "
    "unbound (state 0), strong site bound only (state 1), weak site bound only (state 2), and both "
    "sites bound (state 3). The model comprises 12 reactions: 8 operator state transitions "
    "(binding and unbinding of Mce3R at each site from each accessible state) plus transcription, "
    "translation, mRNA degradation, and protein degradation."
)

add_body(
    "Kinetic parameters are derived from the literature. The on-rate for transcription factor binding "
    "is k_on = 0.0167 per nM per minute, based on Stormo and Zhao (2010). The off-rate for each site "
    "is computed as k_off = Kd times k_on, yielding k_off_strong = 0.04 per minute and k_off_weak = "
    "0.82 per minute. The translation rate is 0.5 proteins per mRNA per minute (Taniguchi et al. "
    "2010). The mRNA half-life is 9.5 minutes (Rustad et al. 2013), corresponding to a degradation "
    "rate of 0.073 per minute. The protein half-life is 1,500 minutes, corresponding to a degradation "
    "rate of 0.000462 per minute."
)

add_body(
    "State-specific transcription rates reflect the degree of repression: state 0 = 0.150 mRNA/min "
    "(full expression), state 1 = 0.0225 mRNA/min (85% blocked by strong site), state 2 = 0.075 "
    "mRNA/min (50% blocked by weak site), state 3 = 0.01125 mRNA/min (92.5% blocked by both sites). "
    "The blocking fractions are derived from the structural model of Mce3R-DNA interaction."
)

# 3.6
add_subsection_header("3.6 Gillespie SSA Implementation")

add_body(
    "The exact Gillespie stochastic simulation algorithm (Gillespie 1977) will be implemented in "
    "Python with Numba JIT compilation for computational efficiency. Each simulation runs 50,000 "
    "independent cells per condition, with a maximum simulation time of 30,000 minutes and a burn-in "
    "period of 15,000 minutes to ensure steady-state sampling. Four primary conditions will be "
    "simulated: Condition A (asymmetric, wild-type Kd values), Condition B (symmetric, both sites "
    "set to the geometric mean Kd = 10.84 nM), Condition C (single-site, strong site only), and "
    "Condition D (unregulated, constitutive expression). An additional Condition E performs an "
    "asymmetry sweep with 10 Kd ratios ranging from 1-fold to 50-fold, holding the geometric mean "
    "constant. All random number generation uses deterministic seeding (master_seed = 42 + offset) "
    "to ensure reproducibility."
)

# 3.7
add_subsection_header("3.7 Berg-von Hippel Calibration")

add_body(
    "Binding energies will be computed using the Berg-von Hippel framework (Stormo & Zhao 2010). "
    "Per-position binding energy contributions are calculated as dG_i = -kT times ln(f_i / p_i), "
    "where f_i is the observed nucleotide frequency at position i in the PWM, p_i is the background "
    "frequency, k is Boltzmann's constant, and T = 310 K (physiological temperature). The total "
    "binding free energy for each site is the sum of per-position contributions."
)

add_body(
    "The energy scale is calibrated against the experimentally measured Kd values: dG_strong = "
    "-12.23 kcal/mol corresponding to Kd = 2.4 nM, and dG_weak = -10.37 kcal/mol corresponding "
    "to Kd = 49 nM. The thermodynamic partition function for the four-state operator is: "
    "Z = 1 + exp(-dG_s/kT)[R] + exp(-dG_w/kT)[R] + omega times exp(-(dG_s + dG_w + dG_sp)/kT)[R]^2, "
    "where [R] is the free Mce3R concentration, omega is the cooperativity parameter, and dG_sp is "
    "the spacer energy contribution."
)

# 3.8
add_subsection_header("3.8 MCMC Inference")

add_body(
    "Cooperativity parameters will be inferred using the emcee ensemble MCMC sampler with 32 walkers "
    "running for 5,000 steps with 1,000 steps of burn-in. The free parameters are ln(omega) and "
    "dG_spacer, with flat priors on ln(omega) from -5 to 5 and dG_spacer from -5 to 5 kcal/mol. "
    "Convergence will be assessed by the Gelman-Rubin R-hat statistic (requiring R-hat < 1.01) and "
    "acceptance fraction (target range 0.2-0.8, achieved value 0.70). The posterior distribution will "
    "be propagated to noise predictions: 50 posterior samples, each simulated with 2,000 cells, will "
    "test whether CV_asymmetric exceeds CV_symmetric across the full parameter uncertainty."
)

# 3.9
add_subsection_header("3.9 Environmental Extensions")

add_body(
    "Four environmental conditions will be modeled by modifying Mce3R concentration and noise "
    "parameters: baseline (standard parameters), cholesterol-rich (2-fold increase in Mce3R "
    "concentration, reflecting mce3 operon induction during cholesterol utilization), acidic pH "
    "(0.7-fold Mce3R with 1.3-fold noise increase, reflecting phagosomal conditions), and host-like "
    "(combined cholesterol and pH effects). A two-species autoregulatory model extends the single-gene "
    "model by including Mce3R autoregulation (negative feedback on its own expression) plus repression "
    "of the target gene, yielding 18 reactions and dynamic free Mce3R concentration. The full "
    "simulation matrix comprises 10,000 cells times 3 architectures times 4 environments times "
    "2 models, totaling 24 conditions."
)

# 3.10
add_subsection_header("3.10 Temporal Noise Dynamics (NEW)")

add_body(
    "Full expression traces will be recorded at 10-minute intervals for 15,000 minutes post-burn-in, "
    "yielding 1,500 time points per cell. A total of 1,000 cells per architecture per environment "
    "will be simulated with full trace recording. The autocorrelation function C(tau) = "
    "<delta_x(t) times delta_x(t+tau)> / <delta_x^2> will be computed via fast Fourier transform "
    "for computational efficiency, where delta_x(t) = x(t) - <x> represents the deviation from the "
    "mean expression level."
)

add_body(
    "Exponential decay models C(tau) approximately equal to exp(-tau/tau_c) will be fitted to "
    "extract the autocorrelation time tau_c for each architecture. Dwell time distributions will "
    "be computed as the durations of contiguous intervals during which protein levels remain below "
    "the persistence threshold. Power spectral density S(f) will be computed as the Fourier transform "
    "of C(tau), identifying characteristic switching frequencies between high-expression and "
    "low-expression states. All temporal metrics will be compared across asymmetric, symmetric, "
    "and single-site architectures using bootstrap confidence intervals."
)

# 3.11
add_subsection_header("3.11 Multi-Generational Lineage Tracking (NEW)")

add_body(
    "Cell division will be simulated at 1,500-minute intervals, reflecting the standard M. tuberculosis "
    "doubling time of approximately 25 hours. At each division event, the daughter cell inherits "
    "protein molecules via binomial partitioning: each protein molecule in the parent cell is assigned "
    "to one of the two daughter cells with probability 0.5, modeled as a Binomial(parent_protein, 0.5) "
    "draw. The operator binding state is inherited directly by the daughter cell, reflecting the "
    "assumption that operator state equilibrates fast relative to cell division."
)

add_body(
    "Full lineage trees will be tracked from single ancestor cells for 7 generations, producing 128 "
    "terminal cells per lineage. A total of 100 lineages per architecture will be simulated, each "
    "initiated from a random ancestor state drawn from the steady-state distribution. The primary "
    "metric is persistence lineage duration, defined as the number of generations until all "
    "descendants of a persister ancestor have exited the persister state. We will also compute the "
    "fraction of terminal cells at generation 7 that remain in the persister state given that the "
    "founding ancestor was in the persister state. These metrics will be compared across architectures "
    "to test whether operator asymmetry extends the duration of trans-generational persistence memory."
)

# 3.12
add_subsection_header("3.12 Noise Quenching Dose-Response (NEW)")

add_body(
    "Two systematic parameter sweeps will generate dose-response curves. The symmetrization sweep "
    "varies the weak site Kd from 49 nM to 2.4 nM in 15 evenly spaced steps, while the strong site "
    "Kd remains fixed at 2.4 nM. This progressively eliminates operator asymmetry, converting the "
    "bipartite operator into a symmetric one. A total of 5,000 cells will be simulated per parameter "
    "point, recording CV, persister fraction, and mean protein level."
)

add_body(
    "The concentration sweep varies Mce3R concentration from 10 to 10,000 nM in 15 logarithmically "
    "spaced steps for both the asymmetric and symmetric architectures, mapping how repressor "
    "abundance modulates noise and persistence across the full dynamic range."
)

add_body(
    "From the symmetrization sweep data, we will extract two therapeutic metrics: IC50_noise, "
    "defined as the Kd_weak value at which CV drops to the midpoint between CV_asymmetric and "
    "CV_symmetric; and IC50_persistence, defined as the Kd_weak value at which the persister "
    "fraction drops by 50 percent. Sigmoidal dose-response curves will be fitted using the "
    "four-parameter model: CV(Kd_weak) = CV_min + (CV_max - CV_min) / (1 + (Kd_weak / IC50)^n), "
    "where n is the Hill coefficient of the dose-response."
)

# 3.13
add_subsection_header("3.13 Mutual Information")

add_body(
    "Mutual information between environmental state and gene expression will be computed as "
    "MI = H(expression) - H(expression | environment), where H denotes Shannon entropy. Four "
    "environmental states will be defined (baseline, cholesterol, acidic, host-like), with 5,000 "
    "cells simulated per environment. Entropy will be estimated using a k-nearest neighbor (k-NN) "
    "estimator to avoid binning artifacts. Mutual information will be computed for each of 3 "
    "architectures (asymmetric, symmetric, single-site) across 10 Mce3R concentrations spanning "
    "the physiological range, testing how architecture modulates the trade-off between noise and "
    "environmental information capacity."
)

# 3.14
add_subsection_header("3.14 Persistence Threshold")

add_body(
    "The persister threshold will be defined using three complementary methods to ensure robustness. "
    "Method 1 (tail_fraction): the threshold is set at the 1st percentile of the unregulated "
    "expression distribution, identifying cells in the extreme low-expression tail. Method 2 "
    "(absolute_calibrated): the threshold is calibrated to match the published wild-type persister "
    "frequency of approximately 0.1 percent in unstressed M. tuberculosis populations. Method 3 "
    "(fold_change): the threshold is set at the population median divided by 1.5, identifying cells "
    "substantially below average expression. Sensitivity analysis will be performed across multiple "
    "percentile cutoffs to verify that the qualitative conclusions (asymmetric produces more persisters "
    "than symmetric) are robust to threshold choice."
)

# 3.15
add_subsection_header("3.15 Evolutionary Simulation (STRETCH)")

add_body(
    "The evolutionary simulation uses a fast surrogate fitness model to avoid running full Gillespie "
    "simulations at each generation. The surrogate is constructed by fitting a polynomial function "
    "CV(Kd_ratio, spacer, block) from Condition E sweep data supplemented by additional grid points "
    "spanning the parameter space. The genetic algorithm maintains a population of 200 virtual "
    "operators, each defined by three mutable parameters: Kd ratio (range 1-100), spacer length "
    "(10-100 bp), and blocking fraction (0.3-0.95)."
)

add_body(
    "Mutation is modeled as Gaussian perturbation with standard deviation equal to 5 percent of "
    "the parameter range per parameter per generation. Selection retains the top 50 percent of the "
    "population by fitness; the remaining slots are filled by crossover (uniform crossover between "
    "two randomly selected parents from the surviving population) plus mutation."
)

add_body(
    "Fitness function A (persistence selection) defines fitness as the survival fraction, computed "
    "as the fraction of cells with expression below the persister threshold after a simulated "
    "antibiotic pulse that kills all cells above the threshold. Fitness function B (growth-rate "
    "selection) defines fitness as mean expression level, where higher expression corresponds to "
    "faster growth. Both functions are run for 1,000 generations, and the evolved Kd ratio "
    "distributions are compared between the two selection regimes."
)

# 3.16
add_subsection_header("3.16 Statistical Analysis")

add_body(
    "Bootstrap resampling with 1,000 iterations will be used to construct 95 percent confidence "
    "intervals for all noise metrics. Distribution comparisons will use the Kolmogorov-Smirnov (KS) "
    "test with Cohen's d for effect size estimation. Genome-wide asymmetry enrichment will be tested "
    "using the Mann-Whitney U test (comparing asymmetry scores between persistence and non-persistence "
    "genes) and Fisher's exact test (comparing the proportion of high-asymmetry sites near persistence "
    "genes versus expected by chance). Gaussian mixture models (GMM) with 1-3 components will be "
    "fitted to protein distributions and compared by Bayesian Information Criterion (BIC) to test for "
    "bimodality. All p-values are two-sided with alpha = 0.05."
)

# 3.17
add_subsection_header("3.17 Software")

add_body(
    "All analyses will be implemented in Python 3.10 using the following packages: NumPy and SciPy "
    "for numerical computation and statistical testing, Pandas for data management, Matplotlib and "
    "Seaborn for visualization, BioPython for genome acquisition and sequence analysis, scikit-learn "
    "for Gaussian mixture modeling, statsmodels for statistical tests, Numba for JIT compilation of "
    "Gillespie simulations, and MEME Suite 5.5.9 for motif discovery and scanning. All random number "
    "generators are seeded deterministically (master seed = 42) to ensure complete reproducibility."
)

page_break()

# =============================================================================
# 4. EXPECTED RESULTS AND SIGNIFICANCE
# =============================================================================
add_section_header("4. EXPECTED RESULTS AND SIGNIFICANCE")

# 4.1
add_subsection_header("4.1 Aim 1: Binding Sites and Genome-Wide Asymmetry")

add_body(
    "The FIMO genome-wide scan is expected to validate the known Mce3R operator as the top-scoring "
    "hit with a p-value of 1.1e-11, confirming the accuracy of our motif discovery pipeline. Across "
    "the H37Rv genome, we expect to identify 1,442 candidate binding sites that will be scored for "
    "the full suite of asymmetry features described in the methods."
)

add_body(
    "The central prediction of Aim 1 is that persistence genes, as defined by published Tn-seq "
    "screens (DeJesus et al. 2017), will have significantly higher median asymmetry scores than "
    "non-persistence genes, tested by Mann-Whitney U (predicted p < 0.05). If this prediction is "
    "confirmed, it constitutes evidence that M. tuberculosis has evolved a persistence regulatory "
    "code: an asymmetry code in which high-asymmetry operators are systematically deployed at genes "
    "that benefit from noise-driven phenotypic diversification. This would represent a fundamentally "
    "new organizing principle for understanding how M. tuberculosis encodes its persistence strategy "
    "in DNA regulatory architecture."
)

add_body(
    "If the prediction is not confirmed, meaning that persistence genes show no enrichment for "
    "high-asymmetry operators, this indicates that the asymmetry observed at Mce3R is a "
    "system-specific feature rather than a genome-wide principle. This negative result is still "
    "informative, as it constrains the scope of operator asymmetry as a noise-generating mechanism "
    "to specific regulatory systems rather than a global strategy."
)

add_body(
    "We also expect to identify novel candidate targets near genes of known significance, including "
    "tgs1 (triacylglycerol synthase, involved in dormancy lipid body formation), espC (ESX-1 "
    "secretion system component, critical for virulence), and mbtG (mycobactin biosynthesis, "
    "essential for iron acquisition). These candidates suggest an expanded Mce3R regulon beyond "
    "the canonical mce3 operon. Conservation across three species supports functional binding, "
    "although the close evolutionary relationship between M. tuberculosis and M. bovis limits the "
    "power of conservation analysis to detect asymmetry-specific constraint."
)

# 4.2
add_subsection_header("4.2 Aim 2: Thermodynamic Model")

add_body(
    "The thermodynamic model is expected to reproduce the experimentally measured Kd values exactly, "
    "validating the Berg-von Hippel calibration. The MCMC inference will yield omega = 1.08 with a "
    "95 percent credible interval of 0.20-5.89, indicating that the two Mce3R binding sites exhibit "
    "no significant cooperativity, meaning they bind independently. This finding has important "
    "mechanistic implications: independent binding creates more distinct intermediate operator states "
    "(states 1 and 2), each with its own transcription rate, generating more graded noise than a "
    "cooperative system that switches sharply between fully bound and fully unbound states."
)

add_body(
    "All 20 classified operators from FIMO sites are expected to function as graded repressors "
    "with Hill coefficients below 2, consistent with the absence of strong cooperativity. The "
    "asymmetric operator will show a broader repression transition region compared to symmetric "
    "architectures, meaning that repression is achieved gradually over a wider range of Mce3R "
    "concentrations rather than switching sharply at a critical threshold."
)

# 4.3
add_subsection_header("4.3 Aim 3a: Noise Amplitude")

add_body(
    "The asymmetric operator is expected to produce a coefficient of variation of 0.187, compared to "
    "0.157 for the symmetric operator, yielding a delta-CV of 0.031 (p < 0.001 by bootstrap test). "
    "This 19 percent increase in noise is robust across 100 percent of MCMC posterior samples, the "
    "two-species autoregulatory model, and all four environmental conditions."
)

add_body(
    "Corrected persister fractions using the calibrated threshold show 0.28 percent for the asymmetric "
    "architecture versus 0.09 percent for the symmetric architecture under baseline conditions, "
    "representing a 3.1-fold enrichment. Under acidic pH conditions mimicking the phagosomal "
    "environment, these fractions increase to 2.65 percent versus 1.23 percent, a 2.15-fold ratio "
    "that translates to a 29.4-fold enrichment in absolute persister excess over the symmetric "
    "architecture."
)

add_body(
    "Mutual information analysis reveals that the asymmetric operator has MI = 0.28 bits, compared "
    "to 0.44 bits for the symmetric operator. This noise-information trade-off is consistent with "
    "a bet-hedging strategy (Veening et al. 2008): the asymmetric operator sacrifices the ability to "
    "faithfully encode environmental information in exchange for phenotypic diversification that "
    "provides a survival advantage under unpredictable antibiotic exposure."
)

add_body(
    "Context is important for interpreting these results. A CV of 0.187 represents moderate noise, "
    "approximately 25 percent above typical regulated genes and 42 percent of the noise level "
    "observed in hipA toxin-antitoxin modules, which are known to be primary noise-driven persistence "
    "switches. Operator asymmetry provides a constitutive noise floor that is always present, "
    "rather than serving as the dominant noise source. Other mechanisms, including toxin-antitoxin "
    "modules, metabolic fluctuations, and stochastic partitioning of growth regulators, likely "
    "contribute more absolute noise. The significance of operator asymmetry is its architectural "
    "permanence: it is encoded in DNA and does not require activation."
)

# 4.4
add_subsection_header("4.4 Aim 3b: Temporal Dynamics (KEY NOVEL RESULT)")

add_body(
    "We predict that the autocorrelation time tau_c will be significantly longer for the asymmetric "
    "architecture than for the symmetric architecture. The mechanism underlying this prediction is "
    "the kinetic asymmetry between the two binding sites: the strong site, with k_off = 0.04 per "
    "minute, remains bound for an average of 25 minutes per binding event, while the weak site, with "
    "k_off = 0.82 per minute, remains bound for an average of only 1.2 minutes. This creates "
    "asymmetric switching dynamics in which the operator spends extended periods in partially "
    "repressed states (particularly state 1, where only the strong site is bound), generating "
    "correlated fluctuations that persist longer than in the symmetric case where both sites have "
    "similar kinetics."
)

add_body(
    "The clinical implication of longer autocorrelation time is profound: if noise fluctuations "
    "persist longer, then cells that stochastically enter the persister state remain dormant longer "
    "before returning to the active state. This means that operator asymmetry affects not just the "
    "probability of entering the persister state (reflected in the persister fraction), but the "
    "duration of persister dormancy, which directly determines treatment duration. A small increase "
    "in noise amplitude combined with a large increase in noise autocorrelation time could have "
    "disproportionate effects on treatment outcomes."
)

add_body(
    "We acknowledge a significant risk: the autocorrelation time may be dominated by the protein "
    "half-life of 1,500 minutes, which sets a fundamental timescale for expression fluctuations "
    "regardless of operator architecture. If protein half-life dominates, tau_c will be similar "
    "across architectures, and the temporal dynamics claim will be falsified. This negative result "
    "would still be informative, establishing that operator architecture modulates noise amplitude "
    "but not noise memory at the protein level."
)

# 4.5
add_subsection_header("4.5 Aim 3c: Lineage Memory (KEY NOVEL RESULT)")

add_body(
    "We predict that the asymmetric operator produces longer persistence lineages, meaning that more "
    "generations are required for all descendants of a persister ancestor to exit the persister state. "
    "The mechanism is that binomial partitioning at cell division creates variable protein inheritance: "
    "each daughter cell receives approximately half of the parent's protein molecules, but with "
    "substantial binomial variation. Because the asymmetric operator generates a wider protein "
    "distribution (higher CV), some daughters of a persister cell are born with very low protein "
    "levels, deep in the persister zone, and require more time and more cell divisions to accumulate "
    "sufficient protein to exit the persister state."
)

add_body(
    "This analysis connects operator architecture to trans-generational persistence memory, explaining "
    "how a regulatory DNA feature could contribute to the extraordinary duration of tuberculosis "
    "latency. While latency is clearly a multi-factorial phenomenon involving immune containment, "
    "granuloma formation, and metabolic adaptation, the contribution of regulatory architecture to "
    "the intrinsic duration of persister states has not been previously modeled. Our lineage tracking "
    "provides the first quantitative framework for this dimension of persistence."
)

add_body(
    "The primary risk is that partitioning noise may dominate over operator-driven noise in "
    "determining lineage memory duration. If binomial partitioning creates sufficient variation to "
    "mask the operator architecture effect, the asymmetric and symmetric architectures may show "
    "similar lineage memory. Both outcomes are publishable and informative."
)

# 4.6
add_subsection_header("4.6 Aim 3d: Noise Quenching (KEY TRANSLATIONAL RESULT)")

add_body(
    "We predict a sigmoidal dose-response relationship between weak-site Kd and both CV and persister "
    "fraction, with the persister fraction declining more steeply than CV due to the nonlinear "
    "relationship between distribution width and tail probability. The IC50_noise value will "
    "identify the specific Kd_weak at which noise is reduced to the midpoint between asymmetric "
    "and symmetric levels, providing a quantitative therapeutic target: any intervention that "
    "shifts the weak-site affinity to this value would achieve half-maximal noise reduction."
)

add_body(
    "The translational significance is grounded in existing experimental work. Flentie et al. (2019) "
    "demonstrated that small molecules targeting the Mce3R system enhance antibiotic efficacy by "
    "16 to 50-fold in M. tuberculosis. Our model provides a mechanistic explanation for this "
    "enhancement: compounds that increase Mce3R affinity for the weak site (symmetrizing the "
    "operator) would reduce noise, decrease the persister subpopulation, and thereby sensitize the "
    "remaining population to antibiotics. The dose-response curves generated by our model provide "
    "quantitative predictions for the required degree of symmetrization to achieve clinically "
    "meaningful persistence reduction."
)

add_body(
    "This represents the first quantitative prediction of noise-quenching dose-response in a "
    "tuberculosis regulatory system, bridging computational noise theory with drug development."
)

# 4.7
add_subsection_header("4.7 Aim 4: Evolutionary Simulation (STRETCH)")

add_body(
    "Three outcomes are possible from the evolutionary simulation. First, if operator asymmetry "
    "increases under persistence selection but not under growth-rate selection, this supports the "
    "hypothesis that asymmetry is positively selected specifically as a noise generator for "
    "persistence. This would provide evolutionary evidence that the 20.4-fold asymmetry observed "
    "at Mce3R is not accidental but has been actively maintained by natural selection."
)

add_body(
    "Second, if asymmetry increases under both selection regimes, this suggests that high-asymmetry "
    "operators have additional functional advantages beyond noise generation, such as enabling more "
    "graded environmental responses. This outcome would weaken the specificity of the "
    "noise-for-persistence argument but would still demonstrate that asymmetry is under positive "
    "selection."
)

add_body(
    "Third, if no directional change in asymmetry is observed under either selection regime, this "
    "suggests that the observed asymmetry at Mce3R is selectively neutral, possibly arising from "
    "mutational drift in the weak binding site. This outcome would weaken the adaptive argument "
    "but would not affect the core findings from Aims 1-3 regarding the quantitative relationship "
    "between asymmetry and noise. As a stretch goal, this aim strengthens the project if positive "
    "but does not weaken it if inconclusive."
)

page_break()

# =============================================================================
# 5. LIMITATIONS AND FALSIFICATION CRITERIA
# =============================================================================
add_section_header("5. LIMITATIONS AND FALSIFICATION CRITERIA")

add_subsection_header("5.1 Limitations")

add_body(
    "1. All results are computational predictions and have not been experimentally validated. "
    "The quantitative predictions generated by our models are hypothesis-generating, not "
    "hypothesis-confirming. Experimental validation through single-cell fluorescence reporters, "
    "operator mutagenesis, and persister assays is required."
)

add_body(
    "2. The stochastic model represents a single gene (or at most two genes in the autoregulatory "
    "model) from a 14-gene operon. The mce3 operon contains multiple co-transcribed genes whose "
    "expression may be subject to additional regulatory mechanisms including translational coupling, "
    "mRNA secondary structure, and post-translational regulation. Our model does not capture these "
    "complexities."
)

add_body(
    "3. The MCMC posterior for the cooperativity parameter omega spans a wide range (0.20-5.89), "
    "reflecting limited calibration data. While the noise predictions are robust across this full "
    "posterior range, the precise value of omega and its mechanistic interpretation should be treated "
    "as approximate. Additional experimental data on cooperative binding would substantially narrow "
    "this posterior."
)

add_body(
    "4. Environmental conditions are modeled as scalar multipliers on Mce3R concentration and noise "
    "amplitude, not as dynamic signals that change over time. In reality, the phagosomal environment "
    "fluctuates as macrophages respond to bacterial cues, creating a time-varying regulatory "
    "landscape that our steady-state analysis does not capture."
)

add_body(
    "5. The symmetric control is a theoretical construct (geometric mean of the two Kd values) that "
    "does not correspond to any known natural operator. The TetR/tetO2 system provides a biological "
    "benchmark for symmetric binding, but it differs from Mce3R in many other respects. Direct "
    "comparison is therefore approximate."
)

add_body(
    "6. The conservation argument is limited by the close evolutionary relationship between "
    "M. tuberculosis and M. bovis (which share over 99.95 percent sequence identity). M. marinum "
    "provides a more distant outgroup, but the power to detect asymmetry-specific conservation "
    "remains limited with only three genomes."
)

add_body(
    "7. The CV difference of 0.031 (delta-CV = 19 percent relative increase) is moderate. Other "
    "molecular mechanisms, including toxin-antitoxin modules such as hipBA and vapBC, metabolic "
    "fluctuations in central carbon metabolism, and stochastic partitioning of growth-rate "
    "determinants at cell division, likely contribute more absolute noise to the persistence "
    "decision. Our work characterizes one component of a multi-component noise landscape."
)

add_body(
    "8. The autocorrelation time tau_c and lineage memory duration may be dominated by the protein "
    "half-life (1,500 minutes) or by binomial partitioning noise at cell division, respectively. "
    "If these non-operator mechanisms dominate, the temporal and lineage memory predictions will "
    "show minimal architecture dependence, limiting these analyses to negative results."
)

add_body(
    "9. The genome-wide enrichment test depends on the quality of both the Tn-seq persistence gene "
    "annotations (which include some false positives and false negatives) and the accuracy of FIMO "
    "site pairing (which assumes that paired sites within 150 bp represent bipartite operators). "
    "Errors in either input will reduce the power or specificity of the enrichment analysis."
)

add_subsection_header("5.2 Falsification Criteria")

add_body(
    "1. If experimentally symmetrizing the Mce3R operator (e.g., by mutating the weak site to match "
    "the strong site sequence) shows no measurable change in expression noise or persister frequency "
    "in single-cell assays, the core claim that asymmetry generates functionally relevant noise "
    "would be falsified."
)

add_body(
    "2. If adding upstream noise sources (such as fluctuations in RNA polymerase availability or "
    "growth rate variation) to the model eliminates the CV difference between asymmetric and "
    "symmetric architectures, this would indicate that operator-level noise is overwhelmed by "
    "extrinsic noise sources and is not functionally relevant."
)

add_body(
    "3. If strong cooperativity (omega substantially greater than 1) is established with improved "
    "experimental data, the four-state model would need fundamental revision, as strong cooperativity "
    "would effectively reduce the four-state operator to a two-state switch, potentially eliminating "
    "the noise-generating intermediate states."
)

add_body(
    "4. If the two-species autoregulatory model shows CV_asymmetric less than or equal to "
    "CV_symmetric, meaning that autoregulatory feedback reverses the noise difference, then the "
    "simple model's prediction would not extend to the more realistic regulatory context. This has "
    "already been tested and the difference persists, but it remains a falsifiable prediction for "
    "experimental systems."
)

add_body(
    "5. If the autocorrelation time tau_c shows no significant architecture dependence, meaning "
    "that the asymmetric and symmetric operators produce noise with similar temporal correlation "
    "structure, the temporal memory claim of Aim 3b would be falsified. However, the amplitude claim "
    "from Aim 3a would remain intact."
)

add_body(
    "6. If persistence genes show no enrichment for high-asymmetry operators in the genome-wide "
    "analysis, the persistence regulatory code hypothesis would be falsified as a genome-wide "
    "principle. However, the single-system claims about the Mce3R operator would remain valid."
)

page_break()

# =============================================================================
# 6. INNOVATION
# =============================================================================
add_section_header("6. INNOVATION")

add_body("This project makes five novel contributions to the field:")

add_body(
    "First, it provides the first genome-wide analysis of operator binding asymmetry as a predictor "
    "of persistence gene function in M. tuberculosis. By scoring all predicted operator sites for "
    "asymmetry and testing for enrichment at persistence-critical loci, we test whether a persistence "
    "regulatory code is encoded in operator DNA architecture. No previous study has examined whether "
    "the arrangement of binding site affinities within regulatory regions correlates with persistence "
    "function at a genomic scale."
)

add_body(
    "Second, it delivers the first quantitative noise characterization, spanning amplitude, temporal "
    "dynamics, and quenching response, of a structurally resolved asymmetric bacterial operator with "
    "real cryo-EM-derived dissociation constants. Previous noise studies have relied on hypothetical "
    "parameter values; our model is anchored to the actual biophysical measurements from Panagoda "
    "et al. (2024)."
)

add_body(
    "Third, it provides the first prediction of noise temporal memory and multi-generational lineage "
    "persistence from operator architecture. By tracking how noise propagates through cell division "
    "lineages, we connect regulatory DNA structure to the duration of persister dormancy, a dimension "
    "of persistence that has not been previously modeled at the operator level."
)

add_body(
    "Fourth, it generates the first therapeutic dose-response curve for noise quenching in a "
    "tuberculosis regulatory system. The IC50 values for both noise reduction and persister fraction "
    "reduction provide quantitative symmetrization targets that can guide the development of "
    "anti-noise therapeutics, complementing the existing experimental work of Flentie et al. (2019)."
)

add_body(
    "Fifth, as a stretch goal, it presents the first evolutionary simulation testing whether operator "
    "asymmetry is positively selected under persistence pressure, addressing the fundamental question "
    "of whether the 20.4-fold asymmetry at Mce3R is an adaptive feature or a neutral byproduct."
)

add_body(
    "Together, these contributions bridge structural biology (Panagoda et al. 2024), persistence "
    "microbiology (Pandey et al. 2023, Quigley & Lewis 2022), and noise theory (Chowdhury et al. "
    "2021, Balazsi et al. 2011) into an integrated framework that generates testable, quantitative "
    "predictions."
)

page_break()

# =============================================================================
# 7. TIMELINE
# =============================================================================
add_section_header("7. TIMELINE")

# Build table
table = doc.add_table(rows=9, cols=4)
table.style = 'Table Grid'

# Header row
headers = ['Month', 'Phase', 'Tasks', 'Deliverables']
for i, h in enumerate(headers):
    cell = table.rows[0].cells[i]
    cell.text = ''
    p = cell.paragraphs[0]
    run = p.add_run(h)
    run.bold = True
    run.font.size = Pt(10)
    run.font.name = 'Calibri'
    p.paragraph_format.line_spacing = 1.5

rows_data = [
    ('1', 'Aim 1a',
     'Motif discovery, genome scanning, conservation analysis',
     '1,442 sites validated'),
    ('2', 'Aim 2a-2b',
     'Gillespie simulation, thermodynamic calibration',
     '50K cells x 4 conditions, energy parameters'),
    ('3', 'Aim 2c-2d + 1b',
     'MCMC inference, operator classification, genome-wide asymmetry scoring',
     'Omega posterior, asymmetry scores for all sites'),
    ('4', 'Aim 1b-1c + 3a',
     'Persistence enrichment test, noise simulation of top candidates, environmental extensions',
     'Enrichment p-value, 24-condition noise data'),
    ('5', 'Aim 3b',
     'Temporal dynamics and lineage tracking',
     'Tau_c values, dwell time distributions, lineage memory'),
    ('6', 'Aim 3c-3d',
     'Noise quenching and mutual information',
     'IC50 values, dose-response curves, MI trade-off'),
    ('7', 'Aim 4 (stretch)',
     'Evolutionary simulation',
     'Selection curves (if time permits)'),
    ('8', 'Manuscript',
     'Figures (19-24 total), manuscript, competition materials',
     'Complete submission package'),
]

for row_idx, (month, phase, tasks, deliverables) in enumerate(rows_data, start=1):
    for col_idx, text in enumerate([month, phase, tasks, deliverables]):
        cell = table.rows[row_idx].cells[col_idx]
        cell.text = ''
        p = cell.paragraphs[0]
        run = p.add_run(text)
        run.font.size = Pt(10)
        run.font.name = 'Calibri'
        p.paragraph_format.line_spacing = 1.5

# Set column widths
for row in table.rows:
    row.cells[0].width = Inches(0.7)
    row.cells[1].width = Inches(1.3)
    row.cells[2].width = Inches(3.0)
    row.cells[3].width = Inches(2.0)

page_break()

# =============================================================================
# 8. REFERENCES
# =============================================================================
add_section_header("8. REFERENCES")

references = [
    "1. WHO (2024). Global Tuberculosis Report. World Health Organization, Geneva.",
    "2. Balaban NQ, Helaine S, Lewis K, et al. (2019). Definitions and guidelines for research on antibiotic persistence. Nature Reviews Microbiology, 17:441-448.",
    "3. Elowitz MB, Levine AJ, Siggia ED, Swain PS (2002). Stochastic gene expression in a single cell. Science, 297:1183-1186.",
    "4. Chowdhury D, Bhattacherjee A, Bhargava R (2021). Modeling gene expression noise in the context of cis-regulatory architecture. Frontiers in Genetics, 12:698910.",
    "5. Lengyel IM, Morelli LG (2017). Multiple binding sites for transcription factors tune the expression noise. Physical Review E, 95:042412.",
    "6. Balazsi G, van Oudenaarden A, Collins JJ (2011). Cellular decision making and biological noise: from microbes to mammals. Cell, 144:910-925.",
    "7. Quigley J, Lewis K (2022). Noise in a metabolic pathway leads to persister formation in Mycobacterium tuberculosis. Microbiology Spectrum, 10:e0094822.",
    "8. Santangelo MP, Goldstein J, Alito A, et al. (2002). Negative transcriptional regulation of the mce3 operon in Mycobacterium tuberculosis. Microbiology, 148:2997-3006.",
    "9. Dunphy KY, Senaratne RH, Masuzawa M, Kendall LV, Riley LW (2010). Attenuation of Mycobacterium tuberculosis functionally disrupted in a fatty acyl-CoA synthetase gene fadD5. Journal of Infectious Diseases, 201:1232-1239.",
    "10. Pandey AK, Sassetti CM (2023). Mycobacterial persistence requires the utilization of host cholesterol. Research in Microbiology, 174:104082.",
    "11. Panagoda GJ, DeJesus MA, Engholm-Keller K, et al. (2024). Structural basis of asymmetric DNA recognition by the Mce3R repressor of Mycobacterium tuberculosis. ACS Chemical Biology, 19:2580-2592.",
    "12. Cock PJA, Antao T, Chang JT, et al. (2009). Biopython: freely available Python tools for computational molecular biology and bioinformatics. Bioinformatics, 25:1422-1423.",
    "13. Bailey TL, Boden M, Buske FA, et al. (2009). MEME Suite: tools for motif discovery and searching. Nucleic Acids Research, 37:W202-W208.",
    "14. Grant CE, Bailey TL, Noble WS (2011). FIMO: scanning for occurrences of a given motif. Bioinformatics, 27:1017-1018.",
    "15. DeJesus MA, Gerrick ER, Xu W, et al. (2017). Comprehensive essentiality analysis of the Mycobacterium tuberculosis genome via saturating transposon mutagenesis. mBio, 8:e02133-16.",
    "16. Gillespie DT (1977). Exact stochastic simulation of coupled chemical reactions. Journal of Physical Chemistry, 81:2340-2361.",
    "17. Stormo GD, Zhao Y (2010). Determining the specificity of protein-DNA interactions. Nature Reviews Genetics, 11:751-760.",
    "18. Taniguchi Y, Choi PJ, Li GW, et al. (2010). Quantifying E. coli proteome and transcriptome with single-molecule sensitivity in single cells. Science, 329:533-538.",
    "19. Rustad TR, Minch KJ, Ma S, et al. (2013). Mapping and manipulating the Mycobacterium tuberculosis transcriptome using a transcription factor overexpression-derived regulatory network. Nucleic Acids Research, 41:509-517.",
    "20. Santangelo MP, Blanco FC, Bianco MV, et al. (2009). Study of the role of Mce3R on the transcription of mce genes of Mycobacterium tuberculosis. Microbiology, 155:882-891.",
    "21. Flentie K, Harrison GA, Tukenmez H, et al. (2019). Chemical disarming of isoniazid resistance in Mycobacterium tuberculosis. ACS Infectious Diseases, 5:1239-1254.",
    "22. Veening JW, Smits WK, Kuipers OP (2008). Bistability, epigenetics, and bet-hedging in bacteria. Annual Review of Microbiology, 62:193-210.",
    "23. Rosenfeld N, Young JW, Alon U, Swain PS, Elowitz MB (2005). Gene regulation at the single-cell level. Science, 307:1962-1965.",
    "24. Farquhar KS, Charlebois DA, Szenk M, et al. (2019). Role of network-mediated stochasticity in mammalian drug resistance. Nature Communications, 10:2766.",
    "25. Rotem E, Loinger A, Ronin I, et al. (2010). Regulation of phenotypic variability by a threshold-based mechanism underlies bacterial persistence. Proceedings of the National Academy of Sciences, 107:12541-12546.",
    "26. Wakamoto Y, Dhar N, Chait R, et al. (2013). Dynamic persistence of antibiotic-stressed mycobacteria. Science, 339:91-95.",
    "27. Sureka K, Ghosh B, Dasgupta A, et al. (2008). Positive feedback and noise activate the stringent response regulator Rel in mycobacteria. PLoS ONE, 3:e1771.",
]

for ref in references:
    p = doc.add_paragraph()
    p.paragraph_format.line_spacing = 2.0
    p.paragraph_format.left_indent = Inches(0.5)
    p.paragraph_format.first_line_indent = Inches(-0.5)
    run = p.add_run(ref)
    run.font.size = Pt(11)
    run.font.name = 'Calibri'

# =============================================================================
# SAVE
# =============================================================================
doc.save(OUTPUT_PATH)
print(f"Saved to: {OUTPUT_PATH}")
print(f"File size: {os.path.getsize(OUTPUT_PATH):,} bytes")
