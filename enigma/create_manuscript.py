#!/usr/bin/env python3
"""Generate ENIGMA manuscript as a Word document."""
import docx
from docx import Document
from docx.shared import Pt, Inches, RGBColor
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.enum.table import WD_TABLE_ALIGNMENT
import os

doc = Document()

# Page setup
for section in doc.sections:
    section.top_margin = Inches(1)
    section.bottom_margin = Inches(1)
    section.left_margin = Inches(1.25)
    section.right_margin = Inches(1.25)

style = doc.styles['Normal']
font = style.font
font.name = 'Times New Roman'
font.size = Pt(12)
font.color.rgb = RGBColor(0, 0, 0)

# Line spacing
from docx.shared import Pt as PtShared
from docx.oxml.ns import qn
style.paragraph_format.space_after = Pt(6)
style.paragraph_format.line_spacing = 1.15

def h1(text):
    p = doc.add_heading(text, level=1)
    for run in p.runs:
        run.font.name = 'Times New Roman'
        run.font.size = Pt(14)
        run.font.color.rgb = RGBColor(0, 0, 0)
        run.bold = True
    return p

def h2(text):
    p = doc.add_heading(text, level=2)
    for run in p.runs:
        run.font.name = 'Times New Roman'
        run.font.size = Pt(12)
        run.font.color.rgb = RGBColor(0, 0, 0)
        run.bold = True
        run.italic = True
    return p

def para(text, bold=False, italic=False, indent=False):
    p = doc.add_paragraph()
    run = p.add_run(text)
    run.font.name = 'Times New Roman'
    run.font.size = Pt(12)
    run.bold = bold
    run.italic = italic
    if indent:
        p.paragraph_format.first_line_indent = Inches(0.5)
    return p

def mixed_para(*segments, indent=False):
    """segments = list of (text, bold, italic) tuples"""
    p = doc.add_paragraph()
    for text, bold, italic in segments:
        run = p.add_run(text)
        run.font.name = 'Times New Roman'
        run.font.size = Pt(12)
        run.bold = bold
        run.italic = italic
    if indent:
        p.paragraph_format.first_line_indent = Inches(0.5)
    return p

def add_table(headers, rows, col_widths=None):
    table = doc.add_table(rows=1 + len(rows), cols=len(headers))
    table.style = 'Table Grid'
    table.alignment = WD_TABLE_ALIGNMENT.CENTER
    for i, h in enumerate(headers):
        cell = table.rows[0].cells[i]
        cell.text = h
        for p in cell.paragraphs:
            p.alignment = WD_ALIGN_PARAGRAPH.CENTER
            for run in p.runs:
                run.bold = True
                run.font.name = 'Times New Roman'
                run.font.size = Pt(10)
    for r_idx, row in enumerate(rows):
        for c_idx, val in enumerate(row):
            cell = table.rows[r_idx + 1].cells[c_idx]
            cell.text = str(val)
            for p in cell.paragraphs:
                p.alignment = WD_ALIGN_PARAGRAPH.CENTER
                for run in p.runs:
                    run.font.name = 'Times New Roman'
                    run.font.size = Pt(10)
    if col_widths:
        for i, w in enumerate(col_widths):
            for row in table.rows:
                row.cells[i].width = Inches(w)
    return table

# ══════════════════════════════════════════════════════════════
# TITLE
# ══════════════════════════════════════════════════════════════
p = doc.add_paragraph()
p.alignment = WD_ALIGN_PARAGRAPH.CENTER
run = p.add_run("Asymmetric Operator Architecture as a Molecular Mechanism for\nNoise-Driven Antibiotic Persistence in Mycobacterium tuberculosis")
run.font.name = 'Times New Roman'
run.font.size = Pt(16)
run.bold = True

p = doc.add_paragraph()
p.alignment = WD_ALIGN_PARAGRAPH.CENTER
run = p.add_run("A Computational Study of the Mce3R Regulatory System")
run.font.name = 'Times New Roman'
run.font.size = Pt(12)
run.italic = True

p = doc.add_paragraph()
p.alignment = WD_ALIGN_PARAGRAPH.CENTER

doc.add_page_break()

# ══════════════════════════════════════════════════════════════
# ABSTRACT
# ══════════════════════════════════════════════════════════════
h1("Abstract")

para(
    "Tuberculosis remains the leading cause of death from a single infectious agent worldwide, "
    "killing over 1.2 million people annually. A major barrier to effective treatment is antibiotic "
    "persistence, in which genetically susceptible bacteria enter a drug-tolerant state that survives "
    "therapy and drives relapse. Unlike genetic resistance, persistence arises from phenotypic "
    "heterogeneity: genetically identical cells randomly adopt different expression states, and rare "
    "outlier cells survive. What creates this heterogeneity at the molecular level remains poorly "
    "understood.",
    indent=True
)

para(
    "Here, we investigate the transcriptional repressor Mce3R (Rv1963c), which controls cholesterol "
    "metabolism genes critical for survival within host macrophages. Mce3R binds an unusual "
    "non-palindromic operator containing two sites of dramatically different affinity (Kd = 2.4 nM "
    "vs. 49 nM), a 20-fold asymmetry unprecedented in the TetR protein family. We develop a "
    "seven-phase computational framework spanning three aims to determine whether this asymmetry is "
    "functionally significant.",
    indent=True
)

mixed_para(
    ("In Aim 1, ", False, False),
    ("we discover and validate conserved Mce3R binding sites ", False, False),
    ("using de novo motif discovery (MEME) and genome-wide scanning (FIMO) across three mycobacterial "
     "species. The top FIMO hit (p = 1.1 x 10", False, False),
    ("-11", False, True),
    (") matches the known operator exactly, and 1,442 candidate sites are identified genome-wide, "
     "including near dormancy (tgs1), virulence (espC), and iron acquisition (mbtG) genes. ", False, False),
    ("In Aim 2, ", False, False),
    ("we model how operator architecture controls repression ", False, False),
    ("by building a four-state stochastic simulation (Gillespie algorithm, 50,000 cells) calibrated "
     "with a thermodynamic partition function (Berg-von Hippel model) and Bayesian cooperativity "
     "inference (MCMC). All 20 top-ranked operators classify as graded repressors (Hill coefficient "
     "< 2), not digital switches. ", False, False),
    ("In Aim 3, ", False, False),
    ("we quantify how asymmetry generates persistence-relevant noise. ", False, False),
    ("The native asymmetric operator produces a coefficient of variation (CV) of 0.187, compared to "
     "0.157 for an equivalent symmetric operator, a statistically significant 19% increase "
     "(bootstrap 95% CI entirely above zero, KS p < 0.001, Cohen's d = -2.29). Under simulated "
     "host-like stress conditions combining cholesterol exposure and acidic pH, the asymmetric operator "
     "maintains higher noise across all four tested environments.", False, False),
    indent=True
)

para(
    "Information-theoretic analysis reveals that the symmetric operator transmits more mutual "
    "information about environmental state (0.44 bits vs. 0.28 bits at physiological concentration), "
    "suggesting that the asymmetric operator sacrifices information capacity for increased stochastic "
    "variation. We propose that this trade-off is an evolved bet-hedging strategy: by generating a "
    "broader expression distribution with a heavier low-expression tail, the asymmetric operator "
    "increases the fraction of cells that stochastically cross a persistence threshold, providing a "
    "molecular mechanism linking DNA-binding site architecture to antibiotic persistence.",
    indent=True
)

doc.add_page_break()

# ══════════════════════════════════════════════════════════════
# INTRODUCTION
# ══════════════════════════════════════════════════════════════
h1("Introduction")

h2("The Problem: Tuberculosis Persistence")

para(
    "Tuberculosis (TB), caused by Mycobacterium tuberculosis, remains the world's deadliest "
    "infectious disease, killing over 1.2 million people per year despite being curable with "
    "existing antibiotics (WHO, 2024). Standard treatment requires 6-9 months of combination "
    "therapy, and even with full compliance, relapse occurs in 5-10% of patients. A central driver "
    "of treatment failure is antibiotic persistence: a subpopulation of genetically susceptible "
    "bacteria enters a drug-tolerant dormant state that survives therapy (Balaban et al., 2019). "
    "Unlike resistance, which arises from genetic mutations, persistence is phenotypic. Genetically "
    "identical cells stochastically adopt different expression states, and rare cells in a "
    "low-metabolic-activity state survive antibiotic exposure. Understanding the molecular origins "
    "of this phenotypic heterogeneity is essential for developing strategies to shorten TB treatment.",
    indent=True
)

h2("The System: Mce3R and Its Unusual Operator")

para(
    "Mce3R (Rv1963c) is a TetR-family transcriptional repressor that controls the mce3 operon "
    "(Rv1964-Rv1977), encoding cholesterol and lipid transport machinery essential for M. tuberculosis "
    "survival within host macrophages (Santangelo et al., 2002; Dunphy et al., 2010). Cholesterol "
    "import is not merely a nutritional convenience; it is required for persistence within the "
    "phagosomal environment where fatty acids serve as the primary carbon source. Deletion of mce3R "
    "increases antibiotic persister frequency (Pandey et al., 2023), directly linking this regulator "
    "to the persistence phenotype.",
    indent=True
)

para(
    "Recent cryo-EM structural work revealed that Mce3R is an unprecedented double TetR-fold repeat "
    "protein that binds a 123-base-pair non-palindromic operator containing two binding sites of "
    "dramatically different affinity (Panagoda, Balazsi & Sampson, 2024). The downstream (strong) "
    "site binds with Kd = 2.4 +/- 0.7 nM, while the upstream (weak) site binds with Kd of "
    "approximately 49 nM, creating a 20-fold affinity asymmetry separated by a 53-base-pair spacer. "
    "This architecture is highly unusual: most TetR-family repressors bind palindromic (symmetric) "
    "operators, and the functional significance of the asymmetry was unknown.",
    indent=True
)

h2("The Theoretical Framework: Noise and Persistence")

para(
    "The connection between gene regulatory architecture and phenotypic heterogeneity has been "
    "established through pioneering work on transcriptional noise. Balazsi, van Oudenaarden, and "
    "Collins (2011) demonstrated that circuit architecture, not just expression level, controls "
    "the magnitude and character of cell-to-cell variation in gene expression. High-noise circuits "
    "produce rare outlier cells that can survive environmental stress, a strategy termed "
    "bet-hedging. Farquhar et al. (2019) showed that network-mediated stochasticity drives drug "
    "resistance in mammalian cells, extending this framework beyond bacteria.",
    indent=True
)

para(
    "The critical gap in our understanding is mechanistic: while it is known that Mce3R deletion "
    "increases persistence (Pandey et al., 2023) and that circuit architecture controls noise "
    "(Balazsi et al., 2011), no one has connected the specific structural feature of the Mce3R "
    "operator (its binding site asymmetry) to the noise theory of persistence. This study addresses "
    "that gap.",
    indent=True
)

h2("Central Hypothesis")

para(
    "We hypothesize that the asymmetric dual-site operator of Mce3R creates four distinct promoter "
    "occupancy states (unbound, strong-only bound, weak-only bound, double-bound) with different "
    "transcription rates. This expanded state space generates more gene expression noise than an "
    "equivalent symmetric operator, increasing the fraction of cells that stochastically cross a "
    "persistence threshold. The asymmetric architecture is not an evolutionary accident; it is a "
    "regulatory strategy that tunes noise to enable phenotypic bet-hedging under host stress.",
    indent=True
)

doc.add_page_break()

# ══════════════════════════════════════════════════════════════
# METHODS
# ══════════════════════════════════════════════════════════════
h1("Methods")

para(
    "All analyses were performed using custom Python code with NumPy, SciPy, and Numba for "
    "acceleration. Parameters were drawn from published experimental measurements (Table 1). "
    "The complete pipeline is organized into seven computational phases, each building on the "
    "outputs of the previous.",
    indent=True
)

h2("Aim 1: Conserved Binding Site Discovery (Phase 1)")

para(
    "To identify Mce3R binding sites without prior bias, we performed de novo motif discovery "
    "using the MEME algorithm (Bailey et al., 2009). We extracted 200 base pairs upstream of the "
    "yrbE3A gene (Rv1964) from three mycobacterial genomes: M. tuberculosis H37Rv (NC_000962.3, "
    "4,411,532 bp, 65.6% GC), M. bovis AF2122/97 (NC_002945.4), and M. marinum M (NC_010612.1). "
    "These three orthologous upstream regions served as input to MEME under the ZOOPS (Zero Or One "
    "Occurrence Per Sequence) model, searching both DNA strands with motif widths from 6 to 110 bp "
    "and a 0-order Markov background model reflecting the high GC composition (A = 0.204, "
    "C = 0.296, G = 0.296, T = 0.204). The top 5 motifs were ranked by E-value.",
    indent=True
)

para(
    "The resulting position weight matrices (PWMs) were then scanned across the complete H37Rv "
    "genome using FIMO (Find Individual Motif Occurrences; Grant et al., 2011) with a p-value "
    "threshold of 1 x 10^-4. Each hit was annotated with the nearest gene, distance to the gene, "
    "and whether it falls in an intergenic region. Cross-species conservation was assessed by "
    "searching for orthologous sequences at corresponding genomic positions in M. bovis and "
    "M. marinum using a sliding-window approach with >80% sequence identity as the conservation "
    "threshold.",
    indent=True
)

h2("Aim 2: Stochastic Simulation of Operator Architecture (Phases 2 and 5)")

para(
    "We modeled the Mce3R operator as a four-state system in which the promoter exists in one of "
    "four configurations: unbound (U), strong-site-only bound (S), weak-site-only bound (W), or "
    "double-bound (D). Each state has a distinct transcription rate determined by the fractional "
    "blockage of each site (block_strong = 0.85, block_weak = 0.50). The system includes 12 "
    "chemical reactions: 8 operator binding/unbinding transitions (kon and koff for each site from "
    "each accessible state), plus transcription, translation, mRNA degradation, and protein "
    "degradation.",
    indent=True
)

para(
    "Binding kinetics were parameterized from experimental measurements: kon = 0.0167 nM^-1 min^-1 "
    "(Stormo & Zhao, 2010), Kd_strong = 2.4 nM and Kd_weak = 49 nM (Panagoda et al., 2024), "
    "yielding koff_strong = 0.040 min^-1 and koff_weak = 0.818 min^-1. Transcription and "
    "translation rates were k_max = 0.15 mRNA/min and k_translation = 0.5 protein/mRNA/min "
    "(Taniguchi et al., 2010). The mRNA half-life was set to 9.5 min (Rustad et al., 2013), giving "
    "gamma_mRNA = 0.073 min^-1. Protein degradation was dominated by dilution with a half-life of "
    "1,500 min (approximately 25 hours), giving gamma_protein = 0.000462 min^-1.",
    indent=True
)

para(
    "We simulated 50,000 independent cells per condition using the Gillespie stochastic simulation "
    "algorithm (SSA; Gillespie, 1977), implemented in Numba-compiled Python for performance. Each "
    "cell was simulated for 30,000 minutes (approximately 25 doubling times) with the first 15,000 "
    "minutes discarded as burn-in to ensure steady-state sampling. Four conditions were compared: "
    "(A) native asymmetric operator, (B) symmetric operator with geometric-mean Kd = 10.84 nM and "
    "geometric-mean block = 0.652 (ensuring equivalent total regulatory capacity), (C) single-site "
    "operator (weak site disabled, Kd_weak = 10^12 nM), and (D) no regulation (kon = 0).",
    indent=True
)

para(
    "To ground the stochastic simulation in physical chemistry, we calibrated binding energies "
    "using the Berg-von Hippel biophysical model (Phase 5). Each position in the MEME-discovered "
    "PWM was converted to a per-position binding energy using DeltaG_i = -kT x ln(f_i / bg_i), "
    "where f_i is the PWM frequency and bg_i is the genome background frequency. The total binding "
    "energy for each site was calibrated to reproduce the experimental Kd values exactly "
    "(DeltaG_strong = -12.23 kcal/mol, DeltaG_weak = -10.37 kcal/mol). A statistical mechanics "
    "partition function was constructed over the four operator states, and the cooperativity "
    "parameter omega and spacer energy penalty DeltaG_spacer were inferred using Markov Chain Monte "
    "Carlo (MCMC) with 32 walkers, 5,000 steps, and 1,000-step burn-in.",
    indent=True
)

h2("Aim 3: Statistical Analysis and Environmental Extensions (Phases 3 and 6)")

para(
    "Noise was quantified using the coefficient of variation (CV = standard deviation / mean), "
    "which normalizes for differences in mean expression across conditions. Bootstrap resampling "
    "(10,000 iterations) was used to construct 95% confidence intervals for the difference in CV "
    "between conditions. The Kolmogorov-Smirnov (KS) test assessed whether protein distributions "
    "differed significantly, and Cohen's d measured effect size. Gaussian mixture models (GMMs) with "
    "1, 2, and 3 components were fit to each condition's protein distribution using expectation-"
    "maximization, and the Bayesian Information Criterion (BIC) was used for model selection. The "
    "Fano factor (variance / mean) was computed as an additional noise metric.",
    indent=True
)

para(
    "To test whether the noise advantage of asymmetry persists under physiologically relevant "
    "conditions, we extended the model in Phase 6 to include four environments: baseline (standard "
    "growth), cholesterol-rich (increased Mce3R effective concentration), acidic pH (simulating "
    "macrophage phagosome conditions), and host-like (combined cholesterol and acidic stress). Each "
    "environment modifies the effective Mce3R concentration via a multiplier and adds a noise "
    "scaling factor to the transcription rate. We also built a two-species model in which Mce3R "
    "autoregulates its own expression (negative feedback) while simultaneously repressing the target "
    "gene, enabling us to capture the anti-correlation between regulator and target levels.",
    indent=True
)

para(
    "Mutual information I(environment; expression) was computed by discretizing protein expression "
    "distributions into bins and calculating the reduction in entropy of the output given knowledge "
    "of the input environment, across 10 Mce3R concentrations spanning 0.01 to 10,000 nM. A "
    "persistence threshold was defined as the 10th percentile of the unregulated condition's protein "
    "distribution (2,057 molecules), and the persister fraction was computed as the proportion of "
    "cells exceeding this threshold under each condition.",
    indent=True
)

doc.add_page_break()

# ══════════════════════════════════════════════════════════════
# TABLE 1: PARAMETERS
# ══════════════════════════════════════════════════════════════
h2("Table 1: Model Parameters")

add_table(
    ["Parameter", "Value", "Source"],
    [
        ["Kd_strong", "2.4 nM", "Panagoda et al. 2024"],
        ["Kd_weak", "49 nM", "Panagoda et al. 2024"],
        ["k_on", "0.0167 nM^-1 min^-1", "Stormo & Zhao 2010"],
        ["k_off_strong", "0.040 min^-1", "Derived (Kd x k_on)"],
        ["k_off_weak", "0.818 min^-1", "Derived (Kd x k_on)"],
        ["k_max (transcription)", "0.15 mRNA/min", "Mtb transcriptomics"],
        ["k_translation", "0.5 protein/mRNA/min", "Taniguchi et al. 2010"],
        ["mRNA half-life", "9.5 min", "Rustad et al. 2013"],
        ["Protein half-life", "1,500 min (~25 hr)", "Dilution-dominated"],
        ["block_strong", "0.85", "Estimated from expression data"],
        ["block_weak", "0.50", "Estimated from expression data"],
        ["Cells per condition", "50,000", "This study"],
        ["Simulation time", "30,000 min", "25 doubling times"],
        ["Burn-in", "15,000 min", "10x protein half-life"],
    ],
    col_widths=[2.0, 2.0, 2.5]
)

doc.add_page_break()

# ══════════════════════════════════════════════════════════════
# RESULTS AND DISCUSSION
# ══════════════════════════════════════════════════════════════
h1("Results and Discussion")

h2("Aim 1: De Novo Motif Discovery Validates the Known Operator and Reveals Novel Targets")

para(
    "MEME analysis of orthologous yrbE3A upstream sequences from three mycobacterial species "
    "identified two conserved motifs: an 8-bp core (ACATCAWA) and a 15-bp extended consensus "
    "(TWTKCATTGYTWTYT). When these motifs were scanned across the complete H37Rv genome using "
    "FIMO, 1,442 candidate binding sites were identified above the significance threshold "
    "(p < 10^-4).",
    indent=True
)

para(
    "The top-ranked FIMO hit (p = 1.1 x 10^-11, score = 34.86 bits) mapped to position 2,207,551 "
    "on the H37Rv chromosome, precisely within the known Mce3R operator region (2,207,477-2,207,699). "
    "The second-ranked hit (p = 3.1 x 10^-11) mapped 23 bp away at position 2,207,528, corresponding "
    "to the second binding site within the operator. This independent computational rediscovery of "
    "the experimentally characterized operator, without any prior structural information as input, "
    "validates the motif discovery pipeline.",
    indent=True
)

para(
    "Among novel candidate sites, several are biologically compelling. A site near Rv3130c (tgs1, "
    "p = 1.6 x 10^-6) is particularly notable because tgs1 encodes triacylglycerol synthase, the "
    "enzyme responsible for lipid body formation, a hallmark of dormant persister cells. If Mce3R "
    "regulates tgs1, it would provide a direct mechanistic link between cholesterol sensing and "
    "dormancy. Additional sites near Rv3615c (espC, an ESX-1 secretion component) and Rv2378c "
    "(mbtG, mycobactin biosynthesis) suggest that Mce3R may coordinate cholesterol metabolism with "
    "virulence factor secretion and iron acquisition. Cross-species conservation analysis confirmed "
    "greater than 80% sequence identity at the primary operator across all three species, indicating "
    "that the asymmetric architecture has been maintained under selective pressure over millions of "
    "years of mycobacterial evolution.",
    indent=True
)

h2("Aim 2: The Operator Functions as a Graded Analog Dimmer, Not a Digital Switch")

para(
    "Thermodynamic calibration using the Berg-von Hippel model reproduced the experimental Kd values "
    "exactly (Kd_strong = 2.4 nM, Kd_weak = 49 nM), corresponding to binding free energies of "
    "DeltaG_strong = -12.23 kcal/mol and DeltaG_weak = -10.37 kcal/mol (DeltaDeltaG = 1.86 kcal/mol). "
    "MCMC inference yielded a cooperativity parameter omega = 1.08 (95% CI: 0.20-5.89) and a spacer "
    "energy penalty DeltaG_spacer = -0.15 kcal/mol, indicating negligible cooperativity between the "
    "two binding sites. The wide posterior distribution for omega reflects an honestly underdetermined "
    "parameter given the available calibration data (two Kd values and population-level mean and CV), "
    "and the near-unity point estimate suggests that the two sites bind independently.",
    indent=True
)

para(
    "Classification of all 20 top-ranked FIMO operators using the Hill coefficient derived from "
    "thermodynamic repression curves revealed that every site functions as a graded repressor "
    "(Hill n_H < 2). None qualify as digital switches (n_H > 2). This means the Mce3R system "
    "operates as a continuous analog dimmer of gene expression rather than a binary on/off switch. "
    "The repression curves show that the native asymmetric operator has a broader transition region "
    "than the symmetric alternative, with more intermediate states accessible at physiological Mce3R "
    "concentrations (approximately 332 nM, corresponding to approximately 200 molecules per cell). "
    "This broad transition region is precisely the regime where small fluctuations in Mce3R "
    "concentration cause large changes in promoter state occupancy, the mechanistic origin of "
    "increased noise.",
    indent=True
)

h2("Aim 3: Asymmetry Generates Significantly More Expression Noise")

para(
    "Gillespie stochastic simulation of 50,000 cells per condition demonstrated that the native "
    "asymmetric operator (Condition A) produces a coefficient of variation of 0.187, compared to "
    "0.157 for the symmetric control (Condition B), 0.149 for the single-site control "
    "(Condition C), and 0.059 for the unregulated control (Condition D). The 19% increase in CV "
    "from asymmetry (A vs. B) is statistically significant by every metric: bootstrap 95% "
    "confidence interval for DeltaCV is entirely above zero, KS test p < 0.001, and Cohen's d = "
    "-2.29 (large effect). The asymmetry sweep (Condition E) confirmed that CV increases "
    "monotonically with the Kd ratio from 0.158 at ratio 1 (symmetric) to 0.185 at ratio 50, while "
    "holding geometric-mean Kd constant to control for total regulatory capacity.",
    indent=True
)

add_table(
    ["Condition", "Architecture", "Mean Protein", "CV", "Fano Factor"],
    [
        ["A", "Asymmetric (native)", "197.3", "0.187", "6.87"],
        ["B", "Symmetric (geom. mean)", "292.6", "0.157", "7.18"],
        ["C", "Single-site only", "341.8", "0.149", "7.63"],
        ["D", "No regulation", "2,224.0", "0.059", "7.77"],
    ],
    col_widths=[1.0, 2.0, 1.2, 0.8, 1.0]
)

para("")  # spacer

para(
    "Gaussian mixture model analysis revealed that a 2-component mixture fits the regulated "
    "conditions (A, B, C) significantly better than a 1-component model by BIC, confirming bimodal "
    "expression. The low-expression component corresponds to cells in a high-repression state "
    "where both operator sites are occupied, consistent with a persister-prone subpopulation with "
    "reduced cholesterol metabolism. The model predicts an 11.3-fold repression of yrbE3A expression, "
    "compared to the experimentally measured 8.5-fold repression reported by Santangelo et al. "
    "(2009), a ratio of 1.33 that is within the expected range given the simplified single-gene "
    "model.",
    indent=True
)

h2("Environmental Extensions Confirm Robustness of the Noise Advantage")

para(
    "The two-species environmental model (Phase 6) confirmed that CV(asymmetric) > CV(symmetric) "
    "in all four tested environments: baseline, cholesterol-rich, acidic pH, and host-like combined "
    "stress. Under host-like conditions, mean expression increased for all architectures "
    "(derepression), consistent with the biological expectation that Mce3R releases the operator "
    "under intracellular stress to allow cholesterol import. The Mce3R and target protein levels "
    "showed a negative Pearson correlation (r = -0.057 to -0.141 depending on condition), "
    "confirming that the autoregulatory negative feedback loop is functional in the two-species "
    "model.",
    indent=True
)

para(
    "Mutual information analysis produced an unexpected and scientifically interesting result: "
    "the symmetric operator transmits more information about environmental state than the "
    "asymmetric operator (MI = 0.44 bits vs. 0.28 bits at the physiological Mce3R concentration "
    "of 332 nM). This implies that the asymmetric operator is not optimized for faithful "
    "environmental sensing. Instead, it appears to sacrifice information fidelity for increased "
    "stochastic variation, generating a broader distribution of expression levels that includes "
    "more cells in the extreme tails. We interpret this as evidence for an evolved trade-off: the "
    "asymmetric operator prioritizes bet-hedging (producing rare persister cells) over accurate "
    "signal transduction (faithfully tracking environmental changes).",
    indent=True
)

h2("Implications for TB Biology and Drug Development")

para(
    "These findings have several implications for understanding and combating TB persistence. "
    "First, they identify the operator architecture itself, not just the regulator protein, as a "
    "determinant of persistence-relevant noise. This shifts the therapeutic target from Mce3R "
    "protein function to the DNA-level logic of the regulatory circuit. Second, the identification "
    "of candidate Mce3R binding sites near tgs1 (dormancy), espC (virulence), and mbtG (iron "
    "acquisition) suggests that Mce3R may coordinate a broader stress response program than "
    "previously recognized. Third, the information-noise trade-off we observe provides a general "
    "principle: asymmetric operators may be a common evolutionary strategy for tuning persistence "
    "in pathogens, and similar architectures should be sought in other TetR-family systems.",
    indent=True
)

para(
    "Flentie et al. (2019) demonstrated that 6-azasteroid compounds targeting the Mce3R pathway "
    "enhance isoniazid activity 16-fold and bedaquiline activity approximately 50-fold. Our model "
    "predicts that drugs which specifically disrupt the weak binding site (rather than both sites "
    "equally) would reduce noise without eliminating regulation, potentially converting the "
    "graded-response system into a more digital, lower-noise switch that produces fewer persister "
    "cells.",
    indent=True
)

h2("Limitations")

para(
    "Several limitations should be noted. First, the model uses a simplified single-gene or "
    "two-gene representation of what is in reality a multi-gene operon. Second, the MCMC-inferred "
    "cooperativity parameter has a wide posterior (95% CI: 0.20-5.89), reflecting limited "
    "calibration data. Third, environmental conditions are modeled as scalar multipliers rather "
    "than dynamic time-varying signals. Fourth, the persistence threshold is defined "
    "computationally rather than measured experimentally. Future work should include experimental "
    "validation of the predicted noise differences using single-cell fluorescence reporters in "
    "M. tuberculosis carrying symmetric vs. asymmetric operator variants.",
    indent=True
)

doc.add_page_break()

# ══════════════════════════════════════════════════════════════
# REFERENCES
# ══════════════════════════════════════════════════════════════
h1("References")

refs = [
    "1. Bailey TL, et al. (2009). MEME Suite: tools for motif discovery and searching. Nucleic Acids Research 37:W202-W208.",
    "2. Balaban NQ, et al. (2019). Definitions and guidelines for research on antibiotic persistence. Nature Reviews Microbiology 17:441-448.",
    "3. Balazsi G, van Oudenaarden A, Collins JJ. (2011). Cellular decision making and biological noise: from microbes to mammals. Cell 144:910-925.",
    "4. Browning DF, Busby SJW. (2004). The regulation of bacterial transcription initiation. Nature Reviews Microbiology 2:57-65.",
    "5. Cock PJA, et al. (2009). Biopython: freely available Python tools for computational molecular biology and bioinformatics. Bioinformatics 25:1422-1423.",
    "6. Cuthbertson L, Nodwell JR. (2013). The TetR family of regulators. Microbiology and Molecular Biology Reviews 77:440-475.",
    "7. Dunphy KY, et al. (2010). Attenuation of Mycobacterium tuberculosis functionally disrupted in a fatty acyl-CoA synthetase gene fadD5. Journal of Infectious Diseases 201:1232-1239.",
    "8. Farquhar KS, et al. (2019). Role of network-mediated stochasticity in mammalian drug resistance. Nature Communications 10:2766.",
    "9. Flentie K, et al. (2019). Chemical disarming of isoniazid resistance in Mycobacterium tuberculosis. ACS Infectious Diseases 5:1239-1254.",
    "10. Gillespie DT. (1977). Exact stochastic simulation of coupled chemical reactions. Journal of Physical Chemistry 81:2340-2361.",
    "11. Grant CE, Bailey TL, Noble WS. (2011). FIMO: scanning for occurrences of a given motif. Bioinformatics 27:1017-1018.",
    "12. Panagoda GDR, Balazsi G, Sampson NS. (2024). Cryo-EM structure of Mce3R bound to the mce3 operator reveals an atypical TetR family member. ACS Chemical Biology 19:2580-2592.",
    "13. Pandey AK, et al. (2023). Deletion of mce3R increases antibiotic persister frequency. Research in Microbiology 174:104082.",
    "14. Pandey AK, Bhatt A. (2023). Mce transporters and their role in Mycobacterium tuberculosis pathogenesis. Tuberculosis.",
    "15. Rustad TR, et al. (2013). Global analysis of mRNA stability in Mycobacterium tuberculosis. Nucleic Acids Research 41:509-517.",
    "16. Santangelo MP, et al. (2002). Characterization of Mce3R, a TetR-type transcriptional repressor of mce3 operon. Microbiology 148:2997-3006.",
    "17. Santangelo MP, et al. (2008). Mce3R, a TetR-type transcriptional repressor, controls the expression of a regulon involved in lipid metabolism. BMC Microbiology 8:38.",
    "18. Santangelo MP, et al. (2009). Study of the role of Mce3R in the transcription of mce genes. Microbiology 155:882-891.",
    "19. Stormo GD, Zhao Y. (2010). Determining the specificity of protein-DNA interactions. Nature Reviews Genetics 11:751-760.",
    "20. Sureka K, et al. (2008). Positive feedback and noise activate the stringent response regulator Rel in mycobacteria. PLoS ONE 3:e1771.",
    "21. Taniguchi Y, et al. (2010). Quantifying E. coli proteome and transcriptome with single-molecule sensitivity in single cells. Science 329:533-538.",
    "22. WHO. (2024). Global Tuberculosis Report 2024. World Health Organization, Geneva.",
]

for ref in refs:
    p = doc.add_paragraph()
    run = p.add_run(ref)
    run.font.name = 'Times New Roman'
    run.font.size = Pt(10)
    p.paragraph_format.space_after = Pt(4)
    p.paragraph_format.left_indent = Inches(0.5)
    p.paragraph_format.first_line_indent = Inches(-0.5)

# ── Save ──
outpath = "/Users/aayanalwani/tb project/mce3r_stochastic/ENIGMA_Manuscript.docx"
doc.save(outpath)
print(f"Saved to {outpath}")
print(f"Size: {os.path.getsize(outpath) / 1024:.1f} KB")
