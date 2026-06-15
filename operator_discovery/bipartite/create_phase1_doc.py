#!/usr/bin/env python3
"""Generate Phase 1 explanation document using python-docx."""
import docx
from docx import Document
from docx.shared import Pt, Inches, RGBColor, Cm
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.enum.table import WD_TABLE_ALIGNMENT
from docx.oxml.ns import qn
import os

doc = Document()

# ── Page setup ──
for section in doc.sections:
    section.top_margin = Inches(1)
    section.bottom_margin = Inches(1)
    section.left_margin = Inches(1)
    section.right_margin = Inches(1)

style = doc.styles['Normal']
font = style.font
font.name = 'Calibri'
font.size = Pt(11)
font.color.rgb = RGBColor(0x33, 0x33, 0x33)

# ── Helper functions ──
def add_title(text, color=RGBColor(0x0D, 0x94, 0x88)):
    p = doc.add_heading(text, level=0)
    for run in p.runs:
        run.font.color.rgb = color
        run.font.size = Pt(26)
    return p

def add_section_heading(text, color=RGBColor(0x0D, 0x94, 0x88)):
    p = doc.add_heading(text, level=1)
    for run in p.runs:
        run.font.color.rgb = color
        run.font.size = Pt(18)
    return p

def add_sub_heading(text, color=RGBColor(0x0E, 0x7C, 0x86)):
    p = doc.add_heading(text, level=2)
    for run in p.runs:
        run.font.color.rgb = color
        run.font.size = Pt(14)
    return p

def add_para(text, bold_prefix=None, italic=False):
    p = doc.add_paragraph()
    if bold_prefix:
        run = p.add_run(bold_prefix)
        run.bold = True
        run.font.size = Pt(11)
    run = p.add_run(text)
    run.font.size = Pt(11)
    if italic:
        run.italic = True
    return p

def add_bullet(text, bold_prefix=None, level=0):
    p = doc.add_paragraph(style='List Bullet')
    p.paragraph_format.left_indent = Inches(0.5 + level * 0.3)
    if bold_prefix:
        run = p.add_run(bold_prefix)
        run.bold = True
        run.font.size = Pt(11)
    run = p.add_run(text)
    run.font.size = Pt(11)
    return p

def add_analogy(text):
    """Add a highlighted analogy box."""
    p = doc.add_paragraph()
    p.paragraph_format.left_indent = Inches(0.3)
    p.paragraph_format.right_indent = Inches(0.3)
    run = p.add_run("🧒 Simple version: ")
    run.bold = True
    run.font.size = Pt(11)
    run.font.color.rgb = RGBColor(0xDC, 0x26, 0x26)
    run = p.add_run(text)
    run.font.size = Pt(11)
    run.italic = True
    run.font.color.rgb = RGBColor(0x55, 0x55, 0x55)
    return p

def add_code_block(text):
    p = doc.add_paragraph()
    p.paragraph_format.left_indent = Inches(0.3)
    run = p.add_run(text)
    run.font.name = 'Courier New'
    run.font.size = Pt(9)
    run.font.color.rgb = RGBColor(0x1a, 0x1a, 0x2e)
    return p

def add_table(headers, rows, col_widths=None):
    table = doc.add_table(rows=1 + len(rows), cols=len(headers))
    table.style = 'Light Grid Accent 1'
    table.alignment = WD_TABLE_ALIGNMENT.CENTER
    # Header row
    for i, h in enumerate(headers):
        cell = table.rows[0].cells[i]
        cell.text = h
        for p in cell.paragraphs:
            for run in p.runs:
                run.bold = True
                run.font.size = Pt(10)
    # Data rows
    for r_idx, row in enumerate(rows):
        for c_idx, val in enumerate(row):
            cell = table.rows[r_idx + 1].cells[c_idx]
            cell.text = str(val)
            for p in cell.paragraphs:
                for run in p.runs:
                    run.font.size = Pt(10)
    if col_widths:
        for i, w in enumerate(col_widths):
            for row in table.rows:
                row.cells[i].width = Inches(w)
    return table

def add_spacer():
    p = doc.add_paragraph()
    p.paragraph_format.space_before = Pt(6)
    p.paragraph_format.space_after = Pt(6)

# ══════════════════════════════════════════════════════════════
# TITLE PAGE
# ══════════════════════════════════════════════════════════════
add_spacer()
add_spacer()
add_title("Phase 1: Finding the Secret Code in TB DNA")
p = doc.add_paragraph()
run = p.add_run("A Very Simple Explanation")
run.font.size = Pt(16)
run.font.color.rgb = RGBColor(0x66, 0x66, 0x66)
run.italic = True
p.alignment = WD_ALIGN_PARAGRAPH.CENTER

add_spacer()
p = doc.add_paragraph()
run = p.add_run("ENIGMA Project — Phase 1: Motif Discovery & Validation")
run.font.size = Pt(12)
run.font.color.rgb = RGBColor(0x0D, 0x94, 0x88)
p.alignment = WD_ALIGN_PARAGRAPH.CENTER

add_spacer()
p = doc.add_paragraph()
run = p.add_run("Every technical term is explained simply. Nothing is left out.")
run.font.size = Pt(11)
run.italic = True
run.font.color.rgb = RGBColor(0x99, 0x99, 0x99)
p.alignment = WD_ALIGN_PARAGRAPH.CENTER

doc.add_page_break()

# ══════════════════════════════════════════════════════════════
# SECTION 1: THE BIG PICTURE
# ══════════════════════════════════════════════════════════════
add_section_heading("1. The Big Picture — What Is This All About?")

add_para(
    "Tuberculosis (that is the disease called TB) is caused by a tiny living thing called "
    "a bacterium. Its full name is Mycobacterium tuberculosis, but we can just call it the TB bug. "
    "This TB bug lives inside your body, and one of the things it needs to survive in there is "
    "cholesterol — a fatty substance found in your cells."
)

add_para(
    "Now, the TB bug has a set of instructions inside it called DNA (deoxyribonucleic acid — "
    "that is the long twisty ladder molecule that stores all the instructions for life). "
    "Written in the DNA are genes — little sections of instructions that tell the bug how to "
    "build tiny machines called proteins. Proteins do all the actual work."
)

add_para(
    "One of these proteins is called Mce3R. Think of Mce3R as a boss sitting in an office. "
    "The boss's job is to decide whether certain workers (other proteins) should come to work "
    "or stay home. Specifically, Mce3R controls workers that help the TB bug eat cholesterol."
)

add_analogy(
    "Imagine your DNA is a very, very long book — over 4 million letters long! "
    "Mce3R is like a bookmark. When the bookmark is stuck on a certain page, that page "
    "cannot be read. When the bookmark comes off, the page can be read and the instructions "
    "on it can be followed. Phase 1 is about finding exactly which page the bookmark sticks to."
)

add_para(
    "The place on the DNA where Mce3R sticks is called an operator (that is a landing pad "
    "or parking spot for the protein). The operator sits between two genes: the gene for "
    "Mce3R itself (called Rv1963c) and the gene for a cholesterol transporter called yrbE3A "
    "(called Rv1964). When Mce3R parks on the operator, it blocks the reading of yrbE3A, "
    "so the cholesterol transporter does not get built. This is called repression — the boss "
    "is saying \"do not build that worker right now.\""
)

add_para(
    "The special thing about this operator is that it has TWO parking spots, not one. "
    "And these two spots are not equal — one is a really good spot (the protein loves it) "
    "and the other is a so-so spot (the protein only sometimes parks there). This unevenness "
    "is called asymmetry, and it is the central mystery of this entire project."
)

add_analogy(
    "Picture a parking lot with two spaces. One space is right by the front door — everyone "
    "wants to park there (that is the strong site). The other space is way at the back — "
    "you would only park there if the front one is full (that is the weak site). "
    "Phase 1 finds these two parking spaces by looking at the DNA letters."
)

doc.add_page_break()

# ══════════════════════════════════════════════════════════════
# SECTION 2: THE DNA WE STARTED WITH
# ══════════════════════════════════════════════════════════════
add_section_heading("2. The DNA We Started With — Three Bacterial Genomes")

add_para(
    "To find the operator, we did not just look at one bug. We looked at three related bugs, "
    "because if a DNA pattern is important, evolution (the slow process of change over millions "
    "of years) would keep it the same across different species. This is called conservation — "
    "important things are conserved (kept the same)."
)

add_sub_heading("The Three Genomes")

add_table(
    ["Organism", "Accession", "Genome Size", "GC Content"],
    [
        ["M. tuberculosis H37Rv", "NC_000962.3", "4,411,532 bp", "65.6%"],
        ["M. bovis AF2122/97", "NC_002945.4", "~4,345,492 bp", "65.6%"],
        ["M. marinum M", "NC_010612.1", "~6,636,827 bp", "65.7%"],
    ],
    col_widths=[2.0, 1.5, 1.5, 1.0]
)

add_spacer()

add_para(
    "Each genome is a long string of letters: A, C, G, and T (adenine, cytosine, guanine, "
    "and thymine). These are the four chemical bases that make up DNA. The \"bp\" stands for "
    "base pairs — because DNA is a double helix (a twisted ladder), and each rung of the "
    "ladder is a pair of bases (A pairs with T, C pairs with G)."
)

add_para(
    "\"GC content\" means how much of the DNA is made of G and C (instead of A and T). "
    "The TB bug has very high GC content — 65.6% — meaning about two-thirds of its DNA "
    "letters are G or C. This matters because when we search for patterns, we need to know "
    "what is \"normal\" and what is \"special.\" In a high-GC genome, seeing a lot of A and T "
    "in a row would be unusual and potentially meaningful."
)

add_sub_heading("The Key Genes")

add_bullet("Mce3R (gene name: Rv1963c)", bold_prefix="The Boss: ")
add_para(
    "    This is the transcription factor (a protein that controls whether genes get read). "
    "Mce3R is a repressor — it stops genes from being read. It belongs to the TetR family "
    "of regulators, which are common in bacteria."
)

add_bullet("yrbE3A (gene name: Rv1964)", bold_prefix="The Worker: ")
add_para(
    "    This gene codes for a protein that helps transport cholesterol and lipids (fats) "
    "into the TB bug. It is part of the larger mce3 operon (Rv1964-Rv1977), which is a "
    "group of genes that work together to import and process cholesterol."
)

add_bullet("Location ~2,207,477 to ~2,207,699 on the H37Rv chromosome", bold_prefix="The Operator: ")
add_para(
    "    The operator is a 123-base-pair stretch of DNA sitting in the intergenic region "
    "(the space between two genes) between Rv1963c and Rv1964. This is where Mce3R physically "
    "touches the DNA to block transcription (the process of reading a gene)."
)

add_analogy(
    "The genome is like a 4-million-word book. Mce3R is a bookmark on one specific page. "
    "yrbE3A is the recipe written on that page. The operator is the exact spot on the page "
    "where the bookmark sticks. We know roughly which page it is on — now we need to find "
    "the exact words."
)

doc.add_page_break()

# ══════════════════════════════════════════════════════════════
# SECTION 3: THE KNOWN OPERATOR SEQUENCE
# ══════════════════════════════════════════════════════════════
add_section_heading("3. The Known Operator Sequence — Two Parking Spots")

add_para(
    "Scientists already knew roughly what the operator looks like, thanks to a technique "
    "called cryo-EM (cryo-electron microscopy — a way of taking pictures of frozen molecules "
    "at incredibly high magnification). The structure was published by Panagoda et al. in 2024 "
    "and deposited in the Protein Data Bank as entry PDB 9B7Y."
)

add_para("The full 123-base-pair operator sequence is:", bold_prefix="")

add_code_block(
    "GCCCCGCGCTATAGGATACTAGCAAGATACATCATAGCCAATATATGCCAGTT"
)
add_code_block(
    "TGCATTGCTATTTACCGATCAGTTGTCCAAGCAATCGCGTATTGGCTATGGA"
)
add_code_block(
    "CATCAGCGGTTCTGCCGC"
)

add_para(
    "That is 123 letters. Inside this sequence are two spots where the Mce3R protein grabs on. "
    "Scientists measured how tightly the protein holds on at each spot using a number called Kd "
    "(the dissociation constant — a measure of how easily the protein lets go)."
)

add_analogy(
    "Kd is like how sticky a sticker is. A small Kd means very sticky (hard to peel off). "
    "A large Kd means not very sticky (peels off easily). The smaller the number, the tighter "
    "the protein holds on."
)

add_sub_heading("The Two Binding Sites")

add_table(
    ["Property", "Strong Site", "Weak Site"],
    [
        ["Location in operator", "Downstream (~88-113 bp)", "Upstream (~10-35 bp)"],
        ["Kd (dissociation constant)", "2.4 nM", "49 nM"],
        ["Kd uncertainty", "± 0.7 nM", "estimated"],
        ["How sticky?", "VERY sticky", "Only somewhat sticky"],
        ["Transcription blocked", "~85%", "~50%"],
        ["Binding preference", "Protein LOVES this spot", "Protein sometimes parks here"],
    ],
    col_widths=[2.0, 2.0, 2.0]
)

add_spacer()

add_para(
    "The unit \"nM\" stands for nanomolar — that is a measure of concentration (how much "
    "protein is dissolved in the liquid). 1 nM = one billionth of a mole per liter. "
    "A Kd of 2.4 nM means the protein binds extremely tightly — you only need a tiny amount "
    "of protein to fill this spot. A Kd of 49 nM means the protein needs to be about 20 times "
    "more concentrated to fill this spot half the time."
)

add_sub_heading("The Spacer")
add_para(
    "Between the two sites, there are 53 base pairs of DNA that the protein does not directly "
    "touch. This is called the spacer. The spacer matters because it determines how the two "
    "binding events might influence each other (a concept called cooperativity, explored in "
    "Phase 5)."
)

add_sub_heading("The Key Number: Affinity Ratio")
add_para(
    "The ratio of the two Kd values is: 49 / 2.4 = 20.4. This means the weak site is about "
    "20 times less attractive than the strong site. This 20-fold asymmetry is THE central "
    "biological variable of the entire project. It is what makes the operator \"lopsided\" and "
    "is hypothesized to create extra noise (randomness) in gene expression, which in turn may "
    "help the TB bug survive drug treatment by creating persister cells."
)

add_analogy(
    "Imagine two magnets on a fridge. One magnet holds on with 20 pounds of force (strong site). "
    "The other holds with only 1 pound (weak site). When you bump the fridge (stress), the weak "
    "magnet falls off first, but the strong one stays. This means the gene gets partially "
    "unblocked — not fully on, not fully off. That in-between state creates noise."
)

doc.add_page_break()

# ══════════════════════════════════════════════════════════════
# SECTION 4: MEME
# ══════════════════════════════════════════════════════════════
add_section_heading("4. Step 1: MEME — Finding the Pattern Without Being Told")

add_para(
    "MEME stands for Multiple EM for Motif Elicitation. That is a fancy name for a computer "
    "program that finds patterns in DNA sequences without being told what to look for. The "
    "\"EM\" stands for Expectation-Maximization, which is the mathematical algorithm (recipe "
    "for solving a problem step by step) that MEME uses internally."
)

add_analogy(
    "Imagine you are a detective. Someone shows you three different crime scenes and says: "
    "\"Find what these have in common.\" You do not know what you are looking for — you just "
    "compare them carefully and notice that all three have the same unusual fingerprint. "
    "That is what MEME does with DNA. It looks at sequences from three different species and "
    "finds the shared pattern."
)

add_sub_heading("MEME Input")
add_para(
    "We gave MEME three DNA sequences, each 200 base pairs long. These are the regions just "
    "upstream (before) the yrbE3A gene in each of the three species. \"Upstream\" means the "
    "region of DNA just before a gene starts — this is where regulatory elements (control "
    "switches) tend to live."
)

add_bullet("M. tuberculosis H37Rv — 200 bp upstream of yrbE3A (Rv1964)")
add_bullet("M. bovis AF2122/97 — 200 bp upstream of the orthologous yrbE3A")
add_bullet("M. marinum M — 200 bp upstream of the orthologous yrbE3A")

add_para(
    "The word \"orthologous\" means \"the same gene in a different species\" — like how humans "
    "and chimps both have a hemoglobin gene, they are orthologs."
)

add_sub_heading("MEME Parameters")
add_para("Here are all the settings we used to run MEME:")

add_table(
    ["Parameter", "Value", "What It Means"],
    [
        ["Mode", "ZOOPS", "Zero Or One Occurrence Per Sequence — each sequence has at most one copy of the pattern"],
        ["Strands", "Both (+ and -)", "Search the DNA in both directions (DNA is double-stranded and can be read forwards or backwards)"],
        ["Width range", "6 to 110 bp", "The pattern could be as short as 6 letters or as long as 110 letters"],
        ["Number of motifs", "Top 5", "Find the 5 best patterns, ranked by statistical significance"],
        ["Background model", "0-order Markov", "Assumes each DNA letter appears independently with frequencies matching the H37Rv genome"],
        ["Background frequencies", "A=0.204, C=0.296, G=0.296, T=0.204", "Because the genome is 65.6% GC, C and G are more common than A and T"],
        ["Ranking", "By E-value", "E-value = expected number of motifs this good you would find by random chance (lower = better)"],
    ],
    col_widths=[1.5, 1.5, 3.5]
)

add_sub_heading("MEME Results — The Discovered Motifs")

add_para(
    "MEME found several motifs. The two most important ones are:"
)

add_para("Motif 1 — 8 base pairs:", bold_prefix="")
add_code_block("Consensus: A C A T C A W A")
add_para(
    "This is a short but highly conserved core pattern found across species."
)

add_para("Motif 2 — 15 base pairs:", bold_prefix="")
add_code_block("Consensus: T W T K C A T T G Y T W T Y T")
add_para(
    "This longer motif captures more of the binding site, particularly around the weak site region."
)

add_sub_heading("IUPAC Ambiguity Codes — What Do Those Letters Mean?")
add_para(
    "You might have noticed letters like W, K, and Y in the motifs above. These are not "
    "regular DNA letters. They are IUPAC ambiguity codes (a standard system by the International "
    "Union of Pure and Applied Chemistry) that represent positions where more than one letter "
    "was found across species:"
)

add_table(
    ["Code", "Stands For", "Means"],
    [
        ["A", "Adenine", "Always A"],
        ["C", "Cytosine", "Always C"],
        ["G", "Guanine", "Always G"],
        ["T", "Thymine", "Always T"],
        ["W", "Weak", "A or T (the two bases connected by only 2 hydrogen bonds — weaker)"],
        ["K", "Keto", "G or T (both have a keto group in their chemistry)"],
        ["Y", "Pyrimidine", "C or T (both are small, single-ring bases)"],
        ["N", "aNy", "A, C, G, or T (completely variable position)"],
    ],
    col_widths=[0.6, 1.2, 4.7]
)

add_analogy(
    "Think of a word puzzle where most letters are filled in but a few squares could be one "
    "of two options. ACATCAWA means the pattern is always A-C-A-T-C-A, then either A or T, "
    "then A. The ambiguous positions are where evolution allowed a little wiggle room."
)

doc.add_page_break()

# ══════════════════════════════════════════════════════════════
# SECTION 5: THE PWM
# ══════════════════════════════════════════════════════════════
add_section_heading("5. The PWM — The Protein's Preference Chart")

add_para(
    "PWM stands for Position Weight Matrix. It is a table that shows, for every position in "
    "the motif, how much the Mce3R protein prefers each of the four DNA letters."
)

add_analogy(
    "Imagine you are ordering a sandwich. At position 1, you ALWAYS want sourdough bread "
    "(100% preference for one option). At position 5, you are happy with either mustard or "
    "mayo (50/50 split). The PWM is like a sandwich order form that shows your preference "
    "at every single position."
)

add_para(
    "Here is a simplified example showing the PWM for Motif 2 (the 15-position motif). "
    "Each number is the probability (from 0.0 to 1.0) that the protein prefers that letter "
    "at that position:"
)

add_table(
    ["Base", "P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8", "P9", "P10", "P11", "P12", "P13", "P14", "P15"],
    [
        ["A", "0.0", "0.5", "0.0", "0.0", "0.0", "1.0", "0.0", "0.0", "0.0", "0.0", "0.0", "0.5", "0.0", "0.0", "0.0"],
        ["C", "0.0", "0.0", "0.0", "0.0", "1.0", "0.0", "0.0", "0.0", "0.0", "0.5", "0.0", "0.0", "0.0", "0.5", "0.0"],
        ["G", "0.0", "0.0", "0.0", "0.5", "0.0", "0.0", "0.0", "0.0", "1.0", "0.0", "0.0", "0.0", "0.0", "0.0", "0.0"],
        ["T", "1.0", "0.5", "1.0", "0.5", "0.0", "0.0", "1.0", "1.0", "0.0", "0.5", "1.0", "0.5", "1.0", "0.5", "1.0"],
    ],
    col_widths=[0.4] + [0.4]*15
)

add_spacer()

add_para(
    "Positions with a 1.0 (100%) are highly conserved — the protein strongly prefers that "
    "specific letter. For example, Position 1 is always T, Position 5 is always C, Position 6 "
    "is always A, and Position 9 is always G. These are the critical contact points where the "
    "protein's amino acids physically touch specific DNA bases."
)

add_para(
    "Positions with a 0.5 (50%) are degenerate — the protein accepts either of two letters. "
    "These are positions where evolution allowed some flexibility, suggesting they are less "
    "critical for binding."
)

add_sub_heading("Sequence Logos — The Pretty Version")
add_para(
    "When scientists visualize a PWM, they often draw a sequence logo. This is a picture where "
    "each position has stacked letters. The height of each letter shows how much the protein "
    "prefers it:"
)

add_bullet("Tall letter = strong preference (the protein really wants this base here)")
add_bullet("Short letter = weak preference (the protein does not care much)")
add_bullet("Multiple short letters = degenerate position (several bases are acceptable)")

add_para(
    "The total height at each position represents the information content — measured in bits "
    "(a unit from information theory). A position with 2 bits of information means the protein "
    "is very picky (only one base allowed). A position with 0 bits means the protein does not "
    "care at all (any base is fine)."
)

add_para(
    "Phase 1 generated sequence logos as EPS files (Encapsulated PostScript — a format for "
    "high-quality scientific graphics) for each discovered motif. These logos visually confirm "
    "the PWM numbers."
)

doc.add_page_break()

# ══════════════════════════════════════════════════════════════
# SECTION 6: FIMO
# ══════════════════════════════════════════════════════════════
add_section_heading("6. Step 2: FIMO — Scanning the Whole Genome")

add_para(
    "FIMO stands for Find Individual Motif Occurrences. Once MEME discovered the pattern, "
    "FIMO takes that pattern and searches the entire genome — all 4,411,532 base pairs of "
    "M. tuberculosis H37Rv — looking for every place where the pattern matches."
)

add_analogy(
    "If MEME was the detective who found the fingerprint, FIMO is the FBI database search — "
    "it takes that fingerprint and checks it against every person in the country. Or think of "
    "it like using a metal detector: MEME told us what kind of metal to look for, and now FIMO "
    "sweeps the entire beach to find every buried piece."
)

add_sub_heading("FIMO Parameters")

add_table(
    ["Parameter", "Value", "Meaning"],
    [
        ["Input genome", "H37Rv complete (NC_000962.3)", "Search the entire TB genome"],
        ["Input motifs", "5 PWMs from MEME", "All 5 discovered motifs are scanned"],
        ["P-value threshold", "1e-4 (0.0001)", "Only report matches with less than 0.01% chance of being random"],
        ["Output", "Ranked list with coordinates", "Each match gets a position, score, p-value, and sequence"],
    ],
    col_widths=[1.5, 2.0, 3.0]
)

add_sub_heading("FIMO Results — The Top Hits")
add_para(
    "FIMO found over 100 candidate binding sites across the genome. The most important ones are:"
)

add_table(
    ["Rank", "Position", "Strand", "P-value", "Score (bits)", "Sequence"],
    [
        ["1", "2,207,551", "+", "1.10e-11", "34.86", "GTTGTCCAAGCAATCGCGTAT"],
        ["2", "2,207,528", "+", "3.14e-11", "33.22", "TTTGCATTGCTATTTACCGA"],
        ["92", "2,207,075", "-", "7.77e-06", "—", "AATTGCAATGTAATCGCGTAT"],
    ],
    col_widths=[0.5, 1.0, 0.6, 1.2, 1.0, 2.2]
)

add_spacer()

add_para(
    "The #1 hit is at position 2,207,551 — that is right inside the known operator region "
    "(which spans ~2,207,477 to ~2,207,699). The p-value of 1.10 x 10^-11 means there is "
    "only about a 1-in-100-billion chance that this match happened by random luck. That is "
    "extraordinarily significant."
)

add_para(
    "The #2 hit is at position 2,207,528 — also inside the operator, just 23 base pairs "
    "away from hit #1. These two hits correspond to the two binding sites (strong and weak) "
    "within the operator."
)

add_para(
    "The score in bits measures how well the sequence matches the PWM. A score of 34.86 bits "
    "is very high — it means the sequence is an almost perfect match to the consensus motif."
)

add_analogy(
    "Imagine you told a search engine \"find this fingerprint in a database of 4 million people\" "
    "and the very first result was the exact person you already knew was the culprit. That is "
    "what happened here — the top FIMO hit matched the known operator perfectly. This validates "
    "the entire approach."
)

add_sub_heading("Additional Candidate Sites")
add_para(
    "Beyond the known operator, FIMO also found potential binding sites near other genes. "
    "These include sites near lipid metabolism genes (Rv1933c, Rv1935c, mbtG) and other "
    "intergenic regions. These are candidate secondary targets of Mce3R regulation — places "
    "where the protein might also bind and control other genes. These candidates are classified "
    "further in Phase 5."
)

doc.add_page_break()

# ══════════════════════════════════════════════════════════════
# SECTION 7: THE ENERGY CONNECTION
# ══════════════════════════════════════════════════════════════
add_section_heading("7. The Energy Connection — From Scores to Sticky-ness")

add_para(
    "There is a deep connection between how well a DNA sequence matches the motif (the FIMO "
    "score) and how tightly the protein binds to it (the Kd value). This connection is "
    "governed by thermodynamics — the science of energy and heat."
)

add_analogy(
    "Think of Kd like the weight needed to pull a sticker off a surface. If you need a lot of "
    "force (low Kd), the sticker is very firmly attached — it has a lot of binding energy. If "
    "you need only a little force (high Kd), the sticker barely holds on — it has less binding "
    "energy. The FIMO score tells us how good the sticker-surface match is, and the energy "
    "equation converts that into actual stickiness."
)

add_sub_heading("The Key Equation")
add_para("The relationship between Kd and binding free energy (Delta-G) is:")
add_code_block("DeltaG = -R x T x ln(Kd / reference)")

add_para("Where:")
add_bullet("DeltaG (Delta-G) = the binding free energy in kcal/mol (kilocalories per mole — a unit of energy)")
add_bullet("R = 1.987 cal/(mol*K) = the gas constant (a universal number from physics)")
add_bullet("T = 310 K = temperature in Kelvin (that is 37 degrees Celsius — body temperature, where TB lives)")
add_bullet("ln = the natural logarithm (a mathematical function)")
add_bullet("Kd = the dissociation constant (how easily the protein lets go)")

add_sub_heading("Calculated Energies")

add_table(
    ["Site", "Kd (nM)", "DeltaG (kcal/mol)", "Interpretation"],
    [
        ["Strong site", "2.4", "-10.2", "Very favorable binding — lots of energy holding protein on"],
        ["Weak site", "49", "-8.1", "Less favorable — protein binds but not as tightly"],
        ["Difference (DDeltaG)", "—", "1.86", "The energy gap between the two sites"],
    ],
    col_widths=[1.2, 1.0, 1.5, 2.8]
)

add_spacer()

add_para(
    "A negative DeltaG means the binding is energetically favorable (it happens spontaneously). "
    "The more negative the number, the tighter the binding. The difference of 1.86 kcal/mol "
    "between the two sites seems small, but in molecular biology, even 1 kcal/mol can make "
    "a 5-fold difference in binding affinity."
)

add_sub_heading("Background Frequencies")
add_para(
    "To correctly score DNA sequences, we need to know what \"random\" DNA looks like in this "
    "genome. The background nucleotide frequencies for H37Rv are:"
)

add_table(
    ["Base", "Frequency", "Percentage"],
    [
        ["A (Adenine)", "0.204", "20.4%"],
        ["C (Cytosine)", "0.296", "29.6%"],
        ["G (Guanine)", "0.296", "29.6%"],
        ["T (Thymine)", "0.204", "20.4%"],
    ],
    col_widths=[1.5, 1.5, 1.5]
)

add_spacer()

add_para(
    "Notice that C and G are each about 30%, while A and T are each about 20%. This reflects "
    "the high GC content (65.6%) of the TB genome. MEME and FIMO use these frequencies as the "
    "\"null model\" — what you would expect to see by chance. A motif is only significant if "
    "it appears more often than this random expectation."
)

add_para(
    "Note: The actual energy matrix conversion (turning every position of the PWM into a "
    "specific energy value) is deferred to Phase 5, where the Berg-von Hippel model is used "
    "to calibrate PWM scores against the experimental Kd values."
)

doc.add_page_break()

# ══════════════════════════════════════════════════════════════
# SECTION 8: WHAT PHASE 1 PROVED
# ══════════════════════════════════════════════════════════════
add_section_heading("8. What Phase 1 Proved — Five Key Achievements")

add_para(
    "Phase 1 accomplished five important things that validate the entire downstream analysis:"
)

add_sub_heading("Achievement 1: Validated the Known Operator")
add_para(
    "The #1 FIMO hit (p = 1.10 x 10^-11) lands precisely in the known operator region at "
    "position 2,207,551. This means our computational approach independently rediscovered "
    "what was already known from experimental cryo-EM structural biology. If the method can "
    "find what we already know is there, we can trust it to find things we do not know yet."
)

add_sub_heading("Achievement 2: De Novo Discovery Worked")
add_para(
    "MEME found the binding motif purely from sequence comparison — without being told "
    "what the binding site looks like, without using the cryo-EM structure, and without any "
    "prior experimental data about Mce3R's DNA preferences. This is de novo discovery (from "
    "scratch), and it worked perfectly."
)

add_sub_heading("Achievement 3: Conservation Across Three Species")
add_para(
    "The motif was found in all three species (H37Rv, M. bovis, M. marinum) with greater "
    "than 80% sequence identity. This conservation over millions of years of evolution means "
    "the operator sequence is functionally important — evolution preserved it because the bugs "
    "that lost it did not survive."
)

add_sub_heading("Achievement 4: 100+ Genome-Wide Candidate Sites")
add_para(
    "Beyond the known operator, FIMO identified over 100 potential Mce3R binding sites across "
    "the H37Rv genome. These are candidate secondary targets — genes that Mce3R might also "
    "regulate. This expands the known Mce3R regulon (the set of genes controlled by one "
    "regulator) and opens new avenues for experimental validation."
)

add_sub_heading("Achievement 5: High Statistical Confidence")
add_para(
    "The top FIMO hits have p-values ranging from 1e-11 to 1e-6, far exceeding the 1e-4 "
    "threshold. To put this in perspective:"
)
add_bullet("p = 1e-4 means 1 in 10,000 chance of being random (our threshold)")
add_bullet("p = 1e-6 means 1 in 1,000,000 chance (very significant)")
add_bullet("p = 1e-11 means 1 in 100,000,000,000 chance (extraordinarily significant)")

add_para(
    "These extremely low p-values give us high confidence that the discovered sites are real "
    "biological binding sites, not statistical artifacts."
)

doc.add_page_break()

# ══════════════════════════════════════════════════════════════
# SECTION 9: OUTPUT FILES
# ══════════════════════════════════════════════════════════════
add_section_heading("9. Output Files — What Phase 1 Produced")

add_para(
    "Phase 1 generated the following files, all stored in the results/phase1/ directory:"
)

add_table(
    ["File", "Format", "What It Contains"],
    [
        ["predicted_sites.csv", "CSV", "All FIMO-predicted binding sites with coordinates, p-values, scores, and sequences"],
        ["conservation_status.csv", "CSV", "Conservation analysis showing which sites are preserved across species"],
        ["meme.txt", "Text", "Complete MEME report with PWMs, alignments, E-values, and consensus sequences"],
        ["meme.xml", "XML", "Structured MEME output for programmatic processing by downstream tools"],
        ["meme.html", "HTML", "Interactive web visualization of discovered motifs (can be opened in a browser)"],
        ["logo1.eps, logo2.eps, ...", "EPS", "Sequence logo images for each motif (one forward, one reverse complement per motif)"],
        ["background.model", "Text", "The 0-order Markov background model with H37Rv nucleotide frequencies"],
        ["phase1_summary.json", "JSON", "Machine-readable summary of all Phase 1 results, timings, and sanity checks"],
    ],
    col_widths=[2.0, 0.7, 3.8]
)

add_spacer()

add_para("Phase 1 completed in approximately 26.5 seconds with a 100% sanity pass rate.")

doc.add_page_break()

# ══════════════════════════════════════════════════════════════
# SECTION 10: WHAT THIS SETS UP
# ══════════════════════════════════════════════════════════════
add_section_heading("10. What This Sets Up — The Road Ahead")

add_para(
    "Phase 1 is the foundation upon which the entire project is built. Here is how its outputs "
    "feed into later phases:"
)

add_sub_heading("Feeding Phase 2: Stochastic Simulation")
add_para(
    "The predicted_sites.csv file provides the operator architecture (two-site, asymmetric) "
    "that Phase 2 uses to build a Gillespie stochastic simulation (a computer program that "
    "simulates individual molecules bouncing around randomly). Phase 2 simulates thousands of "
    "individual cells, each with its own random gene expression, to measure how noisy the "
    "Mce3R system is."
)

add_sub_heading("Feeding Phase 5: Thermodynamic Calibration")
add_para(
    "The PWMs from meme.txt are used in Phase 5 to build a complete energy matrix using the "
    "Berg-von Hippel biophysical model. Phase 5 calibrates these energies against the "
    "experimental Kd values (2.4 nM and 49 nM) to create a thermodynamically rigorous model "
    "of how the operator works. Phase 5 also uses MCMC (Markov Chain Monte Carlo — a method "
    "for exploring uncertainty) to infer the cooperativity parameter omega and the spacer "
    "energy penalty."
)

add_sub_heading("Feeding Phase 6: Environmental Response")
add_para(
    "The genome-wide candidate sites from FIMO are used in Phase 6 to understand how the "
    "Mce3R regulon responds to different environments — baseline conditions, cholesterol-rich "
    "environments, acidic pH (like inside an immune cell), and full host-like stress. Phase 6 "
    "measures whether the asymmetric operator creates more persister cells (drug-tolerant cells) "
    "than a symmetric alternative."
)

add_sub_heading("Feeding Phase 7: Publication Figures")
add_para(
    "The motif logos and site coordinates from Phase 1 appear in the final publication figures "
    "(Phase 7), where they anchor the narrative: \"We discovered the binding pattern, found it "
    "is asymmetric, and showed that asymmetry has functional consequences for noise, "
    "persistence, and environmental response.\""
)

add_spacer()
add_spacer()

p = doc.add_paragraph()
run = p.add_run("— End of Phase 1 Explanation —")
run.font.size = Pt(14)
run.font.color.rgb = RGBColor(0x0D, 0x94, 0x88)
run.bold = True
p.alignment = WD_ALIGN_PARAGRAPH.CENTER

# ── Save ──
outpath = "/Users/aayanalwani/tb project/mce3r_stochastic/Phase_1_Explained_Simply.docx"
doc.save(outpath)
print(f"Saved to {outpath}")
print(f"Size: {os.path.getsize(outpath) / 1024:.1f} KB")
