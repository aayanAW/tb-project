#!/usr/bin/env python3
"""Generate Phase 1 presentation using python-pptx."""
from pptx import Presentation
from pptx.util import Inches, Pt, Emu
from pptx.dml.color import RGBColor
from pptx.enum.text import PP_ALIGN, MSO_ANCHOR
from pptx.enum.shapes import MSO_SHAPE
import os

TEAL = RGBColor(0x0D, 0x94, 0x88)
DARK_TEAL = RGBColor(0x0A, 0x6E, 0x64)
WHITE = RGBColor(0xFF, 0xFF, 0xFF)
DARK = RGBColor(0x1E, 0x29, 0x3B)
GRAY = RGBColor(0x64, 0x74, 0x8B)
LIGHT_BG = RGBColor(0xF8, 0xFA, 0xFC)
ACCENT = RGBColor(0xDC, 0x26, 0x26)

BASE = "/Users/aayanalwani/tb project/mce3r_stochastic"
FIGS = f"{BASE}/results/phase1/slide_figures"

prs = Presentation()
prs.slide_width = Inches(13.333)
prs.slide_height = Inches(7.5)

def add_blank_slide():
    layout = prs.slide_layouts[6]  # blank
    return prs.slides.add_slide(layout)

def add_bg(slide, color=LIGHT_BG):
    bg = slide.background
    fill = bg.fill
    fill.solid()
    fill.fore_color.rgb = color

def add_shape(slide, x, y, w, h, color):
    shape = slide.shapes.add_shape(MSO_SHAPE.RECTANGLE, Inches(x), Inches(y), Inches(w), Inches(h))
    shape.fill.solid()
    shape.fill.fore_color.rgb = color
    shape.line.fill.background()
    return shape

def add_text_box(slide, x, y, w, h, text, size=16, color=DARK, bold=False, align=PP_ALIGN.LEFT, font_name='Calibri'):
    txBox = slide.shapes.add_textbox(Inches(x), Inches(y), Inches(w), Inches(h))
    tf = txBox.text_frame
    tf.word_wrap = True
    p = tf.paragraphs[0]
    p.text = text
    p.font.size = Pt(size)
    p.font.color.rgb = color
    p.font.bold = bold
    p.font.name = font_name
    p.alignment = align
    return txBox

def add_bullets(slide, x, y, w, h, items, size=14, color=DARK, spacing=Pt(8)):
    txBox = slide.shapes.add_textbox(Inches(x), Inches(y), Inches(w), Inches(h))
    tf = txBox.text_frame
    tf.word_wrap = True
    for i, item in enumerate(items):
        if i == 0:
            p = tf.paragraphs[0]
        else:
            p = tf.add_paragraph()
        # Handle bold prefix
        if isinstance(item, tuple):
            run_bold = p.add_run()
            run_bold.text = item[0]
            run_bold.font.bold = True
            run_bold.font.size = Pt(size)
            run_bold.font.color.rgb = color
            run_bold.font.name = 'Calibri'
            run_normal = p.add_run()
            run_normal.text = item[1]
            run_normal.font.size = Pt(size)
            run_normal.font.color.rgb = color
            run_normal.font.name = 'Calibri'
        else:
            p.text = item
            p.font.size = Pt(size)
            p.font.color.rgb = color
            p.font.name = 'Calibri'
        p.space_after = spacing
        p.level = 0
        # Bullet character
        pPr = p._pPr
        if pPr is None:
            from pptx.oxml.ns import qn
            from lxml import etree
            pPr = p._p.get_or_add_pPr()
        from pptx.oxml.ns import qn
        from lxml import etree
        buChar = etree.SubElement(pPr, qn('a:buChar'))
        buChar.set('char', '\u2022')
    return txBox

def add_image(slide, path, x, y, w=None, h=None):
    if w and h:
        slide.shapes.add_picture(path, Inches(x), Inches(y), Inches(w), Inches(h))
    elif w:
        slide.shapes.add_picture(path, Inches(x), Inches(y), width=Inches(w))
    else:
        slide.shapes.add_picture(path, Inches(x), Inches(y))

def teal_header_bar(slide):
    add_shape(slide, 0, 0, 13.333, 0.15, TEAL)

# ══════════════════════════════════════════════════════════════
# SLIDE 1: Title
# ══════════════════════════════════════════════════════════════
slide = add_blank_slide()
add_bg(slide, DARK)
add_shape(slide, 0, 0, 13.333, 0.2, TEAL)
add_shape(slide, 0, 7.3, 13.333, 0.2, TEAL)

add_text_box(slide, 1, 1.8, 11, 1.5, "Phase 1: Binding Motif Discovery",
             size=44, color=WHITE, bold=True, align=PP_ALIGN.CENTER)
add_text_box(slide, 1, 3.3, 11, 0.8, "ENIGMA Project",
             size=28, color=TEAL, bold=True, align=PP_ALIGN.CENTER)
add_text_box(slide, 1, 4.2, 11, 0.8,
             "Finding Where Mce3R Binds DNA in Mycobacterium tuberculosis",
             size=18, color=GRAY, align=PP_ALIGN.CENTER)
add_text_box(slide, 1, 5.5, 11, 0.5, "Panagoda et al. 2024  |  Computational Biology",
             size=12, color=GRAY, align=PP_ALIGN.CENTER)

# ══════════════════════════════════════════════════════════════
# SLIDE 2: The Problem
# ══════════════════════════════════════════════════════════════
slide = add_blank_slide()
add_bg(slide, LIGHT_BG)
teal_header_bar(slide)

add_text_box(slide, 0.8, 0.5, 11, 0.8, "The Problem", size=36, color=DARK, bold=True)

# Big stat callout
add_shape(slide, 0.8, 1.6, 3.5, 1.8, TEAL)
add_text_box(slide, 0.8, 1.7, 3.5, 0.9, "1.23 Million", size=36, color=WHITE, bold=True, align=PP_ALIGN.CENTER)
add_text_box(slide, 0.8, 2.5, 3.5, 0.6, "TB deaths per year", size=16, color=WHITE, align=PP_ALIGN.CENTER)

add_bullets(slide, 5.0, 1.6, 7.5, 4.5, [
    "Mce3R is a transcriptional repressor that controls cholesterol metabolism genes critical for bacterial survival inside the human host",
    "It binds an unusual ASYMMETRIC operator with two binding sites of dramatically different affinity (Kd = 2.4 nM vs 49 nM)",
    "Antibiotic-tolerant persister cells are a major barrier to TB treatment. Deleting mce3R increases persister frequency (Pandey et al. 2023)",
    ("Question: ", "Where exactly does Mce3R bind DNA, and are there other binding sites across the genome?"),
], size=14)

# ══════════════════════════════════════════════════════════════
# SLIDE 3: Pipeline
# ══════════════════════════════════════════════════════════════
slide = add_blank_slide()
add_bg(slide, LIGHT_BG)
teal_header_bar(slide)
add_text_box(slide, 0.8, 0.5, 11, 0.8, "The Pipeline", size=36, color=DARK, bold=True)
add_image(slide, f"{FIGS}/pipeline_flowchart.png", 0.5, 1.5, w=12.3, h=4.5)
add_text_box(slide, 0.8, 6.2, 11, 0.6,
             "3 genomes \u2192 MEME motif discovery \u2192 FIMO genome scan \u2192 1,442 candidate sites",
             size=14, color=GRAY, align=PP_ALIGN.CENTER)

# ══════════════════════════════════════════════════════════════
# SLIDE 4: Three Genomes
# ══════════════════════════════════════════════════════════════
slide = add_blank_slide()
add_bg(slide, LIGHT_BG)
teal_header_bar(slide)
add_text_box(slide, 0.8, 0.5, 11, 0.8, "Three Genomes, One Pattern", size=36, color=DARK, bold=True)

# Three cards
for i, (name, acc, size_bp, note, clr) in enumerate([
    ("M. tuberculosis H37Rv", "NC_000962.3", "4,411,532 bp", "Primary — human TB", TEAL),
    ("M. bovis AF2122/97", "NC_002945.4", "~4,345,492 bp", "Close relative — cattle TB", RGBColor(0x25, 0x63, 0xEB)),
    ("M. marinum M", "NC_010612.1", "~6,636,827 bp", "Distant relative — fish TB", RGBColor(0x7C, 0x3A, 0xED)),
]):
    x = 0.8 + i * 4.0
    add_shape(slide, x, 1.6, 3.6, 2.2, WHITE)
    add_shape(slide, x, 1.6, 3.6, 0.08, clr)
    txBox = add_text_box(slide, x + 0.2, 1.85, 3.2, 0.5, name, size=14, color=DARK, bold=True)
    add_text_box(slide, x + 0.2, 2.35, 3.2, 0.4, acc, size=11, color=GRAY)
    add_text_box(slide, x + 0.2, 2.75, 3.2, 0.4, size_bp, size=13, color=DARK, bold=True)
    add_text_box(slide, x + 0.2, 3.15, 3.2, 0.4, note, size=11, color=GRAY)

add_bullets(slide, 0.8, 4.2, 11.5, 2.5, [
    "If a DNA pattern is conserved across millions of years of evolution, it is functionally important",
    "Extracted 200 bp upstream of yrbE3A (the regulated gene) from each genome",
    "All three genomes have ~65.6% GC content — important for background model",
], size=14)

# ══════════════════════════════════════════════════════════════
# SLIDE 5: GC Content
# ══════════════════════════════════════════════════════════════
slide = add_blank_slide()
add_bg(slide, LIGHT_BG)
teal_header_bar(slide)
add_text_box(slide, 0.8, 0.5, 11, 0.8, "Background: GC Content", size=36, color=DARK, bold=True)
add_image(slide, f"{FIGS}/gc_content.png", 0.5, 1.4, w=5.5, h=4.0)
add_text_box(slide, 6.5, 1.8, 6.0, 3.5,
    "The TB genome is 65.6% GC \u2014 meaning C and G each appear ~30% of the time, "
    "while A and T each appear ~20%.\n\n"
    "MEME uses these frequencies as the \"null model\" so it only reports patterns "
    "that are more structured than random genome composition.\n\n"
    "Without this correction, any stretch of A/T-rich DNA would falsely appear significant "
    "in a GC-rich genome.",
    size=14, color=DARK)

# ══════════════════════════════════════════════════════════════
# SLIDE 6: Operator Architecture
# ══════════════════════════════════════════════════════════════
slide = add_blank_slide()
add_bg(slide, LIGHT_BG)
teal_header_bar(slide)
add_text_box(slide, 0.8, 0.5, 11, 0.8, "The Operator Architecture", size=36, color=DARK, bold=True)
add_image(slide, f"{FIGS}/operator_architecture.png", 0.3, 1.3, w=12.5, h=3.3)

# Key stats
for i, (label, value, clr) in enumerate([
    ("Strong Site", "Kd = 2.4 nM", RGBColor(0x05, 0x96, 0x69)),
    ("Weak Site", "Kd = 49 nM", ACCENT),
    ("Spacer", "53 bp", GRAY),
    ("Affinity Ratio", "20.4\u00d7", RGBColor(0x7C, 0x3A, 0xED)),
]):
    x = 0.8 + i * 3.1
    add_text_box(slide, x, 5.0, 2.8, 0.4, label, size=12, color=GRAY, align=PP_ALIGN.CENTER)
    add_text_box(slide, x, 5.4, 2.8, 0.5, value, size=22, color=clr, bold=True, align=PP_ALIGN.CENTER)

add_text_box(slide, 0.8, 6.3, 11, 0.5,
    "123 bp operator between mce3R (Rv1963c) and yrbE3A (Rv1964) at H37Rv coordinates ~2,207,477\u20132,207,699",
    size=12, color=GRAY, align=PP_ALIGN.CENTER)

# ══════════════════════════════════════════════════════════════
# SLIDE 7: Operator Sequence
# ══════════════════════════════════════════════════════════════
slide = add_blank_slide()
add_bg(slide, DARK)
add_shape(slide, 0, 0, 13.333, 0.15, TEAL)
add_text_box(slide, 0.8, 0.5, 11, 0.8, "The Operator Sequence (123 bp)", size=36, color=WHITE, bold=True)

# Sequence in monospace
seq_lines = [
    "GCCCCGCGCTATAGGATACTAGCAAGATACATCATAGCCAATATATGCCAGTT",
    "TGCATTGCTATTTACCGATCAGTTGTCCAAGCAATCGCGTATTGGCTATGGA",
    "CATCAGCGGTTCTGCCGC"
]
for i, line in enumerate(seq_lines):
    add_text_box(slide, 1.5, 2.0 + i * 0.7, 10, 0.6, line,
                 size=20, color=RGBColor(0x5E, 0xEA, 0xD4), font_name='Courier New', align=PP_ALIGN.CENTER)

add_text_box(slide, 1, 4.5, 11, 0.5, "Source: Panagoda et al. 2024, PDB 9B7Y (cryo-EM structure)",
             size=14, color=GRAY, align=PP_ALIGN.CENTER)

# Color-coded regions
add_text_box(slide, 1, 5.5, 3.5, 0.5, "\u25a0  Positions 10\u201335: Weak site",
             size=13, color=RGBColor(0xFC, 0xA5, 0xA5))
add_text_box(slide, 5, 5.5, 3.5, 0.5, "\u25a0  Positions 35\u201388: Spacer",
             size=13, color=GRAY)
add_text_box(slide, 9, 5.5, 3.5, 0.5, "\u25a0  Positions 88\u2013113: Strong site",
             size=13, color=RGBColor(0x6E, 0xE7, 0xB7))

# ══════════════════════════════════════════════════════════════
# SLIDE 8: MEME
# ══════════════════════════════════════════════════════════════
slide = add_blank_slide()
add_bg(slide, LIGHT_BG)
teal_header_bar(slide)
add_text_box(slide, 0.8, 0.5, 11, 0.8, "Step 1: MEME \u2014 De Novo Motif Discovery", size=32, color=DARK, bold=True)

add_text_box(slide, 0.8, 1.4, 11, 0.5,
    "Multiple EM for Motif Elicitation \u2014 finds DNA patterns without prior knowledge",
    size=14, color=TEAL, bold=True)

add_bullets(slide, 0.8, 2.0, 5.5, 4.5, [
    "Input: 3 orthologous yrbE3A upstream sequences (200 bp each)",
    "Mode: ZOOPS (Zero Or One Occurrence Per Sequence)",
    "Searches both DNA strands (+ and \u2013)",
    "Width range: 6\u2013110 bp",
    "Background: 0-order Markov model using H37Rv GC frequencies",
    "Reports top 5 motifs ranked by E-value (statistical significance)",
], size=13)

# Analogy box
add_shape(slide, 7.0, 2.0, 5.5, 2.5, WHITE)
add_shape(slide, 7.0, 2.0, 5.5, 0.06, TEAL)
add_text_box(slide, 7.2, 2.2, 5.1, 0.4, "The Analogy", size=14, color=TEAL, bold=True)
add_text_box(slide, 7.2, 2.7, 5.1, 1.5,
    "Like a detective examining 3 crime scenes and asking: "
    "\"What do these have in common?\" MEME compares sequences "
    "from 3 species and finds the shared fingerprint \u2014 "
    "without being told what to look for.",
    size=12, color=DARK)

# Parameters table
add_shape(slide, 7.0, 4.8, 5.5, 2.0, WHITE)
add_shape(slide, 7.0, 4.8, 5.5, 0.06, RGBColor(0x25, 0x63, 0xEB))
add_text_box(slide, 7.2, 4.95, 5.1, 0.4, "Background Frequencies", size=13, color=RGBColor(0x25, 0x63, 0xEB), bold=True)
add_text_box(slide, 7.2, 5.4, 5.1, 1.0,
    "A = 20.4%    C = 29.6%\nG = 29.6%    T = 20.4%",
    size=14, color=DARK, font_name='Courier New')

# ══════════════════════════════════════════════════════════════
# SLIDE 9: MEME Results
# ══════════════════════════════════════════════════════════════
slide = add_blank_slide()
add_bg(slide, LIGHT_BG)
teal_header_bar(slide)
add_text_box(slide, 0.8, 0.5, 11, 0.8, "MEME Results \u2014 Discovered Motifs", size=36, color=DARK, bold=True)

# Motif 1 card
add_shape(slide, 0.8, 1.6, 5.5, 2.2, WHITE)
add_shape(slide, 0.8, 1.6, 0.08, 2.2, TEAL)
add_text_box(slide, 1.2, 1.8, 4.8, 0.4, "Motif 1 (8 bp)", size=18, color=TEAL, bold=True)
add_text_box(slide, 1.2, 2.3, 4.8, 0.5, "A C A T C A W A", size=22, color=DARK, font_name='Courier New', bold=True)
add_text_box(slide, 1.2, 2.9, 4.8, 0.5, "Short conserved core pattern", size=13, color=GRAY)

# Motif 2 card
add_shape(slide, 7.0, 1.6, 5.5, 2.2, WHITE)
add_shape(slide, 7.0, 1.6, 0.08, 2.2, RGBColor(0x25, 0x63, 0xEB))
add_text_box(slide, 7.4, 1.8, 4.8, 0.4, "Motif 2 (15 bp)", size=18, color=RGBColor(0x25, 0x63, 0xEB), bold=True)
add_text_box(slide, 7.4, 2.3, 4.8, 0.5, "T W T K C A T T G Y T W T Y T", size=18, color=DARK, font_name='Courier New', bold=True)
add_text_box(slide, 7.4, 2.9, 4.8, 0.5, "Longer motif, captures weak site region", size=13, color=GRAY)

# IUPAC table
add_text_box(slide, 0.8, 4.2, 11, 0.5, "IUPAC Ambiguity Codes", size=18, color=DARK, bold=True)

codes = [
    ["W", "Weak", "A or T", "2 hydrogen bonds"],
    ["K", "Keto", "G or T", "Keto group bases"],
    ["Y", "Pyrimidine", "C or T", "Single-ring bases"],
    ["N", "aNy", "A, C, G, or T", "Completely variable"],
]
for i, (code, name, means, note) in enumerate(codes):
    y_pos = 4.8 + i * 0.5
    bg = WHITE if i % 2 == 0 else LIGHT_BG
    add_shape(slide, 0.8, y_pos, 11.5, 0.45, bg)
    add_text_box(slide, 1.0, y_pos + 0.05, 1.0, 0.35, code, size=16, color=TEAL, bold=True, font_name='Courier New')
    add_text_box(slide, 2.2, y_pos + 0.05, 2.0, 0.35, name, size=13, color=DARK, bold=True)
    add_text_box(slide, 4.5, y_pos + 0.05, 2.5, 0.35, means, size=13, color=DARK)
    add_text_box(slide, 7.5, y_pos + 0.05, 4.0, 0.35, note, size=12, color=GRAY)

# ══════════════════════════════════════════════════════════════
# SLIDE 10: PWM
# ══════════════════════════════════════════════════════════════
slide = add_blank_slide()
add_bg(slide, LIGHT_BG)
teal_header_bar(slide)
add_text_box(slide, 0.8, 0.5, 11, 0.8, "The Position Weight Matrix (PWM)", size=36, color=DARK, bold=True)
add_image(slide, f"{FIGS}/pwm_heatmap.png", 0.5, 1.5, w=8.5, h=2.8)
add_text_box(slide, 9.5, 1.8, 3.3, 3.5,
    "Each cell = probability (0\u20131) that the protein wants that base at that position.\n\n"
    "Dark red (1.0) = absolutely required for binding.\n\n"
    "White (0.0) = never seen at this position.\n\n"
    "This matrix is the \"search template\" that FIMO uses to scan the genome.",
    size=13, color=DARK)

add_text_box(slide, 0.8, 4.8, 11, 1.0,
    "The PWM is derived by aligning the motif instances found across species and counting base frequencies at each position. "
    "It encodes the protein's DNA-binding preferences as a quantitative probability table.",
    size=13, color=GRAY)

# ══════════════════════════════════════════════════════════════
# SLIDE 11: FIMO
# ══════════════════════════════════════════════════════════════
slide = add_blank_slide()
add_bg(slide, LIGHT_BG)
teal_header_bar(slide)
add_text_box(slide, 0.8, 0.5, 11, 0.8, "Step 2: FIMO \u2014 Genome-Wide Scan", size=32, color=DARK, bold=True)
add_text_box(slide, 0.8, 1.3, 11, 0.5,
    "Find Individual Motif Occurrences \u2014 scan all 4,411,532 positions",
    size=14, color=TEAL, bold=True)

add_bullets(slide, 0.8, 1.9, 6.0, 4.0, [
    "Takes the PWM from MEME and slides it across every position in the genome",
    "At each position: computes log-odds score (match vs. random background)",
    "Converts score to p-value (probability of seeing this by chance)",
    "P-value threshold: 1\u00d710\u207b\u2074 (only report if <0.01% chance of random)",
    ("Result: ", "1,442 candidate binding sites across the genome"),
], size=14)

# Analogy box
add_shape(slide, 7.5, 1.9, 5.0, 2.0, WHITE)
add_shape(slide, 7.5, 1.9, 5.0, 0.06, TEAL)
add_text_box(slide, 7.7, 2.1, 4.6, 0.4, "The Analogy", size=14, color=TEAL, bold=True)
add_text_box(slide, 7.7, 2.6, 4.6, 1.0,
    "Like using a metal detector to sweep an entire beach. "
    "MEME told us what metal to look for; FIMO finds every buried piece.",
    size=13, color=DARK)

# Key stat
add_shape(slide, 7.5, 4.3, 5.0, 2.0, TEAL)
add_text_box(slide, 7.5, 4.5, 5.0, 0.8, "1,442", size=48, color=WHITE, bold=True, align=PP_ALIGN.CENTER)
add_text_box(slide, 7.5, 5.3, 5.0, 0.6, "candidate binding sites\nacross the TB genome", size=16, color=WHITE, align=PP_ALIGN.CENTER)

# ══════════════════════════════════════════════════════════════
# SLIDE 12: FIMO Results
# ══════════════════════════════════════════════════════════════
slide = add_blank_slide()
add_bg(slide, LIGHT_BG)
teal_header_bar(slide)
add_text_box(slide, 0.8, 0.5, 11, 0.8, "FIMO Top Hits", size=36, color=DARK, bold=True)
add_image(slide, f"{FIGS}/fimo_results.png", 0.3, 1.3, w=8.5, h=4.3)
add_text_box(slide, 9.2, 1.8, 3.5, 4.0,
    "The #1 and #2 hits land exactly in the known operator region.\n\n"
    "p = 1.1\u00d710\u207b\u00b9\u00b9 means a 1-in-100-billion chance of being random.\n\n"
    "This validates the entire approach: we rediscovered the known binding site without prior bias.\n\n"
    "The tgs1 hit (#21) is biologically exciting \u2014 it links Mce3R to dormancy and lipid body formation.",
    size=13, color=DARK)

# ══════════════════════════════════════════════════════════════
# SLIDE 13: Conservation
# ══════════════════════════════════════════════════════════════
slide = add_blank_slide()
add_bg(slide, LIGHT_BG)
teal_header_bar(slide)
add_text_box(slide, 0.8, 0.5, 11, 0.8, "Conservation Across Species", size=36, color=DARK, bold=True)
add_image(slide, f"{FIGS}/conservation.png", 0.3, 1.4, w=12.5, h=3.0)
add_text_box(slide, 0.8, 4.8, 11, 1.2,
    "Binding sites (red = weak, green = strong) are more conserved than the spacer region. "
    ">80% sequence identity across three species separated by millions of years of evolution. "
    "This confirms the operator is under purifying selection \u2014 organisms that lost this "
    "sequence did not survive.",
    size=14, color=DARK)

# ══════════════════════════════════════════════════════════════
# SLIDE 14: Energy
# ══════════════════════════════════════════════════════════════
slide = add_blank_slide()
add_bg(slide, LIGHT_BG)
teal_header_bar(slide)
add_text_box(slide, 0.8, 0.5, 11, 0.8, "The Energy Connection", size=36, color=DARK, bold=True)
add_image(slide, f"{FIGS}/energy_diagram.png", 0.3, 1.4, w=7.0, h=4.3)

add_text_box(slide, 7.8, 1.6, 5.0, 0.5, "\u0394G = \u2013RT \u00d7 ln(Kd)", size=20, color=DARK, bold=True, font_name='Calibri')
add_text_box(slide, 7.8, 2.3, 5.0, 3.0,
    "PWM scores map directly to physical binding energy.\n\n"
    "R = 1.987 cal/(mol\u00b7K)\n"
    "T = 310 K (37\u00b0C, body temp)\n\n"
    "Strong: \u0394G = \u201310.2 kcal/mol\n"
    "Weak: \u0394G = \u20138.1 kcal/mol\n\n"
    "\u0394\u0394G = 1.86 kcal/mol\n"
    "= 20\u00d7 affinity difference\n\n"
    "This asymmetry is the key variable of the entire project.",
    size=13, color=DARK)

# ══════════════════════════════════════════════════════════════
# SLIDE 15: Key Achievements
# ══════════════════════════════════════════════════════════════
slide = add_blank_slide()
add_bg(slide, LIGHT_BG)
teal_header_bar(slide)
add_text_box(slide, 0.8, 0.5, 11, 0.8, "Key Achievements", size=36, color=DARK, bold=True)

achievements = [
    ("1", "Validated known operator", "Top FIMO hit matches experimental cryo-EM site with p = 1.1\u00d710\u207b\u00b9\u00b9", RGBColor(0x05, 0x96, 0x69)),
    ("2", "De novo discovery", "Found the motif without prior bias, purely from cross-species comparison", RGBColor(0x25, 0x63, 0xEB)),
    ("3", "Conservation confirmed", ">80% identity across M. tuberculosis, M. bovis, M. marinum", RGBColor(0x7C, 0x3A, 0xED)),
    ("4", "1,442 genome-wide sites", "Potential secondary Mce3R targets including tgs1 dormancy gene", RGBColor(0xD9, 0x77, 0x06)),
    ("5", "High confidence", "P-values from 10\u207b\u00b9\u00b9 to 10\u207b\u2076, far exceeding the 10\u207b\u2074 threshold", ACCENT),
]

for i, (num, title, desc, clr) in enumerate(achievements):
    y = 1.5 + i * 1.1
    add_shape(slide, 0.8, y, 0.6, 0.6, clr)
    add_text_box(slide, 0.8, y + 0.05, 0.6, 0.5, num, size=20, color=WHITE, bold=True, align=PP_ALIGN.CENTER)
    add_text_box(slide, 1.7, y + 0.0, 4.0, 0.5, title, size=16, color=DARK, bold=True)
    add_text_box(slide, 1.7, y + 0.4, 10.0, 0.5, desc, size=13, color=GRAY)

# ══════════════════════════════════════════════════════════════
# SLIDE 16: Output Files
# ══════════════════════════════════════════════════════════════
slide = add_blank_slide()
add_bg(slide, LIGHT_BG)
teal_header_bar(slide)
add_text_box(slide, 0.8, 0.5, 11, 0.8, "Output Files", size=36, color=DARK, bold=True)

files = [
    ("predicted_sites.csv", "1,442 FIMO-predicted binding sites with coordinates, p-values, scores, sequences"),
    ("conservation_status.csv", "Cross-species conservation analysis for each predicted site"),
    ("meme.txt / .xml / .html", "Complete MEME reports with PWMs, alignments, E-values, consensus sequences"),
    ("logo1\u20135.eps", "Sequence logos for each discovered motif (forward and reverse complement)"),
    ("background.model", "H37Rv 0-order Markov nucleotide frequency model"),
    ("phase1_summary.json", "Machine-readable summary: timings, counts, sanity check results"),
]

for i, (fname, desc) in enumerate(files):
    y = 1.5 + i * 0.85
    bg = WHITE if i % 2 == 0 else LIGHT_BG
    add_shape(slide, 0.8, y, 11.5, 0.75, bg)
    add_text_box(slide, 1.0, y + 0.1, 3.5, 0.5, fname, size=14, color=TEAL, bold=True, font_name='Courier New')
    add_text_box(slide, 4.8, y + 0.1, 7.2, 0.5, desc, size=13, color=DARK)

add_text_box(slide, 0.8, 6.8, 11, 0.4, "Phase 1 completed in ~26.5 seconds with 100% sanity pass rate",
             size=12, color=GRAY, align=PP_ALIGN.CENTER)

# ══════════════════════════════════════════════════════════════
# SLIDE 17: What Phase 1 Feeds Into
# ══════════════════════════════════════════════════════════════
slide = add_blank_slide()
add_bg(slide, LIGHT_BG)
teal_header_bar(slide)
add_text_box(slide, 0.8, 0.5, 11, 0.8, "What Phase 1 Feeds Into", size=36, color=DARK, bold=True)

phases = [
    ("Phase 2", "Stochastic Simulation", "Operator architecture \u2192 Gillespie SSA with 50,000 cells", TEAL),
    ("Phase 5", "Thermodynamic Calibration", "PWMs \u2192 Berg\u2013von Hippel energy model + MCMC cooperativity inference", RGBColor(0x25, 0x63, 0xEB)),
    ("Phase 6", "Environmental Response", "Candidate sites \u2192 4 stress environments + persistence quantification", RGBColor(0x7C, 0x3A, 0xED)),
    ("Phase 7", "Publication Figures", "Motif logos + site coordinates \u2192 Figures 9\u201314", RGBColor(0xD9, 0x77, 0x06)),
]

for i, (phase, title, desc, clr) in enumerate(phases):
    y = 1.6 + i * 1.3
    add_shape(slide, 0.8, y, 2.0, 1.0, clr)
    add_text_box(slide, 0.8, y + 0.15, 2.0, 0.7, phase, size=20, color=WHITE, bold=True, align=PP_ALIGN.CENTER)
    add_text_box(slide, 3.2, y + 0.0, 4.0, 0.5, title, size=18, color=DARK, bold=True)
    add_text_box(slide, 3.2, y + 0.45, 9.0, 0.5, desc, size=13, color=GRAY)

# ══════════════════════════════════════════════════════════════
# SLIDE 18: Summary
# ══════════════════════════════════════════════════════════════
slide = add_blank_slide()
add_bg(slide, DARK)
add_shape(slide, 0, 0, 13.333, 0.2, TEAL)
add_shape(slide, 0, 7.3, 13.333, 0.2, TEAL)

add_text_box(slide, 1, 1.5, 11, 0.8, "Summary", size=40, color=WHITE, bold=True, align=PP_ALIGN.CENTER)

add_text_box(slide, 1.5, 3.0, 10, 2.5,
    "Phase 1 computationally rediscovered the Mce3R binding site from scratch, "
    "validated it against experimental cryo-EM data, and found 1,442 new candidate "
    "targets across the TB genome \u2014 including a link to the dormancy gene tgs1.",
    size=20, color=WHITE, align=PP_ALIGN.CENTER)

add_text_box(slide, 1, 5.8, 11, 0.5, "ENIGMA: Expression Noise In Gene-regulatory Mechanisms and Architecture",
             size=14, color=TEAL, bold=True, align=PP_ALIGN.CENTER)

# ── Save ──
outpath = f"{BASE}/Phase1_Presentation.pptx"
prs.save(outpath)
print(f"Saved to {outpath}")
print(f"Size: {os.path.getsize(outpath) / 1024:.1f} KB")
print(f"Slides: {len(prs.slides)}")
