#!/usr/bin/env python3
"""Generate Mce3R Pipeline Presentation (.pptx) — 25 slides, dark scientific theme."""

import os
from pptx import Presentation
from pptx.util import Inches, Pt, Emu
from pptx.dml.color import RGBColor
from pptx.enum.text import PP_ALIGN, MSO_ANCHOR
from pptx.enum.shapes import MSO_SHAPE
from pptx.enum.chart import XL_CHART_TYPE, XL_LEGEND_POSITION

# ── Paths ──────────────────────────────────────────────────────────────
PROJECT = os.path.dirname(os.path.abspath(__file__))
FIG_DIR = os.path.join(PROJECT, "results", "figures")
OUT_PATH = os.path.join(PROJECT, "mce3r_presentation.pptx")

# ── Color palette ──────────────────────────────────────────────────────
BG       = RGBColor(0x0B, 0x1D, 0x3A)  # deep navy
BG_LIGHT = RGBColor(0x10, 0x27, 0x4A)  # slightly lighter panel
TEAL     = RGBColor(0x0D, 0x94, 0x88)  # secondary
ACCENT   = RGBColor(0x14, 0xB8, 0xA6)  # accent teal
TEXT     = RGBColor(0xE2, 0xE8, 0xF0)  # light text
GOLD     = RGBColor(0xF5, 0x9E, 0x0B)  # highlight / accent
WHITE    = RGBColor(0xFF, 0xFF, 0xFF)
DIM      = RGBColor(0x94, 0xA3, 0xB8)  # dimmed text
RED_SOFT = RGBColor(0xEF, 0x44, 0x44)
GREEN_S  = RGBColor(0x22, 0xC5, 0x5E)

HEADER_FONT = "Trebuchet MS"
BODY_FONT   = "Calibri"

# Slide dimensions (16:9)
SLIDE_W = Inches(13.333)
SLIDE_H = Inches(7.5)


def set_slide_bg(slide, color=BG):
    bg = slide.background
    fill = bg.fill
    fill.solid()
    fill.fore_color.rgb = color


def add_shape_rect(slide, left, top, width, height, fill_color, border_color=None):
    shape = slide.shapes.add_shape(MSO_SHAPE.ROUNDED_RECTANGLE, left, top, width, height)
    shape.fill.solid()
    shape.fill.fore_color.rgb = fill_color
    if border_color:
        shape.line.color.rgb = border_color
        shape.line.width = Pt(1.5)
    else:
        shape.line.fill.background()
    return shape


def add_textbox(slide, left, top, width, height, text, font_size=18,
                font_name=BODY_FONT, color=TEXT, bold=False, alignment=PP_ALIGN.LEFT):
    txBox = slide.shapes.add_textbox(left, top, width, height)
    tf = txBox.text_frame
    tf.word_wrap = True
    p = tf.paragraphs[0]
    p.text = text
    p.font.size = Pt(font_size)
    p.font.name = font_name
    p.font.color.rgb = color
    p.font.bold = bold
    p.alignment = alignment
    return tf


def add_para(tf, text, font_size=18, color=TEXT, bold=False, font_name=BODY_FONT,
             space_before=Pt(6), alignment=PP_ALIGN.LEFT, level=0):
    p = tf.add_paragraph()
    p.text = text
    p.font.size = Pt(font_size)
    p.font.name = font_name
    p.font.color.rgb = color
    p.font.bold = bold
    p.space_before = space_before
    p.alignment = alignment
    p.level = level
    return p


def add_speaker_notes(slide, text):
    notes_slide = slide.notes_slide
    notes_slide.notes_text_frame.text = text


def try_add_image(slide, filename, left, top, width=None, height=None):
    path = os.path.join(FIG_DIR, filename)
    if os.path.exists(path):
        kwargs = {"image_file": path, "left": left, "top": top}
        if width:
            kwargs["width"] = width
        if height:
            kwargs["height"] = height
        slide.shapes.add_picture(**kwargs)
        return True
    else:
        # placeholder
        shape = add_shape_rect(slide, left, top,
                               width or Inches(5), height or Inches(3),
                               BG_LIGHT, TEAL)
        shape.text_frame.paragraphs[0].text = f"[{filename}]"
        shape.text_frame.paragraphs[0].font.color.rgb = DIM
        shape.text_frame.paragraphs[0].font.size = Pt(14)
        shape.text_frame.paragraphs[0].alignment = PP_ALIGN.CENTER
        return False


def add_accent_line(slide, left, top, width, color=TEAL):
    shape = slide.shapes.add_shape(MSO_SHAPE.RECTANGLE, left, top, width, Pt(3))
    shape.fill.solid()
    shape.fill.fore_color.rgb = color
    shape.line.fill.background()


def add_circle_number(slide, left, top, number, size=Inches(0.55)):
    shape = slide.shapes.add_shape(MSO_SHAPE.OVAL, left, top, size, size)
    shape.fill.solid()
    shape.fill.fore_color.rgb = TEAL
    shape.line.fill.background()
    tf = shape.text_frame
    tf.word_wrap = False
    p = tf.paragraphs[0]
    p.text = str(number)
    p.font.size = Pt(20)
    p.font.color.rgb = WHITE
    p.font.bold = True
    p.alignment = PP_ALIGN.CENTER
    tf.paragraphs[0].space_before = Pt(0)
    shape.text_frame.margin_top = Pt(4)


# ══════════════════════════════════════════════════════════════════════
#  BUILD PRESENTATION
# ══════════════════════════════════════════════════════════════════════
prs = Presentation()
prs.slide_width = SLIDE_W
prs.slide_height = SLIDE_H
blank_layout = prs.slide_layouts[6]  # blank


# ── SLIDE 1: Title ────────────────────────────────────────────────────
sl = prs.slides.add_slide(blank_layout)
set_slide_bg(sl)
add_accent_line(sl, Inches(1), Inches(2.3), Inches(11.3), TEAL)
add_textbox(sl, Inches(1), Inches(2.5), Inches(11.3), Inches(2.0),
            "Genome-Wide Identification of Mce3R\nBinding Sites in Mycobacterium tuberculosis\nUsing a Dual-Motif Computational Pipeline",
            font_size=32, font_name=HEADER_FONT, color=WHITE, bold=True,
            alignment=PP_ALIGN.LEFT)
add_textbox(sl, Inches(1), Inches(4.6), Inches(11.3), Inches(0.6),
            "Informed by the 2024 Asymmetric Operator Structure",
            font_size=20, color=ACCENT, font_name=HEADER_FONT)
add_textbox(sl, Inches(1), Inches(5.5), Inches(5), Inches(0.5),
            "Aayan Alwani", font_size=22, color=DIM, bold=True)
add_textbox(sl, Inches(1), Inches(6.0), Inches(5), Inches(0.5),
            "Great Neck South High School", font_size=16, color=DIM)
# decorative bar bottom
add_accent_line(sl, Inches(1), Inches(6.8), Inches(3), GOLD)
add_speaker_notes(sl, "Title slide. Introduce yourself and the project topic.")


# ── SLIDE 2: The Problem — Tuberculosis ───────────────────────────────
sl = prs.slides.add_slide(blank_layout)
set_slide_bg(sl)
add_textbox(sl, Inches(0.8), Inches(0.4), Inches(11), Inches(0.7),
            "The Problem — Tuberculosis", font_size=30, font_name=HEADER_FONT,
            color=WHITE, bold=True)
add_accent_line(sl, Inches(0.8), Inches(1.1), Inches(4), TEAL)

bullets = [
    ("1.3 million deaths/year", "WHO 2023 — the deadliest bacterial infection"),
    ("Survives inside macrophages", "M. tuberculosis hides inside the very immune cells meant to kill it"),
    ("6–9 month treatment", "Antibiotic-tolerant 'persister' bacteria require prolonged therapy"),
    ("Gene regulation is key", "Transcription factors control which survival genes are switched on/off"),
]
y = Inches(1.5)
for title, sub in bullets:
    add_shape_rect(sl, Inches(0.8), y, Inches(0.12), Inches(0.12), GOLD)
    add_textbox(sl, Inches(1.2), y - Pt(4), Inches(5.2), Inches(0.4),
                title, font_size=20, color=WHITE, bold=True)
    add_textbox(sl, Inches(1.2), y + Inches(0.3), Inches(5.2), Inches(0.4),
                sub, font_size=14, color=DIM)
    y += Inches(1.2)

# Right side — big stat
stat_box = add_shape_rect(sl, Inches(7), Inches(1.5), Inches(5.5), Inches(4.5), BG_LIGHT, TEAL)
tf = stat_box.text_frame
tf.word_wrap = True
p = tf.paragraphs[0]
p.text = "1.3M"
p.font.size = Pt(72)
p.font.color.rgb = GOLD
p.font.bold = True
p.alignment = PP_ALIGN.CENTER
add_para(tf, "deaths per year", font_size=22, color=TEXT, alignment=PP_ALIGN.CENTER)
add_para(tf, "Making TB the #1 bacterial killer worldwide", font_size=14, color=DIM,
         alignment=PP_ALIGN.CENTER, space_before=Pt(20))

add_speaker_notes(sl, "Frame the clinical urgency. TB kills 1.3 million/year. "
    "The bacterium survives inside macrophages and forms persister populations. "
    "Understanding how it regulates survival genes could reveal new drug targets.")


# ── SLIDE 3: What is Mce3R? ──────────────────────────────────────────
sl = prs.slides.add_slide(blank_layout)
set_slide_bg(sl)
add_textbox(sl, Inches(0.8), Inches(0.4), Inches(11), Inches(0.7),
            "What is Mce3R?", font_size=30, font_name=HEADER_FONT, color=WHITE, bold=True)
add_accent_line(sl, Inches(0.8), Inches(1.1), Inches(3), TEAL)

info = [
    "Mce3R (Rv1963c) — a TetR-family transcriptional repressor",
    "Controls the mce3 operon: a cholesterol/lipid import system",
    "When cholesterol is present inside macrophages, Mce3R releases from DNA",
    "This activates lipid import genes needed for intracellular survival",
    "Mce3R is autoregulated — it represses its own gene",
    "Deleting mce3R increases antibiotic persister frequency",
]
y = Inches(1.5)
for line in info:
    add_shape_rect(sl, Inches(0.8), y + Pt(4), Inches(0.1), Inches(0.1), ACCENT)
    add_textbox(sl, Inches(1.15), y, Inches(5.8), Inches(0.45), line, font_size=16, color=TEXT)
    y += Inches(0.55)

# Right side — simplified schematic
box = add_shape_rect(sl, Inches(7.2), Inches(1.5), Inches(5.3), Inches(4.5), BG_LIGHT, TEAL)
tf = box.text_frame
tf.word_wrap = True
p = tf.paragraphs[0]
p.text = "Mce3R Regulatory Model"
p.font.size = Pt(16)
p.font.color.rgb = ACCENT
p.font.bold = True
p.alignment = PP_ALIGN.CENTER

add_para(tf, "", font_size=8, color=BG_LIGHT)
add_para(tf, "─── DNA ─── [ OPERATOR ] ─── mce3 operon ───►",
         font_size=14, color=WHITE, alignment=PP_ALIGN.CENTER)
add_para(tf, "▼", font_size=18, color=GOLD, alignment=PP_ALIGN.CENTER)
add_para(tf, "Mce3R protein sits on operator", font_size=14, color=TEXT, alignment=PP_ALIGN.CENTER)
add_para(tf, "→ BLOCKS transcription of lipid import genes", font_size=14, color=RED_SOFT,
         alignment=PP_ALIGN.CENTER, space_before=Pt(4))
add_para(tf, "", font_size=8, color=BG_LIGHT)
add_para(tf, "+ Cholesterol signal → Mce3R releases",
         font_size=14, color=GREEN_S, alignment=PP_ALIGN.CENTER)
add_para(tf, "→ Lipid genes ACTIVATED", font_size=14, color=GREEN_S, alignment=PP_ALIGN.CENTER)

add_speaker_notes(sl, "Mce3R is a TetR-family repressor that controls cholesterol/lipid import. "
    "It sits on DNA and blocks the mce3 operon. When cholesterol is present inside macrophages, "
    "a ligand causes Mce3R to release, activating lipid import. "
    "Deletion of mce3R increases antibiotic persister frequency — it's clinically relevant.")


# ── SLIDE 4: Asymmetric Operator Discovery ───────────────────────────
sl = prs.slides.add_slide(blank_layout)
set_slide_bg(sl)
add_textbox(sl, Inches(0.8), Inches(0.4), Inches(11), Inches(0.7),
            "The Key Discovery — Asymmetric Operator (2024)",
            font_size=28, font_name=HEADER_FONT, color=WHITE, bold=True)
add_accent_line(sl, Inches(0.8), Inches(1.1), Inches(5), GOLD)

points = [
    ("Classical TetR assumption:", "palindromic (symmetric) DNA operator", TEXT, DIM),
    ("2024 breakthrough:", "Mce3R uses a NONPALINDROMIC, ASYMMETRIC operator", GOLD, WHITE),
    ("Two ~25 bp half-sites", "separated by ~53 bp spacer", TEXT, DIM),
    ("Downstream site:", "higher affinity — Kd = 2.4 ± 0.7 nM", ACCENT, TEXT),
    ("Unprecedented:", "first known asymmetric TetR-family operator", GOLD, DIM),
]
y = Inches(1.5)
for t1, t2, c1, c2 in points:
    add_textbox(sl, Inches(0.8), y, Inches(5.5), Inches(0.35), t1, font_size=16, color=c1, bold=True)
    add_textbox(sl, Inches(0.8), y + Inches(0.3), Inches(5.5), Inches(0.35), t2, font_size=14, color=c2)
    y += Inches(0.75)

# Operator schematic on right
box = add_shape_rect(sl, Inches(6.8), Inches(1.5), Inches(5.8), Inches(3.2), BG_LIGHT, TEAL)
tf = box.text_frame
tf.word_wrap = True
p = tf.paragraphs[0]
p.text = "Validated Operator Architecture"
p.font.size = Pt(14)
p.font.color.rgb = ACCENT
p.font.bold = True
p.alignment = PP_ALIGN.CENTER
add_para(tf, "", font_size=6, color=BG_LIGHT)
add_para(tf, "  ┌─── Upstream ───┐          ┌── Downstream ──┐",
         font_size=13, color=DIM, alignment=PP_ALIGN.CENTER)
add_para(tf, "  │   ~25 bp site  │◄─ 53 bp ─►│   ~25 bp site  │",
         font_size=13, color=WHITE, alignment=PP_ALIGN.CENTER)
add_para(tf, "  └────────────────┘          └────────────────┘",
         font_size=13, color=DIM, alignment=PP_ALIGN.CENTER)
add_para(tf, "       Lower Kd                  Kd = 2.4 nM",
         font_size=12, color=GOLD, alignment=PP_ALIGN.CENTER, space_before=Pt(8))
add_para(tf, "", font_size=6, color=BG_LIGHT)
add_para(tf, "Panagoda, Balázsi & Sampson (2024) ACS Chem. Biol.",
         font_size=11, color=DIM, alignment=PP_ALIGN.CENTER)

add_textbox(sl, Inches(6.8), Inches(5.2), Inches(5.8), Inches(0.5),
            "Citation: Panagoda, Balázsi & Sampson (2024) ACS Chemical Biology",
            font_size=11, color=DIM)

add_speaker_notes(sl, "This is the central scientific motivation. In 2024, Panagoda et al. solved the "
    "cryo-EM structure and showed Mce3R binds an asymmetric nonpalindromic operator — "
    "unprecedented for TetR-family regulators. Two ~25 bp half-sites, 53 bp apart, "
    "with the downstream site having higher affinity (Kd 2.4 nM). "
    "Each Mce3R monomer uses two distinct HTH motifs — Arg53 and Lys262.")


# ── SLIDE 5: Research Question ────────────────────────────────────────
sl = prs.slides.add_slide(blank_layout)
set_slide_bg(sl)
add_textbox(sl, Inches(1), Inches(2.0), Inches(11.3), Inches(1.5),
            "Where else in the M. tuberculosis genome\ndoes Mce3R bind?",
            font_size=38, font_name=HEADER_FONT, color=GOLD, bold=True,
            alignment=PP_ALIGN.CENTER)
add_accent_line(sl, Inches(4), Inches(3.8), Inches(5.3), TEAL)
add_textbox(sl, Inches(2), Inches(4.2), Inches(9.3), Inches(1.0),
            "Can we use the asymmetric operator model to predict\nnew regulatory targets genome-wide?",
            font_size=22, color=TEXT, alignment=PP_ALIGN.CENTER)
add_textbox(sl, Inches(2), Inches(5.5), Inches(9.3), Inches(0.5),
            "Only one operator has been experimentally validated — are there more?",
            font_size=16, color=DIM, alignment=PP_ALIGN.CENTER)

add_speaker_notes(sl, "Frame the research gap. Only one Mce3R operator is experimentally validated. "
    "Our pipeline uses the 2024 asymmetric model to search genome-wide for new binding sites.")


# ── SLIDE 6: Pipeline Overview ────────────────────────────────────────
sl = prs.slides.add_slide(blank_layout)
set_slide_bg(sl)
add_textbox(sl, Inches(0.8), Inches(0.4), Inches(11), Inches(0.7),
            "Pipeline Overview", font_size=30, font_name=HEADER_FONT, color=WHITE, bold=True)
add_accent_line(sl, Inches(0.8), Inches(1.1), Inches(3.5), TEAL)

steps = [
    ("1", "Sequence\nPreparation", "2,351 promoters\nfrom H37Rv genome"),
    ("2", "MEME\nDiscovery", "Find binding motif\nde novo (unbiased)"),
    ("3", "FIMO\nScanning", "Score all promoters\nfor motif matches"),
    ("4", "Analysis &\nRanking", "Architecture, asymmetry\ncomposite score"),
    ("5", "Visualization", "9 publication-quality\nfigures"),
]
x = Inches(0.6)
for num, title, desc in steps:
    add_circle_number(sl, x + Inches(0.65), Inches(1.8), num)
    box = add_shape_rect(sl, x, Inches(2.5), Inches(2.2), Inches(2.8), BG_LIGHT, TEAL)
    tf = box.text_frame
    tf.word_wrap = True
    tf.paragraphs[0].alignment = PP_ALIGN.CENTER
    p = tf.paragraphs[0]
    p.text = title
    p.font.size = Pt(17)
    p.font.color.rgb = ACCENT
    p.font.bold = True
    p.alignment = PP_ALIGN.CENTER
    add_para(tf, desc, font_size=13, color=DIM, alignment=PP_ALIGN.CENTER, space_before=Pt(12))
    # Arrow between boxes (except last)
    if num != "5":
        arrow = sl.shapes.add_shape(MSO_SHAPE.RIGHT_ARROW, x + Inches(2.3), Inches(3.6),
                                     Inches(0.3), Inches(0.3))
        arrow.fill.solid()
        arrow.fill.fore_color.rgb = GOLD
        arrow.line.fill.background()
    x += Inches(2.5)

add_textbox(sl, Inches(0.8), Inches(5.8), Inches(11), Inches(0.5),
            "Python 3.10  •  MEME Suite 5.5.9  •  BioPython  •  pandas  •  matplotlib",
            font_size=14, color=DIM, alignment=PP_ALIGN.CENTER)

add_speaker_notes(sl, "Walk through each pipeline step. Emphasize that MEME discovers the motif "
    "without being told what to look for, then FIMO scores every promoter. "
    "The analysis step combines multiple metrics for confidence ranking.")


# ── SLIDE 7: Sequence Preparation ─────────────────────────────────────
sl = prs.slides.add_slide(blank_layout)
set_slide_bg(sl)
add_textbox(sl, Inches(0.8), Inches(0.4), Inches(11), Inches(0.7),
            "Step 1 — Sequence Preparation", font_size=30, font_name=HEADER_FONT,
            color=WHITE, bold=True)
add_accent_line(sl, Inches(0.8), Inches(1.1), Inches(4), TEAL)

info_items = [
    "Extracted 200–300 bp upstream of every gene in M. tuberculosis H37Rv",
    "2,351 promoter regions total (trimmed at adjacent CDS boundaries)",
    "Also downloaded orthologous sequences from 3 species:",
    "    • M. tuberculosis H37Rv (NC_000962.3)",
    "    • M. bovis AF2122/97 (NC_002945.4)",
    "    • M. marinum M (NC_010612.1)",
    "M. tuberculosis has 65% GC content — unusually high, affects scoring",
]
y = Inches(1.5)
for line in info_items:
    indent = Inches(1.5) if line.startswith("    ") else Inches(0.8)
    add_textbox(sl, indent, y, Inches(6), Inches(0.4), line.strip(), font_size=15, color=TEXT)
    y += Inches(0.45)

# Stats panel
box = add_shape_rect(sl, Inches(7.5), Inches(1.5), Inches(5), Inches(4.5), BG_LIGHT, TEAL)
tf = box.text_frame
tf.word_wrap = True
p = tf.paragraphs[0]
p.text = "Input Data Summary"
p.font.size = Pt(16)
p.font.color.rgb = ACCENT
p.font.bold = True
p.alignment = PP_ALIGN.CENTER
for label, val in [("Genome", "4,411,532 bp"), ("Genes", "~4,000"),
                    ("Promoters extracted", "2,351"), ("Candidate (Mce3R)", "9 sequences"),
                    ("Orthologues", "3 species"), ("GC content", "65.6%")]:
    add_para(tf, f"{label}:", font_size=14, color=DIM, space_before=Pt(10), alignment=PP_ALIGN.LEFT)
    add_para(tf, f"  {val}", font_size=16, color=WHITE, bold=True, space_before=Pt(2), alignment=PP_ALIGN.LEFT)

## Script reference box
scr_box = add_shape_rect(sl, Inches(7.5), Inches(6.1), Inches(5), Inches(1.0), BG_LIGHT, GOLD)
tf_scr = scr_box.text_frame
tf_scr.word_wrap = True
p_scr = tf_scr.paragraphs[0]
p_scr.text = "Scripts Used"
p_scr.font.size = Pt(12)
p_scr.font.color.rgb = GOLD
p_scr.font.bold = True
add_para(tf_scr, "  download_genome.py — fetch H37Rv GenBank from NCBI Entrez", font_size=9, color=TEXT, space_before=Pt(2))
add_para(tf_scr, "  extract_promoters.py — extract 300 bp upstream of each CDS", font_size=9, color=TEXT, space_before=Pt(1))
add_para(tf_scr, "  fetch_sequences.py — download orthologs + tier-specific seqs", font_size=9, color=TEXT, space_before=Pt(1))

add_speaker_notes(sl, "A promoter is the DNA region upstream of a gene where transcription factors bind. "
    "We extract 300 bp upstream of every gene. The 9 candidate promoters are the known "
    "Mce3R-regulated genes used as input for MEME motif discovery. "
    "Cross-species sequences test whether the motif is conserved across mycobacteria. "
    "Scripts: download_genome.py fetches the full H37Rv GenBank record via NCBI Entrez with "
    "retry logic and validates >4 Mb, >3000 CDS features. extract_promoters.py builds a CDS "
    "interval tree, trims promoters at adjacent gene boundaries, and extracts divergent IGRs. "
    "fetch_sequences.py downloads orthologous yrbE3A upstream from M. bovis and M. marinum.")


# ── SLIDE 8: MEME Motif Discovery ────────────────────────────────────
sl = prs.slides.add_slide(blank_layout)
set_slide_bg(sl)
add_textbox(sl, Inches(0.8), Inches(0.4), Inches(11), Inches(0.7),
            "Step 2 — MEME Motif Discovery (Detail)", font_size=30, font_name=HEADER_FONT,
            color=WHITE, bold=True)
add_accent_line(sl, Inches(0.8), Inches(1.1), Inches(5), TEAL)

# Left column — algorithmic detail
bullets = [
    "MEME finds conserved DNA patterns WITHOUT prior knowledge (de novo)",
    "Input: 9 candidate promoters (known Mce3R-regulated genes)",
    "Width range: 20–30 bp  (set by -minw 20 -maxw 30)",
]
y = Inches(1.5)
for b in bullets:
    add_shape_rect(sl, Inches(0.8), y + Pt(5), Inches(0.1), Inches(0.1), ACCENT)
    add_textbox(sl, Inches(1.15), y, Inches(5.5), Inches(0.4), b, font_size=14, color=TEXT)
    y += Inches(0.42)

# EM algorithm box
em_box = add_shape_rect(sl, Inches(0.8), y + Inches(0.1), Inches(5.8), Inches(2.7), BG_LIGHT, TEAL)
tf = em_box.text_frame
tf.word_wrap = True
p = tf.paragraphs[0]
p.text = "Expectation Maximization Algorithm"
p.font.size = Pt(14)
p.font.color.rgb = ACCENT
p.font.bold = True
em_steps = [
    ("E-step:", "Given current PWM, compute probability each position is a motif site"),
    ("M-step:", "Re-estimate PWM from weighted alignment of probable sites"),
    ("Iterate:", "Repeat E→M until convergence (max ΔLL < threshold)"),
    ("Output:", "Position Weight Matrix (PWM) — P(base | position)"),
]
for label, desc in em_steps:
    add_para(tf, f"  {label}  {desc}", font_size=12, color=TEXT, space_before=Pt(6))

# Two modes box
mode_box = add_shape_rect(sl, Inches(0.8), y + Inches(3.0), Inches(5.8), Inches(1.3), BG_LIGHT, GOLD)
tf2 = mode_box.text_frame
tf2.word_wrap = True
p2 = tf2.paragraphs[0]
p2.text = "Two MEME Modes Run Sequentially"
p2.font.size = Pt(13)
p2.font.color.rgb = GOLD
p2.font.bold = True
add_para(tf2, "  ZOOPS (-mod zoops): Zero or One site Per Sequence — standard", font_size=12, color=TEXT, space_before=Pt(6))
add_para(tf2, "  ANR   (-mod anr):  Any Number of Repeats — captures both half-sites", font_size=12, color=TEXT, space_before=Pt(4))

try_add_image(sl, "motif_logo.png", Inches(7), Inches(1.3), width=Inches(5.8))

# Script reference box on right
scr_box8 = add_shape_rect(sl, Inches(7), Inches(4.8), Inches(5.8), Inches(2.2), BG_LIGHT, GOLD)
tf_s8 = scr_box8.text_frame
tf_s8.word_wrap = True
p_s8 = tf_s8.paragraphs[0]
p_s8.text = "Script: run_meme.py"
p_s8.font.size = Pt(13)
p_s8.font.color.rgb = GOLD
p_s8.font.bold = True
add_para(tf_s8, "  build_meme_command() — constructs MEME CLI with -minw/-maxw/-mod flags", font_size=10, color=TEXT, space_before=Pt(3))
add_para(tf_s8, "  run_meme() — executes MEME subprocess, detects availability via shutil.which", font_size=10, color=TEXT, space_before=Pt(2))
add_para(tf_s8, "  simulate_meme_output() — fallback: writes synthetic PWM in MEME v4 format", font_size=10, color=TEXT, space_before=Pt(2))
add_para(tf_s8, "  simulate_asymmetric_meme_output() — creates asymmetric operator PWM", font_size=10, color=TEXT, space_before=Pt(2))
add_para(tf_s8, "  parse_meme_motif() — parses meme.txt → PWM matrix + metadata dict", font_size=10, color=TEXT, space_before=Pt(2))
add_para(tf_s8, "  Output: meme.txt (2 motifs: DOWNSTREAM + UPSTREAM)", font_size=10, color=DIM, space_before=Pt(3))

add_speaker_notes(sl, "MEME uses Expectation Maximization: start with a random alignment, "
    "compute site probabilities (E-step), re-estimate the motif model (M-step), repeat. "
    "We run two modes: ZOOPS assumes each sequence has at most one motif site — "
    "good for finding the single best motif. ANR allows multiple sites per sequence — "
    "this can capture BOTH half-sites of the asymmetric operator in a single run. "
    "Width range 20-30 bp is set based on the 2024 structural data showing ~25 bp half-sites. "
    "Output is a PWM: a 4×W matrix (A/C/G/T × width) of base probabilities.")


# ── SLIDE 9: FIMO Scanning ───────────────────────────────────────────
sl = prs.slides.add_slide(blank_layout)
set_slide_bg(sl)
add_textbox(sl, Inches(0.8), Inches(0.4), Inches(11), Inches(0.7),
            "Step 3 — FIMO Log-Odds Scanning (Detail)", font_size=30, font_name=HEADER_FONT,
            color=WHITE, bold=True)
add_accent_line(sl, Inches(0.8), Inches(1.1), Inches(5), TEAL)

# Left — formula and algorithm
formula_box = add_shape_rect(sl, Inches(0.8), Inches(1.4), Inches(5.8), Inches(2.2), BG_LIGHT, TEAL)
tf = formula_box.text_frame
tf.word_wrap = True
p = tf.paragraphs[0]
p.text = "Log-Odds Scoring Formula"
p.font.size = Pt(14)
p.font.color.rgb = ACCENT
p.font.bold = True
add_para(tf, "", font_size=4, color=BG_LIGHT)
add_para(tf, "  For each window position w of length W:", font_size=12, color=TEXT, space_before=Pt(4))
add_para(tf, "    Score(w) = Σᵢ log₂( PWM[i, base_i] / bg[base_i] )", font_size=13, color=GOLD, space_before=Pt(6))
add_para(tf, "", font_size=4, color=BG_LIGHT)
add_para(tf, "  PWM[i, b] = probability of base b at motif position i", font_size=11, color=DIM, space_before=Pt(2))
add_para(tf, "  bg[b]     = genome background frequency of base b", font_size=11, color=DIM, space_before=Pt(2))
add_para(tf, "  Positive score → sequence more likely under motif than by chance", font_size=11, color=TEXT, space_before=Pt(4))

# Background frequencies box
bg_box = add_shape_rect(sl, Inches(0.8), Inches(3.8), Inches(5.8), Inches(1.3), BG_LIGHT, GOLD)
tf2 = bg_box.text_frame
tf2.word_wrap = True
p2 = tf2.paragraphs[0]
p2.text = "M. tuberculosis Background Frequencies"
p2.font.size = Pt(13)
p2.font.color.rgb = GOLD
p2.font.bold = True
add_para(tf2, "  A: 0.175    C: 0.325    G: 0.325    T: 0.175", font_size=14, color=WHITE, space_before=Pt(8))
add_para(tf2, "  GC = 65% — scoring must account for this bias", font_size=12, color=DIM, space_before=Pt(4))

# P-value box
pv_box = add_shape_rect(sl, Inches(0.8), Inches(5.3), Inches(5.8), Inches(1.5), BG_LIGHT, TEAL)
tf3 = pv_box.text_frame
tf3.word_wrap = True
p3 = tf3.paragraphs[0]
p3.text = "P-Value & Scanning Parameters"
p3.font.size = Pt(13)
p3.font.color.rgb = ACCENT
p3.font.bold = True
add_para(tf3, "  P-value cutoff: < 0.0001 (1 in 10,000)", font_size=12, color=TEXT, space_before=Pt(6))
add_para(tf3, "  Empirical: fraction of random windows scoring ≥ observed", font_size=12, color=TEXT, space_before=Pt(3))
add_para(tf3, "  Both strands scanned (+ and −)", font_size=12, color=TEXT, space_before=Pt(3))
add_para(tf3, "  Two models: palindromic + asymmetric → merged w/ dedup", font_size=12, color=TEXT, space_before=Pt(3))

# Script reference box on right
scr_box9 = add_shape_rect(sl, Inches(7), Inches(1.3), Inches(5.8), Inches(2.5), BG_LIGHT, GOLD)
tf_s9 = scr_box9.text_frame
tf_s9.word_wrap = True
p_s9 = tf_s9.paragraphs[0]
p_s9.text = "Script: run_fimo.py"
p_s9.font.size = Pt(13)
p_s9.font.color.rgb = GOLD
p_s9.font.bold = True
add_para(tf_s9, "  pwm_to_log_odds() — converts PWM probabilities to log₂ scoring matrix", font_size=10, color=TEXT, space_before=Pt(3))
add_para(tf_s9, "  score_window() — scores one DNA window using log-odds matrix", font_size=10, color=TEXT, space_before=Pt(2))
add_para(tf_s9, "  score_sequence_with_pwm() — scans both strands of full sequence", font_size=10, color=TEXT, space_before=Pt(2))
add_para(tf_s9, "  compute_background_score_distribution() — builds empirical null (10K samples)", font_size=10, color=TEXT, space_before=Pt(2))
add_para(tf_s9, "  compute_empirical_pvalue() — p = fraction of null scores ≥ observed", font_size=10, color=TEXT, space_before=Pt(2))
add_para(tf_s9, "  simulate_fimo_scan() — full Python fallback producing canonical fimo.tsv", font_size=10, color=TEXT, space_before=Pt(2))
add_para(tf_s9, "  run_fimo_asymmetric() — separate scan with asymmetric Mce3R PWM", font_size=10, color=TEXT, space_before=Pt(2))

try_add_image(sl, "score_distribution.png", Inches(7), Inches(4.0), width=Inches(5.8))

add_speaker_notes(sl, "FIMO slides a window across every position in all 2,351 promoters. "
    "At each position, it computes a log-odds score: the sum of log2(motif probability / background "
    "probability) for each base. M. tb has 65% GC content, so background frequencies are NOT 25% each — "
    "this is critical to avoid false positives in GC-rich regions. "
    "P-values are computed empirically: sample 100,000 random windows, compute their scores, "
    "then p-value = fraction of random scores ≥ observed score. Threshold is 1 in 10,000. "
    "We run FIMO twice — once with the palindromic PWM, once with the asymmetric PWM — "
    "then merge results, removing duplicates where both models hit the same position.")


# ── SLIDE 10: Palindrome Scoring ─────────────────────────────────────
sl = prs.slides.add_slide(blank_layout)
set_slide_bg(sl)
add_textbox(sl, Inches(0.8), Inches(0.4), Inches(11), Inches(0.7),
            "Step 4a — Palindrome Symmetry Scoring (Detail)", font_size=28, font_name=HEADER_FONT,
            color=WHITE, bold=True)
add_accent_line(sl, Inches(0.8), Inches(1.1), Inches(5), TEAL)

# Formula box
pal_box = add_shape_rect(sl, Inches(0.8), Inches(1.4), Inches(5.8), Inches(2.0), BG_LIGHT, TEAL)
tf = pal_box.text_frame
tf.word_wrap = True
p = tf.paragraphs[0]
p.text = "Palindrome Score Formula"
p.font.size = Pt(14)
p.font.color.rgb = ACCENT
p.font.bold = True
add_para(tf, "", font_size=4, color=BG_LIGHT)
add_para(tf, "  Split matched sequence into left_half and right_half", font_size=12, color=TEXT, space_before=Pt(4))
add_para(tf, "  right_rc = reverse_complement(right_half)", font_size=12, color=TEXT, space_before=Pt(3))
add_para(tf, "  hamming_dist = # mismatches between left_half and right_rc", font_size=12, color=TEXT, space_before=Pt(3))
add_para(tf, "  palindrome_score = 1 − (hamming_dist / half_length)", font_size=13, color=GOLD, space_before=Pt(6))

# Interpretation box
interp_box = add_shape_rect(sl, Inches(0.8), Inches(3.6), Inches(5.8), Inches(2.5), BG_LIGHT, GOLD)
tf2 = interp_box.text_frame
tf2.word_wrap = True
p2 = tf2.paragraphs[0]
p2.text = "Interpretation Ranges"
p2.font.size = Pt(14)
p2.font.color.rgb = GOLD
p2.font.bold = True
ranges = [
    ("Score ≥ 0.85:", "likely_non_mce3r_tetr flag — classical TetR palindrome", RED_SOFT),
    ("Score 0.6–0.85:", "ambiguous — could be either model", DIM),
    ("Score 0.3–0.65:", "EXPECTED Mce3R range — asymmetric operator", GREEN_S),
    ("Score < 0.3:", "fully asymmetric — strong nonpalindromic signature", TEXT),
]
for label, desc, c in ranges:
    add_para(tf2, f"  {label}", font_size=12, color=c, bold=True, space_before=Pt(6))
    add_para(tf2, f"    {desc}", font_size=11, color=DIM, space_before=Pt(1))

# Key insight
add_textbox(sl, Inches(0.8), Inches(6.3), Inches(5.8), Inches(0.5),
            "KEY: Standard TetR pipelines REQUIRE high palindrome scores.\n"
            "Our pipeline PENALIZES them — the opposite approach.",
            font_size=12, color=GOLD, bold=True)

# Script reference box
scr_box10 = add_shape_rect(sl, Inches(7), Inches(1.3), Inches(5.8), Inches(2.2), BG_LIGHT, GOLD)
tf_s10 = scr_box10.text_frame
tf_s10.word_wrap = True
p_s10 = tf_s10.paragraphs[0]
p_s10.text = "Script: analyze_results.py"
p_s10.font.size = Pt(13)
p_s10.font.color.rgb = GOLD
p_s10.font.bold = True
add_para(tf_s10, "  Loads FIMO TSV via pd.read_csv(sep='\\t', comment='#')", font_size=10, color=TEXT, space_before=Pt(3))
add_para(tf_s10, "  Splits matched_sequence at midpoint into halves", font_size=10, color=TEXT, space_before=Pt(2))
add_para(tf_s10, "  Uses utils.reverse_complement() on right half", font_size=10, color=TEXT, space_before=Pt(2))
add_para(tf_s10, "  Hamming distance → palindrome_score column", font_size=10, color=TEXT, space_before=Pt(2))
add_para(tf_s10, "  Flags likely_non_mce3r_tetr if score ≥ 0.85", font_size=10, color=TEXT, space_before=Pt(2))

try_add_image(sl, "palindrome_scores.png", Inches(7), Inches(3.7), width=Inches(5.8))

add_speaker_notes(sl, "The palindrome score splits each matched sequence in half, reverse-complements "
    "the right half, and counts mismatches (Hamming distance). A perfect palindrome scores 1.0 — "
    "both halves are exact reverse complements. The 2024 structural paper showed the real Mce3R "
    "operator scores 0.3-0.65 because the two half-sites are NOT reverse complements. "
    "This is our key innovation: sites scoring ≥ 0.85 are flagged as likely belonging to a "
    "DIFFERENT classical TetR regulator, not Mce3R. Standard TetR pipelines would promote these; "
    "ours penalizes them. The asymmetry_score = 1 - palindrome_score.")


# ── SLIDE 11: Spacing & Architecture ──────────────────────────────────
sl = prs.slides.add_slide(blank_layout)
set_slide_bg(sl)
add_textbox(sl, Inches(0.8), Inches(0.4), Inches(11), Inches(0.7),
            "Step 4b — Site Architecture Classification (Detail)", font_size=28, font_name=HEADER_FONT,
            color=WHITE, bold=True)
add_accent_line(sl, Inches(0.8), Inches(1.1), Inches(5), TEAL)

# Architecture classification table
arch_box = add_shape_rect(sl, Inches(0.8), Inches(1.3), Inches(6), Inches(4.5), BG_LIGHT, TEAL)
tf = arch_box.text_frame
tf.word_wrap = True
p = tf.paragraphs[0]
p.text = "Architecture Classification Rules"
p.font.size = Pt(14)
p.font.color.rgb = ACCENT
p.font.bold = True

arch_types = [
    ("validated_paired_asymmetric", "Same strand, 40–65 bp gap", "HIGHEST", GOLD),
    ("tandem_asymmetric", "Same strand, 0–39 bp gap", "High", GREEN_S),
    ("direct_repeat", "Same strand, >65 bp gap", "Moderate", TEXT),
    ("convergent", "Opposite strands, head-to-head (→ ←)", "Low", DIM),
    ("divergent", "Opposite strands, tail-to-tail (← →)", "Low", DIM),
    ("inverted_repeat", "Classic palindromic TetR geometry", "LOWEST", RED_SOFT),
    ("isolated", "No nearby partner site detected", "Low", DIM),
]
for name, rule, priority, c in arch_types:
    add_para(tf, f"  {name}", font_size=11, color=c, bold=True, space_before=Pt(5))
    add_para(tf, f"    {rule}  →  Priority: {priority}", font_size=10, color=DIM, space_before=Pt(1))

# Spacing detail
sp_box = add_shape_rect(sl, Inches(0.8), Inches(5.9), Inches(6), Inches(1.0), BG_LIGHT, GOLD)
tf2 = sp_box.text_frame
tf2.word_wrap = True
p2 = tf2.paragraphs[0]
p2.text = "Validated Spacing: 53 bp (Panagoda 2024)"
p2.font.size = Pt(13)
p2.font.color.rgb = GOLD
p2.font.bold = True
add_para(tf2, "  matches_validated_spacing flag: 48–58 bp (±5 bp tolerance)", font_size=11, color=TEXT, space_before=Pt(4))
add_para(tf2, "  Spacing = gap between stop of site 1 and start of site 2", font_size=11, color=DIM, space_before=Pt(2))

# Script reference box
scr_box11 = add_shape_rect(sl, Inches(7), Inches(1.3), Inches(5.8), Inches(2.5), BG_LIGHT, GOLD)
tf_s11 = scr_box11.text_frame
tf_s11.word_wrap = True
p_s11 = tf_s11.paragraphs[0]
p_s11.text = "Script: motif_models.py"
p_s11.font.size = Pt(13)
p_s11.font.color.rgb = GOLD
p_s11.font.bold = True
add_para(tf_s11, "  classify_site_pair_architecture() — determines strand/spacing relationship", font_size=10, color=TEXT, space_before=Pt(3))
add_para(tf_s11, "  classify_site_architectures() — identifies all paired sites within max_gap", font_size=10, color=TEXT, space_before=Pt(2))
add_para(tf_s11, "  compare_models() — palindromic vs asymmetric model fit per site", font_size=10, color=TEXT, space_before=Pt(2))
add_para(tf_s11, "  score_sequence_against_consensus() — consensus match scoring", font_size=10, color=TEXT, space_before=Pt(2))
add_para(tf_s11, "  generate_architecture_comparison_report() — compares to validated yrbE3A", font_size=10, color=TEXT, space_before=Pt(2))

try_add_image(sl, "spacing_histogram.png", Inches(7), Inches(4.0), width=Inches(2.8), height=Inches(2.5))
try_add_image(sl, "site_architecture.png", Inches(9.9), Inches(4.0), width=Inches(2.8), height=Inches(2.5))

add_speaker_notes(sl, "Architecture classification is the most biologically informative metric. "
    "The 2024 paper showed Mce3R uses two same-strand half-sites separated by ~53 bp — "
    "this 'validated_paired_asymmetric' geometry gets the highest priority. We define it as "
    "same-strand sites with 40-65 bp spacing (allowing ±12 bp tolerance around 53 bp). "
    "The matches_validated_spacing flag is stricter: 48-58 bp (±5 bp). "
    "Inverted repeats (classic palindromic TetR) get the LOWEST score because the 2024 "
    "paper explicitly showed Mce3R does NOT use this geometry. "
    "Isolated sites have no detectable partner within 80 bp (max_gap parameter).")


# ── SLIDE 12: Composite Scoring ───────────────────────────────────────
sl = prs.slides.add_slide(blank_layout)
set_slide_bg(sl)
add_textbox(sl, Inches(0.8), Inches(0.4), Inches(11), Inches(0.7),
            "Step 4c — Composite Score Formula (Detail)", font_size=28, font_name=HEADER_FONT,
            color=WHITE, bold=True)
add_accent_line(sl, Inches(0.8), Inches(1.1), Inches(5), TEAL)

# Main formula
fbox = add_shape_rect(sl, Inches(0.8), Inches(1.3), Inches(11.5), Inches(0.8), BG_LIGHT, GOLD)
tf = fbox.text_frame
tf.word_wrap = True
p = tf.paragraphs[0]
p.text = "Score = 0.40 × fimo_norm + 0.25 × arch_bonus + 0.20 × asym_weight + 0.15 × tier_prior"
p.font.size = Pt(16)
p.font.color.rgb = GOLD
p.font.bold = True
p.alignment = PP_ALIGN.CENTER

# Four component detail boxes (2×2 grid)
# Component 1: FIMO
c1 = add_shape_rect(sl, Inches(0.8), Inches(2.3), Inches(5.8), Inches(1.5), BG_LIGHT, TEAL)
tf1 = c1.text_frame
tf1.word_wrap = True
p1 = tf1.paragraphs[0]
p1.text = "40% — fimo_norm (Motif Match Quality)"
p1.font.size = Pt(13)
p1.font.color.rgb = GOLD
p1.font.bold = True
add_para(tf1, "  fimo_norm = (score − min) / (max − min)", font_size=11, color=TEXT, space_before=Pt(4))
add_para(tf1, "  Normalized to [0, 1] across all hits", font_size=11, color=DIM, space_before=Pt(2))
add_para(tf1, "  Based on raw FIMO log-odds from PWM scan", font_size=11, color=DIM, space_before=Pt(2))

# Component 2: Architecture
c2 = add_shape_rect(sl, Inches(6.8), Inches(2.3), Inches(5.5), Inches(1.5), BG_LIGHT, TEAL)
tf2 = c2.text_frame
tf2.word_wrap = True
p2 = tf2.paragraphs[0]
p2.text = "25% — arch_bonus (Operator Geometry)"
p2.font.size = Pt(13)
p2.font.color.rgb = ACCENT
p2.font.bold = True
add_para(tf2, "  validated_paired_asymmetric: 1.0", font_size=11, color=WHITE, space_before=Pt(4))
add_para(tf2, "  tandem_asymmetric: 0.6  |  direct_repeat: 0.3", font_size=11, color=DIM, space_before=Pt(2))
add_para(tf2, "  convergent/divergent: 0.15  |  inverted_repeat: 0.05", font_size=11, color=DIM, space_before=Pt(2))

# Component 3: Asymmetry weight
c3 = add_shape_rect(sl, Inches(0.8), Inches(4.0), Inches(5.8), Inches(1.7), BG_LIGHT, TEAL)
tf3 = c3.text_frame
tf3.word_wrap = True
p3 = tf3.paragraphs[0]
p3.text = "20% — asym_weight (Triangular Function)"
p3.font.size = Pt(13)
p3.font.color.rgb = GREEN_S
p3.font.bold = True
add_para(tf3, "  if 0.3 ≤ pal ≤ 0.65:  asym_weight = 1.0  (peak)", font_size=11, color=WHITE, space_before=Pt(4))
add_para(tf3, "  if pal < 0.3:  asym_weight = pal / 0.3  (linear ramp up)", font_size=11, color=DIM, space_before=Pt(2))
add_para(tf3, "  if 0.65 < pal < 0.85:  asym_weight = (0.85−pal) / 0.20  (decay)", font_size=11, color=DIM, space_before=Pt(2))
add_para(tf3, "  if pal ≥ 0.85:  asym_weight = 0.0  (penalized — wrong TetR)", font_size=11, color=RED_SOFT, space_before=Pt(2))

# Component 4: Tier prior
c4 = add_shape_rect(sl, Inches(6.8), Inches(4.0), Inches(5.5), Inches(1.7), BG_LIGHT, TEAL)
tf4 = c4.text_frame
tf4.word_wrap = True
p4 = tf4.paragraphs[0]
p4.text = "15% — tier_prior (Biological Knowledge)"
p4.font.size = Pt(13)
p4.font.color.rgb = DIM
p4.font.bold = True
add_para(tf4, "  Tier 1 (known operator):  prior = 1.0", font_size=11, color=GOLD, space_before=Pt(4))
add_para(tf4, "  Tier 2 (lipid/cholesterol genes):  prior = 0.6", font_size=11, color=GREEN_S, space_before=Pt(2))
add_para(tf4, "  Tier 3 (all other genes):  prior = 0.2", font_size=11, color=DIM, space_before=Pt(2))
add_para(tf4, "  Assigned BEFORE composite score computation", font_size=11, color=TEXT, space_before=Pt(4))

# Script reference box
scr_box12 = add_shape_rect(sl, Inches(0.8), Inches(5.85), Inches(5.8), Inches(0.9), BG_LIGHT, GOLD)
tf_s12 = scr_box12.text_frame
tf_s12.word_wrap = True
p_s12 = tf_s12.paragraphs[0]
p_s12.text = "Script: motif_models.py → compute_composite_score()"
p_s12.font.size = Pt(11)
p_s12.font.color.rgb = GOLD
p_s12.font.bold = True
add_para(tf_s12, "  apply_model_comparison() adds all columns, then compute_composite_score()", font_size=9, color=TEXT, space_before=Pt(2))
add_para(tf_s12, "  analyze_results.py orchestrates: load → palindrome → architecture → composite → rank", font_size=9, color=TEXT, space_before=Pt(1))

# Tier counts
tbox = add_shape_rect(sl, Inches(0.8), Inches(5.9), Inches(11.5), Inches(0.9), BG_LIGHT, TEAL)
tf5 = tbox.text_frame
tf5.word_wrap = True
p5 = tf5.paragraphs[0]
p5.text = "Result:   Tier 1: 33 sites (known)   •   Tier 2: 39 sites (plausible)   •   Tier 3: 631 sites (exploratory)"
p5.font.size = Pt(14)
p5.font.color.rgb = WHITE
p5.alignment = PP_ALIGN.CENTER
add_para(tf5, "Total: 703 ranked predictions — top-ranked = the known validated operator (composite = 0.942)", font_size=12, color=GOLD, alignment=PP_ALIGN.CENTER, space_before=Pt(4))

add_speaker_notes(sl, "The composite score combines four independent signals into a single ranking. "
    "40% comes from raw motif match quality (FIMO log-odds, normalized 0-1). "
    "25% comes from operator geometry — validated_paired_asymmetric (same strand, 40-65 bp) scores 1.0, "
    "while inverted_repeat (classic palindrome) scores only 0.05. "
    "20% is the asymmetry weight — a triangular function that peaks for palindrome scores 0.3-0.65 "
    "(the Mce3R-expected range) and drops to zero for scores ≥ 0.85 (classical TetR). "
    "15% is a biological prior from published literature: Tier 1 (known validated targets) gets 1.0, "
    "Tier 2 (lipid metabolism genes) gets 0.6, Tier 3 (everything else) gets 0.2. "
    "Critically, tier assignment happens BEFORE composite scoring so the prior can inform the ranking.")


# ── SLIDE 13: Visualization Overview ──────────────────────────────────
sl = prs.slides.add_slide(blank_layout)
set_slide_bg(sl)
add_textbox(sl, Inches(0.8), Inches(0.4), Inches(11), Inches(0.7),
            "Step 5 — Visualization (9 Figures)", font_size=30, font_name=HEADER_FONT,
            color=WHITE, bold=True)
add_accent_line(sl, Inches(0.8), Inches(1.1), Inches(4), TEAL)

figs = [
    "score_distribution.png", "motif_position_distribution.png", "spacing_histogram.png",
    "palindrome_scores.png", "gc_vs_hits.png", "motif_logo.png",
    "model_comparison.png", "site_architecture.png", "operator_architecture_diagram.png",
]
labels = [
    "Score Distribution", "Motif Positions", "Spacing Histogram",
    "Palindrome Scores", "GC vs Hits", "Motif Logo",
    "Model Comparison", "Site Architecture", "Operator Diagram",
]
for i, (fig, label) in enumerate(zip(figs, labels)):
    col = i % 3
    row = i // 3
    x = Inches(0.5) + col * Inches(4.2)
    y = Inches(1.5) + row * Inches(1.9)
    try_add_image(sl, fig, x, y, width=Inches(3.8), height=Inches(1.5))
    add_textbox(sl, x, y + Inches(1.5), Inches(3.8), Inches(0.3),
                label, font_size=11, color=DIM, alignment=PP_ALIGN.CENTER)

add_speaker_notes(sl, "Quick overview of all 9 figures generated by the pipeline. "
    "Each provides a different view of the data. We'll examine the key ones in detail.")


# ── SLIDE 14: Validation Results ──────────────────────────────────────
sl = prs.slides.add_slide(blank_layout)
set_slide_bg(sl)
add_textbox(sl, Inches(0.8), Inches(0.4), Inches(11), Inches(0.7),
            "Results — Pipeline Validation", font_size=30, font_name=HEADER_FONT,
            color=WHITE, bold=True)
add_accent_line(sl, Inches(0.8), Inches(1.1), Inches(4), GOLD)

# Big callout
cbox = add_shape_rect(sl, Inches(0.8), Inches(1.5), Inches(5.5), Inches(1.5), BG_LIGHT, GOLD)
tf = cbox.text_frame
tf.word_wrap = True
p = tf.paragraphs[0]
p.text = "Top Hit: IGR_Rv1963c_Rv1964"
p.font.size = Pt(20)
p.font.color.rgb = GOLD
p.font.bold = True
p.alignment = PP_ALIGN.CENTER
add_para(tf, "Composite Score: 0.942", font_size=28, color=WHITE, bold=True,
         alignment=PP_ALIGN.CENTER, space_before=Pt(8))

stats = [
    ("690/703 sites (98%)", "prefer the asymmetric model"),
    ("0 sites scored ≥ 0.8", "palindrome — confirming non-palindromic operator"),
    ("24 validated_paired_asymmetric", "sites match the 53 bp geometry"),
    ("#1 ranked = known operator", "pipeline correctly identifies the validated site"),
]
y = Inches(3.3)
for val, desc in stats:
    add_textbox(sl, Inches(0.8), y, Inches(3), Inches(0.35), val, font_size=16, color=GOLD, bold=True)
    add_textbox(sl, Inches(3.8), y, Inches(3), Inches(0.35), desc, font_size=14, color=TEXT)
    y += Inches(0.55)

try_add_image(sl, "model_comparison.png", Inches(7), Inches(1.3), width=Inches(5.8))

add_speaker_notes(sl, "The pipeline's top-ranked hit IS the known validated operator — this is our "
    "positive control passing. 98% of hits prefer the asymmetric model, and zero sites "
    "have palindrome scores above 0.8, confirming the 2024 structural finding. "
    "The 24 validated_paired_asymmetric sites are our highest-confidence novel predictions.")


# ── SLIDE 15: Simulation vs Real MEME ─────────────────────────────────
sl = prs.slides.add_slide(blank_layout)
set_slide_bg(sl)
add_textbox(sl, Inches(0.8), Inches(0.4), Inches(11), Inches(0.7),
            "Results — Simulation vs. Real MEME", font_size=30, font_name=HEADER_FONT,
            color=WHITE, bold=True)
add_accent_line(sl, Inches(0.8), Inches(1.1), Inches(4.5), TEAL)

# Comparison table
headers = ["Metric", "Simulation", "Real MEME"]
rows = [
    ("Total hits", "152", "703"),
    ("Tier 1 sites", "4", "33"),
    ("Tier 2 sites", "5", "39"),
    ("Top composite score", "0.365", "0.942"),
    ("Validated paired sites", "0", "24"),
]
# Table header
y_t = Inches(1.5)
for j, h in enumerate(headers):
    x_t = Inches(1) + j * Inches(3.5)
    add_textbox(sl, x_t, y_t, Inches(3.3), Inches(0.45), h,
                font_size=16, color=ACCENT, bold=True, alignment=PP_ALIGN.CENTER)
add_accent_line(sl, Inches(1), y_t + Inches(0.45), Inches(10.5), TEAL)

for i, (metric, sim, real) in enumerate(rows):
    y_r = Inches(2.1) + i * Inches(0.6)
    bg_c = BG_LIGHT if i % 2 == 0 else BG
    add_shape_rect(sl, Inches(1), y_r, Inches(10.5), Inches(0.55), bg_c)
    add_textbox(sl, Inches(1), y_r + Pt(4), Inches(3.3), Inches(0.45), metric,
                font_size=15, color=TEXT, alignment=PP_ALIGN.CENTER)
    add_textbox(sl, Inches(4.5), y_r + Pt(4), Inches(3.3), Inches(0.45), sim,
                font_size=15, color=DIM, alignment=PP_ALIGN.CENTER)
    add_textbox(sl, Inches(8), y_r + Pt(4), Inches(3.3), Inches(0.45), real,
                font_size=15, color=GOLD, bold=True, alignment=PP_ALIGN.CENTER)

# Callout
cbox = add_shape_rect(sl, Inches(2), Inches(5.5), Inches(9), Inches(1.2), BG_LIGHT, GOLD)
tf = cbox.text_frame
tf.word_wrap = True
p = tf.paragraphs[0]
p.text = "Key Insight: Real motif discovery dramatically outperforms hardcoded consensus"
p.font.size = Pt(18)
p.font.color.rgb = GOLD
p.font.bold = True
p.alignment = PP_ALIGN.CENTER
add_para(tf, "De novo MEME discovery captures sequence features that a static consensus misses",
         font_size=14, color=DIM, alignment=PP_ALIGN.CENTER, space_before=Pt(8))

add_speaker_notes(sl, "Compare simulation (hardcoded PWM) vs real MEME discovery. "
    "Real MEME found 4.6x more hits, 8x more Tier 1 sites, and the composite score "
    "jumped from 0.365 to 0.942. Most importantly, real MEME found 24 validated paired "
    "asymmetric sites — the simulation found zero. This justifies de novo motif discovery.")


# ── SLIDE 16: 24 High-Confidence Targets ──────────────────────────────
sl = prs.slides.add_slide(blank_layout)
set_slide_bg(sl)
add_textbox(sl, Inches(0.8), Inches(0.4), Inches(11), Inches(0.7),
            "The 24 High-Confidence Targets", font_size=30, font_name=HEADER_FONT,
            color=WHITE, bold=True)
add_accent_line(sl, Inches(0.8), Inches(1.1), Inches(4), GOLD)

add_textbox(sl, Inches(0.8), Inches(1.3), Inches(6), Inches(0.5),
            "Sites with validated_paired_asymmetric architecture (40–65 bp spacer)",
            font_size=14, color=DIM)

targets = [
    ("1", "IGR_Rv1963c_Rv1964", "Known operator (validation)", "0.942", GOLD),
    ("2", "Rv1964 / yrbE3A", "mce3 operon integral membrane", "0.895", GOLD),
    ("3", "Rv1176c", "Hypothetical (54 bp spacer)", "0.857", ACCENT),
    ("4", "Rv3545c / cyp125", "Steroid C26-monooxygenase", "0.815", ACCENT),
    ("5", "Rv3575c / Rv3576", "LacI regulator / lppH (53 bp spacer)", "0.755", ACCENT),
]
y = Inches(1.9)
for rank, locus, desc, score, c in targets:
    add_shape_rect(sl, Inches(0.8), y, Inches(0.4), Inches(0.4), c)
    add_textbox(sl, Inches(0.85), y + Pt(2), Inches(0.35), Inches(0.35), rank,
                font_size=14, color=BG, bold=True, alignment=PP_ALIGN.CENTER)
    add_textbox(sl, Inches(1.4), y, Inches(2.5), Inches(0.4), locus, font_size=15, color=WHITE, bold=True)
    add_textbox(sl, Inches(4), y, Inches(2.8), Inches(0.4), desc, font_size=13, color=DIM)
    add_textbox(sl, Inches(5.8), y + Pt(1), Inches(1), Inches(0.35), score,
                font_size=14, color=c, bold=True, alignment=PP_ALIGN.RIGHT)
    y += Inches(0.55)

try_add_image(sl, "operator_architecture_diagram.png", Inches(7), Inches(1.3), width=Inches(5.8))

add_speaker_notes(sl, "These 24 sites are our strongest predictions — they have the exact same "
    "geometry as the validated operator: two same-strand sites with 40-65 bp spacing. "
    "The top hit IS the known operator (validation). Novel targets include Rv1176c (54 bp spacer), "
    "cyp125 (cholesterol catabolism), and Rv3575c (53 bp spacer — exact match to reference).")


# ── SLIDE 17: Biological Coherence ────────────────────────────────────
sl = prs.slides.add_slide(blank_layout)
set_slide_bg(sl)
add_textbox(sl, Inches(0.8), Inches(0.4), Inches(11), Inches(0.7),
            "Biological Coherence", font_size=30, font_name=HEADER_FONT, color=WHITE, bold=True)
add_accent_line(sl, Inches(0.8), Inches(1.1), Inches(3.5), TEAL)

add_textbox(sl, Inches(0.8), Inches(1.5), Inches(6), Inches(0.5),
            "Tier 2 targets cluster in cholesterol/lipid metabolism — not random",
            font_size=16, color=GOLD, bold=True)

pathway = [
    ("Cholesterol import", "mce3 operon (Rv1964–Rv1977)", ACCENT),
    ("        ↓", "", DIM),
    ("Side-chain oxidation", "cyp125 (Rv3545c) — steroid C26-monooxygenase", TEXT),
    ("        ↓", "", DIM),
    ("Ring degradation", "hsaE, kstD (Rv3536c–Rv3537) — steroid catabolism", TEXT),
    ("        ↓", "", DIM),
    ("Beta-oxidation", "fadE22, fadA5, fadD11 — fatty acid processing", TEXT),
    ("        ↓", "", DIM),
    ("Carbon metabolism", "icl1 (Rv3053c) — glyoxylate shunt", TEXT),
]
y = Inches(2.2)
for step, genes, c in pathway:
    if "↓" in step:
        add_textbox(sl, Inches(1.5), y, Inches(1), Inches(0.3), step, font_size=14, color=DIM)
        y += Inches(0.25)
    else:
        add_textbox(sl, Inches(0.8), y, Inches(2.5), Inches(0.4), step, font_size=15, color=c, bold=True)
        add_textbox(sl, Inches(3.3), y, Inches(3.5), Inches(0.4), genes, font_size=13, color=DIM)
        y += Inches(0.35)

try_add_image(sl, "gc_vs_hits.png", Inches(7), Inches(1.3), width=Inches(5.8))

add_textbox(sl, Inches(7), Inches(5.5), Inches(5.5), Inches(0.5),
            "GC content plot confirms hits are NOT artifacts of GC bias",
            font_size=12, color=DIM, alignment=PP_ALIGN.CENTER)

add_speaker_notes(sl, "The biological coherence is striking. Mce3R is a cholesterol-responsive regulator, "
    "and the predicted targets map to a connected metabolic pathway: import → oxidation → "
    "degradation → carbon metabolism. This functional enrichment is strong evidence against "
    "false positives. The GC content plot confirms hits aren't GC-biased artifacts.")


# ── SLIDE 18: Negative Controls ───────────────────────────────────────
sl = prs.slides.add_slide(blank_layout)
set_slide_bg(sl)
add_textbox(sl, Inches(0.8), Inches(0.4), Inches(11), Inches(0.7),
            "Negative Controls Pass", font_size=30, font_name=HEADER_FONT, color=WHITE, bold=True)
add_accent_line(sl, Inches(0.8), Inches(1.1), Inches(3.5), GREEN_S)

add_textbox(sl, Inches(0.8), Inches(1.5), Inches(11), Inches(0.6),
            "Mce3R does NOT regulate mce1, mce2, or mce4 operons (Santangelo 2008)",
            font_size=18, color=WHITE, bold=True)

neg_results = [
    ("mce1 operon (Rv0169–Rv0178)", "No high-confidence hits", "✓ PASS"),
    ("mce2 operon (Rv0586–Rv0594)", "No high-confidence hits", "✓ PASS"),
    ("mce4 operon (Rv3499c–Rv3513c)", "No high-confidence hits", "✓ PASS"),
]
y = Inches(2.5)
for operon, result, status in neg_results:
    box = add_shape_rect(sl, Inches(1.5), y, Inches(10), Inches(0.8), BG_LIGHT, GREEN_S)
    add_textbox(sl, Inches(1.8), y + Pt(6), Inches(3.5), Inches(0.35), operon,
                font_size=16, color=TEXT)
    add_textbox(sl, Inches(5.5), y + Pt(6), Inches(3), Inches(0.35), result,
                font_size=14, color=DIM)
    add_textbox(sl, Inches(9), y + Pt(6), Inches(2), Inches(0.35), status,
                font_size=16, color=GREEN_S, bold=True, alignment=PP_ALIGN.CENTER)
    y += Inches(1.1)

cbox = add_shape_rect(sl, Inches(2.5), Inches(5.8), Inches(8), Inches(1.0), BG_LIGHT, TEAL)
tf = cbox.text_frame
tf.word_wrap = True
p = tf.paragraphs[0]
p.text = "Pipeline specificity validated — not just finding random TetR sites"
p.font.size = Pt(16)
p.font.color.rgb = ACCENT
p.alignment = PP_ALIGN.CENTER

add_speaker_notes(sl, "Negative controls are as important as positive results. Santangelo 2008 "
    "showed genetically that Mce3R does NOT regulate mce1, mce2, or mce4. "
    "Our pipeline correctly produces no high-confidence hits in those promoters. "
    "This validates the pipeline's specificity — it's not just finding any TetR-like site.")


# ── SLIDE 19: Computational Validation ────────────────────────────────
sl = prs.slides.add_slide(blank_layout)
set_slide_bg(sl)
add_textbox(sl, Inches(0.8), Inches(0.4), Inches(11), Inches(0.7),
            "Computational Validation Approaches", font_size=30, font_name=HEADER_FONT,
            color=WHITE, bold=True)
add_accent_line(sl, Inches(0.8), Inches(1.1), Inches(5), TEAL)

validations = [
    ("Cross-species conservation", "Are sites conserved in M. bovis / M. marinum?"),
    ("Expression data overlap", "Match predicted targets to Δmce3R microarray data"),
    ("Functional enrichment", "Are lipid/cholesterol genes overrepresented? (Fisher's exact test)"),
    ("Motif vs. structure", "Do conserved PWM positions match Arg53/Lys262 DNA contacts?"),
    ("Shuffle test", "Do real motifs outscore shuffled sequences? (controls for GC bias)"),
]
y = Inches(1.5)
for i, (title, desc) in enumerate(validations):
    add_circle_number(sl, Inches(0.8), y, i + 1)
    add_textbox(sl, Inches(1.6), y - Pt(2), Inches(5), Inches(0.35), title,
                font_size=17, color=ACCENT, bold=True)
    add_textbox(sl, Inches(1.6), y + Inches(0.3), Inches(5.5), Inches(0.35), desc,
                font_size=14, color=DIM)
    y += Inches(0.85)

try_add_image(sl, "motif_position_distribution.png", Inches(7), Inches(1.3), width=Inches(5.8))

add_speaker_notes(sl, "Multiple independent lines of evidence strengthen computational predictions "
    "without wet-lab experiments. Cross-species conservation tests evolutionary pressure. "
    "Microarray overlap tests whether predicted targets are differentially expressed in Δmce3R mutants. "
    "Shuffle tests control for GC composition bias in the M. tb genome.")


# ── SLIDE 20: Implications ───────────────────────────────────────────
sl = prs.slides.add_slide(blank_layout)
set_slide_bg(sl)
add_textbox(sl, Inches(0.8), Inches(0.4), Inches(11), Inches(0.7),
            "Implications — Why This Matters", font_size=30, font_name=HEADER_FONT,
            color=WHITE, bold=True)
add_accent_line(sl, Inches(0.8), Inches(1.1), Inches(4), GOLD)

implications = [
    "First genome-wide map of predicted Mce3R binding sites using the asymmetric model",
    "Identifies 24 high-confidence new regulatory targets",
    "Suggests Mce3R controls a broader metabolic network than previously known",
    "Connects cholesterol import → catabolism → carbon metabolism under one regulator",
    "Therapeutic relevance: targeting the Mce3R regulon could address antibiotic persistence",
]
y = Inches(1.5)
for imp in implications:
    c = GOLD if "First" in imp or "Therapeutic" in imp else TEXT
    box = add_shape_rect(sl, Inches(0.8), y, Inches(11.5), Inches(0.7), BG_LIGHT)
    add_shape_rect(sl, Inches(0.8), y, Inches(0.15), Inches(0.7), ACCENT)
    add_textbox(sl, Inches(1.2), y + Pt(6), Inches(11), Inches(0.55), imp, font_size=16, color=c)
    y += Inches(0.9)

add_textbox(sl, Inches(1), Inches(6.2), Inches(11), Inches(0.5),
            "This changes our understanding of Mce3R from a single-operon repressor to a global metabolic coordinator",
            font_size=15, color=DIM, alignment=PP_ALIGN.CENTER)

add_speaker_notes(sl, "Frame the significance. Before this work, Mce3R was thought to regulate "
    "mainly the mce3 operon. Our results suggest it coordinates a much larger metabolic network. "
    "If Mce3R controls cholesterol catabolism genes, targeting this regulon could address "
    "the persister population problem in TB treatment.")


# ── SLIDE 21: Future Directions ───────────────────────────────────────
sl = prs.slides.add_slide(blank_layout)
set_slide_bg(sl)
add_textbox(sl, Inches(0.8), Inches(0.4), Inches(11), Inches(0.7),
            "Future Directions", font_size=30, font_name=HEADER_FONT, color=WHITE, bold=True)
add_accent_line(sl, Inches(0.8), Inches(1.1), Inches(3), TEAL)

futures = [
    ("Experimental Validation",
     "EMSA on top Tier 2 targets in M. smegmatis (BSL-1 safe surrogate)"),
    ("Stochastic Modeling",
     "How does asymmetric operator architecture drive gene expression noise and persister formation?"),
    ("Comparative Genomics",
     "Extend to drug-resistant TB clinical isolates — do operator mutations correlate with resistance?"),
    ("Structural Biology",
     "Predict Mce3R-DNA contacts for novel targets using AlphaFold / molecular dynamics"),
]
y = Inches(1.5)
for i, (title, desc) in enumerate(futures):
    add_circle_number(sl, Inches(0.8), y, i + 1)
    box = add_shape_rect(sl, Inches(1.6), y - Pt(4), Inches(10.5), Inches(1.1), BG_LIGHT, TEAL)
    tf = box.text_frame
    tf.word_wrap = True
    p = tf.paragraphs[0]
    p.text = title
    p.font.size = Pt(17)
    p.font.color.rgb = ACCENT
    p.font.bold = True
    add_para(tf, desc, font_size=13, color=DIM, space_before=Pt(6))
    y += Inches(1.3)

add_speaker_notes(sl, "Show the project has clear next steps. EMSA is the gold-standard experimental "
    "validation. Stochastic modeling connects to Balázsi lab's expertise in gene expression noise. "
    "Comparative genomics could reveal whether operator mutations drive drug resistance. "
    "AlphaFold modeling could predict new Mce3R-DNA interactions without crystallography.")


# ── SLIDE 22: Methods Summary ─────────────────────────────────────────
sl = prs.slides.add_slide(blank_layout)
set_slide_bg(sl)
add_textbox(sl, Inches(0.8), Inches(0.4), Inches(11), Inches(0.7),
            "Methods Summary — Code Architecture", font_size=28, font_name=HEADER_FONT, color=WHITE, bold=True)
add_accent_line(sl, Inches(0.8), Inches(1.1), Inches(4.5), TEAL)

# Left column — scripts
script_box = add_shape_rect(sl, Inches(0.8), Inches(1.3), Inches(5.8), Inches(5.5), BG_LIGHT, TEAL)
tf = script_box.text_frame
tf.word_wrap = True
p = tf.paragraphs[0]
p.text = "Pipeline Scripts (10 Python modules)"
p.font.size = Pt(14)
p.font.color.rgb = ACCENT
p.font.bold = True

scripts = [
    ("main.py", "Orchestrator — argparse --steps/--simulate/--skip-existing flags"),
    ("download_genome.py", "NCBI Entrez → H37Rv GenBank, validates >4 Mb / >3K CDS"),
    ("extract_promoters.py", "CDS interval tree, 300 bp upstream trim, divergent IGR extraction"),
    ("fetch_sequences.py", "Downloads orthologs (Mtb/bovis/marinum) + tier-specific promoters"),
    ("run_meme.py", "MEME subprocess or simulation fallback, dual-motif PWM output"),
    ("run_fimo.py", "Log-odds scanner, empirical p-values, both strands, fimo.tsv output"),
    ("analyze_results.py", "Palindrome scoring, spacing analysis, tier assignment, ranking"),
    ("motif_models.py", "Composite score, architecture classification, model comparison"),
    ("visualize.py", "9 figures (PNG 300 DPI + SVG): score dist, logo, architecture, etc."),
    ("utils.py", "Shared: logging, path setup, GC content, FASTA loader, reverse_complement"),
]
for name, desc in scripts:
    add_para(tf, f"  {name}", font_size=11, color=GOLD, bold=True, space_before=Pt(5))
    add_para(tf, f"    {desc}", font_size=10, color=DIM, space_before=Pt(1))

# Right column — environment and key parameters
env_box = add_shape_rect(sl, Inches(6.8), Inches(1.3), Inches(5.5), Inches(2.5), BG_LIGHT, TEAL)
tfe = env_box.text_frame
tfe.word_wrap = True
pe = tfe.paragraphs[0]
pe.text = "Environment & Dependencies"
pe.font.size = Pt(14)
pe.font.color.rgb = ACCENT
pe.font.bold = True
env_items = [
    "Python 3.10 (Miniconda conda environment)",
    "MEME Suite 5.5.9 (bioconda channel)",
    "BioPython 1.86 — GenBank parsing, NCBI Entrez",
    "pandas 2.x / NumPy 2.x — data manipulation",
    "matplotlib 3.10 / seaborn 0.13 — visualization",
]
for item in env_items:
    add_para(tfe, f"  {item}", font_size=11, color=TEXT, space_before=Pt(4))

# Key parameters
param_box = add_shape_rect(sl, Inches(6.8), Inches(4.0), Inches(5.5), Inches(2.8), BG_LIGHT, GOLD)
tfp = param_box.text_frame
tfp.word_wrap = True
pp = tfp.paragraphs[0]
pp.text = "Key Parameters"
pp.font.size = Pt(14)
pp.font.color.rgb = GOLD
pp.font.bold = True
params = [
    ("Promoter length:", "300 bp upstream (trimmed at adjacent CDS)"),
    ("MEME width:", "-minw 20 -maxw 30 (match ~25 bp half-site)"),
    ("MEME modes:", "-mod zoops, then -mod anr"),
    ("FIMO p-value:", "< 0.0001 (1 in 10,000)"),
    ("Max gap:", "80 bp (site-pairing search window)"),
    ("Validated spacing:", "48–58 bp (±5 bp around 53 bp reference)"),
    ("Composite weights:", "40% FIMO / 25% arch / 20% asym / 15% tier"),
]
for label, val in params:
    add_para(tfp, f"  {label} {val}", font_size=10, color=TEXT, space_before=Pt(3))

add_speaker_notes(sl, "Full technical specifications for reproducibility. The pipeline is 7 Python modules "
    "totaling ~2,500 lines. Key design choices: MEME width 20-30 bp is informed by the structural "
    "data showing ~25 bp half-sites. FIMO p-value 0.0001 balances sensitivity vs. specificity. "
    "Max gap of 80 bp ensures we don't miss wider-than-expected paired sites. "
    "The composite 40/25/20/15 weights were tuned so the known validated operator ranks #1.")


# ── SLIDE 23: Key Literature ──────────────────────────────────────────
sl = prs.slides.add_slide(blank_layout)
set_slide_bg(sl)
add_textbox(sl, Inches(0.8), Inches(0.4), Inches(11), Inches(0.7),
            "Key Literature", font_size=30, font_name=HEADER_FONT, color=WHITE, bold=True)
add_accent_line(sl, Inches(0.8), Inches(1.1), Inches(3), TEAL)

refs = [
    ("Panagoda, Balázsi & Sampson (2024)",
     "ACS Chemical Biology — Cryo-EM structure of Mce3R, asymmetric operator discovery"),
    ("Santangelo et al. (2009)",
     "Microbiology 155:2245 — Mce3R regulon characterization, second regulatory region"),
    ("Santangelo, Blanco et al. (2008)",
     "BMC Microbiology 8:38 — Genetic characterization, mce1/2/4 not Mce3R-regulated"),
    ("Bailey et al. (2015)",
     "Nucleic Acids Research — MEME Suite: tools for motif discovery and searching"),
    ("Grant, Bailey & Noble (2011)",
     "Bioinformatics — FIMO: scanning for occurrences of a given motif"),
]
y = Inches(1.5)
for i, (authors, desc) in enumerate(refs):
    box = add_shape_rect(sl, Inches(0.8), y, Inches(11.5), Inches(0.85), BG_LIGHT if i % 2 == 0 else BG)
    add_textbox(sl, Inches(1), y + Pt(4), Inches(10.5), Inches(0.35), authors,
                font_size=15, color=ACCENT, bold=True)
    add_textbox(sl, Inches(1), y + Inches(0.35), Inches(10.5), Inches(0.4), desc,
                font_size=13, color=DIM)
    y += Inches(0.95)

add_speaker_notes(sl, "Key references. The 2024 ACS Chem Biol paper is the most important — "
    "it established the asymmetric operator model that our pipeline is built around.")


# ── SLIDE 24: Acknowledgments ─────────────────────────────────────────
sl = prs.slides.add_slide(blank_layout)
set_slide_bg(sl)
add_textbox(sl, Inches(0.8), Inches(0.4), Inches(11), Inches(0.7),
            "Acknowledgments", font_size=30, font_name=HEADER_FONT, color=WHITE, bold=True)
add_accent_line(sl, Inches(0.8), Inches(1.1), Inches(3), TEAL)

acks = [
    "Dr. Gábor Balázsi — Stony Brook University, Mce3R structural work",
    "Dr. Nicole Sampson — Stony Brook / University of Rochester, Mce3R biochemistry",
]
y = Inches(2.0)
for ack in acks:
    add_shape_rect(sl, Inches(3), y + Pt(5), Inches(0.12), Inches(0.12), ACCENT)
    add_textbox(sl, Inches(3.4), y, Inches(7), Inches(0.4), ack, font_size=17, color=TEXT)
    y += Inches(0.7)

add_textbox(sl, Inches(3), Inches(4.5), Inches(7), Inches(0.5),
            "Great Neck South High School", font_size=20, color=DIM,
            alignment=PP_ALIGN.CENTER, bold=True)

add_speaker_notes(sl, "Thank your mentors and acknowledge the scientific community.")


# ── SLIDE 25: Questions / Contact ─────────────────────────────────────
sl = prs.slides.add_slide(blank_layout)
set_slide_bg(sl)
add_accent_line(sl, Inches(4), Inches(2.5), Inches(5.3), TEAL)
add_textbox(sl, Inches(1), Inches(2.8), Inches(11.3), Inches(1.2),
            "Thank You", font_size=48, font_name=HEADER_FONT, color=WHITE,
            bold=True, alignment=PP_ALIGN.CENTER)
add_textbox(sl, Inches(1), Inches(4.2), Inches(11.3), Inches(0.6),
            "Questions?", font_size=28, color=ACCENT, alignment=PP_ALIGN.CENTER)
add_accent_line(sl, Inches(4), Inches(5.0), Inches(5.3), GOLD)
add_textbox(sl, Inches(1), Inches(5.5), Inches(11.3), Inches(0.5),
            "Aayan Alwani  •  Great Neck South High School", font_size=16, color=DIM,
            alignment=PP_ALIGN.CENTER)

add_speaker_notes(sl, "Thank the audience and open for questions.")


# ══════════════════════════════════════════════════════════════════════
#  SAVE
# ══════════════════════════════════════════════════════════════════════
prs.save(OUT_PATH)
print(f"Presentation saved to: {OUT_PATH}")
print(f"Slides: {len(prs.slides)}")
