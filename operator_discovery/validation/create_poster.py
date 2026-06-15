#!/usr/bin/env python3
"""
Create Mce3R research poster board as a single-slide 48x36" PowerPoint.
Tri-fold layout: 12" left | 24" center | 12" right.
"""

from pptx import Presentation
from pptx.util import Inches, Pt, Emu
from pptx.dml.color import RGBColor
from pptx.enum.text import PP_ALIGN, MSO_ANCHOR
from pptx.enum.shapes import MSO_SHAPE
import os

# ── Constants ──────────────────────────────────────────────────────────
W, H = 48, 36  # poster dimensions in inches
NAVY = RGBColor(0x0B, 0x1D, 0x3A)
TEAL = RGBColor(0x0D, 0x94, 0x88)
SEAFOAM = RGBColor(0x14, 0xB8, 0xA6)
CHARCOAL = RGBColor(0x1E, 0x29, 0x3B)
LIGHT_GRAY = RGBColor(0xE2, 0xE8, 0xF0)
WHITE = RGBColor(0xFF, 0xFF, 0xFF)
BG_COLOR = RGBColor(0xF8, 0xF9, 0xFA)
CORAL = RGBColor(0xEF, 0x44, 0x44)
AMBER = RGBColor(0xF5, 0x9E, 0x0B)

FIG_DIR = 'results/figures'
POSTER_FIG_DIR = 'results/figures/poster'

prs = Presentation()
prs.slide_width = Inches(W)
prs.slide_height = Inches(H)
slide = prs.slides.add_slide(prs.slide_layouts[6])  # blank layout

# Set slide background
bg = slide.background
fill = bg.fill
fill.solid()
fill.fore_color.rgb = BG_COLOR


def add_box(left, top, width, height, fill_color=None, border_color=None, border_width=Pt(0)):
    """Add a rectangle shape."""
    shape = slide.shapes.add_shape(MSO_SHAPE.ROUNDED_RECTANGLE, Inches(left), Inches(top),
                                   Inches(width), Inches(height))
    shape.fill.solid()
    shape.fill.fore_color.rgb = fill_color or WHITE
    ln = shape.line
    if border_color:
        ln.color.rgb = border_color
        ln.width = border_width
    else:
        ln.fill.background()
    # Small corner radius
    shape.adjustments[0] = 0.02
    return shape


def add_text_box(left, top, width, height, text, font_size=20, bold=False,
                 color=CHARCOAL, alignment=PP_ALIGN.LEFT, font_name='Calibri',
                 italic=False):
    """Add a text box with styled text."""
    txBox = slide.shapes.add_textbox(Inches(left), Inches(top),
                                      Inches(width), Inches(height))
    tf = txBox.text_frame
    tf.word_wrap = True
    p = tf.paragraphs[0]
    p.text = text
    p.font.size = Pt(font_size)
    p.font.bold = bold
    p.font.color.rgb = color
    p.font.name = font_name
    p.alignment = alignment
    p.font.italic = italic
    return txBox


def add_section_header(left, top, width, text, font_size=30):
    """Add a navy section header box."""
    box = add_box(left, top, width, 0.6, fill_color=NAVY)
    tf = box.text_frame
    tf.word_wrap = True
    p = tf.paragraphs[0]
    p.text = text
    p.font.size = Pt(font_size)
    p.font.bold = True
    p.font.color.rgb = WHITE
    p.font.name = 'Trebuchet MS'
    p.alignment = PP_ALIGN.CENTER
    tf.margin_left = Inches(0.15)
    tf.margin_right = Inches(0.15)
    tf.margin_top = Inches(0.05)
    tf.margin_bottom = Inches(0.05)
    return box


def add_bullet_box(left, top, width, height, bullets, font_size=18, color=CHARCOAL,
                   bold_phrases=None):
    """Add a text box with bullet points."""
    txBox = slide.shapes.add_textbox(Inches(left), Inches(top),
                                      Inches(width), Inches(height))
    tf = txBox.text_frame
    tf.word_wrap = True

    for i, bullet in enumerate(bullets):
        if i == 0:
            p = tf.paragraphs[0]
        else:
            p = tf.add_paragraph()
        p.text = bullet
        p.font.size = Pt(font_size)
        p.font.color.rgb = color
        p.font.name = 'Calibri'
        p.space_after = Pt(4)
        p.level = 0
        # Add bullet character
        p.text = '  ' + bullet

    return txBox


def add_image(left, top, width, path):
    """Add an image, return the shape."""
    full = path if os.path.isabs(path) else os.path.join(os.getcwd(), path)
    if not os.path.exists(full):
        print(f'  WARNING: missing {full}')
        return None
    pic = slide.shapes.add_picture(full, Inches(left), Inches(top), width=Inches(width))
    return pic


def add_caption(left, top, width, text, font_size=14):
    """Add a figure caption."""
    return add_text_box(left, top, width, 0.5, text, font_size=font_size,
                        italic=True, color=CHARCOAL)


# ═══════════════════════════════════════════════════════════════════════
# TITLE BANNER (full width)
# ═══════════════════════════════════════════════════════════════════════
title_box = add_box(0, 0, W, 3.0, fill_color=NAVY)
add_text_box(0.5, 0.3, W - 1, 1.5,
             'Stochastic Gene Expression Noise in the Mce3R Autoregulatory Circuit\n'
             'as a Driver of Antibiotic Persistence in Mycobacterium tuberculosis',
             font_size=52, bold=True, color=WHITE, alignment=PP_ALIGN.CENTER,
             font_name='Trebuchet MS')
add_text_box(0.5, 2.0, W - 1, 0.4,
             'Aayan Alwani',
             font_size=32, bold=True, color=SEAFOAM, alignment=PP_ALIGN.CENTER,
             font_name='Trebuchet MS')
add_text_box(0.5, 2.45, W - 1, 0.4,
             'Great Neck South High School, Great Neck, NY',
             font_size=22, bold=False, color=LIGHT_GRAY, alignment=PP_ALIGN.CENTER,
             font_name='Calibri')

# Panel margins
L_LEFT = 0.4          # left panel starts
L_WIDTH = 11.2        # left panel width
C_LEFT = 12.3         # center panel starts
C_WIDTH = 23.4        # center panel width
R_LEFT = 36.2         # right panel starts
R_WIDTH = 11.4        # right panel width

# Vertical divider lines (subtle)
for x in [12.0, 36.0]:
    line = slide.shapes.add_shape(MSO_SHAPE.RECTANGLE, Inches(x), Inches(3.2),
                                   Inches(0.03), Inches(H - 3.5))
    line.fill.solid()
    line.fill.fore_color.rgb = LIGHT_GRAY
    line.line.fill.background()

y_start = 3.4  # content starts below title

# ═══════════════════════════════════════════════════════════════════════
# LEFT PANEL — Background & Methods
# ═══════════════════════════════════════════════════════════════════════

# ── Background ──
y = y_start
add_section_header(L_LEFT, y, L_WIDTH, 'BACKGROUND')
y += 0.75
add_bullet_box(L_LEFT, y, L_WIDTH, 2.5, [
    'TB kills 1.3 million people/year (WHO 2024)',
    'Treatment: 6-9 months due to persister cells',
    'Persisters survive antibiotics without resistance mutations',
    'Caused by random gene expression noise',
    'Mce3R: TetR-family repressor of cholesterol survival genes',
    'Mce3R deletion increases persister frequency (Srivastava 2023)',
], font_size=18)
y += 2.7

# Figure 1: Repression mechanism
add_image(L_LEFT + 0.3, y, L_WIDTH - 0.6, f'{POSTER_FIG_DIR}/fig1_repression_mechanism.png')
y += 4.8
add_caption(L_LEFT, y, L_WIDTH,
            'Figure 1. Mce3R repression mechanism. Figure created by student author.')
y += 0.55

# ── Asymmetric Operator Discovery ──
y += 0.15
add_section_header(L_LEFT, y, L_WIDTH, 'THE ASYMMETRIC OPERATOR')
y += 0.75
add_bullet_box(L_LEFT, y, L_WIDTH, 1.8, [
    '2024 cryo-EM at 2.51 A (Panagoda, Balazsi & Sampson)',
    'Two 25 bp non-palindromic half-sites, 53 bp spacer',
    'Downstream: Kd = 2.4 nM (strong); upstream: weaker',
    'First asymmetric operator in TetR family (unprecedented)',
], font_size=18)
y += 1.9

# Figure 2: Asymmetric vs symmetric
add_image(L_LEFT + 0.3, y, L_WIDTH - 0.6, f'{POSTER_FIG_DIR}/fig2_asymmetric_vs_symmetric.png')
y += 5.3
add_caption(L_LEFT, y, L_WIDTH,
            'Figure 2. Asymmetric vs. symmetric operator architecture. '
            'Figure created by student author.')
y += 0.55

# ── Hypothesis ──
y += 0.15
add_section_header(L_LEFT, y, L_WIDTH, 'HYPOTHESIS')
y += 0.75

# Hypothesis text in a teal box
hyp_box = add_box(L_LEFT + 0.2, y, L_WIDTH - 0.4, 1.8, fill_color=TEAL)
tf = hyp_box.text_frame
tf.word_wrap = True
tf.margin_left = Inches(0.15)
tf.margin_right = Inches(0.15)
tf.margin_top = Inches(0.1)
p = tf.paragraphs[0]
p.text = ('Asymmetric binding sites create a THREE-STATE gene expression system '
          '(OFF / INTERMEDIATE / ON) that generates more persisters than '
          'symmetric sites (OFF / ON only)')
p.font.size = Pt(20)
p.font.bold = True
p.font.color.rgb = WHITE
p.font.name = 'Calibri'
p.alignment = PP_ALIGN.CENTER
y += 2.0

# Figure 3: Four operator states
add_image(L_LEFT + 0.3, y, L_WIDTH - 0.6, f'{POSTER_FIG_DIR}/fig3_four_operator_states.png')
y += 6.2
add_caption(L_LEFT, y, L_WIDTH,
            'Figure 3. Four operator microstates. Intermediate states (2,3) only exist '
            'because sites have different strengths. Figure created by student author.')


# ═══════════════════════════════════════════════════════════════════════
# CENTER PANEL — Pipeline & Stochastic Modeling
# ═══════════════════════════════════════════════════════════════════════

y = y_start

# ── Computational Pipeline ──
add_section_header(C_LEFT, y, C_WIDTH, 'COMPUTATIONAL PIPELINE (AIM 1 — COMPLETED)')
y += 0.75
add_bullet_box(C_LEFT, y, C_WIDTH * 0.48, 2.0, [
    'Scanned 2,351 promoters across M. tuberculosis H37Rv',
    'MEME Suite: de novo motif discovery with grid-search width optimization (w=20-28)',
    'FIMO: genome-wide scanning with dual motif models (palindromic + asymmetric)',
    'Composite scoring: FIMO + architecture + asymmetry + biological prior',
], font_size=17)

# Score distribution figure (right half)
add_image(C_LEFT + C_WIDTH * 0.50, y - 0.2, C_WIDTH * 0.48,
          f'{FIG_DIR}/score_distribution.png')
add_caption(C_LEFT + C_WIDTH * 0.50, y + 4.0, C_WIDTH * 0.48,
            'Figure 4. Composite score distribution across 451 predicted binding sites.')
y += 4.5

# ── Pipeline Results ──
add_section_header(C_LEFT, y, C_WIDTH, 'PIPELINE RESULTS')
y += 0.75

# Results bullets left side
add_bullet_box(C_LEFT, y, C_WIDTH * 0.48, 2.0, [
    '451 predicted Mce3R binding sites',
    '   6 Tier 1 (known) | 23 Tier 2 (plausible) | 422 Tier 3',
    'Top hit: IGR_Rv1963c_Rv1964 (known operator) — validates pipeline',
    '5 sites with validated paired asymmetric architecture (40-65 bp)',
    '98% of sites prefer asymmetric model over palindromic',
], font_size=17)

# Model comparison figure right side
add_image(C_LEFT + C_WIDTH * 0.50, y - 0.2, C_WIDTH * 0.48,
          f'{FIG_DIR}/model_comparison.png')
add_caption(C_LEFT + C_WIDTH * 0.50, y + 3.8, C_WIDTH * 0.48,
            'Figure 5. 98% of predicted sites prefer the asymmetric binding model.')
y += 4.3

# Spacing histogram + architecture (side by side)
half_w = C_WIDTH * 0.48
add_image(C_LEFT, y, half_w, f'{FIG_DIR}/spacing_histogram.png')
add_caption(C_LEFT, y + 3.6, half_w,
            'Figure 6. Inter-site spacing. Vertical line = validated 53 bp.')

add_image(C_LEFT + C_WIDTH * 0.52, y, half_w, f'{FIG_DIR}/operator_architecture_diagram.png')
add_caption(C_LEFT + C_WIDTH * 0.52, y + 3.6, half_w,
            'Figure 7. Validated paired asymmetric operator architecture.')
y += 4.2

# ── Stochastic Modeling ──
add_section_header(C_LEFT, y, C_WIDTH, 'STOCHASTIC MODELING (AIM 2)')
y += 0.75
add_bullet_box(C_LEFT, y, C_WIDTH * 0.48, 2.5, [
    'Gillespie stochastic simulation (SSA) of Mce3R circuit',
    '6 reactions: transcription, translation, mRNA degradation,',
    '   protein dilution, Mce3R binding, Mce3R unbinding',
    'Parameters from published data:',
    '   Kd = 2.4 nM | mRNA t1/2 = 5 min | doubling = 20 hrs',
    '10,000 virtual cells per condition',
], font_size=17)

# Gillespie flowchart
add_image(C_LEFT + C_WIDTH * 0.50, y - 0.3, C_WIDTH * 0.48,
          f'{POSTER_FIG_DIR}/fig8_gillespie_flowchart.png')
add_caption(C_LEFT + C_WIDTH * 0.50, y + 2.4, C_WIDTH * 0.48,
            'Figure 8. Gillespie algorithm. Figure created by student author.')
y += 3.2

# ── Four Simulations ──
add_section_header(C_LEFT, y, C_WIDTH, 'THE FOUR SIMULATIONS')
y += 0.7
add_image(C_LEFT + 1.0, y, C_WIDTH - 2.0,
          f'{POSTER_FIG_DIR}/fig9_four_simulations.png')
y += 5.5
add_caption(C_LEFT + 1.0, y, C_WIDTH - 2.0,
            'Figure 9. Four operator architectures tested. '
            'Figure created by student author.')


# ═══════════════════════════════════════════════════════════════════════
# RIGHT PANEL — Results, Discussion, Future Directions
# ═══════════════════════════════════════════════════════════════════════

y = y_start

# ── Key Finding ──
add_section_header(R_LEFT, y, R_WIDTH, 'KEY FINDING')
y += 0.75

# Trimodal distributions figure
add_image(R_LEFT + 0.2, y, R_WIDTH - 0.4,
          f'{POSTER_FIG_DIR}/fig10_trimodal_distributions.png')
y += 7.8
add_caption(R_LEFT, y, R_WIDTH,
            'Figure 10. Predicted protein distributions. Asymmetric operators '
            'produce trimodal distributions; intermediate state = persisters. '
            'Figure created by student author.')
y += 0.7

# Key finding bullet points
add_bullet_box(R_LEFT, y, R_WIDTH, 1.5, [
    'Asymmetric operator: 3 peaks (OFF / INTERMEDIATE / ON)',
    'Symmetric operator: 2 peaks only (OFF / ON)',
    'The INTERMEDIATE state is the persister subpopulation',
    'Predicted persister fraction (~0.1-1%) matches experiments',
], font_size=17)
y += 1.8

# ── Validation ──
add_section_header(R_LEFT, y, R_WIDTH, 'VALIDATION')
y += 0.75
add_bullet_box(R_LEFT, y, R_WIDTH, 1.5, [
    'Cross-referenced with Santangelo 2009 microarray',
    'Cross-referenced with persister transcriptomics (Keren 2011)',
    'Negative controls PASS: no hits at mce1/mce2/mce4',
    'Pipeline validates: top hit is the known operator',
], font_size=17)
y += 1.7

# Palindrome scores figure
add_image(R_LEFT + 0.5, y, R_WIDTH - 1.0, f'{FIG_DIR}/palindrome_scores.png')
y += 4.0
add_caption(R_LEFT, y, R_WIDTH,
            'Figure 11. Palindrome score distribution confirms '
            'non-palindromic architecture.')
y += 0.6

# ── Novel Targets ──
add_section_header(R_LEFT, y, R_WIDTH, 'NOVEL TARGETS DISCOVERED')
y += 0.7

# Small table of targets
targets = [
    ('eccA3', 'ESX-3 secretion', '54 bp spacing'),
    ('choD', 'Cholesterol oxidase', 'First step in cholesterol degradation'),
    ('icl1', 'Isocitrate lyase', 'Glyoxylate shunt'),
    ('vapC33', 'Toxin-antitoxin', 'Directly linked to persistence'),
]

for i, (gene, desc, note) in enumerate(targets):
    row_y = y + i * 0.55
    # Gene name
    add_text_box(R_LEFT + 0.2, row_y, 1.5, 0.5, gene, font_size=17, bold=True,
                 color=TEAL, font_name='Calibri')
    # Description
    add_text_box(R_LEFT + 1.8, row_y, 4.5, 0.5, f'{desc} — {note}',
                 font_size=15, color=CHARCOAL)
y += 2.5

# ── Conclusions ──
add_section_header(R_LEFT, y, R_WIDTH, 'CONCLUSIONS')
y += 0.75
add_bullet_box(R_LEFT, y, R_WIDTH, 2.0, [
    'Asymmetric operator architecture produces trimodal noise (novel)',
    'Intermediate expression state is a "persistence dimmer switch"',
    'Mce3R regulates a broader network than previously known',
    'DNA binding site shape controls persister frequency',
], font_size=17)
y += 1.8

# ── Future Directions ──
add_section_header(R_LEFT, y, R_WIDTH, 'FUTURE DIRECTIONS')
y += 0.75
add_bullet_box(R_LEFT, y, R_WIDTH, 1.5, [
    'Experimental: EMSA on top targets in M. smegmatis',
    'Flow cytometry: measure noise with fluorescent reporters',
    'Therapeutic: design noise suppressors to reduce persisters',
], font_size=17)
y += 1.6

# ── Acknowledgments ──
add_section_header(R_LEFT, y, R_WIDTH, 'ACKNOWLEDGMENTS', font_size=24)
y += 0.65
add_bullet_box(R_LEFT, y, R_WIDTH, 0.8, [
    'Dr. Gabor Balazsi, Stony Brook University',
    'Dr. Nicole S. Sampson, University of Rochester',
    'Great Neck South High School',
], font_size=14)
y += 1.0

# ── References ──
add_section_header(R_LEFT, y, R_WIDTH, 'REFERENCES', font_size=24)
y += 0.65
refs = (
    'Panagoda et al. 2024, ACS Chem. Biol. 19:2580\n'
    'Santangelo et al. 2009, Microbiology 155:2245\n'
    'Srivastava et al. 2023, Microbiol. Spectrum 11:e00235\n'
    'Balazsi et al. 2011, Cell 144:910\n'
    'Balaban et al. 2004, Science 305:1622\n'
    'Quigley & Lewis 2022, Microbiol. Spectrum\n'
    'Gillespie 1977, J. Phys. Chem. 81:2340'
)
add_text_box(R_LEFT + 0.2, y, R_WIDTH - 0.4, 2.0, refs, font_size=13,
             color=CHARCOAL)

# ═══════════════════════════════════════════════════════════════════════
# Save
# ═══════════════════════════════════════════════════════════════════════
out_path = 'mce3r_poster.pptx'
prs.save(out_path)
print(f'\nPoster saved to: {out_path}')
print(f'Dimensions: {W}" x {H}" (standard tri-fold)')
