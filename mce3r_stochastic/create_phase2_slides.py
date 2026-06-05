#!/usr/bin/env python3
"""Generate Phase 2 presentation: figures + PPTX."""

import os
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
from pptx import Presentation
from pptx.util import Inches, Pt, Emu
from pptx.dml.color import RGBColor
from pptx.enum.text import PP_ALIGN, MSO_ANCHOR
from pptx.enum.shapes import MSO_SHAPE
from pptx.oxml.ns import qn
from lxml import etree

# ============================================================
# PATHS
# ============================================================
BASE = "/Users/aayanalwani/tb project/mce3r_stochastic"
FIGS = f"{BASE}/results/phase2/slide_figures"
EXISTING_FIGS = f"{BASE}/results/figures"
OUT_PPTX = f"{BASE}/Phase2_Presentation.pptx"
os.makedirs(FIGS, exist_ok=True)

# ============================================================
# COLORS (shared)
# ============================================================
TEAL_HEX = '#0D9488'
DARK_HEX = '#1E293B'
LIGHT_HEX = '#F8FAFC'
GRAY_HEX = '#64748B'
ACCENT_HEX = '#DC2626'
WHITE_HEX = '#FFFFFF'

# ============================================================
# FIGURE 1: operator_states.png
# ============================================================
def make_operator_states():
    fig, ax = plt.subplots(figsize=(10, 7))
    ax.set_xlim(-1, 11)
    ax.set_ylim(-1, 9)
    ax.axis('off')
    fig.patch.set_facecolor('#F8FAFC')

    # State positions
    positions = {
        0: (2, 7),   # State 0 top-left
        1: (8, 7),   # State 1 top-right
        2: (2, 2),   # State 2 bottom-left
        3: (8, 2),   # State 3 bottom-right
    }

    labels = {
        0: "State 0\nUnbound",
        1: "State 1\nStrong Bound",
        2: "State 2\nWeak Bound",
        3: "State 3\nBoth Bound",
    }

    txn_rates = {
        0: "k_txn = k_max\n(0.150 /min)",
        1: "k_txn = k_max(1-0.85)\n(0.0225 /min)",
        2: "k_txn = k_max(1-0.50)\n(0.075 /min)",
        3: "k_txn = k_max(1-0.85)(1-0.50)\n(0.01125 /min)",
    }

    box_colors = {0: '#10B981', 1: '#F59E0B', 2: '#3B82F6', 3: '#EF4444'}

    # Draw state boxes
    for state, (x, y) in positions.items():
        box = FancyBboxPatch((x - 1.4, y - 1.0), 2.8, 2.0,
                             boxstyle="round,pad=0.15",
                             facecolor=box_colors[state], alpha=0.15,
                             edgecolor=box_colors[state], linewidth=2.5)
        ax.add_patch(box)
        ax.text(x, y + 0.35, labels[state], ha='center', va='center',
                fontsize=13, fontweight='bold', color=DARK_HEX)
        ax.text(x, y - 0.45, txn_rates[state], ha='center', va='center',
                fontsize=9, color=GRAY_HEX, style='italic')

    # Arrows: horizontal
    arrow_kw = dict(arrowstyle='->', linewidth=2, mutation_scale=18)

    # 0 -> 1 (bind strong)
    ax.annotate('', xy=(6.4, 7.3), xytext=(3.6, 7.3),
                arrowprops=dict(**arrow_kw, color='#F59E0B'))
    ax.text(5, 7.7, r'$k_{on} \cdot [R]$', ha='center', fontsize=10, color='#F59E0B', fontweight='bold')
    ax.text(5, 7.35, 'bind strong', ha='center', fontsize=8, color=GRAY_HEX)

    # 1 -> 0 (unbind strong)
    ax.annotate('', xy=(3.6, 6.7), xytext=(6.4, 6.7),
                arrowprops=dict(**arrow_kw, color='#F59E0B', linestyle='dashed'))
    ax.text(5, 6.25, r'$k_{off,strong}$', ha='center', fontsize=10, color='#F59E0B')

    # 0 -> 2 (bind weak)
    ax.annotate('', xy=(2.3, 3.2), xytext=(2.3, 5.8),
                arrowprops=dict(**arrow_kw, color='#3B82F6'))
    ax.text(1.2, 4.5, r'$k_{on} \cdot [R]$', ha='center', fontsize=10, color='#3B82F6', fontweight='bold',
            rotation=90)

    # 2 -> 0 (unbind weak)
    ax.annotate('', xy=(1.7, 5.8), xytext=(1.7, 3.2),
                arrowprops=dict(**arrow_kw, color='#3B82F6', linestyle='dashed'))
    ax.text(0.6, 4.5, r'$k_{off,weak}$', ha='center', fontsize=10, color='#3B82F6',
            rotation=90)

    # 2 -> 3 (bind strong when weak occupied)
    ax.annotate('', xy=(6.4, 2.3), xytext=(3.6, 2.3),
                arrowprops=dict(**arrow_kw, color='#EF4444'))
    ax.text(5, 2.7, r'$k_{on} \cdot [R]$', ha='center', fontsize=10, color='#EF4444', fontweight='bold')
    ax.text(5, 2.35, 'bind strong', ha='center', fontsize=8, color=GRAY_HEX)

    # 3 -> 2 (unbind strong)
    ax.annotate('', xy=(3.6, 1.7), xytext=(6.4, 1.7),
                arrowprops=dict(**arrow_kw, color='#EF4444', linestyle='dashed'))
    ax.text(5, 1.25, r'$k_{off,strong}$', ha='center', fontsize=10, color='#EF4444')

    # 1 -> 3 (bind weak when strong occupied)
    ax.annotate('', xy=(8.3, 3.2), xytext=(8.3, 5.8),
                arrowprops=dict(**arrow_kw, color='#EF4444'))
    ax.text(9.5, 4.5, r'$k_{on} \cdot [R]$', ha='center', fontsize=10, color='#EF4444', fontweight='bold',
            rotation=90)

    # 3 -> 1 (unbind weak when both occupied)
    ax.annotate('', xy=(7.7, 5.8), xytext=(7.7, 3.2),
                arrowprops=dict(**arrow_kw, color='#EF4444', linestyle='dashed'))
    ax.text(6.6, 4.5, r'$k_{off,weak}$', ha='center', fontsize=10, color='#EF4444',
            rotation=90)

    ax.set_title('4-State Operator Model', fontsize=18, fontweight='bold',
                 color=DARK_HEX, pad=15)

    plt.tight_layout()
    plt.savefig(f"{FIGS}/operator_states.png", dpi=200, bbox_inches='tight',
                facecolor='#F8FAFC')
    plt.close()
    print("  [OK] operator_states.png")


# ============================================================
# FIGURE 2: gillespie_flowchart.png
# ============================================================
def make_gillespie_flowchart():
    fig, ax = plt.subplots(figsize=(10, 7))
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 10)
    ax.axis('off')
    fig.patch.set_facecolor('#F8FAFC')

    steps = [
        ("Initialize", "Set t=0, mRNA=0,\nprotein=0, op=0", '#10B981'),
        ("Compute\nPropensities", "Calculate a_i for\nall 12 reactions", '#0D9488'),
        ("Draw Waiting\nTime", r"$\tau = -\ln(r_1) / a_0$" + "\nfrom Exponential", '#3B82F6'),
        ("Select\nReaction", r"Draw $r_2$, find $j$:" + "\ncumsum(a) >= r2 * a0", '#7C3AED'),
        ("Update\nState", "Execute reaction j:\nchange op/mRNA/protein", '#F59E0B'),
    ]

    # Terminal
    terminal_y = 1.2
    check_label = "t < t_max?"
    check_y = 2.5

    box_w, box_h = 2.8, 1.3
    start_y = 9.0
    x_center = 5.0

    for i, (title, desc, color) in enumerate(steps):
        y = start_y - i * 1.35
        box = FancyBboxPatch((x_center - box_w/2, y - box_h/2), box_w, box_h,
                             boxstyle="round,pad=0.12",
                             facecolor=color, alpha=0.12,
                             edgecolor=color, linewidth=2)
        ax.add_patch(box)
        ax.text(x_center, y + 0.2, title, ha='center', va='center',
                fontsize=12, fontweight='bold', color=color)
        ax.text(x_center, y - 0.3, desc, ha='center', va='center',
                fontsize=8.5, color=GRAY_HEX)
        # Arrow to next
        if i < len(steps) - 1:
            ax.annotate('', xy=(x_center, y - box_h/2 - 0.05),
                        xytext=(x_center, y - box_h/2 - 0.0),
                        arrowprops=dict(arrowstyle='->', linewidth=1.5,
                                       color=GRAY_HEX, mutation_scale=14))

    # Decision diamond (t < t_max?)
    last_step_y = start_y - (len(steps) - 1) * 1.35
    diamond_y = last_step_y - 1.35
    diamond = plt.Polygon([(x_center, diamond_y + 0.55),
                           (x_center + 1.2, diamond_y),
                           (x_center, diamond_y - 0.55),
                           (x_center - 1.2, diamond_y)],
                          facecolor='#EF4444', alpha=0.12,
                          edgecolor='#EF4444', linewidth=2)
    ax.add_patch(diamond)
    ax.text(x_center, diamond_y, check_label, ha='center', va='center',
            fontsize=11, fontweight='bold', color='#EF4444')

    # Arrow from Update to diamond
    ax.annotate('', xy=(x_center, diamond_y + 0.55),
                xytext=(x_center, last_step_y - box_h/2),
                arrowprops=dict(arrowstyle='->', linewidth=1.5,
                                color=GRAY_HEX, mutation_scale=14))

    # Yes arrow: loop back to step 1 (Compute Propensities)
    step1_y = start_y - 1 * 1.35
    ax.annotate('', xy=(x_center - box_w/2 - 0.1, step1_y),
                xytext=(x_center - 1.2, diamond_y),
                arrowprops=dict(arrowstyle='->', linewidth=2,
                                color='#10B981', mutation_scale=14,
                                connectionstyle='arc3,rad=0.4'))
    ax.text(x_center - 2.2, (diamond_y + step1_y) / 2, "Yes",
            fontsize=12, fontweight='bold', color='#10B981')

    # No arrow: to end
    end_x = x_center + 2.8
    ax.annotate('', xy=(end_x, diamond_y),
                xytext=(x_center + 1.2, diamond_y),
                arrowprops=dict(arrowstyle='->', linewidth=2,
                                color='#EF4444', mutation_scale=14))
    ax.text(end_x + 0.1, diamond_y + 0.15, "No", fontsize=12,
            fontweight='bold', color='#EF4444')
    box_end = FancyBboxPatch((end_x + 0.3, diamond_y - 0.35), 1.8, 0.7,
                             boxstyle="round,pad=0.1",
                             facecolor='#EF4444', alpha=0.15,
                             edgecolor='#EF4444', linewidth=2)
    ax.add_patch(box_end)
    ax.text(end_x + 1.2, diamond_y, "Record\nFinal State", ha='center', va='center',
            fontsize=10, fontweight='bold', color='#EF4444')

    ax.set_title('Gillespie Stochastic Simulation Algorithm (SSA)', fontsize=16,
                 fontweight='bold', color=DARK_HEX, pad=15)

    plt.tight_layout()
    plt.savefig(f"{FIGS}/gillespie_flowchart.png", dpi=200, bbox_inches='tight',
                facecolor='#F8FAFC')
    plt.close()
    print("  [OK] gillespie_flowchart.png")


# ============================================================
# FIGURE 3: four_conditions.png
# ============================================================
def make_four_conditions():
    fig, axes = plt.subplots(1, 4, figsize=(14, 4))
    fig.patch.set_facecolor('#F8FAFC')

    conditions = [
        ("A: Asymmetric\n(Native)", "Kd_s=2.4 nM\nKd_w=49 nM\nblock_s=85%, block_w=50%",
         ['#0D9488', '#F59E0B'], [0.85, 0.15]),
        ("B: Symmetric", "Kd=10.84 nM (both)\nblock=65.2% (both)\nGeometric mean",
         ['#3B82F6', '#3B82F6'], [0.5, 0.5]),
        ("C: Single-Site", "Kd_s=2.4 nM\nWeak site disabled\n(Kd_w -> infinity)",
         ['#7C3AED', '#D1D5DB'], [0.85, 0.0]),
        ("D: No Regulation", "k_on = 0\nAlways unbound\nMax expression",
         ['#D1D5DB', '#D1D5DB'], [0.0, 0.0]),
    ]

    for ax, (title, desc, colors, strengths) in zip(axes, conditions):
        ax.set_xlim(0, 10)
        ax.set_ylim(0, 10)
        ax.axis('off')

        # DNA line
        ax.plot([1, 9], [5.5, 5.5], color=GRAY_HEX, linewidth=4, solid_capstyle='round')

        # Strong site
        alpha_s = max(0.15, strengths[0])
        rect_s = FancyBboxPatch((2, 5.0), 2.5, 1.0,
                                boxstyle="round,pad=0.05",
                                facecolor=colors[0], alpha=alpha_s,
                                edgecolor=colors[0], linewidth=2)
        ax.add_patch(rect_s)
        ax.text(3.25, 5.5, "Strong", ha='center', va='center', fontsize=8,
                fontweight='bold', color='white' if alpha_s > 0.4 else DARK_HEX)

        # Weak site
        alpha_w = max(0.15, strengths[1])
        rect_w = FancyBboxPatch((5.5, 5.0), 2.5, 1.0,
                                boxstyle="round,pad=0.05",
                                facecolor=colors[1], alpha=alpha_w,
                                edgecolor=colors[1], linewidth=2)
        ax.add_patch(rect_w)
        ax.text(6.75, 5.5, "Weak", ha='center', va='center', fontsize=8,
                fontweight='bold', color='white' if alpha_w > 0.4 else DARK_HEX)

        ax.text(5, 8.5, title, ha='center', va='center', fontsize=12,
                fontweight='bold', color=DARK_HEX)
        ax.text(5, 2.5, desc, ha='center', va='center', fontsize=9,
                color=GRAY_HEX)

    fig.suptitle('Four Simulation Conditions', fontsize=16, fontweight='bold',
                 color=DARK_HEX, y=0.98)
    plt.tight_layout(rect=[0, 0, 1, 0.93])
    plt.savefig(f"{FIGS}/four_conditions.png", dpi=200, bbox_inches='tight',
                facecolor='#F8FAFC')
    plt.close()
    print("  [OK] four_conditions.png")


# ============================================================
# FIGURE 4: reaction_table.png
# ============================================================
def make_reaction_table():
    fig, ax = plt.subplots(figsize=(12, 7))
    ax.axis('off')
    fig.patch.set_facecolor('#F8FAFC')

    col_labels = ['#', 'Reaction', 'Transition', 'Propensity']
    rows = [
        ['0', 'Bind strong (from empty)',       '0 -> 1', r'$k_{on} \cdot [R_{free}]$'],
        ['1', 'Unbind strong',                   '1 -> 0', r'$k_{off,strong}$'],
        ['2', 'Bind weak (from empty)',          '0 -> 2', r'$k_{on} \cdot [R_{free}]$'],
        ['3', 'Unbind weak',                     '2 -> 0', r'$k_{off,weak}$'],
        ['4', 'Bind strong (weak occupied)',     '2 -> 3', r'$k_{on} \cdot [R_{free}]$'],
        ['5', 'Unbind strong (both occupied)',   '3 -> 2', r'$k_{off,strong}$'],
        ['6', 'Bind weak (strong occupied)',     '1 -> 3', r'$k_{on} \cdot [R_{free}]$'],
        ['7', 'Unbind weak (both occupied)',     '3 -> 1', r'$k_{off,weak}$'],
        ['8', 'Transcription',                    'mRNA += 1', r'$k_{txn}[op\_state]$'],
        ['9', 'Translation',                      'protein += 1', r'$k_{translation} \cdot mRNA$'],
        ['10', 'mRNA decay',                      'mRNA -= 1', r'$\gamma_{mRNA} \cdot mRNA$'],
        ['11', 'Protein decay',                   'protein -= 1', r'$\gamma_{protein} \cdot protein$'],
    ]

    # Color coding
    row_colors = []
    for i in range(len(rows)):
        if i < 8:
            row_colors.append(['#E0F2F1', '#E8F5E9', '#E0F2F1', '#E8F5E9'][i % 2])
        else:
            row_colors.append('#FFF3E0' if i % 2 == 0 else '#FFF8E1')

    table = ax.table(
        cellText=rows,
        colLabels=col_labels,
        loc='center',
        cellLoc='center',
    )
    table.auto_set_font_size(False)
    table.set_fontsize(11)
    table.scale(1, 1.6)

    # Style header
    for j in range(len(col_labels)):
        cell = table[0, j]
        cell.set_facecolor(TEAL_HEX)
        cell.set_text_props(color='white', fontweight='bold', fontsize=12)
        cell.set_edgecolor('white')

    # Style rows
    for i in range(len(rows)):
        for j in range(len(col_labels)):
            cell = table[i + 1, j]
            cell.set_facecolor(row_colors[i])
            cell.set_edgecolor('#E2E8F0')
            if i < 8:
                cell.set_text_props(color=DARK_HEX)
            else:
                cell.set_text_props(color='#9A3412', fontweight='bold')

    # Add section labels
    ax.text(0.02, 0.75, 'Operator\nTransitions\n(8)', transform=ax.transAxes,
            fontsize=11, fontweight='bold', color=TEAL_HEX, va='center',
            rotation=90)
    ax.text(0.02, 0.2, 'Gene\nExpression\n(4)', transform=ax.transAxes,
            fontsize=11, fontweight='bold', color='#9A3412', va='center',
            rotation=90)

    ax.set_title('All 12 Reactions in the Gillespie Model', fontsize=16,
                 fontweight='bold', color=DARK_HEX, pad=20)

    plt.tight_layout()
    plt.savefig(f"{FIGS}/reaction_table.png", dpi=200, bbox_inches='tight',
                facecolor='#F8FAFC')
    plt.close()
    print("  [OK] reaction_table.png")


# ============================================================
# Generate all figures
# ============================================================
print("Generating Phase 2 slide figures...")
make_operator_states()
make_gillespie_flowchart()
make_four_conditions()
make_reaction_table()
print("All figures generated.\n")


# ============================================================
# PPTX GENERATION
# ============================================================
print("Building Phase 2 PPTX...")

TEAL = RGBColor(0x0D, 0x94, 0x88)
DARK_TEAL = RGBColor(0x0A, 0x6E, 0x64)
WHITE = RGBColor(0xFF, 0xFF, 0xFF)
DARK = RGBColor(0x1E, 0x29, 0x3B)
GRAY = RGBColor(0x64, 0x74, 0x8B)
LIGHT_BG = RGBColor(0xF8, 0xFA, 0xFC)
ACCENT = RGBColor(0xDC, 0x26, 0x26)

prs = Presentation()
prs.slide_width = Inches(13.333)
prs.slide_height = Inches(7.5)


def add_blank_slide():
    layout = prs.slide_layouts[6]
    return prs.slides.add_slide(layout)


def add_bg(slide, color=LIGHT_BG):
    bg = slide.background
    fill = bg.fill
    fill.solid()
    fill.fore_color.rgb = color


def add_shape(slide, x, y, w, h, color):
    shape = slide.shapes.add_shape(MSO_SHAPE.RECTANGLE,
                                    Inches(x), Inches(y), Inches(w), Inches(h))
    shape.fill.solid()
    shape.fill.fore_color.rgb = color
    shape.line.fill.background()
    return shape


def add_text_box(slide, x, y, w, h, text, size=16, color=DARK, bold=False,
                 align=PP_ALIGN.LEFT, font_name='Calibri'):
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
        pPr = p._p.get_or_add_pPr()
        buChar = etree.SubElement(pPr, qn('a:buChar'))
        buChar.set('char', '\u2022')
    return txBox


def add_image(slide, path, x, y, w=None, h=None):
    if not os.path.exists(path):
        # Placeholder text if image missing
        add_text_box(slide, x, y, w or 5, h or 3,
                     f"[Image not found:\n{os.path.basename(path)}]",
                     size=12, color=GRAY, align=PP_ALIGN.CENTER)
        return
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

add_text_box(slide, 1, 1.6, 11, 1.5, "Phase 2: Gillespie Stochastic Simulation",
             size=44, color=WHITE, bold=True, align=PP_ALIGN.CENTER)
add_text_box(slide, 1, 3.2, 11, 0.8, "ENIGMA Project",
             size=28, color=TEAL, bold=True, align=PP_ALIGN.CENTER)
add_text_box(slide, 1, 4.1, 11, 0.8,
             "Simulating Single-Cell Gene Expression Noise",
             size=20, color=GRAY, align=PP_ALIGN.CENTER)
add_text_box(slide, 1, 5.2, 11, 0.5,
             "4-State Operator Model  |  50,000 Cells  |  12 Reactions  |  Numba JIT",
             size=13, color=GRAY, align=PP_ALIGN.CENTER)
add_text_box(slide, 1, 5.8, 11, 0.5,
             "Panagoda et al. 2024  |  Computational Biology",
             size=12, color=GRAY, align=PP_ALIGN.CENTER)
print("  Slide 1: Title")

# ══════════════════════════════════════════════════════════════
# SLIDE 2: Why Stochastic?
# ══════════════════════════════════════════════════════════════
slide = add_blank_slide()
add_bg(slide, LIGHT_BG)
teal_header_bar(slide)
add_text_box(slide, 0.8, 0.5, 11, 0.8, "Why Stochastic Simulation?",
             size=36, color=DARK, bold=True)

# Stat callout
add_shape(slide, 0.8, 1.6, 3.2, 1.8, TEAL)
add_text_box(slide, 0.8, 1.7, 3.2, 0.9, "~2 mRNA",
             size=40, color=WHITE, bold=True, align=PP_ALIGN.CENTER)
add_text_box(slide, 0.8, 2.55, 3.2, 0.6, "molecules at\nsteady state",
             size=14, color=WHITE, align=PP_ALIGN.CENTER)

add_bullets(slide, 4.5, 1.6, 8.0, 4.5, [
    "Gene expression is inherently random -- individual molecules bind and unbind stochastically",
    "With only ~2 mRNA molecules at steady state in TB, each binding/unbinding event creates measurable cell-to-cell variation",
    "Deterministic ODE models compute the average behavior, but cannot capture the noise that drives persistence",
    ("Key insight: ", "We need to simulate each molecule individually to measure how operator architecture shapes expression noise"),
    "The Gillespie algorithm (SSA) gives exact stochastic trajectories -- every reaction is sampled from the correct probability distribution",
], size=14)
print("  Slide 2: Why Stochastic?")

# ══════════════════════════════════════════════════════════════
# SLIDE 3: 4-State Operator Model
# ══════════════════════════════════════════════════════════════
slide = add_blank_slide()
add_bg(slide, LIGHT_BG)
teal_header_bar(slide)
add_text_box(slide, 0.8, 0.5, 11, 0.8, "The 4-State Operator Model",
             size=36, color=DARK, bold=True)
add_image(slide, f"{FIGS}/operator_states.png", 0.3, 1.3, w=7.5, h=5.5)

add_bullets(slide, 8.0, 1.5, 4.8, 5.0, [
    ("State 0: ", "Both sites empty -- full transcription (k_max = 0.15/min)"),
    ("State 1: ", "Strong site bound -- 85% blocked (k_txn = 0.0225/min)"),
    ("State 2: ", "Weak site bound -- 50% blocked (k_txn = 0.075/min)"),
    ("State 3: ", "Both bound -- 92.5% blocked (k_txn = 0.01125/min)"),
    "Each state has a different transcription rate, creating state-dependent bursting",
    "8 transitions between states driven by concentration-dependent binding and first-order unbinding",
], size=12.5)
print("  Slide 3: 4-State Operator")

# ══════════════════════════════════════════════════════════════
# SLIDE 4: All 12 Reactions
# ══════════════════════════════════════════════════════════════
slide = add_blank_slide()
add_bg(slide, LIGHT_BG)
teal_header_bar(slide)
add_text_box(slide, 0.8, 0.5, 11, 0.8, "All 12 Reactions",
             size=36, color=DARK, bold=True)
add_image(slide, f"{FIGS}/reaction_table.png", 0.5, 1.3, w=12.3, h=5.8)
print("  Slide 4: 12 Reactions")

# ══════════════════════════════════════════════════════════════
# SLIDE 5: Gillespie Algorithm
# ══════════════════════════════════════════════════════════════
slide = add_blank_slide()
add_bg(slide, LIGHT_BG)
teal_header_bar(slide)
add_text_box(slide, 0.8, 0.5, 11, 0.8, "The Gillespie Algorithm (SSA)",
             size=36, color=DARK, bold=True)
add_image(slide, f"{FIGS}/gillespie_flowchart.png", 0.3, 1.2, w=7.0, h=5.5)

add_bullets(slide, 7.8, 1.5, 5.0, 5.0, [
    ("Step 1: ", "Compute propensity a_i for each of the 12 reactions based on current state"),
    ("Step 2: ", "Sum all propensities to get a_total"),
    ("Step 3: ", "Draw waiting time tau from Exp(a_total) -- time until next reaction"),
    ("Step 4: ", "Select which reaction fires with probability a_i / a_total"),
    ("Step 5: ", "Update molecule counts (change operator state, or +/-1 mRNA/protein)"),
    ("Repeat: ", "Continue until t reaches t_max = 30,000 min (~21 days)"),
], size=12)
print("  Slide 5: Gillespie Algorithm")

# ══════════════════════════════════════════════════════════════
# SLIDE 6: Key Parameters
# ══════════════════════════════════════════════════════════════
slide = add_blank_slide()
add_bg(slide, LIGHT_BG)
teal_header_bar(slide)
add_text_box(slide, 0.8, 0.5, 11, 0.8, "Key Parameters",
             size=36, color=DARK, bold=True)

params_data = [
    ("Parameter", "Value", "Source"),
    ("k_max (max transcription)", "0.15 /min", "Mtb transcriptomics"),
    ("k_translation", "0.5 protein/mRNA/min", "Taniguchi 2010"),
    ("mRNA half-life", "9.5 min", "Rustad 2013"),
    ("Protein half-life", "1,500 min (~25 hr)", "Taniguchi 2010"),
    ("k_on (binding rate)", "0.0167 nM^-1 min^-1", "Stormo & Zhao 2010"),
    ("Kd_strong", "2.4 nM", "Panagoda 2024"),
    ("Kd_weak", "49 nM", "Panagoda 2024"),
    ("block_strong / block_weak", "85% / 50%", "Estimated"),
    ("t_max", "30,000 min (500 hr)", "--"),
    ("t_burn_in", "15,000 min", "10x protein t_1/2"),
    ("n_cells", "50,000", "--"),
]

# Build a table
left = Inches(0.8)
top = Inches(1.5)
width = Inches(11.5)
height = Inches(5.5)
table_shape = slide.shapes.add_table(len(params_data), 3, left, top, width, height)
table = table_shape.table

# Set column widths
table.columns[0].width = Inches(4.5)
table.columns[1].width = Inches(3.5)
table.columns[2].width = Inches(3.5)

for row_idx, row_data in enumerate(params_data):
    for col_idx, cell_text in enumerate(row_data):
        cell = table.cell(row_idx, col_idx)
        cell.text = cell_text
        for paragraph in cell.text_frame.paragraphs:
            paragraph.font.size = Pt(13)
            paragraph.font.name = 'Calibri'
            if row_idx == 0:
                paragraph.font.bold = True
                paragraph.font.color.rgb = WHITE
            else:
                paragraph.font.color.rgb = DARK
        if row_idx == 0:
            cell.fill.solid()
            cell.fill.fore_color.rgb = TEAL
        else:
            cell.fill.solid()
            cell.fill.fore_color.rgb = RGBColor(0xF1, 0xF5, 0xF9) if row_idx % 2 == 0 else WHITE

print("  Slide 6: Key Parameters")

# ══════════════════════════════════════════════════════════════
# SLIDE 7: Four Simulation Conditions
# ══════════════════════════════════════════════════════════════
slide = add_blank_slide()
add_bg(slide, LIGHT_BG)
teal_header_bar(slide)
add_text_box(slide, 0.8, 0.5, 11, 0.8, "Four Simulation Conditions",
             size=36, color=DARK, bold=True)
add_image(slide, f"{FIGS}/four_conditions.png", 0.3, 1.3, w=12.5, h=3.8)

# Condition cards
cards = [
    ("A: Asymmetric", "Native Mce3R\nKd = 2.4/49 nM", TEAL),
    ("B: Symmetric", "Geometric mean\nKd = 10.84/10.84 nM", RGBColor(0x3B, 0x82, 0xF6)),
    ("C: Single-Site", "Strong site only\nWeak disabled", RGBColor(0x7C, 0x3A, 0xED)),
    ("D: No Regulation", "k_on = 0\nAlways unbound", ACCENT),
]
for i, (title, desc, clr) in enumerate(cards):
    x = 0.8 + i * 3.15
    add_shape(slide, x, 5.3, 2.9, 1.8, WHITE)
    add_shape(slide, x, 5.3, 2.9, 0.08, clr)
    add_text_box(slide, x + 0.15, 5.5, 2.6, 0.5, title, size=14, color=clr, bold=True)
    add_text_box(slide, x + 0.15, 6.05, 2.6, 0.8, desc, size=11, color=GRAY)

print("  Slide 7: Four Conditions")

# ══════════════════════════════════════════════════════════════
# SLIDE 8: Protein Distributions
# ══════════════════════════════════════════════════════════════
slide = add_blank_slide()
add_bg(slide, LIGHT_BG)
teal_header_bar(slide)
add_text_box(slide, 0.8, 0.5, 11, 0.8, "What 50,000 Cells Looks Like",
             size=36, color=DARK, bold=True)

fig2_path = f"{EXISTING_FIGS}/fig2_distributions.png"
add_image(slide, fig2_path, 0.5, 1.3, w=8.0, h=5.5)

add_bullets(slide, 8.8, 1.5, 4.0, 5.0, [
    "Each histogram shows the distribution of protein levels across 50,000 simulated cells",
    ("Condition A: ", "Mean = 197, lowest due to strong repression"),
    ("Condition B: ", "Mean = 293, higher because symmetric sites are weaker individually"),
    ("Condition D: ", "Mean = 2,224, the unregulated maximum"),
    "The WIDTH of each distribution (CV) is the key noise metric",
], size=12)
print("  Slide 8: Distributions")

# ══════════════════════════════════════════════════════════════
# SLIDE 9: Key Result - CV Comparison
# ══════════════════════════════════════════════════════════════
slide = add_blank_slide()
add_bg(slide, DARK)
add_shape(slide, 0, 0, 13.333, 0.2, TEAL)
add_shape(slide, 0, 7.3, 13.333, 0.2, TEAL)

add_text_box(slide, 1, 0.6, 11, 0.8, "The Key Result",
             size=36, color=WHITE, bold=True, align=PP_ALIGN.CENTER)

# Big CV comparison
add_shape(slide, 1.5, 1.8, 4.5, 3.0, TEAL)
add_text_box(slide, 1.5, 2.0, 4.5, 0.6, "Asymmetric (Native)",
             size=16, color=WHITE, bold=True, align=PP_ALIGN.CENTER)
add_text_box(slide, 1.5, 2.7, 4.5, 1.2, "CV = 0.187",
             size=48, color=WHITE, bold=True, align=PP_ALIGN.CENTER)
add_text_box(slide, 1.5, 4.0, 4.5, 0.5, "Mean = 197 proteins",
             size=14, color=RGBColor(0xA7, 0xF3, 0xD0), align=PP_ALIGN.CENTER)

add_shape(slide, 7.3, 1.8, 4.5, 3.0, RGBColor(0x3B, 0x82, 0xF6))
add_text_box(slide, 7.3, 2.0, 4.5, 0.6, "Symmetric (Control)",
             size=16, color=WHITE, bold=True, align=PP_ALIGN.CENTER)
add_text_box(slide, 7.3, 2.7, 4.5, 1.2, "CV = 0.157",
             size=48, color=WHITE, bold=True, align=PP_ALIGN.CENTER)
add_text_box(slide, 7.3, 4.0, 4.5, 0.5, "Mean = 293 proteins",
             size=14, color=RGBColor(0xBF, 0xDB, 0xFE), align=PP_ALIGN.CENTER)

# Delta
add_shape(slide, 4.0, 5.2, 5.3, 1.5, RGBColor(0xDC, 0x26, 0x26))
add_text_box(slide, 4.0, 5.3, 5.3, 0.5,
             "Asymmetric CV is ~19% higher",
             size=22, color=WHITE, bold=True, align=PP_ALIGN.CENTER)
add_text_box(slide, 4.0, 5.9, 5.3, 0.6,
             "The native operator produces more noise than an equivalent symmetric design",
             size=13, color=RGBColor(0xFE, 0xCA, 0xCA), align=PP_ALIGN.CENTER)

print("  Slide 9: Key Result")

# ══════════════════════════════════════════════════════════════
# SLIDE 10: Why Asymmetry Creates More Noise
# ══════════════════════════════════════════════════════════════
slide = add_blank_slide()
add_bg(slide, LIGHT_BG)
teal_header_bar(slide)
add_text_box(slide, 0.8, 0.5, 11, 0.8, "Why Asymmetry Creates More Noise",
             size=36, color=DARK, bold=True)

add_bullets(slide, 0.8, 1.5, 11.5, 5.0, [
    ("The weak site flickers: ", "Kd = 49 nM means rapid binding/unbinding. The weak site is occupied only ~25% of the time, constantly switching between states."),
    ("The strong site anchors: ", "Kd = 2.4 nM means the strong site stays bound most of the time, holding the operator in a partially repressed state."),
    ("Transcriptional bursting: ", "When the weak site briefly unbinds, transcription rate jumps from 0.01125/min (both bound) to 0.0225/min (strong only) -- a 2x burst."),
    ("Wider distribution: ", "This flickering creates cells with very different expression histories. Some cells catch more bursts, others fewer, widening the protein distribution."),
    ("Symmetric operators are calmer: ", "With two equal-affinity sites (Kd = 10.84 nM each), both sites flicker at the same rate. The transitions between states are more balanced, producing a narrower distribution."),
    ("The biological implication: ", "Higher noise means a larger tail of the distribution -- more cells reach the extreme expression levels that may trigger the persister phenotype."),
], size=13.5, spacing=Pt(12))
print("  Slide 10: Why Asymmetry")

# ══════════════════════════════════════════════════════════════
# SLIDE 11: Numba JIT Acceleration
# ══════════════════════════════════════════════════════════════
slide = add_blank_slide()
add_bg(slide, LIGHT_BG)
teal_header_bar(slide)
add_text_box(slide, 0.8, 0.5, 11, 0.8, "Numba JIT Acceleration",
             size=36, color=DARK, bold=True)

# Performance stat callout
add_shape(slide, 0.8, 1.6, 3.5, 1.8, TEAL)
add_text_box(slide, 0.8, 1.7, 3.5, 0.9, "~10 min",
             size=40, color=WHITE, bold=True, align=PP_ALIGN.CENTER)
add_text_box(slide, 0.8, 2.55, 3.5, 0.6, "for 200,000 cells\n(all 4 conditions)",
             size=13, color=WHITE, align=PP_ALIGN.CENTER)

add_bullets(slide, 4.8, 1.6, 7.8, 5.0, [
    "Pure Python Gillespie would take hours to days for 50,000 cells per condition",
    ("@numba.njit: ", "Compiles the inner simulation loop to optimized machine code at first call"),
    "The simulate_cell() function runs ~100x faster than interpreted Python",
    "First call takes ~30 seconds for JIT compilation; subsequent calls are near-instant",
    ("Performance breakdown: ", "Condition A: 42s, Condition B: 36s, Condition C: 390s (weak site flickers constantly), Condition D: 83s"),
    "Total Phase 2 runtime: ~593 seconds (~10 minutes) for all conditions + asymmetry sweep",
], size=13)
print("  Slide 11: Numba JIT")

# ══════════════════════════════════════════════════════════════
# SLIDE 12: Condition D Baseline
# ══════════════════════════════════════════════════════════════
slide = add_blank_slide()
add_bg(slide, LIGHT_BG)
teal_header_bar(slide)
add_text_box(slide, 0.8, 0.5, 11, 0.8, "Condition D: The Unregulated Baseline",
             size=36, color=DARK, bold=True)

# Stat cards
stats_d = [
    ("2,224", "Mean Protein", TEAL),
    ("0.059", "CV (Noise)", RGBColor(0x3B, 0x82, 0xF6)),
    ("7.77", "Fano Factor", RGBColor(0x7C, 0x3A, 0xED)),
    ("11.3x", "Fold Repression\n(D/A)", ACCENT),
]
for i, (val, label, clr) in enumerate(stats_d):
    x = 0.8 + i * 3.15
    add_shape(slide, x, 1.6, 2.9, 2.0, clr)
    add_text_box(slide, x, 1.7, 2.9, 1.0, val,
                 size=36, color=WHITE, bold=True, align=PP_ALIGN.CENTER)
    add_text_box(slide, x, 2.7, 2.9, 0.7, label,
                 size=13, color=WHITE, align=PP_ALIGN.CENTER)

add_bullets(slide, 0.8, 4.0, 11.5, 3.0, [
    "With k_on = 0, no repressor binds -- the operator is always in State 0 (fully active)",
    "Mean protein = 2,224 serves as the maximum expression reference for all conditions",
    ("Fold-change validation: ", "Model predicts 11.3x repression (D/A), compared to experimental ~8.5x from Santangelo 2009"),
    ("Fano factor = 7.77: ", "Close to theoretical prediction of 1 + burst_size = 1 + 6.85 = 7.85 for unregulated two-stage model"),
], size=13.5)
print("  Slide 12: Condition D")

# ══════════════════════════════════════════════════════════════
# SLIDE 13: Parameter Sweep
# ══════════════════════════════════════════════════════════════
slide = add_blank_slide()
add_bg(slide, LIGHT_BG)
teal_header_bar(slide)
add_text_box(slide, 0.8, 0.5, 11, 0.8, "Asymmetry Sweep (Condition E)",
             size=36, color=DARK, bold=True)

fig3_path = f"{EXISTING_FIGS}/fig3_asymmetry_sweep.png"
add_image(slide, fig3_path, 0.3, 1.3, w=8.0, h=5.5)

add_bullets(slide, 8.5, 1.5, 4.3, 5.0, [
    "Sweep Kd_weak/Kd_strong ratio from 1 (symmetric) to 50 (extreme asymmetry)",
    "Geometric mean held constant at 10.84 nM to isolate the effect of asymmetry",
    "CV increases from 0.158 (ratio=1) to 0.191 (ratio=20)",
    "5,000 cells simulated per ratio point (10 ratios total)",
    ("Native Mce3R ratio = 20.4: ", "sits near the peak of the noise curve"),
], size=12)
print("  Slide 13: Asymmetry Sweep")

# ══════════════════════════════════════════════════════════════
# SLIDE 14: Single-Cell Traces
# ══════════════════════════════════════════════════════════════
slide = add_blank_slide()
add_bg(slide, LIGHT_BG)
teal_header_bar(slide)
add_text_box(slide, 0.8, 0.5, 11, 0.8, "Single-Cell Traces",
             size=36, color=DARK, bold=True)

fig6_path = f"{EXISTING_FIGS}/fig6_single_cell_traces.png"
add_image(slide, fig6_path, 0.3, 1.3, w=8.0, h=5.5)

add_bullets(slide, 8.5, 1.5, 4.3, 5.0, [
    "Individual cell protein levels fluctuate over time as molecules are produced and degraded",
    "Each line is one cell's trajectory recorded every 10 minutes after burn-in",
    "Stochastic bursting is visible as sharp peaks when the operator briefly opens",
    "Cells in the same condition can have very different instantaneous protein levels",
    "This cell-to-cell variation IS the noise that we quantify with CV",
], size=12)
print("  Slide 14: Single-Cell Traces")

# ══════════════════════════════════════════════════════════════
# SLIDE 15: Sanity Checks
# ══════════════════════════════════════════════════════════════
slide = add_blank_slide()
add_bg(slide, LIGHT_BG)
teal_header_bar(slide)
add_text_box(slide, 0.8, 0.5, 11, 0.8, "Sanity Checks",
             size=36, color=DARK, bold=True)

checks = [
    ("Mean protein in expected range: ", "Condition A = 197 (expected ~200 for regulated TB gene)"),
    ("CV > 0 for all conditions: ", "Confirmed -- all simulations show measurable noise"),
    ("Fano factor > 1 (super-Poissonian): ", "All conditions show Fano 6.87-7.81, consistent with transcriptional bursting"),
    ("Condition D has highest mean: ", "2,224 >> 342 >> 293 >> 197 (D > C > B > A), confirming repression hierarchy"),
    ("Fano(D) ~ 1 + burst_size: ", "7.77 vs expected 7.85 -- validates the two-stage gene expression model"),
    ("Reproducibility: ", "All conditions use deterministic seeds (master_seed=42 + condition offset)"),
    ("All NPZ files saved: ", "condition_A.npz, condition_B.npz, condition_C.npz, condition_D.npz, condition_E_sweep.npz"),
    ("Numba compilation: ", "JIT compiles successfully on first call, no type errors"),
]

add_bullets(slide, 0.8, 1.5, 11.5, 5.5, checks, size=13, spacing=Pt(10))
print("  Slide 15: Sanity Checks")

# ══════════════════════════════════════════════════════════════
# SLIDE 16: Output Files
# ══════════════════════════════════════════════════════════════
slide = add_blank_slide()
add_bg(slide, LIGHT_BG)
teal_header_bar(slide)
add_text_box(slide, 0.8, 0.5, 11, 0.8, "Output Files",
             size=36, color=DARK, bold=True)

files_data = [
    ("File", "Contents", "Size"),
    ("condition_A.npz", "50,000 protein/mRNA/op_state arrays + 20 traces", "Asymmetric (native)"),
    ("condition_B.npz", "50,000 protein/mRNA/op_state arrays + 20 traces", "Symmetric control"),
    ("condition_C.npz", "50,000 protein/mRNA/op_state arrays + 20 traces", "Single-site control"),
    ("condition_D.npz", "50,000 protein/mRNA/op_state arrays + 20 traces", "Unregulated baseline"),
    ("condition_E_sweep.npz", "5,000 cells x 10 ratios", "Asymmetry sweep"),
    ("phase2_summary.json", "Mean, CV, Fano for all conditions", "Summary statistics"),
]

left = Inches(0.8)
top = Inches(1.5)
width = Inches(11.5)
height = Inches(4.0)
table_shape = slide.shapes.add_table(len(files_data), 3, left, top, width, height)
table = table_shape.table

table.columns[0].width = Inches(3.5)
table.columns[1].width = Inches(5.0)
table.columns[2].width = Inches(3.0)

for row_idx, row_data in enumerate(files_data):
    for col_idx, cell_text in enumerate(row_data):
        cell = table.cell(row_idx, col_idx)
        cell.text = cell_text
        for paragraph in cell.text_frame.paragraphs:
            paragraph.font.size = Pt(13)
            paragraph.font.name = 'Calibri'
            if row_idx == 0:
                paragraph.font.bold = True
                paragraph.font.color.rgb = WHITE
            else:
                paragraph.font.color.rgb = DARK
        if row_idx == 0:
            cell.fill.solid()
            cell.fill.fore_color.rgb = TEAL
        else:
            cell.fill.solid()
            cell.fill.fore_color.rgb = RGBColor(0xF1, 0xF5, 0xF9) if row_idx % 2 == 0 else WHITE

add_text_box(slide, 0.8, 5.8, 11.5, 0.8,
             "All output stored in results/phase2/. NPZ files use numpy compressed format for efficient storage of large arrays.",
             size=12, color=GRAY)
print("  Slide 16: Output Files")

# ══════════════════════════════════════════════════════════════
# SLIDE 17: What Phase 2 Feeds Into
# ══════════════════════════════════════════════════════════════
slide = add_blank_slide()
add_bg(slide, LIGHT_BG)
teal_header_bar(slide)
add_text_box(slide, 0.8, 0.5, 11, 0.8, "What Phase 2 Feeds Into",
             size=36, color=DARK, bold=True)

# Pipeline cards
pipeline = [
    ("Phase 3", "Statistical Analysis", "Bootstrap CIs, Mann-Whitney tests,\nGMM fitting, sensitivity analysis", TEAL),
    ("Phase 4", "Figure Generation", "7 publication-quality figures\nat 300 dpi (matplotlib + seaborn)", RGBColor(0x3B, 0x82, 0xF6)),
    ("Phase 5", "Calibration", "Mean protein from Phase 2 becomes\nthe calibration target for\nthermodynamic model", RGBColor(0x7C, 0x3A, 0xED)),
]

for i, (phase, title, desc, clr) in enumerate(pipeline):
    x = 0.8 + i * 4.15
    add_shape(slide, x, 1.6, 3.8, 3.5, WHITE)
    add_shape(slide, x, 1.6, 3.8, 0.1, clr)
    add_text_box(slide, x + 0.2, 1.9, 3.4, 0.5, phase, size=20, color=clr, bold=True)
    add_text_box(slide, x + 0.2, 2.5, 3.4, 0.5, title, size=16, color=DARK, bold=True)
    add_text_box(slide, x + 0.2, 3.1, 3.4, 1.5, desc, size=12, color=GRAY)
    # Arrow between cards
    if i < len(pipeline) - 1:
        arrow_x = x + 3.8
        add_shape(slide, arrow_x, 3.0, 0.35, 0.08, clr)

add_bullets(slide, 0.8, 5.5, 11.5, 1.5, [
    ("Phase 2 is the engine: ", "It produces the raw stochastic data that all downstream analyses depend on"),
    "The 50,000-cell protein arrays from each condition feed directly into Phase 3 statistical tests and Phase 4 visualizations",
], size=13)
print("  Slide 17: Feeds Into")

# ══════════════════════════════════════════════════════════════
# SLIDE 18: Summary
# ══════════════════════════════════════════════════════════════
slide = add_blank_slide()
add_bg(slide, DARK)
add_shape(slide, 0, 0, 13.333, 0.2, TEAL)
add_shape(slide, 0, 7.3, 13.333, 0.2, TEAL)

add_text_box(slide, 1, 0.8, 11, 0.8, "Summary",
             size=36, color=WHITE, bold=True, align=PP_ALIGN.CENTER)

add_text_box(slide, 1.5, 2.0, 10, 1.5,
             "Phase 2 simulated 50,000 individual TB cells for each of 4 operator "
             "architectures using exact stochastic simulation (Gillespie SSA), "
             "demonstrating that the native asymmetric operator produces ~19% "
             "higher expression noise (CV) than an equivalent symmetric design.",
             size=18, color=WHITE, align=PP_ALIGN.CENTER)

# Summary stats row
summary_items = [
    ("50,000", "cells per\ncondition"),
    ("12", "reactions\nmodeled"),
    ("4", "operator\narchitectures"),
    ("19%", "more noise in\nasymmetric"),
    ("~10 min", "total\nruntime"),
]
for i, (val, label) in enumerate(summary_items):
    x = 1.0 + i * 2.4
    add_shape(slide, x, 4.0, 2.1, 2.0, TEAL)
    add_text_box(slide, x, 4.1, 2.1, 0.8, val,
                 size=30, color=WHITE, bold=True, align=PP_ALIGN.CENTER)
    add_text_box(slide, x, 4.9, 2.1, 0.8, label,
                 size=12, color=RGBColor(0xA7, 0xF3, 0xD0), align=PP_ALIGN.CENTER)

add_text_box(slide, 1, 6.4, 11, 0.5,
             "ENIGMA Project  |  Mce3R Stochastic Simulation  |  Phase 2 Complete",
             size=12, color=GRAY, align=PP_ALIGN.CENTER)
print("  Slide 18: Summary")

# ══════════════════════════════════════════════════════════════
# SAVE
# ══════════════════════════════════════════════════════════════
prs.save(OUT_PPTX)
print(f"\nPresentation saved to: {OUT_PPTX}")
print(f"Total slides: {len(prs.slides)}")
