#!/usr/bin/env python3
"""
Create a publication-quality visual abstract summarizing the ENIGMA project
methods pipeline as a flowchart-style infographic.
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import os

# ── Output path ──────────────────────────────────────────────────────────────
OUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
os.makedirs(OUT_DIR, exist_ok=True)
OUT_PATH = os.path.join(OUT_DIR, "visual_abstract.png")

# ── Colour palette ───────────────────────────────────────────────────────────
BLUE   = "#1E3A5F"
TEAL   = "#0D9488"
CORAL  = "#B91C1C"
ARROW  = "#374151"
DARK   = "#1E293B"
WHITE  = "#FFFFFF"

# ── Figure setup ─────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(16, 10), dpi=300)
fig.patch.set_facecolor(WHITE)
ax.set_facecolor(WHITE)
ax.set_xlim(0, 16)
ax.set_ylim(0, 10)
ax.axis('off')

# ── Helper: draw a rounded box with centred text ─────────────────────────────
def draw_box(x, y, w, h, colour, text, fontsize=8):
    box = FancyBboxPatch(
        (x, y), w, h,
        boxstyle="round,pad=0.15",
        facecolor=colour, edgecolor='none',
        zorder=2,
    )
    ax.add_patch(box)
    ax.text(
        x + w / 2, y + h / 2, text,
        ha='center', va='center',
        fontsize=fontsize, color=WHITE,
        fontfamily='sans-serif',
        linespacing=1.3,
        zorder=3,
    )
    return box


# ── Helper: horizontal arrow with label ──────────────────────────────────────
def h_arrow(x_start, x_end, y, label=""):
    ax.annotate(
        "", xy=(x_end, y), xytext=(x_start, y),
        arrowprops=dict(
            arrowstyle="->,head_width=0.25,head_length=0.15",
            color=ARROW, lw=1.5,
        ),
        zorder=1,
    )
    if label:
        mid = (x_start + x_end) / 2
        ax.text(
            mid, y + 0.22, label,
            ha='center', va='bottom',
            fontsize=9, color=ARROW,
            fontstyle='italic',
            fontweight='bold',
            zorder=3,
        )


# ── Helper: vertical dashed arrow with label ─────────────────────────────────
def v_dashed_arrow(x, y_start, y_end, label=""):
    ax.annotate(
        "", xy=(x, y_end), xytext=(x, y_start),
        arrowprops=dict(
            arrowstyle="->,head_width=0.25,head_length=0.15",
            color=ARROW, lw=1.2,
            linestyle='dashed',
        ),
        zorder=1,
    )
    if label:
        mid_y = (y_start + y_end) / 2
        ax.text(
            x + 0.15, mid_y, label,
            ha='left', va='center',
            fontsize=8, color=ARROW,
            fontstyle='italic',
            rotation=0,
            zorder=3,
        )


# ── Layout constants ─────────────────────────────────────────────────────────
BOX_W = 2.6
BOX_H = 1.6
GAP   = 0.55          # horizontal gap between boxes (arrow space)
LABEL_X = 0.15        # x for row labels
BOX_X0  = 2.4         # x for first box in each row

rows = [
    {"y": 7.5, "colour": BLUE,  "label": "Aim 1\nDiscover & Validate\nBinding Sites"},
    {"y": 4.5, "colour": TEAL,  "label": "Aim 2\nModel Operator\nRepression"},
    {"y": 1.5, "colour": CORAL, "label": "Aim 3\nQuantify Noise\n& Persistence"},
]

box_texts = [
    # Row 1 (Aim 1)
    [
        "3 Mycobacterial\nGenomes\n(H37Rv, M. bovis,\nM. marinum)",
        "De Novo Motif\nDiscovery\n(ZOOPS, both strands,\n6\u2013110 bp width)",
        "FIMO Genome\nScan\n(4.4 M bp,\np < 10\u207B\u2074)",
        "Validated Operator\n#1 hit: p = 1.1\u00d710\u207B\u00b9\u00b9\nKd_strong = 2.4 nM\nKd_weak = 49 nM",
    ],
    # Row 2 (Aim 2)
    [
        "4-State Operator\nModel\n(U, S, W, D states,\n12 reactions)",
        "Stochastic\nSimulation\n(50,000 cells \u00d7\n30,000 min each)",
        "Thermodynamic\nCalibration\n(\u0394G_strong = \u221212.2\n\u0394G_weak = \u221210.4 kcal/mol)",
        "Cooperativity\nInference\n(\u03c9 = 1.08\n95% CI: 0.20\u20135.89)",
    ],
    # Row 3 (Aim 3)
    [
        "4 Architectures \u00d7\n4 Environments\n(24 conditions,\ntwo-species model)",
        "CV(asym) = 0.187\nCV(sym) = 0.157\n(p < 0.001, robust\nacross full posterior)",
        "Persister\nFractions\n(3.1\u00d7 enrichment\nbaseline, 29.4\u00d7 acid)",
        "MI Trade-off\nAsym: 0.28 bits\nSym: 0.44 bits\n(noise vs info)",
    ],
]

arrow_labels = [
    ["MEME", "PWM", "1,442 sites"],
    ["Gillespie SSA", "Berg\u2013von Hippel", "MCMC"],
    ["Statistics", "Threshold", "Information"],
]

# ── Draw rows ────────────────────────────────────────────────────────────────
for i, row in enumerate(rows):
    y = row["y"]
    c = row["colour"]

    # Row label
    ax.text(
        LABEL_X, y + BOX_H / 2, row["label"],
        ha='left', va='center',
        fontsize=11, fontweight='bold',
        color=c, linespacing=1.35,
        zorder=3,
    )

    # Boxes and arrows
    for j in range(4):
        bx = BOX_X0 + j * (BOX_W + GAP)
        draw_box(bx, y, BOX_W, BOX_H, c, box_texts[i][j], fontsize=8)

        if j < 3:
            arr_x0 = bx + BOX_W
            arr_x1 = bx + BOX_W + GAP
            h_arrow(arr_x0, arr_x1, y + BOX_H / 2, arrow_labels[i][j])

# ── Vertical dashed arrows between rows ──────────────────────────────────────
# Row 1 Box 4 centre-x  -->  Row 2 Box 1 centre-x
r1_box4_cx = BOX_X0 + 3 * (BOX_W + GAP) + BOX_W / 2
r2_box1_cx = BOX_X0 + BOX_W / 2

# From bottom of Row1 Box4 to top of Row2 Box1
v_dashed_arrow(
    x=(r1_box4_cx + r2_box1_cx) / 2,
    y_start=rows[0]["y"],
    y_end=rows[1]["y"] + BOX_H,
    label="Kd values feed model",
)

# Row 2 Box 4 --> Row 3 Box 1
r2_box4_cx = BOX_X0 + 3 * (BOX_W + GAP) + BOX_W / 2
r3_box1_cx = BOX_X0 + BOX_W / 2

v_dashed_arrow(
    x=(r2_box4_cx + r3_box1_cx) / 2,
    y_start=rows[1]["y"],
    y_end=rows[2]["y"] + BOX_H,
    label="Calibrated parameters",
)

# ── Bottom conclusion bar ────────────────────────────────────────────────────
conc_y = 0.15
conc_h = 0.95
conc_x = 0.3
conc_w = 15.4
conc_box = FancyBboxPatch(
    (conc_x, conc_y), conc_w, conc_h,
    boxstyle="round,pad=0.2",
    facecolor=DARK, edgecolor='none',
    zorder=2,
)
ax.add_patch(conc_box)
ax.text(
    conc_x + conc_w / 2, conc_y + conc_h / 2,
    "CONCLUSION: Operator binding asymmetry (20.4\u00d7 Kd ratio) may generate "
    "persistence-enabling expression noise \u2014 a potential anti-noise therapeutic target",
    ha='center', va='center',
    fontsize=10, color=WHITE,
    fontweight='bold',
    zorder=3,
)

# ── Title and subtitle ───────────────────────────────────────────────────────
ax.text(
    8, 9.7,
    "ENIGMA: Visual Methods Abstract",
    ha='center', va='center',
    fontsize=18, fontweight='bold',
    color=DARK,
    zorder=3,
)
ax.text(
    8, 9.3,
    "Stochastic Noise from Asymmetric Operator Architecture in Mce3R",
    ha='center', va='center',
    fontsize=12, fontstyle='italic',
    color=ARROW,
    zorder=3,
)

# ── Save ─────────────────────────────────────────────────────────────────────
fig.tight_layout()
fig.savefig(OUT_PATH, dpi=300, bbox_inches='tight', facecolor=WHITE)
plt.close(fig)
print(f"Saved visual abstract to {OUT_PATH}")
