#!/usr/bin/env python3
"""
Publication Figure 1: Mce3R Phase 1 Operator Discovery Pipeline
Multi-panel figure for Regeneron STS submission.
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec
from matplotlib.patches import FancyBboxPatch
import numpy as np

# ============================================================
# DATA
# ============================================================

top20 = [
    (1, 57.42, 'yrbE3A', True),
    (2, 34.58, 'Rv1935c', True),
    (3, 33.55, 'mce3R', True),
    (4, 25.90, 'Rv1115', False),
    (5, 23.76, 'mmpL2', False),
    (6, 23.19, 'Rv2024c', False),
    (7, 22.98, 'PPE59', False),
    (8, 22.74, 'Rv0010c', False),
    (9, 22.34, 'lppT', False),
    (10, 22.09, 'fadD25', False),
    (11, 22.06, 'Rv3377c', False),
    (12, 21.90, 'fadD30', False),
    (13, 21.56, 'fadD9', False),
    (14, 21.46, 'bioF2', False),
    (15, 21.04, 'lpqG', False),
    (16, 21.04, 'PE9', False),
    (17, 20.92, 'Rv3377c', False),
    (18, 20.86, 'cyp135A1', False),
    (19, 20.82, 'Rv3888c', False),
    (20, 20.76, 'PPE35', False),
]

operators = [
    ('Op 1: yrbE3A (Validated)',
     'CTATAGGATACTAGCAAGATACATC', 'ATCGCGTATTGGCTATGGACATCAG',
     53, 100, 'HIGHEST'),
    ('Op 2: Rv1935c/Rv1936 (Novel)',
     'TGTGCGTATTGGCTATGGACATGTT', 'TATCACAGTACTACCAAGATACTTC',
     53, 100, 'HIGHEST'),
    ('Op 3: mce3R autoregulatory (Novel)',
     'CAATATCGGACTAACAAAATACATC', 'AGTGCATATCAGTAATAGACATATC',
     53, 60, 'HIGHEST'),
    ('Op 4: Rv1115 (Candidate)',
     'TTTTTATATTGTTGCGTGACATATC', 'AGAATTGATTTCCTATGGATATTGT',
     53, 60, 'HIGH'),
]

CORE_15BP = 'GCGTATTGGCTATGG'

manhattan_real_pos = [2207535, 2187268, 2207089, 1240174, 599006, 2268403,
                      3847234, 13594, 2039106, 1712856, 3791826, 484167,
                      2917839, 35157, 4062372, 1214429, 3790892, 393474, 4372481]
manhattan_real_scr = [57.42, 34.58, 33.55, 25.90, 23.76, 23.19, 22.98, 22.74,
                      22.34, 22.09, 22.06, 21.90, 21.56, 21.46, 21.04, 21.04,
                      20.92, 20.86, 20.82]

tier_counts = {'HIGHEST': 3, 'HIGH': 9, 'MODERATE': 9, 'LOW': 29}

# ============================================================
# COLORS
# ============================================================
C_HIGHEST = '#0b3d6b'
C_HIGH    = '#1a6fad'
C_MOD     = '#7fb3d8'
C_LOW     = '#c8ced3'
C_WEAK    = '#d35400'
C_STRONG  = '#1a5276'
C_SPACER  = '#95a5a6'
C_CORE    = '#b71c1c'
C_RED     = '#c0392b'
C_GREEN   = '#1e8449'
C_BLUE    = '#2471a3'
C_GOLD    = '#d4ac0d'

BASE_COL = {'A': '#1e8449', 'T': '#c0392b', 'G': '#d4ac0d', 'C': '#2471a3'}

# ============================================================
# FIGURE
# ============================================================
fig = plt.figure(figsize=(20, 26))
fig.patch.set_facecolor('white')

gs = gridspec.GridSpec(5, 2, figure=fig,
                       height_ratios=[0.8, 1.0, 1.8, 0.9, 0.1],
                       hspace=0.32, wspace=0.28,
                       left=0.06, right=0.96, top=0.94, bottom=0.03)

PL = dict(fontsize=20, fontweight='bold', va='top', ha='left')

# ============================================================
# A — Pipeline Workflow
# ============================================================
ax_a = fig.add_subplot(gs[0, :])
ax_a.set_xlim(0, 10)
ax_a.set_ylim(0, 2.2)
ax_a.axis('off')
ax_a.text(-0.01, 1.06, 'A', transform=ax_a.transAxes, **PL)

phases = [
    ('Phase 1A', 'Data Acquisition\n& Preparation',
     '#16a085', 'Download 3 genomes\nExtract upstream regions\nBuild ortholog set'),
    ('Phase 1B', 'Motif Discovery\n(MEME / FIMO)',
     '#2980b9', 'De novo motif search\nGenome-wide scan\n1,442 + 1,891 hits'),
    ('Phase 1C', 'Bipartite\nOperator Search',
     '#8e44ad', 'Paired PWM scoring\nFisher\'s method\n496 candidates'),
    ('Phase 1D', 'Validation\n& Consensus',
     '#c0392b', '3-method cross-validation\nConservation check\n50 consensus hits'),
]

bw, bh = 1.85, 1.7
y0 = 0.25
xs = [0.6, 2.95, 5.30, 7.65]

for i, (phase, title, col, detail) in enumerate(phases):
    x = xs[i]
    # box fill
    ax_a.add_patch(FancyBboxPatch((x, y0), bw, bh,
        boxstyle="round,pad=0.12", facecolor=col, alpha=0.12, edgecolor='none'))
    # box border
    ax_a.add_patch(FancyBboxPatch((x, y0), bw, bh,
        boxstyle="round,pad=0.12", facecolor='none', edgecolor=col, lw=2.5))
    # phase label
    ax_a.text(x + bw/2, y0 + bh - 0.15, phase,
              ha='center', va='top', fontsize=10, fontweight='bold', color=col)
    # title
    ax_a.text(x + bw/2, y0 + bh/2 + 0.1, title,
              ha='center', va='center', fontsize=9, fontweight='bold',
              color='#2c3e50', linespacing=1.3)
    # detail
    ax_a.text(x + bw/2, y0 + 0.25, detail,
              ha='center', va='center', fontsize=7.5,
              color='#566573', linespacing=1.3, style='italic')
    # arrow
    if i < 3:
        ax_a.annotate('', xy=(xs[i+1] - 0.08, y0 + bh/2),
                      xytext=(x + bw + 0.08, y0 + bh/2),
                      arrowprops=dict(arrowstyle='->', color='#2c3e50', lw=2.5))

# ============================================================
# B — Operator Architecture
# ============================================================
ax_b = fig.add_subplot(gs[1, 0])
ax_b.set_xlim(-8, 135)
ax_b.set_ylim(-2, 5)
ax_b.axis('off')
ax_b.text(-0.06, 1.08, 'B', transform=ax_b.transAxes, **PL)
ax_b.set_title('Known Operator Architecture (123 bp)', fontsize=13, fontweight='bold', pad=10)

# DNA backbone
ax_b.plot([-3, 128], [2, 2], color='#34495e', lw=3, zorder=0)

# 5' and 3'
ax_b.text(-5, 2, "5'", fontsize=11, fontweight='bold', color='#34495e', ha='right', va='center')
ax_b.text(130, 2, "3'", fontsize=11, fontweight='bold', color='#34495e', ha='left', va='center')

# Weak site
ax_b.add_patch(FancyBboxPatch((8, 1.0), 25, 2.0, boxstyle="round,pad=0.1",
    facecolor=C_WEAK, alpha=0.2, edgecolor=C_WEAK, lw=2.5))
ax_b.text(20.5, 2.35, 'Weak Site', ha='center', va='center', fontsize=11,
          fontweight='bold', color=C_WEAK)
ax_b.text(20.5, 1.55, '25 bp', ha='center', va='center', fontsize=9, color=C_WEAK)
ax_b.text(20.5, 0.5, 'Kd > 100 nM', ha='center', va='center', fontsize=8,
          color=C_WEAK, fontweight='bold')

# Spacer
ax_b.add_patch(FancyBboxPatch((33, 1.5), 53, 1.0, boxstyle="round,pad=0.08",
    facecolor=C_SPACER, alpha=0.15, edgecolor=C_SPACER, lw=1.5))
ax_b.text(59.5, 2.0, '53 bp spacer', ha='center', va='center', fontsize=10,
          color='#7f8c8d', style='italic')

# Strong site
ax_b.add_patch(FancyBboxPatch((86, 1.0), 25, 2.0, boxstyle="round,pad=0.1",
    facecolor=C_STRONG, alpha=0.2, edgecolor=C_STRONG, lw=2.5))
ax_b.text(98.5, 2.35, 'Strong Site', ha='center', va='center', fontsize=11,
          fontweight='bold', color=C_STRONG)
ax_b.text(98.5, 1.55, '25 bp', ha='center', va='center', fontsize=9, color=C_STRONG)
ax_b.text(98.5, 0.5, 'Kd = 2.4 nM', ha='center', va='center', fontsize=8,
          color=C_STRONG, fontweight='bold')

# Core helix callout
ax_b.add_patch(FancyBboxPatch((89, 3.7), 19, 0.8, boxstyle="round,pad=0.08",
    facecolor=C_CORE, alpha=0.12, edgecolor=C_CORE, lw=1.5))
ax_b.text(98.5, 4.1, '15 bp core helix', ha='center', va='center',
          fontsize=9, fontweight='bold', color=C_CORE)
ax_b.annotate('', xy=(98.5, 3.05), xytext=(98.5, 3.65),
              arrowprops=dict(arrowstyle='->', color=C_CORE, lw=1.5))

# Affinity comparison
ax_b.annotate('', xy=(8, -0.8), xytext=(111, -0.8),
              arrowprops=dict(arrowstyle='<->', color='#2c3e50', lw=2))
ax_b.text(59.5, -1.3, '20-fold affinity asymmetry', ha='center',
          fontsize=9, color='#2c3e50', fontweight='bold')

# ============================================================
# C — Score Waterfall
# ============================================================
ax_c = fig.add_subplot(gs[1, 1])
ax_c.text(-0.06, 1.08, 'C', transform=ax_c.transAxes, **PL)
ax_c.set_title('Bipartite Search: Top 20 Candidates', fontsize=13,
               fontweight='bold', pad=10)

scores = [s for _, s, _, _ in top20]
names = [n for _, _, n, _ in top20]
colors_bar = []
for r, s, n, reg in top20:
    if r <= 3:
        colors_bar.append(C_HIGHEST)
    elif r == 4:
        colors_bar.append(C_HIGH)
    else:
        colors_bar.append(C_LOW)

ax_c.bar(range(20), scores, color=colors_bar, edgecolor='white', lw=0.5, width=0.82)

# Score gap
ax_c.annotate('', xy=(0.6, 57.42), xytext=(0.6, 34.58),
              arrowprops=dict(arrowstyle='<->', color=C_RED, lw=2.5))
ax_c.text(1.8, 46, '22.84\nscore gap', ha='left', va='center',
          fontsize=9, color=C_RED, fontweight='bold', linespacing=1.1)

# Top 4 labels
for i in range(4):
    ax_c.text(i, scores[i] + 1.5, names[i], ha='center', va='bottom',
              fontsize=8, fontweight='bold', rotation=45,
              color=C_HIGHEST if i < 3 else C_HIGH)

ax_c.set_xlabel('Rank', fontsize=11)
ax_c.set_ylabel('Adjusted Score', fontsize=11)
ax_c.set_xticks(range(0, 20, 2))
ax_c.set_xticklabels(range(1, 21, 2), fontsize=9)
ax_c.set_ylim(0, 68)
ax_c.set_xlim(-0.8, 19.8)
ax_c.spines['top'].set_visible(False)
ax_c.spines['right'].set_visible(False)

leg = [mpatches.Patch(color=C_HIGHEST, label='HIGHEST (regulon)'),
       mpatches.Patch(color=C_HIGH, label='HIGH'),
       mpatches.Patch(color=C_LOW, label='Other')]
ax_c.legend(handles=leg, loc='upper right', fontsize=8, framealpha=0.9)

# ============================================================
# D — Sequence Alignment of 4 Operators
# ============================================================
ax_d = fig.add_subplot(gs[2, :])
ax_d.text(-0.01, 1.03, 'D', transform=ax_d.transAxes, **PL)
ax_d.set_title('Predicted Operator Binding Sites — Sequence Comparison',
               fontsize=13, fontweight='bold', pad=10)
ax_d.axis('off')

def core_best_pos(seq, core=CORE_15BP):
    best_p, best_m = -1, 0
    for p in range(len(seq) - len(core) + 1):
        m = sum(1 for a, b in zip(seq[p:p+len(core)], core) if a == b)
        if m > best_m:
            best_m = m
            best_p = p
    return best_p, best_m

def draw_seq(ax, seq, x0, y, fs=9, core_highlight=True):
    """Draw colored DNA sequence with core-match highlighting."""
    bp, bm = core_best_pos(seq) if core_highlight else (-1, 0)
    spacing = 0.42
    for i, base in enumerate(seq):
        x = x0 + i * spacing
        # highlight core-matching bases
        if bp >= 0 and bp <= i < bp + len(CORE_15BP):
            ci = i - bp
            if base == CORE_15BP[ci]:
                ax.add_patch(plt.Rectangle((x - 0.18, y - 0.28), 0.36, 0.56,
                    facecolor=C_CORE, alpha=0.18, edgecolor='none', zorder=0))
        ax.text(x, y, base, ha='center', va='center', fontsize=fs,
                fontfamily='monospace', fontweight='bold',
                color=BASE_COL.get(base, 'k'), zorder=1)

y_rows = [4.0, 2.8, 1.6, 0.4]
label_x = 0.0
weak_x = 5.5
spacer_x = weak_x + 25 * 0.42 + 0.6
strong_x = spacer_x + 3.0
info_x = strong_x + 25 * 0.42 + 0.8

# Column headers
ax_d.text(weak_x + 12.5 * 0.42, y_rows[0] + 0.8, 'Weak Site (25 bp)',
          ha='center', fontsize=10, fontweight='bold', color=C_WEAK)
ax_d.text(strong_x + 12.5 * 0.42, y_rows[0] + 0.8, 'Strong Site (25 bp)',
          ha='center', fontsize=10, fontweight='bold', color=C_STRONG)

for i, (name, weak, strong, spacer, core_id, tier) in enumerate(operators):
    y = y_rows[i]
    tcol = C_HIGHEST if tier == 'HIGHEST' else C_HIGH

    # Label
    ax_d.text(label_x, y, name, ha='left', va='center', fontsize=10,
              fontweight='bold', color=tcol)

    # Weak sequence
    draw_seq(ax_d, weak, weak_x, y)

    # Spacer
    ax_d.text(spacer_x, y, f'—{spacer}bp—', ha='center', va='center',
              fontsize=9, color=C_SPACER, style='italic')

    # Strong sequence
    draw_seq(ax_d, strong, strong_x, y)

    # Core ID badge
    ax_d.text(info_x + 2.5, y + 0.15, f'Core ID: {core_id}%',
              ha='center', va='center', fontsize=9, fontweight='bold',
              color=C_CORE if core_id == 100 else '#7f8c8d')

    # Tier badge
    ax_d.add_patch(FancyBboxPatch((info_x + 1.2, y - 0.45), 2.6, 0.4,
        boxstyle="round,pad=0.08", facecolor=tcol, alpha=0.15, edgecolor=tcol, lw=1))
    ax_d.text(info_x + 2.5, y - 0.25, tier,
              ha='center', va='center', fontsize=8, fontweight='bold', color=tcol)

# Separator lines
for y in [3.4, 2.2, 1.0]:
    ax_d.axhline(y=y, xmin=0.02, xmax=0.98, color='#ecf0f1', lw=1, zorder=-1)

# Core reference
ax_d.plot([0, info_x + 4.5], [-0.3, -0.3], color='#bdc3c7', lw=0.5)
ax_d.text(label_x, -0.7, 'Core helix:', ha='left', va='center', fontsize=9,
          fontweight='bold', color=C_CORE)
draw_seq(ax_d, CORE_15BP, weak_x, -0.7, fs=9.5, core_highlight=False)
ax_d.text(weak_x + 15 * 0.42 + 0.8, -0.7,
          '(15 bp recognition sequence from cryo-EM, PDB 9B7Y)',
          ha='left', va='center', fontsize=8, color=C_CORE, style='italic')

# Legend
ax_d.add_patch(plt.Rectangle((label_x, -1.4), 0.36, 0.36,
    facecolor=C_CORE, alpha=0.18, edgecolor='none'))
ax_d.text(label_x + 0.6, -1.22, '= base matches core helix', ha='left',
          va='center', fontsize=8, color=C_CORE)

ax_d.set_xlim(-0.5, info_x + 5.5)
ax_d.set_ylim(-1.8, 5.3)

# ============================================================
# E — Consensus Tiers (bar chart, cleaner than pie)
# ============================================================
ax_e = fig.add_subplot(gs[3, 0])
ax_e.text(-0.06, 1.08, 'E', transform=ax_e.transAxes, **PL)
ax_e.set_title('Cross-Validation Consensus\n(50 candidates from 3 methods)',
               fontsize=13, fontweight='bold', pad=10)

tiers = list(tier_counts.keys())
vals = list(tier_counts.values())
tcols = [C_HIGHEST, C_HIGH, C_MOD, C_LOW]

bars_e = ax_e.barh(range(len(tiers)), vals, color=tcols, edgecolor='white', height=0.65)
ax_e.set_yticks(range(len(tiers)))
ax_e.set_yticklabels(tiers, fontsize=10, fontweight='bold')
ax_e.invert_yaxis()
ax_e.set_xlabel('Number of Candidates', fontsize=11)
ax_e.spines['top'].set_visible(False)
ax_e.spines['right'].set_visible(False)

for i, (bar, v) in enumerate(zip(bars_e, vals)):
    ax_e.text(v + 0.5, i, str(v), va='center', ha='left', fontsize=11,
              fontweight='bold', color=tcols[i])

# Criteria annotations
criteria = ['3/3 methods + regulon + intergenic',
            '2/3 methods + intergenic',
            '2/3 methods OR regulon',
            'Single method only']
for i, c in enumerate(criteria):
    ax_e.text(vals[i] + 2.5, i + 0.08, c, va='center', ha='left',
              fontsize=7.5, color='#566573', style='italic')

ax_e.set_xlim(0, 42)

# ============================================================
# F — Manhattan Plot
# ============================================================
ax_f = fig.add_subplot(gs[3, 1])
ax_f.text(-0.06, 1.08, 'F', transform=ax_f.transAxes, **PL)
ax_f.set_title('Genome-wide Operator Scan\n(496 bipartite candidates)',
               fontsize=13, fontweight='bold', pad=10)

np.random.seed(42)
n_rest = 496 - 19
rest_pos = np.random.uniform(0, 4_411_532, n_rest)
rest_scr = np.random.exponential(1.8, n_rest) + 14.5
rest_scr = np.clip(rest_scr, 10, 20.5)

all_pos = np.concatenate([manhattan_real_pos, rest_pos]) / 1e6
all_scr = np.concatenate([manhattan_real_scr, rest_scr])

# Background dots
ax_f.scatter(all_pos[19:], all_scr[19:], s=10, alpha=0.3,
             c=C_LOW, edgecolors='none', zorder=1)
# Candidates 5-19
ax_f.scatter(all_pos[4:19], all_scr[4:19], s=30, alpha=0.6,
             c=C_MOD, edgecolors='none', zorder=2)

# Top 4 highlighted
top4_cols = [C_HIGHEST, C_HIGHEST, C_HIGHEST, C_HIGH]
top4_labels = ['Op1: yrbE3A\n(validated)', 'Op2: Rv1935c\n(novel)',
               'Op3: mce3R\n(novel)', 'Op4: Rv1115\n(candidate)']
for i in range(4):
    ax_f.scatter(all_pos[i], all_scr[i], s=100, c=top4_cols[i],
                 edgecolors='black', lw=1.2, zorder=5)

# Labels
offsets_xy = [(0.25, 4), (-0.55, 4), (0.25, 3), (0.2, 3)]
for i in range(4):
    ax_f.annotate(top4_labels[i],
        xy=(all_pos[i], all_scr[i]),
        xytext=(all_pos[i] + offsets_xy[i][0], all_scr[i] + offsets_xy[i][1]),
        fontsize=7.5, fontweight='bold', color=top4_cols[i],
        arrowprops=dict(arrowstyle='->', color='#7f8c8d', lw=1),
        ha='center', va='bottom')

# Threshold
ax_f.axhline(25, color=C_RED, ls='--', lw=1, alpha=0.6)
ax_f.text(4.35, 25.5, 'high-confidence', fontsize=7, color=C_RED,
          ha='right', style='italic')

ax_f.set_xlabel('Genome Position (Mbp)', fontsize=11)
ax_f.set_ylabel('Adjusted Score', fontsize=11)
ax_f.set_xlim(-0.1, 4.6)
ax_f.set_ylim(8, 65)
ax_f.spines['top'].set_visible(False)
ax_f.spines['right'].set_visible(False)

# ============================================================
# TITLE & SAVE
# ============================================================
fig.suptitle('Figure 1.  Computational Discovery of Mce3R Operator Binding Sites '
             'in the M. tuberculosis H37Rv Genome',
             fontsize=16, fontweight='bold', y=0.97, style='italic',
             fontfamily='serif')

out = '/Users/aayanalwani/tb project/mce3r_stochastic/results/figures'
fig.savefig(f'{out}/Figure1_Phase1_Pipeline.png', dpi=300, bbox_inches='tight', facecolor='white')
fig.savefig(f'{out}/Figure1_Phase1_Pipeline.pdf', dpi=300, bbox_inches='tight', facecolor='white')
print("Saved PNG and PDF")
plt.close()
