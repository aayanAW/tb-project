#!/usr/bin/env python3
"""Create schematic figures for the Mce3R research poster."""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, ArrowStyle
import numpy as np
import os

# Color palette matching poster
NAVY = '#0B1D3A'
TEAL = '#0D9488'
SEAFOAM = '#14B8A6'
CHARCOAL = '#1E293B'
LIGHT_GRAY = '#E2E8F0'
WHITE = '#FFFFFF'
CORAL = '#EF4444'
AMBER = '#F59E0B'
BLUE = '#3B82F6'
PURPLE = '#8B5CF6'
GREEN = '#10B981'
PINK = '#EC4899'

OUT_DIR = 'results/figures/poster'
os.makedirs(OUT_DIR, exist_ok=True)

DPI = 200


def fig1_repression_mechanism():
    """Figure 1: Mce3R repression mechanism schematic."""
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.5), facecolor='white')

    for ax in axes:
        ax.set_xlim(0, 10)
        ax.set_ylim(0, 8)
        ax.set_aspect('equal')
        ax.axis('off')

    # LEFT: Repressed state
    ax = axes[0]
    ax.set_title('REPRESSED STATE', fontsize=14, fontweight='bold', color=CORAL, pad=10)

    # DNA line
    ax.plot([1, 9], [2, 2], color=NAVY, linewidth=4, zorder=1)
    ax.text(5, 1.2, 'mce3 operon (cholesterol import)', ha='center', fontsize=9,
            color=CHARCOAL, style='italic')

    # Operator box on DNA
    op = FancyBboxPatch((3.5, 1.6), 3, 0.8, boxstyle="round,pad=0.1",
                        facecolor=TEAL, edgecolor=NAVY, linewidth=1.5, zorder=2)
    ax.add_patch(op)
    ax.text(5, 2, 'Operator DNA', ha='center', va='center', fontsize=8,
            color=WHITE, fontweight='bold')

    # Mce3R protein (circle) sitting ON operator
    protein = plt.Circle((5, 3.8), 1.0, facecolor=CORAL, edgecolor=NAVY,
                         linewidth=2, zorder=3)
    ax.add_patch(protein)
    ax.text(5, 3.8, 'Mce3R\nprotein', ha='center', va='center', fontsize=9,
            color=WHITE, fontweight='bold')

    # Arrow from protein to operator
    ax.annotate('', xy=(5, 2.5), xytext=(5, 2.8),
                arrowprops=dict(arrowstyle='->', color=NAVY, lw=2))

    # X mark over genes
    ax.text(5, 0.5, 'GENES OFF', ha='center', fontsize=12, fontweight='bold',
            color=CORAL)
    ax.plot([4.2, 5.8], [0.7, 0.3], color=CORAL, linewidth=3)
    ax.plot([4.2, 5.8], [0.3, 0.7], color=CORAL, linewidth=3)

    # Block arrow
    ax.annotate('', xy=(8.5, 2), xytext=(7, 2),
                arrowprops=dict(arrowstyle='-|>', color=LIGHT_GRAY, lw=3,
                               mutation_scale=20))

    # RIGHT: De-repressed state
    ax = axes[1]
    ax.set_title('DE-REPRESSED STATE', fontsize=14, fontweight='bold', color=GREEN, pad=10)

    # DNA line
    ax.plot([1, 9], [2, 2], color=NAVY, linewidth=4, zorder=1)
    ax.text(5, 1.2, 'mce3 operon (cholesterol import)', ha='center', fontsize=9,
            color=CHARCOAL, style='italic')

    # Operator box on DNA (empty)
    op = FancyBboxPatch((3.5, 1.6), 3, 0.8, boxstyle="round,pad=0.1",
                        facecolor=SEAFOAM, edgecolor=NAVY, linewidth=1.5,
                        alpha=0.5, zorder=2)
    ax.add_patch(op)
    ax.text(5, 2, 'Operator DNA', ha='center', va='center', fontsize=8,
            color=NAVY, fontweight='bold')

    # Mce3R protein floating away
    protein = plt.Circle((2, 5.5), 0.7, facecolor=CORAL, edgecolor=NAVY,
                         linewidth=1.5, alpha=0.4, zorder=3)
    ax.add_patch(protein)
    ax.text(2, 5.5, 'Mce3R', ha='center', va='center', fontsize=8,
            color=WHITE, fontweight='bold', alpha=0.7)

    # Curved arrow showing release
    ax.annotate('', xy=(2, 4.5), xytext=(4, 3),
                arrowprops=dict(arrowstyle='->', color=CHARCOAL, lw=1.5,
                               connectionstyle='arc3,rad=0.3'))
    ax.text(2.5, 4.0, 'released', fontsize=8, color=CHARCOAL, style='italic', rotation=30)

    # Active transcription arrow
    ax.annotate('', xy=(8.5, 2), xytext=(7, 2),
                arrowprops=dict(arrowstyle='-|>', color=GREEN, lw=3, mutation_scale=20))

    # mRNA squiggle
    x_mrna = np.linspace(7.5, 9, 30)
    y_mrna = 3 + 0.3 * np.sin(10 * x_mrna)
    ax.plot(x_mrna, y_mrna, color=BLUE, linewidth=2)
    ax.text(8.25, 3.6, 'mRNA', fontsize=8, color=BLUE, fontweight='bold')

    ax.text(5, 0.5, 'GENES ON', ha='center', fontsize=12, fontweight='bold',
            color=GREEN)

    plt.tight_layout()
    fig.savefig(f'{OUT_DIR}/fig1_repression_mechanism.png', dpi=DPI, bbox_inches='tight',
                facecolor='white')
    plt.close()
    print('  -> fig1_repression_mechanism.png')


def fig2_asymmetric_vs_symmetric():
    """Figure 2: Asymmetric vs symmetric operator architecture."""
    fig, ax = plt.subplots(figsize=(10, 5), facecolor='white')
    ax.set_xlim(0, 12)
    ax.set_ylim(0, 8)
    ax.axis('off')

    # --- ASYMMETRIC (top) ---
    y_top = 6.0
    ax.text(0.3, y_top + 1.0, 'ASYMMETRIC OPERATOR (Mce3R — novel)', fontsize=13,
            fontweight='bold', color=TEAL)

    # DNA backbone
    ax.plot([1, 11], [y_top, y_top], color=CHARCOAL, linewidth=3, zorder=1)

    # Strong site (big arrow)
    strong = FancyBboxPatch((2, y_top - 0.4), 2.5, 0.8, boxstyle="round,pad=0.05",
                            facecolor=CORAL, edgecolor=NAVY, linewidth=1.5, zorder=2)
    ax.add_patch(strong)
    ax.text(3.25, y_top, 'Site A\n(STRONG)', ha='center', va='center', fontsize=8,
            color=WHITE, fontweight='bold')

    # Spacer label
    ax.annotate('', xy=(7.0, y_top - 0.7), xytext=(4.6, y_top - 0.7),
                arrowprops=dict(arrowstyle='<->', color=CHARCOAL, lw=1.5))
    ax.text(5.8, y_top - 1.1, '53 bp spacer', ha='center', fontsize=9,
            color=CHARCOAL, fontweight='bold')

    # Weak site (smaller arrow, different color)
    weak = FancyBboxPatch((7.0, y_top - 0.3), 2.0, 0.6, boxstyle="round,pad=0.05",
                          facecolor=AMBER, edgecolor=NAVY, linewidth=1.5, zorder=2)
    ax.add_patch(weak)
    ax.text(8.0, y_top, 'Site B\n(WEAK)', ha='center', va='center', fontsize=8,
            color=WHITE, fontweight='bold')

    # Kd labels
    ax.text(3.25, y_top + 0.7, 'Kd = 2.4 nM', ha='center', fontsize=9,
            color=CORAL, fontweight='bold')
    ax.text(8.0, y_top + 0.6, 'Kd ~ 10-50 nM', ha='center', fontsize=9,
            color=AMBER, fontweight='bold')

    # Inequality symbol
    ax.text(5.5, y_top + 0.6, '≠', ha='center', fontsize=20, color=NAVY, fontweight='bold')

    # --- SYMMETRIC (bottom) ---
    y_bot = 2.5
    ax.text(0.3, y_bot + 1.3, 'SYMMETRIC OPERATOR (classical TetR)', fontsize=13,
            fontweight='bold', color=PURPLE)

    # DNA backbone
    ax.plot([1, 11], [y_bot, y_bot], color=CHARCOAL, linewidth=3, zorder=1)

    # Identical sites
    s1 = FancyBboxPatch((2, y_bot - 0.4), 2.5, 0.8, boxstyle="round,pad=0.05",
                        facecolor=PURPLE, edgecolor=NAVY, linewidth=1.5, zorder=2)
    ax.add_patch(s1)
    ax.text(3.25, y_bot, 'Site A', ha='center', va='center', fontsize=9,
            color=WHITE, fontweight='bold')

    s2 = FancyBboxPatch((7.0, y_bot - 0.4), 2.5, 0.8, boxstyle="round,pad=0.05",
                        facecolor=PURPLE, edgecolor=NAVY, linewidth=1.5, zorder=2)
    ax.add_patch(s2)
    ax.text(8.25, y_bot, 'Site B', ha='center', va='center', fontsize=9,
            color=WHITE, fontweight='bold')

    # Equal Kd
    ax.text(3.25, y_bot + 0.7, 'Kd = 2.4 nM', ha='center', fontsize=9,
            color=PURPLE, fontweight='bold')
    ax.text(8.25, y_bot + 0.7, 'Kd = 2.4 nM', ha='center', fontsize=9,
            color=PURPLE, fontweight='bold')
    ax.text(5.8, y_bot + 0.7, '=', ha='center', fontsize=20, color=NAVY, fontweight='bold')

    # "Palindromic" label
    ax.text(5.8, y_bot - 0.8, 'Palindromic (mirror-image) sequence', ha='center',
            fontsize=9, color=PURPLE, style='italic')

    # Divider
    ax.plot([0.5, 11.5], [4.2, 4.2], color=LIGHT_GRAY, linewidth=1.5, linestyle='--')

    fig.savefig(f'{OUT_DIR}/fig2_asymmetric_vs_symmetric.png', dpi=DPI, bbox_inches='tight',
                facecolor='white')
    plt.close()
    print('  -> fig2_asymmetric_vs_symmetric.png')


def fig3_four_operator_states():
    """Figure 3: Four operator microstates."""
    fig, ax = plt.subplots(figsize=(10, 6), facecolor='white')
    ax.set_xlim(0, 12)
    ax.set_ylim(0, 10)
    ax.axis('off')

    states = [
        ('State 1: Both empty', 'FULL ON', GREEN, 8.0, False, False),
        ('State 2: Strong occupied', 'PARTIAL ON', AMBER, 6.0, True, False),
        ('State 3: Weak occupied', 'PARTIAL ON', AMBER, 4.0, False, True),
        ('State 4: Both occupied', 'FULL OFF', CORAL, 2.0, True, True),
    ]

    for label, effect, color, y, strong_bound, weak_bound in states:
        # DNA line
        ax.plot([1.5, 7], [y, y], color=CHARCOAL, linewidth=3, zorder=1)

        # Strong site
        sc = CORAL if not strong_bound else NAVY
        sf = CORAL if not strong_bound else CORAL
        alpha_s = 0.3 if not strong_bound else 1.0
        s_box = FancyBboxPatch((2.5, y - 0.3), 1.2, 0.6, boxstyle="round,pad=0.03",
                               facecolor=CORAL, edgecolor=NAVY, linewidth=1, alpha=alpha_s, zorder=2)
        ax.add_patch(s_box)
        ax.text(3.1, y, 'S', ha='center', va='center', fontsize=8, color=WHITE,
                fontweight='bold', alpha=alpha_s)

        # Protein on strong site
        if strong_bound:
            p = plt.Circle((3.1, y + 0.8), 0.45, facecolor=TEAL, edgecolor=NAVY, linewidth=1.5, zorder=3)
            ax.add_patch(p)
            ax.text(3.1, y + 0.8, 'M', ha='center', va='center', fontsize=8, color=WHITE, fontweight='bold')

        # Weak site
        alpha_w = 0.3 if not weak_bound else 1.0
        w_box = FancyBboxPatch((5.0, y - 0.25), 1.0, 0.5, boxstyle="round,pad=0.03",
                               facecolor=AMBER, edgecolor=NAVY, linewidth=1, alpha=alpha_w, zorder=2)
        ax.add_patch(w_box)
        ax.text(5.5, y, 'W', ha='center', va='center', fontsize=8, color=WHITE,
                fontweight='bold', alpha=alpha_w)

        # Protein on weak site
        if weak_bound:
            p = plt.Circle((5.5, y + 0.7), 0.4, facecolor=TEAL, edgecolor=NAVY, linewidth=1.5, zorder=3)
            ax.add_patch(p)
            ax.text(5.5, y + 0.7, 'M', ha='center', va='center', fontsize=7, color=WHITE, fontweight='bold')

        # State label
        ax.text(0.3, y, label, ha='left', va='center', fontsize=10, color=CHARCOAL, fontweight='bold')

        # Effect box
        effect_box = FancyBboxPatch((8, y - 0.35), 3.5, 0.7, boxstyle="round,pad=0.1",
                                    facecolor=color, edgecolor=NAVY, linewidth=1.5, alpha=0.85, zorder=2)
        ax.add_patch(effect_box)
        ax.text(9.75, y, effect, ha='center', va='center', fontsize=11, color=WHITE, fontweight='bold')

    # Highlight bracket for intermediate states
    ax.annotate('', xy=(11.8, 6.5), xytext=(11.8, 3.5),
                arrowprops=dict(arrowstyle='-', color=TEAL, lw=2))
    ax.plot([11.6, 11.8], [6.5, 6.5], color=TEAL, lw=2)
    ax.plot([11.6, 11.8], [3.5, 3.5], color=TEAL, lw=2)
    ax.text(11.8, 5.0, ' NEW\n states', fontsize=9, color=TEAL, fontweight='bold',
            va='center')

    # Legend
    ax.text(3.1, 1.0, 'S = Strong site   W = Weak site   M = Mce3R protein',
            ha='center', fontsize=9, color=CHARCOAL, style='italic',
            transform=ax.transData)

    fig.savefig(f'{OUT_DIR}/fig3_four_operator_states.png', dpi=DPI, bbox_inches='tight',
                facecolor='white')
    plt.close()
    print('  -> fig3_four_operator_states.png')


def fig8_gillespie_flowchart():
    """Figure 8: Gillespie algorithm flow diagram."""
    fig, ax = plt.subplots(figsize=(10, 4), facecolor='white')
    ax.set_xlim(0, 14)
    ax.set_ylim(0, 4)
    ax.axis('off')

    boxes = [
        (1.5, 2, 'Calculate\nall reaction\nrates', NAVY),
        (4.5, 2, 'Pick WHEN\nnext event\nhappens', TEAL),
        (7.5, 2, 'Pick WHICH\nevent\nhappens', SEAFOAM),
        (10.5, 2, 'Update\nmolecule\ncounts', CORAL),
    ]

    for x, y, text, color in boxes:
        box = FancyBboxPatch((x - 1.2, y - 0.9), 2.4, 1.8, boxstyle="round,pad=0.15",
                             facecolor=color, edgecolor=NAVY, linewidth=2, zorder=2)
        ax.add_patch(box)
        ax.text(x, y, text, ha='center', va='center', fontsize=10, color=WHITE,
                fontweight='bold', zorder=3)

    # Arrows between boxes
    for i in range(3):
        x1 = boxes[i][0] + 1.3
        x2 = boxes[i+1][0] - 1.3
        ax.annotate('', xy=(x2, 2), xytext=(x1, 2),
                    arrowprops=dict(arrowstyle='-|>', color=CHARCOAL, lw=2.5,
                                   mutation_scale=20))

    # Loop-back arrow
    ax.annotate('', xy=(0.3, 2), xytext=(11.7, 2),
                arrowprops=dict(arrowstyle='-|>', color=CHARCOAL, lw=2,
                               connectionstyle='arc3,rad=-0.5', mutation_scale=15))
    ax.text(6, 0.15, 'Repeat millions of times per cell', ha='center', fontsize=10,
            color=CHARCOAL, fontweight='bold', style='italic')

    # Sub-labels
    ax.text(4.5, 0.8, 'τ = (1/R) × ln(1/rand)', ha='center', fontsize=8,
            color=CHARCOAL, family='monospace')
    ax.text(7.5, 0.8, 'weighted by rate', ha='center', fontsize=8,
            color=CHARCOAL, style='italic')

    fig.savefig(f'{OUT_DIR}/fig8_gillespie_flowchart.png', dpi=DPI, bbox_inches='tight',
                facecolor='white')
    plt.close()
    print('  -> fig8_gillespie_flowchart.png')


def fig9_four_simulations():
    """Figure 9: Four operator architectures tested."""
    fig, axes = plt.subplots(2, 2, figsize=(10, 6), facecolor='white')
    fig.subplots_adjust(hspace=0.4, wspace=0.3)

    configs = [
        (axes[0, 0], 'A: Asymmetric (Real Mce3R)', TEAL,
         [('S', 2.5, CORAL, 1.0), ('W', 5.5, AMBER, 0.6)]),
        (axes[0, 1], 'B: Symmetric (Control)', PURPLE,
         [('S', 2.5, PURPLE, 1.0), ('S', 5.5, PURPLE, 1.0)]),
        (axes[1, 0], 'C: Single Site', BLUE,
         [('S', 4.0, CORAL, 1.0)]),
        (axes[1, 1], 'D: No Autoregulation', CHARCOAL,
         []),
    ]

    for ax, title, title_color, sites in configs:
        ax.set_xlim(0, 8)
        ax.set_ylim(0, 4)
        ax.axis('off')

        ax.set_title(title, fontsize=11, fontweight='bold', color=title_color, pad=8)

        # DNA line
        ax.plot([0.5, 7.5], [2, 2], color=CHARCOAL, linewidth=3, zorder=1)

        if not sites:
            # No regulation — just an arrow
            ax.annotate('', xy=(7, 2), xytext=(5, 2),
                        arrowprops=dict(arrowstyle='-|>', color=GREEN, lw=3, mutation_scale=20))
            ax.text(4, 2.8, 'Constitutive\nexpression', ha='center', fontsize=9,
                    color=GREEN, fontweight='bold')
            ax.text(4, 1.0, 'No feedback loop', ha='center', fontsize=9,
                    color=CHARCOAL, style='italic')
        else:
            for label, x, color, size in sites:
                w = 1.2 * size
                h = 0.6 * size
                box = FancyBboxPatch((x - w/2, 2 - h/2), w, h,
                                     boxstyle="round,pad=0.03",
                                     facecolor=color, edgecolor=NAVY, linewidth=1.5, zorder=2)
                ax.add_patch(box)
                ax.text(x, 2, label, ha='center', va='center', fontsize=9,
                        color=WHITE, fontweight='bold')

            if len(sites) == 2:
                kd1 = '2.4 nM'
                kd2 = '10-50 nM' if sites[1][3] < 1 else '2.4 nM'
                ax.text(sites[0][1], 2.8, f'Kd={kd1}', ha='center', fontsize=7, color=sites[0][2])
                ax.text(sites[1][1], 2.8, f'Kd={kd2}', ha='center', fontsize=7, color=sites[1][2])
            elif len(sites) == 1:
                ax.text(sites[0][1], 2.8, 'Kd=2.4 nM', ha='center', fontsize=8, color=sites[0][2])
                ax.text(sites[0][1], 1.0, 'No paired site', ha='center', fontsize=9,
                        color=CHARCOAL, style='italic')

    fig.savefig(f'{OUT_DIR}/fig9_four_simulations.png', dpi=DPI, bbox_inches='tight',
                facecolor='white')
    plt.close()
    print('  -> fig9_four_simulations.png')


def fig10_trimodal_distributions():
    """Figure 10: Predicted protein distributions — trimodal vs bimodal."""
    fig, axes = plt.subplots(2, 2, figsize=(10, 7), facecolor='white')
    fig.subplots_adjust(hspace=0.4, wspace=0.3)

    np.random.seed(42)
    x = np.linspace(0, 500, 1000)

    def gauss(x, mu, sigma, amp):
        return amp * np.exp(-0.5 * ((x - mu) / sigma) ** 2)

    # A: Asymmetric — TRIMODAL
    ax = axes[0, 0]
    y = gauss(x, 30, 20, 0.8) + gauss(x, 180, 40, 0.5) + gauss(x, 380, 50, 0.6)
    ax.fill_between(x, y, alpha=0.3, color=TEAL)
    ax.plot(x, y, color=TEAL, linewidth=2.5)
    ax.axvline(30, color=CORAL, linestyle='--', alpha=0.5, linewidth=1)
    ax.axvline(180, color=AMBER, linestyle='--', alpha=0.5, linewidth=1)
    ax.axvline(380, color=GREEN, linestyle='--', alpha=0.5, linewidth=1)
    ax.text(30, max(y)*1.05, 'OFF', ha='center', fontsize=8, color=CORAL, fontweight='bold')
    ax.text(180, max(y)*1.05, 'INTERMEDIATE\n(persisters)', ha='center', fontsize=7,
            color=AMBER, fontweight='bold')
    ax.text(380, max(y)*1.05, 'ON', ha='center', fontsize=8, color=GREEN, fontweight='bold')
    ax.set_title('A: Asymmetric → TRIMODAL', fontsize=11, fontweight='bold', color=TEAL)
    ax.set_xlabel('Mce3R protein copies', fontsize=9)
    ax.set_ylabel('Cell frequency', fontsize=9)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # B: Symmetric — BIMODAL
    ax = axes[0, 1]
    y = gauss(x, 30, 25, 0.9) + gauss(x, 380, 55, 0.7)
    ax.fill_between(x, y, alpha=0.3, color=PURPLE)
    ax.plot(x, y, color=PURPLE, linewidth=2.5)
    ax.axvline(30, color=CORAL, linestyle='--', alpha=0.5, linewidth=1)
    ax.axvline(380, color=GREEN, linestyle='--', alpha=0.5, linewidth=1)
    ax.text(30, max(y)*1.05, 'OFF', ha='center', fontsize=8, color=CORAL, fontweight='bold')
    ax.text(380, max(y)*1.05, 'ON', ha='center', fontsize=8, color=GREEN, fontweight='bold')
    ax.text(200, max(y)*0.5, 'No intermediate\nstate', ha='center', fontsize=9,
            color=CHARCOAL, style='italic')
    ax.set_title('B: Symmetric → BIMODAL', fontsize=11, fontweight='bold', color=PURPLE)
    ax.set_xlabel('Mce3R protein copies', fontsize=9)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # C: Single site — broad unimodal
    ax = axes[1, 0]
    y = gauss(x, 200, 100, 0.6)
    ax.fill_between(x, y, alpha=0.3, color=BLUE)
    ax.plot(x, y, color=BLUE, linewidth=2.5)
    ax.set_title('C: Single Site → Broad Unimodal', fontsize=11, fontweight='bold', color=BLUE)
    ax.set_xlabel('Mce3R protein copies', fontsize=9)
    ax.set_ylabel('Cell frequency', fontsize=9)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # D: No autoregulation — narrow unimodal
    ax = axes[1, 1]
    y = gauss(x, 200, 40, 0.9)
    ax.fill_between(x, y, alpha=0.3, color=CHARCOAL)
    ax.plot(x, y, color=CHARCOAL, linewidth=2.5)
    ax.set_title('D: No Regulation → Narrow Unimodal', fontsize=11, fontweight='bold', color=CHARCOAL)
    ax.set_xlabel('Mce3R protein copies', fontsize=9)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    fig.savefig(f'{OUT_DIR}/fig10_trimodal_distributions.png', dpi=DPI, bbox_inches='tight',
                facecolor='white')
    plt.close()
    print('  -> fig10_trimodal_distributions.png')


if __name__ == '__main__':
    print('Creating poster figures...')
    fig1_repression_mechanism()
    fig2_asymmetric_vs_symmetric()
    fig3_four_operator_states()
    fig8_gillespie_flowchart()
    fig9_four_simulations()
    fig10_trimodal_distributions()
    print('Done! All figures saved to results/figures/poster/')
