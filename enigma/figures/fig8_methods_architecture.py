"""
phase4_figures/fig8_methods_architecture.py — Pipeline architecture diagram.

Produces: results/figures/fig8_methods_architecture.png
"""

import sys
sys.path.insert(0, '/Users/aayanalwani/tb project/mce3r_stochastic')

import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
from figures.figure_style import COLORS, STYLE, apply_style

PROJECT = '/Users/aayanalwani/tb project/mce3r_stochastic'
OUT_DIR = os.path.join(PROJECT, 'results', 'figures')
OUT_FILE = os.path.join(OUT_DIR, 'fig8_methods_architecture.png')

# Phase colors
PHASE_COLORS = {
    1: '#0D9488',  # Teal (bioinformatics)
    2: '#F97316',  # Orange (simulation)
    3: '#8B5CF6',  # Purple (analysis)
    4: '#6B7280',  # Gray (figures)
}
PHASE_COLORS_LIGHT = {
    1: '#CCFBF1',
    2: '#FFF7ED',
    3: '#EDE9FE',
    4: '#F3F4F6',
}


def _draw_box(ax, x, y, w, h, color, light_color, title, items, outputs=None, result=None):
    """Draw a phase box with title, script list, outputs, and key result."""
    box = FancyBboxPatch(
        (x, y), w, h,
        boxstyle="round,pad=0.02",
        facecolor=light_color, edgecolor=color, linewidth=2
    )
    ax.add_patch(box)

    # Title bar
    title_h = 0.06
    title_box = FancyBboxPatch(
        (x, y + h - title_h), w, title_h,
        boxstyle="round,pad=0.02",
        facecolor=color, edgecolor=color, linewidth=1.5
    )
    ax.add_patch(title_box)
    ax.text(x + w / 2, y + h - title_h / 2, title,
            ha='center', va='center', fontsize=10, fontweight='bold', color='white')

    # Script items
    y_text = y + h - title_h - 0.035
    for item in items:
        ax.text(x + 0.015, y_text, item, ha='left', va='top',
                fontsize=6.5, fontfamily='monospace', color='#374151')
        y_text -= 0.028

    # Outputs section
    if outputs:
        y_text -= 0.01
        ax.plot([x + 0.01, x + w - 0.01], [y_text + 0.01, y_text + 0.01],
                color=color, linewidth=0.5, alpha=0.5)
        y_text -= 0.005
        ax.text(x + 0.015, y_text, 'Outputs:', ha='left', va='top',
                fontsize=6, fontweight='bold', color=color)
        y_text -= 0.022
        for out in outputs:
            ax.text(x + 0.015, y_text, out, ha='left', va='top',
                    fontsize=5.5, color='#6B7280')
            y_text -= 0.02

    # Key result
    if result:
        y_text -= 0.005
        ax.plot([x + 0.01, x + w - 0.01], [y_text + 0.01, y_text + 0.01],
                color=color, linewidth=0.5, alpha=0.5)
        y_text -= 0.005
        ax.text(x + w / 2, y_text, result, ha='center', va='top',
                fontsize=6.5, fontweight='bold', color=color, style='italic')


def _draw_arrow(ax, x1, y1, x2, y2, color='#374151'):
    """Draw a curved arrow between two points."""
    arrow = FancyArrowPatch(
        (x1, y1), (x2, y2),
        arrowstyle='->', mutation_scale=15,
        connectionstyle='arc3,rad=0',
        color=color, linewidth=2
    )
    ax.add_patch(arrow)


def generate_fig8():
    """Generate pipeline architecture diagram."""
    apply_style()
    os.makedirs(OUT_DIR, exist_ok=True)

    fig, ax = plt.subplots(figsize=(16, 9))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis('off')
    ax.set_aspect('equal')

    # Title
    ax.text(0.5, 0.97, 'Mce3R Stochastic Simulation Pipeline',
            ha='center', va='top', fontsize=16, fontweight='bold', color='#1F2937')
    ax.text(0.5, 0.935, 'Phases 1 & 2 run in parallel; Phase 3 depends on both; Phase 4 depends on Phase 3',
            ha='center', va='top', fontsize=9, color='#6B7280')

    # Layout constants
    box_w = 0.20
    box_h_top = 0.42    # phases 1 & 2
    box_h_bot = 0.32    # phases 3 & 4
    y_top = 0.45
    y_bot = 0.04
    gap = 0.04

    # Phase 1 — Bioinformatics (top left)
    x1 = 0.05
    _draw_box(ax, x1, y_top, box_w, box_h_top,
              PHASE_COLORS[1], PHASE_COLORS_LIGHT[1],
              'Phase 1: Bioinformatics',
              items=[
                  'download_genomes.py',
                  '  \u2514 Fetch H37Rv, M. bovis, M. marinum',
                  'extract_upstream.py',
                  '  \u2514 200bp upstream of yrbE3A orthologs',
                  'run_meme.py',
                  '  \u2514 De novo motif discovery (ZOOPS)',
                  'run_fimo.py',
                  '  \u2514 Genome-wide motif scan',
                  'conservation_check.py',
                  '  \u2514 Cross-species conservation',
              ],
              outputs=['predicted_sites.csv', 'conservation_status.csv'],
              result='Genome-wide Mce3R binding sites')

    # Phase 2 — Simulation (top right)
    x2 = 0.28
    _draw_box(ax, x2, y_top, box_w, box_h_top,
              PHASE_COLORS[2], PHASE_COLORS_LIGHT[2],
              'Phase 2: Stochastic Simulation',
              items=[
                  'operator_model.py',
                  '  \u2514 4-state operator (both/strong/weak/none)',
                  'gillespie_engine.py',
                  '  \u2514 Numba-accelerated SSA',
                  'run_conditions.py',
                  '  \u2514 A: Asymmetric (wild-type)',
                  '  \u2514 B: Symmetric (geom. mean Kd)',
                  '  \u2514 C: Single-site (weak disabled)',
                  '  \u2514 D: Unregulated (constitutive)',
                  'asymmetry_sweep.py',
                  '  \u2514 E: Kd ratio sweep (1\u201350\u00d7)',
              ],
              outputs=['condition_[A-D].npz', 'condition_E_sweep.npz'],
              result='50k cells/condition, Gillespie SSA')

    # "PARALLEL" label between phases 1 & 2
    mid_x = (x1 + box_w + x2) / 2
    ax.text(mid_x, y_top + box_h_top / 2, 'PARALLEL',
            ha='center', va='center', fontsize=8, fontweight='bold',
            color='#9CA3AF', rotation=90)

    # Phase 3 — Analysis (bottom left)
    x3 = 0.51
    _draw_box(ax, x3, y_top, box_w, box_h_top,
              PHASE_COLORS[3], PHASE_COLORS_LIGHT[3],
              'Phase 3: Statistical Analysis',
              items=[
                  'noise_metrics.py',
                  '  \u2514 CV, Fano, bimodality per condition',
                  'bootstrap_ci.py',
                  '  \u2514 95% CIs (1000 resamples)',
                  'sensitivity_analysis.py',
                  '  \u2514 7 params \u00d7 10 steps \u00d7 1000 cells',
                  'statistical_tests.py',
                  '  \u2514 KS tests, GMM likelihood ratio',
                  'experimental_comparison.py',
                  '  \u2514 Fold-change vs Santangelo 2009',
                  'negative_controls.py',
                  '  \u2514 Shuffle, symmetric TetR, Poisson',
              ],
              outputs=['noise_metrics.csv', 'bootstrap_results.csv',
                       'sensitivity_data.csv', 'statistical_tests.csv'],
              result='CV(A) > CV(B), p < 0.001')

    # Phase 4 — Figures (bottom right)
    x4 = 0.76
    _draw_box(ax, x4, y_top, box_w, box_h_top,
              PHASE_COLORS[4], PHASE_COLORS_LIGHT[4],
              'Phase 4: Publication Figures',
              items=[
                  'fig1: Genome-wide binding sites',
                  'fig2: Protein distributions (A-D)',
                  'fig3: Asymmetry sweep (CV vs ratio)',
                  'fig4: Sensitivity analysis heatmap',
                  'fig5: Model validation',
                  'fig6: Single-cell traces',
                  'fig7: BIC model comparison',
                  'fig8: Pipeline architecture (this fig)',
                  '',
                  'figure_style.py',
                  '  \u2514 Shared colors, fonts, 300 dpi',
              ],
              outputs=['fig1-8_*.png (300 dpi)', 'figures_summary.json'],
              result='8 publication-quality figures')

    # Arrows: Phase 1 -> Phase 3
    _draw_arrow(ax, x1 + box_w, y_top + box_h_top / 2,
                x3, y_top + box_h_top * 0.7,
                color=PHASE_COLORS[1])

    # Arrows: Phase 2 -> Phase 3
    _draw_arrow(ax, x2 + box_w, y_top + box_h_top / 2,
                x3, y_top + box_h_top * 0.4,
                color=PHASE_COLORS[2])

    # Arrows: Phase 3 -> Phase 4
    _draw_arrow(ax, x3 + box_w, y_top + box_h_top / 2,
                x4, y_top + box_h_top / 2,
                color=PHASE_COLORS[3])

    # Bottom: Key parameters box
    param_y = 0.02
    param_h = 0.38
    param_w = 0.96
    param_x = 0.02
    param_box = FancyBboxPatch(
        (param_x, param_y), param_w, param_h,
        boxstyle="round,pad=0.02",
        facecolor='#F9FAFB', edgecolor='#D1D5DB', linewidth=1.5
    )
    ax.add_patch(param_box)
    ax.text(param_x + param_w / 2, param_y + param_h - 0.02,
            'Key Parameters & Configuration (config/parameters.py)',
            ha='center', va='top', fontsize=10, fontweight='bold', color='#374151')

    # Parameter columns
    col_params = [
        [
            'Binding Kinetics:',
            '  Kd_strong = 2.4 nM',
            '  Kd_weak = 49.0 nM',
            '  k_on = 0.0167 nM\u207b\u00b9min\u207b\u00b9',
            '  Kd ratio = 20.4\u00d7',
        ],
        [
            'Transcription Model:',
            '  k_max = 0.15 mRNA/min',
            '  block_strong = 0.85',
            '  block_weak = 0.50',
            '  4-state operator model',
        ],
        [
            'Translation & Decay:',
            '  k_translation = 0.5 prot/mRNA/min',
            '  t\u00bd_mRNA = 9.5 min',
            '  t\u00bd_protein = 1500 min',
            '  burst_size \u2248 6.85',
        ],
        [
            'Simulation Settings:',
            '  n_cells = 50,000/condition',
            '  n_sweep = 5,000/ratio',
            '  t_max = 30,000 min',
            '  t_burn_in = 15,000 min',
        ],
        [
            'Symmetric Control (B):',
            '  Kd = \u221a(2.4\u00d749.0) \u2248 10.84 nM',
            '  block = \u221a(0.85\u00d70.50) \u2248 0.652',
            '  (geometric mean preserves',
            '   total regulatory capacity)',
        ],
    ]

    col_x_start = param_x + 0.015
    col_width = (param_w - 0.03) / len(col_params)
    for i, col in enumerate(col_params):
        cx = col_x_start + i * col_width
        cy = param_y + param_h - 0.07
        for j, line in enumerate(col):
            weight = 'bold' if j == 0 else 'normal'
            color = '#374151' if j == 0 else '#6B7280'
            ax.text(cx, cy - j * 0.025, line, ha='left', va='top',
                    fontsize=6.5, fontweight=weight, color=color)

    # Source reference
    ax.text(0.5, 0.005, 'Source: Panagoda et al. 2024 (PDB 9B7Y) | Stormo & Zhao 2010 | Rustad et al. 2013 | Taniguchi et al. 2010',
            ha='center', va='bottom', fontsize=6, color='#9CA3AF', style='italic')

    fig.savefig(OUT_FILE, dpi=STYLE['dpi'], facecolor='white')
    plt.close(fig)

    file_size_kb = os.path.getsize(OUT_FILE) / 1024
    return {
        'status': 'success',
        'output': OUT_FILE,
        'file_size_kb': round(file_size_kb, 1),
    }


if __name__ == '__main__':
    apply_style()
    result = generate_fig8()
    print(f"Generated: {result['output']} ({result['file_size_kb']:.1f} KB)")

generate = generate_fig8
