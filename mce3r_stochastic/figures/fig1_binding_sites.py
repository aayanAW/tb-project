"""
phase4_figures/fig1_binding_sites.py — Circular genome plot of H37Rv with predicted binding sites.

Produces:
  results/figures/fig1_binding_sites.png

Panel layout:
  - Main: circular (polar) plot of H37Rv genome with predicted Mce3R binding sites
  - Inset: bar chart of top 15 sites by FIMO score
"""

import sys
sys.path.insert(0, '/Users/aayanalwani/tb project/mce3r_stochastic')

import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch
import matplotlib.cm as cm
import matplotlib.colors as mcolors

from config.parameters import PARAMS
from figures.figure_style import COLORS, STYLE, apply_style, despine, add_panel_label

# ── Paths ─────────────────────────────────────────────────────────────────────
PROJECT = '/Users/aayanalwani/tb project/mce3r_stochastic'
DATA_FILE = os.path.join(PROJECT, 'results', 'phase1', 'predicted_sites.csv')
OUT_DIR = os.path.join(PROJECT, 'results', 'figures')
OUT_FILE = os.path.join(OUT_DIR, 'fig1_binding_sites.png')

# Known regulon genes for labelling
REGULON_GENES = set(PARAMS.get('mce3_operon_genes', []))
REGULON_GENES.update(PARAMS.get('regulon_genes_1', []))
REGULON_GENES.update(PARAMS.get('regulon_genes_2', []))
REGULON_GENES.add(PARAMS.get('mce3R_gene', 'Rv1963c'))

GENOME_SIZE = PARAMS['H37Rv_genome_size']


def _placeholder_figure(msg='DATA NOT AVAILABLE'):
    """Return a placeholder figure when data is missing."""
    fig, ax = plt.subplots(figsize=(8, 8))
    ax.text(0.5, 0.5, msg, ha='center', va='center', fontsize=18,
            color='gray', transform=ax.transAxes)
    ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    ax.axis('off')
    return fig


def generate_fig1():
    """Generate Figure 1: Circular genome plot with binding sites."""
    apply_style()
    os.makedirs(OUT_DIR, exist_ok=True)

    # ── Load data ─────────────────────────────────────────────────────────
    if not os.path.isfile(DATA_FILE):
        fig = _placeholder_figure('Fig 1: predicted_sites.csv not found')
        fig.savefig(OUT_FILE, dpi=STYLE['dpi'])
        plt.close(fig)
        return {'status': 'placeholder', 'reason': 'missing data file'}

    df = pd.read_csv(DATA_FILE)
    required_cols = {'start', 'p_value', 'score'}
    if not required_cols.issubset(df.columns):
        fig = _placeholder_figure('Fig 1: missing required columns')
        fig.savefig(OUT_FILE, dpi=STYLE['dpi'])
        plt.close(fig)
        return {'status': 'placeholder', 'reason': 'missing columns'}

    is_mock = 'MEME_mock' in str(df['motif_id'].iloc[0]) if 'motif_id' in df.columns else False

    # ── Prepare coordinates ───────────────────────────────────────────────
    # Convert genomic positions to radians (0 = origin of replication)
    theta = 2 * np.pi * df['start'].values / GENOME_SIZE
    scores = df['score'].values
    p_values = df['p_value'].values

    # Color by -log10(p_value)
    neg_log_p = -np.log10(p_values.clip(min=1e-20))
    norm = mcolors.Normalize(vmin=neg_log_p.min(), vmax=neg_log_p.max())
    cmap = cm.get_cmap('YlOrRd')

    # ── Figure layout ─────────────────────────────────────────────────────
    fig = plt.figure(figsize=(10, 10))
    # Main polar plot
    ax_polar = fig.add_axes([0.1, 0.1, 0.75, 0.75], polar=True)
    # Inset bar chart (top-right)
    ax_inset = fig.add_axes([0.62, 0.65, 0.32, 0.28])

    # ── Main circular plot ────────────────────────────────────────────────
    ax_polar.set_theta_zero_location('N')
    ax_polar.set_theta_direction(-1)  # clockwise

    # Draw genome backbone (outer ring)
    ring_theta = np.linspace(0, 2 * np.pi, 500)
    ax_polar.plot(ring_theta, np.ones_like(ring_theta) * 1.0,
                  color='#CBD5E1', linewidth=3, zorder=1)

    # Plot binding sites as scatter points on the ring
    radii = np.ones(len(theta)) * 1.0
    colors = cmap(norm(neg_log_p))
    sizes = 40 + 160 * (scores / scores.max())
    ax_polar.scatter(theta, radii, c=colors, s=sizes, edgecolors='black',
                     linewidth=0.5, zorder=3, alpha=0.85)

    # Label sites near known regulon genes
    if 'nearest_gene' in df.columns:
        for i, row in df.iterrows():
            gene = row.get('nearest_gene', '')
            gene_name = row.get('nearest_gene_name', '')
            label = gene_name if pd.notna(gene_name) and gene_name else gene
            if gene in REGULON_GENES or str(gene_name) in ['yrbE3A', 'mce3R']:
                angle = theta[i]
                ax_polar.annotate(
                    label,
                    xy=(angle, 1.0),
                    xytext=(angle, 1.22),
                    fontsize=8, fontweight='bold',
                    ha='center', va='center',
                    arrowprops=dict(arrowstyle='-', color='gray', lw=0.6),
                    color=COLORS['asymmetric'],
                )

    # Mark origin of replication
    ax_polar.annotate('oriC', xy=(0, 0.85), fontsize=8, ha='center', va='center',
                      color='gray')

    # Genome size label at center
    ax_polar.text(0, 0, f'H37Rv\n{GENOME_SIZE/1e6:.1f} Mb',
                  ha='center', va='center', fontsize=11, fontweight='bold',
                  transform=ax_polar.transData)

    # Clean up polar axes
    ax_polar.set_ylim(0, 1.35)
    ax_polar.set_yticks([])
    ax_polar.set_xticks([])
    ax_polar.spines['polar'].set_visible(False)
    ax_polar.grid(False)

    # Mock data watermark
    if is_mock:
        ax_polar.text(0, 0.35, '(mock data)', ha='center', va='center',
                      fontsize=9, color='red', alpha=0.7,
                      transform=ax_polar.transData)

    # ── Colorbar ──────────────────────────────────────────────────────────
    sm = cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar_ax = fig.add_axes([0.08, 0.08, 0.25, 0.015])
    cbar = fig.colorbar(sm, cax=cbar_ax, orientation='horizontal')
    cbar.set_label('$-\\log_{10}(p)$', fontsize=9)
    cbar.ax.tick_params(labelsize=8)

    # ── Inset: top 15 sites bar chart ─────────────────────────────────────
    top_n = min(15, len(df))
    df_top = df.nlargest(top_n, 'score').sort_values('score', ascending=True)
    bar_colors = cmap(norm(-np.log10(df_top['p_value'].values.clip(min=1e-20))))

    # Build bar labels
    bar_labels = []
    for _, row in df_top.iterrows():
        gene = row.get('nearest_gene_name', row.get('nearest_gene', ''))
        if pd.isna(gene) or gene == '':
            gene = row.get('nearest_gene', f"pos {int(row['start'])}")
        bar_labels.append(str(gene))

    y_pos = np.arange(top_n)
    ax_inset.barh(y_pos, df_top['score'].values, color=bar_colors,
                  edgecolor='black', linewidth=0.4, height=0.7)
    ax_inset.set_yticks(y_pos)
    ax_inset.set_yticklabels(bar_labels, fontsize=7)
    ax_inset.set_xlabel('FIMO score', fontsize=9)
    ax_inset.set_title(f'Top {top_n} predicted sites', fontsize=10, fontweight='bold')
    despine(ax_inset)

    # ── Title ─────────────────────────────────────────────────────────────
    fig.suptitle('Predicted Mce3R binding sites across the H37Rv genome',
                 fontsize=STYLE['title_size'], fontweight='bold', y=0.97)

    fig.savefig(OUT_FILE, dpi=STYLE['dpi'])
    plt.close(fig)

    return {
        'status': 'mock_data' if is_mock else 'real_data',
        'n_sites': len(df),
        'top_score': float(df['score'].max()),
        'output': OUT_FILE,
    }


# ── Self-tests ────────────────────────────────────────────────────────────────
if __name__ == '__main__':
    n_pass = 0
    n_fail = 0

    # Generate figure
    result = generate_fig1()

    # SANITY: Output file exists
    try:
        assert os.path.isfile(OUT_FILE), f"Output file not found: {OUT_FILE}"
        fsize = os.path.getsize(OUT_FILE)
        assert fsize > 1000, f"Output file too small: {fsize} bytes"
        print(f"SANITY PASS: fig1 saved ({fsize:,} bytes)")
        n_pass += 1
    except AssertionError as e:
        print(f"SANITY FAIL: {e}")
        n_fail += 1

    # SANITY: Result dict has expected keys
    try:
        assert 'status' in result
        assert 'output' in result or 'reason' in result
        print(f"SANITY PASS: result dict valid — status={result['status']}")
        n_pass += 1
    except AssertionError as e:
        print(f"SANITY FAIL: {e}")
        n_fail += 1

    # SANITY: PNG format check
    try:
        with open(OUT_FILE, 'rb') as f:
            header = f.read(8)
        assert header[:4] == b'\x89PNG', "Not a valid PNG file"
        print("SANITY PASS: Output is valid PNG")
        n_pass += 1
    except (AssertionError, Exception) as e:
        print(f"SANITY FAIL: {e}")
        n_fail += 1

    # SCIENTIFIC: Number of sites
    if result.get('n_sites', 0) >= 5:
        print(f"SCIENTIFIC PASS: {result.get('n_sites')} binding sites plotted")
    else:
        print(f"SCIENTIFIC WARNING: Only {result.get('n_sites', 0)} sites — expected >= 5")

    print(f"\n{'='*50}")
    print(f"Sanity checks: {n_pass} PASS, {n_fail} FAIL")
    if n_fail > 0:
        sys.exit(1)
    else:
        print("All checks passed.")
        sys.exit(0)

generate = generate_fig1
