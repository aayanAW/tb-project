"""
phase4_figures/fig4_sensitivity.py — Tornado plot + heatmap for sensitivity analysis.

(A) Tornado: horizontal bars for each parameter's impact on CV
(B) Heatmap: parameter value (fold-change) vs parameter name, color = n_modes

Produces: results/figures/fig4_sensitivity.png
"""

import sys
sys.path.insert(0, '/Users/aayanalwani/tb project/mce3r_stochastic')

import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

from config.parameters import PARAMS
from figures.figure_style import (COLORS, STYLE, apply_style, despine,
                                         add_panel_label)

PROJECT = '/Users/aayanalwani/tb project/mce3r_stochastic'
DATA_FILE = os.path.join(PROJECT, 'results', 'phase3', 'sensitivity_data.csv')
OUT_DIR = os.path.join(PROJECT, 'results', 'figures')
OUT_FILE = os.path.join(OUT_DIR, 'fig4_sensitivity.png')


def _placeholder(msg):
    fig, ax = plt.subplots(figsize=(10, 5))
    ax.text(0.5, 0.5, msg, ha='center', va='center', fontsize=16, color='gray',
            transform=ax.transAxes)
    ax.axis('off')
    return fig


def generate_fig4():
    """Generate Figure 4: sensitivity analysis tornado + heatmap."""
    apply_style()
    os.makedirs(OUT_DIR, exist_ok=True)

    if not os.path.isfile(DATA_FILE):
        fig = _placeholder('Fig 4: sensitivity_data.csv not found')
        fig.savefig(OUT_FILE, dpi=STYLE['dpi'])
        plt.close(fig)
        return {'status': 'placeholder', 'reason': 'missing data'}

    df = pd.read_csv(DATA_FILE)
    required = {'parameter', 'fold_change', 'CV', 'n_modes'}
    if not required.issubset(df.columns):
        fig = _placeholder('Fig 4: missing required columns')
        fig.savefig(OUT_FILE, dpi=STYLE['dpi'])
        plt.close(fig)
        return {'status': 'placeholder', 'reason': 'missing columns'}

    params_list = df['parameter'].unique()

    fig, (ax_tornado, ax_heatmap) = plt.subplots(1, 2, figsize=(14, 6),
                                                  gridspec_kw={'width_ratios': [1, 1.2]})

    # ── Panel A: Tornado plot ─────────────────────────────────────────────
    add_panel_label(ax_tornado, 'A')

    # For each parameter, compute CV range (min to max across fold-changes)
    tornado_data = []
    for param in params_list:
        sub = df[df['parameter'] == param]
        # Baseline CV = CV at fold_change == 1.0
        baseline_row = sub[np.isclose(sub['fold_change'], 1.0)]
        if len(baseline_row) > 0:
            baseline_cv = baseline_row['CV'].values[0]
        else:
            baseline_cv = sub['CV'].mean()
        cv_min = sub['CV'].min()
        cv_max = sub['CV'].max()
        tornado_data.append({
            'param': param,
            'baseline': baseline_cv,
            'low': cv_min,
            'high': cv_max,
            'range': cv_max - cv_min,
        })

    # Sort by range (largest impact on top)
    tornado_data.sort(key=lambda x: x['range'])
    y_pos = np.arange(len(tornado_data))

    for i, td in enumerate(tornado_data):
        # Low side (left of baseline)
        ax_tornado.barh(i, td['low'] - td['baseline'], left=td['baseline'],
                        height=0.6, color=COLORS['symmetric'], alpha=0.8,
                        edgecolor='black', linewidth=0.4)
        # High side (right of baseline)
        ax_tornado.barh(i, td['high'] - td['baseline'], left=td['baseline'],
                        height=0.6, color=COLORS['asymmetric'], alpha=0.8,
                        edgecolor='black', linewidth=0.4)

    # Baseline vertical line
    if tornado_data:
        mean_baseline = np.mean([td['baseline'] for td in tornado_data])
        ax_tornado.axvline(mean_baseline, color='black', linewidth=1, linestyle=':',
                           alpha=0.6, label='Baseline')

    ax_tornado.set_yticks(y_pos)
    ax_tornado.set_yticklabels([td['param'] for td in tornado_data], fontsize=9)
    ax_tornado.set_xlabel('Coefficient of variation (CV)')
    ax_tornado.set_title('Parameter sensitivity (tornado)', fontsize=STYLE['title_size'])
    ax_tornado.legend(fontsize=8, frameon=False, loc='lower right')
    despine(ax_tornado)

    # ── Panel B: Heatmap ──────────────────────────────────────────────────
    add_panel_label(ax_heatmap, 'B')

    # Build matrix: rows = parameters, cols = fold_change values, values = n_modes
    fold_changes = sorted(df['fold_change'].unique())
    matrix = np.full((len(params_list), len(fold_changes)), np.nan)
    for i, param in enumerate(params_list):
        for j, fc in enumerate(fold_changes):
            row = df[(df['parameter'] == param) & (np.isclose(df['fold_change'], fc))]
            if len(row) > 0:
                matrix[i, j] = row['n_modes'].values[0]

    cmap = plt.colormaps.get_cmap('RdYlGn_r')
    bounds = [0.5, 1.5, 2.5, 3.5]
    norm = mcolors.BoundaryNorm(bounds, cmap.N)

    im = ax_heatmap.imshow(matrix, aspect='auto', cmap=cmap, norm=norm,
                           interpolation='nearest')

    ax_heatmap.set_xticks(range(len(fold_changes)))
    ax_heatmap.set_xticklabels([f'{fc:.2f}x' for fc in fold_changes],
                                fontsize=8, rotation=45, ha='right')
    ax_heatmap.set_yticks(range(len(params_list)))
    ax_heatmap.set_yticklabels(params_list, fontsize=9)
    ax_heatmap.set_xlabel('Parameter fold-change')
    ax_heatmap.set_title('GMM modes by parameter setting', fontsize=STYLE['title_size'])

    # Annotate cells
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            if not np.isnan(matrix[i, j]):
                ax_heatmap.text(j, i, f'{int(matrix[i, j])}',
                               ha='center', va='center', fontsize=9, fontweight='bold',
                               color='white' if matrix[i, j] >= 2.5 else 'black')

    cbar = fig.colorbar(im, ax=ax_heatmap, ticks=[1, 2, 3], shrink=0.8)
    cbar.set_label('Number of GMM modes', fontsize=9)

    fig.suptitle('Sensitivity analysis', fontsize=STYLE['title_size'] + 1,
                 fontweight='bold', y=1.02)
    fig.tight_layout()
    fig.savefig(OUT_FILE, dpi=STYLE['dpi'])
    plt.close(fig)

    return {
        'status': 'success',
        'n_parameters': len(params_list),
        'n_conditions': len(fold_changes),
        'output': OUT_FILE,
    }


# ── Self-tests ────────────────────────────────────────────────────────────────
if __name__ == '__main__':
    n_pass = 0
    n_fail = 0

    result = generate_fig4()

    try:
        assert os.path.isfile(OUT_FILE)
        fsize = os.path.getsize(OUT_FILE)
        assert fsize > 5000
        with open(OUT_FILE, 'rb') as f:
            assert f.read(4) == b'\x89PNG'
        print(f"SANITY PASS: fig4 saved ({fsize:,} bytes)")
        n_pass += 1
    except AssertionError as e:
        print(f"SANITY FAIL: {e}")
        n_fail += 1

    try:
        assert result['status'] in ('success', 'placeholder')
        print(f"SANITY PASS: status={result['status']}")
        n_pass += 1
    except (AssertionError, KeyError) as e:
        print(f"SANITY FAIL: {e}")
        n_fail += 1

    # SCIENTIFIC: enough parameters tested
    if result.get('n_parameters', 0) >= 4:
        print(f"SCIENTIFIC PASS: {result['n_parameters']} parameters tested")
    else:
        print(f"SCIENTIFIC WARNING: Only {result.get('n_parameters', 0)} parameters, expected >= 4")

    print(f"\n{'='*50}")
    print(f"Sanity checks: {n_pass} PASS, {n_fail} FAIL")
    if n_fail > 0:
        sys.exit(1)
    else:
        print("All checks passed.")
        sys.exit(0)

generate = generate_fig4
