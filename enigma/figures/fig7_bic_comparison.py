"""
phase4_figures/fig7_bic_comparison.py — BIC grouped bar chart.

X: condition (A, B, C, D), Y: BIC value
3 bars per condition: 1, 2, 3-component GMM
Star on lowest (best) BIC.

Produces: results/figures/fig7_bic_comparison.png
"""

import sys
sys.path.insert(0, '/Users/aayanalwani/tb project/mce3r_stochastic')

import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from sklearn.mixture import GaussianMixture

from config.parameters import PARAMS
from figures.figure_style import (COLORS, STYLE, CONDITION_COLORS,
                                         CONDITION_LABELS, apply_style,
                                         despine, add_panel_label)

PROJECT = '/Users/aayanalwani/tb project/mce3r_stochastic'
METRICS_FILE = os.path.join(PROJECT, 'results', 'phase3', 'noise_metrics.csv')
OUT_DIR = os.path.join(PROJECT, 'results', 'figures')
OUT_FILE = os.path.join(OUT_DIR, 'fig7_bic_comparison.png')


def _compute_bics(cond):
    """Compute BIC for 1, 2, 3-component GMMs for a condition."""
    path = os.path.join(PROJECT, 'results', 'phase2', f'condition_{cond}.npz')
    if not os.path.isfile(path):
        return None
    d = np.load(path, allow_pickle=True)
    proteins = d['proteins'].astype(float).reshape(-1, 1)
    if len(proteins) < 3:
        return None

    bics = []
    for k in [1, 2, 3]:
        gmm = GaussianMixture(n_components=k, reg_covar=PARAMS['gmm_reg_covar'],
                               n_init=PARAMS['gmm_n_init'],
                               random_state=PARAMS['gmm_random_state'])
        gmm.fit(proteins)
        bics.append(gmm.bic(proteins))
    return bics


def generate_fig7():
    """Generate Figure 7: BIC comparison."""
    apply_style()
    os.makedirs(OUT_DIR, exist_ok=True)

    conditions = ['A', 'B', 'C', 'D']
    all_bics = {}
    for cond in conditions:
        bics = _compute_bics(cond)
        if bics is not None:
            all_bics[cond] = bics

    if not all_bics:
        fig, ax = plt.subplots(figsize=(8, 5))
        ax.text(0.5, 0.5, 'Fig 7: No BIC data available', ha='center',
                va='center', fontsize=16, color='gray', transform=ax.transAxes)
        ax.axis('off')
        fig.savefig(OUT_FILE, dpi=STYLE['dpi'])
        plt.close(fig)
        return {'status': 'placeholder', 'reason': 'no data'}

    fig, ax = plt.subplots(figsize=(9, 5.5))

    n_conds = len(all_bics)
    n_components = 3
    width = 0.22
    x = np.arange(n_conds)

    component_colors = ['#93C5FD', '#3B82F6', '#1E3A5F']  # light to dark blue
    component_labels = ['1 component', '2 components', '3 components']

    cond_list = list(all_bics.keys())

    for k_idx in range(n_components):
        offsets = x + (k_idx - 1) * width
        bic_vals = [all_bics[c][k_idx] for c in cond_list]
        bars = ax.bar(offsets, bic_vals, width, color=component_colors[k_idx],
                      edgecolor='black', linewidth=0.5,
                      label=component_labels[k_idx])

    # Mark best (lowest) BIC with a star
    for ci, cond in enumerate(cond_list):
        bics = all_bics[cond]
        best_k = np.argmin(bics)
        best_bic = bics[best_k]
        star_x = ci + (best_k - 1) * width
        ax.annotate('*', xy=(star_x, best_bic), fontsize=22, fontweight='bold',
                    ha='center', va='bottom', color=COLORS['mce3r_marker'])

    ax.set_xticks(x)
    ax.set_xticklabels([CONDITION_LABELS.get(c, c) for c in cond_list],
                        fontsize=10)
    ax.set_ylabel('BIC')
    ax.set_title('Bayesian Information Criterion: GMM model selection',
                 fontsize=STYLE['title_size'], fontweight='bold')
    ax.legend(fontsize=9, frameon=False)
    despine(ax)

    fig.tight_layout()
    fig.savefig(OUT_FILE, dpi=STYLE['dpi'])
    plt.close(fig)

    # Determine best k per condition
    best_k_map = {c: int(np.argmin(all_bics[c])) + 1 for c in cond_list}

    return {
        'status': 'success',
        'conditions': cond_list,
        'best_k': best_k_map,
        'output': OUT_FILE,
    }


# ── Self-tests ────────────────────────────────────────────────────────────────
if __name__ == '__main__':
    n_pass = 0
    n_fail = 0

    result = generate_fig7()

    try:
        assert os.path.isfile(OUT_FILE)
        fsize = os.path.getsize(OUT_FILE)
        assert fsize > 3000
        with open(OUT_FILE, 'rb') as f:
            assert f.read(4) == b'\x89PNG'
        print(f"SANITY PASS: fig7 saved ({fsize:,} bytes)")
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

    # SCIENTIFIC: best_k for asymmetric (A) should be 2 (bimodal)
    best_k = result.get('best_k', {})
    if best_k.get('A') == 2:
        print(f"SCIENTIFIC PASS: Condition A best k = 2 (bimodal)")
    else:
        print(f"SCIENTIFIC WARNING: Condition A best k = {best_k.get('A')}, expected 2 for bimodality")

    if best_k.get('D') == 1:
        print(f"SCIENTIFIC PASS: Condition D best k = 1 (unimodal)")
    else:
        print(f"SCIENTIFIC WARNING: Condition D best k = {best_k.get('D')}, expected 1 for unregulated")

    print(f"\n{'='*50}")
    print(f"Sanity checks: {n_pass} PASS, {n_fail} FAIL")
    if n_fail > 0:
        sys.exit(1)
    else:
        print("All checks passed.")
        sys.exit(0)

generate = generate_fig7
