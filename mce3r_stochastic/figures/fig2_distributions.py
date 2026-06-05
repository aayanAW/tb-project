"""
phase4_figures/fig2_distributions.py — Protein distributions for all conditions.

2x3 panel layout:
  (A) Condition A histogram + GMM fit
  (B) Condition B histogram + GMM fit
  (C) Condition C histogram + GMM fit
  (D) Condition D histogram + GMM fit
  (E) Box plot of CV / Fano / bimodality coefficient
  (F) Overlay of all conditions (density)

Produces: results/figures/fig2_distributions.png
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
OUT_DIR = os.path.join(PROJECT, 'results', 'figures')
OUT_FILE = os.path.join(OUT_DIR, 'fig2_distributions.png')


def _load_condition(cond):
    """Load protein counts for a condition. Returns array or None."""
    path = os.path.join(PROJECT, 'results', 'phase2', f'condition_{cond}.npz')
    if not os.path.isfile(path):
        return None
    d = np.load(path, allow_pickle=True)
    return d['proteins'].astype(float)


def _load_metrics():
    """Load noise metrics. Returns DataFrame or None."""
    path = os.path.join(PROJECT, 'results', 'phase3', 'noise_metrics.csv')
    if not os.path.isfile(path):
        return None
    return pd.read_csv(path)


def _fit_gmm(data, max_k=3):
    """Fit GMMs with 1..max_k components, return best model."""
    best_bic = np.inf
    best_model = None
    X = data.reshape(-1, 1)
    for k in range(1, max_k + 1):
        gmm = GaussianMixture(n_components=k, reg_covar=1e-3, n_init=5,
                               random_state=PARAMS['gmm_random_state'])
        gmm.fit(X)
        bic = gmm.bic(X)
        if bic < best_bic:
            best_bic = bic
            best_model = gmm
    return best_model


def _plot_histogram_gmm(ax, data, color, label, metrics_row=None):
    """Plot histogram with GMM overlay on a given axes."""
    if data is None or len(data) == 0:
        ax.text(0.5, 0.5, 'DATA NOT\nAVAILABLE', ha='center', va='center',
                fontsize=14, color='gray', transform=ax.transAxes)
        ax.set_xlim(0, 1)
        return

    bins = min(50, max(10, int(np.sqrt(len(data)))))
    ax.hist(data, bins=bins, density=True, alpha=0.5, color=color,
            edgecolor='white', linewidth=0.5, label='Data')

    # GMM overlay
    gmm = _fit_gmm(data)
    x_range = np.linspace(data.min() - 10, data.max() + 10, 300)
    if gmm is not None:
        log_prob = gmm.score_samples(x_range.reshape(-1, 1))
        ax.plot(x_range, np.exp(log_prob), color=color, linewidth=2,
                label=f'GMM (k={gmm.n_components})')

    # Persister threshold
    if metrics_row is not None and pd.notna(metrics_row.get('persister_threshold')):
        thresh = metrics_row['persister_threshold']
        ax.axvline(thresh, color='red', linestyle='--', linewidth=1.2,
                   alpha=0.8, label='Persister\nthreshold')

    ax.set_xlabel('Protein count')
    ax.set_ylabel('Density')
    ax.set_title(label, fontsize=STYLE['title_size'])
    ax.legend(fontsize=7, frameon=False)
    despine(ax)


def generate_fig2():
    """Generate Figure 2: protein distributions for all conditions."""
    apply_style()
    os.makedirs(OUT_DIR, exist_ok=True)

    # Load data
    cond_data = {c: _load_condition(c) for c in ['A', 'B', 'C', 'D']}
    metrics = _load_metrics()

    fig, axes = plt.subplots(2, 3, figsize=(14, 9))

    # ── Panels A-D: histograms with GMM ───────────────────────────────────
    panel_map = [('A', axes[0, 0]), ('B', axes[0, 1]),
                 ('C', axes[0, 2]), ('D', axes[1, 0])]

    for cond, ax in panel_map:
        metrics_row = None
        if metrics is not None:
            m = metrics[metrics['condition'] == cond]
            if len(m) > 0:
                metrics_row = m.iloc[0].to_dict()
        _plot_histogram_gmm(ax, cond_data[cond], CONDITION_COLORS[cond],
                            CONDITION_LABELS[cond], metrics_row)

    # Panel labels
    labels = ['A', 'B', 'C', 'D', 'E', 'F']
    for lbl, ax in zip(labels[:4], [axes[0, 0], axes[0, 1], axes[0, 2], axes[1, 0]]):
        add_panel_label(ax, lbl)

    # ── Panel E: box plot of CV, Fano, bimodality ─────────────────────────
    ax_e = axes[1, 1]
    add_panel_label(ax_e, 'E')

    if metrics is not None and len(metrics) >= 4:
        conditions = metrics['condition'].values
        x = np.arange(len(conditions))
        width = 0.25

        cvs = metrics['CV'].values
        fanos = metrics['Fano'].values
        bimod = metrics['bimodality_coeff'].values

        bar_colors = [CONDITION_COLORS.get(c, 'gray') for c in conditions]

        ax_e.bar(x - width, cvs, width, label='CV', color=bar_colors, alpha=0.7,
                 edgecolor='black', linewidth=0.5)
        ax_e.bar(x, fanos / fanos.max(), width, label='Fano (norm.)',
                 color=bar_colors, alpha=0.5, edgecolor='black', linewidth=0.5,
                 hatch='//')
        ax_e.bar(x + width, bimod, width, label='Bimod. coeff.',
                 color=bar_colors, alpha=0.3, edgecolor='black', linewidth=0.5,
                 hatch='xx')

        ax_e.set_xticks(x)
        ax_e.set_xticklabels([CONDITION_LABELS.get(c, c) for c in conditions],
                              fontsize=8, rotation=20, ha='right')
        ax_e.set_ylabel('Value (normalized)')
        ax_e.set_title('Noise metrics comparison', fontsize=STYLE['title_size'])
        ax_e.legend(fontsize=7, frameon=False)
    else:
        ax_e.text(0.5, 0.5, 'DATA NOT\nAVAILABLE', ha='center', va='center',
                  fontsize=14, color='gray', transform=ax_e.transAxes)
    despine(ax_e)

    # ── Panel F: overlay of all conditions ────────────────────────────────
    ax_f = axes[1, 2]
    add_panel_label(ax_f, 'F')

    has_any = False
    for cond in ['A', 'B', 'C', 'D']:
        data = cond_data[cond]
        if data is not None and len(data) > 2:
            has_any = True
            # KDE-like smooth histogram via fine bins
            bins = min(80, max(15, int(np.sqrt(len(data))) * 2))
            ax_f.hist(data, bins=bins, density=True, alpha=0.35,
                      color=CONDITION_COLORS[cond], label=CONDITION_LABELS[cond],
                      edgecolor='none')
            # Overlay GMM density line
            gmm = _fit_gmm(data)
            if gmm is not None:
                x_range = np.linspace(max(0, data.min() - 20), data.max() + 20, 300)
                log_prob = gmm.score_samples(x_range.reshape(-1, 1))
                ax_f.plot(x_range, np.exp(log_prob), color=CONDITION_COLORS[cond],
                          linewidth=1.8)

    if not has_any:
        ax_f.text(0.5, 0.5, 'DATA NOT\nAVAILABLE', ha='center', va='center',
                  fontsize=14, color='gray', transform=ax_f.transAxes)

    ax_f.set_xlabel('Protein count')
    ax_f.set_ylabel('Density')
    ax_f.set_title('All conditions overlay', fontsize=STYLE['title_size'])
    ax_f.legend(fontsize=7, frameon=False, loc='upper right')
    despine(ax_f)

    fig.suptitle('Protein expression distributions across regulatory conditions',
                 fontsize=STYLE['title_size'] + 1, fontweight='bold', y=1.01)
    fig.tight_layout()
    fig.savefig(OUT_FILE, dpi=STYLE['dpi'])
    plt.close(fig)

    return {
        'status': 'success',
        'conditions_plotted': [c for c, d in cond_data.items() if d is not None],
        'output': OUT_FILE,
    }


# ── Self-tests ────────────────────────────────────────────────────────────────
if __name__ == '__main__':
    n_pass = 0
    n_fail = 0

    result = generate_fig2()

    # SANITY: Output file exists and is valid PNG
    try:
        assert os.path.isfile(OUT_FILE), f"Output not found: {OUT_FILE}"
        fsize = os.path.getsize(OUT_FILE)
        assert fsize > 5000, f"File too small: {fsize} bytes"
        with open(OUT_FILE, 'rb') as f:
            assert f.read(4) == b'\x89PNG', "Not a valid PNG"
        print(f"SANITY PASS: fig2 saved ({fsize:,} bytes)")
        n_pass += 1
    except AssertionError as e:
        print(f"SANITY FAIL: {e}")
        n_fail += 1

    # SANITY: Result dict
    try:
        assert result['status'] == 'success'
        assert len(result['conditions_plotted']) > 0
        print(f"SANITY PASS: Plotted conditions {result['conditions_plotted']}")
        n_pass += 1
    except (AssertionError, KeyError) as e:
        print(f"SANITY FAIL: {e}")
        n_fail += 1

    # SCIENTIFIC: All 4 conditions present
    if set(result.get('conditions_plotted', [])) == {'A', 'B', 'C', 'D'}:
        print("SCIENTIFIC PASS: All 4 conditions plotted")
    else:
        print(f"SCIENTIFIC WARNING: Only {result.get('conditions_plotted')} plotted, expected A-D")

    print(f"\n{'='*50}")
    print(f"Sanity checks: {n_pass} PASS, {n_fail} FAIL")
    if n_fail > 0:
        sys.exit(1)
    else:
        print("All checks passed.")
        sys.exit(0)

generate = generate_fig2
