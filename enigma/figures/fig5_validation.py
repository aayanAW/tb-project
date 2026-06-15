"""
phase4_figures/fig5_validation.py — Experimental comparison & controls.

(A) Scatter: predicted vs published fold-change
(B) Bar: persister fraction WT vs unregulated
(C) Shuffle test histogram
(D) Symmetric TetR vs asymmetric distribution

Produces: results/figures/fig5_validation.png
"""

import sys
sys.path.insert(0, '/Users/aayanalwani/tb project/mce3r_stochastic')

import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from config.parameters import PARAMS
from figures.figure_style import (COLORS, STYLE, CONDITION_COLORS,
                                         CONDITION_LABELS, apply_style,
                                         despine, add_panel_label)

PROJECT = '/Users/aayanalwani/tb project/mce3r_stochastic'
OUT_DIR = os.path.join(PROJECT, 'results', 'figures')
OUT_FILE = os.path.join(OUT_DIR, 'fig5_validation.png')

EXP_FILE = os.path.join(PROJECT, 'results', 'phase3', 'experimental_comparison.csv')
METRICS_FILE = os.path.join(PROJECT, 'results', 'phase3', 'noise_metrics.csv')
NEG_FILE = os.path.join(PROJECT, 'results', 'phase3', 'negative_controls.csv')


def _safe_load(path):
    if os.path.isfile(path):
        return pd.read_csv(path)
    return None


def generate_fig5():
    """Generate Figure 5: validation panels."""
    apply_style()
    os.makedirs(OUT_DIR, exist_ok=True)

    exp = _safe_load(EXP_FILE)
    metrics = _safe_load(METRICS_FILE)
    neg = _safe_load(NEG_FILE)

    fig, axes = plt.subplots(2, 2, figsize=(11, 9))
    ax_a, ax_b, ax_c, ax_d = axes[0, 0], axes[0, 1], axes[1, 0], axes[1, 1]

    # ── Panel A: Predicted vs published fold-change ───────────────────────
    add_panel_label(ax_a, 'A')
    if exp is not None and 'fold_change_D_over_A' in exp.columns:
        pred = exp['fold_change_D_over_A'].values
        pub = exp['published_fold_change'].values

        ax_a.scatter(pub, pred, s=120, color=COLORS['asymmetric'], edgecolors='black',
                     linewidth=0.8, zorder=3)

        # Identity line
        all_vals = np.concatenate([pub, pred])
        lo, hi = all_vals.min() * 0.8, all_vals.max() * 1.2
        ax_a.plot([lo, hi], [lo, hi], '--', color='gray', linewidth=1, alpha=0.6,
                  label='y = x')

        # Annotate ratio
        if 'fold_change_ratio' in exp.columns:
            ratio = exp['fold_change_ratio'].values[0]
            ax_a.annotate(f'Ratio: {ratio:.2f}',
                         xy=(pub[0], pred[0]),
                         xytext=(pub[0] + 1, pred[0] + 1),
                         fontsize=9, arrowprops=dict(arrowstyle='->', color='gray'))

        ax_a.set_xlabel('Published fold-change')
        ax_a.set_ylabel('Predicted fold-change')
        ax_a.set_title('Model vs experimental fold-change', fontsize=STYLE['title_size'])
        ax_a.legend(fontsize=8, frameon=False)
    else:
        ax_a.text(0.5, 0.5, 'DATA NOT\nAVAILABLE', ha='center', va='center',
                  fontsize=14, color='gray', transform=ax_a.transAxes)
    despine(ax_a)

    # ── Panel B: Persister fraction WT vs unregulated ─────────────────────
    add_panel_label(ax_b, 'B')
    if metrics is not None and 'persister_fraction' in metrics.columns:
        conds = ['A', 'D']
        fracs = []
        colors_bar = []
        labels_bar = []
        for c in conds:
            row = metrics[metrics['condition'] == c]
            if len(row) > 0:
                fracs.append(row['persister_fraction'].values[0])
                colors_bar.append(CONDITION_COLORS[c])
                labels_bar.append(CONDITION_LABELS[c])

        if fracs:
            x = np.arange(len(fracs))
            ax_b.bar(x, fracs, color=colors_bar, edgecolor='black', linewidth=0.6,
                     width=0.5)
            ax_b.set_xticks(x)
            ax_b.set_xticklabels(labels_bar, fontsize=9)
            ax_b.set_ylabel('Persister fraction')
            ax_b.set_title('Persister fraction: WT vs unregulated',
                          fontsize=STYLE['title_size'])
            ax_b.set_ylim(0, 1.1)
        else:
            ax_b.text(0.5, 0.5, 'NO DATA', ha='center', va='center',
                      fontsize=14, color='gray', transform=ax_b.transAxes)
    else:
        ax_b.text(0.5, 0.5, 'DATA NOT\nAVAILABLE', ha='center', va='center',
                  fontsize=14, color='gray', transform=ax_b.transAxes)
    despine(ax_b)

    # ── Panel C: Shuffle test histogram ───────────────────────────────────
    add_panel_label(ax_c, 'C')
    if neg is not None:
        shuffle_rows = neg[neg['test'] == 'shuffle']
        if len(shuffle_rows) > 0:
            shuffle_vals = shuffle_rows['value'].values
            ax_c.hist(shuffle_vals, bins=max(5, len(shuffle_vals) // 2),
                      color=COLORS['unregulated'], edgecolor='white',
                      linewidth=0.5, alpha=0.7, label='Shuffle scores')

            # Mark p-value threshold
            alpha_thresh = PARAMS.get('shuffle_alpha', 0.01)
            ax_c.axvline(alpha_thresh, color='red', linestyle='--', linewidth=1.5,
                         label=f'alpha = {alpha_thresh}')

            is_mock = any('mock' in str(v) for v in shuffle_rows['note'].values if pd.notna(v))
            if is_mock:
                ax_c.text(0.05, 0.95, '(mock data)', transform=ax_c.transAxes,
                          fontsize=8, color='red', va='top')

            ax_c.set_xlabel('Shuffle test score')
            ax_c.set_ylabel('Count')
            ax_c.set_title('Negative control: binding site shuffle',
                          fontsize=STYLE['title_size'])
            ax_c.legend(fontsize=8, frameon=False)
        else:
            ax_c.text(0.5, 0.5, 'NO SHUFFLE\nDATA', ha='center', va='center',
                      fontsize=14, color='gray', transform=ax_c.transAxes)
    else:
        ax_c.text(0.5, 0.5, 'DATA NOT\nAVAILABLE', ha='center', va='center',
                  fontsize=14, color='gray', transform=ax_c.transAxes)
    despine(ax_c)

    # ── Panel D: Symmetric vs asymmetric distributions ────────────────────
    add_panel_label(ax_d, 'D')
    data_a = None
    data_b = None
    for cond, var in [('A', 'data_a'), ('B', 'data_b')]:
        path = os.path.join(PROJECT, 'results', 'phase2', f'condition_{cond}.npz')
        if os.path.isfile(path):
            d = np.load(path, allow_pickle=True)
            if cond == 'A':
                data_a = d['proteins'].astype(float)
            else:
                data_b = d['proteins'].astype(float)

    if data_a is not None and data_b is not None:
        bins = max(10, int(np.sqrt(max(len(data_a), len(data_b)))))
        ax_d.hist(data_a, bins=bins, density=True, alpha=0.5,
                  color=COLORS['asymmetric'], label='Asymmetric (Mce3R)',
                  edgecolor='white', linewidth=0.3)
        ax_d.hist(data_b, bins=bins, density=True, alpha=0.5,
                  color=COLORS['symmetric'], label='Symmetric (TetR-like)',
                  edgecolor='white', linewidth=0.3)
        ax_d.set_xlabel('Protein count')
        ax_d.set_ylabel('Density')
        ax_d.set_title('Asymmetric vs symmetric regulation',
                      fontsize=STYLE['title_size'])
        ax_d.legend(fontsize=8, frameon=False)
    else:
        ax_d.text(0.5, 0.5, 'DATA NOT\nAVAILABLE', ha='center', va='center',
                  fontsize=14, color='gray', transform=ax_d.transAxes)
    despine(ax_d)

    fig.suptitle('Model validation and controls',
                 fontsize=STYLE['title_size'] + 1, fontweight='bold', y=1.02)
    fig.tight_layout()
    fig.savefig(OUT_FILE, dpi=STYLE['dpi'])
    plt.close(fig)

    return {
        'status': 'success',
        'has_exp_comparison': exp is not None,
        'has_neg_controls': neg is not None,
        'output': OUT_FILE,
    }


# ── Self-tests ────────────────────────────────────────────────────────────────
if __name__ == '__main__':
    n_pass = 0
    n_fail = 0

    result = generate_fig5()

    try:
        assert os.path.isfile(OUT_FILE)
        fsize = os.path.getsize(OUT_FILE)
        assert fsize > 5000
        with open(OUT_FILE, 'rb') as f:
            assert f.read(4) == b'\x89PNG'
        print(f"SANITY PASS: fig5 saved ({fsize:,} bytes)")
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

    # SCIENTIFIC checks
    if result.get('has_exp_comparison'):
        print("SCIENTIFIC PASS: Experimental comparison data present")
    else:
        print("SCIENTIFIC WARNING: No experimental comparison data")

    if result.get('has_neg_controls'):
        print("SCIENTIFIC PASS: Negative control data present")
    else:
        print("SCIENTIFIC WARNING: No negative control data")

    print(f"\n{'='*50}")
    print(f"Sanity checks: {n_pass} PASS, {n_fail} FAIL")
    if n_fail > 0:
        sys.exit(1)
    else:
        print("All checks passed.")
        sys.exit(0)

generate = generate_fig5
