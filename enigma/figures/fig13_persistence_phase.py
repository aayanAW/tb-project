"""
phase7_extended_figures/fig13_persistence_phase.py
Figure 13: Persistence phase diagram.

2D heatmap with X = asymmetry ratio, Y = environment condition,
color = persister fraction.  Marks the native Mce3R operator position
(ratio ~ 20.4) with a star annotation.
"""

import sys
sys.path.insert(0, '/Users/aayanalwani/tb project/mce3r_stochastic')

import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from figures.figure_style import COLORS, STYLE, apply_style, despine
from config.parameters import PARAMS

# ── Paths ────────────────────────────────────────────────────────────────────
BASE = '/Users/aayanalwani/tb project/mce3r_stochastic'
INPUT_CSV = os.path.join(BASE, 'results', 'phase6', 'phase_diagram.csv')
OUT_DIR = os.path.join(BASE, 'results', 'figures')
OUT_PATH = os.path.join(OUT_DIR, 'fig13_persistence_phase.png')

# ── Environment display order & labels ────────────────────────────────────────
ENV_ORDER = ['baseline', 'cholesterol', 'acidic_pH', 'host_like']
ENV_LABELS = {
    'baseline':    'Baseline',
    'cholesterol': 'Cholesterol',
    'acidic_pH':   'Acidic pH',
    'host_like':   'Host-like',
}


def main():
    apply_style()
    os.makedirs(OUT_DIR, exist_ok=True)

    # ── Load data ─────────────────────────────────────────────────────────────
    df = pd.read_csv(INPUT_CSV)

    # If an architecture column exists, filter to asymmetric only
    if 'architecture' in df.columns:
        df = df[df['architecture'] == 'asymmetric'].copy()

    # ── Pivot to 2D matrix  (rows = environment, cols = ratio) ────────────────
    pivot = df.pivot_table(
        index='environment',
        columns='ratio',
        values='persister_fraction',
        aggfunc='mean',
    )

    # Reorder rows to match ENV_ORDER (only those present)
    env_order = [e for e in ENV_ORDER if e in pivot.index]
    pivot = pivot.loc[env_order]

    ratio_vals = pivot.columns.values.astype(float)
    matrix = pivot.values  # shape (n_envs, n_ratios)

    # ── Native operator ratio position ────────────────────────────────────────
    native_ratio = PARAMS['mce3r_actual_ratio']  # 20.4

    # Find the column index closest to native_ratio
    native_col_idx = int(np.argmin(np.abs(ratio_vals - native_ratio)))

    # ── Plot ──────────────────────────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(8, 4))

    # Use pcolormesh for crisp cells; imshow is simpler for categorical axes
    im = ax.imshow(
        matrix,
        aspect='auto',
        cmap='YlOrRd',
        vmin=0,
        vmax=1,
        origin='upper',
        interpolation='nearest',
    )

    # ── Annotate each cell ────────────────────────────────────────────────────
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            val = matrix[i, j]
            # Use white text on dark cells, black on light
            text_color = 'white' if val > 0.65 else 'black'
            ax.text(
                j, i, f'{val:.2f}',
                ha='center', va='center',
                fontsize=STYLE['tick_size'] - 1,
                color=text_color,
                fontweight='bold',
            )

    # ── Mark native operator position with a star ─────────────────────────────
    for i in range(matrix.shape[0]):
        ax.plot(
            native_col_idx, i,
            marker='*',
            markersize=14,
            color=COLORS['mce3r_marker'],
            markeredgecolor='black',
            markeredgewidth=0.8,
            zorder=5,
        )

    # ── Axes labels & ticks ───────────────────────────────────────────────────
    ax.set_xticks(np.arange(len(ratio_vals)))
    ax.set_xticklabels([f'{r:g}' for r in ratio_vals])
    ax.set_xlabel('Asymmetry ratio ($K_{{d,weak}} / K_{{d,strong}}$)')

    y_labels = [ENV_LABELS.get(e, e) for e in env_order]
    ax.set_yticks(np.arange(len(env_order)))
    ax.set_yticklabels(y_labels)
    ax.set_ylabel('Environment')

    ax.set_title('Figure 13: Persistence phase diagram')

    # ── Colorbar ──────────────────────────────────────────────────────────────
    cbar = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.04)
    cbar.set_label('Persister fraction')

    # Re-enable all spines for heatmap (override apply_style defaults)
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_linewidth(STYLE['spine_width'])

    plt.tight_layout()
    fig.savefig(OUT_PATH, dpi=STYLE['dpi'])
    plt.close(fig)
    print(f'Saved: {OUT_PATH}')

    # ── Self-tests ────────────────────────────────────────────────────────────
    n_pass = 0
    n_fail = 0

    # Test 1: file exists
    if os.path.isfile(OUT_PATH):
        print('SELF-TEST PASS: output file exists')
        n_pass += 1
    else:
        print('SELF-TEST FAIL: output file not found')
        n_fail += 1

    # Test 2: file size > 10 KB
    if os.path.isfile(OUT_PATH):
        sz = os.path.getsize(OUT_PATH)
        if sz > 10_000:
            print(f'SELF-TEST PASS: file size = {sz:,} bytes (> 10 KB)')
            n_pass += 1
        else:
            print(f'SELF-TEST FAIL: file size = {sz:,} bytes (expected > 10 KB)')
            n_fail += 1

    # Test 3: valid PNG header
    if os.path.isfile(OUT_PATH):
        with open(OUT_PATH, 'rb') as f:
            header = f.read(8)
        if header[:4] == b'\x89PNG':
            print('SELF-TEST PASS: valid PNG header')
            n_pass += 1
        else:
            print('SELF-TEST FAIL: invalid PNG header')
            n_fail += 1

    print(f'\n{"=" * 50}')
    print(f'Self-tests: {n_pass} PASS, {n_fail} FAIL')
    if n_fail > 0:
        sys.exit(1)
    else:
        print('All checks passed.')
        sys.exit(0)


if __name__ == '__main__':
    main()


generate = main
