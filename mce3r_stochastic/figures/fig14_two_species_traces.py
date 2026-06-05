"""
phase7_extended_figures/fig14_two_species_traces.py

Figure 14: Two-species stochastic traces and scatter plot.
(A) Single-cell time traces of Mce3R and target protein for 5 cells.
(B) Scatter plot of Mce3R vs target protein colored by operator state.
"""

import sys
sys.path.insert(0, '/Users/aayanalwani/tb project/mce3r_stochastic')

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import os

from figures.figure_style import COLORS, STYLE, apply_style, despine, add_panel_label

# ── Paths ─────────────────────────────────────────────────────────────────────
BASE = '/Users/aayanalwani/tb project/mce3r_stochastic'
TRACE_FILE = os.path.join(BASE, 'results', 'phase6', 'two_species_results.npz')
SCATTER_FILE = os.path.join(BASE, 'results', 'phase6',
                            'env_condition_asymmetric_host_like_two_species.npz')
OUT_DIR = os.path.join(BASE, 'results', 'figures')
OUT_PATH = os.path.join(OUT_DIR, 'fig14_two_species_traces.png')

# ── Operator-state colors and labels ──────────────────────────────────────────
OP_COLORS = {0: '#22C55E', 1: '#3B82F6', 2: '#F59E0B', 3: '#EF4444'}
OP_LABELS = {0: 'Unbound', 1: 'Strong bound', 2: 'Weak bound', 3: 'Both bound'}


def generate():
    apply_style()
    os.makedirs(OUT_DIR, exist_ok=True)

    # ── Load data ─────────────────────────────────────────────────────────────
    trace_data = np.load(TRACE_FILE, allow_pickle=True)
    scatter_data = np.load(SCATTER_FILE, allow_pickle=True)

    # ── Create figure ─────────────────────────────────────────────────────────
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4.5))

    # ── Panel A: single-cell time traces ──────────────────────────────────────
    cmap = plt.cm.tab10
    for i in range(5):
        t = trace_data[f'trace_time_{i}']
        mce3r = trace_data[f'trace_mce3r_{i}']
        target = trace_data[f'trace_target_{i}']
        n = int(trace_data[f'n_trace_{i}'].flat[0])
        color = cmap(i)
        ax1.plot(t[:n], mce3r[:n], '-', color=color, linewidth=STYLE['linewidth'],
                 alpha=0.85, label=f'Cell {i+1}' if i == 0 else None)
        ax1.plot(t[:n], target[:n], '--', color=color, linewidth=STYLE['linewidth'],
                 alpha=0.85)

    # Add proxy legend entries for solid vs dashed
    from matplotlib.lines import Line2D
    proxy_solid = Line2D([], [], color='k', linestyle='-', linewidth=STYLE['linewidth'],
                         label='Mce3R')
    proxy_dash = Line2D([], [], color='k', linestyle='--', linewidth=STYLE['linewidth'],
                        label='Target')
    cell_proxies = [Line2D([], [], color=cmap(i), linestyle='-',
                           linewidth=STYLE['linewidth'], label=f'Cell {i+1}')
                    for i in range(5)]
    ax1.legend(handles=[proxy_solid, proxy_dash] + cell_proxies,
               loc='upper right', frameon=False, fontsize=STYLE['legend_size'])
    ax1.set_xlabel('Time (min)')
    ax1.set_ylabel('Protein count')
    despine(ax1)
    add_panel_label(ax1, 'A')

    # ── Panel B: scatter plot colored by operator state ───────────────────────
    mce3r_all = scatter_data['mce3r_proteins']
    target_all = scatter_data['target_proteins']
    op_states = scatter_data['op_states']

    for state in sorted(OP_COLORS.keys()):
        mask = op_states == state
        if np.any(mask):
            ax2.scatter(mce3r_all[mask], target_all[mask],
                        c=OP_COLORS[state], s=5, alpha=0.3,
                        label=OP_LABELS[state], edgecolors='none', rasterized=True)

    ax2.set_xlabel('Mce3R protein count')
    ax2.set_ylabel('Target protein count')
    ax2.legend(loc='upper right', frameon=False, fontsize=STYLE['legend_size'],
               markerscale=3)
    despine(ax2)
    add_panel_label(ax2, 'B')

    # ── Save ──────────────────────────────────────────────────────────────────
    plt.tight_layout()
    fig.savefig(OUT_PATH, dpi=STYLE['dpi'], format='png')
    plt.close(fig)
    print(f"Saved: {OUT_PATH}")

    # ── Self-tests ────────────────────────────────────────────────────────────
    n_pass = 0
    n_fail = 0

    # Test 1: file exists
    if os.path.isfile(OUT_PATH):
        print("PASS: Output file exists")
        n_pass += 1
    else:
        print("FAIL: Output file does not exist")
        n_fail += 1

    # Test 2: file size > 10 KB
    if os.path.isfile(OUT_PATH):
        size_kb = os.path.getsize(OUT_PATH) / 1024
        if size_kb > 10:
            print(f"PASS: File size = {size_kb:.1f} KB (> 10 KB)")
            n_pass += 1
        else:
            print(f"FAIL: File size = {size_kb:.1f} KB (<= 10 KB)")
            n_fail += 1

    # Test 3: valid PNG header
    if os.path.isfile(OUT_PATH):
        with open(OUT_PATH, 'rb') as f:
            header = f.read(8)
        if header[:4] == b'\x89PNG':
            print("PASS: Valid PNG header")
            n_pass += 1
        else:
            print("FAIL: Invalid PNG header")
            n_fail += 1

    print(f"\n{'='*50}")
    print(f"Self-tests: {n_pass} PASS, {n_fail} FAIL")
    if n_fail > 0:
        sys.exit(1)


if __name__ == '__main__':
    generate()
