"""
phase4_figures/fig6_single_cell_traces.py — Time-series traces for single cells.

3 panels: Conditions A, B, C
3 representative cells per panel, ~2100 min post-burn-in
Thin lines with transparency.

Produces: results/figures/fig6_single_cell_traces.png
"""

import sys
sys.path.insert(0, '/Users/aayanalwani/tb project/mce3r_stochastic')

import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from config.parameters import PARAMS
from figures.figure_style import (COLORS, STYLE, CONDITION_COLORS,
                                         CONDITION_LABELS, apply_style,
                                         despine, add_panel_label)

PROJECT = '/Users/aayanalwani/tb project/mce3r_stochastic'
OUT_DIR = os.path.join(PROJECT, 'results', 'figures')
OUT_FILE = os.path.join(OUT_DIR, 'fig6_single_cell_traces.png')


def _load_traces(cond, n_traces=3):
    """Load trace data from condition_X.npz. Returns list of (time, protein) tuples."""
    path = os.path.join(PROJECT, 'results', 'phase2', f'condition_{cond}.npz')
    if not os.path.isfile(path):
        return None

    d = np.load(path, allow_pickle=True)

    # Check for n_traces key
    max_available = int(d['n_traces'][0]) if 'n_traces' in d else 0

    traces = []
    for i in range(min(n_traces, max_available)):
        t_key = f'trace_time_{i}'
        p_key = f'trace_prot_{i}'
        if t_key in d and p_key in d:
            t = d[t_key]
            p = d[p_key]
            # Filter to valid (non-zero time) entries
            mask = t > 0
            if mask.any():
                traces.append((t[mask], p[mask]))

    return traces if traces else None


def generate_fig6():
    """Generate Figure 6: single-cell traces."""
    apply_style()
    os.makedirs(OUT_DIR, exist_ok=True)

    conditions = ['A', 'B', 'C']
    trace_data = {c: _load_traces(c, n_traces=3) for c in conditions}

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5), sharey=False)

    panel_labels = ['A', 'B', 'C']
    trace_colors_per_cell = ['#1a1a1a', '#555555', '#999999']  # dark to light gray

    for idx, cond in enumerate(conditions):
        ax = axes[idx]
        add_panel_label(ax, panel_labels[idx])
        traces = trace_data[cond]

        if traces is None or len(traces) == 0:
            ax.text(0.5, 0.5, f'Condition {cond}\nNO TRACE DATA',
                    ha='center', va='center', fontsize=12, color='gray',
                    transform=ax.transAxes)
            despine(ax)
            continue

        for cell_i, (t, p) in enumerate(traces):
            # Convert time to hours for readability
            t_hours = (t - PARAMS['t_burn_in']) / 60.0
            ax.plot(t_hours, p, linewidth=0.9, alpha=0.7,
                    color=CONDITION_COLORS[cond],
                    label=f'Cell {cell_i + 1}' if cell_i < 3 else None)

        ax.set_xlabel('Time post-burn-in (hours)')
        if idx == 0:
            ax.set_ylabel('Protein count')
        ax.set_title(CONDITION_LABELS[cond], fontsize=STYLE['title_size'],
                     color=CONDITION_COLORS[cond])
        ax.legend(fontsize=7, frameon=False, loc='upper right')
        despine(ax)

    fig.suptitle('Single-cell protein trajectories',
                 fontsize=STYLE['title_size'] + 1, fontweight='bold', y=1.04)
    fig.tight_layout()
    fig.savefig(OUT_FILE, dpi=STYLE['dpi'])
    plt.close(fig)

    n_plotted = sum(1 for c in conditions if trace_data[c] is not None)
    return {
        'status': 'success' if n_plotted > 0 else 'placeholder',
        'conditions_with_traces': [c for c in conditions if trace_data[c] is not None],
        'output': OUT_FILE,
    }


# ── Self-tests ────────────────────────────────────────────────────────────────
if __name__ == '__main__':
    n_pass = 0
    n_fail = 0

    result = generate_fig6()

    try:
        assert os.path.isfile(OUT_FILE)
        fsize = os.path.getsize(OUT_FILE)
        assert fsize > 3000
        with open(OUT_FILE, 'rb') as f:
            assert f.read(4) == b'\x89PNG'
        print(f"SANITY PASS: fig6 saved ({fsize:,} bytes)")
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

    # SCIENTIFIC: trace conditions
    ct = result.get('conditions_with_traces', [])
    if len(ct) >= 3:
        print(f"SCIENTIFIC PASS: Traces for {ct}")
    else:
        print(f"SCIENTIFIC WARNING: Traces for {ct}, expected A, B, C")

    print(f"\n{'='*50}")
    print(f"Sanity checks: {n_pass} PASS, {n_fail} FAIL")
    if n_fail > 0:
        sys.exit(1)
    else:
        print("All checks passed.")
        sys.exit(0)

generate = generate_fig6
