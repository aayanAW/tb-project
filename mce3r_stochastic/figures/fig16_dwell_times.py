"""Figure 16: Operator state dwell times and occupancy comparison."""

import sys, os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from figures.figure_style import COLORS, STYLE, CONDITION_COLORS, apply_style, despine, add_panel_label

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CSV_PATH = os.path.join(BASE, 'results', 'phase8', 'dwell_times.csv')
OUT_DIR = os.path.join(BASE, 'results', 'figures')
OUT_PATH = os.path.join(OUT_DIR, 'fig16_dwell_times.png')

STATE_COLORS = ['#22C55E', '#3B82F6', '#F59E0B', '#EF4444']
STATE_LABELS = ['Unbound', 'Strong only', 'Weak only', 'Both bound']


def generate():
    apply_style()
    os.makedirs(OUT_DIR, exist_ok=True)

    df = pd.read_csv(CSV_PATH)

    fig, axes = plt.subplots(1, 3, figsize=(14, 4.5))

    conditions = df['condition'].tolist()
    labels_map = {'A': 'Asymmetric', 'B': 'Symmetric', 'C': 'Single-site'}
    cond_colors = {
        'A': CONDITION_COLORS.get('A', COLORS['asymmetric']),
        'B': CONDITION_COLORS.get('B', COLORS['symmetric']),
        'C': CONDITION_COLORS.get('C', COLORS['single_site']),
    }

    # Panel A: Dwell times grouped bar chart
    ax = axes[0]
    x = np.arange(len(conditions))
    width = 0.18
    for s in range(4):
        col = f'dwell_state{s}_min'
        vals = df[col].values
        # Cap very small values for visibility
        vals_plot = np.clip(vals, 0.001, None)
        ax.bar(x + s * width, vals_plot, width, color=STATE_COLORS[s],
               label=STATE_LABELS[s], edgecolor='white')

    ax.set_xticks(x + 1.5 * width)
    ax.set_xticklabels([labels_map.get(c, c) for c in conditions], fontsize=STYLE['tick_size'])
    ax.set_ylabel('Mean dwell time (min)', fontsize=STYLE['label_size'])
    ax.set_yscale('log')
    ax.legend(fontsize=7, loc='upper right', ncol=2)
    despine(ax)
    add_panel_label(ax, 'A')

    # Panel B: Stacked bar of occupancy
    ax = axes[1]
    bottom = np.zeros(len(conditions))
    for s in range(4):
        col = f'P_state{s}'
        vals = df[col].values
        ax.bar(x, vals, 0.5, bottom=bottom, color=STATE_COLORS[s], label=STATE_LABELS[s])
        bottom += vals

    ax.set_xticks(x)
    ax.set_xticklabels([labels_map.get(c, c) for c in conditions], fontsize=STYLE['tick_size'])
    ax.set_ylabel('Steady-state probability', fontsize=STYLE['label_size'])
    ax.set_ylim(0, 1.05)
    ax.legend(fontsize=7, loc='upper left', ncol=2)
    despine(ax)
    add_panel_label(ax, 'B')

    # Panel C: Derepression time comparison
    ax = axes[2]
    t_derep = df['t_derepression_min'].values
    bars = ax.bar(x, t_derep, 0.5, color=[cond_colors.get(c, 'gray') for c in conditions],
                  edgecolor='white')
    ax.set_xticks(x)
    ax.set_xticklabels([labels_map.get(c, c) for c in conditions], fontsize=STYLE['tick_size'])
    ax.set_ylabel('Derepression time (min)', fontsize=STYLE['label_size'])
    for bar, val in zip(bars, t_derep):
        if val > 0.01:
            ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.05,
                    f'{val:.2f}', ha='center', va='bottom', fontsize=STYLE['tick_size'])
    despine(ax)
    add_panel_label(ax, 'C')

    fig.suptitle('Figure 16: Operator state dwell times and occupancy',
                 fontsize=STYLE['title_size'], y=1.02)
    plt.tight_layout()
    fig.savefig(OUT_PATH, dpi=STYLE['dpi'], bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved: {OUT_PATH}")
    return OUT_PATH


if __name__ == '__main__':
    generate()
