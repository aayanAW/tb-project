"""Figure 15: Autocorrelation functions for conditions A, B, C."""

import sys, os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from figures.figure_style import COLORS, STYLE, CONDITION_COLORS, apply_style, despine, add_panel_label

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_PATH = os.path.join(BASE, 'results', 'phase8', 'acf_data.npz')
CSV_PATH = os.path.join(BASE, 'results', 'phase8', 'autocorrelation.csv')
OUT_DIR = os.path.join(BASE, 'results', 'figures')
OUT_PATH = os.path.join(OUT_DIR, 'fig15_autocorrelation.png')


def generate():
    apply_style()
    os.makedirs(OUT_DIR, exist_ok=True)

    data = np.load(DATA_PATH)
    import pandas as pd
    df = pd.read_csv(CSV_PATH)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4.5))

    # Panel A: ACF curves
    conditions = ['A', 'B', 'C']
    labels = {'A': 'Asymmetric (Mce3R)', 'B': 'Symmetric control', 'C': 'Single-site control'}
    colors = {
        'A': CONDITION_COLORS.get('A', COLORS['asymmetric']),
        'B': CONDITION_COLORS.get('B', COLORS['symmetric']),
        'C': CONDITION_COLORS.get('C', COLORS['single_site']),
    }

    for cond in conditions:
        lags = data[f'lags_{cond}']
        acf = data[f'mean_acf_{cond}']
        std = data[f'std_acf_{cond}']
        lags_hr = lags / 60.0

        ax1.plot(lags_hr, acf, color=colors[cond], label=labels[cond], linewidth=STYLE['linewidth'])
        ax1.fill_between(lags_hr, acf - std, acf + std, color=colors[cond], alpha=0.15)

    ax1.axhline(1.0 / np.e, color='gray', linestyle='--', linewidth=0.8, alpha=0.7)
    ax1.text(0.5, 1.0 / np.e + 0.03, r'$1/e$', color='gray', fontsize=STYLE['tick_size'])
    ax1.set_xlabel('Lag (hours)', fontsize=STYLE['label_size'])
    ax1.set_ylabel(r'$C(\tau)$', fontsize=STYLE['label_size'])
    ax1.set_xlim(0, 80)
    ax1.set_ylim(-0.1, 1.05)
    ax1.legend(fontsize=STYLE['legend_size'], loc='upper right')
    despine(ax1)
    add_panel_label(ax1, 'A')

    # Panel B: tau_c bar chart
    conds = df['condition'].tolist()
    tau_cs = df['tau_c_mean_acf_hr'].tolist()
    bar_colors = [colors.get(c, 'gray') for c in conds]

    bars = ax2.bar(range(len(conds)), tau_cs, color=bar_colors, width=0.6, edgecolor='white')
    ax2.set_xticks(range(len(conds)))
    ax2.set_xticklabels([labels.get(c, c) for c in conds], fontsize=STYLE['tick_size'], rotation=15, ha='right')
    ax2.set_ylabel(r'$\tau_c$ (hours)', fontsize=STYLE['label_size'])

    for i, (bar, val) in enumerate(zip(bars, tau_cs)):
        ax2.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.3,
                 f'{val:.1f}h', ha='center', va='bottom', fontsize=STYLE['tick_size'])

    despine(ax2)
    add_panel_label(ax2, 'B')

    fig.suptitle('Figure 15: Autocorrelation of protein expression noise',
                 fontsize=STYLE['title_size'], y=1.02)
    plt.tight_layout()
    fig.savefig(OUT_PATH, dpi=STYLE['dpi'], bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved: {OUT_PATH}")
    return OUT_PATH


if __name__ == '__main__':
    generate()
