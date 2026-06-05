"""Figure 17: Power spectral density of protein expression noise."""

import sys, os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from figures.figure_style import COLORS, STYLE, CONDITION_COLORS, apply_style, despine, add_panel_label

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_PATH = os.path.join(BASE, 'results', 'phase8', 'psd_data.npz')
OUT_DIR = os.path.join(BASE, 'results', 'figures')
OUT_PATH = os.path.join(OUT_DIR, 'fig17_power_spectrum.png')


def generate():
    apply_style()
    os.makedirs(OUT_DIR, exist_ok=True)

    data = np.load(DATA_PATH)

    fig, ax = plt.subplots(1, 1, figsize=(7, 5))

    conditions = ['A', 'B', 'C']
    labels = {'A': 'Asymmetric (Mce3R)', 'B': 'Symmetric control', 'C': 'Single-site control'}
    colors = {
        'A': CONDITION_COLORS.get('A', COLORS['asymmetric']),
        'B': CONDITION_COLORS.get('B', COLORS['symmetric']),
        'C': CONDITION_COLORS.get('C', COLORS['single_site']),
    }

    for cond in conditions:
        freqs = data[f'freqs_{cond}']
        psd = data[f'psd_{cond}']

        # Convert frequency from cycles/min to cycles/hr for readability
        freqs_hr = freqs * 60.0

        # Smooth PSD with log-binning for clarity
        ax.loglog(freqs_hr, psd, color=colors[cond], label=labels[cond],
                  linewidth=STYLE['linewidth'], alpha=0.8)

    # Mark protein half-life frequency
    t_half_prot_hr = 1500.0 / 60.0  # 25 hr
    f_prot = 1.0 / t_half_prot_hr
    ax.axvline(f_prot, color='gray', linestyle='--', linewidth=0.8, alpha=0.6)
    ax.text(f_prot * 1.2, ax.get_ylim()[1] * 0.3, r'$f_{prot}$', color='gray',
            fontsize=STYLE['tick_size'])

    ax.set_xlabel('Frequency (cycles/hour)', fontsize=STYLE['label_size'])
    ax.set_ylabel('Power spectral density (protein$^2$ / (cycles/hr))', fontsize=STYLE['label_size'])
    ax.legend(fontsize=STYLE['legend_size'], loc='upper right')
    despine(ax)

    fig.suptitle('Figure 17: Power spectral density of expression noise',
                 fontsize=STYLE['title_size'], y=1.02)
    plt.tight_layout()
    fig.savefig(OUT_PATH, dpi=STYLE['dpi'], bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved: {OUT_PATH}")
    return OUT_PATH


if __name__ == '__main__':
    generate()
