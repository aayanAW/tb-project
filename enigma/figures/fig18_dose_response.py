"""Figure 18: Noise quenching dose-response curve (symmetrization sweep)."""

import sys, os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from figures.figure_style import COLORS, STYLE, apply_style, despine, add_panel_label

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CSV_PATH = os.path.join(BASE, 'results', 'phase9', 'symmetrization_sweep.csv')
OUT_DIR = os.path.join(BASE, 'results', 'figures')
OUT_PATH = os.path.join(OUT_DIR, 'fig18_dose_response.png')


def generate():
    apply_style()
    os.makedirs(OUT_DIR, exist_ok=True)

    df = pd.read_csv(CSV_PATH)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    frac = df['frac_symmetric'].values
    cv = df['cv'].values
    pf = df['persister_fraction'].values
    kd_ratio = df['Kd_ratio'].values

    # Find IC50
    cv_asym = cv[0]
    cv_sym = cv[-1]
    delta_cv_max = cv_asym - cv_sym
    if 'noise_reduction_frac' in df.columns:
        nrf = df['noise_reduction_frac'].values
    else:
        nrf = 1.0 - ((cv - cv_sym) / delta_cv_max) if delta_cv_max > 0 else np.zeros_like(cv)

    crossings = np.where(np.diff(np.sign(nrf - 0.5)))[0]
    ic50 = np.nan
    if len(crossings) > 0:
        idx = crossings[0]
        f1, f2 = frac[idx], frac[idx + 1]
        n1, n2 = nrf[idx], nrf[idx + 1]
        ic50 = f1 + (0.5 - n1) * (f2 - f1) / (n2 - n1) if n2 != n1 else f1

    # Panel A: CV dose-response
    ax1.plot(frac * 100, cv, 'o-', color=COLORS['asymmetric'], linewidth=STYLE['linewidth'],
             markersize=5)
    ax1.axhline(cv_sym, color='gray', linestyle=':', linewidth=0.8, alpha=0.6)
    ax1.axhline(cv_asym, color='gray', linestyle=':', linewidth=0.8, alpha=0.6)

    # Mark IC50
    if np.isfinite(ic50):
        cv_at_ic50 = cv_sym + delta_cv_max * 0.5
        ax1.axvline(ic50 * 100, color=COLORS['mce3r_marker'], linestyle='--', linewidth=1.2)
        ax1.annotate(f'IC50 = {ic50*100:.0f}%', xy=(ic50 * 100, cv_at_ic50),
                     xytext=(ic50 * 100 + 10, cv_at_ic50 + 0.005),
                     fontsize=STYLE['tick_size'],
                     arrowprops=dict(arrowstyle='->', color=COLORS['mce3r_marker']),
                     color=COLORS['mce3r_marker'])

    # Add secondary x-axis showing Kd ratio
    ax1_top = ax1.twiny()
    ax1_top.set_xlim(ax1.get_xlim())
    ratio_ticks = [0, 25, 50, 75, 100]
    ratio_labels = []
    for rt in ratio_ticks:
        idx_closest = np.argmin(np.abs(frac * 100 - rt))
        ratio_labels.append(f'{kd_ratio[idx_closest]:.1f}')
    ax1_top.set_xticks(ratio_ticks)
    ax1_top.set_xticklabels(ratio_labels, fontsize=STYLE['tick_size'] - 1)
    ax1_top.set_xlabel('Kd ratio (weak/strong)', fontsize=STYLE['tick_size'])

    ax1.set_xlabel('Symmetrization dose (%)', fontsize=STYLE['label_size'])
    ax1.set_ylabel('Coefficient of Variation (CV)', fontsize=STYLE['label_size'])
    despine(ax1)
    add_panel_label(ax1, 'A')

    # Panel B: Persister fraction dose-response
    ax2.plot(frac * 100, pf * 100, 's-', color='#EF4444', linewidth=STYLE['linewidth'],
             markersize=5)

    if np.isfinite(ic50):
        ax2.axvline(ic50 * 100, color=COLORS['mce3r_marker'], linestyle='--', linewidth=1.2)

    ax2.set_xlabel('Symmetrization dose (%)', fontsize=STYLE['label_size'])
    ax2.set_ylabel('Persister fraction (%)', fontsize=STYLE['label_size'])
    despine(ax2)
    add_panel_label(ax2, 'B')

    fig.suptitle('Figure 18: Noise quenching dose-response',
                 fontsize=STYLE['title_size'], y=1.02)
    plt.tight_layout()
    fig.savefig(OUT_PATH, dpi=STYLE['dpi'], bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved: {OUT_PATH}")
    return OUT_PATH


if __name__ == '__main__':
    generate()
