"""Figure 19: Summary dashboard — key results across all phases."""

import sys, os
import numpy as np
import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from figures.figure_style import COLORS, STYLE, CONDITION_COLORS, apply_style, despine

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT_DIR = os.path.join(BASE, 'results', 'figures')
OUT_PATH = os.path.join(OUT_DIR, 'fig19_summary_dashboard.png')


def _load_json(path):
    if os.path.exists(path):
        with open(path) as f:
            return json.load(f)
    return {}


def generate():
    apply_style()
    os.makedirs(OUT_DIR, exist_ok=True)

    # Load summaries
    p8 = _load_json(os.path.join(BASE, 'results', 'phase8', 'phase8_summary.json'))
    p9 = _load_json(os.path.join(BASE, 'results', 'phase9', 'phase9_summary.json'))

    fig, axes = plt.subplots(2, 3, figsize=(16, 9))

    # Panel A: CV comparison across conditions
    ax = axes[0, 0]
    import pandas as pd
    noise_path = os.path.join(BASE, 'results', 'phase3', 'noise_metrics.csv')
    if os.path.exists(noise_path):
        nm = pd.read_csv(noise_path)
        conds = nm['condition'].tolist()
        cvs = nm['CV'].tolist()
        ccolors = [CONDITION_COLORS.get(c, 'gray') for c in conds]
        bars = ax.bar(range(len(conds)), cvs, color=ccolors, width=0.6)
        ax.set_xticks(range(len(conds)))
        ax.set_xticklabels(['Asym', 'Sym', 'Single', 'Unreg'], fontsize=STYLE['tick_size'])
        ax.set_ylabel('CV', fontsize=STYLE['label_size'])
        ax.set_title('Noise by architecture', fontsize=STYLE['label_size'])
        for bar, val in zip(bars, cvs):
            ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.002,
                    f'{val:.3f}', ha='center', va='bottom', fontsize=8)
    despine(ax)

    # Panel B: tau_c comparison
    ax = axes[0, 1]
    tc_a = p8.get('tau_c_A_hr', 0)
    tc_b = p8.get('tau_c_B_hr', 0)
    tc_labels = ['Asymmetric', 'Symmetric']
    tc_vals = [tc_a, tc_b]
    tc_colors = [COLORS['asymmetric'], COLORS['symmetric']]
    bars = ax.bar(range(2), tc_vals, color=tc_colors, width=0.5)
    ax.set_xticks(range(2))
    ax.set_xticklabels(tc_labels, fontsize=STYLE['tick_size'])
    ax.set_ylabel(r'$\tau_c$ (hours)', fontsize=STYLE['label_size'])
    ax.set_title('Noise memory', fontsize=STYLE['label_size'])
    for bar, val in zip(bars, tc_vals):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.2,
                f'{val:.1f}h', ha='center', va='bottom', fontsize=9)
    despine(ax)

    # Panel C: IC50 visualization
    ax = axes[0, 2]
    ic50 = p9.get('ic50_of_noise', None)
    delta_cv = p9.get('delta_cv', 0)
    cv_asym = p9.get('cv_asymmetric', 0)
    cv_sym = p9.get('cv_symmetric', 0)

    sweep_path = os.path.join(BASE, 'results', 'phase9', 'symmetrization_sweep.csv')
    if os.path.exists(sweep_path):
        sw = pd.read_csv(sweep_path)
        ax.plot(sw['frac_symmetric'] * 100, sw['cv'], 'o-', color=COLORS['asymmetric'],
                markersize=4, linewidth=STYLE['linewidth'])
        if ic50 is not None:
            ax.axvline(ic50 * 100, color=COLORS['mce3r_marker'], linestyle='--', linewidth=1.5)
            ax.text(ic50 * 100 + 2, cv_asym - 0.003, f'IC50={ic50*100:.0f}%',
                    color=COLORS['mce3r_marker'], fontsize=9)
    ax.set_xlabel('Symmetrization (%)', fontsize=STYLE['label_size'])
    ax.set_ylabel('CV', fontsize=STYLE['label_size'])
    ax.set_title('Dose-response', fontsize=STYLE['label_size'])
    despine(ax)

    # Panel D: Persister fraction comparison
    ax = axes[1, 0]
    pf_a = p9.get('persister_fraction_asymmetric', 0)
    pf_s = p9.get('persister_fraction_symmetric', 0)
    bars = ax.bar([0, 1], [pf_a * 100, pf_s * 100],
                  color=[COLORS['asymmetric'], COLORS['symmetric']], width=0.5)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(['Asymmetric', 'Symmetric'], fontsize=STYLE['tick_size'])
    ax.set_ylabel('Persister fraction (%)', fontsize=STYLE['label_size'])
    ax.set_title('Persistence phenotype', fontsize=STYLE['label_size'])
    fold = p9.get('persister_fold_reduction', 0)
    if fold:
        ax.text(0.5, max(pf_a, pf_s) * 100 * 0.8,
                f'{fold:.0f}x reduction', ha='center', fontsize=10, fontweight='bold')
    despine(ax)

    # Panel E: Key numbers summary box
    ax = axes[1, 1]
    ax.axis('off')
    text_lines = [
        f"Delta-CV: {delta_cv:.4f} ({delta_cv/cv_sym*100:.1f}% increase)" if cv_sym > 0 else "",
        f"IC50 of noise: {ic50*100:.0f}% symmetrization" if ic50 else "IC50: N/A",
        f"tau_c(A): {tc_a:.1f} hr  |  tau_c(B): {tc_b:.1f} hr",
        f"tau_c ratio: {tc_a/tc_b:.2f}x" if tc_b > 0 else "",
        f"Persister fold-reduction: {fold:.0f}x" if fold else "",
    ]
    y_pos = 0.9
    for line in text_lines:
        if line:
            ax.text(0.1, y_pos, line, fontsize=11, transform=ax.transAxes,
                    fontfamily='monospace', verticalalignment='top')
            y_pos -= 0.18

    ax.set_title('Key metrics', fontsize=STYLE['label_size'])
    box = FancyBboxPatch((0.05, 0.05), 0.9, 0.9, boxstyle="round,pad=0.05",
                          facecolor='#F3F4F6', edgecolor='#D1D5DB',
                          transform=ax.transAxes, linewidth=1)
    ax.add_patch(box)

    # Panel F: Phase completion status
    ax = axes[1, 2]
    ax.axis('off')
    phases = [
        ('Phase 1: Motif Discovery', True),
        ('Phase 2: Gillespie SSA', True),
        ('Phase 3: Statistical Analysis', True),
        ('Phase 4: Figures 1-8', True),
        ('Phase 5: Thermodynamic Model', True),
        ('Phase 6: Environmental', True),
        ('Phase 7: Figures 9-14', True),
        ('Phase 8: Temporal Dynamics', True),
        ('Phase 9: Noise Quenching', True),
        ('Phase 10: Figures 15-19', True),
    ]
    y_pos = 0.95
    for name, done in phases:
        marker = '\u2713' if done else '\u2717'
        color = '#22C55E' if done else '#EF4444'
        ax.text(0.05, y_pos, marker, fontsize=12, color=color,
                transform=ax.transAxes, fontweight='bold')
        ax.text(0.12, y_pos, name, fontsize=9, transform=ax.transAxes)
        y_pos -= 0.095

    ax.set_title('Pipeline status', fontsize=STYLE['label_size'])

    fig.suptitle('Figure 19: ENIGMA Project Summary Dashboard',
                 fontsize=14, fontweight='bold', y=1.02)
    plt.tight_layout()
    fig.savefig(OUT_PATH, dpi=STYLE['dpi'], bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved: {OUT_PATH}")
    return OUT_PATH


if __name__ == '__main__':
    generate()
