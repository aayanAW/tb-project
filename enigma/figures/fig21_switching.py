"""Figure 21: State-switching comparison across 4 architectures (Aim 2)."""
import sys, os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from figures.figure_style import COLORS, STYLE, apply_style, despine, add_panel_label

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT_DIR = os.path.join(BASE, 'results', 'figures')

ARCH_COLORS = {
    'wt_asymmetric': COLORS['asymmetric'],
    'symmetric_equal': COLORS['symmetric'],
    'asymmetric_mean_matched': '#8B5CF6',
    'symmetric_variance_matched': '#EF4444',
}
ARCH_LABELS = {
    'wt_asymmetric': 'WT asymmetric',
    'symmetric_equal': 'Symmetric equal',
    'asymmetric_mean_matched': 'Asym (reduced ratio)',
    'symmetric_variance_matched': 'Sym (variance-matched)',
}

def generate():
    apply_style()
    os.makedirs(OUT_DIR, exist_ok=True)
    df = pd.read_csv(os.path.join(BASE, 'results', 'phase_switching', 'switching_data.csv'))
    low = df[df['threshold_type'] == 'low']

    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(15, 5))

    archs = list(ARCH_LABELS.keys())
    x = np.arange(len(archs))

    # Panel A: Mean dwell time
    dwells = [low[low['architecture'] == a]['mean_dwell_min'].values[0] for a in archs]
    colors = [ARCH_COLORS[a] for a in archs]
    bars = ax1.bar(x, dwells, color=colors, width=0.6)
    ax1.set_xticks(x)
    ax1.set_xticklabels([ARCH_LABELS[a] for a in archs], fontsize=8, rotation=20, ha='right')
    ax1.set_ylabel('Mean dwell time (min)', fontsize=STYLE['label_size'])
    ax1.set_yscale('log')
    despine(ax1)
    add_panel_label(ax1, 'A')

    # Panel B: Tail probability (>6hr excursions)
    tails = [low[low['architecture'] == a]['tail_gt_6hr'].values[0] for a in archs]
    bars = ax2.bar(x, tails, color=colors, width=0.6)
    ax2.set_xticks(x)
    ax2.set_xticklabels([ARCH_LABELS[a] for a in archs], fontsize=8, rotation=20, ha='right')
    ax2.set_ylabel('P(excursion > 6 hr)', fontsize=STYLE['label_size'])
    despine(ax2)
    add_panel_label(ax2, 'B')

    # Panel C: Survival under 3 killing models (12hr exposure)
    models = ['survival_hardthresh_12hr', 'survival_hazard_12hr', 'survival_dwell_12hr']
    model_labels = ['Hard threshold', 'Probabilistic hazard', 'Dwell-time']
    width = 0.2
    for j, (model, mlabel) in enumerate(zip(models, model_labels)):
        vals = [low[low['architecture'] == a][model].values[0] for a in archs]
        ax3.bar(x + j * width, vals, width, label=mlabel, alpha=0.8)
    ax3.set_xticks(x + width)
    ax3.set_xticklabels([ARCH_LABELS[a] for a in archs], fontsize=8, rotation=20, ha='right')
    ax3.set_ylabel('Survival fraction (12hr)', fontsize=STYLE['label_size'])
    ax3.legend(fontsize=7)
    despine(ax3)
    add_panel_label(ax3, 'C')

    fig.suptitle('Figure 21: State-switching logic across operator architectures',
                 fontsize=STYLE['title_size'], y=1.02)
    plt.tight_layout()
    path = os.path.join(OUT_DIR, 'fig21_switching.png')
    fig.savefig(path, dpi=STYLE['dpi'], bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved: {path}")
    return path

if __name__ == '__main__':
    generate()
