"""Figure 22: Extended phase diagram with repeated antibiotic pulses (Aim 3)."""
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

ENV_LABELS = {
    'baseline': 'Baseline',
    'cholesterol': 'Cholesterol',
    'acidic_pH': 'Acidic pH',
    'oxidative': 'Oxidative',
    'combined_macrophage': 'Macrophage-like',
}

def generate():
    apply_style()
    os.makedirs(OUT_DIR, exist_ok=True)
    df = pd.read_csv(os.path.join(BASE, 'results', 'phase6', 'extended_phase_diagram.csv'))

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5.5))

    # Panel A: Heatmap — ratio x environment, color = multi-round survival
    envs = list(ENV_LABELS.keys())
    ratios = sorted(df['Kd_ratio'].unique())

    matrix = np.zeros((len(envs), len(ratios)))
    for i, env in enumerate(envs):
        for j, ratio in enumerate(ratios):
            subset = df[(df['environment'] == env) & (np.abs(df['Kd_ratio'] - ratio) < 0.1)]
            if len(subset) > 0:
                matrix[i, j] = subset.iloc[0]['multi_round_survival']

    im = ax1.imshow(matrix, cmap='YlOrRd', aspect='auto', vmin=0)
    ax1.set_xticks(range(len(ratios)))
    ax1.set_xticklabels([f'{r:.0f}' if r >= 1 else f'{r:.1f}' for r in ratios],
                        fontsize=STYLE['tick_size'])
    ax1.set_yticks(range(len(envs)))
    ax1.set_yticklabels([ENV_LABELS[e] for e in envs], fontsize=STYLE['tick_size'])
    ax1.set_xlabel('Kd ratio (weak/strong)', fontsize=STYLE['label_size'])

    # Annotate cells
    for i in range(len(envs)):
        for j in range(len(ratios)):
            val = matrix[i, j]
            color = 'white' if val > 0.3 else 'black'
            ax1.text(j, i, f'{val:.2f}', ha='center', va='center',
                    fontsize=7, color=color)

    # Mark WT ratio
    wt_idx = np.argmin(np.abs(np.array(ratios) - 20.4))
    for i in range(len(envs)):
        ax1.plot(wt_idx, i, '*', color=COLORS['mce3r_marker'], markersize=12,
                markeredgecolor='black', markeredgewidth=0.5)

    plt.colorbar(im, ax=ax1, label='Multi-round survival', shrink=0.8)
    ax1.set_title('4-round antibiotic survival', fontsize=STYLE['label_size'])
    for spine in ax1.spines.values():
        spine.set_visible(True)
    add_panel_label(ax1, 'A')

    # Panel B: Growth-survival tradeoff by environment
    for env in envs:
        env_data = df[df['environment'] == env].sort_values('Kd_ratio')
        ax2.plot(env_data['growth_score'], env_data['multi_round_survival'],
                'o-', label=ENV_LABELS[env], markersize=5, linewidth=STYLE['linewidth'])

    # Mark WT
    wt_data = df[np.abs(df['Kd_ratio'] - 20.4) < 1.0]
    if len(wt_data) > 0:
        for _, row in wt_data.iterrows():
            ax2.plot(row['growth_score'], row['multi_round_survival'],
                    '*', color=COLORS['mce3r_marker'], markersize=15,
                    markeredgecolor='black', zorder=5)

    ax2.set_xlabel('Growth score', fontsize=STYLE['label_size'])
    ax2.set_ylabel('Multi-round survival', fontsize=STYLE['label_size'])
    ax2.legend(fontsize=STYLE['legend_size'], loc='best')
    despine(ax2)
    add_panel_label(ax2, 'B')

    fig.suptitle('Figure 22: Environment-specific asymmetry benefit across repeated treatment',
                 fontsize=STYLE['title_size'], y=1.02)
    plt.tight_layout()
    path = os.path.join(OUT_DIR, 'fig22_repeated_pulses.png')
    fig.savefig(path, dpi=STYLE['dpi'], bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved: {path}")
    return path

if __name__ == '__main__':
    generate()
