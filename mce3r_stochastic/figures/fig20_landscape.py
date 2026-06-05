"""Figure 20: Fitness landscape Pareto front (Aim 1)."""
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

def generate():
    apply_style()
    os.makedirs(OUT_DIR, exist_ok=True)
    df = pd.read_csv(os.path.join(BASE, 'results', 'phase_landscape', 'landscape_data.csv'))

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5.5))

    # Panel A: Growth vs Persistence tradeoff
    wt = df[df['mutation'] == 'WT'].iloc[0]
    strong = df[df['site'] == 'strong_site']
    weak = df[df['site'] == 'weak_site']
    pareto = df[df['is_pareto'] == True]

    ax1.scatter(strong['growth_score'], strong['persister_fraction'],
                c='#3B82F6', alpha=0.5, s=25, label='Strong-site mutants', zorder=2)
    ax1.scatter(weak['growth_score'], weak['persister_fraction'],
                c='#F59E0B', alpha=0.5, s=25, label='Weak-site mutants', zorder=2)

    # Pareto front
    pf = pareto.sort_values('growth_score')
    ax1.plot(pf['growth_score'], pf['persister_fraction'], 'k--', linewidth=1, alpha=0.5, zorder=3)

    # Wild-type
    ax1.scatter([wt['growth_score']], [wt['persister_fraction']],
                c='#EF4444', s=120, marker='*', zorder=5, label='Wild-type Mce3R', edgecolors='black')

    ax1.set_xlabel('Growth score (mean protein / 500)', fontsize=STYLE['label_size'])
    ax1.set_ylabel('Persister fraction', fontsize=STYLE['label_size'])
    ax1.legend(fontsize=STYLE['legend_size'], loc='upper right')
    despine(ax1)
    add_panel_label(ax1, 'A')

    # Panel B: CV vs Kd ratio colored by site
    ax2.scatter(strong['Kd_ratio'], strong['cv'],
                c='#3B82F6', alpha=0.5, s=25, label='Strong-site mutants')
    ax2.scatter(weak['Kd_ratio'], weak['cv'],
                c='#F59E0B', alpha=0.5, s=25, label='Weak-site mutants')
    ax2.scatter([wt['Kd_ratio']], [wt['cv']],
                c='#EF4444', s=120, marker='*', zorder=5, label='Wild-type', edgecolors='black')
    ax2.set_xlabel('Kd ratio (weak/strong)', fontsize=STYLE['label_size'])
    ax2.set_ylabel('CV', fontsize=STYLE['label_size'])
    ax2.set_xscale('log')
    ax2.axvline(20.4, color='gray', linestyle=':', alpha=0.5)
    ax2.legend(fontsize=STYLE['legend_size'])
    despine(ax2)
    add_panel_label(ax2, 'B')

    fig.suptitle('Figure 20: Local sequence-phenotype tradeoff landscape',
                 fontsize=STYLE['title_size'], y=1.02)
    plt.tight_layout()
    path = os.path.join(OUT_DIR, 'fig20_landscape.png')
    fig.savefig(path, dpi=STYLE['dpi'], bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved: {path}")
    return path

if __name__ == '__main__':
    generate()
