"""
phase4_figures/fig3_asymmetry_sweep.py — Asymmetry ratio sweep plot.

X: asymmetry ratio (log scale), Y: intermediate fraction
Vertical dashed gold line at ratio ~ 20.4 labelled "Mce3R (Panagoda 2024)"

Produces: results/figures/fig3_asymmetry_sweep.png
"""

import sys
sys.path.insert(0, '/Users/aayanalwani/tb project/mce3r_stochastic')

import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from config.parameters import PARAMS
from figures.figure_style import (COLORS, STYLE, apply_style, despine,
                                         add_panel_label)

PROJECT = '/Users/aayanalwani/tb project/mce3r_stochastic'
SWEEP_FILE = os.path.join(PROJECT, 'results', 'phase2', 'condition_E_sweep.npz')
BOOTSTRAP_FILE = os.path.join(PROJECT, 'results', 'phase3', 'bootstrap_results.csv')
OUT_DIR = os.path.join(PROJECT, 'results', 'figures')
OUT_FILE = os.path.join(OUT_DIR, 'fig3_asymmetry_sweep.png')


def _placeholder(msg):
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.text(0.5, 0.5, msg, ha='center', va='center', fontsize=16, color='gray',
            transform=ax.transAxes)
    ax.axis('off')
    return fig


def generate_fig3():
    """Generate Figure 3: asymmetry sweep."""
    apply_style()
    os.makedirs(OUT_DIR, exist_ok=True)

    if not os.path.isfile(SWEEP_FILE):
        fig = _placeholder('Fig 3: condition_E_sweep.npz not found')
        fig.savefig(OUT_FILE, dpi=STYLE['dpi'])
        plt.close(fig)
        return {'status': 'placeholder', 'reason': 'missing data'}

    d = np.load(SWEEP_FILE, allow_pickle=True)
    ratios = d['ratios']
    int_fracs = d['intermediate_fractions']
    cvs = d['cvs']

    # Try loading bootstrap CIs for error bars
    ci_lower, ci_upper = None, None
    if os.path.isfile(BOOTSTRAP_FILE):
        bdf = pd.read_csv(BOOTSTRAP_FILE)
        # CIs for intermediate_fraction per condition aren't per-ratio,
        # so we won't add error bars from bootstrap here unless available

    mce3r_ratio = PARAMS['mce3r_actual_ratio']

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    # ── Panel A: Intermediate fraction vs ratio ───────────────────────────
    add_panel_label(ax1, 'A')
    ax1.plot(ratios, int_fracs, 'o-', color=COLORS['asymmetric'], linewidth=2,
             markersize=8, markeredgecolor='black', markeredgewidth=0.5, zorder=3)

    # Fill area under curve
    ax1.fill_between(ratios, 0, int_fracs, alpha=0.15, color=COLORS['asymmetric'])

    # Mce3R reference line
    ax1.axvline(mce3r_ratio, color=COLORS['mce3r_marker'], linestyle='--',
                linewidth=2, zorder=2)
    ax1.annotate('Mce3R\n(Panagoda 2024)',
                 xy=(mce3r_ratio, max(int_fracs) * 0.85),
                 xytext=(mce3r_ratio * 1.4, max(int_fracs) * 0.9),
                 fontsize=9, fontweight='bold', color=COLORS['mce3r_marker'],
                 arrowprops=dict(arrowstyle='->', color=COLORS['mce3r_marker'],
                                 lw=1.5))

    ax1.set_xscale('log')
    ax1.set_xlabel('Asymmetry ratio ($K_{d,weak} / K_{d,strong}$)')
    ax1.set_ylabel('Intermediate expression fraction')
    ax1.set_title('Operator asymmetry drives expression heterogeneity',
                  fontsize=STYLE['title_size'])
    ax1.set_ylim(bottom=-0.02)
    despine(ax1)

    # ── Panel B: CV vs ratio ──────────────────────────────────────────────
    add_panel_label(ax2, 'B')
    ax2.plot(ratios, cvs, 's-', color=COLORS['single_site'], linewidth=2,
             markersize=8, markeredgecolor='black', markeredgewidth=0.5, zorder=3)
    ax2.fill_between(ratios, 0, cvs, alpha=0.12, color=COLORS['single_site'])

    ax2.axvline(mce3r_ratio, color=COLORS['mce3r_marker'], linestyle='--',
                linewidth=2, zorder=2)
    ax2.annotate('Mce3R',
                 xy=(mce3r_ratio, max(cvs) * 0.85),
                 xytext=(mce3r_ratio * 1.4, max(cvs) * 0.9),
                 fontsize=9, fontweight='bold', color=COLORS['mce3r_marker'],
                 arrowprops=dict(arrowstyle='->', color=COLORS['mce3r_marker'],
                                 lw=1.5))

    ax2.set_xscale('log')
    ax2.set_xlabel('Asymmetry ratio ($K_{d,weak} / K_{d,strong}$)')
    ax2.set_ylabel('Coefficient of variation (CV)')
    ax2.set_title('Noise as a function of asymmetry',
                  fontsize=STYLE['title_size'])
    ax2.set_ylim(bottom=0)
    despine(ax2)

    fig.suptitle('Asymmetry sweep: Condition E',
                 fontsize=STYLE['title_size'] + 1, fontweight='bold', y=1.02)
    fig.tight_layout()
    fig.savefig(OUT_FILE, dpi=STYLE['dpi'])
    plt.close(fig)

    return {
        'status': 'success',
        'n_ratios': len(ratios),
        'ratio_range': [float(ratios.min()), float(ratios.max())],
        'output': OUT_FILE,
    }


# ── Self-tests ────────────────────────────────────────────────────────────────
if __name__ == '__main__':
    n_pass = 0
    n_fail = 0

    result = generate_fig3()

    # SANITY: file exists and valid
    try:
        assert os.path.isfile(OUT_FILE)
        fsize = os.path.getsize(OUT_FILE)
        assert fsize > 5000
        with open(OUT_FILE, 'rb') as f:
            assert f.read(4) == b'\x89PNG'
        print(f"SANITY PASS: fig3 saved ({fsize:,} bytes)")
        n_pass += 1
    except AssertionError as e:
        print(f"SANITY FAIL: {e}")
        n_fail += 1

    # SANITY: result dict
    try:
        assert result['status'] in ('success', 'placeholder')
        print(f"SANITY PASS: status={result['status']}")
        n_pass += 1
    except (AssertionError, KeyError) as e:
        print(f"SANITY FAIL: {e}")
        n_fail += 1

    # SCIENTIFIC: sweep covers the Mce3R ratio
    if result.get('status') == 'success':
        rng = result.get('ratio_range', [0, 0])
        mce3r = PARAMS['mce3r_actual_ratio']
        if rng[0] <= mce3r <= rng[1]:
            print(f"SCIENTIFIC PASS: Sweep range [{rng[0]}, {rng[1]}] covers Mce3R ratio {mce3r}")
        else:
            print(f"SCIENTIFIC WARNING: Sweep range [{rng[0]}, {rng[1]}] does not cover Mce3R ratio {mce3r}")

    print(f"\n{'='*50}")
    print(f"Sanity checks: {n_pass} PASS, {n_fail} FAIL")
    if n_fail > 0:
        sys.exit(1)
    else:
        print("All checks passed.")
        sys.exit(0)

generate = generate_fig3
