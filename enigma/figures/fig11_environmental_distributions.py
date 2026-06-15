"""
phase7_extended_figures/fig11_environmental_distributions.py

Figure 11: 3x4 panel grid showing target protein distributions
under different regulatory architectures and environmental conditions.

Rows:    asymmetric, symmetric, single_site
Columns: baseline, cholesterol, acidic_pH, host_like
"""

import sys
sys.path.insert(0, '/Users/aayanalwani/tb project/mce3r_stochastic')

import os
import numpy as np

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from figures.figure_style import COLORS, STYLE, apply_style, despine, add_panel_label

# ── Configuration ─────────────────────────────────────────────────────────────

ARCHITECTURES = ['asymmetric', 'symmetric', 'single_site']
ENVIRONMENTS = ['baseline', 'cholesterol', 'acidic_pH', 'host_like']

ARCH_LABELS = {
    'asymmetric':  'Asymmetric (Mce3R)',
    'symmetric':   'Symmetric',
    'single_site': 'Single-site',
}

ENV_LABELS = {
    'baseline':    'Baseline',
    'cholesterol': 'Cholesterol',
    'acidic_pH':   'Acidic pH',
    'host_like':   'Host-like',
}

# Environment colors from config/parameters.py
ENV_COLORS = {
    'baseline':    '#6B7280',
    'cholesterol': '#F59E0B',
    'acidic_pH':   '#EF4444',
    'host_like':   '#7C3AED',
}

DATA_DIR = '/Users/aayanalwani/tb project/mce3r_stochastic/results/phase6'
OUT_DIR = '/Users/aayanalwani/tb project/mce3r_stochastic/results/figures'
OUT_FILE = os.path.join(OUT_DIR, 'fig11_environmental_distributions.png')

N_BINS = 50


# ── Main figure ───────────────────────────────────────────────────────────────

def make_figure():
    """Generate the 3x4 environmental-distribution panel figure."""
    apply_style()
    os.makedirs(OUT_DIR, exist_ok=True)

    fig, axes = plt.subplots(3, 4, figsize=(14, 8), sharex=True)

    panel_idx = 0
    panel_letters = 'ABCDEFGHIJKL'

    for row, arch in enumerate(ARCHITECTURES):
        for col, env in enumerate(ENVIRONMENTS):
            ax = axes[row, col]

            # Load data
            fname = f'env_condition_{arch}_{env}_single.npz'
            fpath = os.path.join(DATA_DIR, fname)
            data = np.load(fpath)
            proteins = data['proteins']

            # Compute statistics
            mean_val = np.mean(proteins)
            std_val = np.std(proteins)
            cv_val = std_val / mean_val if mean_val > 0 else 0.0

            # Histogram
            color = ENV_COLORS[env]
            counts, bin_edges = np.histogram(proteins, bins=N_BINS)
            # Normalize to density
            widths = np.diff(bin_edges)
            density = counts / (counts.sum() * widths)

            # Step histogram with alpha fill
            ax.step(bin_edges[:-1], density, where='post',
                    color=color, linewidth=STYLE['linewidth'])
            ax.fill_between(bin_edges[:-1], density, step='post',
                            color=color, alpha=0.25)

            # Annotations: mean and CV
            ax.text(0.95, 0.95,
                    f'mean = {mean_val:.0f}\nCV = {cv_val:.2f}',
                    transform=ax.transAxes, fontsize=STYLE['legend_size'],
                    va='top', ha='right',
                    bbox=dict(facecolor='white', edgecolor='none', alpha=0.8))

            # Panel label
            add_panel_label(ax, panel_letters[panel_idx])
            panel_idx += 1

            despine(ax)

            # Column header (top row only)
            if row == 0:
                ax.set_title(ENV_LABELS[env], fontsize=STYLE['title_size'],
                             fontweight='bold')

            # Row label (left column only)
            if col == 0:
                ax.set_ylabel(ARCH_LABELS[arch], fontsize=STYLE['label_size'],
                              fontweight='bold')

            # Bottom row x-label
            if row == len(ARCHITECTURES) - 1:
                ax.set_xlabel('Target protein count')

    fig.suptitle('Figure 11: Protein distributions across environments',
                 fontsize=STYLE['title_size'] + 1, fontweight='bold', y=1.02)
    fig.tight_layout()
    fig.savefig(OUT_FILE, dpi=STYLE['dpi'], bbox_inches='tight')
    plt.close(fig)
    print(f"Saved: {OUT_FILE}")
    return OUT_FILE


# ── Self-tests ────────────────────────────────────────────────────────────────

def run_tests():
    """Verify the generated figure meets publication requirements."""
    n_pass = 0
    n_fail = 0

    # Test 1: file exists
    if os.path.isfile(OUT_FILE):
        print("PASS: Output file exists")
        n_pass += 1
    else:
        print("FAIL: Output file not found")
        n_fail += 1
        return n_pass, n_fail

    # Test 2: file size > 10 KB
    size_kb = os.path.getsize(OUT_FILE) / 1024
    if size_kb > 10:
        print(f"PASS: File size = {size_kb:.1f} KB (> 10 KB)")
        n_pass += 1
    else:
        print(f"FAIL: File size = {size_kb:.1f} KB (expected > 10 KB)")
        n_fail += 1

    # Test 3: valid PNG header
    with open(OUT_FILE, 'rb') as f:
        header = f.read(8)
    png_sig = b'\x89PNG\r\n\x1a\n'
    if header == png_sig:
        print("PASS: Valid PNG header")
        n_pass += 1
    else:
        print("FAIL: Invalid PNG header")
        n_fail += 1

    return n_pass, n_fail


if __name__ == '__main__':
    make_figure()
    n_pass, n_fail = run_tests()
    print(f"\n{'=' * 50}")
    print(f"Self-tests: {n_pass} PASS, {n_fail} FAIL")
    if n_fail > 0:
        sys.exit(1)
    else:
        print("All checks passed.")
        sys.exit(0)


generate = make_figure
