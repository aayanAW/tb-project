"""
phase7_extended_figures/fig12_mutual_information.py

Figure 12: Mutual Information (bits) vs [Mce3R] concentration
for three promoter architectures (asymmetric, symmetric, single-site).
"""

import sys
sys.path.insert(0, '/Users/aayanalwani/tb project/mce3r_stochastic')

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import os

from figures.figure_style import COLORS, STYLE, apply_style, despine, add_panel_label

# ── Paths ────────────────────────────────────────────────────────────────────
BASE = '/Users/aayanalwani/tb project/mce3r_stochastic'
DATA_PATH = os.path.join(BASE, 'results', 'phase6', 'mutual_information.csv')
OUT_DIR = os.path.join(BASE, 'results', 'figures')
OUT_PATH = os.path.join(OUT_DIR, 'fig12_mutual_information.png')

# ── Physiological concentration ──────────────────────────────────────────────
PHYS_CONC = 332  # nM

# ── Architecture styling ─────────────────────────────────────────────────────
ARCH_STYLE = {
    'asymmetric':  {'color': COLORS['asymmetric'],  'label': 'Asymmetric'},
    'symmetric':   {'color': COLORS['symmetric'],   'label': 'Symmetric'},
    'single_site': {'color': COLORS['single_site'], 'label': 'Single-site'},
}


def main():
    apply_style()

    # Load data
    df = pd.read_csv(DATA_PATH)

    # Create figure
    fig, ax = plt.subplots(figsize=(6, 4))

    # Plot each architecture
    for arch, style in ARCH_STYLE.items():
        sub = df[df['architecture'] == arch].sort_values('concentration')
        ax.plot(sub['concentration'], sub['MI_bits'],
                color=style['color'], label=style['label'],
                linewidth=STYLE['linewidth'], marker='o', markersize=4)

    # Vertical line at physiological concentration
    ax.axvline(PHYS_CONC, color='#6B7280', linestyle='--', linewidth=1.0,
               alpha=0.7, zorder=0)
    y_min, y_max = ax.get_ylim()
    ax.text(PHYS_CONC * 1.1, y_max * 0.95,
            f'Physiological\n({PHYS_CONC} nM)',
            fontsize=STYLE['legend_size'], color='#6B7280',
            va='top', ha='left')

    # Axes
    ax.set_xscale('log')
    ax.set_xlabel('[Mce3R] (nM)')
    ax.set_ylabel('Mutual Information (bits)')

    # Legend and spines
    ax.legend(frameon=False)
    despine(ax)

    # Save
    os.makedirs(OUT_DIR, exist_ok=True)
    fig.savefig(OUT_PATH, dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f"Saved: {OUT_PATH}")


# ── Self-tests ───────────────────────────────────────────────────────────────
def run_tests():
    n_pass = 0
    n_fail = 0

    # TEST 1: output file exists
    if os.path.isfile(OUT_PATH):
        print("PASS: Output file exists")
        n_pass += 1
    else:
        print("FAIL: Output file does not exist")
        n_fail += 1
        return n_pass, n_fail

    # TEST 2: file size > 10 KB
    size_kb = os.path.getsize(OUT_PATH) / 1024
    if size_kb > 10:
        print(f"PASS: File size = {size_kb:.1f} KB (> 10 KB)")
        n_pass += 1
    else:
        print(f"FAIL: File size = {size_kb:.1f} KB (<= 10 KB)")
        n_fail += 1

    # TEST 3: valid PNG header
    with open(OUT_PATH, 'rb') as f:
        header = f.read(8)
    if header[:4] == b'\x89PNG':
        print("PASS: Valid PNG header")
        n_pass += 1
    else:
        print("FAIL: Invalid PNG header")
        n_fail += 1

    return n_pass, n_fail


if __name__ == '__main__':
    main()
    n_pass, n_fail = run_tests()
    print(f"\n{'='*50}")
    print(f"Self-tests: {n_pass} PASS, {n_fail} FAIL")
    if n_fail > 0:
        sys.exit(1)
    else:
        print("All checks passed.")
        sys.exit(0)


generate = main
