"""
phase7_extended_figures/fig9_repression_curves.py

Figure 9: Repression fold vs [Mce3R] concentration for 4 operator
architectures.  Loads Phase 5 repression-curve data and produces a
publication-quality plot saved to results/figures/.
"""

import sys
sys.path.insert(0, '/Users/aayanalwani/tb project/mce3r_stochastic')

import os
import numpy as np
import pandas as pd

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from figures.figure_style import COLORS, STYLE, apply_style, despine, add_panel_label

# ── Paths ────────────────────────────────────────────────────────────────────
BASE_DIR = '/Users/aayanalwani/tb project/mce3r_stochastic'
CSV_PATH = os.path.join(BASE_DIR, 'results', 'phase5', 'repression_curves.csv')
OUT_DIR  = os.path.join(BASE_DIR, 'results', 'figures')
OUT_PATH = os.path.join(OUT_DIR, 'fig9_repression_curves.png')

# ── Architecture styling ─────────────────────────────────────────────────────
# The CSV uses 'strong_only' where the project nomenclature says 'asymmetric'.
# Map CSV labels to display labels and colours.
ARCH_CONFIG = {
    'native_asymmetric': {
        'label': 'Native asymmetric',
        'color': '#DC2626',
    },
    'strong_only': {
        'label': 'Asymmetric (strong only)',
        'color': COLORS['asymmetric'],       # '#0D9488'
    },
    'symmetric': {
        'label': 'Symmetric',
        'color': COLORS['symmetric'],        # '#F97316'
    },
    'single_site': {
        'label': 'Single site',
        'color': COLORS['single_site'],      # '#8B5CF6'
    },
}

# Preferred drawing order (back-to-front)
ARCH_ORDER = ['single_site', 'symmetric', 'strong_only', 'native_asymmetric']

# Physiological Mce3R concentration (nM)
PHYS_CONC_NM = 332.0


# ── Main figure function ─────────────────────────────────────────────────────
def make_figure():
    """Generate and save Figure 9."""
    apply_style()

    # Load data
    df = pd.read_csv(CSV_PATH)

    # Identify available architectures
    available = [a for a in ARCH_ORDER if a in df['architecture'].unique()]
    if not available:
        raise ValueError(
            f"No recognised architectures in CSV. "
            f"Found: {df['architecture'].unique().tolist()}"
        )

    # ── Create figure ────────────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(5.5, 4.0))

    # Determine transition region across all architectures: the concentration
    # range where fold_repression changes most rapidly.  We define it as the
    # central band that brackets the steepest gradient for each architecture.
    transition_lo, transition_hi = _find_transition_region(df, available)

    # Shade transition region
    if transition_lo is not None and transition_hi is not None:
        ax.axvspan(
            np.log10(transition_lo), np.log10(transition_hi),
            color='#E5E7EB', alpha=0.55, zorder=0, label='Transition region',
        )

    # Plot each architecture
    for arch in available:
        sub = df[df['architecture'] == arch].sort_values('concentration_nM')
        cfg = ARCH_CONFIG[arch]
        ax.plot(
            np.log10(sub['concentration_nM']),
            sub['fold_repression'],
            color=cfg['color'],
            linewidth=STYLE['linewidth'],
            label=cfg['label'],
            zorder=2,
        )

    # Physiological concentration line
    ax.axvline(
        np.log10(PHYS_CONC_NM),
        color='#374151', linestyle='--', linewidth=1.0,
        zorder=1, label=f'Physiological [Mce3R]\n({PHYS_CONC_NM:.0f} nM)',
    )

    # ── Axes labels and formatting ───────────────────────────────────────
    ax.set_xlabel(r'log$_{10}$([Mce3R] / nM)')
    ax.set_ylabel('Fold repression')
    ax.set_title('Repression curves across operator architectures')

    despine(ax)

    # Legend outside the data area
    ax.legend(
        loc='upper left',
        frameon=False,
        fontsize=STYLE['legend_size'],
    )

    # Panel label
    add_panel_label(ax, 'I')

    # ── Save ─────────────────────────────────────────────────────────────
    os.makedirs(OUT_DIR, exist_ok=True)
    fig.savefig(OUT_PATH, dpi=STYLE['dpi'], format='png')
    plt.close(fig)
    print(f"Saved: {OUT_PATH}")
    return OUT_PATH


def _find_transition_region(df, architectures):
    """Return (lo_nM, hi_nM) bracketing the steepest-gradient region.

    Strategy: for each architecture compute the numerical derivative of
    fold_repression w.r.t. log10(concentration).  The transition region
    is the concentration band where the derivative exceeds 25 % of its
    peak value, union-ed across architectures.
    """
    lo_vals, hi_vals = [], []
    for arch in architectures:
        sub = df[df['architecture'] == arch].sort_values('concentration_nM')
        if len(sub) < 5:
            continue
        conc = sub['concentration_nM'].values
        fr   = sub['fold_repression'].values
        log_conc = np.log10(conc)
        grad = np.gradient(fr, log_conc)
        peak = np.max(np.abs(grad))
        if peak < 1e-6:
            continue
        mask = np.abs(grad) > 0.25 * peak
        if mask.any():
            lo_vals.append(conc[mask].min())
            hi_vals.append(conc[mask].max())
    if lo_vals and hi_vals:
        return min(lo_vals), max(hi_vals)
    return None, None


# ── Self-tests ───────────────────────────────────────────────────────────────
if __name__ == '__main__':
    out = make_figure()

    n_pass = 0
    n_fail = 0

    # TEST 1: Output file exists
    if os.path.isfile(out):
        print("PASS: Output file exists")
        n_pass += 1
    else:
        print("FAIL: Output file does not exist")
        n_fail += 1

    # TEST 2: File size > 10 KB
    if os.path.isfile(out):
        size_kb = os.path.getsize(out) / 1024
        if size_kb > 10:
            print(f"PASS: File size = {size_kb:.1f} KB (> 10 KB)")
            n_pass += 1
        else:
            print(f"FAIL: File size = {size_kb:.1f} KB (<= 10 KB)")
            n_fail += 1
    else:
        print("FAIL: Cannot check size — file missing")
        n_fail += 1

    # TEST 3: Valid PNG (magic bytes)
    if os.path.isfile(out):
        with open(out, 'rb') as fh:
            header = fh.read(8)
        if header[:4] == b'\x89PNG':
            print("PASS: Valid PNG header")
            n_pass += 1
        else:
            print("FAIL: Invalid PNG header")
            n_fail += 1
    else:
        print("FAIL: Cannot check PNG header — file missing")
        n_fail += 1

    print(f"\n{'='*50}")
    print(f"Self-tests: {n_pass} PASS, {n_fail} FAIL")
    sys.exit(0 if n_fail == 0 else 1)


generate = make_figure
