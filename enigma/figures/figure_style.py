"""
phase4_figures/figure_style.py — Shared style config for publication figures.

Provides STYLE dict, COLORS dict, and apply_style() function.
"""

import sys
sys.path.insert(0, '/Users/aayanalwani/tb project/mce3r_stochastic')

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib as mpl

# ── Color scheme ──────────────────────────────────────────────────────────────
COLORS = {
    'asymmetric':   '#0D9488',   # Teal    — Condition A
    'symmetric':    '#F97316',   # Orange  — Condition B
    'single_site':  '#8B5CF6',   # Purple  — Condition C
    'unregulated':  '#6B7280',   # Gray    — Condition D
    'sweep_cmap':   'viridis',
    'mce3r_marker': '#E9C46A',   # Gold
}

# Condition name → color mapping (handy for loops)
CONDITION_COLORS = {
    'A': COLORS['asymmetric'],
    'B': COLORS['symmetric'],
    'C': COLORS['single_site'],
    'D': COLORS['unregulated'],
}

CONDITION_LABELS = {
    'A': 'Asymmetric (Mce3R)',
    'B': 'Symmetric control',
    'C': 'Single-site control',
    'D': 'Unregulated control',
}

# ── Style parameters ──────────────────────────────────────────────────────────
STYLE = {
    'font_family':   'Arial',
    'font_fallback': 'Helvetica',
    'label_size':    12,
    'title_size':    13,
    'tick_size':     10,
    'legend_size':   9,
    'dpi':           300,
    'fig_format':    'png',
    'linewidth':     1.5,
    'spine_width':   0.8,
}


def apply_style():
    """Apply publication-quality matplotlib style globally."""
    # Try Arial, fall back to Helvetica, then sans-serif default
    try:
        mpl.font_manager.findfont('Arial', fallback_to_default=False)
        font_family = 'Arial'
    except Exception:
        try:
            mpl.font_manager.findfont('Helvetica', fallback_to_default=False)
            font_family = 'Helvetica'
        except Exception:
            font_family = 'sans-serif'

    plt.rcParams.update({
        'font.family':        'sans-serif',
        'font.sans-serif':    [font_family, 'Arial', 'Helvetica', 'DejaVu Sans'],
        'font.size':          STYLE['tick_size'],
        'axes.labelsize':     STYLE['label_size'],
        'axes.titlesize':     STYLE['title_size'],
        'xtick.labelsize':    STYLE['tick_size'],
        'ytick.labelsize':    STYLE['tick_size'],
        'legend.fontsize':    STYLE['legend_size'],
        'axes.spines.top':    False,
        'axes.spines.right':  False,
        'axes.linewidth':     STYLE['spine_width'],
        'figure.dpi':         STYLE['dpi'],
        'savefig.dpi':        STYLE['dpi'],
        'savefig.bbox':       'tight',
        'savefig.pad_inches': 0.1,
    })


def despine(ax):
    """Remove top and right spines from an axes object."""
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)


def add_panel_label(ax, label, x=-0.12, y=1.08, fontsize=16):
    """Add a bold panel label (A, B, C, ...) to an axes."""
    ax.text(x, y, label, transform=ax.transAxes,
            fontsize=fontsize, fontweight='bold', va='top', ha='left')


# ── Self-tests ────────────────────────────────────────────────────────────────
if __name__ == '__main__':
    n_pass = 0
    n_fail = 0

    # SANITY: COLORS has required keys
    required_color_keys = ['asymmetric', 'symmetric', 'single_site',
                           'unregulated', 'sweep_cmap', 'mce3r_marker']
    try:
        for k in required_color_keys:
            assert k in COLORS, f"Missing key: {k}"
        print(f"SANITY PASS: COLORS has all {len(required_color_keys)} required keys")
        n_pass += 1
    except AssertionError as e:
        print(f"SANITY FAIL: {e}")
        n_fail += 1

    # SANITY: STYLE has required keys
    required_style_keys = ['font_family', 'label_size', 'title_size', 'dpi']
    try:
        for k in required_style_keys:
            assert k in STYLE, f"Missing key: {k}"
        print(f"SANITY PASS: STYLE has required keys")
        n_pass += 1
    except AssertionError as e:
        print(f"SANITY FAIL: {e}")
        n_fail += 1

    # SANITY: apply_style runs without error
    try:
        apply_style()
        print("SANITY PASS: apply_style() executed successfully")
        n_pass += 1
    except Exception as e:
        print(f"SANITY FAIL: apply_style() raised {e}")
        n_fail += 1

    # SANITY: CONDITION_COLORS maps A-D
    try:
        for c in ['A', 'B', 'C', 'D']:
            assert c in CONDITION_COLORS
        print("SANITY PASS: CONDITION_COLORS maps A, B, C, D")
        n_pass += 1
    except AssertionError as e:
        print(f"SANITY FAIL: {e}")
        n_fail += 1

    # SCIENTIFIC: DPI is publication quality
    if STYLE['dpi'] >= 300:
        print(f"SCIENTIFIC PASS: DPI = {STYLE['dpi']} (>= 300)")
    else:
        print(f"SCIENTIFIC WARNING: DPI = {STYLE['dpi']} is below 300 for publication")

    print(f"\n{'='*50}")
    print(f"Sanity checks: {n_pass} PASS, {n_fail} FAIL")
    if n_fail > 0:
        sys.exit(1)
    else:
        print("All checks passed.")
        sys.exit(0)
