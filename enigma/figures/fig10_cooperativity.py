"""
phase7_extended_figures/fig10_cooperativity.py — Figure 10: Cooperativity analysis.

Panel A: 2D histogram of MCMC posteriors for log10(omega) and dG_spacer.
Panel B: Heatmap of fold repression vs (omega, dG_spacer) at 332 nM Mce3R.
"""

import sys
sys.path.insert(0, '/Users/aayanalwani/tb project/mce3r_stochastic')

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import json
import os

from figures.figure_style import COLORS, STYLE, apply_style, despine, add_panel_label
from phase5_thermodynamic.partition_function import (
    repression_fold, mu_from_concentration,
)
from config.parameters import kT_KCAL

# ── Paths ────────────────────────────────────────────────────────────────────
BASE = '/Users/aayanalwani/tb project/mce3r_stochastic'
POSTERIOR_PATH = os.path.join(BASE, 'results/phase5/mcmc_posteriors.npz')
SUMMARY_PATH  = os.path.join(BASE, 'results/phase5/mcmc_summary.json')
ENERGY_PATH   = os.path.join(BASE, 'results/phase5/energy_parameters.json')
OUT_DIR       = os.path.join(BASE, 'results/figures')
OUT_PATH      = os.path.join(OUT_DIR, 'fig10_cooperativity.png')


def main():
    apply_style()
    os.makedirs(OUT_DIR, exist_ok=True)

    # ── Load data ────────────────────────────────────────────────────────────
    posteriors = np.load(POSTERIOR_PATH)
    chain = posteriors['chain']  # (n_walkers, n_steps, 2): col0=log10(omega), col1=dG_spacer

    with open(SUMMARY_PATH) as f:
        summary = json.load(f)

    with open(ENERGY_PATH) as f:
        energy = json.load(f)

    # Flatten chain (discard burn-in already handled at save time per summary)
    n_walkers, n_steps, _ = chain.shape
    n_burn = summary.get('n_burn', 1000)
    flat = chain[:, n_burn:, :].reshape(-1, 2)
    log10_omega_samples = flat[:, 0]
    dG_spacer_samples   = flat[:, 1]

    # Median values
    omega_median     = summary['omega_median']
    log10_omega_med  = np.log10(omega_median)
    dG_spacer_med    = summary['dG_spacer_median']

    # Calibrated binding energies
    dG_strong = energy['dG_strong']
    dG_weak   = energy['dG_weak']
    conc_nM   = summary.get('mce3r_conc_nM', 332.0)

    # ── Figure ───────────────────────────────────────────────────────────────
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.5))

    # ── Panel A: 2D histogram of posteriors ──────────────────────────────────
    ax_a = axes[0]
    h = ax_a.hist2d(
        log10_omega_samples,
        dG_spacer_samples,
        bins=60,
        cmap='viridis',
        cmin=1,
    )
    plt.colorbar(h[3], ax=ax_a, label='Samples', shrink=0.85)

    # Mark the median
    ax_a.plot(log10_omega_med, dG_spacer_med, marker='*', ms=14,
              color=COLORS['mce3r_marker'], markeredgecolor='white',
              markeredgewidth=0.8, zorder=5)

    ax_a.set_xlabel(r'log$_{10}$($\omega$)')
    ax_a.set_ylabel(r'$\Delta G_{\mathrm{spacer}}$ (kcal/mol)')
    ax_a.set_title('MCMC posterior')
    despine(ax_a)
    add_panel_label(ax_a, 'A')

    # ── Panel B: Heatmap of fold repression ──────────────────────────────────
    ax_b = axes[1]

    # Build grid
    log10_omega_grid = np.linspace(-2, 3, 80)
    dG_spacer_grid   = np.linspace(-3, 3, 80)
    LOG10_OMEGA, DG_SPACER = np.meshgrid(log10_omega_grid, dG_spacer_grid)

    mu = mu_from_concentration(conc_nM)

    FOLD = np.zeros_like(LOG10_OMEGA)
    for i in range(LOG10_OMEGA.shape[0]):
        for j in range(LOG10_OMEGA.shape[1]):
            omega_val = 10.0 ** LOG10_OMEGA[i, j]
            dGs_val   = DG_SPACER[i, j]
            FOLD[i, j] = repression_fold(
                dG_strong, dG_weak, mu,
                omega=omega_val, dG_spacer=dGs_val, kT=kT_KCAL,
            )

    # Use log scale for fold repression
    FOLD_LOG = np.log10(np.clip(FOLD, 1.0, None))

    im = ax_b.pcolormesh(
        log10_omega_grid, dG_spacer_grid, FOLD_LOG,
        cmap='inferno', shading='auto',
    )
    cbar = plt.colorbar(im, ax=ax_b, label='log$_{10}$(fold repression)', shrink=0.85)

    # Mark MCMC median
    ax_b.plot(log10_omega_med, dG_spacer_med, marker='*', ms=14,
              color=COLORS['mce3r_marker'], markeredgecolor='white',
              markeredgewidth=0.8, zorder=5)

    ax_b.set_xlabel(r'log$_{10}$($\omega$)')
    ax_b.set_ylabel(r'$\Delta G_{\mathrm{spacer}}$ (kcal/mol)')
    ax_b.set_title(f'Fold repression at {conc_nM:.0f} nM Mce3R')
    add_panel_label(ax_b, 'B')

    fig.tight_layout()
    fig.savefig(OUT_PATH, dpi=STYLE['dpi'])
    plt.close(fig)
    print(f"Saved: {OUT_PATH}")

    # ── Self-tests ───────────────────────────────────────────────────────────
    n_pass = 0
    n_fail = 0

    # Test 1: file exists
    if os.path.isfile(OUT_PATH):
        print("PASS: output file exists")
        n_pass += 1
    else:
        print("FAIL: output file does not exist")
        n_fail += 1

    # Test 2: file size > 10 KB
    if os.path.isfile(OUT_PATH):
        size_kb = os.path.getsize(OUT_PATH) / 1024
        if size_kb > 10:
            print(f"PASS: file size = {size_kb:.1f} KB (> 10 KB)")
            n_pass += 1
        else:
            print(f"FAIL: file size = {size_kb:.1f} KB (should be > 10 KB)")
            n_fail += 1

    # Test 3: valid PNG header
    if os.path.isfile(OUT_PATH):
        with open(OUT_PATH, 'rb') as f:
            header = f.read(8)
        if header[:4] == b'\x89PNG':
            print("PASS: valid PNG header")
            n_pass += 1
        else:
            print("FAIL: invalid PNG header")
            n_fail += 1

    print(f"\n{'='*50}")
    print(f"Self-tests: {n_pass} PASS, {n_fail} FAIL")
    if n_fail > 0:
        sys.exit(1)
    else:
        print("All checks passed.")


if __name__ == '__main__':
    main()


generate = main
