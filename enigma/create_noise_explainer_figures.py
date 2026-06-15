#!/usr/bin/env python3
"""
Generate 6 explanatory figures about gene expression noise.
Designed for a high-school audience with clean, modern matplotlib styling.
"""

import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
from pathlib import Path

# ── Global style ──────────────────────────────────────────────────────────────
plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.size': 13,
    'axes.titlesize': 15,
    'axes.titleweight': 'bold',
    'axes.labelsize': 13,
    'axes.spines.top': False,
    'axes.spines.right': False,
    'figure.facecolor': 'white',
    'axes.facecolor': '#FAFAFA',
    'axes.edgecolor': '#CCCCCC',
})

OUTDIR = Path("/Users/aayanalwani/tb project/mce3r_stochastic/results/noise_explainers")
OUTDIR.mkdir(parents=True, exist_ok=True)

TEAL   = '#0D9488'
CORAL  = '#E8795A'
BLUE   = '#3B82F6'
DKBLUE = '#1E40AF'


# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 1 — What Is Noise?
# ══════════════════════════════════════════════════════════════════════════════
def fig1_what_is_noise():
    np.random.seed(42)
    n_cells = 20

    # Low-noise data
    low = np.random.normal(200, 15, n_cells)
    low = np.clip(low, 150, 250)

    # High-noise data — bimodal mixture
    high_main = np.random.normal(220, 60, int(n_cells * 0.7))
    high_low  = np.random.normal(80, 30, n_cells - len(high_main))
    high = np.concatenate([high_main, high_low])
    np.random.shuffle(high)
    high = np.clip(high, 20, 420)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 6), constrained_layout=True)

    cells = np.arange(1, n_cells + 1)

    # Panel A — low noise
    ax1.barh(cells, low, height=0.7, color=BLUE, alpha=0.85, edgecolor='white', linewidth=0.5)
    mean_low = np.mean(low)
    ax1.axvline(mean_low, color=DKBLUE, ls='--', lw=2, label=f'Mean = {mean_low:.0f}')
    ax1.set_title("Low Noise: All cells make\nsimilar amounts", fontsize=14, pad=10)
    ax1.set_xlabel("Protein count per cell", fontsize=13)
    ax1.set_ylabel("Cell #", fontsize=13)
    ax1.set_xlim(0, 450)
    ax1.legend(fontsize=12, loc='lower right')
    ax1.text(0.02, 1.02, "A", transform=ax1.transAxes, fontsize=18,
             fontweight='bold', va='bottom')

    # Panel B — high noise
    ax2.barh(cells, high, height=0.7, color=TEAL, alpha=0.85, edgecolor='white', linewidth=0.5)
    mean_high = np.mean(high)
    ax2.axvline(mean_high, color='#065F56', ls='--', lw=2, label=f'Mean = {mean_high:.0f}')
    ax2.set_title("High Noise: Cells make\nvery different amounts", fontsize=14, pad=10)
    ax2.set_xlabel("Protein count per cell", fontsize=13)
    ax2.set_xlim(0, 450)
    ax2.legend(fontsize=12, loc='lower right')
    ax2.text(0.02, 1.02, "B", transform=ax2.transAxes, fontsize=18,
             fontweight='bold', va='bottom')

    fig.suptitle("")
    fig.text(0.5, -0.02,
             "Same gene, same average  \u2192  but VERY different cell-to-cell variation",
             ha='center', fontsize=14, fontweight='bold', color='#333333',
             bbox=dict(boxstyle='round,pad=0.5', facecolor='#FEF3C7', edgecolor='#F59E0B', lw=1.5))

    path = OUTDIR / "what_is_noise.png"
    fig.savefig(path, dpi=300, bbox_inches='tight')
    plt.close(fig)
    return path


# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 2 — Why Noise Matters: Persistence
# ══════════════════════════════════════════════════════════════════════════════
def fig2_persistence():
    np.random.seed(42)
    n = 50

    # Random cell positions in a circle-ish blob
    theta = np.random.uniform(0, 2 * np.pi, n)
    r = np.sqrt(np.random.uniform(0, 1, n)) * 3
    x = r * np.cos(theta)
    y = r * np.sin(theta)

    # Expression levels: most high, a few low
    expression = np.random.beta(5, 1.5, n)  # skewed toward high
    expression[:6] = np.random.uniform(0.0, 0.25, 6)  # 6 persisters
    np.random.shuffle(expression)

    is_persister = expression < 0.3

    fig, axes = plt.subplots(1, 3, figsize=(10, 5), constrained_layout=True)

    # --- Panel A: Before antibiotics ---
    ax = axes[0]
    colors_a = plt.cm.RdYlGn(expression)
    ax.scatter(x, y, c=expression, cmap='RdYlGn', s=120, edgecolors='#555555',
               linewidth=0.8, vmin=0, vmax=1, zorder=2)
    ax.set_title("Before Antibiotics\n(Healthy Population)", fontsize=13, pad=8)
    ax.set_xlim(-4.5, 4.5); ax.set_ylim(-4.5, 4.5)
    ax.set_aspect('equal')
    ax.axis('off')
    ax.text(0.0, 1.05, "A", transform=ax.transAxes, fontsize=18, fontweight='bold')
    # legend
    ax.text(0.5, -0.08, "Green = active metabolism\nRed = dormant",
            transform=ax.transAxes, ha='center', fontsize=11, color='#555')

    # --- Panel B: Antibiotic treatment ---
    ax = axes[1]
    # Draw surviving (persister) cells normally
    ax.scatter(x[is_persister], y[is_persister], c=expression[is_persister],
               cmap='RdYlGn', s=120, edgecolors='#555555', linewidth=0.8,
               vmin=0, vmax=1, zorder=2)
    # Draw dead cells as gray with X
    ax.scatter(x[~is_persister], y[~is_persister], s=120, color='#D1D5DB',
               edgecolors='#9CA3AF', linewidth=0.8, zorder=1, alpha=0.5)
    ax.scatter(x[~is_persister], y[~is_persister], s=80, marker='x',
               color='#6B7280', linewidth=2, zorder=3)
    ax.set_title("Antibiotic Treatment\n(Drug Kills Active Cells)", fontsize=13, pad=8)
    ax.set_xlim(-4.5, 4.5); ax.set_ylim(-4.5, 4.5)
    ax.set_aspect('equal')
    ax.axis('off')
    ax.text(0.0, 1.05, "B", transform=ax.transAxes, fontsize=18, fontweight='bold')
    ax.text(0.5, -0.08, "X = killed by drug\nRed survivors = persisters",
            transform=ax.transAxes, ha='center', fontsize=11, color='#555')

    # --- Panel C: After treatment — persisters regrow ---
    ax = axes[2]
    np.random.seed(99)
    n2 = 50
    theta2 = np.random.uniform(0, 2 * np.pi, n2)
    r2 = np.sqrt(np.random.uniform(0, 1, n2)) * 3
    x2 = r2 * np.cos(theta2)
    y2 = r2 * np.sin(theta2)
    expression2 = np.random.beta(5, 1.5, n2)
    expression2[:7] = np.random.uniform(0.0, 0.25, 7)
    np.random.shuffle(expression2)
    ax.scatter(x2, y2, c=expression2, cmap='RdYlGn', s=120, edgecolors='#555555',
               linewidth=0.8, vmin=0, vmax=1, zorder=2)
    ax.set_title("After Treatment\n(Persisters Regrow)", fontsize=13, pad=8)
    ax.set_xlim(-4.5, 4.5); ax.set_ylim(-4.5, 4.5)
    ax.set_aspect('equal')
    ax.axis('off')
    ax.text(0.0, 1.05, "C", transform=ax.transAxes, fontsize=18, fontweight='bold')
    ax.text(0.5, -0.08, "Population rebounds\nwith new persisters",
            transform=ax.transAxes, ha='center', fontsize=11, color='#555')

    path = OUTDIR / "why_noise_matters_persistence.png"
    fig.savefig(path, dpi=300, bbox_inches='tight')
    plt.close(fig)
    return path


# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 3 — CV Explained
# ══════════════════════════════════════════════════════════════════════════════
def fig3_cv_explained():
    np.random.seed(42)
    mean_val = 200
    low_std, high_std = 15, 40
    low_data  = np.random.normal(mean_val, low_std,  10000)
    high_data = np.random.normal(mean_val, high_std, 10000)

    fig, ax = plt.subplots(figsize=(8, 5), constrained_layout=True)

    bins = np.linspace(80, 320, 60)
    ax.hist(low_data, bins=bins, density=True, alpha=0.65, color=BLUE,
            edgecolor='white', linewidth=0.5,
            label=f"CV = {low_std/mean_val:.3f}  (Low noise)")
    ax.hist(high_data, bins=bins, density=True, alpha=0.55, color=TEAL,
            edgecolor='white', linewidth=0.5,
            label=f"CV = {high_std/mean_val:.3f}  (High noise)")
    ax.axvline(mean_val, color='#333', ls='--', lw=2, label=f"Mean = {mean_val}")

    ax.set_xlabel("Protein count", fontsize=13)
    ax.set_ylabel("Probability density", fontsize=13)
    ax.set_title("Coefficient of Variation (CV) Measures Noise", fontsize=15,
                 fontweight='bold')
    ax.legend(fontsize=12, loc='upper right')

    # Arrows showing width
    peak_y = 0.025
    ax.annotate('', xy=(mean_val - low_std, peak_y + 0.003),
                xytext=(mean_val + low_std, peak_y + 0.003),
                arrowprops=dict(arrowstyle='<->', color=BLUE, lw=2))
    ax.text(mean_val, peak_y + 0.005, 'narrow', ha='center', fontsize=11,
            color=DKBLUE, fontweight='bold')

    ax.annotate('', xy=(mean_val - high_std, peak_y - 0.004),
                xytext=(mean_val + high_std, peak_y - 0.004),
                arrowprops=dict(arrowstyle='<->', color=TEAL, lw=2))
    ax.text(mean_val, peak_y - 0.006, 'wide', ha='center', fontsize=11,
            color='#065F56', fontweight='bold')

    # Text box
    textstr = ("CV = Std Dev / Mean\n"
               "Same average protein, but the\n"
               "teal distribution is much wider\n"
               "\u2192 more noise!")
    props = dict(boxstyle='round,pad=0.6', facecolor='#F0FDFA', edgecolor=TEAL, lw=1.5)
    ax.text(0.02, 0.95, textstr, transform=ax.transAxes, fontsize=12,
            va='top', bbox=props)

    path = OUTDIR / "cv_explained.png"
    fig.savefig(path, dpi=300, bbox_inches='tight')
    plt.close(fig)
    return path


# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 4 — Operator Noise Mechanism
# ══════════════════════════════════════════════════════════════════════════════
def fig4_operator_mechanism():
    fig, axes = plt.subplots(1, 3, figsize=(12, 4), constrained_layout=True)

    def draw_dna(ax, xlim=(-0.5, 5.5)):
        ax.set_xlim(*xlim)
        ax.set_ylim(-2.5, 3.5)
        ax.axis('off')
        # DNA backbone
        ax.plot([xlim[0] + 0.3, xlim[1] - 0.3], [2.2, 2.2], lw=4,
                color='#A3A3A3', solid_capstyle='round', zorder=1)

    def draw_box(ax, x, w, color, label, y=2.2):
        rect = FancyBboxPatch((x - w/2, y - 0.35), w, 0.7,
                              boxstyle="round,pad=0.08", facecolor=color,
                              edgecolor='#333', lw=1.5, zorder=2)
        ax.add_patch(rect)
        ax.text(x, y, label, ha='center', va='center', fontsize=9,
                fontweight='bold', color='white', zorder=3)

    # --- Panel A: Asymmetric ---
    ax = axes[0]
    draw_dna(ax)
    draw_box(ax, 1.5, 1.8, '#FB923C', 'Weak\nKd=49nM')
    draw_box(ax, 3.8, 1.8, '#991B1B', 'Strong\nKd=2.4nM')
    ax.set_title("Asymmetric (Native)", fontsize=14, fontweight='bold', pad=8)
    ax.text(0.5, 0.15, "4 states \u2192 many transitions\n\u2192 HIGH noise",
            transform=ax.transAxes, ha='center', fontsize=12,
            fontweight='bold', color='#991B1B',
            bbox=dict(boxstyle='round,pad=0.4', facecolor='#FEF2F2',
                      edgecolor='#991B1B', lw=1.2))
    # state labels
    states = ["empty", "strong\nonly", "weak\nonly", "both"]
    for i, s in enumerate(states):
        cx = 0.5 + i * 1.4
        ax.text(cx, 0.5, s, ha='center', va='center', fontsize=8,
                bbox=dict(boxstyle='round,pad=0.25', facecolor='#F5F5F5',
                          edgecolor='#999', lw=0.8))
    # arrows between states
    for i in range(3):
        ax.annotate('', xy=(0.5 + (i+1)*1.4 - 0.45, 0.5),
                    xytext=(0.5 + i*1.4 + 0.45, 0.5),
                    arrowprops=dict(arrowstyle='->', color='#666', lw=1.2))
    ax.text(-0.05, 1.05, "A", transform=ax.transAxes, fontsize=18, fontweight='bold')

    # --- Panel B: Symmetric ---
    ax = axes[1]
    draw_dna(ax)
    draw_box(ax, 1.5, 1.8, '#B91C1C', 'Site 1\nKd=10.8nM')
    draw_box(ax, 3.8, 1.8, '#B91C1C', 'Site 2\nKd=10.8nM')
    ax.set_title("Symmetric", fontsize=14, fontweight='bold', pad=8)
    ax.text(0.5, 0.15, "4 states but balanced\n\u2192 MEDIUM noise",
            transform=ax.transAxes, ha='center', fontsize=12,
            fontweight='bold', color='#B91C1C',
            bbox=dict(boxstyle='round,pad=0.4', facecolor='#FEF2F2',
                      edgecolor='#B91C1C', lw=1.2))
    ax.text(-0.05, 1.05, "B", transform=ax.transAxes, fontsize=18, fontweight='bold')

    # --- Panel C: Single site ---
    ax = axes[2]
    draw_dna(ax)
    draw_box(ax, 1.5, 1.8, '#991B1B', 'Strong\nKd=2.4nM')
    draw_box(ax, 3.8, 1.8, '#D4D4D4', 'Deleted', y=2.2)
    # override text color for deleted
    ax.texts[-1] if False else None  # placeholder
    # redraw label for deleted with dark text
    rect2 = FancyBboxPatch((3.8 - 0.9, 2.2 - 0.35), 1.8, 0.7,
                           boxstyle="round,pad=0.08", facecolor='#E5E5E5',
                           edgecolor='#999', lw=1.5, ls='--', zorder=2)
    ax.add_patch(rect2)
    ax.text(3.8, 2.2, 'Deleted', ha='center', va='center', fontsize=9,
            color='#888', zorder=3)
    ax.set_title("Single Site", fontsize=14, fontweight='bold', pad=8)
    ax.text(0.5, 0.15, "2 states \u2192 simple on/off\n\u2192 LOW noise",
            transform=ax.transAxes, ha='center', fontsize=12,
            fontweight='bold', color='#166534',
            bbox=dict(boxstyle='round,pad=0.4', facecolor='#F0FDF4',
                      edgecolor='#166534', lw=1.2))
    ax.text(-0.05, 1.05, "C", transform=ax.transAxes, fontsize=18, fontweight='bold')

    path = OUTDIR / "operator_noise_mechanism.png"
    fig.savefig(path, dpi=300, bbox_inches='tight')
    plt.close(fig)
    return path


# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 5 — Noise-to-Persistence Pipeline
# ══════════════════════════════════════════════════════════════════════════════
def fig5_pipeline():
    fig, ax = plt.subplots(figsize=(14, 3.5), constrained_layout=True)
    ax.set_xlim(-0.5, 14)
    ax.set_ylim(-1.5, 3.5)
    ax.axis('off')

    boxes = [
        "Asymmetric\nOperator\n(Kd ratio = 20.4\u00d7)",
        "More Operator\nState\nTransitions",
        "Higher Gene\nExpression Noise\n(CV = 0.187)",
        "Bimodal\nProtein\nDistribution",
        "Persister\nSub-\npopulation",
        "Antibiotic\nSurvival",
    ]
    arrow_labels = ["causes", "produces", "creates", "enables", "leads to"]

    n = len(boxes)
    box_w, box_h = 1.9, 2.0
    gap = 0.45
    total_w = n * box_w + (n - 1) * gap
    x_start = (14 - total_w) / 2

    for i, label in enumerate(boxes):
        cx = x_start + i * (box_w + gap) + box_w / 2
        cy = 1.0
        rect = FancyBboxPatch((cx - box_w/2, cy - box_h/2), box_w, box_h,
                              boxstyle="round,pad=0.15", facecolor=TEAL,
                              edgecolor='#065F56', lw=2, zorder=2)
        ax.add_patch(rect)
        ax.text(cx, cy, label, ha='center', va='center', fontsize=11,
                fontweight='bold', color='white', zorder=3, linespacing=1.2)

        # Arrow to next box
        if i < n - 1:
            ax_start = cx + box_w / 2 + 0.03
            ax_end = cx + box_w / 2 + gap - 0.03
            ax.annotate('', xy=(ax_end, cy), xytext=(ax_start, cy),
                        arrowprops=dict(arrowstyle='->', lw=3, color='#334155',
                                        mutation_scale=20))
            ax.text((ax_start + ax_end) / 2, cy + box_h/2 + 0.15,
                    arrow_labels[i], ha='center', fontsize=10, fontstyle='italic',
                    color='#64748B')

    ax.set_title("From Operator Architecture to Antibiotic Survival",
                 fontsize=16, fontweight='bold', pad=12)

    path = OUTDIR / "noise_to_persistence_pipeline.png"
    fig.savefig(path, dpi=300, bbox_inches='tight')
    plt.close(fig)
    return path


# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 6 — Key Insight from Real Simulation Data
# ══════════════════════════════════════════════════════════════════════════════
def fig6_key_insight():
    data_a = np.load("/Users/aayanalwani/tb project/mce3r_stochastic/results/phase2/condition_A.npz")
    data_b = np.load("/Users/aayanalwani/tb project/mce3r_stochastic/results/phase2/condition_B.npz")
    prot_a = data_a['proteins']
    prot_b = data_b['proteins']

    mean_a, std_a = np.mean(prot_a), np.std(prot_a)
    mean_b, std_b = np.mean(prot_b), np.std(prot_b)
    cv_a = std_a / mean_a if mean_a > 0 else 0
    cv_b = std_b / mean_b if mean_b > 0 else 0

    fig, ax = plt.subplots(figsize=(10, 6), constrained_layout=True)

    bins = np.linspace(0, 600, 51)

    ax.hist(prot_a, bins=bins, density=True, histtype='step', lw=2.5,
            color=TEAL, alpha=0.9, label=f'Asymmetric (A): CV = {cv_a:.3f}')
    ax.hist(prot_a, bins=bins, density=True, alpha=0.25, color=TEAL)

    ax.hist(prot_b, bins=bins, density=True, histtype='step', lw=2.5,
            color=CORAL, alpha=0.9, label=f'Symmetric (B): CV = {cv_b:.3f}')
    ax.hist(prot_b, bins=bins, density=True, alpha=0.25, color=CORAL)

    # Mean lines
    ax.axvline(mean_a, color=TEAL, ls='--', lw=2, alpha=0.8)
    ax.axvline(mean_b, color=CORAL, ls='--', lw=2, alpha=0.8)
    ax.text(mean_a + 3, ax.get_ylim()[1] * 0.92, f'Mean A = {mean_a:.0f}',
            fontsize=11, color='#065F56')
    ax.text(mean_b + 3, ax.get_ylim()[1] * 0.82, f'Mean B = {mean_b:.0f}',
            fontsize=11, color='#C0392B')

    # Persister zone shading
    ax.axvspan(0, 50, alpha=0.15, color='red', zorder=0)
    ax.text(25, ax.get_ylim()[1] * 0.70, "PERSISTER\nZONE",
            ha='center', fontsize=12, fontweight='bold', color='#B91C1C',
            rotation=90)

    # CV annotations
    ax.text(0.98, 0.95, f"CV(asym) = {cv_a:.3f}\nCV(sym)  = {cv_b:.3f}",
            transform=ax.transAxes, ha='right', va='top', fontsize=13,
            fontweight='bold',
            bbox=dict(boxstyle='round,pad=0.5', facecolor='#FFF7ED',
                      edgecolor='#F97316', lw=1.5))

    # Tail annotation
    frac_a = np.mean(prot_a < 50) * 100
    frac_b = np.mean(prot_b < 50) * 100
    ax.annotate(f"More cells in\nlow-expression tail\n\u2192 more persisters\n"
                f"({frac_a:.1f}% vs {frac_b:.1f}%)",
                xy=(40, ax.get_ylim()[1] * 0.25),
                xytext=(150, ax.get_ylim()[1] * 0.50),
                fontsize=12, fontweight='bold', color='#991B1B',
                arrowprops=dict(arrowstyle='->', lw=2, color='#991B1B'),
                bbox=dict(boxstyle='round,pad=0.4', facecolor='#FEF2F2',
                          edgecolor='#991B1B', lw=1.2))

    ax.set_xlabel("Protein count per cell", fontsize=13)
    ax.set_ylabel("Probability density", fontsize=13)
    ax.set_title("The Key Finding: Asymmetry Creates a Wider Distribution with More Persisters",
                 fontsize=14, fontweight='bold', pad=10)
    ax.legend(fontsize=12, loc='upper right',
              bbox_to_anchor=(0.99, 0.75))

    path = OUTDIR / "asymmetry_key_insight.png"
    fig.savefig(path, dpi=300, bbox_inches='tight')
    plt.close(fig)
    return path


# ══════════════════════════════════════════════════════════════════════════════
# Main
# ══════════════════════════════════════════════════════════════════════════════
if __name__ == "__main__":
    print("Generating noise explainer figures ...\n")

    generators = [
        ("Figure 1", fig1_what_is_noise),
        ("Figure 2", fig2_persistence),
        ("Figure 3", fig3_cv_explained),
        ("Figure 4", fig4_operator_mechanism),
        ("Figure 5", fig5_pipeline),
        ("Figure 6", fig6_key_insight),
    ]

    for name, func in generators:
        path = func()
        size_kb = os.path.getsize(path) / 1024
        print(f"  {name}: {path.name:<45s} ({size_kb:>7.1f} KB)")

    print(f"\nAll 6 figures saved to:\n  {OUTDIR}")
    print("Done.")
