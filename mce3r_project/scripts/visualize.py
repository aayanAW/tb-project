"""
visualize.py — Generate publication-quality figures for Mce3R motif analysis.

Produces six figures summarizing the motif discovery and analysis results:

1. score_distribution          — Distribution of FIMO log-odds scores by strand
2. motif_position_distribution — Motif positions along each sequence
3. spacing_histogram           — Gap between adjacent binding sites
4. palindrome_scores           — Symmetry scores for each hit
5. gc_vs_hits                  — Sequence GC content vs. number of hits
6. motif_logo                  — Position-frequency logo (stacked bar approximation)

All figures are saved as both PNG (300 DPI for print) and SVG (for editing).

Biological interpretation guidance is included in each plot's title/annotations
to make figures self-explanatory for bioinformatics papers.

Usage:
    python visualize.py [--hits PATH] [--spacing PATH] [--figures-dir PATH]
"""

import sys
import warnings
from pathlib import Path

import matplotlib
matplotlib.use("Agg")  # Non-interactive backend for script mode

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd
import seaborn as sns

sys.path.insert(0, str(Path(__file__).resolve().parent))
from utils import get_project_root, setup_logging

logger = setup_logging(__name__)

# Consistent color palette for all figures
BASE_COLORS = {"A": "#2ecc71", "C": "#3498db", "G": "#f39c12", "T": "#e74c3c"}
STRAND_COLORS = {"+": "#2980b9", "-": "#c0392b"}


def _apply_style():
    """Apply a consistent matplotlib style across all figures."""
    style_options = ["seaborn-v0_8-whitegrid", "seaborn-whitegrid", "ggplot"]
    for style in style_options:
        try:
            plt.style.use(style)
            break
        except OSError:
            continue


def _save_figure(fig, output_path: Path, dpi: int = 300):
    """Save figure as both PNG and SVG."""
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(str(output_path.with_suffix(".png")), dpi=dpi, bbox_inches="tight")
    fig.savefig(str(output_path.with_suffix(".svg")), bbox_inches="tight")
    logger.info(f"Saved figure: {output_path.with_suffix('.png')}")
    plt.close(fig)


def plot_score_distribution(hits_df: pd.DataFrame, output_path: Path) -> None:
    """
    Plot the distribution of FIMO log-odds scores, colored by strand.

    A higher score indicates a better match to the Mce3R binding motif.
    The distribution shape reveals whether the motif model discriminates
    well: a bimodal or right-skewed distribution is desirable.

    Args:
        hits_df: DataFrame with 'score' and 'strand' columns
        output_path: Base path for output (extensions added automatically)
    """
    _apply_style()
    fig, ax = plt.subplots(figsize=(10, 6))

    if hits_df.empty or "score" not in hits_df.columns:
        ax.text(0.5, 0.5, "No FIMO hits to plot", transform=ax.transAxes,
                ha="center", va="center", fontsize=14)
        ax.set_title("FIMO Score Distribution")
        _save_figure(fig, output_path)
        return

    for strand, color in STRAND_COLORS.items():
        subset = hits_df[hits_df["strand"] == strand]["score"].dropna()
        if len(subset) > 0:
            sns.histplot(
                subset,
                ax=ax,
                color=color,
                alpha=0.6,
                label=f"Strand {strand} (n={len(subset)})",
                bins=20,
                kde=len(subset) > 5,
            )

    ax.set_xlabel("FIMO Log-Odds Score", fontsize=13)
    ax.set_ylabel("Count", fontsize=13)
    ax.set_title(
        "Distribution of Mce3R Motif Scores\n"
        "(Higher score = stronger match to binding motif)",
        fontsize=13,
    )
    ax.legend(fontsize=11)

    # Annotate mean score
    mean_score = hits_df["score"].mean()
    ax.axvline(mean_score, color="black", linestyle="--", alpha=0.7, label=f"Mean: {mean_score:.2f}")
    ax.text(mean_score + 0.1, ax.get_ylim()[1] * 0.9, f"Mean={mean_score:.2f}",
            fontsize=10, color="black")

    _save_figure(fig, output_path)


def plot_position_distribution(
    hits_df: pd.DataFrame,
    sequences_df: pd.DataFrame,
    output_path: Path,
) -> None:
    """
    Plot motif hit positions along each sequence (lollipop / scatter plot).

    Shows whether Mce3R binding sites occur preferentially at certain
    positions relative to the start of each promoter region.

    Args:
        hits_df: DataFrame with 'sequence_name', 'start', 'stop', 'strand'
        sequences_df: DataFrame with 'sequence_id', 'length'
        output_path: Base output path
    """
    _apply_style()

    if hits_df.empty:
        fig, ax = plt.subplots(figsize=(10, 6))
        ax.text(0.5, 0.5, "No FIMO hits to plot", transform=ax.transAxes,
                ha="center", va="center", fontsize=14)
        ax.set_title("Motif Position Distribution")
        _save_figure(fig, output_path)
        return

    # Compute relative positions (0-1 scale)
    plot_df = hits_df.copy()
    if not sequences_df.empty and "sequence_id" in sequences_df.columns:
        len_map = sequences_df.set_index("sequence_id")["length"].to_dict()
        plot_df["rel_position"] = plot_df.apply(
            lambda r: r["start"] / len_map.get(r["sequence_name"], r["start"] + 1), axis=1
        )
    else:
        plot_df["rel_position"] = plot_df["start"]

    # Sort sequences for a consistent y-axis order
    seq_order = sorted(plot_df["sequence_name"].unique())
    seq_y = {name: i for i, name in enumerate(seq_order)}
    plot_df["y"] = plot_df["sequence_name"].map(seq_y)

    fig, ax = plt.subplots(figsize=(12, max(6, len(seq_order) * 0.4 + 2)))

    for strand, color in STRAND_COLORS.items():
        subset = plot_df[plot_df["strand"] == strand]
        if subset.empty:
            continue
        ax.scatter(
            subset["rel_position"],
            subset["y"],
            c=color,
            alpha=0.8,
            label=f"Strand {strand}",
            s=60,
            zorder=3,
        )

    ax.set_yticks(list(seq_y.values()))
    ax.set_yticklabels(list(seq_y.keys()), fontsize=8)
    ax.set_xlabel("Relative Position in Sequence (0 = start, 1 = end)", fontsize=12)
    ax.set_ylabel("Sequence", fontsize=12)
    ax.set_title(
        "Mce3R Motif Positions Across Promoter Sequences\n"
        "(Each point = one FIMO-predicted binding site)",
        fontsize=13,
    )
    ax.legend(fontsize=11)
    ax.set_xlim(-0.05, 1.05)

    _save_figure(fig, output_path)


def plot_spacing_histogram(spacing_df: pd.DataFrame, output_path: Path) -> None:
    """
    Histogram of inter-site spacing between adjacent Mce3R operators.

    Shades the 40–65 bp validated spacing range (Panagoda et al. 2024) and
    marks the exact 53 bp spacer from the proven mce3R–yrbE3A operator.

    Args:
        spacing_df: DataFrame from compute_inter_site_spacing()
        output_path: Base output path
    """
    _apply_style()
    fig, ax = plt.subplots(figsize=(10, 6))

    if spacing_df.empty or "spacing_bp" not in spacing_df.columns:
        ax.text(
            0.5, 0.5,
            "No multi-hit sequences found.\n"
            "Each sequence had at most one Mce3R binding site.",
            transform=ax.transAxes, ha="center", va="center", fontsize=12,
        )
        ax.set_title("Inter-Site Spacing Between Adjacent Mce3R Operators")
        _save_figure(fig, output_path)
        return

    spacings = spacing_df["spacing_bp"].dropna()
    ax.hist(spacings, bins=max(5, len(spacings) // 3), color="#8e44ad", alpha=0.75, edgecolor="white")

    # Highlight validated 40–65 bp range (2024 ACS Chemical Biology)
    ax.axvspan(40, 65, alpha=0.15, color="#27ae60",
               label="Validated Mce3R operator spacing (40–65 bp)")

    # Mark the exact 53 bp spacer from the proven mce3R–yrbE3A operator
    ax.axvline(53, color="#27ae60", linestyle="--", linewidth=2.0,
               label="Validated spacer: 53 bp (yrbE3A operator)")

    ax.set_xlabel("Inter-Site Spacing (bp)", fontsize=13)
    ax.set_ylabel("Count", fontsize=13)
    ax.set_title(
        "Spacing Between Adjacent Mce3R Binding Sites\n"
        "Validated operator spacer: 53 bp (2024 ACS Chem. Biol.)",
        fontsize=13,
    )
    ax.legend(fontsize=11)

    if len(spacings) > 0:
        n_validated = ((spacings >= 40) & (spacings <= 65)).sum()
        n_exact = ((spacings >= 48) & (spacings <= 58)).sum()
        ax.text(
            0.98, 0.95,
            f"n={len(spacings)} pairs\nMean={spacings.mean():.1f} bp\n"
            f"Median={spacings.median():.1f} bp\n"
            f"In 40–65 bp window: {n_validated}\nNear 53 bp (±5 bp): {n_exact}",
            transform=ax.transAxes, ha="right", va="top", fontsize=9,
            bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8),
        )

    _save_figure(fig, output_path)


def plot_palindrome_scores(hits_df: pd.DataFrame, output_path: Path) -> None:
    """
    Violin + strip plot of palindrome symmetry scores for all FIMO hits.

    A palindrome score >= 0.8 indicates a site consistent with the
    homodimeric TetR binding model (each monomer binds one half-site
    on opposite strands).

    Args:
        hits_df: DataFrame with 'palindrome_score' and 'motif_id' columns
        output_path: Base output path
    """
    _apply_style()
    fig, ax = plt.subplots(figsize=(9, 6))

    if hits_df.empty or "palindrome_score" not in hits_df.columns:
        ax.text(0.5, 0.5, "No palindrome data available", transform=ax.transAxes,
                ha="center", va="center", fontsize=14)
        ax.set_title("Palindrome Symmetry Scores")
        _save_figure(fig, output_path)
        return

    motif_ids = hits_df["motif_id"].unique().tolist()
    plot_data = [hits_df[hits_df["motif_id"] == mid]["palindrome_score"].dropna().values
                 for mid in motif_ids]

    # Violin plot
    parts = ax.violinplot(plot_data, positions=range(len(motif_ids)), showmedians=True)
    for pc in parts["bodies"]:
        pc.set_facecolor("#2980b9")
        pc.set_alpha(0.6)

    # Strip plot (individual points)
    for i, (mid, data) in enumerate(zip(motif_ids, plot_data)):
        jitter = np.random.default_rng(42).uniform(-0.1, 0.1, size=len(data))
        ax.scatter(np.full(len(data), i) + jitter, data, alpha=0.5, s=25, color="#c0392b", zorder=3)

    # Palindrome threshold line
    ax.axhline(0.8, color="#e67e22", linestyle="--", linewidth=1.5,
               label="Palindrome threshold (score = 0.8)")

    ax.set_xticks(range(len(motif_ids)))
    ax.set_xticklabels(motif_ids, fontsize=11)
    ax.set_xlabel("Motif ID", fontsize=13)
    ax.set_ylabel("Palindrome Score (0=asymmetric, 1=perfect palindrome)", fontsize=12)
    ax.set_title(
        "Palindromic Symmetry of Predicted Mce3R Binding Sites\n"
        "(Mce3R uses ASYMMETRIC operator — expected range 0.3–0.65, NOT ≥ 0.8)",
        fontsize=12,
    )
    ax.set_ylim(-0.05, 1.1)
    ax.legend(fontsize=10)

    # Add shading for the expected Mce3R asymmetric range
    ax.axhspan(0.3, 0.65, alpha=0.08, color="#27ae60")
    ax.text(len(hits_df["motif_id"].unique()) - 0.5, 0.475,
            "Expected\nMce3R range\n(0.3–0.65)",
            ha="right", va="center", fontsize=7, color="#27ae60", alpha=0.9)

    n_palindromes = (hits_df["palindrome_score"] >= 0.8).sum()
    n_likely_non = (hits_df["palindrome_score"] >= 0.85).sum() if "palindrome_score" in hits_df.columns else 0
    n_in_target = ((hits_df["palindrome_score"] >= 0.3) &
                   (hits_df["palindrome_score"] <= 0.65)).sum() if "palindrome_score" in hits_df.columns else 0
    ax.text(
        0.02, 0.05,
        f"In Mce3R range (0.3–0.65): {n_in_target}/{len(hits_df)}\n"
        f"Palindromic (≥ 0.8): {n_palindromes}/{len(hits_df)}\n"
        f"Non-Mce3R TetR flag (≥ 0.85): {n_likely_non}",
        transform=ax.transAxes, fontsize=9,
        bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8),
    )

    _save_figure(fig, output_path)


def plot_gc_vs_hits(annotated_df: pd.DataFrame, output_path: Path) -> None:
    """
    Scatter plot: sequence GC content vs. number of FIMO hits in that sequence.

    Colors points by ground-truth motif status (has_motif). Reveals whether
    the motif scanner is biased toward high-GC sequences (a common artifact
    in GC-rich genomes like M. tuberculosis).

    Args:
        annotated_df: DataFrame from annotate_hits_with_sequence_metadata()
        output_path: Base output path
    """
    _apply_style()
    fig, ax = plt.subplots(figsize=(9, 6))

    if annotated_df.empty or "gc_content" not in annotated_df.columns:
        ax.text(0.5, 0.5, "No annotated data available", transform=ax.transAxes,
                ha="center", va="center", fontsize=14)
        ax.set_title("GC Content vs. Number of Motif Hits")
        _save_figure(fig, output_path)
        return

    # Aggregate to per-sequence level
    agg_dict = {"n_hits": ("sequence_name", "count")}
    if "gc_content" in annotated_df.columns:
        agg_dict["gc_content"] = ("gc_content", "first")
    if "has_motif" in annotated_df.columns:
        agg_dict["has_motif"] = ("has_motif", "first")

    agg = annotated_df.groupby("sequence_name").agg(**agg_dict).reset_index()

    if "gc_content" not in agg.columns:
        logger.warning("gc_content not available for GC vs hits plot — skipping figure")
        fig, ax = plt.subplots(figsize=(9, 6))
        ax.text(0.5, 0.5, "GC content data not available\n(requires promoter metadata join)",
                transform=ax.transAxes, ha="center", va="center", fontsize=12)
        ax.set_title("GC Content vs. Number of Motif Hits")
        _save_figure(fig, output_path)
        return

    if "has_motif" in agg.columns:
        color_map = {True: "#27ae60", False: "#e74c3c", None: "#95a5a6"}
        label_map = {True: "True positive (motif embedded)", False: "False positive", None: "Unknown"}
        for has_motif, color in color_map.items():
            subset = agg[agg["has_motif"] == has_motif]
            if subset.empty:
                continue
            ax.scatter(
                subset["gc_content"], subset["n_hits"],
                c=color, alpha=0.8, s=80, label=label_map[has_motif], zorder=3,
            )
    else:
        # Real genome data: no ground truth — color by hit count intensity
        sc = ax.scatter(
            agg["gc_content"], agg["n_hits"],
            c=agg["n_hits"], cmap="Blues", alpha=0.8, s=80, zorder=3,
        )
        plt.colorbar(sc, ax=ax, label="Number of hits")

    ax.set_xlabel("Sequence GC Content", fontsize=13)
    ax.set_ylabel("Number of FIMO Hits", fontsize=13)
    ax.set_title(
        "Sequence GC Content vs. Number of Predicted Binding Sites\n"
        "(Check for GC bias in motif scanner)",
        fontsize=13,
    )
    if "has_motif" in agg.columns:
        ax.legend(fontsize=11)

    _save_figure(fig, output_path)


def plot_sequence_logo(hits_df: pd.DataFrame, output_path: Path) -> None:
    """
    Approximate sequence logo as a stacked bar chart of base frequencies.

    For each position in the matched sequences, plots the frequency of
    each nucleotide as a colored segment. Height proportional to frequency.
    This approximates a sequence logo without the logomaker dependency.

    Args:
        hits_df: DataFrame with 'matched_sequence' column
        output_path: Base output path
    """
    _apply_style()
    fig, ax = plt.subplots(figsize=(14, 5))

    if hits_df.empty or "matched_sequence" not in hits_df.columns:
        ax.text(0.5, 0.5, "No matched sequences available", transform=ax.transAxes,
                ha="center", va="center", fontsize=14)
        ax.set_title("Mce3R Binding Site Sequence Logo")
        _save_figure(fig, output_path)
        return

    # Filter to non-empty matched sequences
    seqs = [str(s).upper() for s in hits_df["matched_sequence"].dropna() if len(str(s)) > 0]

    if not seqs:
        ax.text(0.5, 0.5, "No matched sequences available", transform=ax.transAxes,
                ha="center", va="center")
        _save_figure(fig, output_path)
        return

    # Pad to max length (truncate if needed)
    max_len = max(len(s) for s in seqs)
    padded = [s[:max_len].ljust(max_len, "N") for s in seqs]
    n_seqs = len(padded)

    # Compute frequency matrix (positions x 4 bases)
    bases = ["A", "C", "G", "T"]
    freq_matrix = np.zeros((max_len, 4))
    for seq in padded:
        for i, base in enumerate(seq[:max_len]):
            idx = bases.index(base) if base in bases else -1
            if idx >= 0:
                freq_matrix[i, idx] += 1
    freq_matrix /= n_seqs

    # Stacked bar chart
    positions = np.arange(1, max_len + 1)
    bottoms = np.zeros(max_len)

    for j, base in enumerate(bases):
        heights = freq_matrix[:, j]
        ax.bar(
            positions,
            heights,
            bottom=bottoms,
            color=BASE_COLORS[base],
            label=base,
            width=0.85,
            edgecolor="none",
        )
        bottoms += heights

    ax.set_xlabel("Position in Binding Site (bp)", fontsize=13)
    ax.set_ylabel("Base Frequency", fontsize=13)
    ax.set_title(
        f"Mce3R Binding Site Sequence Logo (n={n_seqs} sites)\n"
        "Stacked bars show nucleotide frequencies at each position",
        fontsize=13,
    )
    ax.set_xlim(0.5, max_len + 0.5)
    ax.set_ylim(0, 1.05)
    ax.set_xticks(positions)
    ax.legend(title="Base", fontsize=11, loc="upper right")

    # Mark half-sites with vertical lines at midpoint
    mid = max_len / 2
    ax.axvline(mid + 0.5, color="gray", linestyle=":", alpha=0.7, linewidth=1.5)
    ax.text(mid * 0.5, 1.02, "Left half-site", ha="center", fontsize=9, color="gray")
    ax.text(mid + mid * 0.5 + 0.5, 1.02, "Right half-site", ha="center", fontsize=9, color="gray")

    _save_figure(fig, output_path)


def plot_operator_architecture_diagram(hits_df: pd.DataFrame, output_path: Path) -> None:
    """
    Schematic diagram of the validated Mce3R tandem two-site operator architecture.

    Shows the reference geometry (two ~25 bp asymmetric half-sites, 53 bp spacer)
    alongside actual detected architectures from the FIMO scan. This figure
    communicates the 2024 ACS Chemical Biology finding that Mce3R uses an
    asymmetric paired-site operator rather than a palindromic inverted repeat.

    Args:
        hits_df: DataFrame with 'site_architecture' and 'partner_spacing_bp' columns
        output_path: Base output path
    """
    _apply_style()
    fig, axes = plt.subplots(2, 1, figsize=(12, 8),
                             gridspec_kw={"height_ratios": [1, 1.5]})

    # ── Panel 1: Reference operator diagram ──────────────────────────────────
    ax_diag = axes[0]
    ax_diag.set_xlim(0, 250)
    ax_diag.set_ylim(-1, 2)
    ax_diag.axis("off")

    # Draw the two asymmetric half-sites as directional boxes
    # Upstream site (lower affinity): positions 10–60
    up_start, up_end = 20, 65
    down_start, down_end = 120, 165
    spacer_start, spacer_end = 65, 120

    # Upstream half-site box
    rect_up = plt.Rectangle((up_start, 0.3), up_end - up_start, 0.7,
                             facecolor="#3498db", edgecolor="#2980b9",
                             linewidth=1.5, alpha=0.85)
    ax_diag.add_patch(rect_up)
    ax_diag.annotate("", xy=(up_end + 2, 0.65), xytext=(up_end, 0.65),
                     arrowprops=dict(arrowstyle="->", color="#2980b9", lw=2))
    ax_diag.text((up_start + up_end) / 2, 0.65, "Upstream\nhalf-site\n~25 bp",
                 ha="center", va="center", fontsize=8, color="white", fontweight="bold")
    ax_diag.text((up_start + up_end) / 2, 0.1, "Lower affinity", ha="center",
                 va="top", fontsize=7, color="#2980b9")

    # Spacer bracket
    spacer_mid = (spacer_start + spacer_end) / 2
    ax_diag.annotate("", xy=(spacer_end, 1.2), xytext=(spacer_start, 1.2),
                     arrowprops=dict(arrowstyle="<->", color="#7f8c8d", lw=1.5))
    ax_diag.text(spacer_mid, 1.32, "~53 bp spacer", ha="center", va="bottom",
                 fontsize=8, color="#7f8c8d", style="italic")
    ax_diag.text(spacer_mid, -0.05, "Non-contacted\nDNA", ha="center", va="top",
                 fontsize=7, color="#95a5a6")

    # Downstream half-site box (higher affinity)
    rect_down = plt.Rectangle((down_start, 0.3), down_end - down_start, 0.7,
                               facecolor="#e74c3c", edgecolor="#c0392b",
                               linewidth=1.5, alpha=0.85)
    ax_diag.add_patch(rect_down)
    ax_diag.annotate("", xy=(down_end + 2, 0.65), xytext=(down_end, 0.65),
                     arrowprops=dict(arrowstyle="->", color="#c0392b", lw=2))
    ax_diag.text((down_start + down_end) / 2, 0.65, "Downstream\nhalf-site\n~25 bp",
                 ha="center", va="center", fontsize=8, color="white", fontweight="bold")
    ax_diag.text((down_start + down_end) / 2, 0.1, "Higher affinity\nKd ≈ 2.4 nM",
                 ha="center", va="top", fontsize=7, color="#c0392b")

    ax_diag.set_title(
        "Validated Mce3R Operator Architecture (mce3R–yrbE3A, H37Rv)\n"
        "Panagoda, Balázsi & Sampson, ACS Chem. Biol. 2024  |  "
        "Both sites asymmetric (nonpalindromic) — same-strand direct repeat",
        fontsize=10, pad=4,
    )

    # ── Panel 2: Detected architecture distribution + spacing ─────────────────
    ax_bar = axes[1]

    if hits_df.empty or "site_architecture" not in hits_df.columns:
        ax_bar.text(0.5, 0.5, "No architecture data available",
                    transform=ax_bar.transAxes, ha="center", va="center", fontsize=12)
        ax_bar.set_title("Detected Site Architectures in This Scan")
    else:
        arch_counts = hits_df["site_architecture"].value_counts()
        arch_color_map = {
            "validated_paired_asymmetric": "#27ae60",
            "tandem_asymmetric":           "#2ecc71",
            "direct_repeat":               "#f39c12",
            "convergent":                  "#e67e22",
            "inverted_repeat":             "#e74c3c",
            "divergent":                   "#9b59b6",
            "isolated":                    "#95a5a6",
        }
        colors = [arch_color_map.get(a, "#bdc3c7") for a in arch_counts.index]
        bars = ax_bar.bar(arch_counts.index, arch_counts.values, color=colors,
                          edgecolor="white", linewidth=0.8, width=0.6)

        for bar, count in zip(bars, arch_counts.values):
            ax_bar.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.2,
                        str(count), ha="center", va="bottom", fontsize=9, fontweight="bold")

        # Add spacing distribution for paired sites if available
        if "partner_spacing_bp" in hits_df.columns:
            valid_spacings = hits_df["partner_spacing_bp"].dropna()
            valid_spacings = valid_spacings[valid_spacings >= 0]
            if not valid_spacings.empty:
                n_in_window = ((valid_spacings >= 40) & (valid_spacings <= 65)).sum()
                ax_bar.text(
                    0.98, 0.95,
                    f"Paired sites with 40–65 bp spacer: {n_in_window}\n"
                    f"(validated_paired_asymmetric arch)",
                    transform=ax_bar.transAxes, ha="right", va="top", fontsize=9,
                    bbox=dict(boxstyle="round,pad=0.3", facecolor="#27ae60", alpha=0.15),
                )

        ax_bar.set_xlabel("Architecture Type", fontsize=11)
        ax_bar.set_ylabel("Number of Predicted Sites", fontsize=11)
        ax_bar.set_title(
            "Detected Site Architectures in This Scan\n"
            "(validated_paired_asymmetric = 40–65 bp same-strand pair; HIGHEST priority)",
            fontsize=10,
        )
        ax_bar.set_xticks(range(len(arch_counts)))
        ax_bar.set_xticklabels(arch_counts.index, rotation=18, ha="right", fontsize=9)

    fig.tight_layout(pad=1.5)
    _save_figure(fig, output_path)


def plot_model_comparison(hits_df: pd.DataFrame, output_path: Path) -> None:
    """
    Scatter plot of palindrome score vs. FIMO score, colored by model preference.

    Each point is a predicted binding site. Sites consistent with the 2024 ACS
    Chemical Biology asymmetric Mce3R model have low palindrome scores and high
    FIMO scores. Sites consistent with the classic palindromic TetR model have
    high palindrome scores.

    Args:
        hits_df: DataFrame with 'palindrome_score', 'score', 'model_preference' columns
        output_path: Base output path
    """
    _apply_style()
    fig, ax = plt.subplots(figsize=(10, 7))

    required_cols = {"palindrome_score", "score", "model_preference"}
    if hits_df.empty or not required_cols.issubset(hits_df.columns):
        ax.text(
            0.5, 0.5,
            "Model comparison data not available.\n"
            "Run the analyze step with motif_models.py present.",
            transform=ax.transAxes, ha="center", va="center", fontsize=12,
        )
        ax.set_title("Model Comparison: Palindromic vs. Asymmetric Mce3R Binding Sites")
        _save_figure(fig, output_path)
        return

    model_colors = {
        "palindromic": "#e74c3c",
        "asymmetric":  "#27ae60",
        "ambiguous":   "#f39c12",
    }
    model_labels = {
        "palindromic": "Palindromic (classic TetR — less likely for Mce3R)",
        "asymmetric":  "Asymmetric (2024 ACS Chem. Biol. model)",
        "ambiguous":   "Ambiguous (palindrome score 0.6–0.8)",
    }

    for model_pref, color in model_colors.items():
        subset = hits_df[hits_df["model_preference"] == model_pref]
        if subset.empty:
            continue
        ax.scatter(
            subset["palindrome_score"],
            subset["score"],
            c=color,
            alpha=0.75,
            s=70,
            label=f"{model_labels[model_pref]} (n={len(subset)})",
            zorder=3,
        )

    # Reference lines for model thresholds
    ax.axvline(0.8, color="#e74c3c", linestyle="--", alpha=0.5, linewidth=1.2,
               label="Palindrome threshold (0.80)")
    ax.axvline(0.6, color="#27ae60", linestyle="--", alpha=0.5, linewidth=1.2,
               label="Asymmetric threshold (0.60)")

    # Shade asymmetric zone
    ax.axvspan(0, 0.6, alpha=0.05, color="#27ae60", label="")
    ax.axvspan(0.8, 1.0, alpha=0.05, color="#e74c3c", label="")

    ax.set_xlabel("Palindrome Score (0 = fully asymmetric, 1 = perfect palindrome)", fontsize=12)
    ax.set_ylabel("FIMO Log-Odds Score (higher = stronger motif match)", fontsize=12)
    ax.set_title(
        "Palindromic vs. Asymmetric Mce3R Binding Site Predictions\n"
        "Based on 2024 ACS Chemical Biology nonpalindromic operator model\n"
        "[COMPUTATIONAL PREDICTIONS — experimental validation required]",
        fontsize=12,
    )
    ax.set_xlim(-0.05, 1.05)
    ax.legend(fontsize=9, loc="upper left")

    # Annotation box
    n_asym = (hits_df["model_preference"] == "asymmetric").sum()
    n_pal = (hits_df["model_preference"] == "palindromic").sum()
    ax.text(
        0.98, 0.02,
        f"Asymmetric: {n_asym} sites\nPalindromic: {n_pal} sites\n"
        f"Total: {len(hits_df)} sites",
        transform=ax.transAxes, ha="right", va="bottom", fontsize=9,
        bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.85),
    )

    _save_figure(fig, output_path)


def plot_site_architecture(hits_df: pd.DataFrame, output_path: Path) -> None:
    """
    Bar chart of site architecture type distribution.

    Architecture types characterize the spatial relationship between paired
    Mce3R binding sites (tandem, inverted repeat, convergent, isolated, etc.).
    The 'tandem_asymmetric' architecture (two direct-repeat sites, 0–60 bp apart)
    is the strongest computational evidence for the nonpalindromic Mce3R model.

    Args:
        hits_df: DataFrame with 'site_architecture' column
        output_path: Base output path
    """
    _apply_style()
    fig, ax = plt.subplots(figsize=(10, 6))

    if hits_df.empty or "site_architecture" not in hits_df.columns:
        ax.text(
            0.5, 0.5,
            "Site architecture data not available.\n"
            "Run the analyze step with motif_models.py present.",
            transform=ax.transAxes, ha="center", va="center", fontsize=12,
        )
        ax.set_title("Site Architecture Distribution")
        _save_figure(fig, output_path)
        return

    arch_counts = hits_df["site_architecture"].value_counts()

    # Color by biological relevance for Mce3R asymmetric model
    arch_color_map = {
        "tandem_asymmetric": "#27ae60",   # strongest Mce3R model evidence
        "direct_repeat":     "#2ecc71",   # same-strand, wide spacing
        "convergent":        "#f39c12",   # head-to-head
        "inverted_repeat":   "#e74c3c",   # palindromic (less likely for Mce3R)
        "divergent":         "#9b59b6",   # unusual
        "isolated":          "#95a5a6",   # no nearby partner
    }
    colors = [arch_color_map.get(a, "#bdc3c7") for a in arch_counts.index]

    bars = ax.bar(arch_counts.index, arch_counts.values, color=colors, edgecolor="white",
                  linewidth=0.8, width=0.6)

    # Count labels on bars
    for bar, count in zip(bars, arch_counts.values):
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            bar.get_height() + 0.1,
            str(count),
            ha="center", va="bottom", fontsize=10, fontweight="bold",
        )

    # Legend patches
    patches = [
        mpatches.Patch(color="#27ae60", label="tandem_asymmetric — Strongest Mce3R model evidence"),
        mpatches.Patch(color="#2ecc71", label="direct_repeat — Same strand, wide spacing"),
        mpatches.Patch(color="#f39c12", label="convergent — Head-to-head pair"),
        mpatches.Patch(color="#e74c3c", label="inverted_repeat — Palindromic (classic TetR)"),
        mpatches.Patch(color="#9b59b6", label="divergent — Unusual geometry"),
        mpatches.Patch(color="#95a5a6", label="isolated — No nearby partner"),
    ]
    ax.legend(handles=patches, fontsize=8, loc="upper right")

    ax.set_xlabel("Site Architecture Type", fontsize=12)
    ax.set_ylabel("Number of Predicted Sites", fontsize=12)
    ax.set_title(
        "Mce3R Binding Site Architecture Distribution\n"
        "(tandem_asymmetric = direct-repeat pair, 0–60 bp, strongest 2024 model support)\n"
        "[COMPUTATIONAL PREDICTIONS — experimental validation required]",
        fontsize=12,
    )
    ax.set_xticks(range(len(arch_counts)))
    ax.set_xticklabels(arch_counts.index, rotation=20, ha="right", fontsize=10)

    _save_figure(fig, output_path)


def generate_all_figures(
    hits_df: pd.DataFrame,
    sequences_df: pd.DataFrame,
    spacing_df: pd.DataFrame,
    figures_dir: Path,
) -> list:
    """
    Generate all six figures for the Mce3R motif analysis.

    Args:
        hits_df: Annotated FIMO hits DataFrame
        sequences_df: Sequence metadata DataFrame
        spacing_df: Inter-site spacing DataFrame
        figures_dir: Directory to save all figures

    Returns:
        List of paths to generated PNG files
    """
    figures_dir = Path(figures_dir)
    figures_dir.mkdir(parents=True, exist_ok=True)

    generated = []

    figures = [
        (plot_score_distribution, [hits_df], "score_distribution"),
        (plot_position_distribution, [hits_df, sequences_df], "motif_position_distribution"),
        (plot_spacing_histogram, [spacing_df], "spacing_histogram"),
        (plot_palindrome_scores, [hits_df], "palindrome_scores"),
        (plot_gc_vs_hits, [hits_df], "gc_vs_hits"),
        (plot_sequence_logo, [hits_df], "motif_logo"),
        (plot_model_comparison, [hits_df], "model_comparison"),
        (plot_site_architecture, [hits_df], "site_architecture"),
        (plot_operator_architecture_diagram, [hits_df], "operator_architecture_diagram"),
    ]

    for func, args, name in figures:
        output_path = figures_dir / name
        try:
            func(*args, output_path)
            generated.append(figures_dir / f"{name}.png")
        except Exception as e:
            logger.error(f"Failed to generate {name}: {e}")

    logger.info(f"Generated {len(generated)}/9 figures in {figures_dir}")
    return generated


def main():
    import argparse

    root = get_project_root()
    parser = argparse.ArgumentParser(
        description="Generate visualization figures for Mce3R motif analysis"
    )
    parser.add_argument(
        "--hits",
        type=Path,
        default=root / "data" / "processed" / "annotated_hits.csv",
        help="Annotated FIMO hits CSV",
    )
    parser.add_argument(
        "--metadata",
        type=Path,
        default=root / "data" / "processed" / "sequence_metadata.csv",
        help="Sequence metadata CSV",
    )
    parser.add_argument(
        "--spacing",
        type=Path,
        default=root / "data" / "processed" / "spacing_analysis.csv",
        help="Spacing analysis CSV",
    )
    parser.add_argument(
        "--figures-dir",
        type=Path,
        default=root / "results" / "figures",
        help="Output directory for figures",
    )
    args = parser.parse_args()

    # Load DataFrames
    hits_df = pd.read_csv(args.hits) if args.hits.exists() else pd.DataFrame()
    sequences_df = pd.read_csv(args.metadata) if args.metadata.exists() else pd.DataFrame()
    spacing_df = pd.read_csv(args.spacing) if args.spacing.exists() else pd.DataFrame()

    generated = generate_all_figures(hits_df, sequences_df, spacing_df, args.figures_dir)

    print(f"\nGenerated {len(generated)} figures:")
    for path in generated:
        print(f"  {path}")


if __name__ == "__main__":
    main()
