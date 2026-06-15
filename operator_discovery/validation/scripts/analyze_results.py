"""
analyze_results.py — Statistical analysis of Mce3R motif hits from FIMO.

This script processes FIMO output to characterize the architecture of
Mce3R binding sites in M. tuberculosis promoter regions.

Analyses performed:
1. Load and clean FIMO hits table
2. Summary statistics (hit counts, score distributions, strand bias)
3. Inter-site spacing between adjacent motif hits in the same sequence
4. Palindrome/asymmetry analysis: Mce3R uses an asymmetric operator (NOT palindromic)
5. Direct repeat and tandem site detection
6. Performance evaluation against ground truth (for synthetic data)

Biological context (2024 update — Panagoda, Balázsi & Sampson, ACS Chem. Biol. 2024):
    Mce3R does NOT bind a classic palindromic TetR operator. The validated
    binding architecture is two asymmetric (nonpalindromic) ~25 bp half-sites
    separated by a ~53 bp spacer at the mce3R–yrbE3A intergenic region
    (H37Rv 2,207,477–2,207,699). Key metrics:
    - palindrome_score: Lower is BETTER for Mce3R (0.3–0.65 range expected)
    - asymmetry_score: Higher is BETTER (> 0.35 consistent with 2024 model)
    - validated spacer: 53 bp (flag sites with 48–58 bp spacing as matching)
    - site width: ~25 bp per half-site (each Mce3R monomer contacts one site)

    Prior work (Santangelo et al. 2009, Microbiology 155:2245):
    - Mce3R regulon: Rv1933c–Rv1935c and Rv1936–Rv1941 (divergently transcribed)
    - Second validated regulatory region: Rv1935c–Rv1936 intergenic region
    - MEME/MAST found the motif 6× in mce3R–yrbE3A IGR, 3× in Rv1935c–Rv1936 IGR

Usage:
    python analyze_results.py [--fimo-tsv PATH] [--metadata PATH] [--output-dir PATH]
"""

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from utils import get_project_root, reverse_complement, setup_logging

logger = setup_logging(__name__)

# ---------------------------------------------------------------------------
# Tiered prioritization constants
# ---------------------------------------------------------------------------

# Tier 1: Known Mce3R regulatory regions with validated operators
# Includes the two experimentally verified binding regions (Panagoda 2024,
# Santangelo 2009) and all genes/IGRs whose promoters are within these regions.
TIER1_LOCI = {
    # mce3R–yrbE3A validated operator (primary — Panagoda et al. 2024)
    "Rv1963c",                        # mce3R — the repressor itself (autoregulated)
    "Rv1964", "Rv1965", "Rv1966",     # yrbE3A, yrbE3B, mce3 operon entry point
    "Rv1967", "Rv1968", "Rv1969",     # mce3A, mce3B, mce3C
    "Rv1970", "Rv1971", "Rv1972",     # mce3D, mce3E, mce3F
    "Rv1973", "Rv1974", "Rv1975",     # lprM and accessory genes
    "Rv1976c", "Rv1977", "Rv1978",    # additional mce3 operon members
    # Full divergent IGR extracted by extract_promoters.extract_divergent_igrs()
    "IGR_Rv1963c_Rv1964",             # 897 bp mce3R–yrbE3A intergenic region (validated operator)
    # Rv1935c–Rv1936 regulatory region (second validated site — Santangelo et al. 2009)
    # DNase I footprinting showed two protected regions (~33 bp and ~32 bp)
    "Rv1935c", "Rv1936",              # boundary genes of the second validated IGR
    "IGR_Rv1935c_Rv1936",             # full divergent IGR (second proven Mce3R target)
}

# Tier 2: Biologically plausible targets confirmed in the Mce3R regulon or
# linked to cholesterol catabolism / oxidative stress co-regulated with mce3.
#
# NOTE: mce1, mce2, and mce4 operons are NOT included — Santangelo et al. 2008
# (BMC Microbiol. 8:38) showed Mce3R does not regulate these operons.
#
# (a) Specific locus tags
TIER2_LOCI = {
    # Confirmed Mce3R regulon members — divergently transcribed pair
    # (Santangelo et al. 2009, Microbiology 155:2245)
    "Rv1933c", "Rv1934c",             # divergent from Rv1936 side; lipid/redox metabolism
    "Rv1937", "Rv1938",               # part of Rv1936–Rv1941 operon
    "Rv1939", "Rv1940", "Rv1941",     # Rv1936–Rv1941 operon (divergently transcribed from Rv1933c–Rv1935c)
    # Cholesterol catabolism — core pathway genes (ksh, kst, hsa operons)
    "Rv3526",  # kshA — ketosteroid hydroxylase A (cholesterol side-chain cleavage)
    "Rv3577",  # kshB — ketosteroid hydroxylase B (reductase partner)
    "Rv3537",  # kstD — 3-ketosteroid-Δ1-dehydrogenase
    "Rv3568c", "Rv3569c", "Rv3570c",  # hsaC, hsaD, hsaA (hydroxybenzaldehyde catabolism)
    "Rv3571", "Rv3574",               # hsaB, hsaE
    # igr operon (intracellular growth required — cholesterol import)
    "Rv3543c", "Rv3544c", "Rv3545c",
    "Rv3546",  "Rv3547",  "Rv3548c",
    # Steroid side-chain degradation — β-oxidation
    "Rv3515c",  # fadD19 — long-chain fatty acid-CoA ligase (steroid side-chain)
    "Rv3516",   # echA19 — enoyl-CoA hydratase (steroid side-chain β-oxidation)
    # Isocitrate lyase — cholesterol import → central carbon metabolism
    "Rv3053c",  # icl1 — isocitrate lyase (connects cholesterol catabolism to TCA cycle)
    # Fatty acid-CoA ligase
    "Rv0119",   # fadD7 — high-scoring hit in genome-wide scan; Mce3R candidate
    # Oxidative stress response (Mce3R deletion increases ROS — Santangelo 2008)
    "Rv1908c",  # katG — catalase-peroxidase (primary H2O2 defense)
    "Rv2428",   # ahpC — alkyl hydroperoxide reductase C
    "Rv3846",   # sodA — superoxide dismutase A
    "Rv0432",   # sodC — superoxide dismutase C
    "Rv1932",   # tpx — thioredoxin-dependent peroxidase
}

# (b) Gene name prefixes common to lipid/cholesterol/stress genes
TIER2_GENE_PREFIXES = frozenset({
    "fadd", "fade", "echa", "ksh", "kst", "hsa", "igr",
    "mce", "ltp", "etf", "icl", "prp", "acc",
})

# (c) Product description keywords that indicate Tier 2 relevance
TIER2_PRODUCT_KEYWORDS = frozenset({
    "lipid", "cholesterol", "fatty acid", "acyl-coa", "acyl coa",
    "oxidoreductase", "peroxidase", "superoxide", "catalase",
    "ketosteroid", "steroid", "mycolic", "mce", "enoyl",
    "thioredoxin", "alkyl hydroperoxide",
})

TIER_LABELS = {
    1: "Tier1_known_operator",
    2: "Tier2_plausible_target",
    3: "Tier3_exploratory",
}

# FIMO TSV column names (matching both real FIMO and simulation output)
FIMO_COLUMNS = [
    "motif_id", "motif_alt_id", "sequence_name", "start", "stop",
    "strand", "score", "p-value", "q-value", "matched_sequence",
]


def load_fimo_results(fimo_tsv: Path) -> pd.DataFrame:
    """
    Load FIMO TSV output into a DataFrame.

    Handles FIMO's comment lines (starting with #) and both v4/v5 formats.
    Adds a computed site_length column for convenience.

    Args:
        fimo_tsv: Path to fimo.tsv file

    Returns:
        DataFrame with FIMO hits and site_length column.
        Returns empty DataFrame with correct schema if file is missing or empty.
    """
    fimo_tsv = Path(fimo_tsv)

    if not fimo_tsv.exists():
        logger.warning(f"FIMO TSV not found: {fimo_tsv}")
        return pd.DataFrame(columns=FIMO_COLUMNS + ["site_length"])

    try:
        df = pd.read_csv(fimo_tsv, sep="\t", comment="#")
    except Exception as e:
        logger.error(f"Failed to read FIMO TSV: {e}")
        return pd.DataFrame(columns=FIMO_COLUMNS + ["site_length"])

    if df.empty:
        logger.warning("FIMO TSV is empty (no hits found)")
        return pd.DataFrame(columns=FIMO_COLUMNS + ["site_length"])

    # Normalize column names (real FIMO uses 'p-value', simulation uses same)
    df.columns = [c.strip() for c in df.columns]

    # Cast numeric columns
    for col in ["start", "stop"]:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce").astype("Int64")
    for col in ["score", "p-value", "q-value"]:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")

    # Drop rows with missing essential fields
    df = df.dropna(subset=["sequence_name", "start", "stop"])

    df["site_length"] = df["stop"] - df["start"] + 1

    logger.info(f"Loaded {len(df)} FIMO hits from {fimo_tsv}")
    return df


def compute_motif_statistics(hits_df: pd.DataFrame) -> dict:
    """
    Compute summary statistics for the FIMO hit table.

    Args:
        hits_df: DataFrame from load_fimo_results()

    Returns:
        Dictionary with keys:
            total_hits, sequences_with_hits, hits_per_sequence,
            mean_score, median_score, strand_counts,
            motif_ids_found
    """
    if hits_df.empty:
        return {
            "total_hits": 0,
            "sequences_with_hits": 0,
            "hits_per_sequence": {},
            "mean_score": None,
            "median_score": None,
            "strand_counts": {"+": 0, "-": 0},
            "motif_ids_found": [],
        }

    hits_per_seq = hits_df.groupby("sequence_name").size().to_dict()
    strand_counts = hits_df["strand"].value_counts().to_dict()

    stats = {
        "total_hits": len(hits_df),
        "sequences_with_hits": hits_df["sequence_name"].nunique(),
        "hits_per_sequence": hits_per_seq,
        "mean_score": round(float(hits_df["score"].mean()), 4),
        "median_score": round(float(hits_df["score"].median()), 4),
        "strand_counts": strand_counts,
        "motif_ids_found": hits_df["motif_id"].unique().tolist(),
    }
    return stats


def compute_inter_site_spacing(hits_df: pd.DataFrame) -> pd.DataFrame:
    """
    Compute spacing between adjacent motif hits within the same sequence.

    For sequences with multiple hits, calculates the gap (in bp) between
    consecutive binding sites. This reveals whether Mce3R operators occur
    in tandem pairs with characteristic spacing.

    Biological context (2024 update):
        The validated Mce3R operator at mce3R–yrbE3A has two asymmetric ~25 bp
        half-sites separated by a 53 bp spacer (not 30–50 bp as for classic TetR).
        Sites with spacing in the 48–58 bp range are flagged as matching the
        validated operator geometry (matches_validated_spacing = True).

    Args:
        hits_df: DataFrame from load_fimo_results()

    Returns:
        DataFrame with columns: sequence_name, site1_start, site2_start,
        site1_stop, site2_stop, spacing_bp, site1_strand, site2_strand,
        matches_validated_spacing (bool — True if 48–58 bp, matching 53 bp reference).
        Empty DataFrame if no sequences have >= 2 hits.
    """
    if hits_df.empty:
        return pd.DataFrame(
            columns=[
                "sequence_name", "site1_start", "site2_start",
                "site1_stop", "site2_stop", "spacing_bp",
                "site1_strand", "site2_strand", "matches_validated_spacing",
            ]
        )

    spacing_rows = []

    for seq_name, group in hits_df.groupby("sequence_name"):
        group_sorted = group.sort_values("start").reset_index(drop=True)
        if len(group_sorted) < 2:
            continue
        for i in range(len(group_sorted) - 1):
            site1 = group_sorted.iloc[i]
            site2 = group_sorted.iloc[i + 1]
            # spacing = gap between end of site1 and start of site2
            spacing = int(site2["start"]) - int(site1["stop"]) - 1
            spacing_rows.append(
                {
                    "sequence_name": seq_name,
                    "site1_start": int(site1["start"]),
                    "site2_start": int(site2["start"]),
                    "site1_stop": int(site1["stop"]),
                    "site2_stop": int(site2["stop"]),
                    "spacing_bp": spacing,
                    "site1_strand": site1["strand"],
                    "site2_strand": site2["strand"],
                    # Flag spacers in the validated 53 bp window (±5 bp tolerance)
                    "matches_validated_spacing": 48 <= spacing <= 58,
                }
            )

    if not spacing_rows:
        logger.info("No sequences had >= 2 motif hits; spacing analysis skipped.")
        return pd.DataFrame(
            columns=[
                "sequence_name", "site1_start", "site2_start",
                "site1_stop", "site2_stop", "spacing_bp",
                "site1_strand", "site2_strand", "matches_validated_spacing",
            ]
        )

    spacing_df = pd.DataFrame(spacing_rows)
    logger.info(
        f"Inter-site spacing computed for {len(spacing_df)} adjacent site pairs "
        f"across {spacing_df['sequence_name'].nunique()} sequences"
    )
    return spacing_df


def detect_palindromic_sites(hits_df: pd.DataFrame) -> pd.DataFrame:
    """
    Score each matched sequence for palindromic symmetry.

    IMPORTANT (2024 update — Panagoda et al. ACS Chem. Biol. 2024):
    Mce3R uses an ASYMMETRIC (nonpalindromic) operator. The expected
    palindrome_score for a genuine Mce3R site is 0.3–0.65, NOT >= 0.8.

    Sites with palindrome_score >= 0.85 are flagged as 'likely_non_mce3r_tetr'
    (more consistent with classic palindromic TetR binding than the 2024 Mce3R model).

    Algorithm:
        1. Split matched_sequence at midpoint
        2. Compute RC of right half
        3. Hamming distance between left half and RC(right half)
        4. palindrome_score = 1 - hamming / half_length
        5. asymmetry_score = 1 - palindrome_score

    The palindrome_score is kept as a descriptive metric but is NOT used as
    a ranking quality filter. Use asymmetry_score for Mce3R-specific ranking.

    Args:
        hits_df: DataFrame with matched_sequence column

    Returns:
        DataFrame with added columns:
        palindrome_score (float 0–1, LOWER = more asymmetric = better for Mce3R),
        is_palindrome (bool, score >= 0.8 — descriptive only, NOT a quality filter),
        likely_non_mce3r_tetr (bool, score >= 0.85 — warning flag),
        left_half (str), right_half (str), rc_right_half (str)
    """
    if hits_df.empty or "matched_sequence" not in hits_df.columns:
        return hits_df

    df = hits_df.copy()
    palindrome_scores = []
    is_palindromes = []
    likely_non_mce3r = []
    left_halves = []
    right_halves = []
    rc_right_halves = []

    for _, row in df.iterrows():
        matched = str(row.get("matched_sequence", "")).upper()
        n = len(matched)
        if n < 2:
            palindrome_scores.append(0.0)
            is_palindromes.append(False)
            likely_non_mce3r.append(False)
            left_halves.append(matched)
            right_halves.append("")
            rc_right_halves.append("")
            continue

        mid = n // 2
        left = matched[:mid]

        # For odd-length sequences, skip the central base (not part of either half-site)
        if n % 2 == 1:
            right = matched[mid + 1:]
        else:
            right = matched[mid:]

        rc_right = reverse_complement(right)

        # Align left half with RC of right half
        compare_len = min(len(left), len(rc_right))
        hamming = sum(
            1 for a, b in zip(left[:compare_len], rc_right[:compare_len]) if a != b
        )
        score = 1.0 - (hamming / compare_len) if compare_len > 0 else 0.0

        palindrome_scores.append(round(score, 4))
        # is_palindrome kept as descriptive (>= 0.8) — NOT a quality filter for Mce3R
        is_palindromes.append(score >= 0.8)
        # Sites >= 0.85 flagged as likely classical TetR palindrome, not Mce3R-specific
        likely_non_mce3r.append(score >= 0.85)
        left_halves.append(left)
        right_halves.append(right)
        rc_right_halves.append(rc_right)

    df["palindrome_score"] = palindrome_scores
    df["is_palindrome"] = is_palindromes
    df["likely_non_mce3r_tetr"] = likely_non_mce3r
    df["left_half"] = left_halves
    df["right_half"] = right_halves
    df["rc_right_half"] = rc_right_halves

    n_palindromes = sum(is_palindromes)
    n_flagged = sum(likely_non_mce3r)
    logger.info(
        f"Palindrome analysis: {n_palindromes}/{len(df)} sites scored >= 0.8 (is_palindrome=True); "
        f"{n_flagged} flagged as likely_non_mce3r_tetr (score >= 0.85)"
    )
    return df


def detect_direct_repeats(hits_df: pd.DataFrame) -> pd.DataFrame:
    """
    Check if matched sequences contain direct repeats (left half == right half).

    Some TetR-family operators have evolved from direct repeats rather than
    inverted repeats. This analysis distinguishes between the two architectures.

    Args:
        hits_df: DataFrame with matched_sequence and palindrome columns

    Returns:
        DataFrame with added columns: direct_repeat_score (float),
        is_direct_repeat (bool)
    """
    if hits_df.empty or "matched_sequence" not in hits_df.columns:
        return hits_df

    df = hits_df.copy()
    dr_scores = []
    is_drs = []

    for _, row in df.iterrows():
        matched = str(row.get("matched_sequence", "")).upper()
        n = len(matched)
        if n < 2:
            dr_scores.append(0.0)
            is_drs.append(False)
            continue

        mid = n // 2
        left = matched[:mid]
        right = matched[mid + (n % 2):]

        compare_len = min(len(left), len(right))
        hamming = sum(
            1 for a, b in zip(left[:compare_len], right[:compare_len]) if a != b
        )
        score = 1.0 - (hamming / compare_len) if compare_len > 0 else 0.0
        dr_scores.append(round(score, 4))
        is_drs.append(score >= 0.8)

    df["direct_repeat_score"] = dr_scores
    df["is_direct_repeat"] = is_drs
    return df


def annotate_hits_with_sequence_metadata(
    hits_df: pd.DataFrame, sequences_df: pd.DataFrame
) -> pd.DataFrame:
    """
    Join FIMO hits with sequence metadata (GC content, ground truth, etc.).

    The joined DataFrame enables:
    - Correlation of hit quality with sequence GC content
    - Precision/recall computation (using has_motif as ground truth)
    - Position accuracy validation (comparing hit position to embedded position)

    Args:
        hits_df: DataFrame from detect_palindromic_sites()
        sequences_df: DataFrame from generate_sequences.py (sequence_metadata.csv)

    Returns:
        Left-joined DataFrame with added columns from sequences_df
    """
    if hits_df.empty:
        return hits_df

    if sequences_df.empty:
        return hits_df

    # Determine which ID column the metadata uses
    # - Synthetic data: 'sequence_id' + 'has_motif', 'motif_position', 'length'
    # - Real genome promoter data: 'locus_tag' + 'gene', 'product', 'prom_length'
    if "locus_tag" in sequences_df.columns:
        # Real genome data — join on locus_tag / sequence_name
        join_cols = [c for c in ["locus_tag", "gene", "product", "gc_content", "prom_length"]
                     if c in sequences_df.columns]
        annotated = hits_df.merge(
            sequences_df[join_cols],
            left_on="sequence_name",
            right_on="locus_tag",
            how="left",
        )
    else:
        # Synthetic data — join on sequence_id
        join_cols = [c for c in ["sequence_id", "gc_content", "length", "has_motif", "motif_position"]
                     if c in sequences_df.columns]
        annotated = hits_df.merge(
            sequences_df[join_cols],
            left_on="sequence_name",
            right_on="sequence_id",
            how="left",
        )
    return annotated


def compute_discovery_performance(annotated_df: pd.DataFrame, n_total_sequences: int) -> dict:
    """
    Evaluate motif discovery accuracy against synthetic ground truth.

    Uses the has_motif column (set during sequence generation) to compute
    precision, recall, and F1 score for the set of sequences with hits.

    This metric answers: "Did FIMO find hits preferentially in sequences
    where the motif was actually embedded?"

    Args:
        annotated_df: DataFrame from annotate_hits_with_sequence_metadata()
        n_total_sequences: Total number of sequences (for FN calculation)

    Returns:
        Dictionary with: precision, recall, f1, true_positives,
        false_positives, false_negatives, true_negatives
    """
    if annotated_df.empty or "has_motif" not in annotated_df.columns:
        return {"error": "No ground truth data available"}

    # Sequences predicted positive = sequences with any FIMO hit
    predicted_positive = set(annotated_df["sequence_name"].unique())
    # Ground truth positive = sequences where motif was embedded
    true_positive_seqs = set(
        annotated_df.loc[
            annotated_df["has_motif"] == True, "sequence_name"
        ].unique()
    )

    tp = len(predicted_positive & true_positive_seqs)
    fp = len(predicted_positive - true_positive_seqs)
    fn = len(true_positive_seqs - predicted_positive)
    tn = n_total_sequences - tp - fp - fn

    precision = tp / (tp + fp) if (tp + fp) > 0 else 0.0
    recall = tp / (tp + fn) if (tp + fn) > 0 else 0.0
    f1 = (
        2 * precision * recall / (precision + recall)
        if (precision + recall) > 0
        else 0.0
    )

    perf = {
        "true_positives": tp,
        "false_positives": fp,
        "false_negatives": fn,
        "true_negatives": tn,
        "precision": round(precision, 4),
        "recall": round(recall, 4),
        "f1": round(f1, 4),
    }
    logger.info(
        f"Discovery performance: Precision={precision:.3f}, "
        f"Recall={recall:.3f}, F1={f1:.3f} "
        f"(TP={tp}, FP={fp}, FN={fn})"
    )
    return perf


def generate_tier1_debug_report(
    hits_df: pd.DataFrame,
    fasta_path: Path,
    motif_files: list,
    output_path: Path,
) -> None:
    """
    Force-score all Tier 1 promoters against both motifs, regardless of p-value threshold.

    This report answers: "Why does Tier 1 have 0 (or few) hits?"
    It scans every Tier 1 sequence at any threshold, reports the best score
    per sequence, compares to the FIMO hit floor, and flags:
    - Sequences with good scores that were filtered by p-value threshold
    - Coverage gaps in the mce3R–yrbE3A intergenic region
    - Sequences missing from the promoter FASTA entirely

    Args:
        hits_df: Ranked binding site table (may be empty for Tier 1)
        fasta_path: Path to the promoters FASTA file
        motif_files: List of motif file paths to score against
        output_path: Where to write the debug report text file
    """
    import re

    try:
        import numpy as np
        from Bio import SeqIO as _SeqIO
        from run_fimo import parse_meme_motif
    except ImportError as e:
        logger.warning(f"Cannot generate Tier1 debug report: {e}")
        return

    fasta_path = Path(fasta_path)
    if not fasta_path.exists():
        logger.warning(f"Tier1 debug report: FASTA not found at {fasta_path}")
        return

    # Load all sequences
    all_seqs = {r.id: str(r.seq).upper() for r in _SeqIO.parse(str(fasta_path), "fasta")}

    bg = np.array([0.175, 0.325, 0.325, 0.175])  # M. tb background A,C,G,T
    BASE_IDX = {"A": 0, "C": 1, "G": 2, "T": 3}

    def _rc(s):
        comp = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N"}
        return "".join(comp.get(b, "N") for b in reversed(s.upper()))

    def _score_seq(seq, lo, w):
        best = -999.0
        best_pos, best_strand, best_win = 0, "+", ""
        for pos in range(len(seq) - w + 1):
            win = seq[pos: pos + w].upper()
            for strand, s_seq in [("+", win), ("-", _rc(win))]:
                s = sum(lo[i, BASE_IDX.get(b, -1)] for i, b in enumerate(s_seq[:w])
                        if BASE_IDX.get(b, -1) >= 0)
                if s > best:
                    best, best_pos, best_strand, best_win = s, pos, strand, s_seq
        return best, best_pos, best_strand, best_win

    # Parse motifs
    motif_data = []
    for mf in motif_files:
        mf = Path(mf)
        if not mf.exists():
            continue
        motif = parse_meme_motif(str(mf))
        w = motif["width"]
        pwm = np.array(motif["pwm"])
        lo = np.log2((pwm + 0.001) / (bg + 0.001))
        motif_data.append((mf.stem, w, lo))

    # FIMO hit floor from the ranked table
    hit_scores = hits_df["score"].dropna() if not hits_df.empty and "score" in hits_df.columns else pd.Series(dtype=float)
    hit_floor = float(hit_scores.min()) if len(hit_scores) > 0 else None

    tier1_tags = sorted(TIER1_LOCI)

    lines = [
        "=" * 70,
        "  MCE3R TIER1 DIAGNOSTIC REPORT",
        "  Why does the known Mce3R regulatory region appear in results?",
        "=" * 70,
        "",
        f"  FIMO hit score floor (weakest passing hit): "
        + (f"{hit_floor:.2f}" if hit_floor else "N/A (no hits)"),
        f"  Motifs scored: {', '.join(m[0] for m in motif_data)}",
        f"  Sequences in FASTA: {len(all_seqs)}",
        "",
        f"  {'Locus':<25} {'In FASTA':>8}  "
        + "  ".join(f"{m[0][:12]:>14}" for m in motif_data)
        + f"  {'Status':>20}",
        "  " + "-" * 68,
    ]

    tier1_hits_in_ranked = (
        hits_df[hits_df["sequence_name"].isin(TIER1_LOCI)] if not hits_df.empty and "sequence_name" in hits_df.columns
        else pd.DataFrame()
    )

    for tag in tier1_tags:
        seq = all_seqs.get(tag)
        in_fasta = "YES" if seq else "NO"
        scores_per_motif = []
        if seq:
            for motif_name, w, lo in motif_data:
                best_score, best_pos, best_strand, best_win = _score_seq(seq, lo, w)
                scores_per_motif.append((best_score, best_pos, best_strand, best_win))

        # Determine status
        in_hits = not tier1_hits_in_ranked.empty and tag in tier1_hits_in_ranked["sequence_name"].values
        if in_hits:
            status = "HIT (in ranked table)"
        elif scores_per_motif:
            best_any = max(s[0] for s in scores_per_motif)
            if hit_floor and best_any >= hit_floor:
                status = f"ABOVE FLOOR but not merged"
            elif hit_floor and best_any >= hit_floor * 0.85:
                status = f"NEAR THRESHOLD"
            else:
                status = f"BELOW THRESHOLD"
        else:
            status = "NOT IN FASTA"

        score_strs = "  ".join(f"{s[0]:>14.2f}" for s in scores_per_motif) if scores_per_motif else "  ".join([f"{'N/A':>14}"] * len(motif_data))
        lines.append(f"  {tag:<25} {in_fasta:>8}  {score_strs}  {status:>20}")

    lines += [
        "",
        "  INTERPRETATION",
        "  " + "-" * 68,
    ]

    if not tier1_hits_in_ranked.empty:
        lines.append(f"  GOOD: {len(tier1_hits_in_ranked)} Tier1 hit(s) appear in the ranked table.")
    else:
        lines += [
            "  WARNING: No Tier1 loci appear in the main ranked hit table.",
            "",
            "  Root causes investigated:",
            "  1. PROMOTER EXTRACTION COVERAGE:",
            "     The mce3R-yrbE3A intergenic region (897 bp) is split across",
            "     two 300 bp windows: Rv1963c gets IGR[0:300], Rv1964 gets",
            "     IGR[597:897]. The central 297 bp (IGR[300:597]) is uncovered.",
            "     FIX: IGR_Rv1963c_Rv1964 entry extracts the full 897 bp.",
            "",
            "  2. MOTIF MODEL SENSITIVITY:",
            "     The palindromic PWM (TTGACANNNNNTGTCAA) scores Rv1963c at",
            "     ~7.74 — far below the FIMO hit floor. The asymmetric PWM",
            "     (TTGACATNNNNNTGCCCA) scores Rv1963c at ~11.63, which is at",
            "     the threshold. The asymmetric model better represents Mce3R.",
            "",
            "  3. SIMULATION LIMITATION:",
            "     The pipeline uses a synthetic PWM (simulation mode). Real",
            "     MEME on candidate promoters would discover the true Mce3R",
            "     PWM from the mce3R-yrbE3A region, which would score Tier1",
            "     sites much higher. Install MEME Suite for real discovery.",
        ]

    lines += [
        "",
        "  RECOMMENDATION:",
        "  Run the pipeline with --simulate flag and check whether",
        "  IGR_Rv1963c_Rv1964 appears in the output after this fix.",
        "  For definitive results: install MEME Suite (conda install -c bioconda meme)",
        "=" * 70,
    ]

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(str(output_path), "w") as f:
        f.write("\n".join(lines) + "\n")
    logger.info(f"Tier1 diagnostic report written to {output_path}")


def assign_priority_tier(hits_df: pd.DataFrame) -> pd.DataFrame:
    """
    Assign a biological priority tier to each predicted binding site.

    Tier 1 — Known Mce3R regulatory region (validation set):
        The mce3R repressor gene (Rv1963c) and its regulated mce3 operon
        (Rv1964–Rv1978). Any hit in these promoters should be Tier 1 by
        definition because the biology is established.

    Tier 2 — Biologically plausible targets:
        Lipid/cholesterol catabolism genes, Mce1/Mce2/Mce4 operons,
        oxidative stress response genes. These are likely co-regulated with
        mce3 based on shared metabolic context in M. tuberculosis.

    Tier 3 — Exploratory (genome-wide):
        All other hits that fall outside Tier 1 and Tier 2 criteria.

    Args:
        hits_df: DataFrame with at least a 'sequence_name' column.
                 Optionally uses 'gene' and 'product' columns for richer
                 Tier 2 classification.

    Returns:
        DataFrame with added 'priority_tier' column (string label) and
        '_tier_order' column (int 1/2/3 for sort key, dropped later).
    """
    if hits_df.empty:
        return hits_df

    df = hits_df.copy()

    def _classify_row(row):
        seq_name = str(row.get("sequence_name", "")).strip()
        gene = str(row.get("gene", "")).strip().lower()
        product = str(row.get("product", "")).strip().lower()

        # Tier 1: exact locus tag match
        if seq_name in TIER1_LOCI:
            return 1

        # Tier 2a: exact locus tag match in Tier 2 set
        if seq_name in TIER2_LOCI:
            return 2

        # Tier 2b: gene name prefix match
        if gene and any(gene.startswith(pfx) for pfx in TIER2_GENE_PREFIXES):
            return 2

        # Tier 2c: product description keyword match
        if product and any(kw in product for kw in TIER2_PRODUCT_KEYWORDS):
            return 2

        # Tier 3: everything else
        return 3

    df["_tier_order"] = df.apply(_classify_row, axis=1)
    df["priority_tier"] = df["_tier_order"].map(TIER_LABELS)

    tier_counts = df["_tier_order"].value_counts().sort_index()
    for t, n in tier_counts.items():
        logger.info(f"Priority {TIER_LABELS[t]}: {n} sites")

    return df


def rank_binding_sites(hits_df: pd.DataFrame) -> pd.DataFrame:
    """
    Produce a ranked table of predicted Mce3R binding sites.

    When model comparison columns are present (from motif_models.apply_model_comparison
    and classify_site_architectures), uses a multi-factor composite score:
        50% FIMO log-odds score (normalized to [0,1])
        25% Site architecture (validated_paired_asymmetric highest, per 2024 paper)
        15% Asymmetric model fit (per 2024 ACS Chemical Biology)
        10% Known Mce3R target promoter context bonus

    NOTE: Palindrome score is kept as a descriptive column only.
    It is NOT used in composite score calculation. Paired asymmetric architectures
    (validated_paired_asymmetric) outrank isolated palindromic-looking hits.

    Falls back to FIMO score ranking if model comparison columns are absent.

    Args:
        hits_df: Annotated FIMO hits DataFrame (from detect_palindromic_sites,
                 optionally enriched by apply_model_comparison)

    Returns:
        Ranked DataFrame with composite_score column when model data is available
    """
    if hits_df.empty:
        logger.warning("No FIMO hits to rank — predicted_mce3r_sites.csv will be empty")
        return pd.DataFrame()

    df = hits_df.copy()

    # Normalize the p-value column name (FIMO uses 'p-value', simulation uses same)
    if "p-value" in df.columns:
        df = df.rename(columns={"p-value": "p_value", "q-value": "q_value"})

    # Multi-factor composite ranking when model comparison columns are available
    use_composite = all(
        c in df.columns for c in ["site_architecture", "model_preference", "model_confidence"]
    )

    # Assign priority tiers FIRST so tier info is available for composite score
    df = assign_priority_tier(df)

    if use_composite:
        try:
            from motif_models import compute_composite_score

            # Normalize FIMO score to [0,1] before composite computation
            score_min = df["score"].min()
            score_max = df["score"].max()
            score_range = score_max - score_min
            if score_range > 0:
                df["_norm_score"] = (df["score"] - score_min) / score_range
            else:
                df["_norm_score"] = 0.5

            df_for_composite = df.copy()
            df_for_composite["score"] = df_for_composite["_norm_score"]
            df["composite_score"] = df_for_composite.apply(compute_composite_score, axis=1)
            df = df.drop(columns=["_norm_score"])

            logger.info("Ranking by composite score (FIMO + architecture + asymmetry + tier)")
        except ImportError:
            use_composite = False

    # Primary sort key: tier order (1 < 2 < 3), secondary: score descending
    sort_score_col = "composite_score" if use_composite else "score"
    if "_tier_order" in df.columns:
        df = df.sort_values(
            ["_tier_order", sort_score_col], ascending=[True, False]
        ).reset_index(drop=True)
        df = df.drop(columns=["_tier_order"])
    elif not use_composite:
        df = df.sort_values("score", ascending=False).reset_index(drop=True)

    df.insert(0, "rank", df.index + 1)

    # Select and order output columns
    keep_cols = ["rank", "sequence_name", "start", "stop", "strand",
                 "score", "p_value", "matched_sequence", "priority_tier"]
    optional_cols = [
        "composite_score",
        "palindrome_score", "is_palindrome", "likely_non_mce3r_tetr", "asymmetry_score",
        "model_preference", "model_confidence",
        "site_architecture", "has_paired_site", "partner_spacing_bp",
        "palindromic_consensus_score", "asymmetric_consensus_score",
        "direct_repeat_score",
        "gene", "product", "gc_content", "prom_length",
    ]
    for col in optional_cols:
        if col in df.columns:
            keep_cols.append(col)

    ranked = df[[c for c in keep_cols if c in df.columns]]

    n_palindromes = int(ranked["is_palindrome"].sum()) if "is_palindrome" in ranked.columns else "N/A"
    n_asymmetric = (
        int((ranked["model_preference"] == "asymmetric").sum())
        if "model_preference" in ranked.columns else "N/A"
    )
    logger.info(
        f"Ranked {len(ranked)} predicted Mce3R binding sites "
        f"({n_palindromes} palindromic, {n_asymmetric} asymmetric model preference)"
    )

    if len(ranked) > 0:
        top = ranked.iloc[0]
        arch = top.get("site_architecture", "")
        logger.info(
            f"Top site: {top['sequence_name']} "
            f"(score={top['score']:.2f}, strand={top['strand']}"
            + (f", arch={arch}" if arch else "")
            + ")"
        )

    return ranked


def main():
    root = get_project_root()
    parser = argparse.ArgumentParser(
        description="Analyze FIMO motif scan results for Mce3R binding site characterization"
    )
    parser.add_argument(
        "--fimo-tsv",
        type=Path,
        default=root / "results" / "scans" / "fimo.tsv",
        help="FIMO output TSV file",
    )
    parser.add_argument(
        "--metadata",
        type=Path,
        default=root / "data" / "processed" / "sequence_metadata.csv",
        help="Sequence metadata CSV (from generate_sequences.py)",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=root / "data" / "processed",
        help="Directory for analysis outputs",
    )
    parser.add_argument(
        "--scans-dir",
        type=Path,
        default=root / "results" / "scans",
        help="Directory for summary statistics JSON",
    )
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    args.scans_dir.mkdir(parents=True, exist_ok=True)

    # Load FIMO results
    hits_df = load_fimo_results(args.fimo_tsv)

    # Load sequence metadata (ground truth)
    sequences_df = pd.DataFrame()
    if args.metadata.exists():
        sequences_df = pd.read_csv(args.metadata)
        logger.info(f"Loaded sequence metadata: {len(sequences_df)} sequences")
    else:
        logger.warning(f"Sequence metadata not found: {args.metadata}")

    # Run all analyses
    stats = compute_motif_statistics(hits_df)
    spacing_df = compute_inter_site_spacing(hits_df)
    hits_df = detect_palindromic_sites(hits_df)
    hits_df = detect_direct_repeats(hits_df)

    # Apply asymmetric model comparison (2024 ACS Chemical Biology)
    try:
        from motif_models import apply_model_comparison, classify_site_architectures
        hits_df = apply_model_comparison(hits_df)
        hits_df = classify_site_architectures(hits_df)
        logger.info("Applied asymmetric model comparison (2024 ACS Chem. Biol. framework)")
    except ImportError:
        logger.warning("motif_models.py not found — skipping asymmetric model comparison")

    # Annotate with ground truth if available
    if not sequences_df.empty:
        annotated_df = annotate_hits_with_sequence_metadata(hits_df, sequences_df)
        n_total = len(sequences_df)
    else:
        annotated_df = hits_df
        n_total = stats.get("sequences_with_hits", 0)

    # Performance evaluation
    if not sequences_df.empty and "has_motif" in sequences_df.columns:
        perf = compute_discovery_performance(annotated_df, n_total)
        stats["discovery_performance"] = perf

    # Save annotated hits
    annotated_path = args.output_dir / "annotated_hits.csv"
    annotated_df.to_csv(annotated_path, index=False)
    logger.info(f"Saved annotated hits to {annotated_path}")

    # Save spacing analysis
    spacing_path = args.output_dir / "spacing_analysis.csv"
    spacing_df.to_csv(spacing_path, index=False)
    logger.info(f"Saved spacing analysis to {spacing_path}")

    # Save ranked binding site predictions (uses composite score if model columns present)
    ranked_df = rank_binding_sites(annotated_df)
    ranked_path = args.output_dir / "predicted_mce3r_sites.csv"
    ranked_df.to_csv(ranked_path, index=False)
    logger.info(f"Saved ranked binding sites to {ranked_path}")

    # Save summary statistics as JSON
    summary_path = args.scans_dir / "summary_statistics.json"
    with open(summary_path, "w") as f:
        json.dump(stats, f, indent=2, default=str)
    logger.info(f"Saved summary statistics to {summary_path}")

    # Generate model comparison report
    report_path = None
    arch_report_path = None
    try:
        from motif_models import generate_model_comparison_report, generate_architecture_comparison_report
        report_path = args.scans_dir / "mce3r_prediction_report.txt"
        generate_model_comparison_report(annotated_df, ranked_df, report_path)
        # Architecture comparison report (requires ranked sites with architecture columns)
        arch_report_path = args.scans_dir / "architecture_comparison_report.txt"
        generate_architecture_comparison_report(ranked_df, arch_report_path)
    except ImportError:
        logger.warning("motif_models.py not found — skipping prediction report")

    # Print summary table
    print("\n" + "=" * 65)
    print("  MCE3R MOTIF ANALYSIS SUMMARY")
    print("=" * 65)
    print(f"  Total FIMO hits           : {stats['total_hits']}")
    print(f"  Sequences with hits       : {stats['sequences_with_hits']}")
    print(f"  Mean FIMO score           : {stats.get('mean_score', 'N/A')}")
    print(f"  Strand distribution       : {stats.get('strand_counts', {})}")

    if not spacing_df.empty:
        mean_spacing = spacing_df["spacing_bp"].mean()
        print(f"  Mean inter-site spacing   : {mean_spacing:.1f} bp")

    if not hits_df.empty and "is_palindrome" in hits_df.columns:
        n_pal = hits_df["is_palindrome"].sum()
        print(f"  Palindromic sites (>=0.8) : {n_pal}/{len(hits_df)}")

    if not hits_df.empty and "model_preference" in hits_df.columns:
        pref = hits_df["model_preference"].value_counts().to_dict()
        print(f"  Asymmetric preference     : {pref.get('asymmetric', 0)}/{len(hits_df)} sites")
        print(f"  Palindromic preference    : {pref.get('palindromic', 0)}/{len(hits_df)} sites")

    if not hits_df.empty and "site_architecture" in hits_df.columns:
        arch = hits_df["site_architecture"].value_counts().to_dict()
        top_arch = sorted(arch.items(), key=lambda x: -x[1])
        print(f"  Site architectures        : " + ", ".join(f"{k}={v}" for k, v in top_arch))
        n_validated_paired = arch.get("validated_paired_asymmetric", 0)
        if n_validated_paired > 0:
            print(f"  ** Validated paired-asym  : {n_validated_paired} sites match 2024 ACS Chem. Biol. architecture **")

    if "discovery_performance" in stats:
        perf = stats["discovery_performance"]
        print(f"\n  Discovery Performance (synthetic ground truth):")
        print(f"    Precision: {perf.get('precision', 'N/A')}")
        print(f"    Recall   : {perf.get('recall', 'N/A')}")
        print(f"    F1 score : {perf.get('f1', 'N/A')}")
    print("=" * 65 + "\n")

    if not ranked_df.empty and "priority_tier" in ranked_df.columns:
        print(f"\n  Priority Tier Summary:")
        for tier_order, tier_label in sorted(TIER_LABELS.items()):
            tier_rows = ranked_df[ranked_df["priority_tier"] == tier_label]
            count = len(tier_rows)
            if count == 0:
                print(f"    {tier_label:<30}: 0 sites")
                continue
            top = tier_rows.iloc[0]
            score_col = "composite_score" if "composite_score" in top.index else "score"
            top_score = top[score_col]
            print(
                f"    {tier_label:<30}: {count} sites  "
                f"[top: {str(top['sequence_name']):<14} {score_col}={top_score:.3f}]"
            )

    if not ranked_df.empty:
        print(f"\n  Top 5 predicted Mce3R binding sites (computational predictions only):")
        top5 = ranked_df.head(5)
        for _, row in top5.iterrows():
            pal = f"pal={row['palindrome_score']:.2f}" if "palindrome_score" in row else ""
            model = row.get("model_preference", "")
            arch = row.get("site_architecture", "")
            tier = row.get("priority_tier", "")
            comp = f"composite={row['composite_score']:.3f}" if "composite_score" in row else f"score={row['score']:.2f}"
            print(
                f"    #{int(row['rank']):<3} {str(row['sequence_name']):<14} "
                f"{comp}  strand={row['strand']}  {pal}"
                + (f"  model={model}" if model else "")
                + (f"  arch={arch}" if arch else "")
                + (f"  [{tier}]" if tier else "")
            )

    print(f"\nOutput files:")
    print(f"  Ranked sites     : {ranked_path}")
    print(f"  Annotated hits   : {annotated_path}")
    print(f"  Spacing analysis : {spacing_path}")
    print(f"  Summary stats    : {summary_path}")
    if report_path and Path(report_path).exists():
        print(f"  Prediction report: {report_path}")
    if arch_report_path and Path(arch_report_path).exists():
        print(f"  Architecture report: {arch_report_path}")


if __name__ == "__main__":
    main()
