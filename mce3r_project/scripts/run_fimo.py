"""
run_fimo.py — FIMO motif scanning wrapper with Python log-odds fallback.

FIMO (Find Individual Motif Occurrences) scans DNA sequences for matches to
a position weight matrix (PWM). For Mce3R analysis, it identifies all candidate
binding sites in the M. tuberculosis promoter regions.

FIMO command used:
    fimo --oc results/scans/ results/motifs/meme.txt data/raw/sequences.fasta

Output (fimo.tsv) columns:
    motif_id, motif_alt_id, sequence_name, start, stop, strand,
    score, p-value, q-value, matched_sequence

Simulation mode (when FIMO is not installed):
    Implements a Python-based log-odds scanner that:
    1. Parses the PWM from meme.txt
    2. Converts to log-odds scores using M. tb background frequencies
    3. Slides a window across both strands of each sequence
    4. Computes empirical p-values from background score distribution
    5. Writes output in the exact FIMO TSV format

Biological context:
    Log-odds scoring compares the probability of observing a sequence
    under the motif model vs. the background model. A positive score
    indicates the sequence matches the motif better than expected by chance.
    This is identical to what FIMO computes internally.

Usage:
    python run_fimo.py --motif-file PATH --fasta PATH --output-dir PATH
"""

import argparse
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import SeqIO

sys.path.insert(0, str(Path(__file__).resolve().parent))
from utils import check_tool_available, get_project_root, reverse_complement, setup_logging
from run_meme import parse_meme_motif

logger = setup_logging(__name__)

# M. tuberculosis background base frequencies (GC-rich genome)
MTB_BACKGROUND = {"A": 0.175, "C": 0.325, "G": 0.325, "T": 0.175}

# FIMO TSV column headers (matching real FIMO output format)
FIMO_COLUMNS = [
    "motif_id",
    "motif_alt_id",
    "sequence_name",
    "start",
    "stop",
    "strand",
    "score",
    "p-value",
    "q-value",
    "matched_sequence",
]


def pwm_to_log_odds(pwm: np.ndarray, background: dict = None) -> np.ndarray:
    """
    Convert a PWM (probabilities) to a log-odds scoring matrix.

    log-odds[i][b] = log2(PWM[i][b] / background[b])

    A positive log-odds score at position i for base b means that base
    is more likely under the motif model than the background model.

    Args:
        pwm: Shape (width, 4) array of probabilities [A, C, G, T]
        background: Base frequency dict (defaults to M. tb frequencies)

    Returns:
        Shape (width, 4) log-odds matrix
    """
    if background is None:
        background = MTB_BACKGROUND

    bg_freq = np.array([background["A"], background["C"], background["G"], background["T"]])
    # Avoid log(0) by clipping very small probabilities
    pwm_clipped = np.clip(pwm, 1e-9, 1.0)
    log_odds = np.log2(pwm_clipped / bg_freq[np.newaxis, :])
    return log_odds


def score_window(window: str, log_odds: np.ndarray) -> float:
    """
    Score a DNA window using a log-odds matrix.

    Returns the sum of log-odds scores for each position in the window.
    Unknown bases (N) contribute 0 (no information).

    Args:
        window: DNA string of length == log_odds.shape[0]
        log_odds: Shape (width, 4) log-odds matrix

    Returns:
        Float log-odds score
    """
    base_to_idx = {"A": 0, "C": 1, "G": 2, "T": 3}
    score = 0.0
    for i, base in enumerate(window.upper()):
        idx = base_to_idx.get(base)
        if idx is not None:
            score += log_odds[i, idx]
    return score


def score_sequence_with_pwm(
    sequence: str,
    log_odds: np.ndarray,
    motif_width: int,
) -> list:
    """
    Slide the log-odds matrix across a sequence and score every window.

    Scans both the forward (+) and reverse complement (-) strands,
    as TetR operators can occur on either strand.

    Args:
        sequence: DNA string (upstream)
        log_odds: Shape (width, 4) log-odds scoring matrix
        motif_width: Window size in bp

    Returns:
        List of tuples: (position_1based, score, strand, matched_sequence)
        Positions follow FIMO convention: 1-based, inclusive
    """
    results = []
    seq_len = len(sequence)

    # Forward strand
    for i in range(seq_len - motif_width + 1):
        window = sequence[i:i + motif_width]
        if len(window) < motif_width:
            continue
        sc = score_window(window, log_odds)
        results.append((i + 1, sc, "+", window))

    # Reverse complement strand
    rc_seq = reverse_complement(sequence)
    for i in range(len(rc_seq) - motif_width + 1):
        window = rc_seq[i:i + motif_width]
        if len(window) < motif_width:
            continue
        sc = score_window(window, log_odds)
        # Convert RC position back to forward-strand coordinates
        fwd_start = seq_len - (i + motif_width) + 1  # 1-based
        results.append((fwd_start, sc, "-", window))

    return results


def compute_background_score_distribution(
    sequences: list,
    log_odds: np.ndarray,
    motif_width: int,
    n_samples: int = 10000,
    rng: np.random.Generator = None,
) -> np.ndarray:
    """
    Build an empirical background score distribution by sampling random windows.

    This is used to compute approximate p-values for observed scores.
    A score's p-value = fraction of background windows with score >= observed.

    Args:
        sequences: List of DNA strings
        log_odds: Log-odds scoring matrix
        motif_width: Window size
        n_samples: Number of random windows to sample
        rng: Random generator

    Returns:
        Sorted array of background scores (ascending) for use with np.searchsorted
    """
    if rng is None:
        rng = np.random.default_rng(42)

    all_windows = []
    for seq in sequences:
        if len(seq) >= motif_width:
            for i in range(len(seq) - motif_width + 1):
                all_windows.append(seq[i:i + motif_width])

    if not all_windows:
        return np.array([0.0])

    # Sample with replacement for a smooth distribution
    n_available = len(all_windows)
    sample_indices = rng.integers(0, n_available, size=min(n_samples, n_available * 2))
    sampled = [all_windows[idx % n_available] for idx in sample_indices]

    scores = [score_window(w, log_odds) for w in sampled]
    return np.sort(scores)


def compute_empirical_pvalue(score: float, background_scores: np.ndarray) -> float:
    """
    Compute an empirical p-value for an observed score.

    p-value = fraction of background scores >= observed score.

    Args:
        score: Observed log-odds score
        background_scores: Sorted array of background scores

    Returns:
        Empirical p-value (float between 0 and 1)
    """
    n = len(background_scores)
    # Number of background scores >= observed (right tail)
    n_exceeding = n - np.searchsorted(background_scores, score)
    pvalue = n_exceeding / n
    # Minimum p-value = 1/n to avoid reporting exactly 0
    return max(pvalue, 1.0 / n)


def simulate_fimo_scan(
    motif_file: Path,
    fasta_path: Path,
    output_dir: Path,
    pvalue_threshold: float = 1e-3,
    rng: np.random.Generator = None,
) -> pd.DataFrame:
    """
    Python-based FIMO simulation: scan sequences using log-odds scoring.

    Produces a fimo.tsv file in the exact canonical FIMO output format,
    so analyze_results.py works identically whether real FIMO or this
    simulation generated the TSV.

    Args:
        motif_file: Path to meme.txt
        fasta_path: Path to input FASTA sequences
        output_dir: Directory to write fimo.tsv
        pvalue_threshold: Only report hits with p-value <= this threshold
        rng: Random generator for background sampling

    Returns:
        DataFrame of all hits passing the threshold
    """
    if rng is None:
        rng = np.random.default_rng(42)

    # Parse the motif PWM
    motif_data = parse_meme_motif(motif_file)
    pwm = motif_data["pwm"]
    motif_id = motif_data["motif_id"]
    motif_width = motif_data["width"]
    logger.info(
        f"Loaded motif '{motif_id}' (width={motif_width}) from {motif_file}"
    )

    # Convert PWM to log-odds
    log_odds = pwm_to_log_odds(pwm)

    # Load all sequences
    sequences = []
    seq_names = []
    for record in SeqIO.parse(str(fasta_path), "fasta"):
        sequences.append(str(record.seq).upper())
        seq_names.append(record.id)

    logger.info(f"Scanning {len(sequences)} sequences for motif '{motif_id}'...")

    # Build background score distribution for p-value estimation
    bg_scores = compute_background_score_distribution(
        sequences, log_odds, motif_width, rng=rng
    )
    logger.debug(f"Background score distribution: min={bg_scores[0]:.2f}, max={bg_scores[-1]:.2f}")

    # Scan each sequence
    hits = []
    for seq, name in zip(sequences, seq_names):
        window_scores = score_sequence_with_pwm(seq, log_odds, motif_width)
        for pos, sc, strand, matched in window_scores:
            pval = compute_empirical_pvalue(sc, bg_scores)
            if pval <= pvalue_threshold:
                stop = pos + motif_width - 1
                hits.append(
                    {
                        "motif_id": motif_id,
                        "motif_alt_id": "Mce3R_operator",
                        "sequence_name": name,
                        "start": pos,
                        "stop": stop,
                        "strand": strand,
                        "score": round(sc, 4),
                        "p-value": pval,
                        "q-value": pval,  # Simplified: use p-value as q-value estimate
                        "matched_sequence": matched,
                    }
                )

    hits_df = pd.DataFrame(hits, columns=FIMO_COLUMNS) if hits else pd.DataFrame(columns=FIMO_COLUMNS)

    # Write FIMO-compatible TSV (with comment header lines matching real FIMO format)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    fimo_tsv = output_dir / "fimo.tsv"

    with open(fimo_tsv, "w") as f:
        f.write("# FIMO (simulation mode) motif scan results\n")
        f.write(f"# Motif: {motif_id}\n")
        f.write(f"# Sequences: {fasta_path}\n")
        f.write(f"# p-value threshold: {pvalue_threshold}\n")
        if hits_df.empty:
            f.write("\t".join(FIMO_COLUMNS) + "\n")
        else:
            hits_df.to_csv(f, sep="\t", index=False)

    n_hits = len(hits_df)
    n_seqs_with_hits = hits_df["sequence_name"].nunique() if not hits_df.empty else 0
    logger.info(
        f"Simulation scan complete: {n_hits} hits in {n_seqs_with_hits} sequences "
        f"(p-value <= {pvalue_threshold})"
    )
    logger.info(f"Wrote FIMO results to {fimo_tsv}")

    return hits_df


def run_fimo(
    motif_file: Path,
    fasta_path: Path,
    output_dir: Path,
    pvalue_threshold: float = 1e-3,
    simulate: bool = False,
    rng: np.random.Generator = None,
) -> Path:
    """
    Run FIMO or fall back to the Python log-odds simulation.

    Args:
        motif_file: Path to meme.txt output from run_meme
        fasta_path: Input FASTA sequences
        output_dir: Output directory for FIMO results
        pvalue_threshold: p-value cutoff for reporting hits
        simulate: Force simulation mode
        rng: Random generator for simulation

    Returns:
        Path to the fimo.tsv output file
    """
    output_dir = Path(output_dir)
    fimo_tsv = output_dir / "fimo.tsv"

    fimo_available = check_tool_available("fimo") and not simulate

    if fimo_available:
        logger.info("FIMO found. Running real FIMO motif scan...")
        cmd = [
            "fimo",
            "--oc", str(output_dir),
            "--thresh", str(pvalue_threshold),
            str(motif_file),
            str(fasta_path),
        ]
        logger.info(f"Command: {' '.join(cmd)}")
        try:
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=600,
                check=True,
            )
            if result.stdout:
                logger.debug(f"FIMO stdout:\n{result.stdout[:2000]}")
            logger.info(f"FIMO completed. Output in {output_dir}")
        except subprocess.CalledProcessError as e:
            logger.error(f"FIMO failed with exit code {e.returncode}")
            if e.stderr:
                logger.error(f"FIMO stderr:\n{e.stderr[:1000]}")
            raise
    else:
        if not simulate and not check_tool_available("fimo"):
            logger.warning("FIMO not found in PATH. Running Python log-odds simulation.")
        simulate_fimo_scan(
            motif_file=motif_file,
            fasta_path=fasta_path,
            output_dir=output_dir,
            pvalue_threshold=pvalue_threshold,
            rng=rng,
        )

    return fimo_tsv


def run_fimo_asymmetric(
    fasta_path: Path,
    output_dir: Path,
    motifs_dir: Path,
    pvalue_threshold: float = 1e-3,
    rng: np.random.Generator = None,
) -> Path:
    """
    Scan sequences using the asymmetric Mce3R motif (2024 ACS Chemical Biology model).

    Always uses the Python log-odds simulation (the asymmetric model is hardcoded,
    not discovered by MEME, so results are comparable to the palindromic scan
    without tool availability requirements).

    The asymmetric motif file (meme_asymmetric.txt) must exist in motifs_dir.
    This is created by run_meme.py in simulation mode or can be placed manually.

    Args:
        fasta_path: Input FASTA sequences to scan
        output_dir: Directory to write fimo_asymmetric.tsv
        motifs_dir: Directory containing meme_asymmetric.txt
        pvalue_threshold: p-value cutoff for reporting hits
        rng: Random generator for background sampling

    Returns:
        Path to fimo_asymmetric.tsv, or None if asymmetric motif file not found
    """
    asym_motif_file = Path(motifs_dir) / "meme_asymmetric.txt"

    if not asym_motif_file.exists():
        logger.warning(
            f"Asymmetric motif file not found: {asym_motif_file}\n"
            "  Run the meme step first (creates meme_asymmetric.txt in simulation mode)."
        )
        return None

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    asym_tsv = output_dir / "fimo_asymmetric.tsv"

    logger.info("Running asymmetric model FIMO scan (Mce3R 2024 operator)...")

    # Write to a temporary subdirectory so we don't overwrite the primary fimo.tsv
    import tempfile, shutil
    with tempfile.TemporaryDirectory() as tmp_dir:
        asym_df = simulate_fimo_scan(
            motif_file=asym_motif_file,
            fasta_path=fasta_path,
            output_dir=Path(tmp_dir),
            pvalue_threshold=pvalue_threshold,
            rng=rng if rng is not None else np.random.default_rng(43),
        )
        tmp_fimo = Path(tmp_dir) / "fimo.tsv"
        if tmp_fimo.exists():
            shutil.copy2(str(tmp_fimo), str(asym_tsv))

    n_hits = len(asym_df)
    n_seqs = asym_df["sequence_name"].nunique() if not asym_df.empty else 0
    logger.info(
        f"Asymmetric model scan: {n_hits} hits in {n_seqs} sequences "
        f"→ {asym_tsv}"
    )
    return asym_tsv


def main():
    parser = argparse.ArgumentParser(
        description="Run FIMO motif scanning on M. tuberculosis promoter sequences"
    )
    parser.add_argument(
        "--motif-file",
        type=Path,
        default=get_project_root() / "results" / "motifs" / "meme.txt",
        help="MEME motif file (meme.txt)",
    )
    parser.add_argument(
        "--fasta",
        type=Path,
        default=get_project_root() / "data" / "raw" / "sequences.fasta",
        help="Input FASTA sequences",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=get_project_root() / "results" / "scans",
        help="Output directory for FIMO results",
    )
    parser.add_argument(
        "--pvalue-threshold",
        type=float,
        default=1e-3,
        help="p-value threshold for reporting hits (default: 1e-3)",
    )
    parser.add_argument(
        "--simulate", action="store_true", help="Force simulation mode (skip real FIMO)"
    )
    args = parser.parse_args()

    fimo_tsv = run_fimo(
        motif_file=args.motif_file,
        fasta_path=args.fasta,
        output_dir=args.output_dir,
        pvalue_threshold=args.pvalue_threshold,
        simulate=args.simulate,
    )
    print(f"\nFIMO output: {fimo_tsv}")


if __name__ == "__main__":
    main()
