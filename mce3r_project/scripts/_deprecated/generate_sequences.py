"""
generate_sequences.py — Synthetic M. tuberculosis promoter sequence generator.

This script generates realistic synthetic promoter/intergenic sequences for
Mce3R binding motif discovery. Since real ChIP-seq or SELEX data for Mce3R
may not be available, we create sequences that:

1. Have ~65% GC content, matching the M. tuberculosis H37Rv genome
2. Use a first-order Markov chain with dinucleotide bias (more realistic
   than purely random GC-rich sequences)
3. Contain embedded Mce3R-like binding motifs in a subset of sequences

Biological context:
    Mce3R is a TetR-family transcription factor that represses the mce3
    operon (Rv1970-Rv1978), which encodes mammalian cell entry proteins
    important for M. tuberculosis virulence. TetR-family proteins typically
    bind palindromic inverted repeat sequences (~18-20 bp) with the structure:
        [half-site 1] - [spacer] - [half-site 2, reverse complement of half-site 1]

    The consensus motif used here (TTGACANNNNNTGTCAA) is based on the
    canonical TetR operator structure, where TGTCAA is the reverse complement
    of TTGACA.

Usage:
    python generate_sequences.py [--n-sequences N] [--seed S] [--output PATH]
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Add project root to path for imports
sys.path.insert(0, str(Path(__file__).resolve().parent))
from utils import get_project_root, setup_logging, compute_gc_content

logger = setup_logging(__name__)

# The TetR-family consensus motif for Mce3R
# TTGACA: left half-site
# NNNNN:  5-base variable spacer (filled with GC-biased random bases)
# TGTCAA: right half-site (reverse complement of TTGACA)
MCE3R_MOTIF_CONSENSUS = "TTGACANNNNNTGTCAA"

# Dinucleotide transition matrix approximating M. tuberculosis H37Rv genome
# Rows: current nucleotide (A, C, G, T)
# Columns: next nucleotide (A, C, G, T)
# Values: transition probabilities (each row sums to 1.0)
# Source: derived from published M. tb genome composition analyses
MTB_TRANSITION_MATRIX = {
    "A": {"A": 0.12, "C": 0.30, "G": 0.38, "T": 0.20},
    "C": {"A": 0.13, "C": 0.35, "G": 0.35, "T": 0.17},
    "G": {"A": 0.22, "C": 0.30, "G": 0.32, "T": 0.16},
    "T": {"A": 0.18, "C": 0.28, "G": 0.38, "T": 0.16},
}

# Stationary distribution (initial nucleotide probabilities)
# Reflects the ~65% GC content of M. tb genome
MTB_BASE_PROBABILITIES = {"A": 0.175, "C": 0.325, "G": 0.325, "T": 0.175}


def generate_gc_rich_sequence(length: int, rng: np.random.Generator) -> str:
    """
    Generate one GC-rich DNA sequence using a first-order Markov chain.

    The Markov chain uses transition probabilities derived from the
    M. tuberculosis H37Rv genome, producing sequences with realistic
    dinucleotide composition (~65% GC content).

    Args:
        length: Desired sequence length in base pairs
        rng: NumPy random generator for reproducibility

    Returns:
        Uppercase DNA string of specified length
    """
    bases = list(MTB_BASE_PROBABILITIES.keys())
    probs = list(MTB_BASE_PROBABILITIES.values())

    # Draw the first nucleotide from the stationary distribution
    current = rng.choice(bases, p=probs)
    sequence = [current]

    for _ in range(length - 1):
        trans = MTB_TRANSITION_MATRIX[current]
        next_bases = list(trans.keys())
        next_probs = list(trans.values())
        current = rng.choice(next_bases, p=next_probs)
        sequence.append(current)

    return "".join(sequence)


def generate_motif_instance(
    consensus: str = MCE3R_MOTIF_CONSENSUS,
    rng: np.random.Generator = None,
) -> str:
    """
    Instantiate a concrete binding motif from the consensus pattern.

    Each 'N' position in the consensus is filled with a GC-biased random
    nucleotide, reflecting that TetR spacers are not under strict sequence
    conservation but still tend toward the host genome's GC composition.

    Args:
        consensus: Motif consensus string with N at degenerate positions
        rng: NumPy random generator (uses global default if None)

    Returns:
        Concrete DNA string of same length as consensus
    """
    if rng is None:
        rng = np.random.default_rng()

    bases = list(MTB_BASE_PROBABILITIES.keys())
    probs = list(MTB_BASE_PROBABILITIES.values())

    result = []
    for char in consensus:
        if char == "N":
            result.append(rng.choice(bases, p=probs))
        else:
            result.append(char)

    return "".join(result)


def embed_motif(sequence: str, motif: str, position: int) -> tuple:
    """
    Insert a motif into a background sequence at a specified position.

    The motif replaces (overwrites) sequence bases at the given position.
    This simulates a TetR operator that evolved within a promoter region.

    Args:
        sequence: Background DNA sequence
        motif: Concrete motif string to insert
        position: 0-based start position for insertion

    Returns:
        Tuple of (modified_sequence, actual_position)

    Raises:
        ValueError: If position + len(motif) exceeds sequence length
    """
    if position + len(motif) > len(sequence):
        raise ValueError(
            f"Motif (len={len(motif)}) at position {position} "
            f"exceeds sequence length {len(sequence)}"
        )
    modified = sequence[:position] + motif + sequence[position + len(motif):]
    return modified, position


def generate_synthetic_dataset(
    n_sequences: int = 25,
    n_with_motif: int = 14,
    min_length: int = 200,
    max_length: int = 300,
    motif_consensus: str = MCE3R_MOTIF_CONSENSUS,
    seed: int = 42,
    output_path: Path = None,
) -> pd.DataFrame:
    """
    Generate a complete synthetic dataset of M. tb-like promoter sequences.

    Creates n_sequences background sequences and embeds the Mce3R motif
    into n_with_motif of them. The sequences are written to a FASTA file
    and metadata is returned as a DataFrame.

    The 'zoops' distribution (zero or one occurrence per sequence) matches
    the MEME discovery mode used in run_meme.py.

    Args:
        n_sequences: Total number of sequences to generate
        n_with_motif: How many sequences should contain the embedded motif
        min_length: Minimum sequence length in bp
        max_length: Maximum sequence length in bp
        motif_consensus: Consensus pattern for the TetR binding motif
        seed: Random seed for reproducibility
        output_path: Path to write the FASTA file

    Returns:
        DataFrame with columns: sequence_id, sequence, length, gc_content,
        has_motif, motif_position
    """
    rng = np.random.default_rng(seed)
    motif_width = len(motif_consensus.replace("N", "N"))  # = full length including N

    if n_with_motif > n_sequences:
        raise ValueError(f"n_with_motif ({n_with_motif}) cannot exceed n_sequences ({n_sequences})")

    # Randomly choose which sequences will receive a motif
    motif_indices = set(rng.choice(n_sequences, size=n_with_motif, replace=False).tolist())

    records = []
    seq_records = []

    for i in range(n_sequences):
        seq_id = f"Mtb_promoter_{i + 1:03d}"
        length = int(rng.integers(min_length, max_length + 1))

        # Generate GC-rich background sequence
        sequence = generate_gc_rich_sequence(length, rng)

        has_motif = i in motif_indices
        motif_position = None

        if has_motif:
            # Embed motif at a random internal position (away from sequence edges)
            # Margin ensures the motif is well within the sequence for accurate
            # positional analysis in downstream FIMO scanning
            margin = 50
            embed_range_start = margin
            embed_range_end = length - motif_width - margin

            if embed_range_end <= embed_range_start:
                # Sequence too short for safe embedding — extend minimum margin
                embed_range_start = 10
                embed_range_end = max(11, length - motif_width - 10)

            position = int(rng.integers(embed_range_start, embed_range_end + 1))
            motif_instance = generate_motif_instance(motif_consensus, rng)
            sequence, motif_position = embed_motif(sequence, motif_instance, position)

            logger.debug(
                f"{seq_id}: embedded motif '{motif_instance}' at position {motif_position}"
            )

        gc = compute_gc_content(sequence)

        # Build FASTA description with metadata for downstream parsing
        if has_motif:
            description = (
                f"gc_content={gc:.3f} length={length} "
                f"has_motif=True motif_pos={motif_position}"
            )
        else:
            description = f"gc_content={gc:.3f} length={length} has_motif=False"

        records.append(
            {
                "sequence_id": seq_id,
                "sequence": sequence,
                "length": length,
                "gc_content": gc,
                "has_motif": has_motif,
                "motif_position": motif_position,
            }
        )

        # Build BioPython SeqRecord for FASTA output
        seq_record = SeqRecord(
            Seq(sequence),
            id=seq_id,
            description=description,
        )
        seq_records.append(seq_record)

    df = pd.DataFrame(records)

    # Write FASTA file
    if output_path is not None:
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        SeqIO.write(seq_records, str(output_path), "fasta")
        logger.info(f"Wrote {len(seq_records)} sequences to {output_path}")

    n_actual_motif = df["has_motif"].sum()
    mean_gc = df["gc_content"].mean()
    logger.info(
        f"Generated {n_sequences} sequences: "
        f"{n_actual_motif} with motif, "
        f"mean GC={mean_gc:.3f}, "
        f"length range {min_length}-{max_length} bp"
    )

    return df


def main():
    parser = argparse.ArgumentParser(
        description="Generate synthetic M. tuberculosis promoter sequences for Mce3R motif discovery"
    )
    parser.add_argument(
        "--n-sequences", type=int, default=25, help="Total number of sequences (default: 25)"
    )
    parser.add_argument(
        "--n-with-motif",
        type=int,
        default=14,
        help="Sequences containing the Mce3R motif (default: 14)",
    )
    parser.add_argument(
        "--min-length", type=int, default=200, help="Minimum sequence length in bp (default: 200)"
    )
    parser.add_argument(
        "--max-length", type=int, default=300, help="Maximum sequence length in bp (default: 300)"
    )
    parser.add_argument("--seed", type=int, default=42, help="Random seed (default: 42)")
    parser.add_argument(
        "--output",
        type=Path,
        default=get_project_root() / "data" / "raw" / "sequences.fasta",
        help="Output FASTA file path",
    )
    parser.add_argument(
        "--metadata-output",
        type=Path,
        default=get_project_root() / "data" / "processed" / "sequence_metadata.csv",
        help="Output CSV for sequence metadata",
    )
    args = parser.parse_args()

    df = generate_synthetic_dataset(
        n_sequences=args.n_sequences,
        n_with_motif=args.n_with_motif,
        min_length=args.min_length,
        max_length=args.max_length,
        seed=args.seed,
        output_path=args.output,
    )

    args.metadata_output.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.metadata_output, index=False)
    logger.info(f"Saved sequence metadata to {args.metadata_output}")

    print(f"\nSequence summary:")
    print(f"  Total sequences : {len(df)}")
    print(f"  With motif      : {df['has_motif'].sum()}")
    print(f"  Mean GC content : {df['gc_content'].mean():.3f}")
    print(f"  Length range    : {df['length'].min()}-{df['length'].max()} bp")
    print(f"\nOutput files:")
    print(f"  FASTA    : {args.output}")
    print(f"  Metadata : {args.metadata_output}")


if __name__ == "__main__":
    main()
