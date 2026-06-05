"""
utils.py — Shared infrastructure for the Mce3R motif discovery pipeline.

This module provides helper functions used across all pipeline scripts:
- Logging setup
- Tool availability checking (MEME, FIMO)
- Directory management
- FASTA loading into pandas DataFrames

Biological context: Mce3R is a TetR-family transcription factor in
Mycobacterium tuberculosis that represses the mce3 operon, which encodes
proteins involved in mammalian cell entry during infection.
"""

import logging
import shutil
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import SeqIO


def get_project_root() -> Path:
    """Return the absolute path to mce3r_project/ regardless of invocation location."""
    # This file lives at mce3r_project/scripts/utils.py
    return Path(__file__).resolve().parent.parent


def setup_logging(name: str, level: str = "INFO") -> logging.Logger:
    """
    Create a named logger with a consistent format.

    Args:
        name: Logger name (typically the calling module's __name__)
        level: Logging level string ("DEBUG", "INFO", "WARNING", "ERROR")

    Returns:
        Configured Logger instance
    """
    logger = logging.getLogger(name)
    if not logger.handlers:
        handler = logging.StreamHandler(sys.stdout)
        formatter = logging.Formatter(
            "[%(asctime)s] [%(name)s] [%(levelname)s] %(message)s",
            datefmt="%Y-%m-%d %H:%M:%S",
        )
        handler.setFormatter(formatter)
        logger.addHandler(handler)
    logger.setLevel(getattr(logging, level.upper(), logging.INFO))
    return logger


def check_tool_available(tool_name: str) -> bool:
    """
    Check whether a command-line tool is available in the system PATH.

    Uses shutil.which() for cross-platform compatibility.

    Args:
        tool_name: Binary name to search for (e.g., "meme", "fimo")

    Returns:
        True if tool is found in PATH, False otherwise
    """
    return shutil.which(tool_name) is not None


def ensure_directories(paths: list) -> None:
    """
    Create all required directories, including any missing parents.

    Args:
        paths: List of Path objects to create
    """
    for path in paths:
        Path(path).mkdir(parents=True, exist_ok=True)


def compute_gc_content(sequence: str) -> float:
    """
    Calculate GC content of a DNA sequence.

    GC content is biologically important for M. tuberculosis sequences,
    which have an unusually high GC content (~65%) compared to most bacteria.
    High GC content affects DNA melting temperature and transcription factor binding.

    Args:
        sequence: DNA string (uppercase or lowercase, IUPAC bases allowed)

    Returns:
        Fraction of G+C bases (0.0 to 1.0). Returns 0.0 for empty sequences.
    """
    if not sequence:
        return 0.0
    seq_upper = sequence.upper()
    gc = sum(1 for base in seq_upper if base in "GC")
    return gc / len(seq_upper)


def load_fasta_to_dataframe(fasta_path: Path) -> pd.DataFrame:
    """
    Load a FASTA file into a pandas DataFrame.

    Each sequence becomes a row with metadata extracted from the header.
    This function also parses optional metadata fields embedded in the
    FASTA description line (e.g., has_motif=True motif_pos=91).

    Args:
        fasta_path: Path to the FASTA file

    Returns:
        DataFrame with columns:
            sequence_id (str): FASTA record ID
            description (str): Full FASTA header description
            sequence (str): DNA sequence (uppercase)
            length (int): Sequence length in bp
            gc_content (float): GC fraction
            has_motif (bool or None): Parsed from header if present
            motif_position (int or None): Parsed from header if present
    """
    fasta_path = Path(fasta_path)
    if not fasta_path.exists():
        raise FileNotFoundError(f"FASTA file not found: {fasta_path}")

    records = []
    for record in SeqIO.parse(str(fasta_path), "fasta"):
        seq_str = str(record.seq).upper()
        desc = record.description

        # Parse optional metadata fields from the description line
        has_motif = None
        motif_position = None
        for token in desc.split():
            if token.startswith("has_motif="):
                has_motif = token.split("=")[1].lower() == "true"
            elif token.startswith("motif_pos="):
                try:
                    motif_position = int(token.split("=")[1])
                except ValueError:
                    pass

        records.append(
            {
                "sequence_id": record.id,
                "description": desc,
                "sequence": seq_str,
                "length": len(seq_str),
                "gc_content": compute_gc_content(seq_str),
                "has_motif": has_motif,
                "motif_position": motif_position,
            }
        )

    if not records:
        # Return empty DataFrame with correct schema
        return pd.DataFrame(
            columns=[
                "sequence_id",
                "description",
                "sequence",
                "length",
                "gc_content",
                "has_motif",
                "motif_position",
            ]
        )

    return pd.DataFrame(records)


def reverse_complement(sequence: str) -> str:
    """
    Compute the reverse complement of a DNA sequence.

    Used in palindrome detection: TetR-family operators are palindromic
    because the protein binds as a homodimer, one monomer per half-site
    on opposite strands.

    Args:
        sequence: DNA string (uppercase, IUPAC ambiguity codes supported)

    Returns:
        Reverse complement string (uppercase)
    """
    complement_map = str.maketrans(
        "ACGTacgtNnRrYyKkMmSsWwBbDdHhVv",
        "TGCAtgcaNnYyRrMmKkSsWwVvHhDdBb",
    )
    return sequence.translate(complement_map)[::-1]
