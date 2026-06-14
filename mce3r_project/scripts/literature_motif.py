"""
literature_motif.py — Build a MEME-format PWM from INDEPENDENTLY-PUBLISHED operator
sites (external reference for the non-circular G1 corroboration test).

The audit (C1) showed the old G1 — Tomtom of the de novo motif vs the knowledge PWM —
was circular: both motifs are MEME runs over the same operator DNA + the same 8 ortholog
windows, so their agreement is guaranteed and proves nothing about independence.

The fix: compare the knowledge-based PWM against a PWM built from operator half-site
sequences that OTHER labs reported (Santangelo 2009 footprinting; mce3r_biology.
SANTANGELO_2009_OPERATOR_SITES). Those sequences are not produced by this pipeline, so a
significant Tomtom match is genuine external corroboration that the PWM captures the real
operator rather than a self-consistent artifact.

This module only converts published, equal-length site strings into a count-based PWM and
writes it in MEME minimal motif format. No genomic sequence from this pipeline enters it.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import mce3r_biology as bio
from utils import get_project_root, setup_logging

logger = setup_logging(__name__)

_ALPHABET = ("A", "C", "G", "T")
_IDX = {b: i for i, b in enumerate(_ALPHABET)}


def sites_to_pwm(sites: tuple[str, ...], pseudocount: float = 0.25) -> np.ndarray:
    """
    Convert equal-length aligned site strings to a probability matrix (w, 4) over ACGT.

    A small pseudocount is added so no column is degenerate (required for Tomtom scoring).
    Non-ACGT characters are skipped for that column.
    """
    sites = tuple(s.strip().upper() for s in sites if s.strip())
    if not sites:
        raise ValueError("No site sequences supplied.")
    width = len(sites[0])
    if any(len(s) != width for s in sites):
        raise ValueError(
            f"Published sites are not equal length: {[len(s) for s in sites]}. "
            "Align/trim them before building a PWM."
        )
    counts = np.full((width, 4), pseudocount, dtype=float)
    for site in sites:
        for j, base in enumerate(site):
            if base in _IDX:
                counts[j, _IDX[base]] += 1.0
    return counts / counts.sum(axis=1, keepdims=True)


def write_meme_pwm(
    pwm: np.ndarray,
    out_meme: Path,
    motif_id: str,
    nsites: int,
    background: tuple[float, float, float, float] = (0.17, 0.33, 0.33, 0.17),
) -> Path:
    """
    Write a single PWM to MEME minimal motif format (v4).

    The default background reflects M. tuberculosis' ~65% GC composition so the
    external motif is scored on the same footing as the genomic models.
    """
    out_meme = Path(out_meme)
    out_meme.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        "MEME version 4",
        "",
        "ALPHABET= ACGT",
        "",
        "strands: + -",
        "",
        "Background letter frequencies",
        "A {:.5f} C {:.5f} G {:.5f} T {:.5f}".format(*background),
        "",
        f"MOTIF {motif_id}",
        f"letter-probability matrix: alength= 4 w= {pwm.shape[0]} nsites= {nsites} E= 0",
    ]
    for row in pwm:
        lines.append(" " + " ".join(f"{p:.6f}" for p in row))
    lines.append("")
    out_meme.write_text("\n".join(lines))
    logger.info(f"External literature PWM ({pwm.shape[0]} bp) -> {out_meme}")
    return out_meme


def build_literature_motif(out_meme: Path) -> Path:
    """Build the Santangelo-2009 external operator PWM and write it in MEME format."""
    sites = bio.SANTANGELO_2009_OPERATOR_SITES
    pwm = sites_to_pwm(sites)
    return write_meme_pwm(
        pwm,
        out_meme,
        motif_id="Santangelo2009_operator",
        nsites=len(sites),
    )


def main():
    root = get_project_root()
    parser = argparse.ArgumentParser(
        description="Build external (published) operator PWM in MEME format"
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=root / "results" / "motifs" / "external" / "santangelo2009.meme",
    )
    args = parser.parse_args()
    path = build_literature_motif(args.output)
    print(f"\nExternal literature motif: {path}")
    print(f"  source: {bio.SANTANGELO_2009_CITATION}")
    print(f"  sites : {len(bio.SANTANGELO_2009_OPERATOR_SITES)}")


if __name__ == "__main__":
    main()
