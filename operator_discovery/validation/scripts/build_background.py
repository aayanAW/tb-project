"""
build_background.py — Build a genomic background model for MEME and FIMO.

The previous pipeline scanned a 65.6% GC genome against an effectively uniform /
self-sampled background, which manufactures false positives (any GC-rich window scores
high). The correct background reflects actual genome composition. We use the MEME Suite
`fasta-get-markov` tool to estimate an order-k Markov background directly from the
H37Rv genome FASTA, and write it in the format MEME/FIMO consume via `-bfile`.

Usage:
    python build_background.py [--genome PATH] [--order 2] [--output PATH]
"""

import argparse
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from utils import get_project_root, require_meme_tool, setup_logging

logger = setup_logging(__name__)

DEFAULT_ORDER = 2  # order-2 (trinucleotide) captures GC-rich di/tri structure


def build_background(
    genome_fasta: Path, output_path: Path, order: int = DEFAULT_ORDER
) -> Path:
    """
    Estimate an order-k Markov background from a genome FASTA.

    Args:
        genome_fasta: Path to the genome FASTA (M. tb H37Rv NC_000962.3).
        output_path: Where to write the MEME background file.
        order: Markov order (0 = single-nucleotide freqs, 2 = trinucleotide).

    Returns:
        Path to the written background file.
    """
    genome_fasta = Path(genome_fasta)
    if not genome_fasta.exists():
        raise FileNotFoundError(
            f"Genome FASTA not found: {genome_fasta}\nRun the 'download' step first."
        )
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    tool = require_meme_tool("fasta-get-markov")
    cmd = [tool, "-m", str(order), "-dna", str(genome_fasta)]
    logger.info(f"Building order-{order} genomic background: {' '.join(cmd)}")

    with open(output_path, "w") as out:
        result = subprocess.run(cmd, stdout=out, stderr=subprocess.PIPE, text=True)
    if result.returncode != 0:
        logger.error(f"fasta-get-markov failed: {result.stderr[-500:]}")
        raise subprocess.CalledProcessError(result.returncode, cmd)

    n_lines = sum(1 for _ in open(output_path))
    logger.info(f"Wrote genomic background ({n_lines} entries) -> {output_path}")
    return output_path


def main():
    root = get_project_root()
    parser = argparse.ArgumentParser(
        description="Build a genomic background model for MEME/FIMO"
    )
    parser.add_argument(
        "--genome", type=Path, default=root / "data" / "raw" / "Mtb_H37Rv.fasta"
    )
    parser.add_argument("--order", type=int, default=DEFAULT_ORDER)
    parser.add_argument(
        "--output", type=Path, default=root / "results" / "background" / "mtb_genome.bg"
    )
    args = parser.parse_args()

    out = build_background(args.genome, args.output, args.order)
    print(f"\nBackground model: {out}")


if __name__ == "__main__":
    main()
