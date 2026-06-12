"""
build_operator_model.py — Knowledge-based Mce3R operator PWM from REAL sequence.

This replaces the hard-coded consensus (TTGACATNNNNNTGCCCA) that the old pipeline injected.
Here the operator model is built from actual genomic sequence at the experimentally-mapped
operator location:

  positive windows = the mce3R-yrbE3A operator region (real M. tb sequence, sliced from the
                     divergent IGR ~40-170 bp upstream of the yrbE3A start, per the 2024
                     cryo-EM paper) + its orthologous upstream regions across Mycobacterium.

Real MEME (mode 'anr', so it can find BOTH ~25 bp sites per region) is run on these windows
against a genomic background. The resulting PWM is the "knowledge-based" operator model. It
is later cross-checked (Tomtom) against the de novo motif discovered WITHOUT telling MEME
where the operator is; agreement = the operator was recovered, not injected.

Usage:
    python build_operator_model.py --igrs data/raw/divergent_igrs.fasta \
        --orthologs data/raw/ortholog_promoters.fasta --bfile results/background/mtb_genome.bg
"""

import argparse
import sys
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

sys.path.insert(0, str(Path(__file__).resolve().parent))
import mce3r_biology as bio
from run_meme import run_meme
from utils import get_project_root, setup_logging

logger = setup_logging(__name__)


def slice_mtb_operator_window(igrs_fasta: Path) -> SeqRecord | None:
    """
    Slice the M. tb operator window from the mce3R-yrbE3A divergent IGR.

    The IGR is oriented Rv1963c.end -> Rv1964(yrbE3A).start; the yrbE3A start is at the
    high-coordinate (3') end on the + strand, so the operator window is the block closest
    to that end, OPERATOR_WINDOW_UPSTREAM_BP bp upstream of the start. Returns real
    genomic sequence (no consensus hard-coding).
    """
    far, near = bio.OPERATOR_WINDOW_UPSTREAM_BP[1], bio.OPERATOR_WINDOW_UPSTREAM_BP[0]
    for rec in SeqIO.parse(str(igrs_fasta), "fasta"):
        if rec.id == bio.PRIMARY_OPERATOR_IGR:
            seq = str(rec.seq)
            if len(seq) < far:
                window = seq  # short IGR: use whole thing
            else:
                window = seq[len(seq) - far : len(seq) - near]
            return SeqRecord(
                Seq(window),
                id="Mtb_H37Rv_mce3R_yrbE3A_operator",
                description=f"operator_window upstream_yrbE3A {near}-{far}bp source={bio.PRIMARY_OPERATOR_IGR}",
            )
    logger.warning(f"{bio.PRIMARY_OPERATOR_IGR} not found in {igrs_fasta}")
    return None


def assemble_operator_windows(
    igrs_fasta: Path, orthologs_fasta: Path, out_fasta: Path
) -> int:
    """Write the combined operator-window FASTA (M. tb + orthologs). Returns count."""
    records = []
    mtb = slice_mtb_operator_window(igrs_fasta)
    if mtb is not None:
        records.append(mtb)
    if Path(orthologs_fasta).exists():
        for rec in SeqIO.parse(str(orthologs_fasta), "fasta"):
            # Skip a duplicate M. tb ortholog entry if present (we already added the IGR slice).
            if rec.id.startswith("Mtb") and mtb is not None:
                continue
            records.append(rec)
    out_fasta = Path(out_fasta)
    out_fasta.parent.mkdir(parents=True, exist_ok=True)
    SeqIO.write(records, str(out_fasta), "fasta")
    logger.info(f"Assembled {len(records)} operator windows -> {out_fasta}")
    return len(records)


def build_operator_model(
    igrs_fasta: Path,
    orthologs_fasta: Path,
    output_dir: Path,
    bfile: Path | None = None,
    threads: int = 4,
) -> dict:
    """
    Build the knowledge-based operator PWM via anchored MEME on real operator windows.

    Returns dict with: n_windows, windows_fasta, meme_txt, motif (parsed PWM dict).
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    windows_fasta = output_dir / "operator_windows.fasta"
    n = assemble_operator_windows(igrs_fasta, orthologs_fasta, windows_fasta)
    if n < 2:
        raise RuntimeError(
            f"Only {n} operator window(s) assembled — need >=2. "
            "Run 'fetch_orthologs' to add more mycobacterial yrbE3A regions."
        )
    if n < 8:
        logger.warning(
            f"Only {n} operator windows available; the knowledge PWM and the conservation "
            "gate (G5, target >=8 orthologs) will be under-powered. Fetch more orthologs."
        )

    w = bio.OPERATOR_SITE_WIDTH
    meme_dir = output_dir / "meme"
    meme_txt = run_meme(
        windows_fasta,
        meme_dir,
        nmotifs=2,
        minw=max(8, w - 4),
        maxw=w + 4,
        mod="anr",
        threads=threads,
        bfile=bfile,
    )

    # Copy the primary motif file to a stable name for downstream FIMO / Tomtom.
    knowledge_meme = output_dir / "operator_knowledge.meme"
    knowledge_meme.write_text(Path(meme_txt).read_text())

    from run_meme import parse_meme_motif

    motif = parse_meme_motif(knowledge_meme, motif_index=0)
    logger.info(
        f"Knowledge-based operator PWM: id={motif['motif_id']} w={motif['width']} "
        f"nsites={motif['nsites']} E={motif['evalue']}"
    )
    return {
        "n_windows": n,
        "windows_fasta": windows_fasta,
        "meme_txt": knowledge_meme,
        "motif": motif,
    }


def main():
    root = get_project_root()
    parser = argparse.ArgumentParser(
        description="Build the knowledge-based Mce3R operator PWM"
    )
    parser.add_argument(
        "--igrs", type=Path, default=root / "data" / "raw" / "divergent_igrs.fasta"
    )
    parser.add_argument(
        "--orthologs",
        type=Path,
        default=root / "data" / "raw" / "ortholog_promoters.fasta",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=root / "results" / "motifs" / "operator_knowledge",
    )
    parser.add_argument("--bfile", type=Path, default=None)
    parser.add_argument("--threads", type=int, default=4)
    args = parser.parse_args()

    result = build_operator_model(
        args.igrs, args.orthologs, args.output_dir, args.bfile, args.threads
    )
    print(f"\nKnowledge-based operator model: {result['meme_txt']}")
    print(f"  windows used : {result['n_windows']}")
    print(
        f"  motif width  : {result['motif']['width']}  E-value: {result['motif']['evalue']}"
    )


if __name__ == "__main__":
    main()
