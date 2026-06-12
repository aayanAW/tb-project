"""
discover_denovo.py — De novo operator discovery + non-circularity cross-check.

The key test that the previous pipeline failed: can the operator be discovered WITHOUT
telling the tool where it is, and does that match the knowledge-based model?

  positive set = the experimentally Mce3R-bound divergent IGRs (mce3R-yrbE3A, with ~6 motif
                 occurrences per Santangelo 2009; Rv1935c-Rv1936) + orthologous upstream
                 regions. MEME mode 'anr' finds any number of sites per sequence.
  negative set = random divergent IGRs (NOT the operator regions), matched by length, as a
                 discriminative control (-neg). This is what stops MEME from "discovering"
                 generic GC-rich repeats.

Then Tomtom compares the de novo motif to the knowledge-based operator PWM. A significant
match (q < 0.05) means the operator was recovered de novo — Gate G1.
"""

import argparse
import subprocess
import sys
from pathlib import Path

import numpy as np
from Bio import SeqIO

sys.path.insert(0, str(Path(__file__).resolve().parent))
import mce3r_biology as bio
from run_meme import extract_evalues_from_meme, run_meme
from utils import get_project_root, require_meme_tool, setup_logging

logger = setup_logging(__name__)

N_CONTROL_IGRS = 30
CONTROL_SEED = 42


def build_positive_set(igrs_fasta: Path, orthologs_fasta: Path, out_fasta: Path) -> int:
    """Positive set = the two bound operator IGRs + ortholog upstream regions."""
    records = []
    igr_index = {r.id: r for r in SeqIO.parse(str(igrs_fasta), "fasta")}
    for igr_id in bio.OPERATOR_IGRS:
        if igr_id in igr_index:
            records.append(igr_index[igr_id])
        else:
            logger.warning(f"Operator IGR not found in {igrs_fasta.name}: {igr_id}")
    if Path(orthologs_fasta).exists():
        records += list(SeqIO.parse(str(orthologs_fasta), "fasta"))
    Path(out_fasta).parent.mkdir(parents=True, exist_ok=True)
    SeqIO.write(records, str(out_fasta), "fasta")
    logger.info(f"De novo positive set: {len(records)} sequences -> {out_fasta}")
    return len(records)


def build_control_set(
    igrs_fasta: Path, out_fasta: Path, n: int = N_CONTROL_IGRS
) -> int:
    """Negative control = random divergent IGRs excluding the operator regions."""
    candidates = [
        r
        for r in SeqIO.parse(str(igrs_fasta), "fasta")
        if r.id not in bio.OPERATOR_IGRS
    ]
    rng = np.random.default_rng(CONTROL_SEED)
    if len(candidates) > n:
        idx = rng.choice(len(candidates), size=n, replace=False)
        chosen = [candidates[i] for i in sorted(idx)]
    else:
        chosen = candidates
    Path(out_fasta).parent.mkdir(parents=True, exist_ok=True)
    SeqIO.write(chosen, str(out_fasta), "fasta")
    logger.info(f"De novo negative/control set: {len(chosen)} IGRs -> {out_fasta}")
    return len(chosen)


def run_tomtom(
    query_meme: Path, target_meme: Path, output_dir: Path, thresh: float = 0.05
) -> dict:
    """
    Compare a de novo motif to the knowledge-based operator PWM with Tomtom.

    Returns dict: matched (bool), best_qvalue (float|None), tsv (path).
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    tomtom = require_meme_tool("tomtom")
    cmd = [
        tomtom,
        "-no-ssc",
        "-oc",
        str(output_dir),
        "-thresh",
        str(thresh),
        "-dist",
        "pearson",
        str(query_meme),
        str(target_meme),
    ]
    logger.info(f"Tomtom: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
    tsv = output_dir / "tomtom.tsv"
    if not tsv.exists():
        logger.warning(
            f"Tomtom produced no tsv (exit {result.returncode}): {result.stderr[-300:]}"
        )
        return {"matched": False, "best_qvalue": None, "tsv": tsv}

    best_q = None
    with open(tsv) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        try:
            q_idx = header.index("q-value")
        except ValueError:
            q_idx = 5
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) <= q_idx:
                continue
            try:
                q = float(parts[q_idx])
            except ValueError:
                continue
            best_q = q if best_q is None else min(best_q, q)
    matched = best_q is not None and best_q < thresh
    logger.info(f"Tomtom best q-value = {best_q} -> matched={matched}")
    return {"matched": matched, "best_qvalue": best_q, "tsv": tsv}


def discover_denovo(
    igrs_fasta: Path,
    orthologs_fasta: Path,
    knowledge_meme: Path,
    output_dir: Path,
    bfile: Path | None = None,
    threads: int = 4,
) -> dict:
    """
    Full de novo discovery + G1 cross-check. Returns dict with denovo motif + G1 verdict.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    pos = output_dir / "positive_set.fasta"
    neg = output_dir / "control_set.fasta"
    n_pos = build_positive_set(igrs_fasta, orthologs_fasta, pos)
    n_neg = build_control_set(igrs_fasta, neg)

    w = bio.OPERATOR_SITE_WIDTH
    meme_dir = output_dir / "meme"
    meme_txt = run_meme(
        pos,
        meme_dir,
        nmotifs=3,
        minw=max(8, w - 4),
        maxw=w + 4,
        mod="anr",
        threads=threads,
        bfile=bfile,
        neg_fasta=neg,
        objfun="de",
    )
    denovo_meme = output_dir / "operator_denovo.meme"
    denovo_meme.write_text(Path(meme_txt).read_text())

    g1 = run_tomtom(denovo_meme, knowledge_meme, output_dir / "tomtom")
    motifs = extract_evalues_from_meme(denovo_meme)
    return {
        "n_positive": n_pos,
        "n_control": n_neg,
        "denovo_meme": denovo_meme,
        "denovo_motifs": motifs,
        "g1_matched": g1["matched"],
        "g1_best_qvalue": g1["best_qvalue"],
        "g1_tomtom_tsv": str(g1["tsv"]),
    }


def main():
    root = get_project_root()
    parser = argparse.ArgumentParser(
        description="De novo operator discovery + Tomtom cross-check"
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
        "--knowledge",
        type=Path,
        default=root
        / "results"
        / "motifs"
        / "operator_knowledge"
        / "operator_knowledge.meme",
    )
    parser.add_argument(
        "--output-dir", type=Path, default=root / "results" / "motifs" / "denovo"
    )
    parser.add_argument("--bfile", type=Path, default=None)
    parser.add_argument("--threads", type=int, default=4)
    args = parser.parse_args()

    res = discover_denovo(
        args.igrs,
        args.orthologs,
        args.knowledge,
        args.output_dir,
        args.bfile,
        args.threads,
    )
    print(f"\nDe novo discovery: {res['denovo_meme']}")
    print(f"  positive seqs: {res['n_positive']}  control IGRs: {res['n_control']}")
    print(
        f"  G1 (de novo matches knowledge PWM, Tomtom q<0.05): "
        f"{'PASS' if res['g1_matched'] else 'FAIL'}  (best q={res['g1_best_qvalue']})"
    )


if __name__ == "__main__":
    main()
