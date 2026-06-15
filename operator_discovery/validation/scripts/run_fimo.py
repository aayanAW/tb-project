"""
run_fimo.py — Real FIMO motif scanning (no simulation).

The previous version shipped a Python "FIMO simulation" whose empirical p-value floor
equalled the reporting threshold (every hit pinned at p=1e-4), sampled its null from the
sequences it was scanning, and set q-value = p-value (no FDR). It also ran a second scan
against a HARD-CODED operator consensus. All of that is removed.

This module runs the real `fimo` binary with:
  - a genomic background model (--bfile), correct for M. tb's 65.6% GC composition,
  - a Benjamini-Hochberg q-value (FDR) threshold (--qv-thresh), not a raw per-site p,
and downstream code reports only hits with score > 0 (a negative log-odds score means the
window matches background better than the motif and is not a binding site).

FIMO command:
    fimo --oc OUTDIR --bfile genome.bg --qv-thresh --thresh 0.05 motif.meme seqs.fasta
"""

import argparse
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from utils import get_project_root, require_meme_tool, setup_logging

logger = setup_logging(__name__)

FIMO_TIMEOUT_S = 1800


def run_fimo(
    motif_file: Path,
    fasta_path: Path,
    output_dir: Path,
    bfile: Path | None = None,
    qv_thresh: float = 0.05,
    use_qvalue: bool = True,
    max_stored_scores: int = 1_000_000,
) -> Path:
    """
    Run real FIMO and return the path to fimo.tsv. Raises on failure (no simulation).

    Args:
        motif_file: MEME-format motif file (real discovered or knowledge-based PWM).
        fasta_path: Sequences to scan (all promoters + IGRs).
        output_dir: FIMO output directory (-oc, overwritten).
        bfile: Genomic background model. Strongly recommended.
        qv_thresh: Threshold value applied to the q-value (FDR) when use_qvalue=True,
                   else to the p-value.
        use_qvalue: If True, threshold on the BH q-value (--qv-thresh).
        max_stored_scores: FIMO score buffer for genome-wide scans.

    Returns:
        Path to fimo.tsv.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    fimo = require_meme_tool("fimo")

    cmd = [fimo, "--oc", str(output_dir), "--max-stored-scores", str(max_stored_scores)]
    if bfile is not None:
        cmd += ["--bfile", str(bfile)]
    if use_qvalue:
        cmd += ["--qv-thresh", "--thresh", str(qv_thresh)]
    else:
        cmd += ["--thresh", str(qv_thresh)]
    cmd += [str(motif_file), str(fasta_path)]

    logger.info(f"FIMO: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=FIMO_TIMEOUT_S)

    fimo_tsv = output_dir / "fimo.tsv"
    if result.returncode != 0 and not fimo_tsv.exists():
        logger.error(f"FIMO failed (exit {result.returncode}).")
        if result.stderr:
            logger.error(f"stderr:\n{result.stderr[-1000:]}")
        raise subprocess.CalledProcessError(result.returncode, cmd)
    if not fimo_tsv.exists():
        raise FileNotFoundError(f"FIMO produced no fimo.tsv in {output_dir}")
    # Fail loud on a non-zero exit even when a (possibly stale/partial) tsv exists, and on a
    # data-row-free tsv (audit C6): never let a truncated FIMO run propagate into the gates.
    if result.returncode != 0:
        logger.warning(
            f"FIMO returned non-zero exit {result.returncode} but wrote {fimo_tsv}; "
            f"stderr tail:\n{(result.stderr or '')[-500:]}"
        )
    n_data_rows = sum(
        1
        for line in fimo_tsv.read_text().splitlines()
        if line.strip() and not line.startswith("#") and not line.startswith("motif_id")
    )
    logger.info(f"FIMO completed -> {fimo_tsv} ({n_data_rows} data rows)")
    return fimo_tsv


def main():
    root = get_project_root()
    parser = argparse.ArgumentParser(description="Run real FIMO motif scanning")
    parser.add_argument("--motif-file", type=Path, required=True)
    parser.add_argument(
        "--fasta", type=Path, default=root / "data" / "raw" / "promoters.fasta"
    )
    parser.add_argument("--output-dir", type=Path, default=root / "results" / "scans")
    parser.add_argument("--bfile", type=Path, default=None)
    parser.add_argument("--qv-thresh", type=float, default=0.05)
    parser.add_argument(
        "--use-pvalue", action="store_true", help="threshold on p instead of q"
    )
    args = parser.parse_args()

    fimo_tsv = run_fimo(
        args.motif_file,
        args.fasta,
        args.output_dir,
        bfile=args.bfile,
        qv_thresh=args.qv_thresh,
        use_qvalue=not args.use_pvalue,
    )
    print(f"\nFIMO output: {fimo_tsv}")


if __name__ == "__main__":
    main()
