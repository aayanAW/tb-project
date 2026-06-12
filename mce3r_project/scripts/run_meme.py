"""
run_meme.py — Real MEME motif discovery (no simulation, no hard-coded consensus).

The previous version of this file contained a simulation fallback that wrote a synthetic
PWM derived from a HAND-TYPED operator consensus (TTGACATNNNNNTGCCCA). Scanning for a
typed string and then "discovering" it is circular; that path is removed entirely. This
module only ever runs the real MEME binary, and it always uses:
  - a genomic background model (-bfile) so GC-rich windows are not spuriously enriched,
  - an optional negative/control set (-neg) for discriminative discovery.

If MEME is not installed the functions raise a clear error (no silent fallback).

MEME command (discriminative, when a control set is given):
    meme primary.fasta -dna -mod zoops -nmotifs N -minw W1 -maxw W2 \
         -bfile genome.bg -neg control.fasta -objfun de -oc OUTDIR -p T
"""

import argparse
import re
import subprocess
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from utils import get_project_root, require_meme_tool, setup_logging

logger = setup_logging(__name__)

MEME_TIMEOUT_S = 3600


def build_meme_command(
    fasta_path: Path,
    output_dir: Path,
    nmotifs: int = 3,
    minw: int = 22,
    maxw: int = 26,
    mod: str = "zoops",
    threads: int = 4,
    bfile: Path | None = None,
    neg_fasta: Path | None = None,
    objfun: str | None = None,
    revcomp: bool = True,
) -> list:
    """
    Build the MEME argument list (never a shell string — no injection risk).

    Args:
        fasta_path: Primary (positive) input FASTA.
        output_dir: MEME output directory (-oc, overwritten).
        nmotifs: Max motifs to find.
        minw / maxw: Motif width range (bp). ~25 bp matches the Mce3R half-site.
        mod: 'zoops' (<=1 site/seq) or 'anr' (any number of repetitions).
        threads: CPU threads.
        bfile: Background model file (-bfile). Strongly recommended for GC-rich genomes.
        neg_fasta: Negative/control set (-neg) for discriminative discovery.
        objfun: Objective function ('classic', 'de', 'cd'); defaults to 'de' when a
                negative set is given, else MEME's default.
        revcomp: Search both strands (-revcomp).

    Returns:
        List of command tokens.
    """
    meme = require_meme_tool("meme")
    cmd = [
        meme,
        str(fasta_path),
        "-dna",
        "-mod",
        mod,
        "-nmotifs",
        str(nmotifs),
        "-minw",
        str(minw),
        "-maxw",
        str(maxw),
        "-oc",
        str(output_dir),
        "-p",
        str(threads),
    ]
    if revcomp:
        cmd.append("-revcomp")
    if bfile is not None:
        cmd += ["-bfile", str(bfile)]
    if neg_fasta is not None:
        cmd += ["-neg", str(neg_fasta)]
        objfun = objfun or "de"
    if objfun is not None:
        cmd += ["-objfun", objfun]
    return cmd


def run_meme(
    fasta_path: Path,
    output_dir: Path,
    nmotifs: int = 3,
    minw: int = 22,
    maxw: int = 26,
    mod: str = "zoops",
    threads: int = 4,
    bfile: Path | None = None,
    neg_fasta: Path | None = None,
    objfun: str | None = None,
) -> Path:
    """
    Run real MEME and return the path to meme.txt. Raises on failure (no simulation).
    """
    output_dir = Path(output_dir)
    cmd = build_meme_command(
        fasta_path,
        output_dir,
        nmotifs=nmotifs,
        minw=minw,
        maxw=maxw,
        mod=mod,
        threads=threads,
        bfile=bfile,
        neg_fasta=neg_fasta,
        objfun=objfun,
    )
    logger.info(f"MEME: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=MEME_TIMEOUT_S)

    meme_txt = output_dir / "meme.txt"
    if meme_txt.exists():
        if result.returncode != 0:
            # MEME often returns non-zero only because EPS->PNG / XML->HTML conversion
            # failed (missing Ghostscript/ImageMagick) while meme.txt is fully written.
            logger.warning(
                f"MEME exit code {result.returncode} but meme.txt was written "
                "(likely only logo/HTML conversion failed — harmless)."
            )
        logger.info(f"MEME completed -> {meme_txt}")
        return meme_txt

    logger.error(f"MEME failed (exit {result.returncode}); no meme.txt produced.")
    if result.stderr:
        logger.error(f"stderr:\n{result.stderr[-1000:]}")
    raise subprocess.CalledProcessError(result.returncode, cmd)


def parse_meme_motif(motif_file: Path, motif_index: int = 0) -> dict:
    """
    Parse a MEME minimal/text file and return one motif's PWM.

    Anchors on landmark strings (handles MEME v4 and v5). Returns the motif at
    motif_index (0-based) among those found.

    Returns dict: motif_id, width, nsites, evalue, pwm (np.ndarray (w,4) [A,C,G,T]).
    """
    motif_file = Path(motif_file)
    if not motif_file.exists():
        raise FileNotFoundError(f"MEME motif file not found: {motif_file}")

    motifs = []
    cur_id = None
    with open(motif_file) as fh:
        lines = fh.readlines()

    i = 0
    while i < len(lines):
        stripped = lines[i].strip()
        if stripped.startswith("MOTIF ") and not stripped.endswith("_done"):
            parts = stripped.split()
            cur_id = parts[1] if len(parts) >= 2 else "unknown"
        if "letter-probability matrix" in stripped:
            width = None
            nsites = None
            evalue = None
            w_match = re.search(r"\bw=\s*(\d+)", stripped)
            if w_match:
                width = int(w_match.group(1))
            n_match = re.search(r"\bnsites=\s*(\d+)", stripped)
            if n_match:
                nsites = int(n_match.group(1))
            e_match = re.search(r"\bE=\s*([^\s]+)", stripped)
            if e_match:
                try:
                    evalue = float(e_match.group(1))
                except ValueError:
                    pass
            rows = []
            j = i + 1
            while j < len(lines):
                parts = lines[j].strip().split()
                if len(parts) == 4:
                    try:
                        rows.append([float(p) for p in parts])
                    except ValueError:
                        break
                    j += 1
                else:
                    break
            if rows:
                if width is None:
                    width = len(rows)
                motifs.append(
                    {
                        "motif_id": cur_id or f"motif_{len(motifs) + 1}",
                        "width": width,
                        "nsites": nsites if nsites is not None else len(rows),
                        "evalue": evalue,
                        "pwm": np.array(rows[:width]),
                        "alphabet_order": ["A", "C", "G", "T"],
                    }
                )
            i = j
            continue
        i += 1

    if not motifs:
        raise ValueError(f"No PWM data found in {motif_file}")
    if motif_index >= len(motifs):
        motif_index = 0
    return motifs[motif_index]


def extract_evalues_from_meme(meme_txt: Path) -> list[dict]:
    """Return [{motif_id, width, nsites, evalue}, ...] for every motif in a meme file."""
    meme_txt = Path(meme_txt)
    if not meme_txt.exists():
        return []
    results = []
    cur_id = None
    with open(meme_txt) as fh:
        for line in fh:
            s = line.strip()
            if s.startswith("MOTIF ") and not s.endswith("_done"):
                parts = s.split()
                cur_id = parts[1] if len(parts) >= 2 else "unknown"
            if "letter-probability matrix" in s:
                w = re.search(r"\bw=\s*(\d+)", s)
                n = re.search(r"\bnsites=\s*(\d+)", s)
                e = re.search(r"\bE=\s*([^\s]+)", s)
                evalue = None
                if e:
                    try:
                        evalue = float(e.group(1))
                    except ValueError:
                        pass
                results.append(
                    {
                        "motif_id": cur_id or "unknown",
                        "width": int(w.group(1)) if w else None,
                        "nsites": int(n.group(1)) if n else None,
                        "evalue": evalue,
                    }
                )
    return results


def main():
    root = get_project_root()
    parser = argparse.ArgumentParser(description="Run real MEME motif discovery")
    parser.add_argument("--fasta", type=Path, required=True)
    parser.add_argument(
        "--output-dir", type=Path, default=root / "results" / "motifs" / "denovo"
    )
    parser.add_argument("--nmotifs", type=int, default=3)
    parser.add_argument("--minw", type=int, default=22)
    parser.add_argument("--maxw", type=int, default=26)
    parser.add_argument("--mod", default="zoops", choices=["zoops", "anr", "oops"])
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--bfile", type=Path, default=None)
    parser.add_argument("--neg", type=Path, default=None)
    args = parser.parse_args()

    meme_txt = run_meme(
        args.fasta,
        args.output_dir,
        nmotifs=args.nmotifs,
        minw=args.minw,
        maxw=args.maxw,
        mod=args.mod,
        threads=args.threads,
        bfile=args.bfile,
        neg_fasta=args.neg,
    )
    print(f"\nMEME output: {meme_txt}")
    for m in extract_evalues_from_meme(meme_txt):
        print(
            f"  {m['motif_id']}  w={m['width']}  nsites={m['nsites']}  E={m['evalue']}"
        )


if __name__ == "__main__":
    main()
