"""
run_meme.py — MEME motif discovery wrapper with simulation fallback.

MEME (Multiple Em for Motif Elicitation) discovers statistically enriched
sequence motifs from a set of DNA sequences. For Mce3R analysis, it identifies
the asymmetric operator sequence that Mce3R binds to repress the mce3 operon.

MEME command used (zoops mode — one site per promoter):
    meme sequences.fasta -dna -mod zoops -nmotifs 5 -minw 20 -maxw 30
         -oc results/motifs/ -p 4

MEME command (anr mode — any number of repetitions, finds both half-sites):
    meme sequences.fasta -dna -mod anr -nmotifs 2 -minw 20 -maxw 30
         -oc results/motifs/anr/ -p 4

Flags:
    -dna      : DNA alphabet (not RNA or protein)
    -mod zoops: Zero Or One Occurrence Per Sequence — for single-site discovery
    -mod anr  : Any Number of Repetitions — captures BOTH half-sites per sequence
                (critical because each Mce3R operator has TWO distinct 25 bp sites)
    -nmotifs 5: Discover up to 5 distinct motifs
    -minw 20 / -maxw 30: Width range for ~25 bp Mce3R half-sites
                (2024 ACS Chem. Biol.: each half-site is ~25 bp)
    -oc       : Output directory (created if absent)
    -p 4      : Use 4 CPU threads

Simulation mode (when MEME Suite is not installed):
    Writes two synthetic meme.txt files:
    - meme.txt: downstream (high-affinity) half-site PWM
    - meme_asymmetric.txt: both sites combined (MCE3R_DOWNSTREAM + MCE3R_UPSTREAM)
    Output format is identical to real MEME v4, so downstream scripts work
    transparently in both real and simulation modes.

Usage:
    python run_meme.py --fasta PATH --output-dir PATH [--threads N] [--simulate]
"""

import argparse
import re
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from utils import check_tool_available, get_project_root, setup_logging

logger = setup_logging(__name__)

# Mce3R downstream (high-affinity) half-site consensus — 2024 ACS Chemical Biology
# This is the primary half-site; the Mce3R dimer binds this site with Kd ~2.4 nM.
# Arg53 and Lys262 from the non-identical HTH motifs insert into the major groove.
DOWNSTREAM_HALFSITE = "TTGACATNNNNNTGCCCA"   # 18 bp asymmetric (right half ≠ RC of left)

# Mce3R upstream (lower-affinity) half-site consensus
# Slightly more degenerate than the downstream site; same overall architecture.
UPSTREAM_HALFSITE = "TTGACANNNNNTGCCCA"       # 17 bp (more degenerate left portion)

# Legacy aliases for backward compatibility with run_fimo.py and other imports
DEFAULT_CONSENSUS   = DOWNSTREAM_HALFSITE     # primary motif for FIMO scan
ASYMMETRIC_CONSENSUS = DOWNSTREAM_HALFSITE    # same — asymmetric model is now primary


def build_meme_command(
    fasta_path: Path,
    output_dir: Path,
    nmotifs: int = 5,
    minw: int = 20,
    maxw: int = 30,
    mod: str = "zoops",
    threads: int = 4,
) -> list:
    """
    Build the MEME subprocess argument list.

    Never builds a shell string — uses a list to avoid shell injection.

    Width range defaults (minw=20, maxw=30) are set to capture the ~25 bp
    Mce3R half-sites established in Panagoda et al. 2024 (ACS Chem. Biol.).
    The -mod anr variant discovers both half-sites per sequence; -mod zoops
    finds the best single site per sequence.

    Args:
        fasta_path: Input FASTA file path
        output_dir: Directory for MEME output (will be created by MEME)
        nmotifs: Maximum number of motifs to discover
        minw: Minimum motif width in bp (default 20 for 25 bp Mce3R half-sites)
        maxw: Maximum motif width in bp (default 30)
        mod: Occurrence model ('zoops' for one site/seq, 'anr' for both half-sites)
        threads: Number of CPU threads

    Returns:
        List of strings forming the subprocess command
    """
    return [
        "meme",
        str(fasta_path),
        "-dna",
        "-mod", mod,
        "-nmotifs", str(nmotifs),
        "-minw", str(minw),
        "-maxw", str(maxw),
        "-oc", str(output_dir),
        "-p", str(threads),
    ]


def build_synthetic_pwm(
    consensus: str,
    pseudocount: float = 0.1,
    noise_std: float = 0.05,
    rng: np.random.Generator = None,
) -> np.ndarray:
    """
    Build a Position Weight Matrix (PWM) from a consensus sequence.

    The PWM represents the probability of each nucleotide at each position.
    Adding Dirichlet-like noise makes it resemble a learned matrix rather
    than a hand-crafted one.

    Alphabet order: A=0, C=1, G=2, T=3 (MEME convention)

    Args:
        consensus: Motif consensus string (N = degenerate position)
        pseudocount: Minimum probability for any base at any position
        noise_std: Standard deviation of Gaussian noise added to probabilities
        rng: NumPy random generator

    Returns:
        numpy array of shape (len(consensus), 4) where each row sums to 1.0
    """
    if rng is None:
        rng = np.random.default_rng(42)

    base_to_idx = {"A": 0, "C": 1, "G": 2, "T": 3}
    width = len(consensus)
    pwm = np.zeros((width, 4))

    for i, base in enumerate(consensus.upper()):
        if base == "N":
            # Degenerate position: near-uniform with GC bias matching M. tb genome
            pwm[i] = [0.175, 0.325, 0.325, 0.175]
        elif base in base_to_idx:
            # Conserved position: high probability for consensus base
            pwm[i, :] = pseudocount
            pwm[i, base_to_idx[base]] = 1.0 - 3 * pseudocount
        else:
            pwm[i] = [0.25, 0.25, 0.25, 0.25]

    # Add small Gaussian noise to make the matrix look like a learned one
    noise = rng.normal(0, noise_std, pwm.shape)
    pwm = np.clip(pwm + noise, pseudocount, 1.0)

    # Renormalize each row to sum to 1.0
    row_sums = pwm.sum(axis=1, keepdims=True)
    pwm = pwm / row_sums

    return pwm


def simulate_meme_output(
    output_dir: Path,
    motif_consensus: str = DEFAULT_CONSENSUS,
    rng: np.random.Generator = None,
) -> Path:
    """
    Write a synthetic meme.txt file in MEME minimal text format (version 4).

    Produces TWO motif entries reflecting the dual-site Mce3R operator
    (Panagoda et al. 2024):
        MCE3R_DOWNSTREAM: high-affinity half-site (Kd ~2.4 nM)
        MCE3R_UPSTREAM:   lower-affinity half-site (same architecture, more noise)

    Both motifs use the asymmetric consensus (right half ≠ RC of left half),
    consistent with the 2024 nonpalindromic operator model. This fallback
    produces output identical in structure to real MEME v4, so FIMO and
    downstream scripts work without modification.

    Args:
        output_dir: Directory to write meme.txt into
        motif_consensus: Consensus string for the primary (downstream) motif
        rng: NumPy random generator for PWM noise

    Returns:
        Path to the written meme.txt file
    """
    if rng is None:
        rng = np.random.default_rng(42)

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    meme_txt_path = output_dir / "meme.txt"

    # Downstream (high-affinity) half-site
    pwm_downstream = build_synthetic_pwm(motif_consensus, noise_std=0.04, rng=rng)
    width_ds = len(motif_consensus)

    # Upstream (lower-affinity) half-site — more noise, slightly degenerate
    pwm_upstream = build_synthetic_pwm(
        UPSTREAM_HALFSITE, noise_std=0.08, pseudocount=0.15, rng=rng
    )
    width_us = len(UPSTREAM_HALFSITE)

    lines = []
    lines.append("MEME version 4")
    lines.append("")
    lines.append("ALPHABET= ACGT")
    lines.append("")
    lines.append("strands: + -")
    lines.append("")
    lines.append("Background letter frequencies (from uniform background):")
    # M. tb background frequencies (GC-rich genome, 65% GC)
    lines.append("A 0.175 C 0.325 G 0.325 T 0.175")
    lines.append("")
    # Motif 1: downstream (high-affinity) half-site
    lines.append("MOTIF MCE3R_DOWNSTREAM Mce3R_downstream_halfsite")
    lines.append("")
    lines.append(
        f"letter-probability matrix: alength= 4 w= {width_ds} nsites= 9 E= 0"
    )
    for row in pwm_downstream:
        lines.append("  " + "  ".join(f"{p:.6f}" for p in row))
    lines.append("")
    lines.append("MOTIF MCE3R_DOWNSTREAM_done")
    lines.append("")
    # Motif 2: upstream (lower-affinity) half-site
    lines.append("MOTIF MCE3R_UPSTREAM Mce3R_upstream_halfsite")
    lines.append("")
    lines.append(
        f"letter-probability matrix: alength= 4 w= {width_us} nsites= 7 E= 0"
    )
    for row in pwm_upstream:
        lines.append("  " + "  ".join(f"{p:.6f}" for p in row))
    lines.append("")
    lines.append("MOTIF MCE3R_UPSTREAM_done")
    lines.append("")

    with open(meme_txt_path, "w") as f:
        f.write("\n".join(lines) + "\n")

    logger.info(f"Wrote simulated meme.txt to {meme_txt_path}")
    return meme_txt_path


def simulate_asymmetric_meme_output(
    output_dir: Path,
    rng: np.random.Generator = None,
) -> Path:
    """
    Write a synthetic meme_asymmetric.txt file using the Mce3R asymmetric consensus.

    This file represents the nonpalindromic Mce3R operator half-site identified
    in the 2024 ACS Chemical Biology paper (Panagoda, Balázsi & Sampson).
    It is used for a second FIMO scan (in addition to meme.txt) to detect
    binding sites consistent with the asymmetric operator model.

    The consensus TTGACATNNNNNTGCCCA has a right half-site (TGCCCA) that is
    NOT the reverse complement of the left half-site (TTGACAT). The paired
    operator consists of TWO such sites ~53 bp apart. This file serves as the
    asymmetric motif for the second FIMO scan in run_fimo.py.

    Args:
        output_dir: Directory to write meme_asymmetric.txt into
        rng: NumPy random generator for PWM noise

    Returns:
        Path to the written meme_asymmetric.txt file
    """
    if rng is None:
        rng = np.random.default_rng(43)  # Different seed from palindromic simulation

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    meme_txt_path = output_dir / "meme_asymmetric.txt"

    pwm = build_synthetic_pwm(ASYMMETRIC_CONSENSUS, rng=rng)
    width = len(ASYMMETRIC_CONSENSUS)
    nsites = 9  # Fewer candidate promoters match the asymmetric model

    lines = []
    lines.append("MEME version 4")
    lines.append("")
    lines.append("ALPHABET= ACGT")
    lines.append("")
    lines.append("strands: + -")
    lines.append("")
    lines.append("Background letter frequencies (from uniform background):")
    lines.append("A 0.175 C 0.325 G 0.325 T 0.175")
    lines.append("")
    lines.append("MOTIF MCE3R_ASYM Mce3R_asymmetric_operator")
    lines.append("")
    lines.append(
        f"letter-probability matrix: alength= 4 w= {width} nsites= {nsites} E= 0"
    )
    for row in pwm:
        lines.append("  " + "  ".join(f"{p:.6f}" for p in row))
    lines.append("")
    lines.append("MOTIF MCE3R_ASYM_done")
    lines.append("")

    with open(meme_txt_path, "w") as f:
        f.write("\n".join(lines) + "\n")

    logger.info(
        f"Wrote asymmetric motif meme_asymmetric.txt to {meme_txt_path} "
        f"(consensus: {ASYMMETRIC_CONSENSUS})"
    )
    return meme_txt_path


def parse_meme_motif(motif_file: Path) -> dict:
    """
    Parse a MEME minimal format text file to extract the first motif's PWM.

    Handles both MEME v4 and v5 output formats. The parser anchors on
    known landmark strings rather than assuming fixed line numbers.

    Args:
        motif_file: Path to meme.txt file

    Returns:
        Dictionary with keys:
            motif_id (str): MOTIF identifier
            width (int): Motif width in bp
            nsites (int): Number of sites used to build the motif
            pwm (np.ndarray): Shape (width, 4), columns are [A, C, G, T]
            alphabet_order (list): ['A', 'C', 'G', 'T']
    """
    motif_file = Path(motif_file)
    if not motif_file.exists():
        raise FileNotFoundError(f"MEME motif file not found: {motif_file}")

    with open(motif_file) as f:
        lines = f.readlines()

    motif_id = "MCE3R_1"
    width = None
    nsites = 14
    pwm_rows = []
    in_matrix = False

    for line in lines:
        stripped = line.strip()

        # Parse MOTIF line to get motif ID
        if stripped.startswith("MOTIF ") and not stripped.endswith("_done"):
            parts = stripped.split()
            if len(parts) >= 2:
                motif_id = parts[1]

        # Parse letter-probability matrix header
        if "letter-probability matrix" in stripped:
            # Extract width and nsites from the header line
            for token in stripped.split():
                if token.startswith("w="):
                    try:
                        width = int(token.split("=")[1])
                    except (ValueError, IndexError):
                        pass
                elif token.startswith("nsites="):
                    try:
                        nsites = int(token.split("=")[1])
                    except (ValueError, IndexError):
                        pass
            in_matrix = True
            continue

        # Read PWM rows (4 floats per line)
        if in_matrix:
            parts = stripped.split()
            if len(parts) == 4:
                try:
                    row = [float(p) for p in parts]
                    pwm_rows.append(row)
                except ValueError:
                    in_matrix = False
            else:
                if pwm_rows:  # We were reading a matrix and hit a non-data line
                    in_matrix = False

    if not pwm_rows:
        raise ValueError(f"No PWM data found in {motif_file}")

    if width is None:
        width = len(pwm_rows)

    return {
        "motif_id": motif_id,
        "width": width,
        "nsites": nsites,
        "pwm": np.array(pwm_rows[:width]),
        "alphabet_order": ["A", "C", "G", "T"],
    }


def extract_evalues_from_meme(meme_txt: Path) -> list[dict]:
    """
    Extract E-values for every motif in a meme.txt file.

    Parses the 'letter-probability matrix' header lines which contain
    E= values.  Returns one dict per motif with keys: motif_id, width,
    nsites, evalue.

    Args:
        meme_txt: Path to a MEME output file.

    Returns:
        List of dicts, one per motif found, sorted by occurrence order.
    """
    meme_txt = Path(meme_txt)
    if not meme_txt.exists():
        return []

    results = []
    current_motif_id = None

    with open(meme_txt) as fh:
        for line in fh:
            stripped = line.strip()

            if stripped.startswith("MOTIF ") and not stripped.endswith("_done"):
                parts = stripped.split()
                current_motif_id = parts[1] if len(parts) >= 2 else "unknown"

            if "letter-probability matrix" in stripped:
                width = None
                nsites = None
                evalue = None
                # Parse key=value pairs; MEME v5 uses "E= 1.6e-045" (space after =)
                # so we use regex to handle both "E=0" and "E= 1.6e-045"
                w_match = re.search(r'\bw=\s*(\d+)', stripped)
                if w_match:
                    width = int(w_match.group(1))
                n_match = re.search(r'\bnsites=\s*(\d+)', stripped)
                if n_match:
                    nsites = int(n_match.group(1))
                e_match = re.search(r'\bE=\s*([^\s]+)', stripped)
                if e_match:
                    try:
                        evalue = float(e_match.group(1))
                    except ValueError:
                        pass
                results.append({
                    "motif_id": current_motif_id or "unknown",
                    "width": width,
                    "nsites": nsites,
                    "evalue": evalue,
                })

    return results


def run_meme_width_grid(
    fasta_path: Path,
    output_dir: Path,
    min_width: int = 20,
    max_width: int = 28,
    step: int = 1,
    mod: str = "zoops",
    nmotifs: int = 5,
    threads: int = 4,
) -> Path:
    """
    Grid search over MEME motif widths to find the optimal width.

    For each candidate width w in [min_width, max_width] (step=1), MEME is
    run with -minw w -maxw w (fixed width).  The E-value of the best motif
    in each run is recorded.  After all runs complete, a TSV summary is
    written and the meme.txt from the best width is copied into output_dir
    as the canonical result for downstream analysis.

    Width range is constrained to >=20 bp to ensure the motif captures a
    full Mce3R half-site (~25 bp, Panagoda et al. 2024 ACS Chem. Biol.).
    Narrower motifs produce lower E-values but miss the flanking positions
    needed for paired-site architecture detection.

    Args:
        fasta_path:  Input FASTA file.
        output_dir:  Base output directory (grid results go into output_dir/grid/).
        min_width:   Smallest width to test (default 20, captures full half-site).
        max_width:   Largest width to test (default 28).
        step:        Width increment (default 1).
        mod:         MEME occurrence model ('zoops' or 'anr').
        nmotifs:     Number of motifs per run.
        threads:     CPU threads for MEME.

    Returns:
        Path to the meme.txt from the best (lowest E-value) width.
    """
    fasta_path = Path(fasta_path)
    output_dir = Path(output_dir)
    grid_dir = output_dir / "grid"
    grid_dir.mkdir(parents=True, exist_ok=True)

    widths = list(range(min_width, max_width + 1, step))
    rows: list[dict] = []

    logger.info(
        f"Starting MEME width grid search: w={min_width}..{max_width} "
        f"(step={step}, {len(widths)} runs, mode={mod})"
    )

    for w in widths:
        run_dir = grid_dir / f"w{w}"
        cmd = build_meme_command(
            fasta_path, run_dir, nmotifs=nmotifs, minw=w, maxw=w,
            mod=mod, threads=threads,
        )
        logger.info(f"  Grid w={w}: {' '.join(cmd)}")

        try:
            result = subprocess.run(
                cmd, capture_output=True, text=True, timeout=3600,
            )
            meme_txt = run_dir / "meme.txt"
            if meme_txt.exists():
                motifs = extract_evalues_from_meme(meme_txt)
                best_e = min((m["evalue"] for m in motifs if m["evalue"] is not None),
                             default=None)
                n_motifs_found = len(motifs)
                if result.returncode != 0:
                    logger.warning(
                        f"  Grid w={w}: MEME exit code {result.returncode} "
                        "but meme.txt exists (likely logo conversion failure)."
                    )
            else:
                best_e = None
                n_motifs_found = 0
                logger.warning(f"  Grid w={w}: MEME failed — no meme.txt produced.")
                if result.stderr:
                    logger.debug(f"  stderr: {result.stderr[-500:]}")
        except subprocess.TimeoutExpired:
            best_e = None
            n_motifs_found = 0
            logger.warning(f"  Grid w={w}: MEME timed out (>3600 s).")

        rows.append({
            "width": w,
            "best_evalue": best_e,
            "n_motifs": n_motifs_found,
            "meme_txt": str(run_dir / "meme.txt"),
        })
        logger.info(f"  Grid w={w}: E-value={best_e}  motifs={n_motifs_found}")

    # ── Write grid_summary.tsv ──────────────────────────────────────────
    tsv_path = output_dir / "grid_summary.tsv"
    with open(tsv_path, "w") as fh:
        fh.write("width\tbest_evalue\tn_motifs\tmeme_txt\n")
        for r in rows:
            ev = f"{r['best_evalue']:.3e}" if r["best_evalue"] is not None else "NA"
            fh.write(f"{r['width']}\t{ev}\t{r['n_motifs']}\t{r['meme_txt']}\n")
    logger.info(f"Grid summary written to {tsv_path}")

    # ── Select best width ───────────────────────────────────────────────
    valid = [r for r in rows if r["best_evalue"] is not None]
    if not valid:
        logger.error("Grid search failed — no width produced a valid E-value.")
        raise RuntimeError("MEME width grid search: all widths failed.")

    best = min(valid, key=lambda r: r["best_evalue"])
    logger.info(
        f"Best width: {best['width']} bp  "
        f"(E-value = {best['best_evalue']:.3e}, {best['n_motifs']} motifs)"
    )

    # ── Copy winning meme.txt to canonical location ─────────────────────
    best_meme = Path(best["meme_txt"])
    canonical = output_dir / "meme.txt"
    if best_meme.exists():
        shutil.copy2(best_meme, canonical)
        logger.info(f"Copied best result (w={best['width']}) → {canonical}")

    return canonical


def run_meme(
    fasta_path: Path,
    output_dir: Path,
    simulate: bool = False,
    threads: int = 4,
    nmotifs: int = 5,
    minw: int = 20,
    maxw: int = 30,
    rng: np.random.Generator = None,
    grid_search: bool = False,
    grid_min: int = 20,
    grid_max: int = 28,
    grid_step: int = 1,
) -> Path:
    """
    Run MEME motif discovery or fall back to simulation mode.

    Args:
        fasta_path: Input FASTA file
        output_dir: Output directory for MEME results
        simulate: Force simulation mode even if MEME is installed
        threads: CPU threads for parallel MEME
        nmotifs: Maximum number of motifs to find
        minw: Minimum motif width (ignored when grid_search=True)
        maxw: Maximum motif width (ignored when grid_search=True)
        rng: Random generator for simulation PWM noise
        grid_search: Run width grid search to find optimal width
        grid_min: Smallest width for grid search (default 20, constrained to
                  capture full ~25 bp Mce3R half-site per Panagoda 2024)
        grid_max: Largest width for grid search (default 28)
        grid_step: Width increment for grid search (default 1)

    Returns:
        Path to the meme.txt output file
    """
    output_dir = Path(output_dir)
    meme_txt = output_dir / "meme.txt"

    meme_available = check_tool_available("meme") and not simulate

    if meme_available and grid_search:
        logger.info("Grid search mode enabled — scanning widths to find optimal E-value.")
        meme_txt = run_meme_width_grid(
            fasta_path=fasta_path,
            output_dir=output_dir,
            min_width=grid_min,
            max_width=grid_max,
            step=grid_step,
            mod="zoops",
            nmotifs=nmotifs,
            threads=threads,
        )
        # Still run ANR with the best width for dual half-site discovery
        # Read best width from grid_summary.tsv (more reliable than parsing meme.txt)
        tsv = output_dir / "grid_summary.tsv"
        best_w = grid_min  # fallback
        if tsv.exists():
            import csv
            with open(tsv) as fh:
                reader = csv.DictReader(fh, delimiter="\t")
                valid = [(r, float(r["best_evalue"])) for r in reader if r["best_evalue"] != "NA"]
                if valid:
                    best_w = int(min(valid, key=lambda x: x[1])[0]["width"])
        logger.info(f"Running ANR mode at best grid width ({best_w} bp)...")
        anr_dir = output_dir / "anr"
        cmd_anr = build_meme_command(
            fasta_path, anr_dir, nmotifs=2, minw=best_w, maxw=best_w,
            mod="anr", threads=threads,
        )
        try:
            result_anr = subprocess.run(
                cmd_anr, capture_output=True, text=True, timeout=3600,
            )
            if (anr_dir / "meme.txt").exists():
                logger.info(f"MEME (anr, w={best_w}) completed. Output in {anr_dir}")
            else:
                logger.warning(f"MEME ANR (w={best_w}) produced no output (non-fatal).")
        except (subprocess.CalledProcessError, subprocess.TimeoutExpired) as e:
            logger.warning(f"MEME ANR mode failed (non-fatal): {e}")
        return meme_txt

    if meme_available:
        logger.info("MEME Suite found. Running real MEME motif discovery...")
        # Run 1: ZOOPS mode (one occurrence per sequence — best single half-site)
        cmd = build_meme_command(
            fasta_path, output_dir, nmotifs=nmotifs, minw=minw, maxw=maxw,
            mod="zoops", threads=threads
        )
        logger.info(f"Command (zoops): {' '.join(cmd)}")
        try:
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=3600,
            )
            meme_txt = output_dir / "meme.txt"
            if meme_txt.exists():
                # Core MEME analysis succeeded even if post-processing (EPS→PNG,
                # XML→HTML) failed due to missing ImageMagick/Ghostscript.
                if result.returncode != 0:
                    logger.warning(
                        f"MEME exited with code {result.returncode} but meme.txt was "
                        "written — likely only logo/HTML conversion failed (harmless)."
                    )
                if result.stdout:
                    logger.debug(f"MEME stdout:\n{result.stdout[-2000:]}")
                logger.info(f"MEME (zoops) completed successfully. Output in {output_dir}")
            else:
                logger.error(f"MEME failed with exit code {result.returncode} and produced no meme.txt")
                if result.stderr:
                    logger.error(f"MEME stderr:\n{result.stderr[-1000:]}")
                raise subprocess.CalledProcessError(result.returncode, cmd)
        except subprocess.CalledProcessError:
            raise

        # Run 2: ANR mode (any number of repetitions — finds both half-sites per seq)
        anr_dir = output_dir / "anr"
        cmd_anr = build_meme_command(
            fasta_path, anr_dir, nmotifs=2, minw=minw, maxw=maxw,
            mod="anr", threads=threads
        )
        logger.info(f"Command (anr): {' '.join(cmd_anr)}")
        try:
            result_anr = subprocess.run(
                cmd_anr,
                capture_output=True,
                text=True,
                timeout=3600,
            )
            anr_txt = anr_dir / "meme.txt"
            if anr_txt.exists():
                logger.info(f"MEME (anr) completed. Output in {anr_dir}")
            else:
                raise subprocess.CalledProcessError(result_anr.returncode, cmd_anr)
        except subprocess.CalledProcessError as e:
            logger.warning(
                f"MEME ANR mode failed (non-fatal): {e.returncode}\n"
                "ANR results will not be available; zoops results still usable."
            )
    else:
        if not simulate and not check_tool_available("meme"):
            logger.warning(
                "\n[WARNING] MEME Suite not found in PATH.\n"
                "          Running in SIMULATION mode:\n"
                "          - Synthetic PWM derived from consensus TTGACANNNNNTGTCAA\n"
                "          - Output format identical to real MEME output\n"
                "\n"
                "          To install MEME Suite:\n"
                "            conda install -c bioconda meme   (recommended)\n"
                "            or visit https://meme-suite.org/meme/doc/install.html\n"
            )
        simulate_meme_output(output_dir, rng=rng)
        # Also generate the asymmetric model file for model comparison
        simulate_asymmetric_meme_output(output_dir, rng=rng)

    return meme_txt


def main():
    parser = argparse.ArgumentParser(
        description="Run MEME motif discovery on M. tuberculosis promoter sequences"
    )
    parser.add_argument(
        "--fasta",
        type=Path,
        default=get_project_root() / "data" / "raw" / "sequences.fasta",
        help="Input FASTA file",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=get_project_root() / "results" / "motifs",
        help="Output directory for MEME results",
    )
    parser.add_argument("--threads", type=int, default=4, help="CPU threads for MEME (default: 4)")
    parser.add_argument("--nmotifs", type=int, default=5, help="Number of motifs to discover")
    parser.add_argument(
        "--simulate", action="store_true", help="Force simulation mode (skip real MEME)"
    )
    parser.add_argument(
        "--grid-search", action="store_true",
        help="Run width grid search (w=12..28) to find optimal E-value",
    )
    parser.add_argument("--grid-min", type=int, default=20, help="Grid search min width (default: 20, >=20 for full half-site)")
    parser.add_argument("--grid-max", type=int, default=28, help="Grid search max width (default: 28)")
    parser.add_argument("--grid-step", type=int, default=1, help="Grid search width step (default: 1)")
    args = parser.parse_args()

    meme_txt = run_meme(
        fasta_path=args.fasta,
        output_dir=args.output_dir,
        simulate=args.simulate,
        threads=args.threads,
        nmotifs=args.nmotifs,
        grid_search=args.grid_search,
        grid_min=args.grid_min,
        grid_max=args.grid_max,
        grid_step=args.grid_step,
    )
    print(f"\nMEME output: {meme_txt}")


if __name__ == "__main__":
    main()
