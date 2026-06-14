"""
main.py — Corrected Mce3R operator-scan pipeline (valid, non-circular).

This replaces the previous orchestrator whose results were circular (MEME on 9 positive-only
promoters), injected (a hard-coded operator consensus scanned by a self-referential
simulator), and unvalidated (~30% of the genome flagged, no FDR). The corrected pipeline:

  1. ensure-genome   Extract the real H37Rv genome FASTA from the GenBank file.
  2. extract         (Re)build genome-wide promoters + divergent IGRs if missing.
  3. background      Genomic order-2 Markov background for MEME/FIMO.
  4. operator-model  Knowledge-based operator PWM from REAL operator-region sequence
                     (M. tb IGR + orthologs), via anchored MEME. No typed consensus.
  5. denovo          De novo MEME on the bound IGRs vs a genomic control set, then Tomtom
                     vs the knowledge PWM  ->  Gate G1 (operator recovered, not injected).
  6. scan            Real FIMO genome-wide with the knowledge PWM, genomic background,
                     reporting q-values (BH FDR); score>0 enforced downstream.
  7. conservation    FIMO the operator PWM across mycobacterial orthologs  ->  Gate G5.
  8. validate        Gates G2 (site recovery), G3 (specificity), G4 (regulon enrichment).
  9. rank            Ranked predicted-sites table by q-value/score (tier is a label only).

Requires the MEME Suite (meme, fimo, tomtom, fasta-get-markov). On Apple Silicon:
    CONDA_SUBDIR=osx-64 conda create -y -n meme_x64 -c conda-forge -c bioconda meme
(the wrappers locate this env automatically).
"""

import argparse
import json
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / "scripts"))

import pandas as pd

import mce3r_biology as bio
from build_background import build_background
from build_operator_model import build_operator_model
from conservation import leave_one_lineage_out
from discover_denovo import discover_denovo, run_tomtom
from literature_motif import build_literature_motif
from run_fimo import run_fimo
from utils import ensure_directories, get_project_root, require_meme_tool, setup_logging
from validate import evaluate, write_report

# Keep formatter from dropping these as "unused" before their use sites are scanned.
_USED = (leave_one_lineage_out, run_tomtom, build_literature_motif)

logger = setup_logging("main")

ALL_STEPS = [
    "ensure-genome",
    "extract",
    "background",
    "operator-model",
    "denovo",
    "scan",
    "conservation",
    "validate",
    "rank",
]


def step_ensure_genome(root: Path) -> Path:
    """Extract the real genome FASTA from the GenBank file if it is missing/too small."""
    from download_genome import extract_fasta_from_genbank

    gb = root / "data" / "raw" / "Mtb_H37Rv.gb"
    fasta = root / "data" / "raw" / "Mtb_H37Rv.fasta"
    if not gb.exists():
        raise FileNotFoundError(
            f"GenBank file missing: {gb}. Run download_genome.py first."
        )
    if (not fasta.exists()) or fasta.stat().st_size < 1_000_000:
        logger.info("Genome FASTA missing/truncated — extracting from GenBank...")
        extract_fasta_from_genbank(gb, fasta)
    else:
        logger.info(f"Genome FASTA present ({fasta.stat().st_size / 1e6:.1f} MB).")
    return fasta


def step_extract(root: Path) -> None:
    """Ensure promoters.fasta and divergent_igrs.fasta exist (build from GenBank if not)."""
    from extract_promoters import extract_all_promoters, extract_divergent_igrs
    from Bio import SeqIO

    gb = root / "data" / "raw" / "Mtb_H37Rv.gb"
    promoters = root / "data" / "raw" / "promoters.fasta"
    igrs = root / "data" / "raw" / "divergent_igrs.fasta"
    if promoters.exists() and igrs.exists():
        logger.info(
            "promoters.fasta and divergent_igrs.fasta present — skipping base extraction."
        )
    else:
        logger.info("Extracting genome-wide promoters + divergent IGRs...")
        extract_all_promoters(gb, promoters, promoter_length=300)
        igr_df = extract_divergent_igrs(gb, igrs, min_igr_length=400)
        if not igr_df.empty:
            igr_records = list(SeqIO.parse(str(igrs), "fasta"))
            with open(str(promoters), "a") as f:
                SeqIO.write(igr_records, f, "fasta")
    # Always ensure the mapped operator IGRs are present, even when shorter than the
    # genome-wide IGR length cutoff (the Rv1935c-Rv1936 operator region is < 400 bp).
    ensure_operator_igrs(root, gb, igrs, promoters)


def ensure_operator_igrs(root: Path, gb: Path, igrs: Path, promoters: Path) -> None:
    """
    Make sure each mapped operator IGR (mce3r_biology.OPERATOR_IGRS) is present in
    divergent_igrs.fasta and promoters.fasta, extracting it from the genome regardless of
    length. Operator regions are short and would otherwise be filtered out of the scan.

    Idempotent (audit C11 fix): each file is checked independently, and an operator IGR is
    appended to a file only if it is not already in THAT file. Re-running the pipeline no
    longer duplicates the operator IGRs in promoters.fasta (which previously changed the
    universe size and every gate number on each run).
    """
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from extract_promoters import build_cds_interval_tree
    from utils import compute_gc_content

    in_igrs = (
        {r.id for r in SeqIO.parse(str(igrs), "fasta")} if igrs.exists() else set()
    )
    in_prom = (
        {r.id for r in SeqIO.parse(str(promoters), "fasta")}
        if promoters.exists()
        else set()
    )
    need = {
        igr: igr.replace("IGR_", "").split("_")
        for igr in bio.OPERATOR_IGRS
        if igr not in in_igrs or igr not in in_prom
    }
    if not need:
        logger.info("Both operator IGRs already present in both FASTA files.")
        return

    record = SeqIO.read(str(gb), "genbank")
    genome = str(record.seq).upper()
    by_tag = {c["locus_tag"]: c for c in build_cds_interval_tree(record)}

    built = {}
    for igr_id, (t1, t2) in need.items():
        a, b = by_tag.get(t1), by_tag.get(t2)
        if a is None or b is None:
            logger.warning(f"Cannot extract {igr_id}: missing CDS for {t1} or {t2}")
            continue
        lo, hi = sorted([a, b], key=lambda c: c["start"])
        igr_start, igr_end = lo["end"], hi["start"]
        if igr_end <= igr_start:
            logger.warning(f"{igr_id}: genes overlap (no intergenic gap) — skipping")
            continue
        seq = genome[igr_start:igr_end]
        desc = (
            f"divergent_igr gene1={lo['locus_tag']} gene2={hi['locus_tag']} "
            f"igr_start={igr_start} igr_end={igr_end} "
            f"gc_content={compute_gc_content(seq):.3f} length={len(seq)} operator_region=True"
        )
        built[igr_id] = SeqRecord(Seq(seq), id=igr_id, description=desc)
        logger.info(f"Extracted operator region {igr_id}: {len(seq)} bp")

    add_to_igrs = [r for i, r in built.items() if i not in in_igrs]
    add_to_prom = [r for i, r in built.items() if i not in in_prom]
    if add_to_igrs:
        with open(str(igrs), "a") as f:
            SeqIO.write(add_to_igrs, f, "fasta")
    if add_to_prom:
        with open(str(promoters), "a") as f:
            SeqIO.write(add_to_prom, f, "fasta")


def compute_provenance(root: Path, inputs: list[Path]) -> dict:
    """
    Capture run provenance (audit S6): git SHA, MEME-tool versions, and input checksums.
    Stamped into gates.json so a result can be tied to exact code + data + tools.
    """
    import hashlib
    import subprocess

    def _git(*args: str) -> str:
        try:
            return subprocess.run(
                ["git", "-C", str(root), *args],
                capture_output=True,
                text=True,
                timeout=15,
            ).stdout.strip()
        except Exception:  # noqa: BLE001
            return ""

    def _tool_version(tool: str) -> str:
        try:
            path = require_meme_tool(tool)
            for flag in ("--version", "-version"):
                out = subprocess.run(
                    [path, flag], capture_output=True, text=True, timeout=15
                )
                text = (out.stdout or out.stderr).strip()
                if text and "error" not in text.splitlines()[0].lower():
                    return text.splitlines()[0]
            return "unknown"
        except Exception:  # noqa: BLE001
            return "unavailable"

    def _checksum(p: Path) -> str:
        try:
            h = hashlib.sha256()
            with open(p, "rb") as fh:
                for chunk in iter(lambda: fh.read(1 << 20), b""):
                    h.update(chunk)
            return f"sha256:{h.hexdigest()[:16]}"
        except Exception:  # noqa: BLE001
            return "missing"

    return {
        "git_sha": _git("rev-parse", "HEAD"),
        "git_dirty": bool(_git("status", "--porcelain")),
        "tool_versions": {t: _tool_version(t) for t in ("meme", "fimo", "tomtom")},
        "input_checksums": {p.name: _checksum(p) for p in inputs if p},
    }


def external_corroboration_g1(knowledge_meme: Path, out_dir: Path) -> dict:
    """
    G1 — external corroboration against an INDEPENDENTLY-published operator motif.

    Replaces the circular de-novo-vs-knowledge Tomtom (audit C1: both are MEME runs over the
    same operator DNA + orthologs). Here the knowledge PWM is compared by Tomtom against a
    PWM built from Santangelo 2009 footprinted operator sites (literature, not this pipeline).
    A significant match is genuine external evidence that the PWM is the real operator.
    """
    out_dir = Path(out_dir)
    lit_meme = build_literature_motif(out_dir / "santangelo2009.meme")
    res = run_tomtom(knowledge_meme, lit_meme, out_dir / "tomtom_external")
    return {
        "name": "external_corroboration",
        "reference": bio.SANTANGELO_2009_CITATION,
        "tomtom_best_qvalue": res["best_qvalue"],
        "matched": bool(res["matched"]),
        "note": (
            "Knowledge PWM vs a PWM built from independently-published footprinted operator "
            "sites. Methodologically independent of this pipeline's MEME runs."
        ),
        "pass": bool(res["matched"]),
    }


def main():
    root = get_project_root()
    parser = argparse.ArgumentParser(
        description="Corrected Mce3R operator-scan pipeline"
    )
    parser.add_argument("--steps", nargs="+", choices=ALL_STEPS, default=ALL_STEPS)
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument(
        "--qv-thresh",
        type=float,
        default=0.05,
        help="FIMO reporting q-value (strict gate uses q<0.05)",
    )
    parser.add_argument(
        "--scan-pthresh",
        type=float,
        default=1e-3,
        help="FIMO p-value reporting cutoff for the genome scan/ranking",
    )
    args = parser.parse_args()

    # Fail fast if the MEME Suite is unavailable (no silent simulation fallback).
    for tool in ("meme", "fimo", "tomtom", "fasta-get-markov"):
        require_meme_tool(tool)

    ensure_directories(
        [
            root / d
            for d in [
                "data/raw",
                "data/processed",
                "results/background",
                "results/motifs",
                "results/scans",
            ]
        ]
    )

    raw = root / "data" / "raw"
    bg = root / "results" / "background" / "mtb_genome.bg"
    knowledge_meme = (
        root / "results" / "motifs" / "operator_knowledge" / "operator_knowledge.meme"
    )
    promoters = raw / "promoters.fasta"
    igrs = raw / "divergent_igrs.fasta"
    orthologs = raw / "ortholog_promoters.fasta"
    genome_fimo = root / "results" / "scans" / "fimo.tsv"

    print("\n" + "=" * 66)
    print("  MCE3R OPERATOR SCAN — corrected, validated pipeline")
    print("=" * 66)
    print(f"  Steps: {', '.join(args.steps)}")
    print("=" * 66 + "\n")

    t0 = time.perf_counter()
    g1 = None
    g5 = None

    if "ensure-genome" in args.steps:
        step_ensure_genome(root)
    if "extract" in args.steps:
        step_extract(root)
    if "background" in args.steps:
        build_background(raw / "Mtb_H37Rv.fasta", bg, order=2)
    if "operator-model" in args.steps:
        build_operator_model(
            igrs, orthologs, knowledge_meme.parent, bfile=bg, threads=args.threads
        )
    if "denovo" in args.steps:
        # G1 = external corroboration vs a published motif (non-circular, audit C1 fix).
        g1 = external_corroboration_g1(
            knowledge_meme, root / "results" / "motifs" / "external"
        )
        # The de-novo-vs-knowledge Tomtom is kept as a DESCRIPTIVE internal-consistency
        # readout only (it is circular and never gates the result).
        res = discover_denovo(
            igrs,
            orthologs,
            knowledge_meme,
            root / "results" / "motifs" / "denovo",
            bfile=bg,
            threads=args.threads,
        )
        g1["internal_consistency_denovo_tomtom_q"] = res["g1_best_qvalue"]
        g1["internal_consistency_matched"] = bool(res["g1_matched"])
    if "scan" in args.steps:
        # Permissive p-value reporting so every promoter gets a score (for ranking);
        # the genome-wide BH q<0.05 strict-hit definition is recomputed in validate.py.
        run_fimo(
            knowledge_meme,
            promoters,
            root / "results" / "scans",
            bfile=bg,
            qv_thresh=args.scan_pthresh,
            use_qvalue=False,
        )
    if "conservation" in args.steps:
        # G5 = leave-one-lineage-out (non-circular, audit C7/S5 fix).
        g5 = leave_one_lineage_out(
            orthologs,
            root / "results" / "scans" / "_conservation_lolo",
            threads=args.threads,
        )

    gates = {}
    if "validate" in args.steps:
        gates = evaluate(promoters, genome_fimo, knowledge_meme, q_thresh=0.05)
        if g1 is not None:
            gates["G1"] = g1
        if g5 is not None:
            gates["G5"] = g5
        # Preserve previously-computed G1/G5 when this is a partial re-run that skipped
        # the denovo/conservation steps, so a `--steps validate` run does not silently
        # destroy a complete gate set (audit self-check H-1).
        prior_path = root / "results" / "scans" / "gates.json"
        if prior_path.exists():
            try:
                prior = json.loads(prior_path.read_text())
                for gk in ("G1", "G5"):
                    if gk not in gates and isinstance(prior.get(gk), dict):
                        gates[gk] = prior[gk]
                        logger.info(
                            f"Carried over {gk} from previous gates.json (step not re-run)."
                        )
            except (json.JSONDecodeError, OSError):
                pass
        core = ["G1", "G2", "G3", "G4", "G5"]
        gates["all_gates_pass"] = bool(
            all(gates.get(k, {}).get("pass") for k in core if k in gates)
            and all(k in gates for k in core)
        )
        gates["provenance"] = compute_provenance(
            root, [promoters, orthologs, knowledge_meme, genome_fimo]
        )
        write_report(gates, root / "results" / "scans" / "validation_report.txt")
        (root / "results" / "scans" / "gates.json").write_text(
            json.dumps(gates, indent=2, default=str)
        )

    if "rank" in args.steps:
        write_ranked_sites(
            promoters,
            genome_fimo,
            root / "data" / "processed" / "predicted_mce3r_sites.csv",
        )

    print("\n" + "=" * 66)
    print("  GATE SUMMARY")
    print("=" * 66)
    gate_order = ["G1", "G2", "G3", "G4", "G5"]
    if g1 is not None and "G1" not in gates:
        gates["G1"] = g1
    if g5 is not None and "G5" not in gates:
        gates["G5"] = g5
    for gk in gate_order:
        g = gates.get(gk)
        if g is None:
            continue
        status = "PASS" if g.get("pass") else "FAIL"
        print(f"  {gk}  {g.get('name', ''):<32} {status}")
    if "all_gates_pass" in gates:
        print("-" * 66)
        print(f"  ALL GATES PASS: {gates['all_gates_pass']}")
    print("=" * 66)
    print(f"  total time: {time.perf_counter() - t0:.1f}s")
    print("=" * 66 + "\n")


def write_ranked_sites(promoters_fasta: Path, fimo_tsv: Path, out_csv: Path) -> None:
    """
    Ranked predicted-site table: q-value then score. Priority tier is a DESCRIPTIVE label
    (from mce3r_biology), never a ranking input. Only score>0 rows are kept.
    """
    if not Path(fimo_tsv).exists():
        logger.warning(f"No FIMO output to rank: {fimo_tsv}")
        return
    df = pd.read_csv(fimo_tsv, sep="\t", comment="#")
    df.columns = [c.strip() for c in df.columns]
    for c in ("score", "p-value", "q-value", "start", "stop"):
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce")
    df = df[df["score"] > 0].copy()
    # Sort key per column: q-value ascending (smaller = better), score descending
    # (larger = better). Building the ascending flags per-column avoids the audit S10
    # bug where, with the q-value column absent, score was sorted ascending (worst first).
    sort_spec = {"q-value": True, "score": False}
    sort_cols = [c for c in ("q-value", "score") if c in df.columns]
    df = df.sort_values(
        sort_cols, ascending=[sort_spec[c] for c in sort_cols]
    ).reset_index(drop=True)
    df.insert(0, "rank", df.index + 1)
    df["locus_class"] = df["sequence_name"].map(bio.classify_locus)
    df["gene"] = df["sequence_name"].map(lambda s: bio.GENE_NAMES.get(str(s), ""))
    Path(out_csv).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_csv, index=False)
    logger.info(f"Ranked {len(df)} score>0 sites (by q then score) -> {out_csv}")


if __name__ == "__main__":
    main()
