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
from discover_denovo import discover_denovo
from run_fimo import run_fimo
from utils import ensure_directories, get_project_root, require_meme_tool, setup_logging
from validate import evaluate, write_report

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
    """
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from extract_promoters import build_cds_interval_tree
    from utils import compute_gc_content

    present = (
        {r.id for r in SeqIO.parse(str(igrs), "fasta")} if igrs.exists() else set()
    )
    pairs = {igr: igr.replace("IGR_", "").split("_") for igr in bio.OPERATOR_IGRS}
    missing = {igr: gp for igr, gp in pairs.items() if igr not in present}
    if not missing:
        logger.info("Both operator IGRs already present.")
        return

    record = SeqIO.read(str(gb), "genbank")
    genome = str(record.seq).upper()
    by_tag = {c["locus_tag"]: c for c in build_cds_interval_tree(record)}

    new_records = []
    for igr_id, (t1, t2) in missing.items():
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
        new_records.append(SeqRecord(Seq(seq), id=igr_id, description=desc))
        logger.info(f"Extracted operator region {igr_id}: {len(seq)} bp")

    if new_records:
        with open(str(igrs), "a") as f:
            SeqIO.write(new_records, f, "fasta")
        with open(str(promoters), "a") as f:
            SeqIO.write(new_records, f, "fasta")


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
    ortholog_fimo = root / "results" / "scans" / "fimo_orthologs.tsv"

    print("\n" + "=" * 66)
    print("  MCE3R OPERATOR SCAN — corrected, validated pipeline")
    print("=" * 66)
    print(f"  Steps: {', '.join(args.steps)}")
    print("=" * 66 + "\n")

    t0 = time.perf_counter()
    g1 = None

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
        res = discover_denovo(
            igrs,
            orthologs,
            knowledge_meme,
            root / "results" / "motifs" / "denovo",
            bfile=bg,
            threads=args.threads,
        )
        g1 = {
            "name": "operator_recovered_not_injected",
            "matched": res["g1_matched"],
            "best_qvalue": res["g1_best_qvalue"],
            "pass": bool(res["g1_matched"]),
        }
    if "scan" in args.steps:
        # Permissive p-value reporting so every promoter gets a score (for ranking/AUPRC);
        # the strict q<0.05 hit definition is applied in validate.py.
        run_fimo(
            knowledge_meme,
            promoters,
            root / "results" / "scans",
            bfile=bg,
            qv_thresh=args.scan_pthresh,
            use_qvalue=False,
        )
    if "conservation" in args.steps:
        run_fimo(
            knowledge_meme,
            orthologs,
            root / "results" / "scans" / "_orthologs",
            bfile=bg,
            qv_thresh=args.scan_pthresh,
            use_qvalue=False,
        )
        src = root / "results" / "scans" / "_orthologs" / "fimo.tsv"
        if src.exists():
            ortholog_fimo.write_text(src.read_text())

    gates = {}
    if "validate" in args.steps:
        gates = evaluate(
            promoters,
            genome_fimo,
            ortholog_fimo if ortholog_fimo.exists() else None,
            q_thresh=0.05,
        )
        if g1 is not None:
            gates["G1"] = g1
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
    for gk in gate_order:
        g = gates.get(gk)
        if g is None:
            continue
        status = "PASS" if g.get("pass") else "FAIL"
        print(f"  {gk}  {g.get('name', ''):<32} {status}")
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
    sort_cols = [c for c in ("q-value", "score") if c in df.columns]
    df = df.sort_values(
        sort_cols, ascending=[True, False][: len(sort_cols)]
    ).reset_index(drop=True)
    df.insert(0, "rank", df.index + 1)
    df["locus_class"] = df["sequence_name"].map(bio.classify_locus)
    df["gene"] = df["sequence_name"].map(lambda s: bio.GENE_NAMES.get(str(s), ""))
    Path(out_csv).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_csv, index=False)
    logger.info(f"Ranked {len(df)} score>0 sites (by q then score) -> {out_csv}")


if __name__ == "__main__":
    main()
