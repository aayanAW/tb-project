"""
main.py — Mce3R Binding Motif Discovery Pipeline Orchestrator

Runs the complete pipeline for identifying DNA binding sites of Mce3R,
a TetR-family transcription factor in Mycobacterium tuberculosis.

Real-data pipeline (recommended):
    python main.py --steps download extract meme fimo analyze visualize

Synthetic-data pipeline (for testing without internet/MEME Suite):
    python main.py --steps generate meme fimo analyze visualize --simulate

Pipeline steps:
    download   — Download M. tb H37Rv genome from NCBI (NC_000962.3)
    extract    — Extract 300 bp upstream promoter regions for all genes;
                 filter to Mce3R candidate genes for MEME input
    generate   — [Synthetic mode] Generate GC-rich M. tb-like sequences
    meme       — Motif discovery on candidate promoters (MEME or simulation)
    fimo       — Genome-wide motif scan on all promoters (FIMO or simulation)
    analyze    — Statistical analysis: spacing, palindromes, ranked site table
    visualize  — Publication-quality figures

MEME Suite is optional. If not installed, Python-based fallbacks are used.
To install: conda install -c bioconda meme
"""

import argparse
import sys
import time
from pathlib import Path

# Add scripts/ to Python path
project_root = Path(__file__).resolve().parent
sys.path.insert(0, str(project_root / "scripts"))

import numpy as np
import pandas as pd

from utils import (
    check_tool_available,
    ensure_directories,
    get_project_root,
    setup_logging,
)

logger = setup_logging("main")

# Ordered set of all valid pipeline step names
ALL_STEPS = ["download", "extract", "generate", "meme", "fimo", "analyze", "visualize"]

# All directories the pipeline writes into
REQUIRED_DIRS = [
    "data/raw",
    "data/processed",
    "results/motifs",
    "results/scans",
    "results/figures",
]


def check_environment() -> dict:
    """Check availability of all required and optional tools/libraries."""
    availability = {}
    for tool in ["meme", "fimo"]:
        availability[tool] = check_tool_available(tool)
    required_libs = ["Bio", "pandas", "numpy", "matplotlib", "seaborn"]
    for lib in required_libs:
        try:
            __import__(lib)
            availability[lib] = True
        except ImportError:
            availability[lib] = False
    return availability


def print_environment_table(availability: dict) -> None:
    """Print a formatted environment check table."""
    print("\n" + "=" * 60)
    print("  ENVIRONMENT CHECK")
    print("=" * 60)
    print(f"  {'Tool/Library':<20} {'Status':<12} Notes")
    print("  " + "-" * 55)

    notes = {
        "meme": "Optional (simulation fallback available)",
        "fimo": "Optional (simulation fallback available)",
        "Bio": "biopython — required",
        "pandas": "Required",
        "numpy": "Required",
        "matplotlib": "Required for figures",
        "seaborn": "Required for figures",
    }

    all_required_ok = True
    for name, available in availability.items():
        status = "OK" if available else "MISSING"
        note = notes.get(name, "")
        is_required = name not in ["meme", "fimo"]
        if not available and is_required:
            all_required_ok = False
            status = "MISSING *"
        print(f"  {name:<20} {status:<12} {note}")
    print("=" * 60)

    if not all_required_ok:
        print("\n[ERROR] Required Python libraries are missing.")
        print("        Run: pip install -r requirements.txt\n")
        sys.exit(1)

    missing_suite = [t for t in ["meme", "fimo"] if not availability.get(t)]
    if missing_suite:
        print(f"\n[INFO] {', '.join(missing_suite)} not found. Simulation fallbacks will be used.")
        print("       To install: conda install -c bioconda meme\n")
    else:
        print("\n[INFO] MEME Suite found. Real MEME/FIMO will be used.\n")


# ── Step functions ─────────────────────────────────────────────────────────────

def step_download(args, root: Path) -> dict:
    """Step: Download M. tb H37Rv genome from NCBI."""
    from download_genome import download_genome

    logger.info("Step [download]: Fetching M. tuberculosis H37Rv genome from NCBI...")
    result = download_genome(
        accession=args.accession,
        output_dir=root / "data" / "raw",
        email=args.email,
        skip_existing=args.skip_existing,
    )
    return result


def step_extract(args, root: Path) -> dict:
    """Step: Extract promoter regions from GenBank annotation."""
    from extract_promoters import (
        extract_all_promoters,
        extract_candidate_promoters,
        extract_divergent_igrs,
    )

    genbank_path = root / "data" / "raw" / "Mtb_H37Rv.gb"
    all_fasta = root / "data" / "raw" / "promoters.fasta"
    igr_fasta = root / "data" / "raw" / "divergent_igrs.fasta"
    candidate_fasta = root / "data" / "raw" / "mce3r_candidate_promoters.fasta"
    all_metadata = root / "data" / "processed" / "promoter_metadata.csv"
    candidate_metadata = root / "data" / "processed" / "candidate_promoter_metadata.csv"

    if not genbank_path.exists():
        raise FileNotFoundError(
            f"GenBank file not found: {genbank_path}\n"
            "Run with --steps download first."
        )

    # Extract all genome-wide promoters (for FIMO scan)
    if args.skip_existing and all_fasta.exists():
        logger.info(f"Skipping all-promoter extraction: {all_fasta} exists")
        all_df = pd.read_csv(all_metadata) if all_metadata.exists() else pd.DataFrame()
    else:
        logger.info("Step [extract]: Extracting genome-wide promoter regions...")
        all_df = extract_all_promoters(
            genbank_path=genbank_path,
            output_fasta=all_fasta,
            promoter_length=args.promoter_length,
        )

        # Append full divergent IGRs — critical for mce3R-yrbE3A (897 bp) which
        # is split between two 300 bp windows leaving a 297 bp uncovered gap.
        logger.info("Step [extract]: Appending full divergent IGRs (>400 bp)...")
        igr_df = extract_divergent_igrs(
            genbank_path=genbank_path,
            output_fasta=igr_fasta,
            min_igr_length=400,
        )
        if not igr_df.empty:
            from Bio import SeqIO as _SeqIO
            igr_records = list(_SeqIO.parse(str(igr_fasta), "fasta"))
            with open(str(all_fasta), "a") as f_out:
                _SeqIO.write(igr_records, f_out, "fasta")
            # Add IGR rows to metadata
            igr_meta_rows = []
            for _, row in igr_df.iterrows():
                igr_meta_rows.append({
                    "locus_tag": row["igr_id"],
                    "gene": f"{row['gene1_name'] or row['gene1_tag']}/{row['gene2_name'] or row['gene2_tag']}",
                    "product": "divergent intergenic region",
                    "strand": "+",
                    "prom_start": row["igr_start"],
                    "prom_end": row["igr_end"],
                    "prom_length": row["igr_length"],
                    "gc_content": row["gc_content"],
                })
            igr_meta_df = pd.DataFrame(igr_meta_rows)
            all_df = pd.concat([all_df, igr_meta_df], ignore_index=True)
            logger.info(f"Added {len(igr_records)} divergent IGRs to promoters.fasta")

        all_df.to_csv(all_metadata, index=False)

    # Filter to Mce3R candidates (for MEME motif discovery)
    if args.skip_existing and candidate_fasta.exists():
        logger.info(f"Skipping candidate extraction: {candidate_fasta} exists")
        candidate_df = pd.read_csv(candidate_metadata) if candidate_metadata.exists() else pd.DataFrame()
    else:
        logger.info("Step [extract]: Filtering to Mce3R candidate promoters...")
        candidate_df = extract_candidate_promoters(
            all_promoters_fasta=all_fasta,
            candidate_fasta=candidate_fasta,
        )
        candidate_df.to_csv(candidate_metadata, index=False)

    return {
        "all_fasta": all_fasta,
        "candidate_fasta": candidate_fasta,
        "n_all": len(all_df),
        "n_candidates": len(candidate_df),
    }


def step_generate(args, root: Path, rng: np.random.Generator) -> dict:
    """Step: Generate synthetic M. tb-like promoter sequences (testing mode)."""
    from generate_sequences import generate_synthetic_dataset

    fasta_path = root / "data" / "raw" / "sequences.fasta"
    metadata_path = root / "data" / "processed" / "sequence_metadata.csv"

    if args.skip_existing and fasta_path.exists():
        logger.info(f"Skipping generate step: {fasta_path} exists")
        return {"fasta_path": fasta_path, "metadata_path": metadata_path}

    logger.info("Step [generate]: Generating synthetic M. tb promoter sequences...")
    df = generate_synthetic_dataset(
        n_sequences=args.n_sequences,
        n_with_motif=args.n_with_motif,
        seed=args.seed,
        output_path=fasta_path,
    )
    metadata_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(metadata_path, index=False)
    return {"fasta_path": fasta_path, "metadata_path": metadata_path}


def _resolve_meme_input(args, root: Path) -> Path:
    """
    Determine which FASTA to use as MEME input.

    Priority:
    1. Mce3R candidate promoters (real data, from extract step)
    2. Synthetic sequences (from generate step)
    """
    candidate_fasta = root / "data" / "raw" / "mce3r_candidate_promoters.fasta"
    synthetic_fasta = root / "data" / "raw" / "sequences.fasta"

    if candidate_fasta.exists():
        logger.info(f"MEME input: candidate promoters ({candidate_fasta.name})")
        return candidate_fasta
    elif synthetic_fasta.exists():
        logger.info(f"MEME input: synthetic sequences ({synthetic_fasta.name})")
        return synthetic_fasta
    else:
        raise FileNotFoundError(
            "No FASTA found for MEME. Run 'download + extract' (real data) "
            "or 'generate' (synthetic) first."
        )


def _resolve_fimo_input(args, root: Path) -> Path:
    """
    Determine which FASTA to scan with FIMO.

    Priority:
    1. All genome-wide promoters (real data, from extract step)
    2. Synthetic sequences (from generate step)
    """
    all_fasta = root / "data" / "raw" / "promoters.fasta"
    synthetic_fasta = root / "data" / "raw" / "sequences.fasta"

    if all_fasta.exists():
        logger.info(f"FIMO scan target: all promoters ({all_fasta.name})")
        return all_fasta
    elif synthetic_fasta.exists():
        logger.info(f"FIMO scan target: synthetic sequences ({synthetic_fasta.name})")
        return synthetic_fasta
    else:
        raise FileNotFoundError(
            "No FASTA found for FIMO. Run 'download + extract' or 'generate' first."
        )


def step_meme(args, root: Path, rng: np.random.Generator) -> dict:
    """Step: Run MEME motif discovery on candidate promoters."""
    from run_meme import run_meme

    fasta_path = _resolve_meme_input(args, root)
    output_dir = root / "results" / "motifs"
    meme_txt = output_dir / "meme.txt"

    if args.skip_existing and meme_txt.exists():
        logger.info(f"Skipping MEME: {meme_txt} exists")
        return {"meme_txt": meme_txt, "fasta_used": fasta_path}

    logger.info("Step [meme]: Running MEME motif discovery...")
    meme_txt = run_meme(
        fasta_path=fasta_path,
        output_dir=output_dir,
        simulate=args.simulate,
        threads=args.threads,
        nmotifs=5,
        minw=20,  # ~25 bp Mce3R half-sites (Panagoda et al. 2024)
        maxw=30,
        rng=rng,
    )
    return {"meme_txt": meme_txt, "fasta_used": fasta_path}


def step_fimo(args, root: Path, meme_txt: Path, rng: np.random.Generator) -> dict:
    """Step: Scan all promoter regions with FIMO (palindromic + asymmetric models)."""
    from run_fimo import run_fimo, run_fimo_asymmetric

    fasta_path = _resolve_fimo_input(args, root)
    output_dir = root / "results" / "scans"
    motifs_dir = root / "results" / "motifs"
    fimo_tsv = output_dir / "fimo.tsv"

    if args.skip_existing and fimo_tsv.exists():
        logger.info(f"Skipping FIMO: {fimo_tsv} exists")
        return {"fimo_tsv": fimo_tsv}

    logger.info("Step [fimo]: Scanning promoters for Mce3R binding sites (palindromic model)...")
    fimo_tsv = run_fimo(
        motif_file=meme_txt,
        fasta_path=fasta_path,
        output_dir=output_dir,
        pvalue_threshold=args.pvalue_threshold,
        simulate=args.simulate,
        rng=rng,
    )

    # Also run the asymmetric model scan (always Python simulation)
    asym_tsv = run_fimo_asymmetric(
        fasta_path=fasta_path,
        output_dir=output_dir,
        motifs_dir=motifs_dir,
        pvalue_threshold=args.pvalue_threshold,
        rng=rng,
    )

    return {"fimo_tsv": fimo_tsv, "fimo_asymmetric_tsv": asym_tsv}


def step_analyze(args, root: Path) -> dict:
    """Step: Analyze FIMO results — ranking, spacing, palindromes, model comparison."""
    from analyze_results import (
        annotate_hits_with_sequence_metadata,
        compute_inter_site_spacing,
        compute_motif_statistics,
        detect_direct_repeats,
        detect_palindromic_sites,
        load_fimo_results,
        rank_binding_sites,
    )
    import json

    logger.info("Step [analyze]: Analyzing Mce3R motif hit statistics...")

    fimo_tsv = root / "results" / "scans" / "fimo.tsv"
    fimo_asym_tsv = root / "results" / "scans" / "fimo_asymmetric.tsv"
    processed_dir = root / "data" / "processed"
    scans_dir = root / "results" / "scans"

    # Load palindromic model hits (primary scan)
    hits_df = load_fimo_results(fimo_tsv)

    # Merge asymmetric model hits — add sites found ONLY by the asymmetric model.
    # This is critical: the mce3R operator (Rv1963c) is detected at p=1e-4 by the
    # asymmetric model but NOT by the palindromic model (score gap ~4 bits).
    if fimo_asym_tsv.exists():
        asym_df = load_fimo_results(fimo_asym_tsv)
        if not asym_df.empty:
            # Identify sites that appear only in the asymmetric scan:
            # a site is "already covered" if the same sequence has a hit within
            # ±30 bp in the palindromic scan. Otherwise add it.
            if hits_df.empty:
                asym_df["source_model"] = "asymmetric_only"
                hits_df = asym_df
            else:
                hits_df["source_model"] = "palindromic"
                asym_df["source_model"] = "asymmetric_only"
                # Build a set of (seq_name, approx_position) for palindromic hits
                pal_positions = set()
                for _, row in hits_df.iterrows():
                    for offset in range(-30, 31):
                        pal_positions.add((str(row["sequence_name"]), int(row["start"]) + offset))
                # Keep asymmetric hits that don't overlap any palindromic hit
                asym_unique_rows = []
                for _, row in asym_df.iterrows():
                    key = (str(row["sequence_name"]), int(row["start"]))
                    if key not in pal_positions:
                        asym_unique_rows.append(row)
                if asym_unique_rows:
                    asym_unique_df = pd.DataFrame(asym_unique_rows)
                    hits_df = pd.concat([hits_df, asym_unique_df], ignore_index=True)
                    n_added = len(asym_unique_rows)
                    logger.info(
                        f"Merged {n_added} asymmetric-only hits into main hit table "
                        f"(total: {len(hits_df)} hits)"
                    )

    # Load promoter metadata (real or synthetic)
    metadata_candidates = [
        processed_dir / "promoter_metadata.csv",        # real data
        processed_dir / "sequence_metadata.csv",        # synthetic data
    ]
    sequences_df = pd.DataFrame()
    for meta_path in metadata_candidates:
        if meta_path.exists():
            sequences_df = pd.read_csv(meta_path)
            logger.info(f"Loaded metadata from {meta_path.name}: {len(sequences_df)} rows")
            break

    stats = compute_motif_statistics(hits_df)
    spacing_df = compute_inter_site_spacing(hits_df)
    hits_df = detect_palindromic_sites(hits_df)
    hits_df = detect_direct_repeats(hits_df)

    # Apply asymmetric model comparison (2024 ACS Chemical Biology)
    try:
        from motif_models import (
            apply_model_comparison,
            classify_site_architectures,
            generate_model_comparison_report,
            generate_architecture_comparison_report,
        )
        hits_df = apply_model_comparison(hits_df)
        hits_df = classify_site_architectures(hits_df)
        logger.info("Applied asymmetric model comparison (2024 ACS Chem. Biol.)")
    except ImportError:
        generate_model_comparison_report = None
        generate_architecture_comparison_report = None
        logger.warning("motif_models.py not found — skipping model comparison")

    if not sequences_df.empty:
        annotated_df = annotate_hits_with_sequence_metadata(hits_df, sequences_df)
    else:
        annotated_df = hits_df

    # Ranked binding sites table (uses composite score if model columns are present)
    ranked_df = rank_binding_sites(annotated_df)

    # Save outputs
    annotated_df.to_csv(processed_dir / "annotated_hits.csv", index=False)
    spacing_df.to_csv(processed_dir / "spacing_analysis.csv", index=False)
    ranked_df.to_csv(processed_dir / "predicted_mce3r_sites.csv", index=False)
    with open(scans_dir / "summary_statistics.json", "w") as f:
        json.dump(stats, f, indent=2, default=str)

    # Generate plain-text model comparison report + architecture comparison report
    report_path = None
    arch_report_path = None
    if "generate_model_comparison_report" in dir() and generate_model_comparison_report is not None:
        try:
            report_path = scans_dir / "mce3r_prediction_report.txt"
            generate_model_comparison_report(annotated_df, ranked_df, report_path)
        except Exception as e:
            logger.warning(f"Report generation failed: {e}")
            report_path = None
    if "generate_architecture_comparison_report" in dir() and generate_architecture_comparison_report is not None:
        try:
            arch_report_path = scans_dir / "architecture_comparison_report.txt"
            generate_architecture_comparison_report(ranked_df, arch_report_path)
        except Exception as e:
            logger.warning(f"Architecture comparison report failed: {e}")
            arch_report_path = None

    logger.info(f"Ranked {len(ranked_df)} predicted Mce3R binding sites")

    # Tier1 debug report — always generated, explains why Tier1 sites do/don't appear
    try:
        from analyze_results import generate_tier1_debug_report
        tier1_debug_path = scans_dir / "tier1_diagnostic_report.txt"
        promoters_fasta = root / "data" / "raw" / "promoters.fasta"
        motif_files = [
            root / "results" / "motifs" / "meme.txt",
            root / "results" / "motifs" / "meme_asymmetric.txt",
        ]
        generate_tier1_debug_report(ranked_df, promoters_fasta, motif_files, tier1_debug_path)
    except Exception as e:
        logger.warning(f"Tier1 debug report failed: {e}")
        tier1_debug_path = None

    # Print tier summary
    if not ranked_df.empty and "priority_tier" in ranked_df.columns:
        from analyze_results import TIER_LABELS
        score_col = "composite_score" if "composite_score" in ranked_df.columns else "score"
        print("\n  Priority Tier Summary (ranked by tier then score):")
        for tier_order, tier_label in sorted(TIER_LABELS.items()):
            tier_rows = ranked_df[ranked_df["priority_tier"] == tier_label]
            count = len(tier_rows)
            if count == 0:
                print(f"    {tier_label:<30}: 0 sites")
                continue
            top = tier_rows.iloc[0]
            print(
                f"    {tier_label:<30}: {count:>3} sites  "
                f"[top: {str(top['sequence_name']):<14} {score_col}={top[score_col]:.3f}  "
                f"strand={top['strand']}]"
            )

    return {
        "hits_df": annotated_df,
        "sequences_df": sequences_df,
        "spacing_df": spacing_df,
        "ranked_df": ranked_df,
        "stats": stats,
        "report_path": report_path,
        "arch_report_path": arch_report_path,
        "tier1_debug_path": tier1_debug_path,
    }


def step_visualize(
    args,
    root: Path,
    hits_df: pd.DataFrame,
    sequences_df: pd.DataFrame,
    spacing_df: pd.DataFrame,
) -> dict:
    """Step: Generate all publication-quality figures."""
    from visualize import generate_all_figures

    logger.info("Step [visualize]: Generating figures...")
    figures_dir = root / "results" / "figures"
    generated = generate_all_figures(hits_df, sequences_df, spacing_df, figures_dir)
    return {"figures": generated}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Mce3R Binding Motif Discovery Pipeline\n"
            "Identifies DNA binding sites of the TetR-family TF Mce3R\n"
            "in Mycobacterium tuberculosis promoter regions.\n\n"
            "Real data:    python main.py --steps download extract meme fimo analyze visualize\n"
            "Synthetic:    python main.py --steps generate meme fimo analyze visualize --simulate"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--steps",
        nargs="+",
        choices=ALL_STEPS,
        default=["download", "extract", "meme", "fimo", "analyze", "visualize"],
        help=f"Pipeline steps to run. Choices: {ALL_STEPS}",
    )
    # NCBI download options
    parser.add_argument(
        "--accession",
        default="NC_000962.3",
        help="NCBI genome accession (default: NC_000962.3 = M. tb H37Rv)",
    )
    parser.add_argument(
        "--email",
        default="user@example.com",
        help="Email for NCBI Entrez API (required by NCBI policy)",
    )
    # Extraction options
    parser.add_argument(
        "--promoter-length",
        type=int,
        default=300,
        help="Upstream region length in bp (default: 300)",
    )
    # Synthetic data options (for --steps generate ...)
    parser.add_argument(
        "--n-sequences",
        type=int,
        default=25,
        help="[Synthetic mode] Number of sequences to generate (default: 25)",
    )
    parser.add_argument(
        "--n-with-motif",
        type=int,
        default=14,
        help="[Synthetic mode] Sequences with embedded motif (default: 14)",
    )
    # MEME / FIMO options
    parser.add_argument("--threads", type=int, default=4, help="CPU threads for MEME")
    parser.add_argument(
        "--pvalue-threshold",
        type=float,
        default=1e-4,
        help="FIMO p-value cutoff for reporting hits (default: 1e-4)",
    )
    parser.add_argument(
        "--simulate",
        action="store_true",
        help="Force simulation mode (skip real MEME/FIMO even if installed)",
    )
    parser.add_argument(
        "--skip-existing",
        action="store_true",
        help="Skip steps where output files already exist",
    )
    parser.add_argument("--seed", type=int, default=42, help="Random seed (default: 42)")
    return parser.parse_args()


def run_pipeline(args: argparse.Namespace) -> None:
    """Execute the Mce3R motif discovery pipeline."""
    root = get_project_root()
    rng = np.random.default_rng(args.seed)

    # Determine mode: real data or synthetic
    using_real_data = any(s in args.steps for s in ["download", "extract"])

    print("\n" + "=" * 65)
    print("  MCE3R BINDING MOTIF DISCOVERY PIPELINE")
    print("  Mycobacterium tuberculosis TetR-family TF Analysis")
    print("=" * 65)
    print(f"  Project root  : {root}")
    print(f"  Steps         : {', '.join(args.steps)}")
    print(f"  Mode          : {'Real M. tb genome data' if using_real_data else 'Synthetic sequences'}")
    print(f"  Simulate      : {args.simulate}")
    if using_real_data:
        print(f"  Accession     : {args.accession}")
        print(f"  Promoter len  : {args.promoter_length} bp")
    print(f"  p-value cutoff: {args.pvalue_threshold}")
    print("=" * 65 + "\n")

    # Environment check
    availability = check_environment()
    print_environment_table(availability)

    # Create all output directories
    ensure_directories([root / d for d in REQUIRED_DIRS])

    step_results = {}
    step_timing = {}
    pipeline_start = time.perf_counter()

    meme_txt = root / "results" / "motifs" / "meme.txt"

    # ── download ────────────────────────────────────────────────────────────
    if "download" in args.steps:
        t0 = time.perf_counter()
        try:
            result = step_download(args, root)
            step_results["download"] = ("OK", result)
        except Exception as e:
            logger.error(f"Download step failed: {e}")
            step_results["download"] = ("FAILED", str(e))
            sys.exit(1)
        step_timing["download"] = time.perf_counter() - t0

    # ── extract ─────────────────────────────────────────────────────────────
    if "extract" in args.steps:
        t0 = time.perf_counter()
        try:
            result = step_extract(args, root)
            step_results["extract"] = ("OK", result)
        except Exception as e:
            logger.error(f"Extract step failed: {e}")
            step_results["extract"] = ("FAILED", str(e))
            sys.exit(1)
        step_timing["extract"] = time.perf_counter() - t0

    # ── generate (synthetic mode) ────────────────────────────────────────────
    if "generate" in args.steps:
        t0 = time.perf_counter()
        try:
            result = step_generate(args, root, rng)
            step_results["generate"] = ("OK", result)
        except Exception as e:
            logger.error(f"Generate step failed: {e}")
            step_results["generate"] = ("FAILED", str(e))
            sys.exit(1)
        step_timing["generate"] = time.perf_counter() - t0

    # ── meme ─────────────────────────────────────────────────────────────────
    if "meme" in args.steps:
        t0 = time.perf_counter()
        try:
            result = step_meme(args, root, rng)
            step_results["meme"] = ("OK", result)
            meme_txt = result["meme_txt"]
        except Exception as e:
            logger.error(f"MEME step failed: {e}")
            step_results["meme"] = ("FAILED", str(e))
            sys.exit(1)
        step_timing["meme"] = time.perf_counter() - t0

    # ── fimo ─────────────────────────────────────────────────────────────────
    if "fimo" in args.steps:
        t0 = time.perf_counter()
        try:
            result = step_fimo(args, root, meme_txt, rng)
            step_results["fimo"] = ("OK", result)
        except Exception as e:
            logger.error(f"FIMO step failed: {e}")
            step_results["fimo"] = ("FAILED", str(e))
            sys.exit(1)
        step_timing["fimo"] = time.perf_counter() - t0

    # ── analyze ───────────────────────────────────────────────────────────────
    hits_df = pd.DataFrame()
    sequences_df = pd.DataFrame()
    spacing_df = pd.DataFrame()
    report_path = None
    arch_report_path = None
    tier1_debug_path = None
    if "analyze" in args.steps:
        t0 = time.perf_counter()
        try:
            result = step_analyze(args, root)
            step_results["analyze"] = ("OK", result)
            hits_df = result["hits_df"]
            sequences_df = result["sequences_df"]
            spacing_df = result["spacing_df"]
            report_path = result.get("report_path")
            arch_report_path = result.get("arch_report_path")
            tier1_debug_path = result.get("tier1_debug_path")
        except Exception as e:
            logger.error(f"Analyze step failed: {e}")
            step_results["analyze"] = ("FAILED", str(e))
        step_timing["analyze"] = time.perf_counter() - t0

    # ── visualize ─────────────────────────────────────────────────────────────
    if "visualize" in args.steps:
        t0 = time.perf_counter()
        # Load from disk if analyze was skipped
        if hits_df.empty:
            p = root / "data" / "processed" / "annotated_hits.csv"
            hits_df = pd.read_csv(p) if p.exists() else pd.DataFrame()
        if sequences_df.empty:
            for mp in [
                root / "data" / "processed" / "promoter_metadata.csv",
                root / "data" / "processed" / "sequence_metadata.csv",
            ]:
                if mp.exists():
                    sequences_df = pd.read_csv(mp)
                    break
        if spacing_df.empty:
            p = root / "data" / "processed" / "spacing_analysis.csv"
            spacing_df = pd.read_csv(p) if p.exists() else pd.DataFrame()

        try:
            result = step_visualize(args, root, hits_df, sequences_df, spacing_df)
            step_results["visualize"] = ("OK", result)
        except Exception as e:
            logger.error(f"Visualize step failed: {e}")
            step_results["visualize"] = ("FAILED", str(e))
        step_timing["visualize"] = time.perf_counter() - t0

    # ── Final Summary ──────────────────────────────────────────────────────────
    total_time = time.perf_counter() - pipeline_start

    output_map = {
        "download":  str(root / "data" / "raw" / "Mtb_H37Rv.gb"),
        "extract":   str(root / "data" / "raw" / "mce3r_candidate_promoters.fasta"),
        "generate":  str(root / "data" / "raw" / "sequences.fasta"),
        "meme":      str(root / "results" / "motifs" / "meme.txt"),
        "fimo":      str(root / "results" / "scans" / "fimo.tsv"),
        "analyze":   str(root / "data" / "processed" / "predicted_mce3r_sites.csv"),
        "visualize": str(root / "results" / "figures"),
    }

    print("\n" + "=" * 65)
    print("  PIPELINE SUMMARY")
    print("=" * 65)
    print(f"  {'Step':<12} {'Status':<10} {'Time':>8}   Output")
    print("  " + "-" * 60)
    for step in args.steps:
        status, _ = step_results.get(step, ("SKIPPED", ""))
        duration = step_timing.get(step, 0)
        out = output_map.get(step, "")
        print(f"  {step:<12} {status:<10} {duration:>6.1f}s   {out}")
    print("=" * 65)
    print(f"  Total pipeline time: {total_time:.1f}s")
    print("=" * 65)

    print(f"\nKey output files:")
    if (root / "data" / "raw" / "promoters.fasta").exists():
        print(f"  All promoters    : {root / 'data' / 'raw' / 'promoters.fasta'}")
    if (root / "data" / "raw" / "mce3r_candidate_promoters.fasta").exists():
        print(f"  Candidates       : {root / 'data' / 'raw' / 'mce3r_candidate_promoters.fasta'}")
    print(f"  Motif (palindromic): {root / 'results' / 'motifs' / 'meme.txt'}")
    asym_motif = root / "results" / "motifs" / "meme_asymmetric.txt"
    if asym_motif.exists():
        print(f"  Motif (asymmetric): {asym_motif}")
    print(f"  FIMO hits        : {root / 'results' / 'scans' / 'fimo.tsv'}")
    asym_tsv = root / "results" / "scans" / "fimo_asymmetric.tsv"
    if asym_tsv.exists():
        print(f"  FIMO (asymmetric): {asym_tsv}")
    print(f"  Ranked sites     : {root / 'data' / 'processed' / 'predicted_mce3r_sites.csv'}")
    print(f"  Figures          : {root / 'results' / 'figures'}/")
    rpt = report_path or (root / "results" / "scans" / "mce3r_prediction_report.txt")
    if Path(rpt).exists():
        print(f"  Prediction report: {rpt}")
    arch_rpt = arch_report_path or (root / "results" / "scans" / "architecture_comparison_report.txt")
    if Path(arch_rpt).exists():
        print(f"  Architecture cmp : {arch_rpt}")
    t1_rpt = tier1_debug_path or (root / "results" / "scans" / "tier1_diagnostic_report.txt")
    if Path(t1_rpt).exists():
        print(f"  Tier1 diagnostic : {t1_rpt}")
    print()


if __name__ == "__main__":
    args = parse_args()
    run_pipeline(args)
