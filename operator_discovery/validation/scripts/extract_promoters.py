"""
extract_promoters.py — Extract promoter regions from M. tuberculosis genome annotation.

Parses the GenBank file (NC_000962.3) to identify upstream promoter regions for
every annotated gene, then filters to a curated set of Mce3R candidate genes.

Promoter extraction logic:
    For each CDS (coding sequence) on the + strand:
        promoter = [gene_start - 300 bp, gene_start)
    For each CDS on the − strand:
        promoter = [gene_end, gene_end + 300 bp)   (then reverse-complemented)

Overlap handling:
    If the upstream region of gene B overlaps with an adjacent gene A's CDS,
    the promoter is trimmed to begin at the end of gene A's CDS.
    If the trimmed region is < 50 bp, the promoter is discarded as uninformative.

Mce3R candidate genes (literature-curated):
    Primary targets — mce3 operon and autoregulation (high confidence):
        Rv1963c (mce3R): The TF itself; TetR-family regulators often autoregulate
        Rv1964 (lprK) to Rv1978: mce3 locus genes encoding mammalian cell entry proteins

    Extended candidates — lipid/cholesterol metabolism (biologically plausible):
        Rv3575c (kshA): 3-ketosteroid-9α-hydroxylase, cholesterol catabolism
        Rv3537 (kshB): 3-ketosteroid-9α-hydroxylase reductase
        Rv3540c–Rv3545c: igr locus, inner membrane lipid transport
        Rv0762c (fadD3): Fatty-acyl-CoA synthetase
        Rv1176c (fadD4): Fatty-acyl-CoA synthetase

    These candidates are used as MEME input to discover the Mce3R binding motif.
    The motif is then used by FIMO to scan all ~4,000 M. tb promoters genome-wide.

References:
    - Casali et al. (2006) Microbiology: Mce3R autoregulation and mce3 operon repression
    - Forrellad et al. (2013) Virulence: Mce proteins in M. tb virulence
    - Pandey & Sassetti (2008) PNAS: mce operons and cholesterol import

Usage:
    python extract_promoters.py [--genbank PATH] [--promoter-length N]
"""

import argparse
import sys
from pathlib import Path

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

sys.path.insert(0, str(Path(__file__).resolve().parent))
from utils import compute_gc_content, get_project_root, setup_logging

logger = setup_logging(__name__)

# ── Candidate gene locus tags ──────────────────────────────────────────────────
# Curated from published literature on Mce3R regulation in M. tuberculosis H37Rv.
# Format: {locus_tag: (gene_name, confidence, biological_rationale)}

MCE3R_CANDIDATE_GENES = {
    # Primary targets: mce3 locus (high confidence)
    "Rv1963c": ("mce3R", "primary", "Mce3R TF itself; TetR-family autoregulation"),
    "Rv1964":  ("lprK",  "primary", "mce3 locus; lipoprotein, surface-exposed"),
    "Rv1965":  ("Rv1965","primary", "mce3 locus; conserved hypothetical"),
    "Rv1966":  ("Rv1966","primary", "mce3 locus; conserved hypothetical"),
    "Rv1967":  ("Rv1967","primary", "mce3 locus; conserved hypothetical"),
    "Rv1968":  ("Rv1968","primary", "mce3 locus; conserved hypothetical"),
    "Rv1969":  ("Rv1969","primary", "mce3 locus; conserved hypothetical"),
    "Rv1970":  ("mce3A", "primary", "mce3 operon; mammalian cell entry protein A"),
    "Rv1971":  ("mce3B", "primary", "mce3 operon; mammalian cell entry protein B"),
    "Rv1972":  ("mce3C", "primary", "mce3 operon; mammalian cell entry protein C"),
    "Rv1973":  ("mce3D", "primary", "mce3 operon; mammalian cell entry protein D"),
    "Rv1974":  ("mce3E", "primary", "mce3 operon; mammalian cell entry protein E"),
    "Rv1975":  ("mce3F", "primary", "mce3 operon; mammalian cell entry protein F"),
    "Rv1976c": ("lprN",  "primary", "mce3 operon; exported lipoprotein N"),
    "Rv1977":  ("Rv1977","primary", "mce3 locus downstream gene"),
    "Rv1978":  ("Rv1978","primary", "mce3 locus downstream gene"),
    # Extended candidates: lipid/cholesterol metabolism (plausible)
    "Rv3575c": ("kshA",  "extended", "3-ketosteroid-9α-hydroxylase; cholesterol catabolism"),
    "Rv3537":  ("kshB",  "extended", "3-ketosteroid-9α-hydroxylase reductase"),
    "Rv3540c": ("ltp2",  "extended", "igr locus; lipid transfer during cholesterol import"),
    "Rv3541c": ("Rv3541c","extended","igr locus; cholesterol import"),
    "Rv3542c": ("Rv3542c","extended","igr locus; cholesterol import"),
    "Rv3543c": ("Rv3543c","extended","igr locus; cholesterol import"),
    "Rv3544c": ("Rv3544c","extended","igr locus; cholesterol import"),
    "Rv3545c": ("cyp125A1","extended","Cytochrome P450; cholesterol side-chain oxidation"),
    "Rv0762c": ("fadD3",  "extended", "Fatty-acyl-CoA synthetase; lipid activation"),
    "Rv1176c": ("fadD4",  "extended", "Fatty-acyl-CoA synthetase; lipid activation"),
}

# Minimum promoter length (bp) — shorter regions are too small for MEME
MIN_PROMOTER_LENGTH = 50


def build_cds_interval_tree(record) -> list:
    """
    Build a sorted list of all CDS intervals from the GenBank record.

    Used to check for overlaps when extracting promoter regions.

    Args:
        record: BioPython SeqRecord from GenBank parsing

    Returns:
        List of dicts: {start, end, strand, locus_tag, gene, product}
        Sorted by start position (ascending).
    """
    cds_list = []
    for feature in record.features:
        if feature.type not in ("CDS",):
            continue
        start = int(feature.location.start)  # 0-based
        end = int(feature.location.end)      # exclusive
        strand = feature.location.strand
        locus_tag = feature.qualifiers.get("locus_tag", ["unknown"])[0]
        gene = feature.qualifiers.get("gene", [""])[0]
        product = feature.qualifiers.get("product", ["hypothetical protein"])[0]
        cds_list.append(
            {
                "start": start,
                "end": end,
                "strand": strand,
                "locus_tag": locus_tag,
                "gene": gene,
                "product": product,
            }
        )
    cds_list.sort(key=lambda x: x["start"])
    return cds_list


def find_trimmed_promoter_bounds(
    cds_start: int,
    cds_end: int,
    strand: int,
    cds_list: list,
    genome_len: int,
    promoter_length: int = 300,
) -> tuple:
    """
    Determine the valid (non-overlapping) promoter region bounds for a gene.

    For a + strand gene starting at `cds_start`:
        - Raw upstream region: [cds_start - promoter_length, cds_start)
        - Trimmed: begins immediately after the closest upstream CDS that
          would otherwise overlap, if any.

    For a − strand gene ending at `cds_end`:
        - Raw downstream region (which is upstream in genomic terms):
          [cds_end, cds_end + promoter_length)
        - Trimmed: ends at the start of the closest downstream CDS.

    Args:
        cds_start: 0-based start of the gene's CDS
        cds_end: 0-based exclusive end of the gene's CDS
        strand: +1 for plus strand, -1 for minus strand
        cds_list: Sorted list of all CDS intervals
        genome_len: Total genome length in bp
        promoter_length: Desired upstream region length

    Returns:
        Tuple (prom_start, prom_end) in 0-based half-open coordinates,
        or (None, None) if the valid region is too short.
    """
    if strand == 1:
        # Plus strand: promoter is upstream (lower coordinates)
        raw_start = max(0, cds_start - promoter_length)
        raw_end = cds_start

        # Find the closest CDS that ends within our upstream region
        # (i.e., an adjacent gene whose CDS overlaps our promoter window)
        max_upstream_end = raw_start
        for other in cds_list:
            if other["start"] == cds_start:
                continue  # Skip self
            # An adjacent CDS ends within our upstream window
            if raw_start <= other["end"] <= raw_end:
                max_upstream_end = max(max_upstream_end, other["end"])
            # An adjacent CDS completely spans our upstream window (rare but possible)
            elif other["start"] <= raw_start and other["end"] >= raw_end:
                return None, None  # Our entire promoter is inside another CDS

        prom_start = max_upstream_end
        prom_end = raw_end

    else:
        # Minus strand: promoter is downstream (higher coordinates) of CDS end
        raw_start = cds_end
        raw_end = min(genome_len, cds_end + promoter_length)

        # Find the closest CDS that starts within our downstream window
        min_downstream_start = raw_end
        for other in cds_list:
            if other["end"] == cds_end:
                continue  # Skip self
            # An adjacent CDS starts within our downstream window
            if raw_start <= other["start"] <= raw_end:
                min_downstream_start = min(min_downstream_start, other["start"])
            # An adjacent CDS completely spans our downstream window
            elif other["start"] <= raw_start and other["end"] >= raw_end:
                return None, None

        prom_start = raw_start
        prom_end = min_downstream_start

    # Check minimum length requirement
    if prom_end - prom_start < MIN_PROMOTER_LENGTH:
        return None, None

    return prom_start, prom_end


def extract_all_promoters(
    genbank_path: Path,
    output_fasta: Path,
    promoter_length: int = 300,
) -> pd.DataFrame:
    """
    Extract promoter regions for every annotated CDS in the genome.

    This is the full genome-wide promoter set used for FIMO scanning.
    Contains ~3,000–4,000 sequences after filtering overlapping regions.

    Args:
        genbank_path: Path to downloaded M. tb H37Rv GenBank file
        output_fasta: Path to write all promoter sequences (FASTA)
        promoter_length: Upstream region length in bp (default: 300)

    Returns:
        DataFrame with columns: locus_tag, gene, product, strand,
        prom_start, prom_end, prom_length, gc_content
    """
    logger.info(f"Parsing GenBank file: {genbank_path}")
    logger.info("This may take 10-30 seconds for the M. tb genome...")

    record = SeqIO.read(str(genbank_path), "genbank")
    genome_seq = str(record.seq).upper()
    genome_len = len(genome_seq)
    logger.info(f"Genome loaded: {genome_len:,} bp, {len(record.features)} features")

    cds_list = build_cds_interval_tree(record)
    logger.info(f"CDS features: {len(cds_list)}")

    seq_records = []
    metadata_rows = []
    seen_loci = set()
    skipped_overlap = 0
    skipped_duplicate = 0

    for cds in cds_list:
        locus_tag = cds["locus_tag"]

        # Skip duplicate locus tags (some genes have multiple CDS entries)
        if locus_tag in seen_loci:
            skipped_duplicate += 1
            continue
        seen_loci.add(locus_tag)

        prom_start, prom_end = find_trimmed_promoter_bounds(
            cds["start"],
            cds["end"],
            cds["strand"],
            cds_list,
            genome_len,
            promoter_length,
        )

        if prom_start is None:
            skipped_overlap += 1
            continue

        # Extract sequence and reverse-complement for minus strand genes
        prom_seq = genome_seq[prom_start:prom_end]
        if cds["strand"] == -1:
            from Bio.Seq import Seq as BioSeq
            prom_seq = str(BioSeq(prom_seq).reverse_complement())

        actual_length = len(prom_seq)
        gc = compute_gc_content(prom_seq)

        description = (
            f"gene={cds['gene'] or locus_tag} "
            f"strand={'+' if cds['strand'] == 1 else '-'} "
            f"prom_start={prom_start} prom_end={prom_end} "
            f"gc_content={gc:.3f} length={actual_length}"
        )
        seq_records.append(
            SeqRecord(Seq(prom_seq), id=locus_tag, description=description)
        )
        metadata_rows.append(
            {
                "locus_tag": locus_tag,
                "gene": cds["gene"],
                "product": cds["product"],
                "strand": "+" if cds["strand"] == 1 else "-",
                "prom_start": prom_start,
                "prom_end": prom_end,
                "prom_length": actual_length,
                "gc_content": gc,
            }
        )

    output_fasta.parent.mkdir(parents=True, exist_ok=True)
    SeqIO.write(seq_records, str(output_fasta), "fasta")

    logger.info(
        f"Extracted {len(seq_records)} promoters "
        f"(skipped: {skipped_overlap} overlapping, {skipped_duplicate} duplicate) "
        f"→ {output_fasta}"
    )

    return pd.DataFrame(metadata_rows)


def extract_divergent_igrs(
    genbank_path: Path,
    output_fasta: Path,
    min_igr_length: int = 400,
) -> pd.DataFrame:
    """
    Detect and extract full intergenic regions between divergent gene pairs.

    A divergent pair occurs when a minus-strand gene is immediately followed
    (in genomic coordinate order) by a plus-strand gene, with a large gap
    between them. Both genes share this intergenic space as their promoter
    region.

    This is biologically critical for Mce3R because the mce3R–yrbE3A
    intergenic region (Rv1963c[-] → Rv1964[+], 897 bp) is the primary
    Mce3R operator region. Standard 300 bp upstream extraction misses the
    central 297 bp of this 897 bp IGR.

    Args:
        genbank_path: Path to M. tb H37Rv GenBank file
        output_fasta: Path to write divergent IGR sequences (FASTA)
        min_igr_length: Only extract IGRs longer than this (default: 400 bp)

    Returns:
        DataFrame with columns: igr_id, gene1_tag, gene2_tag,
        igr_start, igr_end, igr_length, gc_content
    """
    logger.info(f"Detecting divergent gene pairs with IGR > {min_igr_length} bp...")
    record = SeqIO.read(str(genbank_path), "genbank")
    genome_seq = str(record.seq).upper()

    cds_list = build_cds_interval_tree(record)

    igr_records = []
    metadata_rows = []

    for i in range(len(cds_list) - 1):
        g1 = cds_list[i]
        g2 = cds_list[i + 1]

        # Divergent pair: minus-strand gene immediately before a plus-strand gene
        if g1["strand"] != -1 or g2["strand"] != 1:
            continue

        igr_start = g1["end"]    # end of minus-strand gene CDS
        igr_end = g2["start"]    # start of plus-strand gene CDS
        igr_len = igr_end - igr_start

        if igr_len < min_igr_length:
            continue

        # Skip if another CDS lies entirely within this IGR
        inner_cds = [c for c in cds_list
                     if c["start"] >= igr_start and c["end"] <= igr_end
                     and c["locus_tag"] not in (g1["locus_tag"], g2["locus_tag"])]
        if inner_cds:
            continue

        igr_seq = genome_seq[igr_start:igr_end]
        gc = compute_gc_content(igr_seq)
        igr_id = f"IGR_{g1['locus_tag']}_{g2['locus_tag']}"

        description = (
            f"divergent_igr gene1={g1['locus_tag']} gene2={g2['locus_tag']} "
            f"igr_start={igr_start} igr_end={igr_end} "
            f"gc_content={gc:.3f} length={igr_len}"
        )
        igr_records.append(SeqRecord(Seq(igr_seq), id=igr_id, description=description))
        metadata_rows.append({
            "igr_id": igr_id,
            "gene1_tag": g1["locus_tag"],
            "gene2_tag": g2["locus_tag"],
            "gene1_name": g1["gene"],
            "gene2_name": g2["gene"],
            "igr_start": igr_start,
            "igr_end": igr_end,
            "igr_length": igr_len,
            "gc_content": gc,
        })

    output_fasta.parent.mkdir(parents=True, exist_ok=True)
    SeqIO.write(igr_records, str(output_fasta), "fasta")
    logger.info(f"Extracted {len(igr_records)} divergent IGRs > {min_igr_length} bp → {output_fasta}")

    # Specifically call out the mce3R-yrbE3A IGR
    for row in metadata_rows:
        if row["gene1_tag"] == "Rv1963c" and row["gene2_tag"] == "Rv1964":
            logger.info(
                f"  *** mce3R-yrbE3A IGR: {row['igr_length']} bp "
                f"(genomic {row['igr_start']}-{row['igr_end']}) "
                f"GC={row['gc_content']:.3f} — primary Mce3R operator region ***"
            )

    return pd.DataFrame(metadata_rows)


def extract_candidate_promoters(
    all_promoters_fasta: Path,
    candidate_fasta: Path,
    candidate_genes: dict = None,
) -> pd.DataFrame:
    """
    Filter the full promoter set to Mce3R candidate genes only.

    The candidate promoter FASTA is used as MEME input to discover
    the Mce3R binding motif. A smaller, focused input set improves
    MEME's ability to find the true binding motif against background.

    Args:
        all_promoters_fasta: Path to full genome-wide promoters FASTA
        candidate_fasta: Output path for candidate-only FASTA
        candidate_genes: Dict mapping locus_tag → (gene_name, confidence, rationale)
                         Uses MCE3R_CANDIDATE_GENES if None.

    Returns:
        DataFrame of candidate promoters found, with candidate metadata.
    """
    if candidate_genes is None:
        candidate_genes = MCE3R_CANDIDATE_GENES

    logger.info(
        f"Filtering to {len(candidate_genes)} Mce3R candidate genes..."
    )

    # Index all promoters by locus tag
    all_records = {
        r.id: r
        for r in SeqIO.parse(str(all_promoters_fasta), "fasta")
    }

    found = []
    selected_records = []

    for locus_tag, (gene_name, confidence, rationale) in candidate_genes.items():
        record = all_records.get(locus_tag)
        if record is None:
            logger.debug(f"  Not found in promoters: {locus_tag} ({gene_name})")
            continue

        # Annotate header with candidate metadata
        new_desc = (
            f"{record.description} "
            f"candidate_gene={gene_name} "
            f"confidence={confidence}"
        )
        new_record = SeqRecord(record.seq, id=locus_tag, description=new_desc)
        selected_records.append(new_record)

        found.append(
            {
                "locus_tag": locus_tag,
                "gene": gene_name,
                "confidence": confidence,
                "rationale": rationale,
                "promoter_length": len(record.seq),
                "gc_content": compute_gc_content(str(record.seq)),
            }
        )

    candidate_fasta.parent.mkdir(parents=True, exist_ok=True)
    SeqIO.write(selected_records, str(candidate_fasta), "fasta")

    n_found = len(found)
    n_expected = len(candidate_genes)
    logger.info(
        f"Found {n_found}/{n_expected} candidate promoters → {candidate_fasta}"
    )

    if n_found < n_expected:
        missing = set(candidate_genes) - {r["locus_tag"] for r in found}
        logger.info(
            f"Note: {len(missing)} candidate loci were excluded: {', '.join(sorted(missing))}\n"
            "  These are likely interior operon genes co-transcribed from the first gene's\n"
            "  promoter (e.g., Rv1965-Rv1975 are transcribed from the Rv1964 promoter).\n"
            "  Only genes with ≥50 bp of non-overlapping intergenic upstream space are included."
        )

    df = pd.DataFrame(found)

    # Print candidate summary
    if not df.empty:
        primary = df[df["confidence"] == "primary"]
        extended = df[df["confidence"] == "extended"]
        logger.info(
            f"Candidate set: {len(primary)} primary targets, "
            f"{len(extended)} extended candidates"
        )

    return df


def main():
    root = get_project_root()
    parser = argparse.ArgumentParser(
        description=(
            "Extract M. tuberculosis promoter regions from GenBank annotation.\n"
            "Produces two FASTA files:\n"
            "  1. All promoters genome-wide (for FIMO scanning)\n"
            "  2. Mce3R candidate gene promoters (for MEME motif discovery)"
        )
    )
    parser.add_argument(
        "--genbank",
        type=Path,
        default=root / "data" / "raw" / "Mtb_H37Rv.gb",
        help="M. tb H37Rv GenBank file (from download_genome.py)",
    )
    parser.add_argument(
        "--promoter-length",
        type=int,
        default=300,
        help="Upstream region length in bp (default: 300)",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=root / "data" / "raw",
        help="Output directory for FASTA files",
    )
    parser.add_argument(
        "--metadata-dir",
        type=Path,
        default=root / "data" / "processed",
        help="Directory for metadata CSV files",
    )
    args = parser.parse_args()

    all_fasta = args.output_dir / "promoters.fasta"
    candidate_fasta = args.output_dir / "mce3r_candidate_promoters.fasta"
    igr_fasta = args.output_dir / "divergent_igrs.fasta"
    all_metadata = args.metadata_dir / "promoter_metadata.csv"
    candidate_metadata = args.metadata_dir / "candidate_promoter_metadata.csv"

    # Step 1: Extract all genome-wide promoters
    all_df = extract_all_promoters(
        genbank_path=args.genbank,
        output_fasta=all_fasta,
        promoter_length=args.promoter_length,
    )
    args.metadata_dir.mkdir(parents=True, exist_ok=True)
    all_df.to_csv(all_metadata, index=False)
    logger.info(f"Saved promoter metadata: {all_metadata}")

    # Step 1b: Extract full divergent IGRs (e.g., mce3R-yrbE3A 897 bp region)
    igr_df = extract_divergent_igrs(
        genbank_path=args.genbank,
        output_fasta=igr_fasta,
        min_igr_length=400,
    )
    # Append IGR sequences to the main promoters FASTA so FIMO scans them
    if not igr_df.empty:
        igr_records = list(SeqIO.parse(str(igr_fasta), "fasta"))
        with open(str(all_fasta), "a") as f_out:
            SeqIO.write(igr_records, f_out, "fasta")
        # Add IGR rows to metadata so they get tier labels
        igr_meta_rows = []
        for _, row in igr_df.iterrows():
            igr_meta_rows.append({
                "locus_tag": row["igr_id"],
                "gene": f"{row['gene1_tag']}/{row['gene2_tag']}",
                "product": "divergent intergenic region",
                "strand": "+",
                "prom_start": row["igr_start"],
                "prom_end": row["igr_end"],
                "prom_length": row["igr_length"],
                "gc_content": row["gc_content"],
            })
        igr_meta_df = pd.DataFrame(igr_meta_rows)
        all_df = pd.concat([all_df, igr_meta_df], ignore_index=True)
        all_df.to_csv(all_metadata, index=False)
        logger.info(f"Appended {len(igr_records)} divergent IGRs to {all_fasta}")

    # Step 2: Filter to Mce3R candidate genes
    candidate_df = extract_candidate_promoters(
        all_promoters_fasta=all_fasta,
        candidate_fasta=candidate_fasta,
    )
    candidate_df.to_csv(candidate_metadata, index=False)
    logger.info(f"Saved candidate metadata: {candidate_metadata}")

    # Summary
    print("\n" + "=" * 60)
    print("  PROMOTER EXTRACTION SUMMARY")
    print("=" * 60)
    print(f"  Genome            : M. tuberculosis H37Rv (NC_000962.3)")
    print(f"  Promoter length   : {args.promoter_length} bp upstream")
    print(f"  All promoters     : {len(all_df)} sequences → {all_fasta.name}")
    print(f"  Candidate promoters: {len(candidate_df)} sequences → {candidate_fasta.name}")
    print("=" * 60)

    if not candidate_df.empty:
        print("\nMce3R candidate genes:")
        for _, row in candidate_df.iterrows():
            print(
                f"  {row['locus_tag']:<12} {row['gene']:<12} "
                f"[{row['confidence']}]  {row['promoter_length']} bp"
            )

    print(f"\nOutput files:")
    print(f"  All promoters      : {all_fasta}")
    print(f"  Candidate promoters: {candidate_fasta}")
    print(f"  Promoter metadata  : {all_metadata}")
    print(f"  Candidate metadata : {candidate_metadata}")


if __name__ == "__main__":
    main()
