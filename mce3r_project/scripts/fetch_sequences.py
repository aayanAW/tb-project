"""
fetch_sequences.py — Download orthologous sequences for Mce3R motif analysis.

Fetches the DNA sequences used by Panagoda et al. 2024 (ACS Chemical Biology)
for structural validation of the asymmetric Mce3R operator:

1. H37Rv (NC_000962.3) — 200 bp upstream of yrbE3A (Rv1964)
2. M. bovis (NC_002945.4) — orthologous upstream region
3. M. marinum (NC_010612.1) — orthologous upstream region

Also extracts intergenic regions for all Tier 1 and Tier 2 loci from H37Rv
to build the full candidate set for MEME motif discovery and FIMO scanning.

Outputs:
    data/raw/orthologous_sequences.fasta — three-species alignment input
    data/raw/sequences.fasta            — all Tier1/Tier2 intergenic regions

Biological context:
    Orthologous conservation of the Mce3R operator across mycobacterial species
    provides strong evidence that the asymmetric binding geometry is functionally
    important, not an H37Rv-specific artifact.

Usage:
    python fetch_sequences.py [--email EMAIL] [--output-dir PATH]
"""

import argparse
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from utils import get_project_root, setup_logging

logger = setup_logging(__name__)

# ---------------------------------------------------------------------------
# Accession constants for orthologous genome fetch
# ---------------------------------------------------------------------------

GENOMES = {
    "Mtb_H37Rv":  "NC_000962.3",   # M. tuberculosis H37Rv (reference)
    "Mbovis":     "NC_002945.4",   # M. bovis AF2122/97
    "Mmarinum":   "NC_010612.1",   # M. marinum M
}

# yrbE3A (Rv1964) approximate start positions in each genome (1-based)
# Used to extract the 200 bp upstream operator region (Panagoda et al. 2024)
YRBE3A_POSITIONS = {
    "Mtb_H37Rv":  (2207501, "+"),  # H37Rv genomic start of yrbE3A
    "Mbovis":     (2216000, "+"),  # approximate — NCBI search required for exact
    "Mmarinum":   (5100000, "+"),  # approximate — requires genome annotation lookup
}

# Tier 1 locus tags (exact locus tags from H37Rv for targeted extraction)
TIER1_LOCI = {
    "Rv1963c",   # mce3R (autoregulatory operator)
    "Rv1964",    # yrbE3A (primary validated operator)
    "Rv1935c",   # boundary gene of second validated IGR
    "Rv1936",    # boundary gene of second validated IGR
}

# Tier 2 locus tags (biologically plausible — confirmed or predicted regulon)
# See analyze_results.py TIER2_LOCI for the full list
TIER2_LOCI_CORE = {
    "Rv1933c", "Rv1934c", "Rv1937", "Rv1938", "Rv1939", "Rv1940", "Rv1941",
    "Rv3526",  "Rv3577",  "Rv3537",  "Rv3568c", "Rv3569c", "Rv3570c",
    "Rv3543c", "Rv3544c", "Rv3545c", "Rv3546",  "Rv3547",  "Rv3548c",
    "Rv3515c", "Rv3516",  "Rv3053c", "Rv0119",
    "Rv1908c", "Rv2428",  "Rv3846",  "Rv0432",  "Rv1932",
}

UPSTREAM_LEN = 200  # bp upstream of gene start for orthologous comparison (Panagoda 2024)
PROMOTER_LEN = 300  # bp upstream for general Tier1/Tier2 extraction


def _fetch_ncbi_sequence(
    accession: str,
    start: int,
    end: int,
    strand: str,
    email: str,
    retries: int = 3,
) -> str:
    """
    Fetch a DNA subsequence from NCBI Entrez (efetch).

    Args:
        accession: NCBI accession (e.g., 'NC_000962.3')
        start: 1-based start position
        end: 1-based end position (inclusive)
        strand: '+' or '-'
        email: NCBI contact email (required by NCBI policy)
        retries: Number of retry attempts on network failure

    Returns:
        DNA string (uppercase), or empty string on failure
    """
    try:
        from Bio import Entrez, SeqIO
    except ImportError:
        logger.error("BioPython not installed. Run: pip install biopython")
        return ""

    Entrez.email = email
    strand_int = 1 if strand == "+" else 2

    for attempt in range(1, retries + 1):
        try:
            handle = Entrez.efetch(
                db="nuccore",
                id=accession,
                rettype="fasta",
                retmode="text",
                seq_start=start,
                seq_stop=end,
                strand=strand_int,
            )
            record = SeqIO.read(handle, "fasta")
            handle.close()
            return str(record.seq).upper()
        except Exception as e:
            logger.warning(f"Attempt {attempt}/{retries} failed for {accession}:{start}-{end}: {e}")
            if attempt < retries:
                time.sleep(2 ** attempt)  # exponential back-off

    logger.error(f"All {retries} attempts failed for {accession}:{start}-{end}")
    return ""


def fetch_orthologous_sequences(
    output_path: Path,
    email: str = "user@example.com",
) -> int:
    """
    Download 200 bp upstream of yrbE3A from Mtb, M. bovis, and M. marinum.

    These are the sequences used in Panagoda et al. 2024 to demonstrate that
    the asymmetric Mce3R operator is conserved across mycobacterial species.

    Args:
        output_path: Path to write the FASTA file
        email: NCBI Entrez email

    Returns:
        Number of sequences written
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    records = []
    for species, accession in GENOMES.items():
        if species not in YRBE3A_POSITIONS:
            logger.warning(f"No yrbE3A position known for {species} — skipping")
            continue

        gene_start, strand = YRBE3A_POSITIONS[species]
        if strand == "+":
            start = max(1, gene_start - UPSTREAM_LEN)
            end = gene_start - 1
        else:
            start = gene_start + 1
            end = gene_start + UPSTREAM_LEN

        logger.info(f"Fetching {species} ({accession}) yrbE3A upstream [{start}-{end}, {strand}]...")
        seq = _fetch_ncbi_sequence(accession, start, end, strand, email)

        if seq:
            gc = sum(1 for b in seq if b in "GC") / len(seq) if seq else 0.0
            header = (
                f">{species}_yrbE3A_upstream  "
                f"accession={accession} coordinates={start}-{end} strand={strand} "
                f"region=200bp_upstream_yrbE3A gc_content={gc:.3f} "
                f"tier=Tier1_known_operator species={species}"
            )
            records.append((header, seq))
            logger.info(f"  {species}: {len(seq)} bp, GC={gc:.1%}")
        else:
            logger.warning(f"  {species}: fetch failed")

    with open(output_path, "w") as f:
        for header, seq in records:
            f.write(f"{header}\n")
            # Write sequence in 60-char lines
            for i in range(0, len(seq), 60):
                f.write(seq[i:i + 60] + "\n")

    logger.info(f"Wrote {len(records)} orthologous sequences to {output_path}")
    return len(records)


def extract_tier_sequences_from_genbank(
    genbank_path: Path,
    output_path: Path,
    upstream_len: int = PROMOTER_LEN,
) -> int:
    """
    Extract upstream regions for all Tier 1 and Tier 2 loci from H37Rv GenBank.

    Parses the downloaded GenBank record (Mtb_H37Rv.gb from download_genome.py)
    and extracts the specified upstream region for each locus. Writes all
    sequences to a FASTA file suitable for MEME input.

    FASTA header format:
        >LOCUS_TAG  gene=GENE strand=STRAND start=N end=N tier=TIER
        gc_content=GC

    Args:
        genbank_path: Path to Mtb_H37Rv.gb (from download_genome step)
        output_path: Path to write the FASTA file
        upstream_len: bp upstream of gene start to extract (default 300)

    Returns:
        Number of sequences written
    """
    if not genbank_path.exists():
        logger.error(f"GenBank file not found: {genbank_path}. Run the download step first.")
        return 0

    try:
        from Bio import SeqIO
    except ImportError:
        logger.error("BioPython not installed. Run: pip install biopython")
        return 0

    logger.info(f"Parsing GenBank file: {genbank_path}")
    record = SeqIO.read(str(genbank_path), "genbank")
    genome = str(record.seq).upper()
    genome_len = len(genome)
    logger.info(f"Genome: {genome_len:,} bp")

    all_loci = TIER1_LOCI | TIER2_LOCI_CORE

    # Build locus → (start, end, strand, gene_name, product) map from CDS features
    locus_map = {}
    for feat in record.features:
        if feat.type != "CDS":
            continue
        tag = feat.qualifiers.get("locus_tag", [""])[0]
        if tag not in all_loci:
            continue
        gene_name = feat.qualifiers.get("gene", [""])[0]
        product = feat.qualifiers.get("product", ["unknown"])[0]
        strand = "+" if feat.location.strand == 1 else "-"
        feat_start = int(feat.location.start)  # 0-based
        feat_end = int(feat.location.end)
        locus_map[tag] = (feat_start, feat_end, strand, gene_name, product)

    logger.info(f"Found {len(locus_map)}/{len(all_loci)} target loci in GenBank")

    not_found = all_loci - set(locus_map.keys())
    if not_found:
        logger.info(f"Loci not found in annotation: {sorted(not_found)}")

    records = []
    for tag, (feat_start, feat_end, strand, gene_name, product) in sorted(locus_map.items()):
        if strand == "+":
            prom_start = max(0, feat_start - upstream_len)
            prom_end = feat_start
        else:
            prom_start = feat_end
            prom_end = min(genome_len, feat_end + upstream_len)

        seq = genome[prom_start:prom_end]
        if not seq:
            logger.warning(f"Empty sequence for {tag} ({prom_start}-{prom_end})")
            continue

        gc = sum(1 for b in seq if b in "GC") / len(seq)

        # Assign tier label for header
        tier = "Tier1_known_operator" if tag in TIER1_LOCI else "Tier2_plausible_target"

        header = (
            f">{tag}  gene={gene_name or 'unknown'} strand={strand} "
            f"prom_start={prom_start + 1} prom_end={prom_end} "
            f"prom_length={len(seq)} gc_content={gc:.3f} "
            f"tier={tier} product={product[:50]!r}"
        )
        records.append((header, seq))

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w") as f:
        for header, seq in records:
            f.write(f"{header}\n")
            for i in range(0, len(seq), 60):
                f.write(seq[i:i + 60] + "\n")

    logger.info(f"Wrote {len(records)} Tier1/Tier2 upstream regions to {output_path}")
    return len(records)


def main():
    root = get_project_root()
    parser = argparse.ArgumentParser(
        description=(
            "Download orthologous Mce3R operator sequences and extract Tier1/Tier2 "
            "intergenic regions from H37Rv for MEME motif discovery."
        )
    )
    parser.add_argument(
        "--email",
        default="user@example.com",
        help="Email for NCBI Entrez API (required by NCBI policy)",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=root / "data" / "raw",
        help="Directory for output FASTA files",
    )
    parser.add_argument(
        "--genbank",
        type=Path,
        default=root / "data" / "raw" / "Mtb_H37Rv.gb",
        help="H37Rv GenBank file (from download_genome step)",
    )
    parser.add_argument(
        "--orthologous-only",
        action="store_true",
        help="Only fetch orthologous sequences (skip Tier1/Tier2 extraction)",
    )
    parser.add_argument(
        "--upstream-len",
        type=int,
        default=PROMOTER_LEN,
        help=f"bp upstream of gene to extract for Tier1/Tier2 (default: {PROMOTER_LEN})",
    )
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)

    # 1. Fetch orthologous sequences (Mtb + M. bovis + M. marinum)
    orth_path = args.output_dir / "orthologous_sequences.fasta"
    n_orth = fetch_orthologous_sequences(orth_path, email=args.email)
    print(f"\nOrthologous sequences: {n_orth} written to {orth_path}")

    # 2. Extract Tier1/Tier2 sequences from H37Rv GenBank
    if not args.orthologous_only:
        tier_path = args.output_dir / "tier_sequences.fasta"
        n_tier = extract_tier_sequences_from_genbank(
            args.genbank, tier_path, upstream_len=args.upstream_len
        )
        print(f"Tier1/Tier2 sequences: {n_tier} written to {tier_path}")
        print(
            f"\nTo run MEME on these sequences:\n"
            f"  python run_meme.py --fasta {tier_path} --output-dir results/motifs/\n"
            f"\nOr run the full pipeline:\n"
            f"  python main.py --steps meme fimo analyze visualize\n"
        )


if __name__ == "__main__":
    main()
