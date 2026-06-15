"""
download_genome.py — Download M. tuberculosis H37Rv genome from NCBI.

Downloads the complete genome (GenBank format with annotation + FASTA) for
Mycobacterium tuberculosis H37Rv (NCBI accession NC_000962.3) using the
BioPython Entrez API.

The GenBank file contains:
- Complete nucleotide sequence (4,411,532 bp)
- Annotations for ~4,000 coding sequences (CDS)
- tRNA, rRNA, and regulatory feature annotations
- Gene names, locus tags (Rv numbers), and protein products

These are required by extract_promoters.py to identify upstream regions
for every annotated gene in the M. tb genome.

NCBI Entrez usage policy:
- An email address is required to use Entrez
- Rate limit: 3 requests/second (10/second with API key)
- Large files (>10 MB) should not be retrieved repeatedly

Usage:
    python download_genome.py [--email EMAIL] [--accession ACC] [--output-dir DIR]
"""

import argparse
import sys
import time
from pathlib import Path

from Bio import Entrez, SeqIO

sys.path.insert(0, str(Path(__file__).resolve().parent))
from utils import get_project_root, setup_logging

logger = setup_logging(__name__)

# Default NCBI accession for M. tuberculosis H37Rv complete genome
DEFAULT_ACCESSION = "NC_000962.3"
# Genome size in bp for progress reporting
MTB_GENOME_SIZE_BP = 4_411_532


def download_genbank(
    accession: str,
    output_path: Path,
    email: str,
    api_key: str = None,
    max_retries: int = 3,
) -> Path:
    """
    Download a GenBank record from NCBI Entrez in GenBank text format.

    The GenBank format includes both sequence and annotation, allowing
    extract_promoters.py to work from a single downloaded file.

    Args:
        accession: NCBI accession number (e.g., "NC_000962.3")
        output_path: File path to write the GenBank record
        email: Email for NCBI Entrez (required by NCBI policy)
        api_key: Optional NCBI API key (raises rate limit to 10 req/s)
        max_retries: Number of retry attempts on network failure

    Returns:
        Path to the written GenBank file
    """
    Entrez.email = email
    if api_key:
        Entrez.api_key = api_key

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    logger.info(f"Downloading GenBank record: {accession}")
    logger.info(
        f"  Genome: Mycobacterium tuberculosis H37Rv (~{MTB_GENOME_SIZE_BP:,} bp)"
    )
    logger.info(f"  Output: {output_path}")
    logger.info("  This may take 30-120 seconds depending on network speed...")

    # NC_000962.3 is a RefSeq CON (contig) record that references AL123456.3.
    # 'gbwithparts' expands the CON record, embeds all features + full sequence.
    # Fall back to plain 'gb' if gbwithparts is unavailable.
    for rettype in ("gbwithparts", "gb"):
        for attempt in range(1, max_retries + 1):
            try:
                handle = Entrez.efetch(
                    db="nucleotide",
                    id=accession,
                    rettype=rettype,
                    retmode="text",
                )
                genbank_text = handle.read()
                handle.close()

                # Basic sanity check: a valid GenBank with sequence must have ORIGIN
                if "ORIGIN" not in genbank_text:
                    logger.warning(
                        f"rettype={rettype} did not return sequence data "
                        f"(no ORIGIN section). Trying next format..."
                    )
                    break  # Try next rettype

                with open(output_path, "w") as f:
                    f.write(genbank_text)

                file_size_mb = output_path.stat().st_size / 1_048_576
                logger.info(
                    f"Downloaded GenBank file ({rettype}): "
                    f"{file_size_mb:.1f} MB → {output_path}"
                )
                return output_path

            except Exception as e:
                logger.warning(f"Attempt {attempt}/{max_retries} failed ({rettype}): {e}")
                if attempt < max_retries:
                    wait = 5 * attempt
                    logger.info(f"Retrying in {wait}s...")
                    time.sleep(wait)
                else:
                    logger.warning(f"rettype={rettype} exhausted retries. Trying next format...")
                    break

    raise RuntimeError(
        f"Failed to download {accession} with sequence content after {max_retries} attempts.\n"
        "Check your internet connection and NCBI Entrez access."
    )


def extract_fasta_from_genbank(genbank_path: Path, fasta_path: Path) -> Path:
    """
    Extract the nucleotide sequence from a GenBank file and save as FASTA.

    This creates a standalone FASTA file for use with other tools (e.g., BLAST,
    BEDTools) that do not read GenBank format.

    Args:
        genbank_path: Path to the downloaded GenBank file
        fasta_path: Output FASTA path

    Returns:
        Path to the written FASTA file
    """
    fasta_path = Path(fasta_path)
    fasta_path.parent.mkdir(parents=True, exist_ok=True)

    logger.info(f"Extracting FASTA sequence from {genbank_path}...")
    record = SeqIO.read(str(genbank_path), "genbank")

    with open(fasta_path, "w") as f:
        # Write FASTA with informative header
        f.write(
            f">{record.id} {record.description}\n"
        )
        seq_str = str(record.seq)
        for i in range(0, len(seq_str), 60):
            f.write(seq_str[i:i + 60] + "\n")

    logger.info(
        f"Saved genome FASTA: {len(record.seq):,} bp → {fasta_path}"
    )
    return fasta_path


def validate_genbank_file(genbank_path: Path) -> bool:
    """
    Validate that a downloaded GenBank file is complete and parseable.

    Checks:
    1. File exists and is non-empty
    2. BioPython can parse it as GenBank
    3. It contains CDS feature annotations (not just sequence)
    4. Sequence length is consistent with M. tb H37Rv

    Args:
        genbank_path: Path to GenBank file to validate

    Returns:
        True if valid, False otherwise
    """
    genbank_path = Path(genbank_path)

    if not genbank_path.exists() or genbank_path.stat().st_size < 1_000_000:
        logger.warning(f"GenBank file missing or too small: {genbank_path}")
        return False

    try:
        record = SeqIO.read(str(genbank_path), "genbank")
    except Exception as e:
        logger.warning(f"Could not parse GenBank file: {e}")
        return False

    seq_len = len(record.seq)
    n_features = len(record.features)
    n_cds = sum(1 for f in record.features if f.type == "CDS")

    logger.info(
        f"Validation: {seq_len:,} bp, {n_features} features, {n_cds} CDS"
    )

    # M. tb H37Rv should have ~4.4 Mb sequence and >3,000 CDS
    if seq_len < 4_000_000:
        logger.warning(f"Sequence length {seq_len:,} bp seems too short for M. tb H37Rv")
        return False

    if n_cds < 3000:
        logger.warning(f"Only {n_cds} CDS features — expected >3,000 for M. tb H37Rv")
        return False

    logger.info("GenBank file validation passed.")
    return True


def download_genome(
    accession: str = DEFAULT_ACCESSION,
    output_dir: Path = None,
    email: str = "user@example.com",
    api_key: str = None,
    skip_existing: bool = True,
) -> dict:
    """
    Main download function: fetches GenBank + extracts FASTA.

    Args:
        accession: NCBI accession number
        output_dir: Directory to save genome files (default: data/raw/)
        email: Email address for NCBI Entrez
        api_key: Optional NCBI API key
        skip_existing: If True, skip download if valid file already exists

    Returns:
        Dictionary with 'genbank_path' and 'fasta_path'
    """
    if output_dir is None:
        output_dir = get_project_root() / "data" / "raw"
    output_dir = Path(output_dir)

    genbank_path = output_dir / f"Mtb_H37Rv.gb"
    fasta_path = output_dir / f"Mtb_H37Rv.fasta"

    # Check if valid file already exists
    if skip_existing and genbank_path.exists():
        logger.info(f"GenBank file already exists: {genbank_path}")
        if validate_genbank_file(genbank_path):
            logger.info("Existing file is valid — skipping download.")
        else:
            logger.warning("Existing file is invalid — re-downloading.")
            skip_existing = False

    if not skip_existing or not genbank_path.exists():
        download_genbank(accession, genbank_path, email=email, api_key=api_key)

    # Always re-extract FASTA (fast operation)
    if not fasta_path.exists():
        extract_fasta_from_genbank(genbank_path, fasta_path)

    return {"genbank_path": genbank_path, "fasta_path": fasta_path}


def main():
    parser = argparse.ArgumentParser(
        description="Download M. tuberculosis H37Rv genome from NCBI"
    )
    parser.add_argument(
        "--accession",
        default=DEFAULT_ACCESSION,
        help=f"NCBI accession number (default: {DEFAULT_ACCESSION})",
    )
    parser.add_argument(
        "--email",
        default="user@example.com",
        help="Email address for NCBI Entrez API (required by NCBI policy)",
    )
    parser.add_argument(
        "--api-key",
        default=None,
        help="NCBI API key (optional, increases rate limit to 10 req/s)",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=get_project_root() / "data" / "raw",
        help="Output directory for genome files",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Re-download even if files already exist",
    )
    args = parser.parse_args()

    result = download_genome(
        accession=args.accession,
        output_dir=args.output_dir,
        email=args.email,
        api_key=args.api_key,
        skip_existing=not args.force,
    )

    print(f"\nGenome files:")
    print(f"  GenBank : {result['genbank_path']}")
    print(f"  FASTA   : {result['fasta_path']}")


if __name__ == "__main__":
    main()
