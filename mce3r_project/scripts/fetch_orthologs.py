"""
fetch_orthologs.py — Fetch public mycobacterial yrbE3A upstream (operator) regions.

The de novo discovery, the knowledge-based PWM, and the conservation gate all need more
than the 3 ortholog sequences shipped in the repo. This pulls additional orthologous
yrbE3A promoter/operator regions from NCBI (public data only), so MEME has enough real,
conserved sequence to discover the operator without it being injected.

Method (all public NCBI):
  1. esearch the Gene database for yrbE3A across Mycobacterium.
  2. esummary each gene -> genomic accession + coordinates + strand.
  3. efetch the upstream window (default 250 bp, strand-aware) = the divergent
     mce3R-yrbE3A intergenic region that carries the operator.
  4. De-duplicate to one sequence per species and merge with the existing orthologs.

Usage:
    python fetch_orthologs.py --email you@example.com --max-species 20 --upstream 250
"""

import argparse
import sys
import time
from pathlib import Path

from Bio import Entrez, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

sys.path.insert(0, str(Path(__file__).resolve().parent))
from utils import compute_gc_content, get_project_root, setup_logging

logger = setup_logging(__name__)

SEARCH_TERM = "(yrbE3A[gene] OR yrbE3A[All Fields]) AND Mycobacteriaceae[orgn]"
REQUEST_PAUSE_S = 0.4


def _gene_uids(term: str, retmax: int) -> list[str]:
    handle = Entrez.esearch(db="gene", term=term, retmax=retmax)
    rec = Entrez.read(handle)
    handle.close()
    return list(rec.get("IdList", []))


def _gene_genomic_info(uid: str) -> dict | None:
    """Return {organism, accession, start0, stop0, strand} for a gene UID, or None."""
    handle = Entrez.esummary(db="gene", id=uid)
    summ = Entrez.read(handle)
    handle.close()
    docs = summ.get("DocumentSummarySet", {}).get("DocumentSummary", [])
    if not docs:
        return None
    doc = docs[0]
    organism = doc.get("Organism", {}).get("ScientificName", "") or doc.get("Name", "")
    name = str(doc.get("Name", ""))
    description = str(doc.get("Description", ""))
    ginfo = doc.get("GenomicInfo", [])
    if not ginfo:
        return None
    g = ginfo[0]
    try:
        start0 = int(g["ChrStart"])
        stop0 = int(g["ChrStop"])
        acc = g["ChrAccVer"]
    except (KeyError, ValueError):
        return None
    strand = "+" if start0 <= stop0 else "-"
    return {
        "organism": organism,
        "name": name,
        "description": description,
        "accession": acc,
        "start0": start0,
        "stop0": stop0,
        "strand": strand,
    }


def _is_yrbe3a(info: dict) -> bool:
    """Keep only true yrbE3A orthologs (guard against the broad text search)."""
    name = info.get("name", "").lower()
    desc = info.get("description", "").lower()
    return name.startswith("yrbe3") or "yrbe3a" in name or "yrbe3a" in desc


def _fetch_upstream(info: dict, upstream: int) -> str | None:
    """efetch the strand-aware upstream window (gene 5' is ChrStart)."""
    acc = info["accession"]
    if info["strand"] == "+":
        seq_start = max(1, info["start0"] - upstream + 1)
        seq_stop = info["start0"]  # 1-based inclusive end == 0-based start
        strand_flag = 1
    else:
        seq_start = info["start0"] + 1  # 0-based start -> 1-based
        seq_stop = info["start0"] + upstream
        strand_flag = 2
    if seq_stop <= seq_start:
        return None
    try:
        handle = Entrez.efetch(
            db="nuccore",
            id=acc,
            rettype="fasta",
            retmode="text",
            seq_start=seq_start,
            seq_stop=seq_stop,
            strand=strand_flag,
        )
        rec = SeqIO.read(handle, "fasta")
        handle.close()
    except Exception as e:
        logger.warning(f"efetch failed for {acc}: {e}")
        return None
    seq = str(rec.seq).upper()
    return seq if set(seq) <= set("ACGTN") and len(seq) >= upstream // 2 else None


def fetch_orthologs(
    output_fasta: Path,
    email: str,
    max_species: int = 20,
    upstream: int = 250,
    api_key: str | None = None,
    merge_existing: Path | None = None,
) -> int:
    """
    Fetch yrbE3A upstream regions across Mycobacterium, one per species, and write FASTA.
    If merge_existing is given, its records are kept (de-duplicated by species).
    """
    Entrez.email = email
    if api_key:
        Entrez.api_key = api_key

    records = []
    seen_species = set()

    if merge_existing and Path(merge_existing).exists():
        for rec in SeqIO.parse(str(merge_existing), "fasta"):
            sp = rec.id.split("_")[0]
            seen_species.add(sp.lower())
            records.append(rec)
        logger.info(f"Kept {len(records)} existing ortholog records.")

    uids = _gene_uids(SEARCH_TERM, retmax=max_species * 6)
    logger.info(f"Gene UIDs returned: {len(uids)}")

    for uid in uids:
        if len(records) >= max_species:
            break
        time.sleep(REQUEST_PAUSE_S)
        info = _gene_genomic_info(uid)
        if not info or not info["organism"]:
            continue
        if not _is_yrbe3a(info):
            continue  # broad text search may return non-orthologs
        species_key = "_".join(info["organism"].split()[:2]).lower()
        if species_key in seen_species:
            continue
        time.sleep(REQUEST_PAUSE_S)
        seq = _fetch_upstream(info, upstream)
        if not seq:
            continue
        seen_species.add(species_key)
        sid = "".join(
            c if c.isalnum() else "_" for c in info["organism"].split()[:2][-1]
        )
        rec_id = f"{species_key}_yrbE3A_upstream"
        desc = (
            f"accession={info['accession']} strand={info['strand']} "
            f"organism={info['organism']} region={upstream}bp_upstream_yrbE3A "
            f"gc_content={compute_gc_content(seq):.3f} tier=Tier1_known_operator"
        )
        records.append(SeqRecord(Seq(seq), id=rec_id, description=desc))
        logger.info(f"  + {info['organism']} ({info['accession']})  {len(seq)} bp")

    out = Path(output_fasta)
    out.parent.mkdir(parents=True, exist_ok=True)
    SeqIO.write(records, str(out), "fasta")
    logger.info(f"Wrote {len(records)} ortholog upstream regions -> {out}")
    return len(records)


def main():
    root = get_project_root()
    parser = argparse.ArgumentParser(
        description="Fetch mycobacterial yrbE3A upstream orthologs from NCBI"
    )
    parser.add_argument(
        "--email",
        default="alwaniaayan6@gmail.com",
        help="Email for NCBI Entrez (required by policy)",
    )
    parser.add_argument("--api-key", default=None)
    parser.add_argument("--max-species", type=int, default=20)
    parser.add_argument("--upstream", type=int, default=250)
    parser.add_argument(
        "--output",
        type=Path,
        default=root / "data" / "raw" / "ortholog_promoters.fasta",
    )
    parser.add_argument(
        "--merge-existing",
        type=Path,
        default=root / "data" / "raw" / "ortholog_promoters.fasta",
    )
    args = parser.parse_args()

    n = fetch_orthologs(
        args.output,
        args.email,
        max_species=args.max_species,
        upstream=args.upstream,
        api_key=args.api_key,
        merge_existing=args.merge_existing,
    )
    print(f"\nTotal ortholog records: {n} -> {args.output}")


if __name__ == "__main__":
    main()
