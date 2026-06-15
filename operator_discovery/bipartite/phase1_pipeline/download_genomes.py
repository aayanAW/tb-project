"""
phase1_pipeline/download_genomes.py — Download reference genomes and annotations from NCBI.

Downloads:
  - M. tuberculosis H37Rv genome (NC_000962.3) as GenBank + FASTA
  - M. bovis genome (NC_002945.4) as FASTA
  - M. marinum genome (NC_010612.1) as FASTA
  - H37Rv GFF3 annotation
"""

import sys
import os
import time

sys.path.insert(0, '/Users/aayanalwani/tb project/mce3r_stochastic')
from config.parameters import PARAMS

from Bio import Entrez, SeqIO

Entrez.email = "mce3r_project@example.com"

PROJECT_ROOT = '/Users/aayanalwani/tb project/mce3r_stochastic'
GENOME_DIR = os.path.join(PROJECT_ROOT, 'data', 'genomes')


def ensure_dirs():
    os.makedirs(GENOME_DIR, exist_ok=True)


def download_genome_fasta(accession, outpath, retries=3):
    """Download a genome in FASTA format via Bio.Entrez."""
    if os.path.exists(outpath) and os.path.getsize(outpath) > 1000:
        print(f"  Already exists: {outpath}")
        return
    for attempt in range(retries):
        try:
            print(f"  Downloading {accession} (FASTA) attempt {attempt+1}...")
            handle = Entrez.efetch(db="nucleotide", id=accession,
                                   rettype="fasta", retmode="text")
            data = handle.read()
            handle.close()
            with open(outpath, 'w') as f:
                f.write(data)
            print(f"  Saved: {outpath} ({os.path.getsize(outpath)} bytes)")
            return
        except Exception as e:
            print(f"  Attempt {attempt+1} failed: {e}")
            if attempt < retries - 1:
                time.sleep(3)
    raise RuntimeError(f"Failed to download {accession} after {retries} attempts")


def download_genome_genbank(accession, outpath, retries=3):
    """Download a genome in GenBank format via Bio.Entrez."""
    if os.path.exists(outpath) and os.path.getsize(outpath) > 1000:
        print(f"  Already exists: {outpath}")
        return
    for attempt in range(retries):
        try:
            print(f"  Downloading {accession} (GenBank) attempt {attempt+1}...")
            handle = Entrez.efetch(db="nucleotide", id=accession,
                                   rettype="gbwithparts", retmode="text")
            data = handle.read()
            handle.close()
            with open(outpath, 'w') as f:
                f.write(data)
            print(f"  Saved: {outpath} ({os.path.getsize(outpath)} bytes)")
            return
        except Exception as e:
            print(f"  Attempt {attempt+1} failed: {e}")
            if attempt < retries - 1:
                time.sleep(3)
    raise RuntimeError(f"Failed to download {accession} GenBank after {retries} attempts")


def download_gff(accession, outpath, retries=3):
    """Download GFF3 annotation via Bio.Entrez."""
    if os.path.exists(outpath) and os.path.getsize(outpath) > 1000:
        print(f"  Already exists: {outpath}")
        return
    for attempt in range(retries):
        try:
            print(f"  Downloading {accession} GFF3 attempt {attempt+1}...")
            handle = Entrez.efetch(db="nucleotide", id=accession,
                                   rettype="gff3", retmode="text")
            data = handle.read()
            handle.close()
            with open(outpath, 'w') as f:
                f.write(data)
            print(f"  Saved: {outpath} ({os.path.getsize(outpath)} bytes)")
            return
        except Exception as e:
            print(f"  Attempt {attempt+1} failed: {e}")
            if attempt < retries - 1:
                time.sleep(3)
    raise RuntimeError(f"Failed to download {accession} GFF3 after {retries} attempts")


def run_download():
    """Main download routine. Returns dict of output file paths."""
    ensure_dirs()
    outputs = {}

    # H37Rv FASTA
    h37rv_fasta = os.path.join(GENOME_DIR, 'H37Rv.fasta')
    download_genome_fasta(PARAMS['H37Rv_accession'], h37rv_fasta)
    outputs['h37rv_fasta'] = h37rv_fasta

    # H37Rv GenBank
    h37rv_gb = os.path.join(GENOME_DIR, 'H37Rv.gb')
    download_genome_genbank(PARAMS['H37Rv_accession'], h37rv_gb)
    outputs['h37rv_gb'] = h37rv_gb

    # M. bovis FASTA
    bovis_fasta = os.path.join(GENOME_DIR, 'M_bovis.fasta')
    download_genome_fasta(PARAMS['M_bovis_accession'], bovis_fasta)
    outputs['bovis_fasta'] = bovis_fasta

    # M. marinum FASTA
    marinum_fasta = os.path.join(GENOME_DIR, 'M_marinum.fasta')
    download_genome_fasta(PARAMS['M_marinum_accession'], marinum_fasta)
    outputs['marinum_fasta'] = marinum_fasta

    # H37Rv GFF3
    h37rv_gff = os.path.join(GENOME_DIR, 'H37Rv.gff')
    download_gff(PARAMS['H37Rv_accession'], h37rv_gff)
    outputs['h37rv_gff'] = h37rv_gff

    # Brief delay for NCBI rate limit
    time.sleep(1)

    return outputs


if __name__ == "__main__":
    print("=" * 60)
    print("download_genomes.py — Self-test")
    print("=" * 60)

    n_sanity_pass = 0
    n_sanity_fail = 0
    n_science_expected = 0
    n_science_unexpected = 0

    try:
        outputs = run_download()
    except Exception as e:
        print(f"SANITY FAIL: Download raised exception: {e}")
        sys.exit(1)

    # --- SANITY CHECKS (hard-fail) ---

    # 1. All output files exist
    for key, path in outputs.items():
        if os.path.exists(path) and os.path.getsize(path) > 100:
            print(f"SANITY PASS: {key} exists ({os.path.getsize(path)} bytes)")
            n_sanity_pass += 1
        else:
            print(f"SANITY FAIL: {key} missing or empty at {path}")
            n_sanity_fail += 1

    # 2. H37Rv FASTA parseable with correct accession
    try:
        rec = SeqIO.read(outputs['h37rv_fasta'], 'fasta')
        assert len(rec.seq) > 0, "Empty sequence"
        print(f"SANITY PASS: H37Rv FASTA parseable, length={len(rec.seq)}")
        n_sanity_pass += 1
    except Exception as e:
        print(f"SANITY FAIL: H37Rv FASTA parse error: {e}")
        n_sanity_fail += 1

    # 3. GFF file contains feature lines
    try:
        gff_lines = 0
        with open(outputs['h37rv_gff']) as f:
            for line in f:
                if not line.startswith('#') and '\tCDS\t' in line:
                    gff_lines += 1
        assert gff_lines > 100, f"Only {gff_lines} CDS lines found"
        print(f"SANITY PASS: GFF contains {gff_lines} CDS lines")
        n_sanity_pass += 1
    except Exception as e:
        print(f"SANITY FAIL: GFF check error: {e}")
        n_sanity_fail += 1

    # --- SCIENTIFIC EXPECTATIONS (warn only) ---

    # 4. H37Rv genome size ~4,411,532 bp
    try:
        rec = SeqIO.read(outputs['h37rv_fasta'], 'fasta')
        genome_len = len(rec.seq)
        expected = PARAMS['H37Rv_genome_size']
        if abs(genome_len - expected) < 100:
            print(f"SCIENCE EXPECTED: H37Rv length={genome_len} (expected ~{expected})")
            n_science_expected += 1
        else:
            print(f"SCIENCE WARN: H37Rv length={genome_len}, expected ~{expected}")
            n_science_unexpected += 1
    except Exception as e:
        print(f"SCIENCE WARN: Could not check genome size: {e}")
        n_science_unexpected += 1

    # 5. H37Rv GC content ~65.6%
    try:
        seq_str = str(rec.seq).upper()
        gc = (seq_str.count('G') + seq_str.count('C')) / len(seq_str)
        if abs(gc - PARAMS['H37Rv_gc_content']) < 0.01:
            print(f"SCIENCE EXPECTED: H37Rv GC={gc:.4f} (expected ~{PARAMS['H37Rv_gc_content']})")
            n_science_expected += 1
        else:
            print(f"SCIENCE WARN: H37Rv GC={gc:.4f}, expected ~{PARAMS['H37Rv_gc_content']}")
            n_science_unexpected += 1
    except Exception as e:
        print(f"SCIENCE WARN: Could not check GC content: {e}")
        n_science_unexpected += 1

    # --- Summary ---
    print(f"\n{'=' * 50}")
    print(f"SANITY:  {n_sanity_pass} pass, {n_sanity_fail} fail")
    print(f"SCIENCE: {n_science_expected} expected, {n_science_unexpected} unexpected")

    if n_sanity_fail > 0:
        print("STOPPING: Sanity check failed.")
        sys.exit(1)
    else:
        print("All sanity checks passed.")
        sys.exit(0)
