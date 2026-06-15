"""
phase1_pipeline/extract_upstream.py — Extract upstream regions from H37Rv genome.

- Parses GFF3, converts 1-based inclusive coords to 0-based Python indices immediately.
- Extracts 200 bp upstream of every CDS.
- Handles circular genome wrapping at BOTH origin AND terminus.
- Extracts yrbE3A ortholog upstream regions from 3 species for MEME input.
- Writes known_operator.fasta with the experimentally determined operator sequence.
"""

import sys
import os
import re

sys.path.insert(0, '/Users/aayanalwani/tb project/mce3r_stochastic')
from config.parameters import PARAMS

from Bio import SeqIO
from Bio.Seq import Seq

PROJECT_ROOT = '/Users/aayanalwani/tb project/mce3r_stochastic'
GENOME_DIR = os.path.join(PROJECT_ROOT, 'data', 'genomes')
SEQ_DIR = os.path.join(PROJECT_ROOT, 'data', 'sequences')


def ensure_dirs():
    os.makedirs(SEQ_DIR, exist_ok=True)


def load_genome(fasta_path):
    """Load a genome FASTA and return the SeqRecord."""
    rec = SeqIO.read(fasta_path, 'fasta')
    return rec


def extract_circular(seq_str, start, length, genome_len):
    """
    Extract a substring from a circular genome.
    start: 0-based start position
    length: number of bases to extract
    genome_len: total genome length
    Handles wrapping at both origin and terminus.
    """
    start = start % genome_len
    end = start + length
    if end <= genome_len:
        return seq_str[start:end]
    else:
        # Wraps around origin
        return seq_str[start:] + seq_str[:end - genome_len]


def reverse_complement(seq_str):
    """Return reverse complement of a DNA string."""
    return str(Seq(seq_str).reverse_complement())


def parse_gff_cds(gff_path):
    """
    Parse GFF3 file for CDS features.
    Converts 1-based inclusive GFF coords to 0-based Python indices immediately.
    Returns list of dicts with keys: seqid, start0, end0, strand, gene_id, locus_tag
    """
    features = []
    with open(gff_path) as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue
            if parts[2] != 'CDS':
                continue
            seqid = parts[0]
            # GFF is 1-based inclusive: convert to 0-based
            start_1based = int(parts[3])
            end_1based = int(parts[4])
            start0 = start_1based - 1  # 0-based start
            end0 = end_1based          # 0-based exclusive end (== 1-based end)
            strand = parts[6]

            # Parse attributes
            attrs = parts[8]
            gene_id = ''
            locus_tag = ''
            # Try to get gene name
            m = re.search(r'gene=([^;]+)', attrs)
            if m:
                gene_id = m.group(1)
            m = re.search(r'locus_tag=([^;]+)', attrs)
            if m:
                locus_tag = m.group(1)
            # Also try old_locus_tag
            old_locus = ''
            m = re.search(r'old_locus_tag=([^;]+)', attrs)
            if m:
                old_locus = m.group(1)

            features.append({
                'seqid': seqid,
                'start0': start0,
                'end0': end0,
                'strand': strand,
                'gene_id': gene_id,
                'locus_tag': locus_tag,
                'old_locus_tag': old_locus,
            })
    return features


def deduplicate_cds(features):
    """
    Deduplicate CDS features that share the same start position and strand.
    Keep the first occurrence (which typically has the most annotation).
    """
    seen = set()
    deduped = []
    for f in features:
        key = (f['start0'], f['end0'], f['strand'])
        if key not in seen:
            seen.add(key)
            deduped.append(f)
    return deduped


def extract_upstream_regions(features, genome_seq, upstream_len, genome_len):
    """
    Extract upstream_len bp upstream of each CDS.
    For + strand: upstream is BEFORE the start.
    For - strand: upstream is AFTER the end (reverse complement).
    Returns list of (header, sequence) tuples.
    """
    results = []
    for feat in features:
        if feat['strand'] == '+':
            # Upstream region is before start0
            up_start = feat['start0'] - upstream_len
            if up_start < 0:
                # Wraps around origin
                up_start = up_start % genome_len
            seq = extract_circular(genome_seq, up_start, upstream_len, genome_len)
        else:
            # Minus strand: upstream is after end0, on reverse strand
            up_start = feat['end0']
            seq = extract_circular(genome_seq, up_start, upstream_len, genome_len)
            seq = reverse_complement(seq)

        label = feat['locus_tag'] or feat['gene_id'] or f"CDS_{feat['start0']}"
        header = f"{label}_upstream_{upstream_len}bp"
        results.append((header, seq))
    return results


def write_fasta(sequences, outpath):
    """Write list of (header, seq) tuples to FASTA file."""
    with open(outpath, 'w') as f:
        for header, seq in sequences:
            f.write(f">{header}\n{seq}\n")


def find_ortholog_upstream(genome_seq, genome_len, genbank_path, target_locus_tags,
                           upstream_len):
    """
    Find an ortholog by locus_tag in a GenBank file and extract upstream region.
    target_locus_tags: list of locus_tag strings to try matching.
    Returns (header, seq) or None.
    """
    rec = SeqIO.read(genbank_path, 'genbank')
    seq_str = str(rec.seq).upper()
    glen = len(seq_str)

    for feature in rec.features:
        if feature.type != 'CDS':
            continue
        lt = feature.qualifiers.get('locus_tag', [''])[0]
        old_lt = feature.qualifiers.get('old_locus_tag', [''])[0]
        gene = feature.qualifiers.get('gene', [''])[0]

        matched = False
        for target in target_locus_tags:
            # Use substring matching to handle prefix/suffix variations
            # e.g. 'Mb1997' matches 'BQ2027_MB1997C'
            target_lower = target.lower()
            if (target_lower in lt.lower() or target_lower in old_lt.lower()
                    or target_lower in gene.lower()
                    or target in (lt, old_lt, gene)):
                matched = True
                break
        if not matched:
            continue

        # Found it
        strand = feature.location.strand  # 1 or -1
        if strand == 1:
            start0 = int(feature.location.start)
            up_start = (start0 - upstream_len) % glen
            seq = extract_circular(seq_str, up_start, upstream_len, glen)
        else:
            end0 = int(feature.location.end)
            seq = extract_circular(seq_str, end0, upstream_len, glen)
            seq = reverse_complement(seq)

        return (f"{lt}_upstream_{upstream_len}bp", seq)

    return None


def run_extract_upstream():
    """
    Main extraction routine.
    Returns dict with output file paths and counts.
    """
    ensure_dirs()
    upstream_len = PARAMS['upstream_length']

    # Load H37Rv genome
    h37rv_fasta = os.path.join(GENOME_DIR, 'H37Rv.fasta')
    h37rv_rec = load_genome(h37rv_fasta)
    genome_seq = str(h37rv_rec.seq).upper()
    genome_len = len(genome_seq)
    print(f"  H37Rv genome loaded: {genome_len} bp")

    # Parse GFF
    gff_path = os.path.join(GENOME_DIR, 'H37Rv.gff')
    features = parse_gff_cds(gff_path)
    print(f"  GFF parsed: {len(features)} CDS features")

    # Deduplicate
    features = deduplicate_cds(features)
    print(f"  After dedup: {len(features)} unique CDS features")

    # Extract upstream regions
    upstream_seqs = extract_upstream_regions(features, genome_seq, upstream_len, genome_len)
    print(f"  Extracted {len(upstream_seqs)} upstream regions")

    # Write all upstream sequences
    all_upstream_path = os.path.join(SEQ_DIR, 'all_upstream_200bp.fasta')
    write_fasta(upstream_seqs, all_upstream_path)
    print(f"  Written: {all_upstream_path}")

    # --- yrbE3A ortholog upstream for MEME input ---
    meme_seqs = []

    # H37Rv yrbE3A (Rv1964) upstream - find it
    for header, seq in upstream_seqs:
        if PARAMS['yrbE3A_gene'] in header or 'Rv1964' in header or 'yrbE3A' in header.lower():
            meme_seqs.append((f"H37Rv_{header}", seq))
            break
    else:
        # Try harder: search by locus tag pattern
        for feat in features:
            tags = [feat['locus_tag'], feat['old_locus_tag'], feat['gene_id']]
            if any('1964' in t for t in tags if t):
                if feat['strand'] == '+':
                    up_start = (feat['start0'] - upstream_len) % genome_len
                    seq = extract_circular(genome_seq, up_start, upstream_len, genome_len)
                else:
                    seq = extract_circular(genome_seq, feat['end0'], upstream_len, genome_len)
                    seq = reverse_complement(seq)
                meme_seqs.append((f"H37Rv_Rv1964_upstream_{upstream_len}bp", seq))
                break

    # M. bovis ortholog
    bovis_gb = os.path.join(GENOME_DIR, 'M_bovis.fasta')
    # We need GenBank for annotation - but we only downloaded FASTA.
    # Use H37Rv GenBank to find bovis ortholog position by homology
    # Alternative: download bovis GenBank
    bovis_genbank = os.path.join(GENOME_DIR, 'M_bovis.gb')
    if not os.path.exists(bovis_genbank):
        # Download GenBank for M. bovis
        from Bio import Entrez
        import time
        Entrez.email = "mce3r_project@example.com"
        try:
            print("  Downloading M. bovis GenBank for ortholog lookup...")
            handle = Entrez.efetch(db="nucleotide", id=PARAMS['M_bovis_accession'],
                                   rettype="gbwithparts", retmode="text")
            data = handle.read()
            handle.close()
            with open(bovis_genbank, 'w') as f:
                f.write(data)
            time.sleep(1)
        except Exception as e:
            print(f"  WARNING: Could not download M. bovis GenBank: {e}")

    if os.path.exists(bovis_genbank):
        bovis_rec = SeqIO.read(bovis_genbank, 'genbank')
        bovis_seq = str(bovis_rec.seq).upper()
        result = find_ortholog_upstream(
            bovis_seq, len(bovis_seq), bovis_genbank,
            [PARAMS['yrbE3A_bovis_locus'], 'Mb1997', 'yrbE3A'],
            upstream_len
        )
        if result:
            meme_seqs.append((f"M_bovis_{result[0]}", result[1]))
            print(f"  Found M. bovis yrbE3A ortholog upstream")
        else:
            print(f"  WARNING: M. bovis yrbE3A ortholog not found")

    # M. marinum ortholog
    marinum_genbank = os.path.join(GENOME_DIR, 'M_marinum.gb')
    if not os.path.exists(marinum_genbank):
        from Bio import Entrez
        import time
        Entrez.email = "mce3r_project@example.com"
        try:
            print("  Downloading M. marinum GenBank for ortholog lookup...")
            handle = Entrez.efetch(db="nucleotide", id=PARAMS['M_marinum_accession'],
                                   rettype="gbwithparts", retmode="text")
            data = handle.read()
            handle.close()
            with open(marinum_genbank, 'w') as f:
                f.write(data)
            time.sleep(1)
        except Exception as e:
            print(f"  WARNING: Could not download M. marinum GenBank: {e}")

    if os.path.exists(marinum_genbank):
        marinum_rec = SeqIO.read(marinum_genbank, 'genbank')
        marinum_seq = str(marinum_rec.seq).upper()
        result = find_ortholog_upstream(
            marinum_seq, len(marinum_seq), marinum_genbank,
            [PARAMS['yrbE3A_marinum_locus'], 'MMAR_2522', 'yrbE3A'],
            upstream_len
        )
        if result:
            meme_seqs.append((f"M_marinum_{result[0]}", result[1]))
            print(f"  Found M. marinum yrbE3A ortholog upstream")
        else:
            print(f"  WARNING: M. marinum yrbE3A ortholog not found")

    # Write MEME input (ortholog upstream regions)
    meme_input_path = os.path.join(SEQ_DIR, 'meme_input_orthologs.fasta')
    if meme_seqs:
        write_fasta(meme_seqs, meme_input_path)
        print(f"  MEME input written: {meme_input_path} ({len(meme_seqs)} sequences)")
    else:
        print("  WARNING: No MEME input sequences found")

    # Write known operator sequence
    operator_path = os.path.join(SEQ_DIR, 'known_operator.fasta')
    write_fasta([('Mce3R_operator_Panagoda2024', PARAMS['operator_sequence'])],
                operator_path)
    print(f"  Known operator written: {operator_path}")

    return {
        'all_upstream_fasta': all_upstream_path,
        'meme_input_fasta': meme_input_path,
        'known_operator_fasta': operator_path,
        'n_upstream_seqs': len(upstream_seqs),
        'n_meme_seqs': len(meme_seqs),
    }


if __name__ == "__main__":
    print("=" * 60)
    print("extract_upstream.py — Self-test")
    print("=" * 60)

    n_sanity_pass = 0
    n_sanity_fail = 0
    n_science_expected = 0
    n_science_unexpected = 0

    try:
        result = run_extract_upstream()
    except Exception as e:
        import traceback
        traceback.print_exc()
        print(f"SANITY FAIL: Exception during extraction: {e}")
        sys.exit(1)

    # --- SANITY CHECKS (hard-fail) ---

    # 1. Output files exist
    for key in ['all_upstream_fasta', 'known_operator_fasta']:
        path = result[key]
        if os.path.exists(path) and os.path.getsize(path) > 100:
            print(f"SANITY PASS: {key} exists ({os.path.getsize(path)} bytes)")
            n_sanity_pass += 1
        else:
            print(f"SANITY FAIL: {key} missing or empty")
            n_sanity_fail += 1

    # 2. Number of upstream sequences in range 3500-4500
    n_seqs = result['n_upstream_seqs']
    if 3500 <= n_seqs <= 4500:
        print(f"SANITY PASS: {n_seqs} upstream sequences (expected 3500-4500)")
        n_sanity_pass += 1
    else:
        print(f"SANITY FAIL: {n_seqs} upstream sequences (expected 3500-4500)")
        n_sanity_fail += 1

    # 3. All sequences exactly 200 bp
    from Bio import SeqIO as _SeqIO
    all_200 = True
    count_checked = 0
    for rec in _SeqIO.parse(result['all_upstream_fasta'], 'fasta'):
        count_checked += 1
        if len(rec.seq) != PARAMS['upstream_length']:
            print(f"SANITY FAIL: {rec.id} has length {len(rec.seq)}, expected {PARAMS['upstream_length']}")
            all_200 = False
            n_sanity_fail += 1
            break
    if all_200:
        print(f"SANITY PASS: All {count_checked} sequences are exactly {PARAMS['upstream_length']} bp")
        n_sanity_pass += 1

    # 4. Circular wrapping test - extract at origin
    h37rv_rec = load_genome(os.path.join(GENOME_DIR, 'H37Rv.fasta'))
    genome_seq = str(h37rv_rec.seq).upper()
    genome_len = len(genome_seq)

    # Test wrapping at origin (start near 0)
    wrap_origin = extract_circular(genome_seq, genome_len - 50, 100, genome_len)
    if len(wrap_origin) == 100:
        expected = genome_seq[-50:] + genome_seq[:50]
        if wrap_origin == expected:
            print("SANITY PASS: Circular wrapping at origin correct")
            n_sanity_pass += 1
        else:
            print("SANITY FAIL: Circular wrapping at origin produced wrong sequence")
            n_sanity_fail += 1
    else:
        print(f"SANITY FAIL: Origin wrap length {len(wrap_origin)}, expected 100")
        n_sanity_fail += 1

    # Test wrapping at terminus (position near genome_len)
    wrap_term = extract_circular(genome_seq, genome_len - 10, 20, genome_len)
    if len(wrap_term) == 20:
        expected = genome_seq[-10:] + genome_seq[:10]
        if wrap_term == expected:
            print("SANITY PASS: Circular wrapping at terminus correct")
            n_sanity_pass += 1
        else:
            print("SANITY FAIL: Circular wrapping at terminus produced wrong sequence")
            n_sanity_fail += 1
    else:
        print(f"SANITY FAIL: Terminus wrap length {len(wrap_term)}, expected 20")
        n_sanity_fail += 1

    # 5. Known operator FASTA parseable
    try:
        op_rec = _SeqIO.read(result['known_operator_fasta'], 'fasta')
        assert len(op_rec.seq) == PARAMS['operator_length']
        print(f"SANITY PASS: Known operator FASTA correct length ({len(op_rec.seq)} bp)")
        n_sanity_pass += 1
    except Exception as e:
        print(f"SANITY FAIL: Known operator FASTA error: {e}")
        n_sanity_fail += 1

    # --- SCIENTIFIC EXPECTATIONS (warn only) ---

    # 6. Known operator region should appear in upstream sequences
    op_seq = PARAMS['operator_sequence'][:25]  # First 25 bp as probe
    found_operator = False
    for rec in _SeqIO.parse(result['all_upstream_fasta'], 'fasta'):
        if op_seq in str(rec.seq).upper() or reverse_complement(op_seq) in str(rec.seq).upper():
            found_operator = True
            print(f"SCIENCE EXPECTED: Known operator found in upstream of {rec.id}")
            n_science_expected += 1
            break
    if not found_operator:
        print("SCIENCE WARN: Known operator NOT found in any upstream region")
        n_science_unexpected += 1

    # 7. MEME input has 2-3 sequences
    n_meme = result['n_meme_seqs']
    if 2 <= n_meme <= 3:
        print(f"SCIENCE EXPECTED: MEME input has {n_meme} ortholog sequences")
        n_science_expected += 1
    else:
        print(f"SCIENCE WARN: MEME input has {n_meme} sequences (expected 2-3)")
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
