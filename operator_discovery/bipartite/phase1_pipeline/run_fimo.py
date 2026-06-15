"""
phase1_pipeline/run_fimo.py — Run FIMO motif scanning or generate mock output.

- Checks shutil.which('fimo') for MEME suite availability.
- If present: run FIMO with --thresh 1e-4, parse fimo.tsv, annotate nearest gene,
  filter intergenic hits, rank by p-value.
- If missing: generate_mock_fimo_output() with 10 plausible 25-bp hits.
- Writes predicted_sites.csv
"""

import sys
import os
import shutil
import subprocess
import csv
import re

sys.path.insert(0, '/Users/aayanalwani/tb project/mce3r_stochastic')

# Ensure conda env bin is on PATH so shutil.which() can find fimo
_conda_bin = os.path.join(sys.prefix, 'bin')
if _conda_bin not in os.environ.get('PATH', ''):
    os.environ['PATH'] = _conda_bin + os.pathsep + os.environ.get('PATH', '')

from config.parameters import PARAMS

from Bio import SeqIO

PROJECT_ROOT = '/Users/aayanalwani/tb project/mce3r_stochastic'
GENOME_DIR = os.path.join(PROJECT_ROOT, 'data', 'genomes')
SEQ_DIR = os.path.join(PROJECT_ROOT, 'data', 'sequences')
RESULTS_DIR = os.path.join(PROJECT_ROOT, 'results', 'phase1')


def ensure_dirs():
    os.makedirs(RESULTS_DIR, exist_ok=True)


def parse_gff_for_annotation(gff_path):
    """
    Parse GFF to build gene annotation index.
    Returns list of (start0, end0, strand, locus_tag, gene_name).
    Coords are 0-based.
    """
    genes = []
    with open(gff_path) as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue
            if parts[2] not in ('gene', 'CDS'):
                continue
            start0 = int(parts[3]) - 1
            end0 = int(parts[4])
            strand = parts[6]
            attrs = parts[8]
            locus_tag = ''
            gene_name = ''
            m = re.search(r'locus_tag=([^;]+)', attrs)
            if m:
                locus_tag = m.group(1)
            m = re.search(r'gene=([^;]+)', attrs)
            if m:
                gene_name = m.group(1)
            m = re.search(r'old_locus_tag=([^;]+)', attrs)
            old_lt = m.group(1) if m else ''
            genes.append((start0, end0, strand, locus_tag, gene_name, old_lt))
    return genes


def find_nearest_gene(pos, genes):
    """Find the nearest gene to a genomic position."""
    best_dist = float('inf')
    best_gene = ('', '', '')
    for start0, end0, strand, locus_tag, gene_name, old_lt in genes:
        # Distance to gene
        if start0 <= pos <= end0:
            dist = 0
        else:
            dist = min(abs(pos - start0), abs(pos - end0))
        if dist < best_dist:
            best_dist = dist
            tag = locus_tag or old_lt or gene_name
            best_gene = (tag, gene_name, strand)
    return best_gene, best_dist


def is_intergenic(pos, genes):
    """Check if a position falls in an intergenic region."""
    for start0, end0, strand, *_ in genes:
        if start0 <= pos <= end0:
            return False
    return True


def run_real_fimo(meme_txt, input_fasta, output_dir):
    """Run FIMO with the given motif file against the genome."""
    import shutil as _shutil
    if os.path.exists(output_dir):
        _shutil.rmtree(output_dir)
    os.makedirs(output_dir, exist_ok=True)
    bg_path = os.path.join(RESULTS_DIR, 'background.model')
    cmd = [
        'fimo',
        '--thresh', str(PARAMS['fimo_thresh']),
        '--oc', output_dir,
    ]
    if os.path.exists(bg_path):
        cmd.extend(['--bfile', bg_path])
    cmd.extend([meme_txt, input_fasta])
    print(f"  Running FIMO: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
    if result.returncode != 0:
        print(f"  FIMO stderr: {result.stderr}")
        raise RuntimeError(f"FIMO failed with return code {result.returncode}")
    print(f"  FIMO completed, output in {output_dir}")
    return os.path.join(output_dir, 'fimo.tsv')


def parse_fimo_tsv(fimo_tsv_path):
    """Parse FIMO TSV output into list of dicts."""
    hits = []
    with open(fimo_tsv_path) as f:
        reader = csv.DictReader((line for line in f if not line.startswith('#')),
                                delimiter='\t')
        for row in reader:
            if not row.get('start'):
                continue
            try:
                hits.append({
                    'motif_id': row.get('motif_id', ''),
                    'sequence_name': row.get('sequence_name', ''),
                    'start': int(row['start']),
                    'stop': int(row['stop']),
                    'strand': row.get('strand', '+'),
                    'score': float(row.get('score', 0)),
                    'p_value': float(row.get('p-value', 1)),
                    'q_value': float(row.get('q-value', 1)),
                    'matched_sequence': row.get('matched_sequence', ''),
                })
            except (ValueError, KeyError):
                continue
    return hits


def generate_mock_fimo_output():
    """
    Generate mock FIMO output with 10 plausible 25-bp hits.
    Includes the known operator site and other biologically plausible locations.
    """
    output_dir = os.path.join(RESULTS_DIR, 'fimo_output')
    os.makedirs(output_dir, exist_ok=True)

    # The operator sequence strong site (first 25 bp)
    strong_site = PARAMS['operator_sequence'][:25]
    # Known operator region
    op_start = PARAMS['operator_region_h37rv_start']
    op_end = PARAMS['operator_region_h37rv_end']

    # Generate 10 mock hits with plausible properties
    # Hit 1: Known operator site (should be best hit)
    mock_hits = [
        {
            'motif_id': '1-MEME_mock_strong_site',
            'sequence_name': 'NC_000962.3',
            'start': op_start,
            'stop': op_start + 24,
            'strand': '+',
            'score': 25.0,
            'p_value': 1.2e-10,
            'q_value': 5.0e-8,
            'matched_sequence': strong_site,
        },
    ]

    # Additional plausible hits near known mce regulon genes
    # Rv1933c-Rv1935c regulon (lipid metabolism)
    import random
    random.seed(PARAMS['master_seed'])

    other_positions = [
        (1_933_000, 'GCCCCGCGCTATAGGACACTAGCAA', 1.5e-7, '+'),  # Near Rv1933c
        (2_160_000, 'GCCCCGCACTATAGGATTCTAGCAA', 3.2e-7, '+'),  # Near mce3 operon
        (2_207_500, 'GCCCCGCGCTATAGGATACTAGCAA', 8.5e-11, '+'), # In operator region (overlap)
        (1_936_500, 'GCCCTGCGCTATAGGATCCTAGCGA', 5.1e-6, '+'),  # Near Rv1936
        (2_208_000, 'GCCCCGCACTATAGGTTACTAGCAA', 2.3e-7, '-'),  # Downstream of operator
        (3_500_000, 'GCCCAGCGCTATAGGATACTCGCAA', 8.7e-6, '+'),  # Distant potential site
        (1_000_000, 'GCTCCGCGCTATAGGATCCTAGCAA', 4.2e-5, '-'),  # Another potential site
        (750_000, 'GCCCCGCGCTAAAGGATACAAGCAA', 7.8e-5, '+'),    # Weaker match
        (4_000_000, 'GCCCCGCACTATCGGATACTAGCGA', 9.1e-5, '-'),  # Weak distant site
    ]

    for pos, seq, pval, strand in other_positions:
        mock_hits.append({
            'motif_id': '1-MEME_mock_strong_site',
            'sequence_name': 'NC_000962.3',
            'start': pos,
            'stop': pos + 24,
            'strand': strand,
            'score': -10 * __import__('math').log10(pval),
            'p_value': pval,
            'q_value': pval * 10,
            'matched_sequence': seq,
        })

    # Sort by p-value
    mock_hits.sort(key=lambda x: x['p_value'])

    # Write fimo.tsv
    fimo_tsv_path = os.path.join(output_dir, 'fimo.tsv')
    with open(fimo_tsv_path, 'w', newline='') as f:
        f.write("# FIMO mock output\n")
        writer = csv.DictWriter(f, fieldnames=[
            'motif_id', 'sequence_name', 'start', 'stop', 'strand',
            'score', 'p-value', 'q-value', 'matched_sequence'
        ], delimiter='\t')
        writer.writeheader()
        for hit in mock_hits:
            writer.writerow({
                'motif_id': hit['motif_id'],
                'sequence_name': hit['sequence_name'],
                'start': hit['start'],
                'stop': hit['stop'],
                'strand': hit['strand'],
                'score': f"{hit['score']:.4f}",
                'p-value': f"{hit['p_value']:.2e}",
                'q-value': f"{hit['q_value']:.2e}",
                'matched_sequence': hit['matched_sequence'],
            })

    print(f"  Mock FIMO output generated: {fimo_tsv_path} ({len(mock_hits)} hits)")
    return fimo_tsv_path, mock_hits


def annotate_and_filter_hits(hits, gff_path):
    """
    Annotate each hit with nearest gene and filter for intergenic hits.
    Returns annotated hits sorted by p-value.
    """
    genes = parse_gff_for_annotation(gff_path)
    annotated = []
    for hit in hits:
        pos = hit['start']
        nearest, dist = find_nearest_gene(pos, genes)
        hit['nearest_gene'] = nearest[0]
        hit['nearest_gene_name'] = nearest[1]
        hit['nearest_gene_strand'] = nearest[2]
        hit['distance_to_gene'] = dist
        hit['is_intergenic'] = is_intergenic(pos, genes)
        annotated.append(hit)

    # Sort by p-value
    annotated.sort(key=lambda x: x['p_value'])
    return annotated


def write_predicted_sites(hits, outpath):
    """Write predicted binding sites to CSV."""
    fieldnames = [
        'rank', 'motif_id', 'sequence_name', 'start', 'stop', 'strand',
        'score', 'p_value', 'q_value', 'matched_sequence',
        'nearest_gene', 'nearest_gene_name', 'distance_to_gene', 'is_intergenic'
    ]
    with open(outpath, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for i, hit in enumerate(hits, 1):
            writer.writerow({
                'rank': i,
                'motif_id': hit.get('motif_id', ''),
                'sequence_name': hit.get('sequence_name', ''),
                'start': hit['start'],
                'stop': hit['stop'],
                'strand': hit.get('strand', ''),
                'score': f"{hit.get('score', 0):.4f}",
                'p_value': f"{hit['p_value']:.2e}",
                'q_value': f"{hit.get('q_value', 0):.2e}",
                'matched_sequence': hit.get('matched_sequence', ''),
                'nearest_gene': hit.get('nearest_gene', ''),
                'nearest_gene_name': hit.get('nearest_gene_name', ''),
                'distance_to_gene': hit.get('distance_to_gene', ''),
                'is_intergenic': hit.get('is_intergenic', ''),
            })
    print(f"  Predicted sites written: {outpath} ({len(hits)} rows)")


def run_fimo_step():
    """
    Main FIMO step. Returns dict with output paths and metadata.
    """
    ensure_dirs()

    meme_txt = os.path.join(RESULTS_DIR, 'meme_output', 'meme.txt')
    h37rv_fasta = os.path.join(GENOME_DIR, 'H37Rv.fasta')
    gff_path = os.path.join(GENOME_DIR, 'H37Rv.gff')

    if not os.path.exists(meme_txt):
        raise FileNotFoundError(f"MEME output not found: {meme_txt}")

    mock_used = False
    if shutil.which('fimo'):
        print("  FIMO found, running real FIMO...")
        fimo_tsv = run_real_fimo(meme_txt, h37rv_fasta,
                                  os.path.join(RESULTS_DIR, 'fimo_output'))
        hits = parse_fimo_tsv(fimo_tsv)
    else:
        print("  FIMO NOT found, generating mock output...")
        fimo_tsv, hits = generate_mock_fimo_output()
        mock_used = True

    # Annotate hits with gene information
    if os.path.exists(gff_path):
        hits = annotate_and_filter_hits(hits, gff_path)
    else:
        print("  WARNING: GFF not found, skipping annotation")
        for hit in hits:
            hit['nearest_gene'] = ''
            hit['nearest_gene_name'] = ''
            hit['distance_to_gene'] = -1
            hit['is_intergenic'] = True

    # Write predicted sites
    predicted_sites_path = os.path.join(RESULTS_DIR, 'predicted_sites.csv')
    write_predicted_sites(hits, predicted_sites_path)

    return {
        'fimo_tsv': fimo_tsv,
        'predicted_sites_csv': predicted_sites_path,
        'n_hits': len(hits),
        'mock_used': mock_used,
    }


if __name__ == "__main__":
    print("=" * 60)
    print("run_fimo.py — Self-test")
    print("=" * 60)

    n_sanity_pass = 0
    n_sanity_fail = 0
    n_science_expected = 0
    n_science_unexpected = 0

    try:
        result = run_fimo_step()
    except Exception as e:
        import traceback
        traceback.print_exc()
        print(f"SANITY FAIL: Exception during FIMO step: {e}")
        sys.exit(1)

    # --- SANITY CHECKS (hard-fail) ---

    # 1. Predicted sites CSV exists
    csv_path = result['predicted_sites_csv']
    if os.path.exists(csv_path) and os.path.getsize(csv_path) > 50:
        print(f"SANITY PASS: predicted_sites.csv exists ({os.path.getsize(csv_path)} bytes)")
        n_sanity_pass += 1
    else:
        print(f"SANITY FAIL: predicted_sites.csv missing or empty")
        n_sanity_fail += 1

    # 2. CSV is parseable with correct columns
    try:
        with open(csv_path) as f:
            reader = csv.DictReader(f)
            rows = list(reader)
        required_cols = ['rank', 'start', 'stop', 'p_value', 'matched_sequence']
        for col in required_cols:
            assert col in reader.fieldnames, f"Missing column: {col}"
        print(f"SANITY PASS: CSV parseable with {len(rows)} rows, all required columns present")
        n_sanity_pass += 1
    except Exception as e:
        print(f"SANITY FAIL: CSV parse error: {e}")
        n_sanity_fail += 1

    # 3. Hits count > 0
    if result['n_hits'] > 0:
        print(f"SANITY PASS: {result['n_hits']} hits found")
        n_sanity_pass += 1
    else:
        print(f"SANITY FAIL: 0 hits found")
        n_sanity_fail += 1

    # 4. Rows sorted by p-value
    try:
        pvals = [float(r['p_value']) for r in rows]
        assert pvals == sorted(pvals), "Rows not sorted by p-value"
        print(f"SANITY PASS: Rows sorted by p-value (best={pvals[0]:.2e})")
        n_sanity_pass += 1
    except Exception as e:
        print(f"SANITY FAIL: Sort check error: {e}")
        n_sanity_fail += 1

    # 5. All matched sequences are 25 bp (mock) or 20-30 bp (real)
    try:
        lengths = [len(r['matched_sequence']) for r in rows]
        all_ok = all(20 <= l <= 30 for l in lengths)
        if all_ok:
            print(f"SANITY PASS: All matched sequences 20-30 bp (actual: {set(lengths)})")
            n_sanity_pass += 1
        else:
            print(f"SANITY FAIL: Matched sequence lengths out of range: {set(lengths)}")
            n_sanity_fail += 1
    except Exception as e:
        print(f"SANITY FAIL: Sequence length check error: {e}")
        n_sanity_fail += 1

    # --- SCIENTIFIC EXPECTATIONS (warn only) ---

    # 6. Known operator site present (near position 2207477)
    op_start = PARAMS['operator_region_h37rv_start']
    op_end = PARAMS['operator_region_h37rv_end']
    found_op = False
    for r in rows:
        pos = int(r['start'])
        if op_start - 500 <= pos <= op_end + 500:
            found_op = True
            break
    if found_op:
        print(f"SCIENCE EXPECTED: Known operator region hit found")
        n_science_expected += 1
    else:
        print(f"SCIENCE WARN: No hit near known operator region ({op_start}-{op_end})")
        n_science_unexpected += 1

    # 7. Top hit near operator
    try:
        top_pos = int(rows[0]['start'])
        if op_start - 1000 <= top_pos <= op_end + 1000:
            print(f"SCIENCE EXPECTED: Top hit at position {top_pos} near operator")
            n_science_expected += 1
        else:
            print(f"SCIENCE WARN: Top hit at position {top_pos}, not near operator")
            n_science_unexpected += 1
    except Exception:
        n_science_unexpected += 1

    # --- Summary ---
    print(f"\n{'=' * 50}")
    print(f"SANITY:  {n_sanity_pass} pass, {n_sanity_fail} fail")
    print(f"SCIENCE: {n_science_expected} expected, {n_science_unexpected} unexpected")
    print(f"Mock used: {result['mock_used']}")

    if n_sanity_fail > 0:
        print("STOPPING: Sanity check failed.")
        sys.exit(1)
    else:
        print("All sanity checks passed.")
        sys.exit(0)
