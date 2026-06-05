"""
phase1_pipeline/conservation_check.py — Check conservation of predicted binding sites.

For the top 50 predicted sites: sliding-window percent-identity search against
M. bovis and M. marinum genomes.

Conserved = >=80% identity over >=80% of motif length.
Writes conservation_status.csv.
"""

import sys
import os
import csv

sys.path.insert(0, '/Users/aayanalwani/tb project/mce3r_stochastic')
from config.parameters import PARAMS

from Bio import SeqIO
from Bio.Seq import Seq

PROJECT_ROOT = '/Users/aayanalwani/tb project/mce3r_stochastic'
GENOME_DIR = os.path.join(PROJECT_ROOT, 'data', 'genomes')
RESULTS_DIR = os.path.join(PROJECT_ROOT, 'results', 'phase1')


def load_genome_str(fasta_path):
    """Load a genome FASTA and return the uppercase sequence string."""
    rec = SeqIO.read(fasta_path, 'fasta')
    return str(rec.seq).upper()


def reverse_complement(seq_str):
    """Return reverse complement of a DNA string."""
    return str(Seq(seq_str).reverse_complement())


def sliding_window_search(query, target_genome, min_identity=0.80, min_coverage=0.80):
    """
    Sliding-window percent-identity search using numpy vectorization.
    Searches target_genome for the best match to query on both strands.

    Returns (best_identity, best_position, best_strand).
    """
    import numpy as np

    query_len = len(query)
    target_len = len(target_genome)

    # Convert to numpy byte arrays for fast comparison
    target_arr = np.frombuffer(target_genome.encode('ascii'), dtype=np.uint8)

    best_identity = 0.0
    best_position = -1
    best_strand = '+'

    for strand_label, q_str in [('+', query), ('-', reverse_complement(query))]:
        q_arr = np.frombuffer(q_str.encode('ascii'), dtype=np.uint8)
        # Count matches at each position using a rolling sum
        # Build a match array: 1 where target matches query base at offset
        n_positions = target_len - query_len + 1
        match_counts = np.zeros(n_positions, dtype=np.int32)
        for j in range(query_len):
            match_counts += (target_arr[j:j + n_positions] == q_arr[j]).astype(np.int32)
        identities = match_counts / query_len
        idx = np.argmax(identities)
        if identities[idx] > best_identity:
            best_identity = float(identities[idx])
            best_position = int(idx)
            best_strand = strand_label

    return best_identity, best_position, best_strand


def run_conservation_check():
    """
    Main conservation check routine.
    Returns dict with output path and summary stats.
    """
    os.makedirs(RESULTS_DIR, exist_ok=True)

    # Load predicted sites
    predicted_csv = os.path.join(RESULTS_DIR, 'predicted_sites.csv')
    if not os.path.exists(predicted_csv):
        raise FileNotFoundError(f"Predicted sites not found: {predicted_csv}")

    with open(predicted_csv) as f:
        reader = csv.DictReader(f)
        all_sites = list(reader)

    # Take top 50 by rank (already sorted by p-value)
    top_sites = all_sites[:50]
    print(f"  Checking conservation for {len(top_sites)} sites...")

    # Load comparison genomes
    bovis_fasta = os.path.join(GENOME_DIR, 'M_bovis.fasta')
    marinum_fasta = os.path.join(GENOME_DIR, 'M_marinum.fasta')

    bovis_genome = load_genome_str(bovis_fasta) if os.path.exists(bovis_fasta) else None
    marinum_genome = load_genome_str(marinum_fasta) if os.path.exists(marinum_fasta) else None

    if bovis_genome:
        print(f"  M. bovis genome loaded: {len(bovis_genome)} bp")
    else:
        print("  WARNING: M. bovis genome not available")

    if marinum_genome:
        print(f"  M. marinum genome loaded: {len(marinum_genome)} bp")
    else:
        print("  WARNING: M. marinum genome not available")

    # Run conservation checks
    results = []
    for i, site in enumerate(top_sites):
        query_seq = site['matched_sequence'].upper()
        if not query_seq or len(query_seq) < 10:
            results.append({
                'rank': site['rank'],
                'start': site['start'],
                'stop': site['stop'],
                'strand': site['strand'],
                'p_value': site['p_value'],
                'matched_sequence': query_seq,
                'nearest_gene': site.get('nearest_gene', ''),
                'bovis_identity': 0.0,
                'bovis_position': -1,
                'bovis_strand': '+',
                'bovis_conserved': False,
                'marinum_identity': 0.0,
                'marinum_position': -1,
                'marinum_strand': '+',
                'marinum_conserved': False,
                'both_conserved': False,
            })
            continue

        # Search M. bovis
        if bovis_genome:
            b_id, b_pos, b_str = sliding_window_search(query_seq, bovis_genome)
        else:
            b_id, b_pos, b_str = 0.0, -1, '+'

        # Search M. marinum
        if marinum_genome:
            m_id, m_pos, m_str = sliding_window_search(query_seq, marinum_genome)
        else:
            m_id, m_pos, m_str = 0.0, -1, '+'

        b_conserved = b_id >= 0.80
        m_conserved = m_id >= 0.80

        results.append({
            'rank': site['rank'],
            'start': site['start'],
            'stop': site['stop'],
            'strand': site['strand'],
            'p_value': site['p_value'],
            'matched_sequence': query_seq,
            'nearest_gene': site.get('nearest_gene', ''),
            'bovis_identity': round(b_id, 4),
            'bovis_position': b_pos,
            'bovis_strand': b_str,
            'bovis_conserved': b_conserved,
            'marinum_identity': round(m_id, 4),
            'marinum_position': m_pos,
            'marinum_strand': m_str,
            'marinum_conserved': m_conserved,
            'both_conserved': b_conserved and m_conserved,
        })

        if (i + 1) % 10 == 0:
            print(f"  Processed {i+1}/{len(top_sites)} sites...")

    # Write results
    output_path = os.path.join(RESULTS_DIR, 'conservation_status.csv')
    fieldnames = [
        'rank', 'start', 'stop', 'strand', 'p_value', 'matched_sequence',
        'nearest_gene', 'bovis_identity', 'bovis_position', 'bovis_strand',
        'bovis_conserved', 'marinum_identity', 'marinum_position',
        'marinum_strand', 'marinum_conserved', 'both_conserved'
    ]
    with open(output_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for r in results:
            writer.writerow(r)

    n_bovis_conserved = sum(1 for r in results if r['bovis_conserved'])
    n_marinum_conserved = sum(1 for r in results if r['marinum_conserved'])
    n_both_conserved = sum(1 for r in results if r['both_conserved'])

    print(f"  Conservation results written: {output_path}")
    print(f"  Conserved in M. bovis: {n_bovis_conserved}/{len(results)}")
    print(f"  Conserved in M. marinum: {n_marinum_conserved}/{len(results)}")
    print(f"  Conserved in both: {n_both_conserved}/{len(results)}")

    return {
        'conservation_csv': output_path,
        'n_sites_checked': len(results),
        'n_bovis_conserved': n_bovis_conserved,
        'n_marinum_conserved': n_marinum_conserved,
        'n_both_conserved': n_both_conserved,
    }


if __name__ == "__main__":
    print("=" * 60)
    print("conservation_check.py — Self-test")
    print("=" * 60)

    n_sanity_pass = 0
    n_sanity_fail = 0
    n_science_expected = 0
    n_science_unexpected = 0

    try:
        result = run_conservation_check()
    except Exception as e:
        import traceback
        traceback.print_exc()
        print(f"SANITY FAIL: Exception during conservation check: {e}")
        sys.exit(1)

    # --- SANITY CHECKS (hard-fail) ---

    # 1. Output CSV exists
    csv_path = result['conservation_csv']
    if os.path.exists(csv_path) and os.path.getsize(csv_path) > 50:
        print(f"SANITY PASS: conservation_status.csv exists ({os.path.getsize(csv_path)} bytes)")
        n_sanity_pass += 1
    else:
        print(f"SANITY FAIL: conservation_status.csv missing or empty")
        n_sanity_fail += 1

    # 2. CSV parseable with correct columns
    try:
        with open(csv_path) as f:
            reader = csv.DictReader(f)
            rows = list(reader)
        required_cols = ['rank', 'start', 'bovis_identity', 'bovis_conserved',
                         'marinum_identity', 'marinum_conserved', 'both_conserved']
        for col in required_cols:
            assert col in reader.fieldnames, f"Missing column: {col}"
        print(f"SANITY PASS: CSV has correct columns and {len(rows)} rows")
        n_sanity_pass += 1
    except Exception as e:
        print(f"SANITY FAIL: CSV column check error: {e}")
        n_sanity_fail += 1

    # 3. Row count matches sites checked
    if len(rows) == result['n_sites_checked']:
        print(f"SANITY PASS: Row count matches ({len(rows)} rows)")
        n_sanity_pass += 1
    else:
        print(f"SANITY FAIL: Row count mismatch: {len(rows)} vs {result['n_sites_checked']}")
        n_sanity_fail += 1

    # 4. Identity values in [0, 1]
    try:
        for row in rows:
            b_id = float(row['bovis_identity'])
            m_id = float(row['marinum_identity'])
            assert 0.0 <= b_id <= 1.0, f"bovis_identity {b_id} out of range"
            assert 0.0 <= m_id <= 1.0, f"marinum_identity {m_id} out of range"
        print(f"SANITY PASS: All identity values in [0, 1]")
        n_sanity_pass += 1
    except Exception as e:
        print(f"SANITY FAIL: Identity range error: {e}")
        n_sanity_fail += 1

    # 5. Boolean columns are valid
    try:
        for row in rows:
            assert row['bovis_conserved'] in ('True', 'False')
            assert row['marinum_conserved'] in ('True', 'False')
            assert row['both_conserved'] in ('True', 'False')
        print(f"SANITY PASS: Boolean columns valid")
        n_sanity_pass += 1
    except Exception as e:
        print(f"SANITY FAIL: Boolean column error: {e}")
        n_sanity_fail += 1

    # --- SCIENTIFIC EXPECTATIONS (warn only) ---

    # 6. At least some sites conserved in M. bovis (closely related)
    if result['n_bovis_conserved'] > 0:
        print(f"SCIENCE EXPECTED: {result['n_bovis_conserved']} sites conserved in M. bovis")
        n_science_expected += 1
    else:
        print(f"SCIENCE WARN: No sites conserved in M. bovis")
        n_science_unexpected += 1

    # 7. Fewer sites conserved in M. marinum (more distant)
    if result['n_marinum_conserved'] <= result['n_bovis_conserved']:
        print(f"SCIENCE EXPECTED: Fewer sites conserved in M. marinum ({result['n_marinum_conserved']}) "
              f"than M. bovis ({result['n_bovis_conserved']})")
        n_science_expected += 1
    else:
        print(f"SCIENCE WARN: More sites conserved in M. marinum ({result['n_marinum_conserved']}) "
              f"than M. bovis ({result['n_bovis_conserved']})")
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
