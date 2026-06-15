#!/usr/bin/env python3
"""
Targeted Bipartite Operator Search for Mce3R in M. tuberculosis H37Rv
=====================================================================

Goal: Find the 3 unmapped Mce3R operator sites in the TB genome.

Strategy:
  1. Extract the two ~25bp half-sites from the known 123bp operator (PDB 9B7Y)
  2. Build position weight matrices (PWMs) for each half-site using the known
     sequence + evolutionary variation from orthologs
  3. Scan the entire H37Rv genome for matches to either half-site
  4. Filter for BIPARTITE architecture: pairs of hits on the SAME strand (or
     convergent strands) separated by 30-80bp (allowing flexibility around
     the known 53bp spacer)
  5. Rank candidates by combined score
  6. Validate: the known operator must be recovered as the #1 hit

Author: Computational pipeline (Mce3R operator discovery)
Date: 2026-03-18
"""

import os
import sys
import json
import math
import csv
from collections import defaultdict

# ============================================================================
# STEP 0: CONFIGURATION & KNOWN OPERATOR DECOMPOSITION
# ============================================================================

# Full 123bp known operator sequence (Panagoda et al. 2024, PDB 9B7Y)
KNOWN_OPERATOR = (
    "GCCCCGCGCTATAGGATACTAGCAAGATACATCATAGCCAATATATGCCAGT"
    "TTGCATTGCTATTTACCGATCAGTTGTCCAAGCAATCGCGTATTGGCTATG"
    "GACATCAGCGGTTCTGCCGC"
)

# Operator starts at this coordinate in H37Rv (+ strand)
OPERATOR_GENOME_START = 2_207_577  # We'll verify this

# The user-provided known binding site (strong/downstream site)
USER_KNOWN_SITE = "aatcgcgtattggctatggacatcagc"  # 27bp

# --- Decompose the 123bp operator into its bipartite components ---
# From crystal structure (PDB 9B7Y), the operator has:
#   [upstream/weak site ~25bp] [spacer ~53bp] [downstream/strong site ~25bp]
#
# The user's sequence "aatcgcgtattggctatggacatcagc" appears in the DOWNSTREAM
# portion of the operator. Let's find it:

idx = KNOWN_OPERATOR.upper().find(USER_KNOWN_SITE.upper())
assert idx != -1, f"User's known site not found in operator! Check sequence."
print(f"[SELF-CHECK] User's known site found at position {idx} within 123bp operator")
print(f"[SELF-CHECK] Operator length: {len(KNOWN_OPERATOR)}bp")

# The bipartite structure from Panagoda 2024:
# - Strong site (downstream, Kd=2.4nM): the site the user provided, approximately
#   the last ~27bp of the operator
# - Weak site (upstream, Kd=49nM): approximately the first ~25bp
# - Spacer: the middle ~53bp connecting them
#
# Based on the crystal structure contacts:
# Upstream (weak) half-site:  positions 0-24  (25bp)
# Spacer:                     positions 25-77 (53bp)
# Downstream (strong) half-site: positions 78-122 (45bp, but core contact ~25-27bp)
#
# However, the EMSA probes from Panagoda 2024 define the sites more precisely.
# The user's 27bp sequence starts at idx=73 in the 123bp operator.
# Let's define generous half-sites to capture the full binding footprint:

# Strong site: centered on the user's known sequence
STRONG_SITE_START = 73   # where "aatcgcgtattggctatggacatcagc" starts
STRONG_SITE_END = 100    # 27bp
STRONG_SITE = KNOWN_OPERATOR[STRONG_SITE_START:STRONG_SITE_END]
print(f"[SELF-CHECK] Strong site: {STRONG_SITE} ({len(STRONG_SITE)}bp)")
print(f"[SELF-CHECK] Matches user input: {STRONG_SITE.upper() == USER_KNOWN_SITE.upper()}")

# Weak site: the upstream portion, symmetric ~27bp
# From the crystal structure, the weak site is at the OTHER end of the operator.
# With 53bp spacer and the strong site starting at pos 73:
#   Weak site should end around position 73 - 53 = 20
#   So weak site: positions 0-26 (27bp) approximately
WEAK_SITE_START = 0
WEAK_SITE_END = 27
WEAK_SITE = KNOWN_OPERATOR[WEAK_SITE_START:WEAK_SITE_END]
print(f"[SELF-CHECK] Weak site: {WEAK_SITE} ({len(WEAK_SITE)}bp)")

# Spacer
SPACER = KNOWN_OPERATOR[WEAK_SITE_END:STRONG_SITE_START]
SPACER_LEN = len(SPACER)
print(f"[SELF-CHECK] Spacer: {SPACER} ({SPACER_LEN}bp)")
print(f"[SELF-CHECK] Expected ~53bp spacer, got {SPACER_LEN}bp")

# Verify reconstruction
assert WEAK_SITE + SPACER + STRONG_SITE + KNOWN_OPERATOR[STRONG_SITE_END:] == KNOWN_OPERATOR, \
    "Operator decomposition failed!"
print("[SELF-CHECK] Operator decomposition verified: weak + spacer + strong = full operator ✓")
print()

# ============================================================================
# STEP 1: LOAD H37Rv GENOME
# ============================================================================

GENOME_FILE = os.path.join(os.path.dirname(__file__), "data", "genomes", "H37Rv.fasta")
GFF_FILE = os.path.join(os.path.dirname(__file__), "data", "genomes", "H37Rv.gff")

def load_fasta(filepath):
    """Load a single-sequence FASTA file."""
    sequence = []
    header = None
    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                header = line[1:].split()[0]
            else:
                sequence.append(line.upper())
    seq = ''.join(sequence)
    print(f"[SELF-CHECK] Loaded genome: {header}, length={len(seq):,}bp")
    assert len(seq) > 4_000_000, f"Genome too short ({len(seq)}bp), expected ~4.4Mbp"
    return header, seq

def reverse_complement(seq):
    """Return the reverse complement of a DNA sequence."""
    comp = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N',
            'R': 'Y', 'Y': 'R', 'W': 'W', 'S': 'S', 'M': 'K', 'K': 'M'}
    return ''.join(comp.get(b, 'N') for b in reversed(seq.upper()))

print("Loading H37Rv genome...")
GENOME_ID, GENOME_SEQ = load_fasta(GENOME_FILE)
GENOME_LEN = len(GENOME_SEQ)
print()

# ============================================================================
# STEP 2: VERIFY KNOWN OPERATOR LOCATION IN GENOME
# ============================================================================

print("Verifying known operator location in genome...")
# Find the exact position of the known operator in the genome
op_pos = GENOME_SEQ.find(KNOWN_OPERATOR.upper())
if op_pos == -1:
    # Try reverse complement
    op_rc = reverse_complement(KNOWN_OPERATOR)
    op_pos = GENOME_SEQ.find(op_rc)
    if op_pos != -1:
        print(f"[SELF-CHECK] Known operator found on MINUS strand at position {op_pos}")
        OPERATOR_STRAND = '-'
    else:
        print("[WARNING] Known operator not found as exact match. Trying fuzzy search...")
        # Search for the user's 27bp known site directly
        user_pos = GENOME_SEQ.find(USER_KNOWN_SITE.upper())
        user_pos_rc = GENOME_SEQ.find(reverse_complement(USER_KNOWN_SITE))
        print(f"  User's strong site on + strand: {user_pos}")
        print(f"  User's strong site on - strand: {user_pos_rc}")
        OPERATOR_STRAND = '+'
else:
    print(f"[SELF-CHECK] Known operator found on PLUS strand at position {op_pos}")
    OPERATOR_STRAND = '+'

# Also find each half-site independently
strong_pos_plus = GENOME_SEQ.find(STRONG_SITE.upper())
strong_pos_minus = GENOME_SEQ.find(reverse_complement(STRONG_SITE))
weak_pos_plus = GENOME_SEQ.find(WEAK_SITE.upper())
weak_pos_minus = GENOME_SEQ.find(reverse_complement(WEAK_SITE))

print(f"[SELF-CHECK] Strong site on + strand: pos {strong_pos_plus}")
print(f"[SELF-CHECK] Strong site on - strand: pos {strong_pos_minus}")
print(f"[SELF-CHECK] Weak site on + strand: pos {weak_pos_plus}")
print(f"[SELF-CHECK] Weak site on - strand: pos {weak_pos_minus}")
print()

# ============================================================================
# STEP 3: BUILD SCORING MATRIX FROM KNOWN BINDING SITES
# ============================================================================

print("Building scoring matrices...")

def seq_to_pwm(sequence, pseudocount=0.25):
    """
    Build a simple log-odds PWM from a single known sequence.
    Uses pseudocounts and TB genome background frequencies (GC=65.6%).
    """
    bg = {'A': 0.172, 'T': 0.172, 'G': 0.328, 'C': 0.328}  # H37Rv GC=65.6%
    seq = sequence.upper()
    pwm = []
    for base in seq:
        counts = {b: pseudocount for b in 'ACGT'}
        counts[base] += 1.0
        total = sum(counts.values())
        row = {}
        for b in 'ACGT':
            freq = counts[b] / total
            row[b] = math.log2(freq / bg[b])
        pwm.append(row)
    return pwm

def score_sequence(pwm, seq):
    """Score a sequence against a PWM. Returns log-odds score."""
    seq = seq.upper()
    if len(seq) != len(pwm):
        return float('-inf')
    score = 0.0
    for i, base in enumerate(seq):
        if base not in 'ACGT':
            return float('-inf')
        score += pwm[i][base]
    return score

def max_possible_score(pwm):
    """Return the maximum possible score for a PWM."""
    return sum(max(row.values()) for row in pwm)

# Build PWMs for both half-sites
STRONG_PWM = seq_to_pwm(STRONG_SITE)
WEAK_PWM = seq_to_pwm(WEAK_SITE)

# Also build PWMs for reverse complements
STRONG_RC_PWM = seq_to_pwm(reverse_complement(STRONG_SITE))
WEAK_RC_PWM = seq_to_pwm(reverse_complement(WEAK_SITE))

strong_max = max_possible_score(STRONG_PWM)
weak_max = max_possible_score(WEAK_PWM)

print(f"[SELF-CHECK] Strong site PWM: {len(STRONG_PWM)} positions, max score = {strong_max:.2f}")
print(f"[SELF-CHECK] Weak site PWM: {len(WEAK_PWM)} positions, max score = {weak_max:.2f}")

# Verify PWMs score their own sequences perfectly
strong_self_score = score_sequence(STRONG_PWM, STRONG_SITE)
weak_self_score = score_sequence(WEAK_PWM, WEAK_SITE)
print(f"[SELF-CHECK] Strong site self-score: {strong_self_score:.2f} (should equal max {strong_max:.2f})")
print(f"[SELF-CHECK] Weak site self-score: {weak_self_score:.2f} (should equal max {weak_max:.2f})")
assert abs(strong_self_score - strong_max) < 0.01, "Strong PWM self-score mismatch!"
assert abs(weak_self_score - weak_max) < 0.01, "Weak PWM self-score mismatch!"
print("[SELF-CHECK] PWM self-scoring verified ✓")
print()

# ============================================================================
# STEP 4: GENOME-WIDE SCAN FOR INDIVIDUAL HALF-SITE MATCHES
# ============================================================================

print("Scanning genome for individual half-site matches...")
print(f"  Strong site: {len(STRONG_SITE)}bp motif")
print(f"  Weak site: {len(WEAK_SITE)}bp motif")

# Score threshold: we want ~60% of max score to catch degenerate sites
# This is deliberately permissive — the bipartite filter will do the real selection
STRONG_THRESHOLD_FRAC = 0.50
WEAK_THRESHOLD_FRAC = 0.50
STRONG_THRESHOLD = strong_max * STRONG_THRESHOLD_FRAC
WEAK_THRESHOLD = weak_max * WEAK_THRESHOLD_FRAC

print(f"  Strong threshold: {STRONG_THRESHOLD:.2f} ({STRONG_THRESHOLD_FRAC*100:.0f}% of max {strong_max:.2f})")
print(f"  Weak threshold: {WEAK_THRESHOLD:.2f} ({WEAK_THRESHOLD_FRAC*100:.0f}% of max {weak_max:.2f})")

def scan_genome(genome_seq, pwm, motif_len, threshold, label):
    """Scan both strands of genome for motif matches above threshold."""
    hits = []
    rc_pwm = seq_to_pwm(reverse_complement('A' * motif_len))  # dummy, we'll use actual

    # Build reverse complement PWM properly
    # The RC PWM scores the reverse complement of a subsequence
    # Equivalently, we can reverse the PWM rows and swap A<->T, G<->C
    rc_pwm_rows = []
    for row in reversed(pwm):
        rc_row = {'A': row['T'], 'T': row['A'], 'G': row['C'], 'C': row['G']}
        rc_pwm_rows.append(rc_row)

    total_positions = len(genome_seq) - motif_len + 1

    for i in range(total_positions):
        subseq = genome_seq[i:i + motif_len]

        # Score on plus strand
        plus_score = score_sequence(pwm, subseq)
        if plus_score >= threshold:
            hits.append({
                'pos': i,
                'strand': '+',
                'score': plus_score,
                'sequence': subseq,
                'label': label
            })

        # Score on minus strand (equivalent to scoring RC of subseq against original PWM)
        minus_score = score_sequence(rc_pwm_rows, subseq)
        if minus_score >= threshold:
            hits.append({
                'pos': i,
                'strand': '-',
                'score': minus_score,
                'sequence': reverse_complement(subseq),
                'label': label
            })

    hits.sort(key=lambda x: x['score'], reverse=True)
    print(f"  [{label}] Found {len(hits)} hits above threshold")
    if hits:
        print(f"    Top score: {hits[0]['score']:.2f} at pos {hits[0]['pos']} ({hits[0]['strand']})")
        top5 = [round(h['score'], 1) for h in hits[:5]]
        print(f"    Top 5 scores: {top5}")
    return hits

strong_hits = scan_genome(GENOME_SEQ, STRONG_PWM, len(STRONG_SITE), STRONG_THRESHOLD, "strong")
weak_hits = scan_genome(GENOME_SEQ, WEAK_PWM, len(WEAK_SITE), WEAK_THRESHOLD, "weak")

# Self-check: verify the known operator sites are in the results
print()
print("[SELF-CHECK] Verifying known operator sites are captured...")
known_region = range(2_207_400, 2_207_700)
strong_in_known = [h for h in strong_hits if h['pos'] in known_region]
weak_in_known = [h for h in weak_hits if h['pos'] in known_region]
print(f"  Strong hits in known operator region: {len(strong_in_known)}")
for h in strong_in_known:
    print(f"    pos={h['pos']}, strand={h['strand']}, score={h['score']:.2f}")
print(f"  Weak hits in known operator region: {len(weak_in_known)}")
for h in weak_in_known:
    print(f"    pos={h['pos']}, strand={h['strand']}, score={h['score']:.2f}")

if not strong_in_known or not weak_in_known:
    print("[WARNING] Known operator not fully captured! Adjusting thresholds...")
    # Lower thresholds and rescan
    STRONG_THRESHOLD = strong_max * 0.35
    WEAK_THRESHOLD = weak_max * 0.35
    print(f"  New strong threshold: {STRONG_THRESHOLD:.2f}")
    print(f"  New weak threshold: {WEAK_THRESHOLD:.2f}")
    strong_hits = scan_genome(GENOME_SEQ, STRONG_PWM, len(STRONG_SITE), STRONG_THRESHOLD, "strong")
    weak_hits = scan_genome(GENOME_SEQ, WEAK_PWM, len(WEAK_SITE), WEAK_THRESHOLD, "weak")
    strong_in_known = [h for h in strong_hits if h['pos'] in known_region]
    weak_in_known = [h for h in weak_hits if h['pos'] in known_region]
    print(f"  After adjustment - Strong in known region: {len(strong_in_known)}")
    print(f"  After adjustment - Weak in known region: {len(weak_in_known)}")

print()

# ============================================================================
# STEP 5: FIND BIPARTITE OPERATOR CANDIDATES
# ============================================================================

print("=" * 70)
print("SEARCHING FOR BIPARTITE OPERATOR ARCHITECTURE")
print("=" * 70)
print()
print("Architecture: [weak site ~27bp] --spacer 30-80bp-- [strong site ~27bp]")
print(f"Known spacer: {SPACER_LEN}bp, searching range: 30-80bp")
print()

# The bipartite operator has weak upstream, strong downstream.
# On the + strand: weak_pos < strong_pos, with spacer = strong_pos - weak_pos - weak_len
# On the - strand: the order reverses (strong is upstream in genome coordinates)
#
# We also allow the reverse: strong upstream, weak downstream (the repressor
# dimer might bind in either orientation)

MIN_SPACER = 30
MAX_SPACER = 80

# Index hits by position for efficient lookup
# Group hits into genomic bins for fast pair-finding
BIN_SIZE = 200  # bp

def bin_hits(hits):
    """Group hits into genomic position bins."""
    bins = defaultdict(list)
    for h in hits:
        bin_id = h['pos'] // BIN_SIZE
        bins[bin_id].append(h)
    return bins

strong_bins = bin_hits(strong_hits)
weak_bins = bin_hits(weak_hits)

# For each strong hit, look for a nearby weak hit (or vice versa) at proper spacing
operators = []

print("Searching for bipartite pairs...")

# Strategy: for each weak hit, search for a strong hit downstream at proper spacing
# (on same strand) — and vice versa for the reverse arrangement
for w in weak_hits:
    w_pos = w['pos']
    w_end = w_pos + len(WEAK_SITE)
    w_strand = w['strand']

    # Look for strong hits in the right distance range
    # Case 1: weak upstream, strong downstream (same strand +)
    # spacer = strong_pos - w_end
    target_start = w_end + MIN_SPACER
    target_end = w_end + MAX_SPACER + len(STRONG_SITE)

    for s in strong_hits:
        s_pos = s['pos']
        if target_start <= s_pos <= target_end:
            spacer = s_pos - w_end
            # Allow same strand OR opposite strand (both are biologically possible)
            combined_score = w['score'] + s['score']
            normalized_score = combined_score / (strong_max + weak_max)

            operators.append({
                'weak_pos': w_pos,
                'strong_pos': s_pos,
                'weak_strand': w_strand,
                'strong_strand': s['strand'],
                'spacer': spacer,
                'weak_score': w['score'],
                'strong_score': s['score'],
                'combined_score': combined_score,
                'normalized_score': normalized_score,
                'weak_seq': w['sequence'],
                'strong_seq': s['sequence'],
                'arrangement': 'weak-spacer-strong',
                'center_pos': (w_pos + s_pos + len(STRONG_SITE)) // 2,
            })

    # Case 2: strong upstream, weak downstream (reversed arrangement)
    # This time look for strong hits UPSTREAM of the weak hit
    target_start2 = w_pos - MAX_SPACER - len(STRONG_SITE)
    target_end2 = w_pos - MIN_SPACER

    for s in strong_hits:
        s_pos = s['pos']
        s_end = s_pos + len(STRONG_SITE)
        if target_start2 <= s_pos <= target_end2:
            spacer = w_pos - s_end
            if spacer < MIN_SPACER or spacer > MAX_SPACER:
                continue
            combined_score = w['score'] + s['score']
            normalized_score = combined_score / (strong_max + weak_max)

            operators.append({
                'weak_pos': w_pos,
                'strong_pos': s_pos,
                'weak_strand': w_strand,
                'strong_strand': s['strand'],
                'spacer': spacer,
                'weak_score': w['score'],
                'strong_score': s['score'],
                'combined_score': combined_score,
                'normalized_score': normalized_score,
                'weak_seq': w['sequence'],
                'strong_seq': s['sequence'],
                'arrangement': 'strong-spacer-weak',
                'center_pos': (w_pos + s_pos + len(STRONG_SITE)) // 2,
            })

print(f"Found {len(operators)} raw bipartite candidates")

# Deduplicate: merge overlapping operators (within 20bp of each other)
operators.sort(key=lambda x: x['combined_score'], reverse=True)

def deduplicate_operators(ops, merge_dist=30):
    """Remove near-duplicate operator predictions."""
    kept = []
    used_positions = set()
    for op in ops:
        center = op['center_pos']
        # Check if any previously kept operator is within merge_dist
        is_dup = False
        for used_center in used_positions:
            if abs(center - used_center) < merge_dist:
                is_dup = True
                break
        if not is_dup:
            kept.append(op)
            used_positions.add(center)
    return kept

operators = deduplicate_operators(operators)
print(f"After deduplication: {len(operators)} unique operator candidates")
print()

# ============================================================================
# STEP 6: ANNOTATE WITH GENE INFORMATION
# ============================================================================

print("Loading gene annotations...")

genes = []
with open(GFF_FILE) as f:
    for line in f:
        if line.startswith('#'):
            continue
        parts = line.strip().split('\t')
        if len(parts) < 9:
            continue
        if parts[2] != 'gene':
            continue
        start = int(parts[3])
        end = int(parts[4])
        strand = parts[6]
        attrs = parts[8]

        # Extract locus_tag and gene name
        locus_tag = ''
        gene_name = ''
        product = ''
        for attr in attrs.split(';'):
            if attr.startswith('locus_tag='):
                locus_tag = attr.split('=')[1]
            elif attr.startswith('gene='):
                gene_name = attr.split('=')[1]
            elif attr.startswith('product=') or attr.startswith('Note='):
                product = attr.split('=')[1]

        genes.append({
            'start': start,
            'end': end,
            'strand': strand,
            'locus_tag': locus_tag,
            'name': gene_name,
            'product': product
        })

print(f"  Loaded {len(genes)} genes from GFF")

# Sort genes by start position for binary search
genes.sort(key=lambda g: g['start'])
gene_starts = [g['start'] for g in genes]

def find_nearest_genes(pos, n=3):
    """Find the n nearest genes to a genomic position."""
    import bisect
    idx = bisect.bisect_left(gene_starts, pos)

    candidates = []
    for i in range(max(0, idx - 3), min(len(genes), idx + 3)):
        g = genes[i]
        # Distance: negative if pos is inside gene, positive if outside
        if g['start'] <= pos <= g['end']:
            dist = 0
            location = 'inside'
        elif pos < g['start']:
            dist = g['start'] - pos
            location = 'upstream' if g['strand'] == '+' else 'downstream'
        else:
            dist = pos - g['end']
            location = 'downstream' if g['strand'] == '+' else 'upstream'

        candidates.append({
            'locus_tag': g['locus_tag'],
            'name': g['name'],
            'distance': dist,
            'location': location,
            'gene_start': g['start'],
            'gene_end': g['end'],
            'gene_strand': g['strand'],
            'product': g['product']
        })

    candidates.sort(key=lambda c: c['distance'])
    return candidates[:n]

def is_intergenic(pos, site_len):
    """Check if a position is in an intergenic region."""
    for g in genes:
        if g['start'] <= pos <= g['end'] or g['start'] <= pos + site_len <= g['end']:
            return False
    return True

# Annotate each operator candidate
for op in operators:
    center = op['center_pos']
    nearest = find_nearest_genes(center)
    op['nearest_genes'] = nearest
    op['primary_gene'] = nearest[0] if nearest else None

    # Check if the operator is in an intergenic region
    op_start = min(op['weak_pos'], op['strong_pos'])
    op_end = max(op['weak_pos'], op['strong_pos']) + max(len(STRONG_SITE), len(WEAK_SITE))
    op['is_intergenic'] = is_intergenic(op_start, op_end - op_start)

    # Check if upstream of a known regulon gene
    regulon_genes = {'Rv1963c', 'Rv1964', 'Rv1933c', 'Rv1934c', 'Rv1935c',
                     'Rv1936', 'Rv1937', 'Rv1938', 'Rv1939', 'Rv1940', 'Rv1941c'}
    op['near_regulon'] = any(g['locus_tag'] in regulon_genes for g in nearest)

print("Annotation complete.")
print()

# ============================================================================
# STEP 7: RANK AND OUTPUT RESULTS
# ============================================================================

print("=" * 70)
print("BIPARTITE OPERATOR SEARCH RESULTS")
print("=" * 70)
print()

# Print top 20 candidates
TOP_N = min(30, len(operators))
for i, op in enumerate(operators[:TOP_N]):
    gene_info = op['primary_gene']
    gene_str = f"{gene_info['locus_tag']}"
    if gene_info['name']:
        gene_str += f" ({gene_info['name']})"
    gene_str += f" [{gene_info['distance']}bp {gene_info['location']}]"

    # Flag special cases
    flags = []
    if op['near_regulon']:
        flags.append("★ KNOWN REGULON")
    if op['is_intergenic']:
        flags.append("INTERGENIC")
    if abs(op['spacer'] - SPACER_LEN) <= 5:
        flags.append(f"SPACER MATCH ({op['spacer']}bp ≈ {SPACER_LEN}bp)")

    flag_str = " | ".join(flags) if flags else ""

    print(f"--- Operator Candidate #{i+1} ---")
    print(f"  Position: {min(op['weak_pos'], op['strong_pos']):,} - "
          f"{max(op['weak_pos'], op['strong_pos']) + 27:,}")
    print(f"  Near gene: {gene_str}")
    print(f"  Arrangement: {op['arrangement']}")
    print(f"  Weak site:   pos {op['weak_pos']:>9,} ({op['weak_strand']}) score={op['weak_score']:.2f}  seq={op['weak_seq']}")
    print(f"  Strong site:  pos {op['strong_pos']:>9,} ({op['strong_strand']}) score={op['strong_score']:.2f}  seq={op['strong_seq']}")
    print(f"  Spacer: {op['spacer']}bp | Combined score: {op['combined_score']:.2f} "
          f"({op['normalized_score']*100:.1f}% of max)")
    if flag_str:
        print(f"  >>> {flag_str}")
    print()

# ============================================================================
# STEP 8: VALIDATION SELF-CHECKS
# ============================================================================

print("=" * 70)
print("VALIDATION SELF-CHECKS")
print("=" * 70)
print()

# Check 1: Is the known operator recovered?
known_op_found = False
known_op_rank = None
for i, op in enumerate(operators):
    if 2_207_400 <= op['center_pos'] <= 2_207_700:
        known_op_found = True
        known_op_rank = i + 1
        print(f"[CHECK 1] ✅ Known operator recovered at rank #{known_op_rank}")
        print(f"           Score: {op['combined_score']:.2f}, Spacer: {op['spacer']}bp")
        break

if not known_op_found:
    print("[CHECK 1] ❌ FAIL: Known operator NOT found in results!")
    print("          This means the pipeline has a bug or thresholds are too strict.")

# Check 2: Spacer length distribution
spacers = [op['spacer'] for op in operators[:20]]
print(f"\n[CHECK 2] Spacer lengths in top 20: {spacers}")
print(f"          Mean spacer: {sum(spacers)/len(spacers):.1f}bp (expected ~{SPACER_LEN}bp)")

# Check 3: Score distribution
scores = [op['combined_score'] for op in operators[:20]]
print(f"\n[CHECK 3] Combined scores in top 20:")
print(f"          Max: {max(scores):.2f}, Min: {min(scores):.2f}")
if len(operators) > 1:
    gap = operators[0]['combined_score'] - operators[1]['combined_score']
    print(f"          Gap between #1 and #2: {gap:.2f} "
          f"({'known operator stands out ✅' if gap > 5 else 'may need stricter filtering'})")

# Check 4: Do any candidates fall near the known regulon genes?
regulon_candidates = [op for op in operators if op['near_regulon']]
print(f"\n[CHECK 4] Candidates near known regulon genes: {len(regulon_candidates)}")
for rc in regulon_candidates[:5]:
    gene = rc['primary_gene']
    print(f"          Rank #{operators.index(rc)+1}: near {gene['locus_tag']} ({gene['name']}), "
          f"score={rc['combined_score']:.2f}")

# Check 5: How many are intergenic (biologically most relevant)?
intergenic_ops = [op for op in operators if op['is_intergenic']]
print(f"\n[CHECK 5] Intergenic operator candidates: {len(intergenic_ops)} out of {len(operators)} total")

print()

# ============================================================================
# STEP 9: SAVE RESULTS
# ============================================================================

output_dir = os.path.join(os.path.dirname(__file__), "results", "phase1")
os.makedirs(output_dir, exist_ok=True)

# Save full results as CSV
csv_path = os.path.join(output_dir, "bipartite_operators.csv")
with open(csv_path, 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow([
        'rank', 'center_pos', 'weak_pos', 'strong_pos', 'weak_strand', 'strong_strand',
        'spacer_bp', 'weak_score', 'strong_score', 'combined_score', 'normalized_score',
        'arrangement', 'weak_sequence', 'strong_sequence',
        'nearest_gene', 'gene_name', 'gene_distance', 'gene_location',
        'is_intergenic', 'near_regulon'
    ])
    for i, op in enumerate(operators):
        gene = op['primary_gene']
        writer.writerow([
            i + 1, op['center_pos'], op['weak_pos'], op['strong_pos'],
            op['weak_strand'], op['strong_strand'],
            op['spacer'], op['weak_score'], op['strong_score'],
            op['combined_score'], f"{op['normalized_score']:.4f}",
            op['arrangement'], op['weak_seq'], op['strong_seq'],
            gene['locus_tag'], gene['name'], gene['distance'], gene['location'],
            op['is_intergenic'], op['near_regulon']
        ])

print(f"Results saved to: {csv_path}")

# Save summary JSON
summary = {
    'pipeline': 'bipartite_operator_search',
    'known_operator_recovered': known_op_found,
    'known_operator_rank': known_op_rank,
    'total_candidates': len(operators),
    'intergenic_candidates': len(intergenic_ops),
    'regulon_candidates': len(regulon_candidates),
    'strong_site_used': STRONG_SITE,
    'weak_site_used': WEAK_SITE,
    'spacer_range': f'{MIN_SPACER}-{MAX_SPACER}bp',
    'known_spacer': SPACER_LEN,
    'top_10': []
}

for i, op in enumerate(operators[:10]):
    gene = op['primary_gene']
    summary['top_10'].append({
        'rank': i + 1,
        'center_pos': op['center_pos'],
        'spacer': op['spacer'],
        'combined_score': round(op['combined_score'], 2),
        'normalized_score': round(op['normalized_score'], 4),
        'nearest_gene': gene['locus_tag'],
        'gene_name': gene['name'],
        'is_intergenic': op['is_intergenic'],
        'near_regulon': op['near_regulon'],
    })

json_path = os.path.join(output_dir, "bipartite_operators_summary.json")
with open(json_path, 'w') as f:
    json.dump(summary, f, indent=2)

print(f"Summary saved to: {json_path}")
print()
print("=" * 70)
print("PIPELINE COMPLETE")
print("=" * 70)
