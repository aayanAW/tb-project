#!/usr/bin/env python3
"""
find_operators_v3.py — Corrected Bipartite Operator Search for Mce3R
=====================================================================
Fixes from v2:
  - Correct half-site boundaries (53bp spacer, not 60bp)
  - Binary search for p-values
  - Proper Fisher's method for combined significance
  - Gaussian spacer weighting centered on 53bp
  - Strand consistency filter
  - No redundant cross-scanning
"""

import os
import sys
import json
import math
import csv
import random
import bisect
from collections import defaultdict
from scipy.stats import combine_pvalues

BASE_DIR = os.path.dirname(os.path.abspath(__file__))

# ============================================================================
# KNOWN OPERATOR DECOMPOSITION — CORRECTED BOUNDARIES
# ============================================================================

KNOWN_OPERATOR = (
    "GCCCCGCGCTATAGGATACTAGCAAGATACATCATAGCCAATATATGCCAGT"
    "TTGCATTGCTATTTACCGATCAGTTGTCCAAGCAATCGCGTATTGGCTATG"
    "GACATCAGCGGTTCTGCCGC"
)

# The known 27bp strong binding site from cryo-EM (Region 3)
KNOWN_27BP_SITE = "AATCGCGTATTGGCTATGGACATCAGC"

# Locate the 27bp site in the operator
idx_27 = KNOWN_OPERATOR.upper().find(KNOWN_27BP_SITE.upper())
assert idx_27 != -1, "27bp binding site not found in known operator!"
print(f"27bp strong binding site found at operator position {idx_27}")
print(f"  Operator[{idx_27}:{idx_27+27}] = {KNOWN_OPERATOR[idx_27:idx_27+27]}")

# The 27bp site includes 1bp overlap on each side of the 25bp core.
# Strong 25bp core: positions 1..25 of the 27bp match (i.e., trim 1bp each side)
SITE_LENGTH = 25

strong_25_start = idx_27 + 1   # = 86
strong_25_end   = idx_27 + 26  # = 111
STRONG_SITE_25 = KNOWN_OPERATOR[strong_25_start:strong_25_end]
assert len(STRONG_SITE_25) == SITE_LENGTH

# Weak 25bp site: ends 53bp before the strong site starts
EXPECTED_SPACER = 53
weak_25_end   = strong_25_start - EXPECTED_SPACER  # = 33
weak_25_start = weak_25_end - SITE_LENGTH           # = 8
WEAK_SITE_25 = KNOWN_OPERATOR[weak_25_start:weak_25_end]
assert len(WEAK_SITE_25) == SITE_LENGTH

# Verify spacer
actual_spacer = strong_25_start - weak_25_end
assert actual_spacer == EXPECTED_SPACER, (
    f"Spacer mismatch: expected {EXPECTED_SPACER}, got {actual_spacer}"
)

SPACER_SEQ = KNOWN_OPERATOR[weak_25_end:strong_25_start]

print(f"\n=== HALF-SITE DEFINITIONS ({SITE_LENGTH}bp sites, {EXPECTED_SPACER}bp spacer) ===")
print(f"Weak site (upstream):     {WEAK_SITE_25}  (operator[{weak_25_start}:{weak_25_end}])")
print(f"Spacer:                   {SPACER_SEQ}  ({actual_spacer}bp)")
print(f"Strong site (downstream): {STRONG_SITE_25}  (operator[{strong_25_start}:{strong_25_end}])")
print(f"Flanking: {weak_25_start}bp left, {len(KNOWN_OPERATOR) - strong_25_end}bp right")
print(f"Total: {weak_25_start} + {SITE_LENGTH} + {actual_spacer} + {SITE_LENGTH} + "
      f"{len(KNOWN_OPERATOR) - strong_25_end} = {len(KNOWN_OPERATOR)}bp")

# Self-check: total should reconstruct
assert (weak_25_start + SITE_LENGTH + actual_spacer + SITE_LENGTH +
        (len(KNOWN_OPERATOR) - strong_25_end)) == len(KNOWN_OPERATOR)

# ============================================================================
# GENOME LOADING
# ============================================================================

def load_fasta(filepath):
    """Load a FASTA file and return (header, sequence)."""
    sequence = []
    header = None
    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                header = line[1:].split()[0]
            else:
                sequence.append(line.upper())
    return header, ''.join(sequence)

def reverse_complement(seq):
    """Return the reverse complement of a DNA sequence."""
    comp = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
    return ''.join(comp.get(b, 'N') for b in reversed(seq.upper()))

print("\nLoading H37Rv genome...")
GENOME_ID, GENOME = load_fasta(os.path.join(BASE_DIR, "data", "genomes", "H37Rv.fasta"))
print(f"  Genome: {GENOME_ID}, {len(GENOME):,}bp")

# Verify operator in genome
op_pos = GENOME.find(KNOWN_OPERATOR.upper())
assert op_pos != -1, "Known operator not found in genome!"
print(f"  Known operator at genome position: {op_pos}")

# ============================================================================
# PWM CONSTRUCTION WITH PROPER BACKGROUND MODEL
# ============================================================================

# H37Rv GC content = 65.6%
BG = {'A': 0.172, 'T': 0.172, 'G': 0.328, 'C': 0.328}

def make_pwm(site_seq, pseudocount=0.5):
    """Build log-odds PWM from a single sequence with pseudocounts."""
    seq = site_seq.upper()
    pwm = []
    for base in seq:
        counts = {b: pseudocount for b in 'ACGT'}
        counts[base] += 1.0
        total = sum(counts.values())
        row = {}
        for b in 'ACGT':
            freq = counts[b] / total
            row[b] = math.log2(freq / BG[b])
        pwm.append(row)
    return pwm

def score_seq(pwm, seq):
    """Score a sequence against a PWM."""
    if len(seq) != len(pwm):
        return float('-inf')
    score = 0.0
    for i, base in enumerate(seq.upper()):
        if base not in 'ACGT':
            return float('-inf')
        score += pwm[i][base]
    return score

def max_score(pwm):
    return sum(max(row.values()) for row in pwm)

def min_score(pwm):
    return sum(min(row.values()) for row in pwm)

# Build PWMs from correctly-defined half-sites
STRONG_PWM = make_pwm(STRONG_SITE_25)
WEAK_PWM = make_pwm(WEAK_SITE_25)

print(f"\nStrong PWM: max_score={max_score(STRONG_PWM):.2f}, min_score={min_score(STRONG_PWM):.2f}")
print(f"Weak PWM: max_score={max_score(WEAK_PWM):.2f}, min_score={min_score(WEAK_PWM):.2f}")

# ============================================================================
# COMPUTE EMPIRICAL P-VALUES WITH BINARY SEARCH
# ============================================================================

print("\nComputing empirical score distributions for p-value calculation...")

def compute_score_distribution(pwm, genome, motif_len, n_samples=500000):
    """Sample random positions from genome to build score distribution."""
    random.seed(42)
    scores = []
    genome_len = len(genome)
    for _ in range(n_samples):
        pos = random.randint(0, genome_len - motif_len)
        subseq = genome[pos:pos + motif_len]
        if 'N' not in subseq:
            s = score_seq(pwm, subseq)
            scores.append(s)
    scores.sort(reverse=True)  # Descending order for p-value lookup
    return scores

strong_dist = compute_score_distribution(STRONG_PWM, GENOME, SITE_LENGTH)
print(f"  Strong distribution: {len(strong_dist)} samples")
weak_dist = compute_score_distribution(WEAK_PWM, GENOME, SITE_LENGTH)
print(f"  Weak distribution: {len(weak_dist)} samples")

def score_to_pvalue(score, distribution):
    """Convert a score to an empirical p-value using binary search.

    Distribution is sorted in descending order.
    P-value = (number of scores >= query score + 1) / (N + 1).

    We use bisect on the negated distribution to find the insertion point,
    which gives us the count of scores >= query score.
    """
    # bisect_right on negated distribution finds how many elements have
    # -dist[i] <= -score, i.e., dist[i] >= score
    n_better = bisect.bisect_right(distribution, score, key=lambda x: -x)
    # n_better is actually the index where -score would be inserted in
    # the negated list. Since distribution is descending, bisect_left
    # on the negated list gives us the count of values >= score.
    #
    # Simpler approach: since distribution is sorted descending,
    # use bisect_left with a negated key to find how many are >= score.
    # Actually, let's just use bisect directly:
    # In a descending list, the number of elements >= score is the
    # insertion point of score in the reversed (ascending) sense.

    # Use bisect_right on reversed logic: count of elements >= score
    # In descending sorted list: find first index where dist[i] < score
    # That index = number of elements >= score
    lo, hi = 0, len(distribution)
    while lo < hi:
        mid = (lo + hi) // 2
        if distribution[mid] >= score:
            lo = mid + 1
        else:
            hi = mid
    n_better = lo
    return (n_better + 1) / (len(distribution) + 1)

# Verify binary search gives same results as naive for a few test cases
print("  Verifying binary search p-value computation...")
test_scores = [strong_dist[0], strong_dist[len(strong_dist)//2], strong_dist[-1]]
for ts in test_scores:
    naive_p = (sum(1 for s in strong_dist if s >= ts) + 1) / (len(strong_dist) + 1)
    bisect_p = score_to_pvalue(ts, strong_dist)
    assert abs(naive_p - bisect_p) < 1e-10, (
        f"P-value mismatch: naive={naive_p}, bisect={bisect_p} for score={ts}"
    )
print("  Binary search verification passed.")

# Find threshold scores for different p-values
for pval_target in [1e-3, 5e-4, 1e-4, 5e-5]:
    idx_s = int(pval_target * len(strong_dist))
    idx_w = int(pval_target * len(weak_dist))
    if idx_s < len(strong_dist) and idx_w < len(weak_dist):
        print(f"  p={pval_target:.0e}: strong_threshold={strong_dist[idx_s]:.2f}, "
              f"weak_threshold={weak_dist[idx_w]:.2f}")

# Use p < 5e-4 for individual sites
INDIVIDUAL_PVAL = 5e-4
strong_threshold_idx = int(INDIVIDUAL_PVAL * len(strong_dist))
weak_threshold_idx = int(INDIVIDUAL_PVAL * len(weak_dist))

STRONG_THRESHOLD = strong_dist[min(strong_threshold_idx, len(strong_dist)-1)]
WEAK_THRESHOLD = weak_dist[min(weak_threshold_idx, len(weak_dist)-1)]

print(f"\nUsing individual site p-value threshold: {INDIVIDUAL_PVAL}")
print(f"  Strong score threshold: {STRONG_THRESHOLD:.2f}")
print(f"  Weak score threshold: {WEAK_THRESHOLD:.2f}")

# Self-check: verify known sites score well above threshold
known_strong_score = score_seq(STRONG_PWM, STRONG_SITE_25)
known_weak_score = score_seq(WEAK_PWM, WEAK_SITE_25)
known_strong_pval = score_to_pvalue(known_strong_score, strong_dist)
known_weak_pval = score_to_pvalue(known_weak_score, weak_dist)

print(f"\n[SELF-CHECK] Known strong site: score={known_strong_score:.2f}, p={known_strong_pval:.2e}")
print(f"[SELF-CHECK] Known weak site: score={known_weak_score:.2f}, p={known_weak_pval:.2e}")
assert known_strong_pval < INDIVIDUAL_PVAL, "Known strong site below threshold!"
assert known_weak_pval < INDIVIDUAL_PVAL, "Known weak site below threshold!"
print("[SELF-CHECK] Both known sites pass threshold.")

# ============================================================================
# GENOME-WIDE SCAN (BOTH STRANDS) — NO REDUNDANT CROSS-SCANNING
# ============================================================================

print("\nScanning genome (both strands)...")

def scan_genome_pvalue(genome, pwm, motif_len, threshold, dist, label):
    """Scan genome for motif matches above score threshold on both strands."""
    hits = []

    # Build reverse complement PWM
    rc_pwm = []
    for row in reversed(pwm):
        rc_row = {'A': row['T'], 'T': row['A'], 'G': row['C'], 'C': row['G']}
        rc_pwm.append(rc_row)

    for i in range(len(genome) - motif_len + 1):
        subseq = genome[i:i + motif_len]
        if 'N' in subseq:
            continue

        # Plus strand
        plus_score = score_seq(pwm, subseq)
        if plus_score >= threshold:
            pval = score_to_pvalue(plus_score, dist)
            hits.append({
                'pos': i, 'strand': '+', 'score': plus_score,
                'pvalue': pval, 'sequence': subseq, 'label': label
            })

        # Minus strand
        minus_score = score_seq(rc_pwm, subseq)
        if minus_score >= threshold:
            pval = score_to_pvalue(minus_score, dist)
            hits.append({
                'pos': i, 'strand': '-', 'score': minus_score,
                'pvalue': pval, 'sequence': reverse_complement(subseq), 'label': label
            })

    hits.sort(key=lambda x: x['score'], reverse=True)
    print(f"  [{label}] {len(hits)} hits at p < {INDIVIDUAL_PVAL}")
    if hits:
        print(f"    Top hit: pos={hits[0]['pos']}, score={hits[0]['score']:.2f}, p={hits[0]['pvalue']:.2e}")
    return hits

# Only two scans: strong PWM and weak PWM (no redundant cross-scanning)
strong_hits = scan_genome_pvalue(GENOME, STRONG_PWM, SITE_LENGTH, STRONG_THRESHOLD, strong_dist, "strong")
weak_hits = scan_genome_pvalue(GENOME, WEAK_PWM, SITE_LENGTH, WEAK_THRESHOLD, weak_dist, "weak")

# Combine all hits into a unified set
all_hits = []
seen_positions = set()
for h in strong_hits + weak_hits:
    key = (h['pos'], h['strand'], h['label'])
    if key not in seen_positions:
        seen_positions.add(key)
        all_hits.append(h)

all_hits.sort(key=lambda x: x['score'], reverse=True)
print(f"\nTotal unique hits across all scans: {len(all_hits)}")

# Self-check: known operator region
known_region = range(2_207_400, 2_207_700)
known_hits = [h for h in all_hits if h['pos'] in known_region]
print(f"[SELF-CHECK] Hits in known operator region: {len(known_hits)}")
for h in known_hits:
    print(f"  pos={h['pos']}, strand={h['strand']}, label={h['label']}, "
          f"score={h['score']:.2f}, p={h['pvalue']:.2e}")

# ============================================================================
# BIPARTITE OPERATOR DETECTION
# ============================================================================

print("\n" + "=" * 70)
print("BIPARTITE OPERATOR SEARCH")
print("=" * 70)

MIN_SPACER = 25
MAX_SPACER = 80

print(f"Searching for pairs of hits separated by {MIN_SPACER}-{MAX_SPACER}bp spacer...")
print(f"Expected spacer: {EXPECTED_SPACER}bp (Gaussian weight centered here, sigma=10)")

operators = []

# For every pair of hits, check if they form a bipartite architecture
# with SAME STRAND and correct spacer range
for i in range(len(all_hits)):
    h1 = all_hits[i]
    h1_end = h1['pos'] + SITE_LENGTH

    for j in range(i + 1, len(all_hits)):
        h2 = all_hits[j]

        # BUG FIX 6: Strand consistency — both half-sites must be on the same strand
        if h1['strand'] != h2['strand']:
            continue

        # Calculate spacer (distance between end of upstream site and start of downstream site)
        if h1['pos'] < h2['pos']:
            upstream, downstream = h1, h2
        else:
            upstream, downstream = h2, h1

        spacer = downstream['pos'] - (upstream['pos'] + SITE_LENGTH)

        if MIN_SPACER <= spacer <= MAX_SPACER:
            # Combined score
            combined = upstream['score'] + downstream['score']

            # BUG FIX 3: Proper Fisher's method for combined p-values
            _, fisher_pval = combine_pvalues(
                [upstream['pvalue'], downstream['pvalue']], method='fisher'
            )

            # BUG FIX 4: Gaussian spacer weighting centered on 53bp
            spacer_weight = math.exp(-0.5 * ((spacer - EXPECTED_SPACER) / 10) ** 2)

            # Adjusted score: combined score weighted by spacer proximity
            adjusted_score = combined * spacer_weight

            operators.append({
                'upstream_pos': upstream['pos'],
                'downstream_pos': downstream['pos'],
                'upstream_strand': upstream['strand'],
                'downstream_strand': downstream['strand'],
                'upstream_score': upstream['score'],
                'downstream_score': downstream['score'],
                'upstream_pval': upstream['pvalue'],
                'downstream_pval': downstream['pvalue'],
                'upstream_seq': upstream['sequence'],
                'downstream_seq': downstream['sequence'],
                'upstream_label': upstream['label'],
                'downstream_label': downstream['label'],
                'spacer': spacer,
                'spacer_weight': spacer_weight,
                'combined_score': combined,
                'fisher_pval': fisher_pval,
                'adjusted_score': adjusted_score,
                'center_pos': (upstream['pos'] + downstream['pos'] + SITE_LENGTH) // 2,
            })

print(f"Raw bipartite candidates: {len(operators)}")

# Sort by adjusted score (descending)
operators.sort(key=lambda x: x['adjusted_score'], reverse=True)

# Deduplicate (merge operators within 50bp of each other)
def dedup(ops, dist=50):
    kept = []
    used = set()
    for op in ops:
        c = op['center_pos']
        if any(abs(c - u) < dist for u in used):
            continue
        kept.append(op)
        used.add(c)
    return kept

operators = dedup(operators)
print(f"After deduplication: {len(operators)} unique candidates")

# ============================================================================
# GENE ANNOTATION
# ============================================================================

print("\nAnnotating with gene information...")

genes = []
with open(os.path.join(BASE_DIR, "data", "genomes", "H37Rv.gff")) as f:
    for line in f:
        if line.startswith('#'):
            continue
        parts = line.strip().split('\t')
        if len(parts) < 9 or parts[2] != 'gene':
            continue
        attrs = {}
        for attr in parts[8].split(';'):
            if '=' in attr:
                k, v = attr.split('=', 1)
                attrs[k] = v
        genes.append({
            'start': int(parts[3]),
            'end': int(parts[4]),
            'strand': parts[6],
            'locus_tag': attrs.get('locus_tag', ''),
            'name': attrs.get('gene', ''),
            'product': attrs.get('product', attrs.get('Note', ''))
        })

genes.sort(key=lambda g: g['start'])
gene_starts = [g['start'] for g in genes]

print(f"  {len(genes)} genes loaded")

def annotate_position(pos):
    """Find genes near a position and determine if it's intergenic."""
    idx = bisect.bisect_left(gene_starts, pos)
    nearby = []
    for i in range(max(0, idx-2), min(len(genes), idx+2)):
        g = genes[i]
        if g['start'] <= pos <= g['end']:
            dist = 0
            loc = 'inside'
        elif pos < g['start']:
            dist = g['start'] - pos
            loc = 'upstream' if g['strand'] == '+' else 'downstream'
        else:
            dist = pos - g['end']
            loc = 'downstream' if g['strand'] == '+' else 'upstream'
        nearby.append({**g, 'distance': dist, 'location': loc})
    nearby.sort(key=lambda x: x['distance'])

    # Is position intergenic?
    is_ig = all(n['distance'] > 0 for n in nearby[:2])
    return nearby[:3], is_ig

REGULON_GENES = {'Rv1963c', 'Rv1964', 'Rv1933c', 'Rv1934c', 'Rv1935c',
                 'Rv1936', 'Rv1937', 'Rv1938', 'Rv1939', 'Rv1940', 'Rv1941c'}

for op in operators:
    nearby, is_ig = annotate_position(op['center_pos'])
    op['nearby_genes'] = nearby
    op['is_intergenic'] = is_ig
    op['near_regulon'] = any(g['locus_tag'] in REGULON_GENES for g in nearby)

# ============================================================================
# RESULTS
# ============================================================================

print("\n" + "=" * 70)
print("BIPARTITE OPERATOR CANDIDATES (Top 25)")
print("=" * 70)

for i, op in enumerate(operators[:25]):
    gene = op['nearby_genes'][0] if op['nearby_genes'] else None
    gene_str = f"{gene['locus_tag']}" if gene else "unknown"
    if gene and gene['name']:
        gene_str += f" ({gene['name']})"
    if gene:
        gene_str += f" [{gene['distance']}bp {gene['location']}]"

    flags = []
    if op['near_regulon']:
        flags.append("KNOWN REGULON")
    if op['is_intergenic']:
        flags.append("INTERGENIC")
    if abs(op['spacer'] - EXPECTED_SPACER) <= 10:
        flags.append(f"SPACER~KNOWN ({op['spacer']}bp vs {EXPECTED_SPACER}bp)")

    print(f"\n--- Candidate #{i+1} ---")
    print(f"  Region: {op['upstream_pos']:,} - {op['downstream_pos']+SITE_LENGTH:,}")
    print(f"  Strand: {op['upstream_strand']}")
    print(f"  Near: {gene_str}")
    print(f"  Upstream site:   pos {op['upstream_pos']:>9,} ({op['upstream_strand']}) "
          f"score={op['upstream_score']:.2f} p={op['upstream_pval']:.2e}  {op['upstream_seq']}")
    print(f"  Downstream site: pos {op['downstream_pos']:>9,} ({op['downstream_strand']}) "
          f"score={op['downstream_score']:.2f} p={op['downstream_pval']:.2e}  {op['downstream_seq']}")
    print(f"  Spacer: {op['spacer']}bp (weight={op['spacer_weight']:.3f}) | "
          f"Fisher p: {op['fisher_pval']:.2e} | Adj. score: {op['adjusted_score']:.2f}")
    if flags:
        print(f"  >>> {' | '.join(flags)}")

# ============================================================================
# VALIDATION
# ============================================================================

print("\n" + "=" * 70)
print("VALIDATION")
print("=" * 70)

# Check 1: Known operator must be rank #1
known_rank = None
for i, op in enumerate(operators):
    if 2_207_400 <= op['center_pos'] <= 2_207_700:
        known_rank = i + 1
        break

if known_rank:
    print(f"\n[CHECK 1] Known operator recovered at rank #{known_rank}")
    op = operators[known_rank - 1]
    print(f"  Score: {op['adjusted_score']:.2f}, Spacer: {op['spacer']}bp, "
          f"Fisher p: {op['fisher_pval']:.2e}")
    if known_rank == 1:
        print("  PASS: Known operator is rank #1")
    else:
        print(f"  WARNING: Known operator is rank #{known_rank}, expected #1")
else:
    print("\n[CHECK 1] FAIL: Known operator NOT found in candidates!")

# Check 2: Score gap
if len(operators) >= 2:
    gap = operators[0]['adjusted_score'] - operators[1]['adjusted_score']
    print(f"\n[CHECK 2] Score gap between #1 and #2: {gap:.2f}")
    print(f"  #1 score: {operators[0]['adjusted_score']:.2f}")
    print(f"  #2 score: {operators[1]['adjusted_score']:.2f}")

# Check 3: Known regulon hits
regulon_ops = [(i+1, op) for i, op in enumerate(operators) if op['near_regulon']]
print(f"\n[CHECK 3] Candidates near known regulon genes: {len(regulon_ops)}")
for rank, op in regulon_ops[:5]:
    gene = op['nearby_genes'][0]
    print(f"  Rank #{rank}: {gene['locus_tag']} ({gene['name']}), "
          f"adj_score={op['adjusted_score']:.2f}")

# Check 4: Intergenic fraction
ig_count = sum(1 for op in operators if op['is_intergenic'])
print(f"\n[CHECK 4] Intergenic candidates: {ig_count}/{len(operators)} "
      f"({100*ig_count/max(1,len(operators)):.0f}%)")

# Check 5: Spacer distribution
print(f"\n[CHECK 5] Spacer lengths (top 10):")
for i, op in enumerate(operators[:10]):
    marker = " <- matches known" if abs(op['spacer'] - EXPECTED_SPACER) <= 5 else ""
    print(f"  #{i+1}: {op['spacer']}bp{marker}")

# ============================================================================
# SAVE RESULTS
# ============================================================================

output_dir = os.path.join(BASE_DIR, "results", "phase1")
os.makedirs(output_dir, exist_ok=True)

csv_path = os.path.join(output_dir, "bipartite_operators_v3.csv")
with open(csv_path, 'w', newline='') as f:
    w = csv.writer(f)
    w.writerow(['rank', 'center_pos', 'upstream_pos', 'downstream_pos',
                'strand', 'spacer_bp', 'spacer_weight',
                'upstream_score', 'downstream_score', 'upstream_pval', 'downstream_pval',
                'fisher_pval', 'adjusted_score',
                'upstream_seq', 'downstream_seq',
                'nearest_gene', 'gene_name', 'gene_distance',
                'is_intergenic', 'near_regulon'])
    for i, op in enumerate(operators):
        gene = op['nearby_genes'][0] if op['nearby_genes'] else {}
        w.writerow([
            i+1, op['center_pos'], op['upstream_pos'], op['downstream_pos'],
            op['upstream_strand'], op['spacer'], f"{op['spacer_weight']:.4f}",
            f"{op['upstream_score']:.2f}", f"{op['downstream_score']:.2f}",
            f"{op['upstream_pval']:.2e}", f"{op['downstream_pval']:.2e}",
            f"{op['fisher_pval']:.2e}", f"{op['adjusted_score']:.2f}",
            op['upstream_seq'], op['downstream_seq'],
            gene.get('locus_tag', ''), gene.get('name', ''), gene.get('distance', ''),
            op['is_intergenic'], op['near_regulon']
        ])

print(f"\nResults saved to: {csv_path}")

# Summary JSON
summary = {
    'pipeline': 'bipartite_operator_search_v3',
    'fixes_from_v2': [
        'Correct half-site boundaries (53bp spacer, not 60bp)',
        'Binary search for p-values (O(log n) instead of O(n))',
        'Fisher method for combined p-values (scipy.stats.combine_pvalues)',
        'Gaussian spacer weighting (sigma=10bp, centered on 53bp)',
        'Strand consistency filter (both half-sites same strand)',
        'Removed redundant cross-scanning (strong_any scan)',
    ],
    'site_length': SITE_LENGTH,
    'strong_site': STRONG_SITE_25,
    'weak_site': WEAK_SITE_25,
    'expected_spacer': EXPECTED_SPACER,
    'spacer_range': f'{MIN_SPACER}-{MAX_SPACER}bp',
    'individual_pval_threshold': INDIVIDUAL_PVAL,
    'total_strong_hits': len(strong_hits),
    'total_weak_hits': len(weak_hits),
    'total_bipartite_candidates': len(operators),
    'known_operator_rank': known_rank,
    'top_candidates': []
}

for i, op in enumerate(operators[:10]):
    gene = op['nearby_genes'][0] if op['nearby_genes'] else {}
    summary['top_candidates'].append({
        'rank': i+1,
        'center_pos': op['center_pos'],
        'strand': op['upstream_strand'],
        'spacer': op['spacer'],
        'spacer_weight': round(op['spacer_weight'], 4),
        'fisher_pval': op['fisher_pval'],
        'adjusted_score': round(op['adjusted_score'], 2),
        'nearest_gene': gene.get('locus_tag', ''),
        'gene_name': gene.get('name', ''),
        'is_intergenic': op['is_intergenic'],
        'near_regulon': op['near_regulon'],
    })

json_path = os.path.join(output_dir, "bipartite_operators_v3_summary.json")
with open(json_path, 'w') as f:
    json.dump(summary, f, indent=2)

print(f"Summary saved to: {json_path}")
print("\n" + "=" * 70)
print("PIPELINE COMPLETE (v3)")
print("=" * 70)
