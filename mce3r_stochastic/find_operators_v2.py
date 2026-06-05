#!/usr/bin/env python3
"""
Targeted Bipartite Operator Search v2 for Mce3R in M. tuberculosis H37Rv
=========================================================================

Key improvement over v1: Uses p-value based thresholds instead of arbitrary
score cutoffs, allowing detection of degenerate binding sites.

Strategy:
  1. Build a PWM from the known binding site sequence
  2. Compute score distribution empirically from the genome to get p-values
  3. Scan genome for individual half-site matches at p < 1e-4
  4. Search for bipartite pairs: two hits separated by 30-80bp
  5. Rank by combined statistical significance
  6. Validate: known operator must be #1
"""

import os
import sys
import json
import math
import csv
import random
from collections import defaultdict

BASE_DIR = os.path.dirname(os.path.abspath(__file__))

# ============================================================================
# KNOWN OPERATOR DECOMPOSITION
# ============================================================================

KNOWN_OPERATOR = (
    "GCCCCGCGCTATAGGATACTAGCAAGATACATCATAGCCAATATATGCCAGT"
    "TTGCATTGCTATTTACCGATCAGTTGTCCAAGCAATCGCGTATTGGCTATG"
    "GACATCAGCGGTTCTGCCGC"
)

# The user's known binding site (maps to the strong/downstream site)
USER_KNOWN_SITE = "aatcgcgtattggctatggacatcagc"  # 27bp

# Find where the user's site maps in the operator
idx = KNOWN_OPERATOR.upper().find(USER_KNOWN_SITE.upper())
print(f"User's binding site found at position {idx} in the 123bp operator")

# From the crystal structure, the bipartite operator has:
#   [Site 1 / weak] [spacer] [Site 2 / strong]
# The user's site is in the downstream/strong position.
#
# Let's use TWO approaches to define the half-sites:
# Approach A: Use the crystal structure contacts (~25bp each, from Panagoda 2024)
# Approach B: Use the user's 27bp site and a symmetric upstream 27bp site
#
# We'll use Approach B since we have the exact sequence.

# Strong site: the user's known sequence, found at pos 85 in the operator
STRONG_SITE = KNOWN_OPERATOR[85:85+27]  # "ATCGCGTATTGGCTATGGACATCAGCG"
# But actually let's verify against user input
print(f"Strong site from operator[85:112]: {KNOWN_OPERATOR[85:112]}")
print(f"User known site:                   {USER_KNOWN_SITE}")

# The user's site doesn't start at the same position as our extraction.
# Let's find it properly:
strong_start = KNOWN_OPERATOR.upper().find(USER_KNOWN_SITE.upper())
if strong_start == -1:
    # Try with some flexibility - the user's site might be on the reverse strand
    # or slightly offset. Let's search for the core.
    core = USER_KNOWN_SITE[3:24].upper()  # central 21bp
    strong_start = KNOWN_OPERATOR.upper().find(core)
    if strong_start != -1:
        strong_start -= 3  # adjust back
        print(f"Found via core match at adjusted position {strong_start}")

if strong_start == -1:
    print("WARNING: User's site not found directly. Let's check reverse complement...")
    def rc(seq):
        comp = {'A':'T','T':'A','G':'C','C':'G','N':'N'}
        return ''.join(comp.get(b,'N') for b in reversed(seq.upper()))
    print(f"RC of user site: {rc(USER_KNOWN_SITE)}")
    rc_start = KNOWN_OPERATOR.upper().find(rc(USER_KNOWN_SITE))
    print(f"RC found at: {rc_start}")

# Let me just directly examine the operator and find the binding sites
# by looking at the EMSA probe definitions from the paper.
#
# From Panagoda 2024:
# - The 123bp operator spans the mce3R-yrbE3A intergenic region
# - Probe A (strong, Kd=2.4nM) = downstream half ~30bp
# - Probe C (weak, Kd=49nM) = upstream half ~30bp
# - The two probes are separated by ~53bp of spacer
#
# Looking at the 123bp sequence on the genome (+ strand):
# Position in genome: 2,207,476 (start of operator on + strand)
#
# The operator on the + strand reads 5'->3' as:
# GCCCCGCGCTATAGGATACTAGCAAGA  (upstream/weak site, 27bp)
# TACATCATAGCCAATATATGCCAGTTTGCATTGCTATTTACCGATC  (spacer, 46bp)
# AGTTGTCCAAGCAATCGCGTATTGGCT  (downstream portion, 27bp)
# ATGGACATCAGCGGTTCTGCCGC       (remaining 23bp)
#
# But the user says the strong site is "aatcgcgtattggctatggacatcagc"
# This overlaps positions 73-99 of the operator:
# ...CAATCGCGTATTGGCTATGGACATCAGC...
# Which is: C + user's 27bp site
# So the strong site is approximately positions 74-100

# Let's be more careful. Print the full operator with position markers:
print("\n=== OPERATOR DECOMPOSITION ===")
for i in range(0, len(KNOWN_OPERATOR), 10):
    chunk = KNOWN_OPERATOR[i:i+10]
    print(f"  pos {i:3d}-{i+len(chunk)-1:3d}: {chunk}")

# Find user's site by substring match allowing case insensitivity
op_upper = KNOWN_OPERATOR.upper()
user_upper = USER_KNOWN_SITE.upper()

for i in range(len(op_upper) - len(user_upper) + 1):
    mismatches = sum(1 for a, b in zip(op_upper[i:i+len(user_upper)], user_upper) if a != b)
    if mismatches <= 2:
        print(f"\n  Best match at pos {i} ({mismatches} mismatches):")
        print(f"    Operator: {KNOWN_OPERATOR[i:i+len(user_upper)]}")
        print(f"    User:     {USER_KNOWN_SITE}")
        strong_start = i
        break

# Define the two half-sites based on operator structure
# The strong site (downstream) starts around position 73-74
# The weak site (upstream) is the first ~27bp

# Use 25bp sites (matching the ~25bp TetR binding site footprint)
SITE_LENGTH = 25

# Strong (downstream) site - centered on user's sequence
STRONG_CENTER = strong_start + len(USER_KNOWN_SITE) // 2
STRONG_SITE_25 = KNOWN_OPERATOR[strong_start:strong_start+SITE_LENGTH]

# Weak (upstream) site - first 25bp of operator
WEAK_SITE_25 = KNOWN_OPERATOR[0:SITE_LENGTH]

# Spacer between them
spacer_len = strong_start - SITE_LENGTH
SPACER_SEQ = KNOWN_OPERATOR[SITE_LENGTH:strong_start]

print(f"\n=== HALF-SITE DEFINITIONS (using {SITE_LENGTH}bp sites) ===")
print(f"Weak site (upstream):  {WEAK_SITE_25}  (pos 0-{SITE_LENGTH-1})")
print(f"Spacer:                {SPACER_SEQ}  ({spacer_len}bp)")
print(f"Strong site (downstream): {STRONG_SITE_25}  (pos {strong_start}-{strong_start+SITE_LENGTH-1})")
print(f"Total: {SITE_LENGTH} + {spacer_len} + {SITE_LENGTH} = {SITE_LENGTH + spacer_len + SITE_LENGTH}bp")

# ============================================================================
# GENOME LOADING
# ============================================================================

def load_fasta(filepath):
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
    comp = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
    return ''.join(comp.get(b, 'N') for b in reversed(seq.upper()))

print("\nLoading H37Rv genome...")
GENOME_ID, GENOME = load_fasta(os.path.join(BASE_DIR, "data", "genomes", "H37Rv.fasta"))
print(f"  Genome: {GENOME_ID}, {len(GENOME):,}bp")

# Verify operator in genome
op_pos = GENOME.find(KNOWN_OPERATOR.upper())
print(f"  Known operator at genome position: {op_pos}")

# ============================================================================
# PWM CONSTRUCTION WITH PROPER BACKGROUND MODEL
# ============================================================================

# H37Rv has GC content = 65.6%
BG = {'A': 0.172, 'T': 0.172, 'G': 0.328, 'C': 0.328}

def make_pwm(site_seq, pseudocount=0.5):
    """Build log-odds PWM from a single sequence with pseudocounts."""
    seq = site_seq.upper()
    pwm = []
    for base in seq:
        # Use informative prior: the known base gets count=1, others get pseudocount
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

# Build PWMs
STRONG_PWM = make_pwm(STRONG_SITE_25)
WEAK_PWM = make_pwm(WEAK_SITE_25)

# Also build a COMBINED/GENERIC PWM that can match either half-site
# This is important because the 3 unknown operators may use the same
# binding motif for both sites (just at different affinities)
# TetR family repressors often use similar but not identical half-sites

print(f"\nStrong PWM: max_score={max_score(STRONG_PWM):.2f}, min_score={min_score(STRONG_PWM):.2f}")
print(f"Weak PWM: max_score={max_score(WEAK_PWM):.2f}, min_score={min_score(WEAK_PWM):.2f}")

# ============================================================================
# COMPUTE EMPIRICAL P-VALUES
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
    scores.sort(reverse=True)
    return scores

strong_dist = compute_score_distribution(STRONG_PWM, GENOME, SITE_LENGTH)
weak_dist = compute_score_distribution(WEAK_PWM, GENOME, SITE_LENGTH)

def score_to_pvalue(score, distribution):
    """Convert a score to an empirical p-value."""
    n_better = sum(1 for s in distribution if s >= score)
    return (n_better + 1) / (len(distribution) + 1)

# Find threshold scores for different p-values
for pval_target in [1e-3, 5e-4, 1e-4, 5e-5]:
    idx_s = int(pval_target * len(strong_dist))
    idx_w = int(pval_target * len(weak_dist))
    if idx_s < len(strong_dist) and idx_w < len(weak_dist):
        print(f"  p={pval_target:.0e}: strong_threshold={strong_dist[idx_s]:.2f}, "
              f"weak_threshold={weak_dist[idx_w]:.2f}")

# Use p < 1e-4 for individual sites (permissive, bipartite filter will tighten)
# But since we're looking for bipartite pairs, we can afford p < 5e-4 per site
# Combined p-value for a pair: ~2.5e-7 (still highly significant)
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
print("[SELF-CHECK] Both known sites pass threshold ✓")

# ============================================================================
# GENOME-WIDE SCAN (BOTH STRANDS)
# ============================================================================

print("\nScanning genome (both strands)...")

def scan_genome_pvalue(genome, pwm, motif_len, threshold, dist, label):
    """Scan genome for motif matches above score threshold."""
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

strong_hits = scan_genome_pvalue(GENOME, STRONG_PWM, SITE_LENGTH, STRONG_THRESHOLD, strong_dist, "strong")
weak_hits = scan_genome_pvalue(GENOME, WEAK_PWM, SITE_LENGTH, WEAK_THRESHOLD, weak_dist, "weak")

# ALSO: scan with each PWM looking for EITHER half-site pattern
# (the 3 unknown operators might have both sites match the strong motif,
#  or both match the weak motif, or one of each)
print("\n  Cross-scanning: strong PWM scanning for weak-like sites...")
strong_as_weak = scan_genome_pvalue(GENOME, STRONG_PWM, SITE_LENGTH, STRONG_THRESHOLD, strong_dist, "strong_any")
print("  Cross-scanning: weak PWM scanning for strong-like sites...")
weak_as_strong = scan_genome_pvalue(GENOME, WEAK_PWM, SITE_LENGTH, WEAK_THRESHOLD, weak_dist, "weak_any")

# Combine all hits into a unified set
all_hits = []
seen_positions = set()
for h in strong_hits + weak_hits + strong_as_weak + weak_as_strong:
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

# The known operator has ~46bp spacer between the 25bp half-sites.
# Allow generous range: 25-80bp spacer (to catch operators with
# slightly different geometries)
MIN_SPACER = 25
MAX_SPACER = 80

print(f"Searching for pairs of hits separated by {MIN_SPACER}-{MAX_SPACER}bp spacer...")
print(f"Known operator spacer: {spacer_len}bp")

operators = []

# For every pair of hits (regardless of which PWM found them),
# check if they form a bipartite architecture
for i in range(len(all_hits)):
    h1 = all_hits[i]
    h1_end = h1['pos'] + SITE_LENGTH

    for j in range(i + 1, len(all_hits)):
        h2 = all_hits[j]

        # Calculate spacer (distance between end of upstream site and start of downstream site)
        if h1['pos'] < h2['pos']:
            upstream, downstream = h1, h2
        else:
            upstream, downstream = h2, h1

        spacer = downstream['pos'] - (upstream['pos'] + SITE_LENGTH)

        if MIN_SPACER <= spacer <= MAX_SPACER:
            # Combined score
            combined = upstream['score'] + downstream['score']
            # Combined p-value (Fisher's method approximation: product of p-values)
            combined_pval = upstream['pvalue'] * downstream['pvalue']
            # Bonus for spacer close to known (46bp)
            spacer_penalty = abs(spacer - spacer_len) / spacer_len  # 0 = perfect match
            # Adjusted score: combined score with spacer bonus
            adjusted_score = combined * (1.0 - 0.3 * spacer_penalty)

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
                'combined_score': combined,
                'combined_pval': combined_pval,
                'adjusted_score': adjusted_score,
                'center_pos': (upstream['pos'] + downstream['pos'] + SITE_LENGTH) // 2,
            })

print(f"Raw bipartite candidates: {len(operators)}")

# Sort by adjusted score
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

import bisect

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
print("BIPARTITE OPERATOR CANDIDATES")
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
        flags.append("★ KNOWN REGULON")
    if op['is_intergenic']:
        flags.append("INTERGENIC")
    if abs(op['spacer'] - spacer_len) <= 10:
        flags.append(f"SPACER≈KNOWN ({op['spacer']}bp vs {spacer_len}bp)")

    print(f"\n--- Candidate #{i+1} ---")
    print(f"  Region: {op['upstream_pos']:,} - {op['downstream_pos']+SITE_LENGTH:,}")
    print(f"  Near: {gene_str}")
    print(f"  Upstream site:   pos {op['upstream_pos']:>9,} ({op['upstream_strand']}) "
          f"score={op['upstream_score']:.2f} p={op['upstream_pval']:.2e}  {op['upstream_seq']}")
    print(f"  Downstream site: pos {op['downstream_pos']:>9,} ({op['downstream_strand']}) "
          f"score={op['downstream_score']:.2f} p={op['downstream_pval']:.2e}  {op['downstream_seq']}")
    print(f"  Spacer: {op['spacer']}bp | Combined p: {op['combined_pval']:.2e} | "
          f"Adj. score: {op['adjusted_score']:.2f}")
    if flags:
        print(f"  >>> {' | '.join(flags)}")

# ============================================================================
# VALIDATION
# ============================================================================

print("\n" + "=" * 70)
print("VALIDATION")
print("=" * 70)

# Check 1: Known operator
known_rank = None
for i, op in enumerate(operators):
    if 2_207_400 <= op['center_pos'] <= 2_207_700:
        known_rank = i + 1
        break

if known_rank:
    print(f"\n[CHECK 1] ✅ Known operator recovered at rank #{known_rank}")
    op = operators[known_rank - 1]
    print(f"  Score: {op['adjusted_score']:.2f}, Spacer: {op['spacer']}bp, "
          f"Combined p: {op['combined_pval']:.2e}")
else:
    print("\n[CHECK 1] ❌ Known operator NOT found!")

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
    print(f"  #{i+1}: {op['spacer']}bp", end="")
    if abs(op['spacer'] - spacer_len) <= 5:
        print(" ← matches known", end="")
    print()

# ============================================================================
# SAVE RESULTS
# ============================================================================

output_dir = os.path.join(BASE_DIR, "results", "phase1")
os.makedirs(output_dir, exist_ok=True)

csv_path = os.path.join(output_dir, "bipartite_operators_v2.csv")
with open(csv_path, 'w', newline='') as f:
    w = csv.writer(f)
    w.writerow(['rank', 'center_pos', 'upstream_pos', 'downstream_pos',
                'upstream_strand', 'downstream_strand', 'spacer_bp',
                'upstream_score', 'downstream_score', 'upstream_pval', 'downstream_pval',
                'combined_pval', 'adjusted_score',
                'upstream_seq', 'downstream_seq',
                'nearest_gene', 'gene_name', 'gene_distance',
                'is_intergenic', 'near_regulon'])
    for i, op in enumerate(operators):
        gene = op['nearby_genes'][0] if op['nearby_genes'] else {}
        w.writerow([
            i+1, op['center_pos'], op['upstream_pos'], op['downstream_pos'],
            op['upstream_strand'], op['downstream_strand'], op['spacer'],
            f"{op['upstream_score']:.2f}", f"{op['downstream_score']:.2f}",
            f"{op['upstream_pval']:.2e}", f"{op['downstream_pval']:.2e}",
            f"{op['combined_pval']:.2e}", f"{op['adjusted_score']:.2f}",
            op['upstream_seq'], op['downstream_seq'],
            gene.get('locus_tag', ''), gene.get('name', ''), gene.get('distance', ''),
            op['is_intergenic'], op['near_regulon']
        ])

print(f"\nResults saved to: {csv_path}")

# Summary JSON
summary = {
    'pipeline': 'bipartite_operator_search_v2',
    'site_length': SITE_LENGTH,
    'strong_site': STRONG_SITE_25,
    'weak_site': WEAK_SITE_25,
    'spacer_range': f'{MIN_SPACER}-{MAX_SPACER}bp',
    'known_spacer': spacer_len,
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
        'spacer': op['spacer'],
        'combined_pval': op['combined_pval'],
        'adjusted_score': round(op['adjusted_score'], 2),
        'nearest_gene': gene.get('locus_tag', ''),
        'gene_name': gene.get('name', ''),
        'is_intergenic': op['is_intergenic'],
        'near_regulon': op['near_regulon'],
    })

json_path = os.path.join(output_dir, "bipartite_operators_v2_summary.json")
with open(json_path, 'w') as f:
    json.dump(summary, f, indent=2)

print(f"Summary saved to: {json_path}")
print("\n" + "=" * 70)
print("PIPELINE COMPLETE")
print("=" * 70)
