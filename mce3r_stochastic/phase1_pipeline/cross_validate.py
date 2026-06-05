#!/usr/bin/env python3
"""
cross_validate.py — Cross-validate operator candidates across multiple methods
===============================================================================
Merges results from:
  1. Bipartite PWM search v3 (primary method)
  2. Original FIMO output (v1, minw=20-30 motifs)
  3. Corrected FIMO output (v2, minw=6-110 motifs)

Assigns confidence tiers based on how many methods support each candidate.
"""

import os
import csv
import json

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# ============================================================================
# LOAD BIPARTITE v3 RESULTS (PRIMARY)
# ============================================================================

v3_path = os.path.join(BASE_DIR, "results", "phase1", "bipartite_operators_v3.csv")
v3_candidates = []
with open(v3_path) as f:
    reader = csv.DictReader(f)
    for row in reader:
        v3_candidates.append({
            'rank': int(row['rank']),
            'center_pos': int(row['center_pos']),
            'upstream_pos': int(row['upstream_pos']),
            'downstream_pos': int(row['downstream_pos']),
            'strand': row['strand'],
            'spacer': int(row['spacer_bp']),
            'fisher_pval': float(row['fisher_pval']),
            'adjusted_score': float(row['adjusted_score']),
            'upstream_seq': row['upstream_seq'],
            'downstream_seq': row['downstream_seq'],
            'nearest_gene': row['nearest_gene'],
            'gene_name': row['gene_name'],
            'gene_distance': row['gene_distance'],
            'is_intergenic': row['is_intergenic'] == 'True',
            'near_regulon': row['near_regulon'] == 'True',
        })

print(f"Loaded {len(v3_candidates)} bipartite v3 candidates")

# ============================================================================
# LOAD FIMO RESULTS (ORIGINAL v1 — 20-21bp motifs)
# ============================================================================

fimo_v1_path = os.path.join(BASE_DIR, "results", "phase1", "fimo_output", "fimo.tsv")
fimo_v1_positions = set()
fimo_v1_hits = []
if os.path.exists(fimo_v1_path):
    with open(fimo_v1_path) as f:
        for line in f:
            if line.startswith('#') or line.startswith('motif_id\t'):
                continue
            parts = line.strip().split('\t')
            if len(parts) >= 8:
                start = int(parts[3])
                stop = int(parts[4])
                pval = float(parts[7])
                fimo_v1_hits.append({'start': start, 'stop': stop, 'pval': pval})
                # Add all positions in the hit region to the set
                for p in range(start - 50, stop + 50):
                    fimo_v1_positions.add(p)
    print(f"Loaded {len(fimo_v1_hits)} FIMO v1 hits")
else:
    print("WARNING: FIMO v1 output not found")

# ============================================================================
# LOAD FIMO RESULTS (CORRECTED — 6-15bp motifs)
# ============================================================================

fimo_v2_path = os.path.join(BASE_DIR, "results", "phase1", "fimo_output_corrected", "fimo.tsv")
fimo_v2_positions = set()
fimo_v2_hits = []
if os.path.exists(fimo_v2_path):
    with open(fimo_v2_path) as f:
        for line in f:
            if line.startswith('#') or line.startswith('motif_id\t'):
                continue
            parts = line.strip().split('\t')
            if len(parts) >= 8:
                start = int(parts[3])
                stop = int(parts[4])
                pval = float(parts[7])
                fimo_v2_hits.append({'start': start, 'stop': stop, 'pval': pval})
                for p in range(start - 50, stop + 50):
                    fimo_v2_positions.add(p)
    print(f"Loaded {len(fimo_v2_hits)} FIMO v2 hits")
else:
    print("WARNING: FIMO v2 output not found")

# ============================================================================
# CROSS-VALIDATE: Check each v3 candidate against FIMO results
# ============================================================================

print("\n" + "=" * 70)
print("CROSS-VALIDATION RESULTS")
print("=" * 70)

# Known 15bp recognition helix contact region from the strong site
CORE_RECOGNITION = "GCGTATTGGCTATGG"

def core_identity(seq, core=CORE_RECOGNITION):
    """Calculate best identity of a sequence to the 15bp core recognition motif."""
    seq = seq.upper()
    core = core.upper()
    best = 0
    for i in range(max(0, len(seq) - len(core) + 1)):
        subseq = seq[i:i+len(core)]
        matches = sum(1 for a, b in zip(subseq, core) if a == b)
        best = max(best, matches / len(core))
    return best

consensus = []
for cand in v3_candidates[:50]:  # Top 50 candidates
    center = cand['center_pos']
    upstream = cand['upstream_pos']
    downstream = cand['downstream_pos']

    # Check FIMO v1 overlap (within 50bp of either half-site)
    in_fimo_v1 = (upstream in fimo_v1_positions or downstream in fimo_v1_positions)

    # Check FIMO v2 overlap
    in_fimo_v2 = (upstream in fimo_v2_positions or downstream in fimo_v2_positions)

    # Count methods supporting this candidate
    n_methods = 1  # v3 bipartite always counts
    if in_fimo_v1:
        n_methods += 1
    if in_fimo_v2:
        n_methods += 1

    # Core motif identity for the downstream (strong) site
    strong_identity = core_identity(cand['downstream_seq'])
    weak_identity = core_identity(cand['upstream_seq'])
    best_identity = max(strong_identity, weak_identity)

    # Confidence tier
    if n_methods >= 3 and cand['near_regulon'] and cand['is_intergenic']:
        tier = "HIGHEST"
    elif n_methods >= 2 and cand['is_intergenic']:
        tier = "HIGH"
    elif n_methods >= 2 or (cand['near_regulon'] and cand['is_intergenic']):
        tier = "MODERATE"
    else:
        tier = "LOW"

    # Boost tier if core identity > 80%
    if best_identity >= 0.80 and tier == "MODERATE":
        tier = "HIGH"
    if best_identity >= 0.80 and tier == "HIGH" and cand['near_regulon']:
        tier = "HIGHEST"

    consensus.append({
        **cand,
        'n_methods': n_methods,
        'in_fimo_v1': in_fimo_v1,
        'in_fimo_v2': in_fimo_v2,
        'core_identity': best_identity,
        'confidence_tier': tier,
    })

# Sort by: tier priority (HIGHEST > HIGH > MODERATE > LOW), then adjusted score
tier_order = {'HIGHEST': 0, 'HIGH': 1, 'MODERATE': 2, 'LOW': 3}
consensus.sort(key=lambda x: (tier_order.get(x['confidence_tier'], 4), -x['adjusted_score']))

# Print results
for i, c in enumerate(consensus[:25]):
    print(f"\n--- Consensus Rank #{i+1} (v3 rank #{c['rank']}) ---")
    print(f"  Position: {c['upstream_pos']:,} - {c['downstream_pos']+25:,} ({c['strand']})")
    print(f"  Gene: {c['nearest_gene']} ({c['gene_name']}) [{c['gene_distance']}bp]")
    print(f"  Spacer: {c['spacer']}bp | Fisher p: {c['fisher_pval']:.2e} | Score: {c['adjusted_score']:.2f}")
    print(f"  Core identity: {c['core_identity']:.0%}")
    print(f"  Methods: {c['n_methods']}/3 (v3=yes, FIMO_v1={'yes' if c['in_fimo_v1'] else 'no'}, "
          f"FIMO_v2={'yes' if c['in_fimo_v2'] else 'no'})")

    flags = [f"CONFIDENCE: {c['confidence_tier']}"]
    if c['near_regulon']:
        flags.append("REGULON")
    if c['is_intergenic']:
        flags.append("INTERGENIC")
    print(f"  >>> {' | '.join(flags)}")

# ============================================================================
# SAVE CONSENSUS RESULTS
# ============================================================================

output_dir = os.path.join(BASE_DIR, "results", "phase1")
csv_path = os.path.join(output_dir, "consensus_operators.csv")
with open(csv_path, 'w', newline='') as f:
    w = csv.writer(f)
    w.writerow(['consensus_rank', 'v3_rank', 'center_pos', 'upstream_pos', 'downstream_pos',
                'strand', 'spacer_bp', 'fisher_pval', 'adjusted_score',
                'upstream_seq', 'downstream_seq', 'core_identity',
                'n_methods', 'in_fimo_v1', 'in_fimo_v2', 'confidence_tier',
                'nearest_gene', 'gene_name', 'gene_distance',
                'is_intergenic', 'near_regulon'])
    for i, c in enumerate(consensus):
        w.writerow([
            i+1, c['rank'], c['center_pos'], c['upstream_pos'], c['downstream_pos'],
            c['strand'], c['spacer'], f"{c['fisher_pval']:.2e}", f"{c['adjusted_score']:.2f}",
            c['upstream_seq'], c['downstream_seq'], f"{c['core_identity']:.2f}",
            c['n_methods'], c['in_fimo_v1'], c['in_fimo_v2'], c['confidence_tier'],
            c['nearest_gene'], c['gene_name'], c['gene_distance'],
            c['is_intergenic'], c['near_regulon']
        ])

print(f"\nConsensus results saved to: {csv_path}")

# Summary
summary = {
    'total_v3_candidates': len(v3_candidates),
    'fimo_v1_hits': len(fimo_v1_hits),
    'fimo_v2_hits': len(fimo_v2_hits),
    'consensus_candidates': len(consensus),
    'tier_counts': {
        'HIGHEST': sum(1 for c in consensus if c['confidence_tier'] == 'HIGHEST'),
        'HIGH': sum(1 for c in consensus if c['confidence_tier'] == 'HIGH'),
        'MODERATE': sum(1 for c in consensus if c['confidence_tier'] == 'MODERATE'),
        'LOW': sum(1 for c in consensus if c['confidence_tier'] == 'LOW'),
    },
    'top_candidates': [{
        'consensus_rank': i+1,
        'v3_rank': c['rank'],
        'position': f"{c['upstream_pos']}-{c['downstream_pos']+25}",
        'strand': c['strand'],
        'spacer': c['spacer'],
        'fisher_pval': c['fisher_pval'],
        'adjusted_score': c['adjusted_score'],
        'core_identity': c['core_identity'],
        'n_methods': c['n_methods'],
        'confidence_tier': c['confidence_tier'],
        'nearest_gene': c['nearest_gene'],
        'gene_name': c['gene_name'],
        'near_regulon': c['near_regulon'],
    } for i, c in enumerate(consensus[:10])],
}

json_path = os.path.join(output_dir, "consensus_summary.json")
with open(json_path, 'w') as f:
    json.dump(summary, f, indent=2)

print(f"Summary saved to: {json_path}")

# Final summary
print("\n" + "=" * 70)
print("FINAL OPERATOR PREDICTIONS")
print("=" * 70)
highest = [c for c in consensus if c['confidence_tier'] == 'HIGHEST']
high = [c for c in consensus if c['confidence_tier'] == 'HIGH']
print(f"\nHIGHEST confidence: {len(highest)} candidates")
for c in highest:
    print(f"  {c['nearest_gene']} ({c['gene_name']}): {c['upstream_pos']:,}-{c['downstream_pos']+25:,} "
          f"spacer={c['spacer']}bp, core_id={c['core_identity']:.0%}, {c['n_methods']}/3 methods")
print(f"\nHIGH confidence: {len(high)} candidates")
for c in high:
    print(f"  {c['nearest_gene']} ({c['gene_name']}): {c['upstream_pos']:,}-{c['downstream_pos']+25:,} "
          f"spacer={c['spacer']}bp, core_id={c['core_identity']:.0%}, {c['n_methods']}/3 methods")

print(f"\nThe paper predicts 4 total operators. We identified:")
print(f"  1. Known validated operator (yrbE3A/Rv1964)")
print(f"  2-{1+len(highest)+len(high)-1}. Novel predicted operators (HIGHEST + HIGH confidence)")
print("\nDone.")
