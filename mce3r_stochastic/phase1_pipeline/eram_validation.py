"""
Eram validation: Extract the 5 intergenic regions from Eram's analysis,
run extended MEME searches, check for Eram's and Santangelo's motifs,
and create a DNA map of all motifs relative to the Panagoda 2024 operator.

5 Regions (from Eram Kabir's research summary):
  Region 1: mce3R (Rv1963c) - yrbE3A (Rv1964)     ~225 bp (898 bp in Eram's fig)
  Region 2: echA13 (Rv1935c) - Rv1936              ~210 bp
  Region 3: Rv1941 - Rv1942c                       (not Mce3R/VapBC35 candidate)
  Region 4: Mce3F (Rv1970) - Rv1972                (not Mce3R/VapBC35 candidate)
  Region 5: Rv1944c - Rv1945                        ~55 bp
"""

import sys
import os
import re
import json
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from Bio import SeqIO
from Bio.Seq import Seq
from config.parameters import PARAMS

PROJECT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS_DIR = os.path.join(PROJECT, 'results', 'eram_validation')
GENOME_PATH = os.path.join(PROJECT, 'data', 'genomes', 'H37Rv.gb')

# Eram's motifs (from Figures 6-8)
ERAM_MOTIF1_CONSENSUS = "ACTAGCAAGATACATCATAGCCAATATATGCCAGTTTGCATTGCTATTT"  # 49bp, Fig 6
ERAM_MOTIF2_SITES = [  # 21bp, Fig 7 — individual site sequences
    "ATAGCATACTACAACATACA",
    "ATCGCGTATTGGCTATGGACA",
    "ATATCGGACTAACAAAATACA",
    "AGTGCATATCAGTAATAGACA",
    "ATCGGTAAATAACAATGCAAA",
]
ERAM_MOTIF3_CONSENSUS = "TACAACGAAG"  # 10bp, Fig 8

# Santangelo 2009 motifs (from Eram's appendix, page 17)
SANTANGELO_FORWARD = [
    "gactaacaaaatacat",    # Forward strand, Region 1
    "tatcagtaatagacat",
    "tactagcaagatacat",
    "tattggctatggacat",
]
SANTANGELO_REVERSE = [
    "taaatagcaatgcaaactgg",  # Reverse strand (Santangelo)
    "taaattgcaatgtaatcgcg",
]
SANTANGELO_RC_DIRECT = [
    "gattacattgcaattta",     # RC of Santangelo, direct in Region 1
    "cagtttgcattgctattta",
]


def load_genome():
    """Load H37Rv genome sequence."""
    record = SeqIO.read(GENOME_PATH, 'genbank')
    return str(record.seq).upper(), record


def get_gene_coordinates(record, locus_tag):
    """Find CDS coordinates for a given locus tag."""
    for feature in record.features:
        if feature.type == 'CDS':
            lt = ' '.join(feature.qualifiers.get('locus_tag', []))
            if lt == locus_tag:
                return int(feature.location.start), int(feature.location.end), feature.location.strand
    return None, None, None


def extract_intergenic(genome, record, gene1_locus, gene2_locus, label):
    """Extract intergenic region between two genes."""
    s1, e1, strand1 = get_gene_coordinates(record, gene1_locus)
    s2, e2, strand2 = get_gene_coordinates(record, gene2_locus)

    if s1 is None or s2 is None:
        print(f"  WARNING: Could not find {gene1_locus} or {gene2_locus}")
        return None

    # Intergenic region is between the two genes
    if e1 < s2:
        igr_start = e1
        igr_end = s2
    elif e2 < s1:
        igr_start = e2
        igr_end = s1
    else:
        # Overlapping genes
        igr_start = min(e1, e2)
        igr_end = max(s1, s2)

    seq = genome[igr_start:igr_end]

    print(f"  {label}: {gene1_locus}({s1}-{e1}) -- {gene2_locus}({s2}-{e2})")
    print(f"    Intergenic: {igr_start}-{igr_end} ({len(seq)} bp)")

    return {
        'label': label,
        'gene1': gene1_locus,
        'gene2': gene2_locus,
        'start': igr_start,
        'end': igr_end,
        'length': len(seq),
        'sequence': seq,
    }


def search_motif_in_sequence(sequence, motif, max_mismatches=3):
    """Search for a motif in a sequence, allowing mismatches."""
    seq = sequence.upper()
    motif = motif.upper()
    rc_motif = str(Seq(motif).reverse_complement())
    hits = []

    for search_motif, strand in [(motif, '+'), (rc_motif, '-')]:
        mlen = len(search_motif)
        for i in range(len(seq) - mlen + 1):
            window = seq[i:i + mlen]
            mismatches = sum(a != b for a, b in zip(window, search_motif))
            if mismatches <= max_mismatches:
                hits.append({
                    'position': i,
                    'strand': strand,
                    'mismatches': mismatches,
                    'matched': window,
                    'identity': 1.0 - mismatches / mlen,
                })

    hits.sort(key=lambda h: h['mismatches'])
    return hits


def run_eram_validation():
    """Full Eram validation pipeline."""
    os.makedirs(RESULTS_DIR, exist_ok=True)

    print("[Eram Validation] Loading genome...")
    genome, record = load_genome()

    # === Step 1: Extract all 5 intergenic regions ===
    print("\n[Step 1] Extracting 5 intergenic regions...")
    regions = {}

    region_defs = [
        ('Region 1', 'Rv1963c', 'Rv1964'),     # mce3R - yrbE3A
        ('Region 2', 'Rv1935c', 'Rv1936'),      # echA13 - Rv1936
        ('Region 3', 'Rv1941', 'Rv1942c'),       # Rv1941 - Rv1942c
        ('Region 4', 'Rv1970', 'Rv1972'),        # Mce3F - Rv1972
        ('Region 5', 'Rv1944c', 'Rv1945'),       # Rv1944c - Rv1945
    ]

    for label, g1, g2 in region_defs:
        result = extract_intergenic(genome, record, g1, g2, label)
        if result:
            regions[label] = result

    # Write FASTA files for each region
    for label, reg in regions.items():
        fname = label.replace(' ', '_').lower() + '.fasta'
        fpath = os.path.join(RESULTS_DIR, fname)
        with open(fpath, 'w') as f:
            f.write(f">{label} [{reg['gene1']}-{reg['gene2']}] {reg['start']}-{reg['end']}\n")
            f.write(reg['sequence'] + '\n')

    # Write combined regions 1+2 (Eram's best motif came from this)
    if 'Region 1' in regions and 'Region 2' in regions:
        combined_path = os.path.join(RESULTS_DIR, 'regions_1_2_combined.fasta')
        with open(combined_path, 'w') as f:
            f.write(f">Regions1+2 [organism=Mycobacterium tuberculosis]\n")
            f.write(regions['Region 1']['sequence'] + regions['Region 2']['sequence'] + '\n')

    # Write all regions combined
    all_combined_path = os.path.join(RESULTS_DIR, 'all_regions_combined.fasta')
    with open(all_combined_path, 'w') as f:
        for label, reg in regions.items():
            f.write(f">{label}\n{reg['sequence']}\n")

    # === Step 2: Search for Eram's motifs ===
    print("\n[Step 2] Searching for Eram's motifs in intergenic regions...")

    motif_results = []

    for label, reg in regions.items():
        seq = reg['sequence']

        # Eram Motif 1 (49bp)
        hits = search_motif_in_sequence(seq, ERAM_MOTIF1_CONSENSUS, max_mismatches=8)
        for h in hits[:3]:
            motif_results.append({
                'region': label, 'motif': 'Eram_Motif1_49bp',
                'position': h['position'], 'strand': h['strand'],
                'mismatches': h['mismatches'], 'identity': h['identity'],
                'matched_seq': h['matched'],
            })
            print(f"  {label}: Eram Motif 1 at pos {h['position']} ({h['strand']}) "
                  f"identity={h['identity']:.1%} mismatches={h['mismatches']}")

        # Eram Motif 2 (21bp) — check each site variant
        for i, site in enumerate(ERAM_MOTIF2_SITES):
            hits = search_motif_in_sequence(seq, site, max_mismatches=4)
            for h in hits[:2]:
                motif_results.append({
                    'region': label, 'motif': f'Eram_Motif2_site{i+1}_21bp',
                    'position': h['position'], 'strand': h['strand'],
                    'mismatches': h['mismatches'], 'identity': h['identity'],
                    'matched_seq': h['matched'],
                })

        # Eram Motif 3 (10bp)
        hits = search_motif_in_sequence(seq, ERAM_MOTIF3_CONSENSUS, max_mismatches=2)
        for h in hits[:3]:
            motif_results.append({
                'region': label, 'motif': 'Eram_Motif3_10bp',
                'position': h['position'], 'strand': h['strand'],
                'mismatches': h['mismatches'], 'identity': h['identity'],
                'matched_seq': h['matched'],
            })

    # === Step 3: Search for Santangelo motifs ===
    print("\n[Step 3] Searching for Santangelo 2009 motifs...")

    for label, reg in regions.items():
        seq = reg['sequence']
        for i, motif in enumerate(SANTANGELO_FORWARD):
            hits = search_motif_in_sequence(seq, motif, max_mismatches=3)
            for h in hits[:2]:
                motif_results.append({
                    'region': label, 'motif': f'Santangelo_fwd_{i+1}',
                    'position': h['position'], 'strand': h['strand'],
                    'mismatches': h['mismatches'], 'identity': h['identity'],
                    'matched_seq': h['matched'],
                })
                if h['mismatches'] <= 1:
                    print(f"  {label}: Santangelo fwd motif {i+1} at pos {h['position']} "
                          f"identity={h['identity']:.1%}")

        for i, motif in enumerate(SANTANGELO_REVERSE):
            hits = search_motif_in_sequence(seq, motif, max_mismatches=3)
            for h in hits[:2]:
                motif_results.append({
                    'region': label, 'motif': f'Santangelo_rev_{i+1}',
                    'position': h['position'], 'strand': h['strand'],
                    'mismatches': h['mismatches'], 'identity': h['identity'],
                    'matched_seq': h['matched'],
                })

    # === Step 4: Check against Panagoda 2024 operator ===
    print("\n[Step 4] Mapping relative to Panagoda 2024 operator...")

    panagoda_operator = PARAMS['operator_sequence'].upper()
    panagoda_strong = panagoda_operator[88:123]  # Strong site (last 35bp)
    panagoda_weak = panagoda_operator[0:25]       # Weak site (first 25bp)

    print(f"  Panagoda operator: {len(panagoda_operator)} bp")
    print(f"  Strong site (88-123): {panagoda_strong[:25]}...")
    print(f"  Weak site (0-25):     {panagoda_weak}")

    # Find Panagoda operator in Region 1
    if 'Region 1' in regions:
        reg1_seq = regions['Region 1']['sequence']
        hits_op = search_motif_in_sequence(reg1_seq, panagoda_operator, max_mismatches=5)
        if hits_op:
            print(f"  Panagoda operator found in Region 1 at pos {hits_op[0]['position']} "
                  f"({hits_op[0]['strand']}) identity={hits_op[0]['identity']:.1%}")

        hits_strong = search_motif_in_sequence(reg1_seq, panagoda_strong[:21], max_mismatches=3)
        hits_weak = search_motif_in_sequence(reg1_seq, panagoda_weak[:21], max_mismatches=3)

        if hits_strong:
            print(f"  Panagoda strong site in Region 1 at pos {hits_strong[0]['position']}")
        if hits_weak:
            print(f"  Panagoda weak site in Region 1 at pos {hits_weak[0]['position']}")

    # === Step 5: Check FIMO results for echA13-Rv1936 region ===
    print("\n[Step 5] Checking FIMO hits near echA13-Rv1936 (Region 2)...")

    import pandas as pd
    fimo_path = os.path.join(PROJECT, 'results', 'phase1', 'predicted_sites.csv')
    if os.path.exists(fimo_path):
        fimo = pd.read_csv(fimo_path)
        # Region 2 genes
        region2_genes = ['Rv1935c', 'Rv1936', 'Rv1933c', 'Rv1934c']
        region2_hits = fimo[fimo['nearest_gene'].isin(region2_genes)]
        print(f"  FIMO hits near Region 2 genes: {len(region2_hits)}")
        for _, row in region2_hits.head(10).iterrows():
            print(f"    {row['nearest_gene']:10s} rank={int(row['rank']):4d} "
                  f"p={row['p_value']:.2e} score={row['score']:.2f} "
                  f"intergenic={row['is_intergenic']}")

    # === Save results ===
    import pandas as pd
    motif_df = pd.DataFrame(motif_results)
    motif_csv = os.path.join(RESULTS_DIR, 'motif_search_results.csv')
    motif_df.to_csv(motif_csv, index=False)
    print(f"\n  Saved: {motif_csv}")

    # Region summary
    region_rows = []
    for label, reg in regions.items():
        region_rows.append({
            'region': label,
            'gene1': reg['gene1'],
            'gene2': reg['gene2'],
            'start': reg['start'],
            'end': reg['end'],
            'length': reg['length'],
        })
    region_df = pd.DataFrame(region_rows)
    region_csv = os.path.join(RESULTS_DIR, 'intergenic_regions.csv')
    region_df.to_csv(region_csv, index=False)

    # === Step 6: Create DNA map ===
    print("\n[Step 6] Creating DNA map of motifs in Region 1...")
    create_dna_map(regions, motif_results)

    # Summary
    summary = {
        'n_regions': len(regions),
        'regions': {k: {'length': v['length'], 'gene1': v['gene1'], 'gene2': v['gene2']}
                    for k, v in regions.items()},
        'n_motif_hits': len(motif_results),
        'eram_motif1_found': any(m['motif'] == 'Eram_Motif1_49bp' for m in motif_results),
        'santangelo_found_region1': any(m['region'] == 'Region 1' and 'Santangelo' in m['motif']
                                        and m['mismatches'] <= 2 for m in motif_results),
        'santangelo_found_region2': any(m['region'] == 'Region 2' and 'Santangelo' in m['motif']
                                        and m['mismatches'] <= 2 for m in motif_results),
        'region2_fimo_hits': len(region2_hits) if 'region2_hits' in dir() else 0,
    }

    summary_path = os.path.join(RESULTS_DIR, 'eram_validation_summary.json')
    with open(summary_path, 'w') as f:
        json.dump(summary, f, indent=2, default=str)

    print(f"\n[Eram Validation] Complete. Results in {RESULTS_DIR}")
    return summary


def create_dna_map(regions, motif_results):
    """Create a visual DNA map showing motif positions in Region 1."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import matplotlib.patches as patches

    if 'Region 1' not in regions:
        return

    reg1 = regions['Region 1']
    seq_len = reg1['length']

    fig, ax = plt.subplots(1, 1, figsize=(16, 8))
    ax.set_xlim(-10, seq_len + 10)
    ax.set_ylim(-1, 10)

    # Draw the DNA backbone
    ax.plot([0, seq_len], [5, 5], 'k-', linewidth=3)
    ax.text(seq_len / 2, 5.4, f"Region 1: {reg1['gene1']} - {reg1['gene2']} ({seq_len} bp)",
            ha='center', fontsize=12, fontweight='bold')

    # Gene labels
    ax.annotate(reg1['gene1'], xy=(0, 5), xytext=(-5, 6.5),
                fontsize=10, ha='right', arrowprops=dict(arrowstyle='->', color='gray'))
    ax.annotate(reg1['gene2'], xy=(seq_len, 5), xytext=(seq_len + 5, 6.5),
                fontsize=10, ha='left', arrowprops=dict(arrowstyle='->', color='gray'))

    # Panagoda operator position in Region 1
    panagoda_op = PARAMS['operator_sequence'].upper()
    hits = search_motif_in_sequence(reg1['sequence'], panagoda_op, max_mismatches=10)
    if hits:
        pos = hits[0]['position']
        rect = patches.Rectangle((pos, 4.2), len(panagoda_op), 0.6, linewidth=2,
                                  edgecolor='black', facecolor='#FFD700', alpha=0.7)
        ax.add_patch(rect)
        ax.text(pos + len(panagoda_op) / 2, 3.8, 'Panagoda 2024\noperator (123bp)',
                ha='center', fontsize=8, color='#B8860B')

    # Color scheme for motifs
    motif_colors = {
        'Eram_Motif1_49bp': '#FF4444',
        'Eram_Motif2': '#4444FF',
        'Eram_Motif3_10bp': '#44AA44',
        'Santangelo_fwd': '#FF8800',
        'Santangelo_rev': '#AA44FF',
    }

    y_positions = {
        'Eram_Motif1_49bp': 7,
        'Eram_Motif2': 8,
        'Eram_Motif3_10bp': 8.5,
        'Santangelo_fwd': 2,
        'Santangelo_rev': 1,
    }

    # Plot motifs in Region 1
    reg1_motifs = [m for m in motif_results if m['region'] == 'Region 1']
    plotted = set()

    for m in reg1_motifs:
        # Determine color and y position
        motif_base = m['motif'].rsplit('_', 1)[0] if 'site' in m['motif'] else m['motif']
        for key in motif_colors:
            if key in m['motif']:
                color = motif_colors[key]
                y = y_positions.get(key, 7)
                break
        else:
            color = 'gray'
            y = 9

        pos = m['position']
        width = len(m['matched_seq'])
        plot_key = (m['motif'], pos)
        if plot_key in plotted:
            continue
        plotted.add(plot_key)

        if m['mismatches'] <= 4:  # Only show good matches
            rect = patches.Rectangle((pos, y - 0.2), width, 0.4, linewidth=1,
                                      edgecolor=color, facecolor=color, alpha=0.5)
            ax.add_patch(rect)
            if m['mismatches'] <= 2:
                ax.text(pos + width / 2, y + 0.4, f"{m['motif']}\n({m['strand']}, {m['mismatches']}mm)",
                        ha='center', fontsize=6, color=color)

    # Legend
    legend_items = [
        patches.Patch(color='#FFD700', label='Panagoda 2024 operator (123bp)'),
        patches.Patch(color='#FF4444', label='Eram Motif 1 (49bp)'),
        patches.Patch(color='#4444FF', label='Eram Motif 2 (21bp)'),
        patches.Patch(color='#44AA44', label='Eram Motif 3 (10bp)'),
        patches.Patch(color='#FF8800', label='Santangelo forward motifs'),
        patches.Patch(color='#AA44FF', label='Santangelo reverse motifs'),
    ]
    ax.legend(handles=legend_items, loc='lower right', fontsize=8)

    ax.set_xlabel('Position in intergenic region (bp)', fontsize=11)
    ax.set_title('DNA Map: Motif Locations in mce3R-yrbE3A Intergenic Region\n'
                 'Relative to Panagoda 2024 Operator', fontsize=13, fontweight='bold')
    ax.set_yticks([])

    plt.tight_layout()
    map_path = os.path.join(RESULTS_DIR, 'dna_motif_map.png')
    fig.savefig(map_path, dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved: {map_path}")


if __name__ == '__main__':
    run_eram_validation()
