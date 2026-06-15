"""
Zoomed, annotated DNA map of the Mce3R motif clusters.

Addresses Dr. Balazsi's slideshow feedback (2026-05-18):
  1. "Zoom into all of these regions where motifs seem to exist."
  2. "Add to the annotation all 3 suspected motifs (ABC, DEF and HXG)."

The 99 bp de novo motif (extended MEME, regions_1_2_combined, MEME-1,
E = 3.9e-08) occurs three times. Each occurrence corresponds to one of the
three motif clusters in Balazsi's annotated schematic (reference_docs/
TBmotifs.pdf, slide 2):

  ABC  = Region 1 (mce3R-yrbE3A) proximal occurrence  -> 99 bp @ R1 pos 238 (+)
  DEF  = Region 1 (mce3R-yrbE3A) distal occurrence     -> 99 bp @ R1 pos 684 (+)
  HXG  = Region 2 (echA13-Rv1936) occurrence           -> 99 bp @ R2 pos  62 (-)

Coordinates verified against results/extended_meme/report/meme_report.json
and the MEME block diagram 237_[+1]_347_[+1]_176_[-1]_64
(237+99+347+99+176+99+64 = 1121 = Region1 897 bp + Region2 224 bp).

Sub-site core sequences (A-F, G/X/H) are the Santangelo 2009 / Bigi motif
instances catalogued in eram_validation.py and Eram Kabir's appendix (p.17).
Equivalence per TBmotifs.pdf:  A = C = D = F = rcH = rcG ,  B = D = X .

Every sub-site is LOCATED BY SEQUENCE MATCH (both strands, <=3 mismatches)
inside its 99 bp window -- no hand-coded coordinates. The script prints a
full audit table and fails loudly if any expected site is not found.
"""

import os
import sys
import json

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

from Bio.Seq import Seq
from dna_features_viewer import GraphicFeature, GraphicRecord

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from config.parameters import PARAMS

PROJECT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ERAM_DIR = os.path.join(PROJECT, 'results', 'eram_validation')
REPORT_DIR = os.path.join(PROJECT, 'results', 'extended_meme', 'report')
os.makedirs(REPORT_DIR, exist_ok=True)

# Eram's published 49 bp motif (his Fig 6) -- the gold-standard core.
ERAM_MOTIF1 = 'ACTAGCAAGATACATCATAGCCAATATATGCCAGTTTGCATTGCTATTT'

# ---- Deck-consistent palette (white bg, blue/teal) ------------------------
C_GENE      = '#5b6770'   # flanking genes
C_PROMOTER  = '#8c6bb1'   # echA13 promoter
C_OPERATOR  = '#fdae61'   # Panagoda 2024 123 bp operator
C_CLUSTER   = '#2c7fb8'   # 99 bp de novo motif occurrence (the cluster block)
C_FIRST     = '#1b9e77'   # first-half sub-site (A/C/D/F family)
C_SECOND    = '#7570b3'   # second-half sub-site (B/E/X family)
C_ERAM49    = '#253494'   # Eram's 49 bp core (the 100%-identity match)

# ---- The three 99 bp occurrences (region-relative, 0-based) ----------------
# Verified from meme_report.json (combined 1-based starts 238 / 684 / 959).
# Combined->region: Region1 = 0..896, Region2 = combined-897.
OCC_WIDTH = 99
OCCURRENCES = {
    'ABC': dict(region=1, start0=237, strand=+1, pval='2.3e-44',
                eram='73%', desc='Region 1 proximal (near mce3R)'),
    'DEF': dict(region=1, start0=683, strand=+1, pval='1.6e-45',
                eram='100%', desc='Region 1 distal (near yrbE3A)'),
    'HXG': dict(region=2, start0=61,  strand=-1, pval='4.2e-45',
                eram='67%', desc='Region 2 (echA13-Rv1936)'),
}

# ---- Sub-site cores (Santangelo/Bigi instances; Eram appendix p.17) --------
# cluster -> list of (label, core_sequence, role)
#   role 'first'  = first half of Eram's motif  (A/C/D/F/G/H family)
#   role 'second' = second half of Eram's motif (B/E/X family)
SUBSITES = {
    'ABC': [
        ('A', 'GACTAACAAAATACAT',  'first'),
        ('B', 'GATTACATTGCAATTTA', 'second'),
        ('C', 'TATCAGTAATAGACAT',  'first'),
    ],
    'DEF': [
        ('D', 'TACTAGCAAGATACAT',    'first'),
        ('E', 'CAGTTTGCATTGCTATTTA', 'second'),
        ('F', 'TATTGGCTATGGACAT',    'first'),
    ],
    'HXG': [
        ('H', 'ACTACCAAGATACTT',      'first'),
        ('X', 'CAAATTCCCATGCAAAGAAG', 'second'),
        ('G', 'ATTGGCTATGGACAT',      'first'),
    ],
}


def load_region(fasta_name):
    """Return the bare sequence string (skip the FASTA header)."""
    with open(os.path.join(ERAM_DIR, fasta_name)) as fh:
        lines = [ln.strip() for ln in fh if ln.strip()]
    return ''.join(lines[1:]).upper()


def best_match(haystack, needle, lo, hi):
    """Best (min-mismatch) placement of `needle` (either strand) in
    haystack[lo:hi]. Returns (start0, end0, strand, mismatches, matched_seq)."""
    rc = str(Seq(needle).reverse_complement())
    best = None
    for probe, strand in ((needle, +1), (rc, -1)):
        n = len(probe)
        for i in range(max(0, lo), min(len(haystack), hi) - n + 1):
            window = haystack[i:i + n]
            mm = sum(a != b for a, b in zip(window, probe))
            if best is None or mm < best[3]:
                best = (i, i + n, strand, mm, window)
    return best


def locate(seq, region_no, audit):
    """Locate every 99 bp occurrence and sub-site for one region."""
    placed = {}
    for cluster, occ in OCCURRENCES.items():
        if occ['region'] != region_no:
            continue
        o0 = occ['start0']
        o1 = o0 + OCC_WIDTH
        placed[cluster] = {
            'occ': (o0, o1, occ['strand']),
            'subsites': [],
        }
        # Search sub-sites within the occurrence +/- 12 bp slack.
        lo, hi = o0 - 12, o1 + 12
        for label, core, role in SUBSITES[cluster]:
            hit = best_match(seq, core, lo, hi)
            if hit is None or hit[3] > 3:
                raise RuntimeError(
                    f"Sub-site {label} ({core}) not found in Region {region_no} "
                    f"window [{lo},{hi}] (best mismatches="
                    f"{None if hit is None else hit[3]}). Refusing to draw a "
                    f"figure with an unverified annotation.")
            s0, e0, strand, mm, matched = hit
            placed[cluster]['subsites'].append(
                dict(label=label, role=role, start0=s0, end0=e0,
                     strand=strand, mm=mm, matched=matched))
            audit.append(dict(region=region_no, cluster=cluster, feature=label,
                              start0=s0, end0=e0,
                              strand='+' if strand > 0 else '-',
                              mismatches=mm, seq=matched))
        audit.append(dict(region=region_no, cluster=cluster,
                          feature=f'{cluster} 99bp occ',
                          start0=o0, end0=o1,
                          strand='+' if occ['strand'] > 0 else '-',
                          mismatches=0, seq=f"p={occ['pval']}"))
    return placed


def gene_features(region_no, seqlen):
    """Flanking-gene arrows + promoter/operator landmarks."""
    feats = []
    if region_no == 1:
        feats.append(GraphicFeature(start=0, end=22, strand=-1, color=C_GENE,
                                     label='mce3R (Rv1963c)'))
        feats.append(GraphicFeature(start=seqlen - 22, end=seqlen, strand=+1,
                                     color=C_GENE, label='yrbE3A (Rv1964)'))
    else:
        feats.append(GraphicFeature(start=0, end=18, strand=-1, color=C_GENE,
                                     label='echA13 (Rv1935c)'))
        feats.append(GraphicFeature(start=seqlen - 18, end=seqlen, strand=+1,
                                     color=C_GENE, label='Rv1936'))
    return feats


def cluster_features(placed, eram_span=None, short=False):
    """Build GraphicFeatures for the 99 bp clusters + sub-sites.

    short=True  -> compact labels for the overview panels (no collisions).
    short=False -> full labels for the zoom panels.
    """
    feats = []
    for cluster, info in placed.items():
        o0, o1, ostr = info['occ']
        occ = OCCURRENCES[cluster]
        if short:
            clabel = f"{cluster} · 99 bp · p={occ['pval']}"
        else:
            clabel = (f"{cluster}  (99 bp de novo motif, p={occ['pval']}, "
                      f"Eram {occ['eram']})")
        feats.append(GraphicFeature(start=o0, end=o1, strand=ostr,
                                    color=C_CLUSTER, label=clabel))
        for ss in info['subsites']:
            col = C_FIRST if ss['role'] == 'first' else C_SECOND
            feats.append(GraphicFeature(
                start=ss['start0'], end=ss['end0'], strand=ss['strand'],
                color=col, label=ss['label']))
    if eram_span is not None:
        s0, e0, strand = eram_span
        feats.append(GraphicFeature(start=s0, end=e0, strand=strand,
                                    color=C_ERAM49,
                                    label="Eram 49 bp core (100% identity)"))
    return feats


def plot_panel(ax, record, crop, title, with_seq):
    """Plot a (optionally cropped) GraphicRecord onto ax."""
    rec = record.crop(crop) if crop else record
    if with_seq:
        rec.plot(ax=ax, with_ruler=True)
        rec.plot_sequence(ax=ax)
    else:
        rec.plot(ax=ax, with_ruler=True)
    ax.set_title(title, fontsize=12, fontweight='bold', loc='left', pad=8)


def main():
    seq1 = load_region('region_1.fasta')
    seq2 = load_region('region_2.fasta')
    len1, len2 = len(seq1), len(seq2)
    print(f"Region 1 (mce3R-yrbE3A): {len1} bp")
    print(f"Region 2 (echA13-Rv1936): {len2} bp")
    assert len1 == 897 and len2 == 224, \
        f"Unexpected region lengths {len1}/{len2}; coordinate mapping invalid."

    audit = []
    placed1 = locate(seq1, 1, audit)
    placed2 = locate(seq2, 2, audit)

    # Locate Eram's 49 bp core in Region 1 (must fall inside DEF, + strand).
    eram_hit = best_match(seq1, ERAM_MOTIF1, 0, len1)
    es0, ee0, estr, emm, _ = eram_hit
    print(f"Eram 49 bp core: Region 1 {es0}-{ee0} "
          f"({'+' if estr > 0 else '-'}), mismatches={emm}")
    assert emm == 0, f"Eram 49 bp core not exact in Region 1 (mm={emm})."
    def0 = placed1['DEF']['occ']
    assert def0[0] - 5 <= es0 and ee0 <= def0[1] + 5, \
        "Eram 49 bp core is not inside the DEF 99 bp occurrence."
    eram_span = (es0, ee0, estr)

    # Panagoda 2024 operator in Region 1 (landmark, best-effort).
    operator = PARAMS['operator_sequence'].upper()
    op_hit = best_match(seq1, operator, 0, len1)
    op_feat = None
    if op_hit and op_hit[3] <= max(8, int(0.15 * len(operator))):
        op_feat = GraphicFeature(start=op_hit[0], end=op_hit[1],
                                  strand=op_hit[2], color=C_OPERATOR,
                                  label='Panagoda 2024 operator (123 bp)')
        print(f"Panagoda operator: Region 1 {op_hit[0]}-{op_hit[1]} "
              f"mm={op_hit[3]}")

    # echA13 promoter spans the HXG cluster (Balazsi schematic).
    hxg = placed2['HXG']['occ']
    promoter = GraphicFeature(start=max(0, hxg[0] - 8),
                              end=min(len2, hxg[1] + 8), strand=0,
                              color=C_PROMOTER, label='echA13 promoter')

    # ---- Build per-region GraphicRecords ----
    # Overview: compact labels, no Eram-core / operator clutter.
    rec1_ov = GraphicRecord(
        sequence=seq1, sequence_length=len1,
        features=gene_features(1, len1)
        + cluster_features(placed1, short=True))
    rec2_ov = GraphicRecord(
        sequence=seq2, sequence_length=len2,
        features=gene_features(2, len2)
        + cluster_features(placed2, short=True) + [promoter])

    # Detail (for zoom crops): full labels + Eram core + operator + promoter.
    r1_detail = gene_features(1, len1) + cluster_features(placed1, eram_span)
    if op_feat is not None:
        r1_detail.append(op_feat)
    rec1 = GraphicRecord(sequence=seq1, sequence_length=len1,
                         features=r1_detail)
    rec2 = GraphicRecord(
        sequence=seq2, sequence_length=len2,
        features=gene_features(2, len2)
        + cluster_features(placed2) + [promoter])

    # Zoom windows (occurrence +/- 22 bp, clipped).
    def zoom(placed, cluster, seqlen):
        o0, o1, _ = placed[cluster]['occ']
        return (max(0, o0 - 22), min(seqlen, o1 + 22))

    z_abc = zoom(placed1, 'ABC', len1)
    z_def = zoom(placed1, 'DEF', len1)
    z_hxg = zoom(placed2, 'HXG', len2)

    # ---- Figure: 2 overview panels + 3 zoom panels ----
    fig = plt.figure(figsize=(17, 21), facecolor='white')
    gs = GridSpec(5, 1, height_ratios=[1.15, 1.0, 1.3, 1.3, 1.3],
                  hspace=0.75)

    ax0 = fig.add_subplot(gs[0])
    plot_panel(ax0, rec1_ov, None,
               "Region 1  —  mce3R–yrbE3A intergenic (897 bp): "
               "ABC and DEF clusters", with_seq=False)

    ax1 = fig.add_subplot(gs[1])
    plot_panel(ax1, rec2_ov, None,
               "Region 2  —  echA13–Rv1936 intergenic (224 bp): "
               "HXG cluster (− strand)", with_seq=False)

    ax2 = fig.add_subplot(gs[2])
    plot_panel(ax2, rec1, z_abc,
               f"ZOOM — Cluster ABC  (Region 1, 99 bp motif @ pos 238 +, "
               f"p = 2.3e-44; sites A, B, C)", with_seq=True)

    ax3 = fig.add_subplot(gs[3])
    plot_panel(ax3, rec1, z_def,
               f"ZOOM — Cluster DEF  (Region 1, 99 bp motif @ pos 684 +, "
               f"p = 1.6e-45; contains Eram 49 bp @ 100%; sites D, E, F)",
               with_seq=True)

    ax4 = fig.add_subplot(gs[4])
    plot_panel(ax4, rec2, z_hxg,
               f"ZOOM — Cluster HXG  (Region 2, 99 bp motif @ pos 62 −, "
               f"p = 4.2e-45; sites H, X, G)", with_seq=True)

    fig.suptitle(
        "Mce3R binding-site map: three 99 bp de novo motif clusters "
        "(ABC, DEF, HXG)\nde novo extended MEME, Regions 1+2 combined, "
        "E = 3.9e-08  ·  sub-sites A–F / G·X·H per Balázsi (TBmotifs.pdf)",
        fontsize=14, fontweight='bold', y=0.998)

    import matplotlib.patches as mpatches
    legend_handles = [
        mpatches.Patch(color=C_CLUSTER, label='99 bp de novo motif (cluster)'),
        mpatches.Patch(color=C_FIRST,
                       label="first-half sub-site (A/C/D/F · G/H)"),
        mpatches.Patch(color=C_SECOND, label='second-half sub-site (B/E · X)'),
        mpatches.Patch(color=C_ERAM49, label='Eram 49 bp core (100% identity)'),
        mpatches.Patch(color=C_OPERATOR,
                       label='Panagoda 2024 operator (123 bp)'),
        mpatches.Patch(color=C_PROMOTER, label='echA13 promoter'),
        mpatches.Patch(color=C_GENE, label='flanking gene'),
    ]
    fig.legend(handles=legend_handles, loc='lower center', ncol=4,
               fontsize=10, frameon=True, bbox_to_anchor=(0.5, -0.012))

    out_primary = os.path.join(REPORT_DIR, 'dna_motif_map_zoomed.png')
    out_mirror = os.path.join(ERAM_DIR, 'dna_motif_map_zoomed.png')
    fig.savefig(out_primary, dpi=300, bbox_inches='tight', facecolor='white')
    fig.savefig(out_mirror, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print(f"\nSaved: {out_primary}")
    print(f"Saved: {out_mirror}")

    # ---- Audit table + JSON ----
    print("\n" + "=" * 78)
    print("ANNOTATION AUDIT (region-relative 0-based coords)")
    print("=" * 78)
    print(f"{'Reg':>3} {'Cluster':>7} {'Feat':>12} {'start':>6} {'end':>6} "
          f"{'str':>3} {'mm':>3}  seq/notes")
    for r in audit:
        print(f"{r['region']:>3} {r['cluster']:>7} {r['feature']:>12} "
              f"{r['start0']:>6} {r['end0']:>6} {r['strand']:>3} "
              f"{r['mismatches']:>3}  {r['seq']}")

    # Measured gaps between consecutive sub-sites (data-derived, not copied
    # from Balazsi's schematic which uses different genomic windows).
    gaps = {}
    for cluster in ('ABC', 'DEF', 'HXG'):
        src = placed1 if cluster in placed1 else placed2
        ss = sorted(src[cluster]['subsites'], key=lambda d: d['start0'])
        g = [ss[i + 1]['start0'] - ss[i]['end0'] for i in range(len(ss) - 1)]
        # Report in binding-strand 5'->3' order: for a minus-strand
        # occurrence that is high->low region coordinate (so HXG reads
        # H, X, G -- matching Balazsi's cluster name).
        if OCCURRENCES[cluster]['strand'] < 0:
            ss = ss[::-1]
            g = g[::-1]
        gaps[cluster] = {
            'order': [s['label'] for s in ss],
            'gaps_bp': g,
            'order_along': "5'->3' on binding strand",
        }
    print("\nMeasured inter-site gaps (bp, region-relative):")
    for k, v in gaps.items():
        print(f"  {k}: {' '.join(v['order'])}  gaps={v['gaps_bp']}")

    summary = {
        'description': 'Zoomed, ABC/DEF/HXG-annotated Mce3R motif map '
                       '(Balazsi 2026-05-18 slideshow feedback).',
        'motif': '99 bp de novo (extended MEME, regions_1_2_combined, '
                 'MEME-1, E=3.9e-08)',
        'region_lengths': {'Region1_mce3R_yrbE3A': len1,
                           'Region2_echA13_Rv1936': len2},
        'clusters': {
            c: {
                'region': OCCURRENCES[c]['region'],
                'occurrence_start0': OCCURRENCES[c]['start0'],
                'occurrence_end0': OCCURRENCES[c]['start0'] + OCC_WIDTH,
                'strand': '+' if OCCURRENCES[c]['strand'] > 0 else '-',
                'site_pvalue': OCCURRENCES[c]['pval'],
                'eram_49bp_identity': OCCURRENCES[c]['eram'],
                'description': OCCURRENCES[c]['desc'],
                'subsites': (placed1 if OCCURRENCES[c]['region'] == 1
                             else placed2)[c]['subsites'],
                'inter_site_gaps_bp': gaps[c],
            } for c in OCCURRENCES
        },
        'eram_49bp_core_region1': {'start0': es0, 'end0': ee0,
                                   'strand': '+' if estr > 0 else '-',
                                   'mismatches': emm},
        'outputs': [out_primary, out_mirror],
    }
    out_json = os.path.join(REPORT_DIR, 'zoomed_motif_map.json')
    with open(out_json, 'w') as fh:
        json.dump(summary, fh, indent=2, default=str)
    print(f"\nSummary: {out_json}")
    return summary


if __name__ == '__main__':
    main()
