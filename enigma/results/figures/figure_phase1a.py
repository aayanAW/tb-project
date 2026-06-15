#!/usr/bin/env python3
"""
Figure: Phase 1A — Data Acquisition & Sequence Preparation
Clear workflow diagram showing what data goes in and what comes out.
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import numpy as np

fig, ax = plt.subplots(figsize=(16, 20))
ax.set_xlim(0, 16)
ax.set_ylim(0, 22)
ax.axis('off')
fig.patch.set_facecolor('white')

# ============================================================
# COLORS
# ============================================================
C_TB    = '#1a5276'   # M. tuberculosis
C_BOV   = '#2e86c1'   # M. bovis
C_MAR   = '#5dade2'   # M. marinum
C_NCBI  = '#117a65'   # NCBI green
C_ARROW = '#2c3e50'
C_OUT   = '#8e44ad'   # output files
C_OP    = '#c0392b'   # operator
C_GENE  = '#d4ac0d'   # gene annotations
C_BG    = '#f8f9fa'

def box(ax, x, y, w, h, text, color, fontsize=11, subtext=None, subtextsize=9,
        alpha=0.12, lw=2, textcolor=None, bold=True):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.15",
        facecolor=color, alpha=alpha, edgecolor='none'))
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.15",
        facecolor='none', edgecolor=color, lw=lw))
    tc = textcolor or color
    weight = 'bold' if bold else 'normal'
    if subtext:
        ax.text(x + w/2, y + h/2 + 0.2, text, ha='center', va='center',
                fontsize=fontsize, fontweight=weight, color=tc, linespacing=1.3)
        ax.text(x + w/2, y + h/2 - 0.35, subtext, ha='center', va='center',
                fontsize=subtextsize, color='#566573', style='italic', linespacing=1.3)
    else:
        ax.text(x + w/2, y + h/2, text, ha='center', va='center',
                fontsize=fontsize, fontweight=weight, color=tc, linespacing=1.3)

def arrow_down(ax, x, y1, y2, color=C_ARROW):
    ax.annotate('', xy=(x, y2), xytext=(x, y1),
                arrowprops=dict(arrowstyle='->', color=color, lw=2.5,
                               connectionstyle='arc3,rad=0'))

def arrow_down_curved(ax, x1, y1, x2, y2, color=C_ARROW, rad=0.2):
    ax.annotate('', xy=(x2, y2), xytext=(x1, y1),
                arrowprops=dict(arrowstyle='->', color=color, lw=2,
                               connectionstyle=f'arc3,rad={rad}'))

# ============================================================
# TITLE
# ============================================================
ax.text(8, 21.5, 'Phase 1A: Data Acquisition & Sequence Preparation',
        ha='center', va='center', fontsize=18, fontweight='bold',
        color=C_ARROW, fontfamily='sans-serif')
ax.text(8, 21.0, 'Inputs, processing steps, and outputs for the Mce3R operator discovery pipeline',
        ha='center', va='center', fontsize=11, color='#7f8c8d', style='italic')

# ============================================================
# ROW 1: NCBI Downloads (3 genomes)
# ============================================================
# Section label
ax.text(0.5, 20.2, 'Step 1', fontsize=13, fontweight='bold', color=C_NCBI,
        ha='left', va='center')
ax.text(2.2, 20.2, 'Download Reference Genomes from NCBI', fontsize=13,
        fontweight='bold', color=C_ARROW, ha='left', va='center')
ax.text(14.5, 20.2, 'download_genomes.py', fontsize=9, color='#7f8c8d',
        ha='right', va='center', style='italic', fontfamily='monospace')

# NCBI database
box(ax, 6.0, 18.8, 4.0, 1.1, 'NCBI GenBank', C_NCBI, fontsize=13,
    subtext='Public sequence database')

# Arrow down to 3 genomes
arrow_down(ax, 8.0, 18.8, 18.3)

# 3 genome boxes
gw, gh = 4.2, 1.5
box(ax, 0.8, 16.5, gw, gh,
    'M. tuberculosis H37Rv', C_TB, fontsize=11,
    subtext='NC_000962.3  •  4,411,532 bp  •  65.6% GC')
box(ax, 5.9, 16.5, gw, gh,
    'M. bovis AF2122/97', C_BOV, fontsize=11,
    subtext='NC_002945.4  •  4,345,492 bp')
box(ax, 11.0, 16.5, gw, gh,
    'M. marinum M', C_MAR, fontsize=11,
    subtext='NC_010612.1  •  6,636,827 bp')

# Arrows from NCBI to each genome
arrow_down_curved(ax, 8.0, 18.8, 2.9, 18.0, C_NCBI, rad=0.3)
arrow_down(ax, 8.0, 18.8, 18.0)
arrow_down_curved(ax, 8.0, 18.8, 13.1, 18.0, C_NCBI, rad=-0.3)

# Also download GFF annotation
box(ax, 11.0, 15.0, gw, 0.7, 'H37Rv Gene Annotations (GFF3)', C_GENE,
    fontsize=9, bold=True)
ax.text(13.1, 14.55, '3,978 genes with coordinates & strand', ha='center',
        fontsize=8, color='#7f8c8d', style='italic')

arrow_down_curved(ax, 8.0, 18.8, 13.1, 15.7, C_GENE, rad=-0.4)

# File outputs
ax.text(2.9, 16.15, 'H37Rv.fasta', ha='center', fontsize=8,
        fontfamily='monospace', color=C_TB)
ax.text(8.0, 16.15, 'M_bovis.fasta', ha='center', fontsize=8,
        fontfamily='monospace', color=C_BOV)
ax.text(13.1, 16.15, 'M_marinum.fasta', ha='center', fontsize=8,
        fontfamily='monospace', color=C_MAR)

# ============================================================
# ROW 2: Extract Upstream Regions
# ============================================================
# Divider
ax.plot([0.5, 15.5], [13.8, 13.8], color='#d5d8dc', lw=1)

ax.text(0.5, 13.3, 'Step 2', fontsize=13, fontweight='bold', color='#d35400',
        ha='left', va='center')
ax.text(2.2, 13.3, 'Extract 200 bp Upstream of Every Gene', fontsize=13,
        fontweight='bold', color=C_ARROW, ha='left', va='center')
ax.text(14.5, 13.3, 'extract_upstream.py', fontsize=9, color='#7f8c8d',
        ha='right', va='center', style='italic', fontfamily='monospace')

# Genome diagram showing upstream extraction
# Draw a gene arrow with upstream region
gene_y = 12.0
# Chromosome line
ax.plot([1.5, 14.5], [gene_y, gene_y], color='#2c3e50', lw=3)

# Gene body (arrow shape)
gene_start = 7.0
gene_end = 11.0
gene_h = 0.5
arrow_tip = 0.4
gene_body = plt.Polygon([
    (gene_start, gene_y - gene_h/2),
    (gene_end - arrow_tip, gene_y - gene_h/2),
    (gene_end, gene_y),
    (gene_end - arrow_tip, gene_y + gene_h/2),
    (gene_start, gene_y + gene_h/2),
], facecolor=C_GENE, alpha=0.4, edgecolor=C_GENE, lw=1.5)
ax.add_patch(gene_body)
ax.text(8.8, gene_y, 'Gene (CDS)', ha='center', va='center',
        fontsize=9, fontweight='bold', color='#7d6608')

# Upstream region (highlighted)
up_start = 5.0
up_end = gene_start
ax.add_patch(FancyBboxPatch((up_start, gene_y - 0.4), up_end - up_start, 0.8,
    boxstyle="round,pad=0.05", facecolor='#d35400', alpha=0.2,
    edgecolor='#d35400', lw=2, linestyle='--'))
ax.text((up_start + up_end)/2, gene_y, '200 bp\nupstream', ha='center',
        va='center', fontsize=9, fontweight='bold', color='#d35400')

# Annotation arrows
ax.annotate('', xy=(up_start, gene_y - 0.7), xytext=(up_end, gene_y - 0.7),
            arrowprops=dict(arrowstyle='<->', color='#d35400', lw=1.5))
ax.text((up_start + up_end)/2, gene_y - 1.0, '200 bp', ha='center',
        fontsize=9, color='#d35400', fontweight='bold')

# TSS label
ax.plot([gene_start, gene_start], [gene_y - 0.6, gene_y + 0.7], color='#2c3e50',
        lw=1, ls=':')
ax.text(gene_start, gene_y + 0.85, 'Start', ha='center', fontsize=8,
        color='#2c3e50')

# "Repeat for all 3,978 genes"
ax.text(13.0, gene_y + 0.6, 'Repeat for all\n3,978 genes', ha='center',
        va='center', fontsize=9, color='#566573', style='italic',
        bbox=dict(boxstyle='round,pad=0.3', facecolor='#f2f3f4', edgecolor='#d5d8dc'))

# Circular chromosome note
ax.text(2.5, gene_y + 0.6, 'Handles circular\nchromosome wrapping', ha='center',
        va='center', fontsize=8, color='#566573', style='italic',
        bbox=dict(boxstyle='round,pad=0.3', facecolor='#f2f3f4', edgecolor='#d5d8dc'))

# ============================================================
# ROW 3: Build Ortholog Input for MEME
# ============================================================
ax.plot([0.5, 15.5], [10.0, 10.0], color='#d5d8dc', lw=1)

ax.text(0.5, 9.5, 'Step 3', fontsize=13, fontweight='bold', color='#2980b9',
        ha='left', va='center')
ax.text(2.2, 9.5, 'Build Ortholog Upstream Sequences for MEME Input', fontsize=13,
        fontweight='bold', color=C_ARROW, ha='left', va='center')

# Three species upstream boxes converging
ow, oh = 3.5, 1.0
oy = 7.8

box(ax, 0.8, oy, ow, oh, 'H37Rv upstream\nof yrbE3A (Rv1964)', C_TB, fontsize=9)
box(ax, 6.2, oy, ow, oh, 'M. bovis upstream\nof Mb1997', C_BOV, fontsize=9)
box(ax, 11.6, oy, ow, oh, 'M. marinum upstream\nof MMAR_2522', C_MAR, fontsize=9)

# Arrows converging to output
arrow_down_curved(ax, 2.55, 7.8, 8.0, 7.0, C_TB, rad=0.3)
arrow_down(ax, 7.95, 7.8, 7.0)
arrow_down_curved(ax, 13.35, 7.8, 8.0, 7.0, C_MAR, rad=-0.3)

# "200bp each" labels
ax.text(2.55, 8.95, '200 bp', ha='center', fontsize=8, color=C_TB, fontweight='bold')
ax.text(7.95, 8.95, '200 bp', ha='center', fontsize=8, color=C_BOV, fontweight='bold')
ax.text(13.35, 8.95, '200 bp', ha='center', fontsize=8, color=C_MAR, fontweight='bold')

# ============================================================
# ROW 4: Output Files
# ============================================================
ax.plot([0.5, 15.5], [6.6, 6.6], color='#d5d8dc', lw=1)

ax.text(8.0, 6.2, 'OUTPUT FILES', ha='center', fontsize=14, fontweight='bold',
        color=C_OUT)

# Three output files
ow2, oh2 = 4.2, 1.8

# File 1: all upstream
box(ax, 0.5, 3.8, ow2, oh2, 'all_upstream_200bp.fasta', C_OUT, fontsize=10,
    subtext='~3,978 sequences\n200 bp upstream of every H37Rv gene')

# File 2: MEME input
box(ax, 5.9, 3.8, ow2, oh2, 'meme_input_orthologs.fasta', C_OUT, fontsize=10,
    subtext='3 sequences\nOrtholog upstream regions for MEME')

# File 3: known operator
box(ax, 11.3, 3.8, ow2, oh2, 'known_operator.fasta', C_OP, fontsize=10,
    subtext='1 sequence\n123 bp experimentally verified operator')

# Arrows down to outputs
arrow_down_curved(ax, 8.0, 6.6, 2.6, 5.6, C_OUT, rad=0.2)
arrow_down(ax, 8.0, 6.6, 5.6)
arrow_down_curved(ax, 8.0, 6.6, 13.4, 5.6, C_OP, rad=-0.2)

# ============================================================
# ROW 5: Where outputs go next
# ============================================================
ax.plot([0.5, 15.5], [3.4, 3.4], color='#d5d8dc', lw=1)

# Next steps boxes (faded, showing what consumes these files)
nw, nh = 4.2, 1.2

box(ax, 0.5, 1.6, nw, nh, 'Phase 1C\nBipartite PWM Search', '#8e44ad',
    fontsize=10, alpha=0.08, lw=1.5)
ax.text(2.6, 1.25, 'scans all upstream regions', ha='center', fontsize=8,
        color='#7f8c8d', style='italic')

box(ax, 5.9, 1.6, nw, nh, 'Phase 1B\nMEME Motif Discovery', '#2980b9',
    fontsize=10, alpha=0.08, lw=1.5)
ax.text(8.0, 1.25, 'discovers conserved motifs', ha='center', fontsize=8,
        color='#7f8c8d', style='italic')

box(ax, 11.3, 1.6, nw, nh, 'Phase 1C\nPWM Construction', '#c0392b',
    fontsize=10, alpha=0.08, lw=1.5)
ax.text(13.4, 1.25, 'builds search templates', ha='center', fontsize=8,
        color='#7f8c8d', style='italic')

# Arrows
arrow_down(ax, 2.6, 3.8, 2.8)
arrow_down(ax, 8.0, 3.8, 2.8)
arrow_down(ax, 13.4, 3.8, 2.8)

# "feeds into" labels
ax.text(0.5, 0.8, 'feeds into →', fontsize=9, color='#95a5a6', style='italic')

# ============================================================
# SAVE
# ============================================================
out = '/Users/aayanalwani/tb project/mce3r_stochastic/results/figures'
fig.savefig(f'{out}/Figure_Phase1A_DataAcquisition.png', dpi=300,
            bbox_inches='tight', facecolor='white')
fig.savefig(f'{out}/Figure_Phase1A_DataAcquisition.pdf', dpi=300,
            bbox_inches='tight', facecolor='white')
print("Saved PNG and PDF")
plt.close()
