"""
mce3r_biology.py — Single source of truth for Mce3R biology.

Every script imports the regulon, operator regions, and negative-control sets from
here so the repo cannot contradict itself (the previous version disagreed about the
regulon coordinates across files). All facts are cited to primary sources and were
verified against Mycobrowser / UniProt / the cited papers.

Primary sources
---------------
- Santangelo MP et al. 2008, BMC Microbiology 8:38.
    Mce3R represses the mce3 operon and autoregulates; does NOT regulate mce1/mce2/mce4.
- de la Paz Santangelo M et al. 2009, Microbiology 155:2245 (doi:10.1099/mic.0.027086-0).
    Mce3R also represses the divergent Rv1933c-Rv1935c and Rv1936-Rv1941 transcriptional
    units (lipid metabolism / redox). MEME motif found 6x in the mce3R-yrbE3A IGR and 3x
    in the Rv1935c-Rv1936 IGR; footprinting confirmed protected operator regions.
- Panagoda, Balazsi & Sampson 2024, ACS Chem Biol 19:2580-2592
    (doi:10.1021/acschembio.4c00687; PDB 9B7Y).
    Mce3R is an unusual fused double-TFR. The yrbE3A operator = TWO ~25 bp NONPALINDROMIC
    sites separated by ~53 bp; higher-affinity (downstream) site ~100 bp upstream of the
    yrbE3A start codon; Kd 2.4 +/- 0.7 nM.

Genome: M. tuberculosis H37Rv, NC_000962.3 (AL123456.3), 4,411,532 bp.
"""

from __future__ import annotations

# ── Regulon (the ground-truth positive set) ─────────────────────────────────────

# Mce3R autoregulation.
AUTOREGULATION = {"Rv1963c"}

# mce3 structural operon (yrbE3A, yrbE3B, mce3A-F), all + strand, divergent from mce3R.
MCE3_OPERON = {
    "Rv1964": "yrbE3A",
    "Rv1965": "yrbE3B",
    "Rv1966": "mce3A",
    "Rv1967": "mce3B",
    "Rv1968": "mce3C",
    "Rv1969": "mce3D",
    "Rv1970": "mce3E",
    "Rv1971": "mce3F",
}

# Divergent lipid/redox cluster repressed by Mce3R (Santangelo 2009).
# Two adjacent, divergently-transcribed transcriptional units.
DIVERGENT_CLUSTER = {
    "Rv1933c": "fadE18",
    "Rv1934c": "fadE17",
    "Rv1935c": "echA13",
    "Rv1936": "",
    "Rv1937": "",
    "Rv1938": "ephB",
    "Rv1939": "",
    "Rv1940": "ribA1",
    "Rv1941": "",
}

# Full known Mce3R regulon = autoregulation + mce3 operon + divergent cluster.
KNOWN_REGULON = set(AUTOREGULATION) | set(MCE3_OPERON) | set(DIVERGENT_CLUSTER)

GENE_NAMES = {"Rv1963c": "mce3R", **MCE3_OPERON, **DIVERGENT_CLUSTER}


# ── Negative controls (must NOT be Mce3R targets) ───────────────────────────────
# Santangelo 2008 showed no effect on mce1/mce2/mce4. If the motif scores these like
# real targets, it is non-specific. Used as a pre-registered specificity control.

MCE1_OPERON = {f"Rv{n:04d}" for n in range(169, 179)}  # Rv0169-Rv0178
MCE2_OPERON = {f"Rv{n:04d}" for n in range(586, 595)}  # Rv0586-Rv0594
MCE4_OPERON = {f"Rv{n:04d}c" for n in range(3499, 3504)}  # Rv3499c-Rv3503c
NEGATIVE_CONTROL_GENES = MCE1_OPERON | MCE2_OPERON | MCE4_OPERON


# ── Mapped operator regions (divergent intergenic regions) ──────────────────────
# Identified by the IGR ids produced by extract_promoters.extract_divergent_igrs:
#   IGR_<upstream_minus_gene>_<downstream_plus_gene>
PRIMARY_OPERATOR_IGR = "IGR_Rv1963c_Rv1964"  # mce3R-yrbE3A, primary operator
SECOND_OPERATOR_IGR = "IGR_Rv1935c_Rv1936"  # second footprinted region (2009)
OPERATOR_IGRS = {PRIMARY_OPERATOR_IGR, SECOND_OPERATOR_IGR}

# Sequence identifiers that count as "the known operator" for site-level recovery.
# Includes the IGR ids and the flanking-gene promoters that overlap the operator.
KNOWN_OPERATOR_SEQUENCES = OPERATOR_IGRS | {"Rv1963c", "Rv1964", "Rv1935c", "Rv1936"}

# yrbE3A (Rv1964) translation start on + strand (NC_000962.3). The operator's two
# ~25 bp sites lie within ~[40, 170] bp upstream of this start. Used only to define
# the genomic WINDOW from which the knowledge-based PWM is built from REAL sequence;
# no operator sequence is hard-coded.
YRBE3A_START = 2207700
OPERATOR_WINDOW_UPSTREAM_BP = (40, 170)  # (closest, farthest) bp upstream of start
OPERATOR_SITE_WIDTH = 25  # bp per half-site (2024 paper)
OPERATOR_SITE_SPACER = 53  # bp between the two sites (2024 paper)

# Architecture facts (descriptive only — never used to inflate a "match" verdict).
OPERATOR_IS_PALINDROMIC = False  # explicitly nonpalindromic (2024 paper)


# ── Held-out vs training designation (anti-circularity) ─────────────────────────
# The knowledge-based PWM (build_operator_model.py) is trained ONLY on the
# mce3R-yrbE3A operator window (PRIMARY_OPERATOR_IGR) + ortholog upstream windows.
# The Rv1935c-Rv1936 operator region (SECOND_OPERATOR_IGR) is NEVER used in PWM
# construction, so its recovery is a genuine HELD-OUT generalization test, not the
# tautological re-recovery of training data flagged in the audit (C2/C4).
TRAINING_OPERATOR_IGR = PRIMARY_OPERATOR_IGR
HELDOUT_OPERATOR_IGR = SECOND_OPERATOR_IGR


# ── Downstream co-transcribed regulon genes that carry NO operator ──────────────
# Mce3R binds the two divergent operator IGRs; the downstream genes of each operon
# are repressed because they are co-transcribed, NOT because each has its own
# operator site (Santangelo 2009). A genuine binding-site scan must therefore NOT
# light these up — they are the specificity check the leaky "regulon AUPRC" missed.
# These are the regulon genes EXCLUDING the two operator-flanking promoters.
OPERATOR_FLANKING_GENES = {"Rv1963c", "Rv1964", "Rv1935c", "Rv1936"}
DOWNSTREAM_REGULON_NO_OPERATOR = KNOWN_REGULON - OPERATOR_FLANKING_GENES


# ── Mycobacterial lineage grouping for conservation (anti-pseudo-replication) ────
# The conservation gate previously counted 4 near-clonal M. tuberculosis-complex
# (MTBC) members as 4 independent species (>99.9% identical → pseudo-replication,
# audit C7). Collapse each MTBC member to ONE effective lineage; the M. marinum
# clade (marinum/ulcerans/liflandii) and the M. smegmatis outgroup are independent.
# Keyed by the species_key prefix used in ortholog FASTA ids
# (e.g. "mycobacterium_tuberculosis").
ORTHOLOG_LINEAGE = {
    "mycobacterium_tuberculosis": "MTBC",
    "mycobacterium_bovis": "MTBC",
    "mycobacterium_africanum": "MTBC",
    "mycobacterium_canettii": "MTBC",
    "mycobacterium_marinum": "M_marinum_clade",
    "mycobacterium_ulcerans": "M_marinum_clade",
    "mycobacterium_liflandii": "M_marinum_clade",
    "mycobacterium_smegmatis": "M_smegmatis",
}
# Effective independent lineages available in the frozen ortholog set = 3
# (MTBC, M. marinum clade, M. smegmatis). The conservation gate is judged at the
# lineage level, not the species level.
CONSERVATION_MIN_LINEAGES = 2  # operator conserved in >= 2 independent lineages


def ortholog_lineage(seq_id: str) -> str:
    """Map an ortholog FASTA id to its effective lineage (defaults to the id stem)."""
    key = "_".join(str(seq_id).split("_")[:2]).lower()
    return ORTHOLOG_LINEAGE.get(key, key)


# ── Independently-published operator motif sites (EXTERNAL reference) ────────────
# These site sequences were derived by OTHER laboratories using DNase I footprinting
# and independent MEME analysis — they are NOT outputs of this pipeline. Matching the
# knowledge-based PWM against a PWM built from these sites is therefore a genuine
# external-corroboration test (G1), not two MEME runs over shared DNA agreeing
# (the circularity flagged in audit C1).
#
# Source: de la Paz Santangelo M et al. 2009, Microbiology 155:2245
#   (doi:10.1099/mic.0.027086-0) — operator half-site sequences in the mce3R-yrbE3A
#   and Rv1935c-Rv1936 intergenic regions (Fig. 6 / Table; consensus TANNNCAAAATACAT).
# Cross-referenced against the Eram Kabir 2021 lab-rotation footprint summary
# (mce3r_stochastic/reference_docs/Eram_ResearchSummaryFinalUPD.pdf).
SANTANGELO_2009_OPERATOR_SITES = (
    "GACTAACAAAATACAT",
    "TATCAGTAATAGACAT",
    "TACTAGCAAGATACAT",
    "TATTGGCTATGGACAT",
)
SANTANGELO_2009_CITATION = (
    "Santangelo et al. 2009, Microbiology 155:2245 (doi:10.1099/mic.0.027086-0)"
)


def classify_locus(locus_tag: str) -> str:
    """
    Classify a sequence/locus for validation.

    Returns one of:
        'operator'         — a mapped operator region or its flanking promoter
        'regulon'          — a known Mce3R regulon gene (positive)
        'negative_control' — mce1/mce2/mce4 (must behave like background)
        'other'            — everything else (the genome-wide background)
    """
    tag = str(locus_tag).strip()
    if tag in KNOWN_OPERATOR_SEQUENCES:
        return "operator"
    if tag in KNOWN_REGULON:
        return "regulon"
    if tag in NEGATIVE_CONTROL_GENES:
        return "negative_control"
    return "other"


def is_positive(locus_tag: str) -> bool:
    """True for ground-truth positives (operator regions + regulon genes)."""
    return classify_locus(locus_tag) in ("operator", "regulon")
