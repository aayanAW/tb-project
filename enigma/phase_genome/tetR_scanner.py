"""
Genome-wide TetR-family operator classification.

Parses H37Rv GenBank to find all TetR-family regulators, extracts upstream
operator regions, scans for inverted repeats (palindromes), classifies
architectures, and predicts noise levels.

Key question: How common is Mce3R-like asymmetry among Mtb TetR operators?
"""

import sys
import os
import re
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from config.parameters import PARAMS
from Bio import SeqIO
from Bio.Seq import Seq


PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
GENBANK_PATH = os.path.join(PROJECT_ROOT, 'data', 'genomes', 'H37Rv.gb')
RESULTS_DIR = os.path.join(PROJECT_ROOT, 'results', 'phase_genome')


# ============================================================
# 1. PARSE TETR GENES FROM GENBANK
# ============================================================

TETR_PATTERNS = [
    re.compile(r'tetr', re.IGNORECASE),
    re.compile(r'TetR.family', re.IGNORECASE),
    re.compile(r'TetR-family', re.IGNORECASE),
    re.compile(r'TetR.type', re.IGNORECASE),
    re.compile(r'HTH.type.*transcriptional', re.IGNORECASE),
]

# Known TetR-family members in Mtb by locus_tag (supplement pattern matching)
KNOWN_TETR_LOCI = {
    'Rv1963c',  # Mce3R
    'Rv0078',   # EthR (ethionamide resistance)
    'Rv3066',   # TetR/AcrR family
    'Rv3574',   # kstR (cholesterol catabolism)
    'Rv3249c',  # TetR family
}


def parse_tetR_genes(genbank_path=GENBANK_PATH):
    """
    Find all TetR-family transcriptional regulators in H37Rv.

    Searches CDS feature qualifiers (product, note, function) for TetR patterns.
    Returns DataFrame with gene coordinates and annotations.
    """
    record = SeqIO.read(genbank_path, 'genbank')
    genome_seq = str(record.seq).upper()

    rows = []
    seen_loci = set()

    for feature in record.features:
        if feature.type != 'CDS':
            continue

        qualifiers = feature.qualifiers
        product = ' '.join(qualifiers.get('product', []))
        note = ' '.join(qualifiers.get('note', []))
        function = ' '.join(qualifiers.get('function', []))
        gene = ' '.join(qualifiers.get('gene', ['']))
        locus_tag = ' '.join(qualifiers.get('locus_tag', ['']))
        old_locus_tag = ' '.join(qualifiers.get('old_locus_tag', ['']))

        rv_id = locus_tag or old_locus_tag or gene

        # Check if TetR-family (by annotation or known locus list)
        combined = f"{product} {note} {function}"
        is_tetR = any(p.search(combined) for p in TETR_PATTERNS)
        is_known = rv_id in KNOWN_TETR_LOCI or locus_tag in KNOWN_TETR_LOCI

        if not is_tetR and not is_known:
            continue
        if rv_id in seen_loci:
            continue
        seen_loci.add(rv_id)

        start = int(feature.location.start)
        end = int(feature.location.end)
        strand = '+' if feature.location.strand == 1 else '-'

        rows.append({
            'rv_id': rv_id,
            'gene_name': gene,
            'product': product,
            'start': start,
            'end': end,
            'strand': strand,
        })

    df = pd.DataFrame(rows)
    print(f"  Found {len(df)} TetR-family genes in H37Rv")
    return df, genome_seq


# ============================================================
# 2. EXTRACT UPSTREAM / INTERGENIC REGIONS
# ============================================================

def _reverse_complement(seq):
    return str(Seq(seq).reverse_complement())


def _extract_circular(genome, start, length, genome_len):
    start = start % genome_len
    end = start + length
    if end <= genome_len:
        return genome[start:end]
    else:
        return genome[start:] + genome[:end - genome_len]


def extract_upstream_regions(tetR_df, genome_seq, upstream_len=300):
    """
    Extract upstream region for each TetR gene.

    For + strand genes: upstream is before start position.
    For - strand genes: upstream is after end position (reverse complemented).
    """
    genome_len = len(genome_seq)
    regions = []

    for _, row in tetR_df.iterrows():
        if row['strand'] == '+':
            up_start = row['start'] - upstream_len
            seq = _extract_circular(genome_seq, up_start, upstream_len, genome_len)
        else:
            up_start = row['end']
            seq = _extract_circular(genome_seq, up_start, upstream_len, genome_len)
            seq = _reverse_complement(seq)

        regions.append({
            'rv_id': row['rv_id'],
            'upstream_seq': seq,
            'upstream_len': len(seq),
        })

    return pd.DataFrame(regions)


# ============================================================
# 3. PALINDROME / INVERTED REPEAT SCANNER
# ============================================================

def _hamming_distance(s1, s2):
    return sum(c1 != c2 for c1, c2 in zip(s1, s2))


def scan_for_palindromes(sequence, min_arm=6, max_arm=15, max_spacer=60, max_mismatches=2):
    """
    Find inverted repeats (palindromic binding sites) in a DNA sequence.

    An inverted repeat has structure: ARM1---spacer---RC(ARM1)
    where RC = reverse complement.

    Returns list of palindrome hits sorted by quality score.
    """
    seq = sequence.upper()
    n = len(seq)
    hits = []

    for arm_len in range(max_arm, min_arm - 1, -1):
        for i in range(n - arm_len):
            arm1 = seq[i:i + arm_len]
            arm1_rc = _reverse_complement(arm1)

            # Search for arm1_rc downstream (with spacer)
            search_start = i + arm_len
            search_end = min(i + arm_len + max_spacer + arm_len, n)

            for j in range(search_start, search_end - arm_len + 1):
                arm2 = seq[j:j + arm_len]
                mismatches = _hamming_distance(arm1_rc, arm2)

                if mismatches <= max_mismatches:
                    spacer_len = j - (i + arm_len)
                    symmetry_score = 1.0 - (mismatches / arm_len)
                    quality = arm_len * symmetry_score

                    hits.append({
                        'arm1_start': i,
                        'arm1_end': i + arm_len,
                        'arm2_start': j,
                        'arm2_end': j + arm_len,
                        'arm1_seq': arm1,
                        'arm2_seq': arm2,
                        'spacer_len': spacer_len,
                        'arm_len': arm_len,
                        'mismatches': mismatches,
                        'symmetry_score': symmetry_score,
                        'quality': quality,
                    })

    # Deduplicate overlapping hits: keep best quality per region
    hits.sort(key=lambda h: h['quality'], reverse=True)
    filtered = []
    used_positions = set()
    for h in hits:
        pos_key = (h['arm1_start'] // 5, h['arm2_start'] // 5)  # bin by 5bp
        if pos_key not in used_positions:
            filtered.append(h)
            used_positions.add(pos_key)
        if len(filtered) >= 5:  # top 5 per sequence
            break

    return filtered


# ============================================================
# 4. OPERATOR ARCHITECTURE CLASSIFICATION
# ============================================================

def classify_architecture(palindromes, upstream_seq):
    """
    Classify operator architecture based on palindrome scan results.

    Categories:
    - palindromic: high symmetry, both arms similar affinity
    - asymmetric: inverted repeat but arms differ significantly
    - single_site: only one detectable binding motif
    - tandem: same-orientation repeats (not inverted)
    - unresolved: no clear pattern detected
    """
    if not palindromes:
        return {
            'architecture': 'unresolved',
            'symmetry_score': 0.0,
            'best_arm_len': 0,
            'spacer_len': 0,
            'Kd_ratio_est': 1.0,
            'arm1_seq': '',
            'arm2_seq': '',
        }

    best = palindromes[0]

    # Estimate Kd ratio from arm similarity
    # Each mismatch in the recognition helix weakens binding ~3-5 fold
    # But the heuristic must distinguish truly palindromic (0 mismatches) from asymmetric
    kd_ratio_est = 3.0 ** best['mismatches']

    # Look for PERFECT palindromes (0 mismatches) — the real signal
    perfect = [p for p in palindromes if p['mismatches'] == 0 and p['arm_len'] >= 8]
    near_perfect = [p for p in palindromes if p['mismatches'] <= 1 and p['arm_len'] >= 8]
    imperfect = [p for p in palindromes if p['mismatches'] == 2 and p['arm_len'] >= 10]

    if perfect:
        arch = 'palindromic'
        best = perfect[0]
        kd_ratio_est = 1.0
    elif near_perfect:
        arch = 'palindromic'
        best = near_perfect[0]
        kd_ratio_est = 2.0
    elif imperfect and best['arm_len'] >= 12:
        # Long arms with 2 mismatches: genuine but asymmetric binding
        arch = 'asymmetric'
        kd_ratio_est = 3.0 ** best['mismatches']
    else:
        # Only short or heavily mismatched hits — background noise
        arch = 'single_site'
        kd_ratio_est = 1.0

    return {
        'architecture': arch,
        'symmetry_score': best['symmetry_score'],
        'best_arm_len': best['arm_len'],
        'spacer_len': best['spacer_len'],
        'Kd_ratio_est': kd_ratio_est,
        'arm1_seq': best['arm1_seq'],
        'arm2_seq': best['arm2_seq'],
        'mismatches': best['mismatches'],
    }


# ============================================================
# 5. NOISE PREDICTION
# ============================================================

def predict_noise(Kd_ratio, architecture, n_cells=300, seed=42):
    """
    Predict expression noise (CV) for a given operator architecture.

    Uses the Gillespie engine with estimated Kd values based on the
    Mce3R calibration (geometric mean ~10.84 nM).
    """
    from phase_evo.genetic_algorithm import run_population_fast
    from phase2_simulation.operator_model import OperatorModel

    # Use Mce3R geometric mean as baseline, split by ratio
    Kd_geo = PARAMS['Kd_symmetric']  # 10.84 nM
    Kd_s = Kd_geo / np.sqrt(Kd_ratio)
    Kd_w = Kd_geo * np.sqrt(Kd_ratio)

    if architecture == 'palindromic':
        block_s, block_w = 0.65, 0.65
    elif architecture == 'asymmetric':
        block_s, block_w = 0.85, 0.50
    elif architecture == 'single_site':
        Kd_w = 1e12  # disable weak site
        block_s, block_w = 0.85, 0.50
    else:  # unresolved — use symmetric as default
        block_s, block_w = 0.65, 0.65
        Kd_ratio = 1.0
        Kd_s = Kd_geo
        Kd_w = Kd_geo

    model = OperatorModel(
        Kd_strong=Kd_s, Kd_weak=Kd_w,
        block_strong=block_s, block_weak=block_w,
    )
    model_arrays = model.get_numba_arrays()
    proteins = run_population_fast(model_arrays, n_cells, seed)

    pf = proteins.astype(float)
    mean_prot = np.mean(pf)
    std_prot = np.std(pf, ddof=1) if len(pf) > 1 else 0.0
    cv = std_prot / mean_prot if mean_prot > 0 else 0.0
    fano = np.var(pf, ddof=1) / mean_prot if mean_prot > 0 else 0.0
    persister_frac = float(np.sum(pf < 112.0)) / len(pf)

    return {
        'Kd_strong_est': Kd_s,
        'Kd_weak_est': Kd_w,
        'mean_protein': mean_prot,
        'cv': cv,
        'fano': fano,
        'persister_fraction': persister_frac,
    }


# ============================================================
# 6. FULL GENOME SCAN
# ============================================================

def run_genome_scan(n_cells_per_tf=300, seed=42):
    """
    Complete genome-wide TetR operator classification pipeline.

    Returns three DataFrames: catalog, architectures, noise predictions.
    """
    os.makedirs(RESULTS_DIR, exist_ok=True)

    # Step 1: Parse TetR genes
    print("  Step 1: Parsing TetR genes from H37Rv.gb...")
    tetR_df, genome_seq = parse_tetR_genes()

    if len(tetR_df) == 0:
        print("  WARNING: No TetR genes found. Check GenBank annotations.")
        return tetR_df, pd.DataFrame(), pd.DataFrame()

    # Step 2: Extract upstream regions
    print("  Step 2: Extracting upstream regions...")
    upstream_df = extract_upstream_regions(tetR_df, genome_seq)
    catalog = tetR_df.merge(upstream_df, on='rv_id', how='left')

    # Step 3: Scan for palindromes and classify
    print("  Step 3: Scanning for palindromes and classifying...")
    arch_rows = []
    for _, row in catalog.iterrows():
        seq = row.get('upstream_seq', '')
        if not seq or len(seq) < 20:
            arch_rows.append({
                'rv_id': row['rv_id'],
                'architecture': 'unresolved',
                'symmetry_score': 0.0,
                'Kd_ratio_est': 1.0,
                'best_arm_len': 0,
                'spacer_len': 0,
                'arm1_seq': '',
                'arm2_seq': '',
                'mismatches': 0,
            })
            continue

        palindromes = scan_for_palindromes(seq)
        classification = classify_architecture(palindromes, seq)
        classification['rv_id'] = row['rv_id']
        arch_rows.append(classification)

    arch_df = pd.DataFrame(arch_rows)

    # Step 4: Predict noise levels
    print(f"  Step 4: Predicting noise for {len(arch_df)} TFs...")
    noise_rows = []
    for i, row in arch_df.iterrows():
        pred = predict_noise(
            row['Kd_ratio_est'], row['architecture'],
            n_cells=n_cells_per_tf, seed=seed + i * 100,
        )
        pred['rv_id'] = row['rv_id']
        pred['architecture'] = row['architecture']
        noise_rows.append(pred)

    noise_df = pd.DataFrame(noise_rows)

    return catalog, arch_df, noise_df
