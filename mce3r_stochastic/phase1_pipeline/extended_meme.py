"""
Extended MEME analysis replicating Eram Kabir's 2021 protocol:

1. MEME with maxw up to 120 bp (not just 30)
2. ANR mode (any number of repetitions)
3. Multiple motifs discovered per run
4. Run on individual regions AND on Regions 1+2 combined
5. Try to recover Eram's 49 bp motif de novo

Addresses PI requests:
  #2: Don't treat Panagoda 2024 as the only gold standard
  #3: Allow MEME to find longer sites (up to 120 bp)
  #4: Repeat Eram's analysis and recover his long motif
"""

import sys
import os
import re
import json
import subprocess
import shutil
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from config.parameters import PARAMS
from Bio.Seq import Seq

PROJECT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ERAM_DIR = os.path.join(PROJECT, 'results', 'eram_validation')
OUT_DIR = os.path.join(PROJECT, 'results', 'extended_meme')
os.makedirs(OUT_DIR, exist_ok=True)


# Eram's published 49 bp motif consensus (from his Figure 6)
ERAM_MOTIF1 = 'ACTAGCAAGATACATCATAGCCAATATATGCCAGTTTGCATTGCTATTT'


def run_meme(input_fasta, output_dir, minw=8, maxw=120, nmotifs=10, mode='anr'):
    """
    Run MEME with Eram's protocol: ANR mode, wide range, multiple motifs.

    Parameters
    ----------
    input_fasta : str
        Path to input FASTA
    output_dir : str
        Where MEME should write results (will be deleted if exists)
    minw, maxw : int
        Motif width bounds
    nmotifs : int
        How many motifs to find
    mode : str
        'anr' (any number of repetitions), 'zoops', or 'oops'
    """
    if os.path.exists(output_dir):
        shutil.rmtree(output_dir)

    cmd = [
        'meme', input_fasta,
        '-dna',
        '-mod', mode,
        '-revcomp',
        '-minw', str(minw),
        '-maxw', str(maxw),
        '-nmotifs', str(nmotifs),
        '-oc', output_dir,
    ]

    print(f"  Running: meme {os.path.basename(input_fasta)} "
          f"-mod {mode} -minw {minw} -maxw {maxw} -nmotifs {nmotifs}")

    result = subprocess.run(cmd, capture_output=True, text=True, timeout=1800)

    meme_txt = os.path.join(output_dir, 'meme.txt')
    if not os.path.exists(meme_txt):
        print(f"  [ERROR] MEME produced no output: {result.stderr[:300]}")
        return None

    print(f"  [OK] MEME output in {output_dir}")
    return meme_txt


def parse_meme_motifs(meme_txt_path):
    """Parse discovered motifs from meme.txt — width, E-value, sites, consensus."""
    if not os.path.exists(meme_txt_path):
        return []

    with open(meme_txt_path) as f:
        content = f.read()

    motifs = []
    # Find motif blocks
    motif_pattern = re.compile(
        r'MOTIF\s+(\S+)\s+MEME-(\d+).*?width\s*=\s*(\d+).*?sites\s*=\s*(\d+).*?E-value\s*=\s*(\S+)',
        re.DOTALL
    )

    for match in motif_pattern.finditer(content):
        motif_id = match.group(1)
        meme_num = int(match.group(2))
        width = int(match.group(3))
        sites = int(match.group(4))
        evalue = match.group(5)

        # Extract consensus
        consensus_match = re.search(
            rf'Motif\s+{re.escape(motif_id)}\s+MEME-{meme_num}\s+regular expression\s*\n-+\s*\n(\S+)',
            content
        )
        consensus = consensus_match.group(1) if consensus_match else '?'

        # Extract individual sites
        sites_block = re.search(
            rf'Motif\s+{re.escape(motif_id)}\s+MEME-{meme_num}\s+sites sorted by position p-value.*?\n-+\n(.+?)\n-+',
            content, re.DOTALL
        )
        site_seqs = []
        if sites_block:
            for line in sites_block.group(1).strip().split('\n'):
                parts = line.split()
                if len(parts) >= 6:
                    # Format: seq_name strand start p-value [flanking] SITE [flanking]
                    # The site is typically 4th or 5th column depending on format
                    # Look for uppercase DNA string
                    for tok in parts:
                        if re.match(r'^[ACGTN]+$', tok) and len(tok) >= minw_local(width):
                            site_seqs.append(tok)
                            break

        motifs.append({
            'motif_id': motif_id,
            'meme_number': meme_num,
            'width': width,
            'n_sites': sites,
            'e_value': evalue,
            'consensus_regex': consensus,
            'site_sequences': site_seqs,
        })

    return motifs


def minw_local(width):
    """Minimum accepted site length (use 80% of width)."""
    return int(width * 0.8)


def find_in_sequence(sequence, motif_consensus, max_mm=3):
    """Find motif occurrences in a sequence with mismatch tolerance."""
    seq = sequence.upper()
    # Convert IUPAC regex consensus to plain nucleotides (approximate)
    clean = re.sub(r'\[([ACGT])[ACGT]*\]', r'\1', motif_consensus.upper())
    clean = re.sub(r'[^ACGT]', 'N', clean)

    if len(clean) < 8:
        return []

    rc = str(Seq(clean).reverse_complement())
    hits = []
    for m, strand in [(clean, '+'), (rc, '-')]:
        for i in range(len(seq) - len(m) + 1):
            mm = sum(a != b for a, b in zip(seq[i:i+len(m)], m) if b != 'N')
            if mm <= max_mm:
                hits.append({'pos': i, 'strand': strand, 'mm': mm, 'matched': seq[i:i+len(m)]})
    hits.sort(key=lambda h: h['mm'])
    return hits


def compare_to_eram(motifs, eram_motif=ERAM_MOTIF1):
    """Check if any discovered motif overlaps with Eram's 49 bp consensus."""
    eram_len = len(eram_motif)
    results = []

    for m in motifs:
        # Approach: check if one of MEME's discovered sites appears inside Eram's motif
        # OR if Eram's motif contains MEME's consensus
        for site in m.get('site_sequences', []):
            if len(site) < 8:
                continue
            # Check forward
            mm_fwd = 999
            if len(site) <= len(eram_motif):
                for i in range(len(eram_motif) - len(site) + 1):
                    mm = sum(a != b for a, b in zip(site, eram_motif[i:i+len(site)]))
                    mm_fwd = min(mm_fwd, mm)

            # Check reverse complement
            rc = str(Seq(site).reverse_complement())
            mm_rc = 999
            if len(rc) <= len(eram_motif):
                for i in range(len(eram_motif) - len(rc) + 1):
                    mm = sum(a != b for a, b in zip(rc, eram_motif[i:i+len(rc)]))
                    mm_rc = min(mm_rc, mm)

            best_mm = min(mm_fwd, mm_rc)
            identity = 1.0 - (best_mm / len(site))

            results.append({
                'motif_id': m['motif_id'],
                'motif_width': m['width'],
                'motif_evalue': m['e_value'],
                'site_seq': site,
                'site_len': len(site),
                'best_mm_vs_eram': best_mm,
                'identity_vs_eram': identity,
                'direction': 'forward' if mm_fwd <= mm_rc else 'reverse',
                'recovers_eram': identity >= 0.8,
            })

    return results


def run_all():
    """Run extended MEME on all target inputs."""
    print("=" * 70)
    print("EXTENDED MEME ANALYSIS — Eram's protocol replication")
    print("=" * 70)

    # Use the already-prepared FASTA files from eram_validation
    inputs = [
        ('region_1_only',    os.path.join(ERAM_DIR, 'region_1.fasta')),
        ('region_2_only',    os.path.join(ERAM_DIR, 'region_2.fasta')),
        ('regions_1_2_combined', os.path.join(ERAM_DIR, 'regions_1_2_combined.fasta')),
        ('all_5_regions',    os.path.join(ERAM_DIR, 'all_regions_combined.fasta')),
    ]

    all_results = {}

    for name, path in inputs:
        if not os.path.exists(path):
            print(f"  [SKIP] {name}: input file missing ({path})")
            continue

        print(f"\n--- {name} ---")
        out_dir = os.path.join(OUT_DIR, name)
        meme_txt = run_meme(path, out_dir, minw=8, maxw=120, nmotifs=10, mode='anr')

        if meme_txt is None:
            all_results[name] = {'status': 'failed', 'motifs': []}
            continue

        motifs = parse_meme_motifs(meme_txt)
        print(f"  Found {len(motifs)} motifs")

        for m in motifs[:5]:
            print(f"    Motif {m['motif_id']}: width={m['width']} E={m['e_value']} "
                  f"n_sites={m['n_sites']}")
            print(f"      consensus: {m['consensus_regex'][:60]}")

        # Check for Eram's motif recovery
        eram_check = compare_to_eram(motifs)
        recovered = [r for r in eram_check if r['recovers_eram']]

        if recovered:
            best = max(recovered, key=lambda r: r['identity_vs_eram'])
            print(f"  [ERAM MOTIF RECOVERED] motif {best['motif_id']} at "
                  f"{best['identity_vs_eram']:.0%} identity")
        else:
            # At least report closest match
            if eram_check:
                closest = max(eram_check, key=lambda r: r['identity_vs_eram'])
                print(f"  [closest to Eram: {closest['identity_vs_eram']:.0%} — "
                      f"motif {closest['motif_id']}]")

        all_results[name] = {
            'status': 'ok',
            'n_motifs': len(motifs),
            'motifs': motifs,
            'eram_comparison': eram_check,
            'eram_recovered': len(recovered) > 0,
        }

    # Save summary
    summary = {
        'description': 'Extended MEME analysis replicating Eram Kabir 2021 protocol',
        'parameters': {'minw': 8, 'maxw': 120, 'nmotifs': 10, 'mode': 'anr', 'revcomp': True},
        'inputs': [name for name, _ in inputs],
        'eram_recovered_in': [name for name, r in all_results.items() if r.get('eram_recovered')],
        'results': {name: {
            'status': r['status'],
            'n_motifs': r.get('n_motifs', 0),
            'top_motifs': [{
                'id': m['motif_id'],
                'width': m['width'],
                'e_value': m['e_value'],
                'n_sites': m['n_sites'],
            } for m in r.get('motifs', [])[:5]],
            'eram_recovered': r.get('eram_recovered', False),
        } for name, r in all_results.items()},
    }

    summary_path = os.path.join(OUT_DIR, 'extended_meme_summary.json')
    with open(summary_path, 'w') as f:
        json.dump(summary, f, indent=2, default=str)

    # Save flattened motif table for downstream use
    motif_rows = []
    for input_name, result in all_results.items():
        for m in result.get('motifs', []):
            motif_rows.append({
                'input': input_name,
                'motif_id': m['motif_id'],
                'width': m['width'],
                'e_value': m['e_value'],
                'n_sites': m['n_sites'],
                'consensus': m['consensus_regex'],
                'first_site': m['site_sequences'][0] if m['site_sequences'] else '',
            })
    pd.DataFrame(motif_rows).to_csv(os.path.join(OUT_DIR, 'discovered_motifs.csv'), index=False)

    # Save Eram comparison table
    eram_rows = []
    for input_name, result in all_results.items():
        for r in result.get('eram_comparison', []):
            r_copy = dict(r)
            r_copy['input'] = input_name
            eram_rows.append(r_copy)
    pd.DataFrame(eram_rows).to_csv(os.path.join(OUT_DIR, 'eram_comparison.csv'), index=False)

    print(f"\n{'=' * 70}")
    print(f"Extended MEME complete. Results in {OUT_DIR}")
    print(f"Eram motif recovered in: {summary['eram_recovered_in']}")
    return summary, all_results


if __name__ == '__main__':
    run_all()
