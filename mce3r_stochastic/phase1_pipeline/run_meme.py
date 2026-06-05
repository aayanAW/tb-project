"""
phase1_pipeline/run_meme.py — Run MEME motif discovery or generate mock output.

- Checks shutil.which('meme') for MEME suite availability.
- If present: run MEME with -dna -mod zoops -revcomp -minw 20 -maxw 30 -nmotifs 3
  (NO -pal flag). Generates background model with fasta-get-markov.
- If missing: generate_mock_meme_output() with a 25-bp PWM derived from the strong
  binding site in the operator sequence.
"""

import sys
import os
import shutil
import subprocess
import re

sys.path.insert(0, '/Users/aayanalwani/tb project/mce3r_stochastic')

# Ensure conda env bin is on PATH so shutil.which() can find meme/fimo
_conda_bin = os.path.join(sys.prefix, 'bin')
if _conda_bin not in os.environ.get('PATH', ''):
    os.environ['PATH'] = _conda_bin + os.pathsep + os.environ.get('PATH', '')

from config.parameters import PARAMS

PROJECT_ROOT = '/Users/aayanalwani/tb project/mce3r_stochastic'
SEQ_DIR = os.path.join(PROJECT_ROOT, 'data', 'sequences')
RESULTS_DIR = os.path.join(PROJECT_ROOT, 'results', 'phase1')


def ensure_dirs():
    os.makedirs(RESULTS_DIR, exist_ok=True)


def generate_background_model(input_fasta, bg_path):
    """Generate a 0-order Markov background model using fasta-get-markov."""
    if shutil.which('fasta-get-markov'):
        cmd = ['fasta-get-markov', '-dna', input_fasta, bg_path]
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode == 0:
            print(f"  Background model generated: {bg_path}")
            return bg_path
        else:
            print(f"  WARNING: fasta-get-markov failed: {result.stderr}")
    # Fallback: generate from H37Rv GC content
    gc = PARAMS['H37Rv_gc_content']
    at = 1.0 - gc
    with open(bg_path, 'w') as f:
        f.write("# 0-order Markov background model for H37Rv\n")
        f.write(f"A {at/2:.6f}\n")
        f.write(f"C {gc/2:.6f}\n")
        f.write(f"G {gc/2:.6f}\n")
        f.write(f"T {at/2:.6f}\n")
    print(f"  Background model generated (from GC content): {bg_path}")
    return bg_path


def run_real_meme(input_fasta, output_dir, bg_path):
    """Run MEME motif discovery with proper parameters."""
    import shutil as _shutil
    if os.path.exists(output_dir):
        _shutil.rmtree(output_dir)
    os.makedirs(output_dir, exist_ok=True)
    cmd = [
        'meme', input_fasta,
        '-dna',
        '-mod', PARAMS['meme_mod'],
        '-revcomp',
        '-minw', str(PARAMS['meme_minw']),
        '-maxw', str(PARAMS['meme_maxw']),
        '-nmotifs', str(PARAMS['meme_nmotifs']),
        '-bfile', bg_path,
        '-oc', output_dir,
    ]
    print(f"  Running MEME: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
    meme_txt = os.path.join(output_dir, 'meme.txt')
    if result.returncode != 0:
        print(f"  MEME stderr: {result.stderr[-1000:]}")
        # MEME may return non-zero due to meme_xml_to_html Ghostscript issues
        # Check if core output was still produced
        if os.path.exists(meme_txt) and os.path.getsize(meme_txt) > 100:
            print(f"  WARNING: MEME returncode={result.returncode} but meme.txt was produced — continuing")
        else:
            raise RuntimeError(f"MEME failed with return code {result.returncode}")
    print(f"  MEME completed, output in {output_dir}")
    return output_dir


def generate_mock_meme_output():
    """
    Generate mock MEME output with a 25-bp PWM derived from the strong binding site
    portion of the operator sequence.

    The strong binding site is approximately the first 25 bp of the operator.
    We create a PWM with high probability at the consensus positions and
    low background elsewhere.
    """
    output_dir = os.path.join(RESULTS_DIR, 'meme_output')
    os.makedirs(output_dir, exist_ok=True)

    # Use the first 25 bp of the operator as the strong binding site motif
    strong_site = PARAMS['operator_sequence'][:25]  # 'GCCCCGCGCTATAGGATACTAGCAA'
    motif_width = len(strong_site)  # 25

    # Build a PWM: high probability (0.85) at consensus, distribute rest
    pwm_lines = []
    for base in strong_site.upper():
        row = {'A': 0.05, 'C': 0.05, 'G': 0.05, 'T': 0.05}
        row[base] = 0.85
        pwm_lines.append(row)

    # Also create a second motif from the weak site region (middle of operator)
    weak_site = PARAMS['operator_sequence'][49:74]  # 25 bp from weak site region
    pwm_lines_2 = []
    for base in weak_site.upper():
        # Probability at consensus = 0.80, remaining 0.20 split among other 3 bases
        bg = 0.20 / 3.0
        row = {'A': bg, 'C': bg, 'G': bg, 'T': bg}
        row[base] = 0.80
        pwm_lines_2.append(row)

    # Write MEME minimal format output
    meme_txt_path = os.path.join(output_dir, 'meme.txt')
    with open(meme_txt_path, 'w') as f:
        f.write("MEME version 5.5.5\n\n")
        f.write("ALPHABET= ACGT\n\n")
        f.write("strands: + -\n\n")

        # Background frequencies
        gc = PARAMS['H37Rv_gc_content']
        at = 1.0 - gc
        f.write("Background letter frequencies (from H37Rv):\n")
        f.write(f"A {at/2:.6f} C {gc/2:.6f} G {gc/2:.6f} T {at/2:.6f}\n\n")

        # Motif 1 - strong site
        f.write(f"MOTIF 1-MEME_mock_strong_site\n")
        f.write(f"letter-probability matrix: alength= 4 w= {motif_width} "
                f"nsites= 3 E= 1.0e-010\n")
        for row in pwm_lines:
            f.write(f" {row['A']:.6f}  {row['C']:.6f}  {row['G']:.6f}  {row['T']:.6f}\n")
        f.write("\n")

        # Motif 2 - weak site
        f.write(f"MOTIF 2-MEME_mock_weak_site\n")
        f.write(f"letter-probability matrix: alength= 4 w= {len(weak_site)} "
                f"nsites= 3 E= 1.0e-005\n")
        for row in pwm_lines_2:
            f.write(f" {row['A']:.6f}  {row['C']:.6f}  {row['G']:.6f}  {row['T']:.6f}\n")
        f.write("\n")

    # Write the motif in MEME minimal format for FIMO input
    meme_motif_path = os.path.join(output_dir, 'meme.html')
    with open(meme_motif_path, 'w') as f:
        f.write("<html><body><p>Mock MEME output (tools not installed)</p></body></html>\n")

    print(f"  Mock MEME output generated: {output_dir}")
    print(f"  Motif 1: {motif_width}-bp PWM from strong binding site")
    print(f"  Motif 2: {len(weak_site)}-bp PWM from weak binding site region")

    return output_dir


def run_meme_step():
    """
    Main MEME step. Returns dict with output paths and metadata.
    """
    ensure_dirs()

    input_fasta = os.path.join(SEQ_DIR, 'meme_input_orthologs.fasta')
    if not os.path.exists(input_fasta):
        raise FileNotFoundError(f"Input FASTA not found: {input_fasta}")

    bg_path = os.path.join(RESULTS_DIR, 'background.model')
    generate_background_model(input_fasta, bg_path)

    mock_used = False
    if shutil.which('meme'):
        print("  MEME suite found, running real MEME...")
        output_dir = run_real_meme(input_fasta,
                                    os.path.join(RESULTS_DIR, 'meme_output'),
                                    bg_path)
    else:
        print("  MEME suite NOT found, generating mock output...")
        output_dir = generate_mock_meme_output()
        mock_used = True

    meme_txt = os.path.join(output_dir, 'meme.txt')

    return {
        'meme_output_dir': output_dir,
        'meme_txt': meme_txt,
        'background_model': bg_path,
        'mock_used': mock_used,
    }


def parse_meme_motif_widths(meme_txt_path):
    """Parse motif widths from a MEME text output file."""
    widths = []
    with open(meme_txt_path) as f:
        for line in f:
            m = re.search(r'w=\s*(\d+)', line)
            if m:
                widths.append(int(m.group(1)))
    return widths


if __name__ == "__main__":
    print("=" * 60)
    print("run_meme.py — Self-test")
    print("=" * 60)

    n_sanity_pass = 0
    n_sanity_fail = 0
    n_science_expected = 0
    n_science_unexpected = 0

    try:
        result = run_meme_step()
    except Exception as e:
        import traceback
        traceback.print_exc()
        print(f"SANITY FAIL: Exception during MEME step: {e}")
        sys.exit(1)

    # --- SANITY CHECKS (hard-fail) ---

    # 1. MEME output directory exists
    if os.path.isdir(result['meme_output_dir']):
        print(f"SANITY PASS: MEME output directory exists: {result['meme_output_dir']}")
        n_sanity_pass += 1
    else:
        print(f"SANITY FAIL: MEME output directory missing")
        n_sanity_fail += 1

    # 2. meme.txt exists and is non-empty
    meme_txt = result['meme_txt']
    if os.path.exists(meme_txt) and os.path.getsize(meme_txt) > 100:
        print(f"SANITY PASS: meme.txt exists ({os.path.getsize(meme_txt)} bytes)")
        n_sanity_pass += 1
    else:
        print(f"SANITY FAIL: meme.txt missing or empty")
        n_sanity_fail += 1

    # 3. Background model exists
    if os.path.exists(result['background_model']):
        print(f"SANITY PASS: Background model exists")
        n_sanity_pass += 1
    else:
        print(f"SANITY FAIL: Background model missing")
        n_sanity_fail += 1

    # 4. Motif widths in expected range (20-30 bp)
    widths = parse_meme_motif_widths(meme_txt)
    if widths:
        all_in_range = all(20 <= w <= 30 for w in widths)
        if all_in_range:
            print(f"SANITY PASS: Motif widths {widths} all in range 20-30 bp")
            n_sanity_pass += 1
        else:
            print(f"SANITY FAIL: Motif widths {widths} out of range 20-30 bp")
            n_sanity_fail += 1
    else:
        print(f"SANITY FAIL: No motif widths found in meme.txt")
        n_sanity_fail += 1

    # 5. PWM rows sum to ~1.0
    try:
        with open(meme_txt) as f:
            content = f.read()
        # Find probability matrix lines (4 floats)
        pwm_pattern = re.compile(r'^\s*([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s*$', re.MULTILINE)
        rows = pwm_pattern.findall(content)
        if rows:
            bad_rows = 0
            for row in rows:
                total = sum(float(x) for x in row)
                if abs(total - 1.0) > 0.01:
                    bad_rows += 1
            if bad_rows == 0:
                print(f"SANITY PASS: All {len(rows)} PWM rows sum to ~1.0")
                n_sanity_pass += 1
            else:
                print(f"SANITY FAIL: {bad_rows}/{len(rows)} PWM rows don't sum to 1.0")
                n_sanity_fail += 1
        else:
            print(f"SANITY FAIL: No PWM rows found in meme.txt")
            n_sanity_fail += 1
    except Exception as e:
        print(f"SANITY FAIL: Error parsing PWM: {e}")
        n_sanity_fail += 1

    # --- SCIENTIFIC EXPECTATIONS (warn only) ---

    # 6. Mock vs real
    if result['mock_used']:
        print(f"SCIENCE WARN: Using mock MEME output (MEME not installed)")
        n_science_unexpected += 1
    else:
        print(f"SCIENCE EXPECTED: Real MEME output generated")
        n_science_expected += 1

    # 7. At least 2 motifs found
    if len(widths) >= 2:
        print(f"SCIENCE EXPECTED: {len(widths)} motifs found")
        n_science_expected += 1
    else:
        print(f"SCIENCE WARN: Only {len(widths)} motif(s) found (expected >=2)")
        n_science_unexpected += 1

    # --- Summary ---
    print(f"\n{'=' * 50}")
    print(f"SANITY:  {n_sanity_pass} pass, {n_sanity_fail} fail")
    print(f"SCIENCE: {n_science_expected} expected, {n_science_unexpected} unexpected")
    print(f"Mock used: {result['mock_used']}")

    if n_sanity_fail > 0:
        print("STOPPING: Sanity check failed.")
        sys.exit(1)
    else:
        print("All sanity checks passed.")
        sys.exit(0)
