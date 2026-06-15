"""
phase1_pipeline/pipeline_main.py — Orchestrate the full Phase 1 bioinformatics pipeline.

Runs all steps in sequence:
  1. download_genomes
  2. extract_upstream
  3. run_meme
  4. run_fimo
  5. conservation_check

Collects results, writes phase1_summary.json, and returns structured dict.
"""

import sys
import os
import json
import time
import traceback

sys.path.insert(0, '/Users/aayanalwani/tb project/mce3r_stochastic')
from config.parameters import PARAMS

# Import pipeline steps
from phase1_pipeline.download_genomes import run_download
from phase1_pipeline.extract_upstream import run_extract_upstream
from phase1_pipeline.run_meme import run_meme_step
from phase1_pipeline.run_fimo import run_fimo_step
from phase1_pipeline.conservation_check import run_conservation_check

PROJECT_ROOT = '/Users/aayanalwani/tb project/mce3r_stochastic'
RESULTS_DIR = os.path.join(PROJECT_ROOT, 'results', 'phase1')


def run_pipeline():
    """
    Run the full Phase 1 pipeline.
    Returns a structured dict with:
      status, n_sanity_pass, n_sanity_fail, n_science_expected,
      n_science_unexpected, n_warn, mock_used, wall_time_sec, output_files
    """
    os.makedirs(RESULTS_DIR, exist_ok=True)
    start_time = time.time()

    n_sanity_pass = 0
    n_sanity_fail = 0
    n_science_expected = 0
    n_science_unexpected = 0
    n_warn = 0
    mock_used = False
    output_files = []
    step_results = {}

    steps = [
        ('download_genomes', run_download),
        ('extract_upstream', run_extract_upstream),
        ('run_meme', run_meme_step),
        ('run_fimo', run_fimo_step),
        ('conservation_check', run_conservation_check),
    ]

    for step_name, step_func in steps:
        print(f"\n{'='*60}")
        print(f"STEP: {step_name}")
        print(f"{'='*60}")
        try:
            result = step_func()
            step_results[step_name] = result
            n_sanity_pass += 1
            print(f"  {step_name}: COMPLETED")

            # Collect output files
            for key, val in result.items():
                if isinstance(val, str) and os.path.exists(val):
                    output_files.append(val)

            # Track mock usage
            if result.get('mock_used', False):
                mock_used = True
                n_warn += 1

        except Exception as e:
            print(f"  {step_name}: FAILED — {e}")
            traceback.print_exc()
            n_sanity_fail += 1
            step_results[step_name] = {'error': str(e)}

    wall_time = time.time() - start_time

    # --- Post-pipeline sanity checks ---
    print(f"\n{'='*60}")
    print("POST-PIPELINE CHECKS")
    print(f"{'='*60}")

    # Check key output files exist
    expected_files = [
        ('H37Rv genome', os.path.join(PROJECT_ROOT, 'data', 'genomes', 'H37Rv.fasta')),
        ('All upstream', os.path.join(PROJECT_ROOT, 'data', 'sequences', 'all_upstream_200bp.fasta')),
        ('MEME output', os.path.join(RESULTS_DIR, 'meme_output', 'meme.txt')),
        ('Predicted sites', os.path.join(RESULTS_DIR, 'predicted_sites.csv')),
        ('Conservation', os.path.join(RESULTS_DIR, 'conservation_status.csv')),
    ]

    for name, path in expected_files:
        if os.path.exists(path):
            print(f"  SANITY PASS: {name} exists")
            n_sanity_pass += 1
        else:
            print(f"  SANITY FAIL: {name} missing: {path}")
            n_sanity_fail += 1

    # Scientific checks on final results
    # Check predicted sites has known operator
    predicted_csv = os.path.join(RESULTS_DIR, 'predicted_sites.csv')
    if os.path.exists(predicted_csv):
        import csv
        with open(predicted_csv) as f:
            reader = csv.DictReader(f)
            rows = list(reader)
        op_start = PARAMS['operator_region_h37rv_start']
        op_end = PARAMS['operator_region_h37rv_end']
        found = any(op_start - 500 <= int(r['start']) <= op_end + 500 for r in rows)
        if found:
            print(f"  SCIENCE EXPECTED: Known operator found in predicted sites")
            n_science_expected += 1
        else:
            print(f"  SCIENCE WARN: Known operator NOT in predicted sites")
            n_science_unexpected += 1

    # Check conservation has reasonable results
    conservation_csv = os.path.join(RESULTS_DIR, 'conservation_status.csv')
    if os.path.exists(conservation_csv):
        import csv
        with open(conservation_csv) as f:
            reader = csv.DictReader(f)
            cons_rows = list(reader)
        n_bovis = sum(1 for r in cons_rows if r['bovis_conserved'] == 'True')
        if n_bovis > 0:
            print(f"  SCIENCE EXPECTED: {n_bovis} sites conserved in M. bovis")
            n_science_expected += 1
        else:
            print(f"  SCIENCE WARN: No sites conserved in M. bovis")
            n_science_unexpected += 1

    # Determine status
    if n_sanity_fail > 0:
        status = 'FAILED'
    elif n_warn > 0:
        status = 'COMPLETED_WITH_WARNINGS'
    else:
        status = 'COMPLETED'

    summary = {
        'status': status,
        'n_sanity_pass': n_sanity_pass,
        'n_sanity_fail': n_sanity_fail,
        'n_science_expected': n_science_expected,
        'n_science_unexpected': n_science_unexpected,
        'n_warn': n_warn,
        'mock_used': mock_used,
        'wall_time_sec': round(wall_time, 2),
        'output_files': output_files,
    }

    # Write summary JSON
    summary_path = os.path.join(RESULTS_DIR, 'phase1_summary.json')
    with open(summary_path, 'w') as f:
        json.dump(summary, f, indent=2)
    print(f"\n  Summary written: {summary_path}")

    return summary


if __name__ == "__main__":
    print("=" * 60)
    print("pipeline_main.py — Phase 1 Pipeline Self-test")
    print("=" * 60)

    n_sanity_pass = 0
    n_sanity_fail = 0
    n_science_expected = 0
    n_science_unexpected = 0

    try:
        summary = run_pipeline()
    except Exception as e:
        traceback.print_exc()
        print(f"SANITY FAIL: Pipeline raised exception: {e}")
        sys.exit(1)

    # --- SANITY CHECKS on the summary dict ---

    # 1. Summary dict has all required keys
    required_keys = [
        'status', 'n_sanity_pass', 'n_sanity_fail', 'n_science_expected',
        'n_science_unexpected', 'n_warn', 'mock_used', 'wall_time_sec', 'output_files'
    ]
    missing_keys = [k for k in required_keys if k not in summary]
    if not missing_keys:
        print(f"SANITY PASS: Summary dict has all required keys")
        n_sanity_pass += 1
    else:
        print(f"SANITY FAIL: Missing keys in summary: {missing_keys}")
        n_sanity_fail += 1

    # 2. Summary JSON file exists
    summary_path = os.path.join(RESULTS_DIR, 'phase1_summary.json')
    if os.path.exists(summary_path):
        with open(summary_path) as f:
            loaded = json.load(f)
        if loaded == summary:
            print(f"SANITY PASS: phase1_summary.json matches returned dict")
            n_sanity_pass += 1
        else:
            print(f"SANITY FAIL: phase1_summary.json doesn't match returned dict")
            n_sanity_fail += 1
    else:
        print(f"SANITY FAIL: phase1_summary.json not found")
        n_sanity_fail += 1

    # 3. Pipeline completed (no step failures)
    if summary['n_sanity_fail'] == 0:
        print(f"SANITY PASS: All pipeline steps completed without sanity failures")
        n_sanity_pass += 1
    else:
        print(f"SANITY FAIL: {summary['n_sanity_fail']} pipeline sanity failures")
        n_sanity_fail += 1

    # 4. Wall time is reasonable (< 30 minutes)
    if summary['wall_time_sec'] < 1800:
        print(f"SANITY PASS: Wall time {summary['wall_time_sec']:.1f}s (< 30 min)")
        n_sanity_pass += 1
    else:
        print(f"SANITY FAIL: Wall time {summary['wall_time_sec']:.1f}s exceeds 30 min")
        n_sanity_fail += 1

    # 5. Output files list is populated
    if len(summary['output_files']) >= 3:
        print(f"SANITY PASS: {len(summary['output_files'])} output files tracked")
        n_sanity_pass += 1
    else:
        print(f"SANITY FAIL: Only {len(summary['output_files'])} output files tracked")
        n_sanity_fail += 1

    # --- SCIENTIFIC EXPECTATIONS ---

    # 6. Some scientific expectations were met
    if summary['n_science_expected'] > 0:
        print(f"SCIENCE EXPECTED: {summary['n_science_expected']} scientific checks passed")
        n_science_expected += 1
    else:
        print(f"SCIENCE WARN: No scientific expectations met")
        n_science_unexpected += 1

    # --- Summary ---
    print(f"\n{'='*60}")
    print(f"PIPELINE STATUS: {summary['status']}")
    print(f"Wall time: {summary['wall_time_sec']:.1f} seconds")
    print(f"Mock tools used: {summary['mock_used']}")
    print(f"Output files: {len(summary['output_files'])}")
    print(f"\nSelf-test SANITY:  {n_sanity_pass} pass, {n_sanity_fail} fail")
    print(f"Self-test SCIENCE: {n_science_expected} expected, {n_science_unexpected} unexpected")

    if n_sanity_fail > 0:
        print("\nSTOPPING: Self-test sanity check failed.")
        sys.exit(1)
    else:
        print("\nAll self-test sanity checks passed.")
        sys.exit(0)
