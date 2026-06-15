#!/usr/bin/env python3
"""
main.py — Orchestrator for the Mce3R Asymmetric Operator Project

Runs the entire computational pipeline:
  Phase 1 (bioinformatics) and Phase 2 (simulation) in PARALLEL via subprocess
  Phase 3 (analysis) sequentially after Phase 2
  Phase 4 (figures) sequentially after Phase 3

Usage:
    python main.py              # Run everything
    python main.py --phase 1    # Run only Phase 1
    python main.py --phase 2    # Run only Phase 2
    python main.py --phase 3    # Run only Phase 3 (requires Phase 2 output)
    python main.py --phase 4    # Run only Phase 4 (requires Phase 3 output)
"""

import argparse
import json
import os
import subprocess
import sys
import time
from datetime import datetime
from pathlib import Path

PROJECT_ROOT = Path(__file__).parent.resolve()


def create_directories():
    """Create all required output directories."""
    dirs = [
        'data/genomes', 'data/sequences', 'data/external',
        'results/phase1/meme_output', 'results/phase1/fimo_output',
        'results/phase2/traces',
        'results/phase3',
        'results/figures',
        'logs',
    ]
    for d in dirs:
        (PROJECT_ROOT / d).mkdir(parents=True, exist_ok=True)


def read_phase_summary(phase_num):
    """Read a phase summary JSON file."""
    paths = [
        PROJECT_ROOT / f'results/phase{phase_num}/phase{phase_num}_summary.json',
        PROJECT_ROOT / f'results/figures/figures_summary.json',
    ]
    for p in paths:
        if p.exists():
            with open(p) as f:
                return json.load(f)
    return None


def verify_phase2_outputs_exist():
    """Check that Phase 2 output files exist before running Phase 3."""
    required = ['condition_A.npz', 'condition_B.npz', 'condition_C.npz',
                'condition_D.npz', 'condition_E_sweep.npz']
    phase2_dir = PROJECT_ROOT / 'results' / 'phase2'
    missing = [f for f in required if not (phase2_dir / f).exists()]
    if missing:
        raise FileNotFoundError(
            f"Phase 2 outputs missing: {missing}. Run Phase 2 first."
        )


def verify_phase3_outputs_exist():
    """Check that Phase 3 output files exist before running Phase 4."""
    required = ['noise_metrics.csv', 'bootstrap_results.csv']
    phase3_dir = PROJECT_ROOT / 'results' / 'phase3'
    missing = [f for f in required if not (phase3_dir / f).exists()]
    if missing:
        raise FileNotFoundError(
            f"Phase 3 outputs missing: {missing}. Run Phase 3 first."
        )


def run_phase_subprocess(phase_module, log_file):
    """Run a phase as a subprocess, capturing output to log file."""
    log_path = PROJECT_ROOT / 'logs' / log_file
    with open(log_path, 'w') as log_fh:
        proc = subprocess.Popen(
            [sys.executable, '-m', phase_module],
            stdout=log_fh,
            stderr=subprocess.STDOUT,
            cwd=str(PROJECT_ROOT),
        )
    return proc


def write_build_report(build_status, elapsed_sec):
    """Write BUILD_REPORT.md summarizing all phases."""
    report_path = PROJECT_ROOT / 'BUILD_REPORT.md'
    timestamp = datetime.now().strftime('%Y-%m-%d %H:%M:%S')

    # Determine overall status
    any_failed = False
    for phase_key, info in build_status.items():
        if phase_key == 'error':
            any_failed = True
            break
        if isinstance(info, dict):
            if info.get('status') == 'FAILED' or info.get('returncode', 0) != 0:
                any_failed = True
                break

    overall = 'FAILED' if any_failed else 'COMPLETE'

    lines = [
        '# Build Report\n',
        f'- **Timestamp:** {timestamp}',
        f'- **Status:** {overall}',
        f'- **Total wall time:** {elapsed_sec/60:.1f} minutes',
        '',
    ]

    # Phase summaries from JSON files
    for phase_num in [1, 2, 3]:
        summary = read_phase_summary(phase_num)
        if summary:
            # Phase 2 has a different JSON format — handle it
            if phase_num == 2 and 'status' not in summary:
                predictions = summary.get('predictions', {})
                n_cells = sum(
                    summary['conditions'][c].get('n_cells', 0)
                    for c in ['A', 'B', 'C', 'D'] if c in summary.get('conditions', {})
                )
                wt = summary.get('total_time_seconds', '?')
                cv_ok = predictions.get('cv_A_gt_cv_B', False)
                fano_d = predictions.get('fano_D_approx_expected', 0)
                fano_ok = 5 <= fano_d <= 12
                status = 'PASS' if cv_ok and fano_ok else 'WARN'
                lines.append(
                    f'- **Phase {phase_num}:** {status} | '
                    f'{n_cells} cells simulated | '
                    f'CV(A)>CV(B)={cv_ok} | Fano(D)={fano_d:.2f} | '
                    f'Wall time: {wt}s'
                )
            else:
                status = summary.get('status', 'UNKNOWN')
                sp = summary.get('n_sanity_pass', '?')
                sf = summary.get('n_sanity_fail', '?')
                se = summary.get('n_science_expected', '?')
                su = summary.get('n_science_unexpected', '?')
                mock = summary.get('mock_used', False)
                wt = summary.get('wall_time_sec', '?')
                mock_flag = ' **[MOCK]**' if mock else ''
                lines.append(
                    f'- **Phase {phase_num}:** {status} | '
                    f'Sanity: {sp} pass / {sf} fail | '
                    f'Science: {se} expected / {su} unexpected | '
                    f'Wall time: {wt}s{mock_flag}'
                )
        elif phase_num in [1, 2]:
            # Check subprocess return code
            info = build_status.get(f'phase{phase_num}', {})
            rc = info.get('returncode', '?')
            lines.append(f'- **Phase {phase_num}:** subprocess returncode={rc}')

    # Phase 4 (figures) — JSON uses n_success/n_error/total keys
    fig_summary_path = PROJECT_ROOT / 'results' / 'figures' / 'figures_summary.json'
    if fig_summary_path.exists():
        with open(fig_summary_path) as f:
            fig_summary = json.load(f)
        n_success = fig_summary.get('n_success', 0)
        n_error = fig_summary.get('n_error', 0)
        total = fig_summary.get('total', '?')
        fig_status = 'PASS' if n_error == 0 and n_success > 0 else 'FAIL'
        lines.append(
            f'- **Phase 4:** {fig_status} | '
            f'{n_success}/{total} figures generated, {n_error} errors'
        )

    # Error info
    if 'error' in build_status:
        lines.extend(['', f'- **Error:** {build_status["error"]}'])

    # Output files
    lines.extend(['', '## Output Files', ''])
    for dirpath in ['results/phase1', 'results/phase2', 'results/phase3', 'results/figures']:
        full = PROJECT_ROOT / dirpath
        if full.exists():
            for f in sorted(full.iterdir()):
                if f.is_file():
                    size_kb = f.stat().st_size / 1024
                    lines.append(f'- `{dirpath}/{f.name}` ({size_kb:.1f} KB)')

    lines.append('')

    with open(report_path, 'w') as f:
        f.write('\n'.join(lines))

    print(f"[Orchestrator] BUILD_REPORT.md written to {report_path}")


def main():
    parser = argparse.ArgumentParser(
        description='Mce3R Asymmetric Operator Project — Full Pipeline'
    )
    parser.add_argument(
        '--phase', type=int, choices=[1, 2, 3, 4],
        help='Run specific phase only (default: run all)'
    )
    args = parser.parse_args()

    # Ensure we're in the project root
    os.chdir(str(PROJECT_ROOT))

    # Create directory structure
    create_directories()

    build_status = {}
    start_time = time.time()

    try:
        if args.phase is None:
            # === PARALLEL: Phase 1 + Phase 2 ===
            print("[Orchestrator] Launching Phase 1 + Phase 2 in parallel...")
            p1 = run_phase_subprocess(
                'phase1_pipeline.pipeline_main', 'phase1_pipeline.log'
            )
            p2 = run_phase_subprocess(
                'phase2_simulation.simulation_main', 'phase2_gillespie.log'
            )

            p1.wait()
            p2.wait()

            build_status['phase1'] = {'returncode': p1.returncode}
            build_status['phase2'] = {'returncode': p2.returncode}

            print(f"[Orchestrator] Phase 1 exited with code {p1.returncode}")
            print(f"[Orchestrator] Phase 2 exited with code {p2.returncode}")

            if p2.returncode != 0:
                raise RuntimeError(
                    f"Phase 2 failed (returncode={p2.returncode}). "
                    f"Check logs/phase2_gillespie.log"
                )

            # === SEQUENTIAL: Phase 3 (depends on Phase 2) ===
            print("[Orchestrator] Starting Phase 3: Statistical Validation")
            verify_phase2_outputs_exist()
            from phase3_analysis.analysis_main import run_phase3
            build_status['phase3'] = run_phase3()

            # === SEQUENTIAL: Phase 4 (depends on Phase 3) ===
            print("[Orchestrator] Starting Phase 4: Publication Figures")
            verify_phase3_outputs_exist()
            from figures.figures_main import run_phase4
            build_status['phase4'] = run_phase4()

        elif args.phase == 1:
            print("[Orchestrator] Running Phase 1 only")
            from phase1_pipeline.pipeline_main import run_phase1
            build_status['phase1'] = run_phase1()

        elif args.phase == 2:
            print("[Orchestrator] Running Phase 2 only")
            from phase2_simulation.simulation_main import run_phase2
            build_status['phase2'] = run_phase2()

        elif args.phase == 3:
            print("[Orchestrator] Running Phase 3 only")
            verify_phase2_outputs_exist()
            from phase3_analysis.analysis_main import run_phase3
            build_status['phase3'] = run_phase3()

        elif args.phase == 4:
            print("[Orchestrator] Running Phase 4 only")
            verify_phase3_outputs_exist()
            from figures.figures_main import run_phase4
            build_status['phase4'] = run_phase4()

    except Exception as e:
        build_status['error'] = str(e)
        print(f"[Orchestrator] FATAL ERROR: {e}")
        import traceback
        traceback.print_exc()

    # Write BUILD_REPORT.md
    elapsed = time.time() - start_time
    write_build_report(build_status, elapsed)

    print(f"\n{'='*60}")
    print(f"[Orchestrator] Complete. Total time: {elapsed/60:.1f} minutes")
    print(f"[Orchestrator] See BUILD_REPORT.md for details")
    print(f"{'='*60}")


if __name__ == '__main__':
    main()
