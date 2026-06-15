"""
phase5_thermodynamic/dna_shape.py — Optional DNA shape features via DNAshapeR (R package).

Checks if R is available. If not, returns {'status': 'skipped', 'reason': 'R not installed'}.
If available, computes MGW, ProT, Roll, HelT for the operator sequence.

Functions:
  compute_dna_shape(sequence, temp_dir=None) -> dict
"""

import sys
import os
import subprocess
import tempfile

sys.path.insert(0, '/Users/aayanalwani/tb project/mce3r_stochastic')
from config.parameters import PARAMS


def _check_r_available():
    """Return True if R is installed and accessible."""
    try:
        result = subprocess.run(
            ['Rscript', '--version'],
            capture_output=True, text=True, timeout=10
        )
        return result.returncode == 0
    except (FileNotFoundError, subprocess.TimeoutExpired, OSError):
        return False


def _check_dnashaper_available():
    """Return True if DNAshapeR is available in R."""
    try:
        result = subprocess.run(
            ['Rscript', '-e', 'library(DNAshapeR); cat("OK")'],
            capture_output=True, text=True, timeout=30
        )
        return result.returncode == 0 and 'OK' in result.stdout
    except (FileNotFoundError, subprocess.TimeoutExpired, OSError):
        return False


def compute_dna_shape(sequence=None, temp_dir=None):
    """
    Compute DNA shape features for the operator sequence using DNAshapeR.

    If R or DNAshapeR is not available, returns a skipped status dict.

    Parameters
    ----------
    sequence : str — DNA sequence (default: operator_sequence from PARAMS)
    temp_dir : str — directory for temporary files (default: system temp)

    Returns
    -------
    dict with:
      status : 'ok' or 'skipped'
      reason : reason for skipping (if status == 'skipped')
      MGW    : minor groove width array (if status == 'ok')
      ProT   : propeller twist array
      Roll   : roll array
      HelT   : helix twist array
      sequence : the sequence analyzed
    """
    if sequence is None:
        sequence = PARAMS['operator_sequence']

    # Check R availability
    if not _check_r_available():
        return {
            'status': 'skipped',
            'reason': 'R not installed',
            'sequence': sequence,
        }

    # Check DNAshapeR availability
    if not _check_dnashaper_available():
        return {
            'status': 'skipped',
            'reason': 'DNAshapeR R package not installed',
            'sequence': sequence,
        }

    # R and DNAshapeR are available — compute shape features
    try:
        with tempfile.TemporaryDirectory(dir=temp_dir) as tmpdir:
            # Write sequence to FASTA file
            fasta_path = os.path.join(tmpdir, 'operator.fa')
            with open(fasta_path, 'w') as f:
                f.write(f'>operator\n{sequence}\n')

            # Write R script
            r_script = os.path.join(tmpdir, 'shape.R')
            output_prefix = os.path.join(tmpdir, 'shape_output')
            with open(r_script, 'w') as f:
                f.write(f"""
library(DNAshapeR)
fn <- "{fasta_path}"
pred <- getShape(fn)
cat("MGW:", pred$MGW[1,], "\\n")
cat("ProT:", pred$ProT[1,], "\\n")
cat("Roll:", pred$Roll[1,], "\\n")
cat("HelT:", pred$HelT[1,], "\\n")
""")

            result = subprocess.run(
                ['Rscript', r_script],
                capture_output=True, text=True, timeout=120
            )

            if result.returncode != 0:
                return {
                    'status': 'skipped',
                    'reason': f'DNAshapeR computation failed: {result.stderr[:200]}',
                    'sequence': sequence,
                }

            # Parse output
            import numpy as np
            shape_data = {}
            for line in result.stdout.split('\n'):
                for key in ['MGW', 'ProT', 'Roll', 'HelT']:
                    if line.startswith(f'{key}:'):
                        values_str = line[len(key)+1:].strip()
                        values = []
                        for v in values_str.split():
                            try:
                                values.append(float(v))
                            except ValueError:
                                pass
                        if values:
                            shape_data[key] = np.array(values)

            if not shape_data:
                return {
                    'status': 'skipped',
                    'reason': 'Could not parse DNAshapeR output',
                    'sequence': sequence,
                }

            return {
                'status': 'ok',
                'sequence': sequence,
                'MGW':  shape_data.get('MGW', None),
                'ProT': shape_data.get('ProT', None),
                'Roll': shape_data.get('Roll', None),
                'HelT': shape_data.get('HelT', None),
            }

    except Exception as e:
        return {
            'status': 'skipped',
            'reason': f'Error computing DNA shape: {e}',
            'sequence': sequence,
        }


# ---------------------------------------------------------------------------
# Self-tests
# ---------------------------------------------------------------------------
if __name__ == '__main__':
    import sys

    n_pass = 0
    n_fail = 0

    def sanity(name, cond, msg=""):
        global n_pass, n_fail
        if cond:
            print(f"  SANITY PASS: {name} {msg}")
            n_pass += 1
        else:
            print(f"  SANITY FAIL: {name} {msg}")
            n_fail += 1

    print("=== dna_shape.py self-tests ===")

    r_available = _check_r_available()
    print(f"  R available: {r_available}")

    result = compute_dna_shape()
    print(f"  Status: {result['status']}")
    if result['status'] == 'skipped':
        print(f"  Reason: {result.get('reason', 'unknown')}")

    # Test 1: Module returns dict with 'status' key (either 'ok' or 'skipped')
    sanity("Returns dict with 'status' key",
           isinstance(result, dict) and 'status' in result,
           f"— status = {result.get('status', 'MISSING')}")

    sanity("Status is 'ok' or 'skipped'",
           result['status'] in ['ok', 'skipped'],
           f"— status = {result['status']}")

    if result['status'] == 'ok':
        print(f"  MGW length: {len(result.get('MGW', []))}")
        print(f"  ProT length: {len(result.get('ProT', []))}")
        sanity("MGW array present",
               result.get('MGW') is not None,
               f"— length {len(result.get('MGW', []))}")

    print(f"\n{'='*50}")
    print(f"Sanity: {n_pass} pass, {n_fail} fail")
    if n_fail > 0:
        print("STOPPING: Sanity check failed.")
        sys.exit(1)
    else:
        print("All checks passed.")
        sys.exit(0)
