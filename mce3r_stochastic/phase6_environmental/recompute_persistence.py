#!/usr/bin/env python3
"""
phase6_environmental/recompute_persistence.py
Recompute persister fractions for all 24 phase-6 NPZ files using the three
corrected threshold methods from persistence_threshold_v2.

Usage:
    python -m phase6_environmental.recompute_persistence

Outputs:
    results/phase6/persistence_fractions_v2.csv
    Summary table printed to stdout.
"""

import sys, os, glob, re
import numpy as np
import pandas as pd

# Ensure project root is on sys.path
PROJECT_ROOT = '/Users/aayanalwani/tb project/mce3r_stochastic'
sys.path.insert(0, PROJECT_ROOT)

from phase6_environmental.persistence_threshold_v2 import (
    threshold_tail_fraction,
    threshold_absolute_calibrated,
    threshold_fold_change,
    persister_fraction,
    compute_all_methods,
)

NPZ_DIR = os.path.join(PROJECT_ROOT, 'results', 'phase6')
OUT_CSV  = os.path.join(NPZ_DIR, 'persistence_fractions_v2.csv')


# ── helpers ─────────────────────────────────────────────────────────────────

def _parse_filename(path: str):
    """
    Parse  env_condition_{arch}_{env}_{model}.npz
    Returns (architecture, environment, model_type) or None.
    """
    base = os.path.basename(path).replace('.npz', '')
    if not base.startswith('env_condition_'):
        return None
    rest = base[len('env_condition_'):]

    # model_type is the last token: 'single' or 'two_species'
    if rest.endswith('_two_species'):
        model_type = 'two_species'
        rest = rest[:-len('_two_species')]
    elif rest.endswith('_single'):
        model_type = 'single'
        rest = rest[:-len('_single')]
    else:
        return None

    # Now rest = {architecture}_{environment}
    # Architectures: asymmetric, symmetric, single_site
    for arch in ('asymmetric', 'symmetric', 'single_site'):
        if rest.startswith(arch + '_'):
            env = rest[len(arch) + 1:]
            return arch, env, model_type
    return None


def _get_proteins(npz) -> np.ndarray:
    """Return the target-protein array from an NPZ file."""
    if 'target_proteins' in npz:
        return npz['target_proteins'].astype(np.float64)
    elif 'proteins' in npz:
        return npz['proteins'].astype(np.float64)
    else:
        raise KeyError(f"No protein array found. Keys: {list(npz.keys())}")


# ── main ────────────────────────────────────────────────────────────────────

def main():
    print("=" * 72)
    print("Recomputing persister fractions with corrected thresholds (v2)")
    print("=" * 72)

    # 1. Load all 24 NPZ files --------------------------------------------------
    npz_files = sorted(glob.glob(os.path.join(NPZ_DIR, 'env_condition_*.npz')))
    print(f"\nFound {len(npz_files)} NPZ files in {NPZ_DIR}")

    data = {}  # (arch, env, model) -> np.ndarray of proteins
    for path in npz_files:
        parsed = _parse_filename(path)
        if parsed is None:
            print(f"  SKIP (unparseable): {os.path.basename(path)}")
            continue
        arch, env, model = parsed
        npz = np.load(path)
        proteins = _get_proteins(npz)
        data[(arch, env, model)] = proteins

    print(f"Loaded {len(data)} condition arrays.\n")

    # 2. Gather baseline arrays per architecture ---------------------------------
    baselines = {}  # arch -> proteins (using 'single' model baseline preferentially)
    for (arch, env, model), proteins in data.items():
        if env == 'baseline':
            key = arch
            # Prefer 'single' model, but take two_species if that's all there is
            if key not in baselines or model == 'single':
                baselines[key] = proteins

    print("Baseline distributions per architecture:")
    for arch, p in baselines.items():
        print(f"  {arch:15s}  n={len(p):6d}  mean={p.mean():.1f}  "
              f"p1={np.percentile(p, 1):.1f}  min={p.min():.0f}")

    # 3. Symmetric baseline for Method 2 ----------------------------------------
    sym_base = baselines.get('symmetric')
    if sym_base is None:
        print("WARNING: symmetric baseline not found; falling back to first baseline.")
        sym_base = next(iter(baselines.values()))

    # Pre-compute the absolute_calibrated threshold once (it is global).
    abs_thresh = threshold_absolute_calibrated(sym_base, target_fraction=1e-3)
    print(f"\nMethod 2 (absolute_calibrated) threshold = {abs_thresh:.1f} proteins")
    print(f"  (produces fraction {persister_fraction(sym_base, abs_thresh):.4f} "
          f"in symmetric baseline, target 0.001)\n")

    # 4. Compute persister fractions for every (arch, env, model) ----------------
    # Also compute the symmetric-baseline fraction per method for fold-enrichment.
    sym_base_fracs = {}  # method -> persister_fraction in symmetric baseline single
    sym_key = ('symmetric', 'baseline', 'single')
    if sym_key in data:
        sym_p = data[sym_key]
        t1 = threshold_tail_fraction(baselines['symmetric'], 1.0)
        sym_base_fracs['tail_fraction'] = persister_fraction(sym_p, t1)
        sym_base_fracs['absolute_calibrated'] = persister_fraction(sym_p, abs_thresh)
        t3 = threshold_fold_change(sym_p, 1.5)
        sym_base_fracs['fold_change'] = persister_fraction(sym_p, t3)
    else:
        # fallback
        sym_base_fracs = {'tail_fraction': 0.01, 'absolute_calibrated': 0.001,
                          'fold_change': 0.01}

    print("Symmetric-baseline persister fractions (reference for fold-enrichment):")
    for m, f in sym_base_fracs.items():
        print(f"  {m:25s}  {f:.6f}")
    print()

    rows = []
    for (arch, env, model), proteins in sorted(data.items()):
        baseline_same_arch = baselines.get(arch, sym_base)
        results = compute_all_methods(
            proteins=proteins,
            baseline_proteins_same_arch=baseline_same_arch,
            symmetric_baseline_proteins=sym_base,
            tail_percentile=1.0,
            calibrated_target=1e-3,
            fold=1.5,
        )
        for method, (thresh_val, frac) in results.items():
            ref = sym_base_fracs.get(method, 1e-9)
            fold_enrich = frac / ref if ref > 0 else np.nan
            rows.append({
                'architecture': arch,
                'environment': env,
                'model_type': model,
                'method': method,
                'threshold_value': round(thresh_val, 2),
                'persister_fraction': frac,
                'fold_enrichment_vs_symmetric': round(fold_enrich, 4),
            })

    df = pd.DataFrame(rows)
    df = df.sort_values(['method', 'architecture', 'environment', 'model_type']).reset_index(drop=True)

    # 5. Save to CSV -------------------------------------------------------------
    os.makedirs(os.path.dirname(OUT_CSV), exist_ok=True)
    df.to_csv(OUT_CSV, index=False)
    print(f"Saved {len(df)} rows to {OUT_CSV}\n")

    # 6. Summary table -----------------------------------------------------------
    print("=" * 90)
    print(f"{'Method':<25s} {'Architecture':<15s} {'Environment':<15s} "
          f"{'Model':<14s} {'Thresh':>8s} {'Frac':>10s} {'Fold':>8s}")
    print("-" * 90)
    for _, r in df.iterrows():
        print(f"{r['method']:<25s} {r['architecture']:<15s} {r['environment']:<15s} "
              f"{r['model_type']:<14s} {r['threshold_value']:8.1f} "
              f"{r['persister_fraction']:10.6f} {r['fold_enrichment_vs_symmetric']:8.4f}")
    print("=" * 90)

    # 7. Self-tests --------------------------------------------------------------
    print("\n--- Self-tests ---")
    n_pass = 0
    n_fail = 0

    def check(name, cond, msg=""):
        nonlocal n_pass, n_fail
        if cond:
            print(f"  PASS: {name} {msg}")
            n_pass += 1
        else:
            print(f"  FAIL: {name} {msg}")
            n_fail += 1

    # Test 1: No persister fraction is exactly 1.0 (the old bug)
    any_one = (df['persister_fraction'] == 1.0).any()
    check("No persister fraction is 1.0 (old bug fixed)",
          not any_one,
          f"-- max fraction = {df['persister_fraction'].max():.6f}")

    # Test 2: At least one method has fractions in [1e-4, 0.1]
    sensible = df[(df['persister_fraction'] >= 1e-4) &
                  (df['persister_fraction'] <= 0.1)]
    check("Persister fractions in [1e-4, 0.1] for at least one method",
          len(sensible) > 0,
          f"-- {len(sensible)} rows in range")

    # Test 3: Asymmetric fold_enrichment > 1 for at least one environment
    asym_enriched = df[(df['architecture'] == 'asymmetric') &
                       (df['fold_enrichment_vs_symmetric'] > 1.0)]
    check("Asymmetric fold_enrichment > 1 for at least one env",
          len(asym_enriched) > 0,
          f"-- {len(asym_enriched)} rows with fold > 1")

    print(f"\nSelf-tests: {n_pass} passed, {n_fail} failed.")
    if n_fail > 0:
        sys.exit(1)
    else:
        print("All self-tests passed.")
        sys.exit(0)


if __name__ == '__main__':
    main()
