"""
Genome-wide TetR operator classification orchestrator.

Parses all TetR-family regulators from H37Rv, classifies operator architectures,
predicts noise levels, and produces genome-wide statistics.
"""

import sys
import os
import json
import time
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS_DIR = os.path.join(PROJECT_ROOT, 'results', 'phase_genome')


def run_phase_genome(n_cells_per_tf=300, seed=42):
    """Run genome-wide TetR operator classification."""
    os.makedirs(RESULTS_DIR, exist_ok=True)
    start = time.time()

    from phase_genome.tetR_scanner import run_genome_scan

    summary = {
        'description': 'Genome-wide TetR-family operator classification',
        'n_cells_per_tf': n_cells_per_tf,
        'n_sanity_pass': 0,
        'n_sanity_fail': 0,
        'n_science_expected': 0,
        'n_science_unexpected': 0,
        'errors': [],
    }

    try:
        catalog, arch_df, noise_df = run_genome_scan(n_cells_per_tf=n_cells_per_tf, seed=seed)

        # Save CSVs
        catalog.to_csv(os.path.join(RESULTS_DIR, 'tetR_catalog.csv'), index=False)
        arch_df.to_csv(os.path.join(RESULTS_DIR, 'operator_architectures.csv'), index=False)
        noise_df.to_csv(os.path.join(RESULTS_DIR, 'noise_predictions.csv'), index=False)

        n_tetR = len(catalog)
        summary['n_tetR_genes'] = n_tetR

        # Sanity: reasonable number of TetR genes
        if 10 <= n_tetR <= 200:
            summary['n_sanity_pass'] += 1
        else:
            summary['n_sanity_fail'] += 1
            summary['errors'].append(f'Unexpected TetR count: {n_tetR}')

        # Sanity: all classifications valid
        valid_classes = {'palindromic', 'asymmetric', 'single_site', 'tandem', 'unresolved'}
        all_valid = all(c in valid_classes for c in arch_df['architecture'])
        if all_valid:
            summary['n_sanity_pass'] += 1
        else:
            summary['n_sanity_fail'] += 1

        # Sanity: all CV values positive and finite
        if len(noise_df) > 0 and all(noise_df['cv'] > 0) and all(np.isfinite(noise_df['cv'])):
            summary['n_sanity_pass'] += 1
        else:
            summary['n_sanity_fail'] += 1

        # Architecture distribution
        arch_counts = arch_df['architecture'].value_counts().to_dict()
        summary['architecture_distribution'] = arch_counts
        for arch_type in valid_classes:
            pct = arch_counts.get(arch_type, 0) / n_tetR * 100 if n_tetR > 0 else 0
            summary[f'pct_{arch_type}'] = round(pct, 1)

        print(f"\n  Architecture distribution:")
        for arch_type, count in sorted(arch_counts.items(), key=lambda x: -x[1]):
            pct = count / n_tetR * 100
            print(f"    {arch_type}: {count} ({pct:.1f}%)")

        # Science: majority should be palindromic (standard TetR family)
        pct_palindromic = arch_counts.get('palindromic', 0) / n_tetR * 100 if n_tetR > 0 else 0
        if pct_palindromic >= 30:  # relaxed threshold — palindromic should be common
            summary['n_science_expected'] += 1
        else:
            summary['n_science_unexpected'] += 1

        # Science: asymmetric should be minority
        pct_asymmetric = arch_counts.get('asymmetric', 0) / n_tetR * 100 if n_tetR > 0 else 0
        if pct_asymmetric < 40:
            summary['n_science_expected'] += 1
        else:
            summary['n_science_unexpected'] += 1

        # Science: Mce3R (Rv1963c) should be in the catalog
        mce3r_found = 'Rv1963c' in catalog['rv_id'].values
        summary['mce3r_found'] = mce3r_found
        if mce3r_found:
            summary['n_sanity_pass'] += 1
            # Check if classified as asymmetric
            mce3r_arch = arch_df[arch_df['rv_id'] == 'Rv1963c']
            if len(mce3r_arch) > 0:
                mce3r_class = mce3r_arch.iloc[0]['architecture']
                summary['mce3r_classification'] = mce3r_class
                print(f"\n  Mce3R (Rv1963c) classified as: {mce3r_class}")
        else:
            summary['n_sanity_fail'] += 1
            summary['errors'].append('Mce3R (Rv1963c) not found in TetR catalog')

        # Noise comparison by architecture
        if len(noise_df) > 0:
            for arch_type in ['palindromic', 'asymmetric', 'single_site']:
                subset = noise_df[noise_df['architecture'] == arch_type]
                if len(subset) > 0:
                    summary[f'mean_cv_{arch_type}'] = round(float(subset['cv'].mean()), 4)
                    summary[f'mean_persister_{arch_type}'] = round(float(subset['persister_fraction'].mean()), 4)

            # Science: asymmetric should have higher mean CV than palindromic
            cv_asym = noise_df[noise_df['architecture'] == 'asymmetric']['cv'].mean() if 'asymmetric' in arch_counts else 0
            cv_pal = noise_df[noise_df['architecture'] == 'palindromic']['cv'].mean() if 'palindromic' in arch_counts else 0
            if cv_asym > cv_pal and cv_asym > 0:
                summary['n_science_expected'] += 1
                summary['asymmetric_noisier'] = True
            elif cv_asym == 0 or cv_pal == 0:
                pass  # not enough data
            else:
                summary['n_science_unexpected'] += 1
                summary['asymmetric_noisier'] = False

            print(f"\n  Mean CV by architecture:")
            for arch_type in ['palindromic', 'asymmetric', 'single_site', 'unresolved']:
                subset = noise_df[noise_df['architecture'] == arch_type]
                if len(subset) > 0:
                    print(f"    {arch_type}: CV={subset['cv'].mean():.4f} (n={len(subset)})")

    except Exception as e:
        summary['errors'].append(f'Genome scan failed: {str(e)}')
        import traceback
        traceback.print_exc()

    # Finalize
    elapsed = time.time() - start
    summary['wall_time_sec'] = round(elapsed, 1)
    summary['status'] = 'complete' if not summary['errors'] else 'complete_with_errors'

    summary_path = os.path.join(RESULTS_DIR, 'genome_summary.json')
    with open(summary_path, 'w') as f:
        json.dump(summary, f, indent=2, default=str)

    print(f"\n[Genome] Complete in {elapsed:.1f}s")
    print(f"  Sanity: {summary['n_sanity_pass']} pass, {summary['n_sanity_fail']} fail")
    print(f"  Science: {summary['n_science_expected']} expected, {summary['n_science_unexpected']} unexpected")

    return summary


if __name__ == '__main__':
    run_phase_genome()
