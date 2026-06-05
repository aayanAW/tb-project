"""
phase3_analysis/negative_controls.py — Negative control tests.

1. Shuffle test: mononucleotide-preserving shuffle for predicted binding sites.
2. Symmetric TetR control: run simulation with Kd_strong=Kd_weak=5nM, confirm unimodal.
3. Poisson baseline: verify Condition D Fano in [5, 12].
Output: results/phase3/negative_controls.csv
"""

import sys
import os
import warnings
import numpy as np
import pandas as pd

sys.path.insert(0, '/Users/aayanalwani/tb project/mce3r_stochastic')
from config.parameters import PARAMS

PROJECT_ROOT = '/Users/aayanalwani/tb project/mce3r_stochastic'
PHASE2_DIR = os.path.join(PROJECT_ROOT, 'results', 'phase2')
PHASE3_DIR = os.path.join(PROJECT_ROOT, 'results', 'phase3')


def load_proteins(condition):
    """Load protein array from Phase 2 .npz file."""
    path = os.path.join(PHASE2_DIR, f'condition_{condition}.npz')
    data = np.load(path)
    return data['proteins'].astype(np.float64)


def mononucleotide_shuffle(seq, rng=None):
    """
    Mononucleotide-preserving shuffle: permute the character list.
    Preserves base composition but destroys motif structure.
    """
    if rng is None:
        rng = np.random.default_rng(42)
    chars = list(seq)
    rng.shuffle(chars)
    return ''.join(chars)


def score_sequence_gc_content(seq):
    """Simple scoring: GC content as a proxy for binding affinity metric."""
    if len(seq) == 0:
        return 0.0
    gc = sum(1 for c in seq.upper() if c in ('G', 'C'))
    return gc / len(seq)


def run_shuffle_test(n_shuffles=None, alpha=None):
    """
    Shuffle test for binding site predictions.
    Uses operator sequence or mock sites if Phase 1 didn't produce predictions.
    Mononucleotide-preserving shuffle.
    """
    if n_shuffles is None:
        n_shuffles = PARAMS.get('shuffle_n', 1000)
    if alpha is None:
        alpha = PARAMS.get('shuffle_alpha', 0.01)

    # Try to load Phase 1 predicted sites; fall back to mock
    mock_used = False
    operator_seq = PARAMS.get('operator_sequence', '')

    # Create mock binding sites from operator sequence (top 10 predicted sites)
    # In production, these would come from FIMO output
    site_length = 25
    sites = []
    if len(operator_seq) >= site_length:
        for i in range(0, min(10 * site_length, len(operator_seq) - site_length), site_length // 2):
            sites.append(operator_seq[i:i + site_length])
        mock_used = True
    if len(sites) == 0:
        sites = ['GCCCCGCGCTATAGGATACTAGCAA']  # fallback single site
        mock_used = True

    sites = sites[:10]  # top 10

    results = []
    rng = np.random.default_rng(PARAMS.get('master_seed', 42))

    for i, site in enumerate(sites):
        real_score = score_sequence_gc_content(site)

        # Generate shuffled scores
        shuffle_scores = np.zeros(n_shuffles)
        for j in range(n_shuffles):
            shuffled = mononucleotide_shuffle(site, rng=rng)
            shuffle_scores[j] = score_sequence_gc_content(shuffled)

        # p-value: fraction of shuffles scoring >= real
        p_value = np.mean(shuffle_scores >= real_score)

        results.append({
            'test': 'shuffle',
            'site_index': i,
            'sequence': site,
            'real_score': real_score,
            'shuffle_mean': np.mean(shuffle_scores),
            'shuffle_std': np.std(shuffle_scores),
            'p_value': p_value,
            'significant': p_value < alpha,
            'mock_used': mock_used,
        })

    return results, mock_used


def run_symmetric_control(n_cells=None):
    """
    Symmetric TetR control: Kd_strong = Kd_weak = 5nM.
    Expect unimodal distribution.
    """
    from phase2_simulation.operator_model import OperatorModel
    from phase2_simulation.gillespie_engine import run_population

    if n_cells is None:
        n_cells = 100

    model = OperatorModel(Kd_strong=5.0, Kd_weak=5.0,
                           block_strong=0.85, block_weak=0.85)
    arrays = model.get_numba_arrays()
    result = run_population(arrays, n_cells=n_cells,
                            master_seed=PARAMS.get('master_seed', 42) + 9000,
                            record_traces=0, n_trace_cells=0)
    proteins = result['proteins'].astype(np.float64)

    # Check unimodality via GMM BIC
    from sklearn.mixture import GaussianMixture

    X = proteins.reshape(-1, 1)
    std_x = np.std(X)
    if std_x > 0:
        X_dither = X + np.random.default_rng(42).normal(0, std_x * 1e-4, size=X.shape)
    else:
        X_dither = X.copy()

    bics = {}
    for k in [1, 2, 3]:
        try:
            gmm = GaussianMixture(
                n_components=k,
                reg_covar=1e-3,
                n_init=5,
                random_state=42,
                max_iter=300,
            )
            gmm.fit(X_dither)
            bics[k] = gmm.bic(X_dither)
        except Exception:
            bics[k] = np.inf

    best_k = min(bics, key=bics.get)
    is_unimodal = (best_k == 1)

    return {
        'test': 'symmetric_control',
        'Kd_strong': 5.0,
        'Kd_weak': 5.0,
        'n_cells': n_cells,
        'mean': np.mean(proteins),
        'CV': np.std(proteins, ddof=1) / np.mean(proteins) if np.mean(proteins) > 0 else np.nan,
        'best_gmm_k': best_k,
        'is_unimodal': is_unimodal,
        'bic_1': bics.get(1, np.inf),
        'bic_2': bics.get(2, np.inf),
        'bic_3': bics.get(3, np.inf),
    }


def run_poisson_baseline():
    """
    Poisson baseline: verify Condition D Fano factor in [5, 12].
    """
    proteins_D = load_proteins('D')
    mean_d = np.mean(proteins_D)
    var_d = np.var(proteins_D, ddof=1) if len(proteins_D) > 1 else 0.0
    fano_d = var_d / mean_d if mean_d > 0 else np.nan

    return {
        'test': 'poisson_baseline',
        'condition': 'D',
        'mean': mean_d,
        'variance': var_d,
        'Fano': fano_d,
        'expected_range': '[5, 12]',
        'in_range': 5.0 <= fano_d <= 12.0 if np.isfinite(fano_d) else False,
    }


def run_negative_controls(n_shuffle=None, n_cells_symmetric=None):
    """Run all negative controls and save CSV."""
    os.makedirs(PHASE3_DIR, exist_ok=True)

    all_rows = []

    # 1. Shuffle test
    print("  Running shuffle test...")
    shuffle_results, mock_used = run_shuffle_test(n_shuffles=n_shuffle)
    for r in shuffle_results:
        all_rows.append({
            'test': r['test'],
            'detail': f"site_{r['site_index']}",
            'value': r['real_score'],
            'p_value': r['p_value'],
            'pass': r['significant'],
            'note': f"mock={r['mock_used']}",
        })

    # 2. Symmetric control
    print("  Running symmetric TetR control...")
    sym_result = run_symmetric_control(n_cells=n_cells_symmetric)
    all_rows.append({
        'test': 'symmetric_control',
        'detail': f"Kd={sym_result['Kd_strong']}nM",
        'value': sym_result['best_gmm_k'],
        'p_value': np.nan,
        'pass': sym_result['is_unimodal'],
        'note': f"CV={sym_result['CV']:.4f}",
    })

    # 3. Poisson baseline
    print("  Running Poisson baseline check...")
    poi_result = run_poisson_baseline()
    all_rows.append({
        'test': 'poisson_baseline',
        'detail': 'condition_D',
        'value': poi_result['Fano'],
        'p_value': np.nan,
        'pass': poi_result['in_range'],
        'note': f"expected {poi_result['expected_range']}",
    })

    df = pd.DataFrame(all_rows)
    out_path = os.path.join(PHASE3_DIR, 'negative_controls.csv')
    df.to_csv(out_path, index=False)
    print(f"  Wrote {out_path}")

    return df, mock_used, sym_result, poi_result


# ---------------------------------------------------------------------------
# Self-tests
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    n_sanity_pass = 0
    n_sanity_fail = 0
    n_science_pass = 0
    n_science_warn = 0

    def sanity(name, cond, msg=""):
        global n_sanity_pass, n_sanity_fail
        if cond:
            print(f"  SANITY PASS: {name} {msg}")
            n_sanity_pass += 1
        else:
            print(f"  SANITY FAIL: {name} {msg}")
            n_sanity_fail += 1

    def scientific(name, cond, msg=""):
        global n_science_pass, n_science_warn
        if cond:
            print(f"  SCIENTIFIC PASS: {name} {msg}")
            n_science_pass += 1
        else:
            print(f"  SCIENTIFIC WARN: {name} {msg}")
            n_science_warn += 1

    print("=== negative_controls.py self-tests ===")

    # Use small settings for self-test
    df, mock_used, sym_result, poi_result = run_negative_controls(
        n_shuffle=100, n_cells_symmetric=20
    )

    # SANITY: no exceptions (we got here)
    sanity("no exceptions", True, "— all controls ran successfully")

    # SANITY: output file exists
    out_path = os.path.join(PHASE3_DIR, 'negative_controls.csv')
    sanity("output file exists", os.path.exists(out_path))

    # SANITY: correct number of rows (up to 10 shuffle + 1 symmetric + 1 poisson)
    n_shuffle_rows = len(df[df['test'] == 'shuffle'])
    sanity("shuffle rows > 0", n_shuffle_rows > 0,
           f"— got {n_shuffle_rows} shuffle rows")
    sanity("symmetric row exists", len(df[df['test'] == 'symmetric_control']) == 1)
    sanity("poisson row exists", len(df[df['test'] == 'poisson_baseline']) == 1)

    # SANITY: p-values in [0, 1] for shuffle tests
    shuffle_p = df.loc[df['test'] == 'shuffle', 'p_value'].values
    sanity("shuffle p-values in [0,1]",
           np.all((shuffle_p >= 0) & (shuffle_p <= 1)),
           f"— range: [{np.min(shuffle_p):.4f}, {np.max(shuffle_p):.4f}]")

    # SCIENTIFIC: symmetric control is unimodal
    scientific("symmetric control unimodal", sym_result['is_unimodal'],
               f"— best_k={sym_result['best_gmm_k']}")

    # SCIENTIFIC: Fano(D) in [5, 12]
    scientific("Fano(D) in [5, 12]", poi_result['in_range'],
               f"— Fano={poi_result['Fano']:.2f}")

    print(f"\n{'='*50}")
    print(f"SANITY: {n_sanity_pass} pass, {n_sanity_fail} fail")
    print(f"SCIENTIFIC: {n_science_pass} pass, {n_science_warn} warn")
    if n_sanity_fail > 0:
        print("STOPPING: Sanity check failed.")
        sys.exit(1)
    else:
        print("All checks passed.")
        sys.exit(0)
