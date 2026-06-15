"""
phase6_environmental/persistence_threshold_v2.py
Corrected persister-fraction quantification with three threshold methods.

The original persistence_threshold.py used the 10th percentile of Condition D
(unregulated, ~2057 proteins) as the threshold. Because all regulated conditions
have means of 150-370, every cell fell below 2057, producing persister
fractions of 1.0 --- a useless metric.

This module provides three biologically motivated alternatives:

Method 1 - "tail_fraction":
    For each architecture, define persisters as cells in the bottom X% of that
    architecture's OWN baseline distribution. Default X = 1% (matching real
    persister fractions ~10^{-2} to 10^{-3}). Under stress, does the fraction
    of very-low-expression cells increase?

Method 2 - "absolute_calibrated":
    Back-calculate the protein count that yields a persister fraction of ~10^{-3}
    in the symmetric-baseline distribution (matching Pandey et al. 2023 wild-type
    persister frequency). Apply that single threshold to every condition.

Method 3 - "fold_change":
    Persisters are cells expressing below median / fold of the population median.
    Default fold = 1.5 (cells below 2/3 of the median). This value is calibrated
    to produce a non-trivial persister fraction (~1-2%) in the symmetric baseline
    given the CV ~0.18 of these stochastic gene-expression distributions. The
    original 1/10th fold was too extreme (zero cells qualify). Biology-agnostic;
    captures low-expressors regardless of distribution shape.
"""

import numpy as np
from typing import Dict, Tuple, Optional


# ---------------------------------------------------------------------------
# Method 1 --- tail_fraction
# ---------------------------------------------------------------------------

def threshold_tail_fraction(
    baseline_proteins: np.ndarray,
    tail_percentile: float = 1.0,
) -> float:
    """
    Return the protein count at the *tail_percentile*-th percentile of
    *baseline_proteins*.  Cells below this value in any condition are
    counted as persisters.

    Parameters
    ----------
    baseline_proteins : array of protein counts from the architecture's own
        baseline (no-stress) condition.
    tail_percentile : bottom-fraction percentile (default 1%).

    Returns
    -------
    float -- threshold in protein counts.
    """
    return float(np.percentile(baseline_proteins, tail_percentile))


# ---------------------------------------------------------------------------
# Method 2 --- absolute_calibrated
# ---------------------------------------------------------------------------

def threshold_absolute_calibrated(
    symmetric_baseline_proteins: np.ndarray,
    target_fraction: float = 1e-3,
) -> float:
    """
    Find the protein count that produces *target_fraction* persisters in the
    symmetric-baseline distribution (Condition B equivalent).

    We walk the sorted array to find the value below which exactly
    target_fraction of cells fall.

    Parameters
    ----------
    symmetric_baseline_proteins : protein counts from symmetric + baseline.
    target_fraction : desired baseline persister fraction (default 0.001).

    Returns
    -------
    float -- absolute protein-count threshold.
    """
    sorted_p = np.sort(symmetric_baseline_proteins)
    n = len(sorted_p)
    idx = max(0, int(np.floor(target_fraction * n)) - 1)
    # The threshold is the value at that index (inclusive boundary).
    # If idx == 0 and the target fraction is extremely small, use a value
    # slightly below the minimum so that at most 1 cell qualifies.
    if idx <= 0:
        return float(sorted_p[0]) - 0.5  # half-protein below minimum
    return float(sorted_p[idx])


# ---------------------------------------------------------------------------
# Method 3 --- fold_change
# ---------------------------------------------------------------------------

def threshold_fold_change(
    proteins: np.ndarray,
    fold: float = 1.5,
) -> float:
    """
    Return median(proteins) / fold.  Cells below this are persisters.

    Parameters
    ----------
    proteins : protein counts for *the condition being scored* (not a
        reference distribution).
    fold : fold below median (default 1.5, i.e. cells below 2/3 of median).

    Returns
    -------
    float -- threshold in protein counts.
    """
    med = float(np.median(proteins))
    return med / fold


# ---------------------------------------------------------------------------
# Persister-fraction computation
# ---------------------------------------------------------------------------

def persister_fraction(proteins: np.ndarray, threshold: float) -> float:
    """Fraction of cells with protein count strictly below *threshold*."""
    return float(np.sum(proteins < threshold)) / len(proteins)


# ---------------------------------------------------------------------------
# High-level driver: compute fractions under all three methods
# ---------------------------------------------------------------------------

def compute_all_methods(
    proteins: np.ndarray,
    baseline_proteins_same_arch: np.ndarray,
    symmetric_baseline_proteins: np.ndarray,
    tail_percentile: float = 1.0,
    calibrated_target: float = 1e-3,
    fold: float = 10.0,
) -> Dict[str, Tuple[float, float]]:
    """
    Compute persister fraction under all three methods.

    Returns
    -------
    dict  method_name -> (threshold_value, persister_fraction)
    """
    # Method 1
    t1 = threshold_tail_fraction(baseline_proteins_same_arch, tail_percentile)
    f1 = persister_fraction(proteins, t1)

    # Method 2
    t2 = threshold_absolute_calibrated(symmetric_baseline_proteins, calibrated_target)
    f2 = persister_fraction(proteins, t2)

    # Method 3 -- threshold is condition-specific (uses the scored population)
    t3 = threshold_fold_change(proteins, fold)
    f3 = persister_fraction(proteins, t3)

    return {
        'tail_fraction':       (t1, f1),
        'absolute_calibrated': (t2, f2),
        'fold_change':         (t3, f3),
    }
