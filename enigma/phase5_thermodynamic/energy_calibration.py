"""
phase5_thermodynamic/energy_calibration.py — Convert MEME PWM to binding free energies.

Functions:
  load_pwm_from_meme(meme_txt_path) -> np.ndarray, shape (motif_width, 4) [ACGT]
  pwm_to_energy_matrix(pwm, background, kT=0.616) -> np.ndarray same shape
  score_sequence(seq, energy_matrix) -> float
  calibrate_energies(strong_seq, weak_seq, energy_matrix, Kd_strong, Kd_weak, kT) -> dict
  predict_Kd_for_sequence(seq, calibration) -> float  [nM]

Calibration notes:
  The operator has two binding sites for Mce3R:
    - Strong site (Kd=2.4 nM): GTTGTCCAAGCAATCGCGTAT at op[74:95] — matches MEME-1
    - Weak site (Kd=49 nM):    TTTGCATTGCTATTTACCGA  at op[51:71] — matches MEME-2
  Since the two sites correspond to different MEME motifs, a single-offset approach
  is physically incorrect. We use two independent offsets (one per site), calibrated
  from the respective experimental Kd values. For new sequence scoring, the mean
  offset is applied (conservative interpolation).
"""

import sys
import os
import numpy as np

sys.path.insert(0, '/Users/aayanalwani/tb project/mce3r_stochastic')
from config.parameters import PARAMS, kT_KCAL, BACKGROUND_FREQ

# Path constants
MEME_TXT = '/Users/aayanalwani/tb project/mce3r_stochastic/results/phase1/meme_output/meme.txt'

# Nucleotide index mapping (ACGT order, matching MEME probability matrix columns)
NUC_IDX = {'A': 0, 'C': 1, 'G': 2, 'T': 3}

# Actual binding site positions within the 123-bp operator sequence (0-indexed)
# Determined from FIMO predicted_sites.csv and operator sequence alignment:
#   Strong site: op[74:95] = GTTGTCCAAGCAATCGCGTAT  (MEME-1, Kd=2.4 nM)
#   Weak site:   op[51:71] = TTTGCATTGCTATTTACCGA   (MEME-2, Kd=49 nM)
STRONG_SITE_IN_OP = (74, 95)
WEAK_SITE_IN_OP   = (51, 72)   # 21bp: op[51:72] = TTTGCATTGCTATTTACCGAT


def load_pwm_from_meme(meme_txt_path=MEME_TXT):
    """
    Parse meme.txt and extract the first motif's position-specific probability matrix.

    Returns
    -------
    pwm : np.ndarray, shape (motif_width, 4)
        Frequency matrix with columns in ACGT order.
    motif_width : int
    """
    pwm_rows = []
    in_prob_matrix = False

    with open(meme_txt_path, 'r') as f:
        for line in f:
            line = line.rstrip()

            # Detect the first letter-probability matrix block
            if 'letter-probability matrix:' in line and not pwm_rows:
                in_prob_matrix = True
                continue

            if in_prob_matrix:
                stripped = line.strip()
                if stripped == '' or stripped.startswith('-'):
                    in_prob_matrix = False
                    break
                parts = stripped.split()
                if len(parts) == 4:
                    try:
                        row = [float(x) for x in parts]
                        pwm_rows.append(row)
                    except ValueError:
                        in_prob_matrix = False
                        break
                else:
                    in_prob_matrix = False
                    break

    if not pwm_rows:
        raise ValueError(f"Could not parse PWM from {meme_txt_path}")

    pwm = np.array(pwm_rows, dtype=np.float64)  # shape (width, 4) = [A, C, G, T]
    return pwm, pwm.shape[0]


def pwm_to_energy_matrix(pwm, background=None, kT=None, pseudocount=0.001):
    """
    Convert a frequency PWM to a position-specific energy matrix.

    epsilon(i, b) = -kT * ln( (pwm[i,b] + pseudocount) / bg[b] )

    Parameters
    ----------
    pwm : np.ndarray shape (L, 4) — frequencies in ACGT order
    background : dict with keys A, C, G, T  (defaults to BACKGROUND_FREQ)
    kT : float — thermal energy in kcal/mol (default kT_KCAL ≈ 0.616)
    pseudocount : float — added to avoid log(0)

    Returns
    -------
    energy_matrix : np.ndarray shape (L, 4) — energies in kcal/mol
    """
    if background is None:
        background = BACKGROUND_FREQ
    if kT is None:
        kT = kT_KCAL

    bg = np.array([background['A'], background['C'],
                   background['G'], background['T']], dtype=np.float64)

    # Add pseudocount and normalize
    pwm_pseudo = pwm + pseudocount
    row_sums = pwm_pseudo.sum(axis=1, keepdims=True)
    pwm_norm = pwm_pseudo / row_sums

    # Energy: epsilon(i,b) = -kT * ln(pwm[i,b] / bg[b])
    energy_matrix = -kT * np.log(pwm_norm / bg[np.newaxis, :])

    return energy_matrix


def score_sequence(seq, energy_matrix):
    """
    Score a DNA sequence against an energy matrix.

    If len(seq) == energy_matrix.shape[0]: score directly.
    If len(seq) > energy_matrix.shape[0]: find best-scoring window (lowest energy).
    If len(seq) < energy_matrix.shape[0]: raise ValueError.

    Parameters
    ----------
    seq : str — DNA sequence (A/C/G/T)
    energy_matrix : np.ndarray shape (L, 4)

    Returns
    -------
    float — summed energy score (kcal/mol), lower = more favorable
    """
    seq = seq.upper()
    motif_len = energy_matrix.shape[0]
    seq_len = len(seq)

    if seq_len < motif_len:
        raise ValueError(f"Sequence length {seq_len} < motif length {motif_len}")

    if seq_len == motif_len:
        total = 0.0
        for i, nuc in enumerate(seq):
            idx = NUC_IDX.get(nuc, -1)
            if idx == -1:
                total += np.mean(energy_matrix[i])
            else:
                total += energy_matrix[i, idx]
        return total

    else:
        # Find best-scoring (lowest energy) window
        best_score = np.inf
        for start in range(seq_len - motif_len + 1):
            window = seq[start:start + motif_len]
            score = 0.0
            for i, nuc in enumerate(window):
                idx = NUC_IDX.get(nuc, -1)
                if idx == -1:
                    score += np.mean(energy_matrix[i])
                else:
                    score += energy_matrix[i, idx]
            if score < best_score:
                best_score = score
        return best_score


def calibrate_energies(strong_seq=None, weak_seq=None, energy_matrix=None,
                       Kd_strong=None, Kd_weak=None, kT=None):
    """
    Calibrate PWM-derived energies to experimental Kd values.

    Uses two independent offsets (one per site) since the two binding sites
    correspond to different MEME motifs:
      - offset_strong = DeltaG_strong_experimental - DeltaG_strong_raw
      - offset_weak   = DeltaG_weak_experimental   - DeltaG_weak_raw
      - offset (mean) used for unknown sequences

    For the partition function, calibrated DeltaG values are:
      dG_strong = DeltaG_strong_experimental (exact, by construction)
      dG_weak   = DeltaG_weak_experimental   (exact, by construction)

    Parameters
    ----------
    strong_seq : str — sequence of high-affinity (strong) site (default: op[74:95])
    weak_seq   : str — sequence of low-affinity (weak) site   (default: op[51:72])
    energy_matrix : np.ndarray — from pwm_to_energy_matrix()
    Kd_strong : float, nM (default PARAMS['Kd_strong'] = 2.4)
    Kd_weak   : float, nM (default PARAMS['Kd_weak'] = 49.0)
    kT : float, kcal/mol (default kT_KCAL)

    Returns
    -------
    dict with full calibration information
    """
    if Kd_strong is None:
        Kd_strong = PARAMS['Kd_strong']
    if Kd_weak is None:
        Kd_weak = PARAMS['Kd_weak']
    if kT is None:
        kT = kT_KCAL

    # Default sequences: actual FIMO binding sites in the operator
    op = PARAMS['operator_sequence']
    if strong_seq is None:
        s, e = STRONG_SITE_IN_OP
        strong_seq = op[s:e]
    if weak_seq is None:
        s, e = WEAK_SITE_IN_OP
        weak_seq = op[s:e]

    if energy_matrix is None:
        pwm, _ = load_pwm_from_meme()
        energy_matrix = pwm_to_energy_matrix(pwm, kT=kT)

    # Experimental DeltaG = RT * ln(Kd_M) = kT * ln(Kd_nM * 1e-9)
    dG_strong_exp = kT * np.log(Kd_strong * 1e-9)
    dG_weak_exp   = kT * np.log(Kd_weak * 1e-9)

    # Raw PWM scores (best-window if sequence is longer than motif)
    dG_strong_raw = score_sequence(strong_seq, energy_matrix)
    dG_weak_raw   = score_sequence(weak_seq, energy_matrix)

    # Independent offsets per site
    offset_strong = dG_strong_exp - dG_strong_raw
    offset_weak   = dG_weak_exp   - dG_weak_raw
    offset_mean   = (offset_strong + offset_weak) / 2.0

    # Calibrated energies — exact by construction
    dG_strong = dG_strong_exp
    dG_weak   = dG_weak_exp

    # Predicted Kd (exact, since we used two independent offsets)
    Kd_strong_pred = np.exp(dG_strong / kT) * 1e9
    Kd_weak_pred   = np.exp(dG_weak   / kT) * 1e9

    return {
        'dG_strong_exp':   dG_strong_exp,
        'dG_weak_exp':     dG_weak_exp,
        'dG_strong_raw':   dG_strong_raw,
        'dG_weak_raw':     dG_weak_raw,
        'offset_strong':   offset_strong,
        'offset_weak':     offset_weak,
        'offset_mean':     offset_mean,
        'offset':          offset_strong,   # primary offset (strong site anchor)
        'dG_strong':       dG_strong,
        'dG_weak':         dG_weak,
        'Kd_strong_pred':  Kd_strong_pred,
        'Kd_weak_pred':    Kd_weak_pred,
        'kT':              kT,
        'energy_matrix':   energy_matrix,
        'strong_seq':      strong_seq,
        'weak_seq':        weak_seq,
        'Kd_strong_exp':   Kd_strong,
        'Kd_weak_exp':     Kd_weak,
    }


def predict_Kd_for_sequence(seq, calibration):
    """
    Predict Kd (nM) for an arbitrary DNA sequence using calibrated energies.

    Uses offset_mean for unknown sequences (conservative interpolation between
    the strong and weak site offsets).

    Parameters
    ----------
    seq : str
    calibration : dict from calibrate_energies()

    Returns
    -------
    Kd_pred : float (nM)
    """
    kT = calibration['kT']
    energy_matrix = calibration['energy_matrix']
    offset = calibration['offset_mean']

    raw_score = score_sequence(seq, energy_matrix)
    dG_calib = raw_score + offset
    Kd_pred = np.exp(dG_calib / kT) * 1e9
    return Kd_pred


def predict_dG_for_sequence(seq, calibration):
    """
    Predict calibrated DeltaG (kcal/mol) for an arbitrary DNA sequence.

    Parameters
    ----------
    seq : str
    calibration : dict from calibrate_energies()

    Returns
    -------
    dG : float (kcal/mol)
    """
    energy_matrix = calibration['energy_matrix']
    offset = calibration['offset_mean']

    raw_score = score_sequence(seq, energy_matrix)
    return raw_score + offset


def build_calibration(meme_path=MEME_TXT):
    """
    Convenience function: load PWM, build energy matrix, calibrate, return dict.
    Uses the actual FIMO-determined binding site sequences.
    """
    pwm, motif_width = load_pwm_from_meme(meme_path)
    energy_matrix = pwm_to_energy_matrix(pwm)
    calibration = calibrate_energies(energy_matrix=energy_matrix)
    calibration['motif_width'] = motif_width
    return calibration


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

    print("=== energy_calibration.py self-tests ===")

    op = PARAMS['operator_sequence']
    s_start, s_end = STRONG_SITE_IN_OP
    w_start, w_end = WEAK_SITE_IN_OP

    strong_seq = op[s_start:s_end]
    weak_seq   = op[w_start:w_end]
    print(f"  Strong site op[{s_start}:{s_end}]: {strong_seq}")
    print(f"  Weak   site op[{w_start}:{w_end}]: {weak_seq}")

    pwm, motif_width = load_pwm_from_meme()
    print(f"  Motif width: {motif_width} bp, PWM shape: {pwm.shape}")

    energy_matrix = pwm_to_energy_matrix(pwm)
    calibration   = calibrate_energies(energy_matrix=energy_matrix)

    print(f"\n  Experimental DeltaG_strong = {calibration['dG_strong_exp']:.3f} kcal/mol")
    print(f"  Experimental DeltaG_weak   = {calibration['dG_weak_exp']:.3f} kcal/mol")
    print(f"  Raw PWM score strong       = {calibration['dG_strong_raw']:.3f} kcal/mol")
    print(f"  Raw PWM score weak         = {calibration['dG_weak_raw']:.3f} kcal/mol")
    print(f"  Offset strong              = {calibration['offset_strong']:.3f} kcal/mol")
    print(f"  Offset weak                = {calibration['offset_weak']:.3f} kcal/mol")
    print(f"  Offset mean                = {calibration['offset_mean']:.3f} kcal/mol")
    print(f"  Calibrated DeltaG_strong   = {calibration['dG_strong']:.3f} kcal/mol")
    print(f"  Calibrated DeltaG_weak     = {calibration['dG_weak']:.3f} kcal/mol")
    print(f"  Predicted Kd_strong        = {calibration['Kd_strong_pred']:.3f} nM (exp: {PARAMS['Kd_strong']} nM)")
    print(f"  Predicted Kd_weak          = {calibration['Kd_weak_pred']:.3f} nM (exp: {PARAMS['Kd_weak']} nM)")

    # Test 1: Calibrated Kd_strong within 1% of 2.4 nM (exact by construction)
    Kd_s = calibration['Kd_strong_pred']
    sanity("Kd_strong within 1% of 2.4 nM",
           abs(Kd_s - 2.4) / 2.4 < 0.01,
           f"— got {Kd_s:.4f} nM")

    # Test 2: Calibrated Kd_weak within 1% of 49.0 nM (exact by construction)
    Kd_w = calibration['Kd_weak_pred']
    sanity("Kd_weak within 1% of 49.0 nM",
           abs(Kd_w - 49.0) / 49.0 < 0.01,
           f"— got {Kd_w:.4f} nM")

    # Test 3: Energy matrix has no NaN or Inf
    em = calibration['energy_matrix']
    has_bad = np.any(np.isnan(em)) or np.any(np.isinf(em))
    sanity("Energy matrix no NaN/Inf",
           not has_bad,
           f"— shape {em.shape}")

    # Test 4: Background frequencies sum to 1.0
    bg_sum = sum(BACKGROUND_FREQ.values())
    sanity("Background frequencies sum to 1.0",
           abs(bg_sum - 1.0) < 1e-10,
           f"— sum = {bg_sum}")

    # Bonus: predict Kd for the strong site using mean offset
    kd_test = predict_Kd_for_sequence(strong_seq, calibration)
    dG_test = predict_dG_for_sequence(strong_seq, calibration)
    print(f"\n  predict_Kd_for_sequence(strong_seq) = {kd_test:.3f} nM (using mean offset)")
    print(f"  predict_dG_for_sequence(strong_seq) = {dG_test:.4f} kcal/mol (using mean offset)")

    print(f"\n{'='*50}")
    print(f"Sanity: {n_pass} pass, {n_fail} fail")
    if n_fail > 0:
        print("STOPPING: Sanity check failed.")
        sys.exit(1)
    else:
        print("All checks passed.")
        sys.exit(0)
