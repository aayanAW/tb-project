"""
validate.py — Non-circular, paper-grade validation of the genome-wide Mce3R operator scan.

This is a rewrite that addresses the cross-model audit (audit_report.md). The previous
gates were circular/leaky: they validated the knowledge PWM against the very regions it was
built from, and the q-values were not a genuine genome-wide BH FDR. The corrected gates:

  G2 held-out site recovery — the Rv1935c-Rv1936 operator region is NEVER used to build the
                              knowledge PWM (it is trained on the mce3R-yrbE3A window + the
                              ortholog windows only), so recovering it is a true generalization
                              test, not re-recovery of training data. The mce3R-yrbE3A region
                              is reported as a positive control, not counted toward the gate.
  G3 specificity            — genome-wide promoter hit fraction <= 5% at a real BH q<0.05.
  G4 negative-control       — mce1/mce2/mce4 promoters (Mce3R does NOT bind them; Santangelo
     specificity             2008) must behave like background: no strict hits, not in top 5%.

G1 (external corroboration vs a published motif) and G5 (leave-one-lineage-out conservation)
are computed in main.py, where the MEME/Tomtom tooling lives, and merged into the gate set.

FDR: q-values are recomputed here as a transparent, pooled, genome-wide Benjamini-Hochberg
correction over ALL candidate positions (both motifs, both strands), not FIMO's per-report
estimate. Because FIMO is run with a permissive p<1e-3 reporting threshold, every site that
could possibly reach q<0.05 is present in the report (all reported p-values are strictly
smaller than every un-reported one), so the recomputed BH is exact for the significant set.
"""

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import SeqIO

sys.path.insert(0, str(Path(__file__).resolve().parent))
import mce3r_biology as bio
from utils import get_project_root, setup_logging

logger = setup_logging(__name__)

Q_THRESH = 0.05
SPECIFICITY_CEILING = 0.05
TOP_FRACTION = 0.05
N_STRANDS = 2  # FIMO scans both strands; each is an independent test.


def _load_fimo(fimo_tsv: Path) -> pd.DataFrame:
    fimo_tsv = Path(fimo_tsv)
    cols = [
        "sequence_name",
        "start",
        "stop",
        "strand",
        "score",
        "p-value",
        "q-value",
    ]
    if not fimo_tsv.exists():
        return pd.DataFrame(columns=cols)
    try:
        df = pd.read_csv(fimo_tsv, sep="\t", comment="#")
    except pd.errors.EmptyDataError:
        # FIMO writes an empty/header-only file when a scan yields no hits (e.g. a
        # diverged ortholog with no operator) — a valid "no detection" result.
        return pd.DataFrame(columns=cols)
    df.columns = [c.strip() for c in df.columns]
    for c in ("score", "p-value", "q-value"):
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce")
    return df.dropna(subset=["sequence_name"])


def _motif_widths(motif_file: Path | None) -> list[int]:
    """Return the width of every motif in a MEME file (for the genome-wide test count)."""
    if motif_file is None or not Path(motif_file).exists():
        return []
    from run_meme import extract_evalues_from_meme

    return [
        m["width"] for m in extract_evalues_from_meme(Path(motif_file)) if m["width"]
    ]


def genome_wide_test_count(seq_lengths: list[int], motif_widths: list[int]) -> int:
    """
    Total number of independent FIMO tests across the scanned universe:
        sum over motifs, sum over sequences, of max(0, L - w + 1), times both strands.

    This is the denominator for the pooled genome-wide BH FDR.
    """
    total = 0
    for w in motif_widths:
        for L in seq_lengths:
            if L >= w:
                total += L - w + 1
    return total * N_STRANDS


def benjamini_hochberg_q(pvalues: np.ndarray, n_tests: int) -> np.ndarray:
    """
    Genome-wide BH q-values for a set of reported p-values drawn from n_tests total tests.

    Valid when the reported p-values are exactly the smallest of all n_tests p-values
    (guaranteed here: FIMO reports every site below the p-threshold, so all un-reported
    p-values are larger). q_(k) = min_{j>=k} ( p_(j) * n_tests / j ), capped at 1.
    """
    p = np.asarray(pvalues, dtype=float)
    if p.size == 0:
        return p
    order = np.argsort(p, kind="mergesort")
    ranked = p[order]
    ranks = np.arange(1, p.size + 1)
    raw = ranked * n_tests / ranks
    # Enforce monotonicity from the largest rank downward.
    q_sorted = np.minimum.accumulate(raw[::-1])[::-1]
    q_sorted = np.clip(q_sorted, 0.0, 1.0)
    q = np.empty_like(q_sorted)
    q[order] = q_sorted
    return q


def _per_sequence_stats(
    hits: pd.DataFrame, universe: list[str], n_tests: int, q_thresh: float
) -> dict:
    """
    Collapse FIMO rows to per-sequence stats with a SAME-ROW strict-hit definition.

    A sequence is a strict hit iff it has at least one site whose score > 0 AND whose
    genome-wide BH q-value < q_thresh on the SAME row (audit C5 fix).
    """
    best_score = {sid: float("-inf") for sid in universe}
    best_q = {sid: 1.0 for sid in universe}
    strict = set()

    if not hits.empty and "p-value" in hits.columns:
        bh_q = benjamini_hochberg_q(hits["p-value"].to_numpy(), n_tests)
        hits = hits.assign(bh_q=bh_q)
    else:
        hits = hits.assign(bh_q=pd.Series(dtype=float))

    for _, row in hits.iterrows():
        sid = str(row["sequence_name"])
        if sid not in best_score:
            best_score[sid] = float("-inf")
            best_q[sid] = 1.0
        sc = float(row.get("score", float("-inf")))
        rq = float(row.get("bh_q", 1.0))
        if sc > best_score[sid]:
            best_score[sid] = sc
        if np.isfinite(rq):
            best_q[sid] = min(best_q[sid], rq)
        # SAME-ROW strict hit.
        if sc > 0 and np.isfinite(rq) and rq < q_thresh:
            strict.add(sid)
    return {"best_score": best_score, "best_q": best_q, "strict_hits": strict}


def evaluate(
    promoters_fasta: Path,
    genome_fimo_tsv: Path,
    motif_file: Path | None = None,
    q_thresh: float = Q_THRESH,
) -> dict:
    """Compute gates G2, G3, G4 with a genome-wide BH FDR. G1/G5 are merged by main.py."""
    records = list(SeqIO.parse(str(promoters_fasta), "fasta"))
    universe = [r.id for r in records]
    seq_lengths = [len(r.seq) for r in records]
    n_total = len(universe)
    labels_map = {sid: bio.classify_locus(sid) for sid in universe}

    motif_widths = _motif_widths(motif_file)
    n_tests = genome_wide_test_count(seq_lengths, motif_widths)
    if n_tests == 0:
        logger.warning(
            "Genome-wide test count is 0 (no motif widths?). FDR will be unavailable."
        )

    hits = _load_fimo(genome_fimo_tsv)
    stats = _per_sequence_stats(hits, universe, n_tests, q_thresh)
    best_score, best_q, strict_hits = (
        stats["best_score"],
        stats["best_q"],
        stats["strict_hits"],
    )

    # Genome-wide ranking score per sequence (floor for never-scored sequences).
    finite = [s for s in best_score.values() if np.isfinite(s)]
    floor = (min(finite) - 1.0) if finite else 0.0
    rank_score = np.array(
        [best_score[s] if np.isfinite(best_score[s]) else floor for s in universe]
    )
    # Top 5% cutoff over the FULL universe (audit C10 fix: was over scored-only).
    top_cut = (
        float(np.quantile(rank_score, 1 - TOP_FRACTION)) if n_total else float("inf")
    )

    # ── G3 specificity ──────────────────────────────────────────────────────────
    n_hit = len([s for s in universe if s in strict_hits])
    hit_fraction = n_hit / max(n_total, 1)
    g3 = {
        "name": "specificity",
        "hit_fraction": round(hit_fraction, 5),
        "n_hit_promoters": n_hit,
        "n_total_promoters": n_total,
        "n_genome_wide_tests": n_tests,
        "ceiling": SPECIFICITY_CEILING,
        "fdr_q_threshold": q_thresh,
        "pass": hit_fraction <= SPECIFICITY_CEILING,
    }

    # ── G2 held-out site recovery ─────────────────────────────────────────────────
    def recovery(sid: str) -> dict:
        present = sid in best_score and np.isfinite(best_score[sid])
        return {
            "scanned": bool(present),
            "best_score": round(best_score[sid], 3) if present else None,
            "best_bh_q": round(best_q[sid], 8) if present else None,
            "strict_hit": sid in strict_hits,
            "in_top5pct": bool(present and best_score[sid] >= top_cut),
        }

    heldout = bio.HELDOUT_OPERATOR_IGR
    training = bio.TRAINING_OPERATOR_IGR
    heldout_rec = recovery(heldout)
    training_rec = recovery(training)
    g2_pass = bool(heldout_rec["strict_hit"] and heldout_rec["in_top5pct"])
    g2 = {
        "name": "heldout_site_recovery",
        "heldout_operator_igr": heldout,
        "heldout_recovery": heldout_rec,
        "training_operator_igr": training,
        "training_recovery_positive_control": training_rec,
        "top5pct_score_cutoff": round(top_cut, 3) if np.isfinite(top_cut) else None,
        "note": (
            "Pass depends ONLY on the held-out Rv1935c-Rv1936 operator, which is not used "
            "to build the knowledge PWM. The mce3R-yrbE3A region is a training positive "
            "control (expected to recover; not counted)."
        ),
        "pass": g2_pass,
    }

    # ── G4 negative-control specificity ───────────────────────────────────────────
    neg = [s for s in universe if labels_map[s] == "negative_control"]
    neg_strict = [s for s in neg if s in strict_hits]
    neg_top = [
        s for s in neg if np.isfinite(best_score[s]) and best_score[s] >= top_cut
    ]
    neg_rate = len(neg_strict) / len(neg) if neg else None
    n_declared = len(bio.NEGATIVE_CONTROL_GENES)
    g4 = {
        "name": "negative_control_specificity",
        # Honest coverage: co-transcribed interior genes of mce1/2/4 have no own promoter,
        # so only a subset of declared controls is testable. Reported, not hidden (C-1).
        "n_negative_controls_declared": n_declared,
        "n_negative_controls_tested": len(neg),
        "negative_control_coverage_note": (
            f"{len(neg)} of {n_declared} declared mce1/2/4 genes have a promoter in the "
            "universe and are testable; the rest are co-transcribed interior operon genes."
        ),
        "negative_controls_tested": sorted(neg),
        "negative_control_strict_hits": sorted(neg_strict),
        "negative_control_in_top5pct": sorted(neg_top),
        "negative_control_hit_rate": round(neg_rate, 4)
        if neg_rate is not None
        else None,
        "background_hit_rate": round(hit_fraction, 5),
        "note": (
            "Mce3R does not regulate mce1/mce2/mce4 (Santangelo 2008). A specific operator "
            "model must leave them at background: zero strict hits and none in the top 5%."
        ),
        # Descriptive only (NOT a pass criterion): the leaky regulon AUPRC the audit flagged.
        "descriptive_regulon_auprc": _descriptive_regulon_auprc(
            universe, labels_map, rank_score
        ),
        "pass": bool(neg and len(neg_strict) == 0 and len(neg_top) == 0),
    }

    gates = {"G2": g2, "G3": g3, "G4": g4}
    return gates


def _descriptive_regulon_auprc(universe, labels_map, rank_score) -> dict:
    """
    The old G4 metric, kept ONLY as a descriptive readout with an explicit caveat.

    It is circular as a gate (most positives are the PWM's own training regions), so it is
    never used for pass/fail — reported so the regression vs the old pipeline is visible.
    """
    y = np.array(
        [1 if labels_map[s] in ("operator", "regulon") else 0 for s in universe]
    )
    if y.sum() == 0 or y.sum() == len(y):
        return {"auprc": None, "auroc": None, "caveat": "degenerate label set"}
    try:
        from sklearn.metrics import average_precision_score, roc_auc_score

        ap = float(average_precision_score(y, rank_score))
        roc = float(roc_auc_score(y, rank_score))
    except Exception:
        ap, roc = None, None
    return {
        "auprc": round(ap, 4) if ap is not None else None,
        "auroc": round(roc, 4) if roc is not None else None,
        "n_positives": int(y.sum()),
        "caveat": (
            "CIRCULAR — most positives are the operator-flanking regions used to train the "
            "PWM. Descriptive only; not a validation gate."
        ),
    }


def write_report(gates: dict, output_path: Path) -> None:
    order = ["G1", "G2", "G3", "G4", "G5"]
    lines = [
        "=" * 72,
        "  MCE3R OPERATOR SCAN — VALIDATION REPORT (non-circular)",
        "=" * 72,
        "",
    ]
    for key in order:
        g = gates.get(key)
        if g is None:
            continue
        status = "PASS" if g.get("pass") else "FAIL"
        lines.append(f"[{key}] {g.get('name', ''):<28} {status}")
        for k, v in g.items():
            if k in ("name", "pass"):
                continue
            lines.append(f"      {k}: {v}")
        lines.append("")
    if "all_gates_pass" in gates:
        lines.append(f"ALL GATES PASS (G1-G5): {gates['all_gates_pass']}")
    elif "G2_G3_G4_pass" in gates:
        lines.append(
            f"G2-G4 PASS (partial; run main.py for G1+G5): {gates['G2_G3_G4_pass']}"
        )
    lines.append("=" * 72)
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    Path(output_path).write_text("\n".join(lines) + "\n")
    logger.info(f"Validation report -> {output_path}")


def main():
    root = get_project_root()
    parser = argparse.ArgumentParser(description="Validate the Mce3R genome-wide scan")
    parser.add_argument(
        "--promoters", type=Path, default=root / "data" / "raw" / "promoters.fasta"
    )
    parser.add_argument(
        "--fimo", type=Path, default=root / "results" / "scans" / "fimo.tsv"
    )
    parser.add_argument(
        "--motif-file",
        type=Path,
        default=root
        / "results"
        / "motifs"
        / "operator_knowledge"
        / "operator_knowledge.meme",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=root / "results" / "scans" / "validation_report.txt",
    )
    parser.add_argument(
        "--json-out", type=Path, default=root / "results" / "scans" / "gates.json"
    )
    args = parser.parse_args()

    gates = evaluate(args.promoters, args.fimo, args.motif_file)
    # Standalone validate computes ONLY G2/G3/G4 (G1 external + G5 conservation need the
    # MEME/Tomtom steps in main.py). Use a scoped key so this can never be mistaken for the
    # 5-gate verdict that main.py writes (audit self-check H-2).
    gates["G2_G3_G4_pass"] = bool(
        gates["G2"]["pass"] and gates["G3"]["pass"] and gates["G4"]["pass"]
    )
    gates["note"] = (
        "Partial: run main.py for the full 5-gate verdict (G1 + G5 included)."
    )
    write_report(gates, args.output)
    Path(args.json_out).write_text(json.dumps(gates, indent=2, default=str))
    print(Path(args.output).read_text())


if __name__ == "__main__":
    main()
