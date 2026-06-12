"""
validate.py — Paper-grade validation of the genome-wide Mce3R operator scan.

Computes the pre-registered gates (PREREGISTRATION.md) against the known ground truth in
mce3r_biology.py. This is the part the old pipeline never did: it asks whether the scan
behaves like a specific repressor or like a GC detector.

Gates:
  G2 known sites recovered  — both mapped operator regions have a hit at q<0.05 & score>0,
                              in the top 5% of genome-wide scores.
  G3 specificity            — genome-wide promoter hit fraction <= 5%.
  G4 regulon enrichment     — AUPRC > 0.5 and empirical p < 0.05 (random-set + hypergeometric);
                              negative controls (mce1/2/4) behave like background.
  G5 conservation           — operator detected in >= 8 mycobacterial orthologs.

No sklearn dependency required (AUPRC/AUROC implemented locally; uses sklearn if present).
"""

import argparse
import json
import sys
from math import comb
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
N_RANDOM_NULL = 1000
NULL_SEED = 7
CONSERVATION_MIN_ORTHOLOGS = 8


def _load_fimo(fimo_tsv: Path) -> pd.DataFrame:
    fimo_tsv = Path(fimo_tsv)
    if not fimo_tsv.exists():
        return pd.DataFrame(
            columns=[
                "sequence_name",
                "start",
                "stop",
                "strand",
                "score",
                "p-value",
                "q-value",
            ]
        )
    df = pd.read_csv(fimo_tsv, sep="\t", comment="#")
    df.columns = [c.strip() for c in df.columns]
    for c in ("score", "p-value", "q-value"):
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce")
    return df.dropna(subset=["sequence_name"])


def _auprc_auroc(scores: np.ndarray, labels: np.ndarray) -> tuple[float, float]:
    """Average precision + AUROC. Uses sklearn if available, else a local implementation."""
    try:
        from sklearn.metrics import average_precision_score, roc_auc_score

        return float(average_precision_score(labels, scores)), float(
            roc_auc_score(labels, scores)
        )
    except Exception:
        pass
    # Local fallback.
    order = np.argsort(-scores, kind="mergesort")
    y = labels[order]
    P = y.sum()
    N = len(y) - P
    if P == 0 or N == 0:
        return float("nan"), float("nan")
    tp = np.cumsum(y)
    fp = np.cumsum(1 - y)
    precision = tp / (tp + fp)
    recall = tp / P
    # Average precision = sum over thresholds of (recall change * precision).
    rec_prev = np.concatenate([[0.0], recall[:-1]])
    ap = float(np.sum((recall - rec_prev) * precision))
    # AUROC via rank statistic (Mann-Whitney).
    ranks = np.argsort(np.argsort(scores)) + 1
    auroc = (ranks[labels == 1].sum() - P * (P + 1) / 2) / (P * N)
    return ap, float(auroc)


def evaluate(
    promoters_fasta: Path,
    genome_fimo_tsv: Path,
    ortholog_fimo_tsv: Path | None = None,
    q_thresh: float = Q_THRESH,
) -> dict:
    """Compute all gate metrics. Returns a dict with per-gate results + summary."""
    universe = [r.id for r in SeqIO.parse(str(promoters_fasta), "fasta")]
    n_total = len(universe)
    labels_map = {sid: bio.classify_locus(sid) for sid in universe}

    hits = _load_fimo(genome_fimo_tsv)
    has_q = "q-value" in hits.columns

    # Per-sequence best score and best q over reported FIMO rows.
    best_score = {sid: float("-inf") for sid in universe}
    best_q = {sid: 1.0 for sid in universe}
    for _, row in hits.iterrows():
        sid = str(row["sequence_name"])
        if sid not in best_score:
            best_score[sid] = float("-inf")
            best_q[sid] = 1.0
            labels_map[sid] = bio.classify_locus(sid)
        sc = float(row.get("score", float("-inf")))
        if sc > best_score[sid]:
            best_score[sid] = sc
        if has_q and pd.notna(row.get("q-value")):
            best_q[sid] = min(best_q[sid], float(row["q-value"]))

    # Strict hit = q < thresh AND score > 0.
    def is_strict_hit(sid: str) -> bool:
        return best_score[sid] > 0 and (best_q[sid] < q_thresh if has_q else True)

    strict_hits = {sid for sid in best_score if is_strict_hit(sid)}

    # Ranking score per sequence (floor for never-scored sequences).
    finite = [s for s in best_score.values() if np.isfinite(s)]
    floor = (min(finite) - 1.0) if finite else 0.0
    rank_score = np.array(
        [best_score[s] if np.isfinite(best_score[s]) else floor for s in universe]
    )

    # ── G3 specificity ──────────────────────────────────────────────────────────
    hit_fraction = len([s for s in universe if s in strict_hits]) / max(n_total, 1)
    g3 = {
        "name": "specificity",
        "hit_fraction": round(hit_fraction, 4),
        "n_hit_promoters": len([s for s in universe if s in strict_hits]),
        "n_total_promoters": n_total,
        "ceiling": SPECIFICITY_CEILING,
        "pass": hit_fraction <= SPECIFICITY_CEILING,
    }

    # ── G2 known-site recovery ────────────────────────────────────────────────────
    # Top 5% score cutoff among scored sequences.
    scored = np.array([best_score[s] for s in universe if np.isfinite(best_score[s])])
    top_cut = np.quantile(scored, 1 - TOP_FRACTION) if len(scored) else float("inf")
    operator_recovery = {}
    for sid in sorted(bio.KNOWN_OPERATOR_SEQUENCES):
        present = sid in best_score and np.isfinite(best_score[sid])
        operator_recovery[sid] = {
            "scanned": bool(present),
            "best_score": round(best_score[sid], 3) if present else None,
            "best_q": round(best_q[sid], 6) if present else None,
            "strict_hit": sid in strict_hits,
            "in_top5pct": bool(present and best_score[sid] >= top_cut),
        }
    # Require both mapped operator IGRs to be recovered (strict hit + top 5%).
    igr_ok = all(
        operator_recovery.get(igr, {}).get("strict_hit")
        and operator_recovery.get(igr, {}).get("in_top5pct")
        for igr in bio.OPERATOR_IGRS
    )
    g2 = {
        "name": "known_site_recovery",
        "operator_igrs": list(bio.OPERATOR_IGRS),
        "recovery": operator_recovery,
        "top5pct_score_cutoff": round(float(top_cut), 3)
        if np.isfinite(top_cut)
        else None,
        "pass": bool(igr_ok),
    }

    # ── G4 regulon enrichment ─────────────────────────────────────────────────────
    y = np.array(
        [1 if labels_map[s] in ("operator", "regulon") else 0 for s in universe]
    )
    auprc, auroc = _auprc_auroc(rank_score, y)
    baseline_auprc = y.sum() / len(y)

    n_pos = int(y.sum())
    obs_pos_hits = sum(
        1
        for s in universe
        if labels_map[s] in ("operator", "regulon") and s in strict_hits
    )
    total_hits = len([s for s in universe if s in strict_hits])

    # Random-set empirical null: sample n_pos promoters at random, count hits.
    rng = np.random.default_rng(NULL_SEED)
    hit_mask = np.array([1 if s in strict_hits else 0 for s in universe])
    null_counts = np.array(
        [
            hit_mask[rng.choice(n_total, size=n_pos, replace=False)].sum()
            for _ in range(N_RANDOM_NULL)
        ]
    )
    emp_p = float((np.sum(null_counts >= obs_pos_hits) + 1) / (N_RANDOM_NULL + 1))

    # Hypergeometric p (analytic): P(X >= obs) drawing n_pos from n_total with total_hits successes.
    def hypergeom_sf(k, M, n, N):
        # P(X >= k) for X~Hypergeometric(M population, n successes, N draws).
        if N == 0 or n == 0:
            return 1.0
        denom = comb(M, N)
        total = 0
        upper = min(n, N)
        for x in range(k, upper + 1):
            total += comb(n, x) * comb(M - n, N - x)
        return total / denom if denom else 1.0

    try:
        hg_p = float(hypergeom_sf(obs_pos_hits, n_total, total_hits, n_pos))
    except (ValueError, OverflowError):
        hg_p = float("nan")

    # Negative-control hit rate (mce1/2/4 must look like background).
    neg = [s for s in universe if labels_map[s] == "negative_control"]
    neg_hits = sum(1 for s in neg if s in strict_hits)
    neg_rate = neg_hits / len(neg) if neg else None

    g4 = {
        "name": "regulon_enrichment",
        "auprc": round(auprc, 4) if np.isfinite(auprc) else None,
        "auprc_baseline": round(float(baseline_auprc), 4),
        "auroc": round(auroc, 4) if np.isfinite(auroc) else None,
        "n_positives": n_pos,
        "observed_positive_hits": obs_pos_hits,
        "total_strict_hits": total_hits,
        "empirical_p_random_set": round(emp_p, 5),
        "hypergeometric_p": round(hg_p, 6) if np.isfinite(hg_p) else None,
        "negative_control_hit_rate": round(neg_rate, 4)
        if neg_rate is not None
        else None,
        "background_hit_rate": round(hit_fraction, 4),
        "pass": bool(
            np.isfinite(auprc)
            and auprc > 0.5
            and emp_p < 0.05
            and (neg_rate is None or neg_rate <= hit_fraction + 1e-9)
        ),
    }

    # ── G5 conservation ───────────────────────────────────────────────────────────
    g5 = {
        "name": "conservation",
        "n_orthologs_with_operator": 0,
        "min_required": CONSERVATION_MIN_ORTHOLOGS,
        "pass": False,
    }
    if ortholog_fimo_tsv is not None and Path(ortholog_fimo_tsv).exists():
        odf = _load_fimo(ortholog_fimo_tsv)
        if not odf.empty:
            pos = odf[odf["score"] > 0] if "score" in odf.columns else odf
            n_orth = pos["sequence_name"].nunique()
            g5["n_orthologs_with_operator"] = int(n_orth)
            g5["pass"] = n_orth >= CONSERVATION_MIN_ORTHOLOGS
            g5["underpowered"] = n_orth < CONSERVATION_MIN_ORTHOLOGS

    gates = {"G2": g2, "G3": g3, "G4": g4, "G5": g5}
    gates["all_core_pass"] = bool(
        g2["pass"] and g3["pass"] and g4["pass"]
    )  # G2-G4 don't need extra orthologs
    return gates


def write_report(gates: dict, output_path: Path) -> None:
    lines = ["=" * 72, "  MCE3R OPERATOR SCAN — VALIDATION REPORT", "=" * 72, ""]
    for key in ("G2", "G3", "G4", "G5"):
        g = gates[key]
        status = "PASS" if g.get("pass") else "FAIL"
        lines.append(f"[{key}] {g['name']:<22} {status}")
        for k, v in g.items():
            if k in ("name", "pass", "recovery"):
                continue
            lines.append(f"      {k}: {v}")
        if key == "G2":
            for sid, rec in g["recovery"].items():
                if sid in bio.OPERATOR_IGRS:
                    lines.append(f"      operator {sid}: {rec}")
        lines.append("")
    lines.append(f"CORE GATES (G2-G4) PASS: {gates['all_core_pass']}")
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
    parser.add_argument("--ortholog-fimo", type=Path, default=None)
    parser.add_argument(
        "--output",
        type=Path,
        default=root / "results" / "scans" / "validation_report.txt",
    )
    parser.add_argument(
        "--json-out", type=Path, default=root / "results" / "scans" / "gates.json"
    )
    args = parser.parse_args()

    gates = evaluate(args.promoters, args.fimo, args.ortholog_fimo)
    write_report(gates, args.output)
    Path(args.json_out).write_text(json.dumps(gates, indent=2, default=str))
    print(Path(args.output).read_text())


if __name__ == "__main__":
    main()
