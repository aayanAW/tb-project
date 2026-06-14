"""
conservation.py — Leave-one-lineage-out (LOLO) operator conservation (Gate G5).

The audit (C7, S5) showed the old conservation gate was circular and pseudo-replicated:
it scanned the SAME 8 ortholog windows that trained the knowledge PWM, counted 4 near-
clonal M. tuberculosis-complex members as 4 independent species, and used the H37Rv genome
background for non-H37Rv sequence.

This module fixes all three:
  * Leave-one-lineage-out: to test a lineage, a fresh MEME PWM is built from the windows of
    the OTHER lineages only, then FIMO scans the held-out lineage's windows. No lineage is
    ever scored by a model trained on its own sequence.
  * Lineage collapse: the 4 MTBC near-clones count as ONE lineage (mce3r_biology.
    ORTHOLOG_LINEAGE); independent units are MTBC, the M. marinum clade, and M. smegmatis.
  * Mycobacterial background: an order-0 background built from the ortholog windows
    themselves, not the H37Rv whole-genome model.

Detection of an operator in a held-out window = a FIMO site with p < DETECT_P and score > 0.
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

from Bio import SeqIO

sys.path.insert(0, str(Path(__file__).resolve().parent))
import mce3r_biology as bio
from run_fimo import run_fimo
from run_meme import run_meme
from utils import get_project_root, require_meme_tool, setup_logging
from validate import _load_fimo

logger = setup_logging(__name__)

DETECT_P = 1e-3


def _build_ortholog_background(windows_fasta: Path, out_bg: Path) -> Path | None:
    """Order-0 background from the ortholog windows (mycobacterial composition)."""
    fgm = require_meme_tool("fasta-get-markov")
    out_bg.parent.mkdir(parents=True, exist_ok=True)
    try:
        with open(out_bg, "w") as fh:
            r = subprocess.run(
                [fgm, "-m", "0", str(windows_fasta)],
                stdout=fh,
                stderr=subprocess.PIPE,
                text=True,
                timeout=120,
            )
        if r.returncode != 0 or not out_bg.exists() or out_bg.stat().st_size == 0:
            logger.warning(
                f"fasta-get-markov failed (exit {r.returncode}); no LOLO bg."
            )
            return None
        return out_bg
    except Exception as e:  # noqa: BLE001 - background is optional, log and continue
        logger.warning(f"Could not build ortholog background: {e}")
        return None


def leave_one_lineage_out(
    orthologs_fasta: Path,
    output_dir: Path,
    threads: int = 4,
) -> dict:
    """
    Run LOLO conservation across mycobacterial lineages. Returns a G5 gate dict.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    recs = list(SeqIO.parse(str(orthologs_fasta), "fasta"))
    if not recs:
        return {
            "name": "conservation_leave_one_lineage_out",
            "error": f"no ortholog windows in {orthologs_fasta}",
            "pass": False,
        }

    by_lineage: dict[str, list] = defaultdict(list)
    for r in recs:
        by_lineage[bio.ortholog_lineage(r.id)].append(r)
    lineages = sorted(by_lineage)

    bg = _build_ortholog_background(orthologs_fasta, output_dir / "ortholog_order0.bg")

    w = bio.OPERATOR_SITE_WIDTH
    per_lineage = {}
    for test_lin in lineages:
        train_recs = [r for lin in lineages if lin != test_lin for r in by_lineage[lin]]
        test_recs = by_lineage[test_lin]
        if len(train_recs) < 2:
            per_lineage[test_lin] = {
                "detected": False,
                "skipped": f"only {len(train_recs)} training windows (need >=2)",
                "n_test_windows": len(test_recs),
            }
            continue
        lin_dir = output_dir / f"lolo_{test_lin}"
        lin_dir.mkdir(parents=True, exist_ok=True)
        train_fa = lin_dir / "train.fasta"
        test_fa = lin_dir / "test.fasta"
        SeqIO.write(train_recs, str(train_fa), "fasta")
        SeqIO.write(test_recs, str(test_fa), "fasta")

        meme_txt = run_meme(
            train_fa,
            lin_dir / "meme",
            nmotifs=1,
            minw=max(8, w - 4),
            maxw=w + 4,
            mod="anr",
            threads=threads,
            bfile=bg,
        )
        pwm = lin_dir / "lolo.meme"
        pwm.write_text(Path(meme_txt).read_text())

        fimo_tsv = run_fimo(
            pwm,
            test_fa,
            lin_dir / "fimo",
            bfile=bg,
            qv_thresh=DETECT_P,
            use_qvalue=False,
        )
        df = _load_fimo(fimo_tsv)
        detected_ids = set()
        positional_ids = set()
        if not df.empty and "score" in df.columns and "p-value" in df.columns:
            # Gate criterion (prereg-faithful): positive log-odds AND p < DETECT_P.
            detected_ids = set(
                df[(df["score"] > 0) & (df["p-value"] < DETECT_P)][
                    "sequence_name"
                ].astype(str)
            )
            # Sensitivity readout (NOT the gate): positionally significant regardless of
            # log-odds sign. A diverged ortholog can match the operator position at low p
            # yet score below the GC-rich background; reporting this keeps the nuance
            # visible without relaxing the gate to manufacture a pass.
            positional_ids = set(
                df[df["p-value"] < DETECT_P]["sequence_name"].astype(str)
            )
        per_lineage[test_lin] = {
            "detected": len(detected_ids) > 0,
            "n_test_windows": len(test_recs),
            "n_train_windows": len(train_recs),
            "detected_windows": sorted(detected_ids),
            "positionally_significant_windows": sorted(positional_ids),
        }
        logger.info(
            f"LOLO {test_lin}: trained on {len(train_recs)} windows, "
            f"detected={per_lineage[test_lin]['detected']}"
        )

    n_detected = sum(1 for v in per_lineage.values() if v.get("detected"))
    n_positional = sum(
        1 for v in per_lineage.values() if v.get("positionally_significant_windows")
    )
    return {
        "name": "conservation_leave_one_lineage_out",
        "lineages": lineages,
        "n_lineages": len(lineages),
        "n_lineages_detected": n_detected,
        "n_lineages_positionally_significant": n_positional,
        "min_required_lineages": bio.CONSERVATION_MIN_LINEAGES,
        "detection_p_threshold": DETECT_P,
        "per_lineage": per_lineage,
        "note": (
            "Each lineage is scored by a PWM trained on the OTHER lineages only "
            "(non-circular). MTBC near-clones collapsed to one lineage. 'detected' is the "
            "gate (score>0 AND p<thresh); 'positionally_significant' is a sensitivity "
            "readout (p<thresh, any sign) that is NOT used for pass/fail."
        ),
        "pass": n_detected >= bio.CONSERVATION_MIN_LINEAGES,
    }


def main():
    root = get_project_root()
    parser = argparse.ArgumentParser(
        description="Leave-one-lineage-out conservation (G5)"
    )
    parser.add_argument(
        "--orthologs",
        type=Path,
        default=root / "data" / "raw" / "ortholog_promoters.fasta",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=root / "results" / "scans" / "_conservation_lolo",
    )
    parser.add_argument("--threads", type=int, default=4)
    args = parser.parse_args()
    g5 = leave_one_lineage_out(args.orthologs, args.output_dir, args.threads)
    import json

    print(json.dumps(g5, indent=2, default=str))


if __name__ == "__main__":
    main()
