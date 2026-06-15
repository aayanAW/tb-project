# Cross-Model Audit — Bipartite Operator Discovery (find_operators_v3.py)

**Date:** 2026-06-15
**Scope:** the Phase-1 bipartite operator-discovery method that claims to find Mce3R's 3 unmapped
operators — `find_operators_v3.py`, `phase1_pipeline/cross_validate.py`, `results/METHODOLOGY_DOCUMENT.md`,
`results/phase1/bipartite_operators_v3.csv` (496 candidates).
**Models:** A = Claude (4 parallel finders: statistics, circularity, code, biology). B = Codex / GPT-5.5
(independent + adversarial), run blind to A. **Read-only — no code modified.**

## One-line verdict

**The known operator is recovered as a circular positive control, but there is NO statistically
defensible evidence for novel operators.** The headline Fisher p-values are an empirical-floor
artifact, there is no multiple-testing/FDR control over the ~10⁸ pairs tested (expected chance
bipartite pairs ≈ hundreds — same order as the 496 reported), the search is a single-template
similarity scan whose every "validation" layer rewards resemblance to the one known operator, and
the architecture filter doesn't enforce the real weak→strong arrangement. Of the 4 claimed operators:
#1 is the training template (control), #2 was already footprinted by Santangelo 2009 (re-detection,
not novel), #3 sits in the same IGR as #1 (expected multiplicity / likely not distinct), #4 (Rv1115)
is in the undifferentiated false-positive tail. **"Find the remaining 3 operators" does not hold as written.**

## CONFIRMED by BOTH model families (highest confidence)

| #   | Sev      | Finding                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               | Location                                                              |
| --- | -------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | --------------------------------------------------------------------- |
| B1  | CRITICAL | **No FDR / multiple-testing correction** over the 496 candidates or the ~10⁸ position×strand×spacer pairs evaluated. A correct expected-count gives **≈120–500 chance bipartite pairs genome-wide** — the same order as the 496 reported, i.e. the candidate list ≈ the null expectation. The methodology's "~0.07" formula is dimensionally incoherent and its "1e-16" is a per-locus probability quoted as a genome-wide rate (off by ~10⁸).                                                                                        | `find_operators_v3.py:344-406`; `METHODOLOGY_DOCUMENT.md:412,515,670` |
| B2  | CRITICAL | **Headline Fisher p = 1.09e-10 is a floor artifact.** Empirical p floor = 1/(500000+1) ≈ **2.00e-6**; Fisher-combining two floored inputs deterministically yields 1.08978e-10 (both A and B computed this). Rank 1 (the real operator, score 57.4) and rank 2 (score 34.6) get the **identical** p — the p-value has zero resolving power at the top; ranking is actually driven by the uncalibrated `adjusted_score`. All sub-2e-6 p-values are fabricated precision.                                                               | `find_operators_v3.py:183-216,371-380`; CSV ranks 1-2                 |
| B3  | CRITICAL | **Circular / single-template search.** PWMs are built from the ONE known operator (n=1 sequence, 0.5 pseudocount → fuzzy exact-match template). The scan can only find look-alikes; it is structurally blind to the _diverged_ operators Panagoda anticipates ("different DNA contacts"). Recovering the known operator at rank #1 is tautological (the PWM _is_ that operator), not validation of generalization.                                                                                                                    | `find_operators_v3.py:30-34,119-153,532-549`                          |
| B4  | HIGH     | **Proximity prior gates the confidence tiers.** `HIGHEST` is mathematically impossible without `near_regulon` (`cross_validate.py:144-157`), so flagship "discoveries" #2/#3 are HIGHEST partly _because_ they sit next to known regulon genes. The "3 independent methods" (bipartite + FIMO v1 + FIMO v2) are all seeded from the same known-operator/yrbE3A locus → shared provenance, not independence. The "core-identity" boost re-uses the known operator's own 15bp recognition helix → double-counting the same measurement. | `find_operators_v3.py:478-485`; `cross_validate.py:106-157`           |
| B5  | HIGH     | **Architecture not enforced.** Pairing loops over `strong_hits + weak_hits` and only checks same-strand + spacer; it counts (strong,strong), (weak,weak), and reversed-orientation pairs → ~4× candidate inflation and breaks the FP math. The "swapped architecture" highlighted for operator #2 is an artifact of not requiring weak-upstream/strong-downstream, not a discovered feature.                                                                                                                                          | `find_operators_v3.py:348-401`                                        |
| B6  | HIGH     | **Novelty inflation / biology.** #2 (Rv1935c/Rv1936) was already DNase-footprinted by Santangelo 2009 → re-detection, not novel discovery. #3 is 446 bp from #1 in the SAME intergenic region (doc says "~440bp" — coordinate inconsistency) → expected motif multiplicity / probable sub-feature of one extended/looping operator, not a distinct 4th site. #4 (Rv1115) is out-of-regulon, no cholesterol/lipid/stress link, and sits in a smooth score continuum with the FP tail (#4=25.9, #5=23.8, #6=23.2 … no gap).             | `METHODOLOGY_DOCUMENT.md:486-591`; CSV ranks 2-6                      |
| B7  | HIGH     | **Fisher independence + selection bias.** The two half-site p-values are not independent (adjacent windows correlate in a 65.6%-GC genome; tandem-repeat self-pairs appear at CSV ranks 93-96 combining a sequence with a near-copy of itself), and the pair is _selected_ for passing p<5e-4 at ~53bp then tested on that same selection → anti-conservative.                                                                                                                                                                        | `find_operators_v3.py:284,293,371-373`; CSV ranks 93-96               |

## SINGLE-MODEL / lower severity (real, but didn't change the top-4)

- **MEDIUM (code):** minus-strand pairs use plus-strand coordinates → wrong biological orientation and plus/minus double-counting before dedup.
- **HIGH→localized (code):** `REGULON_GENES` contains a dead `'Rv1941c'` (GFF has `Rv1941`) and omits Rv1965–Rv1977 → systematic under-flagging of `near_regulon` (affects tail/tiers, not the top 4).
- **MEDIUM (code):** `dedup()` 50 bp window can merge genuinely distinct adjacent operators (e.g. closely-spaced looping operators).
- **MEDIUM (code):** `random.seed(42)` re-set inside the sampler → strong/weak null distributions are the _same_ sampled loci (reproducible, but not independent draws).
- **LOW (code):** dead `bisect.bisect_right(..., key=...)` call at line 194 (overwritten by the hand loop); `score_to_pvalue` binary search itself verified correct.

## Six-dimension summary

1. **Statistical validity — FAIL.** No FDR; floor-artifact p-values with no resolving power; uncalibrated ranking heuristic presented as significance (B1,B2,B7).
2. **Circularity / leakage — FAIL.** Single-template PWM, proximity-prior tier gating, non-independent "3 methods," core-identity double-count (B3,B4).
3. **Code correctness — real bugs.** Architecture not enforced (4× inflation), minus-strand orientation, regulon-set errors, dedup merging (B5 + single-model).
4. **Biological novelty — overstated.** #2 known (Santangelo), #3 same-IGR multiplicity, #4 likely FP; defensible _novel_ count ≈ 0 (B6).
5. **Reproducibility — OK-ish.** Seeded and deterministic, but the shared-null and floor make the reported precision meaningless.
6. **Honesty of framing — FAIL as written.** "Discovered 3 novel operators" is not supported; #1 is a control, #2 is confirmatory-of-known.

## What it would take to make a defensible claim

1. **Genome-wide permutation null of the full bipartite statistic** — dinucleotide-preserving genome shuffles, re-run the identical pipeline on each, record the distribution of `adjusted_score` / candidate counts; report empirical FDR per candidate. This is the single most important fix (decides whether #4 survives at all).
2. **Real tail p-values** — analytic PWM score distribution (e.g. TFMPvalue) or EVD fit, not the 2e-6 empirical floor.
3. **Enforce the real architecture** — weak-upstream / strong-downstream; strand-correct minus-strand coordinates.
4. **Remove `near_regulon` and core-identity from significance/confidence** — use only for prioritization, disclosed as a prior; obtain a genuinely independent method (ChIP-seq: Minch 2015 / Turkarslan 2015; or phylogenetic footprinting blind to the template).
5. **Reframe honestly** — #1 control, #2 confirmatory of Santangelo 2009, #3 same-region multiplicity, #4 unvalidated hypothesis.

## Reconciliation with the other pipeline

This matches the independent `mce3r_project` single-motif scan: a rigorous genome-wide FDR finds only
the known operator(s) and no defensible novel sites. **Both pipelines, audited separately, reach the
same conclusion** — the operator is specific and already-characterized; no new high-confidence binding
sites are supported by sequence alone.
