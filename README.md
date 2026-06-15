# Mce3R Operator Discovery

**Goal:** Panagoda et al. 2024 (cryo-EM, PDB 9B7Y) state Mce3R has **4 operator sites in its
regulon but structurally mapped only 1** (the mce3R–yrbE3A site). This project computationally
**locates/maps the 3 additional operator sites** within the regulon.

> Honest status (cross-model audited — see `operator_discovery/bipartite/audit_report_bipartite.md`
> and `operator_discovery/validation/audit_report.md`): the pipelines recover the known operator
> and map candidate sites in the **literature-corroborated regulon regions** (Rv1935c–Rv1936; the
> mce3R autoregulatory region), but the headline statistics (floor-artifact p-values, no
> genome-wide FDR) are **not** defensible as written, and the out-of-regulon candidate (Rv1115)
> is likely a false positive. Frame results as **mapping additional in-regulon sites**, not novel
> discovery. The remaining work is a permutation null + FDR, or experimental validation (EMSA/ChIP).

## Layout

```
tb project/
├── operator_discovery/              ← THE PROJECT
│   ├── bipartite/                   primary method: structure-guided bipartite search
│   │   ├── find_operators_v3.py     paired 25bp half-sites ~53bp apart (Fisher + spacer weight)
│   │   ├── phase1_pipeline/         genome download, MEME/FIMO, conservation, cross-validate
│   │   ├── config/  data/  results/ (genomes + phase1 outputs; relative paths — keep layout)
│   │   ├── reference_docs/          Eram summary, TBmotifs
│   │   ├── paper_fulltext.txt        Panagoda 2024 full text
│   │   ├── METHODOLOGY_DOCUMENT.*    (in results/) the STS methods writeup
│   │   └── audit_report_bipartite.md cross-model audit of the discovery method
│   ├── validation/                  single-motif FIMO + non-circular gates (rigor SUPPORT)
│   │   └── (run: python main.py — see its CLAUDE.md)
│   └── Email_to_Balazsi_draft.txt
│
├── enigma/                          ← ENIGMA noise-modeling project (likely to be scrapped)
│   │                                  Gillespie sims, thermodynamics, persistence — phases 2–9.
│   │                                  Does NOT depend on finding new operators (uses the known
│   │                                  operator's Kd asymmetry). Self-contained archive (has its
│   │                                  own config/). Its main.py couples to phase1 (now in
│   │                                  ../operator_discovery/bipartite) — re-point if revived.
│   └── ...
│
└── poster_tbc2026/                  conference poster assets
```

## Running

```bash
# Bipartite discovery (heavy — scans the genome):
cd operator_discovery/bipartite && python find_operators_v3.py

# Single-motif validation pipeline (gates):
cd operator_discovery/validation && PATH="$HOME/miniforge3/envs/meme_x64/bin:$PATH" python main.py
```

MEME Suite is required by both and is NOT on PATH — it lives in the `meme_x64` conda env
(`~/miniforge3/envs/meme_x64/bin`). See `operator_discovery/validation/CLAUDE.md` for details.
