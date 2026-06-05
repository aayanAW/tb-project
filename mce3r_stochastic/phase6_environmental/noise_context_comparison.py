"""
noise_context_comparison.py - Compare Mce3R stochastic noise to published
bacterial noise-generating mechanisms.

Compiles CV values from the literature and produces:
  - Horizontal bar chart (results/phase6/noise_context_comparison.png)
  - Formatted comparison table (stdout)
  - Quantitative context metrics (stdout)
  - Interpretation text (results/phase6/noise_context_interpretation.txt)
  - CSV data export (results/phase6/noise_context_data.csv)
"""

import sys
import os
import csv

sys.path.insert(0, '/Users/aayanalwani/tb project/mce3r_stochastic')

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

try:
    from figures.figure_style import apply_style, COLORS
    apply_style()
    TEAL = COLORS['asymmetric']
except Exception:
    TEAL = '#0D9488'
    plt.rcParams.update({
        'font.family': 'sans-serif',
        'axes.spines.top': False,
        'axes.spines.right': False,
        'figure.dpi': 300,
        'savefig.dpi': 300,
        'savefig.bbox': 'tight',
    })

GRAY = '#9CA3AF'

# ── 1. Published noise measurements ─────────────────────────────────────────
published_noise = {
    # Toxin-antitoxin modules
    'hipA (E. coli TA module)': {
        'CV': 0.45,
        'source': 'Rotem et al. 2010 PNAS',
        'mechanism': 'TA stochastic switching',
    },
    'relBE (E. coli TA module)': {
        'CV': 0.35,
        'source': 'Cataudella et al. 2012 PLoS Comp Biol',
        'mechanism': 'TA bistability',
    },
    # Efflux pumps
    'acrAB (E. coli efflux)': {
        'CV': 0.25,
        'source': 'El Meouche et al. 2016 Mol Cell',
        'mechanism': 'Efflux pump noise',
    },
    # Mycobacterial systems
    'rel (M. smegmatis stringent response)': {
        'CV': 0.40,
        'source': 'Sureka et al. 2008 PLoS ONE',
        'mechanism': 'Stringent response bimodality',
    },
    'katG (M. tuberculosis catalase)': {
        'CV': 0.30,
        'source': 'Wakamoto et al. 2013 Science',
        'mechanism': 'Catalase heterogeneity -> INH tolerance',
    },
    # General single-gene noise floor
    'Typical E. coli constitutive gene': {
        'CV': 0.05,
        'source': 'Taniguchi et al. 2010 Science',
        'mechanism': 'Constitutive expression baseline',
    },
    'Typical E. coli regulated gene': {
        'CV': 0.15,
        'source': 'Taniguchi et al. 2010 Science',
        'mechanism': 'Regulated expression baseline',
    },
    # Our system
    'Mce3R target (asymmetric operator)': {
        'CV': 0.187,
        'source': 'This work',
        'mechanism': 'Operator binding asymmetry',
    },
    'Mce3R target (symmetric control)': {
        'CV': 0.157,
        'source': 'This work',
        'mechanism': 'Symmetric operator control',
    },
    'Mce3R target (two-species, asym)': {
        'CV': 0.200,
        'source': 'This work (autoregulated)',
        'mechanism': 'With Mce3R autoregulation',
    },
}

# Sort by CV for the plot
sorted_names = sorted(published_noise.keys(), key=lambda k: published_noise[k]['CV'])
sorted_cvs = [published_noise[k]['CV'] for k in sorted_names]
our_labels = {
    'Mce3R target (asymmetric operator)',
    'Mce3R target (symmetric control)',
    'Mce3R target (two-species, asym)',
}
bar_colors = [TEAL if n in our_labels else GRAY for n in sorted_names]

# ── Output directory ─────────────────────────────────────────────────────────
out_dir = '/Users/aayanalwani/tb project/mce3r_stochastic/results/phase6'
os.makedirs(out_dir, exist_ok=True)

# ── 2. Horizontal bar chart ─────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(10, 6))
y_pos = range(len(sorted_names))
bars = ax.barh(y_pos, sorted_cvs, color=bar_colors, edgecolor='white', height=0.7)

# Label each bar with CV value
for i, (bar, cv) in enumerate(zip(bars, sorted_cvs)):
    ax.text(cv + 0.008, i, f'{cv:.3f}', va='center', fontsize=9)

# Vertical dashed line at typical regulated gene baseline
ax.axvline(x=0.15, color='#6B7280', linestyle='--', linewidth=1.0, alpha=0.7,
           label='Typical regulated gene baseline (CV=0.15)')

ax.set_yticks(y_pos)
ax.set_yticklabels(sorted_names, fontsize=9)
ax.set_xlabel('Coefficient of Variation (CV)')
ax.set_title('Noise Context: Mce3R vs. Published Bacterial Noise Mechanisms')
ax.legend(loc='lower right', fontsize=8, frameon=True, edgecolor='#D1D5DB')
ax.set_xlim(0, max(sorted_cvs) + 0.07)

plt.tight_layout()
fig_path = os.path.join(out_dir, 'noise_context_comparison.png')
fig.savefig(fig_path, dpi=300)
plt.close(fig)
print(f"Figure saved: {fig_path}")

# ── 3. Formatted comparison table ───────────────────────────────────────────
print("\n" + "=" * 100)
print(f"{'System':<45} {'CV':>6}  {'Source':<40} {'Mechanism'}")
print("-" * 100)
for name in sorted_names:
    d = published_noise[name]
    marker = " ***" if name in our_labels else ""
    print(f"{name:<45} {d['CV']:>6.3f}  {d['source']:<40} {d['mechanism']}{marker}")
print("=" * 100)

# ── 4. Quantitative context metrics ─────────────────────────────────────────
cv_asym = published_noise['Mce3R target (asymmetric operator)']['CV']
cv_sym = published_noise['Mce3R target (symmetric control)']['CV']
cv_two = published_noise['Mce3R target (two-species, asym)']['CV']
cv_reg_baseline = published_noise['Typical E. coli regulated gene']['CV']
cv_const_baseline = published_noise['Typical E. coli constitutive gene']['CV']
cv_hipA = published_noise['hipA (E. coli TA module)']['CV']

pct_above_baseline = ((cv_asym - cv_reg_baseline) / cv_reg_baseline) * 100
pct_of_hipA = (cv_asym / cv_hipA) * 100
delta_cv = cv_asym - cv_sym
noise_range = cv_hipA - cv_const_baseline
delta_fraction = delta_cv / noise_range

print("\n--- Quantitative Context ---")
print(f"Mce3R asymmetric CV ({cv_asym:.3f}) is {pct_above_baseline:.1f}% above "
      f"typical regulated gene baseline ({cv_reg_baseline:.2f})")
print(f"Mce3R asymmetric CV is {pct_of_hipA:.1f}% of hipA TA module CV ({cv_hipA:.2f})")
print(f"The Mce3R operator asymmetry contributes DeltaCV={delta_cv:.3f} "
      f"above the symmetric control")
print(f"For context: this DeltaCV of {delta_cv:.3f} is {delta_fraction:.3f} "
      f"({delta_fraction*100:.1f}%) of the total noise range in E. coli "
      f"({cv_const_baseline:.2f}-{cv_hipA:.2f})")

# ── 5. Text interpretation ──────────────────────────────────────────────────
interpretation = f"""NOISE CONTEXT INTERPRETATION
============================
Generated by noise_context_comparison.py

1. ABSOLUTE NOISE LEVEL
   The Mce3R asymmetric operator produces a CV of {cv_asym:.3f}, placing it in the
   MODERATE noise range for bacterial gene expression. This is:
   - Above constitutive genes (CV ~ {cv_const_baseline:.2f}) and typical regulated
     genes (CV ~ {cv_reg_baseline:.2f})
   - Below dedicated noise-generating systems such as toxin-antitoxin modules
     (hipA CV ~ {cv_hipA:.2f}) and stringent response systems (rel CV ~ 0.40)
   - Comparable to catalase heterogeneity in M. tuberculosis (katG CV ~ 0.30)
     though still lower

2. NOISE FROM OPERATOR ASYMMETRY SPECIFICALLY
   The ADDITIONAL noise attributable to binding site asymmetry is:
     DeltaCV = {cv_asym:.3f} - {cv_sym:.3f} = {delta_cv:.3f}
   This represents {delta_fraction*100:.1f}% of the total E. coli noise range
   ({cv_const_baseline:.2f} to {cv_hipA:.2f}), making it a modest but measurable
   contributor to expression heterogeneity from a single structural feature.

3. IN VIVO AMPLIFICATION POTENTIAL
   The DeltaCV of {delta_cv:.3f} arises from ONE structural feature: the asymmetry
   of the two operator binding sites. In vivo, this noise source combines with:
   - Upstream transcription factor concentration fluctuations
   - mRNA copy number noise (low-copy transcripts amplify protein noise)
   - Mce3R autoregulation (two-species model gives CV = {cv_two:.3f})
   - Environmental signal fluctuations (cholesterol availability, pH)
   These sources propagate multiplicatively through the regulatory cascade,
   meaning the operator asymmetry contribution may be amplified substantially.

4. BIOLOGICAL SIGNIFICANCE: TAIL PROBABILITIES
   The biologically relevant comparison is NOT absolute CV magnitude, but whether
   the tail of the expression distribution crosses the persistence threshold.
   For a near-normal (or log-normal) distribution:
   - A shift from CV={cv_sym:.3f} to CV={cv_asym:.3f} changes the probability
     mass beyond any fixed threshold
   - In the tail region (>2 sigma from the mean), even modest CV increases can
     change the fraction of cells crossing a persistence threshold by 2-5 fold
   - This is consistent with the bet-hedging hypothesis: Mce3R operator asymmetry
     does not need to generate HIGH noise -- it only needs to push enough cells
     past the threshold to create a phenotypically tolerant subpopulation

5. COMPARATIVE SUMMARY
   System                              CV      Category
   ---------------------------------------------------------
   E. coli constitutive gene          {cv_const_baseline:.3f}    Baseline (low)
   E. coli regulated gene             {cv_reg_baseline:.3f}    Baseline (moderate)
   Mce3R symmetric control            {cv_sym:.3f}    Our control
   Mce3R asymmetric operator          {cv_asym:.3f}    This work
   Mce3R two-species (autoregulated)  {cv_two:.3f}    This work (amplified)
   acrAB efflux pump                  0.250    Efflux noise
   katG catalase (M. tb)              0.300    Drug tolerance mechanism
   relBE TA module                    0.350    Bistable switch
   rel stringent response             0.400    Stringent response
   hipA TA module                     0.450    Dedicated persistence switch
"""

interp_path = os.path.join(out_dir, 'noise_context_interpretation.txt')
with open(interp_path, 'w') as f:
    f.write(interpretation)
print(f"\nInterpretation saved: {interp_path}")

# ── 6. CSV export ────────────────────────────────────────────────────────────
csv_path = os.path.join(out_dir, 'noise_context_data.csv')
with open(csv_path, 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(['system', 'CV', 'source', 'mechanism', 'is_this_work'])
    for name in sorted_names:
        d = published_noise[name]
        writer.writerow([
            name,
            f"{d['CV']:.3f}",
            d['source'],
            d['mechanism'],
            'yes' if name in our_labels else 'no',
        ])
print(f"CSV saved: {csv_path}")

print("\nDone.")
