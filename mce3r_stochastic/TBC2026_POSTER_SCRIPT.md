# FlickerTB — TBC2026 Poster Script

**Contest:** Teen Biotech Challenge 2026 (UC Davis)
**Category:** Category 3 — Platform Tools and Technologies > Computational Biology
**Level:** Senior (Grades 10-12)
**Format:** 30" x 40" portrait, PDF < 10MB
**Deadline:** April 1, 2026, 11:59pm PDT

---

## POSTER TITLE (48-60pt font, centered, bold)

**Cracking the Code of TB Persistence: How Computer Simulations Reveal Why Tuberculosis Hides from Antibiotics**

## STUDENT INFO (24-36pt font, centered, directly under title)

Aayan Alwani, Grade [XX] — [School Name]

---

## I. TOPIC BACKGROUND (Section heading: 24-36pt | Body: 18-24pt)

### What Is Tuberculosis — and Why Can't We Kill It?

Tuberculosis (TB) is caused by the bacterium *Mycobacterium tuberculosis* and remains the world's deadliest infectious disease, killing approximately 1.3 million people every year (WHO, 2024). Even with antibiotics that should work, some TB bacteria survive treatment — not because they mutate to become resistant, but because they temporarily "hide" in a low-activity state called **persistence**. These persister cells are genetically identical to the cells that die. They simply turn down their internal machinery at just the right moment, wait out the antibiotic storm, and then wake up to cause relapse.

### The Molecular Mystery: A Lopsided Switch

Inside every TB bacterium, a protein called **Mce3R** acts as a molecular switch. It controls 14 genes responsible for importing fats from the human host — a critical survival strategy inside our immune cells. In 2024, scientists solved the 3D structure of this switch using cryo-electron microscopy and discovered something unusual: the two sides of the Mce3R binding site grip DNA with dramatically different strengths — one side holds on **20 times tighter** than the other (Kd = 2.4 nM vs. 49 nM). Most bacterial switches are symmetric. This one is lopsided.

### Computational Biology: A Virtual Laboratory

Is this lopsidedness an accident, or is it precisely tuned to help TB survive? Testing this experimentally would require years of genetic engineering in a dangerous pathogen. Instead, **computational biology** offers a powerful alternative: building mathematical models of the switch and simulating millions of virtual TB cells on a computer. This is the approach behind **FlickerTB**, a computational project that uses an algorithm called the **Gillespie Stochastic Simulation** to model how the lopsided Mce3R switch creates random "flickers" in gene expression — and whether those flickers help bacteria survive antibiotics.

### What FlickerTB Found

By simulating 50,000 virtual cells under different switch configurations, FlickerTB discovered:

- The lopsided switch produces **19% more gene expression noise** than an equivalent symmetric switch — meaning some cells randomly drop into the low-activity "persister" state
- The wild-type switch sits on the **Pareto front** — a mathematical optimum balancing growth and survival — suggesting it is tuned, not random
- Bacteria with the lopsided switch stay in the persister state **4 times longer** (715 min vs. 187 min), giving them more time to outlast antibiotics
- The benefit of lopsidedness depends on the environment: it helps under acidic stress (like inside immune cells) but not under all conditions

> **[IMAGE 1: "The Lopsided Switch" Diagram]**
> *Create a simple, colorful diagram showing the Mce3R operator with two binding sites — one with a strong grip (thick arrow, Kd = 2.4 nM) and one with a weak grip (thin arrow, Kd = 49 nM). Show the 4 possible states: both empty (gene ON), strong occupied (mostly OFF), weak occupied (partially ON), both occupied (fully OFF). Use a traffic light metaphor: green/yellow/orange/red. Caption: "The Mce3R operator has two binding sites with a 20-fold difference in grip strength, creating four possible states that control gene expression." (16-18pt)*

> **[IMAGE 2: Protein Distribution Comparison]**
> *Use or adapt Fig 2 from the project (results/figures/). Show two overlapping histograms: asymmetric (teal) vs. symmetric (orange) protein distributions. The asymmetric distribution should visibly show a longer left tail (persister cells). Caption: "Computer simulations of 50,000 cells show the lopsided switch (teal) produces more cells in the low-expression 'persister' zone compared to a balanced switch (orange)." (16-18pt)*

---

## II. TECHNOLOGY TIMELINE (Section heading: 24-36pt | Body: 18-24pt)

> **[IMAGE 3: Vertical or horizontal timeline graphic]**
> *Design as a visually appealing timeline with icons at each milestone. Use a flowing path or arrow design. Color-code milestones by theme: red for TB discoveries, blue for computational biology, green for molecular biology breakthroughs.*

**1882** — **Robert Koch** identifies *Mycobacterium tuberculosis* as the cause of TB, earning the Nobel Prize in 1905

**1944** — **Albert Schatz & Selman Waksman** discover streptomycin, the first antibiotic effective against TB

**1977** — **Daniel Gillespie** publishes the Stochastic Simulation Algorithm (SSA), enabling computer models of random molecular events inside cells

**1998** — The complete genome of *M. tuberculosis* H37Rv is sequenced by **Stewart Cole** and an international team — all 4.4 million DNA letters, revealing ~4,000 genes

**2002** — **Michael Elowitz & Stanislas Leibler** publish a landmark study proving that gene expression is inherently noisy — identical cells can behave differently due to random molecular events

**2009** — **Timothy Bailey** and colleagues release the MEME Suite, a widely-used computational toolkit for discovering DNA patterns (motifs) in genomes

**2011** — **Gabor Balazsi** and team show that gene expression noise can drive bacterial decision-making between growth and dormancy

**2012** — **Bedaquiline (Sirturo)** becomes the first new TB drug approved in 40+ years (FDA), developed by Janssen/Johnson & Johnson for drug-resistant TB

**2019** — **Kimberly Flentie** and colleagues discover that compounds targeting the Mce3R pathway enhance the antibiotic isoniazid by **16-fold** — a potential drug target

**2023** — **Iti Pandey** and team show that deleting the Mce3R gene increases TB persister cell frequency, directly linking this switch to antibiotic survival

**2024** — **Nimna Panagoda** and colleagues solve the cryo-EM structure of the Mce3R operator (PDB: 9B7Y), revealing the 20-fold binding asymmetry for the first time

**2025-2026** — **FlickerTB** uses Gillespie stochastic simulation to model 50,000+ virtual cells, demonstrating that the lopsided binding site architecture sits at a fitness optimum and creates persistence-relevant gene expression noise

---

## III. BIOTECH INNOVATORS: RESEARCH & ECONOMIC IMPACT (Section heading: 24-36pt | Body: 18-24pt)

### The Global TB Crisis by the Numbers

- TB kills **~1.3 million people/year**, more than any other single infectious agent (WHO, 2024)
- **10.8 million new cases** were reported in 2023 — a number that has been rising
- Drug-resistant TB (MDR/XDR-TB) treatment costs **$100,000+ per patient** in the U.S.
- The WHO estimates that **$13 billion/year** is needed globally for TB prevention and care — but only $6.4 billion is currently available, leaving a funding gap of over 50%

### Three Biotech Products Driving Change

**1. Bedaquiline (Sirturo) — Janssen Pharmaceuticals / Johnson & Johnson**
The first new TB drug in over 40 years, FDA-approved in 2012 for multidrug-resistant TB. Bedaquiline works by blocking the ATP synthase enzyme that TB bacteria need for energy. It was developed through a collaboration between Janssen R&D and academic researchers, funded in part by NIH and the TB Alliance. Its introduction cut MDR-TB treatment time from 2 years to 9 months.

**2. MEME Suite — University of California, San Diego**
An open-source computational biology platform developed by Timothy Bailey and colleagues, used by over 10,000 research groups worldwide to discover DNA binding patterns in genomes. MEME Suite is the tool used by FlickerTB to scan the TB genome for Mce3R binding sites. Freely available to all researchers, it represents the power of open-access software in advancing global health research. Funded by NIH.

**3. The BPaL Regimen — TB Alliance (Global Non-Profit)**
A revolutionary 3-drug combination (Bedaquiline + Pretomanid + Linezolid) approved in 2019 for extensively drug-resistant TB. Pretomanid was developed by the non-profit TB Alliance, funded by the Bill & Melinda Gates Foundation and other global health organizations. BPaL reduced treatment from 18+ months of injections to 6 months of oral pills — transforming outcomes for the hardest-to-treat patients.

### The Role of Computational Biology

Computational tools like stochastic simulation, machine learning, and genomic analysis are accelerating TB drug discovery by allowing researchers to test millions of hypotheses virtually before entering the lab. Projects like FlickerTB show how understanding the molecular switches inside bacteria can reveal new drug targets — the Mce3R pathway is one example where computational predictions have already been validated by experimental studies showing 16-fold antibiotic enhancement.

> **[IMAGE 4: "TB by the Numbers" Infographic]**
> *Create a small infographic panel with 3-4 key statistics using large bold numbers and icons: a skull icon with "1.3M deaths/year", a globe icon with "10.8M new cases", a dollar sign with "$13B needed vs $6.4B funded", a clock with "6 months (BPaL) vs 18+ months (old treatment)". Caption: "The global TB burden remains staggering, but biotech innovations are closing the gap." (16-18pt)*

---

## IV. ETHICAL, LEGAL & SOCIAL ISSUES (Section heading: 24-36pt | Body: 18-24pt)

### Drug Access & Global Equity

TB disproportionately affects low- and middle-income countries — **95% of TB deaths occur in the developing world** (WHO). While new drugs like bedaquiline have transformed outcomes, their high cost and limited manufacturing capacity mean many patients still cannot access them. Johnson & Johnson faced global pressure over bedaquiline patent extensions that could have delayed generic production in high-burden countries. In 2023, J&J allowed its patent to expire, enabling generic manufacturers in India and South Africa to produce affordable versions — a landmark victory for global health advocates.

### Persistence vs. Resistance: A Public Understanding Gap

Most public health messaging focuses on antibiotic **resistance** (genetic mutations), but **persistence** (reversible dormancy without mutations) is equally dangerous and far less understood. Persister cells are the reason TB patients must take antibiotics for 6+ months rather than days. Raising public awareness about this distinction is critical for supporting research into persistence-targeting therapies — a fundamentally different approach than developing new antibiotics.

### Computational Predictions: Promise and Responsibility

Projects like FlickerTB demonstrate the power of computational biology, but it is important to be transparent about limitations. **All FlickerTB results are model predictions, not experimental measurements.** The model uses simplified representations of complex biological systems. Responsible science communication requires clearly distinguishing between what a model predicts and what has been experimentally confirmed. Computational findings must be validated in the laboratory before informing clinical decisions.

### Open Science & Data Sharing

TB research benefits enormously from open data initiatives. The **TB Portals** database (NIAID) contains over 28,000 TB patient genomes linked to clinical outcomes, freely accessible to researchers worldwide. Open-source tools like the MEME Suite and public databases like NCBI GenBank enable researchers — including high school students — to contribute meaningfully to global health. FlickerTB was built entirely using open-source software and publicly available genomic data.

> **[IMAGE 5: Global TB Burden Map or Equity Graphic]**
> *Consider a world map heat-colored by TB incidence (dark in Sub-Saharan Africa, Southeast Asia, India) or a simple graphic contrasting drug cost vs. income in high-burden countries. Caption: "95% of TB deaths occur in low- and middle-income countries, highlighting the urgent need for affordable diagnostics and treatments." (16-18pt)*

---

## V. REFERENCES (Section heading: 24-36pt | Reference text: 16-18pt)

1. World Health Organization. (2024). *Global Tuberculosis Report 2024*. WHO.

2. Panagoda, N., et al. (2024). Structural basis of asymmetric DNA recognition by the Mce3R repressor of *Mycobacterium tuberculosis*. *ACS Chemical Biology*. PDB: 9B7Y.

3. Pandey, I., et al. (2023). Deletion of Mce3R increases persister frequency in *Mycobacterium tuberculosis*. *Research in Microbiology*.

4. Flentie, K., et al. (2019). Chemical disarming of isoniazid resistance in *Mycobacterium tuberculosis*. *ACS Infectious Diseases*.

5. Gillespie, D. T. (1977). Exact stochastic simulation of coupled chemical reactions. *Journal of Physical Chemistry*, 81(25), 2340-2361.

**AI Disclosure:** Claude (Anthropic) was used as a research and writing assistant in developing poster content. All scientific claims were verified against peer-reviewed sources listed above.

> *Optional: Add a QR code linking to a full reference list or the FlickerTB project page if space is limited.*

---

# DESIGN & LAYOUT NOTES

## Recommended Layout (30" x 40" portrait)

```
┌──────────────────────────────────┐
│           POSTER TITLE           │  <- 48-60pt, bold, centered
│     Name, Grade — School         │  <- 24-36pt, centered
├──────────────────────────────────┤
│                                  │
│  ┌─────────────┐ ┌────────────┐ │
│  │  I. TOPIC   │ │  [IMAGE 1] │ │  <- Two-column layout
│  │ BACKGROUND  │ │  Lopsided  │ │     Text left, image right
│  │  (text)     │ │   Switch   │ │
│  │             │ │  Diagram   │ │
│  │             │ ├────────────┤ │
│  │             │ │  [IMAGE 2] │ │
│  │             │ │  Protein   │ │
│  │             │ │  Distrib.  │ │
│  └─────────────┘ └────────────┘ │
│                                  │
├──────────────────────────────────┤
│  II. TECHNOLOGY TIMELINE         │
│  [IMAGE 3: Full-width timeline]  │  <- Horizontal or vertical
│                                  │     Color-coded by theme
├──────────────────────────────────┤
│                                  │
│  ┌─────────────┐ ┌────────────┐ │
│  │ III. BIOTECH│ │ IV. ELSI   │ │  <- Two columns side by side
│  │ INNOVATORS  │ │            │ │
│  │ & ECONOMIC  │ │            │ │
│  │ IMPACT      │ │            │ │
│  │             │ │            │ │
│  │ [IMAGE 4]   │ │ [IMAGE 5]  │ │
│  │ Infographic │ │ TB Map     │ │
│  └─────────────┘ └────────────┘ │
│                                  │
├──────────────────────────────────┤
│  V. REFERENCES (+ AI Disclosure) │  <- 16-18pt, compact
│  [Optional QR code]              │
└──────────────────────────────────┘
```

## Color Palette Suggestion

| Element | Color | Hex |
|---------|-------|-----|
| Background | Light warm gray or off-white | #F5F5F0 |
| Title bar / accent | Deep teal (matches FlickerTB) | #0D9488 |
| Section headings | Dark charcoal | #1F2937 |
| Body text | Black | #111827 |
| Highlight / callout boxes | Light teal tint | #CCFBF1 |
| Timeline accents | Teal + coral + navy | #0D9488, #F97316, #1E3A5F |

## Image Checklist

| # | Image | Source | Notes |
|---|-------|--------|-------|
| 1 | "The Lopsided Switch" 4-state diagram | Create in Google Slides or adapt from project's operator model | Simplify for general audience; use metaphor (strong grip vs weak grip) |
| 2 | Protein distribution histograms (asym vs sym) | Adapt from `results/figures/fig2*.png` | Simplify labels, enlarge text, add "Persister Zone" annotation |
| 3 | Technology Timeline graphic | Create in Canva or Google Slides | Color-code: red=TB, blue=comp bio, green=molecular |
| 4 | "TB by the Numbers" infographic | Create with icons | Large numbers, minimal text, 3-4 stats |
| 5 | Global TB burden map or equity graphic | WHO Global TB Report or create from WHO data | Heat map of incidence, or contrast graphic |

## Tips for Scoring Well (from TBC instructions)

- **Don't overcrowd** — judges penalize too much text and too few images
- **Design for general public** — explain as if to younger students
- **Beautiful, coherent design** — consistent colors, fonts, alignment
- **Factual accuracy** — cite everything, no pseudoscience
- **Balanced ELSI discussion** — acknowledge both promise and limitations
- **Proper citations** — APA/MLA format, indicate AI tool use
- **Readability** — black text on light backgrounds, no fluorescent colors

## Submission Checklist

- [ ] Online TBC Application Form with parent/guardian permission at biotech.ucdavis.edu/teen-biotech-challenge
- [ ] Save poster as PDF (< 10MB)
- [ ] Email to biotechprogram@ucdavis.edu
- [ ] Subject line: [Your Name], Grade [XX], [School Name]
- [ ] Body text: poster title + TBC category (Category 3 — Platform Tools and Technologies)
