#!/usr/bin/env python3
"""Generate ENIGMA_Project_Plan.docx — detailed week-by-week execution plan."""

import os
from docx import Document
from docx.shared import Pt, Inches, Cm, RGBColor, Emu
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.enum.table import WD_TABLE_ALIGNMENT
from docx.enum.section import WD_ORIENT
from docx.oxml.ns import qn, nsdecls
from docx.oxml import parse_xml

# ── helpers ──────────────────────────────────────────────────────────────────

def set_cell_shading(cell, color_hex):
    """Apply background shading to a table cell."""
    shading = parse_xml(f'<w:shd {nsdecls("w")} w:fill="{color_hex}" w:val="clear"/>')
    cell._tc.get_or_add_tcPr().append(shading)


def set_cell_borders(cell, color="AAAAAA", sz="4"):
    """Set thin borders on a cell."""
    tc = cell._tc
    tcPr = tc.get_or_add_tcPr()
    borders = parse_xml(
        f'<w:tcBorders {nsdecls("w")}>'
        f'  <w:top w:val="single" w:sz="{sz}" w:space="0" w:color="{color}"/>'
        f'  <w:left w:val="single" w:sz="{sz}" w:space="0" w:color="{color}"/>'
        f'  <w:bottom w:val="single" w:sz="{sz}" w:space="0" w:color="{color}"/>'
        f'  <w:right w:val="single" w:sz="{sz}" w:space="0" w:color="{color}"/>'
        f'</w:tcBorders>'
    )
    tcPr.append(borders)


def styled_table(doc, headers, rows, col_widths=None, header_color="2E5090"):
    """Create a styled table with header row shading and borders."""
    n_cols = len(headers)
    table = doc.add_table(rows=1 + len(rows), cols=n_cols)
    table.alignment = WD_TABLE_ALIGNMENT.CENTER
    table.autofit = False

    # Header row
    for i, h in enumerate(headers):
        cell = table.rows[0].cells[i]
        cell.text = ""
        p = cell.paragraphs[0]
        run = p.add_run(h)
        run.bold = True
        run.font.size = Pt(10)
        run.font.color.rgb = RGBColor(0xFF, 0xFF, 0xFF)
        run.font.name = "Calibri"
        p.alignment = WD_ALIGN_PARAGRAPH.LEFT
        set_cell_shading(cell, header_color)
        set_cell_borders(cell, color="2E5090")

    # Data rows
    for r_idx, row_data in enumerate(rows):
        for c_idx, val in enumerate(row_data):
            cell = table.rows[r_idx + 1].cells[c_idx]
            cell.text = ""
            p = cell.paragraphs[0]
            run = p.add_run(str(val))
            run.font.size = Pt(10)
            run.font.name = "Calibri"
            p.alignment = WD_ALIGN_PARAGRAPH.LEFT
            bg = "F2F2F2" if r_idx % 2 == 0 else "FFFFFF"
            set_cell_shading(cell, bg)
            set_cell_borders(cell, color="CCCCCC")

    # Column widths
    if col_widths:
        for row in table.rows:
            for i, w in enumerate(col_widths):
                row.cells[i].width = Inches(w)

    # Spacing after table
    doc.add_paragraph("")
    return table


def add_title(doc, text, size=16, color="1A1A2E"):
    p = doc.add_paragraph()
    p.alignment = WD_ALIGN_PARAGRAPH.LEFT
    run = p.add_run(text)
    run.bold = True
    run.font.size = Pt(size)
    run.font.name = "Calibri"
    run.font.color.rgb = RGBColor(
        int(color[0:2], 16), int(color[2:4], 16), int(color[4:6], 16)
    )
    return p


def add_section_header(doc, text, size=14, color="2E5090"):
    p = doc.add_paragraph()
    p.alignment = WD_ALIGN_PARAGRAPH.LEFT
    run = p.add_run(text)
    run.bold = True
    run.font.size = Pt(size)
    run.font.name = "Calibri"
    run.font.color.rgb = RGBColor(
        int(color[0:2], 16), int(color[2:4], 16), int(color[4:6], 16)
    )
    p_fmt = p.paragraph_format
    p_fmt.space_before = Pt(18)
    p_fmt.space_after = Pt(6)
    return p


def add_subsection_header(doc, text, size=12, color="2E5090"):
    p = doc.add_paragraph()
    run = p.add_run(text)
    run.bold = True
    run.font.size = Pt(size)
    run.font.name = "Calibri"
    run.font.color.rgb = RGBColor(
        int(color[0:2], 16), int(color[2:4], 16), int(color[4:6], 16)
    )
    p_fmt = p.paragraph_format
    p_fmt.space_before = Pt(12)
    p_fmt.space_after = Pt(4)
    return p


def add_body(doc, text):
    p = doc.add_paragraph()
    run = p.add_run(text)
    run.font.size = Pt(11)
    run.font.name = "Calibri"
    return p


def add_code_block(doc, lines):
    """Add a monospaced block of text (for dependency diagram etc.)."""
    for line in lines:
        p = doc.add_paragraph()
        run = p.add_run(line)
        run.font.size = Pt(9)
        run.font.name = "Consolas"
        run.font.color.rgb = RGBColor(0x33, 0x33, 0x33)
        p_fmt = p.paragraph_format
        p_fmt.space_before = Pt(0)
        p_fmt.space_after = Pt(0)
        p_fmt.line_spacing = Pt(12)


def add_bullet(doc, text, level=0):
    p = doc.add_paragraph(style="List Bullet")
    p.clear()
    run = p.add_run(text)
    run.font.size = Pt(11)
    run.font.name = "Calibri"
    if level > 0:
        p.paragraph_format.left_indent = Inches(0.5 * (level + 1))
    return p


# ── document creation ────────────────────────────────────────────────────────

doc = Document()

# Page setup: 1-inch margins, US Letter
for section in doc.sections:
    section.top_margin = Inches(1)
    section.bottom_margin = Inches(1)
    section.left_margin = Inches(1)
    section.right_margin = Inches(1)
    section.page_width = Inches(8.5)
    section.page_height = Inches(11)

# Default font
style = doc.styles["Normal"]
font = style.font
font.name = "Calibri"
font.size = Pt(11)
style.paragraph_format.line_spacing = Pt(16.5)  # ~1.5 line spacing

# ═══════════════════════════════════════════════════════════════════════════
# TITLE PAGE
# ═══════════════════════════════════════════════════════════════════════════

# Spacer
for _ in range(6):
    doc.add_paragraph("")

p = doc.add_paragraph()
p.alignment = WD_ALIGN_PARAGRAPH.CENTER
run = p.add_run("ENIGMA Project Plan")
run.bold = True
run.font.size = Pt(28)
run.font.name = "Calibri"
run.font.color.rgb = RGBColor(0x1A, 0x1A, 0x2E)

p = doc.add_paragraph()
p.alignment = WD_ALIGN_PARAGRAPH.CENTER
run = p.add_run(
    "Temporal Noise Dynamics and Therapeutic Noise Quenching\n"
    "in the Asymmetric Mce3R Operator of Mycobacterium tuberculosis"
)
run.font.size = Pt(14)
run.font.name = "Calibri"
run.font.color.rgb = RGBColor(0x2E, 0x50, 0x90)

doc.add_paragraph("")

p = doc.add_paragraph()
p.alignment = WD_ALIGN_PARAGRAPH.CENTER
run = p.add_run("Detailed Execution Plan \u2014 Week-by-Week")
run.font.size = Pt(13)
run.font.name = "Calibri"
run.font.color.rgb = RGBColor(0x55, 0x55, 0x55)
run.italic = True

# Page break after title
doc.add_page_break()

# ═══════════════════════════════════════════════════════════════════════════
# 1. PROJECT OVERVIEW
# ═══════════════════════════════════════════════════════════════════════════

add_title(doc, "1. PROJECT OVERVIEW")

add_body(doc,
    "This plan covers the complete execution of the ENIGMA project across 10 phases and "
    "approximately 16 weeks. The project builds on existing completed work (Phases 1\u20137) "
    "and adds two novel analysis modules: temporal noise dynamics (Phase 8) and therapeutic "
    "noise quenching (Phase 9), plus updated figures (Phase 10). The plan is organized into "
    "four tracks: (A) existing completed work that serves as foundation, (B) new code "
    "additions for Phases 8\u201310, (C) analysis execution and figure generation, and "
    "(D) manuscript preparation and competition submission materials. Each week has explicit "
    "deliverables, time estimates, and dependency links to ensure the project remains on "
    "schedule for the project submission deadline."
)

add_subsection_header(doc, "Current Project Status")

status_headers = ["Component", "Status", "Files", "Notes"]
status_rows = [
    ["Phase 1: Motif Discovery", "COMPLETE", "phase1_pipeline/*.py",
     "1,442 sites found, top hit p=1.1e-11"],
    ["Phase 2: Gillespie Simulation", "COMPLETE", "phase2_simulation/*.py",
     "50K cells x 4 conditions, CV(A)=0.187"],
    ["Phase 3: Statistical Analysis", "COMPLETE", "phase3_analysis/*.py",
     "Bootstrap, KS, Cohen\u2019s d, sensitivity"],
    ["Phase 4: Manuscript Figures 1\u20138", "COMPLETE", "phase4_figures/*.py",
     "8 publication-quality PNGs"],
    ["Phase 5: Thermodynamic Calibration", "COMPLETE", "phase5_thermodynamic/*.py",
     "omega=1.08, MCMC converged"],
    ["Phase 6: Environmental Extensions", "COMPLETE", "phase6_environmental/*.py",
     "24 conditions, MI, persistence"],
    ["Phase 7: Extended Figures 9\u201314", "COMPLETE", "phase7_extended_figures/*.py",
     "6 extended PNGs"],
    ["Phase 8: Temporal Dynamics", "NOT STARTED", "\u2014",
     "NEW: autocorrelation, dwell time"],
    ["Phase 9: Noise Quenching", "NOT STARTED", "\u2014",
     "NEW: dose-response, IC50"],
    ["Phase 10: Updated Figures", "NOT STARTED", "\u2014",
     "NEW: figures 15\u201319"],
    ["Manuscript v3", "NOT STARTED", "\u2014",
     "Updated with new framing + results"],
    ["Competition Materials", "NOT STARTED", "\u2014",
     "Poster, slides, abstract"],
]
styled_table(doc, status_headers, status_rows, col_widths=[2.0, 1.0, 1.8, 1.7])

# ═══════════════════════════════════════════════════════════════════════════
# 2. ARCHITECTURE OF NEW CODE
# ═══════════════════════════════════════════════════════════════════════════

doc.add_page_break()
add_title(doc, "2. ARCHITECTURE OF NEW CODE")

# --- 2.1 Phase 8 ---
add_section_header(doc, "2.1  Phase 8: Temporal Noise Dynamics")
add_body(doc, "Directory: phase8_temporal/")
add_body(doc, "Files to create:")

phase8_files = [
    ["__init__.py", "Package init", "\u2014"],
    ["autocorrelation.py",
     "Compute C(tau) from Gillespie traces via FFT method. "
     "Functions: compute_autocorrelation(trace, max_lag), fit_exponential_decay(C_tau), "
     "extract_tau_c(C_tau). Input: protein time traces from Phase 6 trace recordings. "
     "Output: tau_c per architecture per environment.",
     "Phase 6 traces"],
    ["dwell_time.py",
     "Compute dwell time distributions in persister state. "
     "Functions: compute_dwell_times(trace, threshold), fit_dwell_distribution(dwell_times). "
     "Output: dwell time histograms, mean/median dwell time per architecture.",
     "Phase 6 traces + persistence threshold"],
    ["power_spectrum.py",
     "Compute power spectral density from traces. "
     "Functions: compute_psd(trace, dt), identify_peaks(psd, freqs). "
     "Output: PSD plots, characteristic frequencies per architecture.",
     "Phase 6 traces"],
    ["temporal_main.py",
     "Phase 8 orchestrator. Load existing Phase 6 trace data. "
     "Run extended trace simulations (1,000 cells x 3 architectures x 4 environments, "
     "recording full traces at 10-min intervals for 15,000 min post-burn-in = "
     "1,500 time points per cell). Compute autocorrelation, dwell time, PSD for each "
     "condition. Save results to results/phase8/. Self-tests: tau_c > 0, dwell times > 0, "
     "PSD integrates to variance.",
     "All Phase 8 modules"],
]

styled_table(doc,
    ["File", "Description", "Dependencies"],
    phase8_files,
    col_widths=[1.5, 3.5, 1.5])

add_body(doc, "Estimated lines of code: ~400")
add_body(doc, "Estimated runtime: ~30 min (1,000 cells with trace recording x 12 conditions)")

# --- 2.2 Phase 9 ---
add_section_header(doc, "2.2  Phase 9: Noise Quenching Dose-Response")
add_body(doc, "Directory: phase9_quenching/")
add_body(doc, "Files to create:")

phase9_files = [
    ["__init__.py", "Package init", "\u2014"],
    ["symmetrization_sweep.py",
     "Sweep Kd_weak from 49 to 2.4 nM. 15 Kd_weak values: "
     "[49, 40, 35, 30, 25, 20, 15, 10.84, 8, 6, 5, 4, 3, 2.4, 1.5]. "
     "At each: 5,000 cells, record CV, persister fraction, mean protein. "
     "Functions: run_symmetrization_sweep(kd_weak_values, n_cells, master_seed). "
     "Output: symmetrization_sweep.csv.",
     "Gillespie engine"],
    ["concentration_sweep.py",
     "Sweep Mce3R concentration. 15 log-spaced concentrations from 10 to 10,000 nM. "
     "Both asymmetric and symmetric architectures at each. "
     "Functions: run_concentration_sweep(conc_values, architectures, n_cells). "
     "Output: concentration_sweep.csv.",
     "Gillespie engine"],
    ["ic50_analysis.py",
     "Compute IC50 of noise and persistence. Fit sigmoidal dose-response: "
     "CV = CV_min + (CV_max - CV_min) / (1 + (Kd/IC50)^n). "
     "Functions: compute_ic50(kd_values, cv_values), "
     "compute_ic50_persistence(kd_values, persister_fractions). "
     "Output: ic50_results.json.",
     "Symmetrization sweep results"],
    ["quenching_main.py",
     "Phase 9 orchestrator. Run both sweeps, compute IC50 values, "
     "save all results to results/phase9/. Self-tests: IC50 between 2.4 and 49 nM, "
     "CV monotonically decreasing with symmetrization.",
     "All Phase 9 modules"],
]

styled_table(doc,
    ["File", "Description", "Dependencies"],
    phase9_files,
    col_widths=[1.5, 3.5, 1.5])

add_body(doc, "Estimated lines of code: ~500")
add_body(doc,
    "Estimated runtime: ~45 min (15 x 5,000 cells symmetrization + "
    "15 x 2 x 5,000 cells concentration)")

# --- 2.3 Phase 10 ---
add_section_header(doc, "2.3  Phase 10: New Figures")
add_body(doc,
    "Directory: phase10_new_figures/ (or extend phase7_extended_figures/). "
    "Each figure: 300 dpi, publication style matching existing figures "
    "(teal/coral/purple color scheme)."
)

fig_data = [
    ["fig15_autocorrelation.py",
     "C(tau) curves for 3 architectures, fitted exponentials, tau_c values annotated",
     "Phase 8 autocorrelation results"],
    ["fig16_dwell_time.py",
     "Dwell time distributions (histograms or survival curves) for 3 architectures",
     "Phase 8 dwell time results"],
    ["fig17_power_spectrum.py",
     "PSD plots for 3 architectures, characteristic frequencies marked",
     "Phase 8 PSD results"],
    ["fig18_symmetrization_dose_response.py",
     "CV and persister fraction vs Kd_weak, IC50 marked",
     "Phase 9 symmetrization sweep"],
    ["fig19_concentration_sweep.py",
     "CV vs [Mce3R] for asymmetric and symmetric, crossover point marked",
     "Phase 9 concentration sweep"],
]

styled_table(doc,
    ["Figure Script", "Content", "Data Source"],
    fig_data,
    col_widths=[2.0, 2.8, 1.7])

# ═══════════════════════════════════════════════════════════════════════════
# 3. WEEK-BY-WEEK EXECUTION PLAN
# ═══════════════════════════════════════════════════════════════════════════

doc.add_page_break()
add_title(doc, "3. WEEK-BY-WEEK EXECUTION PLAN")

week_headers = ["Week", "Dates", "Phase", "Tasks", "Deliverables", "Dependencies", "Est. Hours"]
week_rows = [
    ["1", "Mar 30 \u2013 Apr 5",  "Setup",
     "Review all existing code and results. Read all cited papers "
     "(Chowdhury 2021, Quigley & Lewis 2022, Lengyel & Morelli 2017). "
     "Write updated project context and handoff docs reflecting new aims.",
     "Updated project documentation",
     "None", "15"],

    ["2", "Apr 6 \u2013 Apr 12", "Phase 8",
     "Design Phase 8 architecture (temporal dynamics). "
     "Write autocorrelation.py with self-tests. "
     "Write dwell_time.py with self-tests.",
     "Two tested analysis modules",
     "Week 1 review complete", "20"],

    ["3", "Apr 13 \u2013 Apr 19", "Phase 8",
     "Write power_spectrum.py with self-tests. "
     "Write temporal_main.py orchestrator. "
     "Run Phase 8 on existing Phase 6 trace data (quick test: 100 cells).",
     "Phase 8 code complete, preliminary results",
     "Week 2 modules", "20"],

    ["4", "Apr 20 \u2013 Apr 26", "Phase 8",
     "Run Phase 8 full production (1,000 cells x 12 conditions). "
     "Analyze results: is tau_c(asymmetric) > tau_c(symmetric)? "
     "Compare dwell time distributions across architectures.",
     "Phase 8 results, initial analysis",
     "Week 3 code complete", "18"],

    ["5", "Apr 27 \u2013 May 3", "Phase 9",
     "Design Phase 9 architecture (noise quenching). "
     "Write symmetrization_sweep.py with self-tests. "
     "Write concentration_sweep.py with self-tests.",
     "Two tested simulation modules",
     "Phase 8 complete", "20"],

    ["6", "May 4 \u2013 May 10", "Phase 9",
     "Write ic50_analysis.py with self-tests. "
     "Write quenching_main.py orchestrator. "
     "Run Phase 9 full production.",
     "Phase 9 results, IC50 values computed",
     "Week 5 modules", "22"],

    ["7", "May 11 \u2013 May 17", "Phase 10",
     "Create fig15_autocorrelation.py. "
     "Create fig16_dwell_time.py. "
     "Create fig17_power_spectrum.py.",
     "3 new figures (temporal dynamics)",
     "Phase 8 results", "15"],

    ["8", "May 18 \u2013 May 24", "Phase 10",
     "Create fig18_symmetrization_dose_response.py. "
     "Create fig19_concentration_sweep.py. "
     "Review all 19 figures for consistency.",
     "5 new figures total, all figures reviewed",
     "Phase 9 results", "15"],

    ["9", "May 25 \u2013 May 31", "Robustness",
     "Posterior propagation for temporal dynamics (50 MCMC samples). "
     "Sensitivity of IC50 to parameter perturbation. "
     "Two-species model temporal dynamics.",
     "Robustness analysis complete",
     "Phases 8\u201310 results", "20"],

    ["10", "Jun 1 \u2013 Jun 7", "Manuscript",
     "Write manuscript Introduction and Methods. "
     "Incorporate all critique responses. "
     "Cite all prior work properly (Chowdhury, Quigley & Lewis, Lengyel & Morelli).",
     "Manuscript draft sections 1\u20132",
     "All results finalized", "25"],

    ["11", "Jun 8 \u2013 Jun 14", "Manuscript",
     "Write manuscript Results and Discussion. "
     "Include temporal dynamics results. "
     "Include noise quenching results. "
     "Honest limitations section with falsification criteria.",
     "Manuscript draft sections 3\u20134",
     "Week 10 draft", "25"],

    ["12", "Jun 15 \u2013 Jun 21", "Manuscript",
     "Write Abstract, References, figure legends. "
     "Internal review and revision.",
     "Complete manuscript draft v3",
     "Week 11 draft", "20"],

    ["13", "Jun 22 \u2013 Jun 28", "Competition",
     "Prepare submission materials: 1-page project abstract, "
     "research report (20-page limit), common application essay edits.",
     "competition abstract draft",
     "Manuscript v3 complete", "20"],

    ["14", "Jun 29 \u2013 Jul 5", "Competition",
     "Prepare presentation materials: 15-slide summary presentation, "
     "poster layout (48x36 inches), 3-minute elevator pitch script.",
     "Presentation materials draft",
     "Week 13 abstract", "18"],

    ["15", "Jul 6 \u2013 Jul 12", "Finalization",
     "Final revisions to all materials. Code cleanup and documentation. "
     "Verify all results reproducible from clean run.",
     "Finalized manuscript, code, materials",
     "Weeks 13\u201314 materials", "20"],

    ["16", "Jul 13 \u2013 Jul 19", "Buffer",
     "Address any remaining issues. Practice presentation. Submit.",
     "Submission-ready package",
     "Week 15 finalization", "10"],
]

styled_table(doc, week_headers, week_rows,
             col_widths=[0.45, 0.95, 0.7, 2.0, 1.2, 0.95, 0.5])

# Total hours
add_body(doc, "Total estimated hours: ~303 hours across 16 weeks (~19 hours/week average)")

# ═══════════════════════════════════════════════════════════════════════════
# 4. DEPENDENCY GRAPH
# ═══════════════════════════════════════════════════════════════════════════

doc.add_page_break()
add_title(doc, "4. DEPENDENCY GRAPH")

add_body(doc,
    "The following diagram shows the dependency structure across all project phases. "
    "Arrows indicate that the target phase requires outputs from the source phase."
)

dep_lines = [
    "Phase 1 (motif discovery) ─────────────────────┐",
    "                                                │",
    "Phase 2 (Gillespie) ────┐                      │",
    "                        ├── Phase 3 (stats) ───├── Phase 5 (thermo) ──┐",
    "Phase 2 (Gillespie) ────┘        │             │                     │",
    "                                 │             │                     │",
    "                         Phase 4 (figs 1-8)    │                     │",
    "                                               │                     │",
    "                                 Phase 6 (environmental) ────────────┤",
    "                                               │                     │",
    "                                 Phase 7 (figs 9-14)                 │",
    "                                               │                     │",
    "                                 Phase 8 (temporal dynamics) <───────┤",
    "                                               │                     │",
    "                                 Phase 9 (noise quenching) <─────────┘",
    "                                               │",
    "                                 Phase 10 (figs 15-19)",
    "                                               │",
    "                                 Manuscript v3",
    "                                               │",
    "                                 Competition materials",
]
add_code_block(doc, dep_lines)

add_body(doc,
    "Critical path: Phase 2 -> Phase 5 -> Phase 6 -> Phase 8 -> Phase 9 -> "
    "Phase 10 -> Manuscript -> Competition. Any delay on this path directly "
    "impacts the final submission date."
)

# ═══════════════════════════════════════════════════════════════════════════
# 5. RISK REGISTER
# ═══════════════════════════════════════════════════════════════════════════

doc.add_page_break()
add_title(doc, "5. RISK REGISTER")

risk_headers = ["#", "Risk", "Probability", "Impact", "Mitigation"]
risk_rows = [
    ["1",
     "tau_c dominated by protein half-life, obscuring architecture-dependent memory",
     "Medium", "High",
     "Still publishable as negative result; report that amplitude not memory is the "
     "relevant noise metric for persister fate decisions"],
    ["2",
     "IC50 falls outside biologically achievable Kd range",
     "Low", "Medium",
     "Report the value honestly; discuss what concentration range would be needed "
     "and whether existing drugs could achieve it"],
    ["3",
     "Phase 8 runtime exceeds estimate due to trace recording overhead",
     "Low", "Low",
     "Reduce to 500 cells; noise statistics converge fast for temporal metrics "
     "given 1,500 time points per trace"],
    ["4",
     "Persister fractions too low to detect architecture differences in dwell time",
     "Medium", "Medium",
     "Use multiple threshold methods (already implemented in Phase 6); report "
     "fold-enrichment rather than absolute fractions"],
    ["5",
     "MCMC posterior too wide to constrain temporal dynamics predictions",
     "Medium", "Medium",
     "Report range of tau_c across posterior; if range is huge, this is itself "
     "informative about which parameters control noise memory"],
    ["6",
     "Figures not reaching publication quality standards",
     "Low", "Low",
     "Use existing figure_style.py infrastructure from Phase 4; iterate on "
     "reviewer feedback before final submission"],
    ["7",
     "Manuscript rejected for lack of experimental validation",
     "High", "Medium",
     "Frame as computational prediction paper with explicit falsification criteria; "
     "target computational biology journals (PLoS Comp Bio, Biophysical Journal)"],
    ["8",
     "competition deadline pressure compresses final preparation time",
     "Medium", "High",
     "Buffer week (Week 16) built in; prioritize manuscript and research report "
     "over presentation materials if needed"],
]

styled_table(doc, risk_headers, risk_rows,
             col_widths=[0.3, 2.0, 0.8, 0.7, 2.7])

# ═══════════════════════════════════════════════════════════════════════════
# 6. FILE INVENTORY
# ═══════════════════════════════════════════════════════════════════════════

doc.add_page_break()
add_title(doc, "6. FILE INVENTORY")

add_body(doc,
    "Comprehensive inventory of all files in the ENIGMA project, organized by "
    "functional group. Items marked (NEW) are to be created in Phases 8\u201310."
)

# Group 1: Configuration
add_subsection_header(doc, "Group 1: Configuration")
styled_table(doc,
    ["File", "Description"],
    [
        ["config/__init__.py", "Package initialization"],
        ["config/parameters.py", "ALL parameters, single source of truth for the entire project"],
    ],
    col_widths=[2.5, 4.0])

# Group 2: Phase 1
add_subsection_header(doc, "Group 2: Phase 1 \u2014 Motif Discovery")
styled_table(doc,
    ["File", "Description"],
    [
        ["phase1_pipeline/__init__.py", "Package initialization"],
        ["phase1_pipeline/download_genomes.py", "Download mycobacterial genomes from NCBI"],
        ["phase1_pipeline/extract_upstream.py", "Extract upstream regions of mce3 operons"],
        ["phase1_pipeline/run_meme.py", "Run MEME motif discovery"],
        ["phase1_pipeline/run_fimo.py", "Run FIMO motif scanning"],
        ["phase1_pipeline/conservation_check.py", "Check motif conservation across species"],
        ["phase1_pipeline/pipeline_main.py", "Phase 1 orchestrator"],
    ],
    col_widths=[2.5, 4.0])

# Group 3: Phase 2
add_subsection_header(doc, "Group 3: Phase 2 \u2014 Gillespie Simulation")
styled_table(doc,
    ["File", "Description"],
    [
        ["phase2_simulation/__init__.py", "Package initialization"],
        ["phase2_simulation/operator_model.py", "Operator state model (4-state asymmetric binding)"],
        ["phase2_simulation/gillespie_engine.py", "Core Gillespie SSA engine with tau-leaping option"],
        ["phase2_simulation/run_conditions.py", "Run simulations across 4 binding conditions"],
        ["phase2_simulation/asymmetry_sweep.py", "Sweep asymmetry parameter space"],
        ["phase2_simulation/simulation_main.py", "Phase 2 orchestrator"],
    ],
    col_widths=[2.5, 4.0])

# Group 4: Phase 3
add_subsection_header(doc, "Group 4: Phase 3 \u2014 Statistical Analysis")
styled_table(doc,
    ["File", "Description"],
    [
        ["phase3_analysis/__init__.py", "Package initialization"],
        ["phase3_analysis/noise_metrics.py", "CV, Fano factor, bimodality index computation"],
        ["phase3_analysis/bootstrap_ci.py", "Bootstrap confidence interval estimation"],
        ["phase3_analysis/statistical_tests.py", "KS test, Mann-Whitney U, permutation tests"],
        ["phase3_analysis/sensitivity_analysis.py", "Parameter sensitivity (Sobol, Morris)"],
        ["phase3_analysis/negative_controls.py", "Symmetric null model controls"],
        ["phase3_analysis/experimental_comparison.py", "Comparison with published experimental data"],
        ["phase3_analysis/analysis_main.py", "Phase 3 orchestrator"],
    ],
    col_widths=[2.5, 4.0])

# Group 5: Phase 4
add_subsection_header(doc, "Group 5: Phase 4 \u2014 Figures 1\u20138")
styled_table(doc,
    ["File", "Description"],
    [
        ["phase4_figures/__init__.py", "Package initialization"],
        ["phase4_figures/figure_style.py", "Shared style: teal/coral/purple palette, fonts, layout"],
        ["phase4_figures/fig1_binding_sites.py", "Figure 1: Operator binding site architecture"],
        ["phase4_figures/fig2_noise_distributions.py", "Figure 2: Protein noise distributions"],
        ["phase4_figures/fig3_cv_comparison.py", "Figure 3: CV comparison across architectures"],
        ["phase4_figures/fig4_asymmetry_sweep.py", "Figure 4: Asymmetry parameter sweep"],
        ["phase4_figures/fig5_bootstrap.py", "Figure 5: Bootstrap confidence intervals"],
        ["phase4_figures/fig6_sensitivity.py", "Figure 6: Sensitivity analysis tornado plot"],
        ["phase4_figures/fig7_experimental.py", "Figure 7: Experimental data comparison"],
        ["phase4_figures/fig8_methods_architecture.py", "Figure 8: Methods and architecture overview"],
        ["phase4_figures/figures_main.py", "Phase 4 orchestrator"],
    ],
    col_widths=[2.5, 4.0])

# Group 6: Phase 5
add_subsection_header(doc, "Group 6: Phase 5 \u2014 Thermodynamic Calibration")
styled_table(doc,
    ["File", "Description"],
    [
        ["phase5_thermodynamic/__init__.py", "Package initialization"],
        ["phase5_thermodynamic/energy_calibration.py", "Free energy calibration from binding data"],
        ["phase5_thermodynamic/partition_function.py", "Statistical mechanical partition function"],
        ["phase5_thermodynamic/cooperativity_inference.py", "MCMC inference of cooperativity (omega)"],
        ["phase5_thermodynamic/operator_classification.py", "Classify operator architectures by energy"],
        ["phase5_thermodynamic/dna_shape.py", "DNA shape parameter analysis"],
        ["phase5_thermodynamic/thermodynamic_main.py", "Phase 5 orchestrator"],
    ],
    col_widths=[2.5, 4.0])

# Group 7: Phase 6
add_subsection_header(doc, "Group 7: Phase 6 \u2014 Environmental Extensions")
styled_table(doc,
    ["File", "Description"],
    [
        ["phase6_environmental/__init__.py", "Package initialization"],
        ["phase6_environmental/environmental_signals.py", "Model environmental signal modulation"],
        ["phase6_environmental/two_species_model.py", "Two-species (mRNA + protein) model"],
        ["phase6_environmental/two_species_gillespie.py", "Gillespie engine for two-species model"],
        ["phase6_environmental/mutual_information.py", "Mutual information between noise and environment"],
        ["phase6_environmental/persistence_threshold.py", "Persister threshold calibration"],
        ["phase6_environmental/run_environmental_conditions.py", "Run 24 environmental conditions"],
        ["phase6_environmental/environmental_main.py", "Phase 6 orchestrator"],
    ],
    col_widths=[2.5, 4.0])

# Group 8: Phase 7
add_subsection_header(doc, "Group 8: Phase 7 \u2014 Extended Figures 9\u201314")
styled_table(doc,
    ["File", "Description"],
    [
        ["phase7_extended_figures/__init__.py", "Package initialization"],
        ["phase7_extended_figures/fig9_*.py", "Figure 9: Thermodynamic landscape"],
        ["phase7_extended_figures/fig10_*.py", "Figure 10: MCMC posterior distributions"],
        ["phase7_extended_figures/fig11_*.py", "Figure 11: Environmental phase diagram"],
        ["phase7_extended_figures/fig12_*.py", "Figure 12: Mutual information heatmap"],
        ["phase7_extended_figures/fig13_*.py", "Figure 13: Persistence fractions across conditions"],
        ["phase7_extended_figures/fig14_*.py", "Figure 14: Two-species model comparison"],
        ["phase7_extended_figures/extended_figures_main.py", "Phase 7 orchestrator"],
    ],
    col_widths=[2.5, 4.0])

# Group 9: Phase 8 (NEW)
add_subsection_header(doc, "Group 9: Phase 8 \u2014 Temporal Dynamics (NEW)")
styled_table(doc,
    ["File", "Description"],
    [
        ["phase8_temporal/__init__.py", "Package initialization"],
        ["phase8_temporal/autocorrelation.py", "FFT-based autocorrelation C(tau) computation"],
        ["phase8_temporal/dwell_time.py", "Dwell time distribution in persister state"],
        ["phase8_temporal/power_spectrum.py", "Power spectral density computation"],
        ["phase8_temporal/temporal_main.py", "Phase 8 orchestrator"],
    ],
    col_widths=[2.5, 4.0])

# Group 10: Phase 9 (NEW)
add_subsection_header(doc, "Group 10: Phase 9 \u2014 Noise Quenching (NEW)")
styled_table(doc,
    ["File", "Description"],
    [
        ["phase9_quenching/__init__.py", "Package initialization"],
        ["phase9_quenching/symmetrization_sweep.py", "Kd_weak sweep from 49 to 2.4 nM"],
        ["phase9_quenching/concentration_sweep.py", "Mce3R concentration sweep (10\u201310,000 nM)"],
        ["phase9_quenching/ic50_analysis.py", "Sigmoidal IC50 fitting for noise and persistence"],
        ["phase9_quenching/quenching_main.py", "Phase 9 orchestrator"],
    ],
    col_widths=[2.5, 4.0])

# Group 11: Phase 10 (NEW)
add_subsection_header(doc, "Group 11: Phase 10 \u2014 New Figures (NEW)")
styled_table(doc,
    ["File", "Description"],
    [
        ["fig15_autocorrelation.py", "C(tau) curves, fitted exponentials, tau_c annotated"],
        ["fig16_dwell_time.py", "Dwell time distributions / survival curves"],
        ["fig17_power_spectrum.py", "PSD plots with characteristic frequencies"],
        ["fig18_symmetrization_dose_response.py", "CV and persister fraction vs Kd_weak, IC50 marked"],
        ["fig19_concentration_sweep.py", "CV vs [Mce3R] for asymmetric/symmetric, crossover marked"],
    ],
    col_widths=[2.5, 4.0])

# Group 12: Results
add_subsection_header(doc, "Group 12: Results Directories")
styled_table(doc,
    ["Directory", "Contents"],
    [
        ["results/phase1/", "MEME/FIMO outputs, predicted_sites.csv, conservation_status.csv"],
        ["results/phase2/", "condition_A/B/C/D.npz, sweep results, phase2_summary.json"],
        ["results/phase3/", "noise_metrics.csv, statistical_tests.csv, bootstrap_results.csv, sensitivity_data.csv"],
        ["results/phase5/", "energy_parameters.json, mcmc_summary.json, mcmc_posteriors.npz, "
         "repression_curves.csv, operator_classifications.csv"],
        ["results/phase6/", "24 condition NPZ files, persistence_fractions.csv, "
         "mutual_information.csv, phase_diagram.csv"],
        ["results/figures/", "fig1\u2013fig8 PNGs (publication quality, 300 dpi)"],
        ["results/extended_figures/", "fig9\u2013fig14 PNGs"],
        ["results/phase8/ (NEW)", "autocorrelation.csv, dwell_times.csv, psd_data.csv, phase8_summary.json"],
        ["results/phase9/ (NEW)", "symmetrization_sweep.csv, concentration_sweep.csv, "
         "ic50_results.json, phase9_summary.json"],
        ["results/new_figures/ (NEW)", "fig15\u2013fig19 PNGs"],
    ],
    col_widths=[2.0, 4.5])

# Group 13: Documentation
add_subsection_header(doc, "Group 13: Documentation")
styled_table(doc,
    ["File", "Description"],
    [
        ["HANDOFF.md", "Session handoff notes and current state"],
        ["PROJECT_PROPOSAL.md", "Original project proposal"],
        ["PHASE2_ARCHITECTURE.md", "Phase 2 design document"],
        ["FULL_PROJECT_ARCHITECTURE.md", "Complete project architecture reference"],
        ["BUILD_REPORT.md", "Build and validation report"],
        ["ENIGMA_Manuscript_v2.docx", "Manuscript version 2 (current)"],
        ["ENIGMA_Project_Proposal_v3.docx", "Project proposal version 3"],
        ["ENIGMA_Project_Plan.docx", "This document"],
    ],
    col_widths=[2.5, 4.0])

# ═══════════════════════════════════════════════════════════════════════════
# 7. SUCCESS CRITERIA
# ═══════════════════════════════════════════════════════════════════════════

doc.add_page_break()
add_title(doc, "7. SUCCESS CRITERIA")

add_body(doc,
    "The following table lists all quantitative success criteria for the project. "
    "Criteria 1\u20135 and 8\u20139 have been validated in completed phases. "
    "Criteria 6\u20137 and 10 are pending completion of new phases."
)

criteria_headers = ["#", "Criterion", "Metric", "Threshold", "Status"]
criteria_rows = [
    ["1", "Motif discovery validates known operator",
     "FIMO #1 hit in operator region", "p < 1e-8",
     "PASS (p=1.1e-11)"],
    ["2", "CV(asymmetric) > CV(symmetric)",
     "Bootstrap 95% CI", "CI entirely above zero",
     "PASS (0.023\u20130.037)"],
    ["3", "Result robust across MCMC posterior",
     "Fraction of samples with CV_asym > CV_sym", ">95%",
     "PASS (100%)"],
    ["4", "Result robust under autoregulation",
     "Two-species CV difference", ">0",
     "PASS (0.200 vs 0.197)"],
    ["5", "Result robust across environments",
     "CV_asym > CV_sym in all 4", "4/4",
     "PASS (4/4)"],
    ["6", "tau_c shows architecture dependence",
     "tau_c(asym) != tau_c(sym)",
     "Difference outside 95% CI of null",
     "PENDING"],
    ["7", "IC50 falls in biologically meaningful range",
     "IC50_noise value", "Between 2.4 and 49 nM",
     "PENDING"],
    ["8", "Persister fractions biologically realistic",
     "Calibrated fractions", "10^-4 to 10^-2 range",
     "PASS (0.0028 baseline)"],
    ["9", "All code passes self-tests",
     "Sanity checks", "0 failures",
     "PASS for phases 1\u20137"],
    ["10", "Manuscript addresses all critiques",
     "9 critiques from review", "All addressed",
     "PENDING"],
]

styled_table(doc, criteria_headers, criteria_rows,
             col_widths=[0.3, 1.8, 1.6, 1.4, 1.4])

# ═══════════════════════════════════════════════════════════════════════════
# 8. COMPETITION TIMELINE (competition)
# ═══════════════════════════════════════════════════════════════════════════

doc.add_page_break()
add_title(doc, "8. COMPETITION TIMELINE (competition)")

add_body(doc,
    "The project has a November 2026 submission deadline. "
    "The following timeline maps project milestones to competition requirements, "
    "ensuring all materials are prepared with adequate buffer time."
)

comp_headers = ["Date", "Milestone", "What\u2019s Needed"]
comp_rows = [
    ["April 2026", "Project plan finalized",
     "This document completed and reviewed; all phases scoped with clear deliverables"],
    ["May 2026", "Phases 8\u20139 code complete",
     "Temporal dynamics and noise quenching code written, tested, and validated"],
    ["June 2026", "All simulations run",
     "All Phase 8\u20139 production runs complete; results analyzed and interpreted"],
    ["July 2026", "All figures complete (19 total)",
     "Figures 1\u201319 finalized at publication quality; consistent style across all"],
    ["August 2026", "Manuscript v3 complete",
     "Full paper draft with Introduction, Methods, Results, Discussion; all critiques addressed"],
    ["September 2026", "Research report draft",
     "20-page competition format research report; adapted from manuscript for competition audience"],
    ["October 2026", "Final revisions",
     "All materials polished; code documented; results verified reproducible from clean run"],
    ["November 2026", "competition submission",
     "Research report, 1-page abstract, transcripts, teacher recommendation, "
     "Common App essay, all uploaded to competition portal"],
]

styled_table(doc, comp_headers, comp_rows,
             col_widths=[1.2, 1.8, 3.5])

add_body(doc,
    "Note: The timeline above assumes a November 2026 deadline. If the deadline is "
    "earlier, compress Weeks 13\u201316 first (competition materials and buffer). "
    "The critical deliverable is the manuscript (Week 12), from which all competition "
    "materials derive."
)

# ── Save ─────────────────────────────────────────────────────────────────

output_path = "/Users/aayanalwani/tb project/mce3r_stochastic/ENIGMA_Project_Plan.docx"
doc.save(output_path)
print(f"Saved: {output_path}")
print(f"Size: {os.path.getsize(output_path):,} bytes")
