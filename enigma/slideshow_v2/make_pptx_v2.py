"""
Build an EDITABLE 13-slide PowerPoint matching the corrected deck
(slideshow_v2/Finding Mce3R Binding Sites.html / SLIDES_CONTENT_UPDATED.md).

Native text + tables (fully editable), white background, blue/teal palette
matching the HTML deck. Figures pulled from slideshow_v2/assets/.
"""

import os
from pptx import Presentation
from pptx.util import Inches, Pt, Emu
from pptx.dml.color import RGBColor
from pptx.enum.text import PP_ALIGN, MSO_ANCHOR

HERE = os.path.dirname(os.path.abspath(__file__))
ASSETS = os.path.join(HERE, 'assets')
OUT = os.path.join(HERE, 'Finding_Mce3R_Binding_Sites_v2_editable.pptx')
DL = os.path.expanduser('~/Downloads/Finding_Mce3R_Binding_Sites_v2_editable.pptx')

BLUE = RGBColor(0x1D, 0x3A, 0xD0)
TEAL = RGBColor(0x1A, 0xBC, 0x9C)
INK = RGBColor(0x1A, 0x22, 0x38)
GREY = RGBColor(0x5B, 0x67, 0x70)
MINT = RGBColor(0xE9, 0xF9, 0xF3)
WHITE = RGBColor(0xFF, 0xFF, 0xFF)

prs = Presentation()
prs.slide_width = Inches(13.333)
prs.slide_height = Inches(7.5)
BLANK = prs.slide_layouts[6]
SW, SH = prs.slide_width, prs.slide_height


def slide():
    s = prs.slides.add_slide(BLANK)
    bg = s.background.fill
    bg.solid()
    bg.fore_color.rgb = WHITE
    return s


def textbox(s, x, y, w, h):
    tb = s.shapes.add_textbox(Inches(x), Inches(y), Inches(w), Inches(h))
    tf = tb.text_frame
    tf.word_wrap = True
    return tb, tf


def set_run(r, text, size, color, bold=False, italic=False, mono=False):
    r.text = text
    r.font.size = Pt(size)
    r.font.color.rgb = color
    r.font.bold = bold
    r.font.italic = italic
    r.font.name = 'Consolas' if mono else 'Calibri'


def heading(s, text, size=40):
    tb, tf = textbox(s, 0.7, 0.45, 12.0, 1.1)
    p = tf.paragraphs[0]
    set_run(p.add_run(), text, size, BLUE, bold=True)
    return tb


def subheading(s, text):
    tb, tf = textbox(s, 0.72, 1.42, 12.0, 0.5)
    p = tf.paragraphs[0]
    set_run(p.add_run(), text, 16, TEAL, mono=True)


def bullets(s, items, x=0.75, y=1.8, w=7.4, h=5.2, size=18, gap=10):
    """items: list of plain strings OR list of (text, segments) where
    segments = [(txt, {'bold','color','mono','italic'}), ...]."""
    tb, tf = textbox(s, x, y, w, h)
    tf.word_wrap = True
    first = True
    for it in items:
        p = tf.paragraphs[0] if first else tf.add_paragraph()
        first = False
        p.space_after = Pt(gap)
        p.space_before = Pt(0)
        dash = p.add_run()
        set_run(dash, '—  ', size, TEAL, bold=True)
        if isinstance(it, str):
            set_run(p.add_run(), it, size, INK)
        else:
            for seg_txt, opt in it:
                set_run(p.add_run(), seg_txt, size,
                        opt.get('color', INK), bold=opt.get('bold', False),
                        italic=opt.get('italic', False),
                        mono=opt.get('mono', False))
    return tb


def picture(s, path, x, y, w=None, h=None):
    if not os.path.exists(path):
        raise FileNotFoundError(path)
    kw = {}
    if w:
        kw['width'] = Inches(w)
    if h:
        kw['height'] = Inches(h)
    return s.shapes.add_picture(path, Inches(x), Inches(y), **kw)


def pic_fit(s, path, x, y, boxw, boxh):
    """Add picture scaled to fit (contain) inside box, centered."""
    from PIL import Image
    iw, ih = Image.open(path).size
    ar = iw / ih
    bw, bh = boxw, boxh
    if bw / bh > ar:
        h = bh
        w = bh * ar
    else:
        w = bw
        h = bw / ar
    px = x + (bw - w) / 2
    py = y + (bh - h) / 2
    return s.shapes.add_picture(path, Inches(px), Inches(py),
                                width=Inches(w), height=Inches(h))


# ---- Slide 1 — Title ----
s = slide()
tb, tf = textbox(s, 1.0, 2.55, 11.33, 1.4)
tf.paragraphs[0].alignment = PP_ALIGN.CENTER
set_run(tf.paragraphs[0].add_run(), 'Finding Mce3R Binding Sites', 54, BLUE,
        bold=True)
tb, tf = textbox(s, 1.0, 3.95, 11.33, 0.7)
tf.paragraphs[0].alignment = PP_ALIGN.CENTER
set_run(tf.paragraphs[0].add_run(),
        'Identifying where the Mce3R repressor binds across '
        'M. tuberculosis', 22, GREY, italic=True)
tb, tf = textbox(s, 1.0, 4.95, 11.33, 0.5)
tf.paragraphs[0].alignment = PP_ALIGN.CENTER
set_run(tf.paragraphs[0].add_run(), 'Aayan Alwani', 18, INK)

# ---- Slide 2 — Background ----
s = slide()
heading(s, 'Background')
bullets(s, [
    'Mce3R (Rv1963c) is a TetR-family transcriptional repressor in '
    'M. tuberculosis H37Rv',
    'Panagoda 2024 (PDB 9B7Y) solved its structure — the binding site is '
    'unusually asymmetric',
    'Strong half-site: Kd = 2.4 nM',
    'Weak half-site: Kd = 49.0 nM',
    '20.4-fold affinity difference',
    'To study how architecture affects regulation, we need to find every '
    'binding site',
], w=11.8)

# ---- Slide 3 — Pipeline ----
s = slide()
heading(s, 'A five-step computational pipeline')
bullets(s, [
    '1 · download_genomes.py → retrieve 3 reference genomes',
    '2 · extract_upstream.py → extract 200 bp upstream of each gene',
    '3 · run_meme.py → discover binding motif',
    '4 · run_fimo.py → scan genome for matches',
    '5 · conservation_check.py → validate across species',
], w=11.8)

# ---- Slide 4 — Step 1 ----
s = slide()
heading(s, 'Step 1 — Download reference genomes')
subheading(s, 'phase1_pipeline/download_genomes.py')
bullets(s, [
    'NCBI Entrez download',
    'H37Rv + M. bovis + M. marinum',
    'FASTA + GenBank + GFF3 per species',
    'Three species enables cross-species validation',
], y=2.0, w=11.8)

# ---- Slide 5 — Step 2 ----
s = slide()
heading(s, 'Step 2 — Extract regulatory regions')
subheading(s, 'phase1_pipeline/extract_upstream.py')
bullets(s, [
    '~4,000 protein-coding genes parsed from GFF3',
    '200 bp upstream extracted per gene',
    'Strand-aware, circular-genome-aware',
    'Five specific intergenic regions also extracted (Kabir 2021 protocol)',
], y=2.0, w=11.8)

# ---- Slide 6 — Step 3 (MEME) + logo ----
s = slide()
heading(s, 'Step 3 — Discover the binding motif')
subheading(s, 'phase1_pipeline/run_meme.py')
bullets(s, [
    'MEME Suite (Bailey 2009)',
    'Run on 3 yrbE3A ortholog sequences',
    'ZOOPS mode, widths 20–30 bp, both strands',
    'No palindrome constraint — Mce3R binds asymmetrically',
], y=2.0, w=6.0, size=17)
pic_fit(s, os.path.join(ASSETS, 'logo_mce3r_fimo.png'), 7.0, 2.1, 5.8, 3.4)

# ---- Slide 7 — Step 4 (FIMO) table ----
s = slide()
heading(s, 'Step 4 — Scan the genome (FIMO)')
bullets(s, [
    '1,442 candidate binding sites identified',
    'Top 2 hits match the cryo-EM binding sites at yrbE3A (Panagoda 2024)',
], y=1.55, w=11.8, size=16, gap=4)
rows = [
    ('Rank', 'Gene', 'p-value', 'Score'),
    ('1', 'Rv1964 (yrbE3A)', '1.1e-11', '34.9'),
    ('2', 'Rv1964 (yrbE3A)', '3.1e-11', '33.2'),
    ('3', 'Rv1962A (vapB35)', '4.8e-11', '32.6'),
    ('4', 'Rv0304c (PPE5)', '1.5e-07', '15.8'),
    ('5', 'Rv0251c (hsp)', '3.0e-07', '14.8'),
    ('6', 'Rv2378c (mbtG)', '3.4e-07', '14.8'),
    ('7', 'Rv0412c', '3.5e-07', '11.9'),
    ('8', 'Rv1393c', '3.5e-07', '11.9'),
]
gt = s.shapes.add_table(len(rows), 4, Inches(2.6), Inches(2.5),
                        Inches(8.1), Inches(4.3)).table
for ci, wdt in enumerate((1.2, 3.5, 2.0, 1.4)):
    gt.columns[ci].width = Inches(wdt)
for ri, row in enumerate(rows):
    for ci, val in enumerate(row):
        c = gt.cell(ri, ci)
        c.text = val
        pr = c.text_frame.paragraphs[0]
        pr.alignment = PP_ALIGN.CENTER
        rn = pr.runs[0]
        rn.font.size = Pt(14)
        rn.font.bold = (ri == 0)
        rn.font.color.rgb = WHITE if ri == 0 else INK
        c.fill.solid()
        c.fill.fore_color.rgb = BLUE if ri == 0 else (
            MINT if ri % 2 else WHITE)

# ---- Slide 8 — Step 5 conservation ----
s = slide()
heading(s, 'Step 5 — Conservation across species')
subheading(s, 'phase1_pipeline/conservation_check.py')
bullets(s, [
    'Top 50 sites checked in M. bovis and M. marinum',
    '50 / 50 conserved in M. bovis',
    '44 / 50 conserved in both species',
    'Conservation across divergent species supports biological function',
], y=2.0, w=11.8)

# ---- Slide 9 — Extended MEME 99 bp + logo ----
s = slide()
heading(s, 'Extended MEME recovers a longer 99 bp motif', size=34)
bullets(s, [
    'Per Dr. Balázsi: the Panagoda 2024 123 bp operator is not the only '
    'gold standard',
    [('Re-ran MEME on ', {}), ('Regions 1 + 2 combined', {'italic': True}),
     (' with ', {}), ('-maxw 120 -mod anr', {'mono': True, 'color': BLUE}),
     (' (Eram Kabir 2021 protocol)', {})],
    [('Result: a ', {}), ('99 bp motif', {'bold': True, 'color': BLUE}),
     (', E = 3.9 × 10⁻⁸ — ~4 orders of magnitude better than Eram’s best '
      '(E = 2.8 × 10⁻⁴)', {})],
    [("Eram’s ", {}), ('49 bp motif', {'bold': True, 'color': BLUE}),
     (' is contained within the 99 bp motif at ', {}),
     ('100% identity', {'bold': True, 'color': BLUE})],
    'Three sites: two in Region 1 (mce3R–yrbE3A), one in Region 2 '
    '(echA13–Rv1936)',
], y=1.7, w=6.4, size=16)
pic_fit(s, os.path.join(ASSETS, 'logo_99bp_denovo.png'), 7.3, 1.9, 5.5, 4.3)

# ---- Slide 10 — Second region (Region 2) table ----
s = slide()
heading(s, 'A second regulatory region confirmed')
subheading(s, 'Region 2 · echA13 (Rv1935c) – Rv1936 · 224 bp')
bullets(s, [
    'Original 21 bp FIMO scan missed it (1 weak hit, rank 1422)',
    'The extended 99 bp MEME search found a strong hit in Region 2',
    'First computational recovery of this region with a matching motif',
], y=1.95, w=11.8, size=16, gap=4)
rows = [
    ('Region', 'Position', 'Strand', 'p-value', 'Eram 49 bp'),
    ('1', '684', '+', '1.6e-45', '100% identity'),
    ('1', '238', '+', '2.3e-44', '73% identity'),
    ('2', '62', '−', '4.2e-45', '67% identity'),
]
gt = s.shapes.add_table(len(rows), 5, Inches(3.0), Inches(4.0),
                        Inches(7.3), Inches(2.3)).table
for ri, row in enumerate(rows):
    for ci, val in enumerate(row):
        c = gt.cell(ri, ci)
        c.text = val
        pr = c.text_frame.paragraphs[0]
        pr.alignment = PP_ALIGN.CENTER
        rn = pr.runs[0]
        rn.font.size = Pt(15)
        rn.font.bold = (ri == 0 or ri == 1)
        rn.font.color.rgb = WHITE if ri == 0 else INK
        c.fill.solid()
        c.fill.fore_color.rgb = BLUE if ri == 0 else (
            MINT if ri == 1 else WHITE)

# ---- Slide 11 — Three clusters ABC/DEF/HXG + zoomed map ----
s = slide()
heading(s, 'Three 99 bp motif clusters — ABC, DEF, HXG', size=32)
bullets(s, [
    'The 99 bp motif occurs three times; each occurrence is one of three '
    'suspected binding clusters in Dr. Balázsi’s annotated schematic '
    '(TBmotifs.pdf)',
    [('ABC ', {'bold': True, 'color': BLUE, 'mono': True}),
     ('mce3R–yrbE3A, proximal (near mce3R): 99 bp at position 238 (+), '
      'p = 2.3 × 10⁻⁴⁴; sub-sites A, B, C', {})],
    [('DEF ', {'bold': True, 'color': BLUE, 'mono': True}),
     ('mce3R–yrbE3A, distal (near yrbE3A): 99 bp at position 684 (+), '
      'p = 1.6 × 10⁻⁴⁵; contains the 49 bp core at 100% identity; '
      'sites D, E, F', {})],
    [('HXG ', {'bold': True, 'color': BLUE, 'mono': True}),
     ('echA13–Rv1936, minus strand: 99 bp at position 62 (−), '
      'p = 4.2 × 10⁻⁴⁵; sub-sites H, X, G', {})],
    [('Every sub-site core matches a published instance exactly '
      '(0 mismatches); equivalence ', {}),
     ('A = C = D = F = rcH = rcG', {'mono': True, 'color': BLUE}),
     (' and ', {}), ('B = D = X', {'mono': True, 'color': BLUE})],
], y=1.55, w=5.5, size=13, gap=8)
pic_fit(s, os.path.join(ASSETS, 'dna_motif_map_zoomed.png'),
        6.2, 1.4, 6.9, 5.7)

# ---- Slide 12 — Analyses converge: 4 logos + zoomed map ----
s = slide()
heading(s, 'Multiple independent analyses converge on the same region',
        size=28)
logos = ['logo_mce3r_fimo.png', 'logo_eram_motif2.png',
         'logo_santangelo.png', 'logo_99bp_denovo.png']
ly = 1.5
for lg in logos:
    pic_fit(s, os.path.join(ASSETS, lg), 0.6, ly, 4.4, 1.32)
    ly += 1.42
pic_fit(s, os.path.join(ASSETS, 'dna_motif_map_zoomed.png'),
        5.4, 1.45, 7.6, 5.6)

# ---- Slide 13 — Novel candidates & next steps ----
s = slide()
heading(s, 'Novel candidates & what’s next')
bullets(s, [
    'Rank 3 — Rv1962A (vapB35): adjacent toxin–antitoxin locus',
    'Rank 4 — Rv0304c (PPE5): PE/PPE family',
    'Rank 6 — Rv2378c (mbtG): mycobactin biosynthesis',
    'Rank 38 — Rv1963c: the mce3R gene itself (autoregulation)',
], w=11.8, y=1.7, size=18, gap=6)
tb, tf = textbox(s, 0.75, 4.3, 11.8, 2.6)
p = tf.paragraphs[0]
set_run(p.add_run(), 'Next steps', 20, BLUE, bold=True)
for t in ['Apply the 99 bp PWM genome-wide via FIMO for a broader search',
          'Examine the remaining three intergenic regions (3, 4, 5)',
          'Cross-validate novel candidates with published ChIP-seq data']:
    pp = tf.add_paragraph()
    pp.space_before = Pt(4)
    d = pp.add_run()
    set_run(d, '—  ', 16, TEAL, bold=True)
    set_run(pp.add_run(), t, 16, INK)

prs.save(OUT)
import shutil
shutil.copyfile(OUT, DL)
print(f"Slides: {len(prs.slides.__iter__.__self__._sldIdLst)}")
print(f"Saved: {OUT}")
print(f"Saved: {DL}")
