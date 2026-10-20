"""Apply the 2026 write-up formatting to the exported thesis draft."""
import sys, docx
from docx.shared import Pt, Cm, RGBColor
from docx.enum.text import WD_ALIGN_PARAGRAPH, WD_BREAK
from docx.oxml.ns import qn
from docx.oxml import OxmlElement

src, out = sys.argv[1], sys.argv[2]
d = docx.Document(src)

def font(style_or_run, size=None, bold=None, italic=None):
    f = style_or_run.font
    f.name = "Times New Roman"
    rpr = style_or_run.element.get_or_add_rPr()
    rf = rpr.find(qn("w:rFonts"))
    if rf is None:
        rf = OxmlElement("w:rFonts"); rpr.append(rf)
    for a in ("w:ascii", "w:hAnsi", "w:cs", "w:eastAsia"):
        rf.set(qn(a), "Times New Roman")
    for a in ("w:asciiTheme", "w:hAnsiTheme", "w:cstheme", "w:eastAsiaTheme"):
        if rf.get(qn(a)) is not None: del rf.attrib[qn(a)]
    if size: f.size = Pt(size)
    if bold is not None: f.bold = bold
    if italic is not None: f.italic = italic
    f.color.rgb = RGBColor(0, 0, 0)

# --- styles -------------------------------------------------------------
n = d.styles["Normal"]; font(n, 12)
n.paragraph_format.line_spacing = 1.5
n.paragraph_format.alignment = WD_ALIGN_PARAGRAPH.JUSTIFY
n.paragraph_format.space_after = Pt(6)
for name in ("List Paragraph",):
    s = d.styles[name]; font(s, 12); s.paragraph_format.line_spacing = 1.5
font(d.styles["Heading 1"], 16, bold=True, italic=False)
font(d.styles["Heading 2"], 14, bold=True, italic=False)
font(d.styles["Heading 3"], 11, bold=False, italic=True)
for h in ("Heading 1", "Heading 2", "Heading 3"):
    pf = d.styles[h].paragraph_format
    pf.alignment = WD_ALIGN_PARAGRAPH.LEFT; pf.line_spacing = 1.5
    pf.space_before = Pt(12); pf.space_after = Pt(6)

# runs: force TNR everywhere (keeps bold/italic as exported)
def all_paragraphs():
    for p in d.paragraphs: yield p
    for t in d.tables:
        for row in t.rows:
            for c in row.cells:
                for p in c.paragraphs: yield p
for p in all_paragraphs():
    for r in p.runs:
        font(r)
for t in d.tables:
    for row in t.rows:
        for c in row.cells:
            for p in c.paragraphs:
                p.paragraph_format.alignment = WD_ALIGN_PARAGRAPH.LEFT
                p.paragraph_format.line_spacing = 1.0
                for r in p.runs: r.font.size = Pt(10)

# --- page setup + page numbers -------------------------------------------
for sec in d.sections:
    sec.page_height, sec.page_width = Cm(29.7), Cm(21.0)
    for m in ("top_margin", "bottom_margin", "left_margin", "right_margin"):
        setattr(sec, m, Cm(2.54))
    fp = sec.footer.paragraphs[0] if sec.footer.paragraphs else sec.footer.add_paragraph()
    fp.alignment = WD_ALIGN_PARAGRAPH.CENTER
    def fld(par, instr):
        r = par.add_run(); font(r, 12)
        a = OxmlElement("w:fldChar"); a.set(qn("w:fldCharType"), "begin")
        b = OxmlElement("w:instrText"); b.set(qn("xml:space"), "preserve"); b.text = instr
        c = OxmlElement("w:fldChar"); c.set(qn("w:fldCharType"), "end")
        r._r.extend([a, b, c])
    fld(fp, "PAGE")

def delete(p):
    p._element.getparent().remove(p._element)

P = d.paragraphs
# drop the doc's working header (title, byline, draft note) and the "Title page" heading
for p in P[:4]:
    delete(p)
P = d.paragraphs
# title page: centred, single block, page break after
i = 0
while P[i].style.name != "Heading 2":
    p = P[i]; p.paragraph_format.alignment = WD_ALIGN_PARAGRAPH.CENTER
    p.paragraph_format.space_after = Pt(18)
    i += 1
P[0].paragraph_format.space_before = Pt(150)
for r in P[0].runs: r.font.size = Pt(16)

# table of contents: replace the hand list with a Word TOC field
P = d.paragraphs
toc_h = next(p for p in P if p.style.name == "Heading 2" and p.text.strip() == "Table of contents")
toc_h.text = "Table of Contents"
j = P.index(toc_h) + 1
while P[j].style.name != "Heading 2":
    delete(P[j]); j += 1
newp = OxmlElement("w:p"); toc_h._p.addnext(newp)
tp = docx.text.paragraph.Paragraph(newp, toc_h._parent)
r = tp.add_run(); font(r, 12)
for typ, txt in (("begin", None), (None, 'TOC \\o "2-3" \\h \\z \\u'), ("separate", None)):
    if typ:
        e = OxmlElement("w:fldChar"); e.set(qn("w:fldCharType"), typ)
    else:
        e = OxmlElement("w:instrText"); e.set(qn("xml:space"), "preserve"); e.text = txt
    r._r.append(e)
r2 = tp.add_run("Right-click here and choose Update Field to build the table of contents with page numbers."); font(r2, 12, italic=True)
r3 = tp.add_run(); e = OxmlElement("w:fldChar"); e.set(qn("w:fldCharType"), "end"); r3._r.append(e)

# each main section starts on a new page
for p in d.paragraphs:
    if p.style.name == "Heading 2":
        p.paragraph_format.page_break_before = True

# references: APA hanging indent, left-aligned
P = d.paragraphs
inref = False
for p in P:
    if p.style.name == "Heading 2":
        inref = p.text.strip() == "References"; continue
    if inref and p.text.strip():
        pf = p.paragraph_format
        pf.left_indent = Cm(1.27); pf.first_line_indent = Cm(-1.27)
        pf.alignment = WD_ALIGN_PARAGRAPH.LEFT

d.core_properties.title = "Investigating Protein Splicing in Mycobacterium tuberculosis Using Proteomics Data"
d.core_properties.author = "Alphonse Hounkpe"
# British English (en-GB) as the proofing language, so Word's spell-check
# follows British spelling (requested 16 Oct 2026)
_rpr = d.styles.element.find(qn("w:docDefaults")).find(qn("w:rPrDefault")).find(qn("w:rPr"))
for _old in _rpr.findall(qn("w:lang")):
    _rpr.remove(_old)
_lang = OxmlElement("w:lang"); _lang.set(qn("w:val"), "en-GB"); _lang.set(qn("w:eastAsia"), "en-GB"); _rpr.append(_lang)
for _el in d.element.body.iter(qn("w:lang")):
    _el.set(qn("w:val"), "en-GB")
d.save(out)
print("saved", out)

# --- figures: replace the "[image: ...]" placeholders the export leaves -----------
import os
from docx.shared import Inches
FIG = {"The four steps of protein splicing": "fig1_splicing_mechanism.png",
       "The three inteins to scale": "fig2_inteins_to_scale.png",
       "Analysis workflow": "fig3_workflow.png",
       "Fragment ions across the SufB junction": "fig4_sufb_junction_ions.png",
       "Intein depletion across datasets": "fig5_depletion.png",
       "Rifampicin versus DMSO": "fig6_rif_vs_dmso.png"}
d = docx.Document(out)
figdir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "figures")
n_img = 0
for p in d.paragraphs:
    t = p.text.strip()
    if t.startswith("[image:"):
        key = t[len("[image:"):].rstrip("]").strip()
        for r in list(p.runs):
            r._r.getparent().remove(r._r)
        p.paragraph_format.alignment = WD_ALIGN_PARAGRAPH.CENTER
        p.paragraph_format.line_spacing = 1.0
        p.add_run().add_picture(os.path.join(figdir, FIG[key]), width=Inches(6.0))
        n_img += 1
    # APA figure/table number lines and notes: left aligned, keep with next
    if (t.startswith("Figure ") or t.startswith("Table ")) and len(t) <= 9:
        p.paragraph_format.alignment = WD_ALIGN_PARAGRAPH.LEFT
        p.paragraph_format.keep_with_next = True
P = d.paragraphs
for i, p in enumerate(P[:-1]):
    if (p.text.startswith("Figure ") or p.text.startswith("Table ")) and len(p.text.strip()) <= 9:
        P[i + 1].paragraph_format.keep_with_next = True   # title line stays with the figure/table
    if p._p.xpath(".//a:blip"):
        # figures: picture first, then "Figure N", title and note UNDER it
        # (requested 15 Oct 2026); the picture stays on the page with its caption
        p.paragraph_format.keep_with_next = True
        p.paragraph_format.alignment = WD_ALIGN_PARAGRAPH.CENTER
d.save(out)
print("figures inserted:", n_img)
