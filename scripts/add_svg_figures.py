"""add_svg_figures.py -- make the thesis figures editable in Word.

Word shows an image from its SVG copy when the picture carries an
<asvg:svgBlip> link (Office 2016+), keeping the PNG as a fallback for older
viewers. After this, in Word: right-click a figure > Convert to Shape, and every
label, line and colour can be edited.

Figures are matched to the pictures in document order (Figure 1 is the first
picture, and so on), and each match is checked by aspect ratio.

RUN  python3 add_svg_figures.py in.docx out.docx figures/
"""
import sys, os, re, zipfile
from PIL import Image
import io

src, out, figdir = sys.argv[1:4]
FIGS = ["fig1_splicing_mechanism", "fig2_inteins_to_scale", "fig3_workflow",
        "fig4_sufb_junction_ions", "fig5_depletion", "fig6_rif_vs_dmso"]
SVGNS = "http://schemas.microsoft.com/office/drawing/2016/SVG/main"

zin = zipfile.ZipFile(src)
files = {n: zin.read(n) for n in zin.namelist()}
doc = files["word/document.xml"].decode("utf8")
rels = files["word/_rels/document.xml.rels"].decode("utf8")
ct = files["[Content_Types].xml"].decode("utf8")

rid_target = dict(re.findall(r'Id="(rId\d+)"[^>]*Target="(media/[^"]+)"', rels))
rid_target.update({a: b for b, a in re.findall(r'Target="(media/[^"]+)"[^>]*Id="(rId\d+)"', rels)})
blips = re.findall(r'<a:blip r:embed="(rId\d+)"/>', doc)
assert len(blips) == len(FIGS), "expected 6 pictures, found %d" % len(blips)

for i, (rid, name) in enumerate(zip(blips, FIGS), 1):
    w, h = Image.open(io.BytesIO(files["word/" + rid_target[rid]])).size
    W, H = Image.open(os.path.join(figdir, name + ".png")).size
    assert abs(w / h - W / H) < 0.02, "picture %d does not look like %s" % (i, name)
    files["word/" + rid_target[rid]] = open(os.path.join(figdir, name + ".png"), "rb").read()
    files["word/media/%s.svg" % name] = open(os.path.join(figdir, name + ".svg"), "rb").read()
    new = "rIdSvg%d" % i
    rels = rels.replace("</Relationships>",
        '<Relationship Id="%s" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/image" '
        'Target="media/%s.svg"/></Relationships>' % (new, name))
    doc = doc.replace('<a:blip r:embed="%s"/>' % rid,
        '<a:blip r:embed="%s"><a:extLst><a:ext uri="{96DAC541-7B7A-43D3-8B79-37D633B846F1}">'
        '<asvg:svgBlip xmlns:asvg="%s" r:embed="%s"/></a:ext></a:extLst></a:blip>' % (rid, SVGNS, new))
    print("Figure %d -> %s.svg" % (i, name))

if 'Extension="svg"' not in ct:
    ct = ct.replace("<Default ", '<Default Extension="svg" ContentType="image/svg+xml"/><Default ', 1)
files["word/document.xml"] = doc.encode("utf8")
files["word/_rels/document.xml.rels"] = rels.encode("utf8")
files["[Content_Types].xml"] = ct.encode("utf8")
with zipfile.ZipFile(out, "w", zipfile.ZIP_DEFLATED) as z:
    for n, b in files.items():
        z.writestr(n, b)
print("saved", out)
