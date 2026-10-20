"""Thesis figures 1-6. All numbers are typed from the analysis outputs recorded in
claude/TB-rif-DIA-results.md and the thesis Results; nothing is simulated."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Rectangle
import numpy as np
import re

plt.rcParams.update({
    "font.family": "Liberation Serif", "font.size": 10, "axes.titlesize": 11,
    "axes.labelsize": 10, "axes.spines.top": False, "axes.spines.right": False,
    "axes.edgecolor": "#52514e", "axes.linewidth": 0.8, "xtick.color": "#52514e",
    "ytick.color": "#52514e", "savefig.dpi": 300, "mathtext.fontset": "custom",
    "mathtext.rm": "Liberation Serif", "mathtext.it": "Liberation Serif:italic", "savefig.bbox": "tight",
})
NX, INT, CX = "#2a78d6", "#eb6834", "#1baf7a"     # N-extein, intein, C-extein
RIF, DMSO, POOL, GREY = "#2a78d6", "#eb6834", "#6250d6", "#9a9890"
INK, INK2 = "#0b0b0b", "#52514e"
OUT = "figures/"


def save(fig, name):
    """Save each figure twice: PNG (picture) and SVG (editable).
    In the SVG, text stays as real text (svg.fonttype none), and the font is
    named Times New Roman, so in Word (right-click > Convert to Shape),
    PowerPoint or Inkscape every label, line and colour can be edited."""
    fig.savefig(OUT + name + ".png")
    with plt.rc_context({"svg.fonttype": "none"}):
        fig.savefig(OUT + name + ".svg")
    svg = open(OUT + name + ".svg").read().replace("Liberation Serif", "Times New Roman")
    # Word's SVG renderer mis-draws quadratic curves (Q), which showed as stray
    # lines across the rounded boxes; rewrite every Q as the equivalent cubic (C)
    svg = re.sub(r'(<path d=")([^"]*)(")', lambda m: m.group(1) + q_to_c(m.group(2)) + m.group(3), svg)
    open(OUT + name + ".svg", "w").write(svg)


def q_to_c(d):
    """Convert absolute quadratic Bezier segments (Q) in an SVG path to cubic (C)."""
    tok = re.findall(r"[MLQCZz]|-?[0-9.]+(?:e-?[0-9]+)?", d)
    out, i, cur, start = [], 0, (0.0, 0.0), (0.0, 0.0)
    cmd = None
    while i < len(tok):
        t = tok[i]
        if t in "MLQCZz":
            cmd = t; i += 1
            if t in "Zz":
                out.append("z"); cur = start
            continue
        if cmd in ("M", "L"):
            x, y = float(tok[i]), float(tok[i + 1]); i += 2
            out.append("%s %s %s" % (cmd, tok[i - 2], tok[i - 1]))
            cur = (x, y)
            if cmd == "M":
                start = cur
        elif cmd == "Q":
            qx, qy, x, y = (float(v) for v in tok[i:i + 4]); i += 4
            c1 = (cur[0] + 2 / 3 * (qx - cur[0]), cur[1] + 2 / 3 * (qy - cur[1]))
            c2 = (x + 2 / 3 * (qx - x), y + 2 / 3 * (qy - y))
            out.append("C %.4f %.4f %.4f %.4f %.4f %.4f" % (c1 + c2 + (x, y)))
            cur = (x, y)
        elif cmd == "C":
            v = tok[i:i + 6]; i += 6
            out.append("C " + " ".join(v)); cur = (float(v[4]), float(v[5]))
        else:
            i += 1
    return " ".join(out)


def seg(ax, x, y, w, h, color, label=None, tc="white", fs=9):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0,rounding_size=0.06",
                                fc=color, ec="white", lw=1.5))
    if label:
        ax.text(x + w / 2, y + h / 2, label, ha="center", va="center", color=tc, fontsize=fs)


# ---------------------------------------------------------------- Figure 1
def fig1():
    fig, ax = plt.subplots(figsize=(7.2, 6.2))
    ax.set_xlim(0, 10); ax.set_ylim(0, 10.6); ax.axis("off")
    rows = [
        (9.4, "Precursor", "N-extein–C(O)–NH–Cys1 (intein) … Asn–NH–Cys+1 (C-extein)"),
        (7.3, "Step 1  N→S acyl shift", "Cys1 thiol attacks the preceding peptide bond: linear thioester"),
        (5.2, "Step 2  Transesterification", "Cys+1 of the C-extein takes over the N-extein: branched intermediate"),
        (3.1, "Step 3  Asn cyclisation", "terminal Asn forms a succinimide, releasing the intein"),
        (1.0, "Step 4  S→N acyl shift", "thioester rearranges to a native peptide bond: spliced protein"),
    ]
    for y, title, note in rows:
        ax.text(0.0, y + 0.62, title, fontsize=10.5, fontweight="bold", color=INK)
        ax.text(0.0, y - 0.62, note, fontsize=8.5, color=INK2, style="italic")
    # precursor
    y = 9.4
    seg(ax, 0.2, y - 0.3, 2.4, 0.6, NX, "N-extein")
    seg(ax, 2.6, y - 0.3, 4.6, 0.6, INT, "Intein")
    seg(ax, 7.2, y - 0.3, 2.4, 0.6, CX, "C-extein")
    ax.text(2.75, y + 0.36, "Cys1", ha="center", fontsize=8, color=INK)
    ax.text(7.2, y + 0.36, "Asn | Cys+1", ha="center", fontsize=8, color=INK)
    # step 1: thioester (S between N-extein and intein)
    y = 7.3
    seg(ax, 0.2, y - 0.3, 2.4, 0.6, NX, "N-extein")
    ax.text(2.73, y, "S", ha="center", va="center", fontsize=9, fontweight="bold", color=INK)
    seg(ax, 2.85, y - 0.3, 4.35, 0.6, INT, "Intein")
    seg(ax, 7.2, y - 0.3, 2.4, 0.6, CX, "C-extein")
    # step 2: branched
    y = 5.2
    seg(ax, 2.6, y - 0.3, 4.6, 0.6, INT, "Intein")
    seg(ax, 7.2, y - 0.3, 2.4, 0.6, CX, "C-extein")
    seg(ax, 4.9, y + 0.35, 2.4, 0.45, NX, "N-extein", fs=8)
    ax.text(7.45, y + 0.57, "S", ha="center", va="center", fontsize=9, fontweight="bold", color=INK)
    # step 3: excised intein + ligated-by-thioester exteins
    y = 3.1
    seg(ax, 0.2, y - 0.3, 3.6, 0.6, INT, "Excised intein")
    ax.text(3.85, y, "succinimide", ha="left", va="center", fontsize=7.5, color=INK2)
    seg(ax, 5.0, y - 0.3, 2.4, 0.6, NX, "N-extein")
    ax.text(7.52, y, "S", ha="center", va="center", fontsize=9, fontweight="bold", color=INK)
    seg(ax, 7.65, y - 0.3, 1.95, 0.6, CX, "C-extein")
    # step 4: mature
    y = 1.0
    seg(ax, 5.0, y - 0.3, 2.3, 0.6, NX, "N-extein")
    seg(ax, 7.3, y - 0.3, 2.3, 0.6, CX, "C-extein")
    ax.annotate("new peptide bond", xy=(7.3, y + 0.3), xytext=(4.6, y + 0.75), fontsize=8,
                color=INK, ha="right", va="center", arrowprops=dict(arrowstyle="-", color=INK2, lw=0.8))
    for y in (8.35, 6.25, 4.15, 2.05):
        ax.add_patch(FancyArrowPatch((8.4, y + 0.38), (8.4, y - 0.38), arrowstyle="-|>",
                                     mutation_scale=10, color=INK2, lw=1))
    save(fig, "fig1_splicing_mechanism"); plt.close(fig)


# ---------------------------------------------------------------- Figure 2
def fig2():
    P = [("RecA", "P9WHJ3", 790, 252, 691, "K251 | C692",
          "type 3: CSPPFK (boundary-start); type 2 needs a missed cleavage"),
         ("SufB", "P9WFP7", 846, 253, 611, "G252 | C612",
          "type 2: TLIIADEGSYVHYVEG|CTAPIYK (fully tryptic)"),
         ("DnaB", "P9WMR3", 874, 400, 815, "G399 | S816",
          "type 2: ESG|SLEQDADVVILLHRPDAFDR (fully tryptic)")]
    fig, axs = plt.subplots(len(P), 1, figsize=(7.2, 4.3), sharex=True)
    for ax, (name, acc, L, a, b, junc, diag) in zip(axs, P):
        ax.set_ylim(-1.2, 1.4); ax.axis("off")
        seg(ax, 0, -0.25, a - 1, 0.5, NX)
        seg(ax, a - 1, -0.25, b - a + 1, 0.5, INT, "intein %d–%d (%d aa)" % (a, b, b - a + 1))
        seg(ax, b, -0.25, L - b, 0.5, CX)
        ax.text(-12, 0, "%s\n%s, %d aa" % (name, acc, L), ha="right", va="center", fontsize=9.5,
                color=INK)
        ax.text(a - 1, 0.45, "C%d" % a, ha="center", fontsize=8, color=INK2)
        ax.text(b, 0.45, "N%d" % b, ha="center", fontsize=8, color=INK2)
        ax.text(L, 0.45, "%d" % L, ha="right", fontsize=8, color=INK2)
        ax.text(0, 0.45, "1", ha="left", fontsize=8, color=INK2)
        ax.text(L / 2, -0.75, "junction after splicing: %s;  %s" % (junc, diag), ha="center",
                fontsize=8, color=INK2, style="italic")
    axs[0].set_xlim(-170, 900)
    h = [Rectangle((0, 0), 1, 1, fc=c) for c in (NX, INT, CX)]
    fig.legend(h, ["N-extein", "intein", "C-extein"], loc="lower center", ncol=3, frameon=False,
               fontsize=9, bbox_to_anchor=(0.55, -0.02))
    save(fig, "fig2_inteins_to_scale"); plt.close(fig)


# ---------------------------------------------------------------- Figure 3
def fig3():
    fig, ax = plt.subplots(figsize=(7.2, 5.6))
    ax.set_xlim(0, 12); ax.set_ylim(0, 9.6); ax.axis("off")

    def box(x, y, w, h, t, c="#eef3fb", ec=NX, fs=8.0):
        ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.02,rounding_size=0.12",
                                    fc=c, ec=ec, lw=1))
        ax.text(x + w / 2, y + h / 2, t, ha="center", va="center", fontsize=fs, color=INK,
                linespacing=1.3)

    def arr(x1, y1, x2, y2):
        ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle="-|>", mutation_scale=9,
                                     color=INK2, lw=0.9))
    box(2.6, 8.5, 6.8, 0.85, "H37Rv reference proteome (UniProt UP000001584)\nand intein coordinates of RecA, SufB and DnaB",
        c="#f4f3ef", ec=INK2)
    box(2.6, 7.1, 6.8, 0.85, "In silico digestion of precursor and spliced forms:\ndiagnostic peptide types 1–4 and expected intein share",
        c="#f4f3ef", ec=INK2)
    arr(6, 8.5, 6, 7.95)
    cols = [(0.1, "Dataset 1: Progenio\nSAWC507 / SAWC5527, DDA\nMaxQuant protein-group export",
             "repaired export;\ntargeted panel, cleavage-agnostic\nand split-match searches"),
            (4.1, "Dataset 2: PXD009315\nQ Exactive, DDA, 16 runs",
             "MS-GF+ semi-tryptic search\nwith spliced sequences;\nCasanovo de novo (2 runs)"),
            (8.1, "Dataset 3: rifampicin vs DMSO\ntimsTOF Pro 2, diaPASEF, 10 runs",
             "diaTracer, then MSFragger\nsemi-tryptic search with\nspliced sequences; DIA-NN")]
    for x, top, mid in cols:
        arr(6, 7.1, x + 1.9, 6.3)
        box(x, 5.3, 3.8, 1.0, top)
        arr(x + 1.9, 5.3, x + 1.9, 4.9)
        box(x, 3.8, 3.8, 1.1, mid)
        arr(x + 1.9, 3.8, x + 1.9, 3.35)
    box(0.1, 2.05, 11.8, 1.3, "Intein analysis, the same for every dataset: peptides by region on RecA, SufB and DnaB;\n"
        "diagnostic types 1–4; non-tryptic termini against background; junction matches at any score;\n"
        "depletion test (exact binomial) with pseudo-intein control",
        c="#fdf0ea", ec=INT)
    arr(6, 2.05, 6, 1.6)
    box(0.1, 0.25, 11.8, 1.35, "Evidence standard: fragment ions bracketing the bond in question.\n"
        "Dataset 3 also: DIA-NN q-values and retention time, decoys of identical composition,\n"
        "junction-to-extein ratio, rifampicin vs DMSO. Datasets compared as patterns, never as numbers.",
        c="#eaf7f1", ec=CX)
    save(fig, "fig3_workflow"); plt.close(fig)


# ---------------------------------------------------------------- Figure 5 (depletion)
def fig4():
    # (dataset, protein, k, n, expected share, p)
    D = [("Dataset 1", "RecA", 2, 12, .587, ".0037"), ("Dataset 1", "SufB", 4, 27, .450, ".00095"),
         ("Dataset 2", "RecA", 0, 6, .608, ".0036"), ("Dataset 2", "SufB", 6, 23, .474, ".031"),
         ("Dataset 2", "DnaB", 0, 2, .537, ".21"),
         ("Dataset 3", "RecA", 0, 11, .587, "5.9×10$^{-5}$"), ("Dataset 3", "SufB", 0, 12, .450, "7.7×10$^{-4}$"),
         ("Dataset 3", "DnaB", 1, 3, .559, ".41")]
    fig, axs = plt.subplots(1, 3, figsize=(7.2, 3.3), sharey=True)
    for ax, ds in zip(axs, ["Dataset 1", "Dataset 2", "Dataset 3"]):
        rows = [r for r in D if r[0] == ds]
        x = np.arange(len(rows))
        ax.bar(x - 0.19, [r[4] * 100 for r in rows], 0.36, color="#d6d4cd", label="expected (digest)")
        ax.bar(x + 0.19, [100 * r[2] / r[3] for r in rows], 0.36, color=INT, label="observed")
        for i, r in enumerate(rows):
            ax.text(x[i] + 0.19, 100 * r[2] / r[3] + 2, "%d/%d" % (r[2], r[3]), ha="center", fontsize=7.5, color=INK)
            ax.text(x[i], -11, "$p$ =\n" + r[5], ha="center", va="top", fontsize=7.5, color=INK2, linespacing=1.2)
        ax.set_xticks(x); ax.set_xticklabels([r[1] for r in rows])
        ax.set_title(ds, fontsize=10)
        ax.set_ylim(0, 72); ax.grid(axis="y", color="#e6e4dd", lw=0.6); ax.set_axisbelow(True)
        ax.tick_params(axis="x", length=0)
    if "Dataset 1":
        axs[0].text(2, 3, "DnaB not\ndetected", ha="center", fontsize=7.5, color=INK2, style="italic")
        axs[0].set_xlim(-0.6, 2.6)
    axs[0].set_ylabel("Intein share of distinct peptides (%)")
    axs[0].legend(frameon=False, fontsize=8, loc="upper right")
    fig.subplots_adjust(bottom=0.24, wspace=0.12)
    save(fig, "fig5_depletion"); plt.close(fig)


# ---------------------------------------------------------------- Figure 4 (junction ions)
def fig5():
    pep = "TLIIADEGSYVHYVEGCTAPIYK"; n = len(pep)
    spectra = [("Pool_1", [3, 5, 9], [4, 8, 9, 13, 16, 18, 19, 20], 11, 5),
               ("DMSO_3", [5, 14, 19], [4, 16, 18, 19, 20], 8, 4),
               ("Rif_2", [3, 7, 10, 14], [5, 12, 16, 17, 19, 20, 21], 11, 4),
               ("Rif_1", [21], [4, 8, 9, 10, 16, 19, 20, 21], 9, 4),
               ("Rif_4", [4, 17], [18, 20, 21], 5, 2)]
    fig = plt.figure(figsize=(7.2, 3.9))
    gs = fig.add_gridspec(1, 2, width_ratios=[3.2, 1], wspace=0.32)
    ax = fig.add_subplot(gs[0]); ax2 = fig.add_subplot(gs[1])
    ns = len(spectra)
    ax.set_xlim(-0.2, n + 0.2); ax.set_ylim(-0.4, ns * 1.3 + 1.6); ax.axis("off")
    ytop = ns * 1.3 + 0.35
    for i, c in enumerate(pep):
        ax.add_patch(Rectangle((i + 0.04, ytop), 0.92, 0.8, fc=NX if i < 16 else CX, ec="white"))
        ax.text(i + 0.5, ytop + 0.4, c, ha="center", va="center", color="white", fontsize=8.5,
                fontweight="bold")
    ax.plot([16, 16], [0, ytop], color=INK, lw=0.9, ls=(0, (3, 2)))
    ax.text(16, ytop + 0.95, "junction G252 | C612", ha="center", fontsize=8, color=INK)
    for k, (name, b, y, _, _) in enumerate(spectra):
        yy = (ns - 1 - k) * 1.3 + 0.2
        ax.add_patch(Rectangle((0, yy - 0.05), n, 0.9, fc="#f4f3ef", ec="none"))
        ax.text(-0.4, yy + 0.4, name, ha="right", va="center", fontsize=8.5, color=INK)
        for i in range(1, n):
            ax.plot([i, i], [yy - 0.05, yy + 0.85], color="white", lw=0.6)
        for j in b:   # b_j: cleavage after residue j, N-terminal fragment
            jc = j >= 17
            ax.plot([j - 0.42, j], [yy + 0.62, yy + 0.62], color=INK if jc else NX, lw=2.2,
                    solid_capstyle="butt")
            ax.plot([j, j], [yy + 0.45, yy + 0.8], color=INK if jc else NX, lw=1.6)
        for j in y:   # y_j: cleavage after residue n-j, C-terminal fragment
            x = n - j; jc = j >= 8
            ax.plot([x, x + 0.42], [yy + 0.18, yy + 0.18], color=INK if jc else CX, lw=2.2,
                    solid_capstyle="butt")
            ax.plot([x, x], [yy + 0.0, yy + 0.35], color=INK if jc else CX, lw=1.6)
    yv = np.arange(ns)[::-1]
    ax2.barh(yv + 0.18, [s[3] for s in spectra], 0.34, color=INK, label="real sequence")
    ax2.barh(yv - 0.18, [s[4] for s in spectra], 0.34, color="#c9c7bf", label="best decoy")
    ax2.set_yticks(yv); ax2.set_yticklabels([s[0] for s in spectra], fontsize=8.5)
    ax2.set_xlabel("b and y ions matched"); ax2.set_xlim(0, 12)
    ax2.legend(frameon=False, fontsize=7.5, loc="upper center", bbox_to_anchor=(0.5, 1.2), ncol=1)
    ax2.grid(axis="x", color="#e6e4dd", lw=0.6); ax2.set_axisbelow(True)
    ax2.tick_params(axis="y", length=0)
    save(fig, "fig4_sufb_junction_ions"); plt.close(fig)


# ---------------------------------------------------------------- Figure 6
def fig6():
    ratio = {"Rif": [0.0509, 0.0553, 0.0623, 0.0531], "DMSO": [0.0526, 0.0305, 0.0425, 0.0567]}
    dnab = {"Rif": [0, 0, 897, 0], "DMSO": [720, 1315, 2211, 1547]}
    fig, axs = plt.subplots(1, 2, figsize=(7.2, 3.1))
    for ax, data, yl, title, p in [
        (axs[0], ratio, "Junction peptide ÷ SufB extein peptides", "SufB junction-to-extein ratio", "Welch's $t$, $p$ = .20"),
        (axs[1], dnab, "DnaB signal (precursor quantity)", "DnaB amount (N-extein peptides)", "Mann–Whitney $U$, $p$ = .057")]:
        for i, (g, c) in enumerate([("DMSO", DMSO), ("Rif", RIF)]):
            v = np.array(data[g]); xj = i + np.array([-0.09, -0.03, 0.03, 0.09])
            ax.scatter(xj, v, s=42, color=c, edgecolor="white", lw=1.2, zorder=3)
            ax.plot([i - 0.22, i + 0.22], [v.mean()] * 2, color=INK, lw=1.4)
        ax.set_xticks([0, 1]); ax.set_xticklabels(["DMSO (n = 4)", "rifampicin (n = 4)"])
        ax.set_xlim(-0.6, 1.6); ax.set_ylabel(yl); ax.set_title(title, fontsize=10)
        top = max(max(data["Rif"]), max(data["DMSO"])) * 1.25
        ax.set_ylim(-0.03 * top, top); ax.grid(axis="y", color="#e6e4dd", lw=0.6); ax.set_axisbelow(True)
        ax.text(0.5, 0.98, p, transform=ax.transAxes, ha="center", va="top",
                fontsize=8, color=INK2)
        ax.tick_params(axis="x", length=0)
    axs[1].text(1.0, 380, "3 of 4 runs not detected", fontsize=7.5, color=INK2, ha="center")
    fig.subplots_adjust(wspace=0.38)
    save(fig, "fig6_rif_vs_dmso"); plt.close(fig)


for f in (fig1, fig2, fig3, fig4, fig5, fig6):
    f()
print("done")
