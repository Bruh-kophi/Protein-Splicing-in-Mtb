#!/usr/bin/env python3
# =============================================================================
# make_synthetic_fixture.py -- builds a SYNTHETIC test set for
# fragpipe_intein_peptides.py. NOTHING HERE IS REAL DATA.
#
# The protein sequences are random (fixed seed), merely placed under the
# RecA / SufB / DnaB accessions so the script finds them. Peptides are cut
# from them at known positions, so the correct answer for every one is known
# in advance. All files are named SYNTHETIC_* and every FASTA header says
# SYNTHETIC. Never mix these files with real data.
#
# Run:  python3 make_synthetic_fixture.py OUTDIR
# =============================================================================
import os
import random
import sys

out = sys.argv[1]
os.makedirs(out, exist_ok=True)
rng = random.Random(20261015)
AA = "ACDEFGHKLMNPQRSTVWY"          # no I, so I/L merging cannot confuse the test

INT = {"RecA": ("P9WHJ3", 790, 252, 691), "SufB": ("P9WFP7", 846, 253, 611),
       "DnaB": ("P9WMR3", 874, 400, 815)}
seqs = {}
for name, (acc, L, a, b) in INT.items():
    s = [rng.choice(AA) for _ in range(L)]
    # put the real boundary residues in place (C...HN|C / S), and the residue
    # before the junction (RecA K, SufB/DnaB G) so tryptic logic is realistic
    s[a - 1] = "C"; s[b - 2] = "H"; s[b - 1] = "N"
    s[b] = "S" if name == "DnaB" else "C"
    s[a - 2] = "K" if name == "RecA" else "G"
    seqs[name] = "".join(s)
ctrl = "".join(rng.choice(AA) for _ in range(500))

with open(os.path.join(out, "SYNTHETIC_db.fasta"), "w") as fh:
    for name, (acc, L, a, b) in INT.items():
        p = seqs[name]
        fh.write(">sp|%s|SYNTHETIC_%s SYNTHETIC TEST SEQUENCE - NOT REAL\n%s\n" % (acc, name, p))
        fh.write(">SPLICED_%s_SYNTHETIC spliced form SYNTHETIC TEST\n%s\n" % (name, p[:a - 1] + p[b:]))
    fh.write(">sp|SYN001|SYNTHETIC_CTRL SYNTHETIC TEST SEQUENCE - NOT REAL\n%s\n" % ctrl)


def cut(name, s, e):
    """1-based inclusive slice of a synthetic precursor, with flanks."""
    p = seqs[name]
    return p[s - 1:e], (p[s - 2] if s > 1 else "-"), (p[e] if e < len(p) else "-")


# (run, peptide spec) -> expected classification, recorded for the test
rows, expect = [], []
runs = ["10984_260506_0045_Rif_1_S1-A2", "10988_260506_0045_DMSO_1_S1-A5",
        "10989_260506_0045_Pool_2_S1-A10"]
plan = [
    (runs[0], "RecA", 100, 112, "N-extein", ""),
    (runs[0], "RecA", 300, 312, "intein", ""),
    (runs[0], "RecA", 680, 691, "intein", "type4 intein-end"),
    (runs[0], "RecA", 692, 703, "C-extein", "type3 boundary-start"),
    (runs[1], "RecA", 100, 112, "N-extein", ""),
    (runs[1], "RecA", 720, 731, "C-extein", ""),
    (runs[1], "SufB", 247, 259, "SPANS intein START (unspliced)", "type1 boundary-spanning"),
    (runs[1], "SufB", 253, 265, "intein", "4N-a intein-start"),
    (runs[2], "SufB", 400, 412, "intein", ""),      # Pool_2: excluded from pooled counts
]
for i, (run, name, s, e, region, flag) in enumerate(plan):
    pep, prev, nxt = cut(name, s, e)
    rows.append(dict(Spectrum="%s_diatracer.%d.%d.2" % (run, 100 + i, 100 + i), Peptide=pep,
                     **{"Modified Peptide": "", "Prev AA": prev, "Next AA": nxt, "Charge": "2",
                        "Probability": "0.99", "Protein Start": str(s), "Protein": name}))
    expect.append((run, name, pep, region, flag))

# a junction peptide for RecA (in the spliced form only), seen in a Rif run
sp = seqs["RecA"][:251] + seqs["RecA"][691:]
jpep = sp[245:257]                      # 6 N-side + 6 C-side residues
rows.append(dict(Spectrum="%s_diatracer.900.900.2" % runs[0], Peptide=jpep,
                 **{"Modified Peptide": "", "Prev AA": sp[244], "Next AA": sp[257], "Charge": "2",
                    "Probability": "0.97", "Protein Start": "246", "Protein": "SPLICED_RecA"}))
expect.append((runs[0], "RecA", jpep, "JUNCTION (spliced only)", "type2 junction-spanning"))

# background PSMs on the control protein: 8 tryptic-looking, 2 non-tryptic C-termini
for i in range(10):
    s = 10 + 20 * i
    pep = ctrl[s:s + 10]
    pep = pep[:-1] + ("K" if i < 8 else "G")
    rows.append(dict(Spectrum="%s_diatracer.%d.%d.2" % (runs[1], 500 + i, 500 + i), Peptide=pep,
                     **{"Modified Peptide": "", "Prev AA": "R", "Next AA": "A", "Charge": "2",
                        "Probability": "0.99", "Protein Start": str(s + 1), "Protein": "SYN001"}))

cols = ["Spectrum", "Peptide", "Modified Peptide", "Prev AA", "Next AA", "Charge",
        "Probability", "Protein Start", "Protein"]
with open(os.path.join(out, "SYNTHETIC_psm.tsv"), "w") as fh:
    fh.write("\t".join(cols) + "\n")
    for r in rows:
        fh.write("\t".join(r[c] for c in cols) + "\n")

# pepXML with one WEAK junction hit (DnaB) that is NOT in psm.tsv
dsp = seqs["DnaB"][:399] + seqs["DnaB"][815:]
weak = dsp[393:405]
with open(os.path.join(out, "SYNTHETIC_interact.pep.xml"), "w") as fh:
    fh.write('<?xml version="1.0"?><msms_pipeline_analysis><msms_run_summary>'
             '<spectrum_query spectrum="%s_diatracer.777.777.3" assumed_charge="3">'
             '<search_result><search_hit hit_rank="1" peptide="%s" protein="SPLICED_DnaB_SYNTHETIC">'
             '<search_score name="expect" value="2.5"/><analysis_result analysis="peptideprophet">'
             '<peptideprophet_result probability="0.12"/></analysis_result></search_hit>'
             '</search_result></spectrum_query></msms_run_summary></msms_pipeline_analysis>\n'
             % (runs[0], weak))

with open(os.path.join(out, "SYNTHETIC_expected.tsv"), "w") as fh:
    fh.write("run\tprotein\tpeptide\tregion\tflag\n")
    for e in expect:
        fh.write("\t".join(e) + "\n")
    fh.write("PEPXML\tDnaB\t%s\tJUNCTION (spliced only)\ttype2 junction-spanning\n" % weak)
print("synthetic fixture written to", out)
