# =============================================================================
# ctrl.py -- pseudo-intein control for the depletion test (PXD009315)
#
# THESIS: Methods 3.7 (control); Results 4.3 (667 control proteins,
#         empirical P = .003 for RecA and .010 for SufB).
#
# QUESTION
#   A P value from the binomial formula assumes the model is right. Is it?
#   If the depletion test easily reports "too few peptides in the middle of
#   a protein" even when nothing has been cut out, RecA and SufB looking
#   depleted would mean little.
#
# HOW
#   1. Read the submitters' 16 MS-GF+ result files (0686Mass_*.mzid) and
#      collect, for every protein, the distinct identified peptide positions
#      (rank-1, non-decoy, MS-GF QValue <= 0.01).
#   2. Run the depletion test (same digest as depl.py, <= 1 missed cleavage)
#      on the real RecA and SufB inteins.
#   3. Choose CONTROL proteins: every other protein of >= 300 residues with
#      >= 6 distinct identified peptides (RecA, SufB, DnaB excluded;
#      contaminant entries NOT removed).
#   4. Give each control a "pseudo-intein" at the SAME RELATIVE POSITION as
#      the real intein (RecA: 31.9-87.5% of length; SufB: 29.9-72.2%),
#      boundaries scaled to the control's length and rounded down.
#   5. Run exactly the same test on each control. Nothing is excised in these
#      proteins, so their results show how the test behaves by chance.
#   6. Report: how many controls reach p <= 0.05; how many are at least as
#      extreme as the real protein (r); how many have zero peptides in the
#      pseudo-intein. Empirical probability for the thesis = (r + 1) / (m + 1),
#      where m = number of controls tested (computed by hand from the output).
#
# IN PLAIN WORDS
#   We pretend hundreds of ordinary proteins have an intein in the same place
#   and run the same test on them. If ordinary proteins rarely look as extreme
#   as RecA or SufB, the result is not a quirk of the test.
#
# INPUT   run inside the folder holding the PXD009315 deposit files:
#           *.pep.all.fasta      the submitters' Ensembl FASTA (exactly one)
#           0686Mass_*.mzid      the submitters' 16 MS-GF+ result files
# OUTPUT  two lines per protein, printed; progress messages on stderr.
# RUN     python3 ctrl.py
#
# CAVEATS
#   * PROLINE RULE NOT APPLIED; <= 1 missed cleavage (j < i + 3 below).
#   * A protein is skipped (returns None) if it has no observed peptides, or
#     its digest predicts no peptide inside or none outside the region.
#   * Before the control result was accepted, this script was required to
#     reproduce the real proteins' binomial values from depl.py exactly.
#   * Most controls carry more than 6 peptides, so the control is not
#     matched to RecA's small count. A count-matched control is a [TODO].
#
# NEEDS   Python 3.8+ standard library only.
# =============================================================================
import sys, re, glob, xml.etree.ElementTree as ET
from math import comb
# seq: {Ensembl ID: sequence}; obs: {Ensembl ID: set of (start, end)} 1-based
seq, name, obs = {}, None, {}
for l in open(glob.glob("*.pep.all.fasta")[0]):
    l = l.strip()
    if l[:1] == ">": name = l[1:].split()[0]; seq[name] = ""
    elif name: seq[name] += l
# --- read every result file, streaming the XML to keep memory low ----------
for path in sorted(glob.glob("0686Mass_*.mzid")):
    db, ev = {}, {}
    for _, el in ET.iterparse(path):
        t = el.tag.rsplit("}", 1)[-1]          # tag name without namespace
        # DBSequence: internal ID -> protein accession
        if t == "DBSequence": db[el.get("id")] = el.get("accession")
        # PeptideEvidence: where a peptide sits on a protein (targets only)
        elif t == "PeptideEvidence" and el.get("isDecoy") != "true":
            ev[el.get("id")] = (db[el.get("dBSequence_ref")], int(el.get("start")), int(el.get("end")))
        # SpectrumIdentificationItem: one PSM. Keep rank-1 at q <= 0.01.
        elif t == "SpectrumIdentificationItem" and el.get("rank") == "1":
            q = [c.get("value") for c in el if c.get("name") == "MS-GF:QValue"]
            if q and float(q[0]) <= 0.01:
                for c in el:
                    r = c.get("peptideEvidence_ref")
                    # record the distinct (start, end) position on that protein
                    if r in ev: obs.setdefault(ev[r][0], set()).add(ev[r][1:])
        # free memory for elements already used
        if t in ("DBSequence", "PeptideEvidence", "SpectrumIdentificationResult"): el.clear()
    sys.stderr.write(path + " read\n")
# Depletion test for region a..b (1-based) of protein acc.
# Returns (k intein peptides, n peptides, expected share p, binomial P),
# or None if the test cannot be computed.
def test(acc, a, b):
    s = seq[acc]
    # cut positions, no proline rule (same as depl.py)
    cuts = [0] + [m.end() for m in re.finditer("[KR]", s) if m.end() < len(s)] + [len(s)]
    # c[0] = predicted peptides outside the region, c[1] = inside;
    # k[0] / k[1] = the same for observed peptides. Spanning peptides skipped.
    c, k = [0, 0], [0, 0]
    for i in range(len(cuts) - 1):
        # j < i + 3  ->  at most 1 missed cleavage
        for j in range(i + 1, min(i + 3, len(cuts))):
            st, en = cuts[i] + 1, cuts[j]
            if 6 <= en - st + 1 <= 40 and (en < a or st > b or (st >= a and en <= b)): c[st >= a and en <= b] += 1
    for st, en in obs.get(acc, ()):
        if en < a or st > b or (st >= a and en <= b): k[st >= a and en <= b] += 1
    n = sum(k)
    if not n or not c[0] or not c[1]: return None
    p = c[1] / sum(c)
    # one-sided exact binomial tail P(X <= k_inside)
    return k[1], n, p, sum(comb(n, x) * p**x * (1 - p)**(n - x) for x in range(k[1] + 1))
# Ensembl IDs of RecA, SufB, DnaB: never used as controls
skip = ("CCP45535", "CCP44220", "CCP42780")
for nm, acc, a, b in (("RecA", skip[0], 252, 691), ("SufB", skip[1], 253, 611)):
    k, n, p, pv = test(acc, a, b)
    L = len(seq[acc])
    # controls: pseudo-intein at the same relative position (a/L .. b/L),
    # scaled to each control's length; >= 300 aa and >= 6 distinct peptides
    res = [r for r in (test(x, int(a / L * len(seq[x])) + 1, int(b / L * len(seq[x])))
           for x in obs if x in seq and x not in skip and len(seq[x]) >= 300 and len(obs[x]) >= 6) if r]
    m = len(res)
    print("%s: %d of %d in intein, expected share %.3f, p = %.4g" % (nm, k, n, p, pv))
    print("  controls: %d proteins; p <= 0.05 in %d (%.1f%%); p <= %s's p in %d (%.1f%%); none in pseudo-intein: %d (%.1f%%)" % (
        m, sum(r[3] <= 0.05 for r in res), 100 * sum(r[3] <= 0.05 for r in res) / m, nm,
        sum(r[3] <= pv for r in res), 100 * sum(r[3] <= pv for r in res) / m,
        sum(r[0] == 0 for r in res), 100 * sum(r[0] == 0 for r in res) / m))
