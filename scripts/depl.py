# =============================================================================
# depl.py -- in silico tryptic digest + binomial depletion test (PXD009315)
#
# THESIS: Methods 3.3 (in silico digestion) and 3.7 (depletion test);
#         Results 4.3 (Table 4.4, Dataset 2 / PXD009315).
#
# QUESTION
#   If an intein is cut out of its host protein and the free intein is then
#   degraded, peptides from the intein should be RARER than peptides from the
#   exteins of the same protein. Are they?
#
# HOW
#   1. Read the protein sequences (Ensembl H37Rv, as deposited with PXD009315).
#   2. Digest RecA and SufB in silico with trypsin (cut after K or R) at
#      0, <=1 and <=2 missed cleavages, keeping peptides of 6-40 residues.
#   3. Classify each predicted peptide as N-extein, intein, C-extein, or
#      boundary-spanning (X). Spanning peptides belong to neither region and
#      are left out of the expectation.
#   4. Expected intein share  p = intein / (N-extein + intein + C-extein)
#   5. Observed: k intein peptides out of n distinct peptides (hard-coded
#      below, taken from mzid_intein_peptides.py on all 16 runs).
#   6. One-sided exact binomial test: P(X <= k) for X ~ Binomial(n, p),
#      i.e. the chance of seeing this few intein peptides or fewer if
#      peptides were detected evenly along the protein.
#
# IN PLAIN WORDS
#   Trypsin always cuts after K and R, so we can list on the computer every
#   peptide each part of the protein COULD give. That tells us what share of
#   peptides "should" come from the intein. We then ask how unlikely it is to
#   see as few intein peptides as were actually seen.
#
# INPUT   the Ensembl FASTA at the relative path in `fa` below (PXD009315
#         deposit, sequence/ folder). Run this script from a folder where
#         that path exists, or edit `fa`.
# OUTPUT  printed to the screen: predicted counts, expected share and P value
#         for each protein at each missed-cleavage setting.
# RUN     python3 depl.py
#
# IMPORTANT SETTINGS / CAVEATS
#   * The PROLINE RULE IS NOT APPLIED here (a K or R followed by P is still
#     treated as a cut site). check_known_inteins.py DOES apply it. On the
#     UniProt reference this changes the expected share by about 2
#     percentage points (Methods 3.6). State the setting with every P value.
#   * The observed counts (k, n) are typed in by hand in T. They are distinct
#     peptides, not PSMs, because repeated spectra of one peptide are not
#     independent observations.
#   * Positions are 1-based and inclusive, matching UniProt/Ensembl numbering.
#   * DnaB is not tested here: in PXD009315 it had only 2 N-extein peptides,
#     too few for the test to have any power (Table 4.4 note). To test it,
#     add a line to T:  "CCP42780": ("DnaB", 400, 815, k, n)
#
# NEEDS   Python 3.8+ standard library only (math.comb).
# =============================================================================
import re
from math import comb
# Ensembl H37Rv proteome + contaminants, exactly as searched by the submitters
fa = "../authors_results/20170518-Ensembl-Mycobacterium_tuberculosis_h37rv-Cntms.ASM19595v2.pep.all.fasta"
# Ensembl ID: (name, intein first residue, intein last residue,
#              observed intein peptides k, observed distinct peptides n)
T = {"CCP45535": ("RecA", 252, 691, 0, 6), "CCP44220": ("SufB", 253, 611, 6, 23)}
# --- read the FASTA into {first header word: sequence} -----------------------
seq, name = {}, None
for l in open(fa):
    l = l.strip()
    if l[:1] == ">": name = l[1:].split()[0]; seq[name] = ""
    elif name: seq[name] += l
for acc, (prot, a, b, k, n) in T.items():
    s = seq[acc]
    # cut positions as 0-based slice boundaries: 0, after every K/R (except the
    # last residue), and the protein end. NO proline rule (see caveats above).
    cuts = [0] + [m.end() for m in re.finditer("[KR]", s) if m.end() < len(s)] + [len(s)]
    print("%s: observed %d of %d distinct peptides in the intein (%.0f%%); intein is %.0f%% of the length"
          % (prot, k, n, 100.0 * k / n, 100.0 * (b - a + 1) / len(s)))
    for mc in (0, 1, 2):
        # c counts predicted peptides by region: N-extein, Intein, C-extein,
        # and X = spanning a boundary (excluded from the expectation)
        c = {"N": 0, "I": 0, "C": 0, "X": 0}
        for i in range(len(cuts) - 1):
            # a peptide with mc missed cleavages spans mc + 1 cut intervals,
            # so j runs from i+1 up to i+mc+1
            for j in range(i + 1, min(i + mc + 2, len(cuts))):
                st, en = cuts[i] + 1, cuts[j]            # 1-based, inclusive
                if 6 <= en - st + 1 <= 40:               # length filter
                    c["N" if en < a else "C" if st > b else "I" if st >= a and en <= b else "X"] += 1
        # expected share of (non-spanning) peptides that fall in the intein
        p = c["I"] / (c["N"] + c["I"] + c["C"])
        # one-sided exact binomial tail: P(X <= k), X ~ Binomial(n, p)
        pv = sum(comb(n, x) * p**x * (1 - p)**(n - x) for x in range(k + 1))
        print("  missed cleavages <= %d: predicted N %d, intein %d, C %d, spanning %d; expected intein share %.3f; p = %.4g"
              % (mc, c["N"], c["I"], c["C"], c["X"], p, pv))
