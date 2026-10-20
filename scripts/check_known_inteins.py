#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# REPOSITORY HEADER (added 2026-10-15; the code below is unchanged)
# Stage:    2 - targeted check of the three known inteins
# Thesis:   Methods 3.3 (diagnostic panel, 30 peptides) and 3.4; Results 4.1,
#           4.2 (Dataset 1) and 4.3 (Dataset 2 de novo)
# Datasets: any peptide CSV with a sequence and a score column
# Note:     its in silico digest APPLIES the proline rule (cut_sites), unlike
#           depl.py / ctrl.py. Section 5 was corrected earlier: for SufB and
#           DnaB (junction after glycine) no tryptic boundary-start peptide
#           exists, and the script now says so instead of predicting one.
# -----------------------------------------------------------------------------
# =============================================================================
# check_known_inteins.py
#
# Tests a de novo peptide list against the three experimentally characterised
# Mycobacterium tuberculosis inteins, using published boundary coordinates.
#
# This is the HYPOTHESIS-DRIVEN check: it asks "did we see the specific
# peptides that splicing of RecA, SufB or DnaB would produce?" It does not
# search for novel inteins - intein_search.py does that.
#
# For each of the three host proteins it reports:
#   1. whether the boundary motifs hold in the reference sequence
#   2. whether trypsin cuts at the ligation junction (this differs per protein
#      and decides whether a junction peptide can exist at all)
#   3. which observed peptides land in the N-extein, intein and C-extein
#   4. whether any peptide runs contiguously across a boundary (evidence the
#      PRECURSOR was present, i.e. unspliced)
#   5. the boundary-start peptide - the one beginning at the C-extein's first
#      residue, which is tryptic in the spliced form but not in the precursor
#   6. the junction-spanning peptides - direct evidence of splicing
#
# -----------------------------------------------------------------------------
# HOW TO RUN
# -----------------------------------------------------------------------------
#   python3 check_known_inteins.py UP000001584_83332.fasta peptides.csv
#   python3 check_known_inteins.py ref.fasta peptides.csv --min-score 0.70
#   python3 check_known_inteins.py ref.fasta peptides.csv --out report.txt
#
# FASTA    : UniProt reference proteome. Holds the PRECURSOR (unspliced) form
#            of every protein, because DNA encodes the precursor and splicing
#            happens after translation.
# PEPTIDES : de novo peptide table (e.g. produced by mztab_to_csv.py from a
#            Casanovo mzTab). Needs a sequence column and a score column; the
#            script finds them by name and says which it used.
#
# Needs only the standard library.
#
# -----------------------------------------------------------------------------
# READING THE REPORT - WHICH SECTION MEANS WHAT
# -----------------------------------------------------------------------------
# 1  motifs        A sanity check on the published coordinates, not a result.
#                  3 of 3 means the residues at those positions are what an
#                  intein boundary must look like.
#
# 2  trypsin       READ THIS BEFORE 5 AND 6. It decides what is findable.
#                  RecA's N-extein ends in lysine, so trypsin cuts exactly at
#                  the new bond and destroys the junction peptide. SufB's and
#                  DnaB's end in glycine, so their junctions survive. Same
#                  biology, completely different detectability.
#
# 3  coverage      THE GATE. If "both exteins observed" is NO, a junction
#                  peptide could not have been detected for that protein even
#                  if it existed, so a null result in 5 and 6 means nothing.
#
# 4  unspliced     Peptides running contiguously across a boundary could only
#                  come from the precursor. Their presence argues against
#                  splicing; their absence is weak support for it.
#
# 5  boundary      INDIRECT evidence. A peptide that is tryptic in the spliced
#                  form but not in the precursor. Only applies where trypsin
#                  cuts at the junction, i.e. RecA. This is where CSPPFK comes
#                  from.
#
# 6  junction      DIRECT evidence, and the point of the whole project. A
#                  peptide spanning the new bond cannot exist in any database.
#                  A hit here is the result; it should then be confirmed by
#                  inspecting the spectrum by hand, not accepted on a string
#                  match alone.
#
# -----------------------------------------------------------------------------
# TWO THINGS THAT TRIP PEOPLE UP
# -----------------------------------------------------------------------------
# I/L normalisation. Isoleucine and leucine have identical mass, so no mass
# spectrometer can tell them apart and no de novo tool can either. Both sides
# are normalised (every I becomes L) before matching. The table prints the DE
# NOVO spelling of the peptide with positions taken from the FASTA, so
# ALELAVAQLEK may be listed at a position where the FASTA reads ALELAVAQIEK.
# That is not an error.
#
# Junction peptides can never come from a database search. The junction peptide
# spans a bond that exists only after splicing; its two halves sit hundreds of
# residues apart in the reference. A database search can only propose
# continuous stretches of the database, so it cannot report one even in
# principle. That is why this check needs de novo input.
# =============================================================================

import argparse
import csv
import re
import sys
from collections import defaultdict

# -----------------------------------------------------------------------------
# published intein coordinates, 1-BASED and inclusive
#   name, UniProt accession, first intein residue, last intein residue
# -----------------------------------------------------------------------------
INTEINS = [
    ("RecA", "P9WHJ3", 252, 691),
    ("SufB", "P9WFP7", 253, 611),
    ("DnaB", "P9WMR3", 400, 815),
]

# intein boundary motif expectations
FIRST_RESIDUE_OK = set("CSTA")   # intein's own first residue
LAST_TWO_OK = ("HN", "HQ")       # intein's last two residues
PLUS_ONE_OK = set("CST")         # first residue of the C-extein

MIN_LEN, MAX_LEN = 6, 40         # peptide lengths a run can realistically see
MISSED_CLEAVAGES = 2

SEQ_COL_CANDIDATES = ["sequence", "peptide", "peptide_sequence", "seq",
                      "modified_sequence", "Sequence", "Peptide"]
SCORE_COL_CANDIDATES = ["score", "search_engine_score", "confidence",
                        "aa_scores_mean", "peptide_score", "Score"]


# -----------------------------------------------------------------------------
# sequence helpers
# -----------------------------------------------------------------------------

def norm(seq):
    """I/L normalisation: isoleucine and leucine are isobaric, so collapse."""
    return seq.upper().replace("I", "L")


def strip_mods(seq):
    """Keep only amino-acid letters.

    De novo output may carry modification notation such as +57.021 or
    M+15.995. Those characters are not residues and would break positional
    matching, so they are removed.
    """
    return re.sub(r"[^A-Za-z]", "", seq).upper()


def read_fasta(path):
    """Return {accession: sequence}, parsing sp|ACC|NAME style headers."""
    seqs, acc, buf = {}, None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if acc:
                    seqs[acc] = "".join(buf)
                parts = line[1:].split("|")
                acc = parts[1] if len(parts) >= 3 else line[1:].split()[0]
                buf = []
            else:
                buf.append(line)
    if acc:
        seqs[acc] = "".join(buf)
    return seqs


def pick_column(fieldnames, candidates, what):
    """Find the first candidate column present, case-insensitively."""
    lower = {f.lower(): f for f in fieldnames}
    for c in candidates:
        if c.lower() in lower:
            return lower[c.lower()]
    sys.exit("Could not find a %s column. Looked for %s. "
             "The file has: %s" % (what, candidates, list(fieldnames)))


def read_peptides(path, min_score):
    """Return {peptide: (best_score, n_psms_above_threshold)}.

    One input row is one PSM (peptide-spectrum match). The same peptide
    usually appears in several rows because an abundant peptide is fragmented
    repeatedly as it elutes, so rows are collapsed to distinct sequences while
    keeping the best score and a count.
    """
    best, count = {}, defaultdict(int)
    total = 0

    with open(path, newline="") as fh:
        rdr = csv.DictReader(fh)
        seq_col = pick_column(rdr.fieldnames, SEQ_COL_CANDIDATES, "sequence")
        score_col = pick_column(rdr.fieldnames, SCORE_COL_CANDIDATES, "score")

        for row in rdr:
            total += 1
            pep = strip_mods(row.get(seq_col) or "")
            if not pep:
                continue
            try:
                sc = float(row.get(score_col) or "nan")
            except ValueError:
                continue
            if sc != sc or sc < min_score:        # NaN or below threshold
                continue
            count[pep] += 1
            if pep not in best or sc > best[pep]:
                best[pep] = sc

    return {p: (best[p], count[p]) for p in best}, total, seq_col, score_col


# -----------------------------------------------------------------------------
# trypsin
# -----------------------------------------------------------------------------

def cut_sites(seq):
    """0-based positions AFTER which trypsin cleaves: after K or R, not before P."""
    return [i for i in range(len(seq) - 1) if seq[i] in "KR" and seq[i + 1] != "P"]


def digest(seq, missed=MISSED_CLEAVAGES):
    """Return [(peptide, start_1based, end_1based, missed_cleavages)]."""
    bounds = [-1] + cut_sites(seq) + [len(seq) - 1]
    out = []
    for i in range(len(bounds) - 1):
        for j in range(i + 1, min(i + 2 + missed, len(bounds))):
            s, e = bounds[i] + 1, bounds[j]
            pep = seq[s:e + 1]
            if MIN_LEN <= len(pep) <= MAX_LEN:
                out.append((pep, s + 1, e + 1, j - i - 1))
    return out


# -----------------------------------------------------------------------------
# the report
# -----------------------------------------------------------------------------

def find_all(haystack, needle):
    """All 0-based start positions of needle in haystack, overlaps included."""
    out, i = [], haystack.find(needle)
    while i != -1:
        out.append(i)
        i = haystack.find(needle, i + 1)
    return out


def analyse(name, acc, a, b, fasta, peptides, emit):
    """Run all six checks for one intein host and emit the report lines."""
    if acc not in fasta:
        emit("=" * 72)
        emit("%s  (%s)" % (name, acc))
        emit("=" * 72)
        emit("  NOT FOUND in the FASTA under this accession - skipped")
        emit("")
        return None

    seq = fasta[acc]
    L = len(seq)
    n_ext, intein, c_ext = seq[:a - 1], seq[a - 1:b], seq[b:]
    spliced = n_ext + c_ext
    nseq = norm(seq)

    emit("=" * 72)
    emit("%s  (%s)" % (name, acc))
    emit("=" * 72)
    emit("  precursor %d aa | N-extein 1-%d | intein %d-%d (%d aa) | C-extein %d-%d"
         % (L, a - 1, a, b, len(intein), b + 1, L))
    emit("")

    # ---- 1. boundary motifs -------------------------------------------------
    first, last2, plus1 = seq[a - 1], seq[b - 2:b], (seq[b] if b < L else "")
    checks = [
        ("intein residue %d = %s" % (a, first), first in FIRST_RESIDUE_OK),
        ("residues %d-%d = %s" % (b - 1, b, last2), last2 in LAST_TWO_OK),
        ("+1 residue %d  = %s" % (b + 1, plus1), plus1 in PLUS_ONE_OK),
    ]
    emit("  1. boundary motifs (FASTA)")
    for text, ok in checks:
        emit("     %-26s %s" % (text, "ok" if ok else "MISMATCH"))
    npass = sum(1 for _, ok in checks if ok)
    emit("     -> %d of 3 motifs (%s)"
         % (npass, "coordinates confirmed" if npass == 3 else "check coordinates"))
    emit("")

    # ---- 2. trypsin at the junction ----------------------------------------
    prev = n_ext[-1] if n_ext else ""
    cuts_here = prev in "KR" and (c_ext[:1] != "P")
    emit("  2. trypsin at the junction (spliced form, %d aa)" % len(spliced))
    emit("     junction  ...%s | %s...  (residue %d joined to residue %d)"
         % (n_ext[-8:], c_ext[:8], a - 1, b + 1))
    if cuts_here:
        emit("     %s%d precedes the junction -> trypsin cuts AT the new bond."
             % (prev, a - 1))
        emit("     Junction-spanning peptides need a missed cleavage; expect them rarely.")
    else:
        emit("     %s%d precedes the junction -> not a cut site; junction peptides can form."
             % (prev, a - 1))
    emit("")

    # ---- 3. coverage --------------------------------------------------------
    placements = []
    for pep, (score, n) in peptides.items():
        np_ = norm(pep)
        for off in find_all(nseq, np_):
            s, e = off + 1, off + len(pep)          # 1-based inclusive
            if e < a:
                region = "N-extein"
            elif s > b:
                region = "C-extein"
            elif s >= a and e <= b:
                region = "intein"
            elif s < a <= e:
                region = "SPANS intein START (unspliced)"
            else:
                region = "SPANS intein END (unspliced)"
            placements.append((s, e, pep, score, n, region))
    placements.sort()

    tally = defaultdict(int)
    for _, _, _, _, _, r in placements:
        tally[r] += 1

    emit("  3. coverage (%d peptide placements)" % len(placements))
    for r in ("N-extein", "intein", "C-extein",
              "SPANS intein START (unspliced)", "SPANS intein END (unspliced)"):
        emit("     %-32s %d" % (r, tally[r]))
    both = tally["N-extein"] > 0 and tally["C-extein"] > 0
    emit("     -> both exteins observed: %s%s"
         % ("YES" if both else "NO",
            "" if both else "  (so a junction peptide could not have been detected)"))

    if placements:
        emit("")
        emit("     %-11s %-28s %-6s %-3s %s"
             % ("residues", "peptide", "best", "n", "region"))
        for s, e, pep, score, n, region in placements:
            emit("     %-11s %-28s %-6.2f %-3d %s"
                 % ("%d-%d" % (s, e), pep, score, n, region))
    emit("")

    # ---- 4. unspliced evidence ---------------------------------------------
    spanning = tally["SPANS intein START (unspliced)"] + tally["SPANS intein END (unspliced)"]
    emit("  4. unspliced evidence: %d peptide(s) run contiguously across a boundary"
         % spanning)
    emit("")

    # ---- 5. boundary-start peptide -----------------------------------------
    # the first tryptic peptide of the spliced form that begins at the
    # C-extein's first residue. Tryptic in the spliced protein; in the
    # precursor the preceding residue is the intein's last, which trypsin
    # does not cut, so the precursor cannot easily yield it.
    jn = len(n_ext)                                  # junction index, 0-based
    bdry = None
    for pep, s, e, mc in digest(spliced, missed=0):
        if s == jn + 1:
            bdry = pep
            break
    emit("  5. boundary-start peptide (begins at C-extein residue %d)" % (b + 1))
    emit("     residue before it in the PRECURSOR: %s%d  (trypsin will not cut here)"
         % (seq[b - 1], b))
    emit("     residue before it in the SPLICED form: %s%d  (%s)"
         % (prev, a - 1,
            "trypsin cuts here" if cuts_here else "trypsin will NOT cut here either"))
    if bdry is None:
        # No cut precedes the C-extein in the spliced form either, so no
        # tryptic peptide can begin at the junction. This is the SufB and DnaB
        # case: their N-exteins end in glycine. The read-through sequence is
        # printed for reference, but it is NOT a tryptic product of either
        # form and its absence is therefore not evidence of anything. For
        # these two proteins the diagnostic is the junction-SPANNING peptide
        # in section 6, which trypsin does produce.
        rt = re.match(r"[^KR]*[KR]?", spliced[jn:])
        readthrough = rt.group(0) if rt else ""
        emit("     -> no tryptic peptide can begin at the junction in either form.")
        if readthrough:
            emit("        read-through to the next K/R: %s" % readthrough)
            emit("        NOT a tryptic product of the precursor OR the spliced")
            emit("        form - shown for reference only. Its absence is not")
            emit("        evidence; see section 6 for this protein's diagnostic.")
        found_bdry = False
    else:
        emit("     predicted tryptic peptide in spliced form: %s" % bdry)
        hit = next((p for p in peptides if norm(p) == norm(bdry)), None)
        found_bdry = hit is not None
        if hit:
            score, n = peptides[hit]
            emit("     FOUND  %s  residues %d-%d  best %.2f  seen %dx"
                 % (hit, b + 1, b + len(hit), score, n))
        else:
            emit("     not observed")
    emit("")

    # ---- 6. junction-spanning peptides -------------------------------------
    jpeps = [(pep, mc) for pep, s, e, mc in digest(spliced)
             if s <= jn < e]
    jpeps.sort(key=lambda x: (x[1], len(x[0])))
    emit("  6. junction-spanning peptides (direct evidence of splicing)")
    emit("     predicted, by missed cleavages:")
    for pep, mc in jpeps:
        emit("       %d missed   %s" % (mc, pep))
    observed = [p for p, _ in jpeps if any(norm(q) == norm(p) for q in peptides)]
    if observed:
        for p in observed:
            emit("     OBSERVED  %s" % p)
    else:
        emit("     none observed")
    emit("")

    return {
        "name": name, "motifs": "%d/3" % npass,
        "cut": "yes" if cuts_here else "no",
        "N": tally["N-extein"], "C": tally["C-extein"],
        "intein": tally["intein"], "unspliced": spanning,
        "bdry": "FOUND" if found_bdry else "-",
        "junction": "FOUND" if observed else "-",
    }


def main():
    ap = argparse.ArgumentParser(
        description="Check a de novo peptide list against the three known "
                    "M. tuberculosis inteins.")
    ap.add_argument("fasta", help="UniProt reference proteome (precursor forms)")
    ap.add_argument("peptides", help="de novo peptide CSV")
    ap.add_argument("--min-score", type=float, default=0.0,
                    help="drop PSMs scoring below this (default: 0.0, keep all)")
    ap.add_argument("--out", help="also write the report to this file")
    args = ap.parse_args()

    fasta = read_fasta(args.fasta)
    peptides, total_psms, seq_col, score_col = read_peptides(
        args.peptides, args.min_score)

    lines = []

    def emit(s=""):
        lines.append(s)
        print(s)

    emit("FASTA:    %s  (%d proteins)" % (args.fasta, len(fasta)))
    emit("peptides: %s  (%d distinct at score >= %.2f, from %d PSMs)"
         % (args.peptides, len(peptides), args.min_score, total_psms))
    emit("          columns used: sequence=%r score=%r" % (seq_col, score_col))
    emit("All positions below are 1-BASED residue numbers.")
    emit("")

    summary = []
    for name, acc, a, b in INTEINS:
        row = analyse(name, acc, a, b, fasta, peptides, emit)
        if row:
            summary.append(row)

    emit("=" * 72)
    emit("SUMMARY")
    emit("=" * 72)
    emit("       %-7s %-9s %-5s %-5s %-6s %-9s %-10s %s"
         % ("motifs", "cut@junc", "N", "C", "intein", "unspliced",
            "bdry-start", "junction"))
    for r in summary:
        emit("%-6s %-7s %-9s %-5d %-5d %-6d %-9d %-10s %s"
             % (r["name"], r["motifs"], r["cut"], r["N"], r["C"], r["intein"],
                r["unspliced"], r["bdry"], r["junction"]))

    if args.out:
        with open(args.out, "w") as fh:
            fh.write("\n".join(lines) + "\n")
        print("\nreport written to: %s" % args.out)

    return 0


if __name__ == "__main__":
    sys.exit(main())
