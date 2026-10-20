#!/usr/bin/env python3
# =============================================================================
# depletion_test.py -- general version of depl.py + ctrl.py. Written 2026-10-15.
#
# THESIS: Methods 3.3 (digest), 3.7 (depletion test and control);
#         Results 4.4 (Dataset 3). Same method as depl.py / ctrl.py, which are
#         kept unchanged as the record of the PXD009315 analysis.
#
# WHAT CHANGED FROM depl.py / ctrl.py (method unchanged, inputs generalised)
#   * observed counts come from a file (fragpipe_intein_peptides.py writes
#     <prefix>_depletion_counts.tsv) or the command line, not typed into code
#   * DnaB is always included; if it has too few peptides the script says the
#     test has no power, and prints the smallest P value its n could ever give
#   * the PROLINE RULE is a stated switch. Default ON (the setting fixed for
#     Dataset 3, Methods 3.3). Results with it OFF are printed too, so the
#     effect of the setting is always disclosed
#   * missed cleavages 0, <=1 and <=2 are all reported (as depl.py)
#   * the control reads observed peptide positions from FragPipe psm.tsv
#     (--psm) instead of MS-GF+ mzid files
#   * the empirical probability (r + 1) / (m + 1) is computed and printed
#
# IN PLAIN WORDS
#   Cut each protein on the computer where trypsin cuts. Count how many of
#   the possible peptides lie in the intein: that is the share we "expect".
#   Then ask how unlikely it is to see as few intein peptides as we did. To
#   check the test is fair, run it on ordinary proteins with a pretend
#   intein in the same position.
#
# RUN
#   python3 depletion_test.py --fasta h37rv_uniprot_plus_spliced.fasta \
#       --counts intein_analysis/rifDIA_depletion_counts.tsv \
#       --psm results/psm.tsv --exclude-run 10989 \
#       > intein_analysis/rifDIA_depletion.txt
#   or, for given counts only:
#   python3 depletion_test.py --fasta ref.fasta --observed RecA:2:12 --observed SufB:4:27
#
# INPUT FASTA: precursor sequences (UniProt sp|ACC|... headers, or bare ACC).
# NEEDS: Python 3.8+ standard library only.
# =============================================================================

import argparse
import csv
import sys
from collections import defaultdict
from math import comb

INTEINS = {"RecA": ("P9WHJ3", 252, 691),
           "SufB": ("P9WFP7", 253, 611),
           "DnaB": ("P9WMR3", 400, 815)}
MIN_LEN, MAX_LEN = 6, 40          # as depl.py and check_known_inteins.py


def read_fasta(path):
    """{accession: sequence}; accession = middle field of sp|ACC|NAME, else first word."""
    seqs, acc, buf = {}, None, []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                if acc:
                    seqs[acc] = "".join(buf).upper()
                head = line[1:].split()[0]
                parts = head.split("|")
                acc = parts[1] if len(parts) >= 3 else head
                buf = []
            elif line:
                buf.append(line)
    if acc:
        seqs[acc] = "".join(buf).upper()
    return seqs


def cut_positions(s, proline_rule):
    """0-based slice boundaries: 0, after each K/R (not the last residue), end.
    With the proline rule, a K/R followed by P is not a cut site."""
    cuts = [0]
    for i, ch in enumerate(s[:-1]):
        if ch in "KR" and not (proline_rule and s[i + 1] == "P"):
            cuts.append(i + 1)
    cuts.append(len(s))
    return cuts


def predicted(s, a, b, missed, proline_rule):
    """Predicted peptides by region: N-extein, intein, C-extein, spanning (X).
    Identical to the loop in depl.py when proline_rule is False."""
    cuts = cut_positions(s, proline_rule)
    c = {"N": 0, "I": 0, "C": 0, "X": 0}
    for i in range(len(cuts) - 1):
        for j in range(i + 1, min(i + missed + 2, len(cuts))):
            st, en = cuts[i] + 1, cuts[j]
            if MIN_LEN <= en - st + 1 <= MAX_LEN:
                c["N" if en < a else "C" if st > b else "I" if st >= a and en <= b else "X"] += 1
    return c


def binom_tail(k, n, p):
    """One-sided exact binomial P(X <= k), X ~ Binomial(n, p)."""
    return sum(comb(n, x) * p ** x * (1 - p) ** (n - x) for x in range(k + 1))


def observed_positions(psm_path, seqs, exclude):
    """{accession: set of (start, end)} from FragPipe psm.tsv, distinct
    positions, decoys and excluded runs dropped. Uses 'Protein ID' (or the
    middle field of 'Protein') with 'Protein Start'/'Protein End'."""
    obs = defaultdict(set)
    with open(psm_path, newline="") as fh:
        rdr = csv.DictReader(fh, delimiter="\t")
        cols = {c.lower(): c for c in rdr.fieldnames or []}
        need = ["spectrum", "protein start", "protein end"]
        if any(n not in cols for n in need) or not ("protein id" in cols or "protein" in cols):
            sys.exit("STOPPED: psm.tsv needs Spectrum, Protein ID (or Protein), "
                     "Protein Start and Protein End for the control.")
        for row in rdr:
            # match the run name only (text before the scan numbers), so a scan
            # number such as .10989. can never exclude a PSM by accident
            if any(x in row[cols["spectrum"]].rsplit(".", 3)[0] for x in exclude):
                continue
            if "protein id" in cols and row[cols["protein id"]]:
                acc = row[cols["protein id"]]
            else:
                p = row[cols["protein"]].split("|")
                acc = p[1] if len(p) >= 3 else p[0]
            if acc.startswith("rev_") or acc not in seqs:
                continue
            obs[acc].add((int(row[cols["protein start"]]), int(row[cols["protein end"]])))
    return obs


def region_test(s, a, b, positions, missed, proline_rule):
    """Depletion test on region a..b of sequence s with observed positions.
    Returns (k, n, p, P) or None. Spanning peptides excluded on both sides."""
    c = predicted(s, a, b, missed, proline_rule)
    inside, outside = c["I"], c["N"] + c["C"]
    k = sum(1 for st, en in positions if st >= a and en <= b)
    n = k + sum(1 for st, en in positions if en < a or st > b)
    if not n or not inside or not outside:
        return None
    p = inside / (inside + outside)
    return k, n, p, binom_tail(k, n, p)


def main():
    ap = argparse.ArgumentParser(description="Intein depletion test with pseudo-intein control.")
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--counts", help="<prefix>_depletion_counts.tsv from fragpipe_intein_peptides.py")
    ap.add_argument("--observed", action="append", default=[],
                    help="NAME:k:n, e.g. RecA:0:6 (repeatable; overrides --counts)")
    ap.add_argument("--no-proline-rule", action="store_true",
                    help="primary setting WITHOUT the proline rule (default: rule applied)")
    ap.add_argument("--psm", help="FragPipe psm.tsv; enables the pseudo-intein control")
    ap.add_argument("--control-missed", type=int, default=2,
                    help="missed cleavages for the control digest (default 2)")
    ap.add_argument("--exclude-run", action="append", default=[])
    ap.add_argument("--min-length", type=int, default=300)
    ap.add_argument("--min-peptides", type=int, default=6)
    args = ap.parse_args()

    seqs = read_fasta(args.fasta)
    obs_counts = {}
    if args.counts:
        for r in csv.DictReader(open(args.counts), delimiter="\t"):
            obs_counts[r["protein"]] = (int(r["k_intein"]), int(r["n_total"]))
    for o in args.observed:
        name, k, n = o.split(":")
        obs_counts[name] = (int(k), int(n))
    primary = not args.no_proline_rule

    print("Depletion test. Digest: trypsin, cut after K/R, peptides %d-%d aa." % (MIN_LEN, MAX_LEN))
    print("PRIMARY setting: proline rule %s. The other setting is shown for disclosure."
          % ("APPLIED" if primary else "NOT applied"))
    print("Observed = distinct peptides; spanning peptides excluded from both sides.")
    if args.exclude_run:
        print("Runs excluded: %s" % ", ".join(args.exclude_run))

    real_P = {}
    for name, (acc, a, b) in INTEINS.items():
        print("\n" + "=" * 72)
        if acc not in seqs:
            print("%s (%s): not in FASTA" % (name, acc))
            continue
        s = seqs[acc]
        k, n = obs_counts.get(name, (0, 0))
        print("%s (%s) %d aa, intein %d-%d (%.1f%% of length); observed %d intein of %d"
              % (name, acc, len(s), a, b, 100.0 * (b - a + 1) / len(s), k, n))
        print("=" * 72)
        if n == 0:
            print("  NOT TESTABLE: no peptides observed on %s (detectability limit)." % name)
        for rule in (primary, not primary):
            tag = "proline rule %s%s" % ("applied" if rule else "not applied",
                                         "  [PRIMARY]" if rule == primary else "")
            print("  %s" % tag)
            for mc in (0, 1, 2):
                c = predicted(s, a, b, mc, rule)
                p = c["I"] / (c["N"] + c["I"] + c["C"])
                line = "    MC<=%d  predicted N %3d / intein %3d / C %3d / span %2d  expected %.3f" % (
                    mc, c["N"], c["I"], c["C"], c["X"], p)
                if n:
                    P = binom_tail(k, n, p)
                    line += "  P = %.4g  (smallest P possible with n=%d: %.3g)" % (P, n, (1 - p) ** n)
                    if rule == primary and mc == args.control_missed:
                        real_P[name] = P
                print(line)
        if 0 < n < 6:
            print("  CAUTION: n = %d peptides gives the test little or no power." % n)

    if args.psm:
        obs = observed_positions(args.psm, seqs, args.exclude_run)
        skip = {acc for acc, _, _ in INTEINS.values()}
        print("\n" + "=" * 72)
        print("PSEUDO-INTEIN CONTROL (primary setting, MC <= %d): proteins >= %d aa with >= %d distinct"
              " peptides, RecA/SufB/DnaB excluded" % (args.control_missed, args.min_length, args.min_peptides))
        print("=" * 72)
        for name, (acc, a, b) in INTEINS.items():
            if name not in real_P or acc not in seqs:
                print("%s: no real P value (not testable), control skipped" % name)
                continue
            L = len(seqs[acc])
            res = []
            for x, pos in obs.items():
                if x in skip or len(seqs[x]) < args.min_length or len(pos) < args.min_peptides:
                    continue
                lx = len(seqs[x])
                r = region_test(seqs[x], int(a / L * lx) + 1, int(b / L * lx), pos,
                                args.control_missed, primary)
                if r:
                    res.append(r)
            m = len(res)
            if not m:
                print("%s: no eligible control proteins" % name)
                continue
            r_ext = sum(1 for r in res if r[3] <= real_P[name])
            print("%s: real P = %.4g; %d controls; P <= 0.05 in %d (%.1f%%); as extreme as %s: %d; "
                  "no peptide in pseudo-intein: %d; empirical P = (r+1)/(m+1) = %.3g"
                  % (name, real_P[name], m, sum(r[3] <= 0.05 for r in res),
                     100.0 * sum(r[3] <= 0.05 for r in res) / m, name, r_ext,
                     sum(r[0] == 0 for r in res), (r_ext + 1) / (m + 1)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
