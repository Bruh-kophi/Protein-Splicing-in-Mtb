#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# REPOSITORY HEADER (added 2026-10-15; the code below is unchanged)
# Stage:    2 - safety net for short junction peptides
# Thesis:   Methods 3.5; Results 4.3 (0 peptides at any junction, S7 and S10)
# Datasets: any peptide CSV with a 'peptide' column (de novo or database)
# Plain words: removes each intein on the computer and lists ANY observed
#   peptide that crosses the new join or starts right after it, with no
#   enzyme rule and no minimum length on either side (catches NKCSPPFK).
# -----------------------------------------------------------------------------
"""
junction_scan.py -- report every de novo peptide that touches a known intein
junction, with NO enzyme rule and NO minimum overhang.

Usage:
    python3 junction_scan.py UP000001584_83332.fasta S7_peptides.csv [--out report.txt]

Complements intein_search.py, which needs at least 4 residues on each side of
a junction and a peptide of at least 8 residues, so it cannot see e.g.
NKCSPPFK (2 residues on the N side).

For RecA, SufB and DnaB the intein is excised in silico. Each distinct
peptide (I/L treated as the same residue) is then placed on the spliced
sequence and reported if it is:

  SPANNING        covers residues on both sides of the junction, and does NOT
                  occur in the precursor (so it needs the splice)
  BOUNDARY-START  begins exactly at the first C-extein residue

For every reported peptide the script also prints how many proteins in the
whole FASTA contain it, because a short peptide can occur elsewhere by
chance, and the best Casanovo score among its PSMs.

A hit here is a lead, not a finding: inspect the spectrum before believing it.
Coordinates are 1-based residue numbers in the precursor.
"""

import argparse
import csv
import re
import sys

# name, accession, intein first residue, intein last residue (1-based)
INTEINS = [("RecA", "P9WHJ3", 252, 691),
           ("SufB", "P9WFP7", 253, 611),
           ("DnaB", "P9WMR3", 400, 815)]

MOD = re.compile(r"\([^)]*\)|\[[^\]]*\]|[^A-Za-z]")


def norm(s):
    return s.upper().replace("I", "L")


def read_fasta(path):
    out, pid, chunks = [], None, []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                if pid is not None:
                    out.append((pid, "".join(chunks).upper()))
                pid, chunks = line[1:].split()[0], []
            elif line:
                chunks.append(line)
    if pid is not None:
        out.append((pid, "".join(chunks).upper()))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("fasta")
    ap.add_argument("peptides")
    ap.add_argument("--out")
    args = ap.parse_args()

    proteins = read_fasta(args.fasta)
    normed = [(pid, norm(seq)) for pid, seq in proteins]

    best, npsm = {}, {}
    with open(args.peptides, newline="") as fh:
        reader = csv.DictReader(fh)
        if "peptide" not in (reader.fieldnames or []):
            sys.exit("STOPPED: peptide CSV has no 'peptide' column.")
        for row in reader:
            pep = norm(MOD.sub("", row["peptide"] or ""))
            if not pep:
                continue
            try:
                sc = float(row.get("score") or "nan")
            except ValueError:
                sc = float("nan")
            npsm[pep] = npsm.get(pep, 0) + 1
            if pep not in best or (sc == sc and not best[pep] >= sc):
                best[pep] = sc

    lines = []
    say = lines.append
    say("FASTA: %s (%d proteins)" % (args.fasta, len(proteins)))
    say("peptides: %s (%d distinct, I/L merged)" % (args.peptides, len(best)))
    say("No enzyme rule, no minimum overhang, no score threshold.")
    total = 0

    for name, acc, i_first, i_last in INTEINS:
        hits = [(pid, seq) for pid, seq in proteins if "|%s|" % acc in pid or pid == acc]
        say("")
        say("=" * 70)
        if len(hits) != 1:
            say("%s (%s): expected 1 FASTA entry, found %d -- skipped" % (name, acc, len(hits)))
            continue
        pid, seq = hits[0]
        n_len = i_first - 1                       # N-extein length
        spliced = norm(seq[:n_len] + seq[i_last:])
        precursor = norm(seq)
        say("%s (%s): junction joins residue %d%s to residue %d%s"
            % (name, acc, n_len, seq[n_len - 1], i_last + 1, seq[i_last]))
        say("  spliced context: ...%s | %s..." % (seq[n_len - 8:n_len], seq[i_last:i_last + 8]))
        found = []
        for pep in best:
            for m in re.finditer("(?=%s)" % re.escape(pep), spliced):
                s, e = m.start(), m.start() + len(pep)   # 0-based on spliced
                if s < n_len < e:
                    if pep in precursor:
                        continue
                    kind = "SPANNING"
                    detail = "%d residue(s) N side + %d C side" % (n_len - s, e - n_len)
                elif s == n_len:
                    kind = "BOUNDARY-START"
                    detail = "starts at residue %d" % (i_last + 1)
                else:
                    continue
                elsewhere = sum(1 for q, qs in normed if pep in qs and q != pid)
                found.append((kind, pep, detail, best[pep], npsm[pep], elsewhere))
        if not found:
            say("  no peptide touches this junction")
        for kind, pep, detail, sc, n, elsewhere in sorted(found):
            say("  %-14s %s" % (kind, pep))
            say("                 %s | best score %.3f | %d PSM(s) | also in %d other protein(s)"
                % (detail, sc, n, elsewhere))
        total += len(found)

    say("")
    say("TOTAL peptides touching a junction: %d" % total)
    text = "\n".join(lines)
    print(text)
    if args.out:
        with open(args.out, "w") as fh:
            fh.write(text + "\n")
        print("report written to: %s" % args.out)


if __name__ == "__main__":
    main()
