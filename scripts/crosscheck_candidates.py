#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# REPOSITORY HEADER (added 2026-10-15; the code below is unchanged)
# Stage:    3 - joins intein_search.py output to the known intein boundaries
# Thesis:   Methods 3.5; Results 4.3 (no candidate at a known boundary)
# CAUTION:  the closing message ("two independent methods agreeing ... a
#           stronger negative") should NOT be quoted when both methods were
#           limited by the same lack of coverage, as in PXD009315 S7/S10.
# -----------------------------------------------------------------------------
# =============================================================================
# crosscheck_candidates.py
#
# Reads intein_search.py's candidates.tsv and asks one question:
#
#     did the hypothesis-free discovery search independently flag any of the
#     three KNOWN M. tuberculosis inteins, at the published coordinates?
#
# This is the join between the two scripts. Running it by hand is error-prone
# because the two tools report coordinates on different bases:
#
#   intein_search.py  writes gap_start / gap_end as PYTHON SLICE indices
#                     (0-based), since gap_seq = protein.seq[gap_start:gap_end]
#   check_known_inteins.py  reports 1-BASED residue numbers
#
# RecA's intein is residues 252-691 (1-based) and appears in candidates.tsv as
# gap_start=251, gap_end=691. Compared naively the start looks off by one.
# This script does the conversion so you don't have to think about it.
#
# -----------------------------------------------------------------------------
# HOW TO RUN
# -----------------------------------------------------------------------------
#   python3 crosscheck_candidates.py candidates.tsv
#   python3 crosscheck_candidates.py candidates.tsv --tolerance 5
#   python3 crosscheck_candidates.py candidates.tsv --max-q 0.01
#
# Needs only the standard library.
# =============================================================================

import argparse
import csv
import sys

# published coordinates, 1-BASED inclusive: name, accession, first, last
INTEINS = [
    ("RecA", "P9WHJ3", 252, 691),
    ("SufB", "P9WFP7", 253, 611),
    ("DnaB", "P9WMR3", 400, 815),
]


def to_one_based(gap_start, gap_end):
    """Convert intein_search's 0-based slice indices to 1-based inclusive."""
    return gap_start + 1, gap_end


def main():
    ap = argparse.ArgumentParser(
        description="Cross-check intein_search candidates against the three "
                    "known M. tuberculosis inteins.")
    ap.add_argument("candidates", help="candidates.tsv from intein_search.py")
    ap.add_argument("--tolerance", type=int, default=3,
                    help="residues of slack allowed at each boundary "
                         "(default: 3; greedy anchor extension can shift a "
                         "split point by a residue or two)")
    ap.add_argument("--max-q", type=float, default=None,
                    help="only consider rows at or below this q-value")
    args = ap.parse_args()

    with open(args.candidates, newline="") as fh:
        rdr = csv.DictReader(fh, delimiter="\t")
        needed = {"protein_id", "gap_start", "gap_end", "target_decoy"}
        missing = needed - set(rdr.fieldnames or [])
        if missing:
            sys.exit("candidates file is missing column(s): %s"
                     % ", ".join(sorted(missing)))
        rows = list(rdr)

    # Decoys are reversed-sequence proteins, prefixed DECOY_ by
    # intein_search.make_decoys. A reversed protein has the same length and
    # amino-acid composition as the real one but no real biology, so every
    # candidate found on a decoy is false by construction. The ratio of
    # targets to decoys is therefore a direct readout of how well the score
    # is separating signal from chance.
    targets = [r for r in rows if r.get("target_decoy") == "target"]
    decoys = [r for r in rows if r.get("target_decoy") == "decoy"]

    print("candidates file : %s" % args.candidates)
    print("rows            : %d  (%d target, %d decoy)"
          % (len(rows), len(targets), len(decoys)))
    if targets:
        ratio = len(targets) / float(len(decoys)) if decoys else float("inf")
        print("target:decoy    : %.2f : 1" % ratio)
        if decoys and ratio < 2.0:
            print("                  NOTE: a ratio below 2:1 means the score is")
            print("                  barely separating real from random. Roughly")
            print("                  %d%% of these targets are expected false."
                  % int(100 * len(decoys) / float(len(targets))))
    print()

    tol = args.tolerance
    any_hit = False

    for name, acc, a, b in INTEINS:
        # Match on the accession appearing anywhere in protein_id, because
        # intein_search.read_fasta takes the first whitespace-delimited token
        # of the header - which may be a bare accession, or the whole
        # sp|P9WHJ3|RECA_MYCTU string, depending on the FASTA. A substring
        # test handles both. Decoys are excluded: DECOY_P9WHJ3 also contains
        # the accession, and a hit there would be meaningless.
        mine = [r for r in rows
                if acc in r["protein_id"] and r.get("target_decoy") == "target"]

        print("=" * 68)
        print("%s  (%s)  published intein %d-%d (1-based)" % (name, acc, a, b))
        print("=" * 68)

        if not mine:
            print("  no candidate rows for this protein at all")
            print()
            continue

        print("  %d candidate row(s) on this protein" % len(mine))

        hits = []
        for r in mine:
            try:
                gs, ge = int(r["gap_start"]), int(r["gap_end"])
            except (ValueError, KeyError):
                continue
            s1, e1 = to_one_based(gs, ge)
            if abs(s1 - a) <= tol and abs(e1 - b) <= tol:
                if args.max_q is not None:
                    try:
                        if float(r.get("qvalue") or 1.0) > args.max_q:
                            continue
                    except ValueError:
                        continue
                hits.append((r, s1, e1))

        if not hits:
            print("  none match the published boundaries within +/-%d residues" % tol)
            # Showing the nearest row is more informative than a bare "no".
            # Off by one or two residues suggests the greedy anchor extension
            # overshot the true splice point, which is a known ambiguity
            # (see find_split_matches in intein_search.py) and means you
            # should widen --tolerance. Off by hundreds means the search found
            # something unrelated on the same protein.
            best = min(mine, key=lambda r: abs(int(r["gap_start"]) + 1 - a)
                       + abs(int(r["gap_end"]) - b))
            s1, e1 = to_one_based(int(best["gap_start"]), int(best["gap_end"]))
            print("  closest row: %d-%d  (off by %+d at the start, %+d at the end)"
                  % (s1, e1, s1 - a, e1 - b))
            print("               peptide=%s motif_count=%s validated=%s q=%s"
                  % (best.get("peptide"), best.get("motif_count"),
                     best.get("validated"), best.get("qvalue")))
        else:
            any_hit = True
            print("  MATCH at the published boundaries:")
            for r, s1, e1 in hits:
                print("     %d-%d  peptide=%s  motifs=%s  validated=%s  "
                      "mass_err=%s  score=%s  q=%s"
                      % (s1, e1, r.get("peptide"), r.get("motif_count"),
                         r.get("validated"), r.get("mass_error_da"),
                         r.get("score"), r.get("qvalue")))
        print()

    print("=" * 68)
    if any_hit:
        print("RESULT: the discovery search independently recovered at least one")
        print("known intein. Check motif_count==3 and validated==1 on those rows")
        print("before treating it as support.")
    else:
        print("RESULT: the discovery search flagged none of the three known")
        print("inteins at their published boundaries. Taken together with a")
        print("negative targeted check, that is two independent methods")
        print("agreeing - which is a stronger negative than either alone.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
