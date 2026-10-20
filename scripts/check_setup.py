#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# REPOSITORY HEADER (added 2026-10-15; the code below is unchanged)
# Stage:    0 - run first, before any long de novo analysis
# Thesis:   Methods 3.9; used in the PXD009315 de novo trial (16 ok, 2 warnings)
# Datasets: any de novo peptide CSV (Dataset 2 trial; Dataset 3 if Cascadia runs)
# Plain words: checks that every file and module the de novo scripts need is
#   present BEFORE you start, so a run does not fail halfway through.
# -----------------------------------------------------------------------------
# =============================================================================
# check_setup.py
#
# Run this FIRST, before the long discovery search. It verifies that every
# file the two analysis scripts need is present and usable, so a run does not
# fail forty minutes in on a missing column or an unimportable module.
#
#   python3 check_setup.py
#   python3 check_setup.py --fasta ref.fasta --peptides peptides.csv
#
# Checks nothing destructively and writes no files.
# =============================================================================

import argparse
import csv
import importlib
import os
import sys

INTEIN_ACCESSIONS = {"P9WHJ3": "RecA", "P9WFP7": "SufB", "P9WMR3": "DnaB"}

# intein_search.py raises and stops if any of these is missing
SEARCH_REQUIRED = {"spectrum_id", "peptide", "mz", "charge"}

# check_known_inteins.py looks for one of each of these
KNOWN_SEQ = ["sequence", "peptide", "peptide_sequence", "seq",
             "modified_sequence"]
KNOWN_SCORE = ["score", "search_engine_score", "confidence",
               "aa_scores_mean", "peptide_score"]

ok_count = 0
fail_count = 0
warn_count = 0


def ok(msg):
    global ok_count
    ok_count += 1
    print("  [ ok ] %s" % msg)


def fail(msg):
    global fail_count
    fail_count += 1
    print("  [FAIL] %s" % msg)


def warn(msg):
    global warn_count
    warn_count += 1
    print("  [warn] %s" % msg)


def check_python():
    """3.6+ is needed for f-string-free but ordered-dict-dependent code."""
    print("Python")
    v = sys.version_info
    if v >= (3, 6):
        ok("version %d.%d.%d" % (v.major, v.minor, v.micro))
    else:
        fail("version %d.%d is too old; need 3.6+" % (v.major, v.minor))


def check_tagmatch():
    """intein_search.py does `import tagmatch` at module level.

    If that import fails the script dies immediately, so catching it here
    saves starting a long run that cannot finish. The four functions checked
    below are the ones intein_search.py actually calls; a tagmatch.py that
    imports but lacks one of them would fail partway through instead.
    """
    print("\ntagmatch.py  (required by intein_search.py only)")
    try:
        tm = importlib.import_module("tagmatch")
    except ImportError as exc:
        fail("cannot import tagmatch: %s" % exc)
        fail("intein_search.py will NOT run. Put tagmatch.py in this folder.")
        return
    ok("imports from %s" % getattr(tm, "__file__", "?"))
    for name in ("mz2mw", "peptide_mass", "TagMatch", "missed_cleavages"):
        if hasattr(tm, name):
            ok("tagmatch.%s present" % name)
        else:
            fail("tagmatch.%s MISSING - intein_search.py uses it" % name)


def check_scripts():
    """Report which scripts are present, and which intein_search version."""
    print("\nanalysis scripts")
    for name in ("check_known_inteins.py", "crosscheck_candidates.py"):
        if os.path.exists(name):
            ok("%s present" % name)
        else:
            warn("%s not in this folder" % name)

    if not os.path.exists("intein_search.py"):
        warn("intein_search.py not in this folder")
        return
    size = os.path.getsize("intein_search.py")
    text = open("intein_search.py", errors="replace").read()
    has_contig = "--keep-contiguous" in text
    has_corrob = "--corroborating-out" in text
    print("  intein_search.py: %d bytes" % size)
    # Detect the version by looking for the flag string, NOT by file size.
    # Size is unreliable: adding comments to a script changes its size without
    # changing what it does.
    if has_contig:
        ok("has the contiguity filter (--keep-contiguous) - current version")
    else:
        warn("NO contiguity filter. This is the earlier version. It will")
        warn("report candidates whose peptide occurs contiguously in the same")
        warn("protein and therefore needs no splice to explain it. On the")
        warn("trial data that was 8,932 pairs and 87 extra targets (628 vs")
        warn("541). If another copy of intein_search.py on your machine")
        warn("contains the string '--keep-contiguous', use that one instead.")
    if has_corrob:
        ok("has --corroborating-out")
    else:
        warn("no --corroborating-out flag (corroborating_evidence() is defined")
        warn("but unreachable from the command line in this version)")


def check_fasta(path):
    """Confirm the FASTA parses and contains the three intein host proteins.

    A FASTA missing one of them is not an error - the targeted check simply
    skips that protein and says so - but it silently removes a third of the
    analysis, which is worth knowing before rather than after.
    """
    print("\nFASTA: %s" % path)
    if not os.path.exists(path):
        fail("not found")
        return
    n = 0
    found = {}
    acc = None
    with open(path, errors="replace") as fh:
        for line in fh:
            if line.startswith(">"):
                n += 1
                parts = line[1:].split("|")
                acc = parts[1] if len(parts) >= 3 else line[1:].split()[0]
                if acc in INTEIN_ACCESSIONS:
                    found[acc] = INTEIN_ACCESSIONS[acc]
    if n:
        ok("%d sequences parsed" % n)
    else:
        fail("no FASTA headers found - is this really a FASTA?")
    for a, name in INTEIN_ACCESSIONS.items():
        if a in found:
            ok("%s (%s) present" % (name, a))
        else:
            warn("%s (%s) NOT in this FASTA - that protein will be skipped"
                 % (name, a))


def check_peptides(path):
    """Check the peptide file against BOTH scripts' requirements.

    They differ. intein_search.py needs spectrum_id, peptide, mz and charge,
    because tagmatch validates by mass; it raises ValueError and stops if any
    is absent. check_known_inteins.py only does string matching, so a
    sequence column and a score column are enough. A file with
    spectrum_id,peptide,mz,charge,score satisfies both.
    """
    print("\npeptides: %s" % path)
    if not os.path.exists(path):
        fail("not found")
        return
    with open(path, newline="", errors="replace") as fh:
        rdr = csv.DictReader(fh)
        cols = rdr.fieldnames or []
        rows = 0
        for _ in rdr:
            rows += 1
    if not cols:
        fail("no header row")
        return
    ok("%d columns, %d data rows" % (len(cols), rows))
    print("  columns: %s" % ", ".join(cols))

    # intein_search.py
    missing = SEARCH_REQUIRED - set(cols)
    if missing:
        fail("intein_search.py needs %s - MISSING: %s"
             % (sorted(SEARCH_REQUIRED), sorted(missing)))
    else:
        ok("intein_search.py: all four required columns present")

    # check_known_inteins.py
    lower = {c.lower(): c for c in cols}
    seq = next((lower[c] for c in KNOWN_SEQ if c in lower), None)
    sc = next((lower[c] for c in KNOWN_SCORE if c in lower), None)
    if seq and sc:
        ok("check_known_inteins.py: will use sequence=%r score=%r" % (seq, sc))
    else:
        fail("check_known_inteins.py: could not find a %s column"
             % ("sequence" if not seq else "score"))


def main():
    ap = argparse.ArgumentParser(
        description="Verify the intein analysis folder before running.")
    ap.add_argument("--fasta", default="UP000001584_83332.fasta")
    ap.add_argument("--peptides", default="peptides.csv")
    args = ap.parse_args()

    print("checking folder: %s\n" % os.getcwd())
    check_python()
    check_tagmatch()
    check_scripts()
    check_fasta(args.fasta)
    check_peptides(args.peptides)

    print("\n" + "=" * 60)
    print("%d ok, %d warnings, %d failures" % (ok_count, warn_count, fail_count))
    if fail_count:
        print("Fix the failures before running. Warnings are safe to proceed")
        print("through but change what the results mean - read them.")
        return 1
    print("Ready. Run check_known_inteins.py first, then intein_search.py.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
