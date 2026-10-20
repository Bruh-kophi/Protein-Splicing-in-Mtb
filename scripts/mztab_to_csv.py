#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# REPOSITORY HEADER (added 2026-10-15; the code below is unchanged)
# Stage:    1 - converts de novo output into the table the next scripts read
# Thesis:   Methods 3.5 (Dataset 2 de novo); Table 4.3 context
# Datasets: Casanovo mzTab (Dataset 2: S7, S10). Cascadia also writes mzTab.
# Plain words: turns the de novo program's output file into a simple table
#   (spectrum, peptide, measured m/z, charge, score). Nothing is filtered.
# Safety: refuses to overwrite an existing output file.
# -----------------------------------------------------------------------------
"""
mztab_to_csv.py -- convert a Casanovo mzTab file to the de novo peptide CSV
read by check_known_inteins.py and intein_search.py.

Usage:
    python3 mztab_to_csv.py S7.mztab S7_peptides.csv

Output columns: spectrum_id, peptide, mz, charge, score

  spectrum_id  <input file stem>:scan=<scan number>, so IDs stay unique when
               several samples are analysed
  peptide      the PLAIN sequence (mzTab column "sequence"), no modification
               text. intein_search.py adds carbamidomethyl-cysteine mass itself.
  mz           exp_mass_to_charge (the measured precursor m/z)
  charge       precursor charge
  score        search_engine_score[1] (Casanovo peptide score)

Every PSM row is written; nothing is filtered here. Apply a score threshold
in the downstream scripts. Columns are located by NAME from the PSH header
line, so the script stops with an error rather than read the wrong column.

Note on modifications: PSMs carrying a modification other than
carbamidomethyl-C (for example oxidised methionine) are still written, as
the plain sequence. Their measured mass will not equal the plain-sequence
mass, so intein_search.py's mass validation will reject them. The count is
printed so the size of that loss is known.
"""

import csv
import os
import re
import sys

NEEDED = ["sequence", "search_engine_score[1]", "charge",
          "exp_mass_to_charge", "spectra_ref", "modifications"]


def main():
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    src, dst = sys.argv[1], sys.argv[2]
    if os.path.abspath(src) == os.path.abspath(dst):
        sys.exit("STOPPED: input and output must be different files.")
    if os.path.exists(dst):
        sys.exit("STOPPED: %s already exists; not overwriting." % dst)

    stem = os.path.splitext(os.path.basename(src))[0]
    col = None
    rows = []
    n_psm = n_mod_other = n_carb_only = 0
    seen = set()

    with open(src) as handle:
        for line in handle:
            parts = line.rstrip("\n\r").split("\t")
            if parts[0] == "PSH":
                missing = [c for c in NEEDED if c not in parts]
                if missing:
                    sys.exit("STOPPED: mzTab header lacks column(s): %s"
                             % ", ".join(missing))
                col = {c: parts.index(c) for c in NEEDED}
            elif parts[0] == "PSM":
                if col is None:
                    sys.exit("STOPPED: PSM row found before the PSH header.")
                n_psm += 1
                ref = parts[col["spectra_ref"]]
                m = re.search(r"scan=(\d+)", ref)
                sid = "%s:scan=%s" % (stem, m.group(1)) if m else "%s:%s" % (stem, ref)
                if sid in seen:
                    sid = "%s:psm%d" % (sid, n_psm)
                seen.add(sid)

                mods = parts[col["modifications"]]
                if mods and mods != "null":
                    kinds = [x for x in mods.split(";") if x.strip()]
                    if all("Carbamidomethyl" in x for x in kinds):
                        n_carb_only += 1
                    else:
                        n_mod_other += 1

                rows.append([sid, parts[col["sequence"]],
                             parts[col["exp_mass_to_charge"]],
                             parts[col["charge"]],
                             parts[col["search_engine_score[1]"]]])

    if not rows:
        sys.exit("STOPPED: no PSM rows found in %s; nothing written." % src)

    with open(dst, "w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["spectrum_id", "peptide", "mz", "charge", "score"])
        writer.writerows(rows)

    print("%s: %d PSMs read, %d rows written to %s" % (src, n_psm, len(rows), dst))
    print("  carbamidomethyl-C only:            %d (handled downstream)" % n_carb_only)
    print("  with another modification:         %d (will fail mass validation)" % n_mod_other)
    print("  unmodified:                        %d" % (n_psm - n_carb_only - n_mod_other))


if __name__ == "__main__":
    main()
