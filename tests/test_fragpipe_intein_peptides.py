#!/usr/bin/env python3
# =============================================================================
# test_fragpipe_intein_peptides.py -- checks fragpipe_intein_peptides.py on
# the SYNTHETIC fixture. Exit code 0 = every check passed.
#
# Run from the repository root:
#   python3 tests/test_fragpipe_intein_peptides.py
# =============================================================================
import csv
import os
import shutil
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPT = os.path.join(HERE, "..", "scripts", "fragpipe_intein_peptides.py")
tmp = tempfile.mkdtemp(prefix="SYNTHETIC_test_")
fails = 0


def check(cond, msg):
    global fails
    print(("  PASS  " if cond else "  FAIL  ") + msg)
    fails += 0 if cond else 1


try:
    subprocess.run([sys.executable, os.path.join(HERE, "make_synthetic_fixture.py"), tmp], check=True)
    r = subprocess.run([sys.executable, SCRIPT, "--fasta", os.path.join(tmp, "SYNTHETIC_db.fasta"),
                        "--psm", os.path.join(tmp, "SYNTHETIC_psm.tsv"),
                        "--out-prefix", os.path.join(tmp, "out", "SYNTHETIC"),
                        "--exclude-run", "10989",
                        "--pepxml", os.path.join(tmp, "SYNTHETIC_interact.pep.xml")],
                       capture_output=True, text=True)
    check(r.returncode == 0, "script runs without error")
    if r.returncode:
        print(r.stderr)
    report = open(os.path.join(tmp, "out", "SYNTHETIC_report.txt")).read()
    got = list(csv.DictReader(open(os.path.join(tmp, "out", "SYNTHETIC_peptides.tsv")), delimiter="\t"))
    exp = [e for e in csv.DictReader(open(os.path.join(tmp, "SYNTHETIC_expected.tsv")), delimiter="\t")]

    for e in exp:
        if e["run"] == "PEPXML":
            continue
        m = [g for g in got if g["run"] == e["run"] and g["peptide"] == e["peptide"]
             and g["protein"] == e["protein"]]
        ok = bool(m) and m[0]["region"] == e["region"] and (e["flag"] in m[0]["flags"])
        check(ok, "%s %s -> %s %s" % (e["protein"], e["peptide"], e["region"], e["flag"] or ""))

    check("type4 intein-end (NON-tryptic end)" in report or "type4 intein-end (tryptic end)" in report,
          "intein-end peptide reported with its terminus status")
    check("4N-a intein-start (NON-tryptic start)" in report,
          "SufB intein-start after glycine flagged NON-tryptic")
    check("NO confident peptide from any region of DnaB" in report,
          "DnaB with no psm.tsv peptides is reported as not assessable, not skipped")

    dep = {d["protein"]: d for d in csv.DictReader(open(os.path.join(tmp, "out",
           "SYNTHETIC_depletion_counts.tsv")), delimiter="\t")}
    # RecA: N 1 distinct (seen twice), intein 2, C 2 -> k 2 of n 5
    check(dep["RecA"]["k_intein"] == "2" and dep["RecA"]["n_total"] == "5",
          "RecA pooled depletion counts k=2 of n=5 (got k=%s n=%s)"
          % (dep["RecA"]["k_intein"], dep["RecA"]["n_total"]))
    # SufB: the Pool_2 intein peptide must be excluded -> k 1 of n 1
    check(dep["SufB"]["k_intein"] == "1" and dep["SufB"]["n_total"] == "1",
          "SufB counts exclude Pool_2 (k=1 of n=1; got k=%s n=%s)"
          % (dep["SufB"]["k_intein"], dep["SufB"]["n_total"]))
    check(dep["DnaB"]["n_total"] == "0", "DnaB depletion row present with n=0")

    # background: 18 PSMs on host proteins + 10 control; 2 control C-non-tryptic
    check("ALL" in report, "background table produced")

    jany = list(csv.DictReader(open(os.path.join(tmp, "out", "SYNTHETIC_junction_anyprob.tsv")),
                               delimiter="\t"))
    weak = [e for e in exp if e["run"] == "PEPXML"][0]
    check(any(j["peptide"] == weak["peptide"] and j["protein"] == "DnaB" and j["probability"] == "0.12"
              for j in jany), "weak DnaB junction hit (prob 0.12) found in pepXML, not hidden")

    again = subprocess.run([sys.executable, SCRIPT, "--fasta", os.path.join(tmp, "SYNTHETIC_db.fasta"),
                            "--psm", os.path.join(tmp, "SYNTHETIC_psm.tsv"),
                            "--out-prefix", os.path.join(tmp, "out", "SYNTHETIC")],
                           capture_output=True, text=True)
    check(again.returncode != 0 and "not overwriting" in again.stderr + again.stdout,
          "refuses to overwrite existing results")
finally:
    shutil.rmtree(tmp, ignore_errors=True)

print("\n%s: %d failure(s)" % ("ALL PASSED" if not fails else "FAILED", fails))
sys.exit(1 if fails else 0)
