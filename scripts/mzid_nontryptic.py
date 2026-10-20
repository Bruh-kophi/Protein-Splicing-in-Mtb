#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# REPOSITORY HEADER (added 2026-10-15; the code below is unchanged)
# Stage:    4 - database-search results (MS-GF+ mzIdentML)
# Thesis:   Methods 3.5; Results 4.3 (0.35% non-tryptic C-termini; SufB
#           586-611 ending at N611, the intein's last residue)
# Datasets: Dataset 2 (PXD009315, Ensembl gene tags). For FragPipe use
#           fragpipe_intein_peptides.py.
# Plain words: trypsin cuts after K and R, so nearly every peptide ends in
#   one. A peptide ending where trypsin cannot cut means the protein was
#   already broken there - for an intein's last residue, by the intein itself.
# Note: the proline rule is not applied; only the FIRST PeptideEvidence of a
#   PSM is used (a peptide shared by two proteins is counted once).
# -----------------------------------------------------------------------------
"""mzid_nontryptic.py -- how many confident peptide identifications have a
non-tryptic end, and which ones are on RecA, SufB and DnaB?
Usage: python3 mzid_nontryptic.py FILE.mzid [FILE.mzid ...]
Counts rank-1, non-decoy PSMs with MS-GF QValue <= 0.01.
An N-terminal end is tryptic if the residue before the peptide is K or R, or
the peptide starts the protein (with or without the first Met removed).
A C-terminal end is tryptic if the peptide ends in K or R, or ends the protein."""
import sys
import xml.etree.ElementTree as ET

GENES = {"gene:Rv2737c ": "RecA", "gene:Rv1461 ": "SufB", "gene:Rv0058 ": "DnaB"}

def scan(path):
    target, pepseq, evid = {}, {}, {}
    counts = {"full": 0, "N": 0, "C": 0, "both": 0}
    cterm, rows = {}, {}
    for _, el in ET.iterparse(path, events=("end",)):
        tag = el.tag.rsplit("}", 1)[-1]
        if tag == "DBSequence":
            desc = " ".join(c.get("value", "") for c in el) + " "
            for key, name in GENES.items():
                if key in desc:
                    target[el.get("id")] = name
            el.clear()
        elif tag == "Peptide":
            for c in el:
                if c.tag.endswith("PeptideSequence"):
                    pepseq[el.get("id")] = c.text
            el.clear()
        elif tag == "PeptideEvidence":
            if el.get("isDecoy") != "true":
                evid[el.get("id")] = (el.get("dBSequence_ref"), int(el.get("start")),
                                      int(el.get("end")), el.get("pre"), el.get("post"),
                                      pepseq.get(el.get("peptide_ref"), "?"))
            el.clear()
        elif tag == "SpectrumIdentificationItem":
            if el.get("rank") == "1":
                q, ref = None, None
                for c in el:
                    if c.tag.endswith("PeptideEvidenceRef") and ref is None:
                        ref = c.get("peptideEvidence_ref")
                    elif c.get("name") == "MS-GF:QValue":
                        q = float(c.get("value"))
                if q is not None and q <= 0.01 and ref in evid:
                    db, s, e, pre, post, seq = evid[ref]
                    n_ok = pre in ("K", "R", "-") or (s == 2 and pre == "M")
                    c_ok = seq[-1] in "KR" or post == "-"
                    kind = "full" if n_ok and c_ok else "C" if n_ok else "N" if c_ok else "both"
                    counts[kind] += 1
                    if not c_ok:
                        cterm[seq[-1]] = cterm.get(seq[-1], 0) + 1
                    if db in target:
                        k = (target[db], s, e, seq, pre, post, kind)
                        rows[k] = rows.get(k, 0) + 1
        elif tag == "SpectrumIdentificationResult":
            el.clear()
    return counts, cterm, rows

def main():
    tot = {"full": 0, "N": 0, "C": 0, "both": 0}
    ct, allrows = {}, {}
    print("%-22s %7s %7s %8s %8s %6s %7s" % ("file", "PSMs", "tryptic", "N-nontr", "C-nontr", "both", "%non"))
    for path in sys.argv[1:]:
        c, cterm, rows = scan(path)
        n = sum(c.values())
        print("%-22s %7d %7d %8d %8d %6d %6.1f%%" % (path, n, c["full"], c["N"], c["C"], c["both"],
              100.0 * (n - c["full"]) / n if n else 0))
        for k in tot: tot[k] += c[k]
        for k, v in cterm.items(): ct[k] = ct.get(k, 0) + v
        for k, v in rows.items():
            r, p = allrows.get(k, (0, 0)); allrows[k] = (r + 1, p + v)
    n = sum(tot.values())
    print("%-22s %7d %7d %8d %8d %6d %6.1f%%" % ("ALL", n, tot["full"], tot["N"], tot["C"], tot["both"],
          100.0 * (n - tot["full"]) / n if n else 0))
    print("\nlast residue of peptides with a non-tryptic C-terminal end (all files):")
    print("  " + "  ".join("%s:%d" % kv for kv in sorted(ct.items(), key=lambda kv: -kv[1])))
    print("\nRecA / SufB / DnaB peptides (pre.PEPTIDE.post), non-tryptic ends flagged:")
    for (name, s, e, seq, pre, post, kind), (r, p) in sorted(allrows.items()):
        flag = {"full": "", "N": "<-- non-tryptic START", "C": "<-- non-tryptic END",
                "both": "<-- non-tryptic BOTH ends"}[kind]
        print("  %-5s %-9s %s.%s.%s  runs %d PSMs %d  %s" % (name, "%d-%d" % (s, e), pre, seq, post, r, p, flag))

if __name__ == "__main__":
    main()
