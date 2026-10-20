#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# REPOSITORY HEADER (added 2026-10-15; the code below is unchanged)
# Stage:    4 - database-search results (MS-GF+ mzIdentML)
# Thesis:   Methods 3.5; Results 4.3 (Table 4.4 observed counts, peptide map)
# Datasets: Dataset 2 (PXD009315). Proteins are found by the Ensembl gene tag
#           in the description (gene:Rv2737c / Rv1461 / Rv0058), so it needs
#           the Ensembl database. For FragPipe (Dataset 3) use
#           fragpipe_intein_peptides.py.
# Plain words: lists every confident peptide on RecA, SufB and DnaB and says
#   whether it lies in the N-extein, intein, C-extein or across a boundary.
# -----------------------------------------------------------------------------
"""mzid_intein_peptides.py -- where do the database-search peptides for RecA,
SufB and DnaB fall relative to each intein?  Reads an MS-GF+ mzIdentML file.
Usage: python3 mzid_intein_peptides.py FILE.mzid [max_qvalue, default 0.01]
Positions are 1-based residue numbers in the precursor."""
import sys
import xml.etree.ElementTree as ET

# gene tag in the Ensembl description: (name, precursor length, intein first, intein last)
TARGETS = {"gene:Rv2737c ": ("RecA", 790, 252, 691),
           "gene:Rv1461 ": ("SufB", 846, 253, 611),
           "gene:Rv0058 ": ("DnaB", 874, 400, 815)}

def region(start, end, first, last):
    if end < first: return "N-extein"
    if start > last: return "C-extein"
    if start >= first and end <= last: return "intein"
    if start < first: return "SPANS intein START (unspliced)"
    return "SPANS intein END (unspliced)"

def main():
    path = sys.argv[1]
    maxq = float(sys.argv[2]) if len(sys.argv) > 2 else 0.01
    dbseq, pepseq, evid, hits = {}, {}, {}, {}
    for _, el in ET.iterparse(path, events=("end",)):
        tag = el.tag.rsplit("}", 1)[-1]
        if tag == "DBSequence":
            desc = " ".join(c.get("value", "") for c in el) + " "
            for key, info in TARGETS.items():
                if key in desc:
                    dbseq[el.get("id")] = info + (el.get("accession"), int(el.get("length", 0)))
            el.clear()
        elif tag == "Peptide":
            for c in el:
                if c.tag.endswith("PeptideSequence"):
                    pepseq[el.get("id")] = c.text
            el.clear()
        elif tag == "PeptideEvidence":
            if el.get("dBSequence_ref") in dbseq and el.get("isDecoy") != "true":
                evid[el.get("id")] = (el.get("dBSequence_ref"), int(el.get("start")),
                                      int(el.get("end")), pepseq.get(el.get("peptide_ref"), "?"))
            el.clear()
        elif tag == "SpectrumIdentificationItem":
            if el.get("rank") == "1":
                q, refs = None, []
                for c in el:
                    if c.tag.endswith("PeptideEvidenceRef"):
                        refs.append(c.get("peptideEvidence_ref"))
                    elif c.get("name") == "MS-GF:QValue":
                        q = float(c.get("value"))
                if q is not None and q <= maxq:
                    for r in refs:
                        if r in evid:
                            n, best = hits.get(r, (0, 1.0))
                            hits[r] = (n + 1, min(best, q))
        elif tag == "SpectrumIdentificationResult":
            el.clear()

    print("file: %s   rank-1 PSMs with MS-GF QValue <= %g" % (path, maxq))
    found = {v[0] for v in dbseq.values()}
    for name, plen, first, last in sorted(TARGETS.values()):
        print("\n" + "=" * 68)
        if name not in found:
            print("%s: no entry in this file's database section" % name); continue
        dbid = [k for k, v in dbseq.items() if v[0] == name][0]
        acc, length = dbseq[dbid][4], dbseq[dbid][5]
        print("%s (%s)  length %d%s | intein %d-%d (%.0f%% of length)"
              % (name, acc, length, "" if length == plen else "  WARNING: expected %d" % plen,
                 first, last, 100.0 * (last - first + 1) / plen))
        rows = sorted((evid[r][1], evid[r][2], evid[r][3], n, q)
                      for r, (n, q) in hits.items() if evid[r][0] == dbid)
        if not rows:
            print("  no peptides at this threshold"); continue
        tally = {}
        print("  %-9s %-34s %4s  %s" % ("residues", "peptide", "PSMs", "region"))
        for s, e, seq, n, q in rows:
            reg = region(s, e, first, last)
            d, p = tally.get(reg, (0, 0)); tally[reg] = (d + 1, p + n)
            print("  %-9s %-34s %4d  %s" % ("%d-%d" % (s, e), seq, n, reg))
        print("  summary (distinct peptides / PSMs):")
        for reg in ("N-extein", "intein", "C-extein", "SPANS intein START (unspliced)",
                    "SPANS intein END (unspliced)"):
            d, p = tally.get(reg, (0, 0))
            print("    %-32s %3d / %3d" % (reg, d, p))

if __name__ == "__main__":
    main()
