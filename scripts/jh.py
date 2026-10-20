# =============================================================================
# jh.py -- junction hits: PSMs that cross a splice junction (PXD009315)
#
# THESIS: Methods 3.5 (database search with spliced sequences);
#         Results 4.3 ("Junction peptides"; two weak VVKNKCSPPFK matches).
#
# QUESTION
#   After the spliced (intein-removed) forms of RecA, SufB and DnaB were
#   added to the search database, did any spectrum match a peptide that
#   crosses the NEW extein-extein bond?
#
# HOW
#   Reads MS-GF+ result files (.mzid) from the search against
#   h37rv_plus_spliced.fasta. Spliced entries are named SPLICED_<protein>_...
#   For each target (non-decoy) match on a SPLICED_ entry, the peptide
#   crosses the junction if it starts at or before the last N-extein residue
#   and ends after it. Every rank-1 crossing PSM is listed AT ANY q-value,
#   so nothing is hidden by the 1% threshold; the q-value is shown so that
#   weak matches can be judged (and, if needed, rejected) by inspection.
#
# IN PLAIN WORDS
#   We gave the search engine the "joined" versions of the three proteins.
#   This script lists every spectrum it matched to a piece of protein that
#   only exists after splicing, however weak the match.
#
# INPUT   one or more .mzid files from the spliced-database search
# OUTPUT  printed: per file, the number of confident PSMs and of
#         junction-crossing PSMs (q <= 0.01 and at any q), then up to 25
#         crossing PSMs, best q first, as
#         q, SpecEValue, protein pre.PEPTIDE.post [residues each side], charge, scan
# RUN     python3 jh.py results/*_spliced.mzid
#
# COORDINATES
#   N = last N-extein residue, which is the junction position in the
#   SPLICED sequence (1-based): RecA 251, SufB 252, DnaB 399. A peptide
#   start..end (1-based) crosses the junction when start <= N < end.
#
# A HIT HERE IS A LEAD, NOT A FINDING. Accept a junction peptide only if its
# spectrum shows b and y ions bracketing the new bond (annot.py, Methods 3.8).
#
# NEEDS   Python 3 standard library only. For FragPipe output (Dataset 3),
#         use fragpipe_intein_peptides.py instead.
# =============================================================================
import sys, xml.etree.ElementTree as ET
# last N-extein residue = junction position in the spliced sequence (1-based)
N = {"RecA": 251, "SufB": 252, "DnaB": 399}
for path in sys.argv[1:]:
    # spl: DBSequence id -> protein name (spliced entries only)
    # pep: Peptide id -> sequence; ev: PeptideEvidence id -> description
    # hits: finished crossing PSMs; pend: crossing PSMs of the current
    # spectrum (scan number added when its SpectrumIdentificationResult ends)
    # conf: rank-1 PSMs at q <= 0.01 in the whole file
    spl, pep, ev, hits, pend, conf = {}, {}, {}, [], [], 0
    for _, el in ET.iterparse(path):
        t = el.tag.rsplit("}", 1)[-1]          # tag name without namespace
        if t == "DBSequence":
            a = el.get("accession", "")
            # e.g. SPLICED_RecA_CCP45535 -> "RecA"
            if a.startswith("SPLICED_"): spl[el.get("id")] = a.split("_")[1]
        elif t == "Peptide":
            pep[el.get("id")] = el[0].text      # first child = PeptideSequence
        elif t == "PeptideEvidence":
            d = el.get("dBSequence_ref")
            if d in spl and el.get("isDecoy") != "true":
                n, s, e = N[spl[d]], int(el.get("start")), int(el.get("end"))
                # crosses the junction: residues on both sides of position n
                if s <= n < e:
                    ev[el.get("id")] = "%s %s.%s.%s [%d N-side + %d C-side]" % (
                        spl[d], el.get("pre"), pep[el.get("peptide_ref")], el.get("post"), n - s + 1, e - n)
        elif t == "SpectrumIdentificationItem" and el.get("rank") == "1":
            # all cvParams/userParams of this PSM: name -> value
            v = {c.get("name"): c.get("value") for c in el if c.get("name")}
            q = float(v.get("MS-GF:QValue", 9))   # 9 = missing, never confident
            conf += q <= 0.01
            for c in el:
                r = c.get("peptideEvidence_ref")
                if r in ev: pend.append((q, float(v["MS-GF:SpecEValue"]), ev[r], el.get("chargeState")))
        elif t == "SpectrumIdentificationResult":
            # attach the scan number (last word of spectrumID) to this spectrum's hits
            hits += [h + (el.get("spectrumID").split()[-1],) for h in pend]; pend = []
        if t in ("DBSequence", "Peptide", "PeptideEvidence", "SpectrumIdentificationResult"): el.clear()
    print("\n%s: %d confident PSMs; junction-crossing: %d at q<=0.01, %d at any q" % (
        path, conf, sum(h[0] <= 0.01 for h in hits), len(hits)))
    # best (lowest q) first; at most 25 shown
    for q, sv, d, z, sc in sorted(hits)[:25]:
        print("  q=%-7.4g SpecE=%-8.2g %s z=%s %s" % (q, sv, d, z, sc))
