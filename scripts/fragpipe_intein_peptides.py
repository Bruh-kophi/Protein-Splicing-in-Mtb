#!/usr/bin/env python3
# =============================================================================
# fragpipe_intein_peptides.py
#
# Intein analysis of FragPipe results (Dataset 3: H37Rv, rifampicin 0.5x MIC
# vs DMSO, timsTOF Pro 2 diaPASEF). Written 2026-10-15.
#
# THESIS: Methods 3.6; Results 4.4.
#
# It applies to FragPipe output the same logic that mzid_intein_peptides.py,
# mzid_nontryptic.py and jh.py applied to MS-GF+ output for PXD009315, so
# the method is identical across datasets; only the file reader differs.
#
# -----------------------------------------------------------------------------
# WHAT IT REPORTS, FOR RecA, SufB AND DnaB (all three, always)
# -----------------------------------------------------------------------------
#   1. Where every confidently identified peptide lies: N-extein, intein,
#      C-extein, or running across an intein boundary (unspliced precursor).
#   2. The four diagnostic peptide types (Methods Table 3.2):
#        type 1  boundary-spanning   extein continuous with intein
#        type 2  junction-spanning   crosses the new extein-extein bond
#        type 3  boundary-start      begins at the first C-extein residue
#        type 4  intein-end          ends at the intein's last residue (N)
#      plus two N-terminal cleavage flags (the step that oxidative stress is
#      reported to favour, Topilina et al. 2015):
#        4N-a    intein-start        begins at the intein's first residue
#        4N-b    N-extein-end        ends at the last N-extein residue
#      Each terminus is labelled tryptic or not. A flagged peptide whose
#      flagged end is TRYPTIC is not diagnostic (e.g. RecA: K251 precedes the
#      intein, so a peptide starting at C252 can be an ordinary tryptic one).
#   3. The background rate of non-tryptic termini across ALL confident PSMs,
#      against which a non-tryptic intein-end peptide is judged.
#   4. Per run and per group (Rif / DMSO / Pool), so treatments can be
#      compared within this dataset.
#   5. Pooled distinct-peptide counts (k intein of n) for the depletion test
#      (depletion_test.py), excluding any runs given with --exclude-run.
#   6. Optional (--pepxml): junction-crossing matches at ANY probability from
#      FragPipe's interact-*.pep.xml files, because psm.tsv holds only 1% FDR
#      results and weak matches must be visible to be judged (Methods 3.8).
#
# -----------------------------------------------------------------------------
# HOW PEPTIDES ARE PLACED (independent of FragPipe's protein assignment)
# -----------------------------------------------------------------------------
# FragPipe's "Protein Start/End" refer to one chosen protein only, and a
# peptide in an extein maps to BOTH the precursor and the SPLICED_ entry.
# So every distinct peptide is placed here by direct string search on the
# sequences in the FASTA, after I/L normalisation (I and L have the same
# mass and cannot be told apart):
#   * found in the precursor  -> region from its precursor coordinates
#   * found ONLY in the spliced form, crossing the junction -> type 2
# Coordinates reported are 1-based, inclusive, in PRECURSOR numbering, except
# junction peptides, which have no precursor coordinates and are given as
# "<N-side residues> + <C-side residues>".
#
# -----------------------------------------------------------------------------
# HOW TO RUN (on the cluster, in $PROJ/fragpipe)
# -----------------------------------------------------------------------------
#   python3 fragpipe_intein_peptides.py \
#       --fasta ../db/h37rv_uniprot_plus_spliced.fasta \
#       --psm results/psm.tsv \
#       --out-prefix intein_analysis/rifDIA \
#       --exclude-run 10989 \
#       --pepxml "results/interact-*.pep.xml"
#
#   --fasta        UniProt H37Rv + the three SPLICED_ entries (decoys optional;
#                  the .fas file with rev_ decoys also works)
#   --psm          FragPipe psm.tsv (combined, all runs)
#   --out-prefix   output path stem; the folder is created if needed
#   --exclude-run  run ID(s) left out of the pooled depletion counts (they
#                  are still listed per run). Pool_2 (10989) failed QC.
#   --pepxml       optional glob of interact-*.pep.xml for any-probability
#                  junction matches
#
# OUTPUTS
#   <prefix>_report.txt         human-readable report (also printed)
#   <prefix>_peptides.tsv       one row per (run, peptide, protein placement)
#   <prefix>_depletion_counts.tsv   protein, k (intein), n, runs used
#   <prefix>_junction_anyprob.tsv   (with --pepxml) weak + strong junction hits
#
# NEEDS: Python 3.8+ standard library only. Reads files; changes nothing.
# A HIT IS A LEAD, NOT A FINDING: check its spectrum (Methods 3.8).
# =============================================================================

import argparse
import csv
import glob
import os
import re
import sys
import xml.etree.ElementTree as ET
from collections import defaultdict

# name, UniProt accession, intein first residue, intein last residue (1-based)
INTEINS = [("RecA", "P9WHJ3", 252, 691),
           ("SufB", "P9WFP7", 253, 611),
           ("DnaB", "P9WMR3", 400, 815)]

GROUP_PATTERN = re.compile(r"_(Rif|DMSO|Pool)_(\d+)_")


def norm(s):
    """I/L normalisation: isoleucine and leucine are isobaric."""
    return s.upper().replace("I", "L")


def find_all(hay, needle):
    """Every 0-based start of needle in hay, overlaps included."""
    out, i = [], hay.find(needle)
    while i != -1:
        out.append(i)
        i = hay.find(needle, i + 1)
    return out


def read_fasta(path):
    """Return {header first word: sequence} and {header first word: full header}."""
    seqs, heads, pid, buf = {}, {}, None, []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                if pid is not None:
                    seqs[pid] = "".join(buf).upper()
                pid, buf = line[1:].split()[0], []
                heads[pid] = line[1:]
            elif line:
                buf.append(line)
    if pid is not None:
        seqs[pid] = "".join(buf).upper()
    return seqs, heads


def load_proteins(fasta):
    """Find each host protein's precursor and SPLICED_ entry in the FASTA.

    Returns a dict name -> info, and a list of warnings. Decoys (rev_...)
    are ignored. The spliced entry is checked against precursor[:a-1] +
    precursor[b:], so a database built from a different reference is caught.
    """
    seqs, _ = read_fasta(fasta)
    info, warn = {}, []
    for name, acc, a, b in INTEINS:
        prec = [k for k in seqs if not k.startswith("rev_")
                and (k == acc or "|%s|" % acc in k)]
        spl = [k for k in seqs if k.startswith("SPLICED_%s_" % name)]
        if len(prec) != 1:
            warn.append("%s: expected 1 precursor entry (%s), found %d"
                        % (name, acc, len(prec)))
            continue
        p = seqs[prec[0]]
        expected_spliced = p[:a - 1] + p[b:]
        s = seqs[spl[0]] if spl else None
        if s is None:
            warn.append("%s: no SPLICED_%s entry; junction peptides cannot be "
                        "matched by the search" % (name, name))
        elif s != expected_spliced:
            warn.append("%s: SPLICED_ entry differs from precursor with intein "
                        "%d-%d removed - check the database" % (name, a, b))
        info[name] = dict(acc=acc, a=a, b=b, prec=p, nprec=norm(p),
                          spliced=expected_spliced, nspl=norm(expected_spliced),
                          prec_id=prec[0], spl_id=spl[0] if spl else None)
    return info, warn


def run_of(spectrum):
    """FragPipe spectrum name -> run stem, e.g.
    10984_260506_0045_Rif_1_S1-A2_diatracer.12345.12345.2 -> 10984_..._S1-A2"""
    stem = spectrum.rsplit(".", 3)[0] if spectrum.count(".") >= 3 else spectrum
    return stem.replace("_diatracer", "")


def group_of(run):
    m = GROUP_PATTERN.search(run)
    return ("%s" % m.group(1)) if m else "?"


def tryptic_start(prev):
    """N-terminus made by trypsin: residue before is K/R, or protein start."""
    return prev in ("K", "R", "-", "")


def tryptic_end(last, nxt):
    """C-terminus made by trypsin: peptide ends in K/R, or protein end."""
    return last in ("K", "R") or nxt in ("-", "")


def place(pep, P):
    """All placements of one peptide on one host protein.

    Returns a list of dicts. Proline rule not applied (as mzid_nontryptic.py).
    """
    out = []
    a, b, L = P["a"], P["b"], len(P["prec"])
    npep = norm(pep)
    n = len(npep)
    for off in find_all(P["nprec"], npep):
        s, e = off + 1, off + n                      # 1-based inclusive
        prev = P["prec"][off - 1] if off > 0 else "-"
        nxt = P["prec"][off + n] if off + n < L else "-"
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
        nt, ct = tryptic_start(prev), tryptic_end(npep[-1], nxt)
        flags = []
        if region.startswith("SPANS"):
            flags.append("type1 boundary-spanning")
        if s == b + 1:
            flags.append("type3 boundary-start")
        if e == b:
            flags.append("type4 intein-end (%s end)" % ("tryptic" if ct else "NON-tryptic"))
        if s == a:
            flags.append("4N-a intein-start (%s start)" % ("tryptic" if nt else "NON-tryptic"))
        if e == a - 1:
            flags.append("4N-b N-extein-end (%s end)" % ("tryptic" if ct else "NON-tryptic"))
        out.append(dict(start=s, end=e, region=region, prev=prev, next=nxt,
                        n_tryptic=nt, c_tryptic=ct, flags=flags, coords="%d-%d" % (s, e)))
    if not out:
        # not in the precursor: does it cross the junction of the spliced form?
        j = a - 1                                    # N-extein length
        for off in find_all(P["nspl"], npep):
            if off < j < off + n:
                out.append(dict(start=None, end=None, region="JUNCTION (spliced only)",
                                prev=P["spliced"][off - 1] if off > 0 else "-",
                                next=P["spliced"][off + n] if off + n < len(P["spliced"]) else "-",
                                n_tryptic=None, c_tryptic=None,
                                flags=["type2 junction-spanning"],
                                coords="%d N-side + %d C-side" % (j - off, off + n - j)))
    return out


def pick(cols, *names):
    """First matching column name (case-insensitive), or None."""
    low = {c.lower(): c for c in cols}
    for n in names:
        if n.lower() in low:
            return low[n.lower()]
    return None


def read_psm(path, info):
    """Read psm.tsv. Returns per-PSM records for the host proteins, and
    background termini counts over ALL PSMs."""
    hits = []
    bg = defaultdict(lambda: defaultdict(int))       # run -> class -> count
    cterm_res = defaultdict(int)
    with open(path, newline="") as fh:
        rdr = csv.DictReader(fh, delimiter="\t")
        cols = rdr.fieldnames or []
        c_spec = pick(cols, "Spectrum")
        c_pep = pick(cols, "Peptide")
        c_prev = pick(cols, "Prev AA")
        c_next = pick(cols, "Next AA")
        c_start = pick(cols, "Protein Start")
        c_prob = pick(cols, "Probability", "PeptideProphet Probability")
        c_mod = pick(cols, "Modified Peptide")
        c_z = pick(cols, "Charge")
        if not c_spec or not c_pep:
            sys.exit("STOPPED: psm.tsv lacks 'Spectrum' or 'Peptide'. Columns: %s" % cols)
        if not (c_prev and c_next):
            sys.stderr.write("WARNING: no Prev AA / Next AA columns; background "
                             "termini rate will not be computed\n")
        for row in rdr:
            pep = (row[c_pep] or "").upper()
            if not pep:
                continue
            run = run_of(row[c_spec])
            # background: every confident PSM in the file
            if c_prev and c_next:
                prev, nxt = row[c_prev] or "-", row[c_next] or "-"
                start = row.get(c_start, "") if c_start else ""
                nt = tryptic_start(prev) or (prev == "M" and start == "2")
                ct = tryptic_end(pep[-1], nxt)
                kind = "full" if nt and ct else "C-nontryptic" if nt else \
                       "N-nontryptic" if ct else "both"
                bg[run][kind] += 1
                if not ct:
                    cterm_res[pep[-1]] += 1
            for name, P in info.items():
                pl = place(pep, P)
                if pl:
                    hits.append(dict(run=run, group=group_of(run), protein=name,
                                     peptide=pep,
                                     modified=row.get(c_mod, "") if c_mod else "",
                                     charge=row.get(c_z, "") if c_z else "",
                                     prob=float(row[c_prob]) if c_prob and row.get(c_prob) else None,
                                     spectrum=row[c_spec], placements=pl))
    return hits, bg, cterm_res


def scan_pepxml(pattern, info):
    """Junction-crossing rank-1 hits in pepXML at ANY probability."""
    rows = []
    for path in sorted(glob.glob(pattern)):
        for _, el in ET.iterparse(path):
            if el.tag.rsplit("}", 1)[-1] != "spectrum_query":
                continue
            spec = el.get("spectrum", "")
            for sh in el.iter():
                if sh.tag.rsplit("}", 1)[-1] != "search_hit" or sh.get("hit_rank") != "1":
                    continue
                pep = (sh.get("peptide") or "").upper()
                decoy = (sh.get("protein") or "").startswith("rev_")
                prob, expect = None, None
                for x in sh.iter():
                    t = x.tag.rsplit("}", 1)[-1]
                    if t.endswith("prophet_result") and prob is None:
                        prob = x.get("probability")
                    if t == "search_score" and x.get("name") == "expect":
                        expect = x.get("value")
                for name, P in info.items():
                    for p in place(pep, P):
                        if p["region"].startswith("JUNCTION"):
                            rows.append(dict(file=os.path.basename(path), spectrum=spec,
                                             run=run_of(spec), protein=name, peptide=pep,
                                             coords=p["coords"], decoy_top_hit=int(decoy),
                                             probability=prob, expect=expect,
                                             charge=el.get("assumed_charge", "")))
            el.clear()
    return rows


def main():
    ap = argparse.ArgumentParser(description="Intein analysis of FragPipe psm.tsv "
                                 "for RecA, SufB and DnaB.")
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--psm", required=True)
    ap.add_argument("--out-prefix", required=True)
    ap.add_argument("--exclude-run", action="append", default=[],
                    help="run ID or name fragment to leave out of pooled "
                         "depletion counts (repeatable)")
    ap.add_argument("--pepxml", default=None,
                    help="glob of interact-*.pep.xml for any-probability junction hits")
    args = ap.parse_args()

    for p in (args.out_prefix + "_report.txt", args.out_prefix + "_peptides.tsv"):
        if os.path.exists(p):
            sys.exit("STOPPED: %s already exists; not overwriting." % p)
    d = os.path.dirname(args.out_prefix)
    if d:
        os.makedirs(d, exist_ok=True)

    info, warns = load_proteins(args.fasta)
    hits, bg, cterm_res = read_psm(args.psm, info)

    lines = []
    say = lines.append
    say("FragPipe intein analysis")
    say("FASTA: %s" % args.fasta)
    say("psm.tsv: %s" % args.psm)
    say("Excluded from pooled depletion counts: %s" % (", ".join(args.exclude_run) or "none"))
    say("All coordinates 1-based, precursor numbering. I/L merged. Proline rule not applied.")
    for w in warns:
        say("WARNING: " + w)

    # ---- background rate of non-tryptic termini ----------------------------
    say("")
    say("=" * 72)
    say("BACKGROUND: termini of all confident PSMs (judge intein-end peptides against this)")
    say("=" * 72)
    say("%-38s %7s %7s %8s %8s %6s %7s" % ("run", "PSMs", "tryptic", "N-nontr", "C-nontr", "both", "%C-non"))
    tot = defaultdict(int)
    for run in sorted(bg):
        c = bg[run]
        n = sum(c.values())
        for k, v in c.items():
            tot[k] += v
        say("%-38s %7d %7d %8d %8d %6d %6.2f%%" % (run, n, c["full"], c["N-nontryptic"],
            c["C-nontryptic"], c["both"], 100.0 * (c["C-nontryptic"] + c["both"]) / n if n else 0))
    n = sum(tot.values())
    if n:
        say("%-38s %7d %7d %8d %8d %6d %6.2f%%" % ("ALL", n, tot["full"], tot["N-nontryptic"],
            tot["C-nontryptic"], tot["both"], 100.0 * (tot["C-nontryptic"] + tot["both"]) / n))
        say("last residue of non-tryptic C-termini: " +
            "  ".join("%s:%d" % kv for kv in sorted(cterm_res.items(), key=lambda kv: -kv[1])))

    # ---- per protein --------------------------------------------------------
    tsv_rows = []
    depl_rows = []
    for name, acc, a, b in INTEINS:
        say("")
        say("=" * 72)
        if name not in info:
            say("%s (%s): not in the FASTA - cannot be assessed" % (name, acc))
            continue
        P = info[name]
        L = len(P["prec"])
        say("%s (%s)  %d aa | intein %d-%d | junction ...%s|%s... (residue %d before it: %s)"
            % (name, acc, L, a, b, P["prec"][a - 9:a - 1], P["prec"][b:b + 8], a - 1, P["prec"][a - 2]))
        say("=" * 72)
        mine = [h for h in hits if h["protein"] == name]
        if not mine:
            say("  NO confident peptide from any region of %s in any run." % name)
            say("  -> not assessable in this dataset (report as a detectability limit).")
        # distinct peptides per (run, peptide, placement)
        agg = {}
        for h in mine:
            for p in h["placements"]:
                key = (h["run"], h["peptide"], p["coords"])
                r = agg.setdefault(key, dict(run=h["run"], group=h["group"], protein=name,
                                             peptide=h["peptide"], coords=p["coords"],
                                             region=p["region"], prev=p["prev"], next=p["next"],
                                             flags="; ".join(p["flags"]), psms=0, best_prob=None))
                r["psms"] += 1
                if h["prob"] is not None and (r["best_prob"] is None or h["prob"] > r["best_prob"]):
                    r["best_prob"] = h["prob"]
        rows = sorted(agg.values(), key=lambda r: (r["run"], r["coords"]))
        tsv_rows += rows

        # per-run region counts (distinct peptides)
        regions = ["N-extein", "intein", "C-extein", "SPANS intein START (unspliced)",
                   "SPANS intein END (unspliced)", "JUNCTION (spliced only)"]
        if mine:
            say("  distinct peptides per run (N-extein / intein / C-extein / spans / junction):")
            for run in sorted({r["run"] for r in rows}):
                cnt = defaultdict(set)
                for r in rows:
                    if r["run"] == run:
                        cnt[r["region"]].add(r["peptide"])
                spans = len(cnt[regions[3]]) + len(cnt[regions[4]])
                say("    %-38s %-5s %2d / %2d / %2d / %2d / %2d" % (run, group_of(run),
                    len(cnt["N-extein"]), len(cnt["intein"]), len(cnt["C-extein"]), spans,
                    len(cnt["JUNCTION (spliced only)"])))

        # diagnostic flags
        say("  diagnostic peptides:")
        flagged = [r for r in rows if r["flags"]]
        if not flagged:
            say("    none of types 1-4 (or 4N) observed")
        for r in flagged:
            say("    %-5s %-38s %-28s %s.%s.%s  PSMs %d  best prob %s  [%s]" % (
                r["group"], r["run"], r["coords"], r["prev"], r["peptide"], r["next"],
                r["psms"], "%.3f" % r["best_prob"] if r["best_prob"] is not None else "-",
                r["flags"]))

        # pooled depletion counts (distinct peptides, excluded runs dropped)
        keep = [r for r in rows if not any(x in r["run"] for x in args.exclude_run)]
        pooled = defaultdict(set)
        for r in keep:
            pooled[r["region"]].add((r["peptide"], r["coords"]))
        k = len(pooled["intein"])
        nn = k + len(pooled["N-extein"]) + len(pooled["C-extein"])
        say("  pooled distinct peptides (excluded runs dropped): N %d, intein %d, C %d, "
            "spanning %d, junction %d  ->  depletion input k = %d of n = %d"
            % (len(pooled["N-extein"]), k, len(pooled["C-extein"]),
               len(pooled[regions[3]]) + len(pooled[regions[4]]),
               len(pooled["JUNCTION (spliced only)"]), k, nn))
        depl_rows.append(dict(protein=name, accession=acc, intein_first=a, intein_last=b,
                              k_intein=k, n_total=nn,
                              runs_used=len({r["run"] for r in keep})))
        # same by group
        for g in ("Rif", "DMSO", "Pool"):
            gp = defaultdict(set)
            for r in keep:
                if r["group"] == g:
                    gp[r["region"]].add((r["peptide"], r["coords"]))
            gk = len(gp["intein"])
            gn = gk + len(gp["N-extein"]) + len(gp["C-extein"])
            say("    %-5s N %d, intein %d, C %d  (k %d of n %d)" % (g, len(gp["N-extein"]), gk,
                len(gp["C-extein"]), gk, gn))

    # ---- optional: junction matches at any probability ----------------------
    jrows = []
    if args.pepxml:
        jrows = scan_pepxml(args.pepxml, info)
        say("")
        say("=" * 72)
        say("JUNCTION-CROSSING rank-1 matches at ANY probability (pepXML: %s)" % args.pepxml)
        say("=" * 72)
        if not jrows:
            say("  none, for any of the three proteins")
        for r in sorted(jrows, key=lambda r: -(float(r["probability"]) if r["probability"] else 0)):
            say("  %-5s %-28s %-24s prob %-6s expect %-9s z %s  %s%s" % (
                r["protein"], r["peptide"], r["coords"], r["probability"], r["expect"],
                r["charge"], r["spectrum"], "  (decoy listed as top protein)" if r["decoy_top_hit"] else ""))

    # ---- write ---------------------------------------------------------------
    text = "\n".join(lines) + "\n"
    print(text)
    with open(args.out_prefix + "_report.txt", "w") as fh:
        fh.write(text)
    cols = ["run", "group", "protein", "peptide", "coords", "region", "prev", "next",
            "flags", "psms", "best_prob"]
    with open(args.out_prefix + "_peptides.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
        w.writeheader()
        w.writerows(tsv_rows)
    with open(args.out_prefix + "_depletion_counts.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(depl_rows[0]) if depl_rows else ["protein"],
                           delimiter="\t")
        w.writeheader()
        w.writerows(depl_rows)
    if args.pepxml:
        with open(args.out_prefix + "_junction_anyprob.tsv", "w", newline="") as fh:
            jc = ["file", "spectrum", "run", "protein", "peptide", "coords", "decoy_top_hit",
                  "probability", "expect", "charge"]
            w = csv.DictWriter(fh, fieldnames=jc, delimiter="\t")
            w.writeheader()
            w.writerows(jrows)
    sys.stderr.write("written: %s_report.txt, _peptides.tsv, _depletion_counts.tsv%s\n"
                     % (args.out_prefix, ", _junction_anyprob.tsv" if args.pepxml else ""))
    return 0


if __name__ == "__main__":
    sys.exit(main())
