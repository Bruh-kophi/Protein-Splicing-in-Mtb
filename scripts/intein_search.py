#!/usr/bin/env python
# -----------------------------------------------------------------------------
# REPOSITORY HEADER (added 2026-10-15; the code below is unchanged)
# Stage:    3 - hypothesis-free discovery of UNKNOWN splice junctions
# Thesis:   Methods 3.4 (Dataset 1) and 3.5 (Dataset 2 de novo);
#           Results 4.2 and 4.3 (SAHH candidate = oxidised Trp misread as AM)
# Datasets: de novo peptides only. Database-search peptides are, by
#           construction, continuous stretches of a database protein, so
#           they cannot produce split matches.
# Requires: tagmatch.py in the SAME folder (imported unmodified), plus
#           pandas, numpy and biopython, which tagmatch imports.
# Version:  the current copy contains '--keep-contiguous' (contiguity filter
#           re-implemented 2026-10-08; see VERSION NOTE in the docstring).
# Coordinates in the output TSV are 0-based slice indices; convert with
#           crosscheck_candidates.py.
# -----------------------------------------------------------------------------
"""
intein_search.py -- discover intein splice junctions from de novo MS/MS peptides.

WHAT THIS DETECTS
-----------------
Inteins excise themselves post-translationally and ligate the flanking exteins
with a native peptide bond.  The FASTA (database) sequence is the precursor:

    N-extein + INTEIN + C-extein

The mature protein in the sample is:

    N-extein + C-extein

A peptide spanning the ligation junction therefore maps to TWO non-contiguous
regions of the database protein.  The residues between those two regions are
the intein.

The gap is a DELETION in PROTEIN coordinate space, not a mass gap inside the
peptide.  The peptide is mass-contiguous: its precursor mass equals the sum of
its own residues.  This is the inverse of an open-modification search, where a
gap means ADDED mass.

PROVENANCE
----------
The split-anchoring in _anchor_prefix / _anchor_suffix is adapted from
lib/sequtils.py class frameshift_peptides in the proteomics-pipelines repo,
which solves the same shape of problem for frameshifts (where the two halves
land on two DIFFERENT ORFs).  Four changes were required for inteins:

  (a) both anchors are matched within a single protein sequence; the
      "not in protseq" exclusivity conditions that force the halves onto
      separate records are removed
  (b) suffix_start > prefix_end is required
  (c) str.find() is replaced with finditer, to enumerate every coordinate
      pair rather than only the first occurrence
  (d) protein[prefix_end:suffix_start] is treated as the intein candidate

tagmatch.py is imported unmodified and used only for mass and enzyme
validation on the spliced sequence.

HOW TO RUN
----------
    python3 intein_search.py --fasta ref.fasta --peptides peptides.csv \
            --out candidates.tsv --unspliced-out unspliced.tsv

    --fasta     UniProt reference proteome. Holds PRECURSOR sequences, because
                DNA encodes the precursor and splicing happens afterwards.
    --peptides  de novo CSV. MUST have columns spectrum_id, peptide, mz,
                charge; score is optional. The script raises and stops if any
                of the four is missing, because tagmatch validates by mass and
                cannot work from sequence alone.
    --out       candidates TSV, 29 columns, one row per
                (spectrum, protein, split point).

tagmatch.py must sit in the same folder or on PYTHONPATH. Run check_setup.py
first to confirm.

WHAT TO EXPECT
--------------
A long list, mostly noise. Decoys are reversed sequences; a target:decoy ratio
near 1:1 means the score is barely separating real from random. Filter on
validated == 1 and motif_count == 3 before believing anything, and treat a
q-value computed from a handful of decoys as an order of magnitude rather than
a number (see assign_qvalues).

Coordinates in the output (gap_start, gap_end) are PYTHON SLICE INDICES, i.e.
0-based, because gap_seq = protein.seq[gap_start:gap_end]. An intein at
residues 252-691 in 1-based numbering appears as gap_start=251, gap_end=691.
crosscheck_candidates.py does this conversion for you.

VERSION NOTE
------------
Contiguity filter RE-ADDED 2026-10-08 to the 37 KB copy. A peptide that occurs contiguously
somewhere in the same protein needs no splice to explain it, and should be
discarded; search() now skips such peptide-protein pairs, and --keep-contiguous
overrides that. In the lost ~48 KB revision the filter removed 8,932 such pairs
and 87 targets (628 -> 541) on the pilot. corroborating_evidence() still has
no CLI flag in this version, so it is unreachable from the command line.

THE FIVE STAGES
---------------
    1. k-mer index      prefilter: which proteins could possibly match
    2. split anchoring  find every way to split a peptide across two regions
    3. splice + validate excise the gap, check the mass with tagmatch
    4. FDR              target-decoy q-values
    5. corroboration    supporting peptide counts, as flags not filters
"""

import argparse
import csv
import re
import sys
from collections import defaultdict

import tagmatch

# ----------------------------------------------------------------------------
# Constants
# ----------------------------------------------------------------------------

STANDARD_AA = set("ACDEFGHIKLMNPQRSTVWY")

# Intein N-terminal nucleophile (residue 1 of the intein).
INTEIN_N_NUCLEOPHILE = set("CSTA")
# Intein C-terminal dipeptide.
INTEIN_C_DIPEPTIDE = ("HN", "HQ")
# C-extein +1 residue (first residue after the intein, retained in the mature
# protein).
CEXTEIN_NUCLEOPHILE = set("CST")

# Default search settings (each can be changed on the command line):
#   min_seg   minimum residues matched on EACH side of a split (anchor length)
#   min_total minimum peptide length considered
#   min_gap / max_gap   accepted intein length range, residues (100-650)
#   prec_tol  precursor mass tolerance used by tagmatch, Da
#   gap_tol   tagmatch flank mass tolerance, Da
#   specificity "semi-specific" = one peptide end may be non-tryptic
DEFAULTS = dict(
    min_seg=4,
    min_total=8,
    min_gap=100,
    max_gap=650,
    prec_tol=0.02,
    gap_tol=0.5,
    enzymes=("Trypsin",),
    specificity="semi-specific",
    max_missed_cleavages=2,
    fixed_modifications=("Carbamidomethylation of C",),
)


# ----------------------------------------------------------------------------
# I/O
#
# Protein  -- one FASTA record. Stores both the original sequence and an
#             I/L-normalised copy, because isoleucine and leucine have
#             identical mass and no mass spectrometer can tell them apart.
#             All matching happens on the normalised copy; every sequence
#             REPORTED is recovered from the original by coordinate, so the
#             output always shows the true residues.
# DenovoPSM -- one peptide-spectrum match. One row of the input CSV. The same
#             peptide normally appears in several rows, because an abundant
#             peptide is fragmented repeatedly as it elutes.
# ----------------------------------------------------------------------------

class Protein(object):
    __slots__ = ("id", "description", "seq", "norm", "is_decoy")

    def __init__(self, pid, description, seq, is_decoy=False):
        self.id = pid
        self.description = description
        self.seq = seq.upper()
        self.norm = il_normalise(self.seq)
        self.is_decoy = is_decoy

    def __repr__(self):
        return "Protein(%s, %d aa)" % (self.id, len(self.seq))


class DenovoPSM(object):
    __slots__ = ("spectrum_id", "peptide", "norm", "mz", "charge", "score", "mw")

    def __init__(self, spectrum_id, peptide, mz, charge, score=None):
        self.spectrum_id = spectrum_id
        self.peptide = peptide.upper()
        self.norm = il_normalise(self.peptide)
        self.mz = float(mz)
        self.charge = int(charge)
        self.score = score
        self.mw = float(tagmatch.mz2mw(self.mz, self.charge))

    def __repr__(self):
        return "DenovoPSM(%s, %s, %.4f)" % (self.spectrum_id, self.peptide, self.mw)


def il_normalise(seq):
    """Isoleucine and leucine are isobaric and indistinguishable by MS/MS.

    TagMatch matches by literal regex and has no I/L equivalence, so both the
    peptide and the protein are normalised I -> L before any matching.  All
    reported sequences are recovered from the ORIGINAL protein by coordinate.
    """
    return seq.replace("I", "L")


def read_fasta(path):
    """Minimal FASTA reader.  Returns a list of Protein."""
    proteins = []
    pid = None
    desc = ""
    chunks = []

    def flush():
        if pid is not None:
            proteins.append(Protein(pid, desc, "".join(chunks)))

    with open(path) as handle:
        for line in handle:
            line = line.rstrip("\n\r")
            if not line:
                continue
            if line.startswith(">"):
                flush()
                header = line[1:].strip()
                parts = header.split(None, 1)
                pid = parts[0] if parts else header
                desc = parts[1] if len(parts) > 1 else ""
                chunks = []
            else:
                chunks.append(line.strip())
    flush()
    return proteins


MOD_PATTERN = re.compile(r"\([^)]*\)|\[[^\]]*\]|[+-]?\d+\.\d+|[^A-Za-z]")


def sanitise_peptide(raw):
    """Strip de novo modification annotations and reject non-standard residues.

    De novo engines emit sequences like ``M(+15.99)DKC(+57.02)ER`` or
    ``M[Oxidation]DK``.  TagMatch interpolates the query straight into a regex
    without escaping, so an unstripped modification annotation would either
    raise or silently match the wrong thing.

    Returns the cleaned peptide, or None if it contains residues we cannot
    score (X, B, Z all return inf from tagmatch.peptide_mass).
    """
    if raw is None:
        return None
    cleaned = MOD_PATTERN.sub("", raw).upper()
    if not cleaned:
        return None
    if not set(cleaned) <= STANDARD_AA:
        return None
    return cleaned


def read_denovo_csv(path):
    """Read the canonical de novo peptide list.

    Required columns: spectrum_id, peptide, mz, charge
    Optional column:  score

    Returns (psms, skipped) where skipped is a list of (row_number, reason).
    """
    psms = []
    skipped = []
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle)
        required = {"spectrum_id", "peptide", "mz", "charge"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(
                "de novo file is missing required column(s): %s"
                % ", ".join(sorted(missing))
            )
        for n, row in enumerate(reader, start=2):
            peptide = sanitise_peptide(row.get("peptide"))
            if peptide is None:
                skipped.append((n, "non-standard residue or empty after "
                                   "modification stripping: %r" % row.get("peptide")))
                continue
            try:
                psms.append(
                    DenovoPSM(
                        row["spectrum_id"],
                        peptide,
                        row["mz"],
                        row["charge"],
                        row.get("score"),
                    )
                )
            except (ValueError, TypeError) as exc:
                skipped.append((n, "unparseable mz/charge: %s" % exc))
    return psms, skipped


def make_decoys(proteins, prefix="DECOY_"):
    """Reversed-sequence decoys for FDR estimation."""
    return [
        Protein(prefix + p.id, "reversed decoy of " + p.id, p.seq[::-1], is_decoy=True)
        for p in proteins
    ]


# ----------------------------------------------------------------------------
# Stage 1 -- k-mer index
#
# The problem: comparing every peptide against every protein is
# O(peptides x proteins). With ~12,000 peptides and ~8,000 proteins that is
# 96 million regex scans.
#
# The fix: index every short substring (k-mer) of every protein, then look up
# only those proteins containing BOTH the peptide's first k residues and its
# last k residues. A protein missing either anchor cannot produce a split
# match, so it is never examined. k is --min-seg, default 4.
# ----------------------------------------------------------------------------

class KmerIndex(object):
    """Maps an I/L-normalised k-mer to the set of protein indices containing it.

    This is the prefilter.  Without it every peptide is regex-scanned against
    every protein, which is O(peptides x proteins) and does not scale.  With
    it, only proteins that contain BOTH the peptide's prefix anchor and its
    suffix anchor are ever examined.
    """

    def __init__(self, proteins, k):
        self.k = k
        self.index = defaultdict(set)
        for i, protein in enumerate(proteins):
            seq = protein.norm
            for pos in range(len(seq) - k + 1):
                self.index[seq[pos:pos + k]].add(i)

    def candidates(self, prefix_kmer, suffix_kmer):
        """Protein indices containing both anchors."""
        return self.index.get(prefix_kmer, set()) & self.index.get(suffix_kmer, set())


# ----------------------------------------------------------------------------
# Stage 2 -- split anchoring
#
# This is the heart of the method. Given a peptide and a protein, find every
# way the peptide could be explained as two separate stretches of that protein
# with a gap between them -- the gap being the putative intein.
#
# Anchor the peptide's first few residues, extend forward as far as it keeps
# matching. Anchor the last few residues, extend backward likewise. If the two
# matched stretches together TILE the whole peptide with nothing left over,
# and the gap between them is an intein-plausible length, that is a candidate.
#
# The tiling requirement is the strongest filter in the whole script: any
# residue of the peptide explained by neither half means the candidate is
# simply wrong.
# ----------------------------------------------------------------------------

def _all_occurrences(haystack, needle):
    """Every start offset of needle in haystack, including overlaps.

    frameshift_peptides used str.find(), which returns only the first
    occurrence.  That is harmless when the two halves live on different
    records, but wrong within a single protein containing internal repeats.
    """
    if not needle:
        return []
    return [m.start() for m in re.finditer("(?=%s)" % re.escape(needle), haystack)]


def _anchor_prefix(peptide_norm, protein_norm, min_seg):
    """Anchor on peptide_norm[:min_seg], then extend forward greedily.

    Returns a list of (protein_start, prefix_len).  prefix_len is the length of
    the maximal peptide prefix matching contiguously at that locus.
    """
    anchor = peptide_norm[:min_seg]
    results = []
    for start in _all_occurrences(protein_norm, anchor):
        length = min_seg
        while (
            length < len(peptide_norm)
            and start + length < len(protein_norm)
            and protein_norm[start + length] == peptide_norm[length]
        ):
            length += 1
        results.append((start, length))
    return results


def _anchor_suffix(peptide_norm, protein_norm, min_seg):
    """Anchor on peptide_norm[-min_seg:], then extend backward greedily.

    Returns a list of (protein_end_exclusive, suffix_len).
    """
    n = len(peptide_norm)
    anchor = peptide_norm[-min_seg:]
    results = []
    for start in _all_occurrences(protein_norm, anchor):
        end = start + min_seg
        length = min_seg
        while (
            length < n
            and start - 1 >= 0
            and protein_norm[start - 1] == peptide_norm[n - length - 1]
        ):
            start -= 1
            length += 1
        results.append((end, length))
    return results


class Candidate(object):
    """One (peptide, protein, split point) intein candidate."""

    __slots__ = (
        "psm", "protein", "split", "prefix_start", "prefix_end",
        "suffix_start", "suffix_end", "gap_start", "gap_end",
        "gap_seq", "gap_mass", "motif_n_nucleophile", "motif_c_dipeptide",
        "motif_cextein_plus1", "motif_count", "mass_error", "missed_cleavages",
        "validated", "score", "qvalue",
    )

    def __init__(self, psm, protein, split, prefix_start, suffix_end):
        n = len(psm.peptide)
        self.psm = psm
        self.protein = protein
        self.split = split
        self.prefix_start = prefix_start
        self.prefix_end = prefix_start + split
        self.suffix_end = suffix_end
        self.suffix_start = suffix_end - (n - split)
        self.gap_start = self.prefix_end
        self.gap_end = self.suffix_start
        self.gap_seq = protein.seq[self.gap_start:self.gap_end]
        self.gap_mass = None
        self.motif_n_nucleophile = False
        self.motif_c_dipeptide = False
        self.motif_cextein_plus1 = False
        self.motif_count = 0
        self.mass_error = None
        self.missed_cleavages = None
        self.validated = False
        self.score = 0.0
        self.qvalue = None

    @property
    def gap_length(self):
        return self.gap_end - self.gap_start

    @property
    def spliced(self):
        """Precursor with the candidate intein excised, in ORIGINAL residues."""
        return self.protein.seq[:self.gap_start] + self.protein.seq[self.gap_end:]

    @property
    def nextein_segment(self):
        return self.protein.seq[self.prefix_start:self.prefix_end]

    @property
    def cextein_segment(self):
        return self.protein.seq[self.suffix_start:self.suffix_end]

    def key(self):
        return (self.psm.spectrum_id, self.protein.id, self.gap_start, self.gap_end)


def find_split_matches(psm, protein, min_seg, min_total, min_gap, max_gap):
    """Enumerate every valid split of psm across two regions of protein.

    A split is valid when:
      * the prefix and suffix anchors both match (Stage 2 anchoring)
      * suffix_start > prefix_end, i.e. the regions are ordered and disjoint
        in PROTEIN coordinates
      * prefix_len + suffix_len >= len(peptide), i.e. the two matched regions
        TILE the peptide with no uncovered middle.  This is the
        covered_length == peptide_length check inherited from
        frameshift_peptides and it is the single strongest filter here: any
        residue explained by neither half means the candidate is wrong.
      * the gap length falls in the intein range

    Greedy extension can overshoot the true splice point when the next peptide
    residue happens to continue matching at the prefix locus.  That makes the
    split point ambiguous rather than wrong, so every valid split point in the
    overlap window is emitted and de-duplicated by gap coordinates downstream.
    """
    n = len(psm.peptide)
    if n < min_total:
        return []

    prefixes = _anchor_prefix(psm.norm, protein.norm, min_seg)
    if not prefixes:
        return []
    suffixes = _anchor_suffix(psm.norm, protein.norm, min_seg)
    if not suffixes:
        return []

    seen = set()
    out = []
    for prefix_start, prefix_len in prefixes:
        for suffix_end, suffix_len in suffixes:
            if prefix_len + suffix_len < n:
                continue  # uncovered middle -- reject
            # Valid split points: peptide[:i] inside the prefix match and
            # peptide[i:] inside the suffix match.
            lo = max(min_seg, n - suffix_len)
            hi = min(n - min_seg, prefix_len)
            for split in range(lo, hi + 1):
                cand = Candidate(psm, protein, split, prefix_start, suffix_end)
                if cand.gap_end <= cand.gap_start:
                    continue  # not ordered / overlapping -- not a split
                if not (min_gap <= cand.gap_length <= max_gap):
                    continue
                k = cand.key()
                if k in seen:
                    continue
                seen.add(k)
                out.append(cand)
    return out


# ----------------------------------------------------------------------------
# Stage 3 -- splice and validate with tagmatch
#
# Stage 2 works on sequence alone, which is cheap and therefore permissive.
# This stage brings in mass, which is the expensive and discriminating test.
#
# The candidate intein is excised from the protein in silico, producing the
# spliced sequence. tagmatch is then asked whether the observed peptide is a
# mass-consistent enzymatic product of THAT sequence. tagmatch cannot see
# across an intein in the precursor, because it builds peptides from
# contiguous residues -- splicing first is what makes it applicable at all.
#
# Boundary motifs (first residue C/S/T/A, last two HN or HQ, next residue
# C/S/T) are scored here but never required. Non-canonical inteins exist, and
# a database error can reproduce any motif by chance.
# ----------------------------------------------------------------------------

def validate_candidate(cand, prec_tol, gap_tol, enzymes, specificity,
                       max_missed_cleavages, fixed_modifications):
    """Run tagmatch.TagMatch on the SPLICED sequence.

    TagMatch builds candidate peptides from contiguous target residues and
    checks their mass against the precursor, so it cannot see across an intein
    in the precursor.  Splicing first is what makes it applicable.

    Both query and target are I/L-normalised because TagMatch matches by
    literal regex.  peptide_mass treats I and L as isobaric, so mass validation
    is unaffected.
    """
    spliced = il_normalise(cand.spliced)

    # tagmatch.peptide_mass returns inf for X/B/Z.  Guard rather than let inf
    # propagate silently into the mass comparison.
    if not set(spliced) <= STANDARD_AA | {"J", "U", "O"}:
        bad = sorted(set(spliced) - (STANDARD_AA | {"J", "U", "O"}))
        if bad:
            return False, None, None

    tm = tagmatch.TagMatch(
        cand.psm.norm,
        [(0.0, cand.psm.mw, 0.0)],
        spliced,
        fixed_modifications=list(fixed_modifications),
        enzymes=list(enzymes),
        specificity=specificity,
        prec_tol=prec_tol,
        gap_tol=gap_tol,
        max_missed_cleavages=max_missed_cleavages,
    )

    if not tm.validated:
        return False, None, None

    validated = set(tm.validated_peptides)
    theoretical = tagmatch.peptide_mass(
        cand.psm.norm, fixed_modifications=list(fixed_modifications)
    )[0]
    mass_error = cand.psm.mw - theoretical

    mc = min(
        tagmatch.missed_cleavages(pep, list(enzymes)) for pep in validated
    ) if validated else None

    return True, mass_error, mc


def score_motifs(cand):
    """Intein boundary motifs.

    Scored, never required.  These are supporting evidence, not proof --
    non-canonical inteins exist and a database error can reproduce any motif
    by chance.
    """
    gap = cand.gap_seq
    protein = cand.protein.seq

    cand.motif_n_nucleophile = bool(gap) and gap[0] in INTEIN_N_NUCLEOPHILE
    cand.motif_c_dipeptide = len(gap) >= 2 and gap[-2:] in INTEIN_C_DIPEPTIDE
    cand.motif_cextein_plus1 = (
        cand.gap_end < len(protein)
        and protein[cand.gap_end] in CEXTEIN_NUCLEOPHILE
    )
    cand.motif_count = sum(
        [cand.motif_n_nucleophile, cand.motif_c_dipeptide, cand.motif_cextein_plus1]
    )
    return cand.motif_count


def compute_score(cand, min_seg):
    """Combine anchor lengths and motif hits into a single ranking score.

    Longer anchors on BOTH sides are exponentially less likely by chance, so
    the length term uses the shorter of the two segments.  Motifs add a fixed
    bonus each.  This is a heuristic for ranking, not a probability.
    """
    n_len = cand.prefix_end - cand.prefix_start
    c_len = cand.suffix_end - cand.suffix_start
    shorter = min(n_len, c_len)
    # e.g. anchors of 5 and 7 residues with all 3 motifs: 5 + 2*3 = 11
    cand.score = float(shorter) + 2.0 * cand.motif_count
    return cand.score


# ----------------------------------------------------------------------------
# Stage 4/5 -- FDR and corroborating evidence
#
# Stage 4: how many of these candidates are expected to be false? Each protein
# is searched alongside a reversed-sequence decoy. A reversed protein has the
# same length and composition but no real biology, so any candidate found on
# one is false by construction. The decoy count at a given score therefore
# estimates the false count among targets at that score. See assign_qvalues
# for why the estimate is coarse when decoys are few.
#
# Stage 5: two kinds of supporting evidence, reported and never used to
# filter --
#   unspliced_junction_evidence  peptides running CONTIGUOUSLY across a
#                                candidate boundary. Those could only come
#                                from the precursor, so they argue AGAINST
#                                splicing at that site.
#   corroborating_evidence       peptide counts inside the candidate intein
#                                and inside each extein. The excised intein
#                                persists in the sample, so peptides from it
#                                are expected if splicing occurred.
# ----------------------------------------------------------------------------

def assign_qvalues(candidates):
    """Concatenated target-decoy FDR.

    Sorted by descending score; at each threshold FDR = decoys / targets.
    q-value is the minimum FDR at or below that score (monotone from the
    bottom up).  With few candidates this estimate is coarse -- treat a
    q-value computed from a handful of decoys as an order of magnitude, not a
    number.
    """
    ordered = sorted(candidates, key=lambda c: -c.score)
    targets = 0
    decoys = 0
    fdrs = []
    for cand in ordered:
        if cand.protein.is_decoy:
            decoys += 1
        else:
            targets += 1
        fdrs.append(decoys / float(targets) if targets else 1.0)

    running = 1.0
    for i in range(len(ordered) - 1, -1, -1):
        running = min(running, fdrs[i])
        ordered[i].qvalue = min(running, 1.0)
    return ordered


def unspliced_junction_evidence(candidates, psms, min_overhang=3):
    """Peptides spanning an UNSPLICED precursor junction.

    If the intein were NOT excised in this sample, peptides would map
    contiguously across gap_start or gap_end in the precursor.  Their presence
    argues against splicing (or for a mixed population); their absence is
    supporting evidence for it.

    Reported as a separate section, never used to filter.
    """
    by_protein = defaultdict(list)
    for cand in candidates:
        if not cand.protein.is_decoy:
            by_protein[cand.protein.id].append(cand)

    rows = []
    seen = set()
    for pid, cands in by_protein.items():
        protein = cands[0].protein
        boundaries = set()
        for cand in cands:
            boundaries.add(("intein_start", cand.gap_start))
            boundaries.add(("intein_end", cand.gap_end))
        for psm in psms:
            for start in _all_occurrences(protein.norm, psm.norm):
                end = start + len(psm.norm)
                for label, pos in boundaries:
                    if start + min_overhang <= pos <= end - min_overhang:
                        key = (psm.spectrum_id, pid, start, pos)
                        if key in seen:
                            continue
                        seen.add(key)
                        rows.append(
                            dict(
                                spectrum_id=psm.spectrum_id,
                                peptide=psm.peptide,
                                protein_id=pid,
                                boundary=label,
                                boundary_pos=pos,
                                peptide_start=start,
                                peptide_end=end,
                            )
                        )
    return rows


def corroborating_evidence(candidates, psms):
    """Flags, not filters.

    within_intein  -- peptides lying wholly inside the candidate intein.  The
                      excised intein persists in the sample, so these are
                      expected if splicing occurred.
    within_nextein / within_cextein -- ordinary peptides from each extein,
                      confirming the protein is present at all.
    """
    counts = defaultdict(lambda: dict(within_intein=0, within_nextein=0,
                                      within_cextein=0))
    for cand in candidates:
        if cand.protein.is_decoy:
            continue
        key = (cand.protein.id, cand.gap_start, cand.gap_end)
        if key in counts:
            continue
        protein = cand.protein
        for psm in psms:
            for start in _all_occurrences(protein.norm, psm.norm):
                end = start + len(psm.norm)
                if start >= cand.gap_start and end <= cand.gap_end:
                    counts[key]["within_intein"] += 1
                elif end <= cand.gap_start:
                    counts[key]["within_nextein"] += 1
                elif start >= cand.gap_end:
                    counts[key]["within_cextein"] += 1
    return counts


# ----------------------------------------------------------------------------
# Driver
#
# search()        runs stages 1-4 over every peptide and returns candidates
#                 plus a statistics dict.
# TSV_COLUMNS /   define the 29-column output. candidate_row() is the single
# candidate_row   place where an internal Candidate becomes an output row, so
#                 that is the function to change if you need another column.
# main()          argument parsing, decoy generation, and the run summary
#                 printed to stderr (so stdout stays clean if you pipe).
# ----------------------------------------------------------------------------

def search(psms, proteins, min_seg=DEFAULTS["min_seg"],
           min_total=DEFAULTS["min_total"], min_gap=DEFAULTS["min_gap"],
           max_gap=DEFAULTS["max_gap"], prec_tol=DEFAULTS["prec_tol"],
           gap_tol=DEFAULTS["gap_tol"], enzymes=DEFAULTS["enzymes"],
           specificity=DEFAULTS["specificity"],
           max_missed_cleavages=DEFAULTS["max_missed_cleavages"],
           fixed_modifications=DEFAULTS["fixed_modifications"],
           require_validation=True, verbose=False, keep_contiguous=False):
    """Run the full pipeline.  Returns (candidates, stats)."""

    index = KmerIndex(proteins, min_seg)
    candidates = []
    stats = dict(
        psms=len(psms), proteins=len(proteins), anchored=0, tiled=0,
        gap_filtered=0, validated=0, contiguous_skipped=0,
    )

    for psm in psms:
        if len(psm.peptide) < min_total:
            continue
        prefix_kmer = psm.norm[:min_seg]
        suffix_kmer = psm.norm[-min_seg:]
        for pi in index.candidates(prefix_kmer, suffix_kmer):
            protein = proteins[pi]
            # Contiguity filter: a peptide found contiguously in this protein
            # needs no splice to explain it, so it cannot be junction evidence.
            if not keep_contiguous and psm.norm in protein.norm:
                stats["contiguous_skipped"] += 1
                continue
            splits = find_split_matches(
                psm, protein, min_seg, min_total, min_gap, max_gap
            )
            stats["gap_filtered"] += len(splits)
            for cand in splits:
                score_motifs(cand)
                ok, mass_error, mc = validate_candidate(
                    cand, prec_tol, gap_tol, enzymes, specificity,
                    max_missed_cleavages, fixed_modifications,
                )
                cand.validated = ok
                cand.mass_error = mass_error
                cand.missed_cleavages = mc
                if ok:
                    stats["validated"] += 1
                if ok or not require_validation:
                    try:
                        gm = tagmatch.peptide_mass(
                            cand.gap_seq, nterm=False, cterm=False,
                            fixed_modifications=list(fixed_modifications),
                        )[0]
                        cand.gap_mass = None if gm == float("inf") else round(gm, 4)
                    except KeyError:
                        cand.gap_mass = None
                    compute_score(cand, min_seg)
                    candidates.append(cand)

    assign_qvalues(candidates)
    return candidates, stats


TSV_COLUMNS = [
    "spectrum_id", "peptide", "peptide_length", "precursor_mz", "charge",
    "precursor_mw", "protein_id", "target_decoy", "split_point",
    "prefix_start", "prefix_end", "suffix_start", "suffix_end",
    "nextein_segment", "cextein_segment",
    "gap_start", "gap_end", "gap_length", "gap_mass", "gap_seq",
    "motif_n_nucleophile", "motif_c_dipeptide", "motif_cextein_plus1",
    "motif_count", "mass_error_da", "missed_cleavages", "validated",
    "score", "qvalue",
]


def candidate_row(cand):
    return dict(
        spectrum_id=cand.psm.spectrum_id,
        peptide=cand.psm.peptide,
        peptide_length=len(cand.psm.peptide),
        precursor_mz=cand.psm.mz,
        charge=cand.psm.charge,
        precursor_mw=round(cand.psm.mw, 4),
        protein_id=cand.protein.id,
        target_decoy="decoy" if cand.protein.is_decoy else "target",
        split_point=cand.split,
        prefix_start=cand.prefix_start,
        prefix_end=cand.prefix_end,
        suffix_start=cand.suffix_start,
        suffix_end=cand.suffix_end,
        nextein_segment=cand.nextein_segment,
        cextein_segment=cand.cextein_segment,
        gap_start=cand.gap_start,
        gap_end=cand.gap_end,
        gap_length=cand.gap_length,
        gap_mass="" if cand.gap_mass is None else cand.gap_mass,
        gap_seq=cand.gap_seq,
        motif_n_nucleophile=int(cand.motif_n_nucleophile),
        motif_c_dipeptide=int(cand.motif_c_dipeptide),
        motif_cextein_plus1=int(cand.motif_cextein_plus1),
        motif_count=cand.motif_count,
        mass_error_da="" if cand.mass_error is None else round(cand.mass_error, 5),
        missed_cleavages="" if cand.missed_cleavages is None else cand.missed_cleavages,
        validated=int(cand.validated),
        score=round(cand.score, 3),
        qvalue="" if cand.qvalue is None else round(cand.qvalue, 5),
    )


def write_tsv(path, candidates):
    ordered = sorted(candidates, key=lambda c: (-c.score, c.protein.id))
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=TSV_COLUMNS, delimiter="\t")
        writer.writeheader()
        for cand in ordered:
            writer.writerow(candidate_row(cand))
    return len(ordered)


def write_unspliced_tsv(path, rows):
    cols = ["spectrum_id", "peptide", "protein_id", "boundary", "boundary_pos",
            "peptide_start", "peptide_end"]
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=cols, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
    return len(rows)


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Discover intein splice junctions from de novo MS/MS peptides.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--fasta", required=True, help="protein FASTA (precursor sequences)")
    parser.add_argument("--peptides", required=True,
                        help="de novo peptide CSV: spectrum_id,peptide,mz,charge[,score]")
    parser.add_argument("--out", required=True, help="output TSV of candidates")
    parser.add_argument("--unspliced-out", default=None,
                        help="optional TSV of unspliced-junction evidence")
    parser.add_argument("--min-seg", type=int, default=DEFAULTS["min_seg"],
                        help="minimum anchor length on each side")
    parser.add_argument("--min-total", type=int, default=DEFAULTS["min_total"],
                        help="minimum peptide length")
    parser.add_argument("--min-gap", type=int, default=DEFAULTS["min_gap"],
                        help="minimum intein length (aa)")
    parser.add_argument("--max-gap", type=int, default=DEFAULTS["max_gap"],
                        help="maximum intein length (aa)")
    parser.add_argument("--prec-tol", type=float, default=DEFAULTS["prec_tol"],
                        help="precursor mass tolerance (Da)")
    parser.add_argument("--gap-tol", type=float, default=DEFAULTS["gap_tol"],
                        help="tagmatch flank mass tolerance (Da)")
    parser.add_argument("--enzyme", default="Trypsin",
                        choices=["Trypsin", "Trypsin, no P rule", "Whole protein"])
    parser.add_argument("--specificity", default=DEFAULTS["specificity"],
                        choices=["specific", "semi-specific", "unspecific"])
    parser.add_argument("--max-missed-cleavages", type=int,
                        default=DEFAULTS["max_missed_cleavages"])
    parser.add_argument("--no-decoys", action="store_true",
                        help="skip decoy search (disables FDR estimation)")
    parser.add_argument("--keep-unvalidated", action="store_true",
                        help="report candidates that fail tagmatch validation")
    parser.add_argument("--keep-contiguous", action="store_true",
                        help="keep candidates whose peptide also occurs "
                             "contiguously in the same protein")
    args = parser.parse_args(argv)

    proteins = read_fasta(args.fasta)
    # Every real protein gets a reversed twin (DECOY_...). Candidates found on
    # decoys are false by construction; their count estimates the false
    # discovery rate among targets (assign_qvalues).
    if not args.no_decoys:
        proteins = proteins + make_decoys(proteins)

    psms, skipped = read_denovo_csv(args.peptides)

    sys.stderr.write("proteins (incl. decoys): %d\n" % len(proteins))
    sys.stderr.write("de novo PSMs:            %d\n" % len(psms))
    if skipped:
        sys.stderr.write("skipped rows:            %d\n" % len(skipped))
        for n, reason in skipped[:5]:
            sys.stderr.write("  line %d: %s\n" % (n, reason))
        if len(skipped) > 5:
            sys.stderr.write("  ... and %d more\n" % (len(skipped) - 5))

    candidates, stats = search(
        psms, proteins,
        min_seg=args.min_seg, min_total=args.min_total,
        min_gap=args.min_gap, max_gap=args.max_gap,
        prec_tol=args.prec_tol, gap_tol=args.gap_tol,
        enzymes=(args.enzyme,), specificity=args.specificity,
        max_missed_cleavages=args.max_missed_cleavages,
        require_validation=not args.keep_unvalidated,
        keep_contiguous=args.keep_contiguous,
    )

    n = write_tsv(args.out, candidates)
    targets = sum(1 for c in candidates if not c.protein.is_decoy)
    decoys = len(candidates) - targets

    sys.stderr.write("\n--- run summary ---\n")
    sys.stderr.write("contiguous pairs skipped:           %d\n"
                     % stats["contiguous_skipped"])
    sys.stderr.write("splits passing tiling + gap-length: %d\n" % stats["gap_filtered"])
    sys.stderr.write("tagmatch-validated:                 %d\n" % stats["validated"])
    sys.stderr.write("candidates written:                 %d (%d target, %d decoy)\n"
                     % (n, targets, decoys))
    if candidates:
        best = min(c.qvalue for c in candidates if not c.protein.is_decoy) \
            if targets else None
        if best is not None:
            sys.stderr.write("best target q-value:                %.4f\n" % best)
    sys.stderr.write("written to: %s\n" % args.out)

    if args.unspliced_out:
        rows = unspliced_junction_evidence(candidates, psms)
        m = write_unspliced_tsv(args.unspliced_out, rows)
        sys.stderr.write("unspliced-junction rows:            %d -> %s\n"
                         % (m, args.unspliced_out))

    return 0


if __name__ == "__main__":
    sys.exit(main())
