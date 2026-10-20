# =============================================================================
# annot.py -- fragment-ion (b/y) check of one spectrum against one peptide
#
# THESIS: Methods 3.8 (evidence standard); Results 4.3 (Table 4.3: SufB
#         586-611 intein-end peptide, bonds covered in five runs).
#
# QUESTION
#   Does the raw spectrum really support the peptide the search engine
#   reported -- in particular its END, which for an intein-end peptide is
#   the whole point?
#
# HOW
#   1. Find the spectrum with the given scan number in an mzML file and
#      decode its m/z and intensity arrays (base64, optionally zlib).
#   2. Calculate every b ion (fragment holding the peptide's START) and y ion
#      (fragment holding its END), charge 1+ and 2+, from monoisotopic
#      residue masses. Cysteine is carbamidomethylated (C = 160.03065).
#   3. Match each ion to the nearest peak within 20 ppm.
#   4. Report: b and y ions found; how many of the n-1 backbone bonds have
#      at least one b or y ion; % of total intensity explained; and what the
#      12 most intense peaks are ("?" = unassigned).
#
# IN PLAIN WORDS
#   A peptide breaks into a ladder of fragments in the instrument. If the
#   fragments that contain the last few amino acids are really there (for
#   SufB 586-611, y2 = His-Asn), the peptide's end is confirmed.
#
# INPUT   argv[1] mzML file; argv[2] scan number; argv[3] plain peptide
#         sequence (upper case, no modification text)
# OUTPUT  printed summary
# RUN     python3 annot.py 0686Mass_S7.mzML 34355 EAHEPVYNLDVENPDSYLAYGFAVHN
#
# LIMITS (state them with any result)
#   * Only b and y ions at 1+ and 2+. No isotope peaks, neutral losses,
#     a ions, internal fragments or higher charges, so the explained
#     intensity is a MINIMUM.
#   * No variable modifications (e.g. oxidised M is not modelled).
#   * This is a count, not a probability score.
#   * Finds the scan by the text 'scan=<N>"' in the <spectrum> line, which
#     suits Thermo-style native IDs (PXD009315). For diaTracer pseudo-spectra
#     (Dataset 3), whose ids are 'frame=<N>', pass the whole id instead,
#     e.g.  python3 annot.py run_diatracer.mzML frame=69390 PEPTIDE
#     (any argument containing '=' is matched as id="<argument>").
#   * The ion series index "b5" means the first 5 residues; "y5" the last 5.
#     A bond after residue i is covered by b_i or by y_(n-i).
#
# NEEDS   Python 3 standard library only.
# =============================================================================
import sys, re, base64, zlib, struct
f, scan, pep = sys.argv[1], sys.argv[2], sys.argv[3]
# monoisotopic residue masses (Da); C includes carbamidomethylation (+57.02146)
M = dict(G=57.02146,A=71.03711,S=87.03203,P=97.05276,V=99.06841,T=101.04768,C=160.03065,L=113.08406,I=113.08406,N=114.04293,D=115.02694,Q=128.05858,K=128.09496,E=129.04259,M=131.04049,H=137.05891,F=147.06841,R=156.10111,Y=163.06333,W=186.07931)
# --- pull out the text block of the requested spectrum ----------------------
blk, on = [], False
for line in open(f):
    if not on and '<spectrum ' in line and (('id="%s"' % scan) if '=' in scan else ('scan=%s"' % scan)) in line: on = True
    if on:
        blk.append(line)
        if '</spectrum>' in line: break
# --- decode the binary arrays: 'mz' (m/z) and 'it' (intensity) ---------------
arr = {}
for a in ''.join(blk).split('<binaryDataArray ')[1:]:
    raw = base64.b64decode(re.search(r'<binary>([^<]*)</binary>', a).group(1))
    if 'zlib' in a: raw = zlib.decompress(raw)
    t = 'd' if '64-bit' in a else 'f'          # 64-bit double or 32-bit float
    arr['mz' if 'm/z array' in a else 'it'] = struct.unpack('<%d%s' % (len(raw) // struct.calcsize(t), t), raw)
# peaks as (m/z, intensity), sorted by m/z
pk = sorted(zip(arr['mz'], arr['it']))
# P = proton mass, W = water (y ions carry the C-terminal OH plus H)
n, P, W = len(pep), 1.007276, 18.010565
# --- theoretical b and y ions, 1+ and 2+ ------------------------------------
ions = {}
for i in range(1, n):
    for z in (1, 2):
        ions['b%d%s' % (i, '+' * z)] = (sum(M[c] for c in pep[:i]) + z * P) / z
        ions['y%d%s' % (i, '+' * z)] = (sum(M[c] for c in pep[-i:]) + W + z * P) / z
# --- match each ion to the closest peak within 20 ppm -----------------------
lab, used = {}, set()          # lab: ion name -> peak index; used: matched peaks
for k, m in ions.items():
    c = [(abs(x - m) / m * 1e6, j) for j, (x, y) in enumerate(pk) if abs(x - m) / m * 1e6 <= 20]
    if c: lab[k] = min(c)[1]; used.add(min(c)[1])
tot = sum(y for x, y in pk)                     # total ion intensity
# ion series numbers found (charge state ignored)
b = sorted({int(re.sub(r'\D', '', k)) for k in lab if k[0] == 'b'})
y = sorted({int(re.sub(r'\D', '', k)) for k in lab if k[0] == 'y'})
# bonds covered: bond i is covered by b_i, or by y_j with j = n - i
cov = len(set(b) | {n - i for i in y})
print("scan %s: %d peaks; b ions %d of %d; y ions %d of %d" % (scan, len(pk), len(b), n - 1, len(y), n - 1))
print("bonds with a b or y ion: %d of %d; intensity explained: %.0f%%" % (cov, n - 1, 100 * sum(pk[j][1] for j in used) / tot))
print("b:", b); print("y:", y)
# --- what are the 12 strongest peaks? ----------------------------------------
inv = {}
for k, j in lab.items(): inv.setdefault(j, []).append(k)
print("12 strongest peaks:")
for j in sorted(range(len(pk)), key=lambda j: -pk[j][1])[:12]:
    print("  %10.4f  %5.1f%%  %s" % (pk[j][0], 100 * pk[j][1] / tot, ','.join(inv.get(j, ['?']))))
