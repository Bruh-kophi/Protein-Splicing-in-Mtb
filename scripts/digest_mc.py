"""Expected intein proportion at different missed-cleavage limits.

Reproduces the depletion-test expectation from depl.py, but sweeps the
missed-cleavage limit so the >=2 vs >=3 question can be answered directly.
"""
import sys, itertools

FASTA = sys.argv[1]
INTEINS = {"RecA": ("P9WHJ3", 252, 691), "SufB": ("P9WFP7", 253, 611), "DnaB": ("P9WMR3", 400, 815)}

def load(path):
    seqs, acc, buf = {}, None, []
    for line in open(path):
        if line.startswith(">"):
            if acc: seqs[acc] = "".join(buf)
            acc = line.split("|")[1] if "|" in line else line[1:].split()[0]
            buf = []
        else:
            buf.append(line.strip())
    if acc: seqs[acc] = "".join(buf)
    return seqs

def cut_sites(seq, proline_rule=True):
    """0-based slice boundaries: 0, after each K/R (not before P), and len."""
    s = [0]
    for i, aa in enumerate(seq):
        if aa in "KR":
            if proline_rule and i + 1 < len(seq) and seq[i + 1] == "P":
                continue
            s.append(i + 1)
    if s[-1] != len(seq): s.append(len(seq))
    return s

def peptides(seq, max_mc, lo=6, hi=40, proline_rule=True):
    s = cut_sites(seq, proline_rule)
    out = set()
    for i in range(len(s) - 1):
        for j in range(i + 1, min(i + 2 + max_mc, len(s))):
            a, b = s[i], s[j]
            if lo <= b - a <= hi:
                out.add((a + 1, b))          # 1-based inclusive
    return out

def classify(start, end, a, b):
    if end < a:   return "N"
    if start > b: return "C"
    if start >= a and end <= b: return "I"
    return "B"                                # boundary-spanning

seqs = load(FASTA)
print(f"{'protein':7} {'MC':>3} {'N':>5} {'I':>5} {'C':>5} {'B':>5} {'expected proportion':>21}")
print("-" * 58)
for name, (acc, a, b) in INTEINS.items():
    seq = seqs[acc]
    for mc in (2, 3):
        peps = peptides(seq, mc)
        t = {"N": 0, "I": 0, "C": 0, "B": 0}
        for st, en in peps: t[classify(st, en, a, b)] += 1
        denom = t["N"] + t["I"] + t["C"]
        print(f"{name:7} {mc:>3} {t['N']:>5} {t['I']:>5} {t['C']:>5} {t['B']:>5} {t['I']/denom:>20.4f}")
    print(f"{'':7} {'len':>3} intein {b-a+1} aa of {len(seq)} aa = length-based {(b-a+1)/len(seq):.4f}")
