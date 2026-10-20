#!/usr/bin/env python3
# =============================================================================
# test_depletion_test.py -- proves depletion_test.py computes exactly what
# depl.py computes, on SYNTHETIC random sequences (no real data involved).
#
# depl.py's digest/classification loop is copied here VERBATIM and compared
# with depletion_test.predicted(..., proline_rule=False) on 300 random
# sequences x 3 missed-cleavage settings x random intein positions. It also
# checks the binomial tail against an independent brute-force calculation.
#
# Run from the repository root:  python3 tests/test_depletion_test.py
# =============================================================================
import os
import random
import re
import sys
from math import comb

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "scripts"))
import depletion_test as D   # noqa: E402

rng = random.Random(7)
fails = 0


def depl_py_counts(s, a, b, mc):
    # --- verbatim from depl.py ---
    cuts = [0] + [m.end() for m in re.finditer("[KR]", s) if m.end() < len(s)] + [len(s)]
    c = {"N": 0, "I": 0, "C": 0, "X": 0}
    for i in range(len(cuts) - 1):
        for j in range(i + 1, min(i + mc + 2, len(cuts))):
            st, en = cuts[i] + 1, cuts[j]
            if 6 <= en - st + 1 <= 40:
                c["N" if en < a else "C" if st > b else "I" if st >= a and en <= b else "X"] += 1
    return c


n_cases = 0
for _ in range(300):
    L = rng.randint(300, 1000)
    s = "".join(rng.choice("ACDEFGHIKLMNPQRSTVWYKKRP") for _ in range(L))  # K/R/P enriched
    a = rng.randint(50, L // 2)
    b = rng.randint(a + 50, L - 20)
    for mc in (0, 1, 2):
        n_cases += 1
        if depl_py_counts(s, a, b, mc) != D.predicted(s, a, b, mc, proline_rule=False):
            fails += 1
print("digest identical to depl.py on %d synthetic cases: %s" % (n_cases, "PASS" if not fails else "FAIL"))

# proline rule must suppress exactly the K/R-P sites
s = "AAKPAAKAARPAAR"
got = D.cut_positions(s, True)
ok = got == [0, 7, 14]  # cuts after K7 only (K3 and R10 are followed by P; R14 is the end)
print("proline rule suppresses K/R before P: %s %s" % ("PASS" if ok else "FAIL", got))
fails += not ok

# binomial tail vs brute force
ok = True
for k, n, p in [(0, 6, 0.608), (6, 23, 0.474), (2, 12, 0.587), (4, 27, 0.45)]:
    brute = sum(comb(n, x) * p ** x * (1 - p) ** (n - x) for x in range(0, k + 1))
    ok &= abs(D.binom_tail(k, n, p) - brute) < 1e-15
# published PXD009315 value: RecA 0 of 6 at expected .608 -> P = .0036
ok &= abs(D.binom_tail(0, 6, 0.608) - 0.0036) < 0.00005
print("binomial tail correct (incl. RecA 0/6 at .608 -> .0036): %s" % ("PASS" if ok else "FAIL"))
fails += not ok

print("\n%s: %d failure(s)" % ("ALL PASSED" if not fails else "FAILED", fails))
sys.exit(1 if fails else 0)
