"""Faithful transcription of the MNSX Appendix A program. NOT my algorithm.

PROVENANCE, stated up front. The algorithm and code structure in this file are
transcribed from Appendix A of

    F. Magniez, A. Nayak, M. Santha, J. Sherman, G. Tardos, D. Xiao,
    "Improved bounds for the randomized decision tree complexity of recursive
    majority", arXiv:1309.7565v1, Appendix A, pages 25-28.

The authors also published it at https://www.dropbox.com/s/wcrdoib5h918p2e/
commented-majority.py, which now returns HTTP 403; the paper's appendix is the
surviving copy. Variable names, the record layout and the recursion are theirs.
What is mine here: exact-Fraction reporting, the published-value assertions, and
the comments marked [note].

WHY TRANSCRIBE IT AT ALL, given search/maj3_fast.py is an independent solver?
Two reasons:

1. Speed. Their representation is drastically more compact than mine. A stable
   configuration is a 6-integer record
       r[0] = number of evaluations giving 0
       r[1] = number of evaluations giving 1
       r[2] = number of queried sensitive bits over evaluations with correct output
       r[3:6] = indices of the three children, sorted
   built by a triple loop over sorted child indices with the constraint that at
   most one child is the special index-0 "fully queried subtree". No tree objects,
   no per-state leaf arrays. My own solver visited 25M+ states at k=4 and was
   still on iteration 0 after 3.4 hours; theirs is a list of 246,792 records.

2. Cross-validation in BOTH directions. My independent solver reproducing their
   published alpha_1..alpha_3 already checks my code against their PAPER. Running
   their code and getting the same numbers as my solver checks the transcription
   against my independent derivation. Agreement between two implementations that
   share no representation is the same discipline as verify.py vs verify2.py.

They work with recursive NOT-3-MAJ (a negated majority gate at every node) rather
than 3-MAJ, which is why the parity terms `depth%2` and `(depth-i)%2` appear: it
makes the minority path have a constant value instead of an alternating one, so
no dual configurations are needed.

Run:
    python search/maj3_mnsx.py --k 4
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time
from fractions import Fraction

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS_DIR = os.path.join(ROOT, "results")

PUBLISHED_ALPHA = {
    1: Fraction(2),
    2: Fraction(24, 7),
    3: Fraction(12231, 2203),
    4: Fraction(2027349, 216164),
}

# globals, as in the original
depth = 0
g: list = []
d: list = []
alpha: list = [0, 1]
ot: list = []
opt: list = []


def prep(dep: int):
    """Build g[i], the list of depth-i stable configurations, and the index maps d[i]."""
    global depth, g, d
    depth = dep
    # two stable configurations at depth 0: a variable set to 1, and an unqueried variable
    g = [[[0, 1, depth % 2], [1, 1, 0]]]
    d = [{}]
    for i in range(1, depth + 1):
        gl = g[-1]
        # index 0 of every level: a fully queried subtree of depth i evaluating to 1
        gg = [[0, 1, ((depth - i) % 2) * 2 ** i]]
        dd = {}
        for a in range(len(gl)):
            for b in range(max(a, 1), len(gl)):
                for c in range(b, len(gl)):
                    # a <= b <= c, and no two fully evaluated siblings
                    dd[(a, b, c)] = len(gg)
                    A, B, C = gl[a], gl[b], gl[c]
                    gg.append([
                        A[0] * B[1] * C[1] + A[1] * B[0] * C[1] + A[1] * B[1] * C[0],
                        A[1] * B[0] * C[0] + A[0] * B[1] * C[0] + A[0] * B[0] * C[1],
                        A[0] * B[1] * C[2] + A[0] * B[2] * C[1]
                        + A[1] * B[0] * C[2] + A[1] * B[2] * C[0]
                        + A[2] * B[0] * C[1] + A[2] * B[1] * C[0],
                        a, b, c,
                    ])
        g.append(gg)
        d.append(dd)


def adj(i: int, a: int, s: list, t):
    """Crawl the configuration tree, updating `opt` when querying a variable beats it."""
    global ot, opt
    if s == []:
        tt = [0, 2 ** i, 1, 0, 0]
    else:
        s1, s2 = g[i][s[0][0]], g[i][s[0][1]]
        nn = s1[1] * s2[1]
        ne = s1[0] * s2[1] + s2[0] * s1[1]
        tt = [nn * t[0] + ne * t[3], nn * t[1] + ne * t[4], nn * t[2]]
        if s[0][0] == 0:
            tt.extend([s2[0] * t[0], s2[0] * t[1]])
        else:
            w, j = 0, i
            for si in s:
                j = j + 1
                if w <= si[0]:
                    triple = (w, si[0], si[1])
                elif w <= si[1]:
                    triple = (si[0], w, si[1])
                else:
                    triple = (si[0], si[1], w)
                w = d[j][triple]
            tt.extend(ot[w])
    if i > 0:
        A = g[i][a]
        if A[3] > 0:
            adj(i - 1, A[3], [A[4:6]] + s, tt)
        adj(i - 1, A[4], [[A[3], A[5]]] + s, tt)
        adj(i - 1, A[5], [A[3:5]] + s, tt)
    else:
        pair = [tt[0] + tt[2] + tt[3], tt[1] + tt[4]]
        if (pair[0] - opt[0]) * alpha[0] * 2 ** depth < (pair[1] - opt[1]) * alpha[1]:
            opt = pair


def check(a0: int, a1: int):
    """For alpha = a0/a1, build ot[] and return the optimum record at the root."""
    global alpha, ot, opt
    alpha = [a0, a1]
    ot = [[0, 0]]  # g[depth][0] is impossible
    for a in range(1, len(g[-1])):
        opt = [0, g[-1][a][2]]
        adj(depth, a, [], 0)
        ot.append(opt)
    return ot[-1]


def compute_alpha(dep: int, verbose: bool = True):
    """Fixed-point iteration on alpha. Returns (Fraction, iterations)."""
    prep(dep)
    p, q = 0, 1
    x = check(p, q)
    f = Fraction(0)
    iters = 0
    history = []
    while x != [0, 0]:
        iters += 1
        f = Fraction(x[1], x[0] * (2 ** dep))
        p, q = f.numerator, f.denominator
        history.append({"iter": iters, "alpha": str(f), "ot_last": list(x)})
        if verbose:
            print(f"  iter {iters}: alpha -> {f}")
        x = check(p, q)
        if iters > 60:
            raise RuntimeError("alpha iteration did not converge")
    return f, history


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--k", type=int, required=True)
    args = ap.parse_args()

    sys.setrecursionlimit(100000)
    t0 = time.time()
    print(f"k={args.k}: MNSX Appendix A transcription")
    prep(args.k)
    print(f"stable configurations built: " +
          ", ".join(f"|g[{i}]|={len(g[i])}" for i in range(len(g))))

    a, history = compute_alpha(args.k)
    elapsed = time.time() - t0

    published = PUBLISHED_ALPHA.get(args.k)
    match = published is not None and a == published
    print(f"alpha_{args.k} = {a}")
    if published is not None:
        print(f"published  = {published}  -> {'MATCH' if match else 'MISMATCH'}")
    base = 2 + float(a) ** (-1.0 / args.k)
    print(f"lower-bound base = {base:.9f}   [float, display only]")
    print(f"elapsed {elapsed:.1f}s")

    os.makedirs(RESULTS_DIR, exist_ok=True)
    path = os.path.join(RESULTS_DIR, f"maj3_mnsx_k{args.k}.json")
    with open(path, "w") as fh:
        json.dump({
            "source": "transcription of arXiv:1309.7565v1 Appendix A",
            "k": args.k,
            "alpha_str": str(a),
            "alpha_numerator": a.numerator,
            "alpha_denominator": a.denominator,
            "published": str(published) if published else None,
            "matches_published": match,
            "stable_config_counts": [len(x) for x in g],
            "base_float_display_only": base,
            "elapsed_seconds": round(elapsed, 2),
            "iterations": history,
        }, fh, indent=2)
    print(f"log: {path}")
    return 0 if (published is None or match) else 1


if __name__ == "__main__":
    raise SystemExit(main())
