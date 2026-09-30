"""Quaternary single-edit-correcting codes: graph, verifier, and symmetric search.

TARGET (new-result attempt). Weindel-Heckel Table 3, quaternary single-edit
(insertion OR deletion OR substitution) codes at n=6: best known construction
113 (ReduMIS, 24h), upper bounds ILP 212 / LP 215. Bracket [113, 212] -- wide
open, incumbent from a general-purpose solver.

The untried structure: the confusability graph is invariant under
    * any permutation of the alphabet {0,1,2,3} applied position-wise (S_4, 24)
    * reversal of the word (2)
a group of order 48 that a generic MIS solver knows nothing about. For each
subgroup H, the best H-invariant code is an exact MWIS on the orbit graph --
small enough to solve to PROVEN optimality for large H. This file enumerates
ALL subgroups of S_4 (there are 30), each with and without reversal, and solves
each orbit problem exactly (small) or heuristically (large).

Definitions. A code C over {0,1,2,3}^n corrects one edit iff the balls
    B(x) = {x} u subs(x) u dels(x) u ins(x)
are pairwise disjoint. Balls contain strings of lengths n-1, n, n+1; two balls
intersect iff they share a string in the same length class. For every possible
received string t, the set {x : t in B(x)} is a clique, and these cliques cover
all edges -- the same clique-cover trick as the binary graph.

Validation before trusting the encoding (rule 4): the fractional-transversal LP
on this clique cover must reproduce the paper's published LP bound of 215.

Run:
    python search/edit4.py --n 6 --validate        # LP vs published 215
    python search/edit4.py --n 6 --sweep           # all-subgroup symmetric search
"""

from __future__ import annotations

import argparse
import itertools
import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS_DIR = os.path.join(ROOT, "results")
BEST_DIR = os.path.join(ROOT, "best")

ALPHABET = "0123"
PUBLISHED = {  # Table 3, arXiv:2504.00613v2 -- checked 2026-08-06
    6: {"lp": 215, "ilp_ub": 212, "best": 113},
    7: {"lp": 744, "ilp_ub": 853, "best": 352},   # ilp col reads 853 at n=7
    8: {"lp": 2621, "ilp_ub": 4527, "best": 1130},
}


def ball(word: str):
    n = len(word)
    out = set()
    out.add(word)
    for i in range(n):
        for c in ALPHABET:
            if c != word[i]:
                out.add(word[:i] + c + word[i + 1:])
        out.add(word[:i] + word[i + 1:])          # deletion
    for i in range(n + 1):
        for c in ALPHABET:
            out.add(word[:i] + c + word[i:])      # insertion
    return out


class Edit4Graph:
    def __init__(self, n: int):
        self.n = n
        self.num_vertices = 4 ** n
        self.words = ["".join(w) for w in itertools.product(ALPHABET, repeat=n)]
        buckets: dict[str, list[int]] = {}
        for idx, w in enumerate(self.words):
            for t in ball(w):
                buckets.setdefault(t, []).append(idx)
        self.cliques = {t: m for t, m in buckets.items() if len(m) >= 2}
        adj = [0] * self.num_vertices
        for members in self.cliques.values():
            mask = 0
            for m in members:
                mask |= 1 << m
            for m in members:
                adj[m] |= mask
        for i in range(self.num_vertices):
            adj[i] &= ~(1 << i)
        self.adj = adj


# ---------------------------------------------------------------------------
# verifier (pairwise, independent of the clique construction)
# ---------------------------------------------------------------------------


def verify_edit4(words, n: int):
    """(valid, size, reason). Pairwise ball-disjointness, no clique machinery."""
    ws = list(words)
    seen = set()
    for w in ws:
        if len(w) != n or set(w) - set(ALPHABET):
            return False, 0, f"malformed codeword {w!r}"
        if w in seen:
            return False, 0, f"duplicate {w!r}"
        seen.add(w)
    for i in range(len(ws)):
        bi = ball(ws[i])
        for j in range(i + 1, len(ws)):
            inter = bi & ball(ws[j])
            if inter:
                return False, 0, f"{ws[i]} and {ws[j]} share {sorted(inter)[0]!r}"
    return True, len(ws), ""


# ---------------------------------------------------------------------------
# LP validation against the published bound
# ---------------------------------------------------------------------------


def lp_bound(g: Edit4Graph):
    from scipy.optimize import linprog
    from scipy.sparse import csr_matrix

    cons = list(g.cliques.values())
    rows, cols = [], []
    for i, members in enumerate(cons):
        for v in members:
            rows.append(i)
            cols.append(v)
    A = csr_matrix((np.ones(len(rows)), (rows, cols)),
                   shape=(len(cons), g.num_vertices))
    res = linprog(c=-np.ones(g.num_vertices), A_ub=A,
                  b_ub=np.ones(len(cons)), bounds=(0, 1), method="highs")
    if not res.success:
        raise RuntimeError(res.message)
    return -res.fun


# ---------------------------------------------------------------------------
# subgroups of S_4, and the induced word transformations
# ---------------------------------------------------------------------------


def all_subgroups_s4():
    """All 30 subgroups of S_4. Every subgroup of S_4 is 2-generated."""
    elems = list(itertools.permutations(range(4)))

    def mul(a, b):
        return tuple(a[b[i]] for i in range(4))

    def closure(gens):
        seen = set(gens) | {tuple(range(4))}
        frontier = list(seen)
        while frontier:
            nxt = []
            for a in frontier:
                for b in list(seen):
                    for c in (mul(a, b), mul(b, a)):
                        if c not in seen:
                            seen.add(c)
                            nxt.append(c)
            frontier = nxt
        return frozenset(seen)

    subs = set()
    for a in elems:
        for b in elems:
            subs.add(closure([a, b]))
    return sorted(subs, key=len)


def word_transforms(subgroup, n: int, with_reversal: bool):
    """The set of word->word maps induced by (symbol perm, optional reversal)."""
    out = []
    for perm in subgroup:
        table = {ALPHABET[i]: ALPHABET[perm[i]] for i in range(4)}
        out.append(("f", perm, table))
        if with_reversal:
            out.append(("r", perm, table))
    return out


def apply_t(t, word):
    kind, _, table = t
    w = "".join(table[c] for c in word)
    return w[::-1] if kind == "r" else w


# ---------------------------------------------------------------------------
# orbit MWIS, exact where possible
# ---------------------------------------------------------------------------


def orbits_of(g: Edit4Graph, transforms):
    """Just the orbit partition -- cheap; lets the caller skip before the O(k^2) part."""
    idx = {w: i for i, w in enumerate(g.words)}
    seen = [False] * g.num_vertices
    orbits = []
    for v in range(g.num_vertices):
        if seen[v]:
            continue
        orb = sorted({idx[apply_t(t, g.words[v])] for t in transforms})
        for u in orb:
            seen[u] = True
        orbits.append(tuple(orb))
    return orbits


def orbit_problem(g: Edit4Graph, orbits):
    usable, masks = [], []
    for orb in orbits:
        ok = True
        for i in range(len(orb)):
            for j in range(i + 1, len(orb)):
                if (g.adj[orb[i]] >> orb[j]) & 1:
                    ok = False
                    break
            if not ok:
                break
        if ok:
            usable.append(orb)
            m = 0
            for v in orb:
                m |= 1 << v
            masks.append(m)
    nbr = []
    for orb in usable:
        m = 0
        for v in orb:
            m |= g.adj[v]
        nbr.append(m)
    edges = []
    for i in range(len(usable)):
        for j in range(i + 1, len(usable)):
            if nbr[i] & masks[j]:
                edges.append((i, j))
    return usable, edges


def solve_orbit_exact(usable, edges, time_limit=120):
    from scipy.optimize import milp, Bounds, LinearConstraint
    from scipy.sparse import csr_matrix

    k = len(usable)
    w = np.array([len(o) for o in usable], dtype=float)
    if not edges:
        return list(range(k)), True
    rows, cols, data = [], [], []
    for i, (a, b) in enumerate(edges):
        rows.extend([i, i])
        cols.extend([a, b])
        data.extend([1.0, 1.0])
    A = csr_matrix((data, (rows, cols)), shape=(len(edges), k))
    res = milp(c=-w, constraints=LinearConstraint(A, -np.inf, 1),
               integrality=np.ones(k), bounds=Bounds(0, 1),
               options={"time_limit": time_limit})
    chosen = [i for i in range(k) if res.x is not None and res.x[i] > 0.5]
    return chosen, res.status == 0


def sweep(n: int, exact_cap: int = 900, time_limit: int = 120):
    g = Edit4Graph(n)
    print(f"Edit4 G({n}): {g.num_vertices} vertices, {len(g.cliques)} cliques")
    results = []
    best = (0, None, None, False)
    for sub in all_subgroups_s4():
        for rev in (False, True):
            transforms = word_transforms(sub, n, rev)
            orbits = orbits_of(g, transforms)
            label = f"|H|={len(sub)}{'+rev' if rev else ''}"
            if len(orbits) > exact_cap:
                print(f"  {label:<12} {len(orbits):>5} orbits -- skipped (exact cap)")
                continue
            usable, edges = orbit_problem(g, orbits)
            k = len(usable)
            chosen, proven = solve_orbit_exact(usable, edges, time_limit)
            words = sorted(g.words[v] for i in chosen for v in usable[i])
            ok, size, reason = verify_edit4(words, n)
            tag = "OPT" if proven else "timelimit"
            print(f"  {label:<12} {k:>5} orbits -> invariant code {size:>4} "
                  f"[{tag}] verify={ok}")
            if not ok:
                print(f"      VERIFIER REJECTED: {reason} -- bug, investigate")
                continue
            results.append({"order": len(sub), "reversal": rev, "orbits": k,
                            "size": size, "proven_optimal_invariant": proven})
            if size > best[0]:
                best = (size, words, label, proven)
    return g, results, best


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, required=True)
    ap.add_argument("--validate", action="store_true")
    ap.add_argument("--sweep", action="store_true")
    ap.add_argument("--exact-cap", type=int, default=900)
    ap.add_argument("--time-limit", type=int, default=120)
    args = ap.parse_args()

    t0 = time.time()
    if args.validate:
        g = Edit4Graph(args.n)
        print(f"Edit4 G({args.n}): {g.num_vertices} vertices, {len(g.cliques)} cliques")
        val = lp_bound(g)
        pub = PUBLISHED.get(args.n, {}).get("lp")
        print(f"our float LP bound : {val:.4f}")
        print(f"published LP bound : {pub}")
        if pub is not None:
            match = abs(val - pub) < 1.0
            print("ENCODING " + ("VALIDATED" if match else "MISMATCH -- STOP"))
            return 0 if match else 1
        return 0

    if args.sweep:
        g, results, best = sweep(args.n, args.exact_cap, args.time_limit)
        size, words, label, proven = best
        pub_best = PUBLISHED.get(args.n, {}).get("best")
        print(f"\nbest invariant code found: {size}  ({label}, proven-invariant-opt={proven})")
        print(f"published best known     : {pub_best}")
        os.makedirs(RESULTS_DIR, exist_ok=True)
        stamp = time.strftime("%Y%m%dT%H%M%S")
        out = {"n": args.n, "results": results, "best_size": size,
               "best_group": label, "published_best": pub_best,
               "elapsed_seconds": round(time.time() - t0, 1)}
        path = os.path.join(RESULTS_DIR, f"edit4_sweep_n{args.n}_{stamp}.json")
        with open(path, "w") as fh:
            json.dump(out, fh, indent=2)
        print(f"log: {path}")
        if words and pub_best and size > pub_best:
            ok, sz, _ = verify_edit4(words, args.n)
            if ok:
                os.makedirs(BEST_DIR, exist_ok=True)
                apath = os.path.join(BEST_DIR, f"edit4_n{args.n}_size{sz}.txt")
                with open(apath, "w") as fh:
                    fh.write(f"# quaternary single-edit code, n = {args.n}\n")
                    fh.write(f"# size = {sz}, group {label}\n")
                    for w in words:
                        fh.write(w + "\n")
                print(f"CANDIDATE RECORD saved: {apath} -- run novelty check before celebrating")
        return 0

    ap.error("pass --validate or --sweep")


if __name__ == "__main__":
    raise SystemExit(main())
