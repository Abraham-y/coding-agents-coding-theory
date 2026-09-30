"""Certified LP + odd-cycle cutting planes: push upper bounds below the published table.

Why this can produce a NEW result. In Weindel-Heckel (arXiv:2504.00613v2) Table 2,
the ILP timed out for s=2, n >= 13, so the published upper bounds there ARE the
plain LP bounds (93 at n=13, 160 at n=14; 24 at s=3, n=13). This repo has already
certified those exact LP values. Lovasz theta at (2,10) showed the LP leaves real
slack (20.47 -> 18.64). Odd-cycle inequalities capture part of that slack at LP
cost and remain certifiable in exact rational arithmetic. Any integer improvement
is a bound the published table does not have.

Valid inequalities used, each with a combinatorial proof:
    clique K:      sum_{v in K} x_v <= 1
    odd cycle C:   sum_{v in C} x_v <= (|C|-1)/2      (C simple, |C| odd >= 5)

Separation of odd cycles is exact via the bipartite double cover: with edge
weight w_uv = 1 - x_u - x_v >= 0, an odd cycle violates its inequality iff its
w-weight is < 1, and shortest odd walks are found by Dijkstra from v-side-0 to
v-side-1. Projected walks may repeat vertices, in which case the inequality
would be INVALID; such walks are discarded rather than repaired.

The float LP only proposes. Certification re-checks, in fractions.Fraction:
  * every cycle constraint is a genuine simple odd cycle of G (edges verified
    against the bitset adjacency, oddness and simplicity re-verified),
  * the dual vector covers every vertex >= 1,
  * the reported bound is floor of the exact dual objective.

Run:
    python search/cutplane.py --n 10 --s 2               # calibration: exact=16, LP=20, theta=18.6
    python search/cutplane.py --n 13 --s 2 --rounds 12   # target: beat published 93
    python search/cutplane.py --verify-cert results/cutcert_n13_s2.json
"""

from __future__ import annotations

import argparse
import json
import math
import os
import sys
import time
from fractions import Fraction

import numpy as np
from scipy.optimize import linprog
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import dijkstra

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from search.graph import load_graph  # noqa: E402
from search.lp_bound import _add_clique_cuts  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS_DIR = os.path.join(ROOT, "results")


# ---------------------------------------------------------------------------
# constraint system
# ---------------------------------------------------------------------------


def base_constraints(g, cliques_max=True):
    names = sorted(g.cliques)
    cons = [("clique", tuple(sorted(g.cliques[t]))) for t in names if len(g.cliques[t]) >= 2]
    if cliques_max:
        raw = [c[1] for c in cons]
        grown, _ = _add_clique_cuts(g, raw, [str(i) for i in range(len(raw))])
        seen = set(raw)
        for k in grown[len(raw):]:
            if k not in seen:
                cons.append(("clique", k))
                seen.add(k)
    return cons


def rhs_of(kind, members):
    if kind == "clique":
        return 1.0
    return (len(members) - 1) / 2.0


def solve_lp(nv, cons):
    rows, cols = [], []
    for i, (_, members) in enumerate(cons):
        for v in members:
            rows.append(i)
            cols.append(v)
    A = csr_matrix((np.ones(len(rows)), (rows, cols)), shape=(len(cons), nv))
    b = np.array([rhs_of(k, m) for k, m in cons])
    res = linprog(c=-np.ones(nv), A_ub=A, b_ub=b, bounds=(0, 1), method="highs")
    if not res.success:
        raise RuntimeError(res.message)
    x = np.asarray(res.x)
    y = -np.asarray(res.ineqlin.marginals)  # dual multipliers, >= 0
    return x, y, -res.fun


# ---------------------------------------------------------------------------
# odd-cycle separation via bipartite double cover
# ---------------------------------------------------------------------------


def edges_of(g):
    E = []
    for v in range(g.num_vertices):
        m = g.adj[v]
        while m:
            low = m & -m
            u = low.bit_length() - 1
            m ^= low
            if u > v:
                E.append((v, u))
    return E


def separate_odd_cycles(g, E, x, max_sources=160, max_cuts=400, tol=1e-6):
    """Return violated simple odd cycles as sorted vertex tuples."""
    nv = g.num_vertices
    w = np.array([max(0.0, 1.0 - x[u] - x[v]) for (u, v) in E])

    # double cover: node (v, p) -> v + p*nv ; edge (u,p)-(v,1-p)
    r, c, d = [], [], []
    for k, (u, v) in enumerate(E):
        for p in (0, 1):
            r.extend([u + p * nv, v + (1 - p) * nv])
            c.extend([v + (1 - p) * nv, u + p * nv])
            d.extend([w[k], w[k]])
    M = csr_matrix((d, (r, c)), shape=(2 * nv, 2 * nv))

    # sources: vertices with largest fractional value (most likely in cuts)
    order = np.argsort(-x)
    sources = [int(v) for v in order[:max_sources] if x[v] > tol]
    if not sources:
        return []

    dist, pred = dijkstra(M, indices=sources, return_predecessors=True)
    cuts = set()
    for si, v in enumerate(sources):
        target = v + nv
        if not np.isfinite(dist[si, target]) or dist[si, target] >= 1.0 - tol:
            continue
        # reconstruct walk v(0) -> v(1), project to G
        path = []
        node = target
        while node != -9999 and node >= 0:
            path.append(node % nv)
            if node == v and len(path) > 1:
                break
            node = pred[si, node]
        cyc = path[:-1]  # last element repeats the source
        if len(cyc) < 5 or len(cyc) % 2 == 0:
            continue
        if len(set(cyc)) != len(cyc):
            continue  # closed walk revisiting a vertex: inequality invalid, discard
        # verify each consecutive pair is an edge (and the wrap-around)
        ok = all((g.adj[cyc[i]] >> cyc[(i + 1) % len(cyc)]) & 1 for i in range(len(cyc)))
        if not ok:
            continue
        if sum(x[u] for u in cyc) <= (len(cyc) - 1) / 2.0 + tol:
            continue  # not actually violated
        cuts.add(tuple(sorted(cyc)) + ("|",) + tuple(cyc))  # keep order for cert
        if len(cuts) >= max_cuts:
            break
    # return ordered cycles, deduped by vertex set
    out, seen = [], set()
    for item in cuts:
        sep = item.index("|")
        key = item[:sep]
        if key in seen:
            continue
        seen.add(key)
        out.append(tuple(item[sep + 1:]))
    return out


# ---------------------------------------------------------------------------
# exact certification
# ---------------------------------------------------------------------------


def certify(g, cons, y_float, denom=10 ** 6):
    """Exact rational check. Returns (bound Fraction, floored int)."""
    nv = g.num_vertices
    y = [Fraction(max(0, int(math.ceil(val * denom))), denom) for val in y_float]

    # re-verify every cycle constraint combinatorially before trusting its RHS
    for (kind, members) in cons:
        if kind != "cycle":
            continue
        L = len(members)
        assert L >= 5 and L % 2 == 1, "cycle length invalid"
        assert len(set(members)) == L, "cycle not simple"
        for i in range(L):
            u, v = members[i], members[(i + 1) % L]
            assert (g.adj[u] >> v) & 1, "cycle edge missing from G"

    cover = [Fraction(0)] * nv
    for i, (_, members) in enumerate(cons):
        if y[i] == 0:
            continue
        for v in set(members):
            cover[v] += y[i]
    worst = min(cover)
    if worst <= 0:
        raise RuntimeError("a vertex has zero dual coverage; cannot repair")
    if worst < 1:
        scale = 1 / worst
        y = [val * scale for val in y]

    bound = Fraction(0)
    for i, (kind, members) in enumerate(cons):
        if y[i] == 0:
            continue
        b = Fraction(1) if kind == "clique" else Fraction(len(members) - 1, 2)
        bound += y[i] * b
    return y, bound, math.floor(bound)


def verify_certificate(path):
    with open(path) as fh:
        cert = json.load(fh)
    g = load_graph(cert["n"], cert["s"])
    cons = [(k, tuple(m)) for k, m in cert["constraints"]]
    y = [Fraction(a, b) for a, b in cert["y"]]
    nv = g.num_vertices
    # combinatorial validity
    for (kind, members) in cons:
        if kind == "cycle":
            L = len(members)
            if L < 5 or L % 2 == 0 or len(set(members)) != L:
                return False, "invalid cycle constraint"
            for i in range(L):
                if not (g.adj[members[i]] >> members[(i + 1) % L]) & 1:
                    return False, "cycle edge missing"
    cover = [Fraction(0)] * nv
    for i, (_, members) in enumerate(cons):
        for v in set(members):
            cover[v] += y[i]
    if min(cover) < 1:
        return False, f"dual coverage {min(cover)} < 1"
    bound = sum(
        y[i] * (Fraction(1) if k == "clique" else Fraction(len(m) - 1, 2))
        for i, (k, m) in enumerate(cons)
    )
    if bound != Fraction(*cert["bound"]):
        return False, "objective mismatch"
    if math.floor(bound) != cert["certified_upper_bound"]:
        return False, "floor mismatch"
    return True, f"OK: A_{cert['s']}({cert['n']}) <= {cert['certified_upper_bound']} (exact {bound})"


# ---------------------------------------------------------------------------


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int)
    ap.add_argument("--s", type=int)
    ap.add_argument("--rounds", type=int, default=10)
    ap.add_argument("--sources", type=int, default=160)
    ap.add_argument("--verify-cert", type=str, default=None)
    args = ap.parse_args()

    if args.verify_cert:
        ok, msg = verify_certificate(args.verify_cert)
        print(("PASS " if ok else "FAIL ") + msg)
        return 0 if ok else 1

    t0 = time.time()
    g = load_graph(args.n, args.s)
    E = edges_of(g)
    cons = base_constraints(g)
    print(f"G({args.n},{args.s}): {g.num_vertices} vertices, {len(E)} edges, "
          f"{len(cons)} base clique constraints")

    history = []
    for rnd in range(args.rounds):
        x, y, obj = solve_lp(g.num_vertices, cons)
        cuts = separate_odd_cycles(g, E, x, max_sources=args.sources)
        history.append({"round": rnd, "lp_float": obj, "new_cuts": len(cuts),
                        "constraints": len(cons)})
        print(f"  round {rnd}: LP={obj:.5f}, violated odd cycles found={len(cuts)}")
        if not cuts:
            break
        for cyc in cuts:
            cons.append(("cycle", cyc))

    x, y, obj = solve_lp(g.num_vertices, cons)
    y_exact, bound, ub = certify(g, cons, y)
    elapsed = time.time() - t0
    print(f"final float LP (NOT a result): {obj:.6f}")
    print(f"certified exact dual objective: {bound} (= {float(bound):.6f})")
    print(f"CERTIFIED: A_{args.s}({args.n}) <= {ub}")
    print(f"elapsed {elapsed:.1f}s")

    os.makedirs(RESULTS_DIR, exist_ok=True)
    path = os.path.join(RESULTS_DIR, f"cutcert_n{args.n}_s{args.s}.json")
    with open(path, "w") as fh:
        json.dump({
            "method": "LP + odd-cycle cutting planes, dual certified in exact rationals",
            "n": args.n, "s": args.s,
            "constraints": [[k, list(m)] for k, m in cons],
            "y": [[f.numerator, f.denominator] for f in y_exact],
            "bound": [bound.numerator, bound.denominator],
            "certified_upper_bound": ub,
            "rounds": history,
            "elapsed_seconds": round(elapsed, 2),
        }, fh)
    print(f"certificate: {path}")
    ok, msg = verify_certificate(path)
    print(("re-check PASS " if ok else "re-check FAIL ") + msg)
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
