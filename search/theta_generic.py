"""Lovasz theta of an arbitrary graph in the mis_ils text format (k m / weights /
edge list), as an UPPER bound on the independence number of an orbit graph.

Use: if theta(orbit graph) < 31 for a free-acting order-4 class at n=6, then no
code in that class has 31 orbits, i.e. 124 is impossible there and search in
that class can stop at 120.

FLOAT ONLY unless --certify. The SDP solver proposes; the number printed without
--certify is guidance for where to spend search time, never a verdict (CLAUDE.md
rule 10). With --certify the dual (t, A) is rounded to rationals and t*I + A - J
is checked PSD by exact rational LDL^T, reusing search/theta.py's routine, and
the certified bound is floor(t).

    python search/theta_generic.py graph.txt [--solver SCS] [--certify]
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time
from fractions import Fraction

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


def read_graph(path):
    with open(path) as fh:
        k, m = map(int, fh.readline().split())
        w = list(map(int, fh.readline().split()))
        E = set()
        for line in fh:
            a, b = map(int, line.split())
            if a != b:
                E.add((min(a, b), max(a, b)))
    return k, w, sorted(E)


def theta_primal(k, E, solver, verbose):
    import cvxpy as cp
    X = cp.Variable((k, k), PSD=True)
    rows = np.array([a for a, _ in E])
    cols = np.array([b for _, b in E])
    cons = [cp.trace(X) == 1, X[rows, cols] == 0]
    prob = cp.Problem(cp.Maximize(cp.sum(X)), cons)
    prob.solve(solver=solver, verbose=verbose)
    return prob.value, X.value, cons


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("graph")
    ap.add_argument("--solver", default="SCS")
    ap.add_argument("--verbose", action="store_true")
    ap.add_argument("--certify", action="store_true")
    ap.add_argument("--denom", type=int, default=10 ** 6)
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    k, w, E = read_graph(a.graph)
    if any(x != w[0] for x in w):
        print("WARNING: weighted graph; theta bounds the UNWEIGHTED independence number")
    t0 = time.time()
    val, X, cons = theta_primal(k, E, a.solver, a.verbose)
    dt = time.time() - t0
    print(f"graph {a.graph}: {k} vertices, {len(E)} edges")
    print(f"theta (float, {a.solver}) = {val:.6f}   [{dt:.0f}s]")
    print(f"=> alpha <= {int(np.floor(val + 1e-6))} (FLOAT, not certified)")
    rec = {"graph": a.graph, "k": k, "m": len(E), "theta_float": float(val),
           "solver": a.solver, "seconds": round(dt, 1)}
    if a.certify:
        from search.theta import is_psd_exact
        # dual: theta = min t s.t. t*I + A - J psd, A supported on edges.
        # cvxpy gives the dual multipliers of the edge constraints as A.
        t_f = float(val)
        dual = cons[1].dual_value
        A = np.zeros((k, k))
        for (i, j), d in zip(E, np.ravel(dual)):
            A[i, j] = A[j, i] = float(d)
        # sign convention check: try both signs, pick the one that is nearly psd
        best = None
        for sgn in (1.0, -1.0):
            M = t_f * np.eye(k) + sgn * A - np.ones((k, k))
            lam = np.linalg.eigvalsh(M)[0]
            if best is None or lam > best[0]:
                best = (lam, sgn)
        lam, sgn = best
        print(f"dual check: min eigenvalue of t*I + A - J = {lam:.3e} (sign {sgn:+.0f})")
        den = a.denom
        Aq = [[Fraction(int(round(sgn * A[i][j] * den)), den) if i != j else Fraction(0)
               for j in range(k)] for i in range(k)]
        t = Fraction(int(np.ceil((t_f + max(0.0, -lam) + 1e-9) * den)), den)
        for attempt in range(12):
            M = [[(t if i == j else Fraction(0)) + Aq[i][j] - 1 for j in range(k)] for i in range(k)]
            ok = is_psd_exact(M)
            print(f"  certify attempt {attempt}: t = {float(t):.6f}  psd = {ok}", flush=True)
            if ok:
                break
            t += Fraction(1, 1000)
        if ok:
            print(f"CERTIFIED: alpha <= floor({float(t):.6f}) = {int(t)}")
            rec["certified_t"] = str(t)
            rec["certified_alpha_ub"] = int(t)
        else:
            print("certification FAILED; float value stands as guidance only")
    if a.out:
        with open(a.out, "w") as fh:
            json.dump(rec, fh, indent=2)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
