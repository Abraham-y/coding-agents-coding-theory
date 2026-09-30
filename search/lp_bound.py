"""Certified upper bounds on A_s(n) via LP duality in exact rational arithmetic.

Why this is the tractable side of the problem for this sprint:

  A better LOWER bound is a code, and a code is checked instantly by verify.py.
  A better UPPER bound from SAT is a nonexistence claim needing a DRAT proof,
  which CLAUDE.md rule 11 puts out of scope. But an upper bound from LP DUALITY
  needs only a rational vector y that a few hundred thousand Fraction additions
  can validate. So the upper-bound side needs no certificate infrastructure at
  all -- it IS the certificate.

Primal (fractional relaxation of the code problem):
    max  sum_v x_v
    s.t. sum_{v : t is a subsequence of v} x_v <= 1   for every word t of length n-s
         x >= 0

Dual:
    min  sum_t y_t
    s.t. sum_{t : t is a subsequence of v} y_t >= 1   for every v in {0,1}^n
         y >= 0

ANY dual-feasible y proves A_s(n) <= sum_t y_t. So the float LP solver is used
only to PROPOSE a y; it is never trusted. The proposal is rounded to exact
rationals, repaired by exact scaling until every one of the 2^n constraints holds
with certainty, and only then does its objective become a reported bound. No
float is compared in any accept/reject decision (rule 10), and the reported
result is the certified value, not the number the optimiser printed (rule 12).

Run:
    python search/lp_bound.py --n 12 --s 2
    python search/lp_bound.py --n 12 --s 2 --cuts        # add maximal-clique cuts
    python search/lp_bound.py --verify-cert results/cert_n12_s2.json
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

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from search.graph import load_graph  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS_DIR = os.path.join(ROOT, "results")


# ---------------------------------------------------------------------------
# constraint system
# ---------------------------------------------------------------------------


def build_system(n: int, s: int, cuts: bool = False):
    """Return (constraints, ball_of_vertex, names).

    constraints[i] = sorted tuple of vertex indices forming clique i.
    ball_of_vertex[v] = list of constraint indices containing v.
    """
    g = load_graph(n, s)
    names = sorted(g.cliques)
    constraints = [tuple(sorted(g.cliques[t])) for t in names]

    if cuts:
        constraints, names = _add_clique_cuts(g, constraints, names)

    ball = [[] for _ in range(g.num_vertices)]
    for i, members in enumerate(constraints):
        for v in members:
            ball[v].append(i)
    return g, constraints, ball, names


def _add_clique_cuts(g, constraints, names):
    """Greedily extend each subsequence-clique to a maximal clique.

    Every clique K of G gives a valid inequality sum_{v in K} x_v <= 1, so any
    strictly larger clique is a strictly stronger constraint. Extending is
    cheap with bitset adjacency and can only tighten the LP.
    """
    seen = set(constraints)
    extra = []
    for members in constraints:
        common = None
        for v in members:
            common = g.adj[v] if common is None else (common & g.adj[v])
            if not common:
                break
        if not common:
            continue
        grown = list(members)
        while common:
            low = common & -common
            u = low.bit_length() - 1
            grown.append(u)
            common &= g.adj[u]
        key = tuple(sorted(grown))
        if key not in seen and len(key) > len(members):
            seen.add(key)
            extra.append(key)
    out = list(constraints) + extra
    out_names = list(names) + [f"cut{i}" for i in range(len(extra))]
    return out, out_names


# ---------------------------------------------------------------------------
# float LP: propose a dual vector (never trusted)
# ---------------------------------------------------------------------------


def propose_dual(num_vertices: int, constraints, ball):
    """Solve the DUAL directly as a float LP; return a float y proposal."""
    m = len(constraints)
    # min 1^T y  s.t.  -B y <= -1,  y >= 0   where B[v, t] = 1 iff t in ball(v)
    rows, cols = [], []
    for v, ids in enumerate(ball):
        for t in ids:
            rows.append(v)
            cols.append(t)
    data = np.full(len(rows), -1.0)
    B = csr_matrix((data, (rows, cols)), shape=(num_vertices, m))
    res = linprog(
        c=np.ones(m),
        A_ub=B,
        b_ub=np.full(num_vertices, -1.0),
        bounds=(0, None),
        method="highs",
    )
    if not res.success:
        raise RuntimeError(f"LP failed: {res.message}")
    return np.asarray(res.x, dtype=float), float(res.fun)


# ---------------------------------------------------------------------------
# exact certification
# ---------------------------------------------------------------------------


def certify(y_float, ball, denom: int):
    """Round to rationals, repair by exact scaling, verify all constraints exactly.

    Returns (y_exact, bound_fraction, min_slack). Raises if repair is impossible.
    """
    # round UP so most constraints survive rounding without repair
    y = [Fraction(math.ceil(max(0.0, val) * denom), denom) for val in y_float]

    # exact evaluation of every dual constraint
    def coverage():
        return [sum((y[t] for t in ids), Fraction(0)) for ids in ball]

    cov = coverage()
    if any(c == 0 for c in cov):
        raise RuntimeError("proposal leaves a vertex with zero coverage; cannot repair")

    worst = min(cov)
    if worst < 1:
        scale = 1 / worst  # exact Fraction
        y = [val * scale for val in y]
        cov = coverage()

    min_slack = min(cov)
    if min_slack < 1:
        raise RuntimeError(f"repair failed, min coverage {min_slack}")

    bound = sum(y, Fraction(0))
    return y, bound, min_slack


def verify_certificate(cert: dict) -> tuple[bool, str]:
    """Standalone re-check of a saved certificate. Pure exact arithmetic."""
    n, s = cert["n"], cert["s"]
    g, constraints, ball, _ = build_system(n, s, cuts=cert.get("cuts", False))
    if len(constraints) != len(cert["y"]):
        return False, "certificate length does not match rebuilt constraint system"
    y = [Fraction(a, b) for a, b in cert["y"]]
    for v, ids in enumerate(ball):
        total = sum((y[t] for t in ids), Fraction(0))
        if total < 1:
            return False, f"dual constraint violated at vertex {v}: {total} < 1"
    bound = sum(y, Fraction(0))
    claimed = Fraction(*cert["bound"])
    if bound != claimed:
        return False, f"objective mismatch: recomputed {bound}, claimed {claimed}"
    if math.floor(bound) != cert["certified_upper_bound"]:
        return False, "floor of objective does not match the claimed integer bound"
    return True, f"OK: A_{s}({n}) <= {cert['certified_upper_bound']} (exact dual objective {bound})"


# ---------------------------------------------------------------------------


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int)
    ap.add_argument("--s", type=int)
    ap.add_argument("--cuts", action="store_true")
    ap.add_argument("--denom", type=int, default=100000)
    ap.add_argument("--verify-cert", type=str, default=None)
    args = ap.parse_args()

    if args.verify_cert:
        with open(args.verify_cert) as fh:
            cert = json.load(fh)
        ok, msg = verify_certificate(cert)
        print(("PASS " if ok else "FAIL ") + msg)
        return 0 if ok else 1

    if args.n is None or args.s is None:
        ap.error("--n and --s are required unless --verify-cert is given")

    t0 = time.time()
    g, constraints, ball, names = build_system(args.n, args.s, args.cuts)
    print(f"system: {g.num_vertices} vertices, {len(constraints)} constraints"
          f"{' (with maximal-clique cuts)' if args.cuts else ''}")

    y_float, lp_val = propose_dual(g.num_vertices, constraints, ball)
    print(f"float LP objective (NOT a result): {lp_val:.6f}")

    y, bound, min_slack = certify(y_float, ball, args.denom)
    ub = math.floor(bound)
    elapsed = time.time() - t0

    print(f"certified exact dual objective : {bound}  (= {float(bound):.6f})")
    print(f"min dual coverage (must be >=1): {min_slack}")
    print(f"CERTIFIED: A_{args.s}({args.n}) <= {ub}")
    print(f"elapsed: {elapsed:.1f}s")

    os.makedirs(RESULTS_DIR, exist_ok=True)
    suffix = "_cuts" if args.cuts else ""
    path = os.path.join(RESULTS_DIR, f"cert_n{args.n}_s{args.s}{suffix}.json")
    cert = {
        "method": "LP duality, dual proposed by HiGHS then rounded and repaired exactly",
        "n": args.n,
        "s": args.s,
        "cuts": args.cuts,
        "denominator": args.denom,
        "float_lp_objective": lp_val,
        "bound": [bound.numerator, bound.denominator],
        "certified_upper_bound": ub,
        "y": [[val.numerator, val.denominator] for val in y],
    }
    with open(path, "w") as fh:
        json.dump(cert, fh)
    print(f"certificate: {path}")

    ok, msg = verify_certificate(cert)
    print(("re-check PASS " if ok else "re-check FAIL ") + msg)
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
