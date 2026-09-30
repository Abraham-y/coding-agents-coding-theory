"""Lovasz theta (SDP) upper bounds on A_s(n), with exact rational certification.

Motivation, direct from reading OpenAI's Astra formalisation (2026-08-01). Their
binary-codes result chains

    binaryRate  <=  combinedVariationalRate  <  mrrwRate

i.e. it beats MRRW, which IS the Delsarte linear programming bound, by moving to
Gram-matrix / positive-semidefinite arguments (`projectionGrams`, Jacobi
matrices, Johnson spheres). Their setting is the Hamming metric in the asymptotic
rate regime, so it does not touch the deletion cells here — but the STRUCTURAL
lesson transfers: where an LP bound is loose, the semidefinite strengthening is
the natural next rung.

This repo's existing bound (search/lp_bound.py) is exactly an LP: a fractional
relaxation over the clique cover. The analogue of "go semidefinite" for an
independent-set problem is the Lovasz theta function, which satisfies

    alpha(G)  <=  theta(G)  <=  chi_f(complement)  <=  LP clique-cover bound

so theta is never worse than what we already have, and is often strictly better.

Formulation used (Lovasz's original):
    theta(G) = max  sum_ij X_ij
               s.t. trace(X) = 1
                    X_ij = 0   whenever i ~ j in G
                    X psd

CERTIFICATION. The SDP solver is floating point and is NEVER trusted for a
verdict (CLAUDE.md rule 10). It only proposes. The reported number is obtained
from the DUAL: any symmetric Y with Y_ij <= 0 off-diagonal on non-edges... more
precisely we use the standard dual

    theta(G) = min  t
               s.t. t*I + A - J  psd,  A supported on edges of G

Any feasible (t, A) proves theta(G) <= t, hence alpha(G) <= t. We take the
solver's proposal, round it to exact rationals, repair by inflating t, and then
verify positive semidefiniteness EXACTLY via rational LDL^T. The reported bound
is floor(t) as an exact rational (rule 12).

Run:
    python search/theta.py --n 8 --s 2
    python search/theta.py --n 10 --s 2 --no-certify   # float only, fast
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

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from search.graph import load_graph  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS_DIR = os.path.join(ROOT, "results")


def edges_of(g):
    out = []
    for v in range(g.num_vertices):
        m = g.adj[v]
        while m:
            low = m & -m
            u = low.bit_length() - 1
            m ^= low
            if u > v:
                out.append((v, u))
    return out


def solve_theta(g, solver=None, verbose=False):
    """Float SDP. Returns (theta_value, t, A_dict) from the dual. Never trusted."""
    import cvxpy as cp

    n = g.num_vertices
    E = edges_of(g)
    t = cp.Variable()
    # A is supported on edges only
    a = cp.Variable(len(E))
    M = t * np.eye(n) - np.ones((n, n))
    # build t*I + A - J  with A_ij = A_ji = a_k on edge k
    rows = [e[0] for e in E] + [e[1] for e in E]
    cols = [e[1] for e in E] + [e[0] for e in E]
    A = cp.Constant(np.zeros((n, n)))
    if E:
        Amat = cp.reshape(
            cp.sum([a[k] * _unit(n, E[k]) for k in range(len(E))]), (n, n), order="C"
        ) if False else None
    # the list-comprehension form above is O(|E|) expressions and is far too slow;
    # build a sparse selector instead
    import scipy.sparse as sp

    P = sp.csr_matrix(
        (np.ones(2 * len(E)), (np.array(rows) * n + np.array(cols), np.arange(2 * len(E)) % len(E))),
        shape=(n * n, len(E)),
    )
    Avec = P @ a
    Amat = cp.reshape(Avec, (n, n), order="C")
    constraints = [t * np.eye(n) + Amat - np.ones((n, n)) >> 0]
    prob = cp.Problem(cp.Minimize(t), constraints)
    prob.solve(solver=solver, verbose=verbose)
    if prob.status not in ("optimal", "optimal_inaccurate"):
        raise RuntimeError(f"SDP status {prob.status}")
    return float(t.value), float(t.value), (E, np.asarray(a.value).ravel())


def _unit(n, e):
    Z = np.zeros((n, n))
    Z[e[0], e[1]] = 1.0
    Z[e[1], e[0]] = 1.0
    return Z


# ---------------------------------------------------------------------------
# exact certification: rational LDL^T psd check
# ---------------------------------------------------------------------------


def is_psd_exact(M: list[list[Fraction]]) -> bool:
    """Exact PSD test by symmetric Gaussian elimination (rational LDL^T).

    Works for singular matrices: a zero pivot is acceptable only if the whole
    remaining row is zero, otherwise the matrix is indefinite.
    """
    n = len(M)
    A = [row[:] for row in M]
    for k in range(n):
        if A[k][k] < 0:
            return False
        if A[k][k] == 0:
            for j in range(k + 1, n):
                if A[k][j] != 0:
                    return False
            continue
        piv = A[k][k]
        for i in range(k + 1, n):
            if A[i][k] == 0:
                continue
            f = A[i][k] / piv
            for j in range(k, n):
                A[i][j] -= f * A[k][j]
    return True


def certify(g, t_float: float, E, a_float, denom: int = 10 ** 6):
    """Round the dual to rationals, inflate t until exactly psd, return the bound."""
    n = g.num_vertices
    a_r = [Fraction(int(round(x * denom)), denom) for x in a_float]

    def build(t: Fraction):
        M = [[Fraction(-1) for _ in range(n)] for _ in range(n)]
        for i in range(n):
            M[i][i] = t - 1
        for k, (i, j) in enumerate(E):
            M[i][j] += a_r[k]
            M[j][i] += a_r[k]
        return M

    t = Fraction(int(math.ceil(t_float * denom)), denom)
    for _ in range(60):
        if is_psd_exact(build(t)):
            return t
        t += Fraction(1, 100)
    raise RuntimeError("could not repair the dual to an exactly psd certificate")


def verify_certificate(path: str):
    """Re-check a stored theta certificate from disk, exact arithmetic only.

    Rebuilds t*I + A - J from the stored (t, dual_a, dual_edges), re-runs the
    exact rational LDL^T PSD check, and confirms floor(t) matches the claimed
    bound. Requires the record to contain dual_a/dual_edges (records written
    before 2026-08-07 predate the schema and are not re-checkable)."""
    import json as _json

    with open(path) as fh:
        cert = _json.load(fh)
    if "dual_a" not in cert:
        return False, "record predates the re-checkable schema (no dual stored)"
    n, s = cert["n"], cert["s"]
    g = load_graph(n, s)
    E = [tuple(e) for e in cert["dual_edges"]]
    a_r = [Fraction(p, q) for p, q in cert["dual_a"]]
    if len(a_r) != len(E):
        return False, f"malformed certificate: {len(a_r)} dual values for {len(E)} edges"
    t = Fraction(*cert["certified_t"])
    nv = g.num_vertices
    M = [[Fraction(-1)] * nv for _ in range(nv)]
    for i in range(nv):
        M[i][i] = t - 1
    for k, (i, j) in enumerate(E):
        if not (g.adj[i] >> j) & 1:
            return False, f"stored dual edge ({i},{j}) is not an edge of G({n},{s})"
        M[i][j] += a_r[k]
        M[j][i] += a_r[k]
    if not is_psd_exact(M):
        return False, "stored dual matrix is not PSD"
    if math.floor(t) != cert["certified_upper_bound"]:
        return False, "floor(t) does not match the claimed bound"
    return True, f"OK: A_{s}({n}) <= {cert['certified_upper_bound']} (exact t = {t})"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int)
    ap.add_argument("--s", type=int)
    ap.add_argument("--solver", type=str, default=None)
    ap.add_argument("--no-certify", action="store_true")
    ap.add_argument("--denom", type=int, default=10 ** 6)
    ap.add_argument("--verify-cert", type=str, default=None)
    args = ap.parse_args()

    if args.verify_cert:
        ok, msg = verify_certificate(args.verify_cert)
        print(("PASS " if ok else "FAIL ") + msg)
        return 0 if ok else 1
    if args.n is None or args.s is None:
        ap.error("--n and --s are required unless --verify-cert is given")

    t0 = time.time()
    g = load_graph(args.n, args.s)
    E = edges_of(g)
    print(f"G({args.n},{args.s}): {g.num_vertices} vertices, {len(E)} edges")

    theta, t_float, (E2, a) = solve_theta(g, solver=args.solver)
    # DO NOT print a floor of this. The float is routinely a hair BELOW the true
    # theta -- e.g. G(6,2) returns 3.999997 where theta = 4 -- so flooring it
    # yields 3 and claims A_2(6) <= 3, which is false. This is exactly the
    # failure CLAUDE.md rule 10 forbids. Only the exactly certified dual below
    # may be turned into an integer bound.
    print(f"float theta (NOT a result, NOT floorable): {theta:.6f}")

    record = {
        "method": "Lovasz theta SDP, dual certified in exact rational arithmetic",
        "n": args.n, "s": args.s,
        "vertices": g.num_vertices, "edges": len(E),
        "theta_float_untrusted": theta,
    }

    if not args.no_certify:
        t_exact = certify(g, t_float, E2, a, args.denom)
        ub = math.floor(t_exact)
        print(f"certified exact dual value : {t_exact} (= {float(t_exact):.6f})")
        print(f"CERTIFIED: A_{args.s}({args.n}) <= {ub}")
        # Store the FULL dual so the certificate is re-checkable from disk
        # (audit 2026-08-07 flagged that certified_t alone is not a certificate:
        # without the edge multipliers, is_psd_exact cannot be re-run). The
        # rounded rationals a_r are reconstructed exactly as certify() built
        # them, so (t, denom, a_r) determines the certified matrix.
        a_r = [Fraction(int(round(x * args.denom)), args.denom) for x in a]
        record.update({
            "certified_t": [t_exact.numerator, t_exact.denominator],
            "certified_upper_bound": ub,
            "denominator": args.denom,
            "dual_edges": E2,
            "dual_a": [[f.numerator, f.denominator] for f in a_r],
        })

    record["elapsed_seconds"] = round(time.time() - t0, 2)
    os.makedirs(RESULTS_DIR, exist_ok=True)
    path = os.path.join(RESULTS_DIR, f"theta_n{args.n}_s{args.s}.json")
    with open(path, "w") as fh:
        json.dump(record, fh, indent=2)
    print(f"elapsed {record['elapsed_seconds']}s  log {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
