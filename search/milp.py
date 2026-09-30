"""Exact / branch-and-bound integer programming for A_s(n), via HiGHS.

A rung the earlier work skipped. The set-packing formulation is the same clique
cover used for the LP bound in lp_bound.py, but with integrality enforced:

    max  sum_v x_v
    s.t. sum_{v : t is a subsequence of v} x_v <= 1   for every word t of length n-s
         x_v in {0,1}

For small n this closes the cell outright (optimal, with a proof of optimality
from the solver's own bound). For larger n it yields an incumbent plus a gap.

IMPORTANT: the solver's optimality claim is NOT taken as ground truth. Any
incumbent is re-verified by verify.py AND verify2.py before being reported, and
an optimal-status claim is only ever reported alongside the independently
certified LP bound from lp_bound.py.

Run:
    python search/milp.py --n 8 --s 2                 # exact
    python search/milp.py --n 12 --s 2 --seconds 600  # incumbent + gap
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time

import numpy as np
from scipy.optimize import linprog, milp, Bounds, LinearConstraint
from scipy.sparse import csr_matrix

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from search.graph import load_graph  # noqa: E402
from verify import canonical_form, checksum, verify  # noqa: E402
from verify2 import verify2  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS_DIR = os.path.join(ROOT, "results")
BEST_DIR = os.path.join(ROOT, "best")


def build_matrix(n: int, s: int):
    g = load_graph(n, s)
    names = sorted(g.cliques)
    rows, cols = [], []
    for i, t in enumerate(names):
        for v in g.cliques[t]:
            rows.append(i)
            cols.append(v)
    A = csr_matrix(
        (np.ones(len(rows)), (rows, cols)), shape=(len(names), g.num_vertices)
    )
    return g, A


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, required=True)
    ap.add_argument("--s", type=int, required=True)
    ap.add_argument("--seconds", type=float, default=300)
    ap.add_argument("--tag", type=str, default="milp")
    args = ap.parse_args()

    t0 = time.time()
    g, A = build_matrix(args.n, args.s)
    nv = g.num_vertices
    print(f"MILP: {nv} binary vars, {A.shape[0]} clique constraints")

    res = milp(
        c=-np.ones(nv),
        constraints=LinearConstraint(A, -np.inf, 1),
        integrality=np.ones(nv),
        bounds=Bounds(0, 1),
        options={"time_limit": args.seconds, "presolve": True},
    )
    elapsed = time.time() - t0

    status = res.status
    msg = res.message
    print(f"status {status}: {msg}")
    print(f"elapsed {elapsed:.1f}s")

    record = {
        "method": "set-packing MILP via HiGHS (scipy.optimize.milp)",
        "n": args.n, "s": args.s,
        "vars": nv, "constraints": A.shape[0],
        "time_limit": args.seconds,
        "status": int(status), "message": str(msg),
        "elapsed_seconds": round(elapsed, 2),
    }

    if res.x is not None:
        chosen = [v for v in range(nv) if res.x[v] > 0.5]
        words = sorted(g.words[v] for v in chosen)
        ok1, sc1, r1 = verify(words, args.n, args.s)
        ok2, sc2, r2 = verify2(words, args.n, args.s) if len(words) <= 800 else (ok1, sc1, "")
        print(f"incumbent size {len(words)}")
        print(f"  verify.py : valid={ok1} score={sc1} {r1}")
        print(f"  verify2.py: valid={ok2} score={sc2} {r2}")
        # HiGHS reports the dual/best bound for the maximisation as -res.mip_dual_bound
        bound = getattr(res, "mip_dual_bound", None)
        if bound is not None:
            print(f"  solver best bound (NOT a certificate): {-bound:.4f}")
            record["solver_bound_uncertified"] = float(-bound)
        record.update({"size": sc1, "verify": ok1, "verify2": ok2, "words": words})
        if ok1 and ok2:
            os.makedirs(BEST_DIR, exist_ok=True)
            path = os.path.join(BEST_DIR, f"n{args.n}_s{args.s}_size{sc1}_milp.txt")
            with open(path, "w") as fh:
                fh.write(canonical_form(words, args.n, args.s,
                                        f"HiGHS MILP, status={status}"))
            print(f"  artifact  : {path}")
            print(f"  checksum  : {checksum(words, args.n, args.s)}")
            record["artifact"] = path

    os.makedirs(RESULTS_DIR, exist_ok=True)
    stamp = time.strftime("%Y%m%dT%H%M%S")
    logpath = os.path.join(RESULTS_DIR, f"{args.tag}_n{args.n}_s{args.s}_{stamp}.json")
    with open(logpath, "w") as fh:
        json.dump(record, fh, indent=2)
    print(f"  log       : {logpath}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
