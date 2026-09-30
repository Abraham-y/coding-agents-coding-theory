"""Lovasz theta for a DIMACS-format graph (Sloane challenge instances).

Same solver and certification machinery as search/theta.py, but the graph
comes from an edge-list file instead of a deletion-code construction. Built for
the 1zc.* single-asymmetric-error challenge graphs, whose published upper
bounds predate modern SDP entirely (1zc.1024: [112, 117] since 1998, marked
"Most wanted!" on the challenge page).

The graph file is trusted only because it was independently rebuilt from the
asymmetric-distance definition (z <= x with |x\\z| <= 1 reachable from both
words) and matched edge-for-edge at n = 7, 8, 9, 10 — see notes.

Run:
    python search/theta_dimacs.py --graph <path.txt> --no-certify     # float probe
    python search/theta_dimacs.py --graph <path.txt>                  # + exact rational certificate
    python search/theta_dimacs.py --verify-cert results/theta_1zc1024.json --graph <path.txt>
"""

from __future__ import annotations

import argparse
import json
import math
import os
import sys
import time
from fractions import Fraction

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from search.theta import certify, is_psd_exact, solve_theta  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS_DIR = os.path.join(ROOT, "results")


class DimacsGraph:
    def __init__(self, path: str):
        adj: dict[int, int] = {}
        nv = 0
        for line in open(path):
            t = line.split()
            if not t:
                continue
            if t[0] == "p":
                nv = int(t[2])
            elif t[0] == "e":
                a, b = int(t[1]) - 1, int(t[2]) - 1
                if a == b:
                    continue
                adj[a] = adj.get(a, 0) | (1 << b)
                adj[b] = adj.get(b, 0) | (1 << a)
        self.num_vertices = nv
        self.adj = [adj.get(v, 0) for v in range(nv)]
        self.words = [str(v) for v in range(nv)]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--graph", type=str, required=True)
    ap.add_argument("--no-certify", action="store_true")
    ap.add_argument("--denom", type=int, default=10 ** 6)
    ap.add_argument("--solver", type=str, default=None)
    ap.add_argument("--verify-cert", type=str, default=None)
    args = ap.parse_args()

    g = DimacsGraph(args.graph)
    name = os.path.basename(args.graph).replace(".txt", "").replace(".", "_")
    edges = sum(bin(a).count("1") for a in g.adj) // 2
    print(f"{name}: {g.num_vertices} vertices, {edges} edges")

    if args.verify_cert:
        with open(args.verify_cert) as fh:
            cert = json.load(fh)
        E = [tuple(e) for e in cert["dual_edges"]]
        a_r = [Fraction(p, q) for p, q in cert["dual_a"]]
        t = Fraction(*cert["certified_t"])
        nv = g.num_vertices
        M = [[Fraction(-1)] * nv for _ in range(nv)]
        for i in range(nv):
            M[i][i] = t - 1
        for k, (i, j) in enumerate(E):
            if not (g.adj[i] >> j) & 1:
                print(f"FAIL stored dual edge ({i},{j}) not in graph")
                return 1
            M[i][j] += a_r[k]
            M[j][i] += a_r[k]
        if not is_psd_exact(M):
            print("FAIL dual matrix not PSD")
            return 1
        if math.floor(t) != cert["certified_upper_bound"]:
            print("FAIL floor mismatch")
            return 1
        print(f"PASS alpha({name}) <= {cert['certified_upper_bound']} (exact t = {t})")
        return 0

    t0 = time.time()
    theta, t_float, (E2, a) = solve_theta(g, solver=args.solver)
    print(f"float theta (NOT a result, NOT floorable): {theta:.6f}")
    record = {
        "method": "Lovasz theta SDP on DIMACS graph, dual certified in exact rationals",
        "graph": os.path.basename(args.graph),
        "vertices": g.num_vertices, "edges": edges,
        "theta_float_untrusted": theta,
    }
    if not args.no_certify:
        t_exact = certify(g, t_float, E2, a, args.denom)
        ub = math.floor(t_exact)
        a_r = [Fraction(int(round(x * args.denom)), args.denom) for x in a]
        print(f"certified exact dual value : {t_exact} (= {float(t_exact):.6f})")
        print(f"CERTIFIED: alpha({name}) <= {ub}")
        record.update({
            "certified_t": [t_exact.numerator, t_exact.denominator],
            "certified_upper_bound": ub,
            "denominator": args.denom,
            "dual_edges": E2,
            "dual_a": [[f.numerator, f.denominator] for f in a_r],
        })
    record["elapsed_seconds"] = round(time.time() - t0, 2)
    os.makedirs(RESULTS_DIR, exist_ok=True)
    path = os.path.join(RESULTS_DIR, f"theta_{name}.json")
    with open(path, "w") as fh:
        json.dump(record, fh)
    print(f"elapsed {record['elapsed_seconds']}s  log {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
