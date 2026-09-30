"""Cells this repo closed EXACTLY, by two independent methods that share no code.

Method 1: set-packing MILP solved to proven optimality by HiGHS.
Method 2: exact maximum clique on the complement graph, via networkx.

These are optimality proofs, not best-known constructions. A_2(8) = 7 here
CORRECTS a value this repo originally transcribed as 8; see
notes/external-data-findings.md.

Run:  python tests/test_exact_cells.py
"""
from __future__ import annotations
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np
import networkx as nx
from scipy.optimize import milp, Bounds, LinearConstraint
from scipy.sparse import csr_matrix

from search.graph import load_graph
from verify import verify
from verify2 import verify2
from vt import EXACTLY_SETTLED_HERE

FAIL = []

def milp_opt(n, s):
    g = load_graph(n, s)
    names = sorted(g.cliques)
    rows, cols = [], []
    for i, t in enumerate(names):
        for v in g.cliques[t]:
            rows.append(i); cols.append(v)
    A = csr_matrix((np.ones(len(rows)), (rows, cols)), shape=(len(names), g.num_vertices))
    r = milp(c=-np.ones(g.num_vertices), constraints=LinearConstraint(A, -np.inf, 1),
             integrality=np.ones(g.num_vertices), bounds=Bounds(0, 1),
             options={"time_limit": 300})
    assert r.status == 0, f"MILP did not prove optimality for ({s},{n}): {r.message}"
    chosen = [v for v in range(g.num_vertices) if r.x[v] > 0.5]
    return sorted(g.words[v] for v in chosen)

def clique_opt(n, s):
    g = load_graph(n, s)
    G = nx.Graph(); G.add_nodes_from(range(g.num_vertices))
    for v in range(g.num_vertices):
        m = g.adj[v]
        while m:
            low = m & -m; u = low.bit_length() - 1; m ^= low
            if u > v: G.add_edge(v, u)
    clique, _ = nx.max_weight_clique(nx.complement(G), weight=None)
    return sorted(g.words[v] for v in clique)

def main() -> int:
    checks = 0
    for (s, n), expected in sorted(EXACTLY_SETTLED_HERE.items()):
        w1 = milp_opt(n, s); w2 = clique_opt(n, s)
        for label, w in (("MILP", w1), ("max-clique", w2)):
            ok1, sc1, r1 = verify(w, n, s); ok2, sc2, _ = verify2(w, n, s)
            checks += 1
            if not (ok1 and ok2):
                FAIL.append(f"A_{s}({n}) {label}: artifact rejected ({r1})")
            elif sc1 != expected:
                FAIL.append(f"A_{s}({n}) {label}: size {sc1}, expected {expected}")
        if len(w1) != len(w2):
            FAIL.append(f"A_{s}({n}): MILP {len(w1)} != max-clique {len(w2)}")
        checks += 1
    print(f"checks run : {checks}")
    print(f"failures   : {len(FAIL)}")
    for f in FAIL: print("  FAIL " + f)
    return 1 if FAIL else 0

if __name__ == "__main__":
    raise SystemExit(main())
