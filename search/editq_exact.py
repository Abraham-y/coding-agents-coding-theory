"""Exact maximum-weight independent set on a cached orbit graph, by MaxSAT
(RC2, python-sat).  Used to CLOSE a symmetry class: if the optimum of the
H-invariant problem is W, no H-invariant code is larger than W.  This says
nothing about codes without that symmetry.

The verdict is only trusted when the returned set is expanded to codewords and
re-verified; the optimum value itself is the solver's word, and is recorded as
"class optimum (RC2)" rather than as a bound on the cell.

Instance = the mis_ils text file (k m / weights / edges).  Hard clauses:
(-a v -b) per edge.  Soft clauses: (v) with weight w_v.  RC2 minimises the
weight of falsified soft clauses, i.e. maximises the selected weight.

    python search/editq_exact.py graph.txt [--time 600] [--upper U]
"""
from __future__ import annotations

import argparse
import json
import sys
import time


def read_graph(path):
    with open(path) as fh:
        k, m = map(int, fh.readline().split())
        w = list(map(int, fh.readline().split()))
        E = []
        for line in fh:
            a, b = map(int, line.split())
            if a != b:
                E.append((a, b))
    return k, w, E


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("graph")
    ap.add_argument("--out", default=None)
    ap.add_argument("--solver", default="g4")
    a = ap.parse_args()
    from pysat.formula import WCNF
    from pysat.examples.rc2 import RC2

    k, w, E = read_graph(a.graph)
    wcnf = WCNF()
    for x, y in E:
        wcnf.append([-(x + 1), -(y + 1)])
    for v in range(k):
        wcnf.append([v + 1], weight=w[v])
    t0 = time.time()
    with RC2(wcnf, solver=a.solver, adapt=True, exhaust=True, minz=True) as rc2:
        model = rc2.compute()
        cost = rc2.cost
    dt = time.time() - t0
    sel = [v for v in range(k) if model[v] > 0]
    total = sum(w[v] for v in sel)
    # independent-set sanity check on the returned model
    adj = set(E) | {(y, x) for x, y in E}
    bad = sum(1 for i in sel for j in sel if i < j and (i, j) in adj)
    print(f"graph {a.graph}: {k} vertices, {len(E)} edges, total weight {sum(w)}")
    print(f"RC2 optimum: selected weight {total} ({len(sel)} vertices), cost {cost}, "
          f"{dt:.1f}s, violated edges in model: {bad}")
    print("SELECTED :", " ".join(map(str, sel)))
    if a.out:
        json.dump({"graph": a.graph, "k": k, "m": len(E), "optimum_weight": total,
                   "selected": sel, "seconds": round(dt, 1), "violations": bad},
                  open(a.out, "w"), indent=2)
    return 0 if bad == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
