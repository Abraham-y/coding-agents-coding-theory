"""Encoding validation for the quaternary edit-metric graph.

Produces two things the paper claims and the repo previously did not record:

1. an EXHAUSTIVE differential test of the confusability graph against the
   definition (all C(4096,2) = 8,386,560 vertex pairs, not a sample), and
2. exact values for the small cells where the published table gives a single
   number, which is the closest available substitute for ingesting an incumbent
   whose codeword list was never released.

Run:  python search/encoding_check.py [--exhaustive] [--small-cells]
"""
from __future__ import annotations
import argparse, json, os, sys, time
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from search.edit4 import Edit4Graph

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PUBLISHED_EXACT = {(3, 3): 4, (4, 3): 16}   # Houghten's table, single-valued cells


def lev(a: str, b: str) -> int:
    prev = list(range(len(b) + 1))
    for i, ca in enumerate(a, 1):
        cur = [i]
        for j, cb in enumerate(b, 1):
            cur.append(min(prev[j] + 1, cur[j - 1] + 1, prev[j - 1] + (ca != cb)))
        prev = cur
    return prev[-1]


def exhaustive_check(n: int = 6):
    g = Edit4Graph(n)
    # assert the vertex set itself, not just adjacency over it
    assert len(set(g.words)) == 4 ** n, "vertex set is not the 4^n distinct words"
    assert all(len(w) == n and set(w) <= set("0123") for w in g.words)
    N, mism, pairs = g.num_vertices, 0, 0
    for i in range(N):
        wi, mi = g.words[i], g.adj[i]
        for j in range(i + 1, N):
            pairs += 1
            if bool((mi >> j) & 1) != (lev(wi, g.words[j]) <= 2):
                mism += 1
    edges = sum(bin(a).count("1") for a in g.adj) // 2
    return {"n": n, "vertices": N, "edges": edges, "pairs_compared": pairs,
            "mismatches": mism}


def exact_small_cell(n: int):
    import networkx as nx
    g = Edit4Graph(n)
    G = nx.Graph(); G.add_nodes_from(range(g.num_vertices))
    for i in range(g.num_vertices):
        m = g.adj[i]
        for j in range(i + 1, g.num_vertices):
            if (m >> j) & 1:
                G.add_edge(i, j)
    clique, _ = nx.max_weight_clique(nx.complement(G), weight=None)
    return len(clique)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--exhaustive", action="store_true")
    ap.add_argument("--small-cells", action="store_true")
    a = ap.parse_args()
    if not (a.exhaustive or a.small_cells):
        a.exhaustive = a.small_cells = True
    out, t0 = {}, time.time()
    if a.small_cells:
        rows = []
        for (n, d), pub in sorted(PUBLISHED_EXACT.items()):
            got = exact_small_cell(n)
            ok = got == pub
            print(f"  E_4({n},{d}) exact = {got:>3}   published = {pub:>3}   "
                  f"{'MATCH' if ok else 'MISMATCH'}")
            rows.append({"n": n, "d": d, "computed": got, "published": pub, "match": ok})
        out["small_cells"] = rows
    if a.exhaustive:
        r = exhaustive_check(6)
        print(f"  exhaustive n=6: {r['pairs_compared']:,} pairs, "
              f"{r['edges']:,} edges, mismatches = {r['mismatches']}")
        out["exhaustive"] = r
    out["elapsed_seconds"] = round(time.time() - t0, 1)
    p = os.path.join(ROOT, "results", "encoding_check.json")
    json.dump(out, open(p, "w"), indent=2)
    print(f"  wrote {p}")
