"""Exact maximum-weight code in an orbit class whose conflict graph is nearly
complete: solve max-weight CLIQUE on the sparse compatibility graph written by
search/complement_edges.c (weights = orbit sizes), expand to words, verify.

    python search/compat_clique.py COMPAT.txt ORBITS.json n d [--out artifact]
"""
from __future__ import annotations

import hashlib
import itertools
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from verify_editq import check as verify_editq  # noqa: E402
from verify_editq2 import validate as verify_editq2  # noqa: E402


def main() -> int:
    import networkx as nx
    compat, orbits_path, n, d = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4])
    out = sys.argv[sys.argv.index("--out") + 1] if "--out" in sys.argv else None
    with open(compat) as fh:
        k, m = map(int, fh.readline().split())
        w = list(map(int, fh.readline().split()))
        G = nx.Graph()
        G.add_nodes_from(range(k))
        for line in fh:
            a, b = map(int, line.split())
            G.add_edge(a, b)
    for v in range(k):
        G.nodes[v]["weight"] = w[v]
    print(f"{k} orbits, {m} compatible pairs; max weight clique (exact, networkx)...", flush=True)
    # isolated vertices are cliques of size 1 (a single orbit); include them
    clique, weight = nx.max_weight_clique(G, weight="weight")
    best_single = max(w)
    if best_single > weight:
        clique, weight = [max(range(k), key=lambda v: w[v])], best_single
    orbits = json.load(open(orbits_path))
    words = ["".join("0123"[(idx // 4 ** i) % 4] for i in reversed(range(n))) for o in clique for idx in orbits[o]]
    words = sorted(words)
    ok1, size, why = verify_editq(words, n, d)
    ok2, why2 = verify_editq2(words, n, d, 4)
    print(f"clique of {len(clique)} orbits, weight {weight} -> {size} words; verify_editq={ok1} {why} verify_editq2={ok2} {why2}")
    if out and ok1 and ok2:
        body = "".join(x + "\n" for x in words)
        h = hashlib.sha256(body.encode()).hexdigest()
        with open(out, "w") as fh:
            fh.write(f"# quaternary edit-metric code, minimum Levenshtein distance {d}\n# n = {n}\n# d = {d}\n# q = 4\n# size = {size}\n")
            fh.write("# EXACT optimum within its symmetry class (max clique of the compatibility graph, networkx)\n")
            fh.write(f"# sha256(codewords) = {h}\n{body}")
        print("saved", out, h)
    return 0 if (ok1 and ok2) else 1


if __name__ == "__main__":
    raise SystemExit(main())
