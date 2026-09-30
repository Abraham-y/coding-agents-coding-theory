"""Export the BINARY s-deletion conflict graph G(n,s) to METIS for KaMIS.

Motivation: the published best-known codes at (2,14), (2,15), (2,16) are held
by OnlineMIS — measured on this machine to be strictly weaker than ReduMIS on
this graph family — and the whole table's incumbents come from a solver that
converges in minutes, not the 24h the timeout suggests. Multi-seed ReduMIS on
these cells is the best remaining shot at a new record.

Round trip:
    python search/binary_metis.py --n 14 --s 2 --export
    <redumis graphs/del_n14_s2.metis ...>
    python search/binary_metis.py --n 14 --s 2 --check result.mis
"""

from __future__ import annotations

import argparse
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from search.graph import load_graph  # noqa: E402
from verify import canonical_form, checksum, verify  # noqa: E402
from verify2 import verify2  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
GRAPH_DIR = os.path.join(ROOT, "graphs")

# Best known constructions, Weindel-Heckel arXiv:2504.00613v2 Table 2 (checked
# 2026-08-06). Engine attribution matters: OnlineMIS-held cells are softest.
PUBLISHED_BEST = {
    (12, 2): 37, (13, 2): 56, (14, 2): 87, (15, 2): 136, (16, 2): 215,
    (12, 3): 12, (13, 3): 16, (14, 3): 22,
}


def export(n: int, s: int):
    os.makedirs(GRAPH_DIR, exist_ok=True)
    t0 = time.time()
    g = load_graph(n, s)
    print(f"G({n},{s}): {g.num_vertices} vertices ({time.time()-t0:.0f}s to build)")
    path = os.path.join(GRAPH_DIR, f"del_n{n}_s{s}.metis")
    edges = 0
    with open(path, "w") as fh:
        lines = []
        for v in range(g.num_vertices):
            m = g.adj[v]
            nbrs = []
            while m:
                low = m & -m
                u = low.bit_length() - 1
                m ^= low
                nbrs.append(u + 1)  # METIS is 1-indexed; bit order is ascending
                if u > v:
                    edges += 1
            lines.append(" ".join(map(str, nbrs)))
        fh.write(f"{g.num_vertices} {edges}\n")
        fh.write("\n".join(lines) + "\n")
    print(f"wrote {path}: {edges} edges ({time.time()-t0:.0f}s)")
    return path


def check(n: int, s: int, result_path: str):
    g = load_graph(n, s)
    with open(result_path) as fh:
        flags = [line.strip() for line in fh if line.strip() in ("0", "1")]
    chosen = [i for i, f in enumerate(flags) if f == "1"]
    words = sorted(g.words[i] for i in chosen)
    ok1, sc1, r1 = verify(words, n, s)
    ok2, sc2, r2 = verify2(words, n, s) if len(words) <= 700 else (ok1, sc1, "(skipped, large)")
    pub = PUBLISHED_BEST.get((n, s))
    print(f"selected {len(chosen)}; verify={ok1} ({sc1}) verify2={ok2} ({sc2})")
    print(f"published best known: {pub}")
    if ok1 and ok2 and pub and sc1 > pub:
        out = os.path.join(ROOT, "best", f"n{n}_s{s}_size{sc1}_kamis.txt")
        os.makedirs(os.path.dirname(out), exist_ok=True)
        with open(out, "w") as fh:
            fh.write(canonical_form(words, n, s, "found by ReduMIS (KaMIS), local build"))
        print(f"checksum: {checksum(words, n, s)}")
        print(f"CANDIDATE RECORD saved: {out} -- novelty check required before any claim")
    return 0 if ok1 else 1


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, required=True)
    ap.add_argument("--s", type=int, required=True)
    ap.add_argument("--export", action="store_true")
    ap.add_argument("--check", type=str, default=None)
    args = ap.parse_args()
    if args.export:
        export(args.n, args.s)
        return 0
    if args.check:
        return check(args.n, args.s, args.check)
    ap.error("pass --export or --check FILE")


if __name__ == "__main__":
    raise SystemExit(main())
