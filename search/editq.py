"""Prescribed-symmetry search for fixed-length quaternary edit-metric codes
E_4(n, d) at general minimum distance d, using the same 33 conjugacy classes of
S_4 x C_2 as the d = 3 campaign and the C engine search/mis_ils.c.

WHY. The n=6, d=3 cell went 114 -> 120 by searching orbit graphs of subgroups of
the alphabet-permutation x reversal group. Nothing in that method is specific
to d = 3: symbol permutations and reversal preserve Levenshtein distance, so
they are automorphisms of every conflict graph "lev(u, v) <= d - 1". Houghten's
table has open cells at d = 4, 5, 6, 7, 8 with much smaller answers (e.g.
E_4(6,4) in [28, 32], E_4(7,5) in [18, 24]), where orbit graphs are tiny.

ENCODING GATE (rule 4, in the only form available: the table's lists were never
published). The cells the table marks as EXACT must be reproduced -- never
exceeded -- before any open cell is believed: E_4(5,3) = 44, E_4(5,4) = 11,
E_4(6,5) = 8, E_4(6,6) = 4, E_4(7,6) = 8, E_4(8,7) = 5. A search that beats an
exact cell has the wrong graph.

Records below are Houghten's table as fetched 2026-08-09 (notes/edit-metric-
lists-hunt/), lower bound and upper bound; (6,3) is this repo's 120.

    cc -O2 -o search/lev_graph search/lev_graph.c
    python search/editq.py --n 6 --d 4 --list
    python search/editq.py --n 6 --d 4 --classes all --seeds 4 --seconds 30
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import subprocess
import sys
import time
from collections import Counter
from concurrent.futures import ThreadPoolExecutor

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from search.edit4 import ALPHABET, orbit_problem, orbits_of  # noqa: E402
from search.edit4_diag import attainable_sizes, conjugacy_classes, transforms_of  # noqa: E402
from search.edit4_cmis import ensure_engine, run_engine, write_graph  # noqa: E402
from verify_editq import check as verify_editq  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS_DIR = os.environ.get("FUNMATH_RESULTS_DIR", os.path.join(ROOT, "results"))
BEST_DIR = os.path.join(ROOT, "best")
GRAPH_DIR = os.path.join(ROOT, "graphs")
LEVBIN = os.path.join(ROOT, "search", "lev_graph")
SCRATCH = os.environ.get("FUNMATH_SCRATCH", os.path.join(RESULTS_DIR, "cmis_tmp"))

# (n, d): (lower, upper) from Houghten's table; equal means the cell is exact.
TABLE = {
    (5, 3): (44, 44), (5, 4): (11, 11), (5, 5): (4, 4),
    (6, 3): (120, 176), (6, 4): (28, 32), (6, 5): (8, 8), (6, 6): (4, 4),
    (7, 3): (356, 614), (7, 4): (65, 128), (7, 5): (18, 24), (7, 6): (8, 8), (7, 7): (4, 4),
    (8, 3): (1132, 2340), (8, 4): (176, 512), (8, 5): (38, 96), (8, 6): (14, 32),
    (8, 7): (5, 5), (8, 8): (4, 4),
    (9, 3): (3451, 9360), (9, 4): (495, 2048), (9, 5): (97, 384), (9, 6): (27, 128),
    (9, 7): (11, 20), (9, 8): (5, 5), (9, 9): (4, 4),
    (10, 4): (1463, 8192), (10, 5): (249, 1536), (10, 6): (60, 496), (10, 7): (20, 80),
    (10, 8): (10, 16), (10, 9): (5, 5), (10, 10): (4, 4),
}


class LevGraph:
    """Same interface as Edit4Graph (words, num_vertices, adj bitmasks, n) so the
    orbit machinery of search/edit4.py applies verbatim."""

    def __init__(self, n: int, d: int, q: int = 4):
        self.n, self.d, self.q = n, d, q
        self.num_vertices = q ** n
        import itertools
        self.words = ["".join(w) for w in itertools.product(ALPHABET[:q], repeat=n)]
        os.makedirs(GRAPH_DIR, exist_ok=True)
        path = os.path.join(GRAPH_DIR, f"lev_q{q}_n{n}_d{d}.txt")
        if not os.path.exists(path):
            if not os.path.exists(LEVBIN) or os.path.getmtime(LEVBIN) < os.path.getmtime(LEVBIN + ".c"):
                subprocess.run(["cc", "-O2", "-o", LEVBIN, LEVBIN + ".c"], check=True)
            with open(path, "w") as fh:
                subprocess.run([LEVBIN, str(q), str(n), str(d)], stdout=fh, check=True)
        self.path = path
        adj = [0] * self.num_vertices
        with open(path) as fh:
            k, m = map(int, fh.readline().split())
            assert k == self.num_vertices
            fh.readline()
            for line in fh:
                a, b = map(int, line.split())
                adj[a] |= 1 << b
                adj[b] |= 1 << a
        self.adj = adj
        self.num_edges = m


def cell_instance(g, H):
    orbits = orbits_of(g, transforms_of(H))
    usable, edges = orbit_problem(g, orbits)
    k = len(usable)
    adj = [[] for _ in range(k)]
    for a, b in edges:
        adj[a].append(b)
        adj[b].append(a)
    return usable, adj, [len(o) for o in usable], orbits


def save_artifact(n, d, words, size, label, how, ci=None):
    body = "".join(w + "\n" for w in sorted(words))
    digest = hashlib.sha256(body.encode()).hexdigest()
    os.makedirs(BEST_DIR, exist_ok=True)
    tag = f"_class{ci}" if ci is not None else ""
    path = os.path.join(BEST_DIR, f"editq_n{n}_d{d}_size{size}{tag}.txt")
    with open(path, "w") as fh:
        fh.write(f"# quaternary edit-metric code, minimum Levenshtein distance d\n")
        fh.write(f"# n = {n}\n# d = {d}\n# q = 4\n# size = {size}\n")
        fh.write(f"# invariant under {label}\n# found by {how}\n")
        fh.write(f"# sha256(codewords) = {digest}\n")
        fh.write(body)
    return path, digest


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, required=True)
    ap.add_argument("--d", type=int, required=True)
    ap.add_argument("--list", action="store_true")
    ap.add_argument("--classes", type=str, default="", help="indices, or 'all'")
    ap.add_argument("--min-usable", type=int, default=1)
    ap.add_argument("--max-usable", type=int, default=10 ** 9)
    ap.add_argument("--seeds", type=int, default=4)
    ap.add_argument("--seconds", type=float, default=30.0)
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--tag", type=str, default="")
    a = ap.parse_args()

    ensure_engine()
    os.makedirs(SCRATCH, exist_ok=True)
    t0 = time.time()
    g = LevGraph(a.n, a.d)
    lb, ub = TABLE.get((a.n, a.d), (0, None))
    print(f"E_4({a.n},{a.d}): {g.num_vertices} words, {g.num_edges} conflict edges "
          f"(lev <= {a.d-1}); table [{lb}, {ub}]  ({time.time()-t0:.0f}s)", flush=True)
    classes = conjugacy_classes()
    classes.sort(key=lambda c: (-c["order"], c["diagonal"], c["name"]))

    rows = []
    for i, c in enumerate(classes):
        usable, adj, weights, orbits = cell_instance(g, c["H"])
        sizes = Counter(weights)
        reach = attainable_sizes(sizes, lb + 8) if lb else 0
        nxt = next((t for t in range(lb + 1, lb + 9) if (reach >> t) & 1), None)
        rows.append((i, c, usable, adj, weights, sizes, nxt))
        if a.list:
            print(f"{i:>2} |H|={c['order']:>2} {'DIAG' if c['diagonal'] else 'prod'} "
                  f"usable {len(usable):>6} edges {sum(map(len, adj))//2:>8} "
                  f"sizes {dict(sizes)} next>{lb}: {nxt}  {c['name']}")
    if a.list:
        return 0

    if a.classes == "all":
        want = [r[0] for r in rows if a.min_usable <= len(r[2]) <= a.max_usable]
    else:
        want = [int(x) for x in a.classes.split(",") if x.strip()]
    stamp = time.strftime("%Y%m%dT%H%M%S")
    log = {"n": a.n, "d": a.d, "table_lb": lb, "table_ub": ub, "engine": "search/mis_ils.c",
           "seeds": a.seeds, "seconds": a.seconds, "cells": []}
    best_overall = 0
    for ci in want:
        i, c, usable, adj, weights, sizes, nxt = rows[ci]
        if not usable:
            continue
        tb = time.time()
        gpath = os.path.join(SCRATCH, f"editq_n{a.n}_d{a.d}_class{ci}.txt")
        m = write_graph(gpath, weights, adj)
        secs = max(2.0, min(a.seconds, len(usable) / 20.0))
        seeds = list(range(1, a.seeds + 1))
        with ThreadPoolExecutor(max_workers=a.workers) as ex:
            res = list(ex.map(lambda s: run_engine(gpath, s, secs), seeds))
        by_seed = [r[0] for r in res]
        bi = max(range(len(res)), key=lambda j: res[j][0])
        wt, sel, iters = res[bi]
        words = sorted(g.words[v] for j in sel for v in usable[j])
        ok, size, reason = verify_editq(words, a.n, a.d)
        tag = "DIAG" if c["diagonal"] else "prod"
        flag = ""
        if ok and ub is not None and size > ub:
            flag = "  !!! EXCEEDS TABLE UPPER BOUND -- ENCODING BUG"
        elif ok and size > lb:
            flag = f"  *** beats table lower bound {lb}"
        print(f"  {ci:>2} {tag} |H|={c['order']:>2} usable {len(usable):>6} edges {m:>8} "
              f"{secs:>4.0f}s x{a.seeds} by seed {by_seed} best {wt:>5} "
              f"verified={ok} {reason}{flag}", flush=True)
        rec = {"class_index": ci, "name": c["name"], "order": c["order"],
               "diagonal": c["diagonal"], "orbits_usable": len(usable), "orbit_edges": m,
               "orbit_sizes": dict(sizes), "by_seed": by_seed, "best": wt,
               "verified": ok, "size": size if ok else 0,
               "next_attainable_above_lb": nxt, "seconds_per_seed": secs,
               "wall_seconds": round(time.time() - tb, 1)}
        if ok and size > lb and (ub is None or size <= ub):
            path, digest = save_artifact(a.n, a.d, words, size, c["name"],
                                         "search/editq.py (C ARW local search on the orbit graph)", ci)
            rec["artifact"], rec["sha256"] = path, digest
            print(f"      saved {path} sha256={digest} -- rule 6: bug until re-verified fresh",
                  flush=True)
        best_overall = max(best_overall, size if ok else 0)
        log["cells"].append(rec)
        log["elapsed_seconds"] = round(time.time() - t0, 1)
        os.makedirs(RESULTS_DIR, exist_ok=True)
        lp = os.path.join(RESULTS_DIR, f"editq_n{a.n}_d{a.d}{('_' + a.tag) if a.tag else ''}_{stamp}.json")
        with open(lp, "w") as fh:
            json.dump(log, fh, indent=2)
    print(f"\nbest overall {best_overall}; table [{lb}, {ub}]; log {lp}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
