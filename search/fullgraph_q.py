"""Full-conflict-graph MIS search for fixed-length q-ary edit-metric codes E_q(n, d),
q = 2 or 3, where the whole graph (q^n words) is small enough to hand to
search/mis_ils.c directly with no symmetry assumption.

WHY. The binary and ternary tables were last changed in 2005-2008 (and once in
2024 for binary), and every open cell with q^n <= 3^9 is a graph of at most
19,683 vertices: brute-force-sized for a modern ILS.  This is rung 1 of the
ladder for those tables.  The cells the tables mark EXACT are run too, as an
encoding gate: a result above an exact value means the graph is wrong.

Nothing printed here is a result.  Each best set is expanded to codewords,
written to results/fullgraph_q/ as an edit artifact, and must pass verify.py
in a fresh process before it goes anywhere else.

    cc -O2 -o search/trie_graph search/trie_graph.c   (same graph as lev_graph.c, checked edge for edge)
    python search/fullgraph_q.py --cells 3:8:3,3:9:3 --seeds 1 --seconds 120 --workers 4
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import subprocess
import sys
import time
from multiprocessing import Pool

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "results", "fullgraph_q")
LOG = os.path.join(ROOT, "results", "fullgraph_q_runs.jsonl")


def word(idx, q, n):
    s = []
    for _ in range(n):
        s.append(str(idx % q))
        idx //= q
    return "".join(reversed(s))


def lev(a, b):
    prev = list(range(len(b) + 1))
    for i, ca in enumerate(a, 1):
        cur = [i]
        for j, cb in enumerate(b, 1):
            cur.append(min(prev[j] + 1, cur[j - 1] + 1, prev[j - 1] + (ca != cb)))
        prev = cur
    return prev[-1]


def graph_path(q, n, d):
    path = os.path.join(OUT, f"graph_q{q}_n{n}_d{d}.txt")
    if not os.path.exists(path):
        tmp = path + ".tmp"
        with open(tmp, "w") as fh:
            subprocess.run([os.path.join(ROOT, "search", "trie_graph"), str(q), str(n), str(d)],
                           stdout=fh, check=True)
        os.replace(tmp, path)
    return path


def run(job):
    q, n, d, seed, secs = job
    g = graph_path(q, n, d)
    t0 = time.time()
    p = subprocess.run([os.path.join(ROOT, "search", "mis_ils"), g, str(seed), str(secs), "--quiet"],
                       capture_output=True, text=True, check=True)
    wall = time.time() - t0
    best = [l for l in p.stdout.splitlines() if l.startswith("BEST")][-1]
    verts = [int(x) for x in best.split(":", 1)[1].split()]
    words = sorted(word(v, q, n) for v in verts)
    # candidate check only (rule 10 style: search-side), verify.py decides
    bad = next(((a, b) for i, a in enumerate(words) for b in words[i + 1:] if lev(a, b) < d), None)
    digest = hashlib.sha256("".join(w + "\n" for w in words).encode()).hexdigest()
    art = os.path.join(OUT, f"edit{q}_n{n}_d{d}_fullgraph_seed{seed}_size{len(words)}.txt")
    with open(art, "w") as fh:
        fh.write(f"# edit-metric code, q = {q}, fixed length\n# n = {n}\n# d = {d}\n# q = {q}\n"
                 f"# method: search/fullgraph_q.py, mis_ils on the full conflict graph, seed {seed}, {secs}s\n"
                 f"# size = {len(words)}\n# sha256(codewords) = {digest}\n")
        fh.write("".join(w + "\n" for w in words))
    rec = {"q": q, "n": n, "d": d, "seed": seed, "seconds": secs, "size": len(words),
           "pairwise_ok": bad is None, "wall": round(wall, 1), "artifact": os.path.relpath(art, ROOT),
           "sha256": digest, "date": time.strftime("%Y-%m-%d")}
    return rec


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cells", required=True, help="q:n:d,q:n:d,...")
    ap.add_argument("--seeds", type=int, default=1)
    ap.add_argument("--seed0", type=int, default=1)
    ap.add_argument("--seconds", type=int, default=120)
    ap.add_argument("--workers", type=int, default=4)
    a = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    cells = [tuple(int(x) for x in c.split(":")) for c in a.cells.split(",")]
    for q, n, d in cells:      # build graphs serially, once
        graph_path(q, n, d)
    jobs = [(q, n, d, s, a.seconds) for q, n, d in cells for s in range(a.seed0, a.seed0 + a.seeds)]
    with Pool(a.workers) as pool, open(LOG, "a") as log:
        for rec in pool.imap_unordered(run, jobs):
            log.write(json.dumps(rec) + "\n")
            log.flush()
            print(json.dumps(rec), flush=True)


if __name__ == "__main__":
    sys.exit(main())
