"""Prescribed-symmetry search for fixed-length q-ary edit-metric codes, q = 2 or 3:
codes that are unions of orbits of a group H of symbol permutations x reversal.

WHY. Rung 2 of the ladder for the binary and ternary tables.  At q = 3, n >= 11 the
full conflict graph (177k / 531k vertices) is past mis_ils's V^2 adjacency bitset,
while an orbit graph under |H| = 4..12 is 15k-130k vertices.  The same method moved
the quaternary cells (paper, arXiv:2609.39081).  Symbol permutations and reversal
preserve Levenshtein distance, so every H acts by automorphisms.

    cc -O2 -o search/trie_orbit search/trie_orbit.c
    python search/orbit_q.py --q 3 --n 11 --d 3 --groups s3rev,c3rev --seeds 1 --seconds 600

Results go to results/orbit_q/ (artifact per run) and results/orbit_q_runs.jsonl.
Nothing printed here is a result until verify.py accepts the artifact.
"""
from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import os
import subprocess
import time
from multiprocessing import Pool

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "results", "orbit_q")
LOG = os.path.join(ROOT, "results", "orbit_q_runs.jsonl")


def gens(q, name):
    """Generators (perm, rev) of the named subgroup of S_q x C_2."""
    idp = tuple(range(q))
    sw = tuple([1, 0] + list(range(2, q)))
    c3 = (1, 2, 0)
    table = {
        "rev": [(idp, 1)],
        "swap": [(sw, 0)],
        "swaprev": [(sw, 1)],
        "swap+rev": [(sw, 0), (idp, 1)],
        "c3": [(c3, 0)],
        "c3rev": [(c3, 0), (idp, 1)],
        "s3": [(c3, 0), (sw, 0)],
        "s3rev": [(c3, 0), (sw, 0), (idp, 1)],
        "d3tw": [(c3, 0), (sw, 1)],
    }
    g = table[name]
    if q == 2 and any(p == c3 for p, _ in g):
        raise ValueError(f"group {name} needs q >= 3")
    return g


def closure(q, generators):
    idp = tuple(range(q))
    els = {(idp, 0)}
    frontier = list(els)
    while frontier:
        new = []
        for p, r in frontier:
            for gp, gr in generators:
                e = (tuple(gp[p[i]] for i in range(q)), r ^ gr)
                if e not in els:
                    els.add(e)
                    new.append(e)
        frontier = new
    return sorted(els)


def build(q, n, d, group):
    els = closure(q, gens(q, group))
    tag = f"q{q}_n{n}_d{d}_{group}"
    gfile = os.path.join(OUT, f"group_{tag}.txt")
    graph = os.path.join(OUT, f"graph_{tag}.txt")
    orbs = os.path.join(OUT, f"orbits_{tag}.txt")
    if not os.path.exists(graph):
        with open(gfile, "w") as fh:
            for p, r in els:
                fh.write(" ".join(map(str, p)) + f" {r}\n")
        t0 = time.time()
        p = subprocess.run([os.path.join(ROOT, "search", "trie_orbit"), str(q), str(n), str(d),
                            gfile, graph + ".tmp", orbs + ".tmp"], capture_output=True, text=True, check=True)
        os.replace(orbs + ".tmp", orbs)
        os.replace(graph + ".tmp", graph)
        print(f"built {tag}: {p.stderr.strip()} in {time.time() - t0:.0f}s", flush=True)
    return graph, orbs, len(els)


def word(idx, q, n):
    s = []
    for _ in range(n):
        s.append(str(idx % q))
        idx //= q
    return "".join(reversed(s))


def run(job):
    q, n, d, group, seed, secs = job
    graph, orbs, order = build(q, n, d, group)
    t0 = time.time()
    p = subprocess.run([os.path.join(ROOT, "search", "mis_ils"), graph, str(seed), str(secs), "--quiet"],
                       capture_output=True, text=True, check=True)
    wall = time.time() - t0
    best = [l for l in p.stdout.splitlines() if l.startswith("BEST")][-1]
    chosen = [int(x) for x in best.split(":", 1)[1].split()]
    with open(orbs) as fh:
        orbit_words = [list(map(int, l.split())) for l in fh]
    words = sorted(word(v, q, n) for o in chosen for v in orbit_words[o])
    digest = hashlib.sha256("".join(w + "\n" for w in words).encode()).hexdigest()
    os.makedirs(OUT, exist_ok=True)
    art = os.path.join(OUT, f"edit{q}_n{n}_d{d}_{group}_seed{seed}_size{len(words)}.txt")
    with open(art, "w") as fh:
        fh.write(f"# edit-metric code, q = {q}, fixed length\n# n = {n}\n# d = {d}\n# q = {q}\n"
                 f"# invariant under {group} (|H| = {order}, symbol permutations x reversal)\n"
                 f"# method: search/orbit_q.py, mis_ils on the orbit graph, seed {seed}, {secs}s\n"
                 f"# size = {len(words)}\n# sha256(codewords) = {digest}\n")
        fh.write("".join(w + "\n" for w in words))
    return {"q": q, "n": n, "d": d, "group": group, "order": order, "seed": seed, "seconds": secs,
            "size": len(words), "orbits": len(chosen), "wall": round(wall, 1),
            "artifact": os.path.relpath(art, ROOT), "sha256": digest, "date": time.strftime("%Y-%m-%d")}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--q", type=int, required=True)
    ap.add_argument("--n", type=int, required=True)
    ap.add_argument("--d", type=str, required=True, help="one d or a comma list")
    ap.add_argument("--groups", required=True)
    ap.add_argument("--seeds", type=int, default=1)
    ap.add_argument("--seed0", type=int, default=1)
    ap.add_argument("--seconds", type=int, default=600)
    ap.add_argument("--workers", type=int, default=3)
    a = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    ds = [int(x) for x in a.d.split(",")]
    groups = a.groups.split(",")
    for d, g in itertools.product(ds, groups):
        build(a.q, a.n, d, g)
    jobs = [(a.q, a.n, d, g, s, a.seconds) for d in ds for g in groups
            for s in range(a.seed0, a.seed0 + a.seeds)]
    with Pool(a.workers) as pool, open(LOG, "a") as log:
        for rec in pool.imap_unordered(run, jobs):
            log.write(json.dumps(rec) + "\n")
            log.flush()
            print(json.dumps(rec), flush=True)


if __name__ == "__main__":
    main()
