"""Drive the C iterated-local-search engine (search/mis_ils.c) on orbit graphs of
quaternary single-edit codes, expand the result to a code, and verify it.

WHY. The V_4-invariant 120 at n=6 (best/edit4_n6_V4_size120.txt) was found by an
ARW-style operator that the repo's own Python searcher does not implement fast
enough to reach. This driver builds exactly the same orbit instances as
search/edit4_diag.py (so the orbit-reduction lemma applies unchanged) and hands
them to a C engine that runs ~1000x more local-search iterations per second.

Targets it is written for (n=6, record 120):
  * the three free-acting order-4 classes: all orbits have size 4, so the next
    attainable size is 124 = 31 orbits;
  * the (0 1)(2 3) class: all 2046 orbits have size 2, next above 120 is 122;
  * the trivial class (raw 4096-vertex graph) seeded from a 120 code: 121.

Every candidate is verified by `verify_edit4` (ball disjointness) AND by the
independent pairwise Levenshtein check before anything is written to best/.

Run:
    cc -O2 -o search/mis_ils search/mis_ils.c
    python search/edit4_cmis.py --n 6 --list                    # class census
    python search/edit4_cmis.py --n 6 --classes 20,21,22 --seeds 8 --seconds 120
    python search/edit4_cmis.py --n 6 --classes 0 --init best/edit4_n6_V4_size120.txt
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

from search.edit4 import Edit4Graph, verify_edit4  # noqa: E402
from search.edit4_diag import (  # noqa: E402
    attainable_sizes, cell_instance, conjugacy_classes,
)
from search.edit4_orbit_ls import independent_check  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS_DIR = os.environ.get("FUNMATH_RESULTS_DIR", os.path.join(ROOT, "results"))
BEST_DIR = os.path.join(ROOT, "best")
ENGINE = os.path.join(ROOT, "search", "mis_ils")
SCRATCH = os.environ.get("FUNMATH_SCRATCH", os.path.join(RESULTS_DIR, "cmis_tmp"))

RECORD = {6: 120, 7: 356, 8: 1132, 9: 3451}   # best known lower bounds (ours at 6; table elsewhere)


def ensure_engine():
    src = ENGINE + ".c"
    if not os.path.exists(ENGINE) or os.path.getmtime(ENGINE) < os.path.getmtime(src):
        subprocess.run(["cc", "-O2", "-o", ENGINE, src], check=True)


def ordered_classes():
    """Conjugacy classes in a stable order: by usable-orbit count is not known
    until built, so order by (order, diagonal, name) which is deterministic."""
    cls = conjugacy_classes()
    cls.sort(key=lambda c: (-c["order"], c["diagonal"], c["name"]))
    return cls


def write_graph(path, weights, adj):
    edges = [(a, b) for a in range(len(adj)) for b in adj[a] if a < b]
    with open(path, "w") as fh:
        fh.write(f"{len(weights)} {len(edges)}\n")
        fh.write(" ".join(str(x) for x in weights) + "\n")
        fh.write("".join(f"{a} {b}\n" for a, b in edges))
    return len(edges)


def run_engine(graph_path, seed, seconds, target=None, init_path=None):
    cmd = [ENGINE, graph_path, str(seed), str(seconds), "--quiet"]
    if target:
        cmd += ["--target", str(target)]
    if init_path:
        cmd += ["--init", init_path]
    out = subprocess.run(cmd, capture_output=True, text=True, check=True).stdout
    best_line = next(l for l in out.splitlines() if l.startswith("BEST"))
    head, _, tail = best_line.partition(":")
    wt, cnt = int(head.split()[1]), int(head.split()[2])
    sel = [int(x) for x in tail.split()]
    iters = next((int(l.split()[1]) for l in out.splitlines() if l.startswith("iters")), 0)
    assert len(sel) == cnt
    return wt, sel, iters


def read_graph(path):
    with open(path) as fh:
        k, m = map(int, fh.readline().split())
        weights = list(map(int, fh.readline().split()))
        adj = [[] for _ in range(k)]
        for line in fh:
            x, y = map(int, line.split())
            adj[x].append(y)
            adj[y].append(x)
    return weights, adj, m


def cached_instance(g, ci, H, cache_stem=None):
    """Orbit graph + usable-orbit list for class ci, cached under graphs/ because
    the O(k^2) build is minutes at n=7. The cache is keyed by (n, class index in
    ordered_classes()); the orbit list is what turns a vertex set back into
    codewords, so it is stored next to the graph."""
    gdir = os.path.join(ROOT, "graphs")
    os.makedirs(gdir, exist_ok=True)
    gpath = os.path.join(gdir, f"{cache_stem or f'edit4_n{g.n}_orbit_class'}{ci}.txt")
    opath = gpath[:-4] + ".orbits.json"
    if os.path.exists(gpath) and os.path.exists(opath):
        # header + weights only: the engine reads the edge list itself, and a
        # dense n=8 instance (8e7 edges) does not fit in Python lists
        with open(gpath) as fh:
            k, m = map(int, fh.readline().split())
            weights = list(map(int, fh.readline().split()))
        usable = [tuple(o) for o in json.load(open(opath))]
        assert len(usable) == len(weights) == k
        return usable, None, weights, gpath, m
    usable, adj, weights, orbits = cell_instance(g, H)
    m = write_graph(gpath, weights, adj)
    with open(opath, "w") as fh:
        json.dump([list(o) for o in usable], fh)
    return usable, adj, weights, gpath, m


def code_to_orbit_init(words, usable, g):
    """Map a code to orbit indices; only orbits fully contained are kept."""
    idx = {w: i for i, w in enumerate(g.words)}
    have = {idx[w] for w in words}
    sel = [i for i, orb in enumerate(usable) if all(v in have for v in orb)]
    covered = sum(len(usable[i]) for i in sel)
    return sel, covered


def save_artifact(n, words, size, label, dist=3, ci=None):
    body = "".join(w + "\n" for w in sorted(words))
    d = hashlib.sha256(body.encode()).hexdigest()
    os.makedirs(BEST_DIR, exist_ok=True)
    tag = f"_class{ci}" if ci is not None else ""
    stem = f"edit4_n{n}_cmis" if dist == 3 else f"editq_n{n}_d{dist}"
    path = os.path.join(BEST_DIR, f"{stem}_size{size}{tag}.txt")
    # never overwrite a different code of the same size and class (this
    # silently replaced two artifacts on 2026-09-08); add an _alt suffix
    k = 2
    while os.path.exists(path):
        existing = [l.strip() for l in open(path) if l.strip() and not l.startswith("#")]
        if sorted(existing) == sorted(words):
            break
        path = os.path.join(BEST_DIR, f"{stem}_size{size}{tag}_alt{k}.txt")
        k += 1
    with open(path, "w") as fh:
        fh.write(f"# quaternary edit-metric code, minimum Levenshtein distance {dist}\n")
        fh.write(f"# n = {n}\n# d = {dist}\n# q = 4\n# size = {size}\n")
        fh.write(f"# invariant under {label}\n")
        fh.write(f"# found by search/edit4_cmis.py (C ARW local search on the orbit graph)\n")
        fh.write(f"# sha256(codewords) = {d}\n")
        fh.write(body)
    return path, d


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=6)
    ap.add_argument("--d", type=int, default=3,
                    help="minimum edit distance; d != 3 requires cached instances "
                         "built by search/orbit_graph.c (graphs/editq_n{n}_d{d}_orbit_class{ci}.*)")
    ap.add_argument("--list", action="store_true")
    ap.add_argument("--classes", type=str, default="",
                    help="comma-separated class indices from --list")
    ap.add_argument("--seeds", type=int, default=8)
    ap.add_argument("--seed0", type=int, default=1)
    ap.add_argument("--seconds", type=float, default=120.0)
    ap.add_argument("--workers", type=int, default=6)
    ap.add_argument("--init", type=str, default=None,
                    help="code file to seed from (orbits fully covered are used)")
    ap.add_argument("--target", type=int, default=None,
                    help="stop a seed early once this weight is reached")
    ap.add_argument("--tag", type=str, default="")
    ap.add_argument("--export", action="store_true",
                    help="build and cache the orbit graphs of the classes, then exit")
    a = ap.parse_args()

    ensure_engine()
    os.makedirs(SCRATCH, exist_ok=True)
    t0 = time.time()
    classes = ordered_classes()
    if a.d == 3:
        record = RECORD[a.n]
        cache_stem = f"edit4_n{a.n}_orbit_class"
    else:
        from search.editq import TABLE
        record = TABLE[(a.n, a.d)][0]
        cache_stem = f"editq_n{a.n}_d{a.d}_orbit_class"
    # At n = 8 the Python word graph (65,536 big-int bitmasks) is too heavy to
    # build; when every requested class is cached (built by search/orbit_graph.c)
    # only the word list is needed for expansion, and verify_edit4 is pairwise.
    want_pre = [int(x) for x in a.classes.split(",") if x.strip()]
    gdir = os.path.join(ROOT, "graphs")
    all_cached = bool(want_pre) and not a.list and all(
        os.path.exists(os.path.join(gdir, f"{cache_stem}{ci}.txt")) and
        os.path.exists(os.path.join(gdir, f"{cache_stem}{ci}.orbits.json"))
        for ci in want_pre)
    if a.d != 3 and not all_cached:
        ap.error("d != 3 needs cached instances for every requested class")
    if all_cached:
        import itertools
        class _Words:
            pass
        g = _Words()
        g.n = a.n
        g.words = ["".join(w) for w in itertools.product("0123", repeat=a.n)]
    else:
        g = Edit4Graph(a.n)

    if a.list:
        for i, c in enumerate(classes):
            usable, adj, weights, orbits = cell_instance(g, c["H"])
            sizes = Counter(weights)
            reach = attainable_sizes(sizes, record + 12)
            nxt = next((t for t in range(record + 1, record + 13) if (reach >> t) & 1), None)
            print(f"{i:>2} |H|={c['order']:>2} {'DIAG' if c['diagonal'] else 'prod'} "
                  f"usable {len(usable):>5} edges {sum(map(len, adj))//2:>7} "
                  f"sizes {dict(sizes)} next>{record}: {nxt}  {c['name']}")
        return 0

    want = [int(x) for x in a.classes.split(",") if x.strip()]
    if not want:
        ap.error("--classes required")
    stamp = time.strftime("%Y%m%dT%H%M%S")
    log = {"n": a.n, "engine": "search/mis_ils.c", "seeds": a.seeds, "seed0": a.seed0,
           "seconds": a.seconds, "init": a.init, "target": a.target, "cells": []}
    init_words = None
    if a.init:
        init_words = [l.strip() for l in open(a.init) if l.strip() and not l.startswith("#")]
        ok, sz, why = verify_edit4(init_words, a.n)
        print(f"init code {a.init}: valid={ok} size={sz} {why}")

    overall_best = 0
    lp = None
    for ci in want:
        c = classes[ci]
        tb = time.time()
        usable, adj, weights, gpath, m = cached_instance(g, ci, c["H"], cache_stem)
        sizes = Counter(weights)
        reach = attainable_sizes(sizes, record + 12)
        nxt = next((t for t in range(record + 1, record + 13) if (reach >> t) & 1), None)
        if a.export:
            print(f"cached class {ci} ({len(usable)} vertices, {m} edges) at {gpath}", flush=True)
            continue
        init_path = None
        if init_words:
            sel0, covered = code_to_orbit_init(init_words, usable, g)
            init_path = os.path.join(SCRATCH, f"n{a.n}_class{ci}_init.txt")
            with open(init_path, "w") as fh:
                fh.write("\n".join(map(str, sel0)) + "\n")
            print(f"  init covers {len(sel0)} orbits = {covered} words in class {ci}")
        print(f"class {ci}: |H|={c['order']} {'DIAG' if c['diagonal'] else 'prod'} "
              f"{len(usable)} usable, {m} edges, sizes {dict(sizes)}, "
              f"next>{record}: {nxt}, build {time.time()-tb:.0f}s  {c['name']}", flush=True)
        seeds = list(range(a.seed0, a.seed0 + a.seeds))
        with ThreadPoolExecutor(max_workers=a.workers) as ex:
            res = list(ex.map(lambda s: run_engine(gpath, s, a.seconds, a.target, init_path), seeds))
        by_seed = [r[0] for r in res]
        bi = max(range(len(res)), key=lambda i: res[i][0])
        wt, sel, iters = res[bi]
        words = sorted(g.words[v] for i in sel for v in usable[i])
        if a.d == 3:
            ok, size, reason = verify_edit4(words, a.n)
            bad = independent_check(words) if ok else -1
        else:
            from verify_editq import check as _vq
            from verify_editq2 import validate as _vq2
            ok, size, reason = _vq(words, a.n, a.d)
            bad = 0 if (ok and _vq2(words, a.n, a.d, 4)[0]) else -1
        print(f"  by seed {by_seed}  (iters best seed {iters})", flush=True)
        print(f"  best {wt}  verify_edit4={ok} size={size} {reason}  lev-bad-pairs={bad}", flush=True)
        rec = {"class_index": ci, "name": c["name"], "order": c["order"],
               "diagonal": c["diagonal"], "orbits_usable": len(usable), "orbit_edges": m,
               "orbit_sizes": dict(sizes), "by_seed": by_seed, "seeds": seeds,
               "best": wt, "verified": ok and bad == 0, "size": size,
               "next_attainable_above_record": nxt,
               "wall_seconds": round(time.time() - tb, 1)}
        if ok and bad == 0 and size > record:
            path, d = save_artifact(a.n, words, size, c["name"], a.d, ci)
            rec["artifact"] = path
            rec["sha256"] = d
            print(f"  *** {size} > record {record}: saved {path} sha256={d}\n"
                  f"      CLAUDE.md rule 6: treat as a bug until verify.py passes in a fresh process",
                  flush=True)
        elif ok and bad == 0 and size == record:
            path, d = save_artifact(a.n, words, size, c["name"] + " (ties record)")
            rec["artifact"] = path
            rec["sha256"] = d
        overall_best = max(overall_best, size if ok else 0)
        log["cells"].append(rec)
        log["elapsed_seconds"] = round(time.time() - t0, 1)
        os.makedirs(RESULTS_DIR, exist_ok=True)
        lp = os.path.join(RESULTS_DIR, f"edit4_cmis_n{a.n}{('_' + a.tag) if a.tag else ''}_{stamp}.json")
        with open(lp, "w") as fh:
            json.dump(log, fh, indent=2)
    print(f"\nbest overall {overall_best}  record {record}  log {lp}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
