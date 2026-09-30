"""Diagonal-symmetry sweep at n = 7, targeting E4(7,3) > 356.

Same method that produced E4(6,3) >= 116: prescribed symmetry over all 33
conjugacy classes of subgroups of S_4 x C_2, including the 11 diagonal ones no
earlier sweep could express. Differs from edit4_diag.py only in plumbing: the
orbit instance is built ONCE per class in the parent and shared with the seed
workers via a Pool initializer, because at n=7 the O(k^2) orbit-graph build is
minutes rather than the ~1s it was at n=6.

ARITHMETIC NOTE (matters for target selection): in a class where the group acts
freely all orbits have size |H|, so attainable code sizes are multiples of |H|.
356 = 4*89, so in a free order-4 class the record is attainable and beating it
needs 360 (+4). In an order-2 class with fixed points, 358 (+2) is reachable.
"""
from __future__ import annotations
import json, os, sys, time, hashlib
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from multiprocessing import Pool
from search.edit4 import Edit4Graph, verify_edit4
from search.edit4_diag import conjugacy_classes, cell_instance, attainable_sizes
from search.edit4_orbit_ls import ils
from collections import Counter

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
_W = _A = None

def _init(weights, adj):
    global _W, _A
    _W, _A = weights, adj

def _run(job):
    seed, seconds = job
    sel, wt, iters = ils(_W, _A, seed=seed, time_limit=seconds)
    return {"seed": seed, "weight": wt, "iters": iters, "sel": sel}

def lev(a, b):
    prev = list(range(len(b) + 1))
    for i, ca in enumerate(a, 1):
        cur = [i]
        for j, cb in enumerate(b, 1):
            cur.append(min(prev[j] + 1, cur[j - 1] + 1, prev[j - 1] + (ca != cb)))
        prev = cur
    return prev[-1]

def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=7)
    ap.add_argument("--target", type=int, default=356)
    ap.add_argument("--seeds", type=int, default=4)
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--max-seconds", type=float, default=600.0)
    ap.add_argument("--max-usable", type=int, default=99999,
                    help="skip classes whose orbit graph exceeds this (build cost)")
    a = ap.parse_args()

    t0 = time.time()
    g = Edit4Graph(a.n)
    print(f"Edit4 G({a.n}): {g.num_vertices} vertices, {len(g.cliques)} cliques "
          f"({time.time()-t0:.0f}s)", flush=True)
    classes = conjugacy_classes()
    # biggest groups first: cheapest builds, and they tell us fastest whether
    # the method transfers at all
    classes.sort(key=lambda c: -c["order"])
    out = []
    for c in classes:
        tb = time.time()
        usable, adj, weights, orbits = cell_instance(g, c["H"])
        if len(usable) > a.max_usable:
            print(f"  SKIP {c['name']}  ({len(usable)} usable > cap)", flush=True)
            continue
        sizes = Counter(weights)
        reach = attainable_sizes(sizes, a.target + 8)
        nxt = next((t for t in range(a.target + 1, a.target + 9) if (reach >> t) & 1), None)
        secs = max(30.0, min(a.max_seconds, len(usable) / 6.0))
        with Pool(a.workers, initializer=_init, initargs=(weights, adj)) as pool:
            res = pool.map(_run, [(s, secs) for s in range(1, a.seeds + 1)])
        res.sort(key=lambda r: -r["weight"])
        best = res[0]
        words = sorted(g.words[v] for i in best["sel"] for v in usable[i])
        ok, size, reason = verify_edit4(words, a.n)
        tag = "DIAG" if c["diagonal"] else "    "
        print(f"  {tag} |H|={c['order']:>2} {len(usable):>5} usable  build {time.time()-tb:.0f}s  "
              f"{secs:.0f}s x{a.seeds}  best {best['weight']:>5}  next>{a.target}: {nxt}  "
              f"verified={ok} {reason}", flush=True)
        rec = {"name": c["name"], "order": c["order"], "diagonal": c["diagonal"],
               "orbits_total": len(orbits), "orbits_usable": len(usable),
               "orbit_sizes": dict(sizes), "by_seed": [r["weight"] for r in res],
               "best": best["weight"], "verified": ok,
               "next_attainable_above_target": nxt, "seconds_per_seed": secs}
        out.append(rec)
        if ok and best["weight"] > a.target:
            mind = min(lev(words[i], words[j]) for i in range(len(words))
                       for j in range(i + 1, len(words)))
            print(f"       *** EXCEEDS {a.target}: {size}, min edit distance {mind} ***", flush=True)
            if mind >= 3:
                body = "".join(w + "\n" for w in words)
                d = hashlib.sha256(body.encode()).hexdigest()
                p = os.path.join(ROOT, "best", f"edit4_n{a.n}_size{size}.txt")
                with open(p, "w") as fh:
                    fh.write(f"# quaternary single-edit-correcting code (min edit distance 3)\n")
                    fh.write(f"# n = {a.n}, size = {size}, invariant under {c['name']}\n")
                    fh.write(f"# sha256(codewords) = {d}\n")
                    fh.write(body)
                print(f"       saved {p}  sha256={d}", flush=True)
        json.dump({"n": a.n, "target": a.target, "cells": out,
                   "elapsed_seconds": round(time.time()-t0, 1)},
                  open(os.path.join(ROOT, "results", f"edit4_n{a.n}_diag.json"), "w"), indent=2)
    print(f"\nbest overall: {max((r['best'] for r in out), default=0)}  target {a.target}")

if __name__ == "__main__":
    main()
