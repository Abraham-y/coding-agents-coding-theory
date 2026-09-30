"""Sweep every invariance class of S_4 x C_2 at a given n with the full ARW
operator, a wide pair-candidate cap and a real budget.

Written after an adversarial review found a 120-codeword code at n=6 in a class
this repository had recorded as topping out at 112. That result was a seed-budget
artifact: four seeds per class is a statement about the search, not the cell. This
script exists to stop making that mistake at n=7 and n=8.

Run:  python search/sweep_all.py --n 7 --target 356 --seeds 8 --seconds 600
"""
from __future__ import annotations
import argparse, hashlib, json, os, sys, time
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from multiprocessing import Pool
import search.edit4_orbit_ls as L
from search.edit4 import Edit4Graph, verify_edit4
from search.edit4_diag import conjugacy_classes, cell_instance

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
_W = _A = None

def _init(w, a, cap):
    global _W, _A
    _W, _A = w, a
    L.WeightedSearcher.PAIR_CANDS = cap

def _run(job):
    seed, secs = job
    sel, wt, _ = L.ils(_W, _A, seed=seed, time_limit=secs)
    return {"seed": seed, "weight": wt, "sel": sel}

def lev(a, b):
    prev = list(range(len(b) + 1))
    for i, ca in enumerate(a, 1):
        cur = [i]
        for j, cb in enumerate(b, 1):
            cur.append(min(prev[j] + 1, cur[j - 1] + 1, prev[j - 1] + (ca != cb)))
        prev = cur
    return prev[-1]

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, required=True)
    ap.add_argument("--target", type=int, required=True, help="record to beat")
    ap.add_argument("--seeds", type=int, default=8)
    ap.add_argument("--seconds", type=float, default=600.0)
    ap.add_argument("--workers", type=int, default=5)
    ap.add_argument("--cap", type=int, default=200)
    ap.add_argument("--max-usable", type=int, default=100000)
    a = ap.parse_args()

    t0 = time.time()
    g = Edit4Graph(a.n)
    print(f"Edit4 G({a.n}): {g.num_vertices} vertices; target to beat {a.target}", flush=True)
    classes = sorted(conjugacy_classes(), key=lambda c: -c["order"])
    out, best_overall = [], 0
    for c in classes:
        usable, adj, w, orbits = cell_instance(g, c["H"])
        if len(usable) > a.max_usable or not usable:
            continue
        with Pool(a.workers, initializer=_init, initargs=(w, adj, a.cap)) as pool:
            res = pool.map(_run, [(s, a.seconds) for s in range(1, a.seeds + 1)])
        res.sort(key=lambda r: -r["weight"])
        top = res[0]
        words = sorted(g.words[v] for i in top["sel"] for v in usable[i])
        ok, size, _ = verify_edit4(words, a.n)
        flag = "DIAG" if c["diagonal"] else "    "
        print(f"  {flag} |H|={c['order']:>2} {len(usable):>5} usable  best {top['weight']:>5}"
              f"  seeds {[r['weight'] for r in res]}  verified={ok}", flush=True)
        out.append({"name": c["name"], "order": c["order"], "diagonal": c["diagonal"],
                    "usable": len(usable), "best": top["weight"],
                    "by_seed": [r["weight"] for r in res], "verified": ok})
        if ok and top["weight"] > best_overall:
            best_overall = top["weight"]
        if ok and top["weight"] > a.target:
            mind = min(lev(words[i], words[j]) for i in range(len(words))
                       for j in range(i + 1, len(words)))
            print(f"       *** {size} BEATS {a.target}; min edit distance {mind} ***", flush=True)
            if mind >= 3:
                body = "".join(x + "\n" for x in words)
                d = hashlib.sha256(body.encode()).hexdigest()
                p = os.path.join(ROOT, "best", f"edit4_n{a.n}_size{size}.txt")
                with open(p, "w") as fh:
                    fh.write("# quaternary single-edit-correcting code (min edit distance 3)\n")
                    fh.write(f"# n = {a.n}, alphabet = 0123, size = {size}\n")
                    fh.write(f"# invariant under {c['name']}\n")
                    fh.write(f"# sha256(codewords) = {d}\n{body}")
                print(f"       saved {p}", flush=True)
        json.dump({"n": a.n, "target": a.target, "seeds": a.seeds, "seconds": a.seconds,
                   "pair_cands": a.cap, "cells": out, "best": best_overall,
                   "elapsed_seconds": round(time.time() - t0, 1)},
                  open(os.path.join(ROOT, "results", f"sweep_all_n{a.n}.json"), "w"), indent=2)
    print(f"\nbest over all classes: {best_overall}   target {a.target}", flush=True)

if __name__ == "__main__":
    main()
