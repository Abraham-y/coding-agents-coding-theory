"""Iterated local search for maximum independent set in G(n, s).

Rung 3 of the escalation ladder. Structure follows Andrade-Resende-Werneck:
maximise greedily, apply (1,2)-swaps until no improvement, perturb, repeat.

  1-improvement : a vertex with no neighbour in S can just be added.
  2-improvement : if u in S has two NON-ADJACENT "1-tight" neighbours v,w
                  (each having u as its only S-neighbour), then
                  S - u + v + w is one larger.
  plateau move  : swap u out for a single 1-tight v. Size unchanged, but it
                  moves the search across a plateau, which is where MIS local
                  search spends most of its life.

Every result is a set of codewords and is worthless until verify.py and
verify2.py both accept it. This module never decides anything.

Run:
    python search/ils.py --n 12 --s 2 --seconds 60 --seed 1
    python search/ils.py --n 12 --s 2 --seconds 600 --workers 10
"""

from __future__ import annotations

import argparse
import json
import os
import random
import sys
import time
from multiprocessing import Pool

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from search.graph import load_graph  # noqa: E402

RESULTS_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "results")


class Searcher:
    def __init__(self, graph, seed: int):
        self.g = graph
        self.rng = random.Random(seed)
        self.seed = seed
        self.nv = graph.num_vertices
        self.adj = graph.adj
        # neighbour lists for incremental tightness updates
        self.nbrs = []
        for v in range(self.nv):
            m = self.adj[v]
            lst = []
            while m:
                low = m & -m
                lst.append(low.bit_length() - 1)
                m ^= low
            self.nbrs.append(lst)
        self.in_S = bytearray(self.nv)
        self.tight = [0] * self.nv
        self.size = 0
        # S and the free set (tight == 0, not in S) are maintained incrementally;
        # scanning all nv vertices per move was the bottleneck.
        self.S = set()
        self.free = set(range(self.nv))

    # -- primitive moves --------------------------------------------------
    def add(self, v: int):
        self.in_S[v] = 1
        self.size += 1
        self.S.add(v)
        self.free.discard(v)
        t = self.tight
        free = self.free
        for w in self.nbrs[v]:
            t[w] += 1
            if t[w] == 1:
                free.discard(w)

    def remove(self, v: int):
        self.in_S[v] = 0
        self.size -= 1
        self.S.discard(v)
        t = self.tight
        free = self.free
        if t[v] == 0:
            free.add(v)
        for w in self.nbrs[v]:
            t[w] -= 1
            if t[w] == 0 and not self.in_S[w]:
                free.add(w)

    def clear(self):
        for v in range(self.nv):
            self.in_S[v] = 0
            self.tight[v] = 0
        self.size = 0
        self.S = set()
        self.free = set(range(self.nv))

    def current(self) -> list[int]:
        return sorted(self.S)

    # -- construction -----------------------------------------------------
    def maximize(self):
        """Add free vertices until none remain."""
        while self.free:
            v = self.rng.choice(tuple(self.free))
            self.add(v)

    # -- local search -----------------------------------------------------
    def two_improve(self) -> bool:
        """One pass of 2-improvements. Returns True if the set grew."""
        improved = False
        members = list(self.S)
        self.rng.shuffle(members)
        for u in members:
            if not self.in_S[u]:
                continue
            cands = [v for v in self.nbrs[u] if not self.in_S[v] and self.tight[v] == 1]
            if len(cands) < 2:
                continue
            found = None
            for i in range(len(cands)):
                ai = self.adj[cands[i]]
                for j in range(i + 1, len(cands)):
                    if not (ai >> cands[j]) & 1:
                        found = (cands[i], cands[j])
                        break
                if found:
                    break
            if found:
                self.remove(u)
                self.add(found[0])
                self.add(found[1])
                improved = True
        return improved

    def plateau(self, moves: int):
        """Size-preserving 1-for-1 swaps, to drift across the plateau."""
        for _ in range(moves):
            if not self.S:
                return
            u = self.rng.choice(tuple(self.S))
            cands = [v for v in self.nbrs[u] if not self.in_S[v] and self.tight[v] == 1]
            if not cands:
                continue
            v = self.rng.choice(cands)
            self.remove(u)
            self.add(v)
            self.maximize()

    def local_search(self):
        while True:
            self.maximize()
            if not self.two_improve():
                break

    # -- perturbation -----------------------------------------------------
    def perturb(self, strength: int):
        for _ in range(strength):
            v = self.rng.randrange(self.nv)
            if self.in_S[v]:
                continue
            for w in self.nbrs[v]:
                if self.in_S[w]:
                    self.remove(w)
            self.add(v)

    # -- driver -----------------------------------------------------------
    def run(self, seconds: float, target: int | None = None, report=None):
        t0 = time.time()
        self.clear()
        self.local_search()
        best = self.current()
        best_size = len(best)
        iters = 0
        since_improve = 0

        while time.time() - t0 < seconds:
            iters += 1
            self.perturb(strength=1 + self.rng.randrange(3))
            self.local_search()
            self.plateau(moves=self.rng.randrange(8))
            self.local_search()
            if self.size > best_size:
                best = self.current()
                best_size = self.size
                since_improve = 0
                if report:
                    report(best_size, time.time() - t0)
                if target is not None and best_size >= target:
                    break
            else:
                since_improve += 1
                # restart from the incumbent when stuck
                if since_improve > 400:
                    self.clear()
                    for v in best:
                        self.add(v)
                    since_improve = 0

        return best, best_size, iters, time.time() - t0


def _worker(args):
    n, s, seed, seconds, target = args
    g = load_graph(n, s)
    searcher = Searcher(g, seed)
    best, size, iters, elapsed = searcher.run(seconds, target)
    return {
        "seed": seed,
        "size": size,
        "iters": iters,
        "seconds": round(elapsed, 2),
        "words": g.to_words(sorted(best)),
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, required=True)
    ap.add_argument("--s", type=int, required=True)
    ap.add_argument("--seconds", type=float, default=60)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--workers", type=int, default=1)
    ap.add_argument("--target", type=int, default=None)
    ap.add_argument("--tag", type=str, default="ils")
    args = ap.parse_args()

    load_graph(args.n, args.s)  # warm the cache before forking

    jobs = [
        (args.n, args.s, args.seed + i, args.seconds, args.target)
        for i in range(args.workers)
    ]
    t0 = time.time()
    if args.workers == 1:
        results = [_worker(jobs[0])]
    else:
        with Pool(args.workers) as pool:
            results = pool.map(_worker, jobs)
    wall = time.time() - t0

    results.sort(key=lambda r: -r["size"])
    best = results[0]

    os.makedirs(RESULTS_DIR, exist_ok=True)
    stamp = time.strftime("%Y%m%dT%H%M%S")
    logpath = os.path.join(RESULTS_DIR, f"{args.tag}_n{args.n}_s{args.s}_{stamp}.json")
    record = {
        "method": "iterated local search (ARW-style 1,2-swaps + plateau)",
        "n": args.n,
        "s": args.s,
        "seconds_per_worker": args.seconds,
        "workers": args.workers,
        "base_seed": args.seed,
        "wall_seconds": round(wall, 2),
        "best_size": best["size"],
        "sizes_by_seed": {str(r["seed"]): r["size"] for r in results},
        "total_iters": sum(r["iters"] for r in results),
        "best_words": best["words"],
    }
    with open(logpath, "w") as fh:
        json.dump(record, fh, indent=2)

    print(f"n={args.n} s={args.s} workers={args.workers} wall={wall:.1f}s")
    print(f"best size : {best['size']}")
    print(f"by seed   : {record['sizes_by_seed']}")
    print(f"log       : {logpath}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
