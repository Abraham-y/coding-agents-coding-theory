"""Prescribed-automorphism search: build codes as unions of symmetry orbits.

G(n, s) has the Klein four-group {id, complement, reversal, complement*reversal}
as automorphisms, because LCS is invariant under applying the same bijection to
both words and under reversing both words.

If we DEMAND that the code be invariant under a subgroup H, the code becomes a
union of H-orbits, and the search collapses from "pick vertices out of 2^n" to
"pick orbits out of 2^n/|H|", with each orbit carrying weight |orbit|. That is a
maximum WEIGHT independent set on a graph 2x to 4x smaller.

This is a restriction: the optimal code need not be H-invariant, so this can only
ever find codes, never prove bounds. But it is exactly how most record covering
codes and designs were found, and it is the rung the original ladder omitted.

An orbit is usable only if it is internally conflict-free; orbits O, O' conflict
if any member of O conflicts with any member of O'.

Run:
    python search/orbit.py --n 12 --s 2 --group complement --seconds 120 --workers 10
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

from search.graph import load_graph, orbits  # noqa: E402

RESULTS_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "results")


class OrbitGraph:
    def __init__(self, n: int, s: int, group: str):
        self.n = n
        self.s = s
        self.group = group
        g = load_graph(n, s)
        self.base = g

        all_orbits = orbits(n, group)
        # keep only orbits with no internal conflict
        usable = []
        for orb in all_orbits:
            ok = True
            for i in range(len(orb)):
                for j in range(i + 1, len(orb)):
                    if (g.adj[orb[i]] >> orb[j]) & 1:
                        ok = False
                        break
                if not ok:
                    break
            if ok:
                usable.append(orb)

        self.orbits = usable
        self.num_vertices = len(usable)
        self.weight = [len(o) for o in usable]
        self.dropped = len(all_orbits) - len(usable)

        # orbit masks, then orbit-level adjacency
        masks = []
        for orb in usable:
            m = 0
            for v in orb:
                m |= 1 << v
            masks.append(m)
        self.masks = masks

        # neighbourhood of an orbit = union of member neighbourhoods
        nbr_union = []
        for orb in usable:
            m = 0
            for v in orb:
                m |= g.adj[v]
            nbr_union.append(m)

        adj = [0] * self.num_vertices
        for i in range(self.num_vertices):
            ni = nbr_union[i]
            for j in range(i + 1, self.num_vertices):
                if ni & masks[j]:
                    adj[i] |= 1 << j
                    adj[j] |= 1 << i
        self.adj = adj

    def words_of(self, chosen) -> list[str]:
        out = []
        for i in chosen:
            for v in self.orbits[i]:
                out.append(self.base.words[v])
        return out

    def total_weight(self, chosen) -> int:
        return sum(self.weight[i] for i in chosen)


class WeightedSearcher:
    """ILS for maximum WEIGHT independent set on the orbit graph."""

    def __init__(self, og: OrbitGraph, seed: int):
        self.og = og
        self.rng = random.Random(seed)
        self.nv = og.num_vertices
        self.adj = og.adj
        self.w = og.weight
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
        self.S = set()
        self.free = set(range(self.nv))
        self.wt = 0

    def add(self, v):
        self.in_S[v] = 1
        self.S.add(v)
        self.wt += self.w[v]
        self.free.discard(v)
        t, free = self.tight, self.free
        for u in self.nbrs[v]:
            t[u] += 1
            if t[u] == 1:
                free.discard(u)

    def remove(self, v):
        self.in_S[v] = 0
        self.S.discard(v)
        self.wt -= self.w[v]
        t, free = self.tight, self.free
        if t[v] == 0:
            free.add(v)
        for u in self.nbrs[v]:
            t[u] -= 1
            if t[u] == 0 and not self.in_S[u]:
                free.add(u)

    def clear(self):
        for v in range(self.nv):
            self.in_S[v] = 0
            self.tight[v] = 0
        self.S = set()
        self.free = set(range(self.nv))
        self.wt = 0

    def maximize(self):
        while self.free:
            # prefer heavy orbits, break ties randomly
            best, bw = None, -1
            sample = self.free if len(self.free) <= 24 else self.rng.sample(tuple(self.free), 24)
            for v in sample:
                if self.w[v] > bw:
                    best, bw = v, self.w[v]
            self.add(best)

    def two_improve(self) -> bool:
        improved = False
        for u in list(self.S):
            if not self.in_S[u]:
                continue
            cands = [v for v in self.nbrs[u] if not self.in_S[v] and self.tight[v] == 1]
            if not cands:
                continue
            found = None
            for i in range(len(cands)):
                ai = self.adj[cands[i]]
                for j in range(i + 1, len(cands)):
                    if not (ai >> cands[j]) & 1:
                        if self.w[cands[i]] + self.w[cands[j]] > self.w[u]:
                            found = (cands[i], cands[j])
                            break
                if found:
                    break
            if found is None:
                # single heavier swap
                for v in cands:
                    if self.w[v] > self.w[u]:
                        found = (v,)
                        break
            if found:
                self.remove(u)
                for v in found:
                    self.add(v)
                improved = True
        return improved

    def local_search(self):
        while True:
            self.maximize()
            if not self.two_improve():
                break

    def perturb(self, strength):
        for _ in range(strength):
            v = self.rng.randrange(self.nv)
            if self.in_S[v]:
                continue
            for u in self.nbrs[v]:
                if self.in_S[u]:
                    self.remove(u)
            self.add(v)

    def run(self, seconds):
        t0 = time.time()
        self.clear()
        self.local_search()
        best, best_wt = set(self.S), self.wt
        stall = 0
        while time.time() - t0 < seconds:
            self.perturb(1 + self.rng.randrange(3))
            self.local_search()
            if self.wt > best_wt:
                best, best_wt = set(self.S), self.wt
                stall = 0
            else:
                stall += 1
                if stall > 300:
                    self.clear()
                    for v in best:
                        self.add(v)
                    stall = 0
        return best, best_wt


def _worker(args):
    n, s, group, seed, seconds = args
    og = OrbitGraph(n, s, group)
    ws = WeightedSearcher(og, seed)
    chosen, wt = ws.run(seconds)
    return {
        "seed": seed,
        "size": wt,
        "orbits_used": len(chosen),
        "words": sorted(og.words_of(chosen)),
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, required=True)
    ap.add_argument("--s", type=int, required=True)
    ap.add_argument("--group", choices=["complement", "reversal", "both"], default="complement")
    ap.add_argument("--seconds", type=float, default=60)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--workers", type=int, default=1)
    ap.add_argument("--tag", type=str, default="orbit")
    args = ap.parse_args()

    og = OrbitGraph(args.n, args.s, args.group)
    print(
        f"orbit graph: {og.num_vertices} usable orbits "
        f"({og.dropped} dropped for internal conflict), "
        f"weights {sorted(set(og.weight))}"
    )

    jobs = [(args.n, args.s, args.group, args.seed + i, args.seconds) for i in range(args.workers)]
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
    logpath = os.path.join(
        RESULTS_DIR, f"{args.tag}_{args.group}_n{args.n}_s{args.s}_{stamp}.json"
    )
    record = {
        "method": f"prescribed-symmetry orbit MWIS, group={args.group}",
        "n": args.n,
        "s": args.s,
        "group": args.group,
        "usable_orbits": og.num_vertices,
        "dropped_orbits": og.dropped,
        "seconds_per_worker": args.seconds,
        "workers": args.workers,
        "base_seed": args.seed,
        "wall_seconds": round(wall, 2),
        "best_size": best["size"],
        "sizes_by_seed": {str(r["seed"]): r["size"] for r in results},
        "best_words": best["words"],
    }
    with open(logpath, "w") as fh:
        json.dump(record, fh, indent=2)

    print(f"n={args.n} s={args.s} group={args.group} wall={wall:.1f}s")
    print(f"best size : {best['size']}  (from {best['orbits_used']} orbits)")
    print(f"by seed   : {record['sizes_by_seed']}")
    print(f"log       : {logpath}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
