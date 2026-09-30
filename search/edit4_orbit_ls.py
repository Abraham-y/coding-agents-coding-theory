"""Weighted local search on the orbit graph of quaternary single-edit codes.

WHY THIS FILE EXISTS. The all-subgroup sweep in `search/edit4.py` solves each
orbit problem with an exact MILP guarded by `--exact-cap` (default 900 orbits).
At n = 6 that cap silently SKIPS 15 of the 60 (subgroup x reversal) cells, and
they are precisely the weakest-symmetry ones:

    |H| = 1, rev = False   (4096 orbits)   the unrestricted problem
    |H| = 1, rev = True    (2080 orbits)   reversal-invariant codes
    |H| = 2, rev = False   (9 cells)       a single alphabet involution
    |H| = 3, rev = False   (4 cells)       a single 3-cycle

The recorded "best invariant = 112, proven optimal in the A_4 class" therefore
says nothing about these classes. They are the interesting ones: a mild
symmetry constrains the code little while halving or thirding the search space,
which is the prescribed-automorphism rung CLAUDE.md puts right after brute
force. Where the MILP was pushed past the cap it also degrades visibly -- one
1021-orbit cell in results/edit4_sweep_n6_20260806T025537.json returns an
incumbent of 8 -- so "MILP incumbent" is not evidence about these cells either.

WHAT IT DOES. Max-WEIGHT independent set on the orbit graph by iterated local
search (ARW-style (omega,1)-swaps + forced-insertion perturbation), which does
not care about graph size the way the MILP does. Vertices are usable orbits
(orbits with no internal conflict), weights are orbit sizes, and orbits are
adjacent iff any member conflicts with any member of the other. By the orbit
reduction lemma (formalized sorry-free in lean/OrbitReduction, part_c) the
maximum H-invariant code size equals the maximum weight here -- so a solution
expands to a genuine code, and is verified as one before being reported.

VALIDATION GATE (rule 4 in spirit): --validate re-solves the A_4 cell, whose
optimum is independently certified at 112 (DRUP, drat-trim). If the local
search does not reach 112 there, no number it produces elsewhere is trusted.

Run:
    python search/edit4_orbit_ls.py --n 6 --validate
    python search/edit4_orbit_ls.py --n 6 --missing --seeds 5 --time-limit 300
"""

from __future__ import annotations

import argparse
import json
import os
import random
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from search.edit4 import (  # noqa: E402
    PUBLISHED,
    Edit4Graph,
    all_subgroups_s4,
    orbit_problem,
    orbits_of,
    verify_edit4,
    word_transforms,
)

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS_DIR = os.path.join(ROOT, "results")
BEST_DIR = os.path.join(ROOT, "best")


# ---------------------------------------------------------------------------
# instance construction
# ---------------------------------------------------------------------------


def build_instance(g: Edit4Graph, transforms):
    """(usable orbits, adjacency lists over orbit indices, weights, total orbits)."""
    orbits = orbits_of(g, transforms)
    usable, edges = orbit_problem(g, orbits)
    k = len(usable)
    adj: list[list[int]] = [[] for _ in range(k)]
    for a, b in edges:
        adj[a].append(b)
        adj[b].append(a)
    weights = [len(o) for o in usable]
    return usable, adj, weights, len(orbits)


def adjacency_masks(adj):
    """Bitmask per vertex, for O(1) adjacency tests inside (1,2)-swaps."""
    masks = []
    for lst in adj:
        m = 0
        for u in lst:
            m |= 1 << u
        masks.append(m)
    return masks


# ---------------------------------------------------------------------------
# max-weight independent set by iterated local search
# ---------------------------------------------------------------------------


class MWIS:
    """Incremental MWIS state.

    tight_w[v] = total weight of SELECTED neighbours of v
    tight_c[v] = number of SELECTED neighbours of v

    A non-selected v is an improving insertion exactly when tight_w[v] < w[v]:
    insert v, drop its selected neighbours, gain w[v] - tight_w[v] > 0. This one
    rule subsumes (0,1) insertions and all (omega,1) swaps.
    """

    def __init__(self, weights, adj):
        self.w = weights
        self.adj = adj
        self.k = len(weights)
        self.sel = bytearray(self.k)
        self.tight_w = [0] * self.k
        self.tight_c = [0] * self.k
        self.total = 0
        # candidate queue: vertices whose tightness dropped and so might now be
        # improving. Rescanning all k vertices instead of tracking these is what
        # made the first version of this file useless above ~1000 vertices.
        self.inq = bytearray(self.k)
        self.queue: list[int] = []

    def insert(self, v):
        self.sel[v] = 1
        self.total += self.w[v]
        wv = self.w[v]
        tw, tc = self.tight_w, self.tight_c
        for u in self.adj[v]:
            tw[u] += wv
            tc[u] += 1

    def remove(self, v):
        self.sel[v] = 0
        self.total -= self.w[v]
        wv = self.w[v]
        tw, tc = self.tight_w, self.tight_c
        for u in self.adj[v]:
            tw[u] -= wv
            tc[u] -= 1

    def push(self, v):
        if not self.inq[v]:
            self.inq[v] = 1
            self.queue.append(v)

    def force_insert(self, v):
        """Insert v, evicting selected neighbours regardless of gain.

        Evicting u lowers tight_w for every neighbour of u, so exactly those
        vertices can newly become improving; they go on the queue. Neighbours of
        v only get tighter and are therefore never worth re-examining.
        """
        if self.sel[v]:
            return
        sel, adj, inq, q = self.sel, self.adj, self.inq, self.queue
        for u in adj[v]:
            if sel[u]:
                self.remove(u)
                for x in adj[u]:
                    if not sel[x] and not inq[x]:
                        inq[x] = 1
                        q.append(x)
        self.insert(v)

    def resolve(self):
        """Drain the candidate queue, applying improving insertions."""
        sel, tw, w, adj, inq, q = (self.sel, self.tight_w, self.w,
                                   self.adj, self.inq, self.queue)
        while q:
            v = q.pop()
            inq[v] = 0
            if sel[v] or tw[v] >= w[v]:
                continue
            for u in adj[v]:
                if sel[u]:
                    self.remove(u)
                    for x in adj[u]:
                        if not sel[x] and not inq[x]:
                            inq[x] = 1
                            q.append(x)
            self.insert(v)

    def push_all(self, order):
        for v in order:
            self.push(v)

    def selection(self):
        return [v for v in range(self.k) if self.sel[v]]

    def load(self, chosen):
        for v in range(self.k):
            self.sel[v] = 0
            self.tight_w[v] = 0
            self.tight_c[v] = 0
            self.inq[v] = 0
        self.queue.clear()
        self.total = 0
        for v in chosen:
            self.insert(v)


class WeightedSearcher:
    """Weighted port of `search.ils.Searcher`, the engine that reached 111 on the
    raw n=6 graph.

    The (omega,1) insertion rule alone (insert v when it outweighs every selected
    neighbour) is NOT enough: it can never perform "remove one, add two", which
    is the move maximum-independent-set local search actually lives on. Measured
    on the trivial cell, rule-only search tops out at 92-99 where this move set
    reaches 111. The weighted analogues used here:

      1-improvement  : a free vertex (no selected neighbour) is always added.
      (1,1)-swap     : u in S out, a 1-tight neighbour v in, when w[v] > w[u].
      (1,2)-swap     : u out, two NON-ADJACENT 1-tight neighbours v,x in,
                       when w[v] + w[x] > w[u].
      (omega,1)      : insert v evicting all selected neighbours, when
                       w[v] > total weight of those neighbours.
      plateau move   : u out, one 1-tight v in with w[v] == w[u]. Weight
                       unchanged; this is what crosses plateaus.
    """

    PAIR_CANDS = 30  # cap on 1-tight candidates examined per (1,2)-swap

    def __init__(self, weights, adj, masks, seed: int):
        self.w = weights
        self.nbrs = adj
        self.mask = masks
        self.k = len(weights)
        self.rng = random.Random(seed)
        self.in_S = bytearray(self.k)
        self.tight = [0] * self.k
        self.tight_w = [0] * self.k
        self.total = 0
        self.S: set[int] = set()
        self.free: set[int] = set(range(self.k))
        # (omega,1) can only fire at v when 0 < tight_w[v] < w[v], so it needs
        # w[v] >= 2. On an unweighted cell this list is empty and the rule is
        # skipped entirely -- scanning all k vertices for it was costing more
        # than the move set it was serving.
        self.heavy = [v for v in range(self.k) if weights[v] >= 2]

    def add(self, v):
        self.in_S[v] = 1
        self.total += self.w[v]
        self.S.add(v)
        self.free.discard(v)
        wv, t, tw, free = self.w[v], self.tight, self.tight_w, self.free
        for u in self.nbrs[v]:
            t[u] += 1
            tw[u] += wv
            if t[u] == 1:
                free.discard(u)

    def remove(self, v):
        self.in_S[v] = 0
        self.total -= self.w[v]
        self.S.discard(v)
        wv, t, tw, free = self.w[v], self.tight, self.tight_w, self.free
        if t[v] == 0:
            free.add(v)
        for u in self.nbrs[v]:
            t[u] -= 1
            tw[u] -= wv
            if t[u] == 0 and not self.in_S[u]:
                free.add(u)

    def clear(self):
        for v in range(self.k):
            self.in_S[v] = 0
            self.tight[v] = 0
            self.tight_w[v] = 0
        self.total = 0
        self.S = set()
        self.free = set(range(self.k))

    def maximize(self):
        """Add free vertices, heaviest-of-a-sample first, until none remain."""
        rng, w = self.rng, self.w
        while self.free:
            if len(self.free) <= 8:
                v = max(self.free, key=lambda x: (w[x], rng.random()))
            else:
                sample = rng.sample(tuple(self.free), 8)
                v = max(sample, key=lambda x: w[x])
            self.add(v)

    def improve(self) -> bool:
        """One pass of (1,1)-, (1,2)- and (omega,1)-improvements."""
        improved = False
        w, t, tw, in_S, mask = self.w, self.tight, self.tight_w, self.in_S, self.mask
        members = list(self.S)
        self.rng.shuffle(members)
        for u in members:
            if not in_S[u]:
                continue
            cands = [v for v in self.nbrs[u] if not in_S[v] and t[v] == 1]
            if not cands:
                continue
            cands.sort(key=lambda v: -w[v])
            # (1,1)
            if w[cands[0]] > w[u]:
                self.remove(u)
                self.add(cands[0])
                improved = True
                continue
            # (1,2): heaviest non-adjacent pair beating w[u]
            top = cands[:self.PAIR_CANDS]
            found = None
            for i in range(len(top)):
                vi = top[i]
                mi = mask[vi]
                for j in range(i + 1, len(top)):
                    vj = top[j]
                    if w[vi] + w[vj] <= w[u]:
                        break  # sorted desc: no later j can help
                    if not (mi >> vj) & 1:
                        found = (vi, vj)
                        break
                if found:
                    break
            if found:
                self.remove(u)
                self.add(found[0])
                self.add(found[1])
                improved = True
        # (omega,1): any non-selected vertex outweighing its selected neighbours
        for v in self.heavy:
            if not in_S[v] and tw[v] < w[v]:
                for u in self.nbrs[v]:
                    if in_S[u]:
                        self.remove(u)
                self.add(v)
                improved = True
        return improved

    def plateau(self, moves: int):
        w, t, in_S = self.w, self.tight, self.in_S
        for _ in range(moves):
            if not self.S:
                return
            u = self.rng.choice(tuple(self.S))
            cands = [v for v in self.nbrs[u]
                     if not in_S[v] and t[v] == 1 and w[v] == w[u]]
            if not cands:
                continue
            self.remove(u)
            self.add(self.rng.choice(cands))
            self.maximize()

    def local_search(self):
        while True:
            self.maximize()
            if not self.improve():
                break

    def perturb(self, strength: int):
        for _ in range(strength):
            v = self.rng.randrange(self.k)
            if self.in_S[v]:
                continue
            for u in self.nbrs[v]:
                if self.in_S[u]:
                    self.remove(u)
            self.add(v)

    def run(self, time_limit: float):
        t0 = time.time()
        self.clear()
        self.local_search()
        best, best_w = sorted(self.S), self.total
        iters, stall = 0, 0
        while time.time() - t0 < time_limit:
            iters += 1
            self.perturb(1 + self.rng.randrange(3))
            self.local_search()
            self.plateau(self.rng.randrange(8))
            self.local_search()
            if self.total > best_w:
                best, best_w, stall = sorted(self.S), self.total, 0
            else:
                stall += 1
                if stall > 400:
                    self.clear()
                    for v in best:
                        self.add(v)
                    stall = 0
        return best, best_w, iters


def ils(weights, adj, seed: int, time_limit: float, perturb: int = 3):
    """Iterated local search. Returns (best_selection, best_weight, iterations)."""
    masks = adjacency_masks(adj)
    return WeightedSearcher(weights, adj, masks, seed).run(time_limit)


def ils_rulesonly(weights, adj, seed: int, time_limit: float, perturb: int = 3):
    """The original (omega,1)-only search, kept to document why it is not enough."""
    rng = random.Random(seed)
    k = len(weights)
    st = MWIS(weights, adj)
    order = list(range(k))
    # greedy start: heavy, low-degree first, randomised
    order.sort(key=lambda v: (-weights[v] / (1.0 + len(adj[v])), rng.random()))
    st.push_all(order)
    st.resolve()
    best = st.selection()
    best_w = st.total

    t_end = time.time() + time_limit
    iters = 0
    stall = 0
    randrange = rng.randrange
    # check the clock every so often rather than every iteration -- with the
    # incremental queue an iteration is microseconds and time() would dominate
    CHECK = 256
    while True:
        for _ in range(CHECK):
            iters += 1
            for _ in range(perturb):
                st.force_insert(randrange(k))
            st.resolve()
            if st.total > best_w:
                best_w = st.total
                best = st.selection()
                stall = 0
            else:
                stall += 1
                if stall >= 2000:
                    st.load(best)
                    stall = 0
        if time.time() >= t_end:
            break
    return best, best_w, iters


# ---------------------------------------------------------------------------
# independent check: pairwise Levenshtein >= 3, no shared code with edit4.ball
# ---------------------------------------------------------------------------


def lev(a: str, b: str) -> int:
    prev = list(range(len(b) + 1))
    for i, ca in enumerate(a, 1):
        cur = [i]
        for j, cb in enumerate(b, 1):
            cur.append(min(prev[j] + 1, cur[j - 1] + 1, prev[j - 1] + (ca != cb)))
        prev = cur
    return prev[-1]


def independent_check(words):
    """Valid single-edit code <=> pairwise edit distance >= 3. Returns bad pairs."""
    bad = 0
    for i in range(len(words)):
        for j in range(i + 1, len(words)):
            if lev(words[i], words[j]) < 3:
                bad += 1
    return bad


# ---------------------------------------------------------------------------
# cells
# ---------------------------------------------------------------------------


def cell_label(sub, rev, orbits):
    return f"|H|={len(sub)}{'+rev' if rev else ''} ({orbits} orbits)"


def _cell_work(job):
    """Worker for the parallel driver. Module-level: macOS spawn must import it."""
    n, gens, rev, seed, seconds = job
    g = Edit4Graph(n)
    sub = [tuple(p) for p in gens]
    usable, adj, weights, n_orbits = build_instance(g, word_transforms(sub, n, rev))
    sel, wt, iters = ils(weights, adj, seed=seed, time_limit=seconds)
    words = sorted(g.words[v] for i in sel for v in usable[i])
    return {"seed": seed, "weight": wt, "iters": iters, "words": words,
            "orbits_total": n_orbits, "orbits_usable": len(usable),
            "orbit_edges": sum(len(a) for a in adj) // 2}


def run_cell_parallel(n, sub, rev, seeds, seconds, workers):
    """Best over `seeds` independent runs, `workers` at a time."""
    from multiprocessing import Pool
    gens = sorted(tuple(p) for p in sub)
    jobs = [(n, gens, rev, s, seconds) for s in seeds]
    with Pool(min(workers, len(jobs))) as pool:
        out = pool.map(_cell_work, jobs)
    out.sort(key=lambda r: -r["weight"])
    return out


def missing_cells(exact_cap: int = 900):
    """The (subgroup, reversal) cells the MILP sweep skipped, weakest first."""
    out = []
    for sub in all_subgroups_s4():
        for rev in (False, True):
            # the sweep skipped a cell iff its orbit count exceeded the cap;
            # orbit count is recomputed exactly in run_cell, this is the filter
            if len(sub) <= 3 and not rev:
                out.append((sub, rev))
            elif len(sub) == 1 and rev:
                out.append((sub, rev))
    return out


def run_cell(g, sub, rev, seeds, time_limit, verbose=True):
    transforms = word_transforms(sub, g.n, rev)
    t0 = time.time()
    usable, adj, weights, n_orbits = build_instance(g, transforms)
    build_s = time.time() - t0
    n_edges = sum(len(a) for a in adj) // 2

    best_w, best_sel, best_seed = 0, None, None
    for s in seeds:
        sel, w, iters = ils(weights, adj, seed=s, time_limit=time_limit)
        if verbose:
            print(f"      seed {s:<4} weight {w:>4}  ({iters} ILS iterations)")
        if w > best_w:
            best_w, best_sel, best_seed = w, sel, s

    words = sorted(g.words[v] for i in best_sel for v in usable[i])
    ok, size, reason = verify_edit4(words, g.n)
    return {
        "order": len(sub),
        "reversal": rev,
        "orbits_total": n_orbits,
        "orbits_usable": len(usable),
        "orbit_edges": n_edges,
        "build_seconds": round(build_s, 1),
        "size": size if ok else 0,
        "weight": best_w,
        "best_seed": best_seed,
        "verified": ok,
        "reason": reason,
        "generators": sorted(tuple(p) for p in sub),
    }, words, ok, reason


def conjugacy_representatives(subs):
    """One representative per conjugacy class among the never-attempted cells.

    If H' = s H s^-1 then applying s symbol-wise is an automorphism of the edit
    graph carrying H-invariant codes bijectively onto H'-invariant ones, so the
    two cells have equal optima. The 15 skipped cells therefore hold only 4
    distinct problems. (Corroborated by the censuses: all 6 transposition cells
    build 2080/1359/96182 orbit graphs, all 3 double-transposition cells build
    2048/2046/195224.)
    """
    def order2(perm):
        return next(s for s in subs if len(s) == 2 and tuple(perm) in s)

    triv = next(s for s in subs if len(s) == 1)
    return [
        ("reversal only",                    triv,                     True),
        ("transposition (0 1)",              order2((1, 0, 2, 3)),     False),
        ("double transposition (0 1)(2 3)",  order2((1, 0, 3, 2)),     False),
        ("3-cycle (0 1 2)",
         next(s for s in subs if len(s) == 3 and (1, 2, 0, 3) in s),   False),
    ]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=6)
    ap.add_argument("--batch", action="store_true",
                    help="calibration cell + the 4 distinct never-attempted cells")
    ap.add_argument("--workers", type=int, default=5)
    ap.add_argument("--seconds", type=float, default=900.0)
    ap.add_argument("--validate", action="store_true",
                    help="re-solve the A_4 cell, whose optimum is certified 112")
    ap.add_argument("--missing", action="store_true",
                    help="run the cells the MILP sweep never attempted")
    ap.add_argument("--seeds", type=int, default=3)
    ap.add_argument("--time-limit", type=float, default=120.0)
    args = ap.parse_args()

    t0 = time.time()
    g = Edit4Graph(args.n)
    print(f"Edit4 G({args.n}): {g.num_vertices} vertices, {len(g.cliques)} cliques "
          f"({time.time() - t0:.1f}s)")
    seeds = list(range(1, args.seeds + 1))

    if args.validate:
        # A_4 = the 12 even permutations; certified optimum 112 without reversal
        a4 = next(s for s in all_subgroups_s4() if len(s) == 12)
        print(f"\nVALIDATION GATE: A_4 cell, certified optimum 112")
        rec, words, ok, reason = run_cell(g, a4, False, seeds, args.time_limit)
        print(f"   orbits {rec['orbits_total']} ({rec['orbits_usable']} usable), "
              f"{rec['orbit_edges']} orbit edges")
        print(f"   local search best : {rec['size']}   verify_edit4 = {ok} {reason}")
        print(f"   certified optimum : 112")
        if not ok:
            print("   FAIL: produced set is not a valid code")
            return 1
        if rec["size"] != 112:
            print(f"   FAIL: expected 112, got {rec['size']} -- do not trust other cells")
            return 1
        bad = independent_check(words)
        print(f"   independent Levenshtein check: {bad} violating pairs")
        print("   GATE PASSED" if bad == 0 else "   FAIL: independent check disagrees")
        return 0 if bad == 0 else 1

    if args.batch:
        subs = all_subgroups_s4()
        triv = next(s for s in subs if len(s) == 1)
        seeds = list(range(1, args.seeds + 1))
        cells = [("CALIBRATION: trivial (= unrestricted problem)", triv, False)]
        cells += conjugacy_representatives(subs)

        report, calib = [], None
        for name, sub, rev in cells:
            print(f"\n=== {name} ===", flush=True)
            t = time.time()
            out = run_cell_parallel(args.n, sub, rev, seeds, args.seconds, args.workers)
            best = out[0]
            words = sorted(best["words"])
            ok, size, reason = verify_edit4(words, args.n)
            print(f"  {best['orbits_total']} orbits ({best['orbits_usable']} usable), "
                  f"{best['orbit_edges']} edges, wall {time.time() - t:.0f}s", flush=True)
            print(f"  by seed: {[r['weight'] for r in out]}", flush=True)
            print(f"  best {best['weight']}  verified={ok} {reason}", flush=True)
            rec = {"cell": name, "order": len(sub), "reversal": rev,
                   "orbits_total": best["orbits_total"],
                   "orbits_usable": best["orbits_usable"],
                   "orbit_edges": best["orbit_edges"],
                   "by_seed": {str(r["seed"]): r["weight"] for r in out},
                   "best": best["weight"], "verified": ok}
            report.append(rec)
            if calib is None:
                calib = best["weight"]
                print(f"  bar for calibration: 111 (repo ILS) / 113 (KaMIS)", flush=True)
                if calib < 110:
                    print(f"  GATE FAILED at {calib}: solver still under-powered; "
                          f"downstream cells are not interpretable", flush=True)
            if ok and best["weight"] >= 113:
                bad = independent_check(words)
                print(f"  >= 113 -- independent Levenshtein: {bad} violating pairs",
                      flush=True)
                if bad == 0:
                    os.makedirs(BEST_DIR, exist_ok=True)
                    ap_ = os.path.join(BEST_DIR,
                                       f"edit4_n{args.n}_orbit_size{size}.txt")
                    with open(ap_, "w") as fh:
                        fh.write(f"# quaternary single-edit code, n = {args.n}\n")
                        fh.write(f"# size = {size}, invariant under {name}\n")
                        for w in words:
                            fh.write(w + "\n")
                    print(f"  saved {ap_} -- rule 6: bug until re-verified fresh",
                          flush=True)

        os.makedirs(RESULTS_DIR, exist_ok=True)
        stamp = time.strftime("%Y%m%dT%H%M%S")
        path = os.path.join(RESULTS_DIR, f"edit4_orbit_batch_n{args.n}_{stamp}.json")
        with open(path, "w") as fh:
            json.dump({"n": args.n, "seeds": seeds, "seconds": args.seconds,
                       "workers": args.workers, "calibration_best": calib,
                       "calibration_bar": 111, "cells": report,
                       "published_best": PUBLISHED.get(args.n, {}).get("best"),
                       "literature_best": 114,
                       "elapsed_seconds": round(time.time() - t0, 1)}, fh, indent=2)
        print(f"\nlog: {path}")
        return 0

    if args.missing:
        cells = missing_cells()
        pub = PUBLISHED.get(args.n, {})
        print(f"\n{len(cells)} cells never attempted by the MILP sweep "
              f"(exact_cap=900). W-H best known {pub.get('best')}, "
              f"literature best 114.\n")
        records, best = [], (0, None, None)
        for sub, rev in cells:
            transforms = word_transforms(sub, args.n, rev)
            label = f"|H|={len(sub)}{'+rev' if rev else ''}"
            print(f"  {label}  gens={sorted(tuple(p) for p in sub)}")
            rec, words, ok, reason = run_cell(g, sub, rev, seeds, args.time_limit)
            print(f"      {rec['orbits_total']} orbits ({rec['orbits_usable']} usable), "
                  f"{rec['orbit_edges']} edges, build {rec['build_seconds']}s")
            print(f"      -> invariant code {rec['size']}   verified={ok} {reason}")
            if not ok:
                print(f"      VERIFIER REJECTED -- bug, investigate")
            records.append(rec)
            if ok and rec["size"] > best[0]:
                best = (rec["size"], words, label)
            print()

        size, words, label = best
        print(f"best over never-attempted cells : {size}  ({label})")
        print(f"W-H Table 3 best known          : {pub.get('best')}")
        print(f"literature best known           : 114 (Ashlock-Houghten 2012)")

        os.makedirs(RESULTS_DIR, exist_ok=True)
        stamp = time.strftime("%Y%m%dT%H%M%S")
        path = os.path.join(RESULTS_DIR, f"edit4_orbit_ls_n{args.n}_{stamp}.json")
        with open(path, "w") as fh:
            json.dump({"n": args.n, "seeds": seeds, "time_limit": args.time_limit,
                       "cells": records, "best_size": size, "best_group": label,
                       "published_best": pub.get("best"), "literature_best": 114,
                       "elapsed_seconds": round(time.time() - t0, 1)}, fh, indent=2)
        print(f"log: {path}")

        if words and size >= 113:
            bad = independent_check(words)
            print(f"\nsize {size} >= 113 -- independent Levenshtein check: "
                  f"{bad} violating pairs")
            if bad == 0:
                os.makedirs(BEST_DIR, exist_ok=True)
                apath = os.path.join(BEST_DIR, f"edit4_n{args.n}_orbit_size{size}.txt")
                with open(apath, "w") as fh:
                    fh.write(f"# quaternary single-edit code, n = {args.n}\n")
                    fh.write(f"# size = {size}, invariant under {label}\n")
                    for w in words:
                        fh.write(w + "\n")
                print(f"saved: {apath} -- CLAUDE.md rule 6: treat as a bug until "
                      f"re-verified in a fresh process")
        return 0

    ap.error("pass --validate or --missing")


if __name__ == "__main__":
    raise SystemExit(main())
