"""Search over generalised Varshamov-Tenengolts codes: weighted congruence classes.

Motivation. Plain MIS local search reaches 138 on the (s=1, n=11) cell where the
true optimum is 172, because that optimum is ALGEBRAIC -- it is the VT code
  VT_a(n) = { v : sum_i i*v_i == a  (mod n+1) }
and no amount of vertex swapping stumbles onto a global congruence condition.

So search the family directly. For a weight vector w in Z^n, modulus m, residue a:
  C(w, m, a) = { v in {0,1}^n : sum_i w_i v_i == a  (mod m) }
VT is the special case w = (1,2,...,n), m = n+1. Helberg-Ferreira is another
special case with a recurrence-defined w. Nobody appears to have swept the whole
family for s >= 2, where the best known codes came from general-purpose MIS
solvers rather than from structure.

Each candidate is scored by: take the largest residue class that is a valid
s-deletion code; if no class is valid, greedily extract the largest valid subset
of the best class. That last step is what makes this a hybrid rather than a pure
algebraic sweep.

Nothing here decides anything. Winners go to verify.py and verify2.py.

Run:
    python search/congruence.py --n 12 --s 1 --seconds 60          # should find 316
    python search/congruence.py --n 12 --s 2 --seconds 600 --workers 10
"""

from __future__ import annotations

import argparse
import json
import os
import random
import sys
import time
from multiprocessing import Pool

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from search.graph import load_graph  # noqa: E402

RESULTS_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "results")


class Family:
    """Precomputed structures for evaluating congruence classes fast."""

    def __init__(self, n: int, s: int):
        self.n = n
        self.s = s
        self.graph = load_graph(n, s)
        self.nv = self.graph.num_vertices
        self.words = self.graph.words

        # bit matrix: M[x, i] = symbol at position i of word x
        self.M = np.zeros((self.nv, n), dtype=np.int64)
        for x, w in enumerate(self.words):
            for i, ch in enumerate(w):
                if ch == "1":
                    self.M[x, i] = 1

        # deletion ball as clique ids, padded to a rectangle for vectorised checks
        clique_id = {t: i for i, t in enumerate(sorted(self.graph.cliques))}
        balls = []
        for x, w in enumerate(self.words):
            ids = sorted(clique_id[t] for t in _ball(w, s))
            balls.append(ids)
        self.maxball = max(len(b) for b in balls)
        self.ball_mat = np.full((self.nv, self.maxball), -1, dtype=np.int32)
        for x, ids in enumerate(balls):
            self.ball_mat[x, : len(ids)] = ids
        self.balls = balls

    # -- scoring -----------------------------------------------------------
    def class_is_valid(self, members: np.ndarray) -> bool:
        if members.size < 2:
            return True
        vals = self.ball_mat[members].ravel()
        vals = vals[vals >= 0]
        return np.unique(vals).size == vals.size

    def greedy_subset(self, members, rng) -> list[int]:
        """Largest valid subset found by one randomised greedy pass."""
        order = list(members)
        rng.shuffle(order)
        used: set[int] = set()
        out = []
        for v in order:
            ids = self.balls[v]
            if used.isdisjoint(ids):
                used.update(ids)
                out.append(v)
        return out

    def evaluate(self, w: np.ndarray, m: int, rng, greedy: bool = True):
        """Return (size, residue, members) for the best class of (w, m)."""
        res = (self.M @ w) % m
        best_size = 0
        best_members: list[int] = []
        # order residues by class size so the promising ones are tried first
        counts = np.bincount(res, minlength=m)
        for a in np.argsort(-counts):
            if counts[a] <= best_size:
                break
            members = np.flatnonzero(res == a)
            if self.class_is_valid(members):
                if members.size > best_size:
                    best_size = int(members.size)
                    best_members = members.tolist()
            elif greedy:
                sub = self.greedy_subset(members.tolist(), rng)
                if len(sub) > best_size:
                    best_size = len(sub)
                    best_members = sub
        return best_size, best_members


def _ball(word: str, s: int) -> set[str]:
    from itertools import combinations

    n = len(word)
    if s == 0:
        return {word}
    return {"".join(word[i] for i in keep) for keep in combinations(range(n), n - s)}


# ---------------------------------------------------------------------------


def hill_climb(fam: Family, rng: random.Random, m: int, seconds: float, w0=None):
    """Mutate one weight at a time, keep improvements. Returns (size, members, w)."""
    n = fam.n
    w = np.array(w0 if w0 is not None else [rng.randrange(m) for _ in range(n)], dtype=np.int64)
    best_size, best_members = fam.evaluate(w, m, rng)
    best_w = w.copy()
    t0 = time.time()
    stall = 0
    while time.time() - t0 < seconds:
        cand = best_w.copy()
        # mutate one or two coordinates
        for _ in range(1 if rng.random() < 0.8 else 2):
            cand[rng.randrange(n)] = rng.randrange(m)
        size, members = fam.evaluate(cand, m, rng)
        if size > best_size:
            best_size, best_members, best_w = size, members, cand
            stall = 0
        else:
            stall += 1
            if stall > 60:  # restart this modulus from a fresh weight vector
                cand = np.array([rng.randrange(m) for _ in range(n)], dtype=np.int64)
                size, members = fam.evaluate(cand, m, rng)
                if size > best_size:
                    best_size, best_members, best_w = size, members, cand
                stall = 0
    return best_size, best_members, best_w.tolist()


def _worker(args):
    n, s, seed, seconds, m_lo, m_hi = args
    rng = random.Random(seed)
    fam = Family(n, s)
    t0 = time.time()
    best = (0, [], None, None)
    # sweep moduli round-robin so no single m eats the budget
    moduli = list(range(m_lo, m_hi + 1))
    rng.shuffle(moduli)
    per_m = max(1.0, seconds / max(1, len(moduli)))
    for m in moduli:
        if time.time() - t0 >= seconds:
            break
        size, members, w = hill_climb(fam, rng, m, per_m)
        if size > best[0]:
            best = (size, members, w, m)
    return {
        "seed": seed,
        "size": best[0],
        "weights": best[2],
        "modulus": best[3],
        "words": [fam.words[v] for v in sorted(best[1])],
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, required=True)
    ap.add_argument("--s", type=int, required=True)
    ap.add_argument("--seconds", type=float, default=60)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--workers", type=int, default=1)
    ap.add_argument("--m-lo", type=int, default=None)
    ap.add_argument("--m-hi", type=int, default=None)
    ap.add_argument("--tag", type=str, default="congruence")
    args = ap.parse_args()

    # A valid class has size about 2^n / m, so choose the modulus window around
    # the size we are chasing rather than sweeping blindly.
    if args.m_lo is None or args.m_hi is None:
        target_lo = {1: 0.5, 2: 0.5, 3: 0.5}.get(args.s, 0.5)
        full = 1 << args.n
        args.m_lo = args.m_lo or max(2, int(full / (full * target_lo)) + 1)
        args.m_hi = args.m_hi or max(args.m_lo + 1, 220)

    load_graph(args.n, args.s)
    jobs = [
        (args.n, args.s, args.seed + i, args.seconds, args.m_lo, args.m_hi)
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
        "method": "generalised VT: weighted congruence classes + hill climbing on weights",
        "n": args.n,
        "s": args.s,
        "modulus_range": [args.m_lo, args.m_hi],
        "seconds_per_worker": args.seconds,
        "workers": args.workers,
        "base_seed": args.seed,
        "wall_seconds": round(wall, 2),
        "best_size": best["size"],
        "best_weights": best["weights"],
        "best_modulus": best["modulus"],
        "sizes_by_seed": {str(r["seed"]): r["size"] for r in results},
        "best_words": best["words"],
    }
    with open(logpath, "w") as fh:
        json.dump(record, fh, indent=2)

    print(f"n={args.n} s={args.s} m in [{args.m_lo},{args.m_hi}] wall={wall:.1f}s")
    print(f"best size : {best['size']}")
    print(f"weights   : {best['weights']}  mod {best['modulus']}")
    print(f"by seed   : {record['sizes_by_seed']}")
    print(f"log       : {logpath}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
