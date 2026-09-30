"""alpha_k with the MNSX forced-action rules. Third independent solver.

maj3_alpha.py  keys states on sets of consistent hard inputs   -> dies at k=3
maj3_struct.py keys states on canonical recursive tree objects -> k=3 in 222110 states
maj3_fast.py   adds the two forced-action rules                -> target: k=4

The rules, from MNSX section 3.2, both restated and re-derived here rather than
transcribed on faith:

RULE A (stop when the root value is pinned).
  If the configuration is inconsistent with every 1-hard input, then flipping any
  UNQUERIED leaf yields an input still consistent with the configuration, hence
  still of value 0. So no unqueried leaf is sensitive, every remaining query has
  immediate contribution -alpha*Pr[minority] <= 0, and stopping is optimal.
  Test: c_1(state) == 0  =>  V(state) = 0.

RULE B (force queries that cannot be the absolute minority).
  On-path nodes alternate in value: the minority path has value (depth mod 2) at
  each depth, since the root is 0 and each step moves to the minority child.
  Therefore, for a node whose value is DETERMINED by the configuration:
    - value != depth mod 2  =>  that node is off the minority path, so every leaf
      in its subtree is off-path;
    - value == depth mod 2  =>  it is the unique sibling that can be on the path,
      so both its SIBLINGS' subtrees are off-path.
  An off-path leaf has minority-count 0, so querying it contributes
  2^-k * Pr[sensitive] >= 0 and strictly more information. Querying it first is
  without loss of generality, so the max over actions collapses to a single
  forced action -- branching factor 1 instead of (number of unqueried leaves).

Rule B rests on an exchange argument that this file does not prove; it is
asserted by MNSX. It is therefore treated as a HYPOTHESIS and tested the only way
that matters: the solver must still reproduce alpha_1 = 2, alpha_2 = 24/7 and
alpha_3 = 12231/2203 exactly. If the rule were wrong, those would move.

Exact rational arithmetic throughout (rules 10, 12).

Run:
    python search/maj3_fast.py --k 3     # must give 12231/2203
    python search/maj3_fast.py --k 4     # target 2027349/216164
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time
from fractions import Fraction

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from search.maj3_struct import (  # noqa: E402
    ONE,
    PUBLISHED_ALPHA,
    U,
    ZERO,
    canonical,
    clear_caches,
    counts,
    fresh_state,
    leaf_stats,
    leaf_value,
    set_leaf,
)

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS_DIR = os.path.join(ROOT, "results")


# ---------------------------------------------------------------------------
# forced-action detection
# ---------------------------------------------------------------------------


def _determined(state):
    """0 or 1 if the subtree's value is pinned by the configuration, else None."""
    c0, c1 = counts(state)
    if c0 == 0:
        return 1
    if c1 == 0:
        return 0
    return None


def _first_unqueried(state, h: int, offset: int):
    if h == 0:
        return offset if state == U else None
    w = 3 ** (h - 1)
    for i in range(3):
        got = _first_unqueried(state[i], h - 1, offset + i * w)
        if got is not None:
            return got
    return None


def find_forced(state, h: int, depth: int = 0, offset: int = 0):
    """Return one unqueried leaf provably off the minority path, or None."""
    det = _determined(state)
    onpath = depth % 2
    if det is not None and det != onpath:
        return _first_unqueried(state, h, offset)

    if h == 0:
        return None

    w = 3 ** (h - 1)
    child_onpath = (depth + 1) % 2
    dets = [_determined(c) for c in state]
    special = [i for i in range(3) if dets[i] is not None and dets[i] == child_onpath]
    if special:
        i0 = special[0]
        for j in range(3):
            if j == i0:
                continue
            got = _first_unqueried(state[j], h - 1, offset + j * w)
            if got is not None:
                return got

    for i in range(3):
        got = find_forced(state[i], h - 1, depth + 1, offset + i * w)
        if got is not None:
            return got
    return None


# ---------------------------------------------------------------------------


class FastSolver:
    def __init__(self, k: int, progress: int = 0, forced_cap: int = 1_500_000):
        self.k = k
        self.nleaves = 3 ** k
        self.weight = Fraction(1, 2) ** k
        self.states = 0
        self.forced_hits = 0
        self.progress = progress
        # About 99% of visited states are FORCED (single legal action). Those are
        # exactly the "unstable" configurations MNSX collapse away, and they are
        # what exhausts memory. They still deserve SOME caching -- forced chains
        # are reachable from several parents -- so they live in a separate table
        # that is dropped wholesale when it gets too big. Stable states, which
        # are the genuinely reusable ones, stay in an uncapped table.
        self.forced_cap = forced_cap
        self.forced_flushes = 0

    def solve(self, alpha: Fraction):
        self.cache = {}
        self.forced_cache = {}
        self.states = 0
        self.forced_hits = 0
        self.forced_flushes = 0
        # module caches are keyed on states from the PREVIOUS iteration too;
        # drop them so memory does not accumulate across the fixed point.
        clear_caches()
        return self._value(canonical(fresh_state(self.k)), alpha, root=True)

    def _query(self, state, v: int, alpha: Fraction, c0: int, sens, minr):
        """Value, Pq, Pm of querying leaf v at `state`."""
        p_sens = Fraction(sens[v], c0)
        p_min = Fraction(minr[v], c0)
        total = self.weight * p_sens - alpha * p_min
        pq = p_sens
        pm = p_min
        for a in (ZERO, ONE):
            # Canonicalise IMMEDIATELY. Otherwise every raw permutation of a
            # configuration gets its own entry in the counts/canon caches, and
            # raw states vastly outnumber canonical ones -- that is what
            # exhausted memory on the first k=4 attempt.
            child = canonical(set_leaf(state, self.k, v, a))
            c0a = counts(child)[0]
            if c0a == 0:
                continue
            w = Fraction(c0a, c0)
            sub_val, sub_pq, sub_pm = self._value(child, alpha, False)
            total += w * sub_val
            pq += w * sub_pq
            pm += w * sub_pm
        return total, pq, pm

    def _value(self, state, alpha: Fraction, root: bool = False):
        # `state` is already canonical: _query canonicalises before recursing.
        key = (state, root)
        cached = self.cache.get(key)
        if cached is not None:
            return cached
        cached = self.forced_cache.get(key)
        if cached is not None:
            return cached

        zero = (Fraction(0), Fraction(0), Fraction(0))
        c0, c1 = counts(state)
        if c0 == 0:
            self.cache[key] = zero
            return zero

        # RULE A
        if c1 == 0 and not root:
            self.cache[key] = zero
            return zero

        self.states += 1
        if self.progress and self.states % self.progress == 0:
            print(f"    ... {self.states} states, {self.forced_hits} forced, "
                  f"{len(self.cache)} stable cached, {self.forced_flushes} flushes",
                  flush=True)
        # cache=False: this top-level array is never reused (each canonical state
        # is visited once), and caching it at k=4 wants tens of GB. Sub-states
        # inside the recursion are still cached, which is where the reuse is.
        sens, minr = leaf_stats(state, self.k, 0, cache=False)

        # RULE B
        forced = find_forced(state, self.k)
        if forced is not None:
            self.forced_hits += 1
            out = self._query(state, forced, alpha, c0, sens, minr)
            if len(self.forced_cache) >= self.forced_cap:
                self.forced_cache.clear()
                self.forced_flushes += 1
            self.forced_cache[key] = out
            return out

        best = None if root else zero
        for v in range(self.nleaves):
            if leaf_value(state, self.k, v) != U:
                continue
            cand = self._query(state, v, alpha, c0, sens, minr)
            if best is None or cand[0] > best[0]:
                best = cand
        if best is None:
            best = zero
        self.cache[key] = best
        return best

    def find_alpha(self, alpha0: Fraction = Fraction(0), max_iter: int = 40, verbose=False):
        alpha = alpha0
        history = []
        for it in range(max_iter):
            t0 = time.time()
            value, pq, pm = self.solve(alpha)
            if pm == 0:
                raise RuntimeError("optimal tree never queries the minority")
            new_alpha = self.weight * pq / pm
            dt = time.time() - t0
            history.append({
                "iter": it, "alpha_in": str(alpha), "rho": str(value),
                "Pq": str(pq), "Pm": str(pm), "alpha_out": str(new_alpha),
                "states": self.states, "forced": self.forced_hits,
                "seconds": round(dt, 1),
            })
            if verbose:
                print(f"  iter {it}: alpha={alpha}")
                print(f"           rho={value}")
                print(f"           Pq={pq}  Pm={pm}")
                print(f"           -> alpha={new_alpha}")
                print(f"           {self.states} states ({self.forced_hits} forced), {dt:.1f}s")
            if new_alpha == alpha:
                return alpha, history
            alpha = new_alpha
        raise RuntimeError("alpha iteration did not converge")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--k", type=int, required=True)
    ap.add_argument("--quiet", action="store_true")
    ap.add_argument("--progress", type=int, default=0,
                    help="print a line every N states (0 = off)")
    ap.add_argument("--forced-cap", type=int, default=1_500_000,
                    help="max entries in the forced-state cache before it is dropped")
    args = ap.parse_args()

    t0 = time.time()
    solver = FastSolver(args.k, progress=args.progress, forced_cap=args.forced_cap)
    print(f"k={args.k}: {solver.nleaves} leaves, forced-action rules ON")

    alpha, history = solver.find_alpha(verbose=not args.quiet)
    elapsed = time.time() - t0

    published = PUBLISHED_ALPHA.get(args.k)
    match = published is not None and alpha == published
    print(f"alpha_{args.k} = {alpha}")
    if published is not None:
        print(f"published  = {published}  -> {'MATCH' if match else 'MISMATCH'}")
    else:
        print("published  = (none -- this value is NOT in the literature)")
    base = 2 + float(alpha) ** (-1.0 / args.k)
    print(f"lower-bound base = {base:.9f}   [float, display only]")
    print(f"elapsed {elapsed:.1f}s, {solver.states} states")

    os.makedirs(RESULTS_DIR, exist_ok=True)
    path = os.path.join(RESULTS_DIR, f"maj3_fast_k{args.k}.json")
    with open(path, "w") as fh:
        json.dump({
            "k": args.k, "alpha_str": str(alpha),
            "alpha_numerator": alpha.numerator, "alpha_denominator": alpha.denominator,
            "published": str(published) if published else None,
            "matches_published": match,
            "base_float_display_only": base,
            "states": solver.states,
            "elapsed_seconds": round(elapsed, 2),
            "iterations": history,
        }, fh, indent=2)
    print(f"log: {path}")
    return 0 if (published is None or match) else 1


if __name__ == "__main__":
    raise SystemExit(main())
