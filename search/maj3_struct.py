"""Structured DP for alpha_k: configurations as recursive tree states, not input sets.

The brute-force DP in maj3_alpha.py keys each state on the SET of consistent hard
inputs. That dies at k=3, where |H_3^0| = 3^13 = 1,594,323. This module never
enumerates an input. It represents a configuration as a recursive object

    state ::= 'U' | '0' | '1' | (state, state, state)

('U' = unqueried leaf) and computes every quantity it needs by counting
recursions over that object.

Structure of a b-hard input of height h >= 1: exactly one child is the MINORITY
and carries value 1-b; the other two are majority and carry value b. Hence

    c_b(node)        = sum_i  c_{1-b}(ch_i) * prod_{j != i} c_b(ch_j)
    sens_b(node, v)  = sens_b(ch_i, v) * sum_{j != i} c_{1-b}(ch_j) * prod_{l not in {i,j}} c_b(ch_l)
    min_b(node, v)   = min_{1-b}(ch_i, v) * prod_{j != i} c_b(ch_j)          [v inside child i]

with the base case h = 0: a lone leaf is sensitive in every consistent input (its
flip flips the value) and is trivially the end of the minority path, so
sens_b = min_b = c_b.

Values are invariant under permuting the three children of any node, so states
are memoised on a canonical form obtained by recursively sorting children. That
is the symmetry reduction; the forced-action ("unstable configuration") rules of
MNSX section 3.2 are NOT hard-coded -- the DP rediscovers them, which is slower
but removes a whole class of transcription error.

Exact rational arithmetic throughout (rules 10 and 12).

Run:
    python search/maj3_struct.py --k 2          # must give 24/7
    python search/maj3_struct.py --k 3          # target 12231/2203
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time
from fractions import Fraction

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS_DIR = os.path.join(ROOT, "results")

PUBLISHED_ALPHA = {
    1: Fraction(2),
    2: Fraction(24, 7),
    3: Fraction(12231, 2203),
    4: Fraction(2027349, 216164),
}

U, ZERO, ONE = "U", "0", "1"


# ---------------------------------------------------------------------------
# state construction and canonicalisation
# ---------------------------------------------------------------------------


def fresh_state(h: int):
    if h == 0:
        return U
    sub = fresh_state(h - 1)
    return (sub, sub, sub)


def num_leaves(h: int) -> int:
    return 3 ** h


def set_leaf(state, h: int, v: int, value: str):
    """Return a new state with leaf index v set to `value`."""
    if h == 0:
        return value
    w = 3 ** (h - 1)
    i, rest = divmod(v, w)
    children = list(state)
    children[i] = set_leaf(children[i], h - 1, rest, value)
    return tuple(children)


def leaf_value(state, h: int, v: int):
    if h == 0:
        return state
    w = 3 ** (h - 1)
    i, rest = divmod(v, w)
    return leaf_value(state[i], h - 1, rest)


_canon_cache: dict = {}


def canonical(state):
    """Recursively sort children. Values are invariant under this.

    Sorting uses the natural tuple order, NOT key=repr. All three children of a
    node have the same height, so they are all strings or all tuples and compare
    directly; building a repr string per comparison was pure overhead.
    """
    if not isinstance(state, tuple):
        return state
    cached = _canon_cache.get(state)
    if cached is not None:
        return cached
    out = tuple(sorted(canonical(c) for c in state))
    _canon_cache[state] = out
    return out


def clear_caches():
    """Drop the module-level memo tables. Called between alpha iterations so the
    raw-state keyed caches do not accumulate across the whole fixed-point run."""
    _canon_cache.clear()
    _counts_cache.clear()
    _stats_cache.clear()


# ---------------------------------------------------------------------------
# counting recursions
# ---------------------------------------------------------------------------


_counts_cache: dict = {}


def counts(state):
    """(c_0, c_1): number of 0-hard and 1-hard inputs consistent with `state`."""
    cached = _counts_cache.get(state)
    if cached is not None:
        return cached
    if state == U:
        out = (1, 1)
    elif state == ZERO:
        out = (1, 0)
    elif state == ONE:
        out = (0, 1)
    else:
        ch = [counts(c) for c in state]
        c0 = 0
        c1 = 0
        for i in range(3):
            others = [j for j in range(3) if j != i]
            c0 += ch[i][1] * ch[others[0]][0] * ch[others[1]][0]
            c1 += ch[i][0] * ch[others[0]][1] * ch[others[1]][1]
        out = (c0, c1)
    _counts_cache[state] = out
    return out


_stats_cache: dict = {}


def leaf_stats(state, h: int, b: int, cache: bool = True):
    """Arrays over the 3^h leaves: (sens[v], minr[v]) counts for b-hard inputs.

    `cache=False` is for the TOP-level call only. The DP visits each canonical
    state once, so caching a 2 x 3^k array for every top-level state is pure
    waste -- at k=4 that alone wants ~32 GB. Sub-states are heavily shared and
    are always cached, which is where the real speedup lives.
    """
    key = (state, b)
    if cache:
        cached = _stats_cache.get(key)
        if cached is not None:
            return cached

    if h == 0:
        cb = counts(state)[b]
        out = ([cb], [cb])
    else:
        ch_counts = [counts(c) for c in state]
        sens = []
        minr = []
        for i in range(3):
            others = [j for j in range(3) if j != i]
            # multiplier for sensitivity: child i must be a MAJORITY child
            sens_mult = 0
            for j in others:
                l = [x for x in others if x != j][0]
                sens_mult += ch_counts[j][1 - b] * ch_counts[l][b]
            # multiplier for minority: child i must be THE minority child
            min_mult = ch_counts[others[0]][b] * ch_counts[others[1]][b]

            sub_sens_b, _ = leaf_stats(state[i], h - 1, b)
            _, sub_min_nb = leaf_stats(state[i], h - 1, 1 - b)
            for v in range(3 ** (h - 1)):
                sens.append(sub_sens_b[v] * sens_mult)
                minr.append(sub_min_nb[v] * min_mult)
        out = (sens, minr)

    _stats_cache[key] = out
    return out


# ---------------------------------------------------------------------------
# the DP
# ---------------------------------------------------------------------------


class StructSolver:
    def __init__(self, k: int):
        self.k = k
        self.nleaves = 3 ** k
        self.weight = Fraction(1, 2) ** k
        self.states = 0

    def solve(self, alpha: Fraction):
        self.cache = {}
        self.states = 0
        return self._value(fresh_state(self.k), alpha, root=True)

    def _value(self, state, alpha: Fraction, root: bool = False):
        key = (canonical(state), root)
        cached = self.cache.get(key)
        if cached is not None:
            return cached
        self.states += 1

        c0 = counts(state)[0]
        if c0 == 0:
            out = (Fraction(0), Fraction(0), Fraction(0))
            self.cache[key] = out
            return out

        sens, minr = leaf_stats(state, self.k, 0)
        best = None if root else (Fraction(0), Fraction(0), Fraction(0))

        for v in range(self.nleaves):
            if leaf_value(state, self.k, v) != U:
                continue
            p_sens = Fraction(sens[v], c0)
            p_min = Fraction(minr[v], c0)
            total = self.weight * p_sens - alpha * p_min
            pq = p_sens
            pm = p_min

            for a in (ZERO, ONE):
                child = set_leaf(state, self.k, v, a)
                c0a = counts(child)[0]
                if c0a == 0:
                    continue
                w = Fraction(c0a, c0)
                sub_val, sub_pq, sub_pm = self._value(child, alpha, False)
                total += w * sub_val
                pq += w * sub_pq
                pm += w * sub_pm

            if best is None or total > best[0]:
                best = (total, pq, pm)

        if best is None:
            best = (Fraction(0), Fraction(0), Fraction(0))
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
            history.append({
                "iter": it, "alpha_in": str(alpha), "rho": str(value),
                "Pq": str(pq), "Pm": str(pm), "alpha_out": str(new_alpha),
                "states": self.states, "seconds": round(time.time() - t0, 1),
            })
            if verbose:
                print(f"  iter {it}: alpha={alpha} -> rho={value}")
                print(f"           Pq={pq} Pm={pm} -> alpha={new_alpha}")
                print(f"           {self.states} states, {time.time()-t0:.1f}s")
            if new_alpha == alpha:
                return alpha, history
            alpha = new_alpha
        raise RuntimeError("alpha iteration did not converge")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--k", type=int, required=True)
    ap.add_argument("--verbose", action="store_true", default=True)
    args = ap.parse_args()

    t0 = time.time()
    solver = StructSolver(args.k)
    print(f"k={args.k}: {solver.nleaves} leaves, structured state representation")

    alpha, history = solver.find_alpha(verbose=args.verbose)
    elapsed = time.time() - t0

    published = PUBLISHED_ALPHA.get(args.k)
    match = published is not None and alpha == published
    print(f"alpha_{args.k} = {alpha}")
    if published is not None:
        print(f"published  = {published}  -> {'MATCH' if match else 'MISMATCH'}")
    base = 2 + float(alpha) ** (-1.0 / args.k)
    print(f"lower-bound base = {base:.8f}   [float, display only]")
    print(f"elapsed {elapsed:.1f}s, {solver.states} canonical states")

    os.makedirs(RESULTS_DIR, exist_ok=True)
    path = os.path.join(RESULTS_DIR, f"maj3_struct_k{args.k}.json")
    with open(path, "w") as fh:
        json.dump({
            "k": args.k, "alpha_str": str(alpha),
            "alpha_numerator": alpha.numerator, "alpha_denominator": alpha.denominator,
            "published": str(published) if published else None,
            "matches_published": match,
            "base_float_display_only": base,
            "canonical_states": solver.states,
            "elapsed_seconds": round(elapsed, 2),
            "iterations": history,
        }, fh, indent=2)
    print(f"log: {path}")
    return 0 if (published is None or match) else 1


if __name__ == "__main__":
    raise SystemExit(main())
