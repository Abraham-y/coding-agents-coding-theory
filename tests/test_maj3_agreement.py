"""Three-solver agreement for alpha_k, plus a direct test of the forced-action rule.

Same discipline as verify.py vs verify2.py: three implementations that share no
state representation must produce identical exact rationals.

  maj3_alpha.AlphaSolver   state = set of consistent hard inputs   (k <= 2)
  maj3_struct.StructSolver state = canonical recursive tree object (k <= 3)
  maj3_fast.FastSolver     the above + MNSX forced-action rules    (k <= 4)

The third one is the interesting case. Rule B (force a query on any leaf that
provably cannot be the absolute minority) rests on an exchange argument that
MNSX assert and this repo does not prove. Two things test it here:

  1. INDIRECTLY: FastSolver must still return alpha_1 = 2, alpha_2 = 24/7 and
     alpha_3 = 12231/2203. A wrong rule moves those.
  2. DIRECTLY: every leaf that find_forced() returns must have minority-count
     exactly 0 in the current configuration. That is the rule's own premise, and
     it is checkable without any exchange argument. Tested over states reached by
     an actual DP walk, not over hand-picked examples.

Run:  python tests/test_maj3_agreement.py
"""

from __future__ import annotations

import os
import random
import sys
from fractions import Fraction

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from search.maj3_alpha import AlphaSolver  # noqa: E402
from search.maj3_fast import FastSolver, find_forced  # noqa: E402
from search.maj3_struct import (  # noqa: E402
    ONE,
    PUBLISHED_ALPHA,
    U,
    ZERO,
    counts,
    fresh_state,
    leaf_stats,
    leaf_value,
    set_leaf,
)

SEED = 20260803
FAILURES: list[str] = []
CHECKS = 0


def check(label, got, want):
    global CHECKS
    CHECKS += 1
    if got != want:
        FAILURES.append(f"{label}\n     got : {got}\n     want: {want}")


# ---------------------------------------------------------------------------


def test_three_solvers_agree():
    for k in (1, 2):
        a1, _ = AlphaSolver(k).find_alpha()
        a2, _ = __import__(
            "search.maj3_struct", fromlist=["StructSolver"]
        ).StructSolver(k).find_alpha()
        a3, _ = FastSolver(k).find_alpha()
        check(f"k={k}: brute force == published", a1, PUBLISHED_ALPHA[k])
        check(f"k={k}: structured == brute force", a2, a1)
        check(f"k={k}: forced-rules == structured", a3, a2)

    # k=3 is out of reach for the brute-force solver; compare the other two.
    a2, _ = __import__(
        "search.maj3_struct", fromlist=["StructSolver"]
    ).StructSolver(3).find_alpha()
    a3, _ = FastSolver(3).find_alpha()
    check("k=3: structured == published", a2, PUBLISHED_ALPHA[3])
    check("k=3: forced-rules == structured", a3, a2)


def test_forced_leaves_really_cannot_be_the_minority():
    """The direct test of Rule B's premise, over states from a real DP walk."""
    rng = random.Random(SEED)
    for k in (2, 3):
        state = fresh_state(k)
        checked = 0
        # walk down a few hundred random query paths, testing every state on the way
        for _ in range(150):
            state = fresh_state(k)
            for _ in range(3 ** k):
                c0, c1 = counts(state)
                if c0 == 0:
                    break
                sens, minr = leaf_stats(state, k, 0)
                forced = find_forced(state, k)
                if forced is not None:
                    checked += 1
                    if leaf_value(state, k, forced) != U:
                        FAILURES.append(
                            f"k={k}: find_forced returned an ALREADY QUERIED leaf {forced}")
                        return
                    if minr[forced] != 0:
                        FAILURES.append(
                            f"k={k}: forced leaf {forced} has minority-count "
                            f"{minr[forced]} != 0 -- RULE B PREMISE VIOLATED")
                        return
                # advance: query a random unqueried leaf, random consistent value
                unq = [v for v in range(3 ** k) if leaf_value(state, k, v) == U]
                if not unq:
                    break
                v = rng.choice(unq)
                opts = []
                for a in (ZERO, ONE):
                    if counts(set_leaf(state, k, v, a))[0] > 0:
                        opts.append(a)
                if not opts:
                    break
                state = set_leaf(state, k, v, rng.choice(opts))
        check(f"k={k}: every forced leaf had minority-count 0 ({checked} tested)",
              True, True)
        if checked == 0:
            FAILURES.append(f"k={k}: no forced leaves were ever produced -- test vacuous")


def test_rule_a_premise():
    """If c_1 == 0 then no unqueried leaf is sensitive. Rule A's premise, checked."""
    rng = random.Random(SEED + 1)
    for k in (2, 3):
        hits = 0
        for _ in range(200):
            state = fresh_state(k)
            for _ in range(3 ** k):
                c0, c1 = counts(state)
                if c0 == 0:
                    break
                if c1 == 0:
                    hits += 1
                    sens, _ = leaf_stats(state, k, 0)
                    bad = [v for v in range(3 ** k)
                           if leaf_value(state, k, v) == U and sens[v] != 0]
                    if bad:
                        FAILURES.append(
                            f"k={k}: c_1==0 but unqueried leaves {bad[:3]} are sensitive "
                            f"-- RULE A PREMISE VIOLATED")
                        return
                    break
                unq = [v for v in range(3 ** k) if leaf_value(state, k, v) == U]
                if not unq:
                    break
                v = rng.choice(unq)
                opts = [a for a in (ZERO, ONE)
                        if counts(set_leaf(state, k, v, a))[0] > 0]
                if not opts:
                    break
                state = set_leaf(state, k, v, rng.choice(opts))
        check(f"k={k}: Rule A premise held wherever it fired ({hits} times)", True, True)


def test_counts_match_direct_enumeration():
    """counts() is combinatorial; cross-check it against literal enumeration at k<=2."""
    from search.maj3_alpha import build_hard_inputs

    for k in (1, 2):
        inputs, values, _, _ = build_hard_inputs(k)
        state = fresh_state(k)
        c0, c1 = counts(state)
        check(f"k={k}: counts c_0 matches enumeration",
              c0, sum(1 for v in values if v == 0))
        check(f"k={k}: counts c_1 matches enumeration",
              c1, sum(1 for v in values if v == 1))

        # and after pinning a leaf
        for leaf in range(3 ** k):
            for a in (ZERO, ONE):
                st = set_leaf(state, k, leaf, a)
                got0 = counts(st)[0]
                want0 = sum(1 for i, x in enumerate(inputs)
                            if values[i] == 0 and str(x[leaf]) == a)
                if got0 != want0:
                    FAILURES.append(
                        f"k={k} leaf {leaf}={a}: counts {got0} != enumeration {want0}")
                    return
        check(f"k={k}: counts match enumeration after pinning any single leaf", True, True)


# ---------------------------------------------------------------------------


def main() -> int:
    test_counts_match_direct_enumeration()
    test_rule_a_premise()
    test_forced_leaves_really_cannot_be_the_minority()
    test_three_solvers_agree()

    print(f"checks run : {CHECKS}")
    print(f"failures   : {len(FAILURES)}")
    for f in FAILURES:
        print("  FAIL " + f)
    return 1 if FAILURES else 0


if __name__ == "__main__":
    raise SystemExit(main())
