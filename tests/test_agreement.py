"""Differential test: verify.py and verify2.py must agree on every object.

This is the rule-3 check. The two verifiers use different characterisations of
the same definition (deletion-ball disjointness vs. pairwise LCS threshold), so
agreement across a large randomised object set is evidence that the shared
PROBLEM ENCODING is right -- which is the failure that would silently poison
the whole sprint.

The randomisation here is seeded and fully reproducible. It tests verifier
agreement, never an artifact's verdict, so rule 10 is not in play.

Run:  python tests/test_agreement.py [--full]
"""

from __future__ import annotations

import os
import random
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from verify import verify  # noqa: E402
from verify2 import verify2  # noqa: E402
from vt import vt_code  # noqa: E402

SEED = 20260803
FAILURES: list[str] = []
CHECKS = 0


def agree(label: str, words, n: int, s: int):
    """Both verifiers must return the same (valid, score)."""
    global CHECKS
    CHECKS += 1
    v1, sc1, r1 = verify(words, n, s)
    v2, sc2, r2 = verify2(words, n, s)
    if (v1, sc1) != (v2, sc2):
        FAILURES.append(
            f"{label}\n     verify : valid={v1} score={sc1} reason={r1!r}"
            f"\n     verify2: valid={v2} score={sc2} reason={r2!r}"
        )
    return v1


# ---------------------------------------------------------------------------


def test_literature_objects(max_n: int):
    for n in range(6, max_n + 1):
        code = vt_code(n, 0)
        got = agree(f"VT_0({n}) at s=1", code, n, 1)
        if not got:
            FAILURES.append(f"VT_0({n}) rejected at s=1 by both -- encoding is wrong")
        agree(f"VT_0({n}) at s=2", code, n, 2)
        agree(f"VT_0({n}) at s=3", code, n, 3)


def test_handmade():
    cases = [
        (["000000", "000001"], 6, 1),
        (["000000", "111111"], 6, 1),
        (["000000", "000011"], 6, 2),
        ([], 6, 1),
        (["010101"], 6, 1),
        (["01010"], 6, 1),
        (["01010x"], 6, 1),
        (["010101", "010101"], 6, 1),
        (["0101"], 4, 4),
        (["0101"], 4, -1),
        (["0000", "0001", "0011", "0111", "1111"], 4, 1),
        (["00", "01", "10", "11"], 2, 1),
    ]
    for words, n, s in cases:
        agree(f"handmade {words} n={n} s={s}", words, n, s)


def test_random_differential(rounds: int):
    """The heart of the check: random codes, most invalid, some valid."""
    rng = random.Random(SEED)
    for n in range(4, 10):
        for s in (1, 2, 3):
            if s >= n:
                continue
            universe = [format(x, f"0{n}b") for x in range(1 << n)]
            for r in range(rounds):
                k = rng.randint(1, min(len(universe), 12))
                words = rng.sample(universe, k)
                agree(f"random n={n} s={s} round={r} k={k}", words, n, s)


def test_random_subsets_of_vt(rounds: int):
    """Subsets of a valid code must be judged valid by both."""
    rng = random.Random(SEED + 1)
    for n in (8, 10, 12):
        code = vt_code(n, 0)
        for r in range(rounds):
            k = rng.randint(1, min(len(code), 40))
            words = rng.sample(code, k)
            got = agree(f"VT_0({n}) subset k={k} round={r}", words, n, 1)
            if not got:
                FAILURES.append(f"subset of VT_0({n}) judged invalid -- impossible")


def test_random_near_misses(rounds: int):
    """Valid code plus one outside word: the interesting boundary."""
    rng = random.Random(SEED + 2)
    for n in (8, 10, 12):
        code = vt_code(n, 0)
        inside = set(code)
        outside = [format(x, f"0{n}b") for x in range(1 << n) if format(x, f"0{n}b") not in inside]
        for r in range(rounds):
            k = rng.randint(2, min(len(code), 30))
            words = rng.sample(code, k) + [rng.choice(outside)]
            agree(f"VT_0({n}) subset+intruder round={r}", words, n, 1)


# ---------------------------------------------------------------------------


def main() -> int:
    full = "--full" in sys.argv
    test_handmade()
    test_literature_objects(13 if full else 12)
    test_random_differential(40 if full else 15)
    test_random_subsets_of_vt(30 if full else 12)
    test_random_near_misses(30 if full else 12)

    print(f"agreement checks : {CHECKS}")
    print(f"disagreements    : {len(FAILURES)}")
    for f in FAILURES:
        print("  FAIL " + f)
    return 1 if FAILURES else 0


if __name__ == "__main__":
    raise SystemExit(main())
