"""Tests for the recursive-majority alpha_k computation.

Accept set = the published values of Magniez, Nayak, Santha, Sherman, Tardos and
Xiao (arXiv:1309.7565). Reproducing them exactly, in exact rational arithmetic,
is the rule-4 ingestion check for this track.

Cross-checks that are independent of each other:
  1. alpha_1 = 2 and alpha_2 = 24/7                      (their Theorem 3.8)
  2. rho_alpha for the k=2 optimum equals (48 - 14a)/81  (their Figure 2, by hand)
  3. N_k recurrence gives 2, 7, 112, 246792, 2505258478767772  (their section 3.2)
  4. hard-input counts |H_k^0| follow 3, 81, 3^13         (derivable independently)

Run:  python tests/test_maj3.py
"""

from __future__ import annotations

import os
import sys
from fractions import Fraction

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from search.maj3_alpha import (  # noqa: E402
    PUBLISHED_ALPHA,
    PUBLISHED_N,
    AlphaSolver,
    build_hard_inputs,
    maj3,
    stable_config_counts,
)

FAILURES: list[str] = []
CHECKS = 0


def check(label, got, want):
    global CHECKS
    CHECKS += 1
    if got != want:
        FAILURES.append(f"{label}\n     got : {got}\n     want: {want}")


# ---------------------------------------------------------------------------


def test_maj3_truth_table():
    check("maj3(0,0,0)", maj3(0, 0, 0), 0)
    check("maj3(1,0,0)", maj3(1, 0, 0), 0)
    check("maj3(1,1,0)", maj3(1, 1, 0), 1)
    check("maj3(1,1,1)", maj3(1, 1, 1), 1)


def test_hard_input_structure():
    """|H_k| = 6 * (|H_{k-1}|/2)^3, and every hard input has exactly 2^k sensitive bits."""
    expected_hard = {1: 6, 2: 162, 3: 3188646}
    for k in (1, 2, 3):
        inputs, values, minority, sensitive = build_hard_inputs(k)
        check(f"|H_{k}|", len(inputs), expected_hard[k])
        zeros = sum(1 for v in values if v == 0)
        check(f"|H_{k}^0| is half of |H_{k}|", zeros, expected_hard[k] // 2)
        # every hard input has exactly 2^k sensitive leaves
        bad = [i for i, s in enumerate(sensitive) if len(s) != 2 ** k]
        check(f"every hard input of height {k} has exactly 2^{k} sensitive bits",
              bad[:3], [])
        if k <= 2:
            # flipping a sensitive bit really does flip the value; flipping a
            # non-sensitive bit really does not. Checked directly, no shortcuts.
            for i, x in enumerate(inputs):
                for pos in range(3 ** k):
                    flipped = list(x)
                    flipped[pos] ^= 1
                    newval = _eval(tuple(flipped), k)
                    should_flip = pos in sensitive[i]
                    if (newval != values[i]) != should_flip:
                        FAILURES.append(
                            f"sensitivity wrong at k={k} input {i} pos {pos}")
                        return


def _eval(x, k):
    if k == 0:
        return x[0]
    w = 3 ** (k - 1)
    return maj3(_eval(x[:w], k - 1), _eval(x[w:2 * w], k - 1), _eval(x[2 * w:], k - 1))


def test_minority_is_consistent():
    """The absolute minority bit must itself be a bit whose flip is NOT sensitive."""
    for k in (1, 2):
        inputs, values, minority, sensitive = build_hard_inputs(k)
        for i in range(len(inputs)):
            if minority[i] in sensitive[i]:
                FAILURES.append(f"k={k} input {i}: minority bit marked sensitive")
                return
        check(f"k={k}: minority never sensitive", True, True)


def test_stable_config_counts():
    got = stable_config_counts(5)
    for k, want in PUBLISHED_N.items():
        check(f"N_{k}", got[k], want)


def test_alpha_1():
    solver = AlphaSolver(1)
    alpha, _ = solver.find_alpha()
    check("alpha_1 == 2", alpha, PUBLISHED_ALPHA[1])


def test_alpha_2():
    solver = AlphaSolver(2)
    alpha, hist = solver.find_alpha()
    check("alpha_2 == 24/7", alpha, PUBLISHED_ALPHA[2])

    # independent cross-check against the paper's hand-derived Figure 2:
    # rho_alpha(C') = (48 - 14*alpha)/81, i.e. Pq/4 = 16/27 and Pm = 14/81.
    final = hist[-1]
    pq = Fraction(final["Pq"])
    pm = Fraction(final["Pm"])
    check("k=2 optimum: 2^-2 * Pq == 16/27 (Figure 2 constant term)",
          Fraction(1, 4) * pq, Fraction(16, 27))
    check("k=2 optimum: Pm == 14/81 (Figure 2 alpha coefficient)",
          pm, Fraction(14, 81))
    # and therefore rho_alpha == (48 - 14 alpha)/81 identically in alpha
    for a in (Fraction(0), Fraction(3), Fraction(24, 7), Fraction(7, 2)):
        lhs = Fraction(1, 4) * pq - a * pm
        rhs = (Fraction(48) - 14 * a) / 81
        check(f"k=2: rho_{a} matches Figure 2 formula", lhs, rhs)


def test_alpha_iteration_converges_fast():
    """MNSX report at most four iterations for k = 2,3,4. Ours should match at k<=2."""
    for k in (1, 2):
        solver = AlphaSolver(k)
        _, hist = solver.find_alpha()
        if len(hist) > 4:
            FAILURES.append(f"k={k}: took {len(hist)} iterations, paper says <= 4")
        check(f"k={k} converges in <= 4 iterations", len(hist) <= 4, True)


# ---------------------------------------------------------------------------


def main() -> int:
    test_maj3_truth_table()
    test_hard_input_structure()
    test_minority_is_consistent()
    test_stable_config_counts()
    test_alpha_1()
    test_alpha_2()
    test_alpha_iteration_converges_fast()

    print(f"checks run : {CHECKS}")
    print(f"failures   : {len(FAILURES)}")
    for f in FAILURES:
        print("  FAIL " + f)
    return 1 if FAILURES else 0


if __name__ == "__main__":
    raise SystemExit(main())
