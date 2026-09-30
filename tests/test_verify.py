"""Two-directional test suite for verify.py (CLAUDE.md rule 2).

ACCEPT SET  -- known-good objects from the literature. verify.py must accept
               them AND report exactly the published size (rule 4).
REJECT SET  -- deliberately corrupted objects. verify.py must refuse them.

Run:  python tests/test_verify.py [--full]
      --full also rebuilds VT_0(19) and VT_0(20), which is slow.
"""

from __future__ import annotations

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from verify import deletion_ball, verify  # noqa: E402
from vt import PUBLISHED_VT0_SIZES, SETTLED_OPTIMA, vt_code  # noqa: E402

FAILURES: list[str] = []
FINDINGS: list[str] = []
CHECKS = 0


def check(label: str, got, want):
    global CHECKS
    CHECKS += 1
    if got != want:
        FAILURES.append(f"{label}\n     got : {got!r}\n     want: {want!r}")


def owner_map(words, n, s):
    owner = {}
    for i, w in enumerate(words):
        for sub in deletion_ball(w, s):
            owner[sub] = i
    return owner


# ---------------------------------------------------------------------------
# ACCEPT SET
# ---------------------------------------------------------------------------


def test_accept_vt_codes(max_n: int):
    for n in range(6, max_n + 1):
        code = vt_code(n, 0)
        valid, score, reason = verify(code, n, 1)
        check(f"accept: VT_0({n}) is a valid 1-deletion code ({reason})", valid, True)
        # rule 4: reproduce the published number exactly, not approximately
        check(f"ingest: |VT_0({n})| matches published table", score, PUBLISHED_VT0_SIZES[n])


def test_accept_other_vt_classes():
    # Levenshtein: every residue class a gives a single-deletion-correcting code.
    for n in (7, 9, 11, 12):
        for a in range(n + 1):
            code = vt_code(n, a)
            valid, _, reason = verify(code, n, 1)
            check(f"accept: VT_{a}({n}) valid ({reason})", valid, True)


def test_accept_settled_optima_are_achievable():
    # Where the literature says A_1(n) is settled, VT_0(n) must attain it.
    for n, optimum in SETTLED_OPTIMA[1].items():
        code = vt_code(n, 0)
        _, score, _ = verify(code, n, 1)
        check(f"ingest: VT_0({n}) attains settled optimum A_1({n})={optimum}", score, optimum)


def test_accept_subsets_stay_valid():
    # Any subset of a valid code is a valid code. Cheap structural sanity.
    code = vt_code(12, 0)
    for cut in (1, 7, 100, 315):
        valid, score, _ = verify(code[:cut], 12, 1)
        check(f"accept: VT_0(12) truncated to {cut} stays valid", (valid, score), (True, cut))


# ---------------------------------------------------------------------------
# REJECT SET
# ---------------------------------------------------------------------------


def test_reject_handmade_pairs():
    # Each of these is checkable by eye.
    check("reject: 000000/000001 share 00000", verify(["000000", "000001"], 6, 1)[0], False)
    check("reject: 010101/010101 duplicate", verify(["010101", "010101"], 6, 1)[0], False)
    check("reject: wrong length", verify(["0101010"], 6, 1)[0], False)
    check("reject: non-binary symbol", verify(["01010z"], 6, 1)[0], False)
    check("reject: s >= n", verify(["0101"], 4, 4)[0], False)
    check("reject: negative s", verify(["0101"], 4, -1)[0], False)
    check("reject: 000000/000011 share 0000 at s=2", verify(["000000", "000011"], 6, 2)[0], False)


def test_reject_vt_at_higher_s():
    # A single-deletion code is not automatically a two-deletion code.
    for n in (8, 10, 12):
        code = vt_code(n, 0)
        valid, _, _ = verify(code, n, 2)
        check(f"reject: VT_0({n}) is NOT a 2-deletion code", valid, False)


def test_reject_union_of_two_vt_classes():
    # For n where A_1(n) is settled, VT_0 u VT_1 exceeds the proven optimum,
    # so it CANNOT be a valid code. This reject is forced by the literature.
    for n, optimum in SETTLED_OPTIMA[1].items():
        union = vt_code(n, 0) + vt_code(n, 1)
        if len(union) <= optimum:
            continue  # argument does not apply; skip rather than assert wrongly
        valid, _, _ = verify(union, n, 1)
        check(
            f"reject: VT_0({n}) u VT_1({n}) has {len(union)} > A_1({n})={optimum}, must be invalid",
            valid,
            False,
        )


def test_reject_corrupted_vt():
    """Swap one codeword for an outside word that provably collides."""
    for n in (8, 10, 12):
        code = vt_code(n, 0)
        inside = set(code)
        owner = owner_map(code, n, 1)
        # find the first outside word that collides with the code
        intruder = None
        for x in range(1 << n):
            w = format(x, f"0{n}b")
            if w in inside:
                continue
            if any(sub in owner for sub in deletion_ball(w, 1)):
                intruder = w
                break
        if intruder is None:
            FINDINGS.append(f"VT_0({n}) has NO colliding outside word -- probe further")
            continue
        corrupted = code + [intruder]
        valid, _, _ = verify(corrupted, n, 1)
        check(f"reject: VT_0({n}) + colliding word {intruder}", valid, False)

        # and a true swap, same size as the original
        swapped = code[:-1] + [intruder]
        valid2, _, _ = verify(swapped, n, 1)
        check(f"reject: VT_0({n}) with last word swapped for {intruder}", valid2, False)


# ---------------------------------------------------------------------------
# MAXIMALITY PROBE  (a reject test that doubles as the cheapest possible search)
# ---------------------------------------------------------------------------


def test_vt_maximality_probe(max_n: int):
    """Is VT_0(n) a MAXIMAL independent set, i.e. can no word be added?

    If any word can be added, |VT_0(n)|+1 beats the best known size for the
    open cells n >= 12 and that is a reportable result. Expect: none.
    """
    for n in range(6, max_n + 1):
        code = vt_code(n, 0)
        inside = set(code)
        owner = owner_map(code, n, 1)
        addable = []
        for x in range(1 << n):
            w = format(x, f"0{n}b")
            if w in inside:
                continue
            if not any(sub in owner for sub in deletion_ball(w, 1)):
                addable.append(w)
        if addable:
            FINDINGS.append(
                f"VT_0({n}) is NOT maximal: {len(addable)} addable words, "
                f"e.g. {addable[:3]} -> code of size {len(code) + 1}"
            )
        check(f"probe: VT_0({n}) is a maximal independent set", addable == [], True)


# ---------------------------------------------------------------------------


def main() -> int:
    full = "--full" in sys.argv
    max_n = 20 if full else 16

    test_accept_vt_codes(18 if full else 16)
    test_accept_other_vt_classes()
    test_accept_settled_optima_are_achievable()
    test_accept_subsets_stay_valid()

    test_reject_handmade_pairs()
    test_reject_vt_at_higher_s()
    test_reject_union_of_two_vt_classes()
    test_reject_corrupted_vt()

    test_vt_maximality_probe(max_n if full else 14)

    print(f"checks run : {CHECKS}")
    print(f"failures   : {len(FAILURES)}")
    for f in FAILURES:
        print("  FAIL " + f)
    if FINDINGS:
        print("FINDINGS (investigate, do not celebrate):")
        for f in FINDINGS:
            print("  * " + f)
    return 1 if FAILURES else 0


if __name__ == "__main__":
    raise SystemExit(main())
