"""Two-directional suite for the general-distance verifiers verify_editq.py and
verify_editq2.py (CLAUDE.md rules 2 and 3).

Accept set : the d = 3 artifacts in best/ (which verify.py already accepts),
             and hand-made codes at d = 4, 5, 6.
Reject set : each accepted object with one codeword replaced by a word at
             distance d-1 from another; the same object with d raised by one;
             a duplicate codeword; a codeword of the wrong length; a symbol
             outside the alphabet.
Agreement  : both verifiers, plus the d = 3 ball-disjointness verifier, must
             agree on random codes over several (n, d).

Run:  python tests/test_editq_verify.py
"""
from __future__ import annotations

import glob
import os
import random
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import verify_editq as v1  # noqa: E402
import verify_editq2 as v2  # noqa: E402
from search.edit4 import verify_edit4  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FAILURES: list[str] = []
CHECKS = 0


def check(label, got, want):
    global CHECKS
    CHECKS += 1
    if got != want:
        FAILURES.append(f"{label}: got {got}, want {want}")


def both(words, n, d):
    a = v1.check(list(words), n, d)[0]
    b = v2.validate(list(words), n, d, 4)[0]
    check(f"agreement n={n} d={d} on {len(words)} words", a, b)
    return a


def test_distance_functions():
    rng = random.Random(7)
    cases = [("", "", 0), ("abc", "abc", 0), ("kitten", "sitting", 3),
             ("0123", "1230", 2), ("000000", "000001", 1), ("0000", "1111", 4),
             ("012", "0312", 1), ("0123", "023", 1)]
    for a, b, want in cases:
        check(f"lev1({a},{b})", v1.levenshtein(a, b), want)
        check(f"lev2({a},{b})", v2.edit_distance_bitparallel(a, b), want)
    for _ in range(3000):
        n1, n2 = rng.randint(0, 8), rng.randint(0, 8)
        a = "".join(rng.choice("0123") for _ in range(n1))
        b = "".join(rng.choice("0123") for _ in range(n2))
        check(f"lev agreement {a},{b}", v1.levenshtein(a, b), v2.edit_distance_bitparallel(a, b))


def test_accept_artifacts():
    files = sorted(glob.glob(os.path.join(ROOT, "best", "edit4_n6_*size1*.txt")))
    check("have d=3 artifacts", len(files) > 0, True)
    for f in files:
        words = [l.strip() for l in open(f) if l.strip() and not l.startswith("#")]
        ok3 = verify_edit4(words, 6)[0]
        okq = both(words, 6, 3)
        check(f"accept {os.path.basename(f)} at d=3", okq, True)
        check(f"ball-verifier agrees on {os.path.basename(f)}", ok3, okq)
        # the same code is NOT a d=4 code (rule 2: shift the parameter by one)
        check(f"reject {os.path.basename(f)} at d=4", both(words, 6, 4), False)


def test_handmade():
    # repetition code: min distance n (all n symbols differ, and any n-1 edits
    # cannot map 0^n to 1^n)
    for n in (4, 5, 6):
        code = [c * n for c in "0123"]
        check(f"repetition n={n} at d={n}", both(code, n, n), True)
        check(f"repetition n={n} at d={n+1}", both(code, n, n + 1), False)
    # a d=4 code at n=6 by greedy construction, then re-check both ways
    rng = random.Random(1)
    words = ["".join(rng.choice("0123") for _ in range(6)) for _ in range(4000)]
    code = []
    for w in words:
        if all(v1.levenshtein(w, c) >= 4 for c in code):
            code.append(w)
    check("greedy d=4 code nonempty", len(code) >= 10, True)
    check("greedy d=4 code accepted", both(code, 6, 4), True)
    check("greedy d=4 code rejected at d=5", both(code, 6, 5), False)


def test_reject_corrupted():
    rng = random.Random(3)
    f = os.path.join(ROOT, "best", "edit4_n6_V4_size120.txt")
    words = [l.strip() for l in open(f) if l.strip() and not l.startswith("#")]
    # one codeword replaced by a neighbour (distance 1) of another codeword
    bad = list(words)
    base = bad[5]
    bad[10] = base[:2] + ("0" if base[2] != "0" else "1") + base[3:]
    check("substitution-corrupted rejected", both(bad, 6, 3), False)
    # a word at distance exactly 2 from a codeword (delete one symbol, insert
    # another elsewhere) -- must be rejected at d=3 and accepted at d=2 unless
    # it happens to collide with something else
    bad = list(words)
    base = bad[7]
    cand = base[1:] + "2"
    if v1.levenshtein(cand, base) == 2 and cand not in bad:
        bad[20] = cand
        check("distance-2 word rejected at d=3", both(bad, 6, 3), False)
    # duplicate
    dup = list(words) + [words[0]]
    check("duplicate rejected", both(dup, 6, 3), False)
    # wrong length
    check("wrong length rejected", both(words + ["01234"], 6, 3), False)
    # out-of-alphabet symbol
    check("bad symbol rejected", both(words + ["01234x"], 6, 3), False)
    # deleting a codeword keeps validity (subset of a code is a code)
    fewer = list(words)
    fewer.pop(rng.randrange(len(fewer)))
    check("subset accepted", both(fewer, 6, 3), True)


def test_differential_random():
    rng = random.Random(20260907)
    for _ in range(300):
        n = rng.randint(3, 7)
        d = rng.randint(2, n)
        k = rng.randint(2, 12)
        words = list({"".join(rng.choice("0123") for _ in range(n)) for _ in range(k)})
        okq = both(words, n, d)
        if d == 3:
            check(f"ball verifier agrees n={n}", verify_edit4(words, n)[0], okq)


def main() -> int:
    for t in (test_distance_functions, test_accept_artifacts, test_handmade,
              test_reject_corrupted, test_differential_random):
        t()
    print(f"{CHECKS} checks, {len(FAILURES)} failures")
    for f in FAILURES[:20]:
        print("  FAIL", f)
    return 1 if FAILURES else 0


if __name__ == "__main__":
    raise SystemExit(main())
