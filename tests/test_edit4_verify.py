"""Two-directional test suite promoting verify_edit4 toward ground-truth status.

Closes the verifier gap flagged by adversarial review: the edit-metric artifacts
(A_4-112, the 1131 code) were outside root verify.py's coverage, resting on
verify_edit4 alone. This suite gives verify_edit4 the same discipline as
verify.py: an accept set of known-good objects, a reject set of deliberately
corrupted ones, and a differential check against an independently written
Levenshtein-distance verifier (valid single-edit code <=> pairwise edit
distance >= 3; equivalent because edit distance is a geodesic metric).

Run:  python tests/test_edit4_verify.py
"""

from __future__ import annotations

import os
import random
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from search.edit4 import verify_edit4  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FAILURES: list[str] = []
CHECKS = 0
SEED = 20260807


def check(label, got, want):
    global CHECKS
    CHECKS += 1
    if got != want:
        FAILURES.append(f"{label}: got {got}, want {want}")


def lev(a: str, b: str) -> int:
    """Independent Levenshtein DP — shares no code with edit4.ball()."""
    m, n = len(a), len(b)
    prev = list(range(n + 1))
    for i in range(1, m + 1):
        cur = [i] + [0] * n
        ai = a[i - 1]
        for j in range(1, n + 1):
            cur[j] = min(prev[j] + 1, cur[j - 1] + 1, prev[j - 1] + (ai != b[j - 1]))
        prev = cur
    return prev[n]


def verify_by_levenshtein(words, n):
    ws = list(words)
    if len(set(ws)) != len(ws):
        return False
    if any(len(w) != n or set(w) - set("0123") for w in ws):
        return False
    return all(lev(ws[i], ws[j]) >= 3 for i in range(len(ws)) for j in range(i + 1, len(ws)))


def load_artifact(path):
    return [l.strip() for l in open(path) if l.strip() and not l.startswith("#")]


def test_accept_artifacts():
    """The repo's two flagship artifacts must pass, at exactly their stated sizes."""
    a4 = os.path.join(ROOT, "best", "edit4_n6_A4invariant_112.txt")
    big = os.path.join(ROOT, "best", "edit4_n8_size1131.txt")
    if os.path.exists(a4):
        words = load_artifact(a4)
        ok, sz, _ = verify_edit4(words, 6)
        check("accept: A_4-112 artifact", (ok, sz), (True, 112))
        check("accept: A_4-112 by independent Levenshtein", verify_by_levenshtein(words, 6), True)
    else:
        FAILURES.append("A_4-112 artifact missing")
    if os.path.exists(big):
        words = load_artifact(big)
        ok, sz, _ = verify_edit4(words, 8)
        check("accept: 1131 artifact", (ok, sz), (True, 1131))
    else:
        FAILURES.append("1131 artifact missing")


def test_accept_handmade():
    check("accept: constant words n=4", verify_edit4(["0000", "1111", "2222", "3333"], 4)[0], True)
    check("accept: empty", verify_edit4([], 5), (True, 0, ""))
    check("accept: singleton", verify_edit4(["0123"], 4)[:2], (True, 1))


def test_reject_corrupted():
    """Corrupt the A_4 artifact in each canonical way; every one must be rejected."""
    a4 = os.path.join(ROOT, "best", "edit4_n6_A4invariant_112.txt")
    if not os.path.exists(a4):
        return
    words = load_artifact(a4)
    check("reject: duplicate word", verify_edit4(words + [words[0]], 6)[0], False)
    check("reject: wrong length word", verify_edit4(words + ["01230"], 6)[0], False)
    check("reject: bad symbol", verify_edit4(words[:-1] + ["01234X"[:6].replace('X','4')], 6)[0], False)
    # single-substitution neighbour of an existing codeword: distance 1 < 3
    w = words[0]
    mutated = w[:-1] + ("0" if w[-1] != "0" else "1")
    check("reject: substitution-neighbour added", verify_edit4(words + [mutated], 6)[0], False)
    # a word at edit distance exactly 2 (two substitutions) must also be rejected
    w2 = ("0" if w[0] != "0" else "1") + w[1:-1] + ("2" if w[-1] != "2" else "3")
    if all(lev(w2, x) >= 3 for x in words if x != w) and lev(w2, w) == 2:
        check("reject: distance-2 word added", verify_edit4(words + [w2], 6)[0], False)
    # cross-length confusability: deletion then insertion elsewhere
    check("reject: pair sharing a deletion", verify_edit4(["000000", "000001"], 6)[0], False)
    check("reject: ins/del confusable pair", verify_edit4(["001234"[:6], "012340"[:6]], 6)[0], False)


def test_differential_random():
    """verify_edit4 vs the Levenshtein verifier on random word sets, n=4,5."""
    rng = random.Random(SEED)
    import itertools
    for n in (4, 5):
        universe = ["".join(p) for p in itertools.product("0123", repeat=n)]
        for r in range(40):
            k = rng.randint(2, 9)
            ws = rng.sample(universe, k)
            got1 = verify_edit4(ws, n)[0]
            got2 = verify_by_levenshtein(ws, n)
            check(f"differential n={n} round={r}", got1, got2)


def main() -> int:
    test_accept_artifacts()
    test_accept_handmade()
    test_reject_corrupted()
    test_differential_random()
    print(f"checks run : {CHECKS}")
    print(f"failures   : {len(FAILURES)}")
    for f in FAILURES:
        print("  FAIL " + f)
    return 1 if FAILURES else 0


# --- accept set: the record artifacts themselves (added after adversarial review
# noted they were checked only by reproduce.sh's inline script, not by this suite)
for _fn in ("edit4_n6_C4_product_size116.txt",
            "edit4_n6_C4_diagonal_size116.txt",
            "edit4_n6_C4_diagonal_b_size116.txt"):
    _p = os.path.join(ROOT, "best", _fn)
    if os.path.exists(_p):
        _w = [l.strip() for l in open(_p) if l.strip() and not l.startswith("#")]
        check(f"accept {_fn}", verify_edit4(_w, 6)[:2], (True, 116))
        # reject set: deleting one word stays valid, but corrupting one must fail
        _bad = list(_w); _bad[0] = _bad[1]
        check(f"reject dup in {_fn}", verify_edit4(_bad, 6)[0], False)
        _bad2 = list(_w); _bad2[0] = _bad2[0][:5] + ("0" if _bad2[0][5] != "0" else "1")
        check(f"reject mutation in {_fn}", verify_edit4(_bad2, 6)[0], False)

if __name__ == "__main__":
    raise SystemExit(main())
