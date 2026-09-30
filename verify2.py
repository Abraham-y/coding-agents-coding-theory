"""verify2.py -- INDEPENDENT second verifier (CLAUDE.md rule 3).

Deliberately shares no code and no data structure with verify.py.

    verify.py  : expands each codeword's deletion ball and checks that the
                 balls are pairwise disjoint, via a hash map of subsequences.
    verify2.py : never materialises a deletion ball. It computes the length of
                 the longest common subsequence of each PAIR of codewords by
                 dynamic programming, and rejects if any pair has LCS >= n - s.

The two are equivalent because two words share a common subsequence of length
n - s exactly when their LCS is at least n - s. Getting the same verdict from
both is evidence the PROBLEM ENCODING is right, which is the failure mode
rule 3 exists to catch.

Cost is O(|C|^2 * n^2), quadratic in the code size, so this is the slow
reference. It is used on every test object and on every artifact we report;
it is not used inside search loops.

CAVEAT, recorded honestly: rule 3 asks for this to be written in a fresh
session with no shared context. This file was written in the same session as
verify.py, so it is weaker evidence than the rule intends. PROMPTS.md prompt
4b regenerates it cleanly; do that before reporting any record externally.

CLI:
    python verify2.py <artifact.txt>
    python verify2.py --selftest
"""

from __future__ import annotations

import re

import sys

# --------------------------------------------------------------------------


def lcs_length(u: str, v: str) -> int:
    """Length of the longest common subsequence of u and v. Plain DP, no tricks."""
    m, n = len(u), len(v)
    prev = [0] * (n + 1)
    for i in range(1, m + 1):
        cur = [0] * (n + 1)
        ui = u[i - 1]
        for j in range(1, n + 1):
            if ui == v[j - 1]:
                cur[j] = prev[j - 1] + 1
            else:
                a = prev[j]
                b = cur[j - 1]
                cur[j] = a if a >= b else b
        prev = cur
    return prev[n]


def verify2(words, n: int, s: int):
    """Return (valid, score, reason). Same contract as verify.verify."""
    if not isinstance(n, int) or not isinstance(s, int):
        return (False, 0, "n and s must be ints")
    if n <= 0:
        return (False, 0, f"n must be positive, got {n}")
    if s < 0:
        return (False, 0, f"s must be non-negative, got {s}")
    if s >= n:
        return (False, 0, f"s must be < n, got s={s} n={n}")

    words = list(words)

    for w in words:
        if not isinstance(w, str):
            return (False, 0, f"codeword is not a str: {w!r}")
        if len(w) != n:
            return (False, 0, f"codeword {w!r} has length {len(w)}, expected {n}")
        for ch in w:
            if ch != "0" and ch != "1":
                return (False, 0, f"codeword {w!r} is not binary")

    for i in range(len(words)):
        for j in range(i + 1, len(words)):
            if words[i] == words[j]:
                return (False, 0, f"duplicate codeword {words[i]!r}")

    threshold = n - s
    for i in range(len(words)):
        for j in range(i + 1, len(words)):
            if lcs_length(words[i], words[j]) >= threshold:
                return (
                    False,
                    0,
                    f"codewords {words[i]!r} and {words[j]!r} have LCS >= {threshold}",
                )

    return (True, len(words), "")


# --------------------------------------------------------------------------


def levenshtein(u: str, v: str) -> int:
    """Edit distance by DP. Deliberately a different method from verify.py's
    ball-intersection test, so the two files remain independent checks of the
    same criterion (a code corrects one edit iff min pairwise distance >= 3)."""
    prev = list(range(len(v) + 1))
    for i, cu in enumerate(u, 1):
        cur = [i]
        for j, cv in enumerate(v, 1):
            cur.append(min(prev[j] + 1, cur[j - 1] + 1, prev[j - 1] + (cu != cv)))
        prev = cur
    return prev[-1]


def verify2_edit(words, n: int, q: int = 4):
    """(valid, size, reason) for a single-edit-correcting q-ary code of length n."""
    alpha = set("0123"[:q])
    ws = list(words)
    seen = set()
    for w in ws:
        if len(w) != n or set(w) - alpha:
            return (False, 0, f"malformed codeword {w!r}")
        if w in seen:
            return (False, 0, f"duplicate codeword {w!r}")
        seen.add(w)
    for i in range(len(ws)):
        for j in range(i + 1, len(ws)):
            d = levenshtein(ws[i], ws[j])
            if d < 3:
                return (False, 0, f"{ws[i]!r} and {ws[j]!r} at edit distance {d} < 3")
    return (True, len(ws), "")


def _selftest() -> int:
    failures = []

    def check(label, got, want):
        if got != want:
            failures.append(f"{label}: got {got}, want {want}")

    check("lcs(abc-ish) 010/010", lcs_length("010", "010"), 3)
    check("lcs 000000/000001", lcs_length("000000", "000001"), 5)
    check("lcs 000000/111111", lcs_length("000000", "111111"), 0)
    check("lcs 0101/1010", lcs_length("0101", "1010"), 3)
    check("lcs empty", lcs_length("", "0101"), 0)

    check("reject shared-subseq pair", verify2(["000000", "000001"], 6, 1)[0], False)
    check("accept 000000/111111", verify2(["000000", "111111"], 6, 1)[0], True)
    check("empty", verify2([], 6, 1), (True, 0, ""))
    check("singleton", verify2(["010101"], 6, 1), (True, 1, ""))
    check("wrong length", verify2(["01010"], 6, 1)[0], False)
    check("non binary", verify2(["01010x"], 6, 1)[0], False)
    check("duplicate", verify2(["010101", "010101"], 6, 1)[0], False)
    check("reject s=2 pair", verify2(["000000", "000011"], 6, 2)[0], False)
    check("accept same pair at s=1", verify2(["000000", "000011"], 6, 1)[0], True)

    check("edit accept 3-apart", verify2_edit(["000000", "000111"], 6)[0], True)
    check("edit reject 2-apart", verify2_edit(["000000", "000011"], 6)[0], False)
    check("edit reject 1-apart", verify2_edit(["000000", "000001"], 6)[0], False)
    check("edit reject duplicate", verify2_edit(["012301", "012301"], 6)[0], False)
    check("edit reject bad symbol", verify2_edit(["01234x"], 6)[0], False)
    check("lev sanity 011100/001110", levenshtein("011100", "001110"), 2)
    if failures:
        for f in failures:
            print("FAIL", f)
        return 1
    print("verify2 selftest ok")
    return 0


def main(argv) -> int:
    if len(argv) == 2 and argv[1] == "--selftest":
        return _selftest()
    if len(argv) != 2:
        print(__doc__)
        return 2

    # parse the artifact without importing verify.py's parser
    n = None
    s = None
    words = []
    with open(argv[1]) as fh:
        for raw in fh:
            line = raw.strip()
            if not line:
                continue
            if line.startswith("#"):
                body = line[1:].strip()
                # whole-field match only; prose must not be able to set a parameter
                for field in re.split(r"[,;]", body):
                    m = re.fullmatch(r"\s*\(?\s*([ns])\s*=\s*(\d+)\s*\)?\s*", field)
                    if not m:
                        continue
                    k, v = m.group(1), int(m.group(2))
                    cur = n if k == "n" else s
                    if cur is not None and cur != v:
                        print(f"artifact declares {k} twice with different values")
                        return 2
                    if k == "n":
                        n = v
                    else:
                        s = v
                continue
            words.append(line)
    if n is None:
        print("artifact must declare '# n = ...'")
        return 2

    if s is None:
        # The metric is never inferred from the alphabet: a binary artifact that
        # loses its '# s =' and gains a corrupted symbol would otherwise be
        # re-classified as an edit code and accepted.
        hdr = "\n".join(l for l in open(argv[1]) if l.strip().startswith("#"))
        if not re.search(r"edit", hdr, re.I):
            print("no '# s = ...' and no edit metric declared; refusing to guess")
            return 2
        valid, score, reason = verify2_edit(words, n)
    else:
        valid, score, reason = verify2(words, n, s)
    print(f"file   : {argv[1]}")
    print(f"n      : {n}")
    print(f"s      : {s}")
    print(f"valid  : {valid}")
    print(f"score  : {score}")
    if not valid:
        print(f"reason : {reason}")
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
