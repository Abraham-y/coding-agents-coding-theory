"""Verifier for fixed-length q-ary edit-metric codes E_q(n, d).

A file is a valid artifact iff every codeword has length n over the alphabet
{0,...,q-1}, no duplicates, and every pair has Levenshtein (insert / delete /
substitute) distance >= d.  This is the definition used by Houghten's table and
by Orth's thesis ("minimum edit distance d"), and it is NOT the ball-disjointness
criterion of verify.py: that one is the d = 3 special case only.

File format: '#' comment lines, of which two must read  '# n = <int>'  and
'# d = <int>'  (an optional '# q = <int>' defaults to 4); then one codeword per
line.  The reported checksum is sha256 of the sorted codewords joined by '\n'.

Boring by design: plain Wagner-Fischer, no banding, no early exit.

    python verify_editq.py best/<artifact>
exit 0 iff valid.
"""
from __future__ import annotations

import hashlib
import re
import sys


def levenshtein(a: str, b: str) -> int:
    prev = list(range(len(b) + 1))
    for i in range(1, len(a) + 1):
        cur = [i] + [0] * len(b)
        for j in range(1, len(b) + 1):
            cur[j] = min(prev[j] + 1, cur[j - 1] + 1, prev[j - 1] + (a[i - 1] != b[j - 1]))
        prev = cur
    return prev[len(b)]


def parse(path: str):
    n = d = None
    q = 4
    words = []
    with open(path) as fh:
        for line in fh:
            s = line.strip()
            if not s:
                continue
            if s.startswith("#"):
                m = re.match(r"#\s*n\s*=\s*(\d+)", s)
                if m:
                    n = int(m.group(1))
                m = re.match(r"#\s*d\s*=\s*(\d+)", s)
                if m:
                    d = int(m.group(1))
                m = re.match(r"#\s*q\s*=\s*(\d+)", s)
                if m:
                    q = int(m.group(1))
                continue
            words.append(s)
    return n, d, q, words


def check(words, n: int, d: int, q: int = 4):
    """(valid, size, reason)."""
    alphabet = set(str(i) for i in range(q))
    seen = set()
    for w in words:
        if len(w) != n:
            return False, 0, f"codeword {w!r} has length {len(w)} != {n}"
        if set(w) - alphabet:
            return False, 0, f"codeword {w!r} uses a symbol outside 0..{q-1}"
        if w in seen:
            return False, 0, f"duplicate codeword {w!r}"
        seen.add(w)
    ws = list(words)
    for i in range(len(ws)):
        for j in range(i + 1, len(ws)):
            dist = levenshtein(ws[i], ws[j])
            if dist < d:
                return False, 0, f"{ws[i]} and {ws[j]} are at edit distance {dist} < {d}"
    return True, len(ws), ""


def main(argv) -> int:
    if len(argv) != 2:
        print("usage: python verify_editq.py <artifact>")
        return 2
    n, d, q, words = parse(argv[1])
    if n is None or d is None:
        print("file must declare '# n = <int>' and '# d = <int>'")
        return 2
    ok, size, reason = check(words, n, d, q)
    body = "\n".join(sorted(words)) + "\n"
    print(f"file     : {argv[1]}")
    print(f"n, d, q  : {n}, {d}, {q}")
    print(f"claimed  : {len(words)} codewords")
    print(f"valid    : {ok}{('  (' + reason + ')') if reason else ''}")
    print(f"score    : {size}")
    print(f"checksum : {hashlib.sha256(body.encode()).hexdigest()}")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
