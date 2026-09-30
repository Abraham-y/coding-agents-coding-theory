"""Second, independent verifier for E_q(n, d) artifacts (CLAUDE.md rule 3).

Shares no code with verify_editq.py.  Different representation of the same
criterion: the pairwise edit distance is computed by Myers' bit-parallel
algorithm (Hyyro's formulation), where the DP column is held as bit-vectors of
vertical deltas rather than as an integer array.  If both verifiers agree on the
accept set and the reject set, the encoding of "minimum edit distance d" is at
least consistent across two formulations.

Same file format as verify_editq.py.  exit 0 iff valid.
"""
from __future__ import annotations

import hashlib
import sys


def edit_distance_bitparallel(pattern: str, text: str) -> int:
    """Levenshtein distance via Myers (1999) / Hyyro (2003).  Works for any
    pattern length since Python ints are unbounded."""
    m = len(pattern)
    if m == 0:
        return len(text)
    peq = {}
    for i, ch in enumerate(pattern):
        peq[ch] = peq.get(ch, 0) | (1 << i)
    mask = (1 << m) - 1
    pv = mask
    mv = 0
    score = m
    top = 1 << (m - 1)
    for ch in text:
        eq = peq.get(ch, 0)
        xv = eq | mv
        xh = (((eq & pv) + pv) & mask) ^ pv | eq
        ph = mv | ~(xh | pv) & mask
        mh = pv & xh
        if ph & top:
            score += 1
        if mh & top:
            score -= 1
        ph = ((ph << 1) | 1) & mask
        mh = (mh << 1) & mask
        pv = mh | ~(xv | ph) & mask
        mv = ph & xv
    return score


def read_artifact(path):
    n = d = None
    q = 4
    words = []
    for raw in open(path):
        line = raw.strip()
        if not line:
            continue
        if line[0] == "#":
            body = line[1:].replace(" ", "")
            if body.startswith("n=") and body[2:].isdigit():
                n = int(body[2:])
            elif body.startswith("d=") and body[2:].isdigit():
                d = int(body[2:])
            elif body.startswith("q=") and body[2:].isdigit():
                q = int(body[2:])
            continue
        words.append(line)
    return n, d, q, words


def validate(words, n, d, q):
    symbols = {chr(ord("0") + i) for i in range(q)}
    if len(set(words)) != len(words):
        return False, "duplicate codeword"
    for w in words:
        if len(w) != n or any(c not in symbols for c in w):
            return False, f"bad codeword {w}"
    for i, a in enumerate(words):
        for b in words[i + 1:]:
            if edit_distance_bitparallel(a, b) < d:
                return False, f"pair {a},{b} too close"
    return True, ""


def main() -> int:
    if len(sys.argv) != 2:
        print("usage: python verify_editq2.py <artifact>")
        return 2
    n, d, q, words = read_artifact(sys.argv[1])
    if n is None or d is None:
        print("missing '# n =' or '# d =' header")
        return 2
    ok, why = validate(words, n, d, q)
    digest = hashlib.sha256(("\n".join(sorted(words)) + "\n").encode()).hexdigest()
    print(f"file    : {sys.argv[1]}")
    print(f"n,d,q   : {n},{d},{q}")
    print(f"valid   : {ok}{(' (' + why + ')') if why else ''}")
    print(f"score   : {len(words) if ok else 0}")
    print(f"sha256  : {digest}")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
