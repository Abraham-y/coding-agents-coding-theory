"""verify.py -- ground truth checker for binary s-deletion-correcting codes.

Definition (the one this file implements, no other):
    A set C of binary strings, all of length n, is s-deletion-correcting iff
    no two DISTINCT codewords share a common subsequence of length n - s.

Equivalently: the deletion balls D_s(v) = {all length-(n-s) subsequences of v}
are pairwise disjoint across distinct codewords.

This file is ground truth (CLAUDE.md rule 1). It is deliberately boring.
No floating point appears anywhere in an accept/reject decision (rule 10).

CLI:
    python verify.py <artifact.txt>          # verify one artifact
    python verify.py --selftest              # tiny built-in sanity pairs

Artifact format (see best/README for the canonical form):
    # n = 12
    # s = 1
    # anything else after '#' is a free-text comment
    010110101101
    ...one codeword per line, whitespace ignored...
"""

from __future__ import annotations

import hashlib
import re
import sys
from itertools import combinations

# --------------------------------------------------------------------------
# core
# --------------------------------------------------------------------------


def deletion_ball(word: str, s: int) -> set[str]:
    """Every distinct subsequence of `word` obtained by deleting exactly s symbols.

    Deleting s symbols from a length-n word always yields length n-s. We
    enumerate the C(n, s) index sets and dedupe, because different index sets
    frequently produce the same subsequence (e.g. '000' -> '00' three ways).
    """
    n = len(word)
    if s < 0:
        raise ValueError("s must be >= 0")
    if s > n:
        raise ValueError("cannot delete more symbols than the word has")
    if s == 0:
        return {word}
    keep_count = n - s
    out: set[str] = set()
    for keep in combinations(range(n), keep_count):
        out.add("".join(word[i] for i in keep))
    return out


def verify(words, n: int, s: int):
    """Return (valid, score, reason).

    valid  -- True iff `words` is a well-formed s-deletion-correcting code of length n
    score  -- number of codewords (0 when invalid)
    reason -- '' when valid, else a specific human-readable cause
    """
    if not isinstance(n, int) or not isinstance(s, int):
        return (False, 0, "n and s must be ints")
    if n <= 0:
        return (False, 0, f"n must be positive, got {n}")
    if s < 0:
        return (False, 0, f"s must be non-negative, got {s}")
    if s >= n:
        return (False, 0, f"s must be < n, got s={s} n={n}")

    words = list(words)

    # --- well-formedness -------------------------------------------------
    seen_words = set()
    for w in words:
        if not isinstance(w, str):
            return (False, 0, f"codeword is not a str: {w!r}")
        if len(w) != n:
            return (False, 0, f"codeword {w!r} has length {len(w)}, expected {n}")
        if set(w) - {"0", "1"}:
            return (False, 0, f"codeword {w!r} is not binary")
        if w in seen_words:
            return (False, 0, f"duplicate codeword {w!r}")
        seen_words.add(w)

    # --- the actual code condition ---------------------------------------
    # Map each length-(n-s) subsequence to the index of the codeword that
    # produced it. A subsequence claimed by two distinct codewords is exactly
    # a violation of the definition above.
    owner: dict[str, int] = {}
    for idx, w in enumerate(words):
        for sub in deletion_ball(w, s):
            prev = owner.get(sub)
            if prev is not None and prev != idx:
                return (
                    False,
                    0,
                    f"codewords {words[prev]!r} and {w!r} share subsequence {sub!r}",
                )
            owner[sub] = idx

    return (True, len(words), "")


# --------------------------------------------------------------------------
# artifact I/O
# --------------------------------------------------------------------------


def parse_artifact(text: str):
    """Parse an artifact file. Returns (words, n, s). Raises ValueError.

    A parameter must occupy an ENTIRE comma/semicolon-separated field of a comment.
    An earlier version searched for `n =` anywhere in the line, so prose such as
    "# for context, s = 1 was the baseline" silently set s and re-classified the
    artifact. First declaration wins; a conflicting redeclaration is an error.
    """
    n = None
    s = None
    words = []
    for raw in text.splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith("#"):
            body = line[1:].strip()
            for field in re.split(r"[,;]", body):
                m = re.fullmatch(r"\s*\(?\s*([ns])\s*=\s*(\d+)\s*\)?\s*", field)
                if not m:
                    continue
                key, value = m.group(1), int(m.group(2))
                cur = n if key == "n" else s
                if cur is not None and cur != value:
                    raise ValueError(
                        f"artifact declares {key} twice with different values "
                        f"({cur} and {value})")
                if key == "n":
                    n = value
                else:
                    s = value
            continue
        words.append(line)
    if n is None:
        raise ValueError("artifact must declare '# n = ...'")
    if s is None and set("".join(words)) <= set("01"):
        raise ValueError("binary artifact must declare '# s = ...'")
    return words, n, s


ALPHABET_Q = "0123"


def edit_ball(word: str):
    """All strings within one insertion, deletion or substitution of `word`.

    Self-contained on purpose. The edit-metric artifacts in best/ previously had
    no ground-truth verifier outside the search module, which made CLAUDE.md
    rule 5 ("python verify.py best/<artifact> must work") false for that family.
    """
    out = {word}
    for i in range(len(word)):
        out.add(word[:i] + word[i + 1:])
        for c in ALPHABET_Q:
            if c != word[i]:
                out.add(word[:i] + c + word[i + 1:])
    for i in range(len(word) + 1):
        for c in ALPHABET_Q:
            out.add(word[:i] + c + word[i:])
    return out


def verify_edit(words, n: int):
    """(valid, size, reason) for a single-edit-correcting q=4 code of length n."""
    ws = list(words)
    seen = set()
    for w in ws:
        if len(w) != n or set(w) - set(ALPHABET_Q):
            return (False, 0, f"malformed codeword {w!r}")
        if w in seen:
            return (False, 0, f"duplicate codeword {w!r}")
        seen.add(w)
    balls = [edit_ball(w) for w in ws]
    for i in range(len(ws)):
        for j in range(i + 1, len(ws)):
            shared = balls[i] & balls[j]
            if shared:
                return (False, 0,
                        f"codewords {ws[i]!r} and {ws[j]!r} share {sorted(shared)[0]!r}")
    return (True, len(ws), "")


def canonical_form(words, n: int, s: int, note: str = "") -> str:
    """Canonical serialisation: sorted codewords, declared params (CLAUDE.md rule 5)."""
    lines = [f"# n = {n}", f"# s = {s}"]
    if note:
        for chunk in note.splitlines():
            lines.append(f"# {chunk}")
    lines.append(f"# size = {len(words)}")
    lines.extend(sorted(words))
    return "\n".join(lines) + "\n"


def checksum(words, n=None, s=None) -> str:
    """SHA-256 of the sorted, newline-terminated codewords.

    Same convention the artifact headers use. It previously differed, so
    `python verify.py <artifact>` printed a digest contradicting the file's own
    header while nothing compared them -- under which a substituted (relabelled,
    still valid) code passed unnoticed.
    """
    return hashlib.sha256("".join(w + "\n" for w in sorted(words)).encode()).hexdigest()


def header_checksum(text: str):
    """The sha256 an artifact claims for itself, if it declares one."""
    m = re.search(r"sha256\(codewords\)\s*=\s*([0-9a-f]{64})", text)
    return m.group(1) if m else None


def verify_file(path: str):
    with open(path) as fh:
        text = fh.read()
    header = "\n".join(l for l in text.splitlines() if l.strip().startswith("#"))
    words, n, s = parse_artifact(text)
    if s is None:
        # Never INFER the metric from the alphabet: dropping '# s =' from a binary
        # artifact and corrupting one symbol used to re-classify it as an edit code,
        # which then passed. The metric has to be stated.
        if not re.search(r"edit", header, re.I):
            raise ValueError("no '# s = ...' and no edit metric declared; "
                             "refusing to guess the problem family")
        claimed = header_checksum(header)
        if claimed and claimed != checksum(words):
            raise ValueError(f"checksum mismatch: header says {claimed[:16]}..., "
                             f"codewords hash to {checksum(words)[:16]}...")
        valid, score, reason = verify_edit(words, n)
        return valid, score, reason, n, None, words
    valid, score, reason = verify(words, n, s)
    return valid, score, reason, n, s, words


# --------------------------------------------------------------------------
# self test: a handful of pairs checkable by hand
# --------------------------------------------------------------------------


def _selftest() -> int:
    failures = []

    def check(label, got, want):
        if got != want:
            failures.append(f"{label}: got {got}, want {want}")

    # deletion_ball, verifiable by eye
    check("ball('000000',1)", deletion_ball("000000", 1), {"00000"})
    check("ball('010',1)", deletion_ball("010", 1), {"10", "00", "01"})
    check("ball('01',0)", deletion_ball("01", 0), {"01"})
    check("ball('0011',2)", deletion_ball("0011", 2), {"00", "01", "11"})

    # '000000' and '000001' both contain '00000' -> not 1-deletion-correcting
    check("reject shared-subseq pair", verify(["000000", "000001"], 6, 1)[0], False)
    # distance is not the criterion: these two share no length-5 subsequence
    check("accept 000000/111111", verify(["000000", "111111"], 6, 1)[0], True)
    # empty and singleton codes are trivially valid
    check("empty", verify([], 6, 1), (True, 0, ""))
    check("singleton", verify(["010101"], 6, 1), (True, 1, ""))
    # malformed input is rejected, not silently coerced
    check("wrong length", verify(["01010"], 6, 1)[0], False)
    check("non binary", verify(["01010x"], 6, 1)[0], False)
    check("duplicate", verify(["010101", "010101"], 6, 1)[0], False)
    check("s>=n", verify(["01"], 2, 2)[0], False)
    # 2 deletions: '000000' and '000011' both contain '0000'
    check("reject s=2 pair", verify(["000000", "000011"], 6, 2)[0], False)
    # but they ARE fine for s=1: no shared length-5 subsequence
    check("accept same pair at s=1", verify(["000000", "000011"], 6, 1)[0], True)

    # edit-metric path, both directions
    check("edit accept singleton", verify_edit(["000000"], 6), (True, 1, ""))
    check("edit reject 1-apart", verify_edit(["000000", "000001"], 6)[0], False)
    check("edit reject 2-apart", verify_edit(["000000", "000011"], 6)[0], False)
    check("edit accept 3-apart", verify_edit(["000000", "000111"], 6)[0], True)
    check("edit reject duplicate", verify_edit(["012301", "012301"], 6)[0], False)
    check("edit reject wrong length", verify_edit(["01230"], 6)[0], False)
    check("edit reject bad symbol", verify_edit(["01234x"], 6)[0], False)

    if failures:
        for f in failures:
            print("FAIL", f)
        return 1
    print(f"selftest ok ({12} checks)")
    return 0


def main(argv) -> int:
    if len(argv) == 2 and argv[1] == "--selftest":
        return _selftest()
    if len(argv) != 2:
        print(__doc__)
        return 2
    valid, score, reason, n, s, words = verify_file(argv[1])
    print(f"file      : {argv[1]}")
    print(f"n         : {n}")
    print(f"s         : {s}")
    print(f"claimed   : {len(words)} codewords")
    print(f"valid     : {valid}")
    print(f"score     : {score}")
    print(f"checksum  : {checksum(words, n, s)}")
    if not valid:
        print(f"reason    : {reason}")
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
