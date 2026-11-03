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

Edit-metric artifacts (Houghten's tables): valid iff every codeword is a word of
length n over 0 .. q-1 and every pair of distinct codewords has Levenshtein
distance >= d, where d and q are the artifact's '# d = ...' and '# q = ...'
fields.  Legacy single-edit artifacts that declare no d are checked at d = 3,
which is what they claim, and artifacts that declare no q are quaternary, which
all of them were until 2026-10-30.  q = 2 must be declared explicitly, so that a
binary deletion artifact with its '# s' line removed is still refused.  (Until
2026-10-10 this path ignored '# d' and always checked d = 3, so a d >= 4 claim
was never actually verified here.)

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
    if s is None and set("".join(words)) <= set("01") and declared_q(text) != 2:
        raise ValueError("binary artifact must declare '# s = ...' "
                         "(or '# q = 2' for a binary edit-metric code)")
    return words, n, s


ALPHABET_Q = "0123"


def edit_ball(word: str, alphabet: str = ALPHABET_Q):
    """All strings over `alphabet` within one insertion, deletion or substitution
    of `word`.

    Self-contained on purpose. The edit-metric artifacts in best/ previously had
    no ground-truth verifier outside the search module, which made CLAUDE.md
    rule 5 ("python verify.py best/<artifact> must work") false for that family.
    """
    out = {word}
    for i in range(len(word)):
        out.add(word[:i] + word[i + 1:])
        for c in alphabet:
            if c != word[i]:
                out.add(word[:i] + c + word[i + 1:])
    for i in range(len(word) + 1):
        for c in alphabet:
            out.add(word[:i] + c + word[i:])
    return out


def edit_ball_r(word: str, r: int, alphabet: str = ALPHABET_Q):
    """All strings over `alphabet` within r insertions, deletions or
    substitutions of `word`, built by applying edit_ball r times."""
    ball, frontier = {word}, {word}
    for _ in range(r):
        reached = set()
        for w in frontier:
            reached |= edit_ball(w, alphabet)
        frontier = reached - ball
        ball |= frontier
    return ball


def verify_edit(words, n: int, d: int = 3, q: int = 4):
    """(valid, size, reason) for a q-ary edit-metric code (q = 2, 3 or 4) of
    length n, minimum Levenshtein distance d.

    Edit distance is a path metric on strings of all lengths.  With
    tr = (d-1)//2 and mr = d-1-tr, two codewords a, b are within distance d-1
    exactly when B_tr(a) and B_mr(b) share a string.

    Phase 1 marks every string of each codeword's B_tr ball in a bitset, and
    fails on a string that is already marked: it then lies in another
    codeword's B_tr ball, so the two are within 2*tr <= d-1.  For d = 3 this is
    the original criterion, radius-1 balls pairwise disjoint.  Phase 2 (only
    when mr > tr, i.e. even d) looks up each codeword's B_mr ball minus its own
    B_tr ball; a marked string there lies in some OTHER codeword's B_tr ball,
    so the two are within tr + mr = d-1.

    The balls are taken over the q symbols only.  That loses nothing: an
    optimal edit path from a to b only substitutes or inserts symbols of b, so
    every string on it is q-ary, and the string tr steps along it lies in
    B_tr(a) and B_mr(b) whenever lev(a, b) <= d-1.

    Memory is one bit per possible string of each length n-tr .. n+tr (a
    bytearray of q**L bits), not one dict entry per ball string, so a
    500,000-word code at n = 13 needs about 40 MB rather than about 6 GB.  On a
    failure a second, slow scan names the other codeword."""
    if q not in (2, 3, 4):
        return (False, 0, f"q = {q} is not supported (2, 3 or 4)")
    alphabet = ALPHABET_Q[:q]
    ws = list(words)
    seen = set()
    for w in ws:
        if len(w) != n or set(w) - set(alphabet):
            return (False, 0, f"malformed codeword {w!r}")
        if w in seen:
            return (False, 0, f"duplicate codeword {w!r}")
        seen.add(w)
    if d < 1:
        return (False, 0, f"d = {d} is not a valid minimum distance")
    tr = (d - 1) // 2
    mr = d - 1 - tr
    marks = {L: bytearray(q ** L // 8 + 1) for L in range(max(0, n - tr), n + tr + 1)}

    def index(x):
        return int(x, q) if x else 0

    def other_owner(x, me):
        """Slow path, failures only: a codeword other than `me` whose B_tr ball holds x."""
        for j, w in enumerate(ws):
            if j != me and x in edit_ball_r(w, tr, alphabet):
                return j
        raise AssertionError(f"string {x!r} is marked but no other codeword reaches it")

    for i, w in enumerate(ws):
        for x in edit_ball_r(w, tr, alphabet):
            k = index(x)
            row = marks[len(x)]
            if row[k >> 3] >> (k & 7) & 1:
                j = other_owner(x, i)
                return (False, 0, f"codewords {ws[j]!r} and {w!r} share {x!r}")
            row[k >> 3] |= 1 << (k & 7)
    if mr > tr:
        for i, w in enumerate(ws):
            for x in edit_ball_r(w, mr, alphabet) - edit_ball_r(w, tr, alphabet):
                row = marks.get(len(x))
                if row is None:
                    continue          # outside lengths n-tr .. n+tr: in no B_tr ball
                k = index(x)
                if row[k >> 3] >> (k & 7) & 1:
                    j = other_owner(x, i)
                    return (False, 0, f"codewords {ws[j]!r} and {w!r} are within "
                                      f"distance {d - 1} (both reach {x!r})")
    return (True, len(ws), "")


def _declared(text: str, key: str):
    """The value an artifact declares for `key` ('d' or 'q'), or None.  Same rule
    as parse_artifact: '# d = 4' must fill a whole comma/semicolon-separated
    field, and conflicting declarations are an error."""
    val = None
    for raw in text.splitlines():
        line = raw.strip()
        if not line.startswith("#"):
            continue
        for field in re.split(r"[,;]", line[1:].strip()):
            m = re.fullmatch(r"\s*\(?\s*" + key + r"\s*=\s*(\d+)\s*\)?\s*", field)
            if not m:
                continue
            value = int(m.group(1))
            if val is not None and val != value:
                raise ValueError(f"artifact declares {key} twice with different values ({val} and {value})")
            val = value
    return val


def declared_d(text: str):
    """The minimum distance an edit-metric artifact declares, or None."""
    return _declared(text, "d")


def declared_q(text: str):
    """The alphabet size an edit-metric artifact declares, or None (quaternary)."""
    return _declared(text, "q")


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
        d = declared_d(text)
        q = declared_q(text)
        if q is not None and q not in (2, 3, 4):
            raise ValueError(f"'# q = {q}': edit-metric artifacts must have q = 2, 3 or 4")
        valid, score, reason = verify_edit(words, n, 3 if d is None else d, 4 if q is None else q)
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
    # declared d is enforced: lev('000000','000222') = 3
    check("edit d=3 accept 3-apart", verify_edit(["000000", "000222"], 6, 3)[0], True)
    check("edit d=4 reject 3-apart", verify_edit(["000000", "000222"], 6, 4)[0], False)
    check("edit d=6 reject 3-apart", verify_edit(["000000", "000222"], 6, 6)[0], False)
    # lev('000000','002222') = 4 (even d uses the asymmetric B_1 / B_2 test)
    check("edit d=4 accept 4-apart", verify_edit(["000000", "002222"], 6, 4)[0], True)
    check("edit d=5 reject 4-apart", verify_edit(["000000", "002222"], 6, 5)[0], False)
    # an indel-only neighbour: '012301' vs '123010' (one deletion + one insertion, lev 2)
    check("edit d=3 reject indel pair", verify_edit(["012301", "123010"], 6, 3)[0], False)
    check("edit d=2 accept indel pair", verify_edit(["012301", "123010"], 6, 2)[0], True)
    # the failure names the offending pair, not a bystander
    r = verify_edit(["000000", "111111", "000222"], 6, 4)
    check("edit d=4 reject names pair", (r[0], "'000000'" in r[2] and "'000222'" in r[2]), (False, True))
    r = verify_edit(["111111", "000000", "000022"], 6, 3)
    check("edit d=3 reject names pair", (r[0], "'000000'" in r[2] and "'000022'" in r[2]), (False, True))

    # q = 2 and q = 3: the alphabet is enforced, distances are the same
    check("edit q=3 accept 000/111/222", verify_edit(["000", "111", "222"], 3, 3, 3)[0], True)
    check("edit q=3 reject at d=4", verify_edit(["000", "111", "222"], 3, 4, 3)[0], False)
    check("edit q=3 reject symbol 3", verify_edit(["000", "333"], 3, 3, 3)[0], False)
    check("edit q=2 accept 0000/1111 d=4", verify_edit(["0000", "1111"], 4, 4, 2)[0], True)
    check("edit q=2 reject 0000/1111 d=5", verify_edit(["0000", "1111"], 4, 5, 2)[0], False)
    check("edit q=2 reject symbol 2", verify_edit(["0000", "2222"], 4, 3, 2)[0], False)
    # lev('0101','1010') = 2 via a deletion and an insertion; binary balls must see it
    check("edit q=2 reject indel pair", verify_edit(["0101", "1010"], 4, 3, 2)[0], False)
    check("edit q=2 accept indel pair d=2", verify_edit(["0101", "1010"], 4, 2, 2)[0], True)
    # lev('0011','1100') = 4; d=4 accept exercises the binary B_1 / B_2 phase
    check("edit q=2 d=4 accept 0011/1100", verify_edit(["0011", "1100"], 4, 4, 2)[0], True)
    check("edit q=3 d=6 reject 3-apart", verify_edit(["000000", "000222"], 6, 6, 3)[0], False)
    check("edit q=5 refused", verify_edit(["0000"], 4, 3, 5)[0], False)
    # a binary edit artifact must say q = 2; without it, it is refused as a deletion artifact
    try:
        parse_artifact("# edit metric\n# n = 4\n# d = 4\n0000\n1111\n")
        check("binary edit without q refused", False, True)
    except ValueError:
        pass
    check("binary edit with q=2 parses",
          parse_artifact("# edit metric\n# n = 4\n# d = 4\n# q = 2\n0000\n1111\n")[1:], (4, None))

    if failures:
        for f in failures:
            print("FAIL", f)
        return 1
    print("selftest ok")
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
    if s is None:
        d = declared_d(open(argv[1]).read())
        print(f"d         : {3 if d is None else d}" + ("  (legacy single-edit artifact, no '# d')" if d is None else ""))
        q = declared_q(open(argv[1]).read())
        print(f"q         : {4 if q is None else q}" + ("  (no '# q', quaternary)" if q is None else ""))
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
