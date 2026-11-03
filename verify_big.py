"""verify_big.py -- scalable verifier for quaternary edit-metric codes E_4(n, d).

Valid iff: every codeword has length n over {0,1,2,3}, there are no duplicates,
and every pair of distinct codewords has Levenshtein distance >= d.

Method (different from the pairwise verifiers): edit distance is a path metric on
strings of ALL lengths, so with tr = (d-1)//2 and mr = d-1-tr,
    lev(a, b) <= d-1   iff   B_tr(a) meets B_mr(b),
where B_r(w) is every string reachable from w by at most r single-symbol
insertions, deletions or substitutions.  Balls are generated explicitly, each
string is encoded as an integer key (value << 5 | length), and a collision
between the balls of two DIFFERENT codewords is found by sorting.  Integer keys
and exact set membership only: no floating point, no sampling (CLAUDE.md rule 10).

The pairwise verifiers (verify_editq.py, verify_editq2.py) are quadratic in the
code size and too slow for codes of ~10^5 words; verify_big2.c is the
independent second verifier for this scale (rule 3).

    .venv/bin/python verify_big.py FILE [--d D]
exit 0 iff valid, 1 if invalid, 2 if malformed input or parameters refused.
"""
from __future__ import annotations

import hashlib
import re
import sys

import numpy as np

ROWS = 2_000_000  # target rows per batch when expanding balls (bounds memory)


def chunk_for(r: int) -> int:
    """Codewords per batch for radius-r balls (each has up to ~100**r strings)."""
    return max(1, ROWS // (100 ** r))


def parse(path: str, d_override: int | None):
    n = d = q = None
    words = []
    size_hdr = sha_hdr = None
    for raw in open(path):
        line = raw.strip()
        if not line:
            continue
        if line.startswith("#"):
            # fields are read only from key = value lines, at the start or after a
            # comma: '# n = 6' and '# n = 6, alphabet = 0123, size = 120' parse;
            # prose such as '# quaternary code, n=6, size 112' does NOT (refuse, don't guess)
            def field(name):
                if not re.match(r"#\s*[A-Za-z_]+\s*=", line):
                    return None          # prose line: never parsed
                m = re.search(r"(?:^#\s*|,\s*)" + name + r"\s*=\s*(\d+)\s*(?:,|$)", line)
                return int(m.group(1)) if m else None
            n = field("n") if field("n") is not None else n
            d = field("d") if field("d") is not None else d
            q = field("q") if field("q") is not None else q
            size_hdr = field("size") if field("size") is not None else size_hdr
            m = re.search(r"sha256\(codewords\)\s*=\s*([0-9a-f]{64})", line)
            if m: sha_hdr = m.group(1)
            continue
        words.append(line)
    if d_override is not None:
        if d is not None and d != d_override:
            raise ValueError(f"--d {d_override} contradicts '# d = {d}' in the file")
        d = d_override
    if n is None:
        raise ValueError("missing '# n = ...'")
    if d is None:
        raise ValueError("missing '# d = ...' and no --d given; refusing to guess")
    if q not in (None, 4):
        raise ValueError(f"only q = 4 is supported, file says q = {q}")
    return words, n, d, size_hdr, sha_hdr


def keys_of(digits: np.ndarray) -> np.ndarray:
    """Integer key (value << 5 | length) for each row of a (M, L) digit array."""
    M, L = digits.shape
    v = np.zeros(M, dtype=np.int64)
    for j in range(L):
        v = (v << 2) | digits[:, j].astype(np.int64)
    return (v << 5) | L


def ball1(batches):
    """One edit step from every row: {L: (digits, owners)} -> same, unioned."""
    out = {}

    def put(L, dg, ow):
        if dg.shape[0] == 0 or L < 1:
            return
        out.setdefault(L, []).append((dg, ow))

    for L, (dg, ow) in batches.items():
        put(L, dg, ow)                                       # zero edits
        for i in range(L):                                   # deletions
            put(L - 1, np.delete(dg, i, axis=1), ow)
        for i in range(L):                                   # substitutions
            for c in range(4):
                keep = dg[:, i] != c
                t = dg[keep].copy()
                t[:, i] = c
                put(L, t, ow[keep])
        for i in range(L + 1):                               # insertions
            for c in range(4):
                put(L + 1, np.insert(dg, i, c, axis=1), ow)
    return {L: (np.concatenate([a for a, _ in parts]), np.concatenate([b for _, b in parts]))
            for L, parts in out.items()}


def dedupe(batches):
    """Drop repeated (owner, string) rows; returns (keys, owners) flat arrays.

    (owner, key) is packed into one int64 as owner << kb | key, where kb is the
    key width (2*L + 5 bits), so a 1-D np.unique does the work.  If the packing
    would not fit in 63 bits, fall back to the slower 2-column unique."""
    ks, os = [], []
    for L, (dg, ow) in batches.items():
        k = keys_of(dg)
        ow = ow.astype(np.int64)
        kb = 2 * L + 5
        if ow.size and int(ow.max()).bit_length() + kb <= 63:
            packed = np.unique((ow << kb) | k)
            os.append(packed >> kb); ks.append(packed & ((1 << kb) - 1))
        else:
            pair = np.unique(np.stack([ow, k], axis=1), axis=0)
            os.append(pair[:, 0]); ks.append(pair[:, 1])
    return np.concatenate(ks), np.concatenate(os)


def ball_keys(digits: np.ndarray, owners: np.ndarray, r: int):
    batches = {digits.shape[1]: (digits, owners)}
    for _ in range(r):
        batches = ball1(batches)
        # dedupe between steps so radius-2 and -3 balls stay manageable
        k, o = dedupe(batches)
        batches = {}
        for L in np.unique(k & 31):
            sel = (k & 31) == L
            vals = k[sel] >> 5
            dg = np.zeros((sel.sum(), int(L)), dtype=np.uint8)
            for j in range(int(L) - 1, -1, -1):
                dg[:, j] = vals & 3
                vals = vals >> 2
            batches[int(L)] = (dg, o[sel])
    return dedupe(batches)


def check(words, n: int, d: int):
    seen = set()
    for w in words:
        if len(w) != n or set(w) - set("0123"):
            return False, f"malformed codeword {w!r}"
        if w in seen:
            return False, f"duplicate codeword {w!r}"
        seen.add(w)
    if d <= 1 or len(words) < 2:
        return True, ""
    digits = np.array([[int(c) for c in w] for w in words], dtype=np.uint8)
    tr = (d - 1) // 2
    mr = d - 1 - tr

    # A = B_tr of every codeword.  A key shared by two owners => lev <= 2*tr <= d-1.
    ka, oa = [], []
    ca = chunk_for(tr)
    for s in range(0, len(words), ca):
        k, o = ball_keys(digits[s:s + ca], np.arange(s, min(s + ca, len(words))), tr)
        ka.append(k); oa.append(o)
    ka = np.concatenate(ka); oa = np.concatenate(oa)
    order = np.argsort(ka, kind="stable")
    ka, oa = ka[order], oa[order]
    clash = np.nonzero(ka[1:] == ka[:-1])[0]
    if clash.size:
        i = clash[0]
        a, b = int(oa[i]), int(oa[i + 1])
        return False, f"codewords {words[a]!r} and {words[b]!r} are within distance {2 * tr}"
    if mr == tr:
        return True, ""

    # B = B_mr of every codeword; any element lying in another codeword's B_tr
    # means lev <= tr + mr = d-1.
    cb = chunk_for(mr)
    for s in range(0, len(words), cb):
        e = min(s + cb, len(words))
        kb, ob = ball_keys(digits[s:e], np.arange(s, e), mr)
        pos = np.searchsorted(ka, kb)
        pos[pos >= ka.size] = ka.size - 1
        hit = (ka[pos] == kb) & (oa[pos] != ob)
        if hit.any():
            i = np.nonzero(hit)[0][0]
            a, b = int(oa[pos[i]]), int(ob[i])
            return False, f"codewords {words[a]!r} and {words[b]!r} are within distance {d - 1}"
    return True, ""


def main(argv) -> int:
    args = argv[1:]
    d_over = None
    if "--d" in args:
        i = args.index("--d")
        d_over = int(args[i + 1])
        del args[i:i + 2]
    if len(args) != 1:
        print(__doc__)
        return 2
    try:
        words, n, d, size_hdr, sha_hdr = parse(args[0], d_over)
    except (ValueError, OSError) as e:
        print(f"refused: {e}")
        return 2
    sha = hashlib.sha256("".join(w + "\n" for w in sorted(words)).encode()).hexdigest()
    ok, why = check(words, n, d)
    if ok and size_hdr is not None and size_hdr != len(words):
        ok, why = False, f"header says size {size_hdr}, file has {len(words)} codewords"
    if ok and sha_hdr is not None and sha_hdr != sha:
        ok, why = False, "header sha256 does not match the codewords"
    print(f"file      : {args[0]}")
    print(f"n         : {n}")
    print(f"d         : {d}")
    print(f"size      : {len(words)}")
    print(f"valid     : {ok}")
    if not ok:
        print(f"reason    : {why}")
    print(f"checksum  : {sha}")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv))
