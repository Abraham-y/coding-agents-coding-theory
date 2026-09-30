"""Export the quaternary single-edit conflict graph to METIS format for KaMIS,
and verify KaMIS results back through our own verifier.

Bucket-based construction (never builds bitset adjacency), so it scales to
n = 8 (65,536 vertices) and beyond: for every received string t, the codewords
whose ball contains t form a clique; the union of in-bucket pairs is the edge
set. Edges are deduplicated with numpy on integer codes.

Round trip:
    python search/edit4_metis.py --n 7 --export           # writes graphs/edit4_n7.metis
    <run redumis on it>
    python search/edit4_metis.py --n 7 --check result.mis # verify + report size

The .metis vertex on line i (1-indexed) corresponds to the word with base-4
index i-1; the result file is KaMIS's one-0/1-per-line format.
"""

from __future__ import annotations

import argparse
import itertools
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from search.edit4 import ALPHABET, PUBLISHED, ball, verify_edit4  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
GRAPH_DIR = os.path.join(ROOT, "graphs")


def word_of(idx: int, n: int) -> str:
    out = []
    for _ in range(n):
        out.append(ALPHABET[idx & 3])
        idx >>= 2
    return "".join(reversed(out))


def idx_of(word: str) -> int:
    v = 0
    for c in word:
        v = (v << 2) | int(c)
    return v


def build_edges(n: int):
    """All conflict edges as a deduplicated numpy array of (u, v), u < v."""
    t0 = time.time()
    buckets: dict[str, list[int]] = {}
    for idx in range(4 ** n):
        w = word_of(idx, n)
        for t in ball(w):
            buckets.setdefault(t, []).append(idx)
    print(f"  buckets built: {len(buckets)} ({time.time()-t0:.0f}s)")

    chunks = []
    for members in buckets.values():
        if len(members) < 2:
            continue
        m = sorted(members)
        for u, v in itertools.combinations(m, 2):
            chunks.append(u * (4 ** n) + v)
    codes = np.unique(np.array(chunks, dtype=np.int64))
    u = codes // (4 ** n)
    v = codes % (4 ** n)
    print(f"  edges: {len(codes)} unique ({time.time()-t0:.0f}s)")
    return u.astype(np.int64), v.astype(np.int64)


def export(n: int):
    os.makedirs(GRAPH_DIR, exist_ok=True)
    nv = 4 ** n
    u, v = build_edges(n)
    m = len(u)
    # adjacency lists via CSR on the symmetrised edge list
    src = np.concatenate([u, v])
    dst = np.concatenate([v, u])
    order = np.lexsort((dst, src))  # KaMIS requires each neighbor list sorted
    src, dst = src[order], dst[order]
    starts = np.searchsorted(src, np.arange(nv + 1))
    path = os.path.join(GRAPH_DIR, f"edit4_n{n}.metis")
    t0 = time.time()
    with open(path, "w") as fh:
        fh.write(f"{nv} {m}\n")
        for i in range(nv):
            nbrs = dst[starts[i]:starts[i + 1]] + 1  # METIS is 1-indexed
            fh.write(" ".join(map(str, nbrs)) + "\n")
    print(f"  wrote {path} ({time.time()-t0:.0f}s)")
    return path


def check(n: int, result_path: str):
    with open(result_path) as fh:
        flags = [line.strip() for line in fh if line.strip() in ("0", "1")]
    chosen = [i for i, f in enumerate(flags) if f == "1"]
    words = [word_of(i, n) for i in chosen]
    ok, size, reason = verify_edit4(words, n)
    pub = PUBLISHED.get(n, {}).get("best")
    print(f"result: {len(chosen)} vertices selected")
    print(f"verify_edit4: valid={ok} size={size} {reason}")
    print(f"published best known: {pub}")
    if ok and pub and size > pub:
        out = os.path.join(ROOT, "best", f"edit4_n{n}_size{size}.txt")
        os.makedirs(os.path.dirname(out), exist_ok=True)
        with open(out, "w") as fh:
            fh.write(f"# quaternary single-edit code, n={n}, size {size} (KaMIS)\n")
            for w in sorted(words):
                fh.write(w + "\n")
        print(f"CANDIDATE RECORD saved: {out} -- novelty check required before any claim")
    return 0 if ok else 1


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, required=True)
    ap.add_argument("--export", action="store_true")
    ap.add_argument("--check", type=str, default=None)
    args = ap.parse_args()
    if args.export:
        export(args.n)
        return 0
    if args.check:
        return check(args.n, args.check)
    ap.error("pass --export or --check FILE")


if __name__ == "__main__":
    raise SystemExit(main())
