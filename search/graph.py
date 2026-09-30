"""Conflict graph for binary s-deletion-correcting codes.

G(n, s) has vertex set {0,1}^n, with u ~ v iff u and v share a common
subsequence of length n - s. An s-deletion-correcting code is exactly an
independent set in G(n, s), so "find a big code" == "find a big independent set".

The graph is built from its natural CLIQUE COVER rather than pairwise: for each
word t of length n-s, the set of vertices containing t as a subsequence is a
clique, and every edge lies in at least one such clique. Building it this way is
O(2^n * C(n,s)) instead of O(4^n * n^2), and it hands us the clique cover for
free -- which is exactly the constraint matrix of the fractional-transversal LP
that gives the published upper bounds.

Adjacency is stored as Python big-int bitsets, one per vertex.

Nothing here is trusted. Every code this module helps find goes through
verify.py and verify2.py before it counts.
"""

from __future__ import annotations

import os
import pickle
from itertools import combinations

CACHE_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), ".graphcache")


def deletion_ball_idx(word: str, s: int) -> set[str]:
    n = len(word)
    if s == 0:
        return {word}
    return {"".join(word[i] for i in keep) for keep in combinations(range(n), n - s)}


class ConflictGraph:
    """Vertices are indices 0..2^n-1; vertex i is format(i, '0nb')."""

    def __init__(self, n: int, s: int):
        self.n = n
        self.s = s
        self.num_vertices = 1 << n
        self.words = [format(x, f"0{n}b") for x in range(self.num_vertices)]
        self.adj: list[int] = []
        self.cliques: dict[str, list[int]] = {}
        self._build()

    def _build(self):
        n, s = self.n, self.s
        buckets: dict[str, list[int]] = {}
        for idx, w in enumerate(self.words):
            for t in deletion_ball_idx(w, s):
                buckets.setdefault(t, []).append(idx)
        self.cliques = buckets

        adj = [0] * self.num_vertices
        for members in buckets.values():
            if len(members) < 2:
                continue
            mask = 0
            for m in members:
                mask |= 1 << m
            for m in members:
                adj[m] |= mask
        # drop self-loops
        for i in range(self.num_vertices):
            adj[i] &= ~(1 << i)
        self.adj = adj

    def degree(self, v: int) -> int:
        return bin(self.adj[v]).count("1")

    def is_independent(self, vertices) -> bool:
        seen = 0
        for v in vertices:
            if self.adj[v] & seen:
                return False
            seen |= 1 << v
        return True

    def to_words(self, vertices) -> list[str]:
        return [self.words[v] for v in vertices]

    def stats(self) -> dict:
        degs = [self.degree(v) for v in range(self.num_vertices)]
        return {
            "n": self.n,
            "s": self.s,
            "vertices": self.num_vertices,
            "edges": sum(degs) // 2,
            "min_degree": min(degs),
            "max_degree": max(degs),
            "mean_degree": sum(degs) / len(degs),
            "cliques": len(self.cliques),
            "max_clique_size": max(len(v) for v in self.cliques.values()),
        }


def load_graph(n: int, s: int, use_cache: bool = True) -> ConflictGraph:
    os.makedirs(CACHE_DIR, exist_ok=True)
    path = os.path.join(CACHE_DIR, f"g_{n}_{s}.pkl")
    if use_cache and os.path.exists(path):
        with open(path, "rb") as fh:
            return pickle.load(fh)
    g = ConflictGraph(n, s)
    if use_cache:
        with open(path, "wb") as fh:
            pickle.dump(g, fh, protocol=pickle.HIGHEST_PROTOCOL)
    return g


# ---------------------------------------------------------------------------
# symmetries of G(n, s)
#
# LCS is invariant under applying the SAME bijection to both words, so both of
# these are graph automorphisms:
#     complement : 0 <-> 1
#     reversal   : read the word backwards
# They commute, generating a Klein four-group. Prescribing invariance under a
# subgroup collapses the search space onto orbits (CLAUDE.md escalation ladder).
# ---------------------------------------------------------------------------


def complement_idx(v: int, n: int) -> int:
    return (~v) & ((1 << n) - 1)


def reverse_idx(v: int, n: int) -> int:
    out = 0
    for i in range(n):
        if (v >> i) & 1:
            out |= 1 << (n - 1 - i)
    return out


def orbits(n: int, group: str) -> list[tuple[int, ...]]:
    """Partition vertices into orbits under the named subgroup.

    group in {'complement', 'reversal', 'both'}
    """
    size = 1 << n
    seen = [False] * size
    out = []
    for v in range(size):
        if seen[v]:
            continue
        members = {v}
        if group in ("complement", "both"):
            members.add(complement_idx(v, n))
        if group in ("reversal", "both"):
            members.add(reverse_idx(v, n))
        if group == "both":
            members.add(complement_idx(reverse_idx(v, n), n))
        for m in members:
            seen[m] = True
        out.append(tuple(sorted(members)))
    return out


if __name__ == "__main__":
    import json
    import sys

    n = int(sys.argv[1]) if len(sys.argv) > 1 else 12
    s = int(sys.argv[2]) if len(sys.argv) > 2 else 2
    g = load_graph(n, s, use_cache=False)
    print(json.dumps(g.stats(), indent=2))
    for grp in ("complement", "reversal", "both"):
        o = orbits(n, grp)
        sizes = {}
        for orb in o:
            sizes[len(orb)] = sizes.get(len(orb), 0) + 1
        print(f"orbits under {grp:11s}: {len(o):6d}  size histogram {sizes}")
