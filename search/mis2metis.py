"""Convert a mis_ils text graph (k m / weights / edges) to METIS format for
ReduMIS (1-indexed, sorted adjacency, header "k m"), ignoring weights: use on
UNWEIGHTED orbit graphs (all weights equal) only, and multiply the result by
the common weight.  usage: python search/mis2metis.py in.txt out.metis"""
import sys
src, dst = sys.argv[1], sys.argv[2]
with open(src) as fh:
    k, m = map(int, fh.readline().split())
    w = list(map(int, fh.readline().split()))
    assert len(set(w)) == 1, "weighted graph; ReduMIS is unweighted"
    adj = [[] for _ in range(k)]
    for line in fh:
        a, b = map(int, line.split())
        adj[a].append(b + 1); adj[b].append(a + 1)
with open(dst, "w") as fh:
    fh.write(f"{k} {m}\n")
    for lst in adj:
        fh.write(" ".join(map(str, sorted(set(lst)))) + "\n")
print(f"{dst}: {k} vertices, {m} edges, common weight {w[0]}")
