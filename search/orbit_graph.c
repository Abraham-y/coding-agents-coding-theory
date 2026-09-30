/* orbit_graph.c -- build the orbit graph of E_q(n,d) under a subgroup H of
 * S_q x C_2 (symbol permutations x word reversal), in C, so that n = 8
 * (65,536 words) and beyond are reachable.  The Python builder
 * (search/edit4.py:orbit_problem) is O(k^2) over Python big-int bitmasks and
 * takes minutes at n = 7; here the conflict test is pairwise banded
 * Levenshtein over all word pairs (2.1e9 pairs at n = 8, a few minutes) and
 * orbit edges are deduplicated with a bitset over orbit pairs.
 *
 * usage: orbit_graph q n d GROUPFILE OUT_GRAPH OUT_ORBITS
 *   GROUPFILE: one element per line: q integers (the permutation image of
 *              symbol i) then 0/1 for reversal, e.g. "1 0 3 2 0".
 *              Must contain the identity and be closed (not checked).
 *   OUT_GRAPH : mis_ils text format (k m / weights / edges), vertices = USABLE
 *               orbits (no internal conflict), weight = orbit size.
 *   OUT_ORBITS: one usable orbit per line, word indices (base-q, MSB first)
 *               space-separated, in the same order as OUT_GRAPH vertices.
 * The definition of conflict is lev(u,v) <= d-1, the same as lev_graph.c.
 * Nothing here is trusted: the driver expands and verifies every result.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

typedef uint64_t u64;
static int q, n, t;

static int lev_le(const unsigned char *a, const unsigned char *b) {
    int prev[32], cur[32];
    for (int j = 0; j <= n; j++) prev[j] = j;
    for (int i = 1; i <= n; i++) {
        int lo = i - t < 1 ? 1 : i - t, hi = i + t > n ? n : i + t, rowmin = 1 << 20;
        for (int j = 0; j < lo; j++) cur[j] = 1 << 20;
        cur[0] = i;
        for (int j = lo; j <= hi; j++) {
            int c = prev[j - 1] + (a[i - 1] != b[j - 1]);
            if (prev[j] + 1 < c) c = prev[j] + 1;
            if (cur[j - 1] + 1 < c) c = cur[j - 1] + 1;
            cur[j] = c; if (c < rowmin) rowmin = c;
        }
        for (int j = hi + 1; j <= n; j++) cur[j] = 1 << 20;
        if (rowmin > t) return 0;
        memcpy(prev, cur, sizeof(int) * (n + 1));
    }
    return prev[n] <= t;
}

int main(int argc, char **argv) {
    if (argc < 7) { fprintf(stderr, "usage: %s q n d GROUPFILE OUT_GRAPH OUT_ORBITS\n", argv[0]); return 2; }
    q = atoi(argv[1]); n = atoi(argv[2]); t = atoi(argv[3]) - 1;
    long k = 1; for (int i = 0; i < n; i++) k *= q;
    /* group */
    int G[64][17]; int ng = 0;
    FILE *gf = fopen(argv[4], "r"); if (!gf) { perror("group"); return 2; }
    while (ng < 64) {
        int ok = 1;
        for (int i = 0; i < q; i++) if (fscanf(gf, "%d", &G[ng][i]) != 1) { ok = 0; break; }
        if (!ok) break;
        if (fscanf(gf, "%d", &G[ng][16]) != 1) break;
        ng++;
    }
    fclose(gf);
    unsigned char *words = malloc((size_t)k * n);
    for (long idx = 0; idx < k; idx++) { long x = idx; for (int i = n - 1; i >= 0; i--) { words[idx * n + i] = (unsigned char)(x % q); x /= q; } }
    /* orbits: orbit id per word = index of its minimal image */
    long *orb = malloc(sizeof(long) * k); for (long i = 0; i < k; i++) orb[i] = -1;
    long norb = 0; long *rep = malloc(sizeof(long) * k); int *osize = malloc(sizeof(int) * k);
    unsigned char img[32];
    long *members = malloc(sizeof(long) * k); long *mstart = malloc(sizeof(long) * (k + 1));
    long mfill = 0;
    for (long idx = 0; idx < k; idx++) {
        if (orb[idx] >= 0) continue;
        long ids[64]; int cnt = 0;
        for (int g = 0; g < ng; g++) {
            for (int i = 0; i < n; i++) img[i] = (unsigned char)G[g][words[idx * n + i]];
            if (G[g][16]) for (int i = 0; i < n / 2; i++) { unsigned char c = img[i]; img[i] = img[n - 1 - i]; img[n - 1 - i] = c; }
            long v = 0; for (int i = 0; i < n; i++) v = v * q + img[i];
            int dup = 0; for (int j = 0; j < cnt; j++) if (ids[j] == v) { dup = 1; break; }
            if (!dup) ids[cnt++] = v;
        }
        rep[norb] = idx; osize[norb] = cnt; mstart[norb] = mfill;
        for (int j = 0; j < cnt; j++) { orb[ids[j]] = norb; members[mfill++] = ids[j]; }
        norb++;
    }
    mstart[norb] = mfill;
    fprintf(stderr, "%ld words, %ld orbits under |H|=%d\n", k, norb, ng);
    /* usable: no internal conflict */
    char *usable = malloc(norb); memset(usable, 1, norb);
    for (long o = 0; o < norb; o++)
        for (long a = mstart[o]; a < mstart[o + 1] && usable[o]; a++)
            for (long b = a + 1; b < mstart[o + 1]; b++)
                if (lev_le(words + members[a] * n, words + members[b] * n)) { usable[o] = 0; break; }
    long *newid = malloc(sizeof(long) * norb); long nu = 0;
    for (long o = 0; o < norb; o++) newid[o] = usable[o] ? nu++ : -1;
    fprintf(stderr, "%ld usable orbits\n", nu);
    /* orbit adjacency as a bitset over usable-orbit pairs (upper triangle) */
    size_t nbits = (size_t)nu * (size_t)nu;
    u64 *bits = calloc((nbits + 63) / 64, sizeof(u64));
    if (!bits) { fprintf(stderr, "too many orbits for a bitset (%ld)\n", nu); return 3; }
    long m = 0;
    /* only pairs (u,v) with u a representative: an orbit edge exists iff some
       member of orbit A conflicts with some member of B, and by invariance it
       suffices to test all members of A against ONE representative... no: it
       suffices to test one representative of A against all members of B. */
    for (long o = 0; o < norb; o++) {
        if (!usable[o]) continue;
        long a = rep[o]; long ia = newid[o];
        const unsigned char *wa = words + a * n;
        for (long v = 0; v < k; v++) {
            long ov = orb[v]; if (ov == o || !usable[ov]) continue;
            long iv = newid[ov];
            size_t bit = (size_t)ia * nu + iv;
            if ((bits[bit >> 6] >> (bit & 63)) & 1) continue;
            if (lev_le(wa, words + v * n)) {
                bits[bit >> 6] |= 1ULL << (bit & 63);
                size_t bit2 = (size_t)iv * nu + ia;
                if (!((bits[bit2 >> 6] >> (bit2 & 63)) & 1)) { bits[bit2 >> 6] |= 1ULL << (bit2 & 63); }
            }
        }
        if ((o & 1023) == 0) fprintf(stderr, "\r%ld/%ld orbits", o, norb);
    }
    for (long i = 0; i < nu; i++) for (long j = i + 1; j < nu; j++) { size_t bit = (size_t)i * nu + j; if ((bits[bit >> 6] >> (bit & 63)) & 1) m++; }
    fprintf(stderr, "\n%ld orbit edges\n", m);
    FILE *fo = fopen(argv[5], "w");
    fprintf(fo, "%ld %ld\n", nu, m);
    for (long o = 0; o < norb; o++) if (usable[o]) fprintf(fo, "%d ", osize[o]);
    fprintf(fo, "\n");
    for (long i = 0; i < nu; i++) for (long j = i + 1; j < nu; j++) { size_t bit = (size_t)i * nu + j; if ((bits[bit >> 6] >> (bit & 63)) & 1) fprintf(fo, "%ld %ld\n", i, j); }
    fclose(fo);
    FILE *fr = fopen(argv[6], "w");
    for (long o = 0; o < norb; o++) {
        if (!usable[o]) continue;
        for (long a = mstart[o]; a < mstart[o + 1]; a++) fprintf(fr, "%ld%s", members[a], a + 1 < mstart[o + 1] ? " " : "\n");
    }
    fclose(fr);
    return 0;
}
