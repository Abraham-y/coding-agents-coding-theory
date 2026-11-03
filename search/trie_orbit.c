/* trie_orbit.c -- orbit graph of E_q(n,d) under a group H of symbol permutations x
 * word reversal, built with the trie walk of trie_graph.c instead of all pairs.
 *
 * WHY. orbit_graph.c tests every word pair (1.4e11 pairs at q = 3, n = 12).  H acts by
 * automorphisms of the conflict graph, so it is enough to walk the neighbours of one
 * representative per orbit: orbits A, B are adjacent iff rep(A) has a neighbour in B,
 * and A is unusable (a conflict inside the orbit) iff rep(A) has a neighbour in A.
 *
 * usage: trie_orbit q n d GROUPFILE OUT_GRAPH OUT_ORBITS
 *   GROUPFILE : one element per line, q integers (image of symbol i) then 0/1 for
 *               reversal.  Must contain the identity and be closed (not checked; the
 *               driver generates it).
 *   OUT_GRAPH : mis_ils text format, vertices = usable orbits, weight = orbit size.
 *   OUT_ORBITS: one usable orbit per line, word indices (base q, MSB first).
 * Same contract as orbit_graph.c.  Nothing here is trusted; the driver expands and
 * verifies every result.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

typedef uint64_t u64;

static int q, n, t, G;
static int perm[64][4], rev[64];
static long long N;
static int32_t *orb;          /* word -> orbit id */
static long long *rep;        /* orbit id -> representative word */
static int *osize;
static char *bad;             /* orbit has an internal conflict */
static long long *stamp;      /* stamp[o] == cur_o: edge (cur_o, o) already pushed */
static long long nob;
static unsigned char w[32];
static long long cur_o;
static u64 *E; static size_t nE, capE;

static void push_edge(long long a, long long b) {
    if (a > b) { long long x = a; a = b; b = x; }
    if (nE == capE) { capE = capE ? capE * 2 : 1 << 20; E = realloc(E, capE * sizeof(u64)); if (!E) { fprintf(stderr, "oom\n"); exit(1); } }
    E[nE++] = ((u64)a << 32) | (u64)b;
}

static void walk(int k, const int *prev, long long idx) {
    if (k == n) {
        if (prev[n] <= t) {
            long long o = orb[idx];
            if (o == cur_o) { if (idx != rep[cur_o]) bad[cur_o] = 1; }
            /* the walk from rep(B) finds A whenever the walk from rep(A) finds B, so
               push each edge once, from its smaller end; raw edges were 12x the unique
               count for swap+rev at d = 5, n = 11 (5.6 GB) */
            else if (o > cur_o && stamp[o] != cur_o) { stamp[o] = cur_o; push_edge(cur_o, o); }
        }
        return;
    }
    int cur[33];
    for (int c = 0; c < q; c++) {
        int rowmin;
        cur[0] = k + 1; rowmin = cur[0];
        for (int j = 1; j <= n; j++) {
            int v = prev[j - 1] + (w[j - 1] != c);
            if (prev[j] + 1 < v) v = prev[j] + 1;
            if (cur[j - 1] + 1 < v) v = cur[j - 1] + 1;
            cur[j] = v;
            if (v < rowmin) rowmin = v;
        }
        if (rowmin > t) continue;
        walk(k + 1, cur, idx * q + c);
    }
}

static long long image(long long x, int g) {
    unsigned char s[32], u[32];
    for (int i = n - 1; i >= 0; i--) { s[i] = (unsigned char)(x % q); x /= q; }
    for (int i = 0; i < n; i++) u[i] = (unsigned char)perm[g][s[rev[g] ? n - 1 - i : i]];
    long long y = 0;
    for (int i = 0; i < n; i++) y = y * q + u[i];
    return y;
}

static int cmp64(const void *a, const void *b) {
    u64 x = *(const u64 *)a, y = *(const u64 *)b;
    return x < y ? -1 : x > y;
}

int main(int argc, char **argv) {
    if (argc != 7) { fprintf(stderr, "usage: trie_orbit q n d GROUPFILE OUT_GRAPH OUT_ORBITS\n"); return 2; }
    q = atoi(argv[1]); n = atoi(argv[2]); t = atoi(argv[3]) - 1;
    if (q < 2 || q > 4 || n < 1 || n > 16 || t < 0) { fprintf(stderr, "bad q/n/d\n"); return 2; }
    FILE *gf = fopen(argv[4], "r");
    if (!gf) { perror(argv[4]); return 2; }
    G = 0;
    while (G < 64) {
        int ok = 1;
        for (int i = 0; i < q && ok; i++) ok = fscanf(gf, "%d", &perm[G][i]) == 1;
        if (!ok || fscanf(gf, "%d", &rev[G]) != 1) break;
        G++;
    }
    fclose(gf);
    N = 1;
    for (int i = 0; i < n; i++) N *= q;
    orb = malloc(sizeof(int32_t) * N);
    for (long long x = 0; x < N; x++) orb[x] = -1;
    rep = malloc(sizeof(long long) * N); osize = calloc(N, sizeof(int)); nob = 0;
    for (long long x = 0; x < N; x++) {
        if (orb[x] >= 0) continue;
        rep[nob] = x;
        for (int g = 0; g < G; g++) {
            long long y = image(x, g);
            if (orb[y] < 0) { orb[y] = (int32_t)nob; osize[nob]++; }
        }
        nob++;
    }
    bad = calloc(nob, 1);
    stamp = malloc(sizeof(long long) * nob);
    for (long long o = 0; o < nob; o++) stamp[o] = -1;
    int row0[33];
    for (int j = 0; j <= n; j++) row0[j] = j;
    for (cur_o = 0; cur_o < nob; cur_o++) {
        long long x = rep[cur_o];
        for (int i = n - 1; i >= 0; i--) { w[i] = (unsigned char)(x % q); x /= q; }
        walk(0, row0, 0);
    }
    qsort(E, nE, sizeof(u64), cmp64);
    long long *newid = malloc(sizeof(long long) * nob), k = 0;
    for (long long o = 0; o < nob; o++) newid[o] = bad[o] ? -1 : k++;
    size_t m = 0;
    for (size_t i = 0; i < nE; i++) {
        if (i && E[i] == E[i - 1]) continue;
        long long a = (long long)(E[i] >> 32), b = (long long)(E[i] & 0xffffffffu);
        if (newid[a] >= 0 && newid[b] >= 0) E[m++] = E[i];
    }
    FILE *og = fopen(argv[5], "w"), *oo = fopen(argv[6], "w");
    if (!og || !oo) { perror("output"); return 2; }
    fprintf(og, "%lld %zu\n", k, m);
    for (long long o = 0, c = 0; o < nob; o++) if (!bad[o]) fprintf(og, ++c < k ? "%d " : "%d\n", osize[o]);
    for (size_t i = 0; i < m; i++)
        fprintf(og, "%lld %lld\n", newid[E[i] >> 32], newid[E[i] & 0xffffffffu]);
    /* orbit members: scan words once, bucket by orbit */
    long long *start = calloc(nob + 1, sizeof(long long)), *fillp = malloc(sizeof(long long) * nob);
    for (long long x = 0; x < N; x++) start[orb[x] + 1]++;
    for (long long o = 0; o < nob; o++) start[o + 1] += start[o];
    long long *mem = malloc(sizeof(long long) * N);
    memcpy(fillp, start, sizeof(long long) * nob);
    for (long long x = 0; x < N; x++) mem[fillp[orb[x]]++] = x;
    for (long long o = 0; o < nob; o++) {
        if (bad[o]) continue;
        for (long long j = start[o]; j < start[o + 1]; j++)
            fprintf(oo, j + 1 < start[o + 1] ? "%lld " : "%lld\n", mem[j]);
    }
    fclose(og); fclose(oo);
    fprintf(stderr, "words %lld orbits %lld usable %lld edges %zu (raw %zu)\n", N, nob, k, m, nE);
    return 0;
}
