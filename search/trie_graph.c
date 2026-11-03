/* trie_graph.c -- the same conflict graph as lev_graph.c (vertices = all q^n words,
 * edge iff Levenshtein distance <= d-1), built by walking the trie of target words
 * instead of testing all pairs.
 *
 * WHY. lev_graph tests every pair: 1.7e9 pairs at q = 3, n = 10, over an hour per
 * graph on the laptop.  Here, for each source word w, a depth-first walk over target
 * prefixes v[0..k) carries the DP row lev(w[0..j), v[0..k)) for j = 0..n and prunes a
 * branch as soon as the whole row exceeds t = d-1 (row minima never decrease with k).
 * The work per source word is about the number of prefixes within distance t.
 *
 * usage: trie_graph q n d > graph.txt          (mis_ils text format, as lev_graph)
 * Word index = base-q integer, most significant symbol first.  Edges a < b, each once.
 * Search-side only: tests/ compare its edge sets with lev_graph's on small cases, and
 * every code found with it is verified on its own.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int q, n, t;
static unsigned char w[32];
static long long src;
static long long m;
static FILE *out;
static int emit;

static void walk(int k, const int *prev, long long idx) {
    if (k == n) {
        if (prev[n] <= t && idx > src) {
            if (emit) fprintf(out, "%lld %lld\n", src, idx);
            m++;
        }
        return;
    }
    int cur[33];
    for (int c = 0; c < q; c++) {
        int rowmin;
        cur[0] = k + 1;
        rowmin = cur[0];
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

static void pass(void) {
    long long N = 1;
    for (int i = 0; i < n; i++) N *= q;
    int row0[33];
    for (int j = 0; j <= n; j++) row0[j] = j;
    m = 0;
    for (src = 0; src < N; src++) {
        long long x = src;
        for (int i = n - 1; i >= 0; i--) { w[i] = (unsigned char)(x % q); x /= q; }
        walk(0, row0, 0);
    }
}

int main(int argc, char **argv) {
    if (argc != 4) { fprintf(stderr, "usage: trie_graph q n d\n"); return 2; }
    q = atoi(argv[1]); n = atoi(argv[2]); t = atoi(argv[3]) - 1;
    if (q < 2 || q > 4 || n < 1 || n > 16 || t < 0) { fprintf(stderr, "bad q/n/d\n"); return 2; }
    out = stdout;
    long long N = 1;
    for (int i = 0; i < n; i++) N *= q;
    emit = 0; pass();                    /* count edges for the header */
    printf("%lld %lld\n", N, m);
    for (long long i = 0; i < N; i++) printf(i + 1 < N ? "1 " : "1\n");
    emit = 1; pass();
    return 0;
}
