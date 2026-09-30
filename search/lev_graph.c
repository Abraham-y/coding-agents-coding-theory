/* lev_graph.c -- conflict graph of the fixed-length q-ary edit-metric code
 * problem E_q(n, d): vertices are all q^n words, an edge joins two words whose
 * Levenshtein distance is <= d-1 (they cannot both be in a code of minimum
 * edit distance d).  Output is the mis_ils text format:
 *     k m
 *     1 1 ... 1
 *     a b   (m lines)
 * Word index = base-q integer, most significant symbol first, which is the
 * order of itertools.product in the Python side.
 *
 * usage: lev_graph q n d > graph.txt
 *
 * Distance test: banded Wagner-Fischer with threshold t = d-1 and early exit
 * when every entry of a row exceeds t.  Plain, no bit tricks, on purpose.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int q, n, t;

static int lev_le(const unsigned char *a, const unsigned char *b) {
    /* returns 1 iff lev(a,b) <= t; both length n */
    int prev[32], cur[32];
    for (int j = 0; j <= n; j++) prev[j] = j;
    for (int i = 1; i <= n; i++) {
        int lo = i - t < 1 ? 1 : i - t;
        int hi = i + t > n ? n : i + t;
        int rowmin = 1 << 20;
        cur[0] = i;
        for (int j = 0; j < lo; j++) cur[j] = (1 << 20);
        cur[0] = i;
        for (int j = lo; j <= hi; j++) {
            int c = prev[j - 1] + (a[i - 1] != b[j - 1]);
            int del = prev[j] + 1, ins = cur[j - 1] + 1;
            if (del < c) c = del;
            if (ins < c) c = ins;
            cur[j] = c;
            if (c < rowmin) rowmin = c;
        }
        for (int j = hi + 1; j <= n; j++) cur[j] = (1 << 20);
        if (rowmin > t) return 0;
        memcpy(prev, cur, sizeof(int) * (n + 1));
    }
    return prev[n] <= t;
}

int main(int argc, char **argv) {
    if (argc < 4) { fprintf(stderr, "usage: %s q n d\n", argv[0]); return 2; }
    q = atoi(argv[1]); n = atoi(argv[2]); int d = atoi(argv[3]); t = d - 1;
    if (n > 30 || q > 16) return 2;
    long k = 1; for (int i = 0; i < n; i++) k *= q;
    unsigned char *words = malloc((size_t)k * n);
    for (long idx = 0; idx < k; idx++) {
        long x = idx;
        for (int i = n - 1; i >= 0; i--) { words[idx * n + i] = (unsigned char)(x % q); x /= q; }
    }
    /* two passes: count, then print (so the header is exact) */
    long m = 0;
    for (long a = 0; a < k; a++)
        for (long b = a + 1; b < k; b++)
            if (lev_le(words + a * n, words + b * n)) m++;
    printf("%ld %ld\n", k, m);
    for (long i = 0; i < k; i++) printf(i + 1 < k ? "1 " : "1\n");
    for (long a = 0; a < k; a++)
        for (long b = a + 1; b < k; b++)
            if (lev_le(words + a * n, words + b * n)) printf("%ld %ld\n", a, b);
    return 0;
}
