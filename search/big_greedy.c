/* big_greedy.c -- greedy (lexicode-style) construction of quaternary edit-metric
 * codes for larger n (scoping the n >= 10 cells of Houghten's table).
 *
 * A word w is accepted iff its Levenshtein distance to every accepted codeword
 * is >= d.  Edit distance is a path metric on strings of all lengths, so with
 *   tr = floor((d-1)/2)  and  mr = d-1-tr,
 * lev(w, c) <= d-1  iff  B_tr(w) meets B_mr(c).  We keep one bitset per string
 * length marking the union of B_mr over accepted codewords, and accept w iff no
 * element of B_tr(w) is marked.  The test is exact, not a heuristic; the only
 * heuristic is the ORDER in which candidates are offered.
 *
 * Orders:  lex | rand | runs   (runs = fewest runs first, lex within)
 * Groups:  none | rev | v4 | v4rev | s4 | s4rev
 *   With a group, candidates are whole orbits (alphabet permutations, optionally
 *   with word reversal).  An orbit is accepted iff every member passes the mark
 *   test AND members are pairwise at distance >= d (checked by direct DP).
 *
 * usage: big_greedy n d order seed group  > words.txt   (summary on stderr)
 * build: cc -O3 -o search/big_greedy search/big_greedy.c
 * n <= 16 for lex order (no candidate array); n <= 14 for rand/runs.
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <time.h>

#define MAXL 20
static int N, D, TR, MR;
static uint64_t *bits[MAXL + 1];

static inline uint64_t enc(const uint8_t *s, int L) {
    uint64_t v = 0;
    for (int i = 0; i < L; i++) v = (v << 2) | s[i];
    return v;
}
static inline void dec(uint64_t v, int L, uint8_t *s) {
    for (int i = L - 1; i >= 0; i--) { s[i] = v & 3; v >>= 2; }
}

typedef int (*cb_t)(int L, uint64_t v);
static int cb_test(int L, uint64_t v) { return (int)((bits[L][v >> 6] >> (v & 63)) & 1); }
static int cb_mark(int L, uint64_t v) { bits[L][v >> 6] |= 1ULL << (v & 63); return 0; }

/* Calls cb on every string within r edits of s (with repeats).  Returns 1 as soon
 * as cb returns 1. */
static int ball(const uint8_t *s, int L, int r, cb_t cb) {
    if (r == 0) return cb(L, enc(s, L));
    uint8_t t[MAXL + 2];
    if (ball(s, L, r - 1, cb)) return 1;
    for (int i = 0; i < L; i++) {                       /* deletions */
        int k = 0;
        for (int j = 0; j < L; j++) if (j != i) t[k++] = s[j];
        if (L - 1 >= 1 && ball(t, L - 1, r - 1, cb)) return 1;
    }
    for (int i = 0; i < L; i++) {                       /* substitutions */
        memcpy(t, s, L);
        for (uint8_t c = 0; c < 4; c++) if (c != s[i]) {
            t[i] = c;
            if (ball(t, L, r - 1, cb)) return 1;
        }
    }
    for (int i = 0; i <= L; i++)                        /* insertions */
        for (uint8_t c = 0; c < 4; c++) {
            int k = 0;
            for (int j = 0; j < i; j++) t[k++] = s[j];
            t[k++] = c;
            for (int j = i; j < L; j++) t[k++] = s[j];
            if (ball(t, L + 1, r - 1, cb)) return 1;
        }
    return 0;
}

static int lev(const uint8_t *a, const uint8_t *b, int n) {
    int prev[MAXL + 1], cur[MAXL + 1];
    for (int j = 0; j <= n; j++) prev[j] = j;
    for (int i = 1; i <= n; i++) {
        cur[0] = i;
        for (int j = 1; j <= n; j++) {
            int x = prev[j] + 1, y = cur[j - 1] + 1, z = prev[j - 1] + (a[i - 1] != b[j - 1]);
            cur[j] = x < y ? (x < z ? x : z) : (y < z ? y : z);
        }
        memcpy(prev, cur, sizeof(int) * (n + 1));
    }
    return prev[n];
}

/* group: list of symbol permutations, each optionally combined with reversal */
static int G;                        /* number of group elements */
static uint8_t gperm[96][4];
static uint8_t grev[96];

static void add_elem(const uint8_t p[4], int rev) {
    memcpy(gperm[G], p, 4); grev[G] = (uint8_t)rev; G++;
}
static void build_group(const char *name) {
    static const uint8_t id[4] = {0, 1, 2, 3};
    static const uint8_t v4[4][4] = {{0,1,2,3},{1,0,3,2},{2,3,0,1},{3,2,1,0}};
    int withrev = strstr(name, "rev") != NULL;
    G = 0;
    if (!strncmp(name, "none", 4) || !strcmp(name, "rev")) {
        add_elem(id, 0);
    } else if (!strncmp(name, "v4", 2)) {
        for (int i = 0; i < 4; i++) add_elem(v4[i], 0);
    } else if (!strncmp(name, "s4", 2)) {
        uint8_t p[4];
        for (p[0] = 0; p[0] < 4; p[0]++) for (p[1] = 0; p[1] < 4; p[1]++)
        for (p[2] = 0; p[2] < 4; p[2]++) for (p[3] = 0; p[3] < 4; p[3]++) {
            if (p[0]==p[1]||p[0]==p[2]||p[0]==p[3]||p[1]==p[2]||p[1]==p[3]||p[2]==p[3]) continue;
            add_elem(p, 0);
        }
    } else { fprintf(stderr, "unknown group %s\n", name); exit(2); }
    if (withrev) { int g0 = G; for (int i = 0; i < g0; i++) add_elem(gperm[i], 1); }
}
static uint64_t apply(int g, uint64_t v) {
    uint8_t s[MAXL], t[MAXL];
    dec(v, N, s);
    for (int i = 0; i < N; i++) t[grev[g] ? N - 1 - i : i] = gperm[g][s[i]];
    return enc(t, N);
}

static uint64_t rng_state;
static uint64_t splitmix(void) {
    uint64_t z = (rng_state += 0x9E3779B97F4A7C15ULL);
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    return z ^ (z >> 31);
}

int main(int argc, char **argv) {
    if (argc != 6) {
        fprintf(stderr, "usage: big_greedy n d order(lex|rand|runs) seed group\n");
        return 2;
    }
    N = atoi(argv[1]); D = atoi(argv[2]);
    const char *order = argv[3];
    uint64_t seed = strtoull(argv[4], NULL, 10);
    const char *gname = argv[5];
    if (N < 2 || N > 16 || D < 2 || D > 7) { fprintf(stderr, "bad n or d\n"); return 2; }
    if (strcmp(order, "lex") && N > 14) {
        fprintf(stderr, "orders other than lex need a 4^n candidate array; n <= 14 only\n");
        return 2;
    }
    TR = (D - 1) / 2; MR = D - 1 - TR;
    if (N + MR > MAXL) { fprintf(stderr, "n + mark radius too large\n"); return 2; }
    for (int L = (N - MR < 1 ? 1 : N - MR); L <= N + MR; L++) {
        uint64_t words = (1ULL << (2 * L)) / 64 + 1;
        bits[L] = calloc(words, sizeof(uint64_t));
        if (!bits[L]) { fprintf(stderr, "out of memory at L=%d\n", L); return 2; }
    }
    build_group(gname);

    uint64_t total = 1ULL << (2 * N);
    /* lex order needs no candidate array: candidate idx is the word idx itself */
    uint32_t *cand = NULL;
    if (strcmp(order, "lex")) {
        cand = malloc(sizeof(uint32_t) * total);
        if (!cand) { fprintf(stderr, "out of memory (candidates)\n"); return 2; }
        for (uint64_t v = 0; v < total; v++) cand[v] = (uint32_t)v;
    }
    if (!strcmp(order, "rand")) {
        rng_state = seed;
        for (uint64_t i = total - 1; i > 0; i--) {
            uint64_t j = splitmix() % (i + 1);
            uint32_t t = cand[i]; cand[i] = cand[j]; cand[j] = t;
        }
    } else if (!strcmp(order, "runs")) {
        /* stable counting sort by number of runs (1..N), lex within */
        uint64_t cnt[MAXL + 2] = {0};
        uint8_t s[MAXL];
        for (uint64_t v = 0; v < total; v++) {
            dec(v, N, s); int r = 1;
            for (int i = 1; i < N; i++) r += s[i] != s[i - 1];
            cnt[r + 1]++;
        }
        for (int r = 1; r <= N + 1; r++) cnt[r] += cnt[r - 1];
        for (uint64_t v = 0; v < total; v++) {
            dec(v, N, s); int r = 1;
            for (int i = 1; i < N; i++) r += s[i] != s[i - 1];
            cand[cnt[r]++] = (uint32_t)v;
        }
    } else if (strcmp(order, "lex")) { fprintf(stderr, "unknown order\n"); return 2; }

    clock_t t0 = clock();
    uint64_t size = 0, orbits = 0;
    uint8_t members[96][MAXL];
    uint64_t orb[96];
    for (uint64_t idx = 0; idx < total; idx++) {
        uint64_t v = cand ? cand[idx] : idx;
        int m = 0, rep = 1;
        for (int g = 0; g < G && rep; g++) {
            uint64_t u = apply(g, v);
            if (u < v) rep = 0;                       /* only orbit representatives */
            int dup = 0;
            for (int k = 0; k < m; k++) if (orb[k] == u) { dup = 1; break; }
            if (!dup) orb[m++] = u;
        }
        if (!rep) continue;
        int ok = 1;
        for (int k = 0; k < m && ok; k++) {
            dec(orb[k], N, members[k]);
            if (ball(members[k], N, TR, cb_test)) ok = 0;
        }
        for (int a = 0; a < m && ok; a++)
            for (int b = a + 1; b < m && ok; b++)
                if (lev(members[a], members[b], N) < D) ok = 0;
        if (!ok) continue;
        for (int k = 0; k < m; k++) {
            ball(members[k], N, MR, cb_mark);
            for (int i = 0; i < N; i++) putchar('0' + members[k][i]);
            putchar('\n');
        }
        size += m; orbits++;
    }
    double secs = (double)(clock() - t0) / CLOCKS_PER_SEC;
    fprintf(stderr, "n=%d d=%d order=%s seed=%llu group=%s |G|=%d size=%llu orbits=%llu cpu_s=%.1f\n",
            N, D, order, (unsigned long long)seed, gname, G,
            (unsigned long long)size, (unsigned long long)orbits, secs);
    return 0;
}
