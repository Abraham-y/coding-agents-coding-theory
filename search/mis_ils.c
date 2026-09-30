/* mis_ils.c -- ARW-style iterated local search for maximum (weight) independent
 * set, written for the orbit graphs of quaternary edit-metric codes.
 *
 * WHY C. The repo's Python WeightedSearcher (search/edit4_orbit_ls.py) tops out
 * at 116 on the free-acting order-4 orbit graphs at n=6, and a stronger
 * operator (2-improvements + plateau moves) found 120 there. The move set is
 * not exotic; what matters is iterations per second. This file does the same
 * ARW moves (Andrade, Resende, Werneck 2012) with bitset adjacency and O(1)
 * tightness updates, and runs three orders of magnitude more iterations in the
 * same wall time.
 *
 * Moves (weighted; the unweighted case is the same with all weights equal):
 *   maximize     : add free vertices (no selected neighbour), heaviest of a
 *                  random sample first.
 *   (1,1)-swap   : x out, a 1-tight neighbour v in, when w[v] > w[x].
 *   (1,2)-swap   : x out, two NON-adjacent 1-tight neighbours in, when their
 *                  weight sum exceeds w[x].  This is the "2-improvement".
 *   (omega,1)    : v in, all its selected neighbours out, when w[v] exceeds
 *                  their total weight.  Only fires on weighted instances.
 *   plateau      : x out, an equal-weight 1-tight neighbour in.
 * Perturbation   : force-insert k random non-solution vertices (k = 1 + the
 *                  number of consecutive coin heads), pin them for the next
 *                  local search, and accept the result ARW-style (always if not
 *                  worse; else with probability 1/(1 + delta * delta_best)).
 *
 * Input graph format (plain text):
 *   k m
 *   w_0 w_1 ... w_{k-1}
 *   a b            (m lines, 0-indexed, undirected, no duplicates needed)
 * Usage:
 *   mis_ils graph.txt SEED SECONDS [--target T] [--init FILE] [--quiet] [--restart N]
 * Output: progress lines "improved <t> <weight> <count>" and one final line
 *   BEST <weight> <count> : v1 v2 ...
 * The caller expands the vertex list to codewords and verifies it; nothing
 * printed here is trusted on its own.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <time.h>

typedef uint64_t u64;

static int n, m, W;
static int *deg, *astart, *adj;
static u64 *bits;
static long *w;
static char *inS, *pinned;
static int *tight;
static long *tightw;
static int *S, *posS, nS;
static long curw;
static int *F, *posF, nF;
static int weighted;
static u64 rng_state;

static inline u64 rnd(void) {
    u64 x = rng_state;
    x ^= x << 13; x ^= x >> 7; x ^= x << 17;
    return rng_state = x;
}
static inline int rnd_int(int k) { return (int)(rnd() % (u64)k); }
static inline int isadj(int a, int b) {
    return (int)((bits[(size_t)a * W + (b >> 6)] >> (b & 63)) & 1ULL);
}
static double now(void) {
    struct timespec ts; clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + ts.tv_nsec * 1e-9;
}

static inline void free_add(int v) { posF[v] = nF; F[nF++] = v; }
static inline void free_del(int v) {
    int p = posF[v], last = F[--nF];
    F[p] = last; posF[last] = p; posF[v] = -1;
}

static void add(int v) {
    inS[v] = 1; posS[v] = nS; S[nS++] = v; curw += w[v];
    if (posF[v] >= 0) free_del(v);
    long wv = w[v];
    for (int i = astart[v]; i < astart[v + 1]; i++) {
        int u = adj[i];
        tight[u]++; tightw[u] += wv;
        if (tight[u] == 1 && posF[u] >= 0) free_del(u);
    }
}

static void removev(int v) {
    inS[v] = 0;
    int p = posS[v], last = S[--nS];
    S[p] = last; posS[last] = p; posS[v] = -1;
    curw -= w[v];
    free_add(v);   /* S independent => tight[v] == 0 */
    long wv = w[v];
    for (int i = astart[v]; i < astart[v + 1]; i++) {
        int u = adj[i];
        tight[u]--; tightw[u] -= wv;
        if (tight[u] == 0 && !inS[u]) free_add(u);
    }
}

static void maximize(void) {
    while (nF > 0) {
        int v;
        if (!weighted || nF <= 8) {
            v = F[rnd_int(nF)];
            if (weighted) {
                for (int i = 0; i < nF; i++) if (w[F[i]] > w[v]) v = F[i];
            }
        } else {
            v = F[rnd_int(nF)];
            for (int t = 1; t < 8; t++) { int u = F[rnd_int(nF)]; if (w[u] > w[v]) v = u; }
        }
        add(v);
    }
}

static int *tmpS, *L;

static int improve(void) {
    int improved = 0;
    int cnt0 = nS;
    memcpy(tmpS, S, sizeof(int) * nS);
    for (int i = cnt0 - 1; i > 0; i--) { int j = rnd_int(i + 1); int t = tmpS[i]; tmpS[i] = tmpS[j]; tmpS[j] = t; }
    for (int i = 0; i < cnt0; i++) {
        int x = tmpS[i];
        if (!inS[x] || pinned[x]) continue;
        int cnt = 0;
        for (int j = astart[x]; j < astart[x + 1]; j++) {
            int u = adj[j];
            if (!inS[u] && tight[u] == 1) L[cnt++] = u;
        }
        if (cnt == 0) continue;
        if (weighted) {
            int best = -1;
            for (int a = 0; a < cnt; a++) if (w[L[a]] > w[x] && (best < 0 || w[L[a]] > w[best])) best = L[a];
            if (best >= 0) { removev(x); add(best); maximize(); improved = 1; continue; }
        }
        /* (1,2): random start so different pairs get tried on repeated passes */
        int found = 0, fa = -1, fb = -1;
        int off = rnd_int(cnt);
        for (int a0 = 0; a0 < cnt && !found; a0++) {
            int a = L[(a0 + off) % cnt];
            if (w[a] + w[L[0]] <= 0) continue;
            for (int b0 = a0 + 1; b0 < cnt; b0++) {
                int b = L[(b0 + off) % cnt];
                if (w[a] + w[b] > w[x] && !isadj(a, b)) { found = 1; fa = a; fb = b; break; }
            }
        }
        if (found) { removev(x); add(fa); add(fb); maximize(); improved = 1; }
    }
    if (weighted) {
        for (int v = 0; v < n; v++) {
            if (!inS[v] && tight[v] > 0 && tightw[v] < w[v]) {
                for (int j = astart[v]; j < astart[v + 1]; j++) { int u = adj[j]; if (inS[u]) removev(u); }
                add(v); maximize(); improved = 1;
            }
        }
    }
    return improved;
}

static void local_search(void) { maximize(); while (improve()) ; }

static void plateau_move(void) {
    if (nS == 0) return;
    int x = S[rnd_int(nS)];
    if (pinned[x]) return;
    int cnt = 0;
    for (int j = astart[x]; j < astart[x + 1]; j++) {
        int u = adj[j];
        if (!inS[u] && tight[u] == 1 && w[u] == w[x]) L[cnt++] = u;
    }
    if (cnt == 0) return;
    int v = L[rnd_int(cnt)];
    removev(x); add(v); maximize();
}

static void clear_all(void) {
    while (nS > 0) removev(S[nS - 1]);
}
static void load(const int *sol, int cnt) {
    clear_all();
    for (int i = 0; i < cnt; i++) if (tight[sol[i]] == 0 && !inS[sol[i]]) add(sol[i]);
}

int main(int argc, char **argv) {
    if (argc < 4) { fprintf(stderr, "usage: %s graph seed seconds [--target T] [--init FILE] [--quiet]\n", argv[0]); return 2; }
    const char *gpath = argv[1];
    u64 seed = strtoull(argv[2], NULL, 10);
    double secs = atof(argv[3]);
    long target = -1; const char *init = NULL; int quiet = 0; long restart_after = 20000;
    for (int i = 4; i < argc; i++) {
        if (!strcmp(argv[i], "--target") && i + 1 < argc) target = atol(argv[++i]);
        else if (!strcmp(argv[i], "--init") && i + 1 < argc) init = argv[++i];
        else if (!strcmp(argv[i], "--quiet")) quiet = 1;
        else if (!strcmp(argv[i], "--restart") && i + 1 < argc) restart_after = atol(argv[++i]);
    }
    rng_state = seed * 0x9E3779B97F4A7C15ULL + 0x1234567ULL;
    for (int i = 0; i < 20; i++) rnd();

    FILE *f = fopen(gpath, "r");
    if (!f) { perror("graph"); return 2; }
    if (fscanf(f, "%d %d", &n, &m) != 2) return 2;
    w = malloc(sizeof(long) * n);
    weighted = 0;
    for (int i = 0; i < n; i++) { if (fscanf(f, "%ld", &w[i]) != 1) return 2; if (w[i] != w[0]) weighted = 1; }
    int *ea = malloc(sizeof(int) * m), *eb = malloc(sizeof(int) * m);
    deg = calloc(n, sizeof(int));
    W = (n + 63) / 64;
    bits = calloc((size_t)n * W, sizeof(u64));
    int mm = 0;
    for (int i = 0; i < m; i++) {
        int a, b; if (fscanf(f, "%d %d", &a, &b) != 2) return 2;
        if (a == b || isadj(a, b)) continue;
        bits[(size_t)a * W + (b >> 6)] |= 1ULL << (b & 63);
        bits[(size_t)b * W + (a >> 6)] |= 1ULL << (a & 63);
        ea[mm] = a; eb[mm] = b; mm++; deg[a]++; deg[b]++;
    }
    fclose(f);
    m = mm;
    astart = malloc(sizeof(int) * (n + 1));
    astart[0] = 0; for (int i = 0; i < n; i++) astart[i + 1] = astart[i] + deg[i];
    adj = malloc(sizeof(int) * (2 * (size_t)m + 1));
    int *fill = calloc(n, sizeof(int));
    for (int i = 0; i < m; i++) { adj[astart[ea[i]] + fill[ea[i]]++] = eb[i]; adj[astart[eb[i]] + fill[eb[i]]++] = ea[i]; }
    free(fill); free(ea); free(eb);

    inS = calloc(n, 1); pinned = calloc(n, 1);
    tight = calloc(n, sizeof(int)); tightw = calloc(n, sizeof(long));
    S = malloc(sizeof(int) * n); posS = malloc(sizeof(int) * n);
    F = malloc(sizeof(int) * n); posF = malloc(sizeof(int) * n);
    tmpS = malloc(sizeof(int) * n); L = malloc(sizeof(int) * n);
    nS = 0; nF = 0; curw = 0;
    for (int v = 0; v < n; v++) { posS[v] = -1; posF[v] = -1; }
    for (int v = 0; v < n; v++) free_add(v);

    if (init) {
        FILE *fi = fopen(init, "r");
        if (!fi) { perror("init"); return 2; }
        int v; while (fscanf(fi, "%d", &v) == 1) if (v >= 0 && v < n && tight[v] == 0 && !inS[v]) add(v);
        fclose(fi);
    }
    double t0 = now();
    local_search();
    int *best = malloc(sizeof(int) * n), nbest = nS; long bestw = curw;
    memcpy(best, S, sizeof(int) * nS);
    int *prev = malloc(sizeof(int) * n), nprev; long prevw;
    if (!quiet) printf("improved %.2f %ld %d\n", now() - t0, bestw, nbest);
    fflush(stdout);
    long iters = 0, stall = 0;
    while (1) {
        iters++;
        if ((iters & 63) == 0) { if (now() - t0 > secs) break; }
        if (target > 0 && bestw >= target) break;
        nprev = nS; prevw = curw; memcpy(prev, S, sizeof(int) * nS);
        /* plateau walk */
        int pm = rnd_int(4);
        for (int i = 0; i < pm; i++) plateau_move();
        /* perturbation */
        int k = 1; while ((rnd() & 1) && k < 8) k++;
        int npin = 0; int pins[8];
        for (int j = 0; j < k; j++) {
            int v, tries = 0;
            do { v = rnd_int(n); tries++; } while (inS[v] && tries < 50);
            if (inS[v]) continue;
            for (int i = astart[v]; i < astart[v + 1]; i++) { int u = adj[i]; if (inS[u]) removev(u); }
            add(v); pinned[v] = 1; pins[npin++] = v;
        }
        local_search();
        for (int j = 0; j < npin; j++) pinned[pins[j]] = 0;
        local_search();
        if (curw > bestw) {
            bestw = curw; nbest = nS; memcpy(best, S, sizeof(int) * nS); stall = 0;
            if (!quiet) { printf("improved %.2f %ld %d\n", now() - t0, bestw, nbest); fflush(stdout); }
        } else stall++;
        if (curw < prevw) {
            long d = prevw - curw, db = bestw - curw;
            /* accept worse with probability 1/(1 + d*db) */
            u64 r = rnd() % 1000000ULL;
            double p = 1.0 / (1.0 + (double)d * (double)db);
            if (r >= (u64)(p * 1000000.0)) load(prev, nprev);
        }
        if (stall > restart_after) {
            /* diversify: half the time return to the incumbent, half the time
             * start over from a fresh random maximal solution. Long runs that
             * only ever return to the incumbent sit in one basin (observed:
             * 10-minute seeds worse than 1-minute seeds on the same graph). */
            if (rnd_int(2) == 0) load(best, nbest);
            else { clear_all(); local_search(); }
            stall = 0;
        }
    }
    printf("iters %ld time %.1f\n", iters, now() - t0);
    printf("BEST %ld %d :", bestw, nbest);
    for (int i = 0; i < nbest; i++) printf(" %d", best[i]);
    printf("\n");
    return 0;
}
