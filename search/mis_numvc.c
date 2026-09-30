/* mis_numvc.c -- a second, deliberately different local-search operator for
 * maximum independent set on the orbit graphs: NuMVC (Cai, Su, Luo, Sattar,
 * JAIR 2013), i.e. minimum vertex cover local search with edge weighting and
 * two-stage configuration checking.  An independent set is the complement of
 * a vertex cover.  UNWEIGHTED only (the free-acting classes, where every orbit
 * has the same size); a weighted instance is rejected.
 *
 * Outline (standard NuMVC):
 *   - keep a cover C; "uncovered" edges are those with neither endpoint in C;
 *   - dscore[v] = change in total uncovered-edge weight if v is flipped;
 *   - each step: remove the vertex of C with highest dscore (i.e. least
 *     harmful), then add an endpoint of a random uncovered edge chosen by
 *     dscore with configuration checking (a vertex may re-enter only if a
 *     neighbour changed since it left), tie-break by age;
 *   - every step the weights of uncovered edges grow by 1; weights are scaled
 *     down periodically (forgetting);
 *   - whenever C covers everything, record |V|-|C| as an independent set and
 *     shrink the target by removing one more vertex.
 *
 * Same input format and output line as mis_ils.c:
 *   mis_numvc graph.txt SEED SECONDS [--target T] [--quiet]
 *   BEST <weight> <count> : v1 v2 ...
 * The driver expands and verifies; nothing here is trusted.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <time.h>

typedef uint64_t u64;
static int n, m;
static int *astart, *adj, *eadj;   /* eadj: edge id parallel to adj */
static int *eu, *ev;               /* edge endpoints */
static long *ew;                   /* edge weights */
static char *inC, *conf;
static long *dscore;
static long *age_;                 /* last step vertex was flipped */
static int *uncov, *uncov_pos, nuncov;
static int *Clist, *Cpos, nC;
static u64 rng_state;
static long wcommon = 1;

static inline u64 rnd(void) { u64 x = rng_state; x ^= x << 13; x ^= x >> 7; x ^= x << 17; return rng_state = x; }
static inline int rnd_int(int k) { return (int)(rnd() % (u64)k); }
static double now(void) { struct timespec ts; clock_gettime(CLOCK_MONOTONIC, &ts); return ts.tv_sec + ts.tv_nsec * 1e-9; }

static inline void uncov_add(int e) { uncov_pos[e] = nuncov; uncov[nuncov++] = e; }
static inline void uncov_del(int e) { int p = uncov_pos[e], last = uncov[--nuncov]; uncov[p] = last; uncov_pos[last] = p; uncov_pos[e] = -1; }

static void add_to_cover(int v, long step) {
    inC[v] = 1; Cpos[v] = nC; Clist[nC++] = v; age_[v] = step;
    dscore[v] = -dscore[v];
    for (int i = astart[v]; i < astart[v + 1]; i++) {
        int u = adj[i], e = eadj[i];
        conf[u] = 1;
        if (!inC[u]) { uncov_del(e); dscore[u] -= ew[e]; }
        else dscore[u] += ew[e];
    }
}
static void remove_from_cover(int v, long step) {
    inC[v] = 0; int p = Cpos[v], last = Clist[--nC]; Clist[p] = last; Cpos[last] = p; Cpos[v] = -1;
    age_[v] = step; conf[v] = 0;
    dscore[v] = -dscore[v];
    for (int i = astart[v]; i < astart[v + 1]; i++) {
        int u = adj[i], e = eadj[i];
        conf[u] = 1;
        if (!inC[u]) { uncov_add(e); dscore[u] += ew[e]; }
        else dscore[u] -= ew[e];
    }
}

int main(int argc, char **argv) {
    if (argc < 4) { fprintf(stderr, "usage: %s graph seed seconds [--target T] [--quiet]\n", argv[0]); return 2; }
    u64 seed = strtoull(argv[2], NULL, 10); double secs = atof(argv[3]);
    long target = -1; int quiet = 0;
    for (int i = 4; i < argc; i++) { if (!strcmp(argv[i], "--target") && i + 1 < argc) target = atol(argv[++i]); else if (!strcmp(argv[i], "--quiet")) quiet = 1; }
    rng_state = seed * 0x9E3779B97F4A7C15ULL + 0x2545F491ULL; for (int i = 0; i < 20; i++) rnd();
    FILE *f = fopen(argv[1], "r"); if (!f) { perror("graph"); return 2; }
    if (fscanf(f, "%d %d", &n, &m) != 2) return 2;
    long *w = malloc(sizeof(long) * n);
    for (int i = 0; i < n; i++) { if (fscanf(f, "%ld", &w[i]) != 1) return 2; if (w[i] != w[0]) { fprintf(stderr, "weighted instance: NuMVC variant is unweighted only\n"); return 3; } }
    wcommon = w[0];
    eu = malloc(sizeof(int) * m); ev = malloc(sizeof(int) * m);
    int *deg = calloc(n, sizeof(int)); int mm = 0;
    for (int i = 0; i < m; i++) { int a, b; if (fscanf(f, "%d %d", &a, &b) != 2) return 2; if (a == b) continue; eu[mm] = a; ev[mm] = b; mm++; deg[a]++; deg[b]++; }
    fclose(f); m = mm;
    astart = malloc(sizeof(int) * (n + 1)); astart[0] = 0; for (int i = 0; i < n; i++) astart[i + 1] = astart[i] + deg[i];
    adj = malloc(sizeof(int) * 2 * (size_t)m); eadj = malloc(sizeof(int) * 2 * (size_t)m);
    int *fill = calloc(n, sizeof(int));
    for (int e = 0; e < m; e++) { int a = eu[e], b = ev[e]; adj[astart[a] + fill[a]] = b; eadj[astart[a] + fill[a]++] = e; adj[astart[b] + fill[b]] = a; eadj[astart[b] + fill[b]++] = e; }
    free(fill);
    ew = malloc(sizeof(long) * m); for (int e = 0; e < m; e++) ew[e] = 1;
    inC = calloc(n, 1); conf = malloc(n); memset(conf, 1, n);
    dscore = calloc(n, sizeof(long)); age_ = calloc(n, sizeof(long));
    uncov = malloc(sizeof(int) * m); uncov_pos = malloc(sizeof(int) * m); nuncov = 0;
    Clist = malloc(sizeof(int) * n); Cpos = malloc(sizeof(int) * n); nC = 0;
    for (int v = 0; v < n; v++) Cpos[v] = -1;
    /* initial: all edges uncovered, dscore[v] = deg(v); greedy cover by max dscore */
    for (int e = 0; e < m; e++) uncov_add(e);
    for (int v = 0; v < n; v++) dscore[v] = deg[v];
    long step = 0;
    while (nuncov > 0) {
        int e = uncov[rnd_int(nuncov)]; int a = eu[e], b = ev[e];
        int v = dscore[a] >= dscore[b] ? a : b;
        add_to_cover(v, step++);
    }
    /* drop redundant cover vertices (dscore == 0) */
    for (int i = nC - 1; i >= 0; i--) { int v = Clist[i]; if (inC[v] && dscore[v] == 0) remove_from_cover(v, step++); }
    double t0 = now();
    int *best = malloc(sizeof(int) * n); int nbest = 0;
    long bestsize = n - nC;   /* independent set size */
    for (int v = 0; v < n; v++) if (!inC[v]) best[nbest++] = v;
    if (!quiet) { printf("improved %.2f %ld %d\n", now() - t0, bestsize * wcommon, nbest); fflush(stdout); }
    /* NuMVC main loop: try to find a cover of size nC-1 */
    {
        /* remove the best-dscore vertex to start the search for a smaller cover */
        int bv = Clist[0]; for (int i = 1; i < nC; i++) { int v = Clist[i]; if (dscore[v] > dscore[bv] || (dscore[v] == dscore[bv] && age_[v] < age_[bv])) bv = v; }
        remove_from_cover(bv, step++);
    }
    long iters = 0; const long GAMMA = (long)m / 2 + 1; const double RHO = 0.3;
    while (1) {
        iters++;
        if ((iters & 255) == 0 && now() - t0 > secs) break;
        if (nuncov == 0) {
            /* found a cover of size nC: record, then remove one more */
            if ((long)(n - nC) > bestsize) {
                bestsize = n - nC; nbest = 0; for (int v = 0; v < n; v++) if (!inC[v]) best[nbest++] = v;
                if (!quiet) { printf("improved %.2f %ld %d\n", now() - t0, bestsize * wcommon, nbest); fflush(stdout); }
                if (target > 0 && bestsize * wcommon >= target) break;
            }
            int bv = Clist[0]; for (int i = 1; i < nC; i++) { int v = Clist[i]; if (dscore[v] > dscore[bv] || (dscore[v] == dscore[bv] && age_[v] < age_[bv])) bv = v; }
            remove_from_cover(bv, step++);
            continue;
        }
        /* choose u in C with highest dscore (oldest on ties): full scan of C */
        int u = Clist[0]; long bu = dscore[u]; long au = age_[u];
        for (int i = 1; i < nC; i++) {
            int v = Clist[i];
            if (dscore[v] > bu || (dscore[v] == bu && age_[v] < au)) { u = v; bu = dscore[v]; au = age_[v]; }
        }
        remove_from_cover(u, step++);
        /* pick a random uncovered edge, add the endpoint with higher dscore that is conf-allowed, ties by age */
        int e = uncov[rnd_int(nuncov)]; int a = eu[e], b = ev[e]; int v;
        if (conf[a] && !conf[b]) v = a; else if (conf[b] && !conf[a]) v = b;
        else if (dscore[a] > dscore[b] || (dscore[a] == dscore[b] && age_[a] < age_[b])) v = a; else v = b;
        add_to_cover(v, step++);
        /* edge weighting on remaining uncovered edges */
        for (int i = 0; i < nuncov; i++) { int ue = uncov[i]; ew[ue]++; dscore[eu[ue]]++; dscore[ev[ue]]++; }
        if (nuncov > 0 && (iters % 50) == 0) {
            long tot = 0; for (int i = 0; i < m; i++) tot += ew[i];
            if (tot >= GAMMA * 4) {
                for (int i = 0; i < m; i++) ew[i] = (long)(ew[i] * RHO); if (1) {
                    /* recompute dscores from scratch after forgetting */
                    for (int vv = 0; vv < n; vv++) dscore[vv] = 0;
                    for (int ee = 0; ee < m; ee++) {
                        int a2 = eu[ee], b2 = ev[ee];
                        if (!inC[a2] && !inC[b2]) { dscore[a2] += ew[ee]; dscore[b2] += ew[ee]; }
                        else if (inC[a2] && !inC[b2]) dscore[a2] -= ew[ee];
                        else if (!inC[a2] && inC[b2]) dscore[b2] -= ew[ee];
                    }
                }
            }
        }
    }
    printf("iters %ld time %.1f\n", iters, now() - t0);
    printf("BEST %ld %d :", bestsize * wcommon, nbest);
    for (int i = 0; i < nbest; i++) printf(" %d", best[i]);
    printf("\n");
    return 0;
}
