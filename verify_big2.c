/*
 * verify_big2.c -- independent, scalable verifier for quaternary edit-metric
 * codes (one of two verifiers for the n = 12 cells; CLAUDE.md rule 3).
 *
 * Definition checked (Houghten's table of bounds on quaternary edit-metric
 * codes): C is a set of distinct words of length n over {0,1,2,3}.  C has
 * minimum edit distance d iff every pair of distinct codewords has Levenshtein
 * distance >= d, where Levenshtein distance is the least number of
 * single-symbol insertions, deletions and substitutions (each of cost 1)
 * turning one word into the other.
 *
 * Method.  Every unordered pair {i, j} of codewords is checked, one pair at a
 * time, with the bit-parallel global edit distance of Myers (1999) in Hyyro's
 * (2003) formulation; one 64-bit machine word holds a whole DP column because
 * n <= 32.  There is no pruning by symmetry, balls or distance lower bounds:
 * the only early exit is to stop once a pair below d has been found, at which
 * point the verdict "invalid" is fixed.  The pair loop is split over at most
 * two POSIX threads (rows handed out one at a time).  The reported violating
 * pair is always the lexicographically least (i, j), i < j, in file order,
 * whatever the thread timing.  All arithmetic in the verdict is integer.
 *
 * Input format (as the artifacts in best/):
 *   - a line whose first byte is '#' is a comment;
 *   - '# n = <int>' is required (or the option --n);
 *   - '# d = <int>' is optional; if it is absent, --d is required.  d is
 *     never guessed.  A header and an option that disagree are refused;
 *   - '# q = <int>', if present, must be 4;
 *   - every other non-empty line is one codeword.  Trailing whitespace
 *     (including CR) is stripped; nothing else is.
 *
 * Output (stdout): n, d, size, valid: True/False, the first violation if
 * invalid, and sha256 of the byte-sorted codewords, each followed by '\n',
 * concatenated.  Timing statistics go to stderr.
 *
 * Exit: 0 iff valid, 1 if invalid, 2 on malformed input or refused
 * parameters (and on any failure that prevents a verdict).
 *
 * Build:  cc -O3 -o verify_big2 verify_big2.c -lpthread
 * Usage:  ./verify_big2 [--d D] [--n N] [--threads 1|2] FILE
 *         ./verify_big2 --selftest
 */
#define _POSIX_C_SOURCE 200809L
#include <errno.h>
#include <pthread.h>
#include <stdarg.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define MAX_N 32
#define MAX_D 1000000L
#define MAX_THREADS 2

/* ------------------------------------------------------------------------ */
/* refusals                                                                 */
/* ------------------------------------------------------------------------ */

static void refuse(const char *fmt, ...)
{
    va_list ap;
    fprintf(stderr, "verify_big2: refused: ");
    va_start(ap, fmt);
    vfprintf(stderr, fmt, ap);
    va_end(ap);
    fputc('\n', stderr);
    exit(2);
}

static void *xmalloc(size_t nbytes)
{
    void *p = malloc(nbytes ? nbytes : 1);
    if (!p)
        refuse("out of memory (%zu bytes)", nbytes);
    return p;
}

/* ------------------------------------------------------------------------ */
/* SHA-256 (FIPS 180-4), self-contained                                     */
/* ------------------------------------------------------------------------ */

typedef struct {
    uint32_t h[8];
    uint64_t nbytes;
    unsigned char buf[64];
    size_t nbuf;
} sha256_ctx;

static const uint32_t SHA_K[64] = {
    0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5, 0x3956c25b, 0x59f111f1,
    0x923f82a4, 0xab1c5ed5, 0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3,
    0x72be5d74, 0x80deb1fe, 0x9bdc06a7, 0xc19bf174, 0xe49b69c1, 0xefbe4786,
    0x0fc19dc6, 0x240ca1cc, 0x2de92c6f, 0x4a7484aa, 0x5cb0a9dc, 0x76f988da,
    0x983e5152, 0xa831c66d, 0xb00327c8, 0xbf597fc7, 0xc6e00bf3, 0xd5a79147,
    0x06ca6351, 0x14292967, 0x27b70a85, 0x2e1b2138, 0x4d2c6dfc, 0x53380d13,
    0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85, 0xa2bfe8a1, 0xa81a664b,
    0xc24b8b70, 0xc76c51a3, 0xd192e819, 0xd6990624, 0xf40e3585, 0x106aa070,
    0x19a4c116, 0x1e376c08, 0x2748774c, 0x34b0bcb5, 0x391c0cb3, 0x4ed8aa4a,
    0x5b9cca4f, 0x682e6ff3, 0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208,
    0x90befffa, 0xa4506ceb, 0xbef9a3f7, 0xc67178f2};

#define ROR32(x, k) (((x) >> (k)) | ((x) << (32 - (k))))

static void sha256_block(uint32_t h[8], const unsigned char *p)
{
    uint32_t w[64], a, b, c, d, e, f, g, hh;
    int t;
    for (t = 0; t < 16; t++)
        w[t] = ((uint32_t)p[4 * t] << 24) | ((uint32_t)p[4 * t + 1] << 16) |
               ((uint32_t)p[4 * t + 2] << 8) | (uint32_t)p[4 * t + 3];
    for (t = 16; t < 64; t++) {
        uint32_t s0 = ROR32(w[t - 15], 7) ^ ROR32(w[t - 15], 18) ^ (w[t - 15] >> 3);
        uint32_t s1 = ROR32(w[t - 2], 17) ^ ROR32(w[t - 2], 19) ^ (w[t - 2] >> 10);
        w[t] = w[t - 16] + s0 + w[t - 7] + s1;
    }
    a = h[0]; b = h[1]; c = h[2]; d = h[3];
    e = h[4]; f = h[5]; g = h[6]; hh = h[7];
    for (t = 0; t < 64; t++) {
        uint32_t S1 = ROR32(e, 6) ^ ROR32(e, 11) ^ ROR32(e, 25);
        uint32_t ch = (e & f) ^ (~e & g);
        uint32_t t1 = hh + S1 + ch + SHA_K[t] + w[t];
        uint32_t S0 = ROR32(a, 2) ^ ROR32(a, 13) ^ ROR32(a, 22);
        uint32_t mj = (a & b) ^ (a & c) ^ (b & c);
        uint32_t t2 = S0 + mj;
        hh = g; g = f; f = e; e = d + t1;
        d = c; c = b; b = a; a = t1 + t2;
    }
    h[0] += a; h[1] += b; h[2] += c; h[3] += d;
    h[4] += e; h[5] += f; h[6] += g; h[7] += hh;
}

static void sha256_init(sha256_ctx *s)
{
    static const uint32_t H0[8] = {0x6a09e667, 0xbb67ae85, 0x3c6ef372,
                                   0xa54ff53a, 0x510e527f, 0x9b05688c,
                                   0x1f83d9ab, 0x5be0cd19};
    memcpy(s->h, H0, sizeof H0);
    s->nbytes = 0;
    s->nbuf = 0;
}

static void sha256_update(sha256_ctx *s, const void *data, size_t len)
{
    const unsigned char *p = data;
    s->nbytes += len;
    while (len > 0) {
        size_t take = 64 - s->nbuf;
        if (take > len)
            take = len;
        memcpy(s->buf + s->nbuf, p, take);
        s->nbuf += take;
        p += take;
        len -= take;
        if (s->nbuf == 64) {
            sha256_block(s->h, s->buf);
            s->nbuf = 0;
        }
    }
}

static void sha256_final(sha256_ctx *s, char hex[65])
{
    uint64_t bits = s->nbytes * 8;
    unsigned char pad = 0x80, zero = 0, lenbe[8];
    int i;
    sha256_update(s, &pad, 1);
    while (s->nbuf != 56)
        sha256_update(s, &zero, 1);
    for (i = 0; i < 8; i++)
        lenbe[i] = (unsigned char)(bits >> (56 - 8 * i));
    sha256_update(s, lenbe, 8);
    for (i = 0; i < 8; i++)
        sprintf(hex + 8 * i, "%08x", s->h[i]);
    hex[64] = '\0';
}

/* ------------------------------------------------------------------------ */
/* edit distance                                                            */
/* ------------------------------------------------------------------------ */

/* A word of length n over {0,1,2,3} is packed as sum_k sym[k] << (2k). */

static void build_peq(uint64_t word, int n, uint64_t peq[4])
{
    int k;
    peq[0] = peq[1] = peq[2] = peq[3] = 0;
    for (k = 0; k < n; k++)
        peq[(word >> (2 * k)) & 3] |= (uint64_t)1 << k;
}

/*
 * Global Levenshtein distance between the pattern whose match masks are peq
 * (length n) and the packed text word (length n), Myers/Hyyro bit-vectors.
 * Bit k of Pv/Mv is +1/-1 vertical delta D[k+1][j] - D[k][j]; score tracks
 * D[n][j].  The top row is D[0][j] = j, hence the '| 1' on the horizontal
 * positive delta.  Bits above n-1 hold garbage that never reaches the low n
 * bits (only additions and left shifts move information, and only upwards).
 */
#if defined(__GNUC__) || defined(__clang__)
#define ALWAYS_INLINE inline __attribute__((always_inline))
#else
#define ALWAYS_INLINE inline
#endif

static ALWAYS_INLINE int lev_myers(const uint64_t peq[4], uint64_t text, int n)
{
    uint64_t pv = ~(uint64_t)0, mv = 0;
    int score = n, k;
    for (k = 0; k < n; k++) {
        uint64_t eq = peq[(text >> (2 * k)) & 3];
        uint64_t xv = eq | mv;
        uint64_t xh = (((eq & pv) + pv) ^ pv) | eq;
        uint64_t ph = mv | ~(xh | pv);
        uint64_t mh = pv & xh;
        score += (int)((ph >> (n - 1)) & 1) - (int)((mh >> (n - 1)) & 1);
        ph = (ph << 1) | 1;
        mh = mh << 1;
        pv = mh | ~(xv | ph);
        mv = ph & xv;
    }
    return score;
}

/*
 * The same recurrence as lev_myers, run on LANES independent texts against
 * one pattern in lockstep.  Nothing is shared between lanes; interleaving
 * only gives the CPU independent instruction streams (one pair alone is a
 * serial dependency chain).
 */
#ifndef LANES
#define LANES 16
#endif
static ALWAYS_INLINE void lev_myers_block(const uint64_t peq[4],
                                          const uint64_t *text, int n,
                                          int out[LANES])
{
    uint64_t pv[LANES], mv[LANES];
    int score[LANES], k, l;
    for (l = 0; l < LANES; l++) {
        pv[l] = ~(uint64_t)0;
        mv[l] = 0;
        score[l] = n;
    }
    for (k = 0; k < n; k++) {
        for (l = 0; l < LANES; l++) {
            uint64_t eq = peq[(text[l] >> (2 * k)) & 3];
            uint64_t xv = eq | mv[l];
            uint64_t xh = (((eq & pv[l]) + pv[l]) ^ pv[l]) | eq;
            uint64_t ph = mv[l] | ~(xh | pv[l]);
            uint64_t mh = pv[l] & xh;
            score[l] += (int)((ph >> (n - 1)) & 1) - (int)((mh >> (n - 1)) & 1);
            ph = (ph << 1) | 1;
            mh = mh << 1;
            pv[l] = mh | ~(xv | ph);
            mv[l] = ph & xv;
        }
    }
    for (l = 0; l < LANES; l++)
        out[l] = score[l];
}

/* Plain Wagner-Fischer, used only by --selftest to check lev_myers. */
static int lev_wf(const unsigned char *a, int la, const unsigned char *b, int lb)
{
    int prev[MAX_N + 1], cur[MAX_N + 1], i, j;
    for (j = 0; j <= lb; j++)
        prev[j] = j;
    for (i = 1; i <= la; i++) {
        cur[0] = i;
        for (j = 1; j <= lb; j++) {
            int best = prev[j] + 1;
            if (cur[j - 1] + 1 < best)
                best = cur[j - 1] + 1;
            if (prev[j - 1] + (a[i - 1] != b[j - 1]) < best)
                best = prev[j - 1] + (a[i - 1] != b[j - 1]);
            cur[j] = best;
        }
        memcpy(prev, cur, sizeof(int) * (size_t)(lb + 1));
    }
    return prev[lb];
}

/*
 * Scan row i: words j = i+1 .. N-1.  Returns the least j with distance < d,
 * or N if there is none.  *pairs counts the pairs examined.
 */
static ALWAYS_INLINE size_t scan_row_impl(const uint64_t *w, size_t i, size_t N,
                                          int n, int d, int *dist_out,
                                          uint64_t *pairs)
{
    uint64_t peq[4];
    size_t j;
    build_peq(w[i], n, peq);
    for (j = i + 1; j + LANES <= N; j += LANES) {
        int dist[LANES], any = 0, l;
        lev_myers_block(peq, w + j, n, dist);
        for (l = 0; l < LANES; l++)
            any |= dist[l] < d;
        if (any) {
            for (l = 0; dist[l] >= d; l++)
                ;
            *dist_out = dist[l];
            *pairs += (uint64_t)(j + (size_t)l - i);
            return j + (size_t)l;
        }
    }
    for (; j < N; j++) {
        int dist = lev_myers(peq, w[j], n);
        if (dist < d) {
            *dist_out = dist;
            *pairs += (uint64_t)(j - i);
            return j;
        }
    }
    *pairs += (uint64_t)(N - 1 - i);
    return N;
}

/* Same computation, with the word length a compile-time constant where it
 * matters for speed; every case runs the identical scan_row_impl. */
static size_t scan_row(const uint64_t *w, size_t i, size_t N, int n, int d,
                       int *dist_out, uint64_t *pairs)
{
    switch (n) {
    case 8:  return scan_row_impl(w, i, N, 8, d, dist_out, pairs);
    case 9:  return scan_row_impl(w, i, N, 9, d, dist_out, pairs);
    case 10: return scan_row_impl(w, i, N, 10, d, dist_out, pairs);
    case 11: return scan_row_impl(w, i, N, 11, d, dist_out, pairs);
    case 12: return scan_row_impl(w, i, N, 12, d, dist_out, pairs);
    case 13: return scan_row_impl(w, i, N, 13, d, dist_out, pairs);
    case 14: return scan_row_impl(w, i, N, 14, d, dist_out, pairs);
    case 15: return scan_row_impl(w, i, N, 15, d, dist_out, pairs);
    case 16: return scan_row_impl(w, i, N, 16, d, dist_out, pairs);
    default: return scan_row_impl(w, i, N, n, d, dist_out, pairs);
    }
}

/* ------------------------------------------------------------------------ */
/* threaded pair loop                                                       */
/* ------------------------------------------------------------------------ */

typedef struct {
    const uint64_t *w;
    size_t N;
    int n, d;
    pthread_mutex_t mu;
    size_t next_row;        /* next row to hand out                   */
    int found;              /* a violating pair has been recorded     */
    size_t best_i, best_j;  /* least violating (i, j) recorded so far */
    int best_dist;
    uint64_t pairs;         /* statistics only                        */
} pairjob;

static void *pair_worker(void *arg)
{
    pairjob *J = arg;
    uint64_t local_pairs = 0;
    for (;;) {
        size_t i, j;
        int stop, dist = -1;
        pthread_mutex_lock(&J->mu);
        i = J->next_row++;
        /* rows beyond the best violating row cannot give a smaller pair */
        stop = (i + 1 >= J->N) || (J->found && i > J->best_i);
        pthread_mutex_unlock(&J->mu);
        if (stop)
            break;
        j = scan_row(J->w, i, J->N, J->n, J->d, &dist, &local_pairs);
        if (j < J->N) {
            pthread_mutex_lock(&J->mu);
            if (!J->found || i < J->best_i ||
                (i == J->best_i && j < J->best_j)) {
                J->found = 1;
                J->best_i = i;
                J->best_j = j;
                J->best_dist = dist;
            }
            pthread_mutex_unlock(&J->mu);
        }
    }
    pthread_mutex_lock(&J->mu);
    J->pairs += local_pairs;
    pthread_mutex_unlock(&J->mu);
    return NULL;
}

/* Returns 1 and fills (*bi, *bj, *bdist) with the least violating pair in
 * index order if some pair has distance < d, else 0. */
static int check_pairs(const uint64_t *w, size_t N, int n, int d, int nthreads,
                       size_t *bi, size_t *bj, int *bdist, uint64_t *pairs)
{
    pairjob J;
    pthread_t th;
    int spawned = 0;
    J.w = w;
    J.N = N;
    J.n = n;
    J.d = d;
    J.next_row = 0;
    J.found = 0;
    J.best_i = J.best_j = 0;
    J.best_dist = -1;
    J.pairs = 0;
    if (pthread_mutex_init(&J.mu, NULL) != 0)
        refuse("pthread_mutex_init failed");
    if (nthreads >= 2 && N > 2) {
        if (pthread_create(&th, NULL, pair_worker, &J) != 0)
            refuse("pthread_create failed");
        spawned = 1;
    }
    pair_worker(&J); /* the main thread is the other worker */
    if (spawned && pthread_join(th, NULL) != 0)
        refuse("pthread_join failed");
    pthread_mutex_destroy(&J.mu);
    *pairs = J.pairs;
    if (J.found) {
        *bi = J.best_i;
        *bj = J.best_j;
        *bdist = J.best_dist;
    }
    return J.found;
}

/* ------------------------------------------------------------------------ */
/* input                                                                    */
/* ------------------------------------------------------------------------ */

typedef struct {
    const char *s;  /* not NUL-terminated in general; may contain any byte */
    size_t len;
    long lineno;    /* 1-based physical line number */
} wordline;

static char *read_all(const char *path, size_t *len_out)
{
    FILE *f = fopen(path, "rb");
    size_t cap = (size_t)1 << 16, len = 0;
    char *buf;
    if (!f)
        refuse("cannot open %s: %s", path, strerror(errno));
    buf = xmalloc(cap + 1);
    for (;;) {
        size_t r;
        if (len == cap) {
            char *nb;
            cap *= 2;
            nb = realloc(buf, cap + 1);
            if (!nb)
                refuse("out of memory reading %s", path);
            buf = nb;
        }
        r = fread(buf + len, 1, cap - len, f);
        len += r;
        if (r == 0) {
            if (ferror(f))
                refuse("read error on %s", path);
            break;
        }
    }
    fclose(f);
    buf[len] = '\0';
    *len_out = len;
    return buf;
}

static int is_space(char c)
{
    return c == ' ' || c == '\t' || c == '\r' || c == '\v' || c == '\f';
}

/*
 * If the comment line (NUL-terminated, starting with '#') has the form
 * '#' [ws] key [ws] '=' [ws] <digits> (end | ws | ','), store the value and
 * return 1.  Return 0 if the line is not a header for this key.  A line that
 * starts like a header for this key but has no well-formed non-negative
 * integer is refused.
 */
static int header_value(const char *line, char key, long lineno, long *out)
{
    const char *p = line + 1;
    long v = 0;
    while (*p == ' ' || *p == '\t')
        p++;
    if (*p != key)
        return 0;
    p++;
    while (*p == ' ' || *p == '\t')
        p++;
    if (*p != '=')
        return 0;
    p++;
    while (*p == ' ' || *p == '\t')
        p++;
    if (*p < '0' || *p > '9')
        refuse("line %ld: '# %c =' is not followed by a non-negative integer",
               lineno, key);
    while (*p >= '0' && *p <= '9') {
        v = v * 10 + (*p - '0');
        if (v > MAX_D)
            refuse("line %ld: value of '# %c' is too large", lineno, key);
        p++;
    }
    if (*p != '\0' && *p != ' ' && *p != '\t' && *p != ',')
        refuse("line %ld: '# %c = <int>' followed by unexpected text", lineno,
               key);
    *out = v;
    return 1;
}

static void set_header(long *slot, long v, char key, long lineno)
{
    if (*slot >= 0 && *slot != v)
        refuse("line %ld: '# %c = %ld' conflicts with an earlier '# %c = %ld'",
               lineno, key, v, key, *slot);
    *slot = v;
}

static long parse_option_int(const char *opt, const char *val)
{
    char *end;
    long v;
    if (!val)
        refuse("option %s needs a value", opt);
    errno = 0;
    v = strtol(val, &end, 10);
    if (errno || end == val || *end != '\0' || v < 0 || v > MAX_D)
        refuse("option %s: '%s' is not a non-negative integer", opt, val);
    return v;
}

/* Print a word for humans: printable ASCII as-is, other bytes as '?'. */
static void print_word(const char *s, size_t len)
{
    size_t k, lim = len > 80 ? 80 : len;
    for (k = 0; k < lim; k++) {
        unsigned char c = (unsigned char)s[k];
        putchar(c >= 33 && c < 127 ? c : '?');
    }
    if (len > lim)
        printf("...");
}

static int cmp_wordline(const void *A, const void *B)
{
    const wordline *a = A, *b = B;
    size_t m = a->len < b->len ? a->len : b->len;
    int c = memcmp(a->s, b->s, m);
    if (c)
        return c;
    if (a->len != b->len)
        return a->len < b->len ? -1 : 1;
    return (a->lineno > b->lineno) - (a->lineno < b->lineno);
}

static double now_seconds(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + 1e-9 * (double)ts.tv_nsec;
}

/* ------------------------------------------------------------------------ */
/* self-test                                                                */
/* ------------------------------------------------------------------------ */

static uint64_t rng_state;
static uint64_t splitmix64(void)
{
    uint64_t z = (rng_state += 0x9e3779b97f4a7c15ULL);
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
    return z ^ (z >> 31);
}

static uint64_t pack(const unsigned char *sym, int n)
{
    uint64_t w = 0;
    int k;
    for (k = 0; k < n; k++)
        w |= (uint64_t)sym[k] << (2 * k);
    return w;
}

/* b := a with `edits` random length-preserving edits (substitution, or a
 * deletion paired with an insertion). */
static void mutate(const unsigned char *a, unsigned char *b, int n, int edits)
{
    int e, k;
    memcpy(b, a, (size_t)n);
    for (e = 0; e < edits; e++) {
        if (splitmix64() & 1) {
            b[splitmix64() % (uint64_t)n] = (unsigned char)(splitmix64() & 3);
        } else {
            int del = (int)(splitmix64() % (uint64_t)n);
            int ins = (int)(splitmix64() % (uint64_t)n);
            unsigned char tmp[MAX_N];
            int t = 0;
            for (k = 0; k < n; k++)
                if (k != del)
                    tmp[t++] = b[k];
            /* tmp has n-1 symbols; insert one at position ins */
            for (k = 0, t = 0; k < n; k++) {
                if (k == ins)
                    b[k] = (unsigned char)(splitmix64() & 3);
                else
                    b[k] = tmp[t++];
            }
        }
    }
}

static int selftest(void)
{
    static const struct { const char *msg; const char *hex; } vec[] = {
        {"", "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"},
        {"abc", "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad"},
        {"abcdbcdecdefdefgefghfghighijhijkijkljklmklmnlmnomnopnopq",
         "248d6a61d20638b8e5c026930c3e6039a33ce45964ff2167f6ecedd419db06c1"},
    };
    sha256_ctx s;
    char hex[65];
    unsigned char a[MAX_N], b[MAX_N];
    uint64_t peq[4];
    int failures = 0, L, v;
    long checked = 0;

    for (v = 0; v < 3; v++) {
        sha256_init(&s);
        sha256_update(&s, vec[v].msg, strlen(vec[v].msg));
        sha256_final(&s, hex);
        if (strcmp(hex, vec[v].hex) != 0) {
            printf("selftest FAIL: sha256 vector %d\n", v);
            failures++;
        }
    }
    {   /* one million 'a', fed in uneven pieces */
        static char chunk[1000];
        long done = 0;
        memset(chunk, 'a', sizeof chunk);
        sha256_init(&s);
        while (done < 1000000) {
            long take = 1 + (done % 997);
            if (take > 1000000 - done)
                take = 1000000 - done;
            sha256_update(&s, chunk, (size_t)take);
            done += take;
        }
        sha256_final(&s, hex);
        if (strcmp(hex, "cdc76e5c9914fb9281a1c7e284d73e67f1809a48a497200e046d39ccc7112cd0")) {
            printf("selftest FAIL: sha256 of 10^6 'a'\n");
            failures++;
        }
    }

    /* Myers (single and LANES-wide block) vs Wagner-Fischer: all pairs of
     * words of length 1..6 */
    for (L = 1; L <= 6; L++) {
        uint64_t total = (uint64_t)1 << (2 * L), x, y;
        for (x = 0; x < total; x++) {
            int k, blk[LANES], have_blk = 0;
            for (k = 0; k < L; k++)
                a[k] = (unsigned char)((x >> (2 * k)) & 3);
            build_peq(x, L, peq);
            for (y = 0; y < total; y++) {
                int ref;
                if (y % LANES == 0) {
                    have_blk = (y + LANES <= total);
                    if (have_blk) {
                        uint64_t tx[LANES];
                        int l;
                        for (l = 0; l < LANES; l++)
                            tx[l] = y + (uint64_t)l;
                        lev_myers_block(peq, tx, L, blk);
                    }
                }
                for (k = 0; k < L; k++)
                    b[k] = (unsigned char)((y >> (2 * k)) & 3);
                ref = lev_wf(a, L, b, L);
                if (lev_myers(peq, y, L) != ref ||
                    (have_blk && blk[y % LANES] != ref)) {
                    if (failures++ < 10)
                        printf("selftest FAIL: L=%d x=%llu y=%llu\n", L,
                               (unsigned long long)x, (unsigned long long)y);
                }
                checked++;
            }
        }
    }
    /* random pairs and near pairs for lengths 7..32 (fixed seed) */
    rng_state = 20261009;
    for (L = 7; L <= MAX_N; L++) {
        int t, k, l;
        for (t = 0; t < 10000; t++) {
            unsigned char bb[LANES][MAX_N];
            uint64_t tx[LANES];
            int blk[LANES];
            for (k = 0; k < L; k++)
                a[k] = (unsigned char)(splitmix64() & 3);
            for (l = 0; l < LANES; l++) {
                if (l & 1)
                    mutate(a, bb[l], L, 1 + (int)(splitmix64() % 4));
                else
                    for (k = 0; k < L; k++)
                        bb[l][k] = (unsigned char)(splitmix64() & 3);
                tx[l] = pack(bb[l], L);
            }
            build_peq(pack(a, L), L, peq);
            lev_myers_block(peq, tx, L, blk);
            for (l = 0; l < LANES; l++) {
                int ref = lev_wf(a, L, bb[l], L);
                if (lev_myers(peq, tx[l], L) != ref || blk[l] != ref) {
                    if (failures++ < 10)
                        printf("selftest FAIL: random pair at L=%d\n", L);
                }
                checked++;
            }
        }
    }

    /* the threaded scan vs a brute-force first violating pair */
    {
        static const int lens[] = {3, 5, 8, 12, 16};
        int trial;
        for (trial = 0; trial < 600; trial++) {
            int n = lens[trial % 5];
            int N = 2 + (int)(splitmix64() % 60);
            int d = 1 + (int)(splitmix64() % (uint64_t)(n + 1));
            unsigned char sym[64][MAX_N];
            uint64_t w[64], pairs;
            int i, j, k, ri = -1, rj = -1, rdist = -1, th;
            for (i = 0; i < N; i++) {
                int fresh;
                do {
                    if (i > 0 && (splitmix64() % 3))
                        mutate(sym[splitmix64() % (uint64_t)i], sym[i], n,
                               1 + (int)(splitmix64() % 3));
                    else
                        for (k = 0; k < n; k++)
                            sym[i][k] = (unsigned char)(splitmix64() & 3);
                    w[i] = pack(sym[i], n);
                    fresh = 1;
                    for (j = 0; j < i; j++)
                        if (w[j] == w[i])
                            fresh = 0;
                } while (!fresh);
            }
            for (i = 0; i < N && ri < 0; i++)
                for (j = i + 1; j < N; j++) {
                    int dd = lev_wf(sym[i], n, sym[j], n);
                    if (dd < d) {
                        ri = i; rj = j; rdist = dd;
                        break;
                    }
                }
            for (th = 1; th <= 2; th++) {
                size_t bi = 0, bj = 0;
                int bdist = -1;
                int found = check_pairs(w, (size_t)N, n, d, th, &bi, &bj, &bdist, &pairs);
                int ok = (ri < 0) ? !found
                                  : (found && (int)bi == ri && (int)bj == rj && bdist == rdist);
                if (!ok && failures++ < 10)
                    printf("selftest FAIL: scan trial %d (n=%d N=%d d=%d threads=%d)\n",
                           trial, n, N, d, th);
            }
        }
    }
    if (failures) {
        printf("selftest: FAIL (%d failures)\n", failures);
        return 1;
    }
    printf("selftest: PASS (sha256 vectors, %ld Myers/Wagner-Fischer pairs, 600 scan trials x 2 thread counts)\n",
           checked);
    return 0;
}

/* ------------------------------------------------------------------------ */
/* main                                                                     */
/* ------------------------------------------------------------------------ */

static void usage_refuse(void)
{
    refuse("usage: verify_big2 [--d D] [--n N] [--threads 1|2] FILE | --selftest");
}

int main(int argc, char **argv)
{
    const char *path = NULL;
    long n_opt = -1, d_opt = -1, threads = MAX_THREADS;
    long n_hdr = -1, d_hdr = -1, q_hdr = -1;
    long n, d, lineno = 0;
    char *buf;
    size_t len, pos = 0, nwords = 0, capwords = 1024, t;
    wordline *words, *sorted;
    uint64_t *packed, pairs = 0;
    sha256_ctx sha;
    char hex[65];
    int a, valid = 1;
    double t0, t1;

    for (a = 1; a < argc; a++) {
        if (strcmp(argv[a], "--selftest") == 0) {
            if (argc != 2)
                usage_refuse();
            return selftest();
        } else if (strcmp(argv[a], "--d") == 0) {
            if (d_opt >= 0)
                refuse("--d given twice");
            d_opt = parse_option_int("--d", a + 1 < argc ? argv[a + 1] : NULL);
            a++;
        } else if (strcmp(argv[a], "--n") == 0) {
            if (n_opt >= 0)
                refuse("--n given twice");
            n_opt = parse_option_int("--n", a + 1 < argc ? argv[a + 1] : NULL);
            a++;
        } else if (strcmp(argv[a], "--threads") == 0) {
            threads = parse_option_int("--threads", a + 1 < argc ? argv[a + 1] : NULL);
            if (threads < 1 || threads > MAX_THREADS)
                refuse("--threads must be 1 or 2");
            a++;
        } else if (argv[a][0] == '-' && argv[a][1] != '\0') {
            refuse("unknown option %s", argv[a]);
        } else {
            if (path)
                usage_refuse();
            path = argv[a];
        }
    }
    if (!path)
        usage_refuse();

    t0 = now_seconds();
    buf = read_all(path, &len);
    words = xmalloc(capwords * sizeof *words);

    /* split into lines; headers are read wherever they appear */
    while (pos < len) {
        size_t end = pos, l;
        while (end < len && buf[end] != '\n')
            end++;
        lineno++;
        l = end - pos;
        while (l > 0 && is_space(buf[pos + l - 1]))
            l--;
        buf[pos + l] = '\0'; /* overwrites '\n', stripped space, or final NUL */
        if (l > 0) {
            if (buf[pos] == '#') {
                long v;
                if (header_value(buf + pos, 'n', lineno, &v))
                    set_header(&n_hdr, v, 'n', lineno);
                if (header_value(buf + pos, 'd', lineno, &v))
                    set_header(&d_hdr, v, 'd', lineno);
                if (header_value(buf + pos, 'q', lineno, &v))
                    set_header(&q_hdr, v, 'q', lineno);
            } else {
                if (nwords == capwords) {
                    wordline *nw;
                    capwords *= 2;
                    nw = realloc(words, capwords * sizeof *words);
                    if (!nw)
                        refuse("out of memory");
                    words = nw;
                }
                words[nwords].s = buf + pos;
                words[nwords].len = l;
                words[nwords].lineno = lineno;
                nwords++;
            }
        }
        pos = end + 1;
    }

    /* resolve parameters; never guess */
    if (n_hdr < 0 && n_opt < 0)
        refuse("no '# n = <int>' header and no --n option");
    if (n_hdr >= 0 && n_opt >= 0 && n_hdr != n_opt)
        refuse("'# n = %ld' in the file conflicts with --n %ld", n_hdr, n_opt);
    n = n_hdr >= 0 ? n_hdr : n_opt;
    if (n < 1 || n > MAX_N)
        refuse("n = %ld is outside the supported range 1..%d", n, MAX_N);
    if (d_hdr < 0 && d_opt < 0)
        refuse("no '# d = <int>' header and no --d option; d is never guessed");
    if (d_hdr >= 0 && d_opt >= 0 && d_hdr != d_opt)
        refuse("'# d = %ld' in the file conflicts with --d %ld", d_hdr, d_opt);
    d = d_hdr >= 0 ? d_hdr : d_opt;
    if (d < 1)
        refuse("d = %ld; minimum distance must be at least 1", d);
    if (q_hdr >= 0 && q_hdr != 4)
        refuse("'# q = %ld': this verifier is for the quaternary alphabet 0-3 only", q_hdr);

    printf("file: %s\n", path);
    printf("n: %ld\n", n);
    printf("d: %ld\n", d);
    printf("size: %zu\n", nwords);

    /* checksum and duplicate detection share one byte-wise sort */
    sorted = xmalloc(nwords * sizeof *sorted);
    if (nwords)
        memcpy(sorted, words, nwords * sizeof *sorted);
    qsort(sorted, nwords, sizeof *sorted, cmp_wordline);
    sha256_init(&sha);
    for (t = 0; t < nwords; t++) {
        sha256_update(&sha, sorted[t].s, sorted[t].len);
        sha256_update(&sha, "\n", 1);
    }
    sha256_final(&sha, hex);

    /* 1. length and alphabet, first offending line in file order */
    for (t = 0; t < nwords && valid; t++) {
        size_t k;
        if (words[t].len != (size_t)n) {
            valid = 0;
            printf("valid: False\n");
            printf("violation: length\n");
            printf("first_bad_line: %ld\n", words[t].lineno);
            printf("first_bad_word: ");
            print_word(words[t].s, words[t].len);
            printf("\nreason: word has length %zu, expected n = %ld\n",
                   words[t].len, n);
            break;
        }
        for (k = 0; k < words[t].len; k++) {
            unsigned char c = (unsigned char)words[t].s[k];
            if (c < '0' || c > '3') {
                valid = 0;
                printf("valid: False\n");
                printf("violation: symbol\n");
                printf("first_bad_line: %ld\n", words[t].lineno);
                printf("first_bad_word: ");
                print_word(words[t].s, words[t].len);
                printf("\nreason: byte 0x%02x at position %zu is not a symbol 0-3\n",
                       c, k);
                break;
            }
        }
    }

    /* 2. duplicates: equal words are adjacent in the sort */
    for (t = 1; t < nwords && valid; t++) {
        if (sorted[t - 1].len == sorted[t].len &&
            memcmp(sorted[t - 1].s, sorted[t].s, sorted[t].len) == 0) {
            valid = 0;
            printf("valid: False\n");
            printf("violation: duplicate\n");
            printf("first_violating_pair: ");
            print_word(sorted[t - 1].s, sorted[t - 1].len);
            putchar(' ');
            print_word(sorted[t].s, sorted[t].len);
            printf("\nfirst_violating_lines: %ld %ld\n", sorted[t - 1].lineno,
                   sorted[t].lineno);
            printf("first_violating_distance: 0\n");
            printf("reason: duplicate codeword\n");
        }
    }

    /* 3. every pair at Levenshtein distance >= d */
    if (valid) {
        size_t bi = 0, bj = 0;
        int bdist = -1;
        packed = xmalloc(nwords * sizeof *packed);
        for (t = 0; t < nwords; t++) {
            size_t k;
            uint64_t w = 0;
            for (k = 0; k < (size_t)n; k++)
                w |= (uint64_t)(words[t].s[k] - '0') << (2 * k);
            packed[t] = w;
        }
        if (nwords >= 2 &&
            check_pairs(packed, nwords, (int)n, (int)d, (int)threads, &bi, &bj,
                        &bdist, &pairs)) {
            valid = 0;
            printf("valid: False\n");
            printf("violation: distance\n");
            printf("first_violating_pair: ");
            print_word(words[bi].s, words[bi].len);
            putchar(' ');
            print_word(words[bj].s, words[bj].len);
            printf("\nfirst_violating_lines: %ld %ld\n", words[bi].lineno,
                   words[bj].lineno);
            printf("first_violating_distance: %d\n", bdist);
            printf("reason: Levenshtein distance %d < d = %ld\n", bdist, d);
        }
        free(packed);
    }
    if (valid)
        printf("valid: True\n");
    printf("sha256: %s\n", hex);

    t1 = now_seconds();
    fprintf(stderr,
            "verify_big2: threads=%ld pairs_checked=%llu seconds=%.3f pairs_per_second=%.4g\n",
            threads, (unsigned long long)pairs, t1 - t0,
            (t1 > t0 && pairs) ? (double)pairs / (t1 - t0) : 0.0);
    free(sorted);
    free(words);
    free(buf);
    return valid ? 0 : 1;
}
