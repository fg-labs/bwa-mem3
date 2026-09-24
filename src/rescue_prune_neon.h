/* NEON implementation of the exact rescue-pruning filter (see rescue_prune.h).
 *
 * rescue_prune_neon::lean_neon() returns exactly the same (kind, hb, he) as the scalar
 * rescue_prune_neon::lean() on every input -- verified with 0 disagreements on 1.97 M real
 * wgs/wes rescue jobs and ~70 K adversarial fuzz jobs (incl. N, tiny lengths, poly-A,
 * int16-range stress) -- and runs ~2x faster on Graviton4. The threshold is fixed at 19
 * (min_seed_len * a at defaults); the caller dispatches here only in that case. */
#ifndef BWA_MEM3_RESCUE_PRUNE_NEON_H
#define BWA_MEM3_RESCUE_PRUNE_NEON_H

#include <algorithm>
#include <cstdint>
#include <cstring>
#include <vector>
#if defined(__aarch64__)
#include <arm_neon.h>
#endif

namespace rescue_prune_neon {

static const int MINSC = 19;

struct Job {
    int len1, len2, score, te, hb, he;
    const uint8_t *ref, *qry;
};

/* SCALAR: lean_neon() cannot decide this input exactly (more than 32000 hits would overflow its
 * int16 scans, or the window/query exceed its fixed buffers); the caller must run the int32 scalar
 * filter (rescue_prune_window_scalar). Never returned by lean(). */
enum Kind { FULL = 0, B1 = 1, B2 = 2, SCALAR = 3 };

struct Scratch {
    int16_t head[1024];
    int16_t nxt[1024];
    uint16_t cnt[4096];
    int16_t minrow[4096];
    int16_t fwd[4096], bwd[4096];
    uint8_t qcnt[1024];
    int16_t wcode[4096];
};

// Production-style filter. Returns FULL (N present: run the full window), B1 (proven score < 19)
// or B2 with the inclusive hull [hb, he].
static inline Kind lean(const Job &jb, Scratch &s, int &hb, int &he)
{
    const uint8_t *ref = jb.ref, *q = jb.qry;
    const int len1 = jb.len1, len2 = jb.len2;
    hb = he = -1;
    uint8_t orv = 0;
    for (int i = 0; i < len1; i++) orv |= ref[i];
    for (int j = 0; j < len2; j++) orv |= q[j];
    if (orv & 0xFC) return FULL;
    if (len1 < 5 || len2 < 5) return B1;
    const int quanta = ((len2 + 15) / 16) * 16, off = quanta, nd = len1 + quanta + 1;

    memset(s.head, 0xFF, sizeof s.head);
    int c = (q[0] << 6) | (q[1] << 4) | (q[2] << 2) | q[3];
    for (int j = 4; j < len2; j++) {
        c = ((c << 2) | q[j]) & 1023;
        s.nxt[j] = s.head[c];
        s.head[c] = (int16_t)j;
    }
    memset(s.cnt, 0, nd * sizeof(uint16_t));
    c = (ref[0] << 6) | (ref[1] << 4) | (ref[2] << 2) | ref[3];
    for (int i = 4; i < len1; i++) {
        c = ((c << 2) | ref[i]) & 1023;
        for (int j = s.head[c]; j >= 0; j = s.nxt[j]) {
            const int d = i - j + off;
            if (!s.cnt[d]) s.minrow[d] = (int16_t)(i - 4);
            s.cnt[d]++;
        }
    }
    int best = -1, f = 0;
    for (int d = 0; d < nd; d++) {
        f = (int)s.cnt[d] - 1 + (f > 0 ? f : 0);
        s.fwd[d] = (int16_t)f;
        if (f > best) best = f;
    }
    if (5 + best < MINSC) return B1;
    int b = 0;
    for (int d = nd - 1; d >= 0; d--) {
        b = (int)s.cnt[d] - 1 + (b > 0 ? b : 0);
        s.bwd[d] = (int16_t)b;
    }
    int lo = len1, hi = -1;
    for (int d = 0; d < nd;) {
        auto bnd = [&](int x) { return 5 + s.fwd[x] + s.bwd[x] - ((int)s.cnt[x] - 1); };
        if (bnd(d) < MINSC) { d++; continue; }
        int ub = 0, i0 = len1, dmax = -1;
        while (d < nd && bnd(d) >= MINSC) {
            ub = std::max(ub, bnd(d));
            if (s.cnt[d]) { i0 = std::min(i0, (int)s.minrow[d]); dmax = d; }
            d++;
        }
        if (dmax < 0) continue;
        lo = std::min(lo, i0);
        hi = std::max(hi, (dmax - off) + quanta - 1 + std::max(0, ub - 21));
    }
    if (hi < 0) return B1;
    hb = std::max(0, lo);
    he = std::min(len1 - 1, hi);
    return B2;
}

// ---- lean_neon: NEON rewrite of lean() with identical (Kind, hb, he) on every input ----
//
// Design (aarch64):
//  1. N check (vector OR) while copying ref and query into zero-padded buffers, so the vector
//     5-mer code passes never read past either sequence.
//  2. Query table tab[code] = (off - j_last) | MULTI (0x2000 if the code occurs more than once in
//     the query), or ABSENT (0x4040, 2 KB memset). nxt[j] chains earlier occurrences (as lean).
//  3. Per ref row, 16 rows at a time: 5-mer codes in NEON, lanes moved to GPRs as 64-bit words,
//     four scalar table loads merged per word, one vector add/store per 8 rows:
//     dd[r] = r + tab[code[r]] = diagonal of the row's representative hit (latest query
//     occurrence); >= 0x4000 when the code is absent, bit 0x2000 when it is multi.
//  4. Vector pass over dd: normalise (absent -> DUMMY), find starts/ends of runs of equal diagonal
//     (a run of L consecutive rows on diagonal d = L hits on d) and multi rows; left-pack their row
//     positions into lists with a TBL shuffle table. No data-dependent branches so far.
//  5. Per run: cnt[d] += L, minrow[d] = first touch (select). Remaining occurrences of multi codes
//     in "layers": layer L holds each row whose code occurs >= L times, with its L-th latest query
//     position; within a layer consecutive rows on one diagonal collapse to runs again (vector
//     boundary detection + left-pack), minrow by min. There are no per-hit read-modify-write
//     chains through cnt[] (lean's cost), and no per-row branches (lean's mispredicts).
//  6. Kadane as scans: P = prefix sum of cnt-1, PM(d) = min(0, P[0..d-1]) (log-step int16 scans,
//     block-local, plus a carry), fwd = P - PM, max fwd < 14 -> B1; SX(d) = max_{b>=d} P[b]
//     (backward scan); bnd(d) = 5 + SX(d) - PM(d) (== 5 + fwd + bwd - v). Exact while
//     |P| < 32768, guaranteed by hits <= 32000; above that we fall back to lean().
//  7. Components: bitsets of (bnd >= 19) and (bnd >= 19 && cnt > 0); lo = masked vector min of
//     minrow; per component (few) ub = vector max of bnd, dmax = highest hit bit.
#if defined(__aarch64__)
struct NeonScratch {
    static const int CAP = 4096 + 64;
    alignas(16) uint8_t rbuf[CAP + 64];   // 4 zero bytes, ref, zero padding
    alignas(16) uint8_t qbuf[CAP + 64];   // 4 zero bytes, query, zero padding
    alignas(16) uint16_t qcode[CAP];
    alignas(16) uint16_t RA[CAP], RB[CAP], D[CAP], B[CAP];  // multi-occurrence layers
    alignas(16) int16_t JA[CAP], JB[CAP], Jn[CAP];
    alignas(16) int16_t tab[1024];
    alignas(16) int16_t nxt[CAP];
    alignas(16) uint16_t dd[CAP];
    alignas(16) uint16_t S[CAP], E[CAP];                    // run starts / ends (rows)
    alignas(16) uint16_t cnt[CAP];
    alignas(16) int16_t minrow[CAP];
    alignas(16) int16_t P[CAP], PM[CAP], bnd[CAP];
    alignas(16) uint64_t mw[CAP / 64 + 2], hw[CAP / 64 + 2];  // bitsets (written bytewise)
    alignas(16) uint8_t shuf[256][16];   // left-pack shuffles for 8 x u16 lanes
    uint8_t pc[256];
    // Query cache: the query table (qbuf/qcode/tab/nxt/qn_at) depends only on the oriented mate,
    // which is the same for every anchor rescued with that mate and strand. Rebuilt only when the
    // query bytes differ from the last call's.
    alignas(16) uint16_t qcnt[1024];
    alignas(16) uint16_t qn_at[CAP];      // occurrences of code(j) at or before j (total at the last)
    uint8_t qcache[CAP];
    int qlen_c = -1;
    bool q_has_n = false;
    /* Set to nd when the last lean_neon() call returned B2 from the NEON path, so cnt / minrow /
     * bnd / mw hold that call's per-diagonal arrays (rescue_band.h reads them to derive diagonal
     * components at any threshold); -1 otherwise (other kinds, or the scalar lean() fallback). */
    int view_nd = -1;
    Scratch fallback;
    NeonScratch()
    {
        for (int m = 0; m < 256; m++) {
            int n = 0;
            for (int k = 0; k < 8; k++)
                if (m >> k & 1) { shuf[m][2 * n] = (uint8_t)(2 * k); shuf[m][2 * n + 1] = (uint8_t)(2 * k + 1); n++; }
            for (int k = n; k < 8; k++) shuf[m][2 * k] = shuf[m][2 * k + 1] = 0xFF;
            pc[m] = (uint8_t)n;
        }
        memset(rbuf, 0, sizeof rbuf);
        memset(qbuf, 0, sizeof qbuf);
        memset(cnt, 0, sizeof cnt);
        memset(minrow, 0, sizeof minrow);
    }
};

static inline uint16_t *neon_pack(uint16_t *out, uint16x8_t pos, unsigned m, const NeonScratch &s)
{
    const uint8x16_t sh = vld1q_u8(s.shuf[m]);
    vst1q_u16(out, vreinterpretq_u16_u8(vqtbl1q_u8(vreinterpretq_u8_u16(pos), sh)));
    return out + s.pc[m];
}

// First index >= from with bit set (inv=0) / clear (inv=~0) in the bitset w, or n if none.
static inline int neon_next(const uint64_t *w, int from, int n, uint64_t inv)
{
    if (from >= n) return n;
    int k = from >> 6;
    uint64_t x = (w[k] ^ inv) & (~0ull << (from & 63));
    const int nw = (n + 63) >> 6;
    while (!x) { if (++k >= nw) return n; x = w[k] ^ inv; }
    const int r = (k << 6) + __builtin_ctzll(x);
    return r < n ? r : n;
}

// Highest set bit in [a, b) of bitset w, or -1.
static inline int neon_last(const uint64_t *w, int a, int b)
{
    if (b <= a) return -1;
    int k = (b - 1) >> 6;
    uint64_t x = w[k] & (~0ull >> (63 - ((b - 1) & 63)));
    for (;;) {
        if (x) { const int r = (k << 6) + 63 - __builtin_clzll(x); return r >= a ? r : -1; }
        if (--k < (a >> 6)) return -1;
        x = w[k];
    }
}

// 5-mer codes of rows b..b+15, where buf[k + 4] is sequence position k (buf = padded + b).
static inline void neon_codes16(const uint8_t *buf, uint16_t *out)
{
    const uint8x16_t a0 = vld1q_u8(buf), a1 = vld1q_u8(buf + 1), a2 = vld1q_u8(buf + 2),
                     a3 = vld1q_u8(buf + 3), a4 = vld1q_u8(buf + 4);
    uint8x16_t lo = vorrq_u8(vshlq_n_u8(a1, 6), vshlq_n_u8(a2, 4));
    lo = vorrq_u8(lo, vorrq_u8(vshlq_n_u8(a3, 2), a4));
    vst1q_u8((uint8_t *)out, vzip1q_u8(lo, a0));
    vst1q_u8((uint8_t *)out + 16, vzip2q_u8(lo, a0));
}

// max_hits: return FULL when the window and query share more 5-mer hits than this (the gate;
// same count as summing per-code query occurrences over the window rows), decided before the
// repeated-code layers and the Kadane scans.
static inline Kind lean_neon(const Job &jb, NeonScratch &s, int &hb, int &he, int max_hits = 1 << 30)
{
    const uint8_t *ref = jb.ref, *q = jb.qry;
    const int len1 = jb.len1, len2 = jb.len2;
    hb = he = -1;
    s.view_nd = -1;
    // ---- 1. N check (+ copy ref and query into zero-padded buffers) ----
    if (__builtin_expect(len1 + 64 > NeonScratch::CAP || len2 + 64 > NeonScratch::CAP, 0))
        return SCALAR;
    uint8x16_t ov = vdupq_n_u8(0);
    uint8_t orv = 0;
    {
        uint8_t *rb = s.rbuf + 4;
        int i = 0;
        for (; i + 16 <= len1; i += 16) { const uint8x16_t x = vld1q_u8(ref + i); ov = vorrq_u8(ov, x); vst1q_u8(rb + i, x); }
        for (; i < len1; i++) { orv |= ref[i]; rb[i] = ref[i]; }
        vst1q_u8(rb + len1, vdupq_n_u8(0)); vst1q_u8(rb + len1 + 16, vdupq_n_u8(0));
    }
    const int quanta = ((len2 + 15) / 16) * 16, off = quanta, nd = len1 + quanta + 1;
    // ---- 2. query table (cached per oriented query): copy + N flag, tab[code] = (off - j_last)
    //      | MULTI, nxt[] = earlier occurrences, qn_at[j] = occurrences of code(j) up to j ----
    if (len2 != s.qlen_c || memcmp(q, s.qcache, (size_t)len2) != 0) {
        uint8_t *qb = s.qbuf + 4;
        uint8x16_t qv = vdupq_n_u8(0);
        uint8_t qor = 0;
        int j = 0;
        for (; j + 16 <= len2; j += 16) { const uint8x16_t x = vld1q_u8(q + j); qv = vorrq_u8(qv, x); vst1q_u8(qb + j, x); }
        for (; j < len2; j++) { qor |= q[j]; qb[j] = q[j]; }
        vst1q_u8(qb + len2, vdupq_n_u8(0)); vst1q_u8(qb + len2 + 16, vdupq_n_u8(0));
        s.q_has_n = ((qor | vmaxvq_u8(qv)) & 0xFC) != 0;
        memcpy(s.qcache, q, (size_t)len2);
        s.qlen_c = len2;
        if (!s.q_has_n && len2 >= 5) {
            memset(s.tab, 0x40, sizeof s.tab);           // ABSENT = 0x4040
            memset(s.qcnt, 0, sizeof s.qcnt);
            for (int b = 0; b < len2; b += 16) neon_codes16(s.qbuf + b, s.qcode + b);
            for (int jj = 4; jj < len2; jj++) {
                const int c = s.qcode[jj];
                const int old = (uint16_t)s.tab[c];
                const bool absent = old >= 0x4000;
                s.nxt[jj] = absent ? (int16_t)-1 : (int16_t)(off - (old & 0x1FFF));
                s.tab[c] = (int16_t)((off - jj) | (absent ? 0 : 0x2000));
                s.qn_at[jj] = ++s.qcnt[c];
            }
        }
    }
    orv |= vmaxvq_u8(ov);
    if ((orv & 0xFC) || s.q_has_n) return FULL;
    if (len1 < 5 || len2 < 5) return B1;
    if (__builtin_expect(nd + 32 > NeonScratch::CAP, 0)) return SCALAR;

    static const uint16_t wl[8] = {1, 2, 4, 8, 16, 32, 64, 128};
    static const uint16_t wh[8] = {256, 512, 1024, 2048, 4096, 8192, 16384, 32768};
    static const uint16_t io[8] = {0, 1, 2, 3, 4, 5, 6, 7};
    const uint16x8_t WL = vld1q_u16(wl), WH = vld1q_u16(wh), IO = vld1q_u16(io), k8 = vdupq_n_u16(8);

    // ---- 3+4. per ref row: representative diagonal, runs of equal diagonal, multi rows ----
    // (a) dd[r] = r + tab[code[r]] for 16 rows at a time: 5-mer codes in NEON, lanes moved to
    //     GPRs as 64-bit words, four table loads merged per word, one vector add/store per 8
    //     rows (>= 0x4000: code absent from the query; bit 0x2000: code occurs more than once).
    {
        const uint16_t *__restrict tab = (const uint16_t *)s.tab;
        static const uint16_t io16[16] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15};
        const uint16x8_t io0 = vld1q_u16(io16), io1 = vld1q_u16(io16 + 8);
        for (int b = 0; b < len1; b += 16) {
            const uint8_t *buf = s.rbuf + b;
            const uint8x16_t a0 = vld1q_u8(buf), a1 = vld1q_u8(buf + 1), a2 = vld1q_u8(buf + 2),
                             a3 = vld1q_u8(buf + 3), a4 = vld1q_u8(buf + 4);
            uint8x16_t lo = vorrq_u8(vshlq_n_u8(a1, 6), vshlq_n_u8(a2, 4));
            lo = vorrq_u8(lo, vorrq_u8(vshlq_n_u8(a3, 2), a4));
            const uint64x2_t c0 = vreinterpretq_u64_u8(vzip1q_u8(lo, a0));
            const uint64x2_t c1 = vreinterpretq_u64_u8(vzip2q_u8(lo, a0));
            const uint64_t w[4] = {vgetq_lane_u64(c0, 0), vgetq_lane_u64(c0, 1),
                                   vgetq_lane_u64(c1, 0), vgetq_lane_u64(c1, 1)};
            uint64_t o[4];
            for (int k = 0; k < 4; k++) {
                const uint64_t x = w[k];
                o[k] = (uint64_t)tab[x & 1023] | ((uint64_t)tab[(x >> 16) & 1023] << 16) |
                       ((uint64_t)tab[(x >> 32) & 1023] << 32) | ((uint64_t)tab[x >> 48] << 48);
            }
            const uint16x8_t bb = vdupq_n_u16((uint16_t)b);
            const uint16x8_t v0 = vreinterpretq_u16_u64(vcombine_u64(vcreate_u64(o[0]), vcreate_u64(o[1])));
            const uint16x8_t v1 = vreinterpretq_u16_u64(vcombine_u64(vcreate_u64(o[2]), vcreate_u64(o[3])));
            vst1q_u16(s.dd + b, vaddq_u16(v0, vaddq_u16(io0, bb)));
            vst1q_u16(s.dd + b + 8, vaddq_u16(v1, vaddq_u16(io1, bb)));
        }
        uint16_t *dd = s.dd;
        dd[0] = dd[1] = dd[2] = dd[3] = 0x4000;             // rows 0..3 have no 5-mer
        vst1q_u16(dd + len1, vdupq_n_u16(0x4000)); vst1q_u16(dd + len1 + 8, vdupq_n_u16(0x4000));
    }
    // (b) normalise (absent -> DUMMY), then left-pack run starts, run ends (one past) and multi
    //     rows into position lists: a run of L rows on diagonal d is L hits on d.
    const uint16_t DUMMY = 0x1FFF;
    uint16_t *sp = s.S, *ep = s.E, *mp = s.RA;
    {
        const uint16x8_t k1fff = vdupq_n_u16(0x1FFF), k4000 = vdupq_n_u16(0x4000), k2000 = vdupq_n_u16(0x2000);
        uint16x8_t pos = IO;
        uint16x8_t prev = vdupq_n_u16(DUMMY);
        for (int b = 0; b <= len1; b += 8, pos = vaddq_u16(pos, k8)) {
            const uint16x8_t x = vld1q_u16(s.dd + b);
            const uint16x8_t absent = vcgeq_u16(x, k4000);
            const uint16x8_t multi = vtstq_u16(x, k2000);
            const uint16x8_t xn = vorrq_u16(vandq_u16(x, k1fff), vandq_u16(absent, k1fff));
            vst1q_u16(s.dd + b, xn);
            const uint16x8_t xp = vextq_u16(prev, xn, 7);
            prev = xn;
            const uint16x8_t ne = vmvnq_u16(vceqq_u16(xn, xp));
            const uint16x8_t st = vbicq_u16(ne, absent);                     // new non-dummy run
            const uint16x8_t en = vbicq_u16(ne, vceqq_u16(xp, k1fff));       // one past a run's end
            const unsigned se = vaddvq_u16(vorrq_u16(vandq_u16(st, WL), vandq_u16(en, WH)));
            const unsigned mm = vaddvq_u16(vandq_u16(multi, WL));
            sp = neon_pack(sp, pos, se & 0xFF, s);
            ep = neon_pack(ep, pos, se >> 8, s);
            mp = neon_pack(mp, pos, mm, s);
        }
    }

    // ---- 5. accumulate hits: runs of representative hits, then multi-occurrence layers ----
    const int ndr = (nd + 7) & ~7;
    memset(s.cnt, 0, (ndr + 8) * sizeof(uint16_t));
    int hits = 0;
    {
        const int nr = (int)(sp - s.S);
        const uint16_t *__restrict S = s.S, *__restrict E = s.E, *__restrict xn = s.dd;
        uint16_t *__restrict cnt = s.cnt;
        int16_t *__restrict minrow = s.minrow;
        for (int k = 0; k < nr; k++) {
            const int s0 = S[k], len = E[k] - s0, d = xn[s0];
            const int c = cnt[d];
            const int16_t mr = minrow[d];
            minrow[d] = c ? mr : (int16_t)(s0 - 4);           // runs arrive in row order
            cnt[d] = (uint16_t)(c + len);
            hits += len;
        }
        // Layer L >= 2 holds, for every row whose code occurs >= L times in the query, the L-th
        // latest occurrence J. Within a layer, consecutive rows on one diagonal (tandem repeats,
        // repeated k-mers inside a matching segment) again collapse to runs.
        int n = (int)(mp - s.RA);
        // Gate: the total hit count is known now (each repeated-code row adds its code's
        // remaining occurrences), before the costly layers and scans.
        if (max_hits < (1 << 30)) {
            long total = hits;
            for (int k = 0; k < n; k++) { const int r = s.RA[k]; total += s.qn_at[r + off - xn[r]] - 1; }
            if (total > max_hits) return FULL;
        }
        uint16_t *R = s.RA, *R2 = s.RB;
        int16_t *J = s.JA, *J2 = s.JB, *Jn = s.Jn;
        const int16_t *__restrict nxt = s.nxt;
        for (int k = 0; k < n; k++) { const int r = R[k]; J[k] = nxt[r + off - xn[r]]; }
        const uint16x8_t koff = vdupq_n_u16((uint16_t)off), one16 = vdupq_n_u16(1);
        while (n > 0) {
            hits += n;
            if (__builtin_expect(hits > 32000, 0)) return SCALAR;
            // diagonals and run boundaries (as list indices; k = 0 is always a boundary)
            uint16_t *bp = s.B;
            uint16x8_t pR = vdupq_n_u16(0), pD = vdupq_n_u16(0xFFFF), idx = IO;
            const uint16x8_t vn = vdupq_n_u16((uint16_t)n);
            for (int k0 = 0; k0 < n; k0 += 8, idx = vaddq_u16(idx, k8)) {
                const uint16x8_t r = vld1q_u16(R + k0);
                const uint16x8_t d = vsubq_u16(vaddq_u16(r, koff), vreinterpretq_u16_s16(vld1q_s16(J + k0)));
                vst1q_u16(s.D + k0, d);
                const uint16x8_t rp = vextq_u16(pR, r, 7), dp = vextq_u16(pD, d, 7);
                pR = r; pD = d;
                const uint16x8_t cont = vandq_u16(vceqq_u16(r, vaddq_u16(rp, one16)), vceqq_u16(d, dp));
                const uint16x8_t bd = vandq_u16(vmvnq_u16(cont), vcltq_u16(idx, vn));
                bp = neon_pack(bp, idx, vaddvq_u16(vandq_u16(bd, WL)), s);
            }
            const int nbd = (int)(bp - s.B);
            s.B[nbd] = (uint16_t)n;
            for (int t = 0; t < nbd; t++) {
                const int b0 = s.B[t], len = s.B[t + 1] - b0, d = s.D[b0], row = R[b0];
                const int c = cnt[d], mr = minrow[d];
                minrow[d] = (int16_t)(c ? std::min(mr, row - 4) : row - 4);
                cnt[d] = (uint16_t)(c + len);
            }
            // next layer: rows whose occurrence chain continues
            for (int k = 0; k < n; k++) Jn[k] = nxt[J[k]];
            uint16_t *rp2 = R2, *jp2 = (uint16_t *)J2;
            for (int k0 = 0; k0 < n; k0 += 8) {
                const int16x8_t jn = vld1q_s16(Jn + k0);
                const uint16x8_t lanes = vcltq_u16(vaddq_u16(IO, vdupq_n_u16((uint16_t)k0)), vn);
                const unsigned m = vaddvq_u16(vandq_u16(vandq_u16(vcgezq_s16(jn), lanes), WL));
                const uint8x16_t sh = vld1q_u8(s.shuf[m]);
                vst1q_u16(rp2, vreinterpretq_u16_u8(vqtbl1q_u8(vreinterpretq_u8_u16(vld1q_u16(R + k0)), sh)));
                vst1q_u16(jp2, vreinterpretq_u16_u8(vqtbl1q_u8(vreinterpretq_u8_s16(jn), sh)));
                rp2 += s.pc[m]; jp2 += s.pc[m];
            }
            n = (int)(rp2 - R2);
            std::swap(R, R2); std::swap(J, J2);
        }
    }
    if (hits == 0) return B1;
    if (__builtin_expect(hits > 32000, 0)) return SCALAR;

    // ---- 6a. forward scans: P (inclusive prefix of cnt-1), PM (exclusive running min, with 0) ----
    const int16x8_t zero = vdupq_n_s16(0), one = vdupq_n_s16(1), big = vdupq_n_s16(32767),
                    small = vdupq_n_s16(-32768);
    {
        // Block-local scans are independent of the carries; only one add (P) and one min (PM)
        // sit on the loop-carried chains.
        int16x8_t carryP = zero, carryM = zero, bestv = small;
        for (int d0 = 0; d0 < ndr; d0 += 8) {
            int16x8_t l = vsubq_s16(vreinterpretq_s16_u16(vld1q_u16(s.cnt + d0)), one);
            l = vaddq_s16(l, vextq_s16(zero, l, 7));
            l = vaddq_s16(l, vextq_s16(zero, l, 6));
            l = vaddq_s16(l, vextq_s16(zero, l, 4));            // local inclusive prefix
            int16x8_t lm = vminq_s16(l, vextq_s16(big, l, 7));
            lm = vminq_s16(lm, vextq_s16(big, lm, 6));
            lm = vminq_s16(lm, vextq_s16(big, lm, 4));          // local running min of l
            const int16x8_t p = vaddq_s16(l, carryP);
            const int16x8_t m = vminq_s16(vaddq_s16(lm, carryP), carryM);  // incl. running min
            const int16x8_t pm = vextq_s16(carryM, m, 7);
            carryM = vminq_s16(carryM, vaddq_s16(carryP, vdupq_laneq_s16(lm, 7)));
            carryP = vaddq_s16(carryP, vdupq_laneq_s16(l, 7));
            bestv = vmaxq_s16(bestv, vsubq_s16(p, pm));
            vst1q_s16(s.P + d0, p);
            vst1q_s16(s.PM + d0, pm);
        }
        if (5 + vmaxvq_s16(bestv) < MINSC) return B1;
    }

    // ---- 6b. backward scan: SX = suffix max of P; bnd = 5 + SX - PM; bitsets; lo ----
    const int nb = ndr >> 3;
    uint8_t *mb = (uint8_t *)s.mw, *hbb = (uint8_t *)s.hw;
    memset(mb + nb, 0, 16); memset(hbb + nb, 0, 16);
    int lo;
    {
        const int16x8_t five = vdupq_n_s16(5), minsc = vdupq_n_s16(MINSC);
        int16x8_t carryS = small, lov = big;
        for (int d0 = ndr - 8; d0 >= 0; d0 -= 8) {
            const int16x8_t p = vld1q_s16(s.P + d0);
            int16x8_t sx = vmaxq_s16(p, vextq_s16(p, small, 1));
            sx = vmaxq_s16(sx, vextq_s16(sx, small, 2));
            sx = vmaxq_s16(sx, vextq_s16(sx, small, 4));        // local suffix max
            const int16x8_t sxl = sx;
            sx = vmaxq_s16(sx, carryS);
            carryS = vmaxq_s16(carryS, vdupq_laneq_s16(sxl, 0));
            const int16x8_t bd = vaddq_s16(five, vsubq_s16(sx, vld1q_s16(s.PM + d0)));
            vst1q_s16(s.bnd + d0, bd);
            const uint16x8_t in = vcgeq_s16(bd, minsc);
            const uint16x8_t c = vld1q_u16(s.cnt + d0);
            const uint16x8_t ih = vandq_u16(in, vtstq_u16(c, c));
            lov = vminq_s16(lov, vbslq_s16(ih, vld1q_s16(s.minrow + d0), big));
            const unsigned bits = vaddvq_u16(vorrq_u16(vandq_u16(in, WL), vandq_u16(ih, WH)));
            mb[d0 >> 3] = (uint8_t)bits;
            hbb[d0 >> 3] = (uint8_t)(bits >> 8);
        }
        if (nd & 7) mb[nd >> 3] &= (uint8_t)((1u << (nd & 7)) - 1);  // drop pad lanes
        lo = vminvq_s16(lov);
    }

    // ---- 7. components ----
    int hi = -1;
    for (int d = neon_next(s.mw, 0, nd, 0); d < nd;) {
        const int a = d, b = neon_next(s.mw, a, nd, ~0ull);
        const int dmax = neon_last(s.hw, a, b);
        if (dmax >= 0) {
            int ub = -32768, k = a;
            for (; k < b && (k & 7); k++) ub = std::max(ub, (int)s.bnd[k]);
            int16x8_t mv = small;
            for (; k + 8 <= b; k += 8) mv = vmaxq_s16(mv, vld1q_s16(s.bnd + k));
            ub = std::max(ub, (int)vmaxvq_s16(mv));
            for (; k < b; k++) ub = std::max(ub, (int)s.bnd[k]);
            hi = std::max(hi, (dmax - off) + quanta - 1 + std::max(0, ub - 21));
        }
        d = neon_next(s.mw, b, nd, 0);
    }
    if (hi < 0) return B1;
    hb = std::max(0, lo);
    he = std::min(len1 - 1, hi);
    s.view_nd = nd;
    return B2;
}
#else
struct NeonScratch { Scratch fallback; };
static inline Kind lean_neon(const Job &jb, NeonScratch &s, int &hb, int &he) { return lean(jb, s.fallback, hb, he); }
#endif

}  // namespace rescue_prune_neon

#endif
