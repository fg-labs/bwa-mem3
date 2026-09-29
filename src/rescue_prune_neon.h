/* NEON implementation of the exact rescue-pruning filter (see rescue_prune.h).
 *
 * rescue_prune_neon::lean_neon() returns exactly the same (kind, hb, he) as the scalar
 * rescue_prune_neon::lean() on every input -- verified with 0 disagreements on 1.97 M real
 * wgs/wes rescue jobs and ~70 K adversarial fuzz jobs (incl. N, tiny lengths, poly-A,
 * int16-range stress) -- and runs ~2x faster on Graviton4. The threshold is fixed at 19
 * (min_seed_len * a at defaults); the caller dispatches here only in that case. lean() is the
 * int16 reference the NEON rewrite was derived from, kept as a test oracle.
 *
 * Overview and gates: docs/src/developer-guide/rescue-pruning.md. */
#ifndef BWA_MEM3_RESCUE_PRUNE_NEON_H
#define BWA_MEM3_RESCUE_PRUNE_NEON_H

#include <algorithm>
#include <cstdint>
#include <cstring>
#include <vector>
#include "kswv_quantum.h"
#if defined(__aarch64__)
#include <arm_neon.h>
#include "neon_transpose.h"
#endif

namespace rescue_prune_neon {

static const int MINSC = 19;
// The NEON filter's hit capacity: every int16 sum it forms stays in range at or below this many
// hits (steps 6 and 7); above it (reachable only with the hit gate opened past it) it returns
// FALLBACK and the int32 scalar filter decides.
static const int HIT_CAP = 28000;

struct Job {
    int len1, len2, score, te, hb, he;
    const uint8_t *ref, *qry;
};

// FALLBACK: the job exceeds the NEON path's capacity (window > ~4000 rows, or more than HIT_CAP
// hits); the caller must use the int32 scalar filter. (lean() below is int16 and
// fixed-size, so it must NOT be used as the fallback.)
enum Kind { FULL = 0, B1 = 1, B2 = 2, FALLBACK = 3 };

struct Scratch {
    int16_t head[1024];
    int16_t nxt[1024];
    uint16_t cnt[4096];
    int16_t minrow[4096];
    int16_t fwd[4096], bwd[4096];
    uint8_t qcnt[1024];
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
    const int quanta = kswv_query_quantum8(len2), off = quanta, nd = len1 + quanta + 1;

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
        hi = std::max(hi, (dmax - off) + quanta - 1 + std::max(0, ub - MINSC - 2));
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
//  2. Query table, cached per oriented query: tab[code] = (off - j_last) | MULTI (0x2000 if the
//     code occurs more than once) | occurrences << 16, a 1024-bit presence bitmap pres[], and
//     nxt[j] chaining earlier occurrences (as lean).
//  3. Hit rows, 16 at a time: 5-mer codes and a presence-bitmap lookup (two 64-byte TBLs) in NEON;
//     the rows whose code occurs in the query (~1/3 of a window) and their codes are left-packed
//     into lists with a TBL shuffle table. No scalar work per window row.
//  4. Layer 1, 8 hit rows at a time: four scalar table loads per 64-bit word of codes give each
//     row's representative hit (latest query occurrence), diagonal d = row + off - j_last; runs of
//     consecutive rows on one diagonal (a run of L rows on d = L hits on d), the multi rows, and
//     the exact hit total for the gate (sum of occurrences) come out of the same vector pass.
//  5. Per run: cnt[d] += L, minrow[d] = first touch (layer-1 runs are applied last to first, so a
//     plain store suffices). Remaining occurrences of multi codes in "layers": layer L holds each
//     row whose code occurs >= L times, with its L-th latest query position; within a layer
//     consecutive rows on one diagonal collapse to runs again, minrow by min, and the next layer
//     is gathered in the same pass. There are no per-hit read-modify-write chains through cnt[]
//     (lean's cost), and no per-row branches (lean's mispredicts).
//  6. Kadane as scans: P = prefix sum of cnt-1, PM(d) = min(0, P[0..d-1]), fwd = P - PM, max fwd
//     < 14 -> B1; SX(d) = max_{b>=d} P[b]; bnd(d) = 5 + SX(d) - PM(d) (== 5 + fwd + bwd - v).
//     The diagonals are cut into 8 segments (one per lane) and transposed 8x8 in registers, so
//     each scan step is one vertical op for 8 diagonals; segment carries join them. Exact while
//     |P| < 32768, guaranteed by hits <= HIT_CAP; above that it returns FALLBACK and the caller
//     runs the int32 scalar filter (rescue_prune_window_scalar), never the int16 lean().
//  7. Components: bitsets of (bnd >= 19) and (bnd >= 19 && cnt > 0) built in the backward scan;
//     per component with a hit (few) ub = vector max of bnd, dmax = highest hit bit, and
//     lo = masked vector min of minrow.
#if defined(__aarch64__)
/* A diagonal component at MINSC with a hit, as step 7 finds it (shifted diagonal indices x = d + off):
 * the mw run [a, b), ub = max bnd over it, i0 = min minrow over its hit diagonals, dmax = its
 * highest hit diagonal. rescue_band.cpp reuses the list instead of rescanning the view at MINSC. */
struct Comp { int32_t a, b, ub, i0, dmax; };
struct NeonScratch {
    static const int CAP = 4096 + 64;
    alignas(16) uint8_t rbuf[CAP + 64];   // 4 zero bytes, ref, zero padding
    alignas(16) uint8_t qbuf[CAP + 64];   // 4 zero bytes, query, zero padding
    alignas(16) uint16_t qcode[CAP];
    alignas(16) uint16_t RA[CAP], RB[CAP];                  // multi-occurrence layers: rows
    alignas(16) int16_t JA[CAP], JB[CAP];                   // ... and query positions
    alignas(16) uint16_t B[CAP], BD[CAP], BR[CAP];          // runs: list index, diagonal, row
    alignas(16) uint32_t tab[1024];       // per code: (off - j_last) | MULTI, occurrences << 16
    alignas(16) uint8_t pres[128];        // bitmap of the 5-mer codes present in the query
    alignas(16) int16_t nxt[CAP];
    alignas(16) uint16_t PR[CAP], PC[CAP];  // hit rows and their 5-mer codes
    alignas(16) uint16_t cnt[CAP + 8];
    alignas(16) int16_t minrow[CAP];
    // P (after 8 zeros: P before each segment's first diagonal), PM: segment-transposed (step 6)
    alignas(16) int16_t P[8 + CAP], PM[CAP];
    // bnd is readable VIEW_PAD entries past any nd the filter accepts (nd + 32 <= CAP):
    // rescue_band.cpp reads it in whole 64-diagonal words (static_assert there).
    static const int VIEW_PAD = 64;
    alignas(16) int16_t bnd[CAP + VIEW_PAD];
    alignas(16) uint64_t mw[CAP / 64 + 3], hw[CAP / 64 + 3];  // bitsets
    alignas(16) uint8_t MH[(CAP / 64 + 8) * 16];   // per backward step: mw bytes (lane = segment), hw bytes
    alignas(16) uint8_t shuf[256][16];   // left-pack shuffles for 8 x u16 lanes
    /* Step 7's components, in diagonal order: the first min(ncomp, COMP_CAP) are stored; ncomp
     * counts all of them. Valid after a B2 return, with the rest of the view. */
    static const int COMP_CAP = 128;
    Comp comps[COMP_CAP];
    int ncomp = 0;
    uint8_t pc[256];
    // Query cache: the query table (qbuf/qcode/tab/pres/nxt) depends only on the oriented mate,
    // which is the same for every anchor rescued with that mate and strand. Rebuilt only when the
    // query bytes differ from the last call's.
    uint8_t qcache[CAP];
    int qlen_c = -1;
    bool q_has_n = false;
    /* Result of the last lean_neon() call (memo_ok: a NEON decision, not FALLBACK). A call
     * repeating that job exactly -- same window bytes (still in rbuf), query bytes (qcache),
     * lengths and gate -- returns it again; the per-diagonal arrays behind a B2 (cnt, minrow, bnd,
     * mw, hw, comps) were not touched in between, since only lean_neon writes them. */
    bool memo_ok = false;
    int memo_len1 = -1, memo_mh = 0, memo_hb = -1, memo_he = -1;
    Kind memo_kind = FULL;
    uint64_t memo_hits = 0;   // calls answered from the memo (BWA3_RESCUE_PRUNE_STATS, tests)
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
        memset(P, 0, sizeof P);
        memset(MH, 0, sizeof MH);
    }
};

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

// Left-pack store: the u16 lanes of v selected by the shuffle sh (a row of NeonScratch::shuf),
// moved to the front and stored at out (all 8 lanes are written; the caller advances by pc[]).
static inline void neon_pack_store(uint16_t *out, uint16x8_t v, uint8x16_t sh)
{
    vst1q_u16(out, vreinterpretq_u16_u8(vqtbl1q_u8(vreinterpretq_u8_u16(v), sh)));
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
static inline Kind lean_neon_core(const Job &jb, NeonScratch &s, int &hb, int &he, int max_hits)
{
    const uint8_t *ref = jb.ref, *q = jb.qry;
    const int len1 = jb.len1, len2 = jb.len2;
    hb = he = -1;
    // ---- 1. N check (+ copy ref and query into zero-padded buffers) ----
    if (__builtin_expect(len1 + 64 > NeonScratch::CAP || len2 + 64 > NeonScratch::CAP, 0))
        return FALLBACK;
    uint8x16_t ov = vdupq_n_u8(0);
    uint8_t orv = 0;
    {
        uint8_t *rb = s.rbuf + 4;
        int i = 0;
        for (; i + 16 <= len1; i += 16) { const uint8x16_t x = vld1q_u8(ref + i); ov = vorrq_u8(ov, x); vst1q_u8(rb + i, x); }
        for (; i < len1; i++) { orv |= ref[i]; rb[i] = ref[i]; }
        vst1q_u8(rb + len1, vdupq_n_u8(0)); vst1q_u8(rb + len1 + 16, vdupq_n_u8(0));
    }
    const int quanta = kswv_query_quantum8(len2), off = quanta, nd = len1 + quanta + 1;
    // ---- 2. query table (cached per oriented query): copy + N flag, tab[code] = (off - j_last)
    //      | MULTI | occurrences << 16, pres[] = codes present, nxt[] = earlier occurrences ----
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
            memset(s.tab, 0, sizeof s.tab);
            memset(s.pres, 0, sizeof s.pres);
            for (int b = 0; b < len2; b += 16) neon_codes16(s.qbuf + b, s.qcode + b);
            for (int jj = 4; jj < len2; jj++) {
                const int c = s.qcode[jj];
                const uint32_t old = s.tab[c], occ = (old >> 16) + 1;
                s.nxt[jj] = occ == 1 ? (int16_t)-1 : (int16_t)(off - (old & 0x1FFF));
                s.tab[c] = occ << 16 | (uint32_t)(off - jj) | (occ == 1 ? 0 : 0x2000);
                s.pres[c >> 3] |= (uint8_t)(1u << (c & 7));
            }
        }
    }
    orv |= vmaxvq_u8(ov);
    if ((orv & 0xFC) || s.q_has_n) return FULL;
    if (len1 < 5 || len2 < 5) return B1;
    if (__builtin_expect(nd + 32 > NeonScratch::CAP, 0)) return FALLBACK;

    static const uint16_t wl[8] = {1, 2, 4, 8, 16, 32, 64, 128};
    static const uint16_t wh[8] = {256, 512, 1024, 2048, 4096, 8192, 16384, 32768};
    static const uint16_t io[8] = {0, 1, 2, 3, 4, 5, 6, 7};
    const uint16x8_t WL = vld1q_u16(wl), WH = vld1q_u16(wh), IO = vld1q_u16(io), k8 = vdupq_n_u16(8);

    // ---- 3. rows whose 5-mer occurs in the query (hit rows): 16 rows at a time, 5-mer codes and a
    //      presence-bitmap lookup (two 64-byte TBLs) in NEON; left-pack the hit rows and their codes.
    //      Only hit rows (~1/3 of a window) reach the scalar table loads below. ----
    int np;
    {
        const uint8x16x4_t plo = vld1q_u8_x4(s.pres), phi = vld1q_u8_x4(s.pres + 64);
        static const uint8_t p2[16] = {1, 2, 4, 8, 16, 32, 64, 128, 0, 0, 0, 0, 0, 0, 0, 0};
        static const uint8_t bw[16] = {1, 2, 4, 8, 16, 32, 64, 128, 1, 2, 4, 8, 16, 32, 64, 128};
        const uint8x16_t P2 = vld1q_u8(p2), BW = vld1q_u8(bw), k64 = vdupq_n_u8(64), k7 = vdupq_n_u8(7);
        uint16_t *rp = s.PR, *cp = s.PC;
        // rows b..b+15, of which those with bit set in keep are inside the window
        auto block = [&](int b, unsigned keep) {
            const uint8_t *buf = s.rbuf + b;
            const uint8x16_t a0 = vld1q_u8(buf), a1 = vld1q_u8(buf + 1), a2 = vld1q_u8(buf + 2),
                             a3 = vld1q_u8(buf + 3), a4 = vld1q_u8(buf + 4);
            // code & 255 = a1 a2 a3 a4 (2 bits each; bases are 0..3, so shift-insert is exact);
            // code >> 8 = a0
            const uint8x16_t lo = vsliq_n_u8(vsliq_n_u8(vsliq_n_u8(a4, a3, 2), a2, 4), a1, 6);
            const uint8x16_t idx = vsraq_n_u8(vshlq_n_u8(a0, 5), lo, 3);             // code >> 3
            // idx < 64 -> plo; idx >= 64 -> phi (idx - 64 wraps to >= 192 below 64: TBX keeps plo's)
            const uint8x16_t byte = vqtbx4q_u8(vqtbl4q_u8(plo, idx), phi, vsubq_u8(idx, k64));
            const uint8x16_t hit = vtstq_u8(byte, vqtbl1q_u8(P2, vandq_u8(lo, k7)));
            uint8x16_t m = vandq_u8(hit, BW);
            m = vpaddq_u8(m, m); m = vpaddq_u8(m, m); m = vpaddq_u8(m, m);
            const unsigned m16 = vgetq_lane_u16(vreinterpretq_u16_u8(m), 0) & keep;
            const uint16x8_t c0 = vreinterpretq_u16_u8(vzip1q_u8(lo, a0));
            const uint16x8_t c1 = vreinterpretq_u16_u8(vzip2q_u8(lo, a0));
            const uint16x8_t pos = vaddq_u16(IO, vdupq_n_u16((uint16_t)b));
            const unsigned ml = m16 & 0xFF, mh = m16 >> 8;
            const uint8x16_t shl = vld1q_u8(s.shuf[ml]), shh = vld1q_u8(s.shuf[mh]);
            neon_pack_store(rp, pos, shl);
            neon_pack_store(cp, c0, shl);
            rp += s.pc[ml]; cp += s.pc[ml];
            neon_pack_store(rp, vaddq_u16(pos, k8), shh);
            neon_pack_store(cp, c1, shh);
            rp += s.pc[mh]; cp += s.pc[mh];
        };
        // rows 0..3 have no 5-mer; the last block may extend past the window
        const int nfull = len1 >> 4;
        const unsigned tail = (1u << (len1 & 15)) - 1;
        if (nfull == 0) {
            block(0, 0xFFF0u & tail);
        } else {
            block(0, 0xFFF0u);
            for (int b = 16; b < nfull * 16; b += 16) block(b, 0xFFFFu);
            if (tail) block(nfull * 16, tail);
        }
        np = (int)(rp - s.PR);
    }
    if (np == 0) return B1;
    // ---- 4+5. accumulate hits. The table value of a hit row's code is tab[code] = (off - j_last) |
    //      MULTI, where j_last is the code's latest query occurrence and MULTI (0x2000) marks a code
    //      occurring more than once. Layer 1 = every hit row with its latest occurrence: diagonal
    //      d = row + off - j_last; consecutive rows on one diagonal collapse to runs (a run of L
    //      rows on d = L hits on d). Layer L >= 2 holds each row whose code occurs >= L times, with
    //      its L-th latest occurrence, and again collapses to runs. No per-hit read-modify-write
    //      chains through cnt[] and no per-row branches. ----
    int hits = np;
    const uint16x8_t koff = vdupq_n_u16((uint16_t)off), one16 = vdupq_n_u16(1);
    // Run boundaries of a layer's (row, diagonal) list are left-packed as (list index, diagonal,
    // row) into B / BD / BR; k = 0 is always a boundary and B[nbd] = n closes the last run.
    auto pack_runs = [&](uint16_t *&bp, uint16_t *&dp, uint16_t *&rp, uint16x8_t idx, uint16x8_t d,
                         uint16x8_t r, unsigned m) {
        const uint8x16_t sh = vld1q_u8(s.shuf[m]);
        neon_pack_store(bp, idx, sh);
        neon_pack_store(dp, d, sh);
        neon_pack_store(rp, r, sh);
        bp += s.pc[m]; dp += s.pc[m]; rp += s.pc[m];
    };
    int n;   // rows in the current layer >= 2 (R = s.RA, J = s.JA)
    {
        uint16_t *bp = s.B, *bdp = s.BD, *brp = s.BR, *mp = s.RA, *jp = (uint16_t *)s.JA;
        uint16x8_t pR = vdupq_n_u16(0), pD = vdupq_n_u16(0xFFFF), idx = IO;
        const uint16x8_t vn = vdupq_n_u16((uint16_t)np), k1fff = vdupq_n_u16(0x1FFF), k2000 = vdupq_n_u16(0x2000);
        const uint32_t *__restrict tab = s.tab;
        // table entries of four codes (lanes past np hold stale codes: masked below)
        auto look4 = [tab](const uint16_t *c) {
            uint64_t x;
            memcpy(&x, c, 8);
            const uint64_t a = tab[x & 1023] | (uint64_t)tab[(x >> 16) & 1023] << 32;
            const uint64_t b = tab[(x >> 32) & 1023] | (uint64_t)tab[(x >> 48) & 1023] << 32;
            return vreinterpretq_u16_u64(vcombine_u64(vcreate_u64(a), vcreate_u64(b)));
        };
        uint32x4_t extra = vdupq_n_u32(0);   // gate: occurrences - 1 summed over the hit rows
        for (int k0 = 0; k0 < np; k0 += 8, idx = vaddq_u16(idx, k8)) {
            const uint16x8_t r = vld1q_u16(s.PR + k0);
            const uint16x8_t e0 = look4(s.PC + k0), e1 = look4(s.PC + k0 + 4);
            const uint16x8_t t = vuzp1q_u16(e0, e1), occ = vuzp2q_u16(e0, e1);
            const uint16x8_t tv = vandq_u16(t, k1fff);                  // off - j_last
            const uint16x8_t d = vaddq_u16(r, tv);
            const uint16x8_t rp = vextq_u16(pR, r, 7), dp = vextq_u16(pD, d, 7);
            pR = r; pD = d;
            const uint16x8_t valid = vcltq_u16(idx, vn);
            const uint16x8_t cont = vandq_u16(vceqq_u16(r, vaddq_u16(rp, one16)), vceqq_u16(d, dp));
            const uint16x8_t bd = vbicq_u16(valid, cont);
            const uint16x8_t mu = vandq_u16(vtstq_u16(t, k2000), valid);
            extra = vpadalq_u16(extra, vandq_u16(vsubq_u16(occ, one16), valid));
            const unsigned m = vaddvq_u16(vorrq_u16(vandq_u16(bd, WL), vandq_u16(mu, WH)));
            pack_runs(bp, bdp, brp, idx, d, r, m & 0xFF);
            const uint8x16_t sh = vld1q_u8(s.shuf[m >> 8]);
            neon_pack_store(mp, r, sh);
            neon_pack_store(jp, vsubq_u16(koff, tv), sh);
            mp += s.pc[m >> 8]; jp += s.pc[m >> 8];
        }
        n = (int)(mp - s.RA);
        // Gate: the total hit count is known now (each hit row adds its code's occurrences),
        // before anything is accumulated.
        if ((long)np + vaddvq_u32(extra) > max_hits) return FULL;
        // cnt is scanned up to the 8-segment round-up of nd (< nd + 64)
        memset(s.cnt, 0, (((nd + 63) & ~63) + 8) * sizeof(uint16_t));
        // Layer 1 runs, last to first: rows descend, so the last store to minrow[d] is the first
        // touch and needs no read.
        const int nbd = (int)(bp - s.B);
        s.B[nbd] = (uint16_t)np;
        uint16_t *__restrict cnt = s.cnt;
        int16_t *__restrict minrow = s.minrow;
        for (int k = nbd - 1; k >= 0; k--) {
            const int d = s.BD[k];
            cnt[d] = (uint16_t)(cnt[d] + s.B[k + 1] - s.B[k]);
            minrow[d] = (int16_t)(s.BR[k] - 4);
        }
    }
    // Layers >= 2: diagonals and runs as in layer 1, then minrow by min; the next layer (rows
    // whose occurrence chain continues, J = nxt[J]) is gathered and left-packed in the same pass.
    {
        uint16_t *R = s.RA, *R2 = s.RB;
        int16_t *J = s.JA, *J2 = s.JB;
        const int16_t *__restrict nxt = s.nxt;
        for (int k = 0; k < n; k++) J[k] = nxt[J[k]];
        // J < len2 <= 4095; lanes past n hold stale values, so the index is masked to stay in nxt[]
        auto nxt4 = [nxt](const int16_t *j) {
            uint64_t x;
            memcpy(&x, j, 8);
            return (uint64_t)(uint16_t)nxt[x & 4095] | ((uint64_t)(uint16_t)nxt[(x >> 16) & 4095] << 16) |
                   ((uint64_t)(uint16_t)nxt[(x >> 32) & 4095] << 32) | ((uint64_t)(uint16_t)nxt[(x >> 48) & 4095] << 48);
        };
        while (n > 0) {
            hits += n;
            if (__builtin_expect(hits > HIT_CAP, 0)) return FALLBACK;
            uint16_t *bp = s.B, *bdp = s.BD, *brp = s.BR, *rp2 = R2, *jp2 = (uint16_t *)J2;
            uint16x8_t pR = vdupq_n_u16(0), pD = vdupq_n_u16(0xFFFF), idx = IO;
            const uint16x8_t vn = vdupq_n_u16((uint16_t)n);
            for (int k0 = 0; k0 < n; k0 += 8, idx = vaddq_u16(idx, k8)) {
                const uint16x8_t r = vld1q_u16(R + k0);
                const uint16x8_t d = vsubq_u16(vaddq_u16(r, koff), vreinterpretq_u16_s16(vld1q_s16(J + k0)));
                const uint16x8_t rp = vextq_u16(pR, r, 7), dp = vextq_u16(pD, d, 7);
                pR = r; pD = d;
                const uint16x8_t valid = vcltq_u16(idx, vn);
                const uint16x8_t cont = vandq_u16(vceqq_u16(r, vaddq_u16(rp, one16)), vceqq_u16(d, dp));
                const int16x8_t jn = vreinterpretq_s16_u64(vcombine_u64(vcreate_u64(nxt4(J + k0)), vcreate_u64(nxt4(J + k0 + 4))));
                const uint16x8_t keep = vandq_u16(vcgezq_s16(jn), valid);
                const unsigned m = vaddvq_u16(vorrq_u16(vandq_u16(vbicq_u16(valid, cont), WL), vandq_u16(keep, WH)));
                pack_runs(bp, bdp, brp, idx, d, r, m & 0xFF);
                const uint8x16_t sh = vld1q_u8(s.shuf[m >> 8]);
                neon_pack_store(rp2, r, sh);
                neon_pack_store(jp2, vreinterpretq_u16_s16(jn), sh);
                rp2 += s.pc[m >> 8]; jp2 += s.pc[m >> 8];
            }
            const int nbd = (int)(bp - s.B);
            s.B[nbd] = (uint16_t)n;
            uint16_t *__restrict cnt = s.cnt;
            int16_t *__restrict minrow = s.minrow;
            for (int k = 0; k < nbd; k++) {
                const int d = s.BD[k], row = s.BR[k], c = cnt[d], mr = minrow[d];
                minrow[d] = (int16_t)(c ? std::min(mr, row - 4) : row - 4);
                cnt[d] = (uint16_t)(c + s.B[k + 1] - s.B[k]);
            }
            n = (int)(rp2 - R2);
            std::swap(R, R2); std::swap(J, J2);
        }
    }

    // ---- 6. Kadane bounds as scans over the diagonals: P = inclusive prefix sum of cnt - 1,
    //      PM(d) = min(0, P[0..d-1]), SX(d) = max(P[d..]); bnd(d) = 5 + SX(d) - PM(d)
    //      (== 5 + fwd + bwd - (cnt - 1)), max fwd = max(P - PM). Exact while |P| < 32768,
    //      guaranteed by hits <= HIT_CAP.
    //      The diagonals are split into 8 segments of L (a multiple of 8), one per lane, and an 8x8
    //      transpose turns 8 consecutive diagonals of every segment into 8 vectors, so each scan
    //      step is one vertical op for 8 diagonals (instead of log-step shuffles within a vector).
    //      Segments are scanned from zero, then combined through per-segment carries. ----
    const int L = ((nd + 63) >> 6) << 3, NT = 8 * L;
    const int16x8_t big = vdupq_n_s16(32767), small = vdupq_n_s16(-32768), zero = vdupq_n_s16(0);
    int16x8_t CP, Mb, Xa;   // per segment: P before it; min(0, P before it); max P after it
    {
        const int16x8_t one = vdupq_n_s16(1);
        int16x8_t p = zero, rm = big, mx = small, bestl = small;
        int16_t *__restrict TP = s.P + 8, *__restrict TM = s.PM;
        for (int t = 0; t < L; t += 8) {
            int16x8_t c[8];
            for (int k = 0; k < 8; k++) c[k] = vreinterpretq_s16_u16(vld1q_u16(s.cnt + k * L + t));
            neon_transpose8x8_s16(c);
            for (int k = 0; k < 8; k++) {
                p = vaddq_s16(p, vsubq_s16(c[k], one));
                vst1q_s16(TP + (t + k) * 8, p);
                vst1q_s16(TM + (t + k) * 8, rm);            // segment-local exclusive running min
                bestl = vmaxq_s16(bestl, vqsubq_s16(p, rm));
                rm = vminq_s16(rm, p);
                mx = vmaxq_s16(mx, p);
            }
        }
        // carries across lanes (segments): exclusive prefix sum of the totals; exclusive prefix
        // min (with 0) of the segment minima; exclusive suffix max of the segment maxima
        int16x8_t x = p;
        x = vaddq_s16(x, vextq_s16(zero, x, 7));
        x = vaddq_s16(x, vextq_s16(zero, x, 6));
        x = vaddq_s16(x, vextq_s16(zero, x, 4));
        CP = vextq_s16(zero, x, 7);
        int16x8_t g = vaddq_s16(CP, rm);
        g = vminq_s16(g, vextq_s16(big, g, 7));
        g = vminq_s16(g, vextq_s16(big, g, 6));
        g = vminq_s16(g, vextq_s16(big, g, 4));
        Mb = vextq_s16(zero, g, 7);                            // lane 0: min(0, nothing) = 0
        Mb = vminq_s16(Mb, zero);
        const int16x8_t h = vaddq_s16(CP, mx);
        int16x8_t hs = vmaxq_s16(h, vextq_s16(h, small, 1));
        hs = vmaxq_s16(hs, vextq_s16(hs, small, 2));
        hs = vmaxq_s16(hs, vextq_s16(hs, small, 4));
        Xa = vextq_s16(hs, small, 1);
        // max fwd over segment s = max(local best, max P in s - Mb)
        const int16x8_t best = vmaxq_s16(bestl, vsubq_s16(h, Mb));
        if (5 + vmaxvq_s16(best) < MINSC) return B1;
    }
    // Backward: SX, PM and bnd per segment step, and the bitsets mw = {bnd >= 19} and
    // hw = {bnd >= 19 and cnt > 0} as one byte per segment and step (bit k = step k); bnd is
    // transposed back to diagonal order. cnt > 0 is P(d) >= P(d - 1) (their difference is
    // cnt - 1; P before a segment's first diagonal is the zero guard in front of TP).
    //   The per-segment constants move out of the loop:
    //   SX = max(CP + sx, Xa) = CP + max(sx, Xa - CP), PM = min(Mb, CP + TM) = CP + min(Mb - CP, TM),
    //   so bnd = 5 + max(sx, Xa - CP) - min(Mb - CP, TM), with sx started at Xa - CP. Exact in
    //   int16: |P| <= max(nd, hits) and nd < 4200, so Xa - CP and Mb - CP lie in
    //   [-(HIT_CAP + 4200), HIT_CAP + 4200] (Xa = -32768, no later segment, saturates
    //   below every local P >= -L, as the original max with -32768 ignores it); TM = 32767 (no
    //   earlier diagonal in the segment) gives min = Mb - CP as the saturated original gives Mb.
    uint8_t *mb = (uint8_t *)s.mw, *hbb = (uint8_t *)s.hw;
    {
        const int16x8_t five = vdupq_n_s16(5), minsc = vdupq_n_s16(MINSC);
        const int16_t *__restrict TP = s.P + 8, *__restrict TM = s.PM;
        const int sb = L >> 3;   // bitset bytes per segment
        const int16x8_t Xl = vqsubq_s16(Xa, CP), Ml = vsubq_s16(Mb, CP);
        int16x8_t sx = Xl;
        for (int t = L - 8; t >= 0; t -= 8) {
            int16x8_t bd[8];
            uint16x8_t am = vdupq_n_u16(0), ah = vdupq_n_u16(0);
            int16x8_t pn = vld1q_s16(TP + (t + 8) * 8 - 8);   // P at diagonal t + 7 (k = 7)
            for (int k = 7; k >= 0; k--) {
                const int o = (t + k) * 8;
                const int16x8_t pc = pn;
                pn = vld1q_s16(TP + o - 8);                          // P at the diagonal before
                sx = vmaxq_s16(sx, pc);
                bd[k] = vaddq_s16(five, vsubq_s16(sx, vminq_s16(Ml, vld1q_s16(TM + o))));
                const uint16x8_t in = vcgeq_s16(bd[k], minsc);
                am = vsliq_n_u16(in, am, 1);                        // (am << 1) | bit k
                ah = vsliq_n_u16(vandq_u16(in, vcgeq_s16(pc, pn)), ah, 1);
            }
            neon_transpose8x8_s16(bd);
            for (int k = 0; k < 8; k++) vst1q_s16(s.bnd + k * L + t, bd[k]);
            vst1q_u8(s.MH + (t >> 3) * 16, vcombine_u8(vmovn_u16(am), vmovn_u16(ah)));
        }
        // Scatter MH to the bitsets: byte j of segment s is MH row j, lane s (mw) / 8 + s (hw). Eight
        // rows at a time, transposed 8x8 in registers (both halves at once), one 8-byte store per
        // segment and bitset. Groups go last to first: the last group's stores run past its segment
        // (rows >= sb are stale) into the next segment's first bytes, which group 0 rewrites later,
        // or, for segment 7, into the pad bytes cleared below.
        for (int j0 = (sb - 1) & ~7; j0 >= 0; j0 -= 8) {
            const uint8_t *row = s.MH + j0 * 16;
            uint8x16_t r[8], b[8];
            for (int k = 0; k < 8; k++) r[k] = vld1q_u8(row + 16 * k);
            for (int k = 0; k < 8; k += 2) { b[k] = vtrn1q_u8(r[k], r[k + 1]); b[k + 1] = vtrn2q_u8(r[k], r[k + 1]); }
            uint16x8_t c[8];
            for (int k = 0; k < 8; k += 4)
                for (int e = 0; e < 2; e++) {
                    const uint16x8_t x = vreinterpretq_u16_u8(b[k + e]), y = vreinterpretq_u16_u8(b[k + e + 2]);
                    c[k + e] = vtrn1q_u16(x, y); c[k + e + 2] = vtrn2q_u16(x, y);
                }
            uint8x16_t d[8];   // d[s]: segment s (low half mw, high half hw)
            for (int k = 0; k < 4; k++) {
                const uint32x4_t x = vreinterpretq_u32_u16(c[k]), y = vreinterpretq_u32_u16(c[k + 4]);
                d[k] = vreinterpretq_u8_u32(vtrn1q_u32(x, y)); d[k + 4] = vreinterpretq_u8_u32(vtrn2q_u32(x, y));
            }
            for (int k = 0; k < 8; k++) {   // ascending, so a segment's run-over is rewritten by the next
                vst1_u8(mb + k * sb + j0, vget_low_u8(d[k]));
                vst1_u8(hbb + k * sb + j0, vget_high_u8(d[k]));
            }
        }
        // drop the pad diagonals [nd, NT) from mw (hw has none: cnt = 0 there)
        if (nd & 7) mb[nd >> 3] &= (uint8_t)((1u << (nd & 7)) - 1);
        memset(mb + ((nd + 7) >> 3), 0, (NT >> 3) - ((nd + 7) >> 3) + 16);
        memset(hbb + (NT >> 3), 0, 16);
    }

    // ---- 7. components (runs of mw) with a hit: ub = max bnd; dmax = last hw bit; lo = min over
    //      them of minrow on hit diagonals. 8-wide over [a & ~7, b), lanes outside [a, b) masked. ----
    int hi = -1, lo = 32767, nc = 0;
    for (int d = neon_next(s.mw, 0, nd, 0); d < nd;) {
        const int a = d, b = neon_next(s.mw, a, nd, ~0ull);
        const int dmax = neon_last(s.hw, a, b);
        if (dmax >= 0) {
            const int16x8_t IOs = vreinterpretq_s16_u16(IO), va = vdupq_n_s16((int16_t)a), vb = vdupq_n_s16((int16_t)b);
            int16x8_t mx = small, mn = big;
            // one block of 8 diagonals; in = lanes inside [a, b) (all of them for inner blocks)
            auto block = [&](int x, uint16x8_t in) {
                const uint16x8_t c = vld1q_u16(s.cnt + x);
                mx = vmaxq_s16(mx, vbslq_s16(in, vld1q_s16(s.bnd + x), small));
                mn = vminq_s16(mn, vbslq_s16(vandq_u16(in, vtstq_u16(c, c)), vld1q_s16(s.minrow + x), big));
            };
            auto edge = [&](int x) {
                const int16x8_t pos = vaddq_s16(IOs, vdupq_n_s16((int16_t)x));
                block(x, vandq_u16(vcgeq_s16(pos, va), vcltq_s16(pos, vb)));
            };
            const int x0 = a & ~7, x1 = b & ~7;
            edge(x0);
            for (int x = x0 + 8; x < x1; x += 8) block(x, vdupq_n_u16(0xFFFF));
            if (x1 > x0 && x1 < b) edge(x1);
            const int ub = vmaxvq_s16(mx), i0 = vminvq_s16(mn);
            lo = std::min(lo, i0);
            hi = std::max(hi, (dmax - off) + quanta - 1 + std::max(0, ub - MINSC - 2));
            if (nc < NeonScratch::COMP_CAP) s.comps[nc] = Comp{a, b, ub, i0, dmax};
            nc++;
        }
        d = neon_next(s.mw, b, nd, 0);
    }
    if (hi < 0) return B1;
    hb = std::max(0, lo);
    he = std::min(len1 - 1, hi);
    s.ncomp = nc;
    return B2;
}
// Exact repeats of the previous job (a mate rescued twice in a row against the same window) return
// the previous result without recomputing it: the result is a function of the window bytes, the
// query bytes, the lengths and max_hits only.
static inline Kind lean_neon(const Job &jb, NeonScratch &s, int &hb, int &he, int max_hits = 1 << 30)
{
    if (s.memo_ok && jb.len1 == s.memo_len1 && jb.len2 == s.qlen_c && max_hits == s.memo_mh &&
        memcmp(jb.ref, s.rbuf + 4, (size_t)jb.len1) == 0 && memcmp(jb.qry, s.qcache, (size_t)jb.len2) == 0) {
        hb = s.memo_hb; he = s.memo_he;
        s.memo_hits++;
        return s.memo_kind;
    }
    const Kind k = lean_neon_core(jb, s, hb, he, max_hits);
    s.memo_ok = k != FALLBACK;
    s.memo_len1 = jb.len1; s.memo_mh = max_hits; s.memo_hb = hb; s.memo_he = he; s.memo_kind = k;
    return k;
}
#endif  // __aarch64__ (rescue_prune.h dispatches here only on aarch64)

}  // namespace rescue_prune_neon

#endif
