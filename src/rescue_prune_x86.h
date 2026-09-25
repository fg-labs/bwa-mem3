/* x86 (AVX2 floor) implementation of the exact rescue-pruning filter (see rescue_prune.h).
 *
 * rescue_prune_x86::lean_x86() is a port of rescue_prune_neon::lean_neon() (its design comment,
 * steps 1-7, applies unchanged) to 128-bit SSE4.1 / SSSE3, which every AVX2 host has. It keeps the
 * NEON data layout and scan structure, so it returns the same (Kind, hb, he) on every input AND
 * leaves the same per-diagonal view (cnt, minrow, bnd, the mw / hw bitsets and the component
 * list) that rescue_band.cpp consumes. Differences from the NEON code, none of them in any value:
 *  - step 3's 128-byte presence-bitmap lookup is eight 16-byte PSHUFB lookups joined by a
 *    three-level PBLENDVB tree on the table index bits (PSHUFB reads only index bits 0-3 and 7,
 *    and the index is < 128);
 *  - step 6's bitsets are built after the backward scan from the diagonal-major bnd and cnt with
 *    MOVEMASK (mw = {bnd >= 19}, hw = mw and {cnt > 0}), instead of per-step bytes transposed
 *    and scattered; NEON's hw test P(d) >= P(d - 1) is exactly cnt(d) >= 1.
 * The threshold is fixed at 19 (min_seed_len * a at defaults); the caller dispatches here only in
 * that case. */
#ifndef BWA_MEM3_RESCUE_PRUNE_X86_H
#define BWA_MEM3_RESCUE_PRUNE_X86_H

#include "rescue_prune_neon.h"

#if defined(__AVX2__)
#include <immintrin.h>

namespace rescue_prune_x86 {

using rescue_prune_neon::B1;
using rescue_prune_neon::B2;
using rescue_prune_neon::Comp;
using rescue_prune_neon::FAST_HITS;
using rescue_prune_neon::FULL;
using rescue_prune_neon::Job;
using rescue_prune_neon::Kind;
using rescue_prune_neon::MINSC;
using rescue_prune_neon::SCALAR;
using rescue_prune_neon::neon_last;
using rescue_prune_neon::neon_next;

/* rescue_prune_neon::NeonScratch without the NEON-only bitset staging (MH). */
struct X86Scratch {
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
    alignas(16) int16_t bnd[CAP + 64];   // + 64: rescue_band.cpp reads whole 64-diagonal words
    alignas(16) uint64_t mw[CAP / 64 + 3], hw[CAP / 64 + 3];  // bitsets
    alignas(16) uint8_t shuf[256][16];   // left-pack shuffles for 8 x u16 lanes
    static const int COMP_CAP = 128;
    Comp comps[COMP_CAP];
    int ncomp = 0;
    uint8_t pc[256];
    uint8_t qcache[CAP];
    int qlen_c = -1;
    bool q_has_n = false;
    int view_nd = -1;
    bool memo_ok = false;
    int memo_len1 = -1, memo_mh = 0, memo_hb = -1, memo_he = -1, memo_view_nd = -1;
    Kind memo_kind = FULL;
    X86Scratch()
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
    }
};

static inline __m128i x86_ld(const void *p) { return _mm_loadu_si128((const __m128i *)p); }
static inline void x86_st(void *p, __m128i v) { _mm_storeu_si128((__m128i *)p, v); }

// Horizontal max / min of 8 x int16 (lane 0 only ever combines real lanes, so the zero fill of
// the byte shifts never reaches the result).
static inline int x86_hmax16(__m128i v)
{
    v = _mm_max_epi16(v, _mm_srli_si128(v, 8));
    v = _mm_max_epi16(v, _mm_srli_si128(v, 4));
    v = _mm_max_epi16(v, _mm_srli_si128(v, 2));
    return (int16_t)_mm_extract_epi16(v, 0);
}
static inline int x86_hmin16(__m128i v)
{
    v = _mm_min_epi16(v, _mm_srli_si128(v, 8));
    v = _mm_min_epi16(v, _mm_srli_si128(v, 4));
    v = _mm_min_epi16(v, _mm_srli_si128(v, 2));
    return (int16_t)_mm_extract_epi16(v, 0);
}

// In-register 8x8 transpose of 16-bit lanes: x[k][s] <-> x[s][k] (neon_transpose8).
static inline void x86_transpose8(__m128i *x)
{
    const __m128i t0 = _mm_unpacklo_epi16(x[0], x[1]), t1 = _mm_unpackhi_epi16(x[0], x[1]);
    const __m128i t2 = _mm_unpacklo_epi16(x[2], x[3]), t3 = _mm_unpackhi_epi16(x[2], x[3]);
    const __m128i t4 = _mm_unpacklo_epi16(x[4], x[5]), t5 = _mm_unpackhi_epi16(x[4], x[5]);
    const __m128i t6 = _mm_unpacklo_epi16(x[6], x[7]), t7 = _mm_unpackhi_epi16(x[6], x[7]);
    const __m128i u0 = _mm_unpacklo_epi32(t0, t2), u1 = _mm_unpackhi_epi32(t0, t2);
    const __m128i u2 = _mm_unpacklo_epi32(t1, t3), u3 = _mm_unpackhi_epi32(t1, t3);
    const __m128i u4 = _mm_unpacklo_epi32(t4, t6), u5 = _mm_unpackhi_epi32(t4, t6);
    const __m128i u6 = _mm_unpacklo_epi32(t5, t7), u7 = _mm_unpackhi_epi32(t5, t7);
    x[0] = _mm_unpacklo_epi64(u0, u4); x[1] = _mm_unpackhi_epi64(u0, u4);
    x[2] = _mm_unpacklo_epi64(u1, u5); x[3] = _mm_unpackhi_epi64(u1, u5);
    x[4] = _mm_unpacklo_epi64(u2, u6); x[5] = _mm_unpackhi_epi64(u2, u6);
    x[6] = _mm_unpacklo_epi64(u3, u7); x[7] = _mm_unpackhi_epi64(u3, u7);
}

// Low byte of the 5-mer codes of rows b..b+15 (a1 a2 a3 a4, 2 bits each) from buf = padded + b;
// a0 = the code's high 2 bits. Bases are 0..3, so the 16-bit shifts never carry across bytes.
static inline __m128i x86_code_lo(const uint8_t *buf, __m128i &a0)
{
    a0 = x86_ld(buf);
    const __m128i a1 = x86_ld(buf + 1), a2 = x86_ld(buf + 2), a3 = x86_ld(buf + 3), a4 = x86_ld(buf + 4);
    return _mm_or_si128(_mm_or_si128(_mm_slli_epi16(a1, 6), _mm_slli_epi16(a2, 4)),
                        _mm_or_si128(_mm_slli_epi16(a3, 2), a4));
}

// 5-mer codes of rows b..b+15, where buf[k + 4] is sequence position k (buf = padded + b).
static inline void x86_codes16(const uint8_t *buf, uint16_t *out)
{
    __m128i a0;
    const __m128i lo = x86_code_lo(buf, a0);
    x86_st(out, _mm_unpacklo_epi8(lo, a0));
    x86_st(out + 8, _mm_unpackhi_epi8(lo, a0));
}

static inline Kind lean_x86_core(const Job &jb, X86Scratch &s, int &hb, int &he, int max_hits)
{
    const uint8_t *ref = jb.ref, *q = jb.qry;
    const int len1 = jb.len1, len2 = jb.len2;
    hb = he = -1;
    s.view_nd = -1;
    const __m128i kFC = _mm_set1_epi8((char)0xFC);
    // ---- 1. N check (+ copy ref and query into zero-padded buffers) ----
    if (__builtin_expect(len1 + 64 > X86Scratch::CAP || len2 + 64 > X86Scratch::CAP, 0))
        return SCALAR;
    __m128i ov = _mm_setzero_si128();
    uint8_t orv = 0;
    {
        uint8_t *rb = s.rbuf + 4;
        int i = 0;
        for (; i + 16 <= len1; i += 16) { const __m128i x = x86_ld(ref + i); ov = _mm_or_si128(ov, x); x86_st(rb + i, x); }
        for (; i < len1; i++) { orv |= ref[i]; rb[i] = ref[i]; }
        x86_st(rb + len1, _mm_setzero_si128()); x86_st(rb + len1 + 16, _mm_setzero_si128());
    }
    const int quanta = ((len2 + 15) / 16) * 16, off = quanta, nd = len1 + quanta + 1;
    // ---- 2. query table (cached per oriented query), as lean_neon ----
    if (len2 != s.qlen_c || memcmp(q, s.qcache, (size_t)len2) != 0) {
        uint8_t *qb = s.qbuf + 4;
        __m128i qv = _mm_setzero_si128();
        uint8_t qor = 0;
        int j = 0;
        for (; j + 16 <= len2; j += 16) { const __m128i x = x86_ld(q + j); qv = _mm_or_si128(qv, x); x86_st(qb + j, x); }
        for (; j < len2; j++) { qor |= q[j]; qb[j] = q[j]; }
        x86_st(qb + len2, _mm_setzero_si128()); x86_st(qb + len2 + 16, _mm_setzero_si128());
        s.q_has_n = (qor & 0xFC) != 0 || !_mm_testz_si128(qv, kFC);
        memcpy(s.qcache, q, (size_t)len2);
        s.qlen_c = len2;
        if (!s.q_has_n && len2 >= 5) {
            memset(s.tab, 0, sizeof s.tab);
            memset(s.pres, 0, sizeof s.pres);
            for (int b = 0; b < len2; b += 16) x86_codes16(s.qbuf + b, s.qcode + b);
            for (int jj = 4; jj < len2; jj++) {
                const int c = s.qcode[jj];
                const uint32_t old = s.tab[c], occ = (old >> 16) + 1;
                s.nxt[jj] = occ == 1 ? (int16_t)-1 : (int16_t)(off - (old & 0x1FFF));
                s.tab[c] = occ << 16 | (uint32_t)(off - jj) | (occ == 1 ? 0 : 0x2000);
                s.pres[c >> 3] |= (uint8_t)(1u << (c & 7));
            }
        }
    }
    if ((orv & 0xFC) || !_mm_testz_si128(ov, kFC) || s.q_has_n) return FULL;
    if (len1 < 5 || len2 < 5) return B1;
    if (__builtin_expect(nd + 32 > X86Scratch::CAP, 0)) return SCALAR;

    const __m128i IO = _mm_setr_epi16(0, 1, 2, 3, 4, 5, 6, 7), k8 = _mm_set1_epi16(8);

    // ---- 3. hit rows, 16 at a time; left-pack the hit rows and their codes ----
    int np;
    {
        __m128i T[8];
        for (int t = 0; t < 8; t++) T[t] = x86_ld(s.pres + 16 * t);
        const __m128i P2 = _mm_setr_epi8(1, 2, 4, 8, 16, 32, 64, (char)128, 0, 0, 0, 0, 0, 0, 0, 0);
        const __m128i k1f = _mm_set1_epi8(0x1F), k7 = _mm_set1_epi8(7);
        uint16_t *rp = s.PR, *cp = s.PC;
        auto block = [&](int b, unsigned keep) {
            __m128i a0;
            const __m128i lo = x86_code_lo(s.rbuf + b, a0);
            // code >> 3 (0..127): a0 << 5 | lo >> 3 (the 16-bit shift pulls the next byte's low
            // bits into bits 5-7, masked off)
            const __m128i idx = _mm_or_si128(_mm_slli_epi16(a0, 5), _mm_and_si128(_mm_srli_epi16(lo, 3), k1f));
            // pres[idx]: PSHUFB reads idx bits 0-3 (bit 7 is clear); bits 4, 5, 6 pick the table
            // through PBLENDVB, which tests byte bit 7 (the 16-bit shifts put bit n there)
            const __m128i s4 = _mm_slli_epi16(idx, 3), s5 = _mm_slli_epi16(idx, 2), s6 = _mm_slli_epi16(idx, 1);
            const __m128i l01 = _mm_blendv_epi8(_mm_shuffle_epi8(T[0], idx), _mm_shuffle_epi8(T[1], idx), s4);
            const __m128i l23 = _mm_blendv_epi8(_mm_shuffle_epi8(T[2], idx), _mm_shuffle_epi8(T[3], idx), s4);
            const __m128i l45 = _mm_blendv_epi8(_mm_shuffle_epi8(T[4], idx), _mm_shuffle_epi8(T[5], idx), s4);
            const __m128i l67 = _mm_blendv_epi8(_mm_shuffle_epi8(T[6], idx), _mm_shuffle_epi8(T[7], idx), s4);
            const __m128i byte = _mm_blendv_epi8(_mm_blendv_epi8(l01, l23, s5), _mm_blendv_epi8(l45, l67, s5), s6);
            const __m128i bit = _mm_shuffle_epi8(P2, _mm_and_si128(lo, k7));
            const __m128i hit = _mm_cmpeq_epi8(_mm_and_si128(byte, bit), bit);
            const unsigned m16 = (unsigned)_mm_movemask_epi8(hit) & keep;
            const __m128i c0 = _mm_unpacklo_epi8(lo, a0), c1 = _mm_unpackhi_epi8(lo, a0);
            const __m128i pos = _mm_add_epi16(IO, _mm_set1_epi16((short)b));
            const unsigned ml = m16 & 0xFF, mh = m16 >> 8;
            const __m128i shl = x86_ld(s.shuf[ml]), shh = x86_ld(s.shuf[mh]);
            x86_st(rp, _mm_shuffle_epi8(pos, shl));
            x86_st(cp, _mm_shuffle_epi8(c0, shl));
            rp += s.pc[ml]; cp += s.pc[ml];
            x86_st(rp, _mm_shuffle_epi8(_mm_add_epi16(pos, k8), shh));
            x86_st(cp, _mm_shuffle_epi8(c1, shh));
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
    // ---- 4+5. accumulate hits (layer 1, then layers >= 2), as lean_neon ----
    int hits = np;
    const __m128i koff = _mm_set1_epi16((short)off), one16 = _mm_set1_epi16(1);
    auto pack_runs = [&](uint16_t *&bp, uint16_t *&dp, uint16_t *&rp, __m128i idx, __m128i d, __m128i r, unsigned m) {
        const __m128i sh = x86_ld(s.shuf[m]);
        x86_st(bp, _mm_shuffle_epi8(idx, sh));
        x86_st(dp, _mm_shuffle_epi8(d, sh));
        x86_st(rp, _mm_shuffle_epi8(r, sh));
        bp += s.pc[m]; dp += s.pc[m]; rp += s.pc[m];
    };
    // (valid & ~cont) lanes -> bits 0-7, sel lanes -> bits 8-15 (both are 0 / 0xFFFF per lane)
    auto mask2 = [](__m128i lo8, __m128i hi8) { return (unsigned)_mm_movemask_epi8(_mm_packs_epi16(lo8, hi8)); };
    int n;   // rows in the current layer >= 2 (R = s.RA, J = s.JA)
    {
        uint16_t *bp = s.B, *bdp = s.BD, *brp = s.BR, *mp = s.RA, *jp = (uint16_t *)s.JA;
        __m128i pR = _mm_setzero_si128(), pD = _mm_set1_epi16(-1), idx = IO;
        const __m128i vn = _mm_set1_epi16((short)np), k1fff = _mm_set1_epi16(0x1FFF), k2000 = _mm_set1_epi16(0x2000);
        const __m128i lo16 = _mm_set1_epi32(0xFFFF);
        const uint32_t *__restrict tab = s.tab;
        // table entries of four codes (lanes past np hold stale codes: masked below)
        auto look4 = [tab](const uint16_t *c) {
            uint64_t x;
            memcpy(&x, c, 8);
            const uint64_t a = tab[x & 1023] | (uint64_t)tab[(x >> 16) & 1023] << 32;
            const uint64_t b = tab[(x >> 32) & 1023] | (uint64_t)tab[(x >> 48) & 1023] << 32;
            return _mm_set_epi64x((long long)b, (long long)a);
        };
        __m128i extra = _mm_setzero_si128();   // gate: occurrences - 1 summed over the hit rows (u32)
        for (int k0 = 0; k0 < np; k0 += 8, idx = _mm_add_epi16(idx, k8)) {
            const __m128i r = x86_ld(s.PR + k0);
            const __m128i e0 = look4(s.PC + k0), e1 = look4(s.PC + k0 + 4);
            const __m128i t = _mm_packus_epi32(_mm_and_si128(e0, lo16), _mm_and_si128(e1, lo16));
            const __m128i occ = _mm_packus_epi32(_mm_srli_epi32(e0, 16), _mm_srli_epi32(e1, 16));
            const __m128i tv = _mm_and_si128(t, k1fff);                     // off - j_last
            const __m128i d = _mm_add_epi16(r, tv);
            const __m128i rp = _mm_alignr_epi8(r, pR, 14), dp = _mm_alignr_epi8(d, pD, 14);
            pR = r; pD = d;
            const __m128i valid = _mm_cmplt_epi16(idx, vn);
            const __m128i cont = _mm_and_si128(_mm_cmpeq_epi16(r, _mm_add_epi16(rp, one16)), _mm_cmpeq_epi16(d, dp));
            const __m128i mu = _mm_and_si128(_mm_cmpeq_epi16(_mm_and_si128(t, k2000), k2000), valid);
            extra = _mm_add_epi32(extra, _mm_madd_epi16(_mm_and_si128(_mm_sub_epi16(occ, one16), valid), one16));
            const unsigned m = mask2(_mm_andnot_si128(cont, valid), mu);
            pack_runs(bp, bdp, brp, idx, d, r, m & 0xFF);
            const __m128i sh = x86_ld(s.shuf[m >> 8]);
            x86_st(mp, _mm_shuffle_epi8(r, sh));
            x86_st(jp, _mm_shuffle_epi8(_mm_sub_epi16(koff, tv), sh));
            mp += s.pc[m >> 8]; jp += s.pc[m >> 8];
        }
        n = (int)(mp - s.RA);
        extra = _mm_add_epi32(extra, _mm_srli_si128(extra, 8));
        extra = _mm_add_epi32(extra, _mm_srli_si128(extra, 4));
        if ((long)np + (uint32_t)_mm_cvtsi128_si32(extra) > max_hits) return FULL;
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
        const __m128i m1 = _mm_set1_epi16(-1);
        while (n > 0) {
            hits += n;
            if (__builtin_expect(hits > 32000, 0)) return SCALAR;
            uint16_t *bp = s.B, *bdp = s.BD, *brp = s.BR, *rp2 = R2, *jp2 = (uint16_t *)J2;
            __m128i pR = _mm_setzero_si128(), pD = _mm_set1_epi16(-1), idx = IO;
            const __m128i vn = _mm_set1_epi16((short)n);
            for (int k0 = 0; k0 < n; k0 += 8, idx = _mm_add_epi16(idx, k8)) {
                const __m128i r = x86_ld(R + k0);
                const __m128i d = _mm_sub_epi16(_mm_add_epi16(r, koff), x86_ld(J + k0));
                const __m128i rp = _mm_alignr_epi8(r, pR, 14), dp = _mm_alignr_epi8(d, pD, 14);
                pR = r; pD = d;
                const __m128i valid = _mm_cmplt_epi16(idx, vn);
                const __m128i cont = _mm_and_si128(_mm_cmpeq_epi16(r, _mm_add_epi16(rp, one16)), _mm_cmpeq_epi16(d, dp));
                const __m128i jn = _mm_set_epi64x((long long)nxt4(J + k0 + 4), (long long)nxt4(J + k0));
                const __m128i keep = _mm_and_si128(_mm_cmpgt_epi16(jn, m1), valid);
                const unsigned m = mask2(_mm_andnot_si128(cont, valid), keep);
                pack_runs(bp, bdp, brp, idx, d, r, m & 0xFF);
                const __m128i sh = x86_ld(s.shuf[m >> 8]);
                x86_st(rp2, _mm_shuffle_epi8(r, sh));
                x86_st(jp2, _mm_shuffle_epi8(jn, sh));
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

    // ---- 6. Kadane bounds as segment-transposed scans, as lean_neon ----
    const int L = ((nd + 63) >> 6) << 3, NT = 8 * L;
    const __m128i big = _mm_set1_epi16(32767), small = _mm_set1_epi16(-32768), zero = _mm_setzero_si128();
    __m128i CP, Mb, Xa;   // per segment: P before it; min(0, P before it); max P after it
    {
        __m128i p = zero, rm = big, mx = small, bestl = small;
        int16_t *__restrict TP = s.P + 8, *__restrict TM = s.PM;
        for (int t = 0; t < L; t += 8) {
            __m128i c[8];
            for (int k = 0; k < 8; k++) c[k] = x86_ld(s.cnt + k * L + t);
            x86_transpose8(c);
            for (int k = 0; k < 8; k++) {
                p = _mm_add_epi16(p, _mm_sub_epi16(c[k], one16));
                x86_st(TP + (t + k) * 8, p);
                x86_st(TM + (t + k) * 8, rm);            // segment-local exclusive running min
                bestl = _mm_max_epi16(bestl, _mm_subs_epi16(p, rm));
                rm = _mm_min_epi16(rm, p);
                mx = _mm_max_epi16(mx, p);
            }
        }
        // carries across lanes (segments): exclusive prefix sum of the totals; exclusive prefix
        // min (with 0) of the segment minima; exclusive suffix max of the segment maxima
        __m128i x = p;
        x = _mm_add_epi16(x, _mm_slli_si128(x, 2));
        x = _mm_add_epi16(x, _mm_slli_si128(x, 4));
        x = _mm_add_epi16(x, _mm_slli_si128(x, 8));
        CP = _mm_slli_si128(x, 2);
        __m128i g = _mm_add_epi16(CP, rm);
        g = _mm_min_epi16(g, _mm_alignr_epi8(g, big, 14));
        g = _mm_min_epi16(g, _mm_alignr_epi8(g, big, 12));
        g = _mm_min_epi16(g, _mm_alignr_epi8(g, big, 8));
        Mb = _mm_min_epi16(_mm_slli_si128(g, 2), zero);         // lane 0: min(0, nothing) = 0
        const __m128i h = _mm_add_epi16(CP, mx);
        __m128i hs = _mm_max_epi16(h, _mm_alignr_epi8(small, h, 2));
        hs = _mm_max_epi16(hs, _mm_alignr_epi8(small, hs, 4));
        hs = _mm_max_epi16(hs, _mm_alignr_epi8(small, hs, 8));
        Xa = _mm_alignr_epi8(small, hs, 2);
        // max fwd over segment s = max(local best, max P in s - Mb)
        const __m128i best = _mm_max_epi16(bestl, _mm_sub_epi16(h, Mb));
        if (5 + x86_hmax16(best) < MINSC) return B1;
    }
    // Backward: SX, PM and bnd per segment step (see lean_neon for the FAST_HITS form), bnd
    // transposed back to diagonal order.
    {
        const __m128i five = _mm_set1_epi16(5);
        const int16_t *__restrict TP = s.P + 8, *__restrict TM = s.PM;
        const bool fast = hits <= FAST_HITS;
        const __m128i Xl = _mm_subs_epi16(Xa, CP), Ml = _mm_sub_epi16(Mb, CP);
        __m128i sx = fast ? Xl : small;
        for (int t = L - 8; t >= 0; t -= 8) {
            __m128i bd[8];
            for (int k = 7; k >= 0; k--) {
                const int o = (t + k) * 8;
                sx = _mm_max_epi16(sx, x86_ld(TP + o));
                if (fast) {
                    bd[k] = _mm_add_epi16(five, _mm_sub_epi16(sx, _mm_min_epi16(Ml, x86_ld(TM + o))));
                } else {
                    const __m128i SX = _mm_max_epi16(_mm_add_epi16(CP, sx), Xa);
                    const __m128i PM = _mm_min_epi16(Mb, _mm_adds_epi16(CP, x86_ld(TM + o)));
                    bd[k] = _mm_add_epi16(five, _mm_sub_epi16(SX, PM));
                }
            }
            x86_transpose8(bd);
            for (int k = 0; k < 8; k++) x86_st(s.bnd + k * L + t, bd[k]);
        }
        // Bitsets from the diagonal-major arrays, 16 diagonals per step: mw = {bnd >= 19},
        // hw = mw and {cnt > 0} (cnt <= 32000 here, so the signed compare is exact).
        const __m128i t18 = _mm_set1_epi16(MINSC - 1);
        uint8_t *mb = (uint8_t *)s.mw, *hbb = (uint8_t *)s.hw;
        for (int x0 = 0; x0 < NT; x0 += 16) {
            const unsigned m = (unsigned)_mm_movemask_epi8(_mm_packs_epi16(
                _mm_cmpgt_epi16(x86_ld(s.bnd + x0), t18), _mm_cmpgt_epi16(x86_ld(s.bnd + x0 + 8), t18)));
            const unsigned c = (unsigned)_mm_movemask_epi8(_mm_packs_epi16(
                _mm_cmpgt_epi16(x86_ld(s.cnt + x0), zero), _mm_cmpgt_epi16(x86_ld(s.cnt + x0 + 8), zero)));
            const uint16_t mm = (uint16_t)m, hh = (uint16_t)(m & c);
            memcpy(mb + (x0 >> 3), &mm, 2);
            memcpy(hbb + (x0 >> 3), &hh, 2);
        }
        // drop the pad diagonals [nd, NT) from mw (hw has none: cnt = 0 there), as lean_neon
        if (nd & 7) mb[nd >> 3] &= (uint8_t)((1u << (nd & 7)) - 1);
        memset(mb + ((nd + 7) >> 3), 0, (NT >> 3) - ((nd + 7) >> 3) + 16);
        memset(hbb + (NT >> 3), 0, 16);
    }

    // ---- 7. components (runs of mw) with a hit, as lean_neon ----
    int hi = -1, lo = 32767, nc = 0;
    for (int d = neon_next(s.mw, 0, nd, 0); d < nd;) {
        const int a = d, b = neon_next(s.mw, a, nd, ~0ull);
        const int dmax = neon_last(s.hw, a, b);
        if (dmax >= 0) {
            const __m128i va1 = _mm_set1_epi16((short)(a - 1)), vb = _mm_set1_epi16((short)b);
            __m128i mx = small, mn = big;
            // one block of 8 diagonals; in = lanes inside [a, b) (all of them for inner blocks)
            auto block = [&](int x, __m128i in) {
                const __m128i c = x86_ld(s.cnt + x);
                mx = _mm_max_epi16(mx, _mm_blendv_epi8(small, x86_ld(s.bnd + x), in));
                mn = _mm_min_epi16(mn, _mm_blendv_epi8(big, x86_ld(s.minrow + x), _mm_and_si128(in, _mm_cmpgt_epi16(c, zero))));
            };
            auto edge = [&](int x) {
                const __m128i pos = _mm_add_epi16(IO, _mm_set1_epi16((short)x));
                block(x, _mm_and_si128(_mm_cmpgt_epi16(pos, va1), _mm_cmpgt_epi16(vb, pos)));
            };
            const int x0 = a & ~7, x1 = b & ~7;
            edge(x0);
            for (int x = x0 + 8; x < x1; x += 8) block(x, _mm_set1_epi16(-1));
            if (x1 > x0 && x1 < b) edge(x1);
            const int ub = x86_hmax16(mx), i0 = x86_hmin16(mn);
            lo = std::min(lo, i0);
            hi = std::max(hi, (dmax - off) + quanta - 1 + std::max(0, ub - 21));
            if (nc < X86Scratch::COMP_CAP) s.comps[nc] = Comp{a, b, ub, i0, dmax};
            nc++;
        }
        d = neon_next(s.mw, b, nd, 0);
    }
    if (hi < 0) return B1;
    hb = std::max(0, lo);
    he = std::min(len1 - 1, hi);
    s.ncomp = nc;
    s.view_nd = nd;
    return B2;
}

// Exact repeats of the previous job return the previous result, as lean_neon.
static inline Kind lean_x86(const Job &jb, X86Scratch &s, int &hb, int &he, int max_hits = 1 << 30)
{
    if (s.memo_ok && jb.len1 == s.memo_len1 && jb.len2 == s.qlen_c && max_hits == s.memo_mh &&
        memcmp(jb.ref, s.rbuf + 4, (size_t)jb.len1) == 0 && memcmp(jb.qry, s.qcache, (size_t)jb.len2) == 0) {
        hb = s.memo_hb; he = s.memo_he; s.view_nd = s.memo_view_nd;
        return s.memo_kind;
    }
    const Kind k = lean_x86_core(jb, s, hb, he, max_hits);
    s.memo_ok = k != SCALAR;
    s.memo_len1 = jb.len1; s.memo_mh = max_hits; s.memo_hb = hb; s.memo_he = he;
    s.memo_view_nd = s.view_nd; s.memo_kind = k;
    return k;
}

}  // namespace rescue_prune_x86

#endif  // __AVX2__
#endif
