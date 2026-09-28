/* AVX2 32-lane banded rescue kernel: the x86 counterpart of the NEON kernel section of
 * rescue_band.cpp (included from there, inside its anonymous namespace context). Same band
 * coordinates, SoA layout (one byte per lane, RB_L = 32 lanes per position / row), kernels
 * (rb_dp_core<FScan>, rb_dp_wave2), qe scan and row bookkeeping; every value the caller reads
 * (gmax, te, qe, the row maxima R) is the NEON kernel's for the same lanes. Two encodings differ,
 * because x86 has no NEON TBL / SQADD:
 *
 *  - Query pad code 0x80 (RB_QPAD) instead of 0x40. PSHUFB returns 0 only when bit 7 of the index
 *    is set (and otherwise reads its low nibble), so the pad and the nonexistent code 0xC0 must both
 *    have bit 7 set to XOR any real reference code (0-4) to a zero-scoring index, as TBL does for
 *    any index >= 16 on NEON. An inactive row (reference code 0x80, past the lane's range) against
 *    a real column XORs to a bit-7 index too; against a pad / nonexistent column to a bit-7-clear
 *    one. Either way it scores something, and it does not matter: inactive rows are excluded from
 *    imp (so from gmax, te, the flushes and the qe scan), their row maxima are zeroed when R goes
 *    back to lane-major, and their H / E only feed later rows of the same lane, which are inactive
 *    too (a lane's rows are the prefix [0, nrows)).
 *  - The cell's H + score is the biased form of the x86 kswv u8 kernels:
 *        h = subs_epu8(adds_epu8(Hp, sc), bias),   sc = PSHUFB(tblb, q ^ r),
 *    tblb = the score table + shift (match a + shift, mismatch shift - b, N shift - 1; shift =
 *    max(1, b), all >= 0), bias = shift on real query columns and 0 on pad / nonexistent ones
 *    (where sc is 0, so h = Hp as NEON's +0 gives). This equals NEON's sat(Hp + s) iff
 *    Hp + a + shift <= 255 (an H above that plus a match would saturate at 255 before the
 *    subtract). Every band H is a local alignment score, so H <= ub1 in pass 0 (the K-mer bound of
 *    the parent's best component; plan() refuses ub1 + a + shift > 255 on x86) and H <= S in pass 1
 *    (take_pass1 requires S + a + shift <= 255). At the default scoring: H <= 250.
 * The live-column mask QL (A < 0xC0) is unchanged; the per-position bias BI is built from A at the
 * start of each kernel. */
#ifndef BWA_MEM3_RESCUE_BAND_KERNEL_X86_H
#define BWA_MEM3_RESCUE_BAND_KERNEL_X86_H

#include <immintrin.h>
#include <algorithm>
#include <cstring>
#include <vector>

#define RB_HAVE_SIMD 1
#define RB_X86 1
static const int RB_L = 32;
static const int RB_PAD = 16;
static const uint8_t RB_QPAD = 0x80;

namespace {

struct rb_work {
    std::vector<uint8_t> A, QL, BI, REF, H, E, R, SNAP, ST, RL;
    alignas(32) uint16_t te[32];
    alignas(32) uint8_t gmax[32];
    /* Per-lane early-exit target (pass 1: the pass-0 score S; 0 for empty lanes). */
    alignas(32) uint8_t target[32];
    /* rb_dp_wave2's qe: band index of the lane's qe cell, k = kb_chunk * 255 + kb_off - 1. */
    alignas(32) uint8_t kb_chunk[32], kb_off[32];
    /* The batch scoring (rb_set_scoring): the biased table (both 128-bit halves), the bias, and the
     * gap constants (see the NEON rb_work). */
    alignas(32) uint8_t tblb[32];
    uint8_t shift, oe_del, e_del, oe_ins, e_ins;
    bool asym;   // tblb is indexed (r << 2) | q (REF holds r << 2), not q ^ r
    template <class V> static void fit(V &v, size_t n) { if (v.size() < n) v.resize(n); }
};

static inline void rb_set_scoring(rb_work &w, const rb_scoring &sc)
{
    const int sh = sc.shift();
    w.asym = sc.asym;
    for (int h = 0; h < 32; h += 16) {
        if (sc.asym) {
            for (int i = 0; i < 16; i++) w.tblb[h + i] = (uint8_t)(sc.mat16[i] + sh);
            continue;
        }
        w.tblb[h] = (uint8_t)(sc.a + sh);
        for (int i = 1; i < 4; i++) w.tblb[h + i] = (uint8_t)(sh - sc.b);
        for (int i = 4; i < 16; i++) w.tblb[h + i] = (uint8_t)(sh - 1);
    }
    w.shift = (uint8_t)sh;
    w.oe_del = (uint8_t)std::min(255, sc.o_del + sc.e_del); w.e_del = (uint8_t)std::min(255, sc.e_del);
    w.oe_ins = (uint8_t)std::min(255, sc.o_ins + sc.e_ins); w.e_ins = (uint8_t)std::min(255, sc.e_ins);
}

static inline __m256i rb_ld(const uint8_t *p) { return _mm256_loadu_si256((const __m256i *)p); }
static inline void rb_st(uint8_t *p, __m256i v) { _mm256_storeu_si256((__m256i *)p, v); }
static inline uint32_t rb_mask32(__m256i m) { return (uint32_t)_mm256_movemask_epi8(m); }
/* unsigned byte compares: a >= b, a > b */
static inline __m256i rb_ge_u8(__m256i a, __m256i b) { return _mm256_cmpeq_epi8(_mm256_max_epu8(a, b), a); }
static inline __m256i rb_gt_u8(__m256i a, __m256i b) { return _mm256_andnot_si256(_mm256_cmpeq_epi8(a, b), rb_ge_u8(a, b)); }
/* every lane's gmax has reached its target */
static inline bool rb_all_ge(__m256i gmax, __m256i tgt) { return rb_mask32(rb_ge_u8(gmax, tgt)) == 0xFFFFFFFFu; }
/* te (u16 per lane, two halves of 16 lanes) = row where imp is set */
static inline void rb_te_update(__m256i &te_lo, __m256i &te_hi, __m256i imp, int row)
{
    const __m256i rv = _mm256_set1_epi16((short)row);
    te_lo = _mm256_blendv_epi8(te_lo, rv, _mm256_cvtepi8_epi16(_mm256_castsi256_si128(imp)));
    te_hi = _mm256_blendv_epi8(te_hi, rv, _mm256_cvtepi8_epi16(_mm256_extracti128_si256(imp, 1)));
}
static inline void rb_te_store(rb_work &w, __m256i te_lo, __m256i te_hi)
{
    _mm256_store_si256((__m256i *)w.te, te_lo);
    _mm256_store_si256((__m256i *)(w.te + 16), te_hi);
}

/* In-register 16x16 byte transpose (SSE unpack ladder): r[l] byte t -> r[t] byte l. */
static inline void rb_transpose16_sse(__m128i r[16])
{
    __m128i t[16], u[16], v[16];
    for (int i = 0; i < 8; i++) {
        t[i] = _mm_unpacklo_epi8(r[2 * i], r[2 * i + 1]);
        t[i + 8] = _mm_unpackhi_epi8(r[2 * i], r[2 * i + 1]);
    }
    for (int h = 0; h < 16; h += 8)
        for (int i = 0; i < 4; i++) {
            u[h + i] = _mm_unpacklo_epi16(t[h + 2 * i], t[h + 2 * i + 1]);
            u[h + 4 + i] = _mm_unpackhi_epi16(t[h + 2 * i], t[h + 2 * i + 1]);
        }
    for (int q = 0; q < 16; q += 4)
        for (int i = 0; i < 2; i++) {
            v[q + i] = _mm_unpacklo_epi32(u[q + 2 * i], u[q + 2 * i + 1]);
            v[q + 2 + i] = _mm_unpackhi_epi32(u[q + 2 * i], u[q + 2 * i + 1]);
        }
    for (int q = 0; q < 16; q += 2) {
        r[q] = _mm_unpacklo_epi64(v[q], v[q + 1]);
        r[q + 1] = _mm_unpackhi_epi64(v[q], v[q + 1]);
    }
}

/* ST holds RB_L lane-major rows of stride `stride` (a multiple of 16); write n * RB_L interleaved
 * bytes (position-major, lane-minor) to dst: two 16-lane transposes per 16 positions. */
static inline void rb_transpose_to(const uint8_t *ST, int stride, int n, uint8_t *dst)
{
    __m128i a[16], b[16];
    for (int p = 0; p < n; p += 16) {
        for (int l = 0; l < 16; l++) {
            a[l] = _mm_loadu_si128((const __m128i *)(ST + (size_t)l * stride + p));
            b[l] = _mm_loadu_si128((const __m128i *)(ST + (size_t)(l + 16) * stride + p));
        }
        rb_transpose16_sse(a);
        rb_transpose16_sse(b);
        const int m = std::min(16, n - p);
        for (int t = 0; t < m; t++) {
            _mm_storeu_si128((__m128i *)(dst + (size_t)(p + t) * 32), a[t]);
            _mm_storeu_si128((__m128i *)(dst + (size_t)(p + t) * 32 + 16), b[t]);
        }
    }
}

/* BI[p] = shift on real query columns (A < 0x80), 0 on pad / nonexistent ones, for p in [0, P). */
static inline void rb_build_bias(rb_work &w, int P)
{
    rb_work::fit(w.BI, (size_t)P * 32);
    const __m256i four = _mm256_set1_epi8((char)w.shift), m1 = _mm256_set1_epi8(-1);
    const uint8_t *A = w.A.data();
    uint8_t *BI = w.BI.data();
    for (int p = 0; p < P; p++) rb_st(BI + p * 32, _mm256_and_si256(_mm256_cmpgt_epi8(rb_ld(A + p * 32), m1), four));
}

/* Copy row `row`'s H into SNAP for the lanes in msk (live cells only; dead cells 0). */
static inline void rb_snapshot(rb_work &w, const uint8_t *Hrow, int row, int W, __m256i msk)
{
    uint8_t *SNAP = w.SNAP.data();
    const uint8_t *qp = w.QL.data() + (size_t)row * 32;
    for (int k = W - 1; k >= 0; k--, qp += 32) {
        const __m256i m = _mm256_and_si256(msk, rb_ld(qp));
        const __m256i keep = _mm256_andnot_si256(msk, rb_ld(SNAP + k * 32));   // vbicq(SNAP, msk)
        rb_st(SNAP + k * 32, _mm256_blendv_epi8(keep, rb_ld(Hrow + k * 32), m));
    }
}

/* The biased scoring table (index q ^ r: 0 match, 1-3 mismatch, 4-15 N), both 128-bit halves. */
static inline __m256i rb_tblb(const rb_work &w) { return _mm256_load_si256((const __m256i *)w.tblb); }
/* sat(Hp + score(q, r)) in the biased form (see the file comment). Asym: the --meth table indexed
 * r << 2 | q (REF pre-shifted; the pad 0x80 / nonexistent 0xC0 query codes and the inactive-row
 * code 0x80 still set bit 7, so PSHUFB still gives 0). */
template <bool Asym>
static inline __m256i rb_hs(__m256i hp, __m256i q, __m256i rref, __m256i bi, __m256i tbl)
{
    const __m256i idx = Asym ? _mm256_or_si256(q, rref) : _mm256_xor_si256(q, rref);
    return _mm256_subs_epu8(_mm256_adds_epu8(hp, _mm256_shuffle_epi8(tbl, idx)), bi);
}

/* rb_dp_core of the NEON section, 32 lanes (see there for the row range, the live mask and the
 * FScan cell). */
template <bool FScan, bool Sym, bool Asym>
static long rb_dp_core(rb_work &w, int W, int NR, int omax, int ominq, int omaskq, bool early)
{
    long computed = 0;
    rb_build_bias(w, NR + W - 1);
    const __m256i tbl = rb_tblb(w);
    const __m256i vOI = _mm256_set1_epi8((char)w.oe_ins), vEI = _mm256_set1_epi8((char)w.e_ins);
    const __m256i vEI2 = _mm256_set1_epi8((char)(uint8_t)std::min(255, 2 * w.e_ins));
    const __m256i vOD = _mm256_set1_epi8((char)w.oe_del), vED = _mm256_set1_epi8((char)w.e_del), m1 = _mm256_set1_epi8(-1);
    uint8_t *Hc = w.H.data(), *Hp = Hc + (size_t)W * 32, *E = w.E.data();
    const uint8_t *A = w.A.data(), *QL = w.QL.data(), *BI = w.BI.data(), *REF = w.REF.data();
    uint8_t *Rout = w.R.data();
    memset(Hc, 0, (size_t)W * 64);
    memset(E, 0, (size_t)(W + 1) * 32);
    __m256i gmax = _mm256_setzero_si256(), pend = _mm256_setzero_si256();
    const __m256i tgt = rb_ld(w.target);
    int rlast = NR - 1;
    __m256i te_lo = _mm256_setzero_si256(), te_hi = _mm256_setzero_si256();
    for (int r = 0; r < NR; r++) {
        const __m256i rref = rb_ld(REF + (size_t)r * 32);
        __m256i f = _mm256_setzero_si256(), rmax = _mm256_setzero_si256();
        const int khi = std::min(W - 1, r + omax), klo = std::max(0, r + ominq);
        computed += khi >= klo ? khi - klo + 1 : 0;
        const size_t p0 = (size_t)(r + (W - 1 - khi)) * 32;                    // p = r + W-1-k
        const uint8_t *ap = A + p0, *qp = QL + p0, *bp = BI + p0;
        int k = khi;
        const int kun = std::max(klo, r + omaskq);
#define RB_CELL2(MASK)                                                                          \
        {                                                                                      \
            const __m256i ea = rb_ld(E + k * 32), eb = rb_ld(E + (k - 1) * 32);                \
            const __m256i h0a = _mm256_max_epu8(rb_hs<Asym>(rb_ld(Hp + k * 32), rb_ld(ap), rref, rb_ld(bp), tbl), ea); \
            const __m256i h0b = _mm256_max_epu8(rb_hs<Asym>(rb_ld(Hp + (k - 1) * 32), rb_ld(ap + 32), rref, rb_ld(bp + 32), tbl), eb); \
            const __m256i h07a = _mm256_subs_epu8(h0a, vOI), h07b = _mm256_subs_epu8(h0b, vOI); \
            const __m256i fb = _mm256_max_epu8(h07a, _mm256_subs_epu8(f, vEI));  /* f_in(k-1) */ \
            const __m256i ha = _mm256_max_epu8(h0a, f), hb = _mm256_max_epu8(h0b, fb);         \
            f = _mm256_max_epu8(_mm256_max_epu8(h07b, _mm256_subs_epu8(h07a, vEI)), _mm256_subs_epu8(f, vEI2)); \
            rb_st(Hc + k * 32, ha);                                                            \
            rb_st(Hc + (k - 1) * 32, hb);                                                      \
            const __m256i ra = FScan ? h0a : ha, rb_ = FScan ? h0b : hb;                       \
            if (MASK) rmax = _mm256_max_epu8(rmax, _mm256_max_epu8(_mm256_and_si256(ra, rb_ld(qp)), _mm256_and_si256(rb_, rb_ld(qp + 32)))); \
            else rmax = _mm256_max_epu8(rmax, _mm256_max_epu8(ra, rb_));                       \
            const __m256i hda = FScan ? (Sym ? h07a : _mm256_subs_epu8(h0a, vOD)) : _mm256_subs_epu8(ha, vOD); \
            const __m256i hdb = FScan ? (Sym ? h07b : _mm256_subs_epu8(h0b, vOD)) : _mm256_subs_epu8(hb, vOD); \
            rb_st(E + (k + 1) * 32, _mm256_max_epu8(hda, _mm256_subs_epu8(ea, vED)));          \
            rb_st(E + k * 32, _mm256_max_epu8(hdb, _mm256_subs_epu8(eb, vED)));                \
        }
#define RB_CELL1(MASK)                                                                          \
        {                                                                                      \
            const __m256i e = rb_ld(E + k * 32);                                               \
            const __m256i h0 = _mm256_max_epu8(rb_hs<Asym>(rb_ld(Hp + k * 32), rb_ld(ap), rref, rb_ld(bp), tbl), e); \
            const __m256i h = _mm256_max_epu8(h0, f);                                          \
            rb_st(Hc + k * 32, h);                                                             \
            const __m256i r_ = FScan ? h0 : h;                                                 \
            rmax = _mm256_max_epu8(rmax, (MASK) ? _mm256_and_si256(r_, rb_ld(qp)) : r_);       \
            const __m256i h07 = _mm256_subs_epu8(h0, vOI);                                     \
            const __m256i hd = FScan ? (Sym ? h07 : _mm256_subs_epu8(h0, vOD)) : _mm256_subs_epu8(h, vOD); \
            rb_st(E + (k + 1) * 32, _mm256_max_epu8(hd, _mm256_subs_epu8(e, vED)));            \
            f = _mm256_max_epu8(h07, _mm256_subs_epu8(f, vEI));                                \
        }
        for (; k >= kun + 1; k -= 2, ap += 64, qp += 64, bp += 64) RB_CELL2(false)
        if (k == kun) { RB_CELL1(false) k--; ap += 32; qp += 32; bp += 32; }
        for (; k >= klo + 1; k -= 2, ap += 64, qp += 64, bp += 64) RB_CELL2(true)
        for (; k >= klo; k--, ap += 32, qp += 32, bp += 32) RB_CELL1(true)
#undef RB_CELL2
#undef RB_CELL1
        rb_st(Rout + (size_t)r * 32, rmax);
        const __m256i imp = _mm256_and_si256(rb_gt_u8(rmax, gmax), _mm256_cmpgt_epi8(rref, m1));
        const __m256i flush = _mm256_andnot_si256(imp, pend);
        if (rb_mask32(flush)) rb_snapshot(w, Hp, r - 1, W, flush);
        gmax = _mm256_blendv_epi8(gmax, rmax, imp);
        rb_te_update(te_lo, te_hi, imp, r);
        pend = imp;
        std::swap(Hc, Hp);
        rlast = r;
        if (early && rb_all_ge(gmax, tgt)) break;
    }
    if (rb_mask32(pend)) rb_snapshot(w, Hp, rlast, W, pend);
    rb_te_store(w, te_lo, te_hi);
    rb_st(w.gmax, gmax);
    return computed;
}

/* rb_qe_scan of the NEON section, 32 lanes. */
static inline void rb_qe_scan(rb_work &w, const uint8_t *Hrow, int row, int W, int omax, int ominq,
                              __m256i msk, __m256i gmax)
{
    const int khi = std::min(W - 1, row + omax), klo = std::max(0, row + ominq);
    const __m256i v1 = _mm256_set1_epi8(1), v2 = _mm256_set1_epi8(2), zero = _mm256_setzero_si256();
    __m256i chunk = rb_ld(w.kb_chunk), off = rb_ld(w.kb_off);
    const uint8_t *QL = w.QL.data();
    for (int k0 = klo, c = klo / 255; k0 <= khi; k0 = (c + 1) * 255, c++) {
        const int k1 = std::min(khi, c * 255 + 254);
        __m256i m0 = zero, m1 = zero;
        __m256i kv = _mm256_set1_epi8((char)(uint8_t)(k0 - c * 255 + 1));
        const uint8_t *hp = Hrow + (size_t)k0 * 32, *qp = QL + (size_t)(row + W - 1 - k0) * 32;   // p = row + W-1-k
        int k = k0;
        for (; k + 1 <= k1; k += 2, hp += 64, qp -= 64, kv = _mm256_add_epi8(kv, v2)) {
            const __m256i e0 = _mm256_and_si256(_mm256_cmpeq_epi8(rb_ld(hp), gmax), rb_ld(qp));
            const __m256i e1 = _mm256_and_si256(_mm256_cmpeq_epi8(rb_ld(hp + 32), gmax), rb_ld(qp - 32));
            m0 = _mm256_max_epu8(m0, _mm256_and_si256(e0, kv));
            m1 = _mm256_max_epu8(m1, _mm256_and_si256(e1, _mm256_add_epi8(kv, v1)));
        }
        if (k <= k1) m0 = _mm256_max_epu8(m0, _mm256_and_si256(_mm256_and_si256(_mm256_cmpeq_epi8(rb_ld(hp), gmax), rb_ld(qp)), kv));
        const __m256i mm = _mm256_max_epu8(m0, m1);
        const __m256i none = _mm256_cmpeq_epi8(mm, zero);
        off = _mm256_blendv_epi8(mm, off, none);
        chunk = _mm256_blendv_epi8(_mm256_set1_epi8((char)(uint8_t)c), chunk, none);
    }
    rb_st(w.kb_chunk, _mm256_blendv_epi8(rb_ld(w.kb_chunk), chunk, msk));
    rb_st(w.kb_off, _mm256_blendv_epi8(rb_ld(w.kb_off), off, msk));
}

/* rb_dp_wave2 of the NEON section (the default kernel: the fused cell on two rows per step with
 * the direct qe scan), 32 lanes; see there for the pairing argument. The pair (r, k), (r + 1, k + 1)
 * shares the query slot p, so one A / QL / BI load serves both. */
template <bool Sym, bool Asym>
static long rb_dp_wave2(rb_work &w, int W, int NR, int omax, int ominq, int omaskq, bool early)
{
    long computed = 0;
    rb_build_bias(w, NR + W - 1);
    const __m256i tbl = rb_tblb(w);
    const __m256i vOI = _mm256_set1_epi8((char)w.oe_ins), vEI = _mm256_set1_epi8((char)w.e_ins);
    const __m256i vOD = _mm256_set1_epi8((char)w.oe_del), vED = _mm256_set1_epi8((char)w.e_del), m1 = _mm256_set1_epi8(-1);
    uint8_t *Hp = w.H.data(), *Ha = Hp + (size_t)W * 32, *Hb = Hp + (size_t)W * 64, *E = w.E.data();
    const uint8_t *A = w.A.data(), *QL = w.QL.data(), *BI = w.BI.data(), *REF = w.REF.data();
    uint8_t *Rout = w.R.data();
    memset(Hp, 0, (size_t)W * 96);
    memset(E, 0, (size_t)(W + 2) * 32);
    __m256i gmax = _mm256_setzero_si256(), pend = _mm256_setzero_si256();
    const __m256i tgt = rb_ld(w.target);
    int rlast = NR - 1;
    __m256i te_lo = _mm256_setzero_si256(), te_hi = _mm256_setzero_si256();
    auto row_update = [&](int row, __m256i rmax, __m256i rref) -> __m256i {
        rb_st(Rout + (size_t)row * 32, rmax);
        const __m256i imp = _mm256_and_si256(rb_gt_u8(rmax, gmax), _mm256_cmpgt_epi8(rref, m1));
        __m256i flush = _mm256_andnot_si256(imp, pend);
        if (early) flush = _mm256_and_si256(flush, _mm256_cmpeq_epi8(gmax, tgt));
        gmax = _mm256_blendv_epi8(gmax, rmax, imp);
        rb_te_update(te_lo, te_hi, imp, row);
        pend = imp;
        rlast = row;
        return flush;
    };
    /* The kernel constants are captured by value (as in the NEON rb_dp_wave2): captured by reference
     * they would be address-taken, and every rb_st below (a byte store, which may alias anything)
     * would then force the paired-row loop to reload them from the stack on every step. */
    auto single_row = [&, tbl, vOI, vEI, vOD, vED](int r) -> bool {
        const __m256i rref = rb_ld(REF + (size_t)r * 32);
        __m256i f = _mm256_setzero_si256(), rmax = _mm256_setzero_si256();
        const int khi = std::min(W - 1, r + omax), klo = std::max(0, r + ominq);
        computed += khi >= klo ? khi - klo + 1 : 0;
        const size_t p0 = (size_t)(r + (W - 1 - khi)) * 32;
        const uint8_t *ap = A + p0, *qp = QL + p0, *bp = BI + p0;
        const int kun = std::max(klo, r + omaskq);
        for (int k = khi; k >= klo; k--, ap += 32, qp += 32, bp += 32) {
            const __m256i e = rb_ld(E + k * 32);
            const __m256i h0 = _mm256_max_epu8(rb_hs<Asym>(rb_ld(Hp + k * 32), rb_ld(ap), rref, rb_ld(bp), tbl), e);
            rb_st(Ha + k * 32, _mm256_max_epu8(h0, f));
            rmax = _mm256_max_epu8(rmax, k >= kun ? h0 : _mm256_and_si256(h0, rb_ld(qp)));
            const __m256i h07 = _mm256_subs_epu8(h0, vOI), hd = Sym ? h07 : _mm256_subs_epu8(h0, vOD);
            rb_st(E + (k + 1) * 32, _mm256_max_epu8(hd, _mm256_subs_epu8(e, vED)));
            f = _mm256_max_epu8(h07, _mm256_subs_epu8(f, vEI));
        }
        const __m256i flush = row_update(r, rmax, rref);
        if (rb_mask32(flush)) rb_qe_scan(w, Hp, r - 1, W, omax, ominq, flush, gmax);
        std::swap(Hp, Ha);
        return early && rb_all_ge(gmax, tgt);
    };
    int r = 0;
    bool stop = false;
    for (; !stop && r + 1 < NR; r += 2) {
        const int khi_a = std::min(W - 1, r + omax), klo_a = std::max(0, r + ominq);
        const int khi_b = std::min(W - 1, r + 1 + omax), klo_b = std::max(0, r + 1 + ominq);
        if (khi_a < klo_a || khi_b < klo_b) { stop = single_row(r) || single_row(r + 1); continue; }
        const __m256i rref_a = rb_ld(REF + (size_t)r * 32), rref_b = rb_ld(REF + (size_t)(r + 1) * 32);
        __m256i fa = _mm256_setzero_si256(), fb = _mm256_setzero_si256();
        __m256i rmax_a = _mm256_setzero_si256(), rmax_b = _mm256_setzero_si256();
        __m256i diag = _mm256_setzero_si256();   // H of (r, k + 1); zero above row r's range, never written
        const int kun_a = std::max(klo_a, r + omaskq), kun_b = std::max(klo_b, r + 1 + omaskq);
        const int kun = std::max(kun_a, kun_b - 1);   // both cells of a step unmasked iff k >= kun
        const size_t p0 = (size_t)(r + (W - 1 - khi_a)) * 32;
        const uint8_t *ap = A + p0, *qp = QL + p0, *bp = BI + p0;
        int k = khi_a;
        if (khi_a == W - 1) {   // (r, W - 1) alone: its partner (r + 1, W) does not exist
            const __m256i e = rb_ld(E + k * 32);
            const __m256i h0 = _mm256_max_epu8(rb_hs<Asym>(rb_ld(Hp + k * 32), rb_ld(ap), rref_a, rb_ld(bp), tbl), e);
            diag = _mm256_max_epu8(h0, fa);
            rb_st(Ha + k * 32, diag);
            rmax_a = k >= kun_a ? h0 : _mm256_and_si256(h0, rb_ld(qp));
            fa = _mm256_subs_epu8(h0, vOI);   // max(h0 - oe_ins, sat(0 - e_ins))
            k--; ap += 32; qp += 32; bp += 32;
        }
#define RB_W2_STEP(MASK)                                                                        \
        {                                                                                      \
            const __m256i q = rb_ld(ap), e = rb_ld(E + k * 32), bi = rb_ld(bp);                \
            const __m256i h0a = _mm256_max_epu8(rb_hs<Asym>(rb_ld(Hp + k * 32), q, rref_a, bi, tbl), e); \
            const __m256i h07a = _mm256_subs_epu8(h0a, vOI);                                   \
            const __m256i hda = Sym ? h07a : _mm256_subs_epu8(h0a, vOD);                       \
            const __m256i ha = _mm256_max_epu8(h0a, fa);                                       \
            rb_st(Ha + k * 32, ha);                                                            \
            fa = _mm256_max_epu8(h07a, _mm256_subs_epu8(fa, vEI));                             \
            const __m256i ea = _mm256_max_epu8(hda, _mm256_subs_epu8(e, vED));  /* E in of (r+1, k+1) */ \
            const __m256i h0b = _mm256_max_epu8(rb_hs<Asym>(diag, q, rref_b, bi, tbl), ea);          \
            const __m256i h07b = _mm256_subs_epu8(h0b, vOI);                                   \
            const __m256i hdb = Sym ? h07b : _mm256_subs_epu8(h0b, vOD);                       \
            rb_st(Hb + (k + 1) * 32, _mm256_max_epu8(h0b, fb));                                \
            fb = _mm256_max_epu8(h07b, _mm256_subs_epu8(fb, vEI));                             \
            rb_st(E + (k + 2) * 32, _mm256_max_epu8(hdb, _mm256_subs_epu8(ea, vED)));          \
            diag = ha;                                                                         \
            if (MASK) {                                                                        \
                const __m256i m = rb_ld(qp);                                                   \
                rmax_a = _mm256_max_epu8(rmax_a, _mm256_and_si256(h0a, m));                    \
                rmax_b = _mm256_max_epu8(rmax_b, _mm256_and_si256(h0b, m));                    \
            } else {                                                                           \
                rmax_a = _mm256_max_epu8(rmax_a, h0a);                                         \
                rmax_b = _mm256_max_epu8(rmax_b, h0b);                                         \
            }                                                                                  \
        }
        for (; k >= kun; k--, ap += 32, qp += 32, bp += 32) RB_W2_STEP(false)
        for (; k >= klo_a; k--, ap += 32, qp += 32, bp += 32) RB_W2_STEP(true)
#undef RB_W2_STEP
        if (klo_b == klo_a) {   // (r + 1, 0) alone: its partner (r, -1) does not exist; E in = slot 0 = 0
            const __m256i h0b = rb_hs<Asym>(diag, rb_ld(ap), rref_b, rb_ld(bp), tbl);
            rb_st(Hb + klo_a * 32, _mm256_max_epu8(h0b, fb));
            rb_st(E + (klo_a + 1) * 32, _mm256_subs_epu8(h0b, vOD));
            rmax_b = _mm256_max_epu8(rmax_b, _mm256_and_si256(h0b, rb_ld(qp)));
        }
        computed += khi_a - klo_a + 1;
        const __m256i flush_a = row_update(r, rmax_a, rref_a);
        const __m256i gmax_a = gmax;
        if (early && rb_all_ge(gmax_a, tgt)) {
            if (rb_mask32(flush_a)) rb_qe_scan(w, Hp, r - 1, W, omax, ominq, flush_a, gmax_a);
            std::swap(Hp, Ha);
            stop = true;
            break;
        }
        computed += khi_b - klo_b + 1;
        const __m256i flush_b = row_update(r + 1, rmax_b, rref_b);
        if (rb_mask32(_mm256_or_si256(flush_a, flush_b))) {
            if (rb_mask32(flush_a)) rb_qe_scan(w, Hp, r - 1, W, omax, ominq, flush_a, gmax_a);
            if (rb_mask32(flush_b)) rb_qe_scan(w, Ha, r, W, omax, ominq, flush_b, gmax);
        }
        uint8_t *t = Hp; Hp = Hb; Hb = Ha; Ha = t;
        stop = early && rb_all_ge(gmax, tgt);
    }
    if (!stop && r < NR) single_row(r);
    __m256i fin = pend;
    if (early) fin = _mm256_and_si256(fin, _mm256_cmpeq_epi8(gmax, tgt));
    if (rb_mask32(fin)) rb_qe_scan(w, Hp, rlast, W, omax, ominq, fin, gmax);
    rb_te_store(w, te_lo, te_hi);
    rb_st(w.gmax, gmax);
    return computed;
}

}  // namespace

/* ---- run_jobs helpers (the NEON section's, for 32 lanes) ---- */

/* QL[p] = live mask of A[p] (0xFF unless the query code is 0xC0, nonexistent) for p in [0, P). */
static inline void rb_build_ql(const uint8_t *A, uint8_t *QL, int P)
{
    const __m256i vc0 = _mm256_set1_epi8((char)0xC0), m1 = _mm256_set1_epi8(-1);
    for (int p = 0; p < P; p++) rb_st(QL + p * 32, _mm256_xor_si256(_mm256_cmpeq_epi8(rb_ld(A + p * 32), vc0), m1));
}

/* R (NR rows of RB_L interleaved row maxima) back to lane-major: lane l's rows at RL + l * rls,
 * then zeroed past its own nrows[l] up to rls. */
static inline void rb_rows_to_lane_major(const uint8_t *R, int NR, int nl, const int *nrows, uint8_t *RL, int rls)
{
    __m128i a[16], b[16];
    for (int p = 0; p < NR; p += 16) {
        for (int t = 0; t < 16; t++) {
            a[t] = _mm_loadu_si128((const __m128i *)(R + (size_t)(p + t) * 32));
            b[t] = _mm_loadu_si128((const __m128i *)(R + (size_t)(p + t) * 32 + 16));
        }
        rb_transpose16_sse(a);
        rb_transpose16_sse(b);
        for (int l = 0; l < nl && l < 16; l++) _mm_storeu_si128((__m128i *)(RL + (size_t)l * rls + p), a[l]);
        for (int l = 16; l < nl; l++) _mm_storeu_si128((__m128i *)(RL + (size_t)l * rls + p), b[l - 16]);
    }
    for (int l = 0; l < nl; l++) memset(RL + (size_t)l * rls + nrows[l], 0, rls - nrows[l]);
}

/* buf[r] = max(buf[r], Rl[r]) for r in [0, n) in whole 16-byte chunks: both buffers are readable
 * and writable (buf: zero slack, Rl: zero past n) up to n rounded up to 16. */
static inline void rb_merge_max(uint8_t *buf, const uint8_t *Rl, int n)
{
    for (int r = 0; r < n; r += 16)
        _mm_storeu_si128((__m128i *)(buf + r), _mm_max_epu8(_mm_loadu_si128((const __m128i *)(buf + r)),
                                                            _mm_loadu_si128((const __m128i *)(Rl + r))));
}

/* rescue_band_score2 for a CONTIGUOUS row-max array R[0, n) whose slack R[n, n + 16] is zero, 16
 * rows per vector (the NEON section's rb_score2_vec). */
static void rb_score2_vec(const uint8_t *R, int n, int row0, int Z, int te, int minsc, int *score2, int *te2)
{
    const int low = te - Z, high = te + Z;
    int s2 = -1, t2 = -1, bs = -1, bp = -2;
    const __m128i ms = _mm_set1_epi8((char)(uint8_t)minsc);
    alignas(16) uint8_t vb[16];
    for (int k0 = 0; k0 < n; k0 += 16) {
        const __m128i cur = _mm_loadu_si128((const __m128i *)(R + k0)), nxt = _mm_loadu_si128((const __m128i *)(R + k0 + 1));
        /* rising row -> 0: nxt > cur (unsigned) iff max(nxt, cur) != cur */
        const __m128i rise = _mm_xor_si128(_mm_cmpeq_epi8(_mm_max_epu8(nxt, cur), cur), _mm_set1_epi8(-1));
        const __m128i v = _mm_andnot_si128(rise, cur);
        uint32_t m = (uint32_t)_mm_movemask_epi8(_mm_cmpeq_epi8(_mm_max_epu8(v, ms), v));   // v >= minsc
        if (!m) continue;
        _mm_store_si128((__m128i *)vb, v);
        do {
            const int t = __builtin_ctz(m);
            m &= m - 1;
            const int val = vb[t], i = row0 + k0 + t;
            if (bp + 1 != i) {
                if (bp >= 0 && (bp < low || bp > high) && bs > s2) { s2 = bs; t2 = bp; }
                bs = val; bp = i;
            } else if (bs < val) {
                bs = val; bp = i;
            }
        } while (m);
    }
    if (bp >= 0 && (bp < low || bp > high) && bs > s2) { s2 = bs; t2 = bp; }
    *score2 = s2;
    *te2 = t2;
}

#endif
