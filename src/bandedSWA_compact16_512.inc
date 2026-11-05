/* bandedSWA_compact16_512.inc -- same-row lane compaction for the AVX-512BW 16-bit extension
 * kernel. Included from bandedSWA.cpp right after smithWaterman512_16 (avx512bw tier only).
 * getScores16 runs it for batches of more than one lane group when compact16_groups_ > 0 (see
 * BandedPairWiseSW::set_lane_compaction16).
 *
 * The driver is bandedSWA_compact16.inc's (NEON) on 32-lane groups: superblocks of K groups in
 * block lockstep, R rows per block; at block ends bsw_compact_plan() retires the emptiest group
 * whenever the live lanes fit in one group fewer, moving its live lanes into dead slots of the
 * others. The 16-bit kernel keeps absolute int16 columns, so a moved lane's state is a plain copy
 * of its slot. The group's state stays in registers across a block and goes through the
 * per-group struct only at block ends. See bandedSWA_compact16.inc for why the driver is exact.
 *
 * The row body is smithWaterman512_16's row; any change to one must be mirrored in the other
 * (test/unit/test_bandedswa_compact16.cpp compares the two field by field). */
#if __AVX512BW__

#include <algorithm>
#include <cstring>

#include "bsw_compact.h"

static_assert(SIMD_WIDTH16 == 32, "bandedSWA_compact16_512.inc: Bsw16CVec512 holds 32 lanes");

namespace {
struct alignas(64) Bsw16CVec512 {
    int16_t head[32], tail[32], maxend[32], exit0[32], maxScore[32], x[32], y[32], max_off[32];
    int16_t gscore[32], max_ie[32], myband[32], qlen[32], mlen[32], h0[32];
    int tlenw[32];
    SeqPair *pair[32];
    int nbeg, nend, ncol, minq;
    uint16_t *seq1, *seq2;
    int16_t *H, *F;
};
struct Bsw16CBuf512 {
    uint8_t *mem = nullptr; size_t cap = 0;
    ~Bsw16CBuf512() { if (mem) _mm_free(mem); }
};
}  // namespace
static thread_local Bsw16CBuf512 g_bsw16c_buf512;

template <bool CPG>
void BandedPairWiseSW::smithWatermanBatchWrapper16Compact(SeqPair *pairArray,
                                                          uint8_t *seqBufRef,
                                                          uint8_t *seqBufQer,
                                                          int32_t numPairs,
                                                          int32_t w,
                                                          int K,
                                                          int R)
{
    if (K < 1) K = 1;
    if (R < 1) R = 1;
    if (K > BSW_COMPACT_GROUPS_MAX) K = BSW_COMPACT_GROUPS_MAX;
    const int32_t roundNumPairs = ((numPairs + SIMD_WIDTH16 - 1) / SIMD_WIDTH16) * SIMD_WIDTH16;
    for (int32_t ii = numPairs; ii < roundNumPairs; ii++) {   // padding lanes, as the plain wrapper
        pairArray[ii].id = ii; pairArray[ii].len1 = 0; pairArray[ii].len2 = 0;
        pairArray[ii].idr = 0; pairArray[ii].idq = 0; pairArray[ii].h0 = 0;
    }
    const size_t S = (size_t) MAX_SEQ_LEN16 * SIMD_WIDTH16;   /* elements per slab */
    const size_t need = (size_t) K * (sizeof(Bsw16CVec512) + 4 * S * sizeof(int16_t)) + S * sizeof(int16_t) + 128;
    if (g_bsw16c_buf512.cap < need) {
        if (g_bsw16c_buf512.mem) _mm_free(g_bsw16c_buf512.mem);
        g_bsw16c_buf512.mem = (uint8_t *) _mm_malloc(need, 64);
        if (g_bsw16c_buf512.mem == NULL) { fprintf(stderr, "Error! Mem not allocated!!!\n"); exit(EXIT_FAILURE); }
        g_bsw16c_buf512.cap = need;
    }
    Bsw16CVec512 *Vs = (Bsw16CVec512 *) g_bsw16c_buf512.mem;
    int16_t *slab = (int16_t *) (((uintptr_t) (g_bsw16c_buf512.mem + (size_t) K * sizeof(Bsw16CVec512)) + 63) & ~(uintptr_t) 63);
    for (int v = 0; v < K; v++) {
        Vs[v].seq1 = (uint16_t *) (slab + (size_t) (4 * v + 0) * S);
        Vs[v].seq2 = (uint16_t *) (slab + (size_t) (4 * v + 1) * S);
        Vs[v].H    = slab + (size_t) (4 * v + 2) * S;
        Vs[v].F    = slab + (size_t) (4 * v + 3) * S;
    }
    int16_t *sbt_buf = slab + (size_t) (4 * K) * S;   /* one row's scores, reused by every group */

    /* ---- constants (as smithWatermanBatchWrapper16 + smithWaterman512_16) ---- */
    const int eb = end_bonus;
    int16_t maxsc = 0;
    if (maxsc < w_match) maxsc = w_match;
    if (maxsc < w_mismatch) maxsc = w_mismatch;
    if (maxsc < w_ambig) maxsc = w_ambig;
    const int zdrop = this->zdrop;

    const bool gen_mat = bsw_gen_mat_;
    const BswFreedCell fc = bsw_fc_;
    /* the SoA N codes: the LUT pre-pass (symmetric path) needs the asymmetric AMBR16 / AMBQ16 */
    const uint16_t ambRef = gen_mat ? (uint16_t) 0xFFFF : (uint16_t) AMBR16;
    const uint16_t ambQer = gen_mat ? (uint16_t) 0xFFFF : (uint16_t) AMBQ16;
    const __m512i oe_ins512_s = _mm512_set1_epi16(this->o_ins + this->e_ins);   /* the pack's seed row */
    const __m512i e_ins512_s  = _mm512_set1_epi16(this->e_ins);
    const __m512i zero512_s   = _mm512_setzero_si512();

    uint64_t live[BSW_COMPACT_GROUPS_MAX];
    uint8_t active[BSW_COMPACT_GROUPS_MAX];
    BswMove moves[BSW_COMPACT_GROUPS_MAX * SIMD_WIDTH16];

    auto writeback = [&](Bsw16CVec512 &V, int l) {
        SeqPair *p = V.pair[l];
        if (!p) return;
        p->score   = V.maxScore[l];
        p->tle     = V.x[l];
        p->qle     = V.y[l];
        p->max_off = V.max_off[l];
        p->gscore  = V.gscore[l];
        p->gtle    = V.max_ie[l];
        V.pair[l] = nullptr;
    };

    for (int32_t sb0 = 0; sb0 < roundNumPairs; sb0 += K * SIMD_WIDTH16) {
        const int nv = std::min(K, (int) ((roundNumPairs - sb0) / SIMD_WIDTH16));
        int nrow = 0, ncol = 0;
        compact16_counts_.superblocks++;
        for (int q = sb0; q < sb0 + nv * SIMD_WIDTH16; q++) {
            const SeqPair &sp = pairArray[q];
            if (sp.len1 < 0 || sp.len1 >= MAX_SEQ_LEN16)
                err_fatal(__func__, "bandedSWA: ref window length (len1) %d for pair %d is out of range [0, MAX_SEQ_LEN16=%d)", sp.len1, sp.id, MAX_SEQ_LEN16);
            if (sp.len2 < 0 || sp.len2 >= MAX_SEQ_LEN16)
                err_fatal(__func__, "bandedSWA: query length (len2) %d for pair %d is out of range [0, MAX_SEQ_LEN16=%d)", sp.len2, sp.id, MAX_SEQ_LEN16);
            if (nrow < sp.len1) nrow = sp.len1;
            if (ncol < sp.len2) ncol = sp.len2;
        }

        /* ---- pack and seed every group (the plain wrapper's setup, at superblock extents) ---- */
        /* noinline: kept out of the function that runs the rows, so it cannot disturb their
         * register allocation */
        auto pack_group = [&](const int v) __attribute__((noinline)) {
            Bsw16CVec512 &V = Vs[v];
            SeqPair *p = pairArray + sb0 + v * SIMD_WIDTH16;
            const uint8_t *seq1p[SIMD_WIDTH16], *seq2p[SIMD_WIDTH16];
            int len1a[SIMD_WIDTH16], len2a[SIMD_WIDTH16];
            int vcol = 0, minq = 10000000;
            for (int l = 0; l < SIMD_WIDTH16; l++) {
                {   /* prefetch a later pair's sequences, as the plain wrapper does (PFD16 ahead) */
                    const int32_t q = sb0 + v * SIMD_WIDTH16 + l + PFD16;
                    if (q < roundNumPairs) {
                        const SeqPair &spf = pairArray[q];
                        _mm_prefetch((const char *) seqBufRef + (int64_t) spf.idr, _MM_HINT_NTA);
                        _mm_prefetch((const char *) seqBufRef + (int64_t) spf.idr + 64, _MM_HINT_NTA);
                        _mm_prefetch((const char *) seqBufQer + (int64_t) spf.idq, _MM_HINT_NTA);
                        _mm_prefetch((const char *) seqBufQer + (int64_t) spf.idq + 64, _MM_HINT_NTA);
                    }
                }
                const SeqPair &sp = p[l];
                seq1p[l] = seqBufRef + (int64_t) sp.idr; len1a[l] = sp.len1;
                seq2p[l] = seqBufQer + (int64_t) sp.idq; len2a[l] = sp.len2;
                if (vcol < sp.len2) vcol = sp.len2;
                if (sp.len2 < minq) minq = sp.len2;
                const int reach = sp.len2 * maxsc + eb;
                int max_ins = (int) ((double) (reach - o_ins) / e_ins + 1.0); if (max_ins < 1) max_ins = 1;
                int max_del = (int) ((double) (reach - o_del) / e_del + 1.0); if (max_del < 1) max_del = 1;
                int band = w; if (max_ins < band) band = max_ins; if (max_del < band) band = max_del;
                V.myband[l] = (int16_t) band;
                V.h0[l] = (int16_t) (uint16_t) sp.h0;
                V.qlen[l] = (int16_t) sp.len2;
                V.tlenw[l] = sp.len1;
                {   /* mlen = min_epu16(qlen + myband, tlen) */
                    const uint16_t a = (uint16_t) (sp.len2 + (uint16_t) band), t = (uint16_t) sp.len1;
                    V.mlen[l] = (int16_t) (a < t ? a : t);
                }
                V.pair[l] = &p[l];
                V.head[l] = 0; V.tail[l] = (int16_t) sp.len2; V.maxend[l] = (int16_t) 0x8000;
                V.exit0[l] = (int16_t) 0xFFFF; V.maxScore[l] = V.h0[l]; V.x[l] = 0; V.y[l] = 0;
                V.max_off[l] = 0; V.gscore[l] = -1; V.max_ie[l] = 0;
            }
            /* reference / query bases (N -> ambRef / ambQer), then DUMMY1 / DUMMY2 through the
             * superblock's last row / column, as the plain wrapper's tiled transpose */
            x86_soa_pack_u16<SIMD_WIDTH16>(V.seq1, seq1p, len1a, len1a, nrow + 1, DUMMY1, DUMMY1, AMBIG, ambRef);
            x86_soa_pack_u16<SIMD_WIDTH16>(V.seq2, seq2p, len2a, len2a, ncol + 1, DUMMY2, DUMMY2, AMBIG, ambQer);
            V.nbeg = 0; V.nend = vcol; V.ncol = vcol; V.minq = minq - 1;
            /* insertion seed row (H1) as the plain wrapper, through the superblock's last column */
            const __m512i h0_512 = _mm512_load_si512((const __m512i *) V.h0);
            _mm512_store_si512((__m512i *) V.H, h0_512);
            __m512i t = _mm512_mask_blend_epi16(_mm512_cmpgt_epi16_mask(h0_512, oe_ins512_s), zero512_s,
                                                _mm512_sub_epi16(h0_512, oe_ins512_s));
            _mm512_store_si512((__m512i *) (V.H + SIMD_WIDTH16), t);
            for (int k = 2; k <= ncol; k++) {
                t = _mm512_max_epi16(_mm512_sub_epi16(t, e_ins512_s), zero512_s);
                _mm512_store_si512((__m512i *) (V.H + k * SIMD_WIDTH16), t);
            }
            /* H1[maxLen2] = 0; the plain wrapper writes it first, so H1[0] / H1[1] win over it */
            if (vcol >= 2) _mm512_store_si512((__m512i *) (V.H + vcol * SIMD_WIDTH16), zero512_s);
            for (int j = 0; j <= ncol; j++) _mm512_store_si512((__m512i *) (V.F + j * SIMD_WIDTH16), zero512_s);
            active[v] = 1;
        };
        for (int v = 0; v < nv; v++) pack_group(v);

        /* ---- rows [i0, i1) of one group: smithWaterman512_16's row, with the group's state in
         * registers across the block (loaded once, stored once). always_inline: the hot function
         * must be the block loop with the row in it, not a call per row. Returns the live-lane
         * mask after the block (0 when the group is done: every lane finished, or a row whose
         * every lane's maximum is zero, where the plain kernel leaves the group). ---- */
        auto rows = [&](Bsw16CVec512 &V, const int i0, const int i1) __attribute__((always_inline)) -> uint32_t {
            /* local copies of the constants, under the names the kernel's macros use: the row
             * stores through int16 pointers, which may alias captured vectors */
            const __m512i match512    = _mm512_set1_epi16(this->w_match);
            const __m512i mismatch512 = _mm512_set1_epi16(this->w_mismatch);
            const __m512i w_ambig_512 = _mm512_set1_epi16(this->w_ambig);
            const __m512i e_del512    = _mm512_set1_epi16(this->e_del);
            const __m512i oe_del512   = _mm512_set1_epi16(this->o_del + this->e_del);
            const __m512i e_ins512    = _mm512_set1_epi16(this->e_ins);
            const __m512i oe_ins512   = _mm512_set1_epi16(this->o_ins + this->e_ins);
            const __m512i zero512 = _mm512_setzero_si512();
            const __m512i one512  = _mm512_set1_epi16(1);
            const __m512i two512  = _mm512_set1_epi16(2);
            const __m512i ff512   = _mm512_set1_epi16((int16_t) 0xFFFF);
            const __m512i zdrop512 = _mm512_set1_epi16(zdrop);
            const __mmask32 dmask = 0xFFFFFFFF;
            int16_t *const H_h = V.H, *const F = V.F;
            const uint16_t *seq2SoA = V.seq2;
            const int minq = V.minq, vncol = V.ncol;
            __m512i maxScore512 = _mm512_load_si512((const __m512i *) V.maxScore);
            __m512i x512        = _mm512_load_si512((const __m512i *) V.x);
            __m512i y512        = _mm512_load_si512((const __m512i *) V.y);
            __m512i max_off512  = _mm512_load_si512((const __m512i *) V.max_off);
            __m512i gscore      = _mm512_load_si512((const __m512i *) V.gscore);
            __m512i max_ie512   = _mm512_load_si512((const __m512i *) V.max_ie);
            int nbeg = V.nbeg, nend = V.nend;
            __m512i head512     = _mm512_load_si512((const __m512i *) V.head);
            __m512i tail512     = _mm512_load_si512((const __m512i *) V.tail);
            __m512i maxend512   = _mm512_load_si512((const __m512i *) V.maxend);
            __m512i exit0       = _mm512_load_si512((const __m512i *) V.exit0);
            const __m512i myband512 = _mm512_load_si512((const __m512i *) V.myband);
            const __m512i qlen512   = _mm512_load_si512((const __m512i *) V.qlen);
            const __m512i mlen512   = _mm512_load_si512((const __m512i *) V.mlen);
            const __m512i h0_512    = _mm512_load_si512((const __m512i *) V.h0);
            int j;
            for (int i = i0; i < i1; i++) {
            __m512i e11 = zero512;
            __m512i h00, h11, h10;
            __m512i s10 = _mm512_load_si512((const __m512i *) (V.seq1 + i * SIMD_WIDTH16));
            int beg = nbeg, end = nend;
            if (beg < i - w) beg = i - w;
            if (end > i + w + 1) end = i + w + 1;
            if (end > vncol) end = vncol;
            h10 = zero512;
            if (beg == 0) {   /* H2[i + 1]: the deletion seed column, h0 - o_del - (i+1) e_del in int16, >= 0 */
                const int16_t sub = (int16_t) (this->o_del + (i + 1) * this->e_del);
                h10 = _mm512_max_epi16(_mm512_sub_epi16(h0_512, _mm512_set1_epi16(sub)), zero512);
            }
            __m512i j512 = zero512;
            __m512i maxRS1 = zero512;
            __m512i i1_512 = _mm512_set1_epi16(i + 1);
            __m512i y1_512 = zero512;

            __m512i phead512 = head512, ptail512 = tail512;
            head512 = _mm512_max_epi16(head512, _mm512_sub_epi16(_mm512_set1_epi16(i), myband512));
            tail512 = _mm512_min_epu16(tail512, _mm512_add_epi16(i1_512, myband512));
            tail512 = _mm512_min_epu16(tail512, qlen512);
            {   __mmask32 nm_ = _mm512_cmpeq_epi16_mask(_mm512_sub_epi16(tail512, maxend512), two512) & _mm512_movepi16_mask(exit0);
                if (nm_) {
                    int16_t me_a[SIMD_WIDTH16] __attribute((aligned(64)));
                    _mm512_store_si512((__m512i *) me_a, maxend512);
                    const int oe_ = this->o_ins + this->e_ins;
                    while (nm_) { const int l = __builtin_ctz(nm_); nm_ &= nm_ - 1;
                        const int c = me_a[l] + 1;
                        H_h[c * SIMD_WIDTH16 + l] = (int16_t) bsw_seed_row((int16_t) V.h0[l], oe_, this->e_ins, c);
                        F[c * SIMD_WIDTH16 + l] = 0; }
                }
                maxend512 = _mm512_max_epi16(maxend512, tail512);
            }
            {   /* trimming */
                const __mmask32 cmph = _mm512_cmpeq_epi16_mask(head512, phead512) & _mm512_cmpeq_epi16_mask(tail512, ptail512);
                for (int l = beg; l < end && cmph != dmask; l++) {
                    __m512i h512 = _mm512_load_si512((__m512i *)(H_h + l * SIMD_WIDTH16));
                    __m512i f512 = _mm512_load_si512((__m512i *)(F + l * SIMD_WIDTH16));
                    __mmask32 cmp1 = _mm512_cmpgt_epi16_mask(head512, _mm512_set1_epi16(l));
                    if (cmp1 == 0x00) break;
                    cmp1 = cmp1 | _mm512_cmpgt_epi16_mask(_mm512_set1_epi16(l + 1), tail512);
                    _mm512_store_si512((__m512i *)(F + l * SIMD_WIDTH16), _mm512_mask_blend_epi16(cmp1, f512, zero512));
                    _mm512_store_si512((__m512i *)(H_h + l * SIMD_WIDTH16), _mm512_mask_blend_epi16(cmp1, h512, zero512));
                }
            }
            {
                __mmask32 cmpim = _mm512_cmpgt_epi16_mask(i1_512, mlen512);
                cmpim = cmpim | _mm512_cmpeq_epi16_mask(tail512, head512);
                cmpim = cmpim | _mm512_cmpgt_epi16_mask(head512, tail512);
                exit0 = _mm512_mask_blend_epi16(cmpim, exit0, zero512);
            }
            /* the score pre-pass for [beg, end) */
            if (CPG) {
                const __m512i pmat16_512 = _mm512_loadu_si512((const __m512i *) bsw_pmat16_lut_);
                const __m512i cpgref512 = _mm512_set1_epi16(this->cpg_ref_code_);
                const __m512i cpgread512 = _mm512_set1_epi16(this->cpg_read_);
                const __m512i cpgval512 = _mm512_set1_epi16(this->cpg_val_);
                __mmask32 rowcpg;
                BSW_CPG_ROW16_512(s10, rowcpg, cpgref512);
                for (int jp = beg; jp < end; jp++) {
                    __m512i s2 = _mm512_load_si512((const __m512i *)(seq2SoA + jp * SIMD_WIDTH16));
                    __m512i sbt11;
                    SBT_PREPASS16_LUT(s10, s2, sbt11, pmat16_512);
                    BSW_CPG_OVERRIDE16_512(s2, rowcpg, sbt11, cpgval512, cpgread512);
                    _mm512_store_si512((__m512i *)(sbt_buf + jp * SIMD_WIDTH16), sbt11);
                }
            } else if (!gen_mat) {
                const __m512i pmat16_512 = _mm512_loadu_si512((const __m512i *) bsw_pmat16_lut_);
                for (int jp = beg; jp < end; jp++) {
                    __m512i s2 = _mm512_load_si512((const __m512i *)(seq2SoA + jp * SIMD_WIDTH16));
                    __m512i sbt11;
                    SBT_PREPASS16_LUT(s10, s2, sbt11, pmat16_512);
                    _mm512_store_si512((__m512i *)(sbt_buf + jp * SIMD_WIDTH16), sbt11);
                }
            } else if (fc.rank1) {
                const __m512i frref512 = _mm512_set1_epi16(fc.ref), frread512 = _mm512_set1_epi16(fc.read);
                __m512i alt10 = _mm512_mask_blend_epi16(_mm512_cmpeq_epi16_mask(s10, frref512), s10, frread512);
                for (int jp = beg; jp < end; jp++) {
                    __m512i s2 = _mm512_load_si512((const __m512i *)(seq2SoA + jp * SIMD_WIDTH16));
                    __m512i sbt11;
                    SBT_PREPASS16_RANK1(s10, s2, alt10, sbt11, mismatch512, match512, w_ambig_512);
                    _mm512_store_si512((__m512i *)(sbt_buf + jp * SIMD_WIDTH16), sbt11);
                }
            } else {
                const __m512i amat512 = _mm512_broadcast_i32x4(_mm_load_si128((const __m128i *)bsw_amat_bytes_));
                const __m512i three512 = _mm512_set1_epi16(3);
                for (int jp = beg; jp < end; jp++) {
                    __m512i s2 = _mm512_load_si512((const __m512i *)(seq2SoA + jp * SIMD_WIDTH16));
                    __m512i sbt11;
                    SBT_PREPASS16_AMAT(s10, s2, sbt11, amat512, mismatch512, w_ambig_512, three512);
                    _mm512_store_si512((__m512i *)(sbt_buf + jp * SIMD_WIDTH16), sbt11);
                }
            }
#define C16_CELL512(QCAP) \
            { \
                __m512i f11, f21; \
                __m512i sbt11 = _mm512_load_si512((__m512i *)(sbt_buf + j * SIMD_WIDTH16)); \
                h00 = _mm512_load_si512((__m512i *)(H_h + j * SIMD_WIDTH16)); \
                f11 = _mm512_load_si512((__m512i *)(F + j * SIMD_WIDTH16)); \
                __m512i pj512 = j512; \
                j512 = _mm512_add_epi16(j512, one512); \
                MAIN_CODE16_CORE(sbt11, h00, h11, e11, f11, f21, zero512, \
                                 e_ins512, oe_ins512, e_del512, oe_del512); \
                __mmask32 cmp2 = _mm512_cmpgt_epi16_mask(head512, pj512); \
                __mmask32 cmp1 = _mm512_cmpgt_epi16_mask(pj512, tail512) | cmp2; \
                h10 = _mm512_mask_blend_epi16(cmp1, h10, zero512); \
                __m512i bmaxRS = maxRS1; \
                maxRS1 = _mm512_max_epi16(maxRS1, h11); \
                __mmask32 cmpA = _mm512_cmpeq_epi16_mask(maxRS1, h11); \
                cmp1 = _mm512_cmpgt_epi16_mask(j512, tail512); \
                f21 = _mm512_mask_blend_epi16(cmp1, f21, zero512); \
                cmp1 = cmp1 | cmp2; \
                __m512i blend512 = _mm512_mask_blend_epi16(cmpA, y1_512, j512); \
                y1_512 = _mm512_mask_blend_epi16(cmp1, blend512, y1_512); \
                maxRS1 = _mm512_mask_blend_epi16(cmp1, maxRS1, bmaxRS); \
                _mm512_store_si512((__m512i *)(F + j * SIMD_WIDTH16), f21); \
                _mm512_store_si512((__m512i *)(H_h + j * SIMD_WIDTH16), h10); \
                h10 = h11; \
                if (QCAP) hqe512 = _mm512_mask_mov_epi16(hqe512, _mm512_cmpeq_epi16_mask(j512, qlen512), h11); \
            }
            j512 = _mm512_set1_epi16(beg);
            __m512i hqe512 = zero512;
            {
                const int jq_ = minq < beg ? beg : (minq > end ? end : minq);
                for (j = beg; j < jq_; j++) C16_CELL512(0)
                for (j = jq_; j < end; j++) C16_CELL512(1)
            }
#undef C16_CELL512
            {
                __mmask32 cmp1 = _mm512_cmpgt_epi16_mask(j512, tail512) | _mm512_cmpgt_epi16_mask(head512, j512);
                h10 = _mm512_mask_blend_epi16(cmp1, h10, zero512);
                _mm512_store_si512((__m512i *)(H_h + j * SIMD_WIDTH16), h10);
                _mm512_store_si512((__m512i *)(F + j * SIMD_WIDTH16), zero512);
            }
            {   /* query-end fold, before the zero-row exit (see smithWaterman512_16) */
                __mmask32 qin = _mm512_cmpgt_epi16_mask(qlen512, _mm512_set1_epi16(beg))
                              & ~_mm512_cmpgt_epi16_mask(qlen512, _mm512_set1_epi16(end));
                __mmask32 upd = qin & _mm512_movepi16_mask(exit0)
                              & ~_mm512_cmpgt_epi16_mask(qlen512, tail512);
                __mmask32 cmp_gh = _mm512_cmpgt_epi16_mask(gscore, hqe512);
                __m512i cand_ie = _mm512_mask_blend_epi16(cmp_gh, i1_512, max_ie512);
                max_ie512 = _mm512_mask_blend_epi16(upd, max_ie512, cand_ie);
                gscore = _mm512_mask_blend_epi16(upd, gscore, _mm512_max_epi16(gscore, hqe512));
            }

            __m512i bmaxScore512 = maxScore512;
            __mmask32 tmp = _mm512_cmpeq_epi16_mask(maxRS1, zero512);
            if (tmp == dmask) {
                /* every lane's row max is zero: the plain kernel leaves the group here; every
                 * lane's results are final */
                exit0 = zero512;
                break;
            }
            exit0 = _mm512_mask_blend_epi16(tmp, exit0, zero512);
            __m512i score512 = _mm512_max_epi16(maxScore512, maxRS1);
            maxScore512 = _mm512_mask_blend_epi16(_mm512_movepi16_mask(exit0), maxScore512, score512);
            __mmask32 cmp = _mm512_cmpgt_epi16_mask(maxScore512, bmaxScore512);
            y512 = _mm512_mask_blend_epi16(cmp, y512, y1_512);
            x512 = _mm512_mask_blend_epi16(cmp, x512, i1_512);
            {
                __m512i ind512 = _mm512_abs_epi16(_mm512_sub_epi16(y1_512, i1_512));
                __m512i bmax_off512 = max_off512;
                ind512 = _mm512_max_epi16(max_off512, ind512);
                max_off512 = _mm512_mask_blend_epi16(cmp, bmax_off512, ind512);
            }
            ZSCORE16(i1_512, y1_512);

            {   /* narrowing: smithWaterman512_16's four scans */
                int l;
                for (l = beg; l < end; l++) {
                    __m512i t_ = _mm512_or_si512(_mm512_load_si512((__m512i *)(F + l * SIMD_WIDTH16)),
                                                 _mm512_load_si512((__m512i *)(H_h + l * SIMD_WIDTH16)));
                    if (_mm512_cmpeq_epi16_mask(t_, zero512) == dmask) nbeg = l;
                    else break;
                }
                for (l = end; l >= beg; l--) {
                    __m512i t_ = _mm512_or_si512(_mm512_load_si512((__m512i *)(F + l * SIMD_WIDTH16)),
                                                 _mm512_load_si512((__m512i *)(H_h + l * SIMD_WIDTH16)));
                    if (_mm512_cmpeq_epi16_mask(t_, zero512) != dmask) break;
                }
                nend = l + 2 < vncol ? l + 2 : vncol;
                __m512i exit1 = _mm512_xor_si512(exit0, ff512);
                __mmask32 tmpb = dmask;
                __m512i l512 = _mm512_set1_epi16(beg);
                for (l = beg; l < end; l++) {
                    __m512i t_ = _mm512_or_si512(_mm512_load_si512((__m512i *)(F + l * SIMD_WIDTH16)),
                                                 _mm512_load_si512((__m512i *)(H_h + l * SIMD_WIDTH16)));
                    __mmask32 m_ = _mm512_cmpeq_epi16_mask(_mm512_or_si512(t_, exit1), zero512);
                    if (m_ == 0x00) break;
                    m_ = m_ & tmpb;
                    l512 = _mm512_add_epi16(l512, one512);
                    head512 = _mm512_mask_blend_epi16(m_, head512, l512);
                    tmpb = m_;
                }
                __m512i index512 = tail512;
                tmpb = dmask;
                l512 = _mm512_set1_epi16(end);
                for (l = end; l >= beg; l--) {
                    __m512i t_ = _mm512_or_si512(_mm512_load_si512((__m512i *)(F + l * SIMD_WIDTH16)),
                                                 _mm512_load_si512((__m512i *)(H_h + l * SIMD_WIDTH16)));
                    __mmask32 m_ = _mm512_cmpeq_epi16_mask(_mm512_or_si512(t_, exit1), zero512);
                    if (m_ == 0x00) break;
                    m_ = m_ & tmpb;
                    l512 = _mm512_sub_epi16(l512, one512);
                    index512 = _mm512_mask_blend_epi16(m_, index512, l512);
                    tmpb = m_;
                }
                index512 = _mm512_add_epi16(index512, two512);
                tail512 = _mm512_min_epi16(index512, qlen512);
            }
            if (_mm512_movepi16_mask(exit0) == 0) break;
            }   /* rows */
            _mm512_store_si512((__m512i *) V.head, head512);
            _mm512_store_si512((__m512i *) V.tail, tail512);
            _mm512_store_si512((__m512i *) V.maxend, maxend512);
            _mm512_store_si512((__m512i *) V.exit0, exit0);
            _mm512_store_si512((__m512i *) V.maxScore, maxScore512);
            _mm512_store_si512((__m512i *) V.x, x512);
            _mm512_store_si512((__m512i *) V.y, y512);
            _mm512_store_si512((__m512i *) V.max_off, max_off512);
            _mm512_store_si512((__m512i *) V.gscore, gscore);
            _mm512_store_si512((__m512i *) V.max_ie, max_ie512);
            V.nbeg = nbeg; V.nend = nend;
            return (uint32_t) _mm512_movepi16_mask(exit0);
        };

        /* ---- move one live lane into a dead slot of another group, between rows (next row r);
         * as move_one in bandedSWA_compact16.inc ---- */
        auto move_one = [&](Bsw16CVec512 &Sv, Bsw16CVec512 &D, int sl, int dl, int r) __attribute__((noinline)) {
            writeback(D, dl);
            D.head[dl] = Sv.head[sl]; D.tail[dl] = Sv.tail[sl]; D.maxend[dl] = Sv.maxend[sl];
            D.exit0[dl] = Sv.exit0[sl]; D.maxScore[dl] = Sv.maxScore[sl]; D.x[dl] = Sv.x[sl];
            D.y[dl] = Sv.y[sl]; D.max_off[dl] = Sv.max_off[sl]; D.gscore[dl] = Sv.gscore[sl];
            D.max_ie[dl] = Sv.max_ie[sl]; D.myband[dl] = Sv.myband[sl]; D.qlen[dl] = Sv.qlen[sl];
            D.mlen[dl] = Sv.mlen[sl]; D.h0[dl] = Sv.h0[sl]; D.tlenw[dl] = Sv.tlenw[sl];
            D.pair[dl] = Sv.pair[sl];
            if (Sv.qlen[sl] > D.ncol) D.ncol = Sv.qlen[sl];
            if (Sv.qlen[sl] - 1 < D.minq) D.minq = Sv.qlen[sl] - 1;
            if (Sv.nbeg < D.nbeg) D.nbeg = Sv.nbeg;
            if (Sv.nend > D.nend) D.nend = Sv.nend;
            int ha = Sv.head[sl]; if (ha < 0) ha = 0;
            const int cstart = D.nbeg < 0 ? 0 : D.nbeg;
            for (int c = cstart; c < ha && c <= ncol; c++) {
                D.H[c * SIMD_WIDTH16 + dl] = 0; D.F[c * SIMD_WIDTH16 + dl] = 0;
            }
            for (int c = ha; c <= ncol; c++) {
                D.H[c * SIMD_WIDTH16 + dl] = Sv.H[c * SIMD_WIDTH16 + sl];
                D.F[c * SIMD_WIDTH16 + dl] = Sv.F[c * SIMD_WIDTH16 + sl];
                D.seq2[c * SIMD_WIDTH16 + dl] = Sv.seq2[c * SIMD_WIDTH16 + sl];
            }
            for (int rr = r; rr <= nrow; rr++) D.seq1[rr * SIMD_WIDTH16 + dl] = Sv.seq1[rr * SIMD_WIDTH16 + sl];
            Sv.pair[sl] = nullptr; Sv.exit0[sl] = 0;
        };

        /* ---- block lockstep over the superblock (see bandedSWA_compact16.inc) ---- */
        int prev_live = -1;
        for (int i0 = 0; i0 < nrow; i0 += R) {
            const int i1 = std::min(i0 + R, nrow);
            int total = 0, nact = 0;
            for (int v = 0; v < nv; v++) {
                if (!active[v]) { live[v] = 0; continue; }
                live[v] = rows(Vs[v], i0, i1);
                total += __builtin_popcountll(live[v]);
                if (live[v]) nact++; else active[v] = 0;
            }
            if (nact == 0) break;
            if (nv > 1 && total != prev_live && nact > (total + SIMD_WIDTH16 - 1) / SIMD_WIDTH16) {
                const int nm = bsw_compact_plan(live, active, nv, SIMD_WIDTH16, moves);
                for (int q = 0; q < nm; q++)
                    move_one(Vs[moves[q].sv], Vs[moves[q].dv], moves[q].sl, moves[q].dl, i1);
                int nleft = 0;
                for (int v = 0; v < nv; v++) nleft += active[v];
                compact16_counts_.moves += (uint64_t) nm;
                compact16_counts_.retired += (uint64_t) (nact - nleft);
            }
            prev_live = total;
        }
        for (int v = 0; v < nv; v++)
            for (int l = 0; l < SIMD_WIDTH16; l++) writeback(Vs[v], l);
    }
}

#endif  /* __AVX512BW__ */
