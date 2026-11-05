/* bandedSWA_compact16_256.inc -- same-row lane compaction for the AVX2 16-bit extension kernel.
 * Included from bandedSWA.cpp right after smithWaterman256_16 (avx2 tier only). getScores16 runs
 * it for batches of more than one lane group when compact16_groups_ > 0 (see
 * BandedPairWiseSW::set_lane_compaction16).
 *
 * The driver is bandedSWA_compact16.inc's (NEON) on 16-lane groups: superblocks of K groups in
 * block lockstep, R rows per block; at block ends bsw_compact_plan() retires the emptiest group
 * whenever the live lanes fit in one group fewer, moving its live lanes into dead slots of the
 * others. The 16-bit kernel keeps absolute int16 columns, so a moved lane's state is a plain copy
 * of its slot. The group's state stays in registers across a block and goes through the
 * per-group struct only at block ends. See bandedSWA_compact16.inc for why the driver is exact.
 *
 * The row body is smithWaterman256_16's row; any change to one must be mirrored in the other
 * (test/unit/test_bandedswa_compact16.cpp compares the two field by field). */
#if ((!__AVX512BW__) && (__AVX2__))

#include <algorithm>
#include <cstring>

#include "bsw_compact.h"

static_assert(SIMD_WIDTH16 == 16, "bandedSWA_compact16_256.inc: Bsw16CVec256 holds 16 lanes");

namespace {
struct alignas(64) Bsw16CVec256 {
    int16_t head[16], tail[16], maxend[16], exit0[16], maxScore[16], x[16], y[16], max_off[16];
    int16_t gscore[16], max_ie[16], myband[16], qlen[16], mlen[16], h0[16];
    int tlenw[16];
    SeqPair *pair[16];
    int nbeg, nend, ncol, minq;
    uint16_t *seq1, *seq2;
    int16_t *H, *F;
};
struct Bsw16CBuf256 {
    uint8_t *mem = nullptr; size_t cap = 0;
    ~Bsw16CBuf256() { if (mem) _mm_free(mem); }
};
}  // namespace
static thread_local Bsw16CBuf256 g_bsw16c_buf256;

/* live-lane bits of an int16 all-ones / zero lane mask, lane l at bit l */
static inline uint32_t bsw16c_lane_mask256(__m256i m)
{
    return (uint32_t) _mm_movemask_epi8(_mm_packs_epi16(_mm256_castsi256_si128(m), _mm256_extracti128_si256(m, 1)));
}

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
    const size_t need = (size_t) K * (sizeof(Bsw16CVec256) + 4 * S * sizeof(int16_t)) + S * sizeof(int16_t) + 128;
    if (g_bsw16c_buf256.cap < need) {
        if (g_bsw16c_buf256.mem) _mm_free(g_bsw16c_buf256.mem);
        g_bsw16c_buf256.mem = (uint8_t *) _mm_malloc(need, 64);
        if (g_bsw16c_buf256.mem == NULL) { fprintf(stderr, "Error! Mem not allocated!!!\n"); exit(EXIT_FAILURE); }
        g_bsw16c_buf256.cap = need;
    }
    Bsw16CVec256 *Vs = (Bsw16CVec256 *) g_bsw16c_buf256.mem;
    int16_t *slab = (int16_t *) (((uintptr_t) (g_bsw16c_buf256.mem + (size_t) K * sizeof(Bsw16CVec256)) + 63) & ~(uintptr_t) 63);
    for (int v = 0; v < K; v++) {
        Vs[v].seq1 = (uint16_t *) (slab + (size_t) (4 * v + 0) * S);
        Vs[v].seq2 = (uint16_t *) (slab + (size_t) (4 * v + 1) * S);
        Vs[v].H    = slab + (size_t) (4 * v + 2) * S;
        Vs[v].F    = slab + (size_t) (4 * v + 3) * S;
    }
    int16_t *sbt_buf = slab + (size_t) (4 * K) * S;   /* one row's scores, reused by every group */

    /* ---- constants (as smithWatermanBatchWrapper16 + smithWaterman256_16) ---- */
    const int eb = end_bonus;
    int16_t maxsc = 0;
    if (maxsc < w_match) maxsc = w_match;
    if (maxsc < w_mismatch) maxsc = w_mismatch;
    if (maxsc < w_ambig) maxsc = w_ambig;
    const int zdrop = this->zdrop;

    const bool gen_mat = bsw_gen_mat_;
    const BswFreedCell fc = bsw_fc_;
    const __m256i oe_ins256_s = _mm256_set1_epi16(this->o_ins + this->e_ins);   /* the pack's seed row */
    const __m256i e_ins256_s  = _mm256_set1_epi16(this->e_ins);
    const __m256i zero256_s   = _mm256_setzero_si256();

    uint64_t live[BSW_COMPACT_GROUPS_MAX];
    uint8_t active[BSW_COMPACT_GROUPS_MAX];
    BswMove moves[BSW_COMPACT_GROUPS_MAX * SIMD_WIDTH16];

    auto writeback = [&](Bsw16CVec256 &V, int l) {
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
            Bsw16CVec256 &V = Vs[v];
            SeqPair *p = pairArray + sb0 + v * SIMD_WIDTH16;
            const uint8_t *seq1p[SIMD_WIDTH16], *seq2p[SIMD_WIDTH16];
            int len1a[SIMD_WIDTH16], len2a[SIMD_WIDTH16];
            int vcol = 0, minq = 10000000;
            for (int l = 0; l < SIMD_WIDTH16; l++) {
                {   /* prefetch a later pair's sequences, as the plain wrapper does (PFD ahead) */
                    const int32_t q = sb0 + v * SIMD_WIDTH16 + l + PFD;
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
            /* reference / query bases (N -> 0xFFFF), then DUMMY1 / DUMMY2 through the superblock's
             * last row / column, as the plain wrapper's tiled transpose */
            x86_soa_pack_u16<SIMD_WIDTH16>(V.seq1, seq1p, len1a, len1a, nrow + 1, DUMMY1, DUMMY1, AMBIG, 0xFFFF);
            x86_soa_pack_u16<SIMD_WIDTH16>(V.seq2, seq2p, len2a, len2a, ncol + 1, DUMMY2, DUMMY2, AMBIG, 0xFFFF);
            V.nbeg = 0; V.nend = vcol; V.ncol = vcol; V.minq = minq - 1;
            /* insertion seed row (H1) as the plain wrapper, through the superblock's last column */
            const __m256i h0_256 = _mm256_load_si256((const __m256i *) V.h0);
            _mm256_store_si256((__m256i *) V.H, h0_256);
            __m256i t = _mm256_blendv_epi8(zero256_s, _mm256_sub_epi16(h0_256, oe_ins256_s), _mm256_cmpgt_epi16(h0_256, oe_ins256_s));
            _mm256_store_si256((__m256i *) (V.H + SIMD_WIDTH16), t);
            for (int k = 2; k <= ncol; k++) {
                t = _mm256_max_epi16(_mm256_sub_epi16(t, e_ins256_s), zero256_s);
                _mm256_store_si256((__m256i *) (V.H + k * SIMD_WIDTH16), t);
            }
            /* H1[maxLen2] = 0; the plain wrapper writes it first, so H1[0] / H1[1] win over it */
            if (vcol >= 2) _mm256_store_si256((__m256i *) (V.H + vcol * SIMD_WIDTH16), zero256_s);
            for (int j = 0; j <= ncol; j++) _mm256_store_si256((__m256i *) (V.F + j * SIMD_WIDTH16), zero256_s);
            active[v] = 1;
        };
        for (int v = 0; v < nv; v++) pack_group(v);

        /* ---- rows [i0, i1) of one group: smithWaterman256_16's row, with the group's state in
         * registers across the block (loaded once, stored once). always_inline: the hot function
         * must be the block loop with the row in it, not a call per row. Returns the live-lane
         * mask after the block (0 when the group is done: every lane finished, or a row whose
         * every lane's maximum is zero, where the plain kernel leaves the group). ---- */
        auto rows = [&](Bsw16CVec256 &V, const int i0, const int i1) __attribute__((always_inline)) -> uint32_t {
            /* local copies of the constants, under the names the kernel's macros use: the row
             * stores through int16 pointers, which may alias captured vectors */
            const __m256i match256    = _mm256_set1_epi16(this->w_match);
            const __m256i mismatch256 = _mm256_set1_epi16(this->w_mismatch);
            const __m256i w_ambig_256 = _mm256_set1_epi16(this->w_ambig);
            const __m256i e_del256    = _mm256_set1_epi16(this->e_del);
            const __m256i oe_del256   = _mm256_set1_epi16(this->o_del + this->e_del);
            const __m256i e_ins256    = _mm256_set1_epi16(this->e_ins);
            const __m256i oe_ins256   = _mm256_set1_epi16(this->o_ins + this->e_ins);
            const __m256i zero256 = _mm256_setzero_si256();
            const __m256i one256  = _mm256_set1_epi16(1);
            const __m256i two256  = _mm256_set1_epi16(2);
            const __m256i ff256   = _mm256_set1_epi16(0xFFFF);
            const __m256i zdrop256 = _mm256_set1_epi16(zdrop);
            const __mmask32 dmask32 = 0xAAAAAAAA;
            int16_t *const H_h = V.H, *const F = V.F;
            const uint16_t *seq2SoA = V.seq2;
            const int minq = V.minq, vncol = V.ncol;
            __m256i maxScore256 = _mm256_load_si256((const __m256i *) V.maxScore);
            __m256i x256        = _mm256_load_si256((const __m256i *) V.x);
            __m256i y256        = _mm256_load_si256((const __m256i *) V.y);
            __m256i max_off256  = _mm256_load_si256((const __m256i *) V.max_off);
            __m256i gscore      = _mm256_load_si256((const __m256i *) V.gscore);
            __m256i max_ie256   = _mm256_load_si256((const __m256i *) V.max_ie);
            int nbeg = V.nbeg, nend = V.nend;
            __m256i head256     = _mm256_load_si256((const __m256i *) V.head);
            __m256i tail256     = _mm256_load_si256((const __m256i *) V.tail);
            __m256i maxend256   = _mm256_load_si256((const __m256i *) V.maxend);
            __m256i exit0       = _mm256_load_si256((const __m256i *) V.exit0);
            const __m256i myband256 = _mm256_load_si256((const __m256i *) V.myband);
            const __m256i qlen256   = _mm256_load_si256((const __m256i *) V.qlen);
            const __m256i mlen256   = _mm256_load_si256((const __m256i *) V.mlen);
            const __m256i h0_256    = _mm256_load_si256((const __m256i *) V.h0);
            int j;
            for (int i = i0; i < i1; i++) {
            __m256i e11 = zero256;
            __m256i h00, h11, h10;
            __m256i s10 = _mm256_load_si256((const __m256i *) (V.seq1 + i * SIMD_WIDTH16));
            int beg = nbeg, end = nend;
            if (beg < i - w) beg = i - w;
            if (end > i + w + 1) end = i + w + 1;
            if (end > vncol) end = vncol;
            h10 = zero256;
            if (beg == 0) {   /* H2[i + 1]: the deletion seed column, h0 - o_del - (i+1) e_del in int16, >= 0 */
                const int16_t sub = (int16_t) (this->o_del + (i + 1) * this->e_del);
                h10 = _mm256_max_epi16(_mm256_sub_epi16(h0_256, _mm256_set1_epi16(sub)), zero256);
            }
            __m256i j256 = zero256;
            __m256i maxRS1 = zero256;
            __m256i i1_256 = _mm256_set1_epi16(i + 1);
            __m256i y1_256 = zero256;

            __m256i phead256 = head256, ptail256 = tail256;
            head256 = _mm256_max_epi16(head256, _mm256_sub_epi16(_mm256_set1_epi16(i), myband256));
            tail256 = _mm256_min_epu16(tail256, _mm256_add_epi16(i1_256, myband256));
            tail256 = _mm256_min_epu16(tail256, qlen256);
            {   __m256i need_ = _mm256_and_si256(_mm256_cmpeq_epi16(_mm256_sub_epi16(tail256, maxend256), two256), exit0);
                uint32_t nm_ = (uint32_t) _mm256_movemask_epi8(need_) & 0x55555555u;
                if (nm_) {
                    int16_t me_a[SIMD_WIDTH16] __attribute((aligned(32)));
                    _mm256_store_si256((__m256i *) me_a, maxend256);
                    const int oe_ = this->o_ins + this->e_ins;
                    while (nm_) { const int l = __builtin_ctz(nm_) >> 1; nm_ &= nm_ - 1;
                        const int c = me_a[l] + 1;
                        H_h[c * SIMD_WIDTH16 + l] = (int16_t) bsw_seed_row((int16_t) V.h0[l], oe_, this->e_ins, c);
                        F[c * SIMD_WIDTH16 + l] = 0; }
                }
                maxend256 = _mm256_max_epi16(maxend256, tail256);
            }
            {   /* trimming */
                __m256i cmph = _mm256_and_si256(_mm256_cmpeq_epi16(head256, phead256), _mm256_cmpeq_epi16(tail256, ptail256));
                const uint32_t cmp_ht = (uint32_t) _mm256_movemask_epi8(cmph) & dmask32;
                for (int l = beg; l < end && cmp_ht != dmask32; l++) {
                    __m256i h256 = _mm256_load_si256((__m256i *)(H_h + l * SIMD_WIDTH16));
                    __m256i f256 = _mm256_load_si256((__m256i *)(F + l * SIMD_WIDTH16));
                    __m256i cmp1 = _mm256_cmpgt_epi16(head256, _mm256_set1_epi16(l));
                    if (((uint32_t) _mm256_movemask_epi8(cmp1) & dmask32) == 0) break;
                    cmp1 = _mm256_or_si256(cmp1, _mm256_cmpgt_epi16(_mm256_set1_epi16(l + 1), tail256));
                    _mm256_store_si256((__m256i *)(F + l * SIMD_WIDTH16), _mm256_andnot_si256(cmp1, f256));
                    _mm256_store_si256((__m256i *)(H_h + l * SIMD_WIDTH16), _mm256_andnot_si256(cmp1, h256));
                }
            }
            {
                __m256i cmpim = _mm256_cmpgt_epi16(i1_256, mlen256);
                cmpim = _mm256_or_si256(cmpim, _mm256_cmpeq_epi16(tail256, head256));
                cmpim = _mm256_or_si256(cmpim, _mm256_cmpgt_epi16(head256, tail256));
                exit0 = _mm256_andnot_si256(cmpim, exit0);
            }
            /* the score pre-pass for [beg, end) */
            if (CPG) {
                const __m256i cpgref256 = _mm256_set1_epi16(this->cpg_ref_code_);
                const __m256i cpgread256 = _mm256_set1_epi16(this->cpg_read_);
                const __m256i cpgval256 = _mm256_set1_epi16(this->cpg_val_);
                __m256i rowcpg;
                BSW_CPG_ROW16_256(s10, rowcpg, cpgref256);
                for (int jp = beg; jp < end; jp++) {
                    __m256i s2 = _mm256_load_si256((const __m256i *)(seq2SoA + jp * SIMD_WIDTH16));
                    __m256i sbt11;
                    SBT_PREPASS16_SYM(s10, s2, sbt11, mismatch256, match256, w_ambig_256);
                    BSW_CPG_OVERRIDE16_256(s2, rowcpg, sbt11, cpgval256, cpgread256);
                    _mm256_store_si256((__m256i *)(sbt_buf + jp * SIMD_WIDTH16), sbt11);
                }
            } else if (!gen_mat) {
                for (int jp = beg; jp < end; jp++) {
                    __m256i s2 = _mm256_load_si256((const __m256i *)(seq2SoA + jp * SIMD_WIDTH16));
                    __m256i sbt11;
                    SBT_PREPASS16_SYM(s10, s2, sbt11, mismatch256, match256, w_ambig_256);
                    _mm256_store_si256((__m256i *)(sbt_buf + jp * SIMD_WIDTH16), sbt11);
                }
            } else if (fc.rank1) {
                const __m256i frref256 = _mm256_set1_epi16(fc.ref), frread256 = _mm256_set1_epi16(fc.read);
                __m256i alt10 = _mm256_blendv_epi8(s10, frread256, _mm256_cmpeq_epi16(s10, frref256));
                for (int jp = beg; jp < end; jp++) {
                    __m256i s2 = _mm256_load_si256((const __m256i *)(seq2SoA + jp * SIMD_WIDTH16));
                    __m256i sbt11;
                    SBT_PREPASS16_RANK1(s10, s2, alt10, sbt11, mismatch256, match256, w_ambig_256);
                    _mm256_store_si256((__m256i *)(sbt_buf + jp * SIMD_WIDTH16), sbt11);
                }
            } else {
                const __m256i amat256 = _mm256_broadcastsi128_si256(_mm_load_si128((const __m128i *)bsw_amat_bytes_));
                const __m256i three256 = _mm256_set1_epi16(3);
                for (int jp = beg; jp < end; jp++) {
                    __m256i s2 = _mm256_load_si256((const __m256i *)(seq2SoA + jp * SIMD_WIDTH16));
                    __m256i sbt11;
                    SBT_PREPASS16_AMAT(s10, s2, sbt11, amat256, mismatch256, w_ambig_256, three256);
                    _mm256_store_si256((__m256i *)(sbt_buf + jp * SIMD_WIDTH16), sbt11);
                }
            }
#define C16_CELL256(QCAP) \
            { \
                __m256i f11, f21; \
                __m256i sbt11 = _mm256_load_si256((__m256i *)(sbt_buf + j * SIMD_WIDTH16)); \
                h00 = _mm256_load_si256((__m256i *)(H_h + j * SIMD_WIDTH16)); \
                f11 = _mm256_load_si256((__m256i *)(F + j * SIMD_WIDTH16)); \
                __m256i pj256 = j256; \
                j256 = _mm256_add_epi16(j256, one256); \
                MAIN_CODE16_CORE(sbt11, h00, h11, e11, f11, f21, zero256, \
                                 e_ins256, oe_ins256, e_del256, oe_del256); \
                __m256i cmp2 = _mm256_cmpgt_epi16(head256, pj256); \
                __m256i cmp1 = _mm256_or_si256(_mm256_cmpgt_epi16(pj256, tail256), cmp2); \
                h10 = _mm256_andnot_si256(cmp1, h10); \
                __m256i bmaxRS = maxRS1; \
                maxRS1 = _mm256_max_epi16(maxRS1, h11); \
                __m256i cmpA = _mm256_cmpeq_epi16(maxRS1, h11); \
                cmp1 = _mm256_cmpgt_epi16(j256, tail256); \
                f21 = _mm256_andnot_si256(cmp1, f21); \
                cmp1 = _mm256_or_si256(cmp1, cmp2); \
                cmpA = _mm256_blendv_epi8(y1_256, j256, cmpA); \
                y1_256 = _mm256_blendv_epi8(cmpA, y1_256, cmp1); \
                maxRS1 = _mm256_blendv_epi8(maxRS1, bmaxRS, cmp1); \
                _mm256_store_si256((__m256i *)(F + j * SIMD_WIDTH16), f21); \
                _mm256_store_si256((__m256i *)(H_h + j * SIMD_WIDTH16), h10); \
                h10 = h11; \
                if (QCAP) hqe256 = _mm256_blendv_epi8(hqe256, h11, _mm256_cmpeq_epi16(j256, qlen256)); \
            }
            j256 = _mm256_set1_epi16(beg);
            __m256i hqe256 = zero256;
            {
                const int jq_ = minq < beg ? beg : (minq > end ? end : minq);
                for (j = beg; j < jq_; j++) C16_CELL256(0)
                for (j = jq_; j < end; j++) C16_CELL256(1)
            }
#undef C16_CELL256
            {
                __m256i cmp1 = _mm256_or_si256(_mm256_cmpgt_epi16(j256, tail256), _mm256_cmpgt_epi16(head256, j256));
                h10 = _mm256_andnot_si256(cmp1, h10);
                _mm256_store_si256((__m256i *)(H_h + j * SIMD_WIDTH16), h10);
                _mm256_store_si256((__m256i *)(F + j * SIMD_WIDTH16), zero256);
            }
            {   /* query-end fold, before the zero-row exit (see smithWaterman256_16) */
                __m256i qin = _mm256_andnot_si256(_mm256_cmpgt_epi16(qlen256, _mm256_set1_epi16(end)),
                                                  _mm256_cmpgt_epi16(qlen256, _mm256_set1_epi16(beg)));
                __m256i upd = _mm256_andnot_si256(_mm256_cmpgt_epi16(qlen256, tail256), _mm256_and_si256(qin, exit0));
                __m256i cmp_gh = _mm256_cmpgt_epi16(gscore, hqe256);
                __m256i cand_ie = _mm256_blendv_epi8(i1_256, max_ie256, cmp_gh);
                max_ie256 = _mm256_blendv_epi8(max_ie256, cand_ie, upd);
                gscore = _mm256_blendv_epi8(gscore, _mm256_max_epi16(gscore, hqe256), upd);
            }

            __m256i bmaxScore256 = maxScore256;
            __m256i tmp = _mm256_cmpeq_epi16(maxRS1, zero256);
            if (((uint32_t) _mm256_movemask_epi8(tmp) & dmask32) == dmask32) {
                /* every lane's row max is zero: the plain kernel leaves the group here; every
                 * lane's results are final */
                exit0 = zero256;
                break;
            }
            exit0 = _mm256_andnot_si256(tmp, exit0);
            __m256i score256 = _mm256_max_epi16(maxScore256, maxRS1);
            maxScore256 = _mm256_blendv_epi8(maxScore256, score256, exit0);
            __m256i cmp = _mm256_cmpgt_epi16(maxScore256, bmaxScore256);
            y256 = _mm256_blendv_epi8(y256, y1_256, cmp);
            x256 = _mm256_blendv_epi8(x256, i1_256, cmp);
            {
                tmp = _mm256_abs_epi16(_mm256_sub_epi16(y1_256, i1_256));
                __m256i bmax_off256 = max_off256;
                tmp = _mm256_max_epi16(max_off256, tmp);
                max_off256 = _mm256_blendv_epi8(bmax_off256, tmp, cmp);
            }
            ZSCORE16(i1_256, y1_256);

            {   /* narrowing: smithWaterman256_16's four scans */
                int l;
                for (l = beg; l < end; l++) {
                    __m256i t_ = _mm256_or_si256(_mm256_load_si256((__m256i *)(F + l * SIMD_WIDTH16)),
                                                 _mm256_load_si256((__m256i *)(H_h + l * SIMD_WIDTH16)));
                    if (((uint32_t) _mm256_movemask_epi8(_mm256_cmpeq_epi16(t_, zero256)) & dmask32) == dmask32) nbeg = l;
                    else break;
                }
                for (l = end; l >= beg; l--) {
                    __m256i t_ = _mm256_or_si256(_mm256_load_si256((__m256i *)(F + l * SIMD_WIDTH16)),
                                                 _mm256_load_si256((__m256i *)(H_h + l * SIMD_WIDTH16)));
                    if (((uint32_t) _mm256_movemask_epi8(_mm256_cmpeq_epi16(t_, zero256)) & dmask32) != dmask32) break;
                }
                nend = l + 2 < vncol ? l + 2 : vncol;
                __m256i tmpb = ff256;
                __m256i exit1 = _mm256_xor_si256(exit0, ff256);
                __m256i l256 = _mm256_set1_epi16(beg);
                for (l = beg; l < end; l++) {
                    __m256i t_ = _mm256_or_si256(_mm256_load_si256((__m256i *)(F + l * SIMD_WIDTH16)),
                                                 _mm256_load_si256((__m256i *)(H_h + l * SIMD_WIDTH16)));
                    t_ = _mm256_cmpeq_epi16(_mm256_or_si256(t_, exit1), zero256);
                    if (((uint32_t) _mm256_movemask_epi8(t_) & dmask32) == 0) break;
                    t_ = _mm256_and_si256(t_, tmpb);
                    l256 = _mm256_add_epi16(l256, one256);
                    head256 = _mm256_blendv_epi8(head256, l256, t_);
                    tmpb = t_;
                }
                __m256i index256 = tail256;
                tmpb = ff256;
                l256 = _mm256_set1_epi16(end);
                for (l = end; l >= beg; l--) {
                    __m256i t_ = _mm256_or_si256(_mm256_load_si256((__m256i *)(F + l * SIMD_WIDTH16)),
                                                 _mm256_load_si256((__m256i *)(H_h + l * SIMD_WIDTH16)));
                    t_ = _mm256_cmpeq_epi16(_mm256_or_si256(t_, exit1), zero256);
                    if (((uint32_t) _mm256_movemask_epi8(t_) & dmask32) == 0) break;
                    t_ = _mm256_and_si256(t_, tmpb);
                    l256 = _mm256_sub_epi16(l256, one256);
                    index256 = _mm256_blendv_epi8(index256, l256, t_);
                    tmpb = t_;
                }
                index256 = _mm256_add_epi16(index256, two256);
                tail256 = _mm256_min_epi16(index256, qlen256);
            }
            if (bsw16c_lane_mask256(exit0) == 0) break;
            }   /* rows */
            _mm256_store_si256((__m256i *) V.head, head256);
            _mm256_store_si256((__m256i *) V.tail, tail256);
            _mm256_store_si256((__m256i *) V.maxend, maxend256);
            _mm256_store_si256((__m256i *) V.exit0, exit0);
            _mm256_store_si256((__m256i *) V.maxScore, maxScore256);
            _mm256_store_si256((__m256i *) V.x, x256);
            _mm256_store_si256((__m256i *) V.y, y256);
            _mm256_store_si256((__m256i *) V.max_off, max_off256);
            _mm256_store_si256((__m256i *) V.gscore, gscore);
            _mm256_store_si256((__m256i *) V.max_ie, max_ie256);
            V.nbeg = nbeg; V.nend = nend;
            return bsw16c_lane_mask256(exit0);
        };

        /* ---- move one live lane into a dead slot of another group, between rows (next row r);
         * as move_one in bandedSWA_compact16.inc ---- */
        auto move_one = [&](Bsw16CVec256 &Sv, Bsw16CVec256 &D, int sl, int dl, int r) __attribute__((noinline)) {
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

#endif  /* AVX2 (not AVX-512BW) */
