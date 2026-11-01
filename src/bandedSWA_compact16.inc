/* bandedSWA_compact16.inc -- same-row lane compaction for the NEON 16-bit extension kernel.
 * Included from bandedSWA.cpp right after smithWaterman128_16 (aarch64 only). getScores16 runs
 * it for batches of more than one lane group when compact16_groups_ > 0 (see
 * BandedPairWiseSW::set_lane_compaction16).
 *
 * The plain wrapper runs each group of 8 pairs to its longest lane: a lane whose pair has
 * finished (exit0 cleared) idles until the whole group ends. Here a superblock of K groups runs
 * in block lockstep: every active group runs the same block of R rows, so at the end of a block
 * every lane is at the same target row. There bsw_compact_plan() retires the emptiest group
 * whenever the live lanes fit in one group fewer, moving its live lanes into dead slots of the
 * others. The 16-bit kernel keeps absolute int16 columns, so a moved lane's state is a plain
 * copy of its slot. R is a speed knob only: a lane's result does not depend on when, or
 * whether, it moves; a finished lane just waits up to R - 1 rows longer to free its slot. The
 * group's state stays in registers across a block and goes through the per-group struct only
 * at block ends, which is what keeps this driver as fast as the plain kernel with nothing to
 * compact.
 *
 * The row body is smithWaterman128_16's row; any change to one must be mirrored in the other
 * (test/unit/test_bandedswa_compact16.cpp compares the two field by field). Exactness rests on lane
 * independence (the note at the top of bandedSWA.cpp): a lane's result depends only on its own
 * slot and on the group's column range covering its band. A dead lane (exit0 clear) has frozen
 * results: maxScore / x / y / max_off are only updated under exit0, the query-end fold needs
 * exit0, and the z-drop only clears exit0. */
#if defined(__aarch64__)

#include <algorithm>
#include <cstring>

#include "bsw_compact.h"

static_assert(SIMD_WIDTH16 == 8, "bandedSWA_compact16.inc: Bsw16CVec holds 8 lanes");

namespace {
struct alignas(64) Bsw16CVec {
    int16_t head[8], tail[8], maxend[8], exit0[8], maxScore[8], x[8], y[8], max_off[8];
    int16_t gscore[8], max_ie[8], myband[8], qlen[8], mlen[8], h0[8];
    int tlenw[8];
    SeqPair *pair[8];
    int nbeg, nend, ncol, minq;
    uint16_t *seq1, *seq2;
    int16_t *H, *F;
};
struct Bsw16CBuf {
    uint8_t *mem = nullptr; size_t cap = 0;
    ~Bsw16CBuf() { if (mem) _mm_free(mem); }
};
}  // namespace
static thread_local Bsw16CBuf g_bsw16c_buf;

static inline uint32_t bsw16c_lane_mask(__m128i m)
{
    static const uint16_t bits[8] = {1, 2, 4, 8, 16, 32, 64, 128};
    const uint16x8_t b = vandq_u16(vreinterpretq_u16_m128i(m), vld1q_u16(bits));
    return (uint32_t) vaddvq_u16(b);
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
    const size_t need = (size_t) K * (sizeof(Bsw16CVec) + 4 * S * sizeof(int16_t)) + S * sizeof(int16_t) + 128;
    if (g_bsw16c_buf.cap < need) {
        if (g_bsw16c_buf.mem) _mm_free(g_bsw16c_buf.mem);
        g_bsw16c_buf.mem = (uint8_t *) _mm_malloc(need, 64);
        if (g_bsw16c_buf.mem == NULL) { fprintf(stderr, "Error! Mem not allocated!!!\n"); exit(EXIT_FAILURE); }
        g_bsw16c_buf.cap = need;
    }
    Bsw16CVec *Vs = (Bsw16CVec *) g_bsw16c_buf.mem;
    int16_t *slab = (int16_t *) (((uintptr_t) (g_bsw16c_buf.mem + (size_t) K * sizeof(Bsw16CVec)) + 63) & ~(uintptr_t) 63);
    for (int v = 0; v < K; v++) {
        Vs[v].seq1 = (uint16_t *) (slab + (size_t) (4 * v + 0) * S);
        Vs[v].seq2 = (uint16_t *) (slab + (size_t) (4 * v + 1) * S);
        Vs[v].H    = slab + (size_t) (4 * v + 2) * S;
        Vs[v].F    = slab + (size_t) (4 * v + 3) * S;
    }
    int16_t *sbt_buf = slab + (size_t) (4 * K) * S;   /* one row's scores, reused by every group */

    /* ---- constants (as smithWatermanBatchWrapper16 + smithWaterman128_16) ---- */
    const int eb = end_bonus;
    int16_t maxsc = 0;
    if (maxsc < w_match) maxsc = w_match;
    if (maxsc < w_mismatch) maxsc = w_mismatch;
    if (maxsc < w_ambig) maxsc = w_ambig;
    const int zdrop = this->zdrop;
    const uint16_t ambRef = (uint16_t) 0xFFFF, ambQer = (uint16_t) 0xFFFF;   /* NEON: SYM encoding */

    const bool gen_mat = bsw_gen_mat_;
    const BswFreedCell fc = bsw_fc_;
    const __m128i oe_ins128_s = _mm_set1_epi16(this->o_ins + this->e_ins);   /* the pack's seed row */
    const __m128i e_ins128_s  = _mm_set1_epi16(this->e_ins);
    const __m128i zero128_s   = _mm_setzero_si128();

    uint64_t live[BSW_COMPACT_GROUPS_MAX];
    uint8_t active[BSW_COMPACT_GROUPS_MAX];
    BswMove moves[BSW_COMPACT_GROUPS_MAX * SIMD_WIDTH16];

    auto writeback = [&](Bsw16CVec &V, int l) {
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
            Bsw16CVec &V = Vs[v];
            SeqPair *p = pairArray + sb0 + v * SIMD_WIDTH16;
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
                const uint8_t *s1 = seqBufRef + (int64_t) sp.idr, *s2 = seqBufQer + (int64_t) sp.idq;
                for (int k = 0; k < sp.len1; k++) V.seq1[k * SIMD_WIDTH16 + l] = (s1[k] == AMBIG ? ambRef : s1[k]);
                for (int k = sp.len1; k <= nrow; k++) V.seq1[k * SIMD_WIDTH16 + l] = DUMMY1;
                for (int k = 0; k < sp.len2; k++) V.seq2[k * SIMD_WIDTH16 + l] = (s2[k] == AMBIG ? ambQer : s2[k]);
                for (int k = sp.len2; k <= ncol; k++) V.seq2[k * SIMD_WIDTH16 + l] = DUMMY2;
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
            V.nbeg = 0; V.nend = vcol; V.ncol = vcol; V.minq = minq - 1;
            /* insertion seed row (H1) as the plain wrapper, through the superblock's last column */
            const __m128i h0_128 = _mm_load_si128((const __m128i *) V.h0);
            _mm_store_si128((__m128i *) V.H, h0_128);
            __m128i t = _mm_blendv_epi16(zero128_s, _mm_sub_epi16(h0_128, oe_ins128_s), _mm_cmpgt_epi16(h0_128, oe_ins128_s));
            _mm_store_si128((__m128i *) (V.H + SIMD_WIDTH16), t);
            for (int k = 2; k <= ncol; k++) {
                t = _mm_max_epi16(_mm_sub_epi16(t, e_ins128_s), zero128_s);
                _mm_store_si128((__m128i *) (V.H + k * SIMD_WIDTH16), t);
            }
            /* H1[maxLen2] = 0; the plain wrapper writes it first, so H1[0] / H1[1] win over it */
            if (vcol >= 2) _mm_store_si128((__m128i *) (V.H + vcol * SIMD_WIDTH16), zero128_s);
            for (int j = 0; j <= ncol; j++) _mm_store_si128((__m128i *) (V.F + j * SIMD_WIDTH16), zero128_s);
            active[v] = 1;
        };
        for (int v = 0; v < nv; v++) pack_group(v);

        /* ---- rows [i0, i1) of one group: smithWaterman128_16's row, with the group's state in
         * registers across the block (loaded once, stored once). always_inline: the hot function
         * must be the block loop with the row in it, not a call per row. Returns the live-lane
         * mask after the block (0 when the group is done: every lane finished, or a row whose
         * every lane's maximum is zero, where the plain kernel leaves the group). ---- */
        auto rows = [&](Bsw16CVec &V, const int i0, const int i1) __attribute__((always_inline)) -> uint32_t {
            /* local copies of the constants: the row stores through int16 pointers, which may
             * alias the captured vectors, so captured ones would be reloaded every cell */
            /* Constants from scalars (the compiler may rematerialize these instead of spilling);
             * the arm-specific ones are built inside their arms below. */
            const __m128i match128    = _mm_set1_epi16(this->w_match);
            const __m128i mismatch128 = _mm_set1_epi16(this->w_mismatch);
            const __m128i w_ambig_128 = _mm_set1_epi16(this->w_ambig);
            const __m128i e_del128    = _mm_set1_epi16(this->e_del);
            const __m128i oe_del128   = _mm_set1_epi16(this->o_del + this->e_del);
            const __m128i e_ins128    = _mm_set1_epi16(this->e_ins);
            const __m128i oe_ins128   = _mm_set1_epi16(this->o_ins + this->e_ins);
            const __m128i zero128 = _mm_setzero_si128();
            const __m128i one128  = _mm_set1_epi16(1);
            const __m128i two128  = _mm_set1_epi16(2);
            const __m128i ff128   = _mm_set1_epi16(0xFFFF);
            int16_t *const H_h = V.H, *const F = V.F;
            int16_t *const sbt_row = sbt_buf;
            const int w_ = w;
            const uint16_t *seq2SoA = V.seq2;
            const int minq = V.minq, vncol = V.ncol;
            const __m128i zdrop128 = _mm_set1_epi16(zdrop);
            __m128i maxScore128 = _mm_load_si128((const __m128i *) V.maxScore);
            __m128i x128        = _mm_load_si128((const __m128i *) V.x);
            __m128i y128        = _mm_load_si128((const __m128i *) V.y);
            __m128i max_off128  = _mm_load_si128((const __m128i *) V.max_off);
            __m128i gscore      = _mm_load_si128((const __m128i *) V.gscore);
            __m128i max_ie128   = _mm_load_si128((const __m128i *) V.max_ie);
            int nbeg_ = V.nbeg, nend_ = V.nend;
            __m128i head128     = _mm_load_si128((const __m128i *) V.head);
            __m128i tail128     = _mm_load_si128((const __m128i *) V.tail);
            __m128i maxend128   = _mm_load_si128((const __m128i *) V.maxend);
            __m128i exit0       = _mm_load_si128((const __m128i *) V.exit0);
            const __m128i myband128 = _mm_load_si128((const __m128i *) V.myband);
            const __m128i qlen128   = _mm_load_si128((const __m128i *) V.qlen);
            const __m128i mlen128   = _mm_load_si128((const __m128i *) V.mlen);
            const __m128i h0_128    = _mm_load_si128((const __m128i *) V.h0);
            int j;
            for (int i = i0; i < i1; i++) {
            __m128i e11 = zero128;
            __m128i h00, h11, h10;
            __m128i s10 = _mm_load_si128((const __m128i *) (V.seq1 + i * SIMD_WIDTH16));
            int beg = nbeg_, end = nend_;
            if (beg < i - w) beg = i - w;
            if (end > i + w + 1) end = i + w + 1;
            if (end > vncol) end = vncol;
            h10 = zero128;
            if (beg == 0) {   /* H2[i + 1]: the deletion seed column, h0 - o_del - (i+1) e_del in int16, >= 0 */
                const int16_t sub = (int16_t) (this->o_del + (i + 1) * this->e_del);
                h10 = _mm_max_epi16(_mm_sub_epi16(h0_128, _mm_set1_epi16(sub)), zero128);
            }
            __m128i j128 = zero128;
            __m128i maxRS1 = zero128;
            __m128i i1_128 = _mm_set1_epi16(i + 1);
            __m128i y1_128 = zero128;

            __m128i phead128 = head128, ptail128 = tail128;
            head128 = _mm_max_epi16(head128, _mm_sub_epi16(_mm_set1_epi16(i), myband128));
            tail128 = _mm_min_epu16(tail128, _mm_add_epi16(i1_128, myband128));
            tail128 = _mm_min_epu16(tail128, qlen128);
            {   __m128i need_ = _mm_and_si128(_mm_cmpeq_epi16(_mm_sub_epi16(tail128, maxend128), two128), exit0);
                if (any_lane_set8(need_)) {
                    int16_t need_a[SIMD_WIDTH16] __attribute((aligned(16)));
                    int16_t me_a[SIMD_WIDTH16]   __attribute((aligned(16)));
                    _mm_store_si128((__m128i *) need_a, need_);
                    _mm_store_si128((__m128i *) me_a, maxend128);
                    const int oe_ = this->o_ins + this->e_ins;
                    for (int l = 0; l < SIMD_WIDTH16; l++) if (need_a[l]) {
                        const int c = me_a[l] + 1;
                        H_h[c * SIMD_WIDTH16 + l] = (int16_t) bsw_seed_row((int16_t) V.h0[l], oe_, this->e_ins, c);
                        F[c * SIMD_WIDTH16 + l] = 0;
                    }
                }
                maxend128 = _mm_max_epi16(maxend128, tail128);
            }
            {   /* trim, counted */
                __m128i cmph = _mm_and_si128(_mm_cmpeq_epi16(head128, phead128), _mm_cmpeq_epi16(tail128, ptail128));
                int trim_hi = all_lanes_set8(cmph) ? beg : hmax_epi16_128(head128);
                if (trim_hi > end) trim_hi = end;
                for (int l = beg; l < trim_hi; l++) {
                    __m128i h128 = _mm_load_si128((__m128i *)(H_h + l * SIMD_WIDTH16));
                    __m128i f128 = _mm_load_si128((__m128i *)(F + l * SIMD_WIDTH16));
                    __m128i cmp1 = _mm_or_si128(_mm_cmpgt_epi16(head128, _mm_set1_epi16(l)),
                                                _mm_cmpgt_epi16(_mm_set1_epi16(l + 1), tail128));
                    _mm_store_si128((__m128i *)(F + l * SIMD_WIDTH16), _mm_andnot_si128(cmp1, f128));
                    _mm_store_si128((__m128i *)(H_h + l * SIMD_WIDTH16), _mm_andnot_si128(cmp1, h128));
                }
            }
            {
                __m128i cmpim = _mm_cmpgt_epi16(i1_128, mlen128);
                cmpim = _mm_or_si128(cmpim, _mm_cmpeq_epi16(tail128, head128));
                cmpim = _mm_or_si128(cmpim, _mm_cmpgt_epi16(head128, tail128));
                exit0 = _mm_blendv_epi16(exit0, zero128, cmpim);
            }
            /* the score pre-pass for [beg, end) */
            if (CPG) {
                const __m128i cpgref128 = _mm_set1_epi16(this->cpg_ref_code_);
                const __m128i cpgread128 = _mm_set1_epi16(this->cpg_read_);
                const __m128i cpgval128 = _mm_set1_epi16(this->cpg_val_);
                __m128i rowcpg;
                BSW_CPG_ROW16_128(s10, rowcpg, cpgref128);
                for (int jp = beg; jp < end; jp++) {
                    __m128i s2 = _mm_load_si128((const __m128i *)(seq2SoA + jp * SIMD_WIDTH16));
                    __m128i sbt11;
                    SBT_PREPASS16_SYM(s10, s2, sbt11, mismatch128, match128, w_ambig_128, ff128);
                    BSW_CPG_OVERRIDE16_128(s2, rowcpg, sbt11, cpgval128, cpgread128);
                    _mm_store_si128((__m128i *)(sbt_row + jp * SIMD_WIDTH16), sbt11);
                }
            } else if (!gen_mat) {
                for (int jp = beg; jp < end; jp++) {
                    __m128i s2 = _mm_load_si128((const __m128i *)(seq2SoA + jp * SIMD_WIDTH16));
                    __m128i sbt11;
                    SBT_PREPASS16_SYM(s10, s2, sbt11, mismatch128, match128, w_ambig_128, ff128);
                    _mm_store_si128((__m128i *)(sbt_row + jp * SIMD_WIDTH16), sbt11);
                }
            } else if (fc.rank1) {
                const __m128i frref128 = _mm_set1_epi16(fc.ref), frread128 = _mm_set1_epi16(fc.read);
                __m128i alt10 = _mm_blendv_epi16(s10, frread128, _mm_cmpeq_epi16(s10, frref128));
                for (int jp = beg; jp < end; jp++) {
                    __m128i s2 = _mm_load_si128((const __m128i *)(seq2SoA + jp * SIMD_WIDTH16));
                    __m128i sbt11;
                    SBT_PREPASS16_RANK1(s10, s2, alt10, sbt11, mismatch128, match128, w_ambig_128, ff128);
                    _mm_store_si128((__m128i *)(sbt_row + jp * SIMD_WIDTH16), sbt11);
                }
            } else {
                const __m128i amat128 = _mm_load_si128((const __m128i *)bsw_amat_bytes_);
                const __m128i three128 = _mm_set1_epi16(3);
                for (int jp = beg; jp < end; jp++) {
                    __m128i s2 = _mm_load_si128((const __m128i *)(seq2SoA + jp * SIMD_WIDTH16));
                    __m128i sbt11;
                    SBT_PREPASS16_AMAT(s10, s2, sbt11, amat128, mismatch128, w_ambig_128, ff128, three128);
                    _mm_store_si128((__m128i *)(sbt_row + jp * SIMD_WIDTH16), sbt11);
                }
            }
#define C16_CELL(QCAP) \
            { \
                __m128i f11, f21; \
                __m128i sbt11 = _mm_load_si128((__m128i *)(sbt_row + j * SIMD_WIDTH16)); \
                h00 = _mm_load_si128((__m128i *)(H_h + j * SIMD_WIDTH16)); \
                f11 = _mm_load_si128((__m128i *)(F + j * SIMD_WIDTH16)); \
                __m128i pj128 = j128; \
                j128 = _mm_add_epi16(j128, one128); \
                MAIN_CODE16_CORE(sbt11, h00, h11, e11, f11, f21, zero128, \
                                 e_ins128, oe_ins128, e_del128, oe_del128); \
                __m128i cmp1 = _mm_or_si128(_mm_cmpgt_epi16(head128, pj128), _mm_cmpgt_epi16(pj128, tail128)); \
                h10 = _mm_andnot_si128(cmp1, h10); \
                __m128i bmaxRS = maxRS1; \
                maxRS1 = _mm_max_epi16(maxRS1, h11); \
                __m128i cmpA = _mm_cmpeq_epi16(maxRS1, h11); \
                cmp1 = _mm_cmpgt_epi16(j128, tail128); \
                f21 = _mm_andnot_si128(cmp1, f21); \
                cmpA = _mm_blendv_epi16(y1_128, j128, cmpA); \
                y1_128 = _mm_blendv_epi16(cmpA, y1_128, cmp1); \
                maxRS1 = _mm_blendv_epi16(maxRS1, bmaxRS, cmp1); \
                _mm_store_si128((__m128i *)(F + j * SIMD_WIDTH16), f21); \
                _mm_store_si128((__m128i *)(H_h + j * SIMD_WIDTH16), h10); \
                h10 = h11; \
                if (QCAP) hqe128 = _mm_blendv_epi16(hqe128, h11, _mm_cmpeq_epi16(j128, qlen128)); \
            }
            j128 = _mm_set1_epi16(beg);
            __m128i hqe128 = zero128;
            {
                const int jq_ = minq < beg ? beg : (minq > end ? end : minq);
                for (j = beg; j < jq_; j++) C16_CELL(0)
                for (j = jq_; j < end; j++) C16_CELL(1)
            }
#undef C16_CELL
            {
                __m128i cmp1 = _mm_or_si128(_mm_cmpgt_epi16(head128, j128), _mm_cmpgt_epi16(j128, tail128));
                h10 = _mm_andnot_si128(cmp1, h10);
                _mm_store_si128((__m128i *)(H_h + j * SIMD_WIDTH16), h10);
                _mm_store_si128((__m128i *)(F + j * SIMD_WIDTH16), zero128);
            }
            {   /* query-end fold, before the zero-row exit (see smithWaterman128_16) */
                __m128i qin = _mm_andnot_si128(_mm_cmpgt_epi16(qlen128, _mm_set1_epi16(end)),
                                               _mm_cmpgt_epi16(qlen128, _mm_set1_epi16(beg)));
                __m128i upd = _mm_andnot_si128(_mm_cmpgt_epi16(qlen128, tail128), _mm_and_si128(qin, exit0));
                __m128i cmp_gh = _mm_cmpgt_epi16(gscore, hqe128);
                __m128i cand_ie = _mm_blendv_epi16(i1_128, max_ie128, cmp_gh);
                max_ie128 = _mm_blendv_epi16(max_ie128, cand_ie, upd);
                gscore = _mm_blendv_epi16(gscore, _mm_max_epi16(gscore, hqe128), upd);
            }

            __m128i bmaxScore128 = maxScore128;
            __m128i tmp = _mm_cmpeq_epi16(maxRS1, zero128);
            if (vmaxvq_u16(vreinterpretq_u16_m128i(maxRS1)) == 0) {
                /* every lane's row max is zero: the plain kernel leaves the group here; every
                 * lane's results are final */
                exit0 = zero128;
                break;
            }
            exit0 = _mm_blendv_epi16(exit0, zero128, tmp);
            __m128i score128 = _mm_max_epi16(maxScore128, maxRS1);
            maxScore128 = _mm_blendv_epi16(maxScore128, score128, exit0);
            __m128i cmp = _mm_cmpgt_epi16(maxScore128, bmaxScore128);
            y128 = _mm_blendv_epi16(y128, y1_128, cmp);
            x128 = _mm_blendv_epi16(x128, i1_128, cmp);
            {
                __m128i ab = _mm_subs_epu16(y1_128, i1_128);
                __m128i ba = _mm_subs_epu16(i1_128, y1_128);
                tmp = _mm_or_si128(ab, ba);
                __m128i bmax_off128 = max_off128;
                tmp = _mm_max_epi16(max_off128, tmp);
                max_off128 = _mm_blendv_epi16(bmax_off128, tmp, cmp);
            }
            ZSCORE16(i1_128, y1_128);

            {   /* narrowing: smithWaterman128_16's four scans */
                int l;
                int nbeg = nbeg_;
                for (l = beg; l < end; l++) {
                    __m128i f128 = _mm_load_si128((__m128i *)(F + l * SIMD_WIDTH16));
                    __m128i h128 = _mm_load_si128((__m128i *)(H_h + l * SIMD_WIDTH16));
                    if (vmaxvq_u16(vorrq_u16(vreinterpretq_u16_m128i(f128), vreinterpretq_u16_m128i(h128))) == 0) nbeg = l;
                    else break;
                }
                nbeg_ = nbeg;
                for (l = end; l >= beg; l--) {
                    __m128i f128 = _mm_load_si128((__m128i *)(F + l * SIMD_WIDTH16));
                    __m128i h128 = _mm_load_si128((__m128i *)(H_h + l * SIMD_WIDTH16));
                    if (vmaxvq_u16(vorrq_u16(vreinterpretq_u16_m128i(f128), vreinterpretq_u16_m128i(h128))) != 0) break;
                }
                nend_ = l + 2 < vncol ? l + 2 : vncol;
                __m128i tmpb = ff128;
                __m128i exit1 = _mm_xor_si128(exit0, ff128);
                __m128i l128 = _mm_set1_epi16(beg);
                for (l = beg; l < end; l++) {
                    __m128i f128 = _mm_load_si128((__m128i *)(F + l * SIMD_WIDTH16));
                    __m128i h128 = _mm_load_si128((__m128i *)(H_h + l * SIMD_WIDTH16));
                    __m128i t2 = _mm_cmpeq_epi16(_mm_or_si128(_mm_or_si128(f128, h128), exit1), zero128);
                    if (vmaxvq_u8(vreinterpretq_u8_m128i(t2)) == 0) break;
                    t2 = _mm_and_si128(t2, tmpb);
                    l128 = _mm_add_epi16(l128, one128);
                    head128 = _mm_blendv_epi16(head128, l128, t2);
                    tmpb = t2;
                }
                __m128i index128 = tail128;
                tmpb = ff128;
                l128 = _mm_set1_epi16(end);
                for (l = end; l >= beg; l--) {
                    __m128i f128 = _mm_load_si128((__m128i *)(F + l * SIMD_WIDTH16));
                    __m128i h128 = _mm_load_si128((__m128i *)(H_h + l * SIMD_WIDTH16));
                    __m128i t2 = _mm_cmpeq_epi16(_mm_or_si128(_mm_or_si128(f128, h128), exit1), zero128);
                    if (vmaxvq_u8(vreinterpretq_u8_m128i(t2)) == 0) break;
                    t2 = _mm_and_si128(t2, tmpb);
                    l128 = _mm_sub_epi16(l128, one128);
                    index128 = blendv_fullmask8(index128, l128, t2);
                    tmpb = t2;
                }
                index128 = _mm_add_epi16(index128, two128);
                tail128 = _mm_min_epi16(index128, qlen128);
            }
            if (bsw16c_lane_mask(exit0) == 0) break;
            }   /* rows */
            _mm_store_si128((__m128i *) V.head, head128);
            _mm_store_si128((__m128i *) V.tail, tail128);
            _mm_store_si128((__m128i *) V.maxend, maxend128);
            _mm_store_si128((__m128i *) V.exit0, exit0);
            _mm_store_si128((__m128i *) V.maxScore, maxScore128);
            _mm_store_si128((__m128i *) V.x, x128);
            _mm_store_si128((__m128i *) V.y, y128);
            _mm_store_si128((__m128i *) V.max_off, max_off128);
            _mm_store_si128((__m128i *) V.gscore, gscore);
            _mm_store_si128((__m128i *) V.max_ie, max_ie128);
            V.nbeg = nbeg_; V.nend = nend_;
            return bsw16c_lane_mask(exit0);
        };

        /* ---- move one live lane into a dead slot of another group, between rows (next row r).
         * The lane takes its slot: the per-lane trackers; its DP row (H, F) and query codes from
         * its head to the superblock's last column; its target rows r .. tlen-1. Left of the
         * head the kernel needs zero H/F wherever the group range reaches (the trim/store
         * invariant), so the destination slot is cleared from the group's lowest future column
         * up to the head. The copies are per lane and scalar: moves are rare (a few per superblock). */
        auto move_one = [&](Bsw16CVec &Sv, Bsw16CVec &D, int sl, int dl, int r) __attribute__((noinline)) {
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

        /* ---- row lockstep over the superblock ---- */
        /* Block lockstep: every active group runs rows [i0, i1), then the plan and the moves
         * run with every live lane at row i1. R is a speed knob only: a lane's result does not
         * depend on when it moves (lane independence), only finished lanes wait longer. */
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

#endif  /* __aarch64__ */
