/* bandedSWA_compact.inc -- same-row lane compaction for the NEON 8-bit extension kernel.
 * Included from bandedSWA.cpp right after smithWaterman128_8 (aarch64 only), so it runs that
 * kernel's cell loops and gscore fold through the same macros (EXT13_RUN_SPLIT8_128,
 * EXT13_SBT8_*, DP_CELL_BODY8_128[_FAST], BSW8_GSCORE_NARROW). getScores8 runs it for batches of
 * more than one group at bands w >= compact_min_w_ when compact_groups_ > 0 (see
 * BandedPairWiseSW::set_lane_compaction).
 *
 * The plain wrapper runs each group of 16 pairs to its longest lane: a lane whose pair has
 * finished idles until the whole group ends. Here a superblock of K groups runs in row lockstep
 * (every lane at the same target row i). After a row in which lanes finished, bsw_compact_plan()
 * retires the emptiest group whenever the live lanes fit in one group fewer, moving its live
 * lanes into dead slots of the others. Because every lane is at row i, a lane keeps the
 * absolute-column frame: moving it is a copy of its slot (DP row from its head to ncol, query and
 * remaining target bytes, the per-lane trackers). The row body is the lean (BSW8_ROW_LEAN) row
 * of smithWaterman128_8 with its persistent vectors kept in a per-group struct.
 *
 * Exactness: a lane's result depends only on its own slot and on the group's column range
 * covering its band (lane independence; see the note at the top of bandedSWA.cpp). A move copies
 * the whole slot right of the lane's head, widens the receiving group's column range to cover the
 * moved lane ([nbeg, nend] union with the donor group's, ncol up to and minq down to the lane's
 * query length), and clears the receiving slot left of the head (see move_lane), so every pair's
 * six result fields equal the plain kernel's. That holds for every input the 8-bit kernel is
 * defined on: pairs inside the routing envelope, w <= BSW8_MAX_W (124), gap extends >= 1.
 *
 * Any change to the row body here must be mirrored in smithWaterman128_8's BSW8_ROW_LEAN path
 * (and vice versa); test/unit/test_bandedswa_compact.cpp compares the two field by field. */
#if defined(__aarch64__)

#include "bsw_compact.h"

static_assert(SIMD_WIDTH8 == 16, "bandedSWA_compact.inc: Bsw8CVec and bsw8c_lane_mask hold 16 lanes");

namespace {
struct alignas(64) Bsw8CVec {
    uint8_t head[16], tail[16], maxend[16], exit0[16], maxScore[16], y[16], max_off[16], h0[16], myband[16];
    uint8_t gbest8[16], gset8[16];
    int16_t xrow16[16], ierow16[16];
    int32_t qlenw[16], mlenw[16], tlenw[16], best_abs[16];
    SeqPair *pair[16];
    int nbeg, nend, ncol, minq;
    uint8_t *seq1, *seq2;
    int8_t *H, *F;
};
struct Bsw8CBuf {
    uint8_t *mem = nullptr; size_t cap = 0;
    ~Bsw8CBuf() { if (mem) _mm_free(mem); }
};
}  // namespace
/* Per-thread scratch, grown on demand and kept until thread exit: the pipeline builds a new
 * kernel object for every batch (mem_chain2aln_across_reads_V2), so a per-object buffer would
 * be allocated and faulted in again each time. */
static thread_local Bsw8CBuf g_bsw8c_buf;

static inline uint32_t bsw8c_lane_mask(__m128i m)
{
    static const uint8_t bits[16] = {1, 2, 4, 8, 16, 32, 64, 128, 1, 2, 4, 8, 16, 32, 64, 128};
    const uint8x16_t b = vandq_u8(vreinterpretq_u8_m128i(m), vld1q_u8(bits));
    return (uint32_t) vaddv_u8(vget_low_u8(b)) | ((uint32_t) vaddv_u8(vget_high_u8(b)) << 8);
}

void BandedPairWiseSW::smithWatermanBatchWrapper8Compact(SeqPair *pairArray,
                                                         uint8_t *seqBufRef,
                                                         uint8_t *seqBufQer,
                                                         int32_t numPairs,
                                                         int32_t w,
                                                         int K)
{
    if (K < 1) K = 1;
    if (K > BSW_COMPACT_GROUPS_MAX) K = BSW_COMPACT_GROUPS_MAX;
    const int32_t roundNumPairs = ((numPairs + SIMD_WIDTH8 - 1) / SIMD_WIDTH8) * SIMD_WIDTH8;
    for (int32_t ii = numPairs; ii < roundNumPairs; ii++) {   // padding lanes, as the plain wrapper
        pairArray[ii].id = ii; pairArray[ii].len1 = 0; pairArray[ii].len2 = 0;
        pairArray[ii].idr = 0; pairArray[ii].idq = 0; pairArray[ii].h0 = 0;
    }
    const size_t S = (size_t) MAX_SEQ_LEN8 * SIMD_WIDTH8;
    const size_t need = (size_t) K * (sizeof(Bsw8CVec) + 4 * S) + 128;
    if (g_bsw8c_buf.cap < need) {
        if (g_bsw8c_buf.mem) _mm_free(g_bsw8c_buf.mem);
        g_bsw8c_buf.mem = (uint8_t *) _mm_malloc(need, 64);
        if (g_bsw8c_buf.mem == NULL) { fprintf(stderr, "Error! Mem not allocated!!!\n"); exit(EXIT_FAILURE); }
        g_bsw8c_buf.cap = need;
    }
    Bsw8CVec *Vs = (Bsw8CVec *) g_bsw8c_buf.mem;
    uint8_t *slab = g_bsw8c_buf.mem + (size_t) K * sizeof(Bsw8CVec);
    slab = (uint8_t *) (((uintptr_t) slab + 63) & ~(uintptr_t) 63);
    for (int v = 0; v < K; v++) {
        Vs[v].seq1 = slab + (size_t) (4 * v + 0) * S;
        Vs[v].seq2 = slab + (size_t) (4 * v + 1) * S;
        Vs[v].H    = (int8_t *) (slab + (size_t) (4 * v + 2) * S);
        Vs[v].F    = (int8_t *) (slab + (size_t) (4 * v + 3) * S);
    }

    /* ---- constants (as smithWatermanBatchWrapper8 + smithWaterman128_8) ---- */
    const int eb = end_bonus;
    int8_t maxsc = 0;
    if (maxsc < w_match) maxsc = w_match;
    if (maxsc < w_mismatch) maxsc = w_mismatch;
    if (maxsc < w_ambig) maxsc = w_ambig;
    assert(this->zdrop >= 0 && "8-bit banded SW: negative zdrop");
    const int zdrop = this->zdrop;
#ifdef BSW8_ASSERT_ENVELOPE
    const int BYTE_CEIL = 255 - std::max(1, std::max((int) this->w_match, (int) this->w_ambig));
#endif

    __m128i match128     = _mm_set1_epi8(this->w_match);
    __m128i w_ambig_128  = _mm_set1_epi8(this->w_ambig);
    const bool gen_mat = bsw_gen_mat_;
    const BswFreedCell fc = bsw_fc_;
    __m128i frref128  = _mm_set1_epi8(fc.ref);
    __m128i frread128 = _mm_set1_epi8(fc.read);
    __m128i pmat128 = _mm_load_si128((__m128i *)bsw_pmat_bytes_);
    __m128i pmat_pos128, pmat_neg128;
    SBT_SPLIT8(pmat128, pmat_pos128, pmat_neg128, _mm_setzero_si128());
    __m128i amat128 = _mm_load_si128((__m128i *)bsw_amat_bytes_);
    __m128i three128 = _mm_set1_epi8(3);
    __m128i e_del128    = _mm_set1_epi8(this->e_del);
    __m128i oe_del128   = _mm_set1_epi8(this->o_del + this->e_del);
    __m128i e_ins128    = _mm_set1_epi8(this->e_ins);
    __m128i oe_ins128   = _mm_set1_epi8(this->o_ins + this->e_ins);
    const __m128i zero128 = _mm_setzero_si128();
    const __m128i one128  = _mm_set1_epi8(1);
    const __m128i two128  = _mm_set1_epi8(2);
    const __m128i ff128   = _mm_set1_epi8(0xFF);
    /* column -1 deletion seed H_v[i+1] = subs_epu8 chain h0 - o_del - (i+1)*e_del, with the
     * byte-cast o_del / e_del the plain wrapper's set1_epi8 uses: closed form of the chain. */
    const int o_del8 = (uint8_t) (int8_t) this->o_del, e_del8 = (uint8_t) (int8_t) this->e_del;
    const bool zgate_exact = this->e_del >= 0 && this->e_del <= INT32_MAX / (4 * MAX_SEQ_LEN8) &&
                             this->e_ins >= 0 && this->e_ins <= INT32_MAX / (4 * MAX_SEQ_LEN8);
    const __m128i zgate128 = _mm_set1_epi8((int8_t) (zdrop < 255 ? zdrop : 255));
    /* band-edge write: byte penalties, as the seed row's set1_epi8 / subs_epu8 uses them */
    const int oe_ins8 = (uint8_t) (this->o_ins + this->e_ins), e_ins8 = (uint8_t) this->e_ins;

    uint64_t live[BSW_COMPACT_GROUPS_MAX];
    uint8_t active[BSW_COMPACT_GROUPS_MAX];
    BswMove moves[BSW_COMPACT_GROUPS_MAX * SIMD_WIDTH8];

    auto writeback = [&](Bsw8CVec &V, int l) {
        SeqPair *p = V.pair[l];
        if (!p) return;
        const int xr = V.xrow16[l];
        int sc = V.best_abs[l];
        if (xr != 0 && (int) V.maxScore[l] > sc) sc = V.maxScore[l];
        p->score   = sc;
        p->tle     = xr;
        p->qle     = (xr == 0) ? 0 : ((int) (int8_t) V.y[l] + xr - 1);
        p->max_off = (uint8_t) V.max_off[l];
        p->gscore  = V.gset8[l] ? (int) V.gbest8[l] : -1;
        p->gtle    = V.ierow16[l];
        V.pair[l] = nullptr;
    };

    for (int32_t sb0 = 0; sb0 < roundNumPairs; sb0 += K * SIMD_WIDTH8) {
        const int nv = std::min(K, (int) ((roundNumPairs - sb0) / SIMD_WIDTH8));
        int nrow = 0, ncol = 0;
        compact_counts_.superblocks++;
        for (int q = sb0; q < sb0 + nv * SIMD_WIDTH8; q++) {
            const SeqPair &sp = pairArray[q];
            if (sp.len1 < 0 || sp.len1 >= MAX_SEQ_LEN8)
                err_fatal(__func__, "bandedSWA: ref window length (len1) %d for pair %d is out of range [0, MAX_SEQ_LEN8=%d)", sp.len1, sp.id, MAX_SEQ_LEN8);
            if (sp.len2 < 0 || sp.len2 >= MAX_SEQ_LEN8)
                err_fatal(__func__, "bandedSWA: query length (len2) %d for pair %d is out of range [0, MAX_SEQ_LEN8=%d)", sp.len2, sp.id, MAX_SEQ_LEN8);
            if (nrow < sp.len1) nrow = sp.len1;
            if (ncol < sp.len2) ncol = sp.len2;
        }

        /* ---- pack and seed every group (the plain wrapper's setup, at superblock extents) ---- */
        for (int v = 0; v < nv; v++) {
            Bsw8CVec &V = Vs[v];
            SeqPair *p = pairArray + sb0 + v * SIMD_WIDTH8;
            const uint8_t *seq1p[SIMD_WIDTH8], *seq2p[SIMD_WIDTH8];
            int len1a[SIMD_WIDTH8], len2a[SIMD_WIDTH8];
            int vcol = 0, minq = 10000000;
            for (int l = 0; l < SIMD_WIDTH8; l++) {
                int h0p = p[l].h0; if (h0p < 0) h0p = 0; if (h0p > 255) h0p = 255;
                V.h0[l] = (uint8_t) h0p;
                seq1p[l] = seqBufRef + (int64_t) p[l].idr; len1a[l] = p[l].len1;
                seq2p[l] = seqBufQer + (int64_t) p[l].idq; len2a[l] = p[l].len2;
                if (vcol < p[l].len2) vcol = p[l].len2;
                if (p[l].len2 < minq) minq = p[l].len2;
                const int reach = p[l].len2 * maxsc + eb;
                int max_ins = (int) ((double) (reach - o_ins) / e_ins + 1.0); if (max_ins < 1) max_ins = 1;
                int max_del = (int) ((double) (reach - o_del) / e_del + 1.0); if (max_del < 1) max_del = 1;
                int band = w; if (max_ins < band) band = max_ins; if (max_del < band) band = max_del;
                V.myband[l] = (uint8_t) band;
                V.tlenw[l] = p[l].len1; V.qlenw[l] = p[l].len2;
                int ml = p[l].len2 + V.myband[l]; if (ml > p[l].len1) ml = p[l].len1;
                V.mlenw[l] = ml;
                V.best_abs[l] = p[l].h0;
                V.pair[l] = &p[l];
                V.head[l] = 0; V.tail[l] = 127; V.maxend[l] = (uint8_t) -128; V.exit0[l] = 0xFF;
                V.maxScore[l] = V.h0[l]; V.y[l] = 0; V.max_off[l] = 0;
                V.gbest8[l] = 0; V.gset8[l] = 0; V.xrow16[l] = 0; V.ierow16[l] = 0;
            }
            V.nbeg = 0; V.nend = vcol; V.ncol = vcol; V.minq = minq - 1;
            neon_soa_pack16(V.seq1, seq1p, len1a, len1a, nrow + 1, DUMMY1, DUMMY1, false);
            neon_soa_pack16(V.seq2, seq2p, len2a, len2a, ncol + 1, DUMMY2, DUMMY2, true);
            const __m128i h0_128 = _mm_load_si128((const __m128i *) V.h0);
            _mm_store_si128((__m128i *) (V.H + ncol * SIMD_WIDTH8), zero128);
            _mm_store_si128((__m128i *) V.H, h0_128);
            __m128i t = _mm_subs_epu8(h0_128, oe_ins128);
            _mm_store_si128((__m128i *) (V.H + SIMD_WIDTH8), t);
            for (int k = 2; k < ncol; k++) {
                t = _mm_subs_epu8(t, e_ins128);
                _mm_store_si128((__m128i *) (V.H + k * SIMD_WIDTH8), t);
            }
            for (int j = 0; j <= ncol; j++) _mm_store_si128((__m128i *) (V.F + j * SIMD_WIDTH8), zero128);
            active[v] = 1;
        }

        /* ---- one row of one group: the lean row of smithWaterman128_8 ---- */
        auto row = [&](Bsw8CVec &V, const int i) -> uint32_t {
            int8_t *H_h = V.H, *F = V.F;
            const uint8_t *seq2SoA = V.seq2;
            const int minq = V.minq, vncol = V.ncol;
            __m128i head128     = _mm_load_si128((const __m128i *) V.head);
            __m128i tail128     = _mm_load_si128((const __m128i *) V.tail);
            __m128i maxend128   = _mm_load_si128((const __m128i *) V.maxend);
            __m128i exit0       = _mm_load_si128((const __m128i *) V.exit0);
            __m128i maxScore128 = _mm_load_si128((const __m128i *) V.maxScore);
            __m128i y128        = _mm_load_si128((const __m128i *) V.y);
            __m128i max_off128  = _mm_load_si128((const __m128i *) V.max_off);
            const __m128i myband128 = _mm_load_si128((const __m128i *) V.myband);
            const __m128i h0_128    = _mm_load_si128((const __m128i *) V.h0);
            int j;

            __m128i e11 = zero128;
            __m128i h00, h11, h10;
            __m128i s10 = _mm_load_si128((const __m128i *) (V.seq1 + i * SIMD_WIDTH8));
            int beg = V.nbeg, end = V.nend;
            if (beg < i - w) beg = i - w;
            if (end > i + w + 1) end = i + w + 1;
            if (end > vncol) end = vncol;
            h10 = zero128;
            if (beg == 0) {
                long sub = (long) o_del8 + (long) (i + 1) * e_del8;
                if (sub > 255) sub = 255;
                h10 = _mm_subs_epu8(h0_128, _mm_set1_epi8((char) (uint8_t) sub));
            }
            __m128i j128 = zero128;
            __m128i maxRS1 = zero128;
            __m128i y1_128 = zero128;
            __m128i hqe128 = zero128, qfire128 = zero128;

            int8_t qlen_off_a[SIMD_WIDTH8]   __attribute((aligned(16)));
            int8_t qlen_valid_a[SIMD_WIDTH8] __attribute((aligned(16)));
            int8_t cmpim_a[SIMD_WIDTH8]      __attribute((aligned(16)));
            for (int l = 0; l < SIMD_WIDTH8; l++) {
                int qoff = V.qlenw[l] - i;
                int qoff_sat = qoff;
                if (qoff_sat >  127) qoff_sat =  127;
                if (qoff_sat < -128) qoff_sat = -128;
                qlen_off_a[l]   = (int8_t) qoff_sat;
                qlen_valid_a[l] = (qoff >= -w && qoff <= w + 1) ? (int8_t) 0xFF : 0;
                cmpim_a[l]      = ((i + 1) > V.mlenw[l]) ? (int8_t) 0xFF : 0;
            }
            __m128i qlen_off128   = _mm_load_si128((__m128i *) qlen_off_a);
            __m128i qlen_valid128 = _mm_load_si128((__m128i *) qlen_valid_a);

            __m128i phead128 = head128, ptail128 = tail128;
            head128 = _mm_max_epi8(head128, _mm_sub_epi8(zero128, myband128));
            tail128 = _mm_min_epi8(tail128, _mm_add_epi8(myband128, one128));
            tail128 = _mm_min_epi8(tail128, qlen_off128);
            {   /* band-edge write, as in smithWaterman128_8 */
                __m128i need_ = _mm_and_si128(_mm_cmpeq_epi8(_mm_sub_epi8(tail128, maxend128), two128), exit0);
                if (any_lane_set8(need_)) {
                    int8_t need_a[SIMD_WIDTH8] __attribute((aligned(16)));
                    int8_t me_a[SIMD_WIDTH8]   __attribute((aligned(16)));
                    _mm_store_si128((__m128i *) need_a, need_);
                    _mm_store_si128((__m128i *) me_a, maxend128);
                    for (int l = 0; l < SIMD_WIDTH8; l++) if (need_a[l]) {
                        const int c = i + me_a[l] + 1;
                        H_h[c * SIMD_WIDTH8 + l] = (int8_t) (uint8_t) bsw_seed_row(V.h0[l], oe_ins8, e_ins8, c);
                        F[c * SIMD_WIDTH8 + l] = 0;
                    }
                }
                maxend128 = _mm_max_epi8(maxend128, tail128);
            }
            {   /* trim, counted (BSW8_ROW_LEAN) */
                __m128i cmph = _mm_and_si128(_mm_cmpeq_epi8(head128, phead128), _mm_cmpeq_epi8(tail128, ptail128));
                if (!all_lanes_set8(cmph)) {
                    int trim_hi = i + hmax_epi8(head128);
                    if (trim_hi > end) trim_hi = end;
                    for (int l = beg; l < trim_hi; l++) {
                        __m128i h128 = _mm_load_si128((__m128i *)(H_h + l * SIMD_WIDTH8));
                        __m128i f128 = _mm_load_si128((__m128i *)(F + l * SIMD_WIDTH8));
                        __m128i pj128 = _mm_set1_epi8(l - i);
                        __m128i cmp1 = _mm_or_si128(_mm_cmpgt_epi8(head128, pj128), _mm_cmpgt_epi8(pj128, tail128));
                        h128 = blendv_fullmask8(h128, zero128, cmp1);
                        f128 = blendv_fullmask8(f128, zero128, cmp1);
                        _mm_store_si128((__m128i *)(F + l * SIMD_WIDTH8), f128);
                        _mm_store_si128((__m128i *)(H_h + l * SIMD_WIDTH8), h128);
                    }
                }
            }
            {
                __m128i cmpim = _mm_load_si128((__m128i *) cmpim_a);
                cmpim = _mm_or_si128(cmpim, _mm_cmpeq_epi8(tail128, head128));
                cmpim = _mm_or_si128(cmpim, _mm_cmpgt_epi8(head128, tail128));
                exit0 = blendv_fullmask8(exit0, zero128, cmpim);
            }

            j128 = _mm_set1_epi8(beg - i);
            int fast_lo = beg, fast_hi = beg;
            if (all_lanes_set8(exit0)) {
                fast_lo = i + hmax_epi8(head128);
                fast_hi = i + hmin_epi8(tail128);
                if (fast_lo < beg) fast_lo = beg;
                if (fast_lo > end) fast_lo = end;
                if (fast_hi < fast_lo) fast_hi = fast_lo;
                if (fast_hi > end) fast_hi = end;
            }
            /* the cell loops of smithWaterman128_8, through its own macros */
            if (!gen_mat) {
                EXT13_RUN_SPLIT8_128(EXT13_SBT8_XOR);
            } else if (fc.rank1) {
                __m128i rowfreed = _mm_cmpeq_epi8(s10, frref128);
                EXT13_RUN_SPLIT8_128(EXT13_SBT8_RANK1);
            } else {
                EXT13_RUN_SPLIT8_128(EXT13_SBT8_AMAT);
            }
            {
                __m128i cmp1 = _mm_or_si128(_mm_cmpgt_epi8(head128, j128), _mm_cmpgt_epi8(j128, tail128));
                h10 = blendv_fullmask8(h10, zero128, cmp1);
                _mm_store_si128((__m128i *)(H_h + j * SIMD_WIDTH8), h10);
                _mm_store_si128((__m128i *)(F + j * SIMD_WIDTH8), zero128);
            }

            /* BSW8_GSCORE_NARROW (smithWaterman128_8) folds into these, by name */
            uint8_t *gbest8 = V.gbest8, *gset8 = V.gset8;
            int16_t *ierow16 = V.ierow16;

            __m128i bmaxScore128 = maxScore128;
            __m128i tmp = _mm_cmpeq_epi8(maxRS1, zero128);
            if (vmaxvq_u8(vreinterpretq_u8_m128i(maxRS1)) == 0) {
                /* every lane's row max is zero: the plain kernel folds this row's query-end
                 * capture and leaves the group; every lane is finished. */
                BSW8_GSCORE_NARROW();
                _mm_store_si128((__m128i *) V.exit0, zero128);
                return 0;
            }
            exit0 = blendv_fullmask8(exit0, zero128, tmp);
            __m128i score128 = _mm_max_epu8(maxScore128, maxRS1);
            maxScore128 = blendv_fullmask8(maxScore128, score128, exit0);
            __m128i cmp = _mm_xor_si128(_mm_cmpeq_epi8(maxScore128, bmaxScore128), ff128);
            y128 = blendv_fullmask8(y128, y1_128, cmp);
            {
                __m128i t2 = _mm_abs_epi8(_mm_sub_epi8(y1_128, one128));
                __m128i bmax_off128 = max_off128;
                t2 = _mm_max_epu8(max_off128, t2);
                max_off128 = blendv_fullmask8(bmax_off128, t2, cmp);
            }
            {   /* lean epilogue: xrow16, gscore, gated z-drop */
                const __m128i vip1_16 = _mm_set1_epi16((int16_t) (i + 1));
                __m128i xl = _mm_load_si128((const __m128i *) V.xrow16);
                __m128i xh = _mm_load_si128((const __m128i *) (V.xrow16 + 8));
                xl = blendv_fullmask8(xl, vip1_16, _mm_unpacklo_epi8(cmp, cmp));
                xh = blendv_fullmask8(xh, vip1_16, _mm_unpackhi_epi8(cmp, cmp));
                _mm_store_si128((__m128i *) V.xrow16, xl);
                _mm_store_si128((__m128i *) (V.xrow16 + 8), xh);
                BSW8_GSCORE_NARROW();
                if (zdrop > 0) {
                    const __m128i need_z = _mm_subs_epu8(_mm_subs_epu8(maxScore128, maxRS1), zgate128);
                    if (!zgate_exact || any_lane_set8(need_z)) {
                        int8_t  y1_a[SIMD_WIDTH8]   __attribute((aligned(16)));
                        int8_t  y_a[SIMD_WIDTH8]    __attribute((aligned(16)));
                        int8_t  ms_a[SIMD_WIDTH8]   __attribute((aligned(16)));
                        int8_t  rs_a[SIMD_WIDTH8]   __attribute((aligned(16)));
                        int8_t  exit_a[SIMD_WIDTH8] __attribute((aligned(16)));
                        _mm_store_si128((__m128i *) y1_a, y1_128);
                        _mm_store_si128((__m128i *) y_a, y128);
                        _mm_store_si128((__m128i *) ms_a, maxScore128);
                        _mm_store_si128((__m128i *) rs_a, maxRS1);
                        _mm_store_si128((__m128i *) exit_a, exit0);
                        const __m128i vi   = _mm_set1_epi32(i);
                        const __m128i vip1 = _mm_set1_epi32(i + 1);
                        const __m128i vone = _mm_set1_epi32(1);
                        const __m128i vzd  = _mm_set1_epi32(zdrop);
                        const __m128i vedel = _mm_set1_epi32(this->e_del);
                        const __m128i veins = _mm_set1_epi32(this->e_ins);
                        __m128i die_g[SIMD_WIDTH8 / 4];
                        for (int g = 0; g < SIMD_WIDTH8 / 4; g++) {
                            const int base = g * 4;
                            __m128i msg  = _mm_cvtepu8_epi32(_mm_cvtsi32_si128(*(const int32_t *)(ms_a   + base)));
                            __m128i rsg  = _mm_cvtepu8_epi32(_mm_cvtsi32_si128(*(const int32_t *)(rs_a   + base)));
                            __m128i exitg= _mm_cvtepi8_epi32(_mm_cvtsi32_si128(*(const int32_t *)(exit_a + base)));
                            __m128i y1g  = _mm_cvtepi8_epi32(_mm_cvtsi32_si128(*(const int32_t *)(y1_a   + base)));
                            __m128i yg   = _mm_cvtepi8_epi32(_mm_cvtsi32_si128(*(const int32_t *)(y_a    + base)));
                            __m128i xrg  = _mm_cvtepi16_epi32(_mm_loadl_epi64((const __m128i *)(V.xrow16 + base)));
                            __m128i y1c  = _mm_add_epi32(y1g, vi);
                            __m128i yc   = _mm_add_epi32(yg, _mm_sub_epi32(xrg, vone));
                            yc = _mm_andnot_si128(_mm_cmpeq_epi32(xrg, _mm_setzero_si128()), yc);
                            __m128i tmpi = _mm_sub_epi32(vip1, xrg);
                            __m128i tmpj = _mm_sub_epi32(y1c, yc);
                            __m128i zdelta = _mm_sub_epi32(tmpi, tmpj);
                            __m128i zesel  = blendv_fullmask8(veins, vedel,
                                                 _mm_cmpgt_epi32(zdelta, _mm_setzero_si128()));
                            __m128i dif  = _mm_mullo_epi32(_mm_abs_epi32(zdelta), zesel);
                            __m128i drop = _mm_sub_epi32(msg, rsg);
                            __m128i die  = _mm_cmpgt_epi32(_mm_sub_epi32(drop, dif), vzd);
                            die_g[g] = _mm_and_si128(die, exitg);
                        }
                        __m128i die01 = _mm_packs_epi32(die_g[0], die_g[1]);
                        __m128i die23 = _mm_packs_epi32(die_g[2], die_g[3]);
                        exit0 = _mm_andnot_si128(_mm_packs_epi16(die01, die23), exit0);
                    }
                }
            }
            /* as in smithWaterman128_8: with no band columns maxRS1 stays 0 and the row returned above */
            xassert(beg < end, "lane compaction: empty band reached the narrowing pass");
            {   /* narrowing, one counted pass (BSW8_ROW_LEAN) */
                const __m128i big128 = _mm_set1_epi8(127);
                __m128i off128   = _mm_set1_epi8(beg - i);
                __m128i first128 = _mm_set1_epi8(end - i);
                __m128i last128  = _mm_set1_epi8(beg - 1 - i);
                for (int l = beg; l < end; l++) {
                    __m128i vv = _mm_or_si128(_mm_load_si128((__m128i *)(F + l * SIMD_WIDTH8)),
                                              _mm_load_si128((__m128i *)(H_h + l * SIMD_WIDTH8)));
                    __m128i nz = vreinterpretq_m128i_u8(vtstq_u8(vreinterpretq_u8_m128i(vv), vreinterpretq_u8_m128i(vv)));
                    first128 = _mm_min_epi8(first128, blendv_fullmask8(big128, off128, nz));
                    last128  = blendv_fullmask8(last128, off128, nz);
                    off128   = _mm_add_epi8(off128, one128);
                }
                {
                    __m128i vv = _mm_or_si128(_mm_load_si128((__m128i *)(F + end * SIMD_WIDTH8)),
                                              _mm_load_si128((__m128i *)(H_h + end * SIMD_WIDTH8)));
                    __m128i nz = vreinterpretq_m128i_u8(vtstq_u8(vreinterpretq_u8_m128i(vv), vreinterpretq_u8_m128i(vv)));
                    last128 = blendv_fullmask8(last128, off128, nz);
                }
                const int minfirst = hmin_epi8(first128);
                if (minfirst > beg - i) V.nbeg = i + minfirst - 1;
                const int l = i + hmax_epi8(last128);
                V.nend = l + 2 < vncol ? l + 2 : vncol;
                __m128i upd = _mm_andnot_si128(_mm_cmpeq_epi8(first128, _mm_set1_epi8(beg - i)), exit0);
                head128 = blendv_fullmask8(head128, first128, upd);
                upd = _mm_andnot_si128(_mm_cmpeq_epi8(last128, _mm_set1_epi8(end - i)), exit0);
                __m128i index128 = blendv_fullmask8(tail128, last128, upd);
                index128 = _mm_add_epi8(index128, two128);
                tail128 = _mm_min_epi8(index128, qlen_off128);
            }
            maxend128 = _mm_subs_epi8(maxend128, one128);   /* frame shift, saturating */
            head128 = _mm_sub_epi8(head128, one128);
#ifdef BSW8_ASSERT_ENVELOPE   /* debug only: see the byte-ceiling check in smithWaterman128_8 */
            assert(!_mm_movemask_epi8(_mm_cmpeq_epi8(
                       _mm_subs_epu8(_mm_set1_epi8((int8_t) BYTE_CEIL), maxRS1), zero128)) &&
                   "8-bit banded SW: row max hit the byte ceiling (lane compaction)");
#endif
            tail128 = _mm_sub_epi8(tail128, one128);
            _mm_store_si128((__m128i *) V.head, head128);
            _mm_store_si128((__m128i *) V.tail, tail128);
            _mm_store_si128((__m128i *) V.maxend, maxend128);
            _mm_store_si128((__m128i *) V.exit0, exit0);
            _mm_store_si128((__m128i *) V.maxScore, maxScore128);
            _mm_store_si128((__m128i *) V.y, y128);
            _mm_store_si128((__m128i *) V.max_off, max_off128);
            return bsw8c_lane_mask(exit0);
        };

        /* ---- move one live lane into a dead slot, between rows (next row r) ---- */
        auto move_lane = [&](Bsw8CVec &Sv, int sl, Bsw8CVec &D, int dl, int r) {
            writeback(D, dl);
            D.head[dl] = Sv.head[sl]; D.tail[dl] = Sv.tail[sl]; D.maxend[dl] = Sv.maxend[sl];
            D.exit0[dl] = Sv.exit0[sl]; D.maxScore[dl] = Sv.maxScore[sl]; D.y[dl] = Sv.y[sl];
            D.max_off[dl] = Sv.max_off[sl]; D.h0[dl] = Sv.h0[sl]; D.myband[dl] = Sv.myband[sl];
            D.gbest8[dl] = Sv.gbest8[sl]; D.gset8[dl] = Sv.gset8[sl];
            D.xrow16[dl] = Sv.xrow16[sl]; D.ierow16[dl] = Sv.ierow16[sl];
            D.qlenw[dl] = Sv.qlenw[sl]; D.mlenw[dl] = Sv.mlenw[sl]; D.tlenw[dl] = Sv.tlenw[sl];
            D.best_abs[dl] = Sv.best_abs[sl]; D.pair[dl] = Sv.pair[sl];
            int ha = r + (int8_t) Sv.head[sl];
            if (ha < 0) ha = 0;
            /* The lane's slot from its head to the superblock's last column, and its target
             * rows r to tlen - 1: everything the lane still reads. From row tlen on it is
             * finished (the mlen exit fires before that row's cells), so the target bytes
             * the slot keeps there from its previous lane are never used for it. */
            for (int c = ha; c <= ncol; c++) {
                D.H[c * SIMD_WIDTH8 + dl]    = Sv.H[c * SIMD_WIDTH8 + sl];
                D.F[c * SIMD_WIDTH8 + dl]    = Sv.F[c * SIMD_WIDTH8 + sl];
                D.seq2[c * SIMD_WIDTH8 + dl] = Sv.seq2[c * SIMD_WIDTH8 + sl];
            }
            for (int rr = r; rr < Sv.tlenw[sl]; rr++)
                D.seq1[rr * SIMD_WIDTH8 + dl] = Sv.seq1[rr * SIMD_WIDTH8 + sl];
            if (Sv.nbeg < D.nbeg) D.nbeg = Sv.nbeg;
            if (Sv.nend > D.nend) D.nend = Sv.nend;
            /* Left of a lane's head the kernel needs zero H/F wherever the group range reaches
             * (the trim/store invariant): those cells feed the horizontal carry (h11/e11) into
             * the lane's first in-band column. The slot still holds the previous lane's cells
             * there; clear them from the group's lowest future column (nbeg never decreases
             * except here) up to the head. */
            for (int c = D.nbeg < 0 ? 0 : D.nbeg; c < ha; c++) {
                D.H[c * SIMD_WIDTH8 + dl] = 0;
                D.F[c * SIMD_WIDTH8 + dl] = 0;
            }
            if (Sv.qlenw[sl] > D.ncol) D.ncol = Sv.qlenw[sl];
            if (Sv.qlenw[sl] - 1 < D.minq) D.minq = Sv.qlenw[sl] - 1;
            Sv.pair[sl] = nullptr; Sv.exit0[sl] = 0;
        };

        /* ---- row lockstep over the superblock ---- */
        int prev_live = -1;
        for (int i = 0; i < nrow; i++) {
            int total = 0, nact = 0;
            for (int v = 0; v < nv; v++) {
                if (!active[v]) { live[v] = 0; continue; }
                live[v] = row(Vs[v], i);
                total += __builtin_popcountll(live[v]);
                if (live[v]) nact++; else active[v] = 0;
            }
            if (nact == 0) break;
            if (nv > 1 && total != prev_live && nact > (total + SIMD_WIDTH8 - 1) / SIMD_WIDTH8) {
                const int nm = bsw_compact_plan(live, active, nv, SIMD_WIDTH8, moves);
                for (int k = 0; k < nm; k++)
                    move_lane(Vs[moves[k].sv], moves[k].sl, Vs[moves[k].dv], moves[k].dl, i + 1);
                int nleft = 0;
                for (int v = 0; v < nv; v++) nleft += active[v];
                compact_counts_.moves += (uint64_t) nm;
                compact_counts_.retired += (uint64_t) (nact - nleft);
            }
            prev_live = total;
        }
        for (int v = 0; v < nv; v++)
            for (int l = 0; l < SIMD_WIDTH8; l++) writeback(Vs[v], l);
    }
}

#endif  /* __aarch64__ */
