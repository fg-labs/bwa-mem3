/* bandedSWA_compact256.inc -- same-row lane compaction for the AVX2 8-bit extension kernel.
 * Included from bandedSWA.cpp right after smithWaterman256_8 (avx2 tier only). The driver is
 * bandedSWA_compact512.inc's, on 32-lane groups: superblocks of K groups in row lockstep,
 * bsw_compact_plan() after a row in which lanes finished, moved lanes copied slot to slot in the
 * absolute column frame, the receiving slot cleared left of the lane's head over the group range.
 * The row body is smithWaterman256_8's (its cell loops through the same EXT13_CELL8_256_*
 * macros) with its persistent vectors kept in a per-group struct.
 *
 * Any change to the row body here must be mirrored in smithWaterman256_8 (and vice versa);
 * test/unit/test_bandedswa_compact.cpp compares the two field by field. */
#if ((!__AVX512BW__) && (__AVX2__))

#include "bsw_compact.h"

static_assert(SIMD_WIDTH8 == 32, "bandedSWA_compact256.inc: Bsw8CVec256 holds 32 lanes");

namespace {
struct alignas(32) Bsw8CVec256 {
    uint8_t head[SIMD_WIDTH8], tail[SIMD_WIDTH8], maxend[SIMD_WIDTH8], exit0[SIMD_WIDTH8],
            maxScore[SIMD_WIDTH8], y[SIMD_WIDTH8], max_off[SIMD_WIDTH8], h0[SIMD_WIDTH8], myband[SIMD_WIDTH8];
    int32_t xrow[SIMD_WIDTH8], ierow[SIMD_WIDTH8], best_abs[SIMD_WIDTH8], gbest_abs[SIMD_WIDTH8],
            qlenw[SIMD_WIDTH8], mlenw[SIMD_WIDTH8], tlenw[SIMD_WIDTH8];
    SeqPair *pair[SIMD_WIDTH8];
    int nbeg, nend, ncol, minq;
    uint8_t *seq1, *seq2;
    int8_t *H, *F;
};
struct Bsw8CBuf256 {
    uint8_t *mem = nullptr; size_t cap = 0;
    ~Bsw8CBuf256() { if (mem) _mm_free(mem); }
};
}  // namespace
/* Per-thread scratch, grown on demand and kept until thread exit (see bandedSWA_compact512.inc). */
static thread_local Bsw8CBuf256 g_bsw8c_buf256;

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
    const size_t need = (size_t) K * (sizeof(Bsw8CVec256) + 4 * S) + S + 128;   // + one shared sbt_buf
    if (g_bsw8c_buf256.cap < need) {
        if (g_bsw8c_buf256.mem) _mm_free(g_bsw8c_buf256.mem);
        g_bsw8c_buf256.mem = (uint8_t *) _mm_malloc(need, 64);
        if (g_bsw8c_buf256.mem == NULL) { fprintf(stderr, "Error! Mem not allocated!!!\n"); exit(EXIT_FAILURE); }
        g_bsw8c_buf256.cap = need;
    }
    Bsw8CVec256 *Vs = (Bsw8CVec256 *) g_bsw8c_buf256.mem;
    uint8_t *slab = g_bsw8c_buf256.mem + (size_t) K * sizeof(Bsw8CVec256);
    slab = (uint8_t *) (((uintptr_t) slab + 63) & ~(uintptr_t) 63);
    for (int v = 0; v < K; v++) {
        Vs[v].seq1 = slab + (size_t) (4 * v + 0) * S;
        Vs[v].seq2 = slab + (size_t) (4 * v + 1) * S;
        Vs[v].H    = (int8_t *) (slab + (size_t) (4 * v + 2) * S);
        Vs[v].F    = (int8_t *) (slab + (size_t) (4 * v + 3) * S);
    }
    int8_t *sbt_buf = (int8_t *) (slab + (size_t) (4 * K) * S);   // per-row score pre-pass scratch

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

    /* the constants of smithWaterman256_8, under the names its cell macros use */
    __m256i match256     = _mm256_set1_epi8(this->w_match);
    __m256i w_ambig_256  = _mm256_set1_epi8(this->w_ambig);
    __m128i pmat128 = _mm_load_si128((__m128i *)bsw_pmat_bytes_);
    __m256i pmat256 = _mm256_broadcastsi128_si256(pmat128);
    const bool gen_mat = bsw_gen_mat_;
    const BswFreedCell fc = bsw_fc_;
    __m256i frref256  = _mm256_set1_epi8(fc.ref);
    __m256i frread256 = _mm256_set1_epi8(fc.read);
    __m256i amat256 = _mm256_broadcastsi128_si256(_mm_load_si128((__m128i *)bsw_amat_bytes_));
    __m256i three256_8 = _mm256_set1_epi8(3);
    __m256i e_del256    = _mm256_set1_epi8(this->e_del);
    __m256i oe_del256   = _mm256_set1_epi8(this->o_del + this->e_del);
    __m256i e_ins256    = _mm256_set1_epi8(this->e_ins);
    __m256i oe_ins256   = _mm256_set1_epi8(this->o_ins + this->e_ins);
    const __m256i zero256 = _mm256_setzero_si256();
    const __m256i one256  = _mm256_set1_epi8(1);
    const __m256i two256  = _mm256_set1_epi8(2);
    const __m256i ff256   = _mm256_set1_epi8(0xFF);
    const __m256i zdrop256 = _mm256_set1_epi8(zdrop);
    const uint32_t dmask = 0xFFFFFFFFu;
    const int o_del8 = (uint8_t) (int8_t) this->o_del, e_del8 = (uint8_t) (int8_t) this->e_del;
    /* band-edge write: byte penalties, as the seed row's set1_epi8 / subs_epu8 uses them */
    const int oe_ins8 = (uint8_t) (this->o_ins + this->e_ins), e_ins8 = (uint8_t) this->e_ins;

    uint64_t live[BSW_COMPACT_GROUPS_MAX];
    uint8_t active[BSW_COMPACT_GROUPS_MAX];
    BswMove moves[BSW_COMPACT_GROUPS_MAX * SIMD_WIDTH8];

    auto writeback = [&](Bsw8CVec256 &V, int l) {
        SeqPair *p = V.pair[l];
        if (!p) return;
        p->score   = V.best_abs[l];
        p->tle     = V.xrow[l];
        p->qle     = (V.xrow[l] == 0) ? 0 : ((int) (int8_t) V.y[l] + V.xrow[l] - 1);
        p->max_off = (uint8_t) V.max_off[l];
        p->gscore  = V.gbest_abs[l];
        p->gtle    = V.ierow[l];
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
        for (int v = 0; v < nv; v++) {
            Bsw8CVec256 &V = Vs[v];
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
                /* the plain wrapper's per-lane band clamp, on len2 * max_sc */
                const int reach = p[l].len2 * maxsc + eb;
                int max_ins = (int) ((double) (reach - o_ins) / e_ins + 1.0); if (max_ins < 1) max_ins = 1;
                int max_del = (int) ((double) (reach - o_del) / e_del + 1.0); if (max_del < 1) max_del = 1;
                int band = w; if (max_ins < band) band = max_ins; if (max_del < band) band = max_del;
                V.myband[l] = (uint8_t) band;
                V.tlenw[l] = p[l].len1; V.qlenw[l] = p[l].len2;
                int ml = p[l].len2 + V.myband[l]; if (ml > p[l].len1) ml = p[l].len1;
                V.mlenw[l] = ml;
                V.xrow[l] = 0; V.ierow[l] = 0; V.best_abs[l] = p[l].h0; V.gbest_abs[l] = -1;
                V.pair[l] = &p[l];
                V.head[l] = 0; V.tail[l] = 127; V.maxend[l] = (uint8_t) -128; V.exit0[l] = 0xFF;
                V.maxScore[l] = V.h0[l]; V.y[l] = 0; V.max_off[l] = 0;
            }
            V.nbeg = 0; V.nend = vcol; V.ncol = vcol; V.minq = minq - 1;
            x86_soa_pack<SIMD_WIDTH8>(V.seq1, seq1p, len1a, len1a, nrow + 1, DUMMY1, DUMMY1, false, AMBIG, 8);
            x86_soa_pack<SIMD_WIDTH8>(V.seq2, seq2p, len2a, len2a, ncol + 1, DUMMY2, DUMMY2, true, AMBIG, 8);
            /* the h0-prefix insertion seed row, as the plain wrapper's H1 */
            const __m256i h0_256 = _mm256_load_si256((const __m256i *) V.h0);
            _mm256_store_si256((__m256i *) (V.H + ncol * SIMD_WIDTH8), zero256);
            _mm256_store_si256((__m256i *) V.H, h0_256);
            __m256i t = _mm256_subs_epu8(h0_256, oe_ins256);
            _mm256_store_si256((__m256i *) (V.H + SIMD_WIDTH8), t);
            for (int k = 2; k < ncol; k++) {
                t = _mm256_subs_epu8(t, e_ins256);
                _mm256_store_si256((__m256i *) (V.H + k * SIMD_WIDTH8), t);
            }
            for (int j = 0; j <= ncol; j++) _mm256_store_si256((__m256i *) (V.F + j * SIMD_WIDTH8), zero256);
            active[v] = 1;
        }

        /* ---- one row of one group: the row of smithWaterman256_8 ---- */
        auto row = [&](Bsw8CVec256 &V, const int i) -> uint64_t {
            int8_t *H_h = V.H, *F = V.F;
            const uint8_t *seq2SoA = V.seq2;
            const int minq = V.minq, vncol = V.ncol;
            int32_t *xrow = V.xrow, *ierow = V.ierow, *best_abs = V.best_abs, *gbest_abs = V.gbest_abs;
            __m256i head256     = _mm256_load_si256((const __m256i *) V.head);
            __m256i tail256     = _mm256_load_si256((const __m256i *) V.tail);
            __m256i maxend256   = _mm256_load_si256((const __m256i *) V.maxend);
            __m256i exit0       = _mm256_load_si256((const __m256i *) V.exit0);
            __m256i maxScore256 = _mm256_load_si256((const __m256i *) V.maxScore);
            __m256i y256        = _mm256_load_si256((const __m256i *) V.y);
            __m256i max_off256  = _mm256_load_si256((const __m256i *) V.max_off);
            const __m256i myband256 = _mm256_load_si256((const __m256i *) V.myband);
            const __m256i h0_256    = _mm256_load_si256((const __m256i *) V.h0);
            int j;

            __m256i e11 = zero256;
            __m256i h00, h11, h10;
            __m256i s10 = _mm256_load_si256((const __m256i *) (V.seq1 + i * SIMD_WIDTH8));
            int beg = V.nbeg, end = V.nend;
            if (beg < i - w) beg = i - w;
            if (end > i + w + 1) end = i + w + 1;
            if (end > vncol) end = vncol;
            /* column -1: the h0-prefix deletion seed the plain wrapper keeps in H2 */
            h10 = zero256;
            if (beg == 0) {
                long sub = (long) o_del8 + (long) (i + 1) * e_del8;
                if (sub > 255) sub = 255;
                h10 = _mm256_subs_epu8(h0_256, _mm256_set1_epi8((char) (uint8_t) sub));
            }
            __m256i j256 = zero256;
            __m256i maxRS1 = zero256;
            __m256i y1_256 = zero256;
            __m256i hqe256 = zero256, qfire256 = zero256;

            int8_t qlen_off_a[SIMD_WIDTH8]   __attribute((aligned(32)));
            int8_t qlen_valid_a[SIMD_WIDTH8] __attribute((aligned(32)));
            int8_t cmpim_a[SIMD_WIDTH8]      __attribute((aligned(32)));
            for (int l = 0; l < SIMD_WIDTH8; l++) {
                int qoff = V.qlenw[l] - i;
                int qoff_sat = qoff;
                if (qoff_sat >  127) qoff_sat =  127;
                if (qoff_sat < -128) qoff_sat = -128;
                qlen_off_a[l]   = (int8_t) qoff_sat;
                qlen_valid_a[l] = (qoff >= -w && qoff <= w + 1) ? (int8_t) 0xFF : 0;
                cmpim_a[l]      = ((i + 1) > V.mlenw[l]) ? (int8_t) 0xFF : 0;
            }
            __m256i qlen_off256   = _mm256_load_si256((__m256i *) qlen_off_a);
            __m256i qlen_valid256 = _mm256_load_si256((__m256i *) qlen_valid_a);

            __m256i phead256 = head256, ptail256 = tail256;
            head256 = _mm256_max_epi8(head256, _mm256_sub_epi8(zero256, myband256));
            tail256 = _mm256_min_epi8(tail256, _mm256_add_epi8(myband256, one256));
            tail256 = _mm256_min_epi8(tail256, qlen_off256);
            {   /* band-edge write, as in smithWaterman256_8 */
                __m256i need_ = _mm256_and_si256(_mm256_cmpeq_epi8(_mm256_sub_epi8(tail256, maxend256), two256), exit0);
                uint32_t nm_ = (uint32_t) _mm256_movemask_epi8(need_);
                if (nm_) {
                    int8_t me_a[SIMD_WIDTH8] __attribute((aligned(32)));
                    _mm256_store_si256((__m256i *) me_a, maxend256);
                    while (nm_) { const int l = __builtin_ctz(nm_); nm_ &= nm_ - 1;
                        const int c = i + me_a[l] + 1;
                        H_h[c * SIMD_WIDTH8 + l] = (int8_t) (uint8_t) bsw_seed_row(V.h0[l], oe_ins8, e_ins8, c);
                        F[c * SIMD_WIDTH8 + l] = 0; }
                }
                maxend256 = _mm256_max_epi8(maxend256, tail256);
            }
            {   /* trimming */
                __m256i cmph = _mm256_and_si256(_mm256_cmpeq_epi8(head256, phead256), _mm256_cmpeq_epi8(tail256, ptail256));
                const uint32_t cmp_ht = (uint32_t) _mm256_movemask_epi8(cmph);
                for (int l = beg; l < end && cmp_ht != dmask; l++) {
                    __m256i h256 = _mm256_load_si256((__m256i *)(H_h + l * SIMD_WIDTH8));
                    __m256i f256 = _mm256_load_si256((__m256i *)(F + l * SIMD_WIDTH8));
                    __m256i pj256 = _mm256_set1_epi8(l - i);
                    __m256i cmp1 = _mm256_cmpgt_epi8(head256, pj256);
                    if (_mm256_movemask_epi8(cmp1) == 0) break;
                    cmp1 = _mm256_or_si256(cmp1, _mm256_cmpgt_epi8(pj256, tail256));
                    h256 = _mm256_andnot_si256(cmp1, h256);
                    f256 = _mm256_andnot_si256(cmp1, f256);
                    _mm256_store_si256((__m256i *)(F + l * SIMD_WIDTH8), f256);
                    _mm256_store_si256((__m256i *)(H_h + l * SIMD_WIDTH8), h256);
                }
            }
            {
                __m256i cmpim = _mm256_load_si256((__m256i *) cmpim_a);
                cmpim = _mm256_or_si256(cmpim, _mm256_cmpeq_epi8(tail256, head256));
                cmpim = _mm256_or_si256(cmpim, _mm256_cmpgt_epi8(head256, tail256));
                exit0 = _mm256_andnot_si256(cmpim, exit0);
            }
            if (!gen_mat) {
                for (int jp = beg; jp < end; jp++) {
                    __m256i s2 = _mm256_load_si256((__m256i *)(seq2SoA + jp * SIMD_WIDTH8));
                    __m256i sbt11;
                    SBT_PREPASS8_XOR(s10, s2, sbt11, pmat256);
                    _mm256_store_si256((__m256i *)(sbt_buf + jp * SIMD_WIDTH8), sbt11);
                }
            } else if (fc.rank1) {
                __m256i rowfreed = _mm256_cmpeq_epi8(s10, frref256);
                for (int jp = beg; jp < end; jp++) {
                    __m256i s2 = _mm256_load_si256((__m256i *)(seq2SoA + jp * SIMD_WIDTH8));
                    __m256i sbt11;
                    SBT_PREPASS8_RANK1(s10, s2, rowfreed, sbt11, pmat256, match256, frread256);
                    _mm256_store_si256((__m256i *)(sbt_buf + jp * SIMD_WIDTH8), sbt11);
                }
            } else {
                for (int jp = beg; jp < end; jp++) {
                    __m256i s2 = _mm256_load_si256((__m256i *)(seq2SoA + jp * SIMD_WIDTH8));
                    __m256i sbt11;
                    SBT_PREPASS8_AMAT(s10, s2, sbt11, amat256, w_ambig_256, three256_8);
                    _mm256_store_si256((__m256i *)(sbt_buf + jp * SIMD_WIDTH8), sbt11);
                }
            }
            int fast_lo = beg, fast_hi = beg;
            if (_mm256_movemask_epi8(exit0) == -1) {
                int8_t hh_[SIMD_WIDTH8] __attribute((aligned(32)));
                int8_t tt_[SIMD_WIDTH8] __attribute((aligned(32)));
                _mm256_store_si256((__m256i *) hh_, head256);
                _mm256_store_si256((__m256i *) tt_, tail256);
                int maxhead = -128, mintail = 127;
                for (int l = 0; l < SIMD_WIDTH8; l++) {
                    if (hh_[l] > maxhead) maxhead = hh_[l];
                    if (tt_[l] < mintail) mintail = tt_[l];
                }
                fast_lo = i + maxhead; if (fast_lo < beg) fast_lo = beg; if (fast_lo > end) fast_lo = end;
                fast_hi = i + mintail; if (fast_hi < fast_lo) fast_hi = fast_lo; if (fast_hi > end) fast_hi = end;
            }
            j256 = _mm256_set1_epi8(beg - i);
            /* the cell loops of smithWaterman256_8, through its own macros */
#pragma unroll(4)
            for (j = beg; j < fast_lo; j++)   EXT13_CELL8_256_MASKED
#pragma unroll(4)
            for (j = fast_lo; j < fast_hi; j++) EXT13_CELL8_256_FAST
#pragma unroll(4)
            for (j = fast_hi; j < end; j++)   EXT13_CELL8_256_MASKED
            {
                __m256i cmp1 = _mm256_or_si256(_mm256_cmpgt_epi8(head256, j256), _mm256_cmpgt_epi8(j256, tail256));
                h10 = _mm256_andnot_si256(cmp1, h10);
                _mm256_store_si256((__m256i *)(H_h + j * SIMD_WIDTH8), h10);
                _mm256_store_si256((__m256i *)(F + j * SIMD_WIDTH8), zero256);
            }
            /* exit due to zero score by a row */
            __m256i bmaxScore256 = maxScore256;
            __m256i tmp = _mm256_cmpeq_epi8(maxRS1, zero256);
            if ((uint32_t) _mm256_movemask_epi8(tmp) == dmask) {
                /* finalize this row's query-end capture, as smithWaterman256_8's zero-row break */
                int8_t qf_a_[SIMD_WIDTH8] __attribute((aligned(32)));
                int8_t hq_a_[SIMD_WIDTH8] __attribute((aligned(32)));
                _mm256_store_si256((__m256i *) qf_a_, qfire256);
                _mm256_store_si256((__m256i *) hq_a_, hqe256);
                for (int g = 0; g < SIMD_WIDTH8 / 8; g++) {
                    const int base = g * 8;
                    __m256i hqeg = _mm256_cvtepu8_epi32(_mm_loadl_epi64((const __m128i *)(hq_a_ + base)));
                    __m256i qfg  = _mm256_cvtepi8_epi32(_mm_loadl_epi64((const __m128i *)(qf_a_ + base)));
                    __m256i gba = _mm256_loadu_si256((const __m256i *)(gbest_abs + base));
                    __m256i ge  = _mm256_xor_si256(_mm256_cmpgt_epi32(gba, hqeg), ff256);
                    __m256i gmask = _mm256_and_si256(qfg, ge);
                    gba = _mm256_blendv_epi8(gba, hqeg, gmask);
                    _mm256_storeu_si256((__m256i *)(gbest_abs + base), gba);
                    __m256i ierg = _mm256_loadu_si256((const __m256i *)(ierow + base));
                    ierg = _mm256_blendv_epi8(ierg, _mm256_set1_epi32(i + 1), gmask);
                    _mm256_storeu_si256((__m256i *)(ierow + base), ierg);
                }
                _mm256_store_si256((__m256i *) V.exit0, zero256);
                return 0;
            }
            exit0 = _mm256_andnot_si256(tmp, exit0);
            __m256i score256 = _mm256_max_epu8(maxScore256, maxRS1);
            maxScore256 = _mm256_blendv_epi8(maxScore256, score256, exit0);
            __m256i cmp = _mm256_xor_si256(_mm256_cmpeq_epi8(maxScore256, bmaxScore256), ff256);
            y256 = _mm256_blendv_epi8(y256, y1_256, cmp);
            {
                __m256i ind256 = _mm256_abs_epi8(_mm256_sub_epi8(y1_256, one256));
                __m256i bmax_off256 = max_off256;
                ind256 = _mm256_max_epu8(max_off256, ind256);
                max_off256 = _mm256_blendv_epi8(bmax_off256, ind256, cmp);
            }
            /* the per-lane wide updates of smithWaterman256_8, behind the same gate */
            const __m256i need_z = (zdrop > 0)
                ? _mm256_subs_epu8(_mm256_subs_epu8(maxScore256, maxRS1), zdrop256)
                : zero256;
            const __m256i need_any = _mm256_or_si256(_mm256_or_si256(cmp, qfire256), need_z);
            if (!_mm256_testz_si256(need_any, need_any)) {
                int8_t  cmp_a[SIMD_WIDTH8]      __attribute((aligned(32)));
                int8_t  y1_a[SIMD_WIDTH8]       __attribute((aligned(32)));
                int8_t  y_a[SIMD_WIDTH8]        __attribute((aligned(32)));
                int8_t  ms_a[SIMD_WIDTH8]       __attribute((aligned(32)));
                int8_t  rs_a[SIMD_WIDTH8]       __attribute((aligned(32)));
                int8_t  exit_a[SIMD_WIDTH8]     __attribute((aligned(32)));
                int8_t  qfire_a[SIMD_WIDTH8]    __attribute((aligned(32)));
                int8_t  hqe_a[SIMD_WIDTH8]      __attribute((aligned(32)));
                _mm256_store_si256((__m256i *) cmp_a, cmp);
                _mm256_store_si256((__m256i *) y1_a, y1_256);
                _mm256_store_si256((__m256i *) y_a, y256);
                _mm256_store_si256((__m256i *) ms_a, maxScore256);
                _mm256_store_si256((__m256i *) rs_a, maxRS1);
                _mm256_store_si256((__m256i *) exit_a, exit0);
                _mm256_store_si256((__m256i *) qfire_a, qfire256);
                _mm256_store_si256((__m256i *) hqe_a, hqe256);
                const __m256i vi   = _mm256_set1_epi32(i);
                const __m256i vip1 = _mm256_set1_epi32(i + 1);
                const __m256i vone = _mm256_set1_epi32(1);
                const __m256i vzd  = _mm256_set1_epi32(zdrop);
                const __m256i vedel = _mm256_set1_epi32(this->e_del);
                const __m256i veins = _mm256_set1_epi32(this->e_ins);
                for (int g = 0; g < SIMD_WIDTH8 / 8; g++) {
                    const int base = g * 8;
                    __m256i msg  = _mm256_cvtepu8_epi32(_mm_loadl_epi64((const __m128i *)(ms_a    + base)));
                    __m256i rsg  = _mm256_cvtepu8_epi32(_mm_loadl_epi64((const __m128i *)(rs_a    + base)));
                    __m256i hqeg = _mm256_cvtepu8_epi32(_mm_loadl_epi64((const __m128i *)(hqe_a   + base)));
                    __m256i cmpg = _mm256_cvtepi8_epi32(_mm_loadl_epi64((const __m128i *)(cmp_a   + base)));
                    __m256i exitg= _mm256_cvtepi8_epi32(_mm_loadl_epi64((const __m128i *)(exit_a  + base)));
                    __m256i qfg  = _mm256_cvtepi8_epi32(_mm_loadl_epi64((const __m128i *)(qfire_a + base)));
                    __m256i y1g  = _mm256_cvtepi8_epi32(_mm_loadl_epi64((const __m128i *)(y1_a    + base)));
                    __m256i yg   = _mm256_cvtepi8_epi32(_mm_loadl_epi64((const __m128i *)(y_a     + base)));
                    __m256i xrg = _mm256_loadu_si256((const __m256i *)(xrow + base));
                    xrg = _mm256_blendv_epi8(xrg, vip1, cmpg);
                    _mm256_storeu_si256((__m256i *)(xrow + base), xrg);
                    __m256i bag = _mm256_loadu_si256((const __m256i *)(best_abs + base));
                    bag = _mm256_blendv_epi8(bag, _mm256_max_epi32(bag, msg), cmpg);
                    _mm256_storeu_si256((__m256i *)(best_abs + base), bag);
                    __m256i gba = _mm256_loadu_si256((const __m256i *)(gbest_abs + base));
                    __m256i ge  = _mm256_xor_si256(_mm256_cmpgt_epi32(gba, hqeg), ff256);
                    __m256i gmask = _mm256_and_si256(qfg, ge);
                    gba = _mm256_blendv_epi8(gba, hqeg, gmask);
                    _mm256_storeu_si256((__m256i *)(gbest_abs + base), gba);
                    __m256i ierg = _mm256_loadu_si256((const __m256i *)(ierow + base));
                    ierg = _mm256_blendv_epi8(ierg, vip1, gmask);
                    _mm256_storeu_si256((__m256i *)(ierow + base), ierg);
                    __m256i y1c  = _mm256_add_epi32(y1g, vi);
                    __m256i yc   = _mm256_add_epi32(yg, _mm256_sub_epi32(xrg, vone));
                    yc = _mm256_andnot_si256(_mm256_cmpeq_epi32(xrg, _mm256_setzero_si256()), yc);
                    __m256i tmpi = _mm256_sub_epi32(vip1, xrg);
                    __m256i tmpj = _mm256_sub_epi32(y1c, yc);
                    __m256i zdelta = _mm256_sub_epi32(tmpi, tmpj);
                    __m256i zesel  = _mm256_blendv_epi8(veins, vedel,
                                         _mm256_cmpgt_epi32(zdelta, _mm256_setzero_si256()));
                    __m256i dif  = _mm256_mullo_epi32(_mm256_abs_epi32(zdelta), zesel);
                    __m256i drop = _mm256_sub_epi32(msg, rsg);
                    __m256i die  = _mm256_cmpgt_epi32(_mm256_sub_epi32(drop, dif), vzd);
                    die = _mm256_and_si256(die, exitg);
                    __m128i d16 = _mm_packs_epi32(_mm256_castsi256_si128(die),
                                                  _mm256_extracti128_si256(die, 1));
                    __m128i d8  = _mm_packs_epi16(d16, d16);
                    __m128i ex  = _mm_loadl_epi64((const __m128i *)(exit_a + base));
                    if (zdrop > 0) _mm_storel_epi64((__m128i *)(exit_a + base), _mm_andnot_si128(d8, ex));
                }
                exit0 = _mm256_load_si256((__m256i *) exit_a);
            }
            {   /* narrowing (the four scans) */
                int l;
                for (l = beg; l < end; l++) {
                    __m256i tmp_ = _mm256_or_si256(_mm256_load_si256((__m256i *)(F + l * SIMD_WIDTH8)),
                                                   _mm256_load_si256((__m256i *)(H_h + l * SIMD_WIDTH8)));
                    if ((uint32_t) _mm256_movemask_epi8(_mm256_cmpeq_epi8(tmp_, zero256)) == dmask) V.nbeg = l;
                    else break;
                }
                for (l = end; l >= beg; l--) {
                    __m256i tmp_ = _mm256_or_si256(_mm256_load_si256((__m256i *)(F + l * SIMD_WIDTH8)),
                                                   _mm256_load_si256((__m256i *)(H_h + l * SIMD_WIDTH8)));
                    if ((uint32_t) _mm256_movemask_epi8(_mm256_cmpeq_epi8(tmp_, zero256)) != dmask) break;
                }
                V.nend = l + 2 < vncol ? l + 2 : vncol;
                __m256i tmpb = ff256;
                __m256i exit1 = _mm256_xor_si256(exit0, ff256);
                __m256i l256 = _mm256_set1_epi8(beg - i);
                for (l = beg; l < end; l++) {
                    __m256i tmp_ = _mm256_or_si256(_mm256_load_si256((__m256i *)(F + l * SIMD_WIDTH8)),
                                                   _mm256_load_si256((__m256i *)(H_h + l * SIMD_WIDTH8)));
                    tmp_ = _mm256_or_si256(tmp_, exit1);
                    tmp_ = _mm256_cmpeq_epi8(tmp_, zero256);
                    if (_mm256_movemask_epi8(tmp_) == 0) break;
                    tmp_ = _mm256_and_si256(tmp_, tmpb);
                    l256 = _mm256_add_epi8(l256, one256);
                    head256 = _mm256_blendv_epi8(head256, l256, tmp_);
                    tmpb = tmp_;
                }
                __m256i index256 = tail256;
                tmpb = ff256;
                l256 = _mm256_set1_epi8(end - i);
                for (l = end; l >= beg; l--) {
                    __m256i tmp_ = _mm256_or_si256(_mm256_load_si256((__m256i *)(F + l * SIMD_WIDTH8)),
                                                   _mm256_load_si256((__m256i *)(H_h + l * SIMD_WIDTH8)));
                    tmp_ = _mm256_or_si256(tmp_, exit1);
                    tmp_ = _mm256_cmpeq_epi8(tmp_, zero256);
                    if (_mm256_movemask_epi8(tmp_) == 0) break;
                    tmp_ = _mm256_and_si256(tmp_, tmpb);
                    l256 = _mm256_sub_epi8(l256, one256);
                    index256 = _mm256_blendv_epi8(index256, l256, tmp_);
                    tmpb = tmp_;
                }
                index256 = _mm256_add_epi8(index256, two256);
                tail256 = _mm256_min_epi8(index256, qlen_off256);
            }
            maxend256 = _mm256_subs_epi8(maxend256, one256);   /* frame shift, saturating */
            head256 = _mm256_sub_epi8(head256, one256);
            tail256 = _mm256_sub_epi8(tail256, one256);
#ifdef BSW8_ASSERT_ENVELOPE   /* debug only: see the byte-ceiling check in smithWaterman256_8 */
            assert(!_mm256_movemask_epi8(_mm256_cmpeq_epi8(
                       _mm256_subs_epu8(_mm256_set1_epi8((int8_t) BYTE_CEIL), maxRS1), zero256)) &&
                   "8-bit banded SW: row max hit the byte ceiling (lane compaction)");
#endif
            _mm256_store_si256((__m256i *) V.head, head256);
            _mm256_store_si256((__m256i *) V.tail, tail256);
            _mm256_store_si256((__m256i *) V.maxend, maxend256);
            _mm256_store_si256((__m256i *) V.exit0, exit0);
            _mm256_store_si256((__m256i *) V.maxScore, maxScore256);
            _mm256_store_si256((__m256i *) V.y, y256);
            _mm256_store_si256((__m256i *) V.max_off, max_off256);
            return (uint64_t) (uint32_t) _mm256_movemask_epi8(exit0);
        };

        auto move_lane = [&](Bsw8CVec256 &Sv, int sl, Bsw8CVec256 &D, int dl, int r) {
            writeback(D, dl);
            D.head[dl] = Sv.head[sl]; D.tail[dl] = Sv.tail[sl]; D.maxend[dl] = Sv.maxend[sl];
            D.exit0[dl] = Sv.exit0[sl]; D.maxScore[dl] = Sv.maxScore[sl]; D.y[dl] = Sv.y[sl];
            D.max_off[dl] = Sv.max_off[sl]; D.h0[dl] = Sv.h0[sl]; D.myband[dl] = Sv.myband[sl];
            D.xrow[dl] = Sv.xrow[sl]; D.ierow[dl] = Sv.ierow[sl]; D.best_abs[dl] = Sv.best_abs[sl];
            D.gbest_abs[dl] = Sv.gbest_abs[sl];
            D.qlenw[dl] = Sv.qlenw[sl]; D.mlenw[dl] = Sv.mlenw[sl]; D.tlenw[dl] = Sv.tlenw[sl];
            D.pair[dl] = Sv.pair[sl];
            int ha = r + (int8_t) Sv.head[sl];
            if (ha < 0) ha = 0;
            /* As in bandedSWA_compact.inc: the slot from the head on, target rows r to tlen - 1. */
            for (int c = ha; c <= ncol; c++) {
                D.H[c * SIMD_WIDTH8 + dl]    = Sv.H[c * SIMD_WIDTH8 + sl];
                D.F[c * SIMD_WIDTH8 + dl]    = Sv.F[c * SIMD_WIDTH8 + sl];
                D.seq2[c * SIMD_WIDTH8 + dl] = Sv.seq2[c * SIMD_WIDTH8 + sl];
            }
            for (int rr = r; rr < Sv.tlenw[sl]; rr++)
                D.seq1[rr * SIMD_WIDTH8 + dl] = Sv.seq1[rr * SIMD_WIDTH8 + sl];
            if (Sv.nbeg < D.nbeg) D.nbeg = Sv.nbeg;
            if (Sv.nend > D.nend) D.nend = Sv.nend;
            for (int c = D.nbeg < 0 ? 0 : D.nbeg; c < ha; c++) {   /* see move_lane in bandedSWA_compact.inc */
                D.H[c * SIMD_WIDTH8 + dl] = 0;
                D.F[c * SIMD_WIDTH8 + dl] = 0;
            }
            if (Sv.qlenw[sl] > D.ncol) D.ncol = Sv.qlenw[sl];
            if (Sv.qlenw[sl] - 1 < D.minq) D.minq = Sv.qlenw[sl] - 1;
            Sv.pair[sl] = nullptr; Sv.exit0[sl] = 0;
        };

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

#endif  /* AVX2 (not AVX-512BW) */
