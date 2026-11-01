/* bandedSWA_compact512.inc -- same-row lane compaction for the AVX-512BW 8-bit extension kernel.
 * Included from bandedSWA.cpp right after smithWaterman512_8Impl (avx512bw tier only). The driver is
 * bandedSWA_compact.inc's (NEON), on 64-lane groups: superblocks of K groups in row lockstep,
 * bsw_compact_plan() after a row in which lanes finished, moved lanes copied slot to slot in the
 * absolute column frame, the receiving slot cleared left of the lane's head over the group range.
 * The row body is smithWaterman512_8Impl's lean row (its cell loops through the same EXT13_CELL8_512_*
 * macros, which record the narrowing offsets; the counted trim; the narrow epilogue with the gated
 * z-drop, BSW8_GSCORE_NARROW512 and BSW8_ZDROP_WIDE512) with its persistent vectors kept in a
 * per-group struct.
 *
 * Any change to the row body here must be mirrored in smithWaterman512_8Impl (and vice versa);
 * test/unit/test_bandedswa_compact.cpp compares the two field by field. */
#if __AVX512BW__

#include <algorithm>
#include <cstring>

#include "bsw_compact.h"

static_assert(SIMD_WIDTH8 == 64, "bandedSWA_compact512.inc: Bsw8CVec512 holds 64 lanes");

namespace {
struct alignas(64) Bsw8CVec512 {
    uint8_t head[64], tail[64], maxend[64], exit0[64], maxScore[64], y[64], max_off[64], h0[64], myband[64];
    uint8_t gbest8[64], gset8[64];      /* gscore byte, and 0xFF where it is set (-1 otherwise) */
    int16_t xrow16[64], ierow16[64];    /* best row and best-gscore row (rows < MAX_SEQ_LEN8) */
    int32_t h0_raw[64], qlenw[64], mlenw[64], tlenw[64];   /* h0_raw: the raw seed score (see smithWaterman512_8Impl) */
    SeqPair *pair[64];
    int nbeg, nend, ncol, minq;
    uint8_t *seq1, *seq2;
    int8_t *H, *F;
};
struct Bsw8CBuf512 {
    uint8_t *mem = nullptr; size_t cap = 0;
    ~Bsw8CBuf512() { if (mem) _mm_free(mem); }
};
}  // namespace
/* Per-thread scratch, grown on demand and kept until thread exit: the pipeline builds a new
 * kernel object for every batch (mem_chain2aln_across_reads_V2), so a per-object buffer would
 * be allocated and faulted in again each time. */
static thread_local Bsw8CBuf512 g_bsw8c_buf512;

template <bool SIGNED8, bool CPG>
void BandedPairWiseSW::smithWatermanBatchWrapper8CompactImpl(SeqPair *pairArray,
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
    const size_t need = (size_t) K * (sizeof(Bsw8CVec512) + 4 * S) + S + 128;   // + one shared sbt_buf
    if (g_bsw8c_buf512.cap < need) {
        if (g_bsw8c_buf512.mem) _mm_free(g_bsw8c_buf512.mem);
        g_bsw8c_buf512.mem = (uint8_t *) _mm_malloc(need, 64);
        if (g_bsw8c_buf512.mem == NULL) { fprintf(stderr, "Error! Mem not allocated!!!\n"); exit(EXIT_FAILURE); }
        g_bsw8c_buf512.cap = need;
    }
    Bsw8CVec512 *Vs = (Bsw8CVec512 *) g_bsw8c_buf512.mem;
    uint8_t *slab = g_bsw8c_buf512.mem + (size_t) K * sizeof(Bsw8CVec512);
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

    __m512i match512     = _mm512_set1_epi8(this->w_match);
    __m512i w_ambig_512  = _mm512_set1_epi8(this->w_ambig);
    __m128i pmat128 = _mm_load_si128((__m128i *)bsw_pmat_bytes_);
    __m512i pmat512 = _mm512_broadcast_i32x4(pmat128);
    const bool gen_mat = bsw_gen_mat_;
    const BswFreedCell fc = bsw_fc_;
    __m512i frref512  = _mm512_set1_epi8(fc.ref);
    __m512i frread512 = _mm512_set1_epi8(fc.read);
    __m512i amat512 = _mm512_broadcast_i32x4(_mm_load_si128((__m128i *)bsw_amat_bytes_));
    __m512i three512_8 = _mm512_set1_epi8(3);
    // CpG arm (CPG only; see "CpG arm" in bandedSWA.cpp), as smithWaterman512_8Impl
    __m512i cpgref512  = _mm512_set1_epi8(this->cpg_ref_code_);
    __m512i cpgread512 = _mm512_set1_epi8(this->cpg_read_);
    __m512i cpgval512  = _mm512_set1_epi8(this->cpg_val_);
    __m512i e_del512    = _mm512_set1_epi8(this->e_del);
    __m512i oe_del512   = _mm512_set1_epi8(this->o_del + this->e_del);
    __m512i e_ins512    = _mm512_set1_epi8(this->e_ins);
    __m512i oe_ins512   = _mm512_set1_epi8(this->o_ins + this->e_ins);
    const __m512i zero512 = _mm512_setzero_si512();
    /* score domain (see smithWaterman512_8Impl): the byte of score 0 and the bias of a score byte.
     * The empty asm keeps 0x80 in a register: left as a known constant, clang folds it into the
     * per-cell restart compare of MAIN_CODE8_CORE as a rip-relative memory operand, a fourth load
     * per cell. */
    __m512i x80_512 = _mm512_set1_epi8((char) 0x80);
    __asm__("" : "+v"(x80_512));
    const __m512i hzero512 = SIGNED8 ? x80_512 : zero512;   /* also the XOR bias of a score byte:
                                                              * the seed stores below xor with it */
    const uint8_t hbias8   = SIGNED8 ? 0x80 : 0;
    const __m512i one512  = _mm512_set1_epi8(1);
    const __m512i two512  = _mm512_set1_epi8(2);
    const __m512i ff512   = _mm512_set1_epi8(0xFF);
    /* the epilogue skip below is exact only for gap extends inside this bound */
    const bool zgate_exact = BSW8_ZDROP_GATE_EXACT(this->e_del, this->e_ins);
    /* lean z-drop gate threshold byte: drop <= 255, so for zdrop >= 255 no lane can die */
    const __m512i zgate512 = _mm512_set1_epi8((int8_t) (zdrop < 255 ? zdrop : 255));
    const __mmask64 dmask = 0xFFFFFFFFFFFFFFFFULL;
    const int o_del8 = (uint8_t) (int8_t) this->o_del, e_del8 = (uint8_t) (int8_t) this->e_del;
    /* band-edge write: byte penalties, as the seed row's set1_epi8 / subs_epu8 uses them */
    const int oe_ins8 = (uint8_t) (this->o_ins + this->e_ins), e_ins8 = (uint8_t) this->e_ins;

    uint64_t live[BSW_COMPACT_GROUPS_MAX];
    uint8_t active[BSW_COMPACT_GROUPS_MAX];
    BswMove moves[BSW_COMPACT_GROUPS_MAX * SIMD_WIDTH8];

    auto writeback = [&](Bsw8CVec512 &V, int l) {
        SeqPair *p = V.pair[l];
        if (!p) return;
        const int xr = V.xrow16[l];
        int sc = V.h0_raw[l];   /* raw h0; the score is max(h0, maxScore) once the lane advanced */
        const int ms = (uint8_t) (V.maxScore[l] ^ hbias8);   /* maxScore is kept in the score domain */
        if (xr != 0 && ms > sc) sc = ms;
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
        for (int v = 0; v < nv; v++) {
            Bsw8CVec512 &V = Vs[v];
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
                V.xrow16[l] = 0; V.ierow16[l] = 0; V.h0_raw[l] = p[l].h0; V.gbest8[l] = 0; V.gset8[l] = 0;
                V.pair[l] = &p[l];
                V.head[l] = 0; V.tail[l] = 127; V.maxend[l] = (uint8_t) -128; V.exit0[l] = 0xFF;
                V.maxScore[l] = (uint8_t) (V.h0[l] ^ hbias8); V.y[l] = 0; V.max_off[l] = 0;
            }
            V.nbeg = 0; V.nend = vcol; V.ncol = vcol; V.minq = minq - 1;
            x86_soa_pack<SIMD_WIDTH8>(V.seq1, seq1p, len1a, len1a, nrow + 1, DUMMY1, DUMMY1, false, AMBIG, 8);
            x86_soa_pack<SIMD_WIDTH8>(V.seq2, seq2p, len2a, len2a, ncol + 1, DUMMY2, DUMMY2, true, AMBIG, 8);
            const __m512i h0_512 = _mm512_load_si512((const __m512i *) V.h0);
            /* H and F are kept in the score domain: seeds go in through the bias (hzero512 is the
             * XOR bias: 0x80 signed, 0 unsigned) */
            (SIGNED8 ? domain_counts_.signed_batches : domain_counts_.unsigned_batches)++;
            _mm512_store_si512((__m512i *) (V.H + ncol * SIMD_WIDTH8), hzero512);
            _mm512_store_si512((__m512i *) V.H, _mm512_xor_si512(h0_512, hzero512));
            __m512i t = _mm512_subs_epu8(h0_512, oe_ins512);
            _mm512_store_si512((__m512i *) (V.H + SIMD_WIDTH8), _mm512_xor_si512(t, hzero512));
            for (int k = 2; k < ncol; k++) {
                t = _mm512_subs_epu8(t, e_ins512);
                _mm512_store_si512((__m512i *) (V.H + k * SIMD_WIDTH8), _mm512_xor_si512(t, hzero512));
            }
            for (int j = 0; j <= ncol; j++) _mm512_store_si512((__m512i *) (V.F + j * SIMD_WIDTH8), hzero512);
            active[v] = 1;
        }

        /* ---- one row of one group: the row of smithWaterman512_8Impl ---- */
        auto row = [&](Bsw8CVec512 &V, const int i) -> uint64_t {
            int8_t *H_h = V.H, *F = V.F;
            const uint8_t *seq2SoA = V.seq2;
            const int minq = V.minq, vncol = V.ncol;
            __m512i head512     = _mm512_load_si512((const __m512i *) V.head);
            __m512i tail512     = _mm512_load_si512((const __m512i *) V.tail);
            __m512i maxend512   = _mm512_load_si512((const __m512i *) V.maxend);
            __m512i exit0       = _mm512_load_si512((const __m512i *) V.exit0);
            __m512i maxScore512 = _mm512_load_si512((const __m512i *) V.maxScore);
            __m512i y512        = _mm512_load_si512((const __m512i *) V.y);
            __m512i max_off512  = _mm512_load_si512((const __m512i *) V.max_off);
            const __m512i myband512 = _mm512_load_si512((const __m512i *) V.myband);
            const __m512i h0_512    = _mm512_load_si512((const __m512i *) V.h0);
            int j;

            __m512i e11 = hzero512;
            __m512i h00, h11, h10;
            __m512i s10 = _mm512_load_si512((const __m512i *) (V.seq1 + i * SIMD_WIDTH8));
            int beg = V.nbeg, end = V.nend;
            if (beg < i - w) beg = i - w;
            if (end > i + w + 1) end = i + w + 1;
            if (end > vncol) end = vncol;
            h10 = hzero512;
            if (beg == 0) {
                long sub = (long) o_del8 + (long) (i + 1) * e_del8;
                if (sub > 255) sub = 255;
                h10 = _mm512_xor_si512(_mm512_subs_epu8(h0_512, _mm512_set1_epi8((char) (uint8_t) sub)), hzero512);
            }
            __m512i j512 = zero512;
            __m512i maxRS1 = hzero512;
            __m512i y1_512 = zero512;
            __m512i hqe512 = hzero512, qfire512 = zero512;

            int8_t qlen_off_a[SIMD_WIDTH8]   __attribute((aligned(64)));
            int8_t qlen_valid_a[SIMD_WIDTH8] __attribute((aligned(64)));
            int8_t cmpim_a[SIMD_WIDTH8]      __attribute((aligned(64)));
            for (int l = 0; l < SIMD_WIDTH8; l++) {
                int qoff = V.qlenw[l] - i;
                int qoff_sat = qoff;
                if (qoff_sat >  127) qoff_sat =  127;
                if (qoff_sat < -128) qoff_sat = -128;
                qlen_off_a[l]   = (int8_t) qoff_sat;
                qlen_valid_a[l] = (qoff >= -w && qoff <= w + 1) ? (int8_t) 0xFF : 0;
                cmpim_a[l]      = ((i + 1) > V.mlenw[l]) ? (int8_t) 0xFF : 0;
            }
            __m512i qlen_off512 = _mm512_load_si512((__m512i *) qlen_off_a);
            __mmask64 qlen_valid_k = _mm512_movepi8_mask(_mm512_load_si512((__m512i *) qlen_valid_a));

            __m512i phead512 = head512, ptail512 = tail512;
            head512 = _mm512_max_epi8(head512, _mm512_sub_epi8(zero512, myband512));
            tail512 = _mm512_min_epi8(tail512, _mm512_add_epi8(myband512, one512));
            tail512 = _mm512_min_epi8(tail512, qlen_off512);
            {   /* band-edge write, as in smithWaterman512_8Impl */
                __mmask64 nm_ = _mm512_cmpeq_epi8_mask(_mm512_sub_epi8(tail512, maxend512), two512) & _mm512_movepi8_mask(exit0);
                if (nm_) {
                    int8_t me_a[SIMD_WIDTH8] __attribute((aligned(64)));
                    _mm512_store_si512((__m512i *) me_a, maxend512);
                    while (nm_) { const int l = __builtin_ctzll(nm_); nm_ &= nm_ - 1;
                        const int c = i + me_a[l] + 1;
                        H_h[c * SIMD_WIDTH8 + l] = (int8_t) ((uint8_t) bsw_seed_row(V.h0[l], oe_ins8, e_ins8, c) ^ hbias8);
                        F[c * SIMD_WIDTH8 + l] = (int8_t) hbias8; }
                }
                maxend512 = _mm512_max_epi8(maxend512, tail512);
            }
            {   /* trim, counted (the lean row; see smithWaterman512_8Impl) */
                __mmask64 cmph = _mm512_cmpeq_epi8_mask(head512, phead512) & _mm512_cmpeq_epi8_mask(tail512, ptail512);
                if (cmph != dmask) {
                    int trim_hi = i + bsw8_hmax_epi8_512(head512);
                    if (trim_hi > end) trim_hi = end;
                    for (int l = beg; l < trim_hi; l++) {
                        __m512i h512 = _mm512_load_si512((__m512i *)(H_h + l * SIMD_WIDTH8));
                        __m512i f512 = _mm512_load_si512((__m512i *)(F + l * SIMD_WIDTH8));
                        __m512i pj512 = _mm512_set1_epi8(l - i);
                        __mmask64 cmp1 = _mm512_cmpgt_epi8_mask(head512, pj512) | _mm512_cmpgt_epi8_mask(pj512, tail512);
                        h512 = _mm512_mask_blend_epi8(cmp1, h512, hzero512);
                        f512 = _mm512_mask_blend_epi8(cmp1, f512, hzero512);
                        _mm512_store_si512((__m512i *)(F + l * SIMD_WIDTH8), f512);
                        _mm512_store_si512((__m512i *)(H_h + l * SIMD_WIDTH8), h512);
                    }
                }
            }
            {
                __mmask64 cmpim = _mm512_movepi8_mask(_mm512_load_si512((__m512i *) cmpim_a));
                cmpim = cmpim | _mm512_cmpeq_epi8_mask(tail512, head512) | _mm512_cmpgt_epi8_mask(head512, tail512);
                exit0 = _mm512_mask_blend_epi8(cmpim, exit0, zero512);
            }
            if (CPG) {
                // CpG arm (see "CpG arm"): the row split, then the per-cell score
                __mmask64 rowcpg;
                BSW_CPG_ROW8_512(s10, rowcpg, cpgref512);
                for (int jp = beg; jp < end; jp++) {
                    __m512i s2 = _mm512_load_si512((__m512i *)(seq2SoA + jp * SIMD_WIDTH8));
                    __m512i sbt11;
                    SBT_PREPASS8_RANK1(s10, s2, rowcpg, sbt11, pmat512, cpgval512, cpgread512);
                    _mm512_store_si512((__m512i *)(sbt_buf + jp * SIMD_WIDTH8), sbt11);
                }
            } else if (!gen_mat) {
                for (int jp = beg; jp < end; jp++) {
                    __m512i s2 = _mm512_load_si512((__m512i *)(seq2SoA + jp * SIMD_WIDTH8));
                    __m512i sbt11;
                    SBT_PREPASS8_XOR(s10, s2, sbt11, pmat512);
                    _mm512_store_si512((__m512i *)(sbt_buf + jp * SIMD_WIDTH8), sbt11);
                }
            } else if (fc.rank1) {
                __mmask64 rowfreed = _mm512_cmpeq_epi8_mask(s10, frref512);
                for (int jp = beg; jp < end; jp++) {
                    __m512i s2 = _mm512_load_si512((__m512i *)(seq2SoA + jp * SIMD_WIDTH8));
                    __m512i sbt11;
                    SBT_PREPASS8_RANK1(s10, s2, rowfreed, sbt11, pmat512, match512, frread512);
                    _mm512_store_si512((__m512i *)(sbt_buf + jp * SIMD_WIDTH8), sbt11);
                }
            } else {
                for (int jp = beg; jp < end; jp++) {
                    __m512i s2 = _mm512_load_si512((__m512i *)(seq2SoA + jp * SIMD_WIDTH8));
                    __m512i sbt11;
                    SBT_PREPASS8_AMAT(s10, s2, sbt11, amat512, w_ambig_512, three512_8);
                    _mm512_store_si512((__m512i *)(sbt_buf + jp * SIMD_WIDTH8), sbt11);
                }
            }
            int fast_lo = beg, fast_hi = beg;
            if (_mm512_movepi8_mask(exit0) == dmask) {
                const int maxhead = bsw8_hmax_epi8_512(head512);
                const int mintail = bsw8_hmin_epi8_512(tail512);
                fast_lo = i + maxhead; if (fast_lo < beg) fast_lo = beg; if (fast_lo > end) fast_lo = end;
                fast_hi = i + mintail; if (fast_hi < fast_lo) fast_hi = fast_lo; if (fast_hi > end) fast_hi = end;
            }
            j512 = _mm512_set1_epi8(beg - i);
            __m512i nf_first512 = _mm512_set1_epi8(end - i);       /* see BSW8_NF_TRACK512 */
            __m512i nf_last512  = _mm512_set1_epi8(beg - 1 - i);
            /* the cell loops of smithWaterman512_8Impl, through its own macros */
            for (j = beg; j < fast_lo; j++)   EXT13_CELL8_512_MASKED
            for (j = fast_lo; j < fast_hi; j++) EXT13_CELL8_512_FAST
            for (j = fast_hi; j < end; j++)   EXT13_CELL8_512_MASKED
            {
                __mmask64 cmp1 = _mm512_cmpgt_epi8_mask(head512, j512) | _mm512_cmpgt_epi8_mask(j512, tail512);
                h10 = _mm512_mask_blend_epi8(cmp1, h10, hzero512);
                _mm512_store_si512((__m512i *)(H_h + j * SIMD_WIDTH8), h10);
                _mm512_store_si512((__m512i *)(F + j * SIMD_WIDTH8), hzero512);
                nf_last512 = _mm512_mask_mov_epi8(nf_last512, _mm512_cmpneq_epi8_mask(h10, hzero512), j512);   /* j512 == end - i */
            }
            /* the narrow side channels (BSW8_GSCORE_NARROW512 folds into these) */
            __m512i gbest512 = _mm512_load_si512((const __m512i *) V.gbest8);
            __mmask64 gset64 = _mm512_movepi8_mask(_mm512_load_si512((const __m512i *) V.gset8));
            __m512i ierow_lo = _mm512_load_si512((const __m512i *) V.ierow16);
            __m512i ierow_hi = _mm512_load_si512((const __m512i *) (V.ierow16 + 32));
            auto store_gscore = [&]() {
                _mm512_store_si512((__m512i *) V.gbest8, gbest512);
                _mm512_store_si512((__m512i *) V.gset8, _mm512_movm_epi8(gset64));
                _mm512_store_si512((__m512i *) V.ierow16, ierow_lo);
                _mm512_store_si512((__m512i *) (V.ierow16 + 32), ierow_hi);
            };
            __m512i bmaxScore512 = maxScore512;
            __mmask64 tmp = _mm512_cmpeq_epi8_mask(maxRS1, hzero512);
            if (tmp == dmask) {
                BSW8_GSCORE_NARROW512(gbest512, gset64, ierow_lo, ierow_hi,
                                      bsw8_unbias512<SIGNED8>(hqe512, x80_512),
                                      _mm512_movepi8_mask(qfire512), i + 1);
                store_gscore();
                _mm512_store_si512((__m512i *) V.exit0, zero512);
                return 0;
            }
            exit0 = _mm512_mask_blend_epi8(tmp, exit0, zero512);
            __m512i score512 = bsw8_max512<SIGNED8>(maxScore512, maxRS1);
            __mmask64 mex0 = _mm512_movepi8_mask(exit0);
            maxScore512 = _mm512_mask_blend_epi8(mex0, maxScore512, score512);
            __mmask64 cmp = bsw8_cmpgt512<SIGNED8>(maxScore512, bmaxScore512);
            y512 = _mm512_mask_blend_epi8(cmp, y512, y1_512);
            {
                __m512i ind512 = _mm512_abs_epi8(_mm512_sub_epi8(y1_512, one512));
                __m512i bmax_off512 = max_off512;
                ind512 = _mm512_max_epu8(max_off512, ind512);
                max_off512 = _mm512_mask_blend_epi8(cmp, bmax_off512, ind512);
            }
            {   /* lean epilogue: xrow16, gscore, gated z-drop (see smithWaterman512_8Impl) */
                const __m512i vip1_16 = _mm512_set1_epi16((int16_t) (i + 1));
                __m512i xrow_lo = _mm512_load_si512((const __m512i *) V.xrow16);
                __m512i xrow_hi = _mm512_load_si512((const __m512i *) (V.xrow16 + 32));
                xrow_lo = _mm512_mask_mov_epi16(xrow_lo, (__mmask32) cmp, vip1_16);
                xrow_hi = _mm512_mask_mov_epi16(xrow_hi, (__mmask32) (cmp >> 32), vip1_16);
                _mm512_store_si512((__m512i *) V.xrow16, xrow_lo);
                _mm512_store_si512((__m512i *) (V.xrow16 + 32), xrow_hi);
                BSW8_GSCORE_NARROW512(gbest512, gset64, ierow_lo, ierow_hi,
                                      bsw8_unbias512<SIGNED8>(hqe512, x80_512),
                                      _mm512_movepi8_mask(qfire512), i + 1);
                store_gscore();
                if (zdrop > 0) {   /* live lanes only, as in smithWaterman512_8Impl (mex0 == exit0 here) */
                    const __mmask64 need_z = _mm512_mask_cmpgt_epu8_mask(mex0,
                        _mm512_subs_epu8(bsw8_unbias512<SIGNED8>(maxScore512, x80_512),
                                         bsw8_unbias512<SIGNED8>(maxRS1, x80_512)), zgate512);
                    if (!zgate_exact || need_z != 0)
                        BSW8_ZDROP_WIDE512(xrow_lo, xrow_hi);
                }
            }
            /* as in smithWaterman512_8Impl: with no band columns maxRS1 stays 0 and the row returned above */
            xassert(beg < end, "AVX-512BW lane compaction: empty band reached the narrowing");
            {   /* narrowing from the tracked first / last offsets (see smithWaterman512_8Impl) */
                const int minfirst = bsw8_hmin_epi8_512(nf_first512);
                if (minfirst > beg - i) V.nbeg = i + minfirst - 1;
                const int l = i + bsw8_hmax_epi8_512(nf_last512);
                V.nend = l + 2 < vncol ? l + 2 : vncol;
                const __mmask64 live64 = _mm512_movepi8_mask(exit0);
                __mmask64 upd = live64 & _mm512_cmpneq_epi8_mask(nf_first512, _mm512_set1_epi8(beg - i));
                head512 = _mm512_mask_mov_epi8(head512, upd, nf_first512);
                upd = live64 & _mm512_cmpneq_epi8_mask(nf_last512, _mm512_set1_epi8(end - i));
                __m512i index512 = _mm512_mask_mov_epi8(tail512, upd, nf_last512);
                index512 = _mm512_add_epi8(index512, two512);
                tail512 = _mm512_min_epi8(index512, qlen_off512);
            }
            maxend512 = _mm512_subs_epi8(maxend512, one512);   /* frame shift, saturating */
            head512 = _mm512_sub_epi8(head512, one512);
            tail512 = _mm512_sub_epi8(tail512, one512);
#ifdef BSW8_ASSERT_ENVELOPE   /* debug only: see the byte-ceiling check in smithWaterman512_8Impl */
            assert(_mm512_cmpge_epu8_mask(bsw8_unbias512<SIGNED8>(maxRS1, x80_512), _mm512_set1_epi8((int8_t) BYTE_CEIL)) == 0 &&
                   "8-bit banded SW: row max hit the byte ceiling (lane compaction)");
#endif
            _mm512_store_si512((__m512i *) V.head, head512);
            _mm512_store_si512((__m512i *) V.tail, tail512);
            _mm512_store_si512((__m512i *) V.maxend, maxend512);
            _mm512_store_si512((__m512i *) V.exit0, exit0);
            _mm512_store_si512((__m512i *) V.maxScore, maxScore512);
            _mm512_store_si512((__m512i *) V.y, y512);
            _mm512_store_si512((__m512i *) V.max_off, max_off512);
            return (uint64_t) _mm512_movepi8_mask(exit0);
        };

        /* ---- move live lanes of one group into dead slots of another, between rows (next
         * row r) ----
         * As in bandedSWA_compact.inc, each moved lane takes its slot with it: the per-lane
         * trackers, its H, F and query bytes from its head to the superblock's last column, and
         * its target rows r to tlen - 1; the receiving slot is cleared left of the head, from the
         * group's lowest future column (the kernel needs score 0 there, in the active domain:
         * the byte 0x80 when signed, 0x00 unsigned). The k moves of a run share
         * their source and destination groups and touch distinct lanes, so the copies go a
         * column at a time for all of them, through this tier's cross-lane byte gather: within
         * a column each lane takes the source's byte, is cleared, or keeps its own, a choice that
         * only changes at the lanes' heads and at the cleared range's start, so the columns split
         * into at most k + 1 segments of constant masks; the target rows split at the lanes'
         * target lengths. */
        auto move_run = [&](Bsw8CVec512 &Sv, Bsw8CVec512 &D, const BswMove *mv, int k, int r) {
            /* k < SIMD_WIDTH8 in practice (a destination is an active group, so it has a live
             * lane); the breakpoint lists hold k + 1 */
            assert(k >= 1 && k <= SIMD_WIDTH8 && "move_run: a run moves 1 to SIMD_WIDTH8 lanes");
            uint8_t idx_a[SIMD_WIDTH8] __attribute((aligned(64)));   /* source lane per destination lane */
            uint8_t dl_a[SIMD_WIDTH8];                                /* the destination lanes */
            int ha_a[SIMD_WIDTH8], tl_a[SIMD_WIDTH8];
            memset(idx_a, 0, sizeof idx_a);
            for (int q = 0; q < k; q++) {
                const int sl = mv[q].sl, dl = mv[q].dl;
                writeback(D, dl);
                D.head[dl] = Sv.head[sl]; D.tail[dl] = Sv.tail[sl]; D.maxend[dl] = Sv.maxend[sl];
                D.exit0[dl] = Sv.exit0[sl]; D.maxScore[dl] = Sv.maxScore[sl]; D.y[dl] = Sv.y[sl];
                D.max_off[dl] = Sv.max_off[sl]; D.h0[dl] = Sv.h0[sl]; D.myband[dl] = Sv.myband[sl];
                D.xrow16[dl] = Sv.xrow16[sl]; D.ierow16[dl] = Sv.ierow16[sl]; D.h0_raw[dl] = Sv.h0_raw[sl];
                D.gbest8[dl] = Sv.gbest8[sl]; D.gset8[dl] = Sv.gset8[sl];
                D.qlenw[dl] = Sv.qlenw[sl]; D.mlenw[dl] = Sv.mlenw[sl]; D.tlenw[dl] = Sv.tlenw[sl];
                D.pair[dl] = Sv.pair[sl];
                int ha = r + (int8_t) Sv.head[sl];
                if (ha < 0) ha = 0;
                idx_a[dl] = (uint8_t) sl; dl_a[q] = (uint8_t) dl;
                ha_a[q] = ha; tl_a[q] = Sv.tlenw[sl];
                if (Sv.qlenw[sl] > D.ncol) D.ncol = Sv.qlenw[sl];
                if (Sv.qlenw[sl] - 1 < D.minq) D.minq = Sv.qlenw[sl] - 1;
            }
            if (Sv.nbeg < D.nbeg) D.nbeg = Sv.nbeg;
            if (Sv.nend > D.nend) D.nend = Sv.nend;
            const int cstart = D.nbeg < 0 ? 0 : D.nbeg;   /* the clear runs from here to each head */
            int bp[SIMD_WIDTH8 + 1], nb = 0;   /* column breakpoints: cstart and every head */
            bp[nb++] = cstart;
            for (int q = 0; q < k; q++) bp[nb++] = ha_a[q];
            std::sort(bp, bp + nb);
            nb = (int) (std::unique(bp, bp + nb) - bp);
            int rp[SIMD_WIDTH8 + 1], nr = 0;   /* row breakpoints: r and every target length past it */
            rp[nr++] = r;
            for (int q = 0; q < k; q++) if (tl_a[q] > r) rp[nr++] = tl_a[q];
            std::sort(rp, rp + nr);
            nr = (int) (std::unique(rp, rp + nr) - rp);
            uint8_t *dH = (uint8_t *) D.H, *dF = (uint8_t *) D.F, *dQ = D.seq2, *d1 = D.seq1;
            const uint8_t *sH = (const uint8_t *) Sv.H, *sF = (const uint8_t *) Sv.F, *sQ = Sv.seq2, *s1 = Sv.seq1;
            /* gather: dst byte l = src byte idx[l]. AVX-512BW has no full-width byte permute (that
             * is VBMI), so each 128-bit quarter of the source is broadcast and shuffled in-lane,
             * and the quarter each destination byte reads selects among the four. */
            const __m512i idx = _mm512_load_si512((const __m512i *) idx_a);
            __mmask64 qsel[4] = {0, 0, 0, 0};   /* destination lanes reading each source quarter */
            for (int q = 0; q < k; q++) qsel[mv[q].sl >> 4] |= (__mmask64) 1 << dl_a[q];
            /* quarter 0 is shuffled for every lane; the others only when some lane reads them */
            const bool q1 = qsel[1] != 0, q2 = qsel[2] != 0, q3 = qsel[3] != 0;
            auto gather = [&](const uint8_t *p) {
                const __m512i v = _mm512_load_si512((const __m512i *) p);
                __m512i t = _mm512_shuffle_epi8(_mm512_shuffle_i32x4(v, v, 0x00), idx);
                if (q1) t = _mm512_mask_shuffle_epi8(t, qsel[1], _mm512_shuffle_i32x4(v, v, 0x55), idx);
                if (q2) t = _mm512_mask_shuffle_epi8(t, qsel[2], _mm512_shuffle_i32x4(v, v, 0xAA), idx);
                if (q3) t = _mm512_mask_shuffle_epi8(t, qsel[3], _mm512_shuffle_i32x4(v, v, 0xFF), idx);
                return t;
            };
            /* the lane mask of the destination lanes of the moves q for which pred(q) holds */
            auto lanes_where = [&](auto pred) {
                __mmask64 m = 0;
                for (int q = 0; q < k; q++) if (pred(q)) m |= (__mmask64) 1 << dl_a[q];
                return m;
            };
            const __mmask64 moved = lanes_where([](int) { return true; });
            /* the columns stop at ncol, as the lane-by-lane copy did; a live lane's head is at
             * most its query length, so the clear never needs more */
            for (int j = 0; j < nb; j++) {
                const int c0 = bp[j], c1 = (j + 1 < nb) ? bp[j + 1] : ncol + 1;
                if (c0 > ncol) break;
                const __mmask64 take = lanes_where([&](int q) { return ha_a[q] <= c0; });   /* c >= head */
                const __mmask64 clear = (c0 >= cstart) ? (moved & ~take) : 0;
                for (int c = c0; c < c1 && c <= ncol; c++) {
                    const size_t o = (size_t) c * SIMD_WIDTH8;
                    /* a cleared lane holds score 0 in the score domain */
                    const __m512i h = _mm512_mask_mov_epi8(hzero512, take, gather(sH + o));
                    const __m512i f = _mm512_mask_mov_epi8(hzero512, take, gather(sF + o));
                    _mm512_store_si512((__m512i *) (dH + o), _mm512_mask_mov_epi8(
                                       _mm512_load_si512((const __m512i *) (dH + o)), take | clear, h));
                    _mm512_store_si512((__m512i *) (dF + o), _mm512_mask_mov_epi8(
                                       _mm512_load_si512((const __m512i *) (dF + o)), take | clear, f));
                    if (take)   /* the query row is only copied, never cleared */
                        _mm512_store_si512((__m512i *) (dQ + o), _mm512_mask_mov_epi8(
                                           _mm512_load_si512((const __m512i *) (dQ + o)), take, gather(sQ + o)));
                }
            }
            for (int j = 0; j + 1 < nr; j++) {   /* target rows r to tlen - 1 of each lane */
                const int r0 = rp[j];
                const __mmask64 take = lanes_where([&](int q) { return tl_a[q] > r0; });
                for (int rr = rp[j]; rr < rp[j + 1]; rr++) {
                    const size_t o = (size_t) rr * SIMD_WIDTH8;
                    _mm512_store_si512((__m512i *) (d1 + o), _mm512_mask_mov_epi8(
                                       _mm512_load_si512((const __m512i *) (d1 + o)), take, gather(s1 + o)));
                }
            }
            for (int q = 0; q < k; q++) { Sv.pair[mv[q].sl] = nullptr; Sv.exit0[mv[q].sl] = 0; }
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
                /* in plan order, a run of moves with one source and one destination group at a
                 * time (a later run may move lanes out of an earlier run's destination) */
                for (int k0 = 0; k0 < nm;) {
                    int k1 = k0 + 1;
                    while (k1 < nm && moves[k1].sv == moves[k0].sv && moves[k1].dv == moves[k0].dv) k1++;
                    move_run(Vs[moves[k0].sv], Vs[moves[k0].dv], moves + k0, k1 - k0, i + 1);
                    k0 = k1;
                }
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
// The score domain of a superblock (see the note above MAIN_CODE8_CORE) is decided here, once,
// from the gap penalty bytes; the instantiation seeds and clears H/F in that domain and runs
// the matching kernel body.
void BandedPairWiseSW::smithWatermanBatchWrapper8Compact(SeqPair *pairArray,
                                                         uint8_t *seqBufRef,
                                                         uint8_t *seqBufQer,
                                                         int32_t numPairs,
                                                         int32_t w,
                                                         int K)
{
    if (cpg_) {
        if (bsw8_gap_bytes_fit_int8(o_del, e_del, o_ins, e_ins))
            smithWatermanBatchWrapper8CompactImpl<true, true>(pairArray, seqBufRef, seqBufQer, numPairs, w, K);
        else
            smithWatermanBatchWrapper8CompactImpl<false, true>(pairArray, seqBufRef, seqBufQer, numPairs, w, K);
    } else if (bsw8_gap_bytes_fit_int8(o_del, e_del, o_ins, e_ins))
        smithWatermanBatchWrapper8CompactImpl<true, false>(pairArray, seqBufRef, seqBufQer, numPairs, w, K);
    else
        smithWatermanBatchWrapper8CompactImpl<false, false>(pairArray, seqBufRef, seqBufQer, numPairs, w, K);
}

#endif  /* __AVX512BW__ */
