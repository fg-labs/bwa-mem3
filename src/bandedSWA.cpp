/*************************************************************************************
                           The MIT License

   BWA-MEM2  (Sequence alignment using Burrows-Wheeler Transform),
   Copyright (C) 2019  Intel Corporation, Heng Li.

   Permission is hereby granted, free of charge, to any person obtaining
   a copy of this software and associated documentation files (the
   "Software"), to deal in the Software without restriction, including
   without limitation the rights to use, copy, modify, merge, publish,
   distribute, sublicense, and/or sell copies of the Software, and to
   permit persons to whom the Software is furnished to do so, subject to
   the following conditions:

   The above copyright notice and this permission notice shall be
   included in all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
   BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.

Authors: Vasimuddin Md <vasimuddin.md@intel.com>; Sanchit Misra <sanchit.misra@intel.com>;
*****************************************************************************************/

#include "kernel_dispatch.h"
#include "bandedSWA.h"
#ifdef VTUNE_ANALYSIS
#include <ittnotify.h> 
#endif

#if defined(__clang__) || defined(__GNUC__)
#define __mmask8 uint8_t
#define __mmask16 uint16_t
#define __mmask32 uint32_t
#endif

// ------------------------------------------------------------------------------------
// Sub-slice overshoot guard.
//
// The batched kernels (getScores8 / getScores16) round numPairs up to their
// SIMD width and INITIALIZE the padding lanes pairArray[numPairs .. roundUp)
// in place (len=0 dummy pairs). That is safe when the kernel owns the whole
// array, but a caller may pass a SUB-SLICE of a larger array -- the --meth
// per-hypothesis dispatch (bsw_run_tier) runs three CONTIGUOUS slices
// (OT | OB | SYM) of one pair array back-to-back. Padding writes past numPairs
// would zero the START of the next slice, and that slice's kernel call would
// then score those real pairs as empty (score == h0, gscore == -1) -> a
// spurious soft-clip / confident mismap on --meth reads. This RAII guard saves
// pairArray[numPairs .. roundUp) on construction and restores it on scope exit,
// so getScores{8,16} never leave a caller's pairs past numPairs modified. The
// overshoot is at most (width - 1) <= SIMD_WIDTH8 - 1 pairs; callers of a whole
// array over-allocate by MAX_LINE_LEN so the region is always addressable.
//
// The save/restore is only needed for a SUB-SLICE caller (the --meth OT/OB
// runs, whose overshoot lands in the next slice's live pairs). A whole-array
// caller -- every non-meth extension, plus the last (SYM) --meth slice -- has
// nothing but its own over-allocated padding past numPairs, so guarding it is
// pure overhead on the hot path. It is therefore gated on `active`: only the
// OT/OB kernel objects set it (via set_guard_overshoot); when false the guard
// is a no-op (nOver stays 0, both loops empty).
namespace {
struct BswOvershootGuard {
    SeqPair *pa;
    int base, nOver;
    SeqPair saved[SIMD_WIDTH8];
    BswOvershootGuard(SeqPair *pairArray, int numPairs, int width, bool active)
        : pa(pairArray), base(numPairs), nOver(0) {
        if (!active) return;
        int roundUp = ((numPairs + width - 1) / width) * width;
        nOver = roundUp - numPairs;
        for (int s = 0; s < nOver; s++) saved[s] = pa[base + s];
    }
    ~BswOvershootGuard() {
        for (int s = 0; s < nOver; s++) pa[base + s] = saved[s];
    }
    BswOvershootGuard(const BswOvershootGuard &) = delete;
    BswOvershootGuard &operator=(const BswOvershootGuard &) = delete;
};
}  // namespace

// ------------------------------------------------------------------------------------
// MACROs for vector code
extern uint64_t prof[10][112];
#define AMBIG 4
#define DUMMY1 99
#define DUMMY2 100

// Build the 16-byte LUT used by SBT_PREPASS8_LUT across all 8-bit kernel
// variants (AVX2, AVX-512BW, SSE2/NEON). Layout:
//   [0]      = w_match    (s1 == s2, both ACGT)
//   [1..3]   = w_mismatch (ACGT vs ACGT, XOR ∈ {1,2,3})
//   [4..15]  = w_ambig    (any cell with target=N or query=N; indices
//                           4..7 = target N=4, 8..11 = query N=8, 12 =
//                           both N, 13..15 unreachable but filled for
//                           safety so out-of-band XORs never read junk)
static inline void build_pmat16(int8_t out[16],
                                int8_t w_match, int8_t w_mismatch, int8_t w_ambig)
{
    out[0] = w_match;
    out[1] = w_mismatch;
    out[2] = w_mismatch;
    out[3] = w_mismatch;
    for (int i = 4; i < 16; i++) out[i] = w_ambig;
}

// Build the 16-byte ASYMMETRIC substitution LUT for the generic-matrix scoring
// seam (D3 bisulfite OT/OB matrices). Indexed by (ref<<2)|read over the ACGT
// (0..3) submatrix: out[(i<<2)|j] = mat[i*5 + j] (target-major mat[ref*5+read]).
// N cells (target/query == 4) are NOT in this LUT — callers mask them to ambig.
// For a symmetric matrix this reproduces build_pmat16's scores exactly (the
// diagonal is w_match, every off-diagonal is w_mismatch), so the default path
// stays bit-identical; an asymmetric off-diagonal (e.g. ref-C x read-T freed)
// is now honored where the XOR-LUT could not represent it.
static inline void build_amat16(int8_t out[16], const int8_t *mat)
{
    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++)
            out[(i << 2) | j] = mat[i * 5 + j];
}

// Defensive: numThreads must be ≥ 1 — _mm_malloc(0, 64) is implementation-
// defined and a 0-byte slab would silently OOM downstream writes. Wrappers
// fan out to constructors and per-batch SoA setup; clamping centrally keeps
// every entry point in agreement.
template <typename T>
static inline T effective_threads(T n) { return n < 1 ? T(1) : n; }

// ---------------------------------------------------------------------------
// 8-bit re-baseline floor helpers (shared by wrapper + every SIMD kernel tier)
//
// Initial 8-bit score floor: keep the seed score h0 inside the re-baseline
// keep-window so the stored byte (= H_abs - B) stays within the unsigned [0,255] range.
// Clamping prefixes >REBASE_KEEP below h0 is lossless (REBASE_KEEP=zdrop+1 >
// zdrop => already past the z-drop horizon). MUST be used by every SIMD tier.
// Defined at file scope (before any tier #if block) so every kernel variant
// (AVX2, AVX-512BW, SSE2/NEON) sees the same definition.
static inline int bsw8_rebase_keep(int zdrop)            { return zdrop + 1; }
static inline int bsw8_initial_floor(int h0, int zdrop) {
    int b0 = h0 - bsw8_rebase_keep(zdrop);
    return b0 > 0 ? b0 : 0;
}
// ---------------------------------------------------------------------------

//-----------------------------------------------------------------------------------
// constructor
BandedPairWiseSW::BandedPairWiseSW(const int o_del, const int e_del, const int o_ins,
                                   const int e_ins, const int zdrop,
                                   const int end_bonus, const int8_t *mat_,
                                   const int8_t w_match, const int8_t w_mismatch, int numThreads)
{
    mat = mat_;
    this->m = 5;
    this->end_bonus = end_bonus;
    this->zdrop = zdrop;
    this->o_del = o_del;
    this->o_ins = o_ins;
    this->e_del = e_del;
    this->e_ins = e_ins;
    
    this->w_match    = w_match;
    this->w_mismatch = w_mismatch*-1;
    this->w_open     = o_del;  // redundant, used in vector code.
    this->w_extend   = e_del;  // redundant, used in vector code.
    this->w_ambig    = DEFAULT_AMBIG;
    this->swTicks = 0;
    this->SW_cells = 0;
    setupTicks = 0;
    sort1Ticks = 0;
    swTicks = 0;
    sort2Ticks = 0;
    this->F8_ = this->H8_  = this->H8__ = NULL;
    this->F16_ = this->H16_  = this->H16__ = NULL;
    this->sbt8_ = NULL;
    this->sbt16_ = NULL;
    this->dp_slab_ = NULL;

    // Fold all per-thread scratch into one 64B-aligned slab:
    //   3 × sz8  : F8_, H8_, H8__       (8-bit DP rails)
    //   3 × sz16 : F16_, H16_, H16__    (16-bit DP rails)
    //   1 × sz8  : sbt8_                (8-bit SBT pre-pass scratch)
    //   1 × sz16 : sbt16_               (16-bit SBT pre-pass scratch;
    //                                    multi-MB — moved off the stack)
    // Each partition is a multiple of 64B (MAX_SEQ_LEN8=1088, SIMD_WIDTH8
    // ≥ 16 → 8-bit partition multiple of 17408; 16-bit similarly), so
    // neighbouring views stay 64B-aligned without explicit padding. The
    // static_asserts below trip at compile time if a future MAX_SEQ_LEN*
    // / SIMD_WIDTH* tweak breaks the 64B partition invariant.
    static_assert((MAX_SEQ_LEN8  * SIMD_WIDTH8  * sizeof(int8_t))  % 64 == 0,
                  "8-bit DP partition not a multiple of 64B; tune MAX_SEQ_LEN8 / SIMD_WIDTH8");
    static_assert((MAX_SEQ_LEN16 * SIMD_WIDTH16 * sizeof(int16_t)) % 64 == 0,
                  "16-bit DP partition not a multiple of 64B; tune MAX_SEQ_LEN16 / SIMD_WIDTH16");
    numThreads = effective_threads(numThreads);
    size_t sz8  = (size_t)MAX_SEQ_LEN8  * SIMD_WIDTH8  * numThreads * sizeof(int8_t);
    size_t sz16 = (size_t)MAX_SEQ_LEN16 * SIMD_WIDTH16 * numThreads * sizeof(int16_t);
    size_t total = 4 * sz8 + 4 * sz16;
    dp_slab_ = _mm_malloc(total, 64);
    if (UNLIKELY(dp_slab_ == NULL)) {
        fprintf(stderr, "BSW DP-scratch slab allocation failed (%zu bytes)\n", total);
        exit(EXIT_FAILURE);
    }

    char *base = (char *)dp_slab_;
    F8_   = (int8_t  *)(base + 0       );
    H8_   = (int8_t  *)(base + sz8     );
    H8__  = (int8_t  *)(base + 2 * sz8 );
    sbt8_ = (int8_t  *)(base + 3 * sz8 );
    F16_  = (int16_t *)(base + 4 * sz8);
    H16_  = (int16_t *)(base + 4 * sz8 + sz16);
    H16__ = (int16_t *)(base + 4 * sz8 + 2 * sz16);
    sbt16_ = (int16_t *)(base + 4 * sz8 + 3 * sz16);
}

// destructor
BandedPairWiseSW::~BandedPairWiseSW() {
    _mm_free(dp_slab_);
    F8_ = H8_ = H8__ = NULL;
    F16_ = H16_ = H16__ = NULL;
    sbt8_ = NULL;
    sbt16_ = NULL;
    dp_slab_ = NULL;
}

int64_t BandedPairWiseSW::getTicks()
{
    //printf("oneCount = %ld, totalCount = %ld\n", oneCount, totalCount);
    int64_t totalTicks = sort1Ticks + setupTicks + swTicks + sort2Ticks;
    printf("cost breakup: %lld, %lld, %lld, %lld, %lld\n",
            (long long)sort1Ticks, (long long)setupTicks, (long long)swTicks,
            (long long)sort2Ticks, (long long)totalTicks);

    return totalTicks;
}
// ------------------------------------------------------------------------------------
// Banded SWA - scalar code
// ------------------------------------------------------------------------------------

int BandedPairWiseSW::scalarBandedSWA(int qlen, const uint8_t *query,
                                      int tlen, const uint8_t *target,
                                      int32_t w, int h0, int *_qle, int *_tle,
                                      int *_gtle, int *_gscore,
                                      int *_max_off) {
    
    // uint64_t sw_cells = 0;
    eh_t *eh; // score array
    int8_t *qp; // query profile
    int i, j, k, oe_del = o_del + e_del, oe_ins = o_ins + e_ins, beg, end, max, max_i, max_j, max_ins, max_del, max_ie, gscore, max_off;
    
    // assert(h0 > 0); //check !!!
    
    // allocate memory
    qp = (int8_t *) malloc(qlen * m);
    assert(qp != NULL);
    eh = (eh_t *) calloc(qlen + 1, 8);
    assert(eh != NULL);

    // generate the query profile
    for (k = i = 0; k < m; ++k) {
        const int8_t *p = &mat[k * m];
        //for (j = 0; j < qlen; ++j) qp[i++] = p[query[j]-48];  //sub 48
        for (j = 0; j < qlen; ++j) qp[i++] = p[query[j]];
    }

    // fill the first row
    eh[0].h = h0; eh[1].h = h0 > oe_ins? h0 - oe_ins : 0;
    for (j = 2; j <= qlen && eh[j-1].h > e_ins; ++j)
        eh[j].h = eh[j-1].h - e_ins;

    // adjust $w if it is too large
    k = m * m;
    for (i = 0, max = 0; i < k; ++i) // get the max score
        max = max > mat[i]? max : mat[i];
    max_ins = (int)((double)(qlen * max + end_bonus - o_ins) / e_ins + 1.);
    max_ins = max_ins > 1? max_ins : 1;
    w = w < max_ins? w : max_ins;
    max_del = (int)((double)(qlen * max + end_bonus - o_del) / e_del + 1.);
    max_del = max_del > 1? max_del : 1;
    w = w < max_del? w : max_del; // TODO: is this necessary?

    // DP loop
    max = h0, max_i = max_j = -1; max_ie = -1, gscore = -1;
    max_off = 0;
    beg = 0, end = qlen;
    for (i = 0; (i < tlen); ++i) {
        int t, f = 0, h1, m = 0, mj = -1;
        //int8_t *q = &qp[(target[i]-48) * qlen];   // sub 48
        int8_t *q = &qp[(target[i]) * qlen];
        // apply the band and the constraint (if provided)
        if (beg < i - w) beg = i - w;
        if (end > i + w + 1) end = i + w + 1;
        if (end > qlen) end = qlen;
        // compute the first column
        if (beg == 0) {
            h1 = h0 - (o_del + e_del * (i + 1));
            if (h1 < 0) h1 = 0;
        } else h1 = 0;
        for (j = beg; (j < end); ++j) {
            // At the beginning of the loop: eh[j] = { H(i-1,j-1), E(i,j) }, f = F(i,j) and h1 = H(i,j-1)
            // Similar to SSE2-SW, cells are computed in the following order:
            //   H(i,j)   = max{H(i-1,j-1)+S(i,j), E(i,j), F(i,j)}
            //   E(i+1,j) = max{H(i,j)-gapo, E(i,j)} - gape
            //   F(i,j+1) = max{H(i,j)-gapo, F(i,j)} - gape
            eh_t *p = &eh[j];
            int h, M = p->h, e = p->e; // get H(i-1,j-1) and E(i-1,j)
            p->h = h1;          // set H(i,j-1) for the next row
            M = M? M + q[j] : 0;// separating H and M to disallow a cigar like "100M3I3D20M"
            h = M > e? M : e;   // e and f are guaranteed to be non-negative, so h>=0 even if M<0
            h = h > f? h : f;
            h1 = h;             // save H(i,j) to h1 for the next column
            mj = m > h? mj : j; // record the position where max score is achieved
            m = m > h? m : h;   // m is stored at eh[mj+1]
            // gap-from-H (standard Gotoh): gap-open is subtracted from the cell
            // max H(i,j)=h, not the diagonal-only M. This matches ksw2/minimap2
            // and is the convention the int8 delta recurrence requires. The final
            // reported CIGAR/score still come from ksw_global2 (gap-from-M), so
            // alignment output is unchanged; this only governs the extension
            // boundary scorer. See FINDINGS-delta-recurrence.md (task #16).
            t = h - oe_del;
            t = t > 0? t : 0;
            e -= e_del;
            e = e > t? e : t;   // computed E(i+1,j)
            p->e = e;           // save E(i+1,j) for the next row
            t = h - oe_ins;
            t = t > 0? t : 0;
            f -= e_ins;
            f = f > t? f : t;   // computed F(i,j+1)
            // SW_cells++;
        }
        eh[end].h = h1; eh[end].e = 0;
        if (j == qlen) {
            max_ie = gscore > h1? max_ie : i;
            gscore = gscore > h1? gscore : h1;
        }
        if (m == 0) break;
        if (m > max) {
            max = m, max_i = i, max_j = mj;
            max_off = max_off > abs(mj - i)? max_off : abs(mj - i);
        } else if (zdrop > 0) {
            if (i - max_i > mj - max_j) {
                if (max - m - ((i - max_i) - (mj - max_j)) * e_del > zdrop) break;
            } else {
                if (max - m - ((mj - max_j) - (i - max_i)) * e_ins > zdrop) break;
            }
        }
        // update beg and end for the next round
        for (j = beg; (j < end) && eh[j].h == 0 && eh[j].e == 0; ++j);
        beg = j;
        for (j = end; (j >= beg) && eh[j].h == 0 && eh[j].e == 0; --j);
        end = j + 2 < qlen? j + 2 : qlen;
        //beg = 0; end = qlen; // uncomment this line for debugging
    }
    free(eh); free(qp);
    if (_qle) *_qle = max_j + 1;
    if (_tle) *_tle = max_i + 1;
    if (_gtle) *_gtle = max_ie + 1;
    if (_gscore) *_gscore = gscore;
    if (_max_off) *_max_off = max_off;
    
#if MAXI
    fprintf(stderr, "%d (%d %d) %d %d %d\n", max, max_i+1, max_j+1, gscore, max_off, max_ie+1);
#endif

    // return sw_cells;
    return max;
}

// -------------------------------------------------------------
// Banded SWA, wrapper function
//-------------------------------------------------------------
void BandedPairWiseSW::scalarBandedSWAWrapper(SeqPair *seqPairArray,
                                              uint8_t *seqBufRef,
                                              uint8_t *seqBufQer,
                                              int numPairs,
                                              int nthreads,
                                              int32_t w) {

    for (int i=0; i<numPairs; i++)
    {
        SeqPair *p = seqPairArray + i;
        uint8_t *seq1 = seqBufRef + p->idr;
        uint8_t *seq2 = seqBufQer + p->idq;
        
        p->score = scalarBandedSWA(p->len2, seq2, p->len1,
                                   seq1, w, p->h0, &p->qle, &p->tle,
                                   &p->gtle, &p->gscore, &p->max_off);      
    }

}


#if ((!__AVX512BW__) & (__AVX2__))

//------------------------------------------------------------------------------
// MACROs
// ------------------------ vec-8 ---------------------------------------------



// ------------------------ vec 16 --------------------------------------------------
#define _mm256_blendv_epi16(a,b,c)              \
        _mm256_blendv_epi8(a, b, c);            


#define ZSCORE16(i4_256, y4_256)                                            \
    {                                                                   \
        __m256i tmpi = _mm256_sub_epi16(i4_256, x256);                  \
        __m256i tmpj = _mm256_sub_epi16(y4_256, y256);                  \
        cmp = _mm256_cmpgt_epi16(tmpi, tmpj);                           \
        score256 = _mm256_sub_epi16(maxScore256, maxRS1);               \
        __m256i insdel = _mm256_blendv_epi16(e_ins256, e_del256, cmp);  \
        __m256i sub_a256 = _mm256_sub_epi16(tmpi, tmpj);                    \
        __m256i sub_b256 = _mm256_sub_epi16(tmpj, tmpi);                    \
        tmp = _mm256_blendv_epi16(sub_b256, sub_a256, cmp);             \
        tmp = _mm256_sub_epi16(score256, tmp);                          \
        cmp = _mm256_cmpgt_epi16(tmp, zdrop256);                            \
        exit0 = _mm256_andnot_si256(cmp, exit0);               \
    }



// --- PR 17/16: AVX2 LUT primitive ---
// pmat256 must have the 16-byte LUT broadcast into both 128-bit halves
// (use _mm256_broadcastsi128_si256). 1 xor + 1 shuffle replaces the
// legacy 4-op cmpeq+blendv+max+blendv path. The XOR index is order-insensitive,
// so it cannot represent an asymmetric matrix — that case uses SBT_PREPASS8_AMAT.
#define SBT_PREPASS8_XOR(s1, s2, sbt11_out, pmat256)                    \
    {                                                                   \
        __m256i xor_ = _mm256_xor_si256(s1, s2);                        \
        sbt11_out = _mm256_shuffle_epi8(pmat256, xor_);                 \
    }

// D3 generic-matrix seam (gated on an asymmetric matrix; the default symmetric
// path uses SBT_PREPASS8_XOR above). Index the target-major LUT
// amat[(ref<<2)|read] (broadcast into both 128-bit halves), then mask N cells to
// ambig. s1=target/ref, s2=query/read; ACGT=0..3, N=4(target)/8(query). For a
// symmetric matrix amat == the XOR pmat scores; an asymmetric off-diagonal
// (ref-C x read-T freed) is scored correctly. (s1<<2 via two add_epi8; N detected
// by max_epu8(s1,s2) > 3.) Mirrors the validated NEON SBT_PREPASS8_AMAT.
#define SBT_PREPASS8_AMAT(s1, s2, sbt11_out, amat256, ambig256, three256) \
    {                                                                   \
        __m256i sh_  = _mm256_add_epi8(s1, s1);                         \
        sh_          = _mm256_add_epi8(sh_, sh_);          /* s1 << 2 */ \
        __m256i idx_ = _mm256_or_si256(sh_, s2);           /* (ref<<2)|read */ \
        __m256i acgt_  = _mm256_shuffle_epi8(amat256, idx_);           \
        __m256i nmax_  = _mm256_max_epu8(s1, s2);                       \
        __m256i nmask_ = _mm256_cmpgt_epi8(nmax_, three256); /* N: base > 3 */ \
        sbt11_out = _mm256_or_si256(_mm256_andnot_si256(nmask_, acgt_), \
                                    _mm256_and_si256(nmask_, ambig256)); \
    }

// D3 8-bit rank-1 fast path: symmetric XOR LUT + single freed-to-match override.
// Cheaper than AMAT on the fast 8-bit tier (no index build, no N-mask). Mirrors NEON.
// rowfreed = (ref==fr_ref) hoisted per row (s1 is constant across the band loop).
#define SBT_PREPASS8_RANK1(s1, s2, rowfreed, sbt11_out, pmat256, match256, frread256) \
    {                                                                   \
        __m256i xor_  = _mm256_xor_si256(s1, s2);                       \
        __m256i sbt_  = _mm256_shuffle_epi8(pmat256, xor_);            \
        __m256i freed_ = _mm256_and_si256(rowfreed, _mm256_cmpeq_epi8(s2, frread256)); \
        sbt11_out = _mm256_blendv_epi8(sbt_, match256, freed_);        \
    }

#define MAIN_CODE8_CORE(sbt11, h00, h11, e11, f11, f21, zero256, e_ins256, oe_ins256, e_del256, oe_del256) \
    {                                                                   \
        /* M = max(0, h00 + sbt) computed in UNSIGNED-saturating form so a    \
         * legitimate score in [128,255] is kept (the old signed add_epi8 +   \
         * signed max_epi8 floor wrapped/mis-read >127 as negative). Split    \
         * the signed substitution score into its +bonus and -penalty parts:  \
         * adds_epu8 (no wrap past 255) then subs_epu8 (floors at 0). Mirrors  \
         * the validated NEON smithWaterman128_8 recurrence. */                \
        __m256i sbt_pos = _mm256_max_epi8(sbt11, zero256);              \
        __m256i sbt_neg = _mm256_max_epi8(_mm256_sub_epi8(zero256, sbt11), zero256); \
        __m256i m11 = _mm256_subs_epu8(_mm256_adds_epu8(h00, sbt_pos), sbt_neg); \
        __m256i cmp11 = _mm256_cmpeq_epi8(h00, zero256);                \
        m11 = _mm256_andnot_si256(cmp11, m11);  /* h00==0 -> local restart */ \
        h11 = _mm256_max_epu8(m11, e11);                                \
        h11 = _mm256_max_epu8(h11, f11);                                \
        /* gap-open from H (standard Gotoh), unsigned-saturating */      \
        __m256i temp256 = _mm256_subs_epu8(h11, oe_ins256);            \
        e11 = _mm256_subs_epu8(e11, e_ins256);                          \
        e11 = _mm256_max_epu8(temp256, e11);                            \
        temp256 = _mm256_subs_epu8(h11, oe_del256);                    \
        f21 = _mm256_subs_epu8(f11, e_del256);                          \
        f21 = _mm256_max_epu8(temp256, f21);                            \
    }

// Symmetric (default) fast path. N detected by the high bit of max_epu16 (N=0xFFFF).
#define SBT_PREPASS16_SYM(s1, s2, sbt11_out, mismatch256, match256, w_ambig_256) \
    {                                                                   \
        __m256i cmp_ = _mm256_cmpeq_epi16(s1, s2);                      \
        __m256i sbt_ = _mm256_blendv_epi16(mismatch256, match256, cmp_); \
        __m256i tmp_ = _mm256_max_epu16(s1, s2);                        \
        sbt11_out = _mm256_blendv_epi16(sbt_, w_ambig_256, tmp_);       \
    }

// D3 rank-1 fast path: symmetric + a single off-diagonal freed to a match
// (bisulfite OT/OB). Match when (ref==read) OR (ref==fr_ref AND read==fr_read);
// no LUT, no sign-extend. N(high bit of max_epu16)->ambig. Mirrors NEON.
// alt1 = (ref==fr_ref) ? fr_read : ref, hoisted per row (see NEON variant).
#define SBT_PREPASS16_RANK1(s1, s2, alt1, sbt11_out, mismatch256, match256, w_ambig_256) \
    {                                                                   \
        __m256i eq_  = _mm256_cmpeq_epi16(s2, s1);                      \
        __m256i fr_  = _mm256_cmpeq_epi16(s2, alt1);                    \
        __m256i ism_ = _mm256_or_si256(eq_, fr_);                       \
        __m256i sbt_ = _mm256_blendv_epi16(mismatch256, match256, ism_); \
        __m256i nmax_ = _mm256_max_epu16(s1, s2);                       \
        sbt11_out = _mm256_blendv_epi16(sbt_, w_ambig_256, nmax_);      \
    }

// D3 generic-matrix seam (gated; >=2 freed cells / non-match freed values). LUT
// amat[(ref<<2)|read] for both-ACGT lanes; non-ACGT is mismatch-or-ambig (padding
// pairs are never equal; N->ambig). Byte LUT score in the low byte -> sign-extend.
// acgt mask via cmpeq(max_epu16(maxb,3),3) (unsigned <=3 test). Mirrors NEON.
#define SBT_PREPASS16_AMAT(s1, s2, sbt11_out, amat256, mismatch256, w_ambig_256, three256) \
    {                                                                   \
        __m256i maxb_ = _mm256_max_epu16(s1, s2);                       \
        __m256i base_ = _mm256_blendv_epi16(mismatch256, w_ambig_256, maxb_); /* N high bit */ \
        __m256i acgt_ = _mm256_cmpeq_epi16(_mm256_max_epu16(maxb_, three256), three256); \
        __m256i idx_  = _mm256_or_si256(_mm256_slli_epi16(s1, 2), s2);  \
        __m256i lut_  = _mm256_shuffle_epi8(amat256, idx_);            \
        lut_ = _mm256_srai_epi16(_mm256_slli_epi16(lut_, 8), 8);       \
        sbt11_out = _mm256_blendv_epi16(base_, lut_, acgt_);          \
    }

#define MAIN_CODE16_CORE(sbt11, h00, h11, e11, f11, f21, zero256, e_ins256, oe_ins256, e_del256, oe_del256) \
    {                                                                   \
        __m256i m11 = _mm256_add_epi16(h00, sbt11);                     \
        __m256i cmp11 = _mm256_cmpeq_epi16(h00, zero256);               \
        m11 = _mm256_andnot_si256(cmp11, m11);                 \
        h11 = _mm256_max_epi16(m11, e11);                               \
        h11 = _mm256_max_epi16(h11, f11);                               \
        __m256i temp256 = _mm256_sub_epi16(h11, oe_ins256);             \
        __m256i val256  = _mm256_max_epi16(temp256, zero256);           \
        e11 = _mm256_sub_epi16(e11, e_ins256);                          \
        e11 = _mm256_max_epi16(val256, e11);                            \
        temp256 = _mm256_sub_epi16(h11, oe_del256);                     \
        val256  = _mm256_max_epi16(temp256, zero256);                   \
        f21 = _mm256_sub_epi16(f11, e_del256);                          \
        f21 = _mm256_max_epi16(val256, f21);                            \
    }

// MACROs section ends
// ------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------
// B-SWA - Vector code
// ------------------------- AVX2 - 8 bit SIMD_LANES ---------------------------

inline void sortPairsId(SeqPair *pairArray, int32_t first, int32_t count,
                        SeqPair *tempArray)
{

    int32_t i;

    for(i = 0; i < count; i++)
    {
        SeqPair sp = pairArray[i];
        int32_t pos = sp.id - first;
        tempArray[pos] = sp;
    }

    for(i = 0; i < count; i++)
        pairArray[i] = tempArray[i];    
}

/******************* Vector code, version 2.0 *************************/
#define PFD 2
void BandedPairWiseSW::getScores8(SeqPair *pairArray,
                                  uint8_t *seqBufRef,
                                  uint8_t *seqBufQer,
                                  int32_t numPairs,
                                  uint16_t numThreads,
                                  int32_t w)
{
    int64_t startTick, endTick;

    {
        BswOvershootGuard _g(pairArray, numPairs, SIMD_WIDTH8, guard_overshoot_);
        smithWatermanBatchWrapper8(pairArray, seqBufRef, seqBufQer, numPairs, numThreads, w);
    }

#if MAXI
    printf("AVX2 Vecor code: Writing output..\n");
    for (int l=0; l<numPairs; l++)
    {
        fprintf(stderr, "%d (%d %d) %d %d %d\n",
                pairArray[l].score, pairArray[l].tle, pairArray[l].qle,
                pairArray[l].gscore, pairArray[l].max_off, pairArray[l].gtle);
    }
    printf("Vector code: Writing output completed!!!\n\n");
#endif
    
}

void BandedPairWiseSW::smithWatermanBatchWrapper8(SeqPair *pairArray,
                                                  uint8_t *seqBufRef,
                                                  uint8_t *seqBufQer,
                                                  int32_t numPairs,
                                                  uint16_t numThreads,
                                                  int32_t w)
{
    numThreads = effective_threads(numThreads);
    int64_t st1, st2, st3, st4, st5;
#if RDT
    st1 = ___rdtsc();
#endif

    uint8_t *seq1SoA = NULL;
    seq1SoA = (uint8_t *)_mm_malloc((size_t)MAX_SEQ_LEN8 * SIMD_WIDTH8 * numThreads * sizeof(uint8_t), 64);

    uint8_t *seq2SoA = NULL;
    seq2SoA = (uint8_t *)_mm_malloc((size_t)MAX_SEQ_LEN8 * SIMD_WIDTH8 * numThreads * sizeof(uint8_t), 64);
    
    if (UNLIKELY(seq1SoA == NULL || seq2SoA == NULL)) {
        fprintf(stderr, "Error! Mem not allocated!!!\n");
        exit(EXIT_FAILURE);
    }
    
    int32_t ii;
    int32_t roundNumPairs = ((numPairs + SIMD_WIDTH8 - 1)/SIMD_WIDTH8 ) * SIMD_WIDTH8;
    // assert(roundNumPairs < BATCH_SIZE * SEEDS_PER_READ);
    for(ii = numPairs; ii < roundNumPairs; ii++)
    {
        pairArray[ii].id = ii;
        pairArray[ii].len1 = 0;
        pairArray[ii].len2 = 0;
    }

#if RDT
    st2 = ___rdtsc();
#endif
    

#if RDT
    st3 = ___rdtsc();
#endif
    
    int eb = end_bonus;
//#pragma omp parallel num_threads(numThreads)
    {
        int32_t i;
        // uint16_t tid = omp_get_thread_num();
        uint16_t tid = 0;
        uint8_t *mySeq1SoA = NULL;
        mySeq1SoA = seq1SoA + tid * MAX_SEQ_LEN8 * SIMD_WIDTH8;

        uint8_t *mySeq2SoA = NULL;
        mySeq2SoA = seq2SoA + tid * MAX_SEQ_LEN8 * SIMD_WIDTH8;
        assert(mySeq1SoA != NULL && mySeq2SoA != NULL);
        
        uint8_t *seq1;
        uint8_t *seq2;
        uint8_t h0[SIMD_WIDTH8]   __attribute__((aligned(64)));
        uint8_t band[SIMD_WIDTH8];      
        uint8_t qlen[SIMD_WIDTH8] __attribute__((aligned(64)));
        int32_t bsize = 0;
        
        int8_t *H1 = H8_ + tid * SIMD_WIDTH8 * MAX_SEQ_LEN8;
        int8_t *H2 = H8__ + tid * SIMD_WIDTH8 * MAX_SEQ_LEN8;
        
        __m256i zero256   = _mm256_setzero_si256();
        __m256i e_ins256  = _mm256_set1_epi8(e_ins);
        __m256i oe_ins256 = _mm256_set1_epi8(o_ins + e_ins);
        __m256i o_del256  = _mm256_set1_epi8(o_del);
        __m256i e_del256  = _mm256_set1_epi8(e_del);
        __m256i eb_ins256 = _mm256_set1_epi8(eb - o_ins);
        __m256i eb_del256 = _mm256_set1_epi8(eb - o_del);

        int8_t max = 0;
        if (max < w_match) max = w_match;
        if (max < w_mismatch) max = w_mismatch;
        if (max < w_ambig) max = w_ambig;

        // Initial per-lane score floor B0 = max(0, h0 - REBASE_KEEP_W) (long-read
        // 8-bit, seed score h0 >= 128). Seeding the H bytes from raw h0 overflows
        // the signed-int8 score state on row 0, so we seed from the rebaselined
        // h0' = min(h0, REBASE_KEEP_W) instead (see the per-lane h0 init below).
        // REBASE_KEEP_W MUST equal the kernel's REBASE_KEEP (bsw8_rebase_keep)
        // so the two agree on B0; smithWaterman256_8 recomputes B0 via
        // bsw8_initial_floor(p[].h0, zdrop).
        const int REBASE_KEEP_W = bsw8_rebase_keep(zdrop);
        // The h0-prefix column/row seed below is unsigned-saturating [0,255], so
        // the seeded byte min(h0, REBASE_KEEP_W) only needs to fit a uint8. The
        // routing gate (bsw8_envelope_ok) enforces zdrop + maxStep <= 253 (so the
        // re-baseline window REBASE_HI = 255 - maxStep > REBASE_KEEP holds);
        // assert the byte bound for any direct caller that bypasses the gate.
        assert(REBASE_KEEP_W <= 255 &&
               "8-bit h0-prefix seed byte min(h0,REBASE_KEEP) must fit uint8 (zdrop<=254)");

        int nstart = 0, nend = numPairs;


//#pragma omp for schedule(dynamic, 128)
        for(i = nstart; i < nend; i+=SIMD_WIDTH8)
        {
            int32_t j, k;
            int maxLen1 = 0;
            int maxLen2 = 0;
            bsize = w;
            
            uint64_t tim;
            for(j = 0; j < SIMD_WIDTH8; j++)
            {
                if ((i + j + PFD) < roundNumPairs) { // prefetch block (bounded; see getScores8/16 contract)
                    SeqPair spf = pairArray[i + j + PFD];
                    _mm_prefetch((const char*) seqBufRef + (int64_t)spf.idr, _MM_HINT_NTA);
                    _mm_prefetch((const char*) seqBufRef + (int64_t)spf.idr + 64, _MM_HINT_NTA);
                }

                SeqPair sp = pairArray[i + j];
                // Re-baseline the seed score into the int8 byte frame: seed the H
                // arrays from h0' = h0 - B0 = min(h0, REBASE_KEEP_W) (clamped >= 0)
                // so the initial bytes fit signed-int8 even when the true seed h0
                // is several hundred (long-read). smithWaterman256_8 recomputes the
                // matching B0 = max(0, h0 - REBASE_KEEP) from p[].h0 + zdrop and
                // starts that lane's floor B there, so byte == H_absolute - B.
                {
                    int h0p = sp.h0;
                    if (h0p > REBASE_KEEP_W) h0p = REBASE_KEEP_W;
                    if (h0p < 0) h0p = 0;
                    // Always-on uint8 guard (asserts are compiled out in release):
                    // a direct caller that bypasses the gate with zdrop > 254 would
                    // otherwise truncate the seed byte. In-envelope REBASE_KEEP_W <= 253,
                    // so this never fires on the production path.
                    if (h0p > 255) h0p = 255;
                    h0[j] = (uint8_t) h0p;
                }
                seq1 = seqBufRef + (int64_t)sp.idr;

                for(k = 0; k < sp.len1; k++)
                {
                    mySeq1SoA[k * SIMD_WIDTH8 + j] = seq1[k] /* PR16: N stays 4 */;
                }
                qlen[j] = sp.len2 * max;
                if(maxLen1 < sp.len1) maxLen1 = sp.len1;
            }

            for(j = 0; j < SIMD_WIDTH8; j++)
            {
                SeqPair sp = pairArray[i + j];
                for(k = sp.len1; k <= maxLen1; k++) //removed "="
                {
                    mySeq1SoA[k * SIMD_WIDTH8 + j] = DUMMY1;
                }
            }
            /* B5: only the boundary row H2[maxLen1] survives the h0-prefix
             * deletion seed below (which overwrites rows [0, maxLen1)); write
             * just that row here, before the seed, instead of the dead per-row
             * fills removed above. */
            _mm256_store_si256((__m256i *)(H2 + maxLen1 * SIMD_WIDTH8), _mm256_set1_epi8((char)DUMMY1));

//--------------------
            __m256i h0_256 = _mm256_load_si256((__m256i*) h0);
            _mm256_store_si256((__m256i *) H2, h0_256);
            // h0-prefix deletion seed, unsigned-saturating [0,255]. Was signed
            // sub_epi8 + max_epi8(.,0) floor, which required the seed byte
            // min(h0,REBASE_KEEP) <= 127 and so capped zdrop at 126. subs_epu8
            // floors at 0 inherently and feeds the floored value forward;
            // byte-identical for the monotone-decreasing affine prefix. Mirrors
            // smithWaterman128_8 (the NEON reference, already unsigned here).
            __m256i tmp256 = _mm256_subs_epu8(h0_256, o_del256);

            for(k = 1; k < maxLen1; k++) {
                tmp256 = _mm256_subs_epu8(tmp256, e_del256);
                _mm256_store_si256((__m256i *)(H2 + k* SIMD_WIDTH8), tmp256);
            }
//-------------------
            for(j = 0; j < SIMD_WIDTH8; j++)
            {
                if ((i + j + PFD) < roundNumPairs) { // prefetch block (bounded; see getScores8/16 contract)
                    SeqPair spf = pairArray[i + j + PFD];
                    _mm_prefetch((const char*) seqBufQer + (int64_t)spf.idq, _MM_HINT_NTA);
                    _mm_prefetch((const char*) seqBufQer + (int64_t)spf.idq + 64, _MM_HINT_NTA);
                }

                SeqPair sp = pairArray[i + j];
                seq2 = seqBufQer + (int64_t)sp.idq;
                
                if (sp.len2 > MAX_SEQ_LEN8) fprintf(stderr, "Error !! : %d %d\n", sp.id, sp.len2);
                assert(sp.len2 < MAX_SEQ_LEN8);
                
                for(k = 0; k < sp.len2; k++)
                {
                    mySeq2SoA[k * SIMD_WIDTH8 + j] = (seq2[k]==AMBIG ? 8 : seq2[k]) /* PR16: query N→8 */;
                }
                if(maxLen2 < sp.len2) maxLen2 = sp.len2;
            }

            for(j = 0; j < SIMD_WIDTH8; j++)
            {
                SeqPair sp = pairArray[i + j];
                for(k = sp.len2; k <= maxLen2; k++)
                {
                    mySeq2SoA[k * SIMD_WIDTH8 + j] = DUMMY2;
                }
            }
            /* B5: only the boundary row H1[maxLen2] (value 0) survives the
             * h0-prefix insertion seed below; write just that row, before the
             * seed so its unconditional H1[0]/H1[1] stores still win. */
            _mm256_store_si256((__m256i *)(H1 + maxLen2 * SIMD_WIDTH8), _mm256_setzero_si256());

//------------------------
            _mm256_store_si256((__m256i *) H1, h0_256);
            // h0-prefix insertion seed, unsigned-saturating [0,255]:
            // H1[1] = max(0, h0' - oe_ins), then -e_ins per step. subs_epu8
            // replaces the signed cmpgt+sub+blendv (first gap) and sub+max_epi8
            // (extensions); byte-identical for h0' in [0,255].
            tmp256 = _mm256_subs_epu8(h0_256, oe_ins256);
            _mm256_store_si256((__m256i *) (H1 + SIMD_WIDTH8), tmp256);
            for(k = 2; k < maxLen2; k++)
            {
                tmp256 = _mm256_subs_epu8(tmp256, e_ins256);
                _mm256_store_si256((__m256i *)(H1 + k*SIMD_WIDTH8), tmp256);
            }
//------------------------
            /* Banding calculation in pre-processing */
            uint8_t myband[SIMD_WIDTH8] __attribute__((aligned(64)));
            uint8_t temp[SIMD_WIDTH8] __attribute__((aligned(64)));
            {
                __m256i qlen256 = _mm256_load_si256((__m256i *) qlen);
                __m256i sum256 = _mm256_add_epi8(qlen256, eb_ins256);
                _mm256_store_si256((__m256i *) temp, sum256);               
                for (int l=0; l<SIMD_WIDTH8; l++)
                {
                    double val = temp[l]/e_ins + 1.0;
                    int max_ins = (int) val;
                    max_ins = max_ins > 1? max_ins : 1;
                    myband[l] = min_(bsize, max_ins);
                }
                sum256 = _mm256_add_epi8(qlen256, eb_del256);
                _mm256_store_si256((__m256i *) temp, sum256);               
                for (int l=0; l<SIMD_WIDTH8; l++) {
                    double val = temp[l]/e_del + 1.0;
                    int max_ins = (int) val;
                    max_ins = max_ins > 1? max_ins : 1;
                    myband[l] = min_(myband[l], max_ins);
                    bsize = bsize < myband[l] ? myband[l] : bsize;
                }
            }
            
            smithWaterman256_8(mySeq1SoA,
                               mySeq2SoA,
                               maxLen1,
                               maxLen2,
                               pairArray + i,
                               h0,
                               tid,
                               numPairs,
                               zdrop,
                               bsize,
                               myband);
        }
    }

#if RDT
    st4 = ___rdtsc();
#endif
    

#if RDT
    st5 = ___rdtsc();
    setupTicks = st2 - st1;
    sort1Ticks = st3 - st2;
    swTicks = st4 - st3;
    sort2Ticks = st5 - st4;
#endif
    
    // free mem
    _mm_free(seq1SoA);
    _mm_free(seq2SoA);
    
    return;
}


void BandedPairWiseSW::smithWaterman256_8(uint8_t seq1SoA[],
                                          uint8_t seq2SoA[],
                                          int nrow,
                                          int ncol,
                                          SeqPair *p,
                                          uint8_t h0[],
                                          uint16_t tid,
                                          int32_t numPairs,
                                          int zdrop,
                                          int32_t w,
                                          uint8_t myband[])
{
    __m256i match256     = _mm256_set1_epi8(this->w_match);
    __m256i mismatch256  = _mm256_set1_epi8(this->w_mismatch);
    __m256i w_ambig_256  = _mm256_set1_epi8(this->w_ambig); // ambig penalty

    // PR 16: pmat LUT, broadcast into both 128-bit halves (shuffle_epi8 is
    // lane-wise on AVX2 — each half shuffles against its own half of pmat256).
    int8_t pmat_bytes[16] __attribute__((aligned(16)));
    build_pmat16(pmat_bytes, this->w_match, this->w_mismatch, this->w_ambig);
    __m128i pmat128 = _mm_load_si128((__m128i *)pmat_bytes);
    __m256i pmat256 = _mm256_broadcastsi128_si256(pmat128);
    // D3 generic-matrix seam: symmetric default uses the XOR pmat (fast); an
    // asymmetric matrix (bisulfite OT/OB) uses the target-major amat LUT.
    const bool forced  = bsw_force_generic_matrix();
    const bool gen_mat = bsw_generic_matrix(this->mat, this->w_match, this->w_mismatch)
                         || forced;
    const BswFreedCell fc = bsw_freed_cell(this->mat, this->w_match, this->w_mismatch, forced);
    __m256i frref256  = _mm256_set1_epi8(fc.ref);
    __m256i frread256 = _mm256_set1_epi8(fc.read);
    int8_t amat_bytes[16] __attribute__((aligned(16)));
    build_amat16(amat_bytes, this->mat);
    __m256i amat256 = _mm256_broadcastsi128_si256(_mm_load_si128((__m128i *)amat_bytes));
    __m256i three256_8 = _mm256_set1_epi8(3);

    __m256i e_del256    = _mm256_set1_epi8(this->e_del);
    __m256i oe_del256   = _mm256_set1_epi8(this->o_del + this->e_del);
    __m256i e_ins256    = _mm256_set1_epi8(this->e_ins);
    __m256i oe_ins256   = _mm256_set1_epi8(this->o_ins + this->e_ins);
    
    int8_t  *F  = F8_ + tid * SIMD_WIDTH8 * MAX_SEQ_LEN8;
    int8_t  *H_h    = H8_ + tid * SIMD_WIDTH8 * MAX_SEQ_LEN8;
    int8_t  *H_v = H8__ + tid * SIMD_WIDTH8 * MAX_SEQ_LEN8;

    int i, j;

    uint8_t tail[SIMD_WIDTH8] __attribute((aligned(64)));
    uint8_t head[SIMD_WIDTH8] __attribute((aligned(64)));

    // PR 17: per-row score-vector scratch for fission. Thread-local view
    // into sbt8_ (slab-backed; see constructor) — too large for the stack.
    int8_t *sbt_buf = sbt8_ + tid * SIMD_WIDTH8 * MAX_SEQ_LEN8;

    // --- DIAGONAL-OFFSET POSITION ENCODING (long-read 8-bit, w<=127) ---
    // Every per-cell COLUMN position is tracked as the diagonal offset
    //   d = col - i  in [-w, +w+1]  (fits signed int8 for w <= ~126).
    // ROW quantities (best row, best-gscore row, qlen, tlen, mlen) exceed
    // int8 for long reads, so they live in WIDE per-lane int32 side channels
    // updated O(rows) in the per-row epilogue, not O(cells). Absolute end
    // coordinates are reconstructed at the result store from the wide row.
    // Persistent column-offset state (head256/tail256) is shifted by -1 each
    // row (frame follows i) so the same absolute edge keeps its offset.
    int32_t tlenw[SIMD_WIDTH8];   // raw target length (rows), wide
    int32_t qlenw[SIMD_WIDTH8];   // raw query length (cols), wide
    int32_t mbandw[SIMD_WIDTH8];  // per-lane band width, wide
    int32_t mlenw[SIMD_WIDTH8];   // min(qlen+myband, tlen), wide row bound
    int32_t xrow[SIMD_WIDTH8];    // best row for score (== i+1 at capture)
    int32_t ierow[SIMD_WIDTH8];   // best row for gscore (== i+1 at capture)

    // --- PER-LANE SCORE RE-BASELINING (long-read 8-bit, scores > 127) ---
    // H/F/E and the score maxima are unsigned 8-bit [0,255], but a 1 kb
    // alignment scores ~900. To keep the kernel 32-lane (8-bit) we carry a
    // wide per-lane running FLOOR B[l]: every stored 8-bit score represents
    // (H_absolute - B[l]). When a lane's row-max climbs toward 255 we raise
    // B[l] by `delta` and subtract `delta` (saturating, _mm256_subs_epu8) from
    // every carried score buffer/register for that lane, then add B back when
    // reconstructing the absolute result.
    //
    // Scores live in the UNSIGNED byte range [0,255]. The DP recurrence computes
    // M = max(0, h00 + sbt) with unsigned-saturating arithmetic (adds_epu8 then
    // subs_epu8) and the row/global-max trackers (maxRS1, maxScore256) compare
    // with UNSIGNED order (== after _mm256_max_epu8), so the full [0,255] is usable.
    //
    // CORRECTNESS depends on re-baselining NEVER firing for a pair admitted to
    // this kernel. The saturating-subtract is NOT generally lossless: it can zero
    // a still-positive off-diagonal cell (a cell may sit > zdrop below the ROW max
    // yet still lie on the eventual optimum — z-drop is a row-level early-exit, not
    // a per-cell guarantee), and a zeroed cell is then misread as the h00==0
    // local-restart sentinel, spuriously truncating a valid alignment. So we do
    // NOT rely on re-baseline being lossless; instead it is held inert:
    //   REBASE_HI = 255 - maxStep: a row whose max has not reached REBASE_HI grows
    //     by <= maxStep per row (H(i,j) <= H(i-1,j-1) + maxStep), so its bytes stay
    //     <= 255 with no rebase needed.
    //   The bwamem.cpp routing gate admits a pair only when its MAX ATTAINABLE
    //     score stays < REBASE_HI and h0 <= REBASE_KEEP (so B0 == 0). Under that
    //     gate the row max never reaches REBASE_HI: B stays 0, stored bytes equal
    //     absolute scores, and the kernel is a plain exact unsigned [0,255] SW.
    // The re-baseline machinery below is retained only as a safety net; pairs that
    // could exceed the envelope take the 16-bit path.
    // Precondition REBASE_HI > REBASE_KEEP (a non-empty window) holds for
    // zdrop + maxStep <= 253, which is exactly what the routing gate enforces.
    // The h0-prefix column/row seed (wrapper setup below) is unsigned-saturating
    // [0,255] and imposes no tighter ceiling; it previously used signed int8 ops
    // that required the seed byte <= 127 and capped zdrop at 126.
    int maxStep = (int)this->w_match;
    if ((int)this->w_ambig > maxStep) maxStep = (int)this->w_ambig;
    if (maxStep < 1) maxStep = 1;
    const int REBASE_KEEP = zdrop + 1;
    const int REBASE_HI   = 255 - maxStep;
    assert(REBASE_HI > REBASE_KEEP &&
           "8-bit SW re-baseline: zdrop+maxStep>253, scoring outside safe envelope (route to 16-bit)");
    int32_t B[SIMD_WIDTH8];        // per-lane wide score floor (stored = abs - B)
    int32_t best_abs[SIMD_WIDTH8]; // running best score in ABSOLUTE units
    int32_t gbest_abs[SIMD_WIDTH8];// running gscore (query-end) in ABSOLUTE units

    // --- INITIAL SCORE FLOOR B0 (long-read 8-bit, seed score h0 >= 128) ---
    // The H arrays are seeded from the seed score h0 (H_v[0]=h0, then gap-
    // decremented down column 0 and across row 0). For a long-read seed h0 can
    // be several hundred, which would overflow the unsigned [0,255] byte state
    // before any re-baseline event can fire. So we start each lane's
    // running floor at B0[l] = max(0, h0 - REBASE_KEEP) instead of 0: the
    // wrapper seeds the H bytes from h0' = h0 - B0 = min(h0, REBASE_KEEP)
    // (clamped >= 0, so the initial bytes land in [0, REBASE_KEEP] <= 254), and
    // here we set B[l] = B0[l] so the stored bytes still mean (H_absolute - B).
    //
    // SAFETY (identical argument to the mid-DP re-baseline): the gap-prefix
    // cells in column 0 / row 0 that fall more than REBASE_KEEP below h0 are
    // saturated to byte 0 by the wrapper's unsigned-saturating seed. Because
    // REBASE_KEEP = zdrop+1 > zdrop, any such prefix cell is already beyond the
    // z-drop horizon below the seed, so an alignment routed through it has
    // already z-dropped and cannot be optimal -- clamping it to the floor is
    // lossless. This is the SAME situation that arises mid-DP whenever B>0 after
    // a re-baseline (byte 0 then means H_absolute == B, not 0), which the
    // h00==0 local-restart test in MAIN_CODE8_CORE already handles correctly; B0
    // simply makes B>0 from row 0 instead of from the first re-baseline event.
    int32_t B0[SIMD_WIDTH8];
    for (int l=0; l<SIMD_WIDTH8; l++) {
        B0[l] = bsw8_initial_floor((int)p[l].h0, zdrop);
    }

    int32_t minq = 10000000;
    for (int l=0; l<SIMD_WIDTH8; l++) {
        tlenw[l]  = p[l].len1;
        qlenw[l]  = p[l].len2;
        mbandw[l] = myband[l];
        int ml = qlenw[l] + mbandw[l];
        if (ml > tlenw[l]) ml = tlenw[l];
        mlenw[l]  = ml;
        xrow[l]   = 0;
        ierow[l]  = 0;
        B[l]        = B0[l];   // start the floor at B0 so byte state = abs - B0
        best_abs[l] = p[l].h0; // maxScore256 inits to the h0 seed; record the
                               // ABSOLUTE seed score (wide), not the rebaselined byte
        gbest_abs[l]= -1;      // unset sentinel (-1): gscore=-1 / gtle=0 when no query end is reached, matching scalar
        if (p[l].len2 < minq) minq = p[l].len2;
    }
    minq -= 1; // for gscore

    __m256i myband256 = _mm256_load_si256((__m256i *) myband);
    __m256i zero256 = _mm256_setzero_si256();
    __m256i one256  = _mm256_set1_epi8(1);
    __m256i two256  = _mm256_set1_epi8(2);
    __m256i ff256 = _mm256_set1_epi8(0xFF);

    // gscore query-end capture (see smithWaterman128_8). Reset per row.
    __m256i hqe256   = zero256;   // query-end cell H (rebaselined byte) this row
    __m256i qfire256 = zero256;   // 0xFF where this lane reached its query end this row

    // Offset-frame band edges. head_off starts at 0 (col 0 - row 0); tail_off
    // starts saturated-high (+127) and is immediately clamped by the band-grow
    // min() against (1+myband) and (qlen-i) on the first row.
    __m256i head256 = zero256;
    __m256i tail256 = _mm256_set1_epi8(127);
    _mm256_store_si256((__m256i *) head, head256);
    _mm256_store_si256((__m256i *) tail, tail256);

    __m256i hval = _mm256_load_si256((__m256i *)(H_v));
    __mmask32 dmask = 0xFFFFFFFF;

    __m256i maxScore256 = hval;
    for(j = 0; j < ncol; j++)
        _mm256_store_si256((__m256i *)(F + j * SIMD_WIDTH8), zero256);

    __m256i y256 = zero256;   // best col as diagonal offset: col - i_at_capture (i_at_capture = xrow[l]-1)
    __m256i max_off256 = zero256;
    __m256i exit0 = _mm256_set1_epi8(0xFF);
    __m256i zdrop256 = _mm256_set1_epi8(zdrop);

    int beg = 0, end = ncol;
    int nbeg = beg, nend = end;
    
#if RDT
    uint64_t tim = __rdtsc();
#endif
    
    for(i = 0; i < nrow; i++)
    {
        __m256i e11 = zero256;
        __m256i h00, h11, h10;
        __m256i s10 = _mm256_load_si256((__m256i *)(seq1SoA + (i + 0) * SIMD_WIDTH8));

        beg = nbeg; end = nend;
        // Banding
        if (beg < i - w) beg = i - w;
        if (end > i + w + 1) end = i + w + 1;
        if (end > ncol) end = ncol;

        h10 = zero256;
        if (beg == 0)
            h10 = _mm256_load_si256((__m256i *)(H_v + (i+1) * SIMD_WIDTH8));

        __m256i j256 = zero256;
        __m256i maxRS1 = zero256;

        __m256i y1_256 = zero256;   // row-max column as diagonal offset (col - i)

        // gscore query-end capture resets each row.
        hqe256   = zero256;
        qfire256 = zero256;

        // Per-row diagonal-offset of the query end: qlen_off = qlen - i. Built
        // wide then saturated to int8, with a validity mask so an out-of-band
        // qlen-i (which would wrap and spuriously alias an in-band col offset)
        // never triggers a false gscore/clamp. qlen_off_valid is true only when
        // qlen-i lies within the representable band window [-w, w+1].
        int8_t qlen_off_a[SIMD_WIDTH8]   __attribute((aligned(32)));
        int8_t qlen_valid_a[SIMD_WIDTH8] __attribute((aligned(32)));
        int8_t cmpim_a[SIMD_WIDTH8]      __attribute((aligned(32)));
        for (int l=0; l<SIMD_WIDTH8; l++) {
            int qoff = qlenw[l] - i;                 // wide
            int qoff_sat = qoff;                     // saturate to int8 range
            if (qoff_sat >  127) qoff_sat =  127;
            if (qoff_sat < -128) qoff_sat = -128;
            qlen_off_a[l]   = (int8_t) qoff_sat;
            qlen_valid_a[l] = (qoff >= -w && qoff <= w + 1) ? (int8_t)0xFF : 0;
            // exit when row i+1 has passed the per-lane effective length bound.
            cmpim_a[l]      = ((i + 1) > mlenw[l]) ? (int8_t)0xFF : 0;
        }
        __m256i qlen_off256   = _mm256_load_si256((__m256i *) qlen_off_a);
        __m256i qlen_valid256 = _mm256_load_si256((__m256i *) qlen_valid_a);

#if RDT
        uint64_t tim1 = __rdtsc();
#endif

        // Banding (diagonal-offset frame). head256/tail256 arrive here in row-i's
        // offset frame: the -1 shift at the end of the previous iteration converts
        // (col - (i-1)) -> (col - i), so no additional adjustment is needed here.
        //   abs head-grow: head = max(head, i - myband)  ->  head_off = max(head_off, -myband)
        //   abs tail-clamp: tail = min(tail, (i+1)+myband, qlen)
        //                                          -> tail_off = min(tail_off, 1+myband, qlen-i)
        __m256i cache256;
        __m256i phead256 = head256, ptail256 = tail256;
        __m256i negband256 = _mm256_sub_epi8(zero256, myband256);            // -myband
        head256 = _mm256_max_epi8(head256, negband256);
        cache256 = _mm256_add_epi8(myband256, one256);                       // 1 + myband
        tail256 = _mm256_min_epi8(tail256, cache256);
        tail256 = _mm256_min_epi8(tail256, qlen_off256);

        // NEW, trimming.
        __m256i cmph = _mm256_cmpeq_epi8(head256, phead256);
        __m256i cmpt = _mm256_cmpeq_epi8(tail256, ptail256);
        // cmph &= cmpt;
        cmph = _mm256_and_si256(cmph, cmpt);
        __mmask32 cmp_ht = _mm256_movemask_epi8(cmph);

        for (int l=beg; l<end && cmp_ht != dmask; l++)
        {
            __m256i h256 = _mm256_load_si256((__m256i *)(H_h + l * SIMD_WIDTH8));
            __m256i f256 = _mm256_load_si256((__m256i *)(F + l * SIMD_WIDTH8));

            __m256i pj256 = _mm256_set1_epi8(l - i);   // diagonal offset of column l
            __m256i cmp1 = _mm256_cmpgt_epi8(head256, pj256);
            uint32_t cval = _mm256_movemask_epi8(cmp1);
            if (cval == 0x00) break;
            __m256i cmp2 = _mm256_cmpgt_epi8(pj256, tail256);
            cmp1 = _mm256_or_si256(cmp1, cmp2);
            h256 = _mm256_andnot_si256(cmp1, h256);
            f256 = _mm256_andnot_si256(cmp1, f256);

            _mm256_store_si256((__m256i *)(F + l * SIMD_WIDTH8), f256);
            _mm256_store_si256((__m256i *)(H_h + l * SIMD_WIDTH8), h256);
        }

#if RDT
        prof[DP3][0] += __rdtsc() - tim1;
#endif

        // cmpim: lane exits if row i+1 passed its effective length (precomputed
        // wide into cmpim_a) OR the band collapsed (tail<=head, offset frame).
        __m256i cmpim = _mm256_load_si256((__m256i *) cmpim_a);
        __m256i cmpht = _mm256_cmpeq_epi8(tail256, head256);
        cmpim = _mm256_or_si256(cmpim, cmpht);
        // NEW
        cmpht = _mm256_cmpgt_epi8(head256, tail256);
        cmpim = _mm256_or_si256(cmpim, cmpht);

        exit0 = _mm256_andnot_si256(cmpim, exit0);

#if RDT
        tim1 = __rdtsc();
#endif

        // PR 17: pre-pass — build score vector sbt[j] for all band columns
        // once. Independent per-j, so vectorized loads/stores run at peak
        // throughput without carrying the DP core's dep chain. gen_mat branch is
        // loop-invariant (once per row): symmetric keeps the XOR LUT.
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

        j256 = _mm256_set1_epi8(beg - i);   // diagonal offset of first band column
#pragma unroll(4)
        for(j = beg; j < end; j++)
        {
            __m256i f11, f21, sbt11;
            h00 = _mm256_load_si256((__m256i *)(H_h + j * SIMD_WIDTH8));
            f11 = _mm256_load_si256((__m256i *)(F + j * SIMD_WIDTH8));
            sbt11 = _mm256_load_si256((__m256i *)(sbt_buf + j * SIMD_WIDTH8));

            __m256i pj256 = j256;
            j256 = _mm256_add_epi8(j256, one256);

            MAIN_CODE8_CORE(sbt11, h00, h11, e11, f11, f21, zero256,
                            e_ins256, oe_ins256,
                            e_del256, oe_del256);

            // Masked writing
            __m256i cmp1 = _mm256_cmpgt_epi8(head256, pj256);
            __m256i cmp2 = _mm256_cmpgt_epi8(pj256, tail256);
            cmp1 = _mm256_or_si256(cmp1, cmp2);
            h10 = _mm256_andnot_si256(cmp1, h10);
            f21 = _mm256_andnot_si256(cmp1, f21);

            __m256i bmaxRS = maxRS1;
            maxRS1 =_mm256_max_epu8(maxRS1, h11);
            // "new row-max" argmax mask. maxRS1 = max(bmaxRS,h11), so
            // (maxRS1 != bmaxRS) is a strict subset of (maxRS1 == h11) (both
            // mean h11 >= bmaxRS); the OR was redundant. cmpeq(maxRS1,h11) is
            // the exact combined mask — bit-identical, drops a cmpeq+xor+or.
            __m256i cmpA = _mm256_cmpeq_epi8(maxRS1, h11);
            cmp1 = _mm256_cmpgt_epi8(j256, tail256);
            cmp1 = _mm256_or_si256(cmp1, cmp2);
            cmpA = _mm256_blendv_epi8(y1_256, j256, cmpA);
            y1_256 = _mm256_blendv_epi8(cmpA, y1_256, cmp1);
            maxRS1 = _mm256_blendv_epi8(maxRS1, bmaxRS, cmp1);

            _mm256_store_si256((__m256i *)(F + j * SIMD_WIDTH8), f21);
            _mm256_store_si256((__m256i *)(H_h + j * SIMD_WIDTH8), h10);

            h10 = h11;

            // gscore query-end capture (see smithWaterman128_8 for the full
            // rationale: the re-baseline saturating-subtract zeroes off-diagonal
            // query-end cells, so gscore cannot be reconstructed from the trimmed
            // tail; capture (byte+B) per row and finalize wide). Fire exactly like
            // the byte-identical 16-bit tier: col == qlen-1 (j256 == qlen_off) AND
            // the band-grown tail reached the query end (tail256 == qlen_off ==
            // scalar's end == qlen) AND in-band (qlen_valid) AND lane alive (exit0).
            // gtle CONTRACT (see smithWaterman128_8): exact vs scalar for
            // gscore > 0; may differ only in the unused gscore == 0 tail.
            if (j >= minq)
            {
                __m256i cmp = _mm256_cmpeq_epi8(j256, qlen_off256);
                cmp = _mm256_and_si256(cmp, _mm256_cmpeq_epi8(tail256, qlen_off256));
                cmp = _mm256_and_si256(cmp, qlen_valid256);
                cmp = _mm256_and_si256(cmp, exit0);
                hqe256   = _mm256_blendv_epi8(hqe256, h11, cmp);
                qfire256 = _mm256_blendv_epi8(qfire256, ff256, cmp);
            }
        }
        __m256i cmp1 = _mm256_cmpgt_epi8(head256, j256);
        __m256i cmp2 = _mm256_cmpgt_epi8(j256, tail256);
        cmp1 = _mm256_or_si256(cmp1, cmp2);
        h10 = _mm256_andnot_si256(cmp1, h10);

        _mm256_store_si256((__m256i *)(H_h + j * SIMD_WIDTH8), h10);
        _mm256_store_si256((__m256i *)(F + j * SIMD_WIDTH8), zero256);


        /* exit due to zero score by a row */
        uint32_t cval = 0;
        __m256i bmaxScore256 = maxScore256;
        __m256i tmp = _mm256_cmpeq_epi8(maxRS1, zero256);
        cval = _mm256_movemask_epi8(tmp);
        if (cval == 0xFFFFFFFF) {
            /* Finalize THIS row's query-end (gscore/gtle) capture before exiting
             * (see smithWaterman128_8 zero-row break for the rationale): scalar
             * records the row's gscore then hits its own m==0 break, so the vector
             * must too, else the gscore==0 query-end tail row is dropped -> wrong
             * gtle under asymmetric --meth matrices. Mirror of the epilogue gscore
             * block (>= tie-break: latest query-end row wins). */
            int8_t qf_a_[SIMD_WIDTH8] __attribute((aligned(32)));
            int8_t hq_a_[SIMD_WIDTH8] __attribute((aligned(32)));
            _mm256_store_si256((__m256i *) qf_a_, qfire256);
            _mm256_store_si256((__m256i *) hq_a_, hqe256);
            for (int g = 0; g < SIMD_WIDTH8 / 8; g++) {
                const int base = g * 8;
                __m256i Bg   = _mm256_loadu_si256((const __m256i *)(B + base));
                __m256i hqeg = _mm256_cvtepu8_epi32(_mm_loadl_epi64((const __m128i *)(hq_a_ + base)));
                __m256i qfg  = _mm256_cvtepi8_epi32(_mm_loadl_epi64((const __m128i *)(qf_a_ + base)));
                __m256i gs_abs = _mm256_add_epi32(hqeg, Bg);
                __m256i gba = _mm256_loadu_si256((const __m256i *)(gbest_abs + base));
                __m256i ge  = _mm256_xor_si256(_mm256_cmpgt_epi32(gba, gs_abs), ff256);
                __m256i gmask = _mm256_and_si256(qfg, ge);
                gba = _mm256_blendv_epi8(gba, gs_abs, gmask);
                _mm256_storeu_si256((__m256i *)(gbest_abs + base), gba);
                __m256i ierg = _mm256_loadu_si256((const __m256i *)(ierow + base));
                ierg = _mm256_blendv_epi8(ierg, _mm256_set1_epi32(i + 1), gmask);
                _mm256_storeu_si256((__m256i *)(ierow + base), ierg);
            }
            break;
        }

        exit0 = _mm256_andnot_si256(tmp, exit0);

        __m256i score256 = _mm256_max_epu8(maxScore256, maxRS1);
        maxScore256 = _mm256_blendv_epi8(maxScore256, score256, exit0);

        // UNSIGNED >: maxScore256 (post-update) >= bmaxScore256, so (>u) == (!=).
        // Signed cmpgt_epi8 mis-read scores >127.
        __m256i cmp = _mm256_xor_si256(_mm256_cmpeq_epi8(maxScore256, bmaxScore256), ff256);
        // y256 (best col) stays a diagonal offset captured in the best row's
        // frame; the best row itself moves to the wide xrow[] side channel.
        y256 = _mm256_blendv_epi8(y256, y1_256, cmp);

        // max_off = max running diagonal-distance of the row-max from the main
        // diagonal: |y1col - (i+1)| = |y1_off - 1| in the offset frame.
        // y1_off - 1 is a SIGNED int8 (negative when row-max is left of the
        // sub-diagonal), so use _mm256_abs_epi8.
        __m256i y1_minus1 = _mm256_sub_epi8(y1_256, one256);
        tmp = _mm256_abs_epi8(y1_minus1);                    // |y1_off - 1|
        __m256i bmax_off256 = max_off256;
        tmp = _mm256_max_epu8(max_off256, tmp);
        max_off256 = _mm256_blendv_epi8(bmax_off256, tmp, cmp);

        // Per-lane wide updates (O(rows)): best-score row (xrow), best-gscore
        // row (ierow), absolute-frame score tracking, and the z-drop test — all
        // done in wide scalars so row distances that exceed int8 for long reads
        // are handled exactly.
        {
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
            // VECTORIZED per-lane wide updates (see smithWaterman128_8). 8 int32
            // per __m256i, so SIMD_WIDTH8 lanes -> SIMD_WIDTH8/8 groups of 8. Each
            // group's z-drop die mask is narrowed (lane-crossing-free) to 8 bytes
            // via a 128-bit pack and cleared in exit_a; exit0 reloaded at the end.
            const __m256i vi   = _mm256_set1_epi32(i);
            const __m256i vip1 = _mm256_set1_epi32(i + 1);
            const __m256i vone = _mm256_set1_epi32(1);
            const __m256i vzd  = _mm256_set1_epi32(zdrop);
            const __m256i vedel = _mm256_set1_epi32(this->e_del);
            const __m256i veins = _mm256_set1_epi32(this->e_ins);
            for (int g = 0; g < SIMD_WIDTH8 / 8; g++) {
                const int base = g * 8;
                __m256i Bg   = _mm256_loadu_si256((const __m256i *)(B + base));
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

                __m256i ms_abs = _mm256_add_epi32(msg, Bg);
                __m256i bag = _mm256_loadu_si256((const __m256i *)(best_abs + base));
                bag = _mm256_max_epi32(bag, ms_abs);
                _mm256_storeu_si256((__m256i *)(best_abs + base), bag);

                __m256i gs_abs = _mm256_add_epi32(hqeg, Bg);
                __m256i gba = _mm256_loadu_si256((const __m256i *)(gbest_abs + base));
                __m256i ge  = _mm256_xor_si256(_mm256_cmpgt_epi32(gba, gs_abs), ff256); // gs_abs >= gba
                __m256i gmask = _mm256_and_si256(qfg, ge);
                gba = _mm256_blendv_epi8(gba, gs_abs, gmask);
                _mm256_storeu_si256((__m256i *)(gbest_abs + base), gba);
                __m256i ierg = _mm256_loadu_si256((const __m256i *)(ierow + base));
                ierg = _mm256_blendv_epi8(ierg, vip1, gmask);
                _mm256_storeu_si256((__m256i *)(ierow + base), ierg);

                __m256i y1c  = _mm256_add_epi32(y1g, vi);
                __m256i yc   = _mm256_add_epi32(yg, _mm256_sub_epi32(xrg, vone));
                __m256i tmpi = _mm256_sub_epi32(vip1, xrg);
                __m256i tmpj = _mm256_sub_epi32(y1c, yc);
                // z-drop gap term weighted by gap-extend penalty, matching the
                // scalar reference (drift>0 -> deletion side *e_del, else *e_ins).
                __m256i zdelta = _mm256_sub_epi32(tmpi, tmpj);
                __m256i zesel  = _mm256_blendv_epi8(veins, vedel,
                                     _mm256_cmpgt_epi32(zdelta, _mm256_setzero_si256()));
                __m256i dif  = _mm256_mullo_epi32(_mm256_abs_epi32(zdelta), zesel);
                __m256i drop = _mm256_sub_epi32(msg, rsg);
                __m256i die  = _mm256_cmpgt_epi32(_mm256_sub_epi32(drop, dif), vzd);
                die = _mm256_and_si256(die, exitg);
                // narrow 8 int32 die masks -> 8 bytes (128-bit packs, no lane cross)
                __m128i d16 = _mm_packs_epi32(_mm256_castsi256_si128(die),
                                              _mm256_extracti128_si256(die, 1));
                __m128i d8  = _mm_packs_epi16(d16, d16);
                __m128i ex  = _mm_loadl_epi64((const __m128i *)(exit_a + base));
                _mm_storel_epi64((__m128i *)(exit_a + base), _mm_andnot_si128(d8, ex));
            }
            exit0 = _mm256_load_si256((__m256i *) exit_a);
        }

#if RDT
        prof[DP1][0] += __rdtsc() - tim1;
        tim1 = __rdtsc();
#endif

        /* Narrowing of the band */
        /* From beg */
        int l;
        for (l = beg; l < end; l++)
        {
            __m256i f256 = _mm256_load_si256((__m256i *)(F + l * SIMD_WIDTH8));
            __m256i h256 = _mm256_load_si256((__m256i *)(H_h + l * SIMD_WIDTH8));
            __m256i tmp = _mm256_or_si256(f256, h256);
            tmp = _mm256_cmpeq_epi8(tmp, zero256);
            uint32_t val = _mm256_movemask_epi8(tmp);
            if (val == 0xFFFFFFFF) nbeg = l;
            else
                break;
        }

        /* From end */
        for (l = end; l >= beg; l--)
        {
            __m256i f256 = _mm256_load_si256((__m256i *)(F + l * SIMD_WIDTH8));
            __m256i h256 = _mm256_load_si256((__m256i *)(H_h + l * SIMD_WIDTH8));
            __m256i tmp = _mm256_or_si256(f256, h256);
            tmp = _mm256_cmpeq_epi8(tmp, zero256);
            uint32_t val = _mm256_movemask_epi8(tmp);
            if (val != 0xFFFFFFFF)
                break;
        }
        // int pnend =nend;
        nend = l + 2 < ncol? l + 2: ncol;
        __m256i tmpb = ff256;

        __m256i exit1 = _mm256_xor_si256(exit0, ff256);
        __m256i l256 = _mm256_set1_epi8(beg - i);   // diagonal offset of column beg

        for (l = beg; l < end; l++)
        {
            __m256i f256 = _mm256_load_si256((__m256i *)(F + l * SIMD_WIDTH8));
            __m256i h256 = _mm256_load_si256((__m256i *)(H_h + l * SIMD_WIDTH8));

            __m256i tmp = _mm256_or_si256(f256, h256);
            tmp = _mm256_or_si256(tmp, exit1);
            tmp = _mm256_cmpeq_epi8(tmp, zero256);
            uint32_t val = _mm256_movemask_epi8(tmp);
            if (val == 0x00) {
                break;
            }
            tmp = _mm256_and_si256(tmp,tmpb);
            l256 = _mm256_add_epi8(l256, one256);
            // NEW
            head256 = _mm256_blendv_epi8(head256, l256, tmp);
            tmpb = tmp;
        }

        __m256i  index256 = tail256;
        tmpb = ff256;

        l256 = _mm256_set1_epi8(end - i);   // diagonal offset of column end
        for (l = end; l >= beg; l--)
        {
            __m256i f256 = _mm256_load_si256((__m256i *)(F + l * SIMD_WIDTH8));
            __m256i h256 = _mm256_load_si256((__m256i *)(H_h + l * SIMD_WIDTH8));

            __m256i tmp = _mm256_or_si256(f256, h256);
            tmp = _mm256_or_si256(tmp, exit1);
            tmp = _mm256_cmpeq_epi8(tmp, zero256);
            uint32_t val = _mm256_movemask_epi8(tmp);
            if (val == 0x00)  {
                break;
            }

            tmp = _mm256_and_si256(tmp,tmpb);
            l256 = _mm256_sub_epi8(l256, one256);
            // NEW
            index256 = _mm256_blendv_epi8(index256, l256, tmp);
            tmpb = tmp;
        }
        index256 = _mm256_add_epi8(index256, two256);
        // signed min in the offset frame against qlen-i
        tail256 = _mm256_min_epi8(index256, qlen_off256);

        // Frame shift for the next row: i advances by 1, so the same absolute
        // band edge has its diagonal offset (col - i) decremented by 1. Keep
        // head/tail tracking the same columns as the frame moves.
        head256 = _mm256_sub_epi8(head256, one256);
        tail256 = _mm256_sub_epi8(tail256, one256);

        // --- RE-BASELINE EVENT (per-lane score floor) ---
        // After this row's score state is final, lower any lane whose row-max
        // (maxRS1, rebaselined) has climbed to REBASE_HI back into range. We
        // raise B[l] by delta = maxRS1[l] - REBASE_KEEP and subtract that delta
        // (saturating to 0) from every carried score buffer/register for that
        // lane. Done AFTER band narrowing so this row's trim matches the 16-bit
        // reference; the next row reads the lowered values consistently.
        {
            /* Fast path: one vector test for "any lane >= REBASE_HI?" — see the
             * 128-bit kernel. Under the routing envelope this is ~always false, so
             * the per-lane scalar scan + maxRS1 store below are skipped each row. */
            __m256i hi256 = _mm256_set1_epi8((int8_t) REBASE_HI);
            if (_mm256_movemask_epi8(_mm256_cmpeq_epi8(_mm256_subs_epu8(hi256, maxRS1), zero256))) {
            uint8_t rs_b[SIMD_WIDTH8] __attribute((aligned(32)));
            _mm256_store_si256((__m256i *) rs_b, maxRS1);
            int8_t  delta_a[SIMD_WIDTH8] __attribute((aligned(32)));
            int any = 0;
            for (int l = 0; l < SIMD_WIDTH8; l++) {
                int d = 0;
                if ((int)rs_b[l] >= REBASE_HI) {
                    d = (int)rs_b[l] - REBASE_KEEP;   // keep low REBASE_KEEP, in range
                    B[l] += d;
                    any = 1;
                }
                delta_a[l] = (int8_t) d;              // 0..(REBASE_HI-REBASE_KEEP) small
            }
            if (any) {
                __m256i delta256 = _mm256_load_si256((__m256i *) delta_a);
                // Subtract delta from every carried score buffer over the active
                // band columns. The range [lo,hi] is the union of the current and
                // next band windows; clamping to [0, ncol-1] keeps the loop within
                // the valid column range (a safe superset covering every column the
                // next row will read).
                int lo = nbeg < beg ? nbeg : beg;
                int hi = nend > end ? nend : end;
                if (lo < 0) lo = 0;
                if (hi + 1 > ncol) hi = ncol - 1;
                for (int l = lo; l <= hi; l++) {
                    __m256i h256 = _mm256_load_si256((__m256i *)(H_h + l * SIMD_WIDTH8));
                    __m256i f256 = _mm256_load_si256((__m256i *)(F   + l * SIMD_WIDTH8));
                    h256 = _mm256_subs_epu8(h256, delta256);
                    f256 = _mm256_subs_epu8(f256, delta256);
                    _mm256_store_si256((__m256i *)(H_h + l * SIMD_WIDTH8), h256);
                    _mm256_store_si256((__m256i *)(F   + l * SIMD_WIDTH8), f256);
                }
                // The vertical seed H_v is read only while the band still touches
                // column 0 (beg==0). Lower the rows that may still be consumed.
                if (beg == 0) {
                    for (int r = i + 1; r < nrow; r++) {
                        __m256i v256 = _mm256_load_si256((__m256i *)(H_v + r * SIMD_WIDTH8));
                        v256 = _mm256_subs_epu8(v256, delta256);
                        _mm256_store_si256((__m256i *)(H_v + r * SIMD_WIDTH8), v256);
                    }
                }
                // Carried score maxima (registers) live in the same rebaselined
                // frame and must be lowered too so the next row's z-drop (ms-rs)
                // and the absolute-frame add-back stay consistent. maxScore256 is
                // the running max, so it never saturates to 0 (maxScore256 >=
                // maxRS1 for the triggering row, and maxRS1 >= REBASE_HI >
                // REBASE_KEEP > delta, so it stays positive).
                maxScore256 = _mm256_subs_epu8(maxScore256, delta256);
                // gscore no longer carries a rebaselined byte register: the
                // query-end cell H is captured per row and accumulated wide
                // (gbest_abs, absolute units) in the epilogue, immune to this
                // re-baseline subtract.
            }
            }   /* end "any lane >= REBASE_HI" fast-path guard */
        }

#if RDT
        prof[DP2][0] += __rdtsc() - tim1;
#endif
    }

#if RDT
    prof[DP][0] += __rdtsc() - tim;
#endif

    // Scores are reconstructed in ABSOLUTE units from the per-lane wide
    // best_abs/gbest_abs side channels (= rebaselined byte + B[l], captured each
    // row before any later re-baseline could saturate the byte away). Positions
    // are reconstructed wide from the diagonal-offset lanes plus the per-lane
    // best-row side channels.
    int8_t maxj[SIMD_WIDTH8]  __attribute((aligned(64)));
    _mm256_store_si256((__m256i *) maxj, y256);   // best col as diagonal offset

    int8_t max_off_ar[SIMD_WIDTH8]  __attribute((aligned(64)));
    _mm256_store_si256((__m256i *) max_off_ar, max_off256);

    for(i = 0; i < SIMD_WIDTH8; i++)
    {
        p[i].score   = best_abs[i];                               // absolute score
        p[i].tle     = xrow[i];                                   // best row (target end)
        // qle reconstruction: maxj[l] = col - (xrow-1), so col = maxj[l] + (xrow-1).
        // Guard the unset case (xrow==0 means no update occurred).
        p[i].qle     = (xrow[i] == 0) ? 0 : ((int) maxj[i] + xrow[i] - 1);
        p[i].max_off = (uint8_t) max_off_ar[i];
        p[i].gscore  = gbest_abs[i];                              // absolute gscore (-1 if unset)
        p[i].gtle    = ierow[i];                                  // best gscore row
    }

    return;
}

// ------------------------- AVX2 - 16 bit SIMD_LANES ---------------------------
#define PFD 2
void BandedPairWiseSW::getScores16(SeqPair *pairArray,
                                   uint8_t *seqBufRef,
                                   uint8_t *seqBufQer,
                                   int32_t numPairs,
                                   uint16_t numThreads,
                                   int32_t w)
{
    int64_t startTick, endTick;

    {
        BswOvershootGuard _g(pairArray, numPairs, SIMD_WIDTH16, guard_overshoot_);
        smithWatermanBatchWrapper16(pairArray, seqBufRef, seqBufQer, numPairs, numThreads, w);
    }

#if MAXI
    printf("AVX2 Vecor code: Writing output..\n");
    for (int l=0; l<numPairs; l++)
    {
        fprintf(stderr, "%d (%d %d) %d %d %d\n",
                pairArray[l].score, pairArray[l].x, pairArray[l].y,
                pairArray[l].gscore, pairArray[l].max_off, pairArray[l].max_ie);

    }
    printf("Vector code: Writing output completed!!!\n\n");
#endif
    
}

void BandedPairWiseSW::smithWatermanBatchWrapper16(SeqPair *pairArray,
                                                   uint8_t *seqBufRef,
                                                   uint8_t *seqBufQer,
                                                   int32_t numPairs,
                                                   uint16_t numThreads,
                                                   int32_t w)
{
    numThreads = effective_threads(numThreads);
    int64_t st1, st2, st3, st4, st5;
#if RDT
    st1 = ___rdtsc();
#endif

    uint16_t *seq1SoA = (uint16_t *)_mm_malloc((size_t)MAX_SEQ_LEN16 * SIMD_WIDTH16 * numThreads * sizeof(uint16_t), 64);
    uint16_t *seq2SoA = (uint16_t *)_mm_malloc((size_t)MAX_SEQ_LEN16 * SIMD_WIDTH16 * numThreads * sizeof(uint16_t), 64);

    if (UNLIKELY(seq1SoA == NULL || seq2SoA == NULL)) {
        fprintf(stderr, "Error! Mem not allocated!!!\n");
        exit(EXIT_FAILURE);
    }
    
    int32_t ii;
    int32_t roundNumPairs = ((numPairs + SIMD_WIDTH16 - 1)/SIMD_WIDTH16 ) * SIMD_WIDTH16;
    for(ii = numPairs; ii < roundNumPairs; ii++)
    {
        pairArray[ii].id = ii;
        pairArray[ii].len1 = 0;
        pairArray[ii].len2 = 0;
        pairArray[ii].idr = 0;
        pairArray[ii].idq = 0;
    }

#if RDT 
    st2 = ___rdtsc();
#endif
    

#if RDT 
    st3 = ___rdtsc();
#endif

    int eb = end_bonus;
//#pragma omp parallel num_threads(numThreads)
    {
        int32_t i;
        uint16_t tid = 0; 
        uint16_t *mySeq1SoA = seq1SoA + tid * MAX_SEQ_LEN16 * SIMD_WIDTH16;
        uint16_t *mySeq2SoA = seq2SoA + tid * MAX_SEQ_LEN16 * SIMD_WIDTH16;
        uint8_t *seq1;
        uint8_t *seq2;
        uint16_t h0[SIMD_WIDTH16]   __attribute__((aligned(64)));
        uint16_t band[SIMD_WIDTH16];        
        uint16_t qlen[SIMD_WIDTH16] __attribute__((aligned(64)));
        int32_t bsize = 0;
        
        int16_t *H1 = H16_ + tid * SIMD_WIDTH16 * MAX_SEQ_LEN16;
        int16_t *H2 = H16__ + tid * SIMD_WIDTH16 * MAX_SEQ_LEN16;
        
        __m256i zero256   = _mm256_setzero_si256();
        __m256i e_ins256  = _mm256_set1_epi16(e_ins);
        __m256i oe_ins256 = _mm256_set1_epi16(o_ins + e_ins);
        __m256i o_del256  = _mm256_set1_epi16(o_del);
        __m256i e_del256  = _mm256_set1_epi16(e_del);
        __m256i eb_ins256 = _mm256_set1_epi16(eb - o_ins);
        __m256i eb_del256 = _mm256_set1_epi16(eb - o_del);
        
        int16_t max = 0;
        if (max < w_match) max = w_match;
        if (max < w_mismatch) max = w_mismatch;
        if (max < w_ambig) max = w_ambig;
        
        int nstart = 0, nend = numPairs;

//#pragma omp for schedule(dynamic, 128)
        for(i = nstart; i < nend; i+=SIMD_WIDTH16)
        {
            int32_t j, k;
            uint16_t maxLen1 = 0;
            uint16_t maxLen2 = 0;
            bsize = w;

            uint64_t tim;
            for(j = 0; j < SIMD_WIDTH16; j++)
            {
                if ((i + j + PFD) < roundNumPairs) { // prefetch block (bounded; see getScores8/16 contract)
                    SeqPair spf = pairArray[i + j + PFD];
                    _mm_prefetch((const char*) seqBufRef + (int64_t)spf.idr, _MM_HINT_NTA);
                    _mm_prefetch((const char*) seqBufRef + (int64_t)spf.idr + 64, _MM_HINT_NTA);
                }

                SeqPair sp = pairArray[i + j];
                h0[j] = sp.h0;
                seq1 = seqBufRef + (int64_t)sp.idr;
                
                for(k = 0; k < sp.len1; k++)
                {
                    mySeq1SoA[k * SIMD_WIDTH16 + j] = (seq1[k] == AMBIG?0xFFFF:seq1[k]);
                }
                qlen[j] = sp.len2 * max;
                if(maxLen1 < sp.len1) maxLen1 = sp.len1;
            }

            for(j = 0; j < SIMD_WIDTH16; j++)
            {
                SeqPair sp = pairArray[i + j];
                for(k = sp.len1; k <= maxLen1; k++) //removed "="
                {
                    mySeq1SoA[k * SIMD_WIDTH16 + j] = DUMMY1;
                }
            }
            /* B5: only the boundary row H2[maxLen1] survives the h0-prefix
             * deletion seed below; write just that row, before the seed. */
            _mm256_store_si256((__m256i *)(H2 + maxLen1 * SIMD_WIDTH16), _mm256_set1_epi16((short)DUMMY1));
//--------------------
            __m256i h0_256 = _mm256_load_si256((__m256i*) h0);
            _mm256_store_si256((__m256i *) H2, h0_256);
            __m256i tmp256 = _mm256_sub_epi16(h0_256, o_del256);

            for(k = 1; k < maxLen1; k++) {
                tmp256 = _mm256_sub_epi16(tmp256, e_del256);
                __m256i tmp256_ = _mm256_max_epi16(tmp256, zero256);
                _mm256_store_si256((__m256i *)(H2 + k* SIMD_WIDTH16), tmp256_);
            }
//-------------------
            for(j = 0; j < SIMD_WIDTH16; j++)
            {
                if ((i + j + PFD) < roundNumPairs) { // prefetch block (bounded; see getScores8/16 contract)
                    SeqPair spf = pairArray[i + j + PFD];
                    _mm_prefetch((const char*) seqBufQer + (int64_t)spf.idq, _MM_HINT_NTA);
                    _mm_prefetch((const char*) seqBufQer + (int64_t)spf.idq + 64, _MM_HINT_NTA);
                }
                
                SeqPair sp = pairArray[i + j];
                //seq2 = seqBufQer + (int64_t)sp.id * MAX_SEQ_LEN_QER;
                seq2 = seqBufQer + (int64_t)sp.idq;             
                for(k = 0; k < sp.len2; k++)
                {
                    mySeq2SoA[k * SIMD_WIDTH16 + j] = (seq2[k]==AMBIG?0xFFFF:seq2[k]);
                }
                if(maxLen2 < sp.len2) maxLen2 = sp.len2;
            }
            
            for(j = 0; j < SIMD_WIDTH16; j++)
            {
                SeqPair sp = pairArray[i + j];
                for(k = sp.len2; k <= maxLen2; k++)
                {
                    mySeq2SoA[k * SIMD_WIDTH16 + j] = DUMMY2;
                }
            }
            /* B5: only boundary row H1[maxLen2]=0 survives the seed below. */
            _mm256_store_si256((__m256i *)(H1 + maxLen2 * SIMD_WIDTH16), _mm256_setzero_si256());
//------------------------
            _mm256_store_si256((__m256i *) H1, h0_256);
            __m256i cmp256 = _mm256_cmpgt_epi16(h0_256, oe_ins256);
            tmp256 = _mm256_sub_epi16(h0_256, oe_ins256);

            tmp256 = _mm256_blendv_epi16(zero256, tmp256, cmp256);
            _mm256_store_si256((__m256i *) (H1 + SIMD_WIDTH16), tmp256);
            for(k = 2; k < maxLen2; k++)
            {
                __m256i h1_256 = tmp256;
                tmp256 = _mm256_sub_epi16(h1_256, e_ins256);
                tmp256 = _mm256_max_epi16(tmp256, zero256);
                _mm256_store_si256((__m256i *)(H1 + k*SIMD_WIDTH16), tmp256);
            }
//------------------------
            uint16_t myband[SIMD_WIDTH16] __attribute__((aligned(64)));
            uint16_t temp[SIMD_WIDTH16] __attribute__((aligned(64)));
            {
                __m256i qlen256 = _mm256_load_si256((__m256i *) qlen);
                __m256i sum256 = _mm256_add_epi16(qlen256, eb_ins256);
                _mm256_store_si256((__m256i *) temp, sum256);               
                for (int l=0; l<SIMD_WIDTH16; l++) {
                    double val = temp[l]/e_ins + 1.0;
                    int max_ins = val;
                    max_ins = max_ins > 1? max_ins : 1;
                    myband[l] = min_(bsize, max_ins);
                }
                sum256 = _mm256_add_epi16(qlen256, eb_del256);
                _mm256_store_si256((__m256i *) temp, sum256);               
                for (int l=0; l<SIMD_WIDTH16; l++) {
                    double val = temp[l]/e_del + 1.0;
                    int max_ins = val;
                    max_ins = max_ins > 1? max_ins : 1;
                    myband[l] = min_(myband[l], max_ins);
                    bsize = bsize < myband[l] ? myband[l] : bsize;                  
                }
            }

            smithWaterman256_16(mySeq1SoA,
                                mySeq2SoA,
                                maxLen1,
                                maxLen2,
                                pairArray + i,
                                h0,
                                tid,
                                numPairs,
                                zdrop,
                                bsize, 
                                qlen,
                                myband);
        }
    }

#if RDT
    st4 = ___rdtsc();
#endif
    

#if RDT
    st5 = ___rdtsc();
    setupTicks += st2 - st1;
    sort1Ticks += st3 - st2;
    swTicks += st4 - st3;
    sort2Ticks += st5 - st4;
#endif

    // free mem
    _mm_free(seq1SoA);
    _mm_free(seq2SoA);

    return;
}

void BandedPairWiseSW::smithWaterman256_16(uint16_t seq1SoA[],
                                           uint16_t seq2SoA[],
                                           uint16_t nrow,
                                           uint16_t ncol,
                                           SeqPair *p,
                                           uint16_t h0[],
                                           uint16_t tid,
                                           int32_t numPairs,
                                           int zdrop,
                                           int32_t w,
                                           uint16_t qlen[],
                                           uint16_t myband[])
{
    __m256i match256     = _mm256_set1_epi16(this->w_match);
    __m256i mismatch256  = _mm256_set1_epi16(this->w_mismatch);
    __m256i w_ambig_256  = _mm256_set1_epi16(this->w_ambig);    // ambig penalty

    // D3 generic-matrix seam: symmetric default uses SBT_PREPASS16_SYM; a single
    // freed-to-match cell (bisulfite) uses the rank-1 path; any other asymmetric
    // matrix uses the amat LUT. gen_mat is false on the hot path.
    const bool forced  = bsw_force_generic_matrix();
    const bool gen_mat = bsw_generic_matrix(this->mat, this->w_match, this->w_mismatch)
                         || forced;
    const BswFreedCell fc = bsw_freed_cell(this->mat, this->w_match, this->w_mismatch, forced);
    __m256i frref256  = _mm256_set1_epi16(fc.ref);
    __m256i frread256 = _mm256_set1_epi16(fc.read);
    int8_t amat_bytes[16] __attribute__((aligned(16)));
    build_amat16(amat_bytes, this->mat);
    __m256i amat256   = _mm256_broadcastsi128_si256(_mm_load_si128((__m128i *)amat_bytes));
    __m256i three256  = _mm256_set1_epi16(3);

    __m256i e_del256    = _mm256_set1_epi16(this->e_del);
    __m256i oe_del256   = _mm256_set1_epi16(this->o_del + this->e_del);
    __m256i e_ins256    = _mm256_set1_epi16(this->e_ins);
    __m256i oe_ins256   = _mm256_set1_epi16(this->o_ins + this->e_ins);

    int16_t *F  = F16_ + tid * SIMD_WIDTH16 * MAX_SEQ_LEN16;
    int16_t *H_h    = H16_ + tid * SIMD_WIDTH16 * MAX_SEQ_LEN16;
    int16_t *H_v = H16__ + tid * SIMD_WIDTH16 * MAX_SEQ_LEN16;

    int lane = 0;

    int16_t i, j;

    uint16_t tlen[SIMD_WIDTH16];
    uint16_t tail[SIMD_WIDTH16] __attribute((aligned(64)));
    uint16_t head[SIMD_WIDTH16] __attribute((aligned(64)));

    // PR 17: per-row score-vector scratch for fission (AVX2 16-bit).
    int16_t *sbt_buf = sbt16_ + tid * SIMD_WIDTH16 * MAX_SEQ_LEN16;
    
    int32_t minq = 10000000;
    for (int l=0; l<SIMD_WIDTH16; l++) {
        tlen[l] = p[l].len1;
        qlen[l] = p[l].len2;
        if (p[l].len2 < minq) minq = p[l].len2;
    }
    minq -= 1; // for gscore

    __m256i tlen256 = _mm256_load_si256((__m256i *) tlen);
    __m256i qlen256 = _mm256_load_si256((__m256i *) qlen);
    __m256i myband256 = _mm256_load_si256((__m256i *) myband);
    __m256i zero256 = _mm256_setzero_si256();
    __m256i one256  = _mm256_set1_epi16(1);
    __m256i two256  = _mm256_set1_epi16(2);
    __m256i max_ie256 = zero256;
    __m256i ff256 = _mm256_set1_epi16(0xFFFF);
        
    __m256i tail256 = qlen256, head256 = zero256;
    _mm256_store_si256((__m256i *) head, head256);
    _mm256_store_si256((__m256i *) tail, tail256);

    __m256i mlen256 = _mm256_add_epi16(qlen256, myband256);
    mlen256 = _mm256_min_epu16(mlen256, tlen256);

    uint16_t temp[SIMD_WIDTH16]  __attribute((aligned(64)));
    uint16_t temp1[SIMD_WIDTH16]  __attribute((aligned(64)));
    
    __m256i s00  = _mm256_load_si256((__m256i *)(seq1SoA));
    __m256i hval = _mm256_load_si256((__m256i *)(H_v));
    __mmask16 dmask = 0xFFFF;
    __mmask32 dmask32 = 0xAAAAAAAA;
        
    __m256i maxScore256 = hval;
    for(j = 0; j < ncol; j++)
        _mm256_store_si256((__m256i *)(F + j * SIMD_WIDTH16), zero256);
    
    __m256i x256 = zero256;
    __m256i y256 = zero256;
    __m256i i256 = zero256;
    __m256i gscore = _mm256_set1_epi16(-1);
    __m256i max_off256 = zero256;
    __m256i exit0 = _mm256_set1_epi16(0xFFFF);
    __m256i zdrop256 = _mm256_set1_epi16(zdrop);
    
    int beg = 0, end = ncol;
    int nbeg = beg, nend = end;

#if RDT
    uint64_t tim = __rdtsc();
#endif
    
    for(i = 0; i < nrow; i++)
    {       
        __m256i e11 = zero256;
        __m256i h00, h11, h10;
        __m256i s10 = _mm256_load_si256((__m256i *)(seq1SoA + (i + 0) * SIMD_WIDTH16));

        beg = nbeg; end = nend;
        int pbeg = beg;
        if (beg < i - w) beg = i - w;
        if (end > i + w + 1) end = i + w + 1;
        if (end > ncol) end = ncol;

        h10 = zero256;
        if (beg == 0)
            h10 = _mm256_load_si256((__m256i *)(H_v + (i+1) * SIMD_WIDTH16));

        __m256i j256 = zero256;
        __m256i maxRS1;
        maxRS1 = zero256;

        __m256i i1_256 = _mm256_set1_epi16(i+1);
        __m256i y1_256 = zero256;
        
#if RDT 
        uint64_t tim1 = __rdtsc();
#endif
        
        __m256i i256, cache256;
        __m256i phead256 = head256, ptail256 = tail256;
        i256 = _mm256_set1_epi16(i);
        cache256 = _mm256_sub_epi16(i256, myband256);
        head256 = _mm256_max_epi16(head256, cache256);
        cache256 = _mm256_add_epi16(i1_256, myband256);
        tail256 = _mm256_min_epu16(tail256, cache256);
        tail256 = _mm256_min_epu16(tail256, qlen256);

        // NEW, trimming.
        __m256i cmph = _mm256_cmpeq_epi16(head256, phead256);
        __m256i cmpt = _mm256_cmpeq_epi16(tail256, ptail256);
        // cmph &= cmpt;
        cmph = _mm256_and_si256(cmph, cmpt);
        //__mmask16 cmp_ht = _mm256_movepi16_mask(cmph);
        __mmask32 cmp_ht = _mm256_movemask_epi8(cmph) & dmask32;
        
        for (int l=beg; l<end && cmp_ht != dmask32; l++)
        {
            __m256i h256 = _mm256_load_si256((__m256i *)(H_h + l * SIMD_WIDTH16));
            __m256i f256 = _mm256_load_si256((__m256i *)(F + l * SIMD_WIDTH16));
            
            __m256i pj256 = _mm256_set1_epi16(l);
            __m256i j256 = _mm256_set1_epi16(l+1);
            __m256i cmp1 = _mm256_cmpgt_epi16(head256, pj256);
            //uint16_t cval = _mm256_movepi16_mask(cmp1);
            uint32_t cval = _mm256_movemask_epi8(cmp1) & dmask32;
            if (cval == 0x00) break;
            //__m256i cmp2 = _mm256_cmpgt_epi16(pj256, tail256);
            __m256i cmp2 = _mm256_cmpgt_epi16(j256, tail256);
            cmp1 = _mm256_or_si256(cmp1, cmp2);
            h256 = _mm256_andnot_si256(cmp1, h256);
            f256 = _mm256_andnot_si256(cmp1, f256);
            
            _mm256_store_si256((__m256i *)(F + l * SIMD_WIDTH16), f256);
            _mm256_store_si256((__m256i *)(H_h + l * SIMD_WIDTH16), h256);
        }

#if RDT
        prof[DP3][0] += __rdtsc() - tim1;
#endif

        // beg = nbeg; end = nend;
        __m256i cmp256_1 = _mm256_cmpgt_epi16(i1_256, tlen256);
        
        __m256i cmpim = _mm256_cmpgt_epi16(i1_256, mlen256);
        __m256i cmpht = _mm256_cmpeq_epi16(tail256, head256);
        cmpim = _mm256_or_si256(cmpim, cmpht);

        // NEW
        cmpht = _mm256_cmpgt_epi16(head256, tail256);
        cmpim = _mm256_or_si256(cmpim, cmpht);

        exit0 = _mm256_andnot_si256(cmpim, exit0);

        
#if RDT
        tim1 = __rdtsc();
#endif
        
        // PR 17: AVX2 16-bit score pre-pass (fission). gen_mat branch is
        // loop-invariant: symmetric keeps the cheap SYM prepass.
        if (!gen_mat) {
            for (int jp = beg; jp < end; jp++) {
                __m256i s2 = _mm256_load_si256((__m256i *)(seq2SoA + jp * SIMD_WIDTH16));
                __m256i sbt11;
                SBT_PREPASS16_SYM(s10, s2, sbt11, mismatch256, match256, w_ambig_256);
                _mm256_store_si256((__m256i *)(sbt_buf + jp * SIMD_WIDTH16), sbt11);
            }
        } else if (fc.rank1) {
            __m256i alt10 = _mm256_blendv_epi16(s10, frread256, _mm256_cmpeq_epi16(s10, frref256));
            for (int jp = beg; jp < end; jp++) {
                __m256i s2 = _mm256_load_si256((__m256i *)(seq2SoA + jp * SIMD_WIDTH16));
                __m256i sbt11;
                SBT_PREPASS16_RANK1(s10, s2, alt10, sbt11, mismatch256, match256, w_ambig_256);
                _mm256_store_si256((__m256i *)(sbt_buf + jp * SIMD_WIDTH16), sbt11);
            }
        } else {
            for (int jp = beg; jp < end; jp++) {
                __m256i s2 = _mm256_load_si256((__m256i *)(seq2SoA + jp * SIMD_WIDTH16));
                __m256i sbt11;
                SBT_PREPASS16_AMAT(s10, s2, sbt11, amat256, mismatch256, w_ambig_256, three256);
                _mm256_store_si256((__m256i *)(sbt_buf + jp * SIMD_WIDTH16), sbt11);
            }
        }

        j256 = _mm256_set1_epi16(beg);
        for(j = beg; j < end; j++)
        {
            __m256i f11, f21, sbt11;
            sbt11 = _mm256_load_si256((__m256i *)(sbt_buf + j * SIMD_WIDTH16));
            h00 = _mm256_load_si256((__m256i *)(H_h + j * SIMD_WIDTH16));
            f11 = _mm256_load_si256((__m256i *)(F + j * SIMD_WIDTH16));

            __m256i pj256 = j256;
            j256 = _mm256_add_epi16(j256, one256);

            MAIN_CODE16_CORE(sbt11, h00, h11, e11, f11, f21, zero256,
                             e_ins256, oe_ins256,
                             e_del256, oe_del256);
            
            // Masked writing
            __m256i cmp2 = _mm256_cmpgt_epi16(head256, pj256);
            __m256i cmp1 = _mm256_cmpgt_epi16(pj256, tail256);
            cmp1 = _mm256_or_si256(cmp1, cmp2);
            h10 = _mm256_andnot_si256(cmp1, h10);
            f21 = _mm256_andnot_si256(cmp1, f21);
            
            __m256i bmaxRS = maxRS1;
            maxRS1 =_mm256_max_epi16(maxRS1, h11);
            // "new row-max" argmax mask. maxRS1 = max(bmaxRS,h11), so
            // cmpgt(maxRS1,bmaxRS) is a strict subset of cmpeq(maxRS1,h11)
            // (both mean h11 >= bmaxRS); the OR was redundant. cmpeq(maxRS1,h11)
            // is the exact combined mask — bit-identical, drops a cmpgt+or.
            __m256i cmpA = _mm256_cmpeq_epi16(maxRS1, h11);
            cmp1 = _mm256_cmpgt_epi16(j256, tail256); // change
            cmp1 = _mm256_or_si256(cmp1, cmp2);         // change
            cmpA = _mm256_blendv_epi16(y1_256, j256, cmpA);
            y1_256 = _mm256_blendv_epi16(cmpA, y1_256, cmp1);
            maxRS1 = _mm256_blendv_epi16(maxRS1, bmaxRS, cmp1);                     

            _mm256_store_si256((__m256i *)(F + j * SIMD_WIDTH16), f21);
            _mm256_store_si256((__m256i *)(H_h + j * SIMD_WIDTH16), h10);

            h10 = h11;
            
            //j256 = _mm256_add_epi16(j256, one256);
            
            // gscore calculations
            if (j >= minq)
            {
                __m256i cmp = _mm256_cmpeq_epi16(j256, qlen256);
                // Both blendv pairs below fall back to the same operand under
                // exit0 then cmp, so each collapses to one select against the
                // fused mask (cmp & exit0): saves 2 port-5 blendv per qualifying
                // cell. Bit-identical — lane masks are uniform 0xFFFF/0x0000.
                __m256i sel = _mm256_and_si256(cmp, exit0);
                __m256i max_gh = _mm256_max_epi16(gscore, h11);
                __m256i cmp_gh = _mm256_cmpgt_epi16(gscore, h11);
                __m256i cand_ie = _mm256_blendv_epi16(i1_256, max_ie256, cmp_gh);

                __m256i tmp256_1 = _mm256_blendv_epi16(max_ie256, cand_ie, sel);

                max_gh = _mm256_blendv_epi16(gscore, max_gh, sel);

                cmp = _mm256_cmpgt_epi16(j256, tail256);
                max_gh = _mm256_blendv_epi16(max_gh, gscore, cmp);
                max_ie256 = _mm256_blendv_epi16(tmp256_1, max_ie256, cmp);
                gscore = max_gh;            
            }
        }
        __m256i cmp2 = _mm256_cmpgt_epi16(head256, j256);
        __m256i cmp1 = _mm256_cmpgt_epi16(j256, tail256);
        cmp1 = _mm256_or_si256(cmp1, cmp2);
        h10 = _mm256_andnot_si256(cmp1, h10);
        
        _mm256_store_si256((__m256i *)(H_h + j * SIMD_WIDTH16), h10);
        _mm256_store_si256((__m256i *)(F + j * SIMD_WIDTH16), zero256);
                
        /* exit due to zero score by a row */
        __m256i bmaxScore256 = maxScore256;
        __m256i tmp = _mm256_cmpeq_epi16(maxRS1, zero256);
        uint32_t cval = _mm256_movemask_epi8(tmp) & dmask32;
        if (cval == dmask32) break;

        exit0 = _mm256_andnot_si256(tmp, exit0);

        __m256i score256 = _mm256_max_epi16(maxScore256, maxRS1);
        maxScore256 = _mm256_blendv_epi16(maxScore256, score256, exit0);

        __m256i cmp = _mm256_cmpgt_epi16(maxScore256, bmaxScore256);
        y256 = _mm256_blendv_epi16(y256, y1_256, cmp);
        x256 = _mm256_blendv_epi16(x256, i1_256, cmp);      
        // max_off calculations
        tmp = _mm256_sub_epi16(y1_256, i1_256);
        tmp = _mm256_abs_epi16(tmp);
        __m256i bmax_off256 = max_off256;
        tmp = _mm256_max_epi16(max_off256, tmp);
        max_off256 = _mm256_blendv_epi16(bmax_off256, tmp, cmp);

        // Z-score
        ZSCORE16(i1_256, y1_256);       

#if RDT
        prof[DP1][0] += __rdtsc() - tim1;
        tim1 = __rdtsc();
#endif
        
        /* Narrowing of the band */
        /* From beg */
        int l;
        for (l = beg; l < end; l++) {
            __m256i f256 = _mm256_load_si256((__m256i *)(F + l * SIMD_WIDTH16));
            __m256i h256 = _mm256_load_si256((__m256i *)(H_h + l * SIMD_WIDTH16));
            __m256i tmp = _mm256_or_si256(f256, h256);
            tmp = _mm256_cmpeq_epi16(tmp, zero256);
            //uint16_t val = _mm256_movepi16_mask(tmp);
            uint32_t val = _mm256_movemask_epi8(tmp) & dmask32;
            if (val == dmask32) nbeg = l;
            else
                break;
        }
        
        /* From end */
        bool flg = 1;
        for (l = end; l >= beg; l--)
        {
            __m256i f256 = _mm256_load_si256((__m256i *)(F + l * SIMD_WIDTH16));
            __m256i h256 = _mm256_load_si256((__m256i *)(H_h + l * SIMD_WIDTH16));
            __m256i tmp = _mm256_or_si256(f256, h256);
            tmp = _mm256_cmpeq_epi16(tmp, zero256);
            //uint16_t val = _mm256_movepi16_mask(tmp);
            uint32_t val = _mm256_movemask_epi8(tmp) & dmask32;
            if (val != dmask32 && flg)  
                break;
        }
        nend = l + 2 < ncol? l + 2: ncol;

        __m256i tail256_ = _mm256_sub_epi16(tail256, one256);
        __m256i tmpb = ff256;
        __m256i exit1 = _mm256_xor_si256(exit0, ff256);
        __m256i l256 = _mm256_set1_epi16(beg);
        
        for (l = beg; l < end; l++)
        {
            __m256i f256 = _mm256_load_si256((__m256i *)(F + l * SIMD_WIDTH16));
            __m256i h256 = _mm256_load_si256((__m256i *)(H_h + l * SIMD_WIDTH16));
    
            __m256i tmp = _mm256_or_si256(f256, h256);
            tmp = _mm256_or_si256(tmp, exit1);          
            tmp = _mm256_cmpeq_epi16(tmp, zero256);
            //uint16_t val = _mm256_movepi16_mask(tmp);
            uint32_t val = _mm256_movemask_epi8(tmp) & dmask32;
            if (val == 0x00) {
                break;
            }
            tmp = _mm256_and_si256(tmp,tmpb);
            //__m256i l256 = _mm256_set1_epi16(l+1);
            l256 = _mm256_add_epi16(l256, one256);

            head256 = _mm256_blendv_epi16(head256, l256, tmp);

            tmpb = tmp;         
        }
        // _mm256_store_si256((__m256i *) head, head256);
        
        __m256i  index256 = tail256;
        tmpb = ff256;

        l256 = _mm256_set1_epi16(end);
        for (l = end; l >= beg; l--)
        {
            __m256i f256 = _mm256_load_si256((__m256i *)(F + l * SIMD_WIDTH16));
            __m256i h256 = _mm256_load_si256((__m256i *)(H_h + l * SIMD_WIDTH16));
            
            __m256i tmp = _mm256_or_si256(f256, h256);
            tmp = _mm256_or_si256(tmp, exit1);
            tmp = _mm256_cmpeq_epi16(tmp, zero256);         
            //uint16_t val = _mm256_movepi16_mask(tmp);
            uint32_t val = _mm256_movemask_epi8(tmp) & dmask32;
            if (val == 0x00)  {
                break;
            }
            tmp = _mm256_and_si256(tmp,tmpb);
            l256 = _mm256_sub_epi16(l256, one256);

            // NEW
            index256 = _mm256_blendv_epi8(index256, l256, tmp);

            tmpb = tmp;
        }
        index256 = _mm256_add_epi16(index256, two256);
        tail256 = _mm256_min_epi16(index256, qlen256);
        // _mm256_store_si256((__m256i *) tail, tail256);       

#if RDT
        prof[DP2][0] += __rdtsc() - tim1;
#endif
    }
    
#if RDT
    prof[DP][0] += __rdtsc() - tim;
#endif
    
    int16_t score[SIMD_WIDTH16]  __attribute((aligned(64)));
    _mm256_store_si256((__m256i *) score, maxScore256);

    int16_t maxi[SIMD_WIDTH16]  __attribute((aligned(64)));
    _mm256_store_si256((__m256i *) maxi, x256);

    int16_t maxj[SIMD_WIDTH16]  __attribute((aligned(64)));
    _mm256_store_si256((__m256i *) maxj, y256);

    int16_t max_off_ar[SIMD_WIDTH16]  __attribute((aligned(64)));
    _mm256_store_si256((__m256i *) max_off_ar, max_off256);

    int16_t gscore_ar[SIMD_WIDTH16]  __attribute((aligned(64)));
    _mm256_store_si256((__m256i *) gscore_ar, gscore);

    int16_t maxie_ar[SIMD_WIDTH16]  __attribute((aligned(64)));
    _mm256_store_si256((__m256i *) maxie_ar, max_ie256);
    
    for(i = 0; i < SIMD_WIDTH16; i++)
    {
        p[i].score = score[i];
        p[i].tle = maxi[i];
        p[i].qle = maxj[i];
        p[i].max_off = max_off_ar[i];
        p[i].gscore = gscore_ar[i];
        p[i].gtle = maxie_ar[i];
    }
    
    return;
}

#endif // AVX2



#if __AVX512BW__

// ----------------------------------------------------------------------------------
// AVX512- vec8, vec16 SIMD code
//
// ----------------------------------------------------------------------------------

// ------------------------ vec 8 --------------------------------------------------
// ZSCORE8 is unused in smithWaterman512_8 (z-drop replaced by the wide scalar
// epilogue, as in the NEON/AVX2 twins); retained here in case a future 512-bit
// tier reinstates it.
#define ZSCORE8(i4_512, y4_512)                                         \
    {                                                                   \
        __m512i tmpi = _mm512_sub_epi8(i4_512, x512);                   \
        __m512i tmpj = _mm512_sub_epi8(y4_512, y512);                   \
        cmp = _mm512_cmpgt_epi8_mask(tmpi, tmpj);                       \
        score512 = _mm512_sub_epi8(maxScore512, maxRS1);                \
        __m512i insdel = _mm512_mask_blend_epi8(cmp, e_ins512, e_del512); \
        __m512i sub_a512 = _mm512_sub_epi8(tmpi, tmpj);                 \
        __m512i sub_b512 = _mm512_sub_epi8(tmpj, tmpi);                 \
        __m512i tmp1 = _mm512_mask_blend_epi8(cmp, sub_b512, sub_a512);         \
        tmp1 = _mm512_sub_epi8(score512, tmp1);                         \
        cmp = _mm512_cmpgt_epi8_mask(tmp1, zdrop512);                   \
        exit0 = _mm512_mask_blend_epi8(cmp, exit0, zero512);            \
    }



// ------------------------ vec 16 --------------------------------------------------
#define ZSCORE16(i4_512, y4_512)                                            \
    {                                                                   \
        __m512i tmpi = _mm512_sub_epi16(i4_512, x512);                  \
        __m512i tmpj = _mm512_sub_epi16(y4_512, y512);                  \
        cmp = _mm512_cmpgt_epi16_mask(tmpi, tmpj);                      \
        score512 = _mm512_sub_epi16(maxScore512, maxRS1);               \
        __m512i insdel = _mm512_mask_blend_epi16(cmp, e_ins512, e_del512); \
        __m512i sub_a512 = _mm512_sub_epi16(tmpi, tmpj);                    \
        __m512i sub_b512 = _mm512_sub_epi16(tmpj, tmpi);                    \
        __m512i tmp1 = _mm512_mask_blend_epi16(cmp, sub_b512, sub_a512);            \
        tmp1 = _mm512_sub_epi16(score512, tmp1);                            \
        cmp = _mm512_cmpgt_epi16_mask(tmp1, zdrop512);                  \
        exit0 = _mm512_mask_blend_epi16(cmp, exit0, zero512);           \
    }


// --- PR 17/16: AVX-512BW LUT primitive ---
// Broadcast the 16-byte pmat into all four 128-bit lanes.
// _mm512_shuffle_epi8 is per-128-bit-lane; each lane shuffles against
// its own quarter of the 512-bit pmat. The XOR index is order-insensitive,
// so it cannot represent an asymmetric matrix — that case uses SBT_PREPASS8_AMAT.
#define SBT_PREPASS8_XOR(s1, s2, sbt11_out, pmat512)                    \
    {                                                                   \
        __m512i xor_ = _mm512_xor_si512(s1, s2);                        \
        sbt11_out = _mm512_shuffle_epi8(pmat512, xor_);                 \
    }

// D3 generic-matrix seam (gated; symmetric default uses SBT_PREPASS8_XOR). Index
// the target-major LUT amat[(ref<<2)|read] (broadcast to all 4 lanes), then mask
// N to ambig. ACGT=0..3, N=4(target)/8(query). N detected by max_epu8(s1,s2) > 3.
// Mirrors the validated NEON/AVX2 SBT_PREPASS8_AMAT.
#define SBT_PREPASS8_AMAT(s1, s2, sbt11_out, amat512, ambig512, three512) \
    {                                                                   \
        __m512i sh_  = _mm512_add_epi8(s1, s1);                         \
        sh_          = _mm512_add_epi8(sh_, sh_);          /* s1 << 2 */ \
        __m512i idx_ = _mm512_or_si512(sh_, s2);           /* (ref<<2)|read */ \
        __m512i acgt_  = _mm512_shuffle_epi8(amat512, idx_);           \
        __m512i nmax_  = _mm512_max_epu8(s1, s2);                       \
        __mmask64 nmask_ = _mm512_cmpgt_epi8_mask(nmax_, three512); /* N: base > 3 */ \
        sbt11_out = _mm512_mask_blend_epi8(nmask_, acgt_, ambig512);    \
    }

// D3 8-bit rank-1 fast path: symmetric XOR LUT + single freed-to-match override.
// Cheaper than AMAT on the fast 8-bit tier (no index build, no N-mask). Mirrors NEON.
// rowfreed = cmpeq_epi8_mask(ref, fr_ref) hoisted per row (s1 constant per band).
#define SBT_PREPASS8_RANK1(s1, s2, rowfreed, sbt11_out, pmat512, match512, frread512) \
    {                                                                   \
        __m512i xor_  = _mm512_xor_si512(s1, s2);                       \
        __m512i sbt_  = _mm512_shuffle_epi8(pmat512, xor_);           \
        __mmask64 freed_ = rowfreed & _mm512_cmpeq_epi8_mask(s2, frread512); \
        sbt11_out = _mm512_mask_blend_epi8(freed_, sbt_, match512);    \
    }

#define MAIN_CODE8_CORE(sbt11, h00, h11, e11, f11, f21, zero512, e_ins512, oe_ins512, e_del512, oe_del512) \
    {                                                                   \
        /* M = max(0, h00 + sbt) in UNSIGNED-saturating form so a legitimate  \
         * score in [128,255] is kept (the old signed add_epi8 + signed       \
         * max_epi8 floor wrapped/mis-read >127 as negative). Split the signed \
         * substitution score into +bonus and -penalty parts: adds_epu8 (no   \
         * wrap past 255) then subs_epu8 (floors at 0). Mirrors the validated  \
         * NEON smithWaterman128_8 recurrence. */                              \
        __m512i sbt_pos = _mm512_max_epi8(sbt11, zero512);              \
        __m512i sbt_neg = _mm512_max_epi8(_mm512_sub_epi8(zero512, sbt11), zero512); \
        __m512i m11 = _mm512_subs_epu8(_mm512_adds_epu8(h00, sbt_pos), sbt_neg); \
        __mmask64 cmp11 = _mm512_cmpeq_epi8_mask(h00, zero512);         \
        m11 = _mm512_mask_blend_epi8(cmp11, m11, zero512);  /* h00==0 -> local restart */ \
        h11 = _mm512_max_epu8(m11, e11);                                \
        h11 = _mm512_max_epu8(h11, f11);                                \
        /* gap-open from H (standard Gotoh), unsigned-saturating */      \
        __m512i temp512 = _mm512_subs_epu8(h11, oe_ins512);            \
        e11 = _mm512_subs_epu8(e11, e_ins512);                          \
        e11 = _mm512_max_epu8(temp512, e11);                            \
        temp512 = _mm512_subs_epu8(h11, oe_del512);                    \
        f21 = _mm512_subs_epu8(f11, e_del512);                          \
        f21 = _mm512_max_epu8(temp512, f21);                            \
    }

// Symmetric (default) fast path. N detected by movepi16_mask (N=0xFFFF, high bit).
#define SBT_PREPASS16_SYM(s1, s2, sbt11_out, mismatch512, match512, w_ambig_512) \
    {                                                                   \
        __mmask32 cmp_ = _mm512_cmpeq_epi16_mask(s1, s2);               \
        __m512i sbt_ = _mm512_mask_blend_epi16(cmp_, mismatch512, match512); \
        __m512i tmp_ = _mm512_max_epu16(s1, s2);                        \
        __mmask32 amb_ = _mm512_movepi16_mask(tmp_);                    \
        sbt11_out = _mm512_mask_blend_epi16(amb_, sbt_, w_ambig_512);   \
    }

// D3 rank-1 fast path: symmetric + a single off-diagonal freed to a match
// (bisulfite OT/OB). Match when (ref==read) OR (ref==fr_ref AND read==fr_read);
// no LUT, no sign-extend. Mirrors NEON/AVX2.
// alt1 = (ref==fr_ref) ? fr_read : ref, hoisted per row (see NEON variant).
#define SBT_PREPASS16_RANK1(s1, s2, alt1, sbt11_out, mismatch512, match512, w_ambig_512) \
    {                                                                   \
        __mmask32 eq_  = _mm512_cmpeq_epi16_mask(s2, s1);             \
        __mmask32 fr_  = _mm512_cmpeq_epi16_mask(s2, alt1);          \
        __mmask32 ism_ = eq_ | fr_;                                    \
        __m512i sbt_ = _mm512_mask_blend_epi16(ism_, mismatch512, match512); \
        __mmask32 amb_ = _mm512_movepi16_mask(_mm512_max_epu16(s1, s2)); \
        sbt11_out = _mm512_mask_blend_epi16(amb_, sbt_, w_ambig_512);   \
    }

// D3 generic-matrix seam (gated; >=2 freed cells / non-match freed values). LUT
// amat[(ref<<2)|read] for both-ACGT lanes; non-ACGT is mismatch-or-ambig. Byte LUT
// score in low byte -> sign-extend. acgt mask via cmpeq(max_epu16(maxb,3),3).
#define SBT_PREPASS16_AMAT(s1, s2, sbt11_out, amat512, mismatch512, w_ambig_512, three512) \
    {                                                                   \
        __m512i maxb_ = _mm512_max_epu16(s1, s2);                      \
        __mmask32 nN_ = _mm512_movepi16_mask(maxb_);  /* N high bit */ \
        __m512i base_ = _mm512_mask_blend_epi16(nN_, mismatch512, w_ambig_512); \
        __mmask32 acgt_ = _mm512_cmpeq_epi16_mask(_mm512_max_epu16(maxb_, three512), three512); \
        __m512i idx_  = _mm512_or_si512(_mm512_slli_epi16(s1, 2), s2);  \
        __m512i lut_  = _mm512_shuffle_epi8(amat512, idx_);           \
        lut_ = _mm512_srai_epi16(_mm512_slli_epi16(lut_, 8), 8);      \
        sbt11_out = _mm512_mask_blend_epi16(acgt_, base_, lut_);      \
    }

#define MAIN_CODE16_CORE(sbt11, h00, h11, e11, f11, f21, zero512, e_ins512, oe_ins512, e_del512, oe_del512) \
    {                                                                   \
        __m512i m11 = _mm512_add_epi16(h00, sbt11);                     \
        __mmask32 cmp11 = _mm512_cmpeq_epi16_mask(h00, zero512);        \
        m11 = _mm512_mask_blend_epi16(cmp11, m11, zero512);             \
        h11 = _mm512_max_epi16(m11, e11);                               \
        h11 = _mm512_max_epi16(h11, f11);                               \
        __m512i temp512 = _mm512_sub_epi16(h11, oe_ins512);             \
        __m512i val512  = _mm512_max_epi16(temp512, zero512);           \
        e11 = _mm512_sub_epi16(e11, e_ins512);                          \
        e11 = _mm512_max_epi16(val512, e11);                            \
        temp512 = _mm512_sub_epi16(h11, oe_del512);                     \
        val512  = _mm512_max_epi16(temp512, zero512);                   \
        f21 = _mm512_sub_epi16(f11, e_del512);                          \
        f21 = _mm512_max_epi16(val512, f21);                            \
    }



inline void sortPairsId(SeqPair *pairArray, int32_t first, int32_t count,
                        SeqPair *tempArray)
{
    int32_t i;
    
    for(i = 0; i < count; i++)
    {
        SeqPair sp = pairArray[i];
        int32_t pos = sp.id - first;
        tempArray[pos] = sp;
    }

    for(i = 0; i < count; i++)
        pairArray[i] = tempArray[i];
}

// ____________________________ AVX512 - getScore() _______________________________________
#define PFD8 5
void BandedPairWiseSW::getScores8(SeqPair *pairArray,
                                  uint8_t *seqBufRef,
                                  uint8_t *seqBufQer,
                                  int32_t numPairs,
                                  uint16_t numThreads,
                                  int32_t w)
{
    assert(SIMD_WIDTH8 == 64 && SIMD_WIDTH16 == 32);
    int i;
    int64_t startTick, endTick;

    {
        BswOvershootGuard _g(pairArray, numPairs, SIMD_WIDTH8, guard_overshoot_);
        smithWatermanBatchWrapper8(pairArray, seqBufRef, seqBufQer, numPairs, numThreads, w);
    }

#if MAXI
    printf("AVX512/8 Vecor code: Writing output..\n");
    for (int l=0; l<numPairs; l++)
    {
        fprintf(stderr, "%d (%d %d) %d %d %d\n",
                pairArray[l].score, pairArray[l].tle, pairArray[l].qle,
                pairArray[l].gscore, pairArray[l].max_off, pairArray[l].gtle);

    }
    printf("Vector code: Writing output completed!!!\n\n");
#endif

}

void BandedPairWiseSW::smithWatermanBatchWrapper8(SeqPair *pairArray,
                                                  uint8_t *seqBufRef,
                                                  uint8_t *seqBufQer,
                                                  int32_t numPairs,
                                                  uint16_t numThreads,
                                                  int32_t w)
{
    numThreads = effective_threads(numThreads);
    int64_t st1, st2, st3, st4, st5;
#if RDT
    st1 = ___rdtsc();
#endif
    uint8_t *seq1SoA = (uint8_t *)_mm_malloc((size_t)MAX_SEQ_LEN8 * SIMD_WIDTH8 * numThreads * sizeof(uint8_t), 64);
    uint8_t *seq2SoA = (uint8_t *)_mm_malloc((size_t)MAX_SEQ_LEN8 * SIMD_WIDTH8 * numThreads * sizeof(uint8_t), 64);
    if (UNLIKELY(seq1SoA == NULL || seq2SoA == NULL)) {
        fprintf(stderr, "Error! Mem not allocated!!!\n");
        exit(EXIT_FAILURE);
    }

    int32_t ii;
    int32_t roundNumPairs = ((numPairs + SIMD_WIDTH8 - 1)/SIMD_WIDTH8 ) * SIMD_WIDTH8;
    for(ii = numPairs; ii < roundNumPairs; ii++)
    {
        pairArray[ii].id = ii;
        pairArray[ii].len1 = 0;
        pairArray[ii].len2 = pairArray[numPairs - 1].len2;
        pairArray[ii].idr = 0;
        pairArray[ii].idq = 0;
    }

#if RDT
    st2 = ___rdtsc();
#endif
    
    
#if RDT
    st3 = ___rdtsc();
#endif

    int eb = end_bonus;
//#pragma omp parallel num_threads(numThreads)
    {
        int32_t i;
        uint16_t tid = 0;
        uint8_t *mySeq1SoA = seq1SoA + tid * MAX_SEQ_LEN8 * SIMD_WIDTH8;
        uint8_t *mySeq2SoA = seq2SoA + tid * MAX_SEQ_LEN8 * SIMD_WIDTH8;
        uint8_t *seq1;
        uint8_t *seq2;
        uint8_t h0[SIMD_WIDTH8]   __attribute__((aligned(64)));
        uint8_t band[SIMD_WIDTH8];      
        uint8_t qlen[SIMD_WIDTH8] __attribute__((aligned(64)));
        int32_t bsize = 0;
        
        int8_t *H1 = H8_ + tid * SIMD_WIDTH8 * MAX_SEQ_LEN8;
        int8_t *H2 = H8__ + tid * SIMD_WIDTH8 * MAX_SEQ_LEN8;

        __m512i zero512   = _mm512_setzero_si512();
        __m512i o_ins512  = _mm512_set1_epi8(o_ins);
        __m512i e_ins512  = _mm512_set1_epi8(e_ins);
        __m512i oe_ins512 = _mm512_set1_epi8(o_ins + e_ins);
        __m512i o_del512  = _mm512_set1_epi8(o_del);
        __m512i e_del512  = _mm512_set1_epi8(e_del);
        __m512i eb_ins512 = _mm512_set1_epi8(eb - o_ins);
        __m512i eb_del512 = _mm512_set1_epi8(eb - o_del);
        
        int8_t max = 0;
        if (max < w_match) max = w_match;
        if (max < w_mismatch) max = w_mismatch;
        if (max < w_ambig) max = w_ambig;

        // Initial per-lane score floor B0 = max(0, h0 - REBASE_KEEP_W) (long-read
        // 8-bit, seed score h0 >= 128). Seeding the H bytes from raw h0 overflows
        // the signed-int8 score state on row 0, so we seed from the rebaselined
        // h0' = min(h0, REBASE_KEEP_W) instead (see the per-lane h0 init below).
        // REBASE_KEEP_W MUST equal the kernel's REBASE_KEEP (bsw8_rebase_keep)
        // so the two agree on B0; smithWaterman512_8 recomputes B0 via
        // bsw8_initial_floor(p[].h0, zdrop).
        const int REBASE_KEEP_W = bsw8_rebase_keep(zdrop);
        // The h0-prefix column/row seed below is unsigned-saturating [0,255], so
        // the seeded byte min(h0, REBASE_KEEP_W) only needs to fit a uint8. The
        // routing gate (bsw8_envelope_ok) enforces zdrop + maxStep <= 253 (so the
        // re-baseline window REBASE_HI = 255 - maxStep > REBASE_KEEP holds);
        // assert the byte bound for any direct caller that bypasses the gate.
        assert(REBASE_KEEP_W <= 255 &&
               "8-bit h0-prefix seed byte min(h0,REBASE_KEEP) must fit uint8 (zdrop<=254)");

        int nstart = 0, nend = numPairs;


//#pragma omp for schedule(dynamic, 128)
        for(i = nstart; i < nend; i+=SIMD_WIDTH8)
        {
            int32_t j, k;
            int maxLen1 = 0;
            int maxLen2 = 0;
            bsize = w;
            
            uint64_t tim;
            for(j = 0; j < SIMD_WIDTH8; j++)
            {
                if ((i + j + PFD8) < roundNumPairs) { // prefetch block (bounded; see getScores8/16 contract)
                    SeqPair spf = pairArray[i + j + PFD8];
                    _mm_prefetch((const char*) seqBufRef + (int64_t)spf.idr, _MM_HINT_NTA);
                    _mm_prefetch((const char*) seqBufRef + (int64_t)spf.idr + 64, _MM_HINT_NTA);
                }
                SeqPair sp = pairArray[i + j];
                // Re-baseline the seed score into the int8 byte frame: seed the H
                // arrays from h0' = h0 - B0 = min(h0, REBASE_KEEP_W) (clamped >= 0)
                // so the initial bytes fit signed-int8 even when the true seed h0
                // is several hundred (long-read). smithWaterman512_8 recomputes the
                // matching B0 = max(0, h0 - REBASE_KEEP) from p[].h0 + zdrop and
                // starts that lane's floor B there, so byte == H_absolute - B.
                {
                    int h0p = sp.h0;
                    if (h0p > REBASE_KEEP_W) h0p = REBASE_KEEP_W;
                    if (h0p < 0) h0p = 0;
                    // Always-on uint8 guard (asserts are compiled out in release):
                    // a direct caller that bypasses the gate with zdrop > 254 would
                    // otherwise truncate the seed byte. In-envelope REBASE_KEEP_W <= 253,
                    // so this never fires on the production path.
                    if (h0p > 255) h0p = 255;
                    h0[j] = (uint8_t) h0p;
                }
                seq1 = seqBufRef + (int64_t)sp.idr;

                for(k = 0; k < sp.len1; k++)
                {
                    mySeq1SoA[k * SIMD_WIDTH8 + j] = seq1[k] /* PR16: N stays 4 */;
                }
                qlen[j] = sp.len2 * max;
                if(maxLen1 < sp.len1) maxLen1 = sp.len1;
            }

            for(j = 0; j < SIMD_WIDTH8; j++)
            {
                SeqPair sp = pairArray[i + j];
                for(k = sp.len1; k <= maxLen1; k++)
                {
                    mySeq1SoA[k * SIMD_WIDTH8 + j] = DUMMY1;
                }
            }
            /* B5: only boundary row H2[maxLen1] survives the seed below; write just that row. */
            _mm512_store_si512((__m512i *)(H2 + maxLen1 * SIMD_WIDTH8), _mm512_set1_epi8((char)DUMMY1));
//--------------------
            __m512i h0_512 = _mm512_load_si512((__m512i*) h0);
            _mm512_store_si512((__m512i *) H2, h0_512);
            // h0-prefix deletion seed, unsigned-saturating [0,255] (was signed
            // sub_epi8 + max_epi8 floor; see smithWaterman128_8 / smithWaterman256_8).
            __m512i tmp512 = _mm512_subs_epu8(h0_512, o_del512);

            for(k = 1; k < maxLen1; k++) {
                tmp512 = _mm512_subs_epu8(tmp512, e_del512);
                _mm512_store_si512((__m512i *)(H2 + k* SIMD_WIDTH8), tmp512);
            }
//-------------------
            for(j = 0; j < SIMD_WIDTH8; j++)
            {
                if ((i + j + PFD8) < roundNumPairs) { // prefetch block (bounded; see getScores8/16 contract)
                    SeqPair spf = pairArray[i + j + PFD8];
                    _mm_prefetch((const char*) seqBufQer + (int64_t)spf.idq, _MM_HINT_NTA);
                    _mm_prefetch((const char*) seqBufQer + (int64_t)spf.idq + 64, _MM_HINT_NTA);
                }
                
                SeqPair sp = pairArray[i + j];
                seq2 = seqBufQer + (int64_t)sp.idq;
                for(k = 0; k < sp.len2; k++)
                {
                    mySeq2SoA[k * SIMD_WIDTH8 + j] = (seq2[k]==AMBIG ? 8 : seq2[k]) /* PR16: query N→8 */;
                }
                if(maxLen2 < sp.len2) maxLen2 = sp.len2;
            }
            
            for(j = 0; j < SIMD_WIDTH8; j++)
            {
                SeqPair sp = pairArray[i + j];
                for(k = sp.len2; k <= maxLen2; k++)
                {
                    mySeq2SoA[k * SIMD_WIDTH8 + j] = DUMMY2;
                }
            }
            /* B5: only boundary row H1[maxLen2]=0 survives the seed below. */
            _mm512_store_si512((__m512i *)(H1 + maxLen2 * SIMD_WIDTH8), _mm512_setzero_si512());
//------------------------
            _mm512_store_si512((__m512i *) H1, h0_512);
            // h0-prefix insertion seed, unsigned-saturating [0,255]:
            // H1[1] = max(0, h0' - oe_ins), then -e_ins per step (was signed
            // cmpgt_mask+sub+mask_blend, then sub+max_epi8).
            tmp512 = _mm512_subs_epu8(h0_512, oe_ins512);
            _mm512_store_si512((__m512i *) (H1 + SIMD_WIDTH8), tmp512);

            for(k = 2; k < maxLen2; k++)
            {
                tmp512 = _mm512_subs_epu8(tmp512, e_ins512);
                _mm512_store_si512((__m512i *)(H1 + k*SIMD_WIDTH8), tmp512);
            }
//------------------------
            /* Banding calculation in pre-processing */
            uint8_t myband[SIMD_WIDTH8] __attribute__((aligned(64)));
            uint8_t temp[SIMD_WIDTH8] __attribute__((aligned(64)));
            {
                __m512i qlen512 = _mm512_load_si512((__m512i *) qlen);
                __m512i sum512 = _mm512_add_epi8(qlen512, eb_ins512);
                _mm512_store_si512((__m512i *) temp, sum512);               
                for (int l=0; l<SIMD_WIDTH8; l++) {
                    double val = temp[l]/e_ins + 1.0;
                    int max_ins = val;
                    max_ins = max_ins > 1? max_ins : 1;
                    myband[l] = min_(bsize, max_ins);
                }
                sum512 = _mm512_add_epi8(qlen512, eb_del512);
                _mm512_store_si512((__m512i *) temp, sum512);               
                for (int l=0; l<SIMD_WIDTH8; l++) {
                    double val = temp[l]/e_del + 1.0;
                    int max_ins = val;
                    max_ins = max_ins > 1? max_ins : 1;
                    myband[l] = min_(myband[l], max_ins);
                    bsize = bsize < myband[l] ? myband[l] : bsize;
                }
            }

            smithWaterman512_8(mySeq1SoA,
                               mySeq2SoA,
                               maxLen1,
                               maxLen2,
                               pairArray + i,
                               h0,
                               tid,
                               numPairs,
                               zdrop,
                               bsize,
                               myband);
        }
    }

#if RDT 
    st4 = ___rdtsc();
#endif
    

#if RDT 
    st5 = ___rdtsc();
    setupTicks = st2 - st1;
    sort1Ticks = st3 - st2;
    swTicks = st4 - st3;
    sort2Ticks = st5 - st4;
#endif
    
    // free mem
    _mm_free(seq1SoA);
    _mm_free(seq2SoA);

    return;
}

void BandedPairWiseSW::smithWaterman512_8(uint8_t seq1SoA[],
                                          uint8_t seq2SoA[],
                                          int nrow,
                                          int ncol,
                                          SeqPair *p,
                                          uint8_t h0[],
                                          uint16_t tid,
                                          int32_t numPairs,
                                          int zdrop,
                                          int32_t w,
                                          uint8_t myband[])
{
    __m512i match512     = _mm512_set1_epi8(this->w_match);
    __m512i mismatch512  = _mm512_set1_epi8(this->w_mismatch);
    __m512i gapOpen512   = _mm512_set1_epi8(this->w_open);
    __m512i gapExtend512 = _mm512_set1_epi8(this->w_extend);
    __m512i gapOE512     = _mm512_set1_epi8(this->w_open + this->w_extend);
    __m512i w_ambig_512  = _mm512_set1_epi8(this->w_ambig); // ambig penalty
    __m512i five512      = _mm512_set1_epi8(5);

    // PR 16: pmat LUT broadcast into all 4 128-bit lanes of 512-bit register.
    int8_t pmat_bytes[16] __attribute__((aligned(16)));
    build_pmat16(pmat_bytes, this->w_match, this->w_mismatch, this->w_ambig);
    __m128i pmat128 = _mm_load_si128((__m128i *)pmat_bytes);
    __m512i pmat512 = _mm512_broadcast_i32x4(pmat128);
    // D3 generic-matrix seam: symmetric default uses the XOR pmat; an asymmetric
    // matrix (bisulfite OT/OB) uses the target-major amat LUT.
    const bool forced  = bsw_force_generic_matrix();
    const bool gen_mat = bsw_generic_matrix(this->mat, this->w_match, this->w_mismatch)
                         || forced;
    const BswFreedCell fc = bsw_freed_cell(this->mat, this->w_match, this->w_mismatch, forced);
    __m512i frref512  = _mm512_set1_epi8(fc.ref);
    __m512i frread512 = _mm512_set1_epi8(fc.read);
    int8_t amat_bytes[16] __attribute__((aligned(16)));
    build_amat16(amat_bytes, this->mat);
    __m512i amat512 = _mm512_broadcast_i32x4(_mm_load_si128((__m128i *)amat_bytes));
    __m512i three512_8 = _mm512_set1_epi8(3);

    __m512i e_del512    = _mm512_set1_epi8(this->e_del);
    __m512i oe_del512   = _mm512_set1_epi8(this->o_del + this->e_del);
    __m512i e_ins512    = _mm512_set1_epi8(this->e_ins);
    __m512i oe_ins512   = _mm512_set1_epi8(this->o_ins + this->e_ins);

    int8_t  *F   = F8_ + tid * SIMD_WIDTH8 * MAX_SEQ_LEN8;
    int8_t  *H_h = H8_ + tid * SIMD_WIDTH8 * MAX_SEQ_LEN8;
    int8_t  *H_v = H8__ + tid * SIMD_WIDTH8 * MAX_SEQ_LEN8;

    int i, j;

    uint8_t tail[SIMD_WIDTH8] __attribute((aligned(64)));
    uint8_t head[SIMD_WIDTH8] __attribute((aligned(64)));

    // PR 17: per-row score-vector scratch for fission. Thread-local view
    // into sbt8_ (slab-backed; see constructor) — too large for the stack.
    int8_t *sbt_buf = sbt8_ + tid * SIMD_WIDTH8 * MAX_SEQ_LEN8;

    // --- DIAGONAL-OFFSET POSITION ENCODING (long-read 8-bit, w<=127) ---
    // Every per-cell COLUMN position is tracked as the diagonal offset
    //   d = col - i  in [-w, +w+1]  (fits signed int8 for w <= ~126).
    // ROW quantities (best row, best-gscore row, qlen, tlen, mlen) exceed
    // int8 for long reads, so they live in WIDE per-lane int32 side channels
    // updated O(rows) in the per-row epilogue, not O(cells). Absolute end
    // coordinates are reconstructed at the result store from the wide row.
    // Persistent column-offset state (head512/tail512) is shifted by -1 each
    // row (frame follows i) so the same absolute edge keeps its offset.
    int32_t tlenw[SIMD_WIDTH8];   // raw target length (rows), wide
    int32_t qlenw[SIMD_WIDTH8];   // raw query length (cols), wide
    int32_t mbandw[SIMD_WIDTH8];  // per-lane band width, wide
    int32_t mlenw[SIMD_WIDTH8];   // min(qlen+myband, tlen), wide row bound
    int32_t xrow[SIMD_WIDTH8];    // best row for score (== i+1 at capture)
    int32_t ierow[SIMD_WIDTH8];   // best row for gscore (== i+1 at capture)

    // --- PER-LANE SCORE RE-BASELINING (long-read 8-bit, scores > 127) ---
    // H/F/E and the score maxima are unsigned 8-bit [0,255], but a 1 kb
    // alignment scores ~900. To keep the kernel 64-lane (8-bit) we carry a
    // wide per-lane running FLOOR B[l]: every stored 8-bit score represents
    // (H_absolute - B[l]). When a lane's row-max climbs toward 255 we raise
    // B[l] by `delta` and subtract `delta` (saturating, _mm512_subs_epu8) from
    // every carried score buffer/register for that lane, then add B back when
    // reconstructing the absolute result.
    //
    // Scores live in the UNSIGNED byte range [0,255]. The DP recurrence computes
    // M = max(0, h00 + sbt) with unsigned-saturating arithmetic (adds_epu8 then
    // subs_epu8) and the row/global-max trackers (maxRS1, maxScore512) compare
    // with UNSIGNED order (_mm512_cmpgt_epu8_mask), so the full [0,255] is usable.
    //
    // CORRECTNESS depends on re-baselining NEVER firing for a pair admitted to
    // this kernel. The saturating-subtract is NOT generally lossless: it can zero
    // a still-positive off-diagonal cell (a cell may sit > zdrop below the ROW max
    // yet still lie on the eventual optimum — z-drop is a row-level early-exit, not
    // a per-cell guarantee), and a zeroed cell is then misread as the h00==0
    // local-restart sentinel, spuriously truncating a valid alignment. So we do
    // NOT rely on re-baseline being lossless; instead it is held inert:
    //   REBASE_HI = 255 - maxStep: a row whose max has not reached REBASE_HI grows
    //     by <= maxStep per row (H(i,j) <= H(i-1,j-1) + maxStep), so its bytes stay
    //     <= 255 with no rebase needed.
    //   The bwamem.cpp routing gate admits a pair only when its MAX ATTAINABLE
    //     score stays < REBASE_HI and h0 <= REBASE_KEEP (so B0 == 0). Under that
    //     gate the row max never reaches REBASE_HI: B stays 0, stored bytes equal
    //     absolute scores, and the kernel is a plain exact unsigned [0,255] SW.
    // The re-baseline machinery below is retained only as a safety net; pairs that
    // could exceed the envelope take the 16-bit path.
    // Precondition REBASE_HI > REBASE_KEEP (a non-empty window) holds for
    // zdrop + maxStep <= 253, which is exactly what the routing gate enforces.
    // The h0-prefix column/row seed (wrapper setup below) is unsigned-saturating
    // [0,255] and imposes no tighter ceiling; it previously used signed int8 ops
    // that required the seed byte <= 127 and capped zdrop at 126.
    int maxStep = (int)this->w_match;
    if ((int)this->w_ambig > maxStep) maxStep = (int)this->w_ambig;
    if (maxStep < 1) maxStep = 1;
    const int REBASE_KEEP = zdrop + 1;
    const int REBASE_HI   = 255 - maxStep;
    assert(REBASE_HI > REBASE_KEEP &&
           "8-bit SW re-baseline: zdrop+maxStep>253, scoring outside safe envelope (route to 16-bit)");
    int32_t B[SIMD_WIDTH8];        // per-lane wide score floor (stored = abs - B)
    int32_t best_abs[SIMD_WIDTH8]; // running best score in ABSOLUTE units
    int32_t gbest_abs[SIMD_WIDTH8];// running gscore (query-end) in ABSOLUTE units

    // --- INITIAL SCORE FLOOR B0 (long-read 8-bit, seed score h0 >= 128) ---
    // The H arrays are seeded from the seed score h0 (H_v[0]=h0, then gap-
    // decremented down column 0 and across row 0). For a long-read seed h0 can
    // be several hundred, which would overflow the unsigned [0,255] byte state
    // before any re-baseline event can fire. So we start each lane's
    // running floor at B0[l] = max(0, h0 - REBASE_KEEP) instead of 0: the
    // wrapper seeds the H bytes from h0' = h0 - B0 = min(h0, REBASE_KEEP)
    // (clamped >= 0, so the initial bytes land in [0, REBASE_KEEP] <= 254), and
    // here we set B[l] = B0[l] so the stored bytes still mean (H_absolute - B).
    //
    // SAFETY (identical argument to the mid-DP re-baseline): the gap-prefix
    // cells in column 0 / row 0 that fall more than REBASE_KEEP below h0 are
    // saturated to byte 0 by the wrapper's unsigned-saturating seed. Because
    // REBASE_KEEP = zdrop+1 > zdrop, any such prefix cell is already beyond the
    // z-drop horizon below the seed, so an alignment routed through it has
    // already z-dropped and cannot be optimal -- clamping it to the floor is
    // lossless. This is the SAME situation that arises mid-DP whenever B>0 after
    // a re-baseline (byte 0 then means H_absolute == B, not 0), which the
    // h00==0 local-restart test in MAIN_CODE8_CORE already handles correctly; B0
    // simply makes B>0 from row 0 instead of from the first re-baseline event.
    int32_t B0[SIMD_WIDTH8];
    for (int l=0; l<SIMD_WIDTH8; l++) {
        B0[l] = bsw8_initial_floor((int)p[l].h0, zdrop);
    }

    int32_t minq = 10000000;
    for (int l=0; l<SIMD_WIDTH8; l++) {
        tlenw[l]  = p[l].len1;
        qlenw[l]  = p[l].len2;
        mbandw[l] = myband[l];
        int ml = qlenw[l] + mbandw[l];
        if (ml > tlenw[l]) ml = tlenw[l];
        mlenw[l]  = ml;
        xrow[l]   = 0;
        ierow[l]  = 0;
        B[l]        = B0[l];   // start the floor at B0 so byte state = abs - B0
        best_abs[l] = p[l].h0; // maxScore512 inits to the h0 seed; record the
                               // ABSOLUTE seed score (wide), not the rebaselined byte
        gbest_abs[l]= -1;      // unset sentinel (-1): gscore=-1 / gtle=0 when no query end is reached, matching scalar
        if (p[l].len2 < minq) minq = p[l].len2;
    }
    minq -= 1; // for gscore

    __m512i myband512 = _mm512_load_si512((__m512i *) myband);
    __m512i zero512   = _mm512_setzero_si512();
    __m512i one512    = _mm512_set1_epi8(1);
    __m512i two512    = _mm512_set1_epi8(2);
    __m512i ff512     = _mm512_set1_epi8(0xFF);

    // gscore query-end capture (see smithWaterman128_8). Reset per row.
    __m512i hqe512   = zero512;   // query-end cell H (rebaselined byte) this row
    __m512i qfire512 = zero512;   // 0xFF where this lane reached its query end this row

    // Offset-frame band edges. head_off starts at 0 (col 0 - row 0); tail_off
    // starts saturated-high (+127) and is immediately clamped by the band-grow
    // min() against (1+myband) and (qlen-i) on the first row.
    __m512i head512 = zero512;
    __m512i tail512 = _mm512_set1_epi8(127);
    _mm512_store_si512((__m512i *) head, head512);
    _mm512_store_si512((__m512i *) tail, tail512);

    __m512i hval = _mm512_load_si512((__m512i *)(H_v));
    __mmask64 dmask = 0xFFFFFFFFFFFFFFFFULL;

    __m512i maxScore512 = hval;
    for(j = 0; j < ncol; j++)
        _mm512_store_si512((__m512i *)(F + j * SIMD_WIDTH8), zero512);

    __m512i y512       = zero512;   // best col as diagonal offset: col - i_at_capture (i_at_capture = xrow[l]-1)
    __m512i max_off512 = zero512;
    __m512i exit0      = _mm512_set1_epi8(0xFF);
    __m512i zdrop512   = _mm512_set1_epi8(zdrop);

    int beg = 0, end = ncol;
    int nbeg = beg, nend = end;

#if RDT
    uint64_t tim = __rdtsc();
#endif

    for(i = 0; i < nrow; i++)
    {
        __m512i e11 = zero512;
        __m512i h00, h11, h10;
        __m512i s10 = _mm512_load_si512((__m512i *)(seq1SoA + (i + 0) * SIMD_WIDTH8));

        beg = nbeg; end = nend;
        // Banding
        if (beg < i - w) beg = i - w;
        if (end > i + w + 1) end = i + w + 1;
        if (end > ncol) end = ncol;

        h10 = zero512;
        if (beg == 0)
            h10 = _mm512_load_si512((__m512i *)(H_v + (i+1) * SIMD_WIDTH8));

        __m512i j512 = zero512;
        __m512i maxRS1 = zero512;

        __m512i y1_512 = zero512;   // row-max column as diagonal offset (col - i)

        // gscore query-end capture resets each row.
        hqe512   = zero512;
        qfire512 = zero512;

        // Per-row diagonal-offset of the query end: qlen_off = qlen - i. Built
        // wide then saturated to int8, with a validity mask so an out-of-band
        // qlen-i (which would wrap and spuriously alias an in-band col offset)
        // never triggers a false gscore/clamp. qlen_off_valid is true only when
        // qlen-i lies within the representable band window [-w, w+1].
        int8_t qlen_off_a[SIMD_WIDTH8]   __attribute((aligned(64)));
        int8_t qlen_valid_a[SIMD_WIDTH8] __attribute((aligned(64)));
        int8_t cmpim_a[SIMD_WIDTH8]      __attribute((aligned(64)));
        for (int l=0; l<SIMD_WIDTH8; l++) {
            int qoff = qlenw[l] - i;                 // wide
            int qoff_sat = qoff;                     // saturate to int8 range
            if (qoff_sat >  127) qoff_sat =  127;
            if (qoff_sat < -128) qoff_sat = -128;
            qlen_off_a[l]   = (int8_t) qoff_sat;
            qlen_valid_a[l] = (qoff >= -w && qoff <= w + 1) ? (int8_t)0xFF : 0;
            // exit when row i+1 has passed the per-lane effective length bound.
            cmpim_a[l]      = ((i + 1) > mlenw[l]) ? (int8_t)0xFF : 0;
        }
        __m512i qlen_off512   = _mm512_load_si512((__m512i *) qlen_off_a);
        // qlen_valid as a mask register (high-bit set => in-band query end).
        __mmask64 qlen_valid_k = _mm512_movepi8_mask(
                                    _mm512_load_si512((__m512i *) qlen_valid_a));

#if RDT
        uint64_t tim1 = __rdtsc();
#endif

        // Banding (diagonal-offset frame). head512/tail512 arrive here in row-i's
        // offset frame: the -1 shift at the end of the previous iteration converts
        // (col - (i-1)) -> (col - i), so no additional adjustment is needed here.
        //   abs head-grow: head = max(head, i - myband)  ->  head_off = max(head_off, -myband)
        //   abs tail-clamp: tail = min(tail, (i+1)+myband, qlen)
        //                                          -> tail_off = min(tail_off, 1+myband, qlen-i)
        __m512i cache512;
        __m512i phead512 = head512, ptail512 = tail512;
        __m512i negband512 = _mm512_sub_epi8(zero512, myband512);            // -myband
        head512 = _mm512_max_epi8(head512, negband512);
        cache512 = _mm512_add_epi8(myband512, one512);                       // 1 + myband
        tail512 = _mm512_min_epi8(tail512, cache512);
        tail512 = _mm512_min_epi8(tail512, qlen_off512);

        // NEW, trimming.
        __mmask64 cmph = _mm512_cmpeq_epi8_mask(head512, phead512);
        __mmask64 cmpt = _mm512_cmpeq_epi8_mask(tail512, ptail512);
        cmph &= cmpt;

        for (int l=beg; l<end && cmph != dmask; l++)
        {
            __m512i h512 = _mm512_load_si512((__m512i *)(H_h + l * SIMD_WIDTH8));
            __m512i f512 = _mm512_load_si512((__m512i *)(F + l * SIMD_WIDTH8));

            __m512i pj512 = _mm512_set1_epi8(l - i);   // diagonal offset of column l
            __mmask64 cmp1 = _mm512_cmpgt_epi8_mask(head512, pj512);
            if (cmp1 == 0x00) break;
            __mmask64 cmp2 = _mm512_cmpgt_epi8_mask(pj512, tail512);
            cmp1 = cmp1 | cmp2;
            h512 = _mm512_mask_blend_epi8(cmp1, h512, zero512);
            f512 = _mm512_mask_blend_epi8(cmp1, f512, zero512);

            _mm512_store_si512((__m512i *)(F + l * SIMD_WIDTH8), f512);
            _mm512_store_si512((__m512i *)(H_h + l * SIMD_WIDTH8), h512);
        }

#if RDT
        prof[DP3][0] += __rdtsc() - tim1;
#endif

        // cmpim: lane exits if row i+1 passed its effective length (precomputed
        // wide into cmpim_a) OR the band collapsed (tail<=head, offset frame).
        __mmask64 cmpim = _mm512_movepi8_mask(_mm512_load_si512((__m512i *) cmpim_a));
        __mmask64 cmpht = _mm512_cmpeq_epi8_mask(tail512, head512);
        cmpim = cmpim | cmpht;
        // NEW
        cmpht = _mm512_cmpgt_epi8_mask(head512, tail512);
        cmpim = cmpim |  cmpht;

        exit0 = _mm512_mask_blend_epi8(cmpim, exit0, zero512);

#if RDT
        tim1 = __rdtsc();
#endif

        // PR 17+16: AVX-512 8-bit score pre-pass via LUT. gen_mat branch is
        // loop-invariant: symmetric keeps the cheap XOR LUT.
        if (!gen_mat) {
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

        j512 = _mm512_set1_epi8(beg - i);   // diagonal offset of first band column
        for(j = beg; j < end; j++)
        {
            __m512i f11, f21, f31, f41, f51, jj512, sbt11;
            h00 = _mm512_load_si512((__m512i *)(H_h + j * SIMD_WIDTH8));
            f11 = _mm512_load_si512((__m512i *)(F + j * SIMD_WIDTH8));
            sbt11 = _mm512_load_si512((__m512i *)(sbt_buf + j * SIMD_WIDTH8));

            __m512i pj512 = j512;
            j512 = _mm512_add_epi8(j512, one512);
            
            MAIN_CODE8_CORE(sbt11, h00, h11, e11, f11, f21, zero512,
                            e_ins512, oe_ins512,
                            e_del512, oe_del512);

            // Masked writing
            __mmask64 cmp2 = _mm512_cmpgt_epi8_mask(head512, pj512);
            __mmask64 cmp1 = _mm512_cmpgt_epi8_mask(pj512, tail512);
            cmp1 = cmp1 | cmp2;
            h10 = _mm512_mask_blend_epi8(cmp1, h10, zero512);
            f21 = _mm512_mask_blend_epi8(cmp1, f21, zero512);
            
            /* Part of main code MAIN_CODE */
            __m512i bmaxRS = maxRS1, blend512;
            maxRS1 =_mm512_max_epu8(maxRS1, h11);
            // UNSIGNED >: signed cmpgt_epi8 mis-read scores >127 (long reads).
            __mmask64 cmpA = _mm512_cmpgt_epu8_mask(maxRS1, bmaxRS);
            __mmask64 cmpB =_mm512_cmpeq_epi8_mask(maxRS1, h11);                    
            cmpA = cmpA | cmpB;
            cmp1 = _mm512_cmpgt_epi8_mask(j512, tail512);
            cmp1 = cmp1 | cmp2;
            blend512 = _mm512_mask_blend_epi8(cmpA, y1_512, j512);
            y1_512 = _mm512_mask_blend_epi8(cmp1, blend512, y1_512);
            maxRS1 = _mm512_mask_blend_epi8(cmp1, maxRS1, bmaxRS);                      

            _mm512_store_si512((__m512i *)(F + j * SIMD_WIDTH8), f21);
            _mm512_store_si512((__m512i *)(H_h + j * SIMD_WIDTH8), h10);

            h10 = h11;
                        
            // gscore query-end capture (see smithWaterman128_8). Fire exactly like
            // the byte-identical 16-bit tier: col == qlen-1 (j512 == qlen_off) AND
            // band-grown tail reached the query end (tail512 == qlen_off == scalar's
            // end == qlen) AND in-band (qlen_valid) AND lane alive (exit0 high bit).
            // mask_blend(k, a, b) selects b where k, a where ~k.
            // gtle CONTRACT (see smithWaterman128_8): exact vs scalar for
            // gscore > 0; may differ only in the unused gscore == 0 tail.
            if (j >= minq)
            {
                __mmask64 cmp = _mm512_cmpeq_epi8_mask(j512, qlen_off512);
                cmp = cmp & _mm512_cmpeq_epi8_mask(tail512, qlen_off512);
                cmp = cmp & qlen_valid_k;
                cmp = cmp & _mm512_movepi8_mask(exit0);
                hqe512   = _mm512_mask_blend_epi8(cmp, hqe512, h11);
                qfire512 = _mm512_mask_blend_epi8(cmp, qfire512, ff512);
            }
        }
        __mmask64 cmp1 = _mm512_cmpgt_epi8_mask(head512, j512);
        __mmask64 cmp2 = _mm512_cmpgt_epi8_mask(j512, tail512);
        cmp1 = cmp1 | cmp2;
        h10 = _mm512_mask_blend_epi8(cmp1, h10, zero512);

        _mm512_store_si512((__m512i *)(H_h + j * SIMD_WIDTH8), h10);
        _mm512_store_si512((__m512i *)(F + j * SIMD_WIDTH8), zero512);

        /* exit due to zero score by a row */
        __m512i bmaxScore512 = maxScore512;
        __mmask64 tmp = _mm512_cmpeq_epi8_mask(maxRS1, zero512);
        if (tmp == dmask) {
            /* Finalize THIS row's query-end (gscore/gtle) capture before exiting
             * (see smithWaterman128_8 zero-row break): scalar records the row's
             * gscore then hits its own m==0 break, so the vector must too, else the
             * gscore==0 query-end tail row is dropped -> wrong gtle under the
             * asymmetric --meth matrix. Mirror of the epilogue gscore block. */
            int8_t hq_a_[SIMD_WIDTH8] __attribute((aligned(64)));
            _mm512_store_si512((__m512i *) hq_a_, hqe512);
            const __mmask64 qf64_ = _mm512_movepi8_mask(qfire512);
            for (int g = 0; g < SIMD_WIDTH8 / 16; g++) {
                const int base = g * 16;
                const __mmask16 qfm = (__mmask16)(qf64_ >> (16 * g));
                __m512i Bg   = _mm512_loadu_si512((const void *)(B + base));
                __m512i hqeg = _mm512_cvtepu8_epi32(_mm_loadu_si128((const __m128i *)(hq_a_ + base)));
                __m512i gs_abs = _mm512_add_epi32(hqeg, Bg);
                __m512i gba = _mm512_loadu_si512((const void *)(gbest_abs + base));
                __mmask16 gmask = qfm & _mm512_cmpge_epi32_mask(gs_abs, gba);
                gba = _mm512_mask_mov_epi32(gba, gmask, gs_abs);
                _mm512_storeu_si512((void *)(gbest_abs + base), gba);
                __m512i ierg = _mm512_loadu_si512((const void *)(ierow + base));
                ierg = _mm512_mask_mov_epi32(ierg, gmask, _mm512_set1_epi32(i + 1));
                _mm512_storeu_si512((void *)(ierow + base), ierg);
            }
            break;
        }

        exit0 = _mm512_mask_blend_epi8(tmp, exit0, zero512);

        __m512i score512 = _mm512_max_epu8(maxScore512, maxRS1);
        // exit0: 0xFF=alive; movepi8_mask high-bit => mask bit=1 where alive.
        __mmask64 mex0 = _mm512_movepi8_mask(exit0);
        maxScore512 = _mm512_mask_blend_epi8(mex0, maxScore512, score512);

        // UNSIGNED >: signed cmpgt_epi8 mis-read scores >127 (long reads).
        __mmask64 cmp = _mm512_cmpgt_epu8_mask(maxScore512, bmaxScore512);
        // y512 (best col) stays a diagonal offset captured in the best row's
        // frame; the best row itself moves to the wide xrow[] side channel.
        y512 = _mm512_mask_blend_epi8(cmp, y512, y1_512);

        // max_off = max running diagonal-distance of the row-max from the main
        // diagonal: |y1col - (i+1)| = |y1_off - 1| in the offset frame.
        // y1_off - 1 is a SIGNED int8 (negative when row-max is left of the
        // sub-diagonal), so use _mm512_abs_epi8.
        __m512i y1_minus1 = _mm512_sub_epi8(y1_512, one512);
        __m512i ind512 = _mm512_abs_epi8(y1_minus1);          // |y1_off - 1|
        __m512i bmax_off512 = max_off512;
        ind512 = _mm512_max_epu8(max_off512, ind512);
        max_off512 = _mm512_mask_blend_epi8(cmp, bmax_off512, ind512);

        // Per-lane wide updates (O(rows)): best-score row (xrow), best-gscore
        // row (ierow), absolute-frame score tracking, and the z-drop test — all
        // done in wide scalars so row distances that exceed int8 for long reads
        // are handled exactly.
        {
            // Only the int32 DATA channels need materializing as byte arrays; the
            // cmp/qfire/exit per-lane predicates are read straight from the cmp
            // __mmask64 and movepi8_mask(qfire512)/movepi8_mask(exit0) below.
            int8_t  y1_a[SIMD_WIDTH8]       __attribute((aligned(64)));
            int8_t  y_a[SIMD_WIDTH8]        __attribute((aligned(64)));
            int8_t  ms_a[SIMD_WIDTH8]       __attribute((aligned(64)));
            int8_t  rs_a[SIMD_WIDTH8]       __attribute((aligned(64)));
            int8_t  hqe_a[SIMD_WIDTH8]      __attribute((aligned(64)));
            _mm512_store_si512((__m512i *) y1_a, y1_512);
            _mm512_store_si512((__m512i *) y_a, y512);
            _mm512_store_si512((__m512i *) ms_a, maxScore512);
            _mm512_store_si512((__m512i *) rs_a, maxRS1);
            _mm512_store_si512((__m512i *) hqe_a, hqe512);
            // VECTORIZED per-lane wide updates (see smithWaterman128_8) using
            // AVX-512 mask registers: 16 int32 per __m512i -> SIMD_WIDTH8/16 groups
            // of 16. Per-lane predicates (cmp/qfire/exit) are __mmask16 slices; the
            // z-dropped lanes are accumulated into a combined mask and cleared in
            // exit0. Byte-identical (same >, >= tie-breaks; same abs/z-drop math).
            const __m512i vi   = _mm512_set1_epi32(i);
            const __m512i vip1 = _mm512_set1_epi32(i + 1);
            const __m512i vone = _mm512_set1_epi32(1);
            const __m512i vzd  = _mm512_set1_epi32(zdrop);
            const __m512i vedel = _mm512_set1_epi32(this->e_del);
            const __m512i veins = _mm512_set1_epi32(this->e_ins);
            const __mmask64 qf64 = _mm512_movepi8_mask(qfire512);
            const __mmask64 ex64 = _mm512_movepi8_mask(exit0);
            __mmask64 die64 = 0;
            for (int g = 0; g < SIMD_WIDTH8 / 16; g++) {
                const int base = g * 16;
                const __mmask16 cmpm = (__mmask16)(cmp  >> (16 * g));
                const __mmask16 qfm  = (__mmask16)(qf64 >> (16 * g));
                const __mmask16 exm  = (__mmask16)(ex64 >> (16 * g));
                __m512i Bg   = _mm512_loadu_si512((const void *)(B + base));
                __m512i msg  = _mm512_cvtepu8_epi32(_mm_loadu_si128((const __m128i *)(ms_a  + base)));
                __m512i rsg  = _mm512_cvtepu8_epi32(_mm_loadu_si128((const __m128i *)(rs_a  + base)));
                __m512i hqeg = _mm512_cvtepu8_epi32(_mm_loadu_si128((const __m128i *)(hqe_a + base)));
                __m512i y1g  = _mm512_cvtepi8_epi32(_mm_loadu_si128((const __m128i *)(y1_a  + base)));
                __m512i yg   = _mm512_cvtepi8_epi32(_mm_loadu_si128((const __m128i *)(y_a   + base)));

                __m512i xrg = _mm512_loadu_si512((const void *)(xrow + base));
                xrg = _mm512_mask_mov_epi32(xrg, cmpm, vip1);
                _mm512_storeu_si512((void *)(xrow + base), xrg);

                __m512i bag = _mm512_loadu_si512((const void *)(best_abs + base));
                bag = _mm512_max_epi32(bag, _mm512_add_epi32(msg, Bg));
                _mm512_storeu_si512((void *)(best_abs + base), bag);

                __m512i gs_abs = _mm512_add_epi32(hqeg, Bg);
                __m512i gba = _mm512_loadu_si512((const void *)(gbest_abs + base));
                __mmask16 gmask = qfm & _mm512_cmpge_epi32_mask(gs_abs, gba); // qfire & gs>=gbest
                gba = _mm512_mask_mov_epi32(gba, gmask, gs_abs);
                _mm512_storeu_si512((void *)(gbest_abs + base), gba);
                __m512i ierg = _mm512_loadu_si512((const void *)(ierow + base));
                ierg = _mm512_mask_mov_epi32(ierg, gmask, vip1);
                _mm512_storeu_si512((void *)(ierow + base), ierg);

                __m512i y1c  = _mm512_add_epi32(y1g, vi);
                __m512i yc   = _mm512_add_epi32(yg, _mm512_sub_epi32(xrg, vone));
                __m512i tmpi = _mm512_sub_epi32(vip1, xrg);
                __m512i tmpj = _mm512_sub_epi32(y1c, yc);
                // z-drop gap term weighted by gap-extend penalty, matching the
                // scalar reference (drift>0 -> deletion side *e_del, else *e_ins).
                __m512i zdelta = _mm512_sub_epi32(tmpi, tmpj);
                __mmask16 zdsel = _mm512_cmpgt_epi32_mask(zdelta, _mm512_setzero_si512());
                __m512i zesel  = _mm512_mask_blend_epi32(zdsel, veins, vedel);
                __m512i dif  = _mm512_mullo_epi32(_mm512_abs_epi32(zdelta), zesel);
                __m512i drop = _mm512_sub_epi32(msg, rsg);
                __mmask16 diem = exm & _mm512_cmpgt_epi32_mask(_mm512_sub_epi32(drop, dif), vzd);
                die64 |= ((__mmask64)diem) << (16 * g);
            }
            exit0 = _mm512_mask_mov_epi8(exit0, die64, _mm512_setzero_si512());
        }

#if RDT
        prof[DP1][0] += __rdtsc() - tim1;
#endif

        /* Narrowing of the band */
        /* From beg */
        int l;
        for (l = beg; l < end; l++)
        {
            __m512i f512 = _mm512_load_si512((__m512i *)(F + l * SIMD_WIDTH8));
            __m512i h512 = _mm512_load_si512((__m512i *)(H_h + l * SIMD_WIDTH8));
            __m512i tmp = _mm512_or_si512(f512, h512);
            __mmask64 val = _mm512_cmpeq_epi8_mask(tmp, zero512);
            if (val == dmask) nbeg = l;
            else
                break;
        }

        /* From end */
        for (l = end; l >= beg; l--)
        {
            __m512i f512 = _mm512_load_si512((__m512i *)(F + l * SIMD_WIDTH8));
            __m512i h512 = _mm512_load_si512((__m512i *)(H_h + l * SIMD_WIDTH8));
            __m512i tmp = _mm512_or_si512(f512, h512);
            __mmask64 val = _mm512_cmpeq_epi8_mask(tmp, zero512);
            if (val != dmask)
                break;
        }
        nend = l + 2 < ncol? l + 2: ncol;

#if RDT
        tim1 = __rdtsc();
#endif

        __m512i exit1 = _mm512_xor_si512(exit0, ff512);
        __mmask64 tmpb = dmask;
        __m512i l512 = _mm512_set1_epi8(beg - i);   // diagonal offset of column beg

        for (l = beg; l < end; l++)
        {
            __m512i f512 = _mm512_load_si512((__m512i *)(F + l * SIMD_WIDTH8));
            __m512i h512 = _mm512_load_si512((__m512i *)(H_h + l * SIMD_WIDTH8));
            __m512i tmp_ = _mm512_or_si512(f512, h512);
            tmp_ = _mm512_or_si512(tmp_, exit1);
            __mmask64 tmp = _mm512_cmpeq_epi8_mask(tmp_, zero512);
            if (tmp == 0x00) {
                break;
            }

            tmp = tmp & tmpb;
            l512 = _mm512_add_epi8(l512, one512);
            // NEW
            head512 = _mm512_mask_blend_epi8(tmp, head512, l512);
            tmpb = tmp;
        }

        __m512i  index512 = tail512;
        tmpb = dmask;
        l512 = _mm512_set1_epi8(end - i);   // diagonal offset of column end

        for (l = end; l >= beg; l--)
        {
            __m512i f512 = _mm512_load_si512((__m512i *)(F + l * SIMD_WIDTH8));
            __m512i h512 = _mm512_load_si512((__m512i *)(H_h + l * SIMD_WIDTH8));
            __m512i tmp_ = _mm512_or_si512(f512, h512);
            tmp_ = _mm512_or_si512(tmp_, exit1);
            __mmask64 tmp = _mm512_cmpeq_epi8_mask(tmp_, zero512);
            if (tmp == 0x00)  {
                break;
            }

            tmp = tmp & tmpb;
            l512 = _mm512_sub_epi8(l512, one512);
            // NEW
            index512 = _mm512_mask_blend_epi8(tmp, index512, l512);
            tmpb = tmp;
        }
        index512 = _mm512_add_epi8(index512, two512);
        // signed min in the offset frame against qlen-i
        tail512 = _mm512_min_epi8(index512, qlen_off512);

        // Frame shift for the next row: i advances by 1, so the same absolute
        // band edge has its diagonal offset (col - i) decremented by 1. Keep
        // head/tail tracking the same columns as the frame moves.
        head512 = _mm512_sub_epi8(head512, one512);
        tail512 = _mm512_sub_epi8(tail512, one512);

        // --- RE-BASELINE EVENT (per-lane score floor) ---
        // After this row's score state is final, lower any lane whose row-max
        // (maxRS1, rebaselined) has climbed to REBASE_HI back into range. We
        // raise B[l] by delta = maxRS1[l] - REBASE_KEEP and subtract that delta
        // (saturating to 0) from every carried score buffer/register for that
        // lane. Done AFTER band narrowing so this row's trim matches the 16-bit
        // reference; the next row reads the lowered values consistently.
        {
            /* Fast path: one vector test for "any lane >= REBASE_HI?" — see the
             * 128-bit kernel. Under the routing envelope this is ~always false, so
             * the per-lane scalar scan + maxRS1 store below are skipped each row. */
            if (_mm512_cmpge_epu8_mask(maxRS1, _mm512_set1_epi8((int8_t) REBASE_HI))) {
            uint8_t rs_b[SIMD_WIDTH8] __attribute((aligned(64)));
            _mm512_store_si512((__m512i *) rs_b, maxRS1);
            int8_t  delta_a[SIMD_WIDTH8] __attribute((aligned(64)));
            int any = 0;
            for (int l = 0; l < SIMD_WIDTH8; l++) {
                int d = 0;
                if ((int)rs_b[l] >= REBASE_HI) {
                    d = (int)rs_b[l] - REBASE_KEEP;   // keep low REBASE_KEEP, in range
                    B[l] += d;
                    any = 1;
                }
                delta_a[l] = (int8_t) d;              // 0..(REBASE_HI-REBASE_KEEP) small
            }
            if (any) {
                __m512i delta512 = _mm512_load_si512((__m512i *) delta_a);
                // Subtract delta from every carried score buffer over the active
                // band columns. The range [lo,hi] is the union of the current and
                // next band windows; clamping to [0, ncol-1] keeps the loop within
                // the valid column range (a safe superset covering every column the
                // next row will read).
                int lo = nbeg < beg ? nbeg : beg;
                int hi = nend > end ? nend : end;
                if (lo < 0) lo = 0;
                if (hi + 1 > ncol) hi = ncol - 1;
                for (int l = lo; l <= hi; l++) {
                    __m512i h512 = _mm512_load_si512((__m512i *)(H_h + l * SIMD_WIDTH8));
                    __m512i f512 = _mm512_load_si512((__m512i *)(F   + l * SIMD_WIDTH8));
                    h512 = _mm512_subs_epu8(h512, delta512);
                    f512 = _mm512_subs_epu8(f512, delta512);
                    _mm512_store_si512((__m512i *)(H_h + l * SIMD_WIDTH8), h512);
                    _mm512_store_si512((__m512i *)(F   + l * SIMD_WIDTH8), f512);
                }
                // The vertical seed H_v is read only while the band still touches
                // column 0 (beg==0). Lower the rows that may still be consumed.
                if (beg == 0) {
                    for (int r = i + 1; r < nrow; r++) {
                        __m512i v512 = _mm512_load_si512((__m512i *)(H_v + r * SIMD_WIDTH8));
                        v512 = _mm512_subs_epu8(v512, delta512);
                        _mm512_store_si512((__m512i *)(H_v + r * SIMD_WIDTH8), v512);
                    }
                }
                // Carried score maxima (registers) live in the same rebaselined
                // frame and must be lowered too so the next row's z-drop (ms-rs)
                // and the absolute-frame add-back stay consistent. maxScore512 is
                // the running max, so it never saturates to 0 (maxScore512 >=
                // maxRS1 for the triggering row, and maxRS1 >= REBASE_HI >
                // REBASE_KEEP > delta, so it stays positive).
                maxScore512 = _mm512_subs_epu8(maxScore512, delta512);
                // gscore no longer carries a rebaselined byte register: the
                // query-end cell H is captured per row and accumulated wide
                // (gbest_abs, absolute units) in the epilogue, immune to this
                // re-baseline subtract.
            }
            }   /* end "any lane >= REBASE_HI" fast-path guard */
        }

#if RDT
        prof[DP2][0] += __rdtsc() - tim1;
#endif
    }
    
#if RDT
    prof[DP][0] += __rdtsc() - tim;
#endif

    // Scores are reconstructed in ABSOLUTE units from the per-lane wide
    // best_abs/gbest_abs side channels (= rebaselined byte + B[l], captured each
    // row before any later re-baseline could saturate the byte away). Positions
    // are reconstructed wide from the diagonal-offset lanes plus the per-lane
    // best-row side channels.
    int8_t maxj[SIMD_WIDTH8]  __attribute((aligned(64)));
    _mm512_store_si512((__m512i *) maxj, y512);   // best col as diagonal offset

    int8_t max_off_ar[SIMD_WIDTH8]  __attribute((aligned(64)));
    _mm512_store_si512((__m512i *) max_off_ar, max_off512);

    for(i = 0; i < SIMD_WIDTH8; i++)
    {
        p[i].score   = best_abs[i];                               // absolute score
        p[i].tle     = xrow[i];                                   // best row (target end)
        // qle reconstruction: maxj[l] = col - (xrow-1), so col = maxj[l] + (xrow-1).
        // Guard the unset case (xrow==0 means no update occurred).
        p[i].qle     = (xrow[i] == 0) ? 0 : ((int) maxj[i] + xrow[i] - 1);
        p[i].max_off = (uint8_t) max_off_ar[i];
        p[i].gscore  = gbest_abs[i];                              // absolute gscore (-1 if unset)
        p[i].gtle    = ierow[i];                                  // best gscore row
    }

    return;
}
//----------------------------AVX512 vec 16 bit SIMD lane -------------------------------------
#define PFD16 2
void BandedPairWiseSW::getScores16(SeqPair *pairArray,
                                   uint8_t *seqBufRef,
                                   uint8_t *seqBufQer,
                                   int32_t numPairs,
                                   uint16_t numThreads,
                                   int32_t w)
{
    int i;
    int64_t startTick, endTick;

    {
        BswOvershootGuard _g(pairArray, numPairs, SIMD_WIDTH16, guard_overshoot_);
        smithWatermanBatchWrapper16(pairArray, seqBufRef, seqBufQer, numPairs, numThreads, w);
    }

#if MAXI
    printf("AVX512 Vecor code: Writing output..\n");
    for (int l=0; l<numPairs; l++)
    {
        fprintf(stderr, "%d (%d %d) %d %d %d\n",
                pairArray[l].score, pairArray[l].tle, pairArray[l].qle,
                pairArray[l].gscore, pairArray[l].max_off, pairArray[l].gtle);

    }
    printf("Vector code: Writing output completed!!!\n\n");
#endif

}

void BandedPairWiseSW::smithWatermanBatchWrapper16(SeqPair *pairArray,
                                                   uint8_t *seqBufRef,
                                                   uint8_t *seqBufQer,
                                                   int32_t numPairs,
                                                   uint16_t numThreads,
                                                   int32_t w)
{
    numThreads = effective_threads(numThreads);
    int64_t st1, st2, st3, st4, st5;
#if RDT
    st1 = ___rdtsc();
#endif

    uint16_t *seq1SoA = (uint16_t *)_mm_malloc((size_t)MAX_SEQ_LEN16 * SIMD_WIDTH16 * numThreads * sizeof(uint16_t), 64);
    uint16_t *seq2SoA = (uint16_t *)_mm_malloc((size_t)MAX_SEQ_LEN16 * SIMD_WIDTH16 * numThreads * sizeof(uint16_t), 64);

    if (UNLIKELY(seq1SoA == NULL || seq2SoA == NULL)) {
        fprintf(stderr, "Error! Mem not allocated!!!\n");
        exit(EXIT_FAILURE);
    }
        
    int32_t ii;
    int32_t roundNumPairs = ((numPairs + SIMD_WIDTH16 - 1)/SIMD_WIDTH16 ) * SIMD_WIDTH16;
    for(ii = numPairs; ii < roundNumPairs; ii++)
    {
        pairArray[ii].id = ii;
        pairArray[ii].len1 = 0;
        pairArray[ii].len2 = 0;
        pairArray[ii].idr = 0;
        pairArray[ii].idq = 0;
    }

#if RDT
    st2 = ___rdtsc();
#endif
    
    
#if RDT
    st3 = ___rdtsc();
#endif

    int eb = end_bonus;
//#pragma omp parallel num_threads(numThreads)
    {
        int32_t i;
        uint16_t tid = 0; //omp_get_thread_num();
        uint16_t *mySeq1SoA = seq1SoA + tid * MAX_SEQ_LEN16 * SIMD_WIDTH16;
        uint16_t *mySeq2SoA = seq2SoA + tid * MAX_SEQ_LEN16 * SIMD_WIDTH16;
        uint8_t *seq1;
        uint8_t *seq2;
        uint16_t h0[SIMD_WIDTH16]   __attribute__((aligned(64)));
        uint16_t band[SIMD_WIDTH16];        
        uint16_t qlen[SIMD_WIDTH16] __attribute__((aligned(64)));
        int32_t bsize = 0;
        __mmask16 dmask4 = 0xFFFF;
        
        int16_t *H1 = H16_ + tid * SIMD_WIDTH16 * MAX_SEQ_LEN16;
        int16_t *H2 = H16__ + tid * SIMD_WIDTH16 * MAX_SEQ_LEN16;

        __m512i zero512   = _mm512_setzero_si512();
        __m512i o_ins512  = _mm512_set1_epi16(o_ins);
        __m512i e_ins512  = _mm512_set1_epi16(e_ins);
        __m512i oe_ins512 = _mm512_set1_epi16(o_ins + e_ins);
        __m512i o_del512  = _mm512_set1_epi16(o_del);
        __m512i e_del512  = _mm512_set1_epi16(e_del);
        __m512i eb_ins512 = _mm512_set1_epi16(eb - o_ins);
        __m512i eb_del512 = _mm512_set1_epi16(eb - o_del);
        
        int16_t max = 0;
        if (max < w_match) max = w_match;
        if (max < w_mismatch) max = w_mismatch;
        if (max < w_ambig) max = w_ambig;
        
        int nstart = 0, nend = numPairs;
        
//#pragma omp for schedule(dynamic, 128)
        for(i = nstart; i < nend; i+=SIMD_WIDTH16)
        {
            int32_t j, k;
            uint16_t maxLen1 = 0;
            uint16_t maxLen2 = 0;
            bsize = w;

            for(j = 0; j < SIMD_WIDTH16; j++)
            {
                if ((i + j + PFD16) < roundNumPairs) { // prefetch block (bounded; see getScores8/16 contract)
                    SeqPair spf = pairArray[i + j + PFD16];
                    _mm_prefetch((const char*) seqBufRef + (int64_t)spf.idr, _MM_HINT_NTA);
                    _mm_prefetch((const char*) seqBufRef + (int64_t)spf.idr + 64, _MM_HINT_NTA);
                }
                SeqPair sp = pairArray[i + j];
                h0[j] = sp.h0;

                seq1 = seqBufRef + (int64_t)sp.idr;

                for(k = 0; k < sp.len1; k++)
                {
                    mySeq1SoA[k * SIMD_WIDTH16 + j] = (seq1[k] == AMBIG ? dmask4:seq1[k]);
                }
                
                qlen[j] = sp.len2 * max;
                if(maxLen1 < sp.len1) maxLen1 = sp.len1;
            }

            for(j = 0; j < SIMD_WIDTH16; j++)
            {
                SeqPair sp = pairArray[i + j];
                for(k = sp.len1; k <= maxLen1; k++)
                {
                    mySeq1SoA[k * SIMD_WIDTH16 + j] = DUMMY1;
                }
            }
            /* B5: only boundary row H2[maxLen1] survives the seed below. */
            _mm512_store_si512((__m512i *)(H2 + maxLen1 * SIMD_WIDTH16), _mm512_set1_epi16((short)DUMMY1));
//--------------------
            __m512i h0_512 = _mm512_load_si512((__m512i*) h0);
            _mm512_store_si512((__m512i *) H2, h0_512);
            __m512i tmp512 = _mm512_sub_epi16(h0_512, o_del512);
            
            for(k = 1; k < maxLen1; k++)
            {
                tmp512 = _mm512_sub_epi16(tmp512, e_del512);
                __m512i tmp512_ = _mm512_max_epi16(tmp512, zero512);
                _mm512_store_si512((__m512i *)(H2 + k* SIMD_WIDTH16), tmp512_);
            }
//-------------------
            for(j = 0; j < SIMD_WIDTH16; j++)
            {
                if ((i + j + PFD16) < roundNumPairs) { // prefetch block (bounded; see getScores8/16 contract)
                    SeqPair spf = pairArray[i + j + PFD16];
                    _mm_prefetch((const char*) seqBufQer + (int64_t)spf.idq, _MM_HINT_NTA);
                    _mm_prefetch((const char*) seqBufQer + (int64_t)spf.idq + 64, _MM_HINT_NTA);
                }
                
                SeqPair sp = pairArray[i + j];
                seq2 = seqBufQer + (int64_t)sp.idq;
                for(k = 0; k < sp.len2; k++)
                {
                    mySeq2SoA[k * SIMD_WIDTH16 + j] = (seq2[k]==AMBIG?dmask4:seq2[k]);
                }
                if(maxLen2 < sp.len2) maxLen2 = sp.len2;
            }
            
            for(j = 0; j < SIMD_WIDTH16; j++)
            {
                SeqPair sp = pairArray[i + j];
                for(k = sp.len2; k <= maxLen2; k++)
                {
                    mySeq2SoA[k * SIMD_WIDTH16 + j] = DUMMY2;
                }
            }
            /* B5: only boundary row H1[maxLen2]=0 survives the seed below. */
            _mm512_store_si512((__m512i *)(H1 + maxLen2 * SIMD_WIDTH16), _mm512_setzero_si512());
//------------------------
            _mm512_store_si512((__m512i *) H1, h0_512);
            __mmask32 mask512 = _mm512_cmpgt_epi16_mask(h0_512, oe_ins512);
            tmp512 = _mm512_sub_epi16(h0_512, oe_ins512);
            tmp512 = _mm512_mask_blend_epi16(mask512, zero512, tmp512);
            _mm512_store_si512((__m512i *) (H1 + SIMD_WIDTH16), tmp512);

            for(k = 2; k < maxLen2; k++)
            {
                __m512i h1_512 = tmp512;
                tmp512 = _mm512_sub_epi16(h1_512, e_ins512);
                tmp512 = _mm512_max_epi16(tmp512, zero512);
                _mm512_store_si512((__m512i *)(H1 + k*SIMD_WIDTH16), tmp512);
            }           
//------------------------

            /* Banding calculation in pre-processing */
            uint16_t myband[SIMD_WIDTH16] __attribute__((aligned(64)));
            uint16_t temp[SIMD_WIDTH16] __attribute__((aligned(64)));
            {
                __m512i qlen512 = _mm512_load_si512((__m512i *) qlen);
                __m512i sum512 = _mm512_add_epi16(qlen512, eb_ins512);
                _mm512_store_si512((__m512i *) temp, sum512);               
                for (int l=0; l<SIMD_WIDTH16; l++) {
                    double val = temp[l]/e_ins + 1.0;
                    int max_ins = val;
                    max_ins = max_ins > 1? max_ins : 1;
                    myband[l] = min_(bsize, max_ins);
                }
                sum512 = _mm512_add_epi16(qlen512, eb_del512);
                _mm512_store_si512((__m512i *) temp, sum512);               
                for (int l=0; l<SIMD_WIDTH16; l++) {
                    double val = temp[l]/e_del + 1.0;
                    int max_ins = val;
                    max_ins = max_ins > 1? max_ins : 1;
                    myband[l] = min_(myband[l], max_ins);
                    bsize = bsize < myband[l] ? myband[l] : bsize;
                }               
            }

            smithWaterman512_16(mySeq1SoA,
                                mySeq2SoA,
                                maxLen1,
                                maxLen2,
                                pairArray + i,
                                h0,
                                tid,
                                numPairs,
                                zdrop,
                                bsize,
                                qlen,
                                myband);
        }
    }

#if RDT 
    st4 = ___rdtsc();
#endif
    

#if RDT
    st5 = ___rdtsc();
    setupTicks += st2 - st1;
    sort1Ticks += st3 - st2;
    swTicks += st4 - st3;
    sort2Ticks += st5 - st4;
#endif

    // free mem
    _mm_free(seq1SoA);
    _mm_free(seq2SoA);

    return;
}

void BandedPairWiseSW::smithWaterman512_16(uint16_t seq1SoA[],
                                           uint16_t seq2SoA[],
                                           uint16_t nrow,
                                           uint16_t ncol,
                                           SeqPair *p,
                                           uint16_t *h0,
                                           uint16_t tid,
                                           int32_t numPairs,
                                           int zdrop,
                                           int32_t w,
                                           uint16_t qlen[],
                                           uint16_t myband[])
{   
    __m512i match512     = _mm512_set1_epi16(this->w_match);
    __m512i mismatch512  = _mm512_set1_epi16(this->w_mismatch);
    __m512i gapOpen512   = _mm512_set1_epi16(this->w_open);
    __m512i gapExtend512 = _mm512_set1_epi16(this->w_extend);
    __m512i gapOE512     = _mm512_set1_epi16(this->w_open + this->w_extend);
    __m512i w_ambig_512  = _mm512_set1_epi16(this->w_ambig);    // ambig penalty
    __m512i five512      = _mm512_set1_epi16(5);

    // D3 generic-matrix seam: symmetric default uses SYM; a single freed-to-match
    // cell (bisulfite) uses rank-1; any other asymmetric matrix uses the amat LUT.
    const bool forced  = bsw_force_generic_matrix();
    const bool gen_mat = bsw_generic_matrix(this->mat, this->w_match, this->w_mismatch)
                         || forced;
    const BswFreedCell fc = bsw_freed_cell(this->mat, this->w_match, this->w_mismatch, forced);
    __m512i frref512  = _mm512_set1_epi16(fc.ref);
    __m512i frread512 = _mm512_set1_epi16(fc.read);
    int8_t amat_bytes[16] __attribute__((aligned(16)));
    build_amat16(amat_bytes, this->mat);
    __m512i amat512   = _mm512_broadcast_i32x4(_mm_load_si128((__m128i *)amat_bytes));
    __m512i three512  = _mm512_set1_epi16(3);

    __m512i e_del512    = _mm512_set1_epi16(this->e_del);
    __m512i oe_del512   = _mm512_set1_epi16(this->o_del + this->e_del);
    __m512i e_ins512    = _mm512_set1_epi16(this->e_ins);
    __m512i oe_ins512   = _mm512_set1_epi16(this->o_ins + this->e_ins);

    int16_t *F   = F16_ + tid * SIMD_WIDTH16 * MAX_SEQ_LEN16;
    int16_t *H_h = H16_ + tid * SIMD_WIDTH16 * MAX_SEQ_LEN16;
    int16_t *H_v = H16__ + tid * SIMD_WIDTH16 * MAX_SEQ_LEN16;

    int lane = 0;

    int16_t i, j;

    uint16_t tlen[SIMD_WIDTH16];
    uint16_t tail[SIMD_WIDTH16] __attribute((aligned(64)));
    uint16_t head[SIMD_WIDTH16] __attribute((aligned(64)));

    // PR 17: per-row score-vector scratch for fission (AVX-512 16-bit).
    int16_t *sbt_buf = sbt16_ + tid * SIMD_WIDTH16 * MAX_SEQ_LEN16;
    
    int32_t minq = 10000000;
    for (int l=0; l<SIMD_WIDTH16; l++) {
        tlen[l] = p[l].len1;
        qlen[l] = p[l].len2;
        if (p[l].len2 < minq) minq = p[l].len2;
    }
    minq -= 1; // for gscore

    __m512i tlen512   = _mm512_load_si512((__m512i *) tlen);
    __m512i qlen512   = _mm512_load_si512((__m512i *) qlen);
    __m512i myband512 = _mm512_load_si512((__m512i *) myband);
    __m512i zero512   = _mm512_setzero_si512();
    __m512i one512    = _mm512_set1_epi16(1);
    __m512i two512    = _mm512_set1_epi16(2);
    __m512i i512_1    = _mm512_set1_epi16(1);
    __m512i max_ie512 = zero512;
    __mmask16 dmask4 = 0xFFFF;
    __m512i ff512     = _mm512_set1_epi16(dmask4);
    
    __m512i tail512 = qlen512, head512 = zero512;
    _mm512_store_si512((__m512i *) head, head512);
    _mm512_store_si512((__m512i *) tail, tail512);

    __m512i mlen512 = _mm512_add_epi16(qlen512, myband512);
    mlen512 = _mm512_min_epu16(mlen512, tlen512);
    
    uint16_t temp[SIMD_WIDTH16]  __attribute((aligned(64)));
    uint16_t temp1[SIMD_WIDTH16]  __attribute((aligned(64)));
    
    __m512i s00  = _mm512_load_si512((__m512i *)(seq1SoA));
    __m512i hval = _mm512_load_si512((__m512i *)(H_v));
    __mmask32 dmask = 0xFFFFFFFF;
    
///////
    __m512i maxScore512 = hval;
    for(j = 0; j < ncol; j++)
        _mm512_store_si512((__m512i *)(F + j * SIMD_WIDTH16), zero512);
    
    __m512i x512       = zero512;
    __m512i y512       = zero512;
    __m512i i512       = zero512;
    __m512i gscore     = _mm512_set1_epi16(-1);
    __m512i max_off512 = zero512;
    __m512i exit0      = _mm512_set1_epi16(dmask4);
    __m512i zdrop512   = _mm512_set1_epi16(zdrop);

    int beg = 0, end = ncol;
    int nbeg = beg, nend = end;

#if RDT
    uint64_t tim = __rdtsc();
#endif

    for(i = 0; i < nrow; i++)
    {
        __m512i e11 = zero512;
        __m512i h00, h11, h10;
        __m512i s10 = _mm512_load_si512((__m512i *)(seq1SoA + (i + 0) * SIMD_WIDTH16));

        beg = nbeg; end = nend;
        int pbeg = beg;
        if (beg < i - w) beg = i - w;
        if (end > i + w + 1) end = i + w + 1;
        if (end > ncol) end = ncol;

        h10 = zero512;
        if (beg == 0)
            h10 = _mm512_load_si512((__m512i *)(H_v + (i+1) * SIMD_WIDTH16));

        __m512i j512 = zero512;
        __m512i maxRS1, maxRS2, maxRS3, maxRS4;
        maxRS1 = zero512;
        
        __m512i i1_512 = _mm512_set1_epi16(i+1);
        __m512i y1_512 = zero512;
        
#if RDT 
        uint64_t tim1 = __rdtsc();
#endif
        
        /* Banding */
        __m512i i512, cache512, max512;
        __m512i phead512 = head512, ptail512 = tail512;
        i512 = _mm512_set1_epi16(i);
        cache512 = _mm512_sub_epi16(i512, myband512);
        head512  = _mm512_max_epi16(head512, cache512);
        cache512 = _mm512_add_epi16(i1_512, myband512);
        tail512  = _mm512_min_epu16(tail512, cache512);
        tail512  = _mm512_min_epu16(tail512, qlen512);
        /* Banding ends */
        
        // NEW, trimming.
        __mmask32 cmph = _mm512_cmpeq_epi16_mask(head512, phead512);
        __mmask32 cmpt = _mm512_cmpeq_epi16_mask(tail512, ptail512);
        cmph &= cmpt;
        for (int l=beg; l<end && cmph != dmask; l++)
        {
            __m512i h512 = _mm512_load_si512((__m512i *)(H_h + l * SIMD_WIDTH16));
            __m512i f512 = _mm512_load_si512((__m512i *)(F + l * SIMD_WIDTH16));

            __m512i pj512 = _mm512_set1_epi16(l);
            __m512i j512 = _mm512_set1_epi16(l+1);
            __mmask32 cmp1 = _mm512_cmpgt_epi16_mask(head512, pj512);
            if (cmp1 == 0x00) break;
            // __mmask32 cmp2 = _mm512_cmpgt_epi16_mask(pj512, tail512);
            __mmask32 cmp2 = _mm512_cmpgt_epi16_mask(j512, tail512);
            cmp1 = cmp1 | cmp2;
            h512 = _mm512_mask_blend_epi16(cmp1, h512, zero512);
            f512 = _mm512_mask_blend_epi16(cmp1, f512, zero512);

            _mm512_store_si512((__m512i *)(F + l * SIMD_WIDTH16), f512);
            _mm512_store_si512((__m512i *)(H_h + l * SIMD_WIDTH16), h512);
        }

#if RDT
        prof[DP3][0] += __rdtsc() - tim1;
#endif

        // beg = nbeg; end = nend;
        __mmask32 cmp512_1 = _mm512_cmpgt_epi16_mask(i1_512, tlen512);

        /* Updating row exit status */
        __mmask32 cmpim = _mm512_cmpgt_epi16_mask(i1_512, mlen512);
        __mmask32 cmpht = _mm512_cmpeq_epi16_mask(tail512, head512);
        cmpim = cmpim | cmpht;
        // NEW
        cmpht = _mm512_cmpgt_epi16_mask(head512, tail512);
        cmpim = cmpim |  cmpht;

        exit0 = _mm512_mask_blend_epi16(cmpim, exit0, zero512);
        
#if RDT
        tim1 = __rdtsc();
#endif
        
        // PR 17: AVX-512 16-bit score pre-pass (fission). gen_mat branch is
        // loop-invariant: symmetric keeps the cheap SYM prepass.
        if (!gen_mat) {
            for (int jp = beg; jp < end; jp++) {
                __m512i s2 = _mm512_load_si512((__m512i *)(seq2SoA + jp * SIMD_WIDTH16));
                __m512i sbt11;
                SBT_PREPASS16_SYM(s10, s2, sbt11, mismatch512, match512, w_ambig_512);
                _mm512_store_si512((__m512i *)(sbt_buf + jp * SIMD_WIDTH16), sbt11);
            }
        } else if (fc.rank1) {
            __m512i alt10 = _mm512_mask_blend_epi16(_mm512_cmpeq_epi16_mask(s10, frref512), s10, frread512);
            for (int jp = beg; jp < end; jp++) {
                __m512i s2 = _mm512_load_si512((__m512i *)(seq2SoA + jp * SIMD_WIDTH16));
                __m512i sbt11;
                SBT_PREPASS16_RANK1(s10, s2, alt10, sbt11, mismatch512, match512, w_ambig_512);
                _mm512_store_si512((__m512i *)(sbt_buf + jp * SIMD_WIDTH16), sbt11);
            }
        } else {
            for (int jp = beg; jp < end; jp++) {
                __m512i s2 = _mm512_load_si512((__m512i *)(seq2SoA + jp * SIMD_WIDTH16));
                __m512i sbt11;
                SBT_PREPASS16_AMAT(s10, s2, sbt11, amat512, mismatch512, w_ambig_512, three512);
                _mm512_store_si512((__m512i *)(sbt_buf + jp * SIMD_WIDTH16), sbt11);
            }
        }

        j512 = _mm512_set1_epi16(beg);
        for(j = beg; j < end; j++)
        {
            __m512i f11, f21, f31, f41, f51, jj512, sbt11;
            h00 = _mm512_load_si512((__m512i *)(H_h + j * SIMD_WIDTH16));
            f11 = _mm512_load_si512((__m512i *)(F + j * SIMD_WIDTH16));
            sbt11 = _mm512_load_si512((__m512i *)(sbt_buf + j * SIMD_WIDTH16));

            __m512i pj512 = j512;
            j512 = _mm512_add_epi16(j512, one512);

            MAIN_CODE16_CORE(sbt11, h00, h11, e11, f11, f21, zero512,
                             e_ins512, oe_ins512,
                             e_del512, oe_del512);

            // Masked writing
            __mmask32 cmp2 = _mm512_cmpgt_epi16_mask(head512, pj512);
            __mmask32 cmp1 = _mm512_cmpgt_epi16_mask(pj512, tail512);
            cmp1 = cmp1 | cmp2;
            h10 = _mm512_mask_blend_epi16(cmp1, h10, zero512);
            f21 = _mm512_mask_blend_epi16(cmp1, f21, zero512);
            
            /* Part of main code MAIN_CODE */
            __m512i bmaxRS = maxRS1, blend512;                                      
            maxRS1 =_mm512_max_epi16(maxRS1, h11);                          
            __mmask32 cmpA = _mm512_cmpgt_epi16_mask(maxRS1, bmaxRS);                   
            __mmask32 cmpB =_mm512_cmpeq_epi16_mask(maxRS1, h11);                   
            cmpA = cmpA | cmpB;
            cmp1 = _mm512_cmpgt_epi16_mask(j512, tail512);
            cmp1 = cmp1 | cmp2;         
            blend512 = _mm512_mask_blend_epi16(cmpA, y1_512, j512);
            y1_512 = _mm512_mask_blend_epi16(cmp1, blend512, y1_512);
            maxRS1 = _mm512_mask_blend_epi16(cmp1, maxRS1, bmaxRS);                     

            _mm512_store_si512((__m512i *)(F + j * SIMD_WIDTH16), f21);
            _mm512_store_si512((__m512i *)(H_h + j * SIMD_WIDTH16), h10);

            h10 = h11;
                        
            /* gscore calculations */
            if (j >= minq)
            {
                __mmask32 cmp = _mm512_cmpeq_epi16_mask(j512, qlen512);
                __m512i max_gh = _mm512_max_epi16(gscore, h11);
                __mmask32 cmp_gh = _mm512_cmpgt_epi16_mask(gscore, h11);
                __m512i tmp512_1 = _mm512_mask_blend_epi16(cmp_gh, i1_512, max_ie512);

                tmp512_1 = _mm512_mask_blend_epi16(cmp, max_ie512, tmp512_1);
                __mmask32 mex0 = _mm512_movepi16_mask(exit0);
                tmp512_1 = _mm512_mask_blend_epi16(mex0, max_ie512, tmp512_1);
                
                max_gh = _mm512_mask_blend_epi16(mex0, gscore, max_gh);
                max_gh = _mm512_mask_blend_epi16(cmp, gscore, max_gh);              

                cmp = _mm512_cmpgt_epi16_mask(j512, tail512); 
                max_gh = _mm512_mask_blend_epi16(cmp, max_gh, gscore);
                max_ie512 = _mm512_mask_blend_epi16(cmp, tmp512_1, max_ie512);
                gscore = max_gh;
            }
        }        
        __mmask32 cmp2 = _mm512_cmpgt_epi16_mask(head512, j512);
        __mmask32 cmp1 = _mm512_cmpgt_epi16_mask(j512, tail512);
        cmp1 = cmp1 | cmp2;
        h10 = _mm512_mask_blend_epi16(cmp1, h10, zero512);
        
        _mm512_store_si512((__m512i *)(H_h + j * SIMD_WIDTH16), h10);
        _mm512_store_si512((__m512i *)(F + j * SIMD_WIDTH16), zero512);
                        
        /* exit due to zero score by a row */
        __mmask32 cval = dmask;
        __m512i bmaxScore512 = maxScore512;
        __mmask32 tmp = _mm512_cmpeq_epi16_mask(maxRS1, zero512);
        if (cval == tmp) break;

        exit0 = _mm512_mask_blend_epi16(tmp, exit0, zero512);

        __m512i score512 = _mm512_max_epi16(maxScore512, maxRS1);
        __mmask32 mex0 = _mm512_movepi16_mask(exit0);
        maxScore512 = _mm512_mask_blend_epi16(mex0, maxScore512, score512);

        __mmask32 cmp = _mm512_cmpgt_epi16_mask(maxScore512, bmaxScore512);
        y512 = _mm512_mask_blend_epi16(cmp, y512, y1_512);
        x512 = _mm512_mask_blend_epi16(cmp, x512, i1_512);
        
        /* max_off calculations */
        __m512i ind512 = _mm512_sub_epi16(y1_512, i1_512);
        ind512 = _mm512_abs_epi16(ind512);
        __m512i bmax_off512 = max_off512;
        ind512 = _mm512_max_epi16(max_off512, ind512);
        max_off512 = _mm512_mask_blend_epi16(cmp, bmax_off512, ind512);

        /* Z-score condition for exit */
        ZSCORE16(i1_512, y1_512);       

#if RDT
        prof[DP1][0] += __rdtsc() - tim1;
#endif
        
        /* Narrowing of the band */
        /* Part 1: From beg */
        cval = dmask;
        int l;      
        for (l = beg; l < end; l++)
        {
            __m512i f512 = _mm512_load_si512((__m512i *)(F + l * SIMD_WIDTH16));
            __m512i h512 = _mm512_load_si512((__m512i *)(H_h + l * SIMD_WIDTH16));
            __m512i tmp = _mm512_or_si512(f512, h512);
            __mmask32 val = _mm512_cmpeq_epi16_mask(tmp, zero512);
            if (cval == val) nbeg = l;
            else
                break;
        }
        
        /* From end */
        bool flg = 1;
        for (l = end; l >= beg; l--) {
            __m512i f512 = _mm512_load_si512((__m512i *)(F + l * SIMD_WIDTH16));
            __m512i h512 = _mm512_load_si512((__m512i *)(H_h + l * SIMD_WIDTH16));
            __m512i tmp = _mm512_or_si512(f512, h512);
            __mmask32 val = _mm512_cmpeq_epi16_mask(tmp, zero512);
            if (val != cval && flg)  
                break;
        }
        nend = l + 2 < ncol? l + 2: ncol;

#if RDT
        tim1 = __rdtsc();
#endif
        /* Setting of head and tail for each pair */
        // beg = nbeg; end = l; // keep check on this!!
        beg = nbeg; end = nend; 
        
        __m512i tail512_ = _mm512_sub_epi16(tail512, one512);
        __m512i exit1 = _mm512_xor_si512(exit0, ff512);
        __mmask32 tmpb = dmask;
        __m512i l512 = _mm512_set1_epi16(beg);
        
        for (l = beg; l < end; l++)
        {
            __m512i f512 = _mm512_load_si512((__m512i *)(F + l * SIMD_WIDTH16));
            __m512i h512 = _mm512_load_si512((__m512i *)(H_h + l * SIMD_WIDTH16));  
            __m512i tmp_ = _mm512_or_si512(f512, h512);
            tmp_ = _mm512_or_si512(tmp_, exit1);            
            __mmask32 tmp = _mm512_cmpeq_epi16_mask(tmp_, zero512);
            if (tmp == 0x00) {
                break;
            }
            
            tmp = tmp & tmpb;
            l512 = _mm512_add_epi16(l512, one512);
            // NEW
            head512 = _mm512_mask_blend_epi16(tmp, head512, l512);

            tmpb = tmp;         
        }
        
        __m512i  index512 = tail512;
        tmpb = dmask;
        l512 = _mm512_set1_epi16(end);
        
        for (l = end; l >= beg; l--)
        {
            __m512i f512 = _mm512_load_si512((__m512i *)(F + l * SIMD_WIDTH16));
            __m512i h512 = _mm512_load_si512((__m512i *)(H_h + l * SIMD_WIDTH16));          
            __m512i tmp_ = _mm512_or_si512(f512, h512);
            tmp_ = _mm512_or_si512(tmp_, exit1);
            __mmask32 tmp = _mm512_cmpeq_epi16_mask(tmp_, zero512);         
            if (tmp == 0x00)  {
                break;
            }

            tmp = tmp & tmpb;
            l512 = _mm512_sub_epi16(l512, one512);
            // NEW
            index512 = _mm512_mask_blend_epi16(tmp, index512, l512);

            tmpb = tmp;
        }
        index512 = _mm512_add_epi16(index512, two512);
        tail512 = _mm512_min_epi16(index512, qlen512);

#if RDT
        prof[DP2][0] += __rdtsc() - tim1;
#endif
    }
    
#if RDT
    prof[DP][0] += __rdtsc() - tim;
#endif
    
    int16_t score[SIMD_WIDTH16]  __attribute((aligned(64)));
    _mm512_store_si512((__m512i *) score, maxScore512);

    int16_t maxi[SIMD_WIDTH16]  __attribute((aligned(64)));
    _mm512_store_si512((__m512i *) maxi, x512);

    int16_t maxj[SIMD_WIDTH16]  __attribute((aligned(64)));
    _mm512_store_si512((__m512i *) maxj, y512);

    int16_t max_off_ar[SIMD_WIDTH16]  __attribute((aligned(64)));
    _mm512_store_si512((__m512i *) max_off_ar, max_off512);

    int16_t gscore_ar[SIMD_WIDTH16]  __attribute((aligned(64)));
    _mm512_store_si512((__m512i *) gscore_ar, gscore);

    int16_t maxie_ar[SIMD_WIDTH16]  __attribute((aligned(64)));
    _mm512_store_si512((__m512i *) maxie_ar, max_ie512);
    
    for(i = 0; i < SIMD_WIDTH16; i++)
    {
        p[i].score = score[i];
        p[i].tle = maxi[i];
        p[i].qle = maxj[i];
        p[i].max_off = max_off_ar[i];
        p[i].gscore = gscore_ar[i];
        p[i].gtle = maxie_ar[i];
    }   

    return;
}
#endif  //avx512


/**************** SSE2 code ******************/
/* SBT_PREPASS8_LUT below uses _mm_shuffle_epi8 (PSHUFB, an SSSE3 intrinsic),
 * so this block requires SSSE3. The Makefile passes -mssse3 on all x86 arch
 * targets and sse2neon predefines __SSSE3__=1, so this matches every build
 * that compiles the SSE2/NEON path. */
#if ((!__AVX512BW__) && (!__AVX2__) && (__SSE2__) && (__SSSE3__))

// SSE2/NEON - 16 bit blendv
static inline __m128i
_mm_blendv_epi16(__m128i x, __m128i y, __m128i mask)
{
#if defined(__ARM_NEON) || defined(__aarch64__) || defined(APPLE_SILICON)
    // Use NEON vbsl (bitwise select) - more efficient than AND/OR/ANDNOT
    return vreinterpretq_m128i_s16(
        vbslq_s16(vreinterpretq_u16_m128i(mask),
                  vreinterpretq_s16_m128i(y),
                  vreinterpretq_s16_m128i(x)));
#else
    // x86 SSE2: Replace bit in x with bit in y when matching bit in mask is set
    return _mm_or_si128(_mm_andnot_si128(mask, x), _mm_and_si128(mask, y));
#endif
}

#define ZSCORE16(i4_128, y4_128)                                            \
    {                                                                   \
        __m128i tmpi = _mm_sub_epi16(i4_128, x128);                     \
        __m128i tmpj = _mm_sub_epi16(y4_128, y128);                     \
        cmp = _mm_cmpgt_epi16(tmpi, tmpj);                              \
        score128 = _mm_sub_epi16(maxScore128, maxRS1);                  \
        __m128i insdel = _mm_blendv_epi16(e_ins128, e_del128, cmp);     \
        __m128i sub_a128 = _mm_sub_epi16(tmpi, tmpj);                   \
        __m128i sub_b128 = _mm_sub_epi16(tmpj, tmpi);                   \
        tmp = _mm_blendv_epi16(sub_b128, sub_a128, cmp);                \
        tmp = _mm_sub_epi16(score128, tmp);                             \
        cmp = _mm_cmpgt_epi16(tmp, zdrop128);                           \
        exit0 = _mm_blendv_epi16(exit0, zero128, cmp);                  \
    }



// PR 17: 16-bit fission primitives, mirror of the 8-bit pair above.
// Symmetric (default) fast path: cmpeq match/mismatch + N(0xFFFF)->ambig.
#define SBT_PREPASS16_SYM(s1, s2, sbt11_out, mismatch128, match128, w_ambig_128, ff128) \
    {                                                                   \
        __m128i cmp11_ = _mm_cmpeq_epi16(s1, s2);                       \
        __m128i sbt_ = _mm_blendv_epi16(mismatch128, match128, cmp11_); \
        __m128i tmp_ = _mm_max_epu16(s1, s2);                           \
        tmp_ = _mm_cmpeq_epi16(tmp_, ff128);                            \
        sbt11_out = _mm_blendv_epi16(sbt_, w_ambig_128, tmp_);          \
    }

// D3 rank-1 fast path: symmetric matrix plus a single off-diagonal cell freed
// to a match (bisulfite OT/OB). The freed transition just extends the match
// condition — match when (ref==read) OR (ref==fr_ref AND read==fr_read) — so no
// byte-shuffle LUT and no sign-extend are needed. N(0xFFFF)->ambig exactly as
// the symmetric path; padding (DUMMY) lanes are never equal and never hit the
// freed cell, so they stay mismatch. For the no-op cell (0,0) this is identical
// to SBT_PREPASS16_SYM. The freed-cell test factors into a per-ROW invariant:
// the row base s1 is constant across the band loop, so we hoist
//   alt1 = (s1 == fr_ref) ? fr_read : s1     (computed once per row)
// and the per-cell match condition is just (read==ref) OR (read==alt1) — i.e.
// on a ref-C row, read-T also matches (the freed conversion); on any other row
// alt1==s1 so it degenerates to plain equality. Saves 2 ops/cell vs comparing
// fr_ref/fr_read per cell. ~7 ALU ops/cell, near parity with SYM.
#define SBT_PREPASS16_RANK1(s1, s2, alt1, sbt11_out, mismatch128, match128, w_ambig_128, ff128) \
    {                                                                   \
        __m128i eq_  = _mm_cmpeq_epi16(s2, s1);                         \
        __m128i fr_  = _mm_cmpeq_epi16(s2, alt1);                       \
        __m128i ism_ = _mm_or_si128(eq_, fr_);                          \
        __m128i sbt_ = _mm_blendv_epi16(mismatch128, match128, ism_);   \
        __m128i nN_  = _mm_cmpeq_epi16(_mm_max_epu16(s1, s2), ff128);   \
        sbt11_out = _mm_blendv_epi16(sbt_, w_ambig_128, nN_);          \
    }

// D3 generic-matrix seam (gated on an asymmetric matrix; the default symmetric
// path uses SBT_PREPASS16_SYM above and is untouched). For both-ACGT lanes the
// score comes from a target-major LUT amat[(ref<<2)|read], so an asymmetric
// off-diagonal (e.g. ref-C x read-T freed on bisulfite OT) is scored correctly.
// Every NON-ACGT lane is provably mismatch-or-ambig — padding pairs (DUMMY1=99
// target / DUMMY2=100 query) are never equal and real-vs-padding never matches,
// so the only non-ACGT values are w_ambig (when either base is N=0xFFFF) and
// w_mismatch (otherwise). That lets us drop the symmetric cmpeq match/mismatch
// blend entirely and select acgt ? lut : (N ? ambig : mismatch). For a symmetric
// matrix the LUT diagonal is w_match and off-diagonals are w_mismatch, so this
// is byte-identical to SBT_PREPASS16_SYM. The byte LUT lands the score in each
// 16-bit lane's low byte, so we sign-extend (slli 8, srai 8) before use.
#define SBT_PREPASS16_AMAT(s1, s2, sbt11_out, amat128, mismatch128, w_ambig_128, ff128, three128) \
    {                                                                   \
        __m128i maxb_ = _mm_max_epu16(s1, s2);                          \
        __m128i nN_   = _mm_cmpeq_epi16(maxb_, ff128);  /* either base N */ \
        __m128i base_ = _mm_blendv_epi16(mismatch128, w_ambig_128, nN_); \
        __m128i acgt_ = _mm_cmpeq_epi16(_mm_max_epu16(maxb_, three128), three128); /* both <= 3 */ \
        __m128i idx_  = _mm_or_si128(_mm_slli_epi16(s1, 2), s2);  /* (ref<<2)|read */ \
        __m128i lut_  = _mm_shuffle_epi8(amat128, idx_);  /* score in low byte */ \
        lut_ = _mm_srai_epi16(_mm_slli_epi16(lut_, 8), 8);  /* sign-extend low byte */ \
        sbt11_out = _mm_blendv_epi16(base_, lut_, acgt_);            \
    }

#define MAIN_CODE16_CORE(sbt11, h00, h11, e11, f11, f21, zero128, e_ins128, oe_ins128, e_del128, oe_del128) \
    {                                                                   \
        __m128i m11 = _mm_add_epi16(h00, sbt11);                        \
        __m128i cmp11 = _mm_cmpeq_epi16(h00, zero128);                  \
        m11 = _mm_blendv_epi16(m11, zero128, cmp11);                    \
        h11 = _mm_max_epi16(m11, e11);                                  \
        h11 = _mm_max_epi16(h11, f11);                                  \
        /* Reassociate gap recurrences off h11 (variant B, signed; keep   \
         * the existing max(.,0) floor). me reads OLD e11. */              \
        __m128i mf128 = _mm_max_epi16(m11, f11);                        \
        __m128i me128 = _mm_max_epi16(m11, e11);                        \
        __m128i temp128 = _mm_sub_epi16(mf128, oe_ins128);             \
        __m128i val128  = _mm_max_epi16(temp128, zero128);              \
        e11 = _mm_sub_epi16(e11, e_ins128);                             \
        e11 = _mm_max_epi16(val128, e11);                               \
        temp128 = _mm_sub_epi16(me128, oe_del128);                     \
        val128  = _mm_max_epi16(temp128, zero128);                      \
        f21 = _mm_sub_epi16(f11, e_del128);                             \
        f21 = _mm_max_epi16(val128, f21);                               \
    }



inline void sortPairsId(SeqPair *pairArray, int32_t first,
                        int32_t count, SeqPair *tempArray)
{

    int32_t i;
    
    for(i = 0; i < count; i++)
    {
        SeqPair sp = pairArray[i];
        int32_t pos = sp.id - first;
        tempArray[pos] = sp;
    }

    for(i = 0; i < count; i++)
        pairArray[i] = tempArray[i];
}



// SSE2
#define PFD 2
void BandedPairWiseSW::getScores16(SeqPair *pairArray,
                                   uint8_t *seqBufRef,
                                   uint8_t *seqBufQer,
                                   int32_t numPairs,
                                   uint16_t numThreads,
                                   int32_t w)
{
    {
        BswOvershootGuard _g(pairArray, numPairs, SIMD_WIDTH16, guard_overshoot_);
        smithWatermanBatchWrapper16(pairArray, seqBufRef,
                                    seqBufQer, numPairs,
                                    numThreads, w);
    }

#if MAXI
    for (int l=0; l<numPairs; l++)
    {
        fprintf(stderr, "%d (%d %d) %d %d %d\n",
                pairArray[l].score, pairArray[l].x, pairArray[l].y,
                pairArray[l].gscore, pairArray[l].max_off, pairArray[l].max_ie);

    }
#endif

}

void BandedPairWiseSW::smithWatermanBatchWrapper16(SeqPair *pairArray,
                                                   uint8_t *seqBufRef,
                                                   uint8_t *seqBufQer,
                                                   int32_t numPairs,
                                                   uint16_t numThreads,
                                                   int32_t w)
{
    numThreads = effective_threads(numThreads);
#if RDT
    int64_t st1, st2, st3, st4, st5;
    st1 = ___rdtsc();
#endif

    uint16_t *seq1SoA = (uint16_t *)_mm_malloc((size_t)MAX_SEQ_LEN16 * SIMD_WIDTH16 * numThreads * sizeof(uint16_t), 64);
    uint16_t *seq2SoA = (uint16_t *)_mm_malloc((size_t)MAX_SEQ_LEN16 * SIMD_WIDTH16 * numThreads * sizeof(uint16_t), 64);
    if (UNLIKELY(seq1SoA == NULL || seq2SoA == NULL)) {
        fprintf(stderr, "Error! Mem not allocated!!!\n");
        exit(EXIT_FAILURE);
    }

    int32_t ii;
    int32_t roundNumPairs = ((numPairs + SIMD_WIDTH16 - 1)/SIMD_WIDTH16 ) * SIMD_WIDTH16;
    // assert(roundNumPairs < BATCH_SIZE * SEEDS_PER_READ);
    for(ii = numPairs; ii < roundNumPairs; ii++)
    {
        pairArray[ii].id = ii;
        pairArray[ii].len1 = 0;
        pairArray[ii].len2 = 0;
        pairArray[ii].idr = 0;
        pairArray[ii].idq = 0;
    }

#if RDT
    st2 = ___rdtsc();
#endif
    

#if RDT
    st3 = ___rdtsc();
#endif
    
    int eb = end_bonus;
// #pragma omp parallel num_threads(numThreads)
    {
        int32_t i;
        uint16_t tid = 0; 
        uint16_t *mySeq1SoA = seq1SoA + tid * MAX_SEQ_LEN16 * SIMD_WIDTH16;
        uint16_t *mySeq2SoA = seq2SoA + tid * MAX_SEQ_LEN16 * SIMD_WIDTH16;
        assert(mySeq1SoA != NULL && mySeq2SoA != NULL);
        
        uint8_t *seq1;
        uint8_t *seq2;
        uint16_t h0[SIMD_WIDTH16]  __attribute__((aligned(64)));
        uint16_t qlen[SIMD_WIDTH16] __attribute__((aligned(64)));
        int32_t bsize = 0;
        
        int16_t *H1 = H16_ + tid * SIMD_WIDTH16 * MAX_SEQ_LEN16;
        int16_t *H2 = H16__ + tid * SIMD_WIDTH16 * MAX_SEQ_LEN16;

        __m128i zero128   = _mm_setzero_si128();
        __m128i e_ins128  = _mm_set1_epi16(e_ins);
        __m128i oe_ins128 = _mm_set1_epi16(o_ins + e_ins);
        __m128i o_del128  = _mm_set1_epi16(o_del);
        __m128i e_del128  = _mm_set1_epi16(e_del);
        __m128i eb_ins128 = _mm_set1_epi16(eb - o_ins);
        __m128i eb_del128 = _mm_set1_epi16(eb - o_del);
        
        int16_t max = 0;
        if (max < w_match) max = w_match;
        if (max < w_mismatch) max = w_mismatch;
        if (max < w_ambig) max = w_ambig;
        
        int nstart = 0, nend = numPairs;

// #pragma omp for schedule(dynamic, 128)
        for(i = nstart; i < nend; i+=SIMD_WIDTH16)
        {
            int32_t j, k;
            uint16_t maxLen1 = 0;
            uint16_t maxLen2 = 0;
            bsize = w;

            for(j = 0; j < SIMD_WIDTH16; j++)
            {
                if ((i + j + PFD) < roundNumPairs) { // prefetch block (bounded; see getScores8/16 contract)
                    SeqPair spf = pairArray[i + j + PFD];
                    _mm_prefetch((const char*) seqBufRef + (int64_t)spf.idr, _MM_HINT_NTA);
                    _mm_prefetch((const char*) seqBufRef + (int64_t)spf.idr + 64, _MM_HINT_NTA);
                }
                SeqPair sp = pairArray[i + j];
                h0[j] = sp.h0;
                seq1 = seqBufRef + (int64_t)sp.idr;
                
                for(k = 0; k < sp.len1; k++)
                {
                    mySeq1SoA[k * SIMD_WIDTH16 + j] = (seq1[k] == AMBIG?0xFFFF:seq1[k]);
                }
                qlen[j] = sp.len2 * max;
                if(maxLen1 < sp.len1) maxLen1 = sp.len1;
            }

            for(j = 0; j < SIMD_WIDTH16; j++)
            {
                SeqPair sp = pairArray[i + j];
                for(k = sp.len1; k <= maxLen1; k++) //removed "="
                {
                    mySeq1SoA[k * SIMD_WIDTH16 + j] = DUMMY1;
                }
            }
            /* B5: only the boundary row H2[maxLen1] survives the h0-prefix
             * deletion seed below (which overwrites rows [0, maxLen1)); write
             * just that row here instead of the dead per-row fills removed
             * above, before the seed to preserve store ordering. */
            _mm_store_si128((__m128i *)(H2 + maxLen1 * SIMD_WIDTH16), _mm_set1_epi16((short)DUMMY1));
//--------------------
            __m128i h0_128 = _mm_load_si128((__m128i*) h0);
            _mm_store_si128((__m128i *) H2, h0_128);
            __m128i tmp128 = _mm_sub_epi16(h0_128, o_del128);
            
            for(k = 1; k < maxLen1; k++)
            {
                tmp128 = _mm_sub_epi16(tmp128, e_del128);
                __m128i tmp128_ = _mm_max_epi16(tmp128, zero128);
                _mm_store_si128((__m128i *)(H2 + k* SIMD_WIDTH16), tmp128_);
            }
//-------------------
            for(j = 0; j < SIMD_WIDTH16; j++)
            {
                if ((i + j + PFD) < roundNumPairs) { // prefetch block (bounded; see getScores8/16 contract)
                    SeqPair spf = pairArray[i + j + PFD];
                    _mm_prefetch((const char*) seqBufQer + (int64_t)spf.idq, _MM_HINT_NTA);
                    _mm_prefetch((const char*) seqBufQer + (int64_t)spf.idq + 64, _MM_HINT_NTA);
                }
                
                SeqPair sp = pairArray[i + j];
                seq2 = seqBufQer + (int64_t)sp.idq;             
                for(k = 0; k < sp.len2; k++)
                {
                    mySeq2SoA[k * SIMD_WIDTH16 + j] = (seq2[k]==AMBIG?0xFFFF:seq2[k]);
                }
                if(maxLen2 < sp.len2) maxLen2 = sp.len2;
            }

            for(j = 0; j < SIMD_WIDTH16; j++)
            {
                SeqPair sp = pairArray[i + j];
                for(k = sp.len2; k <= maxLen2; k++)
                {
                    mySeq2SoA[k * SIMD_WIDTH16 + j] = DUMMY2;
                }
            }
            /* B5: only the boundary row H1[maxLen2] (value 0) survives the
             * h0-prefix insertion seed below; write just that row, before the
             * seed so its unconditional H1[0]/H1[1] stores still win. */
            _mm_store_si128((__m128i *)(H1 + maxLen2 * SIMD_WIDTH16), _mm_setzero_si128());
//------------------------
            _mm_store_si128((__m128i *) H1, h0_128);
            __m128i cmp128 = _mm_cmpgt_epi16(h0_128, oe_ins128);
            tmp128 = _mm_sub_epi16(h0_128, oe_ins128);

            tmp128 = _mm_blendv_epi16(zero128, tmp128, cmp128);
            _mm_store_si128((__m128i *) (H1 + SIMD_WIDTH16), tmp128);
            for(k = 2; k < maxLen2; k++)
            {
                __m128i h1_128 = tmp128;
                tmp128 = _mm_sub_epi16(h1_128, e_ins128);
                tmp128 = _mm_max_epi16(tmp128, zero128);
                _mm_store_si128((__m128i *)(H1 + k*SIMD_WIDTH16), tmp128);
            }           
//------------------------
            uint16_t myband[SIMD_WIDTH16] __attribute__((aligned(64)));
            uint16_t temp[SIMD_WIDTH16] __attribute__((aligned(64)));
            {
                __m128i qlen128 = _mm_load_si128((__m128i *) qlen);
                __m128i sum128 = _mm_add_epi16(qlen128, eb_ins128);
                _mm_store_si128((__m128i *) temp, sum128);              
                for (int l=0; l<SIMD_WIDTH16; l++) {
                    double val = temp[l]/e_ins + 1.0;
                    int max_ins = (int) val;
                    max_ins = max_ins > 1? max_ins : 1;
                    myband[l] = min_(bsize, max_ins);
                }
                sum128 = _mm_add_epi16(qlen128, eb_del128);
                _mm_store_si128((__m128i *) temp, sum128);              
                for (int l=0; l<SIMD_WIDTH16; l++) {
                    double val = temp[l]/e_del + 1.0;
                    int max_ins = (int) val;
                    max_ins = max_ins > 1? max_ins : 1;
                    myband[l] = min_(myband[l], max_ins);
                    bsize = bsize < myband[l] ? myband[l] : bsize;                  
                }
            }

            smithWaterman128_16(mySeq1SoA,
                                mySeq2SoA,
                                maxLen1,
                                maxLen2,
                                pairArray + i,
                                h0,
                                tid,
                                numPairs,
                                zdrop,
                                bsize,
                                qlen,
                                myband);
        }
    }

#if RDT 
    st4 = ___rdtsc();
#endif
    

#if RDT
    st5 = ___rdtsc();
    setupTicks += st2 - st1;
    sort1Ticks += st3 - st2;
    swTicks += st4 - st3;
    sort2Ticks += st5 - st4;
#endif

    // free mem
    _mm_free(seq1SoA);
    _mm_free(seq2SoA);

    return;
}

void BandedPairWiseSW::smithWaterman128_16(uint16_t seq1SoA[],
                                           uint16_t seq2SoA[],
                                           uint16_t nrow,
                                           uint16_t ncol,
                                           SeqPair *p,
                                           uint16_t h0[],
                                           uint16_t tid,
                                           int32_t numPairs,
                                           int zdrop,
                                           int32_t w,
                                           uint16_t qlen[],
                                           uint16_t myband[])
{
    
    __m128i match128     = _mm_set1_epi16(this->w_match);
    __m128i mismatch128  = _mm_set1_epi16(this->w_mismatch);
    __m128i w_ambig_128  = _mm_set1_epi16(this->w_ambig);   // ambig penalty

    // D3 generic-matrix seam: when the matrix is asymmetric (bisulfite OT/OB),
    // score the ACGT submatrix from a target-major LUT amat[(ref<<2)|read]
    // overlaid on the symmetric base (SBT_PREPASS16_AMAT). The default aligner
    // matrix is symmetric, so gen_mat is false on the hot path and the kernel
    // uses the original cheap SBT_PREPASS16_SYM — no perf cost when not needed.
    const bool forced  = bsw_force_generic_matrix();
    const bool gen_mat = bsw_generic_matrix(this->mat, this->w_match, this->w_mismatch)
                         || forced;
    const BswFreedCell fc = bsw_freed_cell(this->mat, this->w_match, this->w_mismatch, forced);
    __m128i frref128  = _mm_set1_epi16(fc.ref);
    __m128i frread128 = _mm_set1_epi16(fc.read);
    int8_t amat_bytes[16] __attribute__((aligned(16)));
    build_amat16(amat_bytes, this->mat);
    __m128i amat128 = _mm_load_si128((__m128i *)amat_bytes);
    __m128i three128 = _mm_set1_epi16(3);                    // ACGT mask (base <= 3)

    __m128i e_del128    = _mm_set1_epi16(this->e_del);
    __m128i oe_del128   = _mm_set1_epi16(this->o_del + this->e_del);
    __m128i e_ins128    = _mm_set1_epi16(this->e_ins);
    __m128i oe_ins128   = _mm_set1_epi16(this->o_ins + this->e_ins);

    int16_t *F  = F16_ + tid * SIMD_WIDTH16 * MAX_SEQ_LEN16;
    int16_t *H_h    = H16_ + tid * SIMD_WIDTH16 * MAX_SEQ_LEN16;
    int16_t *H_v = H16__ + tid * SIMD_WIDTH16 * MAX_SEQ_LEN16;

    int16_t i, j;

    uint16_t tlen[SIMD_WIDTH16];
    uint16_t tail[SIMD_WIDTH16] __attribute((aligned(64)));
    uint16_t head[SIMD_WIDTH16] __attribute((aligned(64)));

    // PR 17: per-row score-vector scratch (16-bit variant).
    int16_t *sbt_buf = sbt16_ + tid * SIMD_WIDTH16 * MAX_SEQ_LEN16;
    
    int32_t minq = 10000000;
    for (int l=0; l<SIMD_WIDTH16; l++) {
        tlen[l] = p[l].len1;
        qlen[l] = p[l].len2;
        if (p[l].len2 < minq) minq = p[l].len2;
    }
    minq -= 1; // for gscore

    __m128i tlen128 = _mm_load_si128((__m128i *) tlen);
    __m128i qlen128 = _mm_load_si128((__m128i *) qlen);
    __m128i myband128 = _mm_load_si128((__m128i *) myband);
    __m128i zero128 = _mm_setzero_si128();
    __m128i one128  = _mm_set1_epi16(1);
    __m128i two128  = _mm_set1_epi16(2);
    __m128i max_ie128 = zero128;
    __m128i ff128 = _mm_set1_epi16(0xFFFF);
        
    __m128i tail128 = qlen128, head128 = zero128;
    _mm_store_si128((__m128i *) head, head128);
    _mm_store_si128((__m128i *) tail, tail128);

    __m128i mlen128 = _mm_add_epi16(qlen128, myband128);
    mlen128 = _mm_min_epu16(mlen128, tlen128);
        
    __m128i hval = _mm_load_si128((__m128i *)(H_v));

    __mmask16 dmask16 = 0xAAAA;
    
    __m128i maxScore128 = hval;
    for(j = 0; j < ncol; j++)
        _mm_store_si128((__m128i *)(F + j * SIMD_WIDTH16), zero128);
    
    __m128i x128 = zero128;
    __m128i y128 = zero128;
    __m128i gscore = _mm_set1_epi16(-1);
    __m128i max_off128 = zero128;
    __m128i exit0 = _mm_set1_epi16(0xFFFF);
    __m128i zdrop128 = _mm_set1_epi16(zdrop);
    
    int beg = 0, end = ncol;
    int nbeg = beg, nend = end;

#if RDT
    uint64_t tim = __rdtsc();
#endif
    
    for(i = 0; i < nrow; i++)
    {       
        __m128i e11 = zero128;
        __m128i h00, h11, h10;
        __m128i s10 = _mm_load_si128((__m128i *)(seq1SoA + (i + 0) * SIMD_WIDTH16));

        beg = nbeg; end = nend;
        // Banding
        if (beg < i - w) beg = i - w;
        if (end > i + w + 1) end = i + w + 1;
        if (end > ncol) end = ncol;

        h10 = zero128;
        if (beg == 0)
            h10 = _mm_load_si128((__m128i *)(H_v + (i+1) * SIMD_WIDTH16));

        __m128i j128 = zero128;
        __m128i maxRS1 = zero128;
        
        __m128i i1_128 = _mm_set1_epi16(i+1);
        __m128i y1_128 = zero128;
        
#if RDT 
        uint64_t tim1 = __rdtsc();
#endif
        
        __m128i i128, cache128;
        __m128i phead128 = head128, ptail128 = tail128;
        i128 = _mm_set1_epi16(i);
        cache128 = _mm_sub_epi16(i128, myband128);
        head128 = _mm_max_epi16(head128, cache128);
        cache128 = _mm_add_epi16(i1_128, myband128);
        tail128 = _mm_min_epu16(tail128, cache128);
        tail128 = _mm_min_epu16(tail128, qlen128);
        
        // NEW, trimming.
        __m128i cmph = _mm_cmpeq_epi16(head128, phead128);
        __m128i cmpt = _mm_cmpeq_epi16(tail128, ptail128);
        // cmph &= cmpt;
        cmph = _mm_and_si128(cmph, cmpt);
        //__mmask16 cmp_ht = _mm_movepi16_mask(cmph);
        __mmask16 cmp_ht = _mm_movemask_epi8(cmph) & dmask16;
        
        for (int l=beg; l<end && cmp_ht != dmask16; l++)
        {
            __m128i h128 = _mm_load_si128((__m128i *)(H_h + l * SIMD_WIDTH16));
            __m128i f128 = _mm_load_si128((__m128i *)(F + l * SIMD_WIDTH16));
            
            __m128i pj128 = _mm_set1_epi16(l);
            __m128i j128 = _mm_set1_epi16(l+1);
            __m128i cmp1 = _mm_cmpgt_epi16(head128, pj128);
            // uint32_t cval = _mm_movemask_epi16(cmp1);
            // uint16_t cval = _mm_movepi16_mask(cmp1);
            uint16_t cval = _mm_movemask_epi8(cmp1) & dmask16;          
            if (cval == 0x00) break;
            // __m128i cmp2 = _mm_cmpgt_epi16(pj128, tail128);
            __m128i cmp2 = _mm_cmpgt_epi16(j128, tail128);
            cmp1 = _mm_or_si128(cmp1, cmp2);
            h128 = _mm_blendv_epi16(h128, zero128, cmp1);
            f128 = _mm_blendv_epi16(f128, zero128, cmp1);
            
            _mm_store_si128((__m128i *)(F + l * SIMD_WIDTH16), f128);
            _mm_store_si128((__m128i *)(H_h + l * SIMD_WIDTH16), h128);
        }

#if RDT
        prof[DP3][0] += __rdtsc() - tim1;
#endif

        __m128i cmpim = _mm_cmpgt_epi16(i1_128, mlen128);
        __m128i cmpht = _mm_cmpeq_epi16(tail128, head128);
        cmpim = _mm_or_si128(cmpim, cmpht);
        // NEW
        cmpht = _mm_cmpgt_epi16(head128, tail128);
        cmpim = _mm_or_si128(cmpim, cmpht);

        exit0 = _mm_blendv_epi16(exit0, zero128, cmpim);

#if RDT
        tim1 = __rdtsc();
#endif
        
        // PR 17: 16-bit pre-pass. Same shape as the 8-bit path. The gen_mat
        // branch is loop-invariant across the whole kernel call (once per row,
        // perfectly predicted), so the symmetric default keeps its fast prepass.
        if (!gen_mat) {
            for (int jp = beg; jp < end; jp++) {
                __m128i s2 = _mm_load_si128((__m128i *)(seq2SoA + jp * SIMD_WIDTH16));
                __m128i sbt11;
                SBT_PREPASS16_SYM(s10, s2, sbt11, mismatch128, match128, w_ambig_128, ff128);
                _mm_store_si128((__m128i *)(sbt_buf + jp * SIMD_WIDTH16), sbt11);
            }
        } else if (fc.rank1) {
            // Hoist the loop-invariant freed-cell row test: alt10 = (ref==fr_ref)
            // ? fr_read : ref, computed once per row (s10 is constant across band).
            __m128i alt10 = _mm_blendv_epi16(s10, frread128, _mm_cmpeq_epi16(s10, frref128));
            for (int jp = beg; jp < end; jp++) {
                __m128i s2 = _mm_load_si128((__m128i *)(seq2SoA + jp * SIMD_WIDTH16));
                __m128i sbt11;
                SBT_PREPASS16_RANK1(s10, s2, alt10, sbt11, mismatch128, match128, w_ambig_128, ff128);
                _mm_store_si128((__m128i *)(sbt_buf + jp * SIMD_WIDTH16), sbt11);
            }
        } else {
            for (int jp = beg; jp < end; jp++) {
                __m128i s2 = _mm_load_si128((__m128i *)(seq2SoA + jp * SIMD_WIDTH16));
                __m128i sbt11;
                SBT_PREPASS16_AMAT(s10, s2, sbt11, amat128, mismatch128, w_ambig_128, ff128, three128);
                _mm_store_si128((__m128i *)(sbt_buf + jp * SIMD_WIDTH16), sbt11);
            }
        }

        j128 = _mm_set1_epi16(beg);
        for(j = beg; j < end; j++)
        {
            __m128i f11, f21, sbt11;
            h00 = _mm_load_si128((__m128i *)(H_h + j * SIMD_WIDTH16));
            f11 = _mm_load_si128((__m128i *)(F + j * SIMD_WIDTH16));
            sbt11 = _mm_load_si128((__m128i *)(sbt_buf + j * SIMD_WIDTH16));

            __m128i pj128 = j128;
            j128 = _mm_add_epi16(j128, one128);

            MAIN_CODE16_CORE(sbt11, h00, h11, e11, f11, f21, zero128,
                             e_ins128, oe_ins128,
                             e_del128, oe_del128);

            // Masked writing
            __m128i cmp1 = _mm_cmpgt_epi16(head128, pj128);
            __m128i cmp2 = _mm_cmpgt_epi16(pj128, tail128);
            cmp1 = _mm_or_si128(cmp1, cmp2);
            h10 = _mm_blendv_epi16(h10, zero128, cmp1);
            f21 = _mm_blendv_epi16(f21, zero128, cmp1);
            
            __m128i bmaxRS = maxRS1;                                        
            maxRS1 =_mm_max_epi16(maxRS1, h11);                         
            __m128i cmpA = _mm_cmpgt_epi16(maxRS1, bmaxRS);                 
            __m128i cmpB =_mm_cmpeq_epi16(maxRS1, h11);                 
            cmpA = _mm_or_si128(cmpA, cmpB);
            cmp1 = _mm_cmpgt_epi16(j128, tail128); // change
            cmp1 = _mm_or_si128(cmp1, cmp2);            // change           
            cmpA = _mm_blendv_epi16(y1_128, j128, cmpA);
            y1_128 = _mm_blendv_epi16(cmpA, y1_128, cmp1);
            maxRS1 = _mm_blendv_epi16(maxRS1, bmaxRS, cmp1);                        

            _mm_store_si128((__m128i *)(F + j * SIMD_WIDTH16), f21);
            _mm_store_si128((__m128i *)(H_h + j * SIMD_WIDTH16), h10);

            h10 = h11;
            
            // gscore calculations
            if (j >= minq)
            {
                __m128i cmp = _mm_cmpeq_epi16(j128, qlen128);
                __m128i max_gh = _mm_max_epi16(gscore, h11);
                __m128i cmp_gh = _mm_cmpgt_epi16(gscore, h11);
                __m128i tmp128_1 = _mm_blendv_epi16(i1_128, max_ie128, cmp_gh);

                __m128i tmp128_t = _mm_blendv_epi16(max_ie128, tmp128_1, cmp);
                tmp128_1 = _mm_blendv_epi16(max_ie128, tmp128_t, exit0);                
                
                max_gh = _mm_blendv_epi16(gscore, max_gh, exit0);
                max_gh = _mm_blendv_epi16(gscore, max_gh, cmp);             

                cmp = _mm_cmpgt_epi16(j128, tail128); 
                max_gh = _mm_blendv_epi16(max_gh, gscore, cmp);
                max_ie128 = _mm_blendv_epi16(tmp128_1, max_ie128, cmp);
                gscore = max_gh;
            }
        }
        __m128i cmp1 = _mm_cmpgt_epi16(head128, j128);
        __m128i cmp2 = _mm_cmpgt_epi16(j128, tail128);
        cmp1 = _mm_or_si128(cmp1, cmp2);
        h10 = _mm_blendv_epi16(h10, zero128, cmp1);
            
        _mm_store_si128((__m128i *)(H_h + j * SIMD_WIDTH16), h10);
        _mm_store_si128((__m128i *)(F + j * SIMD_WIDTH16), zero128);
        
        /* exit due to zero score by a row */
        __m128i bmaxScore128 = maxScore128;
        __m128i tmp = _mm_cmpeq_epi16(maxRS1, zero128);
#if defined(__ARM_NEON)
        if (vmaxvq_u16(vreinterpretq_u16_m128i(maxRS1)) == 0) break;
#else
        // uint16_t cval = _mm_movepi16_mask(tmp);
        uint16_t cval = _mm_movemask_epi8(tmp) & dmask16;
        if (cval == dmask16) break;
#endif

        exit0 = _mm_blendv_epi16(exit0, zero128,  tmp);

        __m128i score128 = _mm_max_epi16(maxScore128, maxRS1);
        maxScore128 = _mm_blendv_epi16(maxScore128, score128, exit0);

        __m128i cmp = _mm_cmpgt_epi16(maxScore128, bmaxScore128);
        y128 = _mm_blendv_epi16(y128, y1_128, cmp);
        x128 = _mm_blendv_epi16(x128, i1_128, cmp);     
        // max_off calculations
#if 0
        tmp = _mm_sub_epi16(y1_128, i1_128);
        tmp = _mm_abs_epi16(tmp);
#else
        __m128i ab = _mm_subs_epu16(y1_128, i1_128);
        __m128i ba = _mm_subs_epu16(i1_128, y1_128);
        tmp = _mm_or_si128(ab, ba);
#endif
        __m128i bmax_off128 = max_off128;
        tmp = _mm_max_epi16(max_off128, tmp);
        max_off128 = _mm_blendv_epi16(bmax_off128, tmp, cmp);

        // Z-score
        ZSCORE16(i1_128, y1_128);       

#if RDT
        prof[DP1][0] += __rdtsc() - tim1;
        tim1 = __rdtsc();
#endif
        
        /* Narrowing of the band */
        /* From beg */
        int l;
        for (l = beg; l < end; l++)
        {
            __m128i f128 = _mm_load_si128((__m128i *)(F + l * SIMD_WIDTH16));
            __m128i h128 = _mm_load_si128((__m128i *)(H_h + l * SIMD_WIDTH16));
#if defined(__ARM_NEON)
            if (vmaxvq_u16(vorrq_u16(vreinterpretq_u16_m128i(f128),
                                      vreinterpretq_u16_m128i(h128))) == 0)
                nbeg = l;
            else
                break;
#else
            __m128i tmp = _mm_or_si128(f128, h128);
            tmp = _mm_cmpeq_epi16(tmp, zero128);
            // uint16_t val = _mm_movepi16_mask(tmp);
            uint16_t val = _mm_movemask_epi8(tmp) & dmask16;
            if (val == dmask16) nbeg = l;
            else
                break;
#endif
        }

        /* From end */
        bool flg = 1;
        for (l = end; l >= beg; l--)
        {
            __m128i f128 = _mm_load_si128((__m128i *)(F + l * SIMD_WIDTH16));
            __m128i h128 = _mm_load_si128((__m128i *)(H_h + l * SIMD_WIDTH16));
#if defined(__ARM_NEON)
            if (vmaxvq_u16(vorrq_u16(vreinterpretq_u16_m128i(f128),
                                      vreinterpretq_u16_m128i(h128))) != 0 && flg)
                break;
#else
            __m128i tmp = _mm_or_si128(f128, h128);
            tmp = _mm_cmpeq_epi16(tmp, zero128);
            // uint16_t val = _mm_movepi16_mask(tmp);
            uint16_t val = _mm_movemask_epi8(tmp) & dmask16;
            if (val != dmask16 && flg)
                break;
#endif
        }
        nend = l + 2 < ncol? l + 2: ncol;

        __m128i tmpb = ff128;

        __m128i exit1 = _mm_xor_si128(exit0, ff128);
        __m128i l128 = _mm_set1_epi16(beg);
        for (l = beg; l < end; l++)
        {
            __m128i f128 = _mm_load_si128((__m128i *)(F + l * SIMD_WIDTH16));
            __m128i h128 = _mm_load_si128((__m128i *)(H_h + l * SIMD_WIDTH16));

            __m128i tmp = _mm_or_si128(f128, h128);
            tmp = _mm_or_si128(tmp, exit1);
            tmp = _mm_cmpeq_epi16(tmp, zero128);
#if defined(__ARM_NEON)
            if (vmaxvq_u8(vreinterpretq_u8_m128i(tmp)) == 0) {
                break;
            }
#else
            // uint32_t val = _mm_movemask_epi16(tmp);
            // uint16_t val = _mm_movepi16_mask(tmp);
            uint16_t val = _mm_movemask_epi8(tmp) & dmask16;
            if (val == 0x00) {
                break;
            }
#endif
            tmp = _mm_and_si128(tmp,tmpb);
            //__m128i l128 = _mm_set1_epi16(l+1);
            l128 = _mm_add_epi16(l128, one128);
            // NEW
            head128 = _mm_blendv_epi16(head128, l128, tmp);

            tmpb = tmp;
        }
        // _mm_store_si128((__m128i *) head, head128);

        __m128i  index128 = tail128;
        tmpb = ff128;

        l128 = _mm_set1_epi16(end);
        for (l = end; l >= beg; l--)
        {
            __m128i f128 = _mm_load_si128((__m128i *)(F + l * SIMD_WIDTH16));
            __m128i h128 = _mm_load_si128((__m128i *)(H_h + l * SIMD_WIDTH16));

            __m128i tmp = _mm_or_si128(f128, h128);
            tmp = _mm_or_si128(tmp, exit1);
            tmp = _mm_cmpeq_epi16(tmp, zero128);
#if defined(__ARM_NEON)
            if (vmaxvq_u8(vreinterpretq_u8_m128i(tmp)) == 0) {
                break;
            }
#else
            // uint32_t val = _mm_movemask_epi16(tmp);
            // uint16_t val = _mm_movepi16_mask(tmp);
            uint16_t val = _mm_movemask_epi8(tmp) & dmask16;
            if (val == 0x00)  {
                break;
            }
#endif
            tmp = _mm_and_si128(tmp,tmpb);
            l128 = _mm_sub_epi16(l128, one128);
            // NEW
            index128 = _mm_blendv_epi8(index128, l128, tmp);

            tmpb = tmp;
        }
        index128 = _mm_add_epi16(index128, two128);
        tail128 = _mm_min_epi16(index128, qlen128);

#if RDT
        prof[DP2][0] += __rdtsc() - tim1;
#endif
    }
    
#if RDT
    prof[DP][0] += __rdtsc() - tim;
#endif
    
    int16_t score[SIMD_WIDTH16]  __attribute((aligned(64)));
    _mm_store_si128((__m128i *) score, maxScore128);

    int16_t maxi[SIMD_WIDTH16]  __attribute((aligned(64)));
    _mm_store_si128((__m128i *) maxi, x128);

    int16_t maxj[SIMD_WIDTH16]  __attribute((aligned(64)));
    _mm_store_si128((__m128i *) maxj, y128);

    int16_t max_off_ar[SIMD_WIDTH16]  __attribute((aligned(64)));
    _mm_store_si128((__m128i *) max_off_ar, max_off128);

    int16_t gscore_ar[SIMD_WIDTH16]  __attribute((aligned(64)));
    _mm_store_si128((__m128i *) gscore_ar, gscore);

    int16_t maxie_ar[SIMD_WIDTH16]  __attribute((aligned(64)));
    _mm_store_si128((__m128i *) maxie_ar, max_ie128);
    
    for(i = 0; i < SIMD_WIDTH16; i++)
    {
        p[i].score = score[i];
        p[i].tle = maxi[i];
        p[i].qle = maxj[i];
        p[i].max_off = max_off_ar[i];
        p[i].gscore = gscore_ar[i];
        p[i].gtle = maxie_ar[i];
    }
    
    return;
}

/********************************************************************************/
/* 128-bit - 8 bit version */
// NB: the 128-bit kernel requires SSE4.1 (max_epi8 / min_epi8 / blendv_epi8); the
// lowest x86 build tier is sse41 and arm64 routes through sse2neon, so no
// sub-SSE4.1 polyfill is needed.

// ZSCORE8 is unused in smithWaterman128_8 (z-drop replaced by wide scalar);
// retained here in case a future 128-bit tier reinstates it.
#define ZSCORE8(i4_128, y4_128)                                         \
    {                                                                   \
        __m128i tmpi = _mm_sub_epi8(i4_128, x128);                      \
        __m128i tmpj = _mm_sub_epi8(y4_128, y128);                      \
        cmp = _mm_cmpgt_epi8(tmpi, tmpj);                               \
        score128 = _mm_sub_epi8(maxScore128, maxRS1);                   \
        __m128i insdel = _mm_blendv_epi8(e_ins128, e_del128, cmp);      \
        __m128i sub_a128 = _mm_sub_epi8(tmpi, tmpj);                    \
        __m128i sub_b128 = _mm_sub_epi8(tmpj, tmpi);                    \
        tmp = _mm_blendv_epi8(sub_b128, sub_a128, cmp);                 \
        tmp = _mm_sub_epi8(score128, tmp);                              \
        cmp = _mm_cmpgt_epi8(tmp, zdrop128);                            \
        exit0 = _mm_blendv_epi8(exit0, zero128, cmp);                   \
    }



// --- PR 17/16: SSE2/NEON LUT primitive ---
// Symmetric (default) fast path. With the asymmetric AMBIG encoding (target
// N=4, query N=8), the low 4 bits of (s1 ^ s2) index a 16-byte pmat LUT where:
//   [0]     = w_match   (s1==s2, both ACGT)
//   [1..3]  = w_mismatch (ACGT vs ACGT, XOR ∈ {1,2,3})
//   [4..7]  = w_ambig   (target=N=4, query=ACGT)
//   [8..11] = w_ambig   (target=ACGT, query=N=8)
//   [12]    = w_ambig   (both N)
//   [13..15] = w_ambig   (unreachable; filled for safety)
// Replaces 4-op cmpeq+blendv+max+blendv critical path with 1 XOR + 1
// shuffle_epi8 (TBL on NEON). The XOR index is order-insensitive, so it cannot
// represent an asymmetric matrix — that case uses SBT_PREPASS8_AMAT below.
#define SBT_PREPASS8_XOR(s1, s2, sbt11_out, pmat128)                    \
    {                                                                   \
        __m128i xor_ = _mm_xor_si128(s1, s2);                           \
        sbt11_out = _mm_shuffle_epi8(pmat128, xor_);                    \
    }

// D3 generic-matrix seam (gated on an asymmetric matrix; the default symmetric
// path uses SBT_PREPASS8_XOR above and is untouched). Index the 16-byte
// target-major LUT amat[(ref<<2)|read] (order-sensitive, unlike the symmetric
// XOR), then mask N cells to ambig. s1 = target/ref, s2 = query/read;
// ACGT = 0..3, N = 4 (target) / 8 (query). For a symmetric matrix amat yields
// the same scores as the XOR pmat path; an asymmetric off-diagonal (ref-C x
// read-T freed) is scored correctly. (s1<<2 via two add_epi8 — no byte shift in
// SSE2/NEON; N detected by max_epu8(s1,s2) > 3.)
#define SBT_PREPASS8_AMAT(s1, s2, sbt11_out, amat128, ambig128, three128) \
    {                                                                   \
        __m128i sh_  = _mm_add_epi8(s1, s1);                            \
        sh_          = _mm_add_epi8(sh_, sh_);          /* s1 << 2 */   \
        __m128i idx_ = _mm_or_si128(sh_, s2);           /* (ref<<2)|read */ \
        __m128i acgt_  = _mm_shuffle_epi8(amat128, idx_);              \
        __m128i nmax_  = _mm_max_epu8(s1, s2);                          \
        __m128i nmask_ = _mm_cmpgt_epi8(nmax_, three128); /* N: base > 3 */ \
        sbt11_out = _mm_or_si128(_mm_andnot_si128(nmask_, acgt_),       \
                                 _mm_and_si128(nmask_, ambig128));      \
    }

// D3 8-bit rank-1 fast path: symmetric XOR LUT (folds match/mismatch + N) plus a
// single freed-to-match cell override (bisulfite OT/OB). Cheaper than AMAT on the
// fast 8-bit tier (no index build, no max/cmpgt N-mask): the XOR LUT already
// handles N/ambig, so we only override the one freed ACGT cell to match. frref/
// frread are the 8-bit-broadcast freed (ref,read); padding/N never equal frref.
// rowfreed = (ref==fr_ref) hoisted per row (s1 is constant across the band loop).
#define SBT_PREPASS8_RANK1(s1, s2, rowfreed, sbt11_out, pmat128, match128, frread128) \
    {                                                                   \
        __m128i xor_  = _mm_xor_si128(s1, s2);                          \
        __m128i sbt_  = _mm_shuffle_epi8(pmat128, xor_);                \
        __m128i freed_ = _mm_and_si128(rowfreed, _mm_cmpeq_epi8(s2, frread128)); \
        sbt11_out = _mm_blendv_epi8(sbt_, match128, freed_);            \
    }

// MAIN_CODE8_CORE runs only the cell-update half of MAIN_CODE8 using a
// pre-computed score vector `sbt11` from SBT_PREPASS8.
#define MAIN_CODE8_CORE(sbt11, h00, h11, e11, f11, f21, zero128, e_ins128, oe_ins128, e_del128, oe_del128) \
    {                                                                   \
        /* M = max(0, h00 + sbt) computed in UNSIGNED-saturating form so a    \
         * legitimate score in [128,255] is kept (the old signed `m11 &       \
         * (m11 > 0)` floor mis-read >127 as negative and zeroed it). Split   \
         * the signed substitution score into its +bonus and -penalty parts:  \
         * adds_epu8 (no wrap past 255) then subs_epu8 (floors at 0). */       \
        __m128i sbt_pos = _mm_max_epi8(sbt11, zero128);                 \
        __m128i sbt_neg = _mm_max_epi8(_mm_sub_epi8(zero128, sbt11), zero128); \
        __m128i m11 = _mm_subs_epu8(_mm_adds_epu8(h00, sbt_pos), sbt_neg); \
        __m128i cmp11 = _mm_cmpeq_epi8(h00, zero128);                   \
        m11 = _mm_blendv_epi8(m11, zero128, cmp11);  /* h00==0 -> local restart */ \
        h11 = _mm_max_epu8(m11, e11);                                   \
        h11 = _mm_max_epu8(h11, f11);                                   \
        /* Reassociate gap recurrences off h11 (variant A, saturating u8). \
         * me reads OLD e11 — compute before the e-update. */               \
        __m128i mf128 = _mm_max_epu8(m11, f11);                         \
        __m128i me128 = _mm_max_epu8(m11, e11);                         \
        __m128i temp128 = _mm_subs_epu8(mf128, oe_ins128);             \
        e11 = _mm_subs_epu8(e11, e_ins128);                            \
        e11 = _mm_max_epu8(temp128, e11);                               \
        temp128 = _mm_subs_epu8(me128, oe_del128);                     \
        f21 = _mm_subs_epu8(f11, e_del128);                            \
        f21 = _mm_max_epu8(temp128, f21);                               \
    }


// #define PFD 2 // SSE2
void BandedPairWiseSW::getScores8(SeqPair *pairArray,
                                  uint8_t *seqBufRef,
                                  uint8_t *seqBufQer,
                                  int32_t numPairs,
                                  uint16_t numThreads,
                                  int32_t w)
{
    assert(SIMD_WIDTH8 == 16 && SIMD_WIDTH16 == 8);
    {
        BswOvershootGuard _g(pairArray, numPairs, SIMD_WIDTH8, guard_overshoot_);
        smithWatermanBatchWrapper8(pairArray, seqBufRef, seqBufQer, numPairs, numThreads, w);
    }

#if MAXI
    printf("Vecor code: Writing output..\n");
    for (int l=0; l<numPairs; l++)
    {
        fprintf(stderr, "%d (%d %d) %d %d %d\n",
                pairArray[l].score, pairArray[l].x, pairArray[l].y,
                pairArray[l].gscore, pairArray[l].max_off, pairArray[l].max_ie);

    }
    printf("Vector code: Writing output completed!!!\n\n");
#endif
    
}

void BandedPairWiseSW::smithWatermanBatchWrapper8(SeqPair *pairArray,
                                                  uint8_t *seqBufRef,
                                                  uint8_t *seqBufQer,
                                                  int32_t numPairs,
                                                  uint16_t numThreads,
                                                  int32_t w)
{
    numThreads = effective_threads(numThreads);
#if RDT
    int64_t st1, st2, st3, st4, st5;
    st1 = ___rdtsc();
#endif
    uint8_t *seq1SoA = (uint8_t *)_mm_malloc((size_t)MAX_SEQ_LEN8 * SIMD_WIDTH8 * numThreads * sizeof(uint8_t), 64);
    uint8_t *seq2SoA = (uint8_t *)_mm_malloc((size_t)MAX_SEQ_LEN8 * SIMD_WIDTH8 * numThreads * sizeof(uint8_t), 64);

    if (UNLIKELY(seq1SoA == NULL || seq2SoA == NULL)) {
        fprintf(stderr, "Error! Mem not allocated!!!\n");
        exit(EXIT_FAILURE);
    }
    
    int32_t ii;
    int32_t roundNumPairs = ((numPairs + SIMD_WIDTH8 - 1)/SIMD_WIDTH8 ) * SIMD_WIDTH8;
    // assert(roundNumPairs < BATCH_SIZE * SEEDS_PER_READ);
    for(ii = numPairs; ii < roundNumPairs; ii++)
    {
        pairArray[ii].id = ii;
        pairArray[ii].len1 = 0;
        pairArray[ii].len2 = 0;
    }

#if RDT
    st2 = ___rdtsc();
#endif
    

#if RDT
    st3 = ___rdtsc();
#endif

    int eb = end_bonus;
// #pragma omp parallel num_threads(numThreads)
    {
        int32_t i;
        uint16_t tid =  0; 
        uint8_t *mySeq1SoA = seq1SoA + tid * MAX_SEQ_LEN8 * SIMD_WIDTH8;
        uint8_t *mySeq2SoA = seq2SoA + tid * MAX_SEQ_LEN8 * SIMD_WIDTH8;
        assert(mySeq1SoA != NULL && mySeq2SoA != NULL);     
        uint8_t *seq1;
        uint8_t *seq2;
        uint8_t h0[SIMD_WIDTH8]   __attribute__((aligned(64)));
        uint8_t qlen[SIMD_WIDTH8] __attribute__((aligned(64)));
        int32_t bsize = 0;

        int8_t *H1 = H8_ + tid * SIMD_WIDTH8 * MAX_SEQ_LEN8;
        int8_t *H2 = H8__ + tid * SIMD_WIDTH8 * MAX_SEQ_LEN8;

        __m128i zero128   = _mm_setzero_si128();
        __m128i e_ins128  = _mm_set1_epi8(e_ins);
        __m128i oe_ins128 = _mm_set1_epi8(o_ins + e_ins);
        __m128i o_del128  = _mm_set1_epi8(o_del);
        __m128i e_del128  = _mm_set1_epi8(e_del);
        __m128i eb_ins128 = _mm_set1_epi8(eb - o_ins);
        __m128i eb_del128 = _mm_set1_epi8(eb - o_del);

        int8_t max = 0;
        if (max < w_match) max = w_match;
        if (max < w_mismatch) max = w_mismatch;
        if (max < w_ambig) max = w_ambig;

        // Initial per-lane score floor B0 = max(0, h0 - REBASE_KEEP_W) (long-read
        // 8-bit, seed score h0 >= 128). Seeding the H bytes from raw h0 overflows
        // the signed-int8 score state on row 0, so we seed from the rebaselined
        // h0' = min(h0, REBASE_KEEP_W) instead (see the per-lane h0 init below).
        // REBASE_KEEP_W MUST equal the kernel's REBASE_KEEP (bsw8_rebase_keep)
        // so the two agree on B0; smithWaterman128_8 recomputes B0 via
        // bsw8_initial_floor(p[].h0, zdrop).
        const int REBASE_KEEP_W = bsw8_rebase_keep(zdrop);
        // The h0-prefix column/row seed below is unsigned-saturating [0,255], so
        // the seeded byte min(h0, REBASE_KEEP_W) only needs to fit a uint8. The
        // routing gate (bsw8_envelope_ok) enforces zdrop + maxStep <= 253 (so the
        // re-baseline window REBASE_HI = 255 - maxStep > REBASE_KEEP holds);
        // assert the byte bound for any direct caller that bypasses the gate.
        assert(REBASE_KEEP_W <= 255 &&
               "8-bit h0-prefix seed byte min(h0,REBASE_KEEP) must fit uint8 (zdrop<=254)");

        int nstart = 0, nend = numPairs;
        
// #pragma omp for schedule(dynamic, 128)
        for(i = nstart; i < nend; i+=SIMD_WIDTH8)
        {
            int32_t j, k;
            int maxLen1 = 0;
            int maxLen2 = 0;
            //bsize = 100;
            bsize = w;
            
            for(j = 0; j < SIMD_WIDTH8; j++)
            {
                SeqPair sp = pairArray[i + j];
                // Re-baseline the seed score into the int8 byte frame: seed the H
                // arrays from h0' = h0 - B0 = min(h0, REBASE_KEEP_W) (clamped >= 0)
                // so the initial bytes fit signed-int8 even when the true seed h0
                // is several hundred (long-read). smithWaterman128_8 recomputes the
                // matching B0 = max(0, h0 - REBASE_KEEP) from p[].h0 + zdrop and
                // starts that lane's floor B there, so byte == H_absolute - B.
                {
                    int h0p = sp.h0;
                    if (h0p > REBASE_KEEP_W) h0p = REBASE_KEEP_W;
                    if (h0p < 0) h0p = 0;
                    // Always-on uint8 guard (asserts are compiled out in release):
                    // a direct caller that bypasses the gate with zdrop > 254 would
                    // otherwise truncate the seed byte. In-envelope REBASE_KEEP_W <= 253,
                    // so this never fires on the production path.
                    if (h0p > 255) h0p = 255;
                    h0[j] = (uint8_t) h0p;
                }
                seq1 = seqBufRef + (int64_t)sp.idr;

                for(k = 0; k < sp.len1; k++)
                {
                    mySeq1SoA[k * SIMD_WIDTH8 + j] = seq1[k] /* PR16: N stays 4 */;
                }
                qlen[j] = sp.len2 * max;
                if(maxLen1 < sp.len1) maxLen1 = sp.len1;
            }

            for(j = 0; j < SIMD_WIDTH8; j++)
            {
                SeqPair sp = pairArray[i + j];
                for(k = sp.len1; k <= maxLen1; k++) //removed "="
                {
                    mySeq1SoA[k * SIMD_WIDTH8 + j] = DUMMY1;
                }
            }
            /* B5: the h0-prefix deletion seed below fully overwrites H2 rows
             * [0, maxLen1); only the boundary row H2[maxLen1] survives to be read
             * as the column edge. Write just that row (all lanes = DUMMY1) here
             * instead of the dead per-row 0/DUMMY fills removed above. Placed
             * before the seed so the seed's unconditional H2[0] store still wins,
             * matching the original store ordering for every maxLen1. */
            _mm_store_si128((__m128i *)(H2 + maxLen1 * SIMD_WIDTH8), _mm_set1_epi8((char)DUMMY1));
//--------------------
            __m128i h0_128 = _mm_load_si128((__m128i*) h0);
            _mm_store_si128((__m128i *) H2, h0_128);
            __m128i tmp128 = _mm_subs_epu8(h0_128, o_del128);

            for(k = 1; k < maxLen1; k++)
            {
                tmp128 = _mm_subs_epu8(tmp128, e_del128);
                //__m128i tmp128_ = _mm_max_epi8(tmp128, zero128);    //epi is not present in SSE2
                _mm_store_si128((__m128i *)(H2 + k* SIMD_WIDTH8), tmp128);
            }
//-------------------

            for(j = 0; j < SIMD_WIDTH8; j++)
            {               
                SeqPair sp = pairArray[i + j];
                // seq2 = seqBuf + (2 * (int64_t)sp.id + 1) * MAX_SEQ_LEN;
                seq2 = seqBufQer + (int64_t)sp.idq;
                
                for(k = 0; k < sp.len2; k++)
                {
                    mySeq2SoA[k * SIMD_WIDTH8 + j] = (seq2[k]==AMBIG ? 8 : seq2[k]) /* PR16: query N→8 */;
                }
                if(maxLen2 < sp.len2) maxLen2 = sp.len2;
            }

            //maxLen2 = ((maxLen2  + 3) >> 2) * 4;

            for(j = 0; j < SIMD_WIDTH8; j++)
            {
                SeqPair sp = pairArray[i + j];
                for(k = sp.len2; k <= maxLen2; k++)
                {
                    mySeq2SoA[k * SIMD_WIDTH8 + j] = DUMMY2;
                }
            }
            /* B5: the h0-prefix insertion seed below fully overwrites H1 rows
             * [0, maxLen2); only the boundary row H1[maxLen2] survives as the row
             * edge (value 0). Write just that row here instead of the dead
             * per-row zero fills removed above; before the seed so its
             * unconditional H1[0]/H1[1] stores still win for small maxLen2. */
            _mm_store_si128((__m128i *)(H1 + maxLen2 * SIMD_WIDTH8), _mm_setzero_si128());
//------------------------
            _mm_store_si128((__m128i *) H1, h0_128);
            // h0-prefix insertion seed, unsigned-saturating [0,255]:
            // H1[1] = max(0, h0' - oe_ins), then -e_ins per step. The first gap
            // was the last signed island here (cmpgt_epi8 + sub_epi8 + blendv);
            // the extension loop already used subs_epu8.
            tmp128 = _mm_subs_epu8(h0_128, oe_ins128);
            _mm_store_si128((__m128i *) (H1 + SIMD_WIDTH8), tmp128);
            for(k = 2; k < maxLen2; k++)
            {
                tmp128 = _mm_subs_epu8(tmp128, e_ins128);
                _mm_store_si128((__m128i *)(H1 + k*SIMD_WIDTH8), tmp128);
            }
//------------------------
            uint8_t myband[SIMD_WIDTH8] __attribute__((aligned(64)));
            uint8_t temp[SIMD_WIDTH8] __attribute__((aligned(64)));
            {
                __m128i qlen128 = _mm_load_si128((__m128i *) qlen);
                __m128i sum128 = _mm_add_epi8(qlen128, eb_ins128);
                _mm_store_si128((__m128i *) temp, sum128);              
                for (int l=0; l<SIMD_WIDTH8; l++) {
                    double val = temp[l]/e_ins + 1.0;
                    int max_ins = (int) val;
                    max_ins = max_ins > 1? max_ins : 1;
                    myband[l] = min_(bsize, max_ins);
                }
                sum128 = _mm_add_epi8(qlen128, eb_del128);
                _mm_store_si128((__m128i *) temp, sum128);              
                for (int l=0; l<SIMD_WIDTH8; l++) {
                    double val = temp[l]/e_del + 1.0;
                    int max_ins = (int) val;
                    max_ins = max_ins > 1? max_ins : 1;
                    myband[l] = min_(myband[l], max_ins);
                    bsize = bsize < myband[l] ? myband[l] : bsize;
                }
            }

            smithWaterman128_8(mySeq1SoA,
                               mySeq2SoA,
                               maxLen1,
                               maxLen2,
                               pairArray + i,
                               h0,
                               tid,
                               numPairs,
                               zdrop,
                               bsize,
                               myband);         
        }
    }
#if RDT
     st4 = ___rdtsc();
#endif
     

#if RDT
    st5 = ___rdtsc();
    setupTicks = st2 - st1;
    sort1Ticks = st3 - st2;
    swTicks = st4 - st3;
    sort2Ticks = st5 - st4;
#endif
    
    // free mem
    _mm_free(seq1SoA);
    _mm_free(seq2SoA);
    
    return;
}

void BandedPairWiseSW::smithWaterman128_8(uint8_t seq1SoA[],
                                          uint8_t seq2SoA[],
                                          int nrow,
                                          int ncol,
                                          SeqPair *p,
                                          uint8_t h0[],
                                          uint16_t tid,
                                          int32_t numPairs,
                                          int zdrop,
                                          int32_t w,
                                          uint8_t myband[])
{
    
    __m128i match128     = _mm_set1_epi8(this->w_match);
    __m128i mismatch128  = _mm_set1_epi8(this->w_mismatch);
    __m128i w_ambig_128  = _mm_set1_epi8(this->w_ambig);    // ambig penalty

    // D3 generic-matrix seam: the default symmetric matrix uses the cheap XOR
    // LUT (SBT_PREPASS8_XOR, pmat128); an asymmetric matrix (bisulfite OT/OB)
    // uses the target-major LUT amat[(ref<<2)|read] (SBT_PREPASS8_AMAT). gen_mat
    // is false on the hot path, so symmetric scoring keeps its original speed.
    const bool forced  = bsw_force_generic_matrix();
    const bool gen_mat = bsw_generic_matrix(this->mat, this->w_match, this->w_mismatch)
                         || forced;
    const BswFreedCell fc = bsw_freed_cell(this->mat, this->w_match, this->w_mismatch, forced);
    __m128i frref128  = _mm_set1_epi8(fc.ref);
    __m128i frread128 = _mm_set1_epi8(fc.read);
    int8_t pmat_bytes[16] __attribute__((aligned(16)));
    build_pmat16(pmat_bytes, this->w_match, this->w_mismatch, this->w_ambig);
    __m128i pmat128 = _mm_load_si128((__m128i *)pmat_bytes);
    int8_t amat_bytes[16] __attribute__((aligned(16)));
    build_amat16(amat_bytes, this->mat);
    __m128i amat128 = _mm_load_si128((__m128i *)amat_bytes);
    __m128i three128 = _mm_set1_epi8(3);                     // N threshold (base > 3)

    __m128i e_del128    = _mm_set1_epi8(this->e_del);
    __m128i oe_del128   = _mm_set1_epi8(this->o_del + this->e_del);
    __m128i e_ins128    = _mm_set1_epi8(this->e_ins);
    __m128i oe_ins128   = _mm_set1_epi8(this->o_ins + this->e_ins);

    int8_t  *F   = F8_ + tid * SIMD_WIDTH8 * MAX_SEQ_LEN8;
    int8_t  *H_h = H8_ + tid * SIMD_WIDTH8 * MAX_SEQ_LEN8;
    int8_t  *H_v = H8__ + tid * SIMD_WIDTH8 * MAX_SEQ_LEN8;

    int i, j;

    uint8_t tail[SIMD_WIDTH8] __attribute((aligned(64)));
    uint8_t head[SIMD_WIDTH8] __attribute((aligned(64)));

    // PR 17: per-row score-vector scratch for loop fission. Holds sbt[j]
    // (score for each band column) pre-computed outside the DP core so
    // the core's critical path doesn't carry the cmpeq+blendv+maxu+blendv
    // chain that builds the score. Slab-backed (see constructor sbt8_).
    int8_t *sbt_buf = sbt8_ + tid * SIMD_WIDTH8 * MAX_SEQ_LEN8;

    // --- DIAGONAL-OFFSET POSITION ENCODING (long-read 8-bit, w<=127) ---
    // Every per-cell COLUMN position is tracked as the diagonal offset
    //   d = col - i  in [-w, +w+1]  (fits signed int8 for w <= ~126).
    // ROW quantities (best row, best-gscore row, qlen, tlen, mlen) exceed
    // int8 for long reads, so they live in WIDE per-lane int32 side channels
    // updated O(rows) in the per-row epilogue, not O(cells). Absolute end
    // coordinates are reconstructed at the result store from the wide row.
    // Persistent column-offset state (head128/tail128) is shifted by -1 each
    // row (frame follows i) so the same absolute edge keeps its offset.
    int32_t tlenw[SIMD_WIDTH8];   // raw target length (rows), wide
    int32_t qlenw[SIMD_WIDTH8];   // raw query length (cols), wide
    int32_t mbandw[SIMD_WIDTH8];  // per-lane band width, wide
    int32_t mlenw[SIMD_WIDTH8];   // min(qlen+myband, tlen), wide row bound
    int32_t xrow[SIMD_WIDTH8];    // best row for score (== i+1 at capture)
    int32_t ierow[SIMD_WIDTH8];   // best row for gscore (== i+1 at capture)

    // --- PER-LANE SCORE RE-BASELINING (long-read 8-bit, scores > 127) ---
    // H/F/E and the score maxima are unsigned 8-bit [0,255], but a 1 kb
    // alignment scores ~900. To keep the kernel 16-lane (8-bit) we carry a
    // wide per-lane running FLOOR B[l]: every stored 8-bit score represents
    // (H_absolute - B[l]). When a lane's row-max climbs toward 255 we raise
    // B[l] by `delta` and subtract `delta` (saturating, _mm_subs_epu8) from
    // every carried score buffer/register for that lane, then add B back when
    // reconstructing the absolute result.
    //
    // Scores live in the UNSIGNED byte range [0,255]. The DP recurrence computes
    // M = max(0, h00 + sbt) with unsigned-saturating arithmetic (adds_epu8 then
    // subs_epu8) and the row/global-max trackers (maxRS1, maxScore128) compare
    // with UNSIGNED order (== after _mm_max_epu8), so the full [0,255] is usable.
    //
    // CORRECTNESS depends on re-baselining NEVER firing for a pair admitted to
    // this kernel. The saturating-subtract is NOT generally lossless: it can zero
    // a still-positive off-diagonal cell (a cell may sit > zdrop below the ROW max
    // yet still lie on the eventual optimum — z-drop is a row-level early-exit, not
    // a per-cell guarantee), and a zeroed cell is then misread as the h00==0
    // local-restart sentinel, spuriously truncating a valid alignment. So we do
    // NOT rely on re-baseline being lossless; instead it is held inert:
    //   REBASE_HI = 255 - maxStep: a row whose max has not reached REBASE_HI grows
    //     by <= maxStep per row (H(i,j) <= H(i-1,j-1) + maxStep), so its bytes stay
    //     <= 255 with no rebase needed.
    //   The bwamem.cpp routing gate admits a pair only when its MAX ATTAINABLE
    //     score stays < REBASE_HI and h0 <= REBASE_KEEP (so B0 == 0). Under that
    //     gate the row max never reaches REBASE_HI: B stays 0, stored bytes equal
    //     absolute scores, and the kernel is a plain exact unsigned [0,255] SW.
    // The re-baseline machinery below is retained only as a safety net; pairs that
    // could exceed the envelope take the 16-bit path.
    // Precondition REBASE_HI > REBASE_KEEP (a non-empty window) holds for
    // zdrop + maxStep <= 253, which is exactly what the routing gate enforces.
    // The h0-prefix column/row seed (wrapper setup below) is unsigned-saturating
    // [0,255] and imposes no tighter ceiling; it previously used signed int8 ops
    // that required the seed byte <= 127 and capped zdrop at 126.
    int maxStep = (int)this->w_match;
    if ((int)this->w_ambig > maxStep) maxStep = (int)this->w_ambig;
    if (maxStep < 1) maxStep = 1;
    const int REBASE_KEEP = zdrop + 1;
    const int REBASE_HI   = 255 - maxStep;
    assert(REBASE_HI > REBASE_KEEP &&
           "8-bit SW re-baseline: zdrop+maxStep>253, scoring outside safe envelope (route to 16-bit)");
    int32_t B[SIMD_WIDTH8];        // per-lane wide score floor (stored = abs - B)
    int32_t best_abs[SIMD_WIDTH8]; // running best score in ABSOLUTE units
    int32_t gbest_abs[SIMD_WIDTH8];// running gscore (query-end) in ABSOLUTE units

    // --- INITIAL SCORE FLOOR B0 (long-read 8-bit, seed score h0 >= 128) ---
    // The H arrays are seeded from the seed score h0 (H_v[0]=h0, then gap-
    // decremented down column 0 and across row 0). For a long-read seed h0 can
    // be several hundred, which would overflow the unsigned [0,255] byte state
    // before any re-baseline event can fire. So we start each lane's
    // running floor at B0[l] = max(0, h0 - REBASE_KEEP) instead of 0: the
    // wrapper seeds the H bytes from h0' = h0 - B0 = min(h0, REBASE_KEEP)
    // (clamped >= 0, so the initial bytes land in [0, REBASE_KEEP] <= 254), and
    // here we set B[l] = B0[l] so the stored bytes still mean (H_absolute - B).
    //
    // SAFETY (identical argument to the mid-DP re-baseline): the gap-prefix
    // cells in column 0 / row 0 that fall more than REBASE_KEEP below h0 are
    // saturated to byte 0 by the wrapper's unsigned-saturating seed. Because
    // REBASE_KEEP = zdrop+1 > zdrop, any such prefix cell is already beyond the
    // z-drop horizon below the seed, so an alignment routed through it has
    // already z-dropped and cannot be optimal -- clamping it to the floor is
    // lossless. This is the SAME situation that arises mid-DP whenever B>0 after
    // a re-baseline (byte 0 then means H_absolute == B, not 0), which the
    // h00==0 local-restart test in MAIN_CODE8_CORE already handles correctly; B0
    // simply makes B>0 from row 0 instead of from the first re-baseline event.
    int32_t B0[SIMD_WIDTH8];
    for (int l=0; l<SIMD_WIDTH8; l++) {
        B0[l] = bsw8_initial_floor((int)p[l].h0, zdrop);
    }

    int32_t minq = 10000000;
    for (int l=0; l<SIMD_WIDTH8; l++) {
        tlenw[l]  = p[l].len1;
        qlenw[l]  = p[l].len2;
        mbandw[l] = myband[l];
        int ml = qlenw[l] + mbandw[l];
        if (ml > tlenw[l]) ml = tlenw[l];
        mlenw[l]  = ml;
        xrow[l]   = 0;
        ierow[l]  = 0;
        B[l]        = B0[l];   // start the floor at B0 so byte state = abs - B0
        best_abs[l] = p[l].h0; // maxScore128 inits to the h0 seed; record the
                               // ABSOLUTE seed score (wide), not the rebaselined byte
        gbest_abs[l]= -1;      // unset sentinel (-1): gscore=-1 / gtle=0 when no query end is reached, matching scalar
        if (p[l].len2 < minq) minq = p[l].len2;
    }
    minq -= 1; // for gscore

    __m128i myband128 = _mm_load_si128((__m128i *) myband);
    __m128i zero128   = _mm_setzero_si128();
    __m128i one128    = _mm_set1_epi8(1);
    __m128i two128    = _mm_set1_epi8(2);
    __m128i ff128     = _mm_set1_epi8(0xFF);

    // gscore (query-end score) capture, per row. The diagonal-offset 8-bit kernel
    // CANNOT reconstruct gscore from the rebaselined byte state via the trimmed
    // tail128: re-baselining saturates the off-diagonal query-end cells to 0 (they
    // sit far below the row max), the band masking then trims them, and the gscore
    // is lost. Instead, when a lane's band-grown band reaches the query end we
    // capture that lane's query-end cell H (hqe128) and flag it (qfire128) here in
    // the inner loop -- BEFORE this row's re-baseline -- and finalize a per-lane
    // WIDE running gbest_abs / ierow in the epilogue (value = byte + B). Reset per row.
    __m128i hqe128   = zero128;   // query-end cell H (rebaselined byte) this row
    __m128i qfire128 = zero128;   // 0xFF where this lane reached its query end this row

    // Offset-frame band edges. head_off starts at 0 (col 0 - row 0); tail_off
    // starts saturated-high (+127) and is immediately clamped by the band-grow
    // min() against (1+myband) and (qlen-i) on the first row.
    __m128i head128 = zero128;
    __m128i tail128 = _mm_set1_epi8(127);
    _mm_store_si128((__m128i *) head, head128);
    _mm_store_si128((__m128i *) tail, tail128);

    __m128i hval = _mm_load_si128((__m128i *)(H_v));
    __mmask16 dmask = 0xFFFF;

    __m128i maxScore128 = hval;
    for(j = 0; j < ncol; j++)
        _mm_store_si128((__m128i *)(F + j * SIMD_WIDTH8), zero128);

    __m128i y128 = zero128;   // best col as diagonal offset: col - i_at_capture (i_at_capture = xrow[l]-1)
    __m128i max_off128 = zero128;
    __m128i exit0 = _mm_set1_epi8(0xFF);
    __m128i zdrop128 = _mm_set1_epi8(zdrop);

    int beg = 0, end = ncol;
    int nbeg = beg, nend = end;

#if RDT
    uint64_t tim = __rdtsc();
#endif
    
    for(i = 0; i < nrow; i++)
    {       
        __m128i e11 = zero128;
        __m128i h00, h11, h10;
        __m128i s10 = _mm_load_si128((__m128i *)(seq1SoA + (i + 0) * SIMD_WIDTH8));

        beg = nbeg; end = nend;
        // Banding
        if (beg < i - w) beg = i - w;
        if (end > i + w + 1) end = i + w + 1;
        if (end > ncol) end = ncol;

        h10 = zero128;
        if (beg == 0)
            h10 = _mm_load_si128((__m128i *)(H_v + (i+1) * SIMD_WIDTH8));

        __m128i j128 = zero128;
        __m128i maxRS1 = zero128;

        __m128i y1_128 = zero128;   // row-max column as diagonal offset (col - i)

        // gscore query-end capture resets each row.
        hqe128   = zero128;
        qfire128 = zero128;

        // Per-row diagonal-offset of the query end: qlen_off = qlen - i. Built
        // wide then saturated to int8, with a validity mask so an out-of-band
        // qlen-i (which would wrap and spuriously alias an in-band col offset)
        // never triggers a false gscore/clamp. qlen_off_valid is true only when
        // qlen-i lies within the representable band window [-w, w+1].
        int8_t qlen_off_a[SIMD_WIDTH8]   __attribute((aligned(16)));
        int8_t qlen_valid_a[SIMD_WIDTH8] __attribute((aligned(16)));
        int8_t cmpim_a[SIMD_WIDTH8]      __attribute((aligned(16)));
        for (int l=0; l<SIMD_WIDTH8; l++) {
            int qoff = qlenw[l] - i;                 // wide
            int qoff_sat = qoff;                     // saturate to int8 range
            if (qoff_sat >  127) qoff_sat =  127;
            if (qoff_sat < -128) qoff_sat = -128;
            qlen_off_a[l]   = (int8_t) qoff_sat;
            qlen_valid_a[l] = (qoff >= -w && qoff <= w + 1) ? (int8_t)0xFF : 0;
            // exit when row i+1 has passed the per-lane effective length bound.
            cmpim_a[l]      = ((i + 1) > mlenw[l]) ? (int8_t)0xFF : 0;
        }
        __m128i qlen_off128   = _mm_load_si128((__m128i *) qlen_off_a);
        __m128i qlen_valid128 = _mm_load_si128((__m128i *) qlen_valid_a);

#if RDT
        uint64_t tim1 = __rdtsc();
#endif

        // Banding (diagonal-offset frame). head128/tail128 arrive here in row-i's
        // offset frame: the -1 shift at the end of the previous iteration converts
        // (col - (i-1)) -> (col - i), so no additional adjustment is needed here.
        //   abs head-grow: head = max(head, i - myband)  ->  head_off = max(head_off, -myband)
        //   abs tail-clamp: tail = min(tail, (i+1)+myband, qlen)
        //                                          -> tail_off = min(tail_off, 1+myband, qlen-i)
        __m128i cache128;
        __m128i phead128 = head128, ptail128 = tail128;
        __m128i negband128 = _mm_sub_epi8(zero128, myband128);            // -myband
        head128 = _mm_max_epi8(head128, negband128);
        cache128 = _mm_add_epi8(myband128, one128);                       // 1 + myband
        tail128 = _mm_min_epi8(tail128, cache128);
        tail128 = _mm_min_epi8(tail128, qlen_off128);

        // NEW, trimming.
        __m128i cmph = _mm_cmpeq_epi8(head128, phead128);
        __m128i cmpt = _mm_cmpeq_epi8(tail128, ptail128);
        // cmph &= cmpt;
        cmph = _mm_and_si128(cmph, cmpt);
        // __mmask32 cmp_ht = _mm_movemask_epi8(cmph);
        __mmask16 cmp_ht = _mm_movemask_epi8(cmph);

        for (int l=beg; l<end && cmp_ht != dmask; l++)
        {
            __m128i h128 = _mm_load_si128((__m128i *)(H_h + l * SIMD_WIDTH8));
            __m128i f128 = _mm_load_si128((__m128i *)(F + l * SIMD_WIDTH8));

            __m128i pj128 = _mm_set1_epi8(l - i);   // diagonal offset of column l
            __m128i cmp1 = _mm_cmpgt_epi8(head128, pj128);
            uint32_t cval = _mm_movemask_epi8(cmp1);
            if (cval == 0x00) break;
            __m128i cmp2 = _mm_cmpgt_epi8(pj128, tail128);
            cmp1 = _mm_or_si128(cmp1, cmp2);
            h128 = _mm_blendv_epi8(h128, zero128, cmp1);
            f128 = _mm_blendv_epi8(f128, zero128, cmp1);

            _mm_store_si128((__m128i *)(F + l * SIMD_WIDTH8), f128);
            _mm_store_si128((__m128i *)(H_h + l * SIMD_WIDTH8), h128);
        }

#if RDT
        prof[DP3][0] += __rdtsc() - tim1;
#endif

        // cmpim: lane exits if row i+1 passed its effective length (precomputed
        // wide into cmpim_a) OR the band collapsed (tail<=head, offset frame).
        __m128i cmpim = _mm_load_si128((__m128i *) cmpim_a);
        __m128i cmpht = _mm_cmpeq_epi8(tail128, head128);
        cmpim = _mm_or_si128(cmpim, cmpht);
        // NEW
        cmpht = _mm_cmpgt_epi8(head128, tail128);
        cmpim = _mm_or_si128(cmpim, cmpht);

        exit0 = _mm_blendv_epi8(exit0, zero128, cmpim);

#if RDT
        tim1 = __rdtsc();
#endif
        
        // PR 17: pre-pass — build score vector sbt[j] for all band columns
        // once. Independent per-j, so vectorized loads/stores run at peak
        // throughput without carrying the DP core's dep chain. The gen_mat
        // branch is loop-invariant (once per row, perfectly predicted): the
        // symmetric default keeps the PR16 XOR LUT (1 xor + 1 shuffle); an
        // asymmetric matrix pays the heavier generic LUT only when needed.
        if (!gen_mat) {
            for (int jp = beg; jp < end; jp++) {
                __m128i s2 = _mm_load_si128((__m128i *)(seq2SoA + jp * SIMD_WIDTH8));
                __m128i sbt11;
                SBT_PREPASS8_XOR(s10, s2, sbt11, pmat128);
                _mm_store_si128((__m128i *)(sbt_buf + jp * SIMD_WIDTH8), sbt11);
            }
        } else if (fc.rank1) {
            __m128i rowfreed = _mm_cmpeq_epi8(s10, frref128);
            for (int jp = beg; jp < end; jp++) {
                __m128i s2 = _mm_load_si128((__m128i *)(seq2SoA + jp * SIMD_WIDTH8));
                __m128i sbt11;
                SBT_PREPASS8_RANK1(s10, s2, rowfreed, sbt11, pmat128, match128, frread128);
                _mm_store_si128((__m128i *)(sbt_buf + jp * SIMD_WIDTH8), sbt11);
            }
        } else {
            for (int jp = beg; jp < end; jp++) {
                __m128i s2 = _mm_load_si128((__m128i *)(seq2SoA + jp * SIMD_WIDTH8));
                __m128i sbt11;
                SBT_PREPASS8_AMAT(s10, s2, sbt11, amat128, w_ambig_128, three128);
                _mm_store_si128((__m128i *)(sbt_buf + jp * SIMD_WIDTH8), sbt11);
            }
        }

        j128 = _mm_set1_epi8(beg - i);   // diagonal offset of first band column
        for(j = beg; j < end; j++)
        {
            __m128i f11, f21, sbt11;
            h00 = _mm_load_si128((__m128i *)(H_h + j * SIMD_WIDTH8));
            f11 = _mm_load_si128((__m128i *)(F + j * SIMD_WIDTH8));
            sbt11 = _mm_load_si128((__m128i *)(sbt_buf + j * SIMD_WIDTH8));

            __m128i pj128 = j128;
            j128 = _mm_add_epi8(j128, one128);

            MAIN_CODE8_CORE(sbt11, h00, h11, e11, f11, f21, zero128,
                            e_ins128, oe_ins128,
                            e_del128, oe_del128);

            // Masked writing
            __m128i cmp1 = _mm_cmpgt_epi8(head128, pj128);
            __m128i cmp2 = _mm_cmpgt_epi8(pj128, tail128);
            cmp1 = _mm_or_si128(cmp1, cmp2);
            //__m128i cmpt = _mm_xor_si128(cmp1, ff128);
            h10 = _mm_blendv_epi8(h10, zero128, cmp1);
            f21 = _mm_blendv_epi8(f21, zero128, cmp1);
            
            // got this block out of MAIN_CODE
            __m128i bmaxRS = maxRS1;
            maxRS1 =_mm_max_epu8(maxRS1, h11);   // modif
            // UNSIGNED >: maxRS1 = max_epu8(.,h11) >= bmaxRS always, so
            // (maxRS1 >u bmaxRS) == (maxRS1 != bmaxRS). Signed cmpgt_epi8 here
            // mis-read scores >127 (long reads) as negative.
            __m128i cmpA = _mm_xor_si128(_mm_cmpeq_epi8(maxRS1, bmaxRS), ff128);
            __m128i cmpB =_mm_cmpeq_epi8(maxRS1, h11);                  
            cmpA = _mm_or_si128(cmpA, cmpB);
            cmp1 = _mm_cmpgt_epi8(j128, tail128);
            cmp1 = _mm_or_si128(cmp1, cmp2);
            cmpA = _mm_blendv_epi8(y1_128, j128, cmpA);
            y1_128 = _mm_blendv_epi8(cmpA, y1_128, cmp1);
            maxRS1 = _mm_blendv_epi8(maxRS1, bmaxRS, cmp1);                     

            _mm_store_si128((__m128i *)(F + j * SIMD_WIDTH8), f21);
            _mm_store_si128((__m128i *)(H_h + j * SIMD_WIDTH8), h10);

            h10 = h11;
                        
            // gscore query-end capture. Fire when (a) this column is the lane's
            // query-end column -- (col - i) == (qlen - i) iff col == qlen-1 -- and
            // qlen-i is in the representable in-band window (qlen_valid: qoff in
            // [-w, w+1]), so the lane's band-grown band has reached the query end
            // (== scalar's `end == qlen`), and (b) the lane is still alive at row
            // start (exit0). We capture the raw h11 (rebaselined byte) here; the
            // wide epilogue turns it into an absolute byte+B gscore that survives
            // re-baselining. No dependence on the trimmed tail128 (which the
            // off-diagonal query-end column is, by construction, often beyond).
            if (j >= minq)
            {
                // Fire exactly like the 16-bit tier (which is byte-identical to
                // scalar): the query-end column is reached (j128 == qlen_off, i.e.
                // col == qlen-1) AND the band-grown tail reaches the query end
                // (tail128 == qlen_off -- scalar's `end == qlen`; tail128 is
                // clamped to qlen_off at band-grow, so >= reduces to ==), AND the
                // lane is alive (exit0). qlen_valid guards the wrapped/out-of-band
                // qlen_off alias. Fires regardless of h11 (matches scalar's h1==0
                // firing); the captured value is finalized wide (byte+B) in the
                // epilogue so it survives re-baselining.
                //
                // gtle (= max_ie+1, the captured row) CONTRACT vs scalar: exact
                // when gscore > 0; may differ only in the gscore == 0 query-end
                // tail. The vector row loop ends at nrow = mlenw =
                // min(qlen+myband, tlen) (myband is e-dependent), whereas scalar
                // runs to its dynamic m==0 break and keeps updating max_ie on the
                // trailing end==qlen rows (all h1==0, gscore stays 0) that the
                // vector never reaches. This is harmless: bwa-mem consumes gtle
                // only under `gscore > 0` (see the `if (gscore <= 0 || ...) {...}
                // else {...gtle...}` guards in bwamem.cpp), so the gscore==0
                // divergence never reaches a SAM record. At default -E 1,
                // mlenw == tlen, so there is no divergence at all. Enforced by
                // test/integration/bandedswa_zdrop_eweight_test.cpp (gtle compared
                // iff gscore > 0).
                __m128i cmp = _mm_cmpeq_epi8(j128, qlen_off128);
                cmp = _mm_and_si128(cmp, _mm_cmpeq_epi8(tail128, qlen_off128));
                cmp = _mm_and_si128(cmp, qlen_valid128);
                cmp = _mm_and_si128(cmp, exit0);
                hqe128   = _mm_blendv_epi8(hqe128, h11, cmp);
                qfire128 = _mm_blendv_epi8(qfire128, ff128, cmp);
            }
        }
        __m128i cmp1 = _mm_cmpgt_epi8(head128, j128);
        __m128i cmp2 = _mm_cmpgt_epi8(j128, tail128);
        cmp1 = _mm_or_si128(cmp1, cmp2);
        h10 = _mm_blendv_epi8(h10, zero128, cmp1);
            
        _mm_store_si128((__m128i *)(H_h + j * SIMD_WIDTH8), h10);
        _mm_store_si128((__m128i *)(F + j * SIMD_WIDTH8), zero128);
        
        
        /* exit due to zero score by a row */
        __m128i bmaxScore128 = maxScore128;
        __m128i tmp = _mm_cmpeq_epi8(maxRS1, zero128);
#if defined(__ARM_NEON)
        int zero_row_ = (vmaxvq_u8(vreinterpretq_u8_m128i(maxRS1)) == 0);
#else
        int zero_row_ = (_mm_movemask_epi8(tmp) == 0xFFFF);
#endif
        if (zero_row_) {
            /* Finalize THIS row's query-end (gscore/gtle) capture BEFORE exiting.
             * The capture (hqe128/qfire128) is set in the inner loop above but is
             * folded into the wide gbest_abs/ierow only in the per-row epilogue,
             * which this break would otherwise skip. Scalar processes the row's
             * gscore and then hits its own m==0 break, so it records the query-end
             * row; the vector must too. Skipping it drops the gscore==0 query-end
             * tail row -> wrong gtle. Benign for symmetric scoring (gscore==0 gtle
             * unused), but consumed under the asymmetric --meth (OT/OB) matrix,
             * where it caused soft-clip / placement drift. Mirror of the epilogue's
             * gscore block (>= tie-break: latest query-end row wins, as scalar). */
            int8_t qf_a_[SIMD_WIDTH8] __attribute((aligned(16)));
            int8_t hq_a_[SIMD_WIDTH8] __attribute((aligned(16)));
            _mm_store_si128((__m128i *) qf_a_, qfire128);
            _mm_store_si128((__m128i *) hq_a_, hqe128);
            for (int g = 0; g < SIMD_WIDTH8 / 4; g++) {
                const int base = g * 4;
                __m128i Bg   = _mm_loadu_si128((const __m128i *)(B + base));
                __m128i hqeg = _mm_cvtepu8_epi32(_mm_cvtsi32_si128(*(const int32_t *)(hq_a_ + base)));
                __m128i qfg  = _mm_cvtepi8_epi32(_mm_cvtsi32_si128(*(const int32_t *)(qf_a_ + base)));
                __m128i gs_abs = _mm_add_epi32(hqeg, Bg);
                __m128i gba = _mm_loadu_si128((const __m128i *)(gbest_abs + base));
                __m128i ge  = _mm_xor_si128(_mm_cmpgt_epi32(gba, gs_abs), ff128); // gs_abs >= gba
                __m128i gmask = _mm_and_si128(qfg, ge);
                gba = _mm_blendv_epi8(gba, gs_abs, gmask);
                _mm_storeu_si128((__m128i *)(gbest_abs + base), gba);
                __m128i ierg = _mm_loadu_si128((const __m128i *)(ierow + base));
                ierg = _mm_blendv_epi8(ierg, _mm_set1_epi32(i + 1), gmask);
                _mm_storeu_si128((__m128i *)(ierow + base), ierg);
            }
            break;
        }

        // _mm_store_si128((__m128i *) temp, exit0);
        exit0 = _mm_blendv_epi8(exit0, zero128,  tmp);

        __m128i score128 = _mm_max_epu8(maxScore128, maxRS1);   // epi8 not present, modif
        maxScore128 = _mm_blendv_epi8(maxScore128, score128, exit0);

        // UNSIGNED >: maxScore128 (post-update, = max_epu8 of old & maxRS1 on
        // alive lanes, else unchanged) is >= bmaxScore128, so (>u) == (!=).
        // Signed cmpgt_epi8 mis-read scores >127.
        __m128i cmp = _mm_xor_si128(_mm_cmpeq_epi8(maxScore128, bmaxScore128), ff128);
        // y128 (best col) stays a diagonal offset captured in the best row's
        // frame; the best row itself moves to the wide xrow[] side channel.
        y128 = _mm_blendv_epi8(y128, y1_128, cmp);

        // max_off = max running diagonal-distance of the row-max from the main
        // diagonal: |y1col - (i+1)| = |y1_off - 1| in the offset frame.
        // y1_off - 1 is a SIGNED int8 (negative when row-max is left of the
        // sub-diagonal), so use _mm_abs_epi8 (SSSE3 / sse2neon on arm64).
        __m128i y1_minus1 = _mm_sub_epi8(y1_128, one128);
        tmp = _mm_abs_epi8(y1_minus1);                    // |y1_off - 1|
        __m128i bmax_off128 = max_off128;
        tmp = _mm_max_epu8(max_off128, tmp);  // modif
        max_off128 = _mm_blendv_epi8(bmax_off128, tmp, cmp);

        // Per-lane wide updates (O(rows)): best-score row (xrow), best-gscore
        // row (ierow), and the z-drop test — all done in wide scalars so row
        // distances that exceed int8 for long reads are handled exactly.
        {
            int8_t  cmp_a[SIMD_WIDTH8]      __attribute((aligned(16)));
            int8_t  y1_a[SIMD_WIDTH8]       __attribute((aligned(16)));
            int8_t  y_a[SIMD_WIDTH8]        __attribute((aligned(16)));
            int8_t  ms_a[SIMD_WIDTH8]       __attribute((aligned(16)));
            int8_t  rs_a[SIMD_WIDTH8]       __attribute((aligned(16)));
            int8_t  exit_a[SIMD_WIDTH8]     __attribute((aligned(16)));
            int8_t  qfire_a[SIMD_WIDTH8]    __attribute((aligned(16)));
            int8_t  hqe_a[SIMD_WIDTH8]      __attribute((aligned(16)));
            _mm_store_si128((__m128i *) cmp_a, cmp);
            _mm_store_si128((__m128i *) y1_a, y1_128);
            _mm_store_si128((__m128i *) y_a, y128);
            _mm_store_si128((__m128i *) ms_a, maxScore128);
            _mm_store_si128((__m128i *) rs_a, maxRS1);
            _mm_store_si128((__m128i *) exit_a, exit0);
            _mm_store_si128((__m128i *) qfire_a, qfire128);
            _mm_store_si128((__m128i *) hqe_a, hqe128);
            // VECTORIZED per-lane wide updates (replaces the scalar lane loop;
            // it was the dominant SW-kernel cost on AVX-512 per profiling). int32
            // is 4x int8, so SIMD_WIDTH8 lanes split into SIMD_WIDTH8/4 groups of
            // 4 int32. Each group does the same xrow / best_abs / gbest_abs+ierow /
            // z-drop work as the scalar loop in wide int32 SIMD; byte-identical
            // (same >, >= tie-breaks; same abs/z-drop arithmetic; B[] added per row).
            const __m128i vi   = _mm_set1_epi32(i);
            const __m128i vip1 = _mm_set1_epi32(i + 1);
            const __m128i vone = _mm_set1_epi32(1);
            const __m128i vzd  = _mm_set1_epi32(zdrop);
            const __m128i vedel = _mm_set1_epi32(this->e_del);
            const __m128i veins = _mm_set1_epi32(this->e_ins);
            __m128i die_g[SIMD_WIDTH8 / 4];
            for (int g = 0; g < SIMD_WIDTH8 / 4; g++) {
                const int base = g * 4;
                __m128i Bg   = _mm_loadu_si128((const __m128i *)(B + base));
                __m128i msg  = _mm_cvtepu8_epi32(_mm_cvtsi32_si128(*(const int32_t *)(ms_a    + base)));
                __m128i rsg  = _mm_cvtepu8_epi32(_mm_cvtsi32_si128(*(const int32_t *)(rs_a    + base)));
                __m128i hqeg = _mm_cvtepu8_epi32(_mm_cvtsi32_si128(*(const int32_t *)(hqe_a   + base)));
                __m128i cmpg = _mm_cvtepi8_epi32(_mm_cvtsi32_si128(*(const int32_t *)(cmp_a   + base)));
                __m128i exitg= _mm_cvtepi8_epi32(_mm_cvtsi32_si128(*(const int32_t *)(exit_a  + base)));
                __m128i qfg  = _mm_cvtepi8_epi32(_mm_cvtsi32_si128(*(const int32_t *)(qfire_a + base)));
                __m128i y1g  = _mm_cvtepi8_epi32(_mm_cvtsi32_si128(*(const int32_t *)(y1_a    + base)));
                __m128i yg   = _mm_cvtepi8_epi32(_mm_cvtsi32_si128(*(const int32_t *)(y_a     + base)));

                // (1) best-score row: xrow = cmp ? i+1 : xrow
                __m128i xrg = _mm_loadu_si128((const __m128i *)(xrow + base));
                xrg = _mm_blendv_epi8(xrg, vip1, cmpg);
                _mm_storeu_si128((__m128i *)(xrow + base), xrg);

                // (2) best_abs = max(best_abs, (uint8)ms + B)
                __m128i ms_abs = _mm_add_epi32(msg, Bg);
                __m128i bag = _mm_loadu_si128((const __m128i *)(best_abs + base));
                bag = _mm_max_epi32(bag, ms_abs);
                _mm_storeu_si128((__m128i *)(best_abs + base), bag);

                // (3) gscore: where qfire AND gs_abs >= gbest_abs, take gs_abs / i+1
                __m128i gs_abs = _mm_add_epi32(hqeg, Bg);
                __m128i gba = _mm_loadu_si128((const __m128i *)(gbest_abs + base));
                __m128i ge  = _mm_xor_si128(_mm_cmpgt_epi32(gba, gs_abs), ff128); // gs_abs >= gba
                __m128i gmask = _mm_and_si128(qfg, ge);
                gba = _mm_blendv_epi8(gba, gs_abs, gmask);
                _mm_storeu_si128((__m128i *)(gbest_abs + base), gba);
                __m128i ierg = _mm_loadu_si128((const __m128i *)(ierow + base));
                ierg = _mm_blendv_epi8(ierg, vip1, gmask);
                _mm_storeu_si128((__m128i *)(ierow + base), ierg);

                // (4) z-drop (alive lanes): dif = |((i+1)-xr) - (y1c-yc)|,
                //     drop = (uint8)ms - (uint8)rs; die where drop-dif > zdrop.
                __m128i y1c  = _mm_add_epi32(y1g, vi);
                __m128i yc   = _mm_add_epi32(yg, _mm_sub_epi32(xrg, vone));
                __m128i tmpi = _mm_sub_epi32(vip1, xrg);
                __m128i tmpj = _mm_sub_epi32(y1c, yc);
                // z-drop gap term weighted by gap-extend penalty, matching the
                // scalar reference (drift>0 -> deletion side *e_del, else *e_ins).
                __m128i zdelta = _mm_sub_epi32(tmpi, tmpj);
                __m128i zesel  = _mm_blendv_epi8(veins, vedel,
                                     _mm_cmpgt_epi32(zdelta, _mm_setzero_si128()));
                __m128i dif  = _mm_mullo_epi32(_mm_abs_epi32(zdelta), zesel);
                __m128i drop = _mm_sub_epi32(msg, rsg);
                __m128i die  = _mm_cmpgt_epi32(_mm_sub_epi32(drop, dif), vzd);
                die_g[g] = _mm_and_si128(die, exitg);
            }
            // pack the four int32 die masks back to a byte mask; clear dying lanes.
            __m128i die01 = _mm_packs_epi32(die_g[0], die_g[1]);
            __m128i die23 = _mm_packs_epi32(die_g[2], die_g[3]);
            __m128i die_bytes = _mm_packs_epi16(die01, die23);
            exit0 = _mm_andnot_si128(die_bytes, exit0);
        }

#if RDT
        prof[DP1][0] += __rdtsc() - tim1;
        tim1 = __rdtsc();
#endif
        
        /* Narrowing of the band */
        /* From beg */
        int l;
        for (l = beg; l < end; l++)
        {
            __m128i f128 = _mm_load_si128((__m128i *)(F + l * SIMD_WIDTH8));
            __m128i h128 = _mm_load_si128((__m128i *)(H_h + l * SIMD_WIDTH8));
#if defined(__ARM_NEON)
            if (vmaxvq_u8(vorrq_u8(vreinterpretq_u8_m128i(f128),
                                    vreinterpretq_u8_m128i(h128))) == 0)
                nbeg = l;
            else
                break;
#else
            __m128i tmp = _mm_or_si128(f128, h128);
            tmp = _mm_cmpeq_epi8(tmp, zero128);
            uint16_t val = _mm_movemask_epi8(tmp);
            if (val == 0xFFFF) nbeg = l;
            else
                break;
#endif
        }

        /* From end */
        for (l = end; l >= beg; l--)
        {
            __m128i f128 = _mm_load_si128((__m128i *)(F + l * SIMD_WIDTH8));
            __m128i h128 = _mm_load_si128((__m128i *)(H_h + l * SIMD_WIDTH8));
#if defined(__ARM_NEON)
            if (vmaxvq_u8(vorrq_u8(vreinterpretq_u8_m128i(f128),
                                    vreinterpretq_u8_m128i(h128))) != 0)
                break;
#else
            __m128i tmp = _mm_or_si128(f128, h128);
            tmp = _mm_cmpeq_epi8(tmp, zero128);
            uint16_t val = _mm_movemask_epi8(tmp);
            if (val != 0xFFFF)
                break;
#endif
        }
        // int pnend =nend;
        nend = l + 2 < ncol? l + 2: ncol;
        __m128i tmpb = ff128;

        __m128i exit1 = _mm_xor_si128(exit0, ff128);
        __m128i l128 = _mm_set1_epi8(beg - i);   // diagonal offset of column beg

        for (l = beg; l < end; l++)
        {
            __m128i f128 = _mm_load_si128((__m128i *)(F + l * SIMD_WIDTH8));
            __m128i h128 = _mm_load_si128((__m128i *)(H_h + l * SIMD_WIDTH8));

            __m128i tmp = _mm_or_si128(f128, h128);
            //tmp = _mm_or_si128(tmp, _mm_xor_si128(exit0, ff128));
            tmp = _mm_or_si128(tmp, exit1);
            tmp = _mm_cmpeq_epi8(tmp, zero128);
#if defined(__ARM_NEON)
            if (vmaxvq_u8(vreinterpretq_u8_m128i(tmp)) == 0) {
                break;
            }
#else
            uint32_t val = _mm_movemask_epi8(tmp);
            if (val == 0x00) {
                break;
            }
#endif
            tmp = _mm_and_si128(tmp,tmpb);
            l128 = _mm_add_epi8(l128, one128);
            // NEW
            head128 = _mm_blendv_epi8(head128, l128, tmp);

            tmpb = tmp;
        }

        __m128i  index128 = tail128;
        tmpb = ff128;

        l128 = _mm_set1_epi8(end - i);   // diagonal offset of column end
        for (l = end; l >= beg; l--)
        {
            __m128i f128 = _mm_load_si128((__m128i *)(F + l * SIMD_WIDTH8));
            __m128i h128 = _mm_load_si128((__m128i *)(H_h + l * SIMD_WIDTH8));

            __m128i tmp = _mm_or_si128(f128, h128);
            tmp = _mm_or_si128(tmp, exit1);
            tmp = _mm_cmpeq_epi8(tmp, zero128);
#if defined(__ARM_NEON)
            if (vmaxvq_u8(vreinterpretq_u8_m128i(tmp)) == 0) {
                break;
            }
#else
            uint32_t val = _mm_movemask_epi8(tmp);
            if (val == 0x00)  {
                break;
            }
#endif
            tmp = _mm_and_si128(tmp,tmpb);
            l128 = _mm_sub_epi8(l128, one128);
            // NEW
            index128 = _mm_blendv_epi8(index128, l128, tmp);

            tmpb = tmp;
        }
        index128 = _mm_add_epi8(index128, two128);
        // signed min in the offset frame against qlen-i
        tail128 = _mm_min_epi8(index128, qlen_off128);

        // Frame shift for the next row: i advances by 1, so the same absolute
        // band edge has its diagonal offset (col - i) decremented by 1. Keep
        // head/tail tracking the same columns as the frame moves.
        head128 = _mm_sub_epi8(head128, one128);
        tail128 = _mm_sub_epi8(tail128, one128);

        // --- RE-BASELINE EVENT (per-lane score floor) ---
        // After this row's score state is final, lower any lane whose row-max
        // (maxRS1, rebaselined) has climbed to REBASE_HI back into range. We
        // raise B[l] by delta = maxRS1[l] - REBASE_KEEP and subtract that delta
        // (saturating to 0) from every carried score buffer/register for that
        // lane. Done AFTER band narrowing so this row's trim matches the 16-bit
        // reference; the next row reads the lowered values consistently.
        {
            /* Fast path: one vector test for "any lane >= REBASE_HI?". Under the
             * routing envelope no admitted pair ever reaches REBASE_HI (the gate
             * keeps the max attainable score below it), so this is ~always false and
             * the per-lane scalar scan + maxRS1 store below are skipped each row.
             * Unsigned >=: subs_epu8(REBASE_HI, maxRS1)==0 iff maxRS1 >= REBASE_HI. */
            __m128i hi128 = _mm_set1_epi8((int8_t) REBASE_HI);
            if (_mm_movemask_epi8(_mm_cmpeq_epi8(_mm_subs_epu8(hi128, maxRS1), zero128))) {
            uint8_t rs_b[SIMD_WIDTH8] __attribute((aligned(16)));
            _mm_store_si128((__m128i *) rs_b, maxRS1);
            int8_t  delta_a[SIMD_WIDTH8] __attribute((aligned(16)));
            int any = 0;
            for (int l = 0; l < SIMD_WIDTH8; l++) {
                int d = 0;
                if ((int)rs_b[l] >= REBASE_HI) {
                    d = (int)rs_b[l] - REBASE_KEEP;   // keep low REBASE_KEEP, in range
                    B[l] += d;
                    any = 1;
                }
                delta_a[l] = (int8_t) d;              // 0..(REBASE_HI-REBASE_KEEP) small
            }
            if (any) {
                __m128i delta128 = _mm_load_si128((__m128i *) delta_a);
                // Subtract delta from every carried score buffer over the active
                // band columns. The range [lo,hi] is the union of the current and
                // next band windows; clamping to [0, ncol-1] keeps the loop within
                // the valid column range (a safe superset covering every column the
                // next row will read).
                int lo = nbeg < beg ? nbeg : beg;
                int hi = nend > end ? nend : end;
                if (lo < 0) lo = 0;
                if (hi + 1 > ncol) hi = ncol - 1;
                for (int l = lo; l <= hi; l++) {
                    __m128i h128 = _mm_load_si128((__m128i *)(H_h + l * SIMD_WIDTH8));
                    __m128i f128 = _mm_load_si128((__m128i *)(F   + l * SIMD_WIDTH8));
                    h128 = _mm_subs_epu8(h128, delta128);
                    f128 = _mm_subs_epu8(f128, delta128);
                    _mm_store_si128((__m128i *)(H_h + l * SIMD_WIDTH8), h128);
                    _mm_store_si128((__m128i *)(F   + l * SIMD_WIDTH8), f128);
                }
                // The vertical seed H_v is read only while the band still touches
                // column 0 (beg==0). Lower the rows that may still be consumed.
                if (beg == 0) {
                    for (int r = i + 1; r < nrow; r++) {
                        __m128i v128 = _mm_load_si128((__m128i *)(H_v + r * SIMD_WIDTH8));
                        v128 = _mm_subs_epu8(v128, delta128);
                        _mm_store_si128((__m128i *)(H_v + r * SIMD_WIDTH8), v128);
                    }
                }
                // Carried score maxima (registers) live in the same rebaselined
                // frame and must be lowered too so the next row's z-drop (ms-rs)
                // and the absolute-frame add-back stay consistent. maxScore128 is
                // the running max, so it never saturates to 0 (maxScore128 >=
                // maxRS1 for the triggering row, and maxRS1 >= REBASE_HI >
                // REBASE_KEEP > delta, so it stays positive).
                maxScore128 = _mm_subs_epu8(maxScore128, delta128);
                // gscore no longer carries a rebaselined byte register: the
                // query-end cell H is captured per row and accumulated wide
                // (gbest_abs, in absolute units) in the epilogue, so it is immune
                // to this re-baseline subtract.
            }
            }   /* end "any lane >= REBASE_HI" fast-path guard */
        }

#if RDT
        prof[DP2][0] += __rdtsc() - tim1;
#endif
    }
   
#if RDT
    prof[DP][0] += __rdtsc() - tim;
#endif
    
    // Scores are reconstructed in ABSOLUTE units from the per-lane wide
    // best_abs/gbest_abs side channels (= rebaselined byte + B[l], captured each
    // row before any later re-baseline could saturate the byte away). Positions
    // are reconstructed wide from the diagonal-offset lanes plus the per-lane
    // best-row side channels.
    int8_t maxj[SIMD_WIDTH8]  __attribute((aligned(64)));
    _mm_store_si128((__m128i *) maxj, y128);   // best col as diagonal offset

    int8_t max_off_ar[SIMD_WIDTH8]  __attribute((aligned(64)));
    _mm_store_si128((__m128i *) max_off_ar, max_off128);

    for(i = 0; i < SIMD_WIDTH8; i++)
    {
        p[i].score   = best_abs[i];                               // absolute score
        p[i].tle     = xrow[i];                                   // best row (target end)
        // qle reconstruction: maxj[l] = col - (xrow-1), so col = maxj[l] + (xrow-1).
        // Guard the unset case (xrow==0 means no update occurred).
        p[i].qle     = (xrow[i] == 0) ? 0 : ((int) maxj[i] + xrow[i] - 1);
        p[i].max_off = (uint8_t) max_off_ar[i];
        p[i].gscore  = gbest_abs[i];                              // absolute gscore (-1 if unset)
        p[i].gtle    = ierow[i];                                  // best gscore row
    }

    return;
}

#endif

/* Per-tier factory function. Compiled into each KERNEL_VARIANT build of
 * this TU; the symbol is mangled by kernel_dispatch.h to
 * make_bsw_kernel_<tier>. On arm64 (no KERNEL_VARIANT) this is the
 * unmangled make_bsw_kernel.
 *
 * Returns IBandedPairWiseSW* (not unique_ptr) because extern "C" disallows
 * non-POD return types. The dispatcher in simd_dispatch.cpp wraps it into
 * unique_ptr at the call site. */
extern "C" IBandedPairWiseSW *make_bsw_kernel(
    int o_del, int e_del, int o_ins, int e_ins, int zdrop,
    int end_bonus, const int8_t *mat,
    int8_t w_match, int8_t w_mismatch, int numThreads)
{
    return new BandedPairWiseSW(o_del, e_del, o_ins, e_ins, zdrop, end_bonus,
                                mat, w_match, w_mismatch, numThreads);
}
