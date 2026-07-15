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
         Heng Li <hli@jimmy.harvard.edu>
*****************************************************************************************/

#include "bwamem.h"
#include "FMI_search.h"
#include "smem_dedup.h"
#include "bam_writer.h"
#include "meth_bam.h"
#include "u8vec_scratch.h"
#include "pdqsort_wrap.h"
#include <inttypes.h>      /* PRId64 for int64_t fprintf format strings */

#include "sam_encode.h"
#include "stage_prof.h"
#include "seed_order.h"
/* Not including <htslib/sam.h>: its kstring.h shares the KSTRING_H guard
 * with bwa-mem3's. Opaque bam1_t wrappers live in bam_writer.h. */

/* RAII bracket for the SAM/BAM build (stage_prof --profile). Accumulates the
 * calling (compute) thread's CPU time across all of mem_aln2sam's return paths
 * into the per-thread encode accumulator; harvested by kt_for. Zero cost when off. */
namespace { struct SpEncodeScope {
    double t0; bool on;
    SpEncodeScope() : t0(0.0), on(sp_enabled()) { if (on) t0 = sp_thread_cpu(); }
    ~SpEncodeScope() { if (on) sp_encode_add(sp_thread_cpu() - t0); }
}; }

/* D3 (--meth, PR-5): original-reference pac for XM:Z, set in fastmap.cpp.
 * The meth output path is gated on opt->meth_mode (the chrom map that used to
 * gate it was retired with the f/r output layer). */
const uint8_t *g_meth_orig_pac = NULL;

//----------------
extern uint64_t tprof[LIM_R][LIM_C];
//----------------
#include "kbtree.h"

#define chain_cmp(a, b) (((b).pos < (a).pos) - ((a).pos < (b).pos))
KBTREE_INIT(chn, mem_chain_t, chain_cmp)

#define intv_lt1(a, b) ((((uint64_t)(a).m) <<32 | ((uint64_t)(a).n)) < (((uint64_t)(b).m) <<32 | ((uint64_t)(b).n)))  // trial
KSORT_INIT(mem_intv1, SMEM, intv_lt1)  // debug

#define max_(x, y) ((x)>(y)?(x):(y))
#define min_(x, y) ((x)>(y)?(y):(x))

#define MAX_BAND_TRY  4

/* cap initial band-width at this value. The retry loop
 * doubles w each iteration (up to MAX_BAND_TRY-1 iters), so pairs whose
 * alignment needs a wider band are caught by the retry. Setting this below
 * opt->w forces the kernel to start tight, accept tight-fit pairs early
 * (via sp->max_off heuristic + sp->tight_band), and only expand on demand.
 * BUCKET_MAX_INIT_W=8 with MAX_BAND_TRY=4 gives w-sequence 8, 16, 32, 64. */
#ifndef BUCKET_MAX_INIT_W
#define BUCKET_MAX_INIT_W 8
#endif

/* ------------------------------------------------------------------------
 * 8-bit (16-lane) SW safe-envelope gate.
 *
 * The recovered 8-bit kernel `smithWaterman128_8` (and the AVX2/AVX-512 ports)
 * runs an UNSIGNED [0,255] score recurrence and is byte-exact for long reads
 * ONLY inside the envelope below. The decisive constraint is that the per-lane
 * score re-baselining must NEVER fire for an admitted pair: its saturating
 * subtract is not generally lossless (it can zero a still-positive off-diagonal
 * cell, which the recurrence then misreads as the h00==0 local-restart sentinel
 * and truncates a valid alignment — z-drop is a row-level early-exit, not a
 * per-cell guarantee). Held inert, the kernel is a plain exact unsigned SW.
 *
 *   - len1, len2 < MAX_SEQ_LEN8  : slab width + the wide per-row position
 *     side channel are sized MAX_SEQ_LEN8; rows/cols at or beyond it would
 *     run off the slab.
 *   - len1 >= len2 (target >= query) : the gscore (query-end) capture is
 *     byte-identical to scalar only when the query-end column lies on/left of
 *     the main diagonal. See smithWaterman128_8 gscore capture.
 *   - w <= BSW8_MAX_W (= 127)     : every band/position quantity is encoded
 *     as a diagonal offset d = j - i in [-w, +w], which fits signed int8
 *     only for w <= 127. The default opt->w = 100 qualifies; wide-band
 *     retries (w doubling to 200/400/800) do NOT and fall back to 16-bit.
 *   - zdrop + maxStep <= 253      : keeps the re-baseline window non-empty,
 *     REBASE_HI = 255 - maxStep > REBASE_KEEP = zdrop + 1. The DP body AND the
 *     h0-prefix column/row seed (smithWaterman*_8 setup) are both unsigned-
 *     saturating [0,255], so the seeded byte h0' = min(h0, REBASE_KEEP) only
 *     needs to fit a uint8. maxStep is the largest per-step score increment,
 *     max(w_match, w_ambig, 1) = max(opt->a, 1).
 *   - h0 <= zdrop + 1             : keeps the initial floor B0 = max(0, h0 -
 *     REBASE_KEEP) == 0, so the seed score does not itself force a re-baseline.
 *   - h0 + min(len1,len2)*a < 255 - maxStep : the MAX ATTAINABLE score (seed
 *     plus an all-match diagonal over the shorter sequence) stays below REBASE_HI,
 *     so the running row max never reaches it and re-baseline never fires.
 *
 * NOTE: minval (= h0 + min(len)*a) is still used by sortPairsLenExt as the
 * counting-sort bin index (hist[minval]) and must stay < MAX_SEQ_LEN8 there for
 * the bin to be in range; that is a histogram-sizing constraint, NOT the tier
 * decision. Pairs failing this envelope fall through to the existing 16-bit
 * (then scalar) buckets exactly as before. */
#define BSW8_MAX_W 127                 /* max band: diagonal offset d=j-i must fit signed int8 */
#define BSW8_MAX_ZDROP_STEP 253        /* keep the re-baseline window non-empty:
                                          REBASE_HI=255-max_step > REBASE_KEEP=zdrop+1, i.e.
                                          zdrop+max_step <= 253. The DP body AND the h0-prefix
                                          column/row seed in the 8-bit wrappers are both
                                          unsigned-saturating [0,255], so the seed byte
                                          min(h0,REBASE_KEEP) just needs to fit a uint8.
                                          See smithWaterman128_8 re-baseline + h0-prefix seed. */

static inline int bsw8_envelope_ok(int len1, int len2, int w,
                                   int score_a, int zdrop, int h0)
{
    /* Dev A/B hook: BWAMEM3_DISABLE_BSW8=1 forces every pair off the 8-bit path
     * onto the 16-bit (then scalar) buckets. Used to validate that the 8-bit
     * kernel is byte-identical to the 16-bit reference in the full pipeline
     * (diff the SAM with vs without). Off by default. The function-local static
     * with a non-trivial initializer gets C++11 thread-safe one-time init, so the
     * getenv runs exactly once even under the OpenMP/pthread worker fan-out (a
     * plain `static int x = -1; if (x<0) x = ...;` would be a data race). */
    static const int disable = []() {
        const char *e = getenv("BWAMEM3_DISABLE_BSW8");
        return (e && atoi(e)) ? 1 : 0;
    }();
    if (disable) return 0;
    /* 64-bit math: len < MAX_SEQ_LEN8 and a small score_a keep these tiny in
     * practice, but compute the bounds in int64_t so a pathological scoring
     * parameter cannot wrap an int and admit a pair the guards would reject. */
    int64_t max_step = score_a > 1 ? (int64_t)score_a : 1;   /* max(w_match, w_ambig, 1) */
    int64_t shorter  = len1 < len2 ? (int64_t)len1 : (int64_t)len2;
    /* Max attainable SW score: seed h0 plus an all-match diagonal over the
     * shorter sequence. Must stay strictly below REBASE_HI = 255 - max_step so
     * the running row max never triggers a (lossy) re-baseline; combined with
     * h0 <= zdrop+1 (B0 == 0) the kernel runs as an exact unsigned [0,255] SW. */
    int64_t max_score = (int64_t)h0 + shorter * (int64_t)score_a;
    return len1 < MAX_SEQ_LEN8 && len2 < MAX_SEQ_LEN8 &&
           /* target (rows, len1) >= query (cols, len2). The re-baseline 8-bit
            * kernel's gscore (query-end score) is byte-identical to scalar ONLY
            * when the query-end column (len2-1) lies on or left of the main
            * diagonal, i.e. len1 >= len2. When len2 > len1 the query end is
            * off-diagonal to the right, where the re-baseline saturating-subtract
            * zeroes those low-scoring cells (they sit >zdrop below the row max),
            * the band trims them, and gscore is lost. Such pairs (rare: only when
            * the reference window is shorter than the read segment, e.g. at contig
            * edges) route to the 16-bit tier, which has no re-baseline and is
            * byte-identical for them. See smithWaterman128_8 gscore capture. */
           len1 >= len2 &&
           w <= BSW8_MAX_W &&
           zdrop + max_step <= BSW8_MAX_ZDROP_STEP &&
           h0 <= zdrop + 1 &&
           max_score < 255 - max_step;
}

            int tcnt = 0;
/********************
 * Filtering chains *
 ********************/

#define chn_beg(ch) ((ch).seeds->qbeg)
#define chn_end(ch) ((ch).seeds[(ch).n-1].qbeg + (ch).seeds[(ch).n-1].len)

#define flt_lt(a, b) ((a).w > (b).w)
KSORT_INIT(mem_flt, mem_chain_t, flt_lt)

/* Chain-geometry adaptive band (--adaptive-band). Start the banded-SW tight at
 * opt->band_start on the 16-bit and scalar (long-extension) tiers; carry each
 * pair's chain_band (= its chain's seed diagonal spread, capped at opt->w = the
 * implied inter-seed indel). ACCEPT_PAIR gates the score/max_off "converged"
 * accepts on w >= chain_band, so a pair whose chain implies a larger indel keeps
 * retrying to that band while a diagonal-hugging pair (chain_band=0) accepts in
 * one tight pass. The 8-bit tier stays at opt->w (its int8 diagonal encoding caps
 * at 127, and short extensions are sub-band so gain nothing). tight_band keeps its
 * ungapped-estimate accept-early role. band_start<=0 (default): INIT_W==opt->w and
 * ACCEPT_PAIR reduces to the original condition -> byte-identical. */
#define INIT_W(w) (opt->band_start > 0 ? min_(opt->band_start, (w)) : (w))
#define ACCEPT_PAIR(sc,pv,mo,w,tb,cb,i) \
    ((i)+1==MAX_BAND_TRY || ((tb)>0 && (w)>=(tb)) || \
     (((sc)==(pv) || (mo) < ((w)>>1)+((w)>>2)) && (opt->band_start <= 0 || (w) >= (cb))))
//------------------------------------------------------------------
// Alignment: Construct the alignment from a chain *

static inline int cal_max_gap(const mem_opt_t *opt, int qlen)
{
    int l_del = (int)((double)(qlen * opt->a - opt->o_del) / opt->e_del + 1.);
    int l_ins = (int)((double)(qlen * opt->a - opt->o_ins) / opt->e_ins + 1.);
    //int l_del = (int)((double)(qlen * opt->a - opt->o_del) + 1.);
    //int l_ins = (int)((double)(qlen * opt->a - opt->o_ins) + 1.);

    int l = l_del > l_ins? l_del : l_ins;
    l = l > 1? l : 1;
    return l < opt->w<<1? l : opt->w<<1;
}

//------------------------------------------------------------------
// SMEMs
static smem_aux_t *smem_aux_init()
{
    smem_aux_t *a;
    if ((a = (smem_aux_t *) calloc(BATCH_SIZE, sizeof(smem_aux_t))) == NULL) { fprintf(stderr, "ERROR: out of memory %s\n", __func__); exit(EXIT_FAILURE); }
    for (int i=0; i<BATCH_SIZE; i++)
    {
        a[i].tmpv[0] = (bwtintv_v *) calloc(1, sizeof(bwtintv_v));
        a[i].tmpv[1] = (bwtintv_v *) calloc(1, sizeof(bwtintv_v));
        if (!a[i].tmpv[0] || !a[i].tmpv[1]) { fprintf(stderr, "ERROR: out of memory %s\n", __func__); exit(EXIT_FAILURE); }
    }
    return a;
}

static void smem_aux_destroy(smem_aux_t *a)
{
    for (int i=0; i<BATCH_SIZE; i++)
    {
        free(a[i].tmpv[0]->a);
        free(a[i].tmpv[0]);
        free(a[i].tmpv[1]->a);
        free(a[i].tmpv[1]);
        free(a[i].mem.a);
        free(a[i].mem1.a);
    }
    free(a);
}

mem_opt_t *mem_opt_init()
{
    mem_opt_t *o;
    if ((o = (mem_opt_t *) calloc(1, sizeof(mem_opt_t))) == NULL)  { fprintf(stderr, "ERROR: out of memory\n"); exit(1); }
    o->flag = 0;
    o->a = 1; o->b = 4;
    o->o_del = o->o_ins = 6;
    o->e_del = o->e_ins = 1;
    o->w = 100;
    o->T = 30;
    o->zdrop = 100;
    o->pen_unpaired = 17;
    o->pen_clip5 = o->pen_clip3 = 5;

    o->max_mem_intv = 20;

    o->min_seed_len = 19;
    o->min_ext_len = 0;   // off by default -> byte-identical to baseline
    o->max_extend_chains = 0;   // off by default; opt-in speed lever (--max-extend-chains / --fast)
    o->mate_concordant_window = 0;   // off by default; opt-in (--extend-mate-concordant / --fast --meth)
    o->est_insert_high = 0;          // runtime state; set from data (mem_pestat) or -I during the run
    o->seed_emit_order = SEED_ORDER_OFF;  // byte-identical default
    o->smem_dedup  = 0;   // off by default -> byte-identical to baseline; opt-in via --smem-dedup
    o->skip_contained_ext = 0;   // off by default; opt-in via --skip-contained-ext (byte-identical)
    o->band_start  = 0;   // off by default (adaptive chain-geometry band); opt-in via --adaptive-band
    o->split_width = 10;
    o->max_occ = 500;
    o->max_chain_gap = 10000;
    o->max_ins = 10000;
    o->mask_level = 0.50;
    o->drop_ratio = 0.50;
    o->XA_drop_ratio = 0.80;
    o->split_factor = 1.5;
    o->chunk_size = 10000000;
    o->n_threads = 1;
    o->max_XA_hits = 5;
    o->max_XA_hits_alt = 200;
    o->max_matesw = 50;
    o->mask_level_redun = 0.95;
    o->min_chain_weight = 0;
    o->max_chain_extend = 1<<30;
    o->mapQ_coef_len = 50; o->mapQ_coef_fac = log(o->mapQ_coef_len);
    o->meth_scoring = MEM_METH_SCORING_COLLAPSED;  /* --meth default: bwameth-compatible */
    bwa_fill_scmat(o->a, o->b, o->mat);
    mem_opt_fill_meth_mat(o);
    return o;
}

/* D3 (--meth, PR-4): (re)derive the per-hypothesis ASYMMETRIC matrices from the
 * symmetric `mat` + match score `a`. Target-major mat[ref*5+read], ACGT order
 * A,C,G,T,N. Free exactly ONE off-diagonal cell to a MATCH (+a); leave everything
 * else (incl. the mirror cell, a real variant) at the symmetric value:
 *   OT: free mat[C][T] = mat[1*5+3]  (ref C, read T = unmethylated C→T)
 *   OB: free mat[G][A] = mat[2*5+0]  (ref G, read A = bottom-strand G→A)
 * GENOMIC (opt->meth_scoring): only that one cell is freed ⇒ the mirror stays a
 * real mismatch and the SIMD kernel takes its rank-1 fast path. COLLAPSED: the
 * mirror cell is ALSO freed (two cells) so C/T and G/A are interchangeable
 * (bwameth-compatible), via bandedSWA's general-matrix path. Selected per chain
 * via mem_opt_meth_mat(); valid only under --meth. MUST be called after EVERY
 * rebuild of opt->mat (mem_opt_init AND after main_mem parses -A/-B/-x and
 * re-runs bwa_fill_scmat) — otherwise the meth matrices keep the default
 * match/mismatch and silently ignore the user's scoring options. */
void mem_opt_fill_meth_mat(mem_opt_t *o) {
    memcpy(o->mat_ot, o->mat, sizeof(o->mat));
    memcpy(o->mat_ob, o->mat, sizeof(o->mat));
    o->mat_ot[1 * 5 + 3] = o->a;   /* OT: ref C / read T → match (C→T conversion)   */
    o->mat_ob[2 * 5 + 0] = o->a;   /* OB: ref G / read A → match (G→A conversion)   */
    if (o->meth_scoring == MEM_METH_SCORING_COLLAPSED) {
        /* bwameth-compatible: free the MIRROR cell too so C/T (and G/A) are
         * mutually interchangeable (collapsed 3-letter space). Two freed cells ⇒
         * bandedSWA general path, not the rank-1 fast path. Loses variant vs.
         * conversion discrimination; reproduces bwameth placement. */
        o->mat_ot[3 * 5 + 1] = o->a;   /* OT: ref T / read C → match */
        o->mat_ob[0 * 5 + 2] = o->a;   /* OB: ref A / read G → match */
    }
}

/******************************
 * De-overlap single-end hits *
 ******************************/

/* mem_ars2: primary sort by reference END position (`re`). Used by
 * mem_sort_dedup_patch as the first ordering before its order-sensitive
 * dedup/patch loop. The tie-break keys (rb, score-desc, qb) make the
 * order deterministic at equal `re`; without them, the dedup outcome
 * depended on whichever ordering the underlying sort algorithm happened
 * to produce (klib introsort is unstable, pdqsort is unstable
 * differently, radix is stable) — different orderings produced
 * different surviving alnreg sets, propagating to `sub_n` and MAPQ.
 * With the full key, every sort algorithm produces the same surviving
 * set and the SAM is bit-equal regardless of which sort runs. */
#define alnreg_slt2(a, b) \
    ((a).re < (b).re || ((a).re == (b).re && \
     ((a).rb < (b).rb || ((a).rb == (b).rb && \
      ((a).score > (b).score || ((a).score == (b).score && \
       (a).qb < (b).qb))))))
KSORT_INIT(mem_ars2, mem_alnreg_t, alnreg_slt2)
PDQSORT_INIT(mem_ars2, mem_alnreg_t, alnreg_slt2)

#define alnreg_slt(a, b) ((a).score > (b).score || ((a).score == (b).score && ((a).rb < (b).rb || ((a).rb == (b).rb && (a).qb < (b).qb))))
KSORT_INIT(mem_ars, mem_alnreg_t, alnreg_slt)
PDQSORT_INIT(mem_ars, mem_alnreg_t, alnreg_slt)

#define alnreg_hlt(a, b)  ((a).score > (b).score || ((a).score == (b).score && ((a).is_alt < (b).is_alt || ((a).is_alt == (b).is_alt && (a).hash < (b).hash))))
KSORT_INIT(mem_ars_hash, mem_alnreg_t, alnreg_hlt)

#define alnreg_hlt2(a, b) ((a).is_alt < (b).is_alt || ((a).is_alt == (b).is_alt && ((a).score > (b).score || ((a).score == (b).score && (a).hash < (b).hash))))
KSORT_INIT(mem_ars_hash2, mem_alnreg_t, alnreg_hlt2)

#if MATE_SORT
void sort_alnreg_re(int n, mem_alnreg_t* a) {
    pdqsort_mem_ars2(n, a);
}

void sort_alnreg_score(int n, mem_alnreg_t* a) {
    pdqsort_mem_ars(n, a);
}

#endif

#define PATCH_MAX_R_BW 0.05f
#define PATCH_MIN_SC_RATIO 0.90f

int mem_patch_reg(const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac,
                  uint8_t *query, const mem_alnreg_t *a, const mem_alnreg_t *b,
                  int *_w, const int8_t *mat)
{
    /* D3 (--meth, PR-4): `mat` is the per-hypothesis substitution matrix and
     * `query` the original read bases; the caller passes opt->mat + projected
     * read outside --meth so behavior is unchanged there. */
    int w, score, q_s, r_s;
    double r;
    if (bns == 0 || pac == 0 || query == 0) return 0;
    assert(a->rid == b->rid && a->rb <= b->rb);

    if (a->rb < bns->l_pac && b->rb >= bns->l_pac) return 0; // on different strands
    if (a->qb >= b->qb || a->qe >= b->qe || a->re >= b->re) return 0; // not colinear
    w = (a->re - b->rb) - (a->qe - b->qb); // required bandwidth
    w = w > 0? w : -w; // l = abs(l)
    r = (double)(a->re - b->rb) / (b->re - a->rb) - (double)(a->qe - b->qb) / (b->qe - a->qb); // relative bandwidth
    r = r > 0.? r : -r; // r = fabs(r)

    if (bwa_verbose >= 4)
        fprintf(stderr, "* potential hit merge between [%d,%d)<=>[%ld,%ld) and "
               "[%d,%d)<=>[%ld,%ld), @ %s; w=%d, r=%.4g\n",
               a->qb, a->qe, (long)a->rb, (long)a->re, b->qb, b->qe,
               (long)b->rb, (long)b->re, bns->anns[a->rid].name, w, r);

    if (a->re < b->rb || a->qe < b->qb) // no overlap on query or on ref
    {
        if (w > opt->w<<1 || r >= PATCH_MAX_R_BW) return 0; // the bandwidth or the relative bandwidth is too large
    } else if (w > opt->w<<2 || r >= PATCH_MAX_R_BW*2) return 0; // more permissive if overlapping on both ref and query

    // global alignment
    w += a->w + b->w;
    w = w < opt->w<<2? w : opt->w<<2;

    if (bwa_verbose >= 4)
        fprintf(stderr, "* test potential hit merge with global alignment; w=%d\n", w);

    bwa_gen_cigar2(mat, opt->o_del, opt->e_del, opt->o_ins, opt->e_ins, w,
                   bns->l_pac, pac, b->qe - a->qb, query + a->qb, a->rb, b->re,
                   &score, 0, 0);

    q_s = (int)((double)(b->qe - a->qb) / ((b->qe - b->qb) + (a->qe - a->qb)) *
                (b->score + a->score) + .499); // predicted score from query

    r_s = (int)((double)(b->re - a->rb) / ((b->re - b->rb) + (a->re - a->rb)) *
                (b->score + a->score) + .499); // predicted score from ref

    if (bwa_verbose >= 4)
        fprintf(stderr, "* score=%d;(%d,%d)\n", score, q_s, r_s);

    if ((double)score / (q_s > r_s? q_s : r_s) < PATCH_MIN_SC_RATIO) return 0;
    *_w = w;
    return score;
}
/*********************************
 * Test if a seed is good enough *
 *********************************/

#define MEM_SHORT_EXT 50
#define MEM_SHORT_LEN 200

#define MEM_HSP_COEF 1.1f
#define MEM_MINSC_COEF 5.5f
#define MEM_SEEDSW_COEF 0.05f

#if MATE_SORT
int mem_dedup_patch(const mem_opt_t *opt, const bntseq_t *bns,
                    const uint8_t *pac, uint8_t *query, int n,
                    mem_alnreg_t *a)
{
    int m, i, j;
    if (n <= 1) return n;
    for (i = 0; i < n; ++i) a[i].n_comp = 1;

    for (i = 1; i < n; ++i)
    {
        mem_alnreg_t *p = &a[i];
        if (p->rid != a[i-1].rid || p->rb >= a[i-1].re + opt->max_chain_gap)
            continue; // then no need to go into the loop below

        for (j = i - 1; j >= 0 && p->rid == a[j].rid && p->rb < a[j].re + opt->max_chain_gap; --j) {
            mem_alnreg_t *q = &a[j];
            int64_t or_, oq, mr, mq;
            int score, w;
            if (q->qe == q->qb) continue; // a[j] has been excluded
            or_ = q->re - p->rb; // overlap length on the reference
            oq = q->qb < p->qb? q->qe - p->qb : p->qe - q->qb; // overlap length on the query
            mr = q->re - q->rb < p->re - p->rb? q->re - q->rb : p->re - p->rb; // min ref len in alignment
            mq = q->qe - q->qb < p->qe - p->qb? q->qe - q->qb : p->qe - p->qb; // min qry len in alignment
            if (or_ > opt->mask_level_redun * mr && oq > opt->mask_level_redun * mq) { // one of the hits is redundant
                if (p->score < q->score)
                {
                    p->qe = p->qb;
                    break;
                }
                else q->qe = q->qb;
            }
            else if (q->rb < p->rb && (score = mem_patch_reg(opt, bns, pac, query, q, p, &w, opt->mat)) > 0) { // then merge q into p
                p->n_comp += q->n_comp + 1;
                p->seedcov = p->seedcov > q->seedcov? p->seedcov : q->seedcov;
                p->sub = p->sub > q->sub? p->sub : q->sub;
                p->csub = p->csub > q->csub? p->csub : q->csub;
                p->qb = q->qb, p->rb = q->rb;
                p->truesc = p->score = score;
                p->w = w;
                q->qb = q->qe;
            }
        }
    }
    for (i = 0, m = 0; i < n; ++i) // exclude identical hits
        if (a[i].qe > a[i].qb) {
            if (m != i) a[m++] = a[i];
            else ++m;
        }
    n = m;
    return m;
}
#endif

int mem_sort_dedup_patch(const mem_opt_t *opt, const bntseq_t *bns,
                         const uint8_t *pac, uint8_t *query, int n,
                         mem_alnreg_t *a, const int8_t *mat)
{
    /* D3 (--meth, PR-4): `mat` is the per-read OT/OB matrix and `query` the
     * original read bases (the caller threads both under --meth); outside --meth
     * `mat` is opt->mat and `query` the normal read, so behavior is unchanged.
     * All of a read's alnregs share one hypothesis under directional --meth, so a
     * single per-read matrix is correct for the colinear merge below. mat == NULL
     * (the default used by dedup-only callers, e.g. mate rescue passing 0,0,0)
     * resolves to opt->mat — those callers also pass bns==0 so no patch SW runs. */
    if (mat == NULL) mat = opt->mat;
    int m, i, j;
    if (n <= 1) return n;
    pdqsort_mem_ars2(n, a); // sort by the END position, not START!

    for (i = 0; i < n; ++i) a[i].n_comp = 1;
    for (i = 1; i < n; ++i)
    {
        mem_alnreg_t *p = &a[i];
        if (p->rid != a[i-1].rid || p->rb >= a[i-1].re + opt->max_chain_gap)
            continue; // then no need to go into the loop below

        for (j = i - 1; j >= 0 && p->rid == a[j].rid && p->rb < a[j].re + opt->max_chain_gap; --j) {
            mem_alnreg_t *q = &a[j];
            int64_t or_, oq, mr, mq;
            int score, w;
            if (q->qe == q->qb) continue; // a[j] has been excluded
            or_ = q->re - p->rb; // overlap length on the reference
            oq = q->qb < p->qb? q->qe - p->qb : p->qe - q->qb; // overlap length on the query
            mr = q->re - q->rb < p->re - p->rb? q->re - q->rb : p->re - p->rb; // min ref len in alignment
            mq = q->qe - q->qb < p->qe - p->qb? q->qe - q->qb : p->qe - p->qb; // min qry len in alignment
            if (or_ > opt->mask_level_redun * mr && oq > opt->mask_level_redun * mq) { // one of the hits is redundant
                if (p->score < q->score)
                {
                    p->qe = p->qb;
                    break;
                }
                else q->qe = q->qb;
            }
            else if (q->rb < p->rb && (score = mem_patch_reg(opt, bns, pac, query, q, p, &w, mat)) > 0) { // then merge q into p
                p->n_comp += q->n_comp + 1;
                p->seedcov = p->seedcov > q->seedcov? p->seedcov : q->seedcov;
                p->sub = p->sub > q->sub? p->sub : q->sub;
                p->csub = p->csub > q->csub? p->csub : q->csub;
                p->qb = q->qb, p->rb = q->rb;
                p->truesc = p->score = score;
                p->w = w;
                q->qb = q->qe;
            }
        }
    }
    for (i = 0, m = 0; i < n; ++i) // exclude identical hits
        if (a[i].qe > a[i].qb) {
            if (m != i) a[m++] = a[i];
            else ++m;
        }
    n = m;
    pdqsort_mem_ars(n, a);
    /* The historical post-sort exact-duplicate passes (mark then exclude regions
     * sharing (score, rb, qb)) have been removed: they are provably dead. Any two
     * regions with equal (rb, qb) have the shorter fully contained in the longer,
     * so on the reference or_ == mr and on the query oq >= mq; with the hardcoded
     * mask_level_redun = 0.95 (< 1) the sliding-window merge above always takes
     * its redundant branch (or_ > 0.95*mr && oq > 0.95*mq) and drops one of them.
     * Same-read alignments sharing rb are always within max_chain_gap, so the pair
     * is always in-window. Hence the exact-(score,rb,qb)-duplicate case is a strict
     * subset of what the merge already removes. Confirmed empirically: 0 removals
     * across 337M regions on HG002 WGS, and SAM output byte-identical. The by-score
     * sort is retained — its ordering is relied on downstream by mem_mark_primary_se
     * / mem_pair (removing it changes primary selection). */
    return n;
}


// return 1 if the seed is merged into the chain
static int test_and_merge(const mem_opt_t *opt, int64_t l_pac, mem_chain_t *c,
                          const mem_seed_t *p, int seed_rid, int tid)
{
    int64_t qend, rend, x, y;
    const mem_seed_t *last = &c->seeds[c->n-1];
    qend = last->qbeg + last->len;
    rend = last->rbeg + last->len;

    if (seed_rid != c->rid) return 0; // different chr; request a new chain
    if (p->qbeg >= c->seeds[0].qbeg && p->qbeg + p->len <= qend &&
        p->rbeg >= c->seeds[0].rbeg && p->rbeg + p->len <= rend)
        return 1; // contained seed; do nothing

    if ((last->rbeg < l_pac || c->seeds[0].rbeg < l_pac) &&
        p->rbeg >= l_pac) return 0; // don't chain if on different strand

    x = p->qbeg - last->qbeg; // always non-negtive
    y = p->rbeg - last->rbeg;
    if (y >= 0 && x - y <= opt->w && y - x <= opt->w &&
        x - last->len < opt->max_chain_gap &&
        y - last->len < opt->max_chain_gap) { // grow the chain
        if (c->n == c->m)
        {
            mem_seed_t *auxSeedBuf = NULL;
            int pm = c->m;
            c->m <<= 1;
            if (pm == SEEDS_PER_CHAIN) {  // re-new memory
                if ((auxSeedBuf = (mem_seed_t *) calloc(c->m, sizeof(mem_seed_t))) == NULL) { fprintf(stderr, "ERROR: out of memory auxSeedBuf\n"); exit(1); }
                memcpy((char*) (auxSeedBuf), c->seeds, c->n * sizeof(mem_seed_t));
                c->seeds = auxSeedBuf;
                tprof[PE13][tid]++;
            } else {  // new memory
                // fprintf(stderr, "[%0.4d] re-allocing old seed, m: %d\n", tid, c->m);
                if ((auxSeedBuf = (mem_seed_t *) realloc(c->seeds, c->m * sizeof(mem_seed_t))) == NULL) { fprintf(stderr, "ERROR: out of memory auxSeedBuf\n"); exit(1); }
                c->seeds = auxSeedBuf;
            }
            memset((char*) (c->seeds + c->n), 0, (c->m - c->n) * sizeof(mem_seed_t));
        }
        c->seeds[c->n++] = *p;
        return 1;
    }
    return 0; // request to add a new chain
}

int mem_seed_sw(const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac,
                int l_query, const uint8_t *query, const mem_seed_t *s,
                const int8_t *mat)
{
    /* D3 (--meth, PR-4): `mat` is the (possibly asymmetric, per-hypothesis)
     * substitution matrix for this chain; the caller passes opt->mat outside
     * --meth so behavior is unchanged there. */
    int qb, qe, rid;
    int64_t rb, re, mid, l_pac = bns->l_pac;
    uint8_t *rseq = 0;
    kswr_t x;

    if (s->len >= MEM_SHORT_LEN) return -1; // the seed is longer than the max-extend; no need to do SW
    qb = s->qbeg, qe = s->qbeg + s->len;
    rb = s->rbeg, re = s->rbeg + s->len;
    mid = (rb + re) >> 1;
    qb -= MEM_SHORT_EXT; qb = qb > 0? qb : 0;
    qe += MEM_SHORT_EXT; qe = qe < l_query? qe : l_query;
    rb -= MEM_SHORT_EXT; rb = rb > 0? rb : 0;
    re += MEM_SHORT_EXT; re = re < l_pac<<1? re : l_pac<<1;
    if (rb < l_pac && l_pac < re) {
        if (mid < l_pac) re = l_pac;
        else rb = l_pac;
    }
    if (qe - qb >= MEM_SHORT_LEN || re - rb >= MEM_SHORT_LEN) return -1; // the seed seems good enough; no need to do SW

    {
        static thread_local u8vec_scratch_t t_rseq;
        size_t want = (size_t)((re - rb) + 64);
        if (t_rseq.v.m < want) kv_resize(uint8_t, t_rseq.v, want);
        int64_t rlen;
        bns_fetch_seq_into(bns, pac, &rb, mid, &re, &rid, t_rseq.v.a, &rlen);
        rseq = t_rseq.v.a;
        // No qry-profile cache: each seed slices a different sub-query
        // (query+qb, qe-qb), so ksw_align2 always builds a fresh profile.
        x = ksw_align2(qe - qb, (uint8_t*)query + qb, re - rb, rseq, 5, mat, opt->o_del, opt->e_del, opt->o_ins, opt->e_ins, KSW_XSTART, NULL);
        /* rseq aliases thread-local scratch; do not free. */
    }
    return x.score;
}

int mem_chain_weight(const mem_chain_t *c)
{
    int64_t end;
    int j, w = 0, tmp;
    for (j = 0, end = 0; j < c->n; ++j) {
        const mem_seed_t *s = &c->seeds[j];
        if (s->qbeg >= end) w += s->len;
        else if (s->qbeg + s->len > end) w += s->qbeg + s->len - end;
        end = end > s->qbeg + s->len? end : s->qbeg + s->len;
    }
    tmp = w; w = 0;
    for (j = 0, end = 0; j < c->n; ++j) {
        const mem_seed_t *s = &c->seeds[j];
        if (s->rbeg >= end) w += s->len;
        else if (s->rbeg + s->len > end) w += s->rbeg + s->len - end;
        end = end > s->rbeg + s->len? end : s->rbeg + s->len;
    }
    w = w < tmp? w : tmp;
    return w < 1<<30? w : (1<<30)-1;
}

void mem_print_chain(const bntseq_t *bns, mem_chain_v *chn)
{
    int i, j;
    for (i = 0; i < chn->n; ++i)
    {
        mem_chain_t *p = &chn->a[i];
        fprintf(stderr, "* Found CHAIN(%d): n=%d; weight=%d", i, p->n, mem_chain_weight(p));
        for (j = 0; j < p->n; ++j)
        {
            bwtint_t pos;
            int is_rev;
            pos = bns_depos(bns, p->seeds[j].rbeg, &is_rev);
            if (is_rev) pos -= p->seeds[j].len - 1;
            fprintf(stderr, "\t%d;%d;%d,%ld(%s:%c%ld)",
                       p->seeds[j].score, p->seeds[j].len, p->seeds[j].qbeg,
                       (long)p->seeds[j].rbeg, bns->anns[p->rid].name,
                       "+-"[is_rev], (long)(pos - bns->anns[p->rid].offset) + 1);
        }
        fputc('\n', stderr);
    }
}

// Skip-short-seed extension filter: in a chain that still has a long-enough
// anchor, drop seeds shorter than min_ext_len so they are never extended (their
// banded Smith-Waterman is skipped downstream) -- those short seeds are collinear
// with the anchor and its extension already covers them, so dropping them is
// near output-neutral and pure speed.
//
// A chain with NO seed >= min_ext_len is left untouched: its only evidence is
// short, so dropping it would empty the chain and the read would go unmapped.
// Such all-short chains (common on low-mappability / repetitive reads) therefore
// extend exactly as the default does. This guarantees the filter never empties a
// non-empty chain, so it can only reduce extension work, never lose a read.
//
// Stable compaction -- surviving seeds keep their order. Returns the new seed
// count. min_ext_len <= 0 is a no-op (default), byte-identical to baseline.
int mem_chain_drop_short_seeds(mem_chain_t *c, int min_ext_len)
{
    if (min_ext_len <= 0) return c->n;
    int j, k;
    for (j = 0; j < c->n; ++j)
        if (c->seeds[j].len >= min_ext_len) break;
    if (j == c->n) return c->n;            /* no anchor: all-short chain, leave intact */
    for (j = k = 0; j < c->n; ++j)
        if (c->seeds[j].len >= min_ext_len)
            c->seeds[k++] = c->seeds[j];
    c->n = k;
    return c->n;
}

void mem_flt_chained_seeds(const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac,
                           bseq1_t *seq_, int n_chn, mem_chain_t *a)
{
    //double min_l = opt->min_chain_weight?
    // MEM_HSP_COEF * opt->min_chain_weight : MEM_MINSC_COEF * log(l_query);

    int i, j, k;// min_HSP_score = (int)(opt->a * min_l + .499);
    //if (min_l > MEM_SEEDSW_COEF * l_query) return; // don't run the following for short reads

    /* D3 (--meth, PR-4): the chained-seed SW filter must score the ORIGINAL read
     * against the original ref with the per-chain asymmetric matrix, same as
     * extension — otherwise it over-penalizes bisulfite conversions and drops
     * borderline seeds (the PR-3 placeholder concern). meth_qbuf is the reusable
     * 2-bit encode of meth_orig_seq (same orientation as seq, per bwa.h). */
    uint8_t *meth_qbuf = NULL;
    int      meth_qbuf_cap = 0;

    for (i = 0; i < n_chn; ++i)
    {
        mem_chain_t *c = &a[i];
        const uint8_t *query = (uint8_t*) seq_[c->seqid].seq;
        int l_query = seq_[c->seqid].l_seq;

        // Skip-short-seed extension filter (opt->min_ext_len). Placed BEFORE the
        // MEM_SEEDSW_COEF screen `continue` below, which is skipped for typical
        // read lengths -- dropping seeds inside that loop would no-op on the main
        // workload. Off (min_ext_len==0): no-op, output byte-identical.
        mem_chain_drop_short_seeds(c, opt->min_ext_len);
        if (c->n == 0) continue;  /* Defensive: drop_short_seeds never empties a
                                   * non-empty chain (it leaves all-short chains
                                   * intact), but the meth block below reads
                                   * seeds[0], so guard regardless. */

        /* D3 (--meth, PR-4): swap to the original read + per-hypothesis matrix. */
        const int8_t *mat = opt->mat;
        if (opt->meth_mode && seq_[c->seqid].meth_orig_seq != NULL) {
            if (l_query > meth_qbuf_cap) {
                meth_qbuf_cap = l_query;
                /* CodeRabbit: don't assign realloc() result back to meth_qbuf
                 * directly — a NULL return would leak the old buffer. */
                uint8_t *meth_qbuf_new = (uint8_t *) realloc(meth_qbuf, (size_t) meth_qbuf_cap);
                assert(meth_qbuf_new != NULL);
                meth_qbuf = meth_qbuf_new;
            }
            const char *os = seq_[c->seqid].meth_orig_seq;
            for (int q = 0; q < l_query; ++q) {
                unsigned char ch = (unsigned char) os[q];
                meth_qbuf[q] = (ch < 4) ? ch : nst_nt4_table[ch];
            }
            query = meth_qbuf;
            /* D3 (--meth, fix): strand-adjusted hypothesis — this chain's seeds
             * are single-strand; a reverse chain (seeds[0].rbeg >= l_pac) extends
             * against the RC reference, flipping the conversion's freed cell. */
            int hyp = c->meth_hypothesis;
            if (hyp >= 0 && c->seeds[0].rbeg >= bns->l_pac) hyp ^= 1;
            mat = mem_opt_meth_mat(opt, hyp);
        }

        double min_l = opt->min_chain_weight?
        MEM_HSP_COEF * opt->min_chain_weight : MEM_MINSC_COEF * log(l_query);
        int min_HSP_score = (int)(opt->a * min_l + .499);
        if (min_l > MEM_SEEDSW_COEF * l_query) continue;

        for (j = k = 0; j < c->n; ++j)
        {
            mem_seed_t *s = &c->seeds[j];
            s->score = mem_seed_sw(opt, bns, pac, l_query, query, s, mat);
            if (s->score < 0 || s->score >= min_HSP_score)
            {
                s->score = s->score < 0? s->len * opt->a : s->score;
                c->seeds[k++] = *s;
            }
        }
        c->n = k;
    }
    free(meth_qbuf);  /* NULL outside --meth; free() is NULL-safe */
}

/* --max-extend-chains: cap the number of chains that reach banded-SW extension
 * to the top `max_n` by chain weight (applied after mem_chain_flt). The dropped
 * chains are the lowest-weight secondaries; on holodeck truth (sim-wgs-place)
 * this keeps high-confidence placement accuracy essentially unchanged while
 * cutting extension work -- errors it does add land in low-MAPQ multimappers,
 * not confident calls. Opt-in speed lever (bundled into --fast); NOT
 * byte-identical (drops candidate secondaries, so XS/secondary/MAPQ can move on
 * multi-mapping reads). Always keeps >= 1 chain. Returns the new chain count. */
#define MAX_EXTEND_CHAINS_CAP 4096
static int mem_chain_cap_extend(mem_chain_t *a, int n, int max_n)
{
    if (max_n <= 0 || n <= max_n || n > MAX_EXTEND_CHAINS_CAP) return n;
    int w[MAX_EXTEND_CHAINS_CAP];
    for (int i = 0; i < n; i++) w[i] = a[i].w > 0 ? a[i].w : mem_chain_weight(&a[i]);
    /* keep chain i iff fewer than max_n other chains outrank it by
     * (weight desc, then original index asc) -- a stable top-max_n selection. */
    int k = 0;
    for (int i = 0; i < n; i++) {
        int rank = 0;
        for (int j = 0; j < n; j++) {
            if (j == i) continue;
            if (w[j] > w[i] || (w[j] == w[i] && j < i)) rank++;
        }
        if (rank < max_n) a[k++] = a[i];
        else if (a[i].m > SEEDS_PER_CHAIN) free(a[i].seeds);
    }
    return k;
}

/* --extend-mate-concordant: like mem_chain_cap_extend but additionally RETAINS
 * any chain concordant with a mate chain (same contig, FR orientation, within
 * `win` bp) even if it ranks below max_n. Targets the meth --fast pairing
 * regression: the true chain is often low-weight under bisulfite (collapsed
 * alphabet) and gets capped, yet it anchors the true concordant pair, so PE
 * pairing flips both mates to a wrong concordant locus. chain.pos is in the
 * original 2*l_pac coordinate space, so bns_depos gives the forward coord.
 * `win` is the estimated proper-pair insert high bound (pes[FR].high) when
 * available, else MATE_CONCORDANT_WINDOW_FALLBACK, else a fixed CLI value;
 * matching the aligner's own concordance bound keeps only genuine pair anchors
 * (far/spurious concordant chains stay capped, bounding the extra extension).
 * The mate scan is bounded (MATE_SCAN_MAX) to keep this O(n) in practice; the
 * true anchor is virtually always among a read's higher-ranked chains. */
#define MATE_CONCORDANT_WINDOW_FALLBACK 1000  /* used only in --auto before the insert size is estimated (e.g. first chunk) */
#define MATE_SCAN_MAX 256
static int mem_chain_cap_extend_mate(mem_chain_t *a, int n, int max_n,
        const mem_chain_t *mate, int mate_n, const bntseq_t *bns, int64_t win)
{
    if (max_n <= 0 || n <= max_n || n > MAX_EXTEND_CHAINS_CAP) return n;
    const int MSCAN = MATE_SCAN_MAX;
    int mn = mate_n < MSCAN ? mate_n : MSCAN;
    int64_t mfp[MSCAN]; int mrid[MSCAN], mrev[MSCAN];
    for (int j = 0; j < mn; j++) {
        int is_rev; mfp[j] = bns_depos(bns, mate[j].pos, &is_rev);
        mrid[j] = mate[j].rid; mrev[j] = is_rev;
    }
    int w[MAX_EXTEND_CHAINS_CAP];
    for (int i = 0; i < n; i++) w[i] = a[i].w > 0 ? a[i].w : mem_chain_weight(&a[i]);
    int k = 0;
    for (int i = 0; i < n; i++) {
        int rank = 0;
        for (int j = 0; j < n; j++)
            if (j != i && (w[j] > w[i] || (w[j] == w[i] && j < i))) rank++;
        int keep = rank < max_n;
        if (!keep && mn > 0) {
            int is_rev; int64_t fp = bns_depos(bns, a[i].pos, &is_rev);
            for (int j = 0; j < mn; j++)
                if (mrid[j] == a[i].rid && mrev[j] != is_rev &&
                    /* true FR ("innie"): the forward-strand chain must sit
                     * at/upstream of the reverse-strand chain, else an RF
                     * ("outie") pair within win bp would falsely qualify. */
                    ((!is_rev && fp <= mfp[j]) || (is_rev && mfp[j] <= fp))) {
                    int64_t d = fp > mfp[j] ? fp - mfp[j] : mfp[j] - fp;
                    if (d <= win) { keep = 1; break; }
                }
        }
        if (keep) a[k++] = a[i];
        else if (a[i].m > SEEDS_PER_CHAIN) free(a[i].seeds);
    }
    return k;
}

int mem_chain_flt(const mem_opt_t *opt, int n_chn_, mem_chain_t *a_, int tid)
{
    int i, k, n_numc = 0;
    if (n_chn_ == 0) return 0; // no need to filter
    // compute the weight of each chain and drop chains with small weight
    for (i = k = 0; i < n_chn_; ++i)
    {
        mem_chain_t *c = &a_[i];
        c->first = -1; c->kept = 0;
        c->w = mem_chain_weight(c);
        if (c->w < opt->min_chain_weight)
        {
            if (c->m > SEEDS_PER_CHAIN) {
                tprof[PE11][tid] ++;
                free(c->seeds);
            }
            //free(c->seeds);
        }
        else a_[k++] = *c;
    }
    n_chn_ = k;

    /* Fast path for the dominant single-chain (unique-mapper) case. With one
     * surviving chain the range split is a single [0,1) range and the pairwise-
     * overlap pass reduces to `a_[0].kept = 3` (set unconditionally, kept through
     * compaction) returning 1 — but the general path still builds a std::vector,
     * a kvec, and calls ks_introsort to get there. a_[0].first was already reset
     * to -1 in the weight-filter loop above. (Only n_chn_ == 1 is short-circuited;
     * n_chn_ == 0 falls through to preserve the existing general-path behavior.) */
    if (n_chn_ == 1) {
        a_[0].kept = 3;
        return 1;
    }

    std::vector<std::pair<int, int> > range;
    std::pair<int, int> pr;
    int pseqid = a_[0].seqid;
    pr.first = 0;
    for (i=1; i<n_chn_; i++)
    {
        mem_chain_t *c =&a_[i];
        if (c->seqid != pseqid) {
            //if (flag == -1) {
            //  pr.first = i;
            //flag = 1;
            //}
            //else
            {
                pr.second = i;
                range.push_back(pr);
                pr.first = i;
                // flag = -1;
            }
        }
        pseqid = c->seqid;
    }
    pr.second = i;
    range.push_back(pr);

    int ilag = 0;
    for (int l=0; l<range.size(); l++)
    {
        // this keeps int indices of the non-overlapping chains
        kvec_t(int) chains = {0,0,0};
        mem_chain_t *a =&a_[range[l].first];
        int n_chn = range[l].second - range[l].first;
        // original code block starts
        ks_introsort(mem_flt, n_chn, a);

        // pairwise chain comparisons
        a[0].kept = 3;
        kv_push(int, chains, 0);
        for (i = 1; i < n_chn; ++i)
        {
            int large_ovlp = 0;
            for (k = 0; k < chains.n; ++k)
            {
                int j = chains.a[k];
                int b_max = chn_beg(a[j]) > chn_beg(a[i])? chn_beg(a[j]) : chn_beg(a[i]);
                int e_min = chn_end(a[j]) < chn_end(a[i])? chn_end(a[j]) : chn_end(a[i]);
                if (e_min > b_max && (!a[j].is_alt || a[i].is_alt)) { // have overlap; don't consider ovlp where the kept chain is ALT while the current chain is primary
                    int li = chn_end(a[i]) - chn_beg(a[i]);
                    int lj = chn_end(a[j]) - chn_beg(a[j]);
                    int min_l = li < lj? li : lj;
                    if (e_min - b_max >= min_l * opt->mask_level && min_l < opt->max_chain_gap) { // significant overlap
                        large_ovlp = 1;
                        if (a[j].first < 0) a[j].first = i; // keep the first shadowed hit s.t. mapq can be more accurate
                        if (a[i].w < a[j].w * opt->drop_ratio && a[j].w - a[i].w >= opt->min_seed_len<<1)
                            break;
                    }
                }
            }
            if (k == chains.n)
            {
                kv_push(int, chains, i);
                a[i].kept = large_ovlp? 2 : 3;
            }
        }
        for (i = 0; i < chains.n; ++i)
        {
            mem_chain_t *c = &a[chains.a[i]];
            if (c->first >= 0) a[c->first].kept = 1;
        }
        free(chains.a);
        for (i = k = 0; i < n_chn; ++i) { // don't extend more than opt->max_chain_extend .kept=1/2 chains
            if (a[i].kept == 0 || a[i].kept == 3) continue;
            if (++k >= opt->max_chain_extend) break;
        }

        for (; i < n_chn; ++i)
            if (a[i].kept < 3) a[i].kept = 0;

        for (i = k = 0; i < n_chn; ++i)  // free discarded chains
        {
            mem_chain_t *c = &a[i];
            if (c->kept == 0)
            {
                if (c->m > SEEDS_PER_CHAIN) {
                    tprof[PE11][tid] ++;
                    free(c->seeds);
                }
                //free(c->seeds);
            }
            else a[k++ - ilag] = a[i];
        }
        // original code block ends
        ilag += n_chn - k;
        n_numc += k;
    }

    return n_numc;
}

SMEM *mem_collect_smem(FMI_search *fmi, const mem_opt_t *opt,
                       const bseq1_t *seq_,
                       int nseq,
                       SMEM *matchArray,
                       int32_t *min_intv_ar,
                       int32_t *query_pos_ar,
                       uint8_t *enc_qdb,
                       int32_t *rid,
                       mem_cache *mmc,
                       int64_t &tot_smem,
                       int tid)
{
    int64_t pos = 0;
    int split_len = (int)(opt->min_seed_len * opt->split_factor + .499);
    int64_t num_smem1 = 0, num_smem2 = 0, num_smem3 = 0;
    int max_readlength = -1;

    int32_t *query_cum_len_ar = (int32_t *)_mm_malloc(nseq * sizeof(int32_t), 64);

    int offset = 0;
    for (int l=0; l<nseq; l++)
    {
        min_intv_ar[l] = 1;
        memcpy(enc_qdb + offset, seq_[l].seq, (size_t)seq_[l].l_seq);
        offset += seq_[l].l_seq;
        rid[l] = l;
    }

    max_readlength = seq_[0].l_seq;
    query_cum_len_ar[0] = 0;
    for(int i = 1; i < nseq; i++) {
        query_cum_len_ar[i] = query_cum_len_ar[i - 1] + seq_[i-1].l_seq;
        if (max_readlength < seq_[i].l_seq)
            max_readlength = seq_[i].l_seq;
    }

    fmi->getSMEMsAllPosOneThread(enc_qdb, min_intv_ar, rid, nseq, nseq,
                                 seq_, query_cum_len_ar, max_readlength, opt->min_seed_len,
                                 matchArray, &num_smem1);

    for (int64_t i=0; i<num_smem1; i++)
    {
        SMEM *p = &matchArray[i];
        int start = p->m, end = p->n +1;
        if (end - start < split_len || p->s > opt->split_width)
            continue;

        int len = seq_[p->rid].l_seq;

        rid[pos] = p->rid;

        query_pos_ar[pos] = (end + start)>>1;

        assert(query_pos_ar[pos] < len);

        min_intv_ar[pos] = p->s + 1;
        pos ++;
    }

    #if 1
    // Pre-walk grow: getSMEMsOnePosOneThread emits at most readlength SMEMs
    // per re-walk position (the scalar walker's worst case is one SMEM per
    // backward step plus the final prev[0]), so the safe upper bound on
    // num_smem2 is pos * max_readlength.
    const int64_t safe_smem2_cap = (int64_t)pos * max_readlength;
    if (mmc->wsize_mem[tid] < num_smem1 + safe_smem2_cap) {
        int64_t tmp = mmc->wsize_mem[tid];
        mmc->wsize_mem[tid] = num_smem1 + safe_smem2_cap;

        SMEM* ptr1 = mmc->matchArray[tid];
        // ptr2 = w.mmc.min_intv_ar;
        // ptr3 =
        if (bwa_verbose >= 4) {
            fprintf(stderr, "[%0.4d] REA Re-allocating SMEM after num_smem1: "
                    "%" PRId64 " -> %" PRId64 "\n",
                    tid, tmp, mmc->wsize_mem[tid]);
        }
        mmc->matchArray[tid]    = (SMEM *) _mm_malloc(mmc->wsize_mem[tid] * sizeof(SMEM), 64);
        assert(mmc->matchArray[tid] != NULL);
        matchArray = mmc->matchArray[tid];
        // w.mmc.min_intv_ar[l]   = (int32_t *) malloc(w.mmc.wsize_mem[l] * sizeof(int32_t));
        // w.mmc.query_pos_ar[tid]  = (int16_t *) malloc(w.mmc.wsize_mem[l] * sizeof(int16_t));
        // w.mmc.enc_qdb[l]       = (uint8_t *) malloc(w.mmc.wsize_mem[l] * sizeof(uint8_t));
        // w.mmc.rid[l]           = (int32_t *) malloc(w.mmc.wsize_mem[l] * sizeof(int32_t));
        //w.mmc.lim[l]         = (int32_t *) _mm_malloc((BATCH_SIZE + 32) * sizeof(int32_t), 64);
        for (int i=0; i<num_smem1; i++) {
            mmc->matchArray[tid][i] = ptr1[i];
        }
        _mm_free(ptr1);
    }
    #endif
    
#if SMEM_LOCKSTEP_N > 1
    fmi->getSMEMsOnePosOneThread_lockstep(enc_qdb,
                                          query_pos_ar,
                                          min_intv_ar,
                                          rid,
                                          pos,
                                          pos,
                                          seq_,
                                          query_cum_len_ar,
                                          max_readlength,
                                          opt->min_seed_len,
                                          matchArray + num_smem1,
                                          &num_smem2);
#else
    fmi->getSMEMsOnePosOneThread(enc_qdb,
                                 query_pos_ar,
                                 min_intv_ar,
                                 rid,
                                 pos,
                                 pos,
                                 seq_,
                                 query_cum_len_ar,
                                 max_readlength,
                                 opt->min_seed_len,
                                 matchArray + num_smem1,
                                 &num_smem2);
#endif
    // assert(mmc->wsize_mem[tid] > (num_smem1 + num_smem2));
    // fprintwsize_mem_rf(stderr, "num_smem2: %d\n", num_smem2);
    if (opt->max_mem_intv > 0)
    {
        #if 1
		// Pre-walk grow for the third pass. bwtSeedStrategyAllPosOneThread's
		// documented worst case (FMI_search.h) is nseq * max_seq_length SMEMs
		// — one per query position per read — which dwarfs the prior
		// offset/(min_seed_len+1) heuristic. Use the safe bound so the
		// caller cannot overflow matchArray when this pass writes at
		// (num_smem1 + num_smem2).
		const int64_t safe_smem3_cap = (int64_t)nseq * max_readlength;
		const int64_t mem_intv_cap = num_smem1 + num_smem2 + safe_smem3_cap;
		if (mmc->wsize_mem[tid] < mem_intv_cap) {
            int64_t tmp = mmc->wsize_mem[tid];
		    mmc->wsize_mem[tid] = mem_intv_cap;
		    SMEM* ptr1 = mmc->matchArray[tid];
            if (bwa_verbose >= 4) {
                fprintf(stderr, "[%0.4d] REA Re-allocating SMEM after num_smem2: "
                        "%" PRId64 " -> %" PRId64 "\n",
                        tid, tmp, mmc->wsize_mem[tid]);
            }
		    mmc->matchArray[tid]    = (SMEM *) _mm_malloc(mmc->wsize_mem[tid] * sizeof(SMEM), 64);
		    assert(mmc->matchArray[tid] != NULL);
		    matchArray = mmc->matchArray[tid];
		    // First pass wrote num_smem1 records at offset 0; the second pass
		    // (getSMEMsOnePosOneThread*) just wrote num_smem2 records at offset
		    // num_smem1. Preserve both — bwtSeedStrategyAllPosOneThread() below
		    // appends at offset num_smem1 + num_smem2.
		    int64_t already_written = (int64_t)num_smem1 + num_smem2;
		    for (int64_t i = 0; i < already_written; ++i) {
		        mmc->matchArray[tid][i] = ptr1[i];
		    }
		    _mm_free(ptr1);
		}
        #endif
        for (int l=0; l<nseq; l++)
            min_intv_ar[l] = opt->max_mem_intv;

// Third-pass re-seeding lockstep is gated to arm64. It overlaps the serial
// forward-extension chain's cp_occ misses with N independent reads' chains,
// which is a large seeding win (~-16%, ~-1.7% end-to-end) on non-SMT arm
// (Graviton4, Apple Silicon). On SMT x86 the sibling hyperthread already hides
// that latency, so the lockstep is redundant overhead and regresses ~+2% at
// full-vCPU thread counts — measured across BWTSEED_LOCKSTEP_N in {4,8,16},
// all depths regress x86-SMT. Byte-identical either way (see bwtseed lockstep
// parity harness); the scalar path stays the x86 default.
#if defined(__aarch64__) && BWTSEED_LOCKSTEP_N > 1
        num_smem3 = fmi->bwtSeedStrategyAllPosOneThread_lockstep(enc_qdb, min_intv_ar,
                                                                 nseq, seq_, query_cum_len_ar,
                                                                 opt->min_seed_len + 1,
                                                                 matchArray + num_smem1 + num_smem2,
                                                                 max_readlength);
#else
        num_smem3 = fmi->bwtSeedStrategyAllPosOneThread(enc_qdb, min_intv_ar,
                                                        nseq, seq_, query_cum_len_ar,
                                                        opt->min_seed_len + 1,
                                                        matchArray + num_smem1 + num_smem2);
#endif
    }
    tot_smem = num_smem1 + num_smem2 + num_smem3;
    // assert(mmc->wsize_mem[tid] > (tot_smem));
    // fprintf(stderr, "num_smems: %d %d %d, %d\n", num_smem1, num_smem2, num_smem3, tot_smem);
    fmi->sortSMEMs(matchArray, &tot_smem, nseq, seq_[0].l_seq, 1); // seq_[0].l_seq - only used for blocking when using nthreads

    pos = 0;
    int64_t smem_ptr = 0;
    for (int l=0; l<nseq && pos < tot_smem - 1; l++) {
        pos = smem_ptr - 1;
        do {
            pos++;
        } while (pos < tot_smem - 1 && matchArray[pos].rid == matchArray[pos + 1].rid);
        int64_t n = pos + 1 - smem_ptr;

        if (n > 0)
            ks_introsort(mem_intv1, n, &matchArray[smem_ptr]);
        smem_ptr = pos + 1;
    }

    if (opt->smem_dedup)
        tot_smem = smem_dedup_inplace(matchArray, tot_smem);

    _mm_free(query_cum_len_ar);

    return matchArray;
}

/* D3 (--meth, PR-3): remap one seed hit from SEED-index (f/r-doubled) coordinates
 * to ORIGINAL-reference coordinates + a bisulfite hypothesis label.
 *
 * The seed reference has 2N contigs in r0,f0,r1,f1,... order (rchrK = G→A / OB,
 * fchrK = C→T / OT); both are FORWARD projections of the original contig K, so
 * the seed-local position equals the original-local position. The seed FM-index
 * additionally doubles its own pac (forward + reverse-complement), so a seed hit
 * `seed_rbeg` lives in [0, 2*L_seed) where >= L_seed means the hit is on the RC
 * of the seed reference. That RC strand is the genomic alignment strand and is
 * INDEPENDENT of the OT/OB hypothesis (which is the f/r contig parity).
 *
 * Inputs:
 *   seed_bns   — the seed (f/r-doubled) BNS, used to decode seed_rbeg.
 *   orig_bns   — the original BNS (N contigs), the remap target.
 *   seed_rbeg  — seed-index doubled-pac coordinate of the seed start.
 *   seed_len   — seed length (>0); used to project the start across strands.
 * Outputs (only written on success):
 *   *orig_rbeg   — start coordinate in ORIGINAL doubled-pac space [0, 2*L_orig).
 *   *orig_rid    — original contig id.
 *   *hypothesis  — 1 = OT (odd "f" seed contig), 0 = OB (even "r" seed contig).
 * Returns 0 on success, -1 if the hit cannot be remapped (bridging / out of
 * range), in which case the caller should drop the seed.
 *
 * NOTE: the reverse-strand re-encoding below is the single most placement-
 * critical step in --meth (spec §5.2 / plan §5: "strand bookkeeping can silently
 * corrupt placement"). It is exercised end-to-end and covered by the seed→
 * original remap + OT/OB PE/strand placement regression tests
 * (test/regression/meth_seed_index.sh and the live mem --meth placement tests).
 */
static inline int meth_seed_to_orig(const bntseq_t *seed_bns,
                                    const bntseq_t *orig_bns,
                                    int64_t seed_rbeg, int32_t seed_len,
                                    int64_t *orig_rbeg, int *orig_rid,
                                    int8_t *hypothesis)
{
    if (seed_bns == NULL || orig_bns == NULL) return -1;
    if (seed_rbeg < 0 || seed_rbeg >= (seed_bns->l_pac << 1)) return -1;

    /* Decode the seed start to a forward seed coordinate + RC flag. */
    int seed_is_rev = 0;
    int64_t seed_fwd = bns_depos(seed_bns, seed_rbeg, &seed_is_rev); /* [0,L_seed) */
    int seed_rid = bns_pos2rid(seed_bns, seed_fwd);
    if (seed_rid < 0) return -1;

    /* Seed contig order is r0,f0,r1,f1,...: orig contig = rid/2, OT iff odd. */
    int o_rid = seed_rid >> 1;
    if (o_rid < 0 || o_rid >= orig_bns->n_seqs) return -1;
    *hypothesis = (int8_t)(seed_rid & 1); /* 1=f=OT(C→T), 0=r=OB(G→A) */

    /* f/r seed contigs are forward projections of the original, so the seed-local
     * offset is the original-local offset. Reject seeds that bridge the contig
     * end (parity with bns_intv2rid, which returns -1 on a multi-contig / f-r
     * boundary span). seed_fwd is the FORWARD coordinate of the seed start: for a
     * forward seed the forward span is [seed_local, seed_local+span); for a
     * reverse seed bns_depos returns the rightmost forward base, so the forward
     * span is (seed_local-span, seed_local]. Require that whole span to fit the
     * seed contig (and, since the original shares the layout, the original too). */
    int64_t seed_local = seed_fwd - seed_bns->anns[seed_rid].offset;
    int64_t span = (seed_len > 0) ? (int64_t)seed_len : 1;
    int64_t span_lo = seed_is_rev ? (seed_local - span + 1) : seed_local;
    int64_t span_hi = seed_is_rev ? (seed_local + 1) : (seed_local + span); /* exclusive */
    int64_t seed_contig_len = seed_bns->anns[seed_rid].len;
    int64_t orig_contig_len = orig_bns->anns[o_rid].len;
    if (span_lo < 0 || span_hi > seed_contig_len) return -1;
    if (span_hi > orig_contig_len)                return -1;
    int64_t orig_fwd = orig_bns->anns[o_rid].offset + seed_local; /* [0,L_orig) */

    /* Re-encode the genomic strand into the ORIGINAL doubled-pac space.
     * NOTE: this mirrors the doubled-pac arithmetic used throughout pairing /
     * mate rescue (e.g. b.rb = (l_pac<<1) - (rb + te + 1)).
     * The seed and original references share identical per-contig length and
     * forward layout, so the transform is a pure per-contig forward-coordinate
     * offset shift (seed_offset → orig_offset, already applied to get orig_fwd),
     * re-encoded to the SAME strand. bns_depos is its own inverse on a single
     * position, so the START coordinate is preserved exactly across the shift;
     * extension re-derives the span. */
    if (!seed_is_rev) {
        *orig_rbeg = orig_fwd;                                   /* forward */
    } else {
        *orig_rbeg = (orig_bns->l_pac << 1) - 1 - orig_fwd;      /* reverse */
    }
    /* NOTE: the reverse-strand re-encode above is THE placement-critical line; it
     * is covered by the seed→original remap round-trip and the OT/OB paired-end
     * strand placement regression tests. */
    if (*orig_rbeg < 0 || *orig_rbeg >= (orig_bns->l_pac << 1)) return -1;
    *orig_rid = o_rid;
    return 0;
}

/* Chain one fully-resolved seed into the per-read kbtree: find the closest
 * existing chain and test_and_merge, else open a new chain. Extracted verbatim
 * from mem_chain_seeds's per-seed body so the two seed-emit paths share it: the
 * default (SEED_ORDER_OFF) streams each seed straight here as it is resolved
 * (single pass, no buffering), while the opt-in reorder modes materialize all
 * seeds, order them, then replay them through here. Behaviour is identical to
 * the previous inline body; only the call site differs. */
static inline void chain_add_one_seed(const mem_opt_t *opt, int64_t l_pac,
                                      const bntseq_t *chain_bns, kbtree_t(chn) *tree,
                                      mem_seed_t *seedBuf, int64_t *seedBufCount,
                                      int64_t seedBufSize, int tid, int seqid,
                                      int *num_seqid, const mem_seed_t *seed_in,
                                      int rid, int8_t meth_hyp)
{
    mem_seed_t s = *seed_in;
    mem_chain_t tmp, *lower, *upper;
    int to_add = 0;
    tmp.pos = s.rbeg;

    if (kb_size(tree))
    {
        kb_intervalp(chn, tree, &tmp, &lower, &upper); // find the closest chain

        if (!lower || !test_and_merge(opt, l_pac, lower, &s, rid, tid))
            to_add = 1;
    }
    else to_add = 1;

    if (to_add) // add the seed as a new chain
    {
        tmp.n = 1; tmp.m = SEEDS_PER_CHAIN;
        if((*seedBufCount + tmp.m) > seedBufSize)
        {
            tmp.m += 1;
            tmp.seeds = (mem_seed_t *)calloc (tmp.m, sizeof(mem_seed_t));
            assert(tmp.seeds != NULL);
            tprof[PE13][tid]++;
        }
        else {
            tmp.seeds = seedBuf + *seedBufCount;
            *seedBufCount += tmp.m;
        }
        memset((char*) (tmp.seeds), 0, tmp.m * sizeof(mem_seed_t));
        tmp.seeds[0] = s;
        tmp.rid = rid;
        tmp.seqid = seqid;
        /* is_alt indexes the chain-side bns (original in --meth). */
        tmp.is_alt = !!chain_bns->anns[rid].is_alt;
        /* D3: carry the OT/OB hypothesis on the chain (-1 non-meth). */
        tmp.meth_hypothesis = meth_hyp;
        kb_putp(chn, tree, &tmp);
        (*num_seqid)++;
    }
}

/** NEW ONE **/
void mem_chain_seeds(FMI_search *fmi, const mem_opt_t *opt,
                     const bntseq_t *bns,
                     const bntseq_t *meth_orig_bns,
                     const bseq1_t *seq_,
                     int nseq,
                     int tid,
                     mem_chain_v *chain_ar,
                     mem_seed_t *seedBuf,
                     int64_t seedBufSize,
                     SMEM *matchArray,
                     int64_t num_smem)
{
    int b, e, l_rep, size = 0;
    int64_t i, pos = 0;
    int64_t smem_ptr = 0;
    /* `bns` is the SEED (f/r-doubled) BNS used to decode SA coordinates. In
     * --meth we remap each seed hit into ORIGINAL coordinates before chaining, so
     * the chain-side bns / l_pac (rid lookup, strand-boundary merge test) must be
     * the ORIGINAL ones. Outside --meth chain_bns == bns. */
    const int meth_remap = (opt->meth_mode && meth_orig_bns != NULL);
    const bntseq_t *chain_bns = meth_remap ? meth_orig_bns : bns;
    int64_t l_pac = chain_bns->l_pac;
    /* Seed-emit strategy. OFF (the default) streams each resolved seed straight
     * into the kbtree in resolve order — no recs[] buffer, no order_seeds pass —
     * which is byte-identical to, and as cheap as, the pre-seed-order baseline.
     * Only a selected reorder mode needs the materialize/order/replay path. */
    const int reorder = (opt->seed_emit_order != SEED_ORDER_OFF);

    int num[nseq];
    memset(num, 0, nseq*sizeof(int));
    int smem_buf_size = 6000;
    int64_t *sa_coord = (int64_t *) _mm_malloc(sizeof(int64_t) * opt->max_occ * smem_buf_size, 64);
    int64_t seedBufCount = 0;

    /* Seed-order refactor (Phase A buffer, Spec S4): one fully-resolved seed per
     * SA hit, sized by RESOLVED-seed count (the prefetch ceiling), NOT SMEM count.
     * Persisted across reads in this call and freed at function end with sa_coord. */
    seed_rec_t *recs = NULL;
    int64_t recs_cap = 0;

    for (int l=0; l<nseq; l++)
        kv_init(chain_ar[l]);

    // filter seq at early stage than this!, shifted to collect!!!
    // if (len < opt->min_seed_len) return chain; // if the query is shorter than the seed length, no match

    uint64_t tim = __rdtsc();
    for (int l=0; l<nseq && pos < num_smem - 1; l++)
    {
        // addition, FIX FIX FIX!!!!! THIS!!!!
        //if (aux_->mem.n == 0) continue;

        if (matchArray[smem_ptr].rid > l) continue;
        if (seq_[l].l_seq < opt->min_seed_len) continue;
        assert(matchArray[smem_ptr].rid == l);

        kbtree_t(chn) *tree;
        tree = kb_init(chn, KB_DEFAULT_SIZE + 8); // +8, due to addition of counters in chain
        mem_chain_v *chain = &chain_ar[l];
        size = 0;

        b = e = l_rep = 0;
        pos = smem_ptr - 1;
        //for (i = 0, b = e = l_rep = 0; i < aux_->mem.n; ++i) // compute frac_rep
        do
        {
            pos ++;
            SMEM *p = &matchArray[pos];
            int sb = p->m, se = p->n + 1;
            if (p->s <= opt->max_occ) continue;
            if (sb > e) l_rep += e - b, b = sb, e = se;
            else e = e > se? e : se;
        } while (pos < num_smem - 1 && matchArray[pos].rid == matchArray[pos + 1].rid);
        l_rep += e - b;

        // bwt_sa
        // assert(pos - smem_ptr + 1 < 6000);
        if (pos - smem_ptr + 1 >= smem_buf_size)
        {
            int csize = smem_buf_size;
            smem_buf_size *= 2;
            sa_coord = (int64_t *) _mm_realloc(sa_coord, csize, opt->max_occ * smem_buf_size,
                                               sizeof(int64_t));
            assert(sa_coord != NULL);
        }

        /* Phase A buffer (Spec S4): resolved-seed records for the reorder path.
         * Grow with realloc, persist across reads. Reset nrec per read. Only the
         * opt-in reorder modes populate it; OFF streams straight to chaining. */
        int64_t nrec = 0;
        if (reorder)
        {
            /* Tight ceiling: each SMEM contributes at most min(occ, max_occ)
             * resolved seeds (the resolve loop below stops at the occurrence
             * count and caps at max_occ), so summing that is a far snugger bound
             * than the old (#SMEMs * max_occ) worst case whenever most SMEMs are
             * low-occurrence. Still a strict upper bound on nrec. */
            int64_t recs_need = 0;
            for (int64_t si = smem_ptr; si <= pos; ++si)
            {
                int64_t occ = matchArray[si].s;
                recs_need += occ < opt->max_occ ? occ : opt->max_occ;
            }
            if (recs_need > recs_cap)
            {
                recs_cap = recs_need;
                /* CodeRabbit: temp pointer so a NULL realloc doesn't leak recs. */
                seed_rec_t *recs_new = (seed_rec_t *) realloc(recs, recs_cap * sizeof(seed_rec_t));
                assert(recs_new != NULL);
                recs = recs_new;
            }
        }

        int64_t id = 0, cnt_ = 0, mypos = 0;
        #if SA_COMPRESSION
        uint64_t tim = __rdtsc();
        fmi->get_sa_entries_prefetch(&matchArray[smem_ptr], sa_coord, &cnt_,
                                     pos - smem_ptr + 1, opt->max_occ, tid, id);  // sa compressed prefetch
        tprof[MEM_SA][tid] += __rdtsc() - tim;
        #endif

        for (i = smem_ptr; i <= pos; i++)
        {
            SMEM *p = &matchArray[i];
            int64_t step;
            int32_t count, slen = p->n + 1 - p->m; // seed length
            int64_t k;
            step = p->s > opt->max_occ? p->s / opt->max_occ : 1;

            int cnt = 0;
            #if !SA_COMPRESSION
            uint64_t tim = __rdtsc();
            fmi->get_sa_entries(p, sa_coord, &cnt, 1, opt->max_occ);
            tprof[MEM_SA][tid] += __rdtsc() - tim;
            #endif

            cnt = 0;
            for (k = count = 0; k < p->s && count < opt->max_occ; k += step, ++count)
            {
                mem_seed_t s;
                int rid;

                /* Phase A (resolve). B2/B3: the SA cursor advances BEFORE any drop
                 * (rid<0 / un-remappable meth) — a dropped hit still consumes a
                 * slot. The #if SA_COMPRESSION / #else regimes index different
                 * cursors (global mypos vs per-SMEM cnt); keep both. */
                #if SA_COMPRESSION
                s.rbeg = sa_coord[mypos++];
                #else
                s.rbeg = sa_coord[cnt++];
                #endif

                s.qbeg = p->m;
                s.score= s.len = slen;
                // Propagate SMEM SA-count so chain_n_hits gates --supp-rep-hard-cap.
                s.n_hits = static_cast<int32_t>(p->s);
                if (s.rbeg < 0 || s.len < 0)
                    fprintf(stderr, "rbeg: %lld, slen: %d, cnt: %d, n: %d, m: %d, num_smem: %lld\n",
                            (long long)s.rbeg, s.len, cnt-1, p->n, p->m, (long long)num_smem);

                /* D3 (--meth, PR-3) seed→original remap (spec §5.2). The SA coord
                 * just decoded (s.rbeg) is in SEED (f/r-doubled) coordinates. Map
                 * it to ORIGINAL doubled-pac coords + an OT/OB hypothesis so all
                 * downstream chaining/extension/output sees original space only.
                 * Non-meth: chain_bns == bns and rid is the seed rid as before. */
                int8_t meth_hyp = -1;
                if (meth_remap) {
                    int64_t o_rbeg = 0; int o_rid = -1;
                    if (meth_seed_to_orig(bns, meth_orig_bns, s.rbeg, s.len,
                                          &o_rbeg, &o_rid, &meth_hyp) != 0)
                        continue; /* un-remappable (bridge/range): drop the seed */
                    s.rbeg = o_rbeg;
                    rid = o_rid;
                } else {
                    rid = bns_intv2rid(bns, s.rbeg, s.rbeg + s.len);
                }
                // bridging multiple reference sequences or the
                // forward-reverse boundary; TODO: split the seed;
                // don't discard it!!!
                if (rid < 0) continue;

                if (!reorder)
                {
                    /* OFF (default): chain the seed immediately, in resolve
                     * order — the single-pass streaming path. No recs[] write,
                     * no order_seeds; byte-identical to buffering with the
                     * identity order but without the double memory traffic. */
                    chain_add_one_seed(opt, l_pac, chain_bns, tree, seedBuf,
                                       &seedBufCount, seedBufSize, tid, l, &num[l],
                                       &s, rid, meth_hyp);
                }
                else
                {
                    // Phase A: append the fully-resolved seed (drops above already
                    // consumed the SA slot). B1: copy the FULL mem_seed_t.
                    seed_rec_t *r = &recs[nrec++];
                    r->seed = s;
                    r->rid = rid;
                    r->meth_hyp = meth_hyp; // -1 when !meth_remap
                    r->orig_ix = (uint32_t)(nrec - 1);
                }
            }
        } // seeds

        if (reorder)
        {
            // Phase B (order). Stable on orig_ix.
            order_seeds(recs, nrec, opt->seed_emit_order);

            // Phase C (chain). Replay the ordered buffer through the shared
            // chaining helper. S5: equal-pos insertion order into the kbtree is
            // preserved (no dedup/compact beyond order_seeds).
            for (int64_t ri = 0; ri < nrec; ++ri)
                chain_add_one_seed(opt, l_pac, chain_bns, tree, seedBuf,
                                   &seedBufCount, seedBufSize, tid, l, &num[l],
                                   &recs[ri].seed, recs[ri].rid, recs[ri].meth_hyp);
        } // reorder

        smem_ptr = pos + 1;
        size = kb_size(tree);
        // tprof[PE21][0] += kb_size(tree) * sizeof(mem_chain_t);

        kv_resize(mem_chain_t, *chain, size);

#define traverse_func(p_) (chain->a[chain->n++] = *(p_))
        __kb_traverse(mem_chain_t, tree, traverse_func);
#undef traverse_func

        for (i = 0; i < chain->n; ++i)
            chain->a[i].frac_rep = (float)l_rep / seq_[l].l_seq;

        kb_destroy(chn, tree);

    } // iterations over input reads
    tprof[MEM_SA_BLOCK][tid] += __rdtsc() - tim;

    _mm_free(sa_coord);
    free(recs);
}

int mem_kernel1_core(FMI_search *fmi,
                     const mem_opt_t *opt,
                     bseq1_t *seq_,
                     int nseq,
                     mem_chain_v *chain_ar,
                     mem_seed_t *seedBuf,
                     int64_t seedBufSize,
                     mem_cache *mmc,
                     int tid,
                     const bntseq_t *meth_orig_bns,
                     const uint8_t  *meth_orig_pac)
{
    /* D3 (--meth, PR-3): in --meth chaining/filtering run in ORIGINAL coords;
     * outside --meth (or if a handle is missing) fall back to the seed index so
     * the non-meth path is byte-for-byte unchanged. */
    const bntseq_t *chain_bns = (opt->meth_mode && meth_orig_bns != NULL)
                                ? meth_orig_bns : fmi->idx->bns;
    const uint8_t  *chain_pac = (opt->meth_mode && meth_orig_pac != NULL)
                                ? meth_orig_pac : fmi->idx->pac;
    int i;
    int64_t num_smem = 0, tot_len = 0;
    mem_chain_v *chn;

    uint64_t tim;
    /* convert to 2-bit encoding if we have not done so */
    for (int l=0; l<nseq; l++)
    {
        char *seq = seq_[l].seq;
        int len = seq_[l].l_seq;
        tot_len += len;

        for (i = 0; i < len; ++i)
            seq[i] = seq[i] < 4? seq[i] : nst_nt4_table[(int)seq[i]]; //nst_nt4??
    }
    // Empty-batch fast path: a 0-base batch (e.g. malformed FASTQ parsed as a
    // single zero-length record) tripped the lazy-init grow block below with
    // tot_len(0) >= wsize_mem(0), routing through _mm_realloc(NULL, 0, 0, …)
    // and exiting via the post-collect overflow check. With no bases, there
    // are no SMEMs to collect — initialize chain_ar to empty so mem_kernel2_core
    // (which walks chain_ar unconditionally and free()s chain_ar[l].a) sees
    // a valid empty state, then let downstream stages emit unmapped records.
    if (tot_len == 0) {
        for (int l = 0; l < nseq; ++l)
            kv_init(chain_ar[l]);
        return 1;
    }
    // enc_qdb is sized in *bytes* (= tot_len). Track its capacity via the
    // dedicated wsize_qdb so it stays correct independent of wsize_mem,
    // which mem_collect_smem grows in *SMEM-entry* units. On short-read
    // workloads (≤50 bp aDNA), the post-mem_collect_smem wsize_mem in SMEM
    // units can exceed tot_len in bytes, which would make the legacy
    // `if (tot_len >= wsize_mem)` gate skip the enc_qdb realloc on
    // subsequent batches even when enc_qdb is undersized.
    if (tot_len > mmc->wsize_qdb[tid]) {
        int64_t tmp = mmc->wsize_qdb[tid];
        uint8_t *new_enc_qdb = (uint8_t *) realloc(mmc->enc_qdb[tid],
                                                   tot_len * sizeof(uint8_t));
        if (new_enc_qdb == NULL) err_fatal(__func__, "out of memory growing enc_qdb");
        mmc->enc_qdb[tid] = new_enc_qdb;
        mmc->wsize_qdb[tid] = tot_len;
        if (bwa_verbose >= 4) {
            fprintf(stderr, "[%0.4d] Re-allocating enc_qdb: "
                    "%" PRId64 " -> %" PRId64 "\n",
                    tid, tmp, mmc->wsize_qdb[tid]);
        }
    }
    // tot_len *= N_SMEM_KERNEL;
    // fprintf(stderr, "wsize: %d, tot_len: %d\n", mmc->wsize_mem[tid], tot_len);
    // This covers matchArray and the per-read int arrays. enc_qdb is grown
    // separately above (wsize_qdb).
    if (tot_len >= mmc->wsize_mem[tid])
    {
        int64_t tmp = mmc->wsize_mem[tid];
        mmc->wsize_mem[tid] = tot_len;
        if (bwa_verbose >= 4) {
            fprintf(stderr, "[%0.4d] Re-allocating SMEM scratch: "
                    "%" PRId64 " -> %" PRId64 "\n",
                    tid, tmp, mmc->wsize_mem[tid]);
        }
        mmc->wsize_mem_s[tid] = tot_len;
        mmc->wsize_mem_r[tid] = tot_len;
        mmc->matchArray[tid]   = (SMEM *) _mm_realloc(mmc->matchArray[tid],
                                                      tmp, mmc->wsize_mem[tid], sizeof(SMEM));
            //realloc(mmc->matchArray[tid], mmc->wsize_mem[tid] *   sizeof(SMEM));
        int32_t *new_min_intv_ar = (int32_t *) realloc(mmc->min_intv_ar[tid],
                                                       mmc->wsize_mem[tid] * sizeof(int32_t));
        if (new_min_intv_ar == NULL) err_fatal(__func__, "out of memory growing min_intv_ar");
        mmc->min_intv_ar[tid] = new_min_intv_ar;

        int32_t *new_query_pos_ar = (int32_t *) realloc(mmc->query_pos_ar[tid],
                                                        mmc->wsize_mem[tid] * sizeof(int32_t));
        if (new_query_pos_ar == NULL) err_fatal(__func__, "out of memory growing query_pos_ar");
        mmc->query_pos_ar[tid] = new_query_pos_ar;

        int32_t *new_rid = (int32_t *) realloc(mmc->rid[tid],
                                               mmc->wsize_mem[tid] * sizeof(int32_t));
        if (new_rid == NULL) err_fatal(__func__, "out of memory growing rid");
        mmc->rid[tid] = new_rid;
        // w.mmc.lim[l]        = (int32_t *) _mm_malloc((BATCH_SIZE + 32) * sizeof(int32_t), 64);
    }

    SMEM    *matchArray   = mmc->matchArray[tid];
    int32_t *min_intv_ar  = mmc->min_intv_ar[tid];
    int32_t *query_pos_ar = mmc->query_pos_ar[tid];
    uint8_t *enc_qdb      = mmc->enc_qdb[tid];
    int32_t *rid          = mmc->rid[tid];
    int64_t  *wsize_mem   = &mmc->wsize_mem[tid];

    tim = __rdtsc();
    /********************** Kernel 1: FM+SMEMs *************************/
    printf_(VER, "6. Calling mem_collect_smem.., tid: %d\n", tid);
    matchArray = mem_collect_smem(fmi, opt,
                                  seq_,
                                  nseq,
                                  matchArray,
                                  min_intv_ar,
                                  query_pos_ar,
                                  enc_qdb,
                                  rid,
                                  mmc,
                                  num_smem,
                                  tid);

    // Exact-fit (num_smem == *wsize_mem) is valid: the last write lands at
    // matchArray[*wsize_mem - 1]. Trip only on a true overrun.
    if (num_smem > *wsize_mem){
        fprintf(stderr, "Error [bug]: num_smem: %lld are more than allocated space %lld.\n",
                (long long)num_smem, (long long)*wsize_mem);
        exit(EXIT_FAILURE);
    }
    printf_(VER, "6. Done! mem_collect_smem, num_smem: %lld\n", (long long)num_smem);
    tprof[MEM_COLLECT][tid] += __rdtsc() - tim;


    /********************* Kernel 1.1: SA2REF **********************/
    tim = __rdtsc();
    printf_(VER, "6.1. Calling mem_chain..\n");
    mem_chain_seeds(fmi, opt, fmi->idx->bns,
                    /* remap target (NULL outside --meth) */ meth_orig_bns,
                    seq_, nseq, tid,
                    chain_ar,
                    seedBuf,
                    seedBufSize,
                    matchArray,
                    num_smem);

    printf_(VER, "5. Done mem_chain..\n");
    tprof[MEM_CHAIN][tid] += __rdtsc() - tim;

    /************** Post-processing of collected smems/chains ************/
    // tim = __rdtsc();
    printf_(VER, "6.1. Calling mem_chain_flt..\n");
    /* Pass 1: filter every read first, so the mate's filtered chains are
     * available to the mate-concordance-aware cap (--extend-mate-concordant) below. */
    for (int l=0; l<nseq; l++) {
        chn = &chain_ar[l];
        chn->n = mem_chain_flt(opt, chn->n, chn->a, tid);
    }
    int cap_mate_concordant = opt->mate_concordant_window != 0 && (opt->flag & MEM_F_PE)
                              && opt->max_extend_chains > 0;
    /* Resolve the concordance window: fixed CLI value if >0; else (auto == -1) the
     * estimated proper-pair insert high bound if known, else the fallback. */
    int64_t mate_win = 0;
    if (cap_mate_concordant)
        mate_win = opt->mate_concordant_window > 0 ? opt->mate_concordant_window
                 : (opt->est_insert_high > 0 ? opt->est_insert_high
                                             : MATE_CONCORDANT_WINDOW_FALLBACK);
    /* Pass 2: cap chains (mate-concordance-aware when --extend-mate-concordant). */
    for (int l=0; l<nseq; l++)
    {
        chn = &chain_ar[l];
        if (cap_mate_concordant && (l ^ 1) < nseq) {
            int ml = l ^ 1;  /* interleaved PE: mate is the sibling index */
            chn->n = mem_chain_cap_extend_mate(chn->a, chn->n, opt->max_extend_chains,
                                               chain_ar[ml].a, chain_ar[ml].n,
                                               chain_bns, mate_win);
        } else {
            chn->n = mem_chain_cap_extend(chn->a, chn->n, opt->max_extend_chains);
        }
    }
    printf_(VER, "7. Done mem_chain_flt..\n");
    // tprof[MEM_ALN_M1][tid] += __rdtsc() - tim;


    printf_(VER, "8. Calling mem_flt_chained_seeds..\n");
    for (int l=0; l<nseq; l++) {
        chn = &chain_ar[l];
        /* D3: seeds are now in ORIGINAL coords (remapped in mem_chain_seeds), so
         * the seed-SW filter fetches the ORIGINAL ref window.
         * D3 (--meth, PR-4): mem_flt_chained_seeds now scores the ORIGINAL read
         * (meth_orig_seq) against the original ref with the per-chain asymmetric
         * OT/OB matrix (mem_opt_meth_mat), so bisulfite conversions are no longer
         * mis-penalized at the seed-filter stage. */
        mem_flt_chained_seeds(opt, chain_bns, chain_pac, seq_, chn->n, chn->a);
    }
    printf_(VER, "8. Done mem_flt_chained_seeds..\n");
    // tprof[MEM_ALN_M2][tid] += __rdtsc() - tim;


    return 1;
}

int mem_kernel2_core(FMI_search *fmi,
                     const mem_opt_t *opt,
                     bseq1_t *seq_,
                     mem_alnreg_v *regs,
                     int nseq,
                     mem_chain_v *chain_ar,
                     mem_cache *mmc,
                     uint8_t *ref_string,
                     int tid,
                     const bntseq_t *meth_orig_bns,
                     const uint8_t  *meth_orig_pac)
{
    int i;
    /* D3 (--meth, PR-3): extension/dedup must run in ORIGINAL coords. ref_string
     * is already the original .0123 (worker_aln passes mem_aln_ref_string).
     * Outside --meth these fall back to the seed/normal index unchanged. */
    const bntseq_t *aln_bns = (opt->meth_mode && meth_orig_bns != NULL)
                              ? meth_orig_bns : fmi->idx->bns;
    const uint8_t  *aln_pac = (opt->meth_mode && meth_orig_pac != NULL)
                              ? meth_orig_pac : fmi->idx->pac;
    for (int l=0; l<nseq; l++)
    {
        kv_init(regs[l]);
    }
    /****************** Kernel 2: B-SWA *********************/
    uint64_t tim = __rdtsc();
    printf_(VER, "9. Calling mem_chain2aln...\n");
    /* D3 (--meth, PR-4): mem_chain2aln_across_reads_V2 now extends the ORIGINAL
     * read bases (seq_[].meth_orig_seq, 2-bit-encoded internally) against the
     * ORIGINAL ref window with the per-hypothesis asymmetric matrix selected from
     * chain->meth_hypothesis. The ungapped fast path is disabled under --meth
     * (it cannot express the matrix). a->score / a->truesc — and therefore
     * selection (mem_mark_primary_se) and MAPQ (mem_approx_mapq_se) downstream —
     * are now original-space γ scores. A MIXED-hypothesis batch (interleaved
     * R1/OT + R2/OB — the common paired-end case) is handled by bsw_run_tier's
     * per-hypothesis partition pass inside the function: each tier's pairs are
     * split by meth_hypothesis and scored against a matching OT/OB SW object, so
     * every pair gets its asymmetric matrix. Homogeneous (single-hypothesis)
     * batches collapse to a single partition with no extra cost. */
    mem_chain2aln_across_reads_V2(opt,
                                  aln_bns,
                                  aln_pac,
                                  seq_,
                                  nseq,
                                  chain_ar,
                                  regs,
                                  mmc,
                                  ref_string,
                                  tid);

    printf_(VER, "9. Done mem_chain2aln...\n\n");
    tprof[MEM_ALN2][tid] += __rdtsc() - tim;

    // tim = __rdtsc();
    for (int l=0; l<nseq; l++) {
        mem_chain_v *chain = &chain_ar[l];
        for (int i = 0; i < chain->n; ++i)
        {
            mem_chain_t chn = chain->a[i];
            if (chn.m > SEEDS_PER_CHAIN)
            {
                tprof[PE11][tid] ++;
                free(chn.seeds);
            }
            tprof[PE12][tid]++;
        }
        free(chain_ar[l].a);
    }

    int m = 0;
    for (int l=0; l<nseq; l++)
    {
        mem_alnreg_t *a = regs[l].a;
        int n = regs[l].n;
        for (i = 0, m = 0; i < n; ++i) // exclude identical hits
            if (a[i].qe > a[i].qb) {
                if (m != i) a[m++] = a[i];
                else ++m;
            }
        regs[l].n = m;
    }

    /* D3 (--meth, PR-4): the colinear alnreg merge in mem_sort_dedup_patch
     * recomputes p->score/p->truesc (feeds selection + MAPQ), so it must use the
     * ORIGINAL read + per-read OT/OB matrix like extension. Under directional
     * --meth every alnreg of a read shares one hypothesis, so one matrix per read
     * is correct. Outside --meth (or with no orig seq) it stays the projected
     * read + opt->mat — unchanged. The original read is 2-bit-encoded into a
     * reusable scratch buffer (same orientation as seq, per bwa.h). */
    uint8_t *dd_qbuf = NULL;
    int      dd_qbuf_cap = 0;
    for (int l=0; l<nseq; l++) {
        uint8_t      *dd_query = (uint8_t*) seq_[l].seq;
        const int8_t *dd_mat   = opt->mat;
        if (opt->meth_mode && seq_[l].meth_orig_seq != NULL) {
            int l_query = seq_[l].l_seq;
            if (l_query > dd_qbuf_cap) {
                dd_qbuf_cap = l_query;
                dd_qbuf = (uint8_t *) realloc(dd_qbuf, (size_t) dd_qbuf_cap);
                assert(dd_qbuf != NULL);
            }
            const char *os = seq_[l].meth_orig_seq;
            for (int q = 0; q < l_query; ++q) {
                unsigned char ch = (unsigned char) os[q];
                dd_qbuf[q] = (ch < 4) ? ch : nst_nt4_table[ch];
            }
            dd_query = dd_qbuf;
            /* per-read STRAND-ADJUSTED hypothesis (see meth_strand_hyp); take the
             * first alnreg's (dedup-patch shares one matrix across the read). */
            int hyp = (regs[l].n > 0) ? regs[l].a[0].meth_strand_hyp : -1;
            dd_mat = mem_opt_meth_mat(opt, hyp);
        }
        regs[l].n = mem_sort_dedup_patch(opt, aln_bns,
                                         aln_pac,
                                         dd_query,
                                         regs[l].n, regs[l].a, dd_mat);
    }
    free(dd_qbuf);  /* NULL outside --meth; free() is NULL-safe */

    for (int l=0; l<nseq; l++)
    {
        for (i = 0; i < regs[l].n; ++i)
        {
            mem_alnreg_t *p = &regs[l].a[i];
            if (p->rid >= 0 && aln_bns->anns[p->rid].is_alt)
                p->is_alt = 1;
        }
    }
    // tprof[POST_SWA][tid] += __rdtsc() - tim;

    return 1;
}

static void worker_aln(void *data, int seq_id, int batch_size, int tid)
{
    worker_t *w = (worker_t*) data;

    printf_(VER, "11. Calling mem_kernel2_core..\n");
    /* D3 (--meth): pass the ORIGINAL .0123 ref_string and original bns/pac so
     * extension/dedup fetch original bases for the (now original-coord) seeds.
     * mem_aln_* return the seed/normal handles outside --meth (no-op there). */
    mem_kernel2_core(w->fmi, w->opt,
                     w->seqs + seq_id,
                     w->regs + seq_id,
                     batch_size,
                     w->chain_ar + seq_id,
                     &w->mmc,
                     mem_aln_ref_string(w),
                     tid,
                     w->meth_orig_bns, w->meth_orig_pac);
    printf_(VER, "11. Done mem_kernel2_core....\n");

}

/* Kernel, called by threads */
static void worker_bwt(void *data, int seq_id, int batch_size, int tid)
{
    worker_t *w = (worker_t*) data;
    printf_(VER, "4. Calling mem_kernel1_core..%d %d\n", seq_id, tid);
    int seedBufSz = w->seedBufSize;

    int memSize = w->nreads;
    if (batch_size < BATCH_SIZE) {
        seedBufSz = (memSize - seq_id) * AVG_SEEDS_PER_READ;
        // fprintf(stderr, "[%0.4d] Info: adjusted seedBufSz %d\n", tid, seedBufSz);
    }

    mem_kernel1_core(w->fmi, w->opt,
                     w->seqs + seq_id,
                     batch_size,
                     w->chain_ar + seq_id,
                     w->seedBuf + seq_id * AVG_SEEDS_PER_READ,
                     seedBufSz,
                     &(w->mmc),
                     tid,
                     /* D3 (--meth) remap target; NULL outside --meth */
                     w->meth_orig_bns, w->meth_orig_pac);
    printf_(VER, "4. Done mem_kernel1_core....\n");
}

int64_t sort_classify(mem_cache *mmc, int64_t pcnt, int tid)
{

    SeqPair *seqPairArray = mmc->seqPairArrayLeft128[tid];
    // SeqPair *seqPairArrayAux = mmc->seqPairArrayAux[tid];
    SeqPair *seqPairArrayAux = mmc->seqPairArrayRight128[tid];

    int64_t pos8 = 0, pos16 = 0;
    for (int i=0; i<pcnt; i++)
    {
        SeqPair *s = seqPairArray + i;
        int xtra = s->h0;
        int size = (xtra & KSW_XBYTE)? 1 : 2;
        if (size == 1) // 8
        {
            seqPairArray[pos8++] = seqPairArray[i];
        } else { // 16
            seqPairArrayAux[pos16++] = seqPairArray[i];
        }
    }
    assert(pos8 + pos16 == pcnt);

    for (int i=pos8; i<pcnt; i++) {
        seqPairArray[i] = seqPairArrayAux[i-pos8];
    }

    return pos8;
}

static void worker_sam(void *data, int seqid, int batch_size, int tid)
{
    worker_t *w = (worker_t*) data;

    if (w->opt->flag & MEM_F_PE)
    {
        int64_t pcnt = 0;
        int start = seqid;
        int end = seqid + batch_size;
        int pos = start >> 1;

#if !BWAMEM_BATCHED_MATESW
        // Scalar mem_sam_pe path. Selected when no SIMD kswv kernel is
        // available (e.g. plain SSE2/AVX2 x86 builds, or forced via
        // DISABLE_BATCHED_MATESW=1 for the proto-neon-kswv CI A/B).
        for (int i=start; i< end; i+=2)
        {
            // orig mem_sam_pe() function
            /* D3 (--meth, PR-3): pairing + mate rescue run in ORIGINAL coords
             * (anchors already remapped). mem_aln_bns/pac return the original
             * bns/pac in --meth; the seed/normal index otherwise. This is where
             * the mate-rescue coordinate reconciliation (6a) lands: mem_matesw's
             * window math derives l_pac from bns->l_pac, so passing original
             * bns/pac fixes the doubled-pac↔real-chrom mismatch.
             * D3 (--meth, PR-6): mem_pair_resolve now drives mem_matesw with the
             * ORIGINAL mate bases (s[!i].meth_orig_seq) + the OPPOSITE-strand
             * asymmetric matrix of the anchor; the rescued mate's meth_hypothesis
             * is set to !anchor. (This scalar mem_sam_pe path is selected only
             * when no batched kswv kernel is available.) */
            mem_sam_pe(w->opt, mem_aln_bns(w),
                       mem_aln_pac(w), w->pes,
                       (w->n_processed >> 1) + pos++,   // check!
                       &w->seqs[i],
                       &w->regs[i]);

            free(w->regs[i].a);
            free(w->regs[i+1].a);
        }
#else   // re-structured
        // pre-processing
        // uint64_t tim = __rdtsc();
        int32_t maxRefLen = 0, maxQerLen = 0;
        int32_t gcnt = 0;
        for (int i=start; i< end; i+=2)
        {
            mem_sam_pe_batch_pre(w->opt, mem_aln_bns(w),
                                 mem_aln_pac(w), w->pes,
                                 (w->n_processed >> 1) + pos++,   // check!
                                 &w->seqs[i],
                                 &w->regs[i],
                                 &w->mmc,
                                 pcnt, gcnt,
                                 maxRefLen,
                                 maxQerLen,
                                 tid);
        }

        // tprof[SAM1][tid] += __rdtsc() - tim;
        int64_t pcnt8 = sort_classify(&w->mmc, pcnt, tid);

        kswr_t *aln = (kswr_t *) _mm_malloc ((pcnt + SIMD_WIDTH8) * sizeof(kswr_t), 64);
        assert(aln != NULL);

        // processing
        mem_sam_pe_batch(w->opt, &w->mmc, pcnt, pcnt8, aln, maxRefLen, maxQerLen, tid);

        // post-processing
        // tim = __rdtsc();
        gcnt = 0;
        pos = start >> 1;
        kswr_t *myaln = aln;
        for (int i=start; i< end; i+=2)
        {
            mem_sam_pe_batch_post(w->opt, mem_aln_bns(w),
                                  mem_aln_pac(w), w->pes,
                                  (w->n_processed >> 1) + pos++,   // check!
                                  &w->seqs[i],
                                  &w->regs[i],
                                  &myaln,
                                  &w->mmc,
                                  gcnt,
                                  tid);

            free(w->regs[i].a);
            free(w->regs[i+1].a);
        }
        //tprof[SAM3][tid] += __rdtsc() - tim;
        _mm_free(aln);  // kswr_t
#endif
    }
    else
    {
        for (int i=seqid; i<seqid + batch_size; i++)
        {
            mem_mark_primary_se(w->opt, w->regs[i].n,
                                w->regs[i].a,
                                w->n_processed + i);
#if V17  // Feature from v0.7.17 of orig. bwa-mem
            if (w->opt->flag & MEM_F_PRIMARY5) mem_reorder_primary5(w->opt->T, &w->regs[i]);
#endif
            /* D3 (--meth): emit in ORIGINAL coords/alphabet via mem_aln_bns/pac. */
            mem_reg2sam(w->opt, mem_aln_bns(w), mem_aln_pac(w), &w->seqs[i],
                        &w->regs[i], 0, 0);
            free(w->regs[i].a);
        }
    }
}

void mem_process_seqs(mem_opt_t *opt,
                      int64_t n_processed,
                      int n,
                      bseq1_t *seqs,
                      const mem_pestat_t *pes0,
                      worker_t &w)
{
    mem_pestat_t pes[4];
    double ctime, rtime;

    ctime = cputime(); rtime = realtime();
    w.opt = opt;
    w.seqs = seqs; w.n_processed = n_processed;
    w.pes = &pes[0];

    /* When -I supplied the insert-size distribution (pes0), publish its FR high
     * bound onto opt BEFORE worker_bwt so this chunk's mate-concordant chain cap
     * (--extend-mate-concordant auto) honors it. Otherwise the bound is only set
     * post-pairing at the end of this function, so a single-chunk -I run would
     * cap with the fallback window and never use the user-provided insert size.
     * The inferred-from-data path (pes0 == NULL) still publishes below for the
     * next chunk. pes index 1 = FR orientation. */
    if ((opt->flag & MEM_F_PE) && opt->mate_concordant_window == -1 &&
        pes0 && !pes0[1].failed && pes0[1].high > 0)
        opt->est_insert_high = pes0[1].high;

    //int n_ = (opt->flag & MEM_F_PE) ? n : n;   // this requires n%2==0
    int n_ = n;

    uint64_t tim = __rdtsc();
    fprintf(stderr, "[0000] 1. Calling kt_for - worker_bwt\n");

    kt_for(worker_bwt, &w, n_); // SMEMs (+SAL)

    fprintf(stderr, "[0000] 2. Calling kt_for - worker_aln\n");

    kt_for(worker_aln, &w, n_); // BSW
    tprof[WORKER10][0] += __rdtsc() - tim;


    // PAIRED_END
    if (opt->flag & MEM_F_PE) { // infer insert sizes if not provided
        if (pes0)
            memcpy(pes, pes0, 4 * sizeof(mem_pestat_t)); // if pes0 != NULL, set the insert-size
                                                         // distribution as pes0
        else {
            /* D3 (--meth): insert-size inference operates on alnreg.rb which are
             * now in ORIGINAL doubled-pac space, so use the ORIGINAL l_pac. */
            const bntseq_t *pestat_bns = mem_aln_bns(&w);
            fprintf(stderr, "[0000] Inferring insert size distribution of PE reads from data, "
                    "l_pac: %lld, n: %d\n", (long long)pestat_bns->l_pac, n);
            mem_pestat(opt, pestat_bns->l_pac, n, w.regs, pes); // otherwise, infer the insert size
                                                         // distribution from data
        }
        /* Publish the FR proper-pair insert high bound onto opt so the NEXT
         * chunk's mate-concordant chain cap (which runs before pairing) can use
         * it as its window (--extend-mate-concordant auto). opt is shared and
         * mutable across chunks; this write is single-threaded (between kt_for
         * passes), and the cap only reads it. pes index 1 = FR orientation. */
        if (!pes[1].failed && pes[1].high > 0)
            opt->est_insert_high = pes[1].high;
    }

    tim = __rdtsc();
    fprintf(stderr, "[0000] 3. Calling kt_for - worker_sam\n");

    kt_for(worker_sam, &w,  n_);   // SAM
    tprof[WORKER20][0] += __rdtsc() - tim;

    fprintf(stderr, "\t[0000][ M::%s] Processed %d reads in %.3f "
            "CPU sec, %.3f real sec\n",
            __func__, n, cputime() - ctime, realtime() - rtime);

}

static void mem_mark_primary_se_core(const mem_opt_t *opt, int n, mem_alnreg_t *a, int_v *z)
{ // similar to the loop in mem_chain_flt()
    int i, k, tmp;
    tmp = opt->a + opt->b;
    tmp = opt->o_del + opt->e_del > tmp? opt->o_del + opt->e_del : tmp;
    tmp = opt->o_ins + opt->e_ins > tmp? opt->o_ins + opt->e_ins : tmp;
    z->n = 0;
    kv_push(int, *z, 0);
    for (i = 1; i < n; ++i) {
        for (k = 0; k < z->n; ++k) {
            int j = z->a[k];
            int b_max = a[j].qb > a[i].qb? a[j].qb : a[i].qb;
            int e_min = a[j].qe < a[i].qe? a[j].qe : a[i].qe;
            if (e_min > b_max) { // have overlap
                int min_l = a[i].qe - a[i].qb < a[j].qe - a[j].qb? a[i].qe - a[i].qb : a[j].qe - a[j].qb;
                if (e_min - b_max >= min_l * opt->mask_level) { // significant overlap
                    if (a[j].sub == 0) a[j].sub = a[i].score;
                    if (a[j].score - a[i].score <= tmp && (a[j].is_alt || !a[i].is_alt))
                        ++a[j].sub_n;
                    break;
                }
            }
        }
        if (k == z->n) kv_push(int, *z, i);
        else a[i].secondary = z->a[k];
    }
}

int mem_mark_primary_se(const mem_opt_t *opt, int n, mem_alnreg_t *a, int64_t id)
{
    int i, n_pri;
    int_v z = {0,0,0};
    if (n == 0) return 0;

    /* Fast path for the dominant unique-mapper case. With a single region the
     * general path below still runs two ks_introsort calls and a kv_push in
     * mem_mark_primary_se_core (which heap-allocates z.a) only to leave a[0]
     * with these exact values: mem_mark_primary_se_core touches nothing for
     * n==1 (its inner loop starts at i=1), the first rewrite loop sets
     * secondary_all to the rank 0 then normalizes it to secondary (-1), and
     * sub_n is never touched here. Reproduce that directly. */
    if (n == 1) {
        a[0].sub = a[0].alt_sc = 0;
        a[0].secondary = a[0].secondary_all = -1;
        a[0].hash = hash_64(id);
        return a[0].is_alt ? 0 : 1;
    }

    for (i = n_pri = 0; i < n; ++i)
    {
        a[i].sub = a[i].alt_sc = 0, a[i].secondary = a[i].secondary_all = -1, a[i].hash = hash_64(id+i);
        if (!a[i].is_alt) ++n_pri;
    }
    ks_introsort(mem_ars_hash, n, a);
    mem_mark_primary_se_core(opt, n, a, &z);
    for (i = 0; i < n; ++i)
    {
        mem_alnreg_t *p = &a[i];
        p->secondary_all = i; // keep the rank in the first round
        if (!p->is_alt && p->secondary >= 0 && a[p->secondary].is_alt)
            p->alt_sc = a[p->secondary].score;
    }
    if (n_pri >= 0 && n_pri < n)
    {
        kv_resize(int, z, n);
        if (n_pri > 0) ks_introsort(mem_ars_hash2, n, a);
        for (i = 0; i < n; ++i) z.a[a[i].secondary_all] = i;
        for (i = 0; i < n; ++i)
        {
            if (a[i].secondary >= 0)
            {
                a[i].secondary_all = z.a[a[i].secondary];
                if (a[i].is_alt) a[i].secondary = INT_MAX;
            } else a[i].secondary_all = -1;
        }
        if (n_pri > 0) { // mark primary for hits to the primary assembly only
            for (i = 0; i < n_pri; ++i) a[i].sub = 0, a[i].secondary = -1;
            mem_mark_primary_se_core(opt, n_pri, a, &z);
        }
    }
    else {
        for (i = 0; i < n; ++i)
            a[i].secondary_all = a[i].secondary;
    }
    free(z.a);
    return n_pri;
}

/************************
 * Integrated interface *
 ************************/

int mem_approx_mapq_se(const mem_opt_t *opt, const mem_alnreg_t *a)
{
    int mapq, l, sub = a->sub? a->sub : opt->min_seed_len * opt->a;
    double identity;
    sub = a->csub > sub? a->csub : sub;
    if (sub >= a->score) return 0;
    l = a->qe - a->qb > a->re - a->rb? a->qe - a->qb : a->re - a->rb;
    identity = 1. - (double)(l * opt->a - a->score) / (opt->a + opt->b) / l;
    if (a->score == 0) {
        mapq = 0;
    } else if (opt->mapQ_coef_len > 0) {
        double tmp;
        tmp = l < opt->mapQ_coef_len? 1. : opt->mapQ_coef_fac / log(l);
        tmp *= identity * identity;
        mapq = (int)(6.02 * (a->score - sub) / opt->a * tmp * tmp + .499);
    } else {
        mapq = (int)(MEM_MAPQ_COEF * (1. - (double)sub / a->score) * log(a->seedcov) + .499);
        mapq = identity < 0.95? (int)(mapq * identity * identity + .499) : mapq;
    }
    if (a->sub_n > 0) mapq -= (int)(4.343 * log(a->sub_n+1) + .499);
    if (mapq > 60) mapq = 60;
    if (mapq < 0) mapq = 0;
    mapq = (int)(mapq * (1. - a->frac_rep) + .499);
    return mapq;
}

void mem_reorder_primary5(int T, mem_alnreg_v *a)
{
    int k, n_pri = 0, left_st = INT_MAX, left_k = -1;
    mem_alnreg_t t;
    for (k = 0; k < a->n; ++k)
        if (a->a[k].secondary < 0 && !a->a[k].is_alt && a->a[k].score >= T) ++n_pri;
    if (n_pri <= 1) return; // only one alignment
    for (k = 0; k < a->n; ++k) {
        mem_alnreg_t *p = &a->a[k];
        if (p->secondary >= 0 || p->is_alt || p->score < T) continue;
        if (p->qb < left_st) left_st = p->qb, left_k = k;
    }
    assert(a->a[0].secondary < 0);
    if (left_k == 0) return; // no need to reorder
    t = a->a[0], a->a[0] = a->a[left_k], a->a[left_k] = t;
    for (k = 1; k < a->n; ++k) { // update secondary and secondary_all
        mem_alnreg_t *p = &a->a[k];
        if (p->secondary == 0) p->secondary = left_k;
        else if (p->secondary == left_k) p->secondary = 0;
        if (p->secondary_all == 0) p->secondary_all = left_k;
        else if (p->secondary_all == left_k) p->secondary_all = 0;
    }
}

// TODO (future plan): group hits into a uint64_t[] array. This will be cleaner and more flexible
void mem_reg2sam(const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac,
                 bseq1_t *s, mem_alnreg_v *a, int extra_flag, const mem_aln_t *m)
{
    kstring_t str;
    kvec_t(mem_aln_t) aa;
    int k, l;
    char **XA = 0;
    int *HN = 0;

    if (!(opt->flag & MEM_F_ALL))
        XA = mem_gen_alt(opt, bns, pac, a, s->l_seq, s->seq, &HN);

    kv_init(aa);
    str.l = str.m = 0; str.s = 0;
    for (k = l = 0; k < a->n; ++k)
    {
        mem_alnreg_t *p = &a->a[k];
        mem_aln_t *q;
        if (p->score < opt->T) continue;
        //fprintf(stderr, "%d %d %d %d\n", p->secondary, p->is_alt, opt->flag&MEM_F_ALL,
        //      (p->secondary >= 0 && (p->is_alt || !(opt->flag&MEM_F_ALL))));
        if (p->secondary >= 0 && (p->is_alt || !(opt->flag&MEM_F_ALL))) continue;
        // assert(p->secondary < INT_MAX);
        if (p->secondary >= 0 && p->secondary < INT_MAX && p->score < a->a[p->secondary].score * opt->drop_ratio) continue;
        q = kv_pushp(mem_aln_t, aa);

        *q = mem_reg2aln(opt, bns, pac, s->l_seq, s->seq, p, s->meth_orig_seq);
        assert(q->rid >= 0); // this should not happen with the new code
        q->XA = XA? XA[k] : 0;
        q->HN = HN? HN[k] : -1;
        q->flag |= extra_flag; // flag secondary
        if (p->secondary >= 0) q->sub = -1; // don't output sub-optimal score
        if (l && p->secondary < 0) // if supplementary
            q->flag |= (opt->flag&MEM_F_NO_MULTI)? 0x10000 : 0x800;
#if V17
        if (!(opt->flag & MEM_F_KEEP_SUPP_MAPQ) && l && !p->is_alt && q->mapq > aa.a[0].mapq)
            q->mapq = aa.a[0].mapq; // lower mapq for supplementary mappings, unless -5 or -q is applied
#else
        if (l && !p->is_alt && q->mapq > aa.a[0].mapq) q->mapq = aa.a[0].mapq;
#endif
        // fg-labs: supp alnregs whose chain contains a repetitive seed (many
        // genome occurrences) inherit primary MAPQ despite being ambiguous on
        // their own — upstream issue bwa-mem2/bwa-mem2#260. Opt-in override:
        // force MAPQ=0 when the chain's most-repetitive seed exceeds the cap.
        if (opt->supp_rep_hard_cap > 0 && l && p->secondary < 0 &&
            p->chain_n_hits >= opt->supp_rep_hard_cap)
            q->mapq = 0;
        ++l;
    }
    if (aa.n == 0) { // no alignments good enough; then write an unaligned record
        mem_aln_t t;
        t = mem_reg2aln(opt, bns, pac, s->l_seq, s->seq, 0);
        t.flag |= extra_flag;
        mem_aln2sam(opt, bns, &str, s, 1, &t, 0, m);
    } else {
        for (k = 0; k < aa.n; ++k)
            mem_aln2sam(opt, bns, &str, s, aa.n, aa.a, k, m);
        for (k = 0; k < aa.n; ++k) free(aa.a[k].cigar);
        free(aa.a);
    }
    s->sam = str.s;
    if (XA) {
        for (k = 0; k < a->n; ++k) free(XA[k]);
        free(XA);
    }
    free(HN);
}

/* Caller must have called ks_resize beforehand to leave room for
 * n_cigar*16 bytes (worst-case: each op is 11 digits + 1 letter). */
static inline void add_cigar(const mem_opt_t *opt, mem_aln_t *p, kstring_t *str, int which)
{
    int i;
    if (p->n_cigar) { // aligned
        for (i = 0; i < p->n_cigar; ++i) {
            int c = p->cigar[i]&0xf;
            if (!(opt->flag&MEM_F_SOFTCLIP) && !p->is_alt && (c == 3 || c == 4))
                c = which? 4 : 3; // use hard clipping for supplementary alignments
            kputw_u(p->cigar[i]>>4, str); kputc_u("MIDSH"[c], str);
        }
    } else kputc_u('*', str); // having a coordinate but unaligned (e.g. when copy_mate is true)
}

void mem_aln2sam(const mem_opt_t *opt, const bntseq_t *bns, kstring_t *str,
                 bseq1_t *s, int n, const mem_aln_t *list, int which, const mem_aln_t *m_)
{
    SpEncodeScope _enc;   /* stage_prof: time the SAM/BAM build (all return paths) */
    /* BAM short-circuit: meth_mode applies the bisulfite overlay (chrom
     * consolidation, YD:Z, chimera QC), plain bam_mode uses the generic
     * writer. Either path leaves `str` untouched. */
    if (opt->meth_mode) {
        struct bam1_t *b = bam_writer_alloc();
        if (b == NULL)
            err_fatal(__func__, "out of memory allocating bam1_t for meth record");
        /* `bns` here is the ORIGINAL bns under --meth (mem_aln_bns()); pair
         * it with the original pac global for the XM:Z reference window. */
        if (meth_mem_aln_to_bam(b, opt, bns, g_meth_orig_pac, s, n, list, which, m_) != 0) {
            bam_writer_free(b);
            err_fatal(__func__, "meth BAM conversion failed for read \"%s\"", s->name);
        }
        bam_writer_bseq_push(s, b);
        return;
    }
    if (opt->bam_mode) {
        if (bam_writer_push_aln(s, opt, bns, n, list, which, m_) != 0)
            err_fatal(__func__, "BAM conversion failed for read \"%s\"", s->name);
        return;
    }

    int i, l_name;
    mem_aln_t ptmp = list[which], *p = &ptmp, mtmp, *m = 0; // make a copy of the alignment to convert

    if (m_) mtmp = *m_, m = &mtmp;
    // set flag
    p->flag |= m? 0x1 : 0; // is paired in sequencing
    p->flag |= p->rid < 0? 0x4 : 0; // is mapped
    p->flag |= m && m->rid < 0? 0x8 : 0; // is mate mapped
    if (p->rid < 0 && m && m->rid >= 0) // copy mate to alignment
        p->rid = m->rid, p->pos = m->pos, p->is_rev = m->is_rev, p->n_cigar = 0;
    if (m && m->rid < 0 && p->rid >= 0) // copy alignment to mate
        m->rid = p->rid, m->pos = p->pos, m->is_rev = p->is_rev, m->n_cigar = 0;
    p->flag |= p->is_rev? 0x10 : 0; // is on the reverse strand
    p->flag |= m && m->is_rev? 0x20 : 0; // is mate on the reverse strand

    /* ─── Pre-size kstring once for the worst-case record length, then use
     * unsafe (no-realloc, no-bound-check) writes throughout. The estimate
     * covers: QNAME + FLAG + RNAME + POS + MAPQ + CIGAR + mate fields +
     * SEQ + QUAL + tags + headers; rare overflow falls back to ks_resize.
     */
    l_name = (int)strlen(s->name);
    {
        size_t md_len = 0;
        if (p->n_cigar)
            md_len = strlen((const char*)(p->cigar + p->n_cigar));
        size_t comm_len = s->comment ? strlen(s->comment) : 0;
        size_t anno_len = ((opt->flag & MEM_F_REF_HDR) && p->rid >= 0
                           && bns->anns[p->rid].anno && bns->anns[p->rid].anno[0])
                          ? strlen(bns->anns[p->rid].anno) : 0;
        size_t xa_len   = p->XA ? strlen(p->XA) : 0;
        size_t rg_len   = bwa_rg_id[0] ? strlen(bwa_rg_id) : 0;
        /* RNAME/RNEXT contig names can be very long (T2T-style assemblies).
         * Mate CIGAR (MC:Z) under V17 is also written via add_cigar(). */
        size_t rname_len = (p->rid >= 0) ? strlen(bns->anns[p->rid].name) : 1;
        size_t rnext_len = (m && m->rid >= 0 && m->rid != p->rid)
                           ? strlen(bns->anns[m->rid].name) : 1;
        size_t mate_cigar_len = (m && m->n_cigar) ? (size_t)m->n_cigar * 16 : 0;
        /* SA:Z worst case (only emitted when n>1 and there are other primary
         * hits, but priced unconditionally so the local ks_resize at the
         * SA:Z site can be dropped — single resize keeps the unsafe-fast
         * write invariant for the post-SA:Z tags (XA, HN, comment, XR:Z). */
        size_t sa_need = 16;
        for (int j = 0; j < n; ++j) {
            if (j == which || (list[j].flag&0x100)) continue;
            sa_need += strlen(bns->anns[list[j].rid].name) + 32 + list[j].n_cigar * 16;
        }
        size_t need = (size_t)l_name
                    + (size_t)s->l_seq * 2     /* SEQ + QUAL */
                    + (size_t)p->n_cigar * 16  /* primary CIGAR */
                    + mate_cigar_len           /* MC:Z mate CIGAR */
                    + rname_len + rnext_len    /* RNAME + RNEXT */
                    + md_len + xa_len + rg_len + comm_len + anno_len
                    + sa_need                  /* SA:Z (multi-supp) */
                    + 256;                      /* slack for ints + tags */
        ks_resize(str, str->l + need);
    }

    // print up to CIGAR
    kputsn_u(s->name, l_name, str); kputc_u('\t', str); // QNAME
    kputw_u((p->flag&0xffff) | (p->flag&0x10000? 0x100 : 0), str); kputc_u('\t', str); // FLAG
    if (p->rid >= 0) { // with coordinate
        kputs_u(bns->anns[p->rid].name, str); kputc_u('\t', str); // RNAME
        kputl_u(p->pos + 1, str); kputc_u('\t', str); // POS
        kputw_u(p->mapq, str); kputc_u('\t', str); // MAPQ
#if OLD // from v15
        if (p->n_cigar) { // aligned
            for (i = 0; i < p->n_cigar; ++i) {
                int c = p->cigar[i]&0xf;
                if (!(opt->flag&MEM_F_SOFTCLIP) && !p->is_alt && (c == 3 || c == 4))
                    c = which? 4 : 3; // use hard clipping for supplementary alignments
                kputw_u(p->cigar[i]>>4, str); kputc_u("MIDSH"[c], str);
            }
        } else kputc_u('*', str);
#else
        // modification from v17
        add_cigar(opt, p, str, which);
#endif
    } else kputsn_u("*\t0\t0\t*", 7, str); // without coordinte
    kputc_u('\t', str);

    // print the mate position if applicable
    if (m && m->rid >= 0) {
        if (p->rid == m->rid) kputc_u('=', str);
        else kputs_u(bns->anns[m->rid].name, str);
        kputc_u('\t', str);
        kputl_u(m->pos + 1, str); kputc_u('\t', str);
        if (p->rid == m->rid) {
            int64_t p0 = p->pos + (p->is_rev? get_rlen(p->n_cigar, p->cigar) - 1 : 0);
            int64_t p1 = m->pos + (m->is_rev? get_rlen(m->n_cigar, m->cigar) - 1 : 0);
            if (m->n_cigar == 0 || p->n_cigar == 0) kputc_u('0', str);
            else kputl_u(-(p0 - p1 + (p0 > p1? 1 : p0 < p1? -1 : 0)), str);
        } else kputc_u('0', str);
    } else kputsn_u("*\t0\t0", 5, str);
    kputc_u('\t', str);

    // print SEQ and QUAL — SIMD via NEON tbl1q_u8 byte LUT.
    if (p->flag & 0x100) { // for secondary alignments, don't write SEQ and QUAL
        kputsn_u("*\t*", 3, str);
    } else if (!p->is_rev) { // forward strand
        int qb = 0, qe = s->l_seq;
        if (p->n_cigar && which && !(opt->flag&MEM_F_SOFTCLIP) && !p->is_alt) {
            if ((p->cigar[0]&0xf) == 4 || (p->cigar[0]&0xf) == 3) qb += p->cigar[0]>>4;
            if ((p->cigar[p->n_cigar-1]&0xf) == 4 || (p->cigar[p->n_cigar-1]&0xf) == 3) qe -= p->cigar[p->n_cigar-1]>>4;
        }
        const int n_seq = qe - qb;
        sam_encode_seq_fwd(str->s + str->l, (const uint8_t*)(s->seq + qb), n_seq);
        str->l += n_seq;
        kputc_u('\t', str);
        if (s->qual) {
            memcpy(str->s + str->l, s->qual + qb, (size_t)n_seq);
            str->l += n_seq;
        } else kputc_u('*', str);
    } else { // reverse strand
        int qb = 0, qe = s->l_seq;
        if (p->n_cigar && which && !(opt->flag&MEM_F_SOFTCLIP) && !p->is_alt) {
            if ((p->cigar[0]&0xf) == 4 || (p->cigar[0]&0xf) == 3) qe -= p->cigar[0]>>4;
            if ((p->cigar[p->n_cigar-1]&0xf) == 4 || (p->cigar[p->n_cigar-1]&0xf) == 3) qb += p->cigar[p->n_cigar-1]>>4;
        }
        const int n_seq = qe - qb;
        sam_encode_seq_rev(str->s + str->l, (const uint8_t*)(s->seq + qb), n_seq);
        str->l += n_seq;
        kputc_u('\t', str);
        if (s->qual) {
            sam_encode_qual_rev(str->s + str->l, (const char*)(s->qual + qb), n_seq);
            str->l += n_seq;
        } else kputc_u('*', str);
    }

    // print optional tags
    if (p->n_cigar) {
        kputsn_u("\tNM:i:", 6, str); kputw_u(p->NM, str);
        kputsn_u("\tMD:Z:", 6, str); kputs_u((char*)(p->cigar + p->n_cigar), str);
    }
#if V17
    if (m && m->n_cigar) { kputsn_u("\tMC:Z:", 6, str); add_cigar(opt, m, str, which); }
#endif
    if (m) { kputsn_u("\tMQ:i:", 6, str); kputw_u(m->mapq, str); }
    if (p->score >= 0) { kputsn_u("\tAS:i:", 6, str); kputw_u(p->score, str); }
    if (p->sub >= 0) { kputsn_u("\tXS:i:", 6, str); kputw_u(p->sub, str); }
    if (bwa_rg_id[0]) { kputsn_u("\tRG:Z:", 6, str); kputs_u(bwa_rg_id, str); }
    if (!(p->flag & 0x100)) { // not multi-hit
        for (i = 0; i < n; ++i)
            if (i != which && !(list[i].flag&0x100)) break;
        if (i < n) { // there are other primary hits; output them
            /* SA:Z capacity is reserved in the outer pre-size above. */
            kputsn_u("\tSA:Z:", 6, str);
            for (i = 0; i < n; ++i) {
                const mem_aln_t *r = &list[i];
                int k;
                if (i == which || (r->flag&0x100)) continue;
                kputs_u(bns->anns[r->rid].name, str); kputc_u(',', str);
                kputl_u(r->pos+1, str); kputc_u(',', str);
                kputc_u("+-"[r->is_rev], str); kputc_u(',', str);
                for (k = 0; k < r->n_cigar; ++k) {
                    kputw_u(r->cigar[k]>>4, str); kputc_u("MIDSH"[r->cigar[k]&0xf], str);
                }
                kputc_u(',', str); kputw_u(r->mapq, str);
                kputc_u(',', str); kputw_u(r->NM, str);
                kputc_u(';', str);
            }
        }
        if (p->alt_sc > 0)
            ksprintf(str, "\tpa:f:%.3f", (double)p->score / p->alt_sc);
    }

    if (p->XA) { kputsn_u("\tXA:Z:", 6, str); kputs_u(p->XA, str); }
    if (p->HN >= 0) { kputsn_u("\tHN:i:", 6, str); kputw_u(p->HN, str); }

    if (s->comment) { kputc_u('\t', str); kputs_u(s->comment, str); }
    if ((opt->flag&MEM_F_REF_HDR) && !opt->meth_mode && p->rid >= 0 && bns->anns[p->rid].anno != 0 && bns->anns[p->rid].anno[0] != 0) {
        /* Suppressed under --meth: that mode repurposes XR:Z for the Bismark
         * read-conversion direction (see docs/src/methylation/tags.md). */
        int tmp;
        kputsn_u("\tXR:Z:", 6, str);
        tmp = (int)str->l;
        kputs_u(bns->anns[p->rid].anno, str);
        for (i = tmp; i < (int)str->l; ++i) // replace TAB in the comment to SPACE
            if (str->s[i] == '\t') str->s[i] = ' ';
    }
    kputc_u('\n', str);
    str->s[str->l] = 0;  /* single trailing NUL — unsafe writers skipped this */
}

mem_aln_t mem_reg2aln(const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac, int l_query, const char *query_, const mem_alnreg_t *ar, const char *meth_orig_query)
{
    mem_aln_t a;
    int i, w2, tmp, qb, qe, NM, score, is_rev, last_sc = -(1<<30), l_MD;
    int64_t pos, rb, re;
    uint8_t *query;

    memset(&a, 0, sizeof(mem_aln_t));
    a.HN = -1; // sentinel: HN not computed unless caller fills from mem_gen_alt out_hn
    a.meth_hypothesis = -1; // non-meth default; overwritten below from ar
    if (ar == 0 || ar->rb < 0 || ar->re < 0) { // generate an unmapped record
        a.rid = -1; a.pos = -1; a.flag |= 0x4;
        return a;
    }
    a.meth_hypothesis = ar->meth_hypothesis; // carry hypothesis to the output layer
    qb = ar->qb, qe = ar->qe;
    rb = ar->rb, re = ar->re;
    /* D3 (--meth, PR-5): output CIGAR/NM/MD must reflect the ORIGINAL read vs the
     * ORIGINAL reference in the original alphabet, so a real (non-converting) SNP
     * shows as a mismatch and the bisulfite C→T/G→A conversion shows as read-T at
     * ref-C (the substrate a downstream double-masking step like Revelio expects;
     * plan §4 / spec §6.1). When `meth_orig_query` is supplied, regen from the
     * original bases; else fall back to the projected `query_` (non-meth path,
     * byte-for-byte identical to legacy). */
    const int use_meth_orig = (opt->meth_mode && meth_orig_query != NULL
                               && ar->meth_hypothesis >= 0);
    const char *regen_query = use_meth_orig ? meth_orig_query : query_;
    /* Reuse a per-thread nt4 scratch instead of a malloc/free per region. The
     * buffer is consumed within this call (copied into a.cigar via the DP), so
     * one live instance per thread is safe; grown to the high-water l_query. */
    static thread_local u8vec_scratch_t t_query;
    if (t_query.v.m < (size_t)l_query) kv_resize(uint8_t, t_query.v, (size_t)l_query);
    query = t_query.v.a;
    for (i = 0; i < l_query; ++i) // convert to the nt4 encoding
        query[i] = regen_query[i] < 5? regen_query[i] : nst_nt4_table[(int)regen_query[i]];
    a.mapq = ar->secondary < 0? mem_approx_mapq_se(opt, ar) : 0;
    if (ar->secondary >= 0) a.flag |= 0x100; // secondary alignment
    tmp = infer_bw(qe - qb, re - rb, ar->truesc, opt->a, opt->o_del, opt->e_del);
    w2  = infer_bw(qe - qb, re - rb, ar->truesc, opt->a, opt->o_ins, opt->e_ins);
    w2 = w2 > tmp? w2 : tmp;
    if (bwa_verbose >= 4) fprintf(stderr, "* Band width: inferred=%d, cmd_opt=%d, alnreg=%d\n", w2, opt->w, ar->w);
    if (w2 > opt->w) w2 = w2 < ar->w? w2 : ar->w;
    i = 0; a.cigar = 0;
    /* D3 (--meth, PR-5): native-alphabet output regen.
     *   - `regen_query` is the ORIGINAL read (meth_orig_query) under --meth, else
     *     the projected read for non-meth (legacy, unchanged).
     *   - `regen_mat` is the per-hypothesis ASYMMETRIC matrix under --meth so the
     *     CIGAR/gap shape matches the asymmetric extension that selected this
     *     region (using the symmetric matrix with the ORIGINAL read would
     *     re-penalize every conversion as a mismatch and can shift gaps); for
     *     non-meth it is the symmetric opt->mat (legacy, unchanged).
     *
     * NOTE (--meth NM/MD-vs-conversion policy, RESOLVED — plan §4): bwa_gen_cigar2
     * computes NM/MD by LITERAL base comparison (query[x+i] != rseq[y+i]),
     * independent of `regen_mat`. So a bisulfite C→T (OT) / G→A (OB) at a
     * ref-C/ref-G is COUNTED as one mismatch in NM and emitted in MD, exactly like
     * a real SNP — it is NOT hidden. This is deliberate: it is the substrate a
     * downstream Revelio-style double-mask + conventional caller expects (the
     * conversion is visible as read-T at ref-C, then BQ-masked downstream; the
     * aligner does no masking). The conversion is NOT double-counted against the
     * alignment SCORE: the asymmetric `regen_mat` frees the conversion cell to a
     * match (+a) for the DP, while NM/MD count it once as a literal mismatch — two
     * separate concerns. Hiding conversions from NM/MD (e.g. a Bismark-style
     * XM-only consumer) would be a different, explicit choice and must be made
     * here, not assumed. */
    /* D3 (--meth, fix): the CIGAR regen must use the SAME strand-adjusted matrix
     * as the extension (see mem_alnreg_t.meth_strand_hyp). ar->rb is final here,
     * so derive is_rev directly (rb in doubled-pac; >= l_pac = reverse) and flip
     * the hypothesis; otherwise reverse-strand reads regen with the wrong matrix
     * and their conversions become CIGAR mismatches / soft-clips. */
    int regen_hyp = ar->meth_hypothesis;
    if (regen_hyp >= 0 && ar->rb >= bns->l_pac) regen_hyp ^= 1;
    const int8_t *regen_mat = use_meth_orig
        ? mem_opt_meth_mat(opt, regen_hyp) : opt->mat;
    do {
        free(a.cigar);
        w2 = w2 < opt->w<<2? w2 : opt->w<<2;
        a.cigar = bwa_gen_cigar2(regen_mat, opt->o_del, opt->e_del, opt->o_ins, opt->e_ins, w2, bns->l_pac, pac, qe - qb, (uint8_t*)&query[qb], rb, re, &score, &a.n_cigar, &NM);
        if (bwa_verbose >= 4) fprintf(stderr, "* Final alignment: w2=%d, global_sc=%d, local_sc=%d\n", w2, score, ar->truesc);
        if (score == last_sc || w2 == opt->w<<2) break; // it is possible that global alignment and local alignment give different scores
        last_sc = score;
        w2 <<= 1;
    } while (++i < 3 && score < ar->truesc - opt->a);
    assert(a.cigar != NULL);
    l_MD = strlen((char*)(a.cigar + a.n_cigar)) + 1;
    a.NM = NM;
    pos = bns_depos(bns, rb < bns->l_pac? rb : re - 1, &is_rev);
    a.is_rev = is_rev;
    if (a.n_cigar > 0) { // squeeze out leading or trailing deletions
        assert(a.cigar != NULL);
        if ((a.cigar[0]&0xf) == 2) {
            pos += a.cigar[0]>>4;
            --a.n_cigar;
            memmove(a.cigar, a.cigar + 1, a.n_cigar * 4 + l_MD);
        } else if ((a.cigar[a.n_cigar-1]&0xf) == 2) {
            --a.n_cigar;
            memmove(a.cigar + a.n_cigar, a.cigar + a.n_cigar + 1, l_MD); // MD needs to be moved accordingly
        }
    }
    if (qb != 0 || qe != l_query) { // add clipping to CIGAR
        int clip5, clip3;
        clip5 = is_rev? l_query - qe : qb;
        clip3 = is_rev? qb : l_query - qe;
        a.cigar = (uint32_t*) realloc(a.cigar, 4 * (a.n_cigar + 2) + l_MD);
        if (clip5) {
            memmove(a.cigar+1, a.cigar, a.n_cigar * 4 + l_MD); // make room for 5'-end clipping
            a.cigar[0] = clip5<<4 | 3;
            ++a.n_cigar;
        }
        if (clip3) {
            memmove(a.cigar + a.n_cigar + 1, a.cigar + a.n_cigar, l_MD); // make room for 3'-end clipping
            a.cigar[a.n_cigar++] = clip3<<4 | 3;
        }
    }
    a.rid = bns_pos2rid(bns, pos);
    assert(a.rid == ar->rid);
    a.pos = pos - bns->anns[a.rid].offset;
    a.score = ar->score; a.sub = ar->sub > ar->csub? ar->sub : ar->csub;
    a.is_alt = ar->is_alt; a.alt_sc = ar->alt_sc;
    /* query is the reused thread-local nt4 scratch — no free here. */
    return a;
}

/*****************************
 * Basic hit->SAM conversion *
 *****************************/

static inline int infer_bw(int l1, int l2, int score, int a, int q, int r)
{
    int w;
    if (l1 == l2 && l1 * a - score < (q + r - a)<<1) return 0; // to get equal alignment length, we need at least two gaps
    w = ((double)((l1 < l2? l1 : l2) * a - score - q) / r + 2.);
    if (w < abs(l1 - l2)) w = abs(l1 - l2);
    return w;
}

static inline int get_rlen(int n_cigar, const uint32_t *cigar)
{
    int k, l;
    for (k = l = 0; k < n_cigar; ++k) {
        int op = cigar[k]&0xf;
        if (op == 0 || op == 2)
            l += cigar[k]>>4;
    }
    return l;
}

/************************ New functions, version2*****************************************/
#define _get_pac(pac, l) ((pac)[(l)>>2]>>((~(l)&3)<<1)&3)

void* _mm_realloc(void *ptr, int64_t csize, int64_t nsize, int16_t dsize) {
    if (nsize <= csize)
    {
        fprintf(stderr, "Shringking not supported yet.\n");
        return ptr;
    }
    void *nptr = _mm_malloc(nsize * dsize, 64);
    assert(nptr != NULL);
    /* csize is in elements (matching nsize), not bytes — multiply by dsize so
     * callers passing dsize > 1 (sa_coord with sizeof(int64_t), matchArray
     * with sizeof(SMEM)) get the full pre-grow contents copied across, not
     * just the first csize bytes. Guard the NULL/empty lazy-init path
     * (matchArray[tid] starts NULL with wsize_mem[tid]=0). */
    if (ptr != NULL && csize > 0) {
        memcpy(nptr, ptr, (size_t)csize * (size_t)dsize);
    }
    _mm_free(ptr);

    return nptr;
}

inline void sortPairsLenExt(SeqPair *pairArray, int32_t count, SeqPair *tempArray,
                            int32_t *hist, int &numPairs128, int &numPairs16,
                            int &numPairs1, int score_a, int w, int zdrop)
{
    int32_t i;
    numPairs128 = numPairs16 = numPairs1 = 0;
    if (count <= 0) return;   /* empty extension side: nothing to sort/count */

    int32_t *hist2 = hist + MAX_SEQ_LEN8;
    int32_t *hist3 = hist + MAX_SEQ_LEN8 + MAX_SEQ_LEN16;

    /* The counting-sort bin space is three concatenated regions: 8-bit
     * (hist[0..MAX_SEQ_LEN8)), 16-bit (hist2[0..MAX_SEQ_LEN16)) and scalar
     * (hist3[0]). The 16-bit region alone is 32768 bins, so clearing +
     * cumulating the whole range every call was O(32768) irrespective of the
     * pair count. Window the 16-bit region to the actual minval range: a 16-bit
     * pair's bin is exactly its minval, so every touched bin lies within
     * [min minval, max minval]; bins below the window hold no pairs (they add 0
     * to the cross-region cumulative sum, so starting the sum at that offset is
     * exact) and bins above are never read. The sorted output and the numPairs*
     * tier counts are therefore byte-for-byte identical to the full-range sort.
     * The 8-bit region (1088 bins) is small, so it is still cleared in full to
     * avoid reasoning about its clamped bins. */
    int64_t mn = 0, mx = 0;
    for (i = 0; i < count; i++) {
        SeqPair sp = pairArray[i];
        int64_t minval = (int64_t)sp.h0 + (int64_t)min_(sp.len1, sp.len2) * (int64_t)score_a;
        if (i == 0) { mn = mx = minval; }
        else { if (minval < mn) mn = minval; if (minval > mx) mx = minval; }
    }
    int32_t w_lo = mn < 0 ? 0 : (mn >= MAX_SEQ_LEN16 ? MAX_SEQ_LEN16 - 1 : (int32_t)mn);
    int32_t w_hi = mx < 0 ? 0 : (mx >= MAX_SEQ_LEN16 ? MAX_SEQ_LEN16 - 1 : (int32_t)mx);

    for (i = 0; i < MAX_SEQ_LEN8; i++) hist[i] = 0;
    for (i = w_lo; i <= w_hi; i++)     hist2[i] = 0;
    hist3[0] = 0;

    /* minval is the counting-sort bin (sorts the 8-bit bucket by score for
     * lane packing). The tier decision is the safe envelope, not minval:
     * re-baselining handles arbitrary score, so a pair can be 8-bit even with
     * minval >= MAX_SEQ_LEN8. Clamp the bin to keep it in range — only the
     * ordering, not correctness, depends on the exact key. */
    for(i = 0; i < count; i++)
    {
        SeqPair sp = pairArray[i];
        // int minval = sp.h0 + max_(sp.len1, sp.len2);
        int64_t minval = (int64_t)sp.h0 + (int64_t)min_(sp.len1, sp.len2) * (int64_t)score_a;
        if (bsw8_envelope_ok(sp.len1, sp.len2, w, score_a, zdrop, sp.h0)) {
            int bin = minval < 0 ? 0 : (minval < MAX_SEQ_LEN8 ? (int)minval : MAX_SEQ_LEN8 - 1);
            hist[bin]++;
        }
        else if(sp.len1 < MAX_SEQ_LEN16 && sp.len2 < MAX_SEQ_LEN16 && minval >= 0 && minval < MAX_SEQ_LEN16)
            hist2[(int)minval] ++;
        else
            hist3[0] ++;
    }

    int32_t cumulSum = 0;
    for(i = 0; i < MAX_SEQ_LEN8; i++)
    {
        int32_t cur = hist[i];
        hist[i] = cumulSum;
        cumulSum += cur;
    }
    /* Skipping hist2[0..w_lo) is exact: those bins were cleared to (and hold) no
     * counts, so they would contribute 0 to cumulSum anyway. */
    for(i = w_lo; i <= w_hi; i++)
    {
        int32_t cur = hist2[i];
        hist2[i] = cumulSum;
        cumulSum += cur;
    }
    hist3[0] = cumulSum;

    for(i = 0; i < count; i++)
    {
        SeqPair sp = pairArray[i];
        // int minval = sp.h0 + max_(sp.len1, sp.len2);
        int64_t minval = (int64_t)sp.h0 + (int64_t)min_(sp.len1, sp.len2) * (int64_t)score_a;

        if (bsw8_envelope_ok(sp.len1, sp.len2, w, score_a, zdrop, sp.h0))
        {
            int bin = minval < 0 ? 0 : (minval < MAX_SEQ_LEN8 ? (int)minval : MAX_SEQ_LEN8 - 1);
            int32_t pos = hist[bin];
            tempArray[pos] = sp;
            hist[bin]++;
            numPairs128 ++;
        }
        else if (sp.len1 < MAX_SEQ_LEN16 && sp.len2 < MAX_SEQ_LEN16 && minval >= 0 && minval < MAX_SEQ_LEN16) {
            int32_t pos = hist2[(int)minval];
            tempArray[pos] = sp;
            hist2[(int)minval]++;
            numPairs16 ++;
        }
        else {
            int32_t pos = hist3[0];
            tempArray[pos] = sp;
            hist3[0]++;
            numPairs1 ++;
        }
    }

    for(i = 0; i < count; i++) {
        pairArray[i] = tempArray[i];
    }
}

inline void sortPairsLen(SeqPair *pairArray, int32_t count, SeqPair *tempArray, int32_t *hist)
{

    int32_t i;
    if (count <= 0) return;

    /* Bins are keyed directly by len1, so only [minLen, maxLen] is ever touched;
     * clearing + cumulating the full [0, MAX_SEQ_LEN16] range every call was
     * O(32768) regardless of count. Window to the actual length range — output
     * is identical because every pair falls in [minLen, maxLen], bins below add
     * 0 to the cumulative sum, and bins above are never read. */
    int32_t minLen = pairArray[0].len1, maxLen = pairArray[0].len1;
    for(i = 1; i < count; i++)
    {
        int32_t L = pairArray[i].len1;
        if (L < minLen) minLen = L;
        if (L > maxLen) maxLen = L;
    }
    for(i = minLen; i <= maxLen; i++) hist[i] = 0;

    for(i = 0; i < count; i++)
    {
        SeqPair sp = pairArray[i];
        hist[sp.len1]++;
    }
    int32_t cumulSum = 0;
    for(i = minLen; i <= maxLen; i++)
    {
        int32_t cur = hist[i];
        hist[i] = cumulSum;
        cumulSum += cur;
    }

    for(i = 0; i < count; i++)
    {
        SeqPair sp = pairArray[i];
        int32_t pos = hist[sp.len1];

        tempArray[pos] = sp;
        hist[sp.len1]++;
    }

    for(i = 0; i < count; i++) {
        pairArray[i] = tempArray[i];
    }
}

/* ------------------------------------------------------------------------- *
 * D3 (--meth, A1/PR-4): per-hypothesis banded-SW tier dispatch.
 *
 * A worker batch interleaves R1 (OT, C->T) and R2 (OB, G->A) reads, so a single
 * per-object substitution matrix cannot score every pair in the batch correctly
 * (the PR-3 placeholder fell back to the SYMMETRIC matrix for mixed batches).
 * bsw_run_tier() partitions a tier's pair array in place by the owning alnreg's
 * meth_hypothesis (OT = odd, OB = even, none = <0) and runs the matching SW
 * object on each contiguous run, so each pair is scored in original space with
 * its asymmetric OT/OB matrix. Each run is single-hypothesis, so the kernel's
 * rank-1 fast path still applies. The kernel re-sorts pairs by length internally
 * and the retry loop reads results back via seqid/regid, so the in-place reorder
 * is safe.
 *
 * Outside --meth the partition is skipped and the call is byte-for-byte the
 * original single kernel invocation on the symmetric object.
 * ------------------------------------------------------------------------- */
enum BswMethTier { BSW_TIER_SCALAR, BSW_TIER_16, BSW_TIER_8 };

/* One tier kernel call over pa[0..n). `sort_scratch` must be a SeqPair buffer
 * DISTINCT from `pa` (sortPairsLen is a counting sort that scatters through it;
 * aliasing pa would corrupt the sort). Scalar tier ignores sort_scratch. */
static inline void bsw_tier_kernel(BswMethTier tier, IBandedPairWiseSW *bsw,
        SeqPair *pa, uint8_t *ref, uint8_t *qer, int n, int nthreads, int32_t w,
        SeqPair *sort_scratch, int32_t *hist)
{
    if (n <= 0) return;
    switch (tier) {
    case BSW_TIER_SCALAR:
        bsw->scalarBandedSWAWrapper(pa, ref, qer, n, nthreads, w);
        break;
    case BSW_TIER_16:
#if !HAVE_BSW_VECTOR_8_16
        bsw->scalarBandedSWAWrapper(pa, ref, qer, n, nthreads, w);
#else
        sortPairsLen(pa, n, sort_scratch, hist);
        bsw->getScores16(pa, ref, qer, n, nthreads, w);
#endif
        break;
    case BSW_TIER_8:
#if !HAVE_BSW_VECTOR_8_16
        bsw->scalarBandedSWAWrapper(pa, ref, qer, n, nthreads, w);
#else
        sortPairsLen(pa, n, sort_scratch, hist);
        /* int8 diagonal encoding can't represent a band >= 128; divert to the
         * 16-bit kernel once the doubling retry pushes w past the envelope. */
        if (w <= BSW8_MAX_W) bsw->getScores8(pa, ref, qer, n, nthreads, w);
        else                 bsw->getScores16(pa, ref, qer, n, nthreads, w);
#endif
        break;
    }
}

/* Dispatch one tier of `nump` pairs. Non-meth: a single call on `sym` with the
 * original `base_scratch` (= seqPairArrayAux). Meth: 3-way in-place partition by
 * hypothesis, then one call per non-empty run on the matching object, using
 * `meth_scratch` (= the retry loop's pair_ar_aux, always distinct from pair_ar)
 * as the sort scratch. */
static inline void bsw_run_tier(BswMethTier tier, const mem_opt_t *opt,
        const mem_alnreg_v *av_v,
        IBandedPairWiseSW *sym, IBandedPairWiseSW *ot, IBandedPairWiseSW *ob,
        SeqPair *pair_ar, uint8_t *ref, uint8_t *qer, int nump,
        int nthreads, int32_t w,
        SeqPair *base_scratch, SeqPair *meth_scratch, int32_t *hist)
{
    if (nump <= 0) return;
    if (!opt->meth_mode) {
        bsw_tier_kernel(tier, sym, pair_ar, ref, qer, nump, nthreads, w,
                        base_scratch, hist);
        return;
    }
    /* Dutch-flag partition into [OT (key 0) | OB (key 1) | SYM (key 2)]. Order
     * within/among runs is irrelevant: the kernel re-sorts by length and the
     * caller maps each result back to its alnreg by seqid/regid. */
    int lo = 0, mid = 0, hi = nump - 1;
    while (mid <= hi) {
        const SeqPair *sp = &pair_ar[mid];
        /* strand-adjusted hypothesis (see mem_alnreg_t.meth_strand_hyp): the
         * extension matrix must match the conversion direction AS SEEN against
         * the (possibly RC) reference window, not the raw genome-strand hypothesis. */
        int h = av_v[sp->seqid].a[sp->regid].meth_strand_hyp;
        int key = (h < 0) ? 2 : ((h & 1) ? 0 : 1);
        if (key == 0)      { SeqPair t = pair_ar[lo]; pair_ar[lo]  = pair_ar[mid]; pair_ar[mid] = t; lo++; mid++; }
        else if (key == 1) { mid++; }
        else               { SeqPair t = pair_ar[mid]; pair_ar[mid] = pair_ar[hi]; pair_ar[hi]  = t; hi--; }
    }
    int n_ot = lo, n_ob = mid - lo, n_sym = nump - mid;
    bsw_tier_kernel(tier, ot,  pair_ar,       ref, qer, n_ot,  nthreads, w, meth_scratch, hist);
    bsw_tier_kernel(tier, ob,  pair_ar + lo,  ref, qer, n_ob,  nthreads, w, meth_scratch, hist);
    bsw_tier_kernel(tier, sym, pair_ar + mid, ref, qer, n_sym, nthreads, w, meth_scratch, hist);
}

/* Restructured BSW parent function */
#define FAC 8
#define PFD 2

// ungapped diagonal-extension analyzer.
//
// Walks the N-step diagonal once (SIMD scan + bitmap + scalar trajectory).
// Emits three possible status codes:
//
//   FP_STATUS_HIT      — ungapped is provably optimal (≤ x_threshold
//                        mismatches, no ambig). Caller skips banded SW
//                        and uses sp.score/qle/gscore/gtle directly.
//
//   FP_STATUS_TIGHT    — fast-path fails (too many mismatches) but the
//                        ungapped score bounds the useful SW band:
//                          tight_band = ceil((min(len1,len2)·a - ungapped_score)
//                                            / (o_min + e_min))
//                        Caller runs SW with this tight band instead of
//                        opt->w, and skips MAX_BAND_TRY retries (the
//                        bound is an upper bound on any gapped score).
//
//   FP_STATUS_FALLBACK — ambig base or out-of-range length. Caller uses
//                        opt->w with the full retry loop.
//
// Derivation of tight_band: for any alignment with band offset B from
// diagonal, min B gaps are required; cost ≥ B · (o_min + e_min). Max
// alignment score (all matches) ≤ min(len1, len2) · a. For any gapped
// alignment to beat the observed ungapped max score S:
//   min_len·a - B·(o_min+e_min) > S
//   B < (min_len·a - S) / (o_min + e_min)
// So any band ≥ tight_band is sufficient; narrower bands suffice too
// when the gapped alternative score is < min_len·a. Starting SW at
// tight_band is strictly correct and avoids over-banded DP work.
#define FP_N_MAX 128
#define FP_STATUS_FALLBACK  0
#define FP_STATUS_HIT       1
#define FP_STATUS_TIGHT     2

/* Q2 helper: bin an extension's final alignment score into one of 8 buckets.
 * Bins: {0-10, 11-25, 26-50, 51-75, 76-100, 101-125, 126-150, 151+}. */
static inline int ugp_score_bin(int s)
{
    if      (s <=  10) return 0;
    else if (s <=  25) return 1;
    else if (s <=  50) return 2;
    else if (s <=  75) return 3;
    else if (s <= 100) return 4;
    else if (s <= 125) return 5;
    else if (s <= 150) return 6;
    else               return 7;
}

/* Q3 helper: bin a "delta from perfect extension" into one of 8 buckets.
 * delta = h0 + a * len2 - aln_score, where (h0 + a * len2) is the score of a
 * hypothetical all-match ungapped extension. delta == 0 ⇒ extension was
 * perfect (HIT path with all matches); delta grows with mismatches and gaps.
 * Bins: {0, 1-5, 6-10, 11-25, 26-50, 51-75, 76-100, 101+}. */
static inline int ugp_delta_bin(int d)
{
    if      (d ==   0) return 0;
    else if (d <=   5) return 1;
    else if (d <=  10) return 2;
    else if (d <=  25) return 3;
    else if (d <=  50) return 4;
    else if (d <=  75) return 5;
    else if (d <= 100) return 6;
    else               return 7;
}

/* Q1+Q2+Q3 outcome helpers. Each post-SW dispatch (LEFT scalar / 16-bit /
 * 8-bit; RIGHT scalar / 16-bit / 8-bit) commits a final outcome and bumps
 * the same set of profiling counters. Hoisted so categorization changes
 * (extra category, adjusted bin edges, fourth tier, ...) live in one place
 * instead of drifting across six near-identical inlined blocks. */
static inline void ugp_record_left_outcome(const SeqPair *sp, int a_match, int tid)
{
    int _u = (sp->qle == sp->tle);
    tprof[UGP_OUTCOME_BASE + 0 * 2 + (_u ? 0 : 1)][tid]++;
    tprof[UGP_SCORE_HIST_BASE + 0 * UGP_SCORE_HIST_NBINS
          + ugp_score_bin(sp->score)][tid]++;
    /* Q3: per-category delta-from-perfect histograms. Non-HIT by
     * construction at the post-SW commit. cat0=ALL; cat1 (ungapped final)
     * or cat2 (gapped final) by sp->qle == sp->tle proxy; cat4
     * (ungap_final ∩ ~HIT) when ungapped. */
    int _perfect = sp->h0 + a_match * sp->len2;
    int _bin_ung = ugp_delta_bin(_perfect - sp->ugp_walk_score);
    int _bin_fin = ugp_delta_bin(_perfect - sp->score);
    tprof[UGP_L_CAT_UNG_BASE + 0 * UGP_CAT_NBINS + _bin_ung][tid]++;
    tprof[UGP_L_CAT_FIN_BASE + 0 * UGP_CAT_NBINS + _bin_fin][tid]++;
    if (_u) {
        tprof[UGP_L_CAT_UNG_BASE + 1 * UGP_CAT_NBINS + _bin_ung][tid]++;
        tprof[UGP_L_CAT_FIN_BASE + 1 * UGP_CAT_NBINS + _bin_fin][tid]++;
        tprof[UGP_L_CAT_UNG_BASE + 4 * UGP_CAT_NBINS + _bin_ung][tid]++;
        tprof[UGP_L_CAT_FIN_BASE + 4 * UGP_CAT_NBINS + _bin_fin][tid]++;
    } else {
        tprof[UGP_L_CAT_UNG_BASE + 2 * UGP_CAT_NBINS + _bin_ung][tid]++;
        tprof[UGP_L_CAT_FIN_BASE + 2 * UGP_CAT_NBINS + _bin_fin][tid]++;
    }
}

static inline void ugp_record_right_outcome(const SeqPair *sp, int tid)
{
    int _u = (sp->qle == sp->tle);
    tprof[UGP_OUTCOME_BASE + 1 * 2 + (_u ? 0 : 1)][tid]++;
    tprof[UGP_SCORE_HIST_BASE + 1 * UGP_SCORE_HIST_NBINS
          + ugp_score_bin(sp->score)][tid]++;
}

/* Q3 helper: would-be ungapped extension score for arbitrary N. Mirrors the
 * HIT-path scalar walk in ungapped_analyze (cur with floor at 0; max_sc
 * tracker; ambig terminates the walk to match analyze's FALLBACK semantics).
 * Returned score is what an ungapped extension would produce on this pair
 * regardless of whether the fast-path actually triggered HIT. */
static inline int ungapped_walk_score(const uint8_t *qs, const uint8_t *rs,
                                       int N, int h0, int a, int b)
{
    int cur = h0, max_sc = h0;
    for (int j = 0; j < N; j++) {
        uint8_t qj = qs[j], rj = rs[j];
        if (qj >= 4 || rj >= 4) break;       /* ambig: terminate the walk */
        if (cur == 0) continue;
        if (qj == rj) cur += a;
        else {
            cur -= b;
            if (cur < 0) cur = 0;
        }
        if (cur >= max_sc) max_sc = cur;
    }
    return max_sc;
}

/* optimal ungapped score (no floor) from a per-base mismatch bitmap.
 *
 *     score(i) = h0 + (i − mismatches_i)·a − mismatches_i·b
 *              = h0 + i·a − mismatches_i·(a+b)        for i ∈ [0, N]
 *
 * f(i) is monotone non-decreasing within match runs and drops by b at each
 * mismatch. Local maxima therefore live at one of:
 *   (a) i = 0                                    f(0)  = h0
 *   (b) i = p,  bitmap[p] == 1                   f(p)  = h0 + p·a − k·(a+b)
 *                                                where k = mismatches in [0,p)
 *   (c) i = N                                    f(N)  = h0 + N·a − total·(a+b)
 *
 * Iterating set bits of the bitmap via ctz + clear-low is O(total_mis), not
 * O(N) — typically a handful of operations for HIT/TIGHT candidates.
 *
 * Unlike ungapped_walk_score, there is NO floor: score is allowed to dip
 * below 0 and recover. This yields max_sc ≥ walk's max_sc, giving a strictly
 * tighter (still safe) bound on the SW band in the tight_band proof. The
 * walk's floor exists to mirror SW kernel local-SW semantics for SAM
 * byte-identity on qle/gscore — irrelevant for the band proof. */
static inline int ungapped_max_sc_from_bitmap(uint64_t mis_lo, uint64_t mis_hi,
                                               int N, int h0, int a, int b)
{
    int max_sc = h0;            /* candidate (a): empty extension */
    int gap    = a + b;         /* score drop per mismatch (relative to match) */
    int k      = 0;             /* mismatches counted so far */
    uint64_t lo = mis_lo, hi = mis_hi;
    while (lo) {
        int p = __builtin_ctzll(lo);
        int s = h0 + p * a - k * gap;
        if (s > max_sc) max_sc = s;
        lo &= lo - 1;
        k++;
    }
    while (hi) {
        int p = 64 + __builtin_ctzll(hi);
        int s = h0 + p * a - k * gap;
        if (s > max_sc) max_sc = s;
        hi &= hi - 1;
        k++;
    }
    int s_end = h0 + N * a - k * gap;       /* candidate (c) */
    if (s_end > max_sc) max_sc = s_end;
    return max_sc;
}

static inline int ungapped_analyze(const uint8_t *qs, const uint8_t *rs, int N,
                                    int h0, int a, int b,
                                    int o_min, int e_min,
                                    int x_threshold, int default_w,
                                    int *out_score, int *out_qle,
                                    int *out_gscore, int *out_gtle,
                                    int *out_tight_band)
{
    if (N <= 0 || N > FP_N_MAX || x_threshold < 0) return FP_STATUS_FALLBACK;

    const __m128i v3 = _mm_set1_epi8(3);
    uint64_t mis_lo = 0, mis_hi = 0;
    int total_mis = 0;
    int i = 0;
    for (; i + 16 <= N; i += 16) {
        __m128i qv = _mm_loadu_si128((const __m128i *)(qs + i));
        __m128i rv = _mm_loadu_si128((const __m128i *)(rs + i));
        // AMBIG (code >= 4) detected via max > 3. ACGT = 0..3.
        __m128i mxv = _mm_max_epu8(qv, rv);
        if ((unsigned)_mm_movemask_epi8(_mm_cmpgt_epi8(mxv, v3)))
            return FP_STATUS_FALLBACK;  // ambig: give up on both paths
        __m128i eqv = _mm_cmpeq_epi8(qv, rv);
        unsigned mism_mask = (~(unsigned)_mm_movemask_epi8(eqv)) & 0xFFFFu;
        total_mis += __builtin_popcount(mism_mask);
        if (i < 64) {
            mis_lo |= ((uint64_t)mism_mask) << i;
            if (i + 16 > 64) mis_hi |= ((uint64_t)mism_mask) >> (64 - i);
        } else {
            mis_hi |= ((uint64_t)mism_mask) << (i - 64);
        }
    }
    for (; i < N; i++) {
        uint8_t qi = qs[i], ri = rs[i];
        if (qi >= 4 || ri >= 4) return FP_STATUS_FALLBACK;
        if (qi != ri) {
            total_mis++;
            if (i < 64) mis_lo |= (1ULL << i);
            else        mis_hi |= (1ULL << (i - 64));
        }
    }

    // precise max_sc (no floor) from the mismatch bitmap.
    //
    // An earlier formulation used `S = h0 + first_mis_pos·a` — a loose lower bound on the
    // walk's max_sc — to skip the scalar walk on TIGHT pairs. That preserved
    // walk-time but left the band proof loose: any matching run after the
    // first mismatch could not contribute to S, so the proven band stayed
    // wider than necessary.
    //
    // The bitmap-derived max_sc is the optimal ungapped score over all
    // prefix lengths in [0, N], with no floor (score is allowed to dip and
    // recover). It is ≥ the walk's max_sc, so substituting it into the band
    // proof yields a strictly tighter (still safe) bound. Computing it from
    // the bitmap is O(total_mis), independent of N.
    //
    // Band derivation (see comment above):
    //     B < (min_len·a − S − o_min) / e_min        with S = max_sc − h0
    // For LEFT extensions where the caller gates on len1 ≥ len2,
    // min_len = N. Substituting:
    //     numerator = N·a − (max_sc_proof − h0) − o_min
    //               = N·a + h0 − max_sc_proof − o_min
    int max_sc_proof = ungapped_max_sc_from_bitmap(mis_lo, mis_hi, N, h0, a, b);

    if (total_mis > x_threshold) {
        int64_t numerator = (int64_t)N * a + h0 - max_sc_proof - o_min;
        int band;
        if (numerator <= 0) {
            // Ungapped-optimal (same "tight_band = 0 sentinel" semantic as
            // walk-derived numerator ≤ 0: let SW run with fallback width).
            band = 0;
        } else if (e_min <= 0) {
            band = default_w;
        } else {
            band = (int)((numerator + e_min - 1) / e_min);
            if (band > default_w) band = default_w;
        }
        *out_tight_band = band;
        // Outputs unused on TIGHT; set sane values for any future caller.
        *out_score  = max_sc_proof;
        *out_qle    = 0;
        *out_gscore = 0;
        *out_gtle   = N;
        return FP_STATUS_TIGHT;
    }

    // Closed-form fast path: total_mis == 0 (perfect-match HIT).
    //
    // With no mismatches the scalar walk below is deterministic:
    //   cur starts at h0 and only increases (+a per iter), so it never
    //   touches the `cur == 0` early-skip after the first iteration.
    //   max_sc rises to h0 + N*a; max_i ends at N (the >= tie-break
    //   updates max_i on every step).
    //
    // The h0 > 0 guard preserves bit-identicality with the loop: when
    // h0 == 0 the loop's `if (cur == 0) continue` inhibits all updates,
    // leaving max_sc = 0 / max_i = 0; the closed form would instead
    // compute max_sc = N*a, breaking parity. h0 == 0 is pathological in
    // practice (caller passes a seed score) but the guard is a single
    // compare and free.
    if (total_mis == 0 && h0 > 0) {
        int score = h0 + N * a;
        *out_score      = score;
        *out_qle        = N;
        *out_gscore     = score;
        *out_gtle       = N;
        *out_tight_band = 0;   // unused on HIT (SW skipped)
        return FP_STATUS_HIT;
    }

    // HIT candidate: run the scalar walk for precise qle / gscore / max_sc.
    //
    // MAIN_CODE* local-SW semantics: once cur==0 in the ungapped path it
    // stays 0 (no e/f restart).
    //
    // Tie-break: SW's maxRS tracker (bandedSWA.cpp MAIN_CODE) updates the
    // position on BOTH strictly-greater and tied equal-to-current
    // comparisons — equivalent to "pick the rightmost position where the
    // max was achieved". We must mirror that (use >=) or qle/tle diverge
    // from SW on tied-score walks, breaking byte-identical SAM.
    int cur = h0, max_sc = h0, max_i = 0;
    for (int j = 0; j < N; j++) {
        int is_mis = (j < 64) ? (int)((mis_lo >> j) & 1ULL)
                              : (int)((mis_hi >> (j - 64)) & 1ULL);
        if (cur == 0) continue;
        if (!is_mis) cur += a;
        else {
            cur -= b;
            if (cur < 0) cur = 0;
        }
        if (cur >= max_sc) { max_sc = cur; max_i = j + 1; }
    }

    *out_score      = max_sc;
    *out_qle        = max_i;
    *out_gscore     = cur;
    *out_gtle       = N;
    *out_tight_band = 0;   // unused on HIT (SW skipped)
    return FP_STATUS_HIT;
}

/* ------------------------------------------------------------------------
 * Byte-identical contained-seed extension skip (--skip-contained-ext).
 *
 * A seed s strictly contained (same diagonal, query subinterval) in a longer
 * seed of the same chain is dominated by that longer seed's extension: the
 * longest seed on a diagonal is always extended, its alnreg covers s's region,
 * so s's own alnreg is purged post-extension anyway (the PE18 containment purge
 * later in this function). We detect this before extension and skip building s's
 * banded-SW pair, marking its alnreg purged (qb=qe=-1) exactly as PE18 would --
 * but s stays in c->seeds, so seedcov/MAPQ are untouched.
 *
 * Byte-identity: the skip set is a proven subset of PE18's purge set. (a) seed
 * containment implies containment in the longer seed's (extension-grown) alnreg,
 * so PE18 would purge it; (b) we replicate PE18's interference guard (a
 * comparably-long seed overlapping s on a *different* diagonal forces
 * extension) using seed coordinates; (c) dropping s as a container for other
 * seeds is safe because the longest same-diagonal seed (always extended) is a
 * superset container.
 * ---------------------------------------------------------------------- */
static inline int mem_seed_ext_redundant(const mem_chain_t *c, int si)
{
    const mem_seed_t *s = &c->seeds[si];
    long sd = (long)s->rbeg - s->qbeg;
    int has_container = 0, j;
    for (j = 0; j < c->n; ++j) {
        if (j == si) continue;
        const mem_seed_t *t = &c->seeds[j];
        if (t->len <= s->len) continue;                    /* strictly longer */
        if ((long)t->rbeg - t->qbeg != sd) continue;       /* same diagonal */
        if (s->qbeg >= t->qbeg && s->qbeg + s->len <= t->qbeg + t->len) { has_container = 1; break; }
    }
    if (!has_container) return 0;
    /* interference guard: mirror PE18 (a seed >= .95*len overlapping s on a
     * different diagonal by >= s->len/4 means s could lead to a distinct aln). */
    for (j = 0; j < c->n; ++j) {
        if (j == si) continue;
        const mem_seed_t *u = &c->seeds[j];
        if (u->len < s->len * .95) continue;
        if (s->qbeg <= u->qbeg && s->qbeg + s->len - u->qbeg >= s->len >> 2 &&
            u->qbeg - s->qbeg != u->rbeg - s->rbeg) return 0;
        if (u->qbeg <= s->qbeg && u->qbeg + u->len - s->qbeg >= s->len >> 2 &&
            s->qbeg - u->qbeg != s->rbeg - u->rbeg) return 0;
    }
    return 1;
}


void mem_chain2aln_across_reads_V2(const mem_opt_t *opt, const bntseq_t *bns,
                                   const uint8_t *pac, bseq1_t *seq_, int nseq,
                                   mem_chain_v* chain_ar, mem_alnreg_v *av_v,
                                   mem_cache *mmc, uint8_t *ref_string, int tid)
{
    SeqPair *seqPairArrayAux      = mmc->seqPairArrayAux[tid];
    SeqPair *seqPairArrayLeft128  = mmc->seqPairArrayLeft128[tid];
    SeqPair *seqPairArrayRight128 = mmc->seqPairArrayRight128[tid];
    int64_t *wsize_pair = &(mmc->wsize[tid]);

    uint8_t *seqBufLeftRef  = mmc->seqBufLeftRef[tid*CACHE_LINE];
    uint8_t *seqBufRightRef = mmc->seqBufRightRef[tid*CACHE_LINE];
    uint8_t *seqBufLeftQer  = mmc->seqBufLeftQer[tid*CACHE_LINE];
    uint8_t *seqBufRightQer = mmc->seqBufRightQer[tid*CACHE_LINE];
    int64_t *wsize_buf_ref = &(mmc->wsize_buf_ref[tid*CACHE_LINE]);
    int64_t *wsize_buf_qer = &(mmc->wsize_buf_qer[tid*CACHE_LINE]);

    // int32_t *lim_g = mmc->lim + (BATCH_SIZE + 32) * tid;
    int32_t *lim_g = mmc->lim[tid];

    mem_seed_t *s;
    int64_t l_pac = bns->l_pac, rmax[8] __attribute__((aligned(64)));
    // std::vector<int8_t> nexitv(nseq, 0);

    int numPairsLeft = 0, numPairsRight = 0;
    int numPairsLeft1 = 0, numPairsRight1 = 0;
    int numPairsLeft128 = 0, numPairsRight128 = 0;
    int numPairsLeft16 = 0, numPairsRight16 = 0;

    int64_t leftRefOffset = 0, rightRefOffset = 0;
    int64_t leftQerOffset = 0, rightQerOffset = 0;

    // ungapped fast-path threshold, computed once per call from
    // the scoring model. See ungapped_fastpath_walk docstring.
    const int fp_o_min = opt->o_del < opt->o_ins ? opt->o_del : opt->o_ins;
    const int fp_e_min = opt->e_del < opt->e_ins ? opt->e_del : opt->e_ins;
    const int fp_denom = opt->a + opt->b - fp_e_min;
    const int fp_x_threshold = (fp_denom > 0) ? (fp_o_min / fp_denom) : -1;
    // (fp_o_min, fp_e_min above are reused directly by ungapped_analyze.)

    int srt_size = MAX_SEEDS_PER_READ, fac = FAC;
    uint64_t *srt = (uint64_t *) malloc(srt_size * 8);
    uint32_t *srtgg = (uint32_t*) malloc(nseq * SEEDS_PER_READ * fac * sizeof(uint32_t));

    int spos = 0;

    /* D3 (--meth, PR-4): extension scores the ORIGINAL (unconverted) read bases
     * against the original reference, not the projected (C→T / G→A) read used for
     * seeding. seq_[l].meth_orig_seq holds the original bases as NUL-terminated
     * ASCII in the SAME orientation/order as seq_[l].seq (see the bseq1_t
     * orientation contract in bwa.h), so a 2-bit encode of it indexes identically
     * to the projected query the extension index math already assumes — every
     * qs[i] = query[...] below stays correct with no index change. We 2-bit-encode
     * into this reusable scratch buffer per read. NULL/unused outside --meth. */
    uint8_t *meth_qbuf = NULL;
    int      meth_qbuf_cap = 0;

    // uint64_t timUP = __rdtsc();
    for (int l=0; l<nseq; l++)
    {
        int max = 0;
        uint8_t *rseq = 0;

        uint32_t *srtg = srtgg;
        lim_g[l+1] = 0;

        const uint8_t *query = (uint8_t *) seq_[l].seq;
        int l_query = seq_[l].l_seq;

        /* D3 (--meth, PR-4): swap the extension query to the ORIGINAL read bases.
         * seq_[l].seq is the projected (already C→T/G→A) read, used only for
         * seeding; extension must score the original read against the original
         * ref with the per-hypothesis asymmetric matrix. We 2-bit-encode
         * meth_orig_seq (ASCII, same orientation as seq) into meth_qbuf and point
         * `query` at it. If meth_orig_seq is missing (should not happen under
         * --meth once ingest populates it) we fall back to the projected read so
         * the path still runs. */
        if (opt->meth_mode && seq_[l].meth_orig_seq != NULL) {
            if (l_query > meth_qbuf_cap) {
                meth_qbuf_cap = l_query;
                /* CodeRabbit: temp pointer so a NULL realloc doesn't leak. */
                uint8_t *meth_qbuf_new = (uint8_t *) realloc(meth_qbuf, (size_t) meth_qbuf_cap);
                assert(meth_qbuf_new != NULL);
                meth_qbuf = meth_qbuf_new;
            }
            const char *os = seq_[l].meth_orig_seq;
            for (int i = 0; i < l_query; ++i) {
                unsigned char c = (unsigned char) os[i];
                meth_qbuf[i] = (c < 4) ? c : nst_nt4_table[c];
            }
            query = meth_qbuf;
        }

        mem_chain_v *chn = &chain_ar[l];
        mem_alnreg_v *av = &av_v[l];  // alignment
        mem_chain_t *c;

        _mm_prefetch((const char*) query, _MM_HINT_NTA);

        // aln mem allocation
        av->m = 0;
        for (int j=0; j<chn->n; j++) {
            c = &chn->a[j]; av->m += c->n;
        }
        av->a = (mem_alnreg_t*)calloc(av->m, sizeof(mem_alnreg_t));

        // aln mem allocation ends
        for (int j=0; j<chn->n; j++)
        {
            c = &chn->a[j];
            assert(c->seqid == l);

            int64_t tmp = 0;
            if (c->n == 0) continue;

            _mm_prefetch((const char*) (srtg + spos + 64), _MM_HINT_NTA);
            _mm_prefetch((const char*) (lim_g), _MM_HINT_NTA);

            // get the max possible span
            rmax[0] = l_pac<<1; rmax[1] = 0;

            int chain_max_n_hits = 1;
            int64_t cb_lo=0, cb_hi=0; int cb_f=1, chain_band=0;
            for (int i = 0; i < c->n; ++i) {
                int64_t b, e;
                const mem_seed_t *t = &c->seeds[i];
                b = t->rbeg - (t->qbeg + cal_max_gap(opt, t->qbeg));
                e = t->rbeg + t->len + ((l_query - t->qbeg - t->len) +
                                        cal_max_gap(opt, l_query - t->qbeg - t->len));

                tmp = rmax[0];
                rmax[0] = tmp < b? rmax[0] : b;
                rmax[1] = (rmax[1] > e)? rmax[1] : e;
                if (t->len > max) max = t->len;
                if (t->n_hits > chain_max_n_hits) chain_max_n_hits = t->n_hits;
                if (opt->band_start > 0) { int64_t _d = t->rbeg - t->qbeg; if (cb_f){cb_lo=cb_hi=_d; cb_f=0;} else { if(_d<cb_lo)cb_lo=_d; if(_d>cb_hi)cb_hi=_d; } }
            }
            int64_t cb_span = cb_hi - cb_lo;
            if (cb_span > opt->w) cb_span = opt->w;
            chain_band = (int)cb_span;

            rmax[0] = rmax[0] > 0? rmax[0] : 0;
            rmax[1] = rmax[1] < l_pac<<1? rmax[1] : l_pac<<1;
            if (rmax[0] < l_pac && l_pac < rmax[1])
            {
                if (c->seeds[0].rbeg < l_pac) rmax[1] = l_pac;
                else rmax[0] = l_pac;
            }

            /* retrieve the reference sequence */
            {
                int rid = 0;
                // free rseq
                rseq = bns_fetch_seq_v2(bns, pac, &rmax[0],
                                        c->seeds[0].rbeg,
                                        &rmax[1], &rid, ref_string,
                                        (uint8_t*) seqPairArrayAux);
                assert(c->rid == rid);
            }

            _mm_prefetch((const char*) rseq, _MM_HINT_NTA);
            // _mm_prefetch((const char*) rseq + 64, _MM_HINT_NTA);

            // assert(c->n < MAX_SEEDS_PER_READ);  // temp
            if (c->n > srt_size) {
                srt_size = c->n + 10;
                srt = (uint64_t *) realloc(srt, srt_size * 8);
            }

            for (int i = 0; i < c->n; ++i)
                srt[i] = (uint64_t)c->seeds[i].score<<32 | i;

            if (c->n > 1)
                ks_introsort_64(c->n, srt);

            // assert((spos + c->n) < SEEDS_PER_READ * FAC * nseq);
            if ((spos + c->n) > SEEDS_PER_READ * fac * nseq) {
                fac <<= 1;
                srtgg = (uint32_t *) realloc(srtgg, nseq * SEEDS_PER_READ * fac * sizeof(uint32_t));
            }

            for (int i = 0; i < c->n; ++i)
                srtg[spos++] = srt[i];

            lim_g[l+1] += c->n;

            // uint64_t tim = __rdtsc();
            for (int k=c->n-1; k >= 0; k--)
            {
                s = &c->seeds[(uint32_t)srt[k]];

                mem_alnreg_t *a;
                // a = kv_pushp(mem_alnreg_t, *av);
                a = &av->a[av->n++];
                /* av->a was calloc'd to av->m entries just above; av->n only
                 * ever increments (no slot is reused), so each slot is already
                 * zero on first write — the per-slot memset was redundant. */

                s->aln = av->n-1;

                a->w = opt->w;
                a->score = a->truesc = -1;
                a->rid = c->rid;
                a->frac_rep = c->frac_rep;
                a->seedlen0 = s->len;
                a->c = c; //ptr
                /* D3 (--meth, PR-3): carry the OT/OB hypothesis from chain to
                 * alnreg (-1 outside --meth) so PR-4's asymmetric extension and
                 * the output (XG) layer can select the matrix/strand. */
                a->meth_hypothesis = c->meth_hypothesis;
                /* D3 (--meth, fix): strand-adjusted hypothesis for the extension
                 * matrix. Reverse-strand seeds (s->rbeg >= l_pac, original
                 * doubled-pac) extend against the RC reference window, which
                 * flips the conversion's freed cell — so flip the hypothesis. */
                a->meth_strand_hyp = (c->meth_hypothesis < 0) ? -1
                    : (int8_t)((c->meth_hypothesis ^ (s->rbeg >= l_pac ? 1 : 0)) & 1);
                /* D3 (--meth): the seed matched in projected 3-letter space, but
                 * its 4-letter score against the ORIGINAL read + ref is not
                 * necessarily len*a — a seed-internal variant (e.g. a C/T mirror
                 * that collapsed in seed space) scores as a mismatch under the
                 * per-strand matrix. Recompute the true ungapped seed score so the
                 * alnreg score (the h0 the extension starts from, hence AS, MAPQ,
                 * and mate rescue) is honest rather than optimistically len*a.
                 * --meth-scoring collapsed frees C/T (and G/A) both ways, so this
                 * reproduces len*a (bwameth-like); genomic penalizes the variant.
                 * query is the ORIGINAL read (meth_orig_seq, set above); rseq is
                 * the original 4-letter ref window. mat[ref*5+read] is target-major
                 * (matches mem_opt_fill_meth_mat). */
                int meth_seed_sc = s->len * opt->a;
                if (opt->meth_mode && a->meth_strand_hyp >= 0) {
                    const int8_t *seed_mat = mem_opt_meth_mat(opt, a->meth_strand_hyp);
                    const int64_t roff = s->rbeg - rmax[0];
                    meth_seed_sc = 0;
                    for (int _i = 0; _i < s->len; ++_i)
                        meth_seed_sc += seed_mat[ rseq[roff + _i] * 5 + query[s->qbeg + _i] ];
                }
                a->chain_n_hits = chain_max_n_hits;
                a->rb = a->qb = a->re = a->qe = H0_;

                tprof[PE19][tid] ++;

                /* NB: gated off under --meth. The dominance argument assumes a
                 * longer same-diagonal seed's extension covers the shorter one;
                 * meth's asymmetric C->T matrix breaks this (a contained seed can
                 * extend to a differently-scored aln), so the skip is NOT
                 * byte-identical under --meth (measured ~0.17% of pairs diverge). */
                if (opt->skip_contained_ext && !opt->meth_mode &&
                    mem_seed_ext_redundant(c, (uint32_t)srt[k])) {
                    a->qb = a->qe = -1;   /* pre-purge exactly as PE18 would */
                    continue;
                }

                int flag = 0;
                std::pair<int, int> pr;
                if (s->qbeg)  // left extension
                {
                    SeqPair sp;
                    sp.h0 = meth_seed_sc;   /* D3: true seed score (== len*a outside --meth) */
                    sp.seqid = c->seqid;
                    sp.regid = av->n - 1;

                    if (numPairsLeft >= *wsize_pair) {
                        fprintf(stderr, "[0000][%0.4d] Re-allocating seqPairArrays, in Left\n", tid);
                        *wsize_pair +=  1024;
                        // assert(*wsize_pair > numPairsLeft);
                        *wsize_pair += numPairsLeft + 1024;
                        seqPairArrayAux = (SeqPair *) realloc(seqPairArrayAux,
                                                              (*wsize_pair + MAX_LINE_LEN)
                                                              * sizeof(SeqPair));
                        mmc->seqPairArrayAux[tid] = seqPairArrayAux;
                        seqPairArrayLeft128 = (SeqPair *) realloc(seqPairArrayLeft128,
                                                                  (*wsize_pair + MAX_LINE_LEN)
                                                                  * sizeof(SeqPair));
                        mmc->seqPairArrayLeft128[tid] = seqPairArrayLeft128;
                        seqPairArrayRight128 = (SeqPair *) realloc(seqPairArrayRight128,
                                                                   (*wsize_pair + MAX_LINE_LEN)
                                                                   * sizeof(SeqPair));
                        mmc->seqPairArrayRight128[tid] = seqPairArrayRight128;
                    }


                    sp.idq = leftQerOffset;
                    sp.idr = leftRefOffset;

                    leftQerOffset += s->qbeg;
                    if (leftQerOffset >= *wsize_buf_qer)
                    {
                        fprintf(stderr, "[%0.4d] Re-allocating (doubling) seqBufQers in %s (left)\n",
                                tid, __func__);
                        int64_t tmp = *wsize_buf_qer;
                        *wsize_buf_qer *= 2;
                        assert(*wsize_buf_qer > leftQerOffset);
                        
                        uint8_t *seqBufQer_ = (uint8_t*)
                            _mm_realloc(seqBufLeftQer, tmp, *wsize_buf_qer, sizeof(uint8_t));
                        mmc->seqBufLeftQer[tid*CACHE_LINE] = seqBufLeftQer = seqBufQer_;

                        seqBufQer_ = (uint8_t*)
                            _mm_realloc(seqBufRightQer, tmp, *wsize_buf_qer, sizeof(uint8_t));
                        mmc->seqBufRightQer[tid*CACHE_LINE] = seqBufRightQer = seqBufQer_;
                    }

                    uint8_t *qs = seqBufLeftQer + sp.idq;
                    for (int i = 0; i < s->qbeg; ++i) qs[i] = query[s->qbeg - 1 - i];

                    tmp = s->rbeg - rmax[0];
                    leftRefOffset += tmp;
                    if (leftRefOffset >= *wsize_buf_ref)
                    {
                        fprintf(stderr, "[%0.4d] Re-allocating (doubling) seqBufRefs in %s (left)\n",
                                tid, __func__);
                        int64_t tmp = *wsize_buf_ref;
                        *wsize_buf_ref *= 2;
                        uint8_t *seqBufRef_ = (uint8_t*)
                            _mm_realloc(seqBufLeftRef, tmp, *wsize_buf_ref, sizeof(uint8_t));
                        mmc->seqBufLeftRef[tid*CACHE_LINE] = seqBufLeftRef = seqBufRef_;

                        seqBufRef_ = (uint8_t*)
                            _mm_realloc(seqBufRightRef, tmp, *wsize_buf_ref, sizeof(uint8_t));
                        mmc->seqBufRightRef[tid*CACHE_LINE] = seqBufRightRef = seqBufRef_;
                    }

                    uint8_t *rs = seqBufLeftRef + sp.idr;
                    for (int64_t i = 0; i < tmp; ++i) rs[i] = rseq[tmp - 1 - i]; //seq1

                    sp.len2 = s->qbeg;
                    sp.len1 = tmp;
                    sp.tight_band = 0; sp.chain_band = chain_band;  // 0 sentinel: "no tight bound known"

                    // ungapped analysis.
                    //   HIT      → skip SW; fill a->* from ungapped.
                    //   TIGHT    → save sp.tight_band; SW will use it.
                    //   FALLBACK → use opt->w.
                    /* D3 (--meth, PR-4): the ungapped fast path scores with
                     * opt->a/opt->b hardcoded (ungapped_analyze) — it CANNOT
                     * express the asymmetric OT/OB matrix, so a HIT would
                     * mis-penalize every bisulfite conversion (read T at ref C)
                     * as a mismatch and commit a wrong a->score. Disable it under
                     * --meth so all extension flows through the matrix-aware SW
                     * kernel. (Perf-only fast path; correctness-neutral to skip.) */
                    if (!opt->meth_mode && sp.len1 >= sp.len2 && sp.len2 <= FP_N_MAX) {
                        tprof[UGP_L_ATTEMPT][tid]++;
                        int fp_score, fp_qle, fp_gscore, fp_gtle, fp_band;
                        int fp_st = ungapped_analyze(qs, rs, sp.len2,
                                                     sp.h0, opt->a, opt->b,
                                                     fp_o_min, fp_e_min,
                                                     fp_x_threshold, opt->w,
                                                     &fp_score, &fp_qle,
                                                     &fp_gscore, &fp_gtle,
                                                     &fp_band);
                        if (fp_st == FP_STATUS_HIT) {
#if BWAMEM3_UGP_PROFILE
                            tprof[UGP_L_HIT][tid]++;
                            tprof[UGP_L_UNGAPPED][tid]++;
                            tprof[UGP_SCORE_HIST_BASE + 0 * UGP_SCORE_HIST_NBINS
                                  + ugp_score_bin(fp_score)][tid]++;
                            /* Q3: cat0=ALL, cat1=UNGAP_FINAL, cat3=HIT.
                             * For HIT, fp_score is both the would-be ungapped
                             * score and the committed score, so the delta
                             * (perfect_score - aln_score) is the same for
                             * both histograms. */
                            {
                                int _perfect = sp.h0 + opt->a * sp.len2;
                                int _delta = _perfect - fp_score;
                                int _bin = ugp_delta_bin(_delta);
                                tprof[UGP_L_CAT_UNG_BASE + 0 * UGP_CAT_NBINS + _bin][tid]++;
                                tprof[UGP_L_CAT_UNG_BASE + 1 * UGP_CAT_NBINS + _bin][tid]++;
                                tprof[UGP_L_CAT_UNG_BASE + 3 * UGP_CAT_NBINS + _bin][tid]++;
                                tprof[UGP_L_CAT_FIN_BASE + 0 * UGP_CAT_NBINS + _bin][tid]++;
                                tprof[UGP_L_CAT_FIN_BASE + 1 * UGP_CAT_NBINS + _bin][tid]++;
                                tprof[UGP_L_CAT_FIN_BASE + 3 * UGP_CAT_NBINS + _bin][tid]++;
                            }
#endif
                            // Roll back the qs/rs buffer offsets we just
                            // consumed; the batch won't reference them.
                            leftQerOffset -= s->qbeg;
                            leftRefOffset -= tmp;
                            // Mirror post-SW extraction (~line 2560).
                            a->score = fp_score;
                            if (fp_gscore <= 0 || fp_gscore <= a->score - opt->pen_clip5) {
                                a->qb = s->qbeg - fp_qle;
                                a->rb = s->rbeg - fp_qle; // tle == qle on diagonal
                                a->truesc = a->score;
                            } else {
                                a->qb = 0;
                                a->rb = s->rbeg - fp_gtle;
                                a->truesc = fp_gscore;
                            }
                            a->w = max_(a->w, opt->w);
                            if (a->rb != H0_ && a->qb != H0_ && a->qe != H0_ && a->re != H0_) {
                                int ii;
                                for (ii = 0, a->seedcov = 0; ii < a->c->n; ++ii) {
                                    const mem_seed_t *t = &(a->c->seeds[ii]);
                                    if (t->qbeg >= a->qb && t->qbeg + t->len <= a->qe &&
                                        t->rbeg >= a->rb && t->rbeg + t->len <= a->re)
                                        a->seedcov += t->len;
                                }
                            }
                            goto LEFT_DONE;
                        }
                        if (fp_st == FP_STATUS_TIGHT) {
                            // Fast-path failed but tight_band valid. Pipe
                            // to SW via sp.tight_band.
                            sp.tight_band = fp_band;
                            tprof[UGP_L_TIGHT][tid]++;
                            if (fp_band <= 8)        tprof[UGP_L_TB_1_8][tid]++;
                            else if (fp_band <= 32)  tprof[UGP_L_TB_9_32][tid]++;
                            else                     tprof[UGP_L_TB_33_MAX][tid]++;

                        }
                        // FALLBACK: sp.tight_band stays 0.
                    }

#if BWAMEM3_UGP_PROFILE
                    /* Per-pair tier routing + tight_band histograms. The
                     * numPairsLeft* increments here are dead (sortPairsLenExt
                     * below resets and recomputes all three counters, and
                     * nothing reads them before that) — only t_tier feeds the
                     * UGP histograms, so the whole block is instrumentation. */
                    {
                        int64_t minval = (int64_t)sp.h0 + (int64_t)min_(sp.len1, sp.len2) * (int64_t)opt->a;
                        int t_tier;
                        /* Mirror the routing decision in sortPairsLenExt: 8-bit
                         * iff the safe envelope holds (initial band opt->w). */
                        if (bsw8_envelope_ok(sp.len1, sp.len2, opt->w, opt->a, opt->zdrop, sp.h0)) {
                            numPairsLeft128++; t_tier = 0;
                        }
                        else if (sp.len1 < MAX_SEQ_LEN16 && sp.len2 < MAX_SEQ_LEN16 && minval >= 0 && minval < MAX_SEQ_LEN16){
                            numPairsLeft16++;  t_tier = 1;
                        }
                        else {
                            numPairsLeft1++;   t_tier = 2;  /* scalar; not bucketed */
                        }
                        /* tight_band histogram bin. */
                        int tb_ = sp.tight_band;
                        int fine_bin;
                        if      (tb_ ==  0) fine_bin = 0;
                        else if (tb_ <=  2) fine_bin = 1;
                        else if (tb_ <=  4) fine_bin = 2;
                        else if (tb_ <=  8) fine_bin = 3;
                        else if (tb_ <= 16) fine_bin = 4;
                        else if (tb_ <= 32) fine_bin = 5;
                        else if (tb_ <= 48) fine_bin = 6;
                        else if (tb_ <= 64) fine_bin = 7;
                        else if (tb_ <= 80) fine_bin = 8;
                        else                fine_bin = 9;
                        tprof[UGP_FINE_BASE + 0 * UGP_FINE_NBINS + fine_bin][tid]++;
                        if (t_tier <= 1) {
                            int band_bin;
                            if      (tb_ ==  0) band_bin = 0;
                            else if (tb_ <=  8) band_bin = 1;
                            else if (tb_ <= 32) band_bin = 2;
                            else                band_bin = 3;
                            tprof[UGP_TIER_TB_BASE + 0 * 8 + t_tier * 4 + band_bin][tid]++;
                        }
                    }
#endif

#if BWAMEM3_UGP_PROFILE
                    /* Q3: compute would-be ungapped extension score for this
                     * non-HIT LEFT pair. The walk handles arbitrary N
                     * (including the bypass case len2 > FP_N_MAX where
                     * analyze did not run). Cost: O(len2) scalar; instrumentation
                     * only (feeds ugp_record_left_outcome's tprof histograms). */
                    sp.ugp_walk_score = ungapped_walk_score(qs, rs, sp.len2,
                                                            sp.h0, opt->a, opt->b);
#endif

                    seqPairArrayLeft128[numPairsLeft] = sp;
                    numPairsLeft ++;
                    a->qb = s->qbeg; a->rb = s->rbeg;
                    LEFT_DONE: ;
                }
                else
                {
                    flag = 1;
                    a->score = a->truesc = meth_seed_sc, a->qb = 0, a->rb = s->rbeg;  /* D3: true seed score */
                }

                if (s->qbeg + s->len != l_query)  // right extension
                {
                    int64_t qe = s->qbeg + s->len;
                    int64_t re = s->rbeg + s->len - rmax[0];
                    assert(re >= 0);
                    SeqPair sp;

                    sp.h0 = H0_; //random number
                    sp.seqid = c->seqid;
                    sp.regid = av->n - 1;

                    if (numPairsRight >= *wsize_pair)
                    {
                        fprintf(stderr, "[0000] [%0.4d] Re-allocating seqPairArrays Right\n", tid);
                        *wsize_pair += 1024;
                        // assert(*wsize_pair > numPairsRight);
                        *wsize_pair += numPairsLeft + 1024;
                        seqPairArrayAux = (SeqPair *) realloc(seqPairArrayAux,
                                                              (*wsize_pair + MAX_LINE_LEN)
                                                              * sizeof(SeqPair));
                        mmc->seqPairArrayAux[tid] = seqPairArrayAux;
                        seqPairArrayLeft128 = (SeqPair *) realloc(seqPairArrayLeft128,
                                                                  (*wsize_pair + MAX_LINE_LEN)
                                                                  * sizeof(SeqPair));
                        mmc->seqPairArrayLeft128[tid] = seqPairArrayLeft128;
                        seqPairArrayRight128 = (SeqPair *) realloc(seqPairArrayRight128,
                                                                   (*wsize_pair + MAX_LINE_LEN)
                                                                   * sizeof(SeqPair));
                        mmc->seqPairArrayRight128[tid] = seqPairArrayRight128;
                    }

                    sp.len2 = l_query - qe;
                    sp.len1 = rmax[1] - rmax[0] - re;

                    sp.idq = rightQerOffset;
                    sp.idr = rightRefOffset;

                    rightQerOffset += sp.len2;
                    if (rightQerOffset >= *wsize_buf_qer)
                    {
                        fprintf(stderr, "[%0.4d] Re-allocating (doubling) seqBufQers in %s (right)\n",
                                tid, __func__);
                        int64_t tmp = *wsize_buf_qer;
                        *wsize_buf_qer *= 2;
                        assert(*wsize_buf_qer > rightQerOffset);
                        
                        uint8_t *seqBufQer_ = (uint8_t*)
                            _mm_realloc(seqBufLeftQer, tmp, *wsize_buf_qer, sizeof(uint8_t));
                        mmc->seqBufLeftQer[tid*CACHE_LINE] = seqBufLeftQer = seqBufQer_;

                        seqBufQer_ = (uint8_t*)
                            _mm_realloc(seqBufRightQer, tmp, *wsize_buf_qer, sizeof(uint8_t));
                        mmc->seqBufRightQer[tid*CACHE_LINE] = seqBufRightQer = seqBufQer_;
                    }

                    rightRefOffset += sp.len1;
                    if (rightRefOffset >= *wsize_buf_ref)
                    {
                        fprintf(stderr, "[%0.4d] Re-allocating (doubling) seqBufRefs in %s (right)\n",
                                tid, __func__);
                        int64_t tmp = *wsize_buf_ref;
                        *wsize_buf_ref *= 2;
                        uint8_t *seqBufRef_ = (uint8_t*)
                            _mm_realloc(seqBufLeftRef, tmp, *wsize_buf_ref, sizeof(uint8_t));
                        mmc->seqBufLeftRef[tid*CACHE_LINE] = seqBufLeftRef = seqBufRef_;

                        seqBufRef_ = (uint8_t*)
                            _mm_realloc(seqBufRightRef, tmp, *wsize_buf_ref, sizeof(uint8_t));
                        mmc->seqBufRightRef[tid*CACHE_LINE] = seqBufRightRef = seqBufRef_;
                    }

                    tprof[PE23][tid] += sp.len1 + sp.len2;

                    uint8_t *qs = seqBufRightQer + sp.idq;
                    uint8_t *rs = seqBufRightRef + sp.idr;

                    for (int i = 0; i < sp.len2; ++i) qs[i] = query[qe + i];

                    for (int i = 0; i < sp.len1; ++i) rs[i] = rseq[re + i]; //seq1

                    sp.tight_band = 0; sp.chain_band = chain_band;
                    sp.ugp_r_attempted = 0;

                    // ungapped analysis on right ext.
                    // Precondition a->score != -1 means left either didn't
                    // need extension or was fast-pathed — in both cases
                    // h0 is known. If left was batched, we skip (can't
                    // know h0 yet); SW will run with default band.
                    /* D3 (--meth, PR-4): disable the ungapped fast path under
                     * --meth — see the LEFT-side rationale above (it can't score
                     * the asymmetric OT/OB matrix). */
                    if (!opt->meth_mode && a->score != -1 && sp.len1 >= sp.len2 && sp.len2 <= FP_N_MAX) {
                        sp.ugp_r_attempted = 1;
                        tprof[UGP_R_ATTEMPT][tid]++;
                        int fp_h0 = a->score;  // the real h0 for right ext
                        int fp_score, fp_qle, fp_gscore, fp_gtle, fp_band;
                        int fp_st = ungapped_analyze(qs, rs, sp.len2,
                                                     fp_h0, opt->a, opt->b,
                                                     fp_o_min, fp_e_min,
                                                     fp_x_threshold, opt->w,
                                                     &fp_score, &fp_qle,
                                                     &fp_gscore, &fp_gtle,
                                                     &fp_band);
                        if (fp_st == FP_STATUS_HIT) {
                            tprof[UGP_R_HIT][tid]++;
                            tprof[UGP_R_UNGAPPED][tid]++;
                            tprof[UGP_SCORE_HIST_BASE + 1 * UGP_SCORE_HIST_NBINS
                                  + ugp_score_bin(fp_score)][tid]++;
                            // Roll back the qs/rs buffer offsets.
                            rightQerOffset -= sp.len2;
                            rightRefOffset -= sp.len1;
                            // Mirror post-SW extraction (~line 2777).
                            a->score = fp_score;
                            if (fp_gscore <= 0 || fp_gscore <= a->score - opt->pen_clip3) {
                                a->qe = qe + fp_qle;
                                a->re = rmax[0] + re + fp_qle; // tle == qle on diagonal
                                a->truesc += a->score - fp_h0;
                            } else {
                                a->qe = l_query;
                                a->re = rmax[0] + re + fp_gtle;
                                a->truesc += fp_gscore - fp_h0;
                            }
                            a->w = max_(a->w, opt->w);
                            if (a->rb != H0_ && a->qb != H0_ && a->qe != H0_ && a->re != H0_) {
                                int ii;
                                for (ii = 0, a->seedcov = 0; ii < a->c->n; ++ii) {
                                    const mem_seed_t *t = &a->c->seeds[ii];
                                    if (t->qbeg >= a->qb && t->qbeg + t->len <= a->qe &&
                                        t->rbeg >= a->rb && t->rbeg + t->len <= a->re)
                                        a->seedcov += t->len;
                                }
                            }
                            goto RIGHT_DONE;
                        }
                        if (fp_st == FP_STATUS_TIGHT) {
                            sp.tight_band = fp_band;
                            tprof[UGP_R_TIGHT][tid]++;
                            if (fp_band <= 8)        tprof[UGP_R_TB_1_8][tid]++;
                            else if (fp_band <= 32)  tprof[UGP_R_TB_9_32][tid]++;
                            else                     tprof[UGP_R_TB_33_MAX][tid]++;

                        }
                    }

                    /* The RIGHT-side tier counters (numPairsRight128/16/1) are
                     * recomputed from scratch by sortPairsLenExt below (which
                     * resets them to 0 first), and nothing reads them before
                     * that call — so the per-pair envelope check and routing
                     * that used to live here were dead work. Just stage the
                     * pair; the sort assigns its tier. RIGHT histograms (Groups
                     * A & B) are populated post-left-SW once the right pass has
                     * finalised sp->tight_band. */
                    seqPairArrayRight128[numPairsRight] = sp;
                    numPairsRight ++;
                    a->qe = qe; a->re = rmax[0] + re;
                    RIGHT_DONE: ;
                }
                else
                {
                    a->qe = l_query, a->re = s->rbeg + s->len;
                    // seedcov business, this "if" block should be redundant, check and remove.
                    if (a->rb != H0_ && a->qb != H0_)
                    {
                        int i;
                        for (i = 0, a->seedcov = 0; i < c->n; ++i)
                        {
                            const mem_seed_t *t = &c->seeds[i];
                            if (t->qbeg >= a->qb && t->qbeg + t->len <= a->qe &&
                                t->rbeg >= a->rb && t->rbeg + t->len <= a->re) // seed fully contained
                                a->seedcov += t->len;
                        }
                    }
                }
            }
            // tprof[MEM_ALN2_DOWN1][tid] += __rdtsc() - tim;
        }
    }
    // tprof[MEM_ALN2_UP][tid] += __rdtsc() - timUP;


    int32_t *hist = (int32_t *)_mm_malloc((MAX_SEQ_LEN8 + MAX_SEQ_LEN16 + 32) *
                                          sizeof(int32_t), 64);

    /* Sorting based score is required as that affects the use of SIMD lanes */
    sortPairsLenExt(seqPairArrayLeft128, numPairsLeft, seqPairArrayAux, hist,
                    numPairsLeft128, numPairsLeft16, numPairsLeft1, opt->a,
                    opt->w, opt->zdrop);
    assert(numPairsLeft == (numPairsLeft128 + numPairsLeft16 + numPairsLeft1));

#if BWAMEM3_UGP_PROFILE
    /* instrumentation: per-batch narrow bucket size (Group C, LEFT).
     * Counts pairs in the 8-bit tier with tight_band ∈ [1,8] for THIS worker
     * batch — the size a future narrow bucket would have. */
    {
        int narrow_sz = 0;
        for (int l = 0; l < numPairsLeft128; l++) {
            int tb = seqPairArrayLeft128[l].tight_band;
            if (tb > 0 && tb <= 8) narrow_sz++;
        }
        int bin;
        if      (narrow_sz ==   0) bin = 0;
        else if (narrow_sz <   16) bin = 1;
        else if (narrow_sz <   32) bin = 2;
        else if (narrow_sz <   64) bin = 3;
        else if (narrow_sz <  128) bin = 4;
        else if (narrow_sz <  256) bin = 5;
        else if (narrow_sz <  512) bin = 6;
        else                       bin = 7;
        tprof[UGP_NARROW_SZ_BASE + 0 * UGP_NARROW_SZ_NBINS + bin][tid]++;
    }
#endif


    // SWA
    // uint64_t timL = __rdtsc();
    int nthreads = 1;

    /* D3 (--meth, A1/PR-4): per-batch substitution-matrix objects. The banded-SW
     * kernel reads this->mat once per getScores16/getScores8/scalar call for ALL
     * pairs it processes (the matrix is per-SW-object, not per-pair).
     *
     * The symmetric objects (bswLeft/bswRight, built from opt->mat) drive the
     * entire non-meth pipeline and, under --meth, the non-meth (hypothesis < 0)
     * pairs. Under --meth we additionally build per-hypothesis OT/OB objects;
     * bsw_run_tier() partitions each tier's pairs by hypothesis and scores each
     * partition against its matching object, so a MIXED-hypothesis batch (the
     * common paired-end case: R1 OT + R2 OB interleaved) is scored correctly in
     * original space — the remaining PR-4 work, replacing the PR-3 symmetric
     * fallback. A homogeneous batch (e.g. single-end directional input) yields a
     * single non-empty partition, so there is no extra cost there.
     *
     * Outside --meth only the symmetric objects are built and used, and the path
     * is byte-for-byte identical to baseline. */
    // Now, process all the collected seq-pairs
    // First, left alignment, move out these calls
    auto bswLeft = make_banded_pair_wise_sw(
        opt->o_del, opt->e_del, opt->o_ins, opt->e_ins,
        opt->zdrop, opt->pen_clip5,
        opt->mat, opt->a, opt->b, nthreads);

    auto bswRight = make_banded_pair_wise_sw(
        opt->o_del, opt->e_del, opt->o_ins, opt->e_ins,
        opt->zdrop, opt->pen_clip3,
        opt->mat, opt->a, opt->b, nthreads);

    /* Per-hypothesis OT/OB objects (only under --meth). Left/right differ solely
     * in the clip-extension end bonus (pen_clip5 vs pen_clip3), matching the
     * symmetric objects above. */
    std::unique_ptr<IBandedPairWiseSW> bswLeftOt, bswLeftOb, bswRightOt, bswRightOb;
    if (opt->meth_mode) {
        const int8_t *m_ot = mem_opt_meth_mat(opt, 1);   /* OT: frees C->T */
        const int8_t *m_ob = mem_opt_meth_mat(opt, 0);   /* OB: frees G->A */
        bswLeftOt  = make_banded_pair_wise_sw(opt->o_del, opt->e_del, opt->o_ins, opt->e_ins,
                                              opt->zdrop, opt->pen_clip5, m_ot, opt->a, opt->b, nthreads);
        bswLeftOb  = make_banded_pair_wise_sw(opt->o_del, opt->e_del, opt->o_ins, opt->e_ins,
                                              opt->zdrop, opt->pen_clip5, m_ob, opt->a, opt->b, nthreads);
        bswRightOt = make_banded_pair_wise_sw(opt->o_del, opt->e_del, opt->o_ins, opt->e_ins,
                                              opt->zdrop, opt->pen_clip3, m_ot, opt->a, opt->b, nthreads);
        bswRightOb = make_banded_pair_wise_sw(opt->o_del, opt->e_del, opt->o_ins, opt->e_ins,
                                              opt->zdrop, opt->pen_clip3, m_ob, opt->a, opt->b, nthreads);
        /* OT/OB kernels score a SUB-SLICE of the pair array in bsw_run_tier, so
         * their padding-lane overshoot would corrupt the next slice: enable the
         * guard on exactly these objects. The SYM object (bswLeft/bswRight) is a
         * whole-array / last-slice caller and stays unguarded. */
        bswLeftOt->set_guard_overshoot(true);
        bswLeftOb->set_guard_overshoot(true);
        bswRightOt->set_guard_overshoot(true);
        bswRightOb->set_guard_overshoot(true);
    }

    int i;
    // Left
    SeqPair *pair_ar = seqPairArrayLeft128 + numPairsLeft128 + numPairsLeft16;
    SeqPair *pair_ar_aux = seqPairArrayAux;
    int nump = numPairsLeft1;

    // per-pair tight_band proofs are still piped in via
    // sp->tight_band (and short-circuit the retry loop below once w >=
    // tight_band), but we no longer narrow init_w from opt->w. A batched
    // SW pass shares one band across all pairs in the batch, so narrowing
    // would force FALLBACK pairs (no tight_band proof) to start with a
    // band insufficient for indels their alignment really needs. The
    // heuristic exits (a->score == prev / max_off < 3w/4) can then fire
    // on a suboptimal alignment found within the narrow band, breaking
    // chr22 parity.
    int init_w = INIT_W(opt->w);

    // scalar
    for ( i=0; i<MAX_BAND_TRY; i++)
    {
        int32_t w = init_w << i;
        // uint64_t tim = __rdtsc();
        bsw_run_tier(BSW_TIER_SCALAR, opt, av_v,
                     bswLeft.get(), bswLeftOt.get(), bswLeftOb.get(),
                     pair_ar, seqBufLeftRef, seqBufLeftQer, nump, nthreads, w,
                     pair_ar_aux, pair_ar_aux, hist);  /* CodeRabbit: base_scratch must
                     * track pair_ar_aux (distinct from pair_ar after the retry swap),
                     * not the static seqPairArrayAux which pair_ar can alias. */
        // tprof[PE5][0] += nump;
        // tprof[PE6][0] ++;
        // tprof[MEM_ALN2_B][tid] += __rdtsc() - tim;

        int num = 0;
        for (int l=0; l<nump; l++)
        {
            mem_alnreg_t *a;
            SeqPair *sp = &pair_ar[l];
            a = &(av_v[sp->seqid].a[sp->regid]); // prev
            int prev = a->score;
            a->score = sp->score;

            if (ACCEPT_PAIR(a->score, prev, sp->max_off, w, sp->tight_band, sp->chain_band, i))
            {
#if BWAMEM3_UGP_PROFILE
                ugp_record_left_outcome(sp, opt->a, tid);
#endif
                if (sp->gscore <= 0 || sp->gscore <= a->score - opt->pen_clip5) {
                    a->qb -= sp->qle; a->rb -= sp->tle;
                    a->truesc = a->score;
                } else {
                    a->qb = 0; a->rb -= sp->gtle;
                    a->truesc = sp->gscore;
                }

                a->w = max_(a->w, w);
                if (a->rb != H0_ && a->qb != H0_ && a->qe != H0_ && a->re != H0_)
                {
                    int i = 0;
                    for (i = 0, a->seedcov = 0; i < a->c->n; ++i){
                        const mem_seed_t *t = &(a->c->seeds[i]);
                        if (t->qbeg >= a->qb && t->qbeg + t->len <= a->qe &&
                            t->rbeg >= a->rb && t->rbeg + t->len <= a->re)
                            a->seedcov += t->len;
                    }
                }

            } else {
                pair_ar_aux[num++] = *sp;
            }
        }
        nump = num;
        SeqPair *tmp = pair_ar;
        pair_ar = pair_ar_aux;
        pair_ar_aux = tmp;
    }


    //****************** Left - vector int16 ***********************
    assert(numPairsLeft == (numPairsLeft128 + numPairsLeft16 + numPairsLeft1));

    pair_ar = seqPairArrayLeft128 + numPairsLeft128;
    pair_ar_aux = seqPairArrayAux;

    nump = numPairsLeft16;
    init_w = INIT_W(opt->w);
    for ( i=0; i<MAX_BAND_TRY; i++)
    {
        int32_t w = init_w << i;
        // int64_t tim = __rdtsc();
        bsw_run_tier(BSW_TIER_16, opt, av_v,
                     bswLeft.get(), bswLeftOt.get(), bswLeftOb.get(),
                     pair_ar, seqBufLeftRef, seqBufLeftQer, nump, nthreads, w,
                     pair_ar_aux, pair_ar_aux, hist);  /* CodeRabbit: base_scratch must
                     * track pair_ar_aux (distinct from pair_ar after the retry swap),
                     * not the static seqPairArrayAux which pair_ar can alias. */

        tprof[PE5][0] += nump;
        tprof[PE6][0] ++;
        // tprof[MEM_ALN2_B][tid] += __rdtsc() - tim;

        int num = 0;
        for (int l=0; l<nump; l++)
        {
            mem_alnreg_t *a;
            SeqPair *sp = &pair_ar[l];
            a = &(av_v[sp->seqid].a[sp->regid]); // prev

            int prev = a->score;
            a->score = sp->score;


            if (ACCEPT_PAIR(a->score, prev, sp->max_off, w, sp->tight_band, sp->chain_band, i))
            {
#if BWAMEM3_UGP_PROFILE
                ugp_record_left_outcome(sp, opt->a, tid);
#endif
                if (sp->gscore <= 0 || sp->gscore <= a->score - opt->pen_clip5) {
                    a->qb -= sp->qle; a->rb -= sp->tle;
                    a->truesc = a->score;
                } else {
                    a->qb = 0; a->rb -= sp->gtle;
                    a->truesc = sp->gscore;
                }

                a->w = max_(a->w, w);
                if (a->rb != H0_ && a->qb != H0_ && a->qe != H0_ && a->re != H0_)
                {
                    int i = 0;
                    for (i = 0, a->seedcov = 0; i < a->c->n; ++i){
                        const mem_seed_t *t = &(a->c->seeds[i]);
                        if (t->qbeg >= a->qb && t->qbeg + t->len <= a->qe &&
                            t->rbeg >= a->rb && t->rbeg + t->len <= a->re)
                            a->seedcov += t->len;
                    }
                }
            } else {
                pair_ar_aux[num++] = *sp;
            }
        }
        nump = num;
        SeqPair *tmp = pair_ar;
        pair_ar = pair_ar_aux;
        pair_ar_aux = tmp;
    }

    //****************** Left - vector int8 ***********************
    pair_ar = seqPairArrayLeft128;
    pair_ar_aux = seqPairArrayAux;

    nump = numPairsLeft128;
    init_w = opt->w;
    for ( i=0; i<MAX_BAND_TRY; i++)
    {
        int32_t w = init_w << i;
        // int64_t tim = __rdtsc();

        /* These pairs were bucketed for the 8-bit path at the initial band
         * opt->w (which the envelope proved <= BSW8_MAX_W). Once the doubling
         * retry pushes w >= 128 the diagonal-offset int8 encoding can no longer
         * represent the band, so bsw_run_tier (BSW_TIER_8) diverts that
         * iteration to the 16-bit kernel rather than feeding getScores8 an
         * unrepresentable band. The tight_band early-exit below is unaffected
         * (it short-circuits the retry regardless of which width ran). */
        bsw_run_tier(BSW_TIER_8, opt, av_v,
                     bswLeft.get(), bswLeftOt.get(), bswLeftOb.get(),
                     pair_ar, seqBufLeftRef, seqBufLeftQer, nump, nthreads, w,
                     pair_ar_aux, pair_ar_aux, hist);  /* CodeRabbit: base_scratch must
                     * track pair_ar_aux (distinct from pair_ar after the retry swap),
                     * not the static seqPairArrayAux which pair_ar can alias. */

        tprof[PE1][0] += nump;
        tprof[PE2][0] ++;
        // tprof[MEM_ALN2_D][tid] += __rdtsc() - tim;

        int num = 0;
        for (int l=0; l<nump; l++)
        {
            mem_alnreg_t *a;
            SeqPair *sp = &pair_ar[l];
            a = &(av_v[sp->seqid].a[sp->regid]); // prev

            int prev = a->score;
            a->score = sp->score;

            if (ACCEPT_PAIR(a->score, prev, sp->max_off, w, sp->tight_band, sp->chain_band, i))
            {
#if BWAMEM3_UGP_PROFILE
                ugp_record_left_outcome(sp, opt->a, tid);
#endif
                if (sp->gscore <= 0 || sp->gscore <= a->score - opt->pen_clip5) {
                    a->qb -= sp->qle; a->rb -= sp->tle;
                    a->truesc = a->score;
                } else {
                    a->qb = 0; a->rb -= sp->gtle;
                    a->truesc = sp->gscore;
                }

                a->w = max_(a->w, w);
                if (a->rb != H0_ && a->qb != H0_ && a->qe != H0_ && a->re != H0_)
                {
                    int i = 0;
                    for (i = 0, a->seedcov = 0; i < a->c->n; ++i){
                        const mem_seed_t *t = &(a->c->seeds[i]);
                        if (t->qbeg >= a->qb && t->qbeg + t->len <= a->qe &&
                            t->rbeg >= a->rb && t->rbeg + t->len <= a->re)
                            a->seedcov += t->len;
                    }
                }
            } else {
                pair_ar_aux[num++] = *sp;
            }
        }
        nump = num;
        SeqPair *tmp = pair_ar;
        pair_ar = pair_ar_aux;
        pair_ar_aux = tmp;
    }

    // tprof[CLEFT][tid] += __rdtsc() - timL;

    // uint64_t timR = __rdtsc();
    // **********************************************************
    // h0 fixup: right-extension h0 is the left-extension alnreg score.
    for (int l=0; l<numPairsRight; l++) {
        mem_alnreg_t *a;
        SeqPair *sp = &seqPairArrayRight128[l];
        a = &(av_v[sp->seqid].a[sp->regid]); // prev
        sp->h0 = a->score;
    }

    // post-left-SW right-side ungapped pass. The per-seed fast-
    // path skipped right ext for 93% of pairs because a->score was still
    // -1 at construction time (left SW hadn't run yet). Now that left
    // SW is done, a->score is final and h0 is set above. Run ungapped
    // analysis on every right pair:
    //   HIT    → fill a->* directly, compact the pair out of the array.
    //   TIGHT  → set sp->tight_band for the SW dispatch.
    //   FALL   → leave tight_band at construction default (typically 0).
    {
        int compacted = 0;
        for (int l = 0; l < numPairsRight; l++) {
            SeqPair *sp = &seqPairArrayRight128[l];
            mem_alnreg_t *a = &(av_v[sp->seqid].a[sp->regid]);
            int keep = 1;
            // Skip pairs already analyzed at construction time (a->score
            // was non-(-1) then; UGP_R_ATTEMPT and outcome counters fired
            // there). Without this guard, those pairs would double-count.
            // !opt->meth_mode: ungapped_analyze hardcodes symmetric opt->a/opt->b
            // and cannot honor the per-strand asymmetric matrix, so a HIT here
            // would fill a->score / qe / re from symmetric scoring and compact
            // the pair out, bypassing the asymmetric banded SW. The two
            // construction-time passes are gated the same way (see ~3397/3618).
            if (!opt->meth_mode && !sp->ugp_r_attempted &&
                sp->len1 >= sp->len2 && sp->len2 > 0 && sp->len2 <= FP_N_MAX) {
                tprof[UGP_R_ATTEMPT][tid]++;
                const uint8_t *qs = seqBufRightQer + sp->idq;
                const uint8_t *rs = seqBufRightRef + sp->idr;
                int fp_score, fp_qle, fp_gscore, fp_gtle, fp_band;
                int fp_st = ungapped_analyze(qs, rs, sp->len2, sp->h0,
                                              opt->a, opt->b,
                                              fp_o_min, fp_e_min,
                                              fp_x_threshold, opt->w,
                                              &fp_score, &fp_qle,
                                              &fp_gscore, &fp_gtle,
                                              &fp_band);
                if (fp_st == FP_STATUS_HIT) {
                    tprof[UGP_R_HIT][tid]++;
                    tprof[UGP_R_UNGAPPED][tid]++;
                    tprof[UGP_SCORE_HIST_BASE + 1 * UGP_SCORE_HIST_NBINS
                          + ugp_score_bin(fp_score)][tid]++;
                    // Mirror post-right-SW extraction. a->qe / a->re were
                    // set in the per-seed construction loop to the pre-
                    // extension anchor positions (s->qbeg+s->len / s->rbeg+s->len).
                    int l_query = seq_[sp->seqid].l_seq;
                    a->score = fp_score;
                    if (fp_gscore <= 0 || fp_gscore <= a->score - opt->pen_clip3) {
                        a->qe += fp_qle;
                        a->re += fp_qle; // tle == qle on diagonal
                        a->truesc += a->score - sp->h0;
                    } else {
                        a->qe = l_query;
                        a->re += fp_gtle;
                        a->truesc += fp_gscore - sp->h0;
                    }
                    a->w = max_(a->w, opt->w);
                    if (a->rb != H0_ && a->qb != H0_ && a->qe != H0_ && a->re != H0_) {
                        int ii;
                        for (ii = 0, a->seedcov = 0; ii < a->c->n; ++ii) {
                            const mem_seed_t *t = &a->c->seeds[ii];
                            if (t->qbeg >= a->qb && t->qbeg + t->len <= a->qe &&
                                t->rbeg >= a->rb && t->rbeg + t->len <= a->re)
                                a->seedcov += t->len;
                        }
                    }
                    keep = 0;
                } else if (fp_st == FP_STATUS_TIGHT) {
                    sp->tight_band = fp_band;
                    tprof[UGP_R_TIGHT][tid]++;
                    if (fp_band <= 8)        tprof[UGP_R_TB_1_8][tid]++;
                    else if (fp_band <= 32)  tprof[UGP_R_TB_9_32][tid]++;
                    else                     tprof[UGP_R_TB_33_MAX][tid]++;

                }
            }
            if (keep) {
                if (compacted != l)
                    seqPairArrayRight128[compacted] = *sp;
                compacted++;
            }
        }
        numPairsRight = compacted;
    }

    sortPairsLenExt(seqPairArrayRight128, numPairsRight, seqPairArrayAux,
                    hist, numPairsRight128, numPairsRight16, numPairsRight1, opt->a,
                    opt->w, opt->zdrop);

    assert(numPairsRight == (numPairsRight128 + numPairsRight16 + numPairsRight1));

#if BWAMEM3_UGP_PROFILE
    /* instrumentation (Groups A, B, C — RIGHT). Run after the post-
     * left-SW analyze pass + tier sort: sp->tight_band is now final and the
     * 128 region is contiguous at the head of seqPairArrayRight128. */
    {
        int narrow_sz = 0;
        for (int l = 0; l < numPairsRight; l++) {
            SeqPair *sp_ = &seqPairArrayRight128[l];
            int tb_ = sp_->tight_band;
            int t_tier;
            int64_t minval_ = (int64_t)sp_->h0 + (int64_t)min_(sp_->len1, sp_->len2) * (int64_t)opt->a;
            if (bsw8_envelope_ok(sp_->len1, sp_->len2, opt->w, opt->a, opt->zdrop, sp_->h0))             t_tier = 0;
            else if (sp_->len1 < MAX_SEQ_LEN16 && sp_->len2 < MAX_SEQ_LEN16 && minval_ >= 0 && minval_ < MAX_SEQ_LEN16) t_tier = 1;
            else                                                                                       t_tier = 2;
            int fine_bin;
            if      (tb_ ==  0) fine_bin = 0;
            else if (tb_ <=  2) fine_bin = 1;
            else if (tb_ <=  4) fine_bin = 2;
            else if (tb_ <=  8) fine_bin = 3;
            else if (tb_ <= 16) fine_bin = 4;
            else if (tb_ <= 32) fine_bin = 5;
            else if (tb_ <= 48) fine_bin = 6;
            else if (tb_ <= 64) fine_bin = 7;
            else if (tb_ <= 80) fine_bin = 8;
            else                fine_bin = 9;
            tprof[UGP_FINE_BASE + 1 * UGP_FINE_NBINS + fine_bin][tid]++;
            if (t_tier <= 1) {
                int band_bin;
                if      (tb_ ==  0) band_bin = 0;
                else if (tb_ <=  8) band_bin = 1;
                else if (tb_ <= 32) band_bin = 2;
                else                band_bin = 3;
                tprof[UGP_TIER_TB_BASE + 1 * 8 + t_tier * 4 + band_bin][tid]++;
            }
            if (t_tier == 0 && tb_ > 0 && tb_ <= 8) narrow_sz++;
        }
        int bin;
        if      (narrow_sz ==   0) bin = 0;
        else if (narrow_sz <   16) bin = 1;
        else if (narrow_sz <   32) bin = 2;
        else if (narrow_sz <   64) bin = 3;
        else if (narrow_sz <  128) bin = 4;
        else if (narrow_sz <  256) bin = 5;
        else if (narrow_sz <  512) bin = 6;
        else                       bin = 7;
        tprof[UGP_NARROW_SZ_BASE + 1 * UGP_NARROW_SZ_NBINS + bin][tid]++;
    }
#endif

    pair_ar = seqPairArrayRight128 + numPairsRight128 + numPairsRight16;
    pair_ar_aux = seqPairArrayAux;
    nump = numPairsRight1;
    init_w = INIT_W(opt->w);

    for ( i=0; i<MAX_BAND_TRY; i++)
    {
        int32_t w = init_w << i;
        // tim = __rdtsc();
        bsw_run_tier(BSW_TIER_SCALAR, opt, av_v,
                     bswRight.get(), bswRightOt.get(), bswRightOb.get(),
                     pair_ar, seqBufRightRef, seqBufRightQer, nump, nthreads, w,
                     pair_ar_aux, pair_ar_aux, hist);  /* CodeRabbit: base_scratch must
                     * track pair_ar_aux (distinct from pair_ar after the retry swap),
                     * not the static seqPairArrayAux which pair_ar can alias. */
        // tprof[PE7][0] += nump;
        // tprof[PE8][0] ++;
        // tprof[MEM_ALN2_C][tid] += __rdtsc() - tim;
        int num = 0;

        for (int l=0; l<nump; l++)
        {
            mem_alnreg_t *a;
            SeqPair *sp = &pair_ar[l];
            a = &(av_v[sp->seqid].a[sp->regid]); // prev
            int prev = a->score;
            a->score = sp->score;

            // no further banding
            if (ACCEPT_PAIR(a->score, prev, sp->max_off, w, sp->tight_band, sp->chain_band, i))
            {
#if BWAMEM3_UGP_PROFILE
                ugp_record_right_outcome(sp, tid);
#endif
                if (sp->gscore <= 0 || sp->gscore <= a->score - opt->pen_clip3) {
                    a->qe += sp->qle, a->re += sp->tle;
                    a->truesc += a->score - sp->h0;
                } else {
                    int l_query = seq_[sp->seqid].l_seq;
                    a->qe = l_query, a->re += sp->gtle;
                    a->truesc += sp->gscore - sp->h0;
                }
                a->w = max_(a->w, w);
                if (a->rb != H0_ && a->qb != H0_ && a->qe != H0_ && a->re != H0_)
                {
                    int i = 0;
                    for (i = 0, a->seedcov = 0; i < a->c->n; ++i) {
                        const mem_seed_t *t = &a->c->seeds[i];
                        if (t->qbeg >= a->qb && t->qbeg + t->len <= a->qe &&
                            t->rbeg >= a->rb && t->rbeg + t->len <= a->re)
                            a->seedcov += t->len;
                    }
                }
            } else {
                pair_ar_aux[num++] = *sp;
            }
        }
        nump = num;
        SeqPair *tmp = pair_ar;
        pair_ar = pair_ar_aux;
        pair_ar_aux = tmp;
    }

    // ************************* Right - vector int16 **********************
    pair_ar = seqPairArrayRight128 + numPairsRight128;
    pair_ar_aux = seqPairArrayAux;
    nump = numPairsRight16;
    init_w = INIT_W(opt->w);

    for ( i=0; i<MAX_BAND_TRY; i++)
    {
        int32_t w = init_w << i;
        // uint64_t tim = __rdtsc();
        bsw_run_tier(BSW_TIER_16, opt, av_v,
                     bswRight.get(), bswRightOt.get(), bswRightOb.get(),
                     pair_ar, seqBufRightRef, seqBufRightQer, nump, nthreads, w,
                     pair_ar_aux, pair_ar_aux, hist);  /* CodeRabbit: base_scratch must
                     * track pair_ar_aux (distinct from pair_ar after the retry swap),
                     * not the static seqPairArrayAux which pair_ar can alias. */

        tprof[PE7][0] += nump;
        tprof[PE8][0] ++;
        // tprof[MEM_ALN2_C][tid] += __rdtsc() - tim;

        int num = 0;

        for (int l=0; l<nump; l++)
        {
            mem_alnreg_t *a;
            SeqPair *sp = &pair_ar[l];
            a = &(av_v[sp->seqid].a[sp->regid]); // prev
            //OutScore *o = outScoreArray + l;
            int prev = a->score;
            a->score = sp->score;

            // no further banding
            if (ACCEPT_PAIR(a->score, prev, sp->max_off, w, sp->tight_band, sp->chain_band, i))
            {
#if BWAMEM3_UGP_PROFILE
                ugp_record_right_outcome(sp, tid);
#endif
                if (sp->gscore <= 0 || sp->gscore <= a->score - opt->pen_clip3) {
                    a->qe += sp->qle, a->re += sp->tle;
                    a->truesc += a->score - sp->h0;
                } else {
                    int l_query = seq_[sp->seqid].l_seq;
                    a->qe = l_query, a->re += sp->gtle;
                    a->truesc += sp->gscore - sp->h0;
                }
                a->w = max_(a->w, w);
                if (a->rb != H0_ && a->qb != H0_ && a->qe != H0_ && a->re != H0_)
                {
                    int i = 0;
                    for (i = 0, a->seedcov = 0; i < a->c->n; ++i) {
                        const mem_seed_t *t = &a->c->seeds[i];
                        if (t->qbeg >= a->qb && t->qbeg + t->len <= a->qe &&
                            t->rbeg >= a->rb && t->rbeg + t->len <= a->re)
                            a->seedcov += t->len;
                    }
                }
            } else {
                pair_ar_aux[num++] = *sp;
            }
        }
        nump = num;
        SeqPair *tmp = pair_ar;
        pair_ar = pair_ar_aux;
        pair_ar_aux = tmp;
    }


    // ************************* Right, vector int8 **********************
    pair_ar = seqPairArrayRight128;
    pair_ar_aux = seqPairArrayAux;
    nump = numPairsRight128;
    init_w = opt->w;

    for ( i=0; i<MAX_BAND_TRY; i++)
    {
        int32_t w = init_w << i;
        // uint64_t tim = __rdtsc();

        /* See the LEFT int8 loop: bsw_run_tier (BSW_TIER_8) diverts w >= 128
         * retry iterations to the 16-bit kernel (the diagonal-offset int8 band
         * can't represent it). tight_band early-exit is preserved. */
        bsw_run_tier(BSW_TIER_8, opt, av_v,
                     bswRight.get(), bswRightOt.get(), bswRightOb.get(),
                     pair_ar, seqBufRightRef, seqBufRightQer, nump, nthreads, w,
                     pair_ar_aux, pair_ar_aux, hist);  /* CodeRabbit: base_scratch must
                     * track pair_ar_aux (distinct from pair_ar after the retry swap),
                     * not the static seqPairArrayAux which pair_ar can alias. */

        tprof[PE3][0] += nump;
        tprof[PE4][0] ++;
        // tprof[MEM_ALN2_E][tid] += __rdtsc() - tim;
        int num = 0;

        for (int l=0; l<nump; l++)
        {
            mem_alnreg_t *a;
            SeqPair *sp = &pair_ar[l];
            a = &(av_v[sp->seqid].a[sp->regid]); // prev
            //OutScore *o = outScoreArray + l;
            int prev = a->score;
            a->score = sp->score;
            // no further banding
            if (ACCEPT_PAIR(a->score, prev, sp->max_off, w, sp->tight_band, sp->chain_band, i))
            {
#if BWAMEM3_UGP_PROFILE
                ugp_record_right_outcome(sp, tid);
#endif
                if (sp->gscore <= 0 || sp->gscore <= a->score - opt->pen_clip3) {
                    a->qe += sp->qle, a->re += sp->tle;
                    a->truesc += a->score - sp->h0;
                } else {
                    int l_query = seq_[sp->seqid].l_seq;
                    a->qe = l_query, a->re += sp->gtle;
                    a->truesc += sp->gscore - sp->h0;
                }
                a->w = max_(a->w, w);
                if (a->rb != H0_ && a->qb != H0_ && a->qe != H0_ && a->re != H0_)
                {
                    int i = 0;
                    for (i = 0, a->seedcov = 0; i < a->c->n; ++i) {
                        const mem_seed_t *t = &a->c->seeds[i];
                        if (t->qbeg >= a->qb && t->qbeg + t->len <= a->qe &&
                            t->rbeg >= a->rb && t->rbeg + t->len <= a->re)
                            a->seedcov += t->len;
                    }
                }
            } else {
                pair_ar_aux[num++] = *sp;
            }
        }
        nump = num;
        SeqPair *tmp = pair_ar;
        pair_ar = pair_ar_aux;
        pair_ar_aux = tmp;
    }

    _mm_free(hist);
    // tprof[CRIGHT][tid] += __rdtsc() - timR;

    if (numPairsLeft >= *wsize_pair || numPairsRight >= *wsize_pair)
    {   // refine it!
        fprintf(stderr, "Error: Unexpected behaviour!!!\n");
        fprintf(stderr, "Error: assert failed for seqPair size, "
                "numPairsLeft: %d, numPairsRight %d, lim: %lld\nExiting.\n",
                numPairsLeft, numPairsRight, (long long)*wsize_pair);
        exit(EXIT_FAILURE);
    }
    /* Discard seeds and hence their alignemnts */

    lim_g[0] = 0;
    for (int l=1; l<nseq; l++)
        lim_g[l] += lim_g[l-1];

    // uint64_t tim = __rdtsc();
    // BATCH_SIZE is a compile-time constant (1024 on arm64, 512 elsewhere);
    // sizeof(int)*BATCH_SIZE = 4096/2048 bytes — safe on stack, skips an
    // allocator round-trip per call to this function. nseq is bounded by
    // kt_for's BATCH_SIZE-chunked work distribution (worker_bwt at the
    // call site) — assert it for future-proofing.
    assert(nseq <= BATCH_SIZE);
    int lim[BATCH_SIZE] = {0};

    for (int l=0; l<nseq; l++)
    {
        int s_start = 0, s_end = 0;
        uint32_t *srtg = srtgg + lim_g[l];

        int l_query = seq_[l].l_seq;
        mem_chain_v *chn = &chain_ar[l];
        mem_alnreg_v *av = &av_v[l];  // alignment
        mem_chain_t *c;

        for (int j=0; j<chn->n; j++)
        {
            c = &chn->a[j];
            assert(c->seqid == l);

            s_end = s_start + c->n;

            uint32_t *srt2 = srtg + s_start;
            s_start += c->n;

            int k = 0;
            for (k = c->n-1; k >= 0; k--)
            {
                s = &c->seeds[srt2[k]];
                int i = 0;
                int v = 0;
                for (i = 0; i < av->n && v < lim[l]; ++i)  // test whether extension has been made before
                {
                    mem_alnreg_t *p = &av->a[i];
                    if (p->qb == -1 && p->qe == -1) {
                        continue;
                    }

                    int64_t rd;
                    int qd, w, max_gap;
                    if (s->rbeg < p->rb || s->rbeg + s->len > p->re || s->qbeg < p->qb
                        || s->qbeg + s->len > p->qe) {
                        v++; continue; // not fully contained
                    }

                    if (s->len - p->seedlen0 > .1 * l_query) { v++; continue;}
                    // qd: distance ahead of the seed on query; rd: on reference
                    qd = s->qbeg - p->qb; rd = s->rbeg - p->rb;
                    // the maximal gap allowed in regions ahead of the seed
                    max_gap = cal_max_gap(opt, qd < rd? qd : rd);
                    w = max_gap < p->w? max_gap : p->w; // bounded by the band width
                    if (qd - rd < w && rd - qd < w) break; // the seed is "around" a previous hit
                    // similar to the previous four lines, but this time we look at the region behind
                    qd = p->qe - (s->qbeg + s->len); rd = p->re - (s->rbeg + s->len);
                    max_gap = cal_max_gap(opt, qd < rd? qd : rd);
                    w = max_gap < p->w? max_gap : p->w;
                    if (qd - rd < w && rd - qd < w) break;

                    v++;
                }

                // the seed is (almost) contained in an existing alignment;
                // further testing is needed to confirm it is not leading
                // to a different aln
                if (v < lim[l])
                {
                    for (v = k + 1; v < c->n; ++v)
                    {
                        const mem_seed_t *t;
                        if (srt2[v] == UINT_MAX) continue;
                        t = &c->seeds[srt2[v]];
                        //if (t->done == H0_) continue;  //check for interferences!!!
                        // only check overlapping if t is long enough;
                        // TODO: more efficient by early stopping
                        if (t->len < s->len * .95) continue;
                        if (s->qbeg <= t->qbeg && s->qbeg + s->len - t->qbeg >= s->len>>2 &&
                            t->qbeg - s->qbeg != t->rbeg - s->rbeg) break;
                        if (t->qbeg <= s->qbeg && t->qbeg + t->len - s->qbeg >= s->len>>2 &&
                            s->qbeg - t->qbeg != s->rbeg - t->rbeg) break;
                    }
                    if (v == c->n) {                  // no overlapping seeds; then skip extension
                        mem_alnreg_t *ar = &(av_v[l].a[s->aln]);
                        ar->qb = ar->qe = -1;         // purge the alingment
                        srt2[k] = UINT_MAX;
                        tprof[PE18][tid]++;
                        continue;
                    }
                }
                lim[l]++;
            }
        }
    }
    free(srtgg);
    free(srt);
    free(meth_qbuf);  /* D3 (--meth, PR-4): original-read 2-bit scratch; NULL outside --meth */
    // tprof[MEM_ALN2_DOWN][tid] += __rdtsc() - tim;
}
