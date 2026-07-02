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

#ifndef BWAMEM_HPP
#define BWAMEM_HPP

/* --adaptive-band start band; the chain-geometry retry expands per pair from here */
#define ADAPTIVE_BAND_START 20

#include "bwt.h"
#include "bntseq.h"
#include "bwa.h"
#include "kthread.h"
#include "macro.h"
#include "bandedSWA.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <assert.h>
#include <limits.h>
#include <math.h>
#include <vector>
#include "kstring.h"
#include "ksw.h"
#include "kvec.h"
#include "ksort.h"
#include "utils.h"
#include "macro.h"
#include "profiling.h"
#include "FMI_search.h"

#define MEM_MAPQ_COEF 30.0
#define MEM_MAPQ_MAX  60

struct __smem_i;
typedef struct __smem_i smem_i;

#define MEM_F_PE        0x2
#define MEM_F_NOPAIRING 0x4
#define MEM_F_ALL       0x8
#define MEM_F_NO_MULTI  0x10
#define MEM_F_NO_RESCUE 0x20
#define MEM_F_REF_HDR   0x100
#define MEM_F_SOFTCLIP  0x200
#define MEM_F_SMARTPE   0x400

// V17
#define MEM_F_PRIMARY5  0x800
#define MEM_F_KEEP_SUPP_MAPQ 0x1000
#define MEM_F_XB        0x2000


/* D3 (--meth): bisulfite substitution-matrix mode, selected by --meth-scoring.
 * COLLAPSED (default) frees BOTH conversion directions so C/T (and G/A) are
 * interchangeable — reproduces bwameth's collapsed 3-letter placement (a
 * drop-in). GENOMIC frees only the bisulfite conversion direction, keeping the
 * mirror cell a real mismatch — variant-aware, truthful NM/MD. */
enum mem_meth_scoring { MEM_METH_SCORING_COLLAPSED = 0, MEM_METH_SCORING_GENOMIC = 1 };

typedef enum {
    SEED_ORDER_OFF = 0,
    SEED_ORDER_GLOBAL_LONGEST,
    SEED_ORDER_LOCAL_LONGEST,
    SEED_ORDER_ABSORB_COUNT,
    SEED_ORDER_MOST_ABSORB
} seed_order_t;

typedef struct mem_opt_t {
    int a, b;               // match score and mismatch penalty
    int o_del, e_del;
    int o_ins, e_ins;
    int pen_unpaired;       // phred-scaled penalty for unpaired reads
    int pen_clip5,pen_clip3;// clipping penalty. This score is not deducted from the DP score.
    int w;                  // band width
    int zdrop;              // Z-dropoff

    uint64_t max_mem_intv;

    int T;                  // output score threshold; only affecting output
    int flag;               // see MEM_F_* macros
    int min_seed_len;       // minimum seed length
    int min_ext_len;        // seeds shorter than this are not extended (0 = off)
    int max_extend_chains;  // cap on chains extended per read: keep only the top-N by weight before banded-SW (0 = off). Opt-in speed lever; NOT byte-identical.
    int mate_concordant_window;  // --extend-mate-concordant: when max_extend_chains caps a PE read, also retain chains concordant with a mate chain within this many bp. 0 = off; -1 = auto (use the estimated proper-pair insert high bound); >0 = fixed window. Recovers the true pair's low-weight chain (mainly --meth). NOT byte-identical.
    int est_insert_high;         // runtime state (NOT a user option): upper proper-pair insert bound (pes[FR].high) estimated from data during the run, or from -I; 0 = not yet estimated. Read by the mate-concordant cap when mate_concordant_window == -1 (auto).
    seed_order_t seed_emit_order;  // --seed-order; SEED_ORDER_OFF = byte-identical
    int min_chain_weight;
    int max_chain_extend;
    float split_factor;     // split into a seed if MEM is longer than min_seed_len*split_factor
    int split_width;        // split into a seed if its occurence is smaller than this value
    int max_occ;            // skip a seed if its occurence is larger than this value
    int max_chain_gap;      // do not chain seed if it is max_chain_gap-bp away from the closest seed
    int n_threads;          // number of threads
    int64_t chunk_size;         // process chunk_size-bp sequences in a batch
    float mask_level;       // regard a hit as redundant if the overlap with another better hit is over mask_level times the min length of the two hits
    float drop_ratio;       // drop a chain if its seed coverage is below drop_ratio times the seed coverage of a better chain overlapping with the small chain
    float XA_drop_ratio;    // when counting hits for the XA tag, ignore alignments with score < XA_drop_ratio * max_score; only effective for the XA tag
    float mask_level_redun;
    float mapQ_coef_len;
    int mapQ_coef_fac;
    int max_ins;            // when estimating insert size distribution, skip pairs with insert longer than this value
    int max_matesw;         // perform maximally max_matesw rounds of mate-SW for each end
    int max_XA_hits, max_XA_hits_alt; // if there are max_hits or fewer, output them all
    int8_t mat[25];         // scoring matrix; mat[0] == 0 if unset
    /* D3 (--meth): per-hypothesis substitution matrices, built from `mat` by
     * mem_opt_fill_meth_mat() per opt->meth_scoring (target-major mat[ref*5+read],
     * ACGT order A,C,G,T,N). Each frees the bisulfite conversion cell to a MATCH
     * (+a): mat_ot frees mat[C][T]=mat[1*5+3] (top strand C→T), mat_ob frees
     * mat[G][A]=mat[2*5+0] (bottom strand G→A).
     *   GENOMIC: ONLY that cell is freed; the mirror (mat[T][C], mat[A][G]) STAYS
     *     at −b so genuine variants score as mismatches → one freed cell ⇒ the
     *     SIMD rank-1 fast path. Variant-aware, truthful NM/MD.
     *   COLLAPSED: the mirror cell is ALSO freed (two cells) so C/T and G/A are
     *     interchangeable → reproduces bwameth; uses bandedSWA's general path.
     * Outside --meth these are unused; selection is by mem_chain_t.meth_hypothesis
     * (1 = OT, 0 = OB) via mem_opt_meth_mat(). */
    int8_t mat_ot[25];      // OT (C→T) meth matrix; valid only under --meth
    int8_t mat_ob[25];      // OB (G→A) meth matrix; valid only under --meth
    int    bam_mode;        // 1 = emit BAM instead of SAM text (--bam); meth_mode implies this
    int    bam_level;       // 0..9, BGZF deflate level (0 = uncompressed)
    int    meth_mode;       // 1 = bisulfite mode (--meth); implies bam_mode
    int    meth_scoring;    // bisulfite matrix mode (--meth-scoring): MEM_METH_SCORING_{COLLAPSED,GENOMIC}
    char   meth_set_as_failed;// 'f', 'r', or 0 — flag reads on that strand 0x200
    int    meth_chimera_qc; // 1 to enable bwameth.py-style longest-M <44% chimera heuristic (default off; not in Bismark)
    int    supp_rep_hard_cap; // supp alnregs whose chain's seeds share >=this many genome hits are forced to MAPQ=0; 0 disables
    int    smem_dedup;        // 1 = dedup fully-identical SMEMs before SA expansion (--smem-dedup); 0 = off (default, byte-identical to baseline)
    int    skip_contained_ext; // 1 = skip banded-SW extension of seeds contained (same diagonal) in a longer in-chain seed (--skip-contained-ext); 0 = off. Byte-identical to baseline: the skip set is a subset of the post-extension containment purge (PE18).
    int    band_start;       // >0 = adaptive chain-geometry banding active (start band; set to ADAPTIVE_BAND_START by --adaptive-band); 0 = off (byte-identical). Long-read speed lever; no-op on the 8-bit short-read tier.
} mem_opt_t;


struct mem_alnreg_t;
// * Chaining *
typedef struct abc {
    abc() {
        done = 0;
        rbeg = qbeg = len = score = aln = 0;
        n_hits = 1;
    }
    int64_t rbeg;
    int32_t qbeg;
    int32_t len;
    int32_t score;
    int8_t done;
    int aln;
    int32_t n_hits;  // SMEM SA occurrence count this seed came from; 1 = unique
} mem_seed_t; // unaligned memory

typedef struct {
    int32_t seqid, cseed;
    int32_t n, m, first, rid;
    uint32_t w:29, kept:2, is_alt:1;
    float frac_rep;
    int64_t pos;
    mem_seed_t *seeds;
    /* D3 (--meth, PR-3): per-chain bisulfite hypothesis label, carried from the
     * seed→original remap so PR-4 can pick the asymmetric OT/OB matrix. Encoding
     * matches the seed contig parity (seed_rid & 1): 1 = OT (C→T, odd "f" seed
     * contig), 0 = OB (G→A, even "r" seed contig). -1 = not a meth chain
     * (non-meth runs leave this at -1). For directional libraries (the --meth
     * contract) each read is projected under a SINGLE hypothesis (R1→OT, R2→OB),
     * so all of a read's seeds carry the same label and no cross-hypothesis merge
     * can occur; test_and_merge is therefore NOT hypothesis-guarded.
     * NOTE (--meth directional invariant): if non-directional / dual-hypothesis-
     * per-read support is ever added, test_and_merge MUST gain a hypothesis guard
     * (it is currently safe only because each read is single-hypothesis). */
    int8_t meth_hypothesis;
} mem_chain_t;

typedef struct { size_t n, m, cc; mem_chain_t *a;  } mem_chain_v;

typedef struct mem_alnreg_t {
    // mem_alnreg_t() {c=NULL;}
    int64_t rb, re; // [rb,re): reference sequence in the alignment
    int qb, qe;     // [qb,qe): query sequence in the alignment
    int rid;        // reference seq ID
    mem_chain_t *c;
    int score;      // best local SW score
    int truesc;     // actual score corresponding to the aligned region; possibly smaller than $score
    int sub;        // 2nd best SW score
    int alt_sc;
    int csub;       // SW score of a tandem hit
    int sub_n;      // approximate number of suboptimal hits
    int w;          // actual band width used in extension
    int seedcov;    // length of regions coverged by seeds
    int secondary;  // index of the parent hit shadowing the current hit; <0 if primary
    int secondary_all;
    int seedlen0;   // length of the starting seed
    int chain_n_hits; // max SMEM SA-occurrence count across this chain's seeds (1 = no repetitive seed)
    int n_comp:30, is_alt:2; // number of sub-alignments chained together
    float frac_rep;
    uint64_t hash;
    int flg;
    /* D3 (--meth, PR-3): bisulfite hypothesis (1=OT, 0=OB, -1=non-meth),
     * propagated from the originating chain so PR-4's asymmetric extension and
     * the output (XG strand) layers can select the right matrix/strand. */
    int8_t meth_hypothesis;
    /* D3 (--meth, fix): STRAND-ADJUSTED extension hypothesis = meth_hypothesis
     * XOR is_rev. bwa-mem extends reverse-strand seeds against the reverse-
     * complemented reference window, so a conversion that the OT matrix frees
     * (ref-C x read-T) presents as the OB-freed cell (ref-G x read-A) there, and
     * vice-versa. Extension/regen MUST select the matrix from THIS field, not the
     * raw hypothesis, or reverse-strand reads (the bulk of PE) get the wrong
     * matrix, their real conversions are scored as mismatches, and the 5' end
     * soft-clips. Set at alnreg creation from the seed strand; -1 = non-meth.
     * (Output XG/XM still use the raw meth_hypothesis — the genome strand.) */
    int8_t meth_strand_hyp;
} mem_alnreg_t;

typedef struct { size_t n, m; mem_alnreg_t *a; } mem_alnreg_v;

typedef struct {
    int low, high;   // lower and upper bounds within which a read pair is considered to be properly paired
    int failed;      // non-zero if the orientation is not supported by sufficient data
    double avg, std; // mean and stddev of the insert size distribution
} mem_pestat_t;

typedef struct { // This struct is only used for the convenience of API.
    int64_t pos;     // forward strand 5'-end mapping position
    int rid;         // reference sequence index in bntseq_t; <0 for unmapped
    int flag;        // extra flag
    uint32_t is_rev:1, is_alt:1, mapq:8, NM:22; // is_rev: whether on the reverse strand; mapq: mapping quality; NM: edit distance
    int n_cigar;     // number of CIGAR operations
    uint32_t *cigar; // CIGAR in the BAM encoding: opLen<<4|op; op to integer mapping: MIDSH=>01234
    char *XA;        // alternative mappings
    int HN;          // total # of hits clustered with this primary under XA_drop_ratio; -1 when not computed (e.g., MEM_F_ALL)

    int score, sub, alt_sc;
    /* D3 (--meth, PR-5): bisulfite hypothesis (1=OT, 0=OB, -1=non-meth),
     * copied from the source mem_alnreg_t in mem_reg2aln. The output layer
     * (meth_mem_aln_to_bam) sources the XG strand tag and the XM source
     * strand from THIS, not from the (now retired) f/r contig direction. */
    int8_t meth_hypothesis;
} mem_aln_t;

// struct
typedef struct {
    bwtintv_v mem, mem1, *tmpv[2];
} smem_aux_t;

typedef struct
{
    SeqPair *seqPairArrayAux[MAX_THREADS];
    SeqPair *seqPairArrayLeft128[MAX_THREADS];
    SeqPair *seqPairArrayRight128[MAX_THREADS];
    
    int64_t wsize[MAX_THREADS];

    int64_t wsize_buf_ref[MAX_THREADS*CACHE_LINE]; 
    int64_t wsize_buf_qer[MAX_THREADS*CACHE_LINE];

    uint8_t *seqBufLeftRef[MAX_THREADS*CACHE_LINE];
    uint8_t *seqBufRightRef[MAX_THREADS*CACHE_LINE];
    uint8_t *seqBufLeftQer[MAX_THREADS*CACHE_LINE];
    uint8_t *seqBufRightQer[MAX_THREADS*CACHE_LINE];    

    SMEM *matchArray[MAX_THREADS];
    int32_t *min_intv_ar[MAX_THREADS];
    int32_t *rid[MAX_THREADS];
    int32_t *lim[MAX_THREADS];
    int32_t *query_pos_ar[MAX_THREADS];
    uint8_t *enc_qdb[MAX_THREADS];
    
    int64_t wsize_mem[MAX_THREADS];
    int64_t wsize_mem_s[MAX_THREADS];
    int64_t wsize_mem_r[MAX_THREADS];
    // enc_qdb byte-size tracker. Separate from wsize_mem because
    // mem_collect_smem repurposes wsize_mem as an SMEM-entry count, which
    // breaks the byte-count invariant for enc_qdb on short-read workloads
    // (≤50 bp aDNA), where SMEM count >> tot_len bytes can leave enc_qdb
    // silently undersized between batches.
    int64_t wsize_qdb[MAX_THREADS];

    // Lockstep SMEM batching per-slot state. One pair of contiguous SMEM
    // buffers per thread, each of size SMEM_LOCKSTEP_N * lockstep_buf_cap[tid].
    // Per-slot views are stride-offset into these at lockstep driver entry.
    // Grown on demand when the batch's max_readlength exceeds the per-slot
    // capacity.
    SMEM    *lockstep_prev[MAX_THREADS];
    SMEM    *lockstep_match_buf[MAX_THREADS];
    int64_t  lockstep_buf_cap[MAX_THREADS];

    // Pointer into worker_t::ref_string (the unpacked .0123 reference).
    // Set once in the worker_aln/worker_sam entry points; lets helpers like
    // mem_seed_sw and the mem_matesw_* family invoke bns_fetch_seq_v2 without
    // threading ref_string through 8 function signatures. Read-only and
    // shared across threads — every thread sees the same pointer.
    uint8_t *ref_string;
} mem_cache;

// chain moved to .h
typedef struct worker_t {
    const mem_opt_t      *opt;
    //const bntseq_t         *bns;
    // const uint8_t         *pac;
    const mem_pestat_t   *pes;
    smem_aux_t      **aux;
    bseq1_t          *seqs;
    mem_alnreg_v     *regs;
    int64_t           n_processed;
    mem_chain_v      *chain_ar;
    mem_cache         mmc;
    mem_seed_t       *seedBuf;
    int64_t           seedBufSize;
    mem_seed_t       *auxSeedBuf;
    int64_t           auxSeedBufSize;
    uint8_t          *ref_string;
    int16_t           nthreads;
    int32_t           nreads;
    FMI_search       *fmi;
    /* D3 (--meth, PR-3) coordinate cutover. The seed FM-index in `fmi` lives in
     * f/r-doubled SEED coordinates and is used ONLY for candidate generation;
     * every downstream consumer (chaining merge tests, extension ref fetch,
     * pairing/insert-size, mate rescue, output) must run in ORIGINAL-reference
     * coordinates. These handles are the original (un-converted) bns/pac and the
     * original unpacked `.0123` ref_string; seed hits are remapped into the
     * original-doubled-pac space in mem_chain_seeds before chaining, after which
     * these (NOT fmi->idx->bns/pac) are the bns/pac/ref_string the rest of the
     * pipeline uses. NULL outside --meth. */
    const bntseq_t   *meth_orig_bns;
    const uint8_t    *meth_orig_pac;
    uint8_t          *meth_orig_ref_string;
} worker_t;

/* D3 (--meth, PR-3) helpers. In --meth, this returns the ORIGINAL bns/pac/
 * ref_string the chaining/extension/pairing/output layers must use; outside
 * --meth it returns the single seed/normal index handles. Centralizes the
 * "which coordinate space" decision so every consumer is converted in one place
 * (the coordinate-cutover audit, B4). */
static inline const bntseq_t *mem_aln_bns(const worker_t *w) {
    return (w->opt->meth_mode && w->meth_orig_bns != NULL)
           ? w->meth_orig_bns : w->fmi->idx->bns;
}
static inline const uint8_t *mem_aln_pac(const worker_t *w) {
    return (w->opt->meth_mode && w->meth_orig_pac != NULL)
           ? w->meth_orig_pac : w->fmi->idx->pac;
}
static inline uint8_t *mem_aln_ref_string(const worker_t *w) {
    return (w->opt->meth_mode && w->meth_orig_ref_string != NULL)
           ? w->meth_orig_ref_string : w->ref_string;
}

/* D3 (--meth, PR-4): select the substitution matrix for an extension by the
 * per-chain/per-alnreg OT/OB hypothesis. Outside --meth (or hypothesis < 0,
 * i.e. a non-meth chain) this returns the symmetric `opt->mat` unchanged so the
 * non-meth path is byte-for-byte identical. hypothesis: 1 = OT (mat_ot, frees
 * C→T), 0 = OB (mat_ob, frees G→A). See mem_opt_t::mat_ot/mat_ob. */
static inline const int8_t *mem_opt_meth_mat(const mem_opt_t *opt, int meth_hypothesis) {
    if (!opt->meth_mode || meth_hypothesis < 0) return opt->mat;
    return (meth_hypothesis & 1) ? opt->mat_ot : opt->mat_ob;
}


typedef kvec_t(int) int_v;

smem_i *smem_itr_init(const bwt_t *bwt);
void smem_itr_destroy(smem_i *itr);
void smem_set_query(smem_i *itr, int len, const uint8_t *query);
void smem_config(smem_i *itr, int min_intv, int max_len, uint64_t max_intv);
const bwtintv_v *smem_next(smem_i *itr);

mem_opt_t *mem_opt_init(void);
void mem_fill_scmat(int a, int b, int8_t mat[25]);
/* (Re)derive the --meth per-hypothesis matrices (mat_ot/mat_ob) from opt->mat +
 * opt->a. Call after any rebuild of opt->mat (e.g. CLI -A/-B/-x parsing) so meth
 * scoring tracks the user's options instead of the init-time defaults. */
void mem_opt_fill_meth_mat(mem_opt_t *opt);

// Skip-short-seed extension filter: drop seeds shorter than min_ext_len from a
// chain in place (stable; surviving seeds keep their order). Returns the new
// seed count. min_ext_len <= 0 is a no-op. See mem_opt_t::min_ext_len.
int mem_chain_drop_short_seeds(mem_chain_t *c, int min_ext_len);

void mem_reg2sam(const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac,
                 bseq1_t *s, mem_alnreg_v *a, int extra_flag, const mem_aln_t *m);

int mem_approx_mapq_se(const mem_opt_t *opt, const mem_alnreg_t *a) ;

int mem_mark_primary_se(const mem_opt_t *opt, int n, mem_alnreg_t *a, int64_t id);

static void mem_mark_primary_se_core(const mem_opt_t *opt, int n, mem_alnreg_t *a, int_v *z);

char **mem_gen_alt(const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac,
                   const mem_alnreg_v *a, int l_query, const char *query,
                   int **out_hn); // ONLY work after mem_mark_primary_se(); out_hn may be NULL
void mem_aln2sam(const mem_opt_t *opt, const bntseq_t *bns, kstring_t *str, bseq1_t *s,
                 int n, const mem_aln_t *list, int which, const mem_aln_t *m_);

static inline int get_rlen(int n_cigar, const uint32_t *cigar);
static inline int infer_bw(int l1, int l2, int score, int a, int q, int r);

int mem_kernel1_core(FMI_search *fmi, const mem_opt_t *opt,
                     bseq1_t *seq_,
                     int nseq,
                     mem_chain_v *chain_ar,
                     mem_seed_t *seedBuf,
                     int64_t seedBufSize,
                     mem_cache *mmc,
                     int tid,
                     /* D3 (--meth, PR-3): ORIGINAL bns/pac remap target +
                      * filtering ref. NULL outside --meth → legacy behavior. */
                     const bntseq_t *meth_orig_bns = NULL,
                     const uint8_t  *meth_orig_pac = NULL);

int mem_kernel2_core(FMI_search *fmi, const mem_opt_t *opt,
                     bseq1_t *seq_, mem_alnreg_v *regs, int nseq,
                     mem_chain_v *chain_ar, mem_cache *mmc,
                     uint8_t *ref_string, int tid,
                     /* D3 (--meth, PR-3): ORIGINAL bns/pac for extension/dedup.
                      * NULL outside --meth → legacy behavior. ref_string above
                      * is already the original .0123 in --meth (see worker_aln). */
                     const bntseq_t *meth_orig_bns = NULL,
                     const uint8_t  *meth_orig_pac = NULL);

void* _mm_realloc(void *ptr, int64_t csize, int64_t nsize, int16_t dsize);

void mem_chain2aln_across_reads_V2(const mem_opt_t *opt, const bntseq_t *bns,
                                   const uint8_t *pac, bseq1_t *seq_, int nseq,
                                   mem_chain_v* chain_ar, mem_alnreg_v *av_v,
                                   mem_cache *mmc, uint8_t *ref_string, int tid);

int mem_sam_pe_batch_pre(const mem_opt_t *opt, const bntseq_t *bns,
                         const uint8_t *pac, const mem_pestat_t pes[4],
                         uint64_t id, bseq1_t s[2], mem_alnreg_v a[2],
                         mem_cache *mmc,  int64_t &pcnt, int32_t &gcnt,
                         int32_t&, int32_t&, int tid);

int mem_matesw_batch_pre(const mem_opt_t *opt, const bntseq_t *bns,
                         const uint8_t *pac, const mem_pestat_t pes[4],
                         const mem_alnreg_t *a, int l_ms, const uint8_t *ms,
                         mem_alnreg_v *ma, mem_cache *mmc, int pcnt, int32_t gcnt,
                         int32_t &maxRefLen, int32_t &maxQerLen, int32_t tid);

/* Given two alignment begin positions (rb) on the 2-bit-packed concat-with-
 * reverse-complement index, infer the pair orientation (0=FF, 1=FR, 2=RF,
 * 3=RR) and the distance between 5' ends. Exposed so that external consumers
 * can share the same proper-pair classification used by mem_sam_pe's
 * no_pairing fallback. */
int mem_infer_dir(int64_t l_pac, int64_t b1, int64_t b2, int64_t *dist);

int mem_sam_pe_batch(const mem_opt_t *opt, mem_cache *mmc,
                     int64_t &pcnt, int64_t &pcnt8, kswr_t *aln,
                     int32_t, int32_t, int tid);

int mem_sam_pe_batch_post(const mem_opt_t *opt, const bntseq_t *bns,
                          const uint8_t *pac, const mem_pestat_t pes[4],
                          uint64_t id, bseq1_t s[2], mem_alnreg_v a[2],
                          kswr_t **myaln, mem_cache *mmc,
                          int32_t &gcnt, int tid);

int mem_matesw_batch_post(const mem_opt_t *opt, const bntseq_t *bns,
                          const uint8_t *pac, const mem_pestat_t pes[4],
                          const mem_alnreg_t *a, int l_ms, const uint8_t *ms,
                          mem_alnreg_v *ma, kswr_t **myaln, int32_t gcnt,
                          int32_t *gar, mem_cache *mmc, const char *ms_orig = NULL,
                          const int8_t *mat = NULL);

int mem_sam_pe(const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac,
               const mem_pestat_t pes[4], uint64_t id, bseq1_t s[2],
               mem_alnreg_v a[2]);

// Core pairing decision for a single read pair. Runs mate-rescue SW,
// mem_mark_primary_se, optional MEM_F_PRIMARY5 reorder, mem_pair, is-multi
// sanity check, q_pe/q_se computation, and the secondary<->primary
// secondary_all patch. Does NOT emit SAM/BAM — callers do that.
//
// Output contract:
//   * paired_out: set to 1 iff the paired branch was taken; 0 otherwise.
//     Callers must always inspect *paired_out before reading z[] / q_se[].
//   * z[], q_se[]: valid only when *paired_out == 1. On the no-pairing path
//     they must be treated as undefined — mem_pair() can populate z[]
//     before the is_multi early return, so their contents do not signal
//     anything when *paired_out == 0.
//   * extra_flag_out: on the paired path it is fully assembled (and
//     includes 0x2 iff the paired alignment was preferred); on the
//     no-pairing path it is only partial — the caller (or mem_sam_pe's
//     no_pairing emission block) is responsible for OR-ing in 0x2 itself.
//   * n_pri[]: always populated.
//
// Returns the number of mate-rescue hits produced (same meaning as the
// historical `n` return of mem_sam_pe).
int mem_pair_resolve(const mem_opt_t *opt, const bntseq_t *bns,
                     const uint8_t *pac, const mem_pestat_t pes[4],
                     uint64_t id, bseq1_t s[2], mem_alnreg_v a[2],
                     int n_pri[2], int z[2], int q_se[2],
                     int *extra_flag_out, int *paired_out);
/**
 * Align a batch of sequences and generate the alignments in the SAM format
 *
 * This routine requires $seqs[i].{l_seq,seq,name} and write $seqs[i].sam.
 * Note that $seqs[i].sam may consist of several SAM lines if the
 * corresponding sequence has multiple primary hits.
 *
 * In the paired-end mode (i.e. MEM_F_PE is set in $opt->flag), query
 * sequences must be interleaved: $n must be an even number and the 2i-th
 * sequence and the (2i+1)-th sequence constitute a read pair. In this
 * mode, there should be enough (typically >50) unique pairs for the
 * routine to infer the orientation and insert size.
 *
 * @param opt    alignment parameters
 * @param bwt    FM-index of the reference sequence
 * @param bns    Information of the reference
 * @param pac    2-bit encoded reference
 * @param n      number of query sequences
 * @param seqs   query sequences; $seqs[i].seq/sam to be modified after the call
 * @param pes0   insert-size info; if NULL, infer from data; if not NULL, it should be an array with 4 elements,
 *               corresponding to each FF, FR, RF and RR orientation. See mem_pestat() for more info.
 */
void mem_process_seqs(mem_opt_t *opt, int64_t n_processed,
                      int n, bseq1_t *seqs, const mem_pestat_t *pes0,
                      worker_t &w);


/**
 * Generate CIGAR and forward-strand position from alignment region
 *
 * @param opt    alignment parameters
 * @param bns    Information of the reference
 * @param pac    2-bit encoded reference
 * @param l_seq  length of query sequence
 * @param seq    query sequence
 * @param ar     one alignment region
 *
 * @return       CIGAR, strand, mapping quality and forward-strand position
 */
/* D3 (--meth, PR-5): `meth_orig_query` is the read's ORIGINAL (un-projected)
 * bases (bseq1_t.meth_orig_seq), same orientation/order as `seq`. When non-NULL
 * under --meth, the final CIGAR/NM/MD regen runs the ORIGINAL read against the
 * ORIGINAL reference with the per-hypothesis asymmetric matrix
 * (mem_opt_meth_mat(opt, ar->meth_hypothesis)) so output is native-alphabet.
 * NULL (the default) preserves the legacy symmetric `seq`-vs-ref regen exactly
 * for the non-meth path. */
mem_aln_t mem_reg2aln(const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac,
                      int l_seq, const char *seq, const mem_alnreg_t *ar,
                      const char *meth_orig_query = NULL);


/**
 * Infer the insert size distribution from interleaved alignment regions
 *
 * This function can be called after mem_align1(), as long as paired-end
 * reads are properly interleaved.
 *
 * @param opt    alignment parameters
 * @param l_pac  length of concatenated reference sequence
 * @param n      number of query sequences; must be an even number
 * @param regs   region array of size $n; 2i-th and (2i+1)-th elements constitute a pair
 * @param pes    inferred insert size distribution (output)
 */
void mem_pestat(const mem_opt_t *opt, int64_t l_pac, int n, const mem_alnreg_v *regs,
                mem_pestat_t pes[4]);

void mem_reorder_primary5(int T, mem_alnreg_v *a);

#endif
