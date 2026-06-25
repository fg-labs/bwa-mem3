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

#include <iostream>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include "kstring.h"
#include "bwamem.h"
#include "bam_writer.h"
#include "kvec.h"
#include "u8vec_scratch.h"
#include "utils.h"
#include "ksw.h"
#include "bandedSWA.h"
#include "kswv.h"

#ifdef USE_MALLOC_WRAPPERS
#  include "malloc_wrap.h"
#endif

/* File-scope forward declaration of the dedup/patch entry point defined in
 * src/bwamem.cpp (kept default-free there). The `mat = NULL` default lives
 * here only, so [dcl.fct.default]/4 (no redefining a default in one TU) is
 * satisfied — the historical per-block-scope re-declarations each repeated the
 * default and only compiled under -fpermissive. */
extern int mem_sort_dedup_patch(const mem_opt_t *opt, const bntseq_t *bns,
                                const uint8_t *pac, uint8_t *query, int n,
                                mem_alnreg_t *a, const int8_t *mat = NULL);


#define MIN_RATIO     0.8
#define MIN_DIR_CNT   10
#define MIN_DIR_RATIO 0.05
#define OUTLIER_BOUND 2.0
#define MAPPING_BOUND 3.0
#define MAX_STDDEV    4.0

extern uint64_t tprof[LIM_R][LIM_C];

int mem_infer_dir(int64_t l_pac, int64_t b1, int64_t b2, int64_t *dist)
{
    int64_t p2;
    int r1 = (b1 >= l_pac), r2 = (b2 >= l_pac);
    p2 = r1 == r2? b2 : (l_pac<<1) - 1 - b2; // p2 is the coordinate of read 2 on the read 1 strand
    *dist = p2 > b1? p2 - b1 : b1 - p2;
    return (r1 == r2? 0 : 1) ^ (p2 > b1? 0 : 3);
}

static int cal_sub(const mem_opt_t *opt, mem_alnreg_v *r)
{
    int j;
    for (j = 1; j < r->n; ++j) { // choose unique alignment
        int b_max = r->a[j].qb > r->a[0].qb? r->a[j].qb : r->a[0].qb;
        int e_min = r->a[j].qe < r->a[0].qe? r->a[j].qe : r->a[0].qe;
        if (e_min > b_max) { // have overlap
            int min_l = r->a[j].qe - r->a[j].qb < r->a[0].qe - r->a[0].qb? r->a[j].qe - r->a[j].qb : r->a[0].qe - r->a[0].qb;
            if (e_min - b_max >= min_l * opt->mask_level) break; // significant overlap
        }
    }
    return j < r->n? r->a[j].score : opt->min_seed_len * opt->a;
}

void mem_pestat(const mem_opt_t *opt, int64_t l_pac, int n,
                const mem_alnreg_v *regs, mem_pestat_t pes[4])
{
    int i, d, max;
    uint64_v isize[4];
    
    memset(pes, 0, 4 * sizeof(mem_pestat_t));
    // memset(isize, 0, sizeof(kvec_t(int)) * 4);
    memset(isize, 0, sizeof(uint64_v) * 4);
    for (i = 0; i < n>>1; ++i) {
        int dir;
        int64_t is;
        mem_alnreg_v *r[2];
        r[0] = (mem_alnreg_v*)&regs[i<<1|0];
        r[1] = (mem_alnreg_v*)&regs[i<<1|1];
        if (r[0]->n == 0 || r[1]->n == 0) continue;
        if (cal_sub(opt, r[0]) > MIN_RATIO * r[0]->a[0].score) continue;
        if (cal_sub(opt, r[1]) > MIN_RATIO * r[1]->a[0].score) continue;
        if (r[0]->a[0].rid != r[1]->a[0].rid) continue; // not on the same chr
        dir = mem_infer_dir(l_pac, r[0]->a[0].rb, r[1]->a[0].rb, &is);
        if (is && is <= opt->max_ins) kv_push(uint64_t, isize[dir], is);
    }
    if (bwa_verbose >= 3) fprintf(stderr, "[0000][PE] # candidate unique pairs for (FF, FR, RF, RR): (%ld, %ld, %ld, %ld)\n", isize[0].n, isize[1].n, isize[2].n, isize[3].n);
    for (d = 0; d < 4; ++d) { // TODO: this block is nearly identical to the one in bwtsw2_pair.c. It would be better to merge these two.
        mem_pestat_t *r = &pes[d];
        uint64_v *q = &isize[d];
        int p25, p50, p75, x;
        if (q->n < MIN_DIR_CNT) {
            fprintf(stderr, "[0000][PE] skip orientation %c%c as there are not enough pairs\n", "FR"[d>>1&1], "FR"[d&1]);
            r->failed = 1;
            free(q->a);
            continue;
        } else fprintf(stderr, "[0000][PE] analyzing insert size distribution for orientation %c%c...\n", "FR"[d>>1&1], "FR"[d&1]);
        ks_introsort_64(q->n, q->a);
        p25 = q->a[(int)(.25 * q->n + .499)];
        p50 = q->a[(int)(.50 * q->n + .499)];
        p75 = q->a[(int)(.75 * q->n + .499)];
        r->low  = (int)(p25 - OUTLIER_BOUND * (p75 - p25) + .499);
        if (r->low < 1) r->low = 1;
        r->high = (int)(p75 + OUTLIER_BOUND * (p75 - p25) + .499);
        fprintf(stderr, "[0000][PE] (25, 50, 75) percentile: (%d, %d, %d)\n", p25, p50, p75);
        fprintf(stderr, "[0000][PE] low and high boundaries for computing mean and std.dev: (%d, %d)\n", r->low, r->high);
        for (i = x = 0, r->avg = 0; i < q->n; ++i)
            if (q->a[i] >= r->low && q->a[i] <= r->high)
                r->avg += q->a[i], ++x;
        assert(x != 0);
        r->avg /= x;
        for (i = 0, r->std = 0; i < q->n; ++i)
            if (q->a[i] >= r->low && q->a[i] <= r->high)
                r->std += (q->a[i] - r->avg) * (q->a[i] - r->avg);
        r->std = sqrt(r->std / x);
        fprintf(stderr, "[0000][PE] mean and std.dev: (%.2f, %.2f)\n", r->avg, r->std);
        r->low  = (int)(p25 - MAPPING_BOUND * (p75 - p25) + .499);
        r->high = (int)(p75 + MAPPING_BOUND * (p75 - p25) + .499);
        if (r->low  > r->avg - MAX_STDDEV * r->std) r->low  = (int)(r->avg - MAX_STDDEV * r->std + .499);
        if (r->high < r->avg + MAX_STDDEV * r->std) r->high = (int)(r->avg + MAX_STDDEV * r->std + .499);
        if (r->low < 1) r->low = 1;
        fprintf(stderr, "[0000][PE] low and high boundaries for proper pairs: (%d, %d)\n", r->low, r->high);
        free(q->a);
    }
    for (d = 0, max = 0; d < 4; ++d)
        max = max > isize[d].n? max : isize[d].n;
    for (d = 0; d < 4; ++d)
        if (pes[d].failed == 0 && isize[d].n < max * MIN_DIR_RATIO) {
            pes[d].failed = 1;
            fprintf(stderr, "[0000][PE] skip orientation %c%c\n", "FR"[d>>1&1], "FR"[d&1]);
        }
}

int mem_matesw(const mem_opt_t *opt, const bntseq_t *bns,
               const uint8_t *pac, const mem_pestat_t pes[4],
               const mem_alnreg_t *a, int l_ms, const uint8_t *ms,
               mem_alnreg_v *ma, const char *ms_orig = NULL,
               const int8_t *mat = NULL)
{
    /* D3 (--meth, PR-6): the mate-rescue dedup at the bottom of this function
     * still passes mat=NULL (resolving to opt->mat) — but those calls pass
     * bns/pac/query = 0 (dedup-only, no patch SW), so the matrix is unused
     * there regardless; leaving them on opt->mat is correct. The per-hypothesis
     * asymmetric scorer below uses the `mat` parameter (the caller selects the
     * OPPOSITE-strand matrix of the anchor via mem_opt_meth_mat). */
    /* Outside --meth, mat is NULL and ms_orig is NULL: behavior is byte-for-byte
     * identical to the historical symmetric/projected mate rescue. */
    const int8_t *sw_mat = mat ? mat : opt->mat;
    /* D3 (--meth, PR-6): when ms_orig is set, score the ORIGINAL (unconverted)
     * mate bases against the original ref window instead of the projected mate.
     * meth_orig_seq is ASCII in the SAME orientation as `ms`/`seq` (bwa.h
     * contract), so 2-bit-encode it here and let the existing per-orientation
     * RC below reverse-complement it EXACTLY as it does the projected mate. */
    uint8_t *ms2 = NULL;
    if (ms_orig != NULL) {
        ms2 = (uint8_t*) malloc(l_ms);
        assert(ms2 != NULL);
        for (int k = 0; k < l_ms; ++k) {
            unsigned char c = (unsigned char) ms_orig[k];
            ms2[k] = (c < 4) ? c : nst_nt4_table[c];
        }
        ms = ms2; // alias the projected-mate pointer to the original bases
    }
    #if MATE_SORT
    extern int mem_dedup_patch(const mem_opt_t *opt, const bntseq_t *bns,
                               const uint8_t *pac, uint8_t *query, int n, mem_alnreg_t *a);
    extern void sort_alnreg_re(int n, mem_alnreg_t* a);
    extern void sort_alnreg_score(int n, mem_alnreg_t* a);
    #endif
    
    //int tid = omp_get_thread_num();
    int64_t l_pac = bns->l_pac;
    int i, r, skip[4], n = 0, rid = -1;
    for (r = 0; r < 4; ++r)
        skip[r] = pes[r].failed? 1 : 0;


    for (i = 0; i < ma->n; ++i) { // check which orinentation has been found
        int64_t dist;
        r = mem_infer_dir(l_pac, a->rb, ma->a[i].rb, &dist);
        
        if (dist >= pes[r].low && dist <= pes[r].high) {
            skip[r] = 1;
        }
    }

    if (skip[0] + skip[1] + skip[2] + skip[3] == 4) return 0; // consistent pair exist; no need to perform SW

    for (r = 0; r < 4; ++r) {
        int is_rev, is_larger;
        uint8_t *seq, *rev = 0, *ref = 0;
        int64_t rb, re;
        if (skip[r]) continue;
        is_rev = (r>>1 != (r&1)); // whether to reverse complement the mate
        is_larger = !(r>>1); // whether the mate has larger coordinate
        if (is_rev) {
            rev = (uint8_t*) malloc(l_ms); // this is the reverse complement of $ms
            assert(rev != NULL);
            for (i = 0; i < l_ms; ++i) rev[l_ms - 1 - i] = ms[i] < 4? 3 - ms[i] : 4;
            seq = rev;
        } else seq = (uint8_t*)ms;
        if (!is_rev) {
            rb = is_larger? a->rb + pes[r].low : a->rb - pes[r].high;
            re = (is_larger? a->rb + pes[r].high: a->rb - pes[r].low) + l_ms; // if on the same strand, end position should be larger to make room for the seq length
        } else {
            rb = (is_larger? a->rb + pes[r].low : a->rb - pes[r].high) - l_ms; // similarly on opposite strands
            re = is_larger? a->rb + pes[r].high: a->rb - pes[r].low;
        }
        if (rb < 0) rb = 0;
        if (re > l_pac<<1) re = l_pac<<1;
        static thread_local u8vec_scratch_t t_ref;
        if (rb < re) {
            size_t want = (size_t)((re - rb) + 64);
            if (t_ref.v.m < want) kv_resize(uint8_t, t_ref.v, want);
            int64_t rlen;
            bns_fetch_seq_into(bns, pac, &rb, (rb+re)>>1, &re, &rid, t_ref.v.a, &rlen);
            ref = t_ref.v.a;
        }
        if (a->rid == rid && re - rb >= opt->min_seed_len) { // no funny things happening
            kswr_t aln;
            mem_alnreg_t b;
            int tmp, xtra = KSW_XSUBO | KSW_XSTART | (l_ms * opt->a < 250? KSW_XBYTE : 0) | (opt->min_seed_len * opt->a);

            assert(ref !=0 && re - rb >= 0);
            aln = ksw_align2(l_ms, seq, re - rb, ref, 5,
                             sw_mat, opt->o_del, opt->e_del,
                             opt->o_ins, opt->e_ins, xtra, 0);

            memset(&b, 0, sizeof(mem_alnreg_t));
            if (aln.score >= opt->min_seed_len && aln.qb >= 0) { // something goes wrong if aln.qb < 0
                b.rid = a->rid;
                b.is_alt = a->is_alt;
                /* D3 (--meth, PR-6, B3): the rescued mate's hypothesis is the
                 * OPPOSITE strand of the anchor for directional libraries
                 * (R1→OT, R2→OB): rescuing an OT anchor's mate uses OB and vice
                 * versa. The caller selects sw_mat = mem_opt_meth_mat(opt,
                 * !anchor_hyp) accordingly; record !a->meth_hypothesis here so
                 * the output layer (XG/XM) sources the right strand. -1 anchor
                 * (non-meth) stays -1. Single hypothesis per rescue ⇒ rank-1.
                 * Coordinates are already ORIGINAL (l_pac is the original l_pac
                 * via the original bns), so the 6a coordinate fix holds. */
                b.meth_hypothesis = (a->meth_hypothesis < 0) ? -1
                                                             : !a->meth_hypothesis;
                b.qb = is_rev? l_ms - (aln.qe + 1) : aln.qb;
                b.qe = is_rev? l_ms - aln.qb : aln.qe + 1;
                b.rb = is_rev? (l_pac<<1) - (rb + aln.te + 1) : rb + aln.tb;
                b.re = is_rev? (l_pac<<1) - (rb + aln.tb) : rb + aln.te + 1;
                b.score = aln.score;
                b.csub = aln.score2;
                b.secondary = -1;
                b.seedcov = (b.re - b.rb < b.qe - b.qb? b.re - b.rb : b.qe - b.qb) >> 1;
                b.chain_n_hits = 1; // mate-rescue has no SMEM evidence; treat as unique anchor

                kv_push(mem_alnreg_t, *ma, b); // make room for a new element

                #if !MATE_SORT
                // move b s.t. ma is sorted
                for (i = 0; i < ma->n - 1; ++i) // find the insertion point
                    if (ma->a[i].score < b.score) break;
                tmp = i;
                for (i = ma->n - 1; i > tmp; --i) ma->a[i] = ma->a[i-1];
                ma->a[i] = b;
                
                #else
                
                int resort = 0;
                for (i = 0; i < ma->n - 1; ++i) { // find the insertion point
                    if (ma->a[i].re == b.re) {
                        resort = 1;
                        break;
                    }
                    if (ma->a[i].re > b.re) {
                        break;
                    }
                }
                if (resort) {
                    // Don't know where to put this alignment. So let the scores decide
                    sort_alnreg_score(ma->n - 1, ma->a);
                    for (i = 0; i < ma->n - 1; ++i) { // find the insertion point
                        if (ma->a[i].score < b.score) {
                            break;
                        }
                    }
                    tmp = i;
                    for (i = ma->n - 1; i > tmp; --i) ma->a[i] = ma->a[i-1];
                    ma->a[i] = b;
                    // Now we can sort based on end position
                    sort_alnreg_re(ma->n, ma->a);
                }
                else {
                    tmp = i;
                    for (i = ma->n - 1; i > tmp; --i) ma->a[i] = ma->a[i-1];
                    ma->a[i] = b;
                }
                #endif
                tprof[PE26][0] ++;
            }
            ++n;
        }
        #if !MATE_SORT
        if (n) ma->n = mem_sort_dedup_patch(opt, 0, 0, 0, ma->n, ma->a);
        #else
        if (n) ma->n = mem_dedup_patch(opt, 0, 0, 0, ma->n, ma->a); // sam_improvements
        #endif
        if (rev) free(rev);
        /* ref aliases t_ref thread-local scratch; do not free. */
    }
    if (ms2) free(ms2); // D3 (--meth): original-mate 2-bit scratch
    return n;
}

int mem_pair(const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac, const mem_pestat_t pes[4], bseq1_t s[2], mem_alnreg_v a[2], int id, int *sub, int *n_sub, int z[2], int n_pri[2])
{
    pair64_v v, u;
    int r, i, k, y[4], ret; // y[] keeps the last hit
    int64_t l_pac = bns->l_pac;
    kv_init(v); kv_init(u);
    for (r = 0; r < 2; ++r) { // loop through read number
        for (i = 0; i < n_pri[r]; ++i) {
            pair64_t key;
            mem_alnreg_t *e = &a[r].a[i];
            key.x = e->rb < l_pac? e->rb : (l_pac<<1) - 1 - e->rb; // forward position
            key.x = (uint64_t)e->rid<<32 | (key.x - bns->anns[e->rid].offset);
            key.y = (uint64_t)e->score << 32 | i << 2 | (e->rb >= l_pac)<<1 | r;
            kv_push(pair64_t, v, key);
        }
    }

    ks_introsort_128(v.n, v.a);
    y[0] = y[1] = y[2] = y[3] = -1;
    for (i = 0; i < v.n; ++i) {
        for (r = 0; r < 2; ++r) { // loop through direction
            int dir = r<<1 | (v.a[i].y>>1&1), which;
            if (pes[dir].failed) continue; // invalid orientation
            which = r<<1 | ((v.a[i].y&1)^1);
            if (y[which] < 0) continue; // no previous hits
            for (k = y[which]; k >= 0; --k) { // TODO: this is a O(n^2) solution in the worst case; remember to check if this loop takes a lot of time (I doubt)
                int64_t dist;
                int q;
                double ns;
                pair64_t *p;
                if ((v.a[k].y&3) != which) continue;
                dist = (int64_t)v.a[i].x - v.a[k].x;
                //printf("%d: %lld\n", k, dist);
                if (dist > pes[dir].high) break;
                if (dist < pes[dir].low)  continue;
                ns = (dist - pes[dir].avg) / pes[dir].std;
                q = (int)((v.a[i].y>>32) + (v.a[k].y>>32) + .721 * log(2. * erfc(fabs(ns) * M_SQRT1_2)) * opt->a + .499); // .721 = 1/log(4)
                if (q < 0) q = 0;
                p = kv_pushp(pair64_t, u);
                p->y = (uint64_t)k<<32 | i;
                p->x = (uint64_t)q<<32 | (hash_64(p->y ^ id<<8) & 0xffffffffU);
            }
        }
        y[v.a[i].y&3] = i;
    }
    if (u.n) { // found at least one proper pair
        int tmp = opt->a + opt->b;
        tmp = tmp > opt->o_del + opt->e_del? tmp : opt->o_del + opt->e_del;
        tmp = tmp > opt->o_ins + opt->e_ins? tmp : opt->o_ins + opt->e_ins;
        ks_introsort_128(u.n, u.a);
        i = u.a[u.n-1].y >> 32; k = u.a[u.n-1].y << 32 >> 32;
        z[v.a[i].y&1] = v.a[i].y<<32>>34; // index of the best pair
        z[v.a[k].y&1] = v.a[k].y<<32>>34;
        ret = u.a[u.n-1].x >> 32;
        *sub = u.n > 1? u.a[u.n-2].x>>32 : 0;
        for (i = (long)u.n - 2, *n_sub = 0; i >= 0; --i)
            if (*sub - (int)(u.a[i].x>>32) <= tmp) ++*n_sub;

    } else ret = 0, *sub = 0, *n_sub = 0;
    free(u.a); free(v.a);
    return ret;
}

void mem_aln2sam(const mem_opt_t *opt, const bntseq_t *bns, kstring_t *str, bseq1_t *s, int n, const mem_aln_t *list, int which, const mem_aln_t *m);

#define raw_mapq(diff, a) ((int)(6.02 * (diff) / (a) + .499))


// Core pairing decision for a single read pair.
//
// Performs, in order:
//   1. Mate-rescue SW (unless MEM_F_NO_RESCUE) — may add entries to a[].
//   2. mem_mark_primary_se on both reads — mutates a[] and fills n_pri[].
//   3. mem_reorder_primary5 if MEM_F_PRIMARY5 — mutates a[].
//   4. Pairing attempt via mem_pair, is-multi sanity check, q_pe / q_se
//      computation, and the secondary<->primary secondary_all patch.
//
// On exit:
//   *paired_out = 1 if the paired branch was taken (extra_flag includes 0x2
//                   iff o > score_un, i.e. paired alignment is preferred).
//               = 0 if pairing fell through (MEM_F_NOPAIRING set, neither
//                   read has a primary, mem_pair returned 0, or is_multi).
//                   Callers must treat z[]/q_se[] as undefined in this case:
//                   mem_pair() may already have populated z[] before an
//                   early return (e.g. is_multi), so paired_out is the only
//                   validity check. Callers take the "no_pairing" emission
//                   path here.
//   extra_flag_out carries the common extra_flag bits (always has 0x1).
//                   On the paired path it is fully assembled and includes
//                   0x2 iff the paired alignment was preferred. On the
//                   no-pairing path it is only partial — the caller (or
//                   mem_sam_pe's no_pairing emission block) is responsible
//                   for OR-ing in 0x2 itself.
//   n_pri[] is always populated.
//   z[], q_se[] are valid only when *paired_out == 1.
//
// Returns the number of mate-rescue hits produced (same meaning as the
// historical `n` return from mem_sam_pe — i.e. a caller-visible accounting
// of mate-SW work done).
int mem_pair_resolve(const mem_opt_t *opt, const bntseq_t *bns,
                     const uint8_t *pac, const mem_pestat_t pes[4],
                     uint64_t id, bseq1_t s[2], mem_alnreg_v a[2],
                     int n_pri[2], int z[2], int q_se[2],
                     int *extra_flag_out, int *paired_out)
{
    extern int mem_mark_primary_se(const mem_opt_t *opt, int n, mem_alnreg_t *a, int64_t id);
    extern int mem_approx_mapq_se(const mem_opt_t *opt, const mem_alnreg_t *a);

    #if MATE_SORT
    extern void sort_alnreg_re(int n, mem_alnreg_t* a);
    extern void sort_alnreg_score(int n, mem_alnreg_t* a);
    /* D3 (--meth, PR-6): this dedup call passes bns/pac/query = 0 (dedup-only,
     * no patch SW), so its matrix is never consulted — opt->mat default is
     * correct; the per-hypothesis asymmetric scoring lives in mem_matesw. */
    #endif

    int n = 0, i, j, o, subo, n_sub, extra_flag = 1;

    *paired_out = 0;
    *extra_flag_out = extra_flag;

    if (!(opt->flag & MEM_F_NO_RESCUE)) { // then perform SW for the best alignment

        mem_alnreg_v b[2];
        kv_init(b[0]); kv_init(b[1]);
        for (i = 0; i < 2; ++i)
            for (j = 0; j < a[i].n; ++j)
                if (a[i].a[j].score >= a[i].a[0].score  - opt->pen_unpaired)
                    kv_push(mem_alnreg_t, b[i], a[i].a[j]);

        #if MATE_SORT
        for (i = 0; i < 2; ++i) {
            sort_alnreg_re(a[!i].n, a[!i].a);
            int val = 0, swcount = 0;
            for (j = 0; j < b[i].n && j < opt->max_matesw; ++j) {
                /* D3 (--meth, PR-6, B3): rescue read !i against original ref using
                 * its ORIGINAL bases and the OPPOSITE-strand matrix of anchor
                 * b[i].a[j] (R1→OT / R2→OB ⇒ mate uses the other). Outside --meth
                 * these stay NULL and mem_matesw is identical to before. */
                const char  *ms_orig = opt->meth_mode ? s[!i].meth_orig_seq : NULL;
                const int8_t *rmat    = opt->meth_mode
                    ? mem_opt_meth_mat(opt, !b[i].a[j].meth_hypothesis) : NULL;
                int val = mem_matesw(opt, bns, pac, pes, &b[i].a[j], s[!i].l_seq, (uint8_t*)s[!i].seq, &a[!i], ms_orig, rmat);
                n += val;
                swcount += val;
            }
            if (swcount > 0) {
                mem_alnreg_v* ma = &a[!i];
                ma->n = mem_sort_dedup_patch(opt, 0, 0, 0, ma->n, ma->a);
            }
            else {
                sort_alnreg_score(a[!i].n, a[!i].a);
            }
        }

        #else

        for (i = 0; i < 2; ++i)
            for (j = 0; j < b[i].n && j < opt->max_matesw; ++j) {
                /* D3 (--meth, PR-6, B3): see MATE_SORT branch above — original
                 * mate bases + opposite-strand matrix of the anchor. */
                const char  *ms_orig = opt->meth_mode ? s[!i].meth_orig_seq : NULL;
                const int8_t *rmat    = opt->meth_mode
                    ? mem_opt_meth_mat(opt, !b[i].a[j].meth_hypothesis) : NULL;
                int val = mem_matesw(opt, bns, pac, pes, &b[i].a[j], s[!i].l_seq, (uint8_t*)s[!i].seq, &a[!i], ms_orig, rmat);
                n += val;
            }
        #endif
        free(b[0].a); free(b[1].a);
    }

    n_pri[0] = mem_mark_primary_se(opt, a[0].n, a[0].a, id<<1|0);
    n_pri[1] = mem_mark_primary_se(opt, a[1].n, a[1].a, id<<1|1);

    #if V17
    if (opt->flag & MEM_F_PRIMARY5) {
        mem_reorder_primary5(opt->T, &a[0]);
        mem_reorder_primary5(opt->T, &a[1]);
    }
    #endif

    if (opt->flag & MEM_F_NOPAIRING) {
        *extra_flag_out = extra_flag;
        return n;
    }

    // pairing single-end hits
    if (!(n_pri[0] && n_pri[1] &&
          (o = mem_pair(opt, bns, pac, pes, s, a, id, &subo, &n_sub, z, n_pri)) > 0)) {
        *extra_flag_out = extra_flag;
        return n;
    }

    int is_multi[2], q_pe, score_un;
    for (i = 0; i < 2; ++i) {
        for (j = 1; j < n_pri[i]; ++j)
            if (a[i].a[j].secondary < 0 && a[i].a[j].score >= opt->T) break;
        is_multi[i] = j < n_pri[i] ? 1 : 0;
    }
    if (is_multi[0] || is_multi[1]) { // TODO: in rare cases, the true hit may be long but with low score
        *extra_flag_out = extra_flag;
        return n;
    }

    // compute mapQ for the best SE hit
    score_un = a[0].a[0].score + a[1].a[0].score - opt->pen_unpaired;
    subo = subo > score_un ? subo : score_un;
    q_pe = raw_mapq(o - subo, opt->a);
    if (n_sub > 0) q_pe -= (int)(4.343 * log(n_sub + 1) + .499);
    if (q_pe < 0) q_pe = 0;
    if (q_pe > 60) q_pe = 60;
    q_pe = (int)(q_pe * (1. - .5 * (a[0].a[0].frac_rep + a[1].a[0].frac_rep)) + .499);

    // the following assumes no split hits
    if (o > score_un) { // paired alignment is preferred
        mem_alnreg_t *c[2];
        c[0] = &a[0].a[z[0]]; c[1] = &a[1].a[z[1]];
        for (i = 0; i < 2; ++i) {
            if (c[i]->secondary >= 0)
                c[i]->sub = a[i].a[c[i]->secondary].score, c[i]->secondary = -2;
            q_se[i] = mem_approx_mapq_se(opt, c[i]);
        }
        q_se[0] = q_se[0] > q_pe ? q_se[0] : q_pe < q_se[0] + 40 ? q_pe : q_se[0] + 40;
        q_se[1] = q_se[1] > q_pe ? q_se[1] : q_pe < q_se[1] + 40 ? q_pe : q_se[1] + 40;
        extra_flag |= 2;

        // cap at the tandem repeat score
        q_se[0] = q_se[0] < raw_mapq(c[0]->score - c[0]->csub, opt->a) ? q_se[0] : raw_mapq(c[0]->score - c[0]->csub, opt->a);
        q_se[1] = q_se[1] < raw_mapq(c[1]->score - c[1]->csub, opt->a) ? q_se[1] : raw_mapq(c[1]->score - c[1]->csub, opt->a);
    } else { // the unpaired alignment is preferred
        z[0] = z[1] = 0;
        q_se[0] = mem_approx_mapq_se(opt, &a[0].a[0]);
        q_se[1] = mem_approx_mapq_se(opt, &a[1].a[0]);
    }

    for (i = 0; i < 2; ++i) {
        int k = a[i].a[z[i]].secondary_all;
        if (k >= 0 && k < n_pri[i]) { // switch secondary and primary if both of them are non-ALT
            assert(a[i].a[k].secondary_all < 0);
            for (j = 0; j < a[i].n; ++j)
                if (a[i].a[j].secondary_all == k || j == k)
                    a[i].a[j].secondary_all = z[i];
            a[i].a[z[i]].secondary_all = -1;
        }
    }

    *paired_out = 1;
    *extra_flag_out = extra_flag;
    return n;
}


int mem_sam_pe(const mem_opt_t *opt, const bntseq_t *bns,
               const uint8_t *pac, const mem_pestat_t pes[4],
               uint64_t id, bseq1_t s[2], mem_alnreg_v a[2])
{
    extern void mem_reg2sam(const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac, bseq1_t *s, mem_alnreg_v *a, int extra_flag, const mem_aln_t *m);
    extern char **mem_gen_alt(const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac, const mem_alnreg_v *a, int l_query, const char *query, int **out_hn);

    int i, j, z[2], extra_flag, n_pri[2], q_se[2], n_aa[2], paired;
    kstring_t str;
    mem_aln_t h[2], g[2], aa[2][2];

    str.l = str.m = 0; str.s = 0;
    memset(h, 0, sizeof(mem_aln_t) * 2);
    memset(g, 0, sizeof(mem_aln_t) * 2);
    n_aa[0] = n_aa[1] = 0;

    int n = mem_pair_resolve(opt, bns, pac, pes, id, s, a, n_pri, z, q_se,
                             &extra_flag, &paired);

    if (paired) {
        char **XA[2];
        int *HN[2] = { 0, 0 };
        if (!(opt->flag & MEM_F_ALL)) {
            for (i = 0; i < 2; ++i)
                XA[i] = mem_gen_alt(opt, bns, pac, &a[i], s[i].l_seq, s[i].seq, &HN[i]);
        } else XA[0] = XA[1] = 0;
        // write SAM
        for (i = 0; i < 2; ++i) {
            h[i] = mem_reg2aln(opt, bns, pac, s[i].l_seq, s[i].seq, &a[i].a[z[i]], s[i].meth_orig_seq);
            h[i].mapq = q_se[i];

            h[i].flag |= 0x40<<i | extra_flag;
            h[i].XA = XA[i]? XA[i][z[i]] : 0;
            h[i].HN = HN[i]? HN[i][z[i]] : -1;
            aa[i][n_aa[i]++] = h[i];
            if (n_pri[i] < a[i].n) { // the read has ALT hits
                mem_alnreg_t *p = &a[i].a[n_pri[i]];
                if (p->score < opt->T || p->secondary >= 0 || !p->is_alt) continue;
                g[i] = mem_reg2aln(opt, bns, pac, s[i].l_seq, s[i].seq, p, s[i].meth_orig_seq);
                g[i].flag |= 0x800 | 0x40<<i | extra_flag;
                g[i].XA = XA[i]? XA[i][n_pri[i]] : 0;
                g[i].HN = HN[i]? HN[i][n_pri[i]] : -1;
                if (opt->supp_rep_hard_cap > 0 && p->chain_n_hits >= opt->supp_rep_hard_cap)
                    g[i].mapq = 0; // fg-labs: force repetitive-supp MAPQ to 0
                aa[i][n_aa[i]++] = g[i];
            }
        }
        for (i = 0; i < n_aa[0]; ++i)
            mem_aln2sam(opt, bns, &str, &s[0], n_aa[0], aa[0], i, &h[1]);

        if (opt->bam_mode) {
            /* BAM path (meth or generic): mem_aln2sam short-circuited
             * into s->bams, leaving str untouched. Skip the str.s dance. */
            s[0].sam = NULL;
            str.l = 0;
            for (i = 0; i < n_aa[1]; ++i)
                mem_aln2sam(opt, bns, &str, &s[1], n_aa[1], aa[1], i, &h[0]);
            s[1].sam = NULL;
            free(str.s); str.s = NULL; str.m = 0;
        } else {
            assert(str.s != 0);
            s[0].sam = strdup(str.s); str.l = 0;
            for (i = 0; i < n_aa[1]; ++i)
                mem_aln2sam(opt, bns, &str, &s[1], n_aa[1], aa[1], i, &h[0]); // write read2 hits
            s[1].sam = str.s;
        }
        if (strcmp(s[0].name, s[1].name) != 0) err_fatal(__func__, "paired reads have different names: \"%s\", \"%s\"\n", s[0].name, s[1].name);
        // free
        for (i = 0; i < 2; ++i) {
            free(h[i].cigar); free(g[i].cigar);
            free(HN[i]);
            if (XA[i] == 0) continue;
            for (j = 0; j < a[i].n; ++j) free(XA[i][j]);
            free(XA[i]);
        }
        return n;
    }

    // no_pairing
    int which[2] = { -1, -1 };
    for (i = 0; i < 2; ++i) {
        if (a[i].n) {
            if (a[i].a[0].score >= opt->T) which[i] = 0;
            else if (n_pri[i] < a[i].n && a[i].a[n_pri[i]].score >= opt->T)
                which[i] = n_pri[i];
        }
        if (which[i] >= 0) h[i] = mem_reg2aln(opt, bns, pac, s[i].l_seq, s[i].seq, &a[i].a[which[i]], s[i].meth_orig_seq);
        else h[i] = mem_reg2aln(opt, bns, pac, s[i].l_seq, s[i].seq, 0);
    }
    // Proper-pair flag must be computed from the same alignments that were just
    // emitted via mem_reg2aln — i.e. a[i].a[which[i]]. Using a[i].a[0] here is
    // wrong when which[i] == n_pri[i] (below-T primary + above-T ALT case).
    if (!(opt->flag & MEM_F_NOPAIRING) && which[0] >= 0 && which[1] >= 0 &&
        h[0].rid == h[1].rid && h[0].rid >= 0) {
        int64_t dist;
        int d;
        d = mem_infer_dir(bns->l_pac, a[0].a[which[0]].rb, a[1].a[which[1]].rb, &dist);
        if (!pes[d].failed && dist >= pes[d].low && dist <= pes[d].high) extra_flag |= 2;
    }
    mem_reg2sam(opt, bns, pac, &s[0], &a[0], 0x41|extra_flag, &h[1]);
    mem_reg2sam(opt, bns, pac, &s[1], &a[1], 0x81|extra_flag, &h[0]);
    if (strcmp(s[0].name, s[1].name) != 0)
        err_fatal(__func__, "paired reads have different names: \"%s\", \"%s\"\n",
                  s[0].name, s[1].name);

    free(h[0].cigar); free(h[1].cigar);
    return n;
}

int mem_sam_pe_batch_pre(const mem_opt_t *opt, const bntseq_t *bns,
                         const uint8_t *pac, const mem_pestat_t pes[4],
                         uint64_t id, bseq1_t s[2], mem_alnreg_v a[2],
                         mem_cache *mmc, int64_t &pcnt, int32_t &gcnt,
                         int32_t &maxRefLen, int32_t &maxQerLen,
                         int tid)
{
    //uint8_t *seqBufRef = mmc->seqBufLeftRef[tid*CACHE_LINE];
    //uint8_t *seqBufQer = mmc->seqBufLeftQer[tid*CACHE_LINE];
    // int64_t *wsize_buf = &(mmc->wsize_buf[tid]);

    //SeqPair *seqPairArray = mmc->seqPairArrayLeft128[tid];
    //int32_t *gar = (int32_t*) (mmc->seqPairArrayAux[tid]);
    // int64_t *wsize = &(mmc->wsize[tid]);
    
    int i, j, n_aa[2];
    kstring_t str;
    mem_aln_t h[2], g[2];
    // int tid = omp_get_thread_num();
    
    str.l = str.m = 0; str.s = 0;
    memset(h, 0, sizeof(mem_aln_t) * 2);
    memset(g, 0, sizeof(mem_aln_t) * 2);
    n_aa[0] = n_aa[1] = 0;
    
    if (!(opt->flag & MEM_F_NO_RESCUE)) { // then perform SW for the best alignment
        mem_alnreg_v b[2];
        kv_init(b[0]); kv_init(b[1]);
        for (i = 0; i < 2; ++i)
            for (j = 0; j < a[i].n; ++j)
                if (a[i].a[j].score >= a[i].a[0].score  - opt->pen_unpaired)
                    kv_push(mem_alnreg_t, b[i], a[i].a[j]);
        
        // NEW, batching
        for (i = 0; i < 2; ++i) {
            for (j = 0; j < b[i].n && j < opt->max_matesw; ++j) {
                int64_t val = mem_matesw_batch_pre(opt, bns, pac, pes, &b[i].a[j],
                                                   s[!i].l_seq, (uint8_t*)s[!i].seq,
                                                   &a[!i], mmc, pcnt, gcnt,
                                                   maxRefLen, maxQerLen, tid);

                pcnt = val;
                gcnt += 4;
            }
        }       
        free(b[0].a); free(b[1].a);
    }
    
    return 1;
}

static inline void revseq(int l, uint8_t *s)
{
    int i, t;
    for (i = 0; i < l>>1; ++i)
        t = s[i], s[i] = s[l - 1 - i], s[l - 1 - i] = t;
}

/* Task 5: Read the BWAMEM3_METH_BATCHED_RESCUE escape hatch once. Default ON
 * (batched meth mate rescue). Set to "0" to force the legacy scalar ksw_align2
 * rescue path (used by the regression to A/B the same binary). Only consulted
 * under opt->meth_mode; the non-meth path is unaffected. */
static bool meth_batched_rescue_enabled()
{
    static const bool enabled = []() {
        const char *e = getenv("BWAMEM3_METH_BATCHED_RESCUE");
        return !(e != NULL && strcmp(e, "0") == 0);
    }();
    return enabled;
}

/* Task 5: run the two-phase batched kswv over a contiguous SeqPair slice
 * pairs[0..slice_pcnt) whose 8-bit pairs occupy [0, slice_pcnt8) and 16-bit
 * pairs [slice_pcnt8, slice_pcnt). The slice must have MAX_LINE_LEN trailing
 * headroom (the 16-bit shift writes into pairs[slice_pcnt + MAX_LINE_LEN)).
 * aln is the shared per-batch result array indexed by SeqPair.regid; seqBufRef
 * /seqBufQer are the shared (idr/idq-keyed) sequence buffers. This is the exact
 * body of the pre-Task-5 mem_sam_pe_batch, hoisted verbatim so the non-meth
 * call (one slice = the whole batch) stays byte-identical. */
static void mem_sam_pe_batch_run(Ikswv *pwsw, SeqPair *pairs,
                                 uint8_t *seqBufRef, uint8_t *seqBufQer,
                                 kswr_t *aln, int64_t slice_pcnt,
                                 int64_t slice_pcnt8, int nthreads)
{
    // Shift 16-bit
    for (int i=0; i<slice_pcnt-slice_pcnt8; i++)
        pairs[slice_pcnt + MAX_LINE_LEN - 1 - i] = pairs[slice_pcnt-i-1];

#if BWAMEM_BATCHED_MATESW
    pwsw->getScores8(pairs, seqBufRef, seqBufQer, aln, slice_pcnt8, nthreads, 0);
    pwsw->getScores16(pairs + slice_pcnt8 + MAX_LINE_LEN, seqBufRef, seqBufQer,
                      aln, slice_pcnt-slice_pcnt8, nthreads, 0);
#else
    fprintf(stderr, "Error: mem_sam_pe_batch reached without a batched kswv kernel\n");
    exit(EXIT_FAILURE);
#endif

    // Post-processing
    int pos = 0, pos8 = 0, pos16 = 0;
    for (int i=0; i<slice_pcnt8; i++)
    {
        SeqPair sp = pairs[i];
        int ind = sp.regid;
        kswr_t r = aln[ind];
        int xtra = sp.h0;
        if ((xtra & KSW_XSTART) == 0 || ((xtra & KSW_XSUBO) && r.score < (xtra & 0xffff))) continue;

        sp.h0 = KSW_XSTOP | r.score;
        sp.len2 = r.qe + 1;
        uint8_t *qs = seqBufQer + sp.idq;
        uint8_t *rs = seqBufRef + sp.idr;
        revseq(r.qe + 1, qs); revseq(r.te + 1, rs);
        pairs[pos++] = sp;
        pos8 ++;
    }

    int id = slice_pcnt8 + MAX_LINE_LEN;
    for (int i=0; i<slice_pcnt-slice_pcnt8; i++)
    {
        SeqPair sp = pairs[i + id];
        int ind = sp.regid;
        kswr_t r = aln[ind];
        int xtra = sp.h0;
        if ((xtra & KSW_XSTART) == 0 || ((xtra & KSW_XSUBO) && r.score < (xtra & 0xffff))) continue;

        sp.h0 = KSW_XSTOP | r.score;
        sp.len2 = r.qe + 1;
        uint8_t *qs = seqBufQer + sp.idq;
        uint8_t *rs = seqBufRef + sp.idr;
        revseq(r.qe + 1, qs); revseq(r.te + 1, rs);
        pairs[pos++] = sp;
        pos16 ++;
    }

    int pcnt2 = pos;
    assert(pos8 + pos16 == pcnt2);
    (void) pcnt2;

#if BWAMEM_BATCHED_MATESW
    pwsw->getScores16(pairs + pos8, seqBufRef, seqBufQer, aln, pos16, nthreads, 1);
    pwsw->getScores8(pairs, seqBufRef, seqBufQer, aln, pos8, nthreads, 1);
#else
    fprintf(stderr, "Error: mem_sam_pe_batch reached without a batched kswv kernel\n");
    exit(EXIT_FAILURE);
#endif
}

// This function is equivalent to align2() for axv512
int mem_sam_pe_batch(const mem_opt_t *opt, mem_cache *mmc,
                     int64_t &pcnt, int64_t &pcnt8, kswr_t *aln,
                     int32_t maxRefLen, int32_t maxQerLen, int tid)
{
    uint8_t *seqBufRef = mmc->seqBufLeftRef[tid*CACHE_LINE];
    uint8_t *seqBufQer = mmc->seqBufLeftQer[tid*CACHE_LINE];

    SeqPair *seqPairArray = mmc->seqPairArrayLeft128[tid];

#if DEBUG    // orig function from bwa-mem -- for debugging purpose. Disabled by default.
    // uint64_t tim = __rdtsc();
    SeqPair sp;
    for (int i=0; i<pcnt; i++) {
        sp = seqPairArray[i];
        int xtra = sp.h0;
        uint8_t *qs = seqBufQer + sp.idq;
        uint8_t *rs = seqBufRef + sp.idr;
        aln[i] = ksw_align2(sp.len2, qs, sp.len1, rs, 5,
                            opt->mat, opt->o_del, opt->e_del,
                            opt->o_ins, opt->e_ins, xtra, 0);
    }
    // tprof[SAM2][0] += __rdtsc() - tim;

#else   // avx512, vectorized function

    for (int i=0; i<pcnt; i++) {
        kswr_t *r = &aln[i];
        r->tb = r->qb = -1;
    }

    int nthreads = 1; // no multi-threading here

    /* Task 5: under --meth (batched rescue enabled), partition the enqueued
     * mate-rescue pairs by bisulfite hypothesis (OT/OB) and run the kernel once
     * per group with the matching mat-aware kswv object. mat_ot frees C->T (top
     * strand), mat_ob frees G->A (bottom strand); make_kswv installs the rank-1
     * freed-cell override for each, so the batched score equals the scalar
     * ksw_align2 score with the asymmetric matrix. Non-meth and the env-OFF
     * escape hatch fall through to the single-object path below, byte-identical
     * to the pre-Task-5 code. */
    if (opt->meth_mode && meth_batched_rescue_enabled()) {
        // Scratch slice (Right128 is free here; sort_classify only used it as
        // a transient and the final layout lives in Left128). Sized
        // wsize_pair + MAX_LINE_LEN, which comfortably holds one group plus
        // the 16-bit shift headroom (a group is <= pcnt <= wsize_pair pairs).
        SeqPair *scratch = mmc->seqPairArrayRight128[tid];

        auto pwsw_ot = make_kswv(opt->o_del, opt->e_del, opt->o_ins, opt->e_ins,
                                 opt->a, -1*opt->b, nthreads,
                                 maxRefLen, maxQerLen, opt->mat_ot);
        auto pwsw_ob = make_kswv(opt->o_del, opt->e_del, opt->o_ins, opt->e_ins,
                                 opt->a, -1*opt->b, nthreads,
                                 maxRefLen, maxQerLen, opt->mat_ob);
        // The batched kernel expresses exactly the matrices mem_opt_fill_meth_mat
        // produces: GENOMIC (one freed cell) and COLLAPSED (the conversion cell
        // PLUS its mirror — a symmetric (i,j)/(j,i) pair). needsScalar() is
        // therefore always false here. Guard it at RUNTIME rather than with a
        // bare assert(): under NDEBUG assert is a no-op, which would let a future
        // unsupported matrix run the batched kernel and silently mis-score. Fail
        // loudly instead — this can only fire on a programming error in
        // mem_opt_fill_meth_mat (a matrix that is neither rank-1 nor a mirror pair).
        if (pwsw_ot->needsScalar() || pwsw_ob->needsScalar())
            err_fatal(__func__,
                      "meth mate-rescue matrix not expressible by the batched kernel "
                      "(neither rank-1 genomic nor a collapsed mirror pair)");

        // hyp == 1 -> OT (mat_ot / pwsw_ot); hyp == 0 -> OB (mat_ob / pwsw_ob).
        // -1 (non-meth) cannot occur here: every enqueued pair under meth_mode
        // was tagged with a real hypothesis in mem_matesw_batch_pre. Process
        // each group in the shared scratch buffer one at a time; aln[] is keyed
        // by regid and seqBuf* by idr/idq, so each pair is scored exactly once
        // regardless of grouping.
        int64_t scored = 0;
        for (int hyp = 1; hyp >= 0; --hyp) {
            int64_t n8 = 0, n16 = 0;
            // 8-bit pairs of this hyp first, then 16-bit, mirroring
            // sort_classify's layout so mem_sam_pe_batch_run's split is valid.
            for (int64_t i = 0; i < pcnt; ++i) {
                SeqPair *s = &seqPairArray[i];
                if (s->meth_hyp != hyp) continue;
                if (s->h0 & KSW_XBYTE) scratch[n8++] = *s;
            }
            int64_t group_pcnt8 = n8;
            for (int64_t i = 0; i < pcnt; ++i) {
                SeqPair *s = &seqPairArray[i];
                if (s->meth_hyp != hyp) continue;
                if (!(s->h0 & KSW_XBYTE)) scratch[group_pcnt8 + (n16++)] = *s;
            }
            int64_t group_pcnt = n8 + n16;
            scored += group_pcnt;
            if (group_pcnt == 0) continue;

            Ikswv *pwsw = (hyp == 1) ? pwsw_ot.get() : pwsw_ob.get();
            mem_sam_pe_batch_run(pwsw, scratch, seqBufRef, seqBufQer,
                                 aln, group_pcnt, group_pcnt8, nthreads);
        }
        // Release-safe invariant guard: the OT/OB partition must cover every
        // enqueued pair. A pair tagged with a hypothesis outside {0,1} would be
        // skipped by both groups, leaving its aln[] entry uninitialized for
        // mem_matesw_batch_post to read. The assert in mem_matesw_batch_pre
        // catches this in debug; fail loud under NDEBUG too rather than emit a
        // silently mis-scored rescue.
        if (scored != pcnt)
            err_fatal(__func__,
                      "batched meth rescue partition missed %ld of %ld pairs "
                      "(a pair carried a hypothesis outside {OT,OB})",
                      (long)(pcnt - scored), (long)pcnt);
        return 1;
    }

    auto pwsw = make_kswv(opt->o_del, opt->e_del, opt->o_ins, opt->e_ins,
                          opt->a, -1*opt->b, nthreads,
                          maxRefLen, maxQerLen);

    mem_sam_pe_batch_run(pwsw.get(), seqPairArray, seqBufRef, seqBufQer,
                         aln, pcnt, pcnt8, nthreads);

#endif

    return 1;
}

int mem_sam_pe_batch_post(const mem_opt_t *opt, const bntseq_t *bns,
                          const uint8_t *pac, const mem_pestat_t pes[4],
                          uint64_t id, bseq1_t s[2], mem_alnreg_v a[2],
                          kswr_t **myaln, mem_cache *mmc, 
                          int32_t &gcnt, int tid)
{
    extern int mem_mark_primary_se(const mem_opt_t *opt, int n, mem_alnreg_t *a, int64_t id);
    extern int mem_approx_mapq_se(const mem_opt_t *opt, const mem_alnreg_t *a);
    extern void mem_reg2sam(const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac,
                            bseq1_t *s, mem_alnreg_v *a, int extra_flag, const mem_aln_t *m);
    extern char **mem_gen_alt(const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac,
                              const mem_alnreg_v *a, int l_query, const char *query,
                              int **out_hn);
    #if MATE_SORT
    extern void sort_alnreg_re(int n, mem_alnreg_t* a);
    extern void sort_alnreg_score(int n, mem_alnreg_t* a);
    /* D3 (--meth, PR-6): dedup-only call (query = 0, no patch SW) — matrix
     * unused, opt->mat default correct; asymmetric scoring is in
     * mem_matesw_batch_post's scalar fallback. */
    #endif

    int32_t *gar = (int32_t*) mmc->seqPairArrayAux[tid];
    
    int n = 0, i, j, z[2], o, subo, n_sub, extra_flag = 1, n_pri[2], n_aa[2];
    kstring_t str;
    mem_aln_t h[2], g[2], aa[2][2];
    // int tid = omp_get_thread_num();
    
    str.l = str.m = 0; str.s = 0;
    memset(h, 0, sizeof(mem_aln_t) * 2);
    memset(g, 0, sizeof(mem_aln_t) * 2);
    n_aa[0] = n_aa[1] = 0;
    
    if (!(opt->flag & MEM_F_NO_RESCUE)) { // then perform SW for the best alignment
        mem_alnreg_v b[2];
        kv_init(b[0]); kv_init(b[1]);
        for (i = 0; i < 2; ++i)
            for (j = 0; j < a[i].n; ++j)
                if (a[i].a[j].score >= a[i].a[0].score  - opt->pen_unpaired)
                    kv_push(mem_alnreg_t, b[i], a[i].a[j]);
                
        for (int l=0; l<a[0].n; l++)
            a[0].a[l].flg = 0;
        for (int l=0; l<a[1].n; l++)
            a[1].a[l].flg = 0;

        #if MATE_SORT
        for (i = 0; i < 2; ++i) {
            sort_alnreg_re(a[!i].n, a[!i].a);
            int val = 0, swcount = 0;
            for (j = 0; j < b[i].n && j < opt->max_matesw; ++j) {
                /* D3 (--meth, PR-6/A1, B3): original mate bases + opposite-strand
                 * matrix of anchor b[i].a[j], honored on the scalar ksw_align2
                 * path through which mem_matesw_batch_post scores every meth
                 * rescue (the batched kswv enqueue is skipped under --meth). */
                const char  *ms_orig = opt->meth_mode ? s[!i].meth_orig_seq : NULL;
                const int8_t *rmat    = opt->meth_mode
                    ? mem_opt_meth_mat(opt, !b[i].a[j].meth_hypothesis) : NULL;
                val = mem_matesw_batch_post(opt, bns, pac, pes, &b[i].a[j],
                                                s[!i].l_seq, (uint8_t*)s[!i].seq,
                                                &a[!i], myaln, gcnt, gar, mmc,
                                                ms_orig, rmat);
                n += val;
                swcount += val;
                // ncnt++;
                gcnt += 4;
            }
            if (swcount > 0) {
                mem_alnreg_v* ma = &a[!i];
                ma->n = mem_sort_dedup_patch(opt, 0, 0, 0, ma->n, ma->a);
            }
            else {
                sort_alnreg_score(a[!i].n, a[!i].a);
            }
        }
        #else
        for (i = 0; i < 2; ++i) {
            for (j = 0; j < b[i].n && j < opt->max_matesw; ++j) {
                /* D3 (--meth, PR-6, B3): see MATE_SORT branch above. */
                const char  *ms_orig = opt->meth_mode ? s[!i].meth_orig_seq : NULL;
                const int8_t *rmat    = opt->meth_mode
                    ? mem_opt_meth_mat(opt, !b[i].a[j].meth_hypothesis) : NULL;
                int val = mem_matesw_batch_post(opt, bns, pac, pes, &b[i].a[j],
                                                s[!i].l_seq, (uint8_t*)s[!i].seq,
                                                &a[!i], myaln, gcnt, gar, mmc,
                                                ms_orig, rmat);
                n += val;
                // ncnt++;
                gcnt += 4;
            }
        }
        #endif
        free(b[0].a); free(b[1].a);
    }

    n_pri[0] = mem_mark_primary_se(opt, a[0].n, a[0].a, id<<1|0);
    n_pri[1] = mem_mark_primary_se(opt, a[1].n, a[1].a, id<<1|1);

    #if V17
    if (opt->flag & MEM_F_PRIMARY5) {
        mem_reorder_primary5(opt->T, &a[0]);
        mem_reorder_primary5(opt->T, &a[1]);
    }
    #endif
    
    if (opt->flag&MEM_F_NOPAIRING) goto no_pairing;

    // pairing single-end hits
    if (n_pri[0] && n_pri[1] && (o = mem_pair(opt, bns, pac, pes, s, a, id, &subo, &n_sub, z, n_pri)) > 0)
    {
        int is_multi[2], q_pe, score_un, q_se[2];
        char **XA[2];
        // check if an end has multiple hits even after mate-SW
        for (i = 0; i < 2; ++i) {
            for (j = 1; j < n_pri[i]; ++j)
                if (a[i].a[j].secondary < 0 && a[i].a[j].score >= opt->T) break;
            is_multi[i] = j < n_pri[i]? 1 : 0;
        }
        if (is_multi[0] || is_multi[1]) goto no_pairing; // TODO: in rare cases, the true hit may be long but with low score
        // compute mapQ for the best SE hit
        score_un = a[0].a[0].score + a[1].a[0].score - opt->pen_unpaired;
        //q_pe = o && subo < o? (int)(MEM_MAPQ_COEF * (1. - (double)subo / o) * log(a[0].a[z[0]].seedcov + a[1].a[z[1]].seedcov) + .499) : 0;
        subo = subo > score_un? subo : score_un;
        q_pe = raw_mapq(o - subo, opt->a);

        if (n_sub > 0) q_pe -= (int)(4.343 * log(n_sub+1) + .499);
        if (q_pe < 0) q_pe = 0;
        if (q_pe > 60) q_pe = 60;

        q_pe = (int)(q_pe * (1. - .5 * (a[0].a[0].frac_rep + a[1].a[0].frac_rep)) + .499);

        // the following assumes no split hits
        if (o > score_un) { // paired alignment is preferred
            mem_alnreg_t *c[2];
            c[0] = &a[0].a[z[0]]; c[1] = &a[1].a[z[1]];
            for (i = 0; i < 2; ++i) {
                if (c[i]->secondary >= 0)
                    c[i]->sub = a[i].a[c[i]->secondary].score, c[i]->secondary = -2;
                q_se[i] = mem_approx_mapq_se(opt, c[i]);
            }

            q_se[0] = q_se[0] > q_pe? q_se[0] : q_pe < q_se[0] + 40? q_pe : q_se[0] + 40;
            q_se[1] = q_se[1] > q_pe? q_se[1] : q_pe < q_se[1] + 40? q_pe : q_se[1] + 40;
            extra_flag |= 2;

            // cap at the tandem repeat score
            q_se[0] = q_se[0] < raw_mapq(c[0]->score - c[0]->csub, opt->a)? q_se[0] : raw_mapq(c[0]->score - c[0]->csub, opt->a);
            q_se[1] = q_se[1] < raw_mapq(c[1]->score - c[1]->csub, opt->a)? q_se[1] : raw_mapq(c[1]->score - c[1]->csub, opt->a);

        } else { // the unpaired alignment is preferred
            z[0] = z[1] = 0;
            q_se[0] = mem_approx_mapq_se(opt, &a[0].a[0]);
            q_se[1] = mem_approx_mapq_se(opt, &a[1].a[0]);

        }
        for (i = 0; i < 2; ++i) {
            int k = a[i].a[z[i]].secondary_all;
            if (k >= 0 && k < n_pri[i]) { // switch secondary and primary if both of them are non-ALT
                assert(a[i].a[k].secondary_all < 0);
                for (j = 0; j < a[i].n; ++j)
                    if (a[i].a[j].secondary_all == k || j == k)
                        a[i].a[j].secondary_all = z[i];
                a[i].a[z[i]].secondary_all = -1;
            }
        }
        int *HN[2] = { 0, 0 };
        if (!(opt->flag & MEM_F_ALL)) {
            for (i = 0; i < 2; ++i)
                XA[i] = mem_gen_alt(opt, bns, pac, &a[i], s[i].l_seq, s[i].seq, &HN[i]);
        } else XA[0] = XA[1] = 0;
        // write SAM
        for (i = 0; i < 2; ++i) {
            h[i] = mem_reg2aln(opt, bns, pac, s[i].l_seq, s[i].seq, &a[i].a[z[i]], s[i].meth_orig_seq);
            h[i].mapq = q_se[i];

            h[i].flag |= 0x40<<i | extra_flag;
            h[i].XA = XA[i]? XA[i][z[i]] : 0;
            h[i].HN = HN[i]? HN[i][z[i]] : -1;
            aa[i][n_aa[i]++] = h[i];
            if (n_pri[i] < a[i].n) { // the read has ALT hits
                mem_alnreg_t *p = &a[i].a[n_pri[i]];
                if (p->score < opt->T || p->secondary >= 0 || !p->is_alt) continue;
                g[i] = mem_reg2aln(opt, bns, pac, s[i].l_seq, s[i].seq, p, s[i].meth_orig_seq);
                g[i].flag |= 0x800 | 0x40<<i | extra_flag;
                g[i].XA = XA[i]? XA[i][n_pri[i]] : 0;
                g[i].HN = HN[i]? HN[i][n_pri[i]] : -1;
                if (opt->supp_rep_hard_cap > 0 && p->chain_n_hits >= opt->supp_rep_hard_cap)
                    g[i].mapq = 0; // fg-labs: force repetitive-supp MAPQ to 0
                aa[i][n_aa[i]++] = g[i];
            }
        }
        for (i = 0; i < n_aa[0]; ++i)
            mem_aln2sam(opt, bns, &str, &s[0], n_aa[0], aa[0], i, &h[1]);

        if (opt->bam_mode) {
            /* BAM path (meth or generic): mem_aln2sam short-circuited
             * into s->bams, leaving str untouched. Skip the str.s dance. */
            s[0].sam = NULL;
            str.l = 0;
            for (i = 0; i < n_aa[1]; ++i)
                mem_aln2sam(opt, bns, &str, &s[1], n_aa[1], aa[1], i, &h[0]);
            s[1].sam = NULL;
            free(str.s); str.s = NULL; str.m = 0;
        } else {
            assert(str.s != 0);
            s[0].sam = strdup(str.s); str.l = 0;
            for (i = 0; i < n_aa[1]; ++i)
                mem_aln2sam(opt, bns, &str, &s[1], n_aa[1], aa[1], i, &h[0]);
            s[1].sam = str.s;
        }
        if (strcmp(s[0].name, s[1].name) != 0) err_fatal(__func__, "paired reads have different names: \"%s\", \"%s\"\n", s[0].name, s[1].name);
        // free
        for (i = 0; i < 2; ++i) {
            free(h[i].cigar); free(g[i].cigar);
            free(HN[i]);
            if (XA[i] == 0) continue;
            for (j = 0; j < a[i].n; ++j) free(XA[i][j]);
            free(XA[i]);
        }
    } else goto no_pairing;
    return n;

no_pairing:
    int which[2] = { -1, -1 };
    for (i = 0; i < 2; ++i) {
        if (a[i].n) {
            if (a[i].a[0].score >= opt->T) which[i] = 0;
            else if (n_pri[i] < a[i].n && a[i].a[n_pri[i]].score >= opt->T)
                which[i] = n_pri[i];
        }
        if (which[i] >= 0) h[i] = mem_reg2aln(opt, bns, pac, s[i].l_seq, s[i].seq, &a[i].a[which[i]], s[i].meth_orig_seq);
        else h[i] = mem_reg2aln(opt, bns, pac, s[i].l_seq, s[i].seq, 0);
    }
    // Proper-pair flag must be computed from the same alignments that were just
    // emitted via mem_reg2aln — i.e. a[i].a[which[i]]. Using a[i].a[0] here is
    // wrong when which[i] == n_pri[i] (below-T primary + above-T ALT case).
    if (!(opt->flag & MEM_F_NOPAIRING) && which[0] >= 0 && which[1] >= 0 &&
        h[0].rid == h[1].rid && h[0].rid >= 0) {
        int64_t dist;
        int d;
        d = mem_infer_dir(bns->l_pac, a[0].a[which[0]].rb, a[1].a[which[1]].rb, &dist);
        if (!pes[d].failed && dist >= pes[d].low && dist <= pes[d].high) extra_flag |= 2;
    }
    mem_reg2sam(opt, bns, pac, &s[0], &a[0], 0x41|extra_flag, &h[1]);
    mem_reg2sam(opt, bns, pac, &s[1], &a[1], 0x81|extra_flag, &h[0]);
    if (strcmp(s[0].name, s[1].name) != 0)
        err_fatal(__func__, "paired reads have different names: \"%s\", \"%s\"\n",
                  s[0].name, s[1].name);

    free(h[0].cigar); free(h[1].cigar);
    return n;
}


int mem_matesw_batch_pre(const mem_opt_t *opt, const bntseq_t *bns,
                         const uint8_t *pac, const mem_pestat_t pes[4],
                         const mem_alnreg_t *a, int l_ms, const uint8_t *ms,
                         mem_alnreg_v *ma, mem_cache *mmc, int pcnt, int32_t gcnt,
                         int32_t &maxRefLen, int32_t &maxQerLen, int32_t tid)
{
    /* D3 (--meth, A1/PR-6 → issue 173 Task 5): _pre stages the batched-SIMD
     * mate-rescue query for the kswv getScores8/16 kernel. As of the mat-aware
     * kernels (Tasks 1-4) the batched kernel CAN express the bisulfite OT/OB
     * freed-cell matrices, so under --meth we now stage meth pairs here just
     * like non-meth pairs and tag each with its bisulfite hypothesis
     * (sp.meth_hyp, set near the enqueue below). mem_sam_pe_batch then
     * partitions the enqueued pairs by OT/OB and scores each group with the
     * matching mat-aware kswv object, byte-identical to the per-hypothesis
     * scalar ksw_align2. The legacy scalar path (gar[gcnt+r] = -1 →
     * mem_matesw_batch_post re-runs ksw_align2) is retained as the escape hatch
     * (BWAMEM3_METH_BATCHED_RESCUE=0). The OT/OB matrices mem_opt_fill_meth_mat
     * produces are always expressible (rank-1 genomic or a collapsed mirror
     * pair), so the kswv objects never report needsScalar() here; if a future
     * matrix or tier ever did, mem_sam_pe_batch fails loud via err_fatal rather
     * than silently mis-scoring (see the needsScalar() guard there). The dedup
     * below passes query = 0 (no patch SW), so its mat default (opt->mat) is
     * unused. */

    uint8_t *seqBufRef = mmc->seqBufLeftRef[tid*CACHE_LINE];
    uint8_t *seqBufQer = mmc->seqBufLeftQer[tid*CACHE_LINE];
    SeqPair *seqPairArray = mmc->seqPairArrayLeft128[tid];
    int32_t *gar = (int32_t*) (mmc->seqPairArrayAux[tid]);

    int64_t *wsize_pair = &(mmc->wsize[tid]);
    int64_t *wsize_buf_ref = &(mmc->wsize_buf_ref[tid*CACHE_LINE]);
    int64_t *wsize_buf_qer = &(mmc->wsize_buf_qer[tid*CACHE_LINE]);
    
    int64_t l_pac = bns->l_pac;
    int i, r, skip[4], rid = -1;
    for (r = 0; r < 4; ++r)
        skip[r] = pes[r].failed? 1 : 0;

    for (i = 0; i < ma->n; ++i) { // check which orinentation has been found        
        int64_t dist;
        r = mem_infer_dir(l_pac, a->rb, ma->a[i].rb, &dist);
        if (dist >= pes[r].low && dist <= pes[r].high) {
            skip[r] = 1;
        }
    }

    
    if (skip[0] + skip[1] + skip[2] + skip[3] == 4) //return pcnt; // consistent pair exist; no need to perform SW
    {
        gar[gcnt + 3] = gar[gcnt + 2] = gar[gcnt + 1] = gar[gcnt + 0] = -1;
        return pcnt;
    }
        
    for (r = 0; r < 4; ++r)
    {
        int is_rev, is_larger;
        uint8_t *seq, *rev = 0, *ref = 0;
        int64_t rb, re;
        if (skip[r]) {
            gar[gcnt + r] = -1;
            continue;
        }
        is_rev = (r>>1 != (r&1)); // whether to reverse complement the mate
        is_larger = !(r>>1); // whether the mate has larger coordinate

        if (is_rev) {
            rev = (uint8_t*) malloc(l_ms); // this is the reverse complement of $ms
            assert(rev != NULL);
            for (i = 0; i < l_ms; ++i) rev[l_ms - 1 - i] = ms[i] < 4? 3 - ms[i] : 4;
            seq = rev;
        } else seq = (uint8_t*)ms;

        if (!is_rev) {
            rb = is_larger? a->rb + pes[r].low : a->rb - pes[r].high;
            re = (is_larger? a->rb + pes[r].high: a->rb - pes[r].low) + l_ms; // if on the same strand, end position should be larger to make room for the seq length
        }
        else {
            rb = (is_larger? a->rb + pes[r].low : a->rb - pes[r].high) - l_ms; // similarly on opposite strands
            re = is_larger? a->rb + pes[r].high: a->rb - pes[r].low;
        }

        if (rb < 0) rb = 0;
        if (re > l_pac<<1) re = l_pac<<1;
        // bns_fetch_seq_v2 is zero-copy: it returns a pointer into ref_string
        // (the unpacked .0123 reference). The legacy bns_fetch_seq malloc'd a
        // fresh buffer, which the trailing free(ref) below released; v2 has
        // no allocation, so we drop that free. The seqPairArrayAux scratch
        // is unused by v2 but kept for signature parity with the existing v2
        // caller in bwamem.cpp::mem_chain2aln_across_reads_V2.
        if (rb < re) ref = bns_fetch_seq_v2(bns, pac, &rb, (rb+re)>>1, &re, &rid,
                                            mmc->ref_string,
                                            (uint8_t*) mmc->seqPairArrayAux[tid]);

        if (a->rid == rid && re - rb >= opt->min_seed_len) { // no funny things happening
            /* D3 (--meth, Task 5): meth mate-rescue is now enqueued for the
             * batched kswv kernel just like a non-meth pair. The kernels are
             * mat-aware (Task 2/4): make_kswv(..., mat_ot/mat_ob) installs the
             * rank-1 freed-cell override for the bisulfite OT/OB matrices, so
             * the batched score matches the scalar ksw_align2 score with the
             * asymmetric matrix. The per-pair hypothesis tag (sp.meth_hyp,
             * below) drives the OT/OB partition in mem_sam_pe_batch, which runs
             * getScores8/16 once per hypothesis group with the matching object.
             * The rescued mate uses the OPPOSITE-strand matrix of the anchor
             * `a` (mem_opt_meth_mat(opt, !a->meth_hypothesis)), exactly as
             * mem_matesw_batch_post's scalar fallback does. The legacy scalar
             * path remains as a safety fallback for any pair whose object
             * reports needsScalar() (rank-1 meth never does) and is reachable
             * via BWAMEM3_METH_BATCHED_RESCUE=0 (escape hatch). */
            if (opt->meth_mode && !meth_batched_rescue_enabled()) {
                // Escape hatch (env=0): keep the legacy scalar rescue. Leave
                // gar = -1 so mem_matesw_batch_post re-runs this orientation
                // through ksw_align2 with the asymmetric matrix.
                gar[gcnt + r] = -1;
                if (rev) free(rev);
                continue;
            }
            //kswr_t aln;
            //mem_alnreg_t b;
            int xtra = KSW_XSUBO | KSW_XSTART | (l_ms * opt->a < 250? KSW_XBYTE : 0) | (opt->min_seed_len * opt->a);
            int qerOffset = 0, refOffset = 0;
            if (pcnt != 0)
            {
                SeqPair sp;
                sp = seqPairArray[pcnt - 1];
                refOffset = sp.idr + sp.len1;
                qerOffset = sp.idq + sp.len2;
            }

            SeqPair sp;
            sp.h0 = xtra;
            // assert(pcnt < BATCH_SIZE * SEEDS_PER_READ);
            assert(pcnt < *wsize_pair);
            
            sp.idq = qerOffset;
            sp.idr = refOffset;
            sp.len1 = re - rb;
            sp.len2 = l_ms;
            sp.id = sp.score = sp.seqid = sp.gtle = sp.tle = sp.qle = sp.max_off = sp.gscore = -1; // not needed, remove while code cleaning
            
            assert(sp.len1 >= 0 && sp.len2 >= 0);
            if (refOffset + sp.len1 >= *wsize_buf_ref)
            {
                fprintf(stderr, "[0000][%0.4d] Re-allocating (doubling) seqBufRefs in %s\n",
                        tid, __func__);
                int64_t tmp = *wsize_buf_ref;
                *wsize_buf_ref *= 2;

                uint8_t *seqBufRef_ = (uint8_t*)
                    _mm_realloc(seqBufRef, tmp, *wsize_buf_ref, sizeof(uint8_t)); 
                mmc->seqBufLeftRef[tid*CACHE_LINE] = seqBufRef = seqBufRef_;

                seqBufRef_ = (uint8_t*)
                    _mm_realloc(mmc->seqBufRightRef[tid*CACHE_LINE], tmp,
                                *wsize_buf_ref, sizeof(uint8_t)); 
                mmc->seqBufRightRef[tid*CACHE_LINE] = seqBufRef_;               
            }
            
            if (qerOffset + sp.len2 >= *wsize_buf_qer)
            {
                fprintf(stderr, "[0000][%0.4d] Re-allocating (doubling) seqBufQers in %s\n",
                        tid, __func__);
                int64_t tmp = *wsize_buf_qer;
                *wsize_buf_qer *= 2;

                uint8_t *seqBufQer_ = (uint8_t*)
                    _mm_realloc(seqBufQer, tmp, *wsize_buf_qer, sizeof(uint8_t)); 
                mmc->seqBufLeftQer[tid*CACHE_LINE] = seqBufQer = seqBufQer_;

                seqBufQer_ = (uint8_t*)
                    _mm_realloc(mmc->seqBufRightQer[tid*CACHE_LINE], tmp,
                                *wsize_buf_qer, sizeof(uint8_t)); 
                mmc->seqBufRightQer[tid*CACHE_LINE] = seqBufQer_;               
            }
            
            if (pcnt >= *wsize_pair)
            {
                fprintf(stderr, "[0000][%0.4d] Re-allocating seqPairs in %s\n", tid, __func__);
                *wsize_pair += 1024;
                mmc->seqPairArrayAux[tid] = (SeqPair *) realloc(mmc->seqPairArrayAux[tid],
                                                    (*wsize_pair + MAX_LINE_LEN)
                                                    * sizeof(SeqPair));
                mmc->seqPairArrayLeft128[tid] = (SeqPair *) realloc(mmc->seqPairArrayLeft128[tid],
                                                    (*wsize_pair + MAX_LINE_LEN)
                                                    * sizeof(SeqPair));
                mmc->seqPairArrayRight128[tid] = (SeqPair *) realloc(mmc->seqPairArrayRight128[tid],
                                                    (*wsize_pair + MAX_LINE_LEN)
                                                    * sizeof(SeqPair));
                seqPairArray = mmc->seqPairArrayLeft128[tid];
                gar = (int32_t*) (mmc->seqPairArrayAux[tid]);               
            }

            if (maxRefLen < sp.len1) maxRefLen = sp.len1;
            if (maxQerLen < sp.len2) maxQerLen = sp.len2;
            
            uint8_t *qs = seqBufQer + sp.idq;
            uint8_t *rs = seqBufRef + sp.idr;
            for (int l=0; l<sp.len1; l++) rs[l] = ref[l];
            for (int l=0; l<sp.len2; l++) qs[l] = seq[l];

            /* Task 5: tag the enqueued meth pair with the rescued mate's
             * bisulfite hypothesis = OPPOSITE strand of the anchor `a`
             * (mem_opt_meth_mat(opt, !a->meth_hypothesis) — the exact value
             * mem_matesw_batch_post's scalar fallback already uses to pick
             * rmat). is_rev was already baked into a->meth_hypothesis upstream,
             * so we reuse it verbatim and do NOT re-derive strand logic here.
             *
             * Under --meth every chain that survives to mate rescue carries a
             * real hypothesis (a->meth_hypothesis ∈ {0,1}): meth_seed_to_orig
             * either remaps a seed with a concrete hypothesis or drops it, so a
             * negative hypothesis never reaches an enqueued anchor. The tag is
             * therefore always in {0,1} under --meth, which the OT/OB partition
             * in mem_sam_pe_batch requires. Non-meth pairs keep the SeqPair
             * default (-1), which the kernels ignore. The assert documents the
             * invariant; mem_sam_pe_batch additionally guards it at runtime so a
             * future violation fails loud instead of silently dropping a pair. */
            sp.meth_hyp = opt->meth_mode ? (int8_t)!a->meth_hypothesis : (int8_t)-1;
            assert(!opt->meth_mode || sp.meth_hyp == 0 || sp.meth_hyp == 1);

            gar[gcnt + r] = pcnt;
            sp.regid = pcnt;
            seqPairArray[pcnt++] = sp;
        }
        if (rev) free(rev);
        // ref aliases ref_string (see bns_fetch_seq_v2 above); no free.
    }
    return pcnt;
}

int mem_matesw_batch_post(const mem_opt_t *opt, const bntseq_t *bns,
                          const uint8_t *pac, const mem_pestat_t pes[4],
                          const mem_alnreg_t *a, int l_ms, const uint8_t *ms,
                          mem_alnreg_v *ma, kswr_t **myaln, int32_t gcnt,
                          int32_t *gar, mem_cache *mmc, const char *ms_orig,
                          const int8_t *mat)
{
    extern int mem_sort_dedup_patch_rev(const mem_opt_t *opt, const bntseq_t *bns,
                                        const uint8_t *pac, uint8_t *query, int n,
                                        mem_alnreg_t *a);
    /* The mate-rescue dedup at the bottom passes bns/pac/query = 0 (dedup-only,
     * no patch SW), so its matrix is unused; opt->mat default is correct there. */
    /* D3 (--meth, A1/PR-6, B3): meth mate rescue is now scored by the BATCHED
     * kswv kernel (via the OT/OB hypothesis partition in mem_sam_pe_batch) by
     * default — _pre enqueues each pair and tags it with sp.meth_hyp ∈ {0,1}.
     * The scalar ksw_align2 path below (index == -1) is the ESCAPE HATCH: it
     * is reached only when BWAMEM3_METH_BATCHED_RESCUE=0 forces the legacy
     * scalar path (gar[gcnt+r] == -1 from _pre), or for any matrix where the
     * kswv object reports needsScalar() (a non-rank-1 asymmetric matrix; the
     * current OT/OB matrices are rank-1 and never trigger this). The path is
     * retained as a safety net, not the primary route. */
    const int8_t *sw_mat = mat ? mat : opt->mat;
    #if MATE_SORT    
    extern int mem_dedup_patch(const mem_opt_t *opt, const bntseq_t *bns,
                               const uint8_t *pac, uint8_t *query, int n, mem_alnreg_t *a);
    extern void sort_alnreg_re(int n, mem_alnreg_t* a);
    extern void sort_alnreg_score(int n, mem_alnreg_t* a);
    #endif
    
    int64_t l_pac = bns->l_pac;
    int i, r, skip[4], n = 0, rid = -1;

    /* D3 (--meth, PR-6): for the scalar ksw_align2 fallback (index == -1) point
     * the mate query at the ORIGINAL bases (ASCII meth_orig_seq, same orientation
     * as `ms`/seq). The per-orientation RC below then reverse-complements them
     * exactly as the projected mate. The batched SIMD scores are unaffected (they
     * were computed in mem_sam_pe_batch from _pre's projected query). */
    uint8_t *ms2 = NULL;
    if (ms_orig != NULL) {
        ms2 = (uint8_t*) malloc(l_ms);
        assert(ms2 != NULL);
        for (int k = 0; k < l_ms; ++k) {
            unsigned char c = (unsigned char) ms_orig[k];
            ms2[k] = (c < 4) ? c : nst_nt4_table[c];
        }
        ms = ms2;
    }

    for (r = 0; r < 4; ++r) {
        skip[r] = pes[r].failed? 1 : 0;
    }

    for (i = 0; i < ma->n; ++i) { // check which orinentation has been found
        int64_t dist;
        r = mem_infer_dir(l_pac, a->rb, ma->a[i].rb, &dist);
        if (dist >= pes[r].low && dist <= pes[r].high)
            skip[r] = 1;
    }

    
    if (skip[0] + skip[1] + skip[2] + skip[3] == 4) {
        return 0; // consistent pair exist; no need to perform SW
    }

    for (r = 0; r < 4; ++r) {
        int is_rev, is_larger;
        uint8_t *seq, *rev = 0, *ref = 0;
        int64_t rb, re;
        if (skip[r]) {
                continue;
        }
        is_rev = (r>>1 != (r&1)); // whether to reverse complement the mate
        is_larger = !(r>>1); // whether the mate has larger coordinate
        if (is_rev) {
            rev = (uint8_t*) malloc(l_ms); // this is the reverse complement of $ms
            assert(rev != NULL);
            for (i = 0; i < l_ms; ++i) rev[l_ms - 1 - i] = ms[i] < 4? 3 - ms[i] : 4;
            seq = rev;
        } else seq = (uint8_t*)ms;
        if (!is_rev) {
            rb = is_larger? a->rb + pes[r].low : a->rb - pes[r].high;
            re = (is_larger? a->rb + pes[r].high: a->rb - pes[r].low) + l_ms; // if on the same strand, end position should be larger to make room for the seq length
        } else {
            rb = (is_larger? a->rb + pes[r].low : a->rb - pes[r].high) - l_ms; // similarly on opposite strands
            re = is_larger? a->rb + pes[r].high: a->rb - pes[r].low;
        }
        if (rb < 0) rb = 0;
        if (re > l_pac<<1) re = l_pac<<1;
        // Zero-copy ref slice via bns_fetch_seq_v2 (see mem_matesw_batch_pre
        // for rationale). The scratch arg is unused by v2; pass NULL since
        // mem_matesw_batch_post has no tid in scope to index seqPairArrayAux.
        if (rb < re) ref = bns_fetch_seq_v2(bns, pac, &rb, (rb+re)>>1, &re, &rid,
                                            mmc->ref_string, NULL);

        if (a->rid == rid && re - rb >= opt->min_seed_len) { // no funny things happening
            kswr_t aln;
            mem_alnreg_t b;
            int tmp, xtra = KSW_XSUBO | KSW_XSTART | (l_ms * opt->a < 250? KSW_XBYTE : 0) | (opt->min_seed_len * opt->a);

            //aln = **myaln;
            //(*myaln)++;
            int index = gar[gcnt + r];
            if (index == -1) {
                // fprintf(stderr, "Re-routing: Encountered -ve index for "
                // "gcnt: %d, look into pre.\n", gcnt + r);
                assert(ref != 0);
                // ksw_align2 reverses its target argument in place via
                // revseq (see ksw.cpp:375,381). When mmc->ref_string is
                // shm-backed (PROT_READ mmap of /dev/shm/bwaidx-*), that
                // write SIGSEGVs. Copy the slice into a writable scratch
                // buffer before handing it to ksw_align2.
                int64_t ref_len = re - rb;
                uint8_t *ref_rw = (uint8_t*) malloc((size_t)ref_len);
                assert(ref_rw != NULL);
                memcpy(ref_rw, ref, (size_t)ref_len);
                aln = ksw_align2(l_ms, seq, ref_len, ref_rw, 5,
                                 sw_mat, opt->o_del, opt->e_del,
                                 opt->o_ins, opt->e_ins, xtra, 0);
                free(ref_rw);
            }
            else
                aln = *(*myaln + index);

            memset(&b, 0, sizeof(mem_alnreg_t));
            if (aln.score >= opt->min_seed_len && aln.qb >= 0) { // something goes wrong if aln.qb < 0
                b.rid = a->rid;
                b.is_alt = a->is_alt;
                /* D3 (--meth, PR-6, B3): rescued mate hypothesis = OPPOSITE strand
                 * of the anchor (directional: rescuing an OT anchor's mate uses OB
                 * and vice versa); set from !a->meth_hypothesis so the output
                 * layer (XG/XM) sources the right strand. -1 anchor stays -1.
                 * (Under --meth this score came from the scalar ksw_align2 path
                 * with the asymmetric matrix — see the top of this function — so
                 * it is fully γ-correct, not symmetric/projected.)
                 * Coordinates are already ORIGINAL (l_pac is the original l_pac via
                 * the original bns), so the 6a coordinate fix holds. */
                b.meth_hypothesis = (a->meth_hypothesis < 0) ? -1
                                                             : !a->meth_hypothesis;
                b.qb = is_rev? l_ms - (aln.qe + 1) : aln.qb;
                b.qe = is_rev? l_ms - aln.qb : aln.qe + 1;
                b.rb = is_rev? (l_pac<<1) - (rb + aln.te + 1) : rb + aln.tb;
                b.re = is_rev? (l_pac<<1) - (rb + aln.tb) : rb + aln.te + 1;
                b.score = aln.score;
                b.csub = aln.score2;
                b.secondary = -1;
                b.seedcov = (b.re - b.rb < b.qe - b.qb? b.re - b.rb : b.qe - b.qb) >> 1;
                b.chain_n_hits = 1; // mate-rescue has no SMEM evidence; treat as unique anchor

                kv_push(mem_alnreg_t, *ma, b); // make room for a new element

                #if !MATE_SORT

                // move b s.t. ma is sorted
                for (i = 0; i < ma->n - 1; ++i) // find the insertion point
                    if (ma->a[i].score < b.score) break;
                tmp = i;
                for (i = ma->n - 1; i > tmp; --i) ma->a[i] = ma->a[i-1];
                ma->a[i] = b;

                #else
                int resort = 0;
                // move b s.t. ma is sorted
                for (i = 0; i < ma->n - 1; ++i) { // find the insertion point
                    if (ma->a[i].re == b.re) {
                        resort = 1;
                        break;
                    }
                    if (ma->a[i].re > b.re) {
                        break;
                    }
                }
                if (resort) {
                    // Don't know where to put this alignment. So let the scores decide
                    sort_alnreg_score(ma->n - 1, ma->a);
                    for (i = 0; i < ma->n - 1; ++i) { // find the insertion point
                        if (ma->a[i].score < b.score) {
                            break;
                        }
                    }
                    tmp = i;
                    for (i = ma->n - 1; i > tmp; --i) ma->a[i] = ma->a[i-1];
                    ma->a[i] = b;
                    // Now we can sort based on end position
                    sort_alnreg_re(ma->n, ma->a);
                }
                else {
                    tmp = i;
                    for (i = ma->n - 1; i > tmp; --i) ma->a[i] = ma->a[i-1];
                    ma->a[i] = b;
                }
                #endif
                tprof[PE26][0] ++;
            }
            ++n;
        }
        #if !MATE_SORT
        if (n) ma->n = mem_sort_dedup_patch(opt, 0, 0, 0, ma->n, ma->a);
        #else
        if (n) ma->n = mem_dedup_patch(opt, 0, 0, 0, ma->n, ma->a);
        #endif

        if (rev) free(rev);
        // ref aliases ref_string (see bns_fetch_seq_v2 above); no free.
    }
    if (ms2) free(ms2); // D3 (--meth): original-mate 2-bit scratch
    return n;
}

