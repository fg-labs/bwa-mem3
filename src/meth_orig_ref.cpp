/* SPDX-License-Identifier: MIT */

#include "meth_orig_ref.h"

#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>

struct meth_orig_ref_s {
    /* Borrowed (not owned) — these point at the bwaidx_t already loaded
     * by bwa_idx_load_from_disk for alignment. */
    const bntseq_t *bns;
    const uint8_t  *pac;

    /* Per real-chrom: bns rid for f and r contig at this real chrom, plus
     * the chrom length on the un-doubled forward axis. n_chroms ==
     * cmap->n_output; arrays indexed by real_tid. */
    int      n_chroms;
    int     *f_rid;       /* length n_chroms; -1 if no 'f'-contig at real_tid */
    int     *r_rid;       /* length n_chroms; -1 if no 'r'-contig at real_tid */
    int64_t *chrom_len;   /* length n_chroms; bns->anns[f_rid].len */
};

meth_orig_ref_t *g_meth_orig_ref = NULL;

namespace {

/* (f, r) -> original 2-bit base. Returns 0xff for an illegal pair (treated
 * by caller as N — bwa's random-fill of ambi positions can produce one). */
inline uint8_t recover_pair(uint8_t f, uint8_t r) {
    /* 2-bit: A=0, C=1, G=2, T=3.
     *   (f=A, r=A) -> A    forward original A
     *   (f=T, r=C) -> C    C in f got C->T; r left alone
     *   (f=G, r=A) -> G    f left alone; G in r got G->A
     *   (f=T, r=T) -> T    forward original T
     */
    if (f == 0 && r == 0) return 0;
    if (f == 3 && r == 1) return 1;
    if (f == 2 && r == 0) return 2;
    if (f == 3 && r == 3) return 3;
    return 0xff;
}

}  // namespace

extern "C"
meth_orig_ref_t *meth_orig_ref_load(const bntseq_t *bns, const uint8_t *pac,
                                    const meth_chrom_map_t *cmap)
{
    if (bns == NULL || pac == NULL || cmap == NULL) return NULL;
    meth_orig_ref_t *o = (meth_orig_ref_t *)calloc(1, sizeof(*o));
    if (o == NULL) return NULL;
    o->bns = bns;
    o->pac = pac;
    o->n_chroms = cmap->n_output;
    if (o->n_chroms <= 0) return o;

    o->f_rid     = (int *)    malloc((size_t)o->n_chroms * sizeof(int));
    o->r_rid     = (int *)    malloc((size_t)o->n_chroms * sizeof(int));
    o->chrom_len = (int64_t *)malloc((size_t)o->n_chroms * sizeof(int64_t));
    if (o->f_rid == NULL || o->r_rid == NULL || o->chrom_len == NULL) {
        meth_orig_ref_free(o); return NULL;
    }
    for (int t = 0; t < o->n_chroms; ++t) {
        o->f_rid[t]     = -1;
        o->r_rid[t]     = -1;
        o->chrom_len[t] = cmap->output_lens[t];
    }
    /* Walk cmap to fill in f_rid and r_rid per real_tid. */
    for (int rid = 0; rid < cmap->n_internal; ++rid) {
        int t = cmap->out_tid[rid];
        if (t < 0 || t >= o->n_chroms) continue;
        if (cmap->direction[rid] == 'f' && o->f_rid[t] < 0) o->f_rid[t] = rid;
        if (cmap->direction[rid] == 'r' && o->r_rid[t] < 0) o->r_rid[t] = rid;
    }
    /* Defensive: every real chrom should have both sibling rids. If not
     * (malformed cmap or single-direction index), fail load and let the
     * caller surface a clean error. */
    for (int t = 0; t < o->n_chroms; ++t) {
        if (o->f_rid[t] < 0 || o->r_rid[t] < 0) {
            fprintf(stderr,
                    "ERROR: meth: real chrom %d (\"%s\") missing %s sibling contig\n",
                    t, cmap->output_names[t],
                    o->f_rid[t] < 0 ? "f-prefixed" : "r-prefixed");
            meth_orig_ref_free(o);
            return NULL;
        }
    }
    return o;
}

namespace {

/* bns_iter_ambi callback context: write 'N' over [st, en) of dst, where
 * st/en are global-pac offsets that we translate back to dst[] indices
 * via pac_off (the f-contig's offset for this slice) and st_dst (the
 * caller's chrom-relative start). */
struct n_mask_ctx {
    uint8_t *dst;
    int64_t  pac_off;
    int64_t  st_dst;
};

int n_mask_visit(int64_t amb_st, int64_t amb_en, void *vctx)
{
    n_mask_ctx *c = (n_mask_ctx *)vctx;
    int64_t lo = amb_st - c->pac_off - c->st_dst;
    int64_t hi = amb_en - c->pac_off - c->st_dst;
    for (int64_t j = lo; j < hi; ++j) c->dst[j] = 'N';
    return 0;
}

}  // namespace

extern "C"
void meth_orig_ref_slice(const meth_orig_ref_t *o,
                         int real_tid, int64_t st, int64_t en,
                         uint8_t *dst)
{
    if (o == NULL || dst == NULL || en <= st) return;
    static const char NT_ASCII[4] = {'A', 'C', 'G', 'T'};
    int64_t out_len = en - st;
    if (real_tid < 0 || real_tid >= o->n_chroms) {
        std::memset(dst, 'N', (size_t)out_len);
        return;
    }
    int64_t chrom_len = o->chrom_len[real_tid];
    int64_t f_off = o->bns->anns[o->f_rid[real_tid]].offset;
    int64_t r_off = o->bns->anns[o->r_rid[real_tid]].offset;

    /* Per-position dual decode: read f and r 2-bit bytes, fold via the
     * 5-row recovery table. OOB chrom-relative positions emit 'N'. */
    for (int64_t i = 0; i < out_len; ++i) {
        int64_t p = st + i;
        if (p < 0 || p >= chrom_len) { dst[i] = 'N'; continue; }
        uint8_t f = _get_pac(o->pac, f_off + p);
        uint8_t r = _get_pac(o->pac, r_off + p);
        uint8_t orig = recover_pair(f, r);
        dst[i] = (orig <= 3) ? (uint8_t)NT_ASCII[orig] : 'N';
    }

    /* Mask N over bns->ambs intervals overlapping our slice on EITHER
     * the f-contig or r-contig coordinates — both contigs in the doubled
     * pac correspond to the same real-chrom forward axis, and bwa indexes
     * Ns from both halves into bns->ambs at index time. */
    int64_t st_clipped = (st < 0 ? 0 : st);
    int64_t en_clipped = (en > chrom_len ? chrom_len : en);
    int64_t pac_len = en_clipped - st_clipped;
    if (pac_len > 0) {
        for (int side = 0; side < 2; ++side) {
            int64_t side_off = (side == 0) ? f_off : r_off;
            n_mask_ctx ctx = { dst, side_off, st };
            bns_iter_ambi(o->bns, side_off + st_clipped, (int)pac_len,
                          n_mask_visit, &ctx);
        }
    }
}

extern "C"
void meth_orig_ref_free(meth_orig_ref_t *o)
{
    if (o == NULL) return;
    free(o->f_rid);
    free(o->r_rid);
    free(o->chrom_len);
    free(o);
}
