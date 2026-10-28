/* SPDX-License-Identifier: MIT */

#include "meth_model.h"

void meth_model_fill_mats(mem_opt_t *o)
{
    for (int r = 0; r < METH_M; ++r) {
        for (int c = 0; c < METH_M; ++c) {
            int8_t v = o->mat[meth_code_base((uint8_t)r) * 5 + meth_code_base((uint8_t)c)];
            o->mat7_ot[r * METH_M + c] = v;
            o->mat7_ob[r * METH_M + c] = v;
        }
    }
    /* Same conversion value as mem_opt_fill_meth_mat, on the CpG rows only. */
    int8_t conv = (o->meth_scoring == MEM_METH_SCORING_NEUTRAL) ? 0 : (int8_t)o->a;
    o->mat7_ot[METH_C5 * METH_M + METH_T] = conv;   /* OT: ref CpG C / read T */
    o->mat7_ob[METH_G6 * METH_M + METH_A] = conv;   /* OB: ref CpG G / read A */
}

void meth_unlabel(uint8_t *win, int64_t len)
{
    for (int64_t i = 0; i < len; ++i) win[i] = meth_code_base(win[i]);
}

namespace {

/* Frame geometry of one window. Frame index k (k may be -1 or len, the two
 * outside neighbours) is forward base rb + k on the forward strand and
 * fbase0 - k on the reverse strand, where it reads as the complement. */
struct frame_t {
    int     rev;
    int64_t rb;       /* forward strand: forward coordinate of frame index 0 */
    int64_t fbase0;   /* reverse strand: forward coordinate of frame index 0 */
    int64_t len;
};

inline int64_t frame_to_fwd(const frame_t &f, int64_t k)
{
    return f.rev ? f.fbase0 - k : f.rb + k;
}

inline int64_t fwd_to_frame(const frame_t &f, int64_t p)
{
    return f.rev ? f.fbase0 - p : p - f.rb;
}

struct unlabel_ctx {
    const frame_t *f;
    uint8_t       *win;
};

/* An ambs run [st, en) (forward coordinates): none of its bases, nor the
 * window bases they would pair with, is a CpG. A label at k-1 (C5) depends
 * only on base k, and one at k+1 (G6) only on base k, so undoing those three
 * is exact. */
int unlabel_amb_visit(int64_t st, int64_t en, void *vctx)
{
    unlabel_ctx *c = (unlabel_ctx *)vctx;
    const int64_t len = c->f->len;
    for (int64_t p = st; p < en; ++p) {
        int64_t k = fwd_to_frame(*c->f, p);
        if (k >= 0 && k < len) c->win[k] = meth_code_base(c->win[k]);
        if (k - 1 >= 0 && k - 1 < len && c->win[k - 1] == METH_C5) c->win[k - 1] = METH_C;
        if (k + 1 >= 0 && k + 1 < len && c->win[k + 1] == METH_G6) c->win[k + 1] = METH_G;
    }
    return 0;
}

}  // namespace

void meth_label_cpg(uint8_t *win, int64_t len, int64_t rb, int rid,
                    const bntseq_t *bns, const uint8_t *pac)
{
    if (win == NULL || len <= 0 || bns == NULL || pac == NULL) return;
    if (rid < 0 || rid >= bns->n_seqs) return;
    const int64_t l_pac = bns->l_pac;
    frame_t f;
    f.rev    = rb >= l_pac;
    f.rb     = rb;
    f.fbase0 = (l_pac << 1) - 1 - rb;
    f.len    = len;

    /* Contig rid as a half-open range [klo, khi) of frame indices. */
    const int64_t off = bns->anns[rid].offset, end = off + bns->anns[rid].len;
    const int64_t klo = f.rev ? f.fbase0 - end + 1 : off - rb;
    const int64_t khi = f.rev ? f.fbase0 - off + 1 : end - rb;
    const int64_t lo = klo > 0 ? klo : 0, hi = khi < len ? khi : len;
    if (lo >= hi) return;

    /* Frame base at k in [lo-1, hi], or -1 outside the contig. Inside the
     * window it is the window's own code; outside, the pac base in frame. */
    auto base_at = [&](int64_t k) -> int {
        if (k < klo || k >= khi) return -1;
        if (k >= 0 && k < len) return meth_code_base(win[k]);
        int b = (int)_get_pac(pac, frame_to_fwd(f, k));
        return f.rev ? 3 - b : b;
    };

    /* Each label depends only on plain bases, so labelling left to right in
     * place is safe (base_at strips a label already written at k - 1). */
    int prev = base_at(lo - 1);
    for (int64_t i = lo; i < hi; ++i) {
        const int cur = meth_code_base(win[i]);
        const int next = base_at(i + 1);
        if (cur == METH_C && next == METH_G)      win[i] = METH_C5;
        else if (cur == METH_G && prev == METH_C) win[i] = METH_G6;
        prev = cur;
    }

    /* Undo every label an ambs run touches: frame indices [lo-1, hi], clipped
     * to the contig, as a forward range. */
    int64_t k0 = lo - 1 > klo ? lo - 1 : klo;
    int64_t k1 = hi + 1 < khi ? hi + 1 : khi;   /* exclusive */
    int64_t p0 = frame_to_fwd(f, k0), p1 = frame_to_fwd(f, k1 - 1);
    if (p0 > p1) { int64_t t = p0; p0 = p1; p1 = t; }
    unlabel_ctx ctx = { &f, win };
    bns_iter_ambi(bns, p0, (int)(p1 - p0 + 1), unlabel_amb_visit, &ctx);
}
