/* SPDX-License-Identifier: MIT */

#include "meth_xm.h"

#include <cstdint>
#include <cstdlib>
#include <cstring>

#include "bntseq.h"
#include "kstring.h"

/* Thread-local growable scratch for the per-record ref-window slice and
 * the assembled XM:Z payload. Reused across calls; capacity grows on
 * demand via ks_resize and is freed at thread exit by the RAII wrapper.
 *
 * The kstring_t pattern matches the rest of meth_bam.cpp (XA, SA tag
 * builders), but those use stack-allocated kstrings that malloc + free
 * per record. Promoting them to thread_local removes the per-record
 * allocator round-trip on the BAM-emission hot path. */
namespace {
struct ks_scratch_t {
    kstring_t buf{0, 0, nullptr};
    ~ks_scratch_t() { free(buf.s); }
};
}
static thread_local ks_scratch_t t_ref_scratch;
static thread_local ks_scratch_t t_xm_scratch;

/* Per-source-strand rule table, indexed by is_top_strand (1 = top, 0 = bottom).
 *
 * For top-strand methylation (XG:Z:CT, i.e. f-contig alignment): the C of
 * interest is the forward-strand C; methylated reads as 'C', unmethylated
 * as 'T'; downstream context on top is forward[P+1], forward[P+2]; CpG
 * marker downstream is 'G'. Note the SAM 0x10 (RC) flag is independent —
 * CTOT reads (R2 mapped to top with 0x10 set) have is_top_strand=1.
 *
 * For bottom-strand methylation (XG:Z:GA, i.e. r-contig alignment): the
 * C of interest is a forward G (bottom-strand C); methylated reads as
 * 'G' (after BAM SEQ orientation), unmethylated as 'A'; downstream context
 * on bottom is upstream on forward — forward[P-1], forward[P-2]; CpG
 * marker downstream-on-bottom is 'C' on forward (complement of G).
 *
 * Both strand cases consume seq_text and ref_buf in SEQ orientation
 * (= forward-genome orientation per BAM convention).
 */
static const char C_MARKER[2]       = {'G', 'C'};  /* [bottom, top] */
static const char METH_MARKER[2]    = {'G', 'C'};
static const char UNMETH_MARKER[2]  = {'A', 'T'};
static const char CTX_CPG_MARKER[2] = {'C', 'G'};

/* D3 (--meth, PR-5): slice forward-genome ASCII bases [st, en) of contig
 * `real_tid` directly from the ORIGINAL (un-converted) bns + 2-bit pac,
 * replacing the retired meth_orig_ref f/r fold. `dst` capacity must be
 * >= en - st. OOB chrom-relative positions and bns->ambs overlaps fill
 * 'N'. `real_tid` indexes bns->anns (the alignment's original rid). */
namespace {

struct n_mask_ctx {
    uint8_t *dst;
    int64_t  pac_off;   /* contig's global pac offset */
    int64_t  st_dst;    /* caller's chrom-relative slice start */
    int64_t  out_len;   /* dst length; bounds the write window */
};

int n_mask_visit(int64_t amb_st, int64_t amb_en, void *vctx)
{
    n_mask_ctx *c = (n_mask_ctx *)vctx;
    int64_t lo = amb_st - c->pac_off - c->st_dst;
    int64_t hi = amb_en - c->pac_off - c->st_dst;
    /* bns_iter_ambi already clips [amb_st, amb_en) to its query window (see
     * bns_iter_ambi: st=max(offset,pos_f), en=min(offset+len,window_end)), so
     * lo>=0 and hi<=out_len here. Clamp anyway: this is a buffer write and the
     * guarantee lives in another function — keep the callback self-contained. */
    if (lo < 0) lo = 0;
    if (hi > c->out_len) hi = c->out_len;
    for (int64_t j = lo; j < hi; ++j) c->dst[j] = 'N';
    return 0;
}

void orig_ref_slice(const bntseq_t *bns, const uint8_t *pac,
                    int real_tid, int64_t st, int64_t en, uint8_t *dst)
{
    if (bns == NULL || pac == NULL || dst == NULL || en <= st) return;
    static const char NT_ASCII[4] = {'A', 'C', 'G', 'T'};
    int64_t out_len = en - st;
    if (real_tid < 0 || real_tid >= bns->n_seqs) {
        std::memset(dst, 'N', (size_t)out_len);
        return;
    }
    int64_t chrom_len = bns->anns[real_tid].len;
    int64_t off       = bns->anns[real_tid].offset;
    for (int64_t i = 0; i < out_len; ++i) {
        int64_t p = st + i;
        if (p < 0 || p >= chrom_len) { dst[i] = 'N'; continue; }
        uint8_t b = (uint8_t)_get_pac(pac, off + p);
        dst[i] = (b <= 3) ? (uint8_t)NT_ASCII[b] : 'N';
    }
    int64_t st_clipped = (st < 0 ? 0 : st);
    int64_t en_clipped = (en > chrom_len ? chrom_len : en);
    int64_t pac_len = en_clipped - st_clipped;
    if (pac_len > 0) {
        n_mask_ctx ctx = { dst, off, st, out_len };
        bns_iter_ambi(bns, off + st_clipped, (int)pac_len, n_mask_visit, &ctx);
    }
}

}  // namespace

extern "C"
char *meth_build_xm(const bntseq_t *bns, const uint8_t *pac, int real_tid,
                    int64_t pos, int is_top_strand,
                    const uint32_t *bam_cigar, int n_cigar,
                    const char *seq_text, int l_emit)
{
    if (bns == NULL || pac == NULL || bam_cigar == NULL
            || seq_text == NULL || l_emit <= 0)
        return NULL;

    /* Compute ref length on bam_cigar (M/D/N/=/X consume ref). */
    int64_t rlen = 0;
    for (int i = 0; i < n_cigar; ++i) {
        int op = (int)(bam_cigar[i] & 0xf);
        int len = (int)(bam_cigar[i] >> 4);
        if (op == 0 /*M*/ || op == 2 /*D*/ || op == 3 /*N*/
                || op == 7 /*=*/ || op == 8 /*X*/)
            rlen += len;
    }

    /* Slice forward-genome ref [pos - 2, pos + rlen + 2). The +/-2 covers
     * both is_rev=0 (downstream context lookahead by +2) and is_rev=1
     * (downstream-on-bottom lookahead by -2 on forward). OOB filled with N
     * by orig_ref_slice. Both ref_buf and xm reuse thread-local
     * kstring_t scratch — no per-record malloc/free on the hot path. */
    int64_t buf_len = rlen + 4;
    ks_resize(&t_ref_scratch.buf, (size_t)buf_len);
    if (t_ref_scratch.buf.s == NULL) return NULL;
    uint8_t *ref_buf = (uint8_t *)t_ref_scratch.buf.s;
    orig_ref_slice(bns, pac, real_tid, pos - 2, pos + rlen + 2, ref_buf);

    ks_resize(&t_xm_scratch.buf, (size_t)l_emit + 1);
    if (t_xm_scratch.buf.s == NULL) return NULL;
    char *xm = t_xm_scratch.buf.s;

    int r = 0;        /* read cursor (index into seq_text and into xm) */
    int64_t t = 0;    /* ref cursor (offset from `pos` on forward strand) */
    int idx = is_top_strand ? 1 : 0;
    char c_marker      = C_MARKER[idx];
    char meth_marker   = METH_MARKER[idx];
    char unmeth_marker = UNMETH_MARKER[idx];
    char cpg_marker    = CTX_CPG_MARKER[idx];

    for (int i = 0; i < n_cigar && r < l_emit; ++i) {
        int op = (int)(bam_cigar[i] & 0xf);
        int len = (int)(bam_cigar[i] >> 4);
        switch (op) {
            case 0: /* M */
            case 7: /* = */
            case 8: /* X */
                for (int k = 0; k < len && r < l_emit; ++k) {
                    int64_t ref_idx = t + 2;  /* into ref_buf (which starts -2) */
                    char rb = (char)ref_buf[ref_idx];
                    char read_b = seq_text[r];
                    if (rb != c_marker) {
                        xm[r] = '.';
                    } else {
                        int64_t ctx1_idx = is_top_strand ? ref_idx + 1 : ref_idx - 1;
                        int64_t ctx2_idx = is_top_strand ? ref_idx + 2 : ref_idx - 2;
                        char ctx1 = (ctx1_idx >= 0 && ctx1_idx < buf_len)
                                        ? (char)ref_buf[ctx1_idx] : 'N';
                        char ctx2 = (ctx2_idx >= 0 && ctx2_idx < buf_len)
                                        ? (char)ref_buf[ctx2_idx] : 'N';
                        char meth_state;
                        if (read_b == meth_marker)        meth_state = 'M';
                        else if (read_b == unmeth_marker) meth_state = 'U';
                        else                              meth_state = '.';
                        if (meth_state == '.') {
                            xm[r] = '.';
                        } else if (ctx1 == 'N' || ctx1 == 'X') {
                            xm[r] = (meth_state == 'M') ? 'U' : 'u';
                        } else if (ctx1 == cpg_marker) {
                            xm[r] = (meth_state == 'M') ? 'Z' : 'z';
                        } else if (ctx2 == 'N' || ctx2 == 'X') {
                            xm[r] = (meth_state == 'M') ? 'U' : 'u';
                        } else if (ctx2 == cpg_marker) {
                            xm[r] = (meth_state == 'M') ? 'X' : 'x';
                        } else {
                            xm[r] = (meth_state == 'M') ? 'H' : 'h';
                        }
                    }
                    ++r;
                    ++t;
                }
                break;
            case 1: /* I */
            case 4: /* S */
                for (int k = 0; k < len && r < l_emit; ++k) {
                    xm[r++] = '.';
                }
                break;
            case 2: /* D */
            case 3: /* N */
                t += len;
                break;
            case 5: /* H */
            case 6: /* P */
            default:
                /* No read or ref consume. */
                break;
        }
    }

    /* Defensive padding if cigar walk ended short of l_emit. */
    while (r < l_emit) xm[r++] = '.';
    xm[l_emit] = '\0';
    return xm;
}
