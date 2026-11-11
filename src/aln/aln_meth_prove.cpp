/* SPDX-License-Identifier: MIT */
/* New code for bwa-mem3 (not from upstream bwa): exact refutations for `aln --meth`. See
 * aln_meth_prove.h. */
#include "aln_meth_prove.h"

#include <algorithm>

#include "../bntseq.h"

namespace bwa_aln {

static inline uint8_t comp_mask(uint8_t m)
{
    return (uint8_t)((m & 1) << 3 | (m & 2) << 1 | (m & 4) >> 1 | (m & 8) >> 3);
}

// The masks of converted forward positions [f0, f1) (0 <= f0 <= f1 <= conv l_pac), contig by contig.
static void fwd_masks(const aln_meth_ref_t &ref, int64_t f0, int64_t f1, uint8_t *out)
{
    const bntseq_t *cb = ref.conv_bns, *ob = ref.bns;
    for (int64_t f = f0; f < f1;) {
        const int crid = bns_pos2rid(cb, f);
        const bntann1_t &ca = cb->anns[crid];
        const int64_t cend = std::min<int64_t>(f1, ca.offset + ca.len);
        const int rid = crid >> 1, hyp = crid & 1;
        const int64_t o0 = ob->anns[rid].offset + (f - ca.offset), n = cend - f;
        // C->T copy (hyp 1): an original C may be indexed as T or C; G->A copy: a G as A or G.
        const int from = hyp ? 1 : 2;
        const uint8_t either = hyp ? (uint8_t)(1 << 1 | 1 << 3) : (uint8_t)(1 << 2 | 1 << 0);
        for (int64_t t = 0; t < n; ++t) {
            const int64_t p = o0 + t;
            const int b = ref.pac[p >> 2] >> ((~p & 3) << 1) & 3;
            out[t] = b == from ? either : (uint8_t)(1 << b);
        }
        // ambiguous original bases: any base (the converted index fills them with its own)
        int lo = 0, hi = ob->n_holes; // the first hole ending after o0
        while (lo < hi) {
            const int mid = (lo + hi) >> 1;
            if (ob->ambs[mid].offset + ob->ambs[mid].len <= o0) lo = mid + 1; else hi = mid;
        }
        for (int h = lo; h < ob->n_holes && ob->ambs[h].offset < o0 + n; ++h) {
            const int64_t b0 = std::max<int64_t>(o0, ob->ambs[h].offset);
            const int64_t b1 = std::min<int64_t>(o0 + n, ob->ambs[h].offset + ob->ambs[h].len);
            for (int64_t p = b0; p < b1; ++p) out[p - o0] = 0xf;
        }
        out += n;
        f = cend;
    }
}

void aln_meth_text_masks(const aln_meth_ref_t &ref, int64_t y0, int64_t y1, uint8_t *out)
{
    const int64_t lc = ref.conv_bns->l_pac;
    if (y0 < lc) fwd_masks(ref, y0, std::min(y1, lc), out);
    if (y1 > lc) { // T[p] = complement of forward 2 lc - 1 - p
        const int64_t p0 = std::max(y0, lc);
        uint8_t *o = out + (p0 - y0);
        const int64_t n = y1 - p0;
        fwd_masks(ref, 2 * lc - y1, 2 * lc - p0, o); // forward [2lc - y1, 2lc - p0), reversed below
        std::reverse(o, o + n);
        for (int64_t t = 0; t < n; ++t) o[t] = comp_mask(o[t]);
    }
}

// aln_meth_band_ed with the text read through more(b), which makes t[0, b) available (b <= m) and
// returns t: the masks are fetched only as far as the rows reach, which for a refuted window is
// usually a fraction of it.
template <class More>
static int band_ed_t(const ubyte_t *q, int n, int m, int dlo, int dhi, bool free_start, int kmax,
                     std::vector<int16_t> *scratch, More more)
{
    const int W = dhi - dlo + 1;
    const int16_t INF = 0x3fff;
    if (W <= 0) return kmax + 1;
    scratch->resize(2 * (size_t)W);
    int16_t *cur = scratch->data(), *nxt = cur + W;
    int16_t rmin = INF;
    for (int d = 0; d < W; ++d) {
        const int b = d + dlo;
        cur[d] = b < 0 || b > m ? INF : (int16_t)(free_start ? 0 : std::min<int>(b, INF));
        rmin = std::min(rmin, cur[d]);
    }
    if (rmin > kmax) return kmax + 1;
    for (int a = 1; a <= n; ++a) {
        const int qc = q[a - 1];
        const uint8_t *t = more(std::min(m, std::max(0, a + dhi)));
        rmin = INF;
        for (int d = 0; d < W; ++d) {
            const int b = a + d + dlo;
            if (b < 0 || b > m) { nxt[d] = INF; continue; }
            int v = INF;
            if (b >= 1) v = cur[d] + (qc < 4 && (t[b - 1] >> qc & 1) ? 0 : 1); // diagonal: D[a-1][b-1]
            if (d + 1 < W) v = std::min(v, cur[d + 1] + 1);                    // D[a-1][b]: a read base unaligned
            if (d > 0) v = std::min(v, nxt[d - 1] + 1);                        // D[a][b-1]: a text base unaligned
            nxt[d] = (int16_t)std::min<int>(v, INF);
            rmin = std::min(rmin, nxt[d]);
        }
        if (rmin > kmax) return kmax + 1;
        std::swap(cur, nxt);
    }
    int best = INF;
    for (int d = 0; d < W; ++d) {
        const int b = n + d + dlo;
        if (b >= 0 && b <= m) best = std::min<int>(best, cur[d]);
    }
    return std::min(best, kmax + 1);
}

int aln_meth_band_ed(const ubyte_t *q, int n, const uint8_t *t, int m, int dlo, int dhi, bool free_start, int kmax,
                     std::vector<int16_t> *scratch)
{
    return band_ed_t(q, n, m, dlo, dhi, free_start, kmax, scratch, [t](int) { return t; });
}

int aln_meth_refute_starts(const aln_meth_ref_t &ref, int64_t seq_len, const ubyte_t *seq, int len, int64_t y_lo,
                           int64_t y_hi, int kmax, aln_meth_verify_scratch_t *s)
{
    // An alignment starting at y in [y_lo, y_hi] keeps b - a within kmax of y - y_lo, so in
    // [-kmax, y_hi - y_lo + kmax] of the window from y_lo, and ends by y_hi + len + kmax.
    const int64_t w0 = std::max<int64_t>(0, y_lo), w1 = std::min<int64_t>(seq_len, y_hi + len + (int64_t)kmax);
    if (w1 <= w0 || y_hi < w0) return 1;
    const int c = (int)(w0 - y_lo), span = (int)(y_hi - y_lo);
    const int m = (int)(w1 - w0);
    s->text.resize(m);
    int have = 0;
    auto more = [&](int b) {
        if (b > have) { // in chunks, a few rows ahead
            const int to = std::min(m, std::max(b, have + 32));
            aln_meth_text_masks(ref, w0 + have, w0 + to, s->text.data() + have);
            have = to;
        }
        return (const uint8_t *)s->text.data();
    };
    return band_ed_t(seq, len, m, -kmax - c, span + kmax - c, true, kmax, &s->dp, more) > kmax;
}

int aln_meth_refute_window(const aln_meth_ref_t &ref, int64_t seq_len, const ubyte_t *seq, int len, int64_t y0, int kmax,
                           aln_meth_verify_scratch_t *s)
{
    return aln_meth_refute_starts(ref, seq_len, seq, len, y0, y0 + 2 * (int64_t)kmax, kmax, s);
}

} // namespace bwa_aln
