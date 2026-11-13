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
    if (ref.plain) { // the index text is the original pac; an ambiguous base may be any base
        for (int64_t p = f0; p < f1; ++p) out[p - f0] = (uint8_t)(1 << (ref.pac[p >> 2] >> ((~p & 3) << 1) & 3));
        int lo = 0, hi = ob->n_holes;
        while (lo < hi) {
            const int mid = (lo + hi) >> 1;
            if (ob->ambs[mid].offset + ob->ambs[mid].len <= f0) lo = mid + 1; else hi = mid;
        }
        for (int h = lo; h < ob->n_holes && ob->ambs[h].offset < f1; ++h) {
            const int64_t b0 = std::max<int64_t>(f0, ob->ambs[h].offset);
            const int64_t b1 = std::min<int64_t>(f1, ob->ambs[h].offset + ob->ambs[h].len);
            for (int64_t p = b0; p < b1; ++p) out[p - f0] = 0xf;
        }
        return;
    }
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

namespace bwa_aln {

int aln_meth_pigeonhole(const bwt_t *bwt, const aln_meth_ref_t &ref, const ubyte_t *seq, int len, int kmax,
                        const bwt_width_t *width, aln_meth_ph_scratch_t *s, long *steps, long *rows)
{
    enum { MAX_PIECES = 32 };
    *steps = 0; *rows = 0;
    // n1 pieces allowing one edit and n0 (0 or 1) exact ones: sum(e_j + 1) = 2 n1 + n0 = kmax + 1
    const int n1 = (kmax + 1) / 2, n0 = (kmax + 1) % 2, p = n1 + n0;
    s->why = ALN_METH_PH_SHORT;
    if (kmax < 0 || p > MAX_PIECES || len < p * (int)ALN_METH_PH_MIN_PIECE) return 0;
    int pa[MAX_PIECES + 1], pe[MAX_PIECES], off[MAX_PIECES + 1];
    for (int j = 0; j <= p; ++j) pa[j] = (int)((int64_t)j * len / p);
    for (int j = 0; j < p; ++j) pe[j] = j < n1 ? 1 : 0;
    off[0] = 0;
    for (int j = 0; j < p; ++j) off[j + 1] = off[j] + (pa[j + 1] - pa[j]);
    // Each piece's absent-substring bound: bid[off[j] + t] bounds the edits S[a_j, a_j + t] needs. It
    // is bwt_cal_width's walk from a_j (on the uncomplemented bases: the index holds both strands, so
    // a string and its reverse complement occur equally), the pieces' walks run in step. Piece 0 is
    // the search's own width walk when given.
    s->bid.resize(off[p]);
    struct walk_t { bwtint_t k, l; int t, bid; };
    walk_t w[MAX_PIECES];
    int n_walk = 0;
    for (int j = 0; j < p; ++j) {
        if (j == 0 && width) { for (int t = 0; t < pa[1]; ++t) s->bid[t] = width[t].bid; continue; }
        w[j] = walk_t{0, bwt->seq_len, 0, 0};
        ++n_walk;
    }
    auto vbase = [&](int at) { return seq[at] > 3 ? 4 : 3 - seq[at]; };
    while (n_walk) {
        for (int j = 0; j < p; ++j)
            if (!(j == 0 && width) && w[j].t < pa[j + 1] - pa[j]) bwt_prefetch_2occ(bwt, w[j].k, w[j].l);
        for (int j = 0; j < p; ++j) {
            if ((j == 0 && width) || w[j].t >= pa[j + 1] - pa[j]) continue;
            walk_t &x = w[j];
            const int c = vbase(pa[j] + x.t);
            if (c < 4) {
                bwtint_t ok, ol;
                bwt_2occ(bwt, x.k - 1, x.l, (ubyte_t)c, &ok, &ol);
                ++*steps;
                x.k = bwt->L2[c] + ok + 1;
                x.l = bwt->L2[c] + ol;
            }
            if (x.k > x.l || c > 3) { x.k = 0; x.l = bwt->seq_len; ++x.bid; }
            s->bid[off[j] + x.t] = x.bid;
            if (++x.t == pa[j + 1] - pa[j]) --n_walk;
        }
    }
    // The pieces' occurrences within their allowances: backward search from each piece's last base
    // (S[j] is the next base a node consumes), all nodes of all pieces advanced one lookup per round.
    // A node may still place used < e edits: a mismatch, an unaligned read base (insertion) or an
    // unaligned text base (deletion). Deletions at the piece's right end are not made: dropping that
    // text base leaves an occurrence with fewer edits at the same start. Nor at its left end, after
    // its first base: that occurrence starts one base later with one edit fewer, which the window
    // allows for (aln_meth_refute_window takes starts within kmax of the diagonal either way).
    s->cur.clear();
    s->occ_k.clear(); s->occ_l.clear(); s->occ_a.clear();
    int64_t total = 0;
    for (int j = 0; j < p; ++j)
        if (s->bid[off[j + 1] - 1] <= pe[j]) s->cur.push_back(aln_meth_ph_node_t{0, bwt->seq_len, pa[j + 1] - 1, (uint8_t)j, 0});
    while (!s->cur.empty()) {
        for (const auto &n : s->cur) bwt_prefetch_2occ(bwt, n.k, n.l);
        s->next.clear();
        for (const auto &n : s->cur) {
            const int A = pa[n.piece], e = pe[n.piece], j = n.j, q = seq[j];
            auto bound = [&](int jj) { return jj < A ? 0 : s->bid[off[n.piece] + (jj - A)]; };
            auto emit = [&](bwtint_t k, bwtint_t l) {
                s->occ_k.push_back(k); s->occ_l.push_back(l); s->occ_a.push_back(A);
                total += (int64_t)(l - k + 1);
            };
            auto push = [&](bwtint_t k, bwtint_t l, int jj, int u) {
                s->next.push_back(aln_meth_ph_node_t{k, l, jj, n.piece, (uint8_t)u});
            };
            if (n.used == e) { // no edit left: the rest of the piece exactly
                if (q > 3) continue;
                bwtint_t ok, ol;
                bwt_2occ(bwt, n.k - 1, n.l, (ubyte_t)q, &ok, &ol);
                ++*steps;
                const bwtint_t k = bwt->L2[q] + ok + 1, l = bwt->L2[q] + ol;
                if (k > l) continue;
                if (j == A) emit(k, l);
                else if (bound(j - 1) == 0) push(k, l, j - 1, n.used);
                continue;
            }
            bwtint_t ck[4], cl[4];
            bwt_2occ4(bwt, n.k - 1, n.l, ck, cl);
            ++*steps;
            const bool root = j == pa[n.piece + 1] - 1 && n.used == 0;
            const int left = e - n.used - 1; // edits left after placing one here
            for (int c = 0; c < 4; ++c) {
                const bwtint_t k = bwt->L2[c] + ck[c] + 1, l = bwt->L2[c] + cl[c];
                if (k > l) continue;
                const int u = n.used + (c != q); // a read N is a mismatch with every base
                if (u <= e) {
                    if (j == A) emit(k, l);
                    else if (bound(j - 1) <= e - u) push(k, l, j - 1, u);
                }
                if (!root && left >= 0 && bound(j) <= left) push(k, l, j, n.used + 1); // text base c unaligned
            }
            if (left >= 0) { // read base S[j] unaligned
                if (j == A) emit(n.k, n.l);
                else if (bound(j - 1) <= left) push(n.k, n.l, j - 1, n.used + 1);
            }
        }
        if (total > ALN_METH_PH_ROWS) { s->why = ALN_METH_PH_ROWS_CAP; *rows = total; return 0; }
        if (*steps > ALN_METH_PH_NODES) { s->why = ALN_METH_PH_NODES_CAP; return 0; }
        std::swap(s->cur, s->next);
    }
    // Every occurrence row: the starts S may have there, checked against the text.
    *rows = total;
    s->rows.clear();
    std::vector<int> &row_a = s->occ_a; // reused below as each row's piece start
    std::vector<int64_t> &st = s->starts;
    st.clear();
    for (size_t o = 0; o < s->occ_k.size(); ++o)
        for (bwtint_t r = s->occ_k[o]; r <= s->occ_l[o]; ++r) s->rows.push_back(r);
    bwt_sa_batch(bwt, s->rows.data(), (int64_t)s->rows.size(), s->rows.data());
    size_t at = 0;
    for (size_t o = 0; o < s->occ_k.size(); ++o)
        for (bwtint_t r = s->occ_k[o]; r <= s->occ_l[o]; ++r, ++at) {
            const bwtint_t y = s->rows[at];
            if (y >= bwt->seq_len) { s->why = ALN_METH_PH_NOT_TEXT; return 0; } // not a text position
            st.push_back((int64_t)y - row_a[o] - kmax);
        }
    std::sort(st.begin(), st.end());
    st.erase(std::unique(st.begin(), st.end()), st.end());
    // Rows of one locus (several pieces, or one piece's edit placements) give nearby starts: their
    // start ranges, [y0, y0 + 2 kmax] each, are checked together while they overlap, up to a span
    // of 6 kmax (one banded pass over their union).
    for (size_t i = 0; i < st.size();) {
        const int64_t lo = st[i];
        int64_t hi = st[i] + 2 * (int64_t)kmax;
        size_t j = i + 1;
        while (j < st.size() && st[j] <= hi + 1 && st[j] + 2 * (int64_t)kmax - lo <= 6 * (int64_t)kmax) hi = st[j++] + 2 * (int64_t)kmax;
        if (!aln_meth_refute_starts(ref, (int64_t)bwt->seq_len, seq, len, lo, hi, kmax, &s->vs)) { s->why = ALN_METH_PH_WINDOW; return 0; }
        i = j;
    }
    s->why = ALN_METH_PH_PROVEN;
    return 1;
}

} // namespace bwa_aln
