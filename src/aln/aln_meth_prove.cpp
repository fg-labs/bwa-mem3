/* SPDX-License-Identifier: MIT */
/* New code for bwa-mem3 (not from upstream bwa): exact refutations for `aln --meth`. See
 * aln_meth_prove.h. */
#include "aln_meth_prove.h"

#include <algorithm>
#include <cstring>

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

/* Text masks as four bit planes: bit (p + PAD) of plane c is set when base c may be at text
 * position p. Positions below 0, and those not filled yet, read as 0. */
namespace {
struct mask_planes_t {
    enum { PAD = 64 };
    uint64_t *w = nullptr; // plane c is w[c * nw, (c + 1) * nw)
    int nw = 0;
    // Planes for positions below n, all clear, in *store.
    void init(std::vector<uint64_t> *store, int n)
    {
        nw = (n + PAD) / 64 + 2;
        store->assign(4 * (size_t)nw, 0);
        w = store->data();
    }
    // Masks t[from, to) into the planes; from is a multiple of 8, and t holds zeros from to up to
    // the next multiple of 8.
    void fill(const uint8_t *t, int from, int to)
    {
        for (int p = from; p < to; p += 8) {
            uint64_t x;
            memcpy(&x, t + p, 8);
            const int bi = p + PAD;
            for (int c = 0; c < 4; ++c) {
                // bit c of each byte, gathered into one byte (byte i's bit lands at bit i)
                const uint64_t y = (x >> c) & 0x0101010101010101ULL;
                w[c * nw + (bi >> 6)] |= ((y * 0x0102040810204080ULL) >> 56) << (bi & 63);
            }
        }
    }
    // Bit d: base c may be at position pos + d, for d in [0, 64).
    uint64_t bits(int c, int pos) const
    {
        const int bi = pos + PAD;
        if (bi < 0) return 0; // every position below 0
        const uint64_t *pl = w + c * nw + (bi >> 6);
        const int sh = bi & 63;
        return sh ? pl[0] >> sh | pl[1] << (64 - sh) : pl[0];
    }
};
static_assert(__BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__, "mask_planes_t::fill reads 8 masks as one little-endian word");
} // namespace

// Bits lo..hi of a W-bit word (empty when lo > hi after clamping to [0, W)).
static inline uint64_t bit_range(int lo, int hi, int W)
{
    lo = std::max(lo, 0);
    hi = std::min(hi, W - 1);
    if (lo > hi) return 0;
    const uint64_t up = hi >= 63 ? ~0ULL : (1ULL << (hi + 1)) - 1;
    return up & ~0ULL << lo;
}

/* band_ed_t on bit vectors: the same cells, the same result and the same early exit, for bands of
 * at most 64 diagonals and kmax <= ALN_METH_BITS_KMAX. Bit d of V[t] is set when the cell on
 * diagonal d of the current row, D[a][a + d + dlo], is at most t (Wu and Manber's recurrence, one
 * word per distance up to kmax): a cell is at most t when its diagonal predecessor is at most t and
 * the bases match, or when its diagonal or upper predecessor (row a - 1) or its left neighbour (row
 * a) is at most t - 1. Distances above kmax are never needed: such a cell only feeds cells above
 * kmax. more(b): as for band_ed_t, makes text positions [0, b) of the planes readable. */
template <class More>
static int band_ed_bits(const ubyte_t *q, int n, int m, int dlo, int dhi, bool free_start, int kmax, More more)
{
    const int W = dhi - dlo + 1, K = kmax;
    if (W <= 0) return kmax + 1;
    uint64_t va[ALN_METH_BITS_KMAX + 1], vb[ALN_METH_BITS_KMAX + 1];
    uint64_t *V = va, *N = vb;
    uint64_t valid = bit_range(-dlo, m - dlo, W); // row 0: b = d + dlo in [0, m]
    for (int t = 0; t <= K; ++t) V[t] = free_start ? valid : valid & bit_range(-dlo, t - dlo, W); // D[0][b] = b
    if (!V[K]) return K + 1;
    for (int a = 1; a <= n; ++a) {
        const mask_planes_t &P = more(std::min(m, std::max(0, a + dhi)));
        const int qc = q[a - 1];
        valid = bit_range(-a - dlo, m - a - dlo, W);                  // b = a + d + dlo in [0, m]
        const uint64_t diag = valid & bit_range(1 - a - dlo, W - 1, W); // b >= 1: a diagonal predecessor
        const uint64_t match = qc < 4 ? P.bits(qc, a + dlo - 1) & diag : 0; // t[b - 1] holds base qc
        uint64_t left = 0; // N[t - 1]
        for (int t = 0; t <= K; ++t) {
            uint64_t x = V[t] & match;
            if (t) x |= (V[t - 1] & diag) | V[t - 1] >> 1 | left << 1;
            N[t] = left = x & valid;
        }
        if (!N[K]) return K + 1; // every cell of the row is above kmax
        std::swap(V, N);
    }
    for (int t = 0; t <= K; ++t)
        if (V[t]) return t;
    return K + 1;
}

int aln_meth_band_ed_bits(const ubyte_t *q, int n, const uint8_t *t, int m, int dlo, int dhi, bool free_start, int kmax,
                          aln_meth_verify_scratch_t *s)
{
    if (dhi - dlo + 1 > 64 || kmax < 0 || kmax > ALN_METH_BITS_KMAX) return aln_meth_band_ed(q, n, t, m, dlo, dhi, free_start, kmax, &s->dp);
    s->text.assign((size_t)m + 8, 0);
    std::copy(t, t + m, s->text.begin());
    mask_planes_t P;
    P.init(&s->planes, std::max(m + 8, n + dhi + 64));
    P.fill(s->text.data(), 0, m);
    return band_ed_bits(q, n, m, dlo, dhi, free_start, kmax, [&P](int) -> const mask_planes_t & { return P; });
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
    const int dlo = -kmax - c, dhi = span + kmax - c;
    const bool bits = dhi - dlo + 1 <= 64 && kmax >= 0 && kmax <= ALN_METH_BITS_KMAX;
    s->text.assign((size_t)m + 8, 0);
    mask_planes_t planes;
    if (bits) planes.init(&s->planes, std::max(m + 8, len + dhi + 64));
    int have = 0;
    auto more = [&](int b) {
        if (b > have) { // in chunks, a few rows ahead
            const int to = std::min(m, std::max((b + 31) & ~31, have + 32));
            aln_meth_text_masks(ref, w0 + have, w0 + to, s->text.data() + have);
            if (bits) planes.fill(s->text.data(), have, to);
            have = to;
        }
        return (const uint8_t *)s->text.data();
    };
    if (bits) {
        auto more_planes = [&](int b) -> const mask_planes_t & { more(b); return planes; };
        return band_ed_bits(seq, len, m, dlo, dhi, true, kmax, more_planes) > kmax;
    }
    return band_ed_t(seq, len, m, dlo, dhi, true, kmax, &s->dp, more) > kmax;
}

// Prefetch the original-reference bytes a converted-text window can read, including boundaries.
static void prefetch_window(const aln_meth_ref_t &ref, int64_t lo, int64_t hi)
{
    const int64_t lc = ref.conv_bns->l_pac;
    lo = std::max<int64_t>(0, lo);
    hi = std::min<int64_t>(2 * lc, hi);
    while (lo < hi) {
        const bool reverse = lo >= lc;
        const int64_t f = reverse ? 2 * lc - 1 - lo : lo;
        const int crid = bns_pos2rid(ref.conv_bns, f), rid = crid >> 1;
        const bntann1_t &c = ref.conv_bns->anns[crid];
        const int64_t n = std::min(hi - lo, reverse ? f - c.offset + 1 : c.offset + c.len - f);
        const int64_t p = ref.bns->anns[rid].offset + f - c.offset;
        const int64_t first = reverse ? p - n + 1 : p;
        const int64_t last = reverse ? p : p + n - 1;
        for (int64_t b = first >> 2; b <= last >> 2; b += 64) __builtin_prefetch(ref.pac + b);
        __builtin_prefetch(ref.pac + (last >> 2));
        lo += n;
    }
}

int aln_meth_refute_window(const aln_meth_ref_t &ref, int64_t seq_len, const ubyte_t *seq, int len, int64_t y0, int kmax,
                           aln_meth_verify_scratch_t *s)
{
    return aln_meth_refute_starts(ref, seq_len, seq, len, y0, y0 + 2 * (int64_t)kmax, kmax, s);
}

} // namespace bwa_aln

namespace bwa_aln {

int aln_meth_pigeonhole(const bwt_t *bwt, const aln_meth_ref_t &ref, const ubyte_t *seq, int len, int kmax,
                        const bwt_width_t *width, aln_meth_ph_scratch_t *s, long *steps, long *rows,
                        int64_t max_rows, long max_nodes)
{
    enum { MAX_PIECES = 32 };
    *steps = 0; *rows = 0;
    s->why = ALN_METH_PH_SHORT;
    // kmax >= 2 * MAX_PIECES needs more than MAX_PIECES pieces; refused before kmax + 1 can overflow
    if (kmax < 0 || kmax >= 2 * MAX_PIECES) return 0;
    // n1 pieces allowing one edit and n0 (0 or 1) exact ones: sum(e_j + 1) = 2 n1 + n0 = kmax + 1
    const int n1 = (kmax + 1) / 2, n0 = (kmax + 1) % 2, p = n1 + n0;
    if (p > MAX_PIECES || len < p * (int)ALN_METH_PH_MIN_PIECE) return 0;
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
    // Exact occurrences are a subset of the enumerated rows, even for a one-edit piece. A
    // hopeless attempt can stop before constructing its approximate-search frontier.
    int64_t exact_rows = 0;
    for (int j = 0; j < p; ++j) {
        if (s->bid[off[j + 1] - 1] != 0) continue;
        exact_rows += j == 0 && width ? (int64_t)width[pa[1] - 1].w : (int64_t)(w[j].l - w[j].k + 1);
        if (exact_rows > max_rows) { s->why = ALN_METH_PH_ROWS_CAP; *rows = exact_rows; return 0; }
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
        if (total > max_rows) { s->why = ALN_METH_PH_ROWS_CAP; *rows = total; return 0; }
        if (*steps > max_nodes) { s->why = ALN_METH_PH_NODES_CAP; return 0; }
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
    size_t prefetched = 0;
    for (size_t i = 0; i < st.size();) {
        while (prefetched < std::min(st.size(), i + 8)) {
            const int64_t start = st[prefetched++];
            prefetch_window(ref, start, start + len + 3 * (int64_t)kmax);
        }
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
