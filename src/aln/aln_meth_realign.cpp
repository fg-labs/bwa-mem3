/* SPDX-License-Identifier: MIT */
#include "aln_meth_realign.h"
#include "aln_known_snvs.h"

#include <algorithm>
#include "../bwa.h"
#include "../ksw.h"
#include "../meth_model.h"

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <memory>
#include <new>

namespace bwa_aln {
namespace {

aln_meth_eval_status error(aln_meth_eval_status status, const char *message)
{
    std::fprintf(stderr, "ERROR: --meth candidate realignment: %s\n", message);
    return status;
}

// The refinement's gap model: a gap of length k costs gap_open + k * gap_extend, for both
// insertions and deletions.
constexpr int gap_open = 5, gap_extend = 1;

/* Whether the ungapped alignment of an equal-length query and window is the unique optimum of the
 * global alignment, so ksw_global would return it ([len M]) at any band. *reward gets its reward.
 * With equal lengths, any gapped path has at least one insertion and one deletion (each at least
 * gap_open + gap_extend) and at most len - 1 aligned columns, each worth at most the matrix's
 * largest entry (or 0 when that is negative), so an ungapped reward above that bound wins
 * strictly. */
bool ungapped_is_optimal(const int8_t *mat, int m, const ubyte_t *query, const uint8_t *target, int len,
                         int *reward)
{
    int best_cell = 0, sum = 0;
    for (int i = 0; i < m * m; ++i) best_cell = std::max<int>(best_cell, mat[i]);
    for (int i = 0; i < len; ++i) sum += mat[target[i] * m + query[i]];
    *reward = sum;
    const int gapped_bound = best_cell * (len - 1) - 2 * (gap_open + gap_extend);
    return sum > gapped_bound;
}

struct ambiguity_window_t { uint8_t *bases; int64_t rb; };

int mask_ambiguity(int64_t begin, int64_t end, void *context)
{
    const auto *window = static_cast<ambiguity_window_t *>(context);
    std::fill(window->bases + begin - window->rb, window->bases + end - window->rb, 4);
    return 0;
}

} // namespace

bool aln_meth_window_has_site(const aln_meth_ref_t &ref, int64_t rb, int len, int rid)
{
    if (!ref.known_snvs || rid < 0 || rid >= (int)ref.known_snvs->by_contig.size()) return false;
    const std::vector<ks_site_t> &lst = ref.known_snvs->by_contig[rid];
    const int64_t beg = rb - ref.bns->anns[rid].offset;
    auto it = std::lower_bound(lst.begin(), lst.end(), beg, [](const ks_site_t &x, int64_t p) { return x.pos < p; });
    return it != lst.end() && it->pos < beg + len;
}

void aln_meth_fetch_window(const aln_meth_ref_t &ref, const aln_meth_opt_t &mo,
    int64_t rb, int len, int rid, std::vector<uint8_t> &target, int alt)
{
    target.resize(len);
    for (int i = 0; i < len; ++i) {
        const int64_t p = rb + i;
        target[i] = ref.pac[p >> 2] >> ((~p & 3) << 1) & 3;
    }
    if (alt && ref.known_snvs && rid >= 0 && rid < (int)ref.known_snvs->by_contig.size()) {
        const std::vector<ks_site_t> &lst = ref.known_snvs->by_contig[rid];
        const int64_t beg = rb - ref.bns->anns[rid].offset;
        auto it = std::lower_bound(lst.begin(), lst.end(), beg, [](const ks_site_t &x, int64_t p) { return x.pos < p; });
        for (int64_t last = -1; it != lst.end() && it->pos < beg + len; ++it)
            if (it->pos != last) target[it->pos - beg] = it->alt, last = it->pos;
    }
    ambiguity_window_t window{target.data(), rb};
    bns_iter_ambi(ref.bns, rb, len, mask_ambiguity, &window);
    if (mo.m.meth_context == METH_CTX_CPG && len)
        meth_label_cpg(target.data(), len, rb, rid, ref.bns, ref.pac);
}

aln_meth_eval_status aln_meth_realign_window(const aln_meth_ref_t &ref, const aln_meth_opt_t &mo,
    const aln_meth_candidate_t &c, const ubyte_t *query, int len, int n_ins, int n_del,
    aln_meth_alignment_t *out, int alt, std::vector<uint8_t> *window)
{
    if (out) *out = aln_meth_alignment_t();
    if (!out || !query || !ref.bns || !ref.pac || !ref.bns->anns || len <= 0 || len > 65535 ||
        c.rid < 0 || c.rid >= ref.bns->n_seqs || c.strand < 0 || c.strand > 1 ||
        c.hyp < 0 || c.hyp > 1 || n_ins < 0 || n_del < 0 || c.rb < 0 || c.re <= c.rb ||
        c.re > ref.bns->l_pac)
        return error(aln_meth_eval_status::invalid_input, "invalid query, reference or window");
    const bntann1_t &ann = ref.bns->anns[c.rid];
    if (ann.offset < 0 || ann.len < 0 || c.rb < ann.offset || c.re - ann.offset > ann.len || c.re - c.rb > 65535)
        return error(aln_meth_eval_status::invalid_input, "window crosses a contig or exceeds 65535 bases");
    for (int i = 0; i < len; ++i)
        if (query[i] > 4)
            return error(aln_meth_eval_status::invalid_input, "query code is outside the DNA alphabet");
    const int span = (int)(c.re - c.rb);
    const int m = mo.m.meth_context == METH_CTX_CPG ? METH_M : 5;
    const int8_t *mat = m == 5 ? (c.hyp ? mo.mat_ot : mo.mat_ob) : (c.hyp ? mo.mat7_ot : mo.mat7_ob);
    // Enforce the fixed model's score bound before entering existing int32 kernels.
    for (int i = 0; i < m * m; ++i)
        if (mat[i] < -3 || mat[i] > 1)
            return error(aln_meth_eval_status::invalid_input, "matrix is not the initialized match-1/mismatch-3 model");
    try {
        std::vector<uint8_t> local;
        std::vector<uint8_t> &target = window ? *window : local;
        aln_meth_fetch_window(ref, mo, c.rb, span, c.rid, target, alt);
        const int full = std::max(len, span), delta = std::abs(len - span);
        int band = std::min<int64_t>(full, std::max<int64_t>((int64_t)delta + 3, std::max(n_ins, n_del)));
        const int certificate = m == 5 ? bwa_global_cert_band(mat, gap_open, gap_extend, gap_open, gap_extend,
                                                              len, query, span, target.data()) : -1;
        band = certificate < 0 ? full : std::min(full, std::max(band, certificate));
        // Inputs <=65535 make these products fit uint64_t. This bounds both scalar and wave
        // active arrays conservatively; backend TLS retained from earlier calls is not this budget.
        const uint64_t longest = full, columns = std::min<int64_t>(len, 2LL * band + 1);
        const uint64_t bytes = columns * ((uint64_t)len + span) + 64 + 128 * (longest + 128) + 8 * ((uint64_t)len + span);
        if (bytes > (64ULL << 20))
            return error(aln_meth_eval_status::resource_error, "certified/full window exceeds the 64 MiB active-allocation budget");
        aln_meth_alignment_t result;
        result.rb = c.rb;
        result.re = c.re;
        result.band = band;
        int ungapped = 0;
        if (span == len && ungapped_is_optimal(mat, m, query, target.data(), len, &ungapped)) {
            result.reward = ungapped;
            result.raw_cigar.assign(1, (uint32_t)len << 4);
        } else {
            int n = 0;
            uint32_t *raw = nullptr;
            result.reward = ksw_global(len, query, span, target.data(), m, mat, gap_open, gap_extend, band, &n, &raw);
            std::unique_ptr<uint32_t, decltype(&std::free)> owner(raw, &std::free);
            if (n <= 0 || !raw)
                return error(aln_meth_eval_status::no_alignment, "global alignment returned no CIGAR");
            result.raw_cigar.assign(raw, raw + n);
        }
        result.cigar = result.raw_cigar;
        // Keep stock display normalization separate from the scored raw global alignment.
        if ((result.cigar.back() & 15) == 1) result.cigar.back() = (result.cigar.back() & ~15u) | 4;
        if ((result.cigar.front() & 15) == 1) result.cigar.front() = (result.cigar.front() & ~15u) | 4;
        if ((result.cigar.back() & 15) == 2) {
            result.re -= result.cigar.back() >> 4;
            result.cigar.pop_back();
        }
        if (!result.cigar.empty() && (result.cigar.front() & 15) == 2) {
            result.rb += result.cigar.front() >> 4;
            result.cigar.erase(result.cigar.begin());
        }
        if (result.rb >= result.re || result.cigar.empty())
            return error(aln_meth_eval_status::no_alignment, "terminal normalization left no reference alignment");
        *out = std::move(result);
        return aln_meth_eval_status::ok;
    } catch (const std::bad_alloc &) {
        return error(aln_meth_eval_status::resource_error, "cannot allocate candidate alignment buffers");
    }
}

void aln_meth_fast_init(const aln_meth_opt_t &mo, aln_meth_fast_t *f)
{
    *f = aln_meth_fast_t();
    if (mo.m.meth_context == METH_CTX_CPG) return;
    const int8_t *ot = mo.mat_ot, *ob = mo.mat_ob;
    const int8_t conv = ot[1 * 5 + 3];
    const bool collapsed = ot[3 * 5 + 1] == conv && ot[3 * 5 + 1] != ot[0 * 5 + 1];
    int best = 0;
    for (int h = 0; h < 2; ++h) {
        const int8_t *mat = h ? ot : ob;
        for (int r = 0; r < 5; ++r) for (int q = 0; q < 5; ++q) {
            const int v = mat[r * 5 + q];
            if (v < -3 || v > 1) return; // aln_meth_realign_window's model check: leave it the error
            best = std::max(best, v);
            const bool is_conv = h ? (r == 1 && q == 3) || (collapsed && r == 3 && q == 1)
                                   : (r == 2 && q == 0) || (collapsed && r == 0 && q == 2);
            const int want = (r == 4 || q == 4) ? mat[4 * 5 + 0] : is_conv ? conv : r == q ? mat[0] : mat[0 * 5 + 1];
            if (v != want) return;
        }
    }
    if (ob[2 * 5 + 0] != conv || conv < 0 || ot[0] < 0 || ot[1] >= 0 || ot[4 * 5 + 0] >= 0) return;
    f->collapsed = collapsed;
    f->match = ot[0]; f->mism = ot[1]; f->amb = ot[4 * 5 + 0]; f->conv = conv; f->best_cell = best;
    f->ok = true;
}

void aln_meth_read_bits(const ubyte_t *query, int len, aln_meth_read_bits_t *out)
{
    const size_t nw = ((size_t)len + 31) / 32;
    out->len = len;
    for (auto *v : {&out->code, &out->n, &out->a, &out->c, &out->g, &out->t}) v->assign(nw, 0);
    for (int k = 0; k < len; ++k) {
        const int w = k >> 5, sh = 62 - 2 * (k & 31);
        const int b = query[k];
        if (b > 3) { out->n[w] |= 1ULL << sh; continue; }
        out->code[w] |= (uint64_t)b << sh;
        (b == 0 ? out->a : b == 1 ? out->c : b == 2 ? out->g : out->t)[w] |= 1ULL << sh;
    }
}

namespace {
constexpr uint64_t LANE_LO = 0x5555555555555555ULL;
// 32 reference bases from forward position q, first base in the top bits (lanes past the
// reference end are 0).
inline uint64_t pac_word(const uint8_t *pac, int64_t l_pac, int64_t q)
{
    const int64_t b = q >> 2;
    const int sh = (int)(q & 3) * 2;
    if (b + 9 <= l_pac / 4 + 1) { // the 9 bytes read lie inside the pac allocation (l_pac/4 + 1)
        uint64_t v;
        memcpy(&v, pac + b, 8);
        v = __builtin_bswap64(v);
        if (sh) v = v << sh | (uint64_t)pac[b + 8] >> (8 - sh);
        return v;
    }
    uint64_t v = 0;
    for (int j = 0; j < 32 && q + j < l_pac; ++j) {
        const int64_t p = q + j;
        v |= (uint64_t)(pac[p >> 2] >> ((~p & 3) << 1) & 3) << (62 - 2 * j);
    }
    return v;
}
struct amb_bits_t { uint64_t *bits; int64_t rb; };
int amb_visit(int64_t begin, int64_t end, void *context)
{
    auto *a = static_cast<amb_bits_t *>(context);
    for (int64_t p = begin; p < end; ++p) {
        const int64_t k = p - a->rb;
        a->bits[k >> 5] |= 1ULL << (62 - 2 * (k & 31));
    }
    return 0;
}
} // namespace

namespace {
// 32 lanes of a plane starting at base index b (lanes past the plane's end are 0)
inline uint64_t plane_at(const uint64_t *p, size_t nw, int64_t b)
{
    const int64_t a = b >> 5;
    const int s = (int)(b & 31) * 2;
    const uint64_t hi = (size_t)a < nw ? p[a] : 0;
    if (!s) return hi;
    const uint64_t lo = (size_t)(a + 1) < nw ? p[a + 1] : 0;
    return hi << s | lo >> (64 - s);
}
} // namespace

void aln_meth_amb_bits(const aln_meth_ref_t &ref, int64_t rb, int len, std::vector<uint64_t> *bits)
{
    bits->assign(((size_t)len + 31) / 32, 0);
    amb_bits_t ctx{bits->data(), rb};
    bns_iter_ambi(ref.bns, rb, len, amb_visit, &ctx);
}

void aln_meth_score_block(const aln_meth_ref_t &ref, const aln_meth_fast_t &f, int hyp, int64_t rpos, int qs, int n,
                          const aln_meth_read_bits_t &rd, const uint64_t *amb, int wo, aln_meth_block_counts_t *out)
{
    const size_t nq = rd.code.size(), na = ((size_t)wo + n + 31) / 32;
    const bool seed = !rd.seed.empty();
    for (int w = 0; 32 * w < n; ++w) {
        const int64_t q = qs + 32LL * w;
        const uint64_t R = pac_word(ref.pac, ref.bns->l_pac, rpos + 32LL * w), Q = plane_at(rd.code.data(), nq, q);
        const int rem = n - 32 * w;
        const uint64_t valid = rem >= 32 ? LANE_LO : LANE_LO & (~0ULL << (64 - 2 * rem));
        const uint64_t x = R ^ Q, neq = (x | x >> 1) & LANE_LO;
        const uint64_t lo = R & LANE_LO, hi = R >> 1 & LANE_LO;
        uint64_t conv;
        if (hyp) {
            conv = lo & ~hi & plane_at(rd.t.data(), nq, q);
            if (f.collapsed) conv |= lo & hi & plane_at(rd.c.data(), nq, q);
        } else {
            conv = ~lo & hi & plane_at(rd.a.data(), nq, q);
            if (f.collapsed) conv |= ~lo & ~hi & LANE_LO & plane_at(rd.g.data(), nq, q);
        }
        const uint64_t nn = (plane_at(rd.n.data(), nq, q) | plane_at(amb, na, wo + 32LL * w)) & valid;
        const uint64_t mm = ((neq & ~conv) & valid) | nn;
        out->n_amb += __builtin_popcountll(nn);
        out->n_mm += __builtin_popcountll(mm);
        out->n_conv += __builtin_popcountll(conv & ~nn & valid);
        if (seed) out->n_seed_mm += __builtin_popcountll(mm & plane_at(rd.seed.data(), rd.seed.size(), q));
    }
}

void aln_meth_score_ungapped(const aln_meth_ref_t &ref, const aln_meth_fast_t &f, int hyp, int64_t rb,
                             const aln_meth_read_bits_t &rd, aln_meth_ungapped_t *out)
{
    const int len = rd.len, nw = (len + 31) / 32;
    uint64_t amb_small[8];
    std::vector<uint64_t> amb_big;
    uint64_t *amb = amb_small;
    if (nw > 8) { amb_big.assign(nw, 0); amb = amb_big.data(); }
    else std::fill(amb_small, amb_small + 8, 0);
    amb_bits_t ctx{amb, rb};
    bns_iter_ambi(ref.bns, rb, len, amb_visit, &ctx);
    int n_mm = 0, n_conv = 0, n_amb = 0;
    for (int w = 0; w < nw; ++w) {
        const uint64_t R = pac_word(ref.pac, ref.bns->l_pac, rb + 32LL * w), Q = rd.code[w];
        const int rem = len - 32 * w;
        const uint64_t valid = rem >= 32 ? LANE_LO : LANE_LO & (~0ULL << (64 - 2 * rem));
        const uint64_t x = R ^ Q, neq = (x | x >> 1) & LANE_LO;
        const uint64_t lo = R & LANE_LO, hi = R >> 1 & LANE_LO;
        uint64_t conv;
        if (hyp) { // OT: reference C / read T (collapsed: and reference T / read C)
            conv = lo & ~hi & rd.t[w];
            if (f.collapsed) conv |= lo & hi & rd.c[w];
        } else {   // OB: reference G / read A (collapsed: and reference A / read G)
            conv = ~lo & hi & rd.a[w];
            if (f.collapsed) conv |= ~lo & ~hi & LANE_LO & rd.g[w];
        }
        const uint64_t nn = (rd.n[w] | amb[w]) & valid;
        n_amb += __builtin_popcountll(nn);
        n_mm += __builtin_popcountll(((neq & ~conv) & valid) | nn);
        n_conv += __builtin_popcountll(conv & ~nn & valid);
    }
    out->n_mm = n_mm; out->n_conv = n_conv;
    out->reward = (len - n_mm - n_conv) * f.match + n_conv * f.conv + (n_mm - n_amb) * f.mism + n_amb * f.amb;
}

} // namespace bwa_aln
