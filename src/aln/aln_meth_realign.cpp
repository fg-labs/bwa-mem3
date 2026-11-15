/* SPDX-License-Identifier: MIT */
#include "aln_meth_realign.h"
#include "../bwa.h"
#include "../ksw.h"
#include "../meth_model.h"

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <memory>
#include <new>

namespace bwa_aln {
namespace {

aln_meth_eval_status error(aln_meth_eval_status status, const char *message)
{
    std::fprintf(stderr, "ERROR: --meth candidate realignment: %s\n", message);
    return status;
}

} // namespace

aln_meth_eval_status aln_meth_realign_window(const aln_meth_ref_t &ref, const aln_meth_opt_t &mo,
    const aln_meth_candidate_t &c, const ubyte_t *query, int len, int n_ins, int n_del,
    aln_meth_alignment_t *out)
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
        std::vector<uint8_t> target(span);
        for (int i = 0; i < span; ++i) {
            const int64_t p = c.rb + i;
            target[i] = ref.pac[p >> 2] >> ((~p & 3) << 1) & 3;
        }
        if (m == METH_M) meth_label_cpg(target.data(), span, c.rb, c.rid, ref.bns, ref.pac);
        const int full = std::max(len, span), delta = std::abs(len - span);
        int band = std::min<int64_t>(full, std::max<int64_t>((int64_t)delta + 3, std::max(n_ins, n_del)));
        const int certificate = m == 5 ? bwa_global_cert_band(mat, 5, 1, 5, 1, len, query, span, target.data()) : -1;
        band = certificate < 0 ? full : std::min(full, std::max(band, certificate));
        // Inputs <=65535 make these products fit uint64_t. This bounds both scalar and wave
        // active arrays conservatively; backend TLS retained from earlier calls is not this budget.
        const uint64_t longest = full, columns = std::min<int64_t>(len, 2LL * band + 1);
        const uint64_t bytes = columns * ((uint64_t)len + span) + 64 + 128 * (longest + 128) + 8 * ((uint64_t)len + span);
        if (bytes > (64ULL << 20))
            return error(aln_meth_eval_status::resource_error, "certified/full window exceeds the 64 MiB active-allocation budget");
        int n = 0;
        uint32_t *raw = nullptr;
        const int reward = ksw_global(len, query, span, target.data(), m, mat, 5, 1, band, &n, &raw);
        std::unique_ptr<uint32_t, decltype(&std::free)> owner(raw, &std::free);
        if (n <= 0 || !raw)
            return error(aln_meth_eval_status::no_alignment, "global alignment returned no CIGAR");
        aln_meth_alignment_t result;
        result.reward = reward;
        result.rb = c.rb;
        result.re = c.re;
        result.band = band;
        result.raw_cigar.assign(raw, raw + n);
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

} // namespace bwa_aln
