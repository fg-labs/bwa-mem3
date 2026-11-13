/* See aln_bwt.h. bwt_match_exact_alt is from bwa's Apache2 branch (bwt.c). This file
 * is MODIFIED from the original: the function body is unchanged but moved into
 * namespace bwa_aln, and the occ/SA functions (bwt_sa_batch included) are new code.
 *
 * The MIT License
 *
 * Copyright (c) 2008 Genome Research Ltd (GRL).
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject to
 * the following conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
 * BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
 * ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
 * CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */
#include "aln_bwt.h"
#include "aln_alloc.h"

#include "../FMI_search.h"

#include <stdlib.h>

namespace bwa_aln {

bwt_t *bwt_from_fmi(FMI_search *fmi)
{
    bwt_t *bwt = (bwt_t *)aln_calloc(1, sizeof(bwt_t));
    const int64_t *count = fmi->count_data();
    for (int i = 0; i < 5; ++i) bwt->L2[i] = (bwtint_t)count[i] - 1; // count[] includes the `$`
    bwt->seq_len = (bwtint_t)fmi->reference_seq_len - 1;
    bwt->primary = (bwtint_t)fmi->sentinel_index;
    bwt->fmi = fmi;
    return bwt;
}

void bwt_destroy(bwt_t *bwt)
{
    free(bwt);
}

/* Count of c in bwa-mem3 BWT rows [0, p). The `$` row has no base bit set, so
 * it adds nothing, which is exactly bwa's `$`-removed convention. */
static inline bwtint_t occ_rows_before(const bwt_t *bwt, bwtint_t p, ubyte_t c)
{
    const CP_OCC *cp = (const CP_OCC *)bwt->fmi->cp_occ_data();
    const CP_OCC &blk = cp[p >> CP_SHIFT];
    const unsigned r = (unsigned)(p & CP_MASK);
    const uint64_t mask = r ? ~0ULL << (CP_BLOCK_SIZE - r) : 0;
    return (bwtint_t)blk.cp_count[c] + (bwtint_t)__builtin_popcountll(blk.one_hot_bwt_str[c] & mask);
}

bwtint_t bwt_occ(const bwt_t *bwt, bwtint_t k, ubyte_t c)
{
    if (k == (bwtint_t)(-1)) return 0;
    return occ_rows_before(bwt, k + 1, c);
}

void bwt_2occ(const bwt_t *bwt, bwtint_t k, bwtint_t l, ubyte_t c, bwtint_t *ok, bwtint_t *ol)
{
    *ok = bwt_occ(bwt, k, c);
    *ol = bwt_occ(bwt, l, c);
}

void bwt_2occ4(const bwt_t *bwt, bwtint_t k, bwtint_t l, bwtint_t cntk[4], bwtint_t cntl[4])
{
    for (int c = 0; c < 4; ++c) {
        cntk[c] = bwt_occ(bwt, k, (ubyte_t)c);
        cntl[c] = bwt_occ(bwt, l, (ubyte_t)c);
    }
}

void bwt_prefetch_occ(const bwt_t *bwt, bwtint_t p)
{
    const CP_OCC *cp = (const CP_OCC *)bwt->fmi->cp_occ_data();
    __builtin_prefetch(&cp[p >> CP_SHIFT]);
}

bwtint_t bwt_sa(const bwt_t *bwt, bwtint_t k)
{
    if (k == 0) return (bwtint_t)(-1); // bwa stores -1 for the `$` suffix; bwa-mem3 stores seq_len
    int64_t pos = (int64_t)k, sa_entry = 0, offset = 0;
    // drop_sentinel_offset = 0: bwa's walk reaching the `$` row keeps its offset. This is the walk
    // FMI_search::get_sa_entry_compressed does (same values); call_one_step is used for its
    // branchless NEON symbol select, which the scalar get_sa_entry_compressed chain lacks.
    while (!bwt->fmi->call_one_step(pos, sa_entry, offset, 0)) pos = sa_entry;
    return (bwtint_t)sa_entry;
}

/* Lanes walked at once: mem's resolver's count, SA_RESOLVE_LANES_DEFAULT (FMI_search.h). It only
 * sets how many fetches are in flight; each output depends on its row alone. */

void bwt_sa_batch(const bwt_t *bwt, const bwtint_t *rows, int64_t n, bwtint_t *out)
{
    FMI_search *fmi = bwt->fmi;
    const CP_OCC *cp = (const CP_OCC *)fmi->cp_occ_data();
    const int8_t *ms = fmi->sa_ms_byte_data();
    const uint32_t *ls = fmi->sa_ls_word_data();
    const int64_t compx = fmi->sa_compx, compx_mask = fmi->sa_compx_mask;
    // What the next call_one_step on row p reads: the SA sample for a sampled row, else
    // the row's Occ block (the selection mem's resolver makes).
    auto prefetch_row = [&](int64_t p) {
#ifdef ENABLE_PREFETCH
        if ((p & compx_mask) == 0) {
            __builtin_prefetch(&ms[p >> compx]);
            __builtin_prefetch(&ls[p >> compx]);
        } else {
            __builtin_prefetch(&cp[p >> CP_SHIFT]);
        }
#else
        (void)p, (void)cp, (void)ms, (void)ls, (void)compx, (void)compx_mask;
#endif
    };
    int64_t row[SA_RESOLVE_LANES_DEFAULT], off[SA_RESOLVE_LANES_DEFAULT], dst[SA_RESOLVE_LANES_DEFAULT];
    int64_t next = 0, n_live = 0;
    // Stages the next row needing a walk into lane j; row 0 is bwt_sa's special case.
    // A row is read before its own output is written, so out may alias rows.
    auto stage = [&](int j) -> bool {
        while (next < n) {
            const int64_t i = next++;
            const bwtint_t k = rows[i];
            if (k == 0) {
                out[i] = (bwtint_t)(-1);
                continue;
            }
            row[j] = (int64_t)k, off[j] = 0, dst[j] = i;
            prefetch_row((int64_t)k);
            return true;
        }
        return false;
    };
    int lim = 0;
    while (lim < SA_RESOLVE_LANES_DEFAULT && stage(lim)) ++lim;
    n_live = lim;
    while (n_live > 0) {
        for (int j = 0; j < lim; ++j) {
            if (off[j] < 0) continue;
            int64_t sa_entry = 0;
            // drop_sentinel_offset = 0, as bwt_sa
            if (fmi->call_one_step(row[j], sa_entry, off[j], 0)) {
                out[dst[j]] = (bwtint_t)sa_entry;
                if (!stage(j)) off[j] = -1, --n_live;
            } else {
                row[j] = sa_entry;
                prefetch_row(sa_entry);
            }
        }
    }
}

int bwt_match_exact_alt(const bwt_t *bwt, int len, const ubyte_t *str, bwtint_t *k0, bwtint_t *l0)
{
    int i;
    bwtint_t k, l, ok, ol;
    k = *k0; l = *l0;
    for (i = len - 1; i >= 0; --i) {
        ubyte_t c = str[i];
        if (c > 3) return 0; // there is an N here. no match
        bwt_2occ(bwt, k - 1, l, c, &ok, &ol);
        k = bwt->L2[c] + ok + 1;
        l = bwt->L2[c] + ol;
        if (k > l) return 0; // no match
    }
    *k0 = k; *l0 = l;
    return l - k + 1;
}

} // namespace bwa_aln
