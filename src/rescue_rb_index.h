/* rescue_rb_index.h -- a sorted copy of a mate list's rb values, which answers a rescue anchor's
 * skip set with four binary searches (mem_matesw_skip_probe) instead of a scan of the list.
 * Overview and gates: docs/src/developer-guide/rescue.md.
 *
 * An anchor's rescue skip set depends only on the set of the mate list's rb values (see
 * mem_matesw_skip_probe: duplicates and order do not matter). Two users:
 *
 *   - mem_sam_pe_batch_pre builds the copy once per mate list that many anchors rescue into, and
 *     probes it for each anchor; the list does not change during _pre.
 *   - mem_pair_resolve_batch_post builds it after the first mem_matesw_batch_post call that may
 *     have changed the list, probes it for each later anchor, and after each later call either
 *     adds the rb values the call pushed (mem_rescue_pushed_rb_t) or discards the copy:
 *       - a call never rewrites the rb of a region already in the list (its dedups pass
 *         bns == NULL, so mem_patch_reg merges nothing), and its dedups only drop and reorder
 *         regions;
 *       - so the list after the call holds the old rb values, plus the pushed ones, minus any
 *         dropped;
 *       - n_after == n_before + pushed_rb.n therefore means nothing was dropped, and the copy plus
 *         the pushed values is exactly the new list's rb values. Any other count discards the copy.
 *
 * The BWA_MEM3_DEBUG_RESCUE_SKIP_XCHECK build checks the _post copy against the list after every
 * call (matches), and every probed skip set against the per-region test (matesw_skip_bits). */
#ifndef BWA_MEM3_RESCUE_RB_INDEX_H
#define BWA_MEM3_RESCUE_RB_INDEX_H

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <vector>

#include "bwamem.h" /* mem_alnreg_t, mem_rescue_pushed_rb_t, mem_matesw_skip_query_t, mem_matesw_skip_probe */

struct rescue_rb_index {
    std::vector<int64_t> rb; /* ascending; meaningful only while live */
    bool live = false;

    /* Copy and sort the rb values of a[0, n). */
    void build(const mem_alnreg_t *a, size_t n)
    {
        rb.resize(n);
        for (size_t k = 0; k < n; ++k) rb[k] = a[k].rb;
        std::sort(rb.begin(), rb.end());
        live = true;
    }

    void drop() { live = false; }

    /* Bring a live copy up to date after one mem_matesw_batch_post call that took the list from
     * n_before to n_after regions and reported `pushed_rb`. Returns whether the copy is still live. */
    bool update(size_t n_before, size_t n_after, const mem_rescue_pushed_rb_t &pushed_rb)
    {
        if (!live) return false;
        if (n_after != n_before + (size_t)pushed_rb.n) {
            live = false;
            return false;
        }
        for (int k = 0; k < pushed_rb.n; ++k)
            rb.insert(std::upper_bound(rb.begin(), rb.end(), pushed_rb.rb[k]), pushed_rb.rb[k]);
        return true;
    }

    /* The skip set of query q against the copy (call only while live). */
    unsigned probe(const mem_matesw_skip_query_t *q) const
    {
        return mem_matesw_skip_probe(q, rb.data(), rb.size());
    }

    /* Whether the copy holds exactly the rb values of a[0, n) (the XCHECK test). */
    bool matches(const mem_alnreg_t *a, size_t n) const
    {
        std::vector<int64_t> x(n);
        for (size_t k = 0; k < n; ++k) x[k] = a[k].rb;
        std::sort(x.begin(), x.end());
        return x == rb;
    }
};

#endif /* BWA_MEM3_RESCUE_RB_INDEX_H */
