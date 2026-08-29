#include "smem_dedup.h"
#include <string.h>  /* memcpy */

/* Compact adjacent fully-identical SMEMs in a sorted array; return the new length.
 *
 * Write-pointer compaction: keep a[0], then for each subsequent entry compare
 * to the last *kept* entry on (rid, m, n, k, s). If different, copy to the write
 * position and advance the write pointer. O(n), no allocation.
 *
 * The reverse-complement interval start `l` is deliberately NOT compared: an
 * SMEM's full FM interval (k, l, s) is a deterministic function of the matched
 * substring read[m..n], so (rid, m, n) => (k, l, s). Two entries equal on
 * (rid, m, n) are equal on `l`; entries differing on (rid, m, n) already differ
 * before `l` would be reached — so comparing `l` never changes the result. This
 * also lets the pure-backward SMEM phase skip maintaining the (dead) l-chain
 * (see backwardExt_konly in FMI_search.h).
 *
 * Precondition: identical SMEMs are adjacent in the array (guaranteed after
 * sortSMEMs + per-read ks_introsort(mem_intv1) in mem_collect_smem). */
int64_t smem_dedup_inplace(SMEM *a, int64_t n)
{
    if (n <= 1) return n;

    int64_t w = 0;  /* write pointer — a[w] is the last kept entry */
    for (int64_t i = 1; i < n; ++i) {
        const SMEM *prev = &a[w];
        const SMEM *cur  = &a[i];
        if (cur->rid != prev->rid ||
            cur->m   != prev->m   ||
            cur->n   != prev->n   ||
            cur->k   != prev->k   ||
            cur->s   != prev->s)
        {
            ++w;
            if (w != i)
                memcpy(&a[w], cur, sizeof(SMEM));
        }
        /* else: duplicate — skip */
    }
    return w + 1;
}
