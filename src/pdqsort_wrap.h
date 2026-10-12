/* Wrapper: emit a pdqsort-backed function with the same signature klib's
 * ksort.h emits for ks_introsort. Lets us drop pdqsort into call sites
 * site-by-site (pdqsort_##name(n, a) replaces ks_introsort(name, n, a))
 * without touching klib or the original KSORT_INIT scaffolding.
 *
 * pdqsort (orlp/pdqsort) is the modern industry-standard pattern-
 * defeating quicksort variant (used by Rust stdlib, Boost). Pre-sorted-run
 * detection, faster than klib's introsort at every N >= 16 on the bwa-mem3
 * microbench (measured on a local microbench; the harness is not committed).
 * PDQSORT_INIT uses pdqsort's default partition: with a lambda comparator
 * pdqsort() never selects its branchless one. PDQSORT_BRANCHLESS_INIT below
 * forces the branchless (block) partition, for the sorts that can use it.
 *
 * The CMP argument is the same `__sort_lt(a, b)` macro you pass to
 * KSORT_INIT — they share a strict-weak-order contract, so the existing
 * comparators (alnreg_slt, alnreg_slt2, alnreg_hlt, alnreg_hlt2,
 * pair64_lt, ks_lt_generic, etc.) all work unchanged.
 */
#ifndef BWA_MEM3_PDQSORT_WRAP_H
#define BWA_MEM3_PDQSORT_WRAP_H

#include "../ext/pdqsort/pdqsort.h"

/* PDQSORT_INIT must be invoked from a .cpp file (not a header) since it
 * emits a non-static function definition. Declare the resulting symbol
 * in a header if it needs to be called from another translation unit
 * (declare `void pdqsort_<name>(size_t, type_t *)` in the owning header). */
#define PDQSORT_INIT_WITH_(sortfn, name, type_t, CMP) \
    void pdqsort_##name(size_t n, type_t *a) { \
        sortfn(a, a + n, \
               [](const type_t &__x, const type_t &__y) { return CMP(__x, __y); }); \
    }
#define PDQSORT_INIT(name, type_t, CMP) PDQSORT_INIT_WITH_(pdqsort, name, type_t, CMP)

/* Same, with pdqsort's branchless (block) partition. Only for sorts whose
 * result is unique -- a strict total order over the whole element, or a caller
 * that discards the result on a tie -- because the two partition schemes
 * permute equal elements differently. CMP should itself be branch-free (bitwise
 * & | over the key comparisons) or the partition's comparisons stay branches. */
#define PDQSORT_BRANCHLESS_INIT(name, type_t, CMP) \
    PDQSORT_INIT_WITH_(pdqsort_branchless, name, type_t, CMP)

#endif /* BWA_MEM3_PDQSORT_WRAP_H */
