#ifndef BWA_HUGEPAGES_H
#define BWA_HUGEPAGES_H

#include <stddef.h>

/* Opt-in 1 GB huge pages for the FM-index / suffix-array structures (issue
 * #402, measured in #377). Seeding is memory-latency-bound: the SA/BWT lookups
 * are random accesses over multi-gigabyte arrays, so backing them with 1 GB
 * pages cuts data-TLB misses. bwa-mem3's allocator is the vendored mimalloc,
 * which can reserve 1 GB huge pages at startup; this is the programmatic,
 * auto-sized, safe-by-default form of `MIMALLOC_RESERVE_HUGE_OS_PAGES`. */

#ifdef __cplusplus
extern "C" {
#endif

/* Number of 1 GB huge pages to reserve to cover an index whose two large
 * resident arrays total `bwt_bytes` (the FM-index, `<prefix>.bwt.2bit.64`) and
 * `pac_bytes` (the packed reference, `<prefix>.pac`). Rounds the total up to a
 * whole number of 1 GB pages and adds a small margin for allocator metadata and
 * rounding. Returns 0 when `bwt_bytes <= 0` (size unknown), so the caller
 * declines to reserve. Pure; unit-tested. */
size_t bwamem_huge_pages_needed(long long bwt_bytes, long long pac_bytes);

/* Outcome of a 1 GB huge-page reservation attempt. */
typedef enum {
    BWAMEM_HUGE_FULL,       /* all `need` pages reserved */
    BWAMEM_HUGE_PARTIAL,    /* mimalloc reported success but fewer than `need` landed */
    BWAMEM_HUGE_FAILED,     /* mimalloc reported an error */
    BWAMEM_HUGE_UNVERIFIED  /* mimalloc reported success (`rc == 0`), but the pool
                             * delta could not be measured (a free-count read
                             * failed, or the host-global count did not drop —
                             * e.g. perturbed by another process). The count is
                             * unknown, so we must NOT claim a full reservation. */
} bwamem_huge_outcome_t;

/* Classify a reservation from mimalloc's return code `rc` and the hugetlb
 * free-page count before/after the call. mimalloc's `rc == 0` only guarantees at
 * least one page per NUMA node, so the true count reserved is
 * `free_before - free_after`; a value below `need` is a PARTIAL reservation that
 * must not be reported as a clean success. `*reserved_out` (when non-NULL) gets
 * that measured count, or -1 when it cannot be measured (a free-count read
 * failed, or the host-global count did not drop). The `free_hugepages` count is
 * host-global, so when the delta is unmeasurable the true per-call count is
 * unknown: `rc == 0` in that case is classified UNVERIFIED, never FULL. Pure;
 * unit-tested. */
bwamem_huge_outcome_t bwamem_classify_huge_reservation(
    int rc, long free_before, long free_after, size_t need, long *reserved_out);

/* Opt-in (`--huge-pages`): when the Linux host has enough *free* 1 GB huge pages
 * to cover the index at `idxbase`, reserve them through mimalloc so the
 * FM-index / suffix-array arrays land on 1 GB pages. Safe by default: a no-op —
 * with a one-line `[M::]` note — when the pool is missing or too small, on
 * non-Linux, or on a non-mimalloc build. MUST be called before the index is
 * loaded, since mimalloc only serves later allocations from the reserved arena. */
void bwamem_reserve_huge_pages(const char *idxbase);

#ifdef __cplusplus
}
#endif

#endif /* BWA_HUGEPAGES_H */
