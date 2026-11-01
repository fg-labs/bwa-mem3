#ifndef BWA_HUGEPAGES_H
#define BWA_HUGEPAGES_H

#include <stddef.h>

/* Opt-in 1 GB huge pages for the FM-index / suffix-array structures (issue
 * #402, measured in #377). Seeding is memory-latency-bound: the SA/BWT lookups
 * are random accesses over multi-gigabyte arrays, so 1 GB pages are intended to
 * cut data-TLB misses. bwa-mem3's allocator is the vendored mimalloc, which can
 * reserve 1 GB huge pages at startup; this is the programmatic, auto-sized,
 * safe-by-default form of `MIMALLOC_RESERVE_HUGE_OS_PAGES`.
 *
 * Known limitation: mimalloc serves allocations above its arena object limit
 * (2 GiB) from the OS, not from the reserved arena. The two largest index
 * arrays (the occurrence table and the suffix-array word array on hg38) are
 * above it, so today the reserved pages hold only the smaller allocations and
 * the rest of the reservation sits unused. See "Memory allocator -> Measured
 * effect" in the docs. */

#ifdef __cplusplus
extern "C" {
#endif

/* Number of 1 GB huge pages to reserve to cover an index whose two large
 * resident arrays total `bwt_bytes` (the FM-index, `<prefix>.bwt.2bit.64`) and
 * `pac_bytes` (the packed reference, `<prefix>.pac`). Rounds the total up to a
 * whole number of 1 GB pages and adds a small margin for allocator metadata and
 * rounding. Returns 0 when `bwt_bytes <= 0` (size unknown), so the caller
 * declines to reserve. The count is sized for the whole index on purpose: it
 * does not model which allocations mimalloc can place in its arena (see the
 * known limitation above), so pages for the arrays over 2 GiB are reserved and
 * held for the life of the process but stay unused for now. Pure; unit-tested. */
size_t bwamem_huge_pages_needed(long long bwt_bytes, long long pac_bytes);

/* 1 GB pages a new reservation can still take, from the pool's `free_hugepages`
 * and `resv_hugepages` counters. A hugetlb mapping reserves its pages when it is
 * created but takes them off `free_hugepages` only as they are first touched, so
 * the free count alone overstates what is available. Right after mimalloc's
 * reservation it has dropped only by the pages mimalloc touched (one per NUMA
 * node, for its arena metadata), not by the pages it reserved. Returns
 * `free - resv` (floored at 0), `free` when the reserve count is unreadable
 * (`resv < 0`), or -1 when the free count is unreadable (no 1 GB pool). Pure;
 * unit-tested. */
long bwamem_huge_pages_available(long free_pages, long resv_pages);

/* Read one counter (`name`, e.g. "free_hugepages") from the hugetlb pool
 * directory `dir`, or -1 if it is missing or unparsable. */
long bwamem_read_hugepage_counter(const char *dir, const char *name);

/* Available pages (bwamem_huge_pages_available) of the pool directory `dir`,
 * reading its `free_hugepages` and `resv_hugepages` counters; `resv_hugepages`
 * is read first, then `free_hugepages`, so a page that a concurrent fault moves
 * from reserved to used between the two reads can only make the result one page
 * low, never high. -1 if `free_hugepages` is missing. `dir` is a parameter so a
 * unit test can point it at a temporary directory. */
long bwamem_avail_hugepages_at(const char *dir);

/* Outcome of a 1 GB huge-page reservation attempt. */
typedef enum {
    BWAMEM_HUGE_FULL,       /* all `need` pages reserved */
    BWAMEM_HUGE_PARTIAL,    /* mimalloc reported success but fewer than `need` landed */
    BWAMEM_HUGE_FAILED,     /* mimalloc reported an error */
    BWAMEM_HUGE_UNVERIFIED  /* mimalloc reported success (`rc == 0`), but the pool
                             * delta could not be measured (an available-count
                             * read failed, or the host-global count rose, e.g.
                             * perturbed by another process). The count is
                             * unknown, so we must NOT claim a full reservation. */
} bwamem_huge_outcome_t;

/* Classify a reservation from mimalloc's return code `rc` and the hugetlb
 * available-page count (bwamem_huge_pages_available) before and after the call.
 * mimalloc's `rc == 0` only guarantees at least one page per NUMA node, so the
 * count reserved is `avail_before - avail_after`; a value below `need` is a
 * PARTIAL reservation that must not be reported as a clean success.
 * `*reserved_out` (when non-NULL) gets that count, or -1 when it cannot be
 * measured (an available-count read failed, or the count rose). The counters
 * are host-global, so when the delta is unmeasurable the per-call count is
 * unknown: `rc == 0` in that case is classified UNVERIFIED, never FULL. A
 * concurrent process that reserves or frees 1 GB pages between the two reads
 * can still make the delta wrong in either direction without the count rising,
 * so FULL means "the pool gave up `need` pages", not that every one is ours.
 * Pure; unit-tested. */
bwamem_huge_outcome_t bwamem_classify_huge_reservation(
    int rc, long avail_before, long avail_after, size_t need, long *reserved_out);

/* Opt-in (`--huge-pages`): when the Linux host has enough available 1 GB huge
 * pages (free minus reserved) to cover the index at `idxbase`, reserve them
 * through mimalloc before the index is allocated. Allocations above mimalloc's
 * 2 GiB arena object limit are not served from them (see the known limitation
 * above). Safe by default: a no-op, with a one-line `[M::]` note, when the pool
 * is missing or too small, when `MIMALLOC_RESERVE_HUGE_OS_PAGES` is already
 * set, on non-Linux, or on a non-mimalloc build. MUST be called before the
 * index is loaded, since mimalloc only serves later allocations from the
 * reserved arena. */
void bwamem_reserve_huge_pages(const char *idxbase);

#ifdef __cplusplus
}
#endif

#endif /* BWA_HUGEPAGES_H */
