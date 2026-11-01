#ifndef BWA_HUGEPAGES_H
#define BWA_HUGEPAGES_H

#include <stddef.h>

/* Opt-in 1 GB huge pages for the FM-index occurrence table (issue #402,
 * measured in #377). Seeding is memory-latency-bound: its Occ lookups are
 * random accesses over a multi-gigabyte table, and even on transparent 2 MB
 * pages the page walks for them take a measurable share of the cycles.
 *
 * With `--huge-pages`, FMI_search::load_index maps the occurrence table directly
 * with mmap(MAP_HUGETLB | MAP_HUGE_1GB) from the host's 1 GB hugetlb pool, and
 * faults the pages in over several threads before filling it. The occurrence
 * table's Occ lookups are where nearly all the page walks are; the suffix-array
 * samples always stay on the default (THP-advised) allocation, since on a host
 * where 1 GB pages do reach the TLB, mapping them as well removed few further
 * walks. Every step is best-effort: when the pool is missing or too small, or
 * the mapping or its pre-fault fails (as it does on a kernel before 5.14, which
 * lacks MADV_POPULATE_WRITE), the table falls back to the default allocation. Its
 * contents, and so the output, are the same either way. Linux only; elsewhere
 * every mapping fails with ENOSYS and the default allocation is used. */

#ifdef __cplusplus
extern "C" {
#endif

/* Whether FMI_search::load_index maps the occurrence table on 1 GB pages;
 * passed to load_index by its caller (main_mem passes ON for --huge-pages).
 * Production callers pass only OFF or ON. */
typedef enum {
    BWAMEM_HUGETLB_OFF = 0,        /* default allocation (no --huge-pages) */
    BWAMEM_HUGETLB_ON = 1,         /* try a 1 GB mapping, fall back on failure */

    /* Internal test seams for the unit tests only. A production caller must
     * never pass these: they bypass the size rule below or do not use 1 GB
     * pages at all. */
    BWAMEM_HUGETLB_FORCE_FAIL = 2,  /* skip the mapping, as if it failed (ECANCELED) */
    BWAMEM_HUGETLB_ON_ANY_SIZE = 3, /* ON without the minimum-size rule below */
    BWAMEM_HUGETLB_TEST_PLAIN = 4,  /* ON_ANY_SIZE, but mapped with
                                     * bwamem_hugetlb_map_plain_for_test and not
                                     * pre-faulted: exercises the mapped-table path
                                     * (fill, use, unmap) without a hugetlb pool */
    BWAMEM_HUGETLB_TEST_PREFAULT_FAIL = 5 /* ordinary PROT_NONE mapping: the real
                                          * prefault fails before loading, testing
                                          * fallback without a hugetlb pool */
} bwamem_hugetlb_mode_t;

/* An occurrence table this small or smaller is not put on 1 GB pages: rounding
 * it up to a whole page would more than double it, and a mapped page is pinned,
 * unswappable memory taken from a pool other processes share. 512 MiB keeps the
 * table of a human-sized index (gigabytes) on 1 GB pages and leaves the table of
 * a small genome or amplicon panel on the default allocation. */
#define BWAMEM_HUGETLB_MIN_BYTES ((size_t)512 * 1024 * 1024)

/* `bytes` rounded up to whole 1 GiB pages, or 0 when `bytes` is 0 or the
 * rounding would overflow. Pure; unit-tested. */
size_t bwamem_hugetlb_round_up(size_t bytes);

/* One 1 GB hugetlb mapping. `base == NULL` (and `map_bytes == 0`) means not
 * mapped. */
typedef struct {
    void  *base;
    size_t map_bytes;
} bwamem_hugetlb_buf_t;

/* Map `bytes` (rounded up to whole 1 GiB pages) of private anonymous memory on
 * 1 GB pages. Returns 0 and fills `*out`, or returns an errno value and leaves
 * `*out` unmapped. The pages are reserved from the hugetlb pool when the
 * mapping is created (no MAP_NORESERVE), so a pool that is too small fails
 * here with ENOMEM rather than with SIGBUS on first touch. Returns ENOSYS on
 * non-Linux builds, EINVAL when `bytes` is 0 or too large, and otherwise the
 * mmap errno (EINVAL also when the kernel has no 1 GB huge page size). */
int bwamem_hugetlb_map(size_t bytes, bwamem_hugetlb_buf_t *out);

/* Test seam (BWAMEM_HUGETLB_TEST_PLAIN): the same contract as
 * bwamem_hugetlb_map, including the rounding to whole 1 GiB, but an ordinary
 * private anonymous mapping (no MAP_HUGETLB), so it needs no pool and commits
 * only the pages that are touched. Lets a unit test run the mapped-table path
 * in CI. ENOSYS on non-Linux builds. */
int bwamem_hugetlb_map_plain_for_test(size_t bytes, bwamem_hugetlb_buf_t *out);

/* TEST_PREFAULT_FAIL: ordinary anonymous PROT_NONE mapping, rounded as above.
 * MADV_POPULATE_WRITE must fail on it, exercising the real prefault cleanup
 * and loader fallback without committing a 1 GiB page. ENOSYS off Linux. */
int bwamem_hugetlb_map_protected_for_test(size_t bytes, bwamem_hugetlb_buf_t *out);

/* Unmap `*b` (if mapped) and reset it to unmapped. NULL is a no-op. */
void bwamem_hugetlb_unmap(bwamem_hugetlb_buf_t *b);

/* Fault in every 1 GiB page of the mapping `*buf` before it is filled, using up
 * to `nthreads` threads, one page per task. The kernel zeroes
 * a 1 GB page in the fault that first touches it, on the faulting thread; left
 * to the index readers, several readers share one page and wait on that one
 * zeroing, which made start-up slower than on default pages. Spreading the
 * faults over threads removes that serialization, and the readers then only
 * copy. Faulting uses MADV_POPULATE_WRITE, so it needs Linux 5.14 or later; on
 * an older kernel the pre-fault fails. If any page fails, the mapping is
 * unmapped (and reset), so the caller falls back to the default allocation.
 * Returns 1 when it unmapped the mapping this way, else 0 (also for NULL or an
 * unmapped `*buf`). */
int bwamem_hugetlb_prefault(bwamem_hugetlb_buf_t *buf, int nthreads);

/* Internal test seam: same worker loop and cleanup, but the caller does not
 * claim pages. Reports spawned workers in *workers; a nonzero count plus a
 * prefault failure proves a spawned worker propagated it. Not for loaders. */
int bwamem_hugetlb_prefault_workers_for_test(bwamem_hugetlb_buf_t *buf, int nthreads,
                                           int *workers);

#ifdef __cplusplus
}
#endif

#endif /* BWA_HUGEPAGES_H */
