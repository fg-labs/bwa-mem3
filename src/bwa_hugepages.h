#ifndef BWA_HUGEPAGES_H
#define BWA_HUGEPAGES_H

#include <stddef.h>

/* Opt-in 1 GB huge pages for the FM-index / suffix-array structures (issue
 * #402, measured in #377). Seeding is memory-latency-bound: the SA/BWT lookups
 * are random accesses over multi-gigabyte arrays, and even on transparent 2 MB
 * pages the page walks for them take a measurable share of the cycles.
 *
 * With `--huge-pages`, FMI_search::load_index maps the occurrence table and the
 * suffix-array samples directly with mmap(MAP_HUGETLB | MAP_HUGE_1GB) from the
 * host's 1 GB hugetlb pool, and faults the pages in over several threads before
 * filling them. Every step is best-effort: when the pool is missing or too
 * small, or a mapping or its pre-fault fails (as it does on a kernel before
 * 5.14, which lacks MADV_POPULATE_WRITE), that array falls back to the default
 * allocation. The arrays' contents, and so the output, are the same
 * either way. Linux only; elsewhere every mapping fails with ENOSYS and the
 * default allocation is used. */

#ifdef __cplusplus
extern "C" {
#endif

/* Whether FMI_search::load_index maps its large arrays on 1 GB pages. */
typedef enum {
    BWAMEM_HUGETLB_OFF = 0,        /* default allocation (no --huge-pages) */
    BWAMEM_HUGETLB_ON = 1,         /* try 1 GB mappings, fall back on failure */
    BWAMEM_HUGETLB_FORCE_FAIL = 2  /* tests only: every mapping fails with ECANCELED */
} bwamem_hugetlb_mode_t;

void bwamem_hugetlb_set_index_mode(bwamem_hugetlb_mode_t mode);
bwamem_hugetlb_mode_t bwamem_hugetlb_index_mode(void);

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
 * non-Linux builds, ECANCELED in BWAMEM_HUGETLB_FORCE_FAIL mode, and EINVAL when
 * `bytes` is 0 or too large. */
int bwamem_hugetlb_map(size_t bytes, bwamem_hugetlb_buf_t *out);

/* Unmap `*b` (if mapped) and reset it to unmapped. NULL is a no-op. */
void bwamem_hugetlb_unmap(bwamem_hugetlb_buf_t *b);

/* Fault in every 1 GiB page of the `n` mappings in `bufs` before they are
 * filled, using up to `nthreads` threads, one page per task. The kernel zeroes
 * a 1 GB page in the fault that first touches it, on the faulting thread; left
 * to the index readers, several readers share one page and wait on that one
 * zeroing, which made start-up slower than on default pages. Spreading the
 * faults over threads removes that serialization, and the readers then only
 * copy. Faulting uses MADV_POPULATE_WRITE, so it needs Linux 5.14 or later; on
 * an older kernel every pre-fault fails. A mapping whose pre-fault fails is
 * unmapped (and reset), so the caller falls back for it. Returns the number of
 * mappings unmapped this way. */
int bwamem_hugetlb_prefault(bwamem_hugetlb_buf_t *bufs, int n, int nthreads);

#ifdef __cplusplus
}

/* Sets the index huge-page mode for its lifetime and restores the mode it found
 * when it goes out of scope, so a caller that turns --huge-pages on (main_mem)
 * leaves later index loads in the same process as they were. */
class bwamem_hugetlb_mode_scope {
public:
    explicit bwamem_hugetlb_mode_scope(bwamem_hugetlb_mode_t mode)
        : prev_(bwamem_hugetlb_index_mode()) { bwamem_hugetlb_set_index_mode(mode); }
    ~bwamem_hugetlb_mode_scope() { bwamem_hugetlb_set_index_mode(prev_); }
    bwamem_hugetlb_mode_scope(const bwamem_hugetlb_mode_scope &) = delete;
    bwamem_hugetlb_mode_scope &operator=(const bwamem_hugetlb_mode_scope &) = delete;
private:
    bwamem_hugetlb_mode_t prev_;
};
#endif

#endif /* BWA_HUGEPAGES_H */
