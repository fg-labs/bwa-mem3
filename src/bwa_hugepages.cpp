#include "bwa_hugepages.h"

#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>

#if defined(__linux__) && defined(USE_MIMALLOC)
/* Weak reference to mimalloc's reservation entry point, rather than including
 * <mimalloc.h>. The shipping executable whole-archives libmimalloc, so this
 * binds to the real symbol; a target that links libbwa.a WITHOUT mimalloc (the
 * doctest unit binary does) leaves it NULL, and we fall back at runtime instead
 * of failing to link. The weak attribute must be on our own declaration. */
extern "C" int mi_reserve_huge_os_pages_interleave(size_t pages, size_t numa_nodes,
                                                   size_t timeout_msecs) __attribute__((weak));
#endif

/* 1 GiB in bytes. */
static const unsigned long long BWAMEM_ONE_GB = 1024ULL * 1024ULL * 1024ULL;
/* Extra pages beyond the raw index footprint, covering 1 GB rounding of each
 * array plus mimalloc's own arena metadata. hg38 (~10.5 GiB of .bwt.2bit.64 +
 * .pac) sizes to 11 pages, so this +2 gives 13 to cover the ~10.76 GiB resident
 * index with headroom. The count is sized for the WHOLE index, although
 * mimalloc cannot place allocations above 2 GiB in its arena today (see
 * bwa_hugepages.h), so the pages for those arrays are reserved but unused. The
 * sizing is kept as is on purpose: the array sizes are not known from the file
 * sizes, and a follow-up change places the large arrays explicitly. */
static const size_t BWAMEM_HUGE_MARGIN_PAGES = 2;

size_t bwamem_huge_pages_needed(long long bwt_bytes, long long pac_bytes)
{
    if (bwt_bytes <= 0) return 0;
    unsigned long long total = (unsigned long long)bwt_bytes;
    if (pac_bytes > 0) total += (unsigned long long)pac_bytes;
    /* Ceil division as quotient + (remainder != 0), rather than
     * (total + ONE_GB - 1) / ONE_GB, so a near-SIZE_MAX total cannot overflow
     * the numerator. */
    size_t pages = (size_t)(total / BWAMEM_ONE_GB);
    if (total % BWAMEM_ONE_GB != 0) pages++;
    return pages + BWAMEM_HUGE_MARGIN_PAGES;
}

long bwamem_huge_pages_available(long free_pages, long resv_pages)
{
    if (free_pages < 0) return -1;
    /* Unreadable reserve count: fall back to the free count alone. A reserve
     * count of 0 is a valid reading and falls through to free - 0. */
    if (resv_pages < 0) return free_pages;
    return (resv_pages >= free_pages) ? 0 : free_pages - resv_pages;
}

bwamem_huge_outcome_t bwamem_classify_huge_reservation(
    int rc, long avail_before, long avail_after, size_t need, long *reserved_out)
{
    /* Actual pages reserved = the drop in the hugetlb available count. Only
     * trust it when both reads succeeded and the count did not rise (a rise
     * means another process perturbed the pool, so the delta is not ours to
     * interpret). A process that reserves or frees pages while mimalloc is
     * reserving can still skew the delta without a rise; counting exactly would
     * need mimalloc's per-call arena size, which this binding does not use. */
    long got = (avail_before >= 0 && avail_after >= 0 && avail_after <= avail_before)
                   ? (avail_before - avail_after)
                   : -1;
    if (reserved_out) *reserved_out = got;
    if (rc != 0) return BWAMEM_HUGE_FAILED;
    /* rc == 0 but the delta is unmeasurable: mimalloc only returns a status code
     * (>= 1 page per NUMA node), and the host-global available count did not give
     * us a trustworthy per-call delta, so the true count is unknown — do NOT claim a
     * full reservation. */
    if (got < 0) return BWAMEM_HUGE_UNVERIFIED;
    return (got >= (long)need) ? BWAMEM_HUGE_FULL : BWAMEM_HUGE_PARTIAL;
}

long bwamem_read_hugepage_counter(const char *dir, const char *name)
{
    char path[4096];
    snprintf(path, sizeof(path), "%s/%s", dir, name);
    FILE *f = fopen(path, "r");
    if (f == NULL) return -1;
    long n = -1;
    if (fscanf(f, "%ld", &n) != 1) n = -1;
    fclose(f);
    return n;
}

long bwamem_avail_hugepages_at(const char *dir)
{
    /* Separate statements, so the read order is defined (argument evaluation
     * order is not): resv first, then free (see the header). */
    long resv = bwamem_read_hugepage_counter(dir, "resv_hugepages");
    long freep = bwamem_read_hugepage_counter(dir, "free_hugepages");
    return bwamem_huge_pages_available(freep, resv);
}

#if defined(__linux__) && defined(USE_MIMALLOC)

static long long bwamem_file_size(const char *path)
{
    struct stat st;
    if (stat(path, &st) != 0) return -1;
    return (long long)st.st_size;
}

/* 1 GB pages still available to a new reservation (free minus reserved), or -1
 * if the pool does not exist (no 1 GB hugepage support / not reserved). */
static long bwamem_avail_1g_hugepages(void)
{
    return bwamem_avail_hugepages_at("/sys/kernel/mm/hugepages/hugepages-1048576kB");
}

void bwamem_reserve_huge_pages(const char *idxbase)
{
    char path[4096];
    snprintf(path, sizeof(path), "%s.bwt.2bit.64", idxbase);
    long long bwt = bwamem_file_size(path);
    snprintf(path, sizeof(path), "%s.pac", idxbase);
    long long pac = bwamem_file_size(path);

    /* MIMALLOC_RESERVE_HUGE_OS_PAGES makes mimalloc reserve at its own startup,
     * before this runs, which has already taken pages out of the available
     * count. Reserving again would double-count, and the shortfall check below
     * would misreport the pool as too small. */
    if (getenv("MIMALLOC_RESERVE_HUGE_OS_PAGES") != NULL) {
        fprintf(stderr, "[M::%s] --huge-pages: MIMALLOC_RESERVE_HUGE_OS_PAGES is "
                "set, so mimalloc already reserved its pages at startup; not "
                "reserving again\n", __func__);
        return;
    }

    size_t need = bwamem_huge_pages_needed(bwt, pac);
    if (need == 0) {
        fprintf(stderr, "[M::%s] --huge-pages: cannot size the index at %s "
                "(missing .bwt.2bit.64?); running on default pages\n",
                __func__, idxbase);
        return;
    }

    long avail = bwamem_avail_1g_hugepages();
    if (avail < 0) {
        fprintf(stderr, "[M::%s] --huge-pages: this host has no 1 GB hugepage "
                "pool; reserve some first (e.g. `echo %zu | sudo tee "
                "/sys/kernel/mm/hugepages/hugepages-1048576kB/nr_hugepages`). "
                "Running on default pages\n", __func__, need);
        return;
    }
    if ((size_t)avail < need) {
        fprintf(stderr, "[M::%s] --huge-pages: the index needs ~%zu 1 GB pages "
                "but only %ld are available (free minus already reserved); not "
                "reserving (raise nr_hugepages). Running on default pages\n",
                __func__, need, avail);
        return;
    }

    if (mi_reserve_huge_os_pages_interleave == NULL) {
        fprintf(stderr, "[M::%s] --huge-pages: this binary was built without "
                "mimalloc; running on default pages\n", __func__);
        return;
    }
    /* Interleave across all NUMA nodes (matches MIMALLOC_RESERVE_HUGE_OS_PAGES).
     * Bound the reservation with a finite timeout: a sufficient available-page
     * count was just checked, but reserving contiguous 1 GB blocks can still stall
     * under memory fragmentation, and a 0 (disabled) timeout would let that
     * block startup indefinitely. Budget 1 s per page (floor 5 s); pre-reserved
     * pages normally return at once, and on ETIMEDOUT the rc != 0 branch below
     * falls back to default pages. */
    size_t timeout_msecs = need * 1000;
    if (timeout_msecs < 5000) timeout_msecs = 5000;
    int rc = mi_reserve_huge_os_pages_interleave(need, 0, timeout_msecs);

    /* rc == 0 only means mimalloc reserved at least one page per NUMA node, which
     * can be fewer than `need` under fragmentation — a partial reservation shrinks
     * the arena. Confirm the pool actually gave up `need` pages via the
     * available-count delta rather than trusting rc alone. The counters are
     * host-global, so a concurrent reserver can still skew the delta (see the
     * classifier). */
    long avail_after = bwamem_avail_1g_hugepages();
    long got = -1;
    switch (bwamem_classify_huge_reservation(rc, avail, avail_after, need, &got)) {
    case BWAMEM_HUGE_FULL:
        fprintf(stderr, "[M::%s] --huge-pages: reserved %zu x 1 GB huge pages "
                "via mimalloc; allocations over 2 GiB (the largest index "
                "arrays) are not served from them and stay on default pages\n",
                __func__, need);
        break;
    case BWAMEM_HUGE_PARTIAL:
        fprintf(stderr, "[M::%s] --huge-pages: mimalloc reserved only %ld of %zu "
                "requested 1 GB pages; allocations that do not fit in them "
                "fall back to default pages\n", __func__, got, need);
        break;
    case BWAMEM_HUGE_FAILED:
        fprintf(stderr, "[M::%s] --huge-pages: mimalloc could not reserve %zu x "
                "1 GB pages (rc=%d); running on default pages\n",
                __func__, need, rc);
        break;
    case BWAMEM_HUGE_UNVERIFIED:
        fprintf(stderr, "[M::%s] --huge-pages: mimalloc accepted the 1 GB "
                "reservation but the host pool count could not be confirmed; "
                "the reserved count is unknown\n", __func__);
        break;
    }
}

#else /* not (Linux && mimalloc) */

void bwamem_reserve_huge_pages(const char *idxbase)
{
    (void)idxbase;
    fprintf(stderr, "[M::bwamem_reserve_huge_pages] --huge-pages is only "
            "supported on Linux with the mimalloc build; ignored\n");
}

#endif
