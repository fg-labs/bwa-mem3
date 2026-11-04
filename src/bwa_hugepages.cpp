#include "bwa_hugepages.h"

#include <errno.h>
#include <pthread.h>
#include <sched.h>
#include <stdio.h>
#include <stdlib.h>

#include <atomic>
#include <vector>

#if defined(__linux__)
#include <sys/mman.h>
#endif

/* 1 GiB in bytes. */
static const unsigned long long BWAMEM_ONE_GB = 1024ULL * 1024ULL * 1024ULL;

size_t bwamem_hugetlb_round_up(size_t bytes)
{
    if (bytes == 0) return 0;
    const size_t page = (size_t)BWAMEM_ONE_GB;
    size_t pages = bytes / page + (bytes % page != 0);
    if (pages > (size_t)-1 / page) return 0;
    return pages * page;
}

void bwamem_hugetlb_unmap(bwamem_hugetlb_buf_t *b)
{
    if (b == NULL) return;
#if defined(__linux__)
    if (b->base != NULL) munmap(b->base, b->map_bytes);
#endif
    b->base = NULL;
    b->map_bytes = 0;
}

#if defined(__linux__)

#ifndef MAP_HUGE_SHIFT
#define MAP_HUGE_SHIFT 26
#endif
#ifndef MAP_HUGE_1GB
#define MAP_HUGE_1GB (30 << MAP_HUGE_SHIFT)
#endif
#ifndef MADV_POPULATE_WRITE
#define MADV_POPULATE_WRITE 23
#endif

int bwamem_hugetlb_map(size_t bytes, bwamem_hugetlb_buf_t *out)
{
    out->base = NULL;
    out->map_bytes = 0;
    size_t len = bwamem_hugetlb_round_up(bytes);
    if (len == 0) return EINVAL;
    void *p = mmap(NULL, len, PROT_READ | PROT_WRITE,
                   MAP_PRIVATE | MAP_ANONYMOUS | MAP_HUGETLB | MAP_HUGE_1GB, -1, 0);
    if (p == MAP_FAILED) return errno != 0 ? errno : ENOMEM;
    out->base = p;
    out->map_bytes = len;
    return 0;
}

int bwamem_hugetlb_map_plain_for_test(size_t bytes, bwamem_hugetlb_buf_t *out)
{
    out->base = NULL;
    out->map_bytes = 0;
    size_t len = bwamem_hugetlb_round_up(bytes);
    if (len == 0) return EINVAL;
    void *p = mmap(NULL, len, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
    if (p == MAP_FAILED) return errno != 0 ? errno : ENOMEM;
    out->base = p;
    out->map_bytes = len;
    return 0;
}

int bwamem_hugetlb_map_protected_for_test(size_t bytes, bwamem_hugetlb_buf_t *out)
{
    const int err = bwamem_hugetlb_map_plain_for_test(bytes, out);
    if (err != 0) return err;
    if (mprotect(out->base, out->map_bytes, PROT_NONE) == 0) return 0;
    const int protect_err = errno;
    bwamem_hugetlb_unmap(out);
    return protect_err;
}

namespace {

struct PrefaultJob {
    char               *base = nullptr;  /* the mapping being faulted in */
    size_t              n_pages = 0;     /* its 1 GiB pages */
    std::atomic<size_t> next{0};
    std::atomic<bool>   failed{false};   /* a page failed: stop, unmap after */
    std::atomic<bool>   go{false};       /* set once every worker exists */
};

/* Fault in one 1 GiB page for writing. MADV_POPULATE_WRITE (Linux 5.14) reports
 * a failure as an error instead of a SIGBUS. Any failure, including EINVAL from
 * a kernel without it, is returned rather than retried by touching the page:
 * the mapping reserved its pages, but a hugetlb cgroup's usage limit
 * (hugetlb.1GB.max) is charged at fault time and can still fail the first
 * touch, and that would be a SIGBUS that kills the run instead of a fallback to
 * default pages. */
bool prefault_page(char *page)
{
    return madvise(page, (size_t)BWAMEM_ONE_GB, MADV_POPULATE_WRITE) == 0;
}

void *prefault_worker(void *arg)
{
    PrefaultJob *job = (PrefaultJob *)arg;
    while (!job->go.load(std::memory_order_acquire)) sched_yield();
    for (;;) {
        size_t i = job->next.fetch_add(1, std::memory_order_relaxed);
        if (i >= job->n_pages) break;
        /* Once one page has failed, the whole mapping is unmapped after the
         * join, so zeroing the remaining pages would be wasted. */
        if (job->failed.load(std::memory_order_relaxed)) break;
        if (!prefault_page(job->base + i * (size_t)BWAMEM_ONE_GB))
            job->failed.store(true, std::memory_order_relaxed);
    }
    return NULL;
}

}  // namespace

/* A non-null workers output is the test-only no-caller mode. Production
 * keeps the caller as a worker, including when every pthread_create fails. */
static int hugetlb_prefault(bwamem_hugetlb_buf_t *buf, int nthreads, int *workers)
{
    if (workers != NULL) *workers = 0;
    if (buf == NULL || buf->base == NULL || buf->map_bytes == 0) return 0;
    PrefaultJob job;
    job.base = (char *)buf->base;
    job.n_pages = buf->map_bytes / (size_t)BWAMEM_ONE_GB;
    size_t nt = nthreads > 0 ? (size_t)nthreads : 1;
    if (nt > job.n_pages) nt = job.n_pages;
    /* The calling thread is one of the workers; a thread that cannot be
     * created just leaves its share to the others. No worker starts faulting
     * until every thread exists: creating a thread maps its stack, which takes
     * the address-space lock for writing, and that waits for every fault in
     * flight (each holds the lock for reading while it zeroes a whole 1 GB
     * page) and blocks new ones behind it. Measured on a 16-thread load of
     * hg38 with ten 1 GB pages mapped, starting each worker as it was created
     * left only two or three faulting at a time (0.59 s instead of about
     * 0.15 s). */
    std::vector<pthread_t> tids;
    tids.reserve(nt > 0 ? nt - 1 : 0);
    for (size_t t = 1; t < nt; t++) {
        pthread_t tid;
        if (pthread_create(&tid, NULL, prefault_worker, &job) == 0) tids.push_back(tid);
    }
    if (workers != NULL) *workers = (int)tids.size();
    job.go.store(true, std::memory_order_release);
    if (workers == NULL) prefault_worker(&job);
    for (pthread_t tid : tids) pthread_join(tid, NULL);

    if (!job.failed.load(std::memory_order_relaxed)) return 0;
    bwamem_hugetlb_unmap(buf);
    return 1;
}

int bwamem_hugetlb_prefault(bwamem_hugetlb_buf_t *buf, int nthreads)
{
    return hugetlb_prefault(buf, nthreads, NULL);
}

int bwamem_hugetlb_prefault_workers_for_test(bwamem_hugetlb_buf_t *buf, int nthreads,
                                           int *workers)
{
    return hugetlb_prefault(buf, nthreads, workers);
}

#else /* not Linux */

int bwamem_hugetlb_map(size_t bytes, bwamem_hugetlb_buf_t *out)
{
    out->base = NULL;
    out->map_bytes = 0;
    if (bwamem_hugetlb_round_up(bytes) == 0) return EINVAL;
    return ENOSYS;
}

int bwamem_hugetlb_map_plain_for_test(size_t bytes, bwamem_hugetlb_buf_t *out)
{
    out->base = NULL;
    out->map_bytes = 0;
    if (bwamem_hugetlb_round_up(bytes) == 0) return EINVAL;
    return ENOSYS;
}

int bwamem_hugetlb_map_protected_for_test(size_t bytes, bwamem_hugetlb_buf_t *out)
{
    return bwamem_hugetlb_map_plain_for_test(bytes, out);
}

int bwamem_hugetlb_prefault(bwamem_hugetlb_buf_t *buf, int nthreads)
{
    (void)buf; (void)nthreads;
    return 0;
}

int bwamem_hugetlb_prefault_workers_for_test(bwamem_hugetlb_buf_t *buf, int nthreads,
                                           int *workers)
{
    *workers = 0;
    return bwamem_hugetlb_prefault(buf, nthreads);
}

#endif
