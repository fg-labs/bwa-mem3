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

static std::atomic<int> g_hugetlb_index_mode{BWAMEM_HUGETLB_OFF};

void bwamem_hugetlb_set_index_mode(bwamem_hugetlb_mode_t mode)
{
    g_hugetlb_index_mode.store((int)mode, std::memory_order_relaxed);
}

bwamem_hugetlb_mode_t bwamem_hugetlb_index_mode(void)
{
    return (bwamem_hugetlb_mode_t)g_hugetlb_index_mode.load(std::memory_order_relaxed);
}

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

int bwamem_hugetlb_map(size_t bytes, bwamem_hugetlb_buf_t *out)
{
    out->base = NULL;
    out->map_bytes = 0;
    size_t len = bwamem_hugetlb_round_up(bytes);
    if (len == 0) return EINVAL;
    if (bwamem_hugetlb_index_mode() == BWAMEM_HUGETLB_FORCE_FAIL) return ECANCELED;
    void *p = mmap(NULL, len, PROT_READ | PROT_WRITE,
                   MAP_PRIVATE | MAP_ANONYMOUS | MAP_HUGETLB | MAP_HUGE_1GB, -1, 0);
    if (p == MAP_FAILED) return errno != 0 ? errno : ENOMEM;
    out->base = p;
    out->map_bytes = len;
    return 0;
}

namespace {

struct PrefaultJob {
    std::vector<char *> pages;      /* start of every 1 GiB page to fault in */
    std::vector<int>    owner;      /* index into bufs of each page's mapping */
    std::vector<std::atomic<int>> *failed;  /* per mapping: a fault failed */
    std::atomic<size_t> next{0};
    std::atomic<bool>   go{false};          /* set once every worker exists */
};

/* Fault in one 1 GiB page for writing. MADV_POPULATE_WRITE (Linux 5.14) reports
 * a failure as an error instead of a SIGBUS. Any failure, including EINVAL from
 * a kernel without it, is returned rather than retried by touching the page:
 * the mapping reserved its pages, but a limit charged at fault time (a hugetlb
 * cgroup's 1 GB max) can still fail the first touch, and that is a SIGBUS that
 * kills the run instead of a fallback to default pages. */
bool prefault_page(char *page)
{
#ifndef MADV_POPULATE_WRITE
#define MADV_POPULATE_WRITE 23
#endif
    return madvise(page, (size_t)BWAMEM_ONE_GB, MADV_POPULATE_WRITE) == 0;
}

void *prefault_worker(void *arg)
{
    PrefaultJob *job = (PrefaultJob *)arg;
    while (!job->go.load(std::memory_order_acquire)) sched_yield();
    for (;;) {
        size_t i = job->next.fetch_add(1, std::memory_order_relaxed);
        if (i >= job->pages.size()) break;
        if (!prefault_page(job->pages[i]))
            (*job->failed)[(size_t)job->owner[i]].store(1, std::memory_order_relaxed);
    }
    return NULL;
}

}  // namespace

int bwamem_hugetlb_prefault(bwamem_hugetlb_buf_t *bufs, int n, int nthreads)
{
    if (bufs == NULL || n <= 0) return 0;
    std::vector<std::atomic<int>> failed((size_t)n);
    PrefaultJob job;
    job.failed = &failed;
    for (int b = 0; b < n; b++) {
        failed[(size_t)b].store(0, std::memory_order_relaxed);
        for (size_t off = 0; bufs[b].base != NULL && off < bufs[b].map_bytes;
             off += (size_t)BWAMEM_ONE_GB) {
            job.pages.push_back((char *)bufs[b].base + off);
            job.owner.push_back(b);
        }
    }
    size_t nt = nthreads > 0 ? (size_t)nthreads : 1;
    if (nt > job.pages.size()) nt = job.pages.size();
    /* The calling thread is one of the workers; a thread that cannot be
     * created just leaves its share to the others. No worker starts faulting
     * until every thread exists: creating a thread maps its stack, which takes
     * the address-space lock for writing, and that waits for every fault in
     * flight (each holds the lock for reading while it zeroes a whole 1 GB
     * page) and blocks new ones behind it. Measured on a 16-thread load of
     * hg38, starting each worker as it was created left only two or three
     * faulting at a time (0.59 s for 10 pages instead of about 0.15 s). */
    std::vector<pthread_t> tids;
    for (size_t t = 1; t < nt; t++) {
        pthread_t tid;
        if (pthread_create(&tid, NULL, prefault_worker, &job) == 0) tids.push_back(tid);
    }
    job.go.store(true, std::memory_order_release);
    prefault_worker(&job);
    for (pthread_t tid : tids) pthread_join(tid, NULL);

    int n_failed = 0;
    for (int b = 0; b < n; b++) {
        if (failed[(size_t)b].load(std::memory_order_relaxed)) {
            bwamem_hugetlb_unmap(&bufs[b]);
            n_failed++;
        }
    }
    return n_failed;
}

#else /* not Linux */

int bwamem_hugetlb_map(size_t bytes, bwamem_hugetlb_buf_t *out)
{
    (void)bytes;
    out->base = NULL;
    out->map_bytes = 0;
    return ENOSYS;
}

int bwamem_hugetlb_prefault(bwamem_hugetlb_buf_t *bufs, int n, int nthreads)
{
    (void)bufs; (void)n; (void)nthreads;
    return 0;
}

#endif
