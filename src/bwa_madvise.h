/*************************************************************************************
                           The MIT License

   BWA-MEM2  (Sequence alignment using Burrows-Wheeler Transform),
   Copyright (C) 2019  Intel Corporation, Heng Li.
*****************************************************************************************/
#ifndef BWA_MADVISE_H
#define BWA_MADVISE_H

#include <stddef.h>
#include <stdint.h>

#if defined(__linux__)
#include <errno.h>
#include <stdio.h>
#include <string.h>
#include <sys/mman.h>
#include <unistd.h>
#endif

extern int bwa_verbose;

// Shrink [p, p + sz) to the largest sub-range whose start and end are both
// multiples of `page` (a power of two): start = round_up(p, page),
// end = round_down(p + sz, page). Stores the start in *start_out and returns the
// length, or returns 0 (leaving *start_out untouched) when no whole page fits,
// `page` is not a power of two, or p + sz overflows.
//
// madvise() rejects a start that is not page-aligned with EINVAL. The index
// arrays come from _mm_malloc with 64-byte alignment, and the pack table from
// calloc; glibc serves such large requests from mmap and returns a pointer just
// past its chunk header, so the start is never page-aligned and the whole hint
// would be dropped. Advising only the aligned interior gives up less than one
// page at each end.
static inline size_t bwamem_page_align_inward(uintptr_t p, size_t sz, size_t page,
                                              uintptr_t *start_out)
{
    if (page == 0 || (page & (page - 1)) != 0) return 0;
    const uintptr_t mask = (uintptr_t)page - 1;
    if (p > UINTPTR_MAX - mask) return 0;   // rounding p up would wrap
    const uintptr_t start = (p + mask) & ~mask;
    // If p + sz wraps, the wrapped end is below p and so below start: empty.
    const uintptr_t end = (p + sz) & ~mask;
    if (end <= start) return 0;
    *start_out = start;
    return (size_t)(end - start);
}

// Request transparent huge pages on a large allocation. The bwa-mem3 hg38
// index runtime working set is ~17 GB (cp_occ ~6 GB, sa_ms_byte ~2 GB,
// sa_ls_word ~8 GB) plus an ~800 MB pack table. With 4 KB pages the dTLB
// covers ~256 KB — a rounding error. Promoting to 2 MB THP pages (coverage
// ~128 MB per 64-entry TLB) dramatically cuts TLB pressure in the SMEM
// Occ-lookup, SA expansion, and reference-fetch loops. Best-effort: the
// hint is a no-op on non-Linux builds and on kernels without THP; a failure is
// reported only at -v 4 and never changes output.
static inline void bwamem_madv_hugepage(void *p, size_t sz)
{
#if defined(__linux__) && defined(MADV_HUGEPAGE)
    if (p == NULL || sz == 0) return;
    const long page = sysconf(_SC_PAGESIZE);
    if (page <= 0) return;
    uintptr_t start = 0;
    const size_t len = bwamem_page_align_inward((uintptr_t)p, sz, (size_t)page, &start);
    if (len == 0) return;
    if (madvise((void *)start, len, MADV_HUGEPAGE) != 0 && bwa_verbose >= 4) {
        const int e = errno;
        fprintf(stderr, "[M::%s] madvise(MADV_HUGEPAGE, %zu bytes) failed: %s\n",
                __func__, len, strerror(e));
    }
#else
    (void)p; (void)sz;
#endif
}

#endif // BWA_MADVISE_H
