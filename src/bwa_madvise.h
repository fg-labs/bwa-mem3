/*************************************************************************************
                           The MIT License

   BWA-MEM2  (Sequence alignment using Burrows-Wheeler Transform),
   Copyright (C) 2019  Intel Corporation, Heng Li.
*****************************************************************************************/
#ifndef BWA_MADVISE_H
#define BWA_MADVISE_H

#include <stddef.h>

#if defined(__linux__)
#include <sys/mman.h>
#endif

// Request transparent huge pages on a large allocation. The bwa-mem3 hg38
// index runtime working set is ~17 GB (cp_occ ~6 GB, sa_ms_byte ~2 GB,
// sa_ls_word ~8 GB) plus an ~800 MB pack table. With 4 KB pages the dTLB
// covers ~256 KB — a rounding error. Promoting to 2 MB THP pages (coverage
// ~128 MB per 64-entry TLB) dramatically cuts TLB pressure in the SMEM
// Occ-lookup, SA expansion, and reference-fetch loops. Best-effort: the
// hint is a no-op on non-Linux builds and on kernels without THP.
static inline void bwamem_madv_hugepage(void *p, size_t sz)
{
#if defined(__linux__) && defined(MADV_HUGEPAGE)
    if (p != NULL && sz > 0) (void)madvise(p, sz, MADV_HUGEPAGE);
#else
    (void)p; (void)sz;
#endif
}

#endif // BWA_MADVISE_H
