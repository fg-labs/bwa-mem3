#ifndef SMEM_DEDUP_H
#define SMEM_DEDUP_H

#include <stdint.h>
#include "FMI_search.h"   // SMEM

#ifdef __cplusplus
extern "C" {
#endif

/* Compact adjacent fully-identical SMEMs in a sorted array; return the new length.
 * An entry is kept only if it differs from the previous kept entry on (rid,m,n,k,l,s).
 * Precondition: identical SMEMs are adjacent (true after sortSMEMs + per-read
 * intv_lt1 sort in mem_collect_smem).
 * O(n), no allocation. n <= 1 is handled safely (returns n). */
int64_t smem_dedup_inplace(SMEM *a, int64_t n);

#ifdef __cplusplus
}
#endif

#endif /* SMEM_DEDUP_H */
