/* SPDX-License-Identifier: MIT */
/* New code for bwa-mem3 (not from upstream bwa): allocation that stops with a message on
 * failure. The upstream code dereferences the result of every malloc/calloc/realloc/strdup
 * unchecked; these wrappers keep its call sites and turn a failed allocation into a clean
 * exit instead of a NULL write (and, for realloc, a lost original block). They never change
 * what a successful allocation returns. */
#ifndef BWA_MEM3_ALN_ALLOC_H
#define BWA_MEM3_ALN_ALLOC_H

#include <stdlib.h>
#include <string.h>

#include "../utils.h"

namespace bwa_aln {

static inline void *aln_nonnull(void *p, size_t bytes)
{
	if (!p && bytes) _err_fatal_simple_core("aln", "out of memory");
	return p;
}
static inline void *aln_malloc(size_t n) { return aln_nonnull(malloc(n), n); }
static inline void *aln_calloc(size_t n, size_t sz) { return aln_nonnull(calloc(n, sz), n * sz); }
static inline void *aln_realloc(void *p, size_t n) { return aln_nonnull(realloc(p, n), n); }
static inline char *aln_strdup(const char *s) { return (char*)aln_nonnull(strdup(s), 1); }

} // namespace bwa_aln

#endif
