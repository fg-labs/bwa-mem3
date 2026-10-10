/* ext_reverse_copy.h -- the reversed left-extension windows of stage_seed_extension
 * (bwamem.cpp): dst[i] = src[n - 1 - i] for i in [0, n); n <= 0 copies nothing.
 *
 * n is a by-value parameter, so the loop bound is invariant and the compiler can vectorize the
 * loop (clang 19 on aarch64 reverses 32 bytes per iteration with rev64 + ext); __restrict spares
 * it a runtime overlap check. Read through the window struct instead, the bound must be reloaded
 * after every byte store, which may alias it, and the loop stays scalar. dst and src must not
 * overlap (the staging buffers vs the read or the fetched reference). */
#ifndef BWA_MEM3_EXT_REVERSE_COPY_H
#define BWA_MEM3_EXT_REVERSE_COPY_H

#include <stdint.h>

static inline void ext_reverse_copy(uint8_t *__restrict dst, const uint8_t *__restrict src, int64_t n)
{
    for (int64_t i = 0; i < n; ++i) dst[i] = src[n - 1 - i];
}

#endif
