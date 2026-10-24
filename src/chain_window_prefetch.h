/* chain_window_prefetch.h -- the packed-reference cache lines the extension driver hints before
 * it fetches a chain's reference window (prefetch_chain_window in mem_chain2aln_across_reads_V2,
 * bwamem.cpp).
 *
 * The window a chain's first seed implies is [rbeg - before, rbeg + after) in doubled (both
 * strand) coordinates: the seed's own term of the chain's window bounds. It is clamped to the
 * seed's strand, mapped to the forward-strand bases the packed reference stores (a reverse-strand
 * base x is forward base 2 * l_pac - 1 - x), and then to the 64-byte lines of pac (four bases per
 * byte) it spans. A long window hints its first max_lines lines and its last line, since the
 * fetch's first touches are the window's two ends. Pure hints: nothing reads the result except
 * the prefetch instructions. */
#ifndef BWA_MEM3_CHAIN_WINDOW_PREFETCH_H
#define BWA_MEM3_CHAIN_WINDOW_PREFETCH_H

#include <stdint.h>

/* Write the line-aligned addresses to hint for the window into out[] (room for max_lines + 1)
 * and return how many there are: 0 when the clamped window is empty, else every line the window
 * spans when that is at most max_lines + 1, else the first max_lines lines and the last. A
 * max_lines below 1 hints only the last line (out[] then needs room for 1). */
static inline int chain_window_pac_lines(const uint8_t *pac, int64_t l_pac, int64_t rbeg,
                                         int64_t before, int64_t after, int max_lines,
                                         uintptr_t *out)
{
    const int64_t b = rbeg - before, e = rbeg + after;
    int64_t lo, hi;  /* forward-strand base range [lo, hi] */
    if (rbeg < l_pac) {
        lo = b > 0 ? b : 0;
        hi = (e < l_pac ? e : l_pac) - 1;
    } else {
        const int64_t rb = b > l_pac ? b : l_pac;
        const int64_t re = e < (l_pac << 1) ? e : (l_pac << 1);
        lo = (l_pac << 1) - re;
        hi = (l_pac << 1) - 1 - rb;
    }
    if (lo > hi) return 0;
    const uintptr_t first = ((uintptr_t)pac + (uintptr_t)(lo >> 2)) & ~(uintptr_t)63;
    const uintptr_t last  = ((uintptr_t)pac + (uintptr_t)(hi >> 2)) & ~(uintptr_t)63;
    int n = 0;
    for (uintptr_t at = first; at <= last && n < max_lines; at += 64) out[n++] = at;
    if (n == 0 || out[n - 1] != last) out[n++] = last;
    return n;
}

#endif
