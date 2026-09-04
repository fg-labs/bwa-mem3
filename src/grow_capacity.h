#ifndef GROW_CAPACITY_H
#define GROW_CAPACITY_H

#include <stddef.h>
#include <stdint.h>

/* Geometric growth policy shared by the scratch-buffer growers (BamRecScratch's
 * grow<T> and bwamem's scratch_grow): double the current capacity, but never
 * below `need`, and -- when the capacity is counted in elements of `elem_size`
 * bytes (> 1) -- clamp back to `need` if the doubled value would overflow once
 * scaled by elem_size. Pass elem_size == 1 for a byte-counted buffer (the
 * returned value is itself the byte count, so no scale overflow applies).
 *
 * Callers that allocate `result * elem_size` must still guard `need` itself
 * against `need > SIZE_MAX / elem_size` before calling (this only bounds the
 * doubled value, not a caller-supplied `need`). */
static inline size_t grow_capacity(size_t cur, size_t need, size_t elem_size) {
    size_t want = cur * 2;
    if (want < need) want = need;
    if (elem_size > 1 && want > SIZE_MAX / elem_size) want = need;
    return want;
}

#endif /* GROW_CAPACITY_H */
