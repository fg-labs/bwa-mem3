/* SPDX-License-Identifier: MIT */
#ifndef BWAMEM3_U8VEC_SCRATCH_H
#define BWAMEM3_U8VEC_SCRATCH_H

#include <cstdint>
#include <cstdlib>

#include "kvec.h"

/* Named typedef for a kvec of uint8_t. kvec_t(uint8_t) on its own expands
 * to a fresh anonymous struct at each macro call, so it can't be used
 * across function boundaries. Naming it once here makes it a regular
 * struct type that callers can declare thread_local, pass by pointer,
 * etc. The kv_* macros operate on the public (n, m, a) fields and work
 * with this typedef unchanged. */
typedef struct { size_t n, m; uint8_t *a; } u8vec_t;

/* Thread-local scratch holder. Initializes the kvec on construction
 * (n = m = 0, a = NULL) and frees the underlying buffer on thread exit
 * (or process exit for threads that don't terminate first). Use as:
 *
 *     static thread_local u8vec_scratch_t t_buf;
 *     if (t_buf.v.m < want) kv_resize(uint8_t, t_buf.v, want);
 *     ... write into t_buf.v.a ... no free needed ...
 *
 * One scratch instance per logical buffer per thread; the per-call cost
 * is a single capacity check + amortized realloc-on-grow via kv_resize.
 */
struct u8vec_scratch_t {
    u8vec_t v;
    u8vec_scratch_t() { kv_init(v); }
    ~u8vec_scratch_t() { kv_destroy(v); }
};

#endif
