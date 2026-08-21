/*************************************************************************************
                          The MIT License

   BWA-MEM2  (Sequence alignment using Burrows-Wheeler Transform),
   Copyright (C) 2019  Intel Corporation, Heng Li.

   Permission is hereby granted, free of charge, to any person obtaining
   a copy of this software and associated documentation files (the
   "Software"), to deal in the Software without restriction, including
   without limitation the rights to use, copy, modify, merge, publish,
   distribute, sublicense, and/or sell copies of the Software, and to
   permit persons to whom the Software is furnished to do so, subject to
   the following conditions:

   The above copyright notice and this permission notice shall be
   included in all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
   BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.

   NEON Utility Macros and Functions for Apple Silicon Port
   Native NEON implementations for performance-critical paths
*****************************************************************************************/

#ifndef NEON_UTILS_H
#define NEON_UTILS_H

#if defined(__ARM_NEON) || defined(__aarch64__) || defined(APPLE_SILICON)

#include <arm_neon.h>
#include <stdint.h>

/*
 * NEON vector width constants
 * NEON provides 128-bit vectors (16 bytes, 8 halfwords, 4 words, 2 doublewords)
 */
#define NEON_WIDTH_BYTES 16
#define NEON_WIDTH8  16   /* 16 x 8-bit elements */
#define NEON_WIDTH16 8    /* 8 x 16-bit elements */
#define NEON_WIDTH32 4    /* 4 x 32-bit elements */
#define NEON_WIDTH64 2    /* 2 x 64-bit elements */

/*
 * Helper macros for common NEON operations
 */

/* Zero vector initialization */
#define NEON_ZERO_U8()   vdupq_n_u8(0)
#define NEON_ZERO_U16()  vdupq_n_u16(0)
#define NEON_ZERO_S8()   vdupq_n_s8(0)
#define NEON_ZERO_S16()  vdupq_n_s16(0)

/* Broadcast scalar to all lanes */
#define NEON_SET1_U8(x)  vdupq_n_u8(x)
#define NEON_SET1_U16(x) vdupq_n_u16(x)
#define NEON_SET1_S8(x)  vdupq_n_s8(x)
#define NEON_SET1_S16(x) vdupq_n_s16(x)

/* Load/Store operations */
#define NEON_LOAD_U8(ptr)   vld1q_u8((const uint8_t*)(ptr))
#define NEON_LOAD_U16(ptr)  vld1q_u16((const uint16_t*)(ptr))
#define NEON_LOAD_S8(ptr)   vld1q_s8((const int8_t*)(ptr))
#define NEON_LOAD_S16(ptr)  vld1q_s16((const int16_t*)(ptr))

#define NEON_STORE_U8(ptr, v)   vst1q_u8((uint8_t*)(ptr), v)
#define NEON_STORE_U16(ptr, v)  vst1q_u16((uint16_t*)(ptr), v)
#define NEON_STORE_S8(ptr, v)   vst1q_s8((int8_t*)(ptr), v)
#define NEON_STORE_S16(ptr, v)  vst1q_s16((int16_t*)(ptr), v)

/*
 * Saturating arithmetic (critical for Smith-Waterman)
 * These are direct equivalents to SSE/AVX saturating operations
 */
/* Unsigned saturating add */
#define NEON_ADDS_U8(a, b)   vqaddq_u8(a, b)
#define NEON_ADDS_U16(a, b)  vqaddq_u16(a, b)

/* Unsigned-saturating add of a SIGNED addend into an unsigned accumulator:
 * clamp(a + (int8)b, 0, 255). NEON's USQADD; `b` is passed as uint8x16_t and
 * reinterpreted, matching how the scoring table is laid out. The u8 mate-rescue
 * diagonal uses it to fold the de-bias vqsubq_u8 out of the h00->m11 chain. */
#define NEON_SQADD_U8(a, b)  vsqaddq_u8((a), vreinterpretq_s8_u8(b))

/* Unsigned saturating subtract */
#define NEON_SUBS_U8(a, b)   vqsubq_u8(a, b)
#define NEON_SUBS_U16(a, b)  vqsubq_u16(a, b)

/* Signed saturating add/sub */
#define NEON_ADDS_S8(a, b)   vqaddq_s8(a, b)
#define NEON_ADDS_S16(a, b)  vqaddq_s16(a, b)
#define NEON_SUBS_S8(a, b)   vqsubq_s8(a, b)
#define NEON_SUBS_S16(a, b)  vqsubq_s16(a, b)

/*
 * Max/Min operations
 */
#define NEON_MAX_U8(a, b)   vmaxq_u8(a, b)
#define NEON_MAX_U16(a, b)  vmaxq_u16(a, b)
#define NEON_MAX_S8(a, b)   vmaxq_s8(a, b)
#define NEON_MAX_S16(a, b)  vmaxq_s16(a, b)

#define NEON_MIN_U8(a, b)   vminq_u8(a, b)
#define NEON_MIN_U16(a, b)  vminq_u16(a, b)
#define NEON_MIN_S8(a, b)   vminq_s8(a, b)
#define NEON_MIN_S16(a, b)  vminq_s16(a, b)

/*
 * Comparison operations
 */
#define NEON_CMPEQ_U8(a, b)  vceqq_u8(a, b)
#define NEON_CMPEQ_U16(a, b) vceqq_u16(a, b)
#define NEON_CMPGT_U8(a, b)  vcgtq_u8(a, b)
#define NEON_CMPGT_U16(a, b) vcgtq_u16(a, b)
#define NEON_CMPGT_S8(a, b)  vcgtq_s8(a, b)
#define NEON_CMPGT_S16(a, b) vcgtq_s16(a, b)

/*
 * Bitwise operations
 */
#define NEON_AND(a, b)    vandq_u8(a, b)
#define NEON_OR(a, b)     vorrq_u8(a, b)
#define NEON_XOR(a, b)    veorq_u8(a, b)
#define NEON_NOT(a)       vmvnq_u8(a)

/* Bitwise select: (mask & a) | (~mask & b) - equivalent to _mm_blendv_epi8 */
#define NEON_BLENDV_U8(a, b, mask) vbslq_u8(mask, a, b)
#define NEON_BLENDV_U16(a, b, mask) vbslq_u16(mask, a, b)

/*
 * Shift operations
 */
#define NEON_SHL_U8(a, n)   vshlq_n_u8(a, n)
#define NEON_SHR_U8(a, n)   vshrq_n_u8(a, n)
#define NEON_SHL_U16(a, n)  vshlq_n_u16(a, n)
#define NEON_SHR_U16(a, n)  vshrq_n_u16(a, n)

/*
 * Extract high bit from each byte (movemask equivalent)
 * NEON doesn't have a direct movemask, so we need to implement it
 */
static inline uint16_t neon_movemask_u8(uint8x16_t v) {
    /* Per-byte positional weights (1<<lane), repeated across both 8-byte
     * halves. A compile-time vector literal, so it materializes via MOVI
     * immediates rather than a data-section load. */
    const uint8x16_t bit_mask = {1, 2, 4, 8, 16, 32, 64, 128,
                                 1, 2, 4, 8, 16, 32, 64, 128};

    /* Broadcast each byte's high bit across the whole lane (0xFF/0x00) via an
     * arithmetic right shift, then keep that lane's positional weight with vand.
     * NOTE: a *logical* shift (vshrq_n_u8) yields 0x01, and 0x01 & (1<<lane) is
     * 0 for every lane except 0 and 8 -- that silently collapses the mask. The
     * signed shift is required for vand to be equivalent to the vshl packing. */
    uint8x16_t hi = vreinterpretq_u8_s8(vshrq_n_s8(vreinterpretq_s8_u8(v), 7));
    uint8x16_t bits = vandq_u8(hi, bit_mask);

    /* Sum each 8-byte half: low half -> low byte of result, high -> high. */
    return (uint16_t)vaddv_u8(vget_low_u8(bits)) |
           ((uint16_t)vaddv_u8(vget_high_u8(bits)) << 8);
}

/*
 * Horizontal max within vector
 * Returns the maximum value across all lanes
 */
static inline uint8_t neon_hmax_u8(uint8x16_t v) {
    return vmaxvq_u8(v);
}

static inline uint16_t neon_hmax_u16(uint16x8_t v) {
    return vmaxvq_u16(v);
}

static inline int8_t neon_hmax_s8(int8x16_t v) {
    return vmaxvq_s8(v);
}

static inline int16_t neon_hmax_s16(int16x8_t v) {
    return vmaxvq_s16(v);
}

/*
 * Table lookup (shuffle equivalent)
 * NEON's vtbl/vqtbl can implement byte shuffles like _mm_shuffle_epi8
 */
static inline uint8x16_t neon_shuffle_u8(uint8x16_t tbl, uint8x16_t idx) {
    return vqtbl1q_u8(tbl, idx);
}

/*
 * Prefetch hint (portable)
 */
#define NEON_PREFETCH(addr) __builtin_prefetch((const void*)(addr), 0, 3)
#define NEON_PREFETCH_WRITE(addr) __builtin_prefetch((const void*)(addr), 1, 3)

/*
 * Memory alignment helpers for Apple Silicon
 * Apple Silicon prefers 128-byte alignment for optimal cache performance
 */
#define NEON_ALIGN __attribute__((aligned(128)))
#define NEON_CACHE_LINE 128

static inline void* neon_aligned_alloc(size_t size) {
    void* ptr = NULL;
    if (posix_memalign(&ptr, NEON_CACHE_LINE, size) != 0) {
        return NULL;
    }
    return ptr;
}

static inline void neon_aligned_free(void* ptr) {
    free(ptr);
}

#endif /* __ARM_NEON || __aarch64__ || APPLE_SILICON */

#endif /* NEON_UTILS_H */
