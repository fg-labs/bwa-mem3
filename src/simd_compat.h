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

   SIMD Compatibility Header for Apple Silicon (ARM64/NEON)
   Based on PR #281 from BenjaminDEMAILLE
*****************************************************************************************/

#ifndef SIMD_COMPAT_H
#define SIMD_COMPAT_H

/*
 * Platform detection and SIMD abstraction layer
 *
 * On x86: Use native SSE/AVX intrinsics via immintrin.h
 * On ARM: Use sse2neon for SSE->NEON translation, plus native NEON for hot paths
 */

#if defined(__ARM_NEON) || defined(__aarch64__)
    /* ARM/Apple Silicon */
    #define APPLE_SILICON 1
    #define SIMDE_ENABLE_NATIVE_ALIASES
    #include "sse2neon.h"
    #include <arm_neon.h>

    /* Define SIMD widths for NEON (128-bit) */
    #ifndef SIMD_WIDTH8
    #define SIMD_WIDTH8 16   /* 128-bit / 8-bit = 16 elements */
    #endif
    #ifndef SIMD_WIDTH16
    #define SIMD_WIDTH16 8   /* 128-bit / 16-bit = 8 elements */
    #endif

    /* Memory allocation compatibility */
    #define CACHE_LINE_BYTES 128  /* Apple Silicon uses 128-byte cache lines */

    static inline void* _mm_malloc_compat(size_t size, size_t align) {
        void* ptr = NULL;
        if (posix_memalign(&ptr, align, size) != 0) {
            return NULL;
        }
        return ptr;
    }

    static inline void _mm_free_compat(void* ptr) {
        free(ptr);
    }

    /* Use compatibility functions on ARM
     * Apple Silicon uses 128-byte cache lines, so we enforce minimum 128-byte
     * alignment for all SIMD allocations (vs 64-byte on x86) to avoid false sharing */
    #define _mm_malloc(size, align) _mm_malloc_compat(size, (align) < CACHE_LINE_BYTES ? CACHE_LINE_BYTES : (align))
    #define _mm_free(ptr) _mm_free_compat(ptr)

    /* Prefetch compatibility */
    #ifndef _mm_prefetch
    #define _mm_prefetch(addr, hint) __builtin_prefetch((const void*)(addr), 0, (hint))
    #endif

    /* Mask types for SSE compatibility (sse2neon should provide these) */
    #ifndef __mmask8
    typedef uint8_t __mmask8;
    #endif
    #ifndef __mmask16
    typedef uint16_t __mmask16;
    #endif
    #ifndef __mmask32
    typedef uint32_t __mmask32;
    #endif
    #ifndef __mmask64
    typedef uint64_t __mmask64;
    #endif

    /* __rdtsc compatibility - sse2neon provides _rdtsc */
    #ifndef __rdtsc
    #define __rdtsc _rdtsc
    #endif

    /*
     * Optimized movemask for 16-bit elements (used heavily in bandedSWA.cpp)
     * Instead of _mm_movemask_epi8(v) & 0xAAAA, use _mm_movemask_epi16(v) directly.
     * This extracts the MSB of each 16-bit element into an 8-bit result.
     */
    static inline int _mm_movemask_epi16(__m128i v) {
        /* Shift right to get MSB in LSB position for each 16-bit element */
        uint16x8_t shifted = vshrq_n_u16(vreinterpretq_u16_m128i(v), 15);
        /* Narrow to 8-bit (taking low byte of each 16-bit element) */
        uint8x8_t narrow = vmovn_u16(shifted);
        /* Horizontal add with position weights: bit 0 has weight 1, bit 1 has weight 2, etc. */
        static const uint8_t weights[8] = {1, 2, 4, 8, 16, 32, 64, 128};
        uint8x8_t weighted = vmul_u8(narrow, vld1_u8(weights));
        return vaddv_u8(weighted);
    }

    /*
     * Optimized blendv for 16-bit elements using NEON vbsl (bitwise select).
     * This is more efficient than sse2neon's _mm_blendv_epi8 for 16-bit data.
     */
    static inline __m128i _mm_blendv_epi16_fast(__m128i x, __m128i y, __m128i mask) {
        /* Use vbsl: select y where mask bits are 1, else x */
        return vreinterpretq_m128i_s16(
            vbslq_s16(vreinterpretq_u16_m128i(mask),
                      vreinterpretq_s16_m128i(y),
                      vreinterpretq_s16_m128i(x)));
    }

#elif defined(__AVX512BW__)
    /* x86 with AVX-512 */
    #include <immintrin.h>
    #define CACHE_LINE_BYTES 64

#elif defined(__AVX2__)
    /* x86 with AVX2 */
    #include <immintrin.h>
    #define CACHE_LINE_BYTES 64

#elif defined(__SSE4_1__) || defined(__SSE2__)
    /* x86 with SSE */
    #include <smmintrin.h>
    #include <emmintrin.h>
    #define CACHE_LINE_BYTES 64
    #ifndef __mmask8
    typedef uint8_t __mmask8;
    #endif
    #ifndef __mmask16
    typedef uint16_t __mmask16;
    #endif

#else
    /* Scalar fallback */
    #define CACHE_LINE_BYTES 64
    #warning "No SIMD support detected, using scalar code paths"
#endif

/* Cross-platform aligned allocation macro */
#ifdef APPLE_SILICON
    #define SIMD_ALIGNED_ALLOC(size, align) _mm_malloc_compat(size, (align) < 128 ? 128 : (align))
    #define SIMD_ALIGNED_FREE(ptr) free(ptr)
#else
    #define SIMD_ALIGNED_ALLOC(size, align) _mm_malloc(size, align)
    #define SIMD_ALIGNED_FREE(ptr) _mm_free(ptr)
#endif

#endif /* SIMD_COMPAT_H */
