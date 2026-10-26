/*
 * seq_reverse.h: in-place byte-array reverse.
 *
 * Between the two kswv passes of mate rescue, mem_sam_pe_batch_run reverses the
 * aligned prefix of each surviving query and reference window in place. For
 * the plain swap loop, clang 19.1.7 at -O3 emits an unrolled scalar byte swap
 * (ldrb/strb pairs, no vector code). A Graviton4 cycle profile (5M WGS pairs,
 * -t 16) put those two inlined loops at ~0.48% of total cycles.
 *
 * seq_reverse_inplace swaps 16-byte blocks from both ends, each block reversed
 * with one shuffle (NEON vext+vrev64, SSSE3 pshufb), then 8-byte blocks
 * (bswap64), then single bytes. Every output byte is
 * s_out[k] = s_in[n - 1 - k], exactly what the scalar swap loop produces.
 *
 * When fewer than two blocks' worth remains, the last block pair overlaps in
 * the middle. That is safe: both blocks are loaded before either is stored,
 * and each store writes, at every position it covers, the byte that belongs
 * there in the reversed output, so where the two stores overlap they write the
 * same value. The span left after an overlapping pair is negative, which ends
 * the remaining loops.
 */

#ifndef BWAMEM3_SEQ_REVERSE_H
#define BWAMEM3_SEQ_REVERSE_H

#include <stdint.h>
#include <string.h>

#if defined(__aarch64__) || defined(__ARM_NEON)
#  include <arm_neon.h>
#  define BWA3_SEQ_REV_VEC16 1
#elif defined(__SSSE3__)
#  include <tmmintrin.h>
#  define BWA3_SEQ_REV_VEC16 2
#endif

#if defined(BWA3_SEQ_REV_VEC16) && BWA3_SEQ_REV_VEC16 == 1
typedef uint8x16_t seq_rev_v16;
static inline seq_rev_v16 seq_rev_ld16(const uint8_t *p) { return vld1q_u8(p); }
static inline void seq_rev_st16(uint8_t *p, seq_rev_v16 v) { vst1q_u8(p, v); }
/* Swap the 64-bit halves, then reverse the bytes within each half. A single
 * tbl against a {15..0} index (sam_encode.cpp's form) measured ~0.08% more
 * user CPU on Graviton4 here, so this keeps vext+vrev64. */
static inline seq_rev_v16 seq_rev_rev16(seq_rev_v16 v) { return vrev64q_u8(vextq_u8(v, v, 8)); }
#elif defined(BWA3_SEQ_REV_VEC16) && BWA3_SEQ_REV_VEC16 == 2
typedef __m128i seq_rev_v16;
static inline seq_rev_v16 seq_rev_ld16(const uint8_t *p) { return _mm_loadu_si128((const __m128i *)p); }
static inline void seq_rev_st16(uint8_t *p, seq_rev_v16 v) { _mm_storeu_si128((__m128i *)p, v); }
static inline seq_rev_v16 seq_rev_rev16(seq_rev_v16 v)
{
    return _mm_shuffle_epi8(v, _mm_set_epi8(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15));
}
#endif

static inline uint64_t seq_rev_ld8(const uint8_t *p) { uint64_t x; memcpy(&x, p, 8); return x; }
static inline void seq_rev_st8(uint8_t *p, uint64_t x) { memcpy(p, &x, 8); }

/* Reverse s[i, j) in place with 8-byte bswap blocks, then single bytes. A span
 * j - i < 2 (including a negative one left by an overlapping block pair) is a
 * no-op. Callers keep i, j >= 0. */
static inline void seq_reverse_range_u64(uint8_t *s, int i, int j)
{
    while (j - i >= 8) {
        const uint64_t L = seq_rev_ld8(s + i), R = seq_rev_ld8(s + j - 8);
        seq_rev_st8(s + i, __builtin_bswap64(R));
        seq_rev_st8(s + j - 8, __builtin_bswap64(L));
        i += 8; j -= 8;
    }
    for (--j; i < j; ++i, --j) { const uint8_t t = s[i]; s[i] = s[j]; s[j] = t; }
}

/* Reverse s[0, n) in place using only the 8-byte and single-byte steps; n < 2
 * is a no-op. This is the whole reverse on a target with neither NEON nor
 * SSSE3. Every shipped target has one of them, so this entry point exists for
 * the unit test, to exercise the fallback on every host. */
static inline void seq_reverse_inplace_u64(uint8_t *s, int n)
{
    if (n < 2) return;
    seq_reverse_range_u64(s, 0, n);
}

/* Reverse s[0, n) in place; n < 2 is a no-op. */
static inline void seq_reverse_inplace(uint8_t *s, int n)
{
    if (n < 2) return;
    int i = 0, j = n; /* s[i, j) is the part not yet reversed */
#ifdef BWA3_SEQ_REV_VEC16
    while (j - i >= 16) {
        const seq_rev_v16 L = seq_rev_ld16(s + i), R = seq_rev_ld16(s + j - 16);
        seq_rev_st16(s + i, seq_rev_rev16(R));
        seq_rev_st16(s + j - 16, seq_rev_rev16(L));
        i += 16; j -= 16;
    }
#endif
    seq_reverse_range_u64(s, i, j);
}

#endif /* BWAMEM3_SEQ_REVERSE_H */
