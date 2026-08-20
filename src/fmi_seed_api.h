#ifndef FMI_SEED_API_H
#define FMI_SEED_API_H

#include <stdint.h>

/* ---- Blessed public seeding surface for external consumers (e.g. minibwa). ----
 * Stable: changes here are API-breaking and must be versioned. bwa-mem3's own
 * code and this file share these definitions; do not fork them. */

#define CP_BLOCK_SIZE 64
#define CP_FILENAME_SUFFIX ".bwt.2bit.64"
#define CP_MASK 63
#define CP_SHIFT 6

typedef struct checkpoint_occ_scalar
{
    int64_t cp_count[4];
    uint64_t one_hot_bwt_str[4];
}CP_OCC;

#if defined(__clang__) || defined(__GNUC__)
static inline int _mm_countbits_64(unsigned long x) {
    return __builtin_popcountl(x);
}
#endif

/* One-hot position masks: entry [i] has the top `i` bits set (entry [0] == 0),
 * i.e. one_hot_mask_array[i] == (one_hot_mask_array[i-1] >> 1) | 0x8000...
 * for i in 1..63. The FMI-index lifetime constant is identical on every path,
 * so it lives inline here as a file-scope table rather than a heap allocation
 * reached through a per-object pointer: this removes a dependent load in the
 * (~10^9-call) backwardExt/GET_OCC hot path. Values are byte-identical to the
 * former runtime-computed array. */
static const uint64_t one_hot_mask_array[64] = {
    0x0000000000000000ULL, 0x8000000000000000ULL, 0xc000000000000000ULL, 0xe000000000000000ULL,
    0xf000000000000000ULL, 0xf800000000000000ULL, 0xfc00000000000000ULL, 0xfe00000000000000ULL,
    0xff00000000000000ULL, 0xff80000000000000ULL, 0xffc0000000000000ULL, 0xffe0000000000000ULL,
    0xfff0000000000000ULL, 0xfff8000000000000ULL, 0xfffc000000000000ULL, 0xfffe000000000000ULL,
    0xffff000000000000ULL, 0xffff800000000000ULL, 0xffffc00000000000ULL, 0xffffe00000000000ULL,
    0xfffff00000000000ULL, 0xfffff80000000000ULL, 0xfffffc0000000000ULL, 0xfffffe0000000000ULL,
    0xffffff0000000000ULL, 0xffffff8000000000ULL, 0xffffffc000000000ULL, 0xffffffe000000000ULL,
    0xfffffff000000000ULL, 0xfffffff800000000ULL, 0xfffffffc00000000ULL, 0xfffffffe00000000ULL,
    0xffffffff00000000ULL, 0xffffffff80000000ULL, 0xffffffffc0000000ULL, 0xffffffffe0000000ULL,
    0xfffffffff0000000ULL, 0xfffffffff8000000ULL, 0xfffffffffc000000ULL, 0xfffffffffe000000ULL,
    0xffffffffff000000ULL, 0xffffffffff800000ULL, 0xffffffffffc00000ULL, 0xffffffffffe00000ULL,
    0xfffffffffff00000ULL, 0xfffffffffff80000ULL, 0xfffffffffffc0000ULL, 0xfffffffffffe0000ULL,
    0xffffffffffff0000ULL, 0xffffffffffff8000ULL, 0xffffffffffffc000ULL, 0xffffffffffffe000ULL,
    0xfffffffffffff000ULL, 0xfffffffffffff800ULL, 0xfffffffffffffc00ULL, 0xfffffffffffffe00ULL,
    0xffffffffffffff00ULL, 0xffffffffffffff80ULL, 0xffffffffffffffc0ULL, 0xffffffffffffffe0ULL,
    0xfffffffffffffff0ULL, 0xfffffffffffffff8ULL, 0xfffffffffffffffcULL, 0xfffffffffffffffeULL,
};

typedef struct smem_struct
{
#ifdef DEBUG
    uint64_t info; // for debug
#endif
    uint32_t rid;
    uint32_t m, n;
    int64_t k, l, s;
}SMEM;

#endif /* FMI_SEED_API_H */
