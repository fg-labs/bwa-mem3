/* chain_flat_key.h -- the packed sort key of the flat per-read chaining index (ChainFlat,
 * bwamem.cpp), and the parse of its chain-count cap (BWA3_CHAIN_FLAT_CAP).
 *
 * Each chain is one uint64 key, (pos << CHAIN_FLAT_CI_BITS) | ci: the chain's reference position
 * in the high bits and its index into the insertion-order chain array in the low bits. The flat
 * index keeps the keys ascending, so one sorted insert per new chain replaces the two parallel
 * inserts (positions and chain indices) it would otherwise need, and the predecessor probe reads
 * the position (tie test) and the chain index from the same element.
 *
 * Ordering: every chain index stays below 2^CHAIN_FLAT_CI_BITS (a read holds at most 2^20 chains
 * on the flat index) and 0 <= pos < CHAIN_FLAT_POS_LIMIT, so the shift neither overflows nor
 * lets ci spill into pos; the keys compare exactly as their positions do, with ci breaking only
 * equal positions. "Keys whose pos <= p" is then "keys <= pack(p, CI_MASK)". The caller keeps
 * the chain count at or below the cap, which chain_flat_parse_cap bounds by
 * CHAIN_FLAT_CI_LIMIT, and hands any position outside [0, CHAIN_FLAT_POS_LIMIT) to the kbtree,
 * so every key it builds satisfies both bounds.
 *
 * The constants are macros, so the header is valid C (C11, as macro.h needs) as well as C++. */
#ifndef BWA_MEM3_CHAIN_FLAT_KEY_H
#define BWA_MEM3_CHAIN_FLAT_KEY_H

#include <errno.h>
#include <limits.h>
#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>

#include "macro.h" /* BWA_UNPREDICTABLE */

#define CHAIN_FLAT_CI_BITS 20
#define CHAIN_FLAT_CI_LIMIT ((uint64_t)1 << CHAIN_FLAT_CI_BITS) /* 1048576 chains */
#define CHAIN_FLAT_CI_MASK (CHAIN_FLAT_CI_LIMIT - 1)
/* 2^(64 - 1 - CI_BITS) = 2^43: the key stays below 2^63, so it is also a valid non-negative
 * int64. */
#define CHAIN_FLAT_POS_LIMIT ((int64_t)1 << (63 - CHAIN_FLAT_CI_BITS))

/* The parsed chain-count cap: OK (the value, or the default when unset or empty), MALFORMED
 * (negative, non-numeric, trailing junk or a negative overflow: the default is used) or
 * CLAMPED (above CHAIN_FLAT_CI_LIMIT, including a positive strtol overflow: the limit is used). */
typedef enum {
    CHAIN_FLAT_CAP_OK = 0,
    CHAIN_FLAT_CAP_MALFORMED,
    CHAIN_FLAT_CAP_CLAMPED
} chain_flat_cap_status_t;
typedef struct {
    int cap;
    chain_flat_cap_status_t status;
} chain_flat_cap_t;

/* Parse a BWA3_CHAIN_FLAT_CAP value (NULL when unset) into the effective cap. Pure apart from
 * errno; the caller reports MALFORMED / CLAMPED. The result's cap is always in
 * [0, CHAIN_FLAT_CI_LIMIT] when default_cap is. */
static inline chain_flat_cap_t chain_flat_parse_cap(const char *e, int default_cap)
{
    chain_flat_cap_t r;
    r.cap = default_cap;
    r.status = CHAIN_FLAT_CAP_OK;
    if (e == NULL || e[0] == '\0') return r;
    char *end = NULL;
    errno = 0;
    const long v = strtol(e, &end, 10);
    const int overflow_hi = errno == ERANGE && v == LONG_MAX;
    if (end == e || *end != '\0' || (errno != 0 && !overflow_hi) || v < 0) {
        r.status = CHAIN_FLAT_CAP_MALFORMED;
        return r;
    }
    if (overflow_hi || (unsigned long)v > CHAIN_FLAT_CI_LIMIT) {
        r.cap = (int)CHAIN_FLAT_CI_LIMIT;
        r.status = CHAIN_FLAT_CAP_CLAMPED;
        return r;
    }
    r.cap = (int)v;
    return r;
}

static inline int chain_flat_pos_ok(int64_t pos)
{
    return (uint64_t)pos < (uint64_t)CHAIN_FLAT_POS_LIMIT;
}
/* Precondition: chain_flat_pos_ok(pos) and ci < CHAIN_FLAT_CI_LIMIT; otherwise the key is
 * meaningless. */
static inline uint64_t chain_flat_pack(int64_t pos, uint32_t ci)
{
    return ((uint64_t)pos << CHAIN_FLAT_CI_BITS) | (uint64_t)ci;
}
static inline int64_t chain_flat_pos(uint64_t key)
{
    return (int64_t)(key >> CHAIN_FLAT_CI_BITS);
}
static inline uint32_t chain_flat_ci(uint64_t key)
{
    return (uint32_t)(key & CHAIN_FLAT_CI_MASK);
}

/* Number of keys in the ascending keys[0, n) whose position is <= pos, i.e. the predecessor's
 * rank + 1. The step must stay a conditional move: the comparison is data-dependent (a coin flip
 * per level), and as a branch it mispredicted ~37% of the time on x86 -- clang's cmov-conversion
 * pass had turned the select back into a jump. BWA_UNPREDICTABLE (__builtin_unpredictable under
 * clang) keeps it a cmov there; other compilers get the plain select.
 * Precondition: chain_flat_pos_ok(pos); outside [0, CHAIN_FLAT_POS_LIMIT) the rank is
 * meaningless. */
static inline size_t chain_flat_rank_le(const uint64_t *keys, size_t n, int64_t pos)
{
    if (n == 0) return 0;
    const uint64_t q = chain_flat_pack(pos, (uint32_t)CHAIN_FLAT_CI_MASK);
    const uint64_t *b = keys;
    while (n > 1) { const size_t h = n >> 1; b = BWA_UNPREDICTABLE(b[h] <= q) ? b + h : b; n -= h; }
    return (size_t)(b - keys) + (*b <= q);
}

#endif /* BWA_MEM3_CHAIN_FLAT_KEY_H */
