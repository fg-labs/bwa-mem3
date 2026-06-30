#ifndef SEED_ORDER_H
#define SEED_ORDER_H

#include <stdint.h>
#include "bwamem.h"   // mem_seed_t, seed_order_t

#ifdef __cplusplus
extern "C" {
#endif

// One resolved seed plus the metadata chaining needs. seed is the FULL
// mem_seed_t (all 7 fields), default-constructed then field-filled exactly as
// mem_chain_seeds does today. orig_ix is the strict resolve-order index and is
// every comparator's final stable tiebreak (and the OFF identity key).
typedef struct {
    mem_seed_t seed;
    int32_t    rid;
    int8_t     meth_hyp;
    uint32_t   orig_ix;
} seed_rec_t;

// Permute recs[0..n) in place per mode. OFF is a no-op. Stable on orig_ix.
void order_seeds(seed_rec_t *recs, int64_t n, seed_order_t mode);

// CLI parse: returns (seed_order_t)-1 on an unknown string.
seed_order_t seed_order_from_str(const char *s);
const char  *seed_order_to_str(seed_order_t m);

#ifdef __cplusplus
}
#endif
#endif
