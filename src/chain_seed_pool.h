#ifndef CHAIN_SEED_POOL_H
#define CHAIN_SEED_POOL_H

#include <stdint.h>
#include <stdlib.h>
#include <vector>
#include "bwamem.h"   // mem_chain_t, mem_seed_t
#include "macro.h"    // SEEDS_PER_CHAIN
#include "utils.h"    // xassert

/* The seeds of the read being chained (mem_chain_seeds), each chain's a list
 * linked by index in merge order: mem_chain_t::link holds the indices of its
 * first and last seed. A merge (test_and_merge) appends in O(1) with no
 * allocation; chain_pool_drain copies a finished chain's seeds, contiguous, to
 * the seed window, or to a heap array once the window is full. One per thread,
 * cleared at the start of each read and when the flat chain index hands a read
 * back to the kbtree. The vectors keep their capacity across reads up to
 * KEEP_SEEDS seeds; a read that grew them past that releases the memory at the
 * next clear(), so one seed-heavy read does not pin it for the thread's life. */
struct ChainSeedPool {
    static constexpr size_t KEEP_SEEDS = (size_t) 1 << 18;   /* ~9 MB of seed + next */
    std::vector<mem_seed_t> seed;
    std::vector<int32_t> next;   /* -1 ends a chain's list */

    void clear() {
        if (seed.capacity() > KEEP_SEEDS) {
            std::vector<mem_seed_t>().swap(seed);
            std::vector<int32_t>().swap(next);
        } else {
            seed.clear();
            next.clear();
        }
    }
    /* Add s as a list of its own; its index. */
    inline int32_t add(const mem_seed_t *s) {
        xassert(seed.size() < (size_t) INT32_MAX, "chain seed pool index overflow");
        const int32_t i = (int32_t) seed.size();
        seed.push_back(*s);
        next.push_back(-1);
        return i;
    }
    /* Open chain c on seed s. */
    inline void start(mem_chain_t *c, const mem_seed_t *s) {
        c->n = 1;
        c->link.head = c->link.tail = add(s);
    }
    /* Append seed s to chain c. */
    inline void append(mem_chain_t *c, const mem_seed_t *s) {
        const int32_t i = add(s);
        next[c->link.tail] = i;
        c->link.tail = i;
        c->n++;
    }
};

/* Chain c of the read in `pool`, with its seeds copied contiguous and in merge
 * order into the per-thread seed window (seedBuf, *seedBufCount slots of
 * seedBufSize used), or into a heap array of its own once the window is full.
 * m records which for the free sites: the window holds the array at
 * m == SEEDS_PER_CHAIN, the heap at m == SEEDS_PER_CHAIN + 1 (the ownership
 * marker every free site tests). m is not a capacity: nothing grows a chain
 * after this, and every seed read is bounded by n. */
static inline mem_chain_t chain_pool_drain(const mem_chain_t &c, const ChainSeedPool &pool,
                                           mem_seed_t *seedBuf, int64_t *seedBufCount,
                                           int64_t seedBufSize)
{
    mem_chain_t out = c;
    mem_seed_t *dst;
    if (*seedBufCount + c.n > seedBufSize) {
        out.m = SEEDS_PER_CHAIN + 1;
        dst = (mem_seed_t *) malloc((size_t) c.n * sizeof(mem_seed_t));
        xassert(dst != NULL, "out of memory: chain seeds");
    } else {
        out.m = SEEDS_PER_CHAIN;
        dst = seedBuf + *seedBufCount;
        *seedBufCount += c.n;
    }
    int32_t i = c.link.head;
    for (int k = 0; k < c.n; ++k) {
        dst[k] = pool.seed[i];
        i = pool.next[i];
    }
    /* n hops from head must end exactly past tail: a stale or corrupt link
     * would otherwise copy wrong seeds without any other symptom. */
    xassert(i == -1, "chain seed pool list does not end after n seeds");
    out.seeds = dst;
    return out;
}

#endif /* CHAIN_SEED_POOL_H */
