// test/framework/wang_hash64.h
//
// bwa's hash_64 (Thomas Wang's 64-bit integer mix), written out for the unit
// tests' oracles rather than taken from src/utils.h, so that a change to the
// library's hash, which would reorder equal-score hits and move every pair
// tie-break in the output, is caught rather than followed.

#ifndef BWA_TESTS_WANG_HASH64_H
#define BWA_TESTS_WANG_HASH64_H

#include <cstdint>

namespace bwa_tests {

inline uint64_t wang_hash64(uint64_t key) {
    key += ~(key << 32);
    key ^= (key >> 22);
    key += ~(key << 13);
    key ^= (key >> 8);
    key += (key << 3);
    key ^= (key >> 15);
    key += ~(key << 27);
    key ^= (key >> 31);
    return key;
}

}  // namespace bwa_tests

#endif
