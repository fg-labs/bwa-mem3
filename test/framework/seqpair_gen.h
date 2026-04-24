// test/framework/seqpair_gen.h
//
// Deterministic TestPair generators. Seed behavior is the caller's
// responsibility — pass in an already-seeded std::mt19937 for
// reproducibility across test runs.

#ifndef BWA_TESTS_SEQPAIR_GEN_H
#define BWA_TESTS_SEQPAIR_GEN_H

#include <cstdint>
#include <random>

#include "seqpair.h"

namespace bwa_tests {

// Random ACGT sequences; no N bases. tag = "random".
TestPair gen_random_pair(std::mt19937 &rng, int qlen, int rlen);

// Reference is 2x query length with a deterministic ACGT pattern; query is
// the middle half of the reference (exact match). tag = "exact_match".
TestPair gen_exact_match_pair(int len);

// Reference is all A (base 0), query is all C (base 1). tag = "all_mismatch".
TestPair gen_all_mismatch_pair(int len);

// Reference and query are the same homopolymer. base must be 0..3.
// tag = "homopolymer".
TestPair gen_homopolymer_pair(int len, uint8_t base);

// Random pair, then overwrite `run_len` query bases at `run_pos` with new
// random bases — a substitution cluster, NOT an indel (ref and query stay
// length-matched). Asserts `run_pos + run_len <= qlen` so a mis-parameterized
// call fails loudly instead of silently no-opping. tag = "sub_cluster".
TestPair gen_sub_cluster_pair(std::mt19937 &rng, int qlen, int rlen,
                              int run_pos, int run_len);

// Random pair, then inject N bases (base 4) with probability n_frac_pct/100
// at each position. tag = "with_n_bases".
TestPair gen_with_n_bases_pair(std::mt19937 &rng, int qlen, int rlen,
                               int n_frac_pct);

} // namespace bwa_tests

#endif
