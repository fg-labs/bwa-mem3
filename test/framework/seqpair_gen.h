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

// Reference is a random motif of 18-31 bases tiled to `rlen` with scattered
// substitutions (a tandem repeat); query is a mutated copy of a random window
// of that reference. Unlike gen_random_pair -- whose ref and query are
// independent, so the DP never rises far above zero -- this yields a high
// primary score PLUS many near-threshold secondary maxima against the other
// repeat copies. That plateau structure is what exercises the second-best
// (KSW_XSUBO) machinery: score2/te2 and the b[]-collapse it depends on.
// Asserts `qlen > 0 && rlen > qlen` — the query window is drawn from strictly
// inside the reference — so a mis-parameterized call fails loudly instead of
// sampling an empty window range. tag = "tandem_repeat".
TestPair gen_tandem_repeat_pair(std::mt19937 &rng, int qlen, int rlen);

} // namespace bwa_tests

#endif
