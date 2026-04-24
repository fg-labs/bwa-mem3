// test/framework/seqpair.h
//
// Shared test-only pair type. We deliberately do NOT use src/bandedSWA.h's
// SeqPair here: SeqPair is a flat-buffer descriptor (offsets into seqBufRef
// / seqBufQer) used by the production kernel, while TestPair owns its bytes
// for convenience of parameterization. seqpair_batch.{h,cpp} packs
// TestPairs into the flat SeqPair layout when a kernel needs it.

#ifndef BWA_TESTS_SEQPAIR_H
#define BWA_TESTS_SEQPAIR_H

#include <cstdint>
#include <string>
#include <vector>

namespace bwa_tests {

// A single ref/query pair with a human-readable tag used for diagnostic
// output (CAPTURE / INFO). Sequences are in 2-bit encoding (0=A, 1=C,
// 2=G, 3=T, 4=N), matching what kswv and ksw_align2 consume.
struct TestPair {
    std::vector<uint8_t> ref;
    std::vector<uint8_t> qry;
    std::string          tag;
};

} // namespace bwa_tests

#endif
