// test/framework/ext_result.h
//
// The six result fields of one seed extension -- what scalarBandedSWA and
// ksw_extend2 return through their out-parameters and what getScores8 /
// getScores16 write back into a SeqPair. Shared by the extension-ladder model,
// the batch scorer and the bandedSWA kernel tests that hold a scalar oracle
// per pair.

#ifndef BWA_TESTS_EXT_RESULT_H
#define BWA_TESTS_EXT_RESULT_H

namespace bwa_tests {

struct ExtResult {
    int score, qle, tle, gtle, gscore, max_off;
    bool operator==(const ExtResult &o) const {
        return score == o.score && qle == o.qle && tle == o.tle && gtle == o.gtle &&
               gscore == o.gscore && max_off == o.max_off;
    }
};

} // namespace bwa_tests

#endif
