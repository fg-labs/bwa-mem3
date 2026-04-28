// test/framework/scoring.cpp
#include "scoring.h"

namespace bwa_tests {

ScoringMatrix build_scoring_matrix(int match, int mismatch, int ambig) {
    ScoringMatrix mat{};
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            mat[i * 5 + j] = static_cast<int8_t>(i == j ? match : -mismatch);
        }
        mat[i * 5 + 4] = static_cast<int8_t>(-ambig);
        mat[4 * 5 + i] = static_cast<int8_t>(-ambig);
    }
    mat[24] = static_cast<int8_t>(-ambig);
    return mat;
}

ScoringMatrix default_scoring_matrix() {
    return build_scoring_matrix(1, 4, 1);
}

} // namespace bwa_tests
