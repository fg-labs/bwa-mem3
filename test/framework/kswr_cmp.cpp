// test/framework/kswr_cmp.cpp
#include "kswr_cmp.h"

namespace bwa_tests {

bool kswr_score_eq(const kswr_t &a, const kswr_t &b) {
    return a.score == b.score;
}

bool kswr_coords_eq(const kswr_t &scalar, const kswr_t &batched,
                    int min_seed_len) {
    if (scalar.score == 0) return true;
    if (scalar.score < min_seed_len) return true;
    return scalar.qb == batched.qb && scalar.tb == batched.tb;
}

bool kswr_score2_eq(const kswr_t &scalar, const kswr_t &batched) {
    if (scalar.score == 0) return true;
    int a = (scalar.score2  <= 0) ? -1 : scalar.score2;
    int b = (batched.score2 <= 0) ? -1 : batched.score2;
    return a == b;
}

} // namespace bwa_tests
