#include "meth_scoring.h"

#include <cstdlib>

#include "bwa.h"      // bwa_fill_scmat
#include "bwamem.h"   // mem_opt_init, mem_opt_fill_meth_mat

namespace bwa_tests {

ScoringMatrix meth_scoring_matrix(int scoring, bool ot, int match, int mismatch) {
    mem_opt_t *o = mem_opt_init();
    o->a = match;
    o->b = mismatch;
    o->meth_scoring = scoring;
    bwa_fill_scmat(o->a, o->b, o->mat);
    mem_opt_fill_meth_mat(o);
    const int8_t *src = ot ? o->mat_ot : o->mat_ob;
    ScoringMatrix mat;
    for (int i = 0; i < 25; i++) mat[i] = src[i];
    free(o);
    return mat;
}

} // namespace bwa_tests
