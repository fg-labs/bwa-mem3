// test/integration/test_dedup_n_comp.cpp
//
// n_comp saturation through the real join path. mem_sort_dedup_patch's window
// pass (src/bwamem.cpp) resets every region's n_comp to 1 and then counts each
// join as the two counts plus one (mem_alnreg_add_comp), so k pieces joined
// into one region count 2k - 1. The count can only reach the 16-bit field's
// limit, MEM_N_COMP_MAX (32,767), through real joins: 16,384 of them. A join
// site that adds without the saturating helper wraps past it to a negative
// value, and this case fails.
//
// It is here rather than with the helper's own unit tests
// (test/unit/test_mem_extend_patch.cpp) because it cannot be made fast: every
// join re-aligns the whole span joined so far, so the pass is quadratic in the
// piece count, and its piece count is fixed by the limit (about 160 ms, past
// the 100 ms unit budget). The reference is random bases in a one-contig
// meth_test::OrigRefFixture, and the read is an exact copy of its first k bases.

#include "doctest/doctest.h"

#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

#include "bntseq.h"
#include "bwa.h"
#include "bwamem.h"
#include "../unit/meth_orig_ref_fixture.h"

// Non-static in bwamem.cpp but not declared in bwamem.h.
int mem_sort_dedup_patch(const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac, uint8_t *query, int n,
                         mem_alnreg_t *a, const int8_t *mat);

TEST_CASE("mem_sort_dedup_patch: a join of more pieces than n_comp holds saturates"
          * doctest::test_suite("integration/alnreg_sort_dedup")) {
    // 16,385 pieces would count 2k - 1 = 32,769. Each one-base piece of the
    // exact-match read joins the run before it, in reference-end order (w = r =
    // 0, and the equal-length global alignment is an ungapped scan).
    const int k = MEM_N_COMP_MAX / 2 + 2;
    const std::string ref = meth_test::random_bases(k + 16, 27182u);
    meth_test::OrigRefFixture fx(ref);
    mem_opt_t *opt = mem_opt_init();
    std::vector<uint8_t> q((size_t) k);
    for (int i = 0; i < k; ++i) q[(size_t) i] = meth_test::nt4(ref[(size_t) i]);
    std::vector<mem_alnreg_t> regs((size_t) k);
    for (int i = 0; i < k; ++i) {
        mem_alnreg_t &r = regs[(size_t) i];
        memset(&r, 0, sizeof r);
        r.qb = i, r.qe = i + 1, r.rb = i, r.re = i + 1, r.score = r.truesc = 1;
        r.secondary = -1;
    }
    const int n = mem_sort_dedup_patch(opt, fx.orig_bns, fx.orig_pac, q.data(), k, regs.data(), opt->mat);
    free(opt);
    REQUIRE(n == 1);
    CHECK(regs[0].qb == 0);
    CHECK(regs[0].qe == k);
    CHECK(regs[0].score == k);
    CHECK((int) regs[0].n_comp == MEM_N_COMP_MAX);
}
