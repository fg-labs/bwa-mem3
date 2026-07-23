// test/unit/test_alnreg_sort_dedup.cpp — mem_sort_dedup_patch's two sort modes.
//
// mem_sort_dedup_patch opens with a sort by reference END position whose result
// feeds an ORDER-SENSITIVE dedup loop, so the choice of comparator decides which
// of several equal-`re` regions survives:
//
//   opt->alnreg_sort_fast == 0 (default, set by mem_opt_init)
//       bwa-mem2's comparator -- a PARTIAL order on `re` alone -- plus klib's
//       unstable ks_introsort. Equal-`re` records keep whatever permutation
//       introsort leaves them in, which is exactly what bwa-mem2's output is
//       defined by. Consequence: the survivor depends on the input permutation.
//
//   opt->alnreg_sort_fast == 1 (set by --fast)
//       a STRICT TOTAL order (re, rb, score desc, qb) plus pdqsort. Equal-`re`
//       records are fully ordered, so the survivor is permutation-independent --
//       and necessarily differs from bwa-mem2's on some ties.
//
// The two properties are mutually exclusive, which is why the lever exists. The
// cases below pin both halves: divergence on a tie, agreement when there is no
// tie to resolve, and the permutation (in)dependence that distinguishes them.
//
// Every expectation here is derived from the comparators themselves, not
// recorded from a run:
//   * ks_introsort's n == 2 fast path is a single compare-swap, and the re-only
//     comparator reports `false` for equal `re`, so a 2-element equal-`re` input
//     is returned in its original order.
//   * The strict total order breaks an equal-`re` tie on `rb` ascending.
//   * In the dedup loop's redundancy branch, `p->score < q->score` is false when
//     the scores are equal, so the EARLIER record (q) is dropped -- i.e. the
//     region the first sort placed LAST survives.
//
// Fixtures are built in-memory; no test data files are read.

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <set>
#include <vector>

#include "doctest/doctest.h"
#include "bwamem.h"

// Defined in src/bwamem.cpp. `mat` defaults to NULL there via a declaration in
// src/bwamem_pair.cpp; spell it out here so this TU introduces no new default.
extern int mem_sort_dedup_patch(const mem_opt_t *opt, const bntseq_t *bns,
                                const uint8_t *pac, uint8_t *query, int n,
                                mem_alnreg_t *a, const int8_t *mat);

namespace {

// Region end shared by every tie fixture; the exact value is irrelevant, only
// that the tied regions agree on it.
const int64_t TIE_RE = 1150;

// Build the subset of mem_alnreg_t fields mem_sort_dedup_patch reads: the
// [rb,re) / [qb,qe) spans, rid and score. Everything else is zeroed.
mem_alnreg_t reg(int64_t rb, int64_t re, int qb, int qe, int score, int rid = 0) {
    mem_alnreg_t a;
    memset(&a, 0, sizeof(a));
    a.rb = rb; a.re = re;
    a.qb = qb; a.qe = qe;
    a.rid = rid;
    a.score = a.truesc = score;
    a.meth_hypothesis = a.meth_strand_hyp = -1;
    return a;
}

// A tie member: [TIE_RE - len, TIE_RE) on the reference, [150 - len, 150) on the
// query. Every pair of these is mutually redundant under mask_level_redun = 0.95
// (the reference overlap is the shorter region's full length and the query
// overlap equals the shorter query span), so the dedup loop always takes its
// redundancy branch and mem_patch_reg -- which would dereference the NULL bns we
// pass -- is never reached.
mem_alnreg_t tie_reg(int len, int score) {
    return reg(TIE_RE - len, TIE_RE, 150 - len, 150, score);
}

// Run mem_sort_dedup_patch over a copy of `in` with the given sort mode and
// return the surviving regions.
std::vector<mem_alnreg_t> dedup(const std::vector<mem_alnreg_t> &in, int sort_fast) {
    mem_opt_t *opt = mem_opt_init();
    opt->alnreg_sort_fast = sort_fast;
    std::vector<mem_alnreg_t> a(in);
    int n = mem_sort_dedup_patch(opt, NULL, NULL, NULL,
                                 static_cast<int>(a.size()), a.data(), NULL);
    free(opt);
    REQUIRE(n >= 0);
    REQUIRE(n <= static_cast<int>(a.size()));
    a.resize(static_cast<size_t>(n));
    return a;
}

// The `rb` of every survivor, in output order.
std::vector<int64_t> rbs(const std::vector<mem_alnreg_t> &a) {
    std::vector<int64_t> v;
    for (size_t i = 0; i < a.size(); ++i) v.push_back(a[i].rb);
    return v;
}

// The closing sort's key: score descending, then rb ascending, then qb ascending
// (alnreg_slt in src/bwamem.cpp). mem_mark_primary_se / mem_pair rely on it, so
// it must hold in BOTH modes.
bool ordered_by_score(const std::vector<mem_alnreg_t> &a) {
    for (size_t i = 1; i < a.size(); ++i) {
        const mem_alnreg_t &p = a[i - 1], &q = a[i];
        bool ok = p.score > q.score ||
                  (p.score == q.score &&
                   (p.rb < q.rb || (p.rb == q.rb && p.qb <= q.qb)));
        if (!ok) return false;
    }
    return true;
}

}  // namespace

// ---------------------------------------------------------------------------
// Default contract
// ---------------------------------------------------------------------------

TEST_CASE("mem_opt_init defaults to the bwa-mem2-compatible dedup sort"
          * doctest::test_suite("unit/alnreg_sort_dedup")) {
    mem_opt_t *opt = mem_opt_init();
    CHECK(opt->alnreg_sort_fast == 0);
    free(opt);
}

TEST_CASE("n <= 1 short-circuits identically in both modes"
          * doctest::test_suite("unit/alnreg_sort_dedup")) {
    for (int fast = 0; fast <= 1; ++fast) {
        CAPTURE(fast);
        CHECK(dedup(std::vector<mem_alnreg_t>(), fast).empty());

        std::vector<mem_alnreg_t> one(1, tie_reg(150, 150));
        std::vector<mem_alnreg_t> out = dedup(one, fast);
        REQUIRE(out.size() == 1u);
        CHECK(out[0].rb == one[0].rb);
    }
}

// ---------------------------------------------------------------------------
// Equal-`re` tie: the modes pick different survivors
// ---------------------------------------------------------------------------

TEST_CASE("equal-`re` tie: default keeps bwa-mem2's input-order outcome, --fast does not"
          * doctest::test_suite("unit/alnreg_sort_dedup")) {
    // Two mutually redundant regions ending at TIE_RE with equal scores. `wide`
    // starts earlier (smaller rb), `narrow` later (larger rb).
    mem_alnreg_t wide   = tie_reg(150, 150);  // rb = 1000
    mem_alnreg_t narrow = tie_reg(148, 150);  // rb = 1002
    REQUIRE(wide.rb < narrow.rb);
    REQUIRE(wide.re == narrow.re);
    REQUIRE(wide.score == narrow.score);

    SUBCASE("input already in strict-total-order (narrow last): modes agree") {
        std::vector<mem_alnreg_t> in;
        in.push_back(wide); in.push_back(narrow);
        // Default: equal `re` => no swap => narrow stays last => narrow survives.
        // Fast: rb ascending puts narrow last => narrow survives. Same answer.
        CHECK(rbs(dedup(in, 0)) == std::vector<int64_t>(1, narrow.rb));
        CHECK(rbs(dedup(in, 1)) == std::vector<int64_t>(1, narrow.rb));
    }

    SUBCASE("input reversed (narrow first): the modes diverge") {
        std::vector<mem_alnreg_t> in;
        in.push_back(narrow); in.push_back(wide);
        // Default preserves the given order, so `wide` is last and survives --
        // this is the bwa-mem2 outcome the PR restores.
        CHECK(rbs(dedup(in, 0)) == std::vector<int64_t>(1, wide.rb));
        // Fast re-sorts on the total order, so `narrow` is last and survives.
        CHECK(rbs(dedup(in, 1)) == std::vector<int64_t>(1, narrow.rb));
    }
}

// ---------------------------------------------------------------------------
// Permutation (in)dependence over a larger tie cluster
// ---------------------------------------------------------------------------

TEST_CASE("tie cluster: --fast is permutation-independent, the default is not"
          * doctest::test_suite("unit/alnreg_sort_dedup")) {
    // Six equal-`re`, equal-score, mutually redundant regions. All but one are
    // dropped, so each permutation yields exactly one survivor.
    std::vector<mem_alnreg_t> base;
    for (int len = 145; len <= 150; ++len) base.push_back(tie_reg(len, 150));
    const int64_t max_rb = TIE_RE - 145;  // the largest rb in the cluster

    std::vector<mem_alnreg_t> perm(base);
    std::sort(perm.begin(), perm.end(),
              [](const mem_alnreg_t &x, const mem_alnreg_t &y) { return x.rb < y.rb; });

    std::set<int64_t> default_survivors, fast_survivors;
    do {
        std::vector<mem_alnreg_t> slow_out = dedup(perm, 0);
        std::vector<mem_alnreg_t> fast_out = dedup(perm, 1);
        REQUIRE(slow_out.size() == 1u);
        REQUIRE(fast_out.size() == 1u);
        default_survivors.insert(slow_out[0].rb);
        fast_survivors.insert(fast_out[0].rb);
    } while (std::next_permutation(
                 perm.begin(), perm.end(),
                 [](const mem_alnreg_t &x, const mem_alnreg_t &y) { return x.rb < y.rb; }));

    // The strict total order puts the largest `rb` last among equal `re`, and the
    // redundancy branch keeps the last region, for every one of the 720 inputs.
    CHECK(fast_survivors == std::set<int64_t>{max_rb});
    // The re-only comparator leaves ties to the sort's permutation, so the
    // survivor moves with the input -- the property bwa-mem2's output depends on.
    CHECK(default_survivors.size() > 1u);
}

// ---------------------------------------------------------------------------
// No ties to resolve: the modes must agree
// ---------------------------------------------------------------------------

TEST_CASE("well-separated regions survive intact and identically in both modes"
          * doctest::test_suite("unit/alnreg_sort_dedup")) {
    // Spaced far beyond max_chain_gap (10000), so the dedup loop's window check
    // skips every pair and nothing is merged or dropped. What remains is the
    // closing by-score sort, which downstream primary selection depends on.
    std::vector<mem_alnreg_t> in;
    in.push_back(reg(5000000, 5000150, 0, 150, 120));
    in.push_back(reg(1000000, 1000150, 0, 150, 150));
    in.push_back(reg(3000000, 3000150, 0, 150, 90));
    in.push_back(reg(7000000, 7000150, 0, 150, 150));

    std::vector<mem_alnreg_t> slow_out = dedup(in, 0);
    std::vector<mem_alnreg_t> fast_out = dedup(in, 1);

    REQUIRE(slow_out.size() == in.size());
    REQUIRE(fast_out.size() == in.size());
    CHECK(ordered_by_score(slow_out));
    CHECK(ordered_by_score(fast_out));
    // Score desc, then rb asc: 150@1000000, 150@7000000, 120@5000000, 90@3000000.
    std::vector<int64_t> want;
    want.push_back(1000000); want.push_back(7000000);
    want.push_back(5000000); want.push_back(3000000);
    CHECK(rbs(slow_out) == want);
    CHECK(rbs(fast_out) == want);
}

TEST_CASE("distinct-score redundant regions: the best scoring one wins in both modes"
          * doctest::test_suite("unit/alnreg_sort_dedup")) {
    // Same tie cluster geometry, but the scores are distinct, so the redundancy
    // branch decides on score rather than on position -- no tie for the two
    // comparators to disagree about.
    std::vector<mem_alnreg_t> in;
    in.push_back(tie_reg(150, 130));
    in.push_back(tie_reg(148, 170));  // the unique best score
    in.push_back(tie_reg(146, 110));
    const int64_t best_rb = TIE_RE - 148;

    CHECK(rbs(dedup(in, 0)) == std::vector<int64_t>(1, best_rb));
    CHECK(rbs(dedup(in, 1)) == std::vector<int64_t>(1, best_rb));
}
