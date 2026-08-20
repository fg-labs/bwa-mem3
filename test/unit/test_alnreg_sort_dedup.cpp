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

// ---------------------------------------------------------------------------
// The tie-free fast path must be invisible
// ---------------------------------------------------------------------------
//
// In the default mode each of mem_sort_dedup_patch's two sorts runs pdqsort
// whenever the sorted array has no tied adjacent pair -- the ordering is then
// unique, so every correct comparison sort agrees -- and otherwise restores the
// input and runs ks_introsort over it. Both regimes have to be indistinguishable
// from running ks_introsort unconditionally. That is the whole correctness
// claim, so it is asserted directly on the sorts rather than inferred from the
// dedup output.

// src/bwamem.cpp: each fast path and the exact introsort it must reproduce.
extern void bwamem3_dedup_sort_by_re(int n, mem_alnreg_t *a);
extern void bwamem3_dedup_sort_by_re_exact(int n, mem_alnreg_t *a);
extern void bwamem3_dedup_sort_by_score(int n, mem_alnreg_t *a);
extern void bwamem3_dedup_sort_by_score_exact(int n, mem_alnreg_t *a);

// The permutation-gather sorts mem_sort_dedup_patch uses on the default path:
// they sort a (key,index) permutation and gather the 112-byte records once,
// dropping the per-call save-copy the tie-detection scheme otherwise pays. They
// MUST reproduce the exact same array (permutation included) as the unconditional
// ks_introsort oracle -- the whole correctness claim -- so they run through the
// identical `agrees()` harness. Thin in-place wrappers (perm into a thread-local
// buffer, then copy back) so the existing sort_fn harness can drive them; the
// hot path calls the src->dst core directly and never copies back.
extern void bwamem3_dedup_perm_sort_by_re(int n, mem_alnreg_t *a);
extern void bwamem3_dedup_perm_sort_by_score(int n, mem_alnreg_t *a);

namespace {

// Deterministic 64-bit LCG. std::rand would make the fixtures platform-defined.
struct Rng {
    uint64_t s;  // LCG state; seeded per test so every fixture is reproducible.
    explicit Rng(uint64_t seed) : s(seed) {}
    // Advances the state and returns the high bits, whose period is the long one.
    uint64_t next() { s = s * 6364136223846793005ULL + 1442695040888963407ULL; return s >> 11; }
    // A draw from the inclusive range [lo, hi]; the modulo bias is irrelevant here.
    int in(int lo, int hi) { return lo + static_cast<int>(next() % static_cast<uint64_t>(hi - lo + 1)); }
};

// `n` regions on one reference. re_pool > 0 draws the reference end from that
// many slots, so equal-`re` ties (the fallback) are common; re_pool == 0 gives
// every region a distinct end, so the `re` sort is tie-free (the fast path).
// Scores are drawn from a small range so the closing sort sees ties too.
std::vector<mem_alnreg_t> random_regs(Rng &rng, int n, int re_pool) {
    std::vector<int64_t> ends;
    for (int i = 0; i < n; ++i)
        ends.push_back(1000 + 37 * (re_pool > 0 ? rng.in(0, re_pool - 1) : i));
    for (int i = n - 1; i > 0; --i) std::swap(ends[i], ends[rng.in(0, i)]);  // shuffle

    std::vector<mem_alnreg_t> v;
    for (int i = 0; i < n; ++i) {
        const int len = rng.in(20, 150);
        const int qb = rng.in(0, 150 - len);
        v.push_back(reg(ends[i] - len, ends[i], qb, qb + len, rng.in(1, 12)));
    }
    return v;
}

// Byte-identity over the whole record, not just the sort key: the fast path has
// to reproduce introsort's *permutation*, so two arrays that merely compare equal
// under the comparator are not good enough.
bool same_records(const std::vector<mem_alnreg_t> &x, const std::vector<mem_alnreg_t> &y) {
    if (x.size() != y.size()) return false;
    for (size_t i = 0; i < x.size(); ++i)
        if (memcmp(&x[i], &y[i], sizeof(mem_alnreg_t)) != 0) return false;
    return true;
}

// Sort one fixture both ways and report whether they agree; `tied` reports
// whether the sorted array actually contained a tie, so a test can prove it
// exercised the fallback instead of passing vacuously.
typedef void (*sort_fn)(int, mem_alnreg_t *);
bool agrees(const std::vector<mem_alnreg_t> &in, sort_fn fast, sort_fn exact,
            bool (*tied)(const mem_alnreg_t &, const mem_alnreg_t &), bool *saw_tie) {
    std::vector<mem_alnreg_t> f(in), e(in);
    fast(static_cast<int>(f.size()), f.data());
    exact(static_cast<int>(e.size()), e.data());
    for (size_t i = 1; i < e.size(); ++i)
        if (tied(e[i - 1], e[i])) { *saw_tie = true; break; }
    return same_records(f, e);
}

// The two `tied` predicates `agrees` takes: each says when its comparator ranks
// the pair equally, i.e. when the fast path's tie scan would fire.

// alnreg_slt2_m2 orders on `re` alone, so equal ends are a tie.
bool tied_by_re(const mem_alnreg_t &x, const mem_alnreg_t &y) { return x.re == y.re; }
// alnreg_slt breaks score ties on (rb, qb), so all three must match.
bool tied_by_score(const mem_alnreg_t &x, const mem_alnreg_t &y) {
    return x.score == y.score && x.rb == y.rb && x.qb == y.qb;
}

}  // namespace

TEST_CASE("the `re` fast path is byte-identical to unconditional ks_introsort"
          * doctest::test_suite("unit/alnreg_sort_dedup")) {
    const int TRIALS = 150;

    SUBCASE("tie-free arrays (the pdqsort path itself)") {
        Rng rng(0x51ed0001ULL);
        bool saw_tie = false;
        for (int t = 0; t < TRIALS; ++t) {
            // n spans both sides of the n >= 9 gate.
            CHECK(agrees(random_regs(rng, 1 + t, 0), bwamem3_dedup_sort_by_re,
                         bwamem3_dedup_sort_by_re_exact, tied_by_re, &saw_tie));
        }
        CHECK_FALSE(saw_tie);  // distinct ends by construction
    }

    SUBCASE("tie-heavy arrays (the restore-and-introsort fallback)") {
        Rng rng(0x51ed0002ULL);
        bool saw_tie = false;
        for (int t = 0; t < TRIALS; ++t) {
            CHECK(agrees(random_regs(rng, 1 + t, 3), bwamem3_dedup_sort_by_re,
                         bwamem3_dedup_sort_by_re_exact, tied_by_re, &saw_tie));
        }
        CHECK(saw_tie);  // guards against a vacuous pass
    }

    SUBCASE("mixed arrays") {
        Rng rng(0x51ed0003ULL);
        bool saw_tie = false;
        for (int t = 0; t < TRIALS; ++t) {
            const int n = 1 + t;
            CHECK(agrees(random_regs(rng, n, n / 2 + 1), bwamem3_dedup_sort_by_re,
                         bwamem3_dedup_sort_by_re_exact, tied_by_re, &saw_tie));
        }
        CHECK(saw_tie);
    }
}

TEST_CASE("the by-score fast path is byte-identical to unconditional ks_introsort"
          * doctest::test_suite("unit/alnreg_sort_dedup")) {
    // Real inputs never tie here -- two survivors cannot share (rb, qb) -- but
    // the fast path must not depend on that, so feed it ties anyway. Duplicated
    // regions give equal (score, rb, qb).
    Rng rng(0x51ed0004ULL);
    bool saw_tie = false;
    for (int t = 0; t < 150; ++t) {
        std::vector<mem_alnreg_t> in = random_regs(rng, 1 + t, 4);
        for (size_t i = 1; i < in.size(); i += 3) in[i] = in[i - 1];  // force exact duplicates
        CHECK(agrees(in, bwamem3_dedup_sort_by_score,
                     bwamem3_dedup_sort_by_score_exact, tied_by_score, &saw_tie));
    }
    CHECK(saw_tie);
}

// ---------------------------------------------------------------------------
// The permutation-gather sorts must be byte-identical to the same oracle
// ---------------------------------------------------------------------------
//
// Same contract as the fast paths above, and tested the same way: not merely
// "sorted" but the SAME ARRAY ks_introsort produces, element for element. The
// permutation is what makes this non-trivial -- the (key,index) pairs are sorted
// with klib's OWN ks_introsort, so on tied keys the gathered order can only match
// if the two sorts applied the same permutation. These wrappers have no
// DEDUP_PERM_MIN gate (that gate lives in mem_sort_dedup_patch), so the n-sweep
// drives the permutation core directly at every n in [1, 150].

TEST_CASE("the `re` permutation sort is byte-identical to unconditional ks_introsort"
          * doctest::test_suite("unit/alnreg_sort_dedup")) {
    const int TRIALS = 150;

    SUBCASE("tie-free arrays") {
        Rng rng(0x9e370001ULL);
        bool saw_tie = false;
        for (int t = 0; t < TRIALS; ++t)
            CHECK(agrees(random_regs(rng, 1 + t, 0), bwamem3_dedup_perm_sort_by_re,
                         bwamem3_dedup_sort_by_re_exact, tied_by_re, &saw_tie));
        CHECK_FALSE(saw_tie);
    }
    SUBCASE("tie-heavy arrays (the introsort-on-pairs permutation)") {
        Rng rng(0x9e370002ULL);
        bool saw_tie = false;
        for (int t = 0; t < TRIALS; ++t)
            CHECK(agrees(random_regs(rng, 1 + t, 3), bwamem3_dedup_perm_sort_by_re,
                         bwamem3_dedup_sort_by_re_exact, tied_by_re, &saw_tie));
        CHECK(saw_tie);
    }
    SUBCASE("mixed arrays") {
        Rng rng(0x9e370003ULL);
        bool saw_tie = false;
        for (int t = 0; t < TRIALS; ++t) {
            const int n = 1 + t;
            CHECK(agrees(random_regs(rng, n, n / 2 + 1), bwamem3_dedup_perm_sort_by_re,
                         bwamem3_dedup_sort_by_re_exact, tied_by_re, &saw_tie));
        }
        CHECK(saw_tie);
    }
}

TEST_CASE("the by-score permutation sort is byte-identical to unconditional ks_introsort"
          * doctest::test_suite("unit/alnreg_sort_dedup")) {
    Rng rng(0x9e370004ULL);
    bool saw_tie = false;
    for (int t = 0; t < 150; ++t) {
        std::vector<mem_alnreg_t> in = random_regs(rng, 1 + t, 4);
        // Force a (score, rb, qb) tie, then perturb a NON-key field (seedcov) so
        // the two records stay distinguishable by memcmp. Byte-identical dupes
        // pass under either permutation, so the tie would prove nothing; with a
        // distinct seedcov a wrong permutation on the tied pair is observable.
        for (size_t i = 1; i < in.size(); i += 3) {
            in[i] = in[i - 1];
            in[i].seedcov = in[i - 1].seedcov + 1;
        }
        CHECK(agrees(in, bwamem3_dedup_perm_sort_by_score,
                     bwamem3_dedup_sort_by_score_exact, tied_by_score, &saw_tie));
    }
    CHECK(saw_tie);
}
