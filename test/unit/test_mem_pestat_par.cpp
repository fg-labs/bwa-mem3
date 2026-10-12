// test/unit/test_mem_pestat_par.cpp — the parallel insert-size count equals the
// serial one.
//
// mem_pestat_par runs mem_pestat's candidate loop as a kt_for() pass: each
// worker counts the pairs of its own items into a per-thread histogram, and the
// histograms are summed before the unchanged percentile / mean / std.dev tail.
// The counts are integers, so the result must equal mem_pestat's exactly -- the
// same pes[] down to the bits of avg and std -- for any thread count and any
// way kt_for happens to cut the pairs into items.
//
// The oracle is the serial mem_pestat on the same input, plus one hand-built
// distribution whose answer is known without running either. The inputs are
// synthetic alnreg arrays covering everything the candidate filter looks at:
// unmapped mates, mates on different contigs, ambiguous mates (cal_sub above
// MIN_RATIO of the score), all four orientations, inserts past max_ins, an odd
// trailing read, passes smaller than the thread count (idle workers), and
// repeated calls that change the thread count or max_ins between them.

#include "doctest/doctest.h"
#include "bwamem.h"

#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <random>
#include <vector>

namespace {

const int64_t L_PAC = 50000000;

struct OptGuard {
    mem_opt_t *o;
    OptGuard() : o(mem_opt_init()) {}
    ~OptGuard() { free(o); }
};

// n reads' worth of alnreg vectors; owns the region storage.
struct Regs {
    std::vector<mem_alnreg_v> v;
    std::vector<std::vector<mem_alnreg_t>> store;

    explicit Regs(int n) : v((size_t)n), store((size_t)n) {
        for (auto &r : v) memset(&r, 0, sizeof(r));
    }
    // Give read i one primary region (plus an optional overlapping secondary
    // with score `sub_score`, which makes cal_sub return that score).
    void set(int i, int rid, int64_t rb, int score, int sub_score = -1) {
        mem_alnreg_t a;
        memset(&a, 0, sizeof(a));
        a.rid = rid; a.rb = rb; a.re = rb + 150; a.qb = 0; a.qe = 150; a.score = score;
        store[i].push_back(a);
        if (sub_score >= 0) { a.score = sub_score; a.rb += 7; store[i].push_back(a); }
        v[i].n = store[i].size(); v[i].m = v[i].n; v[i].a = store[i].data();
    }
};

// A pair at forward-strand position `pos` with `ins` between the mates, in
// orientation d (0 FF, 1 FR, 2 RF, 3 RR), using mem_infer_dir's strand coding
// (rb >= l_pac is the reverse strand, mirrored about 2 * l_pac - 1).
void set_pair(Regs &R, int pair, int64_t pos, int64_t ins, int d, int rid = 0)
{
    const int64_t mirror = 2 * L_PAC - 1;
    int64_t b1 = pos, b2 = pos + ins;
    switch (d) {
    case 0: break;                                  // F, F downstream
    case 1: b2 = mirror - (pos + ins); break;       // F, R downstream
    case 2: b1 = pos + ins; b2 = mirror - pos; break;  // R upstream of F
    case 3: b1 = pos + ins; b2 = pos; break;        // F, F upstream
    }
    R.set(2 * pair, rid, b1, 100);
    R.set(2 * pair + 1, rid, b2, 100);
}

// Random cohort exercising every branch of the candidate filter.
Regs random_cohort(int n, unsigned seed)
{
    std::mt19937_64 rng(seed);
    Regs R(n);
    for (int p = 0; p < n / 2; ++p) {
        const int kind = (int)(rng() % 20);
        const int64_t pos = 1000 + (int64_t)(rng() % 1000000);
        if (kind == 0) {                                   // one mate unmapped
            R.set(2 * p, 0, pos, 100);
        } else if (kind == 1) {                            // mates on different contigs
            R.set(2 * p, 0, pos, 100);
            R.set(2 * p + 1, 1, pos + 300, 100);
        } else if (kind == 2) {                            // ambiguous mate
            R.set(2 * p, 0, pos, 100, /*sub_score=*/95);
            R.set(2 * p + 1, 0, pos + 300, 100);
        } else if (kind == 3) {                            // insert past max_ins
            set_pair(R, p, pos, 10000 + (int64_t)(rng() % 5000), 1);
        } else {
            // Mostly FR, some of every other orientation, so all four survive
            // MIN_DIR_CNT / MIN_DIR_RATIO on larger cohorts.
            const int d = kind < 14 ? 1 : (kind - 14) % 4;
            const int64_t ins = 150 + (int64_t)(rng() % 500) + (int64_t)(rng() % 3) * 40;
            set_pair(R, p, pos, ins, d);
        }
    }
    if (n & 1) R.set(n - 1, 0, 12345, 100);               // trailing unpaired read
    return R;
}

void check_same(const mem_pestat_t a[4], const mem_pestat_t b[4])
{
    for (int d = 0; d < 4; ++d) {
        CAPTURE(d);
        CHECK(a[d].failed == b[d].failed);
        CHECK(a[d].low == b[d].low);
        CHECK(a[d].high == b[d].high);
        // Bitwise, not approximately: the summation order must be identical.
        CHECK(memcmp(&a[d].avg, &b[d].avg, sizeof(double)) == 0);
        CHECK(memcmp(&a[d].std, &b[d].std, sizeof(double)) == 0);
    }
}

void check_par_equals_serial(const mem_opt_t *opt, Regs &R, int nthreads)
{
    const int n = (int)R.v.size();
    CAPTURE(n); CAPTURE(nthreads);
    mem_pestat_t serial[4], par[4];
    memset(par, 0xA5, sizeof(par));   // nothing may survive from before the call
    mem_pestat(opt, L_PAC, n, R.v.data(), serial);
    mem_pestat_par(opt, L_PAC, n, R.v.data(), par, nthreads);
    check_same(serial, par);
    // A large default-max_ins cohort must actually fit a distribution, or the
    // comparison above would pass vacuously on all-failed orientations.
    if (n >= 16 * 1024 && opt->max_ins >= 10000) CHECK(serial[1].failed == 0);
}

} // namespace

TEST_CASE("mem_pestat_par equals mem_pestat across thread counts and cohort sizes"
          * doctest::test_suite("unit/pair")) {
    OptGuard g;
    const int threads[] = { 1, 2, 3, 16 };
    const int ns[] = { 0, 1, 2, 3, 41, 2 * 512 + 1, 16 * 1024, 60001 };
    unsigned seed = 1;
    for (int n : ns) {
        Regs R = random_cohort(n, seed++);
        for (int nt : threads) check_par_equals_serial(g.o, R, nt);
    }
}

TEST_CASE("mem_pestat_par: repeated calls carry nothing between them"
          * doctest::test_suite("unit/pair")) {
    OptGuard g;
    // A big pass that every worker counts into, then a tiny pass that only one
    // or two workers touch (the rest must be skipped, not summed). Then more
    // threads (the pool is rebuilt), a larger max_ins, then back down.
    Regs big = random_cohort(40000, 11), tiny = random_cohort(300, 12);
    check_par_equals_serial(g.o, big, 4);
    check_par_equals_serial(g.o, tiny, 4);
    check_par_equals_serial(g.o, big, 16);
    check_par_equals_serial(g.o, tiny, 16);
    g.o->max_ins = 14000;   // now the past-10000 inserts count
    check_par_equals_serial(g.o, big, 16);
    g.o->max_ins = 500;     // and most inserts drop out
    check_par_equals_serial(g.o, big, 8);
    check_par_equals_serial(g.o, tiny, 3);
    kt_pool_destroy();
}

TEST_CASE("mem_pestat_par: a known distribution, independent of mem_pestat"
          * doctest::test_suite("unit/pair")) {
    OptGuard g;
    // 600 FR pairs: 200 each at inserts 300, 310 and 320. Every value is inside
    // the outlier bounds (p25 300, p75 320), so avg is exactly 310 and std is
    // sqrt((100 + 0 + 100) / 3).
    // The proper-pair bounds are p25 - 3 * IQR = 240 and p75 + 3 * IQR = 380
    // (both outside avg -/+ 4 std, so neither is widened). avg is a sum of
    // integers divided by 600, so exactly 310; std goes through sqrt.
    Regs R(1200);
    for (int p = 0; p < 600; ++p) set_pair(R, p, 5000 + 97 * p, 300 + 10 * (p % 3), 1);
    mem_pestat_t pes[4];
    mem_pestat_par(g.o, L_PAC, 1200, R.v.data(), pes, 16);
    CHECK(pes[0].failed == 1);
    CHECK(pes[1].failed == 0);
    CHECK(pes[2].failed == 1);
    CHECK(pes[3].failed == 1);
    CHECK(pes[1].low == 240);
    CHECK(pes[1].high == 380);
    CHECK(pes[1].avg == 310.0);
    CHECK(pes[1].std == doctest::Approx(8.16496580927726));
    kt_pool_destroy();
}

TEST_CASE("mem_pestat_par: an out-of-range max_ins is clamped the same on both paths"
          * doctest::test_suite("unit/pair")) {
    // max_ins is clamped to [0, 1000000] before it sizes any histogram.
    OptGuard g;
    Regs R = random_cohort(5000, 21);
    // Clamped to the cap: 16 threads' histograms (32 MB each) are past
    // mem_pestat_par's memory bound, so this is the serial fallback.
    g.o->max_ins = 2000000000;
    check_par_equals_serial(g.o, R, 16);
    // Clamped to the cap on 2 threads (64 MB in all, inside the bound): the
    // parallel path with the full 1000001-entry histograms, so the candidate
    // filter must use the clamped bound there too.
    g.o->max_ins = 1500000;
    check_par_equals_serial(g.o, R, 2);
    // Negative clamps to 0: the parallel path with one-entry histograms, where
    // no insert passes the filter and every orientation fails.
    g.o->max_ins = -5;
    check_par_equals_serial(g.o, R, 16);
    kt_pool_destroy();
}
