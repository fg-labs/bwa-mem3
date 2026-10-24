// test/unit/test_rescue_kmer_anchor.cpp — the --rescue-kmer anchor scan's k-mer table.
//
// matesw_kmer_anchor (src/bwamem_pair.cpp) indexes the oriented mate's k-mers in a per-thread
// open-addressing table, votes the window's k-mers onto diagonals, and resets only the slots it
// claimed. The table's empty marker must be a value no k-mer code can take. At k = 16 every
// uint32 is a code (the all-T k-mer is 0xFFFFFFFF), so a 32-bit marker made an all-T k-mer
// indistinguishable from an empty slot: the query's all-T k-mer was counted once per occurrence
// and could be overwritten, and the window's all-T k-mers voted with whatever position an
// earlier read on the same thread had left in the first empty slot on their probe.
//
// Checked through mem_matesw_kmer_anchor_counts, at k = 15 and 16:
//   1. an all-T mate indexes exactly one distinct k-mer;
//   2. a mate with no T, scanned against a window that is the mate followed by a long T run,
//      gets exactly one vote per query k-mer on the true diagonal, however many other reads
//      the thread scanned first.
//
// Sequences are generated in-memory; no test data files are read.
#include "doctest/doctest.h"
#include "bwamem.h"
#include <cstdint>
#include <cstdlib>
#include <random>
#include <set>
#include <vector>

namespace {

struct Counts {
    int best_v, n_kmer;
};

Counts scan(const mem_opt_t *opt, const std::vector<uint8_t> &ref, const std::vector<uint8_t> &ms)
{
    Counts c;
    mem_matesw_kmer_anchor_counts(opt, ref.data(), (int)ref.size(), ms.data(), (int)ms.size(), 0, 0,
                                  &c.best_v, &c.n_kmer);
    return c;
}

std::vector<uint8_t> random_bases(std::mt19937 &rng, int n, int alphabet)
{
    std::uniform_int_distribution<int> base(0, alphabet - 1);
    std::vector<uint8_t> v((size_t)n);
    for (auto &b : v) b = (uint8_t)base(rng);
    return v;
}

}  // namespace

TEST_CASE("rescue k-mer anchor: an all-T mate indexes one distinct k-mer"
          * doctest::test_suite("unit/pair")) {
    mem_opt_t *opt = mem_opt_init();
    std::mt19937 rng(16);
    for (int k : {15, 16}) {
        CAPTURE(k);
        opt->rescue_kmer = k;
        const std::vector<uint8_t> ms(40, 3);   // 25 (k = 16) or 26 (k = 15) copies of one k-mer
        const std::vector<uint8_t> ref = random_bases(rng, 200, 4);
        const Counts c = scan(opt, ref, ms);
        CHECK(c.n_kmer == 1);
    }
    free(opt);
}

TEST_CASE("rescue k-mer anchor: a window's T run casts no votes left over from earlier reads"
          * doctest::test_suite("unit/pair")) {
    mem_opt_t *opt = mem_opt_init();
    std::mt19937 rng(2026);
    for (int k : {15, 16}) {
        CAPTURE(k);
        opt->rescue_kmer = k;

        // A, C and G only, so the mate has no all-T k-mer; its k-mers are all distinct.
        const int l_ms = 120;
        std::vector<uint8_t> ms;
        for (;;) {
            ms = random_bases(rng, l_ms, 3);
            std::set<std::vector<uint8_t>> seen;
            for (int i = 0; i + k <= l_ms; i++) seen.emplace(ms.begin() + i, ms.begin() + i + k);
            if ((int)seen.size() == l_ms - k + 1) break;
        }
        // The window is the mate (diagonal 0) followed by 300 T: every T-run k-mer is the
        // all-T code, which the mate does not hold, so it must not vote.
        std::vector<uint8_t> ref = ms;
        ref.insert(ref.end(), 300, 3);
        const int want = l_ms - k + 1;

        const Counts fresh = scan(opt, ref, ms);
        CHECK(fresh.best_v == want);
        CHECK(fresh.n_kmer == want);

        // Earlier reads on this thread leave positions in slots the reset marks empty; a long
        // random read fills about a third of the table each time, at positions past l_ms, so
        // the T run's probes would meet one if they could match an empty slot.
        for (int r = 0; r < 64; r++) {
            CAPTURE(r);
            const std::vector<uint8_t> prior = random_bases(rng, 400, 4);
            const std::vector<uint8_t> prior_ref = random_bases(rng, 600, 4);
            scan(opt, prior_ref, prior);
            const Counts c = scan(opt, ref, ms);
            CHECK(c.best_v == want);
            CHECK(c.n_kmer == want);
        }
    }
    free(opt);
}
