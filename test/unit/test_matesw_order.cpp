// test/unit/test_matesw_order.cpp — the mate-rescue batch order, computed on (len1, idx) keys.
//
// mem_sam_pe_batch scores a rescue batch (or, under --meth, each OT/OB group of it) with the 8-bit
// pairs first, each width length-sorted with std::sort by len1, and the 8-bit banded parents moved
// to the front with a stable partition. It used to move the 72-byte SeqPairs to get there
// (sort_classify's split, or the --meth group copy, then the sort and the partition in place); it
// now orders an 8-byte (len1, idx) key per pair and gathers the SeqPairs once (matesw_order). The
// order of equal lengths decides which pairs share a kswv lane group, so the gathered batch must be
// the same array, every field of every pair.
//
// The oracle below replicates the record-moving code it replaced: the stable 8-bit-first split
// (sort_classify, or the --meth group copy), std::sort over the SeqPairs with the same comparator,
// then RescueBandBatch::partition's stable partition. It runs the same std::sort as the code under
// test, so it pins equivalence with the old code on whichever standard library builds the test.
// Every case compares the two arrays field by field: batch sizes below and well above std::sort's
// insertion-sort cutoff, few distinct lengths (long runs of ties), all-equal lengths, all-8-bit,
// all-16-bit and interleaved batches, every pair (the non-meth path) and each --meth hypothesis,
// no banded table, no banded parents, every 8-bit pair banded, and a regid outside the table.
//
// Fixtures are built in-memory; no test data files are read.
#include "doctest/doctest.h"
#include "bandedSWA.h"
#include "ksw.h"
#include <algorithm>
#include <climits>
#include <cstdint>
#include <random>
#include <vector>

int64_t bwamem3_matesw_order_by_index(const SeqPair *sp, int64_t n, bool any_hyp, int hyp, bool by_len,
                                      const uint8_t *banded, int n_regid, SeqPair *dst, int64_t *n8_out,
                                      int64_t *nb_out);

namespace {

// The test's marker for "every pair" (the non-meth batch); order() passes it to the hook as any_hyp.
const int ANY_HYP = INT_MIN;

// The pairs of hyp, 8-bit first, then length-sorted and partitioned, by moving the records.
std::vector<SeqPair> oracle(const std::vector<SeqPair> &all, int hyp, bool by_len,
                            const std::vector<uint8_t> &banded, int64_t *n8_out, int64_t *nb_out)
{
    auto is_banded = [&](int regid) {
        return regid >= 0 && regid < (int)banded.size() && banded[regid] != 0;
    };
    std::vector<SeqPair> sp;
    for (const SeqPair &p : all)
        if ((hyp == ANY_HYP || p.meth_hyp == hyp) && (p.h0 & KSW_XBYTE)) sp.push_back(p);
    const int64_t n8 = (int64_t)sp.size();
    for (const SeqPair &p : all)
        if ((hyp == ANY_HYP || p.meth_hyp == hyp) && !(p.h0 & KSW_XBYTE)) sp.push_back(p);
    if (by_len) {
        auto by_len1 = [](const SeqPair &x, const SeqPair &y) { return x.len1 < y.len1; };
        std::sort(sp.begin(), sp.begin() + n8, by_len1);
        std::sort(sp.begin() + n8, sp.end(), by_len1);
    }
    std::vector<SeqPair> out(sp);
    int nb = 0;
    if (!banded.empty()) {
        for (int64_t i = 0; i < n8; i++) nb += is_banded(sp[i].regid);
        int a = 0, o = nb;
        for (int64_t i = 0; i < n8; i++) out[is_banded(sp[i].regid) ? a++ : o++] = sp[i];
    }
    *n8_out = n8;
    *nb_out = nb;
    return out;
}

// n pairs, each 8-bit with probability pct8 %, len1 drawn from [lo, hi], meth_hyp drawn from {0, 1}
// when meth (else -1), in random order; every field set to a value unique to the pair so a
// misplaced pair cannot compare equal by accident.
std::vector<SeqPair> make_batch(std::mt19937 &rng, int64_t n, int pct8, int lo, int hi, bool meth)
{
    std::uniform_int_distribution<int> len(lo, hi), pct(0, 99), bit(0, 1);
    std::vector<SeqPair> sp(n);
    for (int64_t i = 0; i < n; i++) {
        SeqPair p;
        p.idr = (int32_t)(i * 7);
        p.idq = (int32_t)(i * 11);
        p.id = (int32_t)i;
        p.len1 = len(rng);
        p.len2 = (int32_t)(i % 151);
        p.h0 = (pct(rng) < pct8 ? KSW_XBYTE : 0) | KSW_XSTART | (int32_t)(i & 0xff);
        p.seqid = (int)(i / 2);
        p.regid = (int)i;
        p.score = (int32_t)(1000 + i);
        p.tle = (int32_t)(2000 + i);
        p.gtle = (int32_t)(3000 + i);
        p.qle = (int32_t)(4000 + i);
        p.gscore = (int32_t)(5000 + i);
        p.max_off = (int32_t)(6000 + i);
        p.tight_band = (int32_t)(i % 97);
        p.chain_band = (int32_t)(i % 89);
        p.ugp_walk_score = (int32_t)(7000 + i);
        p.ugp_r_attempted = (uint8_t)(i & 1);
        p.meth_hyp = meth ? (int8_t)bit(rng) : (int8_t)-1;
        sp[i] = p;
    }
    // The batch arrives in enqueue order, not regid order.
    std::shuffle(sp.begin(), sp.end(), rng);
    return sp;
}

// Every field of a and b, in declaration order (memcmp would also compare the tail padding,
// which a struct copy need not carry).
bool same_pair(const SeqPair &a, const SeqPair &b)
{
    return a.idr == b.idr && a.idq == b.idq && a.id == b.id && a.len1 == b.len1 && a.len2 == b.len2 &&
           a.h0 == b.h0 && a.seqid == b.seqid && a.regid == b.regid && a.score == b.score &&
           a.tle == b.tle && a.gtle == b.gtle && a.qle == b.qle && a.gscore == b.gscore &&
           a.max_off == b.max_off && a.tight_band == b.tight_band && a.chain_band == b.chain_band &&
           a.ugp_walk_score == b.ugp_walk_score && a.ugp_r_attempted == b.ugp_r_attempted &&
           a.meth_hyp == b.meth_hyp;
}

bool same_batch(const std::vector<SeqPair> &a, const std::vector<SeqPair> &b)
{
    if (a.size() != b.size()) return false;
    for (size_t i = 0; i < a.size(); i++)
        if (!same_pair(a[i], b[i])) return false;
    return true;
}

// The hook's gathered batch for hyp (n8 in *n8_out, the banded-parent count in *nb_out).
std::vector<SeqPair> order(const std::vector<SeqPair> &sp, int hyp, bool by_len,
                           const std::vector<uint8_t> &banded, int64_t *n8_out, int64_t *nb_out)
{
    std::vector<SeqPair> got(sp.size());
    const int64_t m = bwamem3_matesw_order_by_index(sp.data(), (int64_t)sp.size(), hyp == ANY_HYP,
                                                    hyp == ANY_HYP ? 0 : hyp, by_len, banded.data(),
                                                    (int)banded.size(), got.data(), n8_out, nb_out);
    got.resize((size_t)m);
    return got;
}

void check_same(const std::vector<SeqPair> &sp, int hyp, bool by_len, const std::vector<uint8_t> &banded)
{
    CAPTURE(hyp);
    CAPTURE(by_len);
    int64_t want8 = -1, got8 = -2, wantb = -1, gotb = -2;
    const std::vector<SeqPair> want = oracle(sp, hyp, by_len, banded, &want8, &wantb);
    const std::vector<SeqPair> got = order(sp, hyp, by_len, banded, &got8, &gotb);
    CHECK(want8 == got8);
    CHECK(wantb == gotb);
    CHECK(same_batch(want, got));
}

std::vector<uint8_t> banded_table(std::mt19937 &rng, int n_regid, int pct)
{
    std::uniform_int_distribution<int> d(0, 99);
    std::vector<uint8_t> b(n_regid);
    for (auto &x : b) x = d(rng) < pct ? 1 : 0;
    return b;
}

} // namespace

TEST_CASE("matesw order by index: random batches match the record split, sort and partition"
          * doctest::test_suite("unit/pair"))
{
    std::mt19937 rng(20261020);
    const int64_t sizes[] = {0, 1, 2, 15, 16, 17, 33, 100, 1100};
    const int len_hi[] = {0, 3, 40, 250};
    const int banded_pct[] = {-1, 0, 10, 50, 100};   // -1: no banded table (no partition)
    const int pct8s[] = {100, 0, 75};
    int interesting = 0;   // cases with ties, some but not all 8-bit pairs banded, and both widths
    for (int64_t n : sizes)
        for (int hi : len_hi)
            for (int bp : banded_pct)
                for (int pct8 : pct8s)
                    for (int meth = 0; meth < 2; meth++) {
                        CAPTURE(n);
                        CAPTURE(hi);
                        CAPTURE(bp);
                        CAPTURE(pct8);
                        CAPTURE(meth);
                        const std::vector<SeqPair> sp = make_batch(rng, n, pct8, 1, hi + 1, meth != 0);
                        const std::vector<uint8_t> banded =
                            bp < 0 ? std::vector<uint8_t>() : banded_table(rng, (int)n, bp);
                        for (int hyp : meth ? std::vector<int>{1, 0} : std::vector<int>{ANY_HYP}) {
                            check_same(sp, hyp, true, banded);
                            check_same(sp, hyp, false, banded);
                        }
                        int n8 = 0, nb = 0;
                        for (const SeqPair &p : sp)
                            if (p.h0 & KSW_XBYTE) { n8++; nb += !banded.empty() && banded[p.regid]; }
                        interesting += hi <= 40 && n > 100 && nb > 0 && nb < n8 && n8 < n;
                    }
    CHECK(interesting > 10);
}

TEST_CASE("matesw order by index: ties keep the record sort's order"
          * doctest::test_suite("unit/pair"))
{
    // Many pairs, two lengths: the order within each run of equal len1 is introsort's, and it
    // decides the kswv lane groups. A stable sort produces a different array here (std::sort is
    // introsort in both libstdc++ and libc++, which is not stable on this input), so this case
    // checks the order of ties, not just the sortedness.
    std::mt19937 rng(7);
    const std::vector<SeqPair> sp = make_batch(rng, 2000, 100, 100, 101, false);
    const std::vector<uint8_t> none;
    check_same(sp, ANY_HYP, true, none);

    std::vector<SeqPair> stable(sp);
    std::stable_sort(stable.begin(), stable.end(),
                     [](const SeqPair &x, const SeqPair &y) { return x.len1 < y.len1; });
    int64_t n8 = 0, nbx = 0;
    CHECK(!same_batch(stable, order(sp, ANY_HYP, true, none, &n8, &nbx)));
}

TEST_CASE("matesw order by index: banded parents come first, in sorted order"
          * doctest::test_suite("unit/pair"))
{
    std::mt19937 rng(11);
    const std::vector<SeqPair> sp = make_batch(rng, 300, 67, 50, 150, false);
    std::vector<uint8_t> banded(300, 0);
    for (int r = 0; r < 300; r += 3) banded[r] = 1;     // includes 16-bit regids, which stay put
    check_same(sp, ANY_HYP, true, banded);

    int64_t n8 = 0, got_nb = -1;
    const std::vector<SeqPair> got = order(sp, ANY_HYP, true, banded, &n8, &got_nb);
    REQUIRE(got.size() == 300);
    REQUIRE(n8 > 0);
    REQUIRE(n8 < 300);
    int nb = 0;
    for (const SeqPair &p : sp) nb += (p.h0 & KSW_XBYTE) && banded[p.regid];
    REQUIRE(nb > 0);
    CHECK(got_nb == nb);
    for (int i = 0; i < n8; i++) CHECK((got[i].h0 & KSW_XBYTE) != 0);
    for (int i = (int)n8; i < 300; i++) CHECK((got[i].h0 & KSW_XBYTE) == 0);
    for (int i = 0; i < n8; i++) CHECK((banded[got[i].regid] != 0) == (i < nb));
    for (int i = 1; i < nb; i++) CHECK(got[i - 1].len1 <= got[i].len1);
    for (int i = nb + 1; i < n8; i++) CHECK(got[i - 1].len1 <= got[i].len1);
    for (int i = (int)n8 + 1; i < 300; i++) CHECK(got[i - 1].len1 <= got[i].len1);

    // A regid past the banded table (or negative) is not a banded parent.
    std::vector<SeqPair> odd(sp);
    odd[0].regid = 5000;
    odd[1].regid = -1;
    check_same(odd, ANY_HYP, true, banded);
}
