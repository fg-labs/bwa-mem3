// test/unit/test_bandedswa_zdrop_wrap.cpp
//
// The 8-bit banded extension kernels test z-drop per row as drop - dif > zdrop,
// with dif = |drift| * e (e the gap extend on the side the drift points to)
// formed in int32. With a gap extend near INT32_MAX, which -E accepts, that
// product wraps, and the step can kill a lane at a drop that does not exceed
// zdrop. Several kernels skip the step on rows where no lane has drop > zdrop:
// the lean NEON row (BSW8_ROW_LEAN), the Apple-silicon epilogue gate, and the
// AVX2 / AVX-512BW epilogue skips and their lane-compaction drivers. That is a
// necessary condition only while dif >= 0, so each gate steps aside for such
// gap extends (BSW8_ZDROP_GATE_EXACT in bandedSWA.cpp) and the step runs on
// every row, as in the ungated kernels.
//
// The pairs below are ones a gate without that guard changes. Their expected
// results are those of the ungated 8-bit kernel (which wraps the same way on
// every target), the reference the gates must reproduce; they were taken from
// the ungated original kernel (BSW8_ROW_LEAN=0 on Linux arm64), which the lean
// row matched on every pair. scalarBandedSWA is not an oracle here: the 8-bit
// kernels broadcast the gap penalties as bytes (set1_epi8), so at gap extends
// this large they and the scalar code can disagree, gate or no gate.
//
// The case runs on whatever getScores8 the build has, so CI checks the lean
// NEON row on Linux arm64, the Apple gate on macOS, and the AVX2 (and, where
// the runner has it, AVX-512BW) skip on x86. Each pair is scored alone, as a
// batch of copies filling every lane of a group, and through the
// lane-compaction driver where the tier has one.

#include <cstdint>
#include <memory>
#include <string>
#include <vector>

#include "doctest/doctest.h"
#include "bandedSWA.h"
#include "bsw_batch.h"
#include "ext_ladder.h"

#if HAVE_BSW_VECTOR_8_16

namespace {

using bwa_tests::ExtPair;
using bwa_tests::ExtResult;
using bwa_tests::ExtScoring;

/// A pair whose z-drop drift wraps int32, with the ungated kernel's result.
struct WrapCase {
    int a, b, o_del, o_ins, e_del, e_ins, pen_clip, zdrop, w, h0;
    const char *ref, *qry;   // bases as digits 0-4
    ExtResult want;          // score, qle, tle, gtle, gscore, max_off
};

// Found by scoring random pairs at huge gap extends with the Apple-silicon gate
// compiled in without the guard, against the ungated kernel; the gate kept each
// of these alive past the row where the wrapped step kills it. Only the
// insertion extend is huge in the fifth pair and only the deletion one in the
// sixth. The positive gap extends are values -E accepts (1..INT_MAX), with
// o + e still inside int (the kernels form o + e in int, so a larger sum would
// make the test itself depend on signed overflow); the last pair's negative
// extend, which makes dif negative without any wrap, is reachable only through
// the kernel API.
const WrapCase kWrapCases[] = {
    {5, 4, 2, 2, 1073741825, 536870912, 6, 58, 105, 40, "3212031220213032022322", "321023122223",
     {58, 5, 4, 0, -1, 1}},
    {2, 2, 1, 1, 268435456, 1073741825, 8, 24, 43, 36, "1032232022133113202222132231", "112222330",
     {39, 4, 5, 0, -1, 1}},
    {4, 4, 0, 0, 1073741825, 1073741824, 5, 3, 46, 207, "1023012111031320302223220200310302312",
     "3003112", {210, 5, 4, 0, -1, 1}},
    {3, 3, 0, 1, 1073741825, 1073741823, 1, 57, 20, 62, "20113023303110103110230321303",
     "20113023301120103110230", {92, 11, 10, 0, -1, 1}},
    {5, 1, 1, 4, 2, 1073741825, 2, 1, 27, 76, "023331223213023220213012202332011133322131312020",
     "000313233322311311123202211101", {95, 8, 9, 0, -1, 1}},
    {5, 1, 6, 4, 1073741825, 5, 0, 7, 47, 83, "10202132233223110202010231300120213122230100311",
     "033322200233210100221332313230330", {104, 19, 18, 0, -1, 1}},
    {3, 2, 3, 6, -3, 2, 5, 3, 75, 113, "20313201230031211021333030", "202230032300322110",
     {119, 2, 2, 0, -1, 0}},
};

std::string show(const ExtResult &x)
{
    return std::to_string(x.score) + "/" + std::to_string(x.qle) + "/" + std::to_string(x.tle) +
           "/" + std::to_string(x.gtle) + "/" + std::to_string(x.gscore) + "/" +
           std::to_string(x.max_off);
}

/// getScores8 on `n` copies of `p`, with compaction off or forced (two groups
/// per superblock, any band); every copy must return `want`. Returns the number
/// of compaction superblocks the kernel ran.
uint64_t check_copies(const ExtScoring &sc, const WrapCase &c, const ExtPair &p, int n, int groups)
{
    std::unique_ptr<IBandedPairWiseSW> k(new BandedPairWiseSW(
        sc.o_del, sc.e_del, sc.o_ins, sc.e_ins, sc.zdrop, sc.pen_clip, sc.mat,
        (int8_t)sc.a, (int8_t)sc.b, 1));
    k->set_lane_compaction(groups, 0);
    const std::vector<const ExtPair *> batch(n, &p);
    const std::vector<ExtResult> got = bwa_tests::score_bsw_batch(*k, batch, 8, c.w);
    for (int q = 0; q < n; q++) {
        INFO("copies=" << n << " groups=" << groups << " copy " << q << ": got " << show(got[q])
             << " want " << show(c.want) << " (score/qle/tle/gtle/gscore/max_off)");
        CHECK(got[q] == c.want);
        if (!(got[q] == c.want)) break;
    }
    return k->lane_compaction_counts().superblocks;
}

} // namespace

TEST_CASE("8-bit extension z-drop gates: exact for gap extends whose drift wraps int32"
          * doctest::test_suite("unit/bandedswa")) {
    uint64_t superblocks = 0;
    for (const WrapCase &c : kWrapCases) {
        ExtScoring sc;
        sc.a = c.a; sc.b = c.b; sc.o_del = c.o_del; sc.o_ins = c.o_ins;
        sc.e_del = c.e_del; sc.e_ins = c.e_ins; sc.pen_clip = c.pen_clip; sc.zdrop = c.zdrop;
        sc.fill_mat();
        ExtPair p;
        p.h0 = c.h0;
        for (const char *q = c.ref; *q; q++) p.target.push_back((uint8_t)(*q - '0'));
        for (const char *q = c.qry; *q; q++) p.query.push_back((uint8_t)(*q - '0'));
        CAPTURE(c.e_del);
        CAPTURE(c.e_ins);
        CAPTURE(c.zdrop);
        CAPTURE(c.w);
        check_copies(sc, c, p, 1, 0);     // alone: the other lanes are padding
        check_copies(sc, c, p, 64, 0);    // every lane of a group (16, 32 or 64 wide)
        superblocks += check_copies(sc, c, p, 256, 2);   // the compaction driver, where there is one
    }
#if defined(__AVX512BW__)
    const std::string tier = "avx512bw";
#elif defined(__AVX2__)
    const std::string tier = "avx2";
#elif defined(__aarch64__)
    const std::string tier = "neon";
#else
    const std::string tier = "128-bit x86";
#endif
    MESSAGE("kernel tier " << tier << "; compaction superblocks run: " << superblocks
            << " (0 on tiers without a driver)");
#if defined(__aarch64__) || defined(__AVX2__)
    CHECK(superblocks > 0);   // this tier has a compaction driver: it must have run
#endif
}

#else

TEST_CASE("8-bit extension z-drop gates: exact for gap extends whose drift wraps int32"
          * doctest::test_suite("unit/bandedswa")) {
    MESSAGE("SKIP: this build has no 8-bit vector extension kernel");
}

#endif
