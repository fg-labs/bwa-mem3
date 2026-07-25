// test/unit/test_kswv_mixed_qlen.cpp
//
// kswv::getScores8 vs scalar ksw_align2 on MIXED-QUERY-LENGTH batches.
//
// kswvBatchWrapper8 pads each lane's query with DUMMY5 up to that lane's own
// 16-base quantum, then fills every column from there to the group's widest
// quantum with 0xFF -- the same high-bit sentinel used for reference rows past
// len1. The kernel must therefore zero the DP diagonal term wherever bit 7 is
// set in the reference base OR the query base.
//
// The NEON and AVX-512BW 8-bit kernels had hoisted that test to the reference
// base alone, on the incorrect premise that a query byte never carries bit 7.
// On a group whose queries straddle a quantum boundary the short lanes' 0xFF
// columns were then left unmasked, carried a diagonal forward, and inflated
// that lane's per-row maximum -- surfacing as a phantom second-best score,
// which mem_matesw_batch_post stores as mem_alnreg_t::csub and the SAM writer
// emits as XS.
//
// test_kswv_correctness.cpp already batches mixed query lengths (50-128), so
// layout alone does not reproduce this: it takes content too. gen_random_pair
// makes ref and query independent, so the DP never rises far enough for a
// spurious row to clear minsc. The tandem-repeat generator used here produces
// a high primary score plus many near-threshold secondary maxima, which is what
// turns an unmasked padding column into an observable score2.
//
// A uniform-length batch cannot expose this IN PHASE 0 -- every lane's quantum
// equals the group's, so no 0xFF query column exists. The lengths below
// alternate per lane so that ANY group width (16, 32 or 64) sees both. Note
// this does NOT make a fixed-length library safe: production phase 1 sets
// len2 = qe + 1 per pair (bwamem_pair.cpp:875), so phase-1 groups have ragged
// quanta whatever the input lengths were. Both phases are checked below.
//
// The --meth cases at the bottom cover the HasFreed=true instantiation, which
// pre-fix was hit far harder than the symmetric one and only on NEON: that
// kernel applied the hoisted rowboundary to BOTH instantiations, where
// AVX-512BW gated it and so kept meth correct, and active_frread's old 0xFF
// marker additionally aliased the 0xFF pad so the freed blend fired on every
// padded column. Symmetric pre-fix loses only score2/te2; meth pre-fix loses
// score, qe/te and phase 1's tb/qb as well.
//
// The pair count is deliberately odd, not a multiple of 64: being odd makes it
// coprime to every 8-bit group width (16/32/64) and every 16-bit one (8/16/32),
// so the final group is always partial. Only a partial group reaches the
// kswvBatchWrapper* prologue that synthesizes len1=len2=0 dummy lanes out to
// roundNumPairs -- exactly the writes BatchBuffers sizes its MAX_SIMD_WIDTH8
// tail padding to absorb. A dummy lane also has query_quantum8(0) == 0, which
// pulls the group's jsplit to 0 and forces the per-cell boundary mask across
// every column, rather than the hoisted per-row form a uniform group takes.
//
// getScores16 carries the identical contract (0xFFFF past query_quantum16) and
// is covered too, so the invariant is locked everywhere it can drift rather
// than only in the two kernels that happened to be wrong. It needs match=14:
// that puts min_seed_len*match = 266 >= 250, so default_xtra_flags drops
// KSW_XBYTE and BOTH sides run word mode. Comparing the 16-bit kernel against
// a byte-mode scalar oracle instead would report spurious score2 differences --
// scalar ksw_u8 is striped and takes its row max before the lazy-F fix-up, so
// the two precisions legitimately collect different b[] entries.

#include <random>
#include <vector>

#include "doctest/doctest.h"

#include "ksw_runner.h"
#include "kswr_cmp.h"
#include "kswv_runner.h"  // BWA_TESTS_HAVE_KSWV
#include "scoring.h"
#include "seqpair_batch.h"
#include "seqpair_gen.h"

#if BWA_TESTS_HAVE_KSWV

#include "bwa.h"            // bwa_fill_scmat
#include "bwamem.h"         // mem_opt_init, mem_opt_fill_meth_mat, MEM_METH_SCORING_*
#include "kswv.h"
#include "simd_dispatch.h"  // make_kswv

namespace {

// Typical mate-rescue reference window (mem_matesw fetches pes.high-ish bp).
constexpr int kRefLen = 460;

// Whether the RUNNING tier has a batched kswv kernel at all. Only AVX2,
// AVX-512BW and NEON do; the sse41/sse42/avx classes are exit() stubs
// (kswv.cpp:3749). These tests reach the kernel through make_kswv, which
// honors BWAMEM3_FORCE_TIER, so a tier sweep down to SSE would otherwise
// take the whole test binary down with it mid-suite rather than skip.
bool batched_kswv_available() {
    bwamem3_simd_init();
    const int tier = bwamem3_simd_tier();
    return tier == BWAMEM3_TIER_AVX2 || tier == BWAMEM3_TIER_AVX512BW ||
           tier == BWAMEM3_TIER_NEON;
}

// Score one batch on the batched kernel and require every forward-pass kswr_t
// field to match the scalar ksw_align2 oracle at the same matrix.
//
// `kernel_mat` is what make_kswv installs: nullptr selects the 9-arg overload
// (symmetric, HasFreed=false); a bisulfite mat25 selects the mat-aware overload
// and, with it, the HasFreed=true kernel instantiation. `mat` is the matrix the
// scalar oracle uses and MUST be the same one, or the comparison is meaningless.
//
// Deliberately uses make_kswv() rather than the framework's run_kswv_batch():
// run_kswv_batch constructs the concrete `kswv` class directly, which binds to
// whichever tier this translation unit was compiled at (BASELINE_ARCH, i.e.
// avx2 on x86). make_kswv dispatches on the RUNTIME tier and honors
// BWAMEM3_FORCE_TIER, so this test reaches the AVX-512BW kernel -- one of the
// two that regressed -- on a capable host, and can be swept across tiers.
// Both phases run: phase 0 for score/score2/te/qe/te2, then phase 1 for tb/qb.
// test_kswv_correctness.cpp also drives both, but only through run_kswv_batch's
// compile-time-tier, symmetric-matrix, independent-ref-and-query configuration,
// which is why neither phase's exposure to the padding contract showed up there.
void compare_batch(uint64_t seed, const std::vector<bwa_tests::TestPair> &pairs,
                   const bwa_tests::ScoringMatrix &mat, const int8_t *kernel_mat,
                   bool use16, int max_qlen) {
    const int xtra = bwa_tests::default_xtra_flags(static_cast<int>(mat[0]));

    std::vector<kswr_t> scalar_aln;
    scalar_aln.reserve(pairs.size());
    for (const auto &p : pairs) {
        scalar_aln.push_back(bwa_tests::run_scalar_ksw(p, mat));
    }

    bwa_tests::BatchBuffers bb(pairs, xtra);
    auto pwsw = kernel_mat
        ? make_kswv(bwa_tests::DEFAULT_GAP_OPEN, bwa_tests::DEFAULT_GAP_EXTEND,
                    bwa_tests::DEFAULT_GAP_OPEN, bwa_tests::DEFAULT_GAP_EXTEND,
                    mat[0], mat[1], 1, kRefLen, max_qlen, kernel_mat)
        : make_kswv(bwa_tests::DEFAULT_GAP_OPEN, bwa_tests::DEFAULT_GAP_EXTEND,
                    bwa_tests::DEFAULT_GAP_OPEN, bwa_tests::DEFAULT_GAP_EXTEND,
                    mat[0], mat[1], 1, kRefLen, max_qlen);
    REQUIRE(pwsw->needsScalar() == false);
    if (use16) {
        pwsw->getScores16(bb.pairs(), bb.ref_buf(), bb.qer_buf(), bb.aln(), bb.n(), 1, 0);
    } else {
        pwsw->getScores8(bb.pairs(), bb.ref_buf(), bb.qer_buf(), bb.aln(), bb.n(), 1, 0);
    }

    // Per-pair CHECKs only on mismatch: at 129 pairs x 5 fields the passing
    // assertions would otherwise bury the signal (matching the convention in
    // test_kswv_correctness.cpp).
    for (int i = 0; i < bb.n(); i++) {
        const kswr_t &want = scalar_aln[i];
        const kswr_t &got  = bb.aln()[bb.pairs()[i].regid];
        CAPTURE(seed);
        CAPTURE(i);
        CAPTURE(use16);
        CAPTURE(bb.pairs()[i].len2);
        if (!bwa_tests::kswr_score_eq(want, got))  CHECK(want.score  == got.score);
        if (!bwa_tests::kswr_score2_eq(want, got)) CHECK(want.score2 == got.score2);
        if (want.te != got.te)                     CHECK(want.te  == got.te);
        if (want.qe != got.qe)                     CHECK(want.qe  == got.qe);
        if (want.te2 != got.te2)                   CHECK(want.te2 == got.te2);
    }

    // Phase 1 (tb/qb recovery) has its own exposure to the padding contract,
    // and a worse-behaved one: production sets len2 = qe + 1 per pair
    // (bwamem_pair.cpp:875), so a phase-1 group's quanta are ragged even when
    // every input read was the same length. A uniform-length library is
    // therefore NOT immune -- phase 0 may be uniform while phase 1 never is.
    //
    // The failure mode differs too. kswv.cpp:887 only writes tb/qb when the
    // phase-1 score reproduces the phase-0 score exactly; a padding-inflated
    // score fails that equality, tb/qb keep the -1 the caller pre-seeded
    // (bwamem_pair.cpp:942), and the `aln.qb >= 0` guard at bwamem_pair.cpp:281
    // then discards the rescue outright. So here the symptom is a SILENTLY
    // DROPPED mate rescue, not a mis-scored one.
    // Guard against a vacuous phase-1 leg: prepare_phase1 keeps only pairs that
    // clear the KSW_XSTART/KSW_XSUBO gate, and a batch that cleared none would
    // run the loop below zero times and still report green.
    const int survivors = bb.prepare_phase1();
    REQUIRE(survivors > 0);
    if (use16) {
        pwsw->getScores16(bb.pairs(), bb.ref_buf(), bb.qer_buf(), bb.aln(), survivors, 1, 1);
    } else {
        pwsw->getScores8(bb.pairs(), bb.ref_buf(), bb.qer_buf(), bb.aln(), survivors, 1, 1);
    }
    for (int i = 0; i < survivors; i++) {
        const int regid = bb.pairs()[i].regid;
        const kswr_t &want = scalar_aln[regid];
        const kswr_t &got  = bb.aln()[regid];
        CAPTURE(seed);
        CAPTURE(regid);
        CAPTURE(use16);
        CAPTURE(bb.pairs()[i].len2);
        if (want.tb != got.tb) CHECK(want.tb == got.tb);
        if (want.qb != got.qb) CHECK(want.qb == got.qb);
    }
}

// Alternating-length tandem-repeat batch: the layout that puts a 0xFF query
// column on the short lanes and the content that makes it observable.
std::vector<bwa_tests::TestPair> mixed_len_batch(std::mt19937 &rng, int n_pairs,
                                                 int len_a, int len_b) {
    std::vector<bwa_tests::TestPair> pairs;
    pairs.reserve(n_pairs);
    for (int i = 0; i < n_pairs; i++) {
        pairs.push_back(bwa_tests::gen_tandem_repeat_pair(
            rng, (i % 2 == 0) ? len_a : len_b, kRefLen));
    }
    return pairs;
}

void check_batch_matches_scalar(uint64_t seed, int n_pairs, int len_a, int len_b,
                                bool use16 = false) {
    // match=1 keeps KSW_XBYTE set (byte mode, both sides); match=14 drops it
    // (word mode, both sides). See the header note.
    const bwa_tests::ScoringMatrix mat =
        use16 ? bwa_tests::build_scoring_matrix(14, 8, 1)
              : bwa_tests::default_scoring_matrix();
    if (!batched_kswv_available()) {
        MESSAGE("tier has no batched kswv kernel; skipping");
        return;
    }
    std::mt19937 rng(static_cast<unsigned>(seed));
    const auto pairs = mixed_len_batch(rng, n_pairs, len_a, len_b);
    compare_batch(seed, pairs, mat, /*kernel_mat=*/nullptr, use16,
                  len_a > len_b ? len_a : len_b);
}

// ---------------------------------------------------------------------------
// Bisulfite (--meth) variants: the HasFreed=true kernel instantiations.
// ---------------------------------------------------------------------------

// Build the production OT/OB matrix for `scoring` at the given match/mismatch
// via the SAME two calls the CLI makes -- bwa_fill_scmat (fastmap.cpp:1745)
// then mem_opt_fill_meth_mat. Restating the freed-cell layout here instead
// would stop testing mem_opt_fill_meth_mat and start testing a copy of it.
bwa_tests::ScoringMatrix meth_matrix(int scoring, bool ot, int match, int mismatch) {
    mem_opt_t *o = mem_opt_init();
    o->a = match;
    o->b = mismatch;
    o->meth_scoring = scoring;
    bwa_fill_scmat(o->a, o->b, o->mat);
    mem_opt_fill_meth_mat(o);
    const int8_t *src = ot ? o->mat_ot : o->mat_ob;
    bwa_tests::ScoringMatrix mat;
    for (int i = 0; i < 25; i++) mat[i] = src[i];
    free(o);
    return mat;
}

// Apply the read-side conversion the matrix frees, so the freed cell actually
// fires: OT frees mat[ref C][read T], OB frees mat[ref G][read A]. Converting
// most (not all) sites mimics a real bisulfite library and leaves both freed
// and unfreed cells in the DP -- converting none would leave HasFreed nominally
// on but never exercised, and converting all would erase the base from the read.
void convert_query_bisulfite(std::mt19937 &rng, bwa_tests::TestPair &p, bool ot) {
    std::uniform_int_distribution<int> pct(0, 99);
    const uint8_t from = ot ? 1 : 2;   // nt4: C=1 (OT), G=2 (OB)
    const uint8_t to   = ot ? 3 : 0;   //      T=3      , A=0
    for (auto &b : p.qry) {
        if (b == from && pct(rng) < 80) b = to;
    }
}

// Returns false (and reports) when the running tier has no freed-cell kernel --
// sse41/sse42/avx report needsScalar() for any freed matrix, and production
// gates the batched meth path off the same property
// (bwamem_pair.cpp: meth_freed_cell_tier_supported).
bool freed_kernel_available(const bwa_tests::ScoringMatrix &mat) {
    auto probe = make_kswv(bwa_tests::DEFAULT_GAP_OPEN, bwa_tests::DEFAULT_GAP_EXTEND,
                           bwa_tests::DEFAULT_GAP_OPEN, bwa_tests::DEFAULT_GAP_EXTEND,
                           mat[0], mat[1], 1, kRefLen, 256, mat.data());
    return !probe->needsScalar();
}

void check_meth_batch_matches_scalar(uint64_t seed, int n_pairs, int len_a, int len_b,
                                     int scoring, bool ot, bool use16 = false) {
    // Same byte-vs-word split as the symmetric case: match=14 puts
    // min_seed_len*match over 250 so default_xtra_flags drops KSW_XBYTE and
    // both sides run word mode.
    const bwa_tests::ScoringMatrix mat =
        use16 ? meth_matrix(scoring, ot, 14, 8) : meth_matrix(scoring, ot, 1, 4);
    // Guard against a vacuous meth leg. If mem_opt_fill_meth_mat ever stopped
    // freeing the conversion cell, the matrix would be symmetric, the mat-aware
    // make_kswv would select HasFreed=false, and every case below would quietly
    // re-run the non-meth path while still reporting green. mat[0*5+1] (ref A x
    // read C) is a real mismatch under both OT and OB, so the freed conversion
    // cell must differ from it in all three scoring modes (+a, or 0 for NEUTRAL).
    const int conv_idx = ot ? (1 * 5 + 3)   // OT: ref C x read T
                            : (2 * 5 + 0);  // OB: ref G x read A
    REQUIRE(mat[conv_idx] != mat[0 * 5 + 1]);

    if (!batched_kswv_available()) {
        MESSAGE("tier has no batched kswv kernel; skipping");
        return;
    }
    if (!freed_kernel_available(mat)) {
        MESSAGE("tier has no freed-cell kswv kernel (needsScalar); skipping");
        return;
    }

    std::mt19937 rng(static_cast<unsigned>(seed));
    auto pairs = mixed_len_batch(rng, n_pairs, len_a, len_b);
    for (auto &p : pairs) convert_query_bisulfite(rng, p, ot);

    compare_batch(seed, pairs, mat, mat.data(), use16,
                  len_a > len_b ? len_a : len_b);
}

} // namespace

// The regression itself. 143 -> quantum 144 and 151 -> quantum 160, so the
// short lanes carry 16 padded query columns; 75 -> quantum 80 widens that to 80.
TEST_CASE("kswv::getScores8 matches scalar on mixed-query-length batches"
          * doctest::test_suite("unit/kswv")) {
    for (uint64_t seed = 1; seed <= 32; ++seed) {
        check_batch_matches_scalar(seed, 129, 143, 151);   // adjacent quanta
        check_batch_matches_scalar(seed, 129, 75, 151);    // wide spread
    }
}

// Same contract on the 16-bit kernels. These were never wrong, so this does not
// fail on origin/main -- it exists so the hoist cannot be introduced there
// later. 143/151 straddle a multiple of 8 as well as of 16.
TEST_CASE("kswv::getScores16 matches scalar on mixed-query-length batches"
          * doctest::test_suite("unit/kswv")) {
    for (uint64_t seed = 1; seed <= 32; ++seed) {
        check_batch_matches_scalar(seed, 129, 143, 151, /*use16=*/true);
        check_batch_matches_scalar(seed, 129, 75, 151, /*use16=*/true);
    }
}

// Control: uniform lengths leave no 0xFF query column on any lane whose result
// is compared, so this passed even before the fix -- the tail group's dummy
// lanes are padded to 0xFF from column 0, but their outputs are never read and
// lanes do not interact. It guards the common production path against a
// regression in the other direction.
TEST_CASE("kswv matches scalar on uniform-query-length batches"
          * doctest::test_suite("unit/kswv")) {
    for (uint64_t seed = 101; seed <= 104; ++seed) {
        check_batch_matches_scalar(seed, 129, 151, 151);
        check_batch_matches_scalar(seed, 129, 151, 151, /*use16=*/true);
    }
}

// The same contract under --meth, i.e. the HasFreed=true instantiation. This is
// the production bisulfite mate-rescue path (bwamem_pair.cpp:963, default ON via
// BWAMEM3_METH_BATCHED_RESCUE), and pre-fix it was worse off than the symmetric
// one on NEON: that kernel applied the hoisted rowboundary to BOTH
// instantiations, while AVX-512BW gated it (`cmp = HasFreed ? per-cell :
// rowboundary512`) and so kept meth correct. On NEON the two defects compounded
// -- active_frread's old 0xFF sentinel aliased the 0xFF query pad, so the freed
// blend fired on every padded column AND the hoisted mask failed to zero it.
// That is what FREED_INACTIVE8 and the restored per-cell test jointly fix, and
// nothing in-repo covered it: every other asymmetric-matrix test targets
// BandedPairWiseSW, not kswv.
//
// All three --meth-scoring modes, both strands: GENOMIC and NEUTRAL free one
// cell (to +a and to 0 respectively), COLLAPSED frees the mirror pair too, and
// the kernel folds all of them into one per-row active_frread target -- so each
// mode is a different fr_val/fr_ref2 configuration through the same blend.
TEST_CASE("kswv::getScores8 matches scalar on mixed-length --meth batches"
          * doctest::test_suite("unit/kswv")) {
    const int modes[] = {MEM_METH_SCORING_GENOMIC, MEM_METH_SCORING_NEUTRAL,
                         MEM_METH_SCORING_COLLAPSED};
    for (int scoring : modes) {
        for (bool ot : {true, false}) {
            for (uint64_t seed = 1; seed <= 8; ++seed) {
                check_meth_batch_matches_scalar(seed, 129, 143, 151, scoring, ot);
                check_meth_batch_matches_scalar(seed, 129, 75, 151, scoring, ot);
            }
        }
    }
}

// 16-bit counterpart. Never wrong here either -- the u16 sentinel 0x7FFF never
// aliased the 0xFFFF pad and all three 16-bit kernels kept the per-cell test --
// so this locks the invariant rather than reproducing a failure.
TEST_CASE("kswv::getScores16 matches scalar on mixed-length --meth batches"
          * doctest::test_suite("unit/kswv")) {
    const int modes[] = {MEM_METH_SCORING_GENOMIC, MEM_METH_SCORING_NEUTRAL,
                         MEM_METH_SCORING_COLLAPSED};
    for (int scoring : modes) {
        for (bool ot : {true, false}) {
            for (uint64_t seed = 1; seed <= 8; ++seed) {
                check_meth_batch_matches_scalar(seed, 129, 143, 151, scoring, ot,
                                                /*use16=*/true);
                check_meth_batch_matches_scalar(seed, 129, 75, 151, scoring, ot,
                                                /*use16=*/true);
            }
        }
    }
}

#endif // BWA_TESTS_HAVE_KSWV
