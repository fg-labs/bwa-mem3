// Unit tests for the two in-memory layouts of the occurrence checkpoints.
//
// The index file stores one 64-byte record per 64 BWT rows: four 64-bit counts
// and one one-hot word per base. The default layout reads those records as they
// are; the opt-in compact layout (--compact-occ) converts each to a 32-byte
// CP_OCC_PLANES record (32-bit counts and two bit planes), where the sentinel row
// reads as an A. These tests check every query the aligner makes against the
// file's own records (test/cp_file_oracle.h), at every row of an index:
//   - compact: occ_at / occ4_at for every base, the sentinel's block included,
//     and row ref_seq_len; the base at every row (the sentinel tested by index,
//     as callers do);
//   - both layouts: backwardExt, backwardExt_konly and forwardExt on intervals
//     that start at, end at and contain the sentinel row, and on every one-row
//     interval (the sentinel-row guards);
//   - both layouts: the SA resolver (get_sa_entry_compressed and call_one_step)
//     at every row, the sentinel's own row included, equal across layouts;
//   - the span-relative ("wide") counts a reference needs when some base
//     occurs 2^32 times or more, forced on this small index with a 512-row span.

#define private public
#include "FMI_search.h"
#undef private

#include "doctest/doctest.h"
#include "scoped_env.h"
#include "temp_index.h"
#include "lockstep_width.h"
#include "run_in_child.h"
#include "temp_dir.h"
#include "../cp_file_oracle.h"
#include "fm_index_writer.h"

#include <climits>
#include <unistd.h>   /* truncate */
#include <cstdio>
#include <algorithm>
#include <cstring>
#include <memory>
#include <tuple>
#include <random>
#include <string>
#include <vector>

namespace {

// Every occurrence query against the file's records; returns the number of
// checks made.
long check_all_rows(FMI_search &f, const std::vector<CP_OCC> &cp) {
    const int64_t L = f.reference_seq_len;
    long checks = 0;
    for (int64_t p = 0; p <= L; p++) {
        const CP_OCC_PLANES &blk = f.cp_planes[p >> CP_SHIFT];
        int64_t o4[4];
        f.occ4_at(blk, p, o4);
        for (int c = 0; c < 4; c++) {
            const int64_t want = cp_file_occ(cp, c, p);
            REQUIRE(f.occ_at(blk, p, c) == want);
            REQUIRE(o4[c] == want);
            checks += 2;
        }
        if (p < L) {
            const int b = cp_file_base(cp, p);
            if (p == f.sentinel_index) {
                REQUIRE(b == 4);
            } else {
                REQUIRE(FMI_search::cp_base_at(blk, p) == b);
            }
            checks++;
        }
    }
    return checks;
}

// backwardExt and backwardExt_konly against the oracle formula.
long check_backward(FMI_search &f, const std::vector<CP_OCC> &cp, int64_t sp, int64_t s) {
    const int64_t ep = sp + s, sent = f.sentinel_index;
    long checks = 0;
    for (uint8_t a = 0; a < 4; a++) {
        SMEM in{};
        in.k = sp; in.l = 1000; in.s = s;
        const int64_t occ_s = cp_file_occ(cp, a, sp);
        const int64_t want_k = f.count[a] + occ_s;
        const int64_t want_s = cp_file_occ(cp, a, ep) - occ_s;
        int64_t want_l = in.l + ((sp <= sent && ep > sent) ? 1 : 0);
        for (int b = a + 1; b < 4; b++) want_l += cp_file_occ(cp, b, ep) - cp_file_occ(cp, b, sp);
        const SMEM got = f.backwardExt(in, a);
        REQUIRE(got.k == want_k);
        REQUIRE(got.s == want_s);
        REQUIRE(got.l == want_l);
        const SMEM ko = f.backwardExt_konly(in, a);
        REQUIRE(ko.s == want_s);
        if (want_s > 0) REQUIRE(ko.k == want_k);   // k is left stale on s' == 0
        // forwardExt reads the interval at l: it must equal swapping k and l
        // around backwardExt by the complement, whose oracle is checked above.
        SMEM fw{};
        fw.k = 2000; fw.l = sp; fw.s = s;
        SMEM sw = fw;
        sw.k = fw.l; sw.l = fw.k;
        const SMEM b = f.backwardExt(sw, (uint8_t) (3 - a));
        const SMEM got_f = f.forwardExt(fw, a);
        REQUIRE(got_f.k == b.l);
        REQUIRE(got_f.l == b.k);
        REQUIRE(got_f.s == b.s);
        checks += 8;
    }
    return checks;
}

void check_index(FMI_search &f, const std::vector<CP_OCC> &cp) {
    const int64_t L = f.reference_seq_len, sent = f.sentinel_index;
    REQUIRE(sent >= 0);
    REQUIRE(sent < L);
    // occ_at / occ4_at / cp_base_at read the compact layout's records
    long checks = f.compact_occ() ? check_all_rows(f, cp) : 0;
    // every one-row interval (the konly one-row path reads the row's base)
    for (int64_t sp = 0; sp < L; sp++) checks += check_backward(f, cp, sp, 1);
    // intervals that start at, end at, or straddle the sentinel row, in and
    // across its block
    for (int64_t s : {2, 3, 7, 63, 64, 65, 130}) {
        for (int64_t d = -s; d <= 1; d++) {
            const int64_t sp = sent + d;
            if (sp < 0 || sp + s > L) continue;
            checks += check_backward(f, cp, sp, s);
        }
    }
    // multi-row intervals away from the sentinel too: s = 2 at every row, the
    // whole BWT, and a few hundred random intervals (sp = 0 and ep = L included)
    for (int64_t sp = 0; sp + 2 <= L; sp++) checks += check_backward(f, cp, sp, 2);
    checks += check_backward(f, cp, 0, L);
    std::mt19937_64 r((uint64_t)L * 1000003u + (uint64_t)sent);
    const int64_t smax = std::min<int64_t>(L, 300);
    for (int i = 0; i < 300; i++) {
        const int64_t s = 2 + (int64_t)(r() % (uint64_t)std::max<int64_t>(smax - 1, 1));
        if (s > L) continue;
        const int64_t sp = i == 0 ? 0 : i == 1 ? L - s : (int64_t)(r() % (uint64_t)(L - s + 1));
        checks += check_backward(f, cp, sp, s);
    }
    MESSAGE("L=" << L << " sentinel=" << sent << " checks=" << checks);
}

std::unique_ptr<FMI_search> load_layout(const std::string &prefix, bool compact) {
    std::unique_ptr<FMI_search> f(new FMI_search(prefix.c_str()));
    f->set_compact_occ(compact);
    f->load_index(/*load_pac=*/false, /*n_threads=*/2, /*build_seed_jump=*/false);
    REQUIRE(f->compact_occ() == compact);
    REQUIRE((f->cp_occ != NULL) == !compact);
    REQUIRE((f->cp_planes != NULL) == compact);
    return f;
}

// A compact load with 512-row spans and span-relative ("wide") counts forced on
// (FMI_search::set_cp_layout_for_test). The override is process-wide, so a guard
// undoes it on every exit, a failed REQUIRE's throw included: it must not stay on
// for the later tests in this binary.
struct ForcedWideLayout {
    ForcedWideLayout() { FMI_search::set_cp_layout_for_test(/*span_shift=*/9, /*force_wide=*/true); }
    ~ForcedWideLayout() { FMI_search::set_cp_layout_for_test(0, false); }
    ForcedWideLayout(const ForcedWideLayout &) = delete;
    ForcedWideLayout &operator=(const ForcedWideLayout &) = delete;
};
std::unique_ptr<FMI_search> load_forced_wide(const std::string &prefix) {
    ForcedWideLayout forced;
    return load_layout(prefix, /*compact=*/true);
}

// The SA resolver walks LF from every row to a sampled row (or the sentinel):
// the two layouts must resolve every row alike, under both sentinel contracts.
void check_sa_parity(FMI_search &dflt, FMI_search &cmp) {
    const int64_t L = dflt.reference_seq_len;
    REQUIRE(cmp.reference_seq_len == L);
    REQUIRE(cmp.sentinel_index == dflt.sentinel_index);
    for (int64_t r = 0; r < L; r++) {
        REQUIRE(cmp.get_sa_entry_compressed(r) == dflt.get_sa_entry_compressed(r));
        for (int drop = 0; drop <= 1; drop++) {
            int64_t a = 0, b = 0, oa = 0, ob = 0;
            const int64_t qa = dflt.call_one_step(r, a, oa, drop);
            const int64_t qb = cmp.call_one_step(r, b, ob, drop);
            REQUIRE(qa == qb);
            REQUIRE(a == b);
            REQUIRE(oa == ob);
        }
    }
    // the sentinel row itself resolves to the walk offset (0 from a fresh walk)
    int64_t v = -1, off = 3;
    if ((dflt.sentinel_index & dflt.sa_compx_mask) != 0) {
        REQUIRE(cmp.call_one_step(cmp.sentinel_index, v, off, 0) == 1);
        REQUIRE(v == 3);
        off = 3;
        REQUIRE(cmp.call_one_step(cmp.sentinel_index, v, off, 1) == 1);
        REQUIRE(v == 0);
    }
}

}  // namespace

TEST_CASE("CP_OCC_PLANES: every occurrence query matches the file's checkpoints"
          * doctest::test_suite("unit/fmi")) {
    std::mt19937 rng(20261116);
    // Several contigs, so the sentinel lands mid-block, and enough rows for a
    // few hundred blocks.
    bwa_tests::TempIndex base({{"c1", bwa_tests::random_acgt(rng, 3001)},
                               {"c2", bwa_tests::random_acgt(rng, 1517)},
                               {"c3", bwa_tests::random_acgt(rng, 2203)}});
    const std::vector<CP_OCC> cp = cp_file_records(base.prefix().c_str(), base.fmi()->reference_seq_len);
    std::unique_ptr<FMI_search> f = load_layout(base.prefix(), /*compact=*/true);
    CHECK_FALSE(f->cp_wide);
    check_index(*f, cp);
}

TEST_CASE("occurrence checkpoints: the default layout passes the same step and sentinel checks"
          * doctest::test_suite("unit/fmi")) {
    std::mt19937 rng(20261118);
    bwa_tests::TempIndex base({{"d1", bwa_tests::random_acgt(rng, 2901)},
                               {"d2", bwa_tests::random_acgt(rng, 1777)}});
    const std::vector<CP_OCC> cp = cp_file_records(base.prefix().c_str(), base.fmi()->reference_seq_len);
    std::unique_ptr<FMI_search> f = load_layout(base.prefix(), /*compact=*/false);
    check_index(*f, cp);
}

TEST_CASE("occurrence checkpoints: both layouts resolve every SA row alike"
          * doctest::test_suite("unit/fmi")) {
    std::mt19937 rng(20261119);
    bwa_tests::TempIndex base({{"s1", bwa_tests::random_acgt(rng, 2311)},
                               {"s2", bwa_tests::random_acgt(rng, 1409)},
                               {"s3", bwa_tests::random_acgt(rng, 997)}});
    std::unique_ptr<FMI_search> dflt = load_layout(base.prefix(), /*compact=*/false);
    std::unique_ptr<FMI_search> cmp  = load_layout(base.prefix(), /*compact=*/true);
    check_sa_parity(*dflt, *cmp);
}

TEST_CASE("CP_OCC_PLANES: span-relative counts (wide layout) match the file's checkpoints"
          * doctest::test_suite("unit/fmi")) {
    std::mt19937 rng(20261117);
    bwa_tests::TempIndex base({{"w1", bwa_tests::random_acgt(rng, 4099)},
                               {"w2", bwa_tests::random_acgt(rng, 2050)}});
    std::unique_ptr<FMI_search> f = load_forced_wide(base.prefix());
    REQUIRE(f->cp_wide);
    REQUIRE(f->cp_span_shift == 9);
    REQUIRE(f->cp_n_spans == (f->reference_seq_len >> 9) + 1);
    const std::vector<CP_OCC> cp = cp_file_records(base.prefix().c_str(), f->reference_seq_len);
    // the stored counts really are span-relative (else this test proves nothing)
    bool rebased = false;
    for (int64_t b = 0; b < (int64_t) cp.size() && !rebased; b++)
        for (int c = 0; c < 4; c++)
            if ((int64_t) f->cp_planes[b].cp_count[c] != cp[(size_t) b].cp_count[c]) rebased = true;
    REQUIRE(rebased);
    check_index(*f, cp);
    std::unique_ptr<FMI_search> dflt = load_layout(base.prefix(), /*compact=*/false);
    check_sa_parity(*dflt, *f);
}

TEST_CASE("occurrence checkpoints: the caller picks the layout; BWA3_COMPACT_OCC is read only on request"
          * doctest::test_suite("unit/fmi")) {
    struct Case { const char *val; bool compact; };
    /* An opt-in switch: on only for a leading '1'; "yes" is reported and ignored. */
    const Case cases[] = {{NULL, false}, {"", false}, {"0", false}, {"1", true}, {"yes", false}};
    bwa_tests::TempDir dir{"cp_env"};
    const std::string err = dir.file("err");
    for (const Case &c : cases) {
        // Restored (or unset again) at the end of each iteration, even when a check throws: the
        // whole unit suite also runs under BWA3_COMPACT_OCC=1.
        const bwa_tests::ScopedEnv env("BWA3_COMPACT_OCC", c.val);
        CAPTURE(c.val ? c.val : "(unset)");
        // In a child: the malformed-value warning is printed once per process.
        const int status = bwa_tests::run_in_child(
            [&] { _exit(FMI_search::compact_occ_env() ? 41 : 40); }, err.c_str());
        REQUIRE(WIFEXITED(status));
        CHECK(WEXITSTATUS(status) == (c.compact ? 41 : 40));
        const bool warned = dir.read("err").find("[W::compact_occ_env] ignoring BWA3_COMPACT_OCC=\"yes\" (not 0 or 1)")
                            != std::string::npos;
        CHECK(warned == (c.val != NULL && std::string(c.val) == "yes"));
        // The library never reads the variable: a new object starts in the default layout.
        FMI_search f("unused-prefix");
        CHECK_FALSE(f.compact_occ());
        f.set_compact_occ(true);   // the mem option, before load
        CHECK(f.compact_occ());
    }
}

TEST_CASE("occurrence checkpoints: set_compact_occ after load_index is refused"
          * doctest::test_suite("unit/fmi")) {
    std::mt19937 rng(20261126);
    bwa_tests::TempIndex base({{"g1", bwa_tests::random_acgt(rng, 701)}});
    bwa_tests::TempDir dir{"cp_guard"};
    const std::string err = dir.file("err");
    for (bool compact : {false, true}) {
        CAPTURE(compact);
        const int status = bwa_tests::run_in_child(
            [&] {
                FMI_search f(base.prefix().c_str());
                f.set_compact_occ(compact);
                f.load_index(/*load_pac=*/false, /*n_threads=*/1, /*build_seed_jump=*/false);
                f.set_compact_occ(!compact);
                _exit(0);
            },
            err.c_str());
        REQUIRE(WIFSIGNALED(status));
        CHECK(WTERMSIG(status) == SIGABRT);
        CHECK(dir.read("err").find("set_compact_occ called after load_index") != std::string::npos);
    }
}

TEST_CASE("CP_OCC_PLANES: cp_is_wide switches at a base count of 2^32"
          * doctest::test_suite("unit/fmi")) {
    for (int b = 0; b < 4; b++) {
        for (int64_t tot : {((int64_t)1 << 32) - 1, (int64_t)1 << 32}) {
            CAPTURE(b);
            CAPTURE(tot);
            int64_t raw[5] = {0, 0, 0, 0, 0};
            for (int i = 0; i < 4; i++) raw[i + 1] = raw[i] + (i == b ? tot : 1000);
            CHECK(FMI_search::cp_is_wide(raw) == (tot >= ((int64_t)1 << 32)));
        }
    }
}

TEST_CASE("CP_OCC_PLANES: the compact load refuses a corrupt count[]"
          * doctest::test_suite("unit/fmi")) {
    // The compact load picks its count width from count[b + 1] - count[b]
    // (cp_is_wide), which a corrupt count[] can overflow, so count[] must be
    // validated first. count[] follows reference_seq_len at the head of the
    // .bwt.2bit.64 file; these values each fail the raw bounds check.
    std::mt19937 rng(20261120);
    bwa_tests::TempIndex base({{"k1", bwa_tests::random_acgt(rng, 1201)}});
    const std::string bwt = base.prefix() + ".bwt.2bit.64";
    bwa_tests::TempDir dir{"cp_count"};
    const std::string err = dir.file("err");
    const int64_t bad[][2] = {{INT64_MIN, INT64_MAX}, {0, (int64_t)1 << 40}, {-1, 0}};
    for (const auto &b : bad) {
        CAPTURE(b[0]);
        CAPTURE(b[1]);
        const int status = bwa_tests::run_in_child(
            [&] {
                FILE *fp = fopen(bwt.c_str(), "r+b");
                if (fp == NULL || fseek(fp, (long)sizeof(int64_t), SEEK_SET) != 0
                    || fwrite(b, sizeof(int64_t), 2, fp) != 2 || fclose(fp) != 0)
                    _exit(90);
                FMI_search f(base.prefix().c_str());
                f.set_compact_occ(true);
                f.load_index(/*load_pac=*/false, /*n_threads=*/2, /*build_seed_jump=*/false);
            },
            err.c_str());
        REQUIRE(WIFSIGNALED(status));
        CHECK(WTERMSIG(status) == SIGABRT);
        CHECK(dir.read("err").find("count[] out of bounds") != std::string::npos);
    }
}

TEST_CASE("CP_OCC_PLANES: the compact load refuses corrupt checkpoint records"
          * doctest::test_suite("unit/fmi")) {
    // Each case corrupts one record of the .bwt.2bit.64 in place, loads it in the
    // compact layout in a child, and restores the bytes. Record i is at
    // FMI_BWT2BIT_HDR_BYTES + 64 i: cp_count[4] (int64), then one_hot[4] (uint64).
    std::mt19937 rng(20261127);
    bwa_tests::TempIndex base({{"x1", bwa_tests::random_acgt(rng, 1301)},
                               {"x2", bwa_tests::random_acgt(rng, 777)}});
    const int64_t L = base.fmi()->reference_seq_len, S = base.fmi()->sentinel_index;
    const std::vector<CP_OCC> cp = cp_file_records(base.prefix().c_str(), L);
    const std::string bwt = base.prefix() + ".bwt.2bit.64";
    // a row away from the sentinel, its block, its bit and its base's plane
    const int64_t row = (S + L / 2) % L == S ? (S + 1) % L : (S + L / 2) % L;
    REQUIRE(row != S);
    const int64_t blk = row >> CP_SHIFT;
    const uint64_t bit = 1ULL << (63 - (row & 63));
    const int c = cp_file_base(cp, row);
    REQUIRE(c < 4);
    const long rec = (long)FMI_BWT2BIT_HDR_BYTES + (long)blk * (long)sizeof(CP_OCC);
    struct Case { const char *what; long off; uint64_t val; const char *msg; };
    const Case cases[] = {
        {"two bases in one row", rec + 32 + 8 * ((c + 1) & 3),
         cp[(size_t)blk].one_hot_bwt_str[(c + 1) & 3] | bit, "not one-hot or out of range"},
        {"a non-sentinel row with no base", rec + 32 + 8 * c,
         cp[(size_t)blk].one_hot_bwt_str[c] & ~bit, "not one-hot or out of range"},
        {"a negative count", rec, (uint64_t)(int64_t)-1, "not one-hot or out of range"},
        {"a count past 32 bits on a narrow index", rec + 8, ((uint64_t)1 << 32) + 5, "not one-hot or out of range"},
    };
    bwa_tests::TempDir dir{"cp_records"};
    const std::string err = dir.file("err");
    for (const Case &k : cases) {
        CAPTURE(k.what);
        uint64_t orig = 0;
        FILE *fp = fopen(bwt.c_str(), "r+b");
        REQUIRE(fp != NULL);
        REQUIRE(fseek(fp, k.off, SEEK_SET) == 0);
        REQUIRE(fread(&orig, sizeof orig, 1, fp) == 1);
        REQUIRE(fseek(fp, k.off, SEEK_SET) == 0);
        REQUIRE(fwrite(&k.val, sizeof k.val, 1, fp) == 1);
        REQUIRE(fclose(fp) == 0);
        const int status = bwa_tests::run_in_child(
            [&] {
                FMI_search f(base.prefix().c_str());
                f.set_compact_occ(true);
                f.load_index(/*load_pac=*/false, /*n_threads=*/2, /*build_seed_jump=*/false);
                _exit(0);
            },
            err.c_str());
        fp = fopen(bwt.c_str(), "r+b");
        REQUIRE(fp != NULL);
        REQUIRE(fseek(fp, k.off, SEEK_SET) == 0);
        REQUIRE(fwrite(&orig, sizeof orig, 1, fp) == 1);
        REQUIRE(fclose(fp) == 0);
        REQUIRE(WIFSIGNALED(status));
        CHECK(WTERMSIG(status) == SIGABRT);
        CHECK(dir.read("err").find(k.msg) != std::string::npos);
    }
    // a record area cut short: the file-size check refuses it before the read
    {
        const std::string cut = dir.file("cut");
        bwa_tests::TempIndex t({{"x3", bwa_tests::random_acgt(rng, 1301)}});
        const std::string tb = t.prefix() + ".bwt.2bit.64";
        REQUIRE(truncate(tb.c_str(), (off_t)FMI_BWT2BIT_HDR_BYTES + 3 * (off_t)sizeof(CP_OCC) + 10) == 0);
        const int status = bwa_tests::run_in_child(
            [&] {
                FMI_search f(t.prefix().c_str());
                f.set_compact_occ(true);
                f.load_index(/*load_pac=*/false, /*n_threads=*/2, /*build_seed_jump=*/false);
                _exit(0);
            },
            cut.c_str());
        REQUIRE(WIFEXITED(status));
        CHECK(WEXITSTATUS(status) == EXIT_FAILURE);
        CHECK(dir.read("cut").find("has an invalid size, header or SA-rate tag") != std::string::npos);
    }
}

TEST_CASE("CP_OCC_PLANES: the compact load refuses counts that are not running totals"
          * doctest::test_suite("unit/fmi")) {
    // A record's count of base c is the c's in the rows before its block: zero
    // for the first, the previous record's plus its one-hot popcount for the
    // rest. The loader's workers each check their own range, slice to slice,
    // and the loader checks where one worker's range meets the next, and the
    // last record against the header's base totals. A
    // synthetic record area long enough for three workers, written out and read
    // back through cp_load_from_file (which returns false, so nothing aborts),
    // has its counts put out of step in turn. The shifted cases also break the
    // final totals; the middle-worker and single-record cases leave them intact,
    // so only the cross-worker or the within-worker check can catch those.
    const int64_t per_worker = (int64_t)(FMI_PREAD_MIN_CHUNK / sizeof(CP_OCC));
    const int64_t n = 3 * per_worker + 5;
    const int T = fmi_pread_worker_count((size_t)n * sizeof(CP_OCC), /*nthreads=*/4);
    REQUIRE(T == 3);
    const int64_t L = n * CP_BLOCK_SIZE - 10;   // a partial last block
    REQUIRE((L >> CP_SHIFT) + 1 == n);
    const int64_t sentinel = L / 3 + 7;
    std::mt19937_64 rng(20261127);
    std::vector<CP_OCC> recs((size_t)n);
    int64_t run[4] = {0, 0, 0, 0};
    for (int64_t b = 0; b < n; b++) {
        CP_OCC &r = recs[(size_t)b];
        const uint64_t hi = rng(), lo = rng();
        uint64_t live = ~0ULL;
        if (L - b * CP_BLOCK_SIZE < CP_BLOCK_SIZE) live = ~0ULL << (CP_BLOCK_SIZE - (L - b * CP_BLOCK_SIZE));
        if (b == sentinel >> CP_SHIFT) live &= ~(1ULL << (63 - (sentinel & 63)));
        r.one_hot_bwt_str[0] = ~hi & ~lo & live;
        r.one_hot_bwt_str[1] = ~hi & lo & live;
        r.one_hot_bwt_str[2] = hi & ~lo & live;
        r.one_hot_bwt_str[3] = hi & lo & live;
        for (int c = 0; c < 4; c++) {
            r.cp_count[c] = run[c];
            run[c] += __builtin_popcountll(r.one_hot_bwt_str[c]);
        }
    }
    // the header's raw C[] array: cumulative base counts, L - 1 in all (every
    // row but the sentinel's)
    const int64_t raw_count[5] = {0, run[0], run[0] + run[1], run[0] + run[1] + run[2],
                                  run[0] + run[1] + run[2] + run[3]};
    REQUIRE(raw_count[4] == L - 1);
    // the first record of each worker, as cp_load_from_file splits the area
    const int64_t worker1 = n / T + (0 < n % T ? 1 : 0);
    bwa_tests::TempDir dir{"cp_running"};
    const std::string path = dir.file("records");
    // CP_OCC_PLANES is 32-byte aligned, which std::vector does not honour
    // under C++14 (no aligned operator new): its buffer can sit 16 bytes off,
    // and code that stores whole records with aligned 32-byte moves (as gcc
    // emits for AVX2) then faults. Allocate it as load_index does.
    struct MmFree { void operator()(void *p) const { _mm_free(p); } };
    std::unique_ptr<CP_OCC_PLANES, MmFree> dst(
        (CP_OCC_PLANES *)_mm_malloc((size_t)n * sizeof(CP_OCC_PLANES), 64));
    REQUIRE(dst != nullptr);
    // Write `v` and load it in the compact layout's narrow form.
    auto load = [&](const std::vector<CP_OCC> &v, int64_t *empty_row) {
        FILE *fp = fopen(path.c_str(), "wb");
        REQUIRE(fp != NULL);
        REQUIRE(fwrite(v.data(), sizeof(CP_OCC), v.size(), fp) == v.size());
        REQUIRE(fclose(fp) == 0);
        FILE *in = fopen(path.c_str(), "rb");
        REQUIRE(in != NULL);
        const bool ok = FMI_search::cp_load_from_file(fileno(in), path.c_str(), 0, n, L, raw_count, dst.get(),
                                                      /*base=*/NULL, /*span_shift=*/32, /*wide=*/false,
                                                      /*nthreads=*/4, empty_row);
        fclose(in);
        return ok;
    };
    int64_t empty_row = -1;
    REQUIRE(load(recs, &empty_row));
    CHECK(empty_row == sentinel);
    // Add 1 to base c's count of every record from `from` on.
    auto shifted = [&](int64_t from, int c) {
        std::vector<CP_OCC> v = recs;
        for (int64_t b = from; b < n; b++) v[(size_t)b].cp_count[c] += 1;
        return v;
    };
    struct Case { const char *what; int64_t from; int c; };
    const Case cases[] = {
        {"the first record's counts are not zero", 0, 1},
        {"a record inside a worker's first slice", 100, 2},
        {"the first record of a worker's second slice", 16384, 0},
        {"the last record", n - 1, 1},
    };
    for (const Case &k : cases) {
        const std::string what = k.what;
        CAPTURE(what);
        CHECK_FALSE(load(shifted(k.from, k.c), &empty_row));
    }
    // Only the second worker's records, [worker1, 2 * worker1): each worker's
    // own range stays consistent and the last record's totals still match, so
    // only the cross-worker checks at its two boundaries can reject this file.
    std::vector<CP_OCC> middle = recs;
    for (int64_t b = worker1; b < 2 * worker1; b++) middle[(size_t)b].cp_count[3] += 1;
    CHECK_FALSE(load(middle, &empty_row));
    // A single record off by one at the first record of a worker's second
    // slice: only the within-worker check carried from slice to slice sees it.
    std::vector<CP_OCC> slice = recs;
    slice[16384].cp_count[0] += 1;
    CHECK_FALSE(load(slice, &empty_row));
    // a single record off by one, inside the third worker's range
    std::vector<CP_OCC> one = recs;
    one[(size_t)(2 * worker1 + 999)].cp_count[2] -= 1;
    CHECK_FALSE(load(one, &empty_row));
    // A row of the last record changed to another base: still one-hot, and no
    // record follows to disagree, so only the base totals can catch it.
    std::vector<CP_OCC> last = recs;
    CP_OCC &lr = last[(size_t)(n - 1)];
    int from_base = 0;
    while (lr.one_hot_bwt_str[from_base] == 0) from_base++;
    const uint64_t bit = lr.one_hot_bwt_str[from_base] & (0 - lr.one_hot_bwt_str[from_base]);
    lr.one_hot_bwt_str[from_base] &= ~bit;
    lr.one_hot_bwt_str[(from_base + 1) & 3] |= bit;
    CHECK_FALSE(load(last, &empty_row));
}

TEST_CASE("TempIndex: seed-jump construction is enabled by default and can be skipped"
          * doctest::test_suite("unit/fmi")) {
    bwa_tests::ScopedEnv env("BWA3_BWTSEED_JUMP", "2");
    std::mt19937 rng(20261127);
    const std::vector<bwa_tests::Contig> ref = {{"jump", bwa_tests::random_acgt(rng, 129)}};
    bwa_tests::TempIndex with_jump(ref);
    bwa_tests::TempIndex without_jump(ref, /*sa_compx=*/3, /*build_seed_jump=*/false);
    REQUIRE(with_jump.fmi()->seed_jump != nullptr);
    CHECK(with_jump.fmi()->seed_jump_depth == 2);
    CHECK(without_jump.fmi()->seed_jump == nullptr);
    CHECK(without_jump.fmi()->seed_jump_depth == 0);
    check_sa_parity(*with_jump.fmi(), *without_jump.fmi());
}

TEST_CASE("CP_OCC_PLANES: the sentinel and the reference end at block edges"
          * doctest::test_suite("unit/fmi")) {
    // The sentinel correction covers the rows after the sentinel in its block, and
    // the last block holds L & 63 real rows; random fixtures reach the edges of
    // both only by luck. Build small indexes until each edge has been checked:
    // L < 64 (the sentinel in block 0, which is also the last block), L & 63 ==
    // 63 and == 1 (L = 2n + 1 is odd, so a reference never ends on a block
    // multiple), and the sentinel as the first (S & 63 == 0) and last (== 63) row
    // of its block.
    std::mt19937 rng(20261125);
    bool small = false, l63 = false, l1 = false, s0 = false, s63 = false;
    int built = 0;
    for (int it = 0; it < 800 && !(small && l63 && l1 && s0 && s63); it++) {
        size_t n;
        if (!small) n = 25;                                        // L = 51
        else if (!l63) n = 31 + 32 * (1 + rng() % 4);              // 2n + 1 == 63 mod 64
        else if (!l1) n = 32 * (2 + rng() % 4);                    // 2n + 1 == 1 mod 64
        else n = 90 + rng() % 120;
        bwa_tests::TempIndex base({{"e", bwa_tests::random_acgt(rng, n)}},
                                 /*sa_compx=*/3, /*build_seed_jump=*/false);
        const int64_t L = base.fmi()->reference_seq_len, S = base.fmi()->sentinel_index;
        const bool hit = (L < 64 && !small) || ((L & 63) == 63 && !l63) || ((L & 63) == 1 && !l1)
                         || ((S & 63) == 0 && !s0) || ((S & 63) == 63 && !s63);
        if (!hit) continue;
        small |= L < 64; l63 |= (L & 63) == 63; l1 |= (L & 63) == 1;
        s0 |= (S & 63) == 0; s63 |= (S & 63) == 63;
        CAPTURE(L);
        CAPTURE(S);
        const std::vector<CP_OCC> cp = cp_file_records(base.prefix().c_str(), L);
        std::unique_ptr<FMI_search> f = load_layout(base.prefix(), /*compact=*/true);
        check_index(*f, cp);
        std::unique_ptr<FMI_search> d = load_layout(base.prefix(), /*compact=*/false);
        check_index(*d, cp);
        check_sa_parity(*d, *f);
        built++;
    }
    MESSAGE("edge fixtures checked: " << built);
    CHECK(small);
    CHECK(l63);
    CHECK(l1);
    CHECK(s0);
    CHECK(s63);
}

namespace {

// The seeding drivers' output on one batch: every SMEM's (rid, m, n, k, s),
// sorted. The interval start l is left out: the backward phase leaves it stale
// by design (backwardExt_konly), and nothing downstream reads it.
typedef std::tuple<uint32_t, uint32_t, uint32_t, int64_t, int64_t> SeedKey;

struct SeedBatch {
    std::vector<std::vector<char>> codes;
    std::vector<bseq1_t> seqs;
    std::vector<uint8_t> enc;
    std::vector<int32_t> cum;
    int32_t max_len = 0;
    explicit SeedBatch(const std::vector<std::vector<uint8_t>> &reads) {
        codes.resize(reads.size());
        seqs.resize(reads.size());
        for (size_t i = 0; i < reads.size(); ++i) {
            codes[i].assign(reads[i].begin(), reads[i].end());
            codes[i].push_back(0);
            std::memset(&seqs[i], 0, sizeof(bseq1_t));
            seqs[i].l_seq = (int)reads[i].size();
            seqs[i].seq = codes[i].data();
            cum.push_back((int32_t)enc.size());
            enc.insert(enc.end(), reads[i].begin(), reads[i].end());
            max_len = std::max(max_len, (int32_t)reads[i].size());
        }
        enc.push_back(4);
    }
    int32_t size() const { return (int32_t)seqs.size(); }
};

std::vector<SeedKey> seed_keys(const SMEM *a, int64_t n) {
    std::vector<SeedKey> v;
    for (int64_t i = 0; i < n; ++i) v.emplace_back(a[i].rid, a[i].m, a[i].n, a[i].k, a[i].s);
    std::sort(v.begin(), v.end());
    return v;
}

// Rounds 1, 2 and 3 of the seeding (getSMEMsAllPosOneThread,
// getSMEMsOnePosOneThread_lockstep at each read's middle with threshold t,
// bwtSeedStrategyAllPosOneThread and its lockstep twin) on `b`.
std::vector<std::vector<SeedKey>> seed_rounds(FMI_search &f, SeedBatch &b, int min_len, int32_t t,
                                              int32_t max_intv) {
    const int32_t nr = b.size();
    const size_t cap = (size_t)nr * (size_t)b.max_len + 1;
    std::vector<std::vector<SeedKey>> out;
    std::vector<int32_t> ones((size_t)nr, 1), rid((size_t)nr), qpos((size_t)nr), tv((size_t)nr, t),
        mi((size_t)nr, max_intv);
    for (int32_t i = 0; i < nr; ++i) {
        rid[(size_t)i] = i;
        qpos[(size_t)i] = b.seqs[(size_t)i].l_seq / 2;
    }
    std::vector<SMEM> buf(cap);
    int64_t n = 0;
    f.getSMEMsAllPosOneThread(b.enc.data(), ones.data(), rid.data(), nr, nr, b.seqs.data(), b.cum.data(),
                              b.max_len, min_len, buf.data(), &n);
    REQUIRE(n <= (int64_t)cap);
    out.push_back(seed_keys(buf.data(), n));
    // The all-position driver compacts rid in place; midpoint inputs are a new round.
    for (int32_t i = 0; i < nr; ++i) rid[(size_t)i] = i;
    n = 0;
    f.getSMEMsOnePosOneThread_lockstep(b.enc.data(), qpos.data(), tv.data(), rid.data(), nr, nr, b.seqs.data(),
                                       b.cum.data(), b.max_len, min_len, buf.data(), &n);
    REQUIRE(n <= (int64_t)cap);
    out.push_back(seed_keys(buf.data(), n));
    n = f.bwtSeedStrategyAllPosOneThread(b.enc.data(), mi.data(), nr, b.seqs.data(), b.cum.data(), min_len,
                                         buf.data());
    REQUIRE(n <= (int64_t)cap);
    out.push_back(seed_keys(buf.data(), n));
    n = f.bwtSeedStrategyAllPosOneThread_lockstep(b.enc.data(), mi.data(), nr, b.seqs.data(), b.cum.data(),
                                                  min_len, buf.data(), b.max_len);
    REQUIRE(n <= (int64_t)cap);
    out.push_back(seed_keys(buf.data(), n));
    return out;
}

char mutate_base(std::mt19937 &rng, char c) {
    static const char b[4] = {'A', 'C', 'G', 'T'};
    char d;
    do d = b[rng() & 3]; while (d == c);
    return d;
}

}  // namespace

// The lockstep widths and switch are process-wide: restored on every exit,
// a failed REQUIRE's throw included.
struct LockstepGlobals {
    const int32_t smem = g_smem_lockstep_n, bwt = g_bwtseed_lockstep_n, on = g_bwtseed_lockstep;
    LockstepGlobals() = default;
    ~LockstepGlobals() { g_smem_lockstep_n = smem; g_bwtseed_lockstep_n = bwt; g_bwtseed_lockstep = on; }
    LockstepGlobals(const LockstepGlobals &) = delete;
    LockstepGlobals &operator=(const LockstepGlobals &) = delete;
};

TEST_CASE("occurrence checkpoints: midpoint inputs are independent of all-position compaction"
          * doctest::test_suite("unit/fmi")) {
    std::mt19937 rng(20261126);
    const std::string ref = bwa_tests::random_acgt(rng, 1400);
    bwa_tests::TempIndex base({{"midpoint", ref}});
    SeedBatch b({bwa_tests::nt4("N"),
                 bwa_tests::nt4(ref.substr(100, 95) + "N" + ref.substr(400, 75)),
                 bwa_tests::nt4(ref.substr(700, 35) + "N" + ref.substr(900, 120))});
    auto dflt = load_layout(base.prefix(), false);
    auto narrow = load_layout(base.prefix(), true);
    auto wide = load_forced_wide(base.prefix());
    REQUIRE_FALSE(narrow->cp_wide);
    REQUIRE(wide->cp_wide);
    LockstepGlobals restore;
    size_t seeds_checked = 0;
    for (int32_t width : {1, 5, 16}) {
        g_smem_lockstep_n = width;
        for (FMI_search *f : {dflt.get(), narrow.get(), wide.get()}) {
            CAPTURE(width);
            CAPTURE(f->compact_occ());
            CAPTURE(f->cp_wide);
            const size_t cap = (size_t)b.size() * (size_t)b.max_len + 1;
            std::vector<SMEM> buf(cap);
            std::vector<int32_t> rid{0, 1, 2}, ones(3, 1);
            const std::vector<int32_t> identity = rid;
            int64_t n = 0;
            f->getSMEMsAllPosOneThread(b.enc.data(), ones.data(), rid.data(), b.size(), b.size(),
                                      b.seqs.data(), b.cum.data(), b.max_len, 8, buf.data(), &n);
            REQUIRE(n <= (int64_t)cap);
            REQUIRE(rid != identity);

            // Supply each original (read id, midpoint) independently: no state
            // from the all-position round participates in this oracle.
            std::vector<SeedKey> want;
            for (int32_t id = 0; id < b.size(); ++id) {
                int32_t read_id = id, midpoint = b.seqs[(size_t)id].l_seq / 2, threshold = 1;
                n = 0;
                f->getSMEMsOnePosOneThread_lockstep(b.enc.data(), &midpoint, &threshold, &read_id,
                                                   1, b.size(), b.seqs.data(), b.cum.data(),
                                                   b.max_len, 8, buf.data(), &n);
                REQUIRE(n <= (int64_t)cap);
                const auto keys = seed_keys(buf.data(), n);
                want.insert(want.end(), keys.begin(), keys.end());
            }
            std::sort(want.begin(), want.end());
            REQUIRE_FALSE(want.empty());
            const auto rounds = seed_rounds(*f, b, 8, 1, 20);
            CHECK(rounds[1] == want);
            seeds_checked += want.size();
        }
    }
    MESSAGE("independent midpoint seeds checked: " << seeds_checked);
    CHECK(seeds_checked > 0);
}

TEST_CASE("CP_OCC_PLANES: the seeding drivers give the same seeds on narrow and span-relative counts"
          * doctest::test_suite("unit/fmi")) {
    // The narrow compact layout (every real genome's, hg38 included) and the
    // wide (span-relative) one, each against the default layout.
    // The wide compact layout is only chosen when some base
    // occurs 2^32 times or more, so no fixture reaches it on its own: force it
    // with a 512-row span and run the three seeding rounds through it -- the
    // compact backward steps specialised for it (ls_backward_general<A, true>,
    // backwardExt_konly_ab<A, true>, cp_occ_t<A, true>), the run-time-base
    // single-interval step (ls_backward_one, backwardExt_konly_rt) and the
    // forward / full steps' wide branches -- comparing every seed with the
    // default layout's on the same index. The reference repeats a block with a
    // few diverged copies so that intervals hold many rows (the numPrev > 1
    // loops), and the reads are drawn from both strands with mismatches and Ns.
    std::mt19937 rng(20261121);
    const std::string block = bwa_tests::random_acgt(rng, 700);
    std::vector<bwa_tests::Contig> ref;
    for (int c = 0; c < 3; c++) {
        std::string seq = bwa_tests::random_acgt(rng, 900);
        for (int copy = 0; copy < 3; copy++) {
            std::string d = block;
            for (size_t i = (size_t)(37 * (c * 3 + copy)) % 61; i < d.size(); i += 61 + c + copy)
                d[i] = mutate_base(rng, d[i]);
            seq += d + bwa_tests::random_acgt(rng, 150 + (size_t)(rng() % 200));
        }
        ref.push_back({"r" + std::to_string(c), seq});
    }
    std::vector<std::vector<uint8_t>> reads;
    for (int r = 0; r < 120; r++) {
        const std::string &src = ref[(size_t)r % ref.size()].seq;
        const size_t len = 30 + (size_t)(rng() % 171);
        std::string q = src.substr((size_t)(rng() % (src.size() - len)), len);
        const int nmis = (int)(rng() % 5);
        for (int k = 0; k < nmis; k++) {
            const size_t at = (size_t)(rng() % len);
            q[at] = (rng() % 8 == 0) ? 'N' : mutate_base(rng, q[at]);
        }
        if (r & 1) q = bwa_tests::revcomp(q);
        reads.push_back(bwa_tests::nt4(q));
    }
    reads.push_back(bwa_tests::nt4(block.substr(100, 150)));   // a read in every copy of the block

    bwa_tests::TempIndex base(ref);
    std::unique_ptr<FMI_search> dflt = load_layout(base.prefix(), /*compact=*/false);
    std::unique_ptr<FMI_search> wide = load_forced_wide(base.prefix());
    REQUIRE(wide->cp_wide);
    REQUIRE(wide->cp_n_spans > 2);

    std::unique_ptr<FMI_search> narrow = load_layout(base.prefix(), /*compact=*/true);
    REQUIRE_FALSE(narrow->cp_wide);

    LockstepGlobals restore;   // a failed REQUIRE must not leak the widths into later tests
    SeedBatch b(reads);
    size_t n_seeds = 0;
    bool multi = false;
    for (int32_t width : {1, 5, 16}) {
        g_smem_lockstep_n = width;
        g_bwtseed_lockstep_n = width;
        g_bwtseed_lockstep = 1;
        for (int min_len : {1, 8, 19}) {
            CAPTURE(width);
            CAPTURE(min_len);
            const auto want = seed_rounds(*dflt, b, min_len, /*t=*/3, /*max_intv=*/20);
            const auto got = seed_rounds(*wide, b, min_len, /*t=*/3, /*max_intv=*/20);
            const auto got_n = seed_rounds(*narrow, b, min_len, /*t=*/3, /*max_intv=*/20);
            REQUIRE(want.size() == got.size());
            REQUIRE(want.size() == got_n.size());
            for (size_t round = 0; round < want.size(); round++) {
                CAPTURE(round);
                CHECK(got[round] == want[round]);
                CHECK(got_n[round] == want[round]);
                n_seeds += want[round].size();
                for (const SeedKey &k : want[round]) multi |= std::get<4>(k) > 1;
            }
        }
    }
    MESSAGE("seeds compared: " << n_seeds);
    CHECK(n_seeds > 0);
    CHECK(multi);   // some interval spans several rows, so the general loops ran
}
