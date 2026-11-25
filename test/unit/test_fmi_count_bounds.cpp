// Unit tests for the bound on the FM-index's count[] (the cumulative C[]
// array): the raw counts leave out the one sentinel row, so a valid raw value
// is at most reference_seq_len - 1, and the +1-adjusted one at most
// reference_seq_len. A larger one lets an interval start at row
// reference_seq_len + 1, whose occurrence block lies one past cp_occ when
// reference_seq_len % 64 == 63. The shm scalar check, the disk loader and
// the stager must all reject it. The stager must also reject counts that
// decrease, which the shm scalar check would then refuse. Every row but the
// sentinel's holds a base, so all three also require the last count to be
// exactly reference_seq_len - 1 raw (reference_seq_len adjusted).

#include "doctest/doctest.h"

#include "FMI_search.h"
#include "bwa_shm.h"
#include "run_in_child.h"
#include "temp_index.h"

#include <stdint.h>
#include <string.h>
#include <sys/wait.h>

#include <cstdio>
#include <random>
#include <string>

namespace {

std::string random_acgt(std::mt19937 &rng, size_t n) {
    static const char kBases[] = "ACGT";
    std::string s(n, 'A');
    for (char &c : s) c = kBases[rng() & 3];
    return s;
}

// Overwrite the raw count[i] in <prefix>.bwt.2bit.64 (after reference_seq_len).
void write_raw_count(const std::string &prefix, int i, int64_t value) {
    FILE *f = fopen((prefix + ".bwt.2bit.64").c_str(), "r+b");
    REQUIRE(f != NULL);
    REQUIRE(fseek(f, (long)(sizeof(int64_t) * (1 + i)), SEEK_SET) == 0);
    REQUIRE(fwrite(&value, sizeof(value), 1, f) == 1);
    REQUIRE(fclose(f) == 0);
}

void write_raw_count4(const std::string &prefix, int64_t value) { write_raw_count(prefix, 4, value); }

// Load `prefix` from disk in a child; true if it loaded. *err (if given) gets
// what the child wrote to stderr.
bool loads_from_disk(const std::string &prefix, std::string *err = NULL) {
    const std::string log = prefix + ".load.log";
    const int status = bwa_tests::run_in_child([&] {
        FMI_search fmi(prefix.c_str());
        fmi.set_compact_occ(FMI_search::compact_occ_env());   // the compact-layout test legs
        fmi.load_index(/*load_pac=*/true, /*n_threads=*/1);
    }, log.c_str());
    if (err != NULL) {
        err->clear();
        if (FILE *f = fopen(log.c_str(), "rb")) {
            char buf[4096];
            size_t got;
            while ((got = fread(buf, 1, sizeof(buf), f)) > 0) err->append(buf, got);
            fclose(f);
        }
    }
    remove(log.c_str());
    return WIFEXITED(status) && WEXITSTATUS(status) == 0;
}

}  // namespace

TEST_CASE("shm_scalars_ok: count[] up to the reference length, not one past it"
          * doctest::test_suite("unit/fmi_count_bounds")) {
    // reference_seq_len, count[0..4] (+1-adjusted), sentinel_index, sa_compx.
    int64_t scalars[8] = {255, 1, 64, 128, 192, 255, 7, 3};
    int64_t sizes[3] = {0, 0, 0};
    char why[160] = "";
    CHECK(FMI_search::shm_scalars_ok((const uint8_t *)scalars, sizes, why, sizeof(why)));
    scalars[5] = 256;   // reference_seq_len + 1: row 256's block is cp_occ[4], one past
    CHECK_FALSE(FMI_search::shm_scalars_ok((const uint8_t *)scalars, sizes, why, sizeof(why)));
    CHECK(std::string(why).find("count[4]=256 out of bounds") != std::string::npos);
}

TEST_CASE("shm_scalars_ok: count[4] must be the reference length"
          * doctest::test_suite("unit/fmi_count_bounds")) {
    // In range and non-decreasing, but one row short of every row but the
    // sentinel's: T's interval would be one row too narrow.
    int64_t scalars[8] = {255, 1, 64, 128, 192, 254, 7, 3};
    int64_t sizes[3] = {0, 0, 0};
    char why[160] = "";
    CHECK_FALSE(FMI_search::shm_scalars_ok((const uint8_t *)scalars, sizes, why, sizeof(why)));
    CHECK(std::string(why).find("count[4]=254 != ref_seq_len=255") != std::string::npos);
}

TEST_CASE("bwa_shm_compute: a raw count[4] short of the reference length is refused"
          * doctest::test_suite("unit/fmi_count_bounds")) {
    std::mt19937 rng(8141);
    bwa_tests::TempIndex idx({{"chrA", random_acgt(rng, 300)}});
    const int64_t ref_len = idx.fmi()->reference_seq_len;
    bwa_shm_layout_t layout;
    write_raw_count4(idx.prefix(), ref_len - 2);
    CHECK(bwa_shm_compute(idx.prefix().c_str(), &layout) != 0);
    write_raw_count4(idx.prefix(), ref_len - 1);
    REQUIRE(bwa_shm_compute(idx.prefix().c_str(), &layout) == 0);
    bwa_shm_layout_free(&layout);
}

TEST_CASE("load_index: a raw count[4] short of the reference length fails closed"
          * doctest::test_suite("unit/fmi_count_bounds")) {
    // In [0, ref_seq_len - 1] and still at least count[3], so only the
    // count[4] == ref_seq_len - 1 check can refuse it.
    std::mt19937 rng(8142);
    bwa_tests::TempIndex idx({{"chrA", random_acgt(rng, 300)}});
    const int64_t ref_len = idx.fmi()->reference_seq_len;
    REQUIRE(idx.fmi()->count_data()[3] - 1 < ref_len - 2);   // the raw count[3]
    write_raw_count4(idx.prefix(), ref_len - 2);
    std::string err;
    CHECK_FALSE(loads_from_disk(idx.prefix(), &err));
    CHECK(err.find("count[4] != reference_seq_len - 1") != std::string::npos);
    write_raw_count4(idx.prefix(), ref_len - 1);
    CHECK(loads_from_disk(idx.prefix()));
}

TEST_CASE("load_index: a raw count[] reaching the reference length fails closed"
          * doctest::test_suite("unit/fmi_count_bounds")) {
    std::mt19937 rng(814);
    bwa_tests::TempIndex idx({{"chrA", random_acgt(rng, 300)}});
    const int64_t ref_len = idx.fmi()->reference_seq_len;
    const int64_t count4 = idx.fmi()->count_data()[4];   // +1-adjusted
    REQUIRE(count4 == ref_len);   // so the raw count[4] is ref_len - 1
    CHECK(loads_from_disk(idx.prefix()));
    write_raw_count4(idx.prefix(), ref_len);   // one row too many
    std::string err;
    CHECK_FALSE(loads_from_disk(idx.prefix(), &err));
    CHECK(err.find("count[] out of bounds") != std::string::npos);
    // Nor will `shm` stage it.
    bwa_shm_layout_t layout;
    CHECK(bwa_shm_compute(idx.prefix().c_str(), &layout) != 0);
    write_raw_count4(idx.prefix(), ref_len - 1);
    CHECK(loads_from_disk(idx.prefix()));
}

TEST_CASE("bwa_shm_compute: decreasing counts are refused, equal ones are not"
          * doctest::test_suite("unit/fmi_count_bounds")) {
    std::mt19937 rng(8140);
    bwa_tests::TempIndex idx({{"chrA", random_acgt(rng, 300)}});
    const int64_t raw2 = idx.fmi()->count_data()[2] - 1;   // undo the +1
    const int64_t raw4 = idx.fmi()->count_data()[4] - 1;
    REQUIRE(raw2 > 0);
    bwa_shm_layout_t layout;
    // Equal adjacent counts: a base with no occurrences (here T, with count[4]
    // still the reference length). Still stageable.
    write_raw_count(idx.prefix(), 3, raw4);
    REQUIRE(bwa_shm_compute(idx.prefix().c_str(), &layout) == 0);
    bwa_shm_layout_free(&layout);
    // count[3] below count[2]: in range, but not cumulative.
    write_raw_count(idx.prefix(), 3, raw2 - 1);
    CHECK(bwa_shm_compute(idx.prefix().c_str(), &layout) != 0);
}
