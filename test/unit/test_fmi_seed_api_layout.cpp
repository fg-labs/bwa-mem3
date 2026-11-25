// Unit tests for the checkpoint layout of the fmi_seed_api.h facade.
//
// fmi_seed_open must give the default layout whatever BWA3_COMPACT_OCC says (that
// variable selects `mem`'s layout, and the facade's contract predates it), so
// fmi_seed_cp_occ() is never NULL on its handle. The compact layout is reached only
// through fmi_seed_open_layout(..., FMI_SEED_LAYOUT_COMPACT); any other layout value
// is refused (NULL).

#include "fmi_seed_api.h"
#include "../cp_file_oracle.h"

#include "doctest/doctest.h"
#include "scoped_env.h"
#include "temp_index.h"

#include <random>
#include <string>

namespace {

struct Handle {
    FmiSeed *h;
    explicit Handle(FmiSeed *x) : h(x) {}
    ~Handle() { fmi_seed_close(h); }
    Handle(const Handle &) = delete;
    Handle &operator=(const Handle &) = delete;
};

void check_layout(FmiSeed *h, bool compact)
{
    REQUIRE(h != NULL);
    CHECK((fmi_seed_cp_occ(h) != NULL) == !compact);
    CHECK((fmi_seed_cp_occ_planes(h) != NULL) == compact);
    REQUIRE(fmi_seed_count(h) != NULL);
}

} // namespace

TEST_CASE("fmi_seed_open: the default layout whatever BWA3_COMPACT_OCC says"
          * doctest::test_suite("unit/fmi_seed_api_layout"))
{
    static_assert(FMI_SEED_API_VERSION >= 2, "fmi_seed_open_layout needs facade version 2");
    std::mt19937 rng(20261117);
    bwa_tests::TempIndex idx({{"c1", bwa_tests::random_acgt(rng, 2003)}, {"c2", bwa_tests::random_acgt(rng, 997)}});
    for (const char *env : {"", "0", "1"}) {
        CAPTURE(env);
        bwa_tests::ScopedEnv e("BWA3_COMPACT_OCC", env);
        Handle d(fmi_seed_open(idx.prefix().c_str()));
        check_layout(d.h, false);
        Handle dl(fmi_seed_open_layout(idx.prefix().c_str(), FMI_SEED_LAYOUT_DEFAULT));
        check_layout(dl.h, false);
        Handle c(fmi_seed_open_layout(idx.prefix().c_str(), FMI_SEED_LAYOUT_COMPACT));
        check_layout(c.h, true);
        /* Both layouts describe the same index. */
        for (int b = 0; b < 5; ++b) CHECK(fmi_seed_count(d.h)[b] == fmi_seed_count(c.h)[b]);
        CHECK(fmi_seed_sentinel(d.h) == fmi_seed_sentinel(c.h));
        /* a narrow index: the compact counts are not span-relative */
        CHECK(fmi_seed_cp_count_base(c.h) == NULL);
        CHECK(fmi_seed_cp_count_base(d.h) == NULL);
        /* an unknown layout fails closed */
        CHECK(fmi_seed_open_layout(idx.prefix().c_str(), 7) == NULL);
        CHECK(fmi_seed_open_layout(idx.prefix().c_str(), -1) == NULL);
    }
}

TEST_CASE("fmi_seed_cp_occ_planes: the compact records hold the file's counts and bases"
          * doctest::test_suite("unit/fmi_seed_api_layout"))
{
    std::mt19937 rng(20261128);
    bwa_tests::TempIndex idx({{"p1", bwa_tests::random_acgt(rng, 1709)}});
    Handle c(fmi_seed_open_layout(idx.prefix().c_str(), FMI_SEED_LAYOUT_COMPACT));
    REQUIRE(c.h != NULL);
    const int64_t L = idx.fmi()->reference_seq_len, S = fmi_seed_sentinel(c.h);
    const std::vector<CP_OCC> cp = cp_file_records(idx.prefix().c_str(), L);
    const CP_OCC_PLANES *pl = fmi_seed_cp_occ_planes(c.h);
    REQUIRE(pl != NULL);
    for (size_t b = 0; b < cp.size(); b++) {
        CAPTURE(b);
        for (int k = 0; k < 4; k++) CHECK((int64_t)pl[b].cp_count[k] == cp[b].cp_count[k]);
        for (int j = 0; j < 64; j++) {
            const int64_t row = (int64_t)b * 64 + j;
            if (row >= L || row == S) continue;   /* the sentinel and padding read as A */
            const int y = 63 - j;
            const int got = (int)(((pl[b].hi >> y) & 1) << 1 | ((pl[b].lo >> y) & 1));
            CHECK(got == cp_file_base(cp, row));
        }
    }
}
