/* rb_group_key (src/rescue_band.h): the sort key RescueBandBatch::run_jobs groups band jobs by.
 *
 * Pass 1 groups by exact width, so two jobs of different widths in one width class must get
 * different keys and sort by width; pass 0 keeps the width-class key, so jobs in one class must tie
 * on width and sort by rows. In both modes the key must stay a non-negative int32_t at any width
 * and row count (job::key is an int32_t and run_jobs compares keys as signed), and the ordering is
 * width (or class) first, rows second. That the grouping leaves every output unchanged, and that
 * run_jobs applies the key to the right pass, is the band harness's to check (rescue_band_harness
 * eq with and without BWA3_RESCUE_BAND_P1_EXACT=0, and the grouping-counter check in CI). */
#include <cstdint>
#include <random>
#include <vector>

#include "doctest/doctest.h"
#include "rescue_band.h"

TEST_CASE("rb_group_key: pass 1 separates widths of one class, pass 0 does not"
          * doctest::test_suite("unit/rescue_band"))
{
    // 65 and 80 are in one class (65-80); 129 and 160 in another (129-160).
    const int pairs[][2] = {{65, 80}, {129, 160}, {17, 24}, {13, 16}};
    for (const auto &p : pairs) {
        CAPTURE(p[0]);
        CAPTURE(p[1]);
        REQUIRE(rb_width_bucket(p[0]) == rb_width_bucket(p[1]));
        // Same rows: the class key ties, the exact key orders by width.
        CHECK(rb_group_key(false, p[0], 100) == rb_group_key(false, p[1], 100));
        CHECK(rb_group_key(true, p[0], 100) < rb_group_key(true, p[1], 100));
        // Fewer rows on the wider job: the class key sorts by rows, the exact key still by width.
        CHECK(rb_group_key(false, p[1], 50) < rb_group_key(false, p[0], 100));
        CHECK(rb_group_key(true, p[0], 100) < rb_group_key(true, p[1], 50));
    }
}

TEST_CASE("rb_width_bucket: classes of 4 diagonals up to 16, of 8 up to 64, of 16 up to 128, of 32 beyond"
          * doctest::test_suite("unit/rescue_band"))
{
    // Class boundaries (each class's widest width), from the documented class sizes.
    std::vector<int> top;
    for (int w = 4; w <= 16; w += 4) top.push_back(w);
    for (int w = 24; w <= 64; w += 8) top.push_back(w);
    for (int w = 80; w <= 128; w += 16) top.push_back(w);
    for (int w = 160; w <= 1024; w += 32) top.push_back(w);
    int cls = 1;
    for (int w = 1; w <= 1024; w++) {
        CAPTURE(w);
        CHECK(rb_width_bucket(w) == cls);
        if (w == top[(size_t)cls - 1]) cls++;
    }
}

TEST_CASE("rb_group_key: bit layout, the class (pass 0) or the exact width (pass 1) over the rows"
          * doctest::test_suite("unit/rescue_band"))
{
    // A pin of the key's bit layout below the clamps; the ordering tests above and below state the
    // contract without the formula.
    for (int w = 1; w <= 400; w++)
        for (int r : {0, 1, 174, 65534, 65535}) {
            CAPTURE(w);
            CAPTURE(r);
            CHECK(rb_group_key(false, w, r) == ((rb_width_bucket(w) << 16) | r));
            CHECK(rb_group_key(true, w, r) == ((w << 16) | r));
        }
}

TEST_CASE("rb_group_key: width first, rows second, non-negative at any size"
          * doctest::test_suite("unit/rescue_band"))
{
    // Clamps: rows past 65535 and widths past 32767 share the top value instead of spilling into
    // the width bits or the sign bit.
    CHECK(rb_group_key(true, 5, 70000) == rb_group_key(true, 5, 65535));
    CHECK(rb_group_key(true, 5, 70000) < rb_group_key(true, 6, 0));
    CHECK(rb_group_key(true, 40000, 3) == rb_group_key(true, 32767, 3));
    CHECK(rb_group_key(true, 1 << 30, 1 << 30) >= 0);
    CHECK(rb_group_key(true, 32767, 65535) == INT32_MAX);
    // The class key too: its class passes 15 bits past a width of about 1M.
    CHECK(rb_group_key(false, 1 << 30, 1 << 30) >= 0);
    CHECK(rb_group_key(false, 2000000, 5) == ((32767 << 16) | 5));
    CHECK(rb_group_key(false, 1000000, 5) < rb_group_key(false, 2000000, 5));

    std::mt19937 rng(7);
    std::uniform_int_distribution<int> wd(1, 40000), rd(0, 70000);
    for (int i = 0; i < 20000; i++) {
        const int w1 = wd(rng), w2 = wd(rng), r1 = rd(rng), r2 = rd(rng);
        CAPTURE(w1);
        CAPTURE(w2);
        CAPTURE(r1);
        CAPTURE(r2);
        const int32_t a = rb_group_key(true, w1, r1), b = rb_group_key(true, w2, r2);
        REQUIRE(a >= 0);
        const int cw1 = w1 < 32767 ? w1 : 32767, cw2 = w2 < 32767 ? w2 : 32767;
        const int cr1 = r1 < 65535 ? r1 : 65535, cr2 = r2 < 65535 ? r2 : 65535;
        const bool lt = cw1 < cw2 || (cw1 == cw2 && cr1 < cr2);
        REQUIRE((a < b) == lt);
    }
}
