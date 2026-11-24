// test/unit/test_kswv_score2_bits.cpp
//
// The u8 kswv kernels' score2/te2 scan visits only the rows a per-lane bitmap
// marks as at or above minsc (kswv_score2.h). These cases check the bitmap
// builder bit for bit and the bitmap scan against the dense per-row walk the
// kernels ran before, at every u8 lane width (16, 32, 64): every minsc
// relation including 0 and 255, lanes whose len1 stops short of the rows the
// batch produced, row counts that end mid-byte, and low-entropy row maxima so
// that runs, ties inside a run and runs straddling the primary region occur.
// The end-to-end score2/te2 against scalar ksw_align2 stays with
// test_kswv_correctness.cpp and test_kswv_mixed_qlen.cpp.

#include <cstdint>
#include <random>
#include <vector>

#include "doctest/doctest.h"
#include "kswv_score2.h"

namespace {

// Reference dense walk: scalar ksw_u8's b[] build and score2 scan over every row
// of the lane, with no bitmap.
template <int W>
void dense_scan(const uint8_t *rowMax, int nrows, int l, int minsc, int low, int high,
                int *score2_out, int *te2_out)
{
    int score2 = -1, te2 = -1, b_score = -1, b_pos = -2;
    for (int i2 = 0; i2 < nrows; i2++) {
        int imax = (int)rowMax[i2 * W + l];
        if (imax < minsc) continue;
        if (b_pos + 1 != i2) {
            if (b_pos >= 0 && (b_pos < low || b_pos > high) && b_score > score2) {
                score2 = b_score;
                te2    = b_pos;
            }
            b_score = imax;
            b_pos   = i2;
        } else if (b_score < imax) {
            b_score = imax;
            b_pos   = i2;
        }
    }
    if (b_pos >= 0 && (b_pos < low || b_pos > high) && b_score > score2) {
        score2 = b_score;
        te2    = b_pos;
    }
    *score2_out = score2;
    *te2_out    = te2;
}

template <int W>
void check_width(uint64_t seed, int trials)
{
    std::mt19937_64 rng(seed);
    auto rnd = [&](int m) { return (int)(rng() % (uint64_t)m); };
    long found = 0, visited_lanes = 0;
    for (int t = 0; t < trials; t++) {
        const int processed = rnd(300);
        // Low-entropy row maxima: a narrow value range makes runs, equal-imax
        // ties and run boundaries common; occasionally the full byte range.
        const int base = rnd(256);
        const int span = rnd(4) == 0 ? 256 : 1 + rnd(6);
        // Extra rows past `processed` hold junk the scan must never visit.
        std::vector<uint8_t> rowMax((size_t)(processed + 9) * W);
        for (auto &v : rowMax) v = (uint8_t)((base + rnd(span)) & 0xFF);

        uint8_t minsc[W];
        for (int l = 0; l < W; l++) {
            const int c = rnd(8);
            minsc[l] = c == 0 ? 0 : c == 1 ? 255 : (uint8_t)((base + rnd(span + 2) - 1) & 0xFF);
        }
        std::vector<uint8_t> bits((size_t)kswv_score2_bits_blocks(processed) * W + W, 0xA5);
        kswv_score2_row_bits_u8<W>(rowMax.data(), processed, minsc, bits.data());

        // The builder: bit k of block b set iff row 8b+k < processed clears minsc;
        // nothing past the last block is written. Mismatches are counted and the
        // first one reported, so a trial costs one assertion.
        long bad = 0;
        int bad_l = -1, bad_r = -1;
        for (int l = 0; l < W; l++) {
            for (int r = 0; r < kswv_score2_bits_blocks(processed) * 8; r++) {
                const bool bit = (bits[(size_t)(r >> 3) * W + l] >> (r & 7)) & 1;
                const bool want = r < processed && rowMax[(size_t)r * W + l] >= minsc[l];
                if (bit != want && bad++ == 0) { bad_l = l; bad_r = r; }
            }
        }
        for (int k = 0; k < W; k++)
            if (bits[(size_t)kswv_score2_bits_blocks(processed) * W + k] != 0xA5 && bad++ == 0) {
                bad_l = k; bad_r = -2;   // wrote past the last block
            }
        REQUIRE_MESSAGE(bad == 0, "W ", W, " builder: ", bad, " bad bits, first lane ", bad_l,
                        " row ", bad_r, " (processed ", processed, ")");

        // A superset bitmap (every row a candidate, as BWA3_RESCUE_SCORE2_BITS=0 builds
        // it) must give the same result: the scan re-tests each visited row against minsc.
        std::vector<uint8_t> all_bits(bits.size(), 0xA5);
        kswv_score2_row_bits_all<W>(processed, all_bits.data());
        const size_t all_n = (size_t)kswv_score2_bits_blocks(processed) * W;
        for (size_t k = 0; k < all_bits.size(); k++)
            if (all_bits[k] != (k < all_n ? 0xFF : 0xA5) && bad++ == 0)
                bad_l = (int)(k % W);   // the all-marked fill missed a byte or overran
        REQUIRE_MESSAGE(bad == 0, "W ", W, " all-marked fill: ", bad, " bad bytes, first lane ",
                        bad_l, " (processed ", processed, ")");

        for (int l = 0; l < W; l++) {
            const int len1 = rnd(processed + 20);
            const int nrows = processed < len1 ? processed : len1;
            const int te = rnd(processed + 1);
            const int val = rnd(12);
            const int low = te - val, high = te + val;
            int s_dense, t_dense, s_bits, t_bits, s_all, t_all;
            dense_scan<W>(rowMax.data(), nrows, l, minsc[l], low, high, &s_dense, &t_dense);
            kswv_score2_scan_bits_u8<W>(rowMax.data(), bits.data(), nrows, l, minsc[l], low,
                                        high, &s_bits, &t_bits);
            kswv_score2_scan_bits_u8<W>(rowMax.data(), all_bits.data(), nrows, l, minsc[l], low,
                                        high, &s_all, &t_all);
            if ((s_bits != s_dense || t_bits != t_dense || s_all != s_dense || t_all != t_dense)
                && bad++ == 0) {
                bad_l = l; bad_r = nrows;
            }
            visited_lanes++;
            found += s_dense >= 0;
        }
        REQUIRE_MESSAGE(bad == 0, "W ", W, " scan: ", bad, " lanes differ from the dense walk, "
                        "first lane ", bad_l, " nrows ", bad_r, " (processed ", processed, ")");
    }
    // The draws must exercise both outcomes, or the comparison proves little.
    CHECK(found > visited_lanes / 10);
    CHECK(found < visited_lanes - visited_lanes / 10);
}

}  // namespace

TEST_CASE("kswv score2 bitmap scan equals the dense row walk at every u8 lane width"
          * doctest::test_suite("unit/kswv"))
{
    check_width<16>(0x5c02e2b175ull, 600);
    check_width<32>(0x5c02e2b176ull, 300);
    check_width<64>(0x5c02e2b177ull, 150);
}

TEST_CASE("kswv score2 bitmap scan: hand-built runs"
          * doctest::test_suite("unit/kswv"))
{
    // One lane (W = 16, lane 3), minsc 10. Rows 0..19:
    //   run A: rows 2-4 (12, 15, 15) -> anchor row 3 (first strict max), 15
    //   run B: rows 7-8 (20, 20)    -> anchor row 7, 20  (inside [6, 9]: excluded)
    //   run C: rows 15-17 (9 row 15 below minsc, 11, 18) -> rows 16-17, anchor 17, 18
    // Expected: score2 = 18 at te2 = 17 (B excluded, C beats A).
    constexpr int W = 16;
    const int nrows = 20, lane = 3;
    std::vector<uint8_t> rowMax((size_t)nrows * W, 0);
    const int vals[nrows] = {0, 9, 12, 15, 15, 0, 1, 20, 20, 0, 0, 0, 0, 0, 0, 9, 11, 18, 0, 0};
    for (int r = 0; r < nrows; r++) rowMax[(size_t)r * W + lane] = (uint8_t)vals[r];
    uint8_t minsc[W];
    for (int l = 0; l < W; l++) minsc[l] = 10;
    std::vector<uint8_t> bits((size_t)kswv_score2_bits_blocks(nrows) * W);
    kswv_score2_row_bits_u8<W>(rowMax.data(), nrows, minsc, bits.data());
    CHECK(bits[0 * W + lane] == 0x9C);  // rows 2,3,4,7
    CHECK(bits[1 * W + lane] == 0x01);  // row 8
    CHECK(bits[2 * W + lane] == 0x03);  // rows 16,17

    int s2, t2;
    kswv_score2_scan_bits_u8<W>(rowMax.data(), bits.data(), nrows, lane, 10, 6, 9, &s2, &t2);
    CHECK(s2 == 18);
    CHECK(t2 == 17);

    // A lane whose len1 stops before run C keeps run A.
    kswv_score2_scan_bits_u8<W>(rowMax.data(), bits.data(), 16, lane, 10, 6, 9, &s2, &t2);
    CHECK(s2 == 15);
    CHECK(t2 == 3);

    // No rows: nothing qualifies.
    kswv_score2_scan_bits_u8<W>(rowMax.data(), bits.data(), 0, lane, 10, 6, 9, &s2, &t2);
    CHECK(s2 == -1);
    CHECK(t2 == -1);
}
