// test/unit/test_chain_window_prefetch.cpp
//
// chain_window_pac_lines (src/chain_window_prefetch.h) picks the packed-reference lines the
// extension driver hints before fetching a chain's reference window. A wrong pick changes no
// output (the lines are only prefetched), so the bounds are checked here against a brute-force
// reference: clamp the window to the seed's strand, map each base to its forward-strand pac byte,
// and collect the distinct 64-byte lines those bytes sit in.

#include <cstdint>
#include <set>
#include <vector>

#include "doctest/doctest.h"
#include "chain_window_prefetch.h"

namespace {

// Every distinct line, in order, that the window's bases occupy in pac.
std::vector<uintptr_t> brute_lines(const uint8_t *pac, int64_t l_pac, int64_t rbeg,
                                   int64_t before, int64_t after)
{
    const bool rev = rbeg >= l_pac;
    const int64_t s_lo = rev ? l_pac : 0, s_hi = rev ? (l_pac << 1) : l_pac;
    std::set<uintptr_t> lines;
    uintptr_t prev = 0;
    bool any = false;
    for (int64_t x = rbeg - before; x < rbeg + after; ++x) {
        if (x < s_lo || x >= s_hi) continue;
        const int64_t fwd = rev ? (l_pac << 1) - 1 - x : x;
        const uintptr_t line = ((uintptr_t)pac + (uintptr_t)(fwd >> 2)) & ~(uintptr_t)63;
        // Consecutive bases mostly share a line; insert only when it changes (the set still
        // collects every distinct line), which keeps the sweep fast in an -O0 coverage build.
        if (!any || line != prev) lines.insert(line);
        prev = line;
        any = true;
    }
    return {lines.begin(), lines.end()};
}

void check_one(const uint8_t *pac, int64_t l_pac, int64_t rbeg, int64_t before, int64_t after,
               int max_lines)
{
    const std::vector<uintptr_t> want_all = brute_lines(pac, l_pac, rbeg, before, after);
    std::vector<uintptr_t> want = want_all;
    if ((int)want.size() > max_lines + 1) {
        want.assign(want_all.begin(), want_all.begin() + max_lines);
        want.push_back(want_all.back());
    }
    std::vector<uintptr_t> got(max_lines + 1, 0);
    const int n = chain_window_pac_lines(pac, l_pac, rbeg, before, after, max_lines, got.data());
    got.resize(n);
    CHECK_MESSAGE(got == want, "l_pac=" << l_pac << " rbeg=" << rbeg << " before=" << before
                                        << " after=" << after << " max_lines=" << max_lines
                                        << " pac%64=" << ((uintptr_t)pac & 63));
}

}  // namespace

TEST_CASE("chain_window_pac_lines: matches a brute-force line set on both strands" * doctest::test_suite("unit/extension")) {
    // pac is only addressed, never read, so a small buffer stands in for any l_pac; offsetting
    // its start covers every alignment of pac's base address.
    std::vector<uint8_t> buf(256);
    const int64_t l_pac = 20000;
    uint32_t x = 2463534242u;
    auto rnd = [&](uint32_t m) { x ^= x << 13; x ^= x >> 17; x ^= x << 5; return x % m; };
    for (int off : {0, 1, 17, 63}) {
        const uint8_t *pac = buf.data() + off;
        // Windows that touch both contig-strand ends, cross the strand boundary, or sit inside.
        for (int64_t rbeg : {int64_t(0), int64_t(5), int64_t(300), l_pac / 2, l_pac - 1, l_pac,
                             l_pac + 3, (l_pac << 1) - 400, (l_pac << 1) - 1}) {
            for (int64_t before : {0, 1, 75, 300, 2000}) {
                for (int64_t after : {1, 2, 150, 260, 3000}) {
                    for (int max_lines : {1, 2, 8}) check_one(pac, l_pac, rbeg, before, after, max_lines);
                }
            }
        }
        for (int t = 0; t < 500; ++t) {
            const int64_t rbeg = rnd((uint32_t)(l_pac << 1));
            check_one(pac, l_pac, rbeg, rnd(1200), 1 + rnd(1200), 1 + (int)rnd(10));
        }
    }
}

TEST_CASE("chain_window_pac_lines: max_lines below 1 hints only the window's last line" * doctest::test_suite("unit/extension")) {
    alignas(64) static uint8_t buf[256];
    const int64_t l_pac = 1000;
    for (int max_lines : {0, -1}) {
        uintptr_t out[2] = {0, 0};
        const int n = chain_window_pac_lines(buf, l_pac, 500, 300, 300, max_lines, out);
        CHECK(n == 1);
        CHECK(out[0] == (((uintptr_t)buf + (uintptr_t)(799 >> 2)) & ~(uintptr_t)63));
        CHECK(out[1] == 0);  // nothing written past the one slot
    }
}

TEST_CASE("chain_window_pac_lines: a window clamped away entirely hints nothing" * doctest::test_suite("unit/extension")) {
    std::vector<uint8_t> buf(64);
    uintptr_t out[9];
    const int64_t l_pac = 1000;
    // A forward seed whose window ends before base 0, and a reverse one that starts past 2*l_pac.
    CHECK(chain_window_pac_lines(buf.data(), l_pac, 10, 100, -50, 8, out) == 0);
    CHECK(chain_window_pac_lines(buf.data(), l_pac, (l_pac << 1) - 1, -60, 100, 8, out) == 0);
}

TEST_CASE("chain_window_pac_lines: a 150 bp read's window covers an end past the seed's neighbours" * doctest::test_suite("unit/extension")) {
    // A seed byte near the start of its line: the window's first byte sits in a line that the
    // seed's line and the lines 64 bytes either side of it do not reach, and it must be hinted.
    alignas(64) static uint8_t buf[4096];
    const int64_t l_pac = 16000;
    const int64_t rbeg = 7960;  // forward; seed byte 1990, window bytes 1915 .. 2064
    uintptr_t out[9];
    const int n = chain_window_pac_lines(buf, l_pac, rbeg, 300, 300, 8, out);
    const uintptr_t lo = ((uintptr_t)buf + (uintptr_t)((rbeg - 300) >> 2)) & ~(uintptr_t)63;
    const uintptr_t hi = ((uintptr_t)buf + (uintptr_t)((rbeg + 299) >> 2)) & ~(uintptr_t)63;
    REQUIRE(n >= 2);
    CHECK(out[0] == lo);
    CHECK(out[n - 1] == hi);
    CHECK(hi - lo == 192);  // four lines: 1856 .. 2048
    CHECK(n == 4);
    // The seed's byte and the bytes 64 either side of it (1926, 1990, 2054) miss the first line.
    const uintptr_t seed = (uintptr_t)buf + (uintptr_t)(rbeg >> 2);
    for (uintptr_t near : {seed - 64, seed, seed + 64}) CHECK((near & ~(uintptr_t)63) != lo);
}
