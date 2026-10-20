// test/unit/test_ksw_global2_scalar.cpp
//
// Differential byte-identity test for the tuned scalar ksw_global2 path
// (ksw_global2_scalar_fast in src/ksw.cpp: per-thread buffers, the row loop,
// the register-window interior rows at w <= 5, and the run-length traceback)
// against the verbatim scalar reference it replaces (ksw_global2_scalar_ref).
//
// The public ksw_global2 sends a wide band to a wavefront kernel on SIMD tiers,
// so the tuned scalar is reached through its own test hook
// (ksw_global2_scalar_tuned) at every width. Every result is compared on score,
// n_cigar and every CIGAR word.
//
// Bands are kept feasible the way bwa_gen_cigar3 keeps them: w >= |tlen-qlen|+3
// for unequal lengths, any w >= 0 for equal lengths. Below that floor a
// traceback can leave the band and read direction bytes the fill never wrote,
// which the reference reads as uninitialized memory; there the tuned path is
// only checked for determinism (its own result, whatever its buffers held).

#include <cstdint>
#include <cstdlib>
#include <sstream>
#include <string>
#include <thread>
#include <vector>

#include "doctest/doctest.h"

#include "ksw.h"  // ksw_global2_scalar_ref, ksw_global2_scalar_tuned (extern "C")

namespace {

extern "C" {
using Ksw2Fn = int (*)(int, const uint8_t *, int, const uint8_t *, int,
                       const int8_t *, int, int, int, int, int, int *, uint32_t **);
}

struct G2 {
    int                   score;
    int                   n_cigar;
    std::vector<uint32_t> cigar;
    bool operator==(const G2 &o) const {
        return score == o.score && n_cigar == o.n_cigar && cigar == o.cigar;
    }
};

// Print a result as "score=S n_cigar=N cigar=10M2I..." so a mismatch shows
// which field differs.
std::ostream &operator<<(std::ostream &os, const G2 &r) {
    os << "score=" << r.score << " n_cigar=" << r.n_cigar << " cigar=";
    for (uint32_t c : r.cigar) os << (c >> 4) << ((c & 0xf) < 3 ? "MID"[c & 0xf] : '?');
    return os;
}

G2 run(Ksw2Fn fn, const std::vector<uint8_t> &q, const std::vector<uint8_t> &t, int m,
       const int8_t *mat, const int gaps[4], int w) {
    int       n_cigar = 0;
    uint32_t *cigar   = nullptr;
    const int score = fn(static_cast<int>(q.size()), q.data(), static_cast<int>(t.size()), t.data(),
                         m, mat, gaps[0], gaps[1], gaps[2], gaps[3], w, &n_cigar, &cigar);
    REQUIRE((n_cigar > 0) == (cigar != nullptr));
    G2 r{score, n_cigar, {}};
    if (n_cigar > 0) r.cigar.assign(cigar, cigar + n_cigar);
    free(cigar);
    return r;
}

// xorshift32, so the case stream is identical on every platform.
struct Rng {
    uint32_t x;
    uint32_t operator()(uint32_t n) {
        x ^= x << 13;
        x ^= x >> 17;
        x ^= x << 5;
        return x % n;
    }
};

// A bwa-style m x m matrix: match a, mismatch -b, the last code ambiguous.
std::vector<int8_t> bwa_mat(int m, int a, int b, int ambig) {
    std::vector<int8_t> mat(static_cast<size_t>(m) * m);
    for (int i = 0; i < m; ++i)
        for (int j = 0; j < m; ++j)
            mat[i * m + j] = static_cast<int8_t>(
                (i == m - 1 || j == m - 1) ? ambig : (i == j ? a : -b));
    return mat;
}

int feasible_floor(const std::vector<uint8_t> &q, const std::vector<uint8_t> &t) {
    return std::abs(static_cast<int>(t.size()) - static_cast<int>(q.size())) + 3;
}

// Check the tuned path against the reference: score, n_cigar and every word.
void check(const std::vector<uint8_t> &q, const std::vector<uint8_t> &t, int m,
           const std::vector<int8_t> &mat, const int gaps[4], int w, int c) {
    const G2 want = run(ksw_global2_scalar_ref, q, t, m, mat.data(), gaps, w);
    const G2 got  = run(ksw_global2_scalar_tuned, q, t, m, mat.data(), gaps, w);
    INFO("case " << c << " qlen=" << q.size() << " tlen=" << t.size() << " w=" << w
                 << " m=" << m << " gaps=" << gaps[0] << "," << gaps[1] << "," << gaps[2]
                 << "," << gaps[3]);
    CHECK(got == want);
}

}  // namespace

TEST_CASE("ksw_global2 tuned scalar matches the reference on random feasible bands"
          * doctest::test_suite("unit/ksw_global2")) {
    Rng rnd{0x5EED1234u};
    const std::vector<int8_t> mats[] = {bwa_mat(5, 1, 4, -1), bwa_mat(5, 127, 4, -1),
                                        bwa_mat(5, 5, 1, 0), bwa_mat(5, 2, 9, -3)};
    const int gaps[][4] = {{6, 1, 6, 1}, {5, 3, 7, 2}, {0, 1, 0, 1}, {40, 9, 25, 30}};
    for (int c = 0; c < 2000; ++c) {
        std::vector<uint8_t> t(1 + rnd(260));
        for (auto &b : t) b = static_cast<uint8_t>(rnd(rnd(30) == 0 ? 5 : 4));
        std::vector<uint8_t> q;
        for (uint8_t b : t) {
            if (rnd(12) == 0) continue;
            q.push_back(rnd(10) == 0 ? static_cast<uint8_t>(rnd(5)) : b);
        }
        if (rnd(7) == 0)
            for (uint32_t n = rnd(8); n > 0; --n) q.push_back(static_cast<uint8_t>(rnd(4)));
        if (q.empty()) q.push_back(0);
        const int fl = feasible_floor(q, t);
        int w;
        switch (rnd(4)) {
            case 0: w = q.size() == t.size() ? static_cast<int>(rnd(4)) : fl; break;
            case 1: w = fl + static_cast<int>(rnd(4)); break;
            case 2: w = fl + static_cast<int>(rnd(60)); break;
            default: w = fl + static_cast<int>(rnd(300)); break;
        }
        check(q, t, 5, mats[c % 4], gaps[(c / 4) % 4], w, c);
    }
}

TEST_CASE("ksw_global2 tuned scalar matches the reference on narrow bands"
          * doctest::test_suite("unit/ksw_global2")) {
    // w = 0..8 straddles the register-window cap (w <= 5); lengths around 2w + 3
    // move the interior row range [w + 1, qlen - w - 1) through empty, one row
    // and many rows; the largest gap costs push E and F below MINUS_INF without
    // overflowing int.
    Rng rnd{0xBAD5EED5u};
    const std::vector<int8_t> mats[] = {bwa_mat(5, 1, 4, -1), bwa_mat(5, 127, 4, -1),
                                        bwa_mat(5, 3, 1, 0)};
    const int gaps[][4] = {{6, 1, 6, 1}, {1, 1, 1, 1}, {5, 3, 7, 2},
                           {100000, 50000, 70000, 60000}, {0x3000000, 0x1000000, 0x2000000, 0x800000}};
    for (int c = 0; c < 8000; ++c) {
        int w = static_cast<int>(rnd(9));
        const uint32_t tl = rnd(3) == 0 ? 1 + rnd(4 * w + 8) : 1 + rnd(220);
        std::vector<uint8_t> t(tl);
        for (auto &b : t) b = static_cast<uint8_t>(rnd(rnd(40) == 0 ? 5 : 4));
        std::vector<uint8_t> q = t;
        for (auto &b : q)
            if (rnd(8) == 0) b = static_cast<uint8_t>(rnd(5));
        if (rnd(3) == 0) {
            // a one-base indel; the band keeps the feasibility floor
            const size_t at = rnd(static_cast<uint32_t>(q.size()));
            if (rnd(2) == 0 && q.size() > 1) q.erase(q.begin() + at);
            else q.insert(q.begin() + at, static_cast<uint8_t>(rnd(4)));
            w = w > feasible_floor(q, t) ? w : feasible_floor(q, t);
        }
        check(q, t, 5, mats[c % 3], gaps[(c / 3) % 5], w, c);
    }
}

TEST_CASE("ksw_global2 tuned scalar matches the reference on other alphabets and edge shapes"
          * doctest::test_suite("unit/ksw_global2")) {
    Rng rnd{0x0DDBA11u};
    const int gaps[4] = {6, 1, 6, 1};
    // A non-5 alphabet changes the query-profile stride.
    for (int m : {3, 4, 6}) {
        const std::vector<int8_t> mat = bwa_mat(m, 2, 3, -1);
        for (int c = 0; c < 500; ++c) {
            std::vector<uint8_t> t(1 + rnd(120));
            for (auto &b : t) b = static_cast<uint8_t>(rnd(m));
            std::vector<uint8_t> q = t;
            for (auto &b : q)
                if (rnd(9) == 0) b = static_cast<uint8_t>(rnd(m));
            check(q, t, m, mat, gaps, static_cast<int>(rnd(12)), c);
        }
    }
    // Single bases and tiny lengths, where n_col < 2w + 1 on every row.
    const std::vector<int8_t> mat = bwa_mat(5, 1, 4, -1);
    for (int ql = 1; ql <= 6; ++ql)
        for (int tl = 1; tl <= 6; ++tl)
            for (int w = 0; w <= 8; ++w) {
                std::vector<uint8_t> q(ql), t(tl);
                for (auto &b : q) b = static_cast<uint8_t>(rnd(5));
                for (auto &b : t) b = static_cast<uint8_t>(rnd(5));
                if (ql != tl && w < feasible_floor(q, t)) continue;
                check(q, t, 5, mat, gaps, w, ql * 100 + tl * 10 + w);
            }
}

TEST_CASE("ksw_global2 tuned scalar is unaffected by buffer reuse and size"
          * doctest::test_suite("unit/ksw_global2")) {
    // The first call's direction store exceeds the per-thread cap, so it
    // allocates its own buffers; the second fills the per-thread buffers with a
    // wide band's bytes and scores, and the small calls after it must not be
    // affected by what it left there.
    Rng rnd{0xC0FFEE1u};
    const std::vector<int8_t> mat = bwa_mat(5, 1, 4, -1);
    const int gaps[4] = {6, 1, 6, 1};
    std::vector<uint8_t> big_t(700), big_q;
    for (auto &b : big_t) b = static_cast<uint8_t>(rnd(4));
    for (uint8_t b : big_t)
        if (rnd(20) != 0) big_q.push_back(rnd(15) == 0 ? static_cast<uint8_t>(rnd(4)) : b);
    check(big_q, big_t, 5, mat, gaps, 400, -1);  // n_col * tlen > the per-thread cap: own buffers
    check(big_q, big_t, 5, mat, gaps, feasible_floor(big_q, big_t) + 20, -2);  // seeds the scratch
    for (int c = 0; c < 2000; ++c) {
        std::vector<uint8_t> t(1 + rnd(160));
        for (auto &b : t) b = static_cast<uint8_t>(rnd(5));
        std::vector<uint8_t> q = t;
        for (auto &b : q)
            if (rnd(10) == 0) b = static_cast<uint8_t>(rnd(5));
        check(q, t, 5, mat, gaps, static_cast<int>(rnd(20)), c);
    }
}

TEST_CASE("ksw_global2 tuned scalar matches the reference on shifted reads"
          * doctest::test_suite("unit/ksw_global2")) {
    // A read shifted against the target opens with s inserted (or deleted)
    // bases, so the alignment starts on the band's first row or first column,
    // where the fill seeds h from the gap costs (-(o_ins + e_ins*j) on row 0,
    // -(o_del + e_del*(i+1)) on column 0).
    Rng rnd{0x5A1F7EDu};
    const std::vector<int8_t> mat = bwa_mat(5, 1, 4, -1);
    const int gaps[][4] = {{6, 1, 6, 1}, {5, 3, 7, 2}, {1, 1, 1, 1}};
    for (int c = 0; c < 3000; ++c) {
        const int n = 8 + static_cast<int>(rnd(150));
        const int sh = 1 + static_cast<int>(rnd(8));
        std::vector<uint8_t> base(n + sh);
        for (auto &b : base) b = static_cast<uint8_t>(rnd(4));
        // Equal lengths: q = s new bases + t[0, n - s) (leading insertion,
        // trailing deletion), or the mirror image.
        std::vector<uint8_t> t(base.begin(), base.begin() + n), q(n);
        for (int j = 0; j < sh; ++j) q[j] = static_cast<uint8_t>(rnd(4));
        for (int j = sh; j < n; ++j) q[j] = t[j - sh];
        if (rnd(2) == 0) std::swap(q, t);
        const int w = sh + static_cast<int>(rnd(4)) - 1;
        check(q, t, 5, mat, gaps[c % 3], w < 0 ? 0 : w, c);
        // Unequal lengths at and just above the floor: the target gains the
        // shift on one end only.
        std::vector<uint8_t> t2(base.begin(), base.end()), q2(base.begin() + sh, base.end());
        if (rnd(2) == 0) std::swap(q2, t2);
        check(q2, t2, 5, mat, gaps[c % 3], feasible_floor(q2, t2) + static_cast<int>(rnd(3)), c);
    }
}

TEST_CASE("ksw_global2 tuned scalar is deterministic on infeasible bands"
          * doctest::test_suite("unit/ksw_global2")) {
    // Below the feasibility floor a traceback can read direction bytes the fill
    // never wrote. The tuned path zeroes them, so its result must not depend on
    // what the per-thread buffers held: run the cases in order on a fresh
    // thread, then again here with a wide-band call before each one that leaves
    // non-zero bytes in the buffers, and compare.
    Rng rnd{0xDE7E2u};
    const std::vector<int8_t> mat = bwa_mat(5, 1, 4, -1);
    const int gaps[4] = {6, 1, 6, 1};
    struct Case { std::vector<uint8_t> q, t; int w; };
    std::vector<Case> cases;
    for (int c = 0; c < 300; ++c) {
        std::vector<uint8_t> t(10 + rnd(120)), q;
        for (auto &b : t) b = static_cast<uint8_t>(rnd(4));
        for (uint8_t b : t)
            if (rnd(3) != 0) q.push_back(b);
        if (q.empty()) q.push_back(0);
        if (rnd(2) == 0) std::swap(q, t);
        const int d = std::abs(static_cast<int>(t.size()) - static_cast<int>(q.size()));
        cases.push_back({q, t, static_cast<int>(rnd(static_cast<uint32_t>(d > 1 ? d - 1 : 1)))});
    }
    auto run_all = [&](bool pollute, std::vector<G2> &out) {
        for (const Case &k : cases) {
            if (pollute) {
                // Fill the scratch with a wide band's direction bytes first.
                std::vector<uint8_t> pt(k.t.size() + 40), pq(k.q.size() + 40);
                for (auto &b : pt) b = static_cast<uint8_t>(rnd(4));
                for (auto &b : pq) b = static_cast<uint8_t>(rnd(4));
                (void) run(ksw_global2_scalar_tuned, pq, pt, 5, mat.data(), gaps, 200);
            }
            out.push_back(run(ksw_global2_scalar_tuned, k.q, k.t, 5, mat.data(), gaps, k.w));
        }
    };
    std::vector<G2> fresh, dirty;
    std::thread th([&] { run_all(false, fresh); });
    th.join();
    run_all(true, dirty);
    REQUIRE(fresh.size() == dirty.size());
    for (size_t c = 0; c < fresh.size(); ++c) {
        INFO("case " << c << " qlen=" << cases[c].q.size() << " tlen=" << cases[c].t.size()
                     << " w=" << cases[c].w);
        CHECK(fresh[c] == dirty[c]);
    }
}
