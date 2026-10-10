// test/unit/test_global_cert_band.cpp
//
// Gate for bwa_global_cert_band (src/bwa.cpp): for two equal-length sequences,
// ksw_global2 must return the same score and CIGAR at the certified band kc as at
// any wider band, and at kc == 0 the result must be the ungapped alignment. That
// is what lets bwa_gen_cigar3 run the CIGAR DP at min(w, kc) (or skip it at 0).
//
// Pairs are equal-length targets and queries built from a random target with
// substitutions, N bases and balanced indel pairs (a deletion and an insertion
// elsewhere, so the lengths stay equal and gapped alignments do occur), scored
// with the default, random symmetric and random asymmetric matrices (an
// off-diagonal cell raised to the match score, as --meth does, or above it, so
// the bound must use the matrix's largest entry, not the match score) and random,
// independently drawn (so asymmetric) gap penalties. Each pair is aligned at
// several bands w, up to the full width; the result at every w must stay within
// kc of the diagonal and score at least the ungapped alignment, and whenever
// kc < w the result at kc (both the dispatched kernel and the scalar reference)
// must equal the result at w. The test also checks that it reaches certified
// bands below w, bands of 1 and 2 (below the wavefront kernels' floor), kc == 0,
// gapped optimal alignments at a narrowed band, and (on tiers that have one) the
// SIMD wavefront kernel at a narrowed band. A second case runs bwa_gen_cigar3
// itself on a packed reference, on both strands, against ksw_global2 at the band
// it computed before the certificate, checking score, CIGAR, NM and MD. A third
// pins the certificate's exact values and the scorings it declines.

#include <algorithm>
#include <cstdlib>
#include <random>
#include <string>
#include <vector>

#include "doctest/doctest.h"

#include "bntseq.h"   // _set_pac
#include "bwa.h"
#include "ksw.h"

namespace {

struct G2 {
    int score = 0, n_cigar = 0;
    std::vector<uint32_t> cigar;
    bool operator==(const G2 &o) const { return score == o.score && n_cigar == o.n_cigar && cigar == o.cigar; }
};

extern "C" {
using Ksw2Fn = int (*)(int, const uint8_t *, int, const uint8_t *, int, const int8_t *, int, int,
                       int, int, int, int *, uint32_t **);
}

G2 run(Ksw2Fn fn, const std::vector<uint8_t> &q, const std::vector<uint8_t> &t, const int8_t *mat,
       int o_del, int e_del, int o_ins, int e_ins, int w)
{
    G2 r;
    uint32_t *cig = nullptr;
    r.score = fn((int) q.size(), q.data(), (int) t.size(), t.data(), 5, mat, o_del, e_del, o_ins, e_ins,
                 w, &r.n_cigar, &cig);
    r.cigar.assign(cig, cig + r.n_cigar);
    free(cig);
    return r;
}

// Largest distance from the diagonal of any cell on the CIGAR's path.
int max_offset(const G2 &r)
{
    int off = 0, mx = 0;
    for (uint32_t c : r.cigar) {
        const int op = c & 0xf, len = (int) (c >> 4);
        if (op == 1) off += len; else if (op == 2) off -= len;
        mx = std::max(mx, std::abs(off));
    }
    return mx;
}

int ungapped_score(const int8_t *mat, const std::vector<uint8_t> &q, const std::vector<uint8_t> &t)
{
    int su = 0;
    for (size_t i = 0; i < q.size(); ++i) su += mat[t[i] * 5 + q[i]];
    return su;
}

// NM and MD of a CIGAR over the query and reference as bwa_gen_cigar3 sees them (after its
// reverse-strand reversal), following its documented contract: MD lists every literal
// mismatch, NM counts those (or, with nm_from_mat, those the matrix penalises) plus the gap
// lengths, and a deletion at either end is neither listed nor counted.
void nm_md(const G2 &r, const std::vector<uint8_t> &q, const std::vector<uint8_t> &t, bool rev,
           const int8_t *mat, int nm_from_mat, int *nm, std::string *md)
{
    const char *base = rev ? "TGCAN" : "ACGTN";
    int x = 0, y = 0, u = 0, n_mm = 0, n_gap = 0;
    md->clear();
    for (int k = 0; k < r.n_cigar; ++k) {
        const int op = r.cigar[k] & 0xf, len = (int) (r.cigar[k] >> 4);
        if (op == 0) {
            for (int i = 0; i < len; ++i) {
                if (q[x + i] != t[y + i]) {
                    *md += std::to_string(u) + base[t[y + i]];
                    u = 0;
                    n_mm += nm_from_mat ? (mat[t[y + i] * 5 + q[x + i]] < 0) : 1;
                } else ++u;
            }
            x += len; y += len;
        } else if (op == 2) {
            if (k > 0 && k < r.n_cigar - 1) {
                *md += std::to_string(u) + "^";
                for (int i = 0; i < len; ++i) *md += base[t[y + i]];
                u = 0; n_gap += len;
            }
            y += len;
        } else if (op == 1) x += len, n_gap += len;
    }
    *md += std::to_string(u);
    *nm = n_mm + n_gap;
}

} // namespace

TEST_CASE("bwa_global_cert_band: ksw_global2 at the certified band equals the wide band"
          * doctest::test_suite("unit/cigar")) {
    std::mt19937_64 g(0xC1A5B4ull);
    auto below = [&](int n) { return (int) (g() % (uint64_t) n); };
    auto unif = [&]() { return (g() >> 11) * (1.0 / 9007199254740992.0); };
    long checked = 0, narrowed = 0, k0 = 0, k12 = 0, gapped = 0, ungapped_dp = 0;
    unsigned long wave_narrow = 0, wave16_narrow = 0;   // wavefront-kernel runs at a narrowed band
    for (int it = 0; it < 2000; ++it) {
        // scoring
        int a = 1, b = 4, o_del = 6, e_del = 1, o_ins = 6, e_ins = 1;
        const int mode = it % 3;
        if (mode > 0) {
            a = 1 + below(3); b = 1 + below(6);
            o_del = below(9); e_del = 1 + below(3); o_ins = below(9); e_ins = 1 + below(3);
        }
        int8_t mat[25];
        bwa_fill_scmat(a, b, mat);
        if (mode == 2) {   // an off-diagonal cell at the match score, or above it
            const int r = below(4), c = (r + 1 + below(3)) & 3;
            mat[r * 5 + c] = (int8_t) (a + (it % 2 ? 0 : 1 + below(2)));
        }
        // an equal-length pair
        const int n = 20 + below(180);
        std::vector<uint8_t> t(n), q;
        for (auto &c : t) c = unif() < 0.02 ? 4 : (uint8_t) below(4);
        const double sub = unif() * 0.12;
        for (int i = 0; i < n; ++i) q.push_back(unif() < sub ? (uint8_t) ((t[i] + 1 + below(3)) & 3) : t[i]);
        for (auto &c : q) if (unif() < 0.01) c = 4;   // N in the query too
        const int npairs = below(3);
        for (int k = 0; k < npairs; ++k) {   // delete a run, insert a run of the same length elsewhere
            const int len = 1 + below(4), p = below(n - len), s = below(n - len);
            q.erase(q.begin() + p, q.begin() + p + len);
            for (int x = 0; x < len; ++x) q.insert(q.begin() + s, (uint8_t) below(4));
        }
        REQUIRE((int) q.size() == n);
        const int kc = bwa_global_cert_band(mat, o_del, e_del, o_ins, e_ins, n, q.data(), t.data());
        REQUIRE(kc >= 0);
        const int su = ungapped_score(mat, q, t);
        for (int w : {3, 1 + below(12), 8 + below(40), n}) {   // n: the full width
            const G2 wide = run(ksw_global2, q, t, mat, o_del, e_del, o_ins, e_ins, w);
            ++checked;
            // the certificate itself: the result at any band stays within kc of the diagonal
            CHECK(wide.score >= su);
            CHECK(max_offset(wide) <= kc);
            if (kc >= w) continue;
            ++narrowed;
            if (kc == 0) {
                ++k0;
                CHECK(wide.score == su);
                REQUIRE(wide.n_cigar == 1);
                CHECK(wide.cigar[0] == ((uint32_t) n << 4 | 0));
                continue;
            }
            if (kc < 3) ++k12;
            if (wide.n_cigar == 1) ++ungapped_dp; else ++gapped;
            const unsigned long wv = ksw_g2_wave_exec_count(), wv16 = ksw_g2_wave16_exec_count();
            const G2 cert = run(ksw_global2, q, t, mat, o_del, e_del, o_ins, e_ins, kc);
            wave_narrow += ksw_g2_wave_exec_count() - wv;
            wave16_narrow += ksw_g2_wave16_exec_count() - wv16;
            const G2 cert_ref = run(ksw_global2_scalar_ref, q, t, mat, o_del, e_del, o_ins, e_ins, kc);
            CHECK(cert == wide);
            CHECK(cert_ref == wide);
        }
    }
    MESSAGE("checked=" << checked << " narrowed=" << narrowed << " k0=" << k0 << " k12=" << k12
            << " gapped=" << gapped << " ungapped_dp=" << ungapped_dp
            << " wave_at_narrowed_band=" << wave_narrow << " wave16_at_narrowed_band=" << wave16_narrow);
    // the cases the test is for all occur
    CHECK(narrowed > checked / 4);
    CHECK(k0 > 100);
    CHECK(k12 > 100);           // bands below the wavefront kernels' floor, scalar only
    CHECK(gapped > 100);        // gapped optima compared at a narrowed band
    CHECK(ungapped_dp > 100);   // ungapped optima that still needed a band > 0
    if (ksw_g2_wave_wmin() > 0 || ksw_g2_wave16_wmin() > 0) CHECK(wave_narrow > 0);
    if (ksw_g2_wave16_wmin() > 0) CHECK(wave16_narrow > 0);
}

TEST_CASE("bwa_global_cert_band: bwa_gen_cigar3 matches ksw_global2 at its uncertified band"
          * doctest::test_suite("unit/cigar")) {
    // The integration: bwa_gen_cigar3 fetches the region from a packed reference (both
    // strands), applies the certificate (or the no-DP path at kc == 0) and must return the
    // score and CIGAR ksw_global2 gives at the band bwa_gen_cigar3 used before the
    // certificate, and the NM and MD those imply. Also exercised: --meth-like asymmetric
    // matrices with nm_from_mat, N bases in the read, a zero band (w_ == 0) and the
    // score-only call (n_cigar == NULL) the hit-merge test makes.
    std::mt19937_64 g(0x5EEDC16ull);
    auto below = [&](int n) { return (int) (g() % (uint64_t) n); };
    const int L = 20000;
    std::vector<uint8_t> ref(L), pac(L / 4 + 64, 0);
    for (int i = 0; i < L; ++i) { ref[i] = (uint8_t) below(4); _set_pac(pac.data(), i, ref[i]); }
    long checked = 0, k0[2] = {0, 0}, narrowed[2] = {0, 0};
    for (int it = 0; it < 2000; ++it) {
        int a = 1, b = 4, o_del = 6, e_del = 1, o_ins = 6, e_ins = 1;
        const int mode = it % 3;
        if (mode > 0) { a = 1 + below(3); b = 1 + below(6); o_del = below(9); e_del = 1 + below(3); o_ins = below(9); e_ins = 1 + below(3); }
        int8_t mat[25];
        bwa_fill_scmat(a, b, mat);
        int nm_from_mat = 0;
        if (mode == 2) {   // a --meth-like conversion cell (scored as a match, or 0), NM from the matrix
            const int r = below(4), c = (r + 1 + below(3)) & 3;
            mat[r * 5 + c] = (int8_t) (below(2) ? a : 0);
            nm_from_mat = 1;
        }
        const int n = 30 + below(170), rb = below(L - n - 1);
        // Forward strand: the read matches the region. Reverse strand (rb >= l_pac): it matches
        // the reverse complement; bwa_gen_cigar3 fetches that and reverses both sequences, so
        // ksw_global2 sees the reversed read against the complemented forward region.
        const bool rev = it % 2;
        std::vector<uint8_t> q(ref.begin() + rb, ref.begin() + rb + n);
        if (rev) {
            std::reverse(q.begin(), q.end());
            for (auto &c : q) c = (uint8_t) (3 - c);
        }
        for (int k = 0, nm = below(12); k < nm; ++k) q[below(n)] = (uint8_t) below(4);
        for (int k = 0, nn = below(3); k < nn; ++k) q[below(n)] = 4;   // N in the read
        if (below(2)) {   // a balanced indel pair
            const int len = 1 + below(3), p = below(n - len), s = below(n - len);
            q.erase(q.begin() + p, q.begin() + p + len);
            for (int x = 0; x < len; ++x) q.insert(q.begin() + s, (uint8_t) below(4));
        }
        const int w_ = it % 50 == 7 ? 0 : 1 + below(100);
        // The band bwa_gen_cigar3 computes before the certificate (its formula in bwa.cpp, at
        // equal lengths; keep in sync).
        int max_ins = (int) ((double) (((n + 1) >> 1) * mat[0] - o_ins) / e_ins + 1.);
        int max_del = (int) ((double) (((n + 1) >> 1) * mat[0] - o_del) / e_del + 1.);
        int max_gap = max_ins > max_del ? max_ins : max_del;
        max_gap = max_gap > 1 ? max_gap : 1;
        int w = (max_gap + 1) >> 1;
        w = w < w_ ? w : w_;
        w = w > 3 ? w : 3;
        std::vector<uint8_t> tw(ref.begin() + rb, ref.begin() + rb + n), qw = q, qb = q;
        if (rev) {
            for (auto &c : tw) c = (uint8_t) (3 - c);
            std::reverse(qw.begin(), qw.end());
        }
        G2 want;
        if (w_ == 0) {   // a zero band is emitted ungapped with no DP, as before the certificate
            want.score = ungapped_score(mat, qw, tw);
            want.n_cigar = 1;
            want.cigar.assign(1, (uint32_t) n << 4 | 0);
        } else {
            want = run(ksw_global2, qw, tw, mat, o_del, e_del, o_ins, e_ins, w);
        }
        int want_nm = 0;
        std::string want_md;
        nm_md(want, qw, tw, rev, mat, nm_from_mat, &want_nm, &want_md);

        int score = 0, n_cigar = 0, NM = 0;
        const int64_t qrb = rev ? 2 * (int64_t) L - (rb + n) : rb;
        uint32_t *cig = bwa_gen_cigar3(mat, o_del, e_del, o_ins, e_ins, w_, L, pac.data(), n, qb.data(),
                                       qrb, qrb + n, &score, &n_cigar, &NM, nm_from_mat);
        REQUIRE(cig != nullptr);
        G2 got;
        got.score = score; got.n_cigar = n_cigar; got.cigar.assign(cig, cig + n_cigar);
        const std::string got_md((const char *) (cig + n_cigar));
        free(cig);
        CHECK(got == want);
        CHECK(NM == want_nm);
        CHECK(got_md == want_md);
        CHECK(qb == q);   // the read is restored after the reverse-strand reversal

        // the score-only call (mem_patch_reg's shape): same score, no CIGAR
        int score_only = 0;
        uint32_t *none = bwa_gen_cigar3(mat, o_del, e_del, o_ins, e_ins, w_, L, pac.data(), n, qb.data(),
                                        qrb, qrb + n, &score_only, nullptr, nullptr, nm_from_mat);
        CHECK(none == nullptr);
        CHECK(score_only == want.score);

        // which path the input mix reaches (recomputed here; the checks above are what prove
        // bwa_gen_cigar3 took it correctly)
        const int kc = bwa_global_cert_band(mat, o_del, e_del, o_ins, e_ins, n, qw.data(), tw.data());
        ++checked;
        if (w_ > 0) { if (kc == 0) ++k0[rev]; else if (kc > 0 && kc < w) ++narrowed[rev]; }
    }
    MESSAGE("checked=" << checked << " k0(fwd,rev)=" << k0[0] << "," << k0[1]
            << " narrowed(fwd,rev)=" << narrowed[0] << "," << narrowed[1]);
    for (int s = 0; s < 2; ++s) {   // both strands reach the zero and the narrowed band
        CHECK(k0[s] > 50);
        CHECK(narrowed[s] > 50);
    }
}

TEST_CASE("bwa_global_cert_band: exact values, and the scorings it declines"
          * doctest::test_suite("unit/cigar")) {
    int8_t mat[25];
    bwa_fill_scmat(1, 4, mat);   // a = 1, b = 4
    const int n = 100;
    std::vector<uint8_t> t(n), q;
    for (int i = 0; i < n; ++i) t[i] = (uint8_t) (i * 7 % 4);
    // With -O 6 -E 1 on both sides the slope is a + e_ins + e_del = 3 and the numerator is
    // a*n - 12 - su = 5*mm - 12 for mm mismatches, so kc = floor((5*mm - 12) / 3), or 0.
    const int expect[] = {0, 0, 0, 1, 2, 4, 6, 7, 9};
    for (int mm = 0; mm <= 8; ++mm) {
        q = t;
        for (int k = 0; k < mm; ++k) q[k * 11] = (uint8_t) ((t[k * 11] + 1) & 3);
        CHECK(bwa_global_cert_band(mat, 6, 1, 6, 1, n, q.data(), t.data()) == expect[mm]);
    }
    q = t;
    // a negative gap open: the bound does not hold, so no certificate
    CHECK(bwa_global_cert_band(mat, -1, 1, 6, 1, n, q.data(), t.data()) == -1);
    CHECK(bwa_global_cert_band(mat, 6, 1, -1, 1, n, q.data(), t.data()) == -1);
    // amax + e_ins + e_del <= 0
    CHECK(bwa_global_cert_band(mat, 6, -1, 6, 0, n, q.data(), t.data()) == -1);
    int8_t neg[25];
    for (auto &v : neg) v = -1;
    CHECK(bwa_global_cert_band(neg, 6, 0, 6, 1, n, q.data(), t.data()) == -1);
    // zero opens are fine
    CHECK(bwa_global_cert_band(mat, 0, 1, 0, 1, n, q.data(), t.data()) == 0);
}
