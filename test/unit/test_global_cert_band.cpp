// test/unit/test_global_cert_band.cpp
//
// Gate for bwa_global_cert_band (src/bwa.cpp): ksw_global2 must return the same
// score and CIGAR at the certified band kc as at any wider band, and for two
// equal-length sequences at kc == 0 the result must be the ungapped alignment.
// That is what lets bwa_gen_cigar3 run the CIGAR DP at min(w, kc) (or skip it at
// 0), and, for unequal lengths, at min(w, max(kc, |tlen - qlen| + 3)). The first
// three cases below cover equal lengths; the last four cover unequal ones (band
// equivalence, the bwa_gen_cigar3 integration on both strands, exact values, and
// the int32 wavefront kernel at a narrowed band).
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
//
// The unequal-length cases build the query from the reference with indels whose
// net length change is d = tlen - qlen (one run or two of the same kind, plus
// sometimes a balanced pair, so optima that leave the corridor between the two
// diagonals occur), under the same scorings. The result at every band must stay
// within kc - |d| of the corridor, and the result at max(kc, |d| + 3) (and, on the
// scalar reference, at kc itself) must equal the result at any wider band. The
// integration case checks bwa_gen_cigar3 on both strands for both gap kinds,
// including, through its -v 4 trace, that it really ran at the narrowed band.
// Exact values cover asymmetric gap opens and gaps at either end, and a case
// with a large match score drives the int32 wavefront kernel at a narrowed band.

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <unistd.h>
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
        const int kc = bwa_global_cert_band(mat, o_del, e_del, o_ins, e_ins, n, q.data(), n, t.data());
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
        const int kc = bwa_global_cert_band(mat, o_del, e_del, o_ins, e_ins, n, qw.data(), n, tw.data());
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
        CHECK(bwa_global_cert_band(mat, 6, 1, 6, 1, n, q.data(), n, t.data()) == expect[mm]);
    }
    q = t;
    // a negative gap open: the bound does not hold, so no certificate
    CHECK(bwa_global_cert_band(mat, -1, 1, 6, 1, n, q.data(), n, t.data()) == -1);
    CHECK(bwa_global_cert_band(mat, 6, 1, -1, 1, n, q.data(), n, t.data()) == -1);
    // amax + e_ins + e_del <= 0
    CHECK(bwa_global_cert_band(mat, 6, -1, 6, 0, n, q.data(), n, t.data()) == -1);
    int8_t neg[25];
    for (auto &v : neg) v = -1;
    CHECK(bwa_global_cert_band(neg, 6, 0, 6, 1, n, q.data(), n, t.data()) == -1);
    // zero opens are fine
    CHECK(bwa_global_cert_band(mat, 0, 1, 0, 1, n, q.data(), n, t.data()) == 0);
}

namespace {

// The band bwa_gen_cigar3 computes before the certificate (its formula in bwa.cpp; keep in
// sync), for a query of length n and a reference region of length m.
int pre_cert_band(const int8_t *mat, int o_del, int e_del, int o_ins, int e_ins, int n, int m, int w_)
{
    const int dl = std::abs(m - n);
    int max_ins = (int) ((double) (((n + 1) >> 1) * mat[0] - o_ins) / e_ins + 1.);
    int max_del = (int) ((double) (((n + 1) >> 1) * mat[0] - o_del) / e_del + 1.);
    int max_gap = max_ins > max_del ? max_ins : max_del;
    max_gap = max_gap > 1 ? max_gap : 1;
    int w = (max_gap + dl + 1) >> 1;
    w = w < w_ ? w : w_;
    return w > dl + 3 ? w : dl + 3;
}

// A random edit of t: substitutions, N bases, and indels whose net length change is dl (an
// extra balanced pair sometimes, so alignments with both an insertion and a deletion occur).
std::vector<uint8_t> edit_unequal(std::mt19937_64 &g, const std::vector<uint8_t> &t, int dl, double sub)
{
    auto below = [&](int n) { return (int) (g() % (uint64_t) n); };
    auto unif = [&]() { return (g() >> 11) * (1.0 / 9007199254740992.0); };
    const int n = (int) t.size();
    std::vector<uint8_t> q;
    for (int i = 0; i < n; ++i) q.push_back(unif() < sub ? (uint8_t) ((t[i] + 1 + below(3)) & 3) : t[i]);
    for (auto &c : q) if (unif() < 0.01) c = 4;
    // the net change: one run, or split into two runs of the same kind
    const int runs = std::abs(dl) > 1 && below(3) == 0 ? 2 : 1;
    for (int r = 0, left = std::abs(dl); r < runs; ++r) {
        const int len = r + 1 == runs ? left : 1 + below(left - 1);
        left -= len;
        if (dl > 0) {   // the reference is longer: delete from the query
            const int p = below((int) q.size() - len);
            q.erase(q.begin() + p, q.begin() + p + len);
        } else {
            const int p = below((int) q.size() + 1);   // at the end too
            for (int x = 0; x < len; ++x) q.insert(q.begin() + p, (uint8_t) below(4));
        }
    }
    if (below(3) == 0) {   // a balanced indel pair on top
        const int len = 1 + below(3), p = below((int) q.size() - len), s = below((int) q.size() - len);
        q.erase(q.begin() + p, q.begin() + p + len);
        for (int x = 0; x < len; ++x) q.insert(q.begin() + s, (uint8_t) below(4));
    }
    return q;
}

// Captures bwa_gen_cigar3's -v 4 trace: band(fn) runs fn with bwa_verbose at 4 and stderr
// redirected to one temporary file (under TMPDIR, as the other stderr-capturing tests do),
// and returns the band from its "* Global bandwidth: <w>" line, or -1 when it printed none
// (no DP ran). A guard restores stderr and bwa_verbose even if fn throws.
class BandTrace {
  public:
    BandTrace()
    {
        const char *tmpdir = getenv("TMPDIR");
        if (tmpdir == NULL || tmpdir[0] == '\0') tmpdir = "/tmp";
        std::string tmpl_s = tmpdir;
        if (tmpl_s[tmpl_s.size() - 1] != '/') tmpl_s += '/';
        tmpl_s += "cert_band_trace_XXXXXX";
        std::vector<char> tmpl(tmpl_s.begin(), tmpl_s.end());
        tmpl.push_back('\0');
        fd_ = mkstemp(tmpl.data());
        REQUIRE(fd_ >= 0);
        unlink(tmpl.data());
    }
    ~BandTrace() { if (fd_ >= 0) close(fd_); }

    template <class F>
    int band(F fn)
    {
        REQUIRE(ftruncate(fd_, 0) == 0);
        REQUIRE(lseek(fd_, 0, SEEK_SET) == 0);   // the redirected stderr shares this offset
        fflush(stderr);
        {
            Redirect r(fd_);
            fn();
            fflush(stderr);
        }
        std::string out;
        char buf[4096];
        ssize_t got;
        for (off_t off = 0; (got = pread(fd_, buf, sizeof buf, off)) > 0; off += got) out.append(buf, (size_t) got);
        int band = -1;
        const std::string key = "* Global bandwidth: ";
        const size_t at = out.find(key);
        if (at != std::string::npos) band = atoi(out.c_str() + at + key.size());
        return band;
    }

  private:
    struct Redirect {   // stderr -> fd and bwa_verbose 4 for its lifetime
        int saved_fd, saved_verbose;
        explicit Redirect(int fd) : saved_fd(dup(fileno(stderr))), saved_verbose(bwa_verbose)
        {
            REQUIRE(saved_fd >= 0);
            REQUIRE(dup2(fd, fileno(stderr)) >= 0);
            bwa_verbose = 4;
        }
        ~Redirect()
        {
            bwa_verbose = saved_verbose;
            fflush(stderr);
            const int ok = dup2(saved_fd, fileno(stderr));
            close(saved_fd);
            if (ok < 0) abort();   // stderr would stay redirected for the rest of the binary
        }
    };
    int fd_ = -1;
};

// Largest excursion of a CIGAR's path beyond the corridor between the diagonal and the end
// diagonal (offsets query - reference in [min(0,-dl), max(0,-dl)], dl = tlen - qlen).
int corridor_excursion(const G2 &r, int dl)
{
    const int lo = std::min(0, -dl), hi = std::max(0, -dl);
    int off = 0, mx = 0;
    for (uint32_t c : r.cigar) {
        const int op = c & 0xf, len = (int) (c >> 4);
        if (op == 1) off += len; else if (op == 2) off -= len;
        mx = std::max(mx, std::max(lo - off, off - hi));
    }
    return mx;
}

} // namespace

TEST_CASE("bwa_global_cert_band: unequal lengths, ksw_global2 at the certified band equals the wide band"
          * doctest::test_suite("unit/cigar")) {
    std::mt19937_64 g(0xD17A5Cull);
    auto below = [&](int n) { return (int) (g() % (uint64_t) n); };
    auto unif = [&]() { return (g() >> 11) * (1.0 / 9007199254740992.0); };
    long checked = 0, narrowed[2] = {0, 0}, both_kinds = 0, tight = 0;
    unsigned long wave_narrow = 0, wave16_narrow = 0;
    for (int it = 0; it < 1000; ++it) {
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
        const int m = 30 + below(170);
        int dl = std::min(1 + (below(4) == 0 ? below(30) : below(6)), m - 20);   // query >= 20 bases
        if (it % 2) dl = -dl;   // dl = tlen - qlen: > 0 the reference is longer, < 0 the query
        std::vector<uint8_t> t(m);
        for (auto &c : t) c = unif() < 0.02 ? 4 : (uint8_t) below(4);
        const std::vector<uint8_t> q = edit_unequal(g, t, dl, unif() * 0.12);
        const int n = (int) q.size();
        REQUIRE(m - n == dl);
        const int kc = bwa_global_cert_band(mat, o_del, e_del, o_ins, e_ins, n, q.data(), m, t.data());
        REQUIRE(kc >= std::abs(dl));
        // Independent of the implementation: the best single-gap alignment by brute force over
        // every gap position, and kc from the proof's U0 and slope.
        int amax = mat[0];
        for (int k = 1; k < 25; ++k) amax = std::max(amax, (int) mat[k]);
        const int L = std::min(n, m), ad = std::abs(dl);
        long s_lb = -(1L << 60);
        for (int p = 0; p <= L; ++p) {
            long s = 0;
            for (int k = 0; k < L; ++k)
                s += k < p ? mat[t[k] * 5 + q[k]] : dl > 0 ? mat[t[k + ad] * 5 + q[k]] : mat[t[k] * 5 + q[k + ad]];
            s_lb = std::max(s_lb, s);
        }
        s_lb -= dl > 0 ? o_del + (long) e_del * ad : o_ins + (long) e_ins * ad;
        const long U0 = (long) amax * L - o_ins - o_del - (long) e_ins * std::max(0, -dl) - (long) e_del * std::max(0, dl);
        CHECK(kc == ad + (U0 < s_lb ? 0 : (U0 - s_lb) / (amax + e_ins + e_del)));
        const int floor_w = std::abs(dl) + 3, kw = std::max(kc, floor_w);
        for (int w : {floor_w, floor_w + below(12), floor_w + 8 + below(40), std::max(n, m)}) {
            const G2 wide = run(ksw_global2, q, t, mat, o_del, e_del, o_ins, e_ins, w);
            ++checked;
            CHECK(wide.score >= s_lb);   // the lower bound is a real alignment inside every band
            // the certificate: the result at any band stays within kc - |dl| of the corridor
            CHECK(corridor_excursion(wide, dl) <= kc - std::abs(dl));
            if (corridor_excursion(wide, dl) > 0) ++both_kinds;
            if (kw >= w) continue;
            ++narrowed[dl > 0];
            const unsigned long wv = ksw_g2_wave_exec_count(), wv16 = ksw_g2_wave16_exec_count();
            const G2 cert = run(ksw_global2, q, t, mat, o_del, e_del, o_ins, e_ins, kw);
            wave_narrow += ksw_g2_wave_exec_count() - wv;
            wave16_narrow += ksw_g2_wave16_exec_count() - wv16;
            CHECK(cert == wide);
            CHECK(run(ksw_global2_scalar_ref, q, t, mat, o_del, e_del, o_ins, e_ins, kw) == wide);
            if (kc < floor_w) {   // below the floor (scalar only; production keeps the floor)
                ++tight;
                CHECK(run(ksw_global2_scalar_ref, q, t, mat, o_del, e_del, o_ins, e_ins, kc) == wide);
            }
        }
    }
    MESSAGE("checked=" << checked << " narrowed(ins,del)=" << narrowed[0] << "," << narrowed[1]
            << " off_corridor=" << both_kinds << " below_floor=" << tight
            << " wave_at_narrowed_band=" << wave_narrow << " wave16_at_narrowed_band=" << wave16_narrow);
    for (int s = 0; s < 2; ++s) CHECK(narrowed[s] > checked / 8);   // both gap kinds narrow
    CHECK(both_kinds > 100);   // optima that leave the corridor (an insertion and a deletion)
    CHECK(tight > 100);
    if (ksw_g2_wave_wmin() > 0 || ksw_g2_wave16_wmin() > 0) CHECK(wave_narrow > 0);
    if (ksw_g2_wave16_wmin() > 0) CHECK(wave16_narrow > 0);
}

TEST_CASE("bwa_global_cert_band: unequal lengths, bwa_gen_cigar3 matches ksw_global2 at its uncertified band"
          * doctest::test_suite("unit/cigar")) {
    // As the equal-length integration case, with reference regions longer or shorter than
    // the read: both strands, --meth-like matrices with nm_from_mat, N in the read, and the
    // score-only call.
    std::mt19937_64 g(0xA11C16ull);
    BandTrace trace;
    auto below = [&](int n) { return (int) (g() % (uint64_t) n); };
    const int L = 20000;
    std::vector<uint8_t> ref(L), pac(L / 4 + 64, 0);
    for (int i = 0; i < L; ++i) { ref[i] = (uint8_t) below(4); _set_pac(pac.data(), i, ref[i]); }
    long checked = 0, narrowed[2][2] = {{0, 0}, {0, 0}};   // [strand][reference longer]
    for (int it = 0; it < 1000; ++it) {
        int a = 1, b = 4, o_del = 6, e_del = 1, o_ins = 6, e_ins = 1;
        const int mode = it % 3;
        if (mode > 0) { a = 1 + below(3); b = 1 + below(6); o_del = below(9); e_del = 1 + below(3); o_ins = below(9); e_ins = 1 + below(3); }
        int8_t mat[25];
        bwa_fill_scmat(a, b, mat);
        int nm_from_mat = 0;
        if (mode == 2) {
            const int r = below(4), c = (r + 1 + below(3)) & 3;
            mat[r * 5 + c] = (int8_t) (below(2) ? a : 0);
            nm_from_mat = 1;
        }
        const int m = 30 + below(170), rb = below(L - m - 1);
        int dl = 1 + (below(4) == 0 ? below(20) : below(5));
        if ((it >> 1) % 2) dl = -dl;
        const bool rev = it % 2;
        // the read, as sequenced: edits of the region (forward) or of its reverse complement
        std::vector<uint8_t> tw(ref.begin() + rb, ref.begin() + rb + m);
        std::vector<uint8_t> src = tw;
        if (rev) {
            std::reverse(src.begin(), src.end());
            for (auto &c : src) c = (uint8_t) (3 - c);
        }
        const std::vector<uint8_t> q = edit_unequal(g, src, dl, (g() % 1000) * 0.00008);
        const int n = (int) q.size();
        const int w_ = 1 + below(100);
        // what ksw_global2 sees: on the reverse strand both sequences reversed
        std::vector<uint8_t> qw = q, qb = q;
        if (rev) {
            for (auto &c : tw) c = (uint8_t) (3 - c);
            std::reverse(qw.begin(), qw.end());
        }
        const int w = pre_cert_band(mat, o_del, e_del, o_ins, e_ins, n, m, w_);
        const G2 want = run(ksw_global2, qw, tw, mat, o_del, e_del, o_ins, e_ins, w);
        int want_nm = 0;
        std::string want_md;
        nm_md(want, qw, tw, rev, mat, nm_from_mat, &want_nm, &want_md);

        int score = 0, n_cigar = 0, NM = 0;
        const int64_t qrb = rev ? 2 * (int64_t) L - (rb + m) : rb;
        uint32_t *cig = nullptr;
        const int used = trace.band([&] {
            cig = bwa_gen_cigar3(mat, o_del, e_del, o_ins, e_ins, w_, L, pac.data(), n, qb.data(),
                                 qrb, qrb + m, &score, &n_cigar, &NM, nm_from_mat);
        });
        REQUIRE(cig != nullptr);
        // the band it ran at: the certified one (never below the |dl| + 3 floor) when narrower.
        // (Recomputed with the function under test; got == want below, with want from the
        // wider pre-certificate band, is the independent check.)
        const int kc = bwa_global_cert_band(mat, o_del, e_del, o_ins, e_ins, n, qw.data(), m, tw.data());
        const int kw = std::max(kc, std::abs(dl) + 3);
        CHECK(used == (kc >= 0 && kw < w ? kw : w));
        if (used < w) ++narrowed[rev][dl > 0];
        G2 got;
        got.score = score; got.n_cigar = n_cigar; got.cigar.assign(cig, cig + n_cigar);
        const std::string got_md((const char *) (cig + n_cigar));
        free(cig);
        CHECK(got == want);
        CHECK(NM == want_nm);
        CHECK(got_md == want_md);
        CHECK(qb == q);

        int score_only = 0;
        uint32_t *none = bwa_gen_cigar3(mat, o_del, e_del, o_ins, e_ins, w_, L, pac.data(), n, qb.data(),
                                        qrb, qrb + m, &score_only, nullptr, nullptr, nm_from_mat);
        CHECK(none == nullptr);
        CHECK(score_only == want.score);
        ++checked;
    }
    MESSAGE("checked=" << checked << " narrowed[fwd](ins,del)=" << narrowed[0][0] << "," << narrowed[0][1]
            << " narrowed[rev](ins,del)=" << narrowed[1][0] << "," << narrowed[1][1]);
    for (int s = 0; s < 2; ++s)
        for (int k = 0; k < 2; ++k) CHECK(narrowed[s][k] > 50);   // both strands, both gap kinds
}

TEST_CASE("bwa_global_cert_band: unequal lengths, exact values"
          * doctest::test_suite("unit/cigar")) {
    int8_t mat[25];
    bwa_fill_scmat(1, 4, mat);   // a = 1, b = 4
    const int n = 100;
    std::vector<uint8_t> t(n);
    for (int i = 0; i < n; ++i) t[i] = (uint8_t) (i * 7 % 4);   // shifted by 1-3 it mismatches everywhere
    // mm mismatches away from the indel, then a 2-base deletion from the read (reference
    // longer, dl = 2) or a 2-base insertion into it (dl = -2), at position 50. The best
    // single-gap alignment is the true one, with L - 5*mm column score (L = min length), so
    // U0 - s_lb = 5*mm - o_other: the open of the gap kind the lower bound does not use.
    auto read = [&](int mm, int dl) {
        std::vector<uint8_t> q = t;
        for (int k = 0; k < mm; ++k) q[1 + k * 11] = (uint8_t) ((t[1 + k * 11] + 1) & 3);
        if (dl > 0) q.erase(q.begin() + 50, q.begin() + 52);
        else q.insert(q.begin() + 50, {(uint8_t) ((t[48] + 1) & 3), (uint8_t) ((t[49] + 1) & 3)});
        return q;
    };
    for (int mm = 0; mm <= 6; ++mm) {
        // -O 6 -E 1 both sides: slope 3, kc = 2 + max(0, floor((5*mm - 6) / 3))
        const int sym[] = {2, 2, 3, 5, 6, 8, 10};
        // -O 3,6 (o_del 3, o_ins 6): a deletion's bound subtracts o_ins = 6, an insertion's o_del = 3
        const int ins_asym[] = {2, 2, 4, 6, 7, 9, 11};
        std::vector<uint8_t> qd = read(mm, 2), qi = read(mm, -2);
        CHECK(bwa_global_cert_band(mat, 6, 1, 6, 1, (int) qd.size(), qd.data(), n, t.data()) == sym[mm]);
        CHECK(bwa_global_cert_band(mat, 6, 1, 6, 1, (int) qi.size(), qi.data(), n, t.data()) == sym[mm]);
        CHECK(bwa_global_cert_band(mat, 3, 1, 6, 1, (int) qd.size(), qd.data(), n, t.data()) == sym[mm]);
        CHECK(bwa_global_cert_band(mat, 3, 1, 6, 1, (int) qi.size(), qi.data(), n, t.data()) == ins_asym[mm]);
    }
    // the gap at either end: the lower bound places it there too
    std::vector<uint8_t> head(t.begin() + 3, t.end()), tail(t.begin(), t.end() - 3);
    CHECK(bwa_global_cert_band(mat, 6, 1, 6, 1, (int) head.size(), head.data(), n, t.data()) == 3);
    CHECK(bwa_global_cert_band(mat, 6, 1, 6, 1, (int) tail.size(), tail.data(), n, t.data()) == 3);
    CHECK(bwa_global_cert_band(mat, 6, 1, 6, 1, n, t.data(), (int) head.size(), head.data()) == 3);
    CHECK(bwa_global_cert_band(mat, 6, 1, 6, 1, n, t.data(), (int) tail.size(), tail.data()) == 3);
    // a matrix entry above the match score (2, in a row t never uses, so the columns score as
    // before): amax = 2, slope 4, and U0 - s_lb = 2*98 - (98 - 5*mm) - 6 = 92 + 5*mm
    int8_t hi[25];
    std::copy(mat, mat + 25, hi);
    hi[4 * 5 + 0] = 2;
    const std::vector<uint8_t> qd0 = read(0, 2), qd2 = read(2, 2);
    CHECK(bwa_global_cert_band(hi, 6, 1, 6, 1, (int) qd0.size(), qd0.data(), n, t.data()) == 2 + 23);
    CHECK(bwa_global_cert_band(hi, 6, 1, 6, 1, (int) qd2.size(), qd2.data(), n, t.data()) == 2 + 25);
    // an empty sequence: no certificate
    CHECK(bwa_global_cert_band(mat, 6, 1, 6, 1, 0, t.data(), n, t.data()) == -1);
    CHECK(bwa_global_cert_band(mat, 6, 1, 6, 1, n, t.data(), 0, t.data()) == -1);
}

TEST_CASE("bwa_global_cert_band: unequal lengths, the int32 wavefront kernel at a narrowed band"
          * doctest::test_suite("unit/cigar")) {
    // A large match score puts qlen*a beyond int16, so on tiers with wavefront kernels a band
    // at or above the int32 kernel's minimum width runs it rather than the int16 one. The
    // certified bands here (|dl| = 60 plus a few) clear every tier's minimum.
    std::mt19937_64 g(0x1E32ull);
    auto below = [&](int n) { return (int) (g() % (uint64_t) n); };
    int8_t mat[25];
    bwa_fill_scmat(40, 60, mat);
    const int o_del = 60, e_del = 10, o_ins = 70, e_ins = 10;
    long narrowed = 0;
    unsigned long wave32_narrow = 0;
    for (int it = 0; it < 10; ++it) {
        const int m = 900 + below(200);
        std::vector<uint8_t> t(m);
        for (auto &c : t) c = (uint8_t) below(4);
        const int dl = it % 2 ? -60 : 60;
        std::vector<uint8_t> q = t;
        for (int k = 0, nm = below(5); k < nm; ++k) { const int p = below(m); q[p] = (uint8_t) ((q[p] + 1) & 3); }
        const int p = 100 + below(m - 300);
        if (dl > 0) q.erase(q.begin() + p, q.begin() + p + dl);
        else for (int x = 0; x < -dl; ++x) q.insert(q.begin() + p, (uint8_t) below(4));
        const int n = (int) q.size();
        const int kc = bwa_global_cert_band(mat, o_del, e_del, o_ins, e_ins, n, q.data(), m, t.data());
        REQUIRE(kc >= 60);
        // ksw_g2_wave16_safe declines int16 when qlen*a + a > 32767 - 64 (its SLACK); the
        // int32 kernel's minimum width is at most 52 on every tier, below kw
        REQUIRE((long) n * 40 + 40 > 32767 - 64);
        const int kw = std::max(kc, 63);
        const G2 wide = run(ksw_global2_scalar_ref, q, t, mat, o_del, e_del, o_ins, e_ins, std::max(n, m));
        const unsigned long wv = ksw_g2_wave_exec_count(), wv16 = ksw_g2_wave16_exec_count();
        const G2 cert = run(ksw_global2, q, t, mat, o_del, e_del, o_ins, e_ins, kw);
        wave32_narrow += (ksw_g2_wave_exec_count() - wv) - (ksw_g2_wave16_exec_count() - wv16);
        ++narrowed;
        CHECK(cert == wide);
        CHECK(run(ksw_global2_scalar_ref, q, t, mat, o_del, e_del, o_ins, e_ins, kw) == wide);
        CHECK(run(ksw_global2, q, t, mat, o_del, e_del, o_ins, e_ins, kw + 40) == wide);
    }
    MESSAGE("narrowed=" << narrowed << " wave32_at_narrowed_band=" << wave32_narrow);
    if (ksw_g2_wave_wmin() > 0) CHECK(wave32_narrow == (unsigned long) narrowed);
}
