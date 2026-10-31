// test/unit/test_mem_extend_local_oracle.cpp
//
// The local Smith-Waterman kernels `mem` runs through ksw.h's ksw_align2 model,
// checked against an oracle written here, on every SIMD tier this host can run:
//
//   - ksw_align2 (src/ksw.cpp), the scalar-dispatched striped kernel: the
//     chained-seed SW screen (mem_seed_sw) and the mate-rescue fallback. Both
//     widths: 16-bit, and 8-bit (KSW_XBYTE) with its overflow contract.
//   - the batched kswv kernels (src/kswv.cpp, getScores8 / getScores16), which
//     score the mate-rescue jobs in two passes (pass 0: score, ends and the
//     second-best score; pass 1 on the reversed prefixes: the starts).
//
// The existing tests are differential (kswv against ksw_align2 on every field,
// test_kswv_correctness.cpp), so a fault the two share would pass them. The
// oracle here is a full-matrix DP written from ksw.h's documentation of
// ksw_align and the padding contract kswv.cpp documents. ksw_align2 is held to
// it on hand-worked and random pairs, which with that comparison covers kswv
// too; kswv is held to it directly on the hand-worked pairs, on every tier:
//
//   - local alignment: H = max(0, diagonal + mat[target*5 + query], E, F); a
//     gap of k bases costs o + k*e (o_del/e_del down the target, o_ins/e_ins
//     along the query);
//   - the query is padded to a whole SIMD segment, 16 columns on the 8-bit
//     path and 8 on the 16-bit path, with columns that score 0 against every
//     base (ksw_qinit's profile padding, which kswv reproduces);
//   - score is the best cell; te the first target row reaching it; qe the
//     smallest column of that row holding it;
//   - KSW_XSUBO: every row whose maximum reaches minsc = xtra & 0xffff is a
//     second-best candidate, a row following the candidate's row merging into
//     it only when its maximum is higher (then it moves the candidate there);
//     score2/te2 is the first highest candidate whose row lies outside
//     te +- ceil(score / max(mat));
//   - KSW_XSTOP: stop at the first row whose maximum reaches xtra & 0xffff;
//   - KSW_XBYTE: when the best score reaches 255 - shift (shift = -min(mat)),
//     the pass stops and reports score 255;
//   - KSW_XSTART: tb/qb come from a second pass over the reversed prefixes
//     query[0..qe] and target[0..te] with KSW_XSTOP at the score.
//
// One thing the documentation leaves open: whether a deletion may directly
// follow an insertion. The batched kernels open both gaps from H; the striped
// kernel's lazy-F loop lets that happen in some cells and not in others. The
// oracle therefore runs twice, once with that transition and once without, and
// a field is only compared when both runs agree on it (all fields agree on the
// large majority of pairs; the test asserts that they did).
//
// What is compared is what `mem` reads of a result: score, te, qe, tb and qb
// (the rescued mate's coordinates, bwamem_pair.cpp mem_matesw_batch_post) and
// score2 (its csub, which enters the MAPQ). te2, the second-best's row, is
// computed by the kernels but read nowhere in `mem`, so it is not compared. With
// KSW_XBYTE, ksw.h documents only the reported score (255) after an overflow;
// the other fields of an overflowed pass are not compared either.
//
// Hand-worked cases at the top pin those six fields on inputs with one optimal
// path; random and adversarial pairs (mutated, homopolymer, tandem repeats with
// a real second-best, N bases, all-mismatch, scores at the 8-bit ceiling) are
// scored in batches of different orders, so a lane's result may not depend on
// the pairs sharing its SIMD group.
//
// Tiers: x86 calls the per-tier entry points by name (ksw_align2_avx2 /
// _avx512bw and make_kswv_kernel_avx2 / _avx512bw; libbwa.a carries every
// tier's objects), arm64 the one NEON build. A tier the host lacks is reported
// with MESSAGE and skipped.
//
// Random inputs use std::mt19937 raw output (rng() % n), which the standard
// fixes.

#include "doctest/doctest.h"
#include "scoring.h"

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <memory>
#include <random>
#include <string>
#include <vector>

#include "ksw.h"
#include "kswv.h"
#include "seqpair.h"
#include "seqpair_batch.h"

#if defined(__x86_64__) || defined(__i386__)
extern "C" {
kswr_t ksw_align2_avx2(int, uint8_t *, int, uint8_t *, int, const int8_t *, int, int, int, int, int, kswq_t **);
kswr_t ksw_align2_avx512bw(int, uint8_t *, int, uint8_t *, int, const int8_t *, int, int, int, int, int,
                           kswq_t **);
}
#endif

namespace {

extern "C" {
using Align2Fn = kswr_t (*)(int, uint8_t *, int, uint8_t *, int, const int8_t *, int, int, int, int, int,
                            kswq_t **);
using KswvFactory = Ikswv *(*)(int, int, int, int, int8_t, int8_t, int, int32_t, int32_t);
}

struct Tier {
    const char *name;
    Align2Fn align2;
    KswvFactory kswv;
};

#if defined(__aarch64__)
Ikswv *make_native_kswv(int o_del, int e_del, int o_ins, int e_ins, int8_t a, int8_t b, int nt, int32_t mr,
                        int32_t mq) {
    return new kswv(o_del, e_del, o_ins, e_ins, a, b, nt, mr, mq);
}
#endif

std::vector<Tier> runnable_tiers() {
    std::vector<Tier> ts;
#if defined(__x86_64__) || defined(__i386__)
    __builtin_cpu_init();
    if (__builtin_cpu_supports("avx2")) ts.push_back({"avx2", ksw_align2_avx2, make_kswv_kernel_avx2});
    if (__builtin_cpu_supports("avx512bw"))
        ts.push_back({"avx512bw", ksw_align2_avx512bw, make_kswv_kernel_avx512bw});
    else
        MESSAGE("local SW oracle: no AVX-512BW on this host; the avx512bw tier is skipped");
#elif defined(__aarch64__)
    ts.push_back({"neon", ksw_align2, make_native_kswv});
#endif
    return ts;
}

struct Scoring {
    int a, b, o_del, e_del, o_ins, e_ins;
    int8_t mat[25];
};

// bwa_fill_scmat's matrix written out: a on the diagonal, -b off it, -1 for N.
Scoring scoring(int a, int b, int o_del, int e_del, int o_ins, int e_ins) {
    Scoring s{a, b, o_del, e_del, o_ins, e_ins, {}};
    const bwa_tests::ScoringMatrix m = bwa_tests::build_scoring_matrix(a, b, 1);
    std::copy(m.begin(), m.end(), s.mat);
    return s;
}

int mat_min(const Scoring &s) {
    int m = 0;
    for (int k = 0; k < 25; ++k) m = std::min(m, (int)s.mat[k]);
    return m;
}

int mat_max(const Scoring &s) {
    int m = 0;
    for (int k = 0; k < 25; ++k) m = std::max(m, (int)s.mat[k]);
    return m;
}

std::string show(const kswr_t &r) {
    return "score " + std::to_string(r.score) + " te " + std::to_string(r.te) + " qe " + std::to_string(r.qe) +
           " score2 " + std::to_string(r.score2) + " te2 " + std::to_string(r.te2) + " tb " +
           std::to_string(r.tb) + " qb " + std::to_string(r.qb);
}

// One pass of the model over query q and target t (see the header). With
// `diverged` non-null (and d_after_i set), it is set when some cell's deletion
// score took the insertion-then-deletion transition, i.e. when the run without
// that transition would differ somewhere.
kswr_t local_pass(const std::vector<uint8_t> &q, const std::vector<uint8_t> &t, const Scoring &s, int xtra,
                  bool d_after_i, bool *diverged = nullptr) {
    const bool byte = (xtra & KSW_XBYTE) != 0;
    const int p = byte ? 16 : 8;
    const int Q = (int)q.size(), Qp = (Q + p - 1) / p * p, T = (int)t.size();
    const int shift = -mat_min(s);
    const int minsc = (xtra & KSW_XSUBO) ? xtra & 0xffff : 0x10000;
    const int endsc = (xtra & KSW_XSTOP) ? xtra & 0xffff : 0x10000;
    std::vector<int> prev_v(Qp, 0), cur_v(Qp, 0), vgap_v(Qp, 0), best_row(Qp, 0);
    std::vector<std::pair<int, int> > cand;  // (row maximum, row)
    // Raw pointers and ternaries: the unit binary is built without
    // optimisation, where std::max and operator[] are calls.
    int *prev = prev_v.data(), *cur = cur_v.data(), *vgap = vgap_v.data();
    const uint8_t *qd = q.data();
    const int oe_del = s.o_del + s.e_del, oe_ins = s.o_ins + s.e_ins, e_del = s.e_del, e_ins = s.e_ins;
    int gmax = 0, te = -1;
    for (int i = 0; i < T; ++i) {
        int f = 0, rmax = 0, hd = 0;
        const int8_t *row = s.mat + t[(size_t)i] * 5;
        for (int j = 0; j < Qp; ++j) {
            const int sc = j < Q ? row[qd[j]] : 0;
            int m = hd + sc;
            if (m < 0) m = 0;
            hd = prev[j];
            const int e = vgap[j];
            int h = m > e ? m : e;
            if (f > h) h = f;
            cur[j] = h;
            if (h > rmax) rmax = h;
            const int me = m > e ? m : e;
            int e_next = e - e_del, t_open = (d_after_i ? h : me) - oe_del;
            if (t_open > e_next) e_next = t_open;
            if (e_next < 0) e_next = 0;
            if (diverged != nullptr) {
                int e_without = e - e_del;
                if (me - oe_del > e_without) e_without = me - oe_del;
                if (e_without < 0) e_without = 0;
                if (e_without != e_next) *diverged = true;
            }
            vgap[j] = e_next;
            int f_next = f - e_ins;
            if (h - oe_ins > f_next) f_next = h - oe_ins;
            f = f_next < 0 ? 0 : f_next;
        }
        if (rmax >= minsc) {
            if (cand.empty() || cand.back().second + 1 != i) cand.push_back(std::make_pair(rmax, i));
            else if (cand.back().first < rmax) cand.back() = std::make_pair(rmax, i);
        }
        bool stop = false;
        if (rmax > gmax) {
            gmax = rmax, te = i;
            std::copy(cur, cur + Qp, best_row.begin());
            stop = (byte && gmax + shift >= 255) || gmax >= endsc;
        }
        std::swap(prev, cur);
        if (stop) break;
    }
    kswr_t r = g_defr;
    r.score = byte && gmax + shift >= 255 ? 255 : gmax;
    r.te = te;
    if (r.score != 255 || !byte) {
        int mx = -1;
        for (int j = 0; j < Qp; ++j)
            if (best_row[(size_t)j] > mx) mx = best_row[(size_t)j], r.qe = j;
        const int c = (r.score + mat_max(s) - 1) / mat_max(s);
        for (const auto &e : cand)
            if ((e.second < te - c || e.second > te + c) && e.first > r.score2) r.score2 = e.first, r.te2 = e.second;
    }
    return r;
}

// ksw_align2's two passes under the model.
kswr_t local_align(const std::vector<uint8_t> &q, const std::vector<uint8_t> &t, const Scoring &s, int xtra,
                   bool d_after_i, bool *diverged = nullptr) {
    kswr_t r = local_pass(q, t, s, xtra, d_after_i, diverged);
    if ((xtra & KSW_XSTART) == 0 || ((xtra & KSW_XSUBO) && r.score < (xtra & 0xffff))) return r;
    std::vector<uint8_t> rq(q.begin(), q.begin() + (r.qe + 1)), rt = t;
    std::reverse(rq.begin(), rq.end());
    std::reverse(rt.begin(), rt.begin() + (r.te + 1));
    const kswr_t rr = local_pass(rq, rt, s, (xtra & KSW_XBYTE) | KSW_XSTOP | r.score, d_after_i, diverged);
    if (r.score == rr.score) r.tb = r.te - rr.te, r.qb = r.qe - rr.qe;
    return r;
}

// The expected result: the model's answer where the two gap conventions agree,
// field by field. `known` marks the fields both agree on.
struct Expected {
    kswr_t r;
    bool known[7];
    bool all() const { return std::all_of(known, known + 7, [](bool k) { return k; }); }
};

Expected expected(const std::vector<uint8_t> &q, const std::vector<uint8_t> &t, const Scoring &s, int xtra) {
    bool diverged = false;
    const kswr_t u = local_align(q, t, s, xtra, true, &diverged);
    // Without a divergent cell the two runs are the same DP.
    const kswr_t l = diverged ? local_align(q, t, s, xtra, false) : u;
    Expected e;
    e.r = u;
    e.known[0] = u.score == l.score;
    e.known[1] = e.known[0] && u.te == l.te;
    e.known[2] = e.known[1] && u.qe == l.qe;
    e.known[3] = e.known[1] && u.score2 == l.score2;
    e.known[4] = e.known[3] && u.te2 == l.te2;
    e.known[5] = e.known[2] && u.tb == l.tb;
    e.known[6] = e.known[2] && u.qb == l.qb;
    return e;
}

// Compare the fields the model pins and `mem` reads (te2 is not read; see the
// header).
void check_fields(const kswr_t &got, const Expected &want) {
    INFO("got " << show(got) << "; want " << show(want.r));
    if (want.known[0]) CHECK(got.score == want.r.score);
    if (want.known[1]) CHECK(got.te == want.r.te);
    if (want.known[2]) CHECK(got.qe == want.r.qe);
    if (want.known[3]) CHECK(got.score2 == want.r.score2);
    if (want.known[5]) CHECK(got.tb == want.r.tb);
    if (want.known[6]) CHECK(got.qb == want.r.qb);
}

kswr_t run_align2(Align2Fn fn, std::vector<uint8_t> q, std::vector<uint8_t> t, const Scoring &s, int xtra) {
    return fn((int)q.size(), q.data(), (int)t.size(), t.data(), 5, s.mat, s.o_del, s.e_del, s.o_ins, s.e_ins,
              xtra, nullptr);
}

// The two-pass batched protocol `mem` runs for mate rescue (bwamem_pair.cpp
// mem_sam_pe_batch_run), through one tier's kswv object. Results in input order.
std::vector<kswr_t> run_kswv(KswvFactory make, const std::vector<bwa_tests::TestPair> &pairs, const Scoring &s,
                             int xtra, bool use16) {
    int32_t max_ref = 1, max_qer = 1;
    for (const auto &p : pairs) {
        max_ref = std::max(max_ref, (int32_t)p.ref.size());
        max_qer = std::max(max_qer, (int32_t)p.qry.size());
    }
    std::unique_ptr<Ikswv> k(make(s.o_del, s.e_del, s.o_ins, s.e_ins, (int8_t)s.a, (int8_t)-s.b, 1, max_ref, max_qer));
    bwa_tests::BatchBuffers bb(pairs, xtra);
    if (use16) k->getScores16(bb.pairs(), bb.ref_buf(), bb.qer_buf(), bb.aln(), bb.n(), 1, 0);
    else k->getScores8(bb.pairs(), bb.ref_buf(), bb.qer_buf(), bb.aln(), bb.n(), 1, 0);
    const int pos = bb.prepare_phase1();
    if (use16) k->getScores16(bb.pairs(), bb.ref_buf(), bb.qer_buf(), bb.aln(), pos, 1, 1);
    else k->getScores8(bb.pairs(), bb.ref_buf(), bb.qer_buf(), bb.aln(), pos, 1, 1);
    return std::vector<kswr_t>(bb.aln(), bb.aln() + bb.n());
}

std::vector<uint8_t> random_seq(std::mt19937 &rng, size_t n) {
    std::vector<uint8_t> v(n);
    for (size_t k = 0; k < n; ++k) v[k] = (uint8_t)(rng() % 4);
    return v;
}

std::vector<uint8_t> digits(const char *s) {
    std::vector<uint8_t> v;
    for (; *s; ++s) v.push_back((uint8_t)(*s - '0'));
    return v;
}

// A rescue-shaped pair: the query planted, mutated, in a longer random target,
// or one of the adversarial shapes.
bwa_tests::TestPair make_pair(std::mt19937 &rng, int kind, int max_q) {
    bwa_tests::TestPair p;
    const size_t ql = 10 + rng() % (uint32_t)(max_q - 10);
    std::vector<uint8_t> q = random_seq(rng, ql), body;
    switch (kind % 6) {
        case 0:  // mutated copy: substitutions and short indels
            for (size_t k = 0; k < q.size(); ++k) {
                const uint32_t r = rng() % 100;
                if (r < 6) body.push_back((uint8_t)((q[k] + 1 + rng() % 3) % 4));
                else if (r < 9) continue;
                else if (r < 12) {
                    for (uint32_t m = 1 + rng() % 4; m > 0; --m) body.push_back((uint8_t)(rng() % 4));
                    body.push_back(q[k]);
                } else body.push_back(q[k]);
            }
            break;
        case 1: {  // homopolymer runs
            const uint8_t base = (uint8_t)(rng() % 4);
            for (auto &c : q) c = rng() % 8 == 0 ? (uint8_t)(rng() % 4) : base;
            body = q;
            break;
        }
        case 2: {  // the query twice, the second copy diverged: a real second-best
            body = q;
            const std::vector<uint8_t> gap = random_seq(rng, 5 + rng() % 40);
            body.insert(body.end(), gap.begin(), gap.end());
            for (uint8_t c : q) body.push_back(rng() % 10 == 0 ? (uint8_t)((c + 1) % 4) : c);
            break;
        }
        case 3:  // N bases on both sides
            for (auto &c : q)
                if (rng() % 12 == 0) c = 4;
            body = q;
            for (auto &c : body)
                if (rng() % 12 == 0) c = rng() % 2 ? (uint8_t)4 : (uint8_t)((c + 2) % 4);
            break;
        case 4:  // every base differs: only short chance alignments
            for (uint8_t c : q) body.push_back((uint8_t)((c + 1 + rng() % 3) % 4));
            break;
        default:  // unrelated
            body = random_seq(rng, ql + rng() % 30);
            break;
    }
    const std::vector<uint8_t> left = random_seq(rng, rng() % 60), right = random_seq(rng, rng() % 60);
    p.ref = left;
    p.ref.insert(p.ref.end(), body.begin(), body.end());
    p.ref.insert(p.ref.end(), right.begin(), right.end());
    if (p.ref.empty()) p.ref.push_back(0);
    p.qry = q;
    return p;
}

}  // namespace

TEST_CASE("local SW oracle: hand-worked ksw_align2 and kswv results on every runnable tier"
          * doctest::test_suite("unit/kswv")) {
    const std::vector<Tier> tiers = runnable_tiers();
    REQUIRE(!tiers.empty());
    const Scoring def = scoring(1, 4, 6, 1, 6, 1);
    // 40 bases, no base equal to its neighbour.
    const std::vector<uint8_t> q = digits("0123102302130120312013201230210310231032");
    const std::vector<uint8_t> flank = digits("3333333333");  // T x 10: q starts and ends in A/C/G
    const int xtra16 = KSW_XSUBO | KSW_XSTART | 19, xtra8 = xtra16 | KSW_XBYTE;

    struct Case {
        const char *what;
        std::vector<uint8_t> qry, ref;
        kswr_t want;
    };
    std::vector<Case> cases;
    {
        // The query planted after 10 flank bases: score 40 at target 10..49.
        std::vector<uint8_t> ref = flank;
        ref.insert(ref.end(), q.begin(), q.end());
        ref.insert(ref.end(), flank.begin(), flank.end());
        cases.push_back({"exact copy at target 10", q, ref, {40, 49, 39, -1, -1, 10, 0}});
        // One substitution at query base 20: 39 - 4 = 35, still end to end
        // (each half alone scores at most 20).
        std::vector<uint8_t> sub = ref;
        sub[30] = (uint8_t)((sub[30] + 2) % 4);
        cases.push_back({"one substitution", q, sub, {35, 49, 39, -1, -1, 10, 0}});
        // Three target bases inserted after query base 19: 40 - (6 + 3) = 31.
        std::vector<uint8_t> del = ref;
        del.insert(del.begin() + 30, {2, 2, 2});
        cases.push_back({"a 3-base deletion", q, del, {31, 52, 39, -1, -1, 10, 0}});
        // The query copied twice, 20 flank bases apart: both copies score 40;
        // the second is the second-best candidate, 60 rows on, outside te +- 40.
        std::vector<uint8_t> two = ref;
        two.insert(two.end(), flank.begin(), flank.end());
        two.insert(two.end(), q.begin(), q.end());
        cases.push_back({"two copies: score2 is the second", q, two, {40, 49, 39, 40, 109, 10, 0}});
    }
    for (const Case &c : cases) {
        CAPTURE(c.what);
        // The oracle agrees with the hand-worked values under both conventions.
        for (bool dai : {true, false}) {
            CAPTURE(dai);
            const kswr_t o = local_align(c.qry, c.ref, def, xtra16, dai);
            INFO("oracle " << show(o) << "; want " << show(c.want));
            CHECK(memcmp(&o, &c.want, sizeof o) == 0);
        }
        Expected want;
        want.r = c.want;
        std::fill(want.known, want.known + 7, true);
        for (const Tier &tier : tiers) {
            CAPTURE(tier.name);
            check_fields(run_align2(tier.align2, c.qry, c.ref, def, xtra16), want);
            check_fields(run_align2(tier.align2, c.qry, c.ref, def, xtra8), want);
            const std::vector<bwa_tests::TestPair> one = {{c.ref, c.qry, ""}};
            check_fields(run_kswv(tier.kswv, one, def, xtra8, false)[0], want);
            check_fields(run_kswv(tier.kswv, one, def, xtra16, true)[0], want);
        }
    }
}

TEST_CASE("local SW oracle: the second-best score's exclusion window, at both edges"
          * doctest::test_suite("unit/kswv")) {
    // ksw.h's KSW_XSUBO second-best ignores candidates within ceil(score /
    // max(mat)) target rows of the best's end te. A 48-base query (a whole
    // segment at both widths, so no padding columns) is placed twice: once
    // exact (48) and once with base 5 changed (43, a candidate whose run of
    // rows at or above minsc 19 starts after the change, so the candidate sits
    // at its copy's last row). With the copies adjacent, the two ends are 48
    // rows apart: exactly ceil(48 / 1), the window's edge, so score2 is -1;
    // with one base between them they are 49 apart and score2 is 43 at the
    // weaker copy's end. Before the best (low edge) the target continues
    // after the best; after it (high edge) the target stops at the weaker
    // copy, so no decaying rows follow it.
    const std::vector<Tier> tiers = runnable_tiers();
    REQUIRE(!tiers.empty());
    const Scoring def = scoring(1, 4, 6, 1, 6, 1);
    std::mt19937 rng(4848u);
    const std::vector<uint8_t> q = random_seq(rng, 48), flank = random_seq(rng, 10);
    std::vector<uint8_t> weak = q;
    weak[5] = (uint8_t)((weak[5] + 1) % 4);
    struct Case {
        const char *what;
        std::vector<uint8_t> ref;
        int score2, te2;
    };
    std::vector<Case> cases;
    {
        std::vector<uint8_t> t = flank;  // weak copy ends at 57, the best at 105
        t.insert(t.end(), weak.begin(), weak.end());
        t.insert(t.end(), q.begin(), q.end());
        t.insert(t.end(), flank.begin(), flank.end());
        cases.push_back({"weaker copy ends 48 rows before the best: excluded", t, -1, -1});
        std::vector<uint8_t> t1 = flank;  // weak copy ends at 57, the best at 106
        t1.insert(t1.end(), weak.begin(), weak.end());
        t1.push_back(0);
        t1.insert(t1.end(), q.begin(), q.end());
        t1.insert(t1.end(), flank.begin(), flank.end());
        cases.push_back({"weaker copy ends 49 rows before the best: counted", t1, 43, 57});
        std::vector<uint8_t> t2 = flank;  // the best ends at 57, the weak copy at 105
        t2.insert(t2.end(), q.begin(), q.end());
        t2.insert(t2.end(), weak.begin(), weak.end());
        cases.push_back({"weaker copy ends 48 rows after the best: excluded", t2, -1, -1});
        std::vector<uint8_t> t3 = flank;  // the best ends at 57, the weak copy at 106
        t3.insert(t3.end(), q.begin(), q.end());
        t3.push_back(0);
        t3.insert(t3.end(), weak.begin(), weak.end());
        cases.push_back({"weaker copy ends 49 rows after the best: counted", t3, 43, 106});
    }
    const int xtra16 = KSW_XSUBO | KSW_XSTART | 19, xtra8 = xtra16 | KSW_XBYTE;
    for (const Case &c : cases) {
        CAPTURE(c.what);
        for (int xtra : {xtra16, xtra8}) {
            CAPTURE(xtra);
            const Expected want = expected(q, c.ref, def, xtra);
            REQUIRE(want.all());
            CHECK(want.r.score == 48);
            CHECK(want.r.score2 == c.score2);
            CHECK(want.r.te2 == c.te2);
            for (const Tier &tier : tiers) {
                CAPTURE(tier.name);
                check_fields(run_align2(tier.align2, q, c.ref, def, xtra), want);
                const std::vector<bwa_tests::TestPair> one = {{c.ref, q, ""}};
                check_fields(run_kswv(tier.kswv, one, def, xtra, xtra == xtra16)[0], want);
            }
        }
    }
}

TEST_CASE("local SW oracle: the 8-bit pass reports 255 once the score reaches 255 - shift"
          * doctest::test_suite("unit/kswv")) {
    // ksw.h: "KSW_XBYTE: use an unsigned byte to store a score. If overflow
    // occurs, kswr_t::score will be set to 255". With -A 1 -B 4, shift is 4,
    // so a 250-base exact match (250 + 4 = 254) is still exact and a 251-base
    // one (251 + 4 = 255) is reported as 255; the 16-bit pass scores both. Only
    // the score is documented after an overflow, so only the score is checked
    // there.
    const std::vector<Tier> tiers = runnable_tiers();
    REQUIRE(!tiers.empty());
    const Scoring def = scoring(1, 4, 6, 1, 6, 1);
    std::mt19937 rng(255u);
    const std::vector<uint8_t> base = random_seq(rng, 300);
    for (int len : {250, 251}) {
        CAPTURE(len);
        const std::vector<uint8_t> qry(base.begin() + 20, base.begin() + 20 + len);
        for (const Tier &tier : tiers) {
            CAPTURE(tier.name);
            const kswr_t r8 = run_align2(tier.align2, qry, base, def, KSW_XBYTE);
            const kswr_t r16 = run_align2(tier.align2, qry, base, def, 0);
            CHECK(r16.score == len);
            CHECK(r16.te == 20 + len - 1);
            CHECK(r16.qe == len - 1);
            if (len == 250) {
                // No overflow: every field as the model gives it (no second-best
                // or start without KSW_XSUBO / KSW_XSTART).
                CHECK(r8.score == 250);
                CHECK(r8.te == 20 + len - 1);
                CHECK(r8.qe == len - 1);
                const kswr_t o8 = local_pass(qry, base, def, KSW_XBYTE, true);
                CHECK(memcmp(&o8, &r8, sizeof o8) == 0);
            } else {
                CHECK(r8.score == 255);
            }
        }
    }
}

TEST_CASE("local SW oracle: ksw_align2 random and adversarial pairs equal the model"
          * doctest::test_suite("unit/kswv")) {
    const std::vector<Tier> tiers = runnable_tiers();
    REQUIRE(!tiers.empty());
    const Scoring scorings[] = {
        scoring(1, 4, 6, 1, 6, 1),    // defaults
        scoring(1, 4, 5, 2, 7, 1),    // -O 5,7 -E 2,1
        scoring(2, 8, 12, 2, 12, 2),  // -A 2
    };
    std::mt19937 rng(31337u);
    int compared = 0, full = 0;
    for (const Scoring &s : scorings) {
        const int minsc = 19 * s.a;
        for (int n = 0; n < 36; ++n) {
            const bwa_tests::TestPair p = make_pair(rng, n, 100);
            // mem's two uses: the rescue job (XSUBO | XSTART | minsc, 8-bit
            // when the query is short enough) and the seed screen (XSTART).
            const bool u8 = (int)p.qry.size() * s.a + -mat_min(s) <= 254;
            const int xtras[] = {KSW_XSUBO | KSW_XSTART | minsc | (u8 ? KSW_XBYTE : 0), KSW_XSTART};
            for (int xtra : xtras) {
                CAPTURE(n);
                CAPTURE(xtra);
                const Expected want = expected(p.qry, p.ref, s, xtra);
                if (want.all()) ++full;
                for (const Tier &tier : tiers) {
                    CAPTURE(tier.name);
                    check_fields(run_align2(tier.align2, p.qry, p.ref, s, xtra), want);
                    ++compared;
                }
            }
        }
    }
    CHECK(compared >= 2 * 3 * 36 * (int)tiers.size());
    // The two gap conventions agreed on every field for most pairs.
    CHECK(full >= 3 * 36 * 2 * 9 / 10);
}
