/* Equivalence check for the SIMD rescue-pruning filters (rescue_prune_neon.h on aarch64,
 * rescue_prune_x86.h on AVX2 hosts) against the scalar filter, through the public
 * rescue_prune_window() / rescue_prune_last_view() API.
 *
 * For every job: the SIMD path must return the scalar path's (kind, hb, he), and on B2 its view
 * must match the scalar arrays: cnt and minrow (where cnt > 0) on [0, nd), bnd16 ==
 * base + fwd + bwd - (a cnt - c), the mw / hw bits (bnd >= minsc, and cnt > 0) on [0, nd) and zero
 * on [nd, nd + 64), and the component list (runs of bnd >= minsc with a hit: a, b, ub, i0, dmax)
 * with its count, and the view's minsc and constants. Every job is run twice in a row half of the
 * time, which exercises the memo. minsc varies per job (dumps: 19, 5, 10, 15, 25, 32 in turn;
 * fuzz: 5-40), and so does the scoring (the default half of the time, else one of SCORINGS: c = 2,
 * a = 2, large gap opens), through rescue_prune_params. At c = 2 a single-hit diagonal has
 * weight a - c < 0, the case the NEON hw test once got wrong. A share of the jobs runs under the
 * --meth relation (set_meth_rel, relx T or A), some of them three times in a row with the other
 * hypothesis in between (the query caches and the memo must key on relx), and the fuzz adds long
 * T- or A-rich mates that exceed the relation's entry cap (FULL on both sides).
 *
 *   rescue_prune_eq dump <files...>     jobs from rescue_band_harness dumps
 *                                       (int32 len1, len2; ref; query)
 *   rescue_prune_eq dump16 <files...>   jobs from the filter-bench dumps
 *                                       (int16 len1, len2, score, te, hb, he; ref; query)
 *   rescue_prune_eq fuzz <n> <seed>     n generated jobs (repeats, low complexity, N, tiny and
 *                                       long windows, gate values)
 *   rescue_prune_eq dumpm <files...>    jobs from the --meth rescue dumps (int32 len1, len2, hyp,
 *                                       xbyte; ref; query): the default scoring under the relation of
 *                                       each job's hypothesis (hyp < 0: none), a quarter of them also
 *                                       under the other hypothesis
 *   rescue_prune_eq timem <mh> <files...>  no checks: the filter's single-thread time over the same
 *                                       jobs at gate mh, relation vs converted copies vs the scalar
 *                                       relation, with the rows and cells (rows x len2) each keeps
 * Exit status 1 on any disagreement. Build: make arch=avx2 rescue-prune-eq (x86), make arch=arm64
 * rescue-prune-eq (aarch64). */
#include "../src/rescue_prune.h"

#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <random>
#include <string>
#include <vector>

namespace {

struct Stats {
    long jobs = 0, simd = 0, fallback = 0, b2 = 0, bad = 0, kind[3] = {0, 0, 0}, rel = 0, rel_b2 = 0;
};

/* Scorings (a, b, o_del, e_del, o_ins, e_ins, k_max): K = 5 ones, and K > 5 ones, which both SIMD
 * filters take with exact matching (the scalar reference is reset in full now and then, which checks
 * the tables' touched-code reset). */
const int SCORINGS[][7] = {
    {1, 4, 8, 2, 8, 2, 5},     // -O 8 -E 2: c = 2
    {1, 4, 6, 2, 6, 2, 5},     // -E 2: c = 2
    {1, 6, 6, 1, 6, 1, 5},     // -B 6
    {1, 9, 16, 1, 16, 1, 5},   // -x intractg: tail offset 12
    {2, 8, 12, 2, 12, 2, 5},   // -A 2 scaled: a = 2, c = 2
    {1, 4, 6, 1, 7, 2, 5},     // -O 6,7 -E 1,2
    {3, 12, 18, 3, 18, 3, 5},  // a = 3
    {1, 6, 6, 1, 6, 1, 7},     // -B 6 at K = 7
    {1, 9, 16, 1, 16, 1, 8},   // -x intractg at K = 8
    {1, 8, 6, 1, 6, 1, 6},     // -B 8 at K = 6
};
const int NSCORINGS = sizeof SCORINGS / sizeof SCORINGS[0];

rescue_prune_params scoring(int which, int minsc)
{
    if (which < 0) return rescue_prune_params::defaults(minsc);
    const int *x = SCORINGS[which];
    return rescue_prune_params::from(x[0], x[1], x[2], x[3], x[4], x[5], minsc * x[0], 5, x[6]);
}

bool check_view(const uint8_t *ref, int len1, const uint8_t *q, int len2, const rescue_prune_params &p,
                const rescue_prune_view &v, rescue_prune_scratch &r, std::string &why)
{
    const int nd = r.view_nd, minsc = p.minsc;
    if (v.minsc != minsc) { why = "view minsc " + std::to_string(v.minsc); return false; }
    if (v.base != p.base() || v.a != p.a || v.c != p.c) { why = "view constants"; return false; }
    if (v.nd != nd) { why = "nd " + std::to_string(v.nd) + " vs " + std::to_string(nd); return false; }
    if (v.off != ((len2 + 15) / 16) * 16) { why = "off"; return false; }
    if (!v.bnd16 || !v.mw || !v.hw || !v.comps || v.ncomp < 0) { why = "SIMD view incomplete"; return false; }
    auto bit = [](const uint64_t *w, int x) { return (int)(w[x >> 6] >> (x & 63) & 1); };
    for (int x = 0; x < nd; x++) {
        const int bnd = p.base() + r.fwd[x] + r.bwd[x] - p.weight(r.cnt[x]);
        if (v.cnt[x] != r.cnt[x]) { why = "cnt at " + std::to_string(x); return false; }
        if (r.cnt[x] && v.minrow[x] != r.minrow[x]) { why = "minrow at " + std::to_string(x); return false; }
        if (v.bnd16[x] != bnd) { why = "bnd at " + std::to_string(x); return false; }
        if (bit(v.mw, x) != (bnd >= minsc)) { why = "mw at " + std::to_string(x); return false; }
        if (bit(v.hw, x) != (bnd >= minsc && r.cnt[x] > 0)) { why = "hw at " + std::to_string(x); return false; }
    }
    for (int x = nd; x < nd + 64; x++)
        if (bit(v.mw, x) || bit(v.hw, x)) { why = "bitset pad at " + std::to_string(x); return false; }
    /* components at minsc with a hit, as the filters' step 7 lists them */
    std::vector<rescue_prune_neon::Comp> want;
    for (int x = 0; x < nd;) {
        auto bnd = [&](int y) { return p.base() + r.fwd[y] + r.bwd[y] - p.weight(r.cnt[y]); };
        if (bnd(x) < minsc) { x++; continue; }
        const int a = x;
        int ub = INT16_MIN, i0 = INT16_MAX, dmax = -1;
        for (; x < nd && bnd(x) >= minsc; x++) {
            ub = std::max(ub, bnd(x));
            if (r.cnt[x]) { i0 = std::min(i0, (int)r.minrow[x]); dmax = x; }
        }
        if (dmax >= 0) want.push_back(rescue_prune_neon::Comp{a, x, ub, i0, dmax});
    }
    if (v.ncomp != (int)want.size()) { why = "ncomp " + std::to_string(v.ncomp) + " vs " + std::to_string(want.size()); return false; }
    for (int k = 0; k < v.ncomp_stored; k++) {
        const rescue_prune_neon::Comp &a = v.comps[k], &b = want[k];
        if (a.a != b.a || a.b != b.b || a.ub != b.ub || a.i0 != b.i0 || a.dmax != b.dmax) {
            why = "comp " + std::to_string(k);
            return false;
        }
    }
    (void)ref; (void)len1; (void)q;
    return true;
}

void run_one(const uint8_t *ref, int len1, const uint8_t *q, int len2, const rescue_prune_params &p, int max_hits,
             bool twice, Stats &st)
{
    static rescue_prune_scratch r;
    if (!p.valid) return;
    const int minsc = p.minsc;
    int h0 = -1, e0 = -1, h1 = -1, e1 = -1;
    static long nk = 0;
    if (p.K != 5 && nk++ % 64 == 0) { r.ntouched = -1; r.qlen_c = -1; }   // full reset of the reference
    const int k0 = rescue_prune_window_scalar(ref, len1, q, len2, p, max_hits, r, &h0, &e0);
    for (int rep = 0; rep < (twice ? 2 : 1); rep++) {
        const int k1 = rescue_prune_window(ref, len1, q, len2, p, max_hits, &h1, &e1);
        const int path = rescue_prune_last_path();
        st.jobs++;
        if (path == 2 || path == 3) st.simd++;
        else if (path == 1) st.fallback++;
        /* the wrapper returns FULL for len1 < 5 / len2 < 5 / oversize before the filters run */
        if (path == 0) continue;
        std::string why;
        bool ok = k1 == k0 && (k1 != RESCUE_PRUNE_B2 || (h1 == h0 && e1 == e0));
        if (!ok) why = "decision (" + std::to_string(k1) + "," + std::to_string(h1) + "," + std::to_string(e1) +
                       ") vs scalar (" + std::to_string(k0) + "," + std::to_string(h0) + "," + std::to_string(e0) + ")";
        if (ok && k1 == RESCUE_PRUNE_B2 && (path == 2 || path == 3))
            ok = check_view(ref, len1, q, len2, p, rescue_prune_last_view(), r, why);
        if (rep == 0) st.kind[k1]++;
        if (k1 == RESCUE_PRUNE_B2) st.b2++;
        if (p.relx >= 0) { st.rel++; st.rel_b2 += k1 == RESCUE_PRUNE_B2; }
        if (!ok) {
            if (st.bad < 20)
                fprintf(stderr, "MISMATCH len1=%d len2=%d minsc=%d a=%d c=%d base=%d relx=%d max_hits=%d rep=%d path=%d: %s\n",
                        len1, len2, minsc, p.a, p.c, p.base(), p.relx, max_hits, rep, path, why.c_str());
            st.bad++;
        }
    }
}

/* run_one under the relation of hypothesis hyp (relx T for OT, A for OB); alias: then under the other
 * hypothesis and this one again, on the identical bytes. */
void run_rel(const uint8_t *ref, int len1, const uint8_t *q, int len2, rescue_prune_params p, int hyp,
             int max_hits, bool twice, bool alias, Stats &st)
{
    p.set_meth_rel(hyp);
    run_one(ref, len1, q, len2, p, max_hits, twice, st);
    if (!alias) return;
    rescue_prune_params o = p;
    o.set_meth_rel(!hyp);
    run_one(ref, len1, q, len2, o, max_hits, false, st);
    run_one(ref, len1, q, len2, p, max_hits, false, st);
}

int run_dumps(int argc, char **argv, bool hdr16, Stats &st)
{
    for (int a = 0; a < argc; a++) {
        FILE *fp = fopen(argv[a], "rb");
        if (!fp) { perror(argv[a]); return 1; }
        std::vector<uint8_t> buf;
        fseek(fp, 0, SEEK_END);
        const long n = ftell(fp);
        fseek(fp, 0, SEEK_SET);
        buf.resize((size_t)n);
        if (fread(buf.data(), 1, (size_t)n, fp) != (size_t)n) { perror("read"); fclose(fp); return 1; }
        fclose(fp);
        long k = 0;
        const size_t hs = hdr16 ? 12 : 8;
        for (size_t p = 0; p + hs <= buf.size(); k++) {
            int h[2];
            if (hdr16) {
                int16_t h16[6];
                memcpy(h16, buf.data() + p, 12);
                h[0] = h16[0]; h[1] = h16[1];
            } else {
                int32_t h32[2];
                memcpy(h32, buf.data() + p, 8);
                h[0] = h32[0]; h[1] = h32[1];
            }
            p += hs;
            if (h[0] < 0 || h[1] < 0 || p + h[0] + h[1] > buf.size()) break;
            const uint8_t *ref = buf.data() + p, *q = ref + h[0];
            p += (size_t)h[0] + h[1];
            /* production gates: 400 (hull path) and 1000 (banded), plus fully open */
            const int mh = k % 3 == 0 ? 400 : k % 3 == 1 ? 1000 : 1 << 30;
            static const int minscs[6] = {19, 5, 10, 15, 25, 32};
            const int which = (k / 18) % 2 == 0 ? -1 : (int)((k / 36) % NSCORINGS);
            const rescue_prune_params sp = scoring(which, minscs[(k / 3) % 6]);
            /* every 4th group of 5 jobs under the relation, hypotheses alternating; each 2nd of them
             * also under the other hypothesis */
            if ((k / 5) % 4 == 3) run_rel(ref, h[0], q, h[1], sp, (int)((k / 20) & 1), mh, (k & 1) == 0, (k / 5) % 8 == 7, st);
            else run_one(ref, h[0], q, h[1], sp, mh, (k & 1) == 0, st);
        }
    }
    return 0;
}

/* --meth dumps: records of int32 len1, len2, hyp, xbyte, then the window and the oriented mate. */
struct MJob { const uint8_t *ref, *q; int len1, len2, hyp; };
bool load_mdumps(int argc, char **argv, std::vector<std::vector<uint8_t>> &bufs, std::vector<MJob> &jobs)
{
    for (int a = 0; a < argc; a++) {
        FILE *fp = fopen(argv[a], "rb");
        if (!fp) { perror(argv[a]); return false; }
        fseek(fp, 0, SEEK_END);
        const long n = ftell(fp);
        fseek(fp, 0, SEEK_SET);
        bufs.emplace_back((size_t)n);
        std::vector<uint8_t> &buf = bufs.back();
        if (fread(buf.data(), 1, (size_t)n, fp) != (size_t)n) { perror("read"); fclose(fp); return false; }
        fclose(fp);
        for (size_t p = 0; p + 16 <= buf.size();) {
            int32_t h[4];
            memcpy(h, buf.data() + p, 16);
            p += 16;
            if (h[0] < 0 || h[1] < 0 || p + h[0] + h[1] > buf.size()) break;
            jobs.push_back(MJob{buf.data() + p, buf.data() + p + h[0], h[0], h[1], h[2]});
            p += (size_t)h[0] + h[1];
        }
    }
    return true;
}

int run_mdumps(int argc, char **argv, Stats &st)
{
    std::vector<std::vector<uint8_t>> bufs;
    std::vector<MJob> jobs;
    if (!load_mdumps(argc, argv, bufs, jobs)) return 1;
    long k = 0;
    for (const MJob &J : jobs) {
        const int mh = k % 3 == 0 ? 400 : k % 3 == 1 ? 1000 : 1 << 30;
        const rescue_prune_params p = rescue_prune_params::defaults(19);
        if (J.hyp < 0) run_one(J.ref, J.len1, J.q, J.len2, p, mh, (k & 1) == 0, st);
        else run_rel(J.ref, J.len1, J.q, J.len2, p, J.hyp, mh, (k & 1) == 0, (k / 5) % 4 == 3, st);
        k++;
    }
    return 0;
}

/* Kill switch for the relation filter: time and yield per variant over the dumped jobs. */
int run_timem(int mh, int argc, char **argv)
{
    std::vector<std::vector<uint8_t>> bufs;
    std::vector<MJob> jobs;
    if (!load_mdumps(argc, argv, bufs, jobs)) return 1;
    static rescue_prune_scratch sr;
    for (int v = 0; v < 4; v++) {   // 0 SIMD relation, 1 converted copies, 2 scalar relation, 3 none
        long kind[3] = {0, 0, 0};
        double rows_in = 0, rows_kept = 0, cells_in = 0, cells_kept = 0;
        const auto t0 = std::chrono::steady_clock::now();
        for (const MJob &J : jobs) {
            rescue_prune_params p = rescue_prune_params::defaults(19);
            int hb = -1, he = -1, k = RESCUE_PRUNE_FULL;
            if (J.hyp >= 0 && v != 3) {
                if (v == 1) p.set_meth(J.hyp);
                else p.set_meth_rel(J.hyp);
                k = v == 2 ? (J.len1 >= 5 && J.len2 >= 5 ? rescue_prune_window_scalar(J.ref, J.len1, J.q, J.len2, p, mh, sr, &hb, &he)
                                                        : RESCUE_PRUNE_FULL)
                           : rescue_prune_window(J.ref, J.len1, J.q, J.len2, p, mh, &hb, &he);
            }
            const int kept = k == RESCUE_PRUNE_B1 ? 0 : k == RESCUE_PRUNE_B2 ? he - hb + 1 : J.len1;
            kind[k]++;
            rows_in += J.len1; rows_kept += kept;
            cells_in += (double)J.len1 * J.len2; cells_kept += (double)kept * J.len2;
        }
        const double sec = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
        static const char *name[4] = {"simd_rel", "conv", "scalar_rel", "none"};
        printf("timem %-10s mh=%d jobs=%zu full=%ld b1=%ld b2=%ld rows_kept=%.4f cells_kept=%.4f sec=%.3f\n", name[v], mh,
               jobs.size(), kind[0], kind[1], kind[2], rows_kept / rows_in, cells_kept / cells_in, sec);
    }
    return 0;
}

void run_fuzz(long n, unsigned seed, Stats &st)
{
    std::mt19937 rng(seed);
    auto rnd = [&](int lo, int hi) { return std::uniform_int_distribution<int>(lo, hi)(rng); };
    std::vector<uint8_t> q, ref;
    for (long it = 0; it < n; it++) {
        const int mode = rnd(0, 10);
        /* relation jobs: relx = T (hyp 1) or A (hyp 0); mode 10 = a long mate rich in relx, up to
         * poly-relx, which can exceed the relation's entry cap */
        const bool rel = mode == 10 || rnd(0, 2) == 0;
        const int hyp = rnd(0, 1), relx = hyp ? 3 : 0;
        const int len2 = mode == 0 ? rnd(1, 12) : mode == 10 ? rnd(300, 1000) : rnd(5, 320);
        q.resize(len2);
        /* query: random, low complexity (period 1-4), or random with repeats of itself */
        const int period = mode >= 7 && mode <= 9 ? rnd(1, 4) : 0;
        const int richp = mode == 10 ? rnd(30, 100) : 0;   // % relx
        for (int j = 0; j < len2; j++)
            q[j] = richp ? (rnd(0, 99) < richp ? (uint8_t)relx : (uint8_t)rnd(0, 3))
                 : period ? (uint8_t)((j % period) * 7 % 4) : (uint8_t)rnd(0, 3);
        if (mode == 6)
            for (int j = rnd(8, 24); j < len2; j++) q[j] = q[j - rnd(1, 8) > 0 ? j - rnd(1, 8) : 0];
        const int len1 = mode == 0 ? rnd(1, 40) : mode == 1 ? rnd(4000, 4100) : rnd(5, 2400);
        ref.resize(len1);
        for (int i = 0; i < len1; i++) ref[i] = period && rnd(0, 3) ? (uint8_t)((i % period) * 7 % 4) : (uint8_t)rnd(0, 3);
        /* embed mutated copies of the query (the mate's true placement and repeats) */
        const int copies = rnd(0, mode == 5 ? 20 : 4);
        for (int c = 0; c < copies && len1 > len2; c++) {
            const int at = rnd(0, len1 - len2);
            const int subs = rnd(0, 8), indel = rnd(-3, 3);
            for (int j = 0; j < len2 && at + j < len1; j++) ref[at + j] = q[j];
            for (int s = 0; s < subs; s++) ref[at + rnd(0, len2 - 1)] = (uint8_t)rnd(0, 3);
            if (indel > 0) {
                const int x = at + rnd(0, len2 - 1);
                for (int i = len1 - 1; i >= x + indel; i--) ref[i] = ref[i - indel];
            } else if (indel < 0) {
                const int x = at + rnd(0, len2 - 1);
                for (int i = x; i + (-indel) < len1; i++) ref[i] = ref[i + (-indel)];
            }
        }
        if (rnd(0, 49) == 0) ref[rnd(0, len1 - 1)] = 4;   // N in the window
        if (rnd(0, 49) == 0) q[rnd(0, len2 - 1)] = 4;     // N in the query
        const int mh = rnd(0, 3) == 0 ? rnd(0, 3000) : rnd(0, 1) ? 1000 : 1 << 30;
        const int minsc = rnd(0, 1) ? 19 : rnd(5, 40);
        const int which = rnd(0, 1) ? -1 : rnd(0, NSCORINGS - 1);
        if (rel) run_rel(ref.data(), len1, q.data(), len2, scoring(which, minsc), hyp, mh, rnd(0, 1) == 0, rnd(0, 3) == 0, st);
        else run_one(ref.data(), len1, q.data(), len2, scoring(which, minsc), mh, rnd(0, 1) == 0, st);
    }
}

}  // namespace

int main(int argc, char **argv)
{
    if (argc < 3) {
        fprintf(stderr, "usage: %s dump|dump16|dumpm <files...> | timem <mh> <files...> | fuzz <n> <seed>\n", argv[0]);
        return 2;
    }
    Stats st;
    const std::string mode = argv[1];
    if (mode == "dump" || mode == "dump16") {
        if (run_dumps(argc - 2, argv + 2, mode == "dump16", st)) return 2;
    } else if (mode == "dumpm") {
        if (run_mdumps(argc - 2, argv + 2, st)) return 2;
    } else if (mode == "timem" && argc >= 4) {
        return run_timem(atoi(argv[2]), argc - 3, argv + 3) ? 2 : 0;
    } else if (mode == "fuzz" && argc >= 4) {
        run_fuzz(atol(argv[2]), (unsigned)atol(argv[3]), st);
    } else {
        fprintf(stderr, "unknown mode %s\n", argv[1]);
        return 2;
    }
    printf("eq %s: jobs=%ld simd=%ld scalar_fallback=%ld full=%ld b1=%ld b2=%ld rel=%ld rel_b2=%ld MISMATCHES=%ld\n",
           mode.c_str(), st.jobs, st.simd, st.fallback, st.kind[0], st.kind[1], st.kind[2], st.rel, st.rel_b2, st.bad);
    return st.bad ? 1 : 0;
}
