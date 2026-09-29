/* Equivalence check for the NEON rescue-pruning filter (src/rescue_prune_neon.h) against the scalar
 * filter, through the public rescue_prune_window() API.
 *
 * For every job at the NEON filter's threshold (rescue_prune_neon::MINSC): the NEON path must return
 * the scalar path's (kind, hb, he), and on a NEON B2 its view must match the scalar arrays -- cnt and
 * minrow (where cnt > 0) on [0, nd), bnd16 == 5 + fwd + bwd - (cnt - 1), the mw / hw bits
 * (bnd >= MINSC, and cnt > 0) on [0, nd) -- and its component list (runs of bnd >= MINSC with a hit:
 * a, b, ub, i0, dmax) with its count. Half of the jobs run twice in a row, the second time from a
 * copy of the bytes in other buffers, which exercises the filter's repeat memo: the repeat must
 * return the same decision and view as the first call.
 *
 *   rescue_prune_eq fuzz <n> <seed>   n generated jobs (repeats, low complexity, N, tiny and long
 *                                     windows, windows with more components than the filter stores,
 *                                     hit-gate values), each sometimes followed by near-repeats the
 *                                     memo must not match (a prefix of the window, its last base
 *                                     changed)
 *   rescue_prune_eq dump <files...>   jobs from rescue_band_harness dumps (int32 len1, len2; ref;
 *                                     query), each at the hull (400) and banded (1000) gates and
 *                                     with the gate open
 *
 * Exit status 1 on any disagreement, and when no job took the NEON path (a vacuous run). aarch64
 * only: elsewhere rescue_prune_window() is the scalar filter itself. Build: make rescue-prune-eq.
 *
 * RPE_NEGATIVE_CONTROL=1 shifts the scalar filter's hb by one on the first B2 job, so a working
 * check must report a MISMATCH and exit 1; CI runs it first to prove the comparison still bites. */
#include "../src/rescue_prune.h"
#include "../src/rescue_env.h"

#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <random>
#include <string>
#include <vector>

#if !defined(__aarch64__)
#error "rescue_prune_eq checks the NEON filter and builds on aarch64 only"
#endif

namespace {

using rescue_prune_neon::Comp;
using rescue_prune_neon::MINSC;

struct Stats {
    long jobs = 0, neon = 0, over_cap = 0, repeats = 0, memo_hits = 0, bad = 0, kind[3] = {0, 0, 0};
};

/* The scalar fields of two views of one job. Their arrays alias the same per-thread scratch, so
 * comparing the pointers would prove nothing; the contents are checked against the scalar filter
 * by check_view on every call, the repeat included. */
bool same_view(const rescue_prune_view &a, const rescue_prune_view &b)
{
    return a.nd == b.nd && a.off == b.off && a.ncomp == b.ncomp && a.ncomp_stored == b.ncomp_stored;
}

/* The NEON view v of a B2 decision against the scalar filter's arrays in r for the same job. */
bool check_view(int len1, int len2, const rescue_prune_view &v, const rescue_prune_scratch &r, std::string &why)
{
    const rescue_prune_view sv = rescue_prune_scalar_view(r, len1, len2);
    const int nd = sv.nd;
    if (v.nd != nd || v.off != sv.off) { why = "nd / off"; return false; }
    if (!v.mw || !v.hw || !v.comps || v.ncomp < 0) { why = "NEON view incomplete"; return false; }
    auto bit = [](const uint64_t *w, int x) { return (int)(w[x >> 6] >> (x & 63) & 1); };
    auto bnd = [&](int x) { return 5 + r.fwd[x] + r.bwd[x] - ((int)r.cnt[x] - 1); };
    for (int x = 0; x < nd; x++) {
        if (v.cnt[x] != r.cnt[x]) { why = "cnt at " + std::to_string(x); return false; }
        if (r.cnt[x] && v.minrow[x] != r.minrow[x]) { why = "minrow at " + std::to_string(x); return false; }
        if (v.bnd16[x] != bnd(x)) { why = "bnd at " + std::to_string(x); return false; }
        if (bit(v.mw, x) != (bnd(x) >= MINSC)) { why = "mw at " + std::to_string(x); return false; }
        if (bit(v.hw, x) != (bnd(x) >= MINSC && r.cnt[x] > 0)) { why = "hw at " + std::to_string(x); return false; }
    }
    std::vector<Comp> want;
    for (int x = 0; x < nd;) {
        if (bnd(x) < MINSC) { x++; continue; }
        const int a = x;
        int ub = INT16_MIN, i0 = INT16_MAX, dmax = -1;
        for (; x < nd && bnd(x) >= MINSC; x++) {
            ub = std::max(ub, bnd(x));
            if (r.cnt[x]) { i0 = std::min(i0, (int)r.minrow[x]); dmax = x; }
        }
        if (dmax >= 0) want.push_back(Comp{a, x, ub, i0, dmax});
    }
    if (v.ncomp != (int)want.size()) {
        why = "ncomp " + std::to_string(v.ncomp) + " vs " + std::to_string(want.size());
        return false;
    }
    if (v.ncomp_stored != std::min(v.ncomp, (int)rescue_prune_neon::NeonScratch::COMP_CAP)) {
        why = "ncomp_stored";
        return false;
    }
    for (int k = 0; k < v.ncomp_stored; k++) {
        const Comp &a = v.comps[k], &b = want[(size_t)k];
        if (a.a != b.a || a.b != b.b || a.ub != b.ub || a.i0 != b.i0 || a.dmax != b.dmax) {
            why = "comp " + std::to_string(k);
            return false;
        }
    }
    return true;
}

void run_one(const uint8_t *ref, int len1, const uint8_t *q, int len2, int max_hits, bool twice, Stats &st)
{
    static rescue_prune_scratch r;
    /* rescue_prune_window's own guards return FULL before either filter runs */
    const bool guarded = len1 < 5 || len2 < 5 || len2 > rescue_prune_scratch::QCAP || len1 > 30000;
    int h0 = -1, e0 = -1;
    const int k0 = guarded ? RESCUE_PRUNE_FULL
                           : rescue_prune_window_scalar(ref, len1, q, len2, MINSC, max_hits, r, &h0, &e0);
    static const bool negative_control = rescue_env_opt_in("RPE_NEGATIVE_CONTROL");
    static bool injected = false;
    if (negative_control && !injected && k0 == RESCUE_PRUNE_B2) { h0++; injected = true; }
    int kind[2] = {-1, -1}, hb[2] = {-1, -1}, he[2] = {-1, -1};
    rescue_prune_view view[2];
    std::vector<uint8_t> ref2, q2;
    for (int rep = 0; rep < (twice ? 2 : 1); rep++) {
        const uint8_t *R = ref, *Q = q;
        if (rep == 1) {   // same bytes, other buffers: only the contents can match the memo
            ref2.assign(ref, ref + len1);
            q2.assign(q, q + len2);
            R = ref2.data(); Q = q2.data();
            st.repeats++;
        }
        const uint64_t hits0 = rescue_prune_memo_hits();
        kind[rep] = rescue_prune_window(R, len1, Q, len2, MINSC, max_hits, &hb[rep], &he[rep], &view[rep]);
        const bool hit = rescue_prune_memo_hits() != hits0;
        st.memo_hits += hit;
        st.jobs++;
        if (rep == 0) st.kind[kind[0]]++;
        const bool neon = kind[rep] == RESCUE_PRUNE_B2 && view[rep].bnd16 != nullptr;
        st.neon += neon;
        st.over_cap += neon && view[rep].ncomp > rescue_prune_neon::NeonScratch::COMP_CAP;
        std::string why;
        bool ok = kind[rep] == k0 && (k0 != RESCUE_PRUNE_B2 || (hb[rep] == h0 && he[rep] == e0));
        if (!ok) why = "decision (" + std::to_string(kind[rep]) + "," + std::to_string(hb[rep]) + "," +
                       std::to_string(he[rep]) + ") vs scalar (" + std::to_string(k0) + "," + std::to_string(h0) +
                       "," + std::to_string(e0) + ")";
        if (ok && neon) ok = check_view(len1, len2, view[rep], r, why);
        if (ok && rep == 1 && !same_view(view[0], view[1])) { ok = false; why = "repeat view differs"; }
        /* The memo may answer only a byte-for-byte repeat of the previous call that reached the
         * NEON filter (the wrapper's guards return before it, leaving its memo alone). */
        static std::vector<uint8_t> last_ref, last_q;
        static int last_mh = -1;
        const bool repeat = max_hits == last_mh && last_ref.size() == (size_t)len1 && last_q.size() == (size_t)len2
                            && std::equal(R, R + len1, last_ref.begin()) && std::equal(Q, Q + len2, last_q.begin());
        if (ok && hit && !repeat) { ok = false; why = "memo answered a job that does not repeat the previous one"; }
        if (!guarded) {
            last_ref.assign(R, R + len1);
            last_q.assign(Q, Q + len2);
            last_mh = max_hits;
        }
        if (!ok) {
            if (st.bad < 20)
                fprintf(stderr, "MISMATCH len1=%d len2=%d max_hits=%d rep=%d: %s\n", len1, len2, max_hits, rep,
                        why.c_str());
            st.bad++;
        }
    }
}

int run_dumps(int argc, char **argv, Stats &st)
{
    for (int a = 0; a < argc; a++) {
        FILE *fp = fopen(argv[a], "rb");
        if (!fp) { perror(argv[a]); return 1; }
        long k = 0;
        int32_t h[2];
        std::vector<uint8_t> ref, q;
        size_t got;
        /* A truncated or corrupt dump is a hard error, as in rescue_band_harness's reader: silently
         * dropping jobs would shrink the check without saying so. */
        while ((got = fread(h, 1, sizeof h, fp)) == sizeof h) {
            if (h[0] < 0 || h[1] < 0 || h[0] > (1 << 24) || h[1] > (1 << 16)) {
                fprintf(stderr, "%s: corrupt record %ld (len1 %d, len2 %d)\n", argv[a], k, h[0], h[1]);
                fclose(fp);
                return 1;
            }
            ref.resize((size_t)h[0]);
            q.resize((size_t)h[1]);
            if (fread(ref.data(), 1, (size_t)h[0], fp) != (size_t)h[0]
                || fread(q.data(), 1, (size_t)h[1], fp) != (size_t)h[1]) {
                fprintf(stderr, "%s: truncated record %ld\n", argv[a], k);
                fclose(fp);
                return 1;
            }
            const int mh = k % 3 == 0 ? 400 : k % 3 == 1 ? 1000 : 1 << 30;
            run_one(ref.data(), h[0], q.data(), h[1], mh, (k & 1) == 0, st);
            k++;
        }
        const bool err = ferror(fp) != 0;
        fclose(fp);
        if (err || got != 0) {
            fprintf(stderr, "%s: %s after record %ld\n", argv[a], err ? "read error" : "trailing partial header", k);
            return 1;
        }
    }
    return 0;
}

void run_fuzz(long n, unsigned seed, Stats &st)
{
    std::mt19937 rng(seed);
    auto rnd = [&](int lo, int hi) { return std::uniform_int_distribution<int>(lo, hi)(rng); };
    std::vector<uint8_t> q, ref;
    for (long it = 0; it < n; it++) {
        const int mode = rnd(0, 9);
        const int len2 = mode == 0 ? rnd(1, 12) : rnd(5, 320);
        q.resize((size_t)len2);
        /* query: random, low complexity (period 1-4), or random with repeats of itself */
        const int period = mode >= 7 ? rnd(1, 4) : 0;
        for (int j = 0; j < len2; j++) q[(size_t)j] = period ? (uint8_t)((j % period) * 7 % 4) : (uint8_t)rnd(0, 3);
        if (mode == 6)
            for (int j = rnd(8, 24); j < len2; j++) q[(size_t)j] = q[(size_t)std::max(0, j - rnd(1, 8))];
        /* windows: tiny, near the NEON filter's capacity (past it the scalar filter decides), long
         * enough for many components, or typical */
        const int len1 = mode == 0 ? rnd(1, 40) : mode == 1 ? rnd(3700, 4200) : mode == 5 ? rnd(2800, 3700)
                                                                                           : rnd(5, 2400);
        ref.resize((size_t)len1);
        for (int i = 0; i < len1; i++)
            ref[(size_t)i] = period && rnd(0, 3) ? (uint8_t)((i % period) * 7 % 4) : (uint8_t)rnd(0, 3);
        /* mode 5: one 20-base query segment planted every `step` rows, one component per copy, which
         * is how a window gets more components than the filter stores (COMP_CAP) */
        if (mode == 5 && len2 >= 20) {
            const int step = rnd(21, 30), qat = rnd(0, len2 - 20);
            for (int at = rnd(0, step - 1); at + 20 <= len1; at += step)
                for (int j = 0; j < 20; j++) ref[(size_t)(at + j)] = q[(size_t)(qat + j)];
        }
        /* mutated copies of the query (the mate's true placement and repeats) */
        const int copies = mode == 5 ? 0 : rnd(0, 4);
        for (int c = 0; c < copies && len1 > len2; c++) {
            const int at = rnd(0, len1 - len2);
            const int subs = rnd(0, 8), indel = rnd(-3, 3);
            const int seg = len2;
            for (int j = 0; j < seg && at + j < len1; j++) ref[(size_t)(at + j)] = q[(size_t)j];
            for (int s = 0; s < subs; s++) ref[(size_t)(at + rnd(0, seg - 1))] = (uint8_t)rnd(0, 3);
            if (indel > 0) {
                const int x = at + rnd(0, seg - 1);
                for (int i = len1 - 1; i >= x + indel; i--) ref[(size_t)i] = ref[(size_t)(i - indel)];
            } else if (indel < 0) {
                const int x = at + rnd(0, seg - 1);
                for (int i = x; i - indel < len1; i++) ref[(size_t)i] = ref[(size_t)(i - indel)];
            }
        }
        if (rnd(0, 49) == 0) ref[(size_t)rnd(0, len1 - 1)] = 4;   // N in the window
        if (rnd(0, 49) == 0) q[(size_t)rnd(0, len2 - 1)] = 4;     // N in the query
        const int mh = rnd(0, 3) == 0 ? rnd(0, 3000) : rnd(0, 1) ? 1000 : 1 << 30;
        run_one(ref.data(), len1, q.data(), len2, mh, rnd(0, 1) == 0, st);
        /* near-repeats straight after, which the repeat memo must not mistake for the job above: a
         * prefix of the window (the memo compares only the new job's bytes), and the window with
         * its last base changed */
        if (len1 > 5 && rnd(0, 3) == 0)
            run_one(ref.data(), len1 - rnd(1, std::min(16, len1 - 5)), q.data(), len2, mh, false, st);
        if (rnd(0, 3) == 0) {
            ref[(size_t)(len1 - 1)] = (uint8_t)((ref[(size_t)(len1 - 1)] + 1) & 3);
            run_one(ref.data(), len1, q.data(), len2, mh, false, st);
        }
    }
}

}  // namespace

int main(int argc, char **argv)
{
    if (argc < 3) {
        fprintf(stderr, "usage: %s fuzz <n> <seed> | dump <files...>\n", argv[0]);
        return 2;
    }
    Stats st;
    const std::string mode = argv[1];
    if (mode == "dump") {
        if (run_dumps(argc - 2, argv + 2, st)) return 2;
    } else if (mode == "fuzz" && argc >= 4) {
        run_fuzz(atol(argv[2]), (unsigned)atol(argv[3]), st);
    } else {
        fprintf(stderr, "unknown mode %s\n", argv[1]);
        return 2;
    }
    printf("eq %s: jobs=%ld repeats=%ld memo_hits=%ld neon_b2=%ld over_comp_cap=%ld full=%ld b1=%ld b2=%ld "
           "MISMATCHES=%ld\n",
           mode.c_str(), st.jobs, st.repeats, st.memo_hits, st.neon, st.over_cap, st.kind[RESCUE_PRUNE_FULL],
           st.kind[RESCUE_PRUNE_B1], st.kind[RESCUE_PRUNE_B2], st.bad);
    if (st.neon == 0) fprintf(stderr, "FAIL: no job took the NEON filter's B2 path; the run compared nothing\n");
    /* Repeats run straight after their job, so a memo that never answers one has stopped working
     * (a correctness no-op, but the repeat checks above would then test nothing). */
    const bool memo_dead = st.repeats > 0 && st.memo_hits == 0;
    if (memo_dead) fprintf(stderr, "FAIL: %ld repeats and no memo hit; the repeat memo never engaged\n", st.repeats);
    return st.bad || st.neon == 0 || memo_dead ? 1 : 0;
}
