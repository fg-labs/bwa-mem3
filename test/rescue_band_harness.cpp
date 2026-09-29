/* Exactness + timing harness for the banded mate rescue (src/rescue_band.{h,cpp}).
 * Build: `make rescue-band-harness` (arch=arm64 on aarch64, arch=avx2 on x86).
 *
 *   rescue_band_harness eq   <gen-jobs> <seed> [dump files...]   exactness (exit 1 on any mismatch)
 *   rescue_band_harness time <reps> <stride> <dump files...>     single-thread pass-0 timing
 *
 * eq: for every job, the PRODUCTION pipeline -- rescue_prune_window, RescueBandBatch::plan/commit,
 * the length sort, partition, kswv phase 0 on the non-banded pairs, run_pass0 on the banded ones,
 * then the mem_sam_pe_batch_run post-processing and kswv phase 1 -- is compared against the
 * pre-pruning pipeline (full window, kswv phase 0 + phase 1), field by field: score (pass/fail
 * when the truth fails), te, qe, score2, tb, qb (offsets applied). The production pipeline includes
 * the banded pass 1 (take_pass1 / run_pass1) of the banded parents and, under the default
 * BWA3_RESCUE_BAND_P1=2, of every other pair (hull or full window, kswv pass 0), so tb / qb check
 * it against kswv phase 1 on the full window; the summary counts how many jobs took it.
 * A scalar full-window DP with kswv's score2 semantics (scalar_dp / scalar_score2 below) plus a scalar
 * reversed-prefix DP for (tb, qb) cross-check the truth on every --scalar-stride-th job.
 * Generated classes: random windows with mutated mate copies, edge copies around te +- S (score2
 * parity / zeroing), tiny windows, N bases (must be FULL), tandem repeats (hit-dense, many
 * components), in-zone secondary copies (round-2 triggers), ragged len2, windows past the NEON
 * filter's capacity (scalar view), more components than the cap, and pass-1 adversaries (class 11:
 * gaps at the band-edge bound, start ties, qe next to the pad columns, te at the window edges).
 * Every B2 job's NEON view is also checked against the scalar filter's view of the same job (same
 * components at several thresholds), and the run fails if nothing was banded. Env knobs are the
 * production ones; the caller sets e.g.
 * BWA3_RESCUE_BAND_COST=100000000 (band every B2 parent) and BWA3_RESCUE_PRUNE_MAX_HITS. RB_MINSC
 * (default 19), RB_SCALAR_STRIDE (default 0: no scalar cross-check) and RB_NEGATIVE_CONTROL are
 * harness-only; RB_NEGATIVE_CONTROL=1 shifts the production te of the first banded passing job by
 * one, so a working comparison must report a MISMATCH and exit 1 (CI runs it first).
 *
 * Dump files: one record per prune-eligible rescue job, as the aligner saw it (full window and
 * oriented mate): int32 len1, int32 len2, len1 reference bytes, len2 query bytes (2-bit codes,
 * 4 = N). The aligner does not write them; capture new ones by instrumenting mem_matesw_batch_pre
 * where it calls rescue_prune_window. */
#include "rescue_band.h"
#include "rescue_env.h"
#include "simd_dispatch.h"

#include <algorithm>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <memory>
#include <random>
#include <string>
#include <vector>

static const int MINSC_DEFAULT = 19;

/* The harness is single-threaded, so one batch stands in for the aligner's per-tid one. */
static RescueBandBatch &rescue_band_batch()
{
    static RescueBandBatch b;
    return b;
}

struct Job { std::vector<uint8_t> ref, q; int cls; };
struct Res { int score, te, qe, score2, tb, qb; };

/* ---------------------------------------------------------------------------------------- */
/* scalar reference                                                                          */
/* ---------------------------------------------------------------------------------------- */
static void scalar_dp(const Job &jb, int &S, int &te, int &qe, std::vector<int> &R)
{
    const int len1 = (int)jb.ref.size(), len2 = (int)jb.q.size(), quanta = kswv_query_quantum8(len2);
    std::vector<int> Hp(quanta, 0), Hc(quanta, 0), E(quanta, 0);
    R.assign(len1, 0);
    S = 0; te = -1; qe = 0;
    for (int i = 0; i < len1; i++) {
        int f = 0, rmax = 0;
        for (int j = 0; j < quanta; j++) {
            const int diag = j > 0 ? Hp[j - 1] : 0;
            const int sc = (j >= len2) ? 0 : (jb.ref[i] >= 4 || jb.q[j] >= 4) ? -1 : (jb.ref[i] == jb.q[j] ? 1 : -4);
            int h = std::max(0, diag + sc);
            h = std::max(h, std::max(E[j], f));
            Hc[j] = h; rmax = std::max(rmax, h);
            E[j] = std::max(0, std::max(E[j] - 1, h - 7));
            f = std::max(0, std::max(f - 1, h - 7));
        }
        R[i] = rmax;
        if (rmax > S) { S = rmax; te = i; for (int j = 0; j < quanta; j++) if (Hc[j] == rmax) { qe = j; break; } }
        std::swap(Hp, Hc);
    }
}
/* kswv phase 1 in scalar form: the DP over ref[te..0] x q[qe..0] plus the pad columns
 * [qe + 1, quanta); first row whose max reaches S, first column holding S there. */
static void scalar_pass1(const Job &jb, int S, int te, int qe, int &tb, int &qb)
{
    const int len2 = qe + 1, quanta = (len2 + 15) / 16 * 16;
    std::vector<int> Hp(quanta, 0), Hc(quanta, 0), E(quanta, 0);
    tb = qb = -1;
    for (int r = 0; r <= te; r++) {
        const int rb = jb.ref[te - r];
        int f = 0, rmax = 0;
        for (int c = 0; c < quanta; c++) {
            const int diag = c > 0 ? Hp[c - 1] : 0;
            const int qc = c < len2 ? jb.q[qe - c] : 0;
            const int sc = (c >= len2) ? 0 : (rb >= 4 || qc >= 4) ? -1 : (rb == qc ? 1 : -4);
            int h = std::max(0, diag + sc);
            h = std::max(h, std::max(E[c], f));
            Hc[c] = h; rmax = std::max(rmax, h);
            E[c] = std::max(0, std::max(E[c] - 1, h - 7));
            f = std::max(0, std::max(f - 1, h - 7));
        }
        if (rmax >= S) {
            for (int c = 0; c < quanta; c++) if (Hc[c] == rmax) { tb = te - r; qb = qe - c; break; }
            return;
        }
        std::swap(Hp, Hc);
    }
}
static int scalar_score2(const std::vector<int> &R, int S, int te, int minsc)
{
    const int n = (int)R.size(), low = te - S, high = te + S;
    int s2 = -1, bs = -1, bp = -2;
    for (int i = 0; i < n; i++) {
        const int v = (i + 1 < n && R[i + 1] > R[i]) ? 0 : R[i];
        if (v < minsc) continue;
        if (bp + 1 != i) { if (bp >= 0 && (bp < low || bp > high) && bs > s2) s2 = bs; bs = v; bp = i; }
        else if (bs < v) { bs = v; bp = i; }
    }
    if (bp >= 0 && (bp < low || bp > high) && bs > s2) s2 = bs;
    return s2;
}

/* ---------------------------------------------------------------------------------------- */
/* generators                                                                                 */
/* ---------------------------------------------------------------------------------------- */
static std::mt19937_64 rng;
static int rnd(int n) { return n <= 1 ? 0 : (int)(rng() % (uint64_t)n); }
static int rndr(int a, int b) { return a + rnd(b - a + 1); }
static double unif() { return (rng() >> 11) * (1.0 / 9007199254740992.0); }

static std::vector<uint8_t> mutate(const std::vector<uint8_t> &s, int a, int b, double mm, double ind)
{
    std::vector<uint8_t> o;
    a = std::max(a, 0);   // callers may ask for a copy longer than a short mate: clip to it
    b = std::min(b, (int)s.size());
    for (int i = a; i < b; i++) {
        if (unif() < ind) {
            const int l = 1 + rnd(4);
            if (rnd(2)) { i += l - 1; continue; }
            for (int k = 0; k < l; k++) o.push_back((uint8_t)rnd(4));
        }
        uint8_t c = s[i];
        if (unif() < mm) c = (uint8_t)((c + 1 + rnd(3)) & 3);
        o.push_back(c);
    }
    return o;
}
static void plant(std::vector<uint8_t> &ref, const std::vector<uint8_t> &m, int pos)
{
    for (size_t t = 0; t < m.size(); t++) {
        const long p = pos + (long)t;
        if (p >= 0 && p < (long)ref.size()) ref[p] = m[t];
    }
}

static Job gen(int cls)
{
    Job J; J.cls = cls;
    int len2 = rnd(4) == 0 ? rndr(17, 240) : rndr(100, 160);
    if (rnd(6) == 0) len2 = 16 * rndr(2, 15);   // exact multiples of 16 (no pad columns)
    std::vector<uint8_t> &q = J.q, &r = J.ref;
    q.resize(len2);
    for (auto &b : q) b = (uint8_t)rnd(4);
    int len1 = rndr(len2 + 5, 1500);
    switch (cls) {
    case 0: {   // random window, 0-4 mutated (partial) copies anywhere, incl. window edges
        r.resize(len1); for (auto &b : r) b = (uint8_t)rnd(4);
        for (int c = rnd(5); c--;) {
            const int a = rnd(len2 / 2 + 1), b = rndr(std::min(len2, a + 20), len2);
            plant(r, mutate(q, a, b, 0.1 * unif() * unif(), 0.02 * unif()), rndr(-40, len1 - 10));
        }
        break;
    }
    case 1: {   // edge fuzz: a full copy + partial copies ending near te-S, te+S 
        len1 = 3 * len2 + 200 + rnd(400);
        r.resize(len1); for (auto &b : r) b = (uint8_t)rnd(4);
        const int P = rnd(len1 - len2), te = P + len2 - 1, S = len2;
        for (int j = 0; j < len2; j++) r[P + j] = q[j];
        for (int c = 1 + rnd(8); c--;) {
            const int sl = 20 + rnd(len2 - 19);
            const int q0 = rnd(2) ? len2 - sl : rnd(len2 - sl + 1);
            int end;
            switch (rnd(3)) {
                case 0: end = te + S + rnd(21) - 10; break;
                case 1: end = te - S + rnd(21) - 10; break;
                default: end = te + 1 + rnd(S + 30); break;
            }
            const int pos = end - sl + 1;
            if (pos < 0 || end >= len1 || (pos <= te && end >= P)) continue;
            plant(r, mutate(q, q0, q0 + sl, 0.01 * rnd(10), 0.01 * rnd(4)), pos);
        }
        break;
    }
    case 2: {   // tiny windows / tiny mates
        len2 = rndr(5, 60); q.resize(len2); for (auto &b : q) b = (uint8_t)rnd(4);
        len1 = rndr(19, 120);
        r.resize(len1); for (auto &b : r) b = (uint8_t)rnd(4);
        if (rnd(3)) plant(r, mutate(q, 0, len2, 0.05 * unif(), 0.0), rndr(-len2 / 2, len1 - 3));
        break;
    }
    case 3: {   // N bases (the filter must bail to FULL) plus a copy
        r.resize(len1); for (auto &b : r) b = (uint8_t)rnd(4);
        plant(r, mutate(q, 0, len2, 0.02, 0.005), rnd(len1));
        if (rnd(2)) for (int t = rndr(1, 4); t--;) r[rnd(len1)] = 4;
        else for (int t = rndr(1, 3); t--;) q[rnd(len2)] = 4;
        break;
    }
    case 4: {   // tandem repeats: repeat-rich window and mate (hit-dense, many components)
        const int per = rndr(1, 8);
        std::vector<uint8_t> unit(per); for (auto &b : unit) b = (uint8_t)rnd(4);
        r.resize(len1);
        for (int i = 0; i < len1; i++) r[i] = unit[i % per];
        for (int i = 0; i < len1; i++) if (unif() < 0.03) r[i] = (uint8_t)rnd(4);
        const int rl = rndr(10, len2);   // repeat stretch inside the mate
        const int rs = rnd(len2 - rl + 1);
        for (int j = rs; j < rs + rl; j++) q[j] = unit[(j + rnd(2)) % per];
        if (rnd(2)) plant(r, mutate(q, 0, len2, 0.02, 0.01), rnd(len1));
        for (int t = rnd(6); t--;) {   // unique stretches between repeats
            const int a = rnd(len1), l = rndr(5, 200);
            for (int i = a; i < std::min(len1, a + l); i++) r[i] = (uint8_t)rnd(4);
        }
        break;
    }
    case 5: {   // round-2 triggers: strong main copy + secondary copies (score ~19..80) in and
                // around the zone (te - S, te + S), so T1 = ub2 / 2 > 19 and score2 < T1 is common
        len1 = std::max(len1, 3 * len2 + 100);
        r.resize(len1); for (auto &b : r) b = (uint8_t)rnd(4);
        const int P = rnd(len1 - len2);
        plant(r, mutate(q, 0, len2, 0.01 * rnd(3), 0.003), P);
        const int te = P + len2 - 1, S = len2;
        for (int c = rndr(1, 5); c--;) {
            const int sl = rndr(20, std::min(len2, 90));
            const int q0 = rnd(len2 - sl + 1);
            const int end = te - S + rnd(3 * S);
            plant(r, mutate(q, q0, q0 + sl, 0.02 * rnd(4), 0.005), end - sl + 1);
        }
        break;
    }
    case 7: {   // round 2 by construction: a mutated main copy, strong secondaries (ub ~40-80, so
                // T1 = ub2/2 > 19) ending INSIDE the zone but on FAR diagonals (separate components),
                // plus weak out-of-zone copies scoring 19..T1 (the exact score2 must come from round 2)
        len2 = rndr(100, 200); q.resize(len2); for (auto &b : q) b = (uint8_t)rnd(4);
        len1 = rndr(4 * len2, 6 * len2);
        r.resize(len1); for (auto &b : r) b = (uint8_t)rnd(4);
        const int P = rndr(len2, len1 - 2 * len2);
        plant(r, mutate(q, 0, len2, 0.02 + 0.06 * unif(), 0.003), P);
        const int te = P + len2 - 1, S = len2;
        for (int c = rndr(1, 3); c--;) {
            const int sl = rndr(60, 90);
            int q0, end;
            // ending in the zone (te - S, te + S) but on a diagonal at least len2 + sl from the
            // primary's: the 5-mer bound falls by one per hit-free diagonal from about len2 + sl,
            // so nearer copies stay inside the primary's component and never raise ub2
            q0 = rnd(len2 - sl + 1);
            const int dp = P, dfar = (int)((len2 + sl) * (1.0 + 0.1 * unif()));
            const int d = rnd(2) ? dp + dfar : dp - dfar;   // diagonal i - j of the copy
            end = d + q0 + sl - 1;
            if (end < te - S + 1 || end > te + S - 1 || end - sl + 1 < 0 || end >= len1) continue;
            plant(r, mutate(q, q0, q0 + sl, 0.01 * rnd(2), 0.0), end - sl + 1);
        }
        for (int c = rnd(4); c--;) {
            const int sl = rndr(19, 27), q0 = rnd(len2 - sl + 1);
            const int end = rnd(2) ? rndr(sl, std::max(sl, te - S - 1)) : rndr(std::min(len1 - 1, te + S + 1), len1 - 1);
            plant(r, mutate(q, q0, q0 + sl, 0.0, 0.0), end - sl + 1);
        }
        break;
    }
    case 8: {   // perfect / near-perfect single copy in a clean window: exercises the tight top band
                // and the "hull inside [te - S, te + S]" termination (score2 = -1), plus partial
                // weak copies that may or may not stay inside the zone
        r.resize(len1); for (auto &b : r) b = (uint8_t)rnd(4);
        const int P = rnd(std::max(1, len1 - len2));
        const int nm = rnd(3) ? 0 : rnd(2);
        std::vector<uint8_t> c(q);
        for (int t = 0; t < nm; t++) { const int x = rnd(len2); c[x] = (uint8_t)((c[x] + 1 + rnd(3)) & 3); }
        plant(r, c, P);
        for (int t = rnd(3); t--;) {
            const int sl = rndr(19, 30), q0 = rnd(len2 - sl + 1);
            plant(r, mutate(q, q0, q0 + sl, 0.0, 0.0), P + rndr(-2 * len2, 3 * len2));
        }
        break;
    }
    case 9: {   // long windows (past the NEON filter's capacity): the scalar filter decides and
                // plan() derives its components from the scalar (fwd / bwd) view
        len1 = rndr(4100, 12000);
        r.resize(len1); for (auto &b : r) b = (uint8_t)rnd(4);
        for (int c = rndr(1, 4); c--;)
            plant(r, mutate(q, 0, len2, 0.02 * rnd(4), 0.005), rnd(len1 - len2));
        for (int c = rnd(20); c--;) {
            const int sl = rndr(19, 40), q0 = rnd(len2 - sl + 1);
            plant(r, mutate(q, q0, q0 + sl, 0.0, 0.0), rnd(len1 - sl));
        }
        break;
    }
    case 10: {  // many short copies (20-120) on separate diagonals: past the component cap once
                // the hit gate is opened (BWA3_RESCUE_PRUNE_MAX_HITS); at the default gate they go FULL
        len1 = rndr(2000, 4000);
        r.resize(len1); for (auto &b : r) b = (uint8_t)rnd(4);
        for (int c = rndr(20, 120); c--;) {
            const int sl = rndr(19, 30), q0 = rnd(len2 - sl + 1);
            plant(r, mutate(q, q0, q0 + sl, 0.0, 0.0), rnd(len1 - sl));
        }
        if (rnd(2)) plant(r, mutate(q, 0, len2, 0.02, 0.003), rnd(len1 - len2));
        break;
    }
    case 11: {  // pass-1 adversaries (banded pass 1, rescue_band.h)
        r.resize(len1); for (auto &b : r) b = (uint8_t)rnd(4);
        const int mode = rnd(5);
        if (mode == 0) {   // qe next to the pad columns: len2 on / just past a multiple of 16
            len2 = 16 * rndr(2, 12) + rnd(2); q.resize(len2); for (auto &b : q) b = (uint8_t)rnd(4);
            len1 = std::max(len1, len2 + 5); r.resize(len1); for (auto &b : r) b = (uint8_t)rnd(4);
        }
        if (mode == 1) {   // start ties: the mate opens with a short tandem unit, so shifted starts tie
            const int per = rndr(1, 3), rl = rndr(4, 12);
            for (int j = 0; j < std::min(rl, len2); j++) q[j] = q[j % per];
        }
        /* The copy: 1-3 gaps whose total length sits at or near the band-edge bound
         * (2I + D <= qe - S - 5), a few mismatches, and sometimes a truncated tail (qe < len2 - 1). */
        std::vector<uint8_t> c;
        const int ng = rndr(mode == 2 ? 1 : 0, 3);
        std::vector<int> gpos(ng), glen(ng);
        for (int t = 0; t < ng; t++) { gpos[t] = rnd(len2); glen[t] = rnd(4) ? rndr(1, 6) : rndr(7, 40); }
        const int qend = rnd(3) ? len2 : rndr(std::min(len2, 25), len2);
        for (int j = 0; j < qend; j++) {
            for (int t = 0; t < ng; t++)
                if (gpos[t] == j) {
                    if (rnd(2)) for (int k = 0; k < glen[t]; k++) c.push_back((uint8_t)rnd(4));   // deletion (ref bases)
                    else { j += glen[t]; break; }                                                  // insertion (query bases)
                }
            if (j >= qend) break;
            uint8_t b = q[j];
            if (rnd(40) == 0) b = (uint8_t)((b + 1 + rnd(3)) & 3);
            c.push_back(b);
        }
        /* te at the window edges: the copy flush with the window start or end, or anywhere. */
        const int where = rnd(3);
        const int pos = where == 0 ? 0 : where == 1 ? len1 - (int)c.size() : rnd(std::max(1, len1 - (int)c.size()));
        plant(r, c, pos);
        if (mode == 3) plant(r, c, pos + rndr(-3, 3) + (rnd(2) ? (int)c.size() : 0));   // a tied second copy
        if (mode == 4) {   // the same copy shifted by a few bases on a nearby diagonal (overlapping ties)
            const int sh = rndr(1, 4);
            plant(r, std::vector<uint8_t>(c.begin(), c.begin() + std::min<size_t>(c.size(), 30)), pos - sh);
        }
        break;
    }
    case 12: {  // the 8-bit ceiling: mates of 241-250 bp (the longest the 8-bit path admits at the
                // default scoring) with a near-perfect copy, so band H reaches 241..250 -- the top of
                // the range the x86 cell's biased add must hold exactly (rescue_band_kernel_x86.h)
        len2 = rndr(241, 250); q.resize(len2); for (auto &b : q) b = (uint8_t)rnd(4);
        len1 = std::max(len1, len2 + rndr(5, 600)); r.resize(len1); for (auto &b : r) b = (uint8_t)rnd(4);
        const int P = rnd(std::max(1, len1 - len2));
        std::vector<uint8_t> c = rnd(3) ? q : mutate(q, 0, len2, 0.0, 0.004);   // perfect, or a 1-bp indel or two
        for (int t = rnd(3); t--;) { const int x = rnd(len2); c[x] = (uint8_t)((c[x] + 1 + rnd(3)) & 3); }
        plant(r, c, P);
        if (rnd(2)) plant(r, mutate(q, rnd(len2 / 2), len2, 0.02, 0.0), rnd(std::max(1, len1 - len2)));   // a weaker second copy (score2)
        break;
    }
    default: {  // two strong copies on nearby diagonals (merged / overlapping components)
        r.resize(len1); for (auto &b : r) b = (uint8_t)rnd(4);
        const int P = rnd(std::max(1, len1 - len2));
        plant(r, mutate(q, 0, len2, 0.01 * rnd(5), 0.01), P);
        plant(r, mutate(q, rnd(len2 / 2), len2, 0.01 * rnd(5), 0.01), P + rndr(-60, 60) + len2 / 2);
        break;
    }
    }
    return J;
}

static void load_dump(const char *path, std::vector<Job> &jobs, int stride)
{
    FILE *fp = fopen(path, "rb");
    if (!fp) { fprintf(stderr, "cannot open %s\n", path); exit(2); }
    long k = 0;
    int32_t h[2];
    if (stride < 1) stride = 1;
    while (fread(h, sizeof h, 1, fp) == 1) {
        /* A truncated or corrupt dump is a hard error: silently dropping jobs would shrink the
         * exactness check without saying so. */
        if (h[0] < 0 || h[1] < 0 || h[0] > (1 << 24) || h[1] > (1 << 16)) {
            fprintf(stderr, "%s: corrupt record %ld (len1 %d, len2 %d)\n", path, k, h[0], h[1]);
            exit(2);
        }
        Job J; J.cls = 100;
        J.ref.resize(h[0]); J.q.resize(h[1]);
        if (fread(J.ref.data(), 1, h[0], fp) != (size_t)h[0] || fread(J.q.data(), 1, h[1], fp) != (size_t)h[1]) {
            fprintf(stderr, "%s: truncated record %ld\n", path, k);
            exit(2);
        }
        if (k++ % stride == 0) jobs.push_back(std::move(J));
    }
    fclose(fp);
}

/* ---------------------------------------------------------------------------------------- */
/* production pipeline replicas                                                               */
/* ---------------------------------------------------------------------------------------- */
static void revseq(int l, uint8_t *s)
{
    for (int i = 0; i < l >> 1; i++) std::swap(s[i], s[l - 1 - i]);
}

struct Batch {
    std::vector<SeqPair> sp;
    std::vector<uint8_t> ref, qer;
    std::vector<kswr_t> aln;
    int n = 0;
    void clear() { sp.clear(); ref.clear(); qer.clear(); n = 0; }
    int add(const uint8_t *r, int len1, const uint8_t *q, int len2, int minsc)
    {
        SeqPair p;
        memset(&p, 0, sizeof p);
        p.idr = (int32_t)ref.size(); p.idq = (int32_t)qer.size();
        p.len1 = len1; p.len2 = len2;
        p.h0 = KSW_XSUBO | KSW_XSTART | KSW_XBYTE | minsc;
        p.regid = n; p.id = n; p.seqid = n;
        ref.insert(ref.end(), r, r + len1);
        qer.insert(qer.end(), q, q + len2);
        sp.push_back(p);
        return n++;
    }
    void finalize()
    {
        sp.resize(n + MAX_LINE_LEN + 64);
        ref.resize(ref.size() + 1024);
        qer.resize(qer.size() + 1024);
        aln.assign(n + 64, kswr_t{0, -1, -1, -1, -1, -1, -1});
    }
};

/* mem_sam_pe_batch_run (8-bit only) with the first nb pairs banded; p1_any offers every other pair
 * to the banded pass 1 too (the harness always uses the default scoring). The truth runs with
 * nb = 0 and p1_any = false: pure kswv. */
static void run_batch(Ikswv *k, Batch &b, int nb, bool p1_any)
{
    SeqPair *pairs = b.sp.data();
    kswr_t *aln = b.aln.data();
    for (int i = 0; i < b.n; i++) aln[i].tb = aln[i].qb = -1;
    k->getScores8(pairs + nb, b.ref.data(), b.qer.data(), aln, b.n - nb, 1, 0);
    if (nb) rescue_band_batch().run_pass0(pairs, nb, b.ref.data(), b.qer.data(), aln, k);
    int pos = 0;
    for (int i = 0; i < b.n; i++) {
        SeqPair sp = pairs[i];
        const kswr_t r = aln[sp.regid];
        const int xtra = sp.h0;
        if ((xtra & KSW_XSTART) == 0 || ((xtra & KSW_XSUBO) && r.score < (xtra & 0xffff))) continue;
        sp.h0 = KSW_XSTOP | r.score;
        sp.len2 = r.qe + 1;
        revseq(r.qe + 1, b.qer.data() + sp.idq);
        revseq(r.te + 1, b.ref.data() + sp.idr);
        if ((i < nb || p1_any) && rescue_band_batch().take_pass1(sp, r, i < nb)) continue;
        pairs[pos++] = sp;
    }
    k->getScores8(pairs, b.ref.data(), b.qer.data(), aln, pos, 1, 1);
    if (nb || p1_any) rescue_band_batch().run_pass1(b.ref.data(), b.qer.data(), aln, k);
}

static std::unique_ptr<Ikswv> make_k(int maxr, int maxq)
{
    return make_kswv(6, 1, 6, 1, 1, -4, 1, maxr + 64, maxq + 64);
}

/* ---------------------------------------------------------------------------------------- */
static int run_eq(std::vector<Job> &jobs, int minsc, int max_hits, int scalar_stride)
{
    long mism = 0, n = 0, npass = 0, nfull = 0, nb1 = 0, nb2 = 0, nband = 0, scal_mm = 0, te2_diff = 0;
    long nview = 0, view_mm = 0, r2 = 0, ccap = 0;
    std::unique_ptr<rescue_prune_scratch> sscratch(new rescue_prune_scratch());
    std::vector<rb_comp> cn, cs;
    long p1_band = 0, p1_guard = 0, band_hi = 0;
    long cls_n[128] = {0}, cls_mm[128] = {0};
    const bool negative_control = rescue_env_opt_in("RB_NEGATIVE_CONTROL");
    bool injected = false;
    const int B = 2048;
    for (size_t base = 0; base < jobs.size(); base += B) {
        const int m = (int)std::min<size_t>(B, jobs.size() - base);
        int maxr = 0, maxq = 0;
        for (int t = 0; t < m; t++) {
            maxr = std::max(maxr, (int)jobs[base + t].ref.size());
            maxq = std::max(maxq, (int)jobs[base + t].q.size() + 16);
        }
        auto k = make_k(maxr, maxq);
        // --- truth: full windows ---
        Batch T;
        for (int t = 0; t < m; t++) {
            const Job &J = jobs[base + t];
            T.add(J.ref.data(), (int)J.ref.size(), J.q.data(), (int)J.q.size(), minsc);
        }
        T.finalize();
        run_batch(k.get(), T, 0, false);
        // --- production pruned + banded pipeline ---
        Batch P;
        std::vector<int> kind(m), off(m, 0), idx(m, -1);
        for (int t = 0; t < m; t++) {
            const Job &J = jobs[base + t];
            const int len1 = (int)J.ref.size(), len2 = (int)J.q.size();
            int hb, he;
            rescue_prune_view view;
            kind[t] = rescue_prune_window(J.ref.data(), len1, J.q.data(), len2, minsc, max_hits, &hb, &he, &view);
            if (kind[t] == RESCUE_PRUNE_B2 && view.bnd16) {
                /* The SIMD and scalar filter views of one job must give plan() the same components
                 * (the scalar view otherwise only serves windows past the SIMD filter's capacity). */
                int shb, she;
                const int sk = rescue_prune_window_scalar(J.ref.data(), len1, J.q.data(), len2, minsc, max_hits,
                                                          *sscratch, &shb, &she);
                bool vok = sk == RESCUE_PRUNE_B2 && shb == hb && she == he;
                if (vok) {
                    const rescue_prune_view sv = rescue_prune_scalar_view(*sscratch, len1, len2, minsc);
                    /* The whole view, then sub-ranges like the ones plan() passes (one component's
                     * diagonals at a higher threshold): an unaligned start, a start just below a
                     * 64-diagonal word boundary, a range shorter than a word, and a random one. */
                    const int nd = view.nd;
                    const int w0 = std::max(0, (nd / 2 & ~63) - 1 - rnd(3));
                    const int ra = rnd(nd), rb = ra + 1 + rnd(nd - ra);
                    const int ranges[5][2] = {{0, nd},
                                              {std::min(nd - 1, 1 + rnd(7)), nd},
                                              {w0, std::min(nd, w0 + 70)},
                                              {ra, std::min(rb, ra + 1 + rnd(63))},
                                              {ra, rb}};
                    for (const int tau : {minsc, minsc + 5, 30, 50, 80})
                        for (const auto &r : ranges) {
                            if (!vok || r[0] >= r[1]) continue;
                            cn.clear(); cs.clear();
                            const bool fn = rescue_band_components(view, tau, r[0], r[1], cn, 1 << 20);
                            const bool fs = rescue_band_components(sv, tau, r[0], r[1], cs, 1 << 20);
                            vok = fn && fs && cn.size() == cs.size();
                            for (size_t c = 0; vok && c < cn.size(); c++)
                                vok = cn[c].ub == cs[c].ub && cn[c].i0 == cs[c].i0 && cn[c].dlo == cs[c].dlo
                                      && cn[c].dhi == cs[c].dhi && cn[c].dmaxhit == cs[c].dmaxhit;
                        }
                }
                nview++;
                if (!vok) {
                    view_mm++;
                    if (view_mm < 5) fprintf(stderr, "VIEW MISMATCH job=%zu cls=%d: NEON and scalar views disagree\n", base + t, J.cls);
                }
            }
            if (kind[t] == RESCUE_PRUNE_B1) { nb1++; continue; }
            if (kind[t] == RESCUE_PRUNE_B2) {
                nb2++;
                if (rescue_band_enabled())
                    rescue_band_batch().plan(view, len1, len2, hb, he, minsc);
                off[t] = hb;
                idx[t] = P.add(J.ref.data() + hb, he - hb + 1, J.q.data(), len2, minsc);
            } else {
                nfull++;
                idx[t] = P.add(J.ref.data(), len1, J.q.data(), len2, minsc);
            }
            if (rescue_band_enabled()) rescue_band_batch().commit(idx[t]);
        }
        const int pn = P.n;
        P.finalize();
        std::stable_sort(P.sp.begin(), P.sp.begin() + pn, [](const SeqPair &x, const SeqPair &y) { return x.len1 < y.len1; });
        const int nbd = rescue_band_enabled() ? rescue_band_batch().partition(P.sp.data(), pn) : 0;
        nband += nbd;
        const uint64_t pb0 = rescue_band_batch().stats().p1_band, pg0 = rescue_band_batch().stats().p1_guard;
        run_batch(k.get(), P, nbd, true);
        std::vector<char> banded(m, 0);
        for (int t = 0; t < m; t++) banded[t] = idx[t] >= 0 && rescue_band_batch().banded(idx[t]);
        r2 += (long)rescue_band_batch().stats().r2_band;
        ccap += (long)rescue_band_batch().stats().comp_cap;
        p1_band += (long)(rescue_band_batch().stats().p1_band - pb0);
        p1_guard += (long)(rescue_band_batch().stats().p1_guard - pg0);
        rescue_band_batch().reset();
        // --- compare ---
        for (int t = 0; t < m; t++) {
            const Job &J = jobs[base + t];
            const kswr_t &a = T.aln[t];
            n++; cls_n[J.cls]++;
            band_hi += banded[t] && a.score >= 241;
            bool ok;
            const bool pass = a.score >= minsc;
            npass += pass;
            if (kind[t] == RESCUE_PRUNE_B1) ok = !pass;
            else {
                kswr_t b = P.aln[idx[t]];
                if (negative_control && !injected && pass && banded[t]) { b.te++; injected = true; }
                if (!pass) ok = b.score < minsc;
                else {
                    ok = a.score == b.score && a.te == b.te + off[t] && a.qe == b.qe && a.score2 == b.score2
                         && a.tb == b.tb + off[t] && a.qb == b.qb;
                    if (ok && a.score2 >= 0 && a.te2 != b.te2 + off[t]) te2_diff++;
                }
                if (!ok && mism < 10) {
                    fprintf(stderr, "MISMATCH job=%zu cls=%d kind=%d len1=%zu len2=%zu off=%d truth{s=%d te=%d qe=%d s2=%d tb=%d qb=%d} got{s=%d te=%d qe=%d s2=%d tb=%d qb=%d} banded=%d\n",
                            base + t, J.cls, kind[t], J.ref.size(), J.q.size(), off[t], a.score, a.te, a.qe, a.score2, a.tb, a.qb,
                            b.score, b.te + off[t], b.qe, b.score2, b.tb + off[t], b.qb, (int)banded[t]);
                }
            }
            if (kind[t] == RESCUE_PRUNE_B1 && !ok && mism < 10)
                fprintf(stderr, "MISMATCH job=%zu cls=%d B1 but truth score=%d\n", base + t, J.cls, a.score);
            if (J.cls == 3 && kind[t] != RESCUE_PRUNE_FULL) { ok = false; if (mism < 10) fprintf(stderr, "N job not FULL\n"); }
            if (!ok) { mism++; cls_mm[J.cls]++; }
            if (scalar_stride > 0 && (n % scalar_stride) == 0) {
                int S, te, qe; std::vector<int> R;
                scalar_dp(J, S, te, qe, R);
                const int s2 = S >= minsc ? scalar_score2(R, S, te, minsc) : -1;
                bool sok = S == a.score;
                if (S > 0) sok = sok && te == a.te && qe == a.qe;
                if (S >= minsc) {
                    int tb, qb;
                    scalar_pass1(J, S, te, qe, tb, qb);
                    sok = sok && s2 == a.score2 && tb == a.tb && qb == a.qb;
                }
                if (!sok) {
                    scal_mm++;
                    if (scal_mm < 5) fprintf(stderr, "SCALAR vs kswv job=%zu cls=%d: scalar{%d %d %d %d} kswv{%d %d %d %d}\n",
                                             base + t, J.cls, S, te, qe, s2, a.score, a.te, a.qe, a.score2);
                }
            }
        }
    }
    printf("eq: jobs=%ld pass=%ld full=%ld b1=%ld b2=%ld banded_parents=%ld banded_s241=%ld pass1_banded=%ld "
           "pass1_guard=%ld round2=%ld comp_cap=%ld views_checked=%ld MISMATCHES=%ld view_mismatch=%ld "
           "scalar_vs_kswv_mismatch=%ld (te2 differs, unconsumed: %ld)\n",
           n, npass, nfull, nb1, nb2, nband, band_hi, p1_band, p1_guard, r2, ccap, nview, mism, view_mm, scal_mm,
           te2_diff);
    for (int c = 0; c < 128; c++)
        if (cls_n[c]) printf("  class %3d: jobs=%ld mismatches=%ld\n", c, cls_n[c], cls_mm[c]);
    /* A pass with nothing banded compares kswv with kswv: fail it where banding runs, so the gate
     * cannot go vacuous (a cost gate or generator change that stops banding). */
    const bool vacuous = rescue_band_enabled() && (nband == 0 || p1_band == 0);
    if (vacuous) fprintf(stderr, "FAIL: no job was banded in pass 0 or pass 1 (open the cost gate: BWA3_RESCUE_BAND_COST=100000000)\n");
    /* The 8-bit ceiling class must reach the top of the band's H range, or the check of the x86
     * cell's biased add there (exact only while H <= 250) went vacuous. */
    const bool ceiling_vacuous = rescue_band_enabled() && cls_n[12] > 0 && band_hi == 0;
    if (ceiling_vacuous) fprintf(stderr, "FAIL: no banded job scored 241 or more (the 8-bit ceiling class stopped reaching it)\n");
    /* Production falls back to kswv when a banded pass-1 max is not S (the guard in rescue_band.h),
     * but that cannot happen, so here any fallback fails the run: it means the band argument broke. */
    if (p1_guard) fprintf(stderr, "FAIL: %ld pass-1 guard fallbacks (banded max != S)\n", p1_guard);
    return mism || scal_mm || view_mm || vacuous || ceiling_vacuous || p1_guard ? 1 : 0;
}

/* ---------------------------------------------------------------------------------------- */
/* timing: single thread, per batch of 2048 dumped jobs, on the SAME parents:                 */
/*   hull  = filter + kswv phase 0 on every non-B1 job (hull or full window) [BAND off]        */
/*   band  = filter + plan + kswv phase 0 on the non-banded + run_pass0 on the banded         */
/* ---------------------------------------------------------------------------------------- */
/* Lanes per 8-bit kswv group at the tier make_kswv dispatched to (its SIMD_WIDTH8). */
static int kswv_lanes8()
{
    bwamem3_simd_init();
    switch (bwamem3_simd_tier()) {
    case BWAMEM3_TIER_AVX512BW: return 64;
    case BWAMEM3_TIER_NEON: return 16;
    default: return 32;   // AVX2, the x86 floor of the batched kswv kernels
    }
}

static void run_time(std::vector<Job> &jobs, int minsc, int max_hits, int reps)
{
    using clk = std::chrono::steady_clock;
    const int kl = kswv_lanes8();
    const int B = 2048;
    double best[2][4];   // [mode][filter, plan, kswv, band]
    double kcells[2] = {0, 0};
    for (auto &x : best) for (double &y : x) y = 1e30;
    for (int rep = 0; rep < reps; rep++) {
        for (int mode = 0; mode < 2; mode++) {
            double acc[4] = {0, 0, 0, 0};
            for (size_t base = 0; base < jobs.size(); base += B) {
                const int m = (int)std::min<size_t>(B, jobs.size() - base);
                int maxr = 0, maxq = 0;
                for (int t = 0; t < m; t++) {
                    maxr = std::max(maxr, (int)jobs[base + t].ref.size());
                    maxq = std::max(maxq, (int)jobs[base + t].q.size() + 16);
                }
                static std::unique_ptr<Ikswv> k;
                static int kr = 0, kq = 0;
                if (!k || maxr > kr || maxq > kq) { kr = std::max(kr, maxr); kq = std::max(kq, maxq); k = make_k(kr, kq); }
                Batch P;
                P.ref.reserve(3000000); P.qer.reserve(400000);
                double tf = 0, tp = 0;
                for (int t = 0; t < m; t++) {
                    const Job &J = jobs[base + t];
                    const int len1 = (int)J.ref.size(), len2 = (int)J.q.size();
                    int hb, he;
                    rescue_prune_view view;
                    auto a = clk::now();
                    const int kd = rescue_prune_window(J.ref.data(), len1, J.q.data(), len2, minsc, max_hits, &hb, &he, &view);
                    auto b = clk::now();
                    tf += std::chrono::duration<double>(b - a).count();
                    if (kd == RESCUE_PRUNE_B1) continue;
                    if (kd == RESCUE_PRUNE_B2 && mode == 1) {
                        a = clk::now();
                        rescue_band_batch().plan(view, len1, len2, hb, he, minsc);
                        tp += std::chrono::duration<double>(clk::now() - a).count();
                    }
                    int id;
                    if (kd == RESCUE_PRUNE_B2) id = P.add(J.ref.data() + hb, he - hb + 1, J.q.data(), len2, minsc);
                    else id = P.add(J.ref.data(), len1, J.q.data(), len2, minsc);
                    rescue_band_batch().commit(id);
                }
                const int pn = P.n;
                P.finalize();
                std::stable_sort(P.sp.begin(), P.sp.begin() + pn, [](const SeqPair &x, const SeqPair &y) { return x.len1 < y.len1; });
                const int nbd = mode == 1 ? rescue_band_batch().partition(P.sp.data(), pn) : 0;
                for (int g = nbd; g < pn; g += kl) {   // kswv padded cells: kl x maxLen1 x maxQuanta
                    int ml = 0, mq = 0;
                    for (int l = g; l < std::min(pn, g + kl); l++) { ml = std::max(ml, (int)P.sp[l].len1); mq = std::max(mq, kswv_query_quantum8((int)P.sp[l].len2)); }
                    kcells[mode] += (double)kl * ml * mq;
                }
                auto a = clk::now();
                k->getScores8(P.sp.data() + nbd, P.ref.data(), P.qer.data(), P.aln.data(), pn - nbd, 1, 0);
                auto b = clk::now();
                if (nbd) rescue_band_batch().run_pass0(P.sp.data(), nbd, P.ref.data(), P.qer.data(), P.aln.data(), k.get());
                auto c = clk::now();
                rescue_band_batch().reset();
                acc[0] += tf; acc[1] += tp;
                acc[2] += std::chrono::duration<double>(b - a).count();
                acc[3] += std::chrono::duration<double>(c - b).count();
            }
            for (int x = 0; x < 4; x++) best[mode][x] = std::min(best[mode][x], acc[x]);
        }
    }
    const double nj = (double)jobs.size();
    for (int mode = 0; mode < 2; mode++) {
        const double tot = best[mode][0] + best[mode][1] + best[mode][2] + best[mode][3];
        printf("%-5s kswv padded cells/job %.0f (%.3f ns/cell)\n", mode ? "band" : "hull", kcells[mode] / reps / nj,
               1e9 * best[mode][2] / (kcells[mode] / reps));
        printf("%-5s us/job: filter %.3f plan %.3f kswv_p0 %.3f band_p0 %.3f total %.3f\n", mode ? "band" : "hull",
               1e6 * best[mode][0] / nj, 1e6 * best[mode][1] / nj, 1e6 * best[mode][2] / nj, 1e6 * best[mode][3] / nj,
               1e6 * tot / nj);
    }
}

int main(int argc, char **argv)
{
    if (argc < 2) {
        fprintf(stderr, "usage: %s eq <ngen> <seed> [dumps...] | time <reps> <stride> <dumps...>\n", argv[0]);
        return 2;
    }
    /* Without banding every mode compares kswv with itself (a build without a SIMD band kernel --
     * x86 below arch=avx2 -- or banding turned off in the environment). */
    if (!rescue_band_enabled()) {
        fprintf(stderr, "FAIL: banded rescue is off (no NEON / AVX2 band kernel in this build, or BWA3_RESCUE_BAND "
                        "turned it off), so there is nothing to compare\n");
        return 1;
    }
    /* The kswv tier is the oracle here and on x86 follows the host (or BWAMEM3_FORCE_TIER), so name
     * it: the same run checks the band kernel against a different kswv on another runner. */
    bwamem3_simd_init();
    printf("kswv tier: %s\n", bwamem3_simd_tier_name(bwamem3_simd_tier()));
    const int minsc = rescue_env_int("RB_MINSC", MINSC_DEFAULT);
    /* The production default (rescue_prune_max_hits in bwamem_pair.cpp). */
    const int max_hits = rescue_env_int("BWA3_RESCUE_PRUNE_MAX_HITS",
                                        rescue_prune_max_hits_default(rescue_band_enabled()));
    std::vector<Job> jobs;
    if (argc < 4) {
        fprintf(stderr, "usage: %s eq <ngen> <seed> [dumps...] | time <reps> <stride> <dumps...>\n", argv[0]);
        return 2;
    }
    if (!strcmp(argv[1], "eq")) {
        const int ngen = atoi(argv[2]);
        rng.seed(strtoull(argv[3], nullptr, 10));
        for (int t = 0; t < ngen; t++) jobs.push_back(gen(t % 13));
        for (int a = 4; a < argc; a++) load_dump(argv[a], jobs, 1);
        const int ss = rescue_env_int("RB_SCALAR_STRIDE", 0);
        return run_eq(jobs, minsc, max_hits, ss);
    }
    if (!strcmp(argv[1], "time")) {
        const int reps = atoi(argv[2]), stride = atoi(argv[3]);
        for (int a = 4; a < argc; a++) load_dump(argv[a], jobs, stride);
        run_time(jobs, minsc, max_hits, reps);
        return 0;
    }
    return 2;
}
