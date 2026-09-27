/* Exactness + timing harness for the banded mate rescue (src/rescue_band.{h,cpp}).
 * Build: `make rescue-band-harness` (arch=arm64 on aarch64, arch=avx2 on x86).
 *
 *   rescue_band_harness eq   <gen-jobs> <seed> [dump files...]   exactness (exit 1 on any mismatch
 *                                                                or pass-1 guard fallback)
 *   rescue_band_harness time <reps> <stride> <dump files...>     single-thread pass-0 timing
 *
 * eq: for every job, the PRODUCTION pipeline -- rescue_prune_window, RescueBandBatch::plan/commit,
 * the length sort, partition, kswv phase 0 on the non-banded pairs, run_pass0 on the banded ones,
 * then the mem_sam_pe_batch_run post-processing and kswv phase 1 -- is compared against the
 * pre-pruning pipeline (full window, kswv phase 0 + phase 1), field by field: score (pass/fail
 * when the truth fails), te, qe, score2, tb, qb (offsets applied). The production pipeline includes
 * the banded pass 1 (take_pass1 / run_pass1) of the banded parents and, under the default
 * BWA3_RESCUE_BAND_P1=2, of every other pair (hull or full window, kswv pass 0), so tb / qb check
 * it against kswv phase 1 on the full window; the summary counts how many jobs took it. A pass-1
 * job whose banded max is not S falls back to kswv in production (the guard in rescue_band.h), but
 * that cannot happen, so here any such fallback fails the run: it means the band argument broke.
 * A scalar full-window DP with kswv's score2 semantics (optsim.cpp dp/score2_of) plus a scalar
 * reversed-prefix DP for (tb, qb) cross-check the truth on every --scalar-stride-th job.
 * Generated classes: random windows with mutated mate copies, edge copies around te +- S (score2
 * parity/zeroing), tiny windows, N bases (must be FULL), tandem repeats (hit-dense, many
 * components), in-zone secondary copies (round-2 triggers), ragged len2, and pass-1 adversaries
 * (class 9: gaps at the band-edge bound, start ties, qe next to the pad columns, te at the
 * window edges). Env knobs are the production ones (src/rescue_env.h); the caller sets e.g.
 * BWA3_RESCUE_BAND_COST=100000000 (band every B2 parent) and BWA3_RESCUE_PRUNE_MAX_HITS. RB_MINSC
 * (the min seed length, default 19; minsc = RB_MINSC * a), RB_SCALAR_STRIDE (default 0: no scalar
 * cross-check; default scoring only) and RB_SCORING (unset: the default scoring; "a,b,o_del,e_del,
 * o_ins,e_ins": fixed; "random": a fresh draw per 2048-job batch, see draw_params) and RB_METH
 * (genomic | neutral | collapsed | random: --meth matrices, see meth_matrix) with RB_METH_REL
 * (conv, the default: the filter matches converted copies; expand: genomic and neutral batches use
 * the relation-expanded filter, rescue_prune_params::set_meth_rel, collapsed stays on conv), RB_16BIT (1: every
 * job runs through kswv's 16-bit kernels, uncut and unbanded, with longer mates) and RB_KMAX (the
 * largest filter K, default 5; 0: a draw in [5, rescue_prune_scratch_kmax] per batch) are harness-only.
 * Banding runs only at the default scoring, as in production.
 *
 * Dump files: one record per prune-eligible rescue job, as the aligner saw it (full window and
 * oriented mate): int32 len1, int32 len2, len1 reference bytes, len2 query bytes (2-bit codes, 4 = N).
 * The aligner does not write them; capture new ones by instrumenting mem_matesw_batch_pre where it
 * calls rescue_prune_window. */
#include "rescue_band.h"
#include "rescue_env.h"

#include <algorithm>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <random>
#include <string>
#include <vector>

static const int MINSC_DEFAULT = 19;

struct Job { std::vector<uint8_t> ref, q; int cls; };
struct Res { int score, te, qe, score2, tb, qb; };

/* ---------------------------------------------------------------------------------------- */
/* scalar reference (optsim.cpp)                                                              */
/* ---------------------------------------------------------------------------------------- */
static void scalar_dp(const Job &jb, int &S, int &te, int &qe, std::vector<int> &R)
{
    const int len1 = (int)jb.ref.size(), len2 = (int)jb.q.size(), quanta = (len2 + 15) / 16 * 16;
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
static bool g_16 = false;   // RB_16BIT
static int rnd(int n) { return n <= 1 ? 0 : (int)(rng() % (uint64_t)n); }
static int rndr(int a, int b) { return a + rnd(b - a + 1); }
static double unif() { return (rng() >> 11) * (1.0 / 9007199254740992.0); }

static std::vector<uint8_t> mutate(const std::vector<uint8_t> &s, int a, int b, double mm, double ind)
{
    std::vector<uint8_t> o;
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
    if (g_16 && rnd(2)) len2 = rnd(2) ? 250 : rndr(17, 500);   // 2x250 and long mates
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
    case 1: {   // edge fuzz: a full copy + partial copies ending near te-S, te+S (edgefuzz.cpp)
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
            // a query PREFIX ending just after te: in the zone, but diagonal ~ len2 - sl away
            if (rnd(4)) { q0 = rnd(6); end = te + rndr(5, len2 / 3); }
            else { q0 = len2 - sl - rnd(10); end = te - S + 1 + rnd(40) + sl; }
            q0 = std::max(0, std::min(q0, len2 - sl));
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
    case 9: {   // pass-1 adversaries (banded pass 1, rescue_band.h)
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
    while (fread(h, sizeof h, 1, fp) == 1) {
        Job J; J.cls = 100;
        J.ref.resize(h[0]); J.q.resize(h[1]);
        if (fread(J.ref.data(), 1, h[0], fp) != (size_t)h[0] || fread(J.q.data(), 1, h[1], fp) != (size_t)h[1]) break;
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
        p.h0 = KSW_XSUBO | KSW_XSTART | (g_16 ? 0 : KSW_XBYTE) | minsc;
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
    if (g_16) {   // RB_16BIT: kswv's 16-bit kernels for both passes, no banding (as production)
        k->getScores16(pairs, b.ref.data(), b.qer.data(), aln, b.n, 1, 0);
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
            pairs[pos++] = sp;
        }
        k->getScores16(pairs, b.ref.data(), b.qer.data(), aln, pos, 1, 1);
        return;
    }
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
        if ((i < nb || p1_any) && rescue_band_batch().take_pass1(sp, r, i < nb, b.ref.data(), b.qer.data())) continue;
        pairs[pos++] = sp;
    }
    k->getScores8(pairs, b.ref.data(), b.qer.data(), aln, pos, 1, 1);
    if (nb || p1_any) rescue_band_batch().run_pass1(b.ref.data(), b.qer.data(), aln, k);
}

static std::unique_ptr<Ikswv> make_k(int maxr, int maxq, const rescue_prune_params &p,
                                     const int8_t *mat25 = nullptr)
{
    return make_kswv(p.o_del, p.e_del, p.o_ins, p.e_ins, (int8_t)p.a, (int8_t)-p.b, 1, maxr + 64, maxq + 64,
                     mat25);
}

/* ---------------------------------------------------------------------------------------- */
/* --meth: RB_METH = genomic | neutral | collapsed | random (per batch). A batch takes one       */
/* hypothesis (hyp 1 = OT: reference C / read T freed; hyp 0 = OB: reference G / read A), its mates */
/* are converted (C->T or G->A) at a per-batch rate, and the truth is the mat-aware kswv built      */
/* with the same matrix as mem_opt_fill_meth_mat.                                                   */
/* ---------------------------------------------------------------------------------------- */
enum { METH_NONE = 0, METH_GENOMIC, METH_NEUTRAL, METH_COLLAPSED, METH_RANDOM };
static int g_meth = METH_NONE;
static bool g_meth_rel = false;   // RB_METH_REL=expand

static void meth_matrix(const rescue_prune_params &p, int kind, int hyp, int8_t mat[25])
{
    for (int i = 0; i < 5; i++)
        for (int j = 0; j < 5; j++)
            mat[i * 5 + j] = (int8_t)(i == 4 || j == 4 ? -1 : i == j ? p.a : -p.b);
    const int fr = hyp ? 1 : 2, fq = hyp ? 3 : 0;
    mat[fr * 5 + fq] = (int8_t)(kind == METH_NEUTRAL ? 0 : p.a);
    if (kind == METH_COLLAPSED) mat[fq * 5 + fr] = (int8_t)p.a;
}

/* ---------------------------------------------------------------------------------------- */
/* scoring: RB_SCORING unset = the default; "a,b,o_del,e_del,o_ins,e_ins" = fixed; "random" = a  */
/* fresh draw per batch (general-scoring plan, section 9), including degenerate parameter sets    */
/* (which must be refused: every job FULL) and the V4 boundary minsc = a (K - 1) + 1 / + 2.        */
/* ---------------------------------------------------------------------------------------- */
struct Scoring { int a = 1, b = 4, o_del = 6, e_del = 1, o_ins = 6, e_ins = 1; };
static bool g_random_scoring = false;
static Scoring g_fixed_scoring;
/* RB_KMAX: the largest filter K (rescue_prune_params::from's k_max; production BWA3_RESCUE_PRUNE_KMAX):
 * 5 (default) to rescue_prune_scratch_kmax, or 0 = a fresh draw in [5, kmax] per batch. */
static int g_kmax = 5;
static int draw_kmax() { return g_kmax > 0 ? g_kmax : rndr(5, rescue_prune_scratch_kmax); }

static rescue_prune_params draw_params(int minsc_default)
{
    Scoring sc = g_fixed_scoring;
    int k = minsc_default;   // min_seed_len; minsc = k * a
    bool fixed_minsc = true;
    if (g_random_scoring) {
        const int mode = rnd(20);
        if (mode == 0) {   // the defaults
            sc = Scoring();
        } else if (mode == 1) {   // degenerate: must be refused before any division
            sc = Scoring();
            switch (rnd(6)) {
                case 0: sc.e_del = 0; break;
                case 1: sc.e_ins = 0; break;
                case 2: sc.e_del = sc.e_ins = 0; break;
                case 3: sc.b = 0; break;
                case 4: sc.a = 0; break;
                default: sc.b = 129; break;
            }
        } else if (mode == 2) {   // one gap type cheaper than (K - 1) a (V2), either side
            sc = Scoring();
            if (rnd(2)) { sc.o_ins = 2; sc.e_ins = 1; } else { sc.o_del = 2; sc.e_del = 1; }
        } else if (mode == 3) {   // split and shared fused cells: -O 6 -E 1,2 / -O 6,7 -E 2,1
            sc = Scoring();
            if (rnd(2)) { sc.e_ins = 2; } else { sc.o_ins = 7; sc.e_del = 2; }
        } else {
            sc.a = 1 + rnd(3);
            sc.b = rndr(1, 12);
            sc.o_del = rndr(0, 20); sc.o_ins = rndr(0, 20);
            sc.e_del = rndr(1, 6); sc.e_ins = rndr(1, 6);
        }
        k = rndr(5, 40);
        fixed_minsc = false;
    }
    int minsc = k * sc.a;
    const int kmax = draw_kmax();
    rescue_prune_params p = rescue_prune_params::from(sc.a, sc.b, sc.o_del, sc.e_del, sc.o_ins, sc.e_ins,
                                                      minsc, 3, kmax);
    if (!fixed_minsc && p.K >= 3 && rnd(8) == 0) {   // the V4 boundary
        minsc = sc.a * (p.K - 1) + 1 + rnd(2);
        p = rescue_prune_params::from(sc.a, sc.b, sc.o_del, sc.e_del, sc.o_ins, sc.e_ins, minsc, 3, kmax);
    }
    if (sc.a < 1 || sc.b < 1 || sc.b > 128 || sc.e_del < 1 || sc.e_ins < 1) {
        /* Degenerate: from() must refuse it (checked: every job FULL); the truth runs with the
         * default scoring, since kswv itself cannot take these values. */
        rescue_prune_params d = rescue_prune_params::defaults(minsc);
        d.valid = p.valid;
        return d;
    }
    /* Keep the raw scoring even when refused: the truth still runs with it. */
    p.a = sc.a; p.b = sc.b; p.o_del = sc.o_del; p.e_del = sc.e_del; p.o_ins = sc.o_ins; p.e_ins = sc.e_ins;
    p.minsc = minsc;
    return p;
}

/* matesw_use_u8 (bwamem_pair.cpp): the job runs in the 8-bit kernels. */
static bool use_u8(int len2, int a, int shift)
{
    const int prod = len2 * a;
    return prod + shift <= 254 && (prod < 250 || len2 % 16 == 0 || len2 % 16 >= 9);
}

/* ---------------------------------------------------------------------------------------- */
static int run_eq(std::vector<Job> &all_jobs, int minsc_default, int max_hits, int scalar_stride)
{
    long mism = 0, n = 0, npass = 0, nfull = 0, nb1 = 0, nb2 = 0, nband = 0, scal_mm = 0, te2_diff = 0;
    long p1_band = 0, p1_guard = 0, nbatch = 0, nvalid = 0, ntrunc = 0, nk[rescue_prune_scratch_kmax + 1] = {0}, nmeth[4] = {0};
    long cls_n[128] = {0}, cls_mm[128] = {0};
    const int B = 2048;
    for (size_t base = 0; base < all_jobs.size(); base += B) {
        const int m = (int)std::min<size_t>(B, all_jobs.size() - base);
        rescue_prune_params pp = draw_params(minsc_default);
        int meth = g_meth, hyp = rnd(2);
        if (meth == METH_RANDOM) meth = METH_GENOMIC + rnd(3);
        int8_t mat[25];
        if (meth) {
            meth_matrix(pp, meth, hyp, mat);
            if (g_meth_rel && meth != METH_COLLAPSED) pp.set_meth_rel(hyp);
            else pp.set_meth(hyp);
            nmeth[meth]++;
        }
        const int minsc = pp.minsc;
        /* As production: pass-0 banding for any scoring the filter accepts, pass-1 banding for any
         * scoring the kernels take (rb_scoring::valid); --meth with the batch's matrix. */
        rb_scoring bsc;
        bsc.a = pp.a; bsc.b = pp.b; bsc.o_del = pp.o_del; bsc.e_del = pp.e_del; bsc.o_ins = pp.o_ins; bsc.e_ins = pp.e_ins;
        if (meth) bsc.set_matrix(mat);
        const bool band = !g_16 && rescue_band_enabled() && pp.valid;
        const bool p1_any = !g_16 && rescue_band_enabled() && bsc.valid();
        if (rescue_band_enabled()) rescue_band_batch().set_scoring(bsc);
        nbatch++; nvalid += pp.valid; if (pp.valid) nk[pp.K]++;
        /* This batch's jobs: mates too long for the 8-bit kernels at this a are cut to a random
         * admissible slice (the planted copies still match it). */
        const int shift = std::max(1, pp.b);
        std::vector<Job> jobs(all_jobs.begin() + base, all_jobs.begin() + base + m);
        for (Job &J : jobs) {
            const int len2 = (int)J.q.size();
            if (g_16 || use_u8(len2, pp.a, shift)) continue;
            int nl = std::min(len2, std::min((254 - shift) / pp.a, 249 / pp.a));
            if (nl < 5) nl = 0;
            if (nl > 17) nl = rndr(17, nl);
            const int s0 = rnd(len2 - nl + 1);
            J.q = std::vector<uint8_t>(J.q.begin() + s0, J.q.begin() + s0 + nl);
            ntrunc++;
        }
        if (meth) {   // convert the mates: the planted copies keep the reference bases
            /* TAPS-like low rates (mostly unconverted reads) are the relation filter's case */
            const double rate = g_meth_rel ? (rnd(2) ? (rnd(2) ? 0.02 : 0.05) : rnd(2) ? 0.5 : 0.95)
                                           : rnd(3) == 0 ? 0.05 : rnd(2) ? 0.5 : 0.95;
            const uint8_t from = hyp ? 1 : 2, to = hyp ? 3 : 0;
            for (Job &J : jobs)
                for (uint8_t &x : J.q)
                    if (x == from && unif() < rate) x = to;
        }
        int maxr = 0, maxq = 0;
        for (int t = 0; t < m; t++) {
            maxr = std::max(maxr, (int)jobs[t].ref.size());
            maxq = std::max(maxq, (int)jobs[t].q.size() + 16);
        }
        auto k = make_k(maxr, maxq, pp, meth ? mat : nullptr);
        if (k->needsScalar()) { fprintf(stderr, "meth matrix needs the scalar kernel\n"); exit(2); }
        // --- truth: full windows ---
        Batch T;
        for (int t = 0; t < m; t++) {
            const Job &J = jobs[t];
            T.add(J.ref.data(), (int)J.ref.size(), J.q.data(), (int)J.q.size(), minsc);
        }
        T.finalize();
        run_batch(k.get(), T, 0, false);
        // --- production pruned + banded pipeline ---
        Batch P;
        std::vector<int> kind(m), off(m, 0), idx(m, -1);
        for (int t = 0; t < m; t++) {
            const Job &J = jobs[t];
            const int len1 = (int)J.ref.size(), len2 = (int)J.q.size();
            int hb, he;
            kind[t] = rescue_prune_window(J.ref.data(), len1, J.q.data(), len2, pp, max_hits, &hb, &he);
            if (kind[t] == RESCUE_PRUNE_B1) { nb1++; continue; }
            if (kind[t] == RESCUE_PRUNE_B2) {
                nb2++;
                if (band)
                    rescue_band_batch().plan(rescue_prune_last_view(), pp, len1, len2, hb, he);
                off[t] = hb;
                idx[t] = P.add(J.ref.data() + hb, he - hb + 1, J.q.data(), len2, minsc);
            } else {
                nfull++;
                idx[t] = P.add(J.ref.data(), len1, J.q.data(), len2, minsc);
            }
            if (band) rescue_band_batch().commit(idx[t]);
        }
        const int pn = P.n;
        P.finalize();
        std::stable_sort(P.sp.begin(), P.sp.begin() + pn, [](const SeqPair &x, const SeqPair &y) { return x.len1 < y.len1; });
        const int nbd = band ? rescue_band_batch().partition(P.sp.data(), pn) : 0;
        nband += nbd;
        const uint64_t pb0 = rescue_band_batch().stats.p1_band, pg0 = rescue_band_batch().stats.p1_guard;
        run_batch(k.get(), P, nbd, p1_any);
        p1_band += (long)(rescue_band_batch().stats.p1_band - pb0);
        p1_guard += (long)(rescue_band_batch().stats.p1_guard - pg0);
        rescue_band_batch().reset();
        // --- compare ---
        for (int t = 0; t < m; t++) {
            const Job &J = jobs[t];
            const kswr_t &a = T.aln[t];
            n++; cls_n[J.cls]++;
            bool ok;
            const bool pass = a.score >= minsc;
            npass += pass;
            if (kind[t] == RESCUE_PRUNE_B1) ok = !pass;
            else {
                const kswr_t &b = P.aln[idx[t]];
                if (!pass) ok = b.score < minsc;
                else {
                    ok = a.score == b.score && a.te == b.te + off[t] && a.qe == b.qe && a.score2 == b.score2
                         && a.tb == b.tb + off[t] && a.qb == b.qb;
                    if (ok && a.score2 >= 0 && a.te2 != b.te2 + off[t]) te2_diff++;
                }
                if (!ok && mism < 10) {
                    fprintf(stderr, "scoring a=%d b=%d o=%d,%d e=%d,%d minsc=%d K=%d c=%d valid=%d\n", pp.a, pp.b,
                            pp.o_del, pp.o_ins, pp.e_del, pp.e_ins, minsc, pp.K, pp.c, pp.valid);
                    fprintf(stderr, "MISMATCH job=%zu cls=%d kind=%d len1=%zu len2=%zu off=%d truth{s=%d te=%d qe=%d s2=%d tb=%d qb=%d} got{s=%d te=%d qe=%d s2=%d tb=%d qb=%d} banded=%d\n",
                            base + t, J.cls, kind[t], J.ref.size(), J.q.size(), off[t], a.score, a.te, a.qe, a.score2, a.tb, a.qb,
                            b.score, b.te + off[t], b.qe, b.score2, b.tb + off[t], b.qb, 0);
                }
            }
            if (kind[t] == RESCUE_PRUNE_B1 && !ok && mism < 10)
                fprintf(stderr, "MISMATCH job=%zu cls=%d B1 but truth score=%d (a=%d b=%d o=%d,%d e=%d,%d minsc=%d K=%d c=%d)\n",
                        base + t, J.cls, a.score, pp.a, pp.b, pp.o_del, pp.o_ins, pp.e_del, pp.e_ins, minsc, pp.K, pp.c);
            /* An N in the window or the mate must give FULL (class 3 plants them; a mate cut for
             * the 8-bit kernels may have lost its N, so test the job itself). */
            const bool has_n = std::any_of(J.ref.begin(), J.ref.end(), [](uint8_t x) { return x > 3; })
                               || std::any_of(J.q.begin(), J.q.end(), [](uint8_t x) { return x > 3; });
            if (has_n && kind[t] != RESCUE_PRUNE_FULL) { ok = false; if (mism < 10) fprintf(stderr, "N job not FULL\n"); }
            if (!pp.valid && kind[t] != RESCUE_PRUNE_FULL) {
                ok = false;
                if (mism < 10) fprintf(stderr, "refused scoring but kind=%d\n", kind[t]);
            }
            if (!ok) { mism++; cls_mm[J.cls]++; }
            /* The scalar model is kswv's 8-bit rowMax (scalar_score2's zeroed rising rows); kswv's
             * 16-bit kernels keep a plain row max, so their score2 can differ from it (both widths
             * as in v0.13.0). RB_16BIT checks pruned against full-window kswv16 only. */
            if (!g_16 && scalar_stride > 0 && pp.default_scoring() && (n % scalar_stride) == 0) {
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
    printf("eq: jobs=%ld pass=%ld full=%ld b1=%ld b2=%ld banded_parents=%ld pass1_banded=%ld pass1_guard=%ld "
           "MISMATCHES=%ld scalar_vs_kswv_mismatch=%ld (te2 differs, unconsumed: %ld)\n",
           n, npass, nfull, nb1, nb2, nband, p1_band, p1_guard, mism, scal_mm, te2_diff);
    printf("  scoring: batches=%ld valid=%ld K3=%ld K4=%ld K5=%ld K6=%ld K7=%ld K8=%ld mates_cut_to_8bit=%ld meth "
           "genomic=%ld neutral=%ld collapsed=%ld\n", nbatch, nvalid, nk[3], nk[4], nk[5], nk[6], nk[7], nk[8], ntrunc,
           nmeth[1], nmeth[2], nmeth[3]);
    for (int c = 0; c < 128; c++)
        if (cls_n[c]) printf("  class %3d: jobs=%ld mismatches=%ld\n", c, cls_n[c], cls_mm[c]);
    if (p1_guard) fprintf(stderr, "FAIL: %ld pass-1 guard fallbacks (banded max != S)\n", p1_guard);
    return mism || scal_mm || p1_guard ? 1 : 0;
}

/* ---------------------------------------------------------------------------------------- */
/* timing: single thread, per batch of 2048 dumped jobs, on the SAME parents:                 */
/*   hull  = filter + kswv phase 0 on every non-B1 job (hull or full window) [BAND off]        */
/*   band  = filter + plan + kswv phase 0 on the non-banded + run_pass0 on the banded         */
/* ---------------------------------------------------------------------------------------- */
static void run_time(std::vector<Job> &jobs, int minsc, int max_hits, int reps)
{
    using clk = std::chrono::steady_clock;
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
                if (!k || maxr > kr || maxq > kq) { kr = std::max(kr, maxr); kq = std::max(kq, maxq); k = make_k(kr, kq, rescue_prune_params::defaults(minsc)); }
                Batch P;
                P.ref.reserve(3000000); P.qer.reserve(400000);
                double tf = 0, tp = 0;
                for (int t = 0; t < m; t++) {
                    const Job &J = jobs[base + t];
                    const int len1 = (int)J.ref.size(), len2 = (int)J.q.size();
                    int hb, he;
                    auto a = clk::now();
                    const int kd = rescue_prune_window(J.ref.data(), len1, J.q.data(), len2, minsc, max_hits, &hb, &he);
                    auto b = clk::now();
                    tf += std::chrono::duration<double>(b - a).count();
                    if (kd == RESCUE_PRUNE_B1) continue;
                    if (kd == RESCUE_PRUNE_B2 && mode == 1) {
                        a = clk::now();
                        rescue_band_batch().plan(rescue_prune_last_view(), rescue_prune_params::defaults(minsc),
                                                 len1, len2, hb, he);
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
                for (int g = nbd; g < pn; g += 16) {   // kswv padded cells: 16 x maxLen1 x maxQuanta
                    int ml = 0, mq = 0;
                    for (int l = g; l < std::min(pn, g + 16); l++) { ml = std::max(ml, (int)P.sp[l].len1); mq = std::max(mq, ((int)P.sp[l].len2 + 15) / 16 * 16); }
                    kcells[mode] += 16.0 * ml * mq;
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
    const int minsc = rescue_env_int("RB_MINSC", MINSC_DEFAULT);
    g_16 = rescue_env_int("RB_16BIT", 0) != 0;
    g_kmax = std::min(rescue_env_int("RB_KMAX", 5), rescue_prune_scratch_kmax);
    if (const char *m = getenv("RB_METH")) {
        g_meth = !strcmp(m, "genomic") ? METH_GENOMIC : !strcmp(m, "neutral") ? METH_NEUTRAL
               : !strcmp(m, "collapsed") ? METH_COLLAPSED : !strcmp(m, "random") ? METH_RANDOM : -1;
        if (g_meth < 0) { fprintf(stderr, "RB_METH: want genomic, neutral, collapsed or random\n"); return 2; }
    }
    if (const char *m = getenv("RB_METH_REL")) {
        if (strcmp(m, "conv") && strcmp(m, "expand")) { fprintf(stderr, "RB_METH_REL: want conv or expand\n"); return 2; }
        g_meth_rel = !strcmp(m, "expand");
    }
    if (const char *sc = getenv("RB_SCORING")) {
        if (!strcmp(sc, "random")) g_random_scoring = true;
        else {
            Scoring &f = g_fixed_scoring;
            if (sscanf(sc, "%d,%d,%d,%d,%d,%d", &f.a, &f.b, &f.o_del, &f.e_del, &f.o_ins, &f.e_ins) != 6) {
                fprintf(stderr, "RB_SCORING: want random or a,b,o_del,e_del,o_ins,e_ins\n");
                return 2;
            }
        }
    }
    /* The production default (rescue_prune_max_hits in bwamem_pair.cpp). */
    const int max_hits = rescue_env_int("BWA3_RESCUE_PRUNE_MAX_HITS",
                                        rescue_band_enabled() && minsc == 19 ? 1000 : 400);
    std::vector<Job> jobs;
    if (!strcmp(argv[1], "eq")) {
        const int ngen = atoi(argv[2]);
        rng.seed(strtoull(argv[3], nullptr, 10));
        for (int t = 0; t < ngen; t++) jobs.push_back(gen(t % 10));
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
