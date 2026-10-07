/* Kernel-level differential test for the banded rescue kernels (src/rescue_band.cpp).
 *
 *   rescue_band_kernel_eq <nconfigs> <seed>     exit 1 on any mismatch
 *
 * Runs the two-row kernel rb_dp_wave2 (kernel 2, the default) and the one-row fused kernel
 * rb_dp_core<true> (kernel 1) on the same synthetic 16-lane groups and compares every per-lane
 * output both produce: the best score (gmax), its first row (te) and every row maximum (R, which
 * feeds the suboptimal score). (Kernel 2 also reads the query end directly, kb_chunk / kb_off,
 * where kernel 1 snapshots rows; rescue_band_harness checks that end result against kswv.) The groups vary
 * the band width (1-96), the row count (1-260), each lane's query offset, quantum and length (so
 * lanes pad from different columns and the masked / unmasked split falls anywhere, including at the
 * last steps of a row pair), N bases, and the early-exit targets, in all three scoring forms
 * (RB_SC_DFLT, RB_SC_SYM, RB_SC_GEN with random tables and gap costs). rescue_band_harness checks
 * end results against kswv; this checks the kernels' intermediate per-row output against each
 * other, which the end results do not always expose. */
#include "../src/rescue_band.cpp"

#include <climits>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <random>

#if defined(__aarch64__)
namespace {

struct group_cfg {
    int W, NR;
    int o[16], quanta[16], len2[16];
};

/* Fill w's inputs for one group: a random reference per lane and a query copied from it with
 * mutpct % substitutions, N bases, and the pad / out-of-band sentinel codes the planner writes. */
void build(rb_work &w, const group_cfg &c, std::mt19937 &rng, int mutpct)
{
    const int W = c.W, NR = c.NR, C = W - 1, P = NR + W - 1;
    const int Pp = (P + 15) & ~15, NRp = (NR + 15) & ~15;
    rb_work::fit(w.A, (size_t)Pp * 16);
    rb_work::fit(w.QL, (size_t)Pp * 16);
    rb_work::fit(w.REF, (size_t)NRp * 16);
    rb_work::fit(w.H, (size_t)W * 3 * 16);
    rb_work::fit(w.E, (size_t)(W + 2) * 16);
    rb_work::fit(w.R, (size_t)NRp * 16);
    rb_work::fit(w.SNAP, (size_t)W * 16);
    std::vector<uint8_t> ref((size_t)NRp * 16);
    for (auto &x : ref) x = rng() & 3;
    for (int r = 0; r < NRp; r++)
        for (int l = 0; l < 16; l++) w.REF[r * 16 + l] = r < NR ? ref[r * 16 + l] : 0x80;
    for (int p = 0; p < Pp; p++)
        for (int l = 0; l < 16; l++) {
            const int j = p - C + c.o[l];
            uint8_t v = 0xC0;
            if (j >= 0 && j < c.len2[l]) {
                const int r = j + W / 2 - c.o[l];
                v = (r >= 0 && r < NR && (int)(rng() % 100) >= mutpct) ? ref[r * 16 + l] : (rng() & 3);
                if (rng() % 200 == 0) v = 8;   // N
            } else if (j >= c.len2[l] && j < c.quanta[l]) {
                v = 0x40;
            }
            w.A[p * 16 + l] = v;
        }
    rb_build_ql(w.A.data(), w.QL.data(), Pp);
}

void bounds(const group_cfg &c, int &omax, int &ominq, int &omaskq)
{
    omax = INT_MIN; ominq = INT_MAX; omaskq = INT_MIN;
    for (int l = 0; l < 16; l++) {
        omax = std::max(omax, c.o[l]);
        ominq = std::min(ominq, c.o[l] - c.quanta[l] + 1);
        omaskq = std::max(omaskq, c.o[l] - c.quanta[l] + 1);
    }
}

template <int SC>
void run(int kern, rb_work &w, const group_cfg &c, bool early)
{
    int omax, ominq, omaskq;
    bounds(c, omax, ominq, omaskq);
    if (kern == 2) rb_dp_wave2<SC>(w, c.W, c.NR, omax, ominq, omaskq, early);
    else rb_dp_core<true, SC>(w, c.W, c.NR, omax, ominq, omaskq, early);
}

void run_form(int form, int kern, rb_work &w, const group_cfg &c, bool early)
{
    if (form == 0) run<RB_SC_DFLT>(kern, w, c, early);
    else if (form == 1) run<RB_SC_SYM>(kern, w, c, early);
    else run<RB_SC_GEN>(kern, w, c, early);
}

}   // namespace

#endif

int main(int argc, char **argv)
{
#if !defined(__aarch64__)
    (void) argc; (void) argv;
    printf("SKIP: rescue_band_kernel_eq compares the NEON (16-lane) kernels; this build is not aarch64\n");
    return 0;
#else
    if (argc != 3) {
        fprintf(stderr, "usage: %s <nconfigs> <seed>\n", argv[0]);
        return 2;
    }
    const long n = atol(argv[1]);
    std::mt19937 rng((uint32_t)strtoul(argv[2], nullptr, 10));
    long bad = 0, per_form[3] = {0, 0, 0};
    for (long t = 0; t < n; t++) {
        group_cfg c;
        c.W = 1 + rng() % 96;
        c.NR = 1 + rng() % 260;
        const int q0 = 16 * (1 + rng() % 12), o0 = (int)(rng() % (c.W + 8)) - 4;
        for (int l = 0; l < 16; l++) {
            c.quanta[l] = std::max(16, q0 - 16 * (int)(rng() % 2));
            c.len2[l] = c.quanta[l] - (int)(rng() % 16);
            c.o[l] = o0 - (int)(rng() % std::min(4, c.W));
        }
        const bool early = rng() & 1;
        const int form = (int)(t % 3);
        const uint32_t seed = rng();
        const int mutpct = (int)(rng() % 30);
        rb_work w2, w1;
        {
            std::mt19937 r(seed);
            build(w2, c, r, mutpct);
        }
        w1.A = w2.A; w1.QL = w2.QL; w1.REF = w2.REF;
        rb_work::fit(w1.H, w2.H.size()); rb_work::fit(w1.E, w2.E.size());
        rb_work::fit(w1.R, w2.R.size()); rb_work::fit(w1.SNAP, w2.SNAP.size());
        for (int l = 0; l < 16; l++) w2.target[l] = w1.target[l] = early ? (uint8_t)(rng() % 120) : 0;
        memset(w2.kb_chunk, 0, 16); memset(w2.kb_off, 0, 16);
        if (form != 0) {   // a random table and gap costs; SYM shares one pair for both gap types
            for (int i = 0; i < 16; i++) w2.tbl[i] = (int8_t)((int)(rng() % 13) - 8);
            w2.tbl[0] = (int8_t)(1 + rng() % 6);
            w2.oe_del = (uint8_t)(2 + rng() % 12); w2.e_del = (uint8_t)(1 + rng() % 4);
            w2.oe_ins = (uint8_t)(2 + rng() % 12); w2.e_ins = (uint8_t)(1 + rng() % 4);
            if (form == 1) { w2.oe_del = w2.oe_ins; w2.e_del = w2.e_ins; }
            memcpy(w1.tbl, w2.tbl, 16);
            w1.oe_del = w2.oe_del; w1.e_del = w2.e_del; w1.oe_ins = w2.oe_ins; w1.e_ins = w2.e_ins;
        }
        run_form(form, 2, w2, c, early);
        run_form(form, 1, w1, c, early);
        const int NRp = (c.NR + 15) & ~15;
        const bool ok = !memcmp(w2.gmax, w1.gmax, 16) && !memcmp(w2.te, w1.te, sizeof w2.te)
                        && !memcmp(w2.R.data(), w1.R.data(), (size_t)NRp * 16);
        if (!ok) {
            if (bad < 5) fprintf(stderr, "MISMATCH config %ld form %d W=%d NR=%d early=%d\n", t, form, c.W, c.NR, early);
            bad++;
            per_form[form]++;
        }
    }
    printf("rescue_band_kernel_eq: %ld configs, MISMATCHES=%ld (DFLT %ld, SYM %ld, GEN %ld)\n", n, bad,
           per_form[0], per_form[1], per_form[2]);
    return bad != 0;
#endif
}
