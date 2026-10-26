// test/unit/test_xa_lazy_slots.cpp — the paired path builds only the XA slots
// its emitter reads (the pair's primary and the emitted ALT record,
// mem_pe_alt_slot). Every slot the emitter reads must equal the full build;
// every other slot is NULL and freeing it is safe. The emitter's read rule is
// restated here rather than taken from mem_pe_alt_slot.

#include "doctest/doctest.h"

#include "bwa.h"
#include "bwamem.h"
#include "meth_orig_ref_fixture.h"

#include <cstdint>
#include <cstdlib>
#include <string>
#include <vector>

using meth_test::flip;
using meth_test::nt4;
using meth_test::random_bases;

TEST_CASE("mem_gen_alt: the paired path builds exactly the XA slots its emitter reads" *
          doctest::test_suite("unit/xa_regen")) {
    const int kRef = 5000, seg_len = 120;
    std::string ref = random_bases(kRef, 99);
    const std::string seg = ref.substr(100, seg_len);
    // copies of seg with 1..3 substitutions, one per hit
    for (int k = 0; k < 12; ++k) {
        std::string c = seg;
        for (int j = 0; j <= k % 3; ++j) c[10 + 37 * j] = flip(c[10 + 37 * j]);
        ref.replace(400 + 350 * k, seg_len, c);
    }
    meth_test::OrigRefFixture fx(ref);
    mem_opt_t *opt = mem_opt_init();
    std::vector<char> seq;
    for (char b : seg) seq.push_back((char)nt4(b));
    // regions 0..2 are primaries, 3 is the first non-primary (n_pri), an ALT
    // primary, and every primary owns XA hits
    const int n_pri = 3;
    struct AltVariant { const char *name; int score, secondary, is_alt; bool emitted; };
    const AltVariant variants[] = {{"emitted ALT", 110, -1, 1, true},
                                   {"ALT at -T", opt->T, -1, 1, true},
                                   {"ALT one below -T", opt->T - 1, -1, 1, false},
                                   {"ALT below -T", 20, -1, 1, false},
                                   {"secondary ALT", 110, 0, 1, false},
                                   {"not ALT", 110, -1, 0, false}};
    for (const AltVariant &v : variants) {
        for (int z = 0; z < n_pri; ++z) {
            CAPTURE(v.name);
            CAPTURE(z);
            std::vector<mem_alnreg_t> regs;
            auto region = [&](int copy, int score, int secondary, int secondary_all, int is_alt) {
                mem_alnreg_t r{};
                r.rb = copy < 0 ? 100 : 400 + 350 * copy;
                r.re = r.rb + seg_len;
                r.qb = 0;
                r.qe = seg_len;
                r.rid = 0;
                r.score = r.truesc = score;
                r.secondary = secondary;
                r.secondary_all = secondary_all;
                r.is_alt = is_alt;
                r.w = opt->w;
                r.seedcov = r.seedlen0 = seg_len;
                r.meth_hypothesis = -1;
                regs.push_back(r);
            };
            region(-1, 120, -1, -1, 0);
            region(0, 115, -1, -1, 0);
            region(1, 110, -1, -1, 0);
            region(2, v.score, v.secondary, -1, v.is_alt);  // n_pri
            for (int k = 3; k < 12; ++k) region(k, 105 - (k % 3), 0, k % 4, 0);  // hits of primaries 0..3
            mem_alnreg_v a{regs.size(), regs.size(), regs.data()};
            int *hn_full = nullptr, *hn_lazy = nullptr;
            char **full = mem_gen_alt(opt, fx.orig_bns, fx.orig_pac, &a, seg_len, seq.data(), &hn_full, NULL);
            char **lazy = mem_gen_alt(opt, fx.orig_bns, fx.orig_pac, &a, seg_len, seq.data(), &hn_lazy, NULL, z,
                                      mem_pe_alt_slot(opt, &a, n_pri));
            REQUIRE(full != nullptr);
            REQUIRE(lazy != nullptr);
            // the emitter's read rule, restated: XA[z] always; XA[n_pri] when that region
            // is emitted as an ALT record (score >= -T, not secondary, is_alt)
            const bool alt_read = regs[n_pri].score >= opt->T && regs[n_pri].secondary < 0 && regs[n_pri].is_alt;
            CHECK(alt_read == v.emitted);
            int n_skipped = 0;
            for (int k = 0; k < (int)regs.size(); ++k) {
                CAPTURE(k);
                CHECK(hn_lazy[k] == hn_full[k]);
                const bool read = k == z || (k == n_pri && alt_read);
                if (read) {
                    REQUIRE(full[k] != nullptr);  // every read slot has hits in this setup
                    REQUIRE(lazy[k] != nullptr);
                    CHECK(std::string(lazy[k]) == full[k]);
                } else {
                    CHECK(lazy[k] == nullptr);
                    n_skipped += full[k] != nullptr;
                }
            }
            CHECK(n_skipped >= 2);  // the filter skipped real work
            for (int k = 0; k < (int)regs.size(); ++k) free(lazy[k]), free(full[k]);  // NULL slots free safely
            free(lazy);
            free(full);
            free(hn_full);
            free(hn_lazy);
        }
    }
    free(opt);
}

TEST_CASE("mem_pe_alt_slot: a read with only primaries has no ALT record" * doctest::test_suite("unit/xa_regen")) {
    mem_opt_t *opt = mem_opt_init();
    std::vector<mem_alnreg_t> regs(3);
    for (int k = 0; k < 3; ++k) {
        mem_alnreg_t &r = regs[k];
        r = mem_alnreg_t{};
        r.rb = 100 + 200 * k;
        r.re = r.rb + 100;
        r.qe = 100;
        r.score = r.truesc = 100;
        r.secondary = r.secondary_all = -1;
        r.is_alt = 1;  // would be an emitted ALT record were it past n_pri
        r.meth_hypothesis = -1;
    }
    mem_alnreg_v a{regs.size(), regs.size(), regs.data()};
    CHECK(mem_pe_alt_slot(opt, &a, (int)a.n) == -1);      // n_pri == a.n
    CHECK(mem_pe_alt_slot(opt, &a, (int)a.n + 1) == -1);
    CHECK(mem_pe_alt_slot(opt, &a, (int)a.n - 1) == (int)a.n - 1);  // control: the same region past n_pri is read
    free(opt);
}
