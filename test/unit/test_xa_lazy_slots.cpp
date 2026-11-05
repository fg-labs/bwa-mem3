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
#include <cstring>
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


TEST_CASE("reg2sam region selection preserves threshold, ALT and secondary boundaries" *
          doctest::test_suite("unit/xa_regen")) {
    mem_opt_t *opt = mem_opt_init();
    mem_alnreg_t regs[2]{};
    regs[0].score = 100;
    mem_alnreg_v a{2, 2, regs};
    struct Case { int score, secondary, alt; bool all, emitted; };
    const Case cases[] = {
        {opt->T - 1, -1, 0, false, false}, {opt->T, -1, 0, false, true},
        {opt->T, -1, 1, false, true}, {100, 0, 0, false, false},
        {100, 0, 1, true, false}, {50, 0, 0, true, true},
        {49, 0, 0, true, false}, {100, INT_MAX, 0, true, true},
        {opt->T - 1, INT_MAX, 1, true, false},
        {opt->T, INT_MAX, 0, false, false}, {100, INT_MAX, 1, true, false},
    };
    for (const Case &c : cases) {
        CAPTURE(c.score); CAPTURE(c.secondary); CAPTURE(c.alt); CAPTURE(c.all);
        regs[1].score = c.score; regs[1].secondary = c.secondary; regs[1].is_alt = c.alt;
        opt->flag = c.all ? MEM_F_ALL : 0;
        CHECK(!!mem_reg2sam_emits_region(opt, &a, 1) == c.emitted);
    }
    free(opt);
}

TEST_CASE("emitted XA owners preserve asymmetric ALT caps and ratio boundaries" *
          doctest::test_suite("unit/xa_regen")) {
    const int len = 100;
    const std::string ref = random_bases(2000, 812);
    meth_test::OrigRefFixture fx(ref);
    std::vector<char> seq;
    for (char b : ref.substr(100, len)) seq.push_back((char)nt4(b));
    mem_opt_t *opt = mem_opt_init();
    opt->XA_drop_ratio = 0.75;  // exactly representable, including the threshold
    std::vector<mem_alnreg_t> regs;
    auto add = [&](int at, int score, int secondary, int owner, bool alt) {
        mem_alnreg_t r{};
        r.rb = at; r.re = at + len; r.qe = len; r.rid = 0;
        r.score = r.truesc = score; r.secondary = secondary;
        r.secondary_all = owner; r.is_alt = alt; r.w = opt->w;
        r.seedcov = r.seedlen0 = len; r.meth_hypothesis = -1;
        regs.push_back(r);
    };
    add(100, 100, -1, -1, false);
    add(300, 100, -1, -1, true);
    add(500, opt->T - 1, -1, -1, false);
    // secondary controls emission, while secondary_all controls XA ownership.
    add(700, 75, 1, 0, false);
    add(900, 76, 1, 0, false);
    add(1100, 75, 0, 1, true);
    add(1300, 76, 0, 1, false);
    add(1500, 74, 0, 1, true);  // below the inclusive ratio boundary
    add(1700, 29, 2, 2, false);
    add(1800, opt->T, -1, -1, false);  // emitted owner at the score boundary
    add(1900, opt->T - 1, 9, 9, false);  // below-T alternate still meets XA_drop_ratio
    const auto before = regs;
    mem_alnreg_v a{regs.size(), regs.size(), regs.data()};
    for (bool xb : {false, true}) for (bool pe : {false, true}) {
        opt->flag = (xb ? MEM_F_XB : 0) | (pe ? MEM_F_PE : 0);
        for (int cap_case = 0; cap_case < 3; ++cap_case) {
            CAPTURE(xb); CAPTURE(pe); CAPTURE(cap_case);
            opt->max_XA_hits = cap_case == 0 ? 1 : 2;
            opt->max_XA_hits_alt = cap_case == 1 ? 1 : 3;
            int *full_hn = nullptr, *selected_hn = nullptr;
            // The original unrestricted generator is the output oracle.
            char **full = mem_gen_alt(opt, fx.orig_bns, fx.orig_pac, &a, len,
                                      seq.data(), &full_hn);
            char **selected = mem_gen_alt(opt, fx.orig_bns, fx.orig_pac, &a, len,
                                          seq.data(), &selected_hn, nullptr, -1, -1, true);
            REQUIRE(full != nullptr); REQUIRE(selected != nullptr);
            const bool expected_slot[] = {cap_case == 2, cap_case != 1};
            for (size_t k = 0; k < regs.size(); ++k) {
                CAPTURE(k);
                const int expected_hn = k < 2 ? 2 : (k == 2 || k == 9) ? 1 : 0;
                CHECK(full_hn[k] == expected_hn);
                CHECK(selected_hn[k] == expected_hn);
                if ((k < 2 && expected_slot[k]) || k == 9) {
                    REQUIRE(full[k] != nullptr); REQUIRE(selected[k] != nullptr);
                    CHECK(std::string(selected[k]) == full[k]);
                } else CHECK(selected[k] == nullptr);
                free(full[k]); free(selected[k]);
            }
            CHECK(std::memcmp(regs.data(), before.data(), regs.size() * sizeof(regs[0])) == 0);
            free(full); free(selected); free(full_hn); free(selected_hn);
        }
    }
    free(opt);
}

// Removing the selective-generation change must fail the reuse-count assertion:
// the unrestricted path regenerates and memoizes hits of a below-T owner.
TEST_CASE("reg2sam: un-emitted owners do not populate the XA memo" *
          doctest::test_suite("unit/xa_regen")) {
    const int len = 120;
    std::string ref = random_bases(2000, 711);
    const std::string segment = ref.substr(100, len);
    for (int at : {400, 700, 1000}) ref.replace(at, len, segment);
    meth_test::OrigRefFixture fx(ref);
    mem_opt_t *opt = mem_opt_init();
    opt->max_XA_hits = opt->max_XA_hits_alt = 10;
    SUBCASE("four-field XA") {}
    SUBCASE("six-field XB") { opt->flag |= MEM_F_XB; }
    std::vector<char> seq;
    for (char b : segment) seq.push_back((char)nt4(b));
    std::vector<mem_alnreg_t> regs;
    auto region = [&](int at, int score, int secondary, int owner, bool alt) {
        mem_alnreg_t r{};
        r.rb = at; r.re = at + len; r.qe = len; r.rid = 0;
        r.score = score; r.truesc = len; r.sub = len;
        r.secondary = secondary; r.secondary_all = owner; r.is_alt = alt;
        r.w = opt->w; r.seedcov = r.seedlen0 = len; r.meth_hypothesis = -1;
        regs.push_back(r);
    };
    region(100, 120, -1, -1, false);  // emitted primary
    region(400, opt->T - 1, -1, -1, false);  // discarded owner
    region(700, opt->T - 1, 1, 1, false);
    region(1000, opt->T - 1, 1, 1, false);
    region(400, 120, -1, -1, true);  // emitted ALT supplementary
    region(700, 110, 4, 4, false);
    region(1000, 110, 4, 4, false);
    mem_alnreg_v a{regs.size(), regs.size(), regs.data()};
    const std::vector<mem_alnreg_t> before = regs;
    char name[] = "xa6";
    bseq1_t read{};
    read.name = name; read.seq = seq.data(); read.l_seq = len;

    // Independent unrestricted-generation oracle: explicitly emit the two
    // constructed owners. Do not call the proposed shared selection helper.
    int *hn = nullptr;
    char **xa = mem_gen_alt(opt, fx.orig_bns, fx.orig_pac, &a, len, seq.data(), &hn);
    REQUIRE(xa != nullptr);
    REQUIRE(xa[1] != nullptr);  // the discarded owner has real generation work
    REQUIRE(xa[4] != nullptr);
    CHECK(hn[0] == 0);
    CHECK(hn[1] == 2);
    CHECK(hn[4] == 2);
    SUBCASE("cap omits strings without losing HN") {
        opt->max_XA_hits = opt->max_XA_hits_alt = 1;
        int *capped_hn = nullptr;
        char **capped = mem_gen_alt(opt, fx.orig_bns, fx.orig_pac, &a, len, seq.data(),
                                    &capped_hn, nullptr, -1, -1, true);
        REQUIRE(capped != nullptr);
        for (size_t k = 0; k < regs.size(); ++k) {
            CHECK(capped_hn[k] == hn[k]);
            CHECK(capped[k] == nullptr);
        }
        free(capped); free(capped_hn);
        opt->max_XA_hits = opt->max_XA_hits_alt = 10;
    }
    int *selected_hn = nullptr;
    char **selected = mem_gen_alt(opt, fx.orig_bns, fx.orig_pac, &a, len, seq.data(),
                                  &selected_hn, nullptr, -1, -1, true);
    REQUIRE(selected != nullptr);
    for (size_t k = 0; k < regs.size(); ++k) {
        CHECK(selected_hn[k] == hn[k]);
        if (k == 4) {
            REQUIRE(selected[k] != nullptr);
            CHECK(std::string(selected[k]) == xa[k]);
        } else CHECK(selected[k] == nullptr);
        free(selected[k]);
    }
    free(selected); free(selected_hn);
    mem_aln_t records[2];
    records[0] = mem_reg2aln(opt, fx.orig_bns, fx.orig_pac, len, seq.data(), &regs[0]);
    records[1] = mem_reg2aln(opt, fx.orig_bns, fx.orig_pac, len, seq.data(), &regs[4]);
    records[1].flag |= 0x800;
    for (int k = 0; k < 2; ++k) {
        const int owner = k == 0 ? 0 : 4;
        records[k].XA = xa[owner]; records[k].HN = hn[owner];
    }
    kstring_t expected{0, 0, nullptr};
    for (int k = 0; k < 2; ++k)
        mem_aln2sam(opt, fx.orig_bns, &expected, &read, 2, records, k, nullptr);
    const long reused_before = mem_gen_alt_n_reused;
    mem_reg2sam(opt, fx.orig_bns, fx.orig_pac, &read, &a, 0, nullptr);
    REQUIRE(read.sam != nullptr);
    CHECK(std::string(read.sam) == expected.s);
    // Only the two retained duplicates should execute: one miss, one reuse.
    CHECK(mem_gen_alt_n_reused - reused_before == 1);
    CHECK(std::memcmp(regs.data(), before.data(), regs.size() * sizeof(regs[0])) == 0);
    for (auto &record : records) free(record.cigar);
    for (size_t k = 0; k < regs.size(); ++k) free(xa[k]);
    free(xa); free(hn); free(expected.s); free(read.sam); free(opt);
}
