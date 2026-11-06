// Unit tests for the shared --meth option family (meth_opts_parse /
// meth_opts_finalize): the --meth-context value and its interaction with the
// chemistry defaults. The exiting error paths are covered by the
// meth_context_cli.sh regression, since they terminate the process.

#include "doctest/doctest.h"
#include "bwamem.h"
#include "meth_opts.h"
#include "aln/aln_meth.h"

#include <cstdlib>
#include <string>

namespace {

struct Fixture {
    mem_opt_t *opt;
    meth_cli_t cli;
    Fixture() : opt(mem_opt_init()), cli() {}
    ~Fixture() { free(opt); }
    void parse(int id, const char *arg) {
        REQUIRE(meth_opts_parse(id, arg, opt, &cli) == 1);
    }
};

}  // namespace

TEST_CASE("meth_opts: 5base + no context -> METH_CTX_ALL, scoring NEUTRAL") {
    Fixture f;
    f.parse(OPT_METH, "5base");
    meth_opts_finalize(f.opt, &f.cli);
    CHECK(f.opt->meth.meth_mode == 1);
    CHECK(f.opt->meth.meth_context == METH_CTX_ALL);
    CHECK(f.opt->meth.meth_scoring == MEM_METH_SCORING_NEUTRAL);
    CHECK(f.opt->meth_seed_prune == MEM_METH_PRUNE_SPEC30);
}

TEST_CASE("meth_opts: taps + cpg -> METH_CTX_CPG") {
    Fixture f;
    f.parse(OPT_METH, "taps");
    f.parse(OPT_METH_CONTEXT, "cpg");
    meth_opts_finalize(f.opt, &f.cli);
    CHECK(f.cli.set_context == 1);
    CHECK(f.opt->meth.meth_context == METH_CTX_CPG);
    CHECK(f.opt->meth.meth_scoring == MEM_METH_SCORING_NEUTRAL);
}

TEST_CASE("meth_opts: 5base + cpg + genomic -> METH_CTX_CPG, scoring GENOMIC") {
    Fixture f;
    f.parse(OPT_METH, "5base");
    f.parse(OPT_METH_CONTEXT, "cpg");
    f.parse(OPT_METH_SCORING, "genomic");
    meth_opts_finalize(f.opt, &f.cli);
    CHECK(f.opt->meth.meth_context == METH_CTX_CPG);
    CHECK(f.opt->meth.meth_scoring == MEM_METH_SCORING_GENOMIC);
}

TEST_CASE("meth_opts: bare --meth + context all stays emseq/collapsed") {
    Fixture f;
    f.parse(OPT_METH, nullptr);
    f.parse(OPT_METH_CONTEXT, "all");
    meth_opts_finalize(f.opt, &f.cli);
    CHECK(f.opt->meth.meth_context == METH_CTX_ALL);
    CHECK(f.opt->meth.meth_scoring == MEM_METH_SCORING_COLLAPSED);
}

TEST_CASE("meth_opts: em-seq with neutral scoring is flagged, every other pairing is not") {
    struct Case { const char *chem; const char *scoring; int want; };
    const Case cases[] = {
        {nullptr, "neutral", 1},   // bare --meth is em-seq
        {"emseq", "neutral", 1},
        {"bisulfite", "neutral", 1},
        {"emseq", nullptr, 0},     // em-seq default: collapsed
        {"emseq", "genomic", 0},
        {"emseq", "collapsed", 0},
        {"taps", nullptr, 0},      // TAPS default: neutral, sparse conversions
        {"taps", "neutral", 0},
        {"5base", "neutral", 0},
    };
    for (const Case &c : cases) {
        CAPTURE(c.chem ? c.chem : "(bare)");
        CAPTURE(c.scoring ? c.scoring : "(default)");
        Fixture f;
        f.parse(OPT_METH, c.chem);
        if (c.scoring) f.parse(OPT_METH_SCORING, c.scoring);
        meth_opts_finalize(f.opt, &f.cli);
        CHECK(meth_opts_emseq_neutral(&f.opt->meth) == c.want);
    }
    Fixture off;  // neutral scoring without --meth: the meth_mode term decides
    off.parse(OPT_METH_SCORING, "neutral");
    meth_opts_finalize(off.opt, &off.cli);
    REQUIRE(off.opt->meth.meth_scoring == MEM_METH_SCORING_NEUTRAL);
    CHECK(meth_opts_emseq_neutral(&off.opt->meth) == 0);
}

TEST_CASE("meth_opts: explicit --meth-seed-prune=off survives finalize") {
    Fixture f;
    f.parse(OPT_METH, nullptr);
    f.parse(OPT_METH_SEED_PRUNE, "off");
    meth_opts_finalize(f.opt, &f.cli);
    CHECK(f.opt->meth_seed_prune == MEM_METH_PRUNE_OFF);
}

TEST_CASE("meth_opts: an option outside the family is not claimed") {
    Fixture f;
    CHECK(meth_opts_parse('k', "19", f.opt, &f.cli) == 0);
    CHECK(f.cli.set_meth == 0);
}

TEST_CASE("meth_opts: stray-value list names cpg and every chemistry") {
    bool cpg = false, taps = false;
    for (const char *const *v = meth_opts_stray_values(); *v != nullptr; ++v) {
        if (std::string(*v) == "cpg") cpg = true;
        if (std::string(*v) == "taps") taps = true;
    }
    CHECK(cpg);
    CHECK(taps);
}

TEST_CASE("meth_opts: aln and mem share named defaults and set-as-failed parsing") {
    Fixture f;
    bwa_aln::aln_meth_opt_t aln;
    CHECK(aln.m.meth_mode == f.opt->meth.meth_mode);
    CHECK(aln.m.meth_scoring == f.opt->meth.meth_scoring);
    CHECK(aln.m.meth_chem == f.opt->meth.meth_chem);
    CHECK(aln.m.meth_tags == f.opt->meth.meth_tags);
    CHECK(aln.m.meth_context == f.opt->meth.meth_context);
    CHECK(aln.m.meth_set_as_failed == f.opt->meth.meth_set_as_failed);
    CHECK(aln.m.meth_set_as_failed == 0);
    for (const char *value : {"f", "r"}) {
        meth_cli_t cli = {};
        REQUIRE(meth_opts_parse_shared(OPT_METH, "taps", &aln.m, &cli) == 1);
        REQUIRE(meth_opts_parse_shared(OPT_METH_SET_AS_FAILED, value, &aln.m, &cli) == 1);
        f.parse(OPT_METH, "taps");
        f.parse(OPT_METH_SET_AS_FAILED, value);
        meth_opts_finalize_shared(&aln.m, &cli);
        meth_opts_finalize(f.opt, &f.cli);
        CHECK(cli.set_failed == 1);
        CHECK(f.cli.set_failed == 1);
        CHECK(aln.m.meth_set_as_failed == value[0]);
        CHECK(aln.m.meth_set_as_failed == f.opt->meth.meth_set_as_failed);
        CHECK(aln.m.meth_scoring == f.opt->meth.meth_scoring);
    }
}
