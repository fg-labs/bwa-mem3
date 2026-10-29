// Unit tests for the shared --meth option family (meth_opts_parse /
// meth_opts_finalize): the --meth-context value and its interaction with the
// chemistry defaults. The exiting error paths are covered by the
// meth_context_cli.sh regression, since they terminate the process.

#include "doctest/doctest.h"
#include "bwamem.h"
#include "meth_opts.h"

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
    CHECK(f.opt->meth_mode == 1);
    CHECK(f.opt->meth_context == METH_CTX_ALL);
    CHECK(f.opt->meth_scoring == MEM_METH_SCORING_NEUTRAL);
    CHECK(f.opt->meth_seed_prune == MEM_METH_PRUNE_SPEC30);
}

TEST_CASE("meth_opts: taps + cpg -> METH_CTX_CPG") {
    Fixture f;
    f.parse(OPT_METH, "taps");
    f.parse(OPT_METH_CONTEXT, "cpg");
    meth_opts_finalize(f.opt, &f.cli);
    CHECK(f.cli.set_context == 1);
    CHECK(f.opt->meth_context == METH_CTX_CPG);
    CHECK(f.opt->meth_scoring == MEM_METH_SCORING_NEUTRAL);
}

TEST_CASE("meth_opts: 5base + cpg + genomic -> METH_CTX_CPG, scoring GENOMIC") {
    Fixture f;
    f.parse(OPT_METH, "5base");
    f.parse(OPT_METH_CONTEXT, "cpg");
    f.parse(OPT_METH_SCORING, "genomic");
    meth_opts_finalize(f.opt, &f.cli);
    CHECK(f.opt->meth_context == METH_CTX_CPG);
    CHECK(f.opt->meth_scoring == MEM_METH_SCORING_GENOMIC);
}

TEST_CASE("meth_opts: bare --meth + context all stays emseq/collapsed") {
    Fixture f;
    f.parse(OPT_METH, nullptr);
    f.parse(OPT_METH_CONTEXT, "all");
    meth_opts_finalize(f.opt, &f.cli);
    CHECK(f.opt->meth_context == METH_CTX_ALL);
    CHECK(f.opt->meth_scoring == MEM_METH_SCORING_COLLAPSED);
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
