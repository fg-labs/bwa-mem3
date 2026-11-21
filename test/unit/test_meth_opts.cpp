// Unit tests for the shared --meth option family (meth_opts_parse /
// meth_opts_check / meth_opts_finalize): the --meth-context value and its
// interaction with the chemistry defaults. The rejections are tested through
// meth_opts_check, which returns them; that finalize prints check's message and
// exits 1 is tested here in a forked child, and through the CLI by the
// meth_context_cli.sh regression.

#include "doctest/doctest.h"
#include "run_in_child.h"
#include "temp_dir.h"
#include "bwamem.h"
#include "meth_opts.h"
#include "meth_xm.h"

#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <initializer_list>
#include <sstream>
#include <string>
#include <sys/wait.h>
#include <unistd.h>
#include <utility>

TEST_SUITE("unit/meth_opts") {

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

TEST_CASE("meth_opts: bare --meth + context all stays emseq/genomic") {
    Fixture f;
    f.parse(OPT_METH, nullptr);
    f.parse(OPT_METH_CONTEXT, "all");
    meth_opts_finalize(f.opt, &f.cli);
    CHECK(f.opt->meth_context == METH_CTX_ALL);
    CHECK(f.opt->meth_scoring == MEM_METH_SCORING_GENOMIC);
}

TEST_CASE("meth_opts: em-seq with neutral scoring is flagged, every other pairing is not") {
    struct Case { const char *chem; const char *scoring; int want; };
    const Case cases[] = {
        {nullptr, "neutral", 1},   // bare --meth is em-seq
        {"emseq", "neutral", 1},
        {"bisulfite", "neutral", 1},
        {"emseq", nullptr, 0},     // em-seq default: genomic
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
        CHECK(meth_opts_emseq_neutral(f.opt) == c.want);
    }
    Fixture off;  // neutral scoring without --meth: the meth_mode term decides
    off.parse(OPT_METH_SCORING, "neutral");
    meth_opts_finalize(off.opt, &off.cli);
    REQUIRE(off.opt->meth_scoring == MEM_METH_SCORING_NEUTRAL);
    CHECK(meth_opts_emseq_neutral(off.opt) == 0);
}

TEST_CASE("meth_opts: check names each cross-option error and accepts what finalize accepts") {
    static const char *const kNoMeth = "--meth-context requires --meth";
    static const char *const kEmseq =
        "--meth-context cpg requires --meth=taps or --meth=5base: under EM-seq/bisulfite every "
        "unmethylated C converts, so conversions are not confined to CpG";
    static const char *const kCollapsed =
        "--meth-context cpg is incompatible with --meth-scoring collapsed: collapsed also frees "
        "reference T against read C, which has no CpG meaning";
    // meth: whether --meth is given; chem: its value (nullptr = bare --meth).
    struct Case { bool meth; const char *chem; const char *context; const char *scoring; const char *want; };
    const Case cases[] = {
        // --meth-context without --meth, whatever the value
        {false, nullptr, "cpg", nullptr, kNoMeth},
        {false, nullptr, "all", nullptr, kNoMeth},
        // cpg under em-seq/bisulfite, explicit scoring or not
        {true, nullptr, "cpg", nullptr, kEmseq},
        {true, "bisulfite", "cpg", "neutral", kEmseq},
        // cpg with collapsed scoring, given explicitly
        {true, "taps", "cpg", "collapsed", kCollapsed},
        {true, "5base", "cpg", "collapsed", kCollapsed},
        // both cpg rejections apply: the em-seq one is reported first
        {true, nullptr, "cpg", "collapsed", kEmseq},
        {true, "bisulfite", "cpg", "collapsed", kEmseq},
        // accepted
        {false, nullptr, nullptr, nullptr, nullptr},
        {true, nullptr, "all", nullptr, nullptr},
        {true, "taps", "cpg", nullptr, nullptr},
        {true, "5base", "cpg", "genomic", nullptr},
        {true, "taps", "all", "collapsed", nullptr},   // collapsed is refused only under cpg
        {true, "5base", nullptr, "collapsed", nullptr},
        {true, "emseq", nullptr, "collapsed", nullptr},
        {true, "emseq", nullptr, "neutral", nullptr},  // warns in finalize, not an error
    };
    for (const Case &c : cases) {
        CAPTURE(c.meth);
        CAPTURE(c.chem ? c.chem : "(bare)");
        CAPTURE(c.context ? c.context : "(default)");
        CAPTURE(c.scoring ? c.scoring : "(default)");
        Fixture f;
        if (c.meth) f.parse(OPT_METH, c.chem);
        if (c.context) f.parse(OPT_METH_CONTEXT, c.context);
        if (c.scoring) f.parse(OPT_METH_SCORING, c.scoring);
        const char *err = meth_opts_check(f.opt, &f.cli);
        if (c.want == nullptr) {
            CHECK(err == nullptr);
        } else {
            CHECK(err != nullptr);
            if (err != nullptr) CHECK(std::string(err) == c.want);  // finalize adds "ERROR: "
        }
    }
}

TEST_CASE("meth_opts: check judges the scoring finalize resolves, not the scoring as given") {
    // TAPS without --meth-scoring resolves to NEUTRAL, so a COLLAPSED value
    // already in the options (here set directly, as a library caller might)
    // is not what finalize would run with: cpg is accepted.
    Fixture f;
    f.parse(OPT_METH, "taps");
    f.parse(OPT_METH_CONTEXT, "cpg");
    f.opt->meth_scoring = MEM_METH_SCORING_COLLAPSED;
    CHECK(meth_opts_check(f.opt, &f.cli) == nullptr);
    meth_opts_finalize(f.opt, &f.cli);
    CHECK(f.opt->meth_scoring == MEM_METH_SCORING_NEUTRAL);
}

TEST_CASE("meth_opts: check judges the options alone when the parser's flags are unset") {
    // An embedder that writes *opt directly leaves meth_cli_t zeroed.
    const meth_cli_t none{};

    Fixture cpg_without_meth;
    cpg_without_meth.opt->meth_context = METH_CTX_CPG;
    const char *err = meth_opts_check(cpg_without_meth.opt, &none);
    CHECK(err != nullptr);
    if (err != nullptr) CHECK(std::string(err) == "--meth-context requires --meth");

    struct Bad { int context, chem, scoring; const char *want; };
    const Bad bad[] = {
        {7, METH_CHEM_EMSEQ, MEM_METH_SCORING_GENOMIC, "--meth-context accepts 'all' (default) or 'cpg'"},
        {METH_CTX_ALL, 9, MEM_METH_SCORING_GENOMIC,
         "--meth accepts 'emseq' (default; aliases em-seq, bisulfite), 'taps', or '5base' (alias 5-base)"},
        {METH_CTX_ALL, METH_CHEM_EMSEQ, -1, "--meth-scoring requires 'collapsed', 'genomic', or 'neutral'"},
    };
    for (const Bad &b : bad) {
        CAPTURE(b.want);
        Fixture f;
        f.opt->meth_mode = 1;
        f.opt->meth_context = b.context;
        f.opt->meth_chem = b.chem;
        f.opt->meth_scoring = b.scoring;
        err = meth_opts_check(f.opt, &none);
        CHECK(err != nullptr);
        if (err != nullptr) CHECK(std::string(err) == b.want);
    }
}

TEST_CASE("meth_opts: check range-checks the stored values, and only under --meth") {
    // The range check reads the stored scoring, before it is resolved: TAPS
    // without --meth-scoring would resolve to NEUTRAL, yet a bad value is
    // still rejected.
    Fixture taps;
    taps.parse(OPT_METH, "taps");
    taps.opt->meth_scoring = -1;
    const char *err = meth_opts_check(taps.opt, &taps.cli);
    CHECK(err != nullptr);
    if (err != nullptr)
        CHECK(std::string(err) == "--meth-scoring requires 'collapsed', 'genomic', or 'neutral'");

    // Without --meth only the --meth-context rule applies: out-of-range
    // chemistry and scoring values are not judged.
    const meth_cli_t none{};
    Fixture off;
    off.opt->meth_chem = 9;
    off.opt->meth_scoring = -1;
    CHECK(meth_opts_check(off.opt, &none) == nullptr);
}

namespace {

std::string slurp(const std::string &path) {
    std::ifstream in(path);
    std::ostringstream ss;
    ss << in.rdbuf();
    return ss.str();
}

}  // namespace

TEST_CASE("meth_opts: check never prints or exits; finalize exits 1 with check's message") {
    bwa_tests::TempDir dir{"meth_opts_check"};
    const std::string out = dir.path() + "/stderr.txt";
    const std::string stdout_path = dir.path() + "/stdout.txt";
    const std::string marker = dir.path() + "/returned.txt";

    // A rejected and an accepted set of options: check returns either way, silently.
    for (const char *context : {"cpg", "all"}) {
        CAPTURE(context);
        std::remove(marker.c_str());
        const int status = bwa_tests::run_in_child(
            [context, &stdout_path, &marker] {
                /* run_in_child captures only stderr: capture stdout too,
                 * unbuffered so the child's _exit(0) cannot drop it. */
                if (freopen(stdout_path.c_str(), "w", stdout) == nullptr) _exit(2);
                setvbuf(stdout, nullptr, _IONBF, 0);
                Fixture f;
                f.parse(OPT_METH, nullptr);  // em-seq: cpg is rejected, all accepted
                f.parse(OPT_METH_CONTEXT, context);
                (void)meth_opts_check(f.opt, &f.cli);
                /* run_in_child's own _exit(0) and an exit(0) inside check look
                 * alike in the wait status; only a return reaches this marker. */
                std::ofstream(marker) << "returned";
            },
            out.c_str());
        CHECK(WIFEXITED(status));
        CHECK(WEXITSTATUS(status) == 0);
        CHECK(slurp(marker) == "returned");
        CHECK(slurp(out).empty());
        CHECK(slurp(stdout_path).empty());
    }

    // finalize on the rejected set: exit 1, and stderr is exactly check's message.
    std::string want;
    {
        Fixture f;
        f.parse(OPT_METH, nullptr);
        f.parse(OPT_METH_CONTEXT, "cpg");
        const char *err = meth_opts_check(f.opt, &f.cli);
        REQUIRE(err != nullptr);
        want = std::string("ERROR: ") + err + "\n";
    }
    const int status = bwa_tests::run_in_child(
        [] {
            Fixture f;
            f.parse(OPT_METH, nullptr);
            f.parse(OPT_METH_CONTEXT, "cpg");
            meth_opts_finalize(f.opt, &f.cli);
        },
        out.c_str());
    CHECK(WIFEXITED(status));
    CHECK(WEXITSTATUS(status) == 1);
    CHECK(slurp(out) == want);
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

namespace {
// The --meth-mapq the options resolve to: argv as (option id, value) pairs, parsed and finalized as
// a subcommand that does (mem) or does not implement it; NULL when meth_opts_check rejects them.
const char *resolved_mapq(std::initializer_list<std::pair<int, const char *>> args, bool supported, int *out) {
    Fixture f;
    f.cli.mapq_supported = supported;
    for (const auto &a : args) f.parse(a.first, a.second);
    const char *err = meth_opts_check(f.opt, &f.cli);
    if (err != NULL) return err;
    meth_opts_finalize(f.opt, &f.cli);
    *out = meth_mapq_conversion_blind(f.opt) ? MEM_METH_MAPQ_CONVERSION_BLIND : f.opt->meth_mapq;
    return NULL;
}
int mapq_of(std::initializer_list<std::pair<int, const char *>> args, bool supported = true) {
    int out = -1;
    const char *err = resolved_mapq(args, supported, &out);
    CHECK(err == nullptr);
    return out;
}
}  // namespace

TEST_CASE("meth_opts: --meth-mapq resolves to conversion-blind only for em-seq genomic where implemented") {
    const int blind = MEM_METH_MAPQ_CONVERSION_BLIND, scoring = MEM_METH_MAPQ_SCORING;
    CHECK(mapq_of({{OPT_METH, nullptr}}) == blind);
    CHECK(mapq_of({{OPT_METH, "emseq"}}) == blind);
    CHECK(mapq_of({{OPT_METH, "bisulfite"}}) == blind);
    CHECK(mapq_of({{OPT_METH, nullptr}, {OPT_METH_SCORING, "collapsed"}}) == scoring);
    CHECK(mapq_of({{OPT_METH, nullptr}, {OPT_METH_SCORING, "neutral"}}) == scoring);
    CHECK(mapq_of({{OPT_METH, "taps"}}) == scoring);
    CHECK(mapq_of({{OPT_METH, "5base"}}) == scoring);
    CHECK(mapq_of({{OPT_METH, "5base"}, {OPT_METH_SCORING, "genomic"}}) == scoring);
    CHECK(mapq_of({{OPT_METH, nullptr}}, false) == scoring); // a subcommand that does not implement it
    CHECK(mapq_of({{OPT_METH, nullptr}, {OPT_METH_MAPQ, "scoring"}}) == scoring);
    CHECK(mapq_of({{OPT_METH_MAPQ, "conversion-blind"}, {OPT_METH, nullptr}}) == blind);
    // collapsed accepts conversion-blind, in either order, as a no-op
    CHECK(mapq_of({{OPT_METH, nullptr}, {OPT_METH_MAPQ, "conversion-blind"}, {OPT_METH_SCORING, "collapsed"}}) == scoring);
    CHECK(mapq_of({{OPT_METH_SCORING, "collapsed"}, {OPT_METH_MAPQ, "conversion-blind"}, {OPT_METH, nullptr}}) == scoring);
}

TEST_CASE("meth_opts: check refuses --meth-mapq without --meth and conversion-blind outside em-seq genomic/collapsed") {
    int out;
    CHECK(std::string(resolved_mapq({{OPT_METH_MAPQ, "scoring"}}, true, &out)) == "--meth-mapq requires --meth");
    // each argv inline: a nested std::initializer_list would leave the inner arrays dangling
    auto refused = [](const char *err) {
        REQUIRE(err != nullptr);
        CHECK(std::string(err).rfind("--meth-mapq conversion-blind applies to em-seq", 0) == 0);
    };
    refused(resolved_mapq({{OPT_METH, "taps"}, {OPT_METH_MAPQ, "conversion-blind"}}, true, &out));
    refused(resolved_mapq({{OPT_METH, "5base"}, {OPT_METH_MAPQ, "conversion-blind"}}, true, &out));
    refused(resolved_mapq({{OPT_METH, nullptr}, {OPT_METH_SCORING, "neutral"}, {OPT_METH_MAPQ, "conversion-blind"}},
                          true, &out));
}

TEST_CASE("meth_opts: the conversion-blind cap applies only under --meth") {
    Fixture f;
    f.opt->meth_mapq = MEM_METH_MAPQ_CONVERSION_BLIND;
    CHECK(!meth_mapq_conversion_blind(f.opt)); // meth_mode unset
    f.opt->meth_mode = 1;
    CHECK(meth_mapq_conversion_blind(f.opt));
    f.opt->meth_mapq = MEM_METH_MAPQ_SCORING;
    CHECK(!meth_mapq_conversion_blind(f.opt));
}

TEST_CASE("meth_opts: --meth-mapq-max-len defaults to 0 (every length), parses a length, and gates by length") {
    Fixture f;
    CHECK(f.opt->meth_mapq_max_len == 0);
    CHECK(meth_mapq_blind_len(f.opt, 1000000));
    f.parse(OPT_METH, nullptr);
    f.parse(OPT_METH_MAPQ_MAX_LEN, "100");
    CHECK(f.cli.set_mapq_max_len == 1);
    CHECK(meth_opts_check(f.opt, &f.cli) == nullptr);
    meth_opts_finalize(f.opt, &f.cli);
    CHECK(f.opt->meth_mapq_max_len == 100);
    CHECK(meth_mapq_blind_len(f.opt, 100));
    CHECK(!meth_mapq_blind_len(f.opt, 101));
    f.parse(OPT_METH_MAPQ_MAX_LEN, "0"); // every length
    CHECK(meth_mapq_blind_len(f.opt, 1000000));
}

TEST_CASE("meth_opts: check refuses --meth-mapq-max-len without --meth") {
    Fixture f;
    f.parse(OPT_METH_MAPQ_MAX_LEN, "50");
    const char *err = meth_opts_check(f.opt, &f.cli);
    REQUIRE(err != nullptr);
    CHECK(std::string(err) == "--meth-mapq-max-len requires --meth");
}

TEST_CASE("meth_opts: the stray-value list routes the --meth-mapq values to --meth-mapq") {
    const char *canon = nullptr;
    CHECK(std::string(meth_opts_stray_flag("conversion-blind", &canon)) == "--meth-mapq");
    CHECK(std::string(canon) == "conversion-blind");
    CHECK(std::string(meth_opts_stray_flag("Scoring", &canon)) == "--meth-mapq");
    CHECK(std::string(canon) == "scoring");
}

}  // TEST_SUITE("unit/meth_opts")
