// test/unit/test_compat_target.cpp — the `--compat` target table.
//
// The table in src/compat_target.cpp is DATA transcribed from upstream
// sources, and one of its rows (bwa-mem) is deliberately not reachable from
// the CLI. Both facts make it prone to silent rot: a transcription error is
// invisible until someone diffs against the real aligner, and an unreachable
// row is exercised by nothing else. So assert every field of every row here,
// against the same citations the table comments carry.
//
// Verified 2026-07-24 against ~/work/git/bwa @ 91fb301 (0.7.19-r1273) and
// bwa-mem2 v2.2.1. If a row changes, the upstream evidence must change with it.

#include <cstring>
#include <string>
#include "doctest/doctest.h"
#include "compat_target.h"

namespace {
// Look a row up by canonical name, failing the test rather than returning NULL
// so a missing row reports as a named failure instead of a segfault.
const compat_target_t *row(const char *name) {
    const compat_target_t *t = compat_target_from_name(name);
    REQUIRE_MESSAGE(t != nullptr, "no compat target row named ", name);
    return t;
}
} // namespace

TEST_CASE("compat target `off` is bwa-mem3's native output") {
    const compat_target_t *t = row("off");
    CHECK(std::string(t->name) == "off");
    CHECK(t->unavailable_reason == nullptr);
    CHECK(t->alias == nullptr);
    // Every output-shaping field is permissive: `off` must be exactly
    // equivalent to not passing --compat at all.
    CHECK(t->emit_hd      == 1);
    CHECK(t->read_sidecar == 1);
    CHECK(t->emit_mq      == 1);
    CHECK(t->emit_hn      == 1);
    // NULL, not a literal: each output path keeps the default it already
    // emits. The SAM text and BAM paths disagree ("VN:1.5 ... GO:query" vs
    // "VN:1.6 SO:unsorted"); pinning a string here would change default output
    // on one of them, which --compat must never do.
    CHECK(t->hd_line == nullptr);
}

TEST_CASE("compat target `bwa-mem2` matches bwa-mem2 v2.2.1") {
    const compat_target_t *t = row("bwa-mem2");
    CHECK(std::string(t->name) == "bwa-mem2");
    CHECK(t->unavailable_reason == nullptr);
    CHECK(std::string(t->alias) == "mem2");
    // bwa-mem2's bwa_print_sam_hdr (src/bwa.cpp:523) has no @HD logic at all:
    // no n_HD counter, no default emission. bwa gained one only in 0.7.18
    // (6b18630), after bwa-mem2 forked at 0.7.17.
    CHECK(t->emit_hd == 0);
    // The <prefix>.hdr / <baseprefix>.dict sidecar is bwa-mem3-only: a port of
    // lh3/bwa#348, which lh3 closed unmerged. grep '\.hdr\|\.dict' is empty in
    // both upstreams, so there is nothing for a compat target to load.
    CHECK(t->read_sidecar == 0);
    // MQ:i arrived in bwa via lh3/bwa#330 (merged 2022-03-06), after the fork.
    CHECK(t->emit_mq == 0);
    // HN:i exists in neither upstream.
    CHECK(t->emit_hn == 0);
}

TEST_CASE("compat target `bwa-mem` is fully specified but not selectable") {
    const compat_target_t *t = row("bwa-mem");
    CHECK(std::string(t->name) == "bwa-mem");
    // Not selectable: bwa-mem3 and bwa still differ on the ALIGNMENTS, so
    // offering the target would advertise a byte-identity we cannot deliver.
    // The row exists because its output shaping is fully determined by the
    // evidence below; when the alignment work lands, this becomes nullptr.
    REQUIRE(t->unavailable_reason != nullptr);
    CHECK(std::string(t->unavailable_reason).find("byte-identity") != std::string::npos);
    // bwa.c:426 — emitted unconditionally when -H supplies no @HD.
    REQUIRE(t->emit_hd == 1);
    REQUIRE(t->hd_line != nullptr);
    CHECK(std::string(t->hd_line) == "@HD\tVN:1.5\tSO:unsorted\tGO:query");
    // lh3/bwa#348 closed unmerged; bwa has no sidecar.
    CHECK(t->read_sidecar == 0);
    // THE field a single compat boolean could not express: bwa DOES emit MQ:i
    // (bwamem.c:935) while bwa-mem2 does not. This asymmetry is why --compat
    // is a target enum rather than a flag bit.
    CHECK(t->emit_mq == 1);
    CHECK(t->emit_hn == 0);
}

TEST_CASE("bwa and bwa-mem2 rows differ exactly where the upstreams do") {
    const compat_target_t *mem  = row("bwa-mem");
    const compat_target_t *mem2 = row("bwa-mem2");
    // Both drop the bwa-mem3-only sidecar and the bwa-mem3-only HN:i tag...
    CHECK(mem->read_sidecar == mem2->read_sidecar);
    CHECK(mem->emit_hn      == mem2->emit_hn);
    // ...and disagree on precisely the two fields the fork point explains:
    // @HD (bwa 0.7.18, 6b18630) and MQ:i (lh3/bwa#330), both post-0.7.17.
    CHECK(mem->emit_hd != mem2->emit_hd);
    CHECK(mem->emit_mq != mem2->emit_mq);
}

TEST_CASE("compat target lookup: aliases, unknown names, NULL") {
    // `mem2` is the documented short spelling and must resolve to the same row.
    CHECK(compat_target_from_name("mem2") == compat_target_from_name("bwa-mem2"));
    // Unknown names are NULL so the caller can say "unknown target"...
    CHECK(compat_target_from_name("bwa-mem4") == nullptr);
    CHECK(compat_target_from_name("") == nullptr);
    CHECK(compat_target_from_name(nullptr) == nullptr);
    // ...while a RECOGNIZED but unavailable name resolves, so the caller can
    // tell the user the real reason instead of "unknown target".
    const compat_target_t *t = compat_target_from_name("bwa-mem");
    REQUIRE(t != nullptr);
    CHECK(t->unavailable_reason != nullptr);
    // Lookup is exact: no prefix or case folding.
    CHECK(compat_target_from_name("bwa") == nullptr);
    CHECK(compat_target_from_name("BWA-MEM2") == nullptr);
    CHECK(compat_target_from_name("bwa-mem2 ") == nullptr);
}

TEST_CASE("COMPAT_TARGET_OFF is the row `off` resolves to") {
    // mem_opt_init() stores &COMPAT_TARGET_OFF, and main_mem tests
    // `opt->compat != &COMPAT_TARGET_OFF` to decide whether a target is
    // active. If --compat=off resolved to a different row with identical
    // fields, that pointer comparison would treat it as active and wrongly
    // reject --compat=off --fast.
    CHECK(compat_target_from_name("off") == &COMPAT_TARGET_OFF);
}

TEST_CASE("every table row is reachable by name, and rows are self-consistent") {
    int n = 0;
    const compat_target_t *const *rows = compat_targets(&n);
    REQUIRE(rows != nullptr);
    CHECK(n >= 3);
    for (int i = 0; i < n; ++i) {
        const compat_target_t *t = rows[i];
        REQUIRE(t != nullptr);
        REQUIRE(t->name != nullptr);
        CHECK_MESSAGE(compat_target_from_name(t->name) == t,
                      "row ", i, " (", t->name, ") is not reachable by name");
        if (t->alias != nullptr)
            CHECK_MESSAGE(compat_target_from_name(t->alias) == t,
                          "row ", t->name, " is not reachable by its alias");
        // A pinned @HD is meaningless unless the row emits one.
        if (t->hd_line != nullptr)
            CHECK_MESSAGE(t->emit_hd == 1,
                          "row ", t->name, " pins hd_line but does not emit @HD");
    }
}

TEST_CASE("the selectable list is derived from the table") {
    // Pinned exactly. A membership loop would be the obvious test, but
    // `find(t->name)` cannot express it: "bwa-mem" is a prefix of "bwa-mem2",
    // so the unavailable bwa-mem row would look listed. Splitting the string
    // to fix that restates, at length, what this one comparison already says --
    // any table change that would break membership breaks equality first, with
    // a clearer failure message.
    CHECK(std::string(compat_target_selectable_list()) == "bwa-mem2 (alias mem2), off");
}
