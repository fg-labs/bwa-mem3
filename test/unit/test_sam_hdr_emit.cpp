// test/unit/test_sam_hdr_emit.cpp — @HD/@SQ precedence in the SAM-text header
// emitter (bwa_print_sam_hdr2, src/bwa.cpp).
//
// This function had NO direct test. test/header_insert_test.cpp covers only how
// -H text is ACCUMULATED (bwa_insert_header/bwa_escape), never how it is
// EMITTED, and the only -H '@HD…' case in the tree was in
// meth_rg_header_test.sh — which passes a LEADING @HD and goes through the
// --meth writer, a different function that was already correct.
//
// That gap hid a spec violation: a user @HD that was not the FIRST record in -H
// suppressed nothing, so the default @HD was emitted as well and the header
// carried two @HD lines (fg-labs/bwa-mem3, fixed alongside this test).
//
// The precedence contract, from lh3/bwa#348 and bwa.c:412-438:
//   @HD : user's (-H) > index sidecar's > the compat target's default,
//         and EXACTLY ONE is emitted, always.
//   @SQ : user's -H block alone if it supplies any; else the sidecar's;
//         else generated from bns.

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>
#include "doctest/doctest.h"
#include "bwa.h"
#include "compat_target.h"

namespace {

// Two contigs, one of them ALT, so the generated @SQ block is non-trivial.
struct TestBns {
    bntann1_t anns[2];
    bntseq_t  bns;
    TestBns() {
        memset(anns, 0, sizeof(anns));
        anns[0].name = const_cast<char *>("chr1");
        anns[0].len = 1000; anns[0].is_alt = 0;
        anns[1].name = const_cast<char *>("chr1_alt");
        anns[1].len = 500;  anns[1].is_alt = 1;
        memset(&bns, 0, sizeof(bns));
        bns.n_seqs = 2;
        bns.anns   = anns;
    }
};

// Render a header to a string. open_memstream gives a real FILE* so the
// function under test is exercised exactly as it is in production.
std::string render(const char *idx_hdr_lines, const char *hdr_line,
                   const compat_target_t *compat) {
    TestBns t;
    char *buf = NULL; size_t len = 0;
    FILE *fp = open_memstream(&buf, &len);
    REQUIRE(fp != NULL);
    // Both inputs are taken as `const char *` and are only ever read from and
    // advanced past internally -- nothing in bwa_print_sam_hdr2 frees or
    // rewrites them -- so the literals can be passed straight through.
    bwa_print_sam_hdr2(&t.bns, idx_hdr_lines, hdr_line, fp, compat);
    fclose(fp);
    std::string out = buf ? std::string(buf, len) : std::string();
    free(buf);
    return out;
}

// Number of lines starting with `prefix`. Pass the record type WITH its
// terminating tab ("@SQ\t", not "@SQ"): the production predicates all match
// "@XX\t", so a bare type would let a malformed "@SQX" line satisfy a count
// that the aligner itself would not have counted.
int count_lines_starting(const std::string &s, const char *prefix) {
    int n = 0;
    size_t pos = 0;
    const size_t plen = strlen(prefix);
    while (pos <= s.size()) {
        const size_t eol = s.find('\n', pos);
        const size_t len = (eol == std::string::npos ? s.size() : eol) - pos;
        if (len >= plen && s.compare(pos, plen, prefix) == 0) ++n;
        if (eol == std::string::npos) break;
        pos = eol + 1;
    }
    return n;
}

// Index of the first line starting with `prefix`, or -1. Same convention as
// count_lines_starting: the prefix carries the record type's tab.
int line_index_of(const std::string &s, const char *prefix) {
    int idx = 0;
    size_t pos = 0;
    const size_t plen = strlen(prefix);
    while (pos <= s.size()) {
        const size_t eol = s.find('\n', pos);
        const size_t len = (eol == std::string::npos ? s.size() : eol) - pos;
        if (len >= plen && s.compare(pos, plen, prefix) == 0) return idx;
        if (eol == std::string::npos) break;
        pos = eol + 1; ++idx;
    }
    return -1;
}

} // namespace

TEST_CASE("@HD: exactly one is emitted, in every -H shape") {
    // THE REGRESSION. Each of these emitted the right count except the third,
    // which emitted two: the default plus the user's.
    struct Case { std::string desc; const char *hdr_line; };
    const std::vector<Case> cases = {
        { "no -H",                        NULL },
        { "-H with a LEADING @HD",        "@HD\tVN:1.6\tSO:coordinate" },
        { "-H with a NON-LEADING @HD",    "@RG\tID:x\tSM:y\n@HD\tVN:1.6\tSO:coordinate" },
        { "-H with @HD last of three",    "@RG\tID:x\n@CO\tnote\n@HD\tVN:1.6" },
        { "-H with no @HD at all",        "@RG\tID:x\tSM:y" },
    };
    for (const Case &c : cases) {
        const std::string out = render(NULL, c.hdr_line, &COMPAT_TARGET_OFF);
        CHECK_MESSAGE(count_lines_starting(out, "@HD\t") == 1,
                      "case: ", c.desc, "\n----\n", out, "----");
    }
}

TEST_CASE("@HD: a leading user @HD wins and is hoisted above @SQ") {
    // Hoisting is what makes the header spec-valid (@HD must come first). It
    // is bwa-mem3-only -- bwa emits -H records after @SQ (bwa.c:438).
    const std::string out = render(NULL, "@HD\tVN:1.6\tSO:coordinate", &COMPAT_TARGET_OFF);
    CHECK(count_lines_starting(out, "@HD\t") == 1);
    CHECK(out.rfind("@HD\tVN:1.6\tSO:coordinate", 0) == 0);   // literally first
    CHECK(line_index_of(out, "@HD\t") < line_index_of(out, "@SQ\t"));
    CHECK(out.find("VN:1.5") == std::string::npos);           // default suppressed
}

TEST_CASE("@HD: a non-leading user @HD suppresses the default but stays inline") {
    // It cannot be hoisted without reordering the user's own records, so it is
    // emitted with the rest of -H after @SQ -- where bwa puts every -H record
    // anyway. What it MUST do is suppress the default; failing to was the bug.
    const std::string out =
        render(NULL, "@RG\tID:x\tSM:y\n@HD\tVN:1.6\tSO:coordinate", &COMPAT_TARGET_OFF);
    CHECK(count_lines_starting(out, "@HD\t") == 1);
    CHECK(out.find("VN:1.5") == std::string::npos);           // default suppressed
    CHECK(out.find("VN:1.6") != std::string::npos);           // the user's survives
    CHECK(line_index_of(out, "@SQ\t") < line_index_of(out, "@HD\t"));   // inline, after @SQ
    CHECK(count_lines_starting(out, "@RG\t") == 1);             // and @RG is not lost
}

TEST_CASE("@HD: the index sidecar supplies one when -H does not") {
    const std::string out =
        render("@HD\tVN:1.4\tSO:queryname\n@CO\tfrom the sidecar", NULL, &COMPAT_TARGET_OFF);
    CHECK(count_lines_starting(out, "@HD\t") == 1);
    CHECK(out.find("VN:1.4") != std::string::npos);
    CHECK(out.find("VN:1.5") == std::string::npos);           // default suppressed
    CHECK(out.find("@CO\tfrom the sidecar") != std::string::npos);
}

TEST_CASE("@HD: the user's -H beats the sidecar's, and only one survives") {
    // Both @HD records lead their own stream here, which is the easy shape: the
    // loser is consumed off the front of the sidecar and simply not printed.
    const std::string out = render("@HD\tVN:1.4\tSO:queryname",
                                   "@HD\tVN:1.6\tSO:coordinate", &COMPAT_TARGET_OFF);
    CHECK(count_lines_starting(out, "@HD\t") == 1);
    CHECK(out.find("VN:1.6") != std::string::npos);           // user's
    CHECK(out.find("VN:1.4") == std::string::npos);           // sidecar's dropped
}

TEST_CASE("@HD: a LOSING sidecar @HD is dropped wherever it sits in the sidecar") {
    // The case above only exercised a LEADING sidecar @HD, which is consumed off
    // the front. A sidecar @HD behind another record cannot be consumed that way,
    // and used to be counted-but-still-printed -- so the stream carried two @HD
    // records whenever -H also supplied one. The BAM writer already filtered
    // these (bam_writer.cpp strips every @HD from the sidecar when -H has one),
    // so this was a SAM-vs-BAM divergence as well as a spec violation.
    struct Case { std::string desc; const char *idx; const char *usr; };
    const std::vector<Case> cases = {
        { "leading user @HD, non-leading sidecar @HD",
          "@RG\tID:z\n@HD\tVN:1.4\tSO:queryname",
          "@HD\tVN:1.6\tSO:coordinate" },
        { "non-leading user @HD, non-leading sidecar @HD",
          "@RG\tID:z\n@HD\tVN:1.4\tSO:queryname",
          "@CO\tnote\n@HD\tVN:1.6\tSO:coordinate" },
        { "non-leading user @HD, leading sidecar @HD",
          "@HD\tVN:1.4\tSO:queryname\n@RG\tID:z",
          "@CO\tnote\n@HD\tVN:1.6\tSO:coordinate" },
    };
    for (const Case &c : cases) {
        const std::string out = render(c.idx, c.usr, &COMPAT_TARGET_OFF);
        CHECK_MESSAGE(count_lines_starting(out, "@HD\t") == 1,
                      "case: ", c.desc, "\n----\n", out, "----");
        CHECK_MESSAGE(out.find("VN:1.6") != std::string::npos,   // the user's wins
                      "case: ", c.desc, "\n----\n", out, "----");
        CHECK_MESSAGE(out.find("VN:1.4") == std::string::npos,   // the sidecar's loses
                      "case: ", c.desc, "\n----\n", out, "----");
        // Dropping the losing @HD must not take the sidecar's other records
        // with it.
        CHECK_MESSAGE(count_lines_starting(out, "@RG\tID:z") == 1,
                      "case: ", c.desc, "\n----\n", out, "----");
    }
}

TEST_CASE("@HD: only the first of several sidecar @HD records survives") {
    // With no -H at all the sidecar supplies the winner, but a sidecar carrying
    // more than one @HD must still yield exactly one record.
    const std::string out =
        render("@HD\tVN:1.4\tSO:queryname\n@RG\tID:z\n@HD\tVN:1.3\tSO:unsorted",
               NULL, &COMPAT_TARGET_OFF);
    CHECK(count_lines_starting(out, "@HD\t") == 1);
    CHECK(out.find("VN:1.4") != std::string::npos);           // first wins
    CHECK(out.find("VN:1.3") == std::string::npos);
    CHECK(count_lines_starting(out, "@RG\tID:z") == 1);
}

TEST_CASE("@HD: a compat target that emits none emits none, but -H still wins") {
    const compat_target_t *mem2 = compat_target_from_name("bwa-mem2");
    REQUIRE(mem2 != nullptr);
    REQUIRE(mem2->emit_hd == 0);

    // No user or sidecar @HD -> nothing at all. bwa-mem2 emits no @HD.
    CHECK(count_lines_starting(render(NULL, NULL, mem2), "@HD\t") == 0);

    // A user @HD is still honored -- the target suppresses only the DEFAULT.
    const std::string out = render(NULL, "@HD\tVN:1.6\tSO:coordinate", mem2);
    CHECK(count_lines_starting(out, "@HD\t") == 1);
    CHECK(out.find("VN:1.6") != std::string::npos);
}

TEST_CASE("@SQ: generated from bns, carrying AH:* on ALT contigs") {
    const std::string out = render(NULL, NULL, &COMPAT_TARGET_OFF);
    CHECK(count_lines_starting(out, "@SQ\t") == 2);
    CHECK(out.find("@SQ\tSN:chr1\tLN:1000\n") != std::string::npos);
    CHECK(out.find("@SQ\tSN:chr1_alt\tLN:500\tAH:*\n") != std::string::npos);
}

TEST_CASE("@SQ: a user -H @SQ block is authoritative and suppresses generation") {
    // Matches both upstreams (bwa.c:429, bwa-mem2 bwa.cpp:534): if -H supplies
    // any @SQ, it is used alone -- no generation, and the sidecar is ignored.
    const std::string out = render("@SQ\tSN:from_sidecar\tLN:7",
                                   "@SQ\tSN:from_user\tLN:9", &COMPAT_TARGET_OFF);
    CHECK(count_lines_starting(out, "@SQ\t") == 1);
    CHECK(out.find("SN:from_user") != std::string::npos);
    CHECK(out.find("SN:from_sidecar") == std::string::npos);
    CHECK(out.find("SN:chr1") == std::string::npos);          // nothing generated
}

TEST_CASE("@SQ: '@SQ' named inside another record does not count as one") {
    // The -H @SQ count decides whether the bns block is generated at all, so a
    // false positive silently drops every contig from the header. Only records
    // that BEGIN with "@SQ\t" count -- a @CO that merely mentions the tag, and a
    // line whose text contains it after the first field, must not.
    const std::string co = render(NULL, "@CO\tsee @SQ\tfor details",
                                  &COMPAT_TARGET_OFF);
    CHECK(count_lines_starting(co, "@SQ\t") == 2);      // still generated from bns
    CHECK(co.find("SN:chr1_alt") != std::string::npos);

    // ... and a real @SQ record that also mentions the tag later on its own
    // line counts exactly once, so generation is still suppressed.
    const std::string one = render(NULL, "@SQ\tSN:from_user\tLN:9\tXX:@SQ\t",
                                   &COMPAT_TARGET_OFF);
    CHECK(count_lines_starting(one, "@SQ\t") == 1);
    CHECK(one.find("SN:chr1") == std::string::npos);
}

TEST_CASE("@SQ: a record type that merely EXTENDS @SQ is not an @SQ") {
    // The record type ends at the tab, so "@SQX" is its own (unknown) type, not
    // an @SQ with a typo. count_SQ requires "@SQ\t", so generation still fires
    // and the malformed record is passed through as an ordinary -H line.
    //
    // This is also what makes the tab mandatory in the assertions above: with a
    // bare "@SQ" prefix, this header counts THREE and the case cannot tell an
    // exact-prefix implementation from a loose one.
    const std::string out = render(NULL, "@SQX\tSN:bogus\tLN:1", &COMPAT_TARGET_OFF);
    CHECK(count_lines_starting(out, "@SQ\t") == 2);       // generated from bns
    CHECK(count_lines_starting(out, "@SQ") == 3);         // ... plus the @SQX line
    CHECK(count_lines_starting(out, "@SQX\t") == 1);      // preserved, not dropped
    CHECK(out.find("SN:chr1_alt") != std::string::npos);
}

TEST_CASE("@SQ: -H text without a trailing newline is still counted and emitted") {
    // -H text arrives unterminated from bwa_insert_header and newline-terminated
    // from a sidecar file; both must yield the same record count. If the final
    // record were dropped, generation would fire and the header would carry the
    // user's @SQ plus all of bns.
    const std::string out = render(NULL, "@SQ\tSN:only\tLN:5", &COMPAT_TARGET_OFF);
    CHECK(count_lines_starting(out, "@SQ\t") == 1);
    CHECK(out.find("SN:only") != std::string::npos);
    CHECK(out.find("SN:chr1") == std::string::npos);
}

TEST_CASE("@HD/@SQ: blank records in -H text are dropped, not emitted") {
    // Header text can carry empty lines (a sidecar with a stray newline). They
    // are not valid records; the emitter drops them rather than writing bare
    // newlines into the header.
    const std::string out = render(NULL, "@RG\tID:a\n\n@RG\tID:b",
                                   &COMPAT_TARGET_OFF);
    CHECK(count_lines_starting(out, "@RG\t") == 2);
    CHECK(out.find("\n\n") == std::string::npos);
}
