// test/unit/test_sq_line_format.cpp — the shared generated-@SQ formatter and
// the shared header-text type predicate (src/bwa.cpp).
//
// bwa_format_sq_line is the one definition of "what a generated @SQ record
// contains", shared by the SAM-text, --bam and --meth writers. It exists
// because those three used to build the record independently and drifted:
// AH:* was correct in only one of them, so --bam and --meth silently dropped
// ALT status for every ALT-aware reference until each copy was fixed
// separately (fg-labs/bwa-mem3#281). Collapsing them to one definition, so
// that class of drift cannot recur, is fg-labs/bwa-mem3#289.
//
// Contents must match what bwa (bwa.c:430-433) and bwa-mem2 (bwa.cpp:535-548)
// emit, so these assertions are against the upstream layout, not against
// whatever the current implementation happens to produce.

#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>
#include "doctest/doctest.h"
#include "bwa.h"

namespace {
// A bntann1_t carrying just the fields the formatter reads. `name` is not
// owned by the struct in bwa, so pointing it at a literal is correct here.
bntann1_t ann(const char *name, int64_t len, int is_alt) {
    bntann1_t a{};
    a.name   = const_cast<char *>(name);
    a.len    = (int32_t)len;
    a.is_alt = is_alt;
    return a;
}

std::string format(const bntann1_t &a) {
    kstring_t s = {0, 0, NULL};
    REQUIRE(bwa_format_sq_line(&s, &a) == 0);
    std::string out = s.s ? s.s : "";
    free(s.s);
    return out;
}
} // namespace

TEST_CASE("bwa_format_sq_line: primary contig is SN + LN, no AH") {
    bntann1_t a = ann("chr1", 248956422, 0);
    CHECK(format(a) == "@SQ\tSN:chr1\tLN:248956422");
}

TEST_CASE("bwa_format_sq_line: ALT contig gains AH:* after LN") {
    // Field ORDER matters, not just presence: bwa appends AH immediately after
    // LN, and --meth appends its M5/UR/AS/SP after that. Asserting the exact
    // string pins the layout the sidecar-enrichment path builds on.
    bntann1_t a = ann("chr1_alt", 5000, 1);
    CHECK(format(a) == "@SQ\tSN:chr1_alt\tLN:5000\tAH:*");
}

TEST_CASE("bwa_format_sq_line: no trailing newline") {
    // Callers add their own terminator -- err_fputc('\n') on the SAM path,
    // sam_hdr_add_lines' length argument on the BAM paths -- and --meth
    // appends more tags after this. A newline here would corrupt all three.
    bntann1_t a = ann("chr2", 100, 0);
    const std::string s = format(a);
    REQUIRE(!s.empty());
    CHECK(s.back() != '\n');
    CHECK(s.find('\n') == std::string::npos);
}

TEST_CASE("bwa_format_sq_line: appends, so a buffer can be reused") {
    // Every caller reuses one kstring across contigs by resetting l to 0.
    // Appending (not overwriting) is what makes --meth's tag enrichment work.
    kstring_t s = {0, 0, NULL};
    bntann1_t a = ann("chr1", 10, 0);
    CHECK(bwa_format_sq_line(&s, &a) == 0);
    const size_t after_first = s.l;
    bntann1_t b = ann("chr2", 20, 0);
    CHECK(bwa_format_sq_line(&s, &b) == 0);
    CHECK(s.l > after_first);
    CHECK(std::string(s.s) == "@SQ\tSN:chr1\tLN:10@SQ\tSN:chr2\tLN:20");
    s.l = 0;                       // the reset every caller performs
    CHECK(bwa_format_sq_line(&s, &b) == 0);
    CHECK(std::string(s.s, s.l) == "@SQ\tSN:chr2\tLN:20");
    free(s.s);
}

TEST_CASE("bwa_format_sq_line: reports success so a caller can refuse to emit") {
    // The formatter reserves the record's worst case up front and returns a
    // status, because ksprintf/kputs swallow their own failures: without it, a
    // caller cannot distinguish a formatted record from a failed append and
    // would fputs a NULL or hand a truncated line to sam_hdr_add_lines. Only
    // the success half is assertable without an allocator hook -- pin it, since
    // all three callers now branch on it.
    kstring_t s = {0, 0, NULL};
    bntann1_t a = ann("chr1_alt", 5000, 1);
    CHECK(bwa_format_sq_line(&s, &a) == 0);
    REQUIRE(s.s != NULL);                     // what the callers rely on
    CHECK(s.l == strlen("@SQ\tSN:chr1_alt\tLN:5000\tAH:*"));
    // The reservation must cover the AH:* branch too, so the record is whole
    // and NUL-terminated within the capacity that was reserved.
    CHECK(s.m > s.l);
    CHECK(s.s[s.l] == '\0');
    free(s.s);
}

TEST_CASE("bwa_format_sq_line: a long contig name is not truncated") {
    // The SAM-text path used to snprintf into char[512] and would silently
    // truncate. kstring grows instead.
    const std::string name(600, 'x');
    bntann1_t a = ann(name.c_str(), 1234, 0);
    CHECK(format(a) == "@SQ\tSN:" + name + "\tLN:1234");
}

TEST_CASE("bwa_hdr_text_has_type matches whole records only") {
    const char *hdr =
        "@HD\tVN:1.5\tSO:unsorted\n"
        "@SQ\tSN:chr1\tLN:100\n"
        "@CO\tSQ is mentioned in this comment\n";
    CHECK(bwa_hdr_text_has_type(hdr, "@HD\t") == 1);   // first line
    CHECK(bwa_hdr_text_has_type(hdr, "@SQ\t") == 1);   // after a newline
    CHECK(bwa_hdr_text_has_type(hdr, "@CO\t") == 1);
    CHECK(bwa_hdr_text_has_type(hdr, "@RG\t") == 0);   // absent

    // A record type named only INSIDE another line must not match -- the
    // comment above says "SQ" but carries no @SQ record of its own.
    CHECK(bwa_hdr_text_has_type("@CO\tmentions @SQ\there\n", "@SQ\t") == 0);

    // A leading match must be at position 0, not merely present.
    CHECK(bwa_hdr_text_has_type("x@HD\tVN:1.5\n", "@HD\t") == 0);

    CHECK(bwa_hdr_text_has_type(NULL, "@HD\t") == 0);
    CHECK(bwa_hdr_text_has_type("", "@HD\t") == 0);
    CHECK(bwa_hdr_text_has_type("@HD\tVN:1.5\n", NULL) == 0);

    // An empty type matches nothing. Returning 1 (every record "starts with"
    // the empty string) would make a caller that built its tag from an empty
    // variable see records that are not there.
    CHECK(bwa_hdr_text_has_type(hdr, "") == 0);

    // A record on the LAST line of unterminated text still matches -- header
    // text arrives without a trailing newline from bwa_insert_header.
    CHECK(bwa_hdr_text_has_type("@SQ\tSN:chr1\tLN:1\n@RG\tID:a", "@RG\t") == 1);
}

TEST_CASE("bwa_hdr_text_has_type has no tag-length ceiling") {
    // This is public API, so the type is whatever the caller passes -- not
    // necessarily a 4-char "@XX\t". A length-capped implementation would
    // truncate a longer type and then match on the truncation, reporting a
    // record that is not present.
    const char *hdr = "@HD\tVN:1.5\n@COMMENTARY\tx\n";
    CHECK(bwa_hdr_text_has_type(hdr, "@COMMENTARY\t") == 1);
    // "@COMMENTX\t" shares its first 8 characters with the record above, so a
    // truncating implementation would wrongly report a match.
    CHECK(bwa_hdr_text_has_type(hdr, "@COMMENTX\t") == 0);
    // A type longer than any record cannot match.
    CHECK(bwa_hdr_text_has_type(hdr, "@COMMENTARY\tx\ty\n") == 0);
}

TEST_CASE("bwa_hdr_next_line iterates records, with or without a trailing newline") {
    auto records = [](const char *text) {
        std::vector<std::string> out;
        const char *p = text, *line; size_t len;
        while (bwa_hdr_next_line(&p, &line, &len)) out.emplace_back(line, len);
        return out;
    };

    // Trailing newline must NOT yield a phantom empty final record -- the
    // sidecar loader returns text without one, htslib emits it with one, and
    // the six call sites this replaced disagreed on which they handled.
    const std::vector<std::string> want{"@HD\tVN:1.5", "@SQ\tSN:chr1\tLN:100"};
    CHECK(records("@HD\tVN:1.5\n@SQ\tSN:chr1\tLN:100\n") == want);
    CHECK(records("@HD\tVN:1.5\n@SQ\tSN:chr1\tLN:100")   == want);

    // A single record, either way.
    CHECK(records("@HD\tVN:1.5")   == std::vector<std::string>{"@HD\tVN:1.5"});
    CHECK(records("@HD\tVN:1.5\n") == std::vector<std::string>{"@HD\tVN:1.5"});

    // Empty records ARE reported; callers that care skip them (both meth_bam
    // helpers do, via `len > 0`). The iterator does not silently drop data.
    CHECK(records("a\n\nb") == std::vector<std::string>({"a", "", "b"}));

    // Degenerate inputs terminate rather than loop or deref.
    CHECK(records("").empty());
    const char *nul = NULL, *line; size_t len;
    CHECK(bwa_hdr_next_line(&nul, &line, &len) == 0);
    CHECK(bwa_hdr_next_line(NULL, &line, &len) == 0);
}

TEST_CASE("bwa_hdr_next_line advances exactly past each record") {
    // The cursor must land on the next record's first byte, so a caller can
    // stop early and resume -- meth_append_sq_extra_tags returns mid-iteration.
    const char *text = "one\ntwo\nthree";
    const char *p = text, *line; size_t len;
    REQUIRE(bwa_hdr_next_line(&p, &line, &len) == 1);
    CHECK(std::string(line, len) == "one");
    CHECK(std::string(p) == "two\nthree");
    REQUIRE(bwa_hdr_next_line(&p, &line, &len) == 1);
    CHECK(std::string(line, len) == "two");
    CHECK(std::string(p) == "three");
    REQUIRE(bwa_hdr_next_line(&p, &line, &len) == 1);
    CHECK(std::string(line, len) == "three");
    CHECK(*p == '\0');
    CHECK(bwa_hdr_next_line(&p, &line, &len) == 0);
}
