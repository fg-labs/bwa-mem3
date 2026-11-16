// test/unit/test_meth_hdr_text.cpp -- @HD placement in the --meth header text
// (meth_hdr_text, src/meth_bam.cpp).
//
// The SAM spec requires @HD, when present, to be the first header line.
// meth_hdr_text appends the user's header text (-R/-r read group, -H lines)
// after the @SQ block, so a user-supplied @HD used to land after @SQ. `mem
// --meth` never showed it, because it hands the text to sam_hdr_add_lines and
// htslib moves a lone @HD to the front; a caller that writes the text verbatim
// (the planned `aln --meth`) would emit it out of order. These tests pin the
// text itself, so they hold for every caller.
//
// A user @HD is treated as the default (non-meth) header builder,
// bwa_build_sam_hdr, treats it: the user's first @HD wins over the default and
// is the first line; every other @HD is dropped.

// htslib headers first: they share the KSTRING_H include guard with
// src/kstring.h, which bwa.h pulls in.
#include "htslib/kstring.h"

#include "doctest/doctest.h"

#include "bwa.h"
#include "compat_target.h"
#include "meth_bam.h"

#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

namespace {

// One plain contig, so the @SQ block is a single known line.
struct TestBns {
    bntann1_t ann;
    bntseq_t  bns;
    TestBns() {
        memset(&ann, 0, sizeof(ann));
        ann.name = const_cast<char *>("chr1");
        ann.len  = 1000;
        memset(&bns, 0, sizeof(bns));
        bns.n_seqs = 1;
        bns.anns   = &ann;
    }
};

// The --meth header text for user header text `hdr_line` (NULL for none),
// with the production default @HD and no sidecar.
std::string render(const char *hdr_line) {
    TestBns t;
    kstring_t out = {0, 0, NULL};
    const int rc = meth_hdr_text(&out, &t.bns, BWAMEM3_DEFAULT_HD_LINE, hdr_line, NULL,
                                 "@PG\tID:bwa-mem3\tPN:bwa-mem3", NULL, "bwa-mem3 test");
    REQUIRE(rc == 0);
    std::string s = out.s ? std::string(out.s, out.l) : std::string();
    free(out.s);
    return s;
}

// The text split into lines (each newline-terminated in the input).
std::vector<std::string> lines_of(const std::string &s) {
    std::vector<std::string> v;
    size_t pos = 0;
    while (pos < s.size()) {
        size_t eol = s.find('\n', pos);
        if (eol == std::string::npos) eol = s.size();
        v.push_back(s.substr(pos, eol - pos));
        pos = eol + 1;
    }
    return v;
}

int count_type(const std::vector<std::string> &lines, const char *type) {
    int n = 0;
    for (const auto &l : lines)
        if (l.compare(0, strlen(type), type) == 0) ++n;
    return n;
}

} // namespace

TEST_CASE("meth_hdr_text: with no user @HD the default @HD is the first line"
          * doctest::test_suite("unit/meth_hdr_text")) {
    for (const char *hdr_line : {(const char *)NULL, "", "@RG\tID:g\tSM:s"}) {
        CAPTURE(hdr_line ? hdr_line : "(null)");
        const auto lines = lines_of(render(hdr_line));
        REQUIRE(!lines.empty());
        CHECK(lines[0] == BWAMEM3_DEFAULT_HD_LINE);
        CHECK(count_type(lines, "@HD\t") == 1);
    }
}

TEST_CASE("meth_hdr_text: a user @HD replaces the default and is the first line"
          * doctest::test_suite("unit/meth_hdr_text")) {
    const std::string user_hd = "@HD\tVN:1.6\tSO:unsorted";
    struct Case { const char *desc; std::string hdr_line; };
    const std::vector<Case> cases = {
        { "leading @HD",            user_hd + "\n@RG\tID:g\tSM:new" },
        { "non-leading @HD",        "@RG\tID:g\tSM:new\n" + user_hd },
        { "@HD between records",    "@RG\tID:g\tSM:new\n" + user_hd + "\n@CO\tnote" },
        { "@HD alone",              user_hd },
        { "trailing newline",       "@RG\tID:g\tSM:new\n" + user_hd + "\n" },
    };
    for (const auto &c : cases) {
        CAPTURE(c.desc);
        const auto lines = lines_of(render(c.hdr_line.c_str()));
        REQUIRE(!lines.empty());
        CHECK(lines[0] == user_hd);
        CHECK(count_type(lines, "@HD\t") == 1);
        // The rest of the user text is still emitted, after the @SQ block.
        if (c.hdr_line.find("@RG\t") != std::string::npos) {
            CHECK(count_type(lines, "@RG\t") == 1);
            CHECK(lines[1] == "@SQ\tSN:chr1\tLN:1000");
        }
    }
}

TEST_CASE("meth_hdr_text: of several user @HD records the first wins, the rest are dropped"
          * doctest::test_suite("unit/meth_hdr_text")) {
    const auto lines = lines_of(render("@RG\tID:g\tSM:s\n@HD\tVN:1.6\tSO:coordinate\n"
                                       "@CO\tc\n@HD\tVN:1.5\tSO:unsorted"));
    REQUIRE(!lines.empty());
    CHECK(lines[0] == "@HD\tVN:1.6\tSO:coordinate");
    CHECK(count_type(lines, "@HD\t") == 1);
    CHECK(count_type(lines, "@RG\t") == 1);
    CHECK(count_type(lines, "@CO\t") == 1);
}

TEST_CASE("meth_hdr_text: the whole text with a hoisted user @HD"
          * doctest::test_suite("unit/meth_hdr_text")) {
    // Every record, in order: the user's @HD, @SQ, the rest of the user text in its own
    // order (the second @HD and the empty line dropped), then the two @PG lines.
    const std::string text =
        render("@RG\tID:a\tSM:s\n@HD\tVN:1.6\n\n@RG\tID:b\tSM:s\n@HD\tVN:1.5\n@CO\tlast");
    // The meth @PG carries the build's version string; everything around it is pinned.
    const std::string head = "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:1000\n"
                             "@RG\tID:a\tSM:s\n@RG\tID:b\tSM:s\n@CO\tlast\n"
                             "@PG\tID:bwa-mem3\tPN:bwa-mem3\n"
                             "@PG\tID:bwa-mem3-meth\tPN:bwa-mem3-meth\tVN:";
    const std::string tail = "-meth\tCL:bwa-mem3 test\n";
    REQUIRE(text.size() > head.size() + tail.size());
    CHECK(text.compare(0, head.size(), head) == 0);
    CHECK(text.compare(text.size() - tail.size(), tail.size(), tail) == 0);
    // Nothing between the two: the version is one token on the last line.
    const std::string version = text.substr(head.size(), text.size() - head.size() - tail.size());
    CHECK(version.find_first_of("\t\n") == std::string::npos);
}

TEST_CASE("meth_hdr_text: a NULL default @HD with no user @HD emits no @HD"
          * doctest::test_suite("unit/meth_hdr_text")) {
    TestBns t;
    kstring_t out = {0, 0, NULL};
    REQUIRE(meth_hdr_text(&out, &t.bns, NULL, "@RG\tID:g\tSM:s", NULL, NULL, NULL, "x") == 0);
    const auto lines = lines_of(std::string(out.s, out.l));
    free(out.s);
    CHECK(count_type(lines, "@HD\t") == 0);
    REQUIRE(!lines.empty());
    CHECK(lines[0] == "@SQ\tSN:chr1\tLN:1000");
}
