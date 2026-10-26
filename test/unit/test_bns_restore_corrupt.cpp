// Unit tests for bns_restore_core's handling of malformed .ann/.amb headers and
// of long .ann comment lines.
//
// The counts and lengths in the .ann ("l_pac n_seqs seed") and .amb ("l_pac
// n_seqs n_holes") headers come straight from files on disk, so a value out of
// range, or two headers that disagree, is a corrupt index rather than a program
// bug: it has to fail with an error in every build, including one compiled with
// -DNDEBUG, rather than reach calloc as a huge size or wrap to a small count.
// Conversely, a comment line that fills the read buffer is valid input that is
// loaded whole, and a longer one is truncated with a warning.
//
// Each case runs bns_restore_core in a forked child (run_in_child.h), because
// err_fatal ends the process. A child that loads the index and passes its
// checks exits 0; err_fatal exits with EXIT_FAILURE; an assert or xassert
// aborts on SIGABRT. The parse-error path (badread) also exits with
// EXIT_FAILURE, so a rejection is identified by its message on the child's
// stderr, not by the status alone.

#include "doctest/doctest.h"
#include "run_in_child.h"

#include "bntseq.h"

#include <climits>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <sstream>
#include <string>
#include <sys/wait.h>
#include <unistd.h>

namespace {

/* A temporary index prefix holding the given .ann and .amb contents and an
 * empty .pac, plus a file for the child's stderr. Removed on destruction. */
struct IndexFixture {
    char dir[PATH_MAX];
    std::string ann, amb, pac, err;

    IndexFixture(const std::string &ann_body, const std::string &amb_body) {
        const char *tmpdir = getenv("TMPDIR");
        if (tmpdir == NULL || tmpdir[0] == '\0') tmpdir = "/tmp";
        const int n = snprintf(dir, sizeof(dir), "%s%sbns_restore_corrupt_XXXXXX",
                               tmpdir, tmpdir[strlen(tmpdir) - 1] == '/' ? "" : "/");
        REQUIRE(n > 0);
        REQUIRE((size_t)n < sizeof(dir));
        REQUIRE(mkdtemp(dir) != NULL);
        ann = std::string(dir) + "/ref.ann";
        amb = std::string(dir) + "/ref.amb";
        pac = std::string(dir) + "/ref.pac";
        err = std::string(dir) + "/stderr.txt";
        write_file(ann, ann_body);
        write_file(amb, amb_body);
        write_file(pac, "");
    }

    ~IndexFixture() {
        unlink(ann.c_str());
        unlink(amb.c_str());
        unlink(pac.c_str());
        unlink(err.c_str());
        rmdir(dir);
    }

    IndexFixture(const IndexFixture &)            = delete;
    IndexFixture &operator=(const IndexFixture &) = delete;

    static void write_file(const std::string &path, const std::string &body) {
        FILE *fp = fopen(path.c_str(), "w");
        REQUIRE(fp != NULL);
        REQUIRE(fwrite(body.data(), 1, body.size(), fp) == body.size());
        REQUIRE(fclose(fp) == 0);
    }

    /* What the last load wrote to stderr. */
    std::string child_stderr() const {
        std::ifstream in(err);
        std::ostringstream text;
        text << in.rdbuf();
        return text.str();
    }

    /* Load the index in a forked child and return its wait status. The child
     * exits with check(bns) after a successful load. */
    int load(int (*check)(const bntseq_t *)) const {
        return bwa_tests::run_in_child([this, check] {
            bntseq_t *bns = bns_restore_core(ann.c_str(), amb.c_str(), pac.c_str());
            const int rc = check(bns);
            bns_destroy(bns);
            _exit(rc);
        }, err.c_str());
    }
};

int never_reached(const bntseq_t *) { return 3; }

/* Whatever the comment's length, the anno keeps at most the first 8191 of its
 * characters (the 8192-byte read limit less the leading space), and the record
 * after it parses. */
int one_contig_8191_char_anno(const bntseq_t *bns) {
    if (bns == NULL || bns->n_seqs != 1) return 4;
    if (strcmp(bns->anns[0].name, "chr1") != 0) return 5;
    if (strlen(bns->anns[0].anno) != 8191) return 6;
    if (bns->anns[0].len != 10 || bns->anns[0].offset != 0) return 7;
    return 0;
}

std::string ann_with_comment(size_t comment_len) {
    return "10 1 11\n0 chr1 " + std::string(comment_len, 'x') + "\n0 10 0\n";
}

}  // namespace

TEST_SUITE("unit/bns_restore") {

TEST_CASE("bns_restore_core rejects an out-of-range or inconsistent header") {
    struct Case {
        const char *what, *ann, *amb;
        bool in_amb;               // which file the error names
        const char *message;       // follows "corrupt index <file>: "
    };
    const Case cases[] = {
        {"negative reference length", "-10 0 11\n", "-10 0 0\n", false,
         "reference length -10 is negative"},
        {"negative sequence count", "0 -1 11\n", "0 -1 0\n", false,
         "sequence count -1 is out of range"},
        {"sequence count above INT_MAX", "0 2147483648 11\n", "0 2147483648 0\n", false,
         "sequence count 2147483648 is out of range"},
        {"negative hole count", "0 0 11\n", "0 0 -1\n", true,
         "hole count -1 is out of range"},
        {"hole count above INT_MAX", "0 0 11\n", "0 0 4294967296\n", true,
         "hole count 4294967296 is out of range"},
        {"headers disagree on the reference length", "0 0 11\n", "5 0 0\n", true,
         "header (5 0) does not match the .ann header (0 0)"},
        {"headers disagree on the sequence count", "10 1 11\n0 chr1 (null)\n0 10 0\n",
         "10 2 0\n", true, "header (10 2) does not match the .ann header (10 1)"},
    };
    for (const Case &c : cases) {
        CAPTURE(c.what);
        IndexFixture fx(c.ann, c.amb);
        const int status = fx.load(never_reached);
        REQUIRE(WIFEXITED(status));
        CHECK(WEXITSTATUS(status) == EXIT_FAILURE);
        const std::string expected =
            "corrupt index " + (c.in_amb ? fx.amb : fx.ann) + ": " + c.message;
        CHECK(fx.child_stderr().find(expected) != std::string::npos);
    }
}

TEST_CASE("bns_restore_core loads an .ann comment that exactly fills the read buffer") {
    IndexFixture fx(ann_with_comment(8191), "10 1 0\n");
    const int status = fx.load(one_contig_8191_char_anno);
    REQUIRE(WIFEXITED(status));
    CHECK(WEXITSTATUS(status) == 0);
    CHECK(fx.child_stderr().find("truncated") == std::string::npos);
}

TEST_CASE("bns_restore_core truncates a longer .ann comment and warns") {
    IndexFixture fx(ann_with_comment(9000), "10 1 0\n");
    const int status = fx.load(one_contig_8191_char_anno);
    REQUIRE(WIFEXITED(status));
    CHECK(WEXITSTATUS(status) == 0);
    CHECK(fx.child_stderr().find(fx.ann + ": the comment of 'chr1' is longer than 8191 "
                                 "characters and was truncated") != std::string::npos);
}

}  // TEST_SUITE
