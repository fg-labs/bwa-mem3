// test/unit/test_temp_index.cpp
//
// Self-checks for the scratch-file helpers in test/framework/ that other unit
// tests build on: TempDir and the temp_path.h calls create their files under
// $TMPDIR, TempDir removes what it made, StderrCapture and QuietStderr
// restore stderr, and TempIndex builds and loads an index whose contig table,
// ambiguity runs and packed sequence match the contigs it was given. A helper
// that drifted would show up here rather than as a wrong expectation elsewhere.

#include "doctest/doctest.h"
#include "run_in_child.h"
#include "scoped_env.h"
#include "stderr_capture.h"
#include "temp_dir.h"
#include "temp_index.h"

#include <dirent.h>
#include <sys/stat.h>
#include <sys/wait.h>
#include <unistd.h>

#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <memory>
#include <random>
#include <string>
#include <vector>

using bwa_tests::Contig;

namespace {

bool exists(const std::string &path) {
    struct stat st;
    return stat(path.c_str(), &st) == 0;
}

// The file stderr refers to, so a case can tell whether it was restored.
struct FileId {
    dev_t dev;
    ino_t ino;
    bool operator==(const FileId &o) const { return dev == o.dev && ino == o.ino; }
};

FileId stderr_id() {
    struct stat st;
    REQUIRE(fstat(fileno(stderr), &st) == 0);
    return {st.st_dev, st.st_ino};
}

std::string slurp(const std::string &path) {
    std::string out;
    if (FILE *fp = fopen(path.c_str(), "rb")) {
        char buf[4096];
        size_t got;
        while ((got = fread(buf, 1, sizeof(buf), fp)) > 0) out.append(buf, got);
        fclose(fp);
    }
    return out;
}

std::string random_bases(std::mt19937 &rng, size_t n) {
    std::string s(n, 'A');
    for (char &c : s) c = "ACGT"[rng() & 3];
    return s;
}

// The number of BWT rows whose suffix-array value, as the index resolves it,
// differs from a suffix array sorted by brute force from the contigs: forward
// strand then reverse complement, an N taking the base pac stored for it.
// Resolving a row reads the BWT (the LF walk from an unsampled row) and the SA
// samples. Row 0 is the sentinel's: the empty suffix, at the text length.
int64_t sa_lookup_mismatches(const bwa_tests::TempIndex &idx, const std::vector<Contig> &ref) {
    const int64_t l_pac = idx.bns()->l_pac;
    std::vector<uint8_t> text;
    for (const Contig &c : ref) {
        const std::vector<uint8_t> enc = bwa_tests::nt4(c.seq);
        text.insert(text.end(), enc.begin(), enc.end());
    }
    for (int64_t p = 0; p < l_pac; ++p)
        if (text[(size_t)p] > 3) text[(size_t)p] = (idx.pac()[p >> 2] >> ((~p & 3) << 1)) & 3;
    for (int64_t p = l_pac - 1; p >= 0; --p) text.push_back((uint8_t)(3 - text[(size_t)p]));
    std::vector<int64_t> sa(text.size());
    for (size_t i = 0; i < sa.size(); ++i) sa[i] = (int64_t)i;
    std::sort(sa.begin(), sa.end(), [&text](int64_t a, int64_t b) {
        return std::lexicographical_compare(text.begin() + a, text.end(), text.begin() + b, text.end());
    });
    int64_t bad = idx.fmi()->get_sa_entry_compressed(0) != 2 * l_pac;
    for (size_t r = 0; r < sa.size(); ++r)
        bad += idx.fmi()->get_sa_entry_compressed((int64_t)r + 1) != sa[r];
    return bad;
}

} // namespace

TEST_CASE("TempDir: makes a fresh directory per object and removes its files"
          * doctest::test_suite("unit/test_framework")) {
    std::string path, file;
    {
        bwa_tests::TempDir a, b;
        CHECK(a.path() != b.path());
        a.write("x.txt", "hello\n");
        path = a.path();
        file = a.file("x.txt");
        CHECK(exists(file));
    }
    CHECK_FALSE(exists(file));
    CHECK_FALSE(exists(path));
}

TEST_CASE("TempDir: a forked child that destroys one leaves the parent's directory"
          * doctest::test_suite("unit/test_framework")) {
    // A child of run_in_child that exits through exit() (err_fatal, say) runs
    // the static destructors it inherited, so a long-lived TempDir (a static
    // test fixture) is destroyed in the child too; it must not remove the
    // directory the parent still uses.
    std::unique_ptr<bwa_tests::TempDir> dir(new bwa_tests::TempDir);
    dir->write("x.txt", "hello\n");
    const std::string file = dir->file("x.txt");
    const int status = bwa_tests::run_in_child([&dir] {
        dir.reset();
        _exit(0);
    });
    CHECK((WIFEXITED(status) && WEXITSTATUS(status) == 0));
    CHECK(exists(file));
    const std::string path = dir->path();
    dir.reset();   // the parent's own destructor still removes it
    CHECK_FALSE(exists(path));
}

TEST_CASE("TempDir: remove() deletes the directory now, and only once"
          * doctest::test_suite("unit/test_framework")) {
    std::string path, file, err;
    {
        bwa_tests::StderrCapture cap;
        {
            bwa_tests::TempDir dir;
            dir.write("x.txt", "hello\n");
            path = dir.path();
            file = dir.file("x.txt");
            dir.remove();
            CHECK_FALSE(exists(file));
            CHECK_FALSE(exists(path));
            CHECK(dir.path().empty());
            dir.remove();   // a second call does nothing
        }                   // nor does the destructor: no warning
        err = cap.text();
    }
    CHECK(err == "");
}

TEST_CASE("TempDir and temp_path.h: scratch paths go under $TMPDIR"
          * doctest::test_suite("unit/test_framework")) {
    bwa_tests::TempDir root("tmpdir_root");
    auto under = [&](const std::string &path, const char *prefix) {
        const std::string want = root.path() + "/" + prefix + ".";
        return path.compare(0, want.size(), want) == 0 && path.size() == want.size() + 6;
    };
    // With and without a trailing slash, every helper names "<TMPDIR>/<prefix>.XXXXXX".
    for (const std::string tmpdir : {root.path(), root.path() + "/"}) {
        CAPTURE(tmpdir);
        bwa_tests::ScopedEnv env("TMPDIR", tmpdir.c_str());
        CHECK(bwa_tests::temp_root() == root.path() + "/");
        {
            bwa_tests::TempDir d("td");
            CHECK(under(d.path(), "td"));
        }
        char buf[BWA_TEST_PATH_MAX];
        REQUIRE(bwa_test_temp_template(buf, sizeof buf, "tt") == 0);
        CHECK(std::string(buf) == root.path() + "/tt.XXXXXX");
        const int fd = bwa_test_mkstemp(buf, sizeof buf, "ms");
        REQUIRE(fd >= 0);
        close(fd);
        CHECK(under(buf, "ms"));
        CHECK(unlink(buf) == 0);
        REQUIRE(bwa_test_mkdtemp(buf, sizeof buf, "md") == 0);
        CHECK(under(buf, "md"));
        CHECK(bwa_test_remove_flat_dir(buf) == 0);
        const int anon = bwa_tests::anon_temp_fd("af");
        CHECK(anon >= 0);
        if (anon >= 0) close(anon);
    }
    // Nothing is left in the root: TempDir removed its own directory,
    // anon_temp_fd unlinked its file, and the loop removed the rest.
    int left = 0;
    if (DIR *d = opendir(root.path().c_str())) {
        while (struct dirent *e = readdir(d))
            if (strcmp(e->d_name, ".") != 0 && strcmp(e->d_name, "..") != 0) ++left;
        closedir(d);
    }
    CHECK(left == 0);

    // A template that does not fit is refused rather than truncated.
    char small[8];
    {
        bwa_tests::QuietStderr quiet;
        CHECK(bwa_test_temp_template(small, sizeof small, "much_too_long") == -1);
    }
}

TEST_CASE("StderrCapture: returns what was written and restores stderr"
          * doctest::test_suite("unit/test_framework")) {
    const FileId before = stderr_id();
    {
        bwa_tests::StderrCapture capture;
        CHECK_FALSE(stderr_id() == before);
        fprintf(stderr, "first %d\n", 1);
        fputs("second\n", stderr);
        CHECK(capture.text() == "first 1\nsecond\n");
        CHECK(stderr_id() == before);   // text() restores
    }
    CHECK(stderr_id() == before);
    {
        // stderr is /dev/null inside (compared with /dev/null itself, since a
        // runner may already have pointed stderr there).
        struct stat null_st;
        REQUIRE(stat("/dev/null", &null_st) == 0);
        bwa_tests::QuietStderr quiet;
        CHECK(stderr_id() == FileId{null_st.st_dev, null_st.st_ino});
    }
    CHECK(stderr_id() == before);
    // A second capture after the first is released sees only its own output.
    bwa_tests::StderrCapture again;
    {
        bwa_tests::QuietStderr quiet;
        fputs("discarded\n", stderr);
    }
    fputs("kept\n", stderr);
    CHECK(again.text() == "kept\n");
}

TEST_CASE("TempIndex: builds and loads the contigs it was given"
          * doctest::test_suite("unit/test_framework")) {
    std::mt19937 rng(20261029);
    std::string c2 = random_bases(rng, 900);
    c2.replace(300, 40, std::string(40, 'N'));
    const std::vector<Contig> ref = {
        {"chrA", random_bases(rng, 1200)}, {"chrB", c2}, {"chrC", random_bases(rng, 150)}};

    bwa_tests::TempIndex idx(ref);
    const bntseq_t *bns = idx.bns();
    REQUIRE(bns != NULL);
    REQUIRE(bns->n_seqs == 3);
    int64_t total = 0;
    for (int i = 0; i < bns->n_seqs; ++i) {
        CHECK(std::string(bns->anns[i].name) == ref[(size_t)i].name);
        CHECK(bns->anns[i].len == (int64_t)ref[(size_t)i].seq.size());
        CHECK(bns->anns[i].offset == total);
        total += bns->anns[i].len;
    }
    CHECK(bns->l_pac == total);
    // The 40 N of chrB are the one ambiguity run.
    REQUIRE(bns->n_holes == 1);
    CHECK(bns->ambs[0].offset == 1200 + 300);
    CHECK(bns->ambs[0].len == 40);
    CHECK(bns->ambs[0].amb == 'N');
    // pac holds the forward strand of every contig back to back, two bits per
    // base, first base in the high bits. An N is stored as a random base, so
    // only ACGT positions are compared.
    REQUIRE(idx.pac() != NULL);
    int64_t pos = 0, compared = 0;
    for (const Contig &c : ref) {
        const std::vector<uint8_t> want = bwa_tests::nt4(c.seq);
        for (size_t j = 0; j < want.size(); ++j, ++pos) {
            if (want[j] > 3) continue;
            CAPTURE(c.name);
            CAPTURE(j);
            CHECK(((idx.pac()[pos >> 2] >> ((~pos & 3) << 1)) & 3) == want[j]);
            ++compared;
        }
    }
    CHECK(pos == bns->l_pac);
    CHECK(compared == bns->l_pac - 40);
    // The FM index: its BWT covers both strands plus the sentinel, and it was
    // sampled at the requested SA rate.
    REQUIRE(idx.fmi() != NULL);
    CHECK(idx.fmi()->reference_seq_len == 2 * bns->l_pac + 1);
    CHECK(idx.fmi()->sa_compx == 3);
    // Every BWT row, the sentinel's included, resolves to its suffix-array value.
    CHECK(sa_lookup_mismatches(idx, ref) == 0);
    // A fixed FASTA name would let a staged shm index of that name be loaded
    // instead; the prefix is the TempDir's unique leaf.
    const std::string &prefix = idx.prefix();
    const std::string base = prefix.substr(prefix.rfind('/') + 1);
    CHECK(base != "ref.fa");
    CHECK(prefix.find("/" + base.substr(0, base.size() - 3) + "/") != std::string::npos);
}

TEST_CASE("TempIndex: honours the SA sampling shift"
          * doctest::test_suite("unit/test_framework")) {
    std::mt19937 rng(7);
    const std::vector<Contig> ref = {{"chrA", random_bases(rng, 500)}};
    bwa_tests::TempIndex idx(ref, /*sa_compx=*/1);
    CHECK(idx.fmi()->sa_compx == 1);
    CHECK(sa_lookup_mismatches(idx, ref) == 0);
}

TEST_CASE("TempIndex: drop_files deletes the index files and the loaded index still answers"
          * doctest::test_suite("unit/test_framework")) {
    std::mt19937 rng(31);
    const std::vector<Contig> ref = {{"chrA", random_bases(rng, 700)}, {"chrB", random_bases(rng, 300)}};
    bwa_tests::TempIndex idx(ref);
    const std::string prefix = idx.prefix();
    const std::string dir = prefix.substr(0, prefix.rfind('/'));
    REQUIRE(exists(prefix));
    idx.drop_files();
    CHECK_FALSE(exists(prefix));
    CHECK_FALSE(exists(dir));
    REQUIRE(idx.bns() != NULL);
    CHECK(idx.bns()->n_seqs == 2);
    CHECK(idx.bns()->l_pac == 1000);
    CHECK(sa_lookup_mismatches(idx, ref) == 0);
}

TEST_CASE("TempIndex: a build that fails and returns throws with its log"
          * doctest::test_suite("unit/test_framework")) {
    // An empty reference: bwa_idx_build reports a malformed .ann and returns
    // nonzero. Run in a child so that a regression to exit() cannot take the
    // runner down.
    bwa_tests::TempDir dir;
    const std::string err = dir.file("child.err");
    const int status = bwa_tests::run_in_child([] {
        try {
            bwa_tests::TempIndex idx({});
        } catch (const std::runtime_error &e) {
            fprintf(stderr, "THROWN: %s\n", e.what());
            _exit(3);
        }
    }, err.c_str());
    const std::string text = slurp(err);
    CAPTURE(text);
    CHECK((WIFEXITED(status) && WEXITSTATUS(status) == 3));
    const size_t thrown = text.find("THROWN: bwa_idx_build failed for ");
    REQUIRE(thrown != std::string::npos);
    // The build's own log, carried in the exception: the builder also writes
    // it to stderr, so only the text after THROWN: shows it was captured.
    CHECK(text.find("malformed", thrown) != std::string::npos);
}

TEST_CASE("LoggedStderr (TempIndex's build log): a process that exits has the log replayed to the real stderr"
          * doctest::test_suite("unit/test_framework")) {
    bwa_tests::TempDir dir;
    const std::string err = dir.file("child.err"), log = dir.file("index.log");
    const int status = bwa_tests::run_in_child([&log] {
        bwa_tests::temp_index_detail::LoggedStderr logged(log);
        fputs("fatal: simulated index failure\n", stderr);
        exit(1);
    }, err.c_str());
    const std::string text = slurp(err);
    CAPTURE(text);
    CHECK((WIFEXITED(status) && WEXITSTATUS(status) == 1));
    CHECK(text.find("TempIndex: index build or load exited the process; its log follows:\n"
                    "fatal: simulated index failure\n") != std::string::npos);
}

TEST_CASE("revcomp and nt4: the encodings mem aligns in"
          * doctest::test_suite("unit/test_framework")) {
    CHECK(bwa_tests::revcomp("AACGTN") == "NACGTT");
    CHECK(bwa_tests::revcomp("acgtR") == "Nacgt");
    CHECK(bwa_tests::nt4("ACGTNacgtR") == std::vector<uint8_t>{0, 1, 2, 3, 4, 0, 1, 2, 3, 4});
}

TEST_CASE("fasta_text: wraps at the given width and refuses a zero width"
          * doctest::test_suite("unit/test_framework")) {
    const std::vector<Contig> contigs = {{"a", "ACGTA"}, {"b", "CC"}};
    CHECK(bwa_tests::fasta_text(contigs, 2) == ">a\nAC\nGT\nA\n>b\nCC\n");
    CHECK_THROWS_AS(bwa_tests::fasta_text(contigs, 0), std::invalid_argument);
}
