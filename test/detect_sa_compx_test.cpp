// Regression test for detect_sa_compx() (src/FMI_search.cpp).
//
// detect_sa_compx() reads the trailing 8-byte SA-rate tag from an FM-index
// file and verifies it reproduces the file's on-disk size via the layout
// formula. It returns default_compx only for a legacy (pre-tag) index or a
// file too small to carry a tail.
//
// The behavior this test pins down: a *read failure* of the tail (I/O error,
// short/EOF read) must NOT be silently converted to default_compx -- doing so
// would mis-size the SA-sample arrays of a genuinely non-default index with no
// diagnostic. Such a failure must abort the load (EXIT_FAILURE). Only a
// fully-read tag that fails the layout test falls back to default_compx.
//
// Usage:
//   detect_sa_compx_test          # runs all cases, exits 0 on success

#include "FMI_search.h"
#include "bwa.h"          // bwa_idx_build(): the real `index -u` build entry point

#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <fcntl.h>        // open/O_RDONLY
#include <sys/stat.h>     // fstat/struct stat
#include <sys/wait.h>
#include <unistd.h>

#define CHECK(cond) do { \
    if (!(cond)) { \
        std::fprintf(stderr, "CHECK failed at %s:%d: %s\n", __FILE__, __LINE__, #cond); \
        std::fflush(stderr); \
        std::abort(); \
    } \
} while (0)

// Mirror detect_sa_compx()'s off_sent_for(): the on-disk offset of the
// sentinel index for a given SA sample rate and reference length.
static int64_t off_sent_for(int64_t ref_seq_len, int64_t compx)
{
    const int64_t hdr_bytes     = 6 * (int64_t)sizeof(int64_t);
    const int64_t cp_occ_cnt    = (ref_seq_len >> CP_SHIFT) + 1;
    const int64_t sa_sample_cnt = (ref_seq_len >> compx) + 1;
    return hdr_bytes + cp_occ_cnt * (int64_t)sizeof(CP_OCC)
                     + sa_sample_cnt * (int64_t)sizeof(int8_t)
                     + sa_sample_cnt * (int64_t)sizeof(uint32_t);
}

// Create a file of exactly `size` bytes with its trailing 8 bytes set to
// `tail`. Returns an open fd (caller closes) via the path in `path_out`.
static int make_file(const char *tmpl_in, int64_t size, int64_t tail,
                     std::string &path_out)
{
    char tmpl[64];
    std::snprintf(tmpl, sizeof(tmpl), "%s", tmpl_in);
    int fd = mkstemp(tmpl);
    CHECK(fd >= 0);
    path_out = tmpl;
    CHECK(ftruncate(fd, (off_t)size) == 0);
    if (size >= (int64_t)sizeof(int64_t)) {
        ssize_t w = pwrite(fd, &tail, sizeof(int64_t),
                           (off_t)(size - (int64_t)sizeof(int64_t)));
        CHECK(w == (ssize_t)sizeof(int64_t));
    }
    return fd;
}

int main()
{
    const int64_t ref_seq_len = 100000;

    // Case 1 (independent oracle): build a REAL index at a non-default SA rate
    // via the production entry point bwa_idx_build() -- the same call chain
    // `bwamem3 index -u` uses (bns_fasta2bntseq -> FMI_search::build_index ->
    // write_fm_index_streaming). detect_sa_compx() then recovers the rate from
    // the file the production writer actually emitted, so a shared offset or
    // byte-count error between the writer's on-disk layout and detect_sa_compx()
    // can no longer agree with a duplicated formula in this test and pass. The
    // off_sent_for() helper above is retained only for the fallback-branch
    // cases 2a/2b below, which exercise control flow, not layout correctness.
    {
        const int chosen_c = 4;  // non-default rate (compile-time default is 3); in [0, 6]

        char dir[] = "/tmp/detect_sa_compx_real_XXXXXX";
        CHECK(mkdtemp(dir) != nullptr);
        const std::string fa_path = std::string(dir) + "/ref.fa";
        const std::string prefix  = std::string(dir) + "/idx";

        // A deterministic ~1 kbp reference (fixed-seed LCG over ACGT) so the
        // build is reproducible and large enough for a well-formed layout.
        {
            FILE *fa = std::fopen(fa_path.c_str(), "w");
            CHECK(fa != nullptr);
            std::fputs(">chr1\n", fa);
            uint64_t s = 0x9e3779b97f4a7c15ULL;  // fixed seed -> reproducible bases
            const char bases[4] = {'A', 'C', 'G', 'T'};
            for (int i = 0; i < 1024; i++) {
                s = s * 6364136223846793005ULL + 1442695040888963407ULL;
                std::fputc(bases[(s >> 33) & 3u], fa);
            }
            std::fputc('\n', fa);
            CHECK(std::fclose(fa) == 0);
        }

        // Build the real FM index at the chosen rate (emit_unpacked_ref = 0).
        CHECK(bwa_idx_build(fa_path.c_str(), prefix.c_str(), 0, chosen_c) == 0);

        const std::string bwt_path = prefix + ".bwt.2bit.64";
        int fd = open(bwt_path.c_str(), O_RDONLY);
        CHECK(fd >= 0);
        struct stat st;
        CHECK(fstat(fd, &st) == 0);

        // reference_seq_len is the file's first header field (offset 0), written
        // verbatim by write_fm_index_streaming -- a pass-through, not part of the
        // disputed layout formula -- so reading it back is not circular. It is
        // exactly the value the production loader feeds to detect_sa_compx().
        int64_t file_ref_seq_len = 0;
        CHECK(pread(fd, &file_ref_seq_len, sizeof(file_ref_seq_len), 0)
              == (ssize_t)sizeof(file_ref_seq_len));

        int64_t got = detect_sa_compx(fd, (int64_t)st.st_size, file_ref_seq_len, SA_COMPX);
        close(fd);
        CHECK(got == chosen_c);

        // Remove every sidecar the writer may have produced, then the tempdir.
        for (const char *ext : {".pac", ".ann", ".amb", ".bwt.2bit.64", ".0123", ".sa"}) {
            const std::string p = prefix + ext;
            unlink(p.c_str());
        }
        unlink(fa_path.c_str());
        rmdir(dir);
    }

    // Case 2a: a legacy-shaped index (no tag: size == off_sent_for(default) +
    // sentinel only) whose trailing 8 bytes are the sentinel index falls back
    // to the default rate.
    {
        const int64_t size = off_sent_for(ref_seq_len, SA_COMPX) + (int64_t)sizeof(int64_t);
        std::string path;
        // Trailing value is a sentinel index, not a valid rate/layout match.
        int fd = make_file("/tmp/detect_sa_compx_legacy_XXXXXX", size, 123456, path);
        int64_t got = detect_sa_compx(fd, size, ref_seq_len, SA_COMPX);
        close(fd);
        unlink(path.c_str());
        CHECK(got == SA_COMPX);
    }

    // Case 2b: a file too small to carry a tail (< 2*int64) falls back to the
    // default without attempting a read.
    {
        std::string path;
        int fd = make_file("/tmp/detect_sa_compx_tiny_XXXXXX", 8, 0, path);
        int64_t got = detect_sa_compx(fd, 8, ref_seq_len, SA_COMPX);
        close(fd);
        unlink(path.c_str());
        CHECK(got == SA_COMPX);
    }

    // Case 3 (the fix): a tail-read failure must abort the load rather than
    // silently returning the default. An invalid fd makes pread() fail with
    // EBADF (not EINTR); the call must exit(EXIT_FAILURE). Before the fix,
    // pread's error was folded into the legacy default and the process
    // returned normally.
    {
        const int fd = -1;  // invalid descriptor -> pread returns EBADF

        pid_t pid = fork();
        CHECK(pid >= 0);
        if (pid == 0) {
            // Child: file_size >= 16 forces the tail-read branch.
            (void)detect_sa_compx(fd, 16, ref_seq_len, SA_COMPX);
            // Reaching here means the read failure was swallowed -> bug.
            _exit(0);
        }
        int status = 0;
        CHECK(waitpid(pid, &status, 0) == pid);
        CHECK(WIFEXITED(status));
        CHECK(WEXITSTATUS(status) == EXIT_FAILURE);
    }

    std::fprintf(stderr, "detect_sa_compx_test: all cases passed\n");
    return 0;
}
