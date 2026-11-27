// Unit tests for the packed SA-sample high parts (FMI_search::sa_hi_*).
//
// On disk each SA sample's high part (bits 32.. of its position) is one byte; in
// memory it is packed at (1 << sa_hi_log2) bits per sample. These tests pin the
// width selection against the largest position a BWT of ref_seq_len rows can
// hold, the round trip through sa_hi_pack for every width (including a chunked
// pack that starts mid-array), and the out-of-range rejection. The production
// accessors (FMI_search::sa_hi_at / sa_hi_addr, with the geometry
// set_sa_hi_geometry derives from ref_seq_len) are checked at every width, and
// the threaded pack, which splits only from 2 x 4096 words up, against the serial
// one.

#define private public
#include "FMI_search.h"
#undef private

#include "doctest/doctest.h"
#include "temp_dir.h"

#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <string>
#include <fcntl.h>
#include <unistd.h>
#include <sys/wait.h>
#include <vector>

namespace {

// Pseudo-random high parts that fit `bits` bits; at 8 bits the whole byte, so
// values 128..255 (negative as int8_t) are packed through their uint8_t value.
std::vector<int8_t> high_parts(int64_t n, int bits, uint64_t seed) {
    std::vector<int8_t> src((size_t)n);
    uint64_t x = seed;
    for (int64_t i = 0; i < n; ++i) {
        x ^= x << 13; x ^= x >> 7; x ^= x << 17;
        src[(size_t)i] = (int8_t)(uint8_t)(x & ((bits >= 8) ? 0xffu : ((1u << bits) - 1)));
    }
    return src;
}

// A reference length whose positions need exactly 1 << log2 high-part bits.
int64_t ref_len_for_log2(int log2) {
    const int64_t G = (int64_t)1 << 32;
    static const int64_t len[4] = {2 * G, 4 * G, 16 * G, 0x7fffffffffLL};
    return len[log2];
}

int64_t unpack(const std::vector<uint64_t> &w, int log2, int64_t i) {
    const uint64_t b = (uint64_t)i << log2;
    const uint64_t mask = (1ull << (1u << log2)) - 1;
    return (int64_t)((w[b >> 6] >> (b & 63)) & mask);
}

// Decode every sample independently of the production packer, with one assertion.
void check_samples(const std::vector<uint64_t> &words, int log2, const std::vector<int8_t> &src) {
    int64_t mismatches = 0;
    for (size_t i = 0; i < src.size(); ++i)
        mismatches += unpack(words, log2, (int64_t)i) != (int64_t)(uint8_t)src[i];
    CHECK(mismatches == 0);
}

void check_zero_padding(const std::vector<uint64_t> &words, int64_t n, int log2) {
    const uint64_t used_bits = ((uint64_t)n << log2) & 63;
    REQUIRE(used_bits != 0);
    const uint64_t unused_mask = ~uint64_t{0} << used_bits;
    CHECK((words.back() & unused_mask) == 0);
}

}  // namespace

TEST_SUITE_BEGIN("unit/fmi_sa_hi_pack");

TEST_CASE("sa_hi_log2_for picks the narrowest width holding every position") {
    const int64_t G = (int64_t)1 << 32;
    CHECK(FMI_search::sa_hi_log2_for(1) == 0);
    CHECK(FMI_search::sa_hi_log2_for(G) == 0);           // max position G-1: high part 0
    CHECK(FMI_search::sa_hi_log2_for(2 * G) == 0);       // high part <= 1: 1 bit
    CHECK(FMI_search::sa_hi_log2_for(6400000000LL) == 0);  // GRCh38-sized (2L+1)
    CHECK(FMI_search::sa_hi_log2_for(2 * G + 1) == 1);   // high part 2: 2 bits
    CHECK(FMI_search::sa_hi_log2_for(4 * G) == 1);
    CHECK(FMI_search::sa_hi_log2_for(4 * G + 1) == 2);
    CHECK(FMI_search::sa_hi_log2_for(16 * G) == 2);
    CHECK(FMI_search::sa_hi_log2_for(16 * G + 1) == 3);
    CHECK(FMI_search::sa_hi_log2_for(0x7fffffffffLL) == 3);
}

TEST_CASE("sa_hi_pack round-trips every width, whole and chunked") {
    for (int log2 = 0; log2 <= 3; ++log2) {
        CAPTURE(log2);
        const int bits = 1 << log2;
        const int64_t n = 1000 + 37;  // not a whole number of words
        const std::vector<int8_t> src = high_parts(n, bits, 0x9e3779b97f4a7c15ull);
        const int64_t words = FMI_search::sa_hi_words(n, log2);
        // The chunked destination starts all ones: the pack assigns whole words, so it
        // needs no zeroing pass, including the unused fields of the partial last word.
        std::vector<uint64_t> whole((size_t)words, 0), chunked((size_t)words, ~0ull);
        REQUIRE(FMI_search::sa_hi_pack(src.data(), 0, n, log2, whole.data(), 4));
        const int64_t step = 2 * (64 >> log2);  // chunk starts stay on word boundaries
        for (int64_t first = 0; first < n; first += step) {
            const int64_t m = (n - first < step) ? n - first : step;
            REQUIRE(FMI_search::sa_hi_pack(src.data() + first, first, m, log2, chunked.data(), 1));
        }
        CHECK(whole == chunked);
        check_zero_padding(whole, n, log2);
        check_zero_padding(chunked, n, log2);
        for (int64_t i = 0; i < n; ++i) CHECK(unpack(whole, log2, i) == (int64_t)(uint8_t)src[(size_t)i]);
    }
}

namespace {
// Restore the actual prior override even when a fatal assertion unwinds the test.
class ScopedSaHiReadChunk {
public:
    explicit ScopedSaHiReadChunk(int64_t n) : saved_(FMI_search::sa_hi_read_chunk()) {
        FMI_search::set_sa_hi_read_chunk(n);
    }
    ~ScopedSaHiReadChunk() { FMI_search::set_sa_hi_read_chunk(saved_); }
    ScopedSaHiReadChunk(const ScopedSaHiReadChunk &) = delete;
    ScopedSaHiReadChunk &operator=(const ScopedSaHiReadChunk &) = delete;
private:
    const int64_t saved_;
};

// Write `pad` junk bytes and then `src` to file `name` in `dir` (removed with the
// TempDir); return its open fd and set `path`.
int bytes_file(const bwa_tests::TempDir &dir, const char *name, const std::vector<int8_t> &src,
               size_t pad, std::string &path) {
    path = dir.file(name);
    const int fd = open(path.c_str(), O_RDWR | O_CREAT | O_TRUNC, 0600);
    REQUIRE(fd >= 0);
    std::vector<char> junk(pad, (char)0x7f);
    REQUIRE(write(fd, junk.data(), pad) == (ssize_t)pad);
    REQUIRE(write(fd, src.data(), src.size()) == (ssize_t)src.size());
    return fd;
}
}  // namespace

TEST_CASE("sa_hi_pread_pack reads in pieces and packs like one whole-array pack") {
    // A piece of 3 x 64 samples and a sample count that leaves a partial last piece and a
    // partial last word, so every piece after the first starts at a nonzero word; the bytes
    // sit after a 13-byte header, as they do after CP_OCC in the index file.
    const ScopedSaHiReadChunk chunk(3 * 64);
    for (int log2 = 0; log2 <= 3; ++log2) {
        CAPTURE(log2);
        const int bits = 1 << log2;
        const int64_t n = 7 * 3 * 64 + 100;
        std::vector<int8_t> src = high_parts(n, bits, 0xa4093822299f31d0ull + (uint64_t)log2);
        const int64_t words = FMI_search::sa_hi_words(n, log2);
        std::vector<uint64_t> whole((size_t)words, 0), got((size_t)words, ~0ull);
        REQUIRE(FMI_search::sa_hi_pack(src.data(), 0, n, log2, whole.data(), 1));
        bwa_tests::TempDir dir("sa_hi_pread");
        std::string path;
        const int fd = bytes_file(dir, "ok", src, 13, path);
        CHECK(FMI_search::sa_hi_pread_pack(fd, path.c_str(), 13, n, log2, got.data(), 4) == 0);
        CHECK(whole == got);
        check_zero_padding(whole, n, log2);
        check_zero_padding(got, n, log2);
        for (int64_t i = 0; i < n; ++i) CHECK(unpack(got, log2, i) == (int64_t)(uint8_t)src[(size_t)i]);
        close(fd);
        if (bits < 8) {   // an over-wide byte in a later piece
            src[(size_t)(4 * 3 * 64 + 5)] = (int8_t)0x80;
            const int bfd = bytes_file(dir, "bad", src, 13, path);
            CHECK(FMI_search::sa_hi_pread_pack(bfd, path.c_str(), 13, n, log2, got.data(), 1)
                  == FMI_search::SA_HI_READ_RANGE);
            close(bfd);
        }
    }
}

TEST_CASE("sa_hi_pread_pack splits over several workers at whole words") {
    // The worker count follows the index load's (one per 8 MiB of bytes), so 3 workers need
    // 24 MiB; the uneven split and the pieces check every worker's range ends.
    const ScopedSaHiReadChunk chunk(4096 * 64);
    const int64_t n = 3 * ((int64_t)8 << 20) + 777;
    std::vector<int8_t> src = high_parts(n, 1, 0x452821e638d01377ull);
    const int64_t words = FMI_search::sa_hi_words(n, 0);
    std::vector<uint64_t> whole((size_t)words, 0), got((size_t)words, ~0ull);
    REQUIRE(FMI_search::sa_hi_pack(src.data(), 0, n, 0, whole.data(), 1));
    bwa_tests::TempDir dir("sa_hi_pread_workers");
    std::string path;
    const int fd = bytes_file(dir, "ok", src, 64, path);
    CHECK(fmi_pread_worker_count((size_t)n, 7) == 3);
    CHECK(FMI_search::sa_hi_pread_pack(fd, path.c_str(), 64, n, 0, got.data(), 7) == 0);
    CHECK(whole == got);
    check_samples(got, 0, src);
    check_zero_padding(whole, n, 0);
    check_zero_padding(got, n, 0);
    close(fd);
    src[(size_t)(n / 2)] = (int8_t)0xff;   // the middle worker's range
    const int bfd = bytes_file(dir, "bad", src, 64, path);
    CHECK(FMI_search::sa_hi_pread_pack(bfd, path.c_str(), 64, n, 0, got.data(), 7)
          == FMI_search::SA_HI_READ_RANGE);
    close(bfd);
}

TEST_CASE("sa_hi_pread_pack returns read failures after joining its workers") {
    for (const int threads : {1, 3}) {
        for (const bool bad_fd : {false, true}) {
            CAPTURE(threads);
            CAPTURE(bad_fd);
            bwa_tests::TempDir dir("sa_hi_read_error");
            std::string path;
            const int fd = bytes_file(dir, "short", std::vector<int8_t>(64, 0), 0, path);
            const int64_t n = threads == 1 ? 128 : 3 * (int64_t)FMI_PREAD_MIN_CHUNK;
            std::vector<uint64_t> got((size_t)FMI_search::sa_hi_words(n, 0));
            fflush(NULL);
            const pid_t pid = fork();
            REQUIRE(pid >= 0);
            if (pid == 0) {
                // An absent path keeps all workers on the deliberately invalid fd.
                const int rc = FMI_search::sa_hi_pread_pack(bad_fd ? -1 : fd,
                    bad_fd ? dir.file("absent").c_str() : path.c_str(), 0, n, 0,
                    got.data(), threads);
                _exit(rc == FMI_search::SA_HI_READ_IO ? 0 : 2);
            }
            int status = 0;
            REQUIRE(waitpid(pid, &status, 0) == pid);
            CHECK(WIFEXITED(status));
            CHECK(WEXITSTATUS(status) == 0);
            close(fd);
        }
    }
}

TEST_CASE("sa_hi_pread_pack reports EOF confined to a spawned worker") {
    bwa_tests::TempDir dir("sa_hi_spawned_eof");
    const std::string path = dir.file("short");
    const int fd = open(path.c_str(), O_RDWR | O_CREAT | O_EXCL, 0600);
    REQUIRE(fd >= 0);
    const int64_t n = 3 * (int64_t)FMI_PREAD_MIN_CHUNK;
    REQUIRE(fmi_pread_worker_count((size_t)n, 3) == 3);
    // Sparse zero-filled bytes cover the first two readers completely; only
    // the last spawned worker encounters EOF, so its status must be collected.
    REQUIRE(ftruncate(fd, (off_t)(2 * FMI_PREAD_MIN_CHUNK)) == 0);
    std::vector<uint64_t> got((size_t)FMI_search::sa_hi_words(n, 0));
    fflush(NULL);
    const pid_t pid = fork();
    REQUIRE(pid >= 0);
    if (pid == 0) {
        const int rc = FMI_search::sa_hi_pread_pack(fd, path.c_str(), 0, n, 0, got.data(), 3);
        _exit(rc == FMI_search::SA_HI_READ_IO ? 0 : 2);
    }
    int status = 0;
    REQUIRE(waitpid(pid, &status, 0) == pid);
    CHECK(WIFEXITED(status));
    CHECK(WEXITSTATUS(status) == 0);
    close(fd);
}

TEST_CASE("ordinary parallel index reads remain fatal on EOF") {
    for (const int threads : {1, 3}) {
        bwa_tests::TempDir dir("ordinary_pread_error");
        std::string path;
        const int fd = bytes_file(dir, "short", std::vector<int8_t>(64, 0), 0, path);
        FILE *fp = fdopen(fd, "rb");
        REQUIRE(fp != NULL);
        const size_t n = threads == 1 ? 128 : 3 * FMI_PREAD_MIN_CHUNK;
        std::vector<char> got(n);
        fflush(NULL);
        const pid_t pid = fork();
        REQUIRE(pid >= 0);
        if (pid == 0) {
            fmi_pread_from_stream(fp, path.c_str(), got.data(), n, threads);
            _exit(2);  // A failed ordinary load must never return to its caller.
        }
        int status = 0;
        REQUIRE(waitpid(pid, &status, 0) == pid);
        CHECK(WIFEXITED(status));
        CHECK(WEXITSTATUS(status) == EXIT_FAILURE);
        fclose(fp);
    }
}

TEST_CASE("ordinary parallel index reads remain fatal on spawned-worker-only EOF") {
    bwa_tests::TempDir dir("ordinary_spawned_eof");
    const std::string path = dir.file("short");
    const int fd = open(path.c_str(), O_RDWR | O_CREAT | O_EXCL, 0600);
    REQUIRE(fd >= 0);
    const size_t n = 3 * FMI_PREAD_MIN_CHUNK;
    REQUIRE(fmi_pread_worker_count(n, 3) == 3);
    // The stream starts at zero: the first two readers succeed and only the
    // last spawned worker can terminate this ordinary index read.
    REQUIRE(ftruncate(fd, (off_t)(2 * FMI_PREAD_MIN_CHUNK)) == 0);
    FILE *fp = fdopen(fd, "rb");
    REQUIRE(fp != NULL);
    std::vector<char> got(n);
    fflush(NULL);
    const pid_t pid = fork();
    REQUIRE(pid >= 0);
    if (pid == 0) {
        fmi_pread_from_stream(fp, path.c_str(), got.data(), n, 3);
        _exit(2);  // Ignoring a spawned reader's failure must not return success.
    }
    int status = 0;
    REQUIRE(waitpid(pid, &status, 0) == pid);
    CHECK(WIFEXITED(status));
    CHECK(WEXITSTATUS(status) == EXIT_FAILURE);
    fclose(fp);
}

TEST_CASE("sa_hi_pack rejects a high part wider than the chosen width") {
    std::vector<int8_t> src = {0, 1, 2, 1};
    std::vector<uint64_t> w((size_t)FMI_search::sa_hi_words((int64_t)src.size(), 0), 0);
    CHECK_FALSE(FMI_search::sa_hi_pack(src.data(), 0, (int64_t)src.size(), 0, w.data(), 1));
    std::vector<uint64_t> w1((size_t)FMI_search::sa_hi_words((int64_t)src.size(), 1), 0);
    CHECK(FMI_search::sa_hi_pack(src.data(), 0, (int64_t)src.size(), 1, w1.data(), 1));
    // Bytes with the sign bit set (negative as int8_t) are over-wide at 1, 2 and 4 bits; a
    // signed comparison against the mask would accept them.
    for (int log2 = 0; log2 <= 2; ++log2) {
        for (int b : {0x80, 0xff}) {
            CAPTURE(log2);
            CAPTURE(b);
            std::vector<int8_t> one = {0, (int8_t)(uint8_t)b, 0};
            std::vector<uint64_t> wo(1, 0);
            CHECK_FALSE(FMI_search::sa_hi_pack(one.data(), 0, 3, log2, wo.data(), 1));
        }
    }
}

TEST_CASE("sa_hi_at and sa_hi_addr read every width's packed high parts") {
    for (int log2 = 0; log2 <= 3; ++log2) {
        CAPTURE(log2);
        const int bits = 1 << log2, per_word = 64 >> log2;
        // Three whole words and a partial fourth, so the last sample sits mid-word.
        const int64_t n = 3 * per_word + per_word / 2 + 1;
        const std::vector<int8_t> src = high_parts(n, bits, 0x243f6a8885a308d3ull + (uint64_t)log2);
        FMI_search f("unused-prefix");
        f.reference_seq_len = ref_len_for_log2(log2);
        f.set_sa_hi_geometry();
        REQUIRE(f.sa_hi_log2 == log2);
        REQUIRE(f.sa_hi_mask == (1ull << bits) - 1);
        std::vector<uint64_t> w((size_t)FMI_search::sa_hi_words(n, log2), 0);
        REQUIRE(w.size() == 4);
        REQUIRE(FMI_search::sa_hi_pack(src.data(), 0, n, log2, w.data(), 1));
        check_zero_padding(w, n, log2);
        f.sa_hi = w.data();
        for (int64_t i = 0; i < n; ++i) {
            CAPTURE(i);
            CHECK(f.sa_hi_at(i) == (int64_t)(uint8_t)src[(size_t)i]);
            CHECK(f.sa_hi_addr(i) == (const char *)&w[(size_t)(i / per_word)]);
        }
        // Word edges: the last sample of a word and the first of the next.
        for (int64_t k = 1; k <= 3; ++k) {
            CHECK(f.sa_hi_addr(k * per_word - 1) == (const char *)&w[(size_t)(k - 1)]);
            CHECK(f.sa_hi_addr(k * per_word) == (const char *)&w[(size_t)k]);
        }
        f.sa_hi = NULL;   // detach before allocations or REQUIRE can unwind
        // All-ones and all-zero samples next to each other read back unmixed.
        std::vector<int8_t> edge((size_t)(2 * per_word), 0);
        for (int64_t i = 0; i < 2 * per_word; i += 2) edge[(size_t)i] = (int8_t)(uint8_t)f.sa_hi_mask;
        std::vector<uint64_t> we(2, 0);
        REQUIRE(FMI_search::sa_hi_pack(edge.data(), 0, 2 * per_word, log2, we.data(), 1));
        f.sa_hi = we.data();
        for (int64_t i = 0; i < 2 * per_word; ++i)
            CHECK(f.sa_hi_at(i) == (i % 2 == 0 ? (int64_t)f.sa_hi_mask : 0));
        f.sa_hi = NULL;   // the test's buffer, not the object's
    }
}

TEST_CASE("sa_hi_pack: the threaded pack equals the serial one at every width") {
    // sa_hi_pack splits over workers only with at least 2 x 4096 words (at 1
    // bit, 524288 samples) and gives each whole words; an uneven word count over
    // three workers checks the split's ends. The chunked leg starts mid-array,
    // as the index load's buffered reads do.
    for (int log2 = 0; log2 <= 3; ++log2) {
        CAPTURE(log2);
        const int bits = 1 << log2, per_word = 64 >> log2;
        const int64_t nw = 3 * 4096 + 5;
        const int64_t n = nw * per_word - per_word / 2;   // last word partial
        const std::vector<int8_t> src = high_parts(n, bits, 0x13198a2e03707344ull + (uint64_t)log2);
        std::vector<uint64_t> serial((size_t)FMI_search::sa_hi_words(n, log2), 0);
        std::vector<uint64_t> threaded(serial.size(), 0), chunked(serial.size(), 0);
        REQUIRE((int64_t)serial.size() == nw);
        REQUIRE(FMI_search::sa_hi_pack(src.data(), 0, n, log2, serial.data(), 1));
        REQUIRE(FMI_search::sa_hi_pack(src.data(), 0, n, log2, threaded.data(), 7));
        CHECK(std::memcmp(serial.data(), threaded.data(), serial.size() * sizeof(uint64_t)) == 0);
        const int64_t first = 4 * per_word;   // word 4; the remaining 12289 words split 3 ways
        REQUIRE(FMI_search::sa_hi_pack(src.data(), 0, first, log2, chunked.data(), 1));
        REQUIRE(FMI_search::sa_hi_pack(src.data() + first, first, n - first, log2, chunked.data(), 7));
        CHECK(std::memcmp(serial.data(), chunked.data(), serial.size() * sizeof(uint64_t)) == 0);
        check_zero_padding(serial, n, log2);
        check_zero_padding(threaded, n, log2);
        check_zero_padding(chunked, n, log2);
        check_samples(serial, log2, src);
        check_samples(threaded, log2, src);
        check_samples(chunked, log2, src);
        // A bad byte in the last worker's range, or a sign-bit byte in the middle worker's,
        // is still reported.
        if (bits < 8) {
            std::vector<int8_t> bad = src;
            bad[(size_t)(n - 2)] = (int8_t)(1 << bits);
            std::vector<uint64_t> wb(serial.size(), 0);
            CHECK_FALSE(FMI_search::sa_hi_pack(bad.data(), 0, n, log2, wb.data(), 7));
            for (int b : {0x80, 0xff}) {
                std::vector<int8_t> mid = src;
                mid[(size_t)(n / 2)] = (int8_t)(uint8_t)b;
                CHECK_FALSE(FMI_search::sa_hi_pack(mid.data(), 0, n, log2, wb.data(), 7));
            }
        }
    }
}

TEST_SUITE_END();
