// Real-file cadence and byte identity across multiple sync boundaries.
// Large writes belong here rather than in the fast unit-test shards.

// htslib headers first: they share the KSTRING_H include guard with
// bwa-mem3's kstring.h.
#include "htslib/sam.h"

#include "doctest/doctest.h"
#include "temp_dir.h"
#include "scoped_env.h"

#include "hts_periodic_sync.h"

#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <string>

namespace {

// One unmapped record carrying a long read, so a few thousand records span
// several sync intervals.
struct Record {
    bam1_t *b;
    explicit Record(size_t len) : b(bam_init1()) {
        std::string seq(len, 'A'), qual(len, 'I');
        const char *qname = "read";
        REQUIRE(bam_set1(b, 4, qname, 4, -1, -1, 0, 0, nullptr, -1, -1, 0, len, seq.c_str(),
                         qual.c_str(), 0) >= 0);
    }
    ~Record() { bam_destroy1(b); }
};

struct Result {
    int64_t n_syncs;
    int64_t bytes;  // block_address before close: compressed bytes written
};

// Write `n` deterministically varied records, calling hts_periodic_sync after
// each write when `sync`. Reusing the seed keeps the two output arms identical.
Result write_file(const std::string &path, const char *mode, const Record &rec, int n,
                  bool sync) {
    htsFile *fp = hts_open(path.c_str(), mode);
    REQUIRE(fp != nullptr);
    sam_hdr_t *hdr = sam_hdr_parse(0, "");
    REQUIRE(hdr != nullptr);
    REQUIRE(sam_hdr_write(fp, hdr) == 0);
    hts_periodic_sync_t s;
    hts_periodic_sync_init(&s);
    uint32_t rng = 123456789;
    for (int i = 0; i < n; ++i) {
        // Vary every record: repeating one random read still compresses well.
        auto next = [&]() {
            rng ^= rng << 13; rng ^= rng >> 17; rng ^= rng << 5;
            return rng;
        };
        uint8_t *seq = bam_get_seq(rec.b), *qual = bam_get_qual(rec.b);
        for (int j = 0; j < rec.b->core.l_qseq; ++j) {
            const uint32_t bits = next();
            bam_set_seqi(seq, j, 1 << (bits & 3));
            qual[j] = (bits >> 8) & 63;
        }
        REQUIRE(sam_write1(fp, hdr, rec.b) >= 0);
        if (sync) REQUIRE(hts_periodic_sync(&s, fp) == 0);
    }
    const int64_t bytes = fp->is_bgzf ? fp->fp.bgzf->block_address : 0;
    REQUIRE(hts_close(fp) == 0);
    sam_hdr_destroy(hdr);
    return {s.n_syncs, bytes};
}

bool same_bytes(const std::string &a, const std::string &b) {
    std::ifstream left(a, std::ios::binary), right(b, std::ios::binary);
    REQUIRE(left.is_open());
    REQUIRE(right.is_open());
    char x[65536], y[65536];
    do {
        left.read(x, sizeof(x));
        right.read(y, sizeof(y));
        if (left.gcount() != right.gcount() ||
            memcmp(x, y, static_cast<size_t>(left.gcount())) != 0) return false;
    } while (left && right);
    REQUIRE_FALSE(left.bad());
    REQUIRE_FALSE(right.bad());
    return left.eof() && right.eof();
}

int count_records(const std::string &path) {
    htsFile *fp = hts_open(path.c_str(), "r");
    REQUIRE(fp != nullptr);
    sam_hdr_t *hdr = sam_hdr_read(fp);
    REQUIRE(hdr != nullptr);
    bam1_t *b = bam_init1();
    int n = 0, r;
    while ((r = sam_read1(fp, hdr, b)) >= 0) ++n;
    CHECK(r == -1);  // clean EOF, not a truncated stream
    bam_destroy1(b);
    sam_hdr_destroy(hdr);
    hts_close(fp);
    return n;
}

const int kRecords = 12000;  // enough high-entropy records for three compressed sync intervals

}  // namespace

TEST_CASE("hts_periodic_sync: one sync per interval on the single-threaded writer"
          * doctest::test_suite("integration/mem_io")) {
    bwa_tests::ScopedEnv env("BWA_MEM3_BAM_PERIODIC_SYNC", "");
    bwa_tests::TempDir dir;
    const Record rec(10000);
    const Result r = write_file(dir.file("out.bam"), "wb0", rec, kRecords, true);
    REQUIRE(r.bytes > 3 * HTS_PERIODIC_SYNC_BYTES);
    // A sync fires on the first write at or past each interval boundary, so
    // the count is the number of whole intervals written.
    CHECK(r.n_syncs == r.bytes / HTS_PERIODIC_SYNC_BYTES);
    CHECK(count_records(dir.file("out.bam")) == kRecords);
}

namespace {

void check_identity(const char *mode) {
    bwa_tests::ScopedEnv env("BWA_MEM3_BAM_PERIODIC_SYNC", "");
    bwa_tests::TempDir dir;
    const Record rec(10000);
    const Result with = write_file(dir.file("with.bam"), mode, rec, kRecords, true);
    write_file(dir.file("without.bam"), mode, rec, kRecords, false);
    CHECK(with.n_syncs >= 3);
    CHECK(same_bytes(dir.file("with.bam"), dir.file("without.bam")));
    CHECK(count_records(dir.file("with.bam")) == kRecords);
}

}  // namespace

TEST_CASE("hts_periodic_sync: uncompressed output bytes are unchanged"
          * doctest::test_suite("integration/mem_io")) {
    check_identity("wb0");
}

TEST_CASE("hts_periodic_sync: compressed output bytes are unchanged"
          * doctest::test_suite("integration/mem_io")) {
    check_identity("wb1");
}
