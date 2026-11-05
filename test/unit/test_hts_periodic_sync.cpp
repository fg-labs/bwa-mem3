// Fast boundary and error checks; real BAM byte identity lives in integration.
#include "htslib/sam.h"
#include "hfile_internal.h"

#include "doctest/doctest.h"
#include "scoped_env.h"
#include "hts_periodic_sync.h"

#include <cerrno>
#include <initializer_list>
#include <string>

namespace {

// Exercise the real hflush frontend without filesystem I/O. Only its backend
// callback counts flushes, independently of the periodic helper's bookkeeping.
struct FlushStream : hFILE {
    int flushes = 0;
    bool fail = false;
    char storage = 0;
    static int flush(hFILE *fp) {
        auto *self = static_cast<FlushStream *>(fp);
        ++self->flushes;
        if (self->fail) { errno = EIO; return -1; }
        return 0;
    }
    FlushStream() : hFILE{} {
        static const hFILE_backend ops = {nullptr, nullptr, nullptr, flush, nullptr};
        buffer = begin = end = limit = &storage;
        backend = &ops;
    }
};

struct Output {
    FlushStream stream;
    BGZF bgzf{};
    htsFile fp{};
    Output() {
        bgzf.fp = &stream;
        fp.is_bgzf = 1;
        fp.fp.bgzf = &bgzf;
    }
};

}  // namespace

TEST_CASE("hts_periodic_sync: flushes at interval boundaries"
          * doctest::test_suite("unit/mem_io")) {
    bwa_tests::ScopedEnv env("BWA_MEM3_BAM_PERIODIC_SYNC", "");
    Output out;
    hts_periodic_sync_t s;
    hts_periodic_sync_init(&s);
    for (int i = 1; i <= 3; ++i) {
        out.bgzf.block_address = i * HTS_PERIODIC_SYNC_BYTES - 1;
        CHECK(hts_periodic_sync(&s, &out.fp) == 0);
        CHECK(out.stream.flushes == i - 1);
        ++out.bgzf.block_address;
        CHECK(hts_periodic_sync(&s, &out.fp) == 0);
        CHECK(out.stream.flushes == i);
        CHECK(s.n_syncs == i);
        CHECK(hts_periodic_sync(&s, &out.fp) == 0);
        CHECK(out.stream.flushes == i);
    }
}

TEST_CASE("hts_periodic_sync: reports backend flush failures"
          * doctest::test_suite("unit/mem_io")) {
    bwa_tests::ScopedEnv env("BWA_MEM3_BAM_PERIODIC_SYNC", "");
    Output out;
    out.stream.fail = true;
    out.bgzf.block_address = HTS_PERIODIC_SYNC_BYTES;
    hts_periodic_sync_t s;
    hts_periodic_sync_init(&s);
    CHECK(hts_periodic_sync(&s, &out.fp) == -1);
    CHECK(out.stream.flushes == 1);
    CHECK(out.stream.has_errno == EIO);
}

TEST_CASE("hts_periodic_sync: stands aside for ineligible outputs"
          * doctest::test_suite("unit/mem_io")) {
    bwa_tests::ScopedEnv env("BWA_MEM3_BAM_PERIODIC_SYNC", "");
    Output out;
    out.bgzf.block_address = HTS_PERIODIC_SYNC_BYTES;
    hts_periodic_sync_t s;
    hts_periodic_sync_init(&s);
    SUBCASE("multithreaded writer") {
        // The helper checks presence only; no pool is dereferenced.
        out.bgzf.mt = reinterpret_cast<decltype(out.bgzf.mt)>(&out.stream);
    }
    SUBCASE("SAM text") { out.fp.is_bgzf = 0; }
    SUBCASE("missing BGZF stream") { out.fp.fp.bgzf = nullptr; }
    CHECK(hts_periodic_sync(&s, &out.fp) == 0);
    CHECK(hts_periodic_sync(&s, nullptr) == 0);
    CHECK(hts_periodic_sync(nullptr, &out.fp) == 0);
    CHECK(out.stream.flushes == 0);
    CHECK(s.n_syncs == 0);
}

TEST_CASE("hts_periodic_sync: environment settings are isolated and restored"
          * doctest::test_suite("unit/mem_io")) {
    bwa_tests::ScopedEnv outer("BWA_MEM3_BAM_PERIODIC_SYNC", "prior-value");
    for (const char *value : {"0", "1", "", "no"}) {
        CAPTURE(value);
        {
            bwa_tests::ScopedEnv env("BWA_MEM3_BAM_PERIODIC_SYNC", value);
            Output out;
            out.bgzf.block_address = HTS_PERIODIC_SYNC_BYTES;
            hts_periodic_sync_t s;
            hts_periodic_sync_init(&s);
            REQUIRE(hts_periodic_sync(&s, &out.fp) == 0);
            CHECK(out.stream.flushes == (value[0] == '0' ? 0 : 1));
        }
        CHECK(std::string(getenv("BWA_MEM3_BAM_PERIODIC_SYNC")) == "prior-value");
    }
}
