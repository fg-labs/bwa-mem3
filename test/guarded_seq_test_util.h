// Shared test helper: a byte sequence copied flush against a PROT_NONE guard page,
// so a read outside [data, data + len) faults instead of passing on stray bytes.
// Used by the SoA-pack tests (test/x86_soa_pack_test.cpp, test/unit/test_neon_soa_pack.cpp)
// and the ungapped fast-path tail test (test/unit/test_ungapped_fastpath.cpp).
// Header-only; include it after doctest/doctest.h (it uses REQUIRE).
#ifndef GUARDED_SEQ_TEST_UTIL_H
#define GUARDED_SEQ_TEST_UTIL_H

#include <cstddef>
#include <cstdint>
#include <cstring>
#include <vector>

#include <sys/mman.h>
#include <unistd.h>

#include "doctest/doctest.h"

namespace bwa_tests {

// The copy sits against the guard page after its last byte (guardAfter) or before its
// first byte. data() is nullptr for an empty sequence.
class GuardedSeq {
  public:
    GuardedSeq(const std::vector<uint8_t> &bytes, bool guardAfter) : len_(bytes.size()) {
        page_ = (size_t)sysconf(_SC_PAGESIZE);
        const size_t body = ((len_ + page_ - 1) / page_) * page_ + page_;   // at least one page
        map_size_ = body + 2 * page_;
        void *m = mmap(nullptr, map_size_, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANON, -1, 0);
        REQUIRE(m != MAP_FAILED);
        map_ = (uint8_t *)m;
        REQUIRE(mprotect(map_, page_, PROT_NONE) == 0);
        REQUIRE(mprotect(map_ + page_ + body, page_, PROT_NONE) == 0);
        data_ = guardAfter ? map_ + page_ + body - len_ : map_ + page_;
        if (len_) memcpy(data_, bytes.data(), len_);
    }
    ~GuardedSeq() { munmap(map_, map_size_); }
    GuardedSeq(const GuardedSeq &) = delete;
    GuardedSeq &operator=(const GuardedSeq &) = delete;
    const uint8_t *data() const { return len_ ? data_ : nullptr; }

  private:
    size_t len_, page_ = 0, map_size_ = 0;
    uint8_t *map_ = nullptr, *data_ = nullptr;
};

} // namespace bwa_tests

#endif // GUARDED_SEQ_TEST_UTIL_H
