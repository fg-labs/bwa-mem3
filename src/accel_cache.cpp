#include "accel_cache.h"

#include <cstring>
#include <cstdio>
#include <cerrno>
#include <cstdlib>
#include <fcntl.h>
#include <unistd.h>
#include <sys/mman.h>
#include <sys/stat.h>

extern "C" {
#include "../ext/sha256/sha-256.h"
}

// -- K-mer encoding ----------------------------------------------------

uint64_t accel_pack_forward(const uint8_t *enc, uint32_t k) {
    uint64_t v = 0;
    for (uint32_t i = 0; i < k; ++i) {
        uint8_t b = enc[i];
        if (b > 3) return UINT64_MAX;
        v |= (uint64_t)b << (2 * i);
    }
    return v;
}

uint64_t accel_revcomp(uint64_t kmer_packed, uint32_t k) {
    uint64_t r = 0;
    for (uint32_t i = 0; i < k; ++i) {
        uint8_t b = (kmer_packed >> (2 * i)) & 0x3u;
        uint8_t c = b ^ 0x3u;
        r |= (uint64_t)c << (2 * (k - 1 - i));
    }
    return r;
}

uint64_t accel_canonicalize(uint64_t fwd, uint32_t k, bool *out_is_rc) {
    uint64_t rc = accel_revcomp(fwd, k);
    bool use_rc = rc < fwd;
    if (out_is_rc) *out_is_rc = use_rc;
    return use_rc ? rc : fwd;
}

// -- SHA-256 -----------------------------------------------------------

bool accel_sha256_file(const char *path, uint8_t out[32]) {
    FILE *f = fopen(path, "rb");
    if (!f) return false;
    struct Sha_256 sha;
    sha_256_init(&sha, out);
    uint8_t buf[1 << 16];
    size_t n;
    while ((n = fread(buf, 1, sizeof(buf), f)) > 0) {
        sha_256_write(&sha, buf, n);
    }
    sha_256_close(&sha);
    fclose(f);
    return true;
}

// -- Hash -------------------------------------------------------------

static inline uint64_t accel_hash64(uint64_t x) {
    x += 0x9E3779B97F4A7C15ULL;
    x = (x ^ (x >> 30)) * 0xBF58476D1CE4E5B9ULL;
    x = (x ^ (x >> 27)) * 0x94D049BB133111EBULL;
    return x ^ (x >> 31);
}

static AccelLookupResult probe_table(const AccelLevelTable &t, uint64_t k) {
    if (t.n_entries == 0 || t.probe == nullptr) return {false, 0, 0, 0};
    uint64_t slot = accel_hash64(k) & t.table_mask;
    // Linear probing until empty slot or match.
    while (true) {
        uint64_t idx1 = t.probe[slot];
        if (idx1 == 0) return {false, 0, 0, 0};
        const AccelEntry &e = t.entries[idx1 - 1];
        if (e.kmer_canonical == k) return {true, e.sa_k, e.sa_l, e.sa_s};
        slot = (slot + 1) & t.table_mask;
    }
}

// -- Cache lifecycle --------------------------------------------------

AccelCache::AccelCache()
    : k_min_(0), k_max_(0), mmap_base_(nullptr), mmap_len_(0) {
    memset(levels_, 0, sizeof(levels_));
    memset(&terminal_, 0, sizeof(terminal_));
}

AccelCache::~AccelCache() {
    if (mmap_base_) {
        ::munmap(mmap_base_, mmap_len_);
        mmap_base_ = nullptr;
        mmap_len_ = 0;
    }
}

bool AccelCache::load(const char *path, const uint8_t ref_sha256[32]) {
    int fd = ::open(path, O_RDONLY);
    if (fd < 0) {
        fprintf(stderr, "[accel] cannot open %s: %s\n", path, strerror(errno));
        return false;
    }
    struct stat st;
    if (::fstat(fd, &st) != 0) {
        ::close(fd);
        return false;
    }
    size_t len = (size_t)st.st_size;
    void *base = ::mmap(nullptr, len, PROT_READ, MAP_PRIVATE, fd, 0);
    ::close(fd);
    if (base == MAP_FAILED) {
        fprintf(stderr, "[accel] mmap failed for %s: %s\n", path, strerror(errno));
        return false;
    }

    const uint8_t *p = (const uint8_t *)base;
    const uint8_t *end = p + len;

#define ACCEL_READ(dst, sz) do { \
        if ((p) + (sz) > (end)) { ::munmap(base, len); return false; } \
        memcpy((dst), (p), (sz)); (p) += (sz); \
    } while (0)

    // File header.
    uint64_t magic;             ACCEL_READ(&magic, 8);
    uint32_t version;           ACCEL_READ(&version, 4);
    uint8_t  sha[32];           ACCEL_READ(sha, 32);
    uint32_t k_min;             ACCEL_READ(&k_min, 4);
    uint32_t k_max;             ACCEL_READ(&k_max, 4);
    uint32_t tier_count;        ACCEL_READ(&tier_count, 4);
    uint32_t collapse_thresh;   ACCEL_READ(&collapse_thresh, 4);
    uint8_t  reserved[16];      ACCEL_READ(reserved, 16);
    (void)collapse_thresh; (void)reserved;

    if (magic != ACCEL_MAGIC) {
        fprintf(stderr, "[accel] bad magic in %s\n", path);
        ::munmap(base, len); return false;
    }
    if (version != ACCEL_VERSION) {
        fprintf(stderr, "[accel] unsupported version %u in %s (expected %u)\n",
                version, path, ACCEL_VERSION);
        ::munmap(base, len); return false;
    }
    if (memcmp(sha, ref_sha256, 32) != 0) {
        fprintf(stderr, "[accel] ref checksum mismatch in %s\n", path);
        ::munmap(base, len); return false;
    }
    if (k_min < ACCEL_K_MIN || k_max > ACCEL_MAX_L || k_min > k_max) {
        fprintf(stderr, "[accel] bad k range [%u, %u] in %s\n", k_min, k_max, path);
        ::munmap(base, len); return false;
    }
    if (tier_count != 1) {
        fprintf(stderr, "[accel] this bwa-mem2 supports tier_count=1 only (got %u)\n",
                tier_count);
        ::munmap(base, len); return false;
    }

    // Tier header (we ignore id for Phase 1).
    uint8_t tier_id;            ACCEL_READ(&tier_id, 1);
    uint8_t tier_pad[7];        ACCEL_READ(tier_pad, 7);
    (void)tier_id; (void)tier_pad;

    k_min_ = k_min;
    k_max_ = k_max;

    // Per-level tables.
    for (uint32_t L = k_min; L <= k_max; ++L) {
        uint64_t n_entries, table_mask;
        ACCEL_READ(&n_entries, 8);
        ACCEL_READ(&table_mask, 8);
        size_t probe_bytes = (table_mask + 1) * sizeof(uint64_t);
        size_t entries_bytes = n_entries * sizeof(AccelEntry);
        if ((size_t)(end - p) < probe_bytes + entries_bytes) {
            fprintf(stderr, "[accel] truncated at level L=%u\n", L);
            ::munmap(base, len); return false;
        }
        levels_[L].n_entries  = n_entries;
        levels_[L].table_mask = table_mask;
        levels_[L].probe      = (const uint64_t *)p;   p += probe_bytes;
        levels_[L].entries    = (const AccelEntry *)p; p += entries_bytes;
    }

    // Terminal table.
    {
        uint64_t n_entries, table_mask;
        ACCEL_READ(&n_entries, 8);
        ACCEL_READ(&table_mask, 8);
        size_t probe_bytes = (table_mask + 1) * sizeof(uint64_t);
        size_t entries_bytes = n_entries * sizeof(AccelEntry);
        if ((size_t)(end - p) < probe_bytes + entries_bytes) {
            fprintf(stderr, "[accel] truncated terminal table\n");
            ::munmap(base, len); return false;
        }
        terminal_.n_entries  = n_entries;
        terminal_.table_mask = table_mask;
        terminal_.probe      = (const uint64_t *)p;   p += probe_bytes;
        terminal_.entries    = (const AccelEntry *)p; p += entries_bytes;
    }

    mmap_base_ = base;
    mmap_len_  = len;
    fprintf(stderr, "[accel] loaded %s: k=[%u..%u], terminal=%llu entries\n",
            path, k_min_, k_max_, (unsigned long long)terminal_.n_entries);
    return true;

#undef ACCEL_READ
}

AccelLookupResult AccelCache::probe(uint32_t L, uint64_t kmer_canonical) const {
    if (!enabled() || L < k_min_ || L > k_max_) return {false, 0, 0, 0};
    return probe_table(levels_[L], kmer_canonical);
}

AccelLookupResult AccelCache::probe_terminal(uint64_t kmer_canonical) const {
    if (!enabled()) return {false, 0, 0, 0};
    return probe_table(terminal_, kmer_canonical);
}
