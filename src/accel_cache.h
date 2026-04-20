#ifndef ACCEL_CACHE_H
#define ACCEL_CACHE_H

#include <cstdint>
#include <cstddef>

// Target-region fast-path accelerator cache.
//
// Multi-level, read-only, canonical-k-mer → SA-interval lookup.
// Mmap'd at load, probed from FMI_search during SMEM extension to
// replace LF-map chains with hash probes while producing identical
// SA intervals (R1 correctness).
//
// See docs/superpowers/specs/2026-04-20-target-region-fast-path-design.md.

static constexpr uint64_t ACCEL_MAGIC   = 0x4C4543434132574DULL; // "MW2ACCEL" LE
static constexpr uint32_t ACCEL_VERSION = 1;
static constexpr uint8_t  ACCEL_TIER_A  = 'A';
static constexpr uint8_t  ACCEL_TIER_B  = 'B';
static constexpr uint8_t  ACCEL_TIER_C  = 'C';
static constexpr uint32_t ACCEL_K_MIN   = 1;
static constexpr uint32_t ACCEL_K_MAX   = 18;
static constexpr uint32_t ACCEL_MAX_L   = 32; // sized upper bound for levels_[] array

// One entry of a per-level table or the terminal table.
struct AccelEntry {
    uint64_t kmer_canonical;   // 2 bits/base packed, low 2*L bits
    int64_t  sa_k;
    int64_t  sa_l;
    int64_t  sa_s;
} __attribute__((packed));
static_assert(sizeof(AccelEntry) == 32, "AccelEntry must be 32 bytes");

struct AccelLookupResult {
    bool    hit;
    int64_t sa_k;
    int64_t sa_l;
    int64_t sa_s;
};

struct AccelLevelTable {
    uint64_t         n_entries;
    uint64_t         table_mask;   // slots = mask + 1 (power of two)
    const uint64_t  *probe;        // mmap'd: mask+1 slots; 0=empty, else 1-based index
    const AccelEntry *entries;     // mmap'd: n_entries records
};

class AccelCache {
public:
    AccelCache();
    ~AccelCache();

    // Load a cache file. Verifies magic/version and that the embedded
    // ref SHA-256 matches the caller-supplied checksum of the current
    // reference's .bwt.2bit.64. Returns false on any mismatch or IO
    // error; the cache is left disabled.
    bool load(const char *path, const uint8_t ref_sha256[32]);

    bool     enabled()  const { return mmap_base_ != nullptr; }
    uint32_t k_min()    const { return k_min_; }
    uint32_t k_max()    const { return k_max_; }

    // Probe the per-length table at level L. Caller must hold L in
    // [k_min(), k_max()]; returns {hit=false} for out-of-range.
    AccelLookupResult probe(uint32_t L, uint64_t kmer_canonical) const;

    // Probe the terminal table at length k_max.
    AccelLookupResult probe_terminal(uint64_t kmer_canonical) const;

private:
    uint32_t k_min_;
    uint32_t k_max_;
    AccelLevelTable levels_[ACCEL_MAX_L + 1]; // indexed 1..k_max
    AccelLevelTable terminal_;
    void   *mmap_base_;
    size_t  mmap_len_;
};

// -- Free helpers ------------------------------------------------------

// Pack the first k 2-bit-encoded bases (A=0..T=3) into a uint64, low
// bits first. Returns UINT64_MAX if any base is > 3 (N).
uint64_t accel_pack_forward(const uint8_t *enc_bases, uint32_t k);

// Reverse-complement a packed k-mer.
uint64_t accel_revcomp(uint64_t kmer_packed, uint32_t k);

// Canonical form = lex-min of (fwd, revcomp(fwd)). Optionally reports
// whether the canonical form was the RC (caller must swap k/l after
// probe if so).
uint64_t accel_canonicalize(uint64_t fwd, uint32_t k, bool *out_is_rc = nullptr);

// SHA-256 of a file's contents.
bool accel_sha256_file(const char *path, uint8_t out[32]);

#endif // ACCEL_CACHE_H
