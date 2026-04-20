// Unit tests for src/accel_cache.{h,cpp}.
// Built via test/Makefile; run with `./accel_cache_test`.

#include <cassert>
#include <cstdio>
#include <cstdint>
#include <cstring>
#include <cstdlib>
#include <unistd.h>

#include "../src/accel_cache.h"

// -- Helpers ----------------------------------------------------------

static void test_pack_forward_and_revcomp() {
    // ACGT = (A=0, C=1, G=2, T=3). Packed low-to-high:
    // bits [1:0]=A(00), [3:2]=C(01), [5:4]=G(10), [7:6]=T(11)
    // => 0b11100100 = 0xE4
    uint8_t acgt[] = {0, 1, 2, 3};
    uint64_t f = accel_pack_forward(acgt, 4);
    assert(f == 0xE4);

    uint64_t rc = accel_revcomp(f, 4);
    // RC of ACGT is ACGT (palindrome).
    assert(rc == f);

    // Sanity: RC of AAAA is TTTT.
    uint8_t aaaa[] = {0, 0, 0, 0};
    uint64_t fa = accel_pack_forward(aaaa, 4);
    assert(fa == 0);
    uint64_t rca = accel_revcomp(fa, 4);
    // TTTT packed: 0b11111111 = 0xFF
    assert(rca == 0xFF);

    printf("test_pack_forward_and_revcomp: PASS\n");
}

static void test_canonical() {
    // CCCC: all 1s. Packed 0b01010101 = 0x55. RC(CCCC) = GGGG (G=2). Packed 0xAA.
    // Canonical = min(0x55, 0xAA) = 0x55.
    uint8_t cccc[] = {1, 1, 1, 1};
    uint64_t f = accel_pack_forward(cccc, 4);
    assert(f == 0x55);
    bool is_rc = true;
    uint64_t c = accel_canonicalize(f, 4, &is_rc);
    assert(c == 0x55);
    assert(!is_rc);

    uint8_t gggg[] = {2, 2, 2, 2};
    f = accel_pack_forward(gggg, 4);
    assert(f == 0xAA);
    c = accel_canonicalize(f, 4, &is_rc);
    assert(c == 0x55); // canonicalized back to CCCC
    assert(is_rc);

    printf("test_canonical: PASS\n");
}

static void test_pack_rejects_n() {
    uint8_t bases[] = {0, 4, 2, 3}; // N at position 1
    uint64_t f = accel_pack_forward(bases, 4);
    assert(f == UINT64_MAX);
    printf("test_pack_rejects_n: PASS\n");
}

static void test_sha256_file() {
    const char *path = "/tmp/accel_sha256_test.bin";
    FILE *f = fopen(path, "wb");
    assert(f != nullptr);
    const char *content = "hello\n";
    fwrite(content, 1, strlen(content), f);
    fclose(f);

    uint8_t digest[32];
    bool ok = accel_sha256_file(path, digest);
    assert(ok);

    // Known SHA-256 of "hello\n"
    uint8_t expected[32] = {
        0x58,0x91,0xb5,0xb5,0x22,0xd5,0xdf,0x08,
        0x6d,0x0f,0xf0,0xb1,0x10,0xfb,0xd9,0xd2,
        0x1b,0xb4,0xfc,0x71,0x63,0xaf,0x34,0xd0,
        0x82,0x86,0xa2,0xe8,0x46,0xf6,0xbe,0x03,
    };
    assert(memcmp(digest, expected, 32) == 0);

    unlink(path);
    printf("test_sha256_file: PASS\n");
}

// -- Minimal cache file builder for load/probe tests ------------------

// Splitmix64 mirror of accel_cache.cpp's accel_hash64.
static inline uint64_t splitmix64(uint64_t x) {
    x += 0x9E3779B97F4A7C15ULL;
    x = (x ^ (x >> 30)) * 0xBF58476D1CE4E5B9ULL;
    x = (x ^ (x >> 27)) * 0x94D049BB133111EBULL;
    return x ^ (x >> 31);
}

struct BuilderEntry { uint64_t kmer; int64_t k, l, s; };

// Write a minimal cache file with one level at L=1 and a terminal table.
// `l1_entries` is used at both the L=1 table AND the terminal table for
// simplicity (so the terminal hit returns the same thing). Uses table
// size = next_pow2(2*n).
static void write_minimal_cache(const char *path,
                                const uint8_t ref_sha256[32],
                                const BuilderEntry *l1_entries,
                                size_t n) {
    // Pow2 table size.
    uint64_t table_size = 4;
    while (table_size < 2 * n) table_size <<= 1;
    uint64_t mask = table_size - 1;

    // Build probe tables (both copies — L=1 and terminal).
    uint64_t *probe = (uint64_t *)calloc(table_size, sizeof(uint64_t));
    for (size_t i = 0; i < n; ++i) {
        uint64_t s = splitmix64(l1_entries[i].kmer) & mask;
        while (probe[s] != 0) s = (s + 1) & mask;
        probe[s] = i + 1;
    }

    FILE *f = fopen(path, "wb");
    assert(f);
    uint64_t magic = ACCEL_MAGIC;              fwrite(&magic, 8, 1, f);
    uint32_t ver   = ACCEL_VERSION;            fwrite(&ver, 4, 1, f);
    fwrite(ref_sha256, 32, 1, f);
    uint32_t kmin = 1;                         fwrite(&kmin, 4, 1, f);
    uint32_t kmax = 1;                         fwrite(&kmax, 4, 1, f);
    uint32_t tier_count = 1;                   fwrite(&tier_count, 4, 1, f);
    uint32_t collapse = 1;                     fwrite(&collapse, 4, 1, f);
    uint8_t reserved[16] = {};                 fwrite(reserved, 16, 1, f);

    uint8_t tier_id = ACCEL_TIER_B;            fwrite(&tier_id, 1, 1, f);
    uint8_t tier_pad[7] = {};                  fwrite(tier_pad, 7, 1, f);

    // Level L=1
    uint64_t n_e = n;                          fwrite(&n_e, 8, 1, f);
    uint64_t t_m = mask;                       fwrite(&t_m, 8, 1, f);
    fwrite(probe, sizeof(uint64_t), table_size, f);
    for (size_t i = 0; i < n; ++i) {
        AccelEntry e = {l1_entries[i].kmer,
                        l1_entries[i].k, l1_entries[i].l, l1_entries[i].s};
        fwrite(&e, sizeof(e), 1, f);
    }

    // Terminal (same contents for test simplicity).
    fwrite(&n_e, 8, 1, f);
    fwrite(&t_m, 8, 1, f);
    fwrite(probe, sizeof(uint64_t), table_size, f);
    for (size_t i = 0; i < n; ++i) {
        AccelEntry e = {l1_entries[i].kmer,
                        l1_entries[i].k, l1_entries[i].l, l1_entries[i].s};
        fwrite(&e, sizeof(e), 1, f);
    }
    fclose(f);
    free(probe);
}

static void test_load_and_probe() {
    const char *path = "/tmp/accel_minimal.cache";
    uint8_t ref[32];
    memset(ref, 0xAB, 32);

    BuilderEntry entries[] = {
        // A (=0), C (=1), G (=2)
        {0, 100, 200, 300},
        {1, 400, 500, 600},
        {2, 700, 800, 900},
    };
    write_minimal_cache(path, ref, entries, 3);

    AccelCache cache;
    bool ok = cache.load(path, ref);
    assert(ok);
    assert(cache.enabled());
    assert(cache.k_min() == 1);
    assert(cache.k_max() == 1);

    AccelLookupResult r = cache.probe(1, 0);
    assert(r.hit); assert(r.sa_k == 100); assert(r.sa_l == 200); assert(r.sa_s == 300);

    r = cache.probe(1, 1);
    assert(r.hit); assert(r.sa_k == 400);

    r = cache.probe(1, 2);
    assert(r.hit); assert(r.sa_k == 700);

    r = cache.probe(1, 99);
    assert(!r.hit);

    r = cache.probe_terminal(0);
    assert(r.hit); assert(r.sa_k == 100);

    // Out-of-range L.
    r = cache.probe(2, 0);
    assert(!r.hit);
    r = cache.probe(0, 0);
    assert(!r.hit);

    unlink(path);
    printf("test_load_and_probe: PASS\n");
}

static void test_load_rejects_bad_magic() {
    const char *path = "/tmp/accel_bad_magic.cache";
    FILE *f = fopen(path, "wb");
    uint64_t bad_magic = 0xDEADBEEFDEADBEEFULL;
    fwrite(&bad_magic, 8, 1, f);
    fclose(f);
    AccelCache c;
    uint8_t ref[32] = {};
    assert(!c.load(path, ref));
    assert(!c.enabled());
    unlink(path);
    printf("test_load_rejects_bad_magic: PASS\n");
}

static void test_load_rejects_checksum_mismatch() {
    const char *path = "/tmp/accel_bad_ref.cache";
    uint8_t good[32]; memset(good, 0xAA, 32);
    BuilderEntry e[] = {{0, 1, 2, 3}};
    write_minimal_cache(path, good, e, 1);

    uint8_t wrong[32]; memset(wrong, 0xBB, 32);
    AccelCache c;
    assert(!c.load(path, wrong));
    assert(!c.enabled());
    unlink(path);
    printf("test_load_rejects_checksum_mismatch: PASS\n");
}

int main() {
    test_pack_forward_and_revcomp();
    test_canonical();
    test_pack_rejects_n();
    test_sha256_file();
    test_load_and_probe();
    test_load_rejects_bad_magic();
    test_load_rejects_checksum_mismatch();
    printf("\nALL TESTS PASS\n");
    return 0;
}
