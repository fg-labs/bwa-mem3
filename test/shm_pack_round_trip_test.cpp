/* Round-trip check: pack the phix index into a self-describing buffer via
 * bwa_shm_pack_from_disk, then verify each section's bytes match what the
 * normal disk loader produces. Run via test/shm_pack_round_trip_test.sh,
 * which builds the phix index first.
 *
 * Usage: shm_pack_round_trip_test <prefix>
 *   <prefix> is the path passed to `bwa-mem3 index`, e.g. test/fixtures/phix.fa
 */

#include "bwa_shm.h"
#include "FMI_search.h"
#include "bntseq.h"
#include "utils.h"

#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#define CHECK(cond) do { \
    if (!(cond)) { \
        std::fprintf(stderr, "CHECK failed at %s:%d: %s\n", __FILE__, __LINE__, #cond); \
        std::fflush(stderr); \
        std::abort(); \
    } \
} while (0)

#define CHECK_EQ(a, b) do { \
    auto _a = (a); auto _b = (b); \
    if (_a != _b) { \
        std::fprintf(stderr, "CHECK_EQ failed at %s:%d: %s (=%lld) != %s (=%lld)\n", \
                     __FILE__, __LINE__, #a, (long long)_a, #b, (long long)_b); \
        std::fflush(stderr); \
        std::abort(); \
    } \
} while (0)

static uint8_t *slurp(const char *path, int64_t *len_out) {
    FILE *fp = std::fopen(path, "rb");
    CHECK(fp != NULL);
    std::fseek(fp, 0, SEEK_END);
    int64_t n = std::ftell(fp);
    std::rewind(fp);
    uint8_t *buf = (uint8_t *)std::malloc((size_t)n);
    CHECK(buf != NULL);
    err_fread_noeof(buf, 1, (size_t)n, fp);
    std::fclose(fp);
    *len_out = n;
    return buf;
}

int main(int argc, char *argv[]) {
    if (argc < 2) {
        std::fprintf(stderr, "Usage: %s <prefix>\n", argv[0]);
        return 2;
    }
    const char *prefix = argv[1];

    /* 1. Pack from disk. */
    uint8_t *buf = NULL;
    uint64_t buf_len = 0;
    CHECK_EQ(bwa_shm_pack_from_disk(prefix, &buf, &buf_len), 0);
    CHECK(buf != NULL);
    CHECK(buf_len > sizeof(bwa_shm_header_t));

    /* 2. Header sanity. */
    bwa_shm_header_t hdr;
    std::memcpy(&hdr, buf, sizeof(hdr));
    CHECK_EQ(hdr.magic,      BWA_SHM_MAGIC);
    CHECK_EQ(hdr.version,    BWA_SHM_VERSION);
    CHECK_EQ(hdr.n_sections, 10u);
    CHECK_EQ(hdr.total_size, buf_len);

    /* 3. Independently load the same index from disk for comparison. */
    FMI_search ref(prefix);
    ref.load_index();

    char zer_path[4096];
    std::snprintf(zer_path, sizeof(zer_path), "%s.0123", prefix);
    int64_t ref_string_len = 0;
    uint8_t *ref_string = slurp(zer_path, &ref_string_len);

    /* 4. Walk sections and compare. */
    auto section = [&](uint32_t kind, uint64_t *off, uint64_t *sz) {
        CHECK_EQ(bwa_shm_section_find(buf, kind, off, sz), 0);
    };

    uint64_t off = 0, sz = 0;

    /* FMI scalars: 8 (refseq_len) + 40 (count[5]) + 8 (sentinel_index) = 56. */
    section(BWA_SHM_SEC_FMI_SCALARS, &off, &sz);
    CHECK_EQ(sz, 56ull);
    int64_t got_refseq = 0, got_count[5] = {0}, got_sentinel = 0;
    std::memcpy(&got_refseq,   buf + off,            sizeof(int64_t));
    std::memcpy(got_count,     buf + off + 8,        sizeof(int64_t) * 5);
    std::memcpy(&got_sentinel, buf + off + 8 + 40,   sizeof(int64_t));
    CHECK_EQ(got_refseq,   ref.reference_seq_len);
    CHECK_EQ(got_sentinel, ref.sentinel_index);
    for (int i = 0; i < 5; ++i) CHECK_EQ(got_count[i], ref.count_data()[i]);

    /* CP_OCC. */
    section(BWA_SHM_SEC_FMI_CP_OCC, &off, &sz);
    CHECK_EQ((int64_t)sz, ref.cp_occ_size_bytes());
    CHECK_EQ(std::memcmp(buf + off, ref.cp_occ_data(), (size_t)sz), 0);

    /* SA_MS. */
    section(BWA_SHM_SEC_FMI_SA_MS, &off, &sz);
    CHECK_EQ((int64_t)sz, ref.sa_ms_byte_size_bytes());
    CHECK_EQ(std::memcmp(buf + off, ref.sa_ms_byte_data(), (size_t)sz), 0);

    /* SA_LS. */
    section(BWA_SHM_SEC_FMI_SA_LS, &off, &sz);
    CHECK_EQ((int64_t)sz, ref.sa_ls_word_size_bytes());
    CHECK_EQ(std::memcmp(buf + off, ref.sa_ls_word_data(), (size_t)sz), 0);

    /* BNS struct. The packed copy contains the SAME bntseq_t bytes as the
     * loaded version, including its (now-zero) fp_pac/fp_dict and its
     * pointer fields. We don't compare the bytes byte-for-byte because the
     * pointer fields will differ (the loader points anns/ambs into freshly-
     * malloc'd memory, while the packed copy serializes them as offsets in
     * the segment, but the bns struct itself was memcpy'd before that
     * patching). Just check size and a couple of stable scalars. */
    section(BWA_SHM_SEC_BNS_STRUCT, &off, &sz);
    CHECK_EQ(sz, sizeof(bntseq_t));
    bntseq_t got_bns;
    std::memcpy(&got_bns, buf + off, sizeof(got_bns));
    CHECK_EQ(got_bns.l_pac,    ref.idx->bns->l_pac);
    CHECK_EQ(got_bns.n_seqs,   ref.idx->bns->n_seqs);
    CHECK_EQ(got_bns.n_holes,  ref.idx->bns->n_holes);
    CHECK_EQ(got_bns.seed,     ref.idx->bns->seed);

    /* AMBS — direct byte compare. */
    section(BWA_SHM_SEC_BNS_AMBS, &off, &sz);
    CHECK_EQ((int64_t)sz, (int64_t)ref.idx->bns->n_holes * (int64_t)sizeof(bntamb1_t));
    CHECK_EQ(std::memcmp(buf + off, ref.idx->bns->ambs, (size_t)sz), 0);

    /* ANNS — byte compare of the struct array. The pointers inside (name/anno)
     * will differ because the loader's anns point at heap strings while the
     * packed anns are direct memcpy of the loader's anns (still pointing at
     * loader heap). Compare scalar fields only. */
    section(BWA_SHM_SEC_BNS_ANNS, &off, &sz);
    CHECK_EQ((int64_t)sz, (int64_t)ref.idx->bns->n_seqs * (int64_t)sizeof(bntann1_t));
    const bntann1_t *got_anns = (const bntann1_t *)(buf + off);
    for (int i = 0; i < ref.idx->bns->n_seqs; ++i) {
        CHECK_EQ(got_anns[i].offset,   ref.idx->bns->anns[i].offset);
        CHECK_EQ(got_anns[i].len,      ref.idx->bns->anns[i].len);
        CHECK_EQ(got_anns[i].n_ambs,   ref.idx->bns->anns[i].n_ambs);
        CHECK_EQ(got_anns[i].gi,       ref.idx->bns->anns[i].gi);
        CHECK_EQ(got_anns[i].is_alt,   ref.idx->bns->anns[i].is_alt);
    }

    /* NAMES — concatenated name\0 anno\0 strings. */
    section(BWA_SHM_SEC_BNS_NAMES, &off, &sz);
    {
        const char *p = (const char *)(buf + off);
        const char *end = p + sz;
        for (int i = 0; i < ref.idx->bns->n_seqs; ++i) {
            CHECK(p < end);
            CHECK_EQ(std::strcmp(p, ref.idx->bns->anns[i].name), 0);
            p += std::strlen(p) + 1;
            CHECK(p < end);
            CHECK_EQ(std::strcmp(p, ref.idx->bns->anns[i].anno), 0);
            p += std::strlen(p) + 1;
        }
        CHECK(p == end);
    }

    /* PAC. */
    section(BWA_SHM_SEC_PAC, &off, &sz);
    CHECK_EQ((int64_t)sz, (int64_t)(ref.idx->bns->l_pac / 4 + 1));
    CHECK_EQ(std::memcmp(buf + off, ref.idx->pac, (size_t)sz), 0);

    /* REF_STRING (.0123). */
    section(BWA_SHM_SEC_REF_STRING, &off, &sz);
    CHECK_EQ((int64_t)sz, ref_string_len);
    CHECK_EQ(std::memcmp(buf + off, ref_string, (size_t)sz), 0);

    /* Cleanup. */
    std::free(ref_string);
    std::free(buf);

    std::printf("shm_pack_round_trip_test: OK\n");
    return 0;
}
