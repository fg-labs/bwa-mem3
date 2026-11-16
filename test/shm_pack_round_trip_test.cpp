/* Round-trip check: pack the phix index into a self-describing buffer via
 * bwa_shm_pack_from_disk, then verify each section's bytes match what the
 * normal disk loader produces. Run via test/shm_pack_round_trip_test.sh,
 * which builds the phix index first.
 *
 * Calls bwa_shm_destroy before and after the orphan check. That drops every
 * registered staged index this user can unlink, plus /bwactl and /bwactl_lock,
 * on this host (unregistered /bwaidx-* objects are left alone). So the test
 * first reads the registry and fails, dropping nothing, if it lists any index
 * other than the fixture. An index staged concurrently after that check can
 * still be dropped; run it where nothing else stages an index.
 *
 * Usage: shm_pack_round_trip_test <prefix>
 *   <prefix> is the path passed to `bwa-mem3 index`, e.g. test/fixtures/phix.fa
 */

#include "bwa_shm.h"
#include "FMI_search.h"
#include "bntseq.h"
#include "utils.h"
#include "framework/shm_registry_guard.h"
#include "framework/temp_dir.h"

#include <cerrno>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fcntl.h>
#include <sys/mman.h>
#include <string>
#include <sys/stat.h>
#include <sys/wait.h>
#include <unistd.h>
#include <vector>

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

/* The .bwt.2bit.64 header before CP_OCC: reference_seq_len + count[5]. */
static const int64_t BWT_2BIT_HEADER_BYTES = (int64_t)sizeof(int64_t) * 6;

/* The raw one-byte-per-sample SA high parts of <prefix>.bwt.2bit.64, which follow
 * the header and CP_OCC. */
static std::vector<uint8_t> read_disk_sa_ms(const char *prefix, int64_t cp_occ_bytes, int64_t n) {
    const std::string path = std::string(prefix) + ".bwt.2bit.64";
    FILE *f = std::fopen(path.c_str(), "rb");
    CHECK(f != NULL);
    CHECK_EQ(fseeko(f, (off_t)(BWT_2BIT_HEADER_BYTES + cp_occ_bytes), SEEK_SET), 0);
    std::vector<uint8_t> v((size_t)n);
    CHECK_EQ(std::fread(v.data(), 1, v.size(), f), v.size());
    std::fclose(f);
    return v;
}

/* Overwrite SA high-part byte i of <prefix>.bwt.2bit.64 in place. */
static void poke_disk_sa_ms(const char *prefix, int64_t cp_occ_bytes, int64_t i, uint8_t b) {
    const std::string path = std::string(prefix) + ".bwt.2bit.64";
    FILE *f = std::fopen(path.c_str(), "r+b");
    CHECK(f != NULL);
    CHECK_EQ(fseeko(f, (off_t)(BWT_2BIT_HEADER_BYTES + cp_occ_bytes + i), SEEK_SET), 0);
    CHECK_EQ(std::fwrite(&b, 1, 1, f), (size_t)1);
    CHECK_EQ(std::fclose(f), 0);
}

/* Sample i of a packed SA_MS array, decoded independently of FMI_search's accessor. */
static int64_t packed_at(const uint8_t *sec, int log2, int64_t i) {
    const int bits = 1 << log2;
    const int64_t per_word = 64 / bits;
    uint64_t w;
    std::memcpy(&w, sec + (size_t)(i / per_word) * sizeof(uint64_t), sizeof(w));
    return (int64_t)((w >> ((i % per_word) * bits)) & ((1ull << bits) - 1));
}

static void copy_file(const std::string &from, const std::string &to) {
    FILE *in = std::fopen(from.c_str(), "rb");
    CHECK(in != NULL);
    FILE *out = std::fopen(to.c_str(), "wb");
    CHECK(out != NULL);
    char b[1 << 16];
    size_t r;
    while ((r = std::fread(b, 1, sizeof(b), in)) > 0) CHECK_EQ(std::fwrite(b, 1, r, out), r);
    std::fclose(in);
    CHECK_EQ(std::fclose(out), 0);
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

    /* 4. Walk sections and compare. */
    auto section = [&](uint32_t kind, uint64_t *off, uint64_t *sz) {
        CHECK_EQ(bwa_shm_section_find(buf, kind, off, sz), 0);
    };

    uint64_t off = 0, sz = 0;

    /* FMI scalars: 8 (refseq_len) + 40 (count[5]) + 8 (sentinel_index)
     * + 8 (sa_compx) = 64. */
    section(BWA_SHM_SEC_FMI_SCALARS, &off, &sz);
    CHECK_EQ(sz, 64ull);
    int64_t got_refseq = 0, got_count[5] = {0}, got_sentinel = 0, got_sa_compx = 0;
    std::memcpy(&got_refseq,   buf + off,            sizeof(int64_t));
    std::memcpy(got_count,     buf + off + 8,        sizeof(int64_t) * 5);
    std::memcpy(&got_sentinel, buf + off + 8 + 40,   sizeof(int64_t));
    std::memcpy(&got_sa_compx, buf + off + 8 + 40 + 8, sizeof(int64_t));
    CHECK_EQ(got_refseq,   ref.reference_seq_len);
    CHECK_EQ(got_sentinel, ref.sentinel_index);
    CHECK_EQ(got_sa_compx, ref.sa_compx);
    for (int i = 0; i < 5; ++i) CHECK_EQ(got_count[i], ref.count_data()[i]);

    /* CP_OCC. */
    section(BWA_SHM_SEC_FMI_CP_OCC, &off, &sz);
    CHECK_EQ((int64_t)sz, ref.cp_occ_size_bytes());
    CHECK_EQ(std::memcmp(buf + off, ref.cp_occ_data(), (size_t)sz), 0);

    /* SA_MS. */
    section(BWA_SHM_SEC_FMI_SA_MS, &off, &sz);
    CHECK_EQ((int64_t)sz, ref.sa_hi_size_bytes());
    CHECK_EQ(std::memcmp(buf + off, ref.sa_hi_data(), (size_t)sz), 0);
    /* Oracle: the loader and the stager share sa_hi_pack, so also check both packed
     * arrays against the raw bytes on disk. */
    const int hi_log2 = FMI_search::sa_hi_log2_for(ref.reference_seq_len);
    {
        const int64_t n = ref.sa_sample_count();
        const std::vector<uint8_t> raw = read_disk_sa_ms(prefix, ref.cp_occ_size_bytes(), n);
        for (int64_t i = 0; i < n; ++i) {
            CHECK_EQ(packed_at((const uint8_t *)ref.sa_hi_data(), hi_log2, i), (int64_t)raw[(size_t)i]);
            CHECK_EQ(packed_at(buf + off, hi_log2, i), (int64_t)raw[(size_t)i]);
        }
    }

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
    /* pos2rid_bucket is derived and process-private: the packer must stage a
     * NULL, not the writer's heap address. The attach path rebuilds it locally.
     * Assert the source table was actually built first — otherwise a loader
     * regression that stopped building it would satisfy the NULL check
     * vacuously and this test would stop guarding the packer at all. */
    CHECK(ref.idx->bns->pos2rid_bucket != NULL);
    CHECK(got_bns.pos2rid_bucket == NULL);

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

    /* REF_STRING (.0123) is NEVER staged: the section slot is present (so
     * n_sections stays 10) but always zero length. `mem` pac-fetches the
     * original reference from `.pac` on demand, so the unpacked `.0123` is
     * neither staged nor required on disk. */
    section(BWA_SHM_SEC_REF_STRING, &off, &sz);
    CHECK_EQ(sz, 0ull);

    /* Cleanup. */
    std::free(buf);

    /* Chunked loads and corrupt high parts, on a scratch copy of the index. phiX needs
     * 1 bit per sample, so its high parts are all zero: set every third one to 1 (in
     * range) so a misplaced sample shows, and lower the read chunk to 128 samples so
     * the loader and the stager both cross many chunk boundaries. The staged segment is
     * then attached (load_index_from_shm) and read back through the resolver. Then set
     * one byte in a later chunk to a value that does not fit in 1 bit (2, and the
     * sign-bit bytes 0x80 and 0xff): the stager, with and without densify, must fail and
     * the disk loader must exit with an error. */
    {
        CHECK_EQ(hi_log2, 0);
        /* Removed with its files when the block ends, whether or not a check failed. */
        const bwa_tests::TempDir tdir("sa_hi_chunk");
        const std::string &dir = tdir.path();
        /* A basename unique to this run: shm segments are keyed by the index basename, so a
         * staged phix.fa must not be attached in place of this modified copy. */
        const std::string tprefix = dir + "/" + dir.substr(dir.rfind('/') + 1) + ".fa";
        static const char *const exts[] = {".amb", ".ann", ".bwt.2bit.64", ".pac"};
        for (const char *e : exts) copy_file(std::string(prefix) + e, tprefix + e);

        const int64_t n = ref.sa_sample_count();
        const int64_t cp_bytes = ref.cp_occ_size_bytes();
        const int64_t chunk = 128;
        CHECK(n > 4 * chunk);
        for (int64_t i = 0; i < n; i += 3) poke_disk_sa_ms(tprefix.c_str(), cp_bytes, i, 1);
        const std::vector<uint8_t> raw = read_disk_sa_ms(tprefix.c_str(), cp_bytes, n);
        FMI_search::set_sa_hi_read_chunk(chunk);

        FMI_search t(tprefix.c_str());
        t.load_index();
        const uint8_t *loaded = (const uint8_t *)t.sa_hi_data();
        for (int64_t i = 0; i < n; ++i) CHECK_EQ(packed_at(loaded, hi_log2, i), (int64_t)raw[(size_t)i]);
        uint8_t *tbuf = NULL;
        uint64_t tlen = 0, toff = 0, tsz = 0;
        CHECK_EQ(bwa_shm_pack_from_disk(tprefix.c_str(), &tbuf, &tlen), 0);
        CHECK_EQ(bwa_shm_section_find(tbuf, BWA_SHM_SEC_FMI_SA_MS, &toff, &tsz), 0);
        CHECK_EQ((int64_t)tsz, t.sa_hi_size_bytes());
        CHECK_EQ(std::memcmp(tbuf + toff, loaded, (size_t)tsz), 0);
        for (int64_t i = 0; i < n; ++i)
            CHECK_EQ(packed_at(tbuf + toff, hi_log2, i), (int64_t)raw[(size_t)i]);
        {
            /* Attach a copy of the staged index as `mem` would (the object unmaps it) and
             * resolve every sampled row: each must match the disk loader's coordinate. */
            void *m = mmap(NULL, (size_t)tlen, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
            CHECK(m != MAP_FAILED);
            std::memcpy(m, tbuf, (size_t)tlen);
            FMI_search a(tprefix.c_str());
            a.load_index_from_shm((uint8_t *)m, (size_t)tlen, (size_t)tlen);
            CHECK_EQ((const uint8_t *)a.sa_hi_data(), (const uint8_t *)m + toff);
            for (int64_t i = 0; i < n; ++i) {
                const int64_t row = i << t.sa_compx;
                CHECK_EQ(a.get_sa_entry_compressed(row, 0), t.get_sa_entry_compressed(row, 0));
                CHECK_EQ(a.get_sa_entry_compressed(row, 0) >> 32, (int64_t)raw[(size_t)i]);
            }
        }
        std::free(tbuf);

        const int64_t bad_at = 2 * chunk + 45;
        const std::string err_path = dir + "/load.err";
        for (int b : {2, 0x80, 0xff}) {
            std::printf("shm_pack_round_trip_test: corrupt high byte 0x%02x\n", b);
            std::fflush(stdout);
            poke_disk_sa_ms(tprefix.c_str(), cp_bytes, bad_at, (uint8_t)b);
            tbuf = NULL;
            CHECK(bwa_shm_pack_from_disk(tprefix.c_str(), &tbuf, &tlen) != 0);
            CHECK(tbuf == NULL);
            /* The densify path (`shm -u`) rejects it as well, by returning an error. */
            bwa_shm_layout_t dl;
            CHECK_EQ(bwa_shm_compute(tprefix.c_str(), &dl, /*bns_only=*/false,
                                     (int)t.sa_compx - 1, /*n_threads=*/2), 0);
            CHECK(dl.sa_compx < dl.disk_sa_compx);
            uint8_t *dbuf2 = (uint8_t *)std::calloc(1, dl.total_size);
            CHECK(dbuf2 != NULL);
            CHECK(bwa_shm_pack_into(&dl, dbuf2) != 0);
            std::free(dbuf2);
            bwa_shm_layout_free(&dl);
            std::fflush(NULL);
            const pid_t pid = fork();
            CHECK(pid >= 0);
            if (pid == 0) {
                /* dup2, not freopen: a reopened stderr may be fully buffered, and the
                 * loader's fatal error aborts without flushing it. */
                const int efd = open(err_path.c_str(), O_WRONLY | O_CREAT | O_TRUNC, 0600);
                if (efd < 0 || dup2(efd, STDERR_FILENO) < 0) _exit(3);
                FMI_search c(tprefix.c_str());
                c.load_index();
                _exit(0);   /* loaded a corrupt index: the parent fails the test */
            }
            int status = 0;
            CHECK_EQ(waitpid(pid, &status, 0), pid);
            CHECK(!(WIFEXITED(status) && WEXITSTATUS(status) == 0));
            FILE *ef = std::fopen(err_path.c_str(), "r");
            CHECK(ef != NULL);
            std::string err;
            char eb[4096];
            size_t r;
            while ((r = std::fread(eb, 1, sizeof(eb), ef)) > 0) err.append(eb, r);
            std::fclose(ef);
            CHECK(err.find("SA sample high part out of range") != std::string::npos);
        }
        FMI_search::set_sa_hi_read_chunk((int64_t)1 << 20);   /* the default */
        std::printf("shm_pack_round_trip_test: chunked load and corrupt high part OK\n");
    }

    /* Densify parity (only when an oracle prefix is given): pack THIS index
     * densified via the `shm -u` path to the oracle's rate, and confirm the
     * staged SA sample arrays are byte-identical to the oracle -- an index
     * physically built at that denser rate (index -u), an INDEPENDENT SA
     * construction, not a self-consistency check. This is the only automated
     * coverage of the shm densify path (bwa_shm_compute sizing +
     * bwa_shm_pack_into's densify branch feeding densify_sa_into). */
    if (argc >= 3) {
        const char *oracle_prefix = argv[2];
        FMI_search oracle(oracle_prefix);
        oracle.load_index();
        CHECK(oracle.sa_compx < ref.sa_compx);   /* oracle must be strictly denser */

        bwa_shm_layout_t layout;
        CHECK_EQ(bwa_shm_compute(prefix, &layout, /*bns_only=*/false,
                                 (int)oracle.sa_compx, /*n_threads=*/1), 0);
        CHECK_EQ(layout.sa_compx, oracle.sa_compx);         /* densify actually took */
        CHECK(layout.sa_compx < layout.disk_sa_compx);
        uint8_t *dbuf = (uint8_t *)std::malloc(layout.total_size);
        CHECK(dbuf != NULL);
        CHECK_EQ(bwa_shm_pack_into(&layout, dbuf), 0);
        bwa_shm_layout_free(&layout);

        uint64_t doff = 0, dsz = 0;
        CHECK_EQ(bwa_shm_section_find(dbuf, BWA_SHM_SEC_FMI_SA_MS, &doff, &dsz), 0);
        CHECK_EQ((int64_t)dsz, oracle.sa_hi_size_bytes());
        CHECK_EQ(std::memcmp(dbuf + doff, oracle.sa_hi_data(), (size_t)dsz), 0);
        CHECK_EQ(bwa_shm_section_find(dbuf, BWA_SHM_SEC_FMI_SA_LS, &doff, &dsz), 0);
        CHECK_EQ((int64_t)dsz, oracle.sa_ls_word_size_bytes());
        CHECK_EQ(std::memcmp(dbuf + doff, oracle.sa_ls_word_data(), (size_t)dsz), 0);
        std::free(dbuf);
        std::printf("shm_pack_round_trip_test: densify parity OK "
                    "(shm -u staged 1/%d == index -u %lld)\n",
                    1 << (int)oracle.sa_compx, (long long)oracle.sa_compx);
    }

    /* A complete-looking segment without a registry entry is not a successful
     * stage. Preserve it rather than adopting it or replacing a live mapping.
     * (A segment whose magic is still 0 is an interrupted stage; the stager
     * reclaims that one, which test/shm_stage_test.cpp covers.) */
    const char *base = std::strrchr(prefix, '/');
    base = base ? base + 1 : prefix;
    char segment[256];
    CHECK(std::snprintf(segment, sizeof(segment), "%s%s", BWA_SHM_IDX_PREFIX, base)
          < (int)sizeof(segment));
    if (!bwa_tests::shm_registry_holds_only(BWA_SHM_CTL_NAME, base)) {
        std::fprintf(stderr, "shm_pack_round_trip_test: FAIL: the shm registry %s lists an "
                     "index other than %s, or could not be read; refusing to run "
                     "bwa_shm_destroy, which would drop every staged index. Run "
                     "`bwa-mem3 shm -l` to inspect it and `bwa-mem3 shm -d` to clear it, "
                     "or use another host, then rerun.\n", BWA_SHM_CTL_NAME, base);
        return 1;
    }
    CHECK_EQ(bwa_shm_destroy(), 0);
    const int fd = shm_open(segment, O_CREAT | O_RDWR | O_EXCL, 0600);
    CHECK(fd >= 0);
    bwa_shm_header_t orphan_hdr;
    std::memset(&orphan_hdr, 0, sizeof(orphan_hdr));
    orphan_hdr.magic      = BWA_SHM_MAGIC;
    orphan_hdr.version    = BWA_SHM_VERSION;
    orphan_hdr.total_size = sizeof(orphan_hdr);
    std::memcpy(&orphan_hdr.reserved[2], "unregistered", sizeof("unregistered"));
    char marker[sizeof(orphan_hdr)];
    std::memcpy(marker, &orphan_hdr, sizeof(marker));
    CHECK_EQ(ftruncate(fd, sizeof(marker)), 0);
    void *orphan = mmap(NULL, sizeof(marker), PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
    CHECK(orphan != MAP_FAILED);
    std::memcpy(orphan, marker, sizeof(marker));
    const int stage_rc = bwa_shm_stage(prefix);
    char saved[sizeof(marker)] = {};
    std::memcpy(saved, orphan, sizeof(saved));
    CHECK_EQ(munmap(orphan, sizeof(marker)), 0);
    CHECK_EQ(close(fd), 0);
    /* Reopen by name and read the bytes back: a regression that replaced the
     * orphan with a new object of the same name must not pass. Check the size
     * before mapping so a shorter replacement fails instead of faulting; use
     * >= because macOS rounds a shm object's size up to the page size. */
    const int preserved_fd = shm_open(segment, O_RDONLY, 0);
    bool preserved = preserved_fd >= 0;
    char reopened[sizeof(marker)] = {};
    if (preserved) {
        struct stat st;
        preserved = fstat(preserved_fd, &st) == 0 && st.st_size >= (off_t)sizeof(marker);
        if (preserved) {
            void *view = mmap(NULL, sizeof(marker), PROT_READ, MAP_SHARED, preserved_fd, 0);
            preserved = view != MAP_FAILED;
            if (preserved) {
                std::memcpy(reopened, view, sizeof(reopened));
                CHECK_EQ(munmap(view, sizeof(marker)), 0);
            }
        }
        CHECK_EQ(close(preserved_fd), 0);
    }
    CHECK_EQ(shm_unlink(segment), 0);
    CHECK_EQ(bwa_shm_destroy(), 0);
    CHECK(stage_rc != 0);
    CHECK(preserved);
    CHECK_EQ(std::memcmp(saved, marker, sizeof(marker)), 0);
    CHECK_EQ(std::memcmp(reopened, marker, sizeof(marker)), 0);

    /* Registered restaging still succeeds, attaches, and is removed by -d. */
    CHECK_EQ(bwa_shm_stage(prefix), 0);
    CHECK_EQ(bwa_shm_test(prefix), 1);
    CHECK_EQ(bwa_shm_stage(prefix), 0);
    size_t attached_len = 0, map_len = 0;
    uint8_t *attached = bwa_shm_attach(prefix, &attached_len, &map_len);
    CHECK(attached != NULL);
    CHECK(attached_len > sizeof(bwa_shm_header_t));
    /* The mapping is the whole POSIX object, whatever size the platform gave it. */
    const int staged_fd = shm_open(segment, O_RDONLY, 0);
    CHECK(staged_fd >= 0);
    struct stat staged_st;
    CHECK_EQ(fstat(staged_fd, &staged_st), 0);
    CHECK_EQ(close(staged_fd), 0);
    CHECK_EQ(map_len, (size_t)staged_st.st_size);
    CHECK(map_len >= attached_len);
    bwa_shm_detach(attached, map_len);
    CHECK_EQ(bwa_shm_destroy(), 0);
    errno = 0;
    const int leftover_fd = shm_open(segment, O_RDONLY, 0);
    const int leftover_errno = errno;
    if (leftover_fd >= 0) close(leftover_fd);
    CHECK(leftover_fd < 0);
    CHECK_EQ(leftover_errno, ENOENT);
    CHECK_EQ(bwa_shm_test(prefix), 0);

    std::printf("shm_pack_round_trip_test: OK\n");
    return 0;
}
