// In-memory FM-index build: the doubled-pac packer, the MemSink, libsais_build_fm_memory and
// FMI_search::adopt_built must produce exactly what the file-based build writes and loads.

#include "doctest/doctest.h"

#include "../framework/temp_dir.h"
#include "../framework/temp_index.h"
#include "../../src/fm_index_writer.h"
#include "../../src/bntseq.h"
#include "../../src/libsais_build.h"
#include "../../src/FMI_search.h"
#include "../../src/aln/aln_bwt.h"

#include <cstdint>
#include <algorithm>
#include <cstdio>
#include <cstring>
#include <memory>
#include <string>
#include <vector>

namespace {

// 2-bit pack bases (0-3), first base in the high bits, unused low bits of a partial final byte
// zero: the .pac byte layout, built here one base at a time as an oracle for pack_doubled.
std::vector<uint8_t> pack_bases(const std::vector<uint8_t> &bases)
{
    std::vector<uint8_t> p((bases.size() + 3) / 4, 0);
    for (size_t i = 0; i < bases.size(); ++i) p[i / 4] |= (uint8_t)(bases[i] << (2 * (3 - i % 4)));
    return p;
}

std::vector<uint8_t> read_file(const std::string &path)
{
    std::vector<uint8_t> v;
    FILE *f = std::fopen(path.c_str(), "rb");
    REQUIRE(f != nullptr);
    uint8_t buf[4096];
    size_t n;
    while ((n = std::fread(buf, 1, sizeof buf, f)) > 0) v.insert(v.end(), buf, buf + n);
    std::fclose(f);
    return v;
}

// Random ACGT bases from a fixed seed.
std::string random_bases(size_t n, uint32_t seed)
{
    static const char kAcgt[] = "ACGT";
    std::string out(n, 'A');
    for (char &c : out) {
        seed = seed * 1664525u + 1013904223u;
        c = kAcgt[(seed >> 24) & 3];
    }
    return out;
}

// A few-kilobase synthetic reference: three contigs, one with an N run (an ambiguity hole,
// which the .pac stores as a random base), with lengths that are not multiples of 4.
std::vector<bwa_tests::Contig> multi_contig_ref()
{
    return {{"a", random_bases(1001, 1u)},
            {"b", random_bases(700, 2u) + "NNNNN" + random_bases(300, 3u)},
            {"c", random_bases(67, 4u)}};
}

} // namespace

TEST_CASE("libsais_mem: pack_doubled packs forward then reverse complement"
          * doctest::test_suite("unit/libsais_mem"))
{
    // ACG -> doubled ACG CGT = 0 1 2 1 2 3 -> 00011001 1011(0000).
    const std::vector<uint8_t> fwd = pack_bases({0, 1, 2});
    REQUIRE(fwd.size() == 1);
    REQUIRE(fwd[0] == 0x18);
    PackedText view(fwd.data(), 3);
    std::unique_ptr<uint8_t[]> mem = pack_doubled(view, 3, 1);
    REQUIRE(mem != nullptr);
    CHECK(mem[0] == 0x19);
    CHECK(mem[1] == 0xB0);
}

TEST_CASE("libsais_mem: pack_doubled and write_doubled_pac equal an independent packing"
          * doctest::test_suite("unit/libsais_mem"))
{
    bwa_tests::TempDir tmp("libsais_mem");
    // Odd lengths end the doubled text in a partial byte; even ones fill the last byte and
    // write_doubled_pac then adds a pad byte before the count byte.
    const int64_t sizes[] = {1, 3, 4, 5, 1003, 1004};
    for (int64_t l_pac : sizes) {
        std::vector<uint8_t> bases((size_t)l_pac);
        uint32_t st = 777u + (uint32_t)l_pac;
        for (int64_t l = 0; l < l_pac; ++l) {
            st = st * 1664525u + 1013904223u;
            bases[(size_t)l] = (uint8_t)((st >> 24) & 3);
        }
        std::vector<uint8_t> doubled(bases);
        for (int64_t l = l_pac - 1; l >= 0; --l) doubled.push_back((uint8_t)(3 - bases[(size_t)l]));
        const std::vector<uint8_t> fwd = pack_bases(bases);
        const std::vector<uint8_t> want = pack_bases(doubled);
        const size_t nbytes = (size_t)((2 * l_pac + 3) / 4);
        REQUIRE(want.size() == nbytes);

        const std::string in = tmp.path() + "/fwd.pac", out = tmp.path() + "/doubled.pac";
        FILE *f = std::fopen(in.c_str(), "wb");
        REQUIRE(f != nullptr);
        REQUIRE(std::fwrite(fwd.data(), 1, fwd.size(), f) == fwd.size());
        std::fclose(f);
        for (int threads : {1, 4}) {
            CAPTURE(l_pac);
            CAPTURE(threads);
            PackedText view(fwd.data(), l_pac);
            std::unique_ptr<uint8_t[]> mem = pack_doubled(view, l_pac, threads);
            REQUIRE(mem != nullptr);
            CHECK(std::vector<uint8_t>(mem.get(), mem.get() + nbytes) == want);

            write_doubled_pac(in.c_str(), l_pac, out.c_str(), threads);
            const std::vector<uint8_t> disk = read_file(out);
            const bool pad = (2 * l_pac) % 4 == 0;
            REQUIRE(disk.size() == nbytes + (pad ? 2 : 1));
            CHECK(std::vector<uint8_t>(disk.data(), disk.data() + nbytes) == want);
            if (pad) CHECK(disk[nbytes] == 0);
            CHECK(disk.back() == (uint8_t)((2 * l_pac) % 4));
        }
    }
}

TEST_CASE("libsais_mem: MemSink holds the same bytes the file writer writes"
          * doctest::test_suite("unit/libsais_mem"))
{
    bwa_tests::TempDir tmp("libsais_mem_sink");
    // N is the doubled-text length; the index holds N + 1 rows. 4099 rows fill no block or
    // sample stride exactly; 4096 rows do, so the layout's last cp_occ entry and last SA sample
    // are reserved but never written, and must read back zero as they do from the file.
    for (const int64_t N : {(int64_t)4098, (int64_t)4095}) {
        std::vector<uint8_t> buf((size_t)N + 1, 0);   // libsais alphabet: 1..4 = ACGT, buf[N] = 0
        uint32_t st = 12345u;
        int64_t count[5] = {0, 0, 0, 0, 0};
        int64_t freq[4] = {0, 0, 0, 0};
        for (int64_t i = 0; i < N; ++i) {
            st = st * 1664525u + 1013904223u;
            const uint8_t b = (uint8_t)((st >> 24) & 3);
            buf[(size_t)i] = (uint8_t)(b + 1);
            ++freq[b];
        }
        for (int c = 1; c <= 4; ++c) count[c] = count[c - 1] + freq[c - 1];

        std::vector<int32_t> sa((size_t)N + 1);
        for (int64_t i = 0; i <= N; ++i) sa[(size_t)i] = (int32_t)i;
        std::sort(sa.begin(), sa.end(), [&](int32_t a, int32_t b) {
            if (a == b) return false;   // irreflexive, as std::sort requires
            while (true) {   // the unique 0 terminator at buf[N] ends every other comparison
                if (buf[(size_t)a] != buf[(size_t)b]) return buf[(size_t)a] < buf[(size_t)b];
                ++a; ++b;
            }
        });

        for (const int sa_compx : {0, 3})
            for (int threads : {1, 4}) {
                CAPTURE(N);
                CAPTURE(sa_compx);
                CAPTURE(threads);
                const std::string path = tmp.path() + "/idx.bwt.2bit.64";
                int64_t sent_file = -1;
                write_fm_index_streaming(path.c_str(), buf.data(), sa.data(), false, N, count,
                                         &sent_file, sa_compx, /*num_threads=*/1);
                const std::vector<uint8_t> file = read_file(path);

                FmArrays arr;
                MemSink sink(N + 1, sa_compx, &arr);
                int64_t sent_mem = -1;
                write_fm_index_streaming(sink, buf.data(), sa.data(), false, N, count, &sent_mem,
                                         sa_compx, threads);

                const Bwt2bitLayout L = fmi_bwt2bit_layout(N + 1, sa_compx);
                REQUIRE(file.size() == (size_t)L.total_bytes);
                CHECK(sent_mem == sent_file);
                CHECK(arr.ref_seq_len == N + 1);
                CHECK(std::memcmp(arr.count, file.data() + 8, sizeof arr.count) == 0);
                CHECK(arr.sentinel_index == sent_file);
                CHECK(arr.sa_compx == sa_compx);
                CHECK(std::memcmp(arr.cp_occ, file.data() + L.off_cp_occ,
                                  (size_t)L.cp_occ_count * sizeof(CP_OCC)) == 0);
                CHECK(std::memcmp(arr.sa_ms_byte, file.data() + L.off_ms_byte,
                                  (size_t)L.sa_sample_count) == 0);
                CHECK(std::memcmp(arr.sa_ls_word, file.data() + L.off_ls_word,
                                  (size_t)L.sa_sample_count * sizeof(uint32_t)) == 0);
                // Independently derive every checkpoint and sample from the naive suffix array,
                // so a shared serial/parallel writer error cannot satisfy the file comparison.
                int64_t running[4] = {0, 0, 0, 0};
                for (int64_t row = 0; row <= N; ++row) {
                    const int64_t block = row / CP_BLOCK_SIZE;
                    if (row % CP_BLOCK_SIZE == 0)
                        for (int c = 0; c < 4; ++c) CHECK(arr.cp_occ[block].cp_count[c] == running[c]);
                    const int64_t pos = sa[(size_t)row];
                    const uint8_t b = buf[(size_t)(pos == 0 ? N : pos - 1)];
                    for (int c = 0; c < 4; ++c) {
                        const uint64_t bit = (arr.cp_occ[block].one_hot_bwt_str[c] >> (63 - row % CP_BLOCK_SIZE)) & 1;
                        CHECK(bit == (uint64_t)(b == c + 1));
                    }
                    if (b != 0) ++running[b - 1];
                    else CHECK(arr.sentinel_index == row);
                    if (row % ((int64_t)1 << sa_compx) == 0) {
                        const int64_t sample = row >> sa_compx;
                        CHECK((uint8_t)arr.sa_ms_byte[sample] == (uint8_t)((uint64_t)pos >> 32));
                        CHECK(arr.sa_ls_word[sample] == (uint32_t)pos);
                    }
                }
                if ((N + 1) % ((int64_t)1 << sa_compx) == 0) {
                    CHECK(arr.sa_ms_byte[L.sa_sample_count - 1] == 0);
                    CHECK(arr.sa_ls_word[L.sa_sample_count - 1] == 0);
                }
            }
    }
}

TEST_CASE("libsais_mem: memory build equals the index file"
          * doctest::test_suite("unit/libsais_mem"))
{
    const std::vector<std::vector<bwa_tests::Contig>> refs = {
        {{"one", random_bases(5000, 7u)}}, multi_contig_ref()};
    for (size_t r = 0; r < refs.size(); ++r) {
        bwa_tests::TempIndex idx(refs[r]);
        const std::string &prefix = idx.prefix();
        const int64_t l_pac = idx.bns()->l_pac;
        const std::vector<uint8_t> pac = read_file(prefix + ".pac");
        REQUIRE(pac.size() >= (size_t)((l_pac + 3) / 4));
        const std::vector<uint8_t> file = read_file(prefix + ".bwt.2bit.64");
        const int sa_compx = 3;
        const Bwt2bitLayout L = fmi_bwt2bit_layout(2 * l_pac + 1, sa_compx);
        REQUIRE(file.size() == (size_t)L.total_bytes);
        for (int threads : {1, 8}) {
            CAPTURE(r);
            CAPTURE(threads);
            LibsaisBuildOpts opts;
            opts.num_threads = threads;
            opts.sa_compx = sa_compx;
            FmArrays arr;
            REQUIRE(libsais_build_fm_memory(pac.data(), l_pac, opts, &arr) == 0);
            CHECK(arr.ref_seq_len == 2 * l_pac + 1);
            CHECK(arr.sa_compx == sa_compx);
            CHECK(std::memcmp(&arr.ref_seq_len, file.data(), 8) == 0);
            CHECK(std::memcmp(arr.count, file.data() + 8, sizeof arr.count) == 0);
            int64_t sent;
            std::memcpy(&sent, file.data() + L.off_sentinel, sizeof sent);
            CHECK(arr.sentinel_index == sent);
            CHECK(std::memcmp(arr.cp_occ, file.data() + L.off_cp_occ,
                              (size_t)L.cp_occ_count * sizeof(CP_OCC)) == 0);
            CHECK(std::memcmp(arr.sa_ms_byte, file.data() + L.off_ms_byte,
                              (size_t)L.sa_sample_count) == 0);
            CHECK(std::memcmp(arr.sa_ls_word, file.data() + L.off_ls_word,
                              (size_t)L.sa_sample_count * sizeof(uint32_t)) == 0);
        }
    }
}

TEST_CASE("libsais_mem: a memory budget refusal allocates nothing"
          * doctest::test_suite("unit/libsais_mem"))
{
    const int64_t l_pac = 1000;
    std::vector<uint8_t> pac((size_t)(l_pac + 3) / 4, 0x1b);
    LibsaisBuildOpts opts;
    opts.max_memory_bytes = 1;   // far below the estimate
    FmArrays arr;
    CHECK(libsais_build_fm_memory(pac.data(), l_pac, opts, &arr) == 3);
    CHECK(arr.cp_occ == nullptr);
    CHECK(arr.sa_ms_byte == nullptr);
    CHECK(arr.sa_ls_word == nullptr);
}

TEST_CASE("libsais_mem: memory budgets include resident output at every sampling rate"
          * doctest::test_suite("unit/libsais_mem"))
{
    const int64_t l_pac = 4096, N = 2 * l_pac;
    std::vector<uint8_t> pac((size_t)(l_pac + 3) / 4, 0x1b);
    for (int sa_compx = 0; sa_compx <= CP_SHIFT; ++sa_compx) {
        CAPTURE(sa_compx);
        LibsaisBuildOpts opts;
        opts.num_threads = 1;
        opts.sa_compx = sa_compx;
        const Bwt2bitLayout L = fmi_bwt2bit_layout(N + 1, sa_compx);
        const int64_t resident_bytes = (N + 3) / 4 + L.cp_occ_count * (int64_t)sizeof(CP_OCC)
                                     + L.sa_sample_count * (int64_t)(sizeof(int8_t) + sizeof(uint32_t));
        const int64_t sa_bytes = libsais_sa_is_64bit(N) ? sizeof(int64_t) : sizeof(int32_t);
        const int64_t need = libsais_estimate_peak_bytes(N) + resident_bytes + 10000 * sa_bytes;
        // Mandatory live allocations, derived directly rather than from the build estimator:
        // packed doubled text, unpacked text, full SA plus scratch, checkpoints and samples.
        const int64_t live_bytes = (N + 3) / 4 + (N + 1) + (N + 1 + 10000) * sa_bytes
                                 + (((N + 1) / 64) + 1) * (int64_t)sizeof(CP_OCC)
                                 + (((N + 1) / ((int64_t)1 << sa_compx)) + 1) * 5;
        REQUIRE(need >= live_bytes);
        opts.max_memory_bytes = libsais_estimate_peak_bytes(N) + resident_bytes;
        FmArrays without_scratch;
        CHECK(libsais_build_fm_memory(pac.data(), l_pac, opts, &without_scratch) == 3);
        CHECK(without_scratch.cp_occ == nullptr);
        CHECK(without_scratch.sa_ms_byte == nullptr);
        CHECK(without_scratch.sa_ls_word == nullptr);
        opts.max_memory_bytes = live_bytes - 1;
        FmArrays below_live;
        CHECK(libsais_build_fm_memory(pac.data(), l_pac, opts, &below_live) == 3);
        CHECK(below_live.cp_occ == nullptr);
        CHECK(below_live.sa_ms_byte == nullptr);
        CHECK(below_live.sa_ls_word == nullptr);
        if (sa_compx <= 1) {
            opts.max_memory_bytes = libsais_estimate_peak_bytes(N);
            FmArrays refused;
            CHECK(libsais_build_fm_memory(pac.data(), l_pac, opts, &refused) == 3);
            CHECK(refused.cp_occ == nullptr);
            CHECK(refused.sa_ms_byte == nullptr);
            CHECK(refused.sa_ls_word == nullptr);
        }
        opts.max_memory_bytes = need - 1;
        FmArrays refused;
        CHECK(libsais_build_fm_memory(pac.data(), l_pac, opts, &refused) == 3);
        CHECK(refused.cp_occ == nullptr);
        CHECK(refused.sa_ms_byte == nullptr);
        CHECK(refused.sa_ls_word == nullptr);
        opts.max_memory_bytes = need;
        FmArrays built;
        REQUIRE(libsais_build_fm_memory(pac.data(), l_pac, opts, &built) == 0);
        CHECK(built.ref_seq_len == N + 1);
        CHECK(built.sa_compx == sa_compx);
    }
}

TEST_CASE("libsais_mem: adopted index equals the loaded index"
          * doctest::test_suite("unit/libsais_mem"))
{
    bwa_tests::TempIndex idx(multi_contig_ref());   // built, then loaded with its pac
    const std::string &prefix = idx.prefix();
    FMI_search &loaded = *idx.fmi();

    bntseq_t *bns0 = bns_restore(prefix.c_str());
    REQUIRE(bns0 != nullptr);
    const int64_t l_pac = bns0->l_pac;
    const std::vector<uint8_t> pac_file = read_file(prefix + ".pac");
    REQUIRE(pac_file.size() >= (size_t)(l_pac / 4 + 1));

    LibsaisBuildOpts opts;
    opts.num_threads = 2;
    opts.sa_compx = 3;
    auto arrays = std::unique_ptr<FmArrays>(new FmArrays());
    REQUIRE(libsais_build_fm_memory(pac_file.data(), l_pac, opts, arrays.get()) == 0);

    // The same bytes bwa_idx_load_ele would have slurped (calloc'd l_pac/4+1).
    uint8_t *pac = (uint8_t *)calloc((size_t)(l_pac / 4 + 1), 1);
    REQUIRE(pac != nullptr);
    std::memcpy(pac, pac_file.data(), (size_t)(l_pac / 4 + 1));

    FMI_search adopted("(memory)");
    adopted.adopt_built(arrays.get(), bns0, pac, 2);
    CHECK(arrays->cp_occ == nullptr);
    CHECK(arrays->sa_ms_byte == nullptr);
    CHECK(arrays->sa_ls_word == nullptr);
    CHECK(adopted.idx->bns == bns0);
    CHECK(adopted.idx->pac == pac);
    CHECK(adopted.idx->mem == nullptr);

    CHECK(adopted.reference_seq_len == loaded.reference_seq_len);
    CHECK(adopted.sentinel_index == loaded.sentinel_index);
    CHECK(adopted.sa_compx == loaded.sa_compx);
    CHECK(adopted.sa_compx_mask == loaded.sa_compx_mask);
    CHECK(std::memcmp(adopted.count_data(), loaded.count_data(), 5 * sizeof(int64_t)) == 0);
    CHECK(std::memcmp(adopted.cp_occ_data(), loaded.cp_occ_data(),
                      (size_t)loaded.cp_occ_size_bytes()) == 0);
    CHECK(std::memcmp(adopted.sa_ms_byte_data(), loaded.sa_ms_byte_data(),
                      (size_t)loaded.sa_ms_byte_size_bytes()) == 0);
    CHECK(std::memcmp(adopted.sa_ls_word_data(), loaded.sa_ls_word_data(),
                      (size_t)loaded.sa_ls_word_size_bytes()) == 0);
    CHECK(adopted.idx->bns->l_pac == loaded.idx->bns->l_pac);
    CHECK(adopted.idx->bns->n_seqs == loaded.idx->bns->n_seqs);
    CHECK(std::memcmp(adopted.idx->pac, loaded.idx->pac, (size_t)(l_pac / 4 + 1)) == 0);

    bwa_aln::bwt_t *bl = bwa_aln::bwt_from_fmi(&loaded);
    bwa_aln::bwt_t *ba = bwa_aln::bwt_from_fmi(&adopted);
    CHECK(ba->primary == bl->primary);
    CHECK(ba->seq_len == bl->seq_len);
    for (int c = 0; c < 5; ++c) CHECK(ba->L2[c] == bl->L2[c]);
    uint32_t st = 99u;
    for (int i = 0; i < 100; ++i) {
        st = st * 1664525u + 1013904223u;
        const int64_t k = (int64_t)(st >> 8) % (2 * l_pac);
        for (int c = 0; c < 4; ++c)
            CHECK(bwa_aln::bwt_occ(ba, k, (uint8_t)c) == bwa_aln::bwt_occ(bl, k, (uint8_t)c));
        const int64_t row = (int64_t)(st >> 4) % loaded.reference_seq_len;
        CHECK(adopted.get_sa_entry_compressed(row) == loaded.get_sa_entry_compressed(row));
    }
    bwa_aln::bwt_destroy(bl);
    bwa_aln::bwt_destroy(ba);
}

TEST_CASE("libsais_mem: adopt_built accepts a null pac (BNS only)"
          * doctest::test_suite("unit/libsais_mem"))
{
    bwa_tests::TempIndex idx(multi_contig_ref());
    const std::string &prefix = idx.prefix();
    bntseq_t *bns = bns_restore(prefix.c_str());
    REQUIRE(bns != nullptr);
    const std::vector<uint8_t> pac_file = read_file(prefix + ".pac");
    LibsaisBuildOpts opts;
    opts.num_threads = 1;
    FmArrays arrays;
    REQUIRE(libsais_build_fm_memory(pac_file.data(), bns->l_pac, opts, &arrays) == 0);
    {
        FMI_search fmi("(memory)");
        fmi.adopt_built(&arrays, bns, nullptr, 1, false);
        CHECK(fmi.idx->pac == nullptr);
        CHECK(fmi.idx->bns == bns);
        CHECK(arrays.cp_occ == nullptr);
    }   // ~FMI_search / ~indexEle free the arrays and bns without error
}
