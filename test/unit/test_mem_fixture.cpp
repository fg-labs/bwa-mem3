// test/unit/test_mem_fixture.cpp
//
// Self-checks for test/framework/mem_fixture.h, which the `mem` path unit tests
// build on: the in-memory synthetic reference must be the one
// test/fixtures/mem_synth_ref.awk writes (same pinned POSIX cksum as
// test/framework/mem_cli.sh checks for the awk output), its planted features
// must sit where the awk script's header says, and TempIndex must build and
// load an index whose contig table matches what it was given. A fixture that
// drifted would move every expected coordinate derived from it.

#include "doctest/doctest.h"
#include "mem_fixture.h"

#include <string>
#include <vector>

using bwa_tests::Contig;

TEST_CASE("mem_fixture: synthetic reference matches the awk generator's pinned cksum"
          * doctest::test_suite("unit/util")) {
    const std::string fa = bwa_tests::mem_synth_fasta();
    const std::pair<uint32_t, uint64_t> sum = bwa_tests::posix_cksum(fa);
    CHECK(sum.first == bwa_tests::kMemSynthCksum);
    CHECK(sum.second == bwa_tests::kMemSynthBytes);
}

TEST_CASE("mem_fixture: posix_cksum matches POSIX cksum on known inputs"
          * doctest::test_suite("unit/util")) {
    // Values from `printf '' | cksum` and `printf 'ACGT\n' | cksum` (POSIX
    // cksum, identical under GNU coreutils and BSD).
    CHECK(bwa_tests::posix_cksum("").first == 4294967295u);
    CHECK(bwa_tests::posix_cksum("ACGT\n").first == 1973602648u);
}

TEST_CASE("mem_fixture: planted repeats sit where mem_synth_ref.awk documents them"
          * doctest::test_suite("unit/util")) {
    const std::vector<Contig> ref = bwa_tests::mem_synth_reference();
    REQUIRE(ref.size() == 4);
    CHECK(ref[0].name == "ctg1");
    CHECK(ref[0].seq.size() == 5000);
    CHECK(ref[1].seq.size() == 4000);
    CHECK(ref[2].seq.size() == 3000);
    CHECK(ref[3].seq.size() == 600);
    // 1-based inclusive coordinates, as in the awk header.
    auto at = [&](int c, int pos, int len) { return ref[(size_t)c].seq.substr((size_t)pos - 1, (size_t)len); };
    const std::string R = at(0, 2001, 300), S = at(0, 3001, 500);
    CHECK(at(1, 1001, 300) == R);
    CHECK(at(3, 101, 300) == bwa_tests::revcomp(R));
    const std::string Sd = at(1, 2001, 500);
    int diffs = 0;
    for (size_t i = 0; i < S.size(); ++i) {
        if (S[i] != Sd[i]) {
            ++diffs;
            CHECK((i + 1) % 20 == 10);
        }
    }
    CHECK(diffs == 25);
    CHECK(at(2, 1501, 40) == std::string(40, 'N'));
    CHECK(at(2, 2201, 50) == std::string(50, 'A'));
    CHECK(at(2, 501, 4) == "ACAC");
}

TEST_CASE("mem_fixture: TempIndex builds and loads the contigs it was given"
          * doctest::test_suite("unit/util")) {
    const std::vector<Contig> ref = bwa_tests::mem_synth_reference();
    bwa_tests::TempIndex idx(ref);
    const bntseq_t *bns = idx.bns();
    REQUIRE(bns != NULL);
    REQUIRE(bns->n_seqs == 4);
    int64_t total = 0;
    for (int i = 0; i < bns->n_seqs; ++i) {
        CHECK(std::string(bns->anns[i].name) == ref[(size_t)i].name);
        CHECK(bns->anns[i].len == (int64_t)ref[(size_t)i].seq.size());
        CHECK(bns->anns[i].offset == total);
        total += bns->anns[i].len;
    }
    CHECK(bns->l_pac == total);
    // The 40 N of ctg3 are its one ambiguity run.
    REQUIRE(bns->n_holes == 1);
    CHECK(bns->ambs[0].offset == 5000 + 4000 + 1500);
    CHECK(bns->ambs[0].len == 40);
    // pac holds the forward strand, two bits per base, first base in the high bits.
    REQUIRE(idx.pac() != NULL);
    const std::vector<uint8_t> want = bwa_tests::nt4(ref[0].seq.substr(0, 64));
    for (int i = 0; i < 64; ++i) {
        CAPTURE(i);
        CHECK(((idx.pac()[i >> 2] >> ((~i & 3) << 1)) & 3) == want[(size_t)i]);
    }
}
