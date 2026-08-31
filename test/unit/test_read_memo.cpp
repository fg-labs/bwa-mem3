// Unit tests for read_memo (--dedup-reads whole-read-pair memoization).
//
// Covers the two pieces of new logic that carry correctness weight in Phase 1:
//   (1) the pre-pass fingerprint + full-byte VERIFY grouping -- in particular
//       that two pairs which collide on the ends-only fingerprint but differ in
//       the middle are NOT grouped (the verify memcmp must catch them), that
//       both mates must match, and that mate ORDER matters; and
//   (2) --dedup-reads / BWAMEM3_DEDUP_READS mode resolution precedence.
//
// Phase 1 does not consume role[]/rep_pair[] in the aligner, so these tests pin
// the pre-pass semantics directly rather than through end-to-end output.

#include "doctest/doctest.h"

#include "read_memo.h"      // read_memo_prepass, read_memo_state, mode API
#include "bwa.h"            // bseq1_t

#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

namespace {

// A set of read-pairs whose bases stay alive for the pre-pass. seqs[] points
// into the owned strings; role/rep_pair come back from the module.
struct PairSet {
    std::vector<std::string> bases;   // 2 per pair (interleaved r1,r2)
    std::vector<bseq1_t>     seqs;

    void add(const std::string &r1, const std::string &r2) {
        bases.push_back(r1);
        bases.push_back(r2);
    }
    // Finalize: build the bseq1_t array pointing at the owned strings.
    read_memo_result run(read_memo_state *st) {
        seqs.assign(bases.size(), bseq1_t{});
        for (size_t i = 0; i < bases.size(); ++i) {
            memset(&seqs[i], 0, sizeof(bseq1_t));
            seqs[i].seq   = const_cast<char *>(bases[i].c_str());
            seqs[i].l_seq = (int)bases[i].size();
        }
        return read_memo_prepass(nullptr, seqs.data(), (int)seqs.size(), st);
    }
};

// A 150bp read whose 32-byte ends are all `end` and whose middle is `mid`.
// Two reads sharing `end` collide on the ends-only fingerprint regardless of
// their middles -- exactly the case the VERIFY memcmp must separate.
std::string longread(char end, char mid) {
    return std::string(40, end) + std::string(70, mid) + std::string(40, end);
}

read_memo_state fresh_state() { return read_memo_state{ 0, nullptr, nullptr, 0 }; }

} // namespace

TEST_CASE("identical pairs are grouped; distinct pairs are not") {
    PairSet ps;
    ps.add("AAAACCCC", "GGGGTTTT");   // pair 0
    ps.add("AAAACCCC", "GGGGTTTT");   // pair 1 == pair 0
    ps.add("AAAACCCC", "GGGGTTTA");   // pair 2: r2 differs -> distinct
    read_memo_state st = fresh_state();
    read_memo_result r = ps.run(&st);

    CHECK(r.pairs == 3);
    CHECK(r.dup_pairs == 1);
    CHECK(st.role[0] == 1);            // REP
    CHECK(st.role[1] == 2);            // DUP
    CHECK(st.rep_pair[1] == 0);        // of pair 0
    CHECK(st.role[2] == 1);            // REP (distinct)
}

TEST_CASE("ends-only fingerprint collision is caught by the byte verify") {
    // Both pairs share identical 32-byte ends on both mates but differ only in
    // the middle -> same fingerprint, must NOT be grouped.
    PairSet ps;
    ps.add(longread('A', 'C'), longread('G', 'T'));
    ps.add(longread('A', 'X'), longread('G', 'T'));   // r1 middle differs
    read_memo_state st = fresh_state();
    read_memo_result r = ps.run(&st);

    CHECK(r.dup_pairs == 0);
    CHECK(st.role[0] == 1);
    CHECK(st.role[1] == 1);            // distinct despite fingerprint collision
}

TEST_CASE("both mates must match, and mate order matters") {
    PairSet ps;
    ps.add("AAAA", "CCCC");            // pair 0
    ps.add("AAAA", "GGGG");            // r2 differs -> distinct
    ps.add("CCCC", "AAAA");            // swapped order -> distinct from pair 0
    ps.add("AAAA", "CCCC");            // == pair 0
    read_memo_state st = fresh_state();
    read_memo_result r = ps.run(&st);

    CHECK(r.dup_pairs == 1);
    CHECK(st.role[1] == 1);
    CHECK(st.role[2] == 1);
    CHECK(st.role[3] == 2);
    CHECK(st.rep_pair[3] == 0);
}

TEST_CASE("length mismatch is not a duplicate") {
    PairSet ps;
    ps.add("AAAACCCC", "GGGG");
    ps.add("AAAACCCC", "GGGGG");       // r2 one base longer
    read_memo_state st = fresh_state();
    read_memo_result r = ps.run(&st);
    CHECK(r.dup_pairs == 0);
}

TEST_CASE("three-way duplicate set collapses to one representative") {
    PairSet ps;
    for (int i = 0; i < 3; ++i) ps.add("ACGTACGT", "TTGGCCAA");
    read_memo_state st = fresh_state();
    read_memo_result r = ps.run(&st);
    CHECK(r.dup_pairs == 2);
    CHECK(st.role[0] == 1);
    CHECK(st.role[1] == 2);
    CHECK(st.role[2] == 2);
    CHECK(st.rep_pair[1] == 0);
    CHECK(st.rep_pair[2] == 0);
}

TEST_CASE("mode resolution: CLI value wins") {
    mem_dedup_reads_configure("off");  CHECK(read_memo_mode() == READMEMO_OFF);
    mem_dedup_reads_configure("on");   CHECK(read_memo_mode() == READMEMO_ON);
    mem_dedup_reads_configure("auto"); CHECK(read_memo_mode() == READMEMO_AUTO);
    mem_dedup_reads_configure("2");    CHECK(read_memo_mode() == READMEMO_AUTO);
    mem_dedup_reads_configure("0");    CHECK(read_memo_mode() == READMEMO_OFF);
}

TEST_CASE("mode resolution: env used when no CLI value, CLI overrides env") {
    setenv("BWAMEM3_DEDUP_READS", "on", 1);
    mem_dedup_reads_configure(nullptr);  CHECK(read_memo_mode() == READMEMO_ON);   // env
    mem_dedup_reads_configure("off");    CHECK(read_memo_mode() == READMEMO_OFF);  // CLI beats env
    unsetenv("BWAMEM3_DEDUP_READS");
}

TEST_CASE("controller latches auto ON when memoizing clearly wins") {
    read_memo_reset_for_testing();      // independent of prior cases / ordering
    mem_dedup_reads_configure("auto");
    REQUIRE(read_memo_mode() == READMEMO_AUTO);
    // Half the pairs are duplicates and the align phase is expensive, so the
    // predicted benefit (dup_pairs * align_ns/pairs) dwarfs the probe cost ->
    // net strongly negative. Vary align_ns so the z-test has nonzero variance.
    const uint64_t align_ns[4] = { 1000000, 1100000, 900000, 1050000 };
    for (int i = 0; i < 4; ++i) {
        read_memo_result r{ /*pairs*/ 1000, /*dup_pairs*/ 500, /*probe_ns*/ 1000 };
        read_memo_controller_observe(r, align_ns[i]);
    }
    CHECK(read_memo_active() == READMEMO_ON);   // auto + latched-on
    mem_dedup_reads_configure("off");
    CHECK(read_memo_active() == READMEMO_OFF);  // explicit off overrides the latch
}
