// test/unit/test_mem_io_header.cpp
//
// The -R read-group line (bwa_set_rg), the index header sidecar
// (bwa_load_hdr_from_index) and the full SAM text header with -H/-R lines and
// @PG (bwa_print_sam_hdr2). Expected values follow the documented contracts:
//
//   - bwa_set_rg: the line must start with "@RG" and carry "\tID:"; "\t",
//     "\n", "\r" and "\\" escapes are expanded; the ID (up to the next TAB or
//     newline, at most 255 characters) becomes bwa_rg_id, the RG:Z value of
//     every record; a bad line returns NULL and clears the ID (bwa's model).
//   - bwa_load_hdr_from_index (src/bwa.cpp header comment): "<prefix>.hdr" if it
//     exists, else "<baseprefix>.dict", where baseprefix drops a trailing ".gz"
//     and then the basename's last dotted suffix; carriage returns and trailing
//     newlines removed; NULL when there is neither file or it is empty.
//   - bwa_print_sam_hdr2: @HD, the @SQ block, the sidecar's other lines, the
//     -H/-R lines, then @PG (test_sam_hdr_emit.cpp covers the @HD/@SQ
//     precedence rules; this checks the whole text with every source present).

#include "doctest/doctest.h"
#include "mem_fixture.h"

#include "bwa.h"

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>

// The @PG line main() builds (src/bwa.cpp defines it; bwa.h does not declare it).
extern char *bwa_pg;

namespace {

// Restores the process globals these functions read or write.
struct GlobalsGuard {
    char rg[sizeof(bwa_rg_id)];
    int verbose;
    char *pg;
    GlobalsGuard() : verbose(bwa_verbose), pg(bwa_pg) { memcpy(rg, bwa_rg_id, sizeof(rg)); }
    ~GlobalsGuard() {
        memcpy(bwa_rg_id, rg, sizeof(rg));
        bwa_verbose = verbose;
        bwa_pg = pg;
    }
};

// bwa_set_rg's returned line ("" for NULL) and the ID it left behind.
std::string set_rg(const std::string &arg, std::string *id) {
    bwa_tests::QuietStderr quiet;  // a refused line is reported on stderr
    char *line = bwa_set_rg(arg.c_str());
    std::string out = line ? line : "<null>";
    free(line);
    *id = bwa_rg_id;
    return out;
}

std::string load_sidecar(const std::string &prefix) {
    char *s = bwa_load_hdr_from_index(prefix.c_str());
    std::string out = s ? s : "<null>";
    free(s);
    return out;
}

}  // namespace

TEST_CASE("mem_io header: -R sets the read-group line and the RG:Z ID"
          * doctest::test_suite("unit/mem_io")) {
    GlobalsGuard guard;
    std::string id;

    SUBCASE("literal TABs") {
        CHECK(set_rg("@RG\tID:grp1\tSM:s1", &id) == "@RG\tID:grp1\tSM:s1");
        CHECK(id == "grp1");
    }
    SUBCASE("escaped TABs and backslashes, as typed on a command line") {
        CHECK(set_rg("@RG\\tID:g2\\tSM:a\\\\b\\tDS:x\\ny", &id) == "@RG\tID:g2\tSM:a\\b\tDS:x\ny");
        CHECK(id == "g2");
    }
    SUBCASE("ID is not the first field, and ends at the line's end") {
        CHECK(set_rg("@RG\\tSM:s\\tID:last", &id) == "@RG\tSM:s\tID:last");
        CHECK(id == "last");
    }
    SUBCASE("a 255-character ID is accepted, 256 is refused") {
        const std::string id255(255, 'x');
        CHECK(set_rg("@RG\tID:" + id255, &id) == "@RG\tID:" + id255);
        CHECK(id == id255);
        CHECK(set_rg("@RG\tID:" + id255 + "y", &id) == "<null>");
    }
    SUBCASE("refused lines clear the previous ID") {
        REQUIRE(set_rg("@RG\tID:keep", &id) == "@RG\tID:keep");
        CHECK(set_rg("RG\tID:x", &id) == "<null>");
        CHECK(id == "");
        REQUIRE(set_rg("@RG\tID:keep", &id) == "@RG\tID:keep");
        CHECK(set_rg("@RG\tSM:no_id", &id) == "<null>");
        CHECK(id == "");
        CHECK(set_rg("@CO\tID:x", &id) == "<null>");
    }
}

TEST_CASE("mem_io header: the index sidecar is <prefix>.hdr or else <baseprefix>.dict"
          * doctest::test_suite("unit/mem_io")) {
    bwa_tests::TempDir dir;
    const std::string d = dir.path() + "/";

    SUBCASE("neither file: NULL") { CHECK(load_sidecar(d + "ref.fa") == "<null>"); }
    SUBCASE(".hdr: CRs and trailing newlines removed") {
        dir.write("ref.fa.hdr", "@HD\tVN:1.6\r\n@SQ\tSN:a\tLN:5\r\n@CO\tc\n\n\n");
        CHECK(load_sidecar(d + "ref.fa") == "@HD\tVN:1.6\n@SQ\tSN:a\tLN:5\n@CO\tc");
    }
    SUBCASE(".dict when there is no .hdr: the basename's last suffix replaced") {
        dir.write("ref.dict", "@SQ\tSN:d\tLN:9\n");
        CHECK(load_sidecar(d + "ref.fa") == "@SQ\tSN:d\tLN:9");
        CHECK(load_sidecar(d + "ref.fasta.gz") == "@SQ\tSN:d\tLN:9");
        // A basename with no suffix keeps it whole: ref -> ref.dict.
        CHECK(load_sidecar(d + "ref") == "@SQ\tSN:d\tLN:9");
    }
    SUBCASE(".hdr wins over .dict") {
        dir.write("ref.fa.hdr", "@CO\tfrom hdr");
        dir.write("ref.dict", "@CO\tfrom dict");
        CHECK(load_sidecar(d + "ref.fa") == "@CO\tfrom hdr");
    }
    SUBCASE("only the last dotted suffix goes") {
        dir.write("GRCh38.dict", "@CO\tG");
        CHECK(load_sidecar(d + "GRCh38.p14") == "@CO\tG");
        CHECK(load_sidecar(d + "GRCh38.p14.fa") == "<null>");  // GRCh38.p14.dict
    }
    SUBCASE("an empty sidecar is none") {
        dir.write("ref.fa.hdr", "\n\n");
        CHECK(load_sidecar(d + "ref.fa") == "<null>");
    }
}

TEST_CASE("mem_io header: the SAM header with every source present"
          * doctest::test_suite("unit/mem_io")) {
    GlobalsGuard guard;
    bntann1_t anns[2];
    memset(anns, 0, sizeof(anns));
    anns[0].len = 1000;
    anns[0].name = const_cast<char *>("chrA");
    anns[0].anno = const_cast<char *>("");
    anns[1].offset = 1000;
    anns[1].len = 500;
    anns[1].name = const_cast<char *>("chrB");
    anns[1].anno = const_cast<char *>("");
    bntseq_t bns;
    memset(&bns, 0, sizeof(bns));
    bns.l_pac = 1500;
    bns.n_seqs = 2;
    bns.anns = anns;

    std::string pg = "@PG\tID:bwa-mem3\tPN:bwa-mem3\tVN:test\n";
    bwa_pg = const_cast<char *>(pg.c_str());
    // -H '@CO\tone' -H '@CO\ttwo' -R '@RG\tID:g\tSM:s', joined as main_mem does.
    char *hdr = bwa_insert_header("@CO\\tone", NULL);
    hdr = bwa_insert_header("@CO\ttwo", hdr);
    hdr = bwa_insert_header("@RG\tID:g\tSM:s", hdr);
    REQUIRE(hdr != NULL);

    bwa_tests::TempDir dir;
    auto emit = [&](const char *sidecar) -> std::string {
        FILE *fp = fopen(dir.file("h.sam").c_str(), "w");
        REQUIRE(fp != NULL);
        bwa_print_sam_hdr2(&bns, sidecar, hdr, fp, NULL);
        REQUIRE(fclose(fp) == 0);
        FILE *in = fopen(dir.file("h.sam").c_str(), "r");
        REQUIRE(in != NULL);
        std::string text;
        char buf[512];
        size_t n;
        while ((n = fread(buf, 1, sizeof(buf), in)) > 0) text.append(buf, n);
        fclose(in);
        return text;
    };
    CHECK(emit(NULL)
          == "@HD\tVN:1.5\tSO:unsorted\tGO:query\n@SQ\tSN:chrA\tLN:1000\n@SQ\tSN:chrB\tLN:500\n"
             "@CO\tone\n@CO\ttwo\n@RG\tID:g\tSM:s\n" + pg);
    // A sidecar with @SQ and another record: its @SQ replace the generated
    // block, its other records follow, before the user's lines.
    CHECK(emit("@SQ\tSN:chrA\tLN:1000\tM5:aa\n@SQ\tSN:chrB\tLN:500\n@CO\tsidecar")
          == "@HD\tVN:1.5\tSO:unsorted\tGO:query\n@SQ\tSN:chrA\tLN:1000\tM5:aa\n@SQ\tSN:chrB\tLN:500\n"
             "@CO\tsidecar\n@CO\tone\n@CO\ttwo\n@RG\tID:g\tSM:s\n" + pg);
    free(hdr);
}
