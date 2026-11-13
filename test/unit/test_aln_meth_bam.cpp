// htslib and bwa-mem3 kstring types must not cross APIs; keep BAM reading in its own TU.
#include "htslib/sam.h"
#include "htslib/kstring.h"
#include "doctest/doctest.h"
#include "aln_test_ref.h"
#include "../../src/aln/aln_main.h"
#include "../../src/main.h"
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>
#include <unistd.h>
#include <getopt.h>

namespace {
int invoke(std::vector<std::string> args, bool index = false) {
    std::vector<char *> argv;
    for (auto &a : args) argv.push_back(&a[0]);
    argv.push_back(nullptr);
    optind = 1;
#if defined(__APPLE__) || defined(__FreeBSD__) || defined(__OpenBSD__) || defined(__NetBSD__)
    optreset = 1;
#endif
    return index ? bwa_index((int)args.size(), argv.data()) : main_aln((int)args.size(), argv.data());
}
std::vector<std::string> records(const std::string &path) {
    samFile *fp = sam_open(path.c_str(), "r");
    REQUIRE(fp != nullptr);
    sam_hdr_t *hdr = sam_hdr_read(fp);
    REQUIRE(hdr != nullptr);
    bam1_t *rec = bam_init1();
    REQUIRE(rec != nullptr);
    std::vector<std::string> out;
    int rc;
    while ((rc = sam_read1(fp, hdr, rec)) >= 0) {
        uint8_t *xr = bam_aux_get(rec, "XR");
        REQUIRE(xr != nullptr);
        CHECK(xr[0] == 'Z');
        CHECK(std::string(bam_aux2Z(xr)) == "CT");
        if (!(rec->core.flag & BAM_FUNMAP)) {
            CHECK(bam_cigar2qlen(rec->core.n_cigar, bam_get_cigar(rec)) == rec->core.l_qseq);
            if (uint8_t *xa = bam_aux_get(rec, "XA")) {
                const std::string text = bam_aux2Z(xa);
                size_t start = 0;
                while (start < text.size()) {
                    size_t first = text.find(',', start), second = text.find(',', first + 1);
                    size_t end = text.find(',', second + 1);
                    REQUIRE(end != std::string::npos);
                    int span = 0, run = 0;
                    for (size_t k = second + 1; k < end; ++k) {
                        const char c = text[k];
                        if (c >= '0' && c <= '9') run = run * 10 + c - '0';
                        else { if (c == 'M' || c == 'I' || c == 'S') span += run; run = 0; }
                    }
                    CHECK(span == rec->core.l_qseq);
                    const size_t terminator = text.find(';', end);
                    REQUIRE(terminator != std::string::npos);
                    start = terminator + 1;
                }
            }
            for (const char *tag : {"XM", "XG", "MD"}) {
                uint8_t *value = bam_aux_get(rec, tag);
                REQUIRE(value != nullptr);
                CHECK(value[0] == 'Z');
            }
            REQUIRE(bam_aux_get(rec, "NM") != nullptr);
        } else {
            CHECK(bam_aux_get(rec, "XM") == nullptr);
            CHECK(bam_aux_get(rec, "XG") == nullptr);
            CHECK((rec->core.flag & BAM_FQCFAIL) == 0);
        }
        kstring_t text = {0, 0, nullptr};
        REQUIRE(sam_format1(hdr, rec, &text) >= 0);
        out.emplace_back(text.s, text.l);
        free(text.s);
    }
    CHECK(rc == -1);
    bam_destroy1(rec);
    sam_hdr_destroy(hdr);
    CHECK(sam_close(fp) == 0);
    return out;
}

// The tab-separated SAM fields of a record.
std::vector<std::string> fields(const std::string &record) {
    std::vector<std::string> f;
    size_t start = 0;
    for (size_t tab; (tab = record.find('\t', start)) != std::string::npos; start = tab + 1)
        f.push_back(record.substr(start, tab - start));
    f.push_back(record.substr(start));
    return f;
}
}

TEST_CASE("aln_meth: trimmed primary and XA CIGARs span the full query in SAM and BAM"
          * doctest::test_suite("unit/aln_meth_bam")) {
    bwa_tests::TempDir tmp{"aln_meth_trim"};
    const std::string &path = tmp.path();
    const std::string fa=path+"/ref.fa", fq=path+"/reads.fq";
    uint32_t seed=713;
    const std::string genome=aln_test::lcg_dna(seed,1000);
    REQUIRE(aln_test::write_test_fasta(fa,{{"one",genome},{"two",genome}}));
    REQUIRE(invoke({"index","--meth",fa},true)==0);
    for(int trim : {5,17000}) {
        FILE *fp=fopen(fq.c_str(),"w"); REQUIRE(fp!=nullptr);
        for(int record=0;record<80;++record) {
            const int reverse=record&1;
            std::string read=genome.substr(200,60);
            if(reverse) read=aln_test::revcomp(read);
            read+=std::string(trim,'A');
            const std::string quality=std::string(60,'I')+std::string(trim,'!');
            fprintf(fp,"@r%d\n%s\n+\n%s\n",record,read.c_str(),quality.c_str());
        }
        REQUIRE(fclose(fp)==0);
        const std::string sam=path+"/out.sam",bam=path+"/out.bam";
        REQUIRE(invoke({"aln","--meth","-q","20","-n","2","-l","1024","-f",sam,fa,fq})==0);
        const auto a=records(sam);
        REQUIRE(a.size()==80);
        // Independent expectation: both reads sit at position 201 of either identical contig, and
        // the trimmed tail is the read's 3' end, so it is clipped after the match on the forward
        // strand and before it on the reverse strand (the record shows the reverse complement).
        const std::string m="60M",clip=std::to_string(trim)+"S";
        for(int i=0;i<80;++i) {
            CAPTURE(i);
            const int reverse=i&1;
            const auto f=fields(a[i]);
            REQUIRE(f.size()>11);
            CHECK(f[1]==(reverse ? "16" : "0"));
            CHECK((f[2]=="one" || f[2]=="two"));
            CHECK(f[3]=="201");
            const std::string cigar=reverse ? clip+m : m+clip;
            CHECK(f[5]==cigar);
            const std::string other=f[2]=="one" ? "two" : "one";
            CHECK(a[i].find(std::string("\tXA:Z:")+other+","+(reverse ? "-" : "+")+"201,"+cigar+",0;")
                  !=std::string::npos);
        }
        REQUIRE(invoke({"aln","--meth","-q","20","-n","2","-l","1024","--bam","-f",bam,fa,fq})==0);
        CHECK(records(bam)==a);
        // Several worker slices repeatedly read the shared stock MAPQ table.
        REQUIRE(invoke({"aln","--meth","-q","20","-n","2","-l","1024","-t","4","-f",sam,fa,fq})==0);
        CHECK(records(sam)==a);
    }
}

TEST_CASE("aln_meth: SAM and BAM agree for both strands unmapped and QC-fail records"
          * doctest::test_suite("unit/aln_meth_bam")) {
    bwa_tests::TempDir tmp{"aln_meth_bam"};
    const std::string &path = tmp.path();
    const std::string fa = path + "/ref.fa", fq = path + "/r.fq";
    uint32_t seed = 1234;
    const std::string genome = aln_test::lcg_dna(seed, 1000);
    REQUIRE(aln_test::write_test_fasta(fa, {{"chr1", genome}}));
    REQUIRE(invoke({"index", "--meth", fa}, true) == 0);
    std::vector<std::string> reads = {genome.substr(200, 50), aln_test::revcomp(genome.substr(400, 50)), std::string(50, 'N')};
    FILE *fp = fopen(fq.c_str(), "w");
    REQUIRE(fp != nullptr);
    for (size_t i = 0; i < reads.size(); ++i) {
        for (char &c : reads[i]) if (c == 'C') c = 'T';
        fprintf(fp, "@r%zu\n%s\n+\n%s\n", i, reads[i].c_str(), std::string(50, 'I').c_str());
    }
    REQUIRE(fclose(fp) == 0);
    for (const char *failed : {"f", "r"}) {
        const std::string sam = path + "/out.sam", bam = path + "/out.bam";
        REQUIRE(invoke({"aln", "--meth", "--set-as-failed", failed, "-f", sam, fa, fq}) == 0);
        REQUIRE(invoke({"aln", "--meth", "--set-as-failed", failed, "--bam", "-f", bam, fa, fq}) == 0);
        const auto a = records(sam), b = records(bam);
        REQUIRE(a.size() == 3);
        REQUIRE(b.size() == a.size());
        CHECK(a == b);
        // SAM flags are parsed independently of htslib's formatting above.
        for (int i = 0; i < 2; ++i) {
            const size_t tab = a[i].find('\t');
            REQUIRE(tab != std::string::npos);
            const int flag = std::stoi(a[i].substr(tab + 1));
            CHECK((flag & BAM_FUNMAP) == 0);
            CHECK(bool(flag & BAM_FQCFAIL) == (failed[0] == (i == 0 ? 'f' : 'r')));
            CHECK(bool(flag & BAM_FREVERSE) == (i == 1));
        }
        // the all-N read is unmapped under either setting
        const size_t tab = a[2].find('\t');
        REQUIRE(tab != std::string::npos);
        CHECK(a[2].substr(0, tab) == "r2");
        CHECK((std::stoi(a[2].substr(tab + 1)) & BAM_FUNMAP) != 0);
    }
}
