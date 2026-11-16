// htslib and bwa-mem3 kstring types must not cross APIs; keep BAM reading in its own TU.
#include "htslib/sam.h"
#include "htslib/kstring.h"
#include "doctest/doctest.h"
#include "aln_test_ref.h"
#include "temp_dir.h"
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
        CHECK(std::string(bam_aux2Z(xr)) == ((rec->core.flag & BAM_FREAD2) ? "GA" : "CT"));
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

TEST_CASE("aln_meth: paired trimmed mates and orphan records agree in SAM and BAM"
          * doctest::test_suite("unit/aln_meth_bam")) {
    bwa_tests::TempDir dir("aln_meth_pe_bam");
    const std::string fa=dir.file("ref.fa"),f1=dir.file("r1.fq"),f2=dir.file("r2.fq");
    uint32_t seed=918;const std::string genome=aln_test::lcg_dna(seed,4000);
    REQUIRE(aln_test::write_test_fasta(fa,{{"one",genome}}));
    REQUIRE(invoke({"index","--meth",fa},true)==0);
    FILE *fp[2]={fopen(f1.c_str(),"w"),fopen(f2.c_str(),"w")};
    REQUIRE(fp[0]!=nullptr);REQUIRE(fp[1]!=nullptr);
    for(int i=0;i<42;++i) {
        const int start=100+60*i,span=180+i%11;
        for(int j=0;j<2;++j) {
            std::string read=j?aln_test::revcomp(genome.substr(start+span-60,60)):genome.substr(start,60);
            if(i==40&&j==1) read=std::string(60,'N');
            if(i==41) read=std::string(60,'N');
            const int trim=i%2?17000:5;
            read+=std::string(trim,'A');
            const std::string qual=std::string(60,'I')+std::string(trim,'!');
            fprintf(fp[j],"@p%d\n%s\n+\n%s\n",i,read.c_str(),qual.c_str());
        }
    }
    REQUIRE(fclose(fp[0])==0);REQUIRE(fclose(fp[1])==0);
    const std::string sam=dir.file("out.sam"),bam=dir.file("out.bam");
    // The planted records, from the layout above: pair i's read 1 is forward at start+1 and read 2
    // reverse at start+span-59, each 60M with its trimmed tail soft-clipped (at the read's 3'
    // end); pair 40's read 2 and both of pair 41 are all N. A pair is proper only when its
    // insert is at least the longest untrimmed read, as sampe requires: so not the 17000-base
    // tails of the odd pairs.
    auto check_planted=[&](const std::vector<std::string> &recs) {
        REQUIRE(recs.size()==84);
        for(int i=0;i<42;++i) for(int j=0;j<2;++j) {
            CAPTURE(i);CAPTURE(j);
            std::vector<std::string> f;
            const std::string &line=recs[2*i+j];
            for(size_t at=0;;) {
                const size_t tab=line.find('\t',at);
                f.push_back(line.substr(at,tab==std::string::npos?tab:tab-at));
                if(tab==std::string::npos) break;
                at=tab+1;
            }
            REQUIRE(f.size()>=11);
            CHECK(f[0]=="p"+std::to_string(i));
            const int flag=std::stoi(f[1]),start=100+60*i,span=180+i%11,trim=i%2?17000:5;
            const int pos[2]={start+1,start+span-59};
            const std::string clip=std::to_string(trim)+"S";
            // SEQ and QUAL: the planted read (as written to the FASTQ, tail included), reverse
            // complemented and reversed on a reverse-strand record
            std::string seq=j?aln_test::revcomp(genome.substr(start+span-60,60)):genome.substr(start,60);
            if((i==40&&j==1)||i==41) seq=std::string(60,'N');
            seq+=std::string(trim,'A');
            std::string qual=std::string(60,'I')+std::string(trim,'!');
            if(flag&BAM_FREVERSE) {seq=aln_test::revcomp(seq);qual.assign(qual.rbegin(),qual.rend());}
            CHECK(f[9]==seq);CHECK(f[10]==qual);
            CHECK((flag&(BAM_FPAIRED|BAM_FREAD1|BAM_FREAD2))==(BAM_FPAIRED|(j?BAM_FREAD2:BAM_FREAD1)));
            if(i<40) {
                CHECK((flag&(BAM_FUNMAP|BAM_FMUNMAP))==0);
                CHECK(bool(flag&BAM_FPROPER_PAIR)==(i%2==0));
                CHECK(bool(flag&BAM_FREVERSE)==(j==1));CHECK(bool(flag&BAM_FMREVERSE)==(j==0));
                CHECK(f[2]=="one");CHECK(f[3]==std::to_string(pos[j]));CHECK(std::stoi(f[4])>0);
                CHECK(f[5]==(j?clip+"60M":"60M"+clip));
                CHECK(f[6]=="=");CHECK(f[7]==std::to_string(pos[1-j]));
                CHECK(f[8]==std::to_string(j?-span:span));
            } else if(i==40) {
                CHECK(f[2]=="one");CHECK(f[3]==std::to_string(pos[0]));CHECK(f[6]=="=");CHECK(f[7]==f[3]);
                CHECK(f[8]=="0");CHECK((flag&BAM_FPROPER_PAIR)==0);
                if(j==0) {
                    CHECK((flag&(BAM_FUNMAP|BAM_FMUNMAP|BAM_FREVERSE))==BAM_FMUNMAP);
                    CHECK(std::stoi(f[4])>0);CHECK(f[5]=="60M"+clip);
                } else {
                    CHECK((flag&(BAM_FUNMAP|BAM_FMUNMAP|BAM_FMREVERSE))==BAM_FUNMAP);
                    CHECK(f[4]=="0");CHECK(f[5]=="*");
                }
            } else {
                CHECK((flag&(BAM_FUNMAP|BAM_FMUNMAP|BAM_FPROPER_PAIR))==(BAM_FUNMAP|BAM_FMUNMAP));
                CHECK(f[2]=="*");CHECK(f[3]=="0");CHECK(f[5]=="*");CHECK(f[6]=="*");CHECK(f[7]=="0");
            }
        }
    };
    const std::vector<std::string> common={"aln","--meth","-q","20","-n","2","-l","1024","-s"};
    auto run=[&](bool binary,int threads) {
        auto args=common;args.insert(args.end(),{"-t",std::to_string(threads),"-f",binary?bam:sam});
        if(binary) args.push_back("--bam");
        args.insert(args.end(),{fa,f1,f2});return invoke(args);
    };
    REQUIRE(run(false,1)==0);const auto a=records(sam);check_planted(a);
    REQUIRE(run(true,3)==0);CHECK(records(bam)==a);
    REQUIRE(run(false,3)==0);CHECK(records(sam)==a);
    auto cpg=common;cpg.insert(cpg.end(),{"--meth=taps","--meth-context","cpg","-f",sam,fa,f1,f2});
    REQUIRE(invoke(cpg)==0);const auto context=records(sam);check_planted(context);
    cpg=common;cpg.insert(cpg.end(),{"--meth=taps","--meth-context","cpg","--bam","-f",bam,fa,f1,f2});
    REQUIRE(invoke(cpg)==0);CHECK(records(bam)==context);
}
