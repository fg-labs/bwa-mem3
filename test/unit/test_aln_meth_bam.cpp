// htslib and bwa-mem3 kstring types must not cross APIs; keep BAM reading in its own TU.
#include "htslib/sam.h"
#include "htslib/kstring.h"
#include "doctest/doctest.h"
#include "aln_test_ref.h"
#include "temp_dir.h"
#include "../../src/aln/aln_main.h"
#include "../../src/main.h"
#include <cstdio>
#include <algorithm>
#include <cstdlib>
#include <set>
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

TEST_SUITE("unit/aln_meth_bam") {

TEST_CASE("aln_meth: neutral conversion guards own opposite-strand CpG CIGARs in SAM and BAM") {
    bwa_tests::TempDir tmp{"aln_meth_conv_bam"};const std::string &path=tmp.path();
    const std::string fa=path+"/ref.fa",fq=path+"/reads.fq",sam=path+"/out.sam",bam=path+"/out.bam";
    uint32_t seed=82731;
    std::string read=aln_test::lcg_dna(seed,50);
    read[20]='T';read[21]='G';
    std::string converted=read;converted[20]='C';
    const std::string spacer=aln_test::lcg_dna(seed,100);
    const std::string genome=spacer+read+spacer+converted+spacer;
    bool best_tie=false;
    SUBCASE("unique best with opposite-strand worse locus") {}
    SUBCASE("two best on opposite strands") {best_tie=true;}
    const std::string reverse_genome=best_tie?genome:spacer+converted+spacer;
    REQUIRE(aln_test::write_test_fasta(fa,{{"forward",genome},{"reverse",aln_test::revcomp(reverse_genome)}}));
    REQUIRE(invoke({"index","--meth",fa},true)==0);
    FILE *fp=fopen(fq.c_str(),"w");REQUIRE(fp!=nullptr);
    fprintf(fp,"@r0\n%s\n+\n%s\n",read.c_str(),std::string(50,'I').c_str());
    REQUIRE(fclose(fp)==0);
    // Opposite-strand alternatives must carry their own normalized alignment and tags.
    for(const char *context : {"all","cpg"}) {
        std::vector<std::string> args={"aln","--meth=5base","--meth-context",context,"-n","2","-l","1024","-f",sam,fa,fq};
        REQUIRE(invoke(args)==0);const auto a=records(sam);REQUIRE(a.size()==1);
        CHECK(a[0].find(best_tie?"\t0\t50M\t":"\t6\t50M\t")!=std::string::npos);
        CHECK(a[0].find(best_tie?"X0:i:2":"X0:i:1")!=std::string::npos);
        CHECK(a[0].find("X1:i:2")!=std::string::npos);
        CHECK(a[0].find("NM:i:0")!=std::string::npos);
        CHECK(a[0].find("MD:Z:50")!=std::string::npos);
        // Every locus but the primary is in XA with its own strand, position and ungapped CIGAR;
        // the converted copies need one free conversion, so their NM is 0 too.
        std::set<std::string> loci={"forward,+101,50M,0","forward,+251,50M,0","reverse,-101,50M,0"};
        if(best_tie) loci.insert("reverse,-251,50M,0");
        const auto f=fields(a[0]);
        const std::string primary=f[2]+","+((std::stoi(f[1])&16)?"-":"+")+f[3]+",50M,0";
        if(best_tie) CHECK((primary=="forward,+101,50M,0" || primary=="reverse,-251,50M,0"));
        else CHECK(primary=="forward,+101,50M,0");
        REQUIRE(loci.erase(primary)==1);
        const size_t xa=a[0].find("XA:Z:");REQUIRE(xa!=std::string::npos);
        const std::string list=a[0].substr(xa+5,a[0].find('\t',xa)-xa-5);
        std::set<std::string> listed;
        for(size_t at=0,end;(end=list.find(';',at))!=std::string::npos;at=end+1) listed.insert(list.substr(at,end-at));
        CHECK(listed==loci);
        args={"aln","--meth=5base","--meth-context",context,"-n","2","-l","1024","--bam","-f",bam,fa,fq};
        REQUIRE(invoke(args)==0);CHECK(records(bam)==a);
    }
}

TEST_CASE("aln_meth: all-resolved pairing beyond 500 preserves trimmed SAM and BAM ownership") {
    bwa_tests::TempDir tmp{"aln_meth_pair601"};const std::string &path=tmp.path();
    const std::string fa=path+"/ref.fa",one=path+"/one.fq",two=path+"/two.fq";
    uint32_t seed=172913;
    const std::string read=aln_test::lcg_dna(seed,36),mate=aln_test::lcg_dna(seed,36);
    std::string repeats=aln_test::lcg_dna(seed,50);
    for(int k=0;k<600;++k) repeats+=read+aln_test::lcg_dna(seed,64);
    const std::string pair=aln_test::lcg_dna(seed,100)+read+aln_test::lcg_dna(seed,150)+mate+aln_test::lcg_dna(seed,100);
    const std::string decoy=aln_test::lcg_dna(seed,100)+mate+aln_test::lcg_dna(seed,100);
    REQUIRE(aln_test::write_test_fasta(fa,{{"repeat",repeats},{"pair",pair},{"decoy",decoy}}));
    REQUIRE(invoke({"index","--meth",fa},true)==0);
    for(int j=0;j<2;++j) {
        FILE *fp=fopen((j?two:one).c_str(),"w");REQUIRE(fp!=nullptr);
        const std::string seq=(j?aln_test::revcomp(mate):read)+"AAAAA";
        const std::string qual=std::string(36,'I')+"!!!!!";
        fprintf(fp,"@r0\n%s\n+\n%s\n",seq.c_str(),qual.c_str());REQUIRE(fclose(fp)==0);
    }
    for(const char *chem:{"--meth","--meth=5base"}) {
        const std::string sam=path+"/out.sam",bam=path+"/out.bam";
        const std::vector<std::string> options={chem,"-q","20","-n","0","-o","0","-l","1024","--max-occ","601"};
        auto args=std::vector<std::string>{"aln"};args.insert(args.end(),options.begin(),options.end());
        args.insert(args.end(),{"-f",sam,fa,one,two});REQUIRE(invoke(args)==0);
        const auto a=records(sam);REQUIRE(a.size()==2);
        CHECK(a[0].find("r0\t99\tpair\t101\t") == 0);
        CHECK(a[1].find("r0\t147\tpair\t287\t") == 0);
        CHECK(a[0].find("\t36M5S\t")!=std::string::npos);
        CHECK(a[1].find("\t5S36M\t")!=std::string::npos);
        CHECK(a[0].find("X0:i:601")!=std::string::npos);
        args={"aln"};args.insert(args.end(),options.begin(),options.end());
        args.insert(args.end(),{"--bam","-f",bam,fa,one,two});REQUIRE(invoke(args)==0);
        CHECK(records(bam)==a);
    }
}

TEST_CASE("aln_meth: cap-bound rescue restores 17000-base tails in paired SAM and BAM") {
    bwa_tests::TempDir tmp{"aln_meth_cap_bam"};const std::string &path=tmp.path();
    const std::string fa=path+"/ref.fa",one=path+"/one.fq",two=path+"/two.fq";
    uint32_t seed=19337;const std::string genome=aln_test::lcg_dna(seed,5000);
    std::string target=aln_test::lcg_dna(seed,50),anchor=aln_test::lcg_dna(seed,50);
    target[20]='C';target[21]='G';
    std::string repeats=aln_test::lcg_dna(seed,50);
    for(int k=0;k<4800;++k) repeats+=target+aln_test::lcg_dna(seed,30);
    const std::string pair=aln_test::lcg_dna(seed,100)+target+aln_test::lcg_dna(seed,150)+anchor+aln_test::lcg_dna(seed,100);
    const std::string second=aln_test::lcg_dna(seed,100)+anchor+aln_test::lcg_dna(seed,100);
    for(bool reversed_reference:{false,true}) {
        auto contigs=std::vector<std::pair<std::string,std::string>>{{"training",genome},{"repeat",repeats},{"pair",pair},{"anchor",second}};
        if(reversed_reference) std::reverse(contigs.begin(),contigs.end());
        if(reversed_reference) REQUIRE(unlink((fa+".meth.fa").c_str())==0);
        REQUIRE(aln_test::write_test_fasta(fa,contigs));REQUIRE(invoke({"index","--meth",fa},true)==0);
        for(const char *chem:{"--meth","--meth=5base"}) for(bool reverse_fragment:{false,true}) {
            FILE *f1=fopen(one.c_str(),"w"),*f2=fopen(two.c_str(),"w");REQUIRE(f1!=nullptr);REQUIRE(f2!=nullptr);
            for(int k=0;k<41;++k) {
                const int start=100+70*k,span=230+(k*7)%31;
                std::string a=k==40?target:genome.substr(start,50);
                std::string b=k==40?aln_test::revcomp(anchor):aln_test::revcomp(genome.substr(start+span-50,50));
                if(reverse_fragment) std::swap(a,b);
                const int trim=k==40?17000:0;
                fprintf(f1,"@r%d\n%s%s\n+\n%s%s\n",k,a.c_str(),std::string(trim,'A').c_str(),std::string(50,'I').c_str(),std::string(trim,'!').c_str());
                fprintf(f2,"@r%d\n%s%s\n+\n%s%s\n",k,b.c_str(),std::string(trim,'A').c_str(),std::string(50,'I').c_str(),std::string(trim,'!').c_str());
            }
            REQUIRE(fclose(f1)==0);REQUIRE(fclose(f2)==0);
            const std::string sam=path+"/out.sam",bam=path+"/out.bam";
            std::vector<std::string> opts={chem,"-q","20","-n","0","-o","0","-l","1024"};
            if(std::string(chem)=="--meth=5base") opts.insert(opts.end(),{"--meth-context","cpg"});
            auto args=std::vector<std::string>{"aln"};args.insert(args.end(),opts.begin(),opts.end());
            args.insert(args.end(),{"-f",sam,fa,one,two});REQUIRE(invoke(args)==0);const auto a=records(sam);REQUIRE(a.size()==82);
            const int rescued=reverse_fragment?81:80,anchored=reverse_fragment?80:81;
            // The rescued end aligns `target` forward with its tail restored as a clip; the
            // anchor aligns reverse, so its record carries the reverse-complemented tail first.
            const auto rescued_rec=fields(a[rescued]),anchored_rec=fields(a[anchored]);
            REQUIRE(rescued_rec.size()>=11);REQUIRE(anchored_rec.size()>=11);
            CHECK(rescued_rec[3]=="101");CHECK(rescued_rec[5]=="50M17000S");
            CHECK(rescued_rec[9]==target+std::string(17000,'A'));
            CHECK(rescued_rec[10]==std::string(50,'I')+std::string(17000,'!'));
            CHECK(anchored_rec[3]=="301");CHECK(anchored_rec[5]=="17000S50M");
            CHECK((std::stoi(anchored_rec[1])&16)!=0);
            CHECK(anchored_rec[9]==std::string(17000,'T')+anchor);
            CHECK(anchored_rec[10]==std::string(17000,'!')+std::string(50,'I'));
            CHECK(a[rescued].find("XT:A:M")!=std::string::npos);
            CHECK(a[rescued].find("\tpair\t101\t0\t")!=std::string::npos);
            CHECK(a[rescued].find("X0:i:")==std::string::npos);CHECK(a[rescued].find("X1:i:")==std::string::npos);
            CHECK(a[rescued].find("XA:Z:")==std::string::npos);
            CHECK(a[anchored].find("XA:Z:")!=std::string::npos);
            const int expected_flag=reverse_fragment?129:65;
            CHECK(a[rescued].find("r40\t"+std::to_string(expected_flag+32)+"\t")==0);
            args={"aln"};args.insert(args.end(),opts.begin(),opts.end());args.insert(args.end(),{"--bam","-f",bam,fa,one,two});
            REQUIRE(invoke(args)==0);CHECK(records(bam)==a);
        }
    }
}

}
