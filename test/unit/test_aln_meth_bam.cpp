// htslib and bwa-mem3 kstring types must not cross APIs; keep BAM reading in its own TU.
#include "htslib/sam.h"
#include "htslib/kstring.h"
#include "doctest/doctest.h"
#include "aln_test_ref.h"
#include "../../src/aln/aln_main.h"
#include "../../src/main.h"
#include <cstdio>
#include <algorithm>
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
}

TEST_CASE("aln_meth: trimmed primary and XA CIGARs span the full query in SAM and BAM") {
    std::string path = aln_test::tmp_template("aln_meth_trim");
    REQUIRE(mkdtemp(&path[0]) != nullptr);
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
        for(int i=0;i<80;++i) {
            CHECK(a[i].find("XA:Z:")!=std::string::npos);
            CHECK(a[i].find((i&1) ? "\t16\t" : "\t0\t")!=std::string::npos);
        }
        REQUIRE(invoke({"aln","--meth","-q","20","-n","2","-l","1024","--bam","-f",bam,fa,fq})==0);
        CHECK(records(bam)==a);
        // Several worker slices repeatedly read the shared stock MAPQ table.
        REQUIRE(invoke({"aln","--meth","-q","20","-n","2","-l","1024","-t","4","-f",sam,fa,fq})==0);
        CHECK(records(sam)==a);
    }
    CHECK(system(("rm -rf '"+path+"'").c_str())==0);
}

TEST_CASE("aln_meth: SAM and BAM agree for both strands unmapped and QC-fail records") {
    std::string path = aln_test::tmp_template("aln_meth_bam");
    REQUIRE(mkdtemp(&path[0]) != nullptr);
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
    }
    const std::string cmd = "rm -rf '" + path + "'";
    CHECK(system(cmd.c_str()) == 0);
}

TEST_CASE("aln_meth: paired trimmed mates and orphan records agree in SAM and BAM") {
    std::string path=aln_test::tmp_template("aln_meth_pe_bam");
    REQUIRE(mkdtemp(&path[0])!=nullptr);
    const std::string fa=path+"/ref.fa",f1=path+"/r1.fq",f2=path+"/r2.fq";
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
    const std::string sam=path+"/out.sam",bam=path+"/out.bam";
    const std::vector<std::string> common={"aln","--meth","-q","20","-n","2","-l","1024","-s"};
    auto run=[&](bool binary,int threads) {
        auto args=common;args.insert(args.end(),{"-t",std::to_string(threads),"-f",binary?bam:sam});
        if(binary) args.push_back("--bam");
        args.insert(args.end(),{fa,f1,f2});return invoke(args);
    };
    REQUIRE(run(false,1)==0);const auto a=records(sam);REQUIRE(a.size()==84);
    REQUIRE(run(true,3)==0);CHECK(records(bam)==a);
    REQUIRE(run(false,3)==0);CHECK(records(sam)==a);
    auto cpg=common;cpg.insert(cpg.end(),{"--meth=taps","--meth-context","cpg","-f",sam,fa,f1,f2});
    REQUIRE(invoke(cpg)==0);const auto context=records(sam);REQUIRE(context.size()==84);
    cpg=common;cpg.insert(cpg.end(),{"--meth=taps","--meth-context","cpg","--bam","-f",bam,fa,f1,f2});
    REQUIRE(invoke(cpg)==0);CHECK(records(bam)==context);
    CHECK(system(("rm -rf '"+path+"'").c_str())==0);
}

TEST_CASE("aln_meth: neutral conversion guards own opposite-strand CpG CIGARs in SAM and BAM") {
    std::string path=aln_test::tmp_template("aln_meth_conv_bam");
    REQUIRE(mkdtemp(&path[0])!=nullptr);
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
        args={"aln","--meth=5base","--meth-context",context,"-n","2","-l","1024","--bam","-f",bam,fa,fq};
        REQUIRE(invoke(args)==0);CHECK(records(bam)==a);
    }
    CHECK(system(("rm -rf '"+path+"'").c_str())==0);
}

TEST_CASE("aln_meth: all-resolved pairing beyond 500 preserves trimmed SAM and BAM ownership") {
    std::string path=aln_test::tmp_template("aln_meth_pair601");REQUIRE(mkdtemp(&path[0])!=nullptr);
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
    CHECK(system(("rm -rf '"+path+"'").c_str())==0);
}

TEST_CASE("aln_meth: cap-bound rescue restores 17000-base tails in paired SAM and BAM") {
    std::string path=aln_test::tmp_template("aln_meth_cap_bam");REQUIRE(mkdtemp(&path[0])!=nullptr);
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
    CHECK(system(("rm -rf '"+path+"'").c_str())==0);
}
