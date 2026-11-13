#include "doctest/doctest.h"
#include "aln/aln_meth_realign.h"
#include "bwa.h"
#include "ksw.h"
#include "meth_model.h"
#include <algorithm>
#include <climits>
#include <cstdlib>
#include <string>
#include <vector>

namespace {
namespace ba = bwa_aln;
struct Ref {
    bntseq_t bns = {};
    bntann1_t ann = {};
    std::vector<uint8_t> pac;
    ba::aln_meth_ref_t ref;
    explicit Ref(const std::vector<uint8_t> &bases) : pac(bases.size()/4+1, 0) {
        ann.len = (int)bases.size(); bns.l_pac = bases.size(); bns.n_seqs = 1; bns.anns = &ann;
        for (size_t i=0;i<bases.size();++i) pac[i/4] |= bases[i] << ((~i&3)*2);
        ref.bns=&bns; ref.pac=pac.data();
    }
};
// Exhaustive legal M/I/D histories, independent of ksw and certificate implementation.
int oracle(const std::vector<uint8_t> &q, const std::vector<uint8_t> &t, const int8_t *mat,
           int r=0, int x=0, int state=0, int m=5) {
    if (r==(int)q.size() && x==(int)t.size()) return 0;
    int best=INT_MIN/4;
    if(r<(int)q.size() && x<(int)t.size())
        best=mat[t[x]*m+q[r]]+oracle(q,t,mat,r+1,x+1,0,m);
    if(r<(int)q.size() && state!=2)
        best=std::max(best,-(state==1?1:6)+oracle(q,t,mat,r+1,x,1,m));
    if(x<(int)t.size() && state!=1)
        best=std::max(best,-(state==2?1:6)+oracle(q,t,mat,r,x+1,2,m));
    return best;
}
int cigar_score(const std::vector<uint8_t>&q,const std::vector<uint8_t>&t,const int8_t*mat,int m,
                const std::vector<uint32_t>&cigar) {
    size_t r=0,x=0;int score=0;
    for(uint32_t c:cigar) {
        int n=c>>4,op=c&15;
        if(op==0) for(int i=0;i<n;++i) { REQUIRE(r<q.size()); REQUIRE(x<t.size()); score+=mat[t[x++]*m+q[r++]]; }
        else { REQUIRE((op==1||op==2)); score-=5+n; if(op==1)r+=n;else x+=n; }
    }
    CHECK(r==q.size());CHECK(x==t.size());return score;
}
}
TEST_SUITE_BEGIN("unit/aln_meth_realign");
TEST_CASE("aln_realign: tiny independent reward oracle and scalar dispatcher parity") {
    for(int mode : {MEM_METH_SCORING_COLLAPSED,MEM_METH_SCORING_GENOMIC,MEM_METH_SCORING_NEUTRAL}) {
        ba::aln_meth_opt_t mo;mo.m.meth_scoring=mode;ba::aln_meth_opt_finish(&mo);
        for(int bits=0;bits<64;++bits) {
            std::vector<uint8_t> q(3),t(3);
            for(int i=0;i<3;++i){q[i]=(bits>>i&1)?3:1;t[i]=(bits>>(i+3)&1)?3:1;}
            Ref ref(t);ba::aln_meth_candidate_t c{0,0,1,0,3};ba::aln_meth_alignment_t a;
            REQUIRE(ba::aln_meth_realign_window(ref.ref,mo,c,q.data(),3,0,0,&a)==ba::aln_meth_eval_status::ok);
            CHECK(a.reward==oracle(q,t,mo.mat_ot));
            CHECK(a.reward==cigar_score(q,t,mo.mat_ot,5,a.raw_cigar));
            int n=0;uint32_t*raw=nullptr;
            CHECK(a.reward==ksw_global2_scalar_ref(3,q.data(),3,t.data(),5,mo.mat_ot,5,1,5,1,3,&n,&raw));
            const std::vector<uint32_t> expected(raw,raw+n);free(raw);CHECK(a.raw_cigar==expected);
        }
    }
}
TEST_CASE("aln_realign: opposing gaps escape zero net delta hint") {
    std::vector<uint8_t> q,t;
    uint32_t seed=3197;
    for(int i=0;i<80;++i){seed=seed*1664525u+1013904223u;t.push_back(seed>>30);}
    q=t;q.insert(q.begin()+15,8,0);q.erase(q.begin()+60,q.begin()+68);
    Ref ref(t);ba::aln_meth_opt_t mo;mo.m.meth_scoring=MEM_METH_SCORING_GENOMIC;ba::aln_meth_opt_finish(&mo);
    ba::aln_meth_alignment_t a;ba::aln_meth_candidate_t c{0,0,1,0,(int64_t)t.size()};
    REQUIRE(ba::aln_meth_realign_window(ref.ref,mo,c,q.data(),(int)q.size(),0,0,&a)==ba::aln_meth_eval_status::ok);
    int n=0;uint32_t*raw=nullptr;
    int exact=ksw_global2_scalar_ref(q.size(),q.data(),t.size(),t.data(),5,mo.mat_ot,5,1,5,1,80,&n,&raw);
    int narrow=ksw_global(q.size(),q.data(),t.size(),t.data(),5,mo.mat_ot,5,1,3,nullptr,nullptr);
    CHECK(exact>narrow);CHECK(a.band>3);
    CHECK(a.reward==exact);CHECK(a.raw_cigar==std::vector<uint32_t>(raw,raw+n));free(raw);
    CHECK(a.reward==cigar_score(q,t,mo.mat_ot,5,a.raw_cigar));
}
TEST_CASE("aln_realign: terminal normalization contig checks and resource refusal") {
    std::vector<uint8_t> t(30,0),q(40,0);Ref ref(t);ba::aln_meth_opt_t mo;ba::aln_meth_opt_finish(&mo);
    ba::aln_meth_candidate_t c{0,0,1,0,30};ba::aln_meth_alignment_t a;
    REQUIRE(ba::aln_meth_realign_window(ref.ref,mo,c,q.data(),40,10,0,&a)==ba::aln_meth_eval_status::ok);
    CHECK(a.reward==30-15);CHECK(a.raw_cigar!=a.cigar);
    CHECK(a.rb==0);CHECK(a.re==30);
    bool clipped=false;for(auto op:a.cigar)clipped|=(op&15)==4;CHECK(clipped);
    c.re=31;CHECK(ba::aln_meth_realign_window(ref.ref,mo,c,q.data(),40,0,0,&a)==ba::aln_meth_eval_status::invalid_input);
    c.re=30;CHECK(ba::aln_meth_realign_window(ref.ref,mo,c,q.data(),INT_MAX,0,0,&a)==ba::aln_meth_eval_status::invalid_input);
    std::vector<uint8_t> huge(20000,0);Ref large(huge);c.re=20000;
    CHECK(ba::aln_meth_realign_window(large.ref,mo,c,huge.data(),20000,0,0,&a)==ba::aln_meth_eval_status::ok);
    mo.m.meth_context=METH_CTX_CPG;
    CHECK(ba::aln_meth_realign_window(large.ref,mo,c,huge.data(),20000,0,0,&a)==ba::aln_meth_eval_status::resource_error);
}

TEST_CASE("aln_realign: unequal tiny windows match independent oracle") {
    ba::aln_meth_opt_t mo;
    mo.m.meth_scoring=MEM_METH_SCORING_NEUTRAL;
    ba::aln_meth_opt_finish(&mo);
    for(int qn=1;qn<=4;++qn) for(int tn=1;tn<=4;++tn) {
        std::vector<uint8_t> q(qn),t(tn);
        for(int i=0;i<qn;++i)q[i]=(i+1)%4;
        for(int i=0;i<tn;++i)t[i]=i%4;
        Ref ref(t);ba::aln_meth_alignment_t a;ba::aln_meth_candidate_t c{0,1,1,0,tn};
        REQUIRE(ba::aln_meth_realign_window(ref.ref,mo,c,q.data(),qn,0,0,&a)==ba::aln_meth_eval_status::ok);
        CHECK(a.reward==oracle(q,t,mo.mat_ot));
        CHECK(a.reward==cigar_score(q,t,mo.mat_ot,5,a.raw_cigar));
    }
}
TEST_CASE("aln_realign: CpG full-band fallback includes outside-window neighbors") {
    std::vector<uint8_t> t={0,1,2,3,1,0,2,0},q={3,2,3,3,0,2};
    Ref ref(t);ba::aln_meth_opt_t mo;mo.m.meth_context=METH_CTX_CPG;
    mo.m.meth_scoring=MEM_METH_SCORING_NEUTRAL;ba::aln_meth_opt_finish(&mo);
    ba::aln_meth_candidate_t c{0,0,1,1,7};ba::aln_meth_alignment_t a;
    REQUIRE(ba::aln_meth_realign_window(ref.ref,mo,c,q.data(),q.size(),0,0,&a)==ba::aln_meth_eval_status::ok);
    std::vector<uint8_t> labeled={METH_C5,METH_G6,3,1,0,2};
    int n=0;uint32_t*raw=nullptr;
    CHECK(a.band==6);
    CHECK(a.reward==ksw_global2_scalar_ref(6,q.data(),6,labeled.data(),7,mo.mat7_ot,5,1,5,1,6,&n,&raw));
    CHECK(a.raw_cigar==std::vector<uint32_t>(raw,raw+n));free(raw);
    CHECK(a.reward==cigar_score(q,labeled,mo.mat7_ot,7,a.raw_cigar));
    // The G just inside the window has its CpG C immediately outside it.
    c.rb=2;c.re=7;q={0,3,1,0,2};
    REQUIRE(ba::aln_meth_realign_window(ref.ref,mo,c,q.data(),q.size(),0,0,&a)==ba::aln_meth_eval_status::ok);
    CHECK(a.reward==1); // four matches and -3 G/A, not freed for OT
}
TEST_CASE("aln_realign: trailing and leading deletions adjust display coordinates") {
    ba::aln_meth_opt_t mo;ba::aln_meth_opt_finish(&mo);
    std::vector<uint8_t> t(40,0),q(30,0);Ref ref(t);
    ba::aln_meth_alignment_t a;ba::aln_meth_candidate_t c{0,1,1,0,40};
    REQUIRE(ba::aln_meth_realign_window(ref.ref,mo,c,q.data(),q.size(),0,10,&a)==ba::aln_meth_eval_status::ok);
    CHECK(a.reward==15);
    CHECK(a.re-a.rb==30);
    CHECK(a.rb==10); // ksw's deterministic tie places the deletion at the leading end
    CHECK(a.re==40);
    CHECK(a.reward==cigar_score(q,t,mo.mat_ot,5,a.raw_cigar));
    std::vector<uint8_t> leading(10,2);leading.insert(leading.end(),30,0);Ref ref2(leading);
    REQUIRE(ba::aln_meth_realign_window(ref2.ref,mo,c,q.data(),q.size(),0,10,&a)==ba::aln_meth_eval_status::ok);
    CHECK(a.rb==10);CHECK(a.re==40);CHECK(a.reward==15);
    std::vector<uint8_t> trailing(30,0);trailing.insert(trailing.end(),10,2);Ref ref3(trailing);
    REQUIRE(ba::aln_meth_realign_window(ref3.ref,mo,c,q.data(),q.size(),0,10,&a)==ba::aln_meth_eval_status::ok);
    CHECK(a.rb==0);CHECK(a.re==30);CHECK(a.reward==15);
}
TEST_CASE("aln_realign: ambiguous query and original packed reference use stock matrix") {
    Ref ref({0,1,2,3});std::vector<uint8_t> q={0,4,2,3};
    ba::aln_meth_opt_t mo;ba::aln_meth_opt_finish(&mo);
    ba::aln_meth_candidate_t c{0,0,1,0,4};ba::aln_meth_alignment_t a;
    REQUIRE(ba::aln_meth_realign_window(ref.ref,mo,c,q.data(),q.size(),0,0,&a)==ba::aln_meth_eval_status::ok);
    CHECK(a.reward==2); // three +1 matches and the existing -1 ambiguous cell
    q[0]=5;
    CHECK(ba::aln_meth_realign_window(ref.ref,mo,c,q.data(),q.size(),0,0,&a)==ba::aln_meth_eval_status::invalid_input);
    c.rb=3;c.re=3;q[0]=0;
    CHECK(ba::aln_meth_realign_window(ref.ref,mo,c,q.data(),q.size(),0,0,&a)==ba::aln_meth_eval_status::invalid_input);
}

TEST_CASE("aln_realign: wide dispatcher calls execute existing SIMD when tier supports it") {
    std::vector<uint8_t> q(150,0),t(150,2);Ref ref(t);
    ba::aln_meth_opt_t mo;mo.m.meth_scoring=MEM_METH_SCORING_NEUTRAL;ba::aln_meth_opt_finish(&mo);
    ba::aln_meth_candidate_t c{0,0,1,0,150};ba::aln_meth_alignment_t a;
    const unsigned long before=ksw_g2_wave_exec_count();
    REQUIRE(ba::aln_meth_realign_window(ref.ref,mo,c,q.data(),q.size(),0,0,&a)==ba::aln_meth_eval_status::ok);
    const int w16=ksw_g2_wave16_wmin(),w32=ksw_g2_wave_wmin();
    if(w16>0&&w32>0&&a.band>=std::max(w16,w32))CHECK(ksw_g2_wave_exec_count()>before);
    int n=0;uint32_t*raw=nullptr;
    CHECK(a.reward==ksw_global2_scalar_ref(150,q.data(),150,t.data(),5,mo.mat_ot,5,1,5,1,150,&n,&raw));
    CHECK(a.raw_cigar==std::vector<uint32_t>(raw,raw+n));free(raw);
}

TEST_CASE("aln_realign: m7 full fallback matches independent oracle without SIMD") {
    for(int mode:{MEM_METH_SCORING_GENOMIC,MEM_METH_SCORING_NEUTRAL})
    for(int hyp:{0,1}) for(int bits=0;bits<64;++bits) {
        std::vector<uint8_t> target={1,2,1,0},query(4);
        for(int k=0;k<4;++k) query[k]=(bits>>(k*2))&3;
        Ref ref(target);ba::aln_meth_opt_t mo;mo.m.meth_context=METH_CTX_CPG;mo.m.meth_scoring=mode;
        ba::aln_meth_opt_finish(&mo);
        ba::aln_meth_candidate_t candidate{0,0,hyp,0,4};ba::aln_meth_alignment_t a;
        const unsigned long before=ksw_g2_wave_exec_count();
        REQUIRE(ba::aln_meth_realign_window(ref.ref,mo,candidate,query.data(),4,0,0,&a)==ba::aln_meth_eval_status::ok);
        CHECK(ksw_g2_wave_exec_count()==before);CHECK(a.band==4);
        const std::vector<uint8_t> labeled={METH_C5,METH_G6,1,0};
        const int8_t *mat=hyp?mo.mat7_ot:mo.mat7_ob;
        CHECK(a.reward==oracle(query,labeled,mat,0,0,0,7));
        CHECK(a.reward==cigar_score(query,labeled,mat,7,a.raw_cigar));
        int n=0;uint32_t *raw=nullptr;
        CHECK(a.reward==ksw_global2_scalar_ref(4,query.data(),4,labeled.data(),7,mat,5,1,5,1,4,&n,&raw));
        CHECK(a.raw_cigar==std::vector<uint32_t>(raw,raw+n));free(raw);
    }
}

TEST_CASE("aln_realign: m7 context cannot borrow contig neighbors or original N runs") {
    Ref ref({1,2,1,2});bntann1_t anns[2]={};anns[0].len=1;anns[1].offset=1;anns[1].len=3;
    ref.bns.anns=anns;ref.bns.n_seqs=2;
    ba::aln_meth_opt_t mo;mo.m.meth_context=METH_CTX_CPG;mo.m.meth_scoring=MEM_METH_SCORING_NEUTRAL;
    ba::aln_meth_opt_finish(&mo);ba::aln_meth_alignment_t a;
    uint8_t t=3,g=0,n=4;
    REQUIRE(ba::aln_meth_realign_window(ref.ref,mo,{0,0,1,0,1},&t,1,0,0,&a)==ba::aln_meth_eval_status::ok);
    CHECK(a.reward==-3); // terminal C cannot borrow the next contig's initial G
    REQUIRE(ba::aln_meth_realign_window(ref.ref,mo,{1,1,0,1,2},&g,1,0,0,&a)==ba::aln_meth_eval_status::ok);
    CHECK(a.reward==-3); // initial G cannot borrow preceding contig C
    REQUIRE(ba::aln_meth_realign_window(ref.ref,mo,{1,0,1,2,3},&t,1,0,0,&a)==ba::aln_meth_eval_status::ok);
    CHECK(a.reward==0); // true same-contig partner beyond window
    bntamb1_t ambiguous={};ambiguous.offset=3;ambiguous.len=1;ambiguous.amb='N';
    ref.bns.ambs=&ambiguous;ref.bns.n_holes=1;
    REQUIRE(ba::aln_meth_realign_window(ref.ref,mo,{1,0,1,2,3},&t,1,0,0,&a)==ba::aln_meth_eval_status::ok);
    CHECK(a.reward==-3); // PAC G is inside an original N run
    REQUIRE(ba::aln_meth_realign_window(ref.ref,mo,{1,0,1,2,3},&n,1,0,0,&a)==ba::aln_meth_eval_status::ok);
    CHECK(a.reward==-1); // original query ambiguity uses the N matrix column
}
