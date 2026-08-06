// test/unit/test_bam_aux_append_propagation.cpp — the BAM writers must not
// report success for a record that lost a tag.
//
// htslib's bam_aux_append() returns -1 (errno ENOMEM or EINVAL) when the
// bam1_t cannot grow. Every call site in bam_writer.cpp and meth_bam.cpp used
// to discard that return, so an allocation failure produced a record that was
// silently missing NM, MD, SA, pa, XA, XM... and a writer that still returned
// 0. Both functions already treat every OTHER allocation failure as `return
// -1` (the bam_cigar / seq_text / qual_bin / ks_resize guards), so the aux
// appends were the outlier.
//
// Forcing a real bam_aux_append() failure from a test is not practical, and
// the reason is worth writing down, because it is why these cases pin the
// contract rather than the failure. htslib's bam_aux_append (ext/htslib/sam.c)
// can only fail in two places:
//
//     new_len = b->l_data + 3 + len;
//     if (new_len > INT32_MAX || new_len < b->l_data) goto nomem;
//     if (realloc_bam_data(b, new_len) < 0) return -1;
//
// Reaching the first means driving `l_data + 3 + len` past INT32_MAX, and all
// three operands are `int` — so a test that manufactured it would be provoking
// signed-overflow UB inside htslib, not exercising our propagation, and at -O3
// the compiler is entitled to delete the very check it was aiming at. The
// second needs a ~2 GB realloc to actually fail, which no unit test can arrange
// portably. Both are true OOM conditions. The only remaining seam is an
// allocator shim or a fault-injecting indirection around every call site, and
// the latter would put a mutable global and an indirect call on a per-tag,
// per-record hot path purely to buy test reachability.
//
// What these cases pin instead is the contract that changed shape, on both of
// the writers this sweep touched:
//
//   - bam_writer_append_generic_aux() is now `int` and reports -1 rather than
//     returning void, so a caller can tell a complete record from a partial one.
//   - the success path still returns 0 and still emits every tag — the risk of
//     a sweep like this is a FALSE -1 that starts dropping good records, and
//     that risk is real on both bam_writer_append_generic_aux() and the
//     mem_aln_to_bam() record writer whose ~10 append sites all grew a
//     `< 0 -> return -1` arm. A false -1 there silently drops a whole record.
//
// htslib headers first: they share the KSTRING_H include guard with
// bwa-mem3's kstring.h, and the BAM writers are compiled against htslib's.
#include "htslib/sam.h"
#include "htslib/kstring.h"

#include "doctest/doctest.h"

#include "bam_writer.h"
#include "bwamem.h"

#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

namespace {

// A one-contig reference carrying an annotation, so the -V / MEM_F_REF_HDR
// XR:Z branch — the one site that owns a heap buffer across the append — is
// reachable.
struct RefFixture {
    bntann1_t ann{};
    bntseq_t  bns{};

    explicit RefFixture(const char *anno = "") {
        ann.offset = 0;
        ann.len    = 1000;
        ann.n_ambs = 0;
        ann.gi     = 0;
        ann.is_alt = 0;
        ann.name   = const_cast<char *>("chr1");
        ann.anno   = const_cast<char *>(anno);
        bns.l_pac  = 1000;
        bns.n_seqs = 1;
        bns.anns   = &ann;
    }
};

// One mapped read, optionally carrying FASTQ comment tags (-C).
struct ReadFixture {
    static const int kLen = 10;

    std::string       name{"read1"};
    std::string       comment;
    std::vector<char> seq;
    std::string       qual;
    bseq1_t           s{};

    explicit ReadFixture(const char *cmt = NULL)
        : comment(cmt ? cmt : ""), seq(kLen), qual(kLen, 'I') {
        for (int i = 0; i < kLen; ++i) seq[i] = (char)(i & 3);
        s.l_seq   = kLen;
        s.id      = 0;
        s.name    = const_cast<char *>(name.c_str());
        s.seq     = seq.data();
        s.qual    = const_cast<char *>(qual.c_str());
        s.comment = cmt ? const_cast<char *>(comment.c_str()) : NULL;
    }
};

// Out-of-class definition for kLen. In C++14 an in-class `static const int`
// initializer is not a definition, so binding it to a reference -- which any
// CHECK(x == ReadFixture::kLen) does -- is an odr-use that fails to link
// without this. Cheaper to define it once than to cast at every use site.
const int ReadFixture::kLen;

// The aux tag names of a bam1_t, in stored order.
std::vector<std::string> bam_tag_names(bam1_t *b) {
    std::vector<std::string> names;
    uint8_t *aux = bam_get_aux(b);
    uint8_t *end = b->data + b->l_data;
    while (aux != NULL && aux + 4 <= end) {
        names.push_back(std::string((const char *)aux, 2));
        uint8_t *next = bam_aux_next(b, aux + 2);
        if (next == NULL) break;
        aux = next - 2;
    }
    return names;
}

bool has_tag(const std::vector<std::string> &names, const char *t) {
    for (const std::string &n : names) if (n == t) return true;
    return false;
}

// One mapped alignment: a full-length match with every optional-tag source
// populated, so a single call reaches as many append sites as one record can.
//
// mem_aln_t keeps the MD string immediately AFTER the CIGAR array in the same
// allocation -- mem_aln_to_bam reads it as `(char *)(p.cigar + p.n_cigar)`,
// the layout mem_reg2aln produces. The fixture has to reproduce that, or the
// MD:Z append reads off the end of the CIGAR buffer.
struct AlnFixture {
    std::vector<uint32_t> cigar_and_md;
    std::string           xa;
    mem_aln_t             aln{};

    explicit AlnFixture(const char *md = "10", const char *xa_ = NULL)
        : xa(xa_ ? xa_ : "") {
        const size_t md_bytes = std::strlen(md) + 1;
        cigar_and_md.assign(1 + (md_bytes + 3) / 4, 0);
        // 10M: op 0 is 'M' in mem_aln_t's MIDSH=>01234 encoding.
        cigar_and_md[0] = ((uint32_t)ReadFixture::kLen << 4) | 0u;
        std::memcpy(&cigar_and_md[1], md, md_bytes);

        aln.pos     = 100;
        aln.rid     = 0;
        aln.flag    = 0;
        aln.is_rev  = 0;
        aln.is_alt  = 0;
        aln.mapq    = 60;
        aln.NM      = 1;
        aln.n_cigar = 1;
        aln.cigar   = cigar_and_md.data();
        aln.XA      = xa_ ? const_cast<char *>(xa.c_str()) : NULL;
        aln.HN      = 3;
        aln.score   = 100;
        aln.sub     = 20;
        // > 0, so the pa:f site is reached. 100/30 is deliberately not exact:
        // pa:f is defined as the float you get by round-tripping "%.3f" of the
        // ratio, so a ratio that quantizes visibly (3.333, not 3.33333...)
        // lets the assertion below check the value independently.
        aln.alt_sc  = 30;
        aln.meth_hypothesis = -1;
    }
};

}  // namespace

TEST_CASE("aux: append_generic_aux rejects null arguments"
          * doctest::test_suite("unit/bam_aux_append_propagation")) {
    // The guard has to answer -1 rather than 0, else a caller that now checks
    // the return would read "record complete" out of a call that did nothing.
    mem_opt_t  *opt = mem_opt_init();
    RefFixture  ref;
    ReadFixture read;

    bam1_t *b = bam_writer_alloc();
    REQUIRE(b != nullptr);

    CHECK(bam_writer_append_generic_aux(NULL, &read.s, opt, &ref.bns, 0) == -1);
    CHECK(bam_writer_append_generic_aux(b, NULL, opt, &ref.bns, 0) == -1);
    CHECK(bam_writer_append_generic_aux(b, &read.s, NULL, &ref.bns, 0) == -1);

    bam_writer_free(b);
    free(opt);
}

TEST_CASE("aux: append_generic_aux reports success and emits the comment tags"
          * doctest::test_suite("unit/bam_aux_append_propagation")) {
    // The sweep's real hazard is a false -1. Pin that the success path still
    // returns 0 AND still emits every tag it used to.
    mem_opt_t  *opt = mem_opt_init();
    RefFixture  ref;
    ReadFixture read("XX:Z:hello\tYY:i:42\tZZ:A:c");

    bam1_t *b = bam_writer_alloc();
    REQUIRE(b != nullptr);

    CHECK(bam_writer_append_generic_aux(b, &read.s, opt, &ref.bns, 0) == 0);

    const std::vector<std::string> names = bam_tag_names(b);
    CHECK(has_tag(names, "XX"));
    CHECK(has_tag(names, "YY"));
    CHECK(has_tag(names, "ZZ"));

    uint8_t *xx = bam_aux_get(b, "XX");
    REQUIRE(xx != nullptr);
    CHECK(std::string(bam_aux2Z(xx)) == "hello");
    uint8_t *yy = bam_aux_get(b, "YY");
    REQUIRE(yy != nullptr);
    CHECK(bam_aux2i(yy) == 42);

    bam_writer_free(b);
    free(opt);
}

TEST_CASE("aux: append_generic_aux succeeds on the heap-owning XR:Z branch"
          * doctest::test_suite("unit/bam_aux_append_propagation")) {
    // -V / MEM_F_REF_HDR is the one branch that mallocs a buffer and holds it
    // across the append, so it is the one where adding a failure return could
    // leak. The free was moved ahead of the return; this pins that the success
    // path still emits XR:Z with the TAB->SPACE rewrite intact.
    mem_opt_t  *opt = mem_opt_init();
    opt->flag |= MEM_F_REF_HDR;
    RefFixture  ref("desc\twith\ttabs");
    ReadFixture read;

    bam1_t *b = bam_writer_alloc();
    REQUIRE(b != nullptr);

    CHECK(bam_writer_append_generic_aux(b, &read.s, opt, &ref.bns, 0) == 0);

    uint8_t *xr = bam_aux_get(b, "XR");
    REQUIRE(xr != nullptr);
    CHECK(std::string(bam_aux2Z(xr)) == "desc with tabs");

    bam_writer_free(b);
    free(opt);
}

TEST_CASE("aux: an empty annotation emits no XR:Z and still succeeds"
          * doctest::test_suite("unit/bam_aux_append_propagation")) {
    // The malloc guard on that branch now returns -1 instead of skipping the
    // tag. That must not turn "there is nothing to emit" into a failure: the
    // branch is gated on a non-empty anno, so an empty one is still a 0.
    mem_opt_t  *opt = mem_opt_init();
    opt->flag |= MEM_F_REF_HDR;
    RefFixture  ref("");
    ReadFixture read;

    bam1_t *b = bam_writer_alloc();
    REQUIRE(b != nullptr);

    CHECK(bam_writer_append_generic_aux(b, &read.s, opt, &ref.bns, 0) == 0);
    CHECK(bam_aux_get(b, "XR") == nullptr);

    bam_writer_free(b);
    free(opt);
}

TEST_CASE("aux: mem_aln_to_bam rejects null arguments"
          * doctest::test_suite("unit/bam_aux_append_propagation")) {
    // Same reasoning as the generic_aux guard: -1, not 0. mem_aln_to_bam
    // already returned int before the sweep, but its callers only started
    // acting on the value once the aux appends could produce one, so the
    // guards are load-bearing now in a way they were not before.
    mem_opt_t  *opt = mem_opt_init();
    RefFixture  ref;
    ReadFixture read;
    AlnFixture  aln;

    bam1_t *b = bam_writer_alloc();
    REQUIRE(b != nullptr);

    CHECK(mem_aln_to_bam(NULL, opt, &ref.bns, &read.s, 1, &aln.aln, 0, NULL) == -1);
    CHECK(mem_aln_to_bam(b, NULL, &ref.bns, &read.s, 1, &aln.aln, 0, NULL) == -1);
    CHECK(mem_aln_to_bam(b, opt, &ref.bns, NULL, 1, &aln.aln, 0, NULL) == -1);
    CHECK(mem_aln_to_bam(b, opt, &ref.bns, &read.s, 1, NULL, 0, NULL) == -1);

    bam_writer_free(b);
    free(opt);
}

TEST_CASE("aux: mem_aln_to_bam reports success and emits every tag"
          * doctest::test_suite("unit/bam_aux_append_propagation")) {
    // The sweep added a `< 0 -> return -1` arm to ~10 append sites in this one
    // function. Its failure mode is not a missed error -- it is a FALSE -1,
    // which drops the whole record silently. This drives one alignment that
    // reaches every site reachable without a mate or a supplementary list, and
    // asserts both halves: the call still returns 0, and each tag is present
    // with the value it carried before.
    mem_opt_t  *opt = mem_opt_init();
    RefFixture  ref;
    ReadFixture read;
    AlnFixture  aln("5A4", "chr1,+200,10M,1;");

    bam1_t *b = bam_writer_alloc();
    REQUIRE(b != nullptr);

    CHECK(mem_aln_to_bam(b, opt, &ref.bns, &read.s, 1, &aln.aln, 0, NULL) == 0);

    uint8_t *nm = bam_aux_get(b, "NM");
    REQUIRE(nm != nullptr);
    CHECK(bam_aux2i(nm) == 1);

    uint8_t *md = bam_aux_get(b, "MD");
    REQUIRE(md != nullptr);
    CHECK(std::string(bam_aux2Z(md)) == "5A4");

    uint8_t *as = bam_aux_get(b, "AS");
    REQUIRE(as != nullptr);
    CHECK(bam_aux2i(as) == 100);

    uint8_t *xs = bam_aux_get(b, "XS");
    REQUIRE(xs != nullptr);
    CHECK(bam_aux2i(xs) == 20);

    uint8_t *xa = bam_aux_get(b, "XA");
    REQUIRE(xa != nullptr);
    CHECK(std::string(bam_aux2Z(xa)) == "chr1,+200,10M,1;");

    // HN:i is gated on the compat target; mem_opt_init() selects "off", which
    // emits it. Pin the gate as well as the value, so a compat-table edit that
    // silently drops HN from the native target fails here.
    REQUIRE(opt->compat != nullptr);
    CHECK(opt->compat->emit_hn == 1);
    uint8_t *hn = bam_aux_get(b, "HN");
    REQUIRE(hn != nullptr);
    CHECK(bam_aux2i(hn) == 3);

    // pa:f rides the same non-secondary block as SA:Z and is the one float
    // append in the function. The expected value is written out rather than
    // taken from bwa_pa_tag_value(): calling the producer to check its own
    // output would pass against a broken producer. score/alt_sc = 100/30
    // = 3.3333..., and pa:f is defined as that ratio round-tripped through
    // "%.3f", so the contract value is 3.333. The epsilon is tight enough to
    // tell the quantized answer from the raw ratio -- they differ by 3.3e-4.
    uint8_t *pa = bam_aux_get(b, "pa");
    REQUIRE(pa != nullptr);
    CHECK(bam_aux2f(pa) == doctest::Approx(3.333).epsilon(1e-6));

    // "Emits every tag" is the claim, so assert the whole set rather than
    // spot-checking members of it: an append site that silently stopped
    // firing, or a new one that appeared, both show up here and in neither
    // of the per-tag checks above. Order is part of it -- bam_writer.cpp
    // appends NM, MD, [MC, MQ,] AS, XS, [RG,] [SA,] pa, XA, HN in that
    // sequence to mirror mem_aln2sam, so the SAM and BAM paths agree field
    // for field. MC/MQ need a mate, RG needs -R, and SA needs a second
    // non-secondary alignment; none are present here.
    const std::vector<std::string> names = bam_tag_names(b);
    CHECK(names == std::vector<std::string>({"NM", "MD", "AS", "XS", "pa", "XA", "HN"}));

    // The record body, not just its tags. A false -1 is one failure mode of
    // the sweep; a 0 returned over a mangled record is the other, and only
    // the core fields catch it.
    CHECK(b->core.tid == 0);
    CHECK(b->core.pos == 100);
    CHECK(b->core.qual == 60);
    CHECK(b->core.l_qseq == ReadFixture::kLen);

    bam_writer_free(b);
    free(opt);
}

TEST_CASE("aux: mem_aln_to_bam still succeeds with a mate and a comment"
          * doctest::test_suite("unit/bam_aux_append_propagation")) {
    // MC:Z and MQ:i are only reachable with a mate, and the trailing
    // bam_writer_append_generic_aux() call is the seam where a -1 from the
    // helper becomes a -1 from the record writer. Driving both together pins
    // that the composed path reports success and carries the FASTQ-comment
    // tags through, which is the only way the helper's return is observable
    // from mem_aln_to_bam at all.
    mem_opt_t  *opt = mem_opt_init();
    RefFixture  ref;
    ReadFixture read("XX:Z:hello");
    AlnFixture  aln;
    AlnFixture  mate;
    mate.aln.pos  = 300;
    mate.aln.mapq = 42;

    bam1_t *b = bam_writer_alloc();
    REQUIRE(b != nullptr);

    CHECK(mem_aln_to_bam(b, opt, &ref.bns, &read.s, 1, &aln.aln, 0, &mate.aln) == 0);

    uint8_t *mc = bam_aux_get(b, "MC");
    REQUIRE(mc != nullptr);
    CHECK(std::string(bam_aux2Z(mc)) == "10M");

    REQUIRE(opt->compat != nullptr);
    CHECK(opt->compat->emit_mq == 1);
    uint8_t *mq = bam_aux_get(b, "MQ");
    REQUIRE(mq != nullptr);
    CHECK(bam_aux2i(mq) == 42);

    uint8_t *xx = bam_aux_get(b, "XX");
    REQUIRE(xx != nullptr);
    CHECK(std::string(bam_aux2Z(xx)) == "hello");

    bam_writer_free(b);
    free(opt);
}
