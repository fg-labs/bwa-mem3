// test/unit/test_xa_pe_emit.cpp — the paired emitter's XA slot selection,
// end to end.
//
// mem_sam_pe_batch_post builds an end's XA strings only for the slots it reads
// (the pair's primary z and the emitted ALT record alt; test_xa_lazy_slots.cpp
// checks mem_gen_alt's side of that against a full build), and reads XA and HN
// of those two slots. This drives the emitter itself, as SAM text and, under
// --meth, as BAM records, and checks the records' FLAG, RNAME, POS, XA and HN
// against values written out by hand:
//   - a pair whose chosen copy is often not the single-end primary (z != 0):
//     end 0's primary carries the other two copies as XA, and HN 2;
//   - a pair with an emitted ALT hit: the ALT record carries its own XA (its
//     two hits) and HN 2, the primary none;
//   - both again with the copies on read 2 (a reverse-strand end), so an
//     attachment that only works for read 1 also fails.
// Building or reading a slot other than z or alt drops or misplaces a tag.

// htslib headers first: they share the KSTRING_H include guard with
// bwa-mem3's kstring.h.
#include "htslib/kstring.h"
#include "htslib/sam.h"

#include "doctest/doctest.h"

#include "bam_writer.h"
#include "bwa.h"
#include "bwamem.h"
#include "meth_bam.h"
#include "meth_opts.h"
#include "meth_orig_ref_fixture.h"

#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

using meth_test::nt4;
using meth_test::random_bases;
using meth_test::revcomp;

namespace {

const int kRef = 6000, kLen = 100;

// A --meth option set scored with the genomic matrix.
void set_meth_genomic(mem_opt_t *opt) {
    opt->meth_mode = 1;
    opt->meth_scoring = MEM_METH_SCORING_GENOMIC;
    mem_opt_fill_meth_mat(opt);
}

// One emitted record's fields under test; xa is "" and hn -1 when absent.
struct Record {
    std::string flag_rname_pos, xa;
    int hn = -1;
};

// One paired-end read run through mem_sam_pe_batch_post, with rescue off and
// an FR model of insert 300 (window 180-420). Owns everything it allocates, so
// a failed REQUIRE leaks nothing. End 0 reads ref[1000, 1100), copied to each
// position in copies0; end 1 is the reverse complement of ref[2200, 2300),
// copied to each position in copies1. regs[1] starts with end 1's region at
// 2200 and regs[0] empty.
struct PairRun {
    std::string ref;
    meth_test::OrigRefFixture *fx = nullptr;
    mem_opt_t *opt = nullptr;
    bool meth;
    std::string r[2], q;
    std::vector<char> nt[2];
    char name[3] = "pp";
    bseq1_t sq[2];
    std::vector<mem_alnreg_t> regs[2];

    PairRun(bool meth_, const std::vector<int> &copies0, const std::vector<int> &copies1 = {})
        : meth(meth_), q(kLen, 'I') {
        ref = random_bases(kRef, 4242);
        const std::string seg = ref.substr(1000, kLen), seg1 = ref.substr(2200, kLen);
        for (int at : copies0) ref.replace(at, kLen, seg);
        for (int at : copies1) ref.replace(at, kLen, seg1);
        fx = new meth_test::OrigRefFixture(ref);
        opt = mem_opt_init();
        opt->flag |= MEM_F_NO_RESCUE | MEM_F_PE;
        if (meth) set_meth_genomic(opt);
        r[0] = seg;
        r[1] = revcomp(ref.substr(2200, kLen));
        memset(sq, 0, sizeof sq);
        for (int e = 0; e < 2; ++e) {
            for (char b : r[e]) nt[e].push_back((char)nt4(b));
            sq[e].name = name;
            sq[e].l_seq = kLen;
            sq[e].seq = nt[e].data();
            sq[e].qual = &q[0];
            if (meth) sq[e].meth_orig_seq = &r[e][0];
        }
        // the mate: reverse, forward span [2200, 2300)
        regs[1].push_back(region_rev(2200, kLen, false));
    }
    ~PairRun() {
        for (int e = 0; e < 2; ++e) {
            free(sq[e].sam);
            for (int k = 0; k < sq[e].n_bams; ++k) bam_writer_free((bam1_t *)sq[e].bams[k]);
            free(sq[e].bams);
        }
        free(opt);
        delete fx;
    }
    PairRun(const PairRun &) = delete;
    PairRun &operator=(const PairRun &) = delete;

    // A reverse-strand region over forward bases [at, at + kLen).
    mem_alnreg_t region_rev(int64_t at, int score, bool is_alt) const {
        return region(2 * (int64_t)kRef - (at + kLen), score, is_alt);
    }
    // A region over forward bases [rb, rb + kLen) of the read's full length.
    mem_alnreg_t region(int64_t rb, int score, bool is_alt) const {
        mem_alnreg_t g{};
        g.rb = rb; g.re = rb + kLen; g.qb = 0; g.qe = kLen; g.rid = 0;
        g.score = g.truesc = score; g.w = opt->w;
        g.seedcov = g.seedlen0 = kLen;
        g.secondary = g.secondary_all = -1;
        g.is_alt = is_alt;
        g.meth_hypothesis = meth ? 1 : -1;
        return g;
    }

    void run(uint64_t id) {
        mem_alnreg_v a[2];
        for (int e = 0; e < 2; ++e) a[e].n = a[e].m = regs[e].size(), a[e].a = regs[e].data(), a[e].capped_w = 0;
        mem_cache *mmc = (mem_cache *)calloc(1, sizeof(mem_cache));
        REQUIRE(mmc != NULL);
        int32_t gcnt = 0;
        const uint8_t *saved_pac = g_meth_orig_pac;
        if (meth) g_meth_orig_pac = fx->orig_pac;
        mem_sam_pe_batch_post(opt, fx->orig_bns, fx->orig_pac, pes(), id, sq, a, NULL, mmc, gcnt, 0);
        g_meth_orig_pac = saved_pac;
        free(mmc);
    }

    static const mem_pestat_t *pes() {
        static mem_pestat_t p[4];
        memset(p, 0, sizeof p);
        for (int d = 0; d < 4; ++d) p[d].failed = 1;
        p[1].failed = 0, p[1].avg = 300, p[1].std = 30, p[1].low = 180, p[1].high = 420;
        return p;
    }

    // End e's records in emission order: from the SAM text, or (--meth) from
    // the BAM records, RNAME read from the record either way.
    std::vector<Record> end(int e) const {
        std::vector<Record> out;
        if (!meth) {
            REQUIRE(sq[e].sam != NULL);
            const std::string text(sq[e].sam);
            for (size_t b = 0, nl; b < text.size(); b = nl + 1) {
                nl = text.find('\n', b);
                if (nl == std::string::npos) nl = text.size();
                const std::string line = text.substr(b, nl - b);
                std::vector<std::string> f;
                size_t c = 0;
                for (size_t t; (t = line.find('\t', c)) != std::string::npos; c = t + 1) f.push_back(line.substr(c, t - c));
                f.push_back(line.substr(c));
                REQUIRE(f.size() > 11);
                Record rec;
                rec.flag_rname_pos = f[1] + " " + f[2] + " " + f[3];
                for (size_t k = 11; k < f.size(); ++k) {
                    if (f[k].compare(0, 5, "XA:Z:") == 0) rec.xa = f[k].substr(5);
                    if (f[k].compare(0, 5, "HN:i:") == 0) rec.hn = std::atoi(f[k].c_str() + 5);
                }
                out.push_back(rec);
            }
        } else {
            for (int k = 0; k < sq[e].n_bams; ++k) {
                const bam1_t *b = (const bam1_t *)sq[e].bams[k];
                REQUIRE(b->core.tid >= 0);
                REQUIRE(b->core.tid < fx->orig_bns->n_seqs);
                Record rec;
                rec.flag_rname_pos = std::to_string(b->core.flag) + " " + fx->orig_bns->anns[b->core.tid].name + " "
                                     + std::to_string(b->core.pos + 1);
                if (const uint8_t *t = bam_aux_get(b, "XA")) rec.xa = bam_aux2Z(t);
                if (const uint8_t *t = bam_aux_get(b, "HN")) rec.hn = (int)bam_aux2i(t);
                out.push_back(rec);
            }
        }
        return out;
    }
};

}  // namespace

TEST_CASE("mem_sam_pe_batch_post: the pair's primary carries its own XA and HN, wherever it ranked" *
          doctest::test_suite("unit/xa_regen")) {
    // End 0 is three exact copies of the segment (at 1000, 2000, 3000); only the
    // copy at 2000 pairs (insert 300). With rescue off the pair's primary is
    // that copy wherever the single-end ranking put it, so z[0] is 0 for some
    // ids and not for others; its XA lists the other two copies in region order
    // and its HN counts them.
    for (bool meth : {false, true}) {
        int nonzero = 0, zero = 0;
        for (uint64_t id = 0; id < 16; ++id) {
            CAPTURE(meth);
            CAPTURE(id);
            PairRun run(meth, {2000, 3000});
            for (int at : {1000, 2000, 3000}) run.regs[0].push_back(run.region(at, kLen, false));
            run.run(id);
            // regs[0] is reordered in place: find the copy at 2000 and the others
            int chosen = -1;
            std::string want;
            for (int k = 0; k < 3; ++k) {
                if (run.regs[0][(size_t)k].rb == 2000) chosen = k;
                else want += "X,+" + std::to_string(run.regs[0][(size_t)k].rb + 1) + ",100M,0;";
            }
            REQUIRE(chosen >= 0);
            const std::vector<Record> recs = run.end(0);
            REQUIRE(recs.size() == 1);
            // proper pair, read 1, mate reverse: 0x1 | 0x2 | 0x20 | 0x40 = 99
            CHECK(recs[0].flag_rname_pos == "99 X 2001");
            CHECK(recs[0].xa == want);
            CHECK(recs[0].hn == 2);
            (chosen != 0 ? nonzero : zero)++;
        }
        CAPTURE(meth);
        CHECK(nonzero > 0);  // the switch ran: z[0] != 0
        CHECK(zero > 0);
    }
}

TEST_CASE("mem_sam_pe_batch_post: an emitted ALT hit carries its own XA and HN" *
          doctest::test_suite("unit/xa_regen")) {
    // End 0: the pairing copy at 2000 (score 100), an ALT hit at 4000 that
    // outscores it (101) and an ALT hit at 5000 (90), all the same read bases.
    // mem_mark_primary_se ranks the 4000 ALT hit first, so the other two are
    // its hits (secondary_all), and it stays a primary: the paired path emits
    // it as the read's ALT record (region n_pri = 1, mem_pe_alt_slot). The
    // pair's primary is the copy at 2000 (z = 0), which owns no hit: no XA,
    // HN 0. The ALT record's XA is its two hits in region order (2000, 5000),
    // HN 2.
    for (bool meth : {false, true}) {
        for (uint64_t id = 0; id < 4; ++id) {
            CAPTURE(meth);
            CAPTURE(id);
            PairRun run(meth, {2000, 4000, 5000});
            run.regs[0].push_back(run.region(2000, kLen, false));
            run.regs[0].push_back(run.region(4000, kLen + 1, true));
            run.regs[0].push_back(run.region(5000, 90, true));
            run.run(id);
            REQUIRE(run.regs[0][0].rb == 2000);
            REQUIRE(run.regs[0][1].rb == 4000);
            const std::vector<Record> recs = run.end(0);
            REQUIRE(recs.size() == 2);
            CHECK(recs[0].flag_rname_pos == "99 X 2001");
            CHECK(recs[0].xa == "");
            CHECK(recs[0].hn == 0);
            // 0x800 | 0x40 | 0x20 | 0x2 | 0x1 = 2147
            CHECK(recs[1].flag_rname_pos == "2147 X 4001");
            CHECK(recs[1].xa == "X,+2001,100M,0;X,+5001,100M,0;");
            CHECK(recs[1].hn == 2);
        }
    }
}

TEST_CASE("mem_sam_pe_batch_post: read 2's pair primary carries its own XA and HN, wherever it ranked" *
          doctest::test_suite("unit/xa_regen")) {
    // The first case mirrored onto read 2 (reverse strand): end 1 is three
    // exact copies (at 2200, 3400, 4400), end 0 one forward hit at 2000. Only
    // the copy at 2200 pairs (insert 300), so z[1] is 0 for some ids and not
    // for others; its XA lists the other two copies (reverse) in region order.
    for (bool meth : {false, true}) {
        int nonzero = 0, zero = 0;
        for (uint64_t id = 0; id < 16; ++id) {
            CAPTURE(meth);
            CAPTURE(id);
            PairRun run(meth, {2000}, {3400, 4400});
            run.regs[0].push_back(run.region(2000, kLen, false));
            run.regs[1].clear();
            for (int at : {2200, 3400, 4400}) run.regs[1].push_back(run.region_rev(at, kLen, false));
            run.run(id);
            int chosen = -1;
            std::string want;
            for (int k = 0; k < 3; ++k) {
                const int64_t at = 2 * (int64_t)kRef - run.regs[1][(size_t)k].rb - kLen;
                if (at == 2200) chosen = k;
                else want += "X,-" + std::to_string(at + 1) + ",100M,0;";
            }
            REQUIRE(chosen >= 0);
            const std::vector<Record> recs = run.end(1);
            REQUIRE(recs.size() == 1);
            // proper pair, read 2, reverse: 0x1 | 0x2 | 0x10 | 0x80 = 147
            CHECK(recs[0].flag_rname_pos == "147 X 2201");
            CHECK(recs[0].xa == want);
            CHECK(recs[0].hn == 2);
            (chosen != 0 ? nonzero : zero)++;
        }
        CAPTURE(meth);
        CHECK(nonzero > 0);  // the switch ran: z[1] != 0
        CHECK(zero > 0);
    }
}

TEST_CASE("mem_sam_pe_batch_post: read 2's emitted ALT hit carries its own XA and HN" *
          doctest::test_suite("unit/xa_regen")) {
    // The ALT case mirrored onto read 2: end 1's pairing copy at 2200 (score
    // 100), an ALT hit at 3400 that outscores it (101) and one at 4400 (90),
    // all reverse; end 0 one forward hit at 2000. The 3400 ALT hit is read 2's
    // ALT record with XA its two hits in region order (2200, 4400) and HN 2;
    // the pair's primary at 2200 has no XA and HN 0.
    for (bool meth : {false, true}) {
        for (uint64_t id = 0; id < 4; ++id) {
            CAPTURE(meth);
            CAPTURE(id);
            PairRun run(meth, {2000}, {3400, 4400});
            run.regs[0].push_back(run.region(2000, kLen, false));
            run.regs[1].clear();
            run.regs[1].push_back(run.region_rev(2200, kLen, false));
            run.regs[1].push_back(run.region_rev(3400, kLen + 1, true));
            run.regs[1].push_back(run.region_rev(4400, 90, true));
            run.run(id);
            REQUIRE(run.regs[1][0].rb == run.region_rev(2200, kLen, false).rb);
            REQUIRE(run.regs[1][1].rb == run.region_rev(3400, kLen, false).rb);
            const std::vector<Record> recs = run.end(1);
            REQUIRE(recs.size() == 2);
            CHECK(recs[0].flag_rname_pos == "147 X 2201");
            CHECK(recs[0].xa == "");
            CHECK(recs[0].hn == 0);
            // 0x800 | 0x80 | 0x10 | 0x2 | 0x1 = 2195
            CHECK(recs[1].flag_rname_pos == "2195 X 3401");
            CHECK(recs[1].xa == "X,-2201,100M,0;X,-4401,100M,0;");
            CHECK(recs[1].hn == 2);
        }
    }
}
