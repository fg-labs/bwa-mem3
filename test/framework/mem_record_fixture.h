// test/framework/mem_record_fixture.h
//
// Hand-built inputs for the record writers (mem_aln2sam, mem_aln_to_bam): a
// two-contig reference, a read and alignments, with no index behind them. The
// writers read only the reference's contig names, lengths and FASTA comments,
// the read's name, codes, quality and comment, and the alignment fields.
//
// Include htslib's kstring.h (and sam.h) before this header: they share the
// KSTRING_H include guard with bwa-mem3's kstring.h, and the BAM writer is
// compiled against htslib's. bwamem.h then declares mem_aln2sam against a
// kstring_t that mangles to a symbol libbwa.a does not define, so this header
// re-declares it against bwa-mem3's struct tag (same layout).

#ifndef BWA_TESTS_MEM_RECORD_FIXTURE_H
#define BWA_TESTS_MEM_RECORD_FIXTURE_H

#include "bntseq.h"
#include "bwamem.h"

#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <deque>
#include <string>
#include <vector>

#include "doctest/doctest.h"

struct __kstring_t {
    size_t l, m;
    char *s;
};
void mem_aln2sam(const mem_opt_t *opt, const bntseq_t *bns, __kstring_t *str, bseq1_t *s,
                 int n, const mem_aln_t *list, int which, const mem_aln_t *m_);

namespace bwa_tests {

// Two contigs, chrA (1000 bp, FASTA comment "first contig<TAB>x") and chrB
// (500 bp, no comment). Change names or annos, then call rebuild().
struct RecordRef {
    std::vector<std::string> names{"chrA", "chrB"};
    std::vector<std::string> annos{"first contig\tx", ""};
    std::vector<bntann1_t> anns;
    bntseq_t bns{};

    RecordRef() { rebuild(); }
    RecordRef(const RecordRef &) = delete;
    RecordRef &operator=(const RecordRef &) = delete;

    void rebuild() {
        anns.assign(names.size(), bntann1_t());
        int64_t off = 0;
        for (size_t i = 0; i < names.size(); ++i) {
            anns[i].offset = off;
            anns[i].len = i == 0 ? 1000 : 500;
            anns[i].name = const_cast<char *>(names[i].c_str());
            anns[i].anno = const_cast<char *>(annos[i].c_str());
            off += anns[i].len;
        }
        bns.l_pac = off;
        bns.n_seqs = (int)names.size();
        bns.anns = anns.data();
    }

    std::string header_text() const {
        std::string h;
        for (size_t i = 0; i < names.size(); ++i)
            h += "@SQ\tSN:" + names[i] + "\tLN:" + std::to_string(anns[i].len) + "\n";
        return h;
    }
};

// A read: bases as given, converted to the 0-4 codes the writers read; an
// empty quality string is FASTA input (no QUAL).
struct RecordRead {
    std::string name, bases, qual, comment;
    std::vector<char> codes;
    bseq1_t s{};

    RecordRead(const std::string &n, const std::string &b, const std::string &q, const char *cmt = NULL)
        : name(n), bases(b), qual(q), comment(cmt ? cmt : "") {
        for (char c : b) codes.push_back((char)nst_nt4_table[(unsigned char)c]);
        s.l_seq = (int)codes.size();
        s.name = const_cast<char *>(name.c_str());
        s.seq = codes.data();
        s.qual = q.empty() ? NULL : const_cast<char *>(qual.c_str());
        s.comment = cmt ? const_cast<char *>(comment.c_str()) : NULL;
        s.meth_base_ot = -1;
    }
    RecordRead(const RecordRead &) = delete;
    RecordRead &operator=(const RecordRead &) = delete;
};

// Builds mem_aln_t records. Each CIGAR lives in a buffer owned here, followed
// by its MD string, the layout mem_reg2aln produces and both writers read.
struct RecordAlns {
    std::deque<std::vector<uint32_t> > bufs;

    // cigar: "6M4S" etc. (M I D S H = 0..4, the mem_aln_t encoding).
    mem_aln_t make(int rid, int64_t pos0, bool rev, const std::string &cigar, const char *md, int nm,
                   int mapq, int score, int sub) {
        static const char kOps[] = "MIDSH";
        std::vector<uint32_t> ops;
        uint32_t len = 0;
        for (char c : cigar) {
            if (c >= '0' && c <= '9') {
                len = len * 10 + (uint32_t)(c - '0');
                continue;
            }
            const char *opc = strchr(kOps, c);
            REQUIRE(opc != NULL);
            ops.push_back(len << 4 | (uint32_t)(opc - kOps));
            len = 0;
        }
        const size_t md_bytes = strlen(md) + 1;
        bufs.push_back(std::vector<uint32_t>(ops.size() + (md_bytes + 3) / 4, 0));
        std::vector<uint32_t> &buf = bufs.back();
        if (!ops.empty()) memcpy(buf.data(), ops.data(), ops.size() * 4);
        memcpy(buf.data() + ops.size(), md, md_bytes);

        mem_aln_t a;
        memset(&a, 0, sizeof(a));
        a.pos = pos0;
        a.rid = rid;
        a.is_rev = rev ? 1 : 0;
        a.mapq = (uint32_t)mapq;
        a.NM = (uint32_t)nm;
        a.n_cigar = (int)ops.size();
        a.cigar = buf.data();
        a.XA = NULL;
        a.HN = -1;
        a.score = score;
        a.sub = sub;
        a.alt_sc = 0;
        a.meth_hypothesis = -1;
        return a;
    }
};

// mem_opt_init()'s options, freed at scope end.
struct OptGuard {
    mem_opt_t *opt;
    OptGuard() : opt(mem_opt_init()) {}
    ~OptGuard() { free(opt); }
    OptGuard(const OptGuard &) = delete;
    OptGuard &operator=(const OptGuard &) = delete;
};

}  // namespace bwa_tests

#endif  // BWA_TESTS_MEM_RECORD_FIXTURE_H
