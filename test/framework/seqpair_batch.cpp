// test/framework/seqpair_batch.cpp
#include "seqpair_batch.h"

#include <algorithm>
#include <cstring>

#include "kswv.h"   // g_defr
#include "macro.h"  // SIMD_WIDTH8

namespace bwa_tests {

namespace {

inline void revseq(int l, uint8_t *s) {
    for (int i = 0; i < l >> 1; ++i) {
        uint8_t t = s[i];
        s[i] = s[l - 1 - i];
        s[l - 1 - i] = t;
    }
}

} // namespace

BatchBuffers::BatchBuffers(const std::vector<TestPair> &pairs, int xtra_flags)
    : n_(static_cast<int>(pairs.size())) {

    int32_t maxRefLen = 0, maxQerLen = 0;
    size_t  ref_total = 0, qer_total = 0;
    for (const auto &p : pairs) {
        ref_total += p.ref.size();
        qer_total += p.qry.size();
        if (static_cast<int32_t>(p.ref.size()) > maxRefLen) {
            maxRefLen = static_cast<int32_t>(p.ref.size());
        }
        if (static_cast<int32_t>(p.qry.size()) > maxQerLen) {
            maxQerLen = static_cast<int32_t>(p.qry.size());
        }
    }
    // SIMD tail padding: the kernel may read past the last pair when
    // lane-filling.
    ref_total += static_cast<size_t>(SIMD_WIDTH8) * maxRefLen;
    qer_total += static_cast<size_t>(SIMD_WIDTH8) * maxQerLen;

    seqBufRef_.assign(ref_total, 0);
    seqBufQer_.assign(qer_total, 0);
    pairs_.assign(n_ + SIMD_WIDTH8, SeqPair{});
    aln_.assign(n_ + SIMD_WIDTH8, g_defr);

    size_t ref_off = 0, qer_off = 0;
    for (int i = 0; i < n_; i++) {
        const auto &p = pairs[i];
        std::memcpy(seqBufRef_.data() + ref_off, p.ref.data(), p.ref.size());
        std::memcpy(seqBufQer_.data() + qer_off, p.qry.data(), p.qry.size());

        SeqPair sp = {};
        sp.idr   = static_cast<int32_t>(ref_off);
        sp.idq   = static_cast<int32_t>(qer_off);
        sp.id    = i;
        sp.len1  = static_cast<int32_t>(p.ref.size());
        sp.len2  = static_cast<int32_t>(p.qry.size());
        sp.h0    = xtra_flags;
        sp.seqid = i;
        sp.regid = i;
        pairs_[i] = sp;

        ref_off += p.ref.size();
        qer_off += p.qry.size();
    }
}

int BatchBuffers::prepare_phase1() {
    int pos = 0;
    for (int i = 0; i < n_; i++) {
        kswr_t  r    = aln_[i];
        SeqPair sp   = pairs_[i];
        int     xtra = sp.h0;
        if ((xtra & KSW_XSTART) == 0 ||
            ((xtra & KSW_XSUBO) && r.score < (xtra & 0xffff))) {
            continue;
        }
        sp.h0   = KSW_XSTOP | r.score;
        sp.len2 = r.qe + 1;
        uint8_t *qs = seqBufQer_.data() + sp.idq;
        uint8_t *rs = seqBufRef_.data() + sp.idr;
        revseq(r.qe + 1, qs);
        revseq(r.te + 1, rs);
        pairs_[pos++] = sp;
    }
    return pos;
}

} // namespace bwa_tests
