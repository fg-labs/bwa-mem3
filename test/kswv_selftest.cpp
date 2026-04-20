// Self-consistency test for kswv::getScores8 versus scalar ksw_align2.
//
// Feeds deterministic random and edge-case SeqPair inputs through both the
// batched SIMD kernel (kswv) and the scalar reference (ksw_align2), then
// diffs kswr_t fields. Exits 0 on full match, nonzero on any divergence.
//
// Purpose: prove that the NEON kswv code path (kswv.cpp:175-628) is bit-
// identical to the scalar reference. NEON kswv is currently unreachable
// from production on ARM (gated by __AVX512BW__), so until we widen the
// gate, this is the only way to exercise it.
//
// Usage:
//   kswv_selftest                  # default: 10k random + ~40 edge cases
//   kswv_selftest <n_random> <seed>

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cstdint>
#include <cstdarg>
#include <random>
#include <vector>
#include <string>

#include "ksw.h"
#include "kswv.h"
#include "macro.h"

// Scoring (matches bwa-mem2 defaults)
static const int8_t SCORE_MATCH    = 1;
static const int8_t SCORE_MISMATCH = 4;
static const int    GAP_OPEN       = 6;
static const int    GAP_EXT        = 1;

// 5x5 substitution matrix for ACGTN encoding 0..4
static int8_t g_mat[25];

static void build_mat() {
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            g_mat[i*5 + j] = (i == j) ? SCORE_MATCH : -SCORE_MISMATCH;
        }
        g_mat[i*5 + 4] = 0;
        g_mat[4*5 + i] = 0;
    }
    g_mat[24] = 0;
}

struct TestPair {
    std::vector<uint8_t> ref;
    std::vector<uint8_t> qry;
    const char *tag;  // label for failure reporting
};

// --- pair generators --------------------------------------------------------

static TestPair gen_random(std::mt19937 &rng, int qlen, int rlen) {
    std::uniform_int_distribution<int> base(0, 3);
    TestPair p;
    p.ref.resize(rlen);
    p.qry.resize(qlen);
    p.tag = "random";
    for (int i = 0; i < rlen; i++) p.ref[i] = base(rng);
    for (int i = 0; i < qlen; i++) p.qry[i] = base(rng);
    return p;
}

static TestPair gen_exact_match(int len) {
    TestPair p;
    p.ref.assign(len * 2, 0);  // AAAAA...
    p.qry.assign(len, 0);
    p.tag = "exact_match";
    for (int i = 0; i < len * 2; i++) p.ref[i] = (i * 7 + 3) & 3;  // deterministic ACGT pattern
    for (int i = 0; i < len;     i++) p.qry[i] = p.ref[len / 2 + i]; // embed query in middle
    return p;
}

static TestPair gen_all_mismatch(int len) {
    TestPair p;
    p.ref.assign(len, 0);  // all A
    p.qry.assign(len, 1);  // all C
    p.tag = "all_mismatch";
    return p;
}

static TestPair gen_homopolymer(int len, uint8_t base) {
    TestPair p;
    p.ref.assign(len, base);
    p.qry.assign(len, base);
    p.tag = "homopolymer";
    return p;
}

static TestPair gen_with_indel(std::mt19937 &rng, int qlen, int rlen, int indel_pos, int indel_len) {
    auto p = gen_random(rng, qlen, rlen);
    p.tag = "with_indel";
    if (indel_pos + indel_len <= (int)p.qry.size()) {
        std::uniform_int_distribution<int> base(0, 3);
        for (int i = 0; i < indel_len; i++) p.qry[indel_pos + i] = base(rng);
    }
    return p;
}

static TestPair gen_with_n_bases(std::mt19937 &rng, int qlen, int rlen, int n_frac_pct) {
    auto p = gen_random(rng, qlen, rlen);
    p.tag = "with_n_bases";
    std::uniform_int_distribution<int> pct(0, 99);
    for (int i = 0; i < qlen; i++) if (pct(rng) < n_frac_pct) p.qry[i] = 4;
    for (int i = 0; i < rlen; i++) if (pct(rng) < n_frac_pct) p.ref[i] = 4;
    return p;
}

// --- runners ----------------------------------------------------------------

static kswr_t run_scalar(const TestPair &p) {
    // ksw_align2 expects non-const pointers, but does not mutate.
    return ksw_align2(p.qry.size(), const_cast<uint8_t*>(p.qry.data()),
                      p.ref.size(), const_cast<uint8_t*>(p.ref.data()),
                      5, g_mat, GAP_OPEN, GAP_EXT, GAP_OPEN, GAP_EXT,
                      KSW_XSTART | KSW_XSUBO, nullptr);
}

// Pack a batch of pairs into kswv's flat seqBuf layout and run getScores8.
// Returns aln array (length == pairs.size()); caller frees nothing.
struct BatchResult {
    std::vector<kswr_t> aln;
};

static BatchResult run_batched(const std::vector<TestPair> &pairs,
                               int32_t maxRefLen, int32_t maxQerLen) {
    BatchResult out;
    int n = (int)pairs.size();
    out.aln.assign(n + SIMD_WIDTH8, g_defr);

    // Flatten into buffers. kswv's internal kernel reads [idr, idr+len1) and
    // [idq, idq+len2) out of seqBufRef / seqBufQer.
    size_t ref_total = 0, qer_total = 0;
    for (const auto &p : pairs) { ref_total += p.ref.size(); qer_total += p.qry.size(); }
    // round up for SIMD tail padding
    ref_total += SIMD_WIDTH8 * maxRefLen;
    qer_total += SIMD_WIDTH8 * maxQerLen;

    std::vector<uint8_t> seqBufRef(ref_total, 0);
    std::vector<uint8_t> seqBufQer(qer_total, 0);
    std::vector<SeqPair> pairArray(n + SIMD_WIDTH8);

    size_t ref_off = 0, qer_off = 0;
    for (int i = 0; i < n; i++) {
        const auto &p = pairs[i];
        std::memcpy(seqBufRef.data() + ref_off, p.ref.data(), p.ref.size());
        std::memcpy(seqBufQer.data() + qer_off, p.qry.data(), p.qry.size());
        SeqPair sp = {};
        sp.idr    = (int32_t)ref_off;
        sp.idq    = (int32_t)qer_off;
        sp.id     = i;
        sp.len1   = (int32_t)p.ref.size();
        sp.len2   = (int32_t)p.qry.size();
        sp.h0     = KSW_XSTART | KSW_XSUBO;
        sp.seqid  = i;
        sp.regid  = i;
        pairArray[i] = sp;
        ref_off += p.ref.size();
        qer_off += p.qry.size();
    }

    kswv *pwsw = new kswv(GAP_OPEN, GAP_EXT, GAP_OPEN, GAP_EXT,
                          SCORE_MATCH, -SCORE_MISMATCH, 1, maxRefLen, maxQerLen);
    pwsw->getScores16(pairArray.data(), seqBufRef.data(), seqBufQer.data(),
                      out.aln.data(), n, 1, 0);
    delete pwsw;
    out.aln.resize(n);
    return out;
}

// --- diffing ----------------------------------------------------------------

struct Mismatch {
    int idx;
    std::string tag;
    kswr_t scalar, batched;
};

static bool kswr_eq_score(const kswr_t &a, const kswr_t &b) {
    return a.score == b.score;  // primary correctness bar
}

int main(int argc, char **argv) {
    int n_random = (argc > 1) ? atoi(argv[1]) : 10000;
    if (n_random < 0) n_random = 0;
    uint32_t seed = (argc > 2) ? (uint32_t)atoi(argv[2]) : 42u;

    build_mat();

    std::mt19937 rng(seed);
    std::vector<TestPair> pairs;
    pairs.reserve((size_t)n_random + 64);

    // Edge cases first (so failures point to them in index-order output).
    pairs.push_back(gen_exact_match(50));
    pairs.push_back(gen_exact_match(100));
    pairs.push_back(gen_all_mismatch(50));
    pairs.push_back(gen_all_mismatch(100));
    pairs.push_back(gen_homopolymer(50, 0));  // all A
    pairs.push_back(gen_homopolymer(50, 1));  // all C
    pairs.push_back(gen_homopolymer(50, 2));  // all G
    pairs.push_back(gen_homopolymer(50, 3));  // all T
    pairs.push_back(gen_with_indel(rng, 100, 150, 40, 3));
    pairs.push_back(gen_with_indel(rng, 100, 150, 40, 10));
    pairs.push_back(gen_with_n_bases(rng, 100, 150, 5));
    pairs.push_back(gen_with_n_bases(rng, 100, 150, 20));

    // Minimum-length sanity
    pairs.push_back(gen_random(rng, 20, 40));
    pairs.push_back(gen_random(rng, 10, 30));

    // Bulk random
    std::uniform_int_distribution<int> qlen_d(50, 128);
    std::uniform_int_distribution<int> rlen_d(100, 250);
    for (int i = 0; i < n_random; i++) {
        pairs.push_back(gen_random(rng, qlen_d(rng), rlen_d(rng)));
    }

    // Compute maxRefLen/maxQerLen for kswv sizing
    int32_t maxRefLen = 0, maxQerLen = 0;
    for (const auto &p : pairs) {
        if ((int32_t)p.ref.size() > maxRefLen) maxRefLen = p.ref.size();
        if ((int32_t)p.qry.size() > maxQerLen) maxQerLen = p.qry.size();
    }

    // Reference (scalar)
    std::vector<kswr_t> scalar_aln(pairs.size());
    for (size_t i = 0; i < pairs.size(); i++) scalar_aln[i] = run_scalar(pairs[i]);

    // Batched (NEON kswv on ARM; AVX-512 kswv on x86)
    BatchResult br = run_batched(pairs, maxRefLen, maxQerLen);

    // Diff
    std::vector<Mismatch> mism;
    for (size_t i = 0; i < pairs.size(); i++) {
        if (!kswr_eq_score(scalar_aln[i], br.aln[i])) {
            Mismatch m{ (int)i, pairs[i].tag, scalar_aln[i], br.aln[i] };
            mism.push_back(m);
        }
    }

    fprintf(stderr, "kswv_selftest: %zu pairs tested (%d random + %zu edge), %zu score mismatch\n",
            pairs.size(), n_random, pairs.size() - (size_t)n_random, mism.size());

    if (!mism.empty()) {
        int printed = 0;
        for (const auto &m : mism) {
            if (printed >= 20) {
                fprintf(stderr, "  ... %zu more\n", mism.size() - (size_t)printed);
                break;
            }
            ++printed;
            fprintf(stderr, "  [%d] tag=%s scalar.score=%d batched.score=%d  "
                            "scalar{tb=%d,te=%d,qb=%d,qe=%d}  batched{tb=%d,te=%d,qb=%d,qe=%d}\n",
                    m.idx, m.tag.c_str(),
                    m.scalar.score, m.batched.score,
                    m.scalar.tb, m.scalar.te, m.scalar.qb, m.scalar.qe,
                    m.batched.tb, m.batched.te, m.batched.qb, m.batched.qe);
        }
        return 1;
    }

    fprintf(stderr, "kswv_selftest: OK\n");
    return 0;
}
