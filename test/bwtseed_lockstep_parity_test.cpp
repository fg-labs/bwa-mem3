/* Parity harness: scalar bwtSeedStrategyAllPosOneThread vs the lockstep
 * variant. Both functions are given identical inputs; outputs (returned
 * numSeed, then every SMEM field but the padding) must be identical. The one
 * exception is a lockstep seed taken from a BsdFirstWalk record that the SMEM
 * walk's unique-match text compare wrote (Case 13h), whose `l` may be stale and
 * is not compared (see l_may_be_stale).
 *
 * Both drivers end a read once no later walk can emit (x_end), so a wrong
 * bound shared by the two would pass the scalar-vs-lockstep compare. Every
 * case therefore also checks both against reference_bwtseed below, the scalar
 * walk as it was before that bound: every outer x up to the read end.
 *
 * Usage: bwtseed_lockstep_parity_test <bwa-mem3 index prefix>
 *   e.g. bwtseed_lockstep_parity_test test/fixtures/phix.fa
 *
 * The checked-in phiX fixture (test/fixtures/phix.fa) is REQUIRED, not just an
 * example: the emission-path cases (10p-12p) assert phiX-specific minimum seed
 * counts, so a run against a different index (e.g. hg38) can legitimately fail
 * those floors. run_unit_tests.sh builds and indexes this fixture.
 */

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <map>
#include <string>
#include <vector>
/* reference_bwtseed needs count[] and backwardExt. */
#define private public
#include "FMI_search.h"
#undef private
#include "seed_jump_test_util.h"
#include "lockstep_parity_test_util.h"
#include "bwa.h"

/* A degenerate N=1 build would invoke bwtSeedStrategyAllPosOneThread_lockstep
 * as a single-slot path, so the cases below would trivially pass without ever
 * exercising slot-interleaving. Fail the build early in that case. */
#if BWTSEED_LOCKSTEP_N <= 1
#error "bwtseed_lockstep_parity_test requires BWTSEED_LOCKSTEP_N > 1 to exercise the lockstep path"
#endif

/* Whether lockstep seed b, compared with the scalar or reference seed a, may
 * carry a stale `l`: l_hint is the hint the lockstep call used (nullptr when it
 * had none or rejected it), and b is that read's first seed (m = 0) taken from
 * a record with s = 1. Such a record may have been written by the SMEM walk's
 * unique-match text compare (ls_s1_compare), which never steps l, so the seed
 * keeps the l of the step that handed off to it. A record from an ordinary
 * forward step with s = 1 cannot be told apart here and is skipped too; every
 * other seed compares l. */
static bool l_may_be_stale(const SMEM &a, const SMEM &b, const BsdFirstWalkHint *l_hint) {
    if (l_hint == nullptr || a.rid != b.rid || a.m != 0 || b.m != 0) return false;
    const BsdFirstWalk &r = l_hint->rec[a.rid];
    return r.x > 0 && r.s == 1;
}

/* Field-by-field SMEM equality. memcmp is unsuitable because the SMEM struct
 * has implicit padding between `n` and `k` that upstream leaves
 * uninitialized — the padding bytes are not part of the semantic contract.
 * Neither driver tags its seeds, so both must leave info 0: SA resolution
 * trusts a tag in info (get_sa_entries_prefetch's trust_s1_tag).
 * `l` is skipped only where l_may_be_stale. */
static bool smem_fields_equal(const SMEM &a, const SMEM &b, const BsdFirstWalkHint *l_hint) {
    return a.rid == b.rid &&
           a.m   == b.m   &&
           a.n   == b.n   &&
           a.k   == b.k   &&
           (a.l == b.l || l_may_be_stale(a, b, l_hint)) &&
           a.s   == b.s   &&
           a.info == 0 && b.info == 0;
}

static int64_t smem_array_first_mismatch(const SMEM *a, const SMEM *b, int64_t n,
                                         const BsdFirstWalkHint *l_hint) {
    for (int64_t i = 0; i < n; i++) {
        if (!smem_fields_equal(a[i], b[i], l_hint)) return i;
    }
    return -1;
}

/* The scalar third-round walk with no early end: every outer x in
 * [0, readlength) is walked, so a read only ends at its last base. It is the
 * scalar driver's walk from before the x_end bound, less the prefetch hint,
 * and is the oracle that bound is checked against. The length test keeps the
 * driver's unsigned compare (SMEM m and n are uint32_t), so a negative
 * minSeedLen emits nothing here, as it does there. */
static int64_t reference_bwtseed(const FMI_search *fmi, const uint8_t *enc_qdb,
                                 const int32_t *max_intv_array, int32_t numReads,
                                 const bseq1_t *seq_, const int32_t *query_cum_len_ar,
                                 int32_t minSeedLen, SMEM *matchArray)
{
    int64_t numTotalSeed = 0;
    for (int32_t i = 0; i < numReads; i++) {
        const int readlength = seq_[i].l_seq;
        const int offset = query_cum_len_ar[i];
        int32_t x = 0;
        while (x < readlength) {
            int next_x = x + 1;
            SMEM smem;
            memset(&smem, 0, sizeof smem);
            smem.rid = i;
            smem.m = x;
            smem.n = x;
            uint8_t a = enc_qdb[offset + x];
            if (a < 4) {
                smem.k = fmi->count[a];
                smem.l = fmi->count[3 - a];
                smem.s = fmi->count[a + 1] - fmi->count[a];
                for (int j = x + 1; j < readlength; j++) {
                    next_x = j + 1;
                    a = enc_qdb[offset + j];
                    if (a >= 4) break;
                    /* Forward extension: backward extension on the swapped interval. */
                    SMEM smem_ = smem;
                    smem_.k = smem.l;
                    smem_.l = smem.k;
                    SMEM newSmem_ = fmi->backwardExt(smem_, 3 - a);
                    SMEM newSmem = newSmem_;
                    newSmem.k = newSmem_.l;
                    newSmem.l = newSmem_.k;
                    newSmem.n = j;
                    smem = newSmem;
                    if (smem.s < max_intv_array[i] && (smem.n - smem.m + 1) >= (uint32_t)minSeedLen) {
                        if (smem.s > 0) matchArray[numTotalSeed++] = smem;
                        break;
                    }
                }
            }
            x = next_x;
        }
    }
    return numTotalSeed;
}

/* Every field (rid, m, n, k, l, s, info) of two SMEM streams of length n, for
 * the SMEM pass (Case 13h), which tags info; returns the first differing index
 * or -1. */
static int64_t smem_all_fields_first_mismatch(const SMEM *a, const SMEM *b, int64_t n) {
    for (int64_t i = 0; i < n; i++) {
        if (a[i].rid != b[i].rid || a[i].m != b[i].m || a[i].n != b[i].n || a[i].k != b[i].k ||
            a[i].l != b[i].l || a[i].s != b[i].s || a[i].info != b[i].info)
            return i;
    }
    return -1;
}

/* Compare two drivers' outputs (count, then every SMEM's fields, `l` skipped
 * only where l_may_be_stale(a, b, l_hint)); prints and returns false on the
 * first difference. */
static bool check_same(const char *case_name,
                       const char *name_a, const SMEM *a_out, int64_t n_a,
                       const char *name_b, const SMEM *b_out, int64_t n_b,
                       const BsdFirstWalkHint *l_hint)
{
    if (n_a != n_b) {
        fprintf(stderr, "[FAIL] %s: numSeed mismatch: %s=%lld %s=%lld\n",
                case_name, name_a, (long long)n_a, name_b, (long long)n_b);
        return false;
    }
    const int64_t mism = smem_array_first_mismatch(a_out, b_out, n_a, l_hint);
    if (mism >= 0) {
        const SMEM &a = a_out[mism];
        const SMEM &b = b_out[mism];
        fprintf(stderr, "[FAIL] %s: SMEM field mismatch at index %lld:\n"
                "  %-9s rid=%d m=%u n=%u k=%lld l=%lld s=%lld info=%#llx\n"
                "  %-9s rid=%d m=%u n=%u k=%lld l=%lld s=%lld info=%#llx\n",
                case_name, (long long)mism,
                name_a, a.rid, a.m, a.n, (long long)a.k, (long long)a.l, (long long)a.s,
                (unsigned long long)a.info,
                name_b, b.rid, b.m, b.n, (long long)b.k, (long long)b.l, (long long)b.s,
                (unsigned long long)b.info);
        return false;
    }
    return true;
}

static int32_t total_cases = 0;
static int32_t passed_cases = 0;

/* A fixed width a case pins, clamped to the build's BWTSEED_LOCKSTEP_N_MAX so a
 * -D build with a small cap still runs every case instead of tripping the
 * driver's width assert. */
static int32_t clamp_width(int32_t w)
{
    return w < BWTSEED_LOCKSTEP_N_MAX ? w : BWTSEED_LOCKSTEP_N_MAX;
}

static bool run_case(FMI_search *fmi,
                     const char *case_name,
                     int32_t numReads,
                     int32_t max_readlength,
                     int32_t minSeedLen,
                     int32_t max_intv,
                     uint8_t *enc_qdb,
                     const bseq1_t *seq_,
                     int32_t *query_cum_len_ar,
                     int32_t expect_min_seeds = 0,
                     bool expect_seed_per_read = false,
                     bool expect_end_seed = false,
                     const BsdFirstWalkHint *bsd_hint = nullptr)
{
    /* Scalar emits at most one SMEM per outer-x ∈ [0, l_seq), so
     * numReads * max_readlength is a tight upper bound for both paths. */
    const int32_t max_out = numReads * max_readlength;
    const size_t out_bytes = (size_t)(max_out > 0 ? max_out : 1) * sizeof(SMEM);
    SMEM *scalar_out   = (SMEM *)_mm_malloc(out_bytes, 64);
    SMEM *lockstep_out = (SMEM *)_mm_malloc(out_bytes, 64);
    SMEM *ref_out      = (SMEM *)_mm_malloc(out_bytes, 64);

    /* max_intv_array per-read; mutated copies for each path so a buggy
     * implementation that writes to its input array would still parity-check. */
    int32_t *mia_s = (int32_t *)malloc(numReads * sizeof(int32_t));
    int32_t *mia_l = (int32_t *)malloc(numReads * sizeof(int32_t));
    int32_t *mia_r = (int32_t *)malloc(numReads * sizeof(int32_t));
    if (!scalar_out || !lockstep_out || !ref_out || !mia_s || !mia_l || !mia_r) {
        fprintf(stderr, "[FAIL] %s: allocation failed\n", case_name);
        _mm_free(scalar_out);
        _mm_free(lockstep_out);
        _mm_free(ref_out);
        free(mia_s); free(mia_l); free(mia_r);
        total_cases++;
        return false;
    }
    for (int i = 0; i < numReads; i++) { mia_s[i] = max_intv; mia_l[i] = max_intv; mia_r[i] = max_intv; }

    int64_t n_scalar = fmi->bwtSeedStrategyAllPosOneThread(
            enc_qdb, mia_s, numReads, seq_, query_cum_len_ar,
            minSeedLen, scalar_out);
    int64_t n_lockstep = fmi->bwtSeedStrategyAllPosOneThread_lockstep(
            enc_qdb, mia_l, numReads, seq_, query_cum_len_ar,
            minSeedLen, lockstep_out, max_readlength, bsd_hint);
    const int64_t n_ref = reference_bwtseed(fmi, enc_qdb, mia_r, numReads, seq_,
                                            query_cum_len_ar, minSeedLen, ref_out);

    total_cases++;
    /* A lockstep seed taken from a BsdFirstWalk record (always the read's
     * first walk, m = 0) carries the SMEM walk's l, which is stale when that
     * walk ended in its text compare (ls_s1_compare, a record with s = 1). No
     * consumer of a seed reads l (SA resolution uses k and s, the sorts key on
     * rid, m, n, the opt-in dedup compares rid, m, n, k, s), so
     * (rid, m, n, k, s) is the output contract for that seed, and l is skipped
     * there (l_may_be_stale). The lockstep driver uses the hint only when its
     * max_intv and min_len match this call's; a rejected hint skips nothing.
     * The scalar driver takes no hint, so reference-vs-scalar compares l in
     * full. */
    const BsdFirstWalkHint *l_hint =
        bsd_hint != nullptr && bsd_hint->max_intv == max_intv && bsd_hint->min_len == minSeedLen
            ? bsd_hint : nullptr;
    /* All three comparisons run, so a failure names the driver that left the
     * reference rather than only reporting that the two drivers differ. */
    const bool same_sl = check_same(case_name, "scalar", scalar_out, n_scalar,
                                    "lockstep", lockstep_out, n_lockstep, l_hint);
    const bool same_rs = check_same(case_name, "reference", ref_out, n_ref,
                                    "scalar", scalar_out, n_scalar, /*l_hint=*/nullptr);
    const bool same_rl = check_same(case_name, "reference", ref_out, n_ref,
                                    "lockstep", lockstep_out, n_lockstep, l_hint);
    bool ok = same_sl && same_rs && same_rl;

    /* Emission-path cases must actually emit seeds. Byte-identity between two
     * paths that both return zero seeds proves nothing about the buffered
     * emit / ordered-flush logic — it would pass even against an index lacking
     * the expected substrings. Require a nonzero floor so such a run fails
     * loudly instead of masquerading as a parity PASS. */
    if (ok && n_scalar < expect_min_seeds) {
        fprintf(stderr, "[FAIL] %s: emission path not exercised: got %lld seeds, "
                "expected >= %d (is the index the checked-in phiX fixture?)\n",
                case_name, (long long)n_scalar, expect_min_seeds);
        ok = false;
    }

    /* A total can be met while some read emits nothing; where every read is an
     * exact phiX substring well above minSeedLen, require each one to emit. */
    if (ok && expect_seed_per_read) {
        bool *seen = (bool *)calloc(numReads, sizeof(bool));
        if (seen == NULL) {
            fprintf(stderr, "[FAIL] %s: allocation failed\n", case_name);
            ok = false;
        } else {
            for (int64_t j = 0; j < n_scalar; j++) {
                const int32_t r = scalar_out[j].rid;
                if (r >= 0 && r < numReads) seen[r] = true;
            }
            for (int32_t r = 0; r < numReads; r++) {
                if (!seen[r]) {
                    fprintf(stderr, "[FAIL] %s: read %d emitted no seed\n", case_name, r);
                    ok = false;
                    break;
                }
            }
            free(seen);
        }
    }

    /* A case built to reach the bound: some read's last walk that can emit,
     * the one from l_seq - max(minSeedLen, 2), must emit at the read's last
     * base. Without that seed the x_end boundary is not exercised. */
    if (ok && expect_end_seed) {
        const int32_t span = minSeedLen > 2 ? minSeedLen : 2;
        bool found = false;
        for (int64_t j = 0; j < n_scalar && !found; j++) {
            const SMEM &a = scalar_out[j];
            const int32_t l = seq_[a.rid].l_seq;
            found = (int32_t)a.n == l - 1 && (int32_t)a.m == l - span;
        }
        if (!found) {
            fprintf(stderr, "[FAIL] %s: no read emitted a seed from its last emitting "
                    "start (l_seq - %d) to its last base\n", case_name, (int)span);
            ok = false;
        }
    }

    if (ok) {
        fprintf(stderr, "[PASS] %s  (n_seed=%lld)\n", case_name, (long long)n_scalar);
        passed_cases++;
    }

    _mm_free(scalar_out);
    _mm_free(lockstep_out);
    _mm_free(ref_out);
    free(mia_s); free(mia_l); free(mia_r);
    return ok;
}

int main(int argc, char *argv[]) {
    if (argc != 2) {
        fprintf(stderr, "Usage: %s <bwa-mem3 index prefix>\n", argv[0]);
        return 2;
    }

    fprintf(stderr, "bwtseed_lockstep_parity_test: BWTSEED_LOCKSTEP_N=%d\n",
            (int)BWTSEED_LOCKSTEP_N);

    FMI_search *fmi = new FMI_search(argv[1]);
    fmi->load_index();

    /* The default max_intv used in bwamem.cpp third-pass is opt->max_mem_intv
     * (default 20). minSeedLen passed to bwtSeed is opt->min_seed_len + 1 = 20.
     * Cases below use those defaults except where varying them exercises a
     * specific path. */
    const int32_t DEFAULT_MAX_INTV   = 20;
    const int32_t DEFAULT_MIN_SEED   = 20;

    /* Case 1: two simple reads. Smoke-test the basic flow. */
    {
        const char *reads[] = {
            "ACGTACGTACGTACGTACGTACGTACGTACGT",
            "TGCATGCATGCATGCATGCATGCATGCATGCA",
        };
        int32_t numReads = 2;
        int32_t max_readlength = 32;
        uint8_t *enc_qdb; bseq1_t *seq_; int32_t *cum_len;
        encode_reads(reads, numReads, &enc_qdb, &seq_, &cum_len);
        run_case(fmi, "Case 1: two simple reads",
                 numReads, max_readlength, DEFAULT_MIN_SEED, DEFAULT_MAX_INTV,
                 enc_qdb, seq_, cum_len);
        free(enc_qdb); free(seq_); free(cum_len);
    }

    /* Case 2: numReads < N — exercises the "unused slots" init path. */
    {
        const char *reads[] = {
            "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT",
            "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT",
        };
        int32_t numReads = 2;
        int32_t max_readlength = 40;
        uint8_t *enc_qdb; bseq1_t *seq_; int32_t *cum_len;
        encode_reads(reads, numReads, &enc_qdb, &seq_, &cum_len);
        run_case(fmi, "Case 2: numReads < N",
                 numReads, max_readlength, DEFAULT_MIN_SEED, DEFAULT_MAX_INTV,
                 enc_qdb, seq_, cum_len);
        free(enc_qdb); free(seq_); free(cum_len);
    }

    /* Case 3: small batch (numReads < N at the default width 24) — every read
     * seeds its own slot, no recycle. Varied read contents. */
    {
        const char *reads[] = {
            "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT",
            "TGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCA",
            "AAAACCCCGGGGTTTTAAAACCCCGGGGTTTTAAAACCCCGGGGTTTTAAAACCCCGGGGTTTT",
            "GATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATC",
        };
        int32_t numReads = 4;
        int32_t max_readlength = 64;
        uint8_t *enc_qdb; bseq1_t *seq_; int32_t *cum_len;
        encode_reads(reads, numReads, &enc_qdb, &seq_, &cum_len);
        run_case(fmi, "Case 3: small batch, no recycle (varied contents)",
                 numReads, max_readlength, DEFAULT_MIN_SEED, DEFAULT_MAX_INTV,
                 enc_qdb, seq_, cum_len);
        free(enc_qdb); free(seq_); free(cum_len);
    }

    /* Case 4: a larger small batch (still < N at the default width 24) — distinct
     * read contents, no recycle here. The recycle path is exercised by Case 11r
     * (width 8) and the default-width recycle case (Case 11s) below. */
    {
        const char *reads[] = {
            "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT",
            "TGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCA",
            "AAAACCCCGGGGTTTTAAAACCCCGGGGTTTTAAAACCCCGGGGTTTTAAAACCCCGGGGTTTT",
            "GATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATC",
            "CAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGT",
            "TTTACCCAGGGCTTTACCCAGGGCTTTACCCAGGGCTTTACCCAGGGCTTTACCCAGGGCTTTA",
            "ATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATAT",
        };
        int32_t numReads = 7;
        int32_t max_readlength = 64;
        uint8_t *enc_qdb; bseq1_t *seq_; int32_t *cum_len;
        encode_reads(reads, numReads, &enc_qdb, &seq_, &cum_len);
        run_case(fmi, "Case 4: small batch, distinct contents (no recycle)",
                 numReads, max_readlength, DEFAULT_MIN_SEED, DEFAULT_MAX_INTV,
                 enc_qdb, seq_, cum_len);
        free(enc_qdb); free(seq_); free(cum_len);
    }

    /* Case 5: first base non-ACGT, also fully-N reads mixed in. Tests the
     * bsd_init_slot seek-past-N path and the all-N → BSD_DONE-with-zero-emits
     * path, with valid neighbors so the flush cursor must still advance. */
    {
        const char *reads[] = {
            "NACGTACGTACGTACGTACGTACGTACGTACGT",  // skip first N then walk
            "TGCATGCATGCATGCATGCATGCATGCATGCA",
            "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN",  // all-N; zero emits, must not stall flush
            "GATCGATCGATCGATCGATCGATCGATCGATC",
        };
        int32_t numReads = 4;
        /* read[0] is 33 bases (leading N + 32), the longest here — the bound
         * passed to the lockstep and the max_out sizing must cover it. */
        int32_t max_readlength = 33;
        uint8_t *enc_qdb; bseq1_t *seq_; int32_t *cum_len;
        encode_reads(reads, numReads, &enc_qdb, &seq_, &cum_len);
        run_case(fmi, "Case 5: non-ACGT cases",
                 numReads, max_readlength, DEFAULT_MIN_SEED, DEFAULT_MAX_INTV,
                 enc_qdb, seq_, cum_len);
        free(enc_qdb); free(seq_); free(cum_len);
    }

    /* Case 6: mid-read N — exercises outer-x re-seed after non-ACGT in inner-j.
     * Scalar: when enc_qdb[offset+j] >= 4 inside inner-j, next_x = j+1; break.
     * Lockstep: must produce identical SMEMs for outer-x positions on either
     * side of the N. */
    {
        const char *reads[] = {
            "ACGTACGTACGTACGTNACGTACGTACGTACG",  // N at pos 16
            "TGCATGCATGCNTGCATGCATGCATGCATGCA",  // N at pos 11
        };
        int32_t numReads = 2;
        int32_t max_readlength = 32;
        uint8_t *enc_qdb; bseq1_t *seq_; int32_t *cum_len;
        encode_reads(reads, numReads, &enc_qdb, &seq_, &cum_len);
        run_case(fmi, "Case 6: mid-read N",
                 numReads, max_readlength, DEFAULT_MIN_SEED, DEFAULT_MAX_INTV,
                 enc_qdb, seq_, cum_len);
        free(enc_qdb); free(seq_); free(cum_len);
    }

    /* Case 7: real 300bp Illumina WGS read pair (SRR6109255.100035 R1/R2).
     * Mirrors the smem_lockstep_parity_test fixture so we cover the same
     * adapter contamination + low-complexity tail patterns in the bwtSeed
     * regime. */
    {
        const char *r1 =
            "ATGACCTCCCTAATATCTTCAGAATAGATCGGAAGAGCACACGTCTGAACTCCAGTCACT"
            "ACAAGATCTCGTATGCCGTCTTCTGCTTGAAAAAAAAAAAAAAAAAATGAGTAACTCTAT"
            "GAATCGAACTCAAACATTCCACACACGTTTATGCCACTACTTCTATCGTGCTTTCATATT"
            "ATACTACTCCCCCTTCCCCCTCCATCTTCACCTCTCCTCCAATCTCCACTCACCTCTACC"
            "TCAACACGTACTTTACACACATCTCTCCCCCACGGACCACAATACCTCTCCTCTCAATTA";
        const char *r2 =
            "ATTCTGAAGATATTAGGGAGGTCATAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAG"
            "ATCTCGGTGGTCGCCGTATCATTAAAAAAAAAAAAAATAAAAACAGAAGCAGTGTCAGTA"
            "GCATAGTGATGAATATAGCAATAAACGCACAGCTTGAAACCTACCGTTTGCGACAGCATC"
            "TCTACAGACGTCTGTTTATCTACTTTGAAAGTGGCTACGTGGAACACTCATAGTCATACC"
            "ACTAATCAATATCATGAAATTACCAGGTTAGTCTGTATACTACGATAAAGACACAGACTT";
        const char *reads[] = { r1, r2 };
        int32_t numReads = 2;
        int32_t max_readlength = 300;
        uint8_t *enc_qdb; bseq1_t *seq_; int32_t *cum_len;
        encode_reads(reads, numReads, &enc_qdb, &seq_, &cum_len);
        run_case(fmi, "Case 7: real 300bp SRR6109255.100035",
                 numReads, max_readlength, DEFAULT_MIN_SEED, DEFAULT_MAX_INTV,
                 enc_qdb, seq_, cum_len);
        free(enc_qdb); free(seq_); free(cum_len);
    }

    /* Case 8: smaller max_intv = 1 — every position must produce a seed
     * (forces the cond `s < max_intv` ≡ `s < 1` ≡ `s == 0`, which excludes
     * the emit because of the `s > 0` guard). Tests the no-emit branch
     * inside the cond. */
    {
        const char *reads[] = {
            "ACGTACGTACGTACGTACGTACGTACGTACGT",
            "TGCATGCATGCATGCATGCATGCATGCATGCA",
        };
        int32_t numReads = 2;
        int32_t max_readlength = 32;
        uint8_t *enc_qdb; bseq1_t *seq_; int32_t *cum_len;
        encode_reads(reads, numReads, &enc_qdb, &seq_, &cum_len);
        run_case(fmi, "Case 8: max_intv = 1",
                 numReads, max_readlength, DEFAULT_MIN_SEED, /*max_intv=*/1,
                 enc_qdb, seq_, cum_len);
        free(enc_qdb); free(seq_); free(cum_len);
    }

    /* Case 9: minSeedLen larger than any plausible match — exercises the
     * (n - m + 1) >= minSeedLen guard. Scalar would break inner-j but skip
     * the emit; lockstep must do the same. */
    {
        const char *reads[] = {
            "ACGTACGTACGTACGTACGTACGTACGTACGT",
            "TGCATGCATGCATGCATGCATGCATGCATGCA",
        };
        int32_t numReads = 2;
        int32_t max_readlength = 32;
        uint8_t *enc_qdb; bseq1_t *seq_; int32_t *cum_len;
        encode_reads(reads, numReads, &enc_qdb, &seq_, &cum_len);
        run_case(fmi, "Case 9: minSeedLen = 1000 (always-fail guard)",
                 numReads, max_readlength, /*minSeed=*/1000, DEFAULT_MAX_INTV,
                 enc_qdb, seq_, cum_len);
        free(enc_qdb); free(seq_); free(cum_len);
    }

    /* Case 10p: phiX subsequence with parameters chosen to actually emit
     * seeds. Cases 1-9 use min_seed_len=20 + max_intv=20 against synthetic
     * reads that don't match phiX, so the inner-j loop collapses to s=0
     * before any emit. To exercise the emission path we feed back the
     * first ~200 bp of phiX itself (genuine substring matches in the
     * index) with min_seed_len=10 and max_intv=5: walks terminate at
     * small-but-nonzero s, satisfying both `s < max_intv` and the
     * `(n - m + 1) >= minSeedLen` guards. This is the case that actually
     * proves parity on non-empty match_buf contents. */
    {
        const char *r1 =
            "GAGTTTTATCGCTTCCATGACGCAGAAGTTAACACTTTCGGATATTTCTGATGAGTCGAAAA";
        const char *r2 =
            "ATTATCTTGATAAAGCAGGAATTACTACTGCTTGTTTACGAATTAAATCGAAGTGGACTGCT";
        const char *r3 =
            "GGCGGAAAATGAGAAAATTCGACCTATCCTTGCGCAGCTCGAGAAGCTCTTACTTTGCGACC";
        const char *r4 =
            "TTTCGCCATCAACTAACGATTCTGTCAAAAACTGACGCGTTGGATGAGGAGAAGTGGCTTAA";
        const char *reads[] = { r1, r2, r3, r4 };
        int32_t numReads = 4;
        int32_t max_readlength = 64;
        uint8_t *enc_qdb; bseq1_t *seq_; int32_t *cum_len;
        encode_reads(reads, numReads, &enc_qdb, &seq_, &cum_len);
        run_case(fmi, "Case 10p: phiX subsequence (real matches)",
                 numReads, max_readlength, /*minSeed=*/10, /*max_intv=*/5,
                 enc_qdb, seq_, cum_len, /*expect_min_seeds=*/10);
        free(enc_qdb); free(seq_); free(cum_len);
    }

    /* Case 11p: phiX subsequence, many reads, exercises the
     * stepping-pass/flush interleave when emitted matches are non-empty.
     * Same parameters as 10p, more reads (still < N at the default width 24, so
     * no recycle) with non-trivial match_buf contents. */
    {
        /* 8 overlapping windows over the phiX start; each will emit
         * several seeds at min_seed_len=10 max_intv=5. */
        const char *r0 = "GAGTTTTATCGCTTCCATGACGCAGAAGTTAACACTTTCGGATA";
        const char *r1 = "TTTATCGCTTCCATGACGCAGAAGTTAACACTTTCGGATATTTC";
        const char *r2 = "TCGCTTCCATGACGCAGAAGTTAACACTTTCGGATATTTCTGAT";
        const char *r3 = "CTTCCATGACGCAGAAGTTAACACTTTCGGATATTTCTGATGAG";
        const char *r4 = "CATGACGCAGAAGTTAACACTTTCGGATATTTCTGATGAGTCGA";
        const char *r5 = "ACGCAGAAGTTAACACTTTCGGATATTTCTGATGAGTCGAAAAA";
        const char *r6 = "AGAAGTTAACACTTTCGGATATTTCTGATGAGTCGAAAAATTAT";
        const char *r7 = "TTAACACTTTCGGATATTTCTGATGAGTCGAAAAATTATCTTGA";
        const char *reads[] = { r0, r1, r2, r3, r4, r5, r6, r7 };
        int32_t numReads = 8;
        int32_t max_readlength = 44;
        uint8_t *enc_qdb; bseq1_t *seq_; int32_t *cum_len;
        encode_reads(reads, numReads, &enc_qdb, &seq_, &cum_len);
        run_case(fmi, "Case 11p: phiX subsequence, many reads",
                 numReads, max_readlength, /*minSeed=*/10, /*max_intv=*/5,
                 enc_qdb, seq_, cum_len, /*expect_min_seeds=*/10);
        free(enc_qdb); free(seq_); free(cum_len);
    }

    /* Case 11r: numReads > N forces slot recycle on non-empty match_buf
     * contents. The width is runtime (g_bwtseed_lockstep_n) and its default (24)
     * exceeds 20, so pin it to 8 here: 20 reads then guarantees
     * the first 8 slots retire while 12 inputs are still pending. The recycle path
     * must reset match_count, ready, and phase for each recycled slot without
     * corrupting in-flight slots' state. Byte-identity is width-invariant, so the
     * pin only sharpens coverage, never the pass/fail. */
    {
        /* 20 sliding 44-bp windows over the phiX start. */
        const char *reads[20];
        const char *src =
            "GAGTTTTATCGCTTCCATGACGCAGAAGTTAACACTTTCGGATATTTCTGATGAGTCGAAAA"
            "ATTATCTTGATAAAGCAGGAATTACTACTGCTTGTTTACGAATTAAATCGAAGTGGACTGCT";
        char buf[20][45];
        for (int i = 0; i < 20; i++) {
            memcpy(buf[i], src + i, 44);
            buf[i][44] = '\0';
            reads[i] = buf[i];
        }
        int32_t numReads = 20;
        int32_t max_readlength = 44;
        uint8_t *enc_qdb; bseq1_t *seq_; int32_t *cum_len;
        encode_reads(reads, numReads, &enc_qdb, &seq_, &cum_len);
        const int32_t saved_n = g_bwtseed_lockstep_n;
        g_bwtseed_lockstep_n = clamp_width(8);   /* force numReads(20) > N to hit the recycle path */
        run_case(fmi, "Case 11r: numReads > N, slot recycle on non-empty match_buf",
                 numReads, max_readlength, /*minSeed=*/10, /*max_intv=*/5,
                 enc_qdb, seq_, cum_len, /*expect_min_seeds=*/20,
                 /*expect_seed_per_read=*/true);
        g_bwtseed_lockstep_n = saved_n;
        free(enc_qdb); free(seq_); free(cum_len);
    }

    /* Case 11s: same slot-recycle path as 11r, but at the compile-time default
     * width (no pin), so numReads must exceed the default. 30 reads > 24 forces
     * the recycle branch at the width non-Zen 5 hosts run. This is the coverage
     * 11r's width-8 pin does not give (the recycle path at N = default); Case 11z
     * covers the AMD Zen 5 host default. */
    {
        /* 30 sliding 44-bp windows over the phiX start (src is 122 bp, so a
         * window at offset i needs i + 44 <= 122, i.e. up to i = 78 > 29). */
        const char *reads[30];
        const char *src =
            "GAGTTTTATCGCTTCCATGACGCAGAAGTTAACACTTTCGGATATTTCTGATGAGTCGAAAA"
            "ATTATCTTGATAAAGCAGGAATTACTACTGCTTGTTTACGAATTAAATCGAAGTGGACTGCT";
        char buf[30][45];
        for (int i = 0; i < 30; i++) {
            memcpy(buf[i], src + i, 44);
            buf[i][44] = '\0';
            reads[i] = buf[i];
        }
        int32_t numReads = 30;
        int32_t max_readlength = 44;
        uint8_t *enc_qdb; bseq1_t *seq_; int32_t *cum_len;
        encode_reads(reads, numReads, &enc_qdb, &seq_, &cum_len);
        if (g_bwtseed_lockstep_n >= numReads)   /* only a -D build raising the default */
            fprintf(stderr, "[NOTE] Case 11s: default width %d >= %d reads, so no slot "
                    "recycles here (Case 11w covers it)\n", (int)g_bwtseed_lockstep_n,
                    (int)numReads);
        run_case(fmi, "Case 11s: numReads > default N, slot recycle at the compile-time default width",
                 numReads, max_readlength, /*minSeed=*/10, /*max_intv=*/5,
                 enc_qdb, seq_, cum_len, /*expect_min_seeds=*/20,
                 /*expect_seed_per_read=*/true);
        free(enc_qdb); free(seq_); free(cum_len);
    }

    /* Case 11z: slot recycle at the AMD Zen 5 host-default width (32 unless a -D
     * override clamps it), the width Zen 5 hosts run by default. 70 sliding
     * 44-bp windows (offsets 0..69 of the 122-bp src) exceed that width twice
     * over, so every slot recycles at least once with non-empty match_buf
     * contents. The width is pinned through the runtime global and restored,
     * as in Case 11r. */
    {
        const int32_t zen5_n = bwa3_bwtseed_host_default_width(1, 0x1a, 0x02);
        enum { NZ = 70 };
        const char *reads[NZ];
        const char *src =
            "GAGTTTTATCGCTTCCATGACGCAGAAGTTAACACTTTCGGATATTTCTGATGAGTCGAAAA"
            "ATTATCTTGATAAAGCAGGAATTACTACTGCTTGTTTACGAATTAAATCGAAGTGGACTGCT";
        char buf[NZ][45];
        for (int i = 0; i < NZ; i++) {
            memcpy(buf[i], src + i, 44);
            buf[i][44] = '\0';
            reads[i] = buf[i];
        }
        if (zen5_n >= NZ) {
            /* Only a -D build raising the widths gets here; Case 11w covers
             * slot recycle at every width up to BWTSEED_LOCKSTEP_N_MAX. */
            fprintf(stderr, "[SKIP] Case 11z: host-default width %d leaves no recycle "
                    "with %d reads (Case 11w covers it)\n", zen5_n, (int)NZ);
        } else {
            int32_t numReads = NZ;
            int32_t max_readlength = 44;
            uint8_t *enc_qdb; bseq1_t *seq_; int32_t *cum_len;
            encode_reads(reads, numReads, &enc_qdb, &seq_, &cum_len);
            const int32_t saved_n = g_bwtseed_lockstep_n;
            g_bwtseed_lockstep_n = zen5_n;
            char name[128];
            snprintf(name, sizeof(name),
                     "Case 11z: numReads(%d) > N(%d), slot recycle at the Zen 5 host-default width",
                     (int)NZ, (int)zen5_n);
            run_case(fmi, name,
                     numReads, max_readlength, /*minSeed=*/10, /*max_intv=*/5,
                     enc_qdb, seq_, cum_len, /*expect_min_seeds=*/NZ,
                     /*expect_seed_per_read=*/true);
            g_bwtseed_lockstep_n = saved_n;
            free(enc_qdb); free(seq_); free(cum_len);
        }
    }

    /* Case 11w: many more reads than slots, at several widths. The flush pass
     * finds the slot holding the next read to emit as input mod N, which relies
     * on slot recycling keeping input i in slot i mod N. 2 * N_MAX + 22 windows
     * (150 at the default N_MAX of 64) of every length from 24 to 44 over the
     * phiX start (so reads retire out of input order, and every slot is reused
     * at least once even at the maximum width) at a small width, the shipped
     * default and the maximum. */
    {
        /* EXTRA keeps NR at 150 for the default N_MAX of 64. */
        enum { EXTRA = 22, NR = 2 * BWTSEED_LOCKSTEP_N_MAX + EXTRA, MAXL = 44 };
        const char *src =
            "GAGTTTTATCGCTTCCATGACGCAGAAGTTAACACTTTCGGATATTTCTGATGAGTCGAAAA"
            "ATTATCTTGATAAAGCAGGAATTACTACTGCTTGTTTACGAATTAAATCGAAGTGGACTGCT";
        static char buf[NR][MAXL + 1];
        const char *reads[NR];
        for (int i = 0; i < NR; i++) {
            /* Both multipliers are coprime to their moduli, so every length and
             * offset occurs and the windows do not repeat with a short period. */
            const int len = 24 + (i * 5) % (MAXL - 23);   /* 24..44 */
            const int off = (i * 7) % (122 - MAXL);        /* 0..77: window stays in src */
            memcpy(buf[i], src + off, len);
            buf[i][len] = '\0';
            reads[i] = buf[i];
        }
        uint8_t *enc_qdb; bseq1_t *seq_; int32_t *cum_len;
        encode_reads(reads, NR, &enc_qdb, &seq_, &cum_len);
        const int32_t saved_n = g_bwtseed_lockstep_n;
        const int32_t widths[] = {clamp_width(3), BWTSEED_LOCKSTEP_N, BWTSEED_LOCKSTEP_N_MAX};
        for (const int32_t w : widths) {
            char name[96];
            snprintf(name, sizeof name, "Case 11w: %d reads > N = %d, slot recycle", (int)NR, (int)w);
            g_bwtseed_lockstep_n = w;
            run_case(fmi, name, NR, MAXL, /*minSeed=*/10, /*max_intv=*/5,
                     enc_qdb, seq_, cum_len, /*expect_min_seeds=*/NR,
                     /*expect_seed_per_read=*/true);
        }
        g_bwtseed_lockstep_n = saved_n;
        free(enc_qdb); free(seq_); free(cum_len);
    }

    /* Case 12p: phiX subsequence with mid-read N — combines real matches
     * (so seeds get emitted) with the non-ACGT outer-x advancement path
     * (so the seek-past-N branch in bsd_advance_step gets exercised on
     * non-empty match_buf state). */
    {
        const char *r1 = "GAGTTTTATCGCTTCCANGACGCAGAAGTTAACACTTTCGGATA";
        const char *r2 = "TTTATCGCTTCCATGACGCAGAANTTAACACTTTCGGATATTTC";
        const char *reads[] = { r1, r2 };
        int32_t numReads = 2;
        int32_t max_readlength = 44;
        uint8_t *enc_qdb; bseq1_t *seq_; int32_t *cum_len;
        encode_reads(reads, numReads, &enc_qdb, &seq_, &cum_len);
        run_case(fmi, "Case 12p: phiX with mid-read N (real matches)",
                 numReads, max_readlength, /*minSeed=*/10, /*max_intv=*/5,
                 enc_qdb, seq_, cum_len, /*expect_min_seeds=*/3);
        free(enc_qdb); free(seq_); free(cum_len);
    }

    /* Case 13j: reads shorter than the seed jump depth, and minSeedLen around
     * it. A walk jumps to length min(depth, minSeedLen - 1) only when that many
     * ACGT bases follow its start, and not at all below length 2 (minSeedLen
     * 2), so these reads take the jump, stop short of it at the read end or an
     * N, or skip it, at every minSeedLen here. */
    {
        const std::vector<std::string> rs = jump_depth_reads();
        std::vector<const char *> reads;
        for (const std::string &r : rs) reads.push_back(r.c_str());
        const int32_t numReads = (int32_t)reads.size();
        uint8_t *enc_qdb; bseq1_t *seq_; int32_t *cum_len;
        encode_reads(reads.data(), numReads, &enc_qdb, &seq_, &cum_len);
        const int32_t min_seeds[] = {2, 3, 4, 5, 11, 12, 13};
        for (const int32_t msl : min_seeds) {
            char name[96];
            snprintf(name, sizeof name, "Case 13j: reads of 1-14 bases, minSeedLen %d", (int)msl);
            run_case(fmi, name, numReads, /*max_readlength=*/14, msl, /*max_intv=*/50,
                     enc_qdb, seq_, cum_len, /*expect_min_seeds=*/1);
        }
        free(enc_qdb); free(seq_); free(cum_len);
    }

    /* Case 13r: the read-end bound. Both drivers end a read at
     * x_end = l_seq - minSeedLen + 1 (the read end for minSeedLen <= 1), where
     * no later walk can emit. Reads of every length from 0 to
     * 3 * max(minSeedLen, 2) + 2 over phiX, and the same reads with an N just
     * before, at and after the last start that can emit and at the last base,
     * put a walk's start, its N or the read end on each side of the bound. At
     * max_intv 1000 every walk long enough emits at its first chance, so the
     * read whose length is a multiple of the seed length emits from its last
     * emitting start to its last base (required). At max_intv 5 walks also run
     * off the read end for minSeedLen 0-3; at 10 and 20 every phiX substring of
     * that length is already under 5, so those runs match max_intv 1000. At
     * max_intv 1 nothing emits (exact phiX substrings never reach s == 0), so
     * every walk from below x_end, through its jump at the larger minSeedLen,
     * runs to an N or the read end and the read retires there. minSeedLen -1
     * emits nothing (the length test is unsigned), 0 and 1 take the read end
     * as the bound, the rest x_end. Each run is checked against the unbounded
     * reference walk. */
    {
        const std::string src =
            "GAGTTTTATCGCTTCCATGACGCAGAAGTTAACACTTTCGGATATTTCTGATGAGTCGAAAA"
            "ATTATCTTGATAAAGCAGGAATTACTACTGCTTGTTTACGAATTAAATCGAAGTGGACTGCT";
        const int32_t min_seeds[] = {-1, 0, 1, 2, 3, 10, 20};
        const int32_t max_intvs[] = {1000, 5, 1};
        for (const int32_t msl : min_seeds) {
            const int32_t span = msl > 2 ? msl : 2;
            const int32_t max_len = 3 * span + 2;
            std::vector<std::string> rs;
            for (int32_t L = 0; L <= max_len; L++) {
                const std::string r = src.substr(3, (size_t)L);
                rs.push_back(r);
                const int32_t ns[] = {L - span - 1, L - span, L - span + 1, L - 1};
                for (const int32_t p : ns) {
                    if (p < 0 || p >= L) continue;
                    std::string rn = r;
                    rn[(size_t)p] = 'N';
                    rs.push_back(rn);
                }
            }
            std::vector<const char *> reads;
            for (const std::string &r : rs) reads.push_back(r.c_str());
            const int32_t numReads = (int32_t)reads.size();
            uint8_t *enc_qdb; bseq1_t *seq_; int32_t *cum_len;
            encode_reads(reads.data(), numReads, &enc_qdb, &seq_, &cum_len);
            for (const int32_t mi : max_intvs) {
                char name[112];
                snprintf(name, sizeof name, "Case 13r: read-end bound, reads of 0-%d bases, "
                         "minSeedLen %d, max_intv %d", (int)max_len, (int)msl, (int)mi);
                run_case(fmi, name, numReads, max_len, msl, mi, enc_qdb, seq_, cum_len,
                         /*expect_min_seeds=*/msl < 0 || mi == 1 ? 0 : 1,
                         /*expect_seed_per_read=*/false,
                         /*expect_end_seed=*/msl >= 0 && mi == 1000);
            }
            free(enc_qdb); free(seq_); free(cum_len);
        }
    }

    /* Case 10: long reads (1000bp, 1500bp) — exercises the heap-allocated
     * match_buf scaling on the lockstep path. The scalar writes directly
     * into matchArray; lockstep buffers per slot then flushes. */
    {
        char r_1000bp[1001];
        char r_1500bp[1501];
        const char *bases = "ACGT";
        for (int i = 0; i < 1000; i++) r_1000bp[i] = bases[(i * 7 + 3) & 3];
        r_1000bp[1000] = '\0';
        for (int i = 0; i < 1500; i++) r_1500bp[i] = bases[(i * 11 + 5) & 3];
        r_1500bp[1500] = '\0';
        const char *reads[] = { r_1000bp, r_1500bp };
        int32_t numReads = 2;
        int32_t max_readlength = 1500;
        uint8_t *enc_qdb; bseq1_t *seq_; int32_t *cum_len;
        encode_reads(reads, numReads, &enc_qdb, &seq_, &cum_len);
        run_case(fmi, "Case 10: long reads (1000bp, 1500bp)",
                 numReads, max_readlength, DEFAULT_MIN_SEED, DEFAULT_MAX_INTV,
                 enc_qdb, seq_, cum_len);
        free(enc_qdb); free(seq_); free(cum_len);
    }

    /* Case 13h: the BsdFirstWalk hand-off. mem_collect_smem runs the SMEM pass
     * (getSMEMsAllPosOneThread, min_intv 1) with a BsdFirstWalkHint, whose first
     * forward walk on each read records where the bwtSeed pass's first walk on
     * it stops; the lockstep bwtSeed pass then emits that walk's seed and starts
     * the read there. The scalar bwtSeedStrategyAllPosOneThread walks every read
     * from x = 0, so it is the oracle: the hinted lockstep stream must equal it
     * field by field and in count, with l skipped only for a read's first seed
     * taken from a record with s = 1 (l_may_be_stale).
     *
     * The records start as garbage with x > 0, so a record the SMEM pass did
     * not write (a read starting with N) would be used if
     * getSMEMsAllPosOneThread did not mark every read unrecorded first; every
     * read whose base 0 is ACGT must then be recorded, and no other. The
     * hinted SMEM pass's own output must equal an unhinted run's, every field
     * (info and l included). Reads:
     * phiX windows of 30-151 bases, forward and reverse-complemented (unique
     * matches with many bases left, so the SMEM walk takes its text compare),
     * windows with an N first, early, mid-read and last (the SMEM walk stops
     * before the bwtSeed walk does), a substituted base, reads at and below the
     * seed length, substitutions and N's around base min_len - 1, palindromes
     * whose walk falls to s = 0 at base min_len - 1 without a unique match,
     * windows whose first 11 bases occur twice in the text (a record with
     * s = 2 at seed length 11, whose l is compared), and non-matching
     * synthetic reads (s = 0). Run at max_intv 1 (no seed: the walk stops only
     * at s = 0), 5 and 20, at seed lengths 20 and 11; then with a hint whose
     * min_len or max_intv differs from the call's, which must be ignored (at
     * max_intv 1 with a hint recorded at 20, using it would emit seeds the
     * scalar does not); with the SMEM pass at min_intv 2, whose walk stops at
     * s = 1 and so must record nothing (a record there would end the bwtSeed
     * walk early); and with a hint max_intv of 0, which must record nothing
     * (at max_intv 0 no walk emits, so only the record count can show it).
     * Each matched run that emits seeds also runs a canary that shifts k in one
     * recorded seed at a time and requires that read's output to change, so a
     * bwtSeed pass that ignored any read's record with s > 0 would fail (a
     * record with s = 0 emits no seed). Every run repeats on the index loaded
     * without its pac, where no record comes from the text compare. */
#if SMEM_LOCKSTEP_N > 1
    {
        /* The last base the windows below take: R.substr(2500 + 97 * 27, 60). */
        const int32_t min_ref_len = 2500 + 97 * 27 + 60;
        int32_t ref_len = 0;
        char *ref = read_fasta_seq(argv[1], &ref_len);
        if (!ref || ref_len < min_ref_len || !fmi->s1_text_available()) {
            total_cases++;
            fprintf(stderr, "[FAIL] Case 13h: needs a reference FASTA of at least %d bases "
                    "at '%s' and its pac loaded (s1_text_available), so the SMEM walk's "
                    "text compare records walks too\n", (int)min_ref_len, argv[1]);
        } else {
            const std::string R(ref, (size_t)ref_len);
            std::vector<std::string> rv;
            for (int i = 0; i < 48; i++) {
                const int len = 30 + (i * 13) % 122;                       /* 30..151 */
                const int off = (int)(((int64_t)i * 389) % (ref_len - 160));
                std::string w = R.substr((size_t)off, (size_t)len);
                if (i & 1) w = revcomp(w);
                switch (i % 8) {
                case 2: w[0] = 'N'; break;                                 /* no first walk */
                case 3: w[3] = 'N'; break;                                 /* SMEM walk stops at 3 */
                case 4: w[w.size() / 2] = 'N'; break;
                case 5: w[w.size() - 1] = 'N'; break;
                case 6: w[22] = (w[22] == 'A') ? 'C' : 'A'; break;         /* a mismatch */
                default: break;
                }
                rv.push_back(w);
            }
            for (int len : {11, 12, 19, 20, 21})                           /* around min_len */
                rv.push_back(R.substr(2000, (size_t)len));
            /* A substitution or an N at, just before and just after base
             * min_len - 1 (10 and 19 for the seed lengths below): the SMEM
             * walk's text compare ends on the step at which the bwtSeed walk
             * could first stop, where an s = 0 step or an N ends that walk
             * without a seed. */
            for (int p : {9, 10, 11, 18, 19, 20})
                for (int kind = 0; kind < 2; kind++) {
                    std::string w = R.substr((size_t)(3000 + 61 * p + 7 * kind), 70);
                    w[(size_t)p] = kind ? 'N' : (w[(size_t)p] == 'A' ? 'C' : 'A');
                    rv.push_back(w);
                }
            /* A 10-base reverse palindrome occurs twice in the doubled text (on
             * both strands), so its walk keeps s >= 2 through base 9 and takes no
             * unique-match hand-off; a next base that neither occurrence
             * continues with then gives s = 0 at base 10, min_len - 1 at seed
             * length 11: at jump depth 12 inside the jump-table steps, else in a
             * forward step. */
            int32_t n_pal = 0;
            for (int32_t a = 1; a + 70 < ref_len && n_pal < 2; a++) {
                const std::string w10 = R.substr((size_t)a, 10);
                if (w10.find_first_not_of("ACGT") != std::string::npos || w10 != revcomp(w10))
                    continue;
                const char f = R[(size_t)a + 10];
                const char r = revcomp(R.substr((size_t)a - 1, 1))[0];
                char c = 'A';
                while (c == f || c == r) c = c == 'A' ? 'C' : c == 'C' ? 'G' : 'T';
                rv.push_back(w10 + c + R.substr((size_t)a + 11, 49));
                n_pal++;
            }
            if (n_pal < 2) {
                total_cases++;
                fprintf(stderr, "[FAIL] Case 13h: found %d of 2 10-base reverse palindromes "
                        "in '%s'\n", (int)n_pal, argv[1]);
            }
            rv.push_back("ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT");       /* s = 0 early */
            rv.push_back("TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT");
            /* A non-matching prefix (s = 0 within a few bases) ahead of a phiX
             * window: the SMEM walk stops in the prefix, before the bwtSeed walk
             * can, and the bwtSeed walk's next start (just past its minimum
             * length, or past an N) decides which phiX seed comes next. */
            for (int pl : {9, 13, 19, 20, 21, 27}) {
                std::string w;
                for (int j = 0; j < pl; j++) w += "CGCGATATCG"[j % 10];
                if (pl == 13) w[11] = 'N';
                rv.push_back(w + R.substr((size_t)(2500 + 97 * pl), 60));
            }
            rv.push_back(R.substr(0, 60));                                 /* the text's first bases */
            rv.push_back(revcomp(R.substr(0, 60)) + "ACGTTGCA");           /* runs past the text end */
            /* Windows whose first 11 bases occur at least twice in the doubled
             * text (no 11-mer is its own reverse complement, so these are
             * true repeats): at seed length 11 their first walk has s >= 2 on
             * every step through base 10 and no unique-match hand-off, so the
             * record (s >= 2) comes from a forward step and its l is exact. */
            {
                std::map<std::string, int32_t> kc;
                const std::string Rrc = revcomp(R);
                for (const std::string *t : {&R, &Rrc})
                    for (size_t a = 0; a + 11 <= t->size(); a++) kc[t->substr(a, 11)]++;
                int32_t n_rep = 0;
                for (int32_t a = 0; a + 70 < ref_len && n_rep < 2; a++) {
                    const std::string w11 = R.substr((size_t)a, 11);
                    if (w11.find_first_not_of("ACGT") != std::string::npos || kc[w11] < 2)
                        continue;
                    rv.push_back(R.substr((size_t)a, 60));
                    n_rep++;
                    a += 60;                                               /* a different repeat */
                }
                if (n_rep < 2) {
                    total_cases++;
                    fprintf(stderr, "[FAIL] Case 13h: found %d of 2 repeated 11-base windows "
                            "in '%s'\n", (int)n_rep, argv[1]);
                }
            }

            const int32_t numReads = (int32_t)rv.size();
            int32_t max_readlength = 0;
            int32_t n_acgt_first = 0;   /* reads the SMEM pass's round 1 walks from base 0 */
            std::vector<bool> acgt_first((size_t)numReads, false);
            std::vector<const char *> reads((size_t)numReads);
            for (int32_t i = 0; i < numReads; i++) {
                const std::string &r = rv[(size_t)i];
                reads[(size_t)i] = r.c_str();
                if ((int32_t)r.size() > max_readlength) max_readlength = (int32_t)r.size();
                if (!r.empty() && strchr("ACGT", r[0]) != NULL) {
                    acgt_first[(size_t)i] = true;
                    n_acgt_first++;
                }
            }
            uint8_t *enc_qdb; bseq1_t *seq_; int32_t *cum_len;
            encode_reads(reads.data(), numReads, &enc_qdb, &seq_, &cum_len);

            /* need_multi: some record must have s >= 2, so a record whose l is
             * compared is exercised (the repeated 11-base windows above). */
            struct HintRun {
                int32_t min_seed, max_intv, rec_intv, rec_min_len_delta, smem_min_intv;
                int32_t floor;
                bool need_multi;
            };
            const HintRun runs[] = {
                {20,  1,  1, 0, 1,  0, false},
                {20,  5,  5, 0, 1, 20, false},
                {20, 20, 20, 0, 1, 20, false},
                {11,  5,  5, 0, 1, 20, true},
                {11, 20, 20, 0, 1, 20, true},
                {20, 20,  5, 0, 1, 20, false},   /* hint max_intv differs: not used */
                {20,  1, 20, 0, 1,  0, false},   /* hint max_intv differs, would emit: not used */
                {20, 20, 20, 1, 1, 20, false},   /* hint min_len differs: not used */
                {20, 20, 20, 0, 2, 20, false},   /* SMEM pass at min_intv 2: no records */
                {20,  0,  0, 0, 1,  0, false},   /* hint max_intv 0: no records */
            };
            const size_t smem_cap = (size_t)numReads * (size_t)max_readlength;
            std::vector<BsdFirstWalk> rec((size_t)numReads);
            std::vector<SMEM> smem_out(smem_cap), smem_plain(smem_cap);
            std::vector<int32_t> mia((size_t)numReads), rid((size_t)numReads);
            std::vector<int32_t> mia_p((size_t)numReads), rid_p((size_t)numReads);
            /* The runs go first on the index with its pac loaded, where the
             * SMEM walk's unique-match text compare writes records, then on the
             * same index loaded without its pac (as the --meth seed index is),
             * where every record comes from a forward step. Each index gets its
             * own record-count and canary checks. */
            FMI_search *fmi_np = new FMI_search(argv[1]);
            fmi_np->load_index(/*load_pac=*/false);
            if (fmi_np->s1_text_available()) {
                total_cases++;
                fprintf(stderr, "[FAIL] Case 13h: the index loaded without its pac still "
                        "reports s1_text_available, so the textless runs would not be textless\n");
            }
            for (FMI_search *idx : {fmi, fmi_np})
            for (const HintRun &hr : runs) {
                const bool textless = idx == fmi_np;
                for (int32_t i = 0; i < numReads; i++) {
                    /* Garbage a correct producer overwrites or marks unrecorded. */
                    rec[(size_t)i].k = 1; rec[(size_t)i].l = 2; rec[(size_t)i].s = 1;
                    rec[(size_t)i].n = 5; rec[(size_t)i].x = 7;
                    mia[(size_t)i] = mia_p[(size_t)i] = hr.smem_min_intv;
                    rid[(size_t)i] = rid_p[(size_t)i] = i;
                }
                BsdFirstWalkHint hint;
                hint.rec      = rec.data();
                hint.max_intv = hr.rec_intv;
                hint.min_len  = hr.min_seed + hr.rec_min_len_delta;
                int64_t n_smem = 0, n_plain = 0;
                /* As mem_collect_smem: the SMEM pass's minSeedLen is one less
                 * than the bwtSeed pass's. The unhinted run gets its own
                 * min_intv and rid arrays, which the driver compacts in place. */
                idx->getSMEMsAllPosOneThread(enc_qdb, mia.data(), rid.data(), numReads, numReads,
                                             seq_, cum_len, max_readlength, hr.min_seed - 1,
                                             smem_out.data(), &n_smem, &hint);
                idx->getSMEMsAllPosOneThread(enc_qdb, mia_p.data(), rid_p.data(), numReads, numReads,
                                             seq_, cum_len, max_readlength, hr.min_seed - 1,
                                             smem_plain.data(), &n_plain, nullptr);
                const bool expect_any = hr.smem_min_intv == 1 && hr.rec_intv >= 1;
                const int32_t expect_rec = expect_any ? n_acgt_first : 0;
                /* wrong_owner counts reads whose recorded state is not the
                 * expected one, so a record written for an N-first read cannot
                 * stand in for a missing one on an ACGT-first read. */
                int32_t recorded = 0, multi = 0, wrong_owner = 0, first_wrong = -1;
                for (int32_t i = 0; i < numReads; i++) {
                    const bool is_rec = rec[(size_t)i].x > 0;
                    recorded += is_rec;
                    multi += is_rec && rec[(size_t)i].s >= 2;
                    if (is_rec != (expect_any && acgt_first[(size_t)i])) {
                        wrong_owner++;
                        if (first_wrong < 0) first_wrong = i;
                    }
                }
                char name[224];
                snprintf(name, sizeof name,
                         "Case 13h%s: BsdFirstWalk hint, minSeed=%d max_intv=%d (hint max_intv=%d "
                         "min_len=%d, SMEM min_intv=%d, %d/%d recorded, %d with s >= 2)",
                         textless ? " (no pac)" : "",
                         (int)hr.min_seed, (int)hr.max_intv, (int)hint.max_intv, (int)hint.min_len,
                         (int)hr.smem_min_intv, (int)recorded, (int)numReads, (int)multi);

                /* The hint must leave the SMEM pass's own output alone. */
                total_cases++;
                const int64_t smism = n_smem != n_plain ? 0
                                    : smem_all_fields_first_mismatch(smem_out.data(),
                                                                     smem_plain.data(), n_smem);
                if (smism >= 0) {
                    fprintf(stderr, "[FAIL] %s: hinted SMEM pass differs from the unhinted one "
                            "(n_smem %lld vs %lld, first difference at %lld)\n", name,
                            (long long)n_smem, (long long)n_plain, (long long)smism);
                } else {
                    fprintf(stderr, "[PASS] %s: SMEM pass unchanged by the hint (n_smem=%lld)\n",
                            name, (long long)n_smem);
                    passed_cases++;
                }

                /* Every read the SMEM pass walks from base 0 at min_intv 1 is
                 * recorded (each exit of that walk writes its record), and no
                 * other read is: checked read by read and in total. The parity
                 * run below goes ahead either way. */
                if (recorded != expect_rec || wrong_owner > 0 || (hr.need_multi && multi == 0)) {
                    total_cases++;
                    fprintf(stderr, "[FAIL] %s: expected %d recorded%s", name, (int)expect_rec,
                            hr.need_multi && multi == 0 ? ", some with s >= 2" : "");
                    if (wrong_owner > 0)
                        fprintf(stderr, "; %d reads have the wrong recorded state (first: read %d, "
                                "%s)", (int)wrong_owner, (int)first_wrong,
                                rec[(size_t)first_wrong].x > 0 ? "recorded but should not be"
                                                               : "not recorded but should be");
                    fprintf(stderr, "\n");
                }
                run_case(idx, name, numReads, max_readlength, hr.min_seed, hr.max_intv,
                         enc_qdb, seq_, cum_len, /*expect_min_seeds=*/hr.floor,
                         /*expect_seed_per_read=*/false, /*expect_end_seed=*/false,
                         &hint);

                /* Canary: parity alone also passes if the bwtSeed pass ignores
                 * the hint and walks every read from x = 0, or uses only some
                 * reads' records. Where the hint matches the call and seeds are
                 * expected, shift k of one recorded seed at a time: that read's
                 * lockstep output must then leave the scalar's, which it can
                 * only do by emitting that read's record. Records with s = 0
                 * emit no seed, so this check does not cover them. */
                if (hr.rec_intv != hr.max_intv || hr.rec_min_len_delta != 0 || hr.floor == 0 ||
                    expect_rec == 0)
                    continue;
                total_cases++;
                std::vector<SMEM> c_scalar(smem_cap), c_lockstep(smem_cap);
                std::vector<int32_t> mia_s((size_t)numReads, hr.max_intv);
                const int64_t c_ns = idx->bwtSeedStrategyAllPosOneThread(
                        enc_qdb, mia_s.data(), numReads, seq_, cum_len, hr.min_seed,
                        c_scalar.data());
                /* The scalar seeds of read r, which the driver emits in read order. */
                auto read_seeds = [](const std::vector<SMEM> &a, int64_t n, int32_t r) {
                    std::vector<SMEM> o;
                    for (int64_t i = 0; i < n; i++)
                        if ((int32_t)a[(size_t)i].rid == r) o.push_back(a[(size_t)i]);
                    return o;
                };
                int32_t tampered = 0, unchanged = 0, first_unchanged = -1;
                for (int32_t r = 0; r < numReads; r++) {
                    if (!(rec[(size_t)r].x > 0 && rec[(size_t)r].s > 0)) continue;
                    std::vector<BsdFirstWalk> bad(rec);
                    bad[(size_t)r].k += 1;
                    tampered++;
                    BsdFirstWalkHint bad_hint = hint;
                    bad_hint.rec = bad.data();
                    std::vector<int32_t> mia_l((size_t)numReads, hr.max_intv);
                    const int64_t c_nl = idx->bwtSeedStrategyAllPosOneThread_lockstep(
                            enc_qdb, mia_l.data(), numReads, seq_, cum_len, hr.min_seed,
                            c_lockstep.data(), max_readlength, &bad_hint);
                    const std::vector<SMEM> s_r = read_seeds(c_scalar, c_ns, r);
                    const std::vector<SMEM> l_r = read_seeds(c_lockstep, c_nl, r);
                    bool differs = s_r.size() != l_r.size();
                    for (size_t i = 0; !differs && i < s_r.size(); i++)
                        differs = !smem_fields_equal(s_r[i], l_r[i], &bad_hint);
                    if (!differs) {
                        unchanged++;
                        if (first_unchanged < 0) first_unchanged = r;
                    }
                }
                if (tampered == 0 || unchanged > 0) {
                    if (tampered == 0)
                        fprintf(stderr, "[FAIL] %s: canary: no recorded seed to tamper with\n", name);
                    else
                        fprintf(stderr, "[FAIL] %s: canary: %d of %d tampered records left their "
                                "read's output unchanged (first: read %d), so the hint is not "
                                "used for every recorded read\n",
                                name, (int)unchanged, (int)tampered, (int)first_unchanged);
                } else {
                    fprintf(stderr, "[PASS] %s: canary (each of %d tampered records changes "
                            "its read's output)\n", name, (int)tampered);
                    passed_cases++;
                }
            }
            delete fmi_np;
            free(enc_qdb); free(seq_); free(cum_len);
        }
        free(ref);
    }
#else
    fprintf(stderr, "[SKIP] Case 13h: needs the lockstep SMEM driver (SMEM_LOCKSTEP_N > 1)\n");
#endif


    fprintf(stderr, "%d / %d cases passed\n", passed_cases, total_cases);
    delete fmi;
    return (passed_cases == total_cases) ? 0 : 1;
}
