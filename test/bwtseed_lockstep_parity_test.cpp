/* Parity harness: scalar bwtSeedStrategyAllPosOneThread vs the lockstep
 * variant. Both functions are given identical inputs; outputs (matchArray
 * bytes, returned numSeed) must be byte-identical.
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
#include "FMI_search.h"
#include "bwa.h"

/* A degenerate N=1 build would invoke bwtSeedStrategyAllPosOneThread_lockstep
 * as a single-slot path, so the cases below would trivially pass without ever
 * exercising slot-interleaving. Fail the build early in that case. */
#if BWTSEED_LOCKSTEP_N <= 1
#error "bwtseed_lockstep_parity_test requires BWTSEED_LOCKSTEP_N > 1 to exercise the lockstep path"
#endif

/* Field-by-field SMEM equality. memcmp is unsuitable because the SMEM struct
 * has implicit padding between `n` and `k` that upstream leaves
 * uninitialized — the padding bytes are not part of the semantic contract. */
static bool smem_fields_equal(const SMEM &a, const SMEM &b) {
    return a.rid == b.rid &&
           a.m   == b.m   &&
           a.n   == b.n   &&
           a.k   == b.k   &&
           a.l   == b.l   &&
           a.s   == b.s;
}

static int64_t smem_array_first_mismatch(const SMEM *a, const SMEM *b, int64_t n) {
    for (int64_t i = 0; i < n; i++) {
        if (!smem_fields_equal(a[i], b[i])) return i;
    }
    return -1;
}

static int32_t total_cases = 0;
static int32_t passed_cases = 0;

static bool run_case(FMI_search *fmi,
                     const char *case_name,
                     int32_t numReads,
                     int32_t max_readlength,
                     int32_t minSeedLen,
                     int32_t max_intv,
                     uint8_t *enc_qdb,
                     const bseq1_t *seq_,
                     int32_t *query_cum_len_ar,
                     int32_t expect_min_seeds = 0)
{
    /* Scalar emits at most one SMEM per outer-x ∈ [0, l_seq), so
     * numReads * max_readlength is a tight upper bound for both paths. */
    const int32_t max_out = numReads * max_readlength;
    SMEM *scalar_out   = (SMEM *)_mm_malloc(max_out * sizeof(SMEM), 64);
    SMEM *lockstep_out = (SMEM *)_mm_malloc(max_out * sizeof(SMEM), 64);

    /* max_intv_array per-read; mutated copies for each path so a buggy
     * implementation that writes to its input array would still parity-check. */
    int32_t *mia_s = (int32_t *)malloc(numReads * sizeof(int32_t));
    int32_t *mia_l = (int32_t *)malloc(numReads * sizeof(int32_t));
    for (int i = 0; i < numReads; i++) { mia_s[i] = max_intv; mia_l[i] = max_intv; }

    int64_t n_scalar = fmi->bwtSeedStrategyAllPosOneThread(
            enc_qdb, mia_s, numReads, seq_, query_cum_len_ar,
            minSeedLen, scalar_out);
    int64_t n_lockstep = fmi->bwtSeedStrategyAllPosOneThread_lockstep(
            enc_qdb, mia_l, numReads, seq_, query_cum_len_ar,
            minSeedLen, lockstep_out, max_readlength);

    total_cases++;
    bool ok = true;
    if (n_scalar != n_lockstep) {
        fprintf(stderr, "[FAIL] %s: numSeed mismatch: scalar=%lld lockstep=%lld\n",
                case_name, (long long)n_scalar, (long long)n_lockstep);
        ok = false;
    } else {
        int64_t mism = smem_array_first_mismatch(scalar_out, lockstep_out, n_scalar);
        if (mism >= 0) {
            const SMEM &a = scalar_out[mism];
            const SMEM &b = lockstep_out[mism];
            fprintf(stderr, "[FAIL] %s: SMEM field mismatch at index %lld:\n"
                    "  scalar   rid=%d m=%u n=%u k=%lld l=%lld s=%lld\n"
                    "  lockstep rid=%d m=%u n=%u k=%lld l=%lld s=%lld\n",
                    case_name, (long long)mism,
                    a.rid, a.m, a.n, (long long)a.k, (long long)a.l, (long long)a.s,
                    b.rid, b.m, b.n, (long long)b.k, (long long)b.l, (long long)b.s);
            ok = false;
        }
    }

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

    if (ok) {
        fprintf(stderr, "[PASS] %s  (n_seed=%lld)\n", case_name, (long long)n_scalar);
        passed_cases++;
    }

    _mm_free(scalar_out);
    _mm_free(lockstep_out);
    free(mia_s); free(mia_l);
    return ok;
}

/* Encode a vector of read strings (ACGT->0123, N->4) into the flat enc_qdb
 * + bseq1_t[] + query_cum_len_ar[] layout the scalar expects. */
static void encode_reads(const char * const *reads, int32_t numReads,
                         uint8_t **out_enc_qdb,
                         bseq1_t **out_seq,
                         int32_t **out_cum_len)
{
    int32_t total_len = 0;
    for (int i = 0; i < numReads; i++) total_len += (int32_t)strlen(reads[i]);
    uint8_t *enc = (uint8_t *)calloc(total_len > 0 ? total_len : 1, 1);
    bseq1_t *seq = (bseq1_t *)calloc(numReads, sizeof(bseq1_t));
    int32_t *cum = (int32_t *)calloc(numReads + 1, sizeof(int32_t));
    int32_t off = 0;
    for (int i = 0; i < numReads; i++) {
        int32_t l = (int32_t)strlen(reads[i]);
        cum[i] = off;
        for (int j = 0; j < l; j++) {
            char c = reads[i][j];
            uint8_t e = 4;
            if (c == 'A' || c == 'a') e = 0;
            else if (c == 'C' || c == 'c') e = 1;
            else if (c == 'G' || c == 'g') e = 2;
            else if (c == 'T' || c == 't') e = 3;
            enc[off + j] = e;
        }
        seq[i].l_seq = l;
        off += l;
    }
    cum[numReads] = off;
    *out_enc_qdb = enc;
    *out_seq = seq;
    *out_cum_len = cum;
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

    /* Case 3: numReads == N — steady state, no slot recycle. */
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
        run_case(fmi, "Case 3: numReads == N (small)",
                 numReads, max_readlength, DEFAULT_MIN_SEED, DEFAULT_MAX_INTV,
                 enc_qdb, seq_, cum_len);
        free(enc_qdb); free(seq_); free(cum_len);
    }

    /* Case 4: numReads > N — exercises the flush-then-pull-next-input recycle. */
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
        run_case(fmi, "Case 4: numReads > N (slot recycle)",
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
     * Same parameters as 10p, more reads to force multiple slot recycles
     * while match_buf contents are non-trivial. */
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
     * contents. With BWTSEED_LOCKSTEP_N=8, 20 reads guarantees the first 8
     * slots retire while there are still 12 pending inputs to pull in. The
     * recycle path must reset match_count, ready, and phase for each
     * recycled slot without corrupting in-flight slots' state. */
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
        run_case(fmi, "Case 11r: numReads > N, slot recycle on non-empty match_buf",
                 numReads, max_readlength, /*minSeed=*/10, /*max_intv=*/5,
                 enc_qdb, seq_, cum_len, /*expect_min_seeds=*/20);
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

    fprintf(stderr, "%d / %d cases passed\n", passed_cases, total_cases);
    delete fmi;
    return (passed_cases == total_cases) ? 0 : 1;
}
