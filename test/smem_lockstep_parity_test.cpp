/* Parity harness: scalar getSMEMsOnePosOneThread vs getSMEMsOnePosOneThread_lockstep.
 * Both functions are given identical inputs; outputs (matchArray bytes,
 * numTotalSmem, query_pos_array write-back) must be byte-identical.
 *
 * Usage: smem_lockstep_parity_test <bwa-mem3 index prefix>
 *   e.g. smem_lockstep_parity_test /path/to/hg38.fa
 */

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include "FMI_search.h"
#include "bwa.h"

/* A degenerate N=1 build would invoke getSMEMsOnePosOneThread_lockstep as a
 * single-slot path, so the cases below would trivially pass without ever
 * exercising slot-interleaving. Fail the build early in that case. */
#if SMEM_LOCKSTEP_N <= 1
#error "smem_lockstep_parity_test requires SMEM_LOCKSTEP_N > 1 to exercise the lockstep path"
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
                     uint8_t *enc_qdb,
                     int16_t *query_pos_array_in,
                     int32_t *min_intv_array_in,
                     int32_t *rid_array_in,
                     const bseq1_t *seq_,
                     int32_t *query_cum_len_ar)
{
    // Tight upper bound: the scalar path emits at most O(readlength) SMEMs per
    // read (one per backward-step decision plus the final prev[0]), so
    // numReads * max_readlength safely caps both scalar_out and lockstep_out.
    const int32_t max_out = numReads * max_readlength;
    SMEM *scalar_out   = (SMEM *)_mm_malloc(max_out * sizeof(SMEM), 64);
    SMEM *lockstep_out = (SMEM *)_mm_malloc(max_out * sizeof(SMEM), 64);

    int16_t *qpa_s = (int16_t *)malloc(numReads * sizeof(int16_t));
    int16_t *qpa_l = (int16_t *)malloc(numReads * sizeof(int16_t));
    int32_t *mia_s = (int32_t *)malloc(numReads * sizeof(int32_t));
    int32_t *mia_l = (int32_t *)malloc(numReads * sizeof(int32_t));
    int32_t *rid_s = (int32_t *)malloc(numReads * sizeof(int32_t));
    int32_t *rid_l = (int32_t *)malloc(numReads * sizeof(int32_t));
    memcpy(qpa_s, query_pos_array_in, numReads * sizeof(int16_t));
    memcpy(qpa_l, query_pos_array_in, numReads * sizeof(int16_t));
    memcpy(mia_s, min_intv_array_in, numReads * sizeof(int32_t));
    memcpy(mia_l, min_intv_array_in, numReads * sizeof(int32_t));
    memcpy(rid_s, rid_array_in, numReads * sizeof(int32_t));
    memcpy(rid_l, rid_array_in, numReads * sizeof(int32_t));

    int64_t n_scalar = 0, n_lockstep = 0;

    // The lockstep variant of getSMEMsOnePosOneThread heap-allocates its
    // per-slot prev/match scratch internally (sized from max_readlength),
    // so the harness no longer pre-allocates them.

    fmi->getSMEMsOnePosOneThread(enc_qdb, qpa_s, mia_s, rid_s,
                                 numReads, numReads, seq_, query_cum_len_ar,
                                 max_readlength, minSeedLen, scalar_out,
                                 &n_scalar);
    fmi->getSMEMsOnePosOneThread_lockstep(enc_qdb, qpa_l, mia_l, rid_l,
                                          numReads, numReads, seq_, query_cum_len_ar,
                                          max_readlength, minSeedLen, lockstep_out,
                                          &n_lockstep);

    total_cases++;
    bool ok = true;
    if (n_scalar != n_lockstep) {
        fprintf(stderr, "[FAIL] %s: numTotalSmem mismatch: scalar=%lld lockstep=%lld\n",
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
        } else if (memcmp(qpa_s, qpa_l, numReads * sizeof(int16_t)) != 0) {
            fprintf(stderr, "[FAIL] %s: query_pos_array write-back differs\n", case_name);
            ok = false;
        } else {
            fprintf(stderr, "[PASS] %s  (n_smem=%lld)\n", case_name, (long long)n_scalar);
            passed_cases++;
        }
    }

    _mm_free(scalar_out);
    _mm_free(lockstep_out);
    free(qpa_s); free(qpa_l); free(mia_s); free(mia_l); free(rid_s); free(rid_l);
    return ok;
}

/* Encode a vector of read strings (ACGT->0123, N->4) into the flat enc_qdb
 * + bseq1_t[] + query_cum_len_ar[] layout getSMEMsOnePosOneThread expects. */
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

    fprintf(stderr, "smem_lockstep_parity_test: SMEM_LOCKSTEP_N=%d\n",
            (int)SMEM_LOCKSTEP_N);

    FMI_search *fmi = new FMI_search(argv[1]);
    fmi->load_index();

    /* Case 1: two simple reads, start at position 0, minSeedLen 19. */
    {
        const char *reads[] = {
            "ACGTACGTACGTACGTACGTACGTACGTACGT",
            "TGCATGCATGCATGCATGCATGCATGCATGCA",
        };
        int32_t numReads = 2;
        int32_t max_readlength = 32;
        uint8_t *enc_qdb; bseq1_t *seq_; int32_t *cum_len;
        encode_reads(reads, numReads, &enc_qdb, &seq_, &cum_len);
        int16_t qpa[] = {0, 0};
        int32_t mia[] = {1, 1};
        int32_t rid[] = {0, 1};
        run_case(fmi, "Case 1: two simple reads",
                 numReads, max_readlength, 19,
                 enc_qdb, qpa, mia, rid, seq_, cum_len);
        free(enc_qdb); free(seq_); free(cum_len);
    }

    /* Case 2: numReads < N (default N=4 → use 2 reads). Exercises the
     * "unused slots" path where the driver must mark slot indices >= numReads
     * as DONE/invalid at entry so the flush pass skips them. */
    {
        const char *reads[] = {
            "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT",
            "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT",
        };
        int32_t numReads = 2;
        int32_t max_readlength = 40;
        uint8_t *enc_qdb; bseq1_t *seq_; int32_t *cum_len;
        encode_reads(reads, numReads, &enc_qdb, &seq_, &cum_len);
        int16_t qpa[] = {0, 0};
        int32_t mia[] = {1, 1};
        int32_t rid[] = {0, 1};
        run_case(fmi, "Case 2: numReads < N",
                 numReads, max_readlength, 19,
                 enc_qdb, qpa, mia, rid, seq_, cum_len);
        free(enc_qdb); free(seq_); free(cum_len);
    }

    /* Case 3: numReads == N (single full batch, no slot recycle).
     * Exercises the steady state where every slot is occupied and no
     * retirement of slots happens before flush. */
    {
        const char *reads[] = {
            "ACGTACGTACGTACGTACGTACGTACGTACGT",
            "TGCATGCATGCATGCATGCATGCATGCATGCA",
            "AAAACCCCGGGGTTTTAAAACCCCGGGGTTTT",
            "GATCGATCGATCGATCGATCGATCGATCGATC",
        };
        int32_t numReads = 4;
        int32_t max_readlength = 32;
        uint8_t *enc_qdb; bseq1_t *seq_; int32_t *cum_len;
        encode_reads(reads, numReads, &enc_qdb, &seq_, &cum_len);
        int16_t qpa[] = {0, 0, 0, 0};
        int32_t mia[] = {1, 1, 1, 1};
        int32_t rid[] = {0, 1, 2, 3};
        run_case(fmi, "Case 3: numReads == N",
                 numReads, max_readlength, 19,
                 enc_qdb, qpa, mia, rid, seq_, cum_len);
        free(enc_qdb); free(seq_); free(cum_len);
    }

    /* Case 4: numReads == N+1 (exactly one slot recycle).
     * Exercises the flush-then-pull-next-input path. */
    {
        const char *reads[] = {
            "ACGTACGTACGTACGTACGTACGTACGTACGT",
            "TGCATGCATGCATGCATGCATGCATGCATGCA",
            "AAAACCCCGGGGTTTTAAAACCCCGGGGTTTT",
            "GATCGATCGATCGATCGATCGATCGATCGATC",
            "CAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGT",
        };
        int32_t numReads = 5;
        int32_t max_readlength = 32;
        uint8_t *enc_qdb; bseq1_t *seq_; int32_t *cum_len;
        encode_reads(reads, numReads, &enc_qdb, &seq_, &cum_len);
        int16_t qpa[] = {0, 0, 0, 0, 0};
        int32_t mia[] = {1, 1, 1, 1, 1};
        int32_t rid[] = {0, 1, 2, 3, 4};
        run_case(fmi, "Case 4: numReads == N+1",
                 numReads, max_readlength, 19,
                 enc_qdb, qpa, mia, rid, seq_, cum_len);
        free(enc_qdb); free(seq_); free(cum_len);
    }

    /* Case 5: first base non-ACGT (N character). The scalar path skips the
     * whole read when enc_qdb[offset + x] >= 4; lockstep must match by
     * going straight to PH_DONE with zero matches, without writing anything
     * to matchArray for this read. Mixed with valid reads to ensure the
     * zero-match case doesn't corrupt the flush cursor. */
    {
        const char *reads[] = {
            "NACGTACGTACGTACGTACGTACGTACGTACG",
            "TGCATGCATGCATGCATGCATGCATGCATGCA",
            "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN",
            "GATCGATCGATCGATCGATCGATCGATCGATC",
        };
        int32_t numReads = 4;
        int32_t max_readlength = 32;
        uint8_t *enc_qdb; bseq1_t *seq_; int32_t *cum_len;
        encode_reads(reads, numReads, &enc_qdb, &seq_, &cum_len);
        int16_t qpa[] = {0, 0, 0, 0};
        int32_t mia[] = {1, 1, 1, 1};
        int32_t rid[] = {0, 1, 2, 3};
        run_case(fmi, "Case 5: first base non-ACGT",
                 numReads, max_readlength, 19,
                 enc_qdb, qpa, mia, rid, seq_, cum_len);
        free(enc_qdb); free(seq_); free(cum_len);
    }

    /* Case 6: numReads not divisible by N (tail partial batch).
     * After filling N=4 slots, 3 more reads get pulled in as slots retire.
     * When only fewer-than-N reads remain and the last ones finish, the
     * driver must correctly mark the slots retired (input_idx = -1) without
     * trying to recycle from an exhausted queue. */
    {
        const char *reads[] = {
            "ACGTACGTACGTACGTACGTACGTACGTACGT",
            "TGCATGCATGCATGCATGCATGCATGCATGCA",
            "AAAACCCCGGGGTTTTAAAACCCCGGGGTTTT",
            "GATCGATCGATCGATCGATCGATCGATCGATC",
            "CAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGT",
            "TTTACCCAGGGCTTTACCCAGGGCTTTACCCA",
            "ATATATATATATATATATATATATATATATAT",
        };
        int32_t numReads = 7;  /* 4 initial + 3 recycles; not divisible by N=4 */
        int32_t max_readlength = 32;
        uint8_t *enc_qdb; bseq1_t *seq_; int32_t *cum_len;
        encode_reads(reads, numReads, &enc_qdb, &seq_, &cum_len);
        int16_t qpa[] = {0, 0, 0, 0, 0, 0, 0};
        int32_t mia[] = {1, 1, 1, 1, 1, 1, 1};
        int32_t rid[] = {0, 1, 2, 3, 4, 5, 6};
        run_case(fmi, "Case 6: numReads not divisible by N",
                 numReads, max_readlength, 19,
                 enc_qdb, qpa, mia, rid, seq_, cum_len);
        free(enc_qdb); free(seq_); free(cum_len);
    }

    /* Case 7: mixed start_pos (some reads start mid-read rather than at 0).
     * Exercises the forward-ext `j = start_pos + 1` initialization and the
     * backward-ext `j = start_pos - 1` initialization for non-zero starts.
     * The scalar path uses `query_pos_array[i]` directly; lockstep carries
     * it in slot.start_pos, which must round-trip correctly. */
    {
        const char *reads[] = {
            "ACGTACGTACGTACGTACGTACGTACGTACGT",
            "TGCATGCATGCATGCATGCATGCATGCATGCA",
            "AAAACCCCGGGGTTTTAAAACCCCGGGGTTTT",
            "GATCGATCGATCGATCGATCGATCGATCGATC",
        };
        int32_t numReads = 4;
        int32_t max_readlength = 32;
        uint8_t *enc_qdb; bseq1_t *seq_; int32_t *cum_len;
        encode_reads(reads, numReads, &enc_qdb, &seq_, &cum_len);
        int16_t qpa[] = {5, 10, 0, 15};   /* mixed start positions */
        int32_t mia[] = {1, 1, 1, 1};
        int32_t rid[] = {0, 1, 2, 3};
        run_case(fmi, "Case 7: mixed start_pos",
                 numReads, max_readlength, 19,
                 enc_qdb, qpa, mia, rid, seq_, cum_len);
        free(enc_qdb); free(seq_); free(cum_len);
    }

    /* Case 8: real 300bp Illumina WGS read pair from SRR6109255
     * (SRR6109255.100035 R1/R2). Exercises long SMEM walks through adapter
     * contamination (AGATCGGAAGAG around base ~22-24) and low-complexity
     * A-run tails — patterns the synthetic Cases 1-7 don't hit. Regression
     * coverage for real-world read composition. */
    {
        /* SRR6109255.100035 R1 (300bp) */
        const char *r1 =
            "ATGACCTCCCTAATATCTTCAGAATAGATCGGAAGAGCACACGTCTGAACTCCAGTCACT"
            "ACAAGATCTCGTATGCCGTCTTCTGCTTGAAAAAAAAAAAAAAAAAATGAGTAACTCTAT"
            "GAATCGAACTCAAACATTCCACACACGTTTATGCCACTACTTCTATCGTGCTTTCATATT"
            "ATACTACTCCCCCTTCCCCCTCCATCTTCACCTCTCCTCCAATCTCCACTCACCTCTACC"
            "TCAACACGTACTTTACACACATCTCTCCCCCACGGACCACAATACCTCTCCTCTCAATTA";
        /* SRR6109255.100035 R2 (300bp) */
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
        int16_t qpa[] = {0, 0};
        int32_t mia[] = {1, 1};
        int32_t rid[] = {0, 1};
        run_case(fmi, "Case 8: real 300bp SRR6109255.100035",
                 numReads, max_readlength, 19,
                 enc_qdb, qpa, mia, rid, seq_, cum_len);
        free(enc_qdb); free(seq_); free(cum_len);
    }

    /* Case 9: same real 300bp pair, driven through the outer
     * `getSMEMsAllPosOneThread` do/while loop. The outer driver compacts
     * active reads and re-calls the inner function with updated
     * query_pos_array[i]; a bug in the inner function's write-back of
     * next_x (slot.next_x via ls_advance_forward_step / ls_init_slot) would
     * show up here even if Case 8 passes. We run the scalar and lockstep
     * outer drivers locally rather than via getSMEMsAllPosOneThread, since
     * that one dispatches based on the SMEM_LOCKSTEP_N macro at compile
     * time — both sides would take the same branch. */
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
        int32_t minSeedLen = 19;
        uint8_t *enc_qdb; bseq1_t *seq_; int32_t *cum_len;
        encode_reads(reads, numReads, &enc_qdb, &seq_, &cum_len);
        const int32_t max_out = numReads * max_readlength;
        SMEM *scalar_out   = (SMEM *)_mm_malloc(max_out * sizeof(SMEM), 64);
        SMEM *lockstep_out = (SMEM *)_mm_malloc(max_out * sizeof(SMEM), 64);
        /* The lockstep variant now owns its per-slot prev/match scratch
         * internally (sized by MAX_READ_LEN_FOR_LOCKSTEP); no external
         * lockstep buffers are needed here. */

        /* The two do/while blocks below intentionally mirror the compaction
         * logic in FMI_search::getSMEMsAllPosOneThread (see src/FMI_search.cpp),
         * in particular the `query_pos_array[head] < seq_[rid].l_seq` predicate
         * that keeps reads active. If that predicate or the surrounding loop
         * structure changes in production, update both drivers here to match. */
        /* Scalar outer driver. */
        int64_t n_s = 0;
        int16_t qpa_s[2] = {0, 0};
        {
            int32_t rid_work[2] = {0, 1};
            int32_t mia_work[2] = {1, 1};
            int32_t numActive = numReads;
            do {
                int32_t head = 0, tail = 0;
                for (head = 0; head < numActive; head++) {
                    int rl = seq_[rid_work[head]].l_seq;
                    if (qpa_s[head] < rl) {
                        rid_work[tail] = rid_work[head];
                        qpa_s[tail] = qpa_s[head];
                        mia_work[tail] = mia_work[head];
                        tail++;
                    }
                }
                fmi->getSMEMsOnePosOneThread(enc_qdb, qpa_s, mia_work, rid_work,
                                             tail, tail, seq_, cum_len,
                                             max_readlength, minSeedLen, scalar_out,
                                             &n_s);
                numActive = tail;
            } while (numActive > 0);
        }

        /* Lockstep outer driver. */
        int64_t n_l = 0;
        int16_t qpa_l[2] = {0, 0};
        {
            int32_t rid_work[2] = {0, 1};
            int32_t mia_work[2] = {1, 1};
            int32_t numActive = numReads;
            do {
                int32_t head = 0, tail = 0;
                for (head = 0; head < numActive; head++) {
                    int rl = seq_[rid_work[head]].l_seq;
                    if (qpa_l[head] < rl) {
                        rid_work[tail] = rid_work[head];
                        qpa_l[tail] = qpa_l[head];
                        mia_work[tail] = mia_work[head];
                        tail++;
                    }
                }
                fmi->getSMEMsOnePosOneThread_lockstep(enc_qdb, qpa_l, mia_work, rid_work,
                                                     tail, tail, seq_, cum_len,
                                                     max_readlength, minSeedLen, lockstep_out,
                                                     &n_l);
                numActive = tail;
            } while (numActive > 0);
        }

        total_cases++;
        if (n_s != n_l) {
            fprintf(stderr, "[FAIL] Case 9 AllPos driver: numTotalSmem mismatch: scalar=%lld lockstep=%lld\n",
                    (long long)n_s, (long long)n_l);
        } else {
            int64_t mism = smem_array_first_mismatch(scalar_out, lockstep_out, n_s);
            if (mism >= 0) {
                const SMEM &a = scalar_out[mism];
                const SMEM &b = lockstep_out[mism];
                fprintf(stderr, "[FAIL] Case 9 AllPos driver: SMEM mismatch at %lld:\n"
                        "  scalar   rid=%d m=%u n=%u k=%lld l=%lld s=%lld\n"
                        "  lockstep rid=%d m=%u n=%u k=%lld l=%lld s=%lld\n",
                        (long long)mism,
                        a.rid, a.m, a.n, (long long)a.k, (long long)a.l, (long long)a.s,
                        b.rid, b.m, b.n, (long long)b.k, (long long)b.l, (long long)b.s);
            } else {
                fprintf(stderr, "[PASS] Case 9 AllPos driver  (n_smem=%lld)\n", (long long)n_s);
                passed_cases++;
            }
        }
        _mm_free(scalar_out); _mm_free(lockstep_out);
        free(enc_qdb); free(seq_); free(cum_len);
    }

    /* Case 10: exercises the forward-ext `newSmem.s < min_intv` exit path
     * by using the real 300bp reads with a large min_intv (mimicking what
     * mem_collect_smem uses for phase-2 re-seeding: min_intv = prev_smem.s
     * + 1). Synthetic Cases 1-7 use min_intv=1, which rarely triggers this
     * exit path. This case guards the post-loop "push smem if smem.s >=
     * min_intv" behavior on the min_intv branch. */
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
        /* Start mid-read (past adapter) and use a large min_intv so the
         * forward walk will terminate via the newSmem.s < min_intv branch. */
        int16_t qpa[] = {100, 100};
        int32_t mia[] = {1000, 1000};
        int32_t rid[] = {0, 1};
        run_case(fmi, "Case 10: forward-ext min_intv break",
                 numReads, max_readlength, 19,
                 enc_qdb, qpa, mia, rid, seq_, cum_len);
        free(enc_qdb); free(seq_); free(cum_len);
    }

    /* Case 11: long reads (>151bp and >512bp) — regression for issue 44.
     * Pre-fix, readlength > 512 would trip the lockstep cap and exit.
     * Post-fix, the lockstep slot buffers are heap-allocated by the
     * driver and sized from max_readlength, so any length works.
     * Exercises a 1000bp and a 1500bp synthetic read. */
    {
        char r_1000bp[1001];
        char r_1500bp[1501];
        // Deterministic non-repetitive synthetic: a 4-base rotation with a
        // per-position twist so the SMEMs don't all collapse into one.
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
        int16_t qpa[] = {0, 0};
        int32_t mia[] = {1, 1};
        int32_t rid[] = {0, 1};
        run_case(fmi, "Case 11: long reads (1000bp, 1500bp)",
                 numReads, max_readlength, 19,
                 enc_qdb, qpa, mia, rid, seq_, cum_len);
        free(enc_qdb); free(seq_); free(cum_len);
    }

    /* Case 12: long reads with mid-read N — exercises the non-ACGT
     * termination path on the long-read code path that Case 11's clean
     * 4-base rotation didn't cover. The forward SMEM walk must split at
     * the N position and re-seed past it, on the heap-allocated lockstep
     * slot buffers. Two reads: 1000bp with N at pos 700, 1500bp with N
     * at pos 1100. */
    {
        char r_1000bp[1001];
        char r_1500bp[1501];
        const char *bases = "ACGT";
        for (int i = 0; i < 1000; i++) r_1000bp[i] = bases[(i * 7 + 3) & 3];
        r_1000bp[700] = 'N';
        r_1000bp[1000] = '\0';
        for (int i = 0; i < 1500; i++) r_1500bp[i] = bases[(i * 11 + 5) & 3];
        r_1500bp[1100] = 'N';
        r_1500bp[1500] = '\0';
        const char *reads[] = { r_1000bp, r_1500bp };
        int32_t numReads = 2;
        int32_t max_readlength = 1500;
        uint8_t *enc_qdb; bseq1_t *seq_; int32_t *cum_len;
        encode_reads(reads, numReads, &enc_qdb, &seq_, &cum_len);
        int16_t qpa[] = {0, 0};
        int32_t mia[] = {1, 1};
        int32_t rid[] = {0, 1};
        run_case(fmi, "Case 12: long reads with mid-read N",
                 numReads, max_readlength, 19,
                 enc_qdb, qpa, mia, rid, seq_, cum_len);
        free(enc_qdb); free(seq_); free(cum_len);
    }

    fprintf(stderr, "%d / %d cases passed\n", passed_cases, total_cases);
    delete fmi;
    return (passed_cases == total_cases) ? 0 : 1;
}
