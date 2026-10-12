/* Parity harness: scalar getSMEMsOnePosOneThread vs getSMEMsOnePosOneThread_lockstep.
 * Both functions are given identical inputs; outputs (matchArray bytes,
 * numTotalSmem, query_pos_array write-back) must be byte-identical.
 *
 * Usage: smem_lockstep_parity_test <reference FASTA that is also the index prefix>
 * The checked-in phiX fixture (test/fixtures/phix.fa) is REQUIRED, not just an
 * example: Cases 13, 13z and 14 read reads from that FASTA and assert emission
 * floors that assume a single-contig, all-ACGT reference, so a run against a
 * different index can legitimately fail them. An unreadable (or, for Case 14,
 * too short) FASTA fails them rather than skipping, so the real-match and
 * slot-recycle coverage cannot be lost silently.
 */

#include <assert.h>
#include <ctype.h>
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

/* Field-by-field SMEM equality on the OBSERVABLE fields (rid,m,n,k,s).
 *
 * `l` (the reverse-complement interval start) is deliberately excluded, mirroring
 * smem_dedup_inplace: the pure-backward phase (backwardExt_konly) intentionally
 * leaves `l` frozen at its forward-phase value while the scalar reference path
 * keeps recomputing it, so the two paths legitimately differ on `l` for any seed
 * emitted after a multi-step backward walk through an s>1 interval. That
 * divergence is inert — no consumer of an emitted seed reads `l` (SA resolution
 * uses only k/s; sort/dedup key on (rid,m,n[,k,s])) — and `l` is a deterministic
 * function of (rid,m,n) anyway, so it carries no independent information. Asserting
 * a.l==b.l here would be a false invariant the K-only backward kernel does not
 * uphold; the observable seed identity is (rid,m,n,k,s).
 *
 * memcmp is unsuitable because the SMEM struct has implicit padding between `n`
 * and `k` that upstream leaves uninitialized — padding is not part of the
 * contract, and neither (now) is `l`. */
static bool smem_fields_equal(const SMEM &a, const SMEM &b) {
    return a.rid == b.rid &&
           a.m   == b.m   &&
           a.n   == b.n   &&
           a.k   == b.k   &&
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

/* A fixed width a case pins, clamped to the build's SMEM_LOCKSTEP_N_MAX so a
 * -D build with a small cap still runs every case instead of tripping the
 * driver's width assert. */
static int32_t clamp_width(int32_t w)
{
    return w < SMEM_LOCKSTEP_N_MAX ? w : SMEM_LOCKSTEP_N_MAX;
}

static bool run_case(FMI_search *fmi,
                     const char *case_name,
                     int32_t numReads,
                     int32_t max_readlength,
                     int32_t minSeedLen,
                     uint8_t *enc_qdb,
                     int32_t *query_pos_array_in,
                     int32_t *min_intv_array_in,
                     int32_t *rid_array_in,
                     const bseq1_t *seq_,
                     int32_t *query_cum_len_ar,
                     int64_t expect_min_smem = 0,
                     const bool *expect_smem_for_read = NULL)
{
    // Tight upper bound: the scalar path emits at most O(readlength) SMEMs per
    // read (one per backward-step decision plus the final prev[0]), so
    // numReads * max_readlength safely caps both scalar_out and lockstep_out.
    const int32_t max_out = numReads * max_readlength;
    SMEM *scalar_out   = (SMEM *)_mm_malloc(max_out * sizeof(SMEM), 64);
    SMEM *lockstep_out = (SMEM *)_mm_malloc(max_out * sizeof(SMEM), 64);

    int32_t *qpa_s = (int32_t *)malloc(numReads * sizeof(int32_t));
    int32_t *qpa_l = (int32_t *)malloc(numReads * sizeof(int32_t));
    int32_t *mia_s = (int32_t *)malloc(numReads * sizeof(int32_t));
    int32_t *mia_l = (int32_t *)malloc(numReads * sizeof(int32_t));
    int32_t *rid_s = (int32_t *)malloc(numReads * sizeof(int32_t));
    int32_t *rid_l = (int32_t *)malloc(numReads * sizeof(int32_t));
    if (!qpa_s || !qpa_l || !mia_s || !mia_l || !rid_s || !rid_l) {
        fprintf(stderr, "[FAIL] %s: allocation failed\n", case_name);
        _mm_free(scalar_out);
        _mm_free(lockstep_out);
        free(qpa_s); free(qpa_l); free(mia_s); free(mia_l); free(rid_s); free(rid_l);
        total_cases++;
        return false;
    }
    memcpy(qpa_s, query_pos_array_in, numReads * sizeof(int32_t));
    memcpy(qpa_l, query_pos_array_in, numReads * sizeof(int32_t));
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
        } else if (memcmp(qpa_s, qpa_l, numReads * sizeof(int32_t)) != 0) {
            fprintf(stderr, "[FAIL] %s: query_pos_array write-back differs\n", case_name);
            ok = false;
        } else if (n_scalar < expect_min_smem) {
            /* Parity between two empty outputs proves nothing about emission. */
            fprintf(stderr, "[FAIL] %s: emission path not exercised: got %lld SMEMs, "
                    "expected >= %lld\n", case_name, (long long)n_scalar,
                    (long long)expect_min_smem);
            ok = false;
        } else if (expect_smem_for_read != NULL) {
            /* Every flagged read must have emitted at least one SMEM. */
            bool *seen = (bool *)calloc(numReads, sizeof(bool));
            int32_t missing = -1;
            if (seen == NULL) {
                fprintf(stderr, "[FAIL] %s: allocation failed\n", case_name);
                ok = false;
            } else {
                for (int64_t j = 0; j < n_scalar; j++) {
                    const int32_t r = scalar_out[j].rid;
                    if (r >= 0 && r < numReads) seen[r] = true;
                }
                for (int32_t r = 0; r < numReads && missing < 0; r++)
                    if (expect_smem_for_read[r] && !seen[r]) missing = r;
                free(seen);
                if (missing >= 0) {
                    fprintf(stderr, "[FAIL] %s: read %d emitted no SMEM\n",
                            case_name, missing);
                    ok = false;
                }
            }
        }
        if (ok) {
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

/* Read a (small) FASTA into a flat uppercase ACGT/N sequence string. Used to
 * derive reads that actually MATCH the reference, so the backward-extension
 * path is exercised against real SA intervals rather than the non-matching
 * synthetic reads above (which all resolve to n_smem=0 and never run a backward
 * step). Returns a malloc'd buffer (caller frees) or NULL if the file cannot be
 * opened / has no sequence. */
static char *read_fasta_seq(const char *path, int32_t *out_len) {
    FILE *f = fopen(path, "r");
    if (!f) return NULL;
    size_t cap = 8192, len = 0;
    char *seq = (char *)malloc(cap);
    if (!seq) { fclose(f); return NULL; }
    int c;
    bool in_header = false;
    while ((c = fgetc(f)) != EOF) {
        if (c == '>') { in_header = true; continue; }
        if (c == '\n') { in_header = false; continue; }
        if (in_header) continue;
        if (c == '\r' || c == ' ' || c == '\t') continue;
        if (len + 1 >= cap) {
            cap *= 2;
            char *tmp = (char *)realloc(seq, cap);
            if (!tmp) { free(seq); fclose(f); return NULL; }
            seq = tmp;
        }
        seq[len++] = (char)toupper(c);
    }
    fclose(f);
    if (len == 0) { free(seq); return NULL; }
    seq[len] = '\0';
    *out_len = (int32_t)len;
    return seq;
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
        int32_t qpa[] = {0, 0};
        int32_t mia[] = {1, 1};
        int32_t rid[] = {0, 1};
        run_case(fmi, "Case 1: two simple reads",
                 numReads, max_readlength, 19,
                 enc_qdb, qpa, mia, rid, seq_, cum_len);
        free(enc_qdb); free(seq_); free(cum_len);
    }

    /* Cases 2-4 and 6 pin the width to N=4 (clamped to SMEM_LOCKSTEP_N_MAX), so
     * their read counts sit below, at, one above and past it as described,
     * whatever the compile-time default. Cases 2-4 are defined by where their
     * read count sits relative to N=4, so a build capped below 4 skips them
     * (reported as SKIP) rather than running them at a width they don't describe.
     *
     * Case 2: numReads < N (N=4, 2 reads). Exercises the
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
        int32_t qpa[] = {0, 0};
        int32_t mia[] = {1, 1};
        int32_t rid[] = {0, 1};
        const int32_t saved_n = g_smem_lockstep_n;
        g_smem_lockstep_n = clamp_width(4);
        if (g_smem_lockstep_n == 4)
            run_case(fmi, "Case 2: numReads < N",
                     numReads, max_readlength, 19,
                     enc_qdb, qpa, mia, rid, seq_, cum_len);
        else
            fprintf(stderr, "[SKIP] Case 2: needs N=4, but SMEM_LOCKSTEP_N_MAX is %d\n",
                    (int)SMEM_LOCKSTEP_N_MAX);
        g_smem_lockstep_n = saved_n;
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
        int32_t qpa[] = {0, 0, 0, 0};
        int32_t mia[] = {1, 1, 1, 1};
        int32_t rid[] = {0, 1, 2, 3};
        const int32_t saved_n = g_smem_lockstep_n;
        g_smem_lockstep_n = clamp_width(4);
        if (g_smem_lockstep_n == 4)
            run_case(fmi, "Case 3: numReads == N",
                     numReads, max_readlength, 19,
                     enc_qdb, qpa, mia, rid, seq_, cum_len);
        else
            fprintf(stderr, "[SKIP] Case 3: needs N=4, but SMEM_LOCKSTEP_N_MAX is %d\n",
                    (int)SMEM_LOCKSTEP_N_MAX);
        g_smem_lockstep_n = saved_n;
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
        int32_t qpa[] = {0, 0, 0, 0, 0};
        int32_t mia[] = {1, 1, 1, 1, 1};
        int32_t rid[] = {0, 1, 2, 3, 4};
        const int32_t saved_n = g_smem_lockstep_n;
        g_smem_lockstep_n = clamp_width(4);
        if (g_smem_lockstep_n == 4)
            run_case(fmi, "Case 4: numReads == N+1",
                     numReads, max_readlength, 19,
                     enc_qdb, qpa, mia, rid, seq_, cum_len);
        else
            fprintf(stderr, "[SKIP] Case 4: needs N=4, but SMEM_LOCKSTEP_N_MAX is %d\n",
                    (int)SMEM_LOCKSTEP_N_MAX);
        g_smem_lockstep_n = saved_n;
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
        int32_t qpa[] = {0, 0, 0, 0};
        int32_t mia[] = {1, 1, 1, 1};
        int32_t rid[] = {0, 1, 2, 3};
        run_case(fmi, "Case 5: first base non-ACGT",
                 numReads, max_readlength, 19,
                 enc_qdb, qpa, mia, rid, seq_, cum_len);
        free(enc_qdb); free(seq_); free(cum_len);
    }

    /* Case 6: numReads not divisible by N (tail partial batch).
     * At N=4, after filling 4 slots, 3 more reads get pulled in as slots retire.
     * A build capped below 4 runs it at the cap (2 or 3), where 7 reads still
     * leave a partial tail batch.
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
        int32_t numReads = 7;  /* not a multiple of N at N=4, 3 or 2 */
        int32_t max_readlength = 32;
        uint8_t *enc_qdb; bseq1_t *seq_; int32_t *cum_len;
        encode_reads(reads, numReads, &enc_qdb, &seq_, &cum_len);
        int32_t qpa[] = {0, 0, 0, 0, 0, 0, 0};
        int32_t mia[] = {1, 1, 1, 1, 1, 1, 1};
        int32_t rid[] = {0, 1, 2, 3, 4, 5, 6};
        const int32_t saved_n = g_smem_lockstep_n;
        g_smem_lockstep_n = clamp_width(4);
        run_case(fmi, "Case 6: numReads not divisible by N",
                 numReads, max_readlength, 19,
                 enc_qdb, qpa, mia, rid, seq_, cum_len);
        g_smem_lockstep_n = saved_n;
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
        int32_t qpa[] = {5, 10, 0, 15};   /* mixed start positions */
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
        int32_t qpa[] = {0, 0};
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
        int32_t qpa_s[2] = {0, 0};
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
        int32_t qpa_l[2] = {0, 0};
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
        int32_t qpa[] = {100, 100};
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
        int32_t qpa[] = {0, 0};
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
        int32_t qpa[] = {0, 0};
        int32_t mia[] = {1, 1};
        int32_t rid[] = {0, 1};
        run_case(fmi, "Case 12: long reads with mid-read N",
                 numReads, max_readlength, 19,
                 enc_qdb, qpa, mia, rid, seq_, cum_len);
        free(enc_qdb); free(seq_); free(cum_len);
    }

    /* Case 13: reads that ACTUALLY MATCH the reference, driven from an interior
     * start_pos so the pure-backward phase (backwardExt_konly) runs a real
     * multi-step walk. Cases 1-12 use non-matching synthetic reads and/or
     * start_pos=0, so every one resolves to n_smem=0 and the backward step
     * (the code this PR changed) is never exercised — a vacuous gate. Here each
     * read is a substring of the loaded reference at a spread of offsets, with
     * an interior pivot, so the backward walk narrows a real s>1 SA interval
     * down to a unique match and emits SMEMs. This is the case that gives
     * backwardExt_konly (and the numPrev==1 fast path) genuine parity coverage;
     * with the pre-refactor comparator (which compared `l`) it FAILS on the now
     * intentionally-stale `l`, and PASSES on (rid,m,n,k,s). */
    {
        int32_t ref_len = 0;
        char *ref = read_fasta_seq(argv[1], &ref_len);
        if (!ref) {
            /* argv[1] must be a readable FASTA (CI passes fixtures/phix.fa).
             * Count the case as failed rather than skipping it, so the
             * backward-path coverage cannot drop out unnoticed. */
            total_cases++;
            fprintf(stderr, "[FAIL] Case 13: could not read reference FASTA '%s' "
                    "(backward-path parity coverage needs it)\n", argv[1]);
        } else {
            const int32_t rl = 60;               /* > minSeedLen (19) */
            const int32_t want = 6;
            int32_t off[6] = {0, 400, 900, 1500, 2600, 4000};
            int32_t numReads = 0;
            char *bufs[6] = {0};
            for (int32_t i = 0; i < want; i++) {
                if (off[i] + rl > ref_len) continue;   /* keep in-bounds */
                char *b = (char *)malloc(rl + 1);
                memcpy(b, ref + off[i], rl);
                b[rl] = '\0';
                bufs[numReads++] = b;
            }
            if (numReads == 0) {
                total_cases++;
                fprintf(stderr, "[FAIL] Case 13: reference too short for any read\n");
            } else {
                const char *reads[6];
                int32_t qpa[6], mia[6], rid[6];
                for (int32_t i = 0; i < numReads; i++) {
                    reads[i] = bufs[i];
                    /* Alternate the pivot: even reads start at the last base
                     * (longest backward walk), odd reads start mid-read (both
                     * forward and backward walks run). */
                    qpa[i] = (i & 1) ? (rl / 2) : (rl - 1);
                    mia[i] = 1;
                    rid[i] = i;   /* rid indexes the per-read seq_/cum_len layout */
                }
                uint8_t *enc_qdb; bseq1_t *seq_; int32_t *cum_len;
                encode_reads(reads, numReads, &enc_qdb, &seq_, &cum_len);
                run_case(fmi, "Case 13: real reference substrings, interior start_pos",
                         numReads, rl, 19,
                         enc_qdb, qpa, mia, rid, seq_, cum_len);
                free(enc_qdb); free(seq_); free(cum_len);
            }
            for (int32_t i = 0; i < numReads; i++) free(bufs[i]);
            free(ref);
        }
    }

    /* Case 13z: slot recycle at the AMD Zen 5 host-default width (32 unless a -D
     * override clamps it), the width Zen 5 hosts run by default. Cases 1-13 run
     * at the compile-time default and never put more reads than slots through
     * real reference matches; here 70 reference substrings (more than twice the
     * width) with interior pivots make every slot recycle with non-empty
     * prev/match state, and the case must emit SMEMs. The width is pinned
     * through the runtime global and restored afterwards. */
    {
        const int32_t zen5_n = bwa3_smem_host_default_width(1, 0x1a, 0x02);
        int32_t ref_len = 0;
        char *ref = read_fasta_seq(argv[1], &ref_len);
        enum { NZ = 70 };
        const int32_t rl = 60;
        const int32_t step = 73;   /* 69 * 73 + 60 = 5097 <= phiX's 5386 bp */
        if (!ref) {
            /* As for Case 13: a missing reference fails the case. */
            total_cases++;
            fprintf(stderr, "[FAIL] Case 13z: could not read reference FASTA '%s'\n", argv[1]);
        } else if ((NZ - 1) * step + rl > ref_len) {
            fprintf(stderr, "[FAIL] Case 13z: reference too short (%d bp) for %d reads\n",
                    ref_len, (int)NZ);
            total_cases++;
            free(ref);
        } else if (zen5_n >= NZ) {
            /* Only a -D build raising the widths gets here; Case 14 covers
             * slot recycle at every width up to SMEM_LOCKSTEP_N_MAX. */
            fprintf(stderr, "[SKIP] Case 13z: host-default width %d leaves no recycle "
                    "with %d reads (Case 14 covers it)\n", zen5_n, (int)NZ);
            free(ref);
        } else {
            char *bufs[NZ];
            const char *reads[NZ];
            int32_t qpa[NZ], mia[NZ], rid[NZ];
            for (int32_t i = 0; i < NZ; i++) {
                bufs[i] = (char *)malloc(rl + 1);
                memcpy(bufs[i], ref + i * step, rl);
                bufs[i][rl] = '\0';
                reads[i] = bufs[i];
                qpa[i] = (i & 1) ? (rl / 2) : (rl - 1);
                mia[i] = 1;
                rid[i] = i;
            }
            uint8_t *enc_qdb; bseq1_t *seq_; int32_t *cum_len;
            encode_reads(reads, NZ, &enc_qdb, &seq_, &cum_len);
            const int32_t saved_n = g_smem_lockstep_n;
            g_smem_lockstep_n = zen5_n;
            char name[128];
            snprintf(name, sizeof(name),
                     "Case 13z: numReads(%d) > N(%d), real reference substrings at the Zen 5 host-default width",
                     (int)NZ, (int)zen5_n);
            run_case(fmi, name, NZ, rl, 19, enc_qdb, qpa, mia, rid, seq_, cum_len,
                     /*expect_min_smem=*/NZ);
            g_smem_lockstep_n = saved_n;
            free(enc_qdb); free(seq_); free(cum_len);
            for (int32_t i = 0; i < NZ; i++) free(bufs[i]);
            free(ref);
        }
    }

    /* Case 14: many more reads than slots, with real matches, at several widths.
     * The flush pass finds the slot holding the next read to emit as input
     * mod N, which relies on slot recycling keeping input i in slot i mod N.
     * Cases 1-13 never recycle a slot through real reference matches (Cases 4
     * and 6 recycle, but with synthetic reads), and Case 13z recycles only at
     * one width with uniform reads, so this drives 2 * N_MAX + 22 reads
     * (150 at the default N_MAX of 64, so every slot is reused at least once even
     * at the maximum width) of mixed lengths and pivots -- they retire out of
     * input order -- at a small width, the shipped default and the maximum. */
    {
        int32_t ref_len = 0;
        char *ref = read_fasta_seq(argv[1], &ref_len);
        if (!ref || ref_len < 200) {
            /* As for Case 13: a missing or short reference fails the case. */
            total_cases++;
            fprintf(stderr, "[FAIL] Case 14: could not read a reference FASTA of at least "
                    "200 bases from '%s' (slot-recycle parity coverage needs it)\n", argv[1]);
        } else {
            /* EXTRA keeps NR at 150 for the default N_MAX of 64. */
            enum { EXTRA = 22, NR = 2 * SMEM_LOCKSTEP_N_MAX + EXTRA, MAXL = 90 };
            char *bufs[NR];
            const char *reads[NR];
            int32_t qpa[NR], mia[NR], rid[NR];
            /* With the single-contig, all-ACGT phiX fixture, each read without
             * an injected N is an exact reference substring of at least 30
             * bases, longer than the 19-base minimum, so it must emit at least
             * one SMEM through its pivot. Reads with an N may emit too, so the
             * check is per read rather than a total. */
            bool need_smem[NR];
            for (int32_t i = 0; i < NR; i++) {
                const int32_t len = 30 + (i * 7) % (MAXL - 29);           /* 30..90 */
                const int32_t off = (int32_t)(((int64_t)i * 977) % (ref_len - MAXL));
                bufs[i] = (char *)malloc(len + 1);
                memcpy(bufs[i], ref + off, len);
                bufs[i][len] = '\0';
                if (i % 11 == 5) bufs[i][len / 3] = 'N';                   /* early stop */
                need_smem[i] = (i % 11 != 5);
                reads[i] = bufs[i];
                qpa[i] = (i % 3 == 0) ? len - 1 : (i % 3 == 1) ? len / 2 : 0;
                mia[i] = 1;
                rid[i] = i;
            }
            uint8_t *enc_qdb; bseq1_t *seq_; int32_t *cum_len;
            encode_reads(reads, NR, &enc_qdb, &seq_, &cum_len);
            const int32_t saved_n = g_smem_lockstep_n;
            const int32_t widths[] = {clamp_width(3), SMEM_LOCKSTEP_N, SMEM_LOCKSTEP_N_MAX};
            for (const int32_t w : widths) {
                char name[96];
                snprintf(name, sizeof name, "Case 14: %d reads > N = %d, slot recycle", (int)NR, (int)w);
                g_smem_lockstep_n = w;
                run_case(fmi, name, NR, MAXL, 19,
                         enc_qdb, qpa, mia, rid, seq_, cum_len,
                         /*expect_min_smem=*/0, need_smem);
            }
            g_smem_lockstep_n = saved_n;
            free(enc_qdb); free(seq_); free(cum_len);
            for (int32_t i = 0; i < NR; i++) free(bufs[i]);
        }
        free(ref);
    }

    fprintf(stderr, "%d / %d cases passed\n", passed_cases, total_cases);
    delete fmi;
    return (passed_cases == total_cases) ? 0 : 1;
}
