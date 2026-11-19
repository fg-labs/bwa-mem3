/* Isolated parity gate for backwardExt_konly vs the full backwardExt.
 *
 * The lockstep parity test (smem_lockstep_parity_test) exercises backwardExt_konly
 * end-to-end through the driver, but only at whatever (sp, s, a) real backward
 * walks happen to produce. This test hammers the kernel DIRECTLY across a
 * systematic sweep, with emphasis on the s==1 fast path — a single BWT-bit test
 * whose MSB-first indexing (bit CP_MASK-pos of one_hot_bwt_str[a]) had a real
 * byte-identity bug during development. Sweeping sp over every block offset
 * pos in [0, CP_MASK] locks that indexing at both block edges (pos==0 reads the
 * MSB with an empty position mask; pos==CP_MASK reads the LSB with a full mask).
 *
 * Contract checked (must hold for every (sp, s, a)):
 *   - backwardExt_konly(smem, a).s == backwardExt(smem, a).s   (always)
 *   - backwardExt_konly(smem, a).k == backwardExt(smem, a).k   (when s' > 0;
 *     konly deliberately leaves k stale on the s'->0 exit, see FMI_search.h)
 *
 * We do NOT require (sp, s) to be a real FM interval — backwardExt's k/s output
 * is pure occ arithmetic over cp_occ/count, well-defined for any in-bounds
 * (sp, sp+s), and both implementations must compute it identically.
 *
 * Under BWA3_COMPACT_OCC=1 (the compact-layout leg) both sides read the same
 * in-memory planes, so this is a self-consistency check of the two compact steps
 * only: a plane-decoding bug common to both would pass here. The independent
 * check of the compact layout is backwardext_parity_test, whose oracle reads the
 * file's own records (cp_file_oracle.h).
 *
 * Usage: backwardext_konly_parity_test <bwa-mem3 index prefix>
 *   e.g. backwardext_konly_parity_test fixtures/phix.fa
 */

/* Reach backwardExt / backwardExt_konly and the FM-index members they read.
 * They are private inline methods with no public shim; widening access in this
 * one test TU is the least-invasive way to exercise them in isolation. */
#define private public
#include "FMI_search.h"
#undef private
#include "bwa.h"

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

int main(int argc, char *argv[]) {
    if (argc != 2) {
        fprintf(stderr, "Usage: %s <bwa-mem3 index prefix>\n", argv[0]);
        return 2;
    }

    FMI_search *fmi = new FMI_search(argv[1]);
    fmi->set_compact_occ(FMI_search::compact_occ_env());   // the compact-layout test legs
    fmi->load_index();
    /* The layout this run loaded, which test/run_unit_tests.sh checks against BWA3_COMPACT_OCC. */
    fprintf(stderr, "occurrence checkpoint layout: %s\n", fmi->compact_occ() ? "compact" : "default");

    const int64_t ref_len = fmi->reference_seq_len;
    fprintf(stderr, "backwardext_konly_parity_test: reference_seq_len=%lld CP_MASK=%d\n",
            (long long)ref_len, (int)CP_MASK);

    /* Block bases spread across the index so the sweep lands in different
     * cp_occ blocks, and offsets 0..CP_MASK so every in-block bit position of
     * the s==1 fast path is exercised (both edges included). */
    const int64_t bases[] = {0, CP_BLOCK_SIZE, 10 * CP_BLOCK_SIZE,
                             ref_len / 2, ref_len - 2 * CP_BLOCK_SIZE};
    const int32_t s_values[] = {1, 2, 3, 5, 8, 17};  /* 1 = fast path; rest = general path */

    int64_t checked = 0, s_mismatch = 0, k_mismatch = 0;
    int64_t s1_checked = 0;
    /* The default layout's column form (backwardExt_konly_col, the SMEM backward
     * step's general loops) must equal backwardExt_konly_cp field for field,
     * the stale k of an s'->0 exit included: same words, same arithmetic. */
    const bool col_on = !fmi->compact_occ();
    int64_t col_checked = 0, col_mismatch = 0;

    for (size_t bi = 0; bi < sizeof(bases) / sizeof(bases[0]); bi++) {
        for (int32_t off = 0; off <= (int32_t)CP_MASK; off++) {
            const int64_t sp = bases[bi] + off;
            if (sp < 0) continue;
            for (size_t si = 0; si < sizeof(s_values) / sizeof(s_values[0]); si++) {
                const int32_t s = s_values[si];
                /* keep sp and ep = sp + s strictly in-bounds for cp_occ reads */
                if (sp + s >= ref_len) continue;
                for (uint8_t a = 0; a < 4; a++) {
                    SMEM smem;
                    smem.rid = 0;
                    smem.m = 0;
                    smem.n = 0;
                    smem.k = sp;
                    smem.l = 12345;  /* sentinel: konly must not depend on l */
                    smem.s = s;

                    SMEM full  = fmi->backwardExt(smem, a);
                    SMEM konly = fmi->backwardExt_konly(smem, a);

                    checked++;
                    if (s == 1) s1_checked++;

                    if (col_on) {
                        const SMEM ref = FMI_search::backwardExt_konly_cp(fmi->cp_occ, fmi->count[a], smem, a);
                        const SMEM col = FMI_search::backwardExt_konly_col(FMI_search::cp_col(fmi->cp_occ, a),
                                                                           fmi->count[a], smem);
                        col_checked++;
                        if (col.k != ref.k || col.s != ref.s || col.l != ref.l || col.m != ref.m ||
                            col.n != ref.n || col.rid != ref.rid) {
                            if (col_mismatch < 10)
                                fprintf(stderr, "[FAIL] column form mismatch: sp=%lld s=%d a=%d "
                                        "cp k=%lld s=%lld col k=%lld s=%lld\n",
                                        (long long)sp, s, a, (long long)ref.k, (long long)ref.s,
                                        (long long)col.k, (long long)col.s);
                            col_mismatch++;
                        }
                    }

                    if (konly.s != full.s) {
                        if (s_mismatch < 10)
                            fprintf(stderr, "[FAIL] s mismatch: sp=%lld pos=%lld s=%d a=%d "
                                    "full.s=%lld konly.s=%lld\n",
                                    (long long)sp, (long long)(sp & CP_MASK), s, a,
                                    (long long)full.s, (long long)konly.s);
                        s_mismatch++;
                        continue;
                    }
                    /* k is only defined on the s' > 0 path (konly leaves it stale on s'->0) */
                    if (full.s > 0 && konly.k != full.k) {
                        if (k_mismatch < 10)
                            fprintf(stderr, "[FAIL] k mismatch: sp=%lld pos=%lld s=%d a=%d "
                                    "full.k=%lld konly.k=%lld (s'=%lld)\n",
                                    (long long)sp, (long long)(sp & CP_MASK), s, a,
                                    (long long)full.k, (long long)konly.k, (long long)full.s);
                        k_mismatch++;
                    }
                }
            }
        }
    }

    fprintf(stderr, "checked=%lld (s==1 fast-path checks=%lld)  s_mismatch=%lld  k_mismatch=%lld\n",
            (long long)checked, (long long)s1_checked,
            (long long)s_mismatch, (long long)k_mismatch);
    fprintf(stderr, "column form: checked=%lld mismatch=%lld%s\n", (long long)col_checked,
            (long long)col_mismatch, col_on ? "" : " (compact layout: not run)");

    delete fmi;

    if (checked == 0) {
        fprintf(stderr, "[FAIL] no comparisons ran (reference too small?)\n");
        return 1;
    }
    if (col_on && col_checked == 0) {
        fprintf(stderr, "[FAIL] no column-form comparisons ran\n");
        return 1;
    }
    if (s_mismatch != 0 || k_mismatch != 0 || col_mismatch != 0) {
        fprintf(stderr, "backwardExt_konly PARITY FAILED\n");
        return 1;
    }
    fprintf(stderr, "backwardExt_konly PARITY PASS\n");
    return 0;
}
