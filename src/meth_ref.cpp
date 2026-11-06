/* SPDX-License-Identifier: MIT */
#include "meth_ref.h"

#include <limits.h>
#include <stdio.h>
#include <stdlib.h>

#include "read_index_ele.h"

/* --meth (mem and aln): load the ORIGINAL reference's bns + pac (un-converted,
 * real contig names) from `prefix` as resident handles, beside the converted seed
 * FM-index: what records are scored, refined and reported against. Mirrors
 * indexEle::bwa_idx_load_ele's disk path (bns_restore, then slurp the full .pac
 * into memory and close fp_pac). On a bns_restore or pac-allocation failure,
 * frees any partial allocation, leaves the out-params NULL, and returns -1. (An
 * I/O error or short read during the .pac slurp itself is fatal --
 * pac_slurp_and_close aborts the process, matching the seed-index loader.)
 * `pread_workers` is the already-resolved worker count for the parallel .pac
 * slurp (as from index_load_threads); 1 reads serially. */
int meth_orig_ref_load(const char *prefix, bntseq_t **bns_out, uint8_t **pac_out, int pread_workers)
{
    *bns_out = NULL;
    *pac_out = NULL;
    bntseq_t *bns = bns_restore(prefix);
    if (bns == NULL) {
        fprintf(stderr,
                "ERROR: --meth could not load the original reference bns from "
                "'%s.{amb,ann,pac}'\n", prefix);
        return -1;
    }
    int64_t pac_bytes = bns->l_pac / 4 + 1;
    uint8_t *pac = (uint8_t*) calloc(pac_bytes, 1);
    if (pac == NULL) {
        fprintf(stderr, "ERROR: --meth failed to allocate %lld bytes for the "
                "original reference pac\n", (long long)pac_bytes);
        bns_destroy(bns);
        return -1;
    }
    /* bns_restore left .pac open in bns->fp_pac; slurp it whole (in parallel,
     * same as the seed index) and close it. Shared with bwa_idx_load_ele. */
    char pac_path[PATH_MAX];
    bns_pac_path(pac_path, sizeof(pac_path), prefix);
    pac_slurp_and_close(&bns->fp_pac, pac_path, pac, pac_bytes, pread_workers);

    *bns_out = bns;
    *pac_out = pac;
    return 0;
}

void meth_orig_ref_free(bntseq_t **bns, uint8_t **pac)
{
    if (*pac != NULL) { free(*pac); *pac = NULL; }
    if (*bns != NULL) { bns_destroy(*bns); *bns = NULL; }
}
