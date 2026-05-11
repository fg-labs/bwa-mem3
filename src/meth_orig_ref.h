/* SPDX-License-Identifier: MIT */
#ifndef BWAMEM3_METH_ORIG_REF_H
#define BWAMEM3_METH_ORIG_REF_H

#include <stdint.h>

#include "bntseq.h"
#include "meth_bam.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Lazy un-converted reference for `bwa-mem3 mem --meth`. Holds borrowed
 * references to the doubled c2t bns and pac (already loaded for alignment)
 * plus a cmap-derived sibling rid table used to find the f/r contig pair
 * for each real chromosome. Read-only after build; safe to share across
 * worker threads.
 *
 * Per-record slice does dual-pac decode: read the f-contig and r-contig
 * windows at the same forward-genome offset from the doubled pac, fold
 * via a 5-row (f, r) -> original table, mask N positions via
 * bns_iter_ambi over the doubled bns->ambs. No separate un-converted pac
 * is materialized — reuses idx->pac in place.
 *
 * Memory cost: zero new pac storage; two int* sibling-rid arrays of
 * length cmap->n_output (a few KB on hg38). */
typedef struct meth_orig_ref_s meth_orig_ref_t;

/* Build from an already-loaded doubled-c2t bns + pac (same handles
 * passed to the rest of bwa-mem3). cmap is consulted to find each
 * real chrom's f-rid and r-rid. Returns NULL on alloc failure or if
 * any real chrom is missing one direction (defensive: bwa-mem3 index
 * --meth always emits both f-<chr> and r-<chr>). */
meth_orig_ref_t *meth_orig_ref_load(const bntseq_t *bns, const uint8_t *pac,
                                    const meth_chrom_map_t *cmap);

/* Slice un-converted forward-strand bases [st, en) of chrom `real_tid`
 * into caller-owned `dst` (capacity must be >= en - st). Encoding is
 * ASCII ('A'/'C'/'G'/'T'/'N'). OOB positions and bns->ambs overlaps
 * fill with 'N'. `real_tid` indexes cmap->output_names. */
void meth_orig_ref_slice(const meth_orig_ref_t *o,
                         int real_tid, int64_t st, int64_t en,
                         uint8_t *dst);

void meth_orig_ref_free(meth_orig_ref_t *o);

/* Global instance, populated once in fastmap.cpp under --meth. */
extern meth_orig_ref_t *g_meth_orig_ref;

#ifdef __cplusplus
}
#endif

#endif
