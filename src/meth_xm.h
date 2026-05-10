/* SPDX-License-Identifier: MIT */
#ifndef BWAMEM3_METH_XM_H
#define BWAMEM3_METH_XM_H

#include <stdint.h>

#include "meth_orig_ref.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Build a Bismark XM:Z payload (length l_emit, NUL-terminated). NULL on
 * alloc failure.
 *
 * Ownership: the returned pointer aliases a thread-local growable scratch
 * buffer owned by the helper. Caller MUST NOT free it. The pointer is
 * valid until the next call to meth_build_xm from the same thread. Pass
 * the bytes directly to bam_aux_append (which copies them into the BAM
 * record) or copy elsewhere before the next call.
 *

 *   o            : un-converted reference (forward-strand bases per real chrom)
 *   real_tid     : index into o's per-chrom array (== cmap->out_tid[p->rid])
 *   pos          : 0-based forward-genome alignment start
 *   is_top_strand: 1 if methylation events are on the top (forward) strand,
 *                  0 if on the bottom strand. This is determined by the
 *                  XG:Z value, NOT by the SAM 0x10 flag: XG:Z:CT (= cmap
 *                  direction 'f') => top strand; XG:Z:GA (= 'r') => bottom.
 *   bam_cigar    : CIGAR in BAM op encoding (0=M, 1=I, 2=D, 3=N, 4=S, 5=H, 6=P, 7==, 8=X)
 *   n_cigar      : number of cigar ops
 *   seq_text     : SEQ-orientation read bases (ASCII A/C/G/T), length l_emit
 *   l_emit       : SEQ length emitted to BAM (post supp soft-clip trim)
 *
 * Output XM matches Bismark's methylation_call: per-base char in SEQ
 * orientation, with z/Z (CpG), x/X (CHG), h/H (CHH), u/U (unknown), '.' for
 * everything else. */
char *meth_build_xm(const meth_orig_ref_t *o, int real_tid,
                    int64_t pos, int is_top_strand,
                    const uint32_t *bam_cigar, int n_cigar,
                    const char *seq_text, int l_emit);

#ifdef __cplusplus
}
#endif

#endif
