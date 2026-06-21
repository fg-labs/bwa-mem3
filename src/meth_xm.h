/* SPDX-License-Identifier: MIT */
#ifndef BWAMEM3_METH_XM_H
#define BWAMEM3_METH_XM_H

#include <stdint.h>

#include "bntseq.h"

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

 *   bns, pac     : the ORIGINAL (un-converted) reference bns + 2-bit pac
 *                  (D3 PR-5: replaces the retired meth_orig_ref f/r fold).
 *                  Forward-strand bases are decoded inline per record.
 *   real_tid     : contig index into bns->anns (== the alignment's original rid)
 *   pos          : 0-based forward-genome alignment start (contig-local)
 *   is_top_strand: 1 if methylation events are on the top (forward) strand,
 *                  0 if on the bottom strand. Sourced from the winning
 *                  hypothesis (= XG:Z), NOT the SAM 0x10 flag: OT (XG:Z:CT)
 *                  => top strand; OB (XG:Z:GA) => bottom.
 *   bam_cigar    : CIGAR in BAM op encoding (0=M, 1=I, 2=D, 3=N, 4=S, 5=H, 6=P, 7==, 8=X)
 *   n_cigar      : number of cigar ops
 *   seq_text     : SEQ-orientation read bases (ASCII A/C/G/T), length l_emit
 *   l_emit       : SEQ length emitted to BAM (post supp soft-clip trim)
 *
 * Output XM matches Bismark's methylation_call: per-base char in SEQ
 * orientation, with z/Z (CpG), x/X (CHG), h/H (CHH), u/U (unknown), '.' for
 * everything else. */
char *meth_build_xm(const bntseq_t *bns, const uint8_t *pac, int real_tid,
                    int64_t pos, int is_top_strand,
                    const uint32_t *bam_cigar, int n_cigar,
                    const char *seq_text, int l_emit);

#ifdef __cplusplus
}
#endif

#endif
