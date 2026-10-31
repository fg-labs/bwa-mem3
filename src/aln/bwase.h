/* Ported from bwa (https://github.com/lh3/bwa), Apache2 branch, commit a08c7c0.
 * Copyright (c) Heng Li and the bwa contributors; licensed under the Apache License,
 * Version 2.0 (see LICENSE in this directory).
 *
 * MODIFIED from the original: moved into namespace bwa_aln, SAM text goes to a
 * kstring_t instead of stdout, the .sai reading is dropped, and the BWT is a
 * parameter (bwa_cal_pac_pos no longer reloads it from disk). */
#ifndef BWA_MEM3_ALN_BWASE_H
#define BWA_MEM3_ALN_BWASE_H

#include "bwtaln.h"
#include "../bntseq.h"
#include "../kstring.h"

namespace bwa_aln {

/* Fills the table bwa_approx_mapQ reads. */
void bwase_initialize(void);
/* Collapses a read's .sai alignments into its best hit (set_main) and up to n_multi
 * extra hits. Consumes drand48(): call in read order, after srand48(bns->seed). */
void bwa_aln2seq_core(int n_aln, const bwt_aln1_t *aln, bwa_seq_t *s, int set_main, int n_multi);
void bwa_aln2seq(int n_aln, const bwt_aln1_t *aln, bwa_seq_t *s);
int bwa_approx_mapQ(const bwa_seq_t *p, int mm);
bwtint_t bwa_sa2pos(const bntseq_t *bns, const bwt_t *bwt, bwtint_t sapos, int ref_len, int *strand);
void bwa_cal_pac_pos_core(const bntseq_t *bns, const bwt_t *bwt, bwa_seq_t *seq, int max_mm, float fnr);
/* bwa_cal_pac_pos_core for every read, then the extra hits' positions. */
void bwa_cal_pac_pos(const bntseq_t *bns, const bwt_t *bwt, int n_seqs, bwa_seq_t *seqs, int max_mm, float fnr);
/* `pacseq` is the 2-bit packed forward reference (l_pac/4 + 1 bytes). Reverses each
 * read's seq in place, as upstream does. */
void bwa_refine_gapped(const bntseq_t *bns, int n_seqs, bwa_seq_t *seqs, ubyte_t *pacseq);
int64_t pos_end(const bwa_seq_t *p);
/* Appends the record(s) stock bwa samse prints for one read, newline included. */
void bwa_print_sam1(kstring_t *out, const bntseq_t *bns, bwa_seq_t *p, const bwa_seq_t *mate,
					int mode, int max_top2);
/* Appends the header stock bwa 0.7.19 samse prints: the default @HD (unless hdr_line has
 * one), the @SQ lines (unless hdr_line has its own), hdr_line, then pg_line. hdr_line and
 * pg_line may be NULL; pg_line gets a newline if it lacks one. */
void aln_print_sam_hdr(kstring_t *out, const bntseq_t *bns, const char *hdr_line, const char *pg_line);

} // namespace bwa_aln

#endif
