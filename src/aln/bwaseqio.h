/* Ported from bwa (https://github.com/lh3/bwa), Apache2 branch, commit a08c7c0.
 * Copyright (c) Heng Li and the bwa contributors; licensed under the Apache License,
 * Version 2.0 (see LICENSE in this directory).
 *
 * MODIFIED from the original: moved into namespace bwa_aln, BAM input is dropped,
 * and C++ casts are added; bwa_seq_close returns an error instead of exiting, and
 * bwa_seq_n_bad is new. */
#ifndef BWA_MEM3_ALN_BWASEQIO_H
#define BWA_MEM3_ALN_BWASEQIO_H

#include "bwtaln.h"

namespace bwa_aln {

/* "-" reads standard input; gzip-compressed input is handled by zlib. */
bwa_seqio_t *bwa_seq_open(const char *fn);
/* Returns 0, or -1 (after a message) if closing the input failed; frees bs either way. */
int bwa_seq_close(bwa_seqio_t *bs);
/* The number of malformed records bwa_read_seq has skipped (each named on stderr). */
int bwa_seq_n_bad(const bwa_seqio_t *bs);
void seq_reverse(int len, ubyte_t *seq, int is_comp);
int bwa_trim_read(int trim_qual, bwa_seq_t *p);
bwa_seq_t *bwa_read_seq(bwa_seqio_t *bs, int n_needed, int *n, int mode, int trim_qual);
void bwa_free_read_seq(int n_seqs, bwa_seq_t *seqs);

} // namespace bwa_aln

#endif
