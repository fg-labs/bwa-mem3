/* Adapter: kseq parser over a fast_reader, producing bseq1_t chunks with the
 * exact contract of bseq_read_orig (same chunking, parity, trim, per-field
 * allocation). Kept in its own translation unit because kseq.h can only be
 * instantiated for one handle type per TU. */
#ifndef FAST_READER_BSEQ_H
#define FAST_READER_BSEQ_H

#include "bwa.h"
#include "fast_reader.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Create a kseq bound to an already-open fast_reader, returned as an opaque
 * kseq_t*. `fr` must be non-NULL and outlive the returned handle (the kseq
 * borrows it; it does not take ownership and fast_kseq_destroy does NOT close
 * the underlying fast_reader). Returns the handle; never returns NULL on the
 * caller-visible path (the underlying kseq_init allocates with calloc and the
 * fast_reader is assumed valid). */
void *fast_kseq_init(fast_reader_t *fr);

/* Destroy a handle from fast_kseq_init. `ks` may be NULL (no-op). Frees only
 * the kseq buffers; the caller still owns and must close the fast_reader the
 * handle was bound to. */
void  fast_kseq_destroy(void *ks);

/* Mirror of bseq_read_orig, reading from fast_reader-backed kseq handles.
 * Reads up to ~chunk_size bytes of sequence (cut on an even record boundary)
 * and returns a malloc'd bseq1_t array of *n_ records; the caller owns the
 * array and every per-record string field (name/comment/seq/qual). `ks1_` is
 * required; `ks2_` is the second mate handle and may be NULL for single-end /
 * interleaved input. At EOF the function sets *n_ = 0 and *s = 0; the returned
 * pointer is then NULL (no records were allocated). Mismatched record counts
 * between the two files are reported to stderr and truncate the batch. */
bseq1_t *bseq_read_fast(int64_t chunk_size, int *n_, void *ks1_, void *ks2_, int64_t *s);

#ifdef __cplusplus
}
#endif

#endif /* FAST_READER_BSEQ_H */
