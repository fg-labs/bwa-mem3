/* Adapter: fr_fastq parser over a fast_reader, producing bseq1_t chunks with
 * the exact contract of bseq_read_orig (same chunking, parity, trim, per-field
 * allocation) — byte-identical to the former kseq path. The handle names below
 * are kept for the pipeline call sites; the opaque handle is an fr_fastq_t*. */
#ifndef FAST_READER_BSEQ_H
#define FAST_READER_BSEQ_H

#include "bwa.h"
#include "fast_reader.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Create a parser bound to an already-open fast_reader, returned as an opaque
 * handle. `fr` must be non-NULL and outlive the returned handle (the parser
 * borrows it; it does not take ownership and fast_kseq_destroy does NOT close
 * the underlying fast_reader). Never returns NULL on the caller-visible path
 * (it aborts on allocation failure). */
void *fast_kseq_init(fast_reader_t *fr);

/* Destroy a handle from fast_kseq_init. `p` may be NULL (no-op). Frees only the
 * parser buffers; the caller still owns and must close the fast_reader the
 * handle was bound to. */
void  fast_kseq_destroy(void *p);

/* Mirror of bseq_read_orig, reading from fast_reader-backed kseq handles.
 * Reads up to ~chunk_size bytes of sequence (cut on an even record boundary)
 * and returns a malloc'd bseq1_t array of *n_ records; the caller owns the
 * array and every per-record string field (name/comment/seq/qual). `ks1_` is
 * required; `ks2_` is the second mate handle and may be NULL for single-end /
 * interleaved input. At EOF the function sets *n_ = 0 and *s = 0; the returned
 * pointer is then NULL (no records were allocated). Mismatched record counts
 * between the two files are reported to stderr and truncate the batch.
 *
 * On success *arena_out receives the per-chunk bump arena backing the returned
 * reads' name/seq/qual fields; the caller owns it and must read_arena_destroy()
 * it after the chunk is fully consumed. At EOF (*n_ == 0, return NULL)
 * *arena_out is set to NULL. See read_arena.h for the ownership contract. */
bseq1_t *bseq_read_fast(int64_t chunk_size, int *n_, void *ks1_, void *ks2_, int64_t *s,
                        read_arena_t **arena_out, int copy_comment);

#ifdef __cplusplus
}
#endif

#endif /* FAST_READER_BSEQ_H */
