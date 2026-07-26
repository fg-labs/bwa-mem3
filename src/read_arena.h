/* Per-chunk bump arena for immutable per-read string fields (PIPE-F6).
 *
 * Motivation: the FASTQ readers used to malloc() every per-read string field
 * (name/seq/qual) individually at ingest, and the output stage free()d each one
 * individually — frequently on a DIFFERENT thread than the one that allocated
 * it (the klib kt_pipeline reads a chunk on one thread and writes it on
 * another). That is ~3 mallocs/read + ~3 frees/read of small buffers, with the
 * frees crossing threads, which stresses the allocator's cross-thread free
 * path across the ~100M reads of a typical run.
 *
 * This arena replaces that with one growable, per-chunk bump allocator: the
 * reader carves each read's name/seq/qual out of the arena, and the output
 * stage frees the WHOLE arena once per chunk instead of per field. Only the
 * allocation strategy changes — field contents, lengths, and NUL-termination
 * are byte-identical to the former per-field malloc/strdup path.
 *
 * Ownership / lifetime contract:
 *   - Exactly one arena per COHORT -- the set of reads whose fields are consumed
 *     together, which is one mem_pestat batch. It is created by the reader on the
 *     first read of the cohort, returned to the pipeline alongside the bseq1_t
 *     array, stored on the chunk (ktp_data_t.read_arena), and destroyed exactly
 *     once, in the write stage, after every field it backs has been consumed.
 *   - A cohort may be READ in several calls (cohort slicing exists so compute can
 *     start before the whole batch has arrived). The `read_arena_t **arena_out`
 *     parameter is therefore in/out: NULL in means "create a fresh arena",
 *     non-NULL means "keep carving from this one". Only the call that CREATED an
 *     arena may destroy it on the empty-batch path; a carried arena still backs
 *     the earlier slices. This is why the unit is the cohort and not the read
 *     call -- the fields outlive the call that carved them.
 *   - The arena is owned by a single cohort and never shared mutably between the
 *     threads processing different cohorts (each carries its own).
 *   - Only allocation grows the arena, and that happens solely on the read
 *     thread during ingest; the process/write stages only READ the fields (the
 *     in-place 2-bit sequence encoding and --meth C->T/G->A projection mutate
 *     the seq bytes in place, same length, which is fine — arena memory is
 *     mutable and those never reallocate the field).
 *
 * Pointer stability: the arena is SEGMENTED (a linked list of blocks). A new
 * allocation that does not fit the current block gets a fresh block; existing
 * blocks are never realloc()'d or moved, so every pointer previously handed out
 * stays valid for the arena's lifetime. This is the property that lets the
 * reader carve many fields and hand out stable char* into the same arena. */
#ifndef READ_ARENA_H
#define READ_ARENA_H

#include <stddef.h>
#include <stdlib.h>
#include <string.h>

#include "utils.h" /* err_fatal — abort loudly on OOM, matching house style */

#ifdef __cplusplus
extern "C" {
#endif

/* Default block size. Read fields are small (tens–hundreds of bytes) and a
 * chunk holds up to ~task_size bytes of sequence, so 1 MiB blocks keep the
 * block count (and thus the malloc/free count) tiny while wasting at most the
 * unused tail of one block per fill. */
#define READ_ARENA_BLOCK_SIZE ((size_t)(1u << 20))

typedef struct read_arena_block_s {
	struct read_arena_block_s *next; /* previous (older) block; NULL terminates */
	char  *data;                     /* block storage */
	size_t cap;                      /* bytes allocated for `data` */
	size_t used;                     /* bytes handed out from this block */
} read_arena_block_t;

typedef struct read_arena_s {
	read_arena_block_t *head; /* current (newest) block; NULL until first alloc */
} read_arena_t;

/* Create an empty arena. Aborts on OOM. */
static inline read_arena_t *read_arena_create(void)
{
	read_arena_t *a = (read_arena_t *)malloc(sizeof(read_arena_t));
	if (a == NULL)
		err_fatal(__func__, "failed to allocate read arena");
	a->head = NULL;
	return a;
}

/* Push a fresh block large enough for at least `need` bytes to the front. */
static inline void read_arena_grow(read_arena_t *a, size_t need)
{
	size_t cap = READ_ARENA_BLOCK_SIZE;
	if (need > cap) cap = need;
	read_arena_block_t *b = (read_arena_block_t *)malloc(sizeof(read_arena_block_t));
	if (b == NULL)
		err_fatal(__func__, "failed to allocate read arena block header");
	b->data = (char *)malloc(cap);
	if (b->data == NULL)
		err_fatal(__func__, "failed to allocate %zu-byte read arena block", cap);
	b->cap = cap;
	b->used = 0;
	b->next = a->head;
	a->head = b;
}

/* Carve `n` bytes from the arena. Returned pointer is stable for the arena's
 * lifetime. Byte alignment is sufficient here — the only callers store NUL-
 * terminated char strings (name/seq/qual), which have no alignment
 * requirement. Aborts on OOM. n == 0 is legal and returns a valid pointer. */
static inline void *read_arena_alloc(read_arena_t *a, size_t n)
{
	read_arena_block_t *b = a->head;
	if (b == NULL || b->used + n > b->cap) {
		read_arena_grow(a, n);
		b = a->head;
	}
	void *p = b->data + b->used;
	b->used += n;
	return p;
}

/* Copy a length-known field into the arena and NUL-terminate it, returning the
 * arena-owned copy. The source need not be NUL-terminated; dst[len] is set
 * explicitly. Produces a byte-identical result to malloc(len+1)+memcpy. */
static inline char *read_arena_dup(read_arena_t *a, const char *src, size_t len)
{
	char *dst = (char *)read_arena_alloc(a, len + 1);
	memcpy(dst, src, len);
	dst[len] = '\0';
	return dst;
}

/* Free the arena and every block it owns. NULL is a no-op. */
static inline void read_arena_destroy(read_arena_t *a)
{
	if (a == NULL) return;
	read_arena_block_t *b = a->head;
	while (b != NULL) {
		read_arena_block_t *prev = b->next;
		free(b->data);
		free(b);
		b = prev;
	}
	free(a);
}

#ifdef __cplusplus
}
#endif

#endif /* READ_ARENA_H */
