/* fr_fastq: a single-copy FASTQ/FASTA parser over a fast_reader byte source.
 *
 * Replaces the kseq layer on the read hot path. kseq copies every field twice
 * (codec buffer -> kstring_t, then kstring_t -> bseq1_t); this parser slices
 * records directly out of its own refillable buffer so the caller copies each
 * field exactly once (or, for #2, defers the copy to a worker thread).
 *
 * The grammar matches kseq's exactly so output is byte-identical:
 *   - record starts at '>' (FASTA) or '@' (FASTQ); leading junk before the
 *     first header is skipped, as are blank lines between records;
 *   - name = bytes up to the first whitespace; comment = the rest of that line
 *     (empty if the name ran to end-of-line), with a single trailing '\r'
 *     dropped when the line has more than one byte;
 *   - seq = the sequence line(s), newlines removed, up to the next header or
 *     the '+' separator; CRLF endings are normalised the same way kseq does;
 *   - qual (FASTQ only) = quality line(s) accumulated until they match the
 *     sequence length.
 *
 * Slice pointers returned by fr_fastq_next() reference parser-internal storage
 * and are valid only until the next fr_fastq_next() / fr_fastq_destroy() call.
 */
#ifndef FR_FASTQ_H
#define FR_FASTQ_H

#include <stddef.h>

#include "fast_reader.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct fr_fastq fr_fastq_t;

/* One parsed record. Lengths are exact; pointers are NOT guaranteed to be
 * NUL-terminated. comment_l == 0 means "no comment" (mirrors kseq's
 * comment.l == 0). qual_l == 0 means no quality (FASTA, or an empty record);
 * downstream should treat comment/qual as absent when their length is 0, which
 * reproduces fr_kseq2bseq1's `comment.l ? ... : 0` / `qual.l ? ... : 0`. */
typedef struct {
    const char *name;    size_t name_l;
    const char *comment; size_t comment_l;
    const char *seq;     size_t seq_l;
    const char *qual;    size_t qual_l;
} fr_fastq_rec_t;

/* Bind a parser to an already-open fast_reader. `fr` is borrowed: it must
 * outlive the parser and is NOT closed by fr_fastq_destroy. Aborts on OOM. */
fr_fastq_t *fr_fastq_init(fast_reader_t *fr);

/* Free the parser's buffers. `p` may be NULL. Does not close the fast_reader. */
void fr_fastq_destroy(fr_fastq_t *p);

/* Parse the next record into *rec.
 *   returns 1  -> a record was parsed; slices in *rec are valid until the next call
 *   returns 0  -> clean end of input
 *   returns -2 -> either a malformed record (matches kseq_read's -2: a '+' line
 *                 with no quality, or a quality string whose length never reaches
 *                 the sequence length) OR an upstream decode/read failure from the
 *                 underlying fast_reader. The caller cannot distinguish the two and
 *                 stops the batch in both cases, as the kseq path did via
 *                 `kseq_read(...) >= 0`. */
int fr_fastq_next(fr_fastq_t *p, fr_fastq_rec_t *rec);

#ifdef __cplusplus
}
#endif

#endif /* FR_FASTQ_H */
