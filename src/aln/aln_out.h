/* Output sink for `bwa-mem3 aln`: SAM text or BAM. New code (not from upstream bwa). */
#ifndef BWA_MEM3_ALN_OUT_H
#define BWA_MEM3_ALN_OUT_H

#include "../kstring.h"

namespace bwa_aln {

struct aln_out_t;

/* Opens `path` (NULL = standard output; "-" too for BAM) and writes the header. hdr_text
 * is the complete SAM header (newline-terminated lines, from aln_print_sam_hdr). For BAM
 * it is parsed with htslib; bam_level is the BGZF level 0-9. Returns NULL on failure
 * after printing why. */
aln_out_t *aln_out_open(const char *path, int bam, int bam_level, const char *hdr_text);
/* Writes complete SAM lines (each ending in a newline). Returns 0, or -1 on failure.
 * The buffer form is the real function: aln_out.cpp sees htslib's kstring_t, every other
 * file bwa-mem3's, and the two have different C++ names. */
int aln_out_write_buf(aln_out_t *o, const char *s, size_t l);
inline int aln_out_write(aln_out_t *o, const kstring_t *sam_lines)
{
	return aln_out_write_buf(o, sam_lines->s, sam_lines->l);
}
/* Flushes and closes, and frees `o`. Returns 0, or -1 if anything failed. */
int aln_out_close(aln_out_t *o);

} // namespace bwa_aln

#endif
