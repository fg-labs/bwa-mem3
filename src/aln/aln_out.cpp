/* SPDX-License-Identifier: MIT */
/* New code for bwa-mem3 (not from upstream bwa): the output sink for `bwa-mem3 aln`. */

/* htslib headers before any bwa-mem3 header that pulls in kstring.h (they share the
 * KSTRING_H include guard); the two kstring_t types have the same layout. */
#include "htslib/sam.h"
#include "htslib/kstring.h"

#include "aln_out.h"

#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

namespace bwa_aln {

struct aln_out_t {
	FILE *fp;       // SAM
	int is_stdout;
	htsFile *hf;    // BAM
	sam_hdr_t *hdr;
	bam1_t *b;
	::kstring_t line;
	int failed;
};

aln_out_t *aln_out_open(const char *path, int bam, int bam_level, const char *hdr_text, int n_threads)
{
	aln_out_t *o = (aln_out_t*)calloc(1, sizeof(aln_out_t));
	if (!o) return NULL;
	if (!hdr_text) hdr_text = "";
	if (!bam) {
		if (path && strcmp(path, "-") != 0) {
			o->fp = fopen(path, "w");
			if (!o->fp) {
				fprintf(stderr, "[E::aln] cannot open '%s' for writing: %s\n", path, strerror(errno));
				free(o);
				return NULL;
			}
		} else {
			o->fp = stdout;
			o->is_stdout = 1;
		}
		if (fputs(hdr_text, o->fp) < 0) {
			fprintf(stderr, "[E::aln] cannot write the header\n");
			if (!o->is_stdout) fclose(o->fp);
			free(o);
			return NULL;
		}
		return o;
	}
	char mode[8];
	if (bam_level < 0) bam_level = 0;
	if (bam_level > 9) bam_level = 9;
	snprintf(mode, sizeof mode, "wb%d", bam_level);
	o->hdr = sam_hdr_parse(strlen(hdr_text), hdr_text);
	o->hf = o->hdr ? hts_open(path ? path : "-", mode) : NULL;
	/* BGZF deflate on htslib's thread pool, as bam_writer does for `mem`: its ordered pool
	 * writes the blocks in order, so the bytes are those of the serial path. */
	if (o->hf && n_threads > 1 && hts_set_threads(o->hf, n_threads) < 0)
		fprintf(stderr, "[W::aln] could not start %d BGZF threads; compressing on one thread\n", n_threads);
	o->b = bam_init1();
	if (!o->hdr || !o->hf || !o->b || sam_hdr_write(o->hf, o->hdr) < 0) {
		fprintf(stderr, "[E::aln] cannot open '%s' for BAM output\n", path ? path : "-");
		if (o->hf) hts_close(o->hf);
		if (o->hdr) sam_hdr_destroy(o->hdr);
		if (o->b) bam_destroy1(o->b);
		free(o);
		return NULL;
	}
	return o;
}

/* The `i`th tab-separated field (from 0) of the line [p, p + n); sets *len. */
static const char *sam_field(const char *p, size_t n, int i, int *len)
{
	const char *end = p + n;
	for (; i > 0 && p < end; --i) {
		const char *t = (const char*)memchr(p, '\t', end - p);
		if (!t) return NULL;
		p = t + 1;
	}
	const char *t = (const char*)memchr(p, '\t', end - p);
	*len = (int)((t ? t : end) - p);
	return p;
}

/* sam_parse1 turns a reference name that is not in the header into tid -1 and marks the record
 * unmapped, with only a warning. That happens when the -r text carries its own @SQ lines (they
 * replace the index's, as in stock samse), and would silently lose every record's reference, so
 * it is an error here. Returns the unknown name's field or NULL. */
static const char *unknown_ref(const bam1_t *b, const char *p, size_t n, int *len)
{
	const char *f;
	if (b->core.tid < 0 && (f = sam_field(p, n, 2, len)) && !(*len == 1 && f[0] == '*')) return f;
	if (b->core.mtid < 0 && (f = sam_field(p, n, 6, len)) && !(*len == 1 && (f[0] == '*' || f[0] == '='))) return f;
	return NULL;
}

int aln_out_write_buf(aln_out_t *o, const char *buf, size_t buf_l)
{
	if (!o || o->failed) return -1;
	if (!buf_l) return 0;
	if (o->fp) {
		if (fwrite(buf, 1, buf_l, o->fp) != buf_l) {
			fprintf(stderr, "[E::aln] write failed: %s\n", strerror(errno));
			o->failed = 1;
			return -1;
		}
		return 0;
	}
	const char *p = buf, *end = buf + buf_l;
	while (p < end) {
		const char *nl = (const char*)memchr(p, '\n', end - p);
		size_t n = nl ? (size_t)(nl - p) : (size_t)(end - p);
		o->line.l = 0;
		kputsn(p, n, &o->line);
		if (sam_parse1(&o->line, o->hdr, o->b) < 0 || sam_write1(o->hf, o->hdr, o->b) < 0) {
			const char *tab = (const char*)memchr(p, '\t', n);
			fprintf(stderr, "[E::aln] cannot convert or write the BAM record of read '%.*s' (a read containing '-' cannot be written as BAM)\n",
					(int)(tab ? tab - p : n), p);
			o->failed = 1;
			return -1;
		}
		int ref_l = 0;
		const char *ref = unknown_ref(o->b, p, n, &ref_l);
		if (ref) {
			const char *tab = (const char*)memchr(p, '\t', n);
			fprintf(stderr, "[E::aln] read '%.*s': reference '%.*s' is not in the BAM header (an -r text with its own @SQ lines replaces the index's)\n",
					(int)(tab ? tab - p : n), p, ref_l, ref);
			o->failed = 1;
			return -1;
		}
		p += n + (nl ? 1 : 0);
	}
	return 0;
}

int aln_out_close(aln_out_t *o)
{
	int ret = 0;
	if (!o) return -1;
	if (o->fp) {
		if (fflush(o->fp) != 0) ret = -1;
		if (!o->is_stdout && fclose(o->fp) != 0) ret = -1;
	} else {
		if (hts_close(o->hf) != 0) ret = -1;
		sam_hdr_destroy(o->hdr);
		bam_destroy1(o->b);
		free(o->line.s);
	}
	if (o->failed) ret = -1;
	free(o);
	return ret;
}

} // namespace bwa_aln
