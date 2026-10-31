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

aln_out_t *aln_out_open(const char *path, int bam, int bam_level, const bntseq_t *bns, const char *hdr_text)
{
	(void)bns;
	aln_out_t *o = (aln_out_t*)calloc(1, sizeof(aln_out_t));
	if (!o) return NULL;
	if (!hdr_text) hdr_text = "";
	if (!bam) {
		if (path && strcmp(path, "-") != 0) {
			o->fp = fopen(path, "w");
			if (!o->fp) {
				fprintf(stderr, "[E::%s] cannot open '%s' for writing: %s\n", __func__, path, strerror(errno));
				free(o);
				return NULL;
			}
		} else {
			o->fp = stdout;
			o->is_stdout = 1;
		}
		if (fputs(hdr_text, o->fp) < 0) {
			fprintf(stderr, "[E::%s] cannot write the header\n", __func__);
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
	o->b = bam_init1();
	if (!o->hdr || !o->hf || !o->b || sam_hdr_write(o->hf, o->hdr) < 0) {
		fprintf(stderr, "[E::%s] cannot open '%s' for BAM output\n", __func__, path ? path : "-");
		if (o->hf) hts_close(o->hf);
		if (o->hdr) sam_hdr_destroy(o->hdr);
		if (o->b) bam_destroy1(o->b);
		free(o);
		return NULL;
	}
	return o;
}

int aln_out_write_buf(aln_out_t *o, const char *buf, size_t buf_l)
{
	if (!o || o->failed) return -1;
	if (!buf_l) return 0;
	if (o->fp) {
		if (fwrite(buf, 1, buf_l, o->fp) != buf_l) {
			fprintf(stderr, "[E::%s] write failed: %s\n", __func__, strerror(errno));
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
			fprintf(stderr, "[E::%s] cannot convert or write a BAM record\n", __func__);
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
