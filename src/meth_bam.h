/* SPDX-License-Identifier: MIT */
#ifndef BWAMEM3_METH_BAM_H
#define BWAMEM3_METH_BAM_H

#include <stdint.h>
#include "bwa.h"
#include "bwamem.h"
#include "bntseq.h"

/* htslib's <htslib/kstring.h> and bwa-mem3's src/kstring.h share the
 * `KSTRING_H` guard, so including htslib/sam.h here would make one of them
 * a no-op. htslib is therefore only included inside src/meth_bam.cpp; this
 * header forward-declares bam1_t and exposes htslib-free wrappers. */
struct bam1_t;

/* Native BAM emission for `bwa-mem3 mem --meth`: consolidates the
 * bwa-meth c2t doubled reference's fchr/rchr contigs into one @SQ per
 * real chromosome and converts mem_aln_t records to bam1_t with the
 * meth-specific tags (YD:Z) and chimera QC. */

#ifdef __cplusplus
extern "C" {
#endif

/* --- Chrom map: internal (fchr/rchr) -> output (chr) ---------------- */

typedef struct meth_chrom_map_s {
    int      n_internal;    /* matches bns->n_seqs */
    int      n_output;      /* after f/r collapse */
    int     *out_tid;       /* length n_internal; index into output_names */
    char    *direction;     /* length n_internal; 'f', 'r', or 0 */
    char   **output_names;  /* length n_output; stripped names (NUL-term) */
    int64_t *output_lens;   /* length n_output */
} meth_chrom_map_t;

meth_chrom_map_t *meth_chrom_map_build_from_bns(const bntseq_t *bns);
void              meth_chrom_map_free(meth_chrom_map_t *m);

/* --- BAM writer lifecycle ------------------------------------------- */

/* Opaque — defined in meth_bam.cpp so this header is htslib-free. */
typedef struct meth_bam_writer_s meth_bam_writer_t;

/* Global chrom map (set once by main_mem when --meth is active). */
extern meth_chrom_map_t *g_meth_cmap;

/* Global BAM writer (same lifecycle). */
extern meth_bam_writer_t *g_meth_bam_writer;

/* Open a BAM writer at path ("-" for stdout). `compression_level` is the
 * BGZF deflate level: 0 = uncompressed, 1..9 = deflate. Returns NULL on
 * failure. */
meth_bam_writer_t *meth_bam_writer_open(const char *path_or_dash,
                                        meth_chrom_map_t *cmap,
                                        const char *bwa_pg,
                                        const char *meth_pg_cl,
                                        int compression_level);

/* Write one bam1_t. Returns 0 on success, -1 on error. */
int meth_bam_writer_write(meth_bam_writer_t *w, struct bam1_t *b);

/* Close the writer and flush the BGZF EOF marker. Frees internal hdr and
 * htsFile; does NOT free the cmap. Returns 0 on success, -1 on error. */
int meth_bam_writer_close(meth_bam_writer_t *w);

/* --- mem_aln_t -> bam1_t --------------------------------------------- */

/* Convert one alignment to a bam1_t with meth transforms applied:
 *   - chrom rewrite: p->rid → cmap->out_tid[p->rid]
 *   - YD:Z:{f,r} tag from cmap->direction[p->rid]
 *   - chimera QC (unless opt->meth_no_chim): if longest M/=/X run < 44%
 *     of l_seq, set 0x200, clear 0x2, cap mapq at 1
 *   - opt->meth_set_as_failed forces 0x200 on the matching strand
 * Caller owns `b`. Returns 0 on success, -1 on error. */
int meth_mem_aln_to_bam(struct bam1_t *b,
                        const mem_opt_t *opt, const bntseq_t *bns,
                        const bseq1_t *s, int n_alns,
                        const mem_aln_t *list, int which,
                        const mem_aln_t *m_,
                        const meth_chrom_map_t *cmap);

/* If any record in `group` has 0x200 set, propagate it to the rest and
 * clear 0x2 on those. */
void meth_bam_group_propagate_qcfail(struct bam1_t **group, int n);

#ifdef __cplusplus
}
#endif

#endif
