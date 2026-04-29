/* SPDX-License-Identifier: MIT */

/* htslib headers must come before any bwa-mem3 header that pulls in
 * bwa-mem3's kstring.h (they share the KSTRING_H include guard). */
#include "htslib/sam.h"
#include "htslib/kstring.h"

#include "meth_bam.h"
#include "bam_writer.h"
#include "cigar_util.h"
#include "version.h"

#include <cstdio>
#include <cstdlib>
#include <cstring>

/* Private writer struct — opaque to callers via meth_bam.h. */
struct meth_bam_writer_s {
    htsFile          *fp;
    sam_hdr_t        *hdr;
    meth_chrom_map_t *cmap; /* non-owning */
};

/* Declared in src/bwa.h (without extern "C"); match that linkage here. */
extern char bwa_rg_id[256];

meth_bam_writer_t *g_meth_bam_writer = NULL;
/* g_meth_cmap lives in bwamem.cpp (so the worker hook can reach it without
 * a link dependency on meth_bam.cpp; both files see the header). */

/* ------------------------------------------------------------------- */
/* Chrom map                                                            */
/* ------------------------------------------------------------------- */

meth_chrom_map_t *meth_chrom_map_build_from_bns(const bntseq_t *bns)
{
    if (bns == NULL) return NULL;
    meth_chrom_map_t *m = (meth_chrom_map_t *)calloc(1, sizeof(*m));
    if (m == NULL) return NULL;
    m->n_internal = bns->n_seqs;
    if (m->n_internal > 0) {
        m->out_tid      = (int *)    calloc((size_t)m->n_internal, sizeof(int));
        m->direction    = (char *)   calloc((size_t)m->n_internal, sizeof(char));
        m->output_names = (char **)  calloc((size_t)m->n_internal, sizeof(char *));
        m->output_lens  = (int64_t *)calloc((size_t)m->n_internal, sizeof(int64_t));
        if (!m->out_tid || !m->direction || !m->output_names || !m->output_lens) {
            meth_chrom_map_free(m);
            return NULL;
        }
    }

    for (int i = 0; i < m->n_internal; ++i) {
        const char *name = bns->anns[i].name;
        char prefix = name[0];
        const char *stripped = name;
        char dir = 0;
        if (prefix == 'f' || prefix == 'r') {
            stripped = name + 1;
            dir = prefix;
        }
        m->direction[i] = dir;

        int existing = -1;
        for (int j = 0; j < m->n_output; ++j) {
            if (strcmp(m->output_names[j], stripped) == 0) { existing = j; break; }
        }
        if (existing >= 0) {
            m->out_tid[i] = existing;
        } else {
            int idx = m->n_output++;
            m->output_names[idx] = strdup(stripped);
            if (m->output_names[idx] == NULL) { meth_chrom_map_free(m); return NULL; }
            m->output_lens[idx] = bns->anns[i].len;
            m->out_tid[i] = idx;
        }
    }
    return m;
}

void meth_chrom_map_free(meth_chrom_map_t *m)
{
    if (m == NULL) return;
    free(m->out_tid);
    free(m->direction);
    if (m->output_names != NULL) {
        for (int i = 0; i < m->n_output; ++i) free(m->output_names[i]);
        free(m->output_names);
    }
    free(m->output_lens);
    free(m);
}

/* ------------------------------------------------------------------- */
/* BAM writer                                                           */
/* ------------------------------------------------------------------- */

/* Replace embedded tabs in a string with spaces (BAM CL tag fields can't
 * contain literal tabs — safer to just sanitize). Caller's buffer. */
static void sanitize_cl(char *s)
{
    for (; *s; ++s) if (*s == '\t') *s = ' ';
}

meth_bam_writer_t *meth_bam_writer_open(const char *path_or_dash,
                                        meth_chrom_map_t *cmap,
                                        const char *bwa_pg,
                                        const char *meth_pg_cl,
                                        int compression_level)
{
    if (path_or_dash == NULL || cmap == NULL) return NULL;
    meth_bam_writer_t *w = (meth_bam_writer_t *)calloc(1, sizeof(*w));
    if (w == NULL) return NULL;
    w->cmap = cmap;

    if (compression_level < 0) compression_level = 0;
    if (compression_level > 9) compression_level = 9;
    char mode[8];
    snprintf(mode, sizeof(mode), "wb%d", compression_level);
    w->fp = hts_open(path_or_dash, mode);
    if (w->fp == NULL) { free(w); return NULL; }

    w->hdr = sam_hdr_init();
    if (w->hdr == NULL) { hts_close(w->fp); free(w); return NULL; }

    /* @HD */
    if (sam_hdr_add_line(w->hdr, "HD", "VN", "1.6", "SO", "unsorted", NULL) < 0) goto fail;

    /* @SQ from consolidated chrom map */
    for (int i = 0; i < cmap->n_output; ++i) {
        char len_buf[32];
        snprintf(len_buf, sizeof(len_buf), "%lld", (long long)cmap->output_lens[i]);
        if (sam_hdr_add_line(w->hdr, "SQ",
                             "SN", cmap->output_names[i],
                             "LN", len_buf,
                             NULL) < 0) goto fail;
    }

    /* Original bwa-mem3 @PG */
    if (bwa_pg != NULL && bwa_pg[0] != '\0') {
        if (sam_hdr_add_lines(w->hdr, bwa_pg, 0) < 0) goto fail;
    }

    /* bwa-mem3-meth @PG — grow as needed so very long CLs aren't truncated. */
    {
        const char *cl_in = (meth_pg_cl && meth_pg_cl[0]) ? meth_pg_cl
                                                          : "bwa-mem3 mem --meth";
        char *cl_copy = strdup(cl_in);
        if (cl_copy == NULL) goto fail;
        sanitize_cl(cl_copy);
        kstring_t pg = {0, 0, NULL};
        ksprintf(&pg, "@PG\tID:bwa-mem3-meth\tPN:bwa-mem3-meth\tVN:%s-meth\tCL:%s\n",
                 PACKAGE_VERSION, cl_copy);
        free(cl_copy);
        int rc = (pg.s != NULL) ? sam_hdr_add_lines(w->hdr, pg.s, 0) : -1;
        free(pg.s);
        if (rc < 0) goto fail;
    }

    if (sam_hdr_write(w->fp, w->hdr) < 0) goto fail;
    return w;

fail:
    sam_hdr_destroy(w->hdr);
    hts_close(w->fp);
    free(w);
    return NULL;
}

int meth_bam_writer_write(meth_bam_writer_t *w, bam1_t *b)
{
    if (w == NULL || b == NULL) return -1;
    return sam_write1(w->fp, w->hdr, b);
}

int meth_bam_writer_close(meth_bam_writer_t *w)
{
    if (w == NULL) return 0;
    int rc = 0;
    if (w->hdr) { sam_hdr_destroy(w->hdr); w->hdr = NULL; }
    if (w->fp) {
        if (hts_close(w->fp) < 0) rc = -1;
        w->fp = NULL;
    }
    free(w);
    return rc;
}

/* ------------------------------------------------------------------- */
/* mem_aln_t -> bam1_t                                                  */
/* ------------------------------------------------------------------- */

/* bwameth.py's chimera heuristic: flag 0x200 + cap mapq when the longest
 * alignment run covers < MIN_LONGEST_M_PCT% of the read. */
static constexpr int MIN_LONGEST_M_PCT = 44;

int meth_mem_aln_to_bam(bam1_t *b,
                        const mem_opt_t *opt, const bntseq_t *bns,
                        const bseq1_t *s, int n_alns,
                        const mem_aln_t *list, int which,
                        const mem_aln_t *m_,
                        const meth_chrom_map_t *cmap)
{
    if (b == NULL || opt == NULL || s == NULL || list == NULL || cmap == NULL) return -1;

    /* Local copies so flag/rid can be mutated without touching the caller's. */
    mem_aln_t p = list[which];
    mem_aln_t m;
    const mem_aln_t *mp = NULL;
    if (m_ != NULL) { m = *m_; mp = &m; }

    p.flag |= mp ? 0x1 : 0;
    p.flag |= p.rid < 0 ? 0x4 : 0;
    p.flag |= mp && mp->rid < 0 ? 0x8 : 0;
    if (p.rid < 0 && mp && mp->rid >= 0) {
        p.rid = mp->rid; p.pos = mp->pos; p.is_rev = mp->is_rev; p.n_cigar = 0;
    }
    if (mp && mp->rid < 0 && p.rid >= 0) {
        m.rid = p.rid; m.pos = p.pos; m.is_rev = p.is_rev; m.n_cigar = 0;
    }
    p.flag |= p.is_rev ? 0x10 : 0;
    p.flag |= mp && mp->is_rev ? 0x20 : 0;

    /* Fold bwa-mem3's high-bit supp flag (0x10000) down into BAM's 0x100. */
    uint16_t flag16 = (uint16_t)((p.flag & 0xffff) | (p.flag & 0x10000 ? 0x100 : 0));

    /* Direction (YD:Z source) comes from the chrom name prefix ('f'/'r')
     * written by `bwa-mem3 index --meth`. */
    int32_t tid = -1, mtid = -1;
    char direction = 0;
    if (p.rid >= 0 && p.rid < cmap->n_internal) {
        tid       = cmap->out_tid[p.rid];
        direction = cmap->direction[p.rid];
    }
    if (mp && mp->rid >= 0 && mp->rid < cmap->n_internal) {
        mtid = cmap->out_tid[mp->rid];
    }

    /* Remap primary CIGAR: bwa-mem3 ops -> BAM ops, + soft->hard for supp */
    uint32_t *bam_cigar = NULL;
    size_t    bam_n_cigar = 0;
    if (p.n_cigar > 0) {
        bam_cigar = (uint32_t *)malloc((size_t)p.n_cigar * sizeof(uint32_t));
        if (bam_cigar == NULL) return -1;
        for (int i = 0; i < p.n_cigar; ++i) {
            int op  = p.cigar[i] & 0xf;
            int len = p.cigar[i] >> 4;
            if (!(opt->flag & MEM_F_SOFTCLIP) && !p.is_alt && (op == 3 || op == 4))
                op = which ? 4 : 3;              /* hard clip for supp */
            uint32_t bam_op = (op >= 0 && op < 5) ? BAM_OP_FROM_MEM[op] : 0;
            bam_cigar[i] = ((uint32_t)len << 4) | bam_op;
        }
        bam_n_cigar = (size_t)p.n_cigar;
    }

    /* TLEN is emitted on the consolidated (post-remap) contigs, so it must be
     * gated on tid/mtid equality — not p.rid/mp->rid. Mates that rescue onto
     * opposite projected strands (f* vs r*) of the same real chromosome have
     * different internal rids but the same output tid, and SAM expects a real
     * TLEN when RNAME==RNEXT. */
    hts_pos_t tlen = 0;
    if (mp && tid >= 0 && mtid >= 0 && tid == mtid
        && p.n_cigar > 0 && m.n_cigar > 0) {
        int64_t p_rlen = cigar_ref_len_mem(p.cigar, p.n_cigar);
        int64_t m_rlen = cigar_ref_len_mem(m.cigar, m.n_cigar);
        int64_t p0 = p.pos + (p.is_rev ? p_rlen - 1 : 0);
        int64_t p1 = m.pos + (m.is_rev ? m_rlen - 1 : 0);
        tlen = -(p0 - p1 + (p0 > p1 ? 1 : p0 < p1 ? -1 : 0));
    }

    /* Compute SEQ/QUAL range with supp soft-clip trim */
    int qb = 0, qe = s->l_seq;
    if (p.n_cigar && which && !(opt->flag & MEM_F_SOFTCLIP) && !p.is_alt) {
        if (!p.is_rev) {
            int c0 = p.cigar[0] & 0xf;
            int cN = p.cigar[p.n_cigar-1] & 0xf;
            if (c0 == 3 || c0 == 4) qb += p.cigar[0] >> 4;
            if (cN == 3 || cN == 4) qe -= p.cigar[p.n_cigar-1] >> 4;
        } else {
            int c0 = p.cigar[0] & 0xf;
            int cN = p.cigar[p.n_cigar-1] & 0xf;
            if (c0 == 3 || c0 == 4) qe -= p.cigar[0] >> 4;
            if (cN == 3 || cN == 4) qb += p.cigar[p.n_cigar-1] >> 4;
        }
    }

    /* Restore pre-c2t bases from YS:Z so MethylDackel sees real C/Ts.
     * FASTQ ingest (fastmap.cpp meth_mode block) builds comments as
     * "YS:Z:<l_seq bytes>\tYC:Z:XX" starting at offset 0 — rely on that. */
    const char *orig_seq = NULL;
    if (s->comment && s->l_seq > 0
        && s->comment[0] == 'Y' && s->comment[1] == 'S' && s->comment[2] == ':') {
        orig_seq = s->comment + 5;
    }

    int emit_seq = !(p.flag & 0x100);
    size_t l_emit = 0;
    char *seq_text = NULL;
    char *qual_bin = NULL;
    if (emit_seq && qe > qb) {
        l_emit = (size_t)(qe - qb);
        seq_text = (char *)malloc(l_emit + 1);
        if (seq_text == NULL) { free(bam_cigar); return -1; }
        if (orig_seq != NULL) {
            /* Apply is_rev RC + supp soft-clip trim over pre-c2t bases. */
            if (!p.is_rev) {
                for (size_t i = 0; i < l_emit; ++i) {
                    char c = orig_seq[qb + (int)i];
                    seq_text[i] = (c >= 'a' && c <= 'z') ? (char)(c - 32) : c;
                }
            } else {
                for (size_t i = 0; i < l_emit; ++i) {
                    unsigned char c = (unsigned char)orig_seq[qe - 1 - (int)i];
                    seq_text[i] = "TGCAN"[nst_nt4_table[c]];
                }
            }
        } else if (!p.is_rev) {
            /* Fallback: YS:Z (orig_seq) was missing. s->seq holds ASCII
             * (post-c2t in meth mode); map ASCII → 0..4 via nst_nt4_table. */
            for (size_t i = 0; i < l_emit; ++i) {
                unsigned char c = (unsigned char)s->seq[qb + (int)i];
                seq_text[i] = "ACGTN"[nst_nt4_table[c]];
            }
        } else {
            for (size_t i = 0; i < l_emit; ++i) {
                unsigned char c = (unsigned char)s->seq[qe - 1 - (int)i];
                seq_text[i] = "TGCAN"[nst_nt4_table[c]];
            }
        }
        seq_text[l_emit] = '\0';
        if (s->qual) {
            qual_bin = (char *)malloc(l_emit);
            if (qual_bin == NULL) { free(seq_text); free(bam_cigar); return -1; }
            if (!p.is_rev) {
                for (size_t i = 0; i < l_emit; ++i) qual_bin[i] = (char)((unsigned char)s->qual[qb + (int)i] - 33);
            } else {
                for (size_t i = 0; i < l_emit; ++i) qual_bin[i] = (char)((unsigned char)s->qual[qe - 1 - (int)i] - 33);
            }
        }
    }

    /* Chimera QC + set-as-failed applied here so group propagation upstream
     * only scans flags. */
    uint8_t mapq = p.mapq;
    int mapped = !(flag16 & 0x4) && direction != 0;
    if (mapped) {
        if (opt->meth_set_as_failed != 0 && opt->meth_set_as_failed == direction) {
            flag16 |= 0x200;
        }
        if (!opt->meth_no_chim && p.n_cigar > 0 && s->l_seq > 0) {
            int lm = cigar_longest_m_mem(p.cigar, p.n_cigar);
            if (100 * lm < MIN_LONGEST_M_PCT * s->l_seq) {
                flag16 |= 0x200;
                flag16 &= ~0x2;
                if (mapq > 1) mapq = 1;
            }
        }
    }

    /* Build the bam1_t. bam_set1 handles 4-bit packing, name storage, etc. */
    int ret = bam_set1(b,
                       strlen(s->name), s->name,
                       flag16,
                       tid,
                       (hts_pos_t)p.pos,
                       mapq,
                       bam_n_cigar, bam_cigar,
                       mtid,
                       mp ? (hts_pos_t)mp->pos : -1,
                       tlen,
                       l_emit, seq_text, qual_bin,
                       /* l_aux */ 0);

    free(bam_cigar);
    free(seq_text);
    free(qual_bin);
    if (ret < 0) return -1;

    /* Aux tags — roughly match mem_aln2sam emission order */
    if (p.n_cigar > 0) {
        int32_t nm = (int32_t)p.NM;
        bam_aux_append(b, "NM", 'i', sizeof(nm), (const uint8_t *)&nm);
        const char *md = (const char *)(p.cigar + p.n_cigar);
        bam_aux_append(b, "MD", 'Z', (int)strlen(md) + 1, (const uint8_t *)md);
    }
    if (mp && mp->n_cigar > 0) {
        /* Growable so long mate CIGARs don't silently truncate. */
        kstring_t mc = {0, 0, NULL};
        for (int i = 0; i < mp->n_cigar; ++i) {
            int op = mp->cigar[i] & 0xf;
            int len = mp->cigar[i] >> 4;
            if (!(opt->flag & MEM_F_SOFTCLIP) && !mp->is_alt && (op == 3 || op == 4))
                op = which ? 4 : 3;
            ksprintf(&mc, "%d%c", len, "MIDSH"[op]);
        }
        if (mc.l > 0)
            bam_aux_append(b, "MC", 'Z', (int)mc.l + 1, (const uint8_t *)mc.s);
        free(mc.s);
    }
    if (p.score >= 0) {
        int32_t as = (int32_t)p.score;
        bam_aux_append(b, "AS", 'i', sizeof(as), (const uint8_t *)&as);
    }
    if (p.sub >= 0) {
        int32_t xs = (int32_t)p.sub;
        bam_aux_append(b, "XS", 'i', sizeof(xs), (const uint8_t *)&xs);
    }
    if (bwa_rg_id[0]) {
        bam_aux_append(b, "RG", 'Z', (int)strlen(bwa_rg_id) + 1, (const uint8_t *)bwa_rg_id);
    }
    if (p.alt_sc > 0) {
        float pa_f = (float)((double)p.score / (double)p.alt_sc);
        bam_aux_append(b, "pa", 'f', sizeof(pa_f), (const uint8_t *)&pa_f);
    }
    /* SA:Z (other primary hits) — mirrors mem_aln2sam, but chrom names are
     * rewritten through cmap so split-alignment linkage references the
     * consolidated output contig rather than the doubled c2t f/r contig. */
    if (!(p.flag & 0x100)) {
        int has_other = 0;
        for (int i = 0; i < n_alns; ++i)
            if (i != which && !(list[i].flag & 0x100)) { has_other = 1; break; }
        if (has_other) {
            kstring_t sa = {0, 0, NULL};
            for (int i = 0; i < n_alns; ++i) {
                const mem_aln_t *r = &list[i];
                if (i == which || (r->flag & 0x100)) continue;
                const char *r_name = NULL;
                if (r->rid >= 0 && r->rid < cmap->n_internal) {
                    int r_out = cmap->out_tid[r->rid];
                    if (r_out >= 0 && r_out < cmap->n_output)
                        r_name = cmap->output_names[r_out];
                }
                if (r_name == NULL) r_name = "*";
                kputs(r_name, &sa); kputc(',', &sa);
                kputl(r->pos + 1, &sa); kputc(',', &sa);
                kputc("+-"[r->is_rev], &sa); kputc(',', &sa);
                for (int k = 0; k < r->n_cigar; ++k) {
                    kputw(r->cigar[k] >> 4, &sa);
                    kputc("MIDSH"[r->cigar[k] & 0xf], &sa);
                }
                kputc(',', &sa); kputw(r->mapq, &sa);
                kputc(',', &sa); kputw(r->NM, &sa);
                kputc(';', &sa);
            }
            if (sa.l > 0)
                bam_aux_append(b, "SA", 'Z', (int)sa.l + 1, (const uint8_t *)sa.s);
            free(sa.s);
        }
    }
    /* XA:Z — like SA, rewrite each entry's contig name through cmap. Entries
     * look like `name,+/-pos,cigar,NM;`; the name up to the first comma is
     * the doubled-ref 'f'/'r'-prefixed contig. Strip the leading f/r so the
     * tag matches the consolidated @SQ names published by meth_bam_writer_open. */
    if (p.XA != NULL) {
        kstring_t xa = {0, 0, NULL};
        for (const char *entry = p.XA; *entry; ) {
            const char *entry_end = strchr(entry, ';');
            size_t entry_len = entry_end ? (size_t)(entry_end - entry) : strlen(entry);
            if (entry_len > 0) {
                const char *name_end = (const char *)memchr(entry, ',', entry_len);
                size_t name_len = name_end ? (size_t)(name_end - entry) : entry_len;
                const char *name = entry;
                if (name_len > 1 && (name[0] == 'f' || name[0] == 'r')) {
                    ++name;
                    --name_len;
                }
                kputsn(name, (int)name_len, &xa);
                if (name_end) {
                    kputsn(name_end, (int)(entry_len - (name_end - entry)), &xa);
                }
                kputc(';', &xa);
            }
            if (!entry_end) break;
            entry = entry_end + 1;
        }
        if (xa.l > 0)
            bam_aux_append(b, "XA", 'Z', (int)xa.l + 1, (const uint8_t *)xa.s);
        free(xa.s);
    }
    /* YD:Z — meth strand hypothesis */
    if (mapped) {
        char yd[2] = { direction, '\0' };
        bam_aux_append(b, "YD", 'Z', 2, (const uint8_t *)yd);
    }

    /* Generic aux shared with the --bam path so --meth -C emits FASTQ tags
     * plus the YS:Z/YC:Z meth tags that `fastmap.cpp` built into s->comment,
     * and --meth -V emits XR:Z. p.rid here is the bwa-mem3 internal contig
     * index (pre-remap), which is what bns uses. */
    bam_writer_append_generic_aux(b, s, opt, bns, p.rid);

    return 0;
}

void meth_bam_group_propagate_qcfail(bam1_t **group, int n)
{
    if (group == NULL || n <= 0) return;
    int any_fail = 0;
    for (int i = 0; i < n; ++i) {
        if (group[i] != NULL && (group[i]->core.flag & 0x200)) { any_fail = 1; break; }
    }
    if (!any_fail) return;
    for (int i = 0; i < n; ++i) {
        if (group[i] == NULL) continue;
        group[i]->core.flag |= 0x200;
        group[i]->core.flag &= (uint16_t)~0x2;
    }
}
