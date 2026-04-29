/* SPDX-License-Identifier: MIT */

/* htslib headers before any bwa-mem3 header that pulls in kstring.h
 * (they share the KSTRING_H include guard). */
#include "htslib/sam.h"
#include "htslib/kstring.h"

#include "bam_writer.h"
#include "cigar_util.h"
#include "utils.h"

#include <cstdio>
#include <cstdlib>
#include <cstring>

struct bam_writer_s {
    htsFile   *fp;
    sam_hdr_t *hdr;
};

extern char bwa_rg_id[256];

extern "C" {
struct bam1_t *bam_writer_alloc(void)             { return bam_init1(); }
void           bam_writer_free(struct bam1_t *b)  { if (b) bam_destroy1(b); }
}

bam_writer_t *bam_writer_open(const char *path, const bntseq_t *bns,
                              const char *idx_hdr_lines,
                              const char *hdr_line, const char *bwa_pg,
                              int compression_level)
{
    if (path == NULL || bns == NULL) return NULL;

    // Detect whether the index .hdr/.dict content already supplies @SQ or
    // @HD records. If @SQ: skip auto-generating @SQ from `bns` — adding
    // both would hit htslib's "duplicate SN" de-dup and fail. If @HD:
    // skip the default @HD — htslib's sam_hdr_add_lines does not de-dup
    // @HD records, so emitting both would produce two @HD lines. The
    // caller (fastmap.cpp) has already nulled `idx_hdr_lines` when the
    // user's -H supplies @SQ, so if we see @SQ here, they're authoritative.
    int idx_has_sq = 0, idx_has_hd = 0;
    if (idx_hdr_lines != NULL && idx_hdr_lines[0] != '\0') {
        if (strncmp(idx_hdr_lines, "@SQ\t", 4) == 0 ||
            strstr(idx_hdr_lines, "\n@SQ\t") != NULL)
            idx_has_sq = 1;
        if (strncmp(idx_hdr_lines, "@HD\t", 4) == 0 ||
            strstr(idx_hdr_lines, "\n@HD\t") != NULL)
            idx_has_hd = 1;
    }
    int user_has_hd = 0;
    if (hdr_line != NULL && hdr_line[0] != '\0') {
        if (strncmp(hdr_line, "@HD\t", 4) == 0 ||
            strstr(hdr_line, "\n@HD\t") != NULL)
            user_has_hd = 1;
    }

    sam_hdr_t *hdr = sam_hdr_init();
    if (hdr == NULL) return NULL;
    // Emit a default @HD only when neither the user's -H nor the index's
    // .hdr/.dict supplies one. Precedence (user > index > default) is
    // enforced by ordering: we add idx_hdr_lines before hdr_line, and the
    // SAM text path uses the same precedence via bwa_print_sam_hdr2.
    if (!idx_has_hd && !user_has_hd &&
        sam_hdr_add_line(hdr, "HD", "VN", "1.6", "SO", "unsorted", NULL) < 0) goto fail;
    if (!idx_has_sq) {
        for (int i = 0; i < bns->n_seqs; ++i) {
            char len_buf[32];
            snprintf(len_buf, sizeof(len_buf), "%lld", (long long)bns->anns[i].len);
            if (sam_hdr_add_line(hdr, "SQ", "SN", bns->anns[i].name, "LN", len_buf, NULL) < 0)
                goto fail;
        }
    }
    // Merge the index's .hdr/.dict records after the default @HD (and
    // auto-generated @SQ, when the index didn't supply its own) but before
    // the user's -H lines. If the user has an @HD in -H, strip the @HD
    // records from idx_hdr_lines first — htslib's sam_hdr_add_lines does
    // not de-dup @HD, so without filtering we'd emit two @HD records. This
    // matches the SAM text path's precedence: user > index > default.
    if (idx_hdr_lines != NULL && idx_hdr_lines[0] != '\0') {
        const char *to_add = idx_hdr_lines;
        char *filtered = NULL;
        if (user_has_hd && idx_has_hd) {
            // Copy idx_hdr_lines, dropping any @HD records.
            size_t n = strlen(idx_hdr_lines);
            filtered = (char *)malloc(n + 1);
            if (filtered == NULL) goto fail;
            size_t w = 0;
            const char *p = idx_hdr_lines;
            while (*p) {
                const char *eol = strchr(p, '\n');
                size_t len = eol ? (size_t)(eol - p) : strlen(p);
                int is_hd = (len >= 4 && strncmp(p, "@HD\t", 4) == 0);
                if (!is_hd && len > 0) {
                    memcpy(filtered + w, p, len); w += len;
                    filtered[w++] = '\n';
                }
                p = eol ? eol + 1 : p + len;
            }
            filtered[w] = '\0';
            to_add = filtered;
        }
        int rc = (to_add[0] != '\0') ? sam_hdr_add_lines(hdr, to_add, 0) : 0;
        free(filtered);
        if (rc < 0) goto fail;
    }
    if (hdr_line != NULL && hdr_line[0] != '\0' && sam_hdr_add_lines(hdr, hdr_line, 0) < 0)
        goto fail;
    if (bwa_pg != NULL && bwa_pg[0] != '\0' && sam_hdr_add_lines(hdr, bwa_pg, 0) < 0)
        goto fail;

    {
        if (compression_level < 0) compression_level = 0;
        if (compression_level > 9) compression_level = 9;
        char mode[8];
        snprintf(mode, sizeof(mode), "wb%d", compression_level);
        htsFile *fp = hts_open(path, mode);
        if (fp == NULL) goto fail;
        if (sam_hdr_write(fp, hdr) < 0) { hts_close(fp); goto fail; }
        bam_writer_t *w = (bam_writer_t *)calloc(1, sizeof(*w));
        if (w == NULL) { hts_close(fp); goto fail; }
        w->fp = fp; w->hdr = hdr;
        return w;
    }
fail:
    sam_hdr_destroy(hdr);
    return NULL;
}

int bam_writer_write(bam_writer_t *w, struct bam1_t *b)
{
    if (w == NULL || b == NULL) return -1;
    return sam_write1(w->fp, w->hdr, b);
}

int bam_writer_close(bam_writer_t *w)
{
    if (w == NULL) return 0;
    int rc = 0;
    if (w->hdr) { sam_hdr_destroy(w->hdr); w->hdr = NULL; }
    if (w->fp)  { if (hts_close(w->fp) < 0) rc = -1; w->fp = NULL; }
    free(w);
    return rc;
}

/* ------------------------------------------------------------------- */
/* mem_aln_t -> bam1_t                                                  */
/* ------------------------------------------------------------------- */

/* Append SAM-text aux fields ("TAG:TYPE:VALUE" tokens, TAB-separated) onto
 * `b` as packed BAM aux. Mirrors the literal s->comment append that
 * mem_aln2sam does on the SAM path so --bam preserves -C output. Only TAB
 * is treated as a separator (matching SAM's grammar); spaces are kept
 * inside Z values, which is exactly how the SAM path treats them when the
 * FASTQ comment is appended verbatim. Only the common scalar types
 * (A/i/f/Z/H) are handled; B (typed array) is skipped because parsing it
 * is non-trivial and FASTQ comments effectively never use it. Malformed
 * tokens are skipped silently. */
static void append_sam_aux_tokens(struct bam1_t *b, const char *s)
{
    if (s == NULL) return;
    while (*s) {
        while (*s == '\t') ++s;
        if (!*s) break;
        const char *tok = s;
        while (*s && *s != '\t') ++s;
        size_t tok_len = (size_t)(s - tok);
        /* Need at least "XX:T:" (5 chars) and a non-empty value. */
        if (tok_len < 6 || tok[2] != ':' || tok[4] != ':') continue;
        char tag[2] = { tok[0], tok[1] };
        char type   = tok[3];
        const char *val  = tok + 5;
        size_t      vlen = tok_len - 5;
        switch (type) {
            case 'A':
                if (vlen >= 1) {
                    char c = val[0];
                    bam_aux_append(b, tag, 'A', 1, (const uint8_t *)&c);
                }
                break;
            case 'i': {
                char *end = NULL;
                int64_t iv = (int64_t)strtoll(val, &end, 10);
                if (end != val) {
                    int32_t i32 = (int32_t)iv;
                    bam_aux_append(b, tag, 'i', sizeof(i32),
                                   (const uint8_t *)&i32);
                }
                break;
            }
            case 'f': {
                char *end = NULL;
                float fv = strtof(val, &end);
                if (end != val) {
                    bam_aux_append(b, tag, 'f', sizeof(fv),
                                   (const uint8_t *)&fv);
                }
                break;
            }
            case 'Z':
            case 'H': {
                char *buf = (char *)malloc(vlen + 1);
                if (buf == NULL) break;
                memcpy(buf, val, vlen);
                buf[vlen] = '\0';
                bam_aux_append(b, tag, type, (int)vlen + 1,
                               (const uint8_t *)buf);
                free(buf);
                break;
            }
            default:
                /* B-type and unknown — leave unconverted rather than emit a
                 * wrong-typed tag. */
                break;
        }
    }
}

int mem_aln_to_bam(struct bam1_t *b,
                   const mem_opt_t *opt, const bntseq_t *bns,
                   const bseq1_t *s, int n_alns,
                   const mem_aln_t *list, int which,
                   const mem_aln_t *m_)
{
    if (b == NULL || opt == NULL || s == NULL || list == NULL) return -1;

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

    uint16_t flag16 = (uint16_t)((p.flag & 0xffff) | (p.flag & 0x10000 ? 0x100 : 0));

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
                op = which ? 4 : 3;
            uint32_t bam_op = (op >= 0 && op < 5) ? BAM_OP_FROM_MEM[op] : 0;
            bam_cigar[i] = ((uint32_t)len << 4) | bam_op;
        }
        bam_n_cigar = (size_t)p.n_cigar;
    }

    hts_pos_t tlen = 0;
    if (mp && mp->rid >= 0 && p.rid == mp->rid && p.n_cigar > 0 && m.n_cigar > 0) {
        /* TLEN uses ref-consumed lengths. bwa-mem3 and BAM both encode
         * M=0, D=2, so we can count directly on the pre-remap CIGARs. */
        int64_t p_rlen = cigar_ref_len_mem(p.cigar, p.n_cigar);
        int64_t m_rlen = cigar_ref_len_mem(m.cigar, m.n_cigar);
        int64_t p0 = p.pos + (p.is_rev ? p_rlen - 1 : 0);
        int64_t p1 = m.pos + (m.is_rev ? m_rlen - 1 : 0);
        tlen = -(p0 - p1 + (p0 > p1 ? 1 : p0 < p1 ? -1 : 0));
    }

    /* SEQ/QUAL range with supp soft-clip trim */
    int qb = 0, qe = s->l_seq;
    if (p.n_cigar && which && !(opt->flag & MEM_F_SOFTCLIP) && !p.is_alt) {
        int c0 = p.cigar[0] & 0xf;
        int cN = p.cigar[p.n_cigar-1] & 0xf;
        if (!p.is_rev) {
            if (c0 == 3 || c0 == 4) qb += p.cigar[0] >> 4;
            if (cN == 3 || cN == 4) qe -= p.cigar[p.n_cigar-1] >> 4;
        } else {
            if (c0 == 3 || c0 == 4) qe -= p.cigar[0] >> 4;
            if (cN == 3 || cN == 4) qb += p.cigar[p.n_cigar-1] >> 4;
        }
    }

    int emit_seq = !(p.flag & 0x100);
    size_t l_emit = 0;
    char *seq_text = NULL;
    char *qual_bin = NULL;
    if (emit_seq && qe > qb) {
        l_emit = (size_t)(qe - qb);
        seq_text = (char *)malloc(l_emit + 1);
        if (seq_text == NULL) { free(bam_cigar); return -1; }
        if (!p.is_rev) {
            for (size_t i = 0; i < l_emit; ++i) seq_text[i] = "ACGTN"[(int)s->seq[qb + (int)i]];
        } else {
            for (size_t i = 0; i < l_emit; ++i) seq_text[i] = "TGCAN"[(int)s->seq[qe - 1 - (int)i]];
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

    int ret = bam_set1(b,
                       strlen(s->name), s->name,
                       flag16, p.rid, (hts_pos_t)p.pos, p.mapq,
                       bam_n_cigar, bam_cigar,
                       mp ? mp->rid : -1, mp ? (hts_pos_t)mp->pos : -1, tlen,
                       l_emit, seq_text, qual_bin,
                       /* l_aux */ 0);
    free(bam_cigar); free(seq_text); free(qual_bin);
    if (ret < 0) return -1;

    if (p.n_cigar > 0) {
        int32_t nm = (int32_t)p.NM;
        bam_aux_append(b, "NM", 'i', sizeof(nm), (const uint8_t *)&nm);
        const char *md = (const char *)(p.cigar + p.n_cigar);
        bam_aux_append(b, "MD", 'Z', (int)strlen(md) + 1, (const uint8_t *)md);
    }
    if (mp && mp->n_cigar > 0) {
        /* Dynamic buffer: CIGARs with >~800 ops would silently truncate a
         * fixed stack buffer. kstring_t grows as needed. */
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
    if (mp) {
        int32_t mq = (int32_t)mp->mapq;
        bam_aux_append(b, "MQ", 'i', sizeof(mq), (const uint8_t *)&mq);
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
    /* SA:Z (other primary hits) — mirrors mem_aln2sam */
    if (!(p.flag & 0x100)) {
        int has_other = 0;
        for (int i = 0; i < n_alns; ++i)
            if (i != which && !(list[i].flag & 0x100)) { has_other = 1; break; }
        if (has_other) {
            kstring_t sa = {0, 0, NULL};
            for (int i = 0; i < n_alns; ++i) {
                const mem_aln_t *r = &list[i];
                if (i == which || (r->flag & 0x100)) continue;
                kputs(bns->anns[r->rid].name, &sa); kputc(',', &sa);
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
    if (p.XA != NULL) {
        bam_aux_append(b, "XA", 'Z', (int)strlen(p.XA) + 1, (const uint8_t *)p.XA);
    }
    if (p.HN >= 0) {
        int32_t hn = (int32_t)p.HN;
        bam_aux_append(b, "HN", 'i', sizeof(hn), (const uint8_t *)&hn);
    }

    bam_writer_append_generic_aux(b, s, opt, bns, p.rid);

    return 0;
}

extern "C" void bam_writer_append_generic_aux(struct bam1_t *b,
                                              const bseq1_t *s,
                                              const mem_opt_t *opt,
                                              const bntseq_t *bns,
                                              int rid)
{
    if (b == NULL || s == NULL || opt == NULL) return;

    /* FASTQ-carried tags (-C; also YS:Z/YC:Z in --meth mode). mem_aln2sam
     * appends s->comment literally to the SAM line; here we parse it as
     * SAM-text aux tokens and emit each as a packed BAM aux so the BAM and
     * SAM paths carry the same tags. */
    if (s->comment != NULL && s->comment[0] != '\0') {
        append_sam_aux_tokens(b, s->comment);
    }

    /* Reference annotation (-V / MEM_F_REF_HDR). Mirrors mem_aln2sam: emit
     * bns->anns[rid].anno as XR:Z, replacing TAB with SPACE. */
    if ((opt->flag & MEM_F_REF_HDR) && rid >= 0 &&
        bns != NULL && bns->anns[rid].anno != NULL &&
        bns->anns[rid].anno[0] != '\0') {
        const char *anno = bns->anns[rid].anno;
        size_t alen = strlen(anno);
        char *xr = (char *)malloc(alen + 1);
        if (xr != NULL) {
            for (size_t i = 0; i < alen; ++i)
                xr[i] = (anno[i] == '\t') ? ' ' : anno[i];
            xr[alen] = '\0';
            bam_aux_append(b, "XR", 'Z', (int)alen + 1, (const uint8_t *)xr);
            free(xr);
        }
    }
}

extern "C" void bam_writer_bseq_push(bseq1_t *s, struct bam1_t *b)
{
    if (s == NULL || b == NULL) return;
    if (s->n_bams == s->cap_bams) {
        int new_cap = s->cap_bams ? s->cap_bams * 2 : 4;
        void **nb = (void **)realloc(s->bams, (size_t)new_cap * sizeof(void *));
        if (nb == NULL) {
            bam_destroy1(b);
            err_fatal(__func__, "out of memory growing per-read BAM list");
        }
        s->bams = nb;
        s->cap_bams = new_cap;
    }
    s->bams[s->n_bams++] = (void *)b;
}

int bam_writer_push_aln(bseq1_t *s,
                        const mem_opt_t *opt, const bntseq_t *bns,
                        int n_alns, const mem_aln_t *list, int which,
                        const mem_aln_t *m)
{
    if (s == NULL) return -1;
    bam1_t *b = bam_init1();
    if (b == NULL) return -1;
    if (mem_aln_to_bam(b, opt, bns, s, n_alns, list, which, m) < 0) {
        bam_destroy1(b);
        return -1;
    }
    bam_writer_bseq_push(s, b);
    return 0;
}
