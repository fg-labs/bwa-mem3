/* Shared by the lockstep parity tests (smem_lockstep_parity_test,
 * bwtseed_lockstep_parity_test): read encoding, a small FASTA loader and a
 * reverse complement. */
#ifndef LOCKSTEP_PARITY_TEST_UTIL_H
#define LOCKSTEP_PARITY_TEST_UTIL_H

#include <ctype.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <string>
#include "bwa.h"

/* Encode a vector of read strings (ACGT->0123, N->4) into the flat enc_qdb
 * + bseq1_t[] + query_cum_len_ar[] layout the seeding drivers expect. The
 * caller frees the three outputs. */
static inline void encode_reads(const char * const *reads, int32_t numReads,
                                uint8_t **out_enc_qdb,
                                bseq1_t **out_seq,
                                int32_t **out_cum_len)
{
    int32_t total_len = 0;
    for (int i = 0; i < numReads; i++) total_len += (int32_t)strlen(reads[i]);
    uint8_t *enc = (uint8_t *)calloc(total_len > 0 ? total_len : 1, 1);
    bseq1_t *seq = (bseq1_t *)calloc(numReads, sizeof(bseq1_t));
    int32_t *cum = (int32_t *)calloc(numReads + 1, sizeof(int32_t));
    int32_t off = 0;
    for (int i = 0; i < numReads; i++) {
        int32_t l = (int32_t)strlen(reads[i]);
        cum[i] = off;
        for (int j = 0; j < l; j++) {
            char c = reads[i][j];
            uint8_t e = 4;
            if (c == 'A' || c == 'a') e = 0;
            else if (c == 'C' || c == 'c') e = 1;
            else if (c == 'G' || c == 'g') e = 2;
            else if (c == 'T' || c == 't') e = 3;
            enc[off + j] = e;
        }
        seq[i].l_seq = l;
        off += l;
    }
    cum[numReads] = off;
    *out_enc_qdb = enc;
    *out_seq = seq;
    *out_cum_len = cum;
}

/* Read a (small) FASTA into a flat uppercase sequence string. Used to derive
 * reads that actually MATCH the reference, so the walks run against real SA
 * intervals rather than non-matching synthetic reads. Returns a malloc'd
 * buffer (caller frees) or NULL if the file cannot be opened / has no
 * sequence. */
static inline char *read_fasta_seq(const char *path, int32_t *out_len) {
    FILE *f = fopen(path, "r");
    if (!f) return NULL;
    size_t cap = 8192, len = 0;
    char *seq = (char *)malloc(cap);
    if (!seq) { fclose(f); return NULL; }
    int c;
    bool in_header = false;
    while ((c = fgetc(f)) != EOF) {
        if (c == '>') { in_header = true; continue; }
        if (c == '\n') { in_header = false; continue; }
        if (in_header) continue;
        if (c == '\r' || c == ' ' || c == '\t') continue;
        if (len + 1 >= cap) {
            cap *= 2;
            char *tmp = (char *)realloc(seq, cap);
            if (!tmp) { free(seq); fclose(f); return NULL; }
            seq = tmp;
        }
        seq[len++] = (char)toupper(c);
    }
    fclose(f);
    if (len == 0) { free(seq); return NULL; }
    seq[len] = '\0';
    *out_len = (int32_t)len;
    return seq;
}

/* Reverse complement of an uppercase sequence; non-ACGT bases (N) are kept. */
static inline std::string revcomp(const std::string &s) {
    std::string r(s.rbegin(), s.rend());
    for (char &c : r)
        c = c == 'A' ? 'T' : c == 'C' ? 'G' : c == 'G' ? 'C' : c == 'T' ? 'A' : c;
    return r;
}

#endif
