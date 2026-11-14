/* SPDX-License-Identifier: MIT */
/*
 * The --meth original reference: the un-converted bns + pac that every step
 * after seeding (mem) or the search (aln) uses, loaded from the original index
 * prefix that `bwa-mem3 index --meth` writes next to the converted one.
 */
#ifndef BWA_MEM3_METH_REF_H
#define BWA_MEM3_METH_REF_H

#include <stdint.h>

#include "bntseq.h"

/* Loads <prefix>.{amb,ann,pac} (bns_restore, then the whole .pac read with
 * `pread_workers` workers; 1 reads serially). Returns 0, or -1 (with an error
 * printed and *bns_out/*pac_out NULL) when bns_restore or the pac allocation
 * fails; an I/O error while reading the .pac is fatal, as for the seed index. */
int meth_orig_ref_load(const char *prefix, bntseq_t **bns_out, uint8_t **pac_out, int pread_workers);

/* Frees what meth_orig_ref_load returned and NULLs the handles (NULL-safe). */
void meth_orig_ref_free(bntseq_t **bns, uint8_t **pac);

#endif
