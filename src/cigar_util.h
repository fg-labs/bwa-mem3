/* SPDX-License-Identifier: MIT */
#ifndef BWAMEM3_CIGAR_UTIL_H
#define BWAMEM3_CIGAR_UTIL_H

#include <stdint.h>

/* bwa-mem3 CIGAR op encoding is 0=M, 1=I, 2=D, 3=S, 4=H.
 * BAM's is 0=M, 1=I, 2=D, 3=N, 4=S, 5=H. Remap on emit. */
static const uint32_t BAM_OP_FROM_MEM[5] = { 0, 1, 2, 4, 5 };

/* Reference length on a bwa-mem3-encoded CIGAR (op 0=M consumes ref + matches;
 * op 2=D consumes ref only). bwa-mem3 never emits =/X so only op 0 and 2 count. */
static inline int64_t cigar_ref_len_mem(const uint32_t *cigar, int n)
{
    int64_t l = 0;
    for (int i = 0; i < n; ++i) {
        int op = cigar[i] & 0xf;
        if (op == 0 || op == 2) l += cigar[i] >> 4;
    }
    return l;
}

/* Longest M-run on a bwa-mem3-encoded CIGAR. */
static inline int cigar_longest_m_mem(const uint32_t *cigar, int n)
{
    int longest = 0;
    for (int i = 0; i < n; ++i) {
        if ((cigar[i] & 0xf) == 0 && (int)(cigar[i] >> 4) > longest)
            longest = (int)(cigar[i] >> 4);
    }
    return longest;
}

#endif
