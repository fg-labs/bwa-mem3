/* SPDX-License-Identifier: MIT */
#ifndef BWA_MEM3_ALN_METH_REALIGN_H
#define BWA_MEM3_ALN_METH_REALIGN_H
#include "aln_meth.h"
namespace bwa_aln {
/* Original forward-reference half-open window. Query is already oriented to this strand. */
struct aln_meth_candidate_t { int rid, strand, hyp; int64_t rb, re; };
enum class aln_meth_eval_status { ok, invalid_input, resource_error, no_alignment };
/* Fixed-reward candidate refiner. Caller cost-ranks its raw CIGAR and owns both CIGARs.
 * Invalid/resource failures print ERROR and must not be interpreted as an absent candidate.
 * Fixed scoring uses the existing meth matrices and ksw_global; hints never certify a band.
 * Active-backend allocation estimate is capped at 64 MiB; existing retained TLS scratch is separate.
 * Query/window lengths are limited to 65535. No silent coordinate clipping or legacy bitfield casts. */
aln_meth_eval_status aln_meth_realign_window(const aln_meth_ref_t &ref, const aln_meth_opt_t &mo,
    const aln_meth_candidate_t &candidate, const ubyte_t *query, int len, int n_ins, int n_del,
    aln_meth_alignment_t *out);
}
#endif
