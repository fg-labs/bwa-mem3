/* SPDX-License-Identifier: MIT */
#ifndef BWA_MEM3_ALN_METH_REALIGN_H
#define BWA_MEM3_ALN_METH_REALIGN_H
#include "aln_meth.h"
namespace bwa_aln {
/* Original forward-reference half-open window. Query is already oriented to this strand. */
struct aln_meth_candidate_t { int rid, strand, hyp; int64_t rb, re; };
enum class aln_meth_eval_status { ok, invalid_input, resource_error, no_alignment };
/* Fetch a validated forward-reference window (length <=65535), masking .amb runs to code N
 * before CpG labelling. Used only by aln's methylation refinement, accounting and emission. */
void aln_meth_fetch_window(const aln_meth_ref_t &ref, const aln_meth_opt_t &mo,
    int64_t rb, int len, int rid, std::vector<uint8_t> &target, int64_t site = -1);
/* --known-snvs: the site alleles (aln_site_key) in [rb, rb + len), in (position, ALT) order, at
 * most KS_RESCUE_MAX_SITES (the first ones; later ones stay at REF); none without sites. With rid
 * set, the range is first clipped to that contig, as mate rescue clips its window to the anchor's
 * contig, so sites in a neighbour cannot use up the cap. With site set to one of them,
 * aln_meth_fetch_window puts that ALT alone in place, every other site at REF, before ambiguity
 * masking and CpG labelling. */
void aln_meth_window_sites(const aln_meth_ref_t &ref, int64_t rb, int len, std::vector<int64_t> *keys,
    int rid = -1);
/* Fixed-reward candidate refiner. Caller cost-ranks its raw CIGAR and owns both CIGARs.
 * band_limit >= 0 replaces the hint- and certificate-derived band with exactly
 * max(band_limit, |len - span|): the best alignment whose indels stay within it (0 with an
 * equal span: ungapped). --known-snvs uses it to re-realign an ALT-window placement whose
 * best-scoring alignment breaks the read's limits.
 * Invalid/resource failures print ERROR and must not be interpreted as an absent candidate.
 * Fixed scoring uses the existing meth matrices and ksw_global; hints never certify a band.
 * Active-backend allocation estimate is capped at 64 MiB; existing retained TLS scratch is separate.
 * Query/window lengths are limited to 65535. No silent coordinate clipping or legacy bitfield casts. */
aln_meth_eval_status aln_meth_realign_window(const aln_meth_ref_t &ref, const aln_meth_opt_t &mo,
    const aln_meth_candidate_t &candidate, const ubyte_t *query, int len, int n_ins, int n_del,
    aln_meth_alignment_t *out, int64_t site = -1, int band_limit = -1);
}
#endif
