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
    int64_t rb, int len, int rid, std::vector<uint8_t> &target, int alt = 0);
/* Whether [rb, rb + len) of contig rid holds a --known-snvs site (false without sites). With
 * alt set, aln_meth_fetch_window puts each site's ALT in place (the first ALT at a multi-allelic
 * position) before ambiguity masking and CpG labelling. */
bool aln_meth_window_has_site(const aln_meth_ref_t &ref, int64_t rb, int len, int rid);
/* Fixed-reward candidate refiner. Caller cost-ranks its raw CIGAR and owns both CIGARs.
 * Invalid/resource failures print ERROR and must not be interpreted as an absent candidate.
 * Fixed scoring uses the existing meth matrices and ksw_global; hints never certify a band.
 * Active-backend allocation estimate is capped at 64 MiB; existing retained TLS scratch is separate.
 * Query/window lengths are limited to 65535. No silent coordinate clipping or legacy bitfield casts.
 * An equal-length window whose ungapped alignment provably beats every gapped one skips the DP; the
 * result (and band) is what ksw_global would have returned. window, when set, gets the fetched
 * window (aln_meth_fetch_window's bases for the candidate and alt), for a caller that reads it
 * again. */
aln_meth_eval_status aln_meth_realign_window(const aln_meth_ref_t &ref, const aln_meth_opt_t &mo,
    const aln_meth_candidate_t &candidate, const ubyte_t *query, int len, int n_ins, int n_del,
    aln_meth_alignment_t *out, int alt = 0, std::vector<uint8_t> *window = nullptr);
}
#endif
