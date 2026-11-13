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

/* The equal-span ungapped scorer: the reward, mismatches and conversions of the [len M] alignment
 * of a query against the reference window [rb, rb + len), as aln_meth_realign_window (ungapped
 * reward) and the search's accounting compute them base by base, but 32 bases per 64-bit word.
 * Only for the 5-letter model (--meth-context all) and the reference haplotype; other settings
 * keep the base-by-base path (aln_meth_fast_t::ok is false). aln_meth_fast_init fills the
 * matrix classes (aln_meth.h) from mo's matrices; aln_meth_opt_finish calls it. */
void aln_meth_fast_init(const aln_meth_opt_t &mo, aln_meth_fast_t *f);
/* The query as 2-bit codes, 32 bases per word, first base in the top bits, and one-bit-per-base
 * planes (the low bit of each 2-bit lane) of its N, A, C, G and T positions. */
struct aln_meth_read_bits_t {
    int len = 0;
    std::vector<uint64_t> code, n, a, c, g, t;
};
void aln_meth_read_bits(const ubyte_t *query, int len, aln_meth_read_bits_t *out);
struct aln_meth_ungapped_t { int reward, n_mm, n_conv; };
/* Scores rd against [rb, rb + rd.len) of the original reference under hypothesis hyp (1: OT,
 * C->T; 0: OB, G->A). Reference .amb runs count as N. Requires f.ok. */
void aln_meth_score_ungapped(const aln_meth_ref_t &ref, const aln_meth_fast_t &f, int hyp, int64_t rb,
                             const aln_meth_read_bits_t &rd, aln_meth_ungapped_t *out);
}
#endif
