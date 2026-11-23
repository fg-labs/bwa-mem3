/* SPDX-License-Identifier: MIT */
#ifndef BWA_MEM3_ALN_METH_REALIGN_H
#define BWA_MEM3_ALN_METH_REALIGN_H
#include "aln_meth.h"
namespace bwa_aln {
/* Original forward-reference half-open window. Query is already oriented to this strand. */
/* cig, when set: the raw CIGAR of a search path to the window (a --known-snvs window hit's), which
 * the main search scores as one of its own paths. */
// blind: a row the ALT-window search's -n (or -k) limit dropped only for its mirror mismatches
// (ks_meth_candidates), which only the --meth-mapq conversion-blind strata count. Brace-initialize
// a candidate: an omitted blind is false.
struct aln_meth_candidate_t { int rid, strand, hyp; int64_t rb, re; std::vector<uint32_t> cig; bool blind; };
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
 * Query/window lengths are limited to 65535. No silent coordinate clipping or legacy bitfield casts.
 * An equal-length window whose ungapped alignment provably beats every gapped one skips the DP; the
 * result (and band) is what ksw_global would have returned. window, when set, gets the fetched
 * window (aln_meth_fetch_window's bases for the candidate and alt), for a caller that reads it
 * again. */
aln_meth_eval_status aln_meth_realign_window(const aln_meth_ref_t &ref, const aln_meth_opt_t &mo,
    const aln_meth_candidate_t &candidate, const ubyte_t *query, int len, int n_ins, int n_del,
    aln_meth_alignment_t *out, int64_t site = -1, std::vector<uint8_t> *window = nullptr,
    int band_limit = -1);

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
    std::vector<uint64_t> seed; // the caller's mask of positions whose differences count toward -k
};
void aln_meth_read_bits(const ubyte_t *query, int len, aln_meth_read_bits_t *out);
struct aln_meth_ungapped_t { int reward, n_mm, n_conv; };
/* Scores rd against [rb, rb + rd.len) of the original reference under hypothesis hyp (1: OT,
 * C->T; 0: OB, G->A). Reference .amb runs count as N. Requires f.ok. */
void aln_meth_score_ungapped(const aln_meth_ref_t &ref, const aln_meth_fast_t &f, int hyp, int64_t rb,
                             const aln_meth_read_bits_t &rd, aln_meth_ungapped_t *out);
/* The path scorer's ungapped block: query [qs, qs + n) against the original reference [rpos, rpos +
 * n), which starts at offset wo of a window whose ambiguity runs are amb (aln_meth_amb_bits; one bit
 * per base in the low bit of a 2-bit lane, as the read's planes). Adds the block's strict
 * mismatches (an N on either side included), the mirror-cell mismatches among them (OT read C over
 * reference T, OB read G over reference A; --meth-mapq conversion-blind), conversions, ambiguous
 * columns and, when rd.seed is set, strict and mirror-cell mismatches at the seed mask's positions.
 * Requires f.ok. */
struct aln_meth_block_counts_t { int n_mm = 0, n_conv = 0, n_amb = 0, n_seed_mm = 0, n_mirror = 0, n_seed_mirror = 0; };
void aln_meth_amb_bits(const aln_meth_ref_t &ref, int64_t rb, int len, std::vector<uint64_t> *bits);
void aln_meth_score_block(const aln_meth_ref_t &ref, const aln_meth_fast_t &f, int hyp, int64_t rpos, int qs, int n,
                          const aln_meth_read_bits_t &rd, const uint64_t *amb, int wo, aln_meth_block_counts_t *out);
}
#endif
