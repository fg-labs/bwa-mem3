/* SPDX-License-Identifier: MIT */
/* `aln --known-snvs`: the known SNV sites a run places reads against. New code. */
#ifndef BWA_MEM3_ALN_KNOWN_SNVS_H
#define BWA_MEM3_ALN_KNOWN_SNVS_H

#include "../bntseq.h"
#include "bwtaln.h"

#include <cstdint>
#include <string>
#include <vector>

class FMI_search;

namespace bwa_aln {

struct bwt_t;

/* One allele at one position: 0-based position on its contig, REF and ALT as 2-bit codes (A=0,
 * C=1, G=2, T=3). A multi-allelic record gives one entry per ALT at the same position; they are
 * alternatives, never combined into one haplotype. */
struct ks_site_t {
	int64_t pos;
	uint8_t ref, alt;
};

struct ks_sites_t {
	std::vector<std::vector<ks_site_t>> by_contig; // bns->n_seqs lists, each sorted by (pos, alt)
	int64_t n_sites = 0;       // alleles kept
	int64_t n_records = 0;     // VCF records read
	int64_t n_not_snv = 0;     // alleles skipped: not a single-base substitution (indel, MNV, symbolic)
	int64_t n_duplicate = 0;   // alleles skipped: the same position and ALT as an earlier record

	bool empty() const { return n_sites == 0; }
};

/* Loads the biallelic and multi-allelic SNVs of a VCF or BCF (optionally compressed) against the
 * reference `bns`/`pac` (the original reference under --meth). FILTER is ignored. A record whose
 * CHROM is not a reference contig, or whose REF base differs from the reference (or falls in an
 * ambiguous N run), is an error naming it. Returns 0, or -1 with an error on stderr. */
int ks_load(const char *path, const bntseq_t *bns, const uint8_t *pac, ks_sites_t *out);

/* Half-width of an ALT window: a read up to KS_WINDOW_HALF - (max gap extensions) bases long that
 * covers a site lies inside its window. Longer reads skip the window search. */
#define KS_WINDOW_HALF 160

/* One ALT window: the reference [beg, beg + len) of contig rid with the ALT of one site in place. */
struct ks_window_t {
	int rid;
	int64_t beg;
	int32_t len;
	int64_t site_pos; // the site, 0-based on rid
	uint8_t alt;
};

/* The ALT windows of a site set and their FM-index, built into a private temporary directory that
 * is removed when the windows are freed. Window i is contig i of the index. */
struct ks_windows_t {
	std::vector<ks_window_t> win;
	std::string dir;                     // the temporary directory (empty once removed)
	::FMI_search *fmi = nullptr;
	bwt_t *bwt = nullptr;
	const bntseq_t *bns = nullptr;       // the window index's contigs (one per window)
};

/* Builds one window per site allele from `pac` (its stored bases, so N runs carry the main index's
 * fill) and indexes them. Returns 0, or -1 with an error on stderr. */
int ks_windows_build(const ks_sites_t *sites, const bntseq_t *bns, const uint8_t *pac, int n_threads,
					 ks_windows_t *out);
void ks_windows_free(ks_windows_t *w);

/* What ks_override did to a read. */
enum ks_outcome_t {
	KS_UNTOUCHED = 0,   // no window hit covers its site: the stock result stands
	KS_SAME,            // window hits cover a site but change no hit's cost: stock result
	KS_CREDITED,        // placement, counts and MAPQ recomputed with the ALT credited
	KS_INCOMPLETE,      // too many hits to expand, or a read too long for its window: stock result
	KS_N_OUTCOMES
};

/* Crediting for one single-end read after bwa_cal_pac_pos (stock position and MAPQ set) and
 * before bwa_refine_gapped. `s->aln` holds the main search's hits; win_aln/win_n the same read's
 * hits on the window index. Window hits that cover their window's site join the main hits; per
 * locus the lower search score wins; then c1/c2, the primary (the stock one when still among the
 * best, else a per-read draw, never drand48), MAPQ (bwa_approx_mapQ) and the alternative hits are
 * recomputed with aln's cuts (score <= best + mismatch penalty, differences <= min(max_diff, best's
 * + 1)), applied before the main intervals are expanded. *alt_site is the packed position of the
 * site whose ALT the new primary's cost assumes, or -1. */
ks_outcome_t ks_override(const ks_windows_t *w, const bntseq_t *bns, const bwt_t *bwt, const gap_opt_t *opt,
						 int max_xa_hits, bwa_seq_t *s, const bwt_aln1_t *win_aln, int win_n, int64_t *alt_site);

/* Start tolerance for "the twin lands on the same locus" (refinement can shift a gapped start). */
#define KS_TWIN_SLOP 5

/* Where a twin read lands, without a draw: its main hits (t->aln) and window hits merged as for
 * ks_override. Returns 1 when it has one best hit at (pos, strand) within KS_TWIN_SLOP, 0 when not
 * (no hit, a tie, or elsewhere), -1 when incomplete; *mapq is its MAPQ. */
int ks_twin_place(const ks_windows_t *w, const bntseq_t *bns, const bwt_t *bwt, const gap_opt_t *opt,
				  const bwa_seq_t *t, const bwt_aln1_t *win_aln, int win_n, int64_t pos, int strand, int *mapq);

/* One site under a read's final placement. */
struct ks_vote_t {
	const ks_site_t *site;
	int rid;
	int qi;        // index in the reference-oriented read
	char allele;   // 'R', 'A' or 'O' (another base)
	int away;      // the placement's differences away from the site
	int max_diff;  // the read's mismatch budget
};

/* The sites under a refined single-end record (s->seq in read order, CIGAR set), in position
 * order. Returns the index of the first REF or ALT vote (the one that gets a twin), or -1. */
int ks_site_votes(const ks_sites_t *sites, const bntseq_t *bns, const gap_opt_t *opt, const bwa_seq_t *s,
				  int64_t alt_site, std::vector<ks_vote_t> *out);

/* The read, in read order, with REF and ALT swapped at the vote's site. */
void ks_twin_seq(const bwa_seq_t *s, const ks_vote_t &v, std::vector<uint8_t> *read_order);

/* The XV:Z value: one `contig:pos:REF>ALT,allele,away,max_diff,twin_same,twin_mapq` per site,
 * `;`-separated; the twin fields are NA except for votes[twin_of] (twin_same >= 0). */
void ks_format_xv(const bntseq_t *bns, const std::vector<ks_vote_t> &votes, int twin_of, int twin_same,
				  int twin_mapq, std::string *out);

} // namespace bwa_aln

#endif
