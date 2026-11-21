/* SPDX-License-Identifier: MIT */
/* `aln --known-snvs`: the known SNV sites a run places reads against. New code. */
#ifndef BWA_MEM3_ALN_KNOWN_SNVS_H
#define BWA_MEM3_ALN_KNOWN_SNVS_H

#include "../bntseq.h"
#include "bwtaln.h"
#include "bwape.h"

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
 * covers a site lies inside its window. Longer reads are searched too, but ordinary placement
 * override falls back to stock; methylation genomic refinement and paired rescue have separate bounds. */
#define KS_WINDOW_HALF 160

/* One ALT window: the reference [beg, beg + len) of contig rid with the ALT of one site in place. */
struct ks_window_t {
	int rid;
	int64_t beg;
	int32_t len;
	int64_t site_pos; // the site, 0-based on rid
	uint8_t alt;
};

/* The ALT windows of a site set and their FM-index, built and held entirely in memory (no files,
 * no temporary directory). Window i is contig i of the index. */
struct ks_windows_t {
	std::vector<ks_window_t> win;
	::FMI_search *fmi = nullptr;
	bwt_t *bwt = nullptr;
	const bntseq_t *bns = nullptr;       // the window index's contigs (one per window)
	struct aln_meth_ref_t *mref = nullptr; // --meth: the windows' converted index as a reference
	bool from_cache = false;               // the FM-index was read from --known-snvs-cache
};

/* Builds one window per site allele from `pac` (its stored bases, so N runs carry the main index's
 * fill) and indexes them in memory: the windows as a forward-only reference (ks_windows_ref), the
 * FM-index from libsais, adopted by an FMI_search. Under `meth` the FM-index is built over the
 * converted windows (ks_windows_meth_convert) and `mref` holds the original windows; no
 * original-alphabet window index is built. With `cache_path` the FM-index is read from that file
 * when it was written for the same indexed text (a 128-bit hash of it, its length, the SA sample
 * rate and the mode) and passes its checksum; otherwise it is built as above and the file
 * rewritten (a failed write is a warning). Either way the index is the same. Returns 0, or -1 with
 * an error on stderr. */
int ks_windows_build(const ks_sites_t *sites, const bntseq_t *bns, const uint8_t *pac, int n_threads,
					 int meth, ks_windows_t *out, const char *cache_path = nullptr);

/* The windows as a forward-only reference: bns (one contig `w<i>` per window, no ambiguity
 * holes, fields as bns_restore gives for today's windows.fa) and pac (l_pac/4+1 bytes, as
 * bwa_idx_load_ele loads). Caller owns both (bns_destroy / free). */
void ks_windows_ref(const std::vector<ks_window_t> &win, const bntseq_t *bns, const uint8_t *pac,
					bntseq_t **out_bns, uint8_t **out_pac);
/* --meth: per window `r<name>` (G->A) then `f<name>` (C->T), the layout aln_meth_ref_load checks.
 * Windows only: valid because window text has no ambiguity holes (the 2-bit codes are exact). */
void ks_windows_meth_convert(const bntseq_t *orig_bns, const uint8_t *orig_pac, bntseq_t **conv_bns,
							 uint8_t **conv_pac);

/* --meth: the candidate windows (original reference, rid = contig) of the hits a read's search on
 * the converted window index found that cover their window's site. */
void ks_meth_candidates(const ks_windows_t *w, const bntseq_t *bns, const struct aln_meth_read_t &wres,
						std::vector<struct aln_meth_candidate_t> *out);
void ks_windows_free(ks_windows_t *w);

/* The options a read's ALT-window search runs with: NONSTOP and no -R cut, so synthetic copies
 * and cross-window hits cannot tighten the genomic discovery bound. The main search and the
 * genomic merge still apply the user's retention policy. */
gap_opt_t ks_discovery_opt(const gap_opt_t &opt);

struct aln_opt_t;
struct aln_dup_state_t;
struct aln_meth_read_t;
/* --meth with --known-snvs (o->ks_windows set; without it, the genomic search alone): reads
 * [0, n) searched on the ALT-window index, then on the converted index `bwt` with the windows'
 * candidates, each read's stop certificate (aln_meth_alt_certificate) and its window-search
 * flags as incoming_flags. Incomplete window discovery taints the genomic result: retained
 * windows may place a read, but cannot certify uniqueness or support pair promotion. Arguments
 * and the return value are aln_meth_search_batch's; window_dup is the window search's
 * duplicate state. */
int ks_meth_search(const aln_opt_t *o, const bwt_t *bwt, int n, bwa_seq_t *seqs, int read_no, int64_t ord,
				   int keep, aln_dup_state_t *dup, aln_dup_state_t *window_dup, aln_meth_read_t *res);

/* What ks_override did to a read. */
enum ks_outcome_t {
	KS_UNTOUCHED = 0,   // no window hit covers its site: the stock result stands
	KS_SAME,            // window hits cover a site but change no hit's cost: stock result
	KS_CREDITED,        // placement, counts and MAPQ recomputed with the ALT credited
	KS_INCOMPLETE,      // too many hits to expand, a read too long for its window, or a genomic or
	                    // window search stopped at -m: stock result
	KS_N_OUTCOMES
};

/* Crediting for one single-end read after bwa_cal_pac_pos (stock position and MAPQ set) and
 * before bwa_refine_gapped. `s->aln` holds the main search's hits; win_aln/win_n the same read's
 * hits on the window index. Window hits that cover their window's site join the main hits; per
 * locus the lower search score wins; then c1/c2, the primary (the stock one when still among the
 * best, else a per-read draw, never drand48), MAPQ (bwa_approx_mapQ) and the alternative hits are
 * recomputed with aln's cuts (score <= best + mismatch penalty, differences <= min(max_diff, best's
 * + 1)), applied before the main intervals are expanded. *alt_site is the site allele
 * (aln_site_key) whose ALT the new primary's cost assumes, or -1. */
/* Paired end: `hits` gets the kept hits, best first, with each selected allele in
 * pe_hit_t::alt_site_plus1; it is cleared unless the read is credited. `stopped`: the read's
 * genomic or window search stopped at -m, so either may lack same-score competitors; a read that
 * would be credited then keeps its stock result (KS_INCOMPLETE). */
ks_outcome_t ks_override(const ks_windows_t *w, const bntseq_t *bns, const bwt_t *bwt, const gap_opt_t *opt,
						 int max_xa_hits, bwa_seq_t *s, const bwt_aln1_t *win_aln, int win_n, int64_t *alt_site,
						 std::vector<pe_hit_t> *hits = nullptr, bool stopped = false);
/* Mate rescue's site lister (pe_alt_sites_t::list) over a ks_rescue_ctx_t: the site alleles in the
 * packed range [beg, end), on every contig it reaches, at most KS_RESCUE_MAX_SITES. */
#define KS_RESCUE_MAX_SITES 32
struct ks_rescue_ctx_t {
	const ks_sites_t *sites;
	const bntseq_t *bns;
};
void ks_rescue_sites(const void *ctx, int64_t beg, int64_t end, std::vector<std::pair<int64_t, uint8_t>> *out);
/* bwa_aln2seq for a twin: the same primary choice with its draws from a per-read stream, so the
 * run's drand48 stream never sees the twin. */
void ks_aln2seq_nodraw(const bntseq_t *bns, bwa_seq_t *s);

/* Start tolerance for "the twin lands on the same locus" (refinement can shift a gapped start). */
#define KS_TWIN_SLOP 5

/* Whether [pos, pos + span) is a positive span inside the contig that holds contig_pos. Both
 * positions are forward-reference packed coordinates. */
bool ks_span_in_contig(const bntseq_t *bns, int64_t pos, int64_t span, int64_t contig_pos);
/* Whether a twin's full positive span stays in the original placement's contig and its start
 * is within KS_TWIN_SLOP of the original's. Both positions are forward-reference packed coordinates. */
bool ks_twin_same_contig(const bntseq_t *bns, int64_t twin_pos, int64_t span, int64_t original_pos);
/* As above, with the twin's refined CIGAR span when it has one, otherwise the search's
 * len + ref_shift. */
bool ks_twin_read_same_contig(const bntseq_t *bns, const bwa_seq_t *twin, int64_t original_pos);
/* A paired twin as ks_twin_read_same_contig, and a mapped supporting mate's span inside the
 * contig its start is on. The mate is only bounds-checked (no slop: it is not compared with any
 * original). An end with a CIGAR spans its CIGAR's reference length. Without one, an end pairing
 * moved spans its new hit's span (twin_span, mate_span: bwa_pe_pair's moved_span, -1 when not
 * moved), since stock pairing leaves such an end's ref_shift as its previous hit's; any other end
 * spans len + ref_shift. */
bool ks_twin_pair_same_contig(const bntseq_t *bns, const bwa_seq_t *twin, const bwa_seq_t *mate,
							  int64_t original_pos, int64_t twin_span = -1, int64_t mate_span = -1);

/* Where a twin read lands, without a draw: its main hits (t->aln) and window hits merged as for
 * ks_override. Returns 1 when it has one best hit at (pos, strand) within KS_TWIN_SLOP whose whole
 * reference span lies in the contig holding pos (ks_twin_same_contig), 0 when not (no hit, a tie,
 * elsewhere, or a span off that contig), -1 when incomplete; *mapq is its MAPQ. */
int ks_twin_place(const ks_windows_t *w, const bntseq_t *bns, const bwt_t *bwt, const gap_opt_t *opt,
				  const bwa_seq_t *t, const bwt_aln1_t *win_aln, int win_n, int64_t pos, int strand, int *mapq);

/* One site under a read's final placement. */
struct ks_vote_t {
	const ks_site_t *site;
	int rid;
	int qi;        // index in the reference-oriented read
	char allele;   // 'R', 'A', 'O' (another base) or, under --meth, 'H' (either allele)
	int away;      // the placement's differences away from the site
	int max_diff;  // the read's mismatch budget
	bool conv;     // under --meth: the strand's conversion applies at the site (see ks_site_votes)
};

/* --meth stores a credited site allele in a hit's 8-bit field as its index from the hit's start:
 * ks_site_rel gives 1 + the index of `key` (aln_site_key) among contig rid's site alleles at or
 * after packed position pos, in (position, ALT) order, or 0 (none, not found, or past 255);
 * ks_site_of_rel turns it back into the key, or -1. */
int ks_site_rel(const ks_sites_t *sites, const bntseq_t *bns, int rid, int64_t pos, int64_t key);
int64_t ks_site_of_rel(const ks_sites_t *sites, const bntseq_t *bns, int rid, int64_t pos, int rel);

/* The sites under a refined single-end record (s->seq in read order, CIGAR set), in position
 * order. alt_site: the site allele (aln_site_key) whose ALT the placement's cost assumes, every
 * other site at REF, or -1 for the REF at every site. A vote's `away` leaves out the site's own
 * base. hyp: -1 for a plain record; under --meth the
 * converted strand (1: C->T, 0: G->A), and an allele the conversion can explain on that strand
 * ('H': hidden) gets no twin, as a vote that cannot tell REF from ALT. With cpg (--meth-context
 * cpg; `pac` is then the original reference) the conversion applies only where the aligner allows
 * it, at a CpG: a C followed by a G (C->T), or a G preceded by a C (G->A), with the placement's
 * credited ALT in place and outside an ambiguity run. Elsewhere only literal allele matches count.
 * The difference budget uses clip_len (the search length); full len still indexes emitted bases
 * and soft clips. collapsed frees the scoring model's mirror conversion cell without changing
 * biological REF/ALT classification. With plain pac, differences are recounted from the final
 * CIGAR and credited allele. Returns the index of the first REF or ALT vote (the one that gets a twin), or -1. */
int ks_site_votes(const ks_sites_t *sites, const bntseq_t *bns, const gap_opt_t *opt, const bwa_seq_t *s,
				  int64_t alt_site, int hyp, std::vector<ks_vote_t> *out, const uint8_t *pac = nullptr,
				  int cpg = 0, int collapsed = 0);

/* The quality-trimmed read (clip_len bases), in read order, with REF and ALT swapped at the vote's
 * site (under --meth, hyp >= 0,
 * the swapped-in base as the conversion leaves it unmethylated, where the conversion applies:
 * v.conv). */
void ks_twin_seq(const bwa_seq_t *s, const ks_vote_t &v, int hyp, std::vector<uint8_t> *read_order);

/* The XV:Z value: one `contig:pos:REF>ALT,allele,away,max_diff,twin_same,twin_mapq` per site,
 * `;`-separated; %, :, , and ; in contig names are percent-escaped with uppercase hex.
 * The twin fields are NA except for votes[twin_of] (twin_same >= 0). */
void ks_format_xv(const bntseq_t *bns, const std::vector<ks_vote_t> &votes, int twin_of, int twin_same,
				  int twin_mapq, std::string *out);

} // namespace bwa_aln

#endif
