/* SPDX-License-Identifier: MIT */
/* New code for bwa-mem3 (not from upstream bwa): methylation-aware `aln --meth`.
 *
 * Converted-space search generates candidate windows within the user's discovery limits.
 * Each exact window is realigned against the original reference with the chemistry matrix.
 * Primary and competitors use stock cost/count accounting on those reward-refined CIGARs.
 * Lookup/work caps mark incomplete results, never an exact unique stratum.
 *
 * Without --meth none of this runs: AlnStock, the stock instantiation, is unchanged. */
#ifndef BWA_MEM3_ALN_METH_H
#define BWA_MEM3_ALN_METH_H

#include <stdint.h>
#include <algorithm>
#include <vector>

#include "../bntseq.h"
#include "../kstring.h"
#include "../meth_opt.h"
#include "aln_bwt.h"
#include "bwtaln.h"
#include "bwape.h"

namespace bwa_aln {

/* The equal-span scorer's matrix classes (aln_meth_realign.h, aln_meth_fast_init): every cell
 * is a match, a mismatch, an ambiguity (either base N) or the conversion cell(s), and each class
 * has one value. ok is false when the scorer does not apply (--meth-context cpg, or another
 * matrix). best_cell is ungapped_is_optimal's largest cell (0 when every cell is negative). */
struct aln_meth_fast_t {
    bool ok = false, collapsed = false;
    int match = 0, mism = 0, amb = 0, conv = 0, best_cell = 0;
};

/* The run's --meth settings. The matrices are at aln's scale (bwa_fill_scmat(1, 3)). */
struct aln_meth_opt_t {
	meth_opt_t m;                  // the shared options (meth_opt.h)
	aln_meth_opt_t() { meth_opt_init(&m); }
	int p_best = 4096;              // legacy budget component; total lookup budget is p_best + p_second
	int p_second = 500;             // second component of the single total budget
	int8_t mat_ot[25], mat_ob[25], mat7_ot[49], mat7_ob[49];
	aln_meth_fast_t fast; // set by aln_meth_opt_finish
};

/* Builds the matrices (meth_fill_mats on aln's base matrix); call once the options are final. */
void aln_meth_opt_finish(aln_meth_opt_t *mo);

/* The converted index's bns (contigs r0,f0,r1,f1,...: rid>>1 is the original contig, rid&1 is 1
 * for the C->T copy, OT, and 0 for the G->A copy, OB) and the original reference's bns and pac. */
struct aln_meth_ref_t {
	const bntseq_t *conv_bns = nullptr;
	bntseq_t *bns = nullptr;
	uint8_t *pac = nullptr;
	const struct ks_sites_t *known_snvs = nullptr; // --known-snvs sites (null without); not owned
	// --known-snvs ALT windows as a reference: contig i's site, as an offset in it (else empty)
	std::vector<int32_t> window_site;
};
int aln_meth_ref_load(aln_meth_ref_t *r, const char *orig_prefix, const bntseq_t *conv_bns, int n_threads);
/* Takes ownership of an in-memory original reference (bns, pac) and checks it pairs with
 * conv_bns exactly as aln_meth_ref_load does. Returns 0, or -1 with the same error (the reference
 * is freed). */
int aln_meth_ref_adopt(aln_meth_ref_t *r, bntseq_t *orig_bns, uint8_t *orig_pac, const bntseq_t *conv_bns);
void aln_meth_ref_free(aln_meth_ref_t *r);

/* Reject user @SQ records before any output/index work; returns 0 or -1 with an error. */
int aln_meth_header_check(const char *rg_line);
/* The --meth header: the shared header text (meth_hdr_text) with aln's @HD, the -r text and
 * aln's @PG, then the second @PG. orig_hdr_lines: the original reference's sidecar, or NULL.
 * Returns 0 on success or -1 on invalid user header text or allocation failure. */
int aln_meth_header(kstring_t *out, const aln_meth_ref_t *ref, const char *rg_line, const char *pg_line,
					const char *orig_hdr_lines);

/* The converted read: read 1 (read_no 1, also single-end) C->T, read 2 (read_no 0) G->A, per base,
 * so any orientation of the read converts the same way. */
void aln_meth_convert(int read_no, const ubyte_t *in, int len, ubyte_t *out);

/* Owned original-reference alignment. Reward includes raw global terminal-gap penalties;
 * display CIGAR uses stock terminal normalization. All CIGAR vectors use BAM op codes.
 * raw_cigar is used while the search ranks and accounts windows, and by mate and cap rescue on
 * the alignments they build. Finishing a read stores its loci as aln_meth_hit_t, without it. */
struct aln_meth_alignment_t {
    int reward = 0;
    int64_t rb = 0, re = 0;
    int band = 0;
    std::vector<uint32_t> raw_cigar, cigar;
};
/* A final locus is (original contig, hypothesis, strand, normalized start). All windows
 * reaching it compete by cost, then neutral TAPS conversion count when enabled,
 * then normalized CIGAR/end and raw CIGAR for deterministic ties.
 * Counts/cost describe the raw scored alignment, before terminal display normalization. */
struct aln_meth_occ_t {
    bwtint_t pos;
    int rid, strand, hyp;
    int n_mm = 0, n_gapo = 0, n_gape = 0, n_ins = 0, n_del = 0, n_conv = 0;
    int n_mirror = 0;      // mismatches in the mirror cell (OT read C x ref T, OB read G x ref A); in n_mm
    int alt = 0;           // --known-snvs: the credited site allele as ks_site_rel gives it (0: none)
    int64_t ks_site = -1;  // --known-snvs: the credited site allele (aln_site_key), or -1
    int64_t cost = 0;
    aln_meth_alignment_t alignment;
};
/* A finished read's stored locus: what pairing, rescue and output read from an aln_meth_occ_t,
 * packed (40 bytes, not 136 plus a CIGAR allocation). pos is the locus start, which is the
 * alignment's start (rb); span is its reference length (re - rb). The display CIGAR is cig
 * itself when it has one operation (nearly every locus), else n_cig operations at offset cig in
 * the read's cig table; see aln_meth_hit_cigar. account_alignment refuses a cost or count that
 * does not fit; storing a locus (meth_hit) is fatal on any field that does not, on pos != rb, and
 * on a cig table past 2^32 entries. */
struct aln_meth_hit_t {
    bwtint_t pos;
    int32_t span, cost, rid, n_conv;
    uint32_t cig, n_cig;
    uint8_t strand, hyp, n_mm, n_gapo, n_gape;
    uint8_t alt; // --known-snvs: the credited site allele, ks_site_rel from pos (0: none; in the padding)
    uint8_t n_mirror; // the locus's mirror-cell mismatches (counted in n_mm; in the padding)
};
static_assert(sizeof(aln_meth_hit_t) == 40, "aln_meth_hit_t is a 40-byte record");
/* --meth-mapq conversion-blind: a locus at its lowest mirror-free row (the row collapsed scoring
 * keeps there), at that row's mirror-free cost and mismatches, its gaps and its span. rid and hyp
 * with strand and pos give collapsed's order of its loci (pairing indexes them, the draw of the
 * main hit walks them). */
struct aln_meth_blind_row_t {
    bwtint_t pos;
    int32_t cost, rid;
    uint16_t span; // a row never spans more than ALN_MAX_READ_LEN (0xffff): longer windows are dropped
    uint8_t strand, hyp, n_mm, n_gapo, n_gape;
    uint8_t alt; // --known-snvs: the row's credited site allele, ks_site_rel from pos (0: none)
};
static_assert(sizeof(aln_meth_blind_row_t) == 24, "aln_meth_blind_row_t is a 24-byte record");
// incomplete: lookup budget, an unrepresentable competitor, CIGAR storage, -m work budget
enum { ALN_METH_BIND_BEST = 1, ALN_METH_BIND_UNREP = 2, ALN_METH_BIND_CIGAR = 4, ALN_METH_BIND_M = 8 };
// A same-read ALT-index lower bound; an unavailable bound must never authorize stopping.
struct aln_meth_alt_bound_t {
    int value = 0;
    bool valid = false;
};
struct aln_meth_read_t {
    // The search's candidate loci; finishing the read moves the retained ones to hits and frees it.
    std::vector<aln_meth_occ_t> occ;
    int main = -1;
    int c1 = 0, c2 = 0;
    // Equal-cost loci displaced from c1 to c2, and their nearest conversion disadvantage.
    int n_conv_moved = 0, conv_gap = 0;
    int flags = 0, hyp = -1;
    // Search work (aln_meth_report); zero on a duplicate that copied another read's search.
    long n_lookups = 0, n_wrong_copy = 0, n_rescored = 0, n_flushes = 0;
    long n_unrepresentable = 0; // refined windows whose counts do not fit a record
    // A finished read's best and retained worse-cost loci, canonical order (main indexes these).
    // An incomplete read (flags != 0) keeps only its best-cost loci, and only its main hit once
    // nothing redraws over them (no duplicate copies left to draw); c1/c2 still count every
    // retained locus. cig holds the CIGARs of more than one operation. Declared last, so the
    // search's fields above keep their offsets.
    std::vector<aln_meth_hit_t> hits;
    std::vector<uint32_t> cig;
    int64_t ks_site = -1; // --known-snvs: the site allele mate rescue credited (aln_site_key), or -1
    aln_meth_alt_bound_t full_width_bound; // captured before search shadowing, only when requested
    // With full_width_bound, on a reference with window_site: the fewest search mismatches of a
    // proposed read-length window that covers its own window's site, and at most one more than
    // the mismatches up to which an unbound search proposes every such window (its difference
    // limit, and the seed's when seeded). Read it through aln_meth_alt_certificate.
    aln_meth_alt_bound_t site_bound;
    long n_over_diff = 0; // rows dropped for more strict differences than the read's -n allows
    long n_over_seed = 0; // rows dropped for more strict differences in the seed than -k allows
    // --meth-mapq conversion-blind: the strata with the mirror cell freed (cost minus s_mm per
    // mirror mismatch), counted over every candidate before the stock cut: the best mirror-free
    // cost, the loci at it, and the loci up to s_mm above it. blind_best is INT64_MAX when unset.
    int64_t blind_best = INT64_MAX;
    int blind_c1 = 0, blind_c2 = 0;
    // The mirror-free strata's loci that are not hits: the genomic cut drops them (a locus that
    // costs two or more mirror mismatches more than the genomic best), or the -n (or -k) limit
    // drops every row of theirs, but not once its mirror mismatches are free (also for a read
    // with no hit at all). For the paired-end shadow pairing, in collapsed's order.
    std::vector<aln_meth_blind_row_t> blind_rows;
    // The hits whose locus's lowest mirror-free row is not the hit's own row (the genomic search
    // keeps the genomically cheapest row of a locus, and on equal cost the lower CIGAR), at that
    // row; pos and strand name the hit, sorted by them. Read through aln_meth_blind_view.
    std::vector<aln_meth_blind_row_t> blind_alt;
    // blind_main: the index, among the mirror-free strata's best-cost loci (blind_best) in
    // collapsed's order (rid, hyp, strand, pos), of the one collapsed scoring draws as the main
    // hit (the same draw over the same loci with the read's seed), and blind_pick that locus. For
    // the paired-end shadow pairing, which starts each end where collapsed does; drawn only for a
    // read that pairs or has duplicate copies (otherwise blind_main is -1). blind_ties holds the
    // loci themselves only while a read's duplicate copies draw again over them, each with its own
    // seed (aln_meth_copy_dup); aln_meth_search_batch frees them once the copies are drawn.
    std::vector<aln_meth_blind_row_t> blind_ties;
    aln_meth_blind_row_t blind_pick = aln_meth_blind_row_t();
    int blind_main = -1;
    // Where collapsed scoring's pairing and mate rescue leave this end, and its MAPQ there (0: it
    // leaves the end unmapped; blind_place_mapq -1: not known, see aln_meth_pe_stage2).
    bwtint_t blind_place_pos = 0;
    int blind_place_strand = 0, blind_place_mapq = -1;
    // The mirror table, while the read is searched: per locus, the lowest mirror-free row among
    // its rows with a mirror mismatch that are within the read's limits once those are free (a
    // read C over a reference T is a mismatch genomically but a conversion mirror-free), CIGARs
    // kept for the tie-break. It only observes the search: nothing the search reads depends on it.
    // Consumed and cleared when the read finishes.
    std::vector<aln_meth_occ_t> mirror_rows;
    // The mirror table outgrew its CIGAR budget: the mirror-free strata are not known, so each
    // end's conversion-blind cap is 0.
    int blind_partial = 0;
    // The MAPQ cap the record's end carries to emission (255: none), and the MAPQ it replaced
    // when the cap lowered it (YQ:i), else -1.
    int blind_cap = 255, yq = -1;
};
/* The stop certificate a read's ALT-window search gives its genomic search: a lower bound on the
 * mismatches of every ungapped placement that covers a site with that site's ALT in place. A
 * search that ran to completion (no flags) proposed every such window within its limits, so its
 * site_bound holds; a bound search may have missed one, and only the index-wide width bound
 * (which also counts copies beside the site) is left.
 * The fallback is defensive. When the window search has flags, ks_meth_search passes the same
 * flags to the genomic search as incoming_flags, and a slot with flags is never stop-eligible
 * (aln_meth_stop_eligible), so the bound returned here is not used for stopping. It then only
 * keys duplicate sharing (copies share a search only with equal bounds), and it is the whole
 * certificate for a window reference without window_site. */
static inline aln_meth_alt_bound_t aln_meth_alt_certificate(const aln_meth_read_t &window)
{
    return !window.flags && window.site_bound.valid ? window.site_bound : window.full_width_bound;
}
// --meth-mapq conversion-blind: hit h of r at its locus's lowest mirror-free row (blind_alt), or,
// when that is h's own row, at h with its mirror mismatches free.
static inline aln_meth_blind_row_t aln_meth_blind_view(const aln_meth_read_t &r, const aln_meth_hit_t &h, int s_mm)
{
    auto before = [](const aln_meth_blind_row_t &x, const aln_meth_blind_row_t &y) {
        return x.pos != y.pos ? x.pos < y.pos : x.strand < y.strand;
    };
    aln_meth_blind_row_t b;
    b.pos = h.pos; b.strand = h.strand; b.rid = h.rid; b.hyp = h.hyp;
    const auto it = std::lower_bound(r.blind_alt.begin(), r.blind_alt.end(), b, before);
    if (it != r.blind_alt.end() && it->pos == h.pos && it->strand == h.strand) return *it;
    b.span = (uint16_t)h.span;
    b.cost = (int32_t)((int64_t)h.cost - (int64_t)s_mm * h.n_mirror); // in [0, h.cost]
    b.n_mm = (uint8_t)(h.n_mm - h.n_mirror); b.n_gapo = h.n_gapo; b.n_gape = h.n_gape;
    return b;
}
// The pointer is into h itself (one operation) or into r.cig, so it is valid while both are.
static inline const uint32_t *aln_meth_hit_cigar(const aln_meth_read_t &r, const aln_meth_hit_t &h)
{
    return h.n_cig == 1 ? &h.cig : r.cig.data() + h.cig;
}

/* --meth's per-slot search state (AlnMeth::ext_t): see aln_meth.cpp. */
struct aln_meth_group_t;

/* Searches every read of a batch under --meth (the counterpart of aln_search_batch: the same
 * threads, read claiming, interleaved slots and duplicate handling). read_no: 1 for single-end
 * and read 1, 0 for read 2. The reads' ordinals in their stream are ord0 + i (each read's main-hit
 * draw is seeded from it). res[i] gets read i's result. keep_strata: keep every locus of both
 * strata (paired-end pairs on them; an incomplete read, which never pairs, keeps only its main
 * hit); otherwise a read keeps its main hit, and the other loci only when they will be printed as
 * XA (at most max_xa_hits + 1 occurrences). extra (--known-snvs, or null): per read, candidate
 * windows the search did not have to propose, resolved (as the reference and with each site
 * allele in place) before the read finishes; identical reads share a search only with equal extras. */
struct aln_meth_candidate_t;
/* Optional alt_bounds are same-read ALT-index certificates, never genomic widths.
 * capture_full_width records width[len-1] before search mutation (and site_bound, on a reference
 * with window_site); both defaults are off.
 * incoming_flags are per-read window-search flags, applied after reset and before searching or
 * finalizing; duplicate sharing requires equal flags. Null means no incoming flags. */
int aln_meth_search_batch(const bwt_t *bwt, const aln_meth_ref_t *ref, const aln_meth_opt_t *mo, int n_seqs,
						  bwa_seq_t *seqs, const gap_opt_t *opt, int n_threads, int read_no, int64_t ord0,
						  int keep_strata, int max_xa_hits, struct aln_dup_state_t *dup, aln_meth_read_t *res,
						  const std::vector<aln_meth_candidate_t> *extra = nullptr,
						  const aln_meth_alt_bound_t *alt_bounds = nullptr, bool capture_full_width = false,
                          const int *incoming_flags = nullptr);

/* samse under --meth for reads [0, n). Caller must run bwase_initialize once
 * before dispatching any finish workers (as aln_run_se_meth does). */
void aln_meth_finish_se(const aln_meth_ref_t *ref, const aln_meth_opt_t *mo, const gap_opt_t *opt,
						int max_xa_hits, int n, bwa_seq_t *seqs, aln_meth_read_t *res);

/* sampe stages: select SE cost/counts, infer insert size, then pair on retained loci.
 * Stage2 preserves the selected owned CIGAR and applies stock pairing MAPQ and mate rescue.
 * Caller initializes the stock MAPQ table before dispatch; explicit incomplete ends are not promoted. */
void aln_meth_pe_stage1(const gap_opt_t *opt, int n, bwa_seq_t *seqs[2], aln_meth_read_t *res[2]);
/* --meth-mapq conversion-blind, after stage 1: the insert size a collapsed run would infer, from
 * each end at collapsed's single-end selection (aln_meth_read_t::blind_main) with collapsed's
 * single-end MAPQ, with last_ii the previous batch's such estimate (as bwa_pe_isize_ref). An end
 * whose mirror-free strata are unknown (blind_partial) counts as MAPQ 0; one with none keeps its
 * genomic selection. The ends are left as they were; nothing is reported. */
void aln_meth_blind_isize(const aln_meth_ref_t *ref, const gap_opt_t *gopt, const pe_opt_t *popt, int n,
                          bwa_seq_t *seqs[2], aln_meth_read_t *const res[2], const uint8_t *excluded,
                          const isize_info_t *last_ii, isize_info_t *ii);
/* blind_ii: under --meth-mapq conversion-blind, the insert size aln_meth_blind_isize gave, for the
 * pairing a collapsed run makes (null: ii). */
int aln_meth_pe_stage2(const aln_meth_ref_t *ref, const aln_meth_opt_t *mo, const gap_opt_t *opt, const pe_opt_t *popt,
					   const isize_info_t *ii, int n, bwa_seq_t *seqs[2], aln_meth_read_t *res[2],
					   uint64_t n_tot[2], uint64_t n_mapped[2], const isize_info_t *blind_ii = nullptr);

/* bwa_print_sam1 under --meth: the stock record without XM:i/XO:i/XG:i, X0/X1 as the strata say,
 * and XR:Z/XG:Z/XM:Z as --meth-tags selects. read_no as for the search. */
void aln_meth_print_sam1(kstring_t *out, const aln_meth_ref_t *ref, const aln_meth_opt_t *mo, bwa_seq_t *p,
						 const aln_meth_read_t *r, int read_no, const bwa_seq_t *mate, int mode, int max_top2);

} // namespace bwa_aln

#endif
