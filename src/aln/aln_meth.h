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
	int cert_budget = 16384; // the row budget of a read verified with saturation certificates
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
	// --known-snvs ALT-window search: its windows are candidates for the main search, which takes no
	// stock cuts under --known-snvs, so neither does this one
	int no_cut = 0;
	// plain (non-converted) index: conv_bns == bns and the index text is pac itself (plain aln's
	// no-hit proof; positions of an .amb run may be any base)
	int plain = 0;
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
    int alt = 0; // --known-snvs: aligned with the window's ALTs in place
    int drawn = 0; // a repeat's drawn primary (or a copy's), not a verified row: kept, not counted
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
    uint8_t alt; // --known-snvs: aligned with the window's ALTs in place (in the padding)
};
static_assert(sizeof(aln_meth_hit_t) == 40, "aln_meth_hit_t is a 40-byte record");
// incomplete: lookup budget, an unrepresentable competitor, CIGAR storage, -m work budget
enum { ALN_METH_BIND_BEST = 1, ALN_METH_BIND_UNREP = 2, ALN_METH_BIND_CIGAR = 4, ALN_METH_BIND_M = 8 };
/* A complete read whose row work ended on a saturation certificate (aln_meth_read_t::cert): REPEAT,
 * more than -R loci at the best cost (X0 is a lower bound, MAPQ 0); CROWDED, a unique best and at
 * least ALN_METH_CROWDED retained worse loci (X1 is a lower bound, MAPQ 0). */
enum { ALN_METH_CERT_REPEAT = 1, ALN_METH_CERT_CROWDED = 2 };
enum { ALN_METH_CROWDED = 255 };
// Best loci that give MAPQ 0 under the TAPS/5base tie-break whichever has the fewest conversions:
// stock MAPQ is 0 from 178 worse loci (bwa_approx_mapQ), and all but one are displaced to X1.
enum { ALN_METH_SAT_ROWS = 179 };
// SA rows a read resolves at the flush that binds its lookup budget (at most what is left of it),
// to pick its best-known placement: that flush makes the read incomplete either way.
enum { ALN_METH_BIND_ROWS = 64 };
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
    int ks_alt = 0; // --known-snvs: mate rescue placed the read against its window's ALTs
    long n_over_diff = 0; // rows dropped for more strict differences than the read's -n allows
    long n_over_seed = 0; // rows dropped for more strict differences in the seed than -k allows
    int cert = 0; // ALN_METH_CERT_*: the certificate that ended a complete read's row work, or 0
    // A repeat certificate's draws for the read's duplicate copies (each its own search's draw):
    // the copy's ordinal and its main locus, which is among hits.
    struct copy_draw_t { int64_t ord; int rid, hyp, strand; bwtint_t pos; };
    std::vector<copy_draw_t> copy_draw;
};
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
 * windows the search did not have to propose, resolved (with and without the ALTs) before the read
 * finishes. */
struct aln_meth_candidate_t;
int aln_meth_search_batch(const bwt_t *bwt, const aln_meth_ref_t *ref, const aln_meth_opt_t *mo, int n_seqs,
						  bwa_seq_t *seqs, const gap_opt_t *opt, int n_threads, int read_no, int64_t ord0,
						  int keep_strata, int max_xa_hits, struct aln_dup_state_t *dup, aln_meth_read_t *res,
						  const std::vector<aln_meth_candidate_t> *extra = nullptr);

/* samse under --meth for reads [0, n). Caller must run bwase_initialize once
 * before dispatching any finish workers (as aln_run_se_meth does). */
void aln_meth_finish_se(const aln_meth_ref_t *ref, const aln_meth_opt_t *mo, const gap_opt_t *opt,
						int max_xa_hits, int n, bwa_seq_t *seqs, aln_meth_read_t *res);

/* sampe stages: select SE cost/counts, infer insert size, then pair on retained loci.
 * Stage2 preserves the selected owned CIGAR and applies stock pairing MAPQ and mate rescue.
 * Caller initializes the stock MAPQ table before dispatch; explicit incomplete ends are not promoted. */
void aln_meth_pe_stage1(const gap_opt_t *opt, int n, bwa_seq_t *seqs[2], aln_meth_read_t *res[2]);
int aln_meth_pe_stage2(const aln_meth_ref_t *ref, const aln_meth_opt_t *mo, const gap_opt_t *opt, const pe_opt_t *popt,
					   const isize_info_t *ii, int n, bwa_seq_t *seqs[2], aln_meth_read_t *res[2],
					   uint64_t n_tot[2], uint64_t n_mapped[2]);

/* bwa_print_sam1 under --meth: the stock record without XM:i/XO:i/XG:i, X0/X1 as the strata say,
 * and XR:Z/XG:Z/XM:Z as --meth-tags selects. read_no as for the search. */
void aln_meth_print_sam1(kstring_t *out, const aln_meth_ref_t *ref, const aln_meth_opt_t *mo, bwa_seq_t *p,
						 const aln_meth_read_t *r, int read_no, const bwa_seq_t *mate, int mode, int max_top2);

} // namespace bwa_aln

#endif
