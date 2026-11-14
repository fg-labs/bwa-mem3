/* alt_lift.h -- ALT-contig -> primary-assembly liftover for `mem --alt-project`.
 *
 * The `.alt` sidecar lists each ALT contig as a SAM record aligned to the primary
 * assembly (bwakit format): QNAME = ALT contig, FLAG 16 = the contig is the reverse
 * complement of the primary sequence it aligns to, RNAME/POS/CIGAR = where. bns
 * only reads the names (to set is_alt); this module keeps the alignment so an
 * interval on an ALT contig can be mapped onto the primary assembly.
 *
 * The table is indexed by bns rid and is separate from bntann1_t, so the index
 * structs and the shared-memory layout are untouched. */
#ifndef ALT_LIFT_H
#define ALT_LIFT_H

#include <stdint.h>
#include <stddef.h>
#include "bntseq.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    int32_t   trid;      /* primary-assembly rid this ALT contig aligns to; -1 = none */
    int32_t   rev;       /* 1 if the contig is the reverse complement (FLAG 16) */
    int64_t   tpos;      /* 0-based start of the alignment on trid */
    int64_t   tlen;      /* length of contig trid (lifted intervals are clamped to it) */
    int32_t   qlen;      /* query length the CIGAR consumes (M/I/S/H/=/X) */
    int32_t   n_cigar;
    uint32_t *cigar;     /* BAM-style: len << 4 | op (0 M,1 I,2 D,3 N,4 S,5 H,7 =,8 X) */
    /* Per-op prefix sums, n_cigar + 1 entries each: query and reference bases
     * before op i. Two binary searches on qoff find the ops an interval starts
     * and ends in; only the ops between them are walked. */
    int64_t  *qoff, *roff;
} alt_lift1_t;

/* The .alt lines alt_lift_load skipped, by reason (the startup warning prints
 * these counts). */
typedef struct {
    long no_ctg;       /* QNAME is not a contig of the index */
    long no_target;    /* RNAME is not a contig of the index */
    long alt_target;   /* RNAME is itself an ALT contig */
    long unmapped;     /* FLAG 0x4, POS 0, RNAME or CIGAR '*' */
    long not_primary;  /* FLAG 0x100 or 0x800 */
    long dup;          /* a later line for a contig that already has one */
    long bad;          /* malformed: too few fields, FLAG/POS not a whole number, POS
                        * past the target's end, a bad CIGAR, a CIGAR length that is
                        * not the contig's, or an alignment that does not end within
                        * the target. A line is counted under the first reason
                        * that applies, tested in this order: too few fields,
                        * no_ctg, FLAG/POS, unmapped, not_primary, no_target,
                        * alt_target, dup, POS past the end, then the CIGAR */
} alt_lift_skips_t;

typedef struct alt_lift_t {
    int32_t          n_seqs;  /* == bns->n_seqs */
    int32_t          n_lift;  /* ALT contigs with a usable alignment */
    alt_lift1_t     *a;       /* n_seqs entries, indexed by rid */
    alt_lift_skips_t skip;    /* skipped lines, by reason */
} alt_lift_t;

/* The result of lifting one interval (see alt_lift_interval). */
typedef struct {
    int     trid;         /* target rid */
    int     rev;          /* 1 if the strand flips */
    int64_t ps, pe;       /* [ps, pe): 0-based forward interval on trid, pe > ps */
    int64_t n_aln;        /* query bases of the interval in M/=/X ops */
    int64_t n_clip;       /* query bases of the interval in S/H ops (not aligned at all) */
    int64_t max_diag;     /* widest |inserted - deleted| bases between the interval's
                           * start and any point inside it: the band an alignment of
                           * the interval to [ps, pe) needs to follow the ALT's indels */
    int     edge_ins;     /* 1 if the interval's first or last base lies in an I op */
} alt_lift_iv_t;

/* Parse `path` against `bns`. Returns NULL (errno set) if the file cannot be
 * opened. Usable lines are primary alignments (no FLAG 0x4/0x100/0x800) of a
 * contig in bns onto a contig in bns that is not itself an ALT, with a CIGAR
 * whose query length equals the contig length and that ends within the target;
 * when a contig has several usable lines the first wins. A malformed line (fewer
 * than 6 fields, a FLAG or POS that is not a non-negative integer, a bad CIGAR, a
 * length mismatch, or an alignment past the target's end) is named in a warning
 * with its reason, the first 10 of them individually, and one summary warning
 * counts every skipped line by reason (also returned in ->skip). */
alt_lift_t *alt_lift_load(const char *path, const bntseq_t *bns);
void alt_lift_destroy(alt_lift_t *l);

/* BED rows and otherwise-usable alignments excluded by restricted loading. */
typedef struct {
    long input_intervals;       /* nonempty validated rows before merge */
    long merged_intervals;
    long zero_width_intervals;
    long excluded;              /* excluded before first-retained selection */
} alt_lift_exclude_stats_t;

/* Load with whole-primary-footprint BED exclusion, preserving bns ALT flags.
 * stats and error are required, with error_size > 0, and initialized on entry.
 * Returns NULL with a bounded diagnostic on input failure. A successful table
 * may have zero retained/excluded entries; the CLI decides whether to refuse it. */
alt_lift_t *alt_lift_load_excluded(const char *path, const bntseq_t *bns,
                                  const char *bed_path,
                                  alt_lift_exclude_stats_t *stats,
                                  char *error, size_t error_size);

/* Lift the half-open interval [s, e) of ALT contig `rid` (0-based, forward contig
 * coordinates) onto the primary assembly. Returns 0 and fills *out, or -1 if rid
 * has no alignment or the interval is empty or outside the contig. Bases inside an
 * ALT insertion or clip map to the adjacent primary position, so a caller that
 * needs the lift to describe the interval checks out->n_aln, out->n_clip and
 * out->edge_ins. */
int alt_lift_interval(const alt_lift_t *l, int rid, int64_t s, int64_t e, alt_lift_iv_t *out);

#ifdef __cplusplus
}
#endif

/* --alt-project: lift the ALT-contig regions of one read onto the primary
 * assembly through opt->alt_lift, drop the regions a lift makes redundant (see
 * alt_lift.cpp), and leave the survivors in by-score order (score descending,
 * rb, qb, as mem_alnreg_score_lt; then the lower contig index before the lift,
 * an ALT contig's for a lifted region). Which of two tied lifted regions at one
 * locus survives is decided first by whether each can replace the
 * primary-assembly regions there (see alt_lift.cpp), and only then by that
 * order. Returns the new region count. No-op when no region lifts. C++ only,
 * like bwamem.h. A positive opt->alt_project_min_gap first refuses a lift whose
 * score advantage over a different-locus, overlapping-query primary or eligible
 * projected hit is <= that gap times -A. Decisions use one projected snapshot;
 * refused hits retain their original ALT region. Zero disables the gate.
 * The developer-only alt_project_pair_gap alternative bypasses this gate and
 * instead caps confident paired selection downstream, retaining all candidates. */
#ifdef __cplusplus
struct mem_opt_t;
struct mem_alnreg_t;
/* Bounded developer diagnostics after rescue/primary marking and pair selection.
 * Each call has its own 64-invocation budget, separate from projection tracing. */
void mem_alt_trace_pair(const bntseq_t *bns, const char *stage,
                        const char *name, int read_id, int n, const mem_alnreg_t *a,
                        int n_primary, int selected = -1, int mapq = -1);
int mem_alt_project(const struct mem_opt_t *opt, const bntseq_t *bns, int n, struct mem_alnreg_t *a,
                    const char *read_name = NULL, int read_id = -1);

/* --alt-project primary: remove every region of a[0..n) still on an ALT contig
 * (those mem_alt_project did not lift), keeping the others in order, and return
 * the new count. Before removing them, set *oa to the best one scoring at least
 * opt->T -- highest score, then lowest rid, rb and qb, so the choice does not
 * depend on where the read sits in the input; when there is none, only
 * oa->rid is written, to -1. *oa is marked primary, with sub and sub_n
 * computed among the read's other ALT regions as primary marking would
 * (mask_level overlap), so mem_reg2aln gives it the MAPQ it would have had if
 * those were its only hits. Only reads whose record ends up unmapped use it
 * (OA:Z). */
int mem_alt_drop_unlifted(const struct mem_opt_t *opt, int n, struct mem_alnreg_t *a,
                          struct mem_alnreg_t *oa);
#endif

#endif
