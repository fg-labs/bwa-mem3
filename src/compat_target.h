/* src/compat_target.h — `--compat` output-compatibility targets.
 *
 * A compat target describes the exact SAM/BAM output shaping needed to
 * reproduce another aligner's byte stream. It is DATA, not control flow: the
 * table in compat_target.cpp holds one row per target, and every consumer reads
 * a field off the selected row instead of testing a flag bit. Adding a target
 * is a new row.
 *
 * A single boolean would not do. bwa and bwa-mem2 diverge from each other on
 * the very fields this shapes -- bwa emits MQ:i and a default @HD, bwa-mem2
 * emits neither -- so "compat" is a choice among targets, not an on/off switch.
 *
 * Compat targets shape OUTPUT ONLY. They never change an alignment, a score, a
 * flag, or any tag's value. See docs/src/whats-different/equivalence.md.
 */

#ifndef BWAMEM3_COMPAT_TARGET_H
#define BWAMEM3_COMPAT_TARGET_H

#ifdef __cplusplus
extern "C" {
#endif

typedef struct compat_target_t {
    /* Canonical spelling, as documented and as reported in diagnostics. */
    const char *name;

    /* One additional accepted spelling, or NULL. */
    const char *alias;

    /* NULL when the target is selectable. Otherwise, WHY it is not: the row is
     * fully specified and unit-tested, but --compat refuses it and prints this
     * instead of a generic "unknown target". For a target whose output shaping
     * we know exactly but whose byte-identity we cannot deliver for other
     * reasons -- see the bwa-mem row. The reason lives here, next to the
     * evidence the row is built from, so the two cannot drift apart. */
    const char *unavailable_reason;

    /* Emit a default @HD when neither -H nor the index sidecar supplies one. */
    int emit_hd;

    /* Exact @HD text (no trailing newline) when emit_hd is set. NULL means
     * "keep whatever default that output path already emits" -- the SAM text
     * and BAM paths currently disagree (fg-labs/bwa-mem3#288), and only the
     * `off` row tolerates that. */
    const char *hd_line;

    /* Honor the bwa-mem3-only <prefix>.hdr / <baseprefix>.dict sidecar. Both
     * upstreams lack the feature entirely (lh3/bwa#348 was closed unmerged), so
     * every non-`off` target must skip it to match. */
    int read_sidecar;

    /* Emit the MQ:i mate mapping quality tag. NOTE: this is not a bwa-mem3
     * invention -- bwa emits it (bwamem.c:935; lh3/bwa#330, merged 2022-03-06)
     * and bwa-mem2 does not, having forked at 0.7.17 before that landed. */
    int emit_mq;

    /* Emit the HN:i hit-count tag. Genuinely bwa-mem3-only: absent from both
     * upstreams. */
    int emit_hn;
} compat_target_t;

/* The `off` row: bwa-mem3's own native output. Never NULL on any mem_opt_t
 * built by mem_opt_init(), so consumers can dereference opt->compat
 * unconditionally. (The `opt0` "was this set explicitly" sentinel in main_mem
 * is memset to zero and is NOT such a struct -- it is only ever read for
 * scalar fields.) */
extern const compat_target_t COMPAT_TARGET_OFF;

/* Resolve a user-supplied --compat value, by canonical name or alias.
 *
 * Returns NULL for an unrecognized name. A name that IS recognized but is
 * unavailable resolves to its row, so the caller can print
 * ->unavailable_reason rather than a generic "unknown target"; callers must
 * check that field before use. */
const compat_target_t *compat_target_from_name(const char *name);

/* Comma-separated list of the selectable target names with their aliases, for
 * usage and error text (e.g. "bwa-mem2 (alias mem2), off"). Derived from the
 * table, so a new row cannot leave it stale. Statically allocated. */
const char *compat_target_selectable_list(void);

/* The table itself, for tests and diagnostics. Sets *n to the row count. */
const compat_target_t *const *compat_targets(int *n);

#ifdef __cplusplus
}
#endif

#endif /* BWAMEM3_COMPAT_TARGET_H */
