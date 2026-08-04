/* src/compat_target.cpp — the `--compat` target table.
 *
 * Every field below is derived from the upstream source cited beside it,
 * verified against lh3/bwa v0.7.19 (b92993c, 0.7.19-r1273) and bwa-mem2 v2.2.1.
 * test/unit/test_compat_target.cpp asserts each field, so a transcription error
 * cannot silently rot into a target that claims a parity it does not deliver.
 */

#include "compat_target.h"

#include <string>
#include <string.h>

/* --- `off`: bwa-mem3's native output -------------------------------------
 *
 * hd_line is the one canonical default (#288). It used to be NULL, meaning
 * "whatever that output path emits" -- which was a different string on the
 * SAM-text and BAM paths. */
const compat_target_t COMPAT_TARGET_OFF = {
    /* .name               */ "off",
    /* .alias              */ NULL,
    /* .unavailable_reason */ NULL,
    /* .emit_hd            */ 1,
    /* .hd_line            */ BWAMEM3_DEFAULT_HD_LINE,
    /* .read_sidecar       */ 1,
    /* .emit_mq            */ 1,
    /* .emit_hn            */ 1,
};

/* --- `bwa-mem2`: bwa-mem2 v2.2.1 -----------------------------------------
 *
 * @HD:      none. bwa-mem2's bwa_print_sam_hdr (src/bwa.cpp:523) has no @HD
 *           logic at all -- no n_HD counter, no default emission. bwa grew one
 *           in 0.7.18 (6b18630), after bwa-mem2 forked at 0.7.17.
 * sidecar:  absent; `grep -n '\.hdr\|\.dict' src/bwa.cpp` is empty.
 * MQ:i:     absent; lh3/bwa#330 post-dates the fork.
 * HN:i:     absent; bwa-mem3-only.
 *
 * Note @SQ AH:* is NOT a field here: bwa-mem2 emits it on generated @SQ
 * (src/bwa.cpp:538, and :545 in the compiled non-ORIG branch), exactly as we
 * do, so it is unconditional rather than a divergence to toggle. */
static const compat_target_t COMPAT_TARGET_BWA_MEM2 = {
    /* .name               */ "bwa-mem2",
    /* .alias              */ "mem2",
    /* .unavailable_reason */ NULL,
    /* .emit_hd            */ 0,
    /* .hd_line            */ NULL,
    /* .read_sidecar       */ 0,
    /* .emit_mq            */ 0,
    /* .emit_hn            */ 0,
};

/* --- `bwa-mem`: bwa 0.7.19 -------------------------------------------------
 *
 * @HD:      "@HD\tVN:1.5\tSO:unsorted\tGO:query" (bwa.c:426) -- byte-identical
 *           to bwa-mem3's own default, so this row shares the constant.
 * sidecar:  absent; lh3/bwa#348 was closed unmerged 2025-03-21.
 * MQ:i:     PRESENT (bwamem.c:935) -- this is the field a single compat
 *           boolean could not express, and the reason this table exists.
 * HN:i:     absent; bwa-mem3-only.
 *
 * Pinned at 0.7.19, and the pin matters: bwa gained both the default @HD
 * (6b18630) and MQ:i (d8dd308) in 0.7.18, so a 0.7.17 target would be this
 * row's opposite on those two fields -- and identical to the bwa-mem2 row.
 *
 * This row was unselectable until 0.9.0, on the strength of a measurement
 * (224/63583 differing records) that was retracted as having been taken
 * against a mis-pinned baseline binary. Re-measured, bwa 0.7.19, bwa-mem2
 * v2.2.1 and bwa-mem3 agree byte-for-byte on alignment records once the
 * additive tags are removed. Auditing 0.7.17..0.7.19 says that is structural
 * rather than incidental: the only output-affecting changes in that range are
 * the two additive items above plus the opt-in -u/-z flags -- nothing touches
 * seeding, chaining, extension, pairing, MAPQ or dedup. So bwa-mem2's
 * "identical to bwa-mem 0.7.17" carries forward to 0.7.19, and bwa-mem2 parity
 * is bwa parity.
 *
 * The limit worth knowing: --compat shapes output and never moves an
 * alignment, so on any record where bwa and bwa-mem2 disagree on an alignment
 * field, at most one of these two targets can be byte-identical. No such
 * record has been observed. */
static const compat_target_t COMPAT_TARGET_BWA_MEM = {
    /* .name               */ "bwa-mem",
    /* .alias              */ NULL,
    /* .unavailable_reason */ NULL,
    /* .emit_hd            */ 1,
    /* .hd_line            */ BWAMEM3_DEFAULT_HD_LINE,   /* == bwa.c:426 */
    /* .read_sidecar       */ 0,
    /* .emit_mq            */ 1,
    /* .emit_hn            */ 0,
};

/* Ordered as the diagnostics should read: real targets first, `off` last. */
static const compat_target_t *const COMPAT_TARGETS[] = {
    &COMPAT_TARGET_BWA_MEM2,
    &COMPAT_TARGET_BWA_MEM,
    &COMPAT_TARGET_OFF,
};
static const int COMPAT_TARGET_N =
    (int)(sizeof(COMPAT_TARGETS) / sizeof(COMPAT_TARGETS[0]));

const compat_target_t *compat_target_from_name(const char *name)
{
    if (name == NULL) return NULL;
    for (int i = 0; i < COMPAT_TARGET_N; ++i) {
        const compat_target_t *t = COMPAT_TARGETS[i];
        if (strcmp(name, t->name) == 0) return t;
        if (t->alias != NULL && strcmp(name, t->alias) == 0) return t;
    }
    return NULL;
}

const char *compat_target_selectable_list(void)
{
    /* Derived from the table rather than written out as a literal, so adding a
     * selectable row cannot leave the usage and error text stale. Function-local
     * static: built once, thread-safely, on first use. The returned pointer is
     * valid for the life of the process and must not be freed -- worth stating
     * because the C linkage on this function invites the opposite assumption. */
    static const std::string list = [] {
        std::string s;
        for (int i = 0; i < COMPAT_TARGET_N; ++i) {
            const compat_target_t *t = COMPAT_TARGETS[i];
            if (t->unavailable_reason != NULL) continue;
            if (!s.empty()) s += ", ";
            s += t->name;
            if (t->alias != NULL) { s += " (alias "; s += t->alias; s += ")"; }
        }
        return s;
    }();
    return list.c_str();
}

const compat_target_t *const *compat_targets(int *n)
{
    if (n != NULL) *n = COMPAT_TARGET_N;
    return COMPAT_TARGETS;
}
