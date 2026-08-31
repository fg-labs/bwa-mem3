#ifndef READ_MEMO_H
#define READ_MEMO_H

#include <cstdint>

/* Whole-read-pair memoization ("--dedup-reads").
 *
 * When two input read-PAIRS are byte-identical in both mates' bases, their
 * post-extension alignment regions (mem_alnreg_v) are byte-identical too (proven
 * by the Phase-0 kernel-invariance audit: the alignment kernels are a pure
 * function of read content + index; every id/name/qual/position dependence lives
 * in the replayed SAM stage). So we can run seed->chain->extend once per distinct
 * pair (the REPRESENTATIVE) and copy its regs to the DUPLICATES, then replay the
 * SAM stage per read for byte-identical output.
 *
 * This header exposes only what bwamem.cpp / fastmap.cpp need. The dedup window
 * is one align invocation (the -K chunk): duplicates are collapsed within a
 * chunk only, because the alignment result depends on chunk-scoped state
 * (est_insert_high from per-chunk mem_pestat), which is constant within a chunk
 * but varies across chunks -- so a cross-chunk cache would not be byte-identical.
 *
 * PHASE 1 (this): the CLI knob, the env config, the fingerprint pre-pass, and the
 * net-cycles controller are all wired, but NOTHING consumes role[]/rep_pair[] --
 * the memoize path is added in Phase 2. Phase 1 is byte-identical by construction
 * and exists to (a) land the surface and (b) confirm the probe cost is ~0 on
 * low-duplicate workloads (WGS).
 *
 * The design mirrors the shipped extension-DP dedup controller (BWAMEM3_DEDUP /
 * --dedup, "#415") deliberately: same off|on|auto surface, same net-cycles
 * self-calibrating controller (no user-facing dup-rate threshold), same
 * ends-only word-fingerprint + full-byte-verify, same env-knob family.
 */

#include "bwa.h"        /* bseq1_t (anonymous-struct typedef; cannot be fwd-declared) */
struct mem_opt_t;

/* Resolved mode. */
enum { READMEMO_OFF = 0, READMEMO_ON = 1, READMEMO_AUTO = 2 };

/* Resolve the mode from --dedup-reads (mode_arg) > BWAMEM3_DEDUP_READS env >
 * default "auto". Fatal (exit 1) on an unrecognized value. An empty CLI value is
 * rejected earlier by the getopt handler; an explicitly-empty env value is fatal
 * (cannot silently fall back to the default). Mirrors mem_dedup_configure. */
void mem_dedup_reads_configure(const char *mode_arg);

/* The resolved mode (lazy env/default resolution for entry points -- unit
 * binaries, library use -- that never call mem_dedup_reads_configure). */
int read_memo_mode(void);

/* Per-align-invocation working state. Owned by the module and grown-in-place
 * across invocations (sized from the pair count). Phase 1 fills role[]/rep_pair[]
 * but no consumer reads them. */
struct read_memo_state {
    int      npairs;       /* n/2 for the current invocation                 */
    uint8_t *role;         /* [npairs] 0=NONE (odd tail), 1=REP, 2=DUP        */
    int32_t *rep_pair;     /* [npairs] DUP pair -> its REP pair index; else -1 */
    int64_t  cap;          /* allocated capacity (pairs)                      */
};

/* Result of one pre-pass, fed to the controller. */
struct read_memo_result {
    int64_t  pairs;        /* pairs examined this invocation      */
    int64_t  dup_pairs;    /* pairs marked DUP (collapsible)       */
    uint64_t probe_ns;     /* measured pre-pass wall cost (ns)     */
};

/* Serial fingerprint pre-pass over the chunk's pairs (seqs[2i], seqs[2i+1]),
 * run on RAW (pre-kernel1) bases. Fills st->role / st->rep_pair (grows st as
 * needed) and returns counts + its own measured cost. Byte-verifies every
 * fingerprint hit (a collision costs a memcmp, never splits a duplicate set).
 * Caller must ensure n is even and > 0. */
read_memo_result read_memo_prepass(const mem_opt_t *opt, const bseq1_t *seqs,
                                   int n, read_memo_state *st);

/* Feed one invocation's pre-pass result + that invocation's measured align-phase
 * cost (WORKER10 delta, ns) into the net-cycles controller, which predicts the
 * net benefit of memoizing (avoided align work - probe/copy overhead) and latches
 * ON/OFF once a z-test clears, with a work-based periodic re-probe. In Phase 1
 * this only drives the STATS/verbose reporting; nothing acts on the latch yet. */
void read_memo_controller_observe(const read_memo_result &r, uint64_t align_ns);

/* Current controller decision for the memoize path (Phase 2 will gate on this):
 * READMEMO_ON while auto is latched-on or mode==on; READMEMO_OFF otherwise. */
int read_memo_active(void);

/* Test-only: reset the process-global controller + stats to their initial
 * (unlatched, MEASURING) state so a unit test does not depend on prior state or
 * case ordering. Not used in production. */
void read_memo_reset_for_testing(void);

#endif /* READ_MEMO_H */
