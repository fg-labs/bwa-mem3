# Correctness Fixes

This page documents bugs present in upstream bwa-mem2 that bwa-mem3 fixes. Each
fix is isolated to a single PR so it can be reviewed independently and dropped
from `main` once upstream merges the equivalent patch.

## `@PG CL:` tab escaping (PR #54)

When a read-group string is passed via `-R '@RG\tID:x\tSM:y'`, the tab
characters in the argument were copied verbatim into the `@PG CL:` SAM header
field. The SAM specification uses tabs as field delimiters, so the resulting
header line appeared to have extra `ID:` and other tag fields embedded inside
`CL:`. Lenient parsers (samtools, htsjdk) tolerated the output; strict parsers
(noodles, some fgbio configurations) rejected the file as malformed.

The fix replaces each tab character with a space when building the `@PG CL:`
value in `src/main.cpp`. The `@RG` line itself is not modified, so the
read-group metadata is preserved correctly. A regression shell test
(`test/pg_cl_escape_test.sh`) asserts that the `@PG` line contains exactly
five tab-separated fields after the fix. Upstream issue reference:
[bwa-mem2#293](https://github.com/bwa-mem2/bwa-mem2/issues/293).

## SMEM buffer overflow on reads longer than 151 bp (PR #55)

bwa-mem2 hardcoded `READ_LEN 151` in `src/macro.h` to size the per-thread
`matchArray` SMEM buffer at compile time. The FMI walk wrote past this buffer
without bounds checking when reads exceeded 151 bp, causing memory corruption
that manifested as segfaults or silent wrong output on 300 bp MiSeq reads,
error-corrected long reads, and any run with a non-default `-k` that extended
seed length.

A second cap, `MAX_READ_LEN_FOR_LOCKSTEP 512`, guarded the lockstep driver's
per-slot stack arrays with a hard assert that aborted on anything longer.

The fix eliminates both compile-time caps. Every per-thread SMEM buffer is now
heap-allocated on the memory management context (`mmc`) and grown on demand
from each batch's observed `max_readlength`. The pre-walk grow in
`mem_collect_smem` sizes `matchArray[tid]` to `BATCH_MUL * BATCH_SIZE *
max_readlength`, and all array writes are bounds-checked with a structured
`smem_overflow_die` on overflow. Regression tests cover 300 bp, 1 kbp, and
3 kbp phiX reads; all three segfaulted before the fix and produce correct
NM:i:0 alignments after. Upstream references:
[bwa-mem2#210](https://github.com/bwa-mem2/bwa-mem2/issues/210) (issue),
[bwa-mem2#238](https://github.com/bwa-mem2/bwa-mem2/pull/238)
(closed unmerged upstream PR).

## `kswv` nrow==0 guard (PR #51)

When a SIMD batch contained only padding pairs (all `len1 == 0`), the DP loop
never executed and `nrow` was zero. The post-loop `rowMax + (i-1) * SIMD_WIDTH`
store still executed, walking SIMD_WIDTH bytes before the beginning of the
`rowMax` allocation. On glibc this produced a `free(): invalid pointer` abort;
on macOS libc it silently corrupted the heap.

The fix wraps the post-loop store in an `if (i > 0)` guard on all five SIMD
kswv kernels: NEON u8, NEON 16, AVX2 u8, AVX-512BW u8, and AVX-512BW 16. The
upstream patch [bwa-mem2#289](https://github.com/bwa-mem2/bwa-mem2/pull/289)
covered only the two AVX-512BW kernels; bwa-mem3 broadens it to the three
additional kernels carried in this fork. A dedicated regression test
(`test/kswv_nrow_zero_test.cpp`) builds all-padding batches and verifies each
kernel is clean under AddressSanitizer.

## kswv score2 plateau series (PRs #26, #27, #28, #29, #30, #31)

The batched mate-rescue Smith-Waterman path (`kswv`) contains a family of
related bugs across its SIMD kernels that inflated the suboptimal score
(`score2` / `XS`) and consequently deflated MAPQ relative to upstream
bwa-mem2.

**AVX-512BW dispatch guard (PR #26).** GCC with `-mavx512bw` automatically
defines `__AVX2__`, so the `#elif __AVX2__` branch in `src/kswv.h` and
`src/kswv.cpp` matched first on every AVX-512BW build. The 256-bit AVX2 kernel
produced only 32-lane results into 64-lane `score[]`/`te1[]`/`qe[]` arrays
sized for AVX-512BW; the upper 32 lanes held uninitialized values.
`mem_matesw_batch_post` read those bogus `te` values, `bwa_gen_cigar2`
returned NULL, and `mem_reg2aln` triggered an `a.cigar != NULL` assertion on
every AVX-512BW dispatch host (AWS c7a, c7i). The fix qualifies the `#elif
__AVX2__` guard with `!__AVX512BW__`, matching the existing pattern in
`bandedSWA.h`. Closes [issue #25](https://github.com/fg-labs/bwa-mem3/issues/25).

**AVX2 score2 plateau fix (PR #27 closed, PR #28 merged).** The AVX2 256-bit
`kswv` kernel added in PR #20 used a dense SIMD max over every `rowMax` row to
compute the suboptimal score. Scalar `ksw_u8` instead collapses consecutive
rows above `minsc` into a single `b[]` entry anchored at the max-score row,
then finds the best anchor outside the primary region. The dense max pulled in
tail rows from a plateau whose anchor sat inside the primary region, inflating
`XS` by 1–4 on a minority of reads and reducing MAPQ by 2–18 on those reads.
PR #27 (closed) temporarily disabled the AVX2 batched path. PR #28 fixes the
kernel itself by replacing the dense scan with a per-lane scalar emulation of
the `b[]` build-and-scan logic.

**NEON and AVX-512BW 8-bit port (PR #29).** The same dense-rowMax score2 scan
existed in `kswv_neon_u8` and `kswv_512_u8`. Confirmed on ARM: rebuilding
smoke-1M on darwin/arm64 pre-fix produced the identical four MAPQ regressions
as the AVX2 case. PR #29 ports the per-lane scalar `b[]`-emulation fix to both
kernels.

**AVX-512BW 16-bit port (PR #30).** `kswv_512_16` carried four bugs: the same
dense-rowMax plateau pattern, aggregate `maxl`/`minh` bounds instead of
per-lane bounds (a gap from PR #21), no `minsc` filter, and no `qe` mask. The
per-lane scalar emulation from PR #29 fixes all four naturally.

**NEON 16-bit rewrite (PR #31).** `kswv_neon_16` was effectively dead code
before this PR. Five interacting bugs produced 20,435 BAM diffs vs scalar
reference on smoke-1M `-A 2`: the score table reinterpreted `int16` xor
indices as `int8` lookups (inflating match scores by ~256 per cell), the table
was too small for the 16-bit SoA encoding, `rowMax` was never written, the
early-exit fired on row 0 for all pairs without a `KSW_XSTOP` target, and all
the fix-3 class bugs from PRs #28–#30 were missing. The PR rewrites the kernel
from scratch against `kswv_neon_u8`'s structure using 32-byte `int8` tables
indexed via `vqtbl2_s8`, per-lane freeze, `exit0` bitmap, and per-lane scalar
score2.

## `kseq2bseq1` zero-initialization (PR #22)

`bseq_read_orig` grows its sequence buffer with `realloc`, leaving tail entries
uninitialized. `kseq2bseq1` populated only `name`, `comment`, `seq`, `qual`,
and `l_seq` for each entry, leaving `sam`, `bams`, `n_bams`, and `cap_bams` at
whatever values `realloc` happened to return. PR #13 added an unconditional
`free(ret->seqs[i].bams)` in the output loop (`fastmap.cpp:571`), which turned
those garbage values into a crash — a `pointer being freed was not allocated`
abort under system malloc and a SIGSEGV under mimalloc — once input exceeded
the initial 256-sequence allocation. The crash was deterministic and
reproducible with `-t1`.

The fix is a single `memset(s, 0, sizeof(*s))` at the top of `kseq2bseq1`.

## Proper-pair flag from emitted alignment (PR #17)

In the `no_pairing` emission path of `mem_sam_pe` and `mem_sam_pe_batch_post`,
the proper-pair bit (`0x2`) was computed from `a[i].a[0].rb` regardless of
which alignment was actually emitted. When the primary's alignment score fell
below the reporting threshold `opt->T` but a non-primary ALT hit cleared it,
`mem_reg2aln` emitted `a[i].a[n_pri[i]]` while `mem_infer_dir` still read the
below-threshold primary. In that case the SAM flag did not reflect the
coordinates in the record.

`#17` stored the selected alignment index per mate in a `which[2]` array and
passed `a[i].a[which[i]].rb` to `mem_infer_dir`, so the proper-pair flag matched
the emitted record. **That derivation is now opt-in**: since
[#362](https://github.com/fg-labs/bwa-mem3/issues/362) the default reads `a[0]`,
as bwa and bwa-mem2 do, and `a[which]` is selected only under
`--proper-pair-from-emitted` — see the note below for why. The `a[0]` derivation
was present in the bwa-mem2 initial commit from 2019, inherited from bwa, which
still has it at 0.7.19 (`bwamem_pair.c:411`). Upstream reference: pre-existing,
no open upstream PR.

> **Opt-in since [#362](https://github.com/fg-labs/bwa-mem3/issues/362) — the
> default now matches both upstreams.** `0x2` is aligner-defined ("properly
> aligned according to the aligner"), so deriving it from `a[which]` is a
> defensible choice rather than a correction — and as a *default* it made
> bwa-mem3's output differ from bwa **and** bwa-mem2 simultaneously. The
> difference is confined to ALT-aware runs, and within those to the records
> where the emitted alignment is not the top-scoring one (`a[which] != a[0]`,
> which needs the read's top primary region to score below `-T` while its top
> ALT region clears it); a run without a `.alt` sidecar is byte-identical
> either way, excluding the `@PG` record, which records the command line.
> Measured at **3,013 of 10,134,006 records** on a 5 M-pair WGS
> slice (HG00096, GRCh38 with the standard `.alt`), on an AWS c6a.4xlarge host
> (AMD EPYC Milan, x86_64, AVX2 tier), all on decoy and ALT/HLA contigs, none
> on the primary assembly.
>
> The behavior is now reached via `--proper-pair-from-emitted` (a hard error
> with `--compat`, overridable via `--compat-allow-divergent`), and the default derives `0x2` from `a[0]` as bwa and
> bwa-mem2 do. This follows the same call [#256](https://github.com/fg-labs/bwa-mem3/pull/256)
> made on gap-open convention: where a fork-carried improvement moves output,
> the drop-in path takes upstream's answer and the improvement becomes opt-in.
>
> The divergence is unreachable without a `.alt` sidecar — `a[which] != a[0]`
> requires `n_pri < a.n`, i.e. the read has ALT hits, and `is_alt` is never set
> without one.

## 8-bit banded-SW envelope cap at `w = 124` (PR #422)

`bsw8_envelope_ok` — the gate that routes a seed-extension pair to the 8-bit
banded-SW kernel — admitted band widths up to `BSW8_MAX_W = 127`. The 8-bit kernel
encodes its band/position quantities as **signed-int8 diagonal offsets** spanning
`[-(w+1), w+3]` (the band-grow term `myband+1` and the tail-trim term `index+2`,
which reaches `w+3`). At `w >= 125` the positive edge `w+3` exceeds `+127` and wraps
negative: the band collapses on row 0, `head > tail`, and lanes die mid-alignment —
the pair then silently returns "no extension" (`score = h0`, `tle = qle = 0`,
`gscore = -1`) where the scalar and 16-bit kernels return the real extension. The
fix drops `BSW8_MAX_W` to **124**, so `w` in {125, 126, 127} route to the 16-bit
tier (exact, byte-identical to scalar) instead.

The overflow begins at `w >= 125`; `-w 127` is the characterized case — the one
user-supplied width that reaches the wrap on essentially every 8-bit-routed
extension (the wide-band retry doubling to 200/400 is already diverted to 16-bit).
`w = 125/126` also cross the boundary but corrupt only exotic pair shapes, so `-w 127`
is the width the tests reproduce. The default `opt->w = 100` never reaches the 8-bit
cap, so the drop-in path is **byte-identical** — verified by a records-only
(header-excluded) whole-aligner md5 on a 2×150 bp WGS workload, where `-w 100` output
is identical to `main`. A `getScores8` characterization test
(`test/unit/test_bandedswa_longread.cpp`, all six result fields vs the scalar
`scalarBandedSWA` oracle across the SIMD tiers exercised in CI) is byte-identical at
the `w = 124` bound and diverges at `w = 127`, locking the reason for the cap; at
`-w 127` this branch's records-only md5 matches an all-16-bit reference while `main`'s
8-bit tier differs.

## Vector banded-SW z-drop gate at `-d 0` (PR #424)

`scalarBandedSWA` gates its z-drop early-exit on `zdrop > 0`, so with z-drop
disabled (`-d 0`) it never truncates. The vector banded-SW kernels (8-bit and
16-bit, every SIMD tier) applied the z-drop test **unconditionally**: at
`zdrop == 0` they compared the running drop against a zero threshold and killed a
lane as soon as its score fell one point below the row max, truncating alignments
the scalar reference runs to completion. The two then disagreed on the query-end
fields (`gscore`/`gtle`), and once the vector stopped early, on the local fields
(`score`/`tle`/`qle`/`max_off`) too. The fix gates all six vector z-drop mask
updates (the three `ZSCORE16` macros and the three 8-bit wide-`epi32` paths) on
`zdrop > 0`, matching the scalar.

The effect is confined to **`-d 0` (z-drop disabled)**, where essentially all
extension pairs previously diverged between the vector tiers and the scalar
reference; disabled z-drop now preserves the alignments the old kernels wrongly
truncated. On the default `-d 100` the guarded branch is always taken, so the
kernels are **byte-identical** there and the drop-in path is unchanged. This was
verified by a controlled records-only (header-excluded) whole-aligner md5 A/B —
control (`main`) versus this branch on identical inputs, default flags, and
batching (`-K` held fixed for both runs) — on a 2×150 bp WGS workload (HG002,
hg38): the default-`-d 100` output matches the pre-fix baseline on both the ARM64
(NEON) and x86_64 (AVX2) tiers. At `-d 0` the fixed kernels' records-only md5
instead equals that same default-flags output — disabling z-drop preserves
exactly the alignments the default z-drop leaves untruncated on these reads —
while the old kernels diverge. A standalone kernel byte-identity gate also passes
at default parameters (extend and rescue). A regression test
(`test/unit/test_bandedswa_zdrop_gate.cpp`) compares `getScores8`/`getScores16` at
`zdrop == 0` to `scalarBandedSWA` on NEON, AVX2, and AVX-512BW: it fails on the old
kernels and passes on the fixed ones.

## 16-bit banded-SW per-lane band clamp (PR #423)

The three 16-bit banded-SW wrappers (`smithWatermanBatchWrapper16`, at 128/256/512
bits) computed the per-lane band-clamp reach with a **16-bit modular add**, reading
the sum back through a `uint16_t`. Because `qlen[l]` already holds `qlen * max_sc`,
the reach `qlen*max_sc + (end_bonus - o)` can go negative (a short query under a
large gap-open, e.g. `-O16`) or exceed 65535; the 16-bit add then wraps and the
clamp evaluates to a huge band that effectively disappears, so `getScores16` runs
the full band `w` where the scalar reference (`scalarBandedSWA`) clamps it as narrow
as 1. The 8-bit wrappers were already fixed for this; this ports the identical wide
`int` computation (mirroring the scalar "adjust `w` if it is too large" block) to
the three 16-bit wrappers.

The effect is scoped to **non-default scoring** (`-B9 -O16 -E3` and similar) on
short-query pairs, where the mis-sized band diverged from the scalar reference on
the query-end fields (`gscore`/`gtle`) and could change `score`, endpoints, and the
emitted alignment records. On bwa's default `-O6 -E1` at normal read lengths the
intermediate provably never wraps, so the wide-int path computes the identical
band and the drop-in path is **byte-identical by construction** — an
integer-arithmetic invariant that holds independent of host, compiler, and thread
count, so the proof (not any single benchmark) is the load-bearing evidence. A
controlled records-only (header-excluded) whole-aligner md5 A/B corroborates it:
control (`main`) versus this branch on a 2×150 bp WGS workload (HG002, hg38) at
default flags, with inputs, flags, thread count (`-t`), and batching (`-K`) all
held fixed across the two runs — so the batch composition is identical — on both
the ARM64 (NEON) and x86_64 (AVX2 and AVX-512BW) tiers. A `getScores16` band-clamp
regression test compares `score`, `tle`, `qle`, and `max_off` to the
`scalarBandedSWA` oracle unconditionally, and `gscore`/`gtle`
only when either side reports `gscore > 0` (a to-end alignment is observable),
across the wrap scoring sets on NEON, AVX2, and AVX-512BW: it fails on the old
kernel and passes on the fixed one.

## Symmetric ref/query length bounds guards on the SW wrappers (PR #467)

Every SIMD banded-SW wrapper (`getScores8`/`getScores16`, all widths) copies each
pair's reference window (`len1`) and query (`len2`) into the SoA lane buffers and
sizes the per-lane band from those lengths. `SeqPair::len1`/`len2` are signed
`int32_t`, and the wrappers previously did almost no length validation: five of the
six checked neither length, and the sixth (the AVX2 8-bit wrapper) checked only
`len2`'s **upper** bound — a non-fatal `stderr` warning plus a debug `assert` — while
no wrapper ever checked `len1` or a lower bound. So a **negative** length, an
out-of-range `len1` in any wrapper, or an over-limit `len2` in any of the other five
would index the SoA copy/padding buffers out of bounds (negative `k`) or wrap the
per-lane `uint8_t`/`uint16_t` band quantity. The guards are now symmetric: each wrapper
rejects both `len1` and `len2` outside `[0, MAX_SEQ_LEN8)` / `[0, MAX_SEQ_LEN16)`
via `err_fatal` **before any use of those lengths** — the SoA copy and padding
loops and the per-lane band-size calculation (`qlen[j] = len2 * max`, whose signed
multiply an unvalidated length could overflow) — and the fatal diagnostic states
the accepted range rather than only "exceeds".

This is an **invalid-input contract**, not a change to any aligned result: valid
alignment pairs always carry `0 <= len < MAX_SEQ_LEN`, so the guard is never taken
on the valid path and the emitted records are **byte-identical by construction** —
the change only converts previously-undefined behavior on malformed input into a
deterministic fatal error. Empirically, a record-level comparison of all result
fields against the pre-change baseline showed zero discordant pairs for the
extension and mate-rescue kernels on arm64 (Apple Silicon, NEON tier), over
synthetic pairs at query/reference lengths 100/150/250 (including ambiguous bases)
at a fixed batch composition; the x86 tiers (SSE4.1, AVX2, AVX-512BW) are built and
exercised on every CI run. That scopes the measured result to this workload, host,
and tier — the by-construction invariant above is what extends it to all valid
inputs.

---

## Banded-SW `qlen` band-clamp reach slot overflow (PR #468)

The batch banded-SW wrappers (`smithWatermanBatchWrapper{8,16}`, all widths)
pre-scale the per-lane band-clamp reach as `len2 * max_sc` and stash it in an SoA
slot the clamp loop reads back. That slot was the narrow per-lane element type —
`uint8_t` on the 8-bit tiers, `uint16_t` on the 16-bit tiers — so the scaled reach
wraps once `len2 * max_sc` exceeds 255 / 65535. `max_sc` is the largest
scoring-matrix entry, i.e. the match reward `A` (mismatch and ambiguous entries are
stored negative), so the wrap trigger is `len2 * A` crossing the slot ceiling and
does **not** depend on `-B`. A wrapped reach makes the vector band **narrower** than
`scalarBandedSWA` computes (which keeps the reach in native `int`), so the extension
can miss an off-diagonal optimum the scalar reference finds. The fix carries the
scaled reach in a wide `int32_t qlen_scaled[]` slot, matching the scalar band
formula exactly. This is the storage-slot sibling of #423, which fixed the reach
*sum* (a `uint16_t` modular add) but still read `len2 * max_sc` back from the narrow
slot.

The reachable divergence is the **16-bit tier only**. On the 8-bit tier the
production routing envelope (`bsw8_envelope_ok`) admits a pair only when
`h0 + len2 * A < 255 - max_step` with `len1 >= len2`, so `len2 * A < 255` for every
routed pair and the narrow `uint8_t` slot never overflows — the 8-bit widening is
**defensive** (no shipped-output change). The 16-bit slot overflows at
non-default `-A >= 3` with a long query segment (`len2 * A > 65535` — e.g.
`len2 ~ 21846 bp` at `-A3`; `MAX_SEQ_LEN16 = 32768` caps `len2 <= 32767`). Output can
diverge only when an off-diagonal optimum requires a band wider than the wrapped
reach; `-A2` cannot reach it (`len2 * 2 <= 65534`). On bwa's default
`-A1` and normal read lengths the product never wraps, so the wide slot computes the
identical band and the emitted **alignment records** are **byte-identical by
construction** for a fixed batch composition — an integer-arithmetic invariant
independent of host, compiler, and thread count; the SAM header block is outside this
claim.

A `getScores{8,16}` band-clamp regression test compares `score`, `tle`, `qle`, and
`max_off` — plus `gscore`/`gtle` under the kernel's query-end contract — to the
`scalarBandedSWA` oracle across the wrap scoring sets, including a case that
constructs the reachable 16-bit long-read overflow: it fails on the old narrow-slot
kernel and passes on the fixed one (18 subcases across 3 doctest cases, 36 assertions).
Its pairs are generated with deterministic, bounded query/reference lengths at a fixed
batch composition, and it runs on the SIMD tier its binary is compiled for — SSE4.1 and
AVX2 on CI's x86 rows, and NEON on the arm64 row and local Apple Silicon builds.
A full-field kernel byte-identity gate over the extension and mate-rescue result
fields — every field
per pair, 60,009 pairs per cell — additionally showed **zero** discordant pairs at
query/reference lengths 100/150/250 for default `-A1` (and a 10%-ambiguous-base
variant) and non-default `-A3` on the forced 16-bit tier, at a fixed batch
composition on arm64 (Apple Silicon M2 Max, NEON tier); the x86 tiers (SSE4.1, AVX2, AVX-512BW) are built and exercised
on every CI run. That scopes the measured result to this workload, host, and tier —
the by-construction invariant above is what extends it to all valid inputs.
bwa-mem2 inherits the same narrow-slot overflow, so on the affected reads this
diverges from bwa-mem2 while converging bwa-mem3's vector tiers onto its own scalar
reference.

---

## `--meth` SEQ restore clamps unsupported bytes to `N` (PR #463)

In `--meth` mode the original read sequence is restored into the emitted BAM `SEQ`
field along one of four branches, chosen by whether the pre-c2t original sequence
(`meth_orig_seq`, else the `YS:Z` comment) is available and by the read's orientation:

- **Forward, original sequence available** — the bytes are copied directly from the
  original sequence (only lowercase folded to uppercase); they do **not** pass through
  `nst_nt4_table`, so a `-` on this path is emitted verbatim as `-`.
- **Reverse, original sequence available**, plus both **fallback paths** (original
  sequence missing, restoring from `s->seq`, forward and reverse) — each ASCII base is
  mapped through `nst_nt4_table` and used to index a five-character string literal
  (`"ACGTN"` / `"TGCAN"`, indices 0–4).

`nst_nt4_table` maps `A/C/G/T/N` (and their lowercase forms) to `0..4` and every
unrecognized byte to `4` — with one exception: `-` maps to `5`. Index `5` selects the
literal's terminating NUL, so on the three lookup-based paths a `-` emitted a `'\0'`
into `SEQ`: an in-bounds read of the string terminator, but corrupt output rather than an
out-of-bounds access. The fix clamps any index `>= 4` to `4`, so on those paths a `-`
now emits `N`.

This is scoped to **`--meth` mode only**. On the lookup-based paths the clamp changes
output **only for a `-` byte** — the single input that maps to `5` — which now renders
as `N` instead of a corrupt NUL; every other byte already resolved to an index in
`0..4`, so its output is unchanged. That equivalence is **by construction** — the clamp
alters the result solely for the `bi == 5` case — not a measured-fixture generalization;
standard fixtures do not exercise the `-` input. The direct-copy forward path (where a
`-` stays `-`) and the default (non-`--meth`) output are unaffected.

## Extension staging copied a chain's windows once per seed (PR #524)

Before banded extension, `stage_seed_extension` copies each seed's extension targets into
per-thread staging buffers: the reversed reference prefix and the reference suffix of the
chain's window, and the reversed read prefix and the read suffix. That per-seed staging is
inherited from bwa-mem2. Every seed got its own copy, and the reference window spans the whole
chain, so a chain cost O(seeds x chain span) bytes (O(seeds x read length) on the read side).

A 150 bp read from a long homopolymer or short tandem repeat has an SMEM with tens of thousands
of hits; seeding keeps a `max_occ` (500) sample of them, and when their spacing is within the
band they all chain into one ~500-seed chain spanning the repeat. A batch of a few hundred such
reads grew the reference staging buffer past its `int32` offset range, and `mem` aborted with
`seqBufRef in stage_seed_extension cannot grow past the int32 offset range`. `--meth` hit it
most readily (it disables the ungapped fast path), but the default path failed the same way.

Every left target of a chain is a suffix of one reversed prefix and every right target a
suffix of one forward suffix, so each is now staged once per chain and each extension points
into the shared copy. Each extension reads exactly the bytes it read before (only their offset
in the staging buffer changes), so on any input that did not abort the emitted records are
**byte-identical by construction**. Inputs that aborted because a many-seed repeat chain copied
its window once per seed, such as the homopolymer reads above, now complete. The `int32` staging
limit itself is unchanged, so a batch whose once-per-chain windows still exceed it (as
kilobase-scale reads can) still aborts with the same error.

Measured scope: SAM minus `@PG` compared against `main` and identical on an AWS c8g.8xlarge
(Graviton4, arm64, NEON kernel tier, Amazon Linux 2023, clang 19.1.7) for hg38 (GRCh38 with
ALT/decoy) with 200,000 real WGS pairs (HG00096) at `-t 32` and `-t 8`; 20,000 simulated pairs on
the 1 Mb synthetic reference, paired (`-t 1`/`-t 8`, `-P`, `-5SP`, `-a`, `-M`, `--bam=0`,
`-P -a`) and single-end; and `--meth` on 3,000 bisulfite-converted pairs, paired with `-C` and
single-end. An earlier revision (reference-side sharing only) was also identical on an Apple M2
Max (arm64, NEON, clang 21.1.8) for 200,000 HG002 WGS pairs on hg38 and 20,000 simulated phiX
pairs. x86 was not compared against `main` for this change; CI's cross-arch check (this branch
on x86-64, AVX-512BW where the runner has it and otherwise AVX2, and on arm64 NEON; ~415K
simulated chr17 pairs) found identical alignment records, and
the staging code is common to every SIMD tier, running before the tier-specific kernels. Pinned by
`test/regression/repeat_chain_extension_window.sh`.

## 8-bit banded-SW score for a negative seed score depended on its SIMD group (PR #528)

The 8-bit banded-SW extension kernels (`smithWaterman128_8`, `smithWaterman256_8`,
`smithWaterman512_8`) keep each lane's best score in a wide per-lane side channel,
`best_abs`, seeded with the raw seed score `h0`, while the byte DP state is seeded with `h0`
clamped to 0. The per-row wide update ran `best_abs = max(best_abs, byte max)` over every
lane of the group whenever *any* lane's maximum advanced that row, on the assumption that
`best_abs` is never below the byte max. A negative `h0` breaks that assumption: a lane whose
cells never score above 0 was lifted from `h0` to 0 exactly when another lane in its group
advanced, and kept `h0` otherwise, so its `score` depended on which pairs shared its group.
The fix updates `best_abs` only on lanes whose maximum advanced that row, as `xrow` already
was, so a negative-`h0` lane now reports what the scalar kernel reports whatever it is
batched with. (On aarch64 outside Apple silicon, `smithWaterman128_8` runs a lean per-row
path, `BSW8_ROW_LEAN` in `src/bandedSWA.cpp`, that does not keep `best_abs` per row. The
byte maximum only rises, and only on the rows where the masked update above fires, so the
path rebuilds the same value at the result store: `max(h0, final byte maximum)` for a lane
whose maximum ever advanced, and `h0` otherwise.)

Negative seed scores occur only under **`--meth`**, where a seed found in converted space is
rescored in original space. There the extension-DP dedup (`--dedup`, default `auto`) latches
on or off from measured wall time and changes which pairs share a group, so identical
`--meth` runs could place a few MAPQ-0 multi-mapper pairs differently (`XA`/`HN` changing,
occasionally the pair moving to another contig). With the fix, `--meth` output no longer
depends on the dedup state, the thread count, or the SIMD tier (AVX2 matches AVX-512BW).
Versus the previous release, `--meth` records whose extension had a negative seed score can
change: on a 5 M-pair EM-seq panel sample, 707 of 10.37 M records under
`--meth-scoring neutral` and 950 under `genomic`; default (collapsed) `--meth` output was
unchanged on that sample. For `h0 >= 0` the masked update is the identity, so **non-`--meth`
output is byte-identical by construction**, corroborated by identical md5s on 5 M-pair WGS
and WES slices with the dedup at `auto` and forced on. Pinned by
`test/bandedswa_negative_h0_test.cpp`, which compares every result field against the scalar
oracle and against the pair scored alone.

## Extension retry ladder stopped early on the ungapped band bound (PR #540)

Before banded SW, the ungapped fast path can prove a *tight band* for an extension: no
alignment reaching a diagonal offset at or beyond it can score above the ungapped walk. The
retry ladder (`w = 100`, then `200`, ...) stopped as soon as the band it had run covered that
bound. Otherwise it stops on the same per-rung test as bwa and bwa-mem2, which have no ungapped
fast path: the score is unchanged from the previous rung, or the alignment's maximum diagonal
offset stayed under three quarters of the band. (bwa-mem3's ladder then allowed up to four
rungs where upstream allows two; see the next entry.) The bound keeps the score right, but when
it stopped the ladder at `w = 100` where the full ladder goes on to `w = 200`, the region
recorded a band width of 100
instead of 200. That width bounds the contained-seed test that decides whether a later seed in
the same region is extended at all, so a seed that bwa skips was extended, and its region's
score could become `XS`. On a synthetic 306 bp read the primary alignment was unchanged
(`143M18I145M`, `AS:i:193`), but `XS` was 91 where bwa 0.7.19 reports 37, as does bwa-mem3
built with the ungapped fast path compiled out (a local build; no option selects it). The bound
also does not pin the kernel's early-termination control flow at the wider rung, so the recorded
extent fields were not provably the full ladder's either.

The ladders that claim byte-identity (the default certified adaptive band, and the full-width
ladder with it disabled, `--no-band-cert`) no longer stop on the bound. At every rung from
`w = 100` up they apply exactly the full-width stop test, so they are **byte-identical to the
full-width ladder by construction**. The bound still feeds the certified narrow probe, which
finalizes on its own certificate. The non-byte-identical `--adaptive-band` ladder (also selected
by `--fast` unless `--no-adaptive-band`) keeps stopping on it on the tiers where it narrows the
starting band; on the 8-bit tier, which starts at `w = 100` and handles the short extensions
that carry a tight band, it now uses the exact ladder too, so `--adaptive-band` stays a no-op on
short reads. The removed early stop almost never fired on its own: on 2 x 1 M 150 bp WGS reads
it never stopped a pair that the ladder's own test would not have stopped, and on 1.86 M
~250 bp single-end SBX HG002 reads it did so for 12 of 24.5 M ladder accepts, none of which
changed a record. Because those pairs now run one more rung, the cost is negligible.

Measured scope: SAM minus `@PG` identical to `main` and to `main` with the ungapped fast path
compiled out, on an Apple M3 Ultra (arm64, NEON, Apple clang 21.0.0), hg38 (GRCh38 with
ALT/decoy), `-t 16`, for 1 M HG002 WGS pairs (paired, single-end, `-A 2 -B 8`, and
`-O 8 -E 2`), 5 M HG00096 WGS pairs, and the 1.86 M SBX reads. Wall time on the 5 M pairs was
32.07-32.12 s against 32.11-32.20 s for `main` (3 interleaved runs each), within run-to-run
noise. `--adaptive-band` output was identical to the default on the 1 M HG002 pairs and to
`main --adaptive-band` on the SBX reads (8.92-8.94 s against 8.95-9.05 s wall). Pinned by `test/tight_band_xs_test.sh`, which checks bwa 0.7.19's record for that read
by default, under `--compat=bwa-mem2`, `--no-band-cert`, and `--adaptive-band`.

---

## Extension retry ladder ran four rungs where upstream runs two (PR #543)

bwa and bwa-mem2 both define `MAX_BAND_TRY 2`: an extension is scored at the band width `-w`
(100), retried once at twice that if the stop test fails there (the score changed and the
alignment's maximum diagonal offset reached three quarters of the band), and the 200-wide
result is then kept whatever the test says. bwa-mem3 had carried `MAX_BAND_TRY 4` since #58,
where it served an initial-band experiment (start at 8 and double) that the same change
reverted, leaving the ladder at `[100, 200, 400, 800]`. A pair that failed the stop test at
200 therefore ran on. That is not a no-op: the wider rungs can find a different alignment when
a gap lies beyond diagonal offset 200, and they always record a wider band width on the
region, which bounds the contained-seed purge, the CIGAR band in `mem_reg2aln` and the band
in `mem_patch_reg`.

When it can happen is bounded exactly. Failing the test at 200 needs a record-setting cell at
diagonal offset 150 or more, which costs a gap of at least `o_min + e_min * 150` and pays back
at most `a` per aligned column, so the extension needs more than 156 query bases at the
default `-w` and scoring (`a*L > o_min + e_min*3w/2` in general). A 150 bp read at the default
`-w` never has such an extension; longer reads can. A synthetic 900 bp read with insertions at
diagonal offsets 80, 150 and 250 reports `209M80I141M70I125M275S`, `AS:i:248` in bwa 0.7.19 and
reported `209M80I141M70I130M100I170M`, `AS:i:267` under the four-rung ladder, in every
exact-ladder mode (default, `--compat=bwa-mem2`, `--no-band-cert`).

The exact ladders now have upstream's two rungs: the default certified ladder runs its narrow
probe and then `[100, 200]`, and the full-width ladder (`--no-band-cert`) runs `[100, 200]`.
`--adaptive-band`'s narrowing tiers keep their four rungs, since that ladder is not
byte-identical by design and needs them to climb from its narrow start; its 8-bit tier follows
the exact ladder, so it stays a no-op on short reads. The ladder is shared code, identical on
every SIMD tier.

Measured scope: on an Apple M3 Ultra (arm64, NEON, Apple clang 21.0.0), hg38 (GRCh38 with
ALT/decoy), `-t 16`, SAM minus `@PG` was identical to the previous ladder on 1 M HG002 WGS
150 bp pairs (paired and single-end), 1.86 M ~250 bp single-end SBX HG002 reads (also under
`--adaptive-band`), and the 1 M pairs at `-A 2 -B 8` and `-O 8 -E 2`: no extension there
reached a third rung. On 100 HG002 HiFi reads under `-x pacbio` the former ladder accepted
1,472 extensions at 400 and 153 at 800. The records of 4 reads changed (2 under
`--compat=bwa-mem`), and under `--compat=bwa-mem` all 109 records now match bwa 0.7.19 exactly
(4 lines differed before). Pinned by `test/ladder_rungs_test.sh`, which
checks bwa 0.7.19's record for the 900 bp read, its reverse complement (a left extension) and
a read shaped for `--adaptive-band`'s narrowing rungs, in the exact-ladder modes and, for the
last, under `--adaptive-band`. `test/unit/test_extension_ladder.cpp` checks the length bound on
a model of the ladder over a parameter grid against both `scalarBandedSWA` and upstream's
`ksw_extend2`, and that past it a third rung changes the committed alignment.

## Ungapped fast path: record tie-break and z-drop (PR #544)

Before banded SW, an extension whose diagonal holds at most `x_threshold` mismatches (one at
default scoring) is committed from an ungapped walk instead, since no gapped alignment can tie
it. The walk must then reproduce every field the extension kernel would have reported, and two
did not match:

- The walk moved its record to the *later* of two tied positions, while the kernels'
  cross-row record test is strict (`m > max`, the earlier row keeps a tie; the rightmost
  tie-break the comment cited is the kernels' within-row rule, which a diagonal walk never
  exercises). A diagonal that returns to its record without exceeding it, a mismatch followed
  by exactly `b/a` matches, therefore reported a longer `qle`. At the default `-L 5` the
  difference is invisible: such an extension ends at its record, so the query-end score equals
  the local score and the clip decision takes the query-end branch on both sides. It shows
  under `-L 0`, where bwa reports `5S145M` and bwa-mem3 reported `150M` with `NM:i:1`, and
  under gap costs that admit a second mismatch after the tie (`-O 20`: `7S143M` against
  `2S148M`).
- The kernels z-drop when the score falls more than `-d` below the running maximum, and one
  mismatch drops it by `b`. With `-d` below `b * x_threshold` the kernels stop at the mismatch
  where the walk went on: under `-d 3` bwa reports `11S139M`, `AS:i:139` and bwa-mem3 reported
  `150M`, `AS:i:145`.

The walk now records a new maximum only on a strict increase, and `ungapped_analyze` falls
back to the kernel when `b * mismatches` exceeds a positive `-d`. The fast path is also off
when `-w` is below 2 (0 or 1): at those widths the ladder cannot accept its first rung, so the
band it records would differ from the fast path's. That also drops the tight-band bound there.
The full argument that a HIT commits exactly the ladder's fields, including the case where the
diagonal falls to zero, is written above `ungapped_analyze` in `src/ungapped_ext.h`. It is
checked in `test/unit/test_ungapped_fastpath.cpp` against both `scalarBandedSWA` and upstream's
`ksw_extend2` over a `-A/-B/-O/-E/-L/-d/-w` grid (asymmetric gap costs and queries up to the
fast path's 512-base limit included) with random and adversarial pairs, and end to end in
`test/ungapped_hit_parity_test.sh`.

At default parameters the output is unchanged: the z-drop guard is inert at `-d 100`, the
`-w` gate cannot fire at `-w 100`, and the tie-break is invisible at `-L 5` with one admissible
mismatch. Measured scope: on an Apple M3 Ultra (arm64, NEON, Apple clang 21.0.0), hg38 (GRCh38
with ALT/decoy), `-t 16`, SAM minus `@PG` was identical to the previous build on 1 M HG002 WGS
150 bp pairs (paired and single-end), 1.86 M ~250 bp single-end SBX HG002 reads, and the 1 M
pairs at `-A 2 -B 8` and `-O 8 -E 2`. On the 1 M pairs, bwa-mem3 under `--compat=bwa-mem` now
matches bwa 0.7.19 line for line at `-L 0`, `-O 20` and `-d 3`, where 8,158, 485 and 210,601
reads changed (29,824, 1,588 and 836,494 lines differed from bwa before), as well as at
defaults. Wall time on 5 M HG00096 pairs was 31.90-32.05 s against 32.07-32.24 s before
(3 interleaved runs each).

## `--meth` mate rescued off a reverse-strand anchor took the other strand (PR #553)

Mate rescue searches for the missing mate in a window next to its anchor. The window lies on one
strand half of the doubled reference: normally the anchor's, but `bns_clamp_window` picks the
half from the window's midpoint, so near the end of the last contig it can be the other one.
Under `--meth` the rescue scores the mate with the matrix that frees its read's conversion,
`(mate_meth_ot ^ is_rev) & 1`, where `is_rev` says only whether the mate lies on the opposite
strand *from the window's half*. That choice is relative to the window's half, and it is the
right matrix for that window. `mem_matesw_batch_post` stored the same value as the rescued
region's `meth_hypothesis`, which the rest of the aligner reads as the genome-strand hypothesis
(OT or OB) and flips itself for regions in the reverse half: chain extension, and `mem_reg2aln`
when it regenerates the CIGAR. The two agree only when the window is on the forward half.

So a mate rescued off a reverse-strand anchor came out with the other strand: an OT read 1
mapped forward was written `XG:Z:GA`, an OB read 2 mapped forward `XG:Z:CT`, and its `XM:Z` was
read off the wrong base. The same happened to a mate rescued off a forward anchor whose window
was clamped onto the reverse half at the end of the last contig. `mem_reg2aln` regenerated the
mate's CIGAR and `NM` under the other conversion, so each C→T (or G→A) counted in `NM` as a
mismatch (`MD`, literal under `--meth`, lists every difference either way and only follows the
CIGAR);
the global regeneration could also place an indel differently, which can move `POS` by a few
bases. The rescued span (`qb`/`qe`, and so any soft clip) and the rescue score were computed with
the right matrix, so the mate's locus, the pairing, `AS` and `XS` did not change, nor did `MAPQ`
except where `--chimera-qc` caps it on the regenerated CIGAR (below); the
anchor's record changes only where it describes the mate (`MC:Z`, and `PNEXT`/`TLEN` when `POS`
moved). A rescued region kept as an alternative hit carries the same hypothesis, so its CIGAR and
edit distance in `XA:Z` can change too. With the opt-in QC filters the flags can change as well:
`--set-as-failed f|r` tests the mate's strand, and `--chimera-qc` its longest match in the
regenerated CIGAR (also capping `MAPQ` at 1), and either result is propagated to every record of
the pair (`0x200` set, `0x2` cleared).

The fix takes the stored hypothesis from the rescued region's own strand, as the extension and
`mem_reg2aln` do: the mate's read-number chemistry XOR whether the region lies in the reverse half.

This is scoped to **`--meth` paired-end mate rescue** (every chemistry and `--meth-scoring` mode,
TAPS included) **whose window lies on the reverse half**. Other records are byte-identical by
construction, apart from a rescued alternative hit in `XA:Z`: the stored value is unchanged
whenever the window is on the forward half, and the hypothesis is `-1` off `--meth`, so mates
placed from their own seeds, single-end runs, `-S` runs and the default (non-`--meth`) path are
unaffected. Pinned by `test/regression/meth_rescue_strand.sh`, which builds one rescued mate per
fragment strand and anchor strand plus one at the end of the contig whose window is clamped onto
the reverse half. Under collapsed and genomic scoring each must keep its position, strand and
`XG`, align full length with `NM` equal to its planted errors, and have exactly its converted
cytosines called in `XM`; under `--meth=taps` its strand and `XG`, with no conversion in `NM`.
The reverse-anchor and contig-end mates fail it without the fix.

## Upstream port divergences: all chains dropped by the weight filter (PR #489) and the SA sentinel offset (PR #469)

bwa-mem2 advertises output identical to bwa. Two records are known where its port is not
faithful, and on both the default path follows bwa; `--compat=bwa-mem2` reproduces bwa-mem2 as
shipped, because parity with that release is its contract (see
[Equivalence → Byte-identical output](equivalence.md#byte-identical-output---compat)).

- **All chains dropped by the weight filter** ([#310](https://github.com/fg-labs/bwa-mem3/issues/310),
  default fixed in [#489](https://github.com/fg-labs/bwa-mem3/pull/489)). With `-W` above zero or
  an `-x pacbio`/`pbref`/`ont2d` preset, when the weight filter drops every chain for a read, bwa
  reports zero survivors and the read goes out unmapped. bwa-mem2's seqid-range bookkeeping
  rebuilds a count of one for the emptied array and extends the chain the filter just rejected —
  whose seeds it has already freed, so what it extends depends on the allocator. Unreachable at
  default settings, since `min_chain_weight` defaults to 0. Pinned by
  `test/regression/chain_flt_empty_compat.sh`.
- **Suffix-array sentinel offset** ([#469](https://github.com/fg-labs/bwa-mem3/pull/469)).
  bwa-mem2's software-pipelined SA lookup drops the accumulated LF-walk offset when the walk
  reaches the sentinel row, placing a hit within the first few bases of the concatenated reference
  up to a few bases too far left; bwa's `bwt_sa` counts the walk. Latent on hg38 (which opens with
  a long `N` run), reachable on small or custom references. Pinned by
  `test/sa_lookup_sentinel_parity_test.cpp` and `test/regression/compat_sa_sentinel.sh`.

## Changes catalog

| Item | bwa-mem3 PR | Upstream PR/issue | Status |
|------|-------------|-------------------|--------|
| `@PG CL:` tab escape | [#54](https://github.com/fg-labs/bwa-mem3/pull/54) | [bwa-mem2#293](https://github.com/bwa-mem2/bwa-mem2/issues/293) | fork-only (open upstream issue) |
| SMEM buffer overflow on >151 bp reads | [#55](https://github.com/fg-labs/bwa-mem3/pull/55) | [bwa-mem2#238](https://github.com/bwa-mem2/bwa-mem2/pull/238), [bwa-mem2#210](https://github.com/bwa-mem2/bwa-mem2/issues/210) | fork-only (upstream PR closed unmerged) |
| kswv nrow==0 guard | [#51](https://github.com/fg-labs/bwa-mem3/pull/51) | [bwa-mem2#289](https://github.com/bwa-mem2/bwa-mem2/pull/289) | fork-only (upstream PR open) |
| AVX-512BW dispatch guard | [#26](https://github.com/fg-labs/bwa-mem3/pull/26) | — | fork-only |
| AVX2 score2 plateau disable (superseded) | [#27](https://github.com/fg-labs/bwa-mem3/issues/27) | — | closed (superseded by #28) |
| AVX2 score2 plateau fix | [#28](https://github.com/fg-labs/bwa-mem3/pull/28) | — | fork-only |
| NEON + AVX-512BW 8-bit score2 fix | [#29](https://github.com/fg-labs/bwa-mem3/pull/29) | — | fork-only |
| AVX-512BW 16-bit score2 fix | [#30](https://github.com/fg-labs/bwa-mem3/pull/30) | — | fork-only |
| NEON 16-bit kernel rewrite | [#31](https://github.com/fg-labs/bwa-mem3/pull/31) | — | fork-only |
| 8-bit banded-SW envelope cap at `w = 124` | [#422](https://github.com/fg-labs/bwa-mem3/pull/422) | — | fork-only (only affects `-w >= 125`; default `-w 100` byte-identical) |
| Vector banded-SW z-drop gate at `-d 0` | [#424](https://github.com/fg-labs/bwa-mem3/pull/424) | — | fork-only (z-drop-disabled only; default `-d 100` byte-identical) |
| 16-bit banded-SW per-lane band clamp (wide arithmetic) | [#423](https://github.com/fg-labs/bwa-mem3/pull/423) | — | fork-only (non-default scoring only; default path unchanged — see the correctness note above) |
| `--meth` SEQ restore clamps unsupported bytes to `N` | [#463](https://github.com/fg-labs/bwa-mem3/pull/463) | — | fork-only (`--meth` only; changes a `-` byte to `N` on the `nst_nt4_table` lookup paths only — the direct-copy forward path, all other bytes, and the default path are unchanged) |
| Symmetric ref/query length bounds guards on the SW wrappers | [#467](https://github.com/fg-labs/bwa-mem3/pull/467) | — | fork-only (invalid-input contract; valid-input records byte-identical by construction — see the PR #467 correctness note above for the measured scope) |
| Banded-SW `qlen` band-clamp reach slot overflow (`int32_t` reach) | [#468](https://github.com/fg-labs/bwa-mem3/pull/468) | — | fork-only (8-bit widening defensive; 16-bit divergence only at non-default `-A >= 3` long-query, `len2 * A > 65535` — default alignment records unchanged for fixed batch composition, see the correctness note above) |
| Extension staging copied a chain's windows once per seed | [#524](https://github.com/fg-labs/bwa-mem3/pull/524) | — | fork-only (inputs that aborted because a repeat chain staged its window once per seed now complete; the `int32` staging limit itself is unchanged; records of inputs that already ran byte-identical by construction — see the PR #524 correctness note above for the measured scope) |
| 8-bit banded-SW score for a negative seed score depended on its SIMD group | [#528](https://github.com/fg-labs/bwa-mem3/pull/528) | — | fork-only (`--meth` only; records whose extension had a negative seed score can change and are now independent of batching; non-`--meth` output byte-identical by construction — see the PR #528 correctness note above) |
| Extension retry ladder stopped early on the ungapped band bound | [#540](https://github.com/fg-labs/bwa-mem3/pull/540) | — | fork-only (a record can change only where the early stop diverged from the full ladder, by recording a narrower band or different extension end fields; none changed on the WGS or SBX sets measured — see the correctness note above) |
| Extension retry ladder ran four rungs where upstream runs two | [#543](https://github.com/fg-labs/bwa-mem3/pull/543) | — | fork-only (only an extension of more than 156 query bases at the default `-w 100` and scoring can reach a third rung, so under those defaults 150 bp reads are unchanged by construction; a narrower `-w` lowers the bound (`a*L > o_min + e_min*3w/2`, more than 36 query bases at `-w 20`); longer reads can change CIGAR, `AS` and the recorded band where the stop test failed at `2w` — see the correctness note above) |
| Ungapped fast path: record tie-break and z-drop guard | [#544](https://github.com/fg-labs/bwa-mem3/pull/544) | — | fork-only (non-default parameters only: `-L 0`, gap costs that admit a second mismatch after a tie, `-d` below `b * x_threshold`, or `-w` below 2; the default `-L 5 -d 100` is byte-identical — see the correctness note above) |
| `--meth` mate rescued off a reverse-strand anchor took the other strand | [#553](https://github.com/fg-labs/bwa-mem3/pull/553) | — | fork-only (`--meth` paired-end only, every chemistry; a mate rescued in a window on the reverse half gets its own `XG`/`XM` and a CIGAR/`NM` that no longer count its conversions (`MD` follows the CIGAR), where an indel can move and shift `POS`; its locus, `AS` and `XS` are unchanged, as is `MAPQ` unless `--chimera-qc` caps it on the regenerated CIGAR, its anchor changes only in `MC`/`PNEXT`/`TLEN`, a rescued alt hit in `XA:Z` can change, and `--set-as-failed`/`--chimera-qc` can flip the pair's `0x200`/`0x2`; every other record is byte-identical by construction — see the correctness note above) |
| kseq2bseq1 zero-initialization | [#22](https://github.com/fg-labs/bwa-mem3/pull/22) | — | fork-only |
| Proper-pair flag from emitted alignment | [#17](https://github.com/fg-labs/bwa-mem3/pull/17) | — | fork-only, **opt-in** (`--proper-pair-from-emitted`; default matches both upstreams, [#362](https://github.com/fg-labs/bwa-mem3/issues/362)) |
| All chains dropped by the weight filter: default follows bwa | [#489](https://github.com/fg-labs/bwa-mem3/pull/489) | — ([#310](https://github.com/fg-labs/bwa-mem3/issues/310)) | fork-only; `--compat=bwa-mem2` reproduces bwa-mem2 |
| Suffix-array sentinel offset in the prefetch lookup | [#469](https://github.com/fg-labs/bwa-mem3/pull/469) | — | fork-only; `--compat=bwa-mem2` reproduces bwa-mem2 |

---

**See also:**
[Performance improvements](performance.md) ·
[Architecture support](arch-support.md) ·
[Fork changes vs. upstream](../reference/pr-catalog.md) ·
[Developer Guide → Regression test framework](../developer-guide/regression-tests.md) ·
[Performance → SIMD dispatch matrix](../performance/simd-dispatch.md)
