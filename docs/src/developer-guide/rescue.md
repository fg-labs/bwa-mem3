# Mate rescue internals

Mate rescue is the batched Smith-Waterman of a read's mate against the reference window its insert-size distribution allows, run when a pair did not place concordantly. It is one of the costliest stages of paired-end alignment, and a stack of shortcuts now runs inside it. Every one of them is exact: it may skip work, never change what the aligner writes (Gate: `rescue_prune_identity.sh`, and the rest under [How exactness is enforced](#how-exactness-is-enforced)). This page is the map; the sub-pages carry the arguments.

- [Exact window pruning](rescue-pruning.md): the 5-mer bound that proves a job cannot reach the rescue threshold, or narrows its window to a hull.
- [Banded rescue DP](rescue-banding.md): the DP of a narrowed job run only inside the diagonal bands that can hold an alignment, for both rescue passes.
- [The 11-op rescue cell](rescue-kswv.md): the kswv cell rebuilt around the score before the in-row gap, on every SIMD tier.

The prose here explains mechanisms and names the tests that enforce them. It deliberately carries no timings: those go stale with every compiler and host, so they live in the pull requests that introduced each change ([#533](https://github.com/fg-labs/bwa-mem3/pull/533), [#535](https://github.com/fg-labs/bwa-mem3/pull/535), [#536](https://github.com/fg-labs/bwa-mem3/pull/536) and their successors) and on the [benchmarks](../performance/benchmarks.md) page.

## The pipeline

A batch of pairs goes through three functions in `src/bwamem_pair.cpp`:

1. `mem_matesw_batch_pre` builds one rescue job per (anchor, orientation): the reference window and the oriented mate, staged into shared sequence buffers as a `SeqPair`. This is where the pruning filter (`rescue_prune_window`) runs. A job proven unable to reach the threshold is not enqueued at all; its slot records `MATESW_GAR_PROVEN_FAIL`. A narrowed job is staged as its hull only, and the hull's row offset is recorded per tid (`matesw_narrow_slot`) for the post step. The band plan of a narrowed job (`RescueBandBatch::plan`) is made here too, from the filter's view, and bound to the pair's `regid` with `RescueBandBatch::commit`.
2. `mem_sam_pe_batch` runs the kswv kernels. `mem_sam_pe_batch_run` scores pass 0 (score, end positions, suboptimal score) with `getScores8` / `getScores16`, then the banded parents' pass 0 (`RescueBandBatch::run_pass0`), then pass 1 (start positions) on the jobs that reached the threshold, banded where `RescueBandBatch::take_pass1` accepts the job and through kswv otherwise.
3. `mem_matesw_batch_post` reads each job's result back, applies the hull offset, and adds a rescued region to the mate's region list, which it then deduplicates.

All per-batch rescue state lives in `mem_cache`, one slot per tid (`rescue_narrow_off`, `rescue_band`), so nothing hides in thread-local statics across batches. The filter's scratch and its query cache are per-thread by design: they depend only on the job and the oriented mate.

## Which shortcut runs where

| Shortcut | Runs when | Turned off by |
|---|---|---|
| Exact pruning | a SIMD filter, so aarch64 or an x86 AVX2 / AVX-512BW build (`rescue_prune_on`), a scoring the lemma holds for (`rescue_prune_params::from`), `--meth` only on aarch64 with EM-seq chemistry, on x86 only where the SIMD filter runs and not at the AVX-512BW tier at `min_seed_len * a` of 25 or more (`-k 25` at `-A 1`; `rescue_prune_runs`), 8-bit job, no `--rescue-kmer` (`rescue_prune_applies`) | `BWA3_RESCUE_PRUNE=0` |
| Banded pass 0 | a NEON or AVX2 band kernel (`rescue_band_enabled`), a pruned job (at any scoring pruning admits, no `--meth`) whose band plan beats the hull (`RescueBandBatch::plan`) | `BWA3_RESCUE_BAND=0`, `BWA3_RESCUE_PRUNE=0` |
| Banded pass 1 | a NEON or AVX2 band kernel (`rescue_band_enabled`), every 8-bit job whose band is cheaper than kswv, narrowed or not (`RescueBandBatch::take_pass1`): at any scoring the band kernels take and any seed length (`rescue_band_runs`, `rb_scoring`), and under `--meth`, with the OT / OB group's matrix, where `--meth` does not prune (`rescue_band_meth_on`) | `BWA3_RESCUE_BAND_P1=0`, `BWA3_RESCUE_BAND=0`, `BWA3_RESCUE_PRUNE=0` |
| 11-op kswv cell | every SIMD kswv body, when the gap costs admit it (`fscan_scoring_ok`) | `BWA3_RESCUE_FSCAN=0` |
| Dedup skip and one-pass insert | every architecture and scoring | `BWA3_RESCUE_DEDUP_SKIP=0` |

Every threshold and cost gate between those switches is a knob too, so an A/B of any single decision is an environment variable away. The knobs, their defaults and the files that read them:

{{#include ../../_generated/rescue/knobs.md}}

The table is generated from the list in `src/rescue_env.h` (`make docs-rescue-knobs`); `rescue_docs_lint.sh` fails when the committed copy is stale, when src/ reads a knob the list lacks, or when a listed default disagrees with the reader's.

## The post-rescue dedup

`mem_matesw_batch_post` used to sort and deduplicate the mate's region list after every rescue orientation, including the ones that added nothing, and the counter that triggers it never resets, so a read with several anchors deduplicated an unchanged list again and again. The caller now carries a `mem_rescue_dedup_state_t` per mate across those calls. A dedup that finds every reference end distinct and drops nothing reports a fixed point (`mem_dedup_only_fixpoint`); repeating it on the same list cannot change a byte, so it is skipped. One new region added to a fixed-point list goes through `mem_dedup_only_insert1`, which places it in one linear pass when its end is distinct, it forms no redundant pair and no (score, rb, qb) tie exists, and runs the full dedup otherwise. The argument is in the comment above `mem_dedup_only_insert1` in `src/bwamem.cpp`; the result is identical to the full dedup on every call. Gates: `the mate-rescue dedup's fixed-point skip and single-insert path are byte-identical to the full dedup`, `rescue_prune_identity.sh`.

## How exactness is enforced

Each shortcut has an argument in its source header and at least one gate that would fail if the argument broke. The gates:

| Gate | What it compares |
|---|---|
| `rescue_prune_identity.sh` | Whole-aligner SAM at `-t 1` and `-t 4` on a generated paired-end fixture, the defaults against a reference with every shortcut and alternative kernel form off (the script's REF_ENV list: `BWA3_RESCUE_PRUNE=0 BWA3_RESCUE_BAND=0 BWA3_RESCUE_DEDUP_SKIP=0 BWA3_RESCUE_REPEAT=0 BWA3_RESCUE_FSCAN=0 BWA3_RESCUE_USQADD=0 BWA3_RESCUE_ROWPAIR=0 BWA3_RESCUE_LAZYQE=0`), on every CI row. Its fixture includes 300 bp mates, so the 16-bit kernels run too, and its stats checks make it fail rather than pass vacuously where a shortcut stops engaging or a leg's switch stops taking effect. Gate: `rescue_prune_identity.sh` |
| The same script under each forced x86 tier | The kswv rescue kernels differ per tier, so the canonical row reruns the reference and default legs under `BWAMEM3_FORCE_TIER` for every tier the runner has. Gate: `Rescue identity under each forced x86 tier (generated PE fixture)` |
| `rescue_prune_eq` fuzz | The SIMD filter's decisions, view, component list and repeat memo against the scalar filter: NEON on the ARM64 rows, the x86 port on the x86 rows (built at AVX2 and, where the runner has it, also run at AVX-512BW). Gate: `SIMD rescue filter == scalar filter (rescue_prune_eq, generated jobs)` |
| `rescue_band_harness` eq | Every banded output field against kswv on the full window, under each band knob: the NEON kernels on the ARM64 rows, the AVX2 kernels on the x86 rows. Gate: `Banded rescue == kswv (rescue_band_harness, generated jobs)` |
| Unit tests | The filter against `ksw_align2` and the scalar filter, the kswv cells against the scalar reference per tier, and the dedup shortcuts against the full dedup. Gates: `rescue prune: B1 and B2 decisions reproduce every consumed ksw_align2 field`, `kswv u8 rescue: BWA3_RESCUE_FSCAN off == on in every u8 body, and FSCAN matches scalar`, `the mate-rescue dedup's fixed-point skip and single-insert path are byte-identical to the full dedup` |

`rescue_docs_lint.sh` holds these pages to the code: every backticked name here must exist in the source, every exactness claim must name a gate that CI runs, and the rescue sources and these pages must name each other.

## Changing rescue code

- Keep the escape hatch. A new shortcut gets an off switch in `src/rescue_env.h`'s knob list and a row in the user-facing table in `docs/src/whats-different/performance.md`, and `rescue_prune_identity.sh`'s reference list (REF_ENV) turns it off.
- Measure non-vacuity, not only equality: an identity test that stopped exercising the shortcut still passes. The `BWA3_RESCUE_PRUNE_STATS` counters exist so the gates can require that each path ran.
- Ask whether the change invalidates these pages; the lint catches renamed symbols and missing gates, not an argument that no longer holds.
