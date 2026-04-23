# bwa-mem2 Performance Porting Plan

Rank-ordered port of performance wins from `ksw2`, `minimap2`, `mm2-fast`, and `sassy` into `bwa-mem2`, one PR at a time. All ports are output-preserving (byte-identical sorted SAM) unless explicitly flag-gated in **Phase K**.

Plan file maintained as a working tracker; updated after each phase (or each PR where relevant). Source audits are in the conversation transcript that produced this plan.

## Status legend

- `PENDING`  — not started
- `ACTIVE`   — in progress, branch open
- `PR`       — PR submitted, awaiting review/merge
- `MERGED`   — merged to `fg-labs/fg-main`
- `DROPPED`  — investigated and rejected (reason in ledger)
- `BLOCKED`  — waiting on prerequisite or external factor

## Execution protocol

Between each PR:
1. **Golden diff test**: SAM md5 must be byte-identical against `fg-main` baseline (except Phase K with flag ON).
2. **Benchmark harness**: full run on the fixed harness input; record wall-clock, max RSS, per-stage `__rdtsc` counters (once PR 1 lands).
3. **Ledger update**: commit delta (`delta_total_s`, `delta_stage_X`, `delta_rss`) to the **Execution ledger** below.
4. **Plan update**: if the PR changed our understanding, update **Plan modifications** with the delta to the plan.
5. **Branch hygiene**: one PR per worktree off `fg-labs/fg-main`, branch name `perf/<short-slug>`, no tracking until first push.

A PR with <0.5% measured gain merges anyway if cheap — compounding is the whole point of the 2% bar — but the expectation miss is noted in the PR body and the ledger.

## Worktree / branch conventions

Per `/Users/nhomer/work/git/bwa-mem2/CLAUDE.md`:
- New worktree per PR: `/Users/nhomer/work/git/bwa-mem2/perf-<slug>/`
- Branch: `perf/<slug>`, off `fg-labs/fg-main`, upstream unset until first push
- First push: `git push -u fg-labs HEAD`

---

## Phase A — Measurement foundation

| Status | # | PR | Effort | Est. gain | Source | Notes |
|---|---|---|---|---|---|---|
| `PR` | 0 | Benchmark harness + golden SAM diff test | Small | 0% (enables all) | — | Branch `perf/bench-harness`. Harness at `bench/` (`run.sh`, `compare.sh`, `normalize_sam.sh`, `README.md`, `config.env.example`, `.gitignore`). Baseline logged below. |
| `PENDING` | 1 | `DISABLE_OUTPUT` compile flag + `__rdtsc` counters around seeding / chaining / DP / format | Trivial | 0% (telemetry) | mm2-fast `map.c:213-249` | Per-stage attribution so later PRs are correctly attributed. |
| `PENDING` | 2 | `(size_t)` cast in large-allocation arithmetic in `bandedSWA.cpp` + `FMI_search.cpp` | Trivial | — | ksw2 commit `ebea754` | Correctness; prevents latent overflow on long reads. |
| `PENDING` | 3 | `-flto=thin` + `-fno-semantic-interposition` + codegen-units=1 equivalent in Makefile | Trivial | 3–6% | sassy `[profile.dist]` | One Makefile PR. PGO deferred to **Phase I** #29. |
| `PENDING` | 4 | Tune default `actual_chunk_size` / batch size (probe 2×, 4×) | Trivial | 1–5% | minimap2 `options.c:188` | Single-line default change + bench sweep. |

## Phase B — Trivial microops

| Status | # | PR | Effort | Est. gain | Source | Notes |
|---|---|---|---|---|---|---|
| `PENDING` | 5 | `LIKELY`/`UNLIKELY` on skewed branches in `bandedSWA.cpp`, `kswv.cpp`, `FMI_search.cpp` | Trivial | 0.5–2% | minimap2 `ksw2_ll_sse.c` | `ksw.cpp` already has them; extend. |
| `PENDING` | 6 | `__builtin_prefetch(&seeds[j-PFD])` in chain-predecessor scans | Trivial | 1–4% | idea | bandedSWA already prefetches; chaining does not. |
| `PENDING` | 7 | Replace small `malloc`/`calloc` top-K arrays with stack VLAs | Trivial | 0.5–2% | sassy `seed.c:61` | Audit `bwamem.cpp` inner paths. |
| `PENDING` | 8 | Fast `mg_log2` float-bit trick for `logf` in chain/score loops | Trivial | 1–3% | minimap2 `mmpriv.h:139-147` | Guard `len >= 2`; validate integer-cast preservation. |
| `PENDING` | 9 | `objdump` check: does compiler already fold striped-profile `qp + target[i]*slen`? Code the hand-table only if not. | Trivial | 0–2% | sassy `tqueries.rs:16-22` | Skip-if-no-win PR. |

## Phase C — Layout / allocator hygiene

| Status | # | PR | Effort | Est. gain | Source | Notes |
|---|---|---|---|---|---|---|
| `PENDING` | 10 | `perf c2c` audit; 64-byte-align + pad per-thread scratch in `kswv.cpp`, `bandedSWA.cpp` to kill false-sharing | Small | 1–5% on many-core | sassy `LaneState` | Diagnosis first, fix follows. |
| `PENDING` | 11 | Per-worker reusable `kswq`/query-profile scratch on `worker_t` | Small | 1–3% | minimap2 `ksw2_ll_sse.c:37-83` | Eliminate per-extension allocator churn. |
| `PENDING` | 12 | Single-block 64-byte-aligned DP scratch arena (fold multi `_mm_malloc` into one slab) | Small | 1–4% | ksw2 `ksw2_extz2_sse.c:84-96` | — |
| `PENDING` | 13 | `KALLOC_POOL_INIT`-style freelist for kbtree nodes / `mem_chain_t` / cigar objects | Trivial | 0–3% | minimap2 `kalloc.h:64-93` | Macro-based, drop-in. |

## Phase D — bandedSWA inner-loop microops (one kernel, four PRs, ordered so each is locally verifiable)

| Status | # | PR | Effort | Est. gain | Source | Notes |
|---|---|---|---|---|---|---|
| `PENDING` | 14 | Pre-reverse query/target once outside DP at extension call sites | Small | 1–3% | minimap2 `ksw2_extz2_sse.c:106-107` | Removes per-row bit-slicing. |
| `PENDING` | 15 | Masked tail-store to fold scalar epilogue on partial tile (`_mm512_mask_storeu_epi8` + AVX2 analog) | Small | 1–3% | mm2-fast `ksw2_extd2_avx.c:104-105,290-299` | Do before #16/#17 — tidies boundaries. |
| `PENDING` | 16 | XOR+shuffle `pmat` score-LUT (`_mm512_shuffle_epi8` of 16-entry table indexed by `q^t`, `N→8` remap) | Small | 2–5% | mm2-fast `ksw2_extd2_avx.c:201-209` | Local swap; trivial to validate. |
| `PENDING` | 17 | Profile-generation loop fission (row score vector in its own pre-pass, DP core after) | Medium | 2–5% (biggest on AVX-512) | ksw2 `ksw2_extz2_sse.c:125-144` | After #16 so fissioned pre-pass uses shuffle-LUT. |

## Phase E — IO / threading / formatting

| Status | # | PR | Effort | Est. gain | Source | Notes |
|---|---|---|---|---|---|---|
| `PENDING` | 18 | Inline integer + fixed-3-decimal SAM-tag formatter (replace `vsnprintf` for `pa:f:%.3f` etc.) | Small | 1–3% | minimap2 `format.c:27-66` | — |
| `PENDING` | 19 | SIMD nibble-LUT base encoding / N-count (`_mm256_shuffle_epi8` + `PACKED_NIBBLES`) | Small | <2% unless hot | sassy `iupac.rs:71-131` | Gate on Phase A profiler data. |
| `PENDING` | 20 | Split `kt_pipeline` into genuine 3-step reader / mapper / writer; `--2-io-threads` flag | Small | 2–6% wall-clock | minimap2 `kthread.c:78-159` | Order preservation already enforced in `fastmap.cpp`. |

## Phase F — Thread-local allocator

| Status | # | PR | Effort | Est. gain | Source | Notes |
|---|---|---|---|---|---|---|
| `PENDING` | 21 | Thread-local `kalloc` arena on `worker_t` + `--cap-kalloc` reset; plumb `void *km` through chain/seed/ksw hot paths | Medium | 2–8% | minimap2 `kalloc.c/.h` | Coexists with mimalloc. Minimal-risk variant first: per-read scope only. |

## Phase G — Seeding / FM-index cache work

| Status | # | PR | Effort | Est. gain | Source | Notes |
|---|---|---|---|---|---|---|
| `PENDING` | 22 | Radix sort for fully-unique 64-bit chain/region keys (`(score<<32)|id`) | Small | ~1% | minimap2 `ksort.h` | Audit tie-stability before each swap. |
| `PENDING` | 23 | k-way heap merge of per-interval-sorted SMEM hit streams (avoid global resort) | Small-medium | 1–3% | minimap2 `collect_seed_hits_heap` | — |
| `PENDING` | 24 | **Batched FM-index Occ/SA lookup with software prefetch** (mm2-fast pattern *without* RMI) | Medium | 2–5% | mm2-fast `seed.c:74-120` | Prior work check: `perf-smem-lockstep` worktree and locked `perf/smem-batch-prefetch` agent worktree may already prototype this. Inspect before starting. |

## Phase H — Extension algorithmic (biggest scalar wins; strong diff-test required)

| Status | # | PR | Effort | Est. gain | Source | Notes |
|---|---|---|---|---|---|---|
| `PENDING` | 25 | Two-pass approximate-then-exact DP in `bandedSWA::getScores8/16` (`KSW_EZ_APPROX_MAX`/`APPROX_DROP` adaptation) | Medium | 3–7% overall | minimap2 `align.c:834-848` | Approximate first pass skips 32-bit H[] tracking; exact only on z-drop/success. |
| `PENDING` | 26 | Ungapped fast-path: when `qe-qs == re-rs`, ungapped scorer before banded DP; fall back if ungapped can't provably beat gapped | Medium | 5–15% on clean reads, 3–7% overall | minimap2 `align.c:773-789` | Biggest single algorithmic win. Validate bit-identical on 10M+ reads. |

## Phase I — Structural refactors & deployment consolidation

| Status | # | PR | Effort | Est. gain | Source | Notes |
|---|---|---|---|---|---|---|
| `PENDING` | 27 | SoA recast of `mem_seed_t` / `mem_chain_t` (with 16-elt head padding); update all scans in `mem_sort_dedup_patch`, `mem_chain_flt`, `mem_chain_weight` | Medium-large | 2–8% | mm2-fast `lchain.c:301-336` | Gate on Phase A showing LLC-bound. Touches many files. |
| `PENDING` | 28 | In-process ISA dispatch via `__attribute__((target_clones(...)))` + IFUNC resolver; retire `runsimd.cpp` `execv` model → single binary | Medium | 5–30% on AVX-512 hosts + eliminates startup reindex cost | sassy | Late; don't juggle build system during inner-loop PRs. |
| `PENDING` | 29 | PGO training + wiring in Makefile (training profile on Phase A harness) | Small + infra | 3–6% on top of LTO | sassy | After #28 so single-binary PGO is meaningful. |

## Phase J — Learned index

| Status | # | PR | Effort | Est. gain | Source | Notes |
|---|---|---|---|---|---|---|
| `PENDING` | 30 | LISA / RMI learned-index for k-mer → SA-range lookup (trained via `learned-systems-rmi` crate, verified by byte-compare) | Large | 5–15% on top of #24 | mm2-fast `ext/TAL/lisa_hash.h` | Pull submodule first. Adds Rust build step to index pipeline. |

## Phase K — Flag-gated, output-changing

| Status | # | PR | Effort | Est. gain | Source | Notes |
|---|---|---|---|---|---|---|
| `PENDING` | 31 | `--fast-rng`: Wang-hash deterministic tie-break RNG (replaces `drand48_r`) | Small | 0.5–2% | minimap2 `map.c:246-248` | Changes output. Default off. |
| `PENDING` | 32 | `--chain-max-skip` / `--chain-max-iter` heuristic bails | Small | 1–3% | minimap2 `lchain.c:168-207` | Changes output on repetitive reads. Default off. |
| `PENDING` | 33 | `--seed-freq-filter`: 2-stage seed selection with top-K heap | Medium | 2–10% on repetitive reads | minimap2 `seed.c:5-132` | Changes chain set. Default off. |

---

## Aggregate expectation

Output-identical (PRs 0–30): **≈ 30–50% wall-clock** on typical short-read workloads, with the bulk from #24, #25, #26, #28, #30.

---

## Execution ledger

Format: `YYYY-MM-DD | PR # | status | wall Δ (harness) | RSS Δ | per-stage Δ | notes`

| Date | PR | Status | Wall Δ | RSS Δ | Notes |
|---|---|---|---|---|---|
| 2026-04-23 | — | plan | — | — | Plan drafted. |
| 2026-04-23 | 0 | `PR` | n/a (baseline) | n/a (baseline) | Harness + baseline on host `fg-nils-00005` (arm64, macOS) against `fg-main` at `61813ef`: 1M pairs hg38, 1-thread golden 28.82s, 4-thread wall median 15.90s (range 14.35–16.23 across 3 trials = 13% shared-laptop noise), golden md5 `e44035cc671c28f42e3ca84f10cf2e5a`, RSS ~16 GB. Representative ranking runs should re-measure on AVX-512 EC2 (`i-09b7c8ae622ba21e6`) with more trials + cold-cache discipline. |

---

## Plan modifications

Log changes to the plan itself (reordering, scope shifts, dropped items). One entry per modification.

| Date | Change | Reason |
|---|---|---|
| 2026-04-23 | Initial plan drafted. | — |
| 2026-04-23 | Plan file moved into `perf/bench-harness` branch so it ships to `fg-main` with PR 0 and is visible from every future branch. | Survives across worktrees. |
| 2026-04-23 | Noted that ranking-quality measurement requires an AVX-512 EC2 host; local arm64 macOS harness is for fast smoke iteration only. | Local host shows 13% wall-clock variance in 3 trials — too coarse for <5% perf deltas. |
