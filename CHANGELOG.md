# Changelog

## [0.6.0](https://github.com/fg-labs/bwa-mem3/compare/v0.5.0...v0.6.0) (2026-07-16)


### Features

* **version:** report whether mimalloc is the active allocator ([#217](https://github.com/fg-labs/bwa-mem3/issues/217)) ([2cd9fe9](https://github.com/fg-labs/bwa-mem3/commit/2cd9fe9c27185af98d2f6a0b59c04431a50eea0a))


### Performance

* **bam:** cut per-record work in the --bam writer (mem_aln_to_bam) ([#212](https://github.com/fg-labs/bwa-mem3/issues/212)) ([837f40e](https://github.com/fg-labs/bwa-mem3/commit/837f40e1e974217ee1c45f7a4ceea6b3b6c00821))
* **bsw:** drop dead per-row H1/H2 setup stores in the SW batch wrappers ([#211](https://github.com/fg-labs/bwa-mem3/issues/211)) ([039b09f](https://github.com/fg-labs/bwa-mem3/commit/039b09f8ecd11e724b098a1729d42703c2ce058d))
* **dedup:** drop provably-dead exact-duplicate passes in mem_sort_dedup_patch ([#205](https://github.com/fg-labs/bwa-mem3/issues/205)) ([c544258](https://github.com/fg-labs/bwa-mem3/commit/c5442582e44da347d8a81e9b58a70a18ed0467a4))
* **fmi:** arm64 lockstep for third-pass reseeding (bwtSeedStrategy) ([#215](https://github.com/fg-labs/bwa-mem3/issues/215)) ([c5f69ec](https://github.com/fg-labs/bwa-mem3/commit/c5f69ec0aaa8d2fae290228c7b24ee774f740d62))
* **kswv:** remove write-only Hmax scratch buffer from batched SW kernels ([#214](https://github.com/fg-labs/bwa-mem3/issues/214)) ([b14ef49](https://github.com/fg-labs/bwa-mem3/commit/b14ef4993ee9a60813b0412b8a5d5d40aa19b66b))
* **mem:** add unique-mapper fast paths to mem_chain_flt and mem_mark_primary_se ([#209](https://github.com/fg-labs/bwa-mem3/issues/209)) ([1700aab](https://github.com/fg-labs/bwa-mem3/commit/1700aabc263a8776dc12b66acbb566bda67a37ba))
* **mem:** drop redundant per-call scratch allocations in mem_gen_alt ([#216](https://github.com/fg-labs/bwa-mem3/issues/216)) ([3418525](https://github.com/fg-labs/bwa-mem3/commit/3418525750b15df76dd0d5c5d3634dab4c8f61d4))
* **mem:** pool per-read scratch in get_sa_entries_prefetch and mem_reg2aln ([#208](https://github.com/fg-labs/bwa-mem3/issues/208)) ([5747a13](https://github.com/fg-labs/bwa-mem3/commit/5747a131cde9c22490c5ace470b0903d395ad130))
* **mem:** remove dead/instrumentation work from the extension hot path ([#206](https://github.com/fg-labs/bwa-mem3/issues/206)) ([439be97](https://github.com/fg-labs/bwa-mem3/commit/439be9709844a0557ad01ec7cdf3b542e9acf66d))
* **seed:** replace sortSMEMs qsort with a counting sort by rid ([#207](https://github.com/fg-labs/bwa-mem3/issues/207)) ([4356b5c](https://github.com/fg-labs/bwa-mem3/commit/4356b5ce50348a2421b2b8e2517d12f642cfb8e8))
* **seed:** vectorize backwardExt occ-counting on arm64 (NEON) + hoist invariants ([#210](https://github.com/fg-labs/bwa-mem3/issues/210)) ([64a55a2](https://github.com/fg-labs/bwa-mem3/commit/64a55a25c18350a66a93438d64f3d93a77259a8a))


### Refactoring

* **mem:** static cleanup — remove dead code, silence int64 format warnings ([#213](https://github.com/fg-labs/bwa-mem3/issues/213)) ([8c57b7e](https://github.com/fg-labs/bwa-mem3/commit/8c57b7ef6b0ff45d3b6d1b17031c9f0c46879199))


### Documentation

* retitle LICENSE for BWA-MEM3 and add Fulcrum copyright ([#203](https://github.com/fg-labs/bwa-mem3/issues/203)) ([089a956](https://github.com/fg-labs/bwa-mem3/commit/089a956375b27372d812a757cc87663526126fc8))

## [0.5.0](https://github.com/fg-labs/bwa-mem3/compare/v0.4.0...v0.5.0) (2026-07-04)


### Features

* add opt-in --seed-order seed reordering (default off, byte-identical) ([#186](https://github.com/fg-labs/bwa-mem3/issues/186)) ([04749a1](https://github.com/fg-labs/bwa-mem3/commit/04749a16fb2e2a4f0bd8a318fc6bdd6439e2689f))
* add opt-in --smem-dedup (dedup identical SMEMs before chaining) ([#187](https://github.com/fg-labs/bwa-mem3/issues/187)) ([1384972](https://github.com/fg-labs/bwa-mem3/commit/1384972417677943b68b8292181017223ad27876))
* **mem:** add --adaptive-band (chain-geometry adaptive banding) for long reads ([#194](https://github.com/fg-labs/bwa-mem3/issues/194)) ([4fe92a6](https://github.com/fg-labs/bwa-mem3/commit/4fe92a6ea55ebe8617cabed79902bceeabecbd2c))
* **mem:** add --extend-mate-concordant; fix --fast --meth placement regression ([#195](https://github.com/fg-labs/bwa-mem3/issues/195)) ([c9ffef1](https://github.com/fg-labs/bwa-mem3/commit/c9ffef1551d13ef22b104a75c5be785857e01efa))
* **mem:** add --fast speed preset ([#189](https://github.com/fg-labs/bwa-mem3/issues/189)) ([a946af8](https://github.com/fg-labs/bwa-mem3/commit/a946af8fcd0309faed952662123f88d152b3b585))
* **mem:** add --max-extend-chains and bundle it into --fast ([#193](https://github.com/fg-labs/bwa-mem3/issues/193)) ([e39b3d4](https://github.com/fg-labs/bwa-mem3/commit/e39b3d497e3e10535b2d9eabb9bb7fa48d1c1379))
* **mem:** add --skip-contained-ext and enable it under --fast ([#192](https://github.com/fg-labs/bwa-mem3/issues/192)) ([2d2b2b4](https://github.com/fg-labs/bwa-mem3/commit/2d2b2b437b70cdc5c445ea83da7cb158af654cd1))


### Bug Fixes

* **bandedSWA:** 8-bit SW drops query-end gscore/gtle on zero-score-row exit ([#198](https://github.com/fg-labs/bwa-mem3/issues/198)) ([611e21b](https://github.com/fg-labs/bwa-mem3/commit/611e21b924d09ccbd7e2e3d18092a90dc798e4b2))
* **bandedSWA:** getScores{8,16} must not scribble padding past numPairs ([#199](https://github.com/fg-labs/bwa-mem3/issues/199)) ([9aae808](https://github.com/fg-labs/bwa-mem3/commit/9aae808c0a721bf88dfc4c486c8bf8b785cdd508))


### Performance

* **bandedSWA:** gate the getScores overshoot guard to sub-slice callers ([#201](https://github.com/fg-labs/bwa-mem3/issues/201)) ([162e909](https://github.com/fg-labs/bwa-mem3/commit/162e90923488fea4634c00a990156794940b9797))


### Refactoring

* **kswv:** drop duplicate F warm-up prefetch in kswv512_16 ([#191](https://github.com/fg-labs/bwa-mem3/issues/191)) ([6e4cf2b](https://github.com/fg-labs/bwa-mem3/commit/6e4cf2b34a527309e8dc675e941a9fd548feadf1))


### Documentation

* **changelog:** backfill the 0.4.0 breaking-change notice ([#183](https://github.com/fg-labs/bwa-mem3/issues/183)) ([e1c381a](https://github.com/fg-labs/bwa-mem3/commit/e1c381ad4b27e7e3bcf26dbbfeea130cff6e2af0))
* **changelog:** render the live changelog, not the frozen NEWS.md ([#184](https://github.com/fg-labs/bwa-mem3/issues/184)) ([12a46d8](https://github.com/fg-labs/bwa-mem3/commit/12a46d8317fd3dabb031f835348ebdf126e76601))
* **contributing:** document breaking-change commit footers ([#181](https://github.com/fg-labs/bwa-mem3/issues/181)) ([4ebd122](https://github.com/fg-labs/bwa-mem3/commit/4ebd1220ffda0a0da544bfe0a9bab94840d73a8c))
* **release:** describe the release-please flow, not manual tagging ([#182](https://github.com/fg-labs/bwa-mem3/issues/182)) ([fc0c4b6](https://github.com/fg-labs/bwa-mem3/commit/fc0c4b65f356e6571c9f33f563a13d9fc6bd3f3e))

## [0.4.0](https://github.com/fg-labs/bwa-mem3/compare/v0.3.0...v0.4.0) (2026-06-27)


### ⚠ BREAKING CHANGES

* **index:** `bwa-mem3 index` no longer writes the unpacked `.0123` reference file by default ([#177](https://github.com/fg-labs/bwa-mem3/pull/177)). `bwa-mem3 mem` now reconstructs reference bases from the packed `.pac` on demand ("pac-fetch") and ignores any `.0123` present, so the file is redundant for bwa-mem3 itself. External tools that read `.0123` directly — most notably sharing a single index with **bwa-mem2** — will break, since the expected file is now absent. To restore the old on-disk layout, re-run indexing with the new opt-in flag `index --emit-unpacked-ref`. Alignment output is byte-for-byte identical; the change is purely to the index artifact set.


### Features

* **mem:** --min-ext-len opt-in filter to skip extension of short seeds ([#169](https://github.com/fg-labs/bwa-mem3/issues/169)) ([13db252](https://github.com/fg-labs/bwa-mem3/commit/13db252ad4c7c39eff18bd3674430aac51570346))
* **meth:** native bisulfite (BS-seq) alignment via --meth (D3) ([#174](https://github.com/fg-labs/bwa-mem3/issues/174)) ([a0296b1](https://github.com/fg-labs/bwa-mem3/commit/a0296b1766f17aab85b720fef1cf7f3984735105))


### Performance

* **mem:** pac-fetch the reference from .pac instead of loading/building .0123 ([#177](https://github.com/fg-labs/bwa-mem3/issues/177)) ([9c4bbf2](https://github.com/fg-labs/bwa-mem3/commit/9c4bbf217d3cbf02143a2fbd9fb75329126d9010))
* **meth:** batched (SIMD) asymmetric mate rescue (closes [#173](https://github.com/fg-labs/bwa-mem3/issues/173)) ([#175](https://github.com/fg-labs/bwa-mem3/issues/175)) ([f146a18](https://github.com/fg-labs/bwa-mem3/commit/f146a18a32232a6ea05b040b10917a464a2d264c))


### Documentation

* **book:** document the libdeflate build prerequisite (incl. AL2023) ([#172](https://github.com/fg-labs/bwa-mem3/issues/172)) ([c2b6ec7](https://github.com/fg-labs/bwa-mem3/commit/c2b6ec78bc8fadb1aa63f99aee4f90671f840e7d))
* **book:** recommend -y 0 (drop 3rd-round seeding) as an opt-in speed knob ([#171](https://github.com/fg-labs/bwa-mem3/issues/171)) ([ca9ac1f](https://github.com/fg-labs/bwa-mem3/commit/ca9ac1f2289c200708193c3f6ebbe094f230d79d))
* collapse the mdbook sidebar into nested, foldable sections ([#180](https://github.com/fg-labs/bwa-mem3/issues/180)) ([c6ac47b](https://github.com/fg-labs/bwa-mem3/commit/c6ac47b9ebc3146144377f4a56e689b3fdf11946))
* deep mdbook cleanup — dedup, consolidate, and tighten ([#179](https://github.com/fg-labs/bwa-mem3/issues/179)) ([54d6d11](https://github.com/fg-labs/bwa-mem3/commit/54d6d1161c374d7eb7fcef0eaac7c1156fe0bfb4))
* **meth:** disclose collapsed-mode placement drift vs bwameth.py ([#178](https://github.com/fg-labs/bwa-mem3/issues/178)) ([96f29e2](https://github.com/fg-labs/bwa-mem3/commit/96f29e257e9acd9072caa246b1589c8943babc26))
* **settings-profiles:** note repeat-aggregating downstream caveat for -m 10 ([#168](https://github.com/fg-labs/bwa-mem3/issues/168)) ([38fd1ec](https://github.com/fg-labs/bwa-mem3/commit/38fd1ec9717e29cace9eac33bd643a7da001da4e))

## [0.3.0](https://github.com/fg-labs/bwa-mem3/compare/v0.2.2...v0.3.0) (2026-06-21)


### Features

* **bsw:** make the 8-bit h0-prefix seed unsigned [0,255] ([#151](https://github.com/fg-labs/bwa-mem3/issues/151)) ([9f51c5f](https://github.com/fg-labs/bwa-mem3/commit/9f51c5f521612a7c3e1c0dae9d71db38be0121b8))
* **bsw:** recover the 8-bit banded Smith–Waterman path for reads ≥128 bp ([#140](https://github.com/fg-labs/bwa-mem3/issues/140)) ([155a916](https://github.com/fg-labs/bwa-mem3/commit/155a91632de602663377197ff1168d200a24f344))
* **kswv:** AVX2 16-bit mate-rescue kernel (kswv256_16) ([#162](https://github.com/fg-labs/bwa-mem3/issues/162)) ([9107b82](https://github.com/fg-labs/bwa-mem3/commit/9107b825189da87822bccac152e4a98daf00aac1))
* **meth:** carry original-reference @SQ M5/UR and @CO/@PG into --meth headers ([#139](https://github.com/fg-labs/bwa-mem3/issues/139)) ([e94ad8b](https://github.com/fg-labs/bwa-mem3/commit/e94ad8b847b6637b95aea6dff3e07f560e22a532))
* **prof:** off-by-default --profile stage-timing instrumentation ([#152](https://github.com/fg-labs/bwa-mem3/issues/152)) ([83cf7ab](https://github.com/fg-labs/bwa-mem3/commit/83cf7ab9d8e4589354e16654726fbdfaa1b4b63c))
* **reader:** content-detecting FASTQ reader fast path (libdeflate BGZF) ([#128](https://github.com/fg-labs/bwa-mem3/issues/128)) ([cdd71bf](https://github.com/fg-labs/bwa-mem3/commit/cdd71bf4590cf2e1a5d908d6e852f780c93c3a8a))


### Bug Fixes

* **bsw:** bound getScores8/16 prefetch reads to the padding contract ([#150](https://github.com/fg-labs/bwa-mem3/issues/150)) ([87ed5d4](https://github.com/fg-labs/bwa-mem3/commit/87ed5d4a78b6b78789bca45be5649373576ee668))
* **fmi:** widen mem_lim to int64 and guard SA-entry allocations ([#156](https://github.com/fg-labs/bwa-mem3/issues/156)) ([2d18c1e](https://github.com/fg-labs/bwa-mem3/commit/2d18c1ed00de5e4a60fa31d15fec04c7b74877a0))
* **kthread:** drive kt_for with a persistent worker pool ([#154](https://github.com/fg-labs/bwa-mem3/issues/154)) ([26b24e7](https://github.com/fg-labs/bwa-mem3/commit/26b24e7bac9fd7f73921b902ba98a628776087b6))
* **meth:** emit -R read group as @RG header in --meth mode ([#137](https://github.com/fg-labs/bwa-mem3/issues/137)) ([ccd1fc5](https://github.com/fg-labs/bwa-mem3/commit/ccd1fc56d26bfd54ab1f9604845332353afd75a6))
* **seeding:** widen SMEM read positions from int16_t to int32_t ([#142](https://github.com/fg-labs/bwa-mem3/issues/142)) ([037c418](https://github.com/fg-labs/bwa-mem3/commit/037c418b3b3c49deeb41bd57e37eacc7895e6fca))
* **test:** make meth layer-2 FAIL diagnostics reachable under set -e ([#133](https://github.com/fg-labs/bwa-mem3/issues/133)) ([d2d6688](https://github.com/fg-labs/bwa-mem3/commit/d2d66880bce736900b81e7610cbae1d1556ba7ea))


### Performance

* **bsw:** AVX2 SIMD tuning for the Smith-Waterman kernels ([#161](https://github.com/fg-labs/bwa-mem3/issues/161)) ([458b216](https://github.com/fg-labs/bwa-mem3/commit/458b21653cb6f84808032360a85a78f291064a4e))
* **bsw:** NEON SIMD tuning for the Smith-Waterman kernels ([#160](https://github.com/fg-labs/bwa-mem3/issues/160)) ([d971ff0](https://github.com/fg-labs/bwa-mem3/commit/d971ff0694d8a90cf836dfff873534c9d3a160a4))
* **bsw:** prefetch next batch's ref/query in the AVX2 8-bit wrapper ([#163](https://github.com/fg-labs/bwa-mem3/issues/163)) ([e6082a0](https://github.com/fg-labs/bwa-mem3/commit/e6082a05fc6af07da26dae90210ea651b601dd91))
* **bsw:** short-circuit the inert per-row re-baseline scan ([#147](https://github.com/fg-labs/bwa-mem3/issues/147)) ([8e284e0](https://github.com/fg-labs/bwa-mem3/commit/8e284e0cf274e8c77c78457e372c1ab92de69f3f))
* **bsw:** vectorize the per-row epilogue side-channel loop ([#149](https://github.com/fg-labs/bwa-mem3/issues/149)) ([403aeb7](https://github.com/fg-labs/bwa-mem3/commit/403aeb76ad8ec2664736998c018a7e033af0c070))
* **fmi:** size SA-entry staging buffers to the exact write count ([#157](https://github.com/fg-labs/bwa-mem3/issues/157)) ([aa0fe33](https://github.com/fg-labs/bwa-mem3/commit/aa0fe3335c58badf97893fae32f67ff02f69217a))
* **read:** vendored zlib-ng inflate + chunk cap + 3rd pipeline worker ([#153](https://github.com/fg-labs/bwa-mem3/issues/153)) ([5cf89e3](https://github.com/fg-labs/bwa-mem3/commit/5cf89e3c60e01bec089183e511efa27327414035))
* **sw:** reassociate affine-gap recurrences on NEON (kswv + bandedSWA) ([#166](https://github.com/fg-labs/bwa-mem3/issues/166)) ([a02fcb4](https://github.com/fg-labs/bwa-mem3/commit/a02fcb446574d5b5d03abdbf73c9b129deead2d4))


### Refactoring

* **bsw:** derive extension gaps from H (standard Gotoh), not M ([#141](https://github.com/fg-labs/bwa-mem3/issues/141)) ([f715fbd](https://github.com/fg-labs/bwa-mem3/commit/f715fbd92b4c4acabe7b7ab8dde35d5ad7fe211c))
* **bsw:** drop the dead qlen[] parameter from the 8-bit kernels ([#143](https://github.com/fg-labs/bwa-mem3/issues/143)) ([42321df](https://github.com/fg-labs/bwa-mem3/commit/42321df6100f04451b8a4d5ed7f214ea28522ec8))
* **bsw:** remove dead SW code paths (SORT_PAIRS, non-CORE macros, SSE2 polyfill) ([#148](https://github.com/fg-labs/bwa-mem3/issues/148)) ([ad8937e](https://github.com/fg-labs/bwa-mem3/commit/ad8937ef05fb48992ca61be70c21d4b6e245a7b3))


### Documentation

* add memory budgeting and data-type tuning guide ([#145](https://github.com/fg-labs/bwa-mem3/issues/145)) ([7127d80](https://github.com/fg-labs/bwa-mem3/commit/7127d80f7204c7951bf74d8191a883c4e0f00a82))
* add situational --supp-rep-hard-cap 20 note for SV-aware pipelines ([#134](https://github.com/fg-labs/bwa-mem3/issues/134)) ([e20bc0e](https://github.com/fg-labs/bwa-mem3/commit/e20bc0e28d3501a732e9b1eba88a9920ec877bae))
* **perf:** refresh reference-architecture table to v0.3.0 (a02fcb4) ([#167](https://github.com/fg-labs/bwa-mem3/issues/167)) ([2cd7bfa](https://github.com/fg-labs/bwa-mem3/commit/2cd7bfa1cfb5392856905df3b88cebeda2a1d659))
* **perf:** what drives the speedup, full perf-PR catalog, and fix the stale RTD build ([#155](https://github.com/fg-labs/bwa-mem3/issues/155)) ([0b01fb7](https://github.com/fg-labs/bwa-mem3/commit/0b01fb73a2e9a8b9e8899f30b2afd94c9e6a1418))
* recommend -s 0 for --meth Pass-2 re-seeding in settings profiles ([#132](https://github.com/fg-labs/bwa-mem3/issues/132)) ([45e02e0](https://github.com/fg-labs/bwa-mem3/commit/45e02e0095f862c99eaf99273b2912ff2229ec9b))
* recommend a recent compiler on aarch64, with measured NEON numbers ([#165](https://github.com/fg-labs/bwa-mem3/issues/165)) ([b591684](https://github.com/fg-labs/bwa-mem3/commit/b59168479ddb0fddba86cd0eb0e61ee0a1c8fcb0))
* settings profiles (bwa drop-in vs recommended) ([#131](https://github.com/fg-labs/bwa-mem3/issues/131)) ([4d845c9](https://github.com/fg-labs/bwa-mem3/commit/4d845c92154d29f944191ec080009d916ce4b72b))

## [0.2.2](https://github.com/fg-labs/bwa-mem3/compare/v0.2.1...v0.2.2) (2026-06-08)


### Bug Fixes

* **lto:** pass explicit -flto=N on GCC to bypass jobserver under sandboxes ([#122](https://github.com/fg-labs/bwa-mem3/issues/122)) ([c6240a7](https://github.com/fg-labs/bwa-mem3/commit/c6240a74959e35a8fcbaa29a6f3b78c6a67ec7c3))
* **smem:** free lockstep SMEM caches at thread exit (closes [#116](https://github.com/fg-labs/bwa-mem3/issues/116)) ([#117](https://github.com/fg-labs/bwa-mem3/issues/117)) ([9454f10](https://github.com/fg-labs/bwa-mem3/commit/9454f106cba06ad3b139a01e5fc2b73a24356ce7))


### Performance

* **sort:** stabilize alnreg tie-breaks + drop in pdqsort at dedup-patch sort sites ([#123](https://github.com/fg-labs/bwa-mem3/issues/123)) ([85f8542](https://github.com/fg-labs/bwa-mem3/commit/85f8542080863c2fda9c198cb8743d18803e83e8))


### Documentation

* **bench:** inject generated divergence catalog + per-release concordance table ([#126](https://github.com/fg-labs/bwa-mem3/issues/126)) ([fea1c94](https://github.com/fg-labs/bwa-mem3/commit/fea1c9407a71b4e6c0efbd12aaa6a3a98b85861a))
* correct concordance claims and document supplementary-alignment divergence ([#125](https://github.com/fg-labs/bwa-mem3/issues/125)) ([8b2dc69](https://github.com/fg-labs/bwa-mem3/commit/8b2dc697c926c779b03be6369e8a08e95283b284))
* document bwa-mem3&lt;-&gt;bwa-mem2 non-bit-identity + auditable PR list ([#124](https://github.com/fg-labs/bwa-mem3/issues/124)) ([bffae5a](https://github.com/fg-labs/bwa-mem3/commit/bffae5a09267877fe514c458d4956b717bcefb8f))

## [0.2.1](https://github.com/fg-labs/bwa-mem3/compare/v0.2.0...v0.2.1) (2026-05-17)


### Bug Fixes

* **changelog:** strip preamble so release-please owns the file ([#112](https://github.com/fg-labs/bwa-mem3/issues/112)) ([56e580c](https://github.com/fg-labs/bwa-mem3/commit/56e580cf4f2515280440b2882a32c3a6a7b6d15c))
* **mapq:** propagate SMEM SA-count to seed n_hits so --supp-rep-hard-cap works ([#101](https://github.com/fg-labs/bwa-mem3/issues/101)) ([cca9d4f](https://github.com/fg-labs/bwa-mem3/commit/cca9d4f41023e501636d2b1ebd2d1e825f95e8e3))
* **smem:** track enc_qdb byte capacity separately from wsize_mem ([#100](https://github.com/fg-labs/bwa-mem3/issues/100)) ([ab922b6](https://github.com/fg-labs/bwa-mem3/commit/ab922b6af531ed4521b6cd9bd3c43a661d51a421))


### Documentation

* **readme:** add bioconda badges and install instructions ([#106](https://github.com/fg-labs/bwa-mem3/issues/106)) ([830276c](https://github.com/fg-labs/bwa-mem3/commit/830276ce01774805c20bfff69c63cfebab239166))
