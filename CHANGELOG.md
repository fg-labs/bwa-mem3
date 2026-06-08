# Changelog

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

## Changelog
