# bwa-mem3

[![CI](https://github.com/fg-labs/bwa-mem3/actions/workflows/ci.yml/badge.svg?branch=main)](https://github.com/fg-labs/bwa-mem3/actions/workflows/ci.yml)
[![codecov](https://codecov.io/gh/fg-labs/bwa-mem3/branch/main/graph/badge.svg)](https://codecov.io/gh/fg-labs/bwa-mem3)
[![Bioconda](https://img.shields.io/conda/vn/bioconda/bwa-mem3.svg?label=bioconda)](https://anaconda.org/bioconda/bwa-mem3)
[![Documentation](https://img.shields.io/readthedocs/bwa-mem3?label=docs)](https://bwa-mem3.readthedocs.io)
[![License](https://img.shields.io/badge/license-MIT-blue.svg)](https://github.com/fg-labs/bwa-mem3/blob/main/LICENSE)

bwa-mem3 is a short-read aligner derived from [bwa-mem2](https://github.com/bwa-mem2/bwa-mem2),
carrying correctness fixes, performance improvements, and new features (methylation alignment,
shared-memory index, mimalloc allocator) maintained by [Fulcrum Genomics](https://fulcrumgenomics.com).

## Performance

Wall-clock speedup of the current release (v0.18.0\*) against `bwa` 0.7.19, `bwa-mem2` v2.2.1, and `minibwa`, on the `wgs-5M` sample. Cells are `stock / --fast`.

| arch | wall_s | vs bwa | vs bwa-mem2 | vs minibwa |
|---|---:|---:|---:|---:|
| ARM | 28.42 / 12.44 | 8.82x / 20.15x | — | 1.43x / 3.27x |
| x86 | 15.44 / 8.95 | 12.85x / 22.16x | 5.72x / 9.87x | 1.94x / 3.36x |

\* v0.18.0, like v0.13.0 through v0.17.0, is measured with a **denser suffix-array index** — stride 2 instead of the stock stride 8 — built with `bwa-mem3 re-sa -u 1` (or staged with `bwa-mem3 shm -u 1`), both new in v0.13.0. Alignment records are byte-identical by construction (each added suffix-array sample is computed from the on-disk ones), and measured so; only the `@PG` header line differs, since its command line names the index path: a full-content `fgumi compare bams` of `wgs-5M` (10,030,558 records) aligned on the stride-2 index against the stock stride-8 index — v0.13.0 against the output-identical v0.12.0 — found 0 differences on both benchmark hosts, Graviton4 m8g (NEON) and AMD m8a (AVX-512), and the change that added it ([#510](https://github.com/fg-labs/bwa-mem3/pull/510)) found 0 differences in a record-by-record comparison of core SAM fields and tags on 5M-pair hg38 WGS and WES slices (that PR does not record the host or SIMD tier). It costs ~12 GB more memory for the human genome. v0.13.0 through v0.18.0 are all measured on it in the tables below, so v0.18.0's `vs prev` is a like-for-like comparison. Every older release and every comparator uses the stock index.

> [!TIP]
> **📈 Full release-history table** — every bwa-mem3 release since v0.2.1, full methodology, and version pins.
>
> <details>
> <summary><strong>Click to expand</strong></summary>
>
> **Graviton4 (m8g, arm64/NEON)**
>
> | release | wall_s | vs bwa | vs bwa-mem2 | vs minibwa | vs prev |
> |---|---:|---:|---:|---:|---:|
> | bwa | 250.66 | 1.00x | — | 0.16x | — |
> | bwa-mem2 | — | — | — | — | — |
> | minibwa | 40.66 | 6.16x | — | 1.00x | — |
> | v0.2.1 | 120.63 | 2.08x | — | 0.34x | — |
> | v0.2.2 | 118.26 | 2.12x | — | 0.34x | 1.020x |
> | v0.3.0 | 102.35 | 2.45x | — | 0.40x | 1.155x |
> | v0.4.0 | 104.10 | 2.41x | — | 0.39x | 0.983x |
> | v0.5.0 | 104.20 / 36.73 | 2.41x / 6.82x | — | 0.39x / 1.11x | 0.999x |
> | v0.6.0 | 97.40 / 36.22 | 2.57x / 6.92x | — | 0.42x / 1.12x | 1.070x / 1.014x |
> | v0.7.0 | 91.45 / 39.33 | 2.74x / 6.37x | — | 0.44x / 1.03x | 1.065x / 0.921x |
> | v0.8.0 | 76.42 / 27.71 | 3.28x / 9.05x | — | 0.53x / 1.47x | 1.197x / 1.419x |
> | v0.9.0 | 76.39 / 27.85 | 3.28x / 9.00x | — | 0.53x / 1.46x | 1.000x / 0.995x |
> | v0.10.0 | 70.85 / 27.04 | 3.54x / 9.27x | — | 0.57x / 1.50x | 1.078x / 1.030x |
> | v0.11.0 | 63.96 / 25.52 | 3.92x / 9.82x | — | 0.64x / 1.59x | 1.108x / 1.060x |
> | v0.12.0 | 58.05 / 18.70 | 4.32x / 13.40x | — | 0.70x / 2.17x | 1.102x / 1.365x |
> | v0.13.0\* | 54.32 / 18.10 | 4.61x / 13.85x | — | 0.75x / 2.25x | 1.069x / 1.033x |
> | v0.14.0\* | 41.50 / 17.17 | 6.04x / 14.60x | — | 0.98x / 2.37x | 1.309x / 1.054x |
> | v0.14.1\* | 35.08 / 15.93 | 7.15x / 15.74x | — | 1.16x / 2.55x | 1.183x / 1.078x |
> | v0.15.0\* | 34.16 / 15.42 | 7.34x / 16.26x | — | 1.19x / 2.64x | 1.027x / 1.033x |
> | v0.15.1\* | 31.40 / 14.12 | 7.98x / 17.75x | — | 1.29x / 2.88x | 1.088x / 1.092x |
> | v0.16.0\* | 28.81 / 12.64 | 8.70x / 19.83x | — | 1.41x / 3.22x | 1.090x / 1.117x |
> | v0.17.0\* | 28.39 / 12.48 | 8.83x / 20.08x | — | 1.43x / 3.26x | 1.015x / 1.013x |
> | **v0.18.0**\* | **28.42 / 12.44** | **8.82x / 20.15x** | **—** | **1.43x / 3.27x** | **0.999x / 1.003x** |
>
> **AMD (m8a, x86-64/AVX-512)**
>
> | release | wall_s | vs bwa | vs bwa-mem2 | vs minibwa | vs prev |
> |---|---:|---:|---:|---:|---:|
> | bwa | 198.33 | 1.00x | 0.45x | 0.15x | — |
> | bwa-mem2 | 88.30 | 2.25x | 1.00x | 0.34x | — |
> | minibwa | 30.03 | 6.60x | 2.94x | 1.00x | — |
> | v0.2.1 | 54.20 | 3.66x | 1.63x | 0.55x | 1.629x |
> | v0.2.2 | 52.18 | 3.80x | 1.69x | 0.58x | 1.039x |
> | v0.3.0 | 49.12 | 4.04x | 1.80x | 0.61x | 1.062x |
> | v0.4.0 | 51.04 | 3.89x | 1.73x | 0.59x | 0.962x |
> | v0.5.0 | 50.92 / 22.08 | 3.89x / 8.98x | 1.73x / 4.00x | 0.59x / 1.36x | 1.002x |
> | v0.6.0 | 49.56 / 21.63 | 4.00x / 9.17x | 1.78x / 4.08x | 0.61x / 1.39x | 1.027x / 1.021x |
> | v0.7.0 | 47.71 / 21.89 | 4.16x / 9.06x | 1.85x / 4.03x | 0.63x / 1.37x | 1.039x / 0.988x |
> | v0.8.0 | 42.78 / 17.76 | 4.64x / 11.17x | 2.06x / 4.97x | 0.70x / 1.69x | 1.115x / 1.233x |
> | v0.9.0 | 42.89 / 17.67 | 4.62x / 11.22x | 2.06x / 5.00x | 0.70x / 1.70x | 0.997x / 1.005x |
> | v0.10.0 | 42.30 / 17.63 | 4.69x / 11.25x | 2.09x / 5.01x | 0.71x / 1.70x | 1.014x / 1.002x |
> | v0.11.0 | 41.07 / 16.80 | 4.83x / 11.81x | 2.15x / 5.26x | 0.73x / 1.79x | 1.030x / 1.049x |
> | v0.12.0 | 32.32 / 13.92 | 6.14x / 14.25x | 2.73x / 6.34x | 0.93x / 2.16x | 1.271x / 1.207x |
> | v0.13.0\* | 29.24 / 13.24 | 6.78x / 14.98x | 3.02x / 6.67x | 1.03x / 2.27x | 1.105x / 1.051x |
> | v0.14.0\* | 23.17 / 12.12 | 8.56x / 16.36x | 3.81x / 7.29x | 1.30x / 2.48x | 1.262x / 1.092x |
> | v0.14.1\* | 21.39 / 11.93 | 9.27x / 16.62x | 4.13x / 7.40x | 1.40x / 2.52x | 1.083x / 1.016x |
> | v0.15.0\* | 20.30 / 11.11 | 9.77x / 17.85x | 4.35x / 7.95x | 1.48x / 2.70x | 1.054x / 1.074x |
> | v0.15.1\* | 18.44 / 10.38 | 10.76x / 19.11x | 4.79x / 8.51x | 1.63x / 2.89x | 1.101x / 1.070x |
> | v0.16.0\* | 15.96 / 9.06 | 12.43x / 21.89x | 5.53x / 9.75x | 1.88x / 3.31x | 1.155x / 1.146x |
> | v0.17.0\* | 15.25 / 9.04 | 13.01x / 21.94x | 5.79x / 9.77x | 1.97x / 3.32x | 1.047x / 1.002x |
> | **v0.18.0**\* | **15.44 / 8.95** | **12.85x / 22.16x** | **5.72x / 9.87x** | **1.94x / 3.36x** | **0.988x / 1.010x** |
>
> `vs prev` is the release-over-release speedup (`prev_wall / this_wall`, `>1` = faster) vs the previous release on this same host, `stock / --fast`. Since v0.16.0 every bwa-mem3 release in these tables, old ones included, writes BAM natively (`--bam=0`) while `bwa`, `bwa-mem2` and `minibwa` pipe SAM through `samtools view -u`; earlier tables sent every arm through SAM, so `vs prev` stays like-for-like but the comparator ratios are not comparable across that change. Since v0.17.0 every bwa-mem3 release whose Makefile supports link-time optimization (LTO), v0.15.1 onward, is built with it, old ones included; releases before v0.15.1 predate it and are built without. Every row was re-measured, so neither absolute times nor ratios are comparable with earlier versions of this table. The first release's predecessor is upstream `bwa-mem2` — bwa-mem3 is its successor — so v0.2.1's `vs prev` is its speedup over bwa-mem2 (blank on ARM, where upstream has no build).
>
> Version pins: `bwa` 0.7.19 · `bwa-mem2` v2.2.1 · `minibwa` commit [`d6d9f87d`](https://github.com/lh3/minibwa) (`minibwa-0.7`). "ARM" = Graviton4 m8g (arm64/NEON, no SMT); "x86" = AMD m8a (x86-64/AVX-512, no SMT). Both are the general-purpose siblings of the c8g/c8a hosts used through v0.12.0 — same CPU family, same core count, no SMT, but 4 GiB/vCPU instead of 2 so every historical arm fits in memory — so absolute times are not comparable with earlier versions of this table, only ratios within it. (The x86 arm replaced an earlier Intel c7i arm, which ran 16 vCPUs over 8 physical cores under 2-way SMT and so wasn't a real core-for-core match for Graviton's 16 real cores.) No ARM `bwa-mem2` build exists, hence the blank cells there. Every arm for a given arch ran interleaved on one fixed on-demand host — 5 rounds for v0.13.0 through v0.18.0, 2 for the rest, median wall-clock shown — so these are same-host comparisons, not medians pooled across separate runs. `—` means the release predates the comparator or predates `--fast`. \* v0.13.0 through v0.18.0 rows use the stride-2 suffix-array index described above; all other rows use the stock stride-8 index. Regenerated at each release; see [Benchmarks](https://bwa-mem3.readthedocs.io/en/latest/performance/benchmarks.html).
>
> </details>

> [!WARNING]
> `--fast` is **not alignment-identical** to the default preset — it trades some sensitivity/specificity at the extremes (repetitive/multi-mapping regions, low-`MAPQ` reads) for the speedup above. See "Three ways to run it" below before switching a production pipeline to it.

## Three ways to run it — plain, `--compat`, `--fast`

bwa-mem3 has three alignment modes that differ in *what alignments come out*, not just in speed:

| mode | where reads align | when to use |
|---|---|---|
| **plain** (default) | bwa-mem2's alignments **plus bonafide correctness fixes**, with two extra tags (`MQ:i`, `HN:i`) and an enriched header. On the cells re-measured for release 0.7.1, the complete alignment-record stream (tags stripped) is byte-identical to bwa-mem2 v2.2.1 on `wgs-5M`/`wes-5M`/`hic-1M` (x86 `c6a` AVX2, with a `c6a`/`c8g` cross-arch check confirming the Arm `c8g` NEON build matches) — differing only by those additive tags and the header. Separately, a 1.07M-record HG00096 WGS slice shows zero diverging **primary** alignments (x86, primary-only; not part of the cross-arch or complete-stream checks). | Migrating a pipeline, validating against bwa/bwa-mem2, or any new pipeline. |
| **`--compat=bwa-mem2` / `--compat=bwa-mem`** | Byte-for-byte identical **alignment records** to a **specific** upstream (bwa-mem2 v2.2.1 or bwa 0.7.19), `@PG` excluded and `-t`/`-K` matched. The two targets are **not** interchangeable. | Diff-clean validation against an existing bwa/bwa-mem2 golden. |
| **`--fast`** | Faster, and **not** record-compatible with the default: it reshuffles the low-confidence tail (~85% of the reads it re-places had `MAPQ 0`; the confident `MAPQ`-60 core moves on ≤0.5%, 0.011% on `wgs-5M`) while staying accuracy-neutral against golden truth (≤0.02 pp across the WGS and methylation sims). Figures from the [benchmark](https://bwa-mem3.readthedocs.io/en/latest/performance/benchmarks.html) release-validation cells (`wgs-5M`/`wes-5M`/`panel-twist-5M` at 5 M reads, `hic-1M`/`sbx-1M` at 1 M) across every SIMD tier (AVX2 `c6a`, AVX-512 `c7a`/`c7i`, NEON `c7g`/`c8g`; meth on `m7i`), each a `.4xlarge` host at `-t 16`, `-K 160000000`. | High-throughput pipelines where you care about the confident, uniquely-mapped calls. |

`--compat` is mutually exclusive with `--fast` (and with `--meth` and `--proper-pair-from-emitted`). See [Alignment modes](https://bwa-mem3.readthedocs.io/en/latest/whats-different/modes.html) for the full side-by-side and [Equivalence with bwa-mem2](https://bwa-mem3.readthedocs.io/en/latest/whats-different/equivalence.html) for the field-by-field audit.

By default bwa-mem3 keeps bwa-mem2's command-line defaults, so it drops into an existing pipeline unchanged. For the fastest configuration — and what each recommended deviation from the bwa defaults trades for speed — see [Settings profiles: bwa drop-in vs recommended](https://bwa-mem3.readthedocs.io/en/latest/best-practices/settings-profiles.html).

**Full documentation:** <https://bwa-mem3.readthedocs.io>

## Install

The recommended way to install bwa-mem3 is via [bioconda](https://bioconda.github.io):

```sh
mamba install -c bioconda bwa-mem3
bwa-mem3 version
```

Prebuilt packages are available for `linux-64`, `linux-aarch64`, and `osx-arm64`.

### Build from source

```sh
git clone --recursive https://github.com/fg-labs/bwa-mem3.git
cd bwa-mem3
make
./bwa-mem3 version
```

See the [installation guide](https://bwa-mem3.readthedocs.io/en/latest/getting-started/installation.html) for prerequisites and architecture-specific notes.

## Quick links

- [Benchmarks](https://bwa-mem3.readthedocs.io/en/latest/performance/benchmarks.html) — published results for every release, and how they are measured
- [Shared-memory index](https://bwa-mem3.readthedocs.io/en/latest/cli/shm.html) — `bwa-mem3 shm`; on Linux, mount `/dev/shm` with `huge=advise` before staging, or alignment uses more CPU than a disk load (hg38, 5M read pairs, `-t 16`: +33 % on Zen 4 x86-64 AVX-512BW, +26 % on Graviton4 arm64 NEON; [measurements](https://bwa-mem3.readthedocs.io/en/latest/cli/shm.html#huge-pages-mount-devshm-with-hugeadvise-or-hugewithin_size))
- [bwa-mem3-rs](https://github.com/fg-labs/bwa-mem3-rs) — Rust bindings for bwa-mem3
- [bioconda recipe](https://github.com/bioconda/bioconda-recipes/tree/master/recipes/bwa-mem3) — conda package on bioconda
- [fgumi](https://github.com/fulcrumgenomics/fgumi) — UMI-aware consensus and deduplication
- [bwa-mem2](https://github.com/bwa-mem2/bwa-mem2) — upstream project

## License

MIT. See the [License page](https://bwa-mem3.readthedocs.io/en/latest/reference/license.html) in the docs.

## Citation

Please cite the bwa-mem2 paper (Vasimuddin Md et al., IPDPS 2019). See the [Citation page](https://bwa-mem3.readthedocs.io/en/latest/reference/citation.html) for BibTeX.

## Issues / contributing

File [issues](https://github.com/fg-labs/bwa-mem3/issues) and [pull requests](https://github.com/fg-labs/bwa-mem3/pulls) on [fg-labs/bwa-mem3](https://github.com/fg-labs/bwa-mem3).
