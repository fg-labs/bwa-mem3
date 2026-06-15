# bwa-mem3

[![CI](https://github.com/fg-labs/bwa-mem3/actions/workflows/ci.yml/badge.svg?branch=main)](https://github.com/fg-labs/bwa-mem3/actions/workflows/ci.yml)
[![codecov](https://codecov.io/gh/fg-labs/bwa-mem3/branch/main/graph/badge.svg)](https://codecov.io/gh/fg-labs/bwa-mem3)
[![Bioconda](https://img.shields.io/conda/vn/bioconda/bwa-mem3.svg?label=bioconda)](https://anaconda.org/bioconda/bwa-mem3)
[![Documentation](https://img.shields.io/readthedocs/bwa-mem3?label=docs)](https://bwa-mem3.readthedocs.io)
[![License](https://img.shields.io/badge/license-MIT-blue.svg)](https://github.com/fg-labs/bwa-mem3/blob/main/LICENSE)

bwa-mem3 is a short-read aligner derived from [bwa-mem2](https://github.com/bwa-mem2/bwa-mem2),
carrying correctness fixes, performance improvements, and new features (methylation alignment,
shared-memory index, mimalloc allocator) maintained by [Fulcrum Genomics](https://fulcrumgenomics.com).

bwa-mem3 is **not** byte-identical to bwa-mem2 — it adds SAM tags, fixes crashes, and changes tie resolution. The core alignment is preserved on the data we have tested, but the SAM byte stream is not. See [Equivalence with bwa-mem2](https://bwa-mem3.readthedocs.io/en/latest/whats-different/equivalence.html) for a field-by-field comparison and a full per-PR audit.

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

- [bwa-mem3-bench](https://github.com/fg-labs/bwa-mem3-bench) — benchmarking harness across CPU architectures
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
