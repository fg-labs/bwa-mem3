# bwa-mem3

> Accelerated short-read alignment, derived from [bwa-mem2](https://github.com/bwa-mem2/bwa-mem2).

[![CI](https://github.com/fg-labs/bwa-mem3/actions/workflows/ci.yml/badge.svg?branch=main)](https://github.com/fg-labs/bwa-mem3/actions/workflows/ci.yml)
[![codecov](https://codecov.io/gh/fg-labs/bwa-mem3/branch/main/graph/badge.svg)](https://codecov.io/gh/fg-labs/bwa-mem3)

## About

bwa-mem3 is a derivative of [bwa-mem2](https://github.com/bwa-mem2/bwa-mem2),
the Burrows-Wheeler Aligner short-read aligner. The original algorithm and
codebase were authored by Vasimuddin Md (@yuk12), Sanchit Misra
(@sanchit-misra), and contributors. bwa-mem3 carries a series of
correctness, performance, and feature improvements maintained by
[Fulcrum Genomics](https://fulcrumgenomics.com). See
[`docs/whats-different.md`](docs/whats-different.md) for the running list
of changes.

## Install

A Bioconda package is coming soon. For now, build from source:

```sh
git clone --recursive https://github.com/fg-labs/bwa-mem3
cd bwa-mem3
make
./bwa-mem3
```

bwa-mem3 bundles [mimalloc](https://github.com/microsoft/mimalloc) and links
it into every binary by default. To build without mimalloc, pass
`USE_MIMALLOC=0`:

```sh
make USE_MIMALLOC=0
```

To confirm mimalloc is active, run `bwa-mem3 version` — a line like
`mimalloc 3.x.x` is printed when the allocator is linked in.

## Usage

Index a reference, then align:

```sh
./bwa-mem3 index ref.fa
./bwa-mem3 mem -t 16 ref.fa r1.fq.gz r2.fq.gz > out.sam
```

### Methylation alignment (`--meth`)

bwa-mem3 is a single-binary drop-in replacement for the `bwameth.py`
pipeline. No Python, no piped preprocessing — one `bwa-mem3 index --meth`
for the reference, one `bwa-mem3 mem --meth` for alignment.

Build a c2t doubled reference once:

```sh
bwa-mem3 index --meth ref.fa
# writes ref.fa.bwameth.c2t and its FMI alongside the input FASTA
```

Align raw FASTQs (paired-end or single-end):

```sh
bwa-mem3 mem --meth -t 16 ref.fa R1.fq.gz R2.fq.gz \
  | samtools sort -o out.bam
samtools index out.bam
```

`--meth` turns on:

- **Inline c2t read conversion** — R1 gets `C→T`, R2 gets `G→A`, with the
  original sequence stashed as a `YS:Z:` comment tag and the conversion
  direction as `YC:Z:` (bwameth.py convention; both pass through to SAM).
- **Auto-append `.bwameth.c2t`** to the reference path so the user passes
  the original FASTA prefix, not the c2t file.
- **bwameth.py-equivalent flag defaults:** `-B 2 -L 10 -U 100 -T 40 -CM`.
  Users can still override any of them.
- **Inline BAM post-processing:** strips `f`/`r` prefix from `@SQ`/`RNAME`
  consolidating to one `@SQ` per real chromosome, emits `YD:Z:{f,r}` on
  mapped records, chimera QC (longest `M`/`=`/`X` run < 44% of read length
  → `0x200`, clear `0x2`, cap MAPQ at 1), pair-level QC-fail propagation,
  and a `@PG ID:bwa-mem3-meth` entry. Output is uncompressed BAM via
  htslib (`wb0`) — near-free CPU cost, fully readable by `samtools`.

Additional options:

- `--set-as-failed {f,r}` — flag alignments aligned to the given strand
  as QC-fail (`0x200`).
- `--do-not-penalize-chimeras` — skip the longest-match < 44% chimera
  heuristic.

For callers that already do their own c2t conversion (e.g. an external
`bwameth.py c2t ... | bwa-mem3 mem ...` pipeline), pass the c2t reference
path directly — `mem --meth` detects the `.bwameth.c2t` suffix and skips
the auto-append:

```sh
bwameth.py c2t R1.fq.gz R2.fq.gz \
  | bwa-mem3 mem --meth -p -t 16 ref.fa.bwameth.c2t /dev/stdin \
  | samtools sort -o out.bam
```

### Shared-memory index (`shm`)

The bwa-mem3 FM-index is large (~28 GB for hg38), and `bwa-mem3 mem`
re-reads it from disk on every invocation. For workloads that align
many small samples back-to-back on the same machine, the index can be
staged once into POSIX shared memory:

```sh
./bwa-mem3 shm ref.fa     # stage; subsequent `mem` runs auto-attach
./bwa-mem3 shm --meth ref.fa  # stage a `--meth` index
./bwa-mem3 shm -l         # list staged indices
./bwa-mem3 shm -d         # drop everything
```

`bwa-mem3 mem ref.fa ...` automatically attaches when a matching segment
is staged; no extra flag is needed.

**Footgun: there is no staleness check.** If you re-run `bwa-mem3 index
ref.fa`, the on-disk files will not match the staged segment, but `mem`
will still attach to the stale segment and silently mis-align. Always
`bwa-mem3 shm -d` before re-indexing.

### Version

```sh
./bwa-mem3 version
```

## What's different from bwa-mem2

bwa-mem3 carries:

- A bundled `mimalloc` allocator linked by default (faster wall-clock
  on multi-threaded alignment; opt out with `USE_MIMALLOC=0`).
- Single-binary methylation alignment (`--meth`), a drop-in replacement
  for the `bwameth.py` pipeline.
- A `bwa-mem3 shm` subcommand for persistent shared-memory indices
  (ported from `bwa shm` in bwa v1).
- Performance and correctness improvements on top of upstream bwa-mem2.

See [`docs/whats-different.md`](docs/whats-different.md) for the full
running list.

## License

MIT, same as upstream bwa-mem2.

## Citation

If you use bwa-mem3, please also cite the original bwa-mem2 paper:

> Vasimuddin Md, Sanchit Misra, Heng Li, Srinivas Aluru.
> *Efficient Architecture-Aware Acceleration of BWA-MEM for Multicore
> Systems.* IEEE Parallel and Distributed Processing Symposium (IPDPS),
> 2019. [10.1109/IPDPS.2019.00041](https://doi.org/10.1109/IPDPS.2019.00041)

## Issues / contributing

File issues and pull requests at
[`fg-labs/bwa-mem3`](https://github.com/fg-labs/bwa-mem3/issues).
