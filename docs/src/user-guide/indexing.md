# Indexing the reference

Before aligning reads, bwa-mem3 builds an FM-index from the reference FASTA.
The index is read back from disk at the start of every `mem` run, so it is
built once and reused indefinitely.

## Basic indexing

```sh
bwa-mem3 index ref.fa
```

The command writes five files alongside the input FASTA:

| File | Contents |
|------|----------|
| `ref.fa.bwt.2bit.64` | Burrows-Wheeler Transform, 2-bit packed, 64-bit offsets |
| `ref.fa.0123` | Forward sequence, 2-bit packed |
| `ref.fa.amb` | Coordinates and counts of ambiguous (N) bases |
| `ref.fa.ann` | Sequence names and lengths |
| `ref.fa.pac` | Forward sequence, 4-bit packed |

The `.bwt.2bit.64` file dominates disk usage. For the human reference (hg38),
expect roughly 28 GB total across all five files.

## Methylation index (`--meth`)

```sh
bwa-mem3 index --meth ref.fa
```

Methylation mode builds a C-to-T doubled reference in addition to the standard
FM-index files. The command writes a `ref.fa.bwameth.c2t` file (the doubled
FASTA) and its own set of five index files with the `.bwameth.c2t` suffix:

```text
ref.fa.bwameth.c2t
ref.fa.bwameth.c2t.bwt.2bit.64
ref.fa.bwameth.c2t.0123
ref.fa.bwameth.c2t.amb
ref.fa.bwameth.c2t.ann
ref.fa.bwameth.c2t.pac
```

The doubled reference is roughly twice the size of the standard one. For hg38,
allow approximately 56 GB of disk space.

> **Tip — Pass the original FASTA to mem, not the c2t file**
>
> When running `bwa-mem3 mem --meth`, pass the original FASTA path (`ref.fa`),
> not `ref.fa.bwameth.c2t`. bwa-mem3 appends `.bwameth.c2t` automatically.
> The auto-append is skipped only when the path already ends in `.bwameth.c2t`,
> which is useful for external-c2t interop pipelines.

## Output file locations

Index files are written to the same directory as the input FASTA by default.
The input path is taken verbatim as a prefix — you can pass an absolute path to
write into a different directory:

```sh
bwa-mem3 index /data/indexes/hg38/hg38.fa
# writes hg38.fa.bwt.2bit.64, etc. into /data/indexes/hg38/
```

## Time and memory

Indexing hg38 takes roughly 60–90 minutes on a single core and requires about
80 GB of RAM during construction. The process is single-threaded; additional
cores do not reduce wall time.

bwa-mem3 uses [libsais](https://github.com/IlyaGrebnov/libsais) to construct
the suffix array, which is faster than the original bwa-mem2 approach. See
[Performance improvements](../whats-different/performance.md) for benchmark
numbers.

> **Warning — Do not index over a live shared-memory segment**
>
> If you have previously staged the index into shared memory with `bwa-mem3 shm`,
> drop the segment first before re-indexing:
>
> ```sh
> bwa-mem3 shm -d
> bwa-mem3 index ref.fa
> ```
>
> There is no staleness check. If `bwa-mem3 mem` finds a matching segment in
> shared memory it will attach to it even when the on-disk index has been updated.
> See [Quick start: shared-memory index](../getting-started/quick-shm.md).

## Arch flags and the index format

The FM-index format is architecture-independent. A single index works
across every SIMD tier and every supported platform: the x86 binary's
AVX2 / AVX-512BW dispatch paths and the arm64 NEON binary all read the
same on-disk layout.

---

**See also:**
[Quick start: align paired-end FASTQs](../getting-started/quick-align.md) ·
[Quick start: methylation alignment](../getting-started/quick-meth.md) ·
[Quick start: shared-memory index](../getting-started/quick-shm.md) ·
[Performance improvements](../whats-different/performance.md) ·
[CLI Reference: index](../cli/index-cmd.md)
