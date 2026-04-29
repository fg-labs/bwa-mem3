# index

`bwa-mem3 index` builds the FM-index (BWT + suffix array) that `bwa-mem3 mem`
requires for alignment. Run it once per reference; the resulting files sit
alongside the input FASTA and are reused for all subsequent alignment jobs.
Pass `--meth` to produce a bwameth-compatible doubled c2t reference for
bisulfite-seq alignment.

## Synopsis

```text
{{#include ../../_generated/cli/index.txt}}
```

## Common usage

Build a standard index using all available cores:

```sh
bwa-mem3 index ref.fa
```

Build a methylation-aware index (required before `bwa-mem3 mem --meth`):

```sh
bwa-mem3 index --meth ref.fa
```

Limit peak RAM to 16 GB and write scratch data to `/scratch`:

```sh
bwa-mem3 index --max-memory 16G --tmp-dir /scratch ref.fa
```

## Flag reference

### `-p STR` — output prefix

By default, index files are written alongside `<in.fasta>` using the FASTA
path as a prefix (e.g. `ref.fa.bwt.2bit.64`, `ref.fa.0123`, etc.). Use `-p`
to write them to a different base path, such as a dedicated index directory:

```sh
bwa-mem3 index -p /idx/hg38 ref.fa
# writes /idx/hg38.bwt.2bit.64, /idx/hg38.0123, …
# align with: bwa-mem3 mem /idx/hg38 R1.fq R2.fq
```

### `-t INT` — worker threads

Controls the number of threads used during index construction. The default
auto-detects available cores and is cgroup-aware, so it behaves correctly
inside containers and on shared cluster nodes. Set explicitly when you want to
cap CPU usage.

### `--max-memory SIZE` — peak memory budget

Limits how much RAM the indexer may use at once. `SIZE` accepts a `G`, `M`, or
`K` suffix (case-insensitive) or a bare byte count. The default is
`min(50% of RAM, 32 GB)`, computed in a cgroup-aware manner.

For large references (hg38 and above) on machines with limited RAM, setting
this to a value lower than the reference size causes the indexer to partition
work and use `--tmp-dir` for intermediate files, at the cost of extra I/O.

### `--tmp-dir PATH` — scratch directory

Scratch directory for intermediate files when memory is partitioned. Defaults
to `$TMPDIR`. Point this at a fast local disk (NVMe or ramdisk) to minimize
wall-clock time when `--max-memory` forces partitioned construction.

### `--meth` — build a methylation (c2t) index

Writes a bwameth-style doubled reference — `<in.fasta>.bwameth.c2t` — and
builds the FM-index over that file rather than the original FASTA. The c2t
file and its index files are placed alongside the original FASTA.

Pass the **original** FASTA prefix (not the `.bwameth.c2t` path) to all three
`index`, `shm`, and `mem` commands. The c2t suffix is appended automatically
when `--meth` is present.

## Notes / Gotchas

> **Tip — Index once, align many times**
>
> Index construction for hg38 takes several minutes and ~28 GB of disk. Build
> the index once and store it on shared storage; all alignment jobs on the same
> reference share the same index files.
>
> **Warning — --meth index is not interchangeable with the standard index**
>
> A `--meth` index is built over the c2t reference and cannot be used for normal
> (non-bisulfite) alignment. Keep separate index directories if you align both
> standard and bisulfite samples to the same reference.

---

**See also:**
[User Guide — Indexing the reference](../user-guide/indexing.md) ·
[CLI Reference — mem](mem.md) ·
[CLI Reference — shm](shm.md) ·
[Getting Started — Quick start: methylation alignment](../getting-started/quick-meth.md) ·
[Methylation Reference — Overview](../methylation/overview.md)
