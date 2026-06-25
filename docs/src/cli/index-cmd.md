# index

`bwa-mem3 index` builds the FM-index (BWT + suffix array) that `bwa-mem3 mem`
requires for alignment. Run it once per reference; the resulting files sit
alongside the input FASTA and are reused for all subsequent alignment jobs.
Pass `--meth` to build a dual index — the normal index plus a converted `.meth`
seed index — for bisulfite-seq alignment.

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

### `--meth` — build a methylation (dual) index

Builds a **dual index**: the normal FM-index over the original FASTA (at the bare
prefix), plus a converted **seed** FM-index under the `.meth` prefix, built over a
per-strand-converted FASTA `<in.fasta>.meth.fa` (`f`-prefixed C→T and `r`-prefixed
G→A doubled contigs). All files are placed alongside the original FASTA.

The seed index omits the unpacked `.meth.0123` reference: `mem --meth` extends
against the original reference, never the seed, so the seed's unpacked bases are
never read. Not building it saves ~13 GB of disk (and RSS) on hg38.

Pass the **original** FASTA prefix to all three `index`, `shm`, and `mem` commands;
the `.meth` seed index is located automatically when `--meth` is present.

## Notes / Gotchas

> **Tip — Index once, align many times**
>
> A standard hg38 index is ~18 GB of index files on disk and takes several
> minutes to build. A `--meth` build adds the seed index on top — the doubled
> seed FM-index (~21 GB) plus its packed `.pac` (~1.6 GB) — for roughly **40 GB**
> of index files (~47 GB including the converted `.meth.fa`). That is a little
> over double the plain footprint, not triple: the seed's unpacked `.0123`
> (~13 GB) is no longer built. Build once and store on shared storage; all
> alignment jobs on the same reference share the files.
>
> **Note — a `--meth` index is a superset, not a separate index**
>
> `index --meth` writes the normal index at the bare prefix *plus* the `.meth`
> seed index. The bare-prefix index is an ordinary index, so `bwa-mem3 mem ref.fa`
> (without `--meth`) works fine for standard alignment against the same files — no
> separate index directory is needed. Only `--meth` runs use the `.meth` seed
> index.

---

**See also:**
[User Guide — Indexing the reference](../user-guide/indexing.md) ·
[CLI Reference — mem](mem.md) ·
[CLI Reference — shm](shm.md) ·
[Getting Started — Quick start: methylation alignment](../getting-started/quick-meth.md) ·
[Methylation Reference — Overview](../methylation/overview.md)
