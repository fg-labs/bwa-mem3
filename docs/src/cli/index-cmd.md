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

```bash
bwa-mem3 index ref.fa
```

Build a methylation-aware index (required before `bwa-mem3 mem --meth`):

```bash
bwa-mem3 index --meth ref.fa
```

Limit peak RAM to 16 GB and write scratch data to `/scratch`:

```bash
bwa-mem3 index --max-memory 16G --tmp-dir /scratch ref.fa
```

## Flag reference

### `-p STR` — output prefix

By default, index files are written alongside `<in.fasta>` using the FASTA
path as a prefix (e.g. `ref.fa.bwt.2bit.64`, `ref.fa.pac`, etc.). Use `-p`
to write them to a different base path, such as a dedicated index directory:

```bash
bwa-mem3 index -p /idx/hg38 ref.fa
# writes /idx/hg38.bwt.2bit.64, /idx/hg38.pac, …
# align with: bwa-mem3 mem /idx/hg38 R1.fq R2.fq
```

### `-t INT` — worker threads

Controls the number of threads used during index construction. The default
auto-detects available cores and is cgroup-aware, so it behaves correctly
inside containers and on shared cluster nodes. Set explicitly when you want to
cap CPU usage.

### `--max-memory SIZE` — peak memory budget

Limits how much RAM the indexer may use at once. `SIZE` accepts a `G`, `M`, or
`K` suffix (case-insensitive) or a bare byte count. The default is the memory
available to the process less a reserve of `min(max(2 GB, 5%), 50%)`, computed in
a cgroup-aware manner — index construction is a one-shot batch job, so the
default lets it use the host it was given. The reserve is the larger of 2 GB and
5% of total, except that it never exceeds half of total, so a host below 4 GB
still resolves to a usable budget (a 1 GB host budgets 512 MB) instead of zero.

Suffix-array construction is **not** currently bounded-memory: the budget is a
precondition that is checked before the build starts, not a target the builder
spills to meet. A reference whose estimated requirement exceeds the budget is
rejected with exit code 3 and a message naming both the `--max-memory` value
that would clear it and the host RAM that would clear it automatically. Human
genomes need roughly 12 bytes per base of doubled text (8 below ~1.07 Gbp, where
the suffix array still fits 32-bit entries) — about 72 GB for hg38 —
so plan on a host with ~76 GB of RAM or more.

### `--tmp-dir PATH` — scratch directory

Scratch directory for intermediate files. Defaults to `$TMPDIR`. The builder
stages a doubled `.pac` here during construction, so point this at a fast local
disk (NVMe or ramdisk) with room for ~2 bits per base of doubled text.

### `--emit-unpacked-ref` — also write `<prefix>.0123`

Off by default. `bwa-mem3 mem` reconstructs reference bases from the packed
`.pac` on demand (*pac-fetch*), so the unpacked `.0123` (~8× the `.pac`; ~6.4 GB
on hg38) is never read and is not built. Enable this flag only when an external
consumer still requires the file — for example, sharing an index with
[bwa-mem2](../related-projects/bwa-mem2.md), which loads `.0123` directly:

```bash
bwa-mem3 index --emit-unpacked-ref ref.fa   # additionally writes ref.fa.0123
```

For `--meth` the flag applies to the **original** index only; the `.meth` seed
index never needs an unpacked reference.

### `--meth` — build a methylation (dual) index

Builds a **dual index**: the normal FM-index over the original FASTA (at the bare
prefix), plus a converted **seed** FM-index under the `.meth` prefix, built over a
per-strand-converted FASTA `<in.fasta>.meth.fa` (`f`-prefixed C→T and `r`-prefixed
G→A doubled contigs). All files are placed alongside the original FASTA.

By default, neither index writes an unpacked `.0123`: `mem --meth` extends against
the original reference (whose bases it pac-fetches from `.pac`), never the seed, so
the original `.0123` (~6.4 GB) is unnecessary and the seed's unpacked bases are
never read at all (~13 GB). Not building either saves ~19 GB of disk on hg38; the
runtime RSS reduction comes from avoiding the original `.0123` load (~6.4 GB).
(`--emit-unpacked-ref` overrides this for the **original** index only; the
`.meth` seed never needs one.)

Pass the **original** FASTA prefix to all three `index`, `shm`, and `mem` commands;
the `.meth` seed index is located automatically when `--meth` is present.

## Notes / Gotchas

> **Tip — Index once, align many times**
>
> A standard hg38 index is ~11 GB of index files on disk and takes several
> minutes to build. A `--meth` build adds the seed index on top — the doubled
> seed FM-index (~21 GB) plus its packed `.pac` (~1.6 GB) — for roughly **34 GB**
> of index files (~37 GB including the converted `.meth.fa`). That is a little
> over double the plain footprint, not triple: by default no unpacked `.0123` is
> built for either index (`--emit-unpacked-ref` adds it for the original only).
> Build once and store on shared storage; all alignment jobs on the same reference
> share the files.
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
