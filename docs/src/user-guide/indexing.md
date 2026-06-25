# Indexing the reference

Before aligning reads, bwa-mem3 builds an FM-index from the reference FASTA.
The index is read back from disk at the start of every `mem` run, so it is
built once and reused indefinitely.

## Basic indexing

```sh
bwa-mem3 index ref.fa
```

The command writes four files alongside the input FASTA:

| File | Contents |
|------|----------|
| `ref.fa.bwt.2bit.64` | Burrows-Wheeler Transform, 2-bit packed, 64-bit offsets |
| `ref.fa.amb` | Coordinates and counts of ambiguous (N) bases |
| `ref.fa.ann` | Sequence names and lengths |
| `ref.fa.pac` | Forward sequence, 2-bit packed (one base per 2 bits) |

The `.bwt.2bit.64` file dominates disk usage. For the human reference (hg38),
expect roughly 11 GB total across all four files.

> **Note — no `.0123` by default**
>
> Earlier releases (and bwa-mem2) also wrote `ref.fa.0123`, an *unpacked*
> forward+reverse reference (~8× the `.pac`; ~6.4 GB on hg38). bwa-mem3 no longer
> builds it: `mem` reconstructs the bases it needs directly from the packed
> `.pac` on demand (*pac-fetch*), so `.0123` is never read. Output is byte-for-byte
> identical. Pass `--emit-unpacked-ref` to `index` if you need the file for an
> external tool that still requires it (e.g. bwa-mem2):
>
> ```sh
> bwa-mem3 index --emit-unpacked-ref ref.fa   # also writes ref.fa.0123
> ```
>
> At run time, `BWAMEM3_REF_PAC_FETCH=0` reloads `.0123` instead of pac-fetching
> (only on an index that still has the file); it exists for A/B verification.

## Methylation index (`--meth`)

```sh
bwa-mem3 index --meth ref.fa
```

Methylation mode builds a **dual index**: the normal FM-index over the original
reference (at the bare prefix), plus a converted **seed** index under the `.meth`
prefix. The seed index is built over a per-strand-converted FASTA
(`ref.fa.meth.fa`) whose contigs are doubled (`f`-prefixed C→T, `r`-prefixed G→A):

```text
# normal index over the original reference (used for scoring/extension)
ref.fa.amb
ref.fa.ann
ref.fa.bwt.2bit.64
ref.fa.pac
# converted seed index (used only for seeding)
ref.fa.meth.fa
ref.fa.meth.amb
ref.fa.meth.ann
ref.fa.meth.bwt.2bit.64
ref.fa.meth.pac
```

> **Note — neither index ships a `.0123`**
>
> By default, neither the original index nor the `.meth` seed index writes an
> unpacked `.0123`. `mem --meth` seeds against the seed FM-index but scores/extends
> against the **original** reference, whose bases it pac-fetches from
> `ref.fa.pac` — so the original `.0123` (~6.4 GB) is unnecessary, and the seed's
> unpacked bases are never read at all (~13 GB). Skipping both saves ~19 GB of
> disk on hg38, while the runtime RSS reduction comes from avoiding the original
> `.0123` load. (`--emit-unpacked-ref`
> applies to the original index only, for bwa-mem2 compatibility; the seed never
> needs one.)

The `.meth` seed FM-index is roughly twice the size of the normal FM-index (its
contigs are doubled), so a `--meth` build is larger than a plain build but well
under 3× (by default no `.0123` is written; `--emit-unpacked-ref` adds it for the
original index only). For hg38, budget on the order of 35 GB
of disk for the combined dual index plus the intermediate `ref.fa.meth.fa`.

> **Tip — Pass the original FASTA to mem, not the seed index**
>
> When running `bwa-mem3 mem --meth`, pass the original FASTA path (`ref.fa`);
> bwa-mem3 finds `ref.fa.meth.*` automatically. A legacy `ref.fa.bwameth.c2t`
> index from an older release is **not** usable — rebuild with
> `bwa-mem3 index --meth` (see
> [Migrating from bwameth.py c2t](../methylation/external-c2t.md)).

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
