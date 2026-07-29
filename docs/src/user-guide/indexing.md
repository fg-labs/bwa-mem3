# Indexing the reference

Before aligning reads, bwa-mem3 builds an FM-index from the reference FASTA.
The index is read back from disk at the start of every `mem` run, so it is
built once and reused indefinitely.

## Basic indexing

```bash
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

> **Note — reusing an existing index**
>
> - **From `bwa-mem2`:** an index built by `bwa-mem2 index` works as-is — no
>   rebuild. bwa-mem3 reads its `.bwt.2bit.64` and `.pac` and ignores the
>   `.0123`.
> - **From `bwa` (v1):** a `bwa index` uses a different FM-index format
>   (`.bwt` / `.sa`) and **cannot** be reused — run `bwa-mem3 index ref.fa` to
>   rebuild (a few minutes; the FASTA is unchanged).
>
> See [Coming from bwa or bwa-mem2](../getting-started/migrating.md).

<!-- Separates adjacent blockquotes: a bare blank line between two of them is
     ambiguous (markdownlint MD028) and some renderers fuse them into one. -->

> **Note — ALT-aware references**
>
> If the reference has a `<prefix>.alt` file *and* a `<prefix>.hdr` /
> `<baseprefix>.dict` header sidecar, `index` warns when the sidecar's `@SQ`
> records lack `AH:*` on the ALT contigs — a Picard/GATK `.dict` never has
> them. See [`@SQ` in Output](output.md#sq) for why, and how to fix it.

<!-- Separates adjacent blockquotes -- see above. -->

> **Note — no `.0123` by default**
>
> Earlier releases (and bwa-mem2) also wrote `ref.fa.0123`, an *unpacked*
> forward+reverse reference (~8× the `.pac`; ~6.4 GB on hg38). bwa-mem3 no longer
> builds it: `mem` reconstructs the bases it needs directly from the packed
> `.pac` on demand (*pac-fetch*), so `.0123` is never read. Output is byte-for-byte
> identical. Pass `--emit-unpacked-ref` to `index` if you need the file for an
> external tool that still requires it (e.g. bwa-mem2):
>
> ```bash
> bwa-mem3 index --emit-unpacked-ref ref.fa   # also writes ref.fa.0123
> ```
>
> `mem` ignores any `.0123` present and always pac-fetches from `.pac`.

## Methylation index (`--meth`)

```bash
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

```bash
bwa-mem3 index /data/indexes/hg38/hg38.fa
# writes hg38.fa.bwt.2bit.64, etc. into /data/indexes/hg38/
```

## Time and memory

Index construction is **multi-threaded**. `index -t` defaults to the detected
core count (cgroup-aware), and `--max-memory` defaults to the memory available
to the process less a reserve of `min(max(2 GB, 5%), 50%)` — a batch build is
allowed to use the host it was given. The reserve never exceeds half of total, so
a host below 4 GB still gets a usable budget (a 1 GB host budgets 512 MB). On a
typical multi-core host, indexing hg38 takes a few minutes (longer if pinned to a
single core).

Peak memory is roughly **12 bytes per base of doubled text**, i.e. ~72 GB for
hg38. Below ~1.07 Gbp the suffix array still fits 32-bit entries and the cost
drops to ~8 bytes per base. `--max-memory` is a precondition checked before the
build starts, not a spill target: a reference that does not fit the budget is
rejected up front (exit code 3) rather than built more slowly. Budget hg38 on a
host with ~76 GB of RAM or more; the rejection message names both the
`--max-memory` override and the host size that would work without one.

`--meth` needs about **twice** that, because it builds a second index over a
per-strand-converted reference whose text is twice as long: ~144 GB estimated for
hg38, so budget a host with ~152 GB of RAM or more. Both builds run in one
process and the first does not fully release its memory before the second starts,
so the per-base figures above already account for that overlap.

bwa-mem3 builds the suffix array with
[libsais](https://github.com/IlyaGrebnov/libsais), whose OpenMP-parallel
construction is faster and leaner than the original bwa-mem2 approach. It is not
a bounded-memory construction — see the per-base cost above. See
[Performance improvements](../whats-different/performance.md) for benchmark
numbers.

> **Warning — Do not index over a live shared-memory segment**
>
> If you have previously staged the index into shared memory with `bwa-mem3 shm`,
> drop the segment first before re-indexing:
>
> ```bash
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
