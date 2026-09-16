# re-sa

`bwa-mem3 re-sa` resamples the **suffix-array (SA) sample table** of an existing
index to a new sample rate, rewriting `<idxbase>.bwt.2bit.64` in place. It lets
you change the space/speed trade-off of an index you already built — **without**
re-running the expensive `index` build from the FASTA. Run it with no `-u` to
just report the index's current SA rate. `re-sa` is a bwa-mem3-specific command
with no bwa-mem2 counterpart, added in
[#510](https://github.com/fg-labs/bwa-mem3/pull/510).

SA resolution turns a BWT row into a reference coordinate by walking back
(LF-mapping) to the nearest stored SA sample and adding the step count. One row
in every `1 << u` is sampled (`u` is the sample-rate *shift*; the default is
`u = 3`, i.e. one row in 8). A **denser** table (smaller `u`) means fewer LF
steps per resolve — faster alignment — at the cost of a larger index and more
resident memory. A **coarser** table (larger `u`) shrinks the index and its
memory footprint, at the cost of more LF steps per resolve.

This is the on-disk counterpart of [`shm -u`](shm.md#-u-int--densify-the-staged-sa-sample-table):
where `shm -u` synthesises a denser table into shared memory once per stage,
`re-sa` writes it to disk once, so every later `mem` that loads from disk — and
every fresh `shm` stage — picks up the new rate with no extra flags. An index
already staged in shared memory keeps serving the *old* rate until it is dropped
and restaged (see the warning below).

## Synopsis

```text
{{#include ../../_generated/cli/re-sa.txt}}
```

## Common usage

Report an index's current SA sample rate (no `-u`):

```bash
bwa-mem3 re-sa ref.fa
# [re-sa] ref.fa.bwt.2bit.64: SA sample rate 1/8 (shift 3), ref_seq_len ..., ... samples
```

Make an existing (default, stride-8) index denser for faster resolves:

```bash
bwa-mem3 re-sa -u 2 -t 16 ref.fa    # rewrite ref.fa.bwt.2bit.64 as a stride-4 table
bwa-mem3 mem -t 16 ref.fa R1.fq R2.fq > out.sam   # picks up the denser rate automatically
```

Shrink an over-dense index back down:

```bash
bwa-mem3 re-sa -u 4 ref.fa          # rewrite as a stride-16 table
```

Keep the original and produce a differently-rated copy — copy the whole index
first, then resample the copy (see the gotcha below on sibling files):

```bash
mkdir -p /idx/dense
cp -a ref.fa.* /idx/dense/             # all index siblings, incl. any .alt/.0123
cp ref.fa  /idx/dense/                 # if you keep the FASTA alongside
bwa-mem3 re-sa -u 1 -t 16 /idx/dense/ref.fa
```

## Flag reference

### `-u INT` — target SA sample-rate shift

Sets the new sample rate to `1 / (1 << INT)`, i.e. one stored sample per
`1 << INT` BWT rows. Must be in `[0, 6]` (the range the on-disk format supports:
the sample period `1 << INT` must divide the checkpoint block size of 64).
**Omit `-u` entirely to run in read-only inspection mode**, which prints the
index's current SA rate, `ref_seq_len`, and sample count, then exits without
touching the file.

- `INT` **smaller** than the current on-disk shift **densifies** the table. Each
  newly added sample is recovered by the same LF-walk the resolver uses, so the
  result is byte-identical to an index built directly at that rate with
  [`index -u INT`](index-cmd.md#-u-int--sa-sample-rate). This is a deterministic
  property of the construction, not a benchmarked measurement — it holds
  independent of host, architecture, or SIMD tier, and is checked against
  `index -u` across sample rates by `test/resa_byte_identity_test.sh`.
- `INT` **larger** than the current shift **coarsens** the table by decimating
  the stored samples (every target-sampled row is already a stored row, so no
  walk is needed).
- `INT` **equal** to the current on-disk rate is a no-op: the command reports
  that nothing needs doing and leaves the file untouched.

### `-t, --threads INT` — worker threads

Number of worker threads for the densify pass, default `1`. The per-sample
LF-walks are independent and write disjoint slots, so densification scales
near-linearly with thread count (on hg38, Graviton4, arm64/NEON tier, the
densify pass runs ~11 min at `-t 1` and ~40 s at `-t 16`). It has **no effect
when coarsening** (a memory-bandwidth-bound decimation, run serially) or **when
inspecting** (no `-u`); passing `-t` in those cases warns and is ignored.

## Performance

A denser on-disk table holds the same SA samples that
[`shm -u`](shm.md#-u-int--densify-the-staged-sa-sample-table) synthesises at the
same rate, so `mem` does the identical reduced resolve work — the same speed-up,
by construction, but permanent and with no per-run staging cost. The magnitude
depends on the workload; representative figures on hg38, a 5M-read-pair WGS
slice (Graviton4, arm64/NEON tier, 16 threads), as relative user-CPU deltas
versus a plain stride-8 index (indicative, not absolute):

| `-u` | stride | user-CPU Δ vs stride-8 | added index / resident size (hg38) |
|---|---|---|---|
| `2` | 4 | **≈ −2.5 %** | +4 GB |
| `1` | 2 | **≈ −4.2 %** | +12 GB |

**`-u 2` (stride-4) is the best speed-per-GB operating point and the recommended
target.** `-u 1` (stride-2) is a real further step (another ≈ −1.7 % for +8 GB
more) where the disk and RAM budget allows. The cost is paid on disk (a larger
`.bwt.2bit.64`) and in resident memory each time `mem` loads the index; the
speed-up scales with how resolve-heavy the workload is and how large the
reference is. Coarsening trades in the other direction — a smaller index and
lower memory for slightly more CPU per resolve.

## Notes / Gotchas

> **Note — output is byte-identical to `index -u`**
>
> A resampled index is bit-for-bit identical to one built directly at the target
> rate: densification recovers exactly the SA value the resolver would compute
> for each added row, and coarsening keeps a subset of the existing samples
> verbatim. Alignments produced against a resampled index therefore match those
> from the original, and from a freshly built index at the same rate. This is a
> deterministic property of the construction, not a benchmarked measurement — it
> holds independent of host, architecture, or SIMD tier, and is checked across
> sample rates by the byte-identity + alignment-equivalence test
> `test/resa_byte_identity_test.sh`.
>
> **Note — only `.bwt.2bit.64` is rewritten; sibling files are shared**
>
> `re-sa` touches only `<idxbase>.bwt.2bit.64`. The BWT/checkpoint-occ block
> inside it is copied through unchanged (it does not depend on the SA rate), and
> the sibling index files (`.pac`, `.ann`, `.amb`, and any `.alt`/`.0123`) are
> left alone. Because those siblings are required for alignment (the `.alt`, when
> present, drives ALT-aware primary/supplementary selection and MAPQ), `re-sa`
> rewrites the index **in place** rather than to a new prefix — a new prefix
> would orphan them. To keep a copy at a different rate, copy the whole index
> first, then `re-sa` the copy.
>
> **Note — the write is atomic**
>
> The new table is written to a **unique** same-directory temporary sibling
> (created with `mkstemp`, an `<idxbase>.bwt.2bit.64.XXXXXX` name) and
> `rename(2)`d over the original on success, so an interrupted run never leaves a
> partial or corrupt index at the canonical path — a concurrent reader sees
> either the old index or the new one. The unique name also means two overlapping
> `re-sa` runs on the same index never share a temporary inode.
>
> **Warning — drop any staged shared-memory segment first**
>
> Shared memory has no staleness check. If the index is currently staged
> (`bwa-mem3 shm <idxbase>`), `mem` keeps attaching to that segment, which still
> serves the **old** SA sample rate — so after `re-sa` you silently keep the
> previous rate and its performance characteristics (not the new rate) until you
> re-stage. Alignments stay correct either way, since `re-sa` changes only the
> sample rate and that is alignment-invariant; it is the speed that fails to
> update. `re-sa` detects a staged index and prints a warning when it rewrites
> one, but it still proceeds — run `bwa-mem3 shm -d` before re-aligning, then
> re-stage to pick up the new rate.

---

**See also:**
[CLI Reference — index](index-cmd.md) ·
[CLI Reference — shm](shm.md) ·
[Getting Started — Quick start: shared-memory index](../getting-started/quick-shm.md) ·
[User Guide — Indexing the reference](../user-guide/indexing.md) ·
[User Guide — Memory budgeting and data-type tuning](../user-guide/memory-and-data-types.md)
