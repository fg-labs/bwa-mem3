# Optimization checklist

The items below are ordered by expected impact for most workloads. Work through them in sequence; there is little point optimizing output format before confirming you are running the right binary for your CPU.

## 1. Confirm the resolved SIMD tier matches your CPU

The default `make` produces a single binary that contains every supported
x86 SIMD tier and selects one in process at startup. Verify which tier
is running:

```bash
bwa-mem3 version
# expect: SIMD floor: <build_floor>; SIMD runtime: <resolved_tier>
```

If the runtime tier is below what your CPU supports, double-check
whether you accidentally built with a lower `BASELINE_ARCH=` or set
`BWAMEM3_FORCE_TIER` in the environment. Set `BWAMEM3_DEBUG_SIMD=1` to
get a startup banner on stderr at the start of a `mem` run.

On ARM / Apple Silicon, the binary has one NEON tier; `bwa-mem3 version`
reports `SIMD runtime: neon`.

See [SIMD dispatch matrix](../performance/simd-dispatch.md) for the full dispatch
logic and the minimum CPU requirements for each tier.

> **Tip — Single-arch deployments**
>
> On a cluster where every node has the same CPU, build with
> `make arch=avx2` (or the appropriate ISA). The runtime dispatch
> overhead is negligible, but a single-arch build trims the binary
> and removes any chance of `BWAMEM3_FORCE_TIER` accidentally
> downgrading throughput in production.

## 2. Build with PGO if you will run repeatedly

For production pipeline nodes that will process many samples against the same reference, a PGO build provides an additional 2–5% throughput at the cost of one extra build pass and a training run:

```bash
make pgo-generate PGO_ARCH=avx2
./bwa-mem3.pgo-instr.avx2 mem -t 16 ref.fa R1.fq.gz R2.fq.gz > /dev/null
make pgo-use PGO_ARCH=avx2
```

See [PGO build](../performance/pgo.md) for the full workflow, including multi-arch and profile portability notes.

## 3. Use shared memory for many small samples

When aligning many samples on one machine against the same reference, loading the index into POSIX shared memory once and reusing it across all `mem` invocations eliminates redundant I/O and reduces per-sample startup time significantly. The benefit grows with the number of samples and the size of the reference.

```bash
# Load the index into shared memory once
bwa-mem3 shm ref.fa

# Align each sample against the in-memory index
bwa-mem3 mem -t 16 ref.fa R1.fq.gz R2.fq.gz | samtools sort -@ 4 -o sample.bam -

# When finished with all samples, drop the shared segment
bwa-mem3 shm -d
```

> **Warning — No staleness check**
>
> `bwa-mem3 shm` does not detect whether the on-disk index has changed after the segment was loaded. Always run `bwa-mem3 shm -d` before re-indexing a reference and re-loading with `bwa-mem3 shm`. Failing to do so results in alignments against a stale index.

See [Getting Started — Shared-memory index](../getting-started/quick-shm.md) and [Best Practices — Multi-sample workflows](multi-sample.md) for complete workflows.

## 4. Emit BAM directly

Use `--bam` (or `--bam=0` for uncompressed BAM) to emit BAM instead of SAM. Uncompressed BAM avoids the text-formatting cost on the aligner side and the text-parsing cost on the downstream side. `samtools sort` reads BAM natively and is fastest when the input is uncompressed:

```bash
bwa-mem3 mem --bam=0 -t 16 ref.fa R1.fq.gz R2.fq.gz \
  | samtools sort -@ 8 -o out.bam -
samtools index out.bam
```

A compression level — `--bam=6`, say — produces BGZF-compressed BAM, useful when writing directly to disk without a downstream piped tool.

See [Best Practices — Output format](output-format.md) for guidance on when SAM is still appropriate.

## 5. Pipe to a multi-threaded sorter

Sorting is typically the bottleneck after alignment. Keep a separate thread budget for `samtools sort`:

```bash
bwa-mem3 mem --bam=0 -t 12 ref.fa R1.fq.gz R2.fq.gz \
  | samtools sort -@ 8 -m 2G -o out.bam -
```

On a 16-core machine, allocating 12 threads to `mem` and 8 to `samtools sort` (with overlap via the pipe) is a common starting point. The aligner is generally CPU-bound; the sorter is I/O-bound during merge. Profile both stages to find the right split for your hardware.

> **Tip — Thread count tuning**
>
> `bwa-mem3 mem` scales well to 16–32 threads on most workloads. Beyond 32 threads the per-thread work unit becomes small enough that synchronization overhead starts to erode gains. See [User Guide — Threading and resource use](../user-guide/threading.md) for thread-scaling data.

## 6. Reorder seeds longest-first (`--seed-order local-longest`)

After SA-interval resolution bwa-mem3 builds chains by scanning seeds in the order they
were emitted. `--seed-order local-longest` re-sorts each read's resolved seeds by
decreasing length before chaining. The longest seed anchors its chain first and absorbs
shorter seeds that are fully contained within it, so contained sub-seeds are never
extended — they contribute nothing a long seed does not already cover.

```bash
bwa-mem3 mem --seed-order local-longest -t 16 ref.fa R1.fq.gz R2.fq.gz | samtools sort -@ 4 -o out.bam -
```

Measured on 50,000 real WGS reads (1000 Genomes HG00096, hg38), `local-longest` cuts
the fraction of seeds that reach banded Smith-Waterman extension from 38.2 % to 43.7 %
absorbed (an ~8.9 % reduction in extended seeds). Since Smith-Waterman extension is
typically the dominant per-read cost in `bwa-mem3 mem`, the gain scales accordingly with
the seed absorption rate for a given dataset.

**Accuracy and byte-identity.** F1 is flat on an easy simulated profile (holodeck, ~94.4 %;
no regression vs `--seed-order off`). Hard-data F1 validation on divergent/indel-rich
reads and GIAB benchmarks is not yet complete, so `--seed-order local-longest` is
**opt-in only** and the default remains `off`. Non-`off` modes are **not byte-identical**
to `off`: they can shift secondary alignments, `XA:Z:`, `XS:i`, and `HN:i` tags, and a
small number of primaries. See
[Equivalence → Seed ordering](../whats-different/equivalence.md#seed-ordering---seed-order-opt-in)
for the full taxonomy.

Additional advanced modes are accepted by the option (`global-longest`, `absorb-count`,
`most-absorb`) but are unadvertised pending further validation.
## 7. Enable SMEM deduplication (opt-in, not byte-identical)

`--smem-dedup` removes duplicate SMEM seeds before SA expansion, cutting SA
lookups by roughly 10 % on typical short-read WGS data. Off by default to preserve
byte-identical output. Enable only when output identity with upstream bwa-mem2 is
not required:

```bash
bwa-mem3 mem --smem-dedup -t 16 ref.fa R1.fq.gz R2.fq.gz | ...
```

The accuracy impact is confined to a small fraction of reads: on 50 k WGS reads vs
hg38 only 2 reads (0.004 %) changed — one XS tag update on a MAPQ-60 read
(primary placement unchanged) and one tie-break shift on a MAPQ-0 equal-score
locus. All uniquely-mapped reads are unaffected.

> **Not byte-identical**
>
> `--smem-dedup` changes the SMEM set seen by chaining for reads that had
> duplicate SMEMs (arising from B-tree duplicate keys in the SA). This can alter
> XS tags and equal-score placements on 0.004 % of reads (2 of 50 k on the WGS
> validation set). Do not enable in pipelines that compare output to a bwa-mem2
> baseline.

## 8. Reserve 1 GB huge pages for the index (Linux, opt-in)

Seeding is memory-latency-bound. On Linux, `--huge-pages` maps the FM-index
occurrence table and the suffix-array samples on explicit 1 GB huge pages, which
removes nearly all of their data-TLB page walks. Measured on one x86 host (Intel Sapphire Rapids, hg38, 5 M WGS read pairs, `-t 16`): about 1 % lower wall time and 0.3–0.4 % lower CPU time on both the avx512bw and avx2 tiers, with index load no slower on that host. The alignment records
stay byte-identical (page size does not change alignments); only the `@PG` `CL:`
header field differs, since it records the `--huge-pages` option on the command
line. Reserve 1 GB hugepages on the host (hg38 needs 10 at the default `index -u 3`; a lower `-u` needs more), then pass
`--huge-pages`; an array that does not fit in the pool falls back to the default
pages. 1 GB HugeTLB support is host-dependent, so reserve only when the pool
exists.

`nr_hugepages` is a **host-wide total**, not a delta: `echo 12` sets the whole
1 GB pool to 12 pages (it does not add 12), so on a shared host writing an
absolute value below the current pool shrinks it under other processes — read the
current count first and pick a non-decreasing target.

```bash
pool=/sys/kernel/mm/hugepages/hugepages-1048576kB
if [ -d "$pool" ]; then
    cur=$(cat "$pool/nr_hugepages")
    echo $(( cur > 12 ? cur : 12 )) | sudo tee "$pool/nr_hugepages"
fi
bwa-mem3 mem --huge-pages ref.fa r1.fq r2.fq > out.sam
```

This is host configuration (privileged, pinned RAM), so the flag is opt-in, not a
default, and needs Linux 5.14 or later. See
[Memory allocator → Large pages for the index](../user-guide/allocator.md#large-pages-for-the-index-linux-deployment-lever).

## 9. Keep transparent huge pages working (Linux)

bwa-mem3 asks the kernel for transparent huge pages (THP, 2 MB pages) on the
index arrays and on its allocator's arenas, which hold the per-batch buffers.
The kernel can only provide them when there is enough free, unfragmented memory
when those pages are first touched. That is roughly the whole working set, about
28 GiB for a human hg38 index. On a bare Linux host whose RAM is full of page cache
(earlier jobs, other indices, a `cat` prewarm of the index on a host smaller than
about 2.2× the index), part of the working set silently falls back to 4 KB pages.

Measured on a 5 M-pair human WGS slice (HG00096, hg38 stride-2 index, `-t 16`,
clang 19) on AWS Graviton4 (arm64, NEON tier) and AMD Zen 5 (x86-64, AVX-512BW
tier) hosts:

| state | cost |
|---|---|
| per GiB on 4 KB pages | Graviton4: ~0.38 s of a ~37 s run (~1 % wall per GiB); Zen 5: ~1.6 % wall per GiB of a ~23 s run |
| page-cache-full host | 1.2–10 GiB of the working set observed on 4 KB pages, costing the per-GiB rate above |
| THP set to `never` | +27.6 % wall |

To keep THP working:

- Keep `/sys/kernel/mm/transparent_hugepage/enabled` at `madvise` or `always`,
  never `never`.
- Start the run with the working set free. If the host has been busy, drop clean
  page cache first (`sync; echo 1 | sudo tee /proc/sys/vm/drop_caches`), or run the job
  and its I/O in a memory-limited cgroup (for example
  `systemd-run --scope -p MemoryMax=…`) so reclaim of its own page cache stays local.
- `echo 1 > /proc/sys/vm/compact_memory` on its own does **not** help when page
  cache fills RAM. There are too few free 2 MB blocks to assemble.
- Do not launch `mem` with THP disabled for the process (`MIMALLOC_ALLOW_THP=0`,
  or a launcher that calls `prctl(PR_SET_THP_DISABLE)`, which children inherit).

`bwa-mem3 mem` checks this for you on Linux. After output batches 1, 2, 4, 8, …
it walks `/proc/self/smaps` and compares `AnonHugePages` with `Anonymous` over the
mappings marked `THPeligible: 1`, which is memory the kernel would back with 2 MB
pages (THP `always`, or `madvise` and advised). Memory the kernel would not back
with THP (under `madvise`, memory that was never advised) is not counted. It
prints a two-line `[W::main_mem]` warning on stderr, at most once per run, when
the process has ≥ 4 GiB of anonymous memory in total and at least 1 GiB, and at
least 4 %, of the memory judged is not on huge pages. There are three warnings:

- eligible memory fell back to 4 KB pages (the page-cache case above; judged over
  the eligible memory);
- THP is `never` (the global `enabled`, or `hugepages-2048kB/enabled` on Linux ≥
  6.8; judged over all anonymous memory);
- THP is disabled for this process (`THP_enabled: 0` in `/proc/self/status`;
  judged over all anonymous memory).

A walk of `/proc/self/smaps` is not a snapshot, so a sample is up to three walks
(re-walked while the walk looks incomplete against `RssAnon`), plus at most one
retry on the next batch. On the hosts above (Linux 6.1) one walk took 2–4 ms
in-process on a fully-THP ~24 GiB process. Read externally (from another process)
the same file took ~8 ms on a fragmented host and 50–190 ms with ~26 GiB on 4 KB
pages, because the kernel visits every page-table entry. A run samples once per
power-of-two batch and stops at the first warning; on a fully-THP Graviton4 WGS
run the check cost −0.007 % ± 0.031 % wall against a build without it (12
interleaved reps per arm). The warning never changes the output. At `-v 4` each
sample is logged as a `[D::main_mem] THP check batch N:` line, or as
`check disabled: <reason>` when the check stops.

The check is skipped silently:

- at `-v 1`;
- on kernels whose base page is not 4 KB (for example 64 KB-page arm64 kernels);
- when the global `enabled` knob is unreadable, or it or a per-size
  `hugepages-<N>kB/enabled` knob holds an unrecognised mode;
- when any multi-size THP (mTHP) size other than 2 MB is enabled (Linux ≥ 6.8;
  out of scope);
- when `/proc/self/smaps` is unreadable.

An unreadable `/proc/self/status` does not skip it, but turns off the
completeness check against `RssAnon`. An unreadable `hugepages-<N>kB/enabled`
knob is treated as absent (2 MB follows the global mode; other sizes are off).
The fallback warning also needs the `THPeligible` field (Linux ≥ 5.0).

The eligibility accounting is validated for the default build (mimalloc). In a
`USE_MIMALLOC=0` build the index arrays are advised over their page-aligned
interior ([PR #589](https://github.com/fg-labs/bwa-mem3/pull/589)), but under THP
`always` its glibc heap is eligible too, so a fallback warning there does not
necessarily mean page-cache pressure.

To check a run by hand, sum `AnonHugePages` and `Anonymous` over the
`THPeligible: 1` mappings in `/proc/<pid>/smaps`, or diff `thp_fault_fallback` in
`/proc/vmstat` before and after. A `bwa-mem3 shm` index (shared memory) and a
`--huge-pages` reservation are not anonymous memory, so the check does not count
them.

## Summary table

| Item | Action | Reference |
|---|---|---|
| Right SIMD tier for CPU | `bwa-mem3 version`; verify `SIMD runtime:` | [SIMD dispatch matrix](../performance/simd-dispatch.md) |
| PGO for production | `pgo-generate` → train → `pgo-use` | [PGO build](../performance/pgo.md) |
| Shared-memory index | `bwa-mem3 shm ref.fa` before batch runs | [Quick start: shm](../getting-started/quick-shm.md) |
| Emit uncompressed BAM | `--bam=0` | [Best Practices — Output format](output-format.md) |
| Multi-threaded sort | `samtools sort -@` with appropriate thread split | [User Guide — Threading](../user-guide/threading.md) |
| Reorder seeds longest-first | `--seed-order local-longest` | [Equivalence](../whats-different/equivalence.md) |
| SMEM deduplication | `--smem-dedup`; ~10 % fewer SA lookups; opt-in, not byte-identical | [Features → --smem-dedup](../whats-different/features.md#--smem-dedup-smem-deduplication) |
| Reproducible output | `-K INT`; pins the batch size so output does not move with `-t` | [Aligning → `-K`](../user-guide/aligning.md) |
| Pipeline overlap at very high `-t` | `--chunk-cap INT`; off by default, opt-in, not byte-identical | [`mem` → `--chunk-cap`](../cli/mem.md) |
| 1 GB huge pages for the index (Linux) | reserve `nr_hugepages` (host-wide total; 10 pages for hg38 at `index -u 3`), run `mem --huge-pages`; ~1 % lower wall, 0.3–0.4 % lower CPU on a 16-thread WGS slice (Sapphire Rapids, avx512bw and avx2); alignment records byte-identical, only the `@PG` `CL:` field differs | [Memory allocator → Large pages](../user-guide/allocator.md#large-pages-for-the-index-linux-deployment-lever) |
| Transparent huge pages (Linux) | THP `madvise`/`always`; free, unfragmented memory at start (drop page cache or a memory-limited cgroup); heed the two-line `[W::]` THP warning | [Keep transparent huge pages working](#9-keep-transparent-huge-pages-working-linux) |

---

**See also:**
[Performance overview](../performance/overview.md) ·
[SIMD dispatch matrix](../performance/simd-dispatch.md) ·
[PGO build](../performance/pgo.md) ·
[Best Practices — Build](build.md) ·
[User Guide — Threading and resource use](../user-guide/threading.md)
