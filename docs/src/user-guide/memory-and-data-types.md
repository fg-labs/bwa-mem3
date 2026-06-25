# Memory budgeting and data-type tuning

Most short-read alignments run comfortably with default settings. A handful of
situations push memory hard enough to matter: tight memory caps (containers,
cgroups, schedulers), very high thread counts, and data types with unusually
wide insert-size distributions such as Hi-C. This page explains how peak memory
is built up, how to fit a run inside a fixed cap, and how to align Hi-C and
similar data without wasting memory.

## How peak memory is built up

Peak resident memory during `bwa-mem3 mem` is the sum of two parts:

```text
peak RSS  ≈  resident index  +  per-batch working set
```

**Resident index.** The FM-index, packed reference (`.pac`), and related
structures are loaded once and shared across all threads (there is no per-thread
copy). For hg38 the resident index baseline is roughly **15 GB**. This is fixed
for a given reference and does not change with `-t` or `-K`.

`mem` reconstructs the reference bases it needs for scoring and extension
directly from the packed `.pac` on demand (*pac-fetch*), so the unpacked
`.0123` reference (~6.4 GB on hg38) is **neither loaded nor required on disk** —
the `.pac` already holds the same bases at one-quarter the size. This is the
default and is byte-for-byte identical to loading `.0123`. The escape hatch
`BWAMEM3_REF_PAC_FETCH=0` restores the old behavior of loading `.0123` (only on
an index that still has the file; see [indexing](indexing.md)), and is intended
only for A/B verification. On hg38 (5M read pairs, `-t 16`) pac-fetch lowered
peak RSS by **~6.2 GB** for a plain alignment and **~6.3 GB** for `--meth`, at
neutral-to-slightly-faster wall time.

**Per-batch working set.** On top of the index, each in-flight batch holds the
reads, their seeds, candidate alignment regions, and the reference windows used
by mate rescue. This scales with the *effective batch size* (below) and with the
data: longer reference windows — for example wide mate-rescue windows — make it
larger.

### Methylation (`--meth`) mode

`--meth` loads a **dual index**: the doubled seed FM-index for seeding plus the
**original** reference's packed `.pac` for scoring/extension. The seed FM-index
is roughly twice the size of a plain FM-index (its contigs are doubled), so the
resident index for hg38 is on the order of **22 GB** (seed FM-index ~21 GB +
original `.pac` ~1 GB), versus ~17 GB for a plain alignment.

As with plain alignment, the original reference's bases are pac-fetched from its
`.pac` — the original unpacked `.0123` (~6.4 GB) is **neither built nor loaded**.
The seed index's own unpacked `.0123` (~13 GB) and packed `.pac` (~1.6 GB) are
likewise **neither built nor loaded** — `mem --meth` extends against the original
reference, never the seed, so the seed's bases are never read. (Earlier `--meth`
builds loaded the original `.0123` plus both seed files, costing ~20 GB of extra
RSS in total.) Under `bwa-mem3 shm --meth` the staged seed segment is likewise
seed-only, omitting the seed PAC and `.0123`. The per-batch levers below (`-K`,
`-t`) apply unchanged.

## Effective batch size: `-K` and the `-t` multiplier

The single most important knob for per-batch memory is the effective batch size,
and it is easy to misjudge because of how it interacts with `-t`:

```text
no -K     effective batch = chunk_size × n_threads   (chunk_size default 10,000,000)
-K INT    effective batch = INT                       (fixed, regardless of n_threads)
```

So with the default `chunk_size` and `-t 16`, each batch is
`10,000,000 × 16 = 160,000,000` bases — **not** 10 M. Setting `-K 1000000` makes
the batch a fixed 1 M bases regardless of thread count: in that example a **160×**
reduction in the per-batch working set, not 10×.

Two consequences worth remembering:

- Raising `-t` raises peak memory by default, because the batch grows with the
  thread count. If you scale threads up on a memory-constrained host, scale `-K`
  down to compensate.
- `-K` also makes output deterministic across different `-t` values (see
  [`-K` in the mem reference](../cli/mem.md)), which is why it is useful for
  reproducibility independent of memory.

## Fitting a run inside a memory cap

Under a hard cap (a `--memory` cgroup limit, a Slurm `--mem`, a container limit),
budget like this:

```text
per-batch budget  ≈  memory cap  −  resident index  −  headroom
```

For hg38 under a 22 GB cap, that is roughly `22 − 15 − a couple GB ≈ a few GB`
for the per-batch working set — not much. The levers, in order of preference:

1. **Lower `-K`.** `-K 1000000` keeps the per-batch working set well under 1 GB
   for typical short reads and costs little throughput. This is the first thing
   to try.
2. **Lower `-t`.** Fewer threads shrink the default batch (since it is
   `chunk_size × n_threads`) and reduce concurrency-related overhead. Prefer
   adjusting `-K` first so you keep your cores busy.
3. **Use `bwa-mem3 shm`** if you run several samples on one host: the index is
   mapped from a shared segment and its pages are shared across processes, so the
   ~15 GB index is paid once rather than per process. See
   [Quick start: shared-memory index](../getting-started/quick-shm.md).

> **Tip — a silent OOM looks like a truncated BAM**
>
> When the kernel OOM-kills `bwa-mem3` (for example, a process exceeding its
> cgroup limit), it receives `SIGKILL` with no error message. A downstream
> `samtools` in the same pipe may still exit 0, leaving a **header-only or
> truncated BAM that masquerades as success**. Always run alignment pipes under
> `set -o pipefail` and verify the output BAM is non-empty (e.g. a record count
> or `samtools quickcheck`) before treating a run as complete.

## Hi-C and other wide-insert data

Mate rescue searches a reference window whose width is derived from the observed
insert-size distribution (roughly the inter-quartile range of inferred insert
sizes, scaled out by several multiples). For standard paired-end libraries this
window is a few hundred bases. For **Hi-C** (and capture-C / 3C-style data) the
"insert size" between mates is not a library property at all — mates are ligation
contacts that can sit anywhere in the genome — so the inferred distribution is
enormously wide and the per-attempt rescue window balloons to tens of kilobases.
The result is both wasted compute and a large per-batch working set, which is a
common cause of OOM on Hi-C.

The fix is to align Hi-C the way Hi-C pipelines expect, with `-5SP`:

```sh
bwa-mem3 mem -5SP -t 16 --bam ref.fa hic_R1.fq.gz hic_R2.fq.gz \
  | samtools sort -n -@ 2 -o hic.namesorted.bam -
```

- **`-S`** skips mate rescue — this removes the wide rescue windows entirely and
  is the bulk of the memory and time saving.
- **`-P`** skips pairing, which is meaningless for Hi-C contacts.
- **`-5`** marks the alignment with the smallest coordinate as primary for split
  reads, the convention expected by Hi-C tools (Juicer, Arima, pairtools, etc.).

`-5SP` is the standard Hi-C invocation for the `bwa mem` family; it is the
*correct* mode for Hi-C, not merely a memory workaround. Because it changes which
alignments are reported, apply it consistently across any runs you intend to
compare. Hi-C output is typically name-sorted (`samtools sort -n`) for the
downstream contact caller rather than coordinate-sorted.

For other data with genuinely large but real inserts (e.g. some long-fragment or
mate-pair libraries), mate rescue is still meaningful — keep it, but lower `-K`
to bound the per-batch working set rather than disabling rescue.

---

**See also:**
[Threading and resource use](threading.md) ·
[Aligning short reads (mem)](aligning.md) ·
[`mem` CLI reference](../cli/mem.md) ·
[Quick start: shared-memory index](../getting-started/quick-shm.md) ·
[Best Practices: anti-patterns](../best-practices/anti-patterns.md)
