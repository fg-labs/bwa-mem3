# Anti-Patterns

This page documents common mistakes that produce incorrect results or
unnecessary failures when using bwa-mem3.

## Re-indexing without dropping the shared-memory segment

> **Warning — Footgun**
>
> `bwa-mem3 shm` does not detect stale segments. If you re-run `bwa-mem3 index`
> after a shared-memory segment is already staged, the on-disk index files will
> not match the in-memory segment. `bwa-mem3 mem` will attach to the stale
> segment and produce incorrect alignments without any warning.
>
> Always run `bwa-mem3 shm -d` before re-indexing:
>
> ```bash
> bwa-mem3 shm -d           # drop all staged segments
> bwa-mem3 index ref.fa     # rebuild the on-disk index
> bwa-mem3 shm ref.fa       # re-stage the new index
> ```
>
> There is no automatic staleness check in the implementation. The segment name
> is derived from the reference basename only; no content hash or modification
> timestamp is stored.

To confirm that no stale segments are staged, use `bwa-mem3 shm -l` before
running any indexing step.

## Forgetting to initialize submodules

bwa-mem3 depends on several submodules (`ext/htslib`, `ext/libsais`,
`ext/mimalloc`, `ext/sse2neon`). A shallow clone or a clone
without `--recursive` will produce a build that fails at the linking step with
missing symbols, or at runtime with missing index files.

> **Warning — Missing submodules**
>
> Always clone with `--recursive`, or initialize submodules after cloning:
>
> ```bash
> git clone --recursive https://github.com/fg-labs/bwa-mem3
> # or, after a bare clone:
> git submodule update --init --recursive
> ```
>
> If `make` reports missing headers (e.g. `htslib/hts.h: No such file or
> directory`), the submodules were not initialized.

## Leaving `BASELINE_ARCH` at the default on a known higher-tier CPU

The default `make` (no `arch=`) builds the multi-tier single binary
with non-kernel TUs compiled at `BASELINE_ARCH=avx2`. On a production
server with a known higher-tier CPU family, this leaves auto-vectorized
non-kernel hot paths at 256-bit width when the host could go wider, or
keeps the host-floor precheck at `avx2` when the deployment surface is
strictly AVX-512. Pass `BASELINE_ARCH=` (or build a single-tier binary
with `arch=`) to align the build with the deployment:

Pick a `BASELINE_ARCH=` (or a single-tier `arch=`) build target that matches the
deployment — see [Build → arch targets](build.md#choose-the-right-arch-target).
The default is correct only when the binary is distributed across multiple CPU
families or the target is genuinely unknown. Note that `BASELINE_ARCH=avx512bw`
does not always beat `avx2` even on AVX-512 hosts — see
[`BASELINE_ARCH=avx512bw` build flag](../whats-different/avx512-baseline.md) for
the empirical characterization.

## Mixing bwa-mem3 and bwa-mem2 outputs in the same pipeline

bwa-mem3 adds several custom SAM tags that bwa-mem2 does not emit: `HN:i`
(total number of primary alignments — both reported and suppressed — that the
aligner found for this read, before the `-h` supplementary cap is applied),
and — in `--meth` mode — the Bismark-compatible `XR:Z` (read conversion
direction), `XG:Z` (genome strand), and `XM:Z` (per-base methylation call
string) tags. In `--meth` mode it also builds each `@SQ` line from the original
reference; the doubled `f`/`r` seed contigs never reach the output.

> **Warning — Header and tag mismatch**
>
> Do not merge BAM files produced by bwa-mem3 and bwa-mem2 without verifying
> that the `@PG` headers and custom tags are handled correctly by the downstream
> tool. In methylation workflows, a bwa-mem2 BAM mixed into a bwa-mem3 `--meth`
> pipeline will be missing the `XR:Z` / `XG:Z` / `XM:Z` Bismark annotations,
> which will cause methylation callers to silently drop or misclassify those
> records.

If you must merge outputs from both tools, run `samtools view -H` on both
files and confirm that `@SQ` lines are consistent and that the downstream tool
can tolerate the tag differences.

## Writing compressed BAM to a pipe

Passing `--bam=1` (compressed BAM) when piping to `samtools sort` compresses
the stream on the bwa-mem3 side and then immediately decompresses it on the
samtools side. This wastes CPU on both ends with no benefit.

Use `--bam=0` (uncompressed BAM) for all pipe-to-sort workflows. See
[Output format](output-format.md) for the full explanation and recommended
pipeline.

## Aligning bisulfite / EM-seq data without `--meth`

Symptoms: alignment several times slower than expected, a low mapping rate, and
a high mismapping rate.

Both chemistries — bisulfite treatment and the enzymatic conversion used by
EM-seq — turn nearly every unmodified C into T, so the reads are effectively a
three-letter alphabet. Queried against the normal four-letter index they match
only marginally — enough to generate abundant weak seeds, not enough to place
confidently. Measured on simulated EM-seq (100,000 pairs vs GRCh38, `holodeck`
truth):

| | placement correct | unmapped | confident (MAPQ ≥ 30) mismappings |
|---|---:|---:|---:|
| plain mode | 14.7% | 30.9% | 46,672 |
| `--meth` | **95.4%** | **0.0%** | **41** |

Plain mode is also roughly 4.5× slower: mate-rescue fan-out rises from 2.55 to
34.7 candidate anchors per read end, because the rescue admission bar is
relative to each read's own best score and therefore *falls* as placement
quality falls.

**The slowness is the visible symptom; the wrong answers are the real cost.**
Do not reach for `-m` or other speed flags to recover the runtime — add
`--meth`, which is both faster and correct.

`--meth` requires a methylation-aware index. Build it once per reference with
`bwa-mem3 index --meth ref.fa` — that writes the normal index at the bare prefix
plus the converted seed index under `ref.fa.meth.*`; at alignment time still
pass the original `ref.fa` path. See
[Indexing → Methylation index](../user-guide/indexing.md#methylation-index---meth)
and [Methylation](methylation.md).

## Tuning performance flags when the reads don't match the reference

The same mechanism applies well beyond bisulfite. Mate-rescue cost rises sharply
on any *marginally mappable* input — a wrong or diverged reference build, an
unexpected species, or contamination — because the admission bar for rescue
candidates is relative to each read's best alignment, so degraded placement
admits more work rather than less.

Note that genuinely *unmappable* reads are cheap: reads with no homology fail at
seeding and never reach rescue. It is the in-between case — homologous but
degraded — that is expensive.

So an unexpectedly slow run **with a low mapping rate** is usually a data or
reference problem, not a tuning problem. Check `samtools flagstat` mapping rate
first. Reaching for speed flags in that situation makes a misconfigured run
finish sooner without making it correct.

---

**See also:**
[Output format](output-format.md) ·
[Multi-sample workflows](multi-sample.md) ·
[Build](build.md) ·
[Quick start: shared-memory index](../getting-started/quick-shm.md) ·
[CLI Reference: shm](../cli/shm.md)
