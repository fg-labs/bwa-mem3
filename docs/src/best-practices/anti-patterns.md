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

bwa-mem3 depends on several submodules (`ext/htslib`, `ext/safestringlib`,
`ext/libsais`, `ext/mimalloc`, `ext/sse2neon`). A shallow clone or a clone
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

> **Warning — Suboptimal build on known hardware**
>
> ```bash
> # Single multi-tier binary with non-kernel TUs at the host's tier:
> make BASELINE_ARCH=avx512bw      # Cascade Lake / Ice Lake / Sapphire Rapids / Zen 4
>
> # Single-tier binary (no dispatch table; smallest install) when the cluster
> # is uniform and you don't need cross-tier portability:
> make arch=avx2                   # Broadwell/Skylake and later x86
> make arch=avx512bw               # Cascade Lake / Sapphire Rapids
> make arch=arm64                  # Apple Silicon / AWS Graviton
> ```
>
> The default (`make` with no overrides) is appropriate when the binary
> will be distributed across multiple CPU families or when the target
> CPU is genuinely unknown. Note that `BASELINE_ARCH=avx512bw` does not
> always win over `avx2` even on AVX-512 hosts — see
> [`BASELINE_ARCH=avx512bw` build flag](../whats-different/avx512-baseline.md)
> for the empirical perf characterization.

See [SIMD dispatch matrix](../performance/simd-dispatch.md) for the full set
of targets and the in-process dispatch architecture.

## Mixing bwa-mem3 and bwa-mem2 outputs in the same pipeline

bwa-mem3 adds several custom SAM tags that bwa-mem2 does not emit: `HN:i`
(total number of primary alignments — both reported and suppressed — that the
aligner found for this read, before the `-h` supplementary cap is applied),
and — in `--meth` mode — the Bismark-compatible `XR:Z` (read conversion
direction), `XG:Z` (genome strand), and `XM:Z` (per-base methylation call
string) tags. It also rewrites `@SQ` header lines in `--meth` mode
(collapsing `f`/`r` strand prefixes back to one entry per chromosome).

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

---

**See also:**
[Output format](output-format.md) ·
[Multi-sample workflows](multi-sample.md) ·
[Build](build.md) ·
[Quick start: shared-memory index](../getting-started/quick-shm.md) ·
[CLI Reference: shm](../cli/shm.md)
