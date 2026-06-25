# Features

This page covers user-facing features added to bwa-mem3 on top of upstream
bwa-mem2. None of these features change default behavior: output produced by
`bwa-mem3 mem` without any of these flags is byte-identical to the
corresponding bwa-mem2 output (except for the `@PG` `ID:` and `PN:` fields
which now read `bwa-mem3`).

## `--meth` bisulfite alignment mode (PR #13)

`--meth` adds native bisulfite/EM-seq alignment to `bwa-mem3 index` and
`bwa-mem3 mem` in a single binary — no Python, no separate post-processing step,
no [bwameth.py](https://github.com/brentp/bwa-meth) dependency. In the default
`--meth-scoring collapsed` mode it reproduces bwameth.py's read *placement* (a
placement drop-in, not a byte-for-byte clone); `--meth-scoring genomic` opts into
variant-aware scoring bwameth cannot produce.

```bash
bwa-mem3 index --meth ref.fa          # once per reference
bwa-mem3 mem --meth ref.fa R1.fq R2.fq | samtools sort -o out.bam
```

`index --meth` builds a dual index: the normal index over the original reference
plus a converted **seed** index `<ref>.meth.*` (over `<ref>.meth.fa`, a doubled
reference with `f`-prefixed C→T and `r`-prefixed G→A contigs). Reads seed in that
3-letter space but are scored against the original 4-letter reference, so the
index layout is **not** the same as bwameth.py's single `.bwameth.c2t` reference.

`mem --meth` projects each read (R1 C→T, R2 G→A) to find seeds in the `.meth`
index, preserving the original bases on the first-class `bseq1_t.meth_orig_seq`
field (a `YS:Z`/`YC:Z` comment carrier is a fallback only; neither reaches the
BAM). It then extends and **scores against the original 4-letter reference** with
a per-strand asymmetric matrix, consolidates the `f`/`r` contig pairs back to one
`@SQ` per real chromosome, emits Bismark-compatible `XR:Z` (read conversion
direction), `XG:Z` (genome strand), and `XM:Z` (per-base methylation call string)
auxiliary tags, restores the original bases into the BAM SEQ field for CpG-calling
tools, optionally applies a chimera QC heuristic (longest M/=/X run < 44% of read
length → set `0x200`, clear proper-pair `0x2`, cap MAPQ at 1) when `--chimera-qc`
is passed (off by default), and writes a `@PG ID:bwa-mem3-meth` entry.

In the default `collapsed` mode this reproduces bwameth.py's read *placement*
(chrom, pos) for the standard case, while scoring against the original reference
rather than in collapsed space — so it is not byte-for-byte identical to bwameth
output. Stacks on PR #12 (`--bam`). See the
[Methylation Reference](../methylation/overview.md) for full details.

## Vendored mimalloc allocator (PR #19)

bwa-mem3 vendors [mimalloc v3.3.0](https://github.com/microsoft/mimalloc) as a
pinned submodule at `ext/mimalloc` and links it into every binary by default
(`USE_MIMALLOC=1`). On Linux, static linkage uses `--whole-archive`; on macOS,
dyld-interposed shared linkage is used.

Measured on AWS c7g.4xlarge (Graviton3, 16 threads, 29M 150 bp paired-end
exome-capture reads vs hg38, page cache dropped between iterations):
**−24.5% wall-clock time** (528.6 s → 424.7 s) compared to the same build
with `USE_MIMALLOC=0`. No user-visible interface change; no runtime
configuration required.

`USE_MIMALLOC=0` is a supported best-effort opt-out and is CI-gated on Linux
x86. `bwa-mem3 version` prints the mimalloc version string when it is active.

## `--supp-rep-hard-cap` supplementary MAPQ rescoring (PR #56)

Supplementary alignments for a split read inherit MAPQ from the full-read
scoring pipeline. Competing repetitive chains for the supplementary fragment
are filtered out during full-read chain scoring (`mem_chain_flt`) before
Smith-Waterman, so they never contribute to `sub`/`sub_n`. A supp fragment
landing in a CCATCC repeat that would map equally well to 50+ locations
standalone can therefore carry MAPQ=60 from its primary.

`--supp-rep-hard-cap INT` opts into rescoring: if any seed in a supplementary
alignment's chain has `>=INT` genome occurrences (from the SMEM SA count), the
supplementary MAPQ is forced to 0. Primary alignment MAPQ and coordinates are
unaffected. Default output (no flag) is byte-identical to upstream bwa-mem2.

The SMEM SA-occurrence count is preserved on each seed as `mem_seed_t.n_hits`
and propagated to `mem_alnreg_t.chain_n_hits` during chain-to-alignment
conversion. Typical values for `INT` are 5–20; lower is more aggressive. The
upstream [bwa-mem2#260](https://github.com/bwa-mem2/bwa-mem2/issues/260)
reporter case drops from MAPQ=60 to MAPQ=0 at `--supp-rep-hard-cap 18`.
Closes [issue #46](https://github.com/fg-labs/bwa-mem3/issues/46).

## Shared-memory index: `bwa-mem3 shm` (PR #65)

`bwa-mem3 mem` reloads the FM-index from disk on every invocation. For hg38
the index is ~18 GB; for short alignment jobs (targeted panels, small sample
batches) this load cost dominates runtime and makes per-invocation IOPS the
bottleneck.

PR #65 ports the `bwa shm` command from bwa-mem v1 to bwa-mem3 with strict v1
CLI parity:

```bash
bwa-mem3 shm <index-prefix>    # load index into shared-memory segment once
bwa-mem3 mem <index-prefix> ...  # subsequent runs attach instead of re-reading
bwa-mem3 shm -d <index-prefix>  # detach and free the segment
```

The index lives in a POSIX shared-memory segment. Multiple `bwa-mem3 mem`
processes on the same host share the same in-memory copy. Closes
[issue #64](https://github.com/fg-labs/bwa-mem3/issues/64).

> **Warning — Stale index**
>
> `bwa-mem3 shm` does not detect when the on-disk index has been rebuilt. Always
> run `bwa-mem3 shm -d <prefix>` before running `bwa-mem3 index` and then
> re-stage with `bwa-mem3 shm <prefix>`. Using a stale shared-memory segment
> produces silently wrong alignments.

## `bwa-mem3 shm --meth` (PR #67)

`bwa-mem3 mem --meth <prefix>` locates the `.meth` seed index built by
`bwa-mem3 index --meth <prefix>` automatically. Before PR #67, staging a
methylation index in shared memory required passing the full suffixed seed-index
path to `shm` while continuing to pass the plain prefix to `mem`. The mismatch was
easy to forget, and the failure mode — a run that silently attached the wrong
segment — was difficult to diagnose.

PR #67 adds `--meth` support to `bwa-mem3 shm` so the same plain-prefix
convention works end-to-end:

```bash
bwa-mem3 shm --meth ref.fa       # stages ref.fa.meth.*
bwa-mem3 mem --meth ref.fa ...   # attaches automatically
bwa-mem3 shm -d --meth ref.fa   # detaches
```

## `HN:i` hit count tag (PR #42)

Every primary SAM/BAM record now carries an `HN:i:<n>` tag reporting the
number of secondary alignment candidates clustered with this primary under
`XA_drop_ratio`. This count is captured before the `-h`/`max_XA_hits` cap
truncates the `XA:Z:` string, so `HN` reports the true number of alternate
loci even when no `XA:Z:` field appears in the record.

This makes it possible to distinguish:

- `HN:i:0` + no `XA:Z:` — genuinely unique mapper.
- `HN:i:N` + `XA:Z:...` (N ≤ `-h`) — multi-mapper with all alternates listed.
- `HN:i:N` + no `XA:Z:` (N > `-h`) — multi-mapper whose alternates were
  suppressed by the cap.

Motivated by [lh3/bwa#438](https://github.com/lh3/bwa/pull/438), which adds
HN to `bwa aln`. HN is emitted in both SAM (`mem_aln2sam`) and BAM
(`mem_aln_to_bam`) paths and is absent when `-a` (`MEM_F_ALL`) is active.

## `--bam=LEVEL` direct BAM output (PR #12)

`bwa-mem3 mem --bam` (or `--bam=0` through `--bam=9`) emits BAM directly via
htslib, bypassing the SAM-text-to-BAM conversion round trip that normally
occurs when the output is piped to `samtools view -bS`.

- `--bam` / `--bam=0`: uncompressed BAM (BGZF framing only) — near-zero CPU
  overhead, smaller than SAM text, fast downstream parsing.
- `--bam=1..9`: BGZF deflate at the specified level.
- No flag: SAM text on stdout (default, unchanged).

The implementation adds `src/bam_writer.{h,cpp}`, a new module that converts
`mem_aln_t` to `bam1_t` via `mem_aln_to_bam`. htslib v1.21 is pulled in as a
submodule at `ext/htslib`. On the bwameth.py example fixture (92,961 records),
`samtools view` of `--bam` output vs SAM text produces a zero-line diff across
all 11 SAM columns and all aux tags. See
[Best Practices → Output format](../best-practices/output-format.md) for the
recommended pipeline.

## `--min-ext-len` short-seed extension filter

`--min-ext-len INT` opts into skipping banded Smith-Waterman extension of seeds
shorter than `INT` bp. Off by default (`0`) → output byte-identical to baseline.

Smith-Waterman extension is ~60 % of `bwa-mem3 mem` CPU, and almost all of it is
spent on short seeds: seeds ≤40 bp hold roughly **90 % of all banded-SW cells**
yet are ~99 % wasted, because long seeds already resolve via the ungapped
fast-path at near-zero cost. `--min-ext-len` drops the short seeds before
extension (`mem_chain_drop_short_seeds`, a stable in-place compaction of each
chain's seeds, called from `mem_flt_chained_seeds` in `src/bwamem.cpp`), so their
extension never runs while seeding and chaining are untouched.

Measured single-thread on hg38:

- **Real HG002 1M PE WGS: ~20 % lower alignment CPU at `--min-ext-len 30`.** The
  accuracy effect is confined to reads that were already low-confidence — 0.40 %
  of reads change (mostly partial, heavily soft-clipped, or repeat/multi-mapping
  reads); **zero** confidently and uniquely mapped reads (MAPQ ≥ 60, NM ≤ 1,
  no soft-clip) regress.
- Simulated clean Illumina and indel-rich data: free (F1 unchanged or slightly
  better, even at larger thresholds).
- The only workload with a real accuracy cost is very high per-base error
  (degraded chemistry, cross-species) — there, every seed is short, so some
  correct alignments depend solely on short seeds. Indels and structural variants
  are *not* a contraindication: an indel leaves two still-long exact segments
  that the fast-path handles.
- Confirmed cross-architecture: Graviton4/Linux reproduces the win (realistic
  +15.8 %, HG002 +20.5 %), matching macOS — the speedup is algorithmic, not
  microarch tuning.

`30` is the recommended opt-in value (accuracy-safe across clean, indel-rich, and
moderately-divergent data); `40`–`50` are also safe on known-clean data but cost
accuracy as divergence rises. Because the speedup is thinning the extension stage
— not faster seeding — the wall-clock gain is smaller than the cell-count
reduction would suggest (seeding, which the filter does not touch, dominates
runtime). See
[CLI → mem `--min-ext-len`](../cli/mem.md#--min-ext-len-int--skip-smith-waterman-extension-of-short-seeds)
and
[Settings profiles](../best-practices/settings-profiles.md#why-we-recommend---min-ext-len-30-standard-error-reads)
for the recommended operating point.

---

## Changes catalog

| Item | bwa-mem3 PR | Upstream PR/issue | Status |
|------|-------------|-------------------|--------|
| `--meth` bisulfite alignment mode | [#13](https://github.com/fg-labs/bwa-mem3/pull/13) | — | fork-only |
| Vendored mimalloc allocator | [#19](https://github.com/fg-labs/bwa-mem3/pull/19) | — | fork-only |
| `--supp-rep-hard-cap` MAPQ rescoring | [#56](https://github.com/fg-labs/bwa-mem3/pull/56) | [bwa-mem2#260](https://github.com/bwa-mem2/bwa-mem2/issues/260) | fork-only (upstream issue open) |
| `bwa-mem3 shm` shared-memory index | [#65](https://github.com/fg-labs/bwa-mem3/pull/65) | — | fork-only |
| `shm --meth` symmetry | [#67](https://github.com/fg-labs/bwa-mem3/pull/67) | — | fork-only |
| `HN:i` hit count tag | [#42](https://github.com/fg-labs/bwa-mem3/pull/42) | [lh3/bwa#438](https://github.com/lh3/bwa/pull/438) | fork-only (analogous to bwa aln) |
| `--bam=LEVEL` direct BAM output | [#12](https://github.com/fg-labs/bwa-mem3/pull/12) | — | fork-only |
| `--min-ext-len` short-seed extension filter | _pending_ | — | fork-only (opt-in, off by default) |

---

**See also:**
[Methylation Reference → Overview](../methylation/overview.md) ·
[User Guide → Memory allocator](../user-guide/allocator.md) ·
[User Guide → Output: SAM/BAM, headers, tags](../user-guide/output.md) ·
[Getting Started → Quick start: shared-memory index](../getting-started/quick-shm.md) ·
[Best Practices → Output format](../best-practices/output-format.md)
