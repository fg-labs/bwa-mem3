# Features

This page covers user-facing features added to bwa-mem3 on top of upstream
bwa-mem2. None of these features change default behavior: output produced by
`bwa-mem3 mem` without any of these flags is byte-identical to the
corresponding bwa-mem2 output (except for the `@PG` `ID:` and `PN:` fields
which now read `bwa-mem3`).

## `--meth` bisulfite alignment mode (PR #13)

`--meth` turns `bwa-mem3 index` and `bwa-mem3 mem` into a single-binary
drop-in replacement for the entire
[bwameth.py](https://github.com/brentp/bwa-meth) pipeline. No Python, no
separate post-processing step, no bwameth.py dependency.

```bash
bwa-mem3 index --meth ref.fa          # once per reference
bwa-mem3 mem --meth ref.fa R1.fq R2.fq | samtools sort -o out.bam
```

`index --meth` writes `<ref>.bwameth.c2t` — a doubled reference with
`f`/`r`-prefixed contigs and C→T / G→A projection, byte-identical to the
index that `bwameth.py index-mem2` produces.

`mem --meth` performs inline C→T conversion of R1 and G→A conversion of R2
before seeding, stashes the original bases in `YS:Z:`, records the conversion
direction in `YC:Z:`, consolidates the `f`/`r` contig pairs back to one `@SQ`
per real chromosome, applies a chimera QC heuristic (longest M/=/X run < 44%
of read length → set `0x200`, clear proper-pair `0x2`, cap MAPQ at 1), copies
`YS:Z:` back into the SEQ field for CpG-calling tools, and writes a `@PG
ID:bwa-mem3-meth` entry.

On the bwameth.py example fixture (92,684 reads), end-to-end output is
byte-identical on chrom, pos, CIGAR, and SEQ vs the bwameth.py oracle. Stacks
on PR #12 (`--bam`). See the
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
the index is ~28 GB; for short alignment jobs (targeted panels, small sample
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

`bwa-mem3 mem --meth <prefix>` auto-appends `.bwameth.c2t` to locate the
methylation index built by `bwa-mem3 index --meth <prefix>`. Before PR #67,
staging a methylation index in shared memory required passing the full
`.bwameth.c2t`-suffixed path to `shm` while continuing to pass the plain
prefix to `mem`. The mismatch was easy to forget, and the failure mode — a
run that silently attached the wrong segment — was difficult to diagnose.

PR #67 adds `--meth` support to `bwa-mem3 shm` so the same plain-prefix
convention works end-to-end:

```bash
bwa-mem3 shm --meth ref.fa       # stages ref.fa.bwameth.c2t
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

---

**See also:**
[Methylation Reference → Overview](../methylation/overview.md) ·
[User Guide → Memory allocator](../user-guide/allocator.md) ·
[User Guide → Output: SAM/BAM, headers, tags](../user-guide/output.md) ·
[Getting Started → Quick start: shared-memory index](../getting-started/quick-shm.md) ·
[Best Practices → Output format](../best-practices/output-format.md)
