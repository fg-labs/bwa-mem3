# aln

`bwa-mem3 aln` aligns short reads with bwa's **backtracking** algorithm (the
algorithm behind `bwa aln`), and writes SAM or BAM in a single step. It exists
for reads that BWA-MEM-style seed-and-extend does not serve well: **ultrashort,
ancient and otherwise degraded DNA**, where reads are 20-40 bp and a few
mismatches or deamination changes matter. The usual settings for such data are
`-l 1024 -n 0.04 -o 2`.

For reads of about 70 bp and longer, use [`mem`](mem.md) instead.

## The fused command and the stock commands it replaces

Stock bwa needs three or four commands and writes `.sai` files in between:

```bash
# single end
bwa aln ref.fa r.fq > r.sai
bwa samse ref.fa r.sai r.fq > out.sam

# paired end
bwa aln ref.fa r1.fq > r1.sai
bwa aln ref.fa r2.fq > r2.sai
bwa sampe ref.fa r1.sai r2.sai r1.fq r2.fq > out.sam
```

`bwa-mem3 aln` does the same in one command, with no `.sai` files:

```bash
bwa-mem3 aln ref.fa r.fq > out.sam            # bwa aln + samse
bwa-mem3 aln ref.fa r1.fq r2.fq > out.sam     # bwa aln (twice) + sampe
```

**Byte-identity.** The tested SAM records and headers matched the output of stock
bwa 0.7.19 `aln` + `samse`/`sampe`, except for the `@PG` line, which names the
program that ran. The comparison uses generated single-end and paired-end reads
(a phiX plus synthetic-contig reference, a few thousand reads per case) with every
option below, at one and at four threads (the input, output and index forms are
checked at one thread). It runs in CI on two GitHub-hosted runners, each at the
one SIMD tier its build selects: `ubuntu-24.04` (Linux x86_64, AVX2 build) and
`ubuntu-24.04-arm` (Linux arm64, NEON). Other inputs and x86_64 tiers are not
compared with stock here. The search is scalar; the gapped refinement and the mate Smith-Waterman use bwa-mem3's
`ksw_global`/`ksw_align`, whose SIMD tiers are checked against the scalar kernels by
the unit tests rather than by this comparison.

The reference is stock bwa built **without floating-point contraction** (fused
multiply-add), which is what stock's own Makefile gives on x86_64 and what the CI
oracle is built with. A stock bwa built with its default flags on arm64 may fuse the
insert-size and pairing arithmetic of `sampe` and then differ from `aln` (and from an
x86_64 stock build) at rare rounding boundaries.

The index is the one `bwa-mem3 index` writes (bwa-mem3's own format, not a
stock `bwa index`). For the generated-read workload above, an index built with
`index --meth` gives the same output as a plain index at one thread on the two CI
runners and SIMD tiers listed above. An index staged with [`shm`](shm.md) is
attached instead of read from disk. `aln` never aligns in methylation mode.

## Synopsis

```text
{{#include ../../_generated/cli/aln.txt}}
```

## Options

The options come from three stock commands. The column "Stock" says which one,
and the name it has there when it differs.

| Option | Stock | Default | Meaning |
| --- | --- | --- | --- |
| `-n NUM` | aln | 0.04 | Maximum differences (integer), or the fraction of reads allowed to be missed at a 2% error rate (float) |
| `-o INT` | aln | 1 | Maximum number of gap opens (a fraction is not accepted, despite stock's help text) |
| `-e INT` | aln | -1 | Maximum gap extensions; -1 disables long gaps |
| `-i INT` | aln | 5 | No indel within INT bp of the read ends |
| `-d INT` | aln | 10 | Maximum occurrences for extending a long deletion |
| `-l INT` | aln | 32 | Seed length |
| `-k INT` | aln | 2 | Maximum differences in the seed |
| `-m INT` | aln | 2000000 | Maximum entries in the search queue |
| `-t INT` | aln | 1 | Threads (at most 256; below 1 uses 1) |
| `-M`, `-O`, `-E INT` | aln | 3, 11, 4 | Mismatch, gap-open and gap-extension penalties |
| `-R INT` | aln | 30 | Once there are more than INT equally best hits, keep collecting equally best hits but stop on the first worse-scoring hit |
| `-q INT` | aln | 0 | Quality threshold for trimming reads down to 35 bp |
| `-B INT` | aln | 0 | Barcode length (at most 63, or 31 with two read files); reads with at most `INT` bases are silently discarded |
| `-L`, `-N`, `-I`, `-Y` | aln | off | Widen the no-gap zone at each read end (`-i`) by `log2(gap opens + gap extensions) / 2 + 1` instead of by the plain sum of gap opens and gap extensions already in the alignment, so further gaps can sit nearer the ends (gap penalties are unchanged); non-iterative search; Illumina 1.3+ qualities; drop Casava-filtered reads |
| `--max-xa-hits INT` | samse, sampe `-n` | 3 | Maximum alternative hits in the `XA` tag |
| `-r STR` | samse, sampe | none | Read-group header line |
| `-f FILE` | samse, sampe | stdout | Output file |
| `-a INT` | sampe | 500 | Maximum insert size when no insert-size estimate is available, including with `-A`; an inferred estimate takes precedence |
| `--max-occ INT` | sampe `-o` | 100000 | Maximum occurrences for one end |
| `--max-xa-hits-discordant INT` | sampe `-N` | 10 | Maximum hits to output for discordant pairs |
| `-c FLOAT` | sampe | 1e-5 | Prior of the chimeric rate (lower bound) |
| `-s`, `-A` | sampe | off | Skip Smith-Waterman for an unmapped mate; skip the insert size estimate (which also skips that Smith-Waterman, as in stock) |
| `-P` | sampe | off | Accepted and ignored: the index is always in memory |
| `--bam` | new | off | Write BAM instead of SAM |
| `--bam-level INT` | new | 1 | BAM compression level (0-9; other values are an error) |

The long names exist because `-n` and `-o` already mean the search options.
`--max-xa-hits` applies to single-end and paired-end output alike. Both alternative-hit
limits must be below 2,147,483,647 (`INT_MAX`) to leave room for the extra candidate
used internally.

The paired-end options are ignored, with a warning, when only one input file is
given.

**Option values.** Numbers must be well formed: where stock's `atoi` reads `abc` as 0 and
`20x` as 20, `aln` stops with an error naming the option. Values that stock accepts but
cannot run are refused with `exit 1` before the index loads: a negative `-l`, `-m`, `-M`,
`-O`, `-E`, `-B` or `-n`, or penalties and `-o` that leave the search no room for one
difference (`-M 0 -O 0 -E 0`, or a `-o` below -1 that cancels the others), a `-n`
probability that is 0 as a float, a `-c` that is not finite (`inf`, `nan`), and an `-r`
that contains an empty line or ends in a newline. On those stock crashes, aborts, drops every read or writes a
malformed header. `-n` is read as a probability when it has a `.` or an exponent (`1e-2`);
stock reads `1e-2` as the integer 1 ([#689](https://github.com/fg-labs/bwa-mem3/pull/689)). Other negative values (`-e -1`, `-o -1`, `-i -1`, and
so on) run as in stock. `-t` below 1 runs one thread, with a warning; above 256
is refused before the index loads or output opens.

During search, an actual score outside the priority stack or its 11-bit packed range
(0-2047) exits 1 with a diagnostic instead of overflowing or indexing beyond the stack.
This check depends on the reads and can occur after output has begun. Large penalties
are not rejected just because their theoretical score bound exceeds that range: searches
whose actual scores fit still run.

## Common usage

Ancient DNA, single end, four threads, written to a file:

```bash
bwa-mem3 aln -l 1024 -n 0.04 -o 2 -t 4 -f out.sam ref.fa reads.fq.gz
```

Paired end, with a read group, to BAM:

```bash
bwa-mem3 aln -l 1024 -n 0.04 -o 2 -r '@RG\tID:lib1\tSM:s1' --bam ref.fa r1.fq.gz r2.fq.gz > out.bam
```

Reads from a pipe (single end; `-` is standard input):

```bash
zcat reads.fq.gz | bwa-mem3 aln ref.fa - > out.sam
```

## Inputs and outputs

- Input is FASTQ or FASTA, plain or gzip-compressed. There is no BAM input.
- Reads with 65,536 bases or more remaining after barcode removal (`-B`) cannot be
  searched. The first one is named on stderr and reading stops there. Single-end
  output includes the preceding reads; paired-end output includes only the common
  batch count available from both input files, discarding unmatched records even
  if they precede the oversized read. The run exits 1. (Stock wraps such a length and can report a false hit; its CIGAR also wraps
  from 16,384 bases, which `aln` keeps whole; [#689](https://github.com/fg-labs/bwa-mem3/pull/689).)
- The reference is given as the index prefix, the same path you gave to `index`.
- Each input file and any existing `-f` path are checked before the index is
  loaded. Absent output paths are not opened by this check: failures to create
  them are reported when the output opens, after loading the index. An `-f`
  that names one of the input files, or a file of the index (`.pac`, `.ann`,
  `.bwt.2bit.64` and the others, including the `.meth.*` seed index that
  `index --meth` writes), is refused (writing it would truncate the reads before
  they are read, or overwrite the index).
- Output goes to standard output or `-f`. With `--bam`, `samtools view` of the
  result gives the SAM records.
- A read that contains `-` cannot be written as BAM: `--bam` stops with an error
  naming the read (SAM output is fine, and matches stock).
- With `--bam`, every reference a record names must be in the header. An `-r`
  text with its own `@SQ` lines replaces the index's (as in stock); SAM output
  then keeps the names as text, as stock does, but `--bam` stops with an error at
  the first record whose reference the header lacks instead of writing it as
  unmapped. An `-r` `@SQ` line for an index contig must also give the index's
  length: `--bam` refuses a different `LN` before writing anything, since records
  could otherwise land past the declared end of the reference.
- Standard input (`-`) can stand for only one input: `aln ref - -` is refused
  with an error before the index loads (the first `-` would take all of standard
  input). For paired reads from pipes, use two named files or process
  substitution (`<(zcat r1.fq.gz) <(zcat r2.fq.gz)`).
- A FASTQ record with a `+` separator but no quality line, or a quality string
  not as long as the sequence, ends that file's batch and is skipped, as in stock.
  Without a `+` separator, the reader can interpret the record as FASTA and emit
  it without qualities instead of reporting malformed FASTQ. In single-end
  mode, a nonempty batch is processed and the reader is called again; an empty
  batch stops reading, leaving later reads unprocessed. Recovery need not retain
  every later record: a short quality string can consume following record lines.
  Paired files end batches independently: only the common count is paired,
  surplus reads are discarded, and a later end of file or name mismatch can stop
  processing. `aln` names the malformed record on standard error and exits 1,
  where stock exits 0. A failure to close an input (for example a truncated gzip stream) is
  reported and exits 1 too, after the output (BAM included) is completed
  ([#689](https://github.com/fg-labs/bwa-mem3/pull/689)).

## Memory and parallel jobs

On hg38 (3.2 Gbp with alternate contigs and decoys) `aln` holds about 10.5 GiB of
index: the occurrence table (6.0 GiB), the suffix-array samples at one row in 8
(3.75 GiB) and the packed reference (0.75 GiB). A run adds about 0.4-1 GiB to that
with the default options at `-t 16`, and more on deep searches such as
`-l 1024 -n 0.04 -o 2`.

To run several `aln` jobs on one machine, for example one sample split into
chunks, stage the index once in shared memory with [`shm`](shm.md) (Linux only).
Every job then attaches to the one staged copy instead of loading its own:

```bash
sudo mount -o remount,huge=advise,size=20g /dev/shm   # room for the segment, 2 MB pages
bwa-mem3 shm ref.fa || exit $?                        # launch no jobs if staging fails
bwa-mem3 aln -t 8 ref.fa chunk1.fq.gz > chunk1.sam & pid1=$!
bwa-mem3 aln -t 8 ref.fa chunk2.fq.gz > chunk2.sam & pid2=$!
rc=0
wait "$pid1" || rc=1
wait "$pid2" || rc=1
exit "$rc"                                            # in a script: exit 1 if either job failed
```

A job attaches only to a complete segment. One started while `shm` is still
staging, or after a stage was interrupted, prints `'ref.fa' is not fully staged`
and loads the index from disk, so start the jobs only after `shm` succeeds.

A successful `shm` does not show that the staged index is the one on disk: if a
complete segment for `ref.fa` is already staged, `shm` prints `already in shm` and
exits 0 without checking its contents or replacing it. If the index files have
changed since that segment was staged, remove it and stage again before launching
the jobs.

The example leaves the segment staged, for later jobs to reuse. Remove it, whether
stale or no longer needed, with `bwa-mem3 shm -d`. That removes every staged
segment on the machine, not only this one, so run it only once no other `aln` or
`mem` job relies on a staged index.

Measured on hg38 with two jobs each aligning the first 5M single-end reads of
ERR5343087 with default options and `-t 8`, on 16 vCPU and 32 GB (Sapphire Rapids,
x86_64; clang 19.1.7 LTO build at its default AVX-512BW runtime tier). For these
runs, SAM output matched after removing the `@PG` line:

| | system memory in use, peak | wall time for both jobs |
| --- | --- | --- |
| each job loads its own index | 23.1 GiB | 53.7 s |
| both attach to one `shm` copy | 12.4 GiB (10.5 GiB of it the segment) | 40.3 s |

Each attached job also costs no index load: a single run attached to the segment
took the same time as one loading the index from the page cache (24.0 s against
24.5 s wall). The private pair was slower because two copies did not fit in the
page cache next to the jobs, so index pages were read again from disk. Each job's
own RSS still shows about 10 GiB, because the shared pages it touches count in
every process's RSS; the saving shows in the machine's total (`free`).

The caveats are those of [`shm`](shm.md):

- **Linux only.** macOS limits shared-memory names and segment sizes, and a
  whole-genome index does not stage there.
- **No staleness check.** A segment is found by the index's file name alone, so
  `ref.fa` and `other/ref.fa` find the same one, and nothing checks the files'
  contents. A job attaches to it even if the index files were rebuilt since, and
  then gives wrong alignments. Run `bwa-mem3 shm -d` before rebuilding an index.
- **Size and page size of `/dev/shm`.** The segment must fit in `/dev/shm`
  (`size=`; containers often have a small one), and without `huge=advise` (or
  `huge=within_size`) it sits on 4 KB pages, which costs CPU on every lookup.
  Re-stage after a remount.
- **The segment stays** until `bwa-mem3 shm -d` or a reboot. `shm -d` removes
  every staged segment on the machine, not only this one; jobs already attached
  keep their mapping and finish.

### Smaller index: sparser suffix array (opt-in)

To save memory per process at a small CPU cost, `aln` can use an index whose
suffix-array samples are one row in 16 instead of the default one in 8. Make a
copy of the index and resample the copy with [`re-sa`](re-sa.md), which rewrites
`<idxbase>.bwt.2bit.64` in place:

```bash
mkdir -p aln-idx || exit $?
cp -L ref.fa.* aln-idx/ || exit $?               # copy every index file, dereferencing symlinks
bwa-mem3 re-sa -u 4 aln-idx/ref.fa || exit $?    # suffix-array samples 1 in 16
bwa-mem3 aln -t 16 aln-idx/ref.fa reads.fq.gz > out.sam
```

The rate is stored in the index file, so the index describes itself: `aln` (and
`mem`, `shm`) read it at load, no option needed, and `bwa-mem3 re-sa aln-idx/ref.fa`
reports it. The one-in-8 and one-in-16 indexes were compared on hg38 at `-t 16`
on 16 vCPU and 128 GB (Sapphire Rapids, x86_64; clang 19.1.7 LTO build at its
default AVX-512BW runtime tier), with five runs per index for each workload:

- The first 5M single-end reads of ERR5343087, with default options.
- The paired-end reads from ERR13475328 and ERR13475329 combined (1,362,990 pairs),
  with default options.
- The same pairs with `-l 1024 -n 0.04 -o 2`.

On these workloads, SAM output matched after removing the `@PG` line. Median peak
memory was 1.87 GiB lower for each workload. Median paired-run CPU increases were
0.2% for single-end reads, 1.3% for default paired-end reads and 0.5% for the deeper
paired-end search; individual runs varied. Keep the copy for `aln` only. These
measurements do not cover `mem`.

The copy has the same file name as the original, so while `ref.fa` is staged with
`shm`, `aln aln-idx/ref.fa` attaches to that staged one-in-8 segment instead. To
share the sparser copy, coordinate with other users of shared memory, run
`bwa-mem3 shm -d` to remove every staged segment, then stage the copy with
`bwa-mem3 shm aln-idx/ref.fa`. Start alignment jobs only after staging succeeds.

## Notes

- **Unequal pair counts.** If the second file has fewer reads than the first,
  `aln` pairs the first reads of each, stops, and prints a warning on standard
  error, as `bwa sampe` does. If the first file has fewer, it pairs only the
  reads that both have (stock `sampe` would read past the end of the first
  file). Check the warning when the counts should be equal: the exit status is
  0 in both cases, as stock `sampe`'s is when the second file is shorter (a
  pipeline that checks only the exit status will not see the dropped reads;
  [#689](https://github.com/fg-labs/bwa-mem3/pull/689)).
- **Read 1 shorter than read 2.** If the first file runs out first, the extra
  reads of the second are ignored, with a warning, and the exit status is 0
  ([#689](https://github.com/fg-labs/bwa-mem3/pull/689)).
- **Mismatched pair names.** A pair whose names differ (after the `/1` and `/2`
  trim) ends the run after that pair is written, with an error and exit 1, as
  `bwa sampe` stops.
- **No `.sai` files.** Nothing is written between the search and the SAM stage.
- **Header.** `@HD` (unless the `-r` text has an `@HD` line), `@SQ` (unless
  the `-r` text has an `@SQ` line), the `-r` line, then `@PG`.

---

**See also:**
[mem](mem.md) · [index](index-cmd.md) · [shm](shm.md)
