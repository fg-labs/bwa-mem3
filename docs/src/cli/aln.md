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

**Byte-identity.** The SAM records are byte-identical to the
output of stock bwa 0.7.19 `aln` + `samse`/`sampe` (the header is too, except
the `@PG` line, which names the program that ran). This is checked against a
stock bwa 0.7.19 binary on generated single-end and paired-end reads (a phiX plus
synthetic-contig reference, a few thousand reads per case) with every option below, at
one and at four threads (the input, output and index forms are checked at one thread).
The check runs in CI on Linux x86_64 (AVX2 build), Linux arm64 and macOS arm64, each
at the one SIMD tier the runner's build selects (NEON on arm64); other x86_64 tiers are
not compared with stock here. The search is scalar; the gapped refinement and the mate Smith-Waterman use bwa-mem3's
`ksw_global`/`ksw_align`, whose SIMD tiers are checked against the scalar kernels by
the unit tests rather than by this comparison.

The reference is stock bwa built **without floating-point contraction** (fused
multiply-add), which is what stock's own Makefile gives on x86_64 and what the CI
oracle is built with. A stock bwa built with its default flags on arm64 may fuse the
insert-size and pairing arithmetic of `sampe` and then differ from `aln` (and from an
x86_64 stock build) at rare rounding boundaries.

The index is the one `bwa-mem3 index` writes (bwa-mem3's own format, not a
stock `bwa index`). An index built with `index --meth` works and gives the same
output as a plain index; an index staged with [`shm`](shm.md) is attached
instead of read from disk. `aln` never aligns in methylation mode.

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
| `-t INT` | aln | 1 | Threads |
| `-M`, `-O`, `-E INT` | aln | 3, 11, 4 | Mismatch, gap-open and gap-extension penalties |
| `-R INT` | aln | 30 | Stop searching when there are more than INT equally best hits |
| `-q INT` | aln | 0 | Quality threshold for trimming reads down to 35 bp |
| `-B INT` | aln | 0 | Barcode length (at most 63, or 31 with two read files) |
| `-L`, `-N`, `-I`, `-Y` | aln | off | Log-scaled gap penalty; non-iterative search; Illumina 1.3+ qualities; drop Casava-filtered reads |
| `--max-xa-hits INT` | samse, sampe `-n` | 3 | Maximum alternative hits in the `XA` tag |
| `-r STR` | samse, sampe | none | Read-group header line |
| `-f FILE` | samse, sampe | stdout | Output file |
| `-a INT` | sampe | 500 | Maximum insert size |
| `--max-occ INT` | sampe `-o` | 100000 | Maximum occurrences for one end |
| `--max-xa-hits-discordant INT` | sampe `-N` | 10 | Maximum hits to output for discordant pairs |
| `-c FLOAT` | sampe | 1e-5 | Prior of the chimeric rate (lower bound) |
| `-s`, `-A` | sampe | off | Skip Smith-Waterman for an unmapped mate; skip the insert size estimate |
| `-P` | sampe | off | Accepted and ignored: the index is always in memory |
| `--bam` | new | off | Write BAM instead of SAM |
| `--bam-level INT` | new | 1 | BAM compression level (0-9; other values are an error) |

The long names exist because `-n` and `-o` already mean the search options.
`--max-xa-hits` applies to single-end and paired-end output alike.

The paired-end options are ignored, with a warning, when only one input file is
given.

**Option values.** Numbers must be well formed: where stock's `atoi` reads `abc` as 0 and
`20x` as 20, `aln` stops with an error naming the option. Values that stock accepts but
cannot run are refused with `exit 1` before the index loads: a negative `-l`, `-m`, `-M`,
`-O`, `-E`, `-B` or `-n`, or penalties and `-o` that leave the search no room for one
difference (`-M 0 -O 0 -E 0`, or a `-o` below -1 that cancels the others). On those stock
crashes, aborts or drops every read. Other negative values (`-e -1`, `-o -1`, `-i -1`, and
so on) run as in stock. `-t` below 1 runs one thread, with a warning.

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
- The reference is given as the index prefix, the same path you gave to `index`.
- Each input file is checked, and `-f` is checked for writing, before the index
  is loaded, so a typo costs no index load and leaves no partial output. An `-f`
  that names one of the input files is refused (writing it would truncate the
  reads before they are read).
- Output goes to standard output or `-f`. With `--bam`, `samtools view` of the
  result gives the SAM records.
- A read that contains `-` cannot be written as BAM: `--bam` stops with an error
  naming the read (SAM output is fine, and matches stock).
- With `--bam`, every reference a record names must be in the header. An `-r`
  text with its own `@SQ` lines replaces the index's (as in stock); SAM output
  then keeps the names as text, as stock does, but `--bam` stops with an error at
  the first record whose reference the header lacks instead of writing it as
  unmapped.
- Standard input (`-`) can stand for only one input: `aln ref - -` is refused
  with an error before the index loads (the first `-` would take all of standard
  input). For paired reads from pipes, use two named files or process
  substitution (`<(zcat r1.fq.gz) <(zcat r2.fq.gz)`).
- A malformed FASTQ record (no quality line, or one not as long as the
  sequence) ends the batch and is skipped, as in stock, so the records written
  are stock's; `aln` names the record on standard error and exits 1, where stock
  exits 0. A failure to close an input (for example a truncated gzip stream) is
  reported and exits 1 too, after the output (BAM included) is completed.

## Notes

- **Unequal pair counts.** If the second file has fewer reads than the first,
  `aln` pairs the first reads of each, stops, and prints a warning on standard
  error, as `bwa sampe` does. If the first file has fewer, it pairs only the
  reads that both have (stock `sampe` would read past the end of the first
  file). Check the warning when the counts should be equal: the exit status is
  0 in both cases, as stock `sampe`'s is when the second file is shorter (a
  pipeline that checks only the exit status will not see the dropped reads).
- **Read 1 shorter than read 2.** If the first file runs out first, the extra
  reads of the second are ignored, with a warning, and the exit status is 0.
- **Mismatched pair names.** A pair whose names differ (after the `/1` and `/2`
  trim) ends the run after that pair is written, with an error and exit 1, as
  `bwa sampe` stops.
- **No `.sai` files.** Nothing is written between the search and the SAM stage.
- **Header.** `@HD` (unless the `-r` text has an `@HD` line), `@SQ` (unless
  the `-r` text has an `@SQ` line), the `-r` line, then `@PG`.

---

**See also:**
[mem](mem.md) · [index](index-cmd.md) · [shm](shm.md)
