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

**Byte-identity without `--meth`, `--short-insert` or `--ref-span-insert`.** The tested SAM records and headers matched the output of stock
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
stock `bwa index`). For the generated-read workload above without `--meth`, an
index built with `index --meth` gives the same output as a plain index at one
thread on the two CI runners and SIMD tiers listed above. An index staged with
[`shm`](shm.md) is attached instead of read from disk.

**Methylation (`--meth`).** With `--meth`, `aln` aligns single-end or paired-end bisulfite,
EM-seq, TAPS or 5-base reads against an index built with `index --meth`: it searches the
converted index and writes records against the original reference. Mapped records
carry the tags selected by `--meth-tags`: `XR:Z`, `XG:Z` and `XM:Z` are all selected
by default. Unmapped records carry only `XR:Z` when selected. Stock bwa has no
methylation mode, so these records have no stock counterpart and the byte-identity
above does not cover them. See [Methylation](#methylation---meth) below.

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
| `--short-insert` | new | off | Estimate short inserts without a batch-wide minimum set by the longest read; paired input required |
| `--ref-span-insert` | new | off | Measure inserts from the alignments' reference spans, not read lengths, so indels count; paired input required |
| `-P` | sampe | off | Accepted and ignored: the index is always in memory |
| `--known-snvs FILE` | new | off | Place reads carrying a known SNV's ALT allele as on a genome holding both alleles, and tag reads over a site with `XV:Z` (see [Known SNVs](#known-snvs---known-snvs)) |
| `-v INT` | new | 3 | 0-2 suppress progress, 3 reports progress, 4 adds lookup/window/incomplete diagnostics; some errors and warnings still print at 0 |
| `--bam` | new | off | Write BAM instead of SAM |
| `--bam-level INT` | new | 1 | BAM compression level (0-9; other values are an error) |

The long names exist because `-n` and `-o` already mean the search options.
`--max-xa-hits` applies to single-end and paired-end output alike. Both alternative-hit
limits must be below 2,147,483,647 (`INT_MAX`) to leave room for the extra candidate
used internally.

The stock paired-end options are ignored, with a warning, when only one input file is
given. `--short-insert` and `--ref-span-insert` instead require two input files and reject
single-end use before loading the index or creating output.

### Short, overlapping fragments

Use `--short-insert` when short fragments overlap and a longer read in the batch would
otherwise exclude them from insert-size estimation, or when both reads often cover the whole
fragment. It works with ordinary and
methylation-aware (`--meth`) paired alignment:

```bash
bwa-mem3 aln --short-insert ref.fa r1.fq r2.fq > out.sam
```

The option sets the estimator's lower bound to the lower interquartile fence, with a
minimum of one base, rather than the longest read in the batch. It keeps the upper
fence, minimum number of high-quality pairs and previous-batch fallback. Failed
estimates use the no-estimate fallback; `-A` overrides the option and disables insert
estimation and mate rescue.

The option also accepts a pair whose reverse read starts at the same base as its forward
read. Stock pairing requires the reverse read to start strictly after the forward one, so
when trimmed reads each cover a whole fragment shorter than the read length, almost no pair
is proper and neither mate's MAPQ benefits from the other. Every other geometry requirement
is kept: the reverse read may not start before the forward read, the outer span must be at
least the longer read, and the insert must fit the estimate (or `-a` without one).

Mate rescue (and, with `--meth`, cap rescue) may also search the anchor's own span. Stock
rescue looks for the mate only beyond the anchor, so it can never place a mate that covers the
same fragment. A rescued mate must stay inside the fragment the anchor implies: the forward read
may neither start nor end after the reverse read. Rescue is a local alignment, so a mate that
runs past its anchor (for example, read 1 starting with adapter sequence ahead of the fragment)
is rescued with that overhang soft-clipped.

This is an opt-in departure from stock output. The inferred distribution, the coincident
starts and the wider rescue windows can change pair selection, proper-pair flags, MAPQ, rescued
mates and their CIGARs. This option does not make every short pair proper. Under this option, pairing and methylation cap rescue
reject candidates when the insert-likelihood calculation underflows or its penalty
cannot be represented as an integer. Existing defaults are unchanged.

### Inserts measured from reference spans

Stock `sampe` measures a pair's insert from the reverse read's start plus its *read length*,
because it pairs reads before their gapped alignments are refined. Where that read's
alignment has indels, the insert is off by their net length: with `-A -a 200`, a 50-base
reverse read whose alignment spans 70 reference bases (a 20-base deletion) is marked proper
at a true insert of 220, while one with a 20-base insertion is refused below 200. `--meth`
pairing follows the same rule. `--ref-span-insert` measures every insert from the
alignments' reference spans (read length plus deletions minus insertions) instead:

```bash
bwa-mem3 aln --ref-span-insert ref.fa r1.fq r2.fq > out.sam
```

With the option:

- insert-size estimation takes each observation to the later read's reference end;
- pairing compares the insert to the estimate (or to `-a`) and scores it from the reverse
  hit's reference end, requires the outer span to cover both hits' reference spans, places a
  reverse hit from its own span rather than the read's best hit's, and gives a read that
  pairing moves to another hit that hit's span (stock keeps the first hit's, which the gapped
  refinement then uses);
- mate rescue places its window from the anchor's reference end.

TLEN is unchanged: it is already taken from each record's CIGAR, so it is the reference-span
insert with or without the option. `--meth` cap rescue already measures inserts from spans;
the option only makes its minimum outer span the two alignments' spans, as in pairing. With
`--short-insert`, a coincident-start pair must cover both hits' spans, and rescue's ceiling for
a reverse anchor is its reference end either way. `-A` does not turn the option off: pairing
still measures the insert it compares to `-a` from the spans.

This is an opt-in departure from stock output: it can change the inferred distribution, pair
selection, proper-pair flags, MAPQ, rescue windows and, for reads pairing moves to a gapped
hit, their position and CIGAR. Without the option the output is unchanged.

### Known SNVs (`--known-snvs`)

Short reads that carry the ALT allele of an SNV have one more mismatch against the reference
than reads carrying REF. For reads of 20-30 bases that mismatch often costs the read its
placement: it ties with a paralog, drops below MAPQ 20, or goes over the mismatch budget. Allele
counts at those sites are then biased toward REF. `--known-snvs FILE` places reads as on a genome
that holds both alleles of every SNV in FILE:

```bash
bwa-mem3 aln --known-snvs sites.vcf.gz ref.fa reads.fq.gz > out.sam
```

FILE is a VCF or BCF (plain or compressed). Its single-base substitutions are used, multi-allelic
ones as alternatives (never combined); other records are skipped, and FILTER is ignored. A record
whose contig is not in the reference, or whose REF differs from it, is an error naming the record.

How it works:

- A window of 160 bases either side of each site, with the ALT in place, is indexed when the run
  starts, in memory, with `-t` threads (no files are written).
- Each read is also searched against the windows. A hit there that covers its window's site is
  a hit at the same reference position with the ALT credited, and the cheaper of the two costs
  wins. Then MAPQ, the primary hit and the alternative hits are recomputed with `aln`'s usual
  cuts. The stock primary is kept while it is still among the best.
- Paired end: the insert size is estimated from stock placements first. Mate rescue also aligns
  the mate against its window with each site's ALT substituted, and keeps the better score.
- `--meth`: the converted search also resolves the windows that the ALT-window search finds, with
  the ALTs in place. A pair whose placement assumes an ALT is left out of the insert-size estimate.
- Reads longer than about 160 bases, and reads with more than 256 candidate loci, keep their
  stock result.

CIGAR, MD and NM are still against the reference, so an ALT read shows its mismatch; `XM:i` is the
mismatch count with the ALT credited. A read that does not touch a site is unchanged.

Every record over a site gets `XV:Z`, one entry per site separated by `;`:

```text
contig:pos:REF>ALT,allele,away,max_diff,twin_same,twin_mapq
```

- `allele`: `R` (REF), `A` (ALT), `O` (another base), or under `--meth` `H` (the conversion can
  explain both alleles on that strand).
- `away`: the read's differences outside the site.
- `max_diff`: the read's mismatch budget.
- `twin_same`, `twin_mapq`: for the first `R` or `A` site only (`NA` otherwise), the read is
  realigned with REF and ALT swapped at that site (as a pair, with its mate, for paired input).
  `twin_same` is 1 when that twin lands on the same strand within 5 bases, and `twin_mapq` is
  the twin's MAPQ. The twin never changes the read's own record or the run's random draws.

To count alleles without reference bias, count a read at a site only when it passes the MAPQ
cut-off, its twin also passes (`twin_same` 1 and `twin_mapq` at or above the cut-off), and
`away` is at most `max_diff - 1`. The last rule matters: a read already at its mismatch budget
places only when an error at the site matches the ALT, so without it errors become ALT counts.
See [Counting alleles at known SNVs](../user-guide/known-snvs.md) for worked examples, paired-end
fragments and the limits.

This is an opt-in departure from stock output. It can change the position, MAPQ, flags and
alternative hits of reads at the sites. Without the option the output is unchanged. With a VCF
that has no SNV, the output is also unchanged.

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
  could otherwise land past the declared end of the reference. Under `--meth`,
  user-supplied `@SQ` lines are refused: the header must describe the original
  reference used for alignment and methylation calls.
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

## Methylation (`--meth`)

`aln --meth[=emseq|taps|5base]` aligns converted reads with the backtracking search. It takes the
same `--meth`, `--meth-context`, `--meth-scoring` and `--meth-tags` options as
[`mem --meth`](mem.md), with the same values, defaults and checks, plus `--set-as-failed`, for
single-end and paired-end input. Single-end support was added in
[#749](https://github.com/fg-labs/bwa-mem3/pull/749); paired-end and CpG support in
[#750](https://github.com/fg-labs/bwa-mem3/pull/750). Every chemistry requires both the original and
converted indices written by `bwa-mem3 index --meth`; give the reference prefix
(`ref.fa`) or the converted index (`ref.fa.meth`). Read 1 is projected C->T and
read 2 G->A for directional discovery. In `aln`, `--set-as-failed` requires
`--meth` (MEM accepts the flag without it).

### Discovery, refinement and accounting

Converted-space backtracking proposes original-reference windows within the user's
search limits (`-n`, `-o`, `-e`, seeding options, `-M/-O/-E`, `-N` and `-m`).
Only the directional copy appropriate to the read is retained. Identical windows
(original contig, conversion hypothesis, strand, start and end) are evaluated once.
`-R` does not prune original-reference competitors. For an unseeded search with
current best cost `B`, let `T = B + M` and `k = floor(T / M)`. When `M > 0`,
`T < O`, `E >= 0`, the configured difference ceiling admits `k` differences,
`-N` is off and discovery has not hit a cap, the search can stop after fully
evaluating converted scores through `M*k`. Every retained raw alignment then has
no gaps and at most `k` strict mismatches, so that stratum already proposes its
exact window. A decreasing best cost can only narrow this bound. Other searches
retain configured generation without a path-cost cutoff. This does not make
discovery exhaustive or guarantee a globally optimal alignment.

Each window is globally refined on the original reference with a fixed
chemistry-aware reward model: match +1, mismatch -3, ambiguity -1 and a gap of
length `k` costing `5+k`. Neutral conversions score 0. **`-M/-O/-E` do not change
this refinement model.** The five-symbol model uses an actual-matrix band
certificate when available and otherwise a full band. The existing dispatcher
may use its SIMD global kernels when their predicates hold. CpG-only scoring
uses seven reference symbols and full-band scalar refinement, not that certificate
or SIMD path. In either alphabet, a window as long as the read whose ungapped reward
exceeds `length - 13` skips the DP: a gapped path needs an insertion and a deletion
(at least 12) and aligns at most `length - 1` columns, so the ungapped alignment is
the unique optimum and the result is the one the DP returns.
Original `.amb` intervals override random packed bases with N at
alignment time, without rebuilding either index. CpG labels use true contig flanks
and respect these ambiguity runs.
`--meth-context cpg` requires TAPS/5base and non-collapsed scoring.

The raw refined CIGAR is then costed using `-M/-O/-E`:

```text
cost = M * mismatches + O * gap_runs + E * gap_extensions
```

A gap run of length `k` has `k-1` extensions. Chemistry-accepted conversions,
including neutral conversions, have mismatch cost 0. An ambiguous original-reference
or query base always counts as a mismatch, including N/N, although refinement keeps
the existing ambiguity reward of -1. Literal MD reports original N bases; XM makes
no methylation call at an ambiguous reference base. These rules affect `aln --meth`
only, not stock aln or MEM. Different window spans
compete before collapsing each final locus (original contig, hypothesis, strand
and normalized start) to its lowest-cost refined CIGAR. Equal costs use canonical
CIGAR/end ordering within a locus, except for the neutral conversion tie-break below. This ranks reward-refined CIGARs by stock cost;
it does not minimize that cost over every possible CIGAR. Stock aln instead
cost-selects search paths before separately reward-refining its output.

`X0` counts best-cost distinct loci. `X1` counts all retained worse-cost loci
through best cost plus `-M`, or all returned worse costs with `-N`, not only the
next distinct cost. Primary selection and MAPQ use stock aln count/cost accounting
and its effective mismatch ceiling. Selection among best loci is seeded per read
and independent of thread count. MAPQ is the stock heuristic, not a calibrated
probability or a certificate that no better undiscovered locus exists. For TAPS/5base neutral scoring (including CpG-only context), equal minimum-cost
loci are ranked by conversion count from the winning raw CIGAR. Same-locus windows
first compete by `(cost, conversions)`, then canonical CIGAR/end ties, keeping the
selected CIGAR and counts together. Conversions still cost 0. Only the fewest-
conversion best-cost loci remain in `X0`; displaced equal-cost loci join `X1`.
The primary is drawn independently per read from that reduced best set. Two
remaining best loci still have MAPQ 0. Otherwise the stock MAPQ is capped at
`min(60, 6 * delta)`, where `delta` is the nearest displaced equal-cost locus's
extra conversions. This is the existing heuristic ceiling, not a new calibration
claim. EM-seq and other scoring modes retain cost-only placement ranking.

Primary and `XA` output use stored CIGARs without a second refinement. Ranking
counts and reward include raw terminal gaps. Display changes terminal insertions
to soft clips and removes terminal deletions with coordinate adjustment. Trimmed
tails are restored for primary and alternative alignments on either strand.
`NM`, literal `MD` and `XM` describe the normalized displayed alignment, so `NM`
can differ from the raw ranking mismatch/gap count. Records use Bismark
`XR:Z`, `XG:Z`, `XM:Z` tags selected by `--meth-tags`, rather than stock
`XM:i`, `XO:i`, `XG:i`. The header includes `@PG ID:bwa-mem3-meth`.

### Pairing and limits

Paired input uses stock sampe cost-based pairing on retained original-reference
loci, restricted to same-contig pairs, and the shared insert-size inference policy
with original reference length, counting only same-contig pairs. Mate rescue is
bounded to the anchor's contig. The insert-size and rescue reports follow `-v`.
Pairing considers every retained distinct locus when neither end has a lookup/work
cutoff and each end has at most `--max-occ` loci (default 100000, stock sampe `-o`).
The fused command's `-o` controls gap opens, not this pairing limit. Equality at
`--max-occ` is eligible. A lookup or `-m` cutoff on either end prevents ordinary
pairing and rescue confidence promotion, without discarding its best-known repeat
placement. There is no separate 500-locus or 255-survivor stop in this path.
If either end exceeds `--max-occ`, candidate pairing is skipped entirely, not run
on a subset. Ordinary mate rescue may still run when enabled (without `-s` or
`-A`), insert-size inference has succeeded, neither end is incomplete, the pair
is not already proper, and at least one end has MAPQ >= 17.
As in stock sampe, the pairing sweep pairs each reverse-strand locus only with the
two nearest preceding forward-strand loci of the mate on the same contig, so an
insert-size-compatible locus farther back is passed over when nearer ones intervene.

Here, complete means that configured discovery finished without those explicit
cutoffs. Pairing still sees only the configured retained cost strata: best through
best plus `-M`, or all returned worse costs under `-N`, including conversion-displaced
equal-cost loci. It does not see every genomic hit or every possible
CIGAR, and neither the absence of cutoff flags nor `-N` proves exhaustive discovery.
Unchanged conversion-selected ends retain that confidence ceiling through pairing
and rescue. An end moved by stock pairing uses stock pairing-derived confidence
(`seQ` is reset to zero), rather than treating its old single-end preference as
support for the new locus. Incomplete ends cannot gain confidence through pairing
or ordinary rescue.

A lookup-cap-bound end can instead be searched near each resolved best-cost mate
locus, provided the mate has no cutoff flags and at most 255 such loci. If there
are more, this rescue is skipped rather than truncating the competitors. Windows
are clipped at both ends to the anchor's contig and identical windows are evaluated
once. Accepted rescue CIGARs must align the trimmed query end to end and cost no
more than the end's best-known cost, when one exists. Raw chemistry-aware counts
and `-M/-O/-E` costs come from the owned realigned CIGAR, with neutral TAPS/5base
conversion ordering among equal pair scores. Pair scores include stock insert-size
penalties, and equal choices use deterministic original-contig/coordinate ordering.

This is placement recovery, not confidence recovery. The rescued end remains
incomplete at MAPQ 0 with `XT:A:M`, no exact `X0`/`X1`/`XA`, and no proper-pair flag
on either end. The former ceiling of 15 cannot raise incomplete MAPQ above 0.
An unchanged unique anchor keeps its own MAPQ and neutral conversion guard; moved
or tied anchors get MAPQ 0. Old rescue accuracy or MAPQ-threshold measurements do
not validate this revised behavior. Rescue uses the same active-allocation and
cumulative CIGAR budgets as candidate refinement and restores trimmed tails before
legacy CIGAR splitting.

The per-read SA lookup budget is a single total of 4596 by default, not separate
best/second-stratum budgets. Lookup or `-m` exhaustion retains a best-known
placement when available, marks it as a repeat with MAPQ 0, and omits exact
`X0`/`X1`/`XA` claims. Candidate rows are resolved in batches; once a batch holds
more rows than the budget has left, the read is incomplete whatever they resolve to,
so that batch resolves at most 64 more rows (in discovery order) to choose the
best-known placement instead of spending the rest of the budget. `-N` does not remove these resource limits. Candidate CIGAR
retention is capped cumulatively at 64 MiB per read. Refinement separately checks
a conservative 64 MiB active-allocation estimate (not total process RSS or retained
backend scratch), and query/window lengths must not exceed 65535 bases. Invalid
results and resource failures exit with an error. A refined candidate with more
than 255 mismatches, gap opens or gap extensions cannot be recorded and is dropped;
if its cost is within `-M` of the best retained cost (or with `-N`), the read is
marked incomplete as above (MAPQ 0, no exact `X0`/`X1`/`XA`). `-v 4` reports
lookup, window, dropped and incomplete counts.

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
- **Header without `--meth`.** `@HD` (unless the `-r` text has an `@HD` line), `@SQ` (unless
  the `-r` text has an `@SQ` line), the `-r` line, then `@PG`.

---

**See also:**
[mem](mem.md) · [index](index-cmd.md) · [shm](shm.md)
