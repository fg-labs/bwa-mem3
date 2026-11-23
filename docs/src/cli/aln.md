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

**Byte-identity without `--meth`, `--short-insert`, `--ref-span-insert`, `--known-snvs` or `--exact-rescue-window`.** The tested SAM records and headers matched the output of stock
bwa 0.7.19 `aln` + `samse`/`sampe`, except for the `@PG` line, which names the
program that ran, and an `-r` text that carries an `@HD` line (see **Header without `--meth`**
under [Notes](#notes)). The comparison uses generated single-end and paired-end reads
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
above does not cover them ([#749](https://github.com/fg-labs/bwa-mem3/pull/749)).
See [Methylation](#methylation---meth) below.

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
| `--exact-rescue-window` | new | off; on under `--meth` | Round the mate-rescue window outward from the insert range instead of truncating it as sampe does (see [Mate-rescue window bounds](#mate-rescue-window-bounds)); paired input required. `--no-exact-rescue-window` turns it off, including under `--meth` |
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
given. `--short-insert`, `--ref-span-insert` and `--exact-rescue-window` instead require two
input files and reject single-end use before loading the index or creating output
(`--no-exact-rescue-window` is accepted with one).

### Short, overlapping fragments

Use `--short-insert` when short fragments overlap and a longer read in the batch would
otherwise exclude them from insert-size estimation, or when both reads often cover the whole
fragment ([#752](https://github.com/fg-labs/bwa-mem3/pull/752)). It works with ordinary and
methylation-aware (`--meth`) paired alignment:

```bash
bwa-mem3 aln --short-insert ref.fa r1.fq r2.fq > out.sam
```

The option sets the estimator's lower bound to the lower interquartile fence, with a
minimum of one base, rather than the longest read in the batch. It keeps the upper
fence, minimum number of high-quality pairs and previous-batch fallback: a failed
estimate reuses the previous batch's estimate when there is one, and otherwise pairing
falls back to `-a` without an estimate. `-A` overrides the option and disables insert
estimation and mate rescue.

The option also accepts a pair whose reverse read starts at the same base as its forward
read. Stock pairing requires the reverse read to start strictly after the forward one, so
when trimmed reads each cover a whole fragment shorter than the read length, almost no pair
is proper and neither mate's MAPQ benefits from the other. The outer span must still be at
least the longer read, but that read is measured after `-q` trimming: trimmed bases need not
lie inside a fragment shorter than the reads (stock, and `-A`, count the whole read). The
other geometry requirements are kept: the reverse read may not start before the forward read,
and the insert must fit the estimate (or `-a` without one).

Mate rescue (and, with `--meth`, cap rescue) may also search the anchor's own span. Stock
rescue looks for the mate only beyond the anchor, so it can never place a mate that covers the
same fragment. A rescued mate must stay inside the fragment the anchor implies: the forward read
may neither start nor end after the reverse read. Rescue is a local alignment, so a mate that
runs past its anchor (for example, read 1 starting with adapter sequence ahead of the fragment)
is rescued with that overhang soft-clipped, also when the overhang runs past the end of the
reference: the window needs room only for the aligned part of the mate (at least 20 bases),
where stock rescue needs room for the whole read.

This is an opt-in departure from stock output. The inferred distribution, the coincident
starts, the trimmed-length outer span and the wider rescue windows can change pair selection, proper-pair flags, MAPQ, rescued
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
alignments' reference spans (read length plus deletions minus insertions) instead
([#765](https://github.com/fg-labs/bwa-mem3/pull/765)):

```bash
bwa-mem3 aln --ref-span-insert ref.fa r1.fq r2.fq > out.sam
```

With the option:

- insert-size estimation takes each observation to the later read's reference end (for mates
  that start at the same base, the reverse read's, as pairing measures it) and floors its lower
  fence at the longest reference span in the batch rather than the longest read;
- pairing compares the insert to the estimate (or to `-a`) and scores it from the reverse
  hit's reference end, requires the outer span to cover both hits' reference spans (without
  `--short-insert`, or with `-A`, each also extended by its `-q`-trimmed bases, as stock counts
  the whole read), places a
  reverse hit from its own span rather than the read's best hit's, and gives a read the hit
  pairing selects, with that hit's span, even when it starts at the same base on the same
  strand (stock keeps the read's first hit, whose span the gapped refinement then uses);
- mate rescue, and `--meth` cap rescue, place a forward mate's window and (without
  `--short-insert`) floor a reverse mate's at the anchor's reference end rather than its
  read-length end.

The TLEN calculation is unchanged: TLEN is taken from each record's CIGAR with or without the
option. Its value can still change, because pairing can select a different hit, with a
different span, for a read. `--meth` cap rescue already measures inserts from spans;
the option also makes its minimum outer span the two alignments' spans, as in pairing. With
`--short-insert`, a coincident-start pair must cover both hits' spans, and rescue's ceiling for
a reverse anchor is its reference end either way. `-A` does not turn the option off: pairing
still measures the insert it compares to `-a` from the spans.

This is an opt-in departure from stock output: it can change the inferred distribution, pair
selection, proper-pair flags, MAPQ, rescue windows and, for reads pairing moves to a gapped
hit, their position and CIGAR. Without the option the output is unchanged.

### Mate-rescue window bounds

Mate rescue searches for an anchor's mate in a window placed from the batch's inferred insert
size: for a forward anchor at `pos`, its reverse mate of `len` bases is searched from
`pos + avg - 3 std - 1.5 len` to `pos + avg + 3 std + 0.5 len` (mirrored for a reverse anchor).
Stock `sampe` truncates the start to an integer and then truncates start + `6 std + 2 len`, so the
end carries the start's discarded fraction and falls up to two bases short of the stated bound.
It also moves with `floor(6 std)` alone: every pair with the same anchor and mate geometry has
the same window offsets, so when the batch's inferred std crosses a multiple of 1/6 (for
example from 12.6685 to 12.6646, across 76/6), the end of all of their windows moves by a base
at once, and it can move back as std grows. For a cohort of pairs whose mate sits at the window's
edge, such as a spike-in with inserts well above the batch's, a 0.004 change in std can change
how every one of them is rescued.

`--exact-rescue-window` computes each bound as its offset from the anchor and rounds it outward
(the start down, the end up), so the window is the smallest one that holds the stated range:

```bash
bwa-mem3 aln --exact-rescue-window ref.fa r1.fq r2.fq > out.sam
```

The start is the same as stock's; the end is up to two bases later. Each bound is monotone in
the inferred mean and std, and it moves only when the stated bound itself crosses a base, which
for the end happens once per 1/3 of a unit of std rather than up to three times as often. A
cohort whose mate sits exactly at a bound still moves with it there: an integer window must
step somewhere. The same window is used by `--meth` cap rescue. The option combines with
`--short-insert` and `--ref-span-insert`, which only move the clamps and the anchor's end.

The option is on by default under `--meth` (whose output has no stock counterpart), and
`--no-exact-rescue-window` turns it off there. Plain `aln` keeps stock's window unless the option
is given. Enabling it can change which mates are rescued and their CIGARs, MAPQ and proper-pair
flags ([#779](https://github.com/fg-labs/bwa-mem3/pull/779)).

### Known SNVs (`--known-snvs`)

Short reads that carry the ALT allele of an SNV have one more mismatch against the reference
than reads carrying REF. For reads of 20-30 bases that mismatch often costs the read its
placement: it ties with a paralog, drops below MAPQ 20, or goes over the mismatch budget. Allele
counts at those sites are then biased toward REF. `--known-snvs FILE` places reads as on a genome
that holds both alleles of every SNV in FILE ([#770](https://github.com/fg-labs/bwa-mem3/pull/770)):

```bash
bwa-mem3 aln --known-snvs sites.vcf.gz ref.fa reads.fq.gz > out.sam
```

FILE is a VCF or BCF (plain or compressed). Its single-base substitutions are used, multi-allelic
ones as alternatives (never combined); other records are skipped, and FILTER is ignored. A record
whose contig is not in the reference, or whose REF differs from it, is an error naming the record.

How it works:

- A window of 160 bases either side of each site, with the ALT in place, is indexed when the run
  starts, in memory, with up to `-t` threads (no files are written).
- Each read is also searched against the windows. A hit there that covers its window's site is
  a hit at the same reference position with the ALT credited, and the cheaper of the two costs
  wins. Then MAPQ, the primary hit and the alternative hits are recomputed with `aln`'s usual
  cuts. The stock primary is kept while it is still among the best.
- Paired end: the insert size is estimated from stock placements first. Mate rescue also aligns
  the mate against its window with each site allele substituted alone, every other site at REF,
  for at most the window's first 32 site alleles in position order (later ones stay at REF), so at
  most 33 alignments per window. It keeps the best score. An alignment that `--short-insert` rules
  out of the fragment never replaces another, and equal best scores at loci that do not overlap
  leave the mate unrescued, whether they come from two alignments or from one window whose best is
  tied: a pick survives only if it scores above every such tie.
- `--meth`: the converted search also resolves the windows that the ALT-window search finds. Each
  candidate window holding a site is scored along the ALT-window search's path to it (its gaps,
  including a gap at either end of the read, as the converted search scores its own rows), against
  the reference and then once per site allele (each ALT of a multi-allelic site separately) with
  only that ALT in place, under the same 32-allele cap. A candidate whose path cannot be recovered,
  and each mate-rescue window holding a site, is realigned instead, the same way (a mate-rescue
  window counts only the sites on its anchor's contig, since rescue never aligns past that
  contig). The locus keeps the cheapest alignment; on a tie the reference wins, then the lower
  (position, ALT). An alignment in a window that only the ALT-window search proposed is kept only
  within the read's difference limit, as the read's own search would be. When a realigned
  candidate's best-scoring alignment breaks the read's `-n`, `-o` or `-e` limits, the window is
  realigned once more with its indels held to the window's own extent (ungapped when the window is
  the read's length), and that alignment is kept if it is within the limits. A pair whose
  placement assumes an ALT is left out of the insert-size estimate.
- Without `--meth`, reads longer than about 160 bases, reads with more than 256 candidate
  loci, and reads whose genomic or ALT-window search stopped at the `-m` limit (it may have
  missed equal-best loci) keep their stock search result. For paired input that end still goes
  through pairing and ALT-aware mate rescue with its mate, so its final record can still differ
  from stock's, for example when its mate is credited. Under `--meth` those override limits do not apply to genomic
  refinement: discovered loci can receive ALT credit under the converted search's own limits
  (its SA lookup budget and `-m`). Full-read ALT-window discovery still requires the read to fit
  inside a fixed window of at most 321 bases; it cannot independently discover an over-budget
  long read that does not fit. Genomic refinement can credit longer reads at discovered loci,
  and paired-end mate rescue separately tries ALT-substituted rescue windows, not these fixed
  discovery windows.

CIGAR is refined against the selected ALT-substituted reference. MD and NM remain against the
original reference, so they include the ALT difference (with the usual conversion accounting
under `--meth`). Without `--meth`, `XM:i` is the mismatch count with the ALT credited; under `--meth`, `XM:Z` is the
methylation call string as usual, and no record carries an ALT-credited mismatch count (`XV:Z`'s
`away` field gives the differences outside each site). A read is unchanged when none of its
candidate placements covers a site, with three exceptions under `--meth`: batch-wide
insert-size effects for paired input; a read whose ALT-window search hit a work limit is reported
as incomplete (MAPQ 0, and for pairs neither end can be promoted by pairing or rescue); and an
unseeded read over 160 bases cannot stop early (see
[what it does not handle](../user-guide/known-snvs.md#what-it-does-not-handle)). For paired input
that must hold for both mates: when one mate is placed, moved or rescued
at a site, its mate's flags, mate fields and MAPQ change with it. Under paired `--meth`, excluding
ALT-credited pairs from the insert-size estimate can also change pairing or mate rescue for other
pairs in the batch, even when neither mate's candidates cover a site.

Every record with a read base aligned to a site gets `XV:Z`, one entry per REF–ALT pair at each
such site, separated by `;`: a multi-allelic site gives one entry for each ALT, all at the same
contig and position. Positions are 1-based. A site the alignment covers only with a deletion has
no read base, so it gets no entry. A record whose reference span crosses a contig boundary
(which `aln` flags unmapped; under `--meth` no placement crosses one) gets no `XV:Z` at all,
even when a read base aligns to a site:

```text
contig:pos:REF>ALT,allele,away,max_diff,twin_same,twin_mapq
```

- `allele`: `R` (REF), `A` (ALT), `O` (another base), or under `--meth` `H` (the conversion can
  explain both alleles on that strand). Under `--meth` a converted base counts as its source
  allele; under `--meth-context cpg` this equivalence applies only at a CpG, as the alignment
  scores it. Literal REF and ALT matches still produce `R` and `A` outside CpGs; other bases are `O`.
- `away`: the read's differences outside the site, counted as the gap search counts them
  (mismatches, gap opens and, unless `-e` is positive, gap extensions).
- `contig`: the reference contig name, with `%`, `:`, `,` and `;` percent-escaped as `%25`,
  `%3A`, `%2C` and `%3B`. Split the entries and fields before decoding the contig once.
- `max_diff`: the read's mismatch budget at its quality-trimmed search length.
- `twin_same`, `twin_mapq`: for the first `R` or `A` site only (`NA` otherwise), the read is
  realigned with REF and ALT swapped at that site (as a pair, with its mate, for paired input).
  For single-end input, `twin_same` is 1 only when the twin has one best hit, on the same strand
  within 5 bases, with its whole reference span on the read's contig; an equal-best hit elsewhere
  makes it 0. For paired input it is 1 when the twin end is mapped on the same strand within 5
  bases after pairing and rescue, with its whole span on the read's contig and, without `--meth`,
  either an unmapped mate or a mapped mate whose span does not cross a contig end. Under `--meth`,
  hits and rescue never cross contigs. `twin_mapq` is the twin's MAPQ. An incomplete plain twin
  search gives `NA,NA` even for the selected vote. The twin uses the original quality-trimmed
  bases and never changes the read's own record or the run's random draws.

For unbiased allele fractions under `--meth`, first exclude both alleles where either is a
conversion-source base in an active strand/context: C for the alignment's CT hypothesis, G for
GA. This includes C>A on CT and G>T on GA, not just confounded C/T or G/A pairs. At T>C on CT,
a literal C can be identifiable as ALT while T is hidden; exclude both from fractions anyway.
Use `XG:Z` for the alignment hypothesis, not read number or SAM orientation, so keep `XG:Z` in
the output: the default `--meth-tags` does, `--meth-tags none` drops it, and `--meth-tags XG`
keeps it alone. Do not infer a record's hypothesis when it lacks `XG:Z`. Under
`--meth-context cpg`, test the allele-substituted original-reference context: the following G
for CT or preceding C for GA must be on the same contig and neither CpG base may be ambiguous.
If either allele context permits source-base conversion, exclude both alleles. Outside that
context, retain the usual rules. For overlapping mates, either covering mate's exclusion removes
the fragment. Numeric twin fields describe placement only, not counting eligibility; this
restriction is for downstream counting and changes no alignment or tag.

After that exclusion, count a read with an `R` or `A` vote only when it passes the MAPQ
cut-off, its twin also passes (`twin_same` 1 and `twin_mapq` at or above the cut-off), and
`away` is at most `max_diff - 1`. The last rule matters: a read already at its mismatch budget
places only when an error at the site matches the ALT, so without it errors become ALT counts.
Only the first `R` or `A` entry is twinned, so at a multi-allelic site the rule is not symmetric
between the ALTs. For per-ALT counts, give each ALT its own biallelic VCF and run once with each;
that works only when the target entry is the read's first `R` or `A` entry. If another site the
read covers can take the twin, run with a VCF holding only the target allele.
See [Counting alleles at known SNVs](../user-guide/known-snvs.md) for worked examples, paired-end
fragments and the limits.

This is an opt-in departure from stock output. It can change the position, MAPQ, flags and
alternative hits of reads at the sites. Without the option the output is unchanged. With a VCF
that has no SNV, the alignment records and every header line but `@PG` are also unchanged; `@PG`
records the command line, which now includes the option.

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
  that names one of the input files (the reads or the `--known-snvs` VCF, also
  when it is the file on standard input for `-`), or a file of the index (`.pac`, `.ann`,
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
df -h /dev/shm                                        # check the current size first
sudo mount -o remount,huge=advise /dev/shm            # allow huge pages; keeps the current size
sudo mount -o remount,size=20g /dev/shm               # only if too small: room for the segment
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
  (`size=`; containers often have a small one). Change `size=` only if
  `df -h /dev/shm` shows it is too small, and never set it below the current size:
  the limit applies to every process that uses `/dev/shm`. A remount that changes
  only `huge=` keeps the current size, so it needs no size check. Without
  `huge=advise` (or `huge=within_size`) the segment stays on base pages (usually
  4 KB), which costs CPU on every lookup. With it, huge pages are allowed but
  best-effort: the kernel backs the segment with them only when it can (see
  [`shm`](shm.md)). Re-stage after a remount.
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

Converted-space backtracking finds hits within the user's search limits (`-n`, `-o`,
`-e`, seeding options, `-M/-O/-E`, `-N` and `-m`). Each SA row of a hit is a candidate
locus. Only the directional copy appropriate to the read is retained.

**Stock cuts.** A row costs at least its path's converted score (below), so the search takes
stock aln's cuts from the best verified row cost `B`. Rows are resolved at the first hit and
at every converted-score boundary. The search stops past `B + M` (stock's score cut). Once
the best stratum is complete (`B` is below the next entry's score), every row that can cost `B`
has been resolved, and the difference ceiling drops to the fewest differences among the stored
rows at `B`, plus one (stock's top-2 rule, at most `-n`); the retained loci obey the same
ceiling. More than `-R` best loci within that ceiling then end the search, as stock's top2b stop does: `X0` is then the full best count and `X1` is
not printed ([#785](https://github.com/fg-labs/bwa-mem3/pull/785)). These cuts apply unless `-N` is given, a penalty is negative or `--known-snvs` is
set. Under `--known-snvs` (whose ALT haplotypes can cost less than the converted path) the
search instead stops early only on the earlier certified bound: for current best cost `B`,
let `T = B + M` and `k = floor(T / M)`. When `M > 0`, `T < O`, `E >= 0`, the configured
difference ceiling admits `k` differences, `-N` is off, discovery has not hit a cap, the read is
at most 160 bases long and the read's ALT-window search gives a lower bound, greater than `k`, on
the mismatches of every ungapped placement over a known site with its ALT in place, the search can stop after fully evaluating converted scores
through `M*k`; a seeded read only when `k` is also at most `-k`
([#778](https://github.com/fg-labs/bwa-mem3/pull/778)). If any condition fails, discovery continues. This does not make discovery exhaustive or guarantee a globally optimal alignment.

**Saturation certificates (single-end).** Under the stock cuts, a single-end read's rows are
verified stratum by stratum in converted score, narrow hits first and wide ones (more than 256
rows) in a pseudo-random order seeded from the read
([#786](https://github.com/fg-labs/bwa-mem3/pull/786)). Once the best stratum is complete (the
best cost `B` is below the current stratum), every row that can cost `B` is verified, so the best
loci and the difference ceiling (their fewest differences plus one) are final. From then on the
search ends as soon as more rows cannot change the stock output; every count below is of the loci
within the ceiling:

- *repeat*: more than `-R` loci at `B` (under the TAPS/5base tie-break, more than `-R` with the
  fewest conversions, or at least 179 best loci, which gives stock MAPQ 0 whichever has the
  fewest, with a best cost below the read's `-n` mismatches). MAPQ is 0 and `X0` is the full
  count (`X1` is not printed, as for stock's `-R` stop). The primary is drawn uniformly over
  every row of the hits that can hold a best row (rejection sampling, seeded from the read),
  among the best-cost rows within the ceiling. After 256 draws without such a row, the primary
  is drawn among the best loci instead.
- *crowded*: a unique best locus (under the TAPS/5base tie-break, a unique one among the best loci
  with the fewest conversions; not using all of `-n`, where stock MAPQ is 25) and at least 255
  retained worse loci: MAPQ 0, `X0` 1 and `X1` the verified count, a lower bound.

Other reads verify every row up to a row budget of 16,384 (4,596 without the certificates: with
`-N`, a negative penalty, `--known-snvs` or before pairing is decided for paired-end input). A
tie-break repeat whose best stratum is still open (more than `-R` best loci so far, all with
conversions) ends incomplete past 4,596 rows, as before. A read that reaches the budget without a certificate is incomplete, as below.

**No-hit proof.** Unless `-e` is given a positive value (which stops counting gap extensions as
differences), a read allowed at least three differences, with no hit
in its exact stratum and a difference to spare over its absent-substring bound, is first checked
for any hit at all. Its difference limit is the read's own: an integer `-n`, or the per-length
limit a fractional `-n` gives. The read is cut into consecutive pieces that allow one edit each,
plus one exact piece when that limit is even, so that the allowances plus one sum to the limit
plus one: any alignment within the limit aligns some piece within its allowance. Every occurrence
of each piece within its allowance is found on the converted index, and the read's alignments
around each are bounded by a banded edit distance on the index text. When none can be within the
limit, the search ends with no hit, where the full search would also end with none; otherwise, or past
4,096 candidate rows, it searches as usual. The ALT-window search of `--known-snvs` is not
checked this way. The output can change: a no-hit read whose search used to stop at the `-m`
work bound, incomplete, is now proven complete; for paired-end input such an end can then be
placed by mate rescue ([#787](https://github.com/fg-labs/bwa-mem3/pull/787)).

**Path cost.** Each row is scored on the original reference along the search path that
found it, as stock aln scores a hit by its path: the path's alignment (its gaps where the
search placed them) is the row's alignment, and its cost is

```text
cost = M * mismatches + O * gap_runs + E * gap_extensions
```

with strict mismatches counted on the original reference. A gap run of length `k` has `k-1`
extensions. Chemistry-accepted conversions, including neutral conversions, cost 0. An
ambiguous original-reference or query base always counts as a mismatch, including N/N.
Original `.amb` intervals override random packed bases with N, without rebuilding either
index. A converted mismatch is always a strict mismatch and the gap terms are the same, so a
row's cost is never below its path's converted score. There is no window realignment in the
search ([#784](https://github.com/fg-labs/bwa-mem3/pull/784)). Two paths that reach the same interval with different gaps (a gap in a tandem repeat,
say) are both scored; the second reuses the first's SA rows. The search entry records at most
4 gaps and counts their extensions in 8 bits, so `--meth` refuses `-o` above 4 and `-e` above 255.
`--meth-context cpg` requires TAPS/5base and non-collapsed scoring; CpG labels use true
contig flanks and respect the ambiguity runs.

Each final locus (original contig, hypothesis, strand and normalized start) keeps its
lowest-cost path. Equal costs use canonical CIGAR/end ordering within a locus, except for
the neutral conversion tie-break below. Literal MD reports original N bases; XM makes no
methylation call at an ambiguous reference base. These rules affect `aln --meth` only,
not stock aln or MEM ([#751](https://github.com/fg-labs/bwa-mem3/pull/751)).

A row whose strict differences exceed the read's `-n` limit (an integer `-n`, or the
per-length limit a fractional `-n` gives), or whose strict differences inside the seed
exceed `-k`, is dropped, as stock aln never reports such a hit. Differences are counted as
stock aln counts them: mismatches plus gap runs, plus gap extensions unless `-e` is positive;
conversions are free. The seed's differences are those the search would count against `-k`.
The converted search cannot see a read C over a reference T (or, for read 2, a read G over a
reference A), so it can reach a row within `-n` (`-k`) converted differences whose
original-reference path has more. A dropped row is not a locus: it does not count toward
`X0`/`X1` or MAPQ and is never paired. Its SA lookup still counts toward the lookup budget,
since the row was resolved to find it ([#782](https://github.com/fg-labs/bwa-mem3/pull/782)).

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
claim. The ceiling belongs to the conversion-selected locus: a paired-end read
that pairing or mate rescue moves to another locus takes the pair's MAPQ instead. EM-seq and other scoring modes retain cost-only placement ranking ([#751](https://github.com/fg-labs/bwa-mem3/pull/751)).

Primary and `XA` output use the stored path CIGARs. Ranking counts include raw terminal
gaps (a path has none unless `-i` is 0). Display changes terminal insertions
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
penalties, and equal choices use deterministic original-contig/coordinate ordering
([#751](https://github.com/fg-labs/bwa-mem3/pull/751)). A candidate pair whose
insert-size penalty is not finite (an insert so far into the tail of a narrow estimate
that the Gaussian tail underflows to 0) is skipped, as cap rescue already skips one;
stock `aln` pairing without `--short-insert` instead scores it as `INT_MIN`, the value stock's x86-64 builds produce
(see **Insert-size penalty in pairing** under [Notes](#notes)).

This is placement recovery, not confidence recovery. The rescued end remains
incomplete at MAPQ 0 with `XT:A:M`, no exact `X0`/`X1`/`XA`, and no proper-pair flag
on either end. The former ceiling of 15 cannot raise incomplete MAPQ above 0.
An unchanged unique anchor keeps its own MAPQ and neutral conversion guard; moved
or tied anchors get MAPQ 0. Old rescue accuracy or MAPQ-threshold measurements do
not validate this revised behavior. Rescue checks the same active-allocation
estimate as candidate refinement and keeps its own 64 MiB cumulative CIGAR budget for each
rescued end; unlike the candidate cap, exceeding it exits with an error. It
restores trimmed tails before legacy CIGAR splitting.

The per-read SA lookup budget is a single total of 4596 by default (16,384 for a read
verified with certificates, above), not separate best/second-stratum budgets. Lookup or `-m` exhaustion retains a best-known
placement when available, marks it as a repeat with MAPQ 0, and omits exact
`X0`/`X1`/`XA` claims. Without certificates, candidate rows are resolved in batches; once a
batch holds more rows than the budget has left, the read is incomplete whatever they resolve to,
so that batch resolves at most 64 more rows (in discovery order) to choose the
best-known placement instead of spending the rest of the budget
([#777](https://github.com/fg-labs/bwa-mem3/pull/777)). A read verified with certificates
resolves its rows 64 at a time, checking the certificates after each chunk, until a certificate
ends it or the 16,384-row budget is spent; a read that spends it is incomplete, placed among the
rows it resolved. `-N` does not remove these resource
limits. Candidate CIGAR retention is capped cumulatively at 64 MiB per read: a read that reaches the cap is
marked incomplete in the same way, stays unmapped if it has no retained placement,
and later reads continue; the run does not exit. Window refinement (mate rescue, and
`--known-snvs` candidates without a search path) separately checks a conservative
64 MiB active-allocation estimate (not total process RSS or retained backend scratch).
A candidate window longer than 65535 bases (only a 65535-base read with a deletion
produces one) is dropped and the read is marked incomplete. Invalid results and other
resource failures exit with an error. A candidate row with more
than 255 mismatches, gap opens or gap extensions cannot be recorded and is dropped;
if it is within the `-n` limit and its cost is within `-M` of the best retained cost
(or with `-N`), the read is marked incomplete as above (MAPQ 0, no exact
`X0`/`X1`/`XA`). A candidate over the `-n` limit is dropped first and never marks
the read incomplete, even with `-N`. `-v 4` reports
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
- **Insert-size penalty in pairing.** Paired-end pairing converts stock's insert-size
  penalty to an integer; when it is not finite, which stock leaves undefined, every build
  uses `INT_MIN`, the value stock's x86-64 builds produce (`--meth` and `--short-insert`
  pairing skip such a candidate instead; see **Pairing and limits**). This preserves x86-64 behavior.
  On aarch64, the changed runner-up cost can change MAPQ at extreme mismatch penalties, even
  with an inferred insert-size estimate and unchanged winning positions (for example,
  `-M 214748364` can change MAPQ from 232 to 3)
  ([correctness note](../whats-different/correctness.md#aln-a-non-finite-insert-size-penalty-in-paired-end-pairing-pr-803);
  [#803](https://github.com/fg-labs/bwa-mem3/pull/803)).
- **Header without `--meth`.** `@HD`, `@SQ` (unless the `-r` text has an `@SQ` line), the `-r`
  line, then `@PG`. The `@HD` is the default `@HD VN:1.5 SO:unsorted GO:query`
  unless the `-r` text has a line that starts `@HD<tab>`: then that line, the
  first if there are several, is written first in its place, and the `-r` text
  after `@SQ` omits its `@HD` lines. This is the one header
  difference from stock: bwa 0.7.19 writes a `-r` `@HD` where it sits in the
  text, after `@SQ`, but the SAM specification requires `@HD` to be the first
  line (bwa 0.7.19 moved only its default `@HD` first;
  [lh3/bwa#442](https://github.com/lh3/bwa/issues/442);
  [#797](https://github.com/fg-labs/bwa-mem3/pull/797)). The input is still
  warned about, as in stock: `[W::aln] please don't include @HD with option -r.
  Continue anyway.` on standard error (stock's text names `-H`, which `aln` does
  not have). A header with no `-r` `@HD` is the same as stock's.

---

**See also:**
[mem](mem.md) · [index](index-cmd.md) · [shm](shm.md)
