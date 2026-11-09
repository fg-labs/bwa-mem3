# Counting alleles at known SNVs

`bwa-mem3 aln --known-snvs FILE` places short reads as if the genome held both alleles of every
SNV in FILE, and tags each read that covers one of those sites with `XV:Z`. This page explains
what the tag records and which reads to count when you want allele counts without reference
bias. The option itself is described in [aln](../cli/aln.md#known-snvs---known-snvs); it was added
in [#770](https://github.com/fg-labs/bwa-mem3/pull/770).

## Why short reads need it

A read that carries the ALT allele of an SNV has one more mismatch against the reference than the
same read carrying REF. For a 100-base read that hardly matters. For a read of 20-30 bases,
`aln` allows only about two differences, so the extra mismatch can:

- push the read over its mismatch budget, leaving it unmapped;
- make it tie with a paralog that matches it equally well, giving MAPQ 0;
- drop its MAPQ below the usual cut-off of 20.

REF reads at the same site are not affected, so a plain count of REF and ALT reads at MAPQ 20
under-counts ALT. `--known-snvs` removes the extra mismatch from an ALT read's cost, so it is
placed as it would be on a genome carrying that ALT.

Placing ALT reads correctly is half the job. The other half is choosing which reads to count, and
that is what `XV:Z` is for.

## Reading `XV:Z`

`XV:Z` has one entry for each REF–ALT pair of each site in FILE where the read's alignment has a
read base, separated by `;`. A multi-allelic site gives one entry per ALT, all at the same contig
and position. A site that the alignment spans with a deletion has no read base and gets no entry:

```text
contig:pos:REF>ALT,allele,away,max_diff,twin_same,twin_mapq
```

| Field | Meaning |
|---|---|
| `contig:pos:REF>ALT` | The site, with a 1-based position. In the contig name, `%`, `:`, `,` and `;` are percent-escaped as `%25`, `%3A`, `%2C` and `%3B`; split entries and fields before decoding the contig once. |
| `allele` | The read's base at the site: `R` (REF), `A` (ALT) or `O` (another base). Under `--meth`, a converted base counts as the allele it came from, and `H` means the conversion could turn either allele into the base the read shows, so the read cannot tell them apart. Under `--meth-context cpg` conversion-based equivalence applies only at a CpG, as in the alignment. Literal REF and ALT matches still produce `R` and `A` outside CpGs; a base matching neither is `O`. |
| `away` | The read's differences outside this site, counted as `aln`'s gap search counts them against `-n`: mismatches, gap opens and, unless `-e` is positive, gap extensions. |
| `max_diff` | The read's mismatch budget, set by its quality-trimmed search length and `-n`. |
| `twin_same` | For single-end input, 1 only if the twin has one best hit on the same strand within 5 bases, with its whole reference span on the read's contig; an equal-best hit elsewhere gives 0. For paired input, 1 if the twin end is mapped on the same strand within 5 bases after pairing and rescue, with its whole span on the read's contig and, without `--meth`, either an unmapped mate or a mapped mate whose span does not cross a contig end; else 0. Under `--meth`, hits and rescue never cross contigs. |
| `twin_mapq` | The twin's MAPQ. |

The **twin** is a copy of the read with REF and ALT swapped at the site. An ALT read's twin
carries REF and a REF read's twin carries ALT. `aln` aligns the twin the same way as the read and
reports where it lands. For paired input the twin is aligned as a pair, with a copy of the mate.
The mate's base is swapped too only when the mate has an `R` or `A` vote for the same REF–ALT
entry; otherwise, including an `O` or `H` vote there, the copy is the mate unchanged. The twin
only adds these two fields:
it never changes the read's own record.

The twin fields are `NA,NA` for:

- every entry after the first `R` or `A` entry, because each read is twinned at one site only;
- `O` and `H` entries;
- the selected entry when its plain twin search is incomplete.

Without `--meth`, a read with more than 256 candidate loci, longer than about 160 bases, or whose
genomic or ALT-window search stopped at the `-m` limit keeps its stock search result, but its first `R` or `A` site still gets a twin. For paired input, pairing and ALT-aware mate rescue still use that end with its mate, so its final record can still differ from stock's, for example when its mate is credited. The twin fields are numeric
when that search completes. Twins use the original quality-trimmed bases, including the mate's
trimmed bases for paired input.

CIGAR is refined against the selected ALT-substituted reference. MD and NM are computed against
the original reference, so they include the ALT difference, with the existing conversion
accounting under `--meth`. In ordinary output `XM:i` gives the read's mismatch count with the ALT
treated as a match.
Under `--meth`, `XM:Z` is the methylation call string and there is no such count; use `away`.

### Examples

Three 24-base reads over the same site, from an ordinary run (no `--meth`) in which each read's
own MAPQ passes the cut-off. Each has a mismatch budget of 2:

```text
XV:Z:chr7:1048576:C>T,A,0,2,1,37
XV:Z:chr7:1048576:C>T,R,1,2,1,37
XV:Z:chr7:1048576:C>T,A,2,2,1,25
```

- The first read carries ALT, has no other differences, and its REF twin lands in the same place at
  MAPQ 37. Count it as ALT.
- The second carries REF with one difference elsewhere, and its ALT twin places too. Count it as
  REF.
- The third carries ALT, but it already has two differences outside the site, its whole budget.
  Do not count it (see the budget rule below).

A read over two sites:

```text
XV:Z:chr1:200100:G>A,R,0,2,1,30;chr1:200112:A>G,A,0,2,NA,NA
```

It can be counted at `chr1:200100` only. It has no twin at `chr1:200112`, so it is not counted
there.

## Which reads to count

For unbiased allele fractions under `--meth`, first exclude both alleles at a site on a strand
and context where either allele can undergo conversion: C in the CT hypothesis, or G in the GA
hypothesis. Apply this exclusion to REF and ALT alike, even when a particular observation is
individually identifiable as `R` or `A`. For example, exclude both C and A at C>A on active CT,
and both G and T at G>T on active GA; excluding only C/T or G/A allele pairs is insufficient.
At T>C on active CT, T can be literal REF or converted ALT (`H`), while literal C can be `A`:
that individually identifiable C must still be excluded from unbiased fractions.

Use the alignment's chemical hypothesis (`XG:Z:CT` or `XG:Z:GA`), not just read number or SAM
orientation. Keep `XG:Z` in the output: the default `--meth-tags` does, `--meth-tags none`
drops it, and `--meth-tags XG` keeps it alone if no other methylation tag is needed. Do not
infer the hypothesis of a record that lacks `XG:Z` from its read number or orientation. With
`--meth-context cpg`, conversion is active only when the following base is G
for CT, or the preceding base is C for GA, in the allele-substituted original-reference context.
The neighbor must be on the same contig; an ambiguity at either CpG base disables conversion.
Evaluate both allele contexts and exclude both alleles if either permits source-base conversion.
Outside an active conversion context, and at sites without that hypothesis's source base, keep
the rules below. For overlapping mates, exclude the fragment if either covering mate triggers
this exclusion; do not retain it by counting the other mate alone.

Numeric twin fields report placement, not counting eligibility. This is a downstream counting
restriction: `aln` does not suppress these records or change their `R`/`A`/`H` or XV fields.

After that exclusion, count a read at a site when all of the following hold:

1. **Allele `R` or `A`.** `O` is neither allele, and `H` cannot tell them apart.
2. **The read has MAPQ at or above your cut-off** (20 is usual).
3. **Its twin also passes:** `twin_same` is 1 and `twin_mapq` is at or above the same cut-off.
4. **`away` is at most `max_diff - 1`.**

Each rule removes one source of bias.

- **The twin rule** makes placement independent of the allele. A read is counted only if it would
  have placed at the same locus, just as confidently, carrying the other allele. This keeps out
  sites where one allele is unique but the other matches a paralog, and it also removes the
  reference bias left in stock placements: REF reads whose ALT twin would not have placed are
  dropped too.
- **The budget rule** stops sequencing errors from becoming ALT counts. A read with `max_diff`
  differences elsewhere is already at its budget. Stock `aln` accepts its REF version but rejects
  its ALT version because that adds another mismatch. With `--known-snvs` both versions can place;
  an ALT base at that depth is often a sequencing error. Its REF twin has the same differences
  elsewhere and places as well, so the twin rule cannot catch it. Requiring one spare difference
  means the base at the site never decides whether the read places.

The conversion-source exclusion applies to both alleles before these four rules; the rules
alone do not guarantee symmetric methylation counts. Multi-allelic sites have a separate
limitation below. The cost is depth. How much the budget rule
removes depends on `-n` and on how far the reads differ from the reference. With `-n 0` every read has
`max_diff` 0, so no read passes and nothing is counted. Sites in repeats, where neither allele
places uniquely, end up with few or no counted reads instead of a biased count.

At a multi-allelic site each ALT has its own entry, but a read is twinned at its first `R` or `A`
entry only, so the twin rule does not treat the ALTs alike. For a `C>A,G` record, a REF read gets
twin fields on the `C>A` entry and `NA,NA` on the `C>G` entry, while a G read is `O` on `C>A` and
gets twin fields on `C>G`. Counting `C>G` with these rules therefore drops every REF read but can
keep G reads. When you need counts for each ALT, split the record into one biallelic VCF per ALT
and run once with each. That works only when the target entry gets the read's twin: if the VCF
keeps another site that the read covers at a lower position, an `R` or `A` vote there takes the
read's only twin and the target entry stays `NA,NA`. When other sites can take that slot, run
with a VCF holding only the target allele.

## Paired-end reads

Each mate has its own `XV:Z`. Count fragments, not mates:

- If only one mate covers the site, the fragment's count is that mate's.
- If both cover it, count the fragment once, and only when both mates pass the rules above. If
  the mates disagree at the site, the fragment is `O`.

The twin of a paired read is a pair, so `twin_mapq` is the twin pair's paired MAPQ.

## What it does not handle

- **Nearby sites.** Each window carries one site's ALT at a time, so a read that carries the ALT at
  two sites is credited for one of them only, and may still be misplaced. This holds under
  `--meth` too. A mate-rescue window, or under `--meth` a candidate window, with more than 32
  site alleles tries only its first 32 (in position order); the rest stay at REF.
- **Long reads.** Without `--meth`, reads longer than about 160 bases keep their stock search
  result; for paired input, pairing and ALT-aware mate rescue can still change the final record,
  for example when the mate is credited.
  Full-read ALT-window discovery, including under `--meth`, requires a read to fit within a window
  of at most 321 bases. It cannot independently discover an over-budget read that does not fit.
  Methylation genomic refinement can still credit longer reads at discovered loci; paired-end
  rescue separately tries ALT-substituted rescue windows and is not limited by that discovery width.
- **Work limits under `--meth`.** If a read's ALT-window search hits a work limit (the lookup
  budget or the `-m` stack cap), the read is reported as incomplete even when none of its
  placements covers a site: it gets MAPQ 0 with no `X0`/`X1`. For paired input, neither end of
  that pair can then be promoted by pairing or mate rescue, so the mate's record can change too.
- **Long unseeded reads under `--meth`.** With known sites loaded, an unseeded read longer than
  160 bases (possible only when `-l` exceeds 160) cannot stop early, so a highly repetitive one
  can use up the lookup budget.
- **Early stopping near sites under `--meth`.** An ALT-window placement that covers a site with at
  most one mismatch more than the best placement found blocks early stopping, even when the same
  locus is also found at lower cost without the ALT. In a repetitive family this includes reads
  that carry the REF base, so a read that matches one locus exactly can end at MAPQ 0.
- **Paralogous sites.** If a site's ALT allele matches a paralog somewhere else in the genome, reads
  from the paralog carry the ALT. Crediting the ALT does not break an exact tie: a paralog read that
  matches the site and the paralog equally well gets MAPQ 0 and fails the MAPQ rule. A paralog read
  that still places at the site and passes the placement and twin rules (for example, one whose
  other differences fit the site better than the paralog) is counted as ALT, and no read-level rule
  removes it, so some paralog bias can remain at such sites.
- **Indels and MNVs.** Only single-base substitutions in FILE are used.

Without `--meth`, at verbosity 3 or higher, `aln` prints how many reads were credited, how many
were at a site with no change, and how many kept their stock result:

```text
[M::aln] --known-snvs: 1203 reads credited, 51877 at a site with no change, 12 kept stock (too many hits, too long, or a search stopped at -m)
```

---

**See also:**
[aln](../cli/aln.md#known-snvs---known-snvs) ·
[Output: SAM/BAM, headers, tags](output.md) ·
[Methylation SAM tags](../methylation/tags.md)
