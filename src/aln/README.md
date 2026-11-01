# src/aln: bwa `aln`/`samse`/`sampe` port

This directory holds the code behind `bwa-mem3 aln`, a fused port of bwa's
short-read `aln` + `samse`/`sampe` pipeline whose output is byte-identical to
stock bwa 0.7.19.

## Provenance

| Item | Value |
| --- | --- |
| Upstream repository | https://github.com/lh3/bwa |
| Branch | `Apache2` |
| Commit | `a08c7c0` (a08c7c05517ab534b1618514507e95cfeaf8cdc7) |
| License | Apache License 2.0 (`LICENSE` in this directory); relicensing statement at https://github.com/lh3/bwa/issues/389 |

Code is copied only from the `Apache2` branch, never from `master` or a release
tag. Upstream's `bwt.c` (in `aln_bwt.*`) and `bntseq.c` (in `bwase.cpp`) carry an MIT
notice, and those files keep it verbatim. Upstream's `bwtaln.c`, `bwtgap.c`, `bwase.c`,
`bwape.c` and `bwaseqio.c` carry no copyright block; the "Heng Li and the bwa
contributors" line in the headers of the files derived from them was added here.

## Files

| File | Origin |
| --- | --- |
| `bwtgap.h`, `bwtgap.cpp` | `bwtgap.h`/`bwtgap.c`: the backtracking search (`bwt_match_gap` and its priority stack) |
| `bwtaln.cpp` | the search half of `bwtaln.c`: `gap_init_opt`, `bwa_cal_maxdiff`, `bwt_cal_width`, `bwa_cal_sa_reg_gap` and the thread fan-out, plus `width_walk_init`/`width_walk_run` (new) |
| `bwase.h`, `bwase.cpp` | `bwase.h`/`bwase.c`: the SAM half of samse (`bwa_aln2seq_core`, `bwa_sa2pos`, `bwa_cal_pac_pos_core`, `bwa_refine_gapped`, `bwa_print_sam1`), plus `bns_cnt_ambi` from `bntseq.c` |
| `bwape.h`, `bwape.cpp` | `bwape.c`: the pairing half of sampe (`bwa_init_pe_opt`, `infer_isize`, `pairing`, `bwa_cal_pac_pos_pe`, `bwa_sw_core`, and `bwa_paired_sw`'s body as `bwa_paired_sw_core`/`bwa_paired_sw_report`) |
| `aln_out.h`, `aln_out.cpp`, `aln_main.h`, `aln_main.cpp` | new code (MIT, as the rest of bwa-mem3; `SPDX-License-Identifier` line): the SAM/BAM sink and the `aln` subcommand driver (`bwa_aln` and `bwa_sai2sam_se` option handling, fused) |
| `aln_alloc.h` | new code (MIT): `aln_malloc`/`aln_calloc`/`aln_realloc`/`aln_strdup`, which exit with a message on allocation failure; every allocation in the files above goes through them (upstream dereferences the result unchecked) |
| `aln_bwt.h`, `aln_bwt.cpp` | bwa's FM-index API (`bwt_occ`, `bwt_2occ`, `bwt_2occ4`, `bwt_sa`) re-implemented over bwa-mem3's own index (`FMI_search`); `bwt_match_exact_alt` is copied from upstream `bwt.c` |

## Modifications

Each change to upstream-derived code adds a line here.

- `aln_bwt.cpp`: `bwt_match_exact_alt` is upstream's body verbatim, moved into
  `namespace bwa_aln`. The occ and SA functions are new code with bwa's
  semantics, computed from bwa-mem3's checkpointed occ table instead of bwa's
  interleaved 2-bit `.bwt` layout.
- `bwtaln.h`: upstream's types and constants unchanged, moved into
  `namespace bwa_aln` (the `bwt_t`/`bwtint_t` types come from `aln_bwt.h`). Dropped:
  `SAI_MAGIC`, the `extern "C"` function declarations and `bwa_bam_open`.
- `bwaseqio.cpp`/`bwaseqio.h`: `bwa_seq_open`, `bwa_seq_close`, `seq_reverse`,
  `bwa_trim_read`, `bwa_read_seq` and `bwa_free_read_seq` are upstream's bodies, in
  `namespace bwa_aln`. Dropped: BAM input (`bwa_bam_open`, `bwa_read_bam`, the `is_bam`
  branches and `bamlite.h`). Input is opened with `xzopen` (zlib `gzopen`, or
  `gzdopen` on stdin for `-`), as bwa does. `nst_nt4_table` is upstream's table, kept
  local because bwa-mem3's differs at `'-'` (5 upstream, 4 here). Parsing uses
  bwa-mem3's `kseq.h`, whose record handling matches upstream's; its differences are
  an OOM guard, a wider `kseq_read` return type and extra sequence-buffer slack.
  `bwa_seq_close` reports a failed `gzclose` and returns -1 instead of exiting (upstream's
  `err_gzclose`), so the driver still completes the output (a BAM gets its EOF block) before
  exiting 1. `bwa_read_seq` still ends a batch at a malformed record (`kseq_read` -2) as
  upstream does, so the batches and the reads are unchanged, but it names the record on stderr
  and counts it (`bwa_seq_n_bad`), and the driver then exits 1 where stock exits 0.
- `bwtgap.cpp`/`bwtgap.h`: `gap_init_stack*`, `gap_destroy_stack`, `bwt_match_gap` and the
  stack helpers are upstream's bodies, in `namespace bwa_aln`. `bwt_match_gap` has no production
  caller: the search runs `gap_job_*` (below), and `bwt_match_gap` is kept as the in-tree
  reference the unit tests check `gap_job_run` against. Dropped: the `extern "C"`
  block, the malloc wrappers and the `bwa_aln2seq` declaration (ported with the SAM
  code). Added to `bwtgap.h`: declarations for `gap_init_opt`, `bwa_cal_maxdiff` and
  `bwt_cal_width` (upstream declares them in `bwtaln.h`), plus `aln_search_batch`,
  `aln_batch_max_gapo` and `width_walk_t`/`width_walk_init`/`width_walk_run`, which are new.
- `bwtaln.cpp`: `gap_init_opt`, `bwa_cal_maxdiff` and `bwt_cal_width` are upstream's
  bodies. `bwa_cal_sa_reg_gap` differs in three ways. (1) The read is complemented into
  a per-thread scratch buffer instead of in place, and `name`, `seq`, `rseq` and `qual`
  are no longer freed, so the SAM stage can use the same reads. (2) The batch-wide
  `max_len`, the `max_diff` derived from it and the `max_gapo` cap are computed once in
  `aln_search_batch` before the threads start (upstream repeats the same computation in
  every thread) and passed down; the per-read `max_diff` is recomputed per read as
  upstream does. (3) The thread count is a parameter rather than `opt->n_threads`, and
  `i % n_threads == tid` is applied for one thread too, so a count of 0 cannot divide by
  zero. Dropped: `bwa_aln_core` (the `.sai` writer and read loop), `bwa_open_reads`,
  `bwa_aln` (option parsing and usage) and the `HAVE_PTHREAD` conditionals (pthreads
  are always used). Every object here is built with `-ffp-contract=off` (Makefile) so the
  floating-point code in `bwape.cpp` (named below) rounds as a stock build without FP
  contraction does; `bwa_cal_maxdiff` itself has no multiply-add to fuse.
- `bwase.cpp`/`bwase.h`: `bwa_aln2seq_core`, `bwa_aln2seq`, `bwa_approx_mapQ`, `bwa_sa2pos`,
  `bwa_cal_pac_pos_core`, `bwa_refine_gapped_core`, `bwa_cal_md1`, `bwa_correct_trimmed`,
  `pos_end`, `bwa_print_sam1` and `bwase_initialize` are upstream's bodies, in
  `namespace bwa_aln`, with these changes: `bwa_print_sam1` and the header writer append to
  a `kstring_t` instead of writing to stdout (`bwa_print_seq` likewise); `bwa_cal_pac_pos`
  takes the `bwt_t` as a parameter instead of reloading `.bwt` and `.sa` from disk, and
  `bwa_refine_gapped` always receives the packed reference (upstream's `NULL` branch that
  reads `.pac` is dropped); the `.sai` reading, `bwa_sai2sam_se_core` and the samse command
  line are not ported (alignments come from the search's `aln` arrays, and the driver is
  `aln_main.cpp`); `bwa_print_sam1`, `bwa_cal_md1` and the header writer format numbers and strings with `kputs`/`kputc`/`kputl` instead of `ksprintf` (same text; libc's `vfprintf` takes a process-wide locale lock that serialises the formatting threads on macOS); and `bwa_print_sam1` / `bwa_refine_gapped_core` have C++ casts. Output
  goes to a sink that writes per batch, in read order.
- `bwase.cpp`: `bns_cnt_ambi` is upstream `bntseq.c`'s body, kept as the file-local
  `aln_bns_cnt_ambi`. bwa-mem3's own `bns_cnt_ambi` counts every overlapping N run (upstream
  stops at the first), which would change the `XN` tag. `bns_get_seq` is not copied:
  `bwa_refine_gapped_core` reads its reference window with bwa-mem3's `bns_get_seq_into`, which
  has upstream's clamp, swap and bridge rules and fills a caller buffer.
- `aln_print_sam_hdr` follows bwa **0.7.19's** header output, not the Apache2 branch's
  `bwa_print_sam_hdr`. The branch prints `@SQ` before `@HD` and finds `@HD` anywhere in the
  `-r` text; the release prints `@HD` first (unless a line of the `-r` text starts `@HD<tab>`),
  then `@SQ`, then the `-r` text. The release is the parity target, so its order is used. The
  `@PG` line is bwa-mem3's own (same builder as `mem`); parity checks drop `@PG` lines.
- `bwa-mem3 aln` differs from running `bwa aln` then `bwa samse` in its command line: it takes
  `aln`'s options, `-f` and `-r` (samse's and sampe's), and `--max-xa-hits` for samse's and
  sampe's `-n`. `-b`, `-0`, `-1` and `-2` (BAM input) are not accepted. `-r` with a literal tab
  is rejected as stock does (bwa-mem3's `bwa_set_rg` alone would accept it). `--bam` and
  `--bam-level` are new; `--bam` compresses on `-t` BGZF threads (the same bytes). Numbers are
  parsed strictly (stock's `atoi` turns `abc` into 0), and values stock accepts but cannot run
  are refused before the index loads: a negative `-l` (stock writes before its seed-width
  buffer, or the port's `resize` threw), `-m`, `-M`, `-O`, `-E`, `-B` (every read dropped) or
  integer `-n`, and penalties with `-o` that leave fewer priority stacks than the scores pushed
  (`-M 0 -O 0 -E 0`, or `-o` below -1 cancelling the others). `aln ref - -` and an `-f` that
  names an input are refused too. None of these changes a run stock completes.
- `srand48(bns->seed)` is called once per run, before the first batch, and `bwa_aln2seq_core` runs
  serially in read order, so the random choice among equally good hits is the same as stock's
  and does not depend on `-t`.
- `bwape.cpp`/`bwape.h`: `bwa_init_pe_opt`, `infer_isize`, `pairing`, `bwa_sw_core` and `bwa_paired_sw` are
  upstream's bodies in `namespace bwa_aln`; `bwa_cal_pac_pos_pe` is upstream's with these changes. Each read's
  alignments come from its own `aln` / `n_aln` (upstream reads them from the `.sai` and keeps a backup copy per
  read, which is not needed here); the BWT is a parameter; `infer_isize` is exported (the unit tests call it).
  Dropped: the dead `if (0 && ...)` hash-table branch and `g_hash` (so `-P` has no effect), the `.sai` reading
  and the sampe command line, the per-step timing lines on stderr, and the branch that loads the reference when
  the caller gives none (it is always the in-memory packed copy). `g_log_n` (bwase.cpp) is now shared with it. The
  mate alignment in `bwa_sw_core` calls bwa-mem3's `ksw_align` and `ksw_global`, which match stock's up to cosmetics;
  `aln_parity.sh` is the check. `infer_isize`'s `erfc`/`sqrt` loop, the pairing score's `log`/`erfc` and the
  `int += double` arithmetic in `bwa_paired_sw_core` are kept verbatim; the objects are built with `-ffp-contract=off`.
- Paired-end driver (`aln_run_pe`, new code). `bwa-mem3 aln [options] <prefix> <R1.fq> <R2.fq>` is `aln` on each
  file followed by `sampe`. Per batch of 0x40000 reads from each file: R1's reads are read, then R2's (stock
  order); each file is searched as its own batch, because stock runs `aln` on each file separately and the
  batch-wide gap-open cap depends on the longest read of the batch; `bwa_pe_aln2seq`, the one stage that draws
  from `drand48`, runs serially in read order; the insert size is estimated (a barrier), with the previous
  batch's estimate carried over when this batch gives none; then pairing (which breaks ties with a hash, not
  the RNG), the mate Smith-Waterman, the gapped refinement and the output run in parallel, each pair printed
  as read 1 then read 2 (see the stage list below). sampe's options are `-a`, `-s`,
  `-A`, `-c` and `-P` (as stock), `--max-occ` (sampe `-o`), `--max-xa-hits` (sampe `-n`, shared with the
  single-end XA limit) and `--max-xa-hits-discordant` (sampe `-N`). Where stock sampe takes the options saved in
  R2's `.sai` for both ends, one option set serves both files here. When the two files give different read
  counts in a batch, sampe's rule is kept (the count of R2) and a warning is written to stderr; if R2 has more
  reads than R1 the count of R1 is used, where sampe would read past R1's reads. A pair whose names differ after
  the `/1` and `/2` trim stops the run with an error after that pair is written, as sampe does.
- The batch-wide `max_gapo` cap does not change any output, so searching the two files of a pair separately or
  together gives the same hits and the same `-m` behaviour. A gap open is tried only from an entry with at least one
  difference still allowed, so its gap opens stay below `min(-o, differences allowed for the file's longest read)`,
  and a larger cap never lets a further entry in. The driver still searches each file as its own batch, as stock does.
- `bwase.cpp`: `aln_print_sam_hdr` is new code that follows the 0.7.19 release's header
  order; its line-start counting is a small helper (`count_line_start`).
- `aln_main.cpp`: new code. Inputs and `-f` are checked before the index loads, so a failed
  run leaves no output; help text lists each option with its stock default and source command.
- `bwape.cpp`/`bwape.h`: `bwa_cal_pac_pos_pe` is cut into its four stages, `bwa_pe_aln2seq` (the
  `drand48` draws), `bwa_pe_se_pos`, `bwa_pe_isize` (`infer_isize` and the `last_ii` fallback) and
  `bwa_pe_pair`, each holding upstream's loop body unchanged; `bwa_cal_pac_pos_pe` calls them in turn
  (the driver calls the stages itself; `bwa_cal_pac_pos_pe` is kept for the unit tests).
  `bwa_paired_sw` is likewise split into `bwa_paired_sw_core` (the loop, adding to the counts) and
  `bwa_paired_sw_report` (the two stderr lines); the wrapper itself is not kept. The `bwa_aln2seq_core` calls in `bwa_pe_pair` use
  `set_main = 0` and never reach the sampling branch (`rest` there equals the hits still to place, so
  no interval exceeds it), so that stage draws nothing.
- The driver runs the stages after the search on `kt_for` (stage order below). Single end: the `drand48`
  pass (`bwa_aln2seq_core`) runs serially over the batch in read order, then position, refinement, MD,
  trim correction and SAM text run in parallel on slices of the batch (`bwa_cal_pac_pos`,
  `bwa_refine_gapped` and `bwa_print_sam1` unchanged, on each slice). Paired end: `bwa_pe_aln2seq`
  serially in stock order (pair by pair, read 1 then read 2), `bwa_pe_se_pos` in parallel,
  `bwa_pe_isize` as a barrier, then pairing, mate Smith-Waterman, refinement and each pair's records in
  parallel. Each slice's text goes to its thread's buffer, and the slices are written in read order; a
  pair with mismatched names still ends the output after that pair.
- Interleaved search. `bwtgap.cpp` gains `gap_job_start` and `gap_job_run` (new code, appended; the
  upstream functions above them are unchanged): `bwt_match_gap`'s loop for one read as a resumable job,
  which stops just before each `bwt_2occ4` lookup after prefetching the two Occ blocks the lookup
  reads (`bwt_prefetch_occ`, new in `aln_bwt.cpp`). The loop body is upstream's, statement for
  statement, with its locals moved into the job. In `bwtaln.cpp`, `bwa_cal_sa_reg_gap` is replaced by
  `bwa_cal_sa_reg_gap_group`, which prepares each read as upstream does (widths, per-read `max_diff`
  and `seed_len`, complemented copy) and keeps 8 of the thread's reads in flight, each with its own
  stack and widths, stepping them in turn so several Occ misses overlap. Each read's hits are exactly
  what `bwt_match_gap` returns for it; `aln_sai_parity.sh` is the check. Reads are still dealt to
  threads round-robin, on fresh pthreads per batch as upstream does (moving the search onto `kt_for_nt`
  would need per-thread slot state and a fresh throughput A/B; see the comment on `aln_search_batch`). A
  slot whose stack grew past 65,536 entries frees its stack buffers after that read, so each thread keeps
  at most 8 x 65,536 (`ALN_STACK_KEEP`) stack entries across its slots. A slot's seed-width buffer is
  sized on the first read longer than `-l` (upstream's single per-thread `calloc` of `-l` + 1 entries
  is never touched for shorter reads; zero-filling one per slot up front cost 8 x 16 bytes x `-l` per
  thread, gigabytes for the large `-l` used to turn seeding off). The group depth of 8 was measured on
  one host (Apple M2 Max); the hits do not depend on it.
- Failure handling (new code, no output change). Allocations go through `aln_alloc.h` instead of the
  bare libc calls. The three `assert`s in `bwa_refine_gapped_core` (window inside the packed reference,
  window length, a CIGAR was produced) are `xassert`, so they hold in a build with `NDEBUG` defined too.
  `aln_search_batch` runs a thread's share of the reads itself if `pthread_create` fails, instead of
  leaving those reads unsearched.
- Interleaved exact-match walks (new code). When an entry has no differences left, `gap_job_run` runs
  `bwt_match_exact_alt`'s loop itself, one base per call: it prefetches the two Occ blocks of the
  next `bwt_2occ` and yields, and the next call does that lookup and updates `k`/`l` as
  `bwt_match_exact_alt` does, with its early exit on `k > l` (its N check is kept but cannot fire
  here: an N forces a width restart, so no entry reaches the tail with an N ahead). The hit, if any, is then
  handled by the same code as before (moved into `gap_job_hit`). `bwt_match_exact_alt` and
  `bwt_match_gap` themselves are unchanged. Likewise, `bwa_cal_sa_reg_gap_group` computes a read's
  two widths (whole read and seed) with `width_walk_run`, `bwt_cal_width`'s loop as a resumable walk
  that yields before each `bwt_2occ`, so they overlap the other slots' lookups; `bwt_cal_width` is
  unchanged and still used by the unit tests as the reference. Per read, the same lookups are done:
  the search's in the same order as before, while the whole-read and seed width walks, which share
  no state, now take turns. Only the interleaving changes; the output does not.
