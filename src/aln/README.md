# src/aln: bwa `aln`/`samse`/`sampe` port

This directory holds the code behind `bwa-mem3 aln`, a fused port of bwa's
short-read `aln` + `samse`/`sampe` pipeline whose SAM records without `--meth` or
`--short-insert` are byte-identical to
stock bwa 0.7.19's (the header is too, except the `@PG` line, which names bwa-mem3), as
checked on the generated reads described in `docs/src/cli/aln.md`.

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
| `bwtgap.h`, `bwtgap.cpp` | `bwtgap.h`/`bwtgap.c`: the backtracking search (`bwt_match_gap` and its priority stack), plus the chunked stack `gap_cstack_t`/`gap_pool_t`, `width_walk_init`/`width_walk_run` and `gap_slots_run` (new). The interleaved search (the job functions, the stack helpers, `gap_slots_run`) is a template over a policy: `AlnStock` is the stock search; `bash test/tools/disasm_identity.sh BASE NEW` compares its instructions against the pre-template build (a manual pre-release check; CI runs only its negative controls, `test/tools/disasm_identity_test.sh`); the batch driver's group loop, thread fan-out and duplicate compaction take the same policy (`gap_group<P>` in `bwtaln.cpp`); `AlnMeth` is `aln --meth`. Both policies use stock-sized entries: 511 per chunk; typed capacity derives from the actual stock chunk allocation |
| `bwtaln.cpp` | the search half of `bwtaln.c`: `gap_init_opt`, `bwa_cal_maxdiff`, `bwt_cal_width`, `bwa_cal_sa_reg_gap` and the thread fan-out |
| `bwase.h`, `bwase.cpp` | `bwase.h`/`bwase.c`: the SAM half of samse (`bwa_aln2seq_core`, `bwa_sa2pos`, `bwa_cal_pac_pos_core`, `bwa_refine_gapped`, `bwa_print_sam1`), plus `bns_cnt_ambi` from `bntseq.c`; `bwa_pos_from_sa` (new) is `bwa_sa2pos`'s tail, split out |
| `bwape.h`, `bwape.cpp` | `bwape.c`: the pairing half of sampe (`bwa_init_pe_opt`, `infer_isize`, `pairing`, `bwa_cal_pac_pos_pe`, `bwa_sw_core`, and `bwa_paired_sw`'s body as `bwa_paired_sw_core`/`bwa_paired_sw_report`) |
| `aln_out.h`, `aln_out.cpp`, `aln_main.h`, `aln_main.cpp` | new code (MIT, as the rest of bwa-mem3; `SPDX-License-Identifier` line): the SAM/BAM sink and the `aln` subcommand driver (`bwa_aln` and `bwa_sai2sam_se` option handling, fused) |
| `aln_meth.h`, `aln_meth.cpp`, `aln_meth_search.h`, `aln_meth_realign.h`, `aln_meth_realign.cpp` | new code (MIT): converted candidate-window discovery, owned fixed-reward global CIGARs, stock chemistry-aware cost/count accounting, and methylation single/paired output. See the ranking and resource limits below |
| `aln_alloc.h` | new code (MIT): `aln_malloc`/`aln_calloc`/`aln_realloc`/`aln_strdup`, which exit with a message on allocation failure; every allocation in the files above goes through them (upstream dereferences the result unchecked) |
| `aln_bwt.h`, `aln_bwt.cpp` | bwa's FM-index API (`bwt_occ`, `bwt_2occ`, `bwt_2occ4`, `bwt_sa`, plus the new `bwt_sa_batch`) re-implemented over bwa-mem3's own index (`FMI_search`); `bwt_match_exact_alt` is copied from upstream `bwt.c` |

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
  `aln_batch_max_gapo`, `read_claim_t`/`ALN_CLAIM_CHUNK`, `ALN_GROUP`, `aln_dup_state_t`, the
  `ALN_DUP_*` constants, `detail::find_dup_reads` and the `detail::dup_sample_*` checks,
  `gap_job_t` with its option copy `gap_job_const_t`,
  `width_walk_t`/`width_walk_init`/`width_walk_run`, `gap_slot_core_t`/`GAP_SLOT_*`/`gap_slots_run`,
  the chunked stack
  (`gap_chunk_t`, `gap_pool_t`, `gap_bucket_t`, `gap_cstack_t` and their init, release and destroy
  functions, `ALN_POOL_KEEP`) and the out-of-line push, pop and reset wrappers the unit tests use
  (`gap_stack_push_ext`, `gap_stack_pop_ext`, `gap_stack_reset_ext`, `gap_cstack_push_ext`,
  `gap_cstack_pop_ext`), which are new. Score arithmetic is widened, the optimistic stack bound
  of both stacks saturates at 2048 (the 11-bit packed score capacity), and an actual score
  outside the stack or packed range exits 1 with a diagnostic before indexing or packing. A
  large optimistic bound alone does not reject a search whose actual scores fit.
- `bwtaln.cpp`: `gap_init_opt`, `bwa_cal_maxdiff` and `bwt_cal_width` are upstream's
  bodies. `bwa_cal_sa_reg_gap` differs in three ways. (1) The read is complemented into
  a per-thread scratch buffer instead of in place, and `name`, `seq`, `rseq` and `qual`
  are no longer freed, so the SAM stage can use the same reads. (2) The batch-wide
  `max_len`, the `max_diff` derived from it and the `max_gapo` cap are computed once in
  `aln_search_batch` before the threads start (upstream repeats the same computation in
  every thread) and passed down; the per-read `max_diff` is recomputed per read as
  upstream does. (3) The thread count is a parameter rather than `opt->n_threads`, and the
  threads claim chunks of consecutive reads from a shared counter instead of taking
  `i % n_threads == tid` (see "Dynamic read claiming" below). Dropped: `bwa_aln_core` (the
  `.sai` writer and read loop), `bwa_open_reads`, `bwa_aln` (option parsing and usage) and the
  `HAVE_PTHREAD` conditionals (pthreads are always used). Every object here is built with
  `-ffp-contract=off` (Makefile) so the floating-point code in `bwape.cpp` (named below) rounds
  as a stock build without FP contraction does; `bwa_cal_maxdiff` itself has no multiply-add to fuse.
- `bwase.cpp`/`bwase.h`: `bwa_aln2seq_core`, `bwa_aln2seq`, `bwa_approx_mapQ`, `bwa_sa2pos` (now
  `bwt_sa` plus its split-out tail `bwa_pos_from_sa`), `bwa_cal_pac_pos_core`, `bwa_refine_gapped_core`, `bwa_cal_md1`, `bwa_correct_trimmed`,
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
  `int += double` arithmetic in `bwa_paired_sw_core` are kept verbatim on the default path; the objects are built with `-ffp-contract=off`.
- `bwtaln.h`, `bwape.cpp`/`bwape.h`: the default-zero `pe_opt_t::short_insert` enables
  `--short-insert` for ordinary and methylation-aware paired alignment. `infer_isize`
  then floors the lower interquartile fence at one base instead of the batch's longest
  read, clears failed estimates, and rejects nonfinite mean/SD before integer conversion.
  Previous-valid-estimate fallback is unchanged and `-A` still overrides estimation.
  Pairing and methylation cap rescue reject insert-likelihood calculations that underflow
  or produce unrepresentable integer penalties under this option. These opt-in changes can alter pair selection, proper-pair flags,
  MAPQ and rescue windows; pairing geometry and the rescue algorithm are unchanged.
  Without the option, the stock inference and scoring behavior is preserved.
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
- `aln_main.cpp`: new code. Inputs and existing `-f` paths are checked before the index loads.
  Absent output paths are left untouched until the real open after loading the index, which
  reports creation failures. Help text lists each option with its stock default and source command.
- `bwape.cpp`/`bwape.h`: `bwa_cal_pac_pos_pe` is cut into its four stages, `bwa_pe_aln2seq` (the
  `drand48` draws), `bwa_pe_se_pos`, `bwa_pe_isize` (`infer_isize` and the `last_ii` fallback) and
  `bwa_pe_pair`, each holding upstream's loop body unchanged apart from the batched SA lookups (see
  "Batched SA lookups" below); `bwa_cal_pac_pos_pe` calls them in turn
  (the driver calls the stages itself; `bwa_cal_pac_pos_pe` is kept for the unit tests).
  `bwa_paired_sw` is likewise split into `bwa_paired_sw_core` (the loop, adding to the counts) and
  `bwa_paired_sw_report` (the two stderr lines); the wrapper itself is not kept. The `bwa_aln2seq_core` calls in `bwa_pe_pair` use
  `set_main = 0` and never reach the sampling branch (`rest` there equals the hits still to place, so
  no interval exceeds it), so that stage draws nothing.
- The driver runs the stages after the search on `kt_for` (stage order below). Single end: the `drand48`
  pass (`bwa_aln2seq_core`) runs serially over the batch in read order, then position, refinement, MD,
  trim correction and SAM text run in parallel on slices of the batch (`bwa_cal_pac_pos`,
  `bwa_refine_gapped` and `bwa_print_sam1` unchanged apart from `bwa_cal_pac_pos`'s batched SA
  lookups, on each slice). Paired end: `bwa_pe_aln2seq`
  serially in stock order (pair by pair, read 1 then read 2), `bwa_pe_se_pos` in parallel,
  `bwa_pe_isize` as a barrier, then pairing, mate Smith-Waterman, refinement and each pair's records in
  parallel. Each slice's text goes to its thread's buffer, and the slices are written in read order; a
  pair with mismatched names still ends the output after that pair.
- Interleaved search. `bwtgap.cpp` gains `gap_job_start` and `gap_job_run` (new code, appended; the
  upstream functions above them are unchanged): `bwt_match_gap`'s loop for one read as a resumable job,
  which stops just before each `bwt_2occ4` lookup after prefetching the two Occ blocks the lookup
  reads (`bwt_prefetch_occ`, new; now inline in `aln_bwt.h`). The loop body is upstream's, statement for
  statement, with its locals moved into the job. In `bwtaln.cpp`, `bwa_cal_sa_reg_gap` is replaced by
  `bwa_cal_sa_reg_gap_group`, which prepares each read as upstream does (widths, per-read
  `max_diff` and `seed_len`, complemented copy) and keeps `ALN_GROUP` of the thread's reads in
  flight, each with its own stack and widths, stepping them in turn so several Occ misses overlap.
  Each read's hits are exactly what `bwt_match_gap` returns for it; `aln_sai_parity.sh` is the
  check. The search still runs on fresh pthreads per batch as upstream does (moving it onto
  `kt_for_nt` would need per-thread slot state and a fresh throughput A/B; see the comment on
  `search_reads`). Each slot's stack is the chunked stack described below. A
  slot's seed-width buffer is sized on the first read longer than `-l` (upstream's single per-thread
  `calloc` of `-l` + 1 entries is never touched for shorter reads; zero-filling one per slot up
  front cost `ALN_GROUP` x 16 bytes x `-l` per thread, gigabytes for the large `-l` used to turn
  seeding off). The hits do not depend on the group depth (`ALN_GROUP`, see the chunked stack
  below).
- Failure handling (new code, no output change). Allocations go through `aln_alloc.h` instead of the
  bare libc calls. The three `assert`s in `bwa_refine_gapped_core` (window inside the packed reference,
  window length, a CIGAR was produced) are `xassert`, so they hold in a build with `NDEBUG` defined too.
  If `pthread_create` fails, the search (`search_reads`) runs a worker on the calling thread in
  place of the missing threads; it claims reads from the shared counter like the others, so no read
  is left unsearched.
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
- Chunked search stack (new code). The interleaved search no longer uses upstream's `gap_stack_t`:
  `gap_cstack_t` (`bwtgap.cpp`) has the same score buckets, the same `best` and entry counts and the
  same push and pop statements (apart from the whole-entry store below), but each bucket is a chain
  of fixed-size chunks (`GAP_CHUNK_BYTES`, 16 KiB, 511 entries) rather than an array that doubles on
  its own. A push that finds its bucket's top chunk full (or no chunk) takes one from a per-thread
  pool; a pop that empties one gives it back; a finished read gives back all of its stack's chunks.
  So the memory the stacks hold follows the entries in use rather than every bucket's high-water
  mark, which lets `bwa_cal_sa_reg_gap_group` keep more reads in flight per thread for less peak
  memory than 8 took with the array stacks. `ALN_GROUP` is 12 on arm64, where it measured faster
  than 8 on Apple M3 Ultra (4-11% wall), Graviton4 and Graviton3; on the Gravitons, against 8 with
  the same chunked stacks, it held up to about 0.1 GiB more at peak at `-t 16` on deep searches.
  Other arm64 cores were not measured and inherit 12. It is 8 on x86-64, where 12 measured slower,
  and on any other target (see the comment at `ALN_GROUP` in `bwtgap.h`; a unit test pins the value
  per target). The pool keeps at most 256 free chunks (`ALN_POOL_KEEP`, 4 MiB) per thread and frees
  the rest; a 1,024-chunk pool measured no faster and held about 0.2 GiB more at peak. Within a
  bucket the last entry pushed is popped first, as in the array, and the bucket popped is chosen by
  the same code, so entries come off in the same order and each read's hits are unchanged.
  `bwt_match_gap` keeps upstream's `gap_stack_t` and is the unit tests' reference.
- Whole-entry push (new code). `gap_cpush` builds the entry in a local `gap_entry_t`, sets every field
  (`n_seed_mm` to 0) and stores it in one piece. Upstream's `gap_push` assigns `n_mm`, `n_gapo`,
  `n_gape` and `state` into the slot's bitfield word, which the compiler turns into a load of that
  word to keep `n_seed_mm`; the slot is often a line not yet in cache, so the load stalled. Nothing
  reads `n_seed_mm`, so the entries popped are the same. `gap_push` is unchanged and still used by
  `bwt_match_gap`.
- Top-of-stack forwarding (new code). After expanding an entry, the job step stores each child but
  the last as before and keeps the last in registers. If its score is at most the stack's best bucket
  once the others are stored, it is the entry the next pop would return (the last entry of the lowest
  non-empty bucket), so the loop continues with it directly; otherwise it is stored and the pop runs.
  The loop top counts a forwarded entry as on the stack (the empty-stack test, `max_entries` and the
  `-m` break see the count a push would give), and a push followed by that pop would leave every
  bucket count, the total and `best` as they are, so the entries searched, their order and every
  break are unchanged.
- Hot-loop layout (new code). `bwt_occ`, `bwt_2occ`, `bwt_2occ4` and `bwt_prefetch_occ` are defined
  inline in `aln_bwt.h` (same bodies), reading the occ table through a pointer cached in `bwt_t`, so
  the job step and the width walks no longer call into `aln_bwt.cpp` for them (`bwt_2occ4` is left
  to the compiler's inlining choice: forcing it inline measured slower on x86-64 and arm64). In
  `bwtgap.cpp` the search's stack helpers are forced inline and its rare paths (a found hit, a stack
  chunk taken from or given back to the pool) are kept out of line. Code placement only; the
  statements run are the same.
- Job layout (new code). `gap_job_t` puts the fields the job step reads on every call first (resume
  state, inputs, then `c`, a copy of the options it reads, made by `gap_job_start`), and two
  expressions are precomputed: `seed_off = len - seed_len` per read, and `score_cut`, the score above
  which the search stops (`best_score + s_mm`, updated where `best_score` is, or the largest packed
  score under `NONSTOP`, where upstream skips the test). The cut is computed in 64 bits and clamped
  to the packed score range, which leaves every comparison with a packed score unchanged.
- Slot loop (new code). `gap_slots_run` (`bwtgap.cpp`) is the loop of `bwa_cal_sa_reg_gap_group` over
  its slots, with the job step (`gap_job_step`, also behind `gap_job_run` for the unit tests) and the
  width-walk step inlined, so stepping a slot is no longer a call; it returns to the driver only
  when a read's widths are done (the driver starts its search) or its search is done (the driver
  takes the hits and refills the slot). `width_walk_init`/`width_walk_run` move from `bwtaln.cpp` to `bwtgap.cpp` for
  this. Only the interleaving of the reads changes; each read's steps and hits are the same.
- Dynamic read claiming (new code). Upstream deals the reads to the threads round-robin (thread `t`
  searches reads `t`, `t + n_threads`, ...); reads differ widely in cost, so threads finished the
  batch unevenly. Each thread now claims `ALN_CLAIM_CHUNK` (64) consecutive reads at a time from a
  shared atomic counter (`read_claim_t`) and keeps its `ALN_GROUP` slots fed across chunks. A read's
  hits do not depend on which thread searches it, so the output is the same at any thread count.
- Duplicate reads (new code, `bwtaln.cpp`). A read's hits depend only on its bases, its length and
  the batch's options (the per-read `max_diff` and `seed_len` follow the length; the batch-wide
  `max_gapo` cap is computed from every read first), so a batch with duplicates is searched once per
  distinct read. `detail::find_dup_reads` points every read at the first read of the batch with the
  same length and bases: a fingerprint of the length and every base, read 8 at a time as 64-bit
  words (computed on the thread pool, `kt_for_nt`), then one table pass on the calling thread with
  a full compare. The distinct reads are moved to the front of the batch (the duplicates set
  aside), searched there by the unchanged search, and every read moved back to its index; every
  later copy gets its own copy of the first one's `aln`. Duplicates are found within one batch
  only. Everything after the search, the `drand48` pass included, still runs per read, so the
  output is unchanged. A batch without duplicates runs the same search code as before; what it
  pays is the check. `detail::dup_sample_check` first samples the batch twice: about 1 read in 32
  chosen by a key of the length and three 8-base words (start, middle, end), so that every copy of
  a read is in or out together, and every 32nd read by position, which catches a read with many
  copies whose key is not sampled. The full pass runs when either sample has at least 1 duplicate
  per 200 reads (`ALN_DUP_MIN_RATIO`), when the batch is under `ALN_DUP_SAMPLE_MIN` reads, or when
  the key sample is far from its expected size. A batch the samples or the full pass find with
  fewer duplicates than that is searched read by read, and turns the check off for the next 7
  batches of its stream (`ALN_DUP_RECHECK`; `aln_dup_state_t`, which the driver keeps once for
  single-end and once per file for paired-end, since one mate's batches need not look like the
  other's); a read searched rather than copied gets the same hits. The full pass holds 8 bytes of
  fingerprint and 8 of indexes per read, a table of 16 to 32 bytes per read, and the duplicates'
  set-aside entries (about 200 bytes each), all freed after the batch; peak RSS measured 30-40 MB
  higher on inputs with many duplicates. `aln_run_se`/`aln_run_pe` report at `-v
  4` how many reads of each batch shared a search; `share_dup_search` (a test hook) turns the
  sharing off.
- Read-ahead (new code, `aln_main.cpp`). `aln_run_se` and `aln_run_pe` read and parse the next
  batch (for paired-end, the next R1 batch and then, if it is non-empty, the next R2 batch) on a
  helper thread, one thread on top of `-t`, while the current batch is searched, paired and
  written, so the threads no longer wait on the input between batches. The batches are processed
  and written in input order, the `drand48` pass and the insert-size carry-over stay serial and in
  order, and each stream is read by one thread at a time, so the output is unchanged; a run that
  stops early (pair names that differ, one file ending first, a write error) stops at the same
  record. A gzip read error of the input does not exit on the helper (stock's `err_gzread` exits
  inside the read): the reader keeps it (`bwa_seq_defer_read_errors`, `bwaseqio.cpp`), and after
  writing the current batch, where stock's read of the next batch would have failed, the driver
  prints the same message and returns 1 (stock exits 1 there; returning lets the caller close the
  output, so a BAM is terminated, as for a failed close of the input). The records written are
  stock's. One extra batch is held in memory. On stderr, messages printed while reading
  (malformed records, trimming) can come before the previous batch's progress line, and after an
  early stop the batch already read ahead (which is never aligned) has been parsed, so its
  messages are printed too, and the run exits only once that read is done (a blocked read cannot
  be interrupted safely; if the helper thread could not be created, that read is skipped instead).
  An allocation failure while reading still exits on the helper thread.
- Batched SA lookups (new code). `bwt_sa_batch` (`aln_bwt.cpp`) gives `bwt_sa` for an array of rows,
  walking up to `SA_RESOLVE_LANES_DEFAULT` (`FMI_search.h`) rows at once with each lane's next Occ
  block or SA sample prefetched, as mem's SA resolver does; each result depends on its row alone
  (row 0 gives `(bwtint_t)-1`, as `bwt_sa`). `bwa_sa2pos` is split into `bwt_sa` and its tail,
  `bwa_pos_from_sa`, unchanged. `bwa_cal_pac_pos` and `bwa_pe_se_pos` collect the rows their loop
  looks up, in the loop's order, resolve them in one call, then run upstream's loop taking each
  result in turn. `bwa_pe_pair` takes the pairs in groups of consecutive pairs whose weight
  (pairing rows, plus an upper bound on the extra hits their bodies can generate, plus one per
  pair) is at most `ALN_PE_GROUP_ROWS` (8,192); a heavier pair is a group alone, so each group
  buffer holds at most that many entries or one pair's rows, which upstream holds for that pair
  too. Per group it resolves the group's pairing rows, runs upstream's body pair by pair up to the
  extra hits' lookups (each pair's array is built in upstream's push order, so `pairing` sees the
  same input), then resolves the group's extra-hit rows and finishes each pair. A pair over
  `max_occ` still skips its extra hits. A pair's body reads and writes only that pair's reads, and
  no stage here draws from `drand48`, so the output is the same as upstream's single pass.
- Pairing sort (new code). `pairing` sorts its (position, k<<2 | strand<<1 | end) rows with
  `pdqsort_128` in place of `ks_introsort_128`. Both order rows by the whole (x, y), so rows equal
  under that order are identical bytes and the sorted array is the same whichever sort makes it.

## Methylation candidate ranking and stock accounting

Under `aln --meth`, converted-space search generates reference windows within the user's discovery limits. Each exact original-contig/hypothesis/strand/start/end window is refined once with the existing chemistry matrix and global reward alignment (match +1, mismatch -3, gap open 5 plus gap length; neutral conversions score 0). The raw chosen CIGAR then supplies chemistry-aware mismatches, gap opens and gap extensions (length minus one per gap run). Candidate cost is `M*mismatches + O*opens + E*extensions`, using `-M/-O/-E`; accepted conversions, including neutral conversions, are not mismatches. At alignment time original `.amb` intervals mask packed random fills to N before refinement/accounting, without rebuilding indexes. Refinement retains the ambiguity reward -1; raw counts and NM count either reference or query N as a mismatch even for N/N. MD reports original N and XM makes no call there. These rules are specific to `aln --meth`.

Different spans compete before collapsing a final locus (original contig, hypothesis, strand, normalized start). Each locus retains its minimum-cost refined alignment, with deterministic CIGAR/end ties. X0 counts all best-cost loci. X1 counts all retained worse-cost loci through best cost plus `-M`, or all returned worse costs under `-N`, not just the next distinct cost. Primary selection, XA and stock `bwa_approx_mapQ` use these counts/costs, with the stock effective mismatch ceiling. Same-locus alternative windows do not inflate placement counts.

This is cost-ranking of reward-refined candidate CIGARs, not minimum-cost optimization over every CIGAR or a globally exhaustive search. Stock aln also separates search-cost selection from fixed-reward output refinement, but refines after selection; meth evaluates each candidate before selection. `-R` and chemistry-dependent width shadows do not prune original-reference competitors. For an unseeded result with current best cost `B`, let `T=B+M` and `k=floor(T/M)`. When `M>0`, `T<O`, `E>=0`, `max_diff>=k`, `-N` is off and no discovery cap has bound, stopping is safe after flushing every converted score through `M*k`: every retained raw alignment is ungapped with at most `k` strict mismatches and has an equivalent proposal in that stratum. A decreasing best cost only narrows the required stratum. Other cases keep configured generation without a path-cost bound. Primary and XA output use the stored winning alignment without a second refinement. Counts/cost and the stored reward include raw terminal gaps; display converts terminal insertions to soft clips and removes terminal deletions with coordinate adjustment. MD/NM/XM describe that normalized display alignment, so its NM can differ from raw ranking counts.

The total per-read SA lookup budget is `p_best + p_second` (4596 by default), not separate rank-dependent budgets. Repeated interval/span proposals and exact windows are deduplicated. Lookup or `-m` exhaustion marks a result incomplete: retain the best-known placement if one exists, classify it as repeat with MAPQ 0, and omit exact X0/X1/XA claims. Cumulative candidate CIGAR bytes are capped at 64 MiB per read, conservatively bounding retention even when a locus is replaced. A candidate that would exceed this cap stops that read's search and marks it incomplete in the same way (`ALN_METH_BIND_CIGAR`); a read with no retained placement stays unmapped, and later reads continue normally. The realignment wrapper separately guards its active allocation estimate. Invalid results and other resource failures are fatal, not skipped competitors. A refined window whose counts do not fit a record (more than 255 mismatches, gap opens or gap extensions) is dropped; when its cost is within the retained window (best cost plus `-M`, or any cost under `-N`) the read is marked incomplete as above. An XA alternative whose edit distance exceeds the legacy 8-bit field drops the XA tag for that read; the primary is unaffected.

MAPQ reuses stock aln's heuristic, not a new reward-gap probability model. Known best-cost ties have MAPQ 0; explicit lookup/work-incomplete reads always have MAPQ 0. No stronger search-completeness or calibration guarantee than stock is implied. Paired-end meth uses stock sampe cost-based pairing and MAPQ on retained original-reference loci, with selected CIGAR/count/hypothesis ownership preserved. Each end may pair on all retained distinct loci up to the fused `--max-occ` limit (sampe `-o`), provided neither end has lookup/work cutoff flags. This is completeness within configured discovery and retained stock cost strata, not exhaustive genomic discovery. Explicit incomplete ends keep best-known repeat output and are not promoted by pairing/rescue. Both paths share paired batch/error handling and the same insert-size inference policy, parameterized by original reference length. Under `--meth-context cpg`, private original-reference windows use the shared seven-symbol CpG labels with true contig flanks and original ambiguity-run masking. Candidate and rescue global refinement uses exact full-band scalar fallback, not the five-symbol certificate or SIMD global path. Shared option validation requires TAPS/5base and non-collapsed scoring. TAPS/5base neutral first compares raw-CIGAR conversion counts within equal-cost windows at a locus and then splits equal minimum-cost loci by the fewest conversions. Displaced loci move from X0 to X1. Independent per-read draws select among the remaining best loci, and stock MAPQ is capped at `min(60, 6 * nearest extra conversions)`. Unchanged ends retain the ceiling through pairing/rescue. Stock pairing-moved ends reset seQ and use pairing confidence instead. Incomplete results never split or gain confidence.

The manual stock instruction gate requires both group and batch drivers and at least four compared functions. Missing symbols, empty disassembly and tool failures are errors. Build the positive control with `make EXTRA_CXXFLAGS=-DALN_DISASM_POSITIVE_CONTROL` and require a difference. This local gate is not CI coverage or cross-platform/LTO validation.

Lookup-cap rescue is placement-only. A cutoff-free mate with at most 255 best-cost loci supplies contig-clipped windows, cached by exact identity. End-to-end raw-CIGAR cost must not exceed the end's known best cost. Candidate pair scores and neutral conversion counts select a deterministic owned alignment. Omitted genomic competitors still force the rescued end to MAPQ 0 and suppress exact strata/XA. Neither end is marked proper solely by this rescue. An unchanged unique anchor retains its confidence, including the neutral guard; moved or tied anchors have MAPQ 0. The historical ceiling of 15 is not a confidence floor and supplies no calibration evidence for this path.
