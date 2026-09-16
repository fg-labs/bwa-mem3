/* The MIT License

   Copyright (c) 2008 Genome Research Ltd (GRL).

   Permission is hereby granted, free of charge, to any person obtaining
   a copy of this software and associated documentation files (the
   "Software"), to deal in the Software without restriction, including
   without limitation the rights to use, copy, modify, merge, publish,
   distribute, sublicense, and/or sell copies of the Software, and to
   permit persons to whom the Software is furnished to do so, subject to
   the following conditions:

   The above copyright notice and this permission notice shall be
   included in all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
   BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.

   Modified Copyright (C) 2019 Intel Corporation, Heng Li.
   Contacts: Vasimuddin Md <vasimuddin.md@intel.com>; Sanchit Misra <sanchit.misra@intel.com>;
   Heng Li <hli@jimmy.harvard.edu> 
*/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <time.h>
#include <getopt.h>
#include <limits.h>
#include <stdint.h>
#include <errno.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <zlib.h>

#include "bntseq.h"
#include "bwa.h"
#include "bwt.h"
#include "utils.h"
#include "FMI_search.h"
#include "fm_index_writer.h"
#include "bwa_shm.h"
#include "kseq.h"
#include "libsais_build.h"
#include "system.h"

KSEQ_DECLARE(gzFile)

/* Writes two converted contigs per input chromosome to the seed FASTA — the
 * G→A-projected reverse-strand target (`>r<name>`) and the C→T-projected
 * forward-strand target (`>f<name>`), wrapped at 100 bp. The caller
 * (meth_index_build) then builds the seed FM-index over it. The per-strand
 * projection matches bwameth.py's `index-mem2` contig layout. */
static void meth_project_and_write(FILE *out, const char *prefix, const char *name,
                                   const char *seq, size_t len, char from, char to)
{
    char buf[65536];
    size_t bl = 0;
    fprintf(out, ">%s%s\n", prefix, name);
    for (size_t i = 0; i < len; ++i) {
        if (bl + 2 > sizeof(buf)) { fwrite(buf, 1, bl, out); bl = 0; }
        char c = seq[i];
        buf[bl++] = (c == from) ? to : c;
        if (((i + 1) % 100) == 0) buf[bl++] = '\n';
    }
    if (len % 100 != 0) {
        if (bl + 1 > sizeof(buf)) { fwrite(buf, 1, bl, out); bl = 0; }
        buf[bl++] = '\n';
    }
    if (bl) fwrite(buf, 1, bl, out);
}

/* Skip the rebuild if the c2t FASTA is already newer than the input. */
static int meth_c2t_is_fresh(const char *in_fa, const char *out_fa)
{
    struct stat a, b;
    if (stat(in_fa, &a) != 0 || stat(out_fa, &b) != 0) return 0;
    return b.st_mtime >= a.st_mtime;
}

/* D3 BS-seq index layout. `index --meth` builds TWO indexes from `fa`:
 *   1. The ORIGINAL-alphabet index `<fa>.{pac,ann,amb,bwt.2bit.64}` — real chrom
 *      names + original bases. This is the extension/scoring reference and the basis
 *      for variant-callable `mem --meth` output. (Identical to a normal `index`.)
 *      `mem` pac-fetches its bases from `.pac`, so no `.0123` is written (see
 *      bwa_idx_build's emit_unpacked_ref default).
 *   2. A converted SEED index `<fa>.meth.{pac,ann,amb,bwt.2bit.64}`, built over a
 *      per-strand-converted FASTA `<fa>.meth.fa` (two contigs per chromosome:
 *      `>r<name>` = G->A reverse-strand target, `>f<name>` = C->T forward-strand
 *      target). 3-letter seeding requires per-strand conversion because C->T and
 *      reverse-complement do not commute. This index is used ONLY to find seed
 *      SA-intervals, which are then remapped to original coordinates for chaining and
 *      extension. The `.meth` separation is by file PREFIX (not by changing the global
 *      CP_FILENAME_SUFFIX, which serves every index). */
/* Sum the reference length the way bns_fasta2bntseq will pack it, so the seed
 * index's cost can be estimated before any index is built. Costs one extra
 * sequential pass over the reference on every --meth run, including the ones
 * that go on to succeed: ~10 s for an uncompressed hg38, longer for a gzipped
 * one where the pass is decompression-bound. That is under 1% of the two SA
 * constructions it guards, and it buys a refusal that leaves no partial index
 * behind. Returns -1 on read error. */
static int64_t fasta_count_bases(const char *fa)
{
    gzFile in = gzopen(fa, "r");
    if (in == NULL) return -1;
    kseq_t *seq = kseq_init(in);
    int64_t total = 0;
    int kr = 0;
    while ((kr = kseq_read(seq)) >= 0) total += (int64_t)seq->seq.l;
    kseq_destroy(seq);
    gzclose(in);
    return kr < -1 ? -1 : total;
}

static int meth_index_build(const char *fa, int emit_unpacked_ref,
                            int64_t max_memory_bytes, int max_memory_user_specified)
{
    /* 0. Preflight the DOMINANT build before writing anything. --meth builds two
     * indexes sequentially in one process: the original over N = 2*l_pac, then a
     * seed over a per-strand-converted FASTA whose text is twice as long again
     * (4*l_pac). The seed therefore costs ~2x the original and is what decides
     * whether the invocation can run at all -- but the per-build check inside
     * libsais_build_fm_index only fires once that build starts, i.e. after the
     * original index has already been built and written. On hg38 that wastes an
     * hour and leaves a half-populated index directory. The seed's size is known
     * from the reference length alone, so refuse up front instead. */
    if (max_memory_bytes > 0) {
        int64_t l_pac = fasta_count_bases(fa);
        if (l_pac < 0) {
            fprintf(stderr, "ERROR: cannot read %s to estimate index memory\n", fa);
            return 2;
        }
        int64_t seed_need = libsais_estimate_peak_bytes(4 * l_pac);
        if (seed_need > max_memory_bytes) {
            fprintf(stderr, "[bwa_index:--meth] the seed index dominates: ~%s vs ~%s for "
                            "the original, because its per-strand-converted text is "
                            "twice as long\n",
                    bwa::fmt_bytes(seed_need).c_str(),
                    bwa::fmt_bytes(libsais_estimate_peak_bytes(2 * l_pac)).c_str());
            libsais_report_budget_shortfall("index --meth", seed_need, max_memory_bytes,
                                            max_memory_user_specified != 0, 2 * l_pac);
            return 3;
        }
    }

    /* 1. Original-alphabet index (extension reference). emit_unpacked_ref applies
     * to the original only (its `.0123` is the legacy bwa-mem2 extension target);
     * the seed index never needs an unpacked ref (step 2b passes false). */
    fprintf(stderr, "[bwa_index:--meth] building original index for %s ...\n", fa);
    if (bwa_idx_build(fa, fa, emit_unpacked_ref) != 0) {
        fprintf(stderr, "ERROR: bwa_idx_build failed on original %s\n", fa);
        return 5;
    }

    /* 2a. Emit the per-strand-converted seed FASTA <fa>.meth.fa. */
    char conv_fa[PATH_MAX];
    int n = snprintf(conv_fa, sizeof(conv_fa), "%s.meth.fa", fa);
    if (n <= 0 || (size_t)n >= sizeof(conv_fa)) {
        fprintf(stderr, "ERROR: reference path too long\n");
        return 1;
    }

    if (meth_c2t_is_fresh(fa, conv_fa)) {
        fprintf(stderr, "[bwa_index:--meth] %s is newer than %s; skipping seed FASTA emission\n",
                conv_fa, fa);
    } else {
        gzFile in = gzopen(fa, "r");
        if (in == NULL) {
            fprintf(stderr, "ERROR: cannot open %s\n", fa);
            return 2;
        }
        FILE *out = fopen(conv_fa, "w");
        if (out == NULL) {
            fprintf(stderr, "ERROR: cannot open %s for writing\n", conv_fa);
            gzclose(in);
            return 3;
        }
        fprintf(stderr, "[bwa_index:--meth] writing seed FASTA %s ...\n", conv_fa);

        kseq_t *seq = kseq_init(in);
        int64_t total_bases = 0, n_seqs = 0;
        int kr = 0;
        while ((kr = kseq_read(seq)) >= 0) {
            /* upper-case before projection so soft-masked ref regions round-trip. */
            for (size_t i = 0; i < seq->seq.l; ++i) {
                char c = seq->seq.s[i];
                if (c >= 'a' && c <= 'z') seq->seq.s[i] = (char)(c - 'a' + 'A');
            }
            meth_project_and_write(out, "r", seq->name.s, seq->seq.s, seq->seq.l, 'G', 'A');
            meth_project_and_write(out, "f", seq->name.s, seq->seq.s, seq->seq.l, 'C', 'T');
            total_bases += (int64_t)seq->seq.l;
            ++n_seqs;
        }
        kseq_destroy(seq);
        gzclose(in);
        /* kseq_read returns -1 on clean EOF; < -1 is a parse/IO error. Don't leave a
         * partial seed FASTA on disk and don't feed it to bwa_idx_build. */
        if (kr < -1) {
            fclose(out);
            unlink(conv_fa);
            fprintf(stderr, "ERROR: failed while reading %s (kseq_read=%d)\n", fa, kr);
            return 4;
        }
        if (fclose(out) != 0) {
            unlink(conv_fa);
            fprintf(stderr, "ERROR: failed to close %s\n", conv_fa);
            return 4;
        }
        fprintf(stderr, "[bwa_index:--meth] emitted %lld seqs, %lld bp (doubled to %lld bp of seed text)\n",
                (long long)n_seqs, (long long)total_bases, (long long)(2 * total_bases));
    }

    /* 2b. Build the converted seed index under the `.meth` prefix. */
    char meth_prefix[PATH_MAX];
    n = snprintf(meth_prefix, sizeof(meth_prefix), "%s.meth", fa);
    if (n <= 0 || (size_t)n >= sizeof(meth_prefix)) {
        fprintf(stderr, "ERROR: reference path too long\n");
        return 1;
    }
    fprintf(stderr, "[bwa_index:--meth] building seed index %s.* ...\n", meth_prefix);
    /* emit_unpacked_ref=false (also the default now): the seed `.0123` is never
     * read by `mem --meth` (extension uses the original reference), so don't
     * write it (~13 GB on hg38). Kept explicit for documentation; the seed
     * `.pac` + `.bwt.2bit.64` + `.ann`/`.amb` are still built. */
    if (bwa_idx_build(conv_fa, meth_prefix, /*emit_unpacked_ref=*/false) != 0) {
        fprintf(stderr, "ERROR: bwa_idx_build failed on seed index %s\n", conv_fa);
        return 5;
    }
    return 0;
}

// Parse a memory spec such as "64G", "512M", "1024K", or a bare integer
// (bytes).  Returns the number of bytes, or -1 on parse error / overflow.
static int64_t parse_memory_spec(const char *s)
{
    char *end;
    errno = 0;
    int64_t v = strtoll(s, &end, 10);
    if (errno != 0 || end == s || v < 0) return -1;
    int shift = 0;
    if      (*end == 'G' || *end == 'g') { shift = 30; ++end; }
    else if (*end == 'M' || *end == 'm') { shift = 20; ++end; }
    else if (*end == 'K' || *end == 'k') { shift = 10; ++end; }
    if (*end != '\0') return -1;
    if (shift && v > (INT64_MAX >> shift)) return -1;
    v <<= shift;
    return v;
}

static void index_usage(void)
{
	fprintf(stderr, "Usage: bwa-mem3 index [-p prefix] [-t N] [--max-memory SIZE] [--tmp-dir PATH] [--meth] <in.fasta>\n");
	fprintf(stderr, "\n"
	        "  -p STR             output prefix (default: <in.fasta>)\n"
	        "  -t INT             worker threads [auto: detected cores, cgroup-aware]\n"
	        "  --max-memory SIZE  peak memory budget; SIZE accepts a G/M/K suffix\n"
	        "                     (case-insensitive) or bare bytes\n"
	        "                     [auto: RAM less a reserve of min(max(2G, 5%%), 50%%),\n"
	        "                     cgroup-aware]\n"
	        "  --tmp-dir PATH     scratch directory [$TMPDIR]\n"
	        "  --meth             build a BS-aware dual index. Writes the original-alphabet\n"
	        "                     index at <in.fasta>.* plus a converted seed FM-index at\n"
	        "                     <in.fasta>.meth.* (used by `bwa-mem3 mem --meth`).\n"
	        "  --emit-unpacked-ref also write the unpacked `<prefix>.0123` reference. Off by\n"
	        "                     default: `mem` pac-fetches bases from `.pac`, so `.0123`\n"
	        "                     is never read. Enable only for an external consumer that\n"
	        "                     still requires it (e.g. bwa-mem2); ~8x the size of `.pac`.\n"
	        "  -u INT             SA sample rate 1/(1<<INT) [3]\n"
	        "  -h, --help         print this help message and exit\n");
}

/* Report, at index time, the ALT/AH gap `mem` also reports -- see
 * bwa_warn_sidecar_missing_AH (bwa.cpp) for why a sidecar's @SQ is emitted
 * verbatim rather than enriched. Worth saying here as well: this is when the
 * reference layout is in front of the user and regenerating the sidecar still
 * costs nothing, whereas `mem` reports it only once an alignment is running.
 *
 * `mem` remains the authoritative check -- the .alt and the sidecar are both
 * optional and either may be dropped in after indexing, so a quiet `index` is
 * not a guarantee.
 *
 * The .alt existence test is what makes this cheap: with no .alt no contig can
 * be ALT, so the check is provably a no-op, and skipping it avoids reading the
 * sidecar and the .ann/.amb (~1 MB on hg38) to reach that conclusion. */
static void warn_if_sidecar_hides_alt(const char *prefix)
{
	char alt_path[PATH_MAX];
	int n = snprintf(alt_path, sizeof(alt_path), "%s.alt", prefix);
	if (n <= 0 || (size_t)n >= sizeof(alt_path)) return;
	if (access(alt_path, F_OK) != 0) return;

	char *idx_hdr = bwa_load_hdr_from_index(prefix);   /* NULL when no sidecar */
	if (idx_hdr == NULL) return;
	bntseq_t *bns = bns_restore(prefix);               /* applies .alt -> is_alt */
	bwa_warn_sidecar_missing_AH(bns, idx_hdr, prefix); /* no-ops on a NULL bns */
	bns_destroy(bns);                                  /* NULL-safe */
	free(idx_hdr);
}

int bwa_index(int argc, char *argv[]) // the "index" command
{
	int c;
	char *prefix = 0;
	int meth = 0;
	int emit_unpacked_ref = 0;     // 0 => don't write <prefix>.0123 (mem pac-fetches)
	int64_t user_max_memory = 0;   // 0 => auto default
	int     user_threads    = 0;   // 0 => auto default
	int     user_sa_compx   = 3;   // SA sample-rate shift; 3 reproduces the historical rate
	static struct option long_opts[] = {
		{"meth",              no_argument,       0, 1000},
		{"max-memory",        required_argument, 0, 1001},
		{"tmp-dir",           required_argument, 0, 1002},
		{"emit-unpacked-ref", no_argument,       0, 1003},
		{"threads",           required_argument, 0, 't'},
		{"help",              no_argument,       0, 'h'},
		{0, 0, 0, 0}
	};
	while ((c = getopt_long(argc, argv, "p:t:u:h", long_opts, NULL)) >= 0) {
		if (c == 'p') prefix = optarg;
		else if (c == 't') {
			// Mirror parse_memory_spec's strict strtol parsing: atoi
			// silently accepts numeric-prefix garbage like "4abc" and
			// has implementation-defined behaviour on overflow.
			char *end = NULL;
			errno = 0;
			long t = strtol(optarg, &end, 10);
			if (errno || end == optarg || *end != '\0' || t <= 0 || t > INT_MAX) {
				fprintf(stderr, "ERROR: invalid -t spec '%s'\n", optarg);
				return 1;
			}
			user_threads = (int)t;
		} else if (c == 'u') {
			// Strict strtol, mirroring -t: reject non-numeric or
			// numeric-prefix-garbage input rather than silently truncating.
			// [0,6] is the writer's supported range (write_fm_index_streaming
			// requires the sample period 1<<u to divide CP_BLOCK_SIZE=64).
			char *end = NULL;
			errno = 0;
			long v = strtol(optarg, &end, 10);
			if (errno || end == optarg || *end != '\0' || v < 0 || v > 6) {
				fprintf(stderr, "[index] -u must be an integer in [0,6]\n");
				return 1;
			}
			user_sa_compx = (int)v;
		} else if (c == 1000) {
			meth = 1;
		} else if (c == 1001) {
			int64_t mem = parse_memory_spec(optarg);
			if (mem <= 0) {
				fprintf(stderr, "ERROR: invalid --max-memory spec '%s'\n", optarg);
				return 1;
			}
			user_max_memory = mem;
		} else if (c == 1002) {
			setenv("BWA_INDEX_TMPDIR", optarg, 1);
		} else if (c == 1003) {
			emit_unpacked_ref = 1;
		} else if (c == 'h') {
			index_usage();
			return 0;
		} else {
			return 1;
		}
	}

	if (optind + 1 > argc) {
		index_usage();
		return 1;
	}

	// Resolve --max-memory and -t: user value wins; otherwise auto from
	// cgroup-aware host detection. Emit a one-line audit per flag.
	int64_t resolved_mem = 0;
	{
		int64_t detected_mem = bwa::detect_total_memory_bytes();
		int     detected_cpu = bwa::detect_cpu_count();

		if (user_max_memory > 0) {
			resolved_mem = user_max_memory;
			fprintf(stderr, "[bwa_index] --max-memory = %s (user-specified)\n",
			        bwa::fmt_bytes(resolved_mem).c_str());
		} else if (detected_mem > 0) {
			// Index construction is a one-shot batch job: give it the host less
			// a reserve, not a fraction. A fractional default refuses hg38 on
			// every machine below ~144 GiB even though the build fits in ~58 GiB.
			resolved_mem = bwa::resolve_batch_memory_budget(detected_mem);
			fprintf(stderr, "[bwa_index] --max-memory = %s (auto: %s detected less a reserve, cgroup-aware)\n",
			        bwa::fmt_bytes(resolved_mem).c_str(),
			        bwa::fmt_bytes(detected_mem).c_str());
		} else {
			resolved_mem = 4LL << 30;
			fprintf(stderr, "[bwa_index] --max-memory = 4.0 GiB (fallback: host detection failed; "
			                "pass --max-memory explicitly to override)\n");
		}

		int resolved_cpu;
		if (user_threads > 0) {
			resolved_cpu = user_threads;
			fprintf(stderr, "[bwa_index] -t = %d (user-specified)\n", resolved_cpu);
		} else if (detected_cpu > 0) {
			resolved_cpu = detected_cpu;
			fprintf(stderr, "[bwa_index] -t = %d (auto: detected cores, cgroup-aware)\n", resolved_cpu);
		} else {
			resolved_cpu = 1;
			fprintf(stderr, "[bwa_index] -t = 1 (fallback: CPU detection failed; "
			                "pass -t explicitly to override)\n");
		}

		char buf[32];
		snprintf(buf, sizeof(buf), "%lld", (long long)resolved_mem);
		setenv("BWA_INDEX_MAX_MEMORY", buf, 1);
		snprintf(buf, sizeof(buf), "%d", resolved_cpu);
		setenv("BWA_INDEX_THREADS", buf, 1);
		// Diagnostics only: lets the preflight suppress a "retry on a bigger
		// host" hint when the budget came from the command line.
		setenv("BWA_INDEX_MAX_MEMORY_USER", user_max_memory > 0 ? "1" : "0", 1);
	}

	if (meth) {
		if (prefix != 0) {
			fprintf(stderr, "ERROR: --meth does not accept -p (outputs <in.fasta>.* and <in.fasta>.meth.*)\n");
			return 1;
		}
		if (user_sa_compx != 3) {
			fprintf(stderr, "ERROR: --meth does not accept -u (SA sampling rate is fixed for meth indexes)\n");
			return 1;
		}
		int rc = meth_index_build(argv[optind], emit_unpacked_ref,
		                         resolved_mem, user_max_memory > 0);
		/* The real reference only. The `.meth` seed index has no .alt of its own
		 * and its contigs are f/r-prefixed, so it can never have ALT status to
		 * lose -- the same reason `mem` skips this check in meth mode
		 * (fastmap.cpp). Stating it here beats relying on the seed prefix
		 * happening not to resolve a sidecar. */
		if (rc == 0) warn_if_sidecar_hides_alt(argv[optind]);
		return rc;
	}
	if (prefix == 0) prefix = argv[optind];
	int rc = bwa_idx_build(argv[optind], prefix, emit_unpacked_ref, user_sa_compx);
	if (rc == 0) warn_if_sidecar_hides_alt(prefix);
	return rc;
}

static void resa_usage(void)
{
	fprintf(stderr, "Usage: bwa-mem3 re-sa [-u INT] [-t INT] <idxbase>\n");
	fprintf(stderr, "\n"
	        "With -u, resample the on-disk SA sample table of an existing index to a\n"
	        "new sample rate, rewriting <idxbase>.bwt.2bit.64 in place. This persists\n"
	        "on disk what `shm -u` synthesises per-stage: a DENSER table (smaller -u)\n"
	        "means fewer LF-walk steps per SA resolve at `mem` time, trading disk +\n"
	        "memory for speed; a COARSER table (larger -u) shrinks the index. Only the\n"
	        "SA-sample layer is touched -- the BWT/cp_occ and the sibling index files\n"
	        "(.pac/.ann/.amb/.0123) are unchanged. Output is byte-identical to an index\n"
	        "built at the target rate with `index -u`.\n\n"
	        "Without -u, report the index's current on-disk SA sample rate and exit.\n\n"
	        "  -u INT             target SA sample rate 1/(1<<INT), in [0,6]. Omit to\n"
	        "                     just report the current rate. e.g. `-u 2` rewrites a\n"
	        "                     stride-8 (default, -u 3) index as a stride-4 table.\n"
	        "                     A target equal to the current rate is a no-op.\n"
	        "  -t, --threads INT  worker threads for the densify pass [1]. The\n"
	        "                     per-sample LF-walks are independent, so this scales\n"
	        "                     near-linearly; has no effect when coarsening or\n"
	        "                     inspecting.\n"
	        "  -h, --help         print this help message and exit\n\n"
	        "The write is atomic (temp file + rename), so a failed run never leaves a\n"
	        "partial index. re-sa warns if the index is currently staged in shared\n"
	        "memory (there is no staleness check): run `bwa-mem3 shm -d` before\n"
	        "re-aligning, or the stale segment will keep serving the old SA rate.\n");
}

// Rewrite an existing index's SA sample table at a new rate. See resa_usage()
// for the user-facing contract. The source .bwt.2bit.64 is mmap'd read-only, so
// the (large) cp_occ block and disk SA samples are borrowed from the page cache
// rather than copied onto the heap -- only the resampled destination arrays are
// allocated. A denser target is filled by FMI_search::densify_sa_into (the same
// LF-walk the resolver uses, so every synthesised sample is byte-identical to
// `index -u`); a coarser target is a straight decimation of the disk table
// (every (target-disk)-th sample), which needs no cp_occ walk.
int main_resa(int argc, char *argv[]) // the "re-sa" command
{
	int c;
	int target_sa_compx = -1;   // < 0 => inspect (report current rate) only
	int n_threads       = 1;

	// Honor `--help`/`-h` before getopt (which would flag -h as unknown). Skip the
	// value token of an argument-taking option so `-u --help` / `-t --help` is a
	// bad-value error from getopt below, not a spurious help request.
	for (int i = 1; i < argc; ++i) {
		const char *a = argv[i];
		if (strcmp(a, "-u") == 0 || strcmp(a, "-t") == 0 ||
		    strcmp(a, "--threads") == 0) {
			++i;   // consume its value
			continue;
		}
		if (strcmp(a, "--help") == 0 || strcmp(a, "-h") == 0) {
			resa_usage();
			return 0;
		}
	}

	optind = 1;
	opterr = 1;
#if defined(__APPLE__) || defined(__FreeBSD__) || defined(__OpenBSD__) || defined(__NetBSD__)
	extern int optreset;
	optreset = 1;
#endif
	static struct option long_opts[] = {
		{"threads", required_argument, 0, 't'},
		{"help",    no_argument,       0, 'h'},
		{0, 0, 0, 0}
	};
	while ((c = getopt_long(argc, argv, "u:t:h", long_opts, NULL)) >= 0) {
		if (c == 'u') {
			// Strict strtol, mirroring `index -u`: [0,6] is the writer-supported
			// range (1<<u must divide CP_BLOCK_SIZE=64).
			char *end = NULL;
			errno = 0;
			long v = strtol(optarg, &end, 10);
			if (errno || end == optarg || *end != '\0' || v < 0 || v > 6) {
				fprintf(stderr, "[re-sa] -u must be an integer in [0,6]\n");
				return 1;
			}
			target_sa_compx = (int)v;
		} else if (c == 't') {
			char *end = NULL;
			errno = 0;
			long v = strtol(optarg, &end, 10);
			if (errno || end == optarg || *end != '\0' || v < 1 || v > 1024) {
				fprintf(stderr, "[re-sa] -t must be an integer in [1,1024]\n");
				return 1;
			}
			n_threads = (int)v;
		} else if (c == 'h') {
			resa_usage();
			return 0;
		} else {
			return 1;   // getopt printed the error
		}
	}

	if (optind + 1 != argc) {
		fprintf(stderr, "[re-sa] expected exactly one <idxbase>\n");
		resa_usage();
		return 1;
	}
	const char *prefix = argv[optind];

	char cp_path[PATH_MAX];
	{
		int n = snprintf(cp_path, sizeof(cp_path), "%s%s", prefix, CP_FILENAME_SUFFIX);
		if (n < 0 || (size_t)n >= sizeof(cp_path)) {
			fprintf(stderr, "[re-sa] index path too long for prefix '%s'\n", prefix);
			return 1;
		}
	}

	int fd = open(cp_path, O_RDONLY);
	if (fd < 0) {
		fprintf(stderr, "[re-sa] cannot open '%s': %s\n", cp_path, strerror(errno));
		return 1;
	}
	struct stat st;
	if (fstat(fd, &st) != 0) {
		fprintf(stderr, "[re-sa] fstat('%s') failed: %s\n", cp_path, strerror(errno));
		close(fd);
		return 1;
	}
	const int64_t file_size = (int64_t)st.st_size;

	// Header: ref_seq_len, then the RAW cumulative C[] (written back verbatim).
	int64_t hdr[6];
	if (file_size < FMI_BWT2BIT_HDR_BYTES ||
	    pread(fd, hdr, sizeof(hdr), 0) != (ssize_t)sizeof(hdr)) {
		fprintf(stderr, "[re-sa] %s: truncated header\n", cp_path);
		close(fd);
		return 1;
	}
	const int64_t ref_seq_len = hdr[0];
	int64_t count_raw[5];
	memcpy(count_raw, &hdr[1], sizeof(count_raw));
	if (ref_seq_len <= 0 || ref_seq_len > 0x7fffffffffLL) {
		fprintf(stderr, "[re-sa] %s: reference_seq_len=%lld out of bounds\n",
		        cp_path, (long long)ref_seq_len);
		close(fd);
		return 1;
	}
	for (int i = 0; i < 5; ++i) {
		if (count_raw[i] < 0 || count_raw[i] > ref_seq_len ||
		    (i > 0 && count_raw[i] < count_raw[i - 1])) {
			fprintf(stderr, "[re-sa] %s: count[] corrupt (bad/non-monotonic)\n", cp_path);
			close(fd);
			return 1;
		}
	}

	// Disk SA rate: the trailing tag if present, else the legacy default.
	const int64_t disk_sa_compx =
		detect_sa_compx(fd, file_size, ref_seq_len, SA_COMPX);

	// -t only affects densification (a target below the current rate). Warn if it
	// was set but no densify will happen (inspection, no-op, or coarsen), mirroring
	// `shm`'s "--threads without -u" warning.
	if (n_threads > 1 && !(target_sa_compx >= 0 && target_sa_compx < disk_sa_compx)) {
		fprintf(stderr,
		        "[re-sa] -t/--threads only affects densification (a -u below the "
		        "current rate 1/%d); ignoring\n", 1 << (int)disk_sa_compx);
	}

	// Inspection mode: no -u given -> report the current rate and exit.
	if (target_sa_compx < 0) {
		int64_t disk_count = (ref_seq_len >> disk_sa_compx) + 1;
		fprintf(stderr,
		        "[re-sa] %s: SA sample rate 1/%d (shift %lld), ref_seq_len %lld, "
		        "%lld samples\n",
		        cp_path, 1 << (int)disk_sa_compx, (long long)disk_sa_compx,
		        (long long)ref_seq_len, (long long)disk_count);
		close(fd);
		return 0;
	}

	if (target_sa_compx == disk_sa_compx) {
		fprintf(stderr,
		        "[re-sa] %s already at SA rate 1/%d (shift %d); nothing to do\n",
		        cp_path, 1 << (int)disk_sa_compx, (int)disk_sa_compx);
		close(fd);
		return 0;
	}

	// From here we will mutate the on-disk index. Shared memory has no staleness
	// check, so a live segment would keep serving the OLD SA rate. Alignments
	// stay correct (the rate change is alignment-invariant); the cost is that the
	// new rate's performance is silently not picked up until the segment is
	// re-staged. Warn (best-effort) but proceed -- consistent with how the
	// `index`/`shm` footgun is handled.
	if (bwa_shm_test(prefix) == 1) {
		fprintf(stderr,
		        "[re-sa] WARNING: '%s' is currently staged in shared memory. This "
		        "rewrite will NOT update the staged segment; run `bwa-mem3 shm -d` "
		        "(and re-stage if desired) or `mem` will keep serving the old SA "
		        "rate.\n", prefix);
	}

	// Section offsets in the DISK layout (disk stride), from the shared layout
	// helper so the reader agrees bit-for-bit with the writer. The target-rate
	// layout gives the destination sample count. sentinel_index sits right after
	// the SA arrays.
	const Bwt2bitLayout disk_layout = fmi_bwt2bit_layout(ref_seq_len, (int)disk_sa_compx);
	const int64_t cp_occ_size = disk_layout.cp_occ_count;
	const int64_t disk_count  = disk_layout.sa_sample_count;
	const int64_t dst_count   = fmi_bwt2bit_layout(ref_seq_len, target_sa_compx).sa_sample_count;
	const int64_t off_ms   = disk_layout.off_ms_byte;
	const int64_t off_ls   = disk_layout.off_ls_word;
	const int64_t off_sent = disk_layout.off_sentinel;
	if (file_size < off_sent + (int64_t)sizeof(int64_t)) {
		fprintf(stderr, "[re-sa] %s: file too small for its declared layout "
		        "(corrupt?)\n", cp_path);
		close(fd);
		return 1;
	}

	// Borrow cp_occ from the page cache instead of copying it onto the heap: it
	// is read-only here and dominates the index size (~6.4 GB on hg38). cp_occ
	// sits at FMI_BWT2BIT_HDR_BYTES (mod 8 == 0) on a page-aligned map, so its
	// 8-byte-aligned int64 fields are safe to read in place.
	void *map = mmap(NULL, (size_t)file_size, PROT_READ, MAP_PRIVATE, fd, 0);
	if (map == MAP_FAILED) {
		fprintf(stderr, "[re-sa] mmap('%s', %lld bytes) failed: %s\n",
		        cp_path, (long long)file_size, strerror(errno));
		close(fd);
		return 1;
	}
	const uint8_t  *base    = (const uint8_t *)map;
	const CP_OCC   *cp_occ  = (const CP_OCC   *)(base + FMI_BWT2BIT_HDR_BYTES);
	const int8_t   *disk_ms = (const int8_t   *)(base + off_ms);   // int8: any alignment
	int64_t sentinel_index;
	memcpy(&sentinel_index, base + off_sent, sizeof(int64_t));
	if (sentinel_index < 0 || sentinel_index >= ref_seq_len) {
		fprintf(stderr, "[re-sa] %s: sentinel_index=%lld out of bounds\n",
		        cp_path, (long long)sentinel_index);
		munmap(map, (size_t)file_size);
		close(fd);
		return 1;
	}

	// sa_ls_word is uint32 but its file offset follows the odd-length sa_ms
	// array, so base+off_ls is not guaranteed 4-byte aligned — reading uint32
	// through it in place is UB (every other reader freads it into an aligned
	// buffer). Copy it into an aligned heap buffer; cp_occ stays borrowed.
	uint32_t *disk_ls = (uint32_t *)malloc((size_t)disk_count * sizeof(uint32_t));
	int8_t   *dst_ms  = (int8_t   *)malloc((size_t)dst_count  * sizeof(int8_t));
	uint32_t *dst_ls  = (uint32_t *)malloc((size_t)dst_count  * sizeof(uint32_t));
	if (!disk_ls || !dst_ms || !dst_ls) {
		fprintf(stderr, "[re-sa] out of memory sizing SA resample buffers\n");
		free(disk_ls); free(dst_ms); free(dst_ls);
		munmap(map, (size_t)file_size);
		close(fd);
		return 1;
	}
	memcpy(disk_ls, base + off_ls, (size_t)disk_count * sizeof(uint32_t));

	const bool densify = (target_sa_compx < disk_sa_compx);
	struct timespec t0, t1;
	clock_gettime(CLOCK_MONOTONIC, &t0);
	if (densify) {
		// Recover each added sample by the resolver's own LF-walk against the disk
		// samples. densify_sa_into needs the +1-adjusted C[].
		fprintf(stderr, "[re-sa] densifying SA samples %lld -> %lld "
		        "(rate 1/%d -> 1/%d, %d thread%s); this can take a while for large "
		        "indices...\n",
		        (long long)disk_count, (long long)dst_count,
		        1 << (int)disk_sa_compx, 1 << target_sa_compx,
		        n_threads, n_threads == 1 ? "" : "s");
		int64_t count_adj[5];
		for (int i = 0; i < 5; ++i) count_adj[i] = count_raw[i] + 1;
		FMI_search fmi(prefix);
		fmi.densify_sa_into(cp_occ, count_adj, ref_seq_len, sentinel_index,
		                    disk_ms, disk_ls, disk_sa_compx,
		                    dst_ms, dst_ls, target_sa_compx, n_threads);
	} else {
		// Coarsen: every target-sampled row is also a disk-sampled row (both
		// periods divide 64 and target > disk), so decimate the disk table --
		// row = j << target, disk index = row >> disk = j << (target - disk).
		// Serial: this is a memory-bandwidth-bound gather that gains nothing from
		// threads, so -t is (as documented) ignored here.
		const int shift = target_sa_compx - (int)disk_sa_compx;
		for (int64_t j = 0; j < dst_count; ++j) {
			int64_t s = j << shift;
			dst_ms[j] = disk_ms[s];
			dst_ls[j] = disk_ls[s];
		}
	}
	clock_gettime(CLOCK_MONOTONIC, &t1);
	const double fill_secs = (t1.tv_sec - t0.tv_sec)
	                       + (t1.tv_nsec - t0.tv_nsec) / 1e9;
	fprintf(stderr, "[re-sa] %s SA samples %lld -> %lld (rate 1/%d -> 1/%d) in %.1fs\n",
	        densify ? "densified" : "coarsened",
	        (long long)disk_count, (long long)dst_count,
	        1 << (int)disk_sa_compx, 1 << target_sa_compx, fill_secs);

	// cp_occ is copied through unchanged (it does not depend on the SA rate); the
	// writer reads it from the mmap and rewrites the file atomically (tmp+rename)
	// before we unmap.
	rewrite_fm_index_resampled_sa(cp_path, ref_seq_len, count_raw,
	                              cp_occ, cp_occ_size, dst_ms, dst_ls, dst_count,
	                              sentinel_index, target_sa_compx);

	free(disk_ls); free(dst_ms); free(dst_ls);
	munmap(map, (size_t)file_size);
	close(fd);
	fprintf(stderr, "[re-sa] wrote %s at SA rate 1/%d (shift %d)\n",
	        cp_path, 1 << target_sa_compx, target_sa_compx);
	return 0;
}

int bwa_idx_build(const char *fa, const char *prefix, int emit_unpacked_ref, int sa_compx)
{
	extern void bwa_pac_rev_core(const char *fn, const char *fn_rev);

	clock_t t;
	int rc = 0;

	{ // nucleotide indexing
		gzFile fp = xzopen(fa, "r");
		t = clock();
		fprintf(stderr, "[bwa_index] Pack FASTA... ");
		bns_fasta2bntseq(fp, prefix, 1);
		fprintf(stderr, "%.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);
		err_gzclose(fp);
        FMI_search *fmi = new FMI_search(prefix);
        rc = fmi->build_index(emit_unpacked_ref, sa_compx);
        delete fmi;
	}
	return rc;
}
