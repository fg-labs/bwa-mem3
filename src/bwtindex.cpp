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
#include <sys/stat.h>
#include <zlib.h>

#include "bntseq.h"
#include "bwa.h"
#include "bwt.h"
#include "utils.h"
#include "FMI_search.h"
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
	static struct option long_opts[] = {
		{"meth",              no_argument,       0, 1000},
		{"max-memory",        required_argument, 0, 1001},
		{"tmp-dir",           required_argument, 0, 1002},
		{"emit-unpacked-ref", no_argument,       0, 1003},
		{"threads",           required_argument, 0, 't'},
		{"help",              no_argument,       0, 'h'},
		{0, 0, 0, 0}
	};
	while ((c = getopt_long(argc, argv, "p:t:h", long_opts, NULL)) >= 0) {
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
	int rc = bwa_idx_build(argv[optind], prefix, emit_unpacked_ref);
	if (rc == 0) warn_if_sidecar_hides_alt(prefix);
	return rc;
}

int bwa_idx_build(const char *fa, const char *prefix, int emit_unpacked_ref)
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
        rc = fmi->build_index(emit_unpacked_ref);
        delete fmi;
	}
	return rc;
}
