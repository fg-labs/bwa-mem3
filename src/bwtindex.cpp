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
#include <string>
#include <vector>

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

/* Scratch is registered with atexit as well as scoped cleanup: index writers
 * below this layer can call exit(). SIGKILL and machine failure cannot clean up. */
static std::string meth_scratch_dir;
static void meth_cleanup_scratch()
{
    if (!meth_scratch_dir.empty()) {
        unlink((meth_scratch_dir + "/seed.fa").c_str());
        rmdir(meth_scratch_dir.c_str());
        meth_scratch_dir.clear();
    }
}
struct MethScratchGuard {
    ~MethScratchGuard() { meth_cleanup_scratch(); }
};

/* Use the same fallback naming as the alignment header loader. */
static std::string meth_dict_path(const std::string& prefix)
{
    char *path = bwa_dict_path(prefix.c_str());
    std::string result(path);
    free(path);
    return result;
}

/* Read and validate the whole source before touching any destination. zlib can
 * report a truncated gzip stream as EOF to kseq, so check both error channels. */
static int64_t fasta_count_bases(const char *fa)
{
    gzFile in = gzopen(fa, "r");
    if (!in) return -1;
    kseq_t *seq = kseq_init(in);
    int64_t total = 0;
    int kr;
    while ((kr = kseq_read(seq)) >= 0) {
        if (seq->seq.l > (size_t)((INT64_MAX / 12 - 1) / 4 - total)) { kr = -2; break; }
        total += (int64_t)seq->seq.l;
    }
    int ze;
    gzerror(in, &ze);
    kseq_destroy(seq);
    int close_rc = gzclose(in);
    return kr < -1 || (ze != Z_OK && ze != Z_STREAM_END) || close_rc != Z_OK ? -1 : total;
}

/* Copy sidecars verbatim, preserving the header reader's .hdr-over-.dict rule. */
static bool meth_copy_file(const std::string& source, const std::string& dest)
{
    FILE *in = fopen(source.c_str(), "rb");
    if (!in) return false;
    FILE *out = fopen(dest.c_str(), "wb");
    if (!out) { fclose(in); return false; }
    char buffer[65536];
    size_t n;
    bool ok = true;
    while ((n = fread(buffer, 1, sizeof(buffer), in)) != 0)
        if (fwrite(buffer, 1, n, out) != n) { ok = false; break; }
    ok = ok && !ferror(in) && !ferror(out);
    if (fclose(in) != 0) ok = false;
    if (fclose(out) != 0) ok = false;
    if (!ok) fprintf(stderr, "ERROR: copying %s to %s failed\n", source.c_str(), dest.c_str());
    return ok;
}

/* Resolve aliases, including a symlinked parent directory, before comparing
 * names. The output parent must already exist; indexing never creates it. */
static std::string meth_absolute_path(const std::string& path)
{
    size_t slash = path.find_last_of('/');
    std::string parent = slash == std::string::npos ? "." : path.substr(0, slash + 1);
    char resolved[PATH_MAX];
    if (!realpath(parent.c_str(), resolved)) return "";
    return std::string(resolved) + "/" + path.substr(slash == std::string::npos ? 0 : slash + 1);
}

static bool meth_same_file(const std::string& a, const std::string& b)
{
    struct stat sa, sb;
    return (stat(a.c_str(), &sa) == 0 && stat(b.c_str(), &sb) == 0 &&
            sa.st_dev == sb.st_dev && sa.st_ino == sb.st_ino) ||
           meth_absolute_path(a) == meth_absolute_path(b);
}

/* Convert into private scratch, checking parser, gzip, write and close errors. */
static bool meth_write_seed_fasta(const char *fa, const std::string& conv)
{
    gzFile in = gzopen(fa, "r");
    FILE *out = fopen(conv.c_str(), "wb");
    if (!in || !out) {
        if (in) gzclose(in);
        if (out) fclose(out);
        fprintf(stderr, "ERROR: cannot open methylation conversion input/output\n");
        return false;
    }
    fprintf(stderr, "[bwa_index:--meth] writing fresh seed FASTA %s\n", conv.c_str());
    kseq_t *seq = kseq_init(in);
    int kr;
    while ((kr = kseq_read(seq)) >= 0) {
        for (size_t i = 0; i < seq->seq.l; ++i) {
            char c = seq->seq.s[i];
            if (c >= 'a' && c <= 'z') seq->seq.s[i] = (char)(c - 'a' + 'A');
        }
        meth_project_and_write(out, "r", seq->name.s, seq->seq.s, seq->seq.l, 'G', 'A');
        meth_project_and_write(out, "f", seq->name.s, seq->seq.s, seq->seq.l, 'C', 'T');
        if (ferror(out)) break;
    }
    int ze;
    gzerror(in, &ze);
    bool ok = kr == -1 && (ze == Z_OK || ze == Z_STREAM_END) && !ferror(out);
    kseq_destroy(seq);
    if (gzclose(in) != Z_OK) ok = false;
    if (fclose(out) != 0) ok = false;
    if (!ok) {
        fprintf(stderr, "ERROR: failed reading or writing methylation conversion\n");
        return false;
    }
    return true;
}

/* Publish retained conversion atomically; scratch and output may be on different
 * filesystems. A failed copy never exposes a partial retained FASTA. */
static bool meth_retain_fasta(const std::string& source, const std::string& dest)
{
    if (rename(source.c_str(), dest.c_str()) == 0) return true;
    if (errno != EXDEV) {
        fprintf(stderr, "ERROR: retaining %s: %s\n", dest.c_str(), strerror(errno));
        return false;
    }
    std::string pattern = dest + ".XXXXXX";
    std::vector<char> name(pattern.begin(), pattern.end());
    name.push_back('\0');
    int fd = mkstemp(name.data());
    if (fd < 0) { fprintf(stderr, "ERROR: retaining %s: %s\n", dest.c_str(), strerror(errno)); return false; }
    close(fd);
    bool ok = meth_copy_file(source, name.data());
    if (ok && rename(name.data(), dest.c_str()) != 0) {
        fprintf(stderr, "ERROR: publishing %s: %s\n", dest.c_str(), strerror(errno));
        ok = false;
    }
    unlink(name.data());
    return ok;
}

/* Only these known artifacts may be replaced. No globbing, recursive removal,
 * or symlink following is permitted, even when --force was requested. */
static int meth_index_build(const char *fa, const char *prefix, bool minimal,
                            bool keep_fasta, bool force, int emit_unpacked_ref,
                            int64_t max_memory_bytes, int max_memory_user_specified)
{
    std::string dest(prefix), source(fa);
    if (dest.empty() || dest.back() == '/' || dest.size() + 32 >= PATH_MAX ||
        (dest.size() >= 5 && dest.compare(dest.size() - 5, 5, ".meth") == 0)) {
        fprintf(stderr, "ERROR: invalid methylation output prefix (must not end in .meth)\n");
        return 1;
    }
    if (meth_absolute_path(dest).empty()) {
        fprintf(stderr, "ERROR: output parent directory does not exist: %s\n", prefix);
        return 1;
    }
    char *header_path = NULL;
    int header_error = 0;
    FILE *header = bwa_open_index_header(fa, 1, &header_path, &header_error);
    std::string hdr(header_path);
    free(header_path);
    if (header_error) {
        fprintf(stderr, "ERROR: cannot read metadata %s: %s\n", hdr.c_str(), strerror(header_error));
        return 1;
    }
    const bool have_hdr = header != NULL;
    if (header) fclose(header);
    std::string alt = source + ".alt";
    struct stat alt_stat;
    bool have_alt = stat(alt.c_str(), &alt_stat) == 0;
    if (!have_alt && (errno != ENOENT || lstat(alt.c_str(), &alt_stat) == 0)) {
        fprintf(stderr, "ERROR: cannot inspect metadata %s\n", alt.c_str());
        return 1;
    }
    if (have_alt) {
        FILE *fp = S_ISREG(alt_stat.st_mode) ? fopen(alt.c_str(), "r") : NULL;
        if (!fp) { fprintf(stderr, "ERROR: cannot read metadata %s\n", alt.c_str()); return 1; }
        fclose(fp);
    }
    // In-place builds must continue to consult the source dictionary on each
    // alignment, not shadow future updates with a newly materialized .hdr.
    const bool relocated = meth_absolute_path(source) != meth_absolute_path(dest);
    std::vector<std::pair<std::string, std::string>> sidecars;
    if (have_alt) sidecars.emplace_back(alt, dest + ".alt");
    if (have_hdr) sidecars.emplace_back(hdr, relocated ? dest + ".hdr" : hdr);
    std::vector<std::string> protected_paths = {source, source + ".alt", source + ".hdr", meth_dict_path(source)};
    std::vector<std::string> remove_paths;
    const char *suffixes[] = {".amb", ".ann", ".pac", ".bwt.2bit.64", ".0123",
        ".meth.amb", ".meth.ann", ".meth.pac", ".meth.bwt.2bit.64", ".meth.0123",
        ".meth.fa", ".alt", ".hdr"};
    for (const char *suffix : suffixes) {
        std::string path = dest + suffix;
        struct stat st;
        if (lstat(path.c_str(), &st) != 0) {
            if (errno == ENOENT) continue;
            fprintf(stderr, "ERROR: cannot inspect output %s: %s\n", path.c_str(), strerror(errno));
            return 1;
        }
        // Sidecars already at their source name are inputs, not replacements.
        bool input_sidecar = false;
        for (const auto& mapping : sidecars)
            if (path == mapping.second &&
                meth_absolute_path(mapping.first) == meth_absolute_path(mapping.second))
                input_sidecar = true;
        if (input_sidecar) continue;
        bool protected_file = false;
        for (const auto& input : protected_paths)
            if (meth_same_file(path, input)) protected_file = true;
        if (!S_ISREG(st.st_mode) || st.st_nlink != 1 || protected_file) {
            fprintf(stderr, "ERROR: unsafe output collision: %s (source, alias, symlink or non-regular file)\n", path.c_str());
            return 1;
        }
        if (!force) {
            fprintf(stderr, "ERROR: output exists: %s; use --force to replace known index artifacts\n", path.c_str());
            return 1;
        }
        remove_paths.push_back(path);
    }
    // Seed metadata would describe the wrong contigs. An unrelated destination
    // dictionary must not silently become the relocated original's header.
    for (const auto& path : {dest + ".meth.alt", dest + ".meth.hdr", meth_dict_path(dest), meth_dict_path(dest + ".meth")}) {
        struct stat st;
        if (lstat(path.c_str(), &st) == 0 &&
            !(path == meth_dict_path(dest) && meth_same_file(path, meth_dict_path(source)))) {
            fprintf(stderr, "ERROR: conflicting metadata: %s; refusing to replace unrelated files\n", path.c_str());
            return 1;
        }
    }

    int64_t length = fasta_count_bases(fa);
    if (length <= 0) {
        fprintf(stderr, "ERROR: cannot read a nonempty reference from %s\n", fa);
        return 2;
    }
    int64_t seed_need = libsais_estimate_peak_bytes(4 * length);
    if (max_memory_bytes > 0 && seed_need > max_memory_bytes) {
        fprintf(stderr, "[bwa_index:--meth] the seed index dominates: ~%s vs ~%s for "
                        "the original, because its per-strand-converted text is twice as long\n",
                bwa::fmt_bytes(seed_need).c_str(),
                bwa::fmt_bytes(libsais_estimate_peak_bytes(2 * length)).c_str());
        libsais_report_budget_shortfall("index --meth", seed_need, max_memory_bytes,
                                        max_memory_user_specified != 0, 2 * length);
        return 3;
    }

    const char *tmp = getenv("BWA_INDEX_TMPDIR");
    if (!tmp) tmp = getenv("TMPDIR");
    if (!tmp) tmp = "/tmp";
    std::string pattern = std::string(tmp) + "/bwa-meth-XXXXXX";
    std::vector<char> name(pattern.begin(), pattern.end());
    name.push_back('\0');
    if (!mkdtemp(name.data())) {
        fprintf(stderr, "ERROR: cannot create methylation scratch in %s: %s\n", tmp, strerror(errno));
        return 1;
    }
    meth_scratch_dir = name.data();
    static bool registered = false;
    if (!registered) { atexit(meth_cleanup_scratch); registered = true; }
    MethScratchGuard guard;
    std::string conv = meth_scratch_dir + "/seed.fa";
    if (!meth_write_seed_fasta(fa, conv)) return 2;
    for (const auto& path : remove_paths) {
        if (unlink(path.c_str()) != 0) {
            fprintf(stderr, "ERROR: cannot remove %s: %s\n", path.c_str(), strerror(errno));
            return 1;
        }
    }
    if (minimal) {
        gzFile fp = xzopen(fa, "r");
        bns_fasta2bntseq(fp, prefix, 1);
        err_gzclose(fp);
    } else if (bwa_idx_build(fa, prefix, emit_unpacked_ref) != 0) return 5;
    if (bwa_idx_build(conv.c_str(), (dest + ".meth").c_str(), false) != 0) return 5;
    for (const auto& mapping : sidecars)
        if (!meth_same_file(mapping.first, mapping.second) &&
            !meth_copy_file(mapping.first, mapping.second)) return 5;
    if (keep_fasta && !meth_retain_fasta(conv, dest + ".meth.fa")) return 5;
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
	        "  --meth             build original <prefix>.* and converted <prefix>.meth.*\n"
	        "                     indexes; -p defaults to <in.fasta>. Use with mem --meth.\n"
	        "  --meth-only        with --meth, pack the original without its FM-index;\n"
	        "                     seven index files, for methylation alignment only\n"
	        "  --keep-meth-fasta  with --meth, retain fresh conversion as <prefix>.meth.fa\n"
	        "                     (otherwise removed from unique scratch on exit)\n"
	        "  --force            with --meth, replace existing known index artifacts;\n"
	        "                     rejects source aliases, symlinks and directories\n"
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
	bool meth_only = false, keep_meth_fasta = false, force = false;
	int emit_unpacked_ref = 0;     // 0 => don't write <prefix>.0123 (mem pac-fetches)
	int64_t user_max_memory = 0;   // 0 => auto default
	int     user_threads    = 0;   // 0 => auto default
	int     user_sa_compx   = 3;   // SA sample-rate shift; 3 reproduces the historical rate
	static struct option long_opts[] = {
		{"meth",              no_argument,       0, 1000},
		{"max-memory",        required_argument, 0, 1001},
		{"tmp-dir",           required_argument, 0, 1002},
		{"emit-unpacked-ref", no_argument,       0, 1003},
		{"meth-only",         no_argument,       0, 1004},
		{"keep-meth-fasta",   no_argument,       0, 1005},
		{"force",             no_argument,       0, 1006},
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
		} else if (c == 1004) {
			meth_only = true;
		} else if (c == 1005) {
			keep_meth_fasta = true;
		} else if (c == 1006) {
			force = true;
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

	if (!meth && (meth_only || keep_meth_fasta || force)) {
		fprintf(stderr, "ERROR: --meth-only, --keep-meth-fasta and --force require --meth\n");
		return 1;
	}
	if (meth_only && emit_unpacked_ref) {
		fprintf(stderr, "ERROR: --meth-only cannot be combined with --emit-unpacked-ref\n");
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
		if (prefix == 0) prefix = argv[optind];
		if (user_sa_compx != 3) {
			fprintf(stderr, "ERROR: --meth does not accept -u (SA sampling rate is fixed for meth indexes)\n");
			return 1;
		}
		int rc = meth_index_build(argv[optind], prefix, meth_only, keep_meth_fasta, force, emit_unpacked_ref,
		                         resolved_mem, user_max_memory > 0);
		/* The real reference only. The `.meth` seed index has no .alt of its own
		 * and its contigs are f/r-prefixed, so it can never have ALT status to
		 * lose -- the same reason `mem` skips this check in meth mode
		 * (fastmap.cpp). Stating it here beats relying on the seed prefix
		 * happening not to resolve a sidecar. */
		if (rc == 0) warn_if_sidecar_hides_alt(prefix);
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
