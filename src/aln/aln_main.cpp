/* New code for bwa-mem3 (not from upstream bwa): the `aln` subcommand driver. The option
 * set is stock `bwa aln`'s plus samse's -n (as --max-xa-hits), -f and -r; see README.md. */
#include "aln_main.h"

#include <getopt.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <errno.h>
#include <sys/stat.h>

#include <algorithm>
#include <vector>

#include "../FMI_search.h"
#include "../bwa.h"
#include "../kthread.h"
#include "../read_index_ele.h"
#include "aln_bwt.h"
#include "bwape.h"
#include "bwase.h"
#include "bwaseqio.h"
#include "bwtgap.h"

extern char *bwa_pg;

namespace bwa_aln {

void aln_opt_init(aln_opt_t *o)
{
	gap_opt_t *g = gap_init_opt();
	o->gap = *g;
	free(g);
	o->max_xa_hits = 3;
	pe_opt_t *pe = bwa_init_pe_opt();
	o->pe = *pe;
	free(pe);
	o->batch_size = 0x40000;
}

/* The stages after the search run in parallel on kt_for (design §6.2): every drand48() draw
 * stays in one serial pass, in stock order, before them; infer_isize is a barrier between the
 * paired-end passes; the records are formatted in parallel and written in read order. */

namespace {

/* One formatted run of consecutive reads (or pairs) [st, st + n): `len` bytes at `off` in
 * thread `tid`'s buffer. `stop` >= 0: the pair at index `stop` has mismatched names, and the
 * run ends with that pair's records. */
struct fmt_seg_t {
	int st, tid, stop;
	size_t off, len;
};

/* Per-thread output buffers and segment lists for one pass, reused across batches. */
struct fmt_bufs_t {
	std::vector<kstring_t> buf;
	std::vector<std::vector<fmt_seg_t>> segs;
	explicit fmt_bufs_t(int n_threads) : buf(n_threads), segs(n_threads)
	{
		for (auto &k : buf) k.l = k.m = 0, k.s = 0;
	}
	~fmt_bufs_t()
	{
		for (auto &k : buf) free(k.s);
	}
	void reset()
	{
		for (auto &k : buf) k.l = 0;
		for (auto &v : segs) v.clear();
	}
	/* Writes the segments in read order, up to and including the first stopped one. Sets
	 * *stop to that segment's stopping index (or -1). Returns 0 or -1 on a write error. */
	int write(aln_out_t *out, int *stop)
	{
		std::vector<fmt_seg_t> all;
		for (auto &v : segs) all.insert(all.end(), v.begin(), v.end());
		std::sort(all.begin(), all.end(), [](const fmt_seg_t &x, const fmt_seg_t &y) { return x.st < y.st; });
		*stop = -1;
		for (const fmt_seg_t &g : all) {
			if (aln_out_write_buf(out, buf[g.tid].s + g.off, g.len) != 0) return -1;
			if (g.stop >= 0) {
				*stop = g.stop;
				break;
			}
		}
		return 0;
	}
};

struct se_ctx_t {
	const aln_opt_t *o;
	const bwt_t *bwt;
	const bntseq_t *bns;
	ubyte_t *pac;
	bwa_seq_t *seqs;
	fmt_bufs_t *fb;
};

/* Single end, reads [st, st + n): SA to position, refinement, MD, trim correction and SAM text;
 * bwa_cal_pac_pos and bwa_refine_gapped are upstream's batch loops run on the slice. */
void se_post_worker(void *data, int st, int n, int tid)
{
	se_ctx_t *c = (se_ctx_t*)data;
	bwa_seq_t *seqs = c->seqs + st;
	bwa_cal_pac_pos(c->bns, c->bwt, n, seqs, c->o->gap.max_diff, c->o->gap.fnr);
	bwa_refine_gapped(c->bns, n, seqs, c->pac);
	kstring_t *k = &c->fb->buf[tid];
	fmt_seg_t g;
	g.st = st, g.tid = tid, g.stop = -1, g.off = k->l;
	for (int i = 0; i < n; ++i)
		bwa_print_sam1(k, c->bns, seqs + i, 0, c->o->gap.mode, c->o->gap.max_top2);
	g.len = k->l - g.off;
	c->fb->segs[tid].push_back(g);
}

struct pe_ctx_t {
	const aln_opt_t *o;
	const bwt_t *bwt;
	const bntseq_t *bns;
	ubyte_t *pac;
	bwa_seq_t *seqs[2];
	const isize_info_t *ii;
	fmt_bufs_t *fb;
	std::vector<long long> cnt_chg;       // per thread
	std::vector<uint64_t> sw_cnt;         // per thread: n_tot[0..1], n_mapped[0..1]
};

/* Paired end, pass 1 (pairs [st, st + n)): single-end mapQ and position. */
void pe_pos_worker(void *data, int st, int n, int tid)
{
	(void)tid;
	pe_ctx_t *c = (pe_ctx_t*)data;
	bwa_seq_t *s[2] = {c->seqs[0] + st, c->seqs[1] + st};
	bwa_pe_se_pos(c->bns, c->bwt, n, s, &c->o->gap);
}

/* Paired end, pass 2 (after infer_isize): pairing and extra hits, mate Smith-Waterman,
 * refinement of both reads, then each pair's records, read 1 then read 2. */
void pe_post_worker(void *data, int st, int n, int tid)
{
	pe_ctx_t *c = (pe_ctx_t*)data;
	bwa_seq_t *s[2] = {c->seqs[0] + st, c->seqs[1] + st};
	c->cnt_chg[tid] += bwa_pe_pair(c->bns, c->bwt, n, s, c->ii, &c->o->pe, &c->o->gap);
	uint64_t *sw = &c->sw_cnt[4 * tid];
	bwa_paired_sw_core(c->bns, c->pac, n, s, &c->o->pe, c->ii, sw, sw + 2);
	for (int j = 0; j < 2; ++j) bwa_refine_gapped(c->bns, n, s[j], c->pac);
	kstring_t *k = &c->fb->buf[tid];
	fmt_seg_t g;
	g.st = st, g.tid = tid, g.stop = -1, g.off = k->l;
	for (int i = 0; i < n; ++i) {
		bwa_seq_t *p0 = s[0] + i, *p1 = s[1] + i;
		if (p0->bc[0] || p1->bc[0]) {
			strcat(p0->bc, p1->bc);
			strcpy(p1->bc, p0->bc);
		}
		bwa_print_sam1(k, c->bns, p0, p1, c->o->gap.mode, c->o->gap.max_top2);
		bwa_print_sam1(k, c->bns, p1, p0, c->o->gap.mode, c->o->gap.max_top2);
		if (strcmp(p0->name, p1->name) != 0) { // as sampe: stop after this pair
			g.stop = st + i;
			break;
		}
	}
	g.len = k->l - g.off;
	c->fb->segs[tid].push_back(g);
}

} // namespace

int aln_run_se(const aln_opt_t *o, const bwt_t *bwt, const bntseq_t *bns, ubyte_t *pac,
			   const char *fn_fq, aln_out_t *out)
{
	bwa_seqio_t *ks = bwa_seq_open(fn_fq);
	if (!ks) {
		fprintf(stderr, "[E::%s] cannot open '%s'\n", __func__, fn_fq);
		return 1;
	}
	bwase_initialize();
	srand48(bns->seed); // once per run, never per batch
	const int n_threads = o->gap.n_threads < 1 ? 1 : o->gap.n_threads;
	fmt_bufs_t fb(n_threads);
	long long tot = 0;
	int n_seqs, ret = 0;
	bwa_seq_t *seqs;
	while ((seqs = bwa_read_seq(ks, o->batch_size, &n_seqs, o->gap.mode, o->gap.trim_qual)) != 0) {
		tot += n_seqs;
		aln_search_batch(bwt, n_seqs, seqs, &o->gap, o->gap.n_threads);
		// serial and in read order: bwa_aln2seq_core draws from drand48
		for (int i = 0; i < n_seqs; ++i)
			bwa_aln2seq_core(seqs[i].n_aln, seqs[i].aln, &seqs[i], 1, o->max_xa_hits);
		se_ctx_t c = {o, bwt, bns, pac, seqs, &fb};
		fb.reset();
		kt_for_nt(se_post_worker, &c, n_seqs, n_threads);
		int stop;
		if (fb.write(out, &stop) != 0) ret = 1;
		bwa_free_read_seq(n_seqs, seqs);
		if (ret) break;
		if (bwa_verbose >= 3) fprintf(stderr, "[M::aln] %lld sequences have been processed.\n", tot);
	}
	bwa_seq_close(ks);
	return ret;
}

int aln_run_pe(const aln_opt_t *o, const bwt_t *bwt, const bntseq_t *bns, ubyte_t *pac,
			   const char *fn_fq1, const char *fn_fq2, aln_out_t *out)
{
	bwa_seqio_t *ks[2];
	ks[0] = bwa_seq_open(fn_fq1);
	if (!ks[0]) {
		fprintf(stderr, "[E::%s] cannot open '%s'\n", __func__, fn_fq1);
		return 1;
	}
	ks[1] = bwa_seq_open(fn_fq2);
	if (!ks[1]) {
		fprintf(stderr, "[E::%s] cannot open '%s'\n", __func__, fn_fq2);
		bwa_seq_close(ks[0]);
		return 1;
	}
	bwase_initialize();
	srand48(bns->seed); // once per run, never per batch
	const int n_threads = o->gap.n_threads < 1 ? 1 : o->gap.n_threads;
	fmt_bufs_t fb(n_threads);
	isize_info_t last_ii;
	last_ii.avg = -1.0;
	long long tot = 0;
	int ret = 0;
	bwa_seq_t *seqs[2];
	int n_in[2];
	while ((seqs[0] = bwa_read_seq(ks[0], o->batch_size, &n_in[0], o->gap.mode, o->gap.trim_qual)) != 0) {
		seqs[1] = bwa_read_seq(ks[1], o->batch_size, &n_in[1], o->gap.mode, o->gap.trim_qual);
		if (!seqs[1]) {
			fprintf(stderr, "[W::aln] read 2 file ended before read 1 file; ignoring the remaining %d read(s) of read 1\n", n_in[0]);
			bwa_free_read_seq(n_in[0], seqs[0]);
			break;
		}
		int n_seqs = n_in[1]; // as sampe: the second file's count
		if (n_in[0] != n_in[1]) {
			if (n_seqs > n_in[0]) n_seqs = n_in[0]; // sampe would read past the first file's reads
			fprintf(stderr, "[W::aln] read counts differ in this batch: R1 %d, R2 %d; using %d%s\n",
					n_in[0], n_in[1], n_seqs, n_seqs == n_in[1]? " (as bwa sampe)" : " (the smaller count; bwa sampe would read past R1)");
		}
		tot += n_seqs;
		// each file is its own batch: the gap-open cap follows that file's longest read
		aln_search_batch(bwt, n_in[0], seqs[0], &o->gap, o->gap.n_threads);
		aln_search_batch(bwt, n_in[1], seqs[1], &o->gap, o->gap.n_threads);
		// bwa_cal_pac_pos_pe, stage by stage: the drand48 pass serially in stock order, then
		// positions in parallel, the insert-size barrier, then the per-pair work in parallel
		bwa_pe_aln2seq(n_seqs, seqs);
		pe_ctx_t c;
		c.o = o, c.bwt = bwt, c.bns = bns, c.pac = pac, c.seqs[0] = seqs[0], c.seqs[1] = seqs[1];
		c.fb = &fb;
		c.cnt_chg.assign(n_threads, 0);
		c.sw_cnt.assign(4 * n_threads, 0);
		kt_for_nt(pe_pos_worker, &c, n_seqs, n_threads);
		isize_info_t ii;
		bwa_pe_isize(bwt, n_seqs, seqs, &ii, &o->pe, &last_ii);
		c.ii = &ii;
		fb.reset();
		kt_for_nt(pe_post_worker, &c, n_seqs, n_threads);
		long long cnt_chg = 0;
		uint64_t n_tot[2] = {0, 0}, n_mapped[2] = {0, 0};
		for (int t = 0; t < n_threads; ++t) {
			cnt_chg += c.cnt_chg[t];
			n_tot[0] += c.sw_cnt[4 * t], n_tot[1] += c.sw_cnt[4 * t + 1];
			n_mapped[0] += c.sw_cnt[4 * t + 2], n_mapped[1] += c.sw_cnt[4 * t + 3];
		}
		if (bwa_verbose >= 3) fprintf(stderr, "[M::aln] changed the coordinates of %lld alignments.\n", cnt_chg);
		if (o->pe.is_sw && ii.avg >= 0.0) bwa_paired_sw_report(n_tot, n_mapped);
		int stop;
		if (fb.write(out, &stop) != 0) ret = 1;
		if (stop >= 0) {
			fprintf(stderr, "[E::aln] paired reads have different names: \"%s\", \"%s\"\n", seqs[0][stop].name, seqs[1][stop].name);
			ret = 1;
		}
		bwa_free_read_seq(n_in[0], seqs[0]);
		bwa_free_read_seq(n_in[1], seqs[1]);
		if (ret) break;
		last_ii = ii;
		if (bwa_verbose >= 3) fprintf(stderr, "[M::aln] %lld read pairs have been processed.\n", tot);
	}
	bwa_seq_close(ks[0]);
	bwa_seq_close(ks[1]);
	return ret;
}

} // namespace bwa_aln

using namespace bwa_aln;

static int aln_usage(const aln_opt_t *o, FILE *fp)
{
	const gap_opt_t *g = &o->gap;
	fprintf(fp, "\n");
	fprintf(fp, "Usage:   bwa-mem3 aln [options] <prefix> <in.fq> [<in2.fq>]\n\n");
	fprintf(fp, "Aligns short reads with bwa's backtracking algorithm and writes SAM (or BAM) in one\n");
	fprintf(fp, "step: bwa aln followed by bwa samse (one FASTQ/FASTA, or '-' for stdin) or bwa sampe\n");
	fprintf(fp, "(two files), with no .sai files. Output matches bwa 0.7.19. <prefix> is an index made\n");
	fprintf(fp, "by `bwa-mem3 index` (also one made with --meth, or one staged by `bwa-mem3 shm`).\n");
	fprintf(fp, "Each option shows its default and the stock command it comes from.\n\n");
	fprintf(fp, "Search options (from bwa aln):\n");
	fprintf(fp, "  -n NUM    max #diff (int) or missing prob under %.2f err rate (float) [%.2f]\n", BWA_AVG_ERR, g->fnr);
	fprintf(fp, "  -o INT    maximum number or fraction of gap opens [%d]\n", g->max_gapo);
	fprintf(fp, "  -e INT    maximum number of gap extensions, -1 for disabling long gaps [-1]\n");
	fprintf(fp, "  -i INT    do not put an indel within INT bp towards the ends [%d]\n", g->indel_end_skip);
	fprintf(fp, "  -d INT    maximum occurrences for extending a long deletion [%d]\n", g->max_del_occ);
	fprintf(fp, "  -l INT    seed length [%d]\n", g->seed_len);
	fprintf(fp, "  -k INT    maximum differences in the seed [%d]\n", g->max_seed_diff);
	fprintf(fp, "  -m INT    maximum entries in the queue [%d]\n", g->max_entries);
	fprintf(fp, "  -t INT    number of threads [%d]\n", g->n_threads);
	fprintf(fp, "  -M INT    mismatch penalty [%d]\n", g->s_mm);
	fprintf(fp, "  -O INT    gap open penalty [%d]\n", g->s_gapo);
	fprintf(fp, "  -E INT    gap extension penalty [%d]\n", g->s_gape);
	fprintf(fp, "  -R INT    stop searching when there are >INT equally best hits [%d]\n", g->max_top2);
	fprintf(fp, "  -q INT    quality threshold for read trimming down to %dbp [%d]\n", BWA_MIN_RDLEN, g->trim_qual);
	fprintf(fp, "  -B INT    length of barcode [0]\n");
	fprintf(fp, "  -L        log-scaled gap penalty for long deletions\n");
	fprintf(fp, "  -N        non-iterative mode: search for all n-difference hits (slooow)\n");
	fprintf(fp, "  -I        the input is in the Illumina 1.3+ FASTQ-like format\n");
	fprintf(fp, "  -Y        filter Casava-filtered sequences\n");
	fprintf(fp, "\nSingle-end output options (from bwa samse):\n");
	fprintf(fp, "  --max-xa-hits INT   maximum number of alternative hits in the XA tag [%d] (samse -n)\n", o->max_xa_hits);
	fprintf(fp, "  -r STR    read group header line such as '@RG\\tID:foo\\tSM:bar' [null] (samse and sampe -r)\n");
	fprintf(fp, "  -f FILE   file to write output to instead of stdout (samse and sampe -f)\n");
	fprintf(fp, "\nPaired-end options (two input files; from bwa sampe):\n");
	fprintf(fp, "  -a INT    maximum insert size [%d]\n", o->pe.max_isize);
	fprintf(fp, "  --max-occ INT   maximum occurrences for one end [%d] (sampe -o)\n", o->pe.max_occ);
	fprintf(fp, "  --max-xa-hits-discordant INT   maximum hits to output for discordant pairs [%d] (sampe -N)\n", o->pe.N_multi);
	fprintf(fp, "  -c FLOAT  prior of chimeric rate (lower bound) [%.1le]\n", o->pe.ap_prior);
	fprintf(fp, "  -s        disable Smith-Waterman for the unmapped mate\n");
	fprintf(fp, "  -A        disable insert size estimate (force -s)\n");
	fprintf(fp, "  -P        accepted and ignored (sampe -P loads the index into memory; always so here)\n");
	fprintf(fp, "\nNew in bwa-mem3 aln:\n");
	fprintf(fp, "  --bam     write BAM instead of SAM\n");
	fprintf(fp, "  --bam-level INT     BAM compression level 0-9 [1]\n");
	fprintf(fp, "  --help    show this message\n");
	fprintf(fp, "\nThe long names stand in for stock samse/sampe -n and -o and sampe -N, whose letters\n");
	fprintf(fp, "are taken by the search options above.\n");
	fprintf(fp, "\n");
	return fp == stdout ? 0 : 1;
}

/* An input that stock would fail to open must fail before any output is created or the index
 * is loaded. "-" is standard input. */
static int aln_input_ok(const char *fn)
{
	struct stat st;
	memset(&st, 0, sizeof st);
	if (strcmp(fn, "-") == 0) return 1;
	if (stat(fn, &st) != 0) {
		fprintf(stderr, "[E::aln] cannot open the input '%s': %s\n", fn, strerror(errno));
		return 0;
	}
	if (S_ISDIR(st.st_mode)) {
		fprintf(stderr, "[E::aln] cannot open the input '%s': is a directory\n", fn);
		return 0;
	}
	if (access(fn, R_OK) != 0) {
		fprintf(stderr, "[E::aln] cannot open the input '%s': %s\n", fn, strerror(errno));
		return 0;
	}
	return 1;
}

/* Checks that -f can be written without leaving a file behind (a file created for the check
 * is removed again). Standard output needs no check. */
static int aln_output_ok(const char *path)
{
	if (!path || strcmp(path, "-") == 0) return 1;
	int existed = access(path, F_OK) == 0;
	FILE *fp = fopen(path, "a");
	if (!fp) {
		fprintf(stderr, "[E::aln] cannot open '%s' for writing: %s\n", path, strerror(errno));
		return 0;
	}
	fclose(fp);
	if (!existed) unlink(path);
	return 1;
}

enum { OPT_MAX_XA = 256, OPT_MAX_XA_DISC, OPT_MAX_OCC, OPT_BAM, OPT_BAM_LEVEL, OPT_HELP };

int main_aln(int argc, char *argv[])
{
	aln_opt_t o;
	aln_opt_init(&o);
	gap_opt_t *opt = &o.gap;
	int c, opte = -1, bam = 0, bam_level = 1, pe_flag = 0;
	const char *out_path = NULL;
	char *rg_line = NULL;
	static const struct option long_opts[] = {
		{"max-xa-hits", required_argument, 0, OPT_MAX_XA},
		{"max-xa-hits-discordant", required_argument, 0, OPT_MAX_XA_DISC},
		{"max-occ", required_argument, 0, OPT_MAX_OCC},
		{"bam", no_argument, 0, OPT_BAM},
		{"bam-level", required_argument, 0, OPT_BAM_LEVEL},
		{"help", no_argument, 0, OPT_HELP},
		{0, 0, 0, 0}
	};

	optind = 1;
	while ((c = getopt_long(argc, argv, "n:o:e:i:d:l:k:LR:m:t:NM:O:E:q:f:r:IYB:a:c:sAP", long_opts, 0)) >= 0) {
		switch (c) {
		case 'n':
			if (strstr(optarg, ".")) opt->fnr = atof(optarg), opt->max_diff = -1;
			else opt->max_diff = atoi(optarg), opt->fnr = -1.0;
			break;
		case 'o': opt->max_gapo = atoi(optarg); break;
		case 'e': opte = atoi(optarg); break;
		case 'M': opt->s_mm = atoi(optarg); break;
		case 'O': opt->s_gapo = atoi(optarg); break;
		case 'E': opt->s_gape = atoi(optarg); break;
		case 'd': opt->max_del_occ = atoi(optarg); break;
		case 'i': opt->indel_end_skip = atoi(optarg); break;
		case 'l': opt->seed_len = atoi(optarg); break;
		case 'k': opt->max_seed_diff = atoi(optarg); break;
		case 'm': opt->max_entries = atoi(optarg); break;
		case 't': opt->n_threads = atoi(optarg); break;
		case 'L': opt->mode |= BWA_MODE_LOGGAP; break;
		case 'R': opt->max_top2 = atoi(optarg); break;
		case 'q': opt->trim_qual = atoi(optarg); break;
		case 'N': opt->mode |= BWA_MODE_NONSTOP; opt->max_top2 = 0x7fffffff; break;
		case 'f': out_path = optarg; break;
		case 'r':
			// stock bwa_set_rg rejects a literal tab; bwa-mem3's bwa_set_rg does not
			if (strchr(optarg, '\t') != NULL) {
				fprintf(stderr, "[E::%s] the read group line contained literal <tab> characters -- replace with escaped tabs: \\t\n", __func__);
				return 1;
			}
			free(rg_line);
			if ((rg_line = bwa_set_rg(optarg)) == 0) return 1;
			break;
		case 'I': opt->mode |= BWA_MODE_IL13; break;
		case 'Y': opt->mode |= BWA_MODE_CFY; break;
		case 'B': opt->mode |= atoi(optarg) << 24; break;
		case 'a': o.pe.max_isize = atoi(optarg); pe_flag = 1; break;
		case 'c': o.pe.ap_prior = atof(optarg); pe_flag = 1; break;
		case 's': o.pe.is_sw = 0; pe_flag = 1; break;
		case 'A': o.pe.force_isize = 1; pe_flag = 1; break;
		case 'P': o.pe.is_preload = 1; pe_flag = 1; break;
		case OPT_MAX_XA: o.max_xa_hits = atoi(optarg); o.pe.n_multi = o.max_xa_hits; break;
		case OPT_MAX_XA_DISC: o.pe.N_multi = atoi(optarg); pe_flag = 1; break;
		case OPT_MAX_OCC: o.pe.max_occ = atoi(optarg); pe_flag = 1; break;
		case OPT_BAM: bam = 1; break;
		case OPT_BAM_LEVEL: bam_level = atoi(optarg); break;
		case OPT_HELP: return aln_usage(&o, stdout);
		default: return 1;
		}
	}
	if (opte > 0) {
		opt->max_gape = opte;
		opt->mode &= ~BWA_MODE_GAPE;
	}
	if (opt->n_threads < 1) opt->n_threads = 1;

	if (optind + 2 > argc) return aln_usage(&o, stderr);
	if (optind + 4 <= argc) return aln_usage(&o, stderr);
	const int paired = optind + 3 == argc;
	if (pe_flag && !paired)
		fprintf(stderr, "[W::%s] paired-end options (-a, -c, -s, -A, -P, --max-occ, --max-xa-hits-discordant) are ignored with one FASTQ\n", __func__);

	const char *prefix = argv[optind], *fn_fq = argv[optind + 1];
	{
		std::string probe = std::string(prefix) + ".bwt.2bit.64";
		if (access(probe.c_str(), F_OK) != 0) {
			fprintf(stderr, "[E::%s] fail to locate the index '%s' (no %s); run `bwa-mem3 index` first\n",
					__func__, prefix, probe.c_str());
			return 1;
		}
	}
	if (!aln_input_ok(fn_fq) || (paired && !aln_input_ok(argv[optind + 2])) || !aln_output_ok(out_path))
		return 1;
	FMI_search *fmi = new FMI_search(prefix);
	fmi->load_index(/*load_pac=*/true, opt->n_threads);
	const bntseq_t *bns = fmi->idx->bns;
	bwa_aln::bwt_t *bwt = bwt_from_fmi(fmi);

	kstring_t hdr = {0, 0, 0};
	aln_print_sam_hdr(&hdr, bns, rg_line, bwa_pg);
	aln_out_t *out = aln_out_open(out_path, bam, bam_level, bns, hdr.s);
	free(hdr.s);
	int ret = 1;
	if (out) {
		ret = paired ? aln_run_pe(&o, bwt, bns, fmi->idx->pac, fn_fq, argv[optind + 2], out)
					 : aln_run_se(&o, bwt, bns, fmi->idx->pac, fn_fq, out);
	}
	if (out && aln_out_close(out) != 0) ret = 1;

	bwt_destroy(bwt);
	delete fmi;
	free(rg_line);
	return ret;
}
