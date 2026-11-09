/* SPDX-License-Identifier: MIT */
/* `aln --known-snvs`: loading and checking the known SNV sites. New code. */
// htslib first: its kstring.h and bwa's share an include guard, and vcf.h needs htslib's.
#include "htslib/hts.h"
#include "htslib/vcf.h"

#include "aln_known_snvs.h"
#include "aln_bwt.h"
#include "aln_meth.h"
#include "aln_meth_realign.h"
#include "aln_alloc.h"
#include "bwase.h"
#include "bwtgap.h"

#include "../FMI_search.h"
#include "../bwa.h"

#include <algorithm>
#include <cerrno>
#include <climits>
#include <cstdlib>
#include <dirent.h>
#include <unistd.h>
#include <cstdio>
#include <cstring>
#include <string>
#include <unordered_map>
#include <getopt.h>

int bwa_index(int argc, char *argv[]); // src/bwtindex.cpp, the `index` command

namespace bwa_aln {

namespace {

int base_code(const char *s)
{
	if (!s || !s[0] || s[1]) return -1; // exactly one base
	switch (s[0]) {
	case 'A': case 'a': return 0;
	case 'C': case 'c': return 1;
	case 'G': case 'g': return 2;
	case 'T': case 't': return 3;
	default: return -1;
	}
}

struct vcf_closer {
	htsFile *fp = nullptr;
	bcf_hdr_t *hdr = nullptr;
	bcf1_t *rec = nullptr;
	~vcf_closer()
	{
		if (rec) bcf_destroy(rec);
		if (hdr) bcf_hdr_destroy(hdr);
		if (fp) hts_close(fp);
	}
};

} // namespace

int ks_load(const char *path, const bntseq_t *bns, const uint8_t *pac, ks_sites_t *out)
{
	*out = ks_sites_t();
	out->by_contig.assign(bns->n_seqs, {});
	std::unordered_map<std::string, int> rid_of;
	rid_of.reserve((size_t)bns->n_seqs * 2);
	for (int i = 0; i < bns->n_seqs; ++i) rid_of.emplace(bns->anns[i].name, i);

	vcf_closer v;
	v.fp = hts_open(path, "r");
	if (!v.fp) {
		fprintf(stderr, "[E::aln] --known-snvs: cannot open '%s'\n", path);
		return -1;
	}
	v.hdr = bcf_hdr_read(v.fp);
	if (!v.hdr) {
		fprintf(stderr, "[E::aln] --known-snvs: '%s' is not a VCF or BCF file\n", path);
		return -1;
	}
	v.rec = bcf_init();
	int r;
	while ((r = bcf_read(v.fp, v.hdr, v.rec)) == 0) {
		bcf1_t *rec = v.rec;
		if (bcf_unpack(rec, BCF_UN_STR) < 0) {
			fprintf(stderr, "[E::aln] --known-snvs: '%s': cannot parse record %lld\n", path,
					(long long)out->n_records + 1);
			return -1;
		}
		++out->n_records;
		const char *chrom = bcf_seqname_safe(v.hdr, rec);
		const long long pos1 = (long long)rec->pos + 1;
		auto it = rid_of.find(chrom);
		if (it == rid_of.end()) {
			fprintf(stderr, "[E::aln] --known-snvs: '%s': %s:%lld: contig '%s' is not in the reference\n",
					path, chrom, pos1, chrom);
			return -1;
		}
		const int rid = it->second;
		const bntann1_t &ann = bns->anns[rid];
		const char *ref_s = rec->n_allele > 0 ? rec->d.allele[0] : nullptr;
		const int ref = base_code(ref_s);
		if (ref < 0) { // REF longer than one base or not ACGT: no allele of it is an SNV
			out->n_not_snv += rec->n_allele > 1 ? rec->n_allele - 1 : 1;
			continue;
		}
		if (rec->pos < 0 || rec->pos >= ann.len) {
			fprintf(stderr, "[E::aln] --known-snvs: '%s': %s:%lld is past the end of the contig (%d bp)\n",
					path, chrom, pos1, ann.len);
			return -1;
		}
		const int64_t pos_f = ann.offset + rec->pos;
		int amb_rid;
		if (bns_cnt_ambi(bns, pos_f, 1, &amb_rid) > 0) {
			fprintf(stderr, "[E::aln] --known-snvs: '%s': %s:%lld: REF %s, but the reference base is N\n",
					path, chrom, pos1, ref_s);
			return -1;
		}
		const int ref_base = (pac[pos_f >> 2] >> ((~pos_f & 3) << 1)) & 3;
		if (ref_base != ref) {
			fprintf(stderr, "[E::aln] --known-snvs: '%s': %s:%lld: REF %s does not match the reference base %c\n",
					path, chrom, pos1, ref_s, "ACGT"[ref_base]);
			return -1;
		}
		for (int a = 1; a < rec->n_allele; ++a) {
			const int alt = base_code(rec->d.allele[a]);
			if (alt < 0 || alt == ref) {
				++out->n_not_snv;
				continue;
			}
			out->by_contig[rid].push_back({rec->pos, (uint8_t)ref, (uint8_t)alt});
		}
	}
	if (r < -1) {
		fprintf(stderr, "[E::aln] --known-snvs: '%s': read error after record %lld\n", path,
				(long long)out->n_records);
		return -1;
	}
	for (std::vector<ks_site_t> &sites : out->by_contig) {
		std::sort(sites.begin(), sites.end(), [](const ks_site_t &x, const ks_site_t &y) {
			return x.pos != y.pos ? x.pos < y.pos : x.alt < y.alt;
		});
		const size_t before = sites.size();
		sites.erase(std::unique(sites.begin(), sites.end(), [](const ks_site_t &x, const ks_site_t &y) {
			return x.pos == y.pos && x.alt == y.alt;
		}), sites.end());
		out->n_duplicate += (int64_t)(before - sites.size());
		out->n_sites += (int64_t)sites.size();
	}
	return 0;
}

namespace {

void remove_dir(const std::string &dir)
{
	DIR *d = opendir(dir.c_str());
	if (d) {
		struct dirent *e;
		while ((e = readdir(d)) != nullptr) {
			if (strcmp(e->d_name, ".") == 0 || strcmp(e->d_name, "..") == 0) continue;
			unlink((dir + "/" + e->d_name).c_str());
		}
		closedir(d);
	}
	rmdir(dir.c_str());
}

} // namespace

int ks_windows_build(const ks_sites_t *sites, const bntseq_t *bns, const uint8_t *pac, int n_threads,
					 int meth, ks_windows_t *out)
{
	ks_windows_free(out);
	for (int rid = 0; rid < (int)sites->by_contig.size(); ++rid) {
		const bntann1_t &ann = bns->anns[rid];
		for (const ks_site_t &s : sites->by_contig[rid]) {
			const int64_t beg = std::max<int64_t>(0, s.pos - KS_WINDOW_HALF);
			const int64_t end = std::min<int64_t>(ann.len, s.pos + KS_WINDOW_HALF + 1);
			out->win.push_back({rid, beg, (int32_t)(end - beg), s.pos, s.alt});
		}
	}
	const char *tmp = getenv("TMPDIR");
	std::string templ = std::string(tmp && *tmp ? tmp : "/tmp") + "/bwa-mem3-known-snvs.XXXXXX";
	std::vector<char> buf(templ.begin(), templ.end());
	buf.push_back('\0');
	if (!mkdtemp(buf.data())) {
		fprintf(stderr, "[E::aln] --known-snvs: cannot create a temporary directory from '%s': %s\n",
				templ.c_str(), strerror(errno));
		return -1;
	}
	out->dir = buf.data();
	const std::string fa = out->dir + "/windows.fa";
	FILE *fp = fopen(fa.c_str(), "w");
	bool ok = fp != nullptr;
	std::string seq;
	for (size_t i = 0; ok && i < out->win.size(); ++i) {
		const ks_window_t &w = out->win[i];
		const int64_t off = bns->anns[w.rid].offset;
		seq.resize((size_t)w.len);
		for (int32_t j = 0; j < w.len; ++j) {
			const int64_t k = off + w.beg + j;
			seq[(size_t)j] = "ACGT"[(pac[k >> 2] >> ((~k & 3) << 1)) & 3];
		}
		seq[(size_t)(w.site_pos - w.beg)] = "ACGT"[w.alt];
		ok = fprintf(fp, ">w%zu\n%s\n", i, seq.c_str()) >= 0;
	}
	if (fp && fclose(fp) != 0) ok = false;
	if (!ok) {
		fprintf(stderr, "[E::aln] --known-snvs: cannot write '%s'\n", fa.c_str());
		ks_windows_free(out);
		return -1;
	}
	if (meth) { // `index --meth`: the converted index, and the original one it is checked against
		std::string a0 = "index", a1 = "--meth", a2 = fa;
		char *argv[] = {&a0[0], &a1[0], &a2[0], nullptr};
#if defined(__APPLE__) || defined(__FreeBSD__) || defined(__OpenBSD__) || defined(__NetBSD__)
		optind = 1;
		optreset = 1;
#else
		optind = 0;
#endif
		if (bwa_index(3, argv) != 0) {
			fprintf(stderr, "[E::aln] --known-snvs: cannot index the ALT windows in '%s'\n", out->dir.c_str());
			ks_windows_free(out);
			return -1;
		}
		out->fmi = new FMI_search((fa + ".meth").c_str());
		out->fmi->load_index(/*load_pac=*/false, n_threads);
		out->bwt = bwt_from_fmi(out->fmi);
		out->mref = new aln_meth_ref_t();
		if (aln_meth_ref_load(out->mref, fa.c_str(), out->fmi->idx->bns, n_threads) != 0) {
			ks_windows_free(out);
			return -1;
		}
		out->bns = out->mref->bns;
	} else {
		if (bwa_idx_build(fa.c_str(), fa.c_str()) != 0) {
			fprintf(stderr, "[E::aln] --known-snvs: cannot index the ALT windows in '%s'\n", out->dir.c_str());
			ks_windows_free(out);
			return -1;
		}
		out->fmi = new FMI_search(fa.c_str());
		out->fmi->load_index(/*load_pac=*/true, n_threads);
		out->bns = out->fmi->idx->bns;
		out->bwt = bwt_from_fmi(out->fmi);
	}
	if (out->bns->n_seqs != (int)out->win.size()) {
		fprintf(stderr, "[E::aln] --known-snvs: the window index has %d contigs for %zu windows\n",
				out->bns->n_seqs, out->win.size());
		ks_windows_free(out);
		return -1;
	}
	return 0;
}

void ks_windows_free(ks_windows_t *w)
{
	if (w->mref) {
		aln_meth_ref_free(w->mref);
		delete w->mref;
	}
	if (w->bwt) bwt_destroy(w->bwt);
	delete w->fmi;
	if (!w->dir.empty()) remove_dir(w->dir);
	*w = ks_windows_t();
}

namespace {

constexpr int KS_MAX_LOCI = 256; // hits expanded per read; above it the stock result stands

struct ks_hit_t {
	int64_t pos;       // forward-strand start in the main index's packed coordinates
	int strand;
	int score, n_mm, n_gapo, n_gape, ref_shift;
	int64_t alt_site;  // packed coordinate of the site whose ALT this hit's cost assumes, or -1
	int diffs() const { return n_mm + n_gapo + n_gape; }
};

uint64_t mix64(uint64_t x)
{
	x += 0x9e3779b97f4a7c15ULL;
	x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
	x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
	return x ^ (x >> 31);
}

uint64_t read_seed(const bntseq_t *bns, const bwa_seq_t *s)
{
	uint64_t h = mix64((uint64_t)bns->seed);
	for (const char *c = s->name; c && *c; ++c) h = mix64(h ^ (uint8_t)*c);
	for (int i = 0; i < s->len; ++i) h = mix64(h ^ s->seq[i]);
	return h;
}

ks_hit_t hit_of(const bwt_aln1_t &q, int64_t pos, int strand, int64_t alt_site)
{
	return {pos, strand, (int)q.score, (int)q.n_mm, (int)q.n_gapo, (int)q.n_gape,
			(int)q.n_del - (int)q.n_ins, alt_site};
}

int read_max_diff(const gap_opt_t *opt, int len)
{
	return opt->fnr > 0.0 ? bwa_cal_maxdiff(len, BWA_AVG_ERR, opt->fnr) : opt->max_diff;
}

/* The merged hit set of one read: kept hits sorted best first, c1 of them at the best score. */
struct ks_merged_t {
	std::vector<ks_hit_t> kept;
	int c1 = 0;
	int max_diff = 0;
};

/* Merges a read's main hits (aln/n_aln, its main search) with its window hits (win_aln/win_n).
 * Returns KS_UNTOUCHED (no window hit covers its site), KS_SAME (they change no hit), KS_INCOMPLETE
 * (too many loci or too long for the window), or KS_CREDITED with `m` filled. */
ks_outcome_t ks_merge(const ks_windows_t *w, const bntseq_t *bns, const bwt_t *bwt, const gap_opt_t *opt,
					  int len, const bwt_aln1_t *aln, int n_aln, const bwt_aln1_t *win_aln, int win_n,
					  ks_merged_t *m)
{
	if (win_n <= 0) return KS_UNTOUCHED;
	// Window hits that cover their window's site, in main-index coordinates.
	std::vector<ks_hit_t> win_hits;
	int64_t n_rows = 0;
	for (int k = 0; k < win_n; ++k) n_rows += (int64_t)(win_aln[k].l - win_aln[k].k + 1);
	if (n_rows > KS_MAX_LOCI) return KS_INCOMPLETE;
	for (int k = 0; k < win_n; ++k) {
		const bwt_aln1_t &q = win_aln[k];
		const int span = len + (int)q.n_del - (int)q.n_ins;
		for (bwtint_t r = q.k; r <= q.l; ++r) {
			int strand;
			const bwtint_t wp = bwa_sa2pos(w->bns, w->bwt, r, span, &strand);
			if (wp == (bwtint_t)-1) continue;
			const int wid = bns_pos2rid(w->bns, (int64_t)wp);
			if (wid < 0) continue;
			const ks_window_t &win = w->win[wid];
			const int64_t off = (int64_t)wp - w->bns->anns[wid].offset;
			if (off < 0 || off + span > win.len) continue; // runs into the next window
			const int64_t beg = win.beg + off;
			if (win.site_pos < beg || win.site_pos >= beg + span) continue; // misses its site
			const int64_t base = bns->anns[win.rid].offset;
			win_hits.push_back(hit_of(q, base + beg, strand, base + win.site_pos));
		}
	}
	if (win_hits.empty()) return KS_UNTOUCHED;
	if (len + opt->max_gape > KS_WINDOW_HALF) return KS_INCOMPLETE;

	// The cuts (aln's: score <= best + mismatch penalty, differences <= min(max_diff, best's + 1))
	// come from the scores alone, so only the main search's intervals that pass them are expanded:
	// crediting often drops hundreds of REF copies that sat one difference behind the stock best.
	const int max_diff = read_max_diff(opt, len);
	int best = INT32_MAX;
	for (const ks_hit_t &h : win_hits) best = std::min(best, h.score);
	for (int k = 0; k < n_aln; ++k) best = std::min(best, (int)aln[k].score);
	int best_diffs = INT32_MAX;
	for (const ks_hit_t &h : win_hits)
		if (h.score == best) best_diffs = std::min(best_diffs, h.diffs());
	for (int k = 0; k < n_aln; ++k)
		if ((int)aln[k].score == best) best_diffs = std::min(best_diffs, (int)(aln[k].n_mm + aln[k].n_gapo + aln[k].n_gape));
	const int cap = std::max(std::min(max_diff, best_diffs + 1), best_diffs);
	auto passes = [&](int score, int diffs) { return score <= best + opt->s_mm && diffs <= cap; };
	std::vector<ks_hit_t> hits;
	n_rows = 0;
	for (int k = 0; k < n_aln; ++k)
		if (passes((int)aln[k].score, (int)(aln[k].n_mm + aln[k].n_gapo + aln[k].n_gape)))
			n_rows += (int64_t)(aln[k].l - aln[k].k + 1);
	if (n_rows > KS_MAX_LOCI) return KS_INCOMPLETE;
	for (int k = 0; k < n_aln; ++k) {
		const bwt_aln1_t &q = aln[k];
		if (!passes((int)q.score, (int)(q.n_mm + q.n_gapo + q.n_gape))) continue;
		const int span = len + (int)q.n_del - (int)q.n_ins;
		for (bwtint_t r = q.k; r <= q.l; ++r) {
			int strand;
			const bwtint_t p = bwa_sa2pos(bns, bwt, r, span, &strand);
			if (p == (bwtint_t)-1) continue;
			hits.push_back(hit_of(q, (int64_t)p, strand, -1));
		}
	}
	// Merge: per (strand, start) the lower score wins.
	bool changed = false;
	for (const ks_hit_t &h : win_hits) {
		if (!passes(h.score, h.diffs())) continue;
		auto it = std::find_if(hits.begin(), hits.end(), [&](const ks_hit_t &x) {
			return x.strand == h.strand && x.pos == h.pos;
		});
		if (it == hits.end()) {
			hits.push_back(h);
			changed = true;
		} else if (h.score < it->score) {
			*it = h;
			changed = true;
		}
	}
	if (!changed) return KS_SAME;
	std::stable_sort(hits.begin(), hits.end(), [](const ks_hit_t &a, const ks_hit_t &b) {
		return a.score != b.score ? a.score < b.score : a.pos != b.pos ? a.pos < b.pos : a.strand < b.strand;
	});
	int c1 = 0;
	while (c1 < (int)hits.size() && hits[c1].score == best) ++c1;
	m->kept = std::move(hits);
	m->c1 = c1;
	m->max_diff = max_diff;
	return KS_CREDITED;
}

int approx_mapq(int c1, int c2, int n_mm, int max_diff)
{
	bwa_seq_t t;
	memset(&t, 0, sizeof t);
	t.c1 = c1; t.c2 = c2; t.n_mm = n_mm;
	return bwa_approx_mapQ(&t, max_diff);
}

} // namespace

ks_outcome_t ks_override(const ks_windows_t *w, const bntseq_t *bns, const bwt_t *bwt, const gap_opt_t *opt,
						 int max_xa_hits, bwa_seq_t *s, const bwt_aln1_t *win_aln, int win_n, int64_t *alt_site)
{
	*alt_site = -1;
	ks_merged_t m;
	const ks_outcome_t r = ks_merge(w, bns, bwt, opt, s->len, s->aln, s->n_aln, win_aln, win_n, &m);
	if (r != KS_CREDITED) return r;
	std::vector<ks_hit_t> &kept = m.kept;
	const int c1 = m.c1;
	// The stock primary stays when it is still among the best; otherwise a per-read draw.
	int pick = -1;
	if (s->type == BWA_TYPE_UNIQUE || s->type == BWA_TYPE_REPEAT)
		for (int i = 0; i < c1 && pick < 0; ++i)
			if (kept[i].pos == (int64_t)s->pos && kept[i].strand == (int)s->strand) pick = i;
	if (pick < 0) pick = (int)(read_seed(bns, s) % (uint64_t)c1);
	const ks_hit_t p = kept[pick];

	s->pos = (bwtint_t)p.pos;
	s->strand = p.strand;
	s->n_mm = p.n_mm; s->n_gapo = p.n_gapo; s->n_gape = p.n_gape;
	s->ref_shift = p.ref_shift;
	s->score = p.score;
	s->c1 = c1;
	s->c2 = kept.size() - c1;
	s->type = c1 > 1 ? BWA_TYPE_REPEAT : BWA_TYPE_UNIQUE;
	s->seQ = s->mapQ = bwa_approx_mapQ(s, m.max_diff);
	*alt_site = p.alt_site;

	if (s->multi) free(s->multi);
	s->multi = 0; s->n_multi = 0;
	if (max_xa_hits > 0 && (int)kept.size() <= max_xa_hits + 1 && kept.size() > 1) {
		s->multi = (bwt_multi1_t *)aln_calloc(kept.size() - 1, sizeof(bwt_multi1_t));
		for (int i = 0; i < (int)kept.size(); ++i) {
			if (i == pick) continue;
			bwt_multi1_t &x = s->multi[s->n_multi++];
			x.pos = (bwtint_t)kept[i].pos;
			x.strand = kept[i].strand;
			x.gap = kept[i].n_gapo + kept[i].n_gape;
			x.mm = kept[i].n_mm;
			x.ref_shift = kept[i].ref_shift;
		}
	}
	return KS_CREDITED;
}

int ks_twin_place(const ks_windows_t *w, const bntseq_t *bns, const bwt_t *bwt, const gap_opt_t *opt,
				  const bwa_seq_t *t, const bwt_aln1_t *win_aln, int win_n, int64_t pos, int strand,
				  int *mapq)
{
	*mapq = 0;
	ks_merged_t m;
	const ks_outcome_t r = ks_merge(w, bns, bwt, opt, t->len, t->aln, t->n_aln, win_aln, win_n, &m);
	if (r == KS_INCOMPLETE) return -1;
	if (r == KS_CREDITED) {
		const ks_hit_t &b = m.kept[0];
		*mapq = approx_mapq(m.c1, (int)m.kept.size() - m.c1, b.n_mm, m.max_diff);
		if (m.c1 > 1) return 0; // MAPQ 0: a tie never passes
		return b.strand == strand && std::llabs(b.pos - pos) <= KS_TWIN_SLOP;
	}
	// The stock result, without a draw: c1 rows at the best score, c2 the rest.
	if (t->n_aln == 0) return 0;
	int64_t c1 = 0, c2 = 0;
	for (int k = 0; k < t->n_aln; ++k) {
		const int64_t n = (int64_t)(t->aln[k].l - t->aln[k].k + 1);
		if (t->aln[k].score == t->aln[0].score) c1 += n; else c2 += n;
	}
	const int max_diff = read_max_diff(opt, t->len);
	*mapq = approx_mapq((int)std::min<int64_t>(c1, 1 << 27), (int)std::min<int64_t>(c2, 1 << 27), t->aln[0].n_mm, max_diff);
	if (c1 != 1) return 0;
	const bwt_aln1_t &q = t->aln[0];
	int s;
	const bwtint_t p = bwa_sa2pos(bns, bwt, q.k, t->len + (int)q.n_del - (int)q.n_ins, &s);
	return p != (bwtint_t)-1 && s == strand && std::llabs((int64_t)p - pos) <= KS_TWIN_SLOP;
}

int ks_site_votes(const ks_sites_t *sites, const bntseq_t *bns, const gap_opt_t *opt, const bwa_seq_t *s,
				  int64_t alt_site, int hyp, std::vector<ks_vote_t> *out)
{
	out->clear();
	if (s->type != BWA_TYPE_UNIQUE && s->type != BWA_TYPE_REPEAT) return -1;
	const int rid = bns_pos2rid(bns, (int64_t)s->pos);
	if (rid < 0) return -1;
	const int64_t cbeg = (int64_t)s->pos - bns->anns[rid].offset;
	const int64_t span = pos_end(s) - (int64_t)s->pos;
	const std::vector<ks_site_t> &lst = sites->by_contig[rid];
	auto it = std::lower_bound(lst.begin(), lst.end(), cbeg, [](const ks_site_t &x, int64_t p) { return x.pos < p; });
	if (it == lst.end() || it->pos >= cbeg + span) return -1;
	// The read in reference orientation (s->seq is in read order after refinement).
	std::vector<uint8_t> rs(s->len);
	for (int i = 0; i < s->len; ++i) rs[i] = s->strand ? (s->seq[s->len - 1 - i] > 3 ? 4 : 3 - s->seq[s->len - 1 - i]) : s->seq[i];
	const int diffs = s->n_mm + s->n_gapo + s->n_gape;
	const int max_diff = read_max_diff(opt, s->len);
	// Under --meth a converted base (C->T, or G->A) also reads as its source.
	const int source = hyp == 1 ? 1 : 2, target = hyp == 1 ? 3 : 0;
	auto in_set = [&](uint8_t b, uint8_t allele) { return b == allele || (hyp >= 0 && allele == source && b == target); };
	int first = -1;
	for (; it != lst.end() && it->pos < cbeg + span; ++it) {
		int qi = -1;
		const int64_t rpos = it->pos;
		if (!s->cigar) {
			qi = (int)(rpos - cbeg);
		} else {
			int64_t r = cbeg;
			int q = 0;
			for (int j = 0; j < s->n_cigar && qi < 0; ++j) {
				const int op = __cigar_op(s->cigar[j]), n = __cigar_len(s->cigar[j]);
				if (op == FROM_M) {
					if (rpos >= r && rpos < r + n) qi = q + (int)(rpos - r);
					q += n; r += n;
				} else if (op == FROM_I || op == FROM_S) {
					q += n;
				} else if (op == FROM_D) {
					if (rpos >= r && rpos < r + n) break;
					r += n;
				}
			}
		}
		if (qi < 0 || qi >= s->len) continue;
		ks_vote_t v;
		v.site = &*it;
		v.rid = rid;
		v.qi = qi;
		const uint8_t b = rs[qi];
		const bool is_ref = in_set(b, it->ref), is_alt = in_set(b, it->alt);
		v.allele = is_ref && is_alt ? 'H' : is_ref ? 'R' : is_alt ? 'A' : 'O';
		const bool alt_assumed = alt_site == KS_ALT_ALL || alt_site == bns->anns[rid].offset + rpos;
		v.away = diffs - (in_set(b, alt_assumed ? it->alt : it->ref) ? 0 : 1);
		v.max_diff = max_diff;
		if (first < 0 && (v.allele == 'R' || v.allele == 'A')) first = (int)out->size();
		out->push_back(v);
	}
	return first;
}

void ks_twin_seq(const bwa_seq_t *s, const ks_vote_t &v, int hyp, std::vector<uint8_t> *read_order)
{
	// s->seq is in read order after refinement; v.qi indexes the reference-oriented read.
	read_order->assign(s->seq, s->seq + s->len);
	uint8_t other = v.allele == 'R' ? v.site->alt : v.site->ref;
	if (hyp == 1 && other == 1) other = 3;      // C reads as T on a C->T strand
	else if (hyp == 0 && other == 2) other = 0; // G reads as A on a G->A strand
	if (s->strand) (*read_order)[s->len - 1 - v.qi] = 3 - other;
	else (*read_order)[v.qi] = other;
}

void ks_format_xv(const bntseq_t *bns, const std::vector<ks_vote_t> &votes, int twin_of, int twin_same,
				  int twin_mapq, std::string *out)
{
	out->clear();
	for (size_t i = 0; i < votes.size(); ++i) {
		const ks_vote_t &v = votes[i];
		char buf[64];
		if (!out->empty()) out->push_back(';');
		out->append(bns->anns[v.rid].name);
		snprintf(buf, sizeof buf, ":%lld:%c>%c,%c,%d,%d,", (long long)v.site->pos + 1, "ACGT"[v.site->ref],
				 "ACGT"[v.site->alt], v.allele, v.away, v.max_diff);
		out->append(buf);
		if ((int)i == twin_of && twin_same >= 0) snprintf(buf, sizeof buf, "%d,%d", twin_same, twin_mapq);
		else snprintf(buf, sizeof buf, "NA,NA");
		out->append(buf);
	}
}

void ks_meth_candidates(const ks_windows_t *w, const bntseq_t *bns, const aln_meth_read_t &wres,
						std::vector<aln_meth_candidate_t> *out)
{
	out->clear();
	for (const aln_meth_hit_t &h : wres.hits) {
		if (h.rid < 0 || h.rid >= (int)w->win.size()) continue;
		const ks_window_t &win = w->win[h.rid];
		const int64_t off = (int64_t)h.pos - w->bns->anns[h.rid].offset;
		if (off < 0 || off + h.span > win.len) continue;
		const int64_t beg = win.beg + off;
		if (win.site_pos < beg || win.site_pos >= beg + h.span) continue; // misses its site
		const int64_t rb = bns->anns[win.rid].offset + beg;
		out->push_back({win.rid, (int)h.strand, (int)h.hyp, rb, rb + h.span});
	}
}

} // namespace bwa_aln
