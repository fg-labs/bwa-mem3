/* SPDX-License-Identifier: MIT */
/* `aln --known-snvs`: loading and checking the known SNV sites. New code. */
// htslib first: its kstring.h and bwa's share an include guard, and vcf.h needs htslib's.
#include "htslib/hts.h"
#include "htslib/vcf.h"

#include "aln_known_snvs.h"
#include "aln_main.h"
#include "aln_bwt.h"
#include "aln_meth.h"
#include "aln_meth_realign.h"
#include "aln_alloc.h"
#include "bwase.h"
#include "bwtgap.h"

#include "../FMI_search.h"
#include "../bwa.h"
#include "../libsais_build.h"

#include <algorithm>
#include <atomic>
#include <climits>
#include <cstdlib>
#include <cerrno>
#include <cstdio>
#include <cstring>
#include <memory>
#include <string>
#include <unordered_map>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>

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
	if (bcf_hdr_set_samples(v.hdr, nullptr, 0) < 0) {
		fprintf(stderr, "[E::aln] --known-snvs: cannot exclude unused samples from '%s'\n", path);
		return -1;
	}
	v.rec = bcf_init();
	if (!v.rec) {
		fprintf(stderr, "[E::aln] --known-snvs: cannot allocate a VCF record for '%s'\n", path);
		return -1;
	}
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

// A restored-state bns (what bns_restore leaves after bwa_idx_load_ele closes fp_pac): seed 11,
// gi 0, empty anno, no holes, fp_pac NULL, pos2rid built.
bntseq_t *ks_new_bns(int n_seqs)
{
	bntseq_t *b = (bntseq_t *)aln_calloc(1, sizeof(bntseq_t));
	b->seed = 11;
	b->n_seqs = n_seqs;
	b->anns = (bntann1_t *)aln_calloc((size_t)std::max(n_seqs, 1), sizeof(bntann1_t));
	return b;
}

void ks_set_ann(bntann1_t *a, const std::string &name, int64_t offset, int32_t len)
{
	a->name = aln_strdup(name.c_str());
	a->anno = aln_strdup("");
	a->offset = offset;
	a->len = len;
}

} // namespace

void ks_windows_ref(const std::vector<ks_window_t> &win, const bntseq_t *bns, const uint8_t *pac,
					bntseq_t **out_bns, uint8_t **out_pac)
{
	bntseq_t *b = ks_new_bns((int)win.size());
	int64_t off = 0;
	for (size_t i = 0; i < win.size(); ++i) {
		ks_set_ann(&b->anns[i], "w" + std::to_string(i), off, win[i].len);
		off += win[i].len;
	}
	b->l_pac = off;
	uint8_t *p = (uint8_t *)aln_calloc((size_t)(off / 4 + 1), 1);
	int64_t o = 0;
	for (const ks_window_t &w : win) {
		const int64_t base = bns->anns[w.rid].offset + w.beg;
		for (int32_t j = 0; j < w.len; ++j, ++o) {
			const int64_t k = base + j;
			// _set_pac ORs into zeroed bytes, so each base is set once: the ALT in place of REF
			_set_pac(p, o, j == w.site_pos - w.beg ? w.alt : (pac[k >> 2] >> ((~k & 3) << 1)) & 3);
		}
	}
	// The windows keep the reference's ambiguity holes (N runs), translated to window coordinates.
	// Their packed bases are arbitrary fill: without the holes a window row would score that fill
	// as matches, cheaper than the same row on the reference, and the window search's cuts would
	// then start from a best the main search can never confirm.
	std::vector<bntamb1_t> holes;
	if (bns->n_holes > 0) {
		const bntamb1_t *h0 = bns->ambs, *h1 = bns->ambs + bns->n_holes; // sorted by offset
		int64_t wo = 0;
		for (size_t i = 0; i < win.size(); ++i) {
			const ks_window_t &w = win[i];
			const int64_t beg = bns->anns[w.rid].offset + w.beg, end = beg + w.len;
			// the first hole that ends after beg
			const bntamb1_t *h = std::lower_bound(h0, h1, beg, [](const bntamb1_t &a, int64_t x) {
				return a.offset + a.len <= x;
			});
			int32_t n = 0;
			for (; h != h1 && h->offset < end; ++h) {
				const int64_t a = std::max(beg, h->offset), z = std::min(end, h->offset + (int64_t)h->len);
				if (z <= a) continue;
				holes.push_back({wo + (a - beg), (int32_t)(z - a), h->amb});
				++n;
			}
			b->anns[i].n_ambs = n;
			wo += w.len;
		}
	}
	if (!holes.empty()) {
		b->n_holes = (int32_t)holes.size();
		b->ambs = (bntamb1_t *)aln_calloc(holes.size(), sizeof(bntamb1_t));
		std::copy(holes.begin(), holes.end(), b->ambs);
	}
	bns_build_pos2rid(b);
	*out_bns = b;
	*out_pac = p;
}

void ks_windows_meth_convert(const bntseq_t *orig_bns, const uint8_t *orig_pac, bntseq_t **conv_bns,
							 uint8_t **conv_pac)
{
	bntseq_t *b = ks_new_bns(2 * orig_bns->n_seqs);
	b->l_pac = 2 * orig_bns->l_pac;
	uint8_t *p = (uint8_t *)aln_calloc((size_t)(b->l_pac / 4 + 1), 1);
	for (int i = 0; i < orig_bns->n_seqs; ++i) {
		const bntann1_t &a = orig_bns->anns[i];
		for (int h = 0; h < 2; ++h) { // r: G->A, then f: C->T (index --meth's order)
			const int64_t off = 2 * a.offset + h * (int64_t)a.len;
			ks_set_ann(&b->anns[2 * i + h], std::string(h ? "f" : "r") + a.name, off, a.len);
			const uint8_t from = h ? 1 : 2, to = h ? 3 : 0;
			for (int32_t j = 0; j < a.len; ++j) {
				const uint8_t c = _get_pac(orig_pac, a.offset + j);
				_set_pac(p, off + j, c == from ? to : c);
			}
		}
	}
	bns_build_pos2rid(b);
	*conv_bns = b;
	*conv_pac = p;
}

namespace {

/* --known-snvs-cache: the ALT-window FM-index on disk. The FM arrays are a function of the text
 * libsais indexes (the packed windows, converted under --meth) and the SA sample rate only, so the
 * file is keyed on a 128-bit hash of exactly that text, its length, the sample rate, the mode and
 * the array layout. A file whose key, sizes or payload checksum do not match is not used: the
 * index is built in memory, as without the option, and the file rewritten. */
constexpr char KS_CACHE_MAGIC[8] = {'B', 'M', '3', 'K', 'S', 'W', 'I', 'X'};
constexpr uint64_t KS_CACHE_VERSION = 1;
constexpr uint64_t KS_CACHE_ENDIAN = 0x0102030405060708ULL;

struct ks_cache_hdr_t {
	char magic[8];
	uint64_t version, endian;
	uint64_t key[2];          // ks_text_key
	int64_t l_text, sa_compx, meth, sizeof_cp_occ;
	int64_t ref_seq_len, count[5], sentinel_index;
	int64_t cp_occ_count, sa_sample_count;
};

inline uint64_t rotl64(uint64_t x, int r) { return (x << r) | (x >> (64 - r)); }

/* A 64-bit hash of n bytes (four multiply-rotate lanes over 8-byte words, then the tail and the
 * length): a checksum against stale or damaged files, not a cryptographic one. */
uint64_t ks_hash(const void *data, size_t n, uint64_t seed)
{
	const uint64_t P1 = 0x9e3779b185ebca87ULL, P2 = 0xc2b2ae3d27d4eb4fULL, P3 = 0x165667b19e3779f9ULL;
	const uint8_t *b = (const uint8_t *)data;
	uint64_t v[4] = {seed + P1 + P2, seed + P2, seed, seed - P1};
	size_t i = 0;
	for (; i + 32 <= n; i += 32)
		for (int k = 0; k < 4; ++k) {
			uint64_t w;
			memcpy(&w, b + i + 8 * k, 8);
			v[k] = rotl64(v[k] + w * P2, 31) * P1;
		}
	uint64_t h = rotl64(v[0], 1) + rotl64(v[1], 7) + rotl64(v[2], 12) + rotl64(v[3], 18);
	for (int k = 0; k < 4; ++k) h = (h ^ (rotl64(v[k] * P2, 31) * P1)) * P1 + P3;
	for (; i < n; ++i) h = rotl64(h ^ (b[i] * P3), 11) * P1;
	h ^= (uint64_t)n;
	h = (h ^ (h >> 30)) * 0xbf58476d1ce4e5b9ULL;
	h = (h ^ (h >> 27)) * 0x94d049bb133111ebULL;
	return h ^ (h >> 31);
}

void ks_text_key(const uint8_t *text, int64_t l_text, int sa_compx, int meth, uint64_t key[2])
{
	const size_t n = (size_t)(l_text + 3) / 4; // the packed bases; later bits are zero
	const int64_t param[4] = {l_text, sa_compx, meth, (int64_t)sizeof(CP_OCC)};
	for (int h = 0; h < 2; ++h) {
		uint64_t k = ks_hash(text, n, h ? 0x6b73636163686532ULL : 0x6b73636163686531ULL);
		key[h] = ks_hash(param, sizeof param, k);
	}
}

/* Reads the cached arrays for the text with this key into *out (as MemSink allocates them).
 * Returns 0, -1 for a cache miss, or -2 when the path must not be rewritten; failures set *why. */
int ks_cache_load(const char *path, const uint64_t key[2], int64_t l_text, int sa_compx, int meth, FmArrays *out,
				  std::string *why)
{
	const int fd = open(path, O_RDONLY | O_NONBLOCK);
	if (fd < 0) {
		const int err = errno;
		*why = err == ENOENT ? "no file yet" : strerror(err);
		return err == ENOENT ? -1 : -2;
	}
	struct stat st;
	if (fstat(fd, &st) != 0 || !S_ISREG(st.st_mode)) {
		*why = "not a readable regular file";
		close(fd);
		return -2;
	}
	FILE *fp = fdopen(fd, "rb");
	if (!fp) { *why = strerror(errno); close(fd); return -2; }
	struct closer { FILE *f; ~closer() { fclose(f); } } c{fp};
	ks_cache_hdr_t h;
	if (fread(&h, sizeof h, 1, fp) != 1) { *why = "short header"; return -1; }
	if (memcmp(h.magic, KS_CACHE_MAGIC, 8) != 0) { *why = "not an ALT-window index"; return -1; }
	if (h.version != KS_CACHE_VERSION || h.endian != KS_CACHE_ENDIAN || h.sizeof_cp_occ != (int64_t)sizeof(CP_OCC)) {
		*why = "another format version";
		return -1;
	}
	if (h.key[0] != key[0] || h.key[1] != key[1] || h.l_text != l_text || h.sa_compx != sa_compx || h.meth != meth) {
		*why = "stale (built for other sites, reference or options)";
		return -1;
	}
	if (h.ref_seq_len != 2 * l_text + 1 || h.sa_compx < 0 || h.sa_compx > CP_SHIFT) { *why = "bad header"; return -1; }
	const Bwt2bitLayout L = fmi_bwt2bit_layout(h.ref_seq_len, (int)h.sa_compx);
	if (h.cp_occ_count != L.cp_occ_count || h.sa_sample_count != L.sa_sample_count) { *why = "bad header"; return -1; }
	const size_t n_cp = (size_t)L.cp_occ_count * sizeof(CP_OCC), n_ms = (size_t)L.sa_sample_count,
				 n_ls = (size_t)L.sa_sample_count * sizeof(uint32_t);
	if (fstat(fileno(fp), &st) != 0 || (uint64_t)st.st_size != sizeof h + n_cp + n_ms + n_ls + sizeof(uint64_t)) {
		*why = "truncated or oversized";
		return -1;
	}
	FmArrays a;
	a.cp_occ = (CP_OCC *)_mm_malloc(n_cp, 64);
	a.sa_ms_byte = (int8_t *)_mm_malloc(n_ms, 64);
	a.sa_ls_word = (uint32_t *)_mm_malloc(n_ls, 64);
	if (!a.cp_occ || !a.sa_ms_byte || !a.sa_ls_word) { *why = "out of memory"; return -1; }
	uint64_t sum = 0;
	if (fread(a.cp_occ, 1, n_cp, fp) != n_cp || fread(a.sa_ms_byte, 1, n_ms, fp) != n_ms ||
		fread(a.sa_ls_word, 1, n_ls, fp) != n_ls || fread(&sum, sizeof sum, 1, fp) != 1) {
		*why = "read error";
		return -1;
	}
	uint64_t want = ks_hash(&h, sizeof h, 0);
	want = ks_hash(a.cp_occ, n_cp, want);
	want = ks_hash(a.sa_ms_byte, n_ms, want);
	want = ks_hash(a.sa_ls_word, n_ls, want);
	if (sum != want) { *why = "checksum mismatch (damaged)"; return -1; }
	// counts and the sentinel are range-checked again by adopt_built
	out->ref_seq_len = h.ref_seq_len;
	for (int i = 0; i < 5; ++i) out->count[i] = h.count[i];
	out->sentinel_index = h.sentinel_index;
	out->sa_compx = (int)h.sa_compx;
	std::swap(out->cp_occ, a.cp_occ);
	std::swap(out->sa_ms_byte, a.sa_ms_byte);
	std::swap(out->sa_ls_word, a.sa_ls_word);
	return 0;
}

/* Writes the arrays to a unique temporary beside path and renames it over path, so a reader sees
 * the old file or the new one. Returns 0, or -1 with the reason in *why (nothing left behind). */
int ks_cache_save(const char *path, const uint64_t key[2], int64_t l_text, int meth, const FmArrays &a,
				  std::string *why)
{
	ks_cache_hdr_t h;
	memset(&h, 0, sizeof h);
	memcpy(h.magic, KS_CACHE_MAGIC, 8);
	h.version = KS_CACHE_VERSION;
	h.endian = KS_CACHE_ENDIAN;
	h.key[0] = key[0];
	h.key[1] = key[1];
	h.l_text = l_text;
	h.sa_compx = a.sa_compx;
	h.meth = meth;
	h.sizeof_cp_occ = (int64_t)sizeof(CP_OCC);
	h.ref_seq_len = a.ref_seq_len;
	for (int i = 0; i < 5; ++i) h.count[i] = a.count[i];
	h.sentinel_index = a.sentinel_index;
	const Bwt2bitLayout L = fmi_bwt2bit_layout(a.ref_seq_len, a.sa_compx);
	h.cp_occ_count = L.cp_occ_count;
	h.sa_sample_count = L.sa_sample_count;
	const size_t n_cp = (size_t)L.cp_occ_count * sizeof(CP_OCC), n_ms = (size_t)L.sa_sample_count,
				 n_ls = (size_t)L.sa_sample_count * sizeof(uint32_t);
	uint64_t sum = ks_hash(&h, sizeof h, 0);
	sum = ks_hash(a.cp_occ, n_cp, sum);
	sum = ks_hash(a.sa_ms_byte, n_ms, sum);
	sum = ks_hash(a.sa_ls_word, n_ls, sum);
	std::string tmp = std::string(path) + ".XXXXXX";
	std::vector<char> tmpl(tmp.begin(), tmp.end());
	tmpl.push_back('\0');
	const int fd = mkstemp(tmpl.data());
	if (fd < 0) { *why = strerror(errno); return -1; }
	tmp = tmpl.data();
	FILE *fp = fdopen(fd, "wb");
	if (!fp) { *why = strerror(errno); close(fd); unlink(tmp.c_str()); return -1; }
	// mkstemp makes the file 0600: give it the mode of the file it replaces, else a new file's mode
	struct stat old;
	mode_t mode;
	if (stat(path, &old) == 0) mode = old.st_mode & 07777;
	else {
		const mode_t mask = umask(0);
		umask(mask);
		mode = 0666 & ~mask;
	}
	bool ok = fchmod(fd, mode) == 0 && fwrite(&h, sizeof h, 1, fp) == 1 && fwrite(a.cp_occ, 1, n_cp, fp) == n_cp &&
			  fwrite(a.sa_ms_byte, 1, n_ms, fp) == n_ms && fwrite(a.sa_ls_word, 1, n_ls, fp) == n_ls &&
			  fwrite(&sum, sizeof sum, 1, fp) == 1 && fflush(fp) == 0;
	if (!ok) *why = strerror(errno);
	if (fclose(fp) != 0 && ok) { ok = false; *why = strerror(errno); }
	if (ok && rename(tmp.c_str(), path) != 0) { ok = false; *why = strerror(errno); }
	if (!ok) { unlink(tmp.c_str()); return -1; }
	return 0;
}

/* The FM arrays of text (l_text bases, forward only): from cache_path when it holds them, else
 * built in memory (then written to cache_path, when given). *loaded: whether they came from the
 * file. Returns libsais_build_fm_memory's status (0 on success). */
int ks_index_text(const uint8_t *text, int64_t l_text, const LibsaisBuildOpts &bo, int meth, const char *cache_path,
				  FmArrays *out, bool *loaded)
{
	*loaded = false;
	uint64_t key[2] = {0, 0};
	if (cache_path) {
		ks_text_key(text, l_text, bo.sa_compx, meth, key);
		std::string why;
		const int cached = ks_cache_load(cache_path, key, l_text, bo.sa_compx, meth, out, &why);
		if (cached == 0) {
			*loaded = true;
			return 0;
		}
		if (bwa_verbose >= 3)
			fprintf(stderr, "[M::aln] --known-snvs-cache: not using '%s' (%s); building the ALT-window index\n",
					cache_path, why.c_str());
		if (cached == -2) cache_path = nullptr; // Ignore special or unreadable files, without replacing them.
	}
	const int rc = libsais_build_fm_memory(text, l_text, bo, out);
	if (rc != 0 || !cache_path) return rc;
	std::string why;
	if (ks_cache_save(cache_path, key, l_text, meth, *out, &why) != 0) {
		if (bwa_verbose >= 2)
			fprintf(stderr, "[W::aln] --known-snvs-cache: cannot write '%s' (%s); the index is used from memory\n",
					cache_path, why.c_str());
	} else if (bwa_verbose >= 3)
		fprintf(stderr, "[M::aln] --known-snvs-cache: wrote the ALT-window index to '%s'\n", cache_path);
	return 0;
}

} // namespace

int ks_windows_build(const ks_sites_t *sites, const bntseq_t *bns, const uint8_t *pac, int n_threads,
					 int meth, ks_windows_t *out, const char *cache_path)
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
	bntseq_t *wbns = nullptr;
	uint8_t *wpac = nullptr;
	ks_windows_ref(out->win, bns, pac, &wbns, &wpac);
	LibsaisBuildOpts bo;
	bo.num_threads = std::max(1, n_threads);
	bo.sa_compx = 3; // bwa_idx_build's default, which the file-built window index used
	FmArrays arrays;
	if (meth) { // the converted index seeds the search; the original windows are the reference
		bntseq_t *cbns = nullptr;
		uint8_t *cpac = nullptr;
		ks_windows_meth_convert(wbns, wpac, &cbns, &cpac);
		bool loaded = false;
		const int rc = ks_index_text(cpac, cbns->l_pac, bo, 1, cache_path, &arrays, &loaded);
		out->from_cache = loaded;
		free(cpac);
		if (rc != 0) {
			fprintf(stderr, "[E::aln] --known-snvs: cannot index the ALT windows (status %d)\n", rc);
			bns_destroy(cbns);
			bns_destroy(wbns);
			free(wpac);
			ks_windows_free(out);
			return -1;
		}
		out->fmi = new FMI_search("(known-snvs windows)");
		out->fmi->adopt_built(&arrays, cbns, nullptr, n_threads, /*build_seed_jump=*/false); // as load_index(false, ...) leaves it; aln never reads the seed jump table
		out->bwt = bwt_from_fmi(out->fmi);
		out->mref = new aln_meth_ref_t();
		if (aln_meth_ref_adopt(out->mref, wbns, wpac, out->fmi->idx->bns) != 0) { // frees wbns, wpac
			ks_windows_free(out);
			return -1;
		}
		out->mref->ks_window = 1;
		out->mref->ks_win = out->win.data();
		out->bns = out->mref->bns;
	} else {
		bool loaded = false;
		const int rc = ks_index_text(wpac, wbns->l_pac, bo, 0, cache_path, &arrays, &loaded);
		out->from_cache = loaded;
		if (rc != 0) {
			fprintf(stderr, "[E::aln] --known-snvs: cannot index the ALT windows (status %d)\n", rc);
			bns_destroy(wbns);
			free(wpac);
			ks_windows_free(out);
			return -1;
		}
		out->fmi = new FMI_search("(known-snvs windows)");
		out->fmi->adopt_built(&arrays, wbns, wpac, n_threads, /*build_seed_jump=*/false); // aln never reads it
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
	*w = ks_windows_t();
}

namespace {

constexpr int KS_MAX_LOCI = 256; // hits expanded per read; above it the stock result stands

// Differences as the gap search counts them against max_diff (bwtgap.cpp): gap extensions only in
// the default mode, not once -e has cleared BWA_MODE_GAPE.
int gap_diffs(const gap_opt_t *opt, int n_mm, int n_gapo, int n_gape)
{
	return n_mm + n_gapo + ((opt->mode & BWA_MODE_GAPE) ? n_gape : 0);
}

struct ks_hit_t {
	int64_t pos;       // forward-strand start in the main index's packed coordinates
	int strand;
	int score, n_mm, n_gapo, n_gape, ref_shift;
	int64_t alt_site;  // the site allele (aln_site_key) this hit's cost assumes, or -1
	int diffs(const gap_opt_t *opt) const { return gap_diffs(opt, n_mm, n_gapo, n_gape); }
};

// A main or window hit's differences, as ks_hit_t::diffs.
int aln_diffs(const gap_opt_t *opt, const bwt_aln1_t &a) { return gap_diffs(opt, a.n_mm, a.n_gapo, a.n_gape); }

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

/* The window hits (win_aln/win_n) that cover their window's site, in main-index coordinates, into
 * *out. False when they hold more than KS_MAX_LOCI rows. */
bool ks_window_hits(const ks_windows_t *w, const bntseq_t *bns, int len, const bwt_aln1_t *win_aln, int win_n,
					std::vector<ks_hit_t> *out)
{
	std::vector<ks_hit_t> &win_hits = *out;
	win_hits.clear();
	int64_t n_rows = 0;
	for (int k = 0; k < win_n; ++k) n_rows += (int64_t)(win_aln[k].l - win_aln[k].k + 1);
	if (n_rows > KS_MAX_LOCI) return false;
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
			win_hits.push_back(hit_of(q, base + beg, strand, aln_site_key(base + win.site_pos, win.alt)));
		}
	}
	return true;
}

/* The merged set's best score and its cap on differences (aln's cuts: score <= best + mismatch
 * penalty, differences <= min(max_diff, best's + 1)), from the main hits and the window hits that
 * cover their site; differences as stock counts them (gap_diffs). *best is INT32_MAX when there is
 * neither. */
void ks_final_cuts(const std::vector<ks_hit_t> &win_hits, const bwt_aln1_t *aln, int n_aln, int max_diff,
				   const gap_opt_t *opt, int *best, int *cap)
{
	*best = INT32_MAX;
	for (const ks_hit_t &h : win_hits) *best = std::min(*best, h.score);
	for (int k = 0; k < n_aln; ++k) *best = std::min(*best, (int)aln[k].score);
	int best_diffs = INT32_MAX;
	for (const ks_hit_t &h : win_hits)
		if (h.score == *best) best_diffs = std::min(best_diffs, h.diffs(opt));
	for (int k = 0; k < n_aln; ++k)
		if ((int)aln[k].score == *best) best_diffs = std::min(best_diffs, aln_diffs(opt, aln[k]));
	*cap = best_diffs == INT32_MAX ? max_diff : std::max(std::min(max_diff, best_diffs + 1), best_diffs);
}

} // namespace

bool ks_cuts_too_tight(const gap_opt_t *opt, int max_diff, const bwt_aln1_t *a, int n, bool top2b, int best, int cap,
					   bool window)
{
	if (n <= 0) return false; // no hit: it ran to -n
	int own = INT32_MAX;
	for (int k = 0; k < n; ++k) own = std::min(own, (int)a[k].score);
	// A window search's -R counts a locus once per window holding it, so its stop can fire on a
	// stratum the merge counts as fewer loci than -R, and the merge then keeps the hits past it:
	// any such stop counts, also when the merged best is its own.
	const bool below = window || best < own;
	// -N takes no top-2 ceiling and no score cut, and the merge keeps every hit within -n; stock's
	// -R stop (which -R turns on under -N too) still drops every hit past the search's own best
	// stratum, which matters once the merged best is below its own.
	if (opt->mode & BWA_MODE_NONSTOP) return top2b && below;
	// A window search's own best can lie below the merged best: the window index holds the windows
	// back to back, and a row that runs from one window into the next (in no window alone, so the
	// merge drops it) can be its first hit. Its score cut then dropped hits the merge keeps.
	if (window && own < best) return true;
	const int first = aln_diffs(opt, a[0]);
	const int ceiling = first + 1 > max_diff ? max_diff : first + 1;
	if (cap > ceiling) return true;
	// With the score cut, the merged set keeps only hits within best + M: a stop past it lost none.
	return top2b && (int64_t)best + opt->s_mm > own && below;
}

namespace {

/* One read's stock search on bwt (its bases seq, as the search was given them) again, with the
 * cuts of the merged set fixed from the start: score cut best + M (none for best INT32_MAX),
 * difference ceiling cap, and no -R stop, so it finds every hit within them. Returns its hits (the
 * caller frees them) and sets *n, or returns null with *n = -1 when the -m bound ends it. */
bwt_aln1_t *ks_research(const bwt_t *bwt, const gap_opt_t *opt, const ubyte_t *seq, int len, int best, int cap, int *n)
{
	gap_job_t j = gap_job_t();
	j.opt = *opt;
	if (opt->fnr > 0.0) j.opt.max_diff = bwa_cal_maxdiff(len, BWA_AVG_ERR, opt->fnr);
	j.opt.seed_len = opt->seed_len < len ? opt->seed_len : 0x7fffffff;
	j.opt.mode |= BWA_MODE_NONSTOP; // its hits do not move the cuts set here
	j.opt.max_top2 = INT_MAX;
	if (cap < j.opt.max_diff) j.opt.max_diff = cap;
	// a gap open is a difference, so this cap on them is the batch's (aln_batch_max_gapo) in effect
	if (j.opt.max_diff < j.opt.max_gapo) j.opt.max_gapo = j.opt.max_diff;
	std::vector<bwt_width_t> width((size_t)len + 1), seed_width;
	bwt_cal_width(bwt, len, seq, width.data());
	if (len > opt->seed_len) {
		seed_width.resize((size_t)opt->seed_len + 1);
		bwt_cal_width(bwt, opt->seed_len, seq + (len - opt->seed_len), seed_width.data());
	}
	std::vector<ubyte_t> cseq((size_t)len + 1);
	for (int k = 0; k < len; ++k) cseq[k] = seq[k] > 3 ? 4 : 3 - seq[k]; // as the search complements it
	gap_pool_t pool;
	gap_init_pool(&pool, 0);
	j.stack = gap_init_cstack(j.opt.max_diff, j.opt.max_gapo, j.opt.max_gape, &j.opt, &pool);
	j.len = len;
	j.seq = cseq.data();
	j.width = width.data();
	j.seed_width = len > opt->seed_len ? seed_width.data() : nullptr;
	int done = gap_job_start(&j, bwt);
	j.score_cut = gap_score_cut(best, opt->s_mm); // saturates (no cut) for best INT32_MAX
	while (!done) done = gap_job_run(&j, bwt);
	const bool bound = j.max_entries > j.opt.max_entries;
	gap_destroy_cstack(j.stack);
	gap_destroy_pool(&pool);
	if (bound) {
		free(j.aln);
		*n = -1;
		return nullptr;
	}
	*n = j.n_aln;
	return j.aln;
}

/* Merges a read's main hits (aln/n_aln, its main search; top2b: it ended on its -R stop) with its
 * window hits (win_aln/win_n; win_top2b likewise). seq: the read's bases as both searches were
 * given them. Each search took stock's cuts from its own first hit; when the merged set's best is
 * lower (or its ceiling looser) than a search's own, or above a window search's own (a row across
 * two windows, which the merge drops), or a window search stopped on -R, those cuts
 * can have left out hits the merged set keeps, and that search is run again with the merged set's
 * cuts (ks_cuts_too_tight, ks_research). Returns KS_UNTOUCHED (no window hit covers its site),
 * KS_SAME (they change no hit), KS_INCOMPLETE (too many loci, too long for the window, or a search
 * run again stopped at -m), or KS_CREDITED with `m` filled. */
ks_outcome_t ks_merge(const ks_windows_t *w, const bntseq_t *bns, const bwt_t *bwt, const gap_opt_t *opt,
					  int len, const ubyte_t *seq, const bwt_aln1_t *aln, int n_aln, bool top2b,
					  const bwt_aln1_t *win_aln, int win_n, bool win_top2b, ks_merged_t *m)
{
	if (win_n <= 0) return KS_UNTOUCHED; // a search with no hit took no cuts
	std::vector<ks_hit_t> win_hits;
	if (!ks_window_hits(w, bns, len, win_aln, win_n, &win_hits)) return KS_INCOMPLETE;
	const int max_diff = read_max_diff(opt, len);
	const bool nonstop = opt->mode & BWA_MODE_NONSTOP;
	int best, cap;
	ks_final_cuts(win_hits, aln, n_aln, max_diff, opt, &best, &cap);
	// A read too long for a window with a covering hit keeps its stock result whatever is run again.
	if (!win_hits.empty() && !ks_fits_window(len, opt)) return KS_INCOMPLETE;
	// The window search first: its cuts may have left out every hit that covers a site. Its second
	// search takes the merged cuts (under -N, every hit within the read's -n); when a dropped row put
	// its own best below the merged one, a covering hit below the merged best may have more
	// differences than the merged cap allows, so it takes the read's -n.
	std::unique_ptr<bwt_aln1_t, decltype(&free)> win_again(nullptr, &free), main_again(nullptr, &free);
	if (ks_cuts_too_tight(opt, max_diff, win_aln, win_n, win_top2b, best, cap, true)) {
		int own = INT32_MAX;
		for (int k = 0; k < win_n; ++k) own = std::min(own, (int)win_aln[k].score);
		int n;
		win_again.reset(ks_research(w->bwt, opt, seq, len, nonstop ? INT32_MAX : best,
									nonstop || own < best ? max_diff : cap, &n));
		if (n < 0) return KS_INCOMPLETE;
		if (!ks_window_hits(w, bns, len, win_again.get(), n, &win_hits)) return KS_INCOMPLETE;
		// it can find a covering hit below the merged best, which lowers the cuts
		ks_final_cuts(win_hits, aln, n_aln, max_diff, opt, &best, &cap);
	}
	if (win_hits.empty()) return KS_UNTOUCHED;
	if (!ks_fits_window(len, opt)) return KS_INCOMPLETE;
	if (ks_cuts_too_tight(opt, max_diff, aln, n_aln, top2b, best, cap, false)) {
		int n;
		main_again.reset(ks_research(bwt, opt, seq, len, nonstop ? INT32_MAX : best, nonstop ? max_diff : cap, &n));
		if (n < 0) return KS_INCOMPLETE;
		aln = main_again.get();
		n_aln = n;
		// A search run again finds every hit within the cuts the merged set had, so these are the same.
		ks_final_cuts(win_hits, aln, n_aln, max_diff, opt, &best, &cap);
	}

	// Only the main search's intervals that pass the cuts are expanded: crediting often drops
	// hundreds of REF copies that sat one difference behind the stock best.
	auto passes = [&](int score, int diffs) {
		return nonstop ? diffs <= max_diff : score <= (int64_t)best + opt->s_mm && diffs <= cap;
	};
	std::vector<ks_hit_t> hits;
	int64_t n_rows = 0;
	for (int k = 0; k < n_aln; ++k)
		if (passes((int)aln[k].score, aln_diffs(opt, aln[k])))
			n_rows += (int64_t)(aln[k].l - aln[k].k + 1);
	if (n_rows > KS_MAX_LOCI) return KS_INCOMPLETE;
	for (int k = 0; k < n_aln; ++k) {
		const bwt_aln1_t &q = aln[k];
		if (!passes((int)q.score, aln_diffs(opt, q))) continue;
		const int span = len + (int)q.n_del - (int)q.n_ins;
		for (bwtint_t r = q.k; r <= q.l; ++r) {
			int strand;
			const bwtint_t p = bwa_sa2pos(bns, bwt, r, span, &strand);
			if (p == (bwtint_t)-1) continue;
			// A row whose span joins two contigs is not a placement (as the twin checks reject it):
			// it neither competes with nor ties a site-covering hit, nor stays the primary.
			if (!ks_span_in_contig(bns, (int64_t)p, span, (int64_t)p)) continue;
			const ks_hit_t h = hit_of(q, (int64_t)p, strand, -1);
			auto it = std::find_if(hits.begin(), hits.end(), [&](const ks_hit_t &x) {
				return x.strand == h.strand && x.pos == h.pos;
			});
			if (it == hits.end()) hits.push_back(h);
			else if (h.score < it->score) *it = h;
		}
	}
	// Merge: per (strand, start) the lower score wins.
	bool changed = false;
	for (const ks_hit_t &h : win_hits) {
		if (!passes(h.score, h.diffs(opt))) continue;
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
	if (hits.size() > KS_MAX_LOCI) return KS_INCOMPLETE;
	std::stable_sort(hits.begin(), hits.end(), [](const ks_hit_t &a, const ks_hit_t &b) {
		return a.score != b.score ? a.score < b.score : a.pos != b.pos ? a.pos < b.pos : a.strand < b.strand;
	});
	int c1 = 0;
	while (c1 < (int)hits.size() && hits[c1].score == best) ++c1;
	// Main-index rows crossing the forward/reverse junction or joining two contigs may have
	// supplied `best` but no usable position. Keep the stock outcome rather than credit an
	// unrepresented best stratum.
	if (c1 == 0) return KS_INCOMPLETE;
	if (!(opt->mode & BWA_MODE_NONSTOP) && c1 > opt->max_top2) hits.resize(c1);
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
						 int max_xa_hits, bwa_seq_t *s, const bwt_aln1_t *win_aln, int win_n, bool win_top2b,
						 int64_t *alt_site, std::vector<pe_hit_t> *hits, bool stopped)
{
	*alt_site = -1;
	if (hits) hits->clear();
	ks_merged_t m;
	const ks_outcome_t r = ks_merge(w, bns, bwt, opt, s->len, s->seq, s->aln, s->n_aln, s->top2b, win_aln, win_n,
										win_top2b, &m);
	if (r != KS_CREDITED) return r;
	if (stopped) return KS_INCOMPLETE; // a search that stopped at -m may lack competitors
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
	if (hits)
		for (const ks_hit_t &h : kept) {
			hits->push_back({h.pos, h.strand, h.score, h.n_mm, h.n_gapo, h.n_gape, h.ref_shift, (uint64_t)(h.alt_site + 1)});
		}

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
			x.alt_site_plus1 = (uint64_t)(kept[i].alt_site + 1);
		}
	}
	return KS_CREDITED;
}

bool ks_span_in_contig(const bntseq_t *bns, int64_t pos, int64_t span, int64_t contig_pos)
{
	if (contig_pos < 0 || contig_pos >= bns->l_pac || pos < 0 || span <= 0) return false;
	const bntann1_t &ann = bns->anns[bns_pos2rid(bns, contig_pos)];
	const int64_t end = ann.offset + ann.len;
	return pos >= ann.offset && pos < end && span <= end - pos;
}

bool ks_twin_same_contig(const bntseq_t *bns, int64_t twin_pos, int64_t span, int64_t original_pos)
{
	return ks_span_in_contig(bns, twin_pos, span, original_pos) && std::llabs(twin_pos - original_pos) <= KS_TWIN_SLOP;
}

namespace {

/* A placed read's reference span: its CIGAR's, or without one the search's len + ref_shift. This
 * differs from pos_end, which takes len alone without a CIGAR: a gapped search hit that was never
 * given a CIGAR still spans its deletions minus insertions (ref_shift), and an ungapped one has
 * ref_shift 0. */
int64_t ks_read_ref_span(const bwa_seq_t *s)
{
	return s->cigar ? bwa_cigar_ref_span(s->cigar, s->n_cigar) : (int64_t)s->len + s->ref_shift;
}

} // namespace

bool ks_twin_read_same_contig(const bntseq_t *bns, const bwa_seq_t *twin, int64_t original_pos)
{
	if (twin->pos > (uint64_t)INT64_MAX) return false;
	return ks_twin_same_contig(bns, (int64_t)twin->pos, ks_read_ref_span(twin), original_pos);
}

bool ks_twin_pair_same_contig(const bntseq_t *bns, const bwa_seq_t *twin, const bwa_seq_t *mate,
							  int64_t original_pos, int64_t twin_span, int64_t mate_span)
{
	// a CIGAR, else the hit pairing moved the end to, else the search's len + ref_shift
	auto span = [](const bwa_seq_t *s, int64_t moved) {
		return !s->cigar && moved >= 0 ? moved : ks_read_ref_span(s);
	};
	if (twin->pos > (uint64_t)INT64_MAX) return false;
	if (!ks_twin_same_contig(bns, (int64_t)twin->pos, span(twin, twin_span), original_pos)) return false;
	if (mate->type == BWA_TYPE_NO_MATCH) return true;
	if (mate->pos > (uint64_t)INT64_MAX) return false;
	return ks_span_in_contig(bns, (int64_t)mate->pos, span(mate, mate_span), (int64_t)mate->pos);
}

int ks_twin_place(const ks_windows_t *w, const bntseq_t *bns, const bwt_t *bwt, const gap_opt_t *opt,
				  const bwa_seq_t *t, const bwt_aln1_t *win_aln, int win_n, bool win_top2b, int64_t pos,
				  int strand, int *mapq)
{
	*mapq = 0;
	ks_merged_t m;
	const ks_outcome_t r = ks_merge(w, bns, bwt, opt, t->len, t->seq, t->aln, t->n_aln, t->top2b, win_aln, win_n,
										win_top2b, &m);
	if (r == KS_INCOMPLETE) return -1;
	if (r == KS_CREDITED) {
		const ks_hit_t &b = m.kept[0];
		*mapq = approx_mapq(m.c1, (int)m.kept.size() - m.c1, b.n_mm, m.max_diff);
		if (m.c1 > 1) return 0; // MAPQ 0: a tie never passes
		return b.strand == strand && ks_twin_same_contig(bns, b.pos, (int64_t)t->len + b.ref_shift, pos);
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
	return p != (bwtint_t)-1 && s == strand
		&& ks_twin_same_contig(bns, (int64_t)p, (int64_t)t->len + q.n_del - q.n_ins, pos);
}

int ks_site_rel(const ks_sites_t *sites, const bntseq_t *bns, int rid, int64_t pos, int64_t key)
{
	if (key < 0 || !sites || rid < 0 || rid >= (int)sites->by_contig.size()) return 0;
	const std::vector<ks_site_t> &lst = sites->by_contig[rid];
	const int64_t off = bns->anns[rid].offset;
	auto it = std::lower_bound(lst.begin(), lst.end(), pos - off, [](const ks_site_t &x, int64_t p) { return x.pos < p; });
	for (int rel = 1; it != lst.end() && rel <= 255; ++it, ++rel)
		if (aln_site_key(off + it->pos, it->alt) == key) return rel;
	return 0;
}

int64_t ks_site_of_rel(const ks_sites_t *sites, const bntseq_t *bns, int rid, int64_t pos, int rel)
{
	if (rel <= 0 || !sites || rid < 0 || rid >= (int)sites->by_contig.size()) return -1;
	const std::vector<ks_site_t> &lst = sites->by_contig[rid];
	const int64_t off = bns->anns[rid].offset;
	auto it = std::lower_bound(lst.begin(), lst.end(), pos - off, [](const ks_site_t &x, int64_t p) { return x.pos < p; });
	if (lst.end() - it < rel) return -1;
	it += rel - 1;
	return aln_site_key(off + it->pos, it->alt);
}

int ks_site_votes(const ks_sites_t *sites, const bntseq_t *bns, const gap_opt_t *opt, const bwa_seq_t *s,
				  int64_t alt_site, int hyp, std::vector<ks_vote_t> *out, const uint8_t *pac, int cpg, int collapsed)
{
	out->clear();
	if (s->type != BWA_TYPE_UNIQUE && s->type != BWA_TYPE_REPEAT && s->type != BWA_TYPE_MATESW) return -1;
	const int rid = bns_pos2rid(bns, (int64_t)s->pos);
	if (rid < 0) return -1;
	const int64_t cbeg = (int64_t)s->pos - bns->anns[rid].offset;
	const int64_t span = pos_end(s) - (int64_t)s->pos;
	if (span > bns->anns[rid].len - cbeg) return -1;
	const std::vector<ks_site_t> &lst = sites->by_contig[rid];
	auto it = std::lower_bound(lst.begin(), lst.end(), cbeg, [](const ks_site_t &x, int64_t p) { return x.pos < p; });
	if (it == lst.end() || it->pos >= cbeg + span) return -1;
	// The read in reference orientation (s->seq is in read order after refinement).
	std::vector<uint8_t> rs(s->len);
	for (int i = 0; i < s->len; ++i) rs[i] = s->strand ? (s->seq[s->len - 1 - i] > 3 ? 4 : 3 - s->seq[s->len - 1 - i]) : s->seq[i];
	int diffs = gap_diffs(opt, s->n_mm, s->n_gapo, s->n_gape); // as counted against max_diff
	// The base the placement's cost assumes at each site: the credited site allele's ALT there,
	// the REF elsewhere.
	const int64_t one_pos = alt_site >= 0 ? aln_site_key_pos(alt_site) : -1;
	const int one_alt = alt_site >= 0 ? aln_site_key_alt(alt_site) : -1;
    // Plain refinement may normalize an adjacent gap. Count the emitted alignment against
    // the selected allele rather than mixing its CIGAR with the search history's counts.
    if (hyp < 0 && pac) {
        diffs = 0;
        int qi = 0;
        int64_t rp = (int64_t)s->pos;
        const int nc = s->cigar ? s->n_cigar : 1;
        for (int j = 0; j < nc; ++j) {
            const int op = s->cigar ? __cigar_op(s->cigar[j]) : FROM_M;
            const int len = s->cigar ? __cigar_len(s->cigar[j]) : s->len;
            if (op == FROM_M) {
                for (int k = 0; k < len; ++k) {
                    const int ref = rp == one_pos ? one_alt : (int)_get_pac(pac, rp);
                    diffs += rs[qi] > 3 || ref != rs[qi]; ++qi; ++rp;
                }
            } else if (op == FROM_I || op == FROM_D) {
                diffs += 1 + ((opt->mode & BWA_MODE_GAPE) ? len - 1 : 0);
                if (op == FROM_I) qi += len; else rp += len;
            } else if (op == FROM_S) qi += len;
        }
    }
	const int max_diff = read_max_diff(opt, s->clip_len > 0 ? s->clip_len : s->len);
	// Under --meth a converted base (C->T, or G->A) also reads as its source.
	const int source = hyp == 1 ? 1 : 2, target = hyp == 1 ? 3 : 0;
	const int64_t c_off = bns->anns[rid].offset, c_len = bns->anns[rid].len;
	// Under --meth-context cpg the aligner converts only at a CpG (meth_label_cpg): the base after
	// a C->T site, or before a G->A one, is the CpG's other half, read with the credited ALT in
	// place; an ambiguity run at either base undoes the label.
	auto converts_at = [&](int64_t cpos) {
		if (hyp < 0) return false;
		if (!cpg || !pac) return true;
		const int64_t npos = hyp == 1 ? cpos + 1 : cpos - 1;
		if (npos < 0 || npos >= c_len) return false;
		const int64_t p0 = c_off + std::min(cpos, npos);
		if (bns_cnt_ambi(bns, p0, 2, nullptr) > 0) return false;
		const int64_t nk = c_off + npos;
		const int nb = one_pos == nk ? one_alt : (int)_get_pac(pac, nk);
		return nb == (hyp == 1 ? 2 : 1);
	};
	bool conv = false; // converts_at the current site
	auto in_set = [&](uint8_t b, uint8_t allele) { return b == allele || (conv && allele == source && b == target); };
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
		conv = converts_at(rpos);
		v.conv = conv;
		v.site = &*it;
		v.rid = rid;
		v.qi = qi;
		const uint8_t b = rs[qi];
		const bool is_ref = in_set(b, it->ref), is_alt = in_set(b, it->alt);
		v.allele = is_ref && is_alt ? 'H' : is_ref ? 'R' : is_alt ? 'A' : 'O';
		const uint8_t assumed = one_pos == bns->anns[rid].offset + rpos ? (uint8_t)one_alt : it->ref;
		const bool scored_match = in_set(b, assumed) || (collapsed && assumed == target && b == source);
		v.away = diffs - (scored_match ? 0 : 1); // the site's own base is not a difference away
		v.max_diff = max_diff;
		if (first < 0 && (v.allele == 'R' || v.allele == 'A')) first = (int)out->size();
		out->push_back(v);
	}
	return first;
}

void ks_twin_seq(const bwa_seq_t *s, const ks_vote_t &v, int hyp, std::vector<uint8_t> *read_order)
{
	// s->seq is in read order after refinement; v.qi indexes the reference-oriented read.
	const int len = s->clip_len > 0 ? s->clip_len : s->len;
	read_order->assign(s->seq, s->seq + len);
	uint8_t other = v.allele == 'R' ? v.site->alt : v.site->ref;
	if (hyp == 1 && v.conv && other == 1) other = 3;      // C reads as T on a C->T strand
	else if (hyp == 0 && v.conv && other == 2) other = 0; // G reads as A on a G->A strand
	const int qi = s->strand ? (int)s->len - 1 - v.qi : v.qi;
	if (qi >= 0 && qi < len) (*read_order)[qi] = s->strand ? 3 - other : other;
}

void ks_format_xv(const bntseq_t *bns, const std::vector<ks_vote_t> &votes, int twin_of, int twin_same,
				  int twin_mapq, std::string *out)
{
	out->clear();
	for (size_t i = 0; i < votes.size(); ++i) {
		const ks_vote_t &v = votes[i];
		char buf[64];
		if (!out->empty()) out->push_back(';');
		for (const unsigned char *p = (const unsigned char *)bns->anns[v.rid].name; *p; ++p) {
			if (*p == '%' || *p == ':' || *p == ',' || *p == ';') {
				out->push_back('%');
				out->push_back("0123456789ABCDEF"[*p >> 4]);
				out->push_back("0123456789ABCDEF"[*p & 15]);
			} else out->push_back((char)*p);
		}
		snprintf(buf, sizeof buf, ":%lld:%c>%c,%c,%d,%d,", (long long)v.site->pos + 1, "ACGT"[v.site->ref],
				 "ACGT"[v.site->alt], v.allele, v.away, v.max_diff);
		out->append(buf);
		if ((int)i == twin_of && twin_same >= 0) snprintf(buf, sizeof buf, "%d,%d", twin_same, twin_mapq);
		else snprintf(buf, sizeof buf, "NA,NA");
		out->append(buf);
	}
}

/* The search path behind an ALT-window hit's display CIGAR c (n operations, reference span span at
 * window offset off): the display soft-clips a terminal insertion and drops a terminal deletion
 * (moving the start past a leading one), so S is restored to I, and the deletions the hit's gap
 * counts hold beyond the display's are restored at its ends. Only the ends' total deletion length
 * is known, but where it sits does not change the path's cost or its display alignment: it goes
 * first before the display (one base at the end when both ends lost one), within the window's
 * len bases. *lead and *trail are the reference bases restored before and after the display. An
 * empty result: no path is recoverable (the hit is refined by the window DP instead). */
static std::vector<uint32_t> ks_hit_path(const uint32_t *c, uint32_t n, const aln_meth_hit_t &h, int64_t off,
										 int64_t len, int *lead, int *trail)
{
	*lead = *trail = 0;
	std::vector<uint32_t> raw(c, c + n);
	int gapo = 0, gape = 0;
	for (uint32_t &op : raw) {
		const int o = op & 15, l = (int)(op >> 4);
		if (o == 4) op = (op & ~15u) | 1; // a soft clip is a terminal insertion
		else if (o != 1 && o != 2) continue;
		++gapo; gape += l - 1;
	}
	const int lost_o = (int)h.n_gapo - gapo, lost_e = (int)h.n_gape - gape;
	if (lost_o < 0 || lost_o > 2 || lost_e < 0 || (lost_o == 0 && lost_e > 0)) return {};
	if (lost_o == 0) return raw;
	const int total = lost_o + lost_e; // the lost deletions' bases
	const int64_t room_after = len - off - h.span;
	if (lost_o == 1) {
		if (off >= total) *lead = total;
		else if (room_after >= total) *trail = total;
		else return {};
	} else if (off >= total - 1 && room_after >= 1) *lead = total - 1, *trail = 1;
	else if (off >= 1 && room_after >= total - 1) *lead = 1, *trail = total - 1;
	else return {};
	if (*lead) raw.insert(raw.begin(), (uint32_t)*lead << 4 | 2);
	if (*trail) raw.push_back((uint32_t)*trail << 4 | 2);
	return raw;
}

void ks_meth_candidates(const ks_windows_t *w, const bntseq_t *bns, const aln_meth_read_t &wres,
						aln_meth_extra_t *out)
{
	out->rows.clear();
	out->horizon = wres.horizon;
	out->flags = wres.flags;
	for (const aln_meth_hit_t &h : wres.hits) {
		if (h.rid < 0 || h.rid >= (int)w->win.size()) continue;
		const ks_window_t &win = w->win[h.rid];
		const int64_t off = (int64_t)h.pos - w->bns->anns[h.rid].offset;
		if (off < 0 || off + h.span > win.len) continue;
		const int64_t beg = win.beg + off;
		if (win.site_pos < beg || win.site_pos >= beg + h.span) continue; // misses its site
		const int64_t rb = bns->anns[win.rid].offset + beg;
		// The candidate is the search path (raw CIGAR and bounds), not the display alignment: the
		// genomic search scores it as it scores its own rows, terminal gaps included.
		int lead, trail;
		std::vector<uint32_t> path = ks_hit_path(aln_meth_hit_cigar(wres, h), h.n_cig, h, off, win.len, &lead, &trail);
		out->rows.push_back({win.rid, (int)h.strand, (int)h.hyp, rb - lead, rb + h.span + trail, std::move(path)});
	}
}

std::atomic<bool> &ks_site_cap_warned()
{
	static std::atomic<bool> warned{false};
	return warned;
}

void ks_rescue_sites(const void *ctx, int64_t beg, int64_t end, std::vector<std::pair<int64_t, uint8_t>> *out)
{
	const ks_rescue_ctx_t *c = (const ks_rescue_ctx_t*)ctx;
	if (end <= beg) return;
	const int rid0 = bns_pos2rid(c->bns, beg);
	if (rid0 < 0) return;
	// The range can run past its first contig (rescue windows are not clipped to one): every
	// contig it reaches, in order.
	for (int rid = rid0; rid < c->bns->n_seqs && c->bns->anns[rid].offset < end; ++rid) {
		const int64_t off = c->bns->anns[rid].offset;
		const std::vector<ks_site_t> &lst = c->sites->by_contig[rid];
		auto it = std::lower_bound(lst.begin(), lst.end(), beg - off, [](const ks_site_t &x, int64_t p) { return x.pos < p; });
		for (; it != lst.end() && off + it->pos < end; ++it) {
			if ((int)out->size() == KS_RESCUE_MAX_SITES) {
				// The later sites stay at REF in this window: say so once per run.
				if (!ks_site_cap_warned().exchange(true))
					fprintf(stderr, "[W::aln] --known-snvs: a window holds more than %d site alleles; the later "
							"ones are scored as the reference there (warned once)\n", KS_RESCUE_MAX_SITES);
				return;
			}
			out->push_back({off + it->pos, it->alt});
		}
	}
}

void ks_aln2seq_nodraw(const bntseq_t *bns, bwa_seq_t *s)
{
	// bwa_aln2seq_core's set_main branch with its draws taken from a per-read stream (erand48)
	uint64_t h = read_seed(bns, s);
	unsigned short x[3] = {(unsigned short)h, (unsigned short)(h >> 16), (unsigned short)(h >> 32)};
	if (s->n_aln == 0) {
		s->type = BWA_TYPE_NO_MATCH;
		s->c1 = s->c2 = 0;
		return;
	}
	const bwt_aln1_t *aln = s->aln;
	int i;
	int64_t cnt = 0;
	const int best = aln[0].score;
	for (i = 0; i < s->n_aln; ++i) {
		const bwt_aln1_t *p = aln + i;
		if (p->score > best) break;
		if (erand48(x) * (double)(p->l - p->k + 1 + cnt) > (double)cnt) {
			s->n_mm = p->n_mm; s->n_gapo = p->n_gapo; s->n_gape = p->n_gape;
			s->ref_shift = (int)p->n_del - (int)p->n_ins;
			s->score = p->score;
			s->sa = p->k + (bwtint_t)((double)(p->l - p->k + 1) * erand48(x));
		}
		cnt += (int64_t)(p->l - p->k + 1);
	}
	s->c1 = cnt;
	for (; i < s->n_aln; ++i) cnt += (int64_t)(aln[i].l - aln[i].k + 1);
	s->c2 = cnt - s->c1;
	s->type = s->c1 > 1 ? BWA_TYPE_REPEAT : BWA_TYPE_UNIQUE;
}

gap_opt_t ks_window_opt(const gap_opt_t &opt)
{
	gap_opt_t window = opt;
	if (opt.mode & BWA_MODE_NONSTOP) window.max_top2 = INT_MAX;
	return window;
}

int ks_meth_search(const aln_opt_t *o, const bwt_t *bwt, int n, bwa_seq_t *seqs, int read_no, int64_t ord,
				   int keep, aln_dup_state_t *dup, aln_dup_state_t *window_dup, aln_meth_read_t *res)
{
	std::vector<aln_meth_extra_t> extra;
	if (o->ks_windows) {
		// The cap: a probe of the main search (aln_meth_probe_batch, its exact stratum and a few
		// verified rows) gives a read with such a row an upper bound U on its final best, and its
		// window search starts from the seed U + M (ks_meth_window_seed), so it stops at U + 2M: a
		// read whose own search stops on -R retains reference rows within B + M that only a window
		// holds, at up to U + 2M there, for pairing and for single-end XA alike. A read with no probe
		// row, a read that does not fit an ALT window, and every read under --meth-context cpg run
		// their window search uncapped. Under -N or a negative penalty the window search takes no
		// score cut, so there is no probe.
		const gap_opt_t &g = o->gap;
		const bool cpg = o->meth->m.meth_context == METH_CTX_CPG;
		std::vector<int64_t> caps;
		if (aln_meth_takes_cuts(&g) && !cpg) {
			caps.assign((size_t)n, INT64_MAX);
			aln_meth_probe_batch(bwt, o->meth_ref, o->meth, n, seqs, &g, g.n_threads, read_no, ord, caps.data());
			for (int i = 0; i < n; ++i) caps[i] = ks_meth_window_seed(caps[i], seqs[i].len, &g, cpg);
		}
		// The ALT-window search runs with the read's own options: as a candidate generator
		// (aln_meth_boundary_gen) it takes only the cuts that hold for any final best.
		std::vector<aln_meth_read_t> windows(n);
		aln_meth_search_batch(o->ks_windows->bwt, o->ks_windows->mref, o->meth, n, seqs, &g, g.n_threads, read_no,
							  ord, 1, o->max_xa_hits, window_dup, windows.data(), nullptr,
							  caps.empty() ? nullptr : caps.data());
		extra.resize(n);
		for (int i = 0; i < n; ++i) ks_meth_candidates(o->ks_windows, o->meth_ref->bns, windows[i], &extra[i]);
	}
	return aln_meth_search_batch(bwt, o->meth_ref, o->meth, n, seqs, &o->gap, o->gap.n_threads, read_no, ord, keep,
								 o->max_xa_hits, dup, res, extra.empty() ? nullptr : extra.data());
}

} // namespace bwa_aln
