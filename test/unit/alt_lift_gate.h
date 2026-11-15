// test/unit/alt_lift_gate.h — shared helpers for the ALT lift tests.
//
// The fixtures (FakeBns, TmpFile, LiftGuard, lift-table comparison) and the
// exclude-BED table gate (restricted loading versus a FLAG-masked sidecar),
// shared by test/unit/test_alt_lift.cpp, which runs the gate on generated
// indexes, and test/integration/test_alt_lift_frozen.cpp, which runs it on a
// real index named by environment variables. Include it from one translation
// unit per test binary: it registers a doctest listener.
#ifndef BWA_TEST_ALT_LIFT_GATE_H
#define BWA_TEST_ALT_LIFT_GATE_H

#include "doctest/doctest.h"
#include "temp_index.h"
#include "temp_path.h"
#include "alt_lift.h"
#include "bntseq.h"
#include "bwamem.h"

#include <algorithm>
#include <atomic>
#include <cerrno>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <functional>
#include <random>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <unistd.h>
#include <vector>
#include <zlib.h>
#include <sys/stat.h>

namespace {
struct Ctg {
    const char *name;
    int len;
    int is_alt;
};

// A minimal bns: names, lengths, offsets and is_alt are all these paths read.
struct FakeBns {
    bntseq_t bns;
    std::vector<bntann1_t> anns;
    explicit FakeBns(const std::vector<Ctg> &ctgs) : anns(ctgs.size()) {
        memset(&bns, 0, sizeof(bns));
        int64_t off = 0;
        for (size_t i = 0; i < ctgs.size(); ++i) {
            memset(&anns[i], 0, sizeof(bntann1_t));
            anns[i].name = const_cast<char *>(ctgs[i].name);
            anns[i].len = ctgs[i].len;
            anns[i].is_alt = ctgs[i].is_alt;
            anns[i].offset = off;
            off += ctgs[i].len;
        }
        bns.n_seqs = (int32_t) ctgs.size();
        bns.anns = anns.data();
        bns.l_pac = off;
    }
};

struct TmpFile {
    std::string path;
    explicit TmpFile(const char *body) {
        char tmpl[4096];
        int fd = bwa_test_mkstemp(tmpl, sizeof(tmpl), "alt_lift_test");
        REQUIRE(fd >= 0);
        path = tmpl;
        FILE *fp = fdopen(fd, "w");
        REQUIRE(fp != nullptr);
        fputs(body, fp);
        fclose(fp);
    }
    ~TmpFile() { unlink(path.c_str()); }
};

struct LiftGuard {
    alt_lift_t *l;
    explicit LiftGuard(alt_lift_t *x) : l(x) {}
    ~LiftGuard() { alt_lift_destroy(l); }
};

// Compare lift data, never allocation addresses or struct padding.
void check_lift_tables(const alt_lift_t *actual, const alt_lift_t *expected) {
    REQUIRE(actual != nullptr);
    REQUIRE(expected != nullptr);
    REQUIRE(actual->n_seqs == expected->n_seqs);
    CHECK(actual->n_lift == expected->n_lift);
    for (int rid = 0; rid < actual->n_seqs; ++rid) {
        CAPTURE(rid);
        const auto &a = actual->a[rid], &b = expected->a[rid];
        REQUIRE(a.trid == b.trid);
        if (a.trid < 0) continue;
        CHECK(a.rev == b.rev);
        CHECK(a.tpos == b.tpos);
        CHECK(a.tlen == b.tlen);
        CHECK(a.qlen == b.qlen);
        REQUIRE(a.n_cigar == b.n_cigar);
        for (int i = 0; i < a.n_cigar; ++i) CHECK(a.cigar[i] == b.cigar[i]);
        for (int i = 0; i <= a.n_cigar; ++i) {
            CHECK(a.qoff[i] == b.qoff[i]);
            CHECK(a.roff[i] == b.roff[i]);
        }
    }
}

void write_gzip(const std::string &path, const std::string &text) {
    gzFile fp = gzopen(path.c_str(), "wb");
    REQUIRE(fp != nullptr);
    CHECK(gzwrite(fp, text.data(), (unsigned) text.size()) == (int) text.size());
    REQUIRE(gzclose(fp) == Z_OK);
}

struct OptGuard {
    mem_opt_t *o;
    OptGuard() : o(mem_opt_init()) {}
    ~OptGuard() { free(o); }
};

// chr1 (primary, 100 kb); alt_fwd and alt_rev: 1500 bases matching chr1:3000-4499
// (0-based), a 3-bp deletion, 497 matching; alt_rev is its reverse complement
// (FLAG 16). alt_clip: 100S200M at chr1:20001 (ALT-only first 100 bases).
// alt_bigdel: 100M500D100M at chr1:30001 (an ALT deletion wider than the band).
// alt_ins: 100M80I100M at chr1:50001 (an ALT insertion shorter than the band).
// alt_wdel: 100M100D100M at chr1:60001 (an ALT deletion of exactly the default
// band, -w 100). alt_bal: 100M5I30M5D100M at chr1:70001 (indels that cancel).
// alt_b35 and alt_b40: 100M{k}I30M{k}D100M at chr1:80001 and 81001, indels that
// cancel 30 bases apart with an excursion of 35 and 40. alt_maxdel and
// alt_hugedel: 100M1022D100M at chr1:90001 and 100M1100D100M at chr1:93001, ALT
// deletions of exactly MEM_ALT_BW_MAX and wider than it. alt_fwd2: a copy of
// alt_fwd, last in the index.
const std::vector<Ctg> kCtgs = {
    {"chr1", 100000, 0}, {"alt_fwd", 1997, 1}, {"alt_rev", 1997, 1},
    {"alt_clip", 300, 1}, {"alt_bigdel", 200, 1}, {"alt_ins", 280, 1},
    {"alt_wdel", 200, 1}, {"alt_bal", 235, 1}, {"alt_b35", 265, 1}, {"alt_b40", 270, 1},
    {"alt_maxdel", 200, 1}, {"alt_hugedel", 200, 1}, {"alt_fwd2", 1997, 1},
};
const char *kAlt =
    "@HD\tVN:1.6\n"
    "alt_fwd\t0\tchr1\t3001\t60\t1500M3D497M\t*\t0\t0\t*\t*\n"
    "alt_rev\t16\tchr1\t3001\t60\t1500M3D497M\t*\t0\t0\t*\t*\n"
    "alt_clip\t0\tchr1\t20001\t60\t100S200M\t*\t0\t0\t*\t*\n"
    "alt_bigdel\t0\tchr1\t30001\t60\t100M500D100M\t*\t0\t0\t*\t*\n"
    "alt_ins\t0\tchr1\t50001\t60\t100M80I100M\t*\t0\t0\t*\t*\n"
    "alt_wdel\t0\tchr1\t60001\t60\t100M100D100M\t*\t0\t0\t*\t*\n"
    "alt_bal\t0\tchr1\t70001\t60\t100M5I30M5D100M\t*\t0\t0\t*\t*\n"
    "alt_b35\t0\tchr1\t80001\t60\t100M35I30M35D100M\t*\t0\t0\t*\t*\n"
    "alt_b40\t0\tchr1\t81001\t60\t100M40I30M40D100M\t*\t0\t0\t*\t*\n"
    "alt_maxdel\t0\tchr1\t90001\t60\t100M1022D100M\t*\t0\t0\t*\t*\n"
    "alt_hugedel\t0\tchr1\t93001\t60\t100M1100D100M\t*\t0\t0\t*\t*\n"
    "alt_fwd2\t0\tchr1\t3001\t60\t1500M3D497M\t*\t0\t0\t*\t*\n";

enum { CHR1 = 0, ALT_FWD = 1, ALT_REV = 2, ALT_CLIP = 3, ALT_BIGDEL = 4, ALT_INS = 5, ALT_WDEL = 6,
       ALT_BAL = 7, ALT_B35 = 8, ALT_B40 = 9, ALT_MAXDEL = 10, ALT_HUGEDEL = 11, ALT_FWD2 = 12 };

// One alignment region: forward or reverse strand of [cs, ce) on contig rid, in
// bwa's doubled coordinate space. `tag` (stored in seedcov) identifies it after
// the rewrite.
mem_alnreg_t reg(const FakeBns &fb, int rid, int64_t cs, int64_t ce, int rev, int score, int tag,
                 int qb = 0, int qe = 150) {
    mem_alnreg_t r;
    memset(&r, 0, sizeof(r));
    const int64_t off = fb.bns.anns[rid].offset, l2 = fb.bns.l_pac << 1;
    r.rid = rid;
    r.rb = rev ? l2 - (off + ce) : off + cs;
    r.re = rev ? l2 - (off + cs) : off + ce;
    r.qb = qb; r.qe = qe;
    r.score = r.truesc = score;
    r.secondary = -1;
    r.is_alt = fb.bns.anns[rid].is_alt;
    r.seedcov = tag;
    return r;
}
// Resolve already-projected candidates without rescue, to isolate pairing/MAPQ.
struct PairOutcome { int paired, flag, q[2], score[2]; int64_t rb[2]; };
PairOutcome project_and_pair(mem_opt_t *opt, const FakeBns &fb,
                             const std::vector<mem_alnreg_t> input[2]) {
    std::vector<mem_alnreg_t> regions[2] = {input[0], input[1]};
    opt->flag |= MEM_F_NO_RESCUE;
    mem_alnreg_v a[2] = {};
    for (int i = 0; i < 2; ++i) {
        regions[i].resize(mem_alt_project(opt, &fb.bns, (int) regions[i].size(), regions[i].data()));
        a[i].a = regions[i].data();
        a[i].n = a[i].m = regions[i].size();
    }
    mem_pestat_t pes[4] = {};
    for (auto &p : pes) p.failed = 1;
    pes[1].failed = 0;
    pes[1].avg = 299; pes[1].std = 30; pes[1].low = 180; pes[1].high = 420;
    bseq1_t reads[2] = {};
    mem_cache *cache = (mem_cache *) calloc(1, sizeof(mem_cache));
    REQUIRE(cache != nullptr);
    int32_t gcnt = 0;
    int n_pri[2], selected[2] = {-1, -1};
    PairOutcome out = {};
    mem_pair_resolve_batch_post(opt, &fb.bns, nullptr, pes, 7, reads, a, nullptr, cache,
                               gcnt, 0, n_pri, selected, out.q, &out.flag, &out.paired);
    for (int i = 0; i < 2; ++i) {
        out.rb[i] = out.paired ? a[i].a[selected[i]].rb : -1;
        out.score[i] = out.paired ? a[i].a[selected[i]].score : -1;
    }
    free(cache);
    return out;
}
} // namespace

namespace {
struct BnsGuard {
    bntseq_t *b;
    explicit BnsGuard(bntseq_t *x) : b(x) {}
    ~BnsGuard() { if (b) bns_destroy(b); }
    BnsGuard(const BnsGuard &) = delete;
    BnsGuard &operator=(const BnsGuard &) = delete;
};

void write_text(const std::string &path, const std::string &text) {
    FILE *fp = fopen(path.c_str(), "wb");
    REQUIRE(fp != nullptr);
    CHECK(fwrite(text.data(), 1, text.size(), fp) == text.size());
    REQUIRE(fclose(fp) == 0);
}

std::vector<std::vector<std::string>> sam_fields(const std::string &text) {
    std::vector<std::vector<std::string>> rows;
    size_t at = 0;
    while (at < text.size()) {
        size_t end = text.find('\n', at);
        if (end == std::string::npos) end = text.size();
        std::vector<std::string> row;
        size_t f = at;
        while (true) {
            const size_t tab = text.find('\t', f);
            if (tab == std::string::npos || tab > end) { row.push_back(text.substr(f, end - f)); break; }
            row.push_back(text.substr(f, tab - f));
            f = tab + 1;
        }
        rows.push_back(row);
        at = end + 1;
    }
    return rows;
}
} // namespace

// ---------------------------------------------------------------------------
// Table gate: restricted loading versus a FLAG-masked sidecar.
//
// The masked sidecar is the source with FLAG 0x4 OR-ed into every record whose
// primary footprint overlaps the BED, as an independent overlap oracle writes
// it. Loading source + BED must give the same lift table, field by field, as
// loading the masked sidecar without a BED, against a bns that each sidecar's
// own names marked. The gate also checks the inputs it relies on to equate
// per-record decisions with table entries: identical QNAME order and non-FLAG
// bytes, unique QNAMEs, and every source record usable when loaded unfiltered.
// ---------------------------------------------------------------------------
namespace {
struct AltGateInputs {
    std::string prefix, source, masked, bed;
};

struct AltGateCounts {
    long n_seqs = 0, alt_contigs = 0, rows = 0, excluded = 0, retained = 0;
    long input_intervals = 0, merged_intervals = 0;
};

struct SidecarRow {
    std::string qname;
    long flag;
    std::string rest;  // every byte after FLAG up to the newline
};

bool file_exists(const std::string &path) {
    struct stat st;
    return stat(path.c_str(), &st) == 0 && S_ISREG(st.st_mode);
}

std::vector<SidecarRow> read_sidecar(const std::string &path) {
    std::ifstream in(path, std::ios::binary);
    REQUIRE(in.good());
    std::vector<SidecarRow> rows;
    std::string line;
    while (std::getline(in, line)) {
        const size_t t1 = line.find('\t');
        const size_t t2 = t1 == std::string::npos ? t1 : line.find('\t', t1 + 1);
        REQUIRE(t2 != std::string::npos);
        REQUIRE(line[0] != '@');
        char *end = nullptr;
        const std::string flag = line.substr(t1 + 1, t2 - t1 - 1);
        const long f = strtol(flag.c_str(), &end, 10);
        REQUIRE((!flag.empty() && *end == '\0'));
        rows.push_back({line.substr(0, t1), f, line.substr(t2)});
    }
    return rows;
}

// bns_restore a private prefix whose .ann/.amb/.pac link to `prefix` and whose
// .alt links to `alt`, so ALT marking comes from that sidecar's names and the
// shared index directory is never written.
bntseq_t *restore_with_alt(const std::string &dir, const char *tag, const std::string &prefix,
                           const std::string &alt) {
    const std::string base = dir + "/" + tag;
    auto link = [](const std::string &target, const std::string &at) {
        char *abs = realpath(target.c_str(), nullptr);
        REQUIRE(abs != nullptr);
        const int rc = symlink(abs, at.c_str());
        free(abs);
        REQUIRE(rc == 0);
    };
    for (const char *ext : {".ann", ".amb", ".pac"}) link(prefix + ext, base + ext);
    link(alt, base + ".alt");
    return bns_restore(base.c_str());
}

// Compare contig dictionaries and ALT bits. Returns the number of differing contigs.
long count_bns_differences(const bntseq_t *a, const bntseq_t *b) {
    long diff = 0;
    for (int rid = 0; rid < a->n_seqs; ++rid) {
        const bntann1_t &x = a->anns[rid], &y = b->anns[rid];
        if (strcmp(x.name, y.name) != 0 || x.len != y.len || x.offset != y.offset || x.is_alt != y.is_alt) {
            if (diff == 0) MESSAGE("first differing contig: rid " << rid << " " << x.name);
            ++diff;
        }
    }
    return diff;
}

// The whole gate. `tamper`, if given, edits the oracle table after loading; the
// generated cases use it to prove a single wrong field is caught.
AltGateCounts run_alt_table_gate(const AltGateInputs &in,
                                 const std::function<void(alt_lift_t *)> &tamper = nullptr) {
    for (const char *ext : {".ann", ".amb", ".pac"}) {
        INFO("missing index file " << in.prefix << ext);
        REQUIRE(file_exists(in.prefix + ext));
    }
    for (const auto &path : {in.source, in.masked, in.bed}) {
        INFO("missing gate input " << path);
        REQUIRE(file_exists(path));
    }
    AltGateCounts counts;

    // Real ALT marking from each sidecar's names, against the same index.
    char dir[4096];
    REQUIRE(bwa_test_mkdtemp(dir, sizeof(dir), "alt_gate") == 0);
    struct DirGuard {
        const char *d;
        ~DirGuard() { bwa_test_remove_flat_dir(d); }
    } dir_guard{dir};
    BnsGuard src_bns(restore_with_alt(dir, "source", in.prefix, in.source));
    BnsGuard mask_bns(restore_with_alt(dir, "masked", in.prefix, in.masked));
    REQUIRE(src_bns.b != nullptr);
    REQUIRE(mask_bns.b != nullptr);
    REQUIRE(src_bns.b->n_seqs == mask_bns.b->n_seqs);
    CHECK(count_bns_differences(src_bns.b, mask_bns.b) == 0);
    // The deployed layout: the index's own sidecar, read by bns_restore(prefix).
    if (file_exists(in.prefix + ".alt")) {
        BnsGuard own_bns(bns_restore(in.prefix.c_str()));
        REQUIRE(own_bns.b != nullptr);
        REQUIRE(own_bns.b->n_seqs == src_bns.b->n_seqs);
        CHECK(count_bns_differences(own_bns.b, src_bns.b) == 0);
    }
    counts.n_seqs = src_bns.b->n_seqs;
    std::vector<int> alt_bits(counts.n_seqs);
    for (int rid = 0; rid < counts.n_seqs; ++rid) {
        alt_bits[rid] = src_bns.b->anns[rid].is_alt;
        counts.alt_contigs += alt_bits[rid] != 0;
    }

    // Row-level oracle shape: same QNAME order and non-FLAG bytes, FLAG either
    // unchanged or exactly FLAG | 4 on a record that lacked 0x4, unique ALT names.
    const auto src_rows = read_sidecar(in.source), mask_rows = read_sidecar(in.masked);
    REQUIRE(src_rows.size() == mask_rows.size());
    std::unordered_map<std::string, int> rid_of;
    rid_of.reserve(counts.n_seqs);
    for (int rid = 0; rid < counts.n_seqs; ++rid) rid_of.emplace(src_bns.b->anns[rid].name, rid);
    std::vector<int> row_rid(src_rows.size());
    std::vector<char> row_masked(src_rows.size());
    std::unordered_set<std::string> seen;
    seen.reserve(src_rows.size());
    long bad_shape = 0, bad_name = 0;
    for (size_t i = 0; i < src_rows.size(); ++i) {
        const auto &s = src_rows[i], &m = mask_rows[i];
        const bool shape_ok = s.qname == m.qname && s.rest == m.rest && !(s.flag & 4) &&
                              (m.flag == s.flag || m.flag == (s.flag | 4));
        if (!shape_ok) {
            if (bad_shape == 0) MESSAGE("first bad masked row " << i << " " << s.qname);
            ++bad_shape;
        }
        const auto it = rid_of.find(s.qname);
        const bool name_ok = it != rid_of.end() && alt_bits[it->second] && seen.insert(s.qname).second;
        if (!name_ok) {
            if (bad_name == 0) MESSAGE("first unknown, non-ALT or repeated QNAME at row " << i << " " << s.qname);
            ++bad_name;
        }
        row_rid[i] = it == rid_of.end() ? -1 : it->second;
        row_masked[i] = m.flag != s.flag;
    }
    CHECK(bad_shape == 0);
    REQUIRE(bad_name == 0);
    counts.rows = (long) src_rows.size();

    // Usability: unfiltered, every record is a table entry, so with unique names
    // each per-record decision is exactly one table entry's presence.
    LiftGuard plain(alt_lift_load(in.source.c_str(), src_bns.b));
    REQUIRE(plain.l != nullptr);
    const alt_lift_skips_t &ps = plain.l->skip;
    CHECK(ps.no_ctg + ps.no_target + ps.alt_target + ps.unmapped + ps.not_primary + ps.dup + ps.bad == 0);
    REQUIRE(plain.l->n_lift == counts.rows);

    alt_lift_exclude_stats_t stats = {};
    char error[1024] = {};
    LiftGuard runtime(alt_lift_load_excluded(in.source.c_str(), src_bns.b, in.bed.c_str(), &stats, error,
                                             sizeof(error)));
    INFO("restricted load error: " << error);
    REQUIRE(runtime.l != nullptr);  // unknown or out-of-bounds BED rows stop here
    LiftGuard oracle(alt_lift_load(in.masked.c_str(), mask_bns.b));
    REQUIRE(oracle.l != nullptr);
    if (tamper) tamper(oracle.l);

    // Per-record decisions against the oracle's FLAG bit.
    long bad_decision = 0;
    for (size_t i = 0; i < src_rows.size(); ++i) {
        const int rid = row_rid[i];
        const bool runtime_excluded = runtime.l->a[rid].trid < 0;
        const bool oracle_excluded = oracle.l->a[rid].trid < 0;
        if (runtime_excluded != (bool) row_masked[i] || oracle_excluded != (bool) row_masked[i]) {
            if (bad_decision == 0) MESSAGE("first differing decision at row " << i << " " << src_rows[i].qname);
            ++bad_decision;
        }
        if (row_masked[i]) ++counts.excluded; else ++counts.retained;
    }
    CHECK(bad_decision == 0);
    CHECK(stats.excluded == counts.excluded);
    CHECK(runtime.l->n_lift == counts.retained);
    CHECK(oracle.l->skip.unmapped == counts.excluded);

    // The full table, every contig, field by field.
    check_lift_tables(runtime.l, oracle.l);

    // Loading never changes ALT marking.
    for (const bntseq_t *b : {(const bntseq_t *) src_bns.b, (const bntseq_t *) mask_bns.b}) {
        long changed = 0;
        for (int rid = 0; rid < counts.n_seqs; ++rid) changed += b->anns[rid].is_alt != alt_bits[rid];
        CHECK(changed == 0);
    }
    counts.input_intervals = stats.input_intervals;
    counts.merged_intervals = stats.merged_intervals;
    return counts;
}

// A generated index with four ALT tiles: a forward tile and a tile whose
// deletion holds the only BED base in its footprint are excluded; a reverse tile
// and a tile with an insertion whose footprint ends where a BED interval starts
// are retained. Names include the alt_/altc_ prefixes and an HLA-style name.
struct AltGateFixture {
    std::string chr1, chr2;
    bwa_tests::TempIndex idx;
    TmpFile source, masked, bed;
    static std::string random_seq(std::mt19937 &rng, size_t n) {
        std::string s(n, 'A');
        for (char &c : s) c = "ACGT"[rng() % 4];
        return s;
    }
    static std::string sidecar(bool masked, int wrong_flag_row) {
        const int flags[4] = {0, 16, 0, 0};
        const bool excluded[4] = {true, false, true, false};
        const char *rest[4] = {
            "\tchr1\t101\t60\t200M\t*\t0\t0\t*\t*\n",
            "\tchr1\t1001\t60\t200M\t*\t0\t0\t*\t*\n",
            "\tchr2\t501\t60\t100M10D100M\t*\t0\t0\t*\t*\n",
            "\tchr1\t2001\t60\t75M5I75M\t*\t0\t0\t*\t*\n",
        };
        const char *names[4] = {"alt_chr1_t1", "alt_chr1_t2", "altc_chr2_t3", "HLA-A*02:01:01:01"};
        std::string out;
        for (int i = 0; i < 4; ++i) {
            const bool mask = masked && (excluded[i] != (i == wrong_flag_row));
            out += std::string(names[i]) + "\t" + std::to_string(flags[i] | (mask ? 4 : 0)) + rest[i];
        }
        return out;
    }
    static std::vector<bwa_tests::Contig> contigs(const std::string &c1, const std::string &c2) {
        return {{"chr1", c1},
                {"chr2", c2},
                {"alt_chr1_t1", c1.substr(100, 200)},
                {"alt_chr1_t2", bwa_tests::revcomp(c1.substr(1000, 200))},
                {"altc_chr2_t3", c2.substr(500, 100) + c2.substr(610, 100)},
                {"HLA-A*02:01:01:01", c1.substr(2000, 75) + "ACGTA" + c1.substr(2075, 75)}};
    }
    explicit AltGateFixture(int wrong_flag_row = -1, std::mt19937 rng = std::mt19937(7754))
        : chr1(random_seq(rng, 3000)), chr2(random_seq(rng, 2000)), idx(contigs(chr1, chr2)),
          source(sidecar(false, -1).c_str()), masked(sidecar(true, wrong_flag_row).c_str()),
          bed("chr1\t250\t260\nchr2\t605\t606\nchr1\t2150\t2200\n") {}
    AltGateInputs inputs() const { return {idx.prefix(), source.path, masked.path, bed.path}; }
};

const char *const kAltGateEnv[4] = {"BWAMEM3_ALT_GATE_PREFIX", "BWAMEM3_ALT_GATE_SOURCE",
                                    "BWAMEM3_ALT_GATE_MASKED", "BWAMEM3_ALT_GATE_BED"};

bool alt_gate_any_env_set() {
    for (const char *name : kAltGateEnv) {
        const char *v = getenv(name);
        if (v && *v) return true;
    }
    return false;
}

// Counts failed CHECK/REQUIRE assertions (not WARNs) over the whole run, so a
// test case can tell whether any assertion failed since a given point. doctest
// reports every assertion to listeners, including non-fatal CHECK failures that
// do not stop the case.
std::atomic<long> g_alt_gate_failed_asserts{0};

struct AltGateFailureCounter : doctest::IReporter {
    explicit AltGateFailureCounter(const doctest::ContextOptions &) {}
    void report_query(const doctest::QueryData &) override {}
    void test_run_start() override {}
    void test_run_end(const doctest::TestRunStats &) override {}
    void test_case_start(const doctest::TestCaseData &) override {}
    void test_case_reenter(const doctest::TestCaseData &) override {}
    void test_case_end(const doctest::CurrentTestCaseStats &) override {}
    void test_case_exception(const doctest::TestCaseException &) override {}
    void subcase_start(const doctest::SubcaseSignature &) override {}
    void subcase_end() override {}
    void log_assert(const doctest::AssertData &ad) override {
        if (ad.m_failed && !(ad.m_at & doctest::assertType::is_warn)) ++g_alt_gate_failed_asserts;
    }
    void log_message(const doctest::MessageData &) override {}
    void test_case_skipped(const doctest::TestCaseData &) override {}
};

// Writes the ALT_GATE_FROZEN_COMPLETE line to `out` only if no assertion has
// failed since `failed_before` was read; returns whether it wrote the line.
bool alt_gate_emit_marker(FILE *out, long failed_before, const AltGateCounts &c) {
    if (g_alt_gate_failed_asserts.load() != failed_before) return false;
    fprintf(out,
            "ALT_GATE_FROZEN_COMPLETE n_seqs=%ld alt_contigs=%ld rows=%ld excluded=%ld retained=%ld "
            "bed_intervals=%ld merged_intervals=%ld\n",
            c.n_seqs, c.alt_contigs, c.rows, c.excluded, c.retained, c.input_intervals, c.merged_intervals);
    fflush(out);
    return true;
}

// Reads back what alt_gate_emit_marker wrote to a temporary file.
std::string alt_gate_marker_text(long failed_before, const AltGateCounts &c, bool *emitted) {
    FILE *f = tmpfile();
    REQUIRE(f != nullptr);
    *emitted = alt_gate_emit_marker(f, failed_before, c);
    rewind(f);
    std::string text;
    char buf[256];
    while (fgets(buf, sizeof(buf), f)) text += buf;
    fclose(f);
    return text;
}
} // namespace

DOCTEST_REGISTER_LISTENER("alt_gate_failure_counter", 1, AltGateFailureCounter);

#endif  // BWA_TEST_ALT_LIFT_GATE_H
