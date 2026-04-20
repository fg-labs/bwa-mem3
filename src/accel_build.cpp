// bwa-mem2 build-accel: construct a target-region fast-path cache file.
//
// Pipeline:
//   1. Parse CLI.
//   2. Load bwa-mem2 FM-index via FMI_search.
//   3. SHA-256 the .bwt.2bit.64 file (stored in cache header).
//   4. Determine the terminal k-mer set:
//        - BED: extract k-mers from reference regions (optionally with flanks).
//        - FASTQ: scan reads, collect k-mers.
//        - Intersect FASTQ k-mers with genome via SA-interval resolution
//          (k-mers with s=0 at length k_max are rejected).
//   5. For each terminal canonical k-mer, walk the FM-index to record
//      (k, l, s) at every length L = 1..k_max.
//   6. Aggregate per-level tables: dedup by canonical L-mer; drop entries
//      whose s <= collapse_threshold.
//   7. Serialize as a cache file (layout matches src/accel_cache.cpp's
//      loader).
//
// Non-goals for this revision: streaming build for huge inputs, hotspots
// bundle, BAM input. See docs/superpowers/plans/... for future work.

#include "accel_build.h"
#include "accel_cache.h"
#include "FMI_search.h"

#include <zlib.h>
#include "kseq.h"
#ifndef KSEQ_INITED_ACCEL_BUILD
#define KSEQ_INITED_ACCEL_BUILD
KSEQ_INIT(gzFile, gzread)
#endif

#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cerrno>
#include <getopt.h>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <unordered_set>
#include <unordered_map>
#include <algorithm>

extern "C" {
#include "../ext/sha256/sha-256.h"
}

namespace {

// -- Options ------------------------------------------------------------

struct BuildOpts {
    std::string ref_path;                       // required
    std::vector<std::string> beds;              // zero or more
    std::vector<std::string> reads;             // zero or more
    std::string out_path;                       // required
    int flanks_bp = 0;
    int k_min = 1;
    int k_max = 18;
    int collapse_thresh = 1;                    // drop entries with s<=thresh from per-level tables
    int n_threads = 1;                          // TODO: parallel trajectory resolve
};

static void usage(FILE *fp) {
    fprintf(fp,
"Usage: bwa-mem2 build-accel --ref REF.fa --out CACHE [options]\n"
"\n"
"Required:\n"
"  --ref PATH             Reference FASTA (must have a bwa-mem2 index\n"
"                         with sibling .bwt.2bit.64 etc.)\n"
"  --out PATH             Output cache file path\n"
"\n"
"At least one of these must be given (can be combined):\n"
"  --bed PATH             BED file of target regions (may be repeated)\n"
"  --reads PATH           FASTQ or FASTQ.gz of reads (may be repeated)\n"
"\n"
"Options:\n"
"  --flanks BP            Extend each BED interval by BP on each side (default 0)\n"
"  --k-min L              Minimum k-mer length to index (default 1)\n"
"  --k-max L              Maximum k-mer length (default 18, max 18)\n"
"  --collapse-threshold S Drop per-level entries with s <= S (default 1)\n"
"  --threads N            Number of threads (default 1; unused in this revision)\n"
"\n"
);
}

static int parse_args(int argc, char **argv, BuildOpts &o) {
    // argv[0] is the subcommand name ("build-accel") per main.cpp's shift,
    // and getopt_long treats argv[0] as the program name, so we don't
    // manually shift here.
    static const struct option longopts[] = {
        {"ref",                 required_argument, 0, 'r'},
        {"out",                 required_argument, 0, 'o'},
        {"bed",                 required_argument, 0, 'b'},
        {"reads",               required_argument, 0, 'R'},
        {"flanks",              required_argument, 0, 'f'},
        {"k-min",               required_argument, 0, 'm'},
        {"k-max",               required_argument, 0, 'M'},
        {"collapse-threshold",  required_argument, 0, 'c'},
        {"threads",             required_argument, 0, 't'},
        {"help",                no_argument,       0, 'h'},
        {0, 0, 0, 0}
    };
    int c;
    int long_index = 0;
    optind = 1;  // reset for use after subcommand dispatch
    while ((c = getopt_long(argc, argv, "r:o:b:R:f:m:M:c:t:h", longopts, &long_index)) != -1) {
        switch (c) {
            case 'r': o.ref_path = optarg; break;
            case 'o': o.out_path = optarg; break;
            case 'b': o.beds.push_back(optarg); break;
            case 'R': o.reads.push_back(optarg); break;
            case 'f': o.flanks_bp = atoi(optarg); break;
            case 'm': o.k_min = atoi(optarg); break;
            case 'M': o.k_max = atoi(optarg); break;
            case 'c': o.collapse_thresh = atoi(optarg); break;
            case 't': o.n_threads = atoi(optarg); break;
            case 'h': usage(stdout); return 2;
            default:  usage(stderr); return 1;
        }
    }
    if (o.ref_path.empty() || o.out_path.empty()) {
        fprintf(stderr, "[build-accel] --ref and --out are required\n");
        return 1;
    }
    if (o.beds.empty() && o.reads.empty()) {
        fprintf(stderr, "[build-accel] provide at least one of --bed or --reads\n");
        return 1;
    }
    if (o.k_min < 1 || o.k_max < o.k_min || o.k_max > (int)ACCEL_K_MAX) {
        fprintf(stderr, "[build-accel] invalid k range [%d, %d]\n", o.k_min, o.k_max);
        return 1;
    }
    return 0;
}

// -- Base encoding ------------------------------------------------------

static inline int base_to_enc(char c) {
    switch (c) {
        case 'A': case 'a': return 0;
        case 'C': case 'c': return 1;
        case 'G': case 'g': return 2;
        case 'T': case 't': return 3;
        default:            return 4; // N or other: reject
    }
}

// -- FASTA loader (whole-reference in memory by contig) ----------------

struct Contig {
    std::string name;
    std::string bases;  // ACGTN characters
};

static bool load_fasta(const std::string &path, std::vector<Contig> &out) {
    gzFile fp = gzopen(path.c_str(), "rb");
    if (!fp) {
        fprintf(stderr, "[build-accel] cannot open FASTA %s: %s\n",
                path.c_str(), strerror(errno));
        return false;
    }
    kseq_t *ks = kseq_init(fp);
    int l;
    while ((l = kseq_read(ks)) >= 0) {
        Contig c;
        c.name.assign(ks->name.s);
        c.bases.assign(ks->seq.s, (size_t)l);
        out.push_back(std::move(c));
    }
    kseq_destroy(ks);
    gzclose(fp);
    return true;
}

// -- BED → k-mers ------------------------------------------------------

static bool load_bed_regions(const std::string &bed_path,
                             std::vector<std::tuple<std::string, long long, long long>> &out) {
    std::ifstream f(bed_path);
    if (!f.is_open()) {
        fprintf(stderr, "[build-accel] cannot open BED %s\n", bed_path.c_str());
        return false;
    }
    std::string line;
    while (std::getline(f, line)) {
        if (line.empty() || line[0] == '#' || line[0] == 't') continue; // skip "track" lines
        std::istringstream s(line);
        std::string chrom; long long start, end;
        if (!(s >> chrom >> start >> end)) continue;
        out.push_back({chrom, start, end});
    }
    return true;
}

// Slide a k-length window across `bases[begin..end)` and insert canonical
// k-mers into the output set.
static void slide_kmers(const std::string &bases, size_t begin, size_t end,
                        int k, std::unordered_set<uint64_t> &out) {
    if ((long long)end - (long long)begin < k) return;
    std::vector<uint8_t> enc((size_t)k);
    for (size_t i = begin; i + (size_t)k <= end; ++i) {
        bool ok = true;
        for (int j = 0; j < k; ++j) {
            int b = base_to_enc(bases[i + (size_t)j]);
            if (b > 3) { ok = false; break; }
            enc[(size_t)j] = (uint8_t)b;
        }
        if (!ok) continue;
        uint64_t f = accel_pack_forward(enc.data(), (uint32_t)k);
        if (f == UINT64_MAX) continue;
        out.insert(accel_canonicalize(f, (uint32_t)k));
    }
}

static bool collect_bed_kmers(const BuildOpts &opt,
                              const std::vector<Contig> &contigs,
                              std::unordered_set<uint64_t> &kmers) {
    std::unordered_map<std::string, const Contig*> by_name;
    for (const auto &c : contigs) by_name[c.name] = &c;

    for (const auto &bed : opt.beds) {
        std::vector<std::tuple<std::string, long long, long long>> regions;
        if (!load_bed_regions(bed, regions)) return false;
        fprintf(stderr, "[build-accel] %s: %zu BED intervals\n",
                bed.c_str(), regions.size());
        for (auto &tpl : regions) {
            const std::string &chrom = std::get<0>(tpl);
            long long s = std::get<1>(tpl);
            long long e = std::get<2>(tpl);
            auto it = by_name.find(chrom);
            if (it == by_name.end()) continue;
            const Contig *c = it->second;
            long long s0 = std::max<long long>(0, s - opt.flanks_bp);
            long long e0 = std::min<long long>((long long)c->bases.size(),
                                               e + opt.flanks_bp);
            if (s0 >= e0) continue;
            slide_kmers(c->bases, (size_t)s0, (size_t)e0, opt.k_max, kmers);
        }
    }
    return true;
}

// -- FASTQ → k-mers ----------------------------------------------------

static bool collect_fastq_kmers(const BuildOpts &opt,
                                std::unordered_set<uint64_t> &kmers) {
    for (const auto &path : opt.reads) {
        gzFile fp = gzopen(path.c_str(), "rb");
        if (!fp) {
            fprintf(stderr, "[build-accel] cannot open %s: %s\n",
                    path.c_str(), strerror(errno));
            return false;
        }
        kseq_t *ks = kseq_init(fp);
        int l;
        size_t n_reads = 0;
        while ((l = kseq_read(ks)) >= 0) {
            if (l < opt.k_max) continue;
            std::string seq(ks->seq.s, (size_t)l);
            slide_kmers(seq, 0, seq.size(), opt.k_max, kmers);
            n_reads++;
        }
        kseq_destroy(ks);
        gzclose(fp);
        fprintf(stderr, "[build-accel] %s: %zu reads\n", path.c_str(), n_reads);
    }
    return true;
}

// -- Trajectory resolution --------------------------------------------

// Per-length record during the walk. We keep (canonical_L_mer, k, l, s).
struct LevelRecord {
    uint64_t canonical;
    int64_t  k, l, s;
};

// Unpack k bases from a packed kmer (low-to-high 2 bits/base).
static void unpack_bases(uint64_t packed, int k, uint8_t *out) {
    for (int i = 0; i < k; ++i) {
        out[i] = (uint8_t)((packed >> (2 * i)) & 0x3u);
    }
}

// Walk one k-mer (bases given) and record the trajectory at each L=1..k_max.
// SA intervals (k, l) are for the forward strand of the INPUT bases.
//
// At each level L, records:
//   - canonical(bases[0..L-1])
//   - the (k, l) that correspond to that canonical form (swapped from
//     the input-forward (k, l) if canonical != input-forward)
// This way, a runtime probe using canonicalize(read[x..x+L-1]) will get
// back (k, l) that the runtime can swap based on its own is_rc flag —
// the same swap semantics as FMI_search.cpp:541-547.
static void walk_bases(FMI_search *fmi,
                       const uint8_t *bases,
                       int k_max,
                       std::vector<LevelRecord> &out) {
    out.clear();
    out.resize((size_t)k_max);

    const int64_t *count = fmi->accel_count();
    SMEM smem;
    smem.rid = 0;
    smem.m = 0;
    smem.n = 0;
    uint8_t a = bases[0];
    smem.k = count[a];
    smem.l = count[3 - a];
    smem.s = count[a + 1] - count[a];

    // Record L=1.
    {
        uint64_t fwd1 = accel_pack_forward(bases, 1);
        bool is_rc = false;
        uint64_t c1 = accel_canonicalize(fwd1, 1, &is_rc);
        // (k_canonical, l_canonical) = swap if is_rc
        int64_t k_c = is_rc ? smem.l : smem.k;
        int64_t l_c = is_rc ? smem.k : smem.l;
        out[0] = {c1, k_c, l_c, smem.s};
    }

    for (int L = 2; L <= k_max; ++L) {
        uint8_t next_base = bases[L - 1];
        SMEM smem_ = smem;
        smem_.k = smem.l;
        smem_.l = smem.k;
        SMEM newSmem_ = fmi->accel_backwardExt(smem_, (uint8_t)(3 - next_base));
        SMEM newSmem = newSmem_;
        newSmem.k = newSmem_.l;
        newSmem.l = newSmem_.k;
        newSmem.n = L - 1;
        smem = newSmem;

        uint64_t fwdL = accel_pack_forward(bases, (uint32_t)L);
        bool is_rc = false;
        uint64_t cL = accel_canonicalize(fwdL, (uint32_t)L, &is_rc);
        int64_t k_c = is_rc ? smem.l : smem.k;
        int64_t l_c = is_rc ? smem.k : smem.l;
        out[(size_t)(L - 1)] = {cL, k_c, l_c, smem.s};

        if (smem.s == 0) {
            for (int L2 = L + 1; L2 <= k_max; ++L2) {
                out[(size_t)(L2 - 1)] = {0, 0, 0, 0};
            }
            return;
        }
    }
}

// For a canonical 18-mer, walk BOTH the forward (canonical bases) and
// reverse-complement (RC of canonical) orientations. Runtime reads may
// arrive in either orientation, and the L-mer prefixes of the two
// orientations canonicalize to DIFFERENT canonical forms in general
// (they only coincide for palindromes). We therefore need both
// trajectories to cover both possible read orientations at the per-L
// tables.
//
// `traj_fwd` and `traj_rc` are filled with the length-k_max trajectories
// for canonical's bases and RC(canonical)'s bases respectively.
static void walk_both_orientations(FMI_search *fmi,
                                   uint64_t canonical_kmer,
                                   int k_max,
                                   std::vector<LevelRecord> &traj_fwd,
                                   std::vector<LevelRecord> &traj_rc) {
    uint8_t bases[64];
    unpack_bases(canonical_kmer, k_max, bases);
    walk_bases(fmi, bases, k_max, traj_fwd);

    // RC bases: reverse order and complement each.
    uint8_t rc_bases[64];
    for (int i = 0; i < k_max; ++i) {
        rc_bases[i] = (uint8_t)(bases[k_max - 1 - i] ^ 0x3u);
    }
    walk_bases(fmi, rc_bases, k_max, traj_rc);
}

// -- Aggregation -------------------------------------------------------

struct LevelTableBuilder {
    // canonical_L_mer → (k, l, s). Last-write-wins is safe because
    // walks are deterministic for a given prefix.
    std::unordered_map<uint64_t, AccelEntry> entries;
};

struct AggregateResult {
    std::vector<LevelTableBuilder> per_level;  // indexed 0..k_max
    LevelTableBuilder terminal;
};

static void aggregate(const std::vector<LevelRecord> &traj,
                      int k_min, int k_max, int collapse_thresh,
                      AggregateResult &agg) {
    for (int L = k_min; L <= k_max; ++L) {
        const LevelRecord &rec = traj[(size_t)(L - 1)];
        if (rec.s <= collapse_thresh) continue;
        AccelEntry e{rec.canonical, rec.k, rec.l, rec.s};
        agg.per_level[(size_t)L].entries[rec.canonical] = e;
    }
    const LevelRecord &term = traj[(size_t)(k_max - 1)];
    // Terminal: keep the terminal state for every k-mer that exists in the
    // reference (s>0). s==0 means the k-mer isn't in the genome; skip.
    if (term.s > 0) {
        AccelEntry e{term.canonical, term.k, term.l, term.s};
        agg.terminal.entries[term.canonical] = e;
    }
}

// -- Serialization ----------------------------------------------------

static inline uint64_t splitmix64(uint64_t x) {
    x += 0x9E3779B97F4A7C15ULL;
    x = (x ^ (x >> 30)) * 0xBF58476D1CE4E5B9ULL;
    x = (x ^ (x >> 27)) * 0x94D049BB133111EBULL;
    return x ^ (x >> 31);
}

static uint64_t next_pow2(uint64_t n) {
    uint64_t p = 1;
    while (p < n) p <<= 1;
    return p;
}

struct SerializedTable {
    std::vector<uint64_t> probe;   // slots
    std::vector<AccelEntry> entries;
    uint64_t table_mask;
};

static SerializedTable serialize_table(const LevelTableBuilder &b) {
    SerializedTable t;
    size_t n = b.entries.size();
    uint64_t table_size = next_pow2(std::max<uint64_t>(4, (uint64_t)(2 * n)));
    t.table_mask = table_size - 1;
    t.probe.assign(table_size, 0);

    // Flatten entries to a vector in deterministic order (sort by canonical).
    t.entries.reserve(n);
    std::vector<std::pair<uint64_t, const AccelEntry*>> sorted;
    sorted.reserve(n);
    for (const auto &kv : b.entries) sorted.push_back({kv.first, &kv.second});
    std::sort(sorted.begin(), sorted.end(),
              [](auto &a, auto &b){ return a.first < b.first; });
    for (size_t i = 0; i < sorted.size(); ++i) {
        t.entries.push_back(*sorted[i].second);
        uint64_t slot = splitmix64(sorted[i].first) & t.table_mask;
        while (t.probe[slot] != 0) slot = (slot + 1) & t.table_mask;
        t.probe[slot] = (uint64_t)(i + 1);
    }
    return t;
}

static bool write_cache(const BuildOpts &opt,
                        const uint8_t ref_sha256[32],
                        const AggregateResult &agg) {
    FILE *f = fopen(opt.out_path.c_str(), "wb");
    if (!f) {
        fprintf(stderr, "[build-accel] cannot write %s: %s\n",
                opt.out_path.c_str(), strerror(errno));
        return false;
    }
    // File header.
    uint64_t magic = ACCEL_MAGIC;              fwrite(&magic, 8, 1, f);
    uint32_t ver   = ACCEL_VERSION;            fwrite(&ver, 4, 1, f);
    fwrite(ref_sha256, 32, 1, f);
    uint32_t kmin = (uint32_t)opt.k_min;       fwrite(&kmin, 4, 1, f);
    uint32_t kmax = (uint32_t)opt.k_max;       fwrite(&kmax, 4, 1, f);
    uint32_t tier_count = 1;                   fwrite(&tier_count, 4, 1, f);
    uint32_t collapse = (uint32_t)opt.collapse_thresh;
                                               fwrite(&collapse, 4, 1, f);
    uint8_t reserved[16] = {};                 fwrite(reserved, 16, 1, f);

    // Tier header.
    uint8_t tier_id = ACCEL_TIER_B;            fwrite(&tier_id, 1, 1, f);
    uint8_t tier_pad[7] = {};                  fwrite(tier_pad, 7, 1, f);

    // Serialize per-level tables.
    for (int L = opt.k_min; L <= opt.k_max; ++L) {
        SerializedTable t = serialize_table(agg.per_level[(size_t)L]);
        uint64_t n_e = (uint64_t)t.entries.size();
        fwrite(&n_e, 8, 1, f);
        fwrite(&t.table_mask, 8, 1, f);
        fwrite(t.probe.data(), sizeof(uint64_t), t.probe.size(), f);
        fwrite(t.entries.data(), sizeof(AccelEntry), t.entries.size(), f);
        fprintf(stderr, "[build-accel] L=%2d: %llu entries, table %llu\n",
                L, (unsigned long long)n_e,
                (unsigned long long)t.probe.size());
    }
    // Terminal table.
    SerializedTable tt = serialize_table(agg.terminal);
    uint64_t tn = (uint64_t)tt.entries.size();
    fwrite(&tn, 8, 1, f);
    fwrite(&tt.table_mask, 8, 1, f);
    fwrite(tt.probe.data(), sizeof(uint64_t), tt.probe.size(), f);
    fwrite(tt.entries.data(), sizeof(AccelEntry), tt.entries.size(), f);
    fprintf(stderr, "[build-accel] terminal: %llu entries, table %llu\n",
            (unsigned long long)tn,
            (unsigned long long)tt.probe.size());

    fclose(f);
    return true;
}

} // namespace

// -- Main entry --------------------------------------------------------

int accel_build_main(int argc, char **argv) {
    BuildOpts opt;
    int pr = parse_args(argc, argv, opt);
    if (pr == 2) return 0;
    if (pr != 0) return pr;

    fprintf(stderr, "[build-accel] ref=%s out=%s k=[%d..%d]\n",
            opt.ref_path.c_str(), opt.out_path.c_str(), opt.k_min, opt.k_max);

    // SHA-256 of the .bwt.2bit.64.
    std::string bwt_path = opt.ref_path + ".bwt.2bit.64";
    uint8_t ref_sha[32];
    if (!accel_sha256_file(bwt_path.c_str(), ref_sha)) {
        fprintf(stderr, "[build-accel] cannot hash %s: %s\n",
                bwt_path.c_str(), strerror(errno));
        return 2;
    }
    {
        char hex[65];
        for (int i = 0; i < 32; ++i) snprintf(hex + 2*i, 3, "%02x", ref_sha[i]);
        hex[64] = 0;
        fprintf(stderr, "[build-accel] ref sha256 = %s\n", hex);
    }

    // Load FM-index.
    FMI_search *fmi = new FMI_search(opt.ref_path.c_str());
    fmi->load_index();

    // Build terminal k-mer set.
    std::unordered_set<uint64_t> kmers;
    if (!opt.beds.empty()) {
        std::vector<Contig> contigs;
        if (!load_fasta(opt.ref_path, contigs)) { delete fmi; return 3; }
        fprintf(stderr, "[build-accel] reference: %zu contigs\n", contigs.size());
        if (!collect_bed_kmers(opt, contigs, kmers)) { delete fmi; return 3; }
    }
    if (!opt.reads.empty()) {
        if (!collect_fastq_kmers(opt, kmers)) { delete fmi; return 3; }
    }
    fprintf(stderr, "[build-accel] %zu distinct canonical %d-mers from inputs\n",
            kmers.size(), opt.k_max);

    // Resolve trajectories + aggregate.
    AggregateResult agg;
    agg.per_level.resize((size_t)(opt.k_max + 1));

    size_t done = 0;
    size_t kept = 0;
    std::vector<LevelRecord> traj_fwd, traj_rc;
    for (uint64_t km : kmers) {
        walk_both_orientations(fmi, km, opt.k_max, traj_fwd, traj_rc);
        const LevelRecord &term = traj_fwd[(size_t)(opt.k_max - 1)];
        if (term.s > 0) {
            aggregate(traj_fwd, opt.k_min, opt.k_max, opt.collapse_thresh, agg);
            aggregate(traj_rc,  opt.k_min, opt.k_max, opt.collapse_thresh, agg);
            kept++;
        }
        done++;
        if ((done % 100000) == 0) {
            fprintf(stderr, "[build-accel]   resolved %zu / %zu kmers (%zu in genome)\n",
                    done, kmers.size(), kept);
        }
    }
    fprintf(stderr, "[build-accel] resolved %zu kmers, %zu in genome\n",
            done, kept);

    // Write.
    if (!write_cache(opt, ref_sha, agg)) { delete fmi; return 4; }

    delete fmi;
    fprintf(stderr, "[build-accel] DONE: wrote %s\n", opt.out_path.c_str());
    return 0;
}
