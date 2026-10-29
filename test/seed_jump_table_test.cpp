/* Exhaustive gate for the seed jump table (FMI_search::seed_jump).
 *
 * The lockstep parity tests only read the entries of strings that occur in
 * their reads. This test checks every entry of every level against the
 * interval the seeding walks reach by stepping: the 1-base interval
 * (count[a], count[3 - a], count[a + 1] - count[a]), then one forward step per
 * base (swap k <-> l, backwardExt with 3 - c, swap back), s == 0 strings and
 * strings whose interval holds the sentinel row included. It does so at the
 * depths on both sides of the builder's serial prefix (2, 3, 4) and at the
 * default depth, each at several builder thread counts, so a wrong level
 * offset, child index or per-prefix partition fails here whatever the reads.
 * At depths 2-4 it also checks every entry's occurrence count s, and the
 * interval starts k and l, against a brute-force count and a naive suffix
 * array of the index text (the .pac forward strand plus its reverse
 * complement), so a defect that backwardExt or the 1-base interval shares
 * with the stepped walk above fails too.
 * It checks that the lockstep SMEM and bwtSeed drivers take their first
 * steps from the table at the default depth: with the table's entries
 * replaced by other valid intervals, their output changes.
 * It also checks the depth BWA3_BWTSEED_JUMP selects (unset, empty, 0,
 * in-range and malformed values), the build_seed_jump=false opt-out on
 * load_index and on fmi_seed_open, the table an index attached from a packed
 * shm segment builds, and load_index's shm-attach path with and without the
 * opt-out (when a phiX segment can be staged on this host).
 *
 * Usage: seed_jump_table_test <bwa-mem3 index prefix>
 */

#define private public
#include "FMI_search.h"
#undef private
#include "bwa.h"
#include "bwa_shm.h"
#include "bntseq.h"
#include "fmi_seed_api.h"
#include "seed_jump_test_util.h"

#include <algorithm>
#include <string>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <sys/mman.h>
#include <fcntl.h>
#include <unistd.h>

static int32_t total_cases = 0;
static int32_t passed_cases = 0;

static void report(bool ok, const char *name)
{
    total_cases++;
    if (ok) passed_cases++;
    fprintf(stderr, "[%s] %s\n", ok ? "PASS" : "FAIL", name);
}

/* Set BWA3_BWTSEED_JUMP (NULL = unset). */
static void set_jump_env(const char *v)
{
    if (v == NULL) unsetenv("BWA3_BWTSEED_JUMP");
    else setenv("BWA3_BWTSEED_JUMP", v, 1);
}

struct Walk { int64_t k, l, s; };

/* One forward step of the seeding walks (bsd_advance_step,
 * ls_advance_forward_step) from w with base c. */
static Walk forward_step(FMI_search *fmi, const Walk &w, int c)
{
    SMEM t;
    memset(&t, 0, sizeof(t));
    t.k = w.l;
    t.l = w.k;
    t.s = w.s;
    const SMEM r = fmi->backwardExt(t, (uint8_t)(3 - c));
    Walk o = { r.l, r.k, r.s };
    return o;
}

/* Depth-first over every string of length L..D extending `code`, comparing
 * each table entry with the stepped interval w. */
static void check_subtree(FMI_search *fmi, int L, uint32_t code, const Walk &w, int D,
                          int64_t *checked, int64_t *bad)
{
    const FMI_search::SeedJumpEntry &e = fmi->seed_jump[FMI_search::seed_jump_level_off(L) + code];
    (*checked)++;
    if (e.k != w.k || e.l != w.l || e.s != w.s) {
        if (*bad < 5)
            fprintf(stderr, "  level %d code %u: table (%lld, %lld, %lld) != stepped (%lld, %lld, %lld)\n",
                    L, code, (long long)e.k, (long long)e.l, (long long)e.s,
                    (long long)w.k, (long long)w.l, (long long)w.s);
        (*bad)++;
    }
    if (L == D) return;
    for (int c = 0; c < 4; c++)
        check_subtree(fmi, L + 1, code * 4 + (uint32_t)c, forward_step(fmi, w, c), D, checked, bad);
}

/* Every entry of fmi's table at depth D. */
static bool check_table(FMI_search *fmi, int D)
{
    if (fmi->seed_jump == NULL || fmi->seed_jump_depth != D) {
        fprintf(stderr, "  table %p, depth %d; expected depth %d\n",
                (void *)fmi->seed_jump, fmi->seed_jump_depth, D);
        return false;
    }
    int64_t checked = 0, bad = 0;
    for (int a = 0; a < 4; a++) {
        Walk w = { fmi->count[a], fmi->count[3 - a], fmi->count[a + 1] - fmi->count[a] };
        check_subtree(fmi, 1, (uint32_t)a, w, D, &checked, &bad);
    }
    const int64_t want = FMI_search::seed_jump_level_off(D + 1);
    if (checked != want) {
        fprintf(stderr, "  checked %lld entries, expected %lld\n", (long long)checked, (long long)want);
        return false;
    }
    return bad == 0;
}

/* ---- Brute-force oracle over the index text -------------------------------
 * The index is built over the .pac forward strand followed by its reverse
 * complement, with a sentinel ($) after it that sorts below every base. */
struct TextOracle {
    std::vector<uint8_t> t;      // 2-bit bases, forward then reverse complement
    std::vector<int32_t> sa;     // suffix array of t + $ (the $ suffix first)
};

static bool load_text(const char *prefix, TextOracle *o)
{
    bntseq_t *bns = bns_restore(prefix);
    if (bns == NULL) return false;
    const int64_t l_pac = bns->l_pac;
    bns_destroy(bns);
    std::string pac_name = std::string(prefix) + ".pac";
    FILE *fp = fopen(pac_name.c_str(), "rb");
    if (fp == NULL) return false;
    std::vector<uint8_t> pac((size_t)(l_pac / 4 + 1));
    const size_t got = fread(pac.data(), 1, pac.size(), fp);
    fclose(fp);
    if ((int64_t)got * 4 < l_pac) return false;
    o->t.resize((size_t)(2 * l_pac));
    for (int64_t i = 0; i < l_pac; i++) {
        const uint8_t b = (pac[(size_t)(i >> 2)] >> ((~i & 3) << 1)) & 3;
        o->t[(size_t)i] = b;
        o->t[(size_t)(2 * l_pac - 1 - i)] = (uint8_t)(3 - b);
    }
    const int32_t n = (int32_t)o->t.size();
    o->sa.resize((size_t)n + 1);
    for (int32_t i = 0; i <= n; i++) o->sa[(size_t)i] = i;
    const uint8_t *t = o->t.data();
    std::sort(o->sa.begin(), o->sa.end(), [t, n](int32_t a, int32_t b) {
        while (a < n && b < n) {
            if (t[a] != t[b]) return t[a] < t[b];
            a++; b++;
        }
        return a == n && b != n;   // the shorter suffix (reaching $) sorts first
    });
    return true;
}

/* Occurrences of p in t, counted position by position. */
static int64_t brute_count(const TextOracle &o, const std::vector<uint8_t> &p)
{
    int64_t c = 0;
    const size_t n = o.t.size(), L = p.size();
    for (size_t i = 0; i + L <= n; i++)
        if (std::equal(p.begin(), p.end(), o.t.begin() + (long)i)) c++;
    return c;
}

/* Rank of the first suffix of t + $ that is not below p (p's SA interval start). */
static int64_t sa_lower(const TextOracle &o, const std::vector<uint8_t> &p)
{
    const int32_t n = (int32_t)o.t.size();
    const uint8_t *t = o.t.data();
    auto it = std::lower_bound(o.sa.begin(), o.sa.end(), 0, [&](int32_t suf, int) {
        for (size_t i = 0; i < p.size(); i++) {
            if (suf + (int32_t)i >= n) return true;   // $ sorts below every base
            const uint8_t b = t[suf + i];
            if (b != p[i]) return b < p[i];
        }
        return false;   // p is a prefix of the suffix: not below p
    });
    return (int64_t)(it - o.sa.begin());
}

/* Every entry of fmi's table (depth D) against the text oracle: s is p's
 * occurrence count; when s > 0, k is p's SA interval start and l is the start
 * of p's reverse complement (both forward-strand and reverse-complement
 * occurrences are in the text, so the two counts agree). */
static bool check_table_text(FMI_search *fmi, const TextOracle &o, int D)
{
    if (fmi->seed_jump == NULL || fmi->seed_jump_depth != D) return false;
    int64_t checked = 0, bad = 0;
    for (int L = 1; L <= D; L++) {
        for (uint32_t code = 0; code < (1u << (2 * L)); code++) {
            std::vector<uint8_t> p((size_t)L), rc((size_t)L);
            for (int i = 0; i < L; i++) p[(size_t)i] = (code >> (2 * (L - 1 - i))) & 3;
            for (int i = 0; i < L; i++) rc[(size_t)i] = (uint8_t)(3 - p[(size_t)(L - 1 - i)]);
            const FMI_search::SeedJumpEntry &e = fmi->seed_jump[FMI_search::seed_jump_level_off(L) + code];
            const int64_t s = brute_count(o, p);
            bool ok = e.s == s;
            if (ok && s > 0) ok = e.k == sa_lower(o, p) && e.l == sa_lower(o, rc);
            checked++;
            if (!ok) {
                if (bad < 5)
                    fprintf(stderr, "  level %d code %u: table (%lld, %lld, %lld), text s %lld k %lld l %lld\n",
                            L, code, (long long)e.k, (long long)e.l, (long long)e.s, (long long)s,
                            (long long)sa_lower(o, p), (long long)sa_lower(o, rc));
                bad++;
            }
        }
    }
    return checked == FMI_search::seed_jump_level_off(D + 1) && bad == 0;
}

/* ---- Does a lockstep driver read the table? --------------------------------
 * Reads: phiX windows (forward and reverse-complement text) and the jump-depth
 * reads, 2-bit encoded, N = 4. */
struct Reads {
    std::vector<uint8_t> enc;
    std::vector<bseq1_t> seq;
    std::vector<int32_t> cum;
    int32_t max_len = 0;
};

static Reads make_reads(const TextOracle &o)
{
    std::vector<std::vector<uint8_t>> rs;
    for (size_t st = 0; st + 60 <= o.t.size() && rs.size() < 40; st += o.t.size() / 40)
        rs.push_back(std::vector<uint8_t>(o.t.begin() + (long)st, o.t.begin() + (long)st + 60));
    for (const std::string &r : jump_depth_reads()) {
        std::vector<uint8_t> e;
        for (char c : r) e.push_back(c == 'A' ? 0 : c == 'C' ? 1 : c == 'G' ? 2 : c == 'T' ? 3 : 4);
        rs.push_back(e);
    }
    Reads R;
    R.seq.resize(rs.size());
    memset(R.seq.data(), 0, rs.size() * sizeof(bseq1_t));
    for (size_t i = 0; i < rs.size(); i++) {
        R.cum.push_back((int32_t)R.enc.size());
        R.enc.insert(R.enc.end(), rs[i].begin(), rs[i].end());
        R.seq[i].l_seq = (int)rs[i].size();
        R.max_len = std::max(R.max_len, (int32_t)rs[i].size());
    }
    R.cum.push_back((int32_t)R.enc.size());
    return R;
}

static bool same_smems(const std::vector<SMEM> &a, const std::vector<SMEM> &b)
{
    if (a.size() != b.size()) return false;
    for (size_t i = 0; i < a.size(); i++)
        if (a[i].rid != b[i].rid || a[i].m != b[i].m || a[i].n != b[i].n
            || a[i].k != b[i].k || a[i].l != b[i].l || a[i].s != b[i].s
            || a[i].info != b[i].info)
            return false;
    return true;
}

/* The SMEM lockstep driver's output: one forward walk from x = 0 per read. */
static std::vector<SMEM> run_smem_lockstep(FMI_search *fmi, Reads &R, int32_t min_intv)
{
    const int32_t n = (int32_t)R.seq.size();
    std::vector<int32_t> qpa((size_t)n, 0), mia((size_t)n, min_intv), rid((size_t)n);
    for (int32_t i = 0; i < n; i++) rid[(size_t)i] = i;
    std::vector<SMEM> out((size_t)n * (size_t)R.max_len);
    int64_t n_out = 0;
    fmi->getSMEMsOnePosOneThread_lockstep(R.enc.data(), qpa.data(), mia.data(), rid.data(),
                                          n, n, R.seq.data(), R.cum.data(), R.max_len,
                                          /*minSeedLen=*/19, out.data(), &n_out);
    out.resize((size_t)n_out);
    return out;
}

/* The bwtSeed lockstep driver's output at the default minSeedLen. */
static std::vector<SMEM> run_bwtseed_lockstep(FMI_search *fmi, Reads &R, int32_t max_intv)
{
    const int32_t n = (int32_t)R.seq.size();
    std::vector<int32_t> mia((size_t)n, max_intv);
    std::vector<SMEM> out((size_t)n * (size_t)R.max_len);
    const int64_t n_out = fmi->bwtSeedStrategyAllPosOneThread_lockstep(
        R.enc.data(), mia.data(), n, R.seq.data(), R.cum.data(), /*minSeedLen=*/19,
        out.data(), R.max_len);
    out.resize((size_t)n_out);
    return out;
}

/* Replace every entry of levels 2..depth with its level's code-0 entry (the
 * interval of a run of As): a valid interval, so a walk that reads it steps on
 * safely, but the wrong one for every other string. */
static void poison_table(FMI_search *fmi)
{
    for (int L = 2; L <= fmi->seed_jump_depth; L++) {
        FMI_search::SeedJumpEntry *lev = fmi->seed_jump + FMI_search::seed_jump_level_off(L);
        for (int64_t c = 1; c < ((int64_t)1 << (2 * L)); c++) lev[c] = lev[0];
    }
}

/* FmiSeed is opaque in fmi_seed_api.h; this is its definition in
 * fmi_seed_api.cpp (identical, as the one-definition rule requires), so the
 * test can reach the handle's FMI_search. */
struct FmiSeed { FMI_search *fmi; };

int main(int argc, char *argv[])
{
    if (argc != 2) {
        fprintf(stderr, "Usage: %s <bwa-mem3 index prefix>\n", argv[0]);
        return 2;
    }
    const char *prefix = argv[1];
    char name[160];

    /* Default load: depth 10, built by load_index. */
    set_jump_env(NULL);
    FMI_search *fmi = new FMI_search(prefix);
    fmi->load_index();
    report(check_table(fmi, 10), "load_index builds the default depth-10 table");

    /* Rebuilds at the depths around the serial prefix (3) and the default,
     * at several thread counts, including more threads than the 64 prefixes
     * the parallel levels are split by (so the builder clamps the count). */
    const struct { int depth; int threads; } builds[] = {
        {2, 1}, {2, 4}, {3, 1}, {3, 4}, {4, 1}, {4, 3}, {5, 64}, {5, 100}, {10, 1}, {10, 7},
    };
    for (const auto &b : builds) {
        char v[8];
        snprintf(v, sizeof v, "%d", b.depth);
        set_jump_env(v);
        fmi->build_seed_jump_table(b.threads);
        snprintf(name, sizeof name, "depth %d, %d builder thread(s): every entry equals the stepped walk",
                 b.depth, b.threads);
        report(check_table(fmi, b.depth), name);
    }

    /* Depths 2-4 against the brute-force text oracle (no backwardExt). */
    TextOracle text;
    const bool have_text = load_text(prefix, &text);
    report(have_text, "index text and naive suffix array loaded from the .pac");
    for (int D = 2; D <= 4; D++) {
        char v[8];
        snprintf(v, sizeof v, "%d", D);
        set_jump_env(v);
        fmi->build_seed_jump_table(2);
        snprintf(name, sizeof name, "depth %d: every entry's s, k and l equal the index text's", D);
        report(have_text && check_table_text(fmi, text, D), name);
    }

    /* The lockstep drivers read the table at the default depth: replacing its
     * entries with other valid intervals must change their output. Without the
     * table both drivers step every walk, and the output is the same either
     * way. */
    if (have_text) {
        set_jump_env(NULL);
        Reads R = make_reads(text);
        fmi->build_seed_jump_table(2);
        const std::vector<SMEM> smem_ok = run_smem_lockstep(fmi, R, 1);
        const std::vector<SMEM> seed_ok = run_bwtseed_lockstep(fmi, R, 1000);
        poison_table(fmi);
        const std::vector<SMEM> smem_bad = run_smem_lockstep(fmi, R, 1);
        const std::vector<SMEM> seed_bad = run_bwtseed_lockstep(fmi, R, 1000);
        report(fmi->seed_jump_depth == 10 && !smem_ok.empty() && !same_smems(smem_ok, smem_bad),
               "SMEM lockstep driver takes its walks' first steps from the depth-10 table");
        report(fmi->seed_jump_depth == 10 && !seed_ok.empty() && !same_smems(seed_ok, seed_bad),
               "bwtSeed lockstep driver takes its walks' first steps from the depth-10 table");

        set_jump_env("0");
        fmi->build_seed_jump_table(2);
        report(fmi->seed_jump == NULL && same_smems(run_smem_lockstep(fmi, R, 1), smem_ok)
                   && same_smems(run_bwtseed_lockstep(fmi, R, 1000), seed_ok),
               "BWA3_BWTSEED_JUMP=0: both lockstep drivers' output equals the table's");
    } else {
        report(false, "lockstep drivers read the table: no index text");
    }

    /* BWA3_BWTSEED_JUMP=0: no table. */
    set_jump_env("0");
    fmi->build_seed_jump_table(4);
    report(fmi->seed_jump == NULL && fmi->seed_jump_depth == 0, "BWA3_BWTSEED_JUMP=0: no table");

    /* Empty, malformed and out-of-range values: the default depth. */
    const char *defaults[] = {"", "1", "13", "-1", "-0", "+5", " 5", "5 ", "10x", "abc",
                              "99999999999999999999"};
    for (const char *v : defaults) {
        set_jump_env(v);
        fmi->build_seed_jump_table(2);
        snprintf(name, sizeof name, "BWA3_BWTSEED_JUMP=\"%s\": default depth 10", v);
        report(fmi->seed_jump != NULL && fmi->seed_jump_depth == 10, name);
    }
    set_jump_env(NULL);

    /* The build_seed_jump=false opt-out, on load_index and on the facade that
     * passes it (fmi_seed_open). */
    {
        FMI_search *np = new FMI_search(prefix);
        np->load_index(/*load_pac=*/false, /*n_threads=*/1, /*build_seed_jump=*/false);
        report(np->seed_jump == NULL && np->seed_jump_depth == 0,
               "load_index(build_seed_jump=false): no table");
        delete np;
        FmiSeed *h = fmi_seed_open(prefix);
        report(h != NULL && h->fmi->shm_base == NULL && h->fmi->seed_jump == NULL
                   && h->fmi->seed_jump_depth == 0,
               "fmi_seed_open (disk load): no table");
        fmi_seed_close(h);
    }

    /* An index attached from a packed segment: load_index_from_shm builds no
     * table, and build_seed_jump_table over the attached arrays builds one
     * equal to the disk load's. The packed segment is copied into an anonymous
     * mapping because the attached FMI_search munmaps its base on destruction. */
    {
        uint8_t *buf = NULL;
        uint64_t len = 0;
        bool ok = bwa_shm_pack_from_disk(prefix, &buf, &len) == 0 && buf != NULL;
        void *map = MAP_FAILED;
        if (ok) {
            map = mmap(NULL, (size_t)len, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANON, -1, 0);
            ok = map != MAP_FAILED;
        }
        bool none = false;
        if (ok) {
            memcpy(map, buf, (size_t)len);
            fmi->build_seed_jump_table(1);   // the disk load's default table again
            FMI_search *at = new FMI_search(prefix);
            at->load_index_from_shm((uint8_t *)map, (size_t)len);
            none = at->seed_jump == NULL && at->seed_jump_depth == 0;
            at->build_seed_jump_table(3);
            const size_t bytes = (size_t)FMI_search::seed_jump_level_off(11) * sizeof(FMI_search::SeedJumpEntry);
            ok = check_table(at, 10) && fmi->seed_jump != NULL
                 && memcmp(at->seed_jump, fmi->seed_jump, bytes) == 0;
            delete at;   // munmaps map
        }
        report(none, "load_index_from_shm alone: no table");
        report(ok, "attached index: table equals the disk load's");
        free(buf);
    }

    /* load_index's shm-attach path, which builds the table itself after the
     * attach and carries the opt-out. Staging needs POSIX shm for the phiX
     * segment. bwa_shm_destroy drops every index staged on the host, so the
     * case runs only when no registry exists yet: then every segment it finds
     * at cleanup is the fixture's own, including after a failed stage. A host
     * with an index already staged skips the case rather than drop it. */
    int ctl_fd = shm_open(BWA_SHM_CTL_NAME, O_RDONLY, 0);
    if (ctl_fd >= 0) {
        close(ctl_fd);
        fprintf(stderr, "[SKIP] load_index's shm-attach path: an index is already staged "
                "in shm on this host (%s exists); not dropping it\n", BWA_SHM_CTL_NAME);
    } else if (bwa_shm_stage(prefix) == 0) {
        FMI_search *at = new FMI_search(prefix);
        at->load_index();
        report(at->shm_base != NULL && check_table(at, 10),
               "load_index (shm attach): builds the default table");
        delete at;
        at = new FMI_search(prefix);
        at->load_index(/*load_pac=*/false, /*n_threads=*/1, /*build_seed_jump=*/false);
        report(at->shm_base != NULL && at->seed_jump == NULL && at->seed_jump_depth == 0,
               "load_index(build_seed_jump=false) (shm attach): no table");
        delete at;
        FmiSeed *h = fmi_seed_open(prefix);
        report(h != NULL && h->fmi->shm_base != NULL && h->fmi->seed_jump == NULL
                   && h->fmi->seed_jump_depth == 0,
               "fmi_seed_open (shm attach): no table");
        fmi_seed_close(h);
        bwa_shm_destroy();
    } else {
        bwa_shm_destroy();   /* the registry, if the failed stage made one, is the fixture's */
        fprintf(stderr, "[SKIP] load_index's shm-attach path: could not stage %s on this host\n", prefix);
    }

    delete fmi;
    fprintf(stderr, "%d / %d cases passed\n", passed_cases, total_cases);
    return passed_cases == total_cases ? 0 : 1;
}
