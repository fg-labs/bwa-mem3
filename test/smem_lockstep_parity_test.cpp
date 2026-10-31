/* Parity harness: scalar getSMEMsOnePosOneThread vs getSMEMsOnePosOneThread_lockstep.
 * Both functions are given identical inputs; outputs (matchArray bytes,
 * numTotalSmem, query_pos_array write-back) must be byte-identical.
 *
 * Usage: smem_lockstep_parity_test <reference FASTA that is also the index prefix>
 * The checked-in phiX fixture (test/fixtures/phix.fa) is REQUIRED, not just an
 * example: Cases 13, 13z, 14, 15, 16, 17 and 18 read reads from that FASTA, and
 * assume a single-contig, all-ACGT reference: Cases 13-14 assert emission
 * floors, and Cases 16-18 build the index's text as the FASTA followed by its
 * reverse complement and expect particular reads to take (or not take) the
 * unique-match text compare. A run against a different index can legitimately
 * fail them; Case 17's text oracle also accepts any read base where the
 * reference is not A/C/G/T (the pac holds a random base there), so N runs
 * weaken it. An unreadable (or too short) FASTA fails them rather than
 * skipping, so the real-match and slot-recycle coverage cannot be lost
 * silently.
 */

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <algorithm>
#include <memory>
#include <string>
#include <vector>
#include "FMI_search.h"
#include "seed_jump_test_util.h"
#include "lockstep_parity_test_util.h"
#include "smem_parked.h"
#include "bwa.h"

/* A degenerate N=1 build would invoke getSMEMsOnePosOneThread_lockstep as a
 * single-slot path, so the cases below would trivially pass without ever
 * exercising slot-interleaving. Fail the build early in that case. */
#if SMEM_LOCKSTEP_N <= 1
#error "smem_lockstep_parity_test requires SMEM_LOCKSTEP_N > 1 to exercise the lockstep path"
#endif

/* Field-by-field SMEM equality on the OBSERVABLE fields (rid,m,n,k,s).
 *
 * `l` (the reverse-complement interval start) is deliberately excluded, mirroring
 * smem_dedup_inplace: the pure-backward phase (backwardExt_konly) intentionally
 * leaves `l` frozen at its forward-phase value while the scalar reference path
 * keeps recomputing it, and the unique-match text compare (ls_s1_compare) never
 * steps it, so the two paths legitimately differ on `l` for a seed emitted after
 * a backward walk through an s>1 interval or finished by the text compare. That
 * divergence is inert — no consumer of an emitted seed reads `l` (SA resolution
 * reads k/s, plus info and m for a tagged s == 1 seed; sort/dedup key on
 * (rid,m,n[,k,s])) — and `l` is a deterministic function of (rid,m,n) anyway, so
 * it carries no independent information. Asserting a.l==b.l here would be a
 * false invariant the lockstep driver does not uphold.
 *
 * `info` is excluded too: the lockstep driver records a unique match's SA
 * position there (a tag) and the scalar walker leaves it 0, so the two differ
 * by design. run_case checks that every untagged SMEM has info 0, and Case 17
 * (check_s1_tags) checks the tags themselves. The seed identity the two paths
 * share is (rid,m,n,k,s).
 *
 * memcmp is unsuitable because the SMEM struct has implicit padding between `n`
 * and `k` that upstream leaves uninitialized — padding is not part of the
 * contract, and neither (now) is `l`. */
static bool smem_fields_equal(const SMEM &a, const SMEM &b) {
    return a.rid == b.rid &&
           a.m   == b.m   &&
           a.n   == b.n   &&
           a.k   == b.k   &&
           a.s   == b.s;
}

static int64_t smem_array_first_mismatch(const SMEM *a, const SMEM *b, int64_t n) {
    for (int64_t i = 0; i < n; i++) {
        if (!smem_fields_equal(a[i], b[i])) return i;
    }
    return -1;
}

static int32_t total_cases = 0;
static int32_t passed_cases = 0;

/* A fixed width a case pins, clamped to the build's SMEM_LOCKSTEP_N_MAX so a
 * -D build with a small cap still runs every case instead of tripping the
 * driver's width assert. */
static int32_t clamp_width(int32_t w)
{
    return w < SMEM_LOCKSTEP_N_MAX ? w : SMEM_LOCKSTEP_N_MAX;
}

static bool run_case(FMI_search *fmi,
                     const char *case_name,
                     int32_t numReads,
                     int32_t max_readlength,
                     int32_t minSeedLen,
                     uint8_t *enc_qdb,
                     int32_t *query_pos_array_in,
                     int32_t *min_intv_array_in,
                     int32_t *rid_array_in,
                     const bseq1_t *seq_,
                     int32_t *query_cum_len_ar,
                     int64_t expect_min_smem = 0,
                     const bool *expect_smem_for_read = NULL)
{
    // Tight upper bound: the scalar path emits at most O(readlength) SMEMs per
    // read (one per backward-step decision plus the final prev[0]), so
    // numReads * max_readlength safely caps both scalar_out and lockstep_out.
    const int32_t max_out = numReads * max_readlength;
    SMEM *scalar_out   = (SMEM *)_mm_malloc(max_out * sizeof(SMEM), 64);
    SMEM *lockstep_out = (SMEM *)_mm_malloc(max_out * sizeof(SMEM), 64);

    int32_t *qpa_s = (int32_t *)malloc(numReads * sizeof(int32_t));
    int32_t *qpa_l = (int32_t *)malloc(numReads * sizeof(int32_t));
    int32_t *mia_s = (int32_t *)malloc(numReads * sizeof(int32_t));
    int32_t *mia_l = (int32_t *)malloc(numReads * sizeof(int32_t));
    int32_t *rid_s = (int32_t *)malloc(numReads * sizeof(int32_t));
    int32_t *rid_l = (int32_t *)malloc(numReads * sizeof(int32_t));
    if (!qpa_s || !qpa_l || !mia_s || !mia_l || !rid_s || !rid_l) {
        fprintf(stderr, "[FAIL] %s: allocation failed\n", case_name);
        _mm_free(scalar_out);
        _mm_free(lockstep_out);
        free(qpa_s); free(qpa_l); free(mia_s); free(mia_l); free(rid_s); free(rid_l);
        total_cases++;
        return false;
    }
    memcpy(qpa_s, query_pos_array_in, numReads * sizeof(int32_t));
    memcpy(qpa_l, query_pos_array_in, numReads * sizeof(int32_t));
    memcpy(mia_s, min_intv_array_in, numReads * sizeof(int32_t));
    memcpy(mia_l, min_intv_array_in, numReads * sizeof(int32_t));
    memcpy(rid_s, rid_array_in, numReads * sizeof(int32_t));
    memcpy(rid_l, rid_array_in, numReads * sizeof(int32_t));

    int64_t n_scalar = 0, n_lockstep = 0;

    // The lockstep variant of getSMEMsOnePosOneThread heap-allocates its
    // per-slot prev/match scratch internally (sized from max_readlength),
    // so the harness no longer pre-allocates them.

    fmi->getSMEMsOnePosOneThread(enc_qdb, qpa_s, mia_s, rid_s,
                                 numReads, numReads, seq_, query_cum_len_ar,
                                 max_readlength, minSeedLen, scalar_out,
                                 &n_scalar);
    fmi->getSMEMsOnePosOneThread_lockstep(enc_qdb, qpa_l, mia_l, rid_l,
                                          numReads, numReads, seq_, query_cum_len_ar,
                                          max_readlength, minSeedLen, lockstep_out,
                                          &n_lockstep);

    total_cases++;
    bool ok = true;
    if (n_scalar != n_lockstep) {
        fprintf(stderr, "[FAIL] %s: numTotalSmem mismatch: scalar=%lld lockstep=%lld\n",
                case_name, (long long)n_scalar, (long long)n_lockstep);
        ok = false;
    } else {
        int64_t mism = smem_array_first_mismatch(scalar_out, lockstep_out, n_scalar);
        if (mism >= 0) {
            const SMEM &a = scalar_out[mism];
            const SMEM &b = lockstep_out[mism];
            fprintf(stderr, "[FAIL] %s: SMEM field mismatch at index %lld:\n"
                    "  scalar   rid=%d m=%u n=%u k=%lld l=%lld s=%lld\n"
                    "  lockstep rid=%d m=%u n=%u k=%lld l=%lld s=%lld\n",
                    case_name, (long long)mism,
                    a.rid, a.m, a.n, (long long)a.k, (long long)a.l, (long long)a.s,
                    b.rid, b.m, b.n, (long long)b.k, (long long)b.l, (long long)b.s);
            ok = false;
        } else if (memcmp(qpa_s, qpa_l, numReads * sizeof(int32_t)) != 0) {
            fprintf(stderr, "[FAIL] %s: query_pos_array write-back differs\n", case_name);
            ok = false;
        } else if (n_scalar < expect_min_smem) {
            /* Parity between two empty outputs proves nothing about emission. */
            fprintf(stderr, "[FAIL] %s: emission path not exercised: got %lld SMEMs, "
                    "expected >= %lld\n", case_name, (long long)n_scalar,
                    (long long)expect_min_smem);
            ok = false;
        } else if (expect_smem_for_read != NULL) {
            /* Every flagged read must have emitted at least one SMEM. */
            bool *seen = (bool *)calloc(numReads, sizeof(bool));
            int32_t missing = -1;
            if (seen == NULL) {
                fprintf(stderr, "[FAIL] %s: allocation failed\n", case_name);
                ok = false;
            } else {
                for (int64_t j = 0; j < n_scalar; j++) {
                    const int32_t r = scalar_out[j].rid;
                    if (r >= 0 && r < numReads) seen[r] = true;
                }
                for (int32_t r = 0; r < numReads && missing < 0; r++)
                    if (expect_smem_for_read[r] && !seen[r]) missing = r;
                free(seen);
                if (missing >= 0) {
                    fprintf(stderr, "[FAIL] %s: read %d emitted no SMEM\n",
                            case_name, missing);
                    ok = false;
                }
            }
        }
        /* SA resolution trusts a tag in info (trust_s1_tag), so every SMEM a
         * producer does not tag must carry info 0. The scalar walker never
         * tags; the lockstep driver tags only s == 1 SMEMs, and only when the
         * index's text is loaded (s1_text_available). */
        if (ok) {
            const char *which = NULL;
            int64_t bad = -1;
            for (int64_t j = 0; j < n_scalar && bad < 0; j++)
                if (scalar_out[j].info != 0) { which = "scalar"; bad = j; }
            const bool text = fmi->s1_text_available();
            for (int64_t j = 0; j < n_lockstep && bad < 0; j++) {
                const SMEM &x = lockstep_out[j];
                if (text && x.s == 1 && (x.info >> 63) != 0) continue;
                if (x.info != 0) { which = "lockstep"; bad = j; }
            }
            if (bad >= 0) {
                const SMEM &x = (which[0] == 's') ? scalar_out[bad] : lockstep_out[bad];
                fprintf(stderr, "[FAIL] %s: %s SMEM %lld (rid=%u m=%u n=%u s=%lld) has info "
                        "%#llx, but an SMEM that is not tagged must have info 0%s\n",
                        case_name, which, (long long)bad, x.rid, x.m, x.n, (long long)x.s,
                        (unsigned long long)x.info,
                        text ? "" : " (the index's text is not loaded, so nothing may be tagged)");
                ok = false;
            }
        }
        if (ok) {
            fprintf(stderr, "[PASS] %s  (n_smem=%lld)\n", case_name, (long long)n_scalar);
            passed_cases++;
        }
    }

    _mm_free(scalar_out);
    _mm_free(lockstep_out);
    free(qpa_s); free(qpa_l); free(mia_s); free(mia_l); free(rid_s); free(rid_l);
    return ok;
}

/* The SA tag the unique-match text compare leaves in SMEM.info: bit 63 set,
 * bits 0-62 the sign-extended q with SA[k] = q + m (FMI_search.cpp,
 * ls_s1_sa_tag). Decoded here independently of the library. */
static const uint64_t S1_TAG = 1ULL << 63;
static bool s1_tagged(const SMEM &a) { return (a.info & S1_TAG) != 0; }
static int64_t s1_coord(const SMEM &a) { return ((int64_t)(a.info << 1) >> 1) + (int64_t)a.m; }

/* Every SMEM a producer does not tag must carry info 0, because SA resolution
 * with trust_s1_tag writes a coordinate for any s == 1 SMEM with bit 63 set.
 * allow_tag: the producer may tag s == 1 SMEMs (the lockstep SMEM driver). */
static bool check_info_zero(const char *case_name, const char *producer,
                            const SMEM *a, int64_t n, bool allow_tag)
{
    for (int64_t i = 0; i < n; i++) {
        if (allow_tag && a[i].s == 1 && s1_tagged(a[i])) continue;
        if (a[i].info != 0) {
            fprintf(stderr, "[FAIL] %s: %s SMEM rid=%u m=%u n=%u s=%lld has info %#llx, "
                    "but an SMEM it does not tag must have info 0\n", case_name, producer,
                    a[i].rid, a[i].m, a[i].n, (long long)a[i].s, (unsigned long long)a[i].info);
            return false;
        }
    }
    return true;
}

/* Restores the process-wide SA-dedup mode (ks_dedup_configure) on every exit
 * path, so a failed case does not leave it forced for the cases after it. */
struct KsDedupModeGuard {
    ~KsDedupModeGuard() { ks_dedup_configure(NULL); }
};

/* The SA tags, checked against two independent oracles:
 *   1. the text itself: a tagged SMEM's coordinate must hold read[m..n] in
 *      `text` (the index's text: forward reference + reverse complement);
 *   2. the LF walk: get_sa_entries_prefetch with trust_s1_tag (the tag written
 *      without a walk) must return the same coordinates as with it off (every
 *      position walked), in the plain and the dedup staging paths, at
 *      max_occ 500 (the default) and 2 (so the sampled step, s / max_occ, of
 *      the SMEMs staged around the tagged ones is exercised, on extra
 *      untagged SMEMs with s >= 2 * max_occ), and must ignore the tag under a
 *      dropped sentinel offset (--compat=bwa-mem2). The case fails unless some
 *      tagged SMEM's walk reaches the sentinel row there (its dropped-offset
 *      coordinate differs from its tag), so that leg cannot pass vacuously.
 * The SMEMs are made as mem_collect_smem makes them: round 1 (the lockstep SMEM
 * driver at min_intv 1 from each read's pivot), round 2 (the re-seed pass, from
 * the middle of each round-1 SMEM of at least split_len bases and s <= 10, at
 * min_intv s + 1), and round 3 (the production bwtSeed driver,
 * bwtSeedStrategyAllPosOneThread_lockstep where built, at min_seed_len + 1 and
 * max_intv 20; the scalar driver's output is checked for info 0 too). Every
 * untagged SMEM of every round must have info 0. They are sorted with
 * smem_sort_rid_mn and tagged by containment with s1_tag_contained_reads, the
 * calls mem_collect_smem makes; this does not cover mem_collect_smem's own
 * gating of that call or smem_dedup_inplace, which runs after it.
 * expect_walk_tag[r] (optional): 1 if read r's round-1 walk must take the text
 * compare (leave a tagged SMEM), -1 if it must not, 0 if either is fine.
 * The case fails if no SMEM is tagged by the driver, or none by containment:
 * either would leave the shortcut untested while the coordinate checks pass. */
static bool check_s1_tags(FMI_search *fmi, const char *case_name, int32_t numReads,
                          int32_t max_readlength, uint8_t *enc_qdb,
                          const int32_t *query_pos_in, const bseq1_t *seq_,
                          int32_t *query_cum_len_ar, const std::string &text,
                          const signed char *expect_walk_tag)
{
    total_cases++;
    if (!fmi->s1_text_available()) {
        fprintf(stderr, "[FAIL] %s: the index's pac is not loaded, so the text compare "
                "cannot run\n", case_name);
        return false;
    }
    const int32_t min_seed_len = 19;                         /* -k default */
    const int32_t split_len = (int32_t)(min_seed_len * 1.5 + .499), split_width = 10;
    const size_t cap = (size_t)numReads * (size_t)max_readlength;

    /* Round 1. */
    std::vector<SMEM> r1(cap);
    std::vector<int32_t> qpa(query_pos_in, query_pos_in + numReads);
    std::vector<int32_t> mia((size_t)numReads, 1), rid((size_t)numReads);
    for (int32_t i = 0; i < numReads; i++) rid[(size_t)i] = i;
    int64_t n1 = 0;
    fmi->getSMEMsOnePosOneThread_lockstep(enc_qdb, qpa.data(), mia.data(), rid.data(),
                                          numReads, numReads, seq_, query_cum_len_ar,
                                          max_readlength, min_seed_len, r1.data(), &n1);
    r1.resize((size_t)n1);
    if (!check_info_zero(case_name, "round-1 lockstep", r1.data(), n1, true)) return false;
    if (expect_walk_tag != NULL) {
        std::vector<char> has_tag((size_t)numReads, 0);
        for (const SMEM &a : r1) if (s1_tagged(a)) has_tag[a.rid] = 1;
        for (int32_t r = 0; r < numReads; r++) {
            if (expect_walk_tag[r] == 0 || (expect_walk_tag[r] > 0) == (has_tag[(size_t)r] != 0))
                continue;
            fprintf(stderr, "[FAIL] %s: read %d (pivot %d) %s the text compare, but its "
                    "input was built to %s\n", case_name, (int)r, (int)query_pos_in[r],
                    has_tag[(size_t)r] ? "took" : "did not take",
                    expect_walk_tag[r] > 0 ? "reach it" : "keep stepping");
            return false;
        }
    }

    /* Round 2: re-seed from the middle of long, low-occurrence round-1 SMEMs. */
    std::vector<int32_t> qpa2, mia2, rid2;
    for (const SMEM &a : r1) {
        const int32_t start = (int32_t)a.m, end = (int32_t)a.n + 1;
        if (end - start < split_len || a.s > split_width) continue;
        qpa2.push_back((end + start) >> 1);
        mia2.push_back((int32_t)a.s + 1);
        rid2.push_back((int32_t)a.rid);
    }
    const int32_t pos2 = (int32_t)qpa2.size();
    std::vector<SMEM> r2((size_t)pos2 * (size_t)max_readlength + 1);
    int64_t n2 = 0;
    if (pos2 > 0)
        fmi->getSMEMsOnePosOneThread_lockstep(enc_qdb, qpa2.data(), mia2.data(), rid2.data(),
                                              pos2, pos2, seq_, query_cum_len_ar,
                                              max_readlength, min_seed_len, r2.data(), &n2);
    r2.resize((size_t)n2);
    if (!check_info_zero(case_name, "round-2 lockstep", r2.data(), n2, true)) return false;

    /* Round 3: both bwtSeed drivers; the production one feeds the oracles. */
    std::vector<int32_t> max_intv((size_t)numReads, 20);
    std::vector<SMEM> r3s(cap);
    const int64_t n3s = fmi->bwtSeedStrategyAllPosOneThread(enc_qdb, max_intv.data(), numReads,
                                                            seq_, query_cum_len_ar,
                                                            min_seed_len + 1, r3s.data());
    r3s.resize((size_t)n3s);
    if (!check_info_zero(case_name, "round-3 scalar bwtSeed", r3s.data(), n3s, false))
        return false;
#if BWTSEED_LOCKSTEP_N > 1
    std::vector<SMEM> r3(cap);
    const int64_t n3 = fmi->bwtSeedStrategyAllPosOneThread_lockstep(
        enc_qdb, max_intv.data(), numReads, seq_, query_cum_len_ar, min_seed_len + 1,
        r3.data(), max_readlength);
    r3.resize((size_t)n3);
    if (!check_info_zero(case_name, "round-3 lockstep bwtSeed", r3.data(), n3, false))
        return false;
#else
    const std::vector<SMEM> &r3 = r3s;
#endif

    std::vector<SMEM> smems;
    smems.reserve(r1.size() + r2.size() + r3.size());
    smems.insert(smems.end(), r1.begin(), r1.end());
    smems.insert(smems.end(), r2.begin(), r2.end());
    smems.insert(smems.end(), r3.begin(), r3.end());
    const int64_t n = (int64_t)smems.size();
    int64_t n_tag_walk = 0, n_tag_contained = 0;
    for (const SMEM &a : smems) n_tag_walk += s1_tagged(a);
    SmemSortScratch scratch;
    memset(&scratch, 0, sizeof scratch);
    smem_sort_rid_mn(smems.data(), n, scratch);
    _mm_free(scratch.cnt);
    _mm_free(scratch.tmp);
    _mm_free(scratch.pairs);
    FMI_search::s1_tag_contained_reads(smems.data(), n);
    for (const SMEM &a : smems) n_tag_contained += s1_tagged(a);
    n_tag_contained -= n_tag_walk;

    /* Oracle 1: the text holds read[m..n] at the tagged coordinate. */
    for (const SMEM &a : smems) {
        if (!s1_tagged(a)) continue;
        const int64_t coord = s1_coord(a);
        const int64_t len = (int64_t)a.n - (int64_t)a.m + 1;
        bool ok = a.s == 1 && coord >= 0 && coord + len <= (int64_t)text.size();
        const uint8_t *q = enc_qdb + query_cum_len_ar[a.rid];
        for (int64_t i = 0; ok && i < len; i++) {
            const char t = text[(size_t)(coord + i)];
            ok = (t != 'A' && t != 'C' && t != 'G' && t != 'T') ||   /* a hole: pac's base is random */
                 "ACGTN"[q[a.m + i]] == t;
        }
        if (!ok) {
            fprintf(stderr, "[FAIL] %s: tagged SMEM rid=%u m=%u n=%u s=%lld: coordinate %lld "
                    "does not hold the read\n", case_name, a.rid, a.m, a.n, (long long)a.s,
                    (long long)coord);
            return false;
        }
    }

    /* Oracle 2: the LF walk, per max_occ, staging path and sentinel mode. */
    KsDedupModeGuard dedup_guard;
    int64_t n_sentinel = 0;   /* tagged SMEMs whose dropped-offset walk differs from the tag */
    for (const int32_t max_occ : {500, 2}) {
        std::vector<SMEM> arr = smems;
        if (max_occ < 500) {
            /* Untagged SMEMs whose s / max_occ step is > 1 (rows of the BWT;
             * nothing here needs them to be some read's SMEMs). */
            const int64_t rows = fmi->reference_seq_len;
            const int64_t wide[][2] = {{rows / 7, 5}, {rows / 3, 9}, {rows / 2 + 1, 37}};
            for (const auto &w : wide) {
                SMEM x;
                memset(&x, 0, sizeof x);
                x.rid = (uint32_t)numReads;
                x.k = w[0];
                x.s = w[1];
                arr.push_back(x);
            }
        }
        const int64_t na = (int64_t)arr.size();
        std::vector<int64_t> slot((size_t)na);
        int64_t total = 0;
        for (int64_t i = 0; i < na; i++) {
            slot[(size_t)i] = total;
            const int64_t s = arr[(size_t)i].s;
            total += s > max_occ ? max_occ : (s < 0 ? 0 : s);
        }
        for (const char *mode : {"off", "on"}) {
            ks_dedup_configure(mode);
            for (int drop : {0, 1}) {
                std::vector<int64_t> walk((size_t)total + 1, -1), tag((size_t)total + 1, -2);
                int64_t cw = 0, ct = 0, iw = 0, it = 0;
                fmi->get_sa_entries_prefetch(arr.data(), walk.data(), &cw, na, max_occ, 0, iw,
                                             drop, false);
                fmi->get_sa_entries_prefetch(arr.data(), tag.data(), &ct, na, max_occ, 0, it,
                                             drop, true);
                if (cw != total || ct != total || iw != it) {
                    fprintf(stderr, "[FAIL] %s: max_occ %d, dedup %s, drop %d: coordinate "
                            "counts differ (expected %lld, walk %lld, tag %lld)\n", case_name,
                            (int)max_occ, mode, drop, (long long)total, (long long)cw,
                            (long long)ct);
                    return false;
                }
                for (int64_t i = 0; i < total; i++) {
                    if (walk[(size_t)i] != tag[(size_t)i]) {
                        fprintf(stderr, "[FAIL] %s: max_occ %d, dedup %s, drop %d: coordinate "
                                "%lld is %lld with the tag, %lld by the LF walk\n", case_name,
                                (int)max_occ, mode, drop, (long long)i,
                                (long long)tag[(size_t)i], (long long)walk[(size_t)i]);
                        return false;
                    }
                }
                if (drop)
                    for (int64_t i = 0; i < na; i++)
                        if (arr[(size_t)i].s == 1 && s1_tagged(arr[(size_t)i]) &&
                            walk[(size_t)slot[(size_t)i]] != s1_coord(arr[(size_t)i]))
                            n_sentinel++;
            }
        }
    }

    if (n_tag_walk == 0 || n_tag_contained == 0 || n_sentinel == 0) {
        fprintf(stderr, "[FAIL] %s: shortcut not exercised: %lld SMEMs tagged by the walk, "
                "%lld by containment, %lld whose walk meets the sentinel row under a dropped "
                "offset (each must be > 0)\n", case_name, (long long)n_tag_walk,
                (long long)n_tag_contained, (long long)n_sentinel);
        return false;
    }
    fprintf(stderr, "[PASS] %s  (n_smem=%lld, tagged: %lld by the walk, %lld by containment; "
            "%lld sentinel-row walks)\n", case_name, (long long)n, (long long)n_tag_walk,
            (long long)n_tag_contained, (long long)n_sentinel);
    passed_cases++;
    return true;
}

/* s1_tag_contained's rules on hand-built SMEMs, through s1_tag_contained_reads
 * (the per-read split mem_collect_smem uses). The tags are arbitrary (no index):
 * each element's info is checked exactly, so a rule that tags an SMEM not
 * inside a tagged one, misses one that is, or crosses reads, fails here even
 * where a real read's single locus would hide it. */
static bool check_s1_tag_contained_rules()
{
    total_cases++;
    const char *name = "Case 19: s1_tag_contained on hand-built SMEMs";
    const uint64_t TA = S1_TAG | 1000, TB = S1_TAG | 5000;
    struct Row { uint32_t rid, m, n; int64_t s; uint64_t info, want; const char *why; };
    const Row rows[] = {
        {0, 0, 10, 1, 0,  0,  "an earlier m run (m < m_T) is not inside T"},
        {0, 2, 15, 1, 0,  TA, "T's m run, before T, inside it: back-filled"},
        {0, 2, 15, 2, 0,  0,  "inside T but s == 2"},
        {0, 2, 30, 1, TA, TA, "T"},
        {0, 5, 20, 1, 0,  TA, "inside T"},
        {0, 5, 30, 1, 0,  TA, "inside T, same end"},
        {0, 6, 20, 3, 0,  0,  "inside T but s == 3"},
        {0, 7, 31, 1, 0,  0,  "ends one base past T"},
        {0, 9, 40, 1, TB, TB, "a second tagged SMEM keeps its own tag"},
        {0, 33, 40, 1, 0, TB, "inside only the second"},
        {0, 35, 41, 1, 0, 0,  "ends past both"},
        {1, 3, 8, 1, 0,   0,  "another read: never takes read 0's tags"},
        {1, 4, 50, 2, 0,  0,  "another read, s == 2"},
    };
    const int64_t n = (int64_t)(sizeof rows / sizeof rows[0]);
    std::vector<SMEM> a((size_t)n);
    for (int64_t i = 0; i < n; i++) {
        memset(&a[(size_t)i], 0, sizeof(SMEM));
        a[(size_t)i].rid = rows[i].rid;
        a[(size_t)i].m = rows[i].m;
        a[(size_t)i].n = rows[i].n;
        a[(size_t)i].s = rows[i].s;
        a[(size_t)i].info = rows[i].info;
    }
    FMI_search::s1_tag_contained_reads(a.data(), n);
    for (int64_t i = 0; i < n; i++) {
        if (a[(size_t)i].info != rows[i].want) {
            fprintf(stderr, "[FAIL] %s: SMEM %lld (rid=%u m=%u n=%u s=%lld, %s) has info %#llx, "
                    "expected %#llx\n", name, (long long)i, rows[i].rid, rows[i].m, rows[i].n,
                    (long long)rows[i].s, rows[i].why, (unsigned long long)a[(size_t)i].info,
                    (unsigned long long)rows[i].want);
            return false;
        }
    }
    fprintf(stderr, "[PASS] %s  (%lld SMEMs)\n", name, (long long)n);
    passed_cases++;
    return true;
}

int main(int argc, char *argv[]) {
    if (argc != 2) {
        fprintf(stderr, "Usage: %s <bwa-mem3 index prefix>\n", argv[0]);
        return 2;
    }

    fprintf(stderr, "smem_lockstep_parity_test: SMEM_LOCKSTEP_N=%d\n",
            (int)SMEM_LOCKSTEP_N);

    FMI_search *fmi = new FMI_search(argv[1]);
    fmi->load_index();

    /* Case 1: two simple reads, start at position 0, minSeedLen 19. */
    {
        const char *reads[] = {
            "ACGTACGTACGTACGTACGTACGTACGTACGT",
            "TGCATGCATGCATGCATGCATGCATGCATGCA",
        };
        int32_t numReads = 2;
        int32_t max_readlength = 32;
        uint8_t *enc_qdb; bseq1_t *seq_; int32_t *cum_len;
        encode_reads(reads, numReads, &enc_qdb, &seq_, &cum_len);
        int32_t qpa[] = {0, 0};
        int32_t mia[] = {1, 1};
        int32_t rid[] = {0, 1};
        run_case(fmi, "Case 1: two simple reads",
                 numReads, max_readlength, 19,
                 enc_qdb, qpa, mia, rid, seq_, cum_len);
        free(enc_qdb); free(seq_); free(cum_len);
    }

    /* Cases 2-4 and 6 pin the width to N=4 (clamped to SMEM_LOCKSTEP_N_MAX), so
     * their read counts sit below, at, one above and past it as described,
     * whatever the compile-time default. Cases 2-4 are defined by where their
     * read count sits relative to N=4, so a build capped below 4 skips them
     * (reported as SKIP) rather than running them at a width they don't describe.
     *
     * Case 2: numReads < N (N=4, 2 reads). Exercises the
     * "unused slots" path where the driver must mark slot indices >= numReads
     * as DONE/invalid at entry so the stepping pass leaves them idle. */
    {
        const char *reads[] = {
            "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT",
            "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT",
        };
        int32_t numReads = 2;
        int32_t max_readlength = 40;
        uint8_t *enc_qdb; bseq1_t *seq_; int32_t *cum_len;
        encode_reads(reads, numReads, &enc_qdb, &seq_, &cum_len);
        int32_t qpa[] = {0, 0};
        int32_t mia[] = {1, 1};
        int32_t rid[] = {0, 1};
        const int32_t saved_n = g_smem_lockstep_n;
        g_smem_lockstep_n = clamp_width(4);
        if (g_smem_lockstep_n == 4)
            run_case(fmi, "Case 2: numReads < N",
                     numReads, max_readlength, 19,
                     enc_qdb, qpa, mia, rid, seq_, cum_len);
        else
            fprintf(stderr, "[SKIP] Case 2: needs N=4, but SMEM_LOCKSTEP_N_MAX is %d\n",
                    (int)SMEM_LOCKSTEP_N_MAX);
        g_smem_lockstep_n = saved_n;
        free(enc_qdb); free(seq_); free(cum_len);
    }

    /* Case 3: numReads == N (single full batch, no slot recycle).
     * Exercises the steady state where every slot is occupied and no
     * slot is refilled. */
    {
        const char *reads[] = {
            "ACGTACGTACGTACGTACGTACGTACGTACGT",
            "TGCATGCATGCATGCATGCATGCATGCATGCA",
            "AAAACCCCGGGGTTTTAAAACCCCGGGGTTTT",
            "GATCGATCGATCGATCGATCGATCGATCGATC",
        };
        int32_t numReads = 4;
        int32_t max_readlength = 32;
        uint8_t *enc_qdb; bseq1_t *seq_; int32_t *cum_len;
        encode_reads(reads, numReads, &enc_qdb, &seq_, &cum_len);
        int32_t qpa[] = {0, 0, 0, 0};
        int32_t mia[] = {1, 1, 1, 1};
        int32_t rid[] = {0, 1, 2, 3};
        const int32_t saved_n = g_smem_lockstep_n;
        g_smem_lockstep_n = clamp_width(4);
        if (g_smem_lockstep_n == 4)
            run_case(fmi, "Case 3: numReads == N",
                     numReads, max_readlength, 19,
                     enc_qdb, qpa, mia, rid, seq_, cum_len);
        else
            fprintf(stderr, "[SKIP] Case 3: needs N=4, but SMEM_LOCKSTEP_N_MAX is %d\n",
                    (int)SMEM_LOCKSTEP_N_MAX);
        g_smem_lockstep_n = saved_n;
        free(enc_qdb); free(seq_); free(cum_len);
    }

    /* Case 4: numReads == N+1 (exactly one slot recycle).
     * Exercises the retire-then-refill path. */
    {
        const char *reads[] = {
            "ACGTACGTACGTACGTACGTACGTACGTACGT",
            "TGCATGCATGCATGCATGCATGCATGCATGCA",
            "AAAACCCCGGGGTTTTAAAACCCCGGGGTTTT",
            "GATCGATCGATCGATCGATCGATCGATCGATC",
            "CAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGT",
        };
        int32_t numReads = 5;
        int32_t max_readlength = 32;
        uint8_t *enc_qdb; bseq1_t *seq_; int32_t *cum_len;
        encode_reads(reads, numReads, &enc_qdb, &seq_, &cum_len);
        int32_t qpa[] = {0, 0, 0, 0, 0};
        int32_t mia[] = {1, 1, 1, 1, 1};
        int32_t rid[] = {0, 1, 2, 3, 4};
        const int32_t saved_n = g_smem_lockstep_n;
        g_smem_lockstep_n = clamp_width(4);
        if (g_smem_lockstep_n == 4)
            run_case(fmi, "Case 4: numReads == N+1",
                     numReads, max_readlength, 19,
                     enc_qdb, qpa, mia, rid, seq_, cum_len);
        else
            fprintf(stderr, "[SKIP] Case 4: needs N=4, but SMEM_LOCKSTEP_N_MAX is %d\n",
                    (int)SMEM_LOCKSTEP_N_MAX);
        g_smem_lockstep_n = saved_n;
        free(enc_qdb); free(seq_); free(cum_len);
    }

    /* Case 5: first base non-ACGT (N character). The scalar path skips the
     * whole read when enc_qdb[offset + x] >= 4; lockstep must match by
     * going straight to PH_DONE with zero matches, without writing anything
     * to matchArray for this read. Mixed with valid reads to ensure the
     * zero-match case doesn't corrupt the emit cursor. */
    {
        const char *reads[] = {
            "NACGTACGTACGTACGTACGTACGTACGTACG",
            "TGCATGCATGCATGCATGCATGCATGCATGCA",
            "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN",
            "GATCGATCGATCGATCGATCGATCGATCGATC",
        };
        int32_t numReads = 4;
        int32_t max_readlength = 32;
        uint8_t *enc_qdb; bseq1_t *seq_; int32_t *cum_len;
        encode_reads(reads, numReads, &enc_qdb, &seq_, &cum_len);
        int32_t qpa[] = {0, 0, 0, 0};
        int32_t mia[] = {1, 1, 1, 1};
        int32_t rid[] = {0, 1, 2, 3};
        run_case(fmi, "Case 5: first base non-ACGT",
                 numReads, max_readlength, 19,
                 enc_qdb, qpa, mia, rid, seq_, cum_len);
        free(enc_qdb); free(seq_); free(cum_len);
    }

    /* Case 6: numReads not divisible by N (tail partial batch).
     * At N=4, after filling 4 slots, 3 more reads get pulled in as slots retire.
     * A build capped below 4 runs it at the cap (2 or 3), where 7 reads still
     * leave a partial tail batch.
     * When only fewer-than-N reads remain and the last ones finish, the
     * driver must correctly mark the slots retired (input_idx = -1) without
     * trying to recycle from an exhausted queue. */
    {
        const char *reads[] = {
            "ACGTACGTACGTACGTACGTACGTACGTACGT",
            "TGCATGCATGCATGCATGCATGCATGCATGCA",
            "AAAACCCCGGGGTTTTAAAACCCCGGGGTTTT",
            "GATCGATCGATCGATCGATCGATCGATCGATC",
            "CAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGT",
            "TTTACCCAGGGCTTTACCCAGGGCTTTACCCA",
            "ATATATATATATATATATATATATATATATAT",
        };
        int32_t numReads = 7;  /* not a multiple of N at N=4, 3 or 2 */
        int32_t max_readlength = 32;
        uint8_t *enc_qdb; bseq1_t *seq_; int32_t *cum_len;
        encode_reads(reads, numReads, &enc_qdb, &seq_, &cum_len);
        int32_t qpa[] = {0, 0, 0, 0, 0, 0, 0};
        int32_t mia[] = {1, 1, 1, 1, 1, 1, 1};
        int32_t rid[] = {0, 1, 2, 3, 4, 5, 6};
        const int32_t saved_n = g_smem_lockstep_n;
        g_smem_lockstep_n = clamp_width(4);
        run_case(fmi, "Case 6: numReads not divisible by N",
                 numReads, max_readlength, 19,
                 enc_qdb, qpa, mia, rid, seq_, cum_len);
        g_smem_lockstep_n = saved_n;
        free(enc_qdb); free(seq_); free(cum_len);
    }

    /* Case 7: mixed start_pos (some reads start mid-read rather than at 0).
     * Exercises the forward-ext `j = start_pos + 1` initialization and the
     * backward-ext `j = start_pos - 1` initialization for non-zero starts.
     * The scalar path uses `query_pos_array[i]` directly; lockstep carries
     * it in slot.start_pos, which must round-trip correctly. */
    {
        const char *reads[] = {
            "ACGTACGTACGTACGTACGTACGTACGTACGT",
            "TGCATGCATGCATGCATGCATGCATGCATGCA",
            "AAAACCCCGGGGTTTTAAAACCCCGGGGTTTT",
            "GATCGATCGATCGATCGATCGATCGATCGATC",
        };
        int32_t numReads = 4;
        int32_t max_readlength = 32;
        uint8_t *enc_qdb; bseq1_t *seq_; int32_t *cum_len;
        encode_reads(reads, numReads, &enc_qdb, &seq_, &cum_len);
        int32_t qpa[] = {5, 10, 0, 15};   /* mixed start positions */
        int32_t mia[] = {1, 1, 1, 1};
        int32_t rid[] = {0, 1, 2, 3};
        run_case(fmi, "Case 7: mixed start_pos",
                 numReads, max_readlength, 19,
                 enc_qdb, qpa, mia, rid, seq_, cum_len);
        free(enc_qdb); free(seq_); free(cum_len);
    }

    /* Case 8: real 300bp Illumina WGS read pair from SRR6109255
     * (SRR6109255.100035 R1/R2). Exercises long SMEM walks through adapter
     * contamination (AGATCGGAAGAG around base ~22-24) and low-complexity
     * A-run tails — patterns the synthetic Cases 1-7 don't hit. Regression
     * coverage for real-world read composition. */
    {
        /* SRR6109255.100035 R1 (300bp) */
        const char *r1 =
            "ATGACCTCCCTAATATCTTCAGAATAGATCGGAAGAGCACACGTCTGAACTCCAGTCACT"
            "ACAAGATCTCGTATGCCGTCTTCTGCTTGAAAAAAAAAAAAAAAAAATGAGTAACTCTAT"
            "GAATCGAACTCAAACATTCCACACACGTTTATGCCACTACTTCTATCGTGCTTTCATATT"
            "ATACTACTCCCCCTTCCCCCTCCATCTTCACCTCTCCTCCAATCTCCACTCACCTCTACC"
            "TCAACACGTACTTTACACACATCTCTCCCCCACGGACCACAATACCTCTCCTCTCAATTA";
        /* SRR6109255.100035 R2 (300bp) */
        const char *r2 =
            "ATTCTGAAGATATTAGGGAGGTCATAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAG"
            "ATCTCGGTGGTCGCCGTATCATTAAAAAAAAAAAAAATAAAAACAGAAGCAGTGTCAGTA"
            "GCATAGTGATGAATATAGCAATAAACGCACAGCTTGAAACCTACCGTTTGCGACAGCATC"
            "TCTACAGACGTCTGTTTATCTACTTTGAAAGTGGCTACGTGGAACACTCATAGTCATACC"
            "ACTAATCAATATCATGAAATTACCAGGTTAGTCTGTATACTACGATAAAGACACAGACTT";
        const char *reads[] = { r1, r2 };
        int32_t numReads = 2;
        int32_t max_readlength = 300;
        uint8_t *enc_qdb; bseq1_t *seq_; int32_t *cum_len;
        encode_reads(reads, numReads, &enc_qdb, &seq_, &cum_len);
        int32_t qpa[] = {0, 0};
        int32_t mia[] = {1, 1};
        int32_t rid[] = {0, 1};
        run_case(fmi, "Case 8: real 300bp SRR6109255.100035",
                 numReads, max_readlength, 19,
                 enc_qdb, qpa, mia, rid, seq_, cum_len);
        free(enc_qdb); free(seq_); free(cum_len);
    }

    /* Case 9: same real 300bp pair, driven through the outer
     * `getSMEMsAllPosOneThread` do/while loop. The outer driver compacts
     * active reads and re-calls the inner function with updated
     * query_pos_array[i]; a bug in the inner function's write-back of
     * next_x (slot.next_x via ls_advance_forward_step / ls_init_slot) would
     * show up here even if Case 8 passes. We run the scalar and lockstep
     * outer drivers locally rather than via getSMEMsAllPosOneThread, since
     * that one dispatches based on the SMEM_LOCKSTEP_N macro at compile
     * time — both sides would take the same branch. */
    {
        const char *r1 =
            "ATGACCTCCCTAATATCTTCAGAATAGATCGGAAGAGCACACGTCTGAACTCCAGTCACT"
            "ACAAGATCTCGTATGCCGTCTTCTGCTTGAAAAAAAAAAAAAAAAAATGAGTAACTCTAT"
            "GAATCGAACTCAAACATTCCACACACGTTTATGCCACTACTTCTATCGTGCTTTCATATT"
            "ATACTACTCCCCCTTCCCCCTCCATCTTCACCTCTCCTCCAATCTCCACTCACCTCTACC"
            "TCAACACGTACTTTACACACATCTCTCCCCCACGGACCACAATACCTCTCCTCTCAATTA";
        const char *r2 =
            "ATTCTGAAGATATTAGGGAGGTCATAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAG"
            "ATCTCGGTGGTCGCCGTATCATTAAAAAAAAAAAAAATAAAAACAGAAGCAGTGTCAGTA"
            "GCATAGTGATGAATATAGCAATAAACGCACAGCTTGAAACCTACCGTTTGCGACAGCATC"
            "TCTACAGACGTCTGTTTATCTACTTTGAAAGTGGCTACGTGGAACACTCATAGTCATACC"
            "ACTAATCAATATCATGAAATTACCAGGTTAGTCTGTATACTACGATAAAGACACAGACTT";
        const char *reads[] = { r1, r2 };
        int32_t numReads = 2;
        int32_t max_readlength = 300;
        int32_t minSeedLen = 19;
        uint8_t *enc_qdb; bseq1_t *seq_; int32_t *cum_len;
        encode_reads(reads, numReads, &enc_qdb, &seq_, &cum_len);
        const int32_t max_out = numReads * max_readlength;
        SMEM *scalar_out   = (SMEM *)_mm_malloc(max_out * sizeof(SMEM), 64);
        SMEM *lockstep_out = (SMEM *)_mm_malloc(max_out * sizeof(SMEM), 64);
        /* The lockstep variant now owns its per-slot prev/match scratch
         * internally (sized by MAX_READ_LEN_FOR_LOCKSTEP); no external
         * lockstep buffers are needed here. */

        /* The two do/while blocks below intentionally mirror the compaction
         * logic in FMI_search::getSMEMsAllPosOneThread (see src/FMI_search.cpp),
         * in particular the `query_pos_array[head] < seq_[rid].l_seq` predicate
         * that keeps reads active. If that predicate or the surrounding loop
         * structure changes in production, update both drivers here to match. */
        /* Scalar outer driver. */
        int64_t n_s = 0;
        int32_t qpa_s[2] = {0, 0};
        {
            int32_t rid_work[2] = {0, 1};
            int32_t mia_work[2] = {1, 1};
            int32_t numActive = numReads;
            do {
                int32_t head = 0, tail = 0;
                for (head = 0; head < numActive; head++) {
                    int rl = seq_[rid_work[head]].l_seq;
                    if (qpa_s[head] < rl) {
                        rid_work[tail] = rid_work[head];
                        qpa_s[tail] = qpa_s[head];
                        mia_work[tail] = mia_work[head];
                        tail++;
                    }
                }
                fmi->getSMEMsOnePosOneThread(enc_qdb, qpa_s, mia_work, rid_work,
                                             tail, tail, seq_, cum_len,
                                             max_readlength, minSeedLen, scalar_out,
                                             &n_s);
                numActive = tail;
            } while (numActive > 0);
        }

        /* Lockstep outer driver. */
        int64_t n_l = 0;
        int32_t qpa_l[2] = {0, 0};
        {
            int32_t rid_work[2] = {0, 1};
            int32_t mia_work[2] = {1, 1};
            int32_t numActive = numReads;
            do {
                int32_t head = 0, tail = 0;
                for (head = 0; head < numActive; head++) {
                    int rl = seq_[rid_work[head]].l_seq;
                    if (qpa_l[head] < rl) {
                        rid_work[tail] = rid_work[head];
                        qpa_l[tail] = qpa_l[head];
                        mia_work[tail] = mia_work[head];
                        tail++;
                    }
                }
                fmi->getSMEMsOnePosOneThread_lockstep(enc_qdb, qpa_l, mia_work, rid_work,
                                                     tail, tail, seq_, cum_len,
                                                     max_readlength, minSeedLen, lockstep_out,
                                                     &n_l);
                numActive = tail;
            } while (numActive > 0);
        }

        total_cases++;
        if (n_s != n_l) {
            fprintf(stderr, "[FAIL] Case 9 AllPos driver: numTotalSmem mismatch: scalar=%lld lockstep=%lld\n",
                    (long long)n_s, (long long)n_l);
        } else {
            int64_t mism = smem_array_first_mismatch(scalar_out, lockstep_out, n_s);
            if (mism >= 0) {
                const SMEM &a = scalar_out[mism];
                const SMEM &b = lockstep_out[mism];
                fprintf(stderr, "[FAIL] Case 9 AllPos driver: SMEM mismatch at %lld:\n"
                        "  scalar   rid=%d m=%u n=%u k=%lld l=%lld s=%lld\n"
                        "  lockstep rid=%d m=%u n=%u k=%lld l=%lld s=%lld\n",
                        (long long)mism,
                        a.rid, a.m, a.n, (long long)a.k, (long long)a.l, (long long)a.s,
                        b.rid, b.m, b.n, (long long)b.k, (long long)b.l, (long long)b.s);
            } else {
                fprintf(stderr, "[PASS] Case 9 AllPos driver  (n_smem=%lld)\n", (long long)n_s);
                passed_cases++;
            }
        }
        _mm_free(scalar_out); _mm_free(lockstep_out);
        free(enc_qdb); free(seq_); free(cum_len);
    }

    /* Case 10: exercises the forward-ext `newSmem.s < min_intv` exit path
     * by using the real 300bp reads with a large min_intv (mimicking what
     * mem_collect_smem uses for phase-2 re-seeding: min_intv = prev_smem.s
     * + 1). Synthetic Cases 1-7 use min_intv=1, which rarely triggers this
     * exit path. This case guards the post-loop "push smem if smem.s >=
     * min_intv" behavior on the min_intv branch. */
    {
        const char *r1 =
            "ATGACCTCCCTAATATCTTCAGAATAGATCGGAAGAGCACACGTCTGAACTCCAGTCACT"
            "ACAAGATCTCGTATGCCGTCTTCTGCTTGAAAAAAAAAAAAAAAAAATGAGTAACTCTAT"
            "GAATCGAACTCAAACATTCCACACACGTTTATGCCACTACTTCTATCGTGCTTTCATATT"
            "ATACTACTCCCCCTTCCCCCTCCATCTTCACCTCTCCTCCAATCTCCACTCACCTCTACC"
            "TCAACACGTACTTTACACACATCTCTCCCCCACGGACCACAATACCTCTCCTCTCAATTA";
        const char *r2 =
            "ATTCTGAAGATATTAGGGAGGTCATAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAG"
            "ATCTCGGTGGTCGCCGTATCATTAAAAAAAAAAAAAATAAAAACAGAAGCAGTGTCAGTA"
            "GCATAGTGATGAATATAGCAATAAACGCACAGCTTGAAACCTACCGTTTGCGACAGCATC"
            "TCTACAGACGTCTGTTTATCTACTTTGAAAGTGGCTACGTGGAACACTCATAGTCATACC"
            "ACTAATCAATATCATGAAATTACCAGGTTAGTCTGTATACTACGATAAAGACACAGACTT";
        const char *reads[] = { r1, r2 };
        int32_t numReads = 2;
        int32_t max_readlength = 300;
        uint8_t *enc_qdb; bseq1_t *seq_; int32_t *cum_len;
        encode_reads(reads, numReads, &enc_qdb, &seq_, &cum_len);
        /* Start mid-read (past adapter) and use a large min_intv so the
         * forward walk will terminate via the newSmem.s < min_intv branch. */
        int32_t qpa[] = {100, 100};
        int32_t mia[] = {1000, 1000};
        int32_t rid[] = {0, 1};
        run_case(fmi, "Case 10: forward-ext min_intv break",
                 numReads, max_readlength, 19,
                 enc_qdb, qpa, mia, rid, seq_, cum_len);
        free(enc_qdb); free(seq_); free(cum_len);
    }

    /* Case 11: long reads (>151bp and >512bp) — regression for issue 44.
     * Pre-fix, readlength > 512 would trip the lockstep cap and exit.
     * Post-fix, the lockstep slot buffers are heap-allocated by the
     * driver and sized from max_readlength, so any length works.
     * Exercises a 1000bp and a 1500bp synthetic read. */
    {
        char r_1000bp[1001];
        char r_1500bp[1501];
        // Deterministic non-repetitive synthetic: a 4-base rotation with a
        // per-position twist so the SMEMs don't all collapse into one.
        const char *bases = "ACGT";
        for (int i = 0; i < 1000; i++) r_1000bp[i] = bases[(i * 7 + 3) & 3];
        r_1000bp[1000] = '\0';
        for (int i = 0; i < 1500; i++) r_1500bp[i] = bases[(i * 11 + 5) & 3];
        r_1500bp[1500] = '\0';
        const char *reads[] = { r_1000bp, r_1500bp };
        int32_t numReads = 2;
        int32_t max_readlength = 1500;
        uint8_t *enc_qdb; bseq1_t *seq_; int32_t *cum_len;
        encode_reads(reads, numReads, &enc_qdb, &seq_, &cum_len);
        int32_t qpa[] = {0, 0};
        int32_t mia[] = {1, 1};
        int32_t rid[] = {0, 1};
        run_case(fmi, "Case 11: long reads (1000bp, 1500bp)",
                 numReads, max_readlength, 19,
                 enc_qdb, qpa, mia, rid, seq_, cum_len);
        free(enc_qdb); free(seq_); free(cum_len);
    }

    /* Case 12: long reads with mid-read N — exercises the non-ACGT
     * termination path on the long-read code path that Case 11's clean
     * 4-base rotation didn't cover. The forward SMEM walk must split at
     * the N position and re-seed past it, on the heap-allocated lockstep
     * slot buffers. Two reads: 1000bp with N at pos 700, 1500bp with N
     * at pos 1100. */
    {
        char r_1000bp[1001];
        char r_1500bp[1501];
        const char *bases = "ACGT";
        for (int i = 0; i < 1000; i++) r_1000bp[i] = bases[(i * 7 + 3) & 3];
        r_1000bp[700] = 'N';
        r_1000bp[1000] = '\0';
        for (int i = 0; i < 1500; i++) r_1500bp[i] = bases[(i * 11 + 5) & 3];
        r_1500bp[1100] = 'N';
        r_1500bp[1500] = '\0';
        const char *reads[] = { r_1000bp, r_1500bp };
        int32_t numReads = 2;
        int32_t max_readlength = 1500;
        uint8_t *enc_qdb; bseq1_t *seq_; int32_t *cum_len;
        encode_reads(reads, numReads, &enc_qdb, &seq_, &cum_len);
        int32_t qpa[] = {0, 0};
        int32_t mia[] = {1, 1};
        int32_t rid[] = {0, 1};
        run_case(fmi, "Case 12: long reads with mid-read N",
                 numReads, max_readlength, 19,
                 enc_qdb, qpa, mia, rid, seq_, cum_len);
        free(enc_qdb); free(seq_); free(cum_len);
    }

    /* Case 13: reads that ACTUALLY MATCH the reference, driven from an interior
     * start_pos so the pure-backward phase (backwardExt_konly) runs a real
     * multi-step walk. Cases 1-12 use non-matching synthetic reads and/or
     * start_pos=0, so every one resolves to n_smem=0 and the backward step
     * (the code this PR changed) is never exercised — a vacuous gate. Here each
     * read is a substring of the loaded reference at a spread of offsets, with
     * an interior pivot, so the backward walk narrows a real s>1 SA interval
     * down to a unique match and emits SMEMs. This is the case that gives
     * backwardExt_konly (and the numPrev==1 fast path) genuine parity coverage;
     * with the pre-refactor comparator (which compared `l`) it FAILS on the now
     * intentionally-stale `l`, and PASSES on (rid,m,n,k,s). */
    {
        int32_t ref_len = 0;
        char *ref = read_fasta_seq(argv[1], &ref_len);
        if (!ref) {
            /* argv[1] must be a readable FASTA (CI passes fixtures/phix.fa).
             * Count the case as failed rather than skipping it, so the
             * backward-path coverage cannot drop out unnoticed. */
            total_cases++;
            fprintf(stderr, "[FAIL] Case 13: could not read reference FASTA '%s' "
                    "(backward-path parity coverage needs it)\n", argv[1]);
        } else {
            const int32_t rl = 60;               /* > minSeedLen (19) */
            const int32_t want = 6;
            int32_t off[6] = {0, 400, 900, 1500, 2600, 4000};
            int32_t numReads = 0;
            char *bufs[6] = {0};
            for (int32_t i = 0; i < want; i++) {
                if (off[i] + rl > ref_len) continue;   /* keep in-bounds */
                char *b = (char *)malloc(rl + 1);
                memcpy(b, ref + off[i], rl);
                b[rl] = '\0';
                bufs[numReads++] = b;
            }
            if (numReads == 0) {
                total_cases++;
                fprintf(stderr, "[FAIL] Case 13: reference too short for any read\n");
            } else {
                const char *reads[6];
                int32_t qpa[6], mia[6], rid[6];
                for (int32_t i = 0; i < numReads; i++) {
                    reads[i] = bufs[i];
                    /* Alternate the pivot: even reads start at the last base
                     * (longest backward walk), odd reads start mid-read (both
                     * forward and backward walks run). */
                    qpa[i] = (i & 1) ? (rl / 2) : (rl - 1);
                    mia[i] = 1;
                    rid[i] = i;   /* rid indexes the per-read seq_/cum_len layout */
                }
                uint8_t *enc_qdb; bseq1_t *seq_; int32_t *cum_len;
                encode_reads(reads, numReads, &enc_qdb, &seq_, &cum_len);
                run_case(fmi, "Case 13: real reference substrings, interior start_pos",
                         numReads, rl, 19,
                         enc_qdb, qpa, mia, rid, seq_, cum_len);
                free(enc_qdb); free(seq_); free(cum_len);
            }
            for (int32_t i = 0; i < numReads; i++) free(bufs[i]);
            free(ref);
        }
    }

    /* Case 13z: slot recycle at the AMD Zen 5 host-default width (32 unless a -D
     * override clamps it), the width Zen 5 hosts run by default. Cases 1-13 run
     * at the compile-time default and never put more reads than slots through
     * real reference matches; here 70 reference substrings (more than twice the
     * width) with interior pivots make every slot recycle with non-empty
     * prev/match state, and the case must emit SMEMs. The width is pinned
     * through the runtime global and restored afterwards. */
    {
        const int32_t zen5_n = bwa3_smem_host_default_width(1, 0x1a, 0x02);
        int32_t ref_len = 0;
        char *ref = read_fasta_seq(argv[1], &ref_len);
        enum { NZ = 70 };
        const int32_t rl = 60;
        const int32_t step = 73;   /* 69 * 73 + 60 = 5097 <= phiX's 5386 bp */
        if (!ref) {
            /* As for Case 13: a missing reference fails the case. */
            total_cases++;
            fprintf(stderr, "[FAIL] Case 13z: could not read reference FASTA '%s'\n", argv[1]);
        } else if ((NZ - 1) * step + rl > ref_len) {
            fprintf(stderr, "[FAIL] Case 13z: reference too short (%d bp) for %d reads\n",
                    ref_len, (int)NZ);
            total_cases++;
            free(ref);
        } else if (zen5_n >= NZ) {
            /* Only a -D build raising the widths gets here; Case 14 covers
             * slot recycle at every width up to SMEM_LOCKSTEP_N_MAX. */
            fprintf(stderr, "[SKIP] Case 13z: host-default width %d leaves no recycle "
                    "with %d reads (Case 14 covers it)\n", zen5_n, (int)NZ);
            free(ref);
        } else {
            char *bufs[NZ];
            const char *reads[NZ];
            int32_t qpa[NZ], mia[NZ], rid[NZ];
            for (int32_t i = 0; i < NZ; i++) {
                bufs[i] = (char *)malloc(rl + 1);
                memcpy(bufs[i], ref + i * step, rl);
                bufs[i][rl] = '\0';
                reads[i] = bufs[i];
                qpa[i] = (i & 1) ? (rl / 2) : (rl - 1);
                mia[i] = 1;
                rid[i] = i;
            }
            uint8_t *enc_qdb; bseq1_t *seq_; int32_t *cum_len;
            encode_reads(reads, NZ, &enc_qdb, &seq_, &cum_len);
            const int32_t saved_n = g_smem_lockstep_n;
            g_smem_lockstep_n = zen5_n;
            char name[128];
            snprintf(name, sizeof(name),
                     "Case 13z: numReads(%d) > N(%d), real reference substrings at the Zen 5 host-default width",
                     (int)NZ, (int)zen5_n);
            run_case(fmi, name, NZ, rl, 19, enc_qdb, qpa, mia, rid, seq_, cum_len,
                     /*expect_min_smem=*/NZ);
            g_smem_lockstep_n = saved_n;
            free(enc_qdb); free(seq_); free(cum_len);
            for (int32_t i = 0; i < NZ; i++) free(bufs[i]);
            free(ref);
        }
    }

    /* Case 14: many more reads than slots, with real matches, at several widths.
     * A finished read's slot is refilled at once, so a read can land in any
     * slot and finished reads are emitted in input order from the parked list.
     * Cases 1-13 never recycle a slot through real reference matches (Cases 4
     * and 6 recycle, but with synthetic reads), and Case 13z recycles only at
     * one width with uniform reads, so this drives 2 * N_MAX + 22 reads
     * (150 at the default N_MAX of 64, so every slot is reused at least once even
     * at the maximum width) of mixed lengths and pivots -- they retire out of
     * input order -- at a small width, the shipped default and the maximum. */
    {
        int32_t ref_len = 0;
        char *ref = read_fasta_seq(argv[1], &ref_len);
        if (!ref || ref_len < 200) {
            /* As for Case 13: a missing or short reference fails the case. */
            total_cases++;
            fprintf(stderr, "[FAIL] Case 14: could not read a reference FASTA of at least "
                    "200 bases from '%s' (slot-recycle parity coverage needs it)\n", argv[1]);
        } else {
            /* EXTRA keeps NR at 150 for the default N_MAX of 64. */
            enum { EXTRA = 22, NR = 2 * SMEM_LOCKSTEP_N_MAX + EXTRA, MAXL = 90 };
            char *bufs[NR];
            const char *reads[NR];
            int32_t qpa[NR], mia[NR], rid[NR];
            /* With the single-contig, all-ACGT phiX fixture, each read without
             * an injected N is an exact reference substring of at least 30
             * bases, longer than the 19-base minimum, so it must emit at least
             * one SMEM through its pivot. Reads with an N may emit too, so the
             * check is per read rather than a total. */
            bool need_smem[NR];
            for (int32_t i = 0; i < NR; i++) {
                const int32_t len = 30 + (i * 7) % (MAXL - 29);           /* 30..90 */
                const int32_t off = (int32_t)(((int64_t)i * 977) % (ref_len - MAXL));
                bufs[i] = (char *)malloc(len + 1);
                memcpy(bufs[i], ref + off, len);
                bufs[i][len] = '\0';
                if (i % 11 == 5) bufs[i][len / 3] = 'N';                   /* early stop */
                need_smem[i] = (i % 11 != 5);
                reads[i] = bufs[i];
                qpa[i] = (i % 3 == 0) ? len - 1 : (i % 3 == 1) ? len / 2 : 0;
                mia[i] = 1;
                rid[i] = i;
            }
            uint8_t *enc_qdb; bseq1_t *seq_; int32_t *cum_len;
            encode_reads(reads, NR, &enc_qdb, &seq_, &cum_len);
            const int32_t saved_n = g_smem_lockstep_n;
            const int32_t widths[] = {clamp_width(3), SMEM_LOCKSTEP_N, SMEM_LOCKSTEP_N_MAX};
            for (const int32_t w : widths) {
                char name[96];
                snprintf(name, sizeof name, "Case 14: %d reads > N = %d, slot recycle", (int)NR, (int)w);
                g_smem_lockstep_n = w;
                run_case(fmi, name, NR, MAXL, 19,
                         enc_qdb, qpa, mia, rid, seq_, cum_len,
                         /*expect_min_smem=*/0, need_smem);
            }
            g_smem_lockstep_n = saved_n;
            free(enc_qdb); free(seq_); free(cum_len);
            for (int32_t i = 0; i < NR; i++) free(bufs[i]);
        }
        free(ref);
    }

    /* Case 15: reads that finish far out of input order, so the driver parks
     * finished reads and emits them later. A slot is refilled as soon as its read
     * finishes, and finished reads are emitted in input order from a parked list;
     * the window admits a read only while it is fewer than 8 reads per slot,
     * rounded up to a power of two, ahead of the oldest unemitted one. A pair of
     * long reads (up to 1500 bases) every 40th input sits among short ones (25-29
     * bases), so short reads finish behind every long read still walking and the
     * window fills, leaving slots idle until the long read is out. About a third
     * of the short reads start with N, which is done the moment its slot is
     * filled, so a fill can finish at once (sometimes several in a row), including
     * the read the cursor is waiting for. There are two windows' worth of reads at
     * the maximum width, so the window can fill there too. Run at widths 1, 2, 3,
     * the default and the maximum. A model of the scheduler on these inputs (each
     * read costing about its length in steps) fills the window at every width but
     * 1, where the read the cursor waits for is the only one in flight. */
    {
        int32_t ref_len = 0;
        char *ref = read_fasta_seq(argv[1], &ref_len);
        enum { NR = 2 * LockstepParked::window_for(SMEM_LOCKSTEP_N_MAX) + 37, LONGL = 1500, SHORTL = 25 };
        if (!ref || ref_len < LONGL + 100) {
            total_cases++;
            fprintf(stderr, "[FAIL] Case 15: could not read a reference FASTA of at least "
                    "%d bases from '%s' (out-of-order parity coverage needs it)\n",
                    LONGL + 100, argv[1]);
        } else {
            std::vector<std::string> bufs(NR);
            const char *reads[NR];
            int32_t qpa[NR], mia[NR], rid[NR];
            bool need_smem[NR];
            for (int32_t i = 0; i < NR; i++) {
                /* Two adjacent long reads every 40th input. */
                const bool is_long = (i % 40 == 0) || (i % 40 == 1);
                const int32_t len = is_long ? LONGL - (i % 7) * 100 : SHORTL + (i % 5);
                const int32_t off = (int32_t)(((int64_t)i * 613) % (ref_len - LONGL));
                bufs[i].assign(ref + off, len);
                const bool first_n = !is_long && (i % 4 == 3 || i % 13 == 6);
                if (first_n) bufs[i][0] = 'N';
                need_smem[i] = !first_n;
                reads[i] = bufs[i].c_str();
                qpa[i] = 0;
                mia[i] = 1;
                rid[i] = i;
            }
            uint8_t *enc_qdb; bseq1_t *seq_; int32_t *cum_len;
            encode_reads(reads, NR, &enc_qdb, &seq_, &cum_len);
            const int32_t saved_n = g_smem_lockstep_n;
            const int32_t widths[] = {1, clamp_width(2), clamp_width(3), SMEM_LOCKSTEP_N, SMEM_LOCKSTEP_N_MAX};
            for (const int32_t w : widths) {
                char name[112];
                snprintf(name, sizeof name, "Case 15: %d reads finishing out of order, N = %d", (int)NR, (int)w);
                g_smem_lockstep_n = w;
                run_case(fmi, name, NR, LONGL, 19,
                         enc_qdb, qpa, mia, rid, seq_, cum_len,
                         /*expect_min_smem=*/0, need_smem);
            }
            g_smem_lockstep_n = saved_n;
            free(enc_qdb); free(seq_); free(cum_len);
        }
        free(ref);
    }

    /* Case 16: the unique-match text compare's exits. Cases 13-15 read only
     * forward-strand reference substrings, so the compare's other exits are
     * reached at most by chance. Here, with the forward text followed by its
     * reverse complement (2 * l_pac bases, the index's text):
     *   - reverse-complement windows (long matches in the second half);
     *   - forward and reverse-complement windows with one substituted base 3-40
     *     bases past the pivot (the text-mismatch exit; at 3 or 5 bases the
     *     match may not be unique yet, at 20 or more the compare has run 8-base
     *     words past it);
     *   - RC(ref[0..60)) followed by extra bases: the match is the text's last
     *     60 bases, so the compare stops at the text end, not at a mismatch;
     *   - ref[l_pac-30..l_pac) followed by the start of the reverse complement:
     *     a match across the forward/reverse junction;
     *   - forward and reverse-complement windows with an N 20-40 bases past the
     *     pivot (the non-ACGT exit, which resumes the next walk past the N);
     *   - short forward windows (24-36 bases from pivot 0), on both sides of the
     *     minimum read rest (LS_S1_MIN_REST, 16) below which the walk keeps
     *     stepping: the 24-base read turns unique with too few bases left;
     *   - forward windows at the start of the text (for Case 17's dropped
     *     sentinel offset leg);
     *   - two distant loci joined (for Case 17's containment checks);
     * at lengths that are not multiples of 8, from pivot 0 and interior pivots.
     * Every read emits at least one SMEM: each holds an exact match of at least
     * 20 bases through its pivot, except a pivot-0 read substituted 3 or 5 bases
     * in, which is not required to.
     * Case 17 then checks the SA tags these walks leave (see check_s1_tags),
     * including that each read built to reach the text compare took it (its
     * round-1 walk left a tagged SMEM) and that the reads built not to (an RC
     * window from its last base, the 24-base short read) did not. Case 18 runs
     * the reads on the index loaded without its pac. */
    {
        int32_t ref_len = 0;
        char *ref = read_fasta_seq(argv[1], &ref_len);
        if (!ref || ref_len < 4260) {
            total_cases++;
            fprintf(stderr, "[FAIL] Case 16: could not read a reference FASTA of at least "
                    "4260 bases from '%s' (text-compare parity coverage needs it)\n", argv[1]);
        } else {
            const std::string R(ref, (size_t)ref_len);
            std::vector<std::string> rv;
            std::vector<int32_t> qv;
            std::vector<char> nv;   /* need_smem, as char: vector<bool> has no data() */
            std::vector<signed char> ev;   /* expect_walk_tag (check_s1_tags) */
            auto add = [&](const std::string &r, int32_t q, bool need, int tag) {
                rv.push_back(r); qv.push_back(q); nv.push_back(need ? 1 : 0);
                ev.push_back((signed char)tag);
            };
            /* Reverse-complement windows. */
            const int32_t rc_len[] = {61, 75, 93, 117, 131};
            for (int32_t i = 0; i < 5; i++) {
                const int32_t len = rc_len[i];
                const std::string w = revcomp(R.substr((size_t)(300 + 701 * i), (size_t)len));
                add(w, 0, true, 1);
                add(w, len / 2, true, 1);
                add(w, len - 1, true, -1);
            }
            /* One substitution d bases past the pivot, forward and RC. */
            const int32_t subs[] = {3, 5, 20, 27, 33, 40};
            for (int32_t i = 0; i < 6; i++) {
                const int32_t d = subs[i];
                for (int strand = 0; strand < 2; strand++) {
                    const int32_t len = strand ? 115 : 101;
                    std::string w = R.substr((size_t)(150 + 611 * i + 37 * strand), (size_t)len);
                    if (strand) w = revcomp(w);
                    for (int32_t q : {0, len / 3}) {
                        std::string m = w;
                        const size_t at = (size_t)(q + d);
                        m[at] = (m[at] == 'A') ? 'G' : 'A';   /* differs from the reference */
                        add(m, q, q > 0 || d >= 20, d >= 20 ? 1 : 0);
                    }
                }
            }
            /* The match runs to the text end (the sentinel). */
            const std::string tail = revcomp(R.substr(0, 60));
            add(tail + "ACGTTGCAGTC", 0, true, 1);
            add(tail + "ACGTTGCAGTC", 10, true, 1);
            add(tail + "GGA", 30, true, 1);
            /* The match crosses the forward/reverse junction at l_pac. */
            const std::string fend = R.substr((size_t)(ref_len - 30), 30);
            add(fend + revcomp(fend), 0, true, 1);
            add(fend + revcomp(fend), 25, true, 1);
            add(fend + revcomp(fend).substr(0, 13), 0, true, 1);
            /* An N past the pivot ends the compare at a non-ACGT read base. */
            const int32_t nsubs[] = {20, 27, 33, 40};
            for (int32_t i = 0; i < 4; i++) {
                const int32_t d = nsubs[i];
                for (int strand = 0; strand < 2; strand++) {
                    const int32_t len = strand ? 109 : 97;
                    std::string w = R.substr((size_t)(2500 + 263 * i + 41 * strand), (size_t)len);
                    if (strand) w = revcomp(w);
                    for (int32_t q : {0, len / 4}) {
                        std::string m = w;
                        m[(size_t)(q + d)] = 'N';
                        add(m, q, true, 1);
                    }
                }
            }
            /* Short reads: the rest after the match turns unique straddles the
             * minimum rest for the text compare. */
            for (int32_t len : {24, 28, 32, 36})
                add(R.substr((size_t)(3500 + 7 * len), (size_t)len), 0, true,
                    len <= 24 ? -1 : 1);
            /* Matches at the start of the text: their SA walks can reach the
             * sentinel row, where a dropped sentinel offset (--compat=bwa-mem2)
             * changes the walked coordinate, so Case 17's dropped-offset leg
             * can tell a trusted tag from a walk. */
            for (int32_t p = 0; p < 16; p++)
                add(R.substr((size_t)p, 60), 0, true, 1);
            /* Two distant loci joined (forward/forward, forward/RC, RC/forward):
             * an s == 1 SMEM of one locus is not inside a tagged SMEM of the
             * other, and a broken containment rule that handed it that tag
             * would give it a coordinate that does not hold it (Case 17). */
            add(R.substr(400, 45) + R.substr(2000, 45), 0, true, 1);
            add(R.substr(400, 45) + R.substr(2000, 45), 30, true, 1);
            add(R.substr(1200, 50) + revcomp(R.substr(3000, 50)), 0, true, 1);
            add(revcomp(R.substr(1700, 40)) + R.substr(4200, 60), 10, true, 1);

            const int32_t NR = (int32_t)rv.size();
            int32_t maxl = 0;
            std::vector<const char *> reads((size_t)NR);
            std::vector<int32_t> mia((size_t)NR, 1), rid((size_t)NR);
            std::unique_ptr<bool[]> need(new bool[(size_t)NR]);
            for (int32_t i = 0; i < NR; i++) {
                reads[(size_t)i] = rv[(size_t)i].c_str();
                rid[(size_t)i] = i;
                need[i] = nv[(size_t)i] != 0;
                if ((int32_t)rv[(size_t)i].size() > maxl) maxl = (int32_t)rv[(size_t)i].size();
            }
            int32_t n_need = 0;
            for (int32_t i = 0; i < NR; i++) n_need += need[i];
            uint8_t *enc_qdb; bseq1_t *seq_; int32_t *cum_len;
            encode_reads(reads.data(), NR, &enc_qdb, &seq_, &cum_len);
            const int32_t saved_n = g_smem_lockstep_n;
            const int32_t widths[] = {clamp_width(3), SMEM_LOCKSTEP_N};
            for (const int32_t w : widths) {
                char name[112];
                snprintf(name, sizeof name, "Case 16: %d text-compare exit reads, N = %d", (int)NR, (int)w);
                g_smem_lockstep_n = w;
                run_case(fmi, name, NR, maxl, 19,
                         enc_qdb, qv.data(), mia.data(), rid.data(), seq_, cum_len,
                         /*expect_min_smem=*/n_need, need.get());
            }
            for (const int32_t w : widths) {
                char name[112];
                snprintf(name, sizeof name, "Case 17: SA tags of %d text-compare reads, N = %d", (int)NR, (int)w);
                g_smem_lockstep_n = w;
                check_s1_tags(fmi, name, NR, maxl, enc_qdb, qv.data(), seq_, cum_len,
                              R + revcomp(R), ev.data());
            }
            /* Case 18: the same reads on the index loaded without its pac
             * (as --meth and the fmi_seed_api facade load it): the text compare
             * is off, the lockstep walk steps every base, and no SMEM may be
             * tagged (run_case checks info 0 on both outputs). */
            {
                FMI_search *fmi_np = new FMI_search(argv[1]);
                fmi_np->load_index(/*load_pac=*/false);
                if (fmi_np->s1_text_available()) {
                    total_cases++;
                    fprintf(stderr, "[FAIL] Case 18: s1_text_available() is true on an index "
                            "loaded without its pac\n");
                } else {
                    for (const int32_t w : widths) {
                        char name[112];
                        snprintf(name, sizeof name, "Case 18: %d text-compare exit reads, no pac, N = %d",
                                 (int)NR, (int)w);
                        g_smem_lockstep_n = w;
                        run_case(fmi_np, name, NR, maxl, 19,
                                 enc_qdb, qv.data(), mia.data(), rid.data(), seq_, cum_len,
                                 /*expect_min_smem=*/n_need, need.get());
                    }
                }
                delete fmi_np;
            }
            g_smem_lockstep_n = saved_n;
            free(enc_qdb); free(seq_); free(cum_len);
        }
        free(ref);
    }

    /* Case 20: reads shorter than the seed jump depth. A forward walk takes its
     * first steps from the table up to the read end or the first N, so walks
     * from every start of reads of 1-14 bases (and of a 14-base read with an N
     * at each position) cover table runs of every length up past the default
     * depth, at min_intv 1 and at min_intv values that end a run early. */
    {
        const std::vector<std::string> rs = jump_depth_reads();
        std::vector<const char *> reads;
        for (const std::string &r : rs) reads.push_back(r.c_str());
        uint8_t *enc_qdb; bseq1_t *seq_; int32_t *cum_len;
        encode_reads(reads.data(), (int32_t)reads.size(), &enc_qdb, &seq_, &cum_len);
        std::vector<int32_t> qpa, rid;
        for (size_t r = 0; r < rs.size(); r++)
            for (size_t x = 0; x < rs[r].size(); x++) {
                qpa.push_back((int32_t)x);
                rid.push_back((int32_t)r);
            }
        const int32_t n = (int32_t)qpa.size();
        const int32_t min_intvs[] = {1, 4, 30};
        for (const int32_t mi : min_intvs) {
            std::vector<int32_t> q = qpa, r = rid, m((size_t)n, mi);
            char name[96];
            snprintf(name, sizeof name, "Case 20: walks from every start of reads of 1-14 bases, min_intv %d",
                     (int)mi);
            run_case(fmi, name, n, /*max_readlength=*/14, /*minSeedLen=*/3,
                     enc_qdb, q.data(), m.data(), r.data(), seq_, cum_len,
                     /*expect_min_smem=*/1);
        }
        free(enc_qdb); free(seq_); free(cum_len);
    }

    /* Case 19: s1_tag_contained's rules (no index needed). */
    check_s1_tag_contained_rules();

    fprintf(stderr, "%d / %d cases passed\n", passed_cases, total_cases);
    delete fmi;
    return (passed_cases == total_cases) ? 0 : 1;
}
