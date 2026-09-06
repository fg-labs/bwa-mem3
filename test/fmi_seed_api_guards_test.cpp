/* Regression test for the fmi_seed_api.h facade's input guards.
 *
 * fmi_seed_sa_prefetch forwards directly to
 * FMI_search::get_sa_entries_prefetch, which divides by max_occ for any
 * SMEM with s > max_occ. Before the facade validated max_occ, a caller
 * passing max_occ <= 0 alongside a non-empty SMEM crashed with a division
 * by zero (or drove an invalid negative staging size for a negative
 * max_occ) instead of getting a documented safe no-op.
 *
 * Usage: fmi_seed_api_guards_test <bwa-mem3 index prefix>
 *   e.g. fmi_seed_api_guards_test test/fixtures/phix.fa
 */

#include "fmi_seed_api.h"
#include <cstdio>
#include <cstdlib>

static void expect(bool cond, const char *what) {
    if (!cond) {
        fprintf(stderr, "FAIL: %s\n", what);
        exit(1);
    }
}

int main(int argc, char **argv) {
    if (argc != 2) {
        fprintf(stderr, "usage: %s <bwa-mem3 index prefix>\n", argv[0]);
        return 1;
    }

    FmiSeed *h = fmi_seed_open(argv[1]);
    expect(h != NULL, "fmi_seed_open returned NULL");

    // A small, in-bounds SMEM interval (any indexed reference has at least
    // this many BWT rows) with s > 0 so the staging loop's step computation
    // actually reaches the max_occ divide this test targets.
    SMEM smem;
    smem.rid = 0; smem.m = 0; smem.n = 0;
    smem.k = 0; smem.l = 0; smem.s = 5;

    int64_t coords[8];
    int64_t coord_count, id;

    // A sentinel distinguishes "the facade left the output untouched" from "the
    // facade wrote a value that happens to equal our pre-set" -- the no-op guard
    // must do the former, so seed every output slot with SENT and assert the
    // slots it must not touch still hold SENT afterwards.
    const int64_t SENT = -999;

    // max_occ == 0: must be a safe no-op, not a forwarded call (which would
    // divide smem.s by max_occ). The guard returns before touching any output,
    // so coord_count, id, and coords must all remain at their sentinels.
    coord_count = SENT; id = SENT;
    for (int i = 0; i < 8; i++) coords[i] = SENT;
    fmi_seed_sa_prefetch(h, &smem, coords, &coord_count, 1, 0, 0, &id);
    expect(coord_count == SENT, "max_occ == 0 must leave coord_count untouched");
    expect(id == SENT, "max_occ == 0 must leave id untouched");
    for (int i = 0; i < 8; i++)
        expect(coords[i] == SENT, "max_occ == 0 must not write any coordinate");

    // Negative max_occ: same guard, same reasoning.
    coord_count = SENT; id = SENT;
    for (int i = 0; i < 8; i++) coords[i] = SENT;
    fmi_seed_sa_prefetch(h, &smem, coords, &coord_count, 1, -3, 0, &id);
    expect(coord_count == SENT, "negative max_occ must leave coord_count untouched");
    expect(id == SENT, "negative max_occ must leave id untouched");
    for (int i = 0; i < 8; i++)
        expect(coords[i] == SENT, "negative max_occ must not write any coordinate");

    // A valid positive max_occ must still resolve entries -- the guard must
    // not swallow legitimate calls. With smem.s == 5 and max_occ == 4, the
    // staging loop steps by 1 and stops at c < max_occ, so it resolves exactly
    // four entries (coord_count) and stages exactly four buffer slots (id).
    // id is the in/out staging base offset, so it starts at 0 here.
    for (int i = 0; i < 8; i++) coords[i] = SENT;
    coord_count = 0; id = 0;
    fmi_seed_sa_prefetch(h, &smem, coords, &coord_count, 1, 4, 0, &id);
    expect(coord_count == 4, "valid max_occ must resolve exactly min(s, max_occ) entries");
    expect(id == 4, "valid max_occ must stage exactly min(s, max_occ) buffer slots");

    // Beyond the count, validate the resolved records themselves so a corrupt
    // coordinate cannot pass on counts alone. The four resolved entries are SA
    // positions into the reference, so each must be a valid (non-negative)
    // offset, all four map distinct BWT rows so they must be pairwise distinct,
    // and exactly four slots may be written -- the fifth must stay at SENT.
    // (These are self-consistent structural invariants, not an exact-value
    // oracle; see the reply on this thread for why an upstream-bwa value oracle
    // is out of scope for this facade regression test.)
    for (int i = 0; i < 4; i++)
        expect(coords[i] >= 0, "resolved SA coordinate must be a valid (non-negative) offset");
    for (int i = 0; i < 4; i++)
        for (int j = i + 1; j < 4; j++)
            expect(coords[i] != coords[j], "resolved SA coordinates must be pairwise distinct");
    expect(coords[4] == SENT, "exactly four coordinates must be staged -- no overrun");

    fmi_seed_close(h);
    printf("PASS\n");
    return 0;
}
