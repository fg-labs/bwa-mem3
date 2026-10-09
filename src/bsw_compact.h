/* bsw_compact.h -- the plan for same-row lane compaction in the 8-bit banded extension kernels
 * (bandedSWA_compact.inc on NEON, bandedSWA_compact256.inc on AVX2, bandedSWA_compact512.inc
 * on AVX-512BW).
 *
 * A batch runs in superblocks of K lane groups (vectors) that advance in row lockstep, so every
 * lane of a superblock is at the same target row i. When pairs finish, the superblock may need
 * fewer vectors than it has active; the plan below then moves the live lanes of the emptiest
 * active vector into dead slots ("holes") of the other active vectors and retires it. Lanes
 * keep their absolute-column frame (every lane is at row i), so a moved lane's DP state is a
 * plain copy of its slot, and a vector's column range stays the union of same-row bands, as in
 * the plain kernel.
 *
 * Exactness rests on lane independence: a lane's result does not depend on which lanes share
 * its vector or on the vector's column range, provided the range covers the lane's own band
 * (see the lane-independence note at the top of bandedSWA.cpp). The plan is pure bookkeeping
 * on live-lane bitmasks; it never looks at DP state.
 */
#ifndef BSW_COMPACT_H
#define BSW_COMPACT_H

#include <assert.h>
#include <stdint.h>

struct BswMove {
    uint8_t sv, sl;   /* source vector, lane */
    uint8_t dv, dl;   /* destination vector, lane */
};

/* live[v]: bit l set iff lane l of vector v is live. active[v]: vector v still runs. L is the
 * lane count (16, 32 or 64), nv <= BSW_COMPACT_GROUPS_MAX (bandedSWA.h).
 *
 * Vectors with no live lane are retired first. Then, while the number of active vectors exceeds
 * ceil(total_live / L), the active vector with the fewest live lanes (lowest index on ties) is
 * evacuated: each of its live lanes, in lane order, goes to the first hole of the other active
 * vectors (lowest vector, then lowest lane). The donor always empties: the other active vectors
 * then hold at least as many holes as it has live lanes. live[] and active[] are updated, the
 * moves are appended to mv (capacity nv * L), and the number of moves is returned.
 *
 * So the moves come donor by donor, and a donor's moves fill one destination vector before the
 * next. The NEON and AVX-512BW drivers apply each such run (one source and one destination
 * vector) in one pass (move_run in bandedSWA_compact.inc and bandedSWA_compact512.inc); that is
 * correct for any order, but emitting the moves in another order would split the runs and cost
 * speed. */
static inline int bsw_compact_plan(uint64_t *live, uint8_t *active, int nv, int L, BswMove *mv)
{
    const uint64_t full = (L >= 64) ? ~(uint64_t) 0 : (((uint64_t) 1 << L) - 1);
    int total = 0, nact = 0;
    for (int v = 0; v < nv; v++) {
        if (active[v] && live[v] == 0) active[v] = 0;
        if (active[v]) { nact++; total += __builtin_popcountll(live[v]); }
    }
    int nm = 0;
    while (nact > (total + L - 1) / L) {
        int d = -1;
        for (int v = 0; v < nv; v++)
            if (active[v] && (d < 0 || __builtin_popcountll(live[v]) < __builtin_popcountll(live[d]))) d = v;
        for (int s = 0; s < L; s++) {
            if (!(live[d] >> s & 1u)) continue;
            int rv = -1, rl = -1;
            for (int v = 0; v < nv && rv < 0; v++) {
                if (!active[v] || v == d || live[v] == full) continue;
                rv = v; rl = __builtin_ctzll(~live[v] & full);
            }
            /* nact > ceil(total / L) leaves the other active vectors at least popcount(live[d])
             * holes, so a hole always exists */
            if (rv < 0) { assert(!"bsw_compact_plan: no hole for an evacuated lane"); return nm; }
            mv[nm].sv = (uint8_t) d; mv[nm].sl = (uint8_t) s;
            mv[nm].dv = (uint8_t) rv; mv[nm].dl = (uint8_t) rl; nm++;
            live[rv] |= (uint64_t) 1 << rl; live[d] &= ~((uint64_t) 1 << s);
        }
        active[d] = 0; nact--;
    }
    return nm;
}

#endif
