/* Thin wrapper implementation of the opaque-handle facade declared in
 * fmi_seed_api.h. Every FMI_search member used here is public (see
 * FMI_search.h); this file exists so external consumers (e.g. minibwa's
 * cp_occ backend) can link against a minimal, stable surface instead of
 * including FMI_search.h and depending on the full class. */

#include "fmi_seed_api.h"
#include "FMI_search.h"

struct FmiSeed { FMI_search *fmi; };

FmiSeed *fmi_seed_open(const char *prefix) {
    // Allocate the handle first (value-initialized, so h->fmi == nullptr), then
    // build the FMI_search inside a try/catch. A throw out of `new FMI_search`,
    // `load_index()`, or the handle allocation itself leaves nothing to leak: on
    // a throw after the handle exists we delete both (delete of a null h->fmi is
    // a no-op), and if `new FmiSeed()` itself throws there is nothing to clean up.
    FmiSeed *h = new FmiSeed();
    try {
        h->fmi = new FMI_search(prefix);
        // load_pac=false: this facade exposes only FM-index seeding data
        // (cp_occ/count/sentinel/SA resolution); the 2-bit packed reference
        // is never read through it, so skip the ~1.6 GB load (see
        // FMI_search::load_index's doc comment).
        h->fmi->load_index(/*load_pac=*/false);
    } catch (...) {
        delete h->fmi;                    // null-safe if `new FMI_search` threw
        delete h;
        throw;
    }
    return h;
}

void fmi_seed_close(FmiSeed *h) { if (!h) return; delete h->fmi; delete h; }

const CP_OCC *fmi_seed_cp_occ(const FmiSeed *h) {
    if (!h) return NULL;                  // precondition: h from a live fmi_seed_open()
    return (const CP_OCC *)h->fmi->cp_occ_data();
}

const int64_t *fmi_seed_count(const FmiSeed *h) {
    if (!h) return NULL;                  // precondition: h from a live fmi_seed_open()
    return h->fmi->count_data();
}

int64_t fmi_seed_sentinel(const FmiSeed *h) {
    if (!h) return 0;                     // precondition: h from a live fmi_seed_open()
    return h->fmi->sentinel_index;
}

void fmi_seed_sa_prefetch(FmiSeed *h, SMEM *smems, int64_t *coords,
        int64_t *coord_counts, int64_t n, int32_t max_occ, int tid, int64_t *id) {
    if (!h) return;                       // precondition: h from a live fmi_seed_open()
    // FMI_search::get_sa_entries_prefetch divides by max_occ for any SMEM
    // with s > max_occ; max_occ <= 0 is a division by zero (or an invalid
    // negative staging size) rather than a meaningful "resolve nothing"
    // request. Reject it as a safe no-op instead of forwarding it.
    if (max_occ <= 0) return;
    h->fmi->get_sa_entries_prefetch(smems, coords, coord_counts, n, max_occ, tid, *id);
}
