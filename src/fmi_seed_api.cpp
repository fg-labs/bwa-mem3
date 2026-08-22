/* Thin wrapper implementation of the opaque-handle facade declared in
 * fmi_seed_api.h. Every FMI_search member used here is public (see
 * FMI_search.h); this file exists so external consumers (e.g. minibwa's
 * cp_occ backend) can link against a minimal, stable surface instead of
 * including FMI_search.h and depending on the full class. */

#include "fmi_seed_api.h"
#include "FMI_search.h"

struct FmiSeed { FMI_search *fmi; };

FmiSeed *fmi_seed_open(const char *prefix) {
    // Build the FMI_search fully before allocating the handle, so a throw out
    // of `new FMI_search` or `load_index()` leaves nothing for the caller to
    // leak -- there is no partially-constructed FmiSeed to clean up.
    FMI_search *fmi = new FMI_search(prefix);
    try {
        fmi->load_index();                // defaults: load_pac=true, n_threads=1
    } catch (...) {
        delete fmi;
        throw;
    }
    FmiSeed *h = new FmiSeed();
    h->fmi = fmi;
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
    h->fmi->get_sa_entries_prefetch(smems, coords, coord_counts, n, max_occ, tid, *id);
}
