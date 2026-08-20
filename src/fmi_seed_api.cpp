/* Thin wrapper implementation of the opaque-handle facade declared in
 * fmi_seed_api.h. Every FMI_search member used here is public (see
 * FMI_search.h); this file exists so external consumers (e.g. minibwa's
 * cp_occ backend) can link against a minimal, stable surface instead of
 * including FMI_search.h and depending on the full class. */

#include "fmi_seed_api.h"
#include "FMI_search.h"

struct FmiSeed { FMI_search *fmi; };

FmiSeed *fmi_seed_open(const char *prefix) {
    FmiSeed *h = new FmiSeed();
    h->fmi = new FMI_search(prefix);
    h->fmi->load_index();                 // defaults: load_pac=true, n_threads=1
    return h;
}

void fmi_seed_close(FmiSeed *h) { if (!h) return; delete h->fmi; delete h; }

const CP_OCC *fmi_seed_cp_occ(const FmiSeed *h) {
    return (const CP_OCC *)h->fmi->cp_occ_data();
}

const int64_t *fmi_seed_count(const FmiSeed *h) { return h->fmi->count_data(); }

int64_t fmi_seed_sentinel(const FmiSeed *h) { return h->fmi->sentinel_index; }

void fmi_seed_sa_prefetch(FmiSeed *h, SMEM *smems, int64_t *coords,
        int64_t *coord_counts, int64_t n, int32_t max_occ, int tid, int64_t *id) {
    h->fmi->get_sa_entries_prefetch(smems, coords, coord_counts, n, max_occ, tid, *id);
}
