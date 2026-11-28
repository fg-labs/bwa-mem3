#ifndef FMI_SEED_API_H
#define FMI_SEED_API_H

#include <stdint.h>
#include <stddef.h>   /* offsetof (the CP_OCC_PLANES static asserts) */

/* ---- Blessed public seeding surface for external consumers (e.g. minibwa). ----
 * Stable: changes here are API-breaking and must be versioned. bwa-mem3's own
 * code and this file share these definitions; do not fork them. */

#define CP_BLOCK_SIZE 64
#define CP_FILENAME_SUFFIX ".bwt.2bit.64"
#define CP_MASK 63
#define CP_SHIFT 6

/* The on-disk occurrence checkpoint of <prefix>.bwt.2bit.64: one per 64 BWT
 * rows, the four bases' counts before the block and one one-hot word per base,
 * bit 63 - j for the block's row j. The sentinel row and the padding past the
 * last row set no bit. */
typedef struct checkpoint_occ_scalar
{
    int64_t cp_count[4];
    uint64_t one_hot_bwt_str[4];
}CP_OCC;

/* The compact in-memory occurrence checkpoint, built from CP_OCC at load when
 * the compact layout is selected (mem --compact-occ or BWA3_COMPACT_OCC=1, or
 * fmi_seed_open_layout with FMI_SEED_LAYOUT_COMPACT; the default keeps the
 * file's CP_OCC records): 32 bytes, two per 64-byte line,
 * none straddling one. `hi` and `lo` hold the two bits of
 * each row's base code (A=00 C=01 G=10 T=11), bit 63 - j for row j, MSB-first
 * like the file: hi = G|T, lo = C|T. The sentinel row and the padding read as A
 * here; the occurrence counts correct for the sentinel (FMI_search::occ_at).
 * cp_count[c] is the count of c before the block, less cp_count_base[span][c]
 * for the block's 2^32-row span on a reference where some base occurs 2^32
 * times or more (FMI_search::cp_wide); otherwise the full count. */
typedef struct checkpoint_occ_planes
{
    uint32_t cp_count[4];
    uint64_t hi;
    uint64_t lo;
} __attribute__((aligned(32))) CP_OCC_PLANES;
/* 32 bytes (two per 64-byte line, none straddling one), and lo right after hi
 * (the NEON path loads [hi, lo] as one vector). */
#ifdef __cplusplus
static_assert(sizeof(CP_OCC_PLANES) == 32, "CP_OCC_PLANES must stay 32 bytes");
static_assert(offsetof(CP_OCC_PLANES, lo) == offsetof(CP_OCC_PLANES, hi) + 8, "CP_OCC_PLANES: lo must follow hi");
#else
_Static_assert(sizeof(CP_OCC_PLANES) == 32, "CP_OCC_PLANES must stay 32 bytes");
_Static_assert(offsetof(CP_OCC_PLANES, lo) == offsetof(CP_OCC_PLANES, hi) + 8, "CP_OCC_PLANES: lo must follow hi");
#endif

#if defined(__clang__) || defined(__GNUC__)
static inline int _mm_countbits_64(unsigned long x) {
    return __builtin_popcountl(x);
}
#endif

/* One-hot position masks: entry [i] has the top `i` bits set (entry [0] == 0),
 * i.e. one_hot_mask_array[i] == (one_hot_mask_array[i-1] >> 1) | 0x8000...
 * for i in 1..63. The FMI-index lifetime constant is identical on every path,
 * so it lives inline here as a file-scope table rather than a heap allocation
 * reached through a per-object pointer: this removes a dependent load in the
 * (~10^10-call) backwardExt/GET_OCC hot path. Values are byte-identical to the
 * former runtime-computed array. */
static const uint64_t one_hot_mask_array[64] = {
    0x0000000000000000ULL, 0x8000000000000000ULL, 0xc000000000000000ULL, 0xe000000000000000ULL,
    0xf000000000000000ULL, 0xf800000000000000ULL, 0xfc00000000000000ULL, 0xfe00000000000000ULL,
    0xff00000000000000ULL, 0xff80000000000000ULL, 0xffc0000000000000ULL, 0xffe0000000000000ULL,
    0xfff0000000000000ULL, 0xfff8000000000000ULL, 0xfffc000000000000ULL, 0xfffe000000000000ULL,
    0xffff000000000000ULL, 0xffff800000000000ULL, 0xffffc00000000000ULL, 0xffffe00000000000ULL,
    0xfffff00000000000ULL, 0xfffff80000000000ULL, 0xfffffc0000000000ULL, 0xfffffe0000000000ULL,
    0xffffff0000000000ULL, 0xffffff8000000000ULL, 0xffffffc000000000ULL, 0xffffffe000000000ULL,
    0xfffffff000000000ULL, 0xfffffff800000000ULL, 0xfffffffc00000000ULL, 0xfffffffe00000000ULL,
    0xffffffff00000000ULL, 0xffffffff80000000ULL, 0xffffffffc0000000ULL, 0xffffffffe0000000ULL,
    0xfffffffff0000000ULL, 0xfffffffff8000000ULL, 0xfffffffffc000000ULL, 0xfffffffffe000000ULL,
    0xffffffffff000000ULL, 0xffffffffff800000ULL, 0xffffffffffc00000ULL, 0xffffffffffe00000ULL,
    0xfffffffffff00000ULL, 0xfffffffffff80000ULL, 0xfffffffffffc0000ULL, 0xfffffffffffe0000ULL,
    0xffffffffffff0000ULL, 0xffffffffffff8000ULL, 0xffffffffffffc000ULL, 0xffffffffffffe000ULL,
    0xfffffffffffff000ULL, 0xfffffffffffff800ULL, 0xfffffffffffffc00ULL, 0xfffffffffffffe00ULL,
    0xffffffffffffff00ULL, 0xffffffffffffff80ULL, 0xffffffffffffffc0ULL, 0xffffffffffffffe0ULL,
    0xfffffffffffffff0ULL, 0xfffffffffffffff8ULL, 0xfffffffffffffffcULL, 0xfffffffffffffffeULL,
};

/* `info` is unconditional (not gated on DEBUG): this struct crosses a
 * library/consumer build boundary through fmi_seed_sa_prefetch, so its
 * layout must be identical regardless of the DEBUG setting the library was
 * built with vs. what a consumer TU defines when it includes this header.
 * A conditional field here would silently desync array-element offsets and
 * sizes between mismatched builds. */
typedef struct smem_struct
{
    uint64_t info; // reserved for bwa-mem3's own seeding, which may record an
                   // SMEM's already-resolved SA position here (bit 63 set) for
                   // its internal SA resolution; fmi_seed_sa_prefetch ignores
                   // info, so external callers may leave it 0
    uint32_t rid;
    uint32_t m, n;
    int64_t k, l, s;
}SMEM;

/* ---- Opaque-handle facade over FMI_search. ----
 * External consumers (e.g. minibwa's cp_occ backend) load an index, resolve
 * SA entries, and read the raw cp_occ/count/sentinel arrays through this
 * facade instead of including FMI_search.h and touching the full class.
 * Stable: signatures here are API-breaking to change. */

struct FmiSeed;                                             /* opaque handle */
typedef struct FmiSeed FmiSeed;               /* lets C consumers name the type */

#ifdef __cplusplus
extern "C" {
#endif

/* The facade's version. 1 was the unversioned surface; 2 adds fmi_seed_open_layout,
 * fmi_seed_cp_occ_planes and fmi_seed_cp_count_base. A consumer that needs one of
 * them checks FMI_SEED_API_VERSION >= 2. */
#define FMI_SEED_API_VERSION 2

/* The checkpoint layouts fmi_seed_open_layout takes. */
#define FMI_SEED_LAYOUT_DEFAULT 0   /* the file's CP_OCC records (fmi_seed_cp_occ) */
#define FMI_SEED_LAYOUT_COMPACT 1   /* CP_OCC_PLANES (fmi_seed_cp_occ_planes), as mem --compact-occ */

/* new FMI_search + load_index(load_pac=false, no seed jump table), in the
 * default checkpoint layout whatever the environment says: BWA3_COMPACT_OCC,
 * which selects `mem`'s layout, is ignored here, so fmi_seed_cp_occ() is never
 * NULL on a handle from this call. Missing, malformed or structurally inconsistent
 * index files terminate the process with a diagnostic, as the underlying loader does;
 * these failures are not reported as NULL. */
FmiSeed       *fmi_seed_open(const char *prefix);
/* fmi_seed_open in the layout the caller names (FMI_SEED_LAYOUT_*); NULL for
 * any other value. The environment is ignored here too.
 * Both open calls attach an index staged in shared memory (`bwa-mem3 shm`) only
 * when it holds the layout the call names; a stage of the other layout (`shm
 * --compact-occ`, or a default stage under FMI_SEED_LAYOUT_COMPACT) is not
 * attached, and the index loads from disk, next to the stage, with a
 * [W::bwa_shm_attach] line on stderr. */
FmiSeed       *fmi_seed_open_layout(const char *prefix, int layout);
void           fmi_seed_close(FmiSeed *h);                   /* delete */
/* The checkpoints, one record per 64 rows ((ref_seq_len >> 6) + 1 of them), in
 * the layout the handle was opened with; the other layout's accessor returns
 * NULL. Default layout: fmi_seed_cp_occ() gives the file's CP_OCC records, and
 * fmi_seed_cp_occ_planes() and fmi_seed_cp_count_base() return NULL. Compact
 * layout (fmi_seed_open_layout with FMI_SEED_LAYOUT_COMPACT only):
 * fmi_seed_cp_occ() returns NULL; read
 * fmi_seed_cp_occ_planes() with fmi_seed_cp_count_base() (NULL unless the
 * counts are span-relative, one int64_t[4] per 2^32 rows). */
const CP_OCC  *fmi_seed_cp_occ(const FmiSeed *h);
const CP_OCC_PLANES *fmi_seed_cp_occ_planes(const FmiSeed *h);
const int64_t *fmi_seed_cp_count_base(const FmiSeed *h);
const int64_t *fmi_seed_count(const FmiSeed *h);              /* count_data() */
int64_t        fmi_seed_sentinel(const FmiSeed *h);           /* sentinel_index */
/* get_sa_entries_prefetch(). max_occ <= 0 is a safe no-op (resolves nothing,
 * does not forward to FMI_search): the underlying implementation divides by
 * max_occ for any SMEM with s > max_occ.
 * *id is accumulated (+=) with the number of resolved coordinates staged into
 * coords -- i.e. sum(min(s, max_occ)) over the n SMEMs. This is invariant to
 * the internal (k,s) dedup mode: cross-read duplicate intervals skip their
 * LF-walks but still occupy their coordinate slots, so *id never shrinks when
 * dedup engages. */
void           fmi_seed_sa_prefetch(FmiSeed *h, SMEM *smems, int64_t *coords,
                    int64_t *coord_counts, int64_t n, int32_t max_occ,
                    int tid, int64_t *id);

#ifdef __cplusplus
}
#endif

/* fmi_seed_sa_intv() (sampling-interval accessor) is reserved for a later PR
 * (B2/M3); do not add it here. */

#endif /* FMI_SEED_API_H */
