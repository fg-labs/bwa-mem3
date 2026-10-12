#ifndef LOCKSTEP_WIDTH_INTERNAL_H
#define LOCKSTEP_WIDTH_INTERNAL_H

/* Pure pieces of the lockstep-width resolution, exposed only to
 * lockstep_width.cpp and its unit test (test/unit/test_lockstep_width.cpp).
 * Not part of the library interface: callers use the initializers and
 * accessors in lockstep_width.h. */

#include "lockstep_width.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Nonzero iff an AMD family 1Ah display model is Zen 5. Ranges from the Linux
 * kernel's Zen generation table (bsp_init_amd, arch/x86/kernel/cpu/amd.c):
 * 00h-2Fh, 40h-4Fh, 60h-7Fh, D0h-D7h. The Zen 6 ranges (50h-5Fh, 80h-AFh,
 * C0h-CFh, D8h-EFh) and unlisted models return 0. */
int bwa3_amd_family_1ah_model_is_zen5(uint32_t model);

/* Nonzero iff the identity is AMD Zen 5: AMD vendor, display family 1Ah, and a
 * Zen 5 model. */
int bwa3_cpu_is_amd_zen5(int is_amd, uint32_t family, uint32_t model);

/* A resolved width and where it came from (BWA3_LOCKSTEP_WIDTH_SRC_*). */
typedef struct {
    int32_t width;
    int     src;
} bwa3_lockstep_width_res_t;

/* Resolve the phase-2 SMEM width from already-classified inputs, WITHOUT
 * running the probe: `pinned` is the bwa3_lockstep_width_parse_env result (> 0
 * pin, 0 unset, -1 invalid, which resolves as unset), `probe_on` the
 * bwa3_lockstep_probe_enabled result, `host_width` the
 * bwa3_smem_host_default_width result. Returns
 *   {pinned, ENV}            for a valid pin;
 *   {0, PROBE}               when the probe is opted in: the caller runs it and
 *                            finishes with bwa3_resolve_smem_width_from_probe;
 *   {host, HOST | DEFAULT}   otherwise, host_width clamped to
 *                            [SMEM_LOCKSTEP_N, SMEM_LOCKSTEP_N_MAX]; HOST when
 *                            that differs from SMEM_LOCKSTEP_N. */
bwa3_lockstep_width_res_t bwa3_resolve_smem_width(int32_t pinned, int probe_on,
                                                  int32_t host_width);

/* Finish a PROBE resolution with the probe's raw result: a measured knee
 * (raw_mlp > 0) gives {bwa3_lockstep_width_from_probe(raw_mlp), PROBE}; a probe
 * that could not measure (0) falls back to the host default exactly as if the
 * probe were off. */
bwa3_lockstep_width_res_t bwa3_resolve_smem_width_from_probe(int32_t raw_mlp,
                                                             int32_t host_width);

/* Resolve the third-pass bwtSeed width: {pinned, ENV} for a valid pin, else
 * host_width clamped to [BWTSEED_LOCKSTEP_N, BWTSEED_LOCKSTEP_N_MAX], HOST when
 * that differs from BWTSEED_LOCKSTEP_N, DEFAULT otherwise. */
bwa3_lockstep_width_res_t bwa3_resolve_bwtseed_width(int32_t pinned,
                                                     int32_t host_width);

#ifdef __cplusplus
}
#endif

#endif /* LOCKSTEP_WIDTH_INTERNAL_H */
