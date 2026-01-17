#ifndef BWA_SAFESTRINGLIB_H
#define BWA_SAFESTRINGLIB_H

/*
 * Wrapper for safestringlib headers.
 *
 * On macOS, memset_s is declared in <string.h> with a 4-argument signature
 * (C11 Annex K), which conflicts with safestringlib's 3-argument version.
 * bwa-mem2 doesn't call memset_s, so we hide safestringlib's declaration
 * by renaming it during inclusion.
 */

#ifdef __APPLE__
#define memset_s _safestringlib_memset_s
#endif

#ifdef __cplusplus
extern "C" {
#endif
#include "safe_mem_lib.h"
#include "safe_str_lib.h"
#ifdef __cplusplus
}
#endif

#ifdef __APPLE__
#undef memset_s
#endif

#endif /* BWA_SAFESTRINGLIB_H */
