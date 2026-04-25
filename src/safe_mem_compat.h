/*************************************************************************************
                          The MIT License

   BWA-MEM2  (Sequence alignment using Burrows-Wheeler Transform),
   Copyright (C) 2019  Intel Corporation, Heng Li.

   Permission is hereby granted, free of charge, to any person obtaining
   a copy of this software and associated documentation files (the
   "Software"), to deal in the Software without restriction, including
   without limitation the rights to use, copy, modify, merge, publish,
   distribute, sublicense, and/or sell copies of the Software, and to
   permit persons to whom the Software is furnished to do so, subject to
   the following conditions:

   The above copyright notice and this permission notice shall be
   included in all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
   BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.

   Safe Memory Compatibility Header for macOS
   Resolves conflicts with macOS built-in memset_s
   Based on PR #281 from BenjaminDEMAILLE
*****************************************************************************************/

#ifndef SAFE_MEM_COMPAT_H
#define SAFE_MEM_COMPAT_H

/*
 * macOS defines memset_s in string.h (C11 Annex K), which conflicts with
 * safestringlib's memset_s. The signatures are different:
 *   macOS:        errno_t memset_s(void *s, rsize_t smax, int c, rsize_t n);
 *   safestringlib: errno_t memset_s(void *dest, rsize_t dmax, uint8_t value);
 *
 * We need to prevent safestringlib from declaring its memset_s on macOS.
 */

#ifdef __APPLE__
    /* Include standard headers */
    #include <string.h>
    #include <stdlib.h>
    #include <stddef.h>
    #include <stdint.h>

    /*
     * Trick to hide safestringlib's memset_s declaration:
     * We redefine 'memset_s' just before safe_mem_lib.h is included,
     * so the extern declaration becomes 'extern errno_t _hidden_safestringlib_memset_s(...)'
     * which doesn't conflict with macOS's memset_s.
     * Then we undef it after so we can use macOS's memset_s normally.
     */
    #define APPLE_SAFESTRINGLIB_COMPAT 1

#endif /* __APPLE__ */

#endif /* SAFE_MEM_COMPAT_H */
