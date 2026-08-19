/* The MIT License

   Copyright (c) 2008, by Attractive Chaos <attractor@live.co.uk>

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
*/

/*
  An example:

#include "kvec.h"
int main() {
	kvec_t(int) array;
	kv_init(array);
	kv_push(int, array, 10); // append
	kv_a(int, array, 20) = 5; // dynamic
	kv_A(array, 20) = 4; // static
	kv_destroy(array);
	return 0;
}
*/

/*
  2008-09-22 (0.1.0):

	* The initial version.

*/

#ifndef AC_KVEC_H
#define AC_KVEC_H

#include <stdio.h>
#include <stdlib.h>

#ifdef USE_MALLOC_WRAPPERS
#  include "malloc_wrap.h"
#endif

/* Checked realloc for the growth macros below.
 *
 * The stock klib macros assigned realloc()'s result straight back onto (v).a:
 * on failure realloc returns NULL, which both leaks the original buffer and is
 * then written through by the immediately-following element store. Routing
 * every resize through this helper closes both holes at the single source
 * point, leaving the ~40 kv_* call sites untouched.
 *
 * OOM is treated as unrecoverable and aborts loudly with a diagnostic, matching
 * the house convention for allocation failure on the read/align paths (see
 * src/fr_fastq.c, src/read_arena.h) rather than threading a failure code
 * through every caller -- several of which sit in the chaining/pairing inner
 * loops. The abort is noreturn, so the common (non-growing) path pays nothing:
 * the check lives only on the realloc branch, which already costs orders of
 * magnitude more than a compare. A zero-size request is honoured rather than
 * mistaken for failure: realloc(p, 0) may legally return NULL. */
static inline void *kv_realloc_or_die(void *p, size_t n)
{
	void *q = realloc(p, n);
	if (q == NULL && n != 0) {
		fprintf(stderr, "[kvec] out of memory: failed to (re)allocate %zu bytes\n", n);
		abort();
	}
	return q;
}

#define kv_roundup32(x) (--(x), (x)|=(x)>>1, (x)|=(x)>>2, (x)|=(x)>>4, (x)|=(x)>>8, (x)|=(x)>>16, ++(x))
typedef struct {
	size_t n,m;
	int *a;
}  kvecint;

#define kvec_t(type) struct { size_t n, m; type *a; }
#define kv_init(v) ((v).n = (v).m = 0, (v).a = 0)
#define kv_destroy(v) free((v).a)
#define kv_A(v, i) ((v).a[(i)])
#define kv_pop(v) ((v).a[--(v).n])
#define kv_size(v) ((v).n)
#define kv_max(v) ((v).m)

#define kv_resize(type, v, s)  ((v).m = (s), (v).a = (type*)kv_realloc_or_die((v).a, sizeof(type) * (v).m))

#define kv_copy(type, v1, v0) do {							\
		if ((v1).m < (v0).n) kv_resize(type, v1, (v0).n);	\
		(v1).n = (v0).n;									\
		memcpy((v1).a, (v0).a, sizeof(type) * (v0).n);		\
	} while (0)												\

#define kv_push(type, v, x) do {									\
		if ((v).n == (v).m) {										\
			(v).m = (v).m? (v).m<<1 : 2;							\
			(v).a = (type*)kv_realloc_or_die((v).a, sizeof(type) * (v).m); \
		}															\
		(v).a[(v).n++] = (x);										\
	} while (0)

#define kv_pushp(type, v) ((((v).n == (v).m)?							\
						   ((v).m = ((v).m? (v).m<<1 : 2),				\
							(v).a = (type*)kv_realloc_or_die((v).a, sizeof(type) * (v).m), 0)	\
						   : 0), &(v).a[(v).n++])

#define kv_a(type, v, i) (((v).m <= (size_t)(i)? \
						  ((v).m = (v).n = (i) + 1, kv_roundup32((v).m), \
						   (v).a = (type*)kv_realloc_or_die((v).a, sizeof(type) * (v).m), 0) \
						  : (v).n <= (size_t)(i)? (v).n = (i) + 1 \
						  : 0), (v).a[(i)])

#endif
