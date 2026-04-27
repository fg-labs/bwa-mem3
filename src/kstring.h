/* The MIT License

   Copyright (c) by Attractive Chaos <attractivechaos@aol.co.uk>

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

   Modified Copyright (C) 2020 Intel Corporation, Heng Li.
   Contacts: Vasimuddin Md <vasimuddin.md@intel.com>; Sanchit Misra <sanchit.misra@intel.com>;
   Heng Li <hli@jimmy.harvard.edu> 
*/

#ifndef KSTRING_H
#define KSTRING_H

#include <stdlib.h>
#include <string.h>
#include "memcpy_bwamem.h"

#ifdef USE_MALLOC_WRAPPERS
#  include "malloc_wrap.h"
#endif

#ifndef kroundup32
#define kroundup32(x) (--(x), (x)|=(x)>>1, (x)|=(x)>>2, (x)|=(x)>>4, (x)|=(x)>>8, (x)|=(x)>>16, ++(x))
#endif

#ifndef KSTRING_T
#define KSTRING_T kstring_t
typedef struct __kstring_t {
	size_t l, m;
	char *s;
} kstring_t;
#endif

static inline void ks_resize(kstring_t *s, size_t size)
{
	if (s->m < size) {
		s->m = size;
		kroundup32(s->m);
		s->s = (char*)realloc(s->s, s->m);
	}
}

static inline int kputsn(const char *p, int l, kstring_t *s)
{
	if (s->l + l + 1 >= s->m) {
		s->m = s->l + l + 2;
		kroundup32(s->m);
		s->s = (char*)realloc(s->s, s->m);
	}
	memcpy_bwamem(s->s + s->l, s->m - s->l, p, l, __FILE__, __LINE__);
	s->l += l;
	s->s[s->l] = 0;
	return l;
}

static inline int kputs(const char *p, kstring_t *s)
{
	return kputsn(p, strlen(p), s);
}

static inline int kputc(int c, kstring_t *s)
{
	if (s->l + 1 >= s->m) {
		s->m = s->l + 2;
		kroundup32(s->m);
		s->s = (char*)realloc(s->s, s->m);
	}
	s->s[s->l++] = c;
	s->s[s->l] = 0;
	return c;
}

static inline int kputw(int c, kstring_t *s)
{
	char buf[16];
	int l, x;
	if (c == 0) return kputc('0', s);
	for (l = 0, x = c < 0? -c : c; x > 0; x /= 10) buf[l++] = x%10 + '0';
	if (c < 0) buf[l++] = '-';
	if (s->l + l + 1 >= s->m) {
		s->m = s->l + l + 2;
		kroundup32(s->m);
		s->s = (char*)realloc(s->s, s->m);
	}
	for (x = l - 1; x >= 0; --x) s->s[s->l++] = buf[x];
	s->s[s->l] = 0;
	return 0;
}

static inline int kputuw(unsigned c, kstring_t *s)
{
	char buf[16];
	int l, i;
	unsigned x;
	if (c == 0) return kputc('0', s);
	for (l = 0, x = c; x > 0; x /= 10) buf[l++] = x%10 + '0';
	if (s->l + l + 1 >= s->m) {
		s->m = s->l + l + 2;
		kroundup32(s->m);
		s->s = (char*)realloc(s->s, s->m);
	}
	for (i = l - 1; i >= 0; --i) s->s[s->l++] = buf[i];
	s->s[s->l] = 0;
	return 0;
}

static inline int kputl(long c, kstring_t *s)
{
	char buf[32];
	long l, x;
	if (c == 0) return kputc('0', s);
	for (l = 0, x = c < 0? -c : c; x > 0; x /= 10) buf[l++] = x%10 + '0';
	if (c < 0) buf[l++] = '-';
	if (s->l + l + 1 >= s->m) {
		s->m = s->l + l + 2;
		kroundup32(s->m);
		s->s = (char*)realloc(s->s, s->m);
	}
	for (x = l - 1; x >= 0; --x) s->s[s->l++] = buf[x];
	s->s[s->l] = 0;
	return 0;
}

/* ─── unsafe variants — caller MUST ks_resize() up-front ──────────────────
 * These skip the per-call length-check + realloc and the trailing NUL. The
 * caller is responsible for: (a) ensuring buffer has room, and (b) setting
 * the trailing s->s[s->l] = 0 once at the end. About 5-10× faster than the
 * safe variants because they're branchless and inlinable. */

static inline void kputc_u(int c, kstring_t *s) {
	s->s[s->l++] = (char)c;
}

static inline void kputsn_u(const char *p, int l, kstring_t *s) {
	memcpy(s->s + s->l, p, (size_t)l);
	s->l += l;
}

static inline void kputs_u(const char *p, kstring_t *s) {
	int l = (int)strlen(p);
	memcpy(s->s + s->l, p, (size_t)l);
	s->l += l;
}

static inline void kputuw_u(unsigned c, kstring_t *s) {
	if (c == 0) { s->s[s->l++] = '0'; return; }
	char buf[12];
	int l = 0;
	while (c > 0) { buf[l++] = (char)(c % 10 + '0'); c /= 10; }
	while (l-- > 0) s->s[s->l++] = buf[l];
}

static inline void kputw_u(int c, kstring_t *s) {
	/* Negate in the unsigned domain to avoid signed-overflow UB on INT_MIN. */
	if (c < 0) { s->s[s->l++] = '-'; kputuw_u(0u - (unsigned)c, s); return; }
	kputuw_u((unsigned)c, s);
}

static inline void kputul_u(unsigned long c, kstring_t *s) {
	if (c == 0) { s->s[s->l++] = '0'; return; }
	char buf[24];
	int l = 0;
	while (c > 0) { buf[l++] = (char)(c % 10 + '0'); c /= 10; }
	while (l-- > 0) s->s[s->l++] = buf[l];
}

static inline void kputl_u(long c, kstring_t *s) {
	/* Negate in the unsigned domain to avoid signed-overflow UB on LONG_MIN. */
	if (c < 0) { s->s[s->l++] = '-'; kputul_u(0ul - (unsigned long)c, s); return; }
	kputul_u((unsigned long)c, s);
}

int ksprintf(kstring_t *s, const char *fmt, ...);

#endif
