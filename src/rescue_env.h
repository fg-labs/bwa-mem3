/* Shared readers for the BWA3_RESCUE_* knobs (kswv.cpp, rescue_prune.h, rescue_band.h).
 *
 * rescue_env_int: an integer tuning knob, read once per process into a function-local static by
 * its caller. It must be a non-negative decimal integer. An unset or empty knob silently takes the
 * default; anything else (a sign, leading or trailing characters, overflow) is reported to stderr
 * and the default used, so a typo cannot silently set a knob.
 *
 * rescue_env_on: an on/off toggle, default on; a value starting with '0' turns it off and
 * anything else leaves it on. Not reported and not cached: kswv reads its toggles on every call so
 * a unit test can flip them in-process, and a per-call report would repeat once per batch;
 * callers that want the value once wrap it in a function-local static. */
#ifndef BWA_MEM3_RESCUE_ENV_H
#define BWA_MEM3_RESCUE_ENV_H

#include <errno.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>

static inline int rescue_env_int(const char *name, int dflt)
{
    const char *e = getenv(name);
    if (!e || !*e) return dflt;
    char *end = NULL;
    errno = 0;
    /* strtol skips leading whitespace and takes a sign, so require a digit first. */
    const long x = (e[0] >= '0' && e[0] <= '9') ? strtol(e, &end, 10) : -1;
    if (x < 0 || errno || *end || x > INT_MAX) {
        fprintf(stderr, "ERROR: %s=\"%s\" is not a non-negative integer; using %d.\n", name, e, dflt);
        return dflt;
    }
    return (int)x;
}

static inline bool rescue_env_on(const char *name)
{
    const char *e = getenv(name);
    return !e || e[0] != '0';
}

#endif
