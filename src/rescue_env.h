/* Shared reader for the integer BWA3_RESCUE_* tuning knobs (rescue_prune.h, rescue_band.h).
 *
 * Every knob is read once per process into a function-local static by its caller. A knob must be
 * a non-negative decimal integer. An unset or empty knob silently takes the default; anything
 * else (a sign, leading or trailing characters, overflow) is reported to stderr and the default
 * used, so a typo cannot silently set a knob. */
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

#endif
