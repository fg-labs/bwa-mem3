/* test/framework/temp_path.h
 *
 * Scratch files and directories for the standalone C tests in test/ (C++
 * tests use temp_dir.h, which builds on this). C99, header-only.
 */

#ifndef BWA_TESTS_TEMP_PATH_H
#define BWA_TESTS_TEMP_PATH_H

#include <dirent.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

/* The scratch root: $TMPDIR when set and non-empty, else /tmp. CI points
 * TMPDIR at a writable scratch area, and /tmp is RAM-backed on several
 * distributions. */
static inline const char *bwa_test_temp_root(void)
{
    const char *t = getenv("TMPDIR");
    return (t != NULL && t[0] != '\0') ? t : "/tmp";
}

/* Writes "<root>/<prefix>.XXXXXX" into buf (size n), ready for mkstemp or
 * mkdtemp. Returns 0, or -1 (with a message on stderr) if it does not fit. */
static inline int bwa_test_temp_template(char *buf, size_t n, const char *prefix)
{
    const char *root = bwa_test_temp_root();
    const char *sep = root[strlen(root) - 1] == '/' ? "" : "/";
    const int len = snprintf(buf, n, "%s%s%s.XXXXXX", root, sep, prefix);
    if (len < 0 || (size_t)len >= n) {
        fprintf(stderr, "FAIL: scratch path too long under %s\n", root);
        return -1;
    }
    return 0;
}

/* bwa_test_temp_template into a fresh malloc'd buffer (the caller frees it);
 * aborts if that fails. */
static inline char *bwa_test_temp_alloc(const char *prefix)
{
    const size_t n = strlen(bwa_test_temp_root()) + strlen(prefix) + 16;
    char *buf = (char *)malloc(n);
    if (buf == NULL || bwa_test_temp_template(buf, n, prefix) != 0) abort();
    return buf;
}

/* Removes every entry of `dir` (it must hold only regular files) and then
 * `dir` itself. Returns rmdir's result: 0, or -1 with errno set. */
static inline int bwa_test_remove_flat_dir(const char *dir)
{
    DIR *d = opendir(dir);
    if (d != NULL) {
        struct dirent *e;
        char path[4096];
        while ((e = readdir(d)) != NULL) {
            if (strcmp(e->d_name, ".") == 0 || strcmp(e->d_name, "..") == 0) continue;
            if (snprintf(path, sizeof path, "%s/%s", dir, e->d_name) < (int)sizeof path) unlink(path);
        }
        closedir(d);
    }
    return rmdir(dir);
}

#endif /* BWA_TESTS_TEMP_PATH_H */
