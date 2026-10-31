/* test/framework/temp_path.h
 *
 * Scratch files and directories under $TMPDIR, as a C99 header-only API, so
 * that the standalone C tests in test/ and any C++ test that wants a plain
 * path in a char buffer (e.g. one it hands to a child process) share one
 * template and cleanup. temp_dir.h's TempDir and anon_temp_fd build on it.
 */

#ifndef BWA_TESTS_TEMP_PATH_H
#define BWA_TESTS_TEMP_PATH_H

#include <dirent.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

/* Buffer size for a scratch path (the template, or a file in a scratch dir). */
#define BWA_TEST_PATH_MAX 4096

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

/* Creates a fresh scratch file and leaves its path in buf (size n). Returns
 * the open fd, or -1 with a message on stderr. */
static inline int bwa_test_mkstemp(char *buf, size_t n, const char *prefix)
{
    if (bwa_test_temp_template(buf, n, prefix) != 0) return -1;
    const int fd = mkstemp(buf);
    if (fd < 0) fprintf(stderr, "FAIL: mkstemp(%s): %s\n", buf, strerror(errno));
    return fd;
}

/* Creates a fresh scratch directory and leaves its path in buf (size n).
 * Returns 0, or -1 with a message on stderr. */
static inline int bwa_test_mkdtemp(char *buf, size_t n, const char *prefix)
{
    if (bwa_test_temp_template(buf, n, prefix) != 0) return -1;
    if (mkdtemp(buf) == NULL) {
        fprintf(stderr, "FAIL: mkdtemp(%s): %s\n", buf, strerror(errno));
        return -1;
    }
    return 0;
}

/* Removes every entry of `dir` (it must hold only regular files) and then
 * `dir` itself. Returns rmdir's result: 0, or -1 with errno set. */
static inline int bwa_test_remove_flat_dir(const char *dir)
{
    DIR *d = opendir(dir);
    if (d != NULL) {
        struct dirent *e;
        char path[BWA_TEST_PATH_MAX];
        while ((e = readdir(d)) != NULL) {
            if (strcmp(e->d_name, ".") == 0 || strcmp(e->d_name, "..") == 0) continue;
            if (snprintf(path, sizeof path, "%s/%s", dir, e->d_name) < (int)sizeof path) unlink(path);
        }
        closedir(d);
    }
    return rmdir(dir);
}

#endif /* BWA_TESTS_TEMP_PATH_H */
