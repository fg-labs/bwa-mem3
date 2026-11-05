/* Test-only stat interposition for input_fail_fast.sh. EOVERFLOW can mean an
 * inode/size/block count is not representable in struct stat, even though
 * open/read work. Fail only the generated fixture's metadata lookup; leave
 * the real input open and reader untouched. Never linked into bwa-mem3. */
/* Older glibc headers expose RTLD_NEXT only with GNU extensions enabled. */
#if defined(__linux__) && !defined(_GNU_SOURCE)
#define _GNU_SOURCE
#endif
#include <errno.h>
#include <stdio.h>
#include <string.h>
#include <sys/stat.h>

#ifdef INPUT_STAT_PROBE
#include <fcntl.h>
#include <unistd.h>

/* Check the compatibility premise independently of main_mem's preflight. */
int main(int argc, char **argv)
{
    struct stat st;
    char buf[4096];
    ssize_t n;
    if (argc != 2 || stat(argv[1], &st) == 0 || errno != EOVERFLOW) return 1;
    int fd = open(argv[1], O_RDONLY);
    if (fd < 0) return 2;
    while ((n = read(fd, buf, sizeof buf)) > 0) {
        if (fwrite(buf, 1, (size_t)n, stdout) != (size_t)n) return 3;
    }
    return n < 0 || close(fd) != 0 || fflush(stdout) != 0;
}
#else
#ifndef __APPLE__
#include <dlfcn.h>
#endif

static int input_stat(const char *path, struct stat *st)
{
    const char *name = strrchr(path, '/');
    name = name ? name + 1 : path;
    if (strcmp(name, "stat-overflow.fa") == 0) {
        fprintf(stderr, "test stat: EOVERFLOW\n");
        errno = EOVERFLOW;
        return -1;
    }
#ifdef __APPLE__
    return stat(path, st);
#else
    int (*real_stat)(const char *, struct stat *) = dlsym(RTLD_NEXT, "stat");
    if (!real_stat) {
        errno = ENOSYS;
        return -1;
    }
    return real_stat(path, st);
#endif
}

#ifdef __APPLE__
static const struct {
    int (*replacement)(const char *, struct stat *);
    int (*original)(const char *, struct stat *);
} stat_interpose __attribute__((used, section("__DATA,__interpose"))) = {
    input_stat, stat
};
#else
int stat(const char *path, struct stat *st)
{
    return input_stat(path, st);
}
#endif
#endif
