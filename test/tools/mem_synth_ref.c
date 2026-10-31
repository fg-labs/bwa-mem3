// test/tools/mem_synth_ref.c
//
// Generator and checker for test/fixtures/mem_synth.fa, the synthetic
// reference the `mem` CLI regression scripts (test/regression/mem_*.sh) align
// against. The committed FASTA is what the tests read; this program is how it
// was made and how it is verified:
//
//     ./test/mem_synth_ref > test/fixtures/mem_synth.fa      # regenerate
//     ./test/mem_synth_ref --check test/fixtures/mem_synth.fa
//
// --check confirms that the committed file is byte-identical to the output,
// that the planted features below sit where this header says, and that every
// 19-mer outside them is unique in the reference on both strands. Every
// expected coordinate, MAPQ and XA in the scripts is worked out from those
// three facts, so a change here must keep --check passing and then revisit
// the scripts' expected values.
//
// Bases come from a Park-Miller generator (x = 48271 * x mod 2^31-1, seed
// 20261024); base i of the stream is "ACGT"[x / 2^29] after the i-th step. The
// contigs are drawn from one stream in order, then the features are written
// over them. Coordinates are 1-based and inclusive, as in SAM.
//
//   ctg1  5000 bp  unique random sequence, except:
//                  2001-2300  block R (300 bp; the original copy)
//                  3001-3500  block S (500 bp; the original copy)
//   ctg2  4000 bp  1001-1300  exact copy of R, same strand
//                  2001-2500  copy of S with a substitution at block offsets
//                             10, 30, 50, ..., 490 (25 sites, 5%): each base
//                             becomes the next one in ACGT order (A>C>G>T>A)
//   ctg3  3000 bp   501-620   (AC)x60 dinucleotide repeat
//                  1501-1540  40 N
//                  2201-2250  poly-A (50 A)
//   ctg4   600 bp   101-400   reverse complement of R
//
// So a read inside R has three equally good hits (ctg1 +, ctg2 +, ctg4 -), and
// a read inside S has a best hit on ctg1 and a 5%-diverged second hit on ctg2.

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define N_CONTIGS 4
#define KMER 19
#define LINE_WIDTH 60
#define MAX_UNIQ_REPORTS 20

static const char *const kNames[N_CONTIGS] = {"ctg1", "ctg2", "ctg3", "ctg4"};
static const int kLens[N_CONTIGS] = {5000, 4000, 3000, 600};

// The planted features, as {contig index, first, last} (1-based, inclusive).
static const int kFeatures[][3] = {
    {0, 2001, 2300}, {0, 3001, 3500}, {1, 1001, 1300}, {1, 2001, 2500},
    {2, 501, 620},   {2, 1501, 1540}, {2, 2201, 2250}, {3, 101, 400},
};
#define N_FEATURES ((int)(sizeof(kFeatures) / sizeof(kFeatures[0])))

static char next_base(char c) {
    return c == 'A' ? 'C' : c == 'C' ? 'G' : c == 'G' ? 'T' : 'A';
}

static char comp(char c) {
    return c == 'A' ? 'T' : c == 'C' ? 'G' : c == 'G' ? 'C' : c == 'T' ? 'A' : 'N';
}

static void revcomp(const char *s, int n, char *out) {
    for (int i = 0; i < n; i++) out[i] = comp(s[n - 1 - i]);
}

// contigs[i] gets kLens[i] + 1 bytes, NUL-terminated.
static void build(char *contigs[N_CONTIGS]) {
    uint64_t x = 20261024u;
    for (int c = 0; c < N_CONTIGS; c++) {
        contigs[c] = malloc((size_t)kLens[c] + 1);
        if (contigs[c] == NULL) {
            perror("malloc");
            exit(2);
        }
        for (int i = 0; i < kLens[c]; i++) {
            x = (x * 48271u) % 2147483647u;
            contigs[c][i] = "ACGT"[x / 536870912u];
        }
        contigs[c][kLens[c]] = '\0';
    }
    const char *R = contigs[0] + 2000, *S = contigs[0] + 3000;
    memcpy(contigs[1] + 1000, R, 300);
    for (int i = 1; i <= 500; i++)
        contigs[1][2000 + i - 1] = i % 20 == 10 ? next_base(S[i - 1]) : S[i - 1];
    for (int i = 0; i < 120; i++) contigs[2][500 + i] = i % 2 == 0 ? 'A' : 'C';
    memset(contigs[2] + 1500, 'N', 40);
    memset(contigs[2] + 2200, 'A', 50);
    revcomp(R, 300, contigs[3] + 100);
}

// The FASTA text, LINE_WIDTH bases per line; *len receives its length.
static char *fasta(char *const contigs[N_CONTIGS], size_t *len) {
    size_t cap = 0;
    for (int c = 0; c < N_CONTIGS; c++)
        cap += strlen(kNames[c]) + 2 + (size_t)kLens[c] + (size_t)kLens[c] / LINE_WIDTH + 1;
    char *out = malloc(cap + 1);
    if (out == NULL) {
        perror("malloc");
        exit(2);
    }
    size_t n = 0;
    for (int c = 0; c < N_CONTIGS; c++) {
        n += (size_t)snprintf(out + n, cap + 1 - n, ">%s\n", kNames[c]);
        for (int i = 0; i < kLens[c]; i += LINE_WIDTH) {
            const int w = kLens[c] - i < LINE_WIDTH ? kLens[c] - i : LINE_WIDTH;
            memcpy(out + n, contigs[c] + i, (size_t)w);
            n += (size_t)w;
            out[n++] = '\n';
        }
    }
    out[n] = '\0';
    *len = n;
    return out;
}

static char *slurp(const char *path, size_t *len) {
    FILE *fp = fopen(path, "rb");
    if (fp == NULL) return NULL;
    size_t cap = 1 << 14, n = 0;
    char *buf = malloc(cap + 1);
    for (;;) {
        if (buf == NULL) break;
        n += fread(buf + n, 1, cap - n, fp);
        if (n < cap) break;
        cap *= 2;
        char *grown = realloc(buf, cap + 1);
        if (grown == NULL) free(buf);
        buf = grown;
    }
    const int err = ferror(fp);
    fclose(fp);
    if (buf == NULL || err) {
        free(buf);
        return NULL;
    }
    buf[n] = '\0';
    *len = n;
    return buf;
}

// Parses `text` into the N_CONTIGS contigs named in kNames, in that order.
// Returns 0 on success; on failure reports and returns 1, with every
// contigs[i] either NULL or allocated.
static int parse(const char *text, char *contigs[N_CONTIGS], int lens[N_CONTIGS]) {
    int c = -1;
    for (int i = 0; i < N_CONTIGS; i++) {
        contigs[i] = NULL;
        lens[i] = 0;
    }
    for (const char *p = text; *p != '\0';) {
        const char *eol = strchr(p, '\n');
        const size_t n = eol != NULL ? (size_t)(eol - p) : strlen(p);
        if (n > 0 && p[0] == '>') {
            if (++c >= N_CONTIGS || n - 1 != strlen(kNames[c]) ||
                strncmp(p + 1, kNames[c], n - 1) != 0) {
                fprintf(stderr, "FAIL: header %d is '%.*s', expected '>%s'\n", c + 1,
                        (int)n, p, c < N_CONTIGS ? kNames[c] : "(none)");
                return 1;
            }
            contigs[c] = malloc(strlen(p) + 1);
            if (contigs[c] == NULL) {
                perror("malloc");
                exit(2);
            }
        } else if (c >= 0) {
            memcpy(contigs[c] + lens[c], p, n);
            lens[c] += (int)n;
            contigs[c][lens[c]] = '\0';
        }
        p += n + (eol != NULL);
    }
    if (c != N_CONTIGS - 1) {
        fprintf(stderr, "FAIL: %d contigs, expected %d\n", c + 1, N_CONTIGS);
        return 1;
    }
    return 0;
}

// The planted features of the header, read back from the parsed file.
static int check_features(char *const contigs[N_CONTIGS], const int lens[N_CONTIGS]) {
    int fails = 0;
    for (int c = 0; c < N_CONTIGS; c++) {
        if (lens[c] != kLens[c]) {
            fprintf(stderr, "FAIL: %s is %d bp, expected %d\n", kNames[c], lens[c], kLens[c]);
            fails++;
        }
    }
    if (fails) return fails;
    const char *R = contigs[0] + 2000, *S = contigs[0] + 3000;
    char rc[300];
    revcomp(R, 300, rc);
    if (memcmp(contigs[1] + 1000, R, 300) != 0) {
        fprintf(stderr, "FAIL: ctg2:1001-1300 is not a copy of R (ctg1:2001-2300)\n");
        fails++;
    }
    if (memcmp(contigs[3] + 100, rc, 300) != 0) {
        fprintf(stderr, "FAIL: ctg4:101-400 is not the reverse complement of R\n");
        fails++;
    }
    for (int i = 1; i <= 500; i++) {
        const char want = i % 20 == 10 ? next_base(S[i - 1]) : S[i - 1];
        if (contigs[1][2000 + i - 1] != want) {
            fprintf(stderr, "FAIL: ctg2:%d is '%c', expected '%c' (S block offset %d)\n",
                    2000 + i, contigs[1][2000 + i - 1], want, i);
            fails++;
        }
    }
    for (int i = 0; i < 120; i++) {
        if (contigs[2][500 + i] != (i % 2 == 0 ? 'A' : 'C')) {
            fprintf(stderr, "FAIL: ctg3:501-620 is not (AC)x60\n");
            fails++;
            break;
        }
    }
    for (int i = 0; i < 40; i++) {
        if (contigs[2][1500 + i] != 'N') {
            fprintf(stderr, "FAIL: ctg3:1501-1540 is not 40 N\n");
            fails++;
            break;
        }
    }
    for (int i = 0; i < 50; i++) {
        if (contigs[2][2200 + i] != 'A') {
            fprintf(stderr, "FAIL: ctg3:2201-2250 is not 50 A\n");
            fails++;
            break;
        }
    }
    return fails;
}

static int code(char b) {
    return b == 'A' ? 0 : b == 'C' ? 1 : b == 'G' ? 2 : b == 'T' ? 3 : -1;
}

// 2 bits per base, first base highest; -1 when the k-mer holds a non-ACGT base.
static int64_t encode(const char *s) {
    int64_t v = 0;
    for (int i = 0; i < KMER; i++) {
        const int b = code(s[i]);
        if (b < 0) return -1;
        v = v << 2 | b;
    }
    return v;
}

static int64_t encode_rc(int64_t v) {
    int64_t r = 0;
    for (int i = 0; i < KMER; i++, v >>= 2) r = r << 2 | (3 - (v & 3));
    return r;
}

static int cmp_i64(const void *a, const void *b) {
    const int64_t x = *(const int64_t *)a, y = *(const int64_t *)b;
    return (x > y) - (x < y);
}

// Occurrences of `key` in the sorted array.
static int count(const int64_t *v, size_t n, int64_t key) {
    size_t lo = 0, hi = n;
    while (lo < hi) {
        const size_t mid = lo + (hi - lo) / 2;
        if (v[mid] < key) lo = mid + 1;
        else hi = mid;
    }
    int k = 0;
    while (lo + (size_t)k < n && v[lo + (size_t)k] == key) k++;
    return k;
}

// Does [first, first + KMER) (1-based) of contig c overlap a planted feature?
static int in_feature(int c, int first) {
    const int last = first + KMER - 1;
    for (int f = 0; f < N_FEATURES; f++)
        if (kFeatures[f][0] == c && first <= kFeatures[f][2] && last >= kFeatures[f][1])
            return 1;
    return 0;
}

// Every 19-mer outside the features occurs once on the forward strand of the
// whole reference and never as a reverse complement.
static int check_unique(char *const contigs[N_CONTIGS], const int lens[N_CONTIGS]) {
    size_t total = 0;
    for (int c = 0; c < N_CONTIGS; c++) total += (size_t)lens[c];
    int64_t *all = malloc(total * sizeof(*all));
    if (all == NULL) {
        perror("malloc");
        exit(2);
    }
    size_t n = 0;
    for (int c = 0; c < N_CONTIGS; c++)
        for (int i = 0; i + KMER <= lens[c]; i++) {
            const int64_t v = encode(contigs[c] + i);
            if (v >= 0) all[n++] = v;
        }
    qsort(all, n, sizeof(*all), cmp_i64);
    int fails = 0;
    for (int c = 0; c < N_CONTIGS; c++)
        for (int i = 0; i + KMER <= lens[c]; i++) {
            const int64_t v = encode(contigs[c] + i);
            if (v < 0 || in_feature(c, i + 1)) continue;
            const int fwd = count(all, n, v), rev = count(all, n, encode_rc(v));
            if (fwd == 1 && rev == 0) continue;
            if (++fails <= MAX_UNIQ_REPORTS)
                fprintf(stderr,
                        "FAIL: 19-mer at %s:%d occurs %d times forward, %d reverse-complement\n",
                        kNames[c], i + 1, fwd, rev);
        }
    if (fails > MAX_UNIQ_REPORTS)
        fprintf(stderr, "FAIL: ... %d non-unique 19-mers in all\n", fails);
    free(all);
    return fails;
}

static int check(const char *path, const char *want, size_t want_len) {
    size_t len = 0;
    char *text = slurp(path, &len);
    if (text == NULL) {
        fprintf(stderr, "FAIL: cannot read %s\n", path);
        return 1;
    }
    int fails = 0;
    if (len != want_len || memcmp(text, want, len) != 0) {
        size_t i = 0;
        while (i < len && i < want_len && text[i] == want[i]) i++;
        fprintf(stderr, "FAIL: %s differs from mem_synth_ref's output at byte %zu (%zu vs %zu bytes)\n",
                path, i, len, want_len);
        fails++;
    }
    char *contigs[N_CONTIGS];
    int lens[N_CONTIGS];
    if (parse(text, contigs, lens) != 0) {
        fails++;
    } else {
        const int feature_fails = check_features(contigs, lens);
        fails += feature_fails;
        if (feature_fails == 0) fails += check_unique(contigs, lens);
    }
    for (int c = 0; c < N_CONTIGS; c++) free(contigs[c]);
    free(text);
    if (fails == 0) printf("ok: %s matches mem_synth_ref\n", path);
    return fails != 0;
}

int main(int argc, char **argv) {
    if (!(argc == 1 || (argc == 3 && strcmp(argv[1], "--check") == 0))) {
        fprintf(stderr, "usage: %s [--check <mem_synth.fa>]\n", argv[0]);
        return 2;
    }
    char *contigs[N_CONTIGS];
    build(contigs);
    size_t len = 0;
    char *out = fasta(contigs, &len);
    int ret = 0;
    if (argc == 1) {
        if (fwrite(out, 1, len, stdout) != len || fflush(stdout) != 0) {
            perror("write");
            ret = 2;
        }
    } else {
        ret = check(argv[2], out, len);
    }
    free(out);
    for (int c = 0; c < N_CONTIGS; c++) free(contigs[c]);
    return ret;
}
