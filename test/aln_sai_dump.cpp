// Writes the body of a stock-format `bwa aln` .sai file for a FASTQ, using the
// ported search (src/aln) over bwa-mem3's own index.
//
// usage: aln_sai_dump [aln options] <idxbase> <in.fq>      .sai body on stdout
//        aln_sai_dump --gap-opt-size                       sizeof(gap_opt_t)
//        aln_sai_dump --write-fasta <path>                 the shared synthetic reference
//        aln_sai_dump --explain <offset> [aln options] <idxbase> <in.fq>
//
// Options are stock `bwa aln`'s: -n -o -e -i -d -l -k -m -t -M -O -E -q -I -Y -B -L -R -N.
// The body is, for each read in order, `int n_aln` followed by n_aln bwt_aln1_t records,
// which is what stock writes after its 4-byte magic and gap_opt_t header. The header is
// not written: it carries -t, so it differs between thread counts while the body does
// not. test/regression/aln_sai_parity.sh strips stock's header and compares bodies.
//
// --explain takes the 1-based byte offset `cmp -l` reports in a body and prints which
// read and field it falls in, then exits.
//
// This is a test tool; the product has no .sai writer.

#include "FMI_search.h"
#include "aln/aln_bwt.h"
#include "aln/bwaseqio.h"
#include "aln/bwtgap.h"
// The synthetic reference is the unit tests' own generator, included rather than copied: the
// regression scripts' expected outputs depend on its exact sequence, so a second copy (in awk,
// say) could drift from the one the unit tests check.
#include "unit/aln_test_ref.h"

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <unistd.h>

using namespace bwa_aln;

namespace {

// Offsets of the fields within a bwt_aln1_t as written by fwrite.
std::string field_of(size_t off_in_rec)
{
    // bits 0..7 n_mm, 8..15 n_gapo, 16..23 n_gape, 24..43 score, 44..53 n_ins, 54..63 n_del
    // in the first (little-endian) 8 bytes; then k and l.
    char buf[96];
    if (off_in_rec < 8) {
        int lo = (int)off_in_rec * 8, hi = lo + 7;
        const char *names[6] = {"n_mm", "n_gapo", "n_gape", "score", "n_ins", "n_del"};
        const int start[6] = {0, 8, 16, 24, 44, 54}, end[6] = {7, 15, 23, 43, 53, 63};
        std::string s;
        for (int f = 0; f < 6; ++f)
            if (start[f] <= hi && end[f] >= lo) s += std::string(s.empty() ? "" : "/") + names[f];
        snprintf(buf, sizeof buf, "packed counts byte %d (%s)", (int)off_in_rec, s.c_str());
    } else if (off_in_rec < 16) {
        snprintf(buf, sizeof buf, "k byte %d", (int)off_in_rec - 8);
    } else {
        snprintf(buf, sizeof buf, "l byte %d", (int)off_in_rec - 16);
    }
    return buf;
}

// Describes where body offset `off` (0-based) falls, given each read's n_aln.
void explain(size_t off, const bwa_seq_t *seqs, int n, size_t base_read)
{
    size_t pos = 0;
    for (int i = 0; i < n; ++i) {
        size_t len = 4 + (size_t)seqs[i].n_aln * sizeof(bwt_aln1_t);
        if (off < pos + len) {
            size_t o = off - pos;
            if (o < 4) {
                printf("read %zu (%s): n_aln byte %zu\n", base_read + i, seqs[i].name, o);
            } else {
                size_t r = (o - 4) / sizeof(bwt_aln1_t);
                printf("read %zu (%s): alignment %zu of %d, %s\n", base_read + i, seqs[i].name, r,
                       seqs[i].n_aln, field_of((o - 4) % sizeof(bwt_aln1_t)).c_str());
            }
            exit(0);
        }
        pos += len;
    }
}

} // namespace

int main(int argc, char *argv[])
{
    for (int i = 1; i < argc; ++i)
        if (strcmp(argv[i], "--gap-opt-size") == 0) {
            printf("%zu\n", sizeof(gap_opt_t));
            return 0;
        }
    if (argc == 3 && strcmp(argv[1], "--write-fasta") == 0) {
        if (!aln_test::write_test_fasta(argv[2])) { fprintf(stderr, "cannot write %s\n", argv[2]); return 1; }
        return 0;
    }
    long long explain_off = -1;
    if (argc > 2 && strcmp(argv[1], "--explain") == 0) {
        explain_off = atoll(argv[2]);
        if (explain_off < 1) { fprintf(stderr, "--explain wants a 1-based offset\n"); return 1; }
        argv += 2;
        argc -= 2;
    }

    int c, opte = -1;
    gap_opt_t *opt = gap_init_opt();
    while ((c = getopt(argc, argv, "n:o:e:i:d:l:k:LR:m:t:NM:O:E:q:IYB:")) >= 0) {
        switch (c) {
        case 'n':
            if (strstr(optarg, ".")) opt->fnr = atof(optarg), opt->max_diff = -1;
            else opt->max_diff = atoi(optarg), opt->fnr = -1.0;
            break;
        case 'o': opt->max_gapo = atoi(optarg); break;
        case 'e': opte = atoi(optarg); break;
        case 'M': opt->s_mm = atoi(optarg); break;
        case 'O': opt->s_gapo = atoi(optarg); break;
        case 'E': opt->s_gape = atoi(optarg); break;
        case 'd': opt->max_del_occ = atoi(optarg); break;
        case 'i': opt->indel_end_skip = atoi(optarg); break;
        case 'l': opt->seed_len = atoi(optarg); break;
        case 'k': opt->max_seed_diff = atoi(optarg); break;
        case 'm': opt->max_entries = atoi(optarg); break;
        case 't': opt->n_threads = atoi(optarg); break;
        case 'L': opt->mode |= BWA_MODE_LOGGAP; break;
        case 'R': opt->max_top2 = atoi(optarg); break;
        case 'q': opt->trim_qual = atoi(optarg); break;
        case 'N': opt->mode |= BWA_MODE_NONSTOP; opt->max_top2 = 0x7fffffff; break;
        case 'I': opt->mode |= BWA_MODE_IL13; break;
        case 'Y': opt->mode |= BWA_MODE_CFY; break;
        case 'B': opt->mode |= atoi(optarg) << 24; break;
        default: return 1;
        }
    }
    if (opte > 0) {
        opt->max_gape = opte;
        opt->mode &= ~BWA_MODE_GAPE;
    }
    if (optind + 2 > argc) {
        fprintf(stderr, "usage: aln_sai_dump [aln options] <idxbase> <in.fq>\n");
        return 1;
    }

    FMI_search fmi(argv[optind]);
    fmi.load_index(false, 1);
    bwa_aln::bwt_t *bwt = bwt_from_fmi(&fmi);
    bwa_seqio_t *ks = bwa_seq_open(argv[optind + 1]);
    if (!ks) { fprintf(stderr, "cannot open %s\n", argv[optind + 1]); return 1; }

    int n_seqs;
    size_t n_done = 0, body_pos = 0;
    bwa_seq_t *seqs;
    while ((seqs = bwa_read_seq(ks, 0x40000, &n_seqs, opt->mode, opt->trim_qual)) != 0) {
        aln_search_batch(bwt, n_seqs, seqs, opt, opt->n_threads);
        if (explain_off > 0) {
            explain((size_t)explain_off - 1 - body_pos, seqs, n_seqs, n_done);
        } else {
            for (int i = 0; i < n_seqs; ++i) {
                fwrite(&seqs[i].n_aln, 4, 1, stdout);
                if (seqs[i].n_aln) fwrite(seqs[i].aln, sizeof(bwt_aln1_t), seqs[i].n_aln, stdout);
            }
        }
        for (int i = 0; i < n_seqs; ++i) {
            body_pos += 4 + (size_t)seqs[i].n_aln * sizeof(bwt_aln1_t); // aln is freed by bwa_free_read_seq
        }
        n_done += n_seqs;
        bwa_free_read_seq(n_seqs, seqs);
    }
    if (explain_off > 0) { fprintf(stderr, "offset %lld is past the end of the body\n", explain_off); return 1; }
    bwa_seq_close(ks);
    bwt_destroy(bwt);
    free(opt);
    return fflush(stdout) == 0 ? 0 : 1;
}
