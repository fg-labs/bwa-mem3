#!/usr/bin/env bash
# test/regression/aln_sai_parity.sh
#
# Regression: the `aln` port's backtracking search returns the same alignments as
# stock `bwa aln` (0.7.19), byte for byte.
#
# A .sai file is a 4-byte magic and a 64-byte gap_opt_t header (which records -t,
# so it differs between thread counts) followed by the body: for each read, its
# alignment count and alignment records. aln_sai_dump writes only the body, from
# the ported search over a bwa-mem3 index; the script strips the header from stock's
# output and compares the bodies with cmp. Both sides read the same FASTQ, so the
# reader's trimming and barcode handling are exercised too.
#
# The reference is phiX plus the synthetic contigs the unit tests use (an N run, a
# reverse-complement copy, a tandem repeat), indexed once by stock bwa and once by
# bwa-mem3. Reads are generated here by a seeded awk program (integer arithmetic
# only, so every awk gives the same reads) from that reference: exact copies,
# substitutions, 1-bp insertions and deletions, ambiguous bases, reads from the
# repeat and from the reverse-complement copy, on both strands.
#
# Each case runs with -t 1 and -t 4. On a mismatch the first differing byte is
# reported with the read and field it falls in (aln_sai_dump --explain).
#
# Inputs (env vars):
#   BWA_MEM3          -- path to the bwa-mem3 binary under test
#   BWA_ORACLE        -- path to a stock bwa 0.7.19 binary
#                        (test/tools/build_bwa_oracle.sh); unset => SKIP
#   FIXTURES          -- directory holding phix.fa (default: test/fixtures)
#   ALN_SAI_DUMP      -- path to aln_sai_dump (default: repo root)

set -Eeuo pipefail
# shellcheck disable=SC2154  # rc is assigned by `rc=$?` at the start of the trap body
trap 'rc=$?; echo "FAIL: unexpected command failure at line $LINENO (exit $rc)" >&2; exit "$rc"' ERR

if [ -z "${BWA_ORACLE:-}" ]; then
    echo "SKIP: BWA_ORACLE unset"
    exit 0
fi
if [ -z "${BWA_MEM3:-}" ]; then
    echo "FAIL: BWA_MEM3 must be set" >&2
    exit 1
fi

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
fixtures="${FIXTURES:-$repo_root/test/fixtures}"
dump="${ALN_SAI_DUMP:-$repo_root/aln_sai_dump}"

work="$(mktemp -d "${TMPDIR:-/tmp}/aln_sai_parity.XXXXXX")"
cleanup() { rm -rf "$work"; }
trap cleanup EXIT

# Reference: phiX followed by the synthetic contigs.
"$dump" --write-fasta "$work/synthetic.fa"
cat "$fixtures/phix.fa" "$work/synthetic.fa" > "$work/ref.fa"
mkdir -p "$work/stock" "$work/mem3"
cp "$work/ref.fa" "$work/stock/ref.fa"
cp "$work/ref.fa" "$work/mem3/ref.fa"
"$BWA_ORACLE" index "$work/stock/ref.fa" > "$work/stock.index.log" 2>&1
"$BWA_MEM3" index "$work/mem3/ref.fa" > "$work/mem3.index.log" 2>&1

# Seeded read generator. usage: gen_reads <seed> <n> <lens> <tail>
#   lens  -- comma-separated read lengths, or MIN-MAX for a uniform range
#   tail  -- 1 gives a fifth of the reads a low-quality 30-base tail
# Park-Miller generator (16807 mod 2^31-1): products stay below 2^46, so the
# arithmetic is exact in any awk's doubles.
cat > "$work/gen.awk" << 'EOF'
function rnd() { st = (st * 16807) % 2147483647; return st }
function ri(n) { return rnd() % n }
function rc(s,   i, r, c) {
    r = ""
    for (i = length(s); i >= 1; i--) {
        c = substr(s, i, 1)
        r = r (c == "A" ? "T" : c == "C" ? "G" : c == "G" ? "C" : c == "T" ? "A" : c)
    }
    return r
}
function other(c,   b, k) {
    b = "ACGT"
    do { k = substr(b, ri(4) + 1, 1) } while (k == c)
    return k
}
/^>/ { nc++; name[nc] = substr($1, 2); next }
{ seq[nc] = seq[nc] toupper($0) }
END {
    st = seed
    for (i = 0; i < 20; i++) rnd()
    nl = split(lens, lv, ",")
    for (r = 0; r < n; r++) {
        if (nl == 1 && index(lv[1], "-")) {
            split(lv[1], mm, "-"); len = mm[1] + ri(mm[2] - mm[1] + 1)
        } else len = lv[ri(nl) + 1]
        # contig: a tenth of the reads come from the tandem repeat, another tenth
        # from the reverse-complement copy, the rest uniformly over contigs
        x = ri(10)
        if (x == 0) c = nc; else if (x == 1) c = nc - 1; else c = ri(nc) + 1
        s = seq[c]
        L = length(s)
        if (c == nc - 1 && x == 1) { lo = 1; span = 300 } else { lo = 1; span = L }
        if (len > span) len = span
        pos = lo + ri(span - len + 1)
        t = substr(s, pos, len)
        if (ri(2)) t = rc(t)
        # edits
        e = ri(100)
        p = ri(len > 0 ? len : 1) + 1
        if (e < 25) t = substr(t, 1, p - 1) other(substr(t, p, 1)) substr(t, p + 1)
        else if (e < 35) {
            t = substr(t, 1, p - 1) other(substr(t, p, 1)) substr(t, p + 1)
            p = ri(len) + 1
            t = substr(t, 1, p - 1) other(substr(t, p, 1)) substr(t, p + 1)
        }
        else if (e < 50) t = substr(t, 1, p - 1) substr(t, p + 1)
        else if (e < 60) t = substr(t, 1, p - 1) substr("ACGT", ri(4) + 1, 1) substr(t, p)
        else if (e < 65) t = substr(t, 1, p - 1) "N" substr(t, p + 1)
        n2 = length(t)
        q = ""
        for (i = 1; i <= n2; i++) q = q ((tail && r % 5 == 0 && i > n2 - 30) ? "#" : "I")
        printf "@r%d\n%s\n+\n%s\n", r, t, q
    }
}
EOF
gen_reads() {
    awk -v seed="$1" -v n="$2" -v lens="$3" -v tail="$4" -f "$work/gen.awk" "$work/ref.fa"
}

hdr=$((4 + $("$dump" --gap-opt-size)))
n_cfg=0

# run_case <name> <seed> <n> <lens> <tail> [aln options...]
run_case() {
    local name="$1" seed="$2" n="$3" lens="$4" tail="$5"
    shift 5
    gen_reads "$seed" "$n" "$lens" "$tail" > "$work/$name.fq"
    local t
    for t in 1 4; do
        n_cfg=$((n_cfg + 1))
        local tag="$name -t $t $*"
        "$BWA_ORACLE" aln "$@" -t "$t" "$work/stock/ref.fa" "$work/$name.fq" \
            > "$work/stock.sai" 2> "$work/stock.aln.log"
        tail -c +$((hdr + 1)) "$work/stock.sai" > "$work/stock.body"
        "$dump" "$@" -t "$t" "$work/mem3/ref.fa" "$work/$name.fq" > "$work/mem3.body" 2> "$work/mem3.dump.log"
        if ! cmp -s "$work/stock.body" "$work/mem3.body"; then
            local off
            off="$(cmp -l "$work/stock.body" "$work/mem3.body" 2> /dev/null | head -n 1 | awk '{ print $1 }' || true)"
            echo "FAIL: aln_sai_parity: $tag differs from stock at body byte ${off:-?}" >&2
            if [ -n "$off" ]; then
                "$dump" --explain "$off" "$@" -t "$t" "$work/mem3/ref.fa" "$work/$name.fq" > "$work/explain.txt" 2> /dev/null || true
                cat "$work/explain.txt" >&2
            fi
            cmp -l "$work/stock.body" "$work/mem3.body" 2> /dev/null | head -5 >&2 || true
            exit 1
        fi
        # A body of all-zero counts would pass trivially: require real alignments.
        if [ "$(wc -c < "$work/stock.body")" -le $((4 * n)) ]; then
            echo "FAIL: aln_sai_parity: $tag found no alignments at all" >&2
            exit 1
        fi
        echo "aln_sai_parity: ok $tag ($(wc -c < "$work/stock.body") body bytes)"
    done
}

run_case default36 101 2000 36 0
run_case adna 102 2000 20-40 0 -l 1024 -n 0.04 -o 2
run_case mixed 103 1000 20,100 0 -n 0.04 -o 2
run_case qtrim 104 500 150 1 -q 15
run_case barcode 105 600 3,5,6,36,50 0 -B 5
run_case nonstop 106 600 36 0 -N
run_case penalties 107 1000 30-60 0 -M 4 -R 10
run_case seedoff 108 1000 36 0 -k 0 -l 16

echo "PASS: aln_sai_parity ($n_cfg configurations)"
