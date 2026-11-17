#!/usr/bin/env bash
# test/regression/compact_occ_byte_identity.sh
#
# Regression: `mem --compact-occ` (and BWA3_COMPACT_OCC=1) holds the FM-index
# occurrence checkpoints in the 32-byte compact layout (CP_OCC_PLANES) instead of
# the index file's 64-byte records. Seeding, the seed jump table and the SA
# resolver read the checkpoints through a separate instantiation per layout, so
# the alignment records must be byte-identical to the default layout's. This
# pins that contract end to end, on paired-end, single-end and long reads, at
# -t 1 and -t 4, and under --compat=bwa-mem2 (whose SA resolver drops the
# sentinel walk offset, the one place the two layouts find the sentinel row
# differently).
#
# The reference has several contigs, so the sentinel row of the BWT lands
# inside a checkpoint block, and a contig made of a repeated segment, so seeds
# have multi-row SA intervals. Fixtures are generated here, not committed
# (repo convention).
#
# Also checked, so the comparisons cannot be vacuous: every compact run reports
# the compact layout on stderr, naming its source (--compact-occ or
# BWA3_COMPACT_OCC), and no default run does; BWA3_COMPACT_OCC=0 leaves the
# default layout, and a malformed value is reported and ignored; only @PG
# differs (it records the flag). Further legs: --meth (its own seed index
# through the same compact load), --huge-pages (on a host without a 1 GB pool
# it falls back to the default allocation; the mapping itself is covered by the
# unit tests' compact leg), and `index` with BWA3_COMPACT_OCC=1 exported writes
# the same files (the library never reads the variable).
#
# Inputs (env vars):
#   BWA_MEM3                 — path to the bwa-mem3 binary under test
#   COMPACT_OCC_SYN_FA       — the 1 Mb synthetic reference (test/fixtures/synthetic_1mb.fa)
#   COMPACT_OCC_PHIX_FA      — phix reference FASTA (test/fixtures/phix.fa)
#   COMPACT_OCC_WORK_DIR     — scratch directory; each run builds its fixtures in a
#                              fresh subdirectory of it (removed on exit), so a rerun
#                              never meets the index files of an earlier run, which
#                              `index` refuses to replace without --force

set -euo pipefail
for v in BWA_MEM3 COMPACT_OCC_SYN_FA COMPACT_OCC_PHIX_FA COMPACT_OCC_WORK_DIR; do
    if [[ -z ${!v:-} ]]; then
        echo "FAIL: $v must be set" >&2
        exit 1
    fi
done

HERE="$(cd "$(dirname "$0")" && pwd)"
# shellcheck source=test/framework/dna.sh
source "$HERE/../framework/dna.sh"
fail() {
    echo "FAIL: $*" >&2
    exit 1
}
ok() { echo "PASS: $*"; }

mkdir -p "$COMPACT_OCC_WORK_DIR" || fail "create work directory"
W=$(mktemp -d "$COMPACT_OCC_WORK_DIR/run.XXXXXX") || fail "create run directory"
trap 'rm -rf "$W"' EXIT

# Reference: the synthetic sequence cut into five 200 kb contigs, phiX, and a
# contig of one 3 kb segment repeated three times.
# shellcheck disable=SC2016 # $0 is awk's record, not a shell expansion
awk '/^>/ { next } { seq = seq $0 }
     END {
         for (i = 0; i < 5; i++) {
             printf ">s%d\n", i + 1
             s = substr(seq, i * 200000 + 1, 200000)
             for (j = 1; j <= length(s); j += 60) print substr(s, j, 60)
         }
         r = substr(seq, 400001, 3000)
         r = r r r
         print ">rep"
         for (j = 1; j <= length(r); j += 60) print substr(r, j, 60)
     }' "$COMPACT_OCC_SYN_FA" > "$W/ref.fa" || fail "build reference"
sed 's/^>.*/>phix/' "$COMPACT_OCC_PHIX_FA" >> "$W/ref.fa" || fail "append phix"
"$BWA_MEM3" index "$W/ref.fa" > /dev/null 2>&1 || fail "index"

# Reads: windows with ~1% substitutions and ~0.1% 1 bp indels from each contig.
# 400 bp windows become pairs (mate 1 = the first 150 bp, mate 2 = the reverse
# complement of the last 150 bp); 150 bp windows are single-end reads; 2-3 kb
# windows are long reads.
: > "$W/w400.fq"
: > "$W/se.fq"
: > "$W/long.fq"
seed=11
for c in s1 s2 s3 s4 s5 rep phix; do
    awk -v NAME="$c" -v SEED=$seed -v NREADS=300 -v LMIN=400 -v LMAX=400 \
        -f "$HERE/../fixtures/make_long_reads.awk" "$W/ref.fa" >> "$W/w400.fq" \
        || fail "generate pair windows ($c)"
    awk -v NAME="$c" -v SEED=$((seed + 1)) -v NREADS=300 -v LMIN=150 -v LMAX=150 \
        -f "$HERE/../fixtures/make_long_reads.awk" "$W/ref.fa" >> "$W/se.fq" \
        || fail "generate single-end reads ($c)"
    awk -v NAME="$c" -v SEED=$((seed + 2)) -v NREADS=4 -v LMIN=2000 -v LMAX=3000 \
        -f "$HERE/../fixtures/make_long_reads.awk" "$W/ref.fa" >> "$W/long.fq" \
        || fail "generate long reads ($c)"
    seed=$((seed + 3))
done
# shellcheck disable=SC2016 # $0 is awk's record, not a shell expansion
awk -v out1="$W/pe_1.fq" -v out2="$W/pe_2.fq" '
    function revcomp(s,    j, c, r) {
        r = ""
        for (j = length(s); j >= 1; j--) {
            c = substr(s, j, 1)
            r = r ((c == "A") ? "T" : (c == "C") ? "G" : (c == "G") ? "C" : (c == "T") ? "A" : c)
        }
        return r
    }
    NR % 4 == 1 { name = substr($1, 2) }
    NR % 4 == 2 {
        q = sprintf("%150s", ""); gsub(/ /, "I", q)
        printf "@%s/1\n%s\n+\n%s\n", name, substr($0, 1, 150), q > out1
        printf "@%s/2\n%s\n+\n%s\n", name, revcomp(substr($0, length($0) - 149)), q > out2
    }' "$W/w400.fq" || fail "split pairs"
# Exact 100 bp reads from the first bases of contig s1, on both strands: their
# SA rows resolve by LF walks that reach the sentinel row (whose offset
# --compat=bwa-mem2 drops), the row the compact layout finds by index.
# shellcheck disable=SC2016 # $0 is awk's record, not a shell expansion
awk '/^>/ { n++; next } n == 1 { seq = seq $0 }
     END { for (k = 0; k < 24; k++) print substr(seq, k + 1, 100) }' "$W/ref.fa" > "$W/edge.seq" \
    || fail "generate edge reads"
q100=$(printf '%100s' '' | tr ' ' 'I')
k=0
: > "$W/edge.fq"
while read -r r; do
    printf '@edge_f_%02d\n%s\n+\n%s\n@edge_r_%02d\n%s\n+\n%s\n' \
        "$k" "$r" "$q100" "$k" "$(dna_revcomp "$r")" "$q100" >> "$W/edge.fq"
    k=$((k + 1))
done < "$W/edge.seq"
[[ $k -eq 24 ]] || fail "generate edge reads ($k of 24)"
for f in pe_1.fq pe_2.fq se.fq long.fq edge.fq; do
    [[ -s "$W/$f" ]] || fail "empty fixture $f"
done

# run <tag> <layout: default|flag|env|env0|envbad> <threads> <input: pe|se|long|edge> [extra mem flags...]
run() {
    local tag=$1 layout=$2 t=$3 in=$4
    shift 4
    local -a reads
    case $in in
        pe) reads=("$W/pe_1.fq" "$W/pe_2.fq") ;;
        se) reads=("$W/se.fq") ;;
        long) reads=("$W/long.fq") ;;
        edge) reads=("$W/edge.fq") ;;
    esac
    local -a flags=()
    local -a envs=(env -u BWA3_COMPACT_OCC)
    case $layout in
        flag) flags=(--compact-occ) ;;
        env) envs=(env BWA3_COMPACT_OCC=1) ;;
        env0) envs=(env BWA3_COMPACT_OCC=0) ;;
        envbad) envs=(env BWA3_COMPACT_OCC=yes) ;;
    esac
    "${envs[@]}" "$BWA_MEM3" mem -t "$t" ${flags[@]+"${flags[@]}"} "$@" "$W/ref.fa" "${reads[@]}" \
        > "$W/$tag.sam" 2> "$W/$tag.err" || {
        cat "$W/$tag.err" >&2
        fail "$tag: mem exited non-zero"
    }
    if [[ $layout == flag ]]; then
        grep -qF "occurrence checkpoints: compact layout (--compact-occ)" "$W/$tag.err" \
            || fail "$tag: --compact-occ run did not report the compact layout from the flag"
    elif [[ $layout == env ]]; then
        grep -qF "occurrence checkpoints: compact layout (BWA3_COMPACT_OCC)" "$W/$tag.err" \
            || fail "$tag: BWA3_COMPACT_OCC=1 run did not report the compact layout from the variable"
    else
        if grep -q "occurrence checkpoints: compact layout" "$W/$tag.err"; then
            fail "$tag: default run reported the compact layout"
        fi
    fi
    grep -v '^@PG' "$W/$tag.sam" > "$W/$tag.body"
    [[ $(grep -vc '^@' "$W/$tag.body") -gt 0 ]] || fail "$tag: no alignment records"
}

same() { # same <tag-a> <tag-b>
    cmp -s "$W/$1.body" "$W/$2.body" && return 0
    echo "FAIL: $2 differs from $1" >&2
    { diff "$W/$1.body" "$W/$2.body" | head -20 >&2; } || true
    exit 1
}

n=0
while read -r name t in extra; do
    # shellcheck disable=SC2086 # extra is a deliberately word-split flag list
    run "$name.def" default "$t" "$in" $extra
    for layout in flag env; do
        # shellcheck disable=SC2086
        run "$name.$layout" "$layout" "$t" "$in" $extra
        same "$name.def" "$name.$layout"
    done
    ok "$name: --compact-occ and BWA3_COMPACT_OCC=1 byte-identical to the default layout ($(grep -vc '^@' "$W/$name.def.body") records)"
    n=$((n + 1))
done << 'EOF'
pe-t1 1 pe
pe-t4 4 pe
se-t4 4 se
long-t1 1 long
pe-compat-t4 4 pe --compat=bwa-mem2
pe-huge-t4 4 pe --huge-pages
edge-t1 1 edge
edge-compat-t1 1 edge --compat=bwa-mem2
EOF
[[ $n -eq 8 ]] || fail "ran $n of 8 configurations"

# The edge reads must really reach the sentinel walk: --compat=bwa-mem2 drops
# its offset, so some edge read is placed at a different position than in the
# default run.
pos() { grep -v '^@' "$1" | cut -f1,3,4 | sort; }
moved=$(($(comm -23 <(pos "$W/edge-t1.def.body") <(pos "$W/edge-compat-t1.def.body") | wc -l)))
[[ $moved -gt 0 ]] \
    || fail "edge reads never reached the sentinel walk (--compat=bwa-mem2 places every edge read as the default does)"
ok "edge reads reach the sentinel walk ($moved placed differently under --compat=bwa-mem2)"

run pe-env0 env0 4 pe
same pe-t4.def pe-env0
ok "BWA3_COMPACT_OCC=0 keeps the default layout"

run pe-envbad envbad 4 pe
same pe-t4.def pe-envbad
grep -qF '[W::compact_occ_env] ignoring BWA3_COMPACT_OCC="yes" (not 0 or 1)' "$W/pe-envbad.err" \
    || fail "BWA3_COMPACT_OCC=yes was not reported"
ok "BWA3_COMPACT_OCC=yes is reported and keeps the default layout"

# --meth: the f/r-doubled seed index goes through the same compact load.
"$BWA_MEM3" index --meth --force "$W/ref.fa" > "$W/index-meth.log" 2>&1 || {
    cat "$W/index-meth.log" >&2
    fail "index --meth"
}
for mode in se pe; do
    run "meth-$mode.def" default 4 "$mode" --meth
    run "meth-$mode.flag" flag 4 "$mode" --meth
    same "meth-$mode.def" "meth-$mode.flag"
done
ok "--meth --compact-occ byte-identical to --meth ($(grep -vc '^@' "$W/meth-pe.def.body") pe records)"

# `index` never reads BWA3_COMPACT_OCC: the files match a build without it.
mkdir -p "$W/idx-env"
cp "$COMPACT_OCC_PHIX_FA" "$W/idx-env/a.fa"
cp "$COMPACT_OCC_PHIX_FA" "$W/idx-env/b.fa"
env -u BWA3_COMPACT_OCC "$BWA_MEM3" index "$W/idx-env/a.fa" > /dev/null 2>&1 || fail "index (no variable)"
env BWA3_COMPACT_OCC=1 "$BWA_MEM3" index "$W/idx-env/b.fa" > /dev/null 2>&1 || fail "index (BWA3_COMPACT_OCC=1)"
# (.0123 is written only with --emit-unpacked-ref)
for ext in bwt.2bit.64 ann amb pac; do
    [[ -s "$W/idx-env/a.fa.$ext" ]] || fail "index did not write .$ext"
    cmp -s "$W/idx-env/a.fa.$ext" "$W/idx-env/b.fa.$ext" || fail "index with BWA3_COMPACT_OCC=1 differs in .$ext"
done
ok "index ignores BWA3_COMPACT_OCC"

# Only @PG differs between a flag run and a default run (it records the flag).
grep -q -- '--compact-occ' "$W/pe-t4.flag.sam" || fail "@PG CL does not record --compact-occ"
ok "@PG CL records --compact-occ"
