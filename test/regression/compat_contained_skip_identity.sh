#!/usr/bin/env bash
# test/regression/compat_contained_skip_identity.sh
#
# Regression: `--compat` keeps the contained-seed extension skip, and under
# both targets the skip is byte-identical to the reference extension path
# (`--keep-contained-ext`).
#
# The compat policy is proof-based: a target keeps every optimization that is
# byte-identical by construction or proven inside a code-enforced envelope, and
# drops only output shaping. The two-wave skip is proven (see
# mem_skip_contained_ext_sound in src/bwamem.cpp), so `--compat` no longer
# forces the reference path. This pins three things:
#
#   1. Under `--compat=bwa-mem2` and `--compat=bwa-mem`, the default (skip on)
#      SAM is byte-identical to `--keep-contained-ext`, at -t 1 and -t 3, on
#      the default scoring and on scorings that move every term of the
#      containment test (-A 2 -B 8 -O 8 -E 2, -L 0, -d 3, -w 20) and on -W 200,
#      which drops every chain of a 150 bp read and so drives the chain-filter
#      resurrection that only the bwa-mem2 target reproduces (#310). The reads
#      also include slices of the first bases of contig 0, whose suffix-array
#      walk reaches the sentinel row (the offset the bwa-mem2 target drops,
#      #469), so both known port divergences run with the skip in force.
#   2. The fixture really exercises both halves of the skip: with
#      BWA3_CHAIN_STATS=1 every compat skip run must report
#      contained_deferred > 0 (except --compat=bwa-mem -W 200 on the PE input,
#      where no chain survives), and over the run set both contained_purged
#      (seeds whose SW was skipped) and contained_extended (seeds the guard
#      sent to the second batch) must be > 0; every --keep-contained-ext run
#      must defer nothing. A fixture that never defers, or never extends, a
#      seed would make (1) vacuous for that half.
#   3. The proof envelope is enforced by the extension driver and reported: a
#      negative match score (-A -1, the one scoring where cal_max_gap is not
#      non-decreasing) prints the `contained-seed extension skip disabled`
#      notice and defers nothing although the flag stays set; -A 0 is inside
#      the envelope, defers seeds and prints no notice.
#
# Fixture (deterministic, no PRNG): a two-contig reference -- the committed
# phix.fa as contig 0, and a contig assembled from phix slices in which short
# units (24, 40, 60 bp) are planted twice at tandem distances 0, 30 and 300 bp,
# once exactly and once with a single substitution. A read spanning a planted
# unit has one long unique seed and, from re-seeding, a shorter seed for the
# unit on the same diagonal (contained -> deferred) whose other copies fall
# inside the chaining band on a shifted diagonal (the interference guard, at
# and around its .95 length boundary) or far enough away to chain separately.
# SE reads are 150 and 250 bp windows, plus 2-3 kb reads from each contig with
# ~1% substitutions and ~0.2% 1 bp indels (test/fixtures/make_long_reads.awk,
# deterministic across awk implementations). Long reads are rescored by
# mem_flt_chained_seeds, so container/contained visit order can invert and the
# guard sends the seed to the second batch. PE reads
# pair a window with the reverse complement of a downstream window.
#
# Inputs (env vars):
#   BWA_MEM3                  — path to the bwa-mem3 binary under test
#   COMPAT_CONTAINED_PHIX_FA  — phix reference FASTA
#   COMPAT_CONTAINED_WORK_DIR — private scratch directory for this script

set -euo pipefail
: "${BWA_MEM3:?BWA_MEM3 must be set}"
: "${COMPAT_CONTAINED_PHIX_FA:?COMPAT_CONTAINED_PHIX_FA must be set}"
: "${COMPAT_CONTAINED_WORK_DIR:?COMPAT_CONTAINED_WORK_DIR must be set}"

HERE="$(cd "$(dirname "$0")" && pwd)"
fail() {
    echo "FAIL: $*" >&2
    exit 1
}
ok() { echo "PASS: $*"; }

W="$COMPAT_CONTAINED_WORK_DIR"
mkdir -p "$W" || fail "create work directory"

# --- Reference: phix + a planted tandem-repeat contig (deterministic). -------
SEQ=$(awk '/^>/{next}{printf "%s",$0}' "$COMPAT_CONTAINED_PHIX_FA" | tr 'acgt' 'ACGT')
[[ ${#SEQ} -ge 5000 ]] || fail "phix reference shorter than expected (${#SEQ} bp)"

# substitute position 1 of a unit with a different base (deterministic swap)
mutate1() {
    local u=$1 b sub
    b=${u:1:1}
    case $b in
        A) sub=C ;;
        C) sub=G ;;
        G) sub=T ;;
        *) sub=A ;;
    esac
    printf '%s%s%s' "${u:0:1}" "$sub" "${u:2}"
}

REP=""
i=0
for U in 24 40 60; do
    for D in 0 30 300; do
        unit=${SEQ:$((1000 + i * 97)):U}
        spacer=${SEQ:3000:D}
        filler=${SEQ:$((4000 + i * 53)):160}
        REP+="${unit}${spacer}${unit}${filler}"
        REP+="${unit}${spacer}$(mutate1 "$unit")${filler}"
        i=$((i + 1))
    done
done
{
    cat "$COMPAT_CONTAINED_PHIX_FA"
    echo ">rep"
    printf '%s\n' "$REP" | fold -w 60
} > "$W/ref.fa"
"$BWA_MEM3" index "$W/ref.fa" > /dev/null 2>&1 || fail "index"

# --- Reads (deterministic slices). --------------------------------------------
revcomp() { printf '%s' "$1" | rev | tr ACGTacgt TGCAtgca; }
QUAL150=$(printf 'I%.0s' $(seq 150))
QUAL250=$(printf 'I%.0s' $(seq 250))

se="$W/se.fq"
: > "$se"
n=0
# 150 bp windows over the repeat contig (stride 23) and phix (stride 37).
for ((off = 0; off + 150 <= ${#REP}; off += 23)); do
    n=$((n + 1))
    printf '@rep150_%04d\n%s\n+\n%s\n' "$n" "${REP:$off:150}" "$QUAL150" >> "$se"
done
for ((off = 0; off + 150 <= ${#SEQ}; off += 37)); do
    n=$((n + 1))
    printf '@phix150_%04d\n%s\n+\n%s\n' "$n" "${SEQ:$off:150}" "$QUAL150" >> "$se"
done
# 250 bp windows over the repeat contig (stride 31); rescoring range.
for ((off = 0; off + 250 <= ${#REP}; off += 31)); do
    n=$((n + 1))
    printf '@rep250_%04d\n%s\n+\n%s\n' "$n" "${REP:$off:250}" "$QUAL250" >> "$se"
done
# 150 bp windows with a 3 bp deletion (gapped extension), reverse strand.
for ((off = 5; off + 153 <= ${#REP}; off += 41)); do
    n=$((n + 1))
    r=${REP:$off:153}
    r="${r:0:70}${r:73}"
    printf '@repdel_%04d\n%s\n+\n%s\n' "$n" "$(revcomp "$r")" "$QUAL150" >> "$se"
done
# Openers of contig 0 (SA sentinel walk) at three offsets and lengths.
for off in 0 1 2; do
    for L in 50 100 150; do
        n=$((n + 1))
        printf '@open_%d_%d\n%s\n+\n%s\n' "$off" "$L" "${SEQ:$off:$L}" "$(printf 'I%.0s' $(seq "$L"))" >> "$se"
    done
done

# 2-3 kb reads with scattered errors, 10 per contig (second-batch extensions).
contig0=$(awk '/^>/{print substr($1, 2); exit}' "$COMPAT_CONTAINED_PHIX_FA")
for c in "$contig0" rep; do
    awk -v NAME="$c" -v SEED=7 -v NREADS=10 -v LMIN=2000 -v LMAX=3000 \
        -f "$HERE/../fixtures/make_long_reads.awk" "$W/ref.fa" >> "$se" || fail "long-read generator"
done

r1="$W/pe_1.fq"
r2="$W/pe_2.fq"
: > "$r1"
: > "$r2"
m=0
for ((off = 0; off + 300 <= ${#REP}; off += 29)); do
    m=$((m + 1))
    printf '@pe_%04d/1\n%s\n+\n%s\n' "$m" "${REP:$off:150}" "$QUAL150" >> "$r1"
    printf '@pe_%04d/2\n%s\n+\n%s\n' "$m" "$(revcomp "${REP:$((off + 150)):150}")" "$QUAL150" >> "$r2"
done
for ((off = 0; off + 300 <= ${#SEQ}; off += 61)); do
    m=$((m + 1))
    printf '@pephix_%04d/1\n%s\n+\n%s\n' "$m" "${SEQ:$off:150}" "$QUAL150" >> "$r1"
    printf '@pephix_%04d/2\n%s\n+\n%s\n' "$m" "$(revcomp "${SEQ:$((off + 150)):150}")" "$QUAL150" >> "$r2"
done
[[ $n -ge 200 && $m -ge 100 ]] || fail "fixture too small (se=$n pe=$m)"

# --- Runs. --------------------------------------------------------------------
# run <tag> <input: se|pe> <threads> [flags...]  -> $W/<tag>.sam (no @PG), $W/<tag>.err
run() {
    local tag=$1 in=$2 t=$3
    shift 3
    local -a reads
    if [[ $in == se ]]; then reads=("$se"); else reads=("$r1" "$r2"); fi
    BWA3_CHAIN_STATS=1 "$BWA_MEM3" mem -t "$t" "$@" "$W/ref.fa" "${reads[@]}" \
        2> "$W/$tag.err" | grep -v '^@PG' > "$W/$tag.sam" \
        || fail "$tag: mem run failed"
    [[ -s "$W/$tag.sam" ]] || fail "$tag: empty output"
}
# counter <err-file> <name>  -> the value of contained_<name> on the [chain-stats] line
counter() {
    local v
    v=$(sed -n 's/.*contained_'"$2"'=\([0-9]*\).*/\1/p' "$1" | tail -n 1)
    [[ -n $v ]] || fail "$1: no contained_$2 counter on the [chain-stats] line"
    echo "$v"
}

params=("" "-A 2 -B 8 -O 8 -E 2" "-L 0" "-d 3" "-w 20" "-W 200")
total_purged=0
total_extended=0
for target in bwa-mem2 bwa-mem; do
    for p in "${params[@]}"; do
        # shellcheck disable=SC2206  # word-split the parameter set on purpose
        flags=($p)
        ptag=${p//[^A-Za-z0-9]/}
        for in in se pe; do
            for t in 1 3; do
                tag="${target}_${ptag:-def}_${in}_t${t}"
                run "skip_$tag" "$in" "$t" --compat="$target" ${flags[@]+"${flags[@]}"}
                run "keep_$tag" "$in" "$t" --compat="$target" ${flags[@]+"${flags[@]}"} --keep-contained-ext
                cmp -s "$W/skip_$tag.sam" "$W/keep_$tag.sam" \
                    || fail "--compat=$target ${p:-defaults} ($in, -t $t): skip != --keep-contained-ext"
                grep -q 'contained-seed extension skip disabled' "$W/skip_$tag.err" \
                    && fail "$tag: the skip must not be disabled inside its envelope"
                d=$(counter "$W/skip_$tag.err" deferred)
                pu=$(counter "$W/skip_$tag.err" purged)
                ex=$(counter "$W/skip_$tag.err" extended)
                kd=$(counter "$W/keep_$tag.err" deferred)
                [[ $kd -eq 0 ]] || fail "$tag: --keep-contained-ext deferred $kd seeds (skip not off)"
                # -W 200 drops every chain of the 150 bp PE reads; only the bwa-mem2
                # target resurrects one, so the bwa-mem target has nothing to defer
                # there (the SE input's longer reads keep chains and still defer).
                if [[ $p != "-W 200" || $target == bwa-mem2 || $in == se ]]; then
                    [[ $d -gt 0 ]] || fail "$tag: no seed deferred -- fixture does not exercise the skip"
                fi
                [[ $((pu + ex)) -eq "$d" ]] || fail "$tag: deferred=$d != purged=$pu + extended=$ex"
                total_purged=$((total_purged + pu))
                total_extended=$((total_extended + ex))
                ok "--compat=$target ${p:-defaults} ($in, -t $t): byte-identical; deferred=$d purged=$pu extended=$ex"
            done
        done
    done
done
[[ $total_purged -gt 0 ]] || fail "no deferred seed was ever purged -- the skip saved nothing"
[[ $total_extended -gt 0 ]] || fail "no deferred seed was ever extended -- the second batch was never run"
ok "skip exercised under --compat: purged=$total_purged extended=$total_extended"

# --- Proof envelope: outside it the driver disables the skip and says so. -----
# A negative -A scales every penalty the command line leaves unset (-B, -E, -d,
# -L, -U ...) by -A, which makes the z-drop negative and aborts the extension
# kernel on both paths, so the penalties are pinned explicitly here. -A -1 is
# outside the envelope; -A 0 is inside it. main_mem only reports: the flag
# stays set, so "deferred nothing" is the driver's own guard at work.
envelope=(-B 4 -O 6 -E 1 -d 100 -L 5 -U 17)
run "env_out" se 1 --compat=bwa-mem2 -A -1 "${envelope[@]}"
grep -q 'contained-seed extension skip disabled' "$W/env_out.err" \
    || fail "-A -1: expected the skip-disabled notice"
od=$(counter "$W/env_out.err" deferred)
[[ $od -eq 0 ]] || fail "-A -1: $od seeds were deferred outside the proof envelope"
run "env_out_keep" se 1 --compat=bwa-mem2 -A -1 "${envelope[@]}" --keep-contained-ext
cmp -s "$W/env_out.sam" "$W/env_out_keep.sam" \
    || fail "-A -1: disabled skip != --keep-contained-ext"
run "env_in" se 1 --compat=bwa-mem2 -A 0 "${envelope[@]}"
grep -q 'contained-seed extension skip disabled' "$W/env_in.err" \
    && fail "-A 0 (inside the envelope) must not disable the skip"
id=$(counter "$W/env_in.err" deferred)
[[ $id -gt 0 ]] || fail "-A 0: inside the envelope, yet no seed was deferred"
run "env_in_keep" se 1 --compat=bwa-mem2 -A 0 "${envelope[@]}" --keep-contained-ext
cmp -s "$W/env_in.sam" "$W/env_in_keep.sam" || fail "-A 0: skip != --keep-contained-ext"
ok "proof envelope enforced by the driver and reported (-A -1 disabled, -A 0 kept)"

# --- Record-count sanity. -----------------------------------------------------
rc=0
cnt=$(grep -vc '^@' "$W/skip_bwa-mem2_def_se_t1.sam") || rc=$?
[[ $rc -gt 1 ]] && fail "grep read error (status $rc)"
[[ $cnt -ge $n ]] || fail "record count $cnt < $n reads"
ok "record count sanity ($cnt records for $n SE reads)"
