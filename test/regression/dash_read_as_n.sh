#!/usr/bin/env bash
# Regression: a read base '-' is aligned and printed exactly as an N.
#
# nst_nt4_table maps '-' to 4, the N code. Upstream bwa/bwa-mem2 map it to a
# code 5 that no scoring-matrix row defines, so each alignment kernel scored it
# by however its N test happened to work: w_ambig in the 8-bit extension
# kernels, w_mismatch in the 16-bit AVX2 and NEON ones, w_ambig in the 16-bit
# AVX-512BW one, a read past the matrix row in the scalar paths. The output for
# a read with a '-' then depended on the SIMD tier and the kernel width.
#
# The oracle is exact and needs no golden file: a fixture with '-' bases must
# give byte-identical SAM (records, @PG excluded) to the same fixture with
# every '-' replaced by N, and the '-' fixture must give the same records at
# every tier this host can run (the unforced host tier plus every x86 tier at or
# below it, forced with BWAMEM3_FORCE_TIER). Checked paired-end and single-end,
# at -t 1 and -t 4, at the default scoring, with BWAMEM3_DISABLE_BSW8=1 (every
# extension on the 16-bit kernels, where the tiers disagreed), at a scaled
# scoring (-A 3 -B 12 -O 18 -E 3), and under both --compat targets (which also
# score '-' as N).
#
# Env:
#   BWA_MEM3        -- the bwa-mem3 binary
#   DASH_PHIX_FA    -- the committed test/fixtures/phix.fa
#   DASH_WORK_DIR   -- scratch directory (created)
set -euo pipefail
: "${BWA_MEM3:?BWA_MEM3 must be set}"
: "${DASH_PHIX_FA:?DASH_PHIX_FA must be set}"
: "${DASH_WORK_DIR:?DASH_WORK_DIR must be set}"
HERE="$(cd "$(dirname "$0")" && pwd)"
fail() {
    echo "FAIL: $*" >&2
    exit 1
}
ok() { echo "PASS: $*"; }
W="$DASH_WORK_DIR"
mkdir -p "$W" || fail "create work directory"
cp "$DASH_PHIX_FA" "$W/phix.fa" || fail "copy reference"
"$BWA_MEM3" index "$W/phix.fa" > /dev/null 2>&1 || fail "index"

# 600 distinct concordant 97 bp pairs (each with a 3 bp deletion, so they reach
# the banded extension kernels rather than the ungapped fast path).
"$HERE/../fixtures/make_dedup_reads.sh" "$W/phix.fa" "$W/base" 600 1 \
    || fail "generate fixture"

# Two bases per read, at positions that vary with the record, become '-' (or N
# in the twin). Seeds cannot span them, so they land in the extended flanks.
mark() { # <in.fq> <out.fq> <char>
    awk -v ch="$3" 'NR % 4 == 2 {
        r = (NR - 2) / 4
        p1 = 3 + (r % 23); p2 = 66 + (r * 7) % 29
        print substr($0, 1, p1 - 1) ch substr($0, p1 + 1, p2 - p1 - 1) ch substr($0, p2 + 1)
        next
    } { print }' "$1" > "$2"
}
for m in 1 2; do
    mark "$W/base_$m.fq" "$W/dash_$m.fq" - || fail "generate '-' R$m"
    mark "$W/base_$m.fq" "$W/n_$m.fq" N || fail "generate N R$m"
done
grep -q -- '-' "$W/dash_1.fq" || fail "the '-' fixture has no '-'"
cmp -s "$W/dash_1.fq" "$W/n_1.fq" && fail "the '-' and N fixtures are identical"

# Tiers: the host's own (unforced) and every x86 tier at or below it that
# simd_sweep_tiers covers (none below avx2, which has no batched kswv kernel).
# shellcheck source=test/framework/simd_tiers.sh
source "$HERE/../framework/simd_tiers.sh"
simd_sweep_tiers "$BWA_MEM3" || fail "could not detect the host SIMD tier from BWAMEM3_DEBUG_SIMD output"
tiers=(host)
for t in ${SIMD_TIERS[@]+"${SIMD_TIERS[@]}"}; do
    [ "$t" = neon ] || tiers+=("$t")
done
echo "Host tier: $SIMD_HOST_TIER; tiers checked: ${tiers[*]}"

# run <out> <tier> <threads> <mode: pe|se> <prefix> <env/flags...>
run() {
    local out=$1 tier=$2 t=$3 mode=$4 p=$5
    shift 5
    local envs=() flags=()
    for a in "$@"; do
        case "$a" in
            *=*) [[ "$a" == --* ]] && flags+=("$a") || envs+=("$a") ;;
            *) flags+=("$a") ;;
        esac
    done
    [ "$tier" = host ] || envs+=("BWAMEM3_FORCE_TIER=$tier")
    local reads=("$W/${p}_1.fq")
    [ "$mode" = pe ] && reads+=("$W/${p}_2.fq")
    env ${envs[@]+"${envs[@]}"} "$BWA_MEM3" mem -t "$t" ${flags[@]+"${flags[@]}"} \
        "$W/phix.fa" "${reads[@]}" 2> "$out.err" | grep -v '^@PG' > "$out.sam" \
        || fail "$out: $(tail -1 "$out.err")"
    if grep -q 'ignoring BWAMEM3_FORCE_TIER' "$out.err"; then
        fail "the dispatcher ignored BWAMEM3_FORCE_TIER=$tier ($out)"
    fi
}

legs=(default bsw16 scaled compat2 compat1)
leg_args() { # leg -> env/flags
    case "$1" in
        default) ;;
        bsw16) echo BWAMEM3_DISABLE_BSW8=1 ;;
        scaled) echo -A 3 -B 12 -O 18 -E 3 ;;
        compat2) echo --compat=bwa-mem2 ;;
        compat1) echo --compat=bwa-mem ;;
    esac
}

n_checks=0
for leg in "${legs[@]}"; do
    read -r -a args <<< "$(leg_args "$leg")"
    for mode in pe se; do
        for t in 1 4; do
            first=""
            for tier in "${tiers[@]}"; do
                tag="$W/$leg.$mode.t$t.$tier"
                run "$tag.dash" "$tier" "$t" "$mode" dash ${args[@]+"${args[@]}"}
                run "$tag.n" "$tier" "$t" "$mode" n ${args[@]+"${args[@]}"}
                cmp -s "$tag.dash.sam" "$tag.n.sam" \
                    || fail "$leg $mode -t $t tier $tier: '-' fixture != N fixture ($(diff "$tag.dash.sam" "$tag.n.sam" | grep -c '^<') records differ)"
                if [ -z "$first" ]; then
                    first="$tag.dash.sam"
                else
                    cmp -s "$first" "$tag.dash.sam" \
                        || fail "$leg $mode -t $t: tier $tier != tier ${tiers[0]} on the '-' fixture"
                fi
                n_checks=$((n_checks + 1))
            done
        done
    done
done

# Non-vacuity: the '-' reads were aligned, through the N in their SEQ.
mapped_n=$(awk -F'\t' '!/^@/ && int($2 / 4) % 2 == 0 && $10 ~ /N/' "$W/default.pe.t1.host.dash.sam" | wc -l | tr -d " ")
[ "$mapped_n" -gt 1000 ] || fail "only $mapped_n mapped records carry an N base; the fixture did not exercise '-' columns"
ok "a read '-' aligns and prints as N: $n_checks fixture pairs equal across legs (${legs[*]}), PE/SE, -t 1/4, tiers ${tiers[*]} ($mapped_n mapped records with an N)"
