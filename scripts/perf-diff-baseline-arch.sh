#!/usr/bin/env bash
# scripts/perf-diff-baseline-arch.sh
#
# A/B comparison harness for BASELINE_ARCH choices on a single host.
# Builds bwa-mem3 with two (or more) different non-kernel ARCH_FLAGS
# values, runs each on the same input, and emits a comparison table.
#
# Both builds use the production `single` target, so the per-tier kernel
# objects (bandedSWA, kswv, ksw, sam_encode at every tier) are identical
# between runs and the runtime dispatcher picks the host's tier
# regardless. The only thing that differs between runs is the
# auto-vectorization width and ISA used by NON-kernel TUs — which is
# what we attribute when investigating BASELINE_ARCH choices.
#
# Modes:
#   time    — wrap each run with `tricorder` (fg-labs/tricord). Same tool
#             the bwa-mem3-bench Snakemake workflow uses, so the wall /
#             cpu_time / max_rss / io numbers are directly comparable to
#             benchmark.db rows. Fast (one invocation per rep). Default.
#   record  — `perf record -F 999 -g --call-graph dwarf` × 1 rep per
#             variant. Used for finding individual hot functions where one
#             variant beats another. For Phase 3 pragma-target candidate
#             selection.
#
# Variants are passed as a comma-separated list of <label>:<flags>
# tokens. Examples:
#   --variants 'avx2:,avx512bw:'
#       avx2 = default ARCH_FLAGS for arch=avx2
#       avx512bw = default ARCH_FLAGS for arch=avx512bw
#   --variants 'avx2:,avx512bw:,avx512bw-pvw256:-mprefer-vector-width=256'
#       three-way comparison; the third variant builds avx512bw with
#       the autovec-width cap appended via EXTRA_CXXFLAGS.
#
# Usage:
#   scripts/perf-diff-baseline-arch.sh \
#       --ref <ref.fa> --r1 <r1.fq.gz> --r2 <r2.fq.gz> \
#       --out <out_dir> [--mode time|record] [--reps N] [--threads T] \
#       [--variants <label>:<extra_flags>,...] [--no-shm]
#
# Output files in <out_dir> (mode=time):
#   bwa-mem3.<label>           — built binary per variant
#   tric-<label>-rep<N>.tsv    — tricorder TSV per rep
#   summary.csv                — median per variant per metric
#   delta.csv                  — variant deltas vs the first variant
# Output files (mode=record):
#   record-<label>.data        — perf.data per variant
#   record-<label>.txt         — perf report --stdio
#   delta.tsv                  — function self-time delta (variant 2 vs 1)
#
# Requires: g++ ≥ 11; tricorder (cargo install tricord) on PATH for
# mode=time; perf for mode=record. Run on the target host (c7a / c7i /
# m7i / etc.). Lower kernel.perf_event_paranoid (-1 on a throwaway box).
# /dev/shm must hold the staged FMI index when --no-shm is not set
# (~17 GB for hg38; remount with `mount -o remount,size=28g /dev/shm`
# if the default tmpfs is too small).

set -euo pipefail

REF=""
R1=""
R2=""
OUT=""
MODE="time"
REPS=3
THREADS="$(nproc 2> /dev/null || echo 4)"
VARIANTS_RAW="avx2:,avx512bw:"
USE_SHM=1

while [ $# -gt 0 ]; do
    case "$1" in
        --ref)
            REF="$2"
            shift 2
            ;;
        --r1)
            R1="$2"
            shift 2
            ;;
        --r2)
            R2="$2"
            shift 2
            ;;
        --out)
            OUT="$2"
            shift 2
            ;;
        --mode)
            MODE="$2"
            shift 2
            ;;
        --reps)
            REPS="$2"
            shift 2
            ;;
        --threads)
            THREADS="$2"
            shift 2
            ;;
        --variants)
            VARIANTS_RAW="$2"
            shift 2
            ;;
        --no-shm)
            USE_SHM=0
            shift
            ;;
        -h | --help)
            sed -n '/^# /,/^$/p' "$0" | sed 's/^# //;s/^#$//' | head -55
            exit 0
            ;;
        *)
            echo "unknown arg: $1" >&2
            exit 2
            ;;
    esac
done

: "${REF:?--ref is required}"
: "${R1:?--r1 is required}"
: "${R2:?--r2 is required}"
: "${OUT:?--out is required}"

case "$MODE" in
    time | record) ;;
    *)
        echo "--mode must be time or record (got: $MODE)" >&2
        exit 2
        ;;
esac

if [ "$MODE" = "time" ] && ! command -v tricorder > /dev/null 2>&1; then
    echo "tricorder not found on PATH; install with: cargo install tricord" >&2
    exit 2
fi

REF="$(readlink -f "$REF")"
R1="$(readlink -f "$R1")"
R2="$(readlink -f "$R2")"
mkdir -p "$OUT"
OUT="$(readlink -f "$OUT")"

REPO_ROOT="$(git rev-parse --show-toplevel)"

# Parse VARIANTS_RAW into parallel arrays LABELS / FLAGS.
LABELS=()
FLAGS=()
IFS=',' read -ra _VARIANTS <<< "$VARIANTS_RAW"
for v in "${_VARIANTS[@]}"; do
    label="${v%%:*}"
    flags="${v#*:}"
    LABELS+=("$label")
    FLAGS+=("$flags")
done
echo "Variants: ${#LABELS[@]}"
for i in "${!LABELS[@]}"; do
    echo "  ${LABELS[$i]}: arch_flags+='${FLAGS[$i]}'"
done

# Cleanup helper: remove only bwa-mem3's own build artifacts.
clean_bwamem3_only() {
    cd "$REPO_ROOT"
    rm -f src/*.o src/*/*.o libbwa.a bwa-mem3 src/version.h
    rm -f src/*.gcno src/*.gcda
}

# Build one variant. The label maps to BASELINE_ARCH (first segment
# before any '-'); the flags string is appended via EXTRA_CXXFLAGS.
build_one() {
    local label="$1" flags="$2"
    cd "$REPO_ROOT"
    clean_bwamem3_only
    echo "=== Building variant=$label flags='$flags' ==="
    local baseline="${label%%-*}"
    EXTRA_CXXFLAGS="$flags" make BASELINE_ARCH="$baseline" -j"$THREADS" >&2
    cp bwa-mem3 "$OUT/bwa-mem3.$label"
    echo "=== Built $OUT/bwa-mem3.$label ==="
}

run_tricord_rep() {
    local label="$1" rep="$2"
    local bin="$OUT/bwa-mem3.$label"
    local tsv="$OUT/tric-$label-rep$rep.tsv"
    echo "=== tricord $label rep $rep ==="
    tricorder --out "$tsv" -- \
        "$bin" mem -t "$THREADS" "$REF" "$R1" "$R2" \
        > /dev/null 2> "$OUT/stderr-$label-rep$rep.log"
}

run_record() {
    local label="$1"
    local bin="$OUT/bwa-mem3.$label"
    local data="$OUT/record-$label.data"
    local txt="$OUT/record-$label.txt"
    echo "=== perf record $label ==="
    perf record -F 999 -g --call-graph dwarf -o "$data" -- \
        "$bin" mem -t "$THREADS" "$REF" "$R1" "$R2" > /dev/null
    perf report -i "$data" --no-children --stdio --percent-limit 0.05 \
        --demangle > "$txt"
}

# Build all variants.
for i in "${!LABELS[@]}"; do
    build_one "${LABELS[$i]}" "${FLAGS[$i]}"
done

# Shared-memory FMI index prewarm/teardown. Mirrors the bench's pattern:
# pull the ~17 GB hg38 index load out of the timed region so wall times
# reflect alignment work only. Index data is binary-independent (same
# source SHA across variants), so one prewarm covers every variant.
shm_cleanup() {
    if [ "$USE_SHM" = "1" ] && [ -n "${SHM_PREWARM_BIN:-}" ]; then
        echo "=== shm teardown ==="
        "$SHM_PREWARM_BIN" shm -d 2> /dev/null || true
    fi
}
trap shm_cleanup EXIT

if [ "$USE_SHM" = "1" ]; then
    SHM_PREWARM_BIN="$OUT/bwa-mem3.${LABELS[0]}"
    echo "=== shm prewarm via $SHM_PREWARM_BIN ==="
    "$SHM_PREWARM_BIN" shm -d 2> /dev/null || true
    "$SHM_PREWARM_BIN" shm "$REF"
fi

# Run.
if [ "$MODE" = "time" ]; then
    for label in "${LABELS[@]}"; do
        for rep in $(seq 1 "$REPS"); do
            run_tricord_rep "$label" "$rep"
        done
    done

    # tricorder TSV columns (header + 1 data row):
    #   s  h:m:s  max_rss  max_vms  max_uss  max_pss  io_in  io_out  mean_load  cpu_time
    # We summarize the numeric columns we care about: s, max_rss, io_in,
    # io_out, mean_load, cpu_time. The metric label "s" is renamed to
    # "wall_seconds" for clarity.
    {
        echo "variant,metric,median"
        for label in "${LABELS[@]}"; do
            for col_idx in 1 3 7 8 9 10; do
                # Map col_idx → metric name
                case $col_idx in
                    1) metric=wall_seconds ;;
                    3) metric=max_rss_mib ;;
                    7) metric=io_in_mib ;;
                    8) metric=io_out_mib ;;
                    9) metric=mean_load_pct_one_core ;;
                    10) metric=cpu_seconds ;;
                esac
                vals=()
                for rep in $(seq 1 "$REPS"); do
                    local_tsv="$OUT/tric-$label-rep$rep.tsv"
                    [ -f "$local_tsv" ] || continue
                    v=$(awk -v c=$col_idx 'NR==2{print $c}' "$local_tsv")
                    [ -n "$v" ] && [ "$v" != "NA" ] && [ "$v" != "-" ] && vals+=("$v")
                done
                if [ ${#vals[@]} -gt 0 ]; then
                    median=$(printf '%s\n' "${vals[@]}" | sort -n \
                        | awk -v n=${#vals[@]} 'NR==int((n+1)/2)')
                    echo "$label,$metric,$median"
                fi
            done
        done
    } > "$OUT/summary.csv"

    BASE_LABEL="${LABELS[0]}"
    {
        echo "variant,metric,base_median,variant_median,abs_delta,rel_delta_pct"
        awk -F',' -v base="$BASE_LABEL" '
            NR==1 {next}
            { val[$1, $2] = $3; seen_var[$1] = 1; seen_m[$2] = 1 }
            END {
                for (var in seen_var) {
                    if (var == base) continue
                    for (m in seen_m) {
                        b = val[base, m] + 0
                        v = val[var, m] + 0
                        if (b == 0) continue
                        absd = v - b
                        rel = (absd / b) * 100
                        printf "%s,%s,%g,%g,%g,%.2f\n", var, m, b, v, absd, rel
                    }
                }
            }
        ' "$OUT/summary.csv"
    } | sort -t',' -k1,1 -k2,2 > "$OUT/delta.csv"

    echo
    echo "=== Summary (tricord median of $REPS reps) ==="
    column -t -s',' < "$OUT/summary.csv"
    echo
    echo "=== Deltas vs $BASE_LABEL ==="
    column -t -s',' < "$OUT/delta.csv"
    echo
    echo "Files: $OUT/summary.csv $OUT/delta.csv"

elif [ "$MODE" = "record" ]; then
    for label in "${LABELS[@]}"; do
        run_record "$label"
    done

    # Build per-function self-time delta (variant 2 vs variant 1).
    extract_self_time() {
        awk '
            /^[[:space:]]*[0-9]/ {
                pct = $1; gsub(/%/, "", pct);
                for (i = 1; i <= NF; i++) {
                    if ($i == "[.]" || $i == "[k]") {
                        sym = "";
                        for (j = i+1; j <= NF; j++) {
                            sym = (sym == "" ? $j : sym " " $j);
                        }
                        print sym "\t" pct;
                        next;
                    }
                }
            }
        ' "$1"
    }

    if [ ${#LABELS[@]} -ge 2 ]; then
        L0="${LABELS[0]}"
        L1="${LABELS[1]}"
        extract_self_time "$OUT/record-$L0.txt" | sort -k1,1 > "$OUT/.self-$L0.tsv"
        extract_self_time "$OUT/record-$L1.txt" | sort -k1,1 > "$OUT/.self-$L1.tsv"
        awk -F'\t' -v l0="$L0" -v l1="$L1" '
            NR==FNR { a[$1] = $2; next }
            { b[$1] = $2 }
            END {
                for (s in a) all[s] = 1
                for (s in b) all[s] = 1
                printf "function\tself_%s_pct\tself_%s_pct\tabs_delta_pct\trel_delta_pct\n", l0, l1
                for (s in all) {
                    av = (s in a) ? a[s] + 0 : 0
                    bv = (s in b) ? b[s] + 0 : 0
                    absd = bv - av
                    rel = (av > 0) ? (absd / av) * 100 : 0
                    printf "%s\t%.3f\t%.3f\t%.3f\t%.2f\n", s, av, bv, absd, rel
                }
            }
        ' "$OUT/.self-$L0.tsv" "$OUT/.self-$L1.tsv" \
            | {
                read -r header
                printf '%s\n' "$header"
                sort -t$'\t' -k4,4 -gr
            } \
                > "$OUT/delta.tsv"
        rm -f "$OUT/.self-$L0.tsv" "$OUT/.self-$L1.tsv"

        echo
        echo "=== Top 30 functions by abs_delta_pct ($L1 vs $L0) ==="
        head -31 "$OUT/delta.tsv" | column -t -s$'\t'
        echo
        echo "Files: $OUT/delta.tsv"
    fi
fi
