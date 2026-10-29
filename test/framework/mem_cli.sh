#!/usr/bin/env bash
# test/framework/mem_cli.sh
#
# Source-only helpers for regression scripts that check `bwa-mem3 mem` output
# field by field (test/regression/mem_cli_basics.sh is the worked example).
# Nothing here runs on its own; a script sources it with
#
#     # shellcheck source=test/framework/mem_cli.sh
#     source "$(cd "$(dirname "${BASH_SOURCE[0]}")/../framework" && pwd)/mem_cli.sh"
#
# calls mem_cli_init first and ends with
#
#     if mem_cli_verdict; then echo "PASS: $MEM_CLI_VERDICT"; else echo "FAIL: $MEM_CLI_VERDICT"; exit 1; fi
#
# The checks print per-check `ok:`, `FAIL:` and `SKIP:` lines and count them;
# the script still prints its own PASS:/FAIL: line and sets its exit status
# (test/regression/README.md).
#
# mem_cli_init sets an EXIT trap (mem_cli_cleanup) that removes the work
# directory; a script that needs its own EXIT handling must call
# mem_cli_cleanup from it rather than replace the trap.
#
# Portability: bash 3.2 (macOS /bin/bash), BSD awk, mawk and gawk, BSD and GNU
# userland. So no mapfile, no associative arrays, no awk extensions (and(),
# gensub, strtonum), no `sed -i`.
#
# Environment read here:
#   BWA_MEM3        the binary under test (required)
#   MEM_CLI_TIERS   optional space-separated SIMD tier list overriding
#                   mem_cli_tiers' detection (e.g. "avx2" to run one tier); each
#                   must be one of the sweepable tiers (avx2 avx512bw neon).
#                   Empty counts as unset (so a CI step can pass an empty
#                   value through); whitespace alone fails as naming no tier
#   MEM_CLI_KEEP    optional; 1 keeps the work directory and prints its path
#   OMP_NUM_THREADS defaults to 4 when unset, so `index` cannot claim a large
#                   host's every core
#   TMPDIR          where the work directory is created (default /tmp)

# shellcheck source=test/framework/simd_tiers.sh
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/simd_tiers.sh"
# shellcheck source=test/framework/dna.sh
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/dna.sh"

# mem_cli_init <script-name>
#   Checks BWA_MEM3, resolves the SIMD tiers to run (mem_cli_tiers), makes the
#   work directory MEM_CLI_WORK (removed on exit unless MEM_CLI_KEEP=1), sets
#   MEM_CLI_REPO and zeroes the check counters. On a host below the avx2 floor,
#   where `mem` refuses to run, the script reports SKIP and exits 0.
mem_cli_init() {
    MEM_CLI_NAME="$1"
    if [[ -z ${BWA_MEM3:-} ]]; then
        echo "FAIL: $MEM_CLI_NAME: BWA_MEM3 must be set" >&2
        exit 1
    fi
    if [[ ! -x $BWA_MEM3 ]]; then
        echo "FAIL: $MEM_CLI_NAME: BWA_MEM3=$BWA_MEM3 is not an executable" >&2
        exit 1
    fi
    BWA_MEM3="$(cd "$(dirname "$BWA_MEM3")" && pwd)/$(basename "$BWA_MEM3")"
    MEM_CLI_REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
    MEM_CLI_TIER_LIST=()
    if [[ -n ${MEM_CLI_TIERS:-} ]]; then
        local forced t
        read -r -a forced <<< "$MEM_CLI_TIERS"
        for t in ${forced[@]+"${forced[@]}"}; do
            case " ${SIMD_SWEEP_TIERS[*]} " in
                *" $t "*)
                    case " ${MEM_CLI_TIER_LIST[*]-} " in
                        *" $t "*)
                            echo "FAIL: $MEM_CLI_NAME: MEM_CLI_TIERS names '$t' more than once" >&2
                            exit 1
                            ;;
                    esac
                    MEM_CLI_TIER_LIST+=("$t")
                    ;;
                *)
                    echo "FAIL: $MEM_CLI_NAME: MEM_CLI_TIERS names '$t', not one of: ${SIMD_SWEEP_TIERS[*]}" >&2
                    exit 1
                    ;;
            esac
        done
        if [[ ${#MEM_CLI_TIER_LIST[@]} -eq 0 ]]; then
            echo "FAIL: $MEM_CLI_NAME: MEM_CLI_TIERS='$MEM_CLI_TIERS' names no tier" >&2
            exit 1
        fi
    else
        simd_sweep_tiers "$BWA_MEM3" || exit $?
        if [[ ${#SIMD_TIERS[@]} -eq 0 ]]; then
            echo "SKIP: $MEM_CLI_NAME: host tier '$SIMD_HOST_TIER' is below the avx2 floor, where mem does not run; nothing was checked. This is NOT a pass."
            exit 0
        fi
        MEM_CLI_TIER_LIST=("${SIMD_TIERS[@]}")
    fi
    local tmp="${TMPDIR:-/tmp}"
    MEM_CLI_WORK="$(mktemp -d "${tmp%/}/$MEM_CLI_NAME.XXXXXX")"
    trap mem_cli_cleanup EXIT
    export OMP_NUM_THREADS="${OMP_NUM_THREADS:-4}"
    MEM_CLI_CHECKS=0
    MEM_CLI_FAILS=0
    MEM_CLI_SKIPS=0
}

mem_cli_cleanup() {
    if [[ ${MEM_CLI_KEEP:-0} == 1 ]]; then
        echo "note: $MEM_CLI_NAME: work directory kept at $MEM_CLI_WORK" >&2
    else
        rm -rf "$MEM_CLI_WORK"
    fi
}

# mem_cli_verdict
#   Sets MEM_CLI_VERDICT to a one-line summary and returns 1 if any check
#   failed or none ran at all, else 0. The summary counts skipped checks
#   separately, so a run that skipped part of its work does not read as a
#   full pass.
# MEM_CLI_VERDICT is read by the calling script, which shellcheck sees separately.
# shellcheck disable=SC2034
mem_cli_verdict() {
    local tiers
    tiers="$(mem_cli_tiers | tr '\n' ' ')"
    if ((MEM_CLI_FAILS > 0)); then
        MEM_CLI_VERDICT="$MEM_CLI_NAME: $MEM_CLI_FAILS of $MEM_CLI_CHECKS checks failed (tiers: $tiers)"
        return 1
    fi
    if ((MEM_CLI_CHECKS == 0)); then
        MEM_CLI_VERDICT="$MEM_CLI_NAME: no check ran (tiers: $tiers)"
        return 1
    fi
    MEM_CLI_VERDICT="$MEM_CLI_NAME: $MEM_CLI_CHECKS checks on tiers: $tiers($MEM_CLI_SKIPS skipped)"
}

# ---------------------------------------------------------------------------
# References and indexes
# ---------------------------------------------------------------------------

# mem_cli_index <ref.fa>
#   `bwa-mem3 index` with its log beside the reference; a failure ends the
#   script with the log's tail.
mem_cli_index() {
    if ! "$BWA_MEM3" index "$1" > "$1.index.log" 2>&1; then
        echo "FAIL: $MEM_CLI_NAME: bwa-mem3 index $1 failed:" >&2
        tail -n 20 "$1.index.log" >&2
        exit 1
    fi
}

# mem_cli_setup_synth
#   Copies the committed synthetic reference (test/fixtures/mem_synth.fa; the
#   header of test/tools/mem_synth_ref.c lists its contigs and planted
#   repeats) into the work directory and indexes it, at MEM_CLI_SYNTH. The
#   copy's basename carries the script name and pid: `mem` first tries to
#   attach an index staged in shared memory (`bwa-mem3 shm`), matched by
#   basename alone, so a fixed name could pick up a staged index instead.
mem_cli_setup_synth() {
    MEM_CLI_SYNTH="$MEM_CLI_WORK/synth/$MEM_CLI_NAME.$$.fa"
    mkdir -p "$MEM_CLI_WORK/synth"
    if ! cp "$MEM_CLI_REPO/test/fixtures/mem_synth.fa" "$MEM_CLI_SYNTH"; then
        echo "FAIL: $MEM_CLI_NAME: cannot copy test/fixtures/mem_synth.fa to $MEM_CLI_SYNTH" >&2
        exit 1
    fi
    mem_cli_index "$MEM_CLI_SYNTH"
}

# ---------------------------------------------------------------------------
# Sequences and reads
# ---------------------------------------------------------------------------

# mem_cli_seq <ref.fa> <contig> <pos> <len>
#   The reference bases [pos, pos+len) of a contig, 1-based, as stored. A
#   range outside the contig (or an unknown contig) prints a FAIL line and
#   returns 1. Assign the result to a variable on its own
#   (`s="$(mem_cli_seq ...)"`), which ends a `set -e` caller; inside a larger
#   expression or as a command's argument the status is lost.
mem_cli_seq() {
    if ! awk -v c="$2" -v p="$3" -v l="$4" '
        /^>/ { on = (substr($1, 2) == c); found = found || on; next }
        on { s = s $0 }
        END {
            if (!found || p < 1 || length(s) < p + l - 1) exit 1
            print substr(s, p, l)
        }' "$1"; then
        echo "FAIL: ${MEM_CLI_NAME:-mem_cli}: mem_cli_seq $2:$3+$4 is outside the reference $1" >&2
        return 1
    fi
}

# mem_cli_subst <seq> <pos> <base>
#   <seq> with the base at 1-based <pos> replaced.
mem_cli_subst() {
    printf '%s%s%s\n' "${1:0:$(($2 - 1))}" "$3" "${1:$2}"
}

# mem_cli_fq <name> <seq> [<qual>]
#   One FASTQ record; the quality string defaults to all 'I' (Q40).
mem_cli_fq() {
    local q="${3:-}"
    if [[ -z $q ]]; then
        q="$(printf '%s\n' "$2" | sed 's/./I/g')"
    fi
    printf '@%s\n%s\n+\n%s\n' "$1" "$2" "$q"
}

# ---------------------------------------------------------------------------
# Running mem
# ---------------------------------------------------------------------------

# mem_cli_tiers
#   The SIMD tiers `mem` can run on this host, one per line, lowest first, as
#   mem_cli_init resolved them: avx2 and, where the host has it, avx512bw on
#   x86; neon on arm64 (simd_sweep_tiers in test/framework/simd_tiers.sh), or
#   MEM_CLI_TIERS when set. Resolved once in mem_cli_init rather than here,
#   because callers read this through $(...), where a failed detection's exit
#   would end only the subshell and leave the script looping over no tiers.
mem_cli_tiers() {
    printf '%s\n' "${MEM_CLI_TIER_LIST[@]}"
}

# mem_cli_mem <tier> <out> [mem args...]
#   Runs `bwa-mem3 mem <args>` with BWAMEM3_FORCE_TIER=<tier>, stdout to
#   <out>, stderr to <out>.err. A nonzero exit, or a forced tier the banner
#   does not report, ends the script.
mem_cli_mem() {
    local tier="$1" out="$2"
    shift 2
    if ! BWAMEM3_FORCE_TIER="$tier" "$BWA_MEM3" mem "$@" > "$out" 2> "$out.err"; then
        echo "FAIL: $MEM_CLI_NAME: bwa-mem3 mem $* (tier $tier) exited nonzero:" >&2
        tail -n 20 "$out.err" >&2
        exit 1
    fi
    if ! grep -q "^Executing in $tier mode" "$out.err"; then
        echo "FAIL: $MEM_CLI_NAME: BWAMEM3_FORCE_TIER=$tier was not honoured:" >&2
        head -n 5 "$out.err" >&2
        exit 1
    fi
}

# ---------------------------------------------------------------------------
# Reading SAM
# ---------------------------------------------------------------------------

# mem_cli_records <sam>
#   The alignment records (no header lines).
mem_cli_records() {
    awk '!/^@/' "$1"
}

# mem_cli_primary <sam> <qname> [1|2]
#   The one primary record (neither 0x100 nor 0x800) of read <qname>; with 1
#   or 2, of that end of a pair (0x40 / 0x80). When there is not exactly one,
#   prints nothing on stdout, so an expectation on it fails instead of reading
#   a neighbour's fields, and says how many it found on stderr.
mem_cli_primary() {
    awk -F '\t' -v q="$2" -v end="${3:-}" -v sam="$1" '
        /^@/ || $1 != q { next }
        {
            f = $2 + 0
            if (int(f / 256) % 2 || int(f / 2048) % 2) next
            if (end == "1" && int(f / 64) % 2 == 0) next
            if (end == "2" && int(f / 128) % 2 == 0) next
            rec = $0
            n++
        }
        END {
            if (n == 1) print rec
            else printf "note: %s: %d primary records for %s%s, expected 1\n", sam, n, q, (end == "" ? "" : "/" end) > "/dev/stderr"
        }' "$1"
}

# mem_cli_col <record> <n>
#   Column <n> (1-based) of a SAM record.
mem_cli_col() {
    printf '%s\n' "$1" | awk -F '\t' -v n="$2" '{ print $n }'
}

# mem_cli_tag <record> <TAG>
#   The value of optional field <TAG> (NM:i:3 -> 3), or nothing when absent.
mem_cli_tag() {
    printf '%s\n' "$1" | awk -F '\t' -v t="$2" '{
        for (i = 12; i <= NF; i++)
            if (substr($i, 1, length(t) + 1) == t ":") { print substr($i, length(t) + 4); exit }
    }'
}

# ---------------------------------------------------------------------------
# Checks
# ---------------------------------------------------------------------------

# mem_cli_expect <label> <expected> <actual>
#   String equality; prints `ok: <label>` or `FAIL: <label>: ...` and counts.
mem_cli_expect() {
    MEM_CLI_CHECKS=$((MEM_CLI_CHECKS + 1))
    if [[ $2 == "$3" ]]; then
        echo "ok: $1"
    else
        echo "FAIL: $1: expected '$2', got '$3'" >&2
        MEM_CLI_FAILS=$((MEM_CLI_FAILS + 1))
    fi
}

# mem_cli_same_file <label> <expected-file> <actual-file>
#   Byte equality of two files, header lines starting with @PG removed from
#   both (its CL: records the argv, which differs between runs by design).
#   A file that is missing or holds no alignment record fails the check.
mem_cli_same_file() {
    MEM_CLI_CHECKS=$((MEM_CLI_CHECKS + 1))
    local f
    for f in "$2" "$3"; do
        # Two missing or header-only files would otherwise compare equal.
        if [[ ! -s $f ]] || ! grep -qv '^@' "$f"; then
            echo "FAIL: $1: $f is missing or has no records" >&2
            MEM_CLI_FAILS=$((MEM_CLI_FAILS + 1))
            return 0
        fi
    done
    if diff <(grep -v '^@PG' "$2") <(grep -v '^@PG' "$3") > "$MEM_CLI_WORK/same_file.diff"; then
        echo "ok: $1"
    else
        echo "FAIL: $1: $2 and $3 differ:" >&2
        head -n 20 "$MEM_CLI_WORK/same_file.diff" >&2
        MEM_CLI_FAILS=$((MEM_CLI_FAILS + 1))
    fi
}

# mem_cli_skip <label> <reason>
#   Records a check that could not run here. Never counted as a pass.
mem_cli_skip() {
    MEM_CLI_SKIPS=$((MEM_CLI_SKIPS + 1))
    echo "SKIP: $1: $2"
}
