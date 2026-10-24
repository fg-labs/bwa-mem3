#!/usr/bin/env bash
# test/framework/mem_cli.sh
#
# Source-only helpers for the `bwa-mem3 mem` CLI regression scripts
# (test/regression/mem_cli_basics.sh and the mem_*_cli.sh scripts beside it).
# Nothing here runs on its own; a script sources it with
#
#     # shellcheck source=test/framework/mem_cli.sh
#     source "$(cd "$(dirname "${BASH_SOURCE[0]}")/../framework" && pwd)/mem_cli.sh"
#
# and calls mem_cli_init first. Each script still owns its PASS:/FAIL: lines
# and its exit status (test/regression/README.md); these helpers only print
# per-check `ok:`, `FAIL:` and `SKIP:` lines and count them.
#
# Portability: bash 3.2 (macOS /bin/bash), BSD awk, mawk and gawk, BSD and GNU
# userland. So no mapfile, no associative arrays, no awk extensions (and(),
# gensub, strtonum), no `sed -i`, no md5sum/sha256sum (cksum is POSIX).
#
# Environment read here:
#   BWA_MEM3        the binary under test (required)
#   MEM_CLI_TIERS   optional space-separated SIMD tier list overriding
#                   mem_cli_tiers' detection (e.g. "avx2" to run one tier)
#   MEM_CLI_KEEP    optional; 1 keeps the work directory and prints its path
#   OMP_NUM_THREADS defaults to 4 when unset, so `index` cannot claim a large
#                   host's every core
#   TMPDIR          where the work directory is created (default /tmp)

# The POSIX cksum ("<crc> <bytes>") of test/fixtures/mem_synth_ref.awk's output.
# Every expected coordinate in the mem_*_cli.sh scripts is derived from that
# reference, so a generator change, or an awk that computes it differently, must
# fail here rather than move the expected values silently.
MEM_CLI_REF_CKSUM="1290712259 12835"

# mem_cli_init <script-name>
#   Checks BWA_MEM3, makes the work directory MEM_CLI_WORK (removed on exit
#   unless MEM_CLI_KEEP=1), sets MEM_CLI_REPO and zeroes the check counters.
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

# mem_cli_synth_ref <out.fa>
#   Writes the synthetic reference (test/fixtures/mem_synth_ref.awk; its header
#   lists the contigs and the planted repeats) and checks the pinned cksum.
mem_cli_synth_ref() {
    awk -f "$MEM_CLI_REPO/test/fixtures/mem_synth_ref.awk" > "$1"
    local sum
    sum="$(cksum < "$1" | awk '{ print $1, $2 }')"
    if [[ $sum != "$MEM_CLI_REF_CKSUM" ]]; then
        echo "FAIL: $MEM_CLI_NAME: synthetic reference cksum is '$sum', pinned '$MEM_CLI_REF_CKSUM'" >&2
        echo "      (test/fixtures/mem_synth_ref.awk changed, or this awk computes it differently)" >&2
        exit 1
    fi
}

# mem_cli_setup_synth
#   Generates and indexes the synthetic reference once per script, at
#   MEM_CLI_SYNTH ($MEM_CLI_WORK/synth/mem_synth.fa).
mem_cli_setup_synth() {
    MEM_CLI_SYNTH="$MEM_CLI_WORK/synth/mem_synth.fa"
    if [[ ! -s $MEM_CLI_SYNTH.bwt.2bit.64 ]]; then
        mkdir -p "$MEM_CLI_WORK/synth"
        mem_cli_synth_ref "$MEM_CLI_SYNTH"
        mem_cli_index "$MEM_CLI_SYNTH"
    fi
}

# mem_cli_setup_phix
#   Copies the committed phiX reference (test/fixtures/phix.fa, one contig
#   "NC_001422.1", 5386 bp) and indexes it once per script, at MEM_CLI_PHIX.
mem_cli_setup_phix() {
    MEM_CLI_PHIX="$MEM_CLI_WORK/phix/phix.fa"
    if [[ ! -s $MEM_CLI_PHIX.bwt.2bit.64 ]]; then
        mkdir -p "$MEM_CLI_WORK/phix"
        cp "$MEM_CLI_REPO/test/fixtures/phix.fa" "$MEM_CLI_PHIX"
        mem_cli_index "$MEM_CLI_PHIX"
    fi
}

# ---------------------------------------------------------------------------
# Sequences and reads
# ---------------------------------------------------------------------------

# mem_cli_seq <ref.fa> <contig> <pos> <len>
#   The reference bases [pos, pos+len) of a contig, 1-based, as stored.
mem_cli_seq() {
    awk -v c="$2" -v p="$3" -v l="$4" '
        /^>/ { on = (substr($1, 2) == c); next }
        on { s = s $0 }
        END {
            if (p < 1 || length(s) < p + l - 1) exit 1
            print substr(s, p, l)
        }' "$1"
}

# mem_cli_revcomp <seq>
mem_cli_revcomp() {
    printf '%s\n' "$1" | rev | tr 'ACGTNacgtn' 'TGCANtgcan'
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
#   The SIMD tiers `mem` can run on this host, one per line, lowest first: on
#   x86 avx2 and, where the host has it, avx512bw (avx2 is the runtime floor;
#   the sub-AVX2 tiers refuse `mem`); on arm64 neon. MEM_CLI_TIERS overrides.
mem_cli_tiers() {
    if [[ -n ${MEM_CLI_TIERS:-} ]]; then
        local forced
        read -r -a forced <<< "$MEM_CLI_TIERS"
        printf '%s\n' "${forced[@]}"
        return 0
    fi
    local runtime
    runtime="$(env -u BWAMEM3_FORCE_TIER "$BWA_MEM3" version 2> /dev/null \
        | sed -n 's/^SIMD runtime: \([a-z0-9]*\).*/\1/p')"
    case "$runtime" in
        avx512bw) printf 'avx2\navx512bw\n' ;;
        avx2 | neon) echo "$runtime" ;;
        *)
            echo "FAIL: $MEM_CLI_NAME: cannot read the runtime SIMD tier from '$BWA_MEM3 version' (got '$runtime')" >&2
            exit 1
            ;;
    esac
}

# mem_cli_mem <tier> <out> [mem args...]
#   Runs `bwa-mem3 mem <args>` with BWAMEM3_FORCE_TIER=<tier> ("" leaves the
#   tier to the host), stdout to <out>, stderr to <out>.err. A nonzero exit,
#   or a forced tier the banner does not report, ends the script.
mem_cli_mem() {
    local tier="$1" out="$2"
    shift 2
    if ! BWAMEM3_FORCE_TIER="$tier" "$BWA_MEM3" mem "$@" > "$out" 2> "$out.err"; then
        echo "FAIL: $MEM_CLI_NAME: bwa-mem3 mem $* (tier ${tier:-host}) exited nonzero:" >&2
        tail -n 20 "$out.err" >&2
        exit 1
    fi
    if [[ -n $tier ]] && ! grep -q "^Executing in $tier mode" "$out.err"; then
        echo "FAIL: $MEM_CLI_NAME: BWAMEM3_FORCE_TIER=$tier was not honoured:" >&2
        head -n 5 "$out.err" >&2
        exit 1
    fi
}

# mem_cli_mem_status <tier> <out> [mem args...]
#   As mem_cli_mem, but for a run that may fail: sets MEM_CLI_STATUS to mem's
#   exit status instead of ending the script.
# MEM_CLI_STATUS is read by the calling script, which shellcheck sees separately.
# shellcheck disable=SC2034
mem_cli_mem_status() {
    local tier="$1" out="$2"
    shift 2
    MEM_CLI_STATUS=0
    BWAMEM3_FORCE_TIER="$tier" "$BWA_MEM3" mem "$@" > "$out" 2> "$out.err" || MEM_CLI_STATUS=$?
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
#   or 2, of that end of a pair (0x40 / 0x80). Prints nothing when there is
#   not exactly one, so an expectation on it fails instead of reading a
#   neighbour's fields.
mem_cli_primary() {
    awk -F '\t' -v q="$2" -v end="${3:-}" '
        /^@/ || $1 != q { next }
        {
            f = $2 + 0
            if (int(f / 256) % 2 || int(f / 2048) % 2) next
            if (end == "1" && int(f / 64) % 2 == 0) next
            if (end == "2" && int(f / 128) % 2 == 0) next
            rec = $0
            n++
        }
        END { if (n == 1) print rec }' "$1"
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

# mem_cli_count <sam> <awk-condition>
#   How many records satisfy an awk condition over the tab-split fields, e.g.
#   'int($2 / 2) % 2 == 1' (proper pairs; int() division, not and(), which BSD
#   awk and mawk lack).
mem_cli_count() {
    awk -F '\t' '!/^@/ && ('"$2"') { n++ } END { print n + 0 }' "$1"
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
mem_cli_same_file() {
    MEM_CLI_CHECKS=$((MEM_CLI_CHECKS + 1))
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
