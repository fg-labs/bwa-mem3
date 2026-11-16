#!/usr/bin/env bash
# test/meth_rg_header_test.sh
#
# Asserts that a `-R '@RG\t...'` read group supplied on the command line is
# emitted as an @RG header line in BOTH the default path and `--meth` mode.
#
# Bug: `mem --meth` stamps every record with the RG:Z tag derived from the
# user's -R argument (via the global bwa_rg_id), but the meth BAM writer
# (meth_bam_writer_open) never received the assembled hdr_line, so no @RG
# header line was written. The result is a malformed BAM whose records
# reference a read group ID that the header never declares — strict
# validators (Picard, noodles) reject it. The default (non-meth) path
# threads hdr_line into its writer and was always correct; this test pins
# both so the two modes can never diverge again.
#
# The assertion is the SAM invariant "every RG:Z tag value must be declared
# by an @RG ID: in the header":
#   * the @RG header line is present and carries ID:<rg>,
#   * every alignment record carries RG:Z:<rg>, and
#   * the RG:Z value matches the declared @RG ID.
#
# The --meth output is read back with samtools (it parses the SAM or BAM container). When
# samtools is absent (minimal dev boxes) the --meth half is skipped with a
# warning rather than failing; CI always has samtools, so the regression
# stays enforced there.
#
# Usage: test/meth_rg_header_test.sh <bwa-mem3-binary> <fixtures-dir>

set -euo pipefail

if [[ $# -ne 2 ]]; then
    echo "usage: $0 <bwa-mem3-binary> <fixtures-dir>" >&2
    exit 2
fi

bin="$1"
fixtures="$2"
ref="$fixtures/phix.fa"
reads="$fixtures/reads.fa"

[[ -x "$bin" ]] || {
    echo "FAIL: bwa-mem3 binary not executable at $bin" >&2
    exit 1
}
[[ -s "$ref" ]] || {
    echo "FAIL: phix.fa missing at $ref" >&2
    exit 1
}
[[ -s "$reads" ]] || {
    echo "FAIL: reads.fa missing at $reads" >&2
    exit 1
}

# A read group whose ID is distinct from every other token so a stray match
# can't accidentally satisfy the assertions.
RG_ID="rg_meth_test"
RG=$'@RG\tID:'"$RG_ID"$'\tSM:sample1\tPL:ILLUMINA'

tmp="$(mktemp -d)"
trap 'rm -rf "$tmp"' EXIT

# assert_rg_consistent <label> <header-text-file> <records-text-file>
# Verifies the SAM invariant: an @RG ID:<RG_ID> line exists, every record
# carries RG:Z:<RG_ID>, and there is at least one record to check.
assert_rg_consistent() {
    local label="$1" hdr="$2" recs="$3"

    # ID:<RG_ID> is always followed by a tab here (SM: trails it in $RG).
    if ! grep -q $'^@RG\t' "$hdr" || ! grep -q $'\tID:'"$RG_ID"$'\t' "$hdr"; then
        echo "FAIL [$label]: no '@RG ... ID:$RG_ID' header line" >&2
        echo "---- header (@RG lines, tabs as <TAB>) ----" >&2
        grep '^@RG' "$hdr" | sed $'s/\t/<TAB>/g' >&2 || echo "(none)" >&2
        echo "-------------------------------------------" >&2
        exit 1
    fi

    local n_recs n_tagged
    n_recs="$(wc -l < "$recs" | tr -d ' ')"
    [[ "$n_recs" -ge 1 ]] || {
        echo "FAIL [$label]: no alignment records emitted" >&2
        exit 1
    }

    # RG:Z:<RG_ID> may be the last tag on a record (no trailing tab) — match
    # a tab-or-end boundary so a record ending in RG:Z still counts.
    n_tagged="$(grep -cE $'\tRG:Z:'"$RG_ID"$'(\t|$)' "$recs" || true)"
    if [[ "$n_tagged" -ne "$n_recs" ]]; then
        echo "FAIL [$label]: $n_tagged/$n_recs records carry RG:Z:$RG_ID" >&2
        exit 1
    fi

    echo "OK:   [$label] @RG ID:$RG_ID declared and all $n_recs records tagged"
}

# --- default (non-meth) path: SAM output -----------------------------------
# This path was always correct; it is the control that proves --meth must
# match it.
if [[ ! -s "$ref.bwt.2bit.64" || ! -s "$ref.amb" ||
    ! -s "$ref.ann" || ! -s "$ref.pac" ]]; then
    "$bin" index "$ref" > /dev/null 2>&1 || {
        echo "FAIL: bwa-mem3 index on phix.fa failed" >&2
        exit 1
    }
fi

"$bin" mem -R "$RG" "$ref" "$reads" > "$tmp/plain.sam" 2> "$tmp/plain.err" || {
    echo "FAIL: bwa-mem3 mem (non-meth) exited non-zero" >&2
    cat "$tmp/plain.err" >&2
    exit 1
}
grep '^@' "$tmp/plain.sam" > "$tmp/plain.hdr" || true
grep -v '^@' "$tmp/plain.sam" > "$tmp/plain.recs" || true
assert_rg_consistent "non-meth" "$tmp/plain.hdr" "$tmp/plain.recs"

# --- --meth path: BAM output (needs samtools to read the header) -----------
if ! command -v samtools > /dev/null 2>&1; then
    echo "SKIP: samtools not found; --meth half not checked" >&2
    echo "PASS: @RG header emitted in non-meth mode (--meth skipped)"
    exit 0
fi

# --meth requires the D3 seed index built with `index --meth`. Check the full
# artifact set (converted seed FASTA + the four seed FMI files), not just one,
# so a partial index left by an interrupted run is rebuilt rather than skipped —
# matching the non-meth block above. The seed has no `.0123` (never read in
# --meth; extension uses the original reference, which itself pac-fetches).
if [[ ! -s "$ref.meth.fa" || ! -s "$ref.meth.bwt.2bit.64" ||
    ! -s "$ref.meth.amb" || ! -s "$ref.meth.ann" ||
    ! -s "$ref.meth.pac" ]]; then
    "$bin" index --meth "$ref" > /dev/null 2>&1 \
        || {
            echo "FAIL: bwa-mem3 index --meth on phix.fa failed" >&2
            exit 1
        }
fi

"$bin" mem --meth -R "$RG" "$ref" "$reads" > "$tmp/meth.bam" 2> "$tmp/meth.err" || {
    echo "FAIL: bwa-mem3 mem --meth exited non-zero" >&2
    cat "$tmp/meth.err" >&2
    exit 1
}
samtools view -H "$tmp/meth.bam" > "$tmp/meth.hdr" 2> /dev/null \
    || {
        echo "FAIL: samtools could not read --meth BAM header" >&2
        exit 1
    }
samtools view "$tmp/meth.bam" > "$tmp/meth.recs" 2> /dev/null \
    || {
        echo "FAIL: samtools could not read --meth BAM records" >&2
        exit 1
    }
assert_rg_consistent "--meth" "$tmp/meth.hdr" "$tmp/meth.recs"

# --- --meth @HD de-dup: a user -H @HD must suppress the default @HD ---------
# The --meth header text (meth_hdr_text) starts with a default @HD
# (BWAMEM3_DEFAULT_HD_LINE, "@HD VN:1.5 SO:unsorted GO:query") UNLESS the
# user's -H/-R text supplies an @HD, whose first one is written there instead.
# Writing both would give two @HD lines (htslib does not de-dup @HD). Pass a
# DISTINCT @HD (SO:coordinate) alongside -R and assert the output carries
# exactly one @HD line and it is the user's — proving the default was
# suppressed and the user's line survived. This pins the policy that the -R
# assertions above never exercise (an @RG line carries no @HD).
USER_HD=$'@HD\tVN:1.6\tSO:coordinate'
"$bin" mem --meth -R "$RG" -H "$USER_HD" "$ref" "$reads" \
    > "$tmp/meth_hd.bam" 2> "$tmp/meth_hd.err" || {
    echo "FAIL: bwa-mem3 mem --meth -H @HD exited non-zero" >&2
    cat "$tmp/meth_hd.err" >&2
    exit 1
}
samtools view -H "$tmp/meth_hd.bam" > "$tmp/meth_hd.hdr" 2> /dev/null \
    || {
        echo "FAIL: samtools could not read --meth -H @HD BAM header" >&2
        exit 1
    }

n_hd="$(grep -c $'^@HD\t' "$tmp/meth_hd.hdr" || true)"
if [[ "$n_hd" -ne 1 ]]; then
    echo "FAIL [--meth @HD dedup]: expected exactly 1 @HD line, got $n_hd" >&2
    echo "---- header @HD lines (tabs as <TAB>) ----" >&2
    grep '^@HD' "$tmp/meth_hd.hdr" | sed $'s/\t/<TAB>/g' >&2 || echo "(none)" >&2
    echo "------------------------------------------" >&2
    exit 1
fi
if ! grep -qE $'^@HD\t.*\tSO:coordinate(\t|$)' "$tmp/meth_hd.hdr"; then
    echo "FAIL [--meth @HD dedup]: the surviving @HD is not the user's (SO:coordinate)" >&2
    grep '^@HD' "$tmp/meth_hd.hdr" | sed $'s/\t/<TAB>/g' >&2
    exit 1
fi
echo "OK:   [--meth @HD dedup] exactly one @HD line, user's SO:coordinate preserved"

# --- --meth with several -H @HD lines: the first wins, as the first line ------
# The SAM spec requires @HD first and allows only one. The --meth header text
# writes the user's first @HD as its first line and drops every later one (the
# default writer's policy). htslib alone would keep a second @HD where it sat
# in the user text, so this pins the text-level de-dup.
"$bin" mem --meth -R "$RG" -H $'@CO\tbefore\n@HD\tVN:1.6\tSO:coordinate\n@HD\tVN:1.5\tSO:unsorted' \
    "$ref" "$reads" > "$tmp/meth_hd2.bam" 2> "$tmp/meth_hd2.err" || {
    echo "FAIL: bwa-mem3 mem --meth with two -H @HD lines exited non-zero" >&2
    cat "$tmp/meth_hd2.err" >&2
    exit 1
}
samtools view -H "$tmp/meth_hd2.bam" > "$tmp/meth_hd2.hdr" 2> /dev/null \
    || {
        echo "FAIL: samtools could not read --meth two-@HD header" >&2
        exit 1
    }
n_hd="$(grep -c $'^@HD\t' "$tmp/meth_hd2.hdr" || true)"
first_line="$(head -n 1 "$tmp/meth_hd2.hdr")"
if [[ "$n_hd" -ne 1 || "$first_line" != $'@HD\tVN:1.6\tSO:coordinate' ]]; then
    echo "FAIL [--meth two @HD]: expected one @HD, the user's first, as line 1" >&2
    grep -n '^@' "$tmp/meth_hd2.hdr" | sed $'s/\t/<TAB>/g' >&2
    exit 1
fi
if [[ "$(grep -c $'^@CO\tbefore$' "$tmp/meth_hd2.hdr" || true)" -ne 1 ]]; then
    echo "FAIL [--meth two @HD]: the user's @CO record before the @HD lines is not kept exactly once" >&2
    exit 1
fi
echo "OK:   [--meth two @HD] one @HD, the user's first, written as line 1"

# --- --meth with an empty line in -H: dropped, as in plain mem -----------------
# An empty line is not a header record. The --meth header text drops it, as the
# default writer does; handed to htslib it was rejected and the run exited 1.
"$bin" mem --meth -H $'@CO	x

@CO	y' "$ref" "$reads" \
    > "$tmp/meth_blank.out" 2> "$tmp/meth_blank.err" || {
    echo "FAIL: bwa-mem3 mem --meth with an empty -H line exited non-zero" >&2
    cat "$tmp/meth_blank.err" >&2
    exit 1
}
samtools view -H "$tmp/meth_blank.out" > "$tmp/meth_blank.hdr" 2> /dev/null \
    || {
        echo "FAIL: samtools could not read --meth empty-line header" >&2
        exit 1
    }
for co in x y; do
    if [[ "$(grep -c "^@CO"$'\t'"$co\$" "$tmp/meth_blank.hdr" || true)" -ne 1 ]]; then
        echo "FAIL [--meth empty -H line]: the @CO record '$co' around it is not kept exactly once" >&2
        exit 1
    fi
done
echo "OK:   [--meth empty -H line] dropped; the records around it kept"

echo "PASS: -R read group emitted as @RG header in both default and --meth modes"
