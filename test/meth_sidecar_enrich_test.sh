#!/usr/bin/env bash
# test/meth_sidecar_enrich_test.sh
#
# Asserts that `mem --meth` enriches its consolidated @SQ block with the
# identity tags (M5/UR/AS/SP) of the *original* (pre-c2t) reference, read from
# that reference's .hdr/.dict sidecar and matched by SN, and forwards the
# sidecar's @CO/@PG/@RG provenance — while NEVER consulting the c2t index's
# own sidecar (which describes the doubled f/r converted contigs).
#
# Why this matters: the default (--bam) and SAM paths forward the index
# sidecar, so their @SQ carries M5/UR that strict validators (GATK, Picard)
# use to bind a BAM to its reference. --meth builds a consolidated @SQ from
# the chrom map (SN/LN only); without this enrichment that provenance is lost.
#
# The test is isolated: it copies phix.fa to a temp dir and builds the c2t
# index there, so the shared fixtures (and other tests' meth headers) are
# untouched.
#
# Usage: test/meth_sidecar_enrich_test.sh <bwa-mem3-binary> <fixtures-dir>

set -euo pipefail

if [[ $# -ne 2 ]]; then
    echo "usage: $0 <bwa-mem3-binary> <fixtures-dir>" >&2
    exit 2
fi

bin="$1"
fixtures="$2"
src_ref="$fixtures/phix.fa"
reads="$fixtures/reads.fa"

[[ -x "$bin" ]]     || { echo "FAIL: bwa-mem3 binary not executable at $bin" >&2; exit 1; }
[[ -s "$src_ref" ]] || { echo "FAIL: phix.fa missing at $src_ref" >&2; exit 1; }
[[ -s "$reads" ]]   || { echo "FAIL: reads.fa missing at $reads" >&2; exit 1; }

fail() { echo "FAIL [meth sidecar]: $*" >&2; exit 1; }

# Reading the BAM header back needs samtools; skip cleanly without it (CI has
# it, so the regression stays enforced there).
if ! command -v samtools >/dev/null 2>&1; then
    echo "SKIP: samtools not found; --meth @SQ enrichment not checked" >&2
    echo "PASS: (skipped, no samtools)"
    exit 0
fi

tmp="$(mktemp -d)"
trap 'rm -rf "$tmp"' EXIT

ref="$tmp/phix.fa"
cp "$src_ref" "$ref"

SN="NC_001422.1"                                  # phiX contig name
M5="0123456789abcdef0123456789abcdef"             # distinctive, not the real md5
UR="file:/refs/original/phix.fa"
AS="phiX-test-asm"
SP="Enterobacteria phage phiX174"
CO="provenance: original-reference sidecar"
PG_ID="orig-sidecar"                              # @PG provenance to forward
RG_ID="orig-rg"                                   # @RG provenance to forward
RG_SM="orig-sample"

LN="5386"                                         # phiX length (= chrom-map LN)

# Build the c2t index once (this is what `mem --meth` actually aligns against).
"$bin" index --meth "$ref" >/dev/null 2>&1 \
    || fail "bwa-mem3 index --meth on phix.fa failed"

# DECOY sidecar on the c2t index. This is the exact path the *wrong* code
# (loading the sidecar from the c2t prefix) would read: "<c2t-prefix>.hdr". It
# carries converted-contig SN and bogus identity tags; if any of it reaches the
# output header, the writer consulted the c2t sidecar — a bug.
printf '%s\n' \
    $'@HD\tVN:1.6\tSO:unsorted' \
    $'@SQ\tSN:fNC_001422.1\tLN:5386\tM5:ffffffffffffffffffffffffffffffff\tUR:file:/WRONG/c2t.fa' \
    $'@CO\tDECOY c2t sidecar must not be forwarded' \
    > "$ref.bwameth.c2t.hdr"

# sq_for_sn <header-file> -> the @SQ line for $SN (robust to field order: SN
# may be the last field, so match a tab-or-end boundary).
sq_for_sn() { grep $'^@SQ\t' "$1" | grep -E $'\tSN:'"$SN"$'(\t|$)' || true; }

# run_meth <dict-LN> -> writes the --meth BAM header for a phix.dict whose @SQ
# carries the given LN. bwa_load_hdr_from_index("$tmp/phix.fa") resolves to
# "$tmp/phix.dict" (it strips the .fa suffix when no .hdr is present).
run_meth() {
    printf '%s\n' \
        $'@HD\tVN:1.6\tSO:unsorted' \
        "@SQ"$'\t'"SN:$SN"$'\t'"LN:$1"$'\t'"M5:$M5"$'\t'"UR:$UR"$'\t'"AS:$AS"$'\t'"SP:$SP" \
        "@CO"$'\t'"$CO" \
        "@PG"$'\t'"ID:$PG_ID"$'\t'"PN:$PG_ID"$'\t'"VN:1.0" \
        "@RG"$'\t'"ID:$RG_ID"$'\t'"SM:$RG_SM" \
        > "$tmp/phix.dict"
    "$bin" mem --meth "$ref" "$reads" >"$tmp/out.bam" 2>"$tmp/err" \
        || { echo "FAIL [meth sidecar]: mem --meth exited non-zero" >&2; cat "$tmp/err" >&2; exit 1; }
    samtools view -H "$tmp/out.bam" >"$tmp/hdr" 2>/dev/null \
        || fail "samtools could not read --meth BAM header"
}

# === positive: sidecar LN matches the chrom map -> @SQ is enriched ==========
run_meth "$LN"
sq="$(sq_for_sn "$tmp/hdr")"
[[ -n "$sq" ]] || { echo "FAIL [meth sidecar]: no @SQ for SN:$SN" >&2
                    grep '^@SQ' "$tmp/hdr" | sed $'s/\t/<TAB>/g' >&2; exit 1; }

want_tag() {
    grep -qF -- "$1" <<<"$sq" \
        || { echo "FAIL [meth sidecar]: @SQ missing '$1'" >&2
             echo "$sq" | sed $'s/\t/<TAB>/g' >&2; exit 1; }
}
want_tag $'\tM5:'"$M5"
want_tag $'\tUR:'"$UR"
want_tag $'\tAS:'"$AS"
want_tag $'\tSP:'"$SP"
grep -qE $'\tLN:'"$LN"$'(\t|$)' <<<"$sq" || fail "@SQ LN is not $LN (chrom map); got: $(sed $'s/\t/<TAB>/g' <<<"$sq")"

# @CO/@PG/@RG provenance from the original sidecar is forwarded (the writer
# passes through every non-@HD/non-@SQ record, not just @CO).
grep -qF -- "$(printf '@CO\t%s' "$CO")" "$tmp/hdr" || fail "forwarded @CO provenance missing"
grep -qE $'^@PG\t.*\tID:'"$PG_ID"$'(\t|$)' "$tmp/hdr" \
    || grep -qE $'^@PG\tID:'"$PG_ID"$'(\t|$)' "$tmp/hdr" \
    || fail "forwarded @PG provenance (ID:$PG_ID) missing"
grep -qE $'^@RG\t.*ID:'"$RG_ID"$'(\t|$)' "$tmp/hdr" || fail "forwarded @RG provenance (ID:$RG_ID) missing"
grep -qE $'^@RG\t.*SM:'"$RG_SM"$'(\t|$)' "$tmp/hdr" || fail "forwarded @RG SM:$RG_SM missing"

# the c2t decoy sidecar must NOT have been consulted
if grep -q "M5:ffffffffffffffffffffffffffffffff" "$tmp/hdr"; then fail "c2t decoy M5 leaked into output"; fi
if grep -q "file:/WRONG/c2t.fa"                  "$tmp/hdr"; then fail "c2t decoy UR leaked into output"; fi
if grep -q "SN:fNC_001422.1"                     "$tmp/hdr"; then fail "c2t converted contig name leaked into output"; fi
if grep -q "DECOY"                               "$tmp/hdr"; then fail "c2t decoy @CO leaked into output"; fi
echo "OK:   [meth sidecar] @SQ enriched with M5/UR/AS/SP from the original reference; @CO forwarded; c2t sidecar ignored"

# === negative: sidecar LN disagrees -> identity tags are NOT trusted ========
# A stale/foreign sidecar (SN matches, LN does not) must not have its M5/UR
# attached to a sequence they no longer identify. @CO pass-through is
# length-independent and still forwarded; the @SQ keeps the chrom-map LN.
run_meth 9999
sq="$(sq_for_sn "$tmp/hdr")"
[[ -n "$sq" ]] || fail "no @SQ for SN:$SN in LN-mismatch case"
grep -qE $'\tLN:'"$LN"$'(\t|$)' <<<"$sq" || fail "@SQ LN is not $LN (chrom map) in LN-mismatch case"
if grep -qF $'\tM5:'"$M5" <<<"$sq"; then fail "M5 carried despite sidecar LN mismatch (should be skipped)"; fi
if grep -qF $'\tUR:'"$UR" <<<"$sq"; then fail "UR carried despite sidecar LN mismatch (should be skipped)"; fi
echo "OK:   [meth sidecar] LN mismatch skips @SQ identity-tag enrichment"

echo "PASS: --meth carries original-reference @SQ identity tags (SN+LN matched) and @CO provenance"
