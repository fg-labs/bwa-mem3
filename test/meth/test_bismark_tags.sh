#!/usr/bin/env bash
# Layer 3: holodeck golden-truth XR/XG/XM diff.
#
# Builds holodeck from a pinned SHA (nh/bisulfite-conversion HEAD as of the
# spec date), simulates paired-end reads with --methylation-mode em-seq, runs
# `bwa-mem3 mem --meth` on the simulated FASTQ, and stream-diffs the
# (qname, flag, XR, XG, XM) tuples on primary alignments. Expected divergence
# on a clean simulator + bwa-mem3 run is 0 records.
#
# Skip if cargo is not on PATH (needed to build holodeck).
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
BWAMEM3="$HERE/../../bwa-mem3"
SAMTOOLS="${SAMTOOLS:-samtools}"

# Pinned holodeck SHA (HEAD of nh/bisulfite-conversion at spec date).
HOLODECK_SHA="d240816c8c2dfbd72bdb3386679d0ad4d4f223e6"
HOLODECK_REPO="https://github.com/fg-labs/holodeck.git"
HOLODECK_DIR="${HOLODECK_DIR:-/tmp/holodeck-${HOLODECK_SHA:0:8}}"
HOLODECK_BIN="${HOLODECK_DIR}/target/release/holodeck"

if [[ ! -x "$BWAMEM3" ]]; then
    echo "ERROR: bwa-mem3 not built at $BWAMEM3 (run 'make' or 'make arm64' first)." >&2
    exit 2
fi
if ! command -v "$SAMTOOLS" >/dev/null 2>&1; then
    echo "ERROR: samtools not found on PATH." >&2
    exit 2
fi
if ! command -v cargo >/dev/null 2>&1; then
    echo "SKIP layer 3: cargo not on PATH (needed to build holodeck)"
    exit 0
fi

if [[ ! -x "$HOLODECK_BIN" ]]; then
    echo "[layer 3] cloning + building holodeck @ $HOLODECK_SHA ..."
    rm -rf "$HOLODECK_DIR"
    git clone --filter=blob:none "$HOLODECK_REPO" "$HOLODECK_DIR" >/dev/null 2>&1
    git -C "$HOLODECK_DIR" checkout "$HOLODECK_SHA" >/dev/null 2>&1
    (cd "$HOLODECK_DIR" && cargo build --release --quiet)
fi

cd "$HERE"

# Refresh the bwa-mem3 c2t index for the existing test fixture.
if [[ ! -f ref.fa.bwameth.c2t.bwt.2bit.64 ]]; then
    "$BWAMEM3" index --meth ref.fa >/dev/null 2>&1
fi

# Index the original ref.fa with samtools faidx (holodeck needs .fai).
[[ -f ref.fa.fai ]] || "$SAMTOOLS" faidx ref.fa

PREFIX="/tmp/holodeck-bismark-tags"
rm -f "$PREFIX".* "$PREFIX".golden.bam

"$HOLODECK_BIN" simulate \
    -r ref.fa \
    -o "$PREFIX" \
    --methylation-mode em-seq \
    --methylation-rate 0.7 \
    --methylation-conversion-rate 0.99 \
    --golden-bam \
    --simple-names \
    --seed 42 \
    -c 30 -l 150 -t 2 \
    >/dev/null 2>&1

GOLDEN_BAM="${PREFIX}.golden.bam"
R1="${PREFIX}.r1.fastq.gz"
R2="${PREFIX}.r2.fastq.gz"

if [[ ! -s "$GOLDEN_BAM" || ! -s "$R1" || ! -s "$R2" ]]; then
    echo "FAIL layer 3: holodeck did not produce expected outputs"
    ls -la "${PREFIX}".* 2>&1 | head -10
    exit 1
fi

# Run bwa-mem3 on the simulated FASTQs.
MINE_BAM="/tmp/bismark-tags-mine.bam"
"$BWAMEM3" mem --meth -t 2 ref.fa "$R1" "$R2" 2>/dev/null > "$MINE_BAM"

# Sort both by qname so per-qname comparison is straightforward.
"$SAMTOOLS" sort -n -@ 2 -O bam -o "${PREFIX}.golden.qsort.bam"  "$GOLDEN_BAM"
"$SAMTOOLS" sort -n -@ 2 -O bam -o "${PREFIX}.mine.qsort.bam"    "$MINE_BAM"

# Extract (qname, flag, pos, cigar, XR, XG, XM) per primary record
# (flag & 0x904 == 0: mapped, not secondary, not supplementary). The pos+cigar
# pair lets us tell which records have identical alignments — bwa-mem3 may
# extend the alignment past where holodeck soft-clipped (or vice versa), and
# in those cases the XM strings will legitimately differ at the ends. We
# require XR/XG to match always, and XM to match exactly only on records
# with matching (pos, cigar).
extract_tags() {
    "$SAMTOOLS" view "$1" \
        | mawk -F '\t' '
            # Bitwise AND for mawk (no built-in band, no hex literals).
            # 0x904 == 2308 (unmapped|secondary|supplementary).
            # 0x40  == 64   (first segment / R1).
            function bit_any(a, mask,    rest, m, abit, mbit) {
                rest = a
                m = mask
                while (m > 0) {
                    abit = rest % 2
                    mbit = m % 2
                    if (abit == 1 && mbit == 1) return 1
                    rest = int(rest / 2)
                    m    = int(m / 2)
                }
                return 0
            }
            {
                if (bit_any($2, 2308)) next
                xr = ""; xg = ""; xm = ""
                for (i = 12; i <= NF; i++) {
                    if (substr($i, 1, 5) == "XR:Z:") xr = substr($i, 6)
                    if (substr($i, 1, 5) == "XG:Z:") xg = substr($i, 6)
                    if (substr($i, 1, 5) == "XM:Z:") xm = substr($i, 6)
                }
                is_r1 = bit_any($2, 64)
                printf "%s\t%d\t%s\t%s\t%s\t%s\t%s\n",
                       $1, is_r1, $4, $6, xr, xg, xm
            }
        ' \
        | sort -k1,1 -k2,2n
}

extract_tags "${PREFIX}.golden.qsort.bam" > "${PREFIX}.golden.tags"
extract_tags "${PREFIX}.mine.qsort.bam"   > "${PREFIX}.mine.tags"

# Join on (qname, is_r1) and assert XR/XG always match, XM matches when
# (pos, cigar) match.
mawk_check=$(mawk -F '\t' '
    BEGIN {
        # Read golden into associative arrays keyed by qname\tis_r1.
        while ((getline line < "'"${PREFIX}.golden.tags"'") > 0) {
            split(line, f, "\t")
            key = f[1] "\t" f[2]
            g_pos[key]   = f[3]
            g_cigar[key] = f[4]
            g_xr[key]    = f[5]
            g_xg[key]    = f[6]
            g_xm[key]    = f[7]
        }
        close("'"${PREFIX}.golden.tags"'")
    }
    {
        key = $1 "\t" $2
        if (!(key in g_pos)) { uniq_mine++; next }
        n_total++
        if ($5 != g_xr[key]) { fail_xr++; if (fail_xr <= 3) print "XR diff " key " mine=" $5 " gold=" g_xr[key] }
        if ($6 != g_xg[key]) { fail_xg++; if (fail_xg <= 3) print "XG diff " key " mine=" $6 " gold=" g_xg[key] }
        if ($3 == g_pos[key] && $4 == g_cigar[key]) {
            n_aln_match++
            if ($7 == g_xm[key]) {
                n_xm_match++
            } else {
                if (fail_xm <= 3) {
                    print "XM diff (same aln) " key " pos=" $3 " cigar=" $4
                    print "  mine: " $7
                    print "  gold: " g_xm[key]
                }
                fail_xm++
            }
        } else {
            n_aln_diverge++
        }
        # Same length always.
        if (length($7) != length(g_xm[key])) fail_xmlen++
    }
    END {
        printf "%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",
               n_total+0, n_aln_match+0, n_xm_match+0,
               fail_xr+0, fail_xg+0, fail_xm+0, fail_xmlen+0, n_aln_diverge+0
    }
' "${PREFIX}.mine.tags" 2>&1)

# Last line of mawk_check is the tab-separated summary; everything before is diags.
SUMMARY=$(echo "$mawk_check" | tail -1)
DIAGS=$(echo "$mawk_check" | sed '$d')
read -r N_TOTAL N_ALN_MATCH N_XM_MATCH FAIL_XR FAIL_XG FAIL_XM FAIL_XMLEN N_ALN_DIVERGE <<<"$(echo "$SUMMARY" | tr '\t' ' ')"

echo "Layer 3 summary:"
echo "  records compared:             $N_TOTAL"
echo "  alignment match (pos+cigar):  $N_ALN_MATCH"
echo "  alignment diverged:           $N_ALN_DIVERGE"
echo "  XM match given matching aln:  $N_XM_MATCH / $N_ALN_MATCH"
echo "  XR mismatches:                $FAIL_XR"
echo "  XG mismatches:                $FAIL_XG"
echo "  XM mismatches (aln match):    $FAIL_XM"
echo "  XM length mismatches:         $FAIL_XMLEN"

if [[ "$FAIL_XR" -ne 0 || "$FAIL_XG" -ne 0 || "$FAIL_XMLEN" -ne 0 ]]; then
    echo "FAIL layer 3: XR/XG/XM-length mismatches"
    echo "$DIAGS" | head -20
    exit 1
fi

if [[ "$FAIL_XM" -ne 0 ]]; then
    echo "FAIL layer 3: XM diverges from holodeck golden on records with matching alignment"
    echo "$DIAGS" | head -20
    exit 1
fi

echo "OK layer 3: XR/XG/XM Bismark golden-truth match (records=$N_TOTAL, aln-match=$N_ALN_MATCH, XM-match=$N_XM_MATCH)"
