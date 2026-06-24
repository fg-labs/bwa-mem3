#!/usr/bin/env bash
# Regression test: bwa-mem3 mem --meth end-to-end.
#
# Layer 1 (always runs):  valid BAM emission.
#   - binary builds and runs with --meth
#   - produces uncompressed BAM readable by samtools
#   - @PG ID:bwa-mem3-meth present
#   - BGZF EOF marker at tail
#   - --set-as-failed / --chimera-qc parse cleanly
#
# Layers 2-3 (bwameth structural / byte equivalence) are RETIRED in D3 — see the
# note near the exit at the bottom of this file. D3 intentionally diverges from
# bwameth (the genomic --meth-scoring mode is variant-aware, not collapsed), so
# byte/structural equivalence to bwameth.py is no longer an invariant. D3
# correctness is covered by the CI-wired whole-aligner regressions
# (test/regression/meth_*.sh) plus the directed unit tests.
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
BWAMEM3="$HERE/../../bwa-mem3"
SAMTOOLS="${SAMTOOLS:-samtools}"
BWAMETH_DIR="${BWAMETH_DIR:-$HOME/work/git/bwa-meth}"
BWAMETH_PY="$BWAMETH_DIR/bwameth.py"

if [[ ! -x "$BWAMEM3" ]]; then
    echo "ERROR: bwa-mem3 binary not found at $BWAMEM3. Run 'make arm64' first."
    exit 2
fi
if ! command -v "$SAMTOOLS" >/dev/null 2>&1; then
    echo "ERROR: samtools not found on PATH."
    exit 2
fi

cd "$HERE"

# ---------------------------------------------------------------------------
# Layer 1: BAM emission smoke test
# ---------------------------------------------------------------------------

if [[ ! -f ref.fa.bwameth.c2t.bwt.2bit.64 ]]; then
    "$BWAMEM3" index --meth ref.fa >/dev/null 2>&1
fi

"$BWAMEM3" mem --meth -t 2 ref.fa t_R1.fastq.gz 2>/dev/null > /tmp/meth_test.bam

EXPECT_EOF="1f8b08040000000000ff0600424302001b0003000000000000000000"
ACTUAL_EOF="$(tail -c 28 /tmp/meth_test.bam | od -An -v -t x1 | tr -d ' \n')"
if [[ "${ACTUAL_EOF%$'\n'}" != "${EXPECT_EOF}" ]]; then
    echo "FAIL: BGZF EOF marker mismatch (actual=$ACTUAL_EOF)"; exit 1
fi

HDR="$("$SAMTOOLS" view -H /tmp/meth_test.bam 2>&1)"
if echo "$HDR" | grep -qi 'truncated\|EOF marker is absent'; then
    echo "FAIL: samtools reports truncated BAM"; echo "$HDR"; exit 1
fi
if ! echo "$HDR" | grep -q 'ID:bwa-mem3-meth'; then
    echo "FAIL: @PG ID:bwa-mem3-meth missing"; exit 1
fi

TOTAL="$("$SAMTOOLS" view -c /tmp/meth_test.bam 2>/dev/null)"
if [[ "$TOTAL" -lt 1 ]]; then echo "FAIL: zero records in output BAM"; exit 1; fi

"$BWAMEM3" mem --meth --set-as-failed f --chimera-qc \
    ref.fa t_R1.fastq.gz 2>/dev/null > /tmp/meth_test2.bam
if [[ ! -s /tmp/meth_test2.bam ]]; then
    echo "FAIL: --set-as-failed + --chimera-qc produced empty output"
    exit 1
fi

echo "OK layer 1: bwa-mem3 mem --meth (records=$TOTAL, BGZF-EOF ok, @PG bwa-mem3-meth ok)"

# --- Bismark XR:Z / XG:Z / XM:Z emission assertions ----------------------
# Every primary mapped record (FLAG & 0x904 == 0) must carry XR:Z:(CT|GA),
# XG:Z:(CT|GA), and XM:Z whose payload length equals SEQ length. Unmapped
# records (FLAG & 0x4) carry XR:Z only. No record emits the legacy Y* tags.

PRIMARY_MAPPED_NO_XR=$("$SAMTOOLS" view -F 0x904 /tmp/meth_test.bam \
    | mawk '!/\tXR:Z:(CT|GA)/{n++} END{print n+0}')
if [[ "$PRIMARY_MAPPED_NO_XR" -ne 0 ]]; then
    echo "FAIL: $PRIMARY_MAPPED_NO_XR primary mapped record(s) missing XR:Z:(CT|GA)"
    exit 1
fi

PRIMARY_MAPPED_NO_XG=$("$SAMTOOLS" view -F 0x904 /tmp/meth_test.bam \
    | mawk '!/\tXG:Z:(CT|GA)/{n++} END{print n+0}')
if [[ "$PRIMARY_MAPPED_NO_XG" -ne 0 ]]; then
    echo "FAIL: $PRIMARY_MAPPED_NO_XG primary mapped record(s) missing XG:Z:(CT|GA)"
    exit 1
fi

XM_LEN_BAD=$("$SAMTOOLS" view -F 0x904 /tmp/meth_test.bam \
    | mawk '
        {
            seq_len = length($10)
            xm_len = -1
            for (i = 12; i <= NF; i++) {
                if (substr($i, 1, 5) == "XM:Z:") { xm_len = length($i) - 5; break }
            }
            if (xm_len < 0) { n++; next }
            if (xm_len != seq_len) n++
        }
        END { print n+0 }')
if [[ "$XM_LEN_BAD" -ne 0 ]]; then
    echo "FAIL: $XM_LEN_BAD record(s) with missing or wrong-length XM:Z"
    exit 1
fi

UNMAPPED_BAD=$("$SAMTOOLS" view -f 0x4 /tmp/meth_test.bam \
    | mawk '
        {
            has_xr = 0; has_xg = 0; has_xm = 0
            for (i = 12; i <= NF; i++) {
                if (substr($i, 1, 5) == "XR:Z:") has_xr = 1
                if (substr($i, 1, 5) == "XG:Z:") has_xg = 1
                if (substr($i, 1, 5) == "XM:Z:") has_xm = 1
            }
            if (!has_xr || has_xg || has_xm) n++
        }
        END { print n+0 }')
if [[ "$UNMAPPED_BAD" -ne 0 ]]; then
    echo "FAIL: $UNMAPPED_BAD unmapped record(s) with wrong tag set (XR required, XG/XM forbidden)"
    exit 1
fi

Y_LEAKS=$("$SAMTOOLS" view /tmp/meth_test.bam \
    | grep -cE '\b(YS|YC|YD):[ZA]:' || true)
if [[ "$Y_LEAKS" -ne 0 ]]; then
    echo "FAIL: $Y_LEAKS record(s) still emit YS/YC/YD legacy tags"
    exit 1
fi

echo "OK layer 1 Bismark tags: XR/XG/XM well-formed, no Y* leak"

# ---------------------------------------------------------------------------
# Layers 2-3 (bwameth equivalence) RETIRED in D3.
# ---------------------------------------------------------------------------
# The pre-D3 --meth was a 3-letter / bwameth-style aligner, so Layers 2-3
# asserted structural (QNAME/FLAG/RNAME/POS/CIGAR) and byte equivalence to
# bwameth.py via the doubled c2t reference. D3 redesigns --meth to seed in
# collapsed space but SCORE/EXTEND against the ORIGINAL 4-letter reference with
# a per-strand asymmetric matrix -- it deliberately diverges from bwameth's
# collapsed-space alignment (it distinguishes real variants from conversions).
# Structural/byte equivalence to bwameth is therefore no longer a goal, and the
# old `mem --meth <c2t-prefix> <pre-converted-reads>` invocation does not exist
# in the new dual-index model (bare prefix + original reads + `.meth` seed).
#
# D3 correctness is covered by the CI-wired whole-aligner regressions
# (test/regression/meth_*.sh: dual index + seed->original remap, mixed-PE
# asymmetric scoring, OT/OB x strand placement, original-alphabet output
# integrity, reverse-strand-conversion) plus the directed gamma unit tests.
echo "OK: Layer 1 passed; Layers 2-3 (bwameth structural/byte equivalence) retired in D3 (see comment)."
exit 0
