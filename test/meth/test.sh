#!/usr/bin/env bash
# Regression test: bwa-mem3 mem --meth end-to-end.
#
# Layer 1 (always runs):  valid BAM emission.
#   - binary builds and runs with --meth
#   - produces uncompressed BAM readable by samtools (--bam; --meth alone emits
#     SAM text, so the BGZF assertions below require it explicitly)
#   - @PG ID:bwa-mem3-meth present
#   - BGZF EOF marker at tail
#   - the run stayed PAIRED: every record flagged 0x1, both mates present
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

if [[ ! -x "$BWAMEM3" ]]; then
    echo "ERROR: bwa-mem3 binary not found at $BWAMEM3. Run 'make arm64' first."
    exit 2
fi
if ! command -v "$SAMTOOLS" > /dev/null 2>&1; then
    echo "ERROR: samtools not found on PATH."
    exit 2
fi

cd "$HERE"

# The fixtures are gitignored and are copied in from the bwa-meth checkout --
# see the "Copy bwa-meth fixtures into test/meth/" step in ci.yml. Checked here
# because the indexing step below sends both streams to /dev/null: without
# this, a missing ref.fa surfaced as a bare `exit 2` with no output at all,
# since `set -e` aborts before any of this script's own diagnostics can run.
for fixture in ref.fa t_R1.fastq.gz t_R2.fastq.gz; do
    if [[ ! -f "$fixture" || ! -s "$fixture" ]]; then
        echo "ERROR: fixture $PWD/$fixture is missing, empty, or not a regular file."
        echo "       These are gitignored. Copy them from a bwa-meth checkout:"
        echo "         cp <bwa-meth>/example/{ref.fa,t_R1.fastq.gz,t_R2.fastq.gz} $PWD/"
        exit 2
    fi
done

# Every BAM this test writes goes into a private per-run directory. Fixed names
# under /tmp can be pre-created as symlinks by another local user, in which case
# the redirections below would follow them and clobber an unrelated file; they
# also stop two runs of this script from coexisting.
SCRATCH="$(mktemp -d)"
trap 'rm -rf "$SCRATCH"' EXIT
BAM="$SCRATCH/meth_test.bam"
BAM2="$SCRATCH/meth_test2.bam"

# ---------------------------------------------------------------------------
# Layer 1: BAM emission smoke test
# ---------------------------------------------------------------------------

if [[ ! -f ref.fa.bwameth.c2t.bwt.2bit.64 ]]; then
    "$BWAMEM3" index --meth ref.fa > /dev/null 2>&1
fi

# PAIRED-END on purpose. This ran paired until the D3 rewrite dropped the mate
# while it was retiring Layers 2-3; nothing about that retirement required Layer
# 1 to go single-end, and going single-end took every assertion below with it.
# The tag rules these assertions check are per-strand -- R1 to OT (XR:CT), R2 to
# OB (XR:GA) -- so an SE run exercises one half of the contract and the mate
# half of the XR/XG assignment is never seen on real reads. The directed PE
# regressions (test/regression/meth_pe_placement.sh, meth_mixed_pe_asym.sh)
# cover the placement logic, but on 60 bp hand-written reads over a synthetic
# reference; this is the only --meth test that sees a real read pair.
"$BWAMEM3" mem --meth --bam -t 2 ref.fa t_R1.fastq.gz t_R2.fastq.gz 2> /dev/null > "$BAM"

EXPECT_EOF="1f8b08040000000000ff0600424302001b0003000000000000000000"
ACTUAL_EOF="$(tail -c 28 "$BAM" | od -An -v -t x1 | tr -d ' \n')"
if [[ "${ACTUAL_EOF%$'\n'}" != "${EXPECT_EOF}" ]]; then
    echo "FAIL: BGZF EOF marker mismatch (actual=$ACTUAL_EOF)"
    exit 1
fi

HDR="$("$SAMTOOLS" view -H "$BAM" 2>&1)"
if echo "$HDR" | grep -qi 'truncated\|EOF marker is absent'; then
    echo "FAIL: samtools reports truncated BAM"
    echo "$HDR"
    exit 1
fi
if ! echo "$HDR" | grep -q 'ID:bwa-mem3-meth'; then
    echo "FAIL: @PG ID:bwa-mem3-meth missing"
    exit 1
fi

TOTAL="$("$SAMTOOLS" view -c "$BAM" 2> /dev/null)"
if [[ "$TOTAL" -lt 1 ]]; then
    echo "FAIL: zero records in output BAM"
    exit 1
fi

# Assert the run was actually paired, so the invocation above cannot quietly
# revert to single-end and leave every assertion below still passing over half
# the contract -- which is exactly how the paired coverage was lost the first
# time. Both mates must be present: 0x1 alone would be satisfied by a batch
# where one side failed to load.
PAIRED="$("$SAMTOOLS" view -c -f 0x1 "$BAM" 2> /dev/null)"
R1_COUNT="$("$SAMTOOLS" view -c -f 0x41 "$BAM" 2> /dev/null)"
R2_COUNT="$("$SAMTOOLS" view -c -f 0x81 "$BAM" 2> /dev/null)"
if [[ "$PAIRED" -ne "$TOTAL" ]]; then
    echo "FAIL: $((TOTAL - PAIRED)) of $TOTAL record(s) are not flagged paired (0x1)"
    echo "      -- the --meth run above went single-end"
    exit 1
fi
if [[ "$R1_COUNT" -lt 1 || "$R2_COUNT" -lt 1 ]]; then
    echo "FAIL: expected both mates, got R1=$R1_COUNT R2=$R2_COUNT"
    exit 1
fi

# --chimera-qc classifies read PAIRS, so single-end input makes it a no-op.
"$BWAMEM3" mem --meth --bam --set-as-failed f --chimera-qc \
    ref.fa t_R1.fastq.gz t_R2.fastq.gz 2> /dev/null > "$BAM2"
if [[ ! -s "$BAM2" ]]; then
    echo "FAIL: --set-as-failed + --chimera-qc produced empty output"
    exit 1
fi
# Guarded the same way as the run above, and for the same reason: this is the
# invocation whose flag NEEDS pairs, so it is the one that must not be allowed
# to quietly go single-end while a non-empty BAM keeps the check green.
TOTAL2="$("$SAMTOOLS" view -c "$BAM2" 2> /dev/null)"
PAIRED2="$("$SAMTOOLS" view -c -f 0x1 "$BAM2" 2> /dev/null)"
if [[ "$TOTAL2" -lt 1 ]]; then
    # The -s check above only proves a header was written. Without this, the
    # paired comparison below is 0 -ne 0 and passes over an empty BAM.
    echo "FAIL: zero records in the --chimera-qc output BAM"
    exit 1
fi
if [[ "$PAIRED2" -ne "$TOTAL2" ]]; then
    echo "FAIL: $((TOTAL2 - PAIRED2)) of $TOTAL2 --chimera-qc record(s) are not flagged paired (0x1)"
    echo "      -- --chimera-qc classifies PAIRS, so that run is a no-op single-end"
    exit 1
fi

echo "OK layer 1: bwa-mem3 mem --meth --bam (records=$TOTAL paired, R1=$R1_COUNT R2=$R2_COUNT, BGZF-EOF ok, @PG bwa-mem3-meth ok)"

# --- Bismark XR:Z / XG:Z / XM:Z emission assertions ----------------------
# Every primary mapped record (FLAG & 0x904 == 0) must carry XR:Z:(CT|GA),
# XG:Z:(CT|GA), and XM:Z whose payload length equals SEQ length. Unmapped
# records (FLAG & 0x4) carry XR:Z only. Every record's XR:Z must agree with its
# mate flag. No record emits the legacy Y* tags.

PRIMARY_MAPPED_NO_XR=$("$SAMTOOLS" view -F 0x904 "$BAM" \
    | mawk '!/\tXR:Z:(CT|GA)/{n++} END{print n+0}')
if [[ "$PRIMARY_MAPPED_NO_XR" -ne 0 ]]; then
    echo "FAIL: $PRIMARY_MAPPED_NO_XR primary mapped record(s) missing XR:Z:(CT|GA)"
    exit 1
fi

# XR:Z is the READ conversion, and --meth derives it from the mate alone: R1 is
# C->T, R2 is G->A (src/fastmap.cpp, `yc = is_r2 ? "GA" : "CT"`). The check
# above only asserts XR is one of the two legal values, which stays green even
# if every record were tagged CT -- so it does not actually cover the per-strand
# assignment that going paired here is meant to exercise. This does: it is the
# assertion that fails if the R1/R2 classification regresses (the parity and
# adjacent-name rules that same function uses to pick the mate). Run over EVERY
# record, not just primary mapped, because XR is emitted on all of them.
XR_MATE_BAD=$("$SAMTOOLS" view "$BAM" \
    | mawk '
        {
            want = (int($2 / 64) % 2) ? "CT" : ((int($2 / 128) % 2) ? "GA" : "")
            xr = ""
            for (i = 12; i <= NF; i++) {
                if (substr($i, 1, 5) == "XR:Z:") { xr = substr($i, 6); break }
            }
            if (want == "" || xr != want) { n++ }
        }
        END { print n+0 }')
if [[ "$XR_MATE_BAD" -ne 0 ]]; then
    echo "FAIL: $XR_MATE_BAD record(s) with XR:Z not matching the mate (R1 must be CT, R2 must be GA)"
    exit 1
fi

PRIMARY_MAPPED_NO_XG=$("$SAMTOOLS" view -F 0x904 "$BAM" \
    | mawk '!/\tXG:Z:(CT|GA)/{n++} END{print n+0}')
if [[ "$PRIMARY_MAPPED_NO_XG" -ne 0 ]]; then
    echo "FAIL: $PRIMARY_MAPPED_NO_XG primary mapped record(s) missing XG:Z:(CT|GA)"
    exit 1
fi

XM_LEN_BAD=$("$SAMTOOLS" view -F 0x904 "$BAM" \
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

UNMAPPED_BAD=$("$SAMTOOLS" view -f 0x4 "$BAM" \
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

# Counted with mawk rather than `grep -c ... || true`: grep exits 1 on "no
# match", so the `|| true` needed to tolerate the (expected) zero-match case
# would equally swallow a `samtools view` failure, leaving Y_LEAKS=0 and the
# check passing without ever reading "$BAM". mawk prints 0 on no match and
# exits 0, so pipefail still propagates a samtools failure.
Y_LEAKS=$("$SAMTOOLS" view "$BAM" \
    | mawk '
        {
            for (i = 12; i <= NF; i++) {
                if ($i ~ /^(YS|YC|YD):[ZA]:/) { n++; break }
            }
        }
        END { print n+0 }')
if [[ "$Y_LEAKS" -ne 0 ]]; then
    echo "FAIL: $Y_LEAKS record(s) still emit YS/YC/YD legacy tags"
    exit 1
fi

echo "OK layer 1 Bismark tags: XR/XG/XM well-formed, XR matches mate, no Y* leak"

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
