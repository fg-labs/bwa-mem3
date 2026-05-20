#!/usr/bin/env bash
# test/regression/supp_rep_hard_cap.sh
#
# Regression: --supp-rep-hard-cap N must force MAPQ=0 on supplementary
# alignments whose chain originates from a repetitive seed.
#
# Before #101, the SMEM->seed materialization loop in bwamem.cpp left
# s.n_hits at its default constructor value (1), so chain_n_hits was
# always 1, the gate `chain_n_hits >= cap` never fired, and the option
# shipped as a silent no-op. No test in test/ exercised the flag against
# a workload with repetitive seeds, which is how the regression slipped
# past CI. This script closes that gap.
#
# Fixture (deterministic, regenerated each run from build_fixture.awk;
# bytes are sliced from committed phix.fa so the output is byte-identical
# across awk implementations — no PRNG):
#   - reference with a 60 bp motif duplicated 8x; one copy is flanked by
#     a 200 bp UNIQUE_C that disambiguates the chain extension
#   - 30 SE reads = [400 bp UNIQUE_A][60 bp MOTIF][200 bp UNIQUE_C]
#   - primary chain anchors on UNIQUE_A (chain_n_hits=1, MAPQ=60)
#   - supp chain spans MOTIF+UNIQUE_C: high natural MAPQ from local score
#     margin, but chain_n_hits=8 from the re-seeded sub-MEM within the
#     MOTIF+UNIQUE_C SMEM (the MOTIF alone has SA-count 8)
#
# With --supp-rep-hard-cap 0 every supp keeps its natural MAPQ (60).
# With --supp-rep-hard-cap 3 the gate fires and every supp goes to MAPQ=0.
# Pre-#101 the two SAMs are identical because chain_n_hits is stuck at 1.
#
# Inputs:
#   BWA_MEM3              — path to bwa-mem3 binary under test
#   SUPP_REP_FIXTURE_AWK  — path to test/fixtures/supp_rep/build_fixture.awk
#   SUPP_REP_PHIX_FA      — path to test/fixtures/phix.fa (the fixture source)
#
# Asserts both directions of the gate so a future regression in either
# branch (cap=0 silently zeroing supps, or cap>0 failing to fire) shows up.

set -euo pipefail

: "${BWA_MEM3:?BWA_MEM3 must be set}"
: "${SUPP_REP_FIXTURE_AWK:?SUPP_REP_FIXTURE_AWK must be set}"
: "${SUPP_REP_PHIX_FA:?SUPP_REP_PHIX_FA must be set}"

OUT_DIR="$(mktemp -d -t bwamem3-supprep-XXXXXX)"
trap 'rm -rf "$OUT_DIR"' EXIT

REF="$OUT_DIR/supp_rep_ref.fa"
READS="$OUT_DIR/supp_rep_reads.fq"
CAP0_SAM="$OUT_DIR/cap0.sam"
CAP_N_SAM="$OUT_DIR/cap_n.sam"
# Cap value: must lie in (1, N_COPIES] where N_COPIES=8 in the fixture.
# 3 leaves a comfortable margin and matches the value documented in #104.
CAP=3
# Slack for any future seed-strategy or chain-merge tweaks that legitimately
# split one of the 30 reads into a non-supp shape. Set well below 30 so the
# test still fails loudly if PR #101's one-liner is reverted (which would
# drop all 30 cap=N MAPQ=0 supps to 0).
MIN_SUPP=20

mawk -v MODE=ref   -f "$SUPP_REP_FIXTURE_AWK" "$SUPP_REP_PHIX_FA" > "$REF"
mawk -v MODE=reads -f "$SUPP_REP_FIXTURE_AWK" "$SUPP_REP_PHIX_FA" > "$READS"

"$BWA_MEM3" index "$REF" > "$OUT_DIR/index.log" 2>&1

"$BWA_MEM3" mem -t 1 "$REF" "$READS" \
    > "$CAP0_SAM" 2>"$OUT_DIR/cap0.log"
"$BWA_MEM3" mem -t 1 --supp-rep-hard-cap "$CAP" "$REF" "$READS" \
    > "$CAP_N_SAM" 2>"$OUT_DIR/cap_n.log"

# Count supplementary alignments (SAM flag bit 0x800 / 2048) and split each
# count by MAPQ==0 vs MAPQ>0.
count_supps()        { samtools view -c -f 2048 "$1"; }
count_supps_mapq0()  { samtools view    -f 2048 "$1" | mawk '$5 == 0' | wc -l | tr -d ' '; }
count_supps_mapqpos(){ samtools view    -f 2048 "$1" | mawk '$5 >  0' | wc -l | tr -d ' '; }

cap0_supps=$(count_supps        "$CAP0_SAM")
cap0_mapq0=$(count_supps_mapq0  "$CAP0_SAM")
cap0_mapqp=$(count_supps_mapqpos "$CAP0_SAM")
capn_supps=$(count_supps        "$CAP_N_SAM")
capn_mapq0=$(count_supps_mapq0  "$CAP_N_SAM")
capn_mapqp=$(count_supps_mapqpos "$CAP_N_SAM")

echo "cap=0:      $cap0_supps supps ($cap0_mapqp at MAPQ>0, $cap0_mapq0 at MAPQ=0)"
echo "cap=$CAP:      $capn_supps supps ($capn_mapqp at MAPQ>0, $capn_mapq0 at MAPQ=0)"

# Sanity: the workload must actually produce supplementaries. If it doesn't,
# the rest of the assertions trivially "pass" with zero counts. A fixture
# regression (e.g. read-length tweak that stops splitting the read) shows
# up here, not as a misleading silent pass.
if [ "$cap0_supps" -lt "$MIN_SUPP" ]; then
    echo "FAIL: cap=0 produced $cap0_supps supplementary records, expected >= $MIN_SUPP" >&2
    echo "      fixture or aligner behavior changed in a way that breaks the test setup" >&2
    exit 1
fi

# Direction 1: with cap disabled, the natural MAPQ on the supp chain is
# high (UNIQUE_C disambiguates the chain to one position). If supps come
# back at MAPQ=0 here, either the fixture stopped triggering the
# disambiguation or supp MAPQ computation regressed unrelatedly.
if [ "$cap0_mapq0" -gt 0 ]; then
    echo "FAIL: cap=0 produced $cap0_mapq0 supps at MAPQ=0 (expected 0; fixture should yield high natural MAPQ)" >&2
    exit 1
fi

# Direction 2: with cap=N, every supp goes to MAPQ=0. This is the assertion
# that fails when PR #101's one-line fix is reverted — without it,
# chain_n_hits stays at 1 and the gate never fires.
if [ "$capn_mapq0" -lt "$MIN_SUPP" ]; then
    echo "FAIL: cap=$CAP produced only $capn_mapq0 supps at MAPQ=0 (expected >= $MIN_SUPP)" >&2
    echo "      probable cause: chain_n_hits is not being propagated from SMEM SA-count" >&2
    echo "      (this is the failure mode #101 fixed at src/bwamem.cpp:944-948)" >&2
    exit 1
fi

# Belt-and-suspenders: the supp count should be the same with and without
# the cap. The gate only changes MAPQ; it should not add or drop records.
if [ "$cap0_supps" != "$capn_supps" ]; then
    echo "FAIL: supp record count differs cap=0 ($cap0_supps) vs cap=$CAP ($capn_supps)" >&2
    exit 1
fi

echo "PASS: --supp-rep-hard-cap forces supp MAPQ=0 on repetitive-seed chains ($capn_mapq0 of $capn_supps flipped)"
