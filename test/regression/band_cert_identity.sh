#!/usr/bin/env bash
# test/regression/band_cert_identity.sh
#
# Regression: the certified adaptive extension band (on by default) must be
# byte-identical to the full-width extension ladder, and `--no-band-cert` must
# select that full-width ladder as an escape hatch / A-B handle.
#
# The default seed-extension path scores each pair at a narrow probe band first
# and finalizes only the pairs it can PROVE are already optimal there (a per-pair
# certificate plus a clip-decision guard); every other pair falls through to the
# exact full-width ceiling ladder. The claim is that this changes no emitted
# field. `--no-band-cert` disables the probe entirely and runs the full-width
# ladder for every pair, so:
#
#   md5(default)  ==  md5(--no-band-cert)
#
# on the same inputs -- byte-for-byte, @PG excluded (its CL: records the argv,
# which differs between the two runs). If the certificate ever certifies a pair
# whose wider-band result would differ (a wrong `d_max`, a flipped clip decision,
# a polluted retry-ladder `prev`), the two streams diverge and this fails.
#
# This is the flag-isolation + contract smoke test. band_cert's effect on the
# emitted record is subtle -- the final CIGAR/POS is recomputed by mem_reg2aln,
# which re-aligns with its own band, so on simple phix reads a narrow-vs-full-width
# extension difference is largely masked. The load-bearing byte-identity gate for
# band_cert on real reads is therefore (a) the whole-genome md5 reported in the PR
# and (b) chr22_parity.sh, which runs the DEFAULT path (band_cert on) against the
# bwa golden on ~50k real chr22 PE reads on every CI row -- a band_cert regression
# that moved a real read would fail it. This script pins the toggle, the flag
# wiring, and the parameter-envelope gate that those real-read gates do not isolate.
# The whole-aligner ~3% throughput win is a separate benchmark, not asserted here.
#
# Fixture (deterministic, no PRNG): PE reads sliced directly from the committed
# phix.fa -- the same construction compat_byte_identical.sh uses, so the FASTQ is
# byte-identical across awk implementations. read1 is a forward window; read2 is
# the reverse complement of a downstream window, so every pair aligns
# concordantly and drives both the extend and the clip-to-query-end finalize
# branches the certificate has to keep invariant.
#
# Inputs (env vars):
#   BWA_MEM3          — path to the bwa-mem3 binary under test
#   BAND_CERT_PHIX_FA — path to test/fixtures/phix.fa (the reference source)
#   BAND_CERT_WORK_DIR — directory for fixture-private intermediates

set -euo pipefail
# Any command that fails unexpectedly (index/aligner/diff/conversion) exits via
# set -e; surface it as a FAIL: record with the failing line and status rather
# than dying silently.
# shellcheck disable=SC2154  # rc is assigned by `rc=$?` at the start of the trap body
trap 'rc=$?; echo "FAIL: unexpected command failure at line $LINENO (exit $rc)" >&2; exit "$rc"' ERR

: "${BWA_MEM3:?BWA_MEM3 must be set}"
: "${BAND_CERT_PHIX_FA:?BAND_CERT_PHIX_FA must be set}"
: "${BAND_CERT_WORK_DIR:?BAND_CERT_WORK_DIR must be set}"

# The non-vacuity deletion-width scan below needs mawk; if it is unavailable,
# skip (per the regression-suite SKIP convention) rather than trip the ERR trap.
if ! command -v mawk > /dev/null 2>&1; then
    echo "SKIP: mawk not on PATH — band_cert byte-identity regression skipped"
    exit 0
fi

mkdir -p "$BAND_CERT_WORK_DIR"

ref="$BAND_CERT_WORK_DIR/phix.fa"
cp "$BAND_CERT_PHIX_FA" "$ref"

# --- Build the PE FASTQs by slicing phix (deterministic, no PRNG). ---
seq="$(grep -v '^>' "$ref" | tr -d '\n' | tr 'acgt' 'ACGT')"
r1="$BAND_CERT_WORK_DIR/r1.fq"
r2="$BAND_CERT_WORK_DIR/r2.fq"
: > "$r1"
: > "$r2"
L=120 # read length
G=300 # inner gap between the two windows (fragment ~ G+L)
n_pairs=0
for off in 200 900 1600 2300 3000 3700; do
    s1="${seq:$off:$L}"
    s2fwd="${seq:$((off + G)):$L}"
    s2="$(printf '%s' "$s2fwd" | rev | tr 'ACGT' 'TGCA')"
    [ "${#s1}" -eq "$L" ] && [ "${#s2}" -eq "$L" ] || continue
    q="$(printf 'I%.0s' $(seq 1 "$L"))"
    printf '@pair%d/1\n%s\n+\n%s\n' "$n_pairs" "$s1" "$q" >> "$r1"
    printf '@pair%d/2\n%s\n+\n%s\n' "$n_pairs" "$s2" "$q" >> "$r2"
    n_pairs=$((n_pairs + 1))
done

if [ "$n_pairs" -lt 3 ]; then
    echo "FAIL: fixture built only $n_pairs pairs (phix slicing regression?)" >&2
    exit 1
fi

# --- Deletion-bearing pairs so the certificate is genuinely EXERCISED. ---
# The exact-match pairs above never need a band wider than the probe (w0=20), so
# their optimal alignment is trivially certifiable and a broken certificate that
# certified everything would still match. These pairs carry a DEL-bp deletion
# (DEL > w0): read1 = ref[off..off+L1] + ref[off+L1+DEL .. off+L1+DEL+L2], so its
# optimal alignment is L1 M / DEL D / L2 M and REQUIRES a band >= DEL to find. At
# the narrow probe the deletion is unreachable, so a correct certificate declines
# the pair (falls through to the exact full-width ladder); a certificate that
# wrongly certified it narrow would emit the band-truncated (clipped) alignment
# instead -- a different CIGAR/POS -- and the default-vs-full-width diff below
# would then FAIL. So these pairs are what makes this regression able to fail.
DEL=24 # > w0 (ADAPTIVE_BAND_START = 20)
L1=55
L2=55
n_del=0
for off in 500 1500 2500 3500; do
    d1="${seq:$off:$L1}${seq:$((off + L1 + DEL)):$L2}"
    d2fwd="${seq:$((off + 700)):120}"
    d2="$(printf '%s' "$d2fwd" | rev | tr 'ACGT' 'TGCA')"
    [ "${#d1}" -eq "$((L1 + L2))" ] && [ "${#d2}" -eq 120 ] || continue
    q1="$(printf 'I%.0s' $(seq 1 "$((L1 + L2))"))"
    q2="$(printf 'I%.0s' $(seq 1 120))"
    printf '@del%d/1\n%s\n+\n%s\n' "$n_del" "$d1" "$q1" >> "$r1"
    printf '@del%d/2\n%s\n+\n%s\n' "$n_del" "$d2" "$q2" >> "$r2"
    n_del=$((n_del + 1))
done
if [ "$n_del" -lt 2 ]; then
    echo "FAIL: built only $n_del deletion pairs (phix too short?) — test would be vacuous" >&2
    exit 1
fi
echo "fixture: $n_pairs exact PE pairs + $n_del deletion-bearing pairs (DEL=$DEL > w0) sliced from phix"

# --- Index the reference. ---
"$BWA_MEM3" index "$ref" > "$BAND_CERT_WORK_DIR/index.log" 2>&1

# --- The load-bearing A/B: default (band_cert on) vs --no-band-cert. ---
# @PG excluded from both sides: its CL: records the actual argv, which differs
# between the two invocations. Everything else must match byte-for-byte.
def_sam="$BAND_CERT_WORK_DIR/default.sam"
raw_sam="$BAND_CERT_WORK_DIR/no_band_cert.sam"
"$BWA_MEM3" mem "$ref" "$r1" "$r2" > "$def_sam" 2> /dev/null
"$BWA_MEM3" mem --no-band-cert "$ref" "$r1" "$r2" > "$raw_sam" 2> /dev/null

# Non-vacuity: both runs must actually emit alignment records.
if ! grep -q '^[^@]' "$def_sam" || ! grep -q '^[^@]' "$raw_sam"; then
    echo "FAIL: one of the runs emitted no alignment records (fixture too weak?)" >&2
    exit 1
fi

def_stripped="$BAND_CERT_WORK_DIR/default.stripped.sam"
raw_stripped="$BAND_CERT_WORK_DIR/no_band_cert.stripped.sam"
grep -v '^@PG' "$def_sam" > "$def_stripped"
grep -v '^@PG' "$raw_sam" > "$raw_stripped"
if ! diff -q "$def_stripped" "$raw_stripped" > /dev/null; then
    echo "FAIL: --no-band-cert SAM differs from the default (certified band is NOT byte-identical):" >&2
    diff "$def_stripped" "$raw_stripped" | head -20 >&2
    exit 1
fi
echo "PASS: SAM-text default (certified band) is byte-identical to --no-band-cert (@PG excluded)"

# Non-vacuity floor: the default output must actually contain a wide-band (> w0)
# deletion alignment, i.e. the fixture really does drive the extension past the
# narrow probe (where the certificate must decline). If no CIGAR carries a deletion
# larger than the probe width, the deletion pairs did not align as intended and the
# A/B above degenerates to the exact-match pairs only.
widest_del="$(mawk '!/^@/ {
    cig = $6
    while (match(cig, /[0-9]+D/)) {
        n = substr(cig, RSTART, RLENGTH - 1) + 0
        if (n > max) max = n
        cig = substr(cig, RSTART + RLENGTH)
    }
} END { print max + 0 }' "$def_sam")"
if [ "$widest_del" -le 20 ]; then
    echo "FAIL: default output has no deletion wider than w0=20 (widest=$widest_del) — the" >&2
    echo "      certificate is never exercised on a decline case, so this test is vacuous" >&2
    exit 1
fi
echo "  ok: default output exercises the certificate (widest deletion = ${widest_del} > w0=20)"

# --- Non-default parameter sweep: default vs --no-band-cert must match under
# scoring/gap/zdrop knobs too. The certificate bounds the optimal score but not the
# extension kernel's early-termination heuristics (zdrop / m==0 / band-edge shrink),
# so band_cert is applied only inside a parameter envelope and falls back to the
# exact full-width ladder outside it. Either way default == --no-band-cert: inside
# the envelope by the certificate, outside it because both run the full-width ladder.
ab_identical() { # <label> <mem args...>
    local label="$1"
    shift
    "$BWA_MEM3" mem "$@" "$ref" "$r1" "$r2" 2> /dev/null | grep -v '^@PG' > "$BAND_CERT_WORK_DIR/ab_def.sam"
    "$BWA_MEM3" mem --no-band-cert "$@" "$ref" "$r1" "$r2" 2> /dev/null | grep -v '^@PG' > "$BAND_CERT_WORK_DIR/ab_raw.sam"
    if ! diff -q "$BAND_CERT_WORK_DIR/ab_def.sam" "$BAND_CERT_WORK_DIR/ab_raw.sam" > /dev/null; then
        echo "FAIL: default != --no-band-cert under '$label' ($*)" >&2
        diff "$BAND_CERT_WORK_DIR/ab_def.sam" "$BAND_CERT_WORK_DIR/ab_raw.sam" | head -20 >&2
        exit 1
    fi
    echo "  ok: default == --no-band-cert under '$label'"
}
ab_identical "small zdrop (-d 24)" -d 24
ab_identical "large clip penalty (-L 8)" -L 8
ab_identical "custom gaps (-O 2 -E 1)" -O 2 -E 1
ab_identical "custom scoring (-A 2 -B 6)" -A 2 -B 6
ab_identical "combined (-d 15 -L 10)" -d 15 -L 10
echo "PASS: default == --no-band-cert across non-default -d/-L/-O/-E/-A/-B"

# --- BAM path: same assertion, if samtools is available. ---
if command -v samtools > /dev/null 2>&1; then
    def_bam="$BAND_CERT_WORK_DIR/default.bam"
    raw_bam="$BAND_CERT_WORK_DIR/no_band_cert.bam"
    "$BWA_MEM3" mem --bam "$ref" "$r1" "$r2" > "$def_bam" 2> /dev/null
    "$BWA_MEM3" mem --bam --no-band-cert "$ref" "$r1" "$r2" > "$raw_bam" 2> /dev/null
    def_rec="$BAND_CERT_WORK_DIR/default.bam.rec"
    raw_rec="$BAND_CERT_WORK_DIR/no_band_cert.bam.rec"
    samtools view --no-PG "$def_bam" > "$def_rec"
    samtools view --no-PG "$raw_bam" > "$raw_rec"
    if ! diff -q "$def_rec" "$raw_rec" > /dev/null; then
        echo "FAIL: --no-band-cert BAM records differ from the default:" >&2
        diff "$def_rec" "$raw_rec" | head -20 >&2
        exit 1
    fi
    echo "PASS: BAM records default (certified band) byte-identical to --no-band-cert"
else
    echo "SKIP: samtools not on PATH — BAM-path assertion skipped"
fi

# --- The flag must actually be recognized (not a fallback no-op). ---
# `mem` with no operands prints usage and exits nonzero, so capture to a file
# rather than pipe: under `set -o pipefail` the nonzero exit would mask grep.
usage_txt="$BAND_CERT_WORK_DIR/usage.txt"
"$BWA_MEM3" mem > "$usage_txt" 2>&1 || true
if ! grep -q -- '--no-band-cert' "$usage_txt"; then
    echo "FAIL: --no-band-cert is not listed in the usage text (help/flag drift)" >&2
    exit 1
fi
# --no-band-cert composes with --fast (band_cert is off under --fast anyway):
# this must not become a hard error.
if ! "$BWA_MEM3" mem --fast --no-band-cert "$ref" "$r1" "$r2" > /dev/null 2> "$BAND_CERT_WORK_DIR/compose.log"; then
    echo "FAIL: 'mem --fast --no-band-cert' exited nonzero, expected success:" >&2
    cat "$BAND_CERT_WORK_DIR/compose.log" >&2
    exit 1
fi
echo "PASS: --no-band-cert is advertised in usage and composes with --fast"

# --- Source guard: the certificate is ON by default, and the aggressive presets
# turn it OFF. This is not observable on phix's short reads (--adaptive-band is a
# no-op there), so pin the wiring structurally, the same way compat_byte_identical.sh
# guards its proper-pair derivation.
BWAMEM_SRC="$(dirname "$0")/../../src/bwamem.cpp"
FASTMAP_SRC="$(dirname "$0")/../../src/fastmap.cpp"
[[ -f "$BWAMEM_SRC" ]] || {
    echo "FAIL: cannot find src/bwamem.cpp at $BWAMEM_SRC" >&2
    exit 1
}
[[ -f "$FASTMAP_SRC" ]] || {
    echo "FAIL: cannot find src/fastmap.cpp at $FASTMAP_SRC" >&2
    exit 1
}

if ! grep -qE 'o(pt)?->band_cert *= *1' "$BWAMEM_SRC"; then
    echo "FAIL: bwamem.cpp no longer initializes band_cert = 1 (default-on regression)" >&2
    exit 1
fi
# --adaptive-band handler and the --fast preset block must both clear band_cert.
adaptive_off=$(grep -c 'OPT_ADAPTIVE_BAND.*band_cert *= *0' "$FASTMAP_SRC" || true)
if [[ "$adaptive_off" -lt 1 ]]; then
    echo "FAIL: --adaptive-band handler no longer sets band_cert = 0" >&2
    exit 1
fi
# Scope this to the --fast preset block specifically (anchored on its opt-out
# comment), so it cannot pass on the strength of the --adaptive-band or
# --no-band-cert assignments alone if --fast stops clearing band_cert.
if ! grep -qE 'band_cert *= *0;.*--fast opts out' "$FASTMAP_SRC"; then
    echo "FAIL: --fast preset no longer sets band_cert = 0 (aggressive-band opt-out lost)" >&2
    exit 1
fi
# --no-band-cert must route to band_cert = 0.
if ! grep -qE 'OPT_NO_BAND_CERT.*band_cert *= *0|band_cert *= *0;.*opt out of the certified' "$FASTMAP_SRC"; then
    echo "FAIL: OPT_NO_BAND_CERT no longer clears band_cert" >&2
    exit 1
fi
# The parameter-envelope gate must exist and be consulted, or non-default -d/-L/-A/-B
# would silently run the certified band outside its proven-safe regime.
if ! grep -q 'int mem_band_cert_params_safe(' "$BWAMEM_SRC"; then
    echo "FAIL: bwamem.cpp lost mem_band_cert_params_safe (certified-band safety envelope)" >&2
    exit 1
fi
if ! grep -q 'mem_band_cert_params_safe(opt)' "$FASTMAP_SRC"; then
    echo "FAIL: fastmap.cpp no longer gates band_cert on mem_band_cert_params_safe" >&2
    exit 1
fi
# The public extension entry point must re-apply that same envelope. main_mem()
# clears band_cert outside the safe regime for the CLI, but a library caller that
# mutates scoring/gap fields after mem_opt_init() and calls
# mem_chain2aln_across_reads_V2() directly would bypass that guard and select
# certified narrowing outside the envelope — changing alignment records. The
# entry point therefore sanitizes a local opt (band_cert && params_safe) so every
# band decision below falls back to the exact full-width ladder when unsafe. This
# is not observable on phix's short CLI reads (the CLI never reaches the entry
# point with an unsafe band_cert=1), so pin it structurally, like the gates above.
if ! grep -qE 'mem_band_cert_params_safe\(opt_in\)' "$BWAMEM_SRC"; then
    echo "FAIL: mem_chain2aln_across_reads_V2 no longer re-applies the band_cert safety" >&2
    echo "      envelope at the entry point (library-caller certified-narrowing regression)" >&2
    exit 1
fi
echo "PASS: band_cert defaults on; --fast/--adaptive-band/--no-band-cert clear it; envelope gate wired at the CLI and the public entry point (source guard)"

# --- Source guard: the tight_band (tb) bound must be derived from a REALIZABLE
# ungapped extension score. tb short-circuits the retry ladder (a band >= tb is
# certified complete), which is sound only if the offset-0 in-band run actually
# ACHIEVES the score S fed into the bound. ungapped_walk_score is the floored
# score under ksw_extend local-truncation semantics (once the running score hits
# 0 it stays 0) -- the value the rung-1 banded DP reaches on the diagonal, i.e. a
# valid lower bound on the in-band optimum. A NO-FLOOR score (one that lets a
# negative prefix recover via later matches) can EXCEED the floor-killed value
# the DP actually reaches; feeding that larger S shrinks tb below soundness and
# lets the ladder skip the wider rung a gapped alignment in (default_w, wider]
# genuinely needs -- a CIGAR/coordinate divergence from the full-width ladder.
# The trigger is an adversarial double-deletion geometry (a diag-0 floor-death
# whose out-of-band gapped optimum sits just past default_w); on single reads it
# is masked by mem_reg2aln recomputing the final CIGAR at its own band, so it
# cannot be pinned as a SAM byte-fixture. Guard the realizable-score choice
# structurally instead, the same way the wiring guards above do.
if ! grep -qE 'max_sc_proof = ungapped_walk_score\(' "$BWAMEM_SRC"; then
    echo "FAIL: tight_band no longer derives max_sc_proof from the realizable" >&2
    echo "      floored ungapped_walk_score -- the band proof requires an" >&2
    echo "      in-band-achievable score; a no-floor bound makes tb unsound" >&2
    exit 1
fi
if grep -qE 'ungapped_max_sc_from_bitmap' "$BWAMEM_SRC"; then
    echo "FAIL: the NO-FLOOR ungapped_max_sc_from_bitmap bound was reintroduced --" >&2
    echo "      it over-estimates the realizable extension optimum and makes the" >&2
    echo "      tight_band bound unsound (breaks byte-identity vs the full-width ladder)" >&2
    exit 1
fi
echo "PASS: tight_band bound uses the realizable floored ungapped score (source guard)"

echo "PASS: band_cert byte-identity regression"
