#!/usr/bin/env bash
set -euo pipefail
# test/regression/meth_seed_prune.sh
#
# Regression for --meth-seed-prune (the 3-letter over-seeding prune, default
# spec30 under --meth). Under --meth the bisulfite (3-letter) alphabet
# over-produces short, high-multiplicity spurious SMEMs; pruning them before SA
# resolution removes ~2/3 of meth seed volume for a large speedup at ~0
# truth-based accuracy cost. This pins the correctness invariants:
#
#  1. num_smem==1 must-fix. The prune's up-front compaction can drive a whole
#     batch's SMEM count to exactly 1 (a lone read, or a ragged trailing batch
#     whose one rid-group prunes to its single longest seed). A pre-existing
#     off-by-one in the resolve loop's guard then skipped the loop entirely and
#     dropped the sole surviving read as UNMAPPED. Forced here deterministically:
#     a single-read --meth alignment with T=1 (always regime B) and L=U=4096
#     (nothing but the always-kept longest SMEM survives) must STILL map.
#  2. Default is spec30 under --meth: a plain --meth run is byte-identical to an
#     explicit --meth-seed-prune=spec30 run.
#  3. High-MAPQ placement is invariant to the prune: --meth-seed-prune=off and
#     the spec30 default agree on POS/MAPQ/CIGAR for confidently-mapped reads
#     (the prune removes only spurious short/repetitive seeds, never the true
#     long/unique seed).
#  4. The pruned path is thread-count deterministic (-t1 == -t4), and the flag
#     and the BWAMEM3_METH_SEED_PRUNE env select the same mode (=off == env off).
#  5. The baseline rule runs end-to-end and places the confident reads.
#  6. Malformed values follow the cohort_ramp_validation contract: a bad FLAG
#     value is a hard error; a bad ENV value (mode or numeric threshold) warns
#     and falls back to the default instead of silently parsing a prefix.
#  7. The prune actually removes seeds (BWAMEM3_METH_SEED_PRUNE_STATS: spec30
#     keeps fewer than the input count -- a silent no-op would fail); regime A
#     (T above every SMEM length) keeps all and is byte-identical to off; and
#     baseline keeps strictly more than spec30, so the two rules genuinely differ.
#
# Inputs:
#   BWA_MEM3 — the bwa-mem3 binary under test: either a bare command name resolved
#              via PATH, or a path to the binary. A slash-containing relative path
#              (e.g. build/bwa-mem3) is resolved to an absolute path up front so it
#              survives the `cd "$WORK"` below; a bare command name is left as-is.
if [ -z "${BWA_MEM3:-}" ]; then
    echo "FAIL: BWA_MEM3 must be set" >&2
    exit 1
fi
# Normalize a slash-containing path to absolute before we change directory --
# otherwise a relative BWA_MEM3 would be unreachable after `cd "$WORK"`. A bare
# command name (no slash) is a PATH lookup and must stay untouched.
case "$BWA_MEM3" in
    */*) BWA_MEM3="$(cd "$(dirname "$BWA_MEM3")" && pwd)/$(basename "$BWA_MEM3")" ;;
esac
command -v samtools > /dev/null 2>&1 || {
    echo "SKIP: samtools not on PATH (--meth emits BAM)"
    exit 0
}
command -v mawk > /dev/null 2>&1 || {
    echo "SKIP: mawk not on PATH (required to inspect alignment records and stats)"
    exit 0
}

WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT
cd "$WORK"
fail() {
    echo "FAIL: $*" >&2
    exit 1
}
# Dump a BAM's records to a file, failing loudly if `samtools view` errors. A
# bare `cmp -s <(samtools view a) <(samtools view b)` swallows each view's exit
# status (only cmp's is seen), so two *failed* views yield equal empty streams
# and the comparison passes falsely. Materialize each side first, gated on the
# view succeeding, then cmp the files.
view_to() {
    samtools view "$1" > "$2" || fail "samtools view $1 failed"
}

# Deterministic 1500 bp reference (PRNG seed 4242), shared with meth_pe_placement.
REF=TCATTGGCTATCCTAACCCGACCCTAGGAGCGGTTGGCGTGTATGCCGTGAATTTTCTCATTTCCGCTAGACATAATCGTTCTGCCTATATCTGGACAACATCCCGGCGACTTAGGCGACCCACAGAATCGTCCCTTCTAACGTAGTTCGCATAGTTCCCGTCCGTAGCCGGACTATTCGAACACCCAGTATTCGATTAACTCGGGCTTGACGTATTAGAGGCGTTAGTGTGCCAGGTAAGATACGCCAACGGAATTAACCTCTGTGACACTCCGCGGAGCCTTCGGACATATAAGTGATCGGGTCTACGTTTGTTAGACTTGAGACGTCTGTTAAGAGTTGGGTCTAATAAATCGCCTACACGTGGAGTCTAACGGGGAAGCGTCGAATCCTGATACATCATATAATGGAGCGTGTTATGAAAAAAGAGCATTCCATTGTACGAGCCGTGCCAGAAACGGCTTGACTACGTGAGCGTAGTGTTAGATAAACAGGAAACACTGACGCGGTTAGAAGGCGGATTGCCGGTAGGTTTTGGAAACATAAATACACACGGTATCATGTTGGGTCACGATTCCTATCACCGCACAGGGCCAACCATAGAAGAACTGAAAGAACTAATCTGGCGGCGGGCTCGGTGCTTATATTTTCCACCCAACATCGTGCACATTAGGCTCACCGCGCCCTACGGGCGAAGGGTGCGTACGGTGTTTATAAGGCGTGACGGCCCCAAGTAGAGGGTAATTCTGTGAAAGAATCTCAGGACGGTGGCATGAATTCAATTCCTTTTAAACCTATCGTTCCGACCTTATGCAATCCTTCAATGAAGATCGTCAACGACCATCGTTCTTCTGCTTTAAGTGTGAGTTCTCTCTTACAAGCTAATACACCCCAGCGTTCTCCGTACTCTTCACTGCCCAAGCGAGGCTAACCTTTTGAAATGTCACAGTCGAAGCATATCTCCCGTACATCTTTTTCGGAGATCGCAGCTCGCGGAGCTATAAGCGACTTAAGCCCTTGTGTCGGTGATCCCAAGGGTCTGACTCCTGTACCAGGGTTACTGTTTCGCTTTACGGAGTAGCCTGTGAGGTGAACTGAAAGGAGCATATTTGAGATCTAAGATAGGGTCCTCCTCTGCGTCTACGTTCTCTCCGTTACGTACGGCTTCGCACCGGAGTGCATCTTGGCCCCGAAACGCACTGTGTGTGCTGATACAGCGTCCCTGGCCGGCCATGGGTTCAGAACTCCCGGGAACGCTTTTCAACTTAGAGGAACCCCGTCATGGAAGTAGATCGCGTCGAATGAGGGAGTTAGTCCTCGTTCCAGCTGGTAATTGTTTTACCGCTTGGGACCACTATAGGCCGCGGGTAGAAGTTGCTGGGTGTTGATTCCAACCCTCGAACCACGATACGACCTGCCATTTATGGCACAGTAAGGTTCAAACAGCATAATGAATACAGTTATAGTAACTTCCTCACGTACGATTAGGACGCAGCCTTG

# Fixtures that place uniquely at known loci (mapq 60), from meth_pe_placement.
FWD_R1=ATTCTGGCGATTTAGGTGACTCACAGAATCGTTCTTTCTAACGTAGTTTGTATAGTTCTC # OT fwd @101
FWD_R2=TAAGCGATTTATTAAACCCAACTCTTAACAAACGTCTCAAATCTAACAAACGTAAACCCG # OB rev @301
REV_R1=AGATTCTTTCACAGAATTACTCTCTATTTGGGGCTGTCACGCTTTATAAATATCGTACGC # OT rev @701
REV_R2=CTAACGCGATTAAAAGACAGATTGCCAGTAAGTTTTAGAAACATAAATACACACAGTATC # OB fwd @501

printf '>chrA\n%s\n' "$REF" > ref.fa
Q=$(printf 'I%.0s' $(seq 1 60))
emit() { printf '@%s\n%s\n+\n%s\n' "$1" "$2" "$3" > "$4"; }
"$BWA_MEM3" index --meth ref.fa > /dev/null 2>&1 || fail "index --meth nonzero exit"

# Primary record (not secondary 0x100, not supplementary 0x800) of a SE BAM,
# printed as "RNAME POS MAPQ FLAG CIGAR". Drains the stream (see meth_pe_placement).
prim() {
    samtools view "$1" | mawk '
      out=="" && (int($2/256)%2)==0 && (int($2/2048)%2)==0 { out = $3" "$4" "$5" "$2" "$6 }
      END { if (out != "") print out }'
}
# One mate (bit 64=READ1, 128=READ2) as "RNAME POS MAPQ CIGAR".
mate() {
    samtools view "$1" | mawk -v bit="$2" '
      out=="" && (int($2/bit)%2)==1 && (int($2/256)%2)==0 && (int($2/2048)%2)==0 { out = $3" "$4" "$5" "$6 }
      END { if (out != "") print out }'
}

# --- 1. num_smem==1 must-fix: single read prunes to one SMEM, still maps -------
# T=1 => M>=T always (regime B for every read); L=U=4096 => in regime B the only
# survivor is the always-kept longest SMEM. With a one-read FASTQ the whole batch
# therefore compacts to num_smem==1 -- the exact condition that used to skip the
# resolve loop and emit the read as unmapped.
emit one "$FWD_R1" "$Q" one.fq
env BWAMEM3_METH_SEED_PRUNE=spec30 BWAMEM3_METH_PRUNE_T=1 \
    BWAMEM3_METH_PRUNE_L=4096 BWAMEM3_METH_PRUNE_U=4096 \
    BWAMEM3_METH_SEED_PRUNE_STATS=1 \
    "$BWA_MEM3" mem --meth -t 1 ref.fa one.fq > one.bam 2> one.stats \
    || fail "single-read forced-prune mem --meth nonzero exit"
# Prove the setup actually reached num_smem==1 rather than just mapping: with
# L=U=4096 in regime B (T=1) only each read's always-kept longest SMEM survives,
# so this one-read batch must compact to kept_smem==1. Without this gate, a run
# that ignored PRUNE_L/PRUNE_U (retaining several SMEMs) would still pass the
# mapping check below and hide the num_smem==1 regression this test exists for.
kept="$(mawk '/^\[meth-seed-prune-stats\]/ {
                for (i = 1; i <= NF; i++) if ($i ~ /^kept_smem=/) { sub(/kept_smem=/, "", $i); print $i }
              }' one.stats)"
[ "$kept" = "1" ] \
    || fail "single-read forced-prune kept_smem='$kept', want 1 (PRUNE_L/PRUNE_U ignored, or the batch did not reach num_smem==1)"
samtools quickcheck one.bam || fail "single-read produced an invalid BAM"
one="$(prim one.bam)" || fail "reading one.bam failed"
[ -n "$one" ] || fail "no primary record for the single read"
read -r orn opos omapq oflag ocig <<< "$one"
[ $((oflag & 4)) -eq 0 ] \
    || fail "single-read --meth pruned to num_smem==1 came back UNMAPPED (must-fix #1 regressed): flag=$oflag"
[ "$orn" = "chrA" ] && [ "$opos" = "101" ] \
    || fail "single-read placement $orn:$opos, want chrA:101 (CIGAR $ocig, MAPQ $omapq)"

# --- 2/3. default == spec30, and off vs spec30 agree on confident placement ----
emit f "$FWD_R1" "$Q" f1.fq # PE mates must share a read name (mem_sam_pe_batch_post)
emit f "$FWD_R2" "$Q" f2.fq
"$BWA_MEM3" mem --meth -t 1 ref.fa f1.fq f2.fq > def.bam 2> /dev/null || fail "default --meth nonzero exit"
"$BWA_MEM3" mem --meth --meth-seed-prune=spec30 -t 1 ref.fa f1.fq f2.fq > spec.bam 2> /dev/null || fail "spec30 --meth nonzero exit"
"$BWA_MEM3" mem --meth --meth-seed-prune=off -t 1 ref.fa f1.fq f2.fq > off.bam 2> /dev/null || fail "off --meth nonzero exit"
for b in def.bam spec.bam off.bam; do samtools quickcheck "$b" || fail "$b invalid"; done

# 2: a plain --meth run defaults to spec30 -> byte-identical to explicit spec30.
view_to def.bam def.sam
view_to spec.bam spec.sam
cmp -s def.sam spec.sam \
    || fail "default --meth output != explicit --meth-seed-prune=spec30 (default is not spec30)"

# 3: the prune keeps the true seed, so confidently-mapped reads land identically
# with the prune off vs on (spec30). Compare each mate's POS/MAPQ/CIGAR.
for bit in 64 128; do
    a="$(mate off.bam "$bit")"
    [ -n "$a" ] || fail "off: no mate $bit"
    c="$(mate spec.bam "$bit")"
    [ -n "$c" ] || fail "spec30: no mate $bit"
    [ "$a" = "$c" ] || fail "mate $bit differs off vs spec30: '$a' vs '$c' (prune moved a confident read)"
    read -r _ _ mapq _ <<< "$c" # only MAPQ is asserted here
    [ "$mapq" = "60" ] || fail "mate $bit MAPQ $mapq, want 60 (fixture must map confidently for the invariance to bite)"
done

# --- 4. determinism (-t1 == -t4) and flag/env equivalence ----------------------
emit r "$REV_R1" "$Q" r1.fq
emit r "$REV_R2" "$Q" r2.fq
cat f1.fq f2.fq r1.fq r2.fq > many.fq # run single-end below, so shared names are fine
"$BWA_MEM3" mem --meth -t 1 ref.fa many.fq > t1.bam 2> /dev/null || fail "spec30 -t1 nonzero exit"
"$BWA_MEM3" mem --meth -t 4 ref.fa many.fq > t4.bam 2> /dev/null || fail "spec30 -t4 nonzero exit"
# Compare the record streams directly (no sort): output order is contractually
# deterministic across thread counts, so a -t1 vs -t4 record-order difference is
# itself the regression this check must catch -- sorting first would mask it.
view_to t1.bam t1.sam
view_to t4.bam t4.sam
cmp -s t1.sam t4.sam \
    || fail "-t1 != -t4 under spec30 (pruned path is thread-nondeterministic)"

env BWAMEM3_METH_SEED_PRUNE=off "$BWA_MEM3" mem --meth -t 1 ref.fa f1.fq f2.fq > envoff.bam 2> /dev/null \
    || fail "env off nonzero exit"
view_to off.bam off.sam
view_to envoff.bam envoff.sam
cmp -s off.sam envoff.sam \
    || fail "BWAMEM3_METH_SEED_PRUNE=off != --meth-seed-prune=off"

# --- 5. baseline mode runs end-to-end and places the confident reads -----------
# baseline is a shipped, user-selectable rule; run it as a real alignment (the
# footgun test only uses the word as an orphan token). For these unique reads it
# places identically to spec30.
"$BWA_MEM3" mem --meth --meth-seed-prune=baseline -t 1 ref.fa f1.fq f2.fq > base.bam 2> /dev/null \
    || fail "--meth-seed-prune=baseline nonzero exit"
samtools quickcheck base.bam || fail "baseline produced an invalid BAM"
for bit in 64 128; do
    b="$(mate base.bam "$bit")"
    [ -n "$b" ] || fail "baseline: no mate $bit"
    s="$(mate spec.bam "$bit")"
    [ "$b" = "$s" ] || fail "mate $bit differs baseline vs spec30 for a confident read: '$b' vs '$s'"
done

# --- 6. malformed flag/env values: flag errors, env warns and falls back -------
# Mirrors the cohort_ramp_validation contract: a bad FLAG value is a hard error
# (nothing is loaded yet), while a bad ENV value warns and falls back so a sweep
# knob typo does not throw away a loaded index. Covers the mode env and a numeric
# threshold env; the flag path already errors in the parser.
"$BWA_MEM3" mem --meth --meth-seed-prune=bogus -t 1 ref.fa f1.fq f2.fq > /dev/null 2> flagerr.txt \
    && fail "--meth-seed-prune=bogus exited 0; a malformed flag value must be a hard error"
grep -q "meth-seed-prune" flagerr.txt \
    || fail "--meth-seed-prune=bogus error does not name the flag: $(head -1 flagerr.txt)"

# Malformed MODE env: warn once, fall back to the flag/default (still maps).
env BWAMEM3_METH_SEED_PRUNE=bogus "$BWA_MEM3" mem --meth -t 1 ref.fa f1.fq f2.fq > modeenv.bam 2> modeenv.txt \
    || fail "malformed BWAMEM3_METH_SEED_PRUNE exited nonzero (must warn and fall back, not abort)"
grep -q "\[W::meth-seed-prune\]" modeenv.txt \
    || fail "malformed BWAMEM3_METH_SEED_PRUNE did not warn: $(head -1 modeenv.txt)"
# Fell back to the flag/default (spec30) -> byte-identical to the default run.
view_to modeenv.bam modeenv.sam
cmp -s modeenv.sam def.sam \
    || fail "malformed mode env did not fall back to the spec30 default"

# Malformed THRESHOLD env: warn and keep the default threshold (byte-identical).
env BWAMEM3_METH_PRUNE_T=bogus "$BWA_MEM3" mem --meth -t 1 ref.fa f1.fq f2.fq > threnv.bam 2> threnv.txt \
    || fail "malformed BWAMEM3_METH_PRUNE_T exited nonzero (must warn and keep default)"
grep -q "\[W::meth-seed-prune\].*BWAMEM3_METH_PRUNE_T" threnv.txt \
    || fail "malformed BWAMEM3_METH_PRUNE_T did not warn by name: $(head -1 threnv.txt)"
view_to threnv.bam threnv.sam
cmp -s threnv.sam def.sam \
    || fail "malformed threshold env did not keep the default threshold"

# --- 7. the prune actually removes seeds (not a silent no-op) ------------------
# The placement-invariance checks above pass even if the prune became a no-op:
# the spurious short seeds it drops score below -T and never surface in the SAM,
# so identical output does not prove the prune ran. BWAMEM3_METH_SEED_PRUNE_STATS
# reports input-vs-kept SMEM counts, which pins that it did.
# kept/total from the stats line, or "MISSING" if no line was emitted.
kept_of() { mawk '/meth-seed-prune-stats/ {for(i=1;i<=NF;i++){split($i,a,"=");v[a[1]]=a[2]}; print v["kept_smem"]"/"v["total_smem"]; f=1} END{if(!f)print "MISSING"}' "$1"; }
BWAMEM3_METH_SEED_PRUNE_STATS=1 "$BWA_MEM3" mem --meth --meth-seed-prune=spec30 -t 1 ref.fa f1.fq f2.fq > /dev/null 2> stat_spec.txt || fail "spec30 stats run nonzero exit"
spec_kt="$(kept_of stat_spec.txt)"
spec_kept="${spec_kt%/*}"
spec_tot="${spec_kt#*/}"
[ "$spec_kt" != "MISSING" ] || fail "no meth-seed-prune-stats line under spec30 (stats plumbing broken)"
[ "$spec_tot" -gt 0 ] || fail "spec30 stats report 0 total SMEMs (prune not reached)"
[ "$spec_kept" -lt "$spec_tot" ] || fail "spec30 kept all $spec_tot SMEMs -- the prune removed nothing (silent no-op?)"

# regime A: T above any SMEM length forces M<T for every read -> keep all. This
# exercises the keep-all branch and must be byte-identical to off.
BWAMEM3_METH_SEED_PRUNE_STATS=1 BWAMEM3_METH_PRUNE_T=4096 "$BWA_MEM3" mem --meth --meth-seed-prune=spec30 -t 1 ref.fa f1.fq f2.fq > rega.bam 2> stat_rega.txt || fail "regime-A run nonzero exit"
rega_kt="$(kept_of stat_rega.txt)"
# Guard against a missing stats line first: with no '/' in "MISSING" both
# expansions below yield "MISSING", so the keep-all comparison would pass
# vacuously if the stats plumbing broke for this run (cf. spec_kt/base_kt).
[ "$rega_kt" != "MISSING" ] || fail "no meth-seed-prune-stats line in regime A"
[ "${rega_kt%/*}" = "${rega_kt#*/}" ] || fail "regime A (T=4096) did not keep all SMEMs: $rega_kt"
view_to rega.bam rega.sam
cmp -s rega.sam off.sam || fail "regime-A keep-all is not byte-identical to off"

# baseline and spec30 are genuinely different rules: on this fixture baseline
# keeps strictly more (it never drops a unique seed; spec30's regime B does).
BWAMEM3_METH_SEED_PRUNE_STATS=1 "$BWA_MEM3" mem --meth --meth-seed-prune=baseline -t 1 ref.fa f1.fq f2.fq > /dev/null 2> stat_base.txt || fail "baseline stats run nonzero exit"
base_kt="$(kept_of stat_base.txt)"
[ "$base_kt" != "MISSING" ] || fail "no stats line under baseline"
[ "${base_kt%/*}" -gt "$spec_kept" ] || fail "baseline kept ${base_kt%/*} <= spec30 kept $spec_kept -- the two rules did not diverge on this fixture"

echo "PASS: meth_seed_prune (num_smem==1 maps; default==spec30; off/spec30 placement-invariant; -t1==-t4; flag==env; baseline maps; malformed flag errors, env warns+falls back; prune removes seeds; regime-A keep-all==off; baseline!=spec30)"
