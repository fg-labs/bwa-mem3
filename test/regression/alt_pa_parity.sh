#!/usr/bin/env bash
# test/regression/alt_pa_parity.sh
#
# Regression: the `pa:f:` tag must carry the same value in SAM text and in
# `--bam` output of the same run (fg-labs/bwa-mem3#365).
#
# `pa` is the ratio of a hit's score to the score of the better overlapping
# hit on an ALT contig. The SAM writer renders it with "%.3f", as bwa and
# bwa-mem2 both do; the two BAM writers used to store the raw unrounded
# quotient, so a `--bam` consumer saw `pa:f:0.806723` where a `samtools view
# -b` of the same run's SAM text stored float32("0.807").
#
# Why this needs its own script. `pa` is unreachable without an ALT sidecar:
# `alt_sc` is only set for a non-ALT region whose better overlapping hit is on
# a contig `is_alt` marks, so no fixture that lacks a `<prefix>.alt` can emit
# the tag at all. bam_roundtrip.sh compares record counts, not tag values, and
# header_parity.sh does build an ALT-marked index but slices its reads from the
# PRIMARY contig, so nothing there ever gets a positive `alt_sc` either. That
# is why hundreds of millions of records of prior byte-identity testing never
# saw this defect, and why the coverage has to be an ALT-aware fixture whose
# reads come from the ALT sequence.
#
# The comparison is on the DECODED float, not on the text `samtools view`
# prints: samtools renders an 'f' aux field with %g, so the identical value
# 0.8 prints as "0.800" from bwa-mem3's SAM writer and "0.8" from a BAM
# decode. Comparing the printed strings would report ~30% spurious
# differences on a fixture where the two agree perfectly.
#
# Fixture is generated here, not committed (repo convention): a primary contig
# and an ALT contig that is a diverged copy of a window of it, plus SE reads
# sampled from the ALT sequence so the primary hit is the worse, non-ALT one
# that carries `pa`. Sequence comes from python3's Mersenne Twister at a fixed
# seed, so the fixture is identical on every run and platform.
#
# `samtools` is a REQUIRED dependency, not an optional one: decoding the BAM is
# the only way to read back the value this script exists to compare, so every
# assertion below depends on it. Skipping them and still exiting zero would
# report a green run for a check that never executed — the same reason
# default_hd_parity.sh and meth_sam_output.sh hard-fail rather than SKIP.
#
# Scope: the DEFAULT path only -- mem_aln2sam vs mem_aln_to_bam. The third
# writer, meth_mem_aln_to_bam (src/meth_bam.cpp), carries its own copy of the
# aux-tag block and is covered by meth_alt_pa_parity.sh next door, whose fixture
# additionally has to be bisulfite-converted to reach it.
#
# Inputs (env vars):
#   BWA_MEM3          — path to the bwa-mem3 binary under test
#   ALT_PA_WORK_DIR   — directory for fixture-private intermediates. This
#                       script writes and deletes files inside it, so give it
#                       a directory of its own.

set -euo pipefail

# `set -e` aborts on an unchecked failure but prints nothing, and this
# directory's contract (test/regression/README.md) is that a script reports
# every outcome with a `PASS:` or `FAIL:` line. The invocations below that carry
# their own `|| { echo "FAIL: ..."; }` block do not reach this trap -- a failure
# that is handled is not an ERR -- so it covers exactly the ones that are not
# individually wrapped: fixture generation and the `samtools view` decode.
trap 'echo "FAIL: alt_pa_parity: line $LINENO: $BASH_COMMAND" >&2' ERR

: "${BWA_MEM3:?BWA_MEM3 must be set}"
: "${ALT_PA_WORK_DIR:?ALT_PA_WORK_DIR must be set}"

command -v samtools > /dev/null 2>&1 || {
    echo "FAIL: samtools not on PATH; the SAM vs --bam pa:f: comparison cannot run" >&2
    exit 1
}

mkdir -p "$ALT_PA_WORK_DIR"
cd "$ALT_PA_WORK_DIR"
rm -f ref.fa ref.fa.* reads.fq out.sam out.bam out.bam.sam

python3 - << 'PY'
import random

# Fixed seed: the fixture must be byte-identical on every host and every run,
# so the PASS/FAIL below is about the writers and never about the input.
random.seed(365)
BASES = "ACGT"

primary = "".join(random.choice(BASES) for _ in range(20000))

# The ALT contig is chr1[5000:13000] with ~2% substitutions. Divergent enough
# that a read sampled from it scores strictly better against the ALT than
# against the primary -- which is the precondition for alt_sc, and therefore
# for pa, being set on the primary record.
window = list(primary[5000:13000])
for i, base in enumerate(window):
    if random.random() < 0.02:
        window[i] = random.choice([b for b in BASES if b != base])
alt = "".join(window)


def wrap(seq, width=60):
    return "\n".join(seq[i:i + width] for i in range(0, len(seq), width))


with open("ref.fa", "w") as fh:
    fh.write(">chr1\n%s\n>chr1_alt\n%s\n" % (wrap(primary), wrap(alt)))

# <prefix>.alt is a SAM file; bwa-mem3 reads only column 1 (QNAME) to set
# is_alt (bntseq.cpp). Keep the full record so the file is well-formed.
with open("ref.fa.alt", "w") as fh:
    fh.write("chr1_alt\t0\tchr1\t5001\t60\t8000M\t*\t0\t0\t*\t*\n")

# Reads come from the ALT sequence, so the hit on chr1 is the non-ALT one an
# ALT hit beats. 150 bp is long enough that the ~2% divergence separates the
# two loci by a wide, deterministic margin.
with open("reads.fq", "w") as fh:
    for i in range(2000):
        start = random.randrange(0, len(alt) - 150)
        fh.write("@r%d\n%s\n+\n%s\n" % (i, alt[start:start + 150], "I" * 150))
PY

[ -s ref.fa ] && [ -s ref.fa.alt ] && [ -s reads.fq ] || {
    echo "FAIL: fixture generation produced empty files" >&2
    exit 1
}

"$BWA_MEM3" index ref.fa > index.log 2>&1 || {
    echo "FAIL: bwa-mem3 index failed" >&2
    cat index.log >&2
    exit 1
}

"$BWA_MEM3" mem -t 1 ref.fa reads.fq > out.sam 2> sam.err || {
    echo "FAIL: mem exited non-zero ($?)" >&2
    tail -5 sam.err >&2
    exit 1
}
"$BWA_MEM3" mem -t 1 --bam=0 ref.fa reads.fq > out.bam 2> bam.err || {
    echo "FAIL: mem --bam exited non-zero ($?)" >&2
    tail -5 bam.err >&2
    exit 1
}
samtools view out.bam > out.bam.sam

python3 - << 'PY' || exit 1
import struct
import sys


def records(path):
    """[((QNAME, FLAG, RNAME, POS), pa-as-float32-or-None, core-fields)]."""
    out = []
    with open(path) as fh:
        for line in fh:
            if line.startswith("@"):
                continue
            f = line.rstrip("\n").split("\t")
            pa = None
            for tag in f[11:]:
                if tag.startswith("pa:f:"):
                    # Narrow through float32: a BAM 'f' field holds a float,
                    # so comparing the SAM token as a double would report a
                    # difference that no consumer can observe.
                    pa = struct.unpack("f", struct.pack("f", float(tag[5:])))[0]
                    break
            out.append(((f[0], f[1], f[2], f[3]), pa, f[:11]))
    return out


sam = records("out.sam")
bam = records("out.bam.sam")

if len(sam) != len(bam):
    print("FAIL: SAM (%d) vs --bam (%d) record count mismatch" % (len(sam), len(bam)),
          file=sys.stderr)
    sys.exit(1)

pa_records = 0
differing = []
core_mismatch = 0
for (skey, spa, score), (bkey, bpa, bcore) in zip(sam, bam):
    if skey != bkey:
        print("FAIL: records out of lockstep at %s vs %s" % (skey, bkey), file=sys.stderr)
        sys.exit(1)
    if score != bcore:
        core_mismatch += 1
    if spa is not None:
        pa_records += 1
    if spa != bpa:
        differing.append((skey[0], spa, bpa))

# The fixture must actually emit pa, else every assertion here passes vacuously
# -- which is precisely how this defect survived until now.
if pa_records < 100:
    print("FAIL: only %d records carry pa:f: (expected >= 100); the .alt sidecar"
          " or the ALT divergence is no longer doing its job" % pa_records,
          file=sys.stderr)
    sys.exit(1)

if core_mismatch:
    print("FAIL: %d records differ in the eleven core SAM fields" % core_mismatch,
          file=sys.stderr)
    sys.exit(1)

if differing:
    print("FAIL: %d of %d pa records differ between SAM and --bam" %
          (len(differing), pa_records), file=sys.stderr)
    for name, spa, bpa in differing[:5]:
        print("      %s: SAM %r vs BAM %r" % (name, spa, bpa), file=sys.stderr)
    sys.exit(1)

# The rounding itself, independent of the two containers agreeing: an unrounded
# quotient would agree across containers if BOTH writers regressed together.
over_precision = [pa for _, pa, _ in bam
                  if pa is not None and abs(pa * 1000 - round(pa * 1000)) > 1e-4]
if over_precision:
    print("FAIL: %d pa values carry more than three decimals (e.g. %r)" %
          (len(over_precision), over_precision[0]), file=sys.stderr)
    sys.exit(1)

print("PASS: pa:f: identical in SAM and --bam (%d records, %d carrying pa)" %
      (len(sam), pa_records))
PY

echo "PASS: alt_pa_parity"
