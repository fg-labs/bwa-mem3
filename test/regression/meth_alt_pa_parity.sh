#!/usr/bin/env bash
# test/regression/meth_alt_pa_parity.sh
#
# Regression: `--meth` must render the `pa:f:` tag from the same definition the
# SAM-text writer uses (fg-labs/bwa-mem3#365).
#
# alt_pa_parity.sh covers the default path — mem_aln2sam vs mem_aln_to_bam. It
# cannot reach `meth_mem_aln_to_bam`, which is a third, independent writer: it
# built `pa` from the raw unrounded quotient exactly as the `--bam` writer did,
# and would have kept doing so after the default path was fixed.
#
# Under `--meth` the text and BAM containers share one bam1_t construction path
# and differ only in the htsFile mode (see meth_sam_output.sh), so a SAM-vs-BAM
# comparison alone is nearly vacuous here — it only catches someone
# reintroducing a second meth formatting path. The load-bearing assertion is
# the third one: every emitted value must be the float32 of a three-decimal
# rendering, which is the contract `"%.3f"` defines and which the raw quotient
# violates. On the pre-fix binary this fixture emits values like 0.8866667
# where the fixed writer emits 0.887.
#
# Note that `--meth` SAM text is rendered by htslib from the bam1_t, not by
# bwa-mem3's own SAM writer, so htslib's `%g` prints the identical value as
# "0.96" where the default path's `"%.3f"` prints "0.960". Every comparison
# below is therefore on the DECODED float, never on the printed text.
#
# Why the fixture is ALT-aware and bisulfite-converted. `alt_sc` is only set on
# a non-ALT region whose better overlapping hit is on a contig `<prefix>.alt`
# marks, so no fixture without that sidecar emits `pa` at all — and the reads
# have to come from the ALT sequence for the primary hit to be the worse, non-
# ALT one that carries the tag. `--meth` reads `is_alt` from the ORIGINAL bns
# (fastmap.cpp), so the sidecar reaches the meth path unchanged. The reads are
# C->T converted, the OT hypothesis, so they align under bisulfite semantics.
#
# `samtools` is a REQUIRED dependency, not an optional one: decoding the BAM is
# the only way to read back the values being compared. Skipping and still
# exiting zero would report a green run for a check that never executed.
#
# Fixture is generated here, not committed (repo convention). Sequence comes
# from python3's Mersenne Twister at a fixed seed, so it is identical on every
# run and platform.
#
# Inputs (env vars):
#   BWA_MEM3              — path to the bwa-mem3 binary under test
#   METH_ALT_PA_WORK_DIR  — directory for fixture-private intermediates. This
#                           script writes and deletes files inside it, so give
#                           it a directory of its own.

set -euo pipefail

# `set -e` aborts on an unchecked failure but prints nothing, and this
# directory's contract (test/regression/README.md) is that a script reports
# every outcome with a `PASS:` or `FAIL:` line. The invocations below that carry
# their own `|| { echo "FAIL: ..."; }` block do not reach this trap -- a failure
# that is handled is not an ERR -- so it covers exactly the ones that are not
# individually wrapped: fixture generation, the plain index build, and the
# `samtools view` decode.
trap 'echo "FAIL: meth_alt_pa_parity: line $LINENO: $BASH_COMMAND" >&2' ERR

: "${BWA_MEM3:?BWA_MEM3 must be set}"
: "${METH_ALT_PA_WORK_DIR:?METH_ALT_PA_WORK_DIR must be set}"

command -v samtools > /dev/null 2>&1 || {
    echo "FAIL: samtools not on PATH; the --meth pa:f: comparison cannot run" >&2
    exit 1
}

mkdir -p "$METH_ALT_PA_WORK_DIR"
cd "$METH_ALT_PA_WORK_DIR"
# Remove only this script's own artifacts (never `rm -rf` a caller-supplied path).
rm -f ref.fa ref.fa.* reads.fq meth.sam meth.bam meth.bam.sam ./*.log ./*.err

python3 - << 'PY'
import random

# Fixed seed: the fixture must be byte-identical on every host and every run,
# so the PASS/FAIL below is about the writer and never about the input.
random.seed(365)
BASES = "ACGT"

primary = "".join(random.choice(BASES) for _ in range(20000))

# The ALT contig is chr1[5000:13000] with ~2% substitutions -- divergent enough
# that a read sampled from it scores strictly better against the ALT than
# against the primary, which is the precondition for alt_sc, and therefore for
# pa, being set on the primary record.
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

# Reads come from the ALT sequence, with every C converted to T: the directional
# (Lister) OT hypothesis, which is what --meth undoes. 150 bp is long enough
# that the ~2% divergence separates the two loci by a wide, deterministic
# margin even after conversion has collapsed one base.
with open("reads.fq", "w") as fh:
    for i in range(2000):
        start = random.randrange(0, len(alt) - 150)
        fh.write("@r%d\n%s\n+\n%s\n"
                 % (i, alt[start:start + 150].replace("C", "T"), "I" * 150))
PY

[ -s ref.fa ] && [ -s ref.fa.alt ] && [ -s reads.fq ] || {
    echo "FAIL: fixture generation produced empty files" >&2
    exit 1
}

"$BWA_MEM3" index ref.fa > index.log 2>&1
"$BWA_MEM3" index --meth ref.fa > index-meth.log 2>&1 || {
    echo "FAIL: bwa-mem3 index --meth failed" >&2
    cat index-meth.log >&2
    exit 1
}

"$BWA_MEM3" mem --meth -t 1 ref.fa reads.fq > meth.sam 2> meth.err || {
    echo "FAIL: mem --meth exited non-zero ($?)" >&2
    tail -5 meth.err >&2
    exit 1
}
"$BWA_MEM3" mem --meth --bam=0 -t 1 ref.fa reads.fq > meth.bam 2> meth-bam.err || {
    echo "FAIL: mem --meth --bam exited non-zero ($?)" >&2
    tail -5 meth-bam.err >&2
    exit 1
}
samtools view meth.bam > meth.bam.sam

# The fixture must actually emit pa, else every assertion below passes
# vacuously -- which is precisely how this defect survived until now.
pa_count=$(grep -c 'pa:f:' meth.sam || true)
if [ "$pa_count" -lt 100 ]; then
    echo "FAIL: --meth fixture emitted only $pa_count pa:f: records (expected >= 100);" >&2
    echo "      the .alt sidecar, the ALT divergence, or the C->T conversion is no" >&2
    echo "      longer doing its job" >&2
    exit 1
fi

python3 - << 'PY' || exit 1
import struct
import sys


def records(path):
    """[((QNAME, FLAG, RNAME, POS), pa float32 or None, tag names)]."""
    out = []
    with open(path) as fh:
        for line in fh:
            if line.startswith("@"):
                continue
            f = line.rstrip("\n").split("\t")
            pa = None
            names = []
            for tag in f[11:]:
                names.append(tag[:2])
                if tag.startswith("pa:f:"):
                    # Narrow through float32: a BAM 'f' field holds a float, so
                    # comparing the token as a double would report a difference
                    # no consumer can observe.
                    pa = struct.unpack("f", struct.pack("f", float(tag[5:])))[0]
            out.append(((f[0], f[1], f[2], f[3]), pa, names))
    return out


sam = records("meth.sam")
bam = records("meth.bam.sam")

if len(sam) != len(bam):
    print("FAIL: --meth SAM (%d) vs --meth --bam (%d) record count mismatch"
          % (len(sam), len(bam)), file=sys.stderr)
    sys.exit(1)

# --- 1. the two containers agree ---
# Near-vacuous by construction today (one bam1_t path, two htsFile modes), and
# deliberately kept: it is what fails if a second meth formatting path is ever
# reintroduced, which is the drift that produced three different writers here.
differing = []
pa_records = 0
for (skey, spa, _), (bkey, bpa, _) in zip(sam, bam):
    if skey != bkey:
        print("FAIL: records out of lockstep at %s vs %s" % (skey, bkey), file=sys.stderr)
        sys.exit(1)
    if spa is not None:
        pa_records += 1
    if spa != bpa:
        differing.append((skey[0], spa, bpa))

if differing:
    print("FAIL: %d of %d pa records differ between --meth SAM and --meth --bam"
          % (len(differing), pa_records), file=sys.stderr)
    for name, spa, bpa in differing[:5]:
        print("      %s: SAM %r vs BAM %r" % (name, spa, bpa), file=sys.stderr)
    sys.exit(1)

# --- 2. the rounding contract ---
# The assertion that actually pins meth_mem_aln_to_bam to the shared
# definition. "%.3f" is what the SAM-text writer emits and what both upstreams
# emit; the raw quotient this writer used to store does not survive it.
over_precision = [pa for _, pa, _ in bam
                  if pa is not None and abs(pa * 1000 - round(pa * 1000)) > 1e-4]
if over_precision:
    print("FAIL: %d --meth pa values carry more than three decimals (e.g. %r);"
          % (len(over_precision), over_precision[0]), file=sys.stderr)
    print("      meth_mem_aln_to_bam is not rendering pa through the shared"
          " SAM definition", file=sys.stderr)
    sys.exit(1)

# --- 3. secondary records carry no pa ---
secondary = [k for k, pa, _ in bam if int(k[1]) & 0x100 and pa is not None]
n_secondary = sum(1 for k, _, _ in bam if int(k[1]) & 0x100)
if not n_secondary:
    print("FAIL: fixture produced no secondary records; the pa-omission check"
          " would pass vacuously", file=sys.stderr)
    sys.exit(1)
if secondary:
    print("FAIL: %d secondary records carry pa (e.g. %s); mem_aln2sam emits it"
          " only for non-secondary records" % (len(secondary), secondary[0]),
          file=sys.stderr)
    sys.exit(1)

# --- 4. SA:Z precedes pa ---
both = [names for _, pa, names in bam if pa is not None and "SA" in names]
if not both:
    print("FAIL: no record carries both SA:Z and pa; the tag-order check would"
          " pass vacuously", file=sys.stderr)
    sys.exit(1)
out_of_order = [n for n in both if n.index("SA") > n.index("pa")]
if out_of_order:
    print("FAIL: %d records emit pa before SA:Z; mem_aln2sam emits SA:Z first"
          % len(out_of_order), file=sys.stderr)
    sys.exit(1)

print("PASS: --meth pa:f: matches the shared SAM definition (%d records, %d"
      " carrying pa, %d secondary, %d with SA:Z)"
      % (len(bam), pa_records, n_secondary, len(both)))
PY

echo "PASS: meth_alt_pa_parity"
