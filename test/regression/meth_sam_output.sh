#!/usr/bin/env bash
# test/regression/meth_sam_output.sh
#
# Regression: `--meth` selects bisulfite ALIGNMENT semantics, not an output
# format. Output format is chosen by `--bam` on every mode alike.
#
# It was not. `--meth` hard-set `opt->bam_mode = 1` in the option parser, and
# the writer opened `hts_open(path, "wb<level>")` unconditionally, so text SAM
# was unreachable under `--meth` — there was no flag that could undo it. That
# made the CLI mean two different things depending on the mode: `--bam` was
# opt-in without `--meth` and forced with it.
#
# This pins the uniform rule:
#   mem --meth        -> SAM text  (same default as non-meth)
#   mem --meth --bam  -> BGZF BAM  (same opt-in as non-meth)
#
# The third assertion is the one that matters: the two must carry byte-identical
# RECORDS. Text and BAM share one bam1_t construction path and differ only in
# the htsFile mode, so a divergence here means someone reintroduced a second
# formatting path for meth — exactly the drift that gave the three writers three
# different default @HD lines (fg-labs/bwa-mem3#288).
#
# `samtools` is a REQUIRED dependency, not an optional one: the record- and
# header-parity assertions below need it, and they are the only checks that
# would catch the drift described above. Skipping them and still exiting zero
# would report a green run for a check that never executed.
#
# Inputs (env vars):
#   BWA_MEM3             — path to the bwa-mem3 binary under test
#   METH_SAM_PHIX_FA     — path to test/fixtures/phix.fa
#   METH_SAM_WORK_DIR    — directory for fixture-private intermediates

set -euo pipefail

: "${BWA_MEM3:?BWA_MEM3 must be set}"
: "${METH_SAM_PHIX_FA:?METH_SAM_PHIX_FA must be set}"
: "${METH_SAM_WORK_DIR:?METH_SAM_WORK_DIR must be set}"

command -v samtools > /dev/null 2>&1 || {
    echo "FAIL: samtools not on PATH; the text/BAM record and header parity checks cannot run" >&2
    exit 1
}

mkdir -p "$METH_SAM_WORK_DIR"
cd "$METH_SAM_WORK_DIR"
# Remove only this script's own artifacts (never `rm -rf` a caller-supplied path).
rm -f phix.fa phix.fa.* phix.dict phix.meth.* r1.fq pe1.fq pe2.fq ./*.sam ./*.bam ./*.log ./*.err ./records-*.txt

cp "$METH_SAM_PHIX_FA" phix.fa
ref=phix.fa

"$BWA_MEM3" index "$ref" > index.log 2>&1
"$BWA_MEM3" index --meth "$ref" > index-meth.log 2>&1 \
    || {
        echo "FAIL: bwa-mem3 index --meth failed" >&2
        cat index-meth.log >&2
        exit 1
    }

# Bisulfite-converted reads sliced from phix. Deterministic, no PRNG.
#
# SE: unmethylated C -> T on the read, which is what the OT hypothesis undoes.
#
# PE: directional (Lister) protocol — convert the FRAGMENT C -> T, then R1 is
# its 5' end and R2 the reverse complement of its 3' end, so R2 shows G -> A in
# its own orientation (what OB expects). Converting both mates C -> T instead
# would make R2 unalignable under OB and every pair would fall out of pairing.
#
# The insert size MUST vary. mem_pair() scores against the inferred FR insert
# distribution, and a single fixed insert gives it zero width — mem_pair then
# returns 0 for every pair and they all take the `no_pairing` path, which
# silently skips the paired emission branch this test exists to cover. Verified:
# with a uniform 300 bp insert the branch is reached 0 times; with the spread
# below, 2468 times.
python3 - << 'PY'
seq = "".join(l.strip() for l in open("phix.fa") if not l.startswith(">")).upper()
comp = str.maketrans("ACGT", "TGCA")

with open("r1.fq", "w") as fh:
    for i, off in enumerate([200, 900, 1600, 2400]):
        read = seq[off:off+120].replace("C", "T")
        fh.write("@p%d\n%s\n+\n%s\n" % (i, read, "I"*120))

inserts = [260, 280, 300, 320, 340, 360, 380, 400]
with open("pe1.fq", "w") as a, open("pe2.fq", "w") as b:
    n = 0
    for off in range(0, len(seq) - 450, 2):
        ins = inserts[n % len(inserts)]
        conv = seq[off:off+ins].replace("C", "T")
        r1, r2 = conv[:120], conv[-120:][::-1].translate(comp)
        if len(r1) < 120 or len(r2) < 120:
            continue
        a.write("@v%d\n%s\n+\n%s\n" % (n, r1, "I"*120))
        b.write("@v%d\n%s\n+\n%s\n" % (n, r2, "I"*120))
        n += 1
PY

# --- default: text SAM, no --bam ---
"$BWA_MEM3" mem --meth "$ref" r1.fq > meth.sam 2> meth-sam.err \
    || {
        echo "FAIL: mem --meth exited non-zero ($?)" >&2
        cat meth-sam.err >&2
        exit 1
    }

# A BGZF/gzip stream starts with 0x1f 0x8b; SAM text starts with '@'. Read the
# first byte rather than trusting `file`, which is not installed everywhere.
first_byte=$(head -c 1 meth.sam | od -An -tx1 | tr -d ' \n')
if [ "$first_byte" != "40" ]; then
    echo "FAIL: mem --meth did not emit SAM text (first byte 0x$first_byte, expected 0x40 '@')" >&2
    exit 1
fi
grep -q '^@HD' meth.sam || {
    echo "FAIL: --meth SAM text has no @HD" >&2
    exit 1
}
grep -q '^@SQ' meth.sam || {
    echo "FAIL: --meth SAM text has no @SQ" >&2
    exit 1
}

# The meth-only tags must survive the text path — they are built into the
# bam1_t aux block, so a text writer that dropped them would be a real bug.
grep -q 'XM:Z:' meth.sam || {
    echo "FAIL: --meth SAM text is missing XM:Z" >&2
    exit 1
}
grep -q 'XG:Z:' meth.sam || {
    echo "FAIL: --meth SAM text is missing XG:Z" >&2
    exit 1
}
grep -q 'XR:Z:' meth.sam || {
    echo "FAIL: --meth SAM text is missing XR:Z" >&2
    exit 1
}

sam_records=$(grep -cv '^@' meth.sam || true)
if [ "$sam_records" -eq 0 ]; then
    echo "FAIL: --meth SAM text has no alignment records" >&2
    exit 1
fi

# --- PE: the paired emission branch must survive the text container ---
# This is the branch that actually broke. mem_aln2sam() short-circuits into
# bam1_t under --meth and never writes the SAM kstring_t, so the paired branch
# in mem_sam_pe_batch_post/mem_sam_pe has to be selected on "records are bam1_t"
# (bam_mode || meth_mode), not on bam_mode alone. Selecting on bam_mode alone
# sends --meth-without--bam into the text branch and aborts on its
# `assert(str.s != 0)`. Verified: that variant exits 134 on exactly this input.
"$BWA_MEM3" mem --meth "$ref" pe1.fq pe2.fq > meth-pe.sam 2> meth-pe.err \
    || {
        echo "FAIL: mem --meth PE exited non-zero ($?)" >&2
        tail -5 meth-pe.err >&2
        exit 1
    }
pe_records=$(grep -cv '^@' meth-pe.sam || true)
if [ "$pe_records" -eq 0 ]; then
    echo "FAIL: --meth PE SAM text has no alignment records" >&2
    exit 1
fi
# Proper pairs (0x2) prove pairing actually resolved rather than every pair
# falling out to the single-end path, which would not exercise the branch.
proper=$(mawk '!/^@/ && int($2/2)%2==1' meth-pe.sam | wc -l | tr -d ' ')
if [ "$proper" -eq 0 ]; then
    echo "FAIL: --meth PE produced no proper pairs — the paired branch was not exercised" >&2
    exit 1
fi

# --- opt-in: --meth --bam is BGZF ---
"$BWA_MEM3" mem --meth --bam "$ref" r1.fq > meth.bam 2> meth-bam.err \
    || {
        echo "FAIL: mem --meth --bam exited non-zero ($?)" >&2
        cat meth-bam.err >&2
        exit 1
    }
"$BWA_MEM3" mem --meth --bam "$ref" pe1.fq pe2.fq > meth-pe.bam 2> meth-pe-bam.err \
    || {
        echo "FAIL: mem --meth --bam PE exited non-zero ($?)" >&2
        tail -5 meth-pe-bam.err >&2
        exit 1
    }
magic=$(head -c 2 meth.bam | od -An -tx1 | tr -d ' \n')
if [ "$magic" != "1f8b" ]; then
    echo "FAIL: mem --meth --bam did not emit BGZF (magic 0x$magic, expected 0x1f8b)" >&2
    exit 1
fi

samtools quickcheck meth.bam meth-pe.bam

# --- the two containers must agree, record for record and header for header ---
# --no-PG so samtools does not stamp its own @PG onto one side only.

# Records: byte-identical, no exceptions. This is the assertion with teeth, and
# it runs on BOTH read layouts — SE and PE take different emission paths
# (mem_reg2sam vs the paired branch in mem_sam_pe_batch_post), and it was the
# paired one that regressed.
for layout in "se:meth.sam:meth.bam" "pe:meth-pe.sam:meth-pe.bam"; do
    name=${layout%%:*}
    rest=${layout#*:}
    txt=${rest%%:*}
    bin=${rest#*:}
    samtools view --no-PG "$txt" > "records-$name-text.txt"
    samtools view --no-PG "$bin" > "records-$name-bam.txt"
    if ! cmp -s "records-$name-text.txt" "records-$name-bam.txt"; then
        echo "FAIL: --meth SAM text and --meth --bam differ in $name records" >&2
        diff "records-$name-text.txt" "records-$name-bam.txt" | head -20 >&2 || true
        exit 1
    fi
done

# Header: identical except the `@PG ... CL:` line, which records the invoking
# command line and therefore MUST differ — one run was given --bam and the other
# was not. Strip only the CL: field, not the whole @PG record, so a drift in
# ID/PN/VN still fails.
strip_cl() { sed 's/\tCL:.*$//'; }
samtools view -H --no-PG meth.sam | strip_cl > text-hdr.sam
samtools view -H --no-PG meth.bam | strip_cl > bam-hdr.sam
if ! cmp -s text-hdr.sam bam-hdr.sam; then
    echo "FAIL: --meth SAM text and --meth --bam headers differ (beyond @PG CL:)" >&2
    diff text-hdr.sam bam-hdr.sam | head -20 >&2 || true
    exit 1
fi

echo "PASS: --meth defaults to SAM text; --meth --bam emits identical records as BGZF (SE $sam_records, PE $pe_records)"
