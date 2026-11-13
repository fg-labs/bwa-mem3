#!/usr/bin/env bash
# test/regression/alt_project.sh
#
# End-to-end behavior of --alt-project.
#
# ALT-aware mapping keeps a primary-assembly hit primary whatever an ALT contig
# scores, so a read whose haplotype exists only as an ALT is placed by the
# primary assembly alone. --alt-project lifts ALT hits onto the primary assembly
# through <idxbase>.alt before pairing and MAPQ, so the ALT decides.
#
# Fixture (generated): copy X on chrA and a paralog Y on chrB differing at 8
# sites. chrA_alt1 is a haplotype of X that carries Y's base at 6 of them (gene
# conversion), one private SNP and a 3-bp deletion. Read pairs come from
# chrA_alt1. Reads covering the private SNP match the ALT best; on the primary
# assembly alone they cannot tell X from Y.
#
# Asserts:
#   * without the flag, two runs give identical output (determinism; flag-off
#     byte-identity to earlier releases is pinned by the rest of the suite), and
#     the SNP reads are not all at chrA with MAPQ > 0 (so the test is not vacuous)
#   * with --alt-project (paired, -I fixed): every SNP-covering read is on chrA
#     at its lifted position (chrA:3001 + its ALT start, from the read name) as
#     150M with MAPQ > 0, both mates of its pair are on chrA as a proper pair, its
#     AS is the ALT score that placed it (150), and its NM is exactly the number
#     of bases where its ALT window differs from chrA; no record (RNAME, RNEXT or
#     any tag) names the ALT contig; and every read spanning the ALT deletion
#     carries the 3D at the lifted position (the read name encodes its ALT start,
#     an oracle independent of the aligner)
#   * a reverse-complemented ALT (.alt FLAG 16) gives the same records
#   * ALTs whose insertion and deletion cancel inside a read (1I...1D, 2I...2D,
#     5I...5D, 40 bp apart; 35I...35D, 30 bp apart, the widest excursion the
#     re-alignment band of a 150 bp read can follow) and one with a wide pair
#     (25I then 30D): every lifted read is written at its lifted position with
#     both indels and NM equal to the indel bases, not as an ungapped run of
#     mismatches; for the 1-5 bp balanced pairs its record equals the default
#     run's chrA record apart from AS and the tags that name the ALT hit
#   * a 40I...40D ALT, past that band: its reads stay on the ALT, with records
#     equal to the default run's, not lifted with a garbled CIGAR
#   * a read that runs off the end of a divergent ALT onto sequence matching
#     chrA: its shorter ALT hit outscores its full chrA hit, but not by more than
#     the 50-54 chrA-aligned bases it leaves out are worth, so the full chrA record
#     is kept, unclipped: equal to the default run's primary record apart from
#     the SA and pa tags, with no ALT supplementary record (the default writes
#     one, so the check is not vacuous)
#   * the other direction: a read whose ALT hit misses 4 bases at its start
#     (4S146M, AS 146) but beats its full chrA hit (AS 120) by more than they are
#     worth is placed there by the lifted hit, not at a paralog elsewhere on chrA
#     (AS 130) where the default places it
#   * a tandem repeat U U U' on chrA whose ALT is U U' U': a read across copies 1
#     and 2 matches the ALT at its own locus and chrA one unit on equally well, so
#     both placements survive, as MAPQ 0 with XS equal to AS and, under -a, a
#     secondary record one unit away
#   * a fixture whose ALT starts with 300 bp of ALT-only sequence (.alt CIGAR
#     300S1000M): reads inside that clip stay on the ALT contig; reads straddling
#     its edge keep their ALT hit on the ALT (records equal to the default run's,
#     no long insertion); reads past it are at chrA or chrB at the lifted position
#     and no record of theirs names the ALT, which the default run's do (so the
#     check is not vacuous); an unmapped read fails the check
#   * records are byte-identical (apart from @PG) at -t 1, 2 and 4; --fast at
#     -t 1 and 4, where it also lifts the SNP reads to their chrA positions and
#     no record names the ALT; and --dedup-reads on and off over input with every
#     pair duplicated, with the memo shown to engage (dup > 0, final=ON) and the
#     BWAMEM3_DEDUP_READS_VERIFY check of duplicate regions passing
#   * matching ALT haplotypes represented at both paralogs preserve low paired
#     confidence for identical ambiguous reads, while a mate in distinct flanks
#     resolves either locus; an asymmetric-reference negative control fails the
#     uncertainty check, and the paired-confidence experiment keeps placement
#   * --alt-project is refused with -j, with --meth, without a .alt and with a
#     .alt that has no usable line, and --compat-allow-divergent admits it under
#     --compat (the --compat refusal itself is checked in compat_byte_identical.sh)
#   * --alt-project-exclude-bed: help lists it; a BED that excludes a .alt line
#     whose later duplicate line is usable is accepted (plain, gzip and with
#     --alt-project-min-gap 0, records identical apart from @PG; and under
#     --compat-allow-divergent) and logs exact counts, while plain --alt-project
#     logs no BED; it is refused, with a message naming the cause and no
#     alignment record, without --alt-project, when repeated, with a missing,
#     empty or stdin argument, with -j, --meth, --compat, a positive gap or
#     BWAMEM3_ALT_PROJECT_PAIR_GAP=1, for a missing, malformed or truncated gzip
#     BED, a BED naming an unknown contig or an ALT contig or reaching past a
#     contig's end (each message carries the option, the BED path and the line),
#     a BED with no nonempty interval, one that excludes nothing, one that
#     excludes every usable line, and a .alt with no usable line (keeping the
#     ordinary error)
#   * mixed exclusion (two independent ALT haplotypes, one excluded by a BED
#     interval inside its primary footprint but away from its informative reads,
#     one retained): every informative end from the retained ALT is lifted to its
#     chrA position with XL:Z, and none from the excluded ALT carries its XL:Z
#     while its records still name that ALT; the check fails on unrestricted
#     --alt-project (for the excluded reads) and without projection (for the
#     retained ones). The runtime BED run equals the run with a FLAG-masked copy
#     of the .alt (decided by the script, beside the same index), every header
#     line, record, field and tag in order apart from @PG CL, single-end and
#     paired at -t 1, 2 and 4, --fast at -t 1 and 4, and --dedup-reads on and off
#     over duplicated pairs, with -K and -I fixed; with projection off the masked
#     .alt maps as the original; and direct BAM decodes to the runtime SAM and to
#     the masked run's BAM (OPEN without samtools, a failure with
#     ALT_PROJECT_REQUIRE_SAMTOOLS=1)
#
# Inputs:
#   BWA_MEM3 — path to the bwa-mem3 binary under test
#   ALT_PROJECT_REQUIRE_SAMTOOLS (optional) — 1 makes a missing samtools a
#     failure instead of an OPEN mixed-exclusion BAM check
set -euo pipefail

if [[ -z ${BWA_MEM3:-} ]]; then
    echo "FAIL: BWA_MEM3 must be set" >&2
    exit 1
fi
command -v python3 > /dev/null 2>&1 || {
    echo "FAIL: python3 not on PATH (it generates the fixture)" >&2
    exit 1
}
case "$BWA_MEM3" in
    /*) BIN="$BWA_MEM3" ;;
    *) BIN="$PWD/$BWA_MEM3" ;;
esac
[[ -x "$BIN" ]] || {
    echo "FAIL: BWA_MEM3 ($BIN) is not executable" >&2
    exit 1
}

WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT
cd "$WORK"
fail() {
    echo "FAIL: $*" >&2
    exit 1
}

python3 - << 'EOF'
import random
random.seed(7)
B = "ACGT"
def rnd(n): return "".join(random.choice(B) for _ in range(n))
def mut(c): return random.choice([x for x in B if x != c])
RC = str.maketrans("ACGT", "TGCA")
def rc(s): return s.translate(RC)[::-1]
X = rnd(2000)
chrA = rnd(3000) + X + rnd(3000)
psv = sorted(random.sample(range(200, 1800), 8))
Y = list(X)
for p in psv: Y[p] = mut(X[p])
chrB = rnd(3000) + "".join(Y) + rnd(3000)
alt = list(X)
for p in psv[:6]: alt[p] = Y[p]
alt[1000] = mut(X[1000])                        # private SNP
# ALT positions (all before the deletion, so equal to X positions) where the ALT
# differs from chrA: the oracle for a lifted read's NM.
open("diffs.txt", "w").write(" ".join(str(p) for p in sorted(psv[:6] + [1000])) + "\n")
alt = "".join(alt[:1500] + alt[1503:])          # 3-bp deletion at X:1500
cig = "1500M3D497M"
for name, seq, flag in (("fwd", alt, 0), ("rev", rc(alt), 16)):
    open(f"{name}.fa", "w").write(f">chrA\n{chrA}\n>chrB\n{chrB}\n>chrA_alt1\n{seq}\n")
    open(f"{name}.fa.alt", "w").write(f"chrA_alt1\t{flag}\tchrA\t3001\t60\t{cig}\t*\t0\t0\t*\t*\n")
open("noalt.fa", "w").write(f">chrA\n{chrA}\n>chrB\n{chrB}\n>chrA_alt1\n{alt}\n")
with open("r1.fq", "w") as r1, open("r2.fq", "w") as r2:
    for i in range(60):
        s = random.randint(0, len(alt) - 400)
        f = alt[s:s + 400]
        q = "I" * 150
        r1.write(f"@p{i}_{s}\n{f[:150]}\n+\n{q}\n")
        r2.write(f"@p{i}_{s}\n{rc(f[-150:])}\n+\n{q}\n")
# Clip fixture: 300 bp of ALT-only sequence, then X[0:1000] at chrA:3001.
alt2 = rnd(300) + X[:1000]
open("clip.fa", "w").write(f">chrA\n{chrA}\n>chrB\n{chrB}\n>chrA_alt2\n{alt2}\n")
open("clip.fa.alt", "w").write("chrA_alt2\t0\tchrA\t3001\t60\t300S1000M\t*\t0\t0\t*\t*\n")
with open("clip.fq", "w") as fq:
    for s in list(range(0, 300, 10)) + list(range(400, 1100, 50)):
        fq.write(f"@c_{s}\n{alt2[s:s + 150]}\n+\n{'I' * 150}\n")
# Indel fixtures: an ALT of Z = rnd(2000) at chrA:3001 with an insertion after
# Z[1000) and a deletion GAP bases later. Single-end reads from the ALT start at
# s < 1000 and cover both indels; their lifted position is chrA:3001 + s.
def indel_fixture(name, ins, gap, dl, seed):
    rng = random.Random(seed)
    r = lambda n: "".join(rng.choice(B) for _ in range(n))
    Z = r(2000)
    a = Z[:1000] + r(ins) + Z[1000:1000 + gap] + Z[1000 + gap + dl:]
    cig = f"1000M{ins}I{gap}M{dl}D{2000 - 1000 - gap - dl}M"
    open(f"{name}.fa", "w").write(f">chrA\n{r(3000) + Z + r(3000)}\n>chrA_alt1\n{a}\n")
    open(f"{name}.fa.alt", "w").write(f"chrA_alt1\t0\tchrA\t3001\t60\t{cig}\t*\t0\t0\t*\t*\n")
    with open(f"{name}.fq", "w") as fq:
        for s in (930, 950, 970):
            fq.write(f"@i_{s}_{ins}_{dl}\n{a[s:s + 150]}\n+\n{'I' * 150}\n")
for k in (1, 2, 5):
    indel_fixture(f"bal{k}", k, 40, k, 100 + k)
indel_fixture("bal35", 35, 30, 35, 135)
indel_fixture("bal40", 40, 30, 40, 140)
indel_fixture("wide", 25, 30, 30, 12)
# End fixture: an ALT of W[0:1000) at chrA:3001 with 12 ALT-only SNPs in its last
# 100 bases. Each read is the ALT's last 100 - k bases then W's next 50 + k, which
# match chrA exactly: 100 - k on the ALT (no mismatch) against 150 on chrA (12
# mismatches).
# Its own generator, so the fixtures after it are as before.
rng = random.Random(11)
rw = lambda n: "".join(rng.choice(B) for _ in range(n))
W = rw(2000)
altw = list(W[:1000])
for p in range(904, 1000, 8): altw[p] = rng.choice([x for x in B if x != altw[p]])
altw = "".join(altw)
open("end.fa", "w").write(f">chrA\n{rw(3000) + W + rw(3000)}\n>chrA_alt1\n{altw}\n")
open("end.fa.alt", "w").write("chrA_alt1\t0\tchrA\t3001\t60\t1000M\t*\t0\t0\t*\t*\n")
with open("end.fq", "w") as fq:
    for k in (0, 2, 4):
        seq = altw[900 + k:1000] + W[1000:1050 + k]
        fq.write(f"@e_{900 + k}\n{seq}\n+\n{'I' * len(seq)}\n")
# Duplicated pairs (new names) for the --dedup-reads check.
with open("r1d.fq", "w") as o1, open("r2d.fq", "w") as o2:
    for src, dst in (("r1.fq", o1), ("r2.fq", o2)):
        lines = open(src).read().splitlines()
        for i in range(0, len(lines), 4):
            dst.write("\n".join(lines[i:i + 4]) + "\n")
            dst.write("\n".join([lines[i] + "_dup"] + lines[i + 1:i + 4]) + "\n")
# Tandem fixture: U U U' on chrA, U U' U' on the ALT (U' = U with one SNP).
U = rnd(100)
U2 = U[:20] + mut(U[20]) + U[21:]
lf, rt = rnd(500), rnd(500)
open("tan.fa", "w").write(f">chrA\n{rnd(3000) + lf + U + U + U2 + rt + rnd(3000)}\n>chrA_alt1\n{lf + U + U2 + U2 + rt}\n")
open("tan.fa.alt", "w").write("chrA_alt1\t0\tchrA\t3001\t60\t1300M\t*\t0\t0\t*\t*\n")
alt3 = lf + U + U2 + U2 + rt
with open("tan.fq", "w") as fq:
    for s in (500, 510, 530, 550):  # across copies 1 and 2, over U2's SNP
        fq.write(f"@t_{s}\n{alt3[s:s + 150]}\n+\n{'I' * 150}\n")
# Paralog fixture (its own generator): an ALT of V[0:1000) at chrA:3001 with
# SNPs at s+1, s+3 and six in s+104..s+144, s = 400. The read is V[s:s+100] then
# the ALT's next 50 bases: on the ALT it aligns from base 4 (4S146M, AS 146), on
# chrA's copy of V over all 150 bases with the six SNPs (AS 120), and on a
# paralog P (the read with 4 mismatches, at chrA:8001) at AS 130.
rng = random.Random(5)
rv = lambda n: "".join(rng.choice(B) for _ in range(n))
mutv = lambda c: rng.choice([x for x in B if x != c])
V = rv(2000)
s = 400
altv = list(V[:1000])
for p in [s + 1, s + 3] + list(range(s + 104, s + 150, 8)): altv[p] = mutv(altv[p])
altv = "".join(altv)
read = V[s:s + 100] + altv[s + 100:s + 150]
P = list(read)
for p in (20, 60, 90, 130): P[p] = mutv(P[p])
open("par.fa", "w").write(f">chrA\n{rv(3000) + V + rv(3000) + ''.join(P) + rv(1000)}\n>chrA_alt1\n{altv}\n")
open("par.fa.alt", "w").write("chrA_alt1\t0\tchrA\t3001\t60\t1000M\t*\t0\t0\t*\t*\n")
open("par.fq", "w").write(f"@v_{s}\n{read}\n+\n{'I' * 150}\n")
# The reader every check below shares: a SAM's records in file order, grouped by
# QNAME, the primary record of each read, and a record's tags.
open("samrecs.py", "w").write('''def records(path):
    for l in open(path):
        if l[0] != "@": yield l.rstrip("\\n").split("\\t")
def recs(path):
    out = {}
    for t in records(path): out.setdefault(t[0], []).append(t)
    return out
def prim(path):
    return {t[0]: t for t in records(path) if not int(t[1]) & 0x900}
def tags(t):
    return {x[:2]: x[5:] for x in t[11:]}
''')
EOF

for ref in fwd rev noalt clip bal1 bal2 bal5 bal35 bal40 wide tan end par; do
    "$BIN" index "$ref.fa" > "idx.$ref.log" 2>&1 || fail "index $ref failed"
done
# An index whose .alt has no usable line (its CIGAR length does not match the ALT).
for f in fwd.fa.*; do
    case "$f" in
        *.alt) ;;
        *) ln -s "$f" "bad.fa${f#fwd.fa}" ;;
    esac
done
ln -s fwd.fa bad.fa
printf 'chrA_alt1\t0\tchrA\t3001\t60\t100M\t*\t0\t0\t*\t*\n' > bad.fa.alt

aln() { # ref out extra...
    local ref=$1 out=$2
    shift 2
    "$BIN" mem -t 2 -I 400,10 "$@" "$ref.fa" r1.fq r2.fq > "$out" 2> "$out.log" || fail "mem $* on $ref failed (see $out.log)"
}
# Summarize primary records: name, mate, chrom, pos, MAPQ, CIGAR, FLAG, AS, NM,
# and the mate's ALT span (R1 = [s, s+150), R2 = [s+250, s+400)).
summ() {
    python3 - "$1" << 'EOF'
import sys
from samrecs import records, tags
for t in records(sys.argv[1]):
    flag = int(t[1])
    if flag & 0x900: continue
    s = int(t[0].split("_")[1]); m = 1 if flag & 0x40 else 2
    a, b = (s, s + 150) if m == 1 else (s + 250, s + 400)
    tg = tags(t)
    print(t[0], m, t[2], t[3], t[4], t[5], flag, tg.get("AS", "NA"), tg.get("NM", "NA"), a, b)
EOF
}

aln fwd def1.sam
aln fwd def2.sam
grep -v '^@PG' def1.sam > def1.nopg
grep -v '^@PG' def2.sam > def2.nopg
cmp -s def1.nopg def2.nopg || fail "default output is not deterministic"
grep -q "alt-project" def1.sam.log && fail "default run reports alt-project"

aln fwd proj.sam --alt-project
grep -q "alt-project: 1 ALT contigs" proj.sam.log || fail "--alt-project did not report the lift table"
aln rev projr.sam --alt-project

# Trace only selected complete read groups, without changing any SAM records.
trace_name=$(
    python3 - << 'EOF'
with open("r1.fq") as fq:
    print(fq.readline().strip()[1:])
EOF
)
BWAMEM3_ALT_PROJECT_TRACE="$trace_name" aln fwd trace.sam --alt-project
python3 - "$trace_name" << 'EOF'
import sys
from samrecs import records
name = sys.argv[1]
assert list(records("trace.sam")) == list(records("proj.sam")), "trace changed SAM records"
rows = [l.rstrip().split("\t") for l in open("trace.sam.log") if l.startswith("[alt-project-trace]\t")]
assert rows, "selected read has no projection trace"
assert {r[1] for r in rows} == {name}, "unselected read leaked into trace"
assert {r[3] for r in rows} >= {"input", "project", "compare", "gate", "group", "summary", "paired_input", "paired_result"}, "missing candidate stages"
assert len({r[2] for r in rows if r[3] == "summary"}) == 2, "trace lost a mate"
assert any("gap=0" in r for r in rows if r[3] == "summary"), "disabled gate was not identified"
for mate in (1, 2):
    with open(f"r{mate}.fq") as fq:
        record = "".join(next(fq) for _ in range(4))
    with open(f"trace{mate}.fq", "w") as fq:
        fq.write(record * 40)
EOF
BWAMEM3_ALT_PROJECT_TRACE="absent,$trace_name" "$BIN" mem -t 4 -I 400,10 --dedup-reads off --alt-project fwd.fa trace1.fq trace2.fq > trace.cap.sam 2> trace.cap.log
python3 - << 'EOF'
rows = [l.rstrip().split("\t") for l in open("trace.cap.log") if l.startswith("[alt-project-trace]\t")]
assert sum(r[3] == "begin" for r in rows) == 64, "read invocation cap failed"
assert sum(r[3] == "summary" for r in rows) == 64, "bounded traces lost summaries"
assert sum(r[3] == "call_limit" and "phase=0" in r for r in rows) == 1, "invocation truncation must be explicit"
EOF

summ def1.sam > def.tsv
summ proj.sam > proj.tsv
summ projr.sam > projr.tsv

# SNP-covering mates (ALT coordinate 1000).
n_snp=$(awk '$10 <= 1000 && 1000 < $11' proj.tsv | wc -l)
[[ $n_snp -ge 5 ]] || fail "fixture has only $((n_snp)) SNP-covering reads"
bad=$(awk '$10 <= 1000 && 1000 < $11 && ($3 != "chrA" || $5 == 0)' proj.tsv | wc -l)
[[ $bad -eq 0 ]] || fail "--alt-project left $((bad)) SNP-covering reads off chrA or at MAPQ 0"
def_ok=$(awk '$10 <= 1000 && 1000 < $11 && $3 == "chrA" && $5 > 0' def.tsv | wc -l)
[[ $def_ok -lt $n_snp ]] || fail "default already places every SNP read at chrA with MAPQ > 0 (test is vacuous)"
# AS is the ALT score that placed the read; the read is 150M at its lifted
# position (chrA:3001 + its ALT start), and NM is exactly the number of bases
# where its ALT window differs from chrA.
cat > snpcheck.py << 'EOF'
import sys
diffs = [int(x) for x in open(sys.argv[2]).read().split()]
n = 0
for l in open(sys.argv[1]):
    name, m, c, pos, q, cig, flag, AS, NM, a, b = l.split()
    a, b = int(a), int(b)
    if not (a <= 1000 < b): continue
    n += 1
    want_nm = sum(a <= d < b for d in diffs)
    if c != "chrA" or int(pos) != 3001 + a or cig != "150M" or AS != "150" or int(NM) != want_nm:
        print(f"{name}/{m}: {c}:{pos} {cig} AS:{AS} NM:{NM}, want chrA:{3001 + a} 150M AS:150 NM:{want_nm}",
              file=sys.stderr)
        sys.exit(1)
if n < 5: print(f"only {n} SNP-covering reads", file=sys.stderr); sys.exit(1)
EOF
python3 snpcheck.py proj.tsv diffs.txt || fail "a lifted SNP read is misplaced or has the wrong AS or NM"
# Both mates of every SNP fragment: on chrA, proper pair.
awk '$10 <= 1000 && 1000 < $11 { print $1 }' proj.tsv | sort -u > snp.names
bad=$(awk 'NR == FNR { n[$1] = 1; next } ($1 in n) && ($3 != "chrA" || int($7 / 2) % 2 != 1)' snp.names proj.tsv | wc -l)
[[ $bad -eq 0 ]] || fail "$((bad)) mates of SNP fragments are off chrA or not a proper pair"
# No alignment field names the ALT contig; XL records historical provenance.
awk -F '\t' '!/^@/ {for(i=1;i<=NF;i++) if($i !~ /^XL:Z:/) print $i}' proj.sam > proj.body
grep -q 'chrA_alt1' proj.body && fail "--alt-project wrote a record that names the ALT contig"
# Deletion: X:1500 is chrA 0-based 4500, so the D starts at 1-based 4501. Every
# mate with at least 10 ALT bases on each side of it must carry it there.
python3 - proj.tsv << 'EOF' || fail "a read spanning the ALT deletion lacks the 3D at chrA:4501"
import re, sys
n = 0
for l in open(sys.argv[1]):
    name, m, c, pos, q, cig, flag, AS, NM, a, b = l.split()
    if not (int(a) <= 1490 and int(b) >= 1510): continue
    n += 1
    ref, ok = 0, False
    for ln, op in re.findall(r"(\d+)([MIDNSHP=X])", cig):
        if op == "D" and int(ln) == 3 and int(pos) + ref == 4501: ok = True
        if op in "MDN=X": ref += int(ln)
    if c != "chrA" or not ok:
        print(f"{name}/{m}: {c}:{pos} {cig}", file=sys.stderr); sys.exit(1)
if n < 3: print(f"only {n} reads span the deletion", file=sys.stderr); sys.exit(1)
EOF
cmp -s proj.tsv projr.tsv || fail "reverse-complemented ALT (.alt FLAG 16) placed reads differently"

# Indel fixtures (single-end): every lifted read keeps both ALT indels at its
# lifted position, with NM equal to the indel bases (the reads are exact ALT
# substrings), rather than an ungapped run of mismatches. For the 1-5 bp
# balanced pairs the record also equals the default run's chrA record, apart
# from AS (the ALT score) and the tags that name the ALT hit (at 35 bp the
# default run soft-clips the read on chrA instead of paying for two gaps).
for ref in bal1 bal2 bal5 bal35 wide; do
    "$BIN" mem -t 2 --alt-project "$ref.fa" "$ref.fq" > "$ref.proj.sam" 2> "$ref.proj.log" \
        || fail "mem --alt-project on $ref failed"
    "$BIN" mem -t 2 "$ref.fa" "$ref.fq" > "$ref.def.sam" 2> "$ref.def.log" || fail "mem on $ref failed"
    python3 - "$ref.proj.sam" "$ref.def.sam" "$ref" << 'EOF' || fail "--alt-project misaligned a read across the ALT indels of $ref"
import re, sys
from samrecs import prim, tags
proj, dflt = prim(sys.argv[1]), prim(sys.argv[2])
if len(proj) != 3: print(f"{len(proj)} primary records, want 3", file=sys.stderr); sys.exit(1)
for name, t in proj.items():
    s, ins, dl = (int(x) for x in name.split("_")[1:])
    tg = tags(t)
    ops = re.findall(r"(\d+)([MIDNSHP=X])", t[5])
    ok = (t[2] == "chrA" and int(t[3]) == 3001 + s and (str(ins), "I") in ops and (str(dl), "D") in ops
          and all(op in "MID" for _, op in ops) and int(tg.get("NM", -1)) == ins + dl and tg.get("AS") == "150")
    if ok and ins == dl and ins <= 5:
        d = dflt[name]
        dtg = tags(d)
        ok = d[2:6] == t[2:6] and dtg.get("MD") == tg.get("MD") and dtg.get("NM") == tg.get("NM")
    if not ok:
        print(f"{name}: {t[2]}:{t[3]} {t[5]} NM:{tg.get('NM')} AS:{tg.get('AS')}", file=sys.stderr)
        sys.exit(1)
EOF
done

# Past the band: a 40I...40D ALT needs a band of 40, more than the 35 the
# re-alignment of a 150 bp read can have, so its reads stay on the ALT and their
# records equal the default run's.
"$BIN" mem -t 2 --alt-project bal40.fa bal40.fq > bal40.proj.sam 2> bal40.proj.log || fail "mem --alt-project on bal40 failed"
"$BIN" mem -t 2 bal40.fa bal40.fq > bal40.def.sam 2> bal40.def.log || fail "mem on bal40 failed"
grep -v '^@PG' bal40.proj.sam > bal40.proj.nopg
grep -v '^@PG' bal40.def.sam > bal40.def.nopg
cmp -s bal40.proj.nopg bal40.def.nopg || fail "--alt-project changed reads whose ALT indels are past the re-alignment band (bal40)"
grep -q $'\tchrA_alt1\t' bal40.def.nopg || fail "bal40: the default run has no ALT record (test is vacuous)"

# End fixture (single-end): for read e_<900 + k> (k = 0, 2, 4) the ALT hit is
# 100 - k bases at AS 100 - k and the chrA hit 150M at AS 90 + k (its first
# mismatch is at read offset 4 - k; AS is the local score, which leaves out the
# bases up to it, while the clipping penalty keeps the CIGAR end to end); the
# 50 + k bases the lifted hit leaves out are worth 50 + k, so it does not
# replace the chrA hit and is dropped. The read's only record is the default
# run's primary chrA record (150M at chrA:3001 + its ALT start),
# field for field and tag for tag, minus what the default derives from the ALT
# hit: its SA tag (to the default's ALT supplementary record, which is gone too)
# and its pa tag. The default run must have that supplementary record and those
# tags, so the comparison is not vacuous.
"$BIN" mem -t 2 --alt-project end.fa end.fq > end.proj.sam 2> end.proj.log || fail "mem --alt-project on end failed"
"$BIN" mem -t 2 end.fa end.fq > end.def.sam 2> end.def.log || fail "mem on end failed"
python3 - end.proj.sam end.def.sam << 'EOF' || fail "--alt-project let a shorter lifted hit replace a full chrA hit"
import sys
from samrecs import recs, tags
got, dflt = recs(sys.argv[1]), recs(sys.argv[2])
if len(got) != 3 or set(got) != set(dflt): sys.exit(f"{len(got)} reads, want the default's 3")
for name, rs in got.items():
    s = int(name.split("_")[1])
    dp = [t for t in dflt[name] if not int(t[1]) & 0x900]
    ds = [t for t in dflt[name] if int(t[1]) & 0x800 and t[2] == "chrA_alt1"]
    if len(dp) != 1 or len(ds) != 1: sys.exit(f"{name}: the default run has no ALT supplementary record (vacuous)")
    d = dp[0]
    if "SA" not in tags(d) or "pa" not in tags(d):
        sys.exit(f"{name}: the default primary record has no SA or pa tag (vacuous)")
    want = d[:11] + [x for x in d[11:] if not x.startswith(("SA:Z:", "pa:f:"))]
    if len(rs) != 1 or rs[0] != want or rs[0][2] != "chrA" or int(rs[0][3]) != 3001 + s or rs[0][5] != "150M":
        got_s = ["\t".join(t[1:6] + t[11:]) for t in rs]
        sys.exit(f"{name}: {got_s}, want one record, the default primary minus SA and pa: {want[1:6] + want[11:]}")
EOF

# Paralog fixture (single-end), the other direction: the lifted hit misses the
# read's first 4 bases (4S146M, AS 146) and the chrA hit at the same locus aligns
# all 150 (AS 120); 146 beats 120 + 4, so the lifted hit replaces it, and the read
# is written at chrA:3405 (3001 + 400 + 4) rather than at the paralog (chrA:8001,
# AS 130), where the default run places it.
"$BIN" mem -t 2 --alt-project par.fa par.fq > par.proj.sam 2> par.proj.log || fail "mem --alt-project on par failed"
"$BIN" mem -t 2 par.fa par.fq > par.def.sam 2> par.def.log || fail "mem on par failed"
python3 - par.proj.sam par.def.sam << 'EOF' || fail "--alt-project did not let a much better lifted hit that clips 4 bases replace a full chrA hit"
import sys
from samrecs import recs, prim, tags
got, d = recs(sys.argv[1]), prim(sys.argv[2])
if set(got) != {"v_400"}: sys.exit(f"reads {sorted(got)}, want v_400")
rs = got["v_400"]
p = [t for t in rs if not int(t[1]) & 0x900]
if not p: sys.exit("no primary record")
p = p[0]
if (p[2], p[3], p[5], tags(p).get("AS")) != ("chrA", "3405", "4S146M", "146") or int(p[4]) == 0:
    sys.exit(f"primary {p[2]}:{p[3]} {p[5]} MAPQ {p[4]} AS:{tags(p).get('AS')}, want chrA:3405 4S146M AS:146, MAPQ > 0")
if any("chrA_alt1" in "\t".join(x for x in t if not x.startswith("XL:Z:")) for t in rs): sys.exit("a record names the ALT")
dv = d["v_400"]
if (dv[2], dv[3]) != ("chrA", "8001"): sys.exit(f"the default places it at {dv[2]}:{dv[3]}, not the paralog (vacuous)")
EOF

# Tandem fixture (single-end): the lifted hit and the chrA placement one repeat
# unit away score the same, so both survive: MAPQ 0, XS equal to AS, and under
# -a a secondary at the other placement.
"$BIN" mem -t 2 -a --alt-project tan.fa tan.fq > tan.sam 2> tan.log || fail "mem -a --alt-project on tan failed"
BWAMEM3_ALT_PROJECT_TRACE=t_500 "$BIN" mem -t 2 -a --alt-project --alt-project-min-gap 2 tan.fa tan.fq > tan.gap.sam 2> tan.gap.log
python3 - << 'EOF'
rows = [l.rstrip().split("\t") for l in open("tan.gap.log") if l.startswith("[alt-project-trace]\t")]
assert rows and {r[1] for r in rows} == {"t_500"}, "trace must match exact names"
assert any(r[3] == "veto" for r in rows), "tandem competitor did not veto the lift"
assert any(r[3] == "gate" and "decision=refuse" in r for r in rows), "missing refused-candidate trace"
EOF
BWAMEM3_ALT_PROJECT_PAIR_GAP=1 BWAMEM3_ALT_PROJECT_TRACE=t_500 "$BIN" mem -t 2 -a --alt-project --alt-project-min-gap 2 tan.fa tan.fq > tan.pair.sam 2> tan.pair.log
python3 - << 'EOF'
from samrecs import records
rows = [l.rstrip().split("\t") for l in open("tan.pair.log") if l.startswith("[alt-project-trace]\t")]
assert any(r[3] == "gate" and "decision=defer_to_pair" in r for r in rows), "pair-confidence selector not enabled"
assert list(records("tan.pair.sam")) == list(records("tan.sam")), "pair-only confidence experiment changed single-end records"
EOF
python3 - tan.sam << 'EOF' || fail "--alt-project dropped a competing placement in a tandem repeat"
import sys
from samrecs import recs, tags
by_name = recs(sys.argv[1])
for rs in by_name.values():
    for t in rs:
        if int(t[1]) & 0x4: print(f"{t[0]}: unmapped", file=sys.stderr); sys.exit(1)
if len(by_name) != 4: print(f"{len(by_name)} reads, want 4", file=sys.stderr); sys.exit(1)
for name, rs in by_name.items():
    s = int(name.split("_")[1])
    p = [t for t in rs if not int(t[1]) & 0x900]
    if not p: print(f"{name}: no primary record", file=sys.stderr); sys.exit(1)
    tg = tags(p[0])
    pos150 = {(t[2], int(t[3])) for t in rs if not int(t[1]) & 0x800 and "AS:i:150" in t[11:]}
    if (len(p) != 1 or p[0][4] != "0" or tg.get("AS") != "150" or tg.get("XS") != "150"
            or pos150 != {("chrA", 3001 + s), ("chrA", 3101 + s)}):
        print(f"{name}: MAPQ {p[0][4]} AS:{tg.get('AS')} XS:{tg.get('XS')} at {sorted(pos150)}", file=sys.stderr)
        sys.exit(1)
EOF

# Clip fixture (single-end). Reads wholly inside the ALT-only 300 bp stay on the
# ALT. Reads straddling its edge keep their ALT hit on the ALT: their records
# equal the default run's and none carries a long insertion (an unguarded lift
# writes the ALT-only bases as one). Reads past it are at chrA or chrB at the
# lifted position, and none of their records names the ALT: the ALT hit went to
# chrA and merged with chrA's own. The default run names the ALT on those reads
# (as XA), so that last check fails there and is not vacuous. An unmapped read is
# a failure.
"$BIN" mem -t 2 --alt-project clip.fa clip.fq > clip.sam 2> clip.log || fail "mem --alt-project on clip failed"
"$BIN" mem -t 2 clip.fa clip.fq > clip.def.sam 2> clip.def.log || fail "mem on clip failed"
cat > clipcheck.py << 'EOF'
import re, sys
from samrecs import recs
got, dflt = recs(sys.argv[1]), recs(sys.argv[2])
n_in = n_edge = n_past = 0
for name, rs in got.items():
    s = int(name.split("_")[1])
    for t in rs:
        if int(t[1]) & 0x4: sys.exit(f"{name}: unmapped")
        if max([int(x) for x in re.findall(r"(\d+)I", t[5])] + [0]) >= 20: sys.exit(f"{name}: {t[5]}")
    p = [t for t in rs if not int(t[1]) & 0x900]
    if not p: sys.exit(f"{name}: no primary record")
    p = p[0]
    if s + 150 <= 300:
        n_in += 1
        if p[2] != "chrA_alt2": sys.exit(f"{name}: {p[2]}, want chrA_alt2")
    elif s < 300:
        n_edge += 1
        if not any(t[2] == "chrA_alt2" for t in rs): sys.exit(f"{name}: its ALT hit left the ALT")
        if rs != dflt[name]: sys.exit(f"{name}: records differ from the default run's")
    else:
        n_past += 1
        if p[2] not in ("chrA", "chrB") or int(p[3]) != 3001 + s - 300: sys.exit(f"{name}: {p[2]}:{p[3]}")
        if any("chrA_alt2" in "\t".join(x for x in t if not x.startswith("XL:Z:")) for t in rs): sys.exit(f"{name}: a record names chrA_alt2")
if n_in < 10 or n_edge < 10 or n_past < 10: sys.exit(f"only {n_in} clip / {n_edge} edge / {n_past} past reads")
EOF
python3 clipcheck.py clip.sam clip.def.sam || fail "--alt-project mishandled reads in, across or past an ALT clip"
# The same check on the default run must fail, and for the reason it exists:
# a past-the-clip read whose records name the ALT.
python3 clipcheck.py clip.def.sam clip.def.sam 2> clipvac.log && fail "the clip check passes on the default run (vacuous)"
grep -q "a record names chrA_alt2" clipvac.log || fail "the clip check fails on the default run for another reason (see clipvac.log)"

# Determinism: records are byte-identical (apart from @PG, whose CL records the
# arguments) across thread counts, under --fast (its pdqsort orders score ties
# differently), and with --dedup-reads on and off over input with every pair
# duplicated, so memoization replays records.
nopg() { grep -v '^@PG' "$1"; }
for t in 1 4; do
    "$BIN" mem -t "$t" -I 400,10 --alt-project fwd.fa r1.fq r2.fq > "det.t$t.sam" 2> "det.t$t.log" || fail "mem -t $t --alt-project failed"
    "$BIN" mem -t "$t" -I 400,10 --fast --alt-project fwd.fa r1.fq r2.fq > "detf.t$t.sam" 2> "detf.t$t.log" || fail "mem -t $t --fast --alt-project failed"
    BWAMEM3_DEDUP_READS_STATS=1 "$BIN" mem -t "$t" -I 400,10 --dedup-reads on --alt-project fwd.fa r1d.fq r2d.fq \
        > "detd.on.t$t.sam" 2> "detd.on.t$t.log" || fail "mem -t $t --dedup-reads on --alt-project failed"
    "$BIN" mem -t "$t" -I 400,10 --dedup-reads off --alt-project fwd.fa r1d.fq r2d.fq > "detd.off.t$t.sam" 2> "detd.off.t$t.log" \
        || fail "mem -t $t --dedup-reads off --alt-project failed"
done
cmp -s <(nopg det.t1.sam) <(nopg det.t4.sam) || fail "--alt-project output differs between -t 1 and -t 4"
cmp -s <(nopg det.t1.sam) <(nopg proj.sam) || fail "--alt-project output differs between -t 1 and -t 2"
cmp -s <(nopg detf.t1.sam) <(nopg detf.t4.sam) || fail "--fast --alt-project output differs between -t 1 and -t 4"
# --fast lifts too: its SNP reads pass the same placement check as the default
# mode's (which also rules out an empty or unmapped output), and no record names
# the ALT (a file, not a pipe into grep -q, whose early exit could SIGPIPE the
# producer and, under pipefail, hide a match).
summ detf.t1.sam > detf.tsv
python3 snpcheck.py detf.tsv diffs.txt || fail "--fast --alt-project did not lift the SNP reads to their chrA positions"
awk -F '\t' '!/^@/ {for(i=1;i<=NF;i++) if($i !~ /^XL:Z:/) print $i}' detf.t1.sam > detf.body
grep -q 'chrA_alt1' detf.body && fail "--fast --alt-project wrote a record that names the ALT contig"
# --dedup-reads: byte-identical on and off (its documented contract), and the
# memo must engage, or on == off is vacuous: its stats line reports duplicate
# pairs and the final state ON.
for t in 1 4; do
    cmp -s <(nopg "detd.on.t$t.sam") <(nopg "detd.off.t$t.sam") || fail "--alt-project output differs with --dedup-reads on and off (-t $t)"
    grep -q -E '^\[dedup-reads-stats\] pairs=[0-9]+ dup=[1-9][0-9]* .*final=ON' "detd.on.t$t.log" \
        || fail "--dedup-reads on --alt-project did not memoize any pair (-t $t, see detd.on.t$t.log)"
done
cmp -s <(nopg detd.on.t1.sam) <(nopg detd.on.t4.sam) || fail "--dedup-reads on --alt-project output differs between -t 1 and -t 4"
# The VERIFY diagnostic aligns each duplicate itself and aborts if its regions
# (lifted ones and their band floor included) differ from its representative's;
# its output must equal the off run's, and it must have compared duplicates.
BWAMEM3_DEDUP_READS_VERIFY=1 BWAMEM3_DEDUP_READS_STATS=1 "$BIN" mem -t 4 -I 400,10 --dedup-reads on --alt-project fwd.fa r1d.fq r2d.fq \
    > detd.verify.sam 2> detd.verify.log || fail "--dedup-reads VERIFY found a duplicate whose regions differ (see detd.verify.log)"
grep -q -E '^\[dedup-reads-stats\] pairs=[0-9]+ dup=[1-9][0-9]* ' detd.verify.log \
    || fail "--dedup-reads VERIFY compared no duplicate pair (see detd.verify.log)"
cmp -s <(nopg detd.verify.sam) <(nopg detd.off.t4.sam) || fail "--dedup-reads VERIFY output differs from the off run's"

# Lift provenance identifies the actual winning ALT, not an overlapping tile.
python3 - << 'EOF'
from samrecs import recs
rs = [r for values in recs("proj.sam").values() for r in values]
lifted = [r for r in rs if any(t == "XL:Z:chrA_alt1" for t in r[11:])]
assert lifted, "no lifted record has source ALT provenance"
assert all(r[2] == "chrA" for r in lifted), "source tag on an unlifted ALT record"
assert not any(t.startswith("XL:") for values in recs("def1.sam").values() for r in values for t in r[11:])
EOF

# Direct BAM must retain the same source tags and alignment content as SAM.
if command -v samtools > /dev/null 2>&1; then
    "$BIN" mem -t 2 -I 400,10 --alt-project --bam=0 fwd.fa r1.fq r2.fq > proj.bam 2> proj.bam.log || fail "direct BAM projection failed"
    samtools quickcheck proj.bam || fail "projected BAM failed quickcheck"
    samtools view proj.bam > proj.bam.sam
    python3 - << 'EOF'
from samrecs import recs
a, b = recs("proj.sam"), recs("proj.bam.sam")
assert a.keys() == b.keys()
for name in a:
    left = sorted((r[:11], sorted(r[11:])) for r in a[name])
    right = sorted((r[:11], sorted(r[11:])) for r in b[name])
    assert left == right, name
EOF
fi

# The gap is explicit, nonnegative, and meaningless without projection.
for gap in -1 1x 2147483648; do
    "$BIN" mem --alt-project --alt-project-min-gap "$gap" fwd.fa r1.fq r2.fq > /dev/null 2> gap.bad.log && fail "invalid gap $gap accepted"
    grep -q 'alt-project-min-gap' gap.bad.log || fail "invalid gap gave unrelated error"
done
"$BIN" mem --alt-project-min-gap 0 fwd.fa r1.fq r2.fq > /dev/null 2> gap.dep.log && fail "gap accepted without projection"
grep -q 'requires --alt-project' gap.dep.log || fail "gap dependency gave wrong error"
"$BIN" mem -t 2 -I 400,10 --alt-project --alt-project-min-gap 0 fwd.fa r1.fq r2.fq > gap.zero.sam 2> gap.zero.log || fail "zero gap failed"
cmp -s <(nopg gap.zero.sam) <(nopg proj.sam) || fail "zero gap changes projection output"
for t in 1 4; do
    BWAMEM3_DEDUP_READS_VERIFY=1 BWAMEM3_DEDUP_READS_STATS=1 "$BIN" mem -t "$t" -I 400,10 --alt-project --alt-project-min-gap 2 --dedup-reads on fwd.fa r1d.fq r2d.fq > "gap.t$t.sam" 2> "gap.t$t.log" || fail "positive gap failed with dedup verification"
    grep -q -E '^\[dedup-reads-stats\] pairs=[0-9]+ dup=[1-9][0-9]* ' "gap.t$t.log" || fail "positive gap compared no duplicate pairs"
done
cmp -s <(nopg gap.t1.sam) <(nopg gap.t4.sam) || fail "positive gap differs across thread counts"

# Confidence-only experiment must leave placement untouched, remain deterministic,
# replay duplicate pairs correctly, and agree between direct SAM and BAM output.
for t in 1 4; do
    BWAMEM3_ALT_PROJECT_PAIR_GAP=1 BWAMEM3_DEDUP_READS_VERIFY=1 BWAMEM3_DEDUP_READS_STATS=1 "$BIN" mem -t "$t" -I 400,10 --alt-project --alt-project-min-gap 2 --dedup-reads on fwd.fa r1d.fq r2d.fq > "pairgap.t$t.sam" 2> "pairgap.t$t.log"
    grep -q -E '^\[dedup-reads-stats\] pairs=[0-9]+ dup=[1-9][0-9]* ' "pairgap.t$t.log" || fail "pair confidence compared no duplicate pairs"
done
cmp -s <(nopg pairgap.t1.sam) <(nopg pairgap.t4.sam) || fail "pair confidence differs across threads"
BWAMEM3_ALT_PROJECT_PAIR_GAP=1 "$BIN" mem -t 2 -I 400,10 --alt-project --alt-project-min-gap 0 fwd.fa r1.fq r2.fq > pairgap.zero.sam 2> pairgap.zero.log
cmp -s <(nopg pairgap.zero.sam) <(nopg proj.sam) || fail "zero pair confidence gap changes records"
python3 - << 'EOF'
from samrecs import records
base, current = list(records("detd.off.t4.sam")), list(records("pairgap.t4.sam"))
assert len(base) == len(current)
for b, r in zip(base, current):
    assert b[:4] + b[5:11] == r[:4] + r[5:11], "pair confidence changed placement"
    assert int(r[4]) <= int(b[4]), "pair confidence raised MAPQ"
    assert [t for t in r[11:] if t[:3] in ("AS:", "XS:", "XL:", "NM:")] == [t for t in b[11:] if t[:3] in ("AS:", "XS:", "XL:", "NM:")]
EOF
if command -v samtools > /dev/null 2>&1; then
    BWAMEM3_ALT_PROJECT_PAIR_GAP=1 "$BIN" mem -t 4 -I 400,10 --alt-project --alt-project-min-gap 2 --dedup-reads on --bam=0 fwd.fa r1d.fq r2d.fq > pairgap.bam 2> pairgap.bam.log
    samtools quickcheck pairgap.bam
    samtools view pairgap.bam > pairgap.bam.sam
    python3 - << 'EOF'
from samrecs import records
left, right = list(records("pairgap.t4.sam")), list(records("pairgap.bam.sam"))
assert len(left) == len(right)
for a, b in zip(left, right):
    assert a[:11] == b[:11] and sorted(a[11:]) == sorted(b[11:])
EOF
fi

# Equal haplotype hypotheses at two paralogs must remain in paired competition.
# Ambiguous pairs have identical bytes but different source labels; a mate in
# distinct flanks resolves the locus. The asymmetric reference is a negative
# control: omitting one hypothesis must make the uncertainty check fail.
python3 - << 'EOF'
import random
rng = random.Random(775119)
rc = str.maketrans("ACGT", "TGCA")
references, placements = [], []
for family in range(8):
    a = "".join(rng.choice("ACGT") for _ in range(3000))
    b = list(a)
    for pos in range(650, 1150, 50):
        b[pos] = rng.choice([base for base in "ACGT" if base != a[pos]])
    for pos in range(1200, 1600, 10):
        b[pos] = rng.choice([base for base in "ACGT" if base != a[pos]])
    b = "".join(b)
    tile = b[600:780] + rng.choice([base for base in "ACGT" if base != b[780]]) + b[781:1200]
    references.append(f">A{family}\n{a}\n>B{family}\n{b}\n>hapA{family}\n{tile}\n")
    placements.append(f"hapA{family}\t0\tA{family}\t601\t60\t600M\t*\t0\t0\t*\t*\n")
    with open("sym.r1.fq", "a") as r1, open("sym.r2.fq", "a") as r2:
        for start in (650, 660, 670, 680):
            for informative, length in ((0, 400), (1, 750)):
                for source, primary in (("A", a), ("B", b)):
                    haplotype = primary[:600] + tile + primary[1200:]
                    end = start + length
                    reads = (haplotype[start:start + 151],
                             haplotype[end - 151:end].translate(rc)[::-1])
                    name = f"sym_{family}_{start}_{informative}_{source}"
                    for handle, read in zip((r1, r2), reads):
                        handle.write(f"@{name}\n{read}\n+\n{'I' * 151}\n")
open("asym.fa", "w").write("".join(references))
open("asym.fa.alt", "w").write("".join(placements))
for family in range(8):
    tile = references[family].splitlines()[-1]
    references.append(f">hapB{family}\n{tile}\n")
    placements.append(f"hapB{family}\t0\tB{family}\t601\t60\t600M\t*\t0\t0\t*\t*\n")
open("sym.fa", "w").write("".join(references))
open("sym.fa.alt", "w").write("".join(placements))
EOF
for ref in asym sym; do
    "$BIN" index "$ref.fa" > "$ref.index.log" 2>&1
    "$BIN" mem -t 2 -I 575,175,1100,100 --alt-project "$ref.fa" sym.r1.fq sym.r2.fq > "$ref.sam" 2> "$ref.log"
done
BWAMEM3_ALT_PROJECT_PAIR_GAP=1 "$BIN" mem -t 2 -I 575,175,1100,100 --alt-project --alt-project-min-gap 2 sym.fa sym.r1.fq sym.r2.fq > sym.pair.sam 2> sym.pair.log
python3 - << 'EOF'
from samrecs import records, tags

def check_symmetric(path):
    observed = {}
    for row in records(path):
        flag = int(row[1])
        if flag & 0x900:
            continue
        _, family, start, informative, source = row[0].split("_")
        mate = 1 if flag & 64 else 2
        key = row[0], mate
        assert key not in observed, "duplicate primary end"
        observed[key] = row
        assert not flag & 4, "unmapped symmetric read"
        if informative == "0":
            assert int(row[4]) < 20, "equally represented paired loci lost uncertainty"
        else:
            expected = int(start) + (0 if mate == 1 else 750 - 151)
            assert row[2] == source + family and int(row[3]) == expected + 1
            assert bool(flag & 16) == (mate == 2) and flag & 2
            assert int(row[4]) >= 20, "informative mate did not resolve the locus"
    assert len(observed) == 256, "missing primary ends"
    return observed

base = check_symmetric("sym.sam")
current = check_symmetric("sym.pair.sam")
assert base.keys() == current.keys()
for key, old in base.items():
    new = current[key]
    assert old[:4] + old[5:11] == new[:4] + new[5:11]
    assert int(new[4]) <= int(old[4])
    for tag in ("AS", "XS", "NM", "XL"):
        assert tags(old).get(tag) == tags(new).get(tag)
try:
    check_symmetric("asym.sam")
except AssertionError as error:
    assert str(error) == "equally represented paired loci lost uncertainty", str(error)
else:
    raise AssertionError("symmetric uncertainty check accepts a missing haplotype hypothesis")
EOF

# CLI restriction fixture: the original line is excluded and a later usable
# alignment for the same ALT is retained. This proves acceptance and fallback;
# mixed excluded/retained mapping is checked on its own fixture below.
for f in fwd.fa.*; do
    case "$f" in
        *.alt) ;;
        *) ln -s "$f" "exclude.fa${f#fwd.fa}" ;;
    esac
done
ln -s fwd.fa exclude.fa
cp fwd.fa.alt exclude.fa.alt
printf 'chrA_alt1\t0\tchrA\t1\t60\t1500M3D497M\t*\t0\t0\t*\t*\n' >> exclude.fa.alt
printf 'chrA\t3000\t3001\n' > exclude.bed
printf 'chrA\t0\t8000\n' > all.bed
printf 'chrA\t7000\t7001\n' > noeffect.bed
printf '# no intervals\nchrA\t1\t1\n' > empty.bed
printf 'chrA\t-1\t2\n' > invalid.bed
python3 - << 'EOF'
import gzip
with gzip.open("exclude.bed.gz", "wb") as out:
    out.write(open("exclude.bed", "rb").read())
raw = open("exclude.bed.gz", "rb").read()
open("truncated.bed.gz", "wb").write(raw[:-8])
EOF
"$BIN" mem --help > exclude.help 2>&1 || fail "mem --help failed"
"$BIN" mem -t 1 -I 400,10 --alt-project --alt-project-exclude-bed exclude.bed exclude.fa r1.fq r2.fq > exclude.sam 2> exclude.log \
    || fail "valid exclusion BED was rejected"
grep -q -- '--alt-project-exclude-bed FILE' exclude.help || fail "help omits exclusion BED option"
grep -q 'exclusion BED exclude.bed: input=1 merged=1 zero_width=0 excluded=1 retained=1$' exclude.log \
    || fail "restriction counts missing or wrong"
grep -q -v '^@' exclude.sam || fail "restricted run emitted no alignment records"
# det.t1.log is a real --alt-project log (its lift line), so the absence check
# below cannot pass on a missing or empty file.
grep -q 'alt-project: 1 ALT contigs' det.t1.log || fail "det.t1.log is not an --alt-project log"
grep -q 'exclusion BED' det.t1.log && fail "unrestricted --alt-project logged an exclusion BED"
"$BIN" mem -t 1 -I 400,10 --alt-project --alt-project-min-gap 0 --alt-project-exclude-bed exclude.bed.gz exclude.fa r1.fq r2.fq > exclude.zero.sam 2> exclude.zero.log \
    || fail "explicit gap zero or gzip BED was rejected"
grep -q 'exclusion BED exclude.bed.gz: input=1 merged=1 zero_width=0 excluded=1 retained=1$' exclude.zero.log \
    || fail "gzip restriction counts missing or wrong"
cmp -s <(nopg exclude.sam) <(nopg exclude.zero.sam) || fail "gap zero/gzip changed restricted records"
"$BIN" mem -t 1 --compat bwa-mem --compat-allow-divergent --alt-project --alt-project-exclude-bed exclude.bed exclude.fa r1.fq r2.fq > exclude.compat.sam 2> exclude.compat.log \
    || fail "compat override did not admit exclusion BED"

exclude_refuse() {
    local label=$1 message=$2
    shift 2
    if "$BIN" mem "$@" > "$label.sam" 2> "$label.log"; then
        fail "$label was accepted"
    fi
    grep -q -e "$message" "$label.log" || fail "$label gave the wrong error"
    if grep -q -v '^@' "$label.sam"; then
        fail "$label emitted alignment records"
    fi
}
exclude_refuse exclude.requires 'requires --alt-project' --alt-project-exclude-bed exclude.bed exclude.fa r1.fq r2.fq
exclude_refuse exclude.repeat 'only once' --alt-project --alt-project-exclude-bed exclude.bed --alt-project-exclude-bed exclude.bed exclude.fa r1.fq r2.fq
# getopt words this differently across C libraries ("option '--x' requires an
# argument", "option requires an argument: x"); both name the option and say
# it needs an argument, and neither says the option is unrecognized.
exclude_refuse exclude.argument 'alt-project-exclude-bed' --alt-project-exclude-bed
grep -q 'requires an argument' exclude.argument.log || fail "exclude.argument did not report the missing argument"
grep -q -i 'unrecognized' exclude.argument.log && fail "exclude.argument reported an unrecognized option"
exclude_refuse exclude.stdin 'requires a named BED file (not stdin)' --alt-project --alt-project-exclude-bed - exclude.fa r1.fq r2.fq
exclude_refuse exclude.emptyarg 'requires a named BED file (not stdin)' --alt-project --alt-project-exclude-bed '' exclude.fa r1.fq r2.fq
exclude_refuse exclude.j 'cannot be combined with -j' -j --alt-project --alt-project-exclude-bed exclude.bed exclude.fa r1.fq r2.fq
exclude_refuse exclude.meth 'cannot be combined with --meth' --meth --alt-project --alt-project-exclude-bed exclude.bed exclude.fa r1.fq r2.fq
exclude_refuse exclude.compat 'mutually exclusive: --alt-project' --compat bwa-mem --alt-project --alt-project-exclude-bed exclude.bed exclude.fa r1.fq r2.fq
exclude_refuse exclude.gap 'positive --alt-project-min-gap' --alt-project --alt-project-min-gap 1 --alt-project-exclude-bed exclude.bed exclude.fa r1.fq r2.fq
BWAMEM3_ALT_PROJECT_PAIR_GAP=1 exclude_refuse exclude.pairgap 'BWAMEM3_ALT_PROJECT_PAIR_GAP=1' --alt-project --alt-project-exclude-bed exclude.bed exclude.fa r1.fq r2.fq
exclude_refuse exclude.missing 'missing.bed' --alt-project --alt-project-exclude-bed missing.bed exclude.fa r1.fq r2.fq
# A mistyped BED path fails before the index is touched: pointing at an index
# that does not exist still yields the BED error, not an index error.
exclude_refuse exclude.missing.early '^\[E::main_mem\] --alt-project-exclude-bed missing.bed: cannot read: No such file' \
    --alt-project --alt-project-exclude-bed missing.bed no-such-index.fa r1.fq r2.fq
exclude_refuse exclude.invalid 'invalid.bed line 1' --alt-project --alt-project-exclude-bed invalid.bed exclude.fa r1.fq r2.fq
# Dictionary errors reach the user through the option's prefix with the BED's
# path and line: a contig the index lacks, an ALT contig, and a span past the
# contig's end.
printf 'chrZ\t0\t10\n' > unknown.bed
printf 'chrA_alt1\t0\t10\n' > altctg.bed
printf 'chrA\t0\t8001\n' > beyond.bed
exclude_refuse exclude.unknown '^\[E::main_mem\] --alt-project-exclude-bed: unknown.bed line 1: unknown reference contig' \
    --alt-project --alt-project-exclude-bed unknown.bed exclude.fa r1.fq r2.fq
exclude_refuse exclude.altctg '^\[E::main_mem\] --alt-project-exclude-bed: altctg.bed line 1: BED contig is an ALT contig' \
    --alt-project --alt-project-exclude-bed altctg.bed exclude.fa r1.fq r2.fq
exclude_refuse exclude.beyond '^\[E::main_mem\] --alt-project-exclude-bed: beyond.bed line 1: interval exceeds reference contig length' \
    --alt-project --alt-project-exclude-bed beyond.bed exclude.fa r1.fq r2.fq
exclude_refuse exclude.truncated 'truncated.bed.gz.*failure' --alt-project --alt-project-exclude-bed truncated.bed.gz exclude.fa r1.fq r2.fq
exclude_refuse exclude.empty 'no nonempty intervals' --alt-project --alt-project-exclude-bed empty.bed exclude.fa r1.fq r2.fq
exclude_refuse exclude.noeffect 'excludes no usable' --alt-project --alt-project-exclude-bed noeffect.bed exclude.fa r1.fq r2.fq
exclude_refuse exclude.all 'excludes all 2 usable ALT alignments' --alt-project --alt-project-exclude-bed all.bed exclude.fa r1.fq r2.fq
grep -q 'ordinary ALT-aware mapping, not primary-only' exclude.all.log || fail "all-excluded error does not explain omitting --alt-project"
# An empty BED over a .alt with no usable line keeps the ordinary zero-usable
# error (the order is: zero-usable, empty BED, no effect, all excluded).
exclude_refuse exclude.empty.unusable 'no ALT contig in bad.fa.alt has a usable alignment' --alt-project --alt-project-exclude-bed empty.bed bad.fa r1.fq r2.fq
grep -q 'no nonempty intervals' exclude.empty.unusable.log && fail "empty BED masked the zero-usable error"
exclude_refuse exclude.unusable 'no ALT contig in bad.fa.alt has a usable alignment' --alt-project --alt-project-exclude-bed exclude.bed bad.fa r1.fq r2.fq
grep -q 'excludes all' exclude.unusable.log && fail "an unusable sidecar was reported as all-excluded"

# Mixed exclusion fixture: two independent loci, each a copy (X1, X2 on chrA)
# with a paralog (Y1, Y2 on chrB) and an ALT haplotype carrying gene-converted
# paralog bases and a private SNP at ALT position 1000, like the fixture above.
# The ALT of X1 (a conventional *_alt name) is retained; the ALT of X2 (an
# HLA-style name with no ALT prefix or suffix, stored reverse-complemented,
# .alt FLAG 16) is excluded by a BED interval 100 bp from its footprint's end,
# away from every informative read: exclusion takes the whole primary footprint,
# not the read's aligned span. Two more intervals touch the retained footprint at
# its start and end without overlapping it. The oracle is a copy of the .alt
# whose excluded line has FLAG 0x4 added, decided here from the CIGAR's
# reference span (M/D/N/=/X) and half-open overlap, not by the aligner's parser,
# and loaded beside the same index files (symlinks, no rebuild) as mixmask.fa.
python3 - << 'EOF'
import random, re
rng = random.Random(77501)
B = "ACGT"
r = lambda n: "".join(rng.choice(B) for _ in range(n))
mut = lambda c: rng.choice([x for x in B if x != c])
RC = str.maketrans("ACGT", "TGCA")
rc = lambda s: s.translate(RC)[::-1]
def locus():
    X = r(2000)
    psv = sorted(rng.sample(range(200, 1800), 8))
    Y, H = list(X), list(X)
    for p in psv: Y[p] = mut(X[p])
    for p in psv[:6]: H[p] = Y[p]
    H[1000] = mut(X[1000])
    return X, "".join(Y), "".join(H)
X1, Y1, HR = locus()
X2, Y2, HE = locus()
chrA = r(3000) + X1 + r(3000) + X2 + r(3000)   # X1 at 3000, X2 at 8000 (0-based)
chrB = r(3000) + Y1 + r(3000) + Y2 + r(3000)
RN, EN = "chrA_KI000001v1_alt", "HLA-Z*01:01"
open("mix.fa", "w").write(f">chrA\n{chrA}\n>chrB\n{chrB}\n>{RN}\n{HR}\n>{EN}\n{rc(HE)}\n")
alt = [[RN, "0", "chrA", "3001", "60", "2000M", "*", "0", "0", "*", "*"],
       [EN, "16", "chrA", "8001", "60", "2000M", "*", "0", "0", "*", "*"]]
open("mix.fa.alt", "w").write("".join("\t".join(x) + "\n" for x in alt))
bed = [("chrA", 9900, 9910), ("chrA", 2990, 3000), ("chrA", 5000, 5010)]
open("mix.bed", "w").write("".join(f"{c}\t{s}\t{e}\n" for c, s, e in bed))
masked = []
for x in alt:
    a = int(x[3]) - 1
    b = a + sum(int(n) for n, op in re.findall(r"(\d+)([MIDNSHP=X])", x[5]) if op in "MDN=X")
    hit = any(c == x[2] and s < b and a < e for c, s, e in bed)
    masked.append([x[0], str(int(x[1]) | (4 if hit else 0))] + x[2:])
assert [x[1] for x in masked] == ["0", "20"], "the oracle must mask exactly the excluded line"
open("mixmask.fa.alt", "w").write("".join("\t".join(x) + "\n" for x in masked))
# Pairs from both haplotypes (R1 = ALT [s, s+150), R2 = [s+250, s+400)), the
# same pairs each followed by a duplicate (new name) for --dedup-reads, and
# single-end reads around the private SNP. Names: <m|s>_<r|e>_<ALT start>.
with open("mix1.fq", "w") as r1, open("mix2.fq", "w") as r2, \
        open("mixd1.fq", "w") as d1, open("mixd2.fq", "w") as d2:
    for tag, h in (("r", HR), ("e", HE)):
        for s in range(0, 1600, 20):
            f, q = h[s:s + 400], "I" * 150
            for o1, o2, sfx in ((r1, r2, ""), (d1, d2, ""), (d1, d2, "_dup")):
                o1.write(f"@m_{tag}_{s}{sfx}\n{f[:150]}\n+\n{q}\n")
                o2.write(f"@m_{tag}_{s}{sfx}\n{rc(f[-150:])}\n+\n{q}\n")
with open("mixse.fq", "w") as fq:
    for tag, h in (("r", HR), ("e", HE)):
        for s in range(700, 1001, 10):
            fq.write(f"@s_{tag}_{s}\n{h[s:s + 150]}\n+\n{'I' * 150}\n")
# The mixed-exclusion oracle. An informative read end has the private SNP at
# least 10 bases inside it. Each retained one is lifted: primary on chrA at
# 3001 + its ALT start, 150M, MAPQ > 0, XL:Z naming the retained ALT. No
# excluded one names the excluded ALT in XL:Z (it was not projected), and its
# records still name that ALT elsewhere (RNAME or XA: the unprojected ALT hit
# stays an ordinary ALT-aware candidate); where it is finally placed is not
# checked, as a primary hit may win. Both ALTs keep AH:* in the header. Prints
# the counts, and fails unless each is at least 10.
open("mixcheck.py", "w").write(f'''import sys
from samrecs import records
RN, EN = {RN!r}, {EN!r}
sq = [l.rstrip("\\n").split("\\t") for l in open(sys.argv[1]) if l.startswith("@SQ")]
if sorted((t[1][3:], "AH:*" in t) for t in sq) != sorted([("chrA", False), ("chrB", False), (RN, True), (EN, True)]):
    sys.exit(f"header ALT marking: {{sq}}")
ends = {{}}
for t in records(sys.argv[1]):
    ends.setdefault((t[0], int(t[1]) & 0xC0), []).append(t)
n = {{"retained": 0, "excluded": 0, "excluded_alt": 0}}
for (name, mate), rs in ends.items():
    tag, s = name.split("_")[1], int(name.split("_")[2])
    a = s + (250 if mate == 0x80 else 0)
    if not a + 10 <= 1000 < a + 140: continue
    p = [t for t in rs if not int(t[1]) & 0x900]
    if len(p) != 1: sys.exit(f"{{name}}/{{mate}}: {{len(p)}} primary records")
    p = p[0]
    xl = [x[5:] for t in rs for x in t[11:] if x.startswith("XL:Z:")]
    body = "\\t".join(x for t in rs for x in t if not x.startswith("XL:Z:"))
    if tag == "r":
        n["retained"] += 1
        if (p[2], int(p[3]), p[5]) != ("chrA", 3001 + a, "150M") or int(p[4]) == 0 or "XL:Z:" + RN not in p[11:]:
            sys.exit(f"{{name}}/{{mate}}: retained ALT read at {{p[2]}}:{{p[3]}} {{p[5]}} MAPQ {{p[4]}}, want chrA:{{3001 + a}} 150M MAPQ > 0 XL:Z:{{RN}}")
    else:
        assert 8000 + a + 150 <= 9900, "an informative excluded read overlaps the BED (fixture)"
        n["excluded"] += 1
        if EN in xl: sys.exit(f"{{name}}/{{mate}}: projected from the excluded ALT ({{p[2]}}:{{p[3]}})")
        if EN in body: n["excluded_alt"] += 1
print(" ".join(f"{{k}}={{v}}" for k, v in n.items()))
if min(n.values()) < 10: sys.exit(f"too few informative records: {{n}}")
''')
# Runtime and oracle outputs must be equal: header lines in order, @PG apart
# from its CL (the argv), and every record, field and tag in order. A decoded
# BAM's samtools @PG is dropped.
open("samecmp.py", "w").write('''import sys
def split(path):
    head, pg, body = [], [], []
    for l in open(path):
        l = l.rstrip("\\n")
        if l.startswith("@PG"):
            if "\\tPN:samtools" not in l: pg.append([x for x in l.split("\\t") if not x.startswith("CL:")])
        elif l.startswith("@"): head.append(l)
        else: body.append(l)
    return head, pg, body
a, b = split(sys.argv[1]), split(sys.argv[2])
if a[0] != b[0]: sys.exit("header dictionaries differ")
if a[1] != b[1]: sys.exit(f"@PG differs beyond its CL: {a[1]} vs {b[1]}")
if not a[2]: sys.exit("no records")
for i, (x, y) in enumerate(zip(a[2], b[2])):
    if x != y: sys.exit(f"record {i + 1} differs:\\n  {x}\\n  {y}")
if len(a[2]) != len(b[2]): sys.exit(f"{len(a[2])} records vs {len(b[2])}")
''')
EOF
"$BIN" index mix.fa > idx.mix.log 2>&1 || fail "index mix failed"
for f in mix.fa.*; do
    case "$f" in
        *.alt) ;;
        *) ln -s "$f" "mixmask.fa${f#mix.fa}" ;;
    esac
done
ln -s mix.fa mixmask.fa
mixrun() { # out ref se|pe|dup extra...
    local out=$1 ref=$2 mode=$3
    shift 3
    case "$mode" in
        se) set -- "$@" "$ref.fa" mixse.fq ;;
        pe) set -- -I 400,10 "$@" "$ref.fa" mix1.fq mix2.fq ;;
        dup) set -- -I 400,10 "$@" "$ref.fa" mixd1.fq mixd2.fq ;;
    esac
    "$BIN" mem -K 20000 "$@" > "$out" 2> "$out.log" || fail "mem $* failed (see $out.log)"
}
# Negative controls: unrestricted projection lifts the excluded ALT's reads, and
# the default (no projection) does not lift the retained ALT's, so the oracle
# check fails on each, for that reason.
for mode in se pe; do
    mixrun "mix.unr.$mode.sam" mix "$mode" -t 2 --alt-project
    python3 mixcheck.py "mix.unr.$mode.sam" > /dev/null 2> "mix.unr.$mode.chk" \
        && fail "the mixed-exclusion check passes on unrestricted --alt-project ($mode, vacuous)"
    grep -q 'projected from the excluded ALT' "mix.unr.$mode.chk" \
        || fail "the mixed-exclusion check fails on unrestricted --alt-project for another reason (see mix.unr.$mode.chk)"
done
mixrun mix.off.pe.sam mix pe -t 2
python3 mixcheck.py mix.off.pe.sam > /dev/null 2> mix.off.pe.chk && fail "the mixed-exclusion check passes without projection (vacuous)"
grep -q 'retained ALT read at' mix.off.pe.chk || fail "the mixed-exclusion check fails without projection for another reason (see mix.off.pe.chk)"
# Runtime BED and FLAG-masked sidecar, each checked against the oracle and equal
# to each other: SE and PE at -t 1, 2 and 4, --fast at -t 1 and 4, and
# --dedup-reads on and off over duplicated pairs (memo shown to engage), all
# with -K fixed (and -I for pairs).
: > mix.counts
mixpair() { # label se|pe|dup extra...
    local label=$1 mode=$2 counts
    shift 2
    mixrun "mix.$label.rt.sam" mix "$mode" "$@" --alt-project --alt-project-exclude-bed mix.bed
    mixrun "mix.$label.mask.sam" mixmask "$mode" "$@" --alt-project
    grep -q 'exclusion BED mix.bed: input=3 merged=3 zero_width=0 excluded=1 retained=1$' "mix.$label.rt.sam.log" \
        || fail "mixed $label: restriction counts missing or wrong"
    grep -q 'exclusion BED' "mix.$label.mask.sam.log" && fail "mixed $label: the masked-sidecar run logged a BED"
    counts=$(python3 mixcheck.py "mix.$label.rt.sam" 2> "mix.$label.chk") \
        || fail "mixed $label: runtime exclusion failed the oracle (see mix.$label.chk)"
    echo "$label $counts" >> mix.counts
    python3 samecmp.py "mix.$label.rt.sam" "mix.$label.mask.sam" 2> "mix.$label.cmp" \
        || fail "mixed $label: runtime BED and FLAG-masked sidecar outputs differ (see mix.$label.cmp)"
}
for mode in se pe; do
    for t in 1 2 4; do
        mixpair "$mode.t$t" "$mode" -t "$t"
    done
    for t in 2 4; do
        python3 samecmp.py "mix.$mode.t1.rt.sam" "mix.$mode.t$t.rt.sam" || fail "mixed $mode: runtime output differs between -t 1 and -t $t"
    done
    for t in 1 4; do
        mixpair "$mode.fast.t$t" "$mode" -t "$t" --fast
    done
done
for t in 1 4; do
    BWAMEM3_DEDUP_READS_STATS=1 mixpair "dup.on.t$t" dup -t "$t" --dedup-reads on
    grep -q -E '^\[dedup-reads-stats\] pairs=[0-9]+ dup=[1-9][0-9]* .*final=ON' "mix.dup.on.t$t.rt.sam.log" \
        || fail "mixed --dedup-reads on did not memoize any pair (-t $t)"
    mixpair "dup.off.t$t" dup -t "$t" --dedup-reads off
    python3 samecmp.py "mix.dup.on.t$t.rt.sam" "mix.dup.off.t$t.rt.sam" || fail "mixed: --dedup-reads on and off differ (-t $t)"
done
# Projection off: the masked sidecar is the original for ALT-aware mapping (only
# QNAMEs mark ALT contigs), and the output names both ALTs (not vacuous).
for mode in se pe; do
    [[ $mode == pe ]] || mixrun mix.off.se.sam mix se -t 2
    mixrun "mixmask.off.$mode.sam" mixmask "$mode" -t 2
    python3 samecmp.py "mix.off.$mode.sam" "mixmask.off.$mode.sam" || fail "masked sidecar changes ALT-aware mapping ($mode)"
    grep -v '^@' "mix.off.$mode.sam" > "mix.off.$mode.body"
    grep -q 'chrA_KI000001v1_alt' "mix.off.$mode.body" || fail "the projection-off run names no retained ALT (vacuous, $mode)"
    grep -q 'HLA-Z\*01:01' "mix.off.$mode.body" || fail "the projection-off run names no excluded ALT (vacuous, $mode)"
done
# Direct BAM: decoded runtime BAM equals the runtime SAM, and the oracle's BAM.
# Without samtools this check is OPEN (reported, not passed);
# ALT_PROJECT_REQUIRE_SAMTOOLS=1 makes that a failure.
if command -v samtools > /dev/null 2>&1; then
    mixrun mix.rt.bam mix pe -t 2 --bam=0 --alt-project --alt-project-exclude-bed mix.bed
    mixrun mix.mask.bam mixmask pe -t 2 --bam=0 --alt-project
    for b in rt mask; do
        samtools quickcheck "mix.$b.bam" || fail "mixed $b BAM failed quickcheck"
        samtools view -h "mix.$b.bam" > "mix.$b.bam.sam" || fail "samtools view mix.$b.bam failed"
    done
    python3 samecmp.py mix.pe.t2.rt.sam mix.rt.bam.sam || fail "mixed runtime SAM and decoded BAM differ"
    python3 samecmp.py mix.rt.bam.sam mix.mask.bam.sam || fail "mixed runtime and masked-sidecar BAMs differ"
    mix_bam="BAM parity checked"
elif [[ ${ALT_PROJECT_REQUIRE_SAMTOOLS:-0} == 1 ]]; then
    fail "samtools not on PATH (ALT_PROJECT_REQUIRE_SAMTOOLS=1)"
else
    echo "OPEN: samtools not on PATH; mixed-exclusion direct BAM parity not checked" >&2
    mix_bam="BAM parity OPEN (no samtools)"
fi
mix_min=$(awk '{ for (i = 2; i <= NF; i++) { split($i, kv, "="); v = kv[2] + 0; if (m == "" || v < m) m = v } } END { print m }' mix.counts)

# Refusals
"$BIN" mem -j --alt-project fwd.fa r1.fq r2.fq > /dev/null 2> j.log && fail "-j --alt-project was accepted"
grep -q -e "cannot be combined with -j (it needs" j.log || fail "-j --alt-project gave the wrong error"
"$BIN" mem --meth --alt-project fwd.fa r1.fq r2.fq > /dev/null 2> m.log && fail "--meth --alt-project was accepted"
grep -q -e "cannot be combined with --meth" m.log || fail "--meth --alt-project gave the wrong error"
"$BIN" mem --alt-project noalt.fa r1.fq r2.fq > /dev/null 2> n.log && fail "--alt-project without a .alt was accepted"
grep -q "needs noalt.fa.alt" n.log || fail "--alt-project without a .alt gave the wrong error"
"$BIN" mem --alt-project bad.fa r1.fq r2.fq > /dev/null 2> b.log && fail "--alt-project with no usable .alt line was accepted"
grep -q "no ALT contig in bad.fa.alt has a usable alignment" b.log || fail "--alt-project with no usable .alt line gave the wrong error"
# The --compat refusal (both targets) is checked in compat_byte_identical.sh,
# with the other levers --compat refuses; only the override is checked here.
"$BIN" mem --compat bwa-mem --compat-allow-divergent --alt-project fwd.fa r1.fq r2.fq > /dev/null 2> ca.log \
    || fail "--compat-allow-divergent did not admit --alt-project (see ca.log)"

echo "PASS: alt_project ($((n_snp)) SNP reads on chrA as proper pairs, FLAG 16 identical, ALT indels kept across a lift, past-band indels kept on the ALT, full primary hits kept, a clearly better clipped lifted hit placed, tandem placements kept, ALT clip kept on the ALT, byte-identical across -t, --fast and --dedup-reads; mixed exclusion equals its FLAG-masked oracle in $(($(wc -l < mix.counts))) runs, at least $((mix_min)) informative records of each kind per run, $mix_bam)"
