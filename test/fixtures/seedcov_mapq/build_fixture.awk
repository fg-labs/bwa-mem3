# Deterministically construct the reference and reads for
# test/regression/seedcov_mapq_legacy.sh. Every base is sliced from the
# committed phix.fa and every substitution is a fixed rotation (A>C>G>T>A), so
# the fixture is byte-identical across awk implementations -- no PRNG.
#
# Run as:
#   mawk -v W=<dir> -f build_fixture.awk test/fixtures/phix.fa
# Writes <dir>/ref.fa, <dir>/se.fq, <dir>/r1.fq and <dir>/r2.fq.
#
# Reference: two contigs.
#   phix -- the whole phiX genome.
#   dup  -- phix[1000, 2700) (0-based, half-open) with every 41st base rotated.
# Every read is sliced from phix[1000, 2700), so it has a best hit on `phix`
# and a near-equal second hit on `dup`. That makes sub > 0 and keeps the
# legacy MAPQ (-Q 0: 30 * (1 - sub/score) * log(seedcov)) below the cap of 60,
# so the emitted MAPQ actually depends on seedcov.
#
# Reads: one per CASES entry, "<offset>:<length>:<edits>[:rc]", where the read
# is phix[1000 + offset, +length) with <edits> applied in read coordinates:
#   DEL<p>  a 2 bp deletion after read base p (applied first)
#   INS<p>  a 2 bp insertion ("GA") after read base p (applied first)
#   <p>     a substitution (rotation) at read base p
# and ":rc" reverse-complements the result. The cases are chosen so that the
# seedcov recount sites a phiX-sized fixture can reach each finish some region:
# exact reads (no extension), substitutions near either end (the right ungapped
# fast path, and the certified-band probe), a middle substitution (two seeds),
# indels (banded SW at both 8- and 16-bit widths; the long reads exceed the
# 8-bit envelope), and a left indel with a right-side substitution (the right
# fast path that runs only after the left SW). Not reached: the left ungapped
# fast-path recount in stage_seed_extension (a->qe/re are still unset there,
# so its H0_ guard never passes) and the scalar-tier recounts (they need
# extensions longer than 32 kb).
#
# se.fq holds the reads. r1.fq/r2.fq pair each read with a 100 bp mate on the
# other strand so every pair is a concordant FR pair with a 1000 bp insert:
# for a forward read starting at x the mate is revcomp(phix[x + 900, +100)),
# and for a ":rc" read ending at y the mate is phix[y - 1000, +100) forward.

/^>/ { next }
     { s = s toupper($0) }

function rot(b) { return b == "A" ? "C" : b == "C" ? "G" : b == "G" ? "T" : "A" }
function comp(b) { return b == "A" ? "T" : b == "C" ? "G" : b == "G" ? "C" : "A" }
function rc(x,    o, i) { o = ""; for (i = length(x); i >= 1; i--) o = o comp(substr(x, i, 1)); return o }
function sub1(x, p) { return substr(x, 1, p) rot(substr(x, p + 1, 1)) substr(x, p + 2) }
function subs(x, list,    n, a, i) {
    n = split(list, a, ",")
    for (i = 1; i <= n; i++) if (a[i] != "") x = sub1(x, a[i] + 0)
    return x
}
function qual(n,    o) { o = ""; while (n-- > 0) o = o "I"; return o }
function fq(f, name, x) { printf "@%s\n%s\n+\n%s\n", name, x, qual(length(x)) > f }

END {
    CASES = "0:100:|60:100:5|120:100:94|180:100:5,94|240:100:50|300:100:DEL40|360:100:INS60" \
            "|420:100:3,30,70,96|480:100:5,94:rc|540:100:DEL40:rc|600:100:12,88" \
            "|700:400:20|750:400:380|800:400:20,380|850:400:DEL30|900:400:INS370" \
            "|950:600:30,570|1000:600:DEL25:rc|700:250:8,240|760:250:DEL20" \
            "|1100:100:DEL12,94|1150:100:DEL8,90|1200:150:INS10,140" \
            "|1250:150:DEL12,60|1300:150:DEL10,45,100|1350:150:DEL70,30"
    ref = W "/ref.fa"; se = W "/se.fq"; r1 = W "/r1.fq"; r2 = W "/r2.fq"
    printf ">phix\n%s\n", s > ref
    d = substr(s, 1001, 1700); dd = ""
    for (i = 1; i <= length(d); i++) { b = substr(d, i, 1); dd = dd ((i % 41) == 20 ? rot(b) : b) }
    printf ">dup\n%s\n", dd > ref
    n = split(CASES, C, "|")
    for (k = 1; k <= n; k++) {
        m = split(C[k], f, ":")
        off = 1000 + f[1]; L = f[2] + 0; ed = f[3]; isrc = (m >= 4 && f[4] == "rc")
        x = substr(s, off + 1, L)
        if (ed ~ /^DEL/) {
            p = substr(ed, 4) + 0
            x = substr(s, off + 1, p) substr(s, off + p + 3, L - p)
            sub(/^DEL[0-9]+,?/, "", ed)
        } else if (ed ~ /^INS/) {
            p = substr(ed, 4) + 0
            x = substr(s, off + 1, p) "GA" substr(s, off + p + 1, L - 2 - p)
            sub(/^INS[0-9]+,?/, "", ed)
        }
        x = subs(x, ed)
        if (isrc) x = rc(x)
        name = "c" (k - 1)
        fq(se, name, x)
        fq(r1, name, x)
        if (isrc) fq(r2, name, substr(s, off + L - 1000 + 1, 100))
        else      fq(r2, name, rc(substr(s, off + 900 + 1, 100)))
    }
}
