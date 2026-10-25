#!/usr/bin/env bash
# test/regression/seedcov_mapq_legacy.sh
#
# Regression: the legacy MAPQ model (-Q 0) must keep reading a correct seedcov,
# and the default model (-Q 50) must not change.
#
# mem_alnreg_t::seedcov (the summed length of the chain's seeds that lie fully
# inside the region) has one consumer: mem_approx_mapq_se's legacy branch,
#
#     mapq = 30 * (1 - sub/score) * log(seedcov)
#
# taken only when mapQ_coef_len <= 0 (-Q 0 or negative). The default (-Q 50)
# never reads it, so every seedcov recount site in the extension driver is
# gated on mem_want_seedcov(opt) -- the same predicate the MAPQ branch tests --
# and the default run skips the seed walk. The risk that gate creates is
# specific to -Q 0: a site whose gate is wrong, or a new finalize path that
# forgets to recount, leaves seedcov at 0 or stale, and -Q 0 MAPQs change while
# every default-Q byte-identity gate stays green. This script pins both.
#
#   1. Goldens: QNAME, FLAG and MAPQ of every record must equal the values
#      below, at -Q 0 and at the default -Q, for single-end reads, for pairs
#      left to infer their insert size (too few pairs, so each end is scored on
#      its own), and for pairs with the insert size pinned by -I (so the
#      paired-end MAPQ path runs; every record is asserted to be a proper pair).
#      The -Q 0 legs also run under --no-band-cert, against the same golden.
#      The goldens were produced by a build without the gate (seedcov recounted
#      unconditionally). The fixture is built so the -Q 0 MAPQ depends on
#      seedcov: every read has a near-equal second hit on a planted `dup` contig,
#      so sub > 0 and the MAPQ sits below the cap of 60.
#
#      Coverage of the 11 recount sites in src/bwamem.cpp: forcing any one of
#      the no-extension, right ungapped fast-path (staging and post-left-SW),
#      certified-band probe, or 8-/16-bit banded-SW accept sites off changes the
#      output (the 16-bit ones under --no-band-cert). Three sites are covered
#      only by the source guard in (3): the left ungapped fast-path recount in
#      stage_seed_extension, which cannot run (a->qe/re are still unset there, so
#      its H0_ guard never passes), and the two scalar-tier recounts, which need
#      an extension longer than 32 kb.
#   2. Non-vacuity: the -Q 0 MAPQs must differ from the default -Q MAPQs (so the
#      legacy branch really ran) and at least half the SE records must have a
#      MAPQ strictly between 0 and 60 (so seedcov is not masked by the cap or by
#      sub >= score).
#   3. Source guard: every seedcov recount in src/bwamem.cpp (`seedcov = 0` /
#      `seedcov +=`) must sit inside the braced block of an `if` whose condition
#      begins with want_seedcov (`if (want_seedcov && ...` or `if (want_seedcov)`),
#      and there must be exactly 11 recount sites (counted by
#      their `seedcov = 0` loop headers). The check tracks brace scope (comments
#      and string / char literals stripped), so a recount after the guarded
#      block closes, or in its `else`, is rejected; a braceless want_seedcov
#      `if` authorizes nothing, and so does an `if` that only mentions it
#      (`!want_seedcov`, `x || want_seedcov`). It requires the file's braces to
#      balance (so a lexing slip cannot silently shift scopes) and self-tests
#      first on an ungated recount, a recount after a closed guarded block, one
#      in the guarded `if`'s `else`, and one under a negated want_seedcov.
#
# Regenerating the goldens (only after a deliberate MAPQ change): build a binary
# without the gate, run the legs below on the fixture at the matching -Q, and
# paste `cut -f1,2,5` of the records, space-joined, into the expect_* strings.
#
# Fixture (deterministic, no PRNG): test/fixtures/seedcov_mapq/build_fixture.awk
# slices the reference and reads from the committed phix.fa; see its header.
#
# Inputs (env vars):
#   BWA_MEM3          -- path to the bwa-mem3 binary under test
#   SEEDCOV_PHIX_FA   -- path to test/fixtures/phix.fa (the reference source)
#   SEEDCOV_WORK_DIR  -- directory for fixture-private intermediates
#   AWK               -- optional; defaults to mawk, then gawk, then awk (the
#                        generator is byte-identical across them)

set -Eeuo pipefail
# shellcheck disable=SC2154  # rc is assigned by `rc=$?` at the start of the trap body
trap 'rc=$?; echo "FAIL: unexpected command failure at line $LINENO (exit $rc)" >&2; exit "$rc"' ERR

# Checked explicitly rather than with ${VAR:?}: that exits before the ERR trap
# runs, so the failure would carry no FAIL: line.
for v in BWA_MEM3 SEEDCOV_PHIX_FA SEEDCOV_WORK_DIR; do
    if [[ -z "${!v:-}" ]]; then
        echo "FAIL: $v must be set" >&2
        exit 1
    fi
done
AWK="${AWK:-$(command -v mawk || command -v gawk || command -v awk)}"

here="$(cd "$(dirname "$0")" && pwd)"
gen="$here/../fixtures/seedcov_mapq/build_fixture.awk"
src="$here/../../src/bwamem.cpp"
[[ -f "$gen" ]] || {
    echo "FAIL: fixture generator not found at $gen" >&2
    exit 1
}
[[ -f "$src" ]] || {
    echo "FAIL: cannot find src/bwamem.cpp at $src" >&2
    exit 1
}

W="$SEEDCOV_WORK_DIR"
mkdir -p "$W"
rm -f "$W"/ref.fa* "$W"/se.fq "$W"/r1.fq "$W"/r2.fq
"$AWK" -v W="$W" -f "$gen" "$SEEDCOV_PHIX_FA"
"$BWA_MEM3" index "$W/ref.fa" > "$W/index.log" 2>&1

# QNAME FLAG MAPQ of every record, space-joined on one line, into <out>.
fields() { # <out> <mem args...>
    local out="$1"
    shift
    if ! "$BWA_MEM3" mem "$@" > "$out.sam" 2> "$out.log"; then
        echo "FAIL: bwa-mem3 mem $* exited nonzero; last lines of its log:" >&2
        tail -n 20 "$out.log" >&2
        exit 1
    fi
    grep -v '^@' "$out.sam" | cut -f1,2,5 | tr '\t' ' ' | paste -sd' ' - > "$out"
}

# Produced by a build without the seedcov gate (recount unconditional).
expect_se_q0="c0 0 14 c1 0 16 c2 0 15 c3 0 16 c4 0 15 c5 0 24 c6 0 16 c7 0 60 c8 16 15"
expect_se_q0+=" c9 16 25 c10 0 14 c11 0 23 c12 0 23 c13 0 21 c14 0 20 c15 0 23 c16 0 25 c17 16 24"
expect_se_q0+=" c18 0 20 c19 0 21 c20 0 15 c21 0 20 c22 0 17 c23 0 17 c24 0 17 c25 0 24"
expect_se_def="c0 0 26 c1 0 27 c2 0 25 c3 0 26 c4 0 25 c5 0 35 c6 0 24 c7 0 60 c8 16 24"
expect_se_def+=" c9 16 35 c10 0 24 c11 0 60 c12 0 60 c13 0 60 c14 0 60 c15 0 60 c16 0 60 c17 16 60"
expect_se_def+=" c18 0 52 c19 0 51 c20 0 22 c21 0 29 c22 0 32 c23 0 30 c24 0 29 c25 0 40"
expect_pe_q0="c0 97 14 c0 145 14 c1 97 16 c1 145 19 c2 97 15 c2 145 14 c3 97 16 c3 145 22 c4 97 15"
expect_pe_q0+=" c4 145 14 c5 97 24 c5 145 22 c6 97 16 c6 145 14 c7 97 60 c7 145 22 c8 81 15 c8 161 60"
expect_pe_q0+=" c9 81 25 c9 161 60 c10 97 14 c10 145 14 c11 97 23 c11 145 14 c12 97 23 c12 145 60 c13 97 21"
expect_pe_q0+=" c13 145 60 c14 97 20 c14 145 60 c15 97 23 c15 145 60 c16 97 25 c16 145 60 c17 81 24 c17 161 14"
expect_pe_q0+=" c18 97 20 c18 145 14 c19 97 21 c19 145 60 c20 97 15 c20 145 60 c21 97 20 c21 145 60 c22 97 17"
expect_pe_q0+=" c22 145 60 c23 97 17 c23 145 60 c24 97 17 c24 145 60 c25 97 24 c25 145 60"
expect_pe_def="c0 97 26 c0 145 26 c1 97 27 c1 145 33 c2 97 25 c2 145 26 c3 97 26 c3 145 38 c4 97 25"
expect_pe_def+=" c4 145 26 c5 97 35 c5 145 38 c6 97 24 c6 145 26 c7 97 60 c7 145 38 c8 81 24 c8 161 60"
expect_pe_def+=" c9 81 35 c9 161 60 c10 97 24 c10 145 26 c11 97 60 c11 145 26 c12 97 60 c12 145 60 c13 97 60"
expect_pe_def+=" c13 145 60 c14 97 60 c14 145 60 c15 97 60 c15 145 60 c16 97 60 c16 145 60 c17 81 60 c17 161 26"
expect_pe_def+=" c18 97 52 c18 145 26 c19 97 51 c19 145 60 c20 97 22 c20 145 60 c21 97 29 c21 145 60 c22 97 32"
expect_pe_def+=" c22 145 60 c23 97 30 c23 145 60 c24 97 29 c24 145 60 c25 97 40 c25 145 60"
expect_pei_q0="c0 99 54 c0 147 54 c1 99 56 c1 147 59 c2 99 55 c2 147 54 c3 99 56 c3 147 60 c4 99 55"
expect_pei_q0+=" c4 147 54 c5 99 60 c5 147 60 c6 99 56 c6 147 54 c7 99 60 c7 147 60 c8 83 55 c8 163 60"
expect_pei_q0+=" c9 83 60 c9 163 60 c10 99 54 c10 147 54 c11 99 60 c11 147 54 c12 99 60 c12 147 60 c13 99 60"
expect_pei_q0+=" c13 147 60 c14 99 60 c14 147 60 c15 99 60 c15 147 60 c16 99 60 c16 147 60 c17 83 60 c17 163 54"
expect_pei_q0+=" c18 99 60 c18 147 54 c19 99 60 c19 147 60 c20 99 55 c20 147 60 c21 99 60 c21 147 60 c22 99 57"
expect_pei_q0+=" c22 147 60 c23 99 57 c23 147 60 c24 99 57 c24 147 60 c25 99 60 c25 147 60"
expect_pei_def="c0 99 60 c0 147 60 c1 99 60 c1 147 60 c2 99 60 c2 147 60 c3 99 60 c3 147 60 c4 99 60"
expect_pei_def+=" c4 147 60 c5 99 60 c5 147 60 c6 99 60 c6 147 60 c7 99 60 c7 147 60 c8 83 60 c8 163 60"
expect_pei_def+=" c9 83 60 c9 163 60 c10 99 60 c10 147 60 c11 99 60 c11 147 60 c12 99 60 c12 147 60 c13 99 60"
expect_pei_def+=" c13 147 60 c14 99 60 c14 147 60 c15 99 60 c15 147 60 c16 99 60 c16 147 60 c17 83 60 c17 163 60"
expect_pei_def+=" c18 99 60 c18 147 60 c19 99 60 c19 147 60 c20 99 60 c20 147 60 c21 99 60 c21 147 60 c22 99 60"
expect_pei_def+=" c22 147 60 c23 99 60 c23 147 60 c24 99 60 c24 147 60 c25 99 60 c25 147 60"

check() { # <label> <expected> <got-file>
    local label="$1" want="$2" got
    got="$(cat "$3")"
    if [[ "$got" != "$want" ]]; then
        echo "FAIL: QNAME/FLAG/MAPQ differ from the golden under '$label' (< expected, > got):" >&2
        diff <(printf '%s\n' "$want" | xargs -n 3) <(printf '%s\n' "$got" | xargs -n 3) | head -20 >&2 || true
        echo "      If this MAPQ change is deliberate, regenerate the goldens from a build" >&2
        echo "      without the seedcov gate at the matching -Q (see the script header)." >&2
        exit 1
    fi
    echo "  ok: QNAME/FLAG/MAPQ match the golden under '$label'"
}

se=("$W/ref.fa" "$W/se.fq")
pe=("$W/ref.fa" "$W/r1.fq" "$W/r2.fq")
pin=(-I "1000,50")

fields "$W/se.q0" -Q 0 "${se[@]}"
fields "$W/se.q0.nbc" -Q 0 --no-band-cert "${se[@]}"
fields "$W/se.def" "${se[@]}"
fields "$W/pe.q0" -Q 0 "${pe[@]}"
fields "$W/pe.q0.nbc" -Q 0 --no-band-cert "${pe[@]}"
fields "$W/pe.def" "${pe[@]}"
fields "$W/pei.q0" -Q 0 "${pin[@]}" "${pe[@]}"
fields "$W/pei.q0.nbc" -Q 0 --no-band-cert "${pin[@]}" "${pe[@]}"
fields "$W/pei.def" "${pin[@]}" "${pe[@]}"

check "SE -Q 0" "$expect_se_q0" "$W/se.q0"
check "SE -Q 0 --no-band-cert" "$expect_se_q0" "$W/se.q0.nbc"
check "SE default -Q" "$expect_se_def" "$W/se.def"
check "PE -Q 0" "$expect_pe_q0" "$W/pe.q0"
check "PE -Q 0 --no-band-cert" "$expect_pe_q0" "$W/pe.q0.nbc"
check "PE default -Q" "$expect_pe_def" "$W/pe.def"
check "PE -I -Q 0" "$expect_pei_q0" "$W/pei.q0"
check "PE -I -Q 0 --no-band-cert" "$expect_pei_q0" "$W/pei.q0.nbc"
check "PE -I default -Q" "$expect_pei_def" "$W/pei.def"
echo "PASS: -Q 0 and default -Q MAPQs match the pre-gate goldens (SE, PE, PE with -I; -Q 0 also under --no-band-cert)"

# --- Non-vacuity. ---
# The pinned-insert legs must really take the paired-end MAPQ path: every
# record a proper pair (FLAG 0x2). int($2/2)%2 rather than and(): and() is a
# gawk builtin absent from mawk and BSD awk.
for f in "$W/pei.q0" "$W/pei.def"; do
    n_rec="$(xargs -n 3 < "$f" | wc -l | tr -d ' ')"
    # shellcheck disable=SC2016  # $AWK runs an awk program; its $N are awk fields, not shell expansions
    n_proper="$(xargs -n 3 < "$f" | "$AWK" 'int($2 / 2) % 2 == 1 { n++ } END { print n + 0 }')"
    if ((n_rec == 0 || n_proper != n_rec)); then
        echo "FAIL: only $n_proper of $n_rec records in $(basename "$f") are proper pairs — the paired-end MAPQ path is not exercised" >&2
        exit 1
    fi
done
for p in se pe pei; do
    if [[ "$(cat "$W/$p.def")" == "$(cat "$W/$p.q0")" ]]; then
        echo "FAIL: default -Q and -Q 0 emit the same MAPQs for $p — the legacy branch is not exercised" >&2
        exit 1
    fi
done
# shellcheck disable=SC2016  # $AWK runs an awk program; its $N are awk fields, not shell expansions
n_mid="$(xargs -n 3 < "$W/se.q0" | "$AWK" '$3 > 0 && $3 < 60 { n++ } END { print n + 0 }')"
n_all="$(xargs -n 3 < "$W/se.q0" | wc -l | tr -d ' ')"
if ((n_all == 0 || 2 * n_mid < n_all)); then
    echo "FAIL: only $n_mid of $n_all SE records have 0 < MAPQ < 60 at -Q 0 — seedcov is masked" >&2
    exit 1
fi
echo "PASS: pinned-insert pairs are all proper; -Q 0 differs from default -Q; $n_mid of $n_all SE MAPQs are seedcov-sensitive (0 < MAPQ < 60)"

# --- Source guard: every recount sits inside a want_seedcov `if` block. ---
# Prints "BAD <line> ..." for each recount (`seedcov = 0` / `seedcov +=`) that is
# not inside the braced block of an `if` whose condition begins with
# want_seedcov (`if (want_seedcov && ...` or `if (want_seedcov)`), then
# "SITES <n>" (the number of recount sites, counted by their `seedcov = 0` loop
# headers) and "DEPTH <d>" (the brace depth at end of input; 0 when the braces
# balance). Braces are counted after stripping comments and string / char
# literals. A want_seedcov `if` guards the block its `{` opens, on the same line
# or the next code line; the block ends at the matching `}`, so an `else`
# branch or code after it is not guarded.
# shellcheck disable=SC2016  # an awk program; its $0 and /.../ are awk, not shell
guard_prog='
    # Strip comments, string and char literals, so only code braces count.
    {
        s = $0
        if (in_cmt) {
            if (!(p = index(s, "*/"))) next
            s = substr(s, p + 2); in_cmt = 0
        }
        gsub(/"([^"\\]|\\.)*"/, "\"\"", s)
        gsub(/\047([^\047\\]|\\.)\047/, "\047\047", s)
        gsub(/\/\*([^*]|\*+[^*\/])*\*+\//, " ", s)
        sub(/\/\/.*/, "", s)
        if ((p = index(s, "/*"))) { s = substr(s, 1, p - 1); in_cmt = 1 }
    }
    s !~ /[^ \t]/ { next }
    # A want_seedcov `if` whose block did not open on its own line opens on the
    # next code line or not at all.
    pending && s !~ /^[ \t]*\{/ { pending = 0 }
    {
        g = match(s, /if \(want_seedcov( &&|\))/)
        r = match(s, /seedcov = 0|seedcov \+=/)
        for (i = 1; i <= length(s); i++) {
            if (i == g) pending = 1
            if (i == r) {
                if (substr(s, r, 11) == "seedcov = 0") n++
                if (!gsp) printf "BAD %d (not inside an if (want_seedcov ...) block): %s\n", NR, $0
            }
            c = substr(s, i, 1)
            if (c == "{") {
                depth++
                if (pending) { gstack[++gsp] = depth; pending = 0 }
            } else if (c == "}") {
                while (gsp && gstack[gsp] == depth) gsp--
                depth--
            }
        }
    }
    END { printf "SITES %d\nDEPTH %d\n", n + 0, depth }'

# Self-test: the guard must flag an ungated recount, a recount after a guarded
# block has closed, a recount in the guarded `if`'s `else`, and one under a
# negated want_seedcov, or it proves nothing.
selftest_bad() { # <label>; the case is on stdin
    local out
    out="$("$AWK" "$guard_prog")"
    if ! printf '%s\n' "$out" | grep -q '^BAD'; then
        echo "FAIL: source-guard self-test: $1 was not reported as BAD" >&2
        exit 1
    fi
}
selftest_bad "an ungated seedcov recount" << 'EOF'
    if (a->rb != H0_) {
        for (i = 0, a->seedcov = 0; i < n; ++i)
            a->seedcov += t->len;
    }
EOF
selftest_bad "a recount after its want_seedcov block closed" << 'EOF'
    if (want_seedcov) {
        a->seedcov = 0;
    }
    a->seedcov += t->len;
EOF
selftest_bad "a recount in the else of a want_seedcov if" << 'EOF'
    if (want_seedcov && a->rb != H0_) {
        x = 0;
    } else {
        for (i = 0, a->seedcov = 0; i < n; ++i)
            a->seedcov += t->len;
    }
EOF
selftest_bad "a recount under a negated want_seedcov" << 'EOF'
    if (!want_seedcov) {
        for (i = 0, a->seedcov = 0; i < n; ++i)
            a->seedcov += t->len;
    }
EOF

guard="$("$AWK" "$guard_prog" "$src")"
n_sites="$(printf '%s\n' "$guard" | sed -n 's/^SITES //p')"
depth="$(printf '%s\n' "$guard" | sed -n 's/^DEPTH //p')"
if [[ "$depth" != 0 ]]; then
    echo "FAIL: source guard: braces in src/bwamem.cpp do not balance after stripping comments and literals (end depth ${depth:-none}); the scope tracking cannot be trusted" >&2
    exit 1
fi
if printf '%s\n' "$guard" | grep -q '^BAD'; then
    echo "FAIL: seedcov recount line(s) in src/bwamem.cpp not inside a want_seedcov block:" >&2
    printf '%s\n' "$guard" | grep '^BAD' >&2
    exit 1
fi
expected_sites=11
if ((n_sites != expected_sites)); then
    echo "FAIL: found $n_sites seedcov recount sites in src/bwamem.cpp, expected $expected_sites." >&2
    echo "      Adding or removing a recount site means updating expected_sites here (and the" >&2
    echo "      site list in this header and in test/regression/README.md)." >&2
    exit 1
fi
echo "PASS: all $n_sites seedcov recount sites in src/bwamem.cpp sit inside a want_seedcov block (scope-aware source guard, self-tested)"

echo "PASS: seedcov -Q 0 / default-Q MAPQ regression"
