#!/usr/bin/env bash
# test/fast_preset_test.sh
#
# Asserts that `bwa-mem3 mem --fast` resolves the characterized speed levers
# (-m 10 -y 0 --min-ext-len 30 --smem-dedup
# --max-extend-chains 20 --extend-mate-concordant --rescue-kmer=6, plus -s 2 and a
# lower --max-extend-chains 10 under --meth), that explicit user flags override the
# preset, and that the default path is untouched when --fast is absent.
# The contained-seed extension skip is the byte-identical DEFAULT, not a --fast
# lever, so `--skip-contained-ext` must NOT appear on either audit line (the line
# records only the output-changing levers --fast turns on); --max-extend-chains
# applies under --meth too but at a lower cap of 10 (non-meth uses 20).
# --extend-mate-concordant recovers the chain-cap pairing regression and is now
# enabled for both non-meth and --meth --fast (fg-labs/bwa-mem3#202), so it must
# be present on both audit lines. The deprecated `--skip-contained-ext` and its
# opt-out `--keep-contained-ext` must both still parse (with and without --fast).
#
# The assertion surface is the audit line main_mem prints to stderr when
# --fast is active: it reports the *resolved* mem_opt_t values, so an explicit
# override is visible in the line itself.
#
# Usage: test/fast_preset_test.sh <bwa-mem3-binary> <fixtures-dir>
set -euo pipefail

[[ $# -eq 2 ]] || {
    echo "usage: $0 <bwa-mem3-binary> <fixtures-dir>" >&2
    exit 2
}
bin="$1"
fixtures="$2"
src_ref="$fixtures/phix.fa"
reads="$fixtures/reads.fa"

[[ -x "$bin" ]] || {
    echo "FAIL: binary not executable: $bin" >&2
    exit 1
}
[[ -s "$src_ref" ]] || {
    echo "FAIL: phix.fa missing: $src_ref" >&2
    exit 1
}
[[ -s "$reads" ]] || {
    echo "FAIL: reads.fa missing: $reads" >&2
    exit 1
}

err="$(mktemp)"
mdir="$(mktemp -d)"
trap 'rm -f "$err"; rm -rf "$mdir"' EXIT

# Index a private copy of phiX in temp space so the test never writes derived
# index files into the (possibly read-only or shared) fixtures tree. The meth
# path below does the same; keep both hermetic. The copy is always fresh, so
# index unconditionally.
ref="$mdir/phix.fa"
cp "$src_ref" "$ref"
"$bin" index "$ref" > /dev/null 2>&1 || {
    echo "FAIL: index phix.fa" >&2
    exit 1
}

# Run mem with extra args; echo the resolved --fast audit line (empty if none).
fast_line() {
    "$bin" mem "$@" "$ref" "$reads" > /dev/null 2> "$err" \
        || {
            echo "FAIL: mem exited nonzero (args: $*)" >&2
            cat "$err" >&2
            exit 1
        }
    grep -E '^\[M::main_mem\] --fast:' "$err" || true
}

# 1. Bundle correctness (non-meth).
line="$(fast_line --fast)"
[[ "$line" == *"-m 10"* && "$line" == *"-y 0"* &&
    "$line" == *"--min-ext-len 30"* && "$line" == *"--smem-dedup"* &&
    "$line" == *"--max-extend-chains 20"* &&
    "$line" == *"--adaptive-band"* &&
    "$line" == *"--extend-mate-concordant"* &&
    "$line" == *"--extend-tie-frac 0.95"* &&
    "$line" == *"--extend-tie-floor 1"* &&
    "$line" == *"--extend-csub"* &&
    "$line" == *"--rescue-kmer=6"* &&
    "$line" == *"alnreg-sort=fast"* ]] \
    || {
        echo "FAIL: --fast bundle wrong: '$line'" >&2
        exit 1
    }
[[ "$line" != *"-s "* ]] \
    || {
        echo "FAIL: non-meth --fast must not set -s: '$line'" >&2
        exit 1
    }
# The contained-seed skip is the default, not a --fast lever: it must not be
# reported as one (a stale audit line would misdescribe what --fast changes).
[[ "$line" != *"--skip-contained-ext"* ]] \
    || {
        echo "FAIL: --skip-contained-ext is the default, not a --fast lever; must not be on the audit line: '$line'" >&2
        exit 1
    }
echo "OK:   --fast bundle resolves -m 10 -y 0 --min-ext-len 30 --smem-dedup --max-extend-chains 20 --adaptive-band --extend-mate-concordant --extend-tie-frac 0.95 --extend-tie-floor 1 --extend-csub --rescue-kmer=6 alnreg-sort=fast (no -s, no --skip-contained-ext)"

# 1b. The gate levers honor an explicit user value: --extend-tie-frac 0 disables the
#     competitiveness gate even under --fast (opt0 wins), and the audit line reflects it.
#     --extend-csub has no opt-out and stays on; --extend-tie-floor is untouched.
line="$(fast_line --fast --extend-tie-frac 0)"
[[ "$line" == *"--extend-tie-frac 0 "* &&
    "$line" == *"--extend-tie-floor 1"* &&
    "$line" == *"--extend-csub"* ]] \
    || {
        echo "FAIL: explicit --extend-tie-frac 0 should override the preset in the audit line: '$line'" >&2
        exit 1
    }
echo "OK:   --fast folds in --extend-tie-frac 0.95 / --extend-tie-floor 1 / --extend-csub; explicit --extend-tie-frac 0 overrides the gate"

# 2. Override precedence: an explicit flag wins *only* for its own field; the
#    rest of the preset (including the unconditional --smem-dedup) must survive.
line="$(fast_line --fast -m 30)"
[[ "$line" == *"-m 30"* && "$line" == *"-y 0"* &&
    "$line" == *"--min-ext-len 30"* && "$line" == *"--smem-dedup"* &&
    "$line" == *"--max-extend-chains 20"* &&
    "$line" == *"--adaptive-band"* && "$line" == *"--extend-mate-concordant"* &&
    "$line" == *"alnreg-sort=fast"* ]] \
    || {
        echo "FAIL: explicit -m 30 should only override -m: '$line'" >&2
        exit 1
    }
line="$(fast_line --fast -y 5)"
[[ "$line" == *"-m 10"* && "$line" == *"-y 5"* &&
    "$line" == *"--min-ext-len 30"* && "$line" == *"--smem-dedup"* &&
    "$line" == *"--max-extend-chains 20"* &&
    "$line" == *"--adaptive-band"* && "$line" == *"--extend-mate-concordant"* &&
    "$line" == *"alnreg-sort=fast"* ]] \
    || {
        echo "FAIL: explicit -y 5 should only override -y: '$line'" >&2
        exit 1
    }
line="$(fast_line --fast --min-ext-len 45)"
[[ "$line" == *"-m 10"* && "$line" == *"-y 0"* &&
    "$line" == *"--min-ext-len 45"* && "$line" == *"--smem-dedup"* &&
    "$line" == *"--max-extend-chains 20"* &&
    "$line" == *"--adaptive-band"* && "$line" == *"--extend-mate-concordant"* &&
    "$line" == *"alnreg-sort=fast"* ]] \
    || {
        echo "FAIL: explicit --min-ext-len 45 should only override min-ext-len: '$line'" >&2
        exit 1
    }
# --max-extend-chains is overridable under --fast (respects an explicit value).
line="$(fast_line --fast --max-extend-chains 8)"
[[ "$line" == *"-m 10"* && "$line" == *"-y 0"* &&
    "$line" == *"--min-ext-len 30"* && "$line" == *"--smem-dedup"* &&
    "$line" == *"--max-extend-chains 8"* &&
    "$line" == *"--adaptive-band"* && "$line" == *"--extend-mate-concordant"* &&
    "$line" == *"alnreg-sort=fast"* ]] \
    || {
        echo "FAIL: explicit --max-extend-chains 8 should only override that lever: '$line'" >&2
        exit 1
    }
echo "OK:   explicit -m/-y/--min-ext-len/--max-extend-chains override only their field; rest of preset (--adaptive-band, --extend-mate-concordant, alnreg-sort=fast, ...) survives"

# --rescue-kmer is overridable both ways under --fast. The `=0` case is the one
# that matters: 0 means "off", so a parser that treated an explicitly-set 0 as
# "unset" would silently re-enable the preset's K=6 -- and since --rescue-kmer
# changes MAPQ on rescued reads, that is an invisible output change. The value
# must be attached with `=`: the option takes an OPTIONAL argument, so a
# space-separated `--rescue-kmer 0` would leave 0 as a positional instead.
line="$(fast_line --fast --rescue-kmer=0)"
[[ "$line" == *"-m 10"* && "$line" == *"--smem-dedup"* &&
    "$line" == *"--rescue-kmer=0"* &&
    "$line" == *"alnreg-sort=fast"* ]] \
    || {
        echo "FAIL: explicit --rescue-kmer=0 should opt out of just that lever: '$line'" >&2
        exit 1
    }
line="$(fast_line --fast --rescue-kmer=11)"
[[ "$line" == *"--rescue-kmer=11"* && "$line" == *"--max-extend-chains 20"* ]] \
    || {
        echo "FAIL: explicit --rescue-kmer=11 should override the preset K: '$line'" >&2
        exit 1
    }
echo "OK:   explicit --rescue-kmer=0/=11 overrides the preset K without disturbing the rest"

# 3. Default contract: no --fast => no audit line at all.
"$bin" mem "$ref" "$reads" > /dev/null 2> "$err" || {
    echo "FAIL: plain mem nonzero" >&2
    exit 1
}
! grep -qE '^\[M::main_mem\] --fast:' "$err" \
    || {
        echo "FAIL: audit line present without --fast" >&2
        exit 1
    }
echo "OK:   no --fast => no audit line (default path untouched)"

# 3b. Contained-seed skip CLI surface. The skip is the byte-identical default;
#     `--keep-contained-ext` is its opt-out and `--skip-contained-ext` is the
#     deprecated spelling of the default. Both must still parse (alone and with
#     --fast), the deprecated one must say so on stderr, and the opt-out must not.
"$bin" mem --keep-contained-ext "$ref" "$reads" > /dev/null 2> "$err" || {
    echo "FAIL: mem --keep-contained-ext nonzero" >&2
    cat "$err" >&2
    exit 1
}
! grep -q 'deprecated' "$err" \
    || {
        echo "FAIL: --keep-contained-ext must not print a deprecation notice" >&2
        exit 1
    }
"$bin" mem --skip-contained-ext "$ref" "$reads" > /dev/null 2> "$err" || {
    echo "FAIL: mem --skip-contained-ext (deprecated) must still parse" >&2
    cat "$err" >&2
    exit 1
}
grep -qE '^\[W::main_mem\] --skip-contained-ext is deprecated' "$err" \
    || {
        echo "FAIL: --skip-contained-ext must print a deprecation notice" >&2
        cat "$err" >&2
        exit 1
    }
line="$(fast_line --fast --keep-contained-ext)"
[[ "$line" == *"--smem-dedup"* && "$line" != *"--skip-contained-ext"* ]] \
    || {
        echo "FAIL: --fast --keep-contained-ext should keep the preset and not report the skip: '$line'" >&2
        exit 1
    }
echo "OK:   --keep-contained-ext parses silently; --skip-contained-ext parses with a deprecation notice; both coexist with --fast"

# 4. Meth path: --fast --meth additionally sets -s 2. Build a meth index on a
#    copy of phiX; if meth indexing is unavailable in this build, SKIP.
cp "$ref" "$mdir/ref.fa"
if "$bin" index --meth "$mdir/ref.fa" > /dev/null 2>&1; then
    "$bin" mem --meth --fast -t 1 "$mdir/ref.fa" "$reads" > /dev/null 2> "$err" \
        || {
            echo "FAIL: mem --meth --fast nonzero" >&2
            cat "$err" >&2
            exit 1
        }
    line="$(grep -E '^\[M::main_mem\] --fast:' "$err" || true)"
    [[ "$line" == *"-s 2"* ]] \
        || {
            echo "FAIL: --fast --meth should resolve -s 2: '$line'" >&2
            exit 1
        }
    [[ "$line" != *"--skip-contained-ext"* ]] \
        || {
            echo "FAIL: --skip-contained-ext is the default, not a --fast lever; must not be on the meth audit line: '$line'" >&2
            exit 1
        }
    [[ "$line" == *"--max-extend-chains 10"* ]] \
        || {
            echo "FAIL: --max-extend-chains applies under --meth (cap 10); must be in audit line: '$line'" >&2
            exit 1
        }
    [[ "$line" == *"--extend-mate-concordant"* ]] \
        || {
            echo "FAIL: --fast --meth must enable --extend-mate-concordant: '$line'" >&2
            exit 1
        }
    # The score-gated chain-extension cap is folded into --fast under --meth too
    # (applies on top of the meth --max-extend-chains 10; validated placement-neutral).
    [[ "$line" == *"--extend-tie-frac 0.95"* &&
        "$line" == *"--extend-tie-floor 1"* &&
        "$line" == *"--extend-csub"* ]] \
        || {
            echo "FAIL: --fast --meth must fold in --extend-tie-frac 0.95/--extend-tie-floor 1/--extend-csub: '$line'" >&2
            exit 1
        }
    # --rescue-kmer applies under --meth as well (the anchor scan collapses to the
    # pair's bisulfite alphabet there), so the meth audit line must report it too.
    [[ "$line" == *"--rescue-kmer=6"* ]] \
        || {
            echo "FAIL: --fast --meth must enable --rescue-kmer=6: '$line'" >&2
            exit 1
        }
    # The dedup-sort lever is not meth-gated; it must be reported under --meth too.
    [[ "$line" == *"alnreg-sort=fast"* ]] \
        || {
            echo "FAIL: --fast --meth must enable alnreg-sort=fast: '$line'" >&2
            exit 1
        }
    echo "OK:   --fast --meth additionally sets -s 2, --extend-mate-concordant, the extend-tie gate and alnreg-sort=fast (--max-extend-chains lowered to 10, no --skip-contained-ext)"
    # Explicit -s wins even under --meth (src/fastmap.cpp: -s 2 is gated on !opt0.split_width).
    "$bin" mem --meth --fast -s 7 -t 1 "$mdir/ref.fa" "$reads" > /dev/null 2> "$err" \
        || {
            echo "FAIL: mem --meth --fast -s 7 nonzero" >&2
            cat "$err" >&2
            exit 1
        }
    line="$(grep -E '^\[M::main_mem\] --fast:' "$err" || true)"
    [[ "$line" == *"-s 7"* ]] \
        || {
            echo "FAIL: explicit -s 7 should win under --meth: '$line'" >&2
            exit 1
        }
    echo "OK:   explicit -s 7 overrides --fast --meth's -s 2"
else
    echo "SKIP: index --meth unavailable; meth -s 2 case not checked" >&2
fi

echo "PASS: fast_preset_test"
