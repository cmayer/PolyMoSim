#!/bin/sh
# SPDX-License-Identifier: GPL-3.0-only
# Tests for PolyMoSim.
#
# Usage:  sh tests/run_tests.sh [path/to/PolyMoSim]      (or: make test)
#
# The simulations use fixed seeds, so the output has to be identical to the reference files.

# awk must read numbers with a decimal point. In e.g. a German locale mawk reads 0.5 as 0.
LC_ALL=C
export LC_ALL

BIN=${1:-./PolyMoSim}
TESTS=$(cd "$(dirname "$0")" && pwd)
ROOT=$(dirname "$TESTS")
EXAMPLE=$ROOT/Examples/Example-for-sitepattern

case $BIN in
    */*) BIN=$(cd "$(dirname "$BIN")" && pwd)/$(basename "$BIN") ;;
esac

TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT

passed=0
failed=0

pass() { passed=$((passed + 1)); echo "ok:     $1"; }
fail() { failed=$((failed + 1)); echo "FAILED: $1"; }

# check NAME EXPECTED_FILE ARGS... : runs PolyMoSim with ARGS and compares its output with EXPECTED_FILE.
check() {
    name=$1
    expected=$2
    shift 2
    "$BIN" -l "$TMP/log.txt" "$@" > "$TMP/out.txt" 2> "$TMP/err.txt"
    if [ -s "$TMP/out.txt" ] && cmp -s "$TMP/out.txt" "$expected"; then
        pass "$name"
    else
        fail "$name"
        diff "$expected" "$TMP/out.txt" | head -n 10
        tail -n 5 "$TMP/err.txt"
    fi
}

# 1. --version prints the version defined in PolyMoSim.h.
version=$(sed -n 's/^#define VERSION *"\(.*\)".*/\1/p' "$ROOT/PolyMoSim.h")
if "$BIN" --version 2>/dev/null | grep -q "version: $version\$"; then
    pass "--version ($version)"
else
    fail "--version: expected $version"
fi

# 2. The example for site patterns: the output files in the example folder.
#    Only fasta and site_pattern_freq_relative: the other files there were made with an older version.
cd "$EXAMPLE" || exit 1
out=output-that-will-be-generated
check "example for site patterns, fasta" "$out/sim-sequences.fas" \
      -s 127 -m Modelfile_basic.txt -t tree1.txt -n 1 -f fasta
check "example for site patterns, site_pattern_freq_relative" "$out/site_pattern_freq_rel.txt" \
      -s 127 -m Modelfile_basic.txt -t tree1.txt -n 1 -f site_pattern_freq_relative

# 3. Nucleotide models JC, HKY and GTR with gamma and invariant sites, and a model change on a branch.
cd "$TESTS" || exit 1
check "nucleotide models" expected/dna.fas -s 11 -m dna-models.txt -t dna-tree.txt -n 1 -f fasta
check "nucleotide site patterns, all 256 patterns" expected/dna-site-patterns-fill.txt \
      -s 11 -m dna-models.txt -t dna-tree.txt -n 1 -f site_pattern_freq_absolute_fill

# 4. All fixed rate amino acid models, and LG and WAG with gamma and invariant sites.
check "amino acid models" expected/aa.fas -s 11 -m aa-models.txt -t aa-tree.txt -n 1 -f fasta
check "amino acid site patterns" expected/aa-site-patterns.txt \
      -s 11 -m aa-models.txt -t aa-tree.txt -n 1 -f site_pattern_freq_absolute

# 5. Amino acid site pattern frequencies with all 20^4 = 160000 patterns, as used by PhyloAI.
"$BIN" -l "$TMP/log.txt" -s 11 -m aa-models.txt -t aa-tree.txt -n 1 \
       -f site_pattern_freq_relative_fill > "$TMP/out.txt" 2> /dev/null
if awk '{ n++; s += $2 } END { exit !(n == 160000 && s > 0.99 && s < 1.01) }' "$TMP/out.txt"; then
    pass "amino acid site patterns, all 160000 patterns, sum 1"
else
    fail "amino acid site patterns, all 160000 patterns, sum 1"
fi

echo "$passed passed, $failed failed"
[ "$failed" -eq 0 ]
