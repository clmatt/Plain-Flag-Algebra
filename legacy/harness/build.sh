#!/usr/bin/env bash
# Build the legacy test harness against the unmodified thesis-era headers.
# flags.h is copied up to (not including) the Mosek-dependent copositive functions.
set -euo pipefail
here="$(cd "$(dirname "$0")" && pwd)"
src="${LEGACY_SRC:-$here/..}"
out="${1:-$here/build}"
mkdir -p "$out"
cut_line=$(grep -n '^void excitingFlagAlgebra' "$src/flags.h" | cut -d: -f1)
head -n $((cut_line - 1)) "$src/flags.h" > "$out/flags_legacy.h"
cp "$src/general.h" "$src/fractions.h" "$out/"
g++ -std=c++17 -O2 -fopenmp -w \
    -I"$out" -I"$here/shim" -I/usr/include/nauty \
    -o "$out/harness" "$here/harness.cpp" -lnauty
echo "built $out/harness"
