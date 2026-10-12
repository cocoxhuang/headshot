#!/usr/bin/env bash
# Compute class groups of all quadratic fields with |D| < XMAX (both signs),
# in chunks of CHUNK discriminants, JOBS chunks at a time.  Re-running the
# script resumes: finished chunks (renamed from *.part on success) are skipped.
#
# usage: run_all.sh [OUTDIR] [XMAX] [CHUNK] [JOBS] [SIGNS] [XMIN]
#   defaults: work 100000000 1000000 $(nproc) "-1 1" 0
#   e.g. continue the real fields from 1.7e8 to 1e9:
#        run_all.sh work 1000000000 1000000 4 "1" 170000000
set -euo pipefail
SRC=$(cd "$(dirname "$0")" && pwd)
OUT=${1:-work}; XMAX=${2:-100000000}; CHUNK=${3:-1000000}
JOBS=${4:-$(nproc)}; SIGNS=${5:-"-1 1"}; XMIN=${6:-0}
mkdir -p "$OUT/imag" "$OUT/real"
for sign in $SIGNS; do
  name=$([ "$sign" = -1 ] && echo imag || echo real)
  for ((lo = XMIN; lo < XMAX; lo += CHUNK)); do
    hi=$((lo + CHUNK)); [ "$hi" -gt "$XMAX" ] && hi=$XMAX
    f=$(printf "%s/%s/%s_%010d_%010d.txt" "$OUT" "$name" "$name" "$lo" "$hi")
    [ -s "$f" ] || echo "$sign $((lo == 0 ? 1 : lo)) $hi $f"
  done
done | xargs -P "$JOBS" -L 1 sh -c \
  'python3 "$0" "$1" "$2" "$3" > "$4.part" && mv "$4.part" "$4"' "$SRC/clgrp_pari.py"
