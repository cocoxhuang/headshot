#!/usr/bin/env bash
# Imaginary quadratic fields with the unconditional C program clgrp_imag,
# written in files of GROUP discriminants each: OUT/imag_<lo>_<hi>.txt
#
# usage: run_imag_fast.sh OUTDIR XMIN XMAX [GROUP] [THREADS]
#   e.g. run_imag_fast.sh work/imag_fast 0 1000000000 10000000 4
set -euo pipefail
SRC=$(cd "$(dirname "$0")" && pwd)
OUT=$1; XMIN=$2; XMAX=$3; GROUP=${4:-10000000}; THREADS=${5:-$(nproc)}
BIN=${BIN:-$OUT/clgrp_imag}
mkdir -p "$OUT"
[ -x "$BIN" ] || gcc -O3 -march=native -fopenmp -o "$BIN" "$SRC/clgrp_imag.c" -lm
for ((lo = XMIN; lo < XMAX; lo += GROUP)); do
  hi=$((lo + GROUP)); [ "$hi" -gt "$XMAX" ] && hi=$XMAX
  f=$(printf "%s/imag_%010d_%010d.txt" "$OUT" "$lo" "$hi")
  [ -s "$f" ] && continue
  OMP_NUM_THREADS=$THREADS "$BIN" "$((lo < 3 ? 3 : lo))" "$hi" > "$f.part"
  mv "$f.part" "$f"
done
