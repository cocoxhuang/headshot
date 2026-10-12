#!/usr/bin/env bash
# Pack class group data files ("D c1,c2,..." lines, files named
# <family>_<lo>_<hi>.txt with zero-padded bounds) into xz-compressed files
# covering GROUP discriminants each:  OUTDIR/<family>_<lo>_<hi>.txt.xz
#
# usage: package.sh SRCDIR FAMILY OUTDIR [XMAX] [GROUP]
#   e.g. package.sh work/real real data/real 100000000 10000000
set -euo pipefail
SRCDIR=$1; FAM=$2; OUT=$3; XMAX=${4:-100000000}; GROUP=${5:-10000000}
mkdir -p "$OUT"
for ((lo = 0; lo < XMAX; lo += GROUP)); do
  hi=$((lo + GROUP))
  files=()
  for f in "$SRCDIR/${FAM}"_*.txt; do
    b=$(basename "$f" .txt); a=${b#*_}; f_lo=$((10#${a%_*}))
    if [ "$f_lo" -ge "$lo" ] && [ "$f_lo" -lt "$hi" ]; then files+=("$f"); fi
  done
  [ ${#files[@]} -eq 0 ] && continue
  dest=$(printf '%s/%s_%010d_%010d.txt.xz' "$OUT" "$FAM" "$lo" "$hi")
  [ -s "$dest" ] && continue                 # already packed
  printf '%s\0' "${files[@]}" | sort -z | xargs -0 cat > "$dest.tmp"
  echo "$dest"
done | xargs -P "$(nproc)" -I{} sh -c 'xz -9 -T1 -c "{}.tmp" > "{}.part" && mv "{}.part" "{}" && rm "{}.tmp"'
