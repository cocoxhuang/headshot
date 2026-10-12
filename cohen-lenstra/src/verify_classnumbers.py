#!/usr/bin/env python3
"""
Compare the class numbers in the class group data (h = product of the
invariant factors) with the unconditional class numbers from countforms.

usage: countforms X | verify_classnumbers.py DATAFILE [DATAFILE ...]

Both streams list the imaginary fundamental discriminants in increasing |D|;
the script checks that the discriminants agree one-to-one and that the class
numbers match, and reports the first mismatches.
"""
import lzma
import sys
from math import prod


def data_lines(paths):
    for path in paths:
        opener = lzma.open if path.endswith(".xz") else open
        with opener(path, "rt") as f:
            for line in f:
                D, cyc = line.split()
                yield int(D), prod(int(c) for c in cyc.split(","))


def main():
    paths = sys.argv[1:]
    n = bad = 0
    ref = (tuple(map(int, line.split())) for line in sys.stdin)
    for (D, h), (D2, h2) in zip(data_lines(paths), ref):
        n += 1
        if D != D2 or h != h2:
            bad += 1
            if bad <= 10:
                print(f"MISMATCH: data D={D} h={h}  countforms D={D2} h={h2}")
            if D != D2:
                print("discriminant lists out of sync; stopping")
                break
    leftover = sum(1 for _ in ref)
    print(f"compared {n} fields: {bad} mismatches; {leftover} extra countforms lines")
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
