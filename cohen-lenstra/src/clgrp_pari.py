#!/usr/bin/env python3
"""
Class groups of the quadratic fields Q(sqrt(D)) for every fundamental
discriminant D of a given sign with  lo <= |D| < hi.

Uses PARI's quadclassunit (Buchmann's subexponential algorithm) through cypari
(or cypari2).  As with PARI itself, the output is guaranteed under GRH: the
factor base is certified up to a Bach-type bound.  The tuning parameter c1
only changes how relations are searched for, not the certification bound, so
it affects speed but never the answer.

Output: one line per field, in increasing |D|:

    D c1,c2,...,ck

where c1, c2, ..., ck are the invariant factors of Cl(D) (c_{i+1} | c_i), so
Cl(D) = Z/c1 x ... x Z/ck and h(D) = c1*...*ck; the trivial group is "1".

Example:
    python3 clgrp_pari.py -1 1 1000000 > imag_0000000-1000000.txt
"""
import argparse
import sys
from math import isqrt

try:
    from cypari import pari
except ImportError:  # pragma: no cover
    from cypari2 import Pari
    pari = Pari()

# Relation-search factor-base constant (PARI's "c1" in quadclassunit's tech
# vector), tuned for speed at 10^6 <= |D| <= 10^8.  It does not affect
# correctness (under GRH).
C1 = {-1: 0.2, 1: 0.1}

_bnf_cyc = pari("D -> bnfinit(quadpoly(D), 1).cyc")


def squarefree_flags(lo, hi):
    """bytearray f with f[n - lo] = 1 iff n is squarefree, lo <= n < hi."""
    f = bytearray(b"\x01") * (hi - lo)
    lim = isqrt(hi) + 1
    sieve = bytearray(b"\x01") * (lim + 1)
    for p in range(2, lim + 1):
        if not sieve[p]:
            continue
        sieve[p * p :: p] = bytearray(len(range(p * p, lim + 1, p)))
        p2 = p * p
        start = -(-lo // p2) * p2
        if start < hi:
            f[start - lo :: p2] = bytearray(len(range(start, hi, p2)))
    return f


def fundamental_abs(sign, lo, hi):
    """Yield d = |D|, lo <= d < hi, in increasing order, for which D = sign*d
    is a fundamental discriminant (D = 1 is excluded)."""
    sf = squarefree_flags(lo, hi)
    lo4, hi4 = lo // 4, (hi - 1) // 4 + 1
    sf4 = squarefree_flags(lo4, hi4) if hi4 > lo4 else b""
    odd_res = 1 if sign > 0 else 3          # D = 1 mod 4  <=>  d = odd_res mod 4
    even_res = (2, 3) if sign > 0 else (1, 2)  # D/4 = 2,3 mod 4 <=> d/4 in even_res
    for d in range(lo, hi):
        if d & 1:
            if (d & 3) == odd_res and sf[d - lo] and d != 1:
                yield d
        elif not d & 3:
            n = d >> 2
            if (n & 3) in even_res and sf4[n - lo4]:
                yield d


def class_group(D, c1):
    """Invariant factors of Cl(Q(sqrt(D))) as a string like '9,3' or '1'."""
    try:
        cyc = pari.quadclassunit(D, 0, [c1])[1]
    except Exception as e:  # fall back on the general bnfinit
        print(f"quadclassunit failed at D={D}: {e}; using bnfinit", file=sys.stderr)
        cyc = _bnf_cyc(D)
    s = str(cyc)[1:-1].replace(" ", "")
    return s or "1"


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("sign", type=int, choices=(-1, 1))
    ap.add_argument("lo", type=int)
    ap.add_argument("hi", type=int)
    ap.add_argument("--c1", type=float, default=None,
                    help="override PARI's relation-search constant c1")
    args = ap.parse_args()
    if args.lo < 1 or args.hi <= args.lo:
        ap.error("need 1 <= lo < hi")
    c1 = C1[args.sign] if args.c1 is None else args.c1
    out = sys.stdout
    buf = []
    for d in fundamental_abs(args.sign, args.lo, args.hi):
        D = args.sign * d
        buf.append(f"{D} {class_group(D, c1)}\n")
        if len(buf) >= 4096:
            out.write("".join(buf))
            buf.clear()
    out.write("".join(buf))


if __name__ == "__main__":
    main()
