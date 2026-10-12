#!/usr/bin/env python3
"""
Cohen-Lenstra predictions for class groups of quadratic fields.

For an odd prime p, Cohen-Lenstra predict that the p-Sylow subgroup Cl(D)_p of
a quadratic field, ordered by |D|, is distributed like a random finite abelian
p-group G drawn with probability

    P_u(G) = c_u(p) / ( |G|^u |Aut(G)| ),      c_u(p) = eta_oo(p) / eta_u(p),

where u is the unit rank (u = 0 for imaginary, u = 1 for real quadratic
fields) and eta_k(p) = prod_{i=1..k} (1 - p^-i).  Equivalently,
E[#Sur(Cl_p, A)] = |A|^-u for every finite abelian p-group A.

For p = 2 genus theory gives Cl[2] exactly (rank = omega(D) - 1); Gerth's
extension predicts that 2*Cl, i.e. the subgroup of squares, has the same law
(with p = 2), which Fouvry-Kluners proved for the 4-rank and Smith proved for
the 2^oo part of imaginary quadratic fields.

All functions return floats; partitions are tuples in decreasing order, e.g.
(2, 1) is Z/p^2 x Z/p.

Run as a script to print the main predictions.
"""
from fractions import Fraction
from functools import lru_cache
from itertools import combinations
from math import log2, prod


def eta(p, k):
    """eta_k(p) = prod_{i=1}^k (1 - p^-i); k = None means k = infinity
    (truncated once p^-i < 2^-64, beyond double precision)."""
    if k is None:
        k = int(64 / log2(p)) + 1
    return prod(1.0 - p ** (-i) for i in range(1, k + 1))


def gauss_binom(n, k, p):
    """Number of k-dimensional subspaces of F_p^n."""
    if k < 0 or k > n:
        return 0
    num = prod(p ** (n - i) - 1 for i in range(k))
    den = prod(p ** (k - i) - 1 for i in range(k))
    return num // den


def aut_order(partition, p):
    """|Aut(G)| for G = sum_i Z/p^{e_i}, partition = (e_1, e_2, ...).
    Formula of Hillar-Rhea (Amer. Math. Monthly 114 (2007))."""
    e = sorted(x for x in partition if x > 0)       # e_1 <= ... <= e_n
    n = len(e)
    if n == 0:
        return 1
    d = [max(l for l in range(n) if e[l] == e[k]) + 1 for k in range(n)]
    c = [min(l for l in range(n) if e[l] == e[k]) + 1 for k in range(n)]
    r = 1
    for k in range(n):                               # k is 0-based here
        r *= p ** d[k] - p ** k
    for j in range(n):
        r *= (p ** e[j]) ** (n - d[j])
    for i in range(n):
        r *= (p ** (e[i] - 1)) ** (n - c[i] + 1)
    return r


def group_order(partition, p):
    return p ** sum(partition)


def prob_group(partition, p, u):
    """Cohen-Lenstra probability that Cl_p = sum Z/p^{e_i} (unit rank u)."""
    den = group_order(partition, p) ** u * aut_order(partition, p)
    return eta(p, None) / eta(p, u) * (1 / den)   # int true division: no overflow


def prob_rank(r, p, u):
    """Probability that the p-rank of Cl equals r."""
    return p ** (-r * (r + u)) * eta(p, None) / (eta(p, r) * eta(p, r + u))


def prob_divides(p, u):
    """Probability that p divides h."""
    return 1.0 - prob_rank(0, p, u)


def prob_p_power_divides(p, k, u):
    """Probability that p^k divides h (sum over groups of order >= p^k)."""
    return 1.0 - sum(prob_group(lam, p, u) for n in range(k)
                     for lam in partitions(n))


def moment_rank(k, p, u):
    """E[ |Cl[p]|^k ] = E[ p^(k * rank) ] = sum_j binom(k, j)_p p^(-j u)."""
    return float(sum(Fraction(gauss_binom(k, j, p)) / Fraction(p) ** (j * u)
                     for j in range(k + 1)))


def mean_surj(order_A, u):
    """E[#Sur(Cl_p, A)] = |A|^-u."""
    return order_A ** (-u)


@lru_cache(maxsize=None)
def partitions(n, maxpart=None):
    """All partitions of n as decreasing tuples."""
    if maxpart is None:
        maxpart = n
    if n == 0:
        return ((),)
    out = []
    for first in range(min(n, maxpart), 0, -1):
        for rest in partitions(n - first, first):
            out.append((first,) + rest)
    return tuple(out)


def prob_cyclic(p, u):
    """Probability that Cl_p is cyclic (p-rank <= 1)."""
    return prob_rank(0, p, u) + prob_rank(1, p, u)


def odd_primes(limit):
    sieve = bytearray(b"\x01") * (limit + 1)
    out = []
    for q in range(2, limit + 1):
        if sieve[q]:
            if q > 2:
                out.append(q)
            sieve[q * q :: q] = bytearray(len(range(q * q, limit + 1, q)))
    return out


def prob_odd_part_cyclic(u, plimit=100000):
    """Probability that the odd part of Cl is cyclic: prod_{p odd} P(rk_p<=1).
    The tail beyond plimit is < 1e-14."""
    return prod(prob_cyclic(q, u) for q in odd_primes(plimit))


def prob_odd_part_trivial_real(plimit=10**6):
    """u = 1: probability that the odd part of Cl is trivial,
    prod_{p odd} prod_{k>=2} (1 - p^-k)  (~ 0.754458, the famous h=1 density
    for real quadratic fields of prime discriminant).  The primes beyond
    plimit contribute a factor exp(-sum_{q > plimit} q^-2) ~ 1 - 1/(L log L),
    which we include."""
    from math import exp, log
    tail = 1.0 / (plimit * log(plimit))
    return prod(prob_rank(0, q, 1) for q in odd_primes(plimit)) * exp(-tail)


def describe(partition):
    if not partition:
        return "1"
    return " x ".join(f"Z/p^{e}" if e > 1 else "Z/p" for e in partition)


def main():
    print("Cohen-Lenstra predictions (u=0 imaginary, u=1 real)\n")
    print(f"{'p':>3} {'P(p|h) imag':>12} {'P(p|h) real':>12} "
          f"{'E|Cl[p]| imag':>14} {'E|Cl[p]| real':>14}")
    for p in (3, 5, 7, 11, 13, 17, 19, 23):
        print(f"{p:>3} {prob_divides(p, 0):12.6f} {prob_divides(p, 1):12.6f} "
              f"{moment_rank(1, p, 0):14.6f} {moment_rank(1, p, 1):14.6f}")
    for p in (3, 5, 7):
        print(f"\np = {p}: p-rank distribution")
        for r in range(5):
            print(f"  rank {r}: imag {prob_rank(r, p, 0):.8f}   "
                  f"real {prob_rank(r, p, 1):.8f}")
        print(f"  p-Sylow distribution (groups of order <= p^4)")
        for n in range(5):
            for lam in partitions(n):
                print(f"  {describe(lam):<28} imag {prob_group(lam, p, 0):.8f}"
                      f"   real {prob_group(lam, p, 1):.8f}")
    print("\nGerth (p = 2, applied to 2*Cl): 4-rank distribution")
    for r in range(5):
        print(f"  4-rank {r}: imag {prob_rank(r, 2, 0):.8f}   "
              f"real {prob_rank(r, 2, 1):.8f}")
    print(f"\nP(odd part of Cl cyclic): imag {prob_odd_part_cyclic(0):.6f}"
          f"   real {prob_odd_part_cyclic(1):.6f}")
    print(f"P(odd part trivial), real: {prob_odd_part_trivial_real():.6f}")


if __name__ == "__main__":
    main()
