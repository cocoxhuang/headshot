#!/usr/bin/env python3
"""
Empirical class group statistics of quadratic fields vs. the Cohen-Lenstra
heuristics (and Gerth's extension to p = 2, and the Davenport-Heilbronn /
Bhargava-Shankar-Tsimerman / Taniguchi-Thorne theorems for p = 3).

Reads class group data files ("D c1,c2,...", one field per line, increasing
|D|; plain text or .xz), computes invariants of every field and accumulates
counts in windows of |D| (log-spaced, 8 per decade, fixed edges 10^(2+k/8)).

usage:
  analyze.py all    FAMILY OUTDIR FILE [FILE ...]   count + report in one go
  analyze.py count  FAMILY WINDOWS.csv FILE [FILE ...]
  analyze.py report FAMILY OUTDIR WINDOWS.csv [WINDOWS.csv ...]

"count" can be run on disjoint sets of files in parallel; "report" adds up
their window counts and writes to OUTDIR:

  <family>_windows.csv      raw counts per window (long format: lo,hi,stat,count)
  <family>_cumulative.csv   empirical frequency of each statistic for |D| < X at
                            each window edge X (and in the window itself), next
                            to the Cohen-Lenstra prediction
  <family>_summary.md       table at the largest X

family is "imag" (D < 0, unit rank u = 0) or "real" (D > 0, u = 1).
"""
import csv
import lzma
import math
import sys
from collections import defaultdict

import numpy as np

sys.path.insert(0, __import__("os").path.dirname(__file__))
import cl_predictions as cl  # noqa: E402

ODD_PRIMES = (3, 5, 7, 11, 13, 17, 19)
TYPE_PRIMES = (3, 5, 7)          # p-Sylow isomorphism types are tabulated
MAX_TYPE_ORDER = 4               # ... for groups of order <= p^4
PAIRS = ((3, 5), (3, 7), (5, 7))  # independence of different primes
LOCAL_PRIMES = (3, 5)            # statistics split by how p behaves in K
EDGES_PER_DECADE = 8


def window_edges(xmax=10 ** 11):
    """Log-spaced window edges 1, 10, 100, round(10^(2 + k/8)), ... >= xmax.
    The edges do not depend on the data, so counts from different runs can be
    added up."""
    edges = [1, 10, 100]
    k = 1
    while edges[-1] < xmax:
        edges.append(int(round(10 ** (2 + k / EDGES_PER_DECADE))))
        k += 1
    return np.array(edges, dtype=np.int64)


def open_any(path):
    return lzma.open(path, "rt") if path.endswith(".xz") else open(path)


def read_chunks(paths, chunk_lines=500000):
    """Yield (D, C) numpy arrays: D (n,), C (n, k) invariant factors padded
    with 1s, in file order."""
    for path in paths:
        with open_any(path) as f:
            while True:
                lines = f.readlines(chunk_lines * 16)
                if not lines:
                    break
                n = len(lines)
                D = np.empty(n, dtype=np.int64)
                rows = []
                for i, line in enumerate(lines):
                    a, b = line.split()
                    D[i] = int(a)
                    rows.append(b.split(","))
                k = max(len(r) for r in rows)
                C = np.ones((n, k), dtype=np.int64)
                for i, r in enumerate(rows):
                    C[i, : len(r)] = [int(x) for x in r]
                yield D, C


def valuation_matrix(C, p):
    V = np.zeros(C.shape, dtype=np.int64)
    X = C.copy()
    while True:
        m = X % p == 0
        if not m.any():
            return V
        V += m
        X = np.where(m, X // p, X)


def partition_key(row):
    return ",".join(str(int(e)) for e in row if e > 0) or "0"


def partition_keys(V):
    """Row-wise partition keys ("2,1", "0" = trivial, "larger" if the group
    has order > p^MAX_TYPE_ORDER) of a valuation matrix V."""
    uniq, inv = np.unique(V, axis=0, return_inverse=True)
    names = np.array([partition_key(row) if row.sum() <= MAX_TYPE_ORDER
                      else "larger" for row in uniq])
    return names[inv.ravel()]


def prime_flags(lo, hi):
    """Boolean array: is n prime, for lo <= n < hi."""
    n = hi - lo
    flags = np.ones(n, dtype=bool)
    lim = math.isqrt(hi) + 1
    small = np.ones(lim + 1, dtype=bool)
    small[:2] = False
    for q in range(2, math.isqrt(lim) + 1):
        if small[q]:
            small[q * q :: q] = False
    for q in np.nonzero(small)[0]:
        q = int(q)
        start = max(q * q, -(-lo // q) * q)
        flags[start - lo :: q] = False
    for x in (0, 1):
        if lo <= x < hi:
            flags[x - lo] = False
    return flags


class Accumulator:
    """counts[stat][window] += ..."""

    def __init__(self, edges):
        self.edges = edges
        self.counts = defaultdict(lambda: np.zeros(len(edges) - 1, dtype=np.float64))
        self.dmax = 0                 # largest |D| seen
        self.range_end = 0            # end of the computed range (from file names)

    def add(self, stat, win, weights=None, mask=None):
        if mask is not None:
            win = win[mask]
            if weights is not None:
                weights = weights[mask]
        c = np.bincount(win, weights=weights, minlength=len(self.edges) - 1)
        self.counts[stat] += c

    def add_keys(self, prefix, win, keys):
        """Count occurrences of string keys per window."""
        for key in np.unique(keys):
            self.add(f"{prefix}{key}", win, mask=(keys == key))


def range_end_of(path):
    """Data files are named <family>_<lo>_<hi>.txt[.xz]: return hi, or 0."""
    import os
    import re
    m = re.search(r"_(\d+)_(\d+)\.txt", os.path.basename(path))
    return int(m.group(2)) if m else 0


def process(family, paths, edges):
    acc = Accumulator(edges)
    acc.range_end = max(range_end_of(p) for p in paths)
    sign = -1 if family == "imag" else 1
    for D, C in read_chunks(paths):
        d = sign * D
        assert (d > 0).all(), "wrong family for these files"
        win = np.searchsorted(edges, d, side="right") - 1
        acc.dmax = max(acc.dmax, int(d.max()))
        h = np.prod(C, axis=1)
        acc.add("N", win)

        odd = C.copy()
        V2 = valuation_matrix(C, 2)
        odd //= 2 ** V2
        # odd part cyclic <=> second invariant factor has trivial odd part
        oddcyc = odd[:, 1] == 1 if odd.shape[1] > 1 else np.ones(len(d), bool)
        oddtriv = odd[:, 0] == 1
        acc.add("odd_part_cyclic", win, mask=oddcyc)
        acc.add("odd_part_trivial", win, mask=oddtriv)
        acc.add("h_eq_1", win, mask=(h == 1))

        # prime discriminants: imaginary -4, -8, -q (q = 3 mod 4); real 8, q
        lo, hi = int(d.min()), int(d.max()) + 1
        isprime = prime_flags(lo, hi)[d - lo]
        primedisc = isprime | (d == 8) | ((d == 4) & (sign < 0))
        acc.add("prime_N", win, mask=primedisc)
        acc.add("prime_h_eq_1", win, mask=primedisc & (h == 1))
        acc.add("prime_odd_part_trivial", win, mask=primedisc & oddtriv)

        # odd primes: divisibility, rank, moments, Sylow types
        divs = {}
        for p in ODD_PRIMES:
            V = valuation_matrix(C, p)
            rank = (V > 0).sum(axis=1)
            divs[p] = rank > 0
            acc.add(f"p{p}_div", win, mask=divs[p])
            acc.add(f"p{p}_div_p2", win, mask=V.sum(axis=1) >= 2)
            for r in range(3):
                acc.add(f"p{p}_rank_{r}", win, mask=(rank == r))
            acc.add(f"p{p}_rank_3+", win, mask=(rank >= 3))
            acc.add(f"p{p}_sum_p^rank", win, weights=(float(p) ** rank))
            acc.add(f"p{p}_sum_p^2rank", win, weights=(float(p) ** (2 * rank)))
            if p in TYPE_PRIMES:
                acc.add_keys(f"p{p}_type_", win, partition_keys(V))
        for p, q in PAIRS:
            acc.add(f"p{p}p{q}_both_div", win, mask=divs[p] & divs[q])

        # local conditions: how p behaves in K = Q(sqrt D) (Legendre symbol (D/p))
        for p in LOCAL_PRIMES:
            Dm = D % p
            squares = np.unique((np.arange(1, p) ** 2) % p)
            kinds = {"ram": Dm == 0, "split": np.isin(Dm, squares),
                     "inert": (Dm != 0) & ~np.isin(Dm, squares)}
            V = valuation_matrix(C, p)
            rank = (V > 0).sum(axis=1)
            for kind, m in kinds.items():
                tag = f"@{p}{kind}"
                acc.add("N" + tag, win, mask=m)
                acc.add(f"p{p}_div" + tag, win, mask=m & (rank > 0))
                acc.add(f"p{p}_sum_p^rank" + tag, win, weights=float(p) ** rank, mask=m)

        # 2-part: 2-rank t, 4-rank, 8-rank, and type of 2*Cl[2^oo] (Gerth)
        rank2 = (V2 >= 1).sum(axis=1)
        rank4 = (V2 >= 2).sum(axis=1)
        rank8 = (V2 >= 3).sum(axis=1)
        for r in range(4):
            acc.add(f"two_rank4_{r}", win, mask=(rank4 == r))
            acc.add(f"two_rank8_{r}", win, mask=(rank8 == r))
        acc.add("two_rank4_4+", win, mask=(rank4 >= 4))
        acc.add("two_rank8_4+", win, mask=(rank8 >= 4))
        acc.add("two_sum_2^rank4", win, weights=2.0 ** rank4)
        for t in range(int(rank2.max()) + 1):
            mt = rank2 == t
            if not mt.any():
                continue
            acc.add(f"two_rank2_{t}", win, mask=mt)
            for r in range(t + 1):
                acc.add(f"two_rank2_{t}_rank4_{r}", win, mask=mt & (rank4 == r))
        V2sq = np.maximum(V2 - 1, 0)        # 2-adic type of 2*Cl
        acc.add_keys("two_2Cl_type_", win, partition_keys(V2sq))
    return acc


# --------------------------------------------------------------------------
# predictions

def corank_dist_square(t, r, q=2):
    """Probability that a uniformly random t x t matrix over F_q has corank r."""
    k = t - r
    if k < 0:
        return 0.0
    num = 1
    for i in range(k):
        num *= (q ** t - q ** i) ** 2
    den = 1
    for i in range(k):
        den *= q ** k - q ** i
    return num / den / q ** (t * t)


def prediction(stat, u):
    """Cohen-Lenstra (Gerth for p = 2) prediction for a statistic, or None.
    A subfamily tag "@..." (local conditions) does not change the prediction."""
    stat = stat.split("@")[0]
    parts = stat.split("_")
    if stat.startswith("p") and parts[0][1:].isdigit():
        p = int(parts[0][1:])
        rest = "_".join(parts[1:])
        if rest == "div":
            return cl.prob_divides(p, u)
        if rest == "div_p2":
            return cl.prob_p_power_divides(p, 2, u)
        if rest.startswith("rank_"):
            r = parts[2]
            if r == "3+":
                return 1 - sum(cl.prob_rank(s, p, u) for s in range(3))
            return cl.prob_rank(int(r), p, u)
        if rest == "sum_p^rank":
            return cl.moment_rank(1, p, u)
        if rest == "sum_p^2rank":
            return cl.moment_rank(2, p, u)
        if rest.startswith("type_"):
            key = parts[2]
            if key == "larger":
                return 1 - sum(cl.prob_group(lam, p, u)
                               for n in range(MAX_TYPE_ORDER + 1)
                               for lam in cl.partitions(n))
            lam = () if key == "0" else tuple(int(x) for x in key.split(","))
            return cl.prob_group(lam, p, u)
    if stat.startswith("p") and "both_div" in stat:
        p, q = (int(x) for x in stat.split("_")[0][1:].split("p"))
        return cl.prob_divides(p, u) * cl.prob_divides(q, u)
    if stat == "odd_part_cyclic":
        return cl.prob_odd_part_cyclic(u)
    if stat in ("odd_part_trivial", "prime_odd_part_trivial", "prime_h_eq_1"):
        if u == 1:
            return cl.prob_odd_part_trivial_real()
        return 0.0
    if u == 0 and stat.startswith("two_"):
        # Gerth: 2*Cl[2^oo] is CL-distributed with p = 2 (imaginary fields;
        # Fouvry-Kluners for the 4-rank, Smith for 2^oo)
        if stat.startswith("two_rank4_"):
            r = parts[2]
            if r == "4+":
                return 1 - sum(cl.prob_rank(s, 2, 0) for s in range(4))
            return cl.prob_rank(int(r), 2, 0)
        if stat == "two_sum_2^rank4":
            return cl.moment_rank(1, 2, 0)
        if stat.startswith("two_2Cl_type_"):
            key = parts[3]
            if key == "larger":
                return 1 - sum(cl.prob_group(lam, 2, 0)
                               for n in range(MAX_TYPE_ORDER + 1)
                               for lam in cl.partitions(n))
            lam = () if key == "0" else tuple(int(x) for x in key.split(","))
            return cl.prob_group(lam, 2, 0)
    return None


def conditional_prediction(stat, u):
    """Gerth's finite-t prediction for two_rank2_t_rank4_r: P(4-rank = r |
    2-rank = t) = P(random t x t F_2-matrix has corank r) (imaginary)."""
    parts = stat.split("_")
    if u == 0 and stat.startswith("two_rank2_") and len(parts) == 5:
        return corank_dist_square(int(parts[2]), int(parts[4]))
    return None


# --------------------------------------------------------------------------
# p = 3: Davenport-Heilbronn main term and BST / Taniguchi-Thorne second term

def dh_constants(u, plimit=10 ** 7):
    """(M, K) with  sum_{fund. D, 0 < +-D < X} |Cl(D)[3]| = M X + K X^{5/6}
    + o(X^{5/6}):  M = 6/pi^2 (imag) or 4/pi^2 (real), and
    K = c * 8 zeta(1/3) / (5 Gamma(2/3)^3) * prod_p (1 - (p^{1/3}+1)/(p(p+1)))
    with c = sqrt(3) (imag) or 1 (real)."""
    from cypari import pari
    zeta13 = float(pari("zeta(1/3)"))
    gamma23 = math.gamma(2.0 / 3.0)
    sieve = np.ones(plimit + 1, dtype=bool)
    sieve[:2] = False
    for q in range(2, math.isqrt(plimit) + 1):
        if sieve[q]:
            sieve[q * q :: q] = False
    ps = np.nonzero(sieve)[0].astype(np.float64)
    logprod = np.sum(np.log1p(-(np.cbrt(ps) + 1) / (ps * (ps + 1))))
    logprod -= 1.5 * plimit ** (-2.0 / 3.0) / math.log(plimit)   # tail
    euler = math.exp(logprod)
    c = math.sqrt(3.0) if u == 0 else 1.0
    M = (6.0 if u == 0 else 4.0) / math.pi ** 2
    K = c * 8 * zeta13 / (5 * gamma23 ** 3) * euler
    return M, K


# --------------------------------------------------------------------------
# output

def gerth_mixture_rows(acc, last):
    """Imaginary fields: the 4-rank distribution predicted by Gerth's finite-t
    law averaged over the *observed* 2-ranks t = omega(D) - 1, i.e.
    sum_t P_emp(2-rank = t) * P(corank of random t x t F_2-matrix = r).
    Compared with the t -> oo limit this isolates the effect of omega(D)
    being small (~ log log X) at accessible X."""
    edges, N = acc.edges, acc.counts["N"][:last]
    ts = sorted(int(s.split("_")[2]) for s in acc.counts
                if s.startswith("two_rank2_") and len(s.split("_")) == 3)
    rows = []
    for r in range(4):
        stat = f"two_rank4_{r}"
        mix_win = np.zeros(last)
        for t in ts:
            mix_win += acc.counts[f"two_rank2_{t}"][:last] * corank_dist_square(t, r)
        mix_cum = np.cumsum(mix_win)
        cumN = np.cumsum(N)
        for row_i in range(last):
            if cumN[row_i] == 0:
                continue
            rows.append(dict(stat=stat + "_vs_mixture", X=int(edges[row_i + 1]),
                             N=int(cumN[row_i]),
                             cum_freq=np.cumsum(acc.counts[stat][:last])[row_i] / cumN[row_i],
                             cl=cl.prob_rank(r, 2, 0),
                             win_lo=int(edges[row_i]), win_N=int(N[row_i]),
                             win_freq=(acc.counts[stat][row_i] / N[row_i]) if N[row_i] else "",
                             gerth_mixture=mix_cum[row_i] / cumN[row_i],
                             gerth_mixture_win=(mix_win[row_i] / N[row_i]) if N[row_i] else ""))
    return rows


def write_outputs(family, acc, outdir):
    import os
    os.makedirs(outdir, exist_ok=True)
    u = 0 if family == "imag" else 1
    edges = acc.edges
    N = acc.counts["N"]
    used = np.nonzero(N)[0]
    last = used[-1] + 1
    stats = sorted(acc.counts)
    # upper end of the data: the last window may be only partly covered.  Use
    # the end of the computed range when the file names give it.
    end = acc.range_end if acc.range_end > acc.dmax else acc.dmax + 1
    xtop = min(int(edges[last]), end) if acc.dmax else int(edges[last])
    edges = edges.copy()
    edges[last] = xtop

    save_counts(acc, os.path.join(outdir, f"{family}_windows.csv"))

    M, K = dh_constants(u)
    cumN = np.cumsum(N[:last])
    rows = []
    for s in stats:
        if s == "N":
            continue
        cum = np.cumsum(acc.counts[s][:last])
        if "@" in s and not s.startswith("N@"):
            sub = acc.counts["N@" + s.split("@")[1]][:last]
            denom, wdenom = np.cumsum(sub), sub
        elif s.startswith("prime_") and s != "prime_N":
            denom = np.cumsum(acc.counts["prime_N"][:last])
            wdenom = acc.counts["prime_N"][:last]
        elif s.startswith("two_rank2_") and len(s.split("_")) == 5:
            t = s.split("_")[2]
            denom = np.cumsum(acc.counts[f"two_rank2_{t}"][:last])
            wdenom = acc.counts[f"two_rank2_{t}"][:last]
        elif s == "prime_N":
            denom, wdenom = cumN, N[:last]
        else:
            denom, wdenom = cumN, N[:last]
        pred = prediction(s, u)
        if pred is None:
            pred = conditional_prediction(s, u)
        win = acc.counts[s][:last]
        for i in range(last):
            if denom[i] == 0:
                continue
            row = dict(stat=s, X=int(edges[i + 1]), N=int(denom[i]),
                       cum_freq=cum[i] / denom[i],
                       win_lo=int(edges[i]),
                       win_freq=(win[i] / wdenom[i]) if wdenom[i] else "",
                       win_N=int(wdenom[i]),
                       cl=pred if pred is not None else "")
            if "both_div" in s:
                p, q = s.split("_")[0][1:].split("p")
                fp = np.cumsum(acc.counts[f"p{p}_div"][:last])[i] / denom[i]
                fq = np.cumsum(acc.counts[f"p{q}_div"][:last])[i] / denom[i]
                row["marginals_product"] = fp * fq
            if s == "p3_sum_p^rank":
                X = float(edges[i + 1])
                lo = float(edges[i])
                row["second_order"] = (M * X + K * X ** (5 / 6)) / cumN[i]
                row["second_order_win"] = ((M * (X - lo) + K * (X ** (5 / 6) - lo ** (5 / 6)))
                                           / N[i] if N[i] else "")
            rows.append(row)
    if u == 0:
        rows += gerth_mixture_rows(acc, last)
    fields = ["stat", "X", "N", "cum_freq", "cl", "win_lo", "win_N", "win_freq",
              "second_order", "second_order_win", "gerth_mixture", "gerth_mixture_win",
              "marginals_product"]
    with open(os.path.join(outdir, f"{family}_cumulative.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for row in rows:
            w.writerow(row)

    # summary at the largest X
    Xmax = int(edges[last])
    lines = [f"# {family}: {int(cumN[-1])} fields, |D| < {Xmax}\n",
             f"Unit rank u = {u}.  Davenport-Heilbronn/BST constants: "
             f"M = {M:.6f}, K = {K:.6f}\n",
             "| statistic | empirical | Cohen-Lenstra | emp/CL - 1 |",
             "|---|---|---|---|"]
    for row in rows:
        if row["X"] != Xmax or row["cl"] == "":
            continue
        emp, pr = row["cum_freq"], row["cl"]
        rel = f"{emp / pr - 1:+.4f}" if pr else ""
        extra = ""
        if row.get("second_order"):
            extra = f" (with X^5/6 term: {row['second_order']:.6f})"
        if row.get("marginals_product"):
            extra = f" (product of observed marginals: {row['marginals_product']:.6f})"
        if row.get("gerth_mixture"):
            extra = f" (finite-t Gerth mixture: {row['gerth_mixture']:.6f})"
        lines.append(f"| {row['stat']} | {emp:.6f} | {pr:.6f}{extra} | {rel} |")
    with open(os.path.join(outdir, f"{family}_summary.md"), "w") as f:
        f.write("\n".join(lines) + "\n")
    write_intervals(family, acc, edges, last, outdir, M, K)
    return rows


def interval_bounds(edges, last):
    """Decades [10^k, 10^(k+1)) from 100 up to the end of the data, as lists
    of window indices."""
    xtop = int(edges[last])
    bounds = [10 ** k for k in range(2, 20) if 10 ** k < xtop] + [xtop]
    out = []
    for lo, hi in zip(bounds[:-1], bounds[1:]):
        idx = [i for i in range(last) if edges[i] >= lo and edges[i + 1] <= hi]
        out.append((lo, hi, idx))
    return out


def write_intervals(family, acc, edges, last, outdir, M, K):
    """Distributions in each decade of |D|: <family>_intervals.md / .csv."""
    import os
    u = 0 if family == "imag" else 1
    ivs = interval_bounds(edges, last)
    tot = lambda st, idx: float(sum(acc.counts[st][i] for i in idx)) if st in acc.counts else 0.0
    groups = []        # (title, [(column label, stat, denominator stat)], prediction fn)

    groups.append(("P(p | h)", [(f"p = {p}", f"p{p}_div", "N") for p in ODD_PRIMES[:5]]))
    groups.append(("Mean of |Cl[p]|", [(f"p = {p}", f"p{p}_sum_p^rank", "N") for p in (3, 5, 7)]))
    groups.append(("3-rank", [(f"rank {r}", f"p3_rank_{r}", "N") for r in ("0", "1", "2", "3+")]))
    for p in (3, 5):
        keys = ["0", "1", "2", "1,1", "3", "2,1", "4"] if p == 3 else ["0", "1", "2", "1,1", "3"]
        cols = []
        for k in keys:
            lam = [] if k == "0" else [int(e) for e in k.split(",")]
            cols.append(("trivial" if not lam else " x ".join(f"Z/{p ** e}" for e in lam),
                         f"p{p}_type_{k}", "N"))
        groups.append((f"{p}-Sylow subgroup", cols))
    if u == 0:
        groups.append(("4-rank (Gerth)", [(f"4-rank {r}", f"two_rank4_{r}", "N") for r in range(4)]))
    else:
        groups.append(("Class number one", [("odd part of h trivial", "odd_part_trivial", "N"),
                                            ("h = 1, prime D", "prime_h_eq_1", "prime_N")]))
    lines = [f"# {family}: class group statistics by interval of |D|\n",
             "Frequencies among the fields with |D| in each interval; the last row is the "
             "Cohen-Lenstra prediction (Gerth's for 4-ranks).\n"]
    csv_rows = []
    for title, cols in groups:
        head = "| interval | fields | " + " | ".join(c[0] for c in cols) + " |"
        lines += [f"\n## {title}\n", head, "|" + "---|" * (len(cols) + 2)]
        for lo, hi, idx in ivs:
            n = tot("N", idx)
            if not n:
                continue
            cells = []
            for label, st, den in cols:
                d = tot(den, idx)
                v = tot(st, idx) / d if d else float("nan")
                cells.append(f"{v:.5f}")
                csv_rows.append([lo, hi, st, tot(st, idx), d, v])
            extra = ""
            if title.startswith("Mean") and idx:
                a, b = float(edges[idx[0]]), float(edges[idx[-1] + 1])
                extra = f" (2-term: {(M * (b - a) + K * (b ** (5 / 6) - a ** (5 / 6))) / n:.5f})"
                cells[0] += extra
            # 3 significant digits, more only where they are needed (1.004e10)
            g = lambda x: f"{x:.3g}" if float(f"{x:.3g}") == x else f"{x:.6g}"
            lines.append(f"| [{g(lo)}, {g(hi)}) | {int(n):,} | " + " | ".join(cells) + " |")
        preds = []
        for label, st, den in cols:
            pr = prediction(st, u)
            preds.append(f"{pr:.5f}" if pr is not None else "")
        lines.append("| **CL** | | " + " | ".join(preds) + " |")
    with open(os.path.join(outdir, f"{family}_intervals.md"), "w") as f:
        f.write("\n".join(lines).replace("e+0", "e").replace("e+", "e") + "\n")
    with open(os.path.join(outdir, f"{family}_intervals.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["lo", "hi", "stat", "count", "denominator", "freq"])
        w.writerows(csv_rows)


def save_counts(acc, path):
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["lo", "hi", "stat", "count"])
        for i in range(len(acc.edges) - 1):
            for st in sorted(acc.counts):
                c = acc.counts[st][i]
                # keep statistics that never occur (one zero row) so that
                # merged runs know about them
                if c or (i == 0 and not acc.counts[st].any()):
                    w.writerow([acc.edges[i], acc.edges[i + 1], st,
                                int(c) if c == int(c) else repr(float(c))])
        for name, v in (("max_absD", acc.dmax), ("range_end", acc.range_end)):
            if v:             # merged with max(), not summed
                i = int(np.searchsorted(acc.edges, v - 1, side="right")) - 1
                w.writerow([acc.edges[i], acc.edges[i + 1], name, v])


def load_counts(paths):
    acc = Accumulator(window_edges())
    index = {int(lo): i for i, lo in enumerate(acc.edges[:-1])}
    for path in paths:
        with open(path) as f:
            for row in csv.DictReader(f):
                if row["stat"] == "max_absD":
                    acc.dmax = max(acc.dmax, int(float(row["count"])))
                    continue
                if row["stat"] == "range_end":
                    acc.range_end = max(acc.range_end, int(float(row["count"])))
                    continue
                acc.counts[row["stat"]][index[int(row["lo"])]] += float(row["count"])
    return acc


def main():
    if len(sys.argv) < 5 or sys.argv[1] not in ("all", "count", "report") \
            or sys.argv[2] not in ("imag", "real"):
        print(__doc__)
        sys.exit(1)
    mode, family, target, paths = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4:]
    if mode == "report":
        write_outputs(family, load_counts(paths), target)
        return
    acc = process(family, paths, window_edges())
    if mode == "count":
        save_counts(acc, target)
    else:
        write_outputs(family, acc, target)


if __name__ == "__main__":
    main()
