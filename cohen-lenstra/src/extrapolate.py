#!/usr/bin/env python3
"""
Extrapolate empirical class group statistics to X = infinity and compare the
limit with Cohen-Lenstra.

For each statistic F(X) (cumulative frequency over |D| < X, from
<family>_cumulative.csv) we fit, over X >= XMIN,

  (a) F(X) = F_oo + a X^(-1/6)               (the proven shape for p = 3)
  (b) F(X) = F_oo + a X^(-theta)              (theta free)

and report F_oo next to the Cohen-Lenstra value.  These are least-squares
fits to log-spaced points (cumulative values are strongly correlated, so treat
the quoted fit errors as indicative only).

usage: extrapolate.py RESULTS_DIR [XMIN]    (writes RESULTS_DIR/extrapolation.md)
"""
import csv
import math
import os
import sys
from collections import defaultdict

import numpy as np
from scipy.optimize import curve_fit

THETA_MIN = 0.05
STATS = ["p3_div", "p5_div", "p7_div", "p11_div", "p3_rank_2", "p5_rank_2",
         "p3_sum_p^rank", "p5_sum_p^rank", "p7_sum_p^rank",
         "p3_div_p2", "odd_part_cyclic", "prime_h_eq_1"]


def load(path):
    d = defaultdict(list)
    with open(path) as f:
        for row in csv.DictReader(f):
            d[row["stat"]].append(row)
    return d


def fit(xs, ys):
    xs, ys = np.asarray(xs, float), np.asarray(ys, float)
    A = np.vstack([np.ones_like(xs), xs ** (-1 / 6)]).T
    (c0, c1), *_ = np.linalg.lstsq(A, ys, rcond=None)
    out = {"lim16": c0, "a16": c1}
    try:
        popt, pcov = curve_fit(lambda x, c, a, t: c + a * x ** (-t), xs, ys,
                               p0=(c0, c1, 1 / 6), maxfev=20000,
                               bounds=([-np.inf, -np.inf, THETA_MIN], [np.inf, np.inf, 1.0]))
        if THETA_MIN + 1e-3 < popt[2] < 1.0 - 1e-3:   # drop fits stuck at a bound
            out.update(limfree=popt[0], afree=popt[1], theta=popt[2],
                       err=math.sqrt(max(pcov[0, 0], 0)))
    except Exception:  # noqa: BLE001
        pass
    return out


def main():
    results = sys.argv[1] if len(sys.argv) > 1 else "results"
    xmin = float(sys.argv[2]) if len(sys.argv) > 2 else 1e5
    lines = ["# Extrapolation to X = infinity\n",
             f"Fits to the cumulative frequency F(X) over |D| < X, using X >= {xmin:.0e}.\n",
             "(a): F = F_oo + a X^(-1/6).  (b): F = F_oo + a X^(-theta), theta free in [0.05, 1] "
             "(left blank when the fit runs into a bound).\n"]
    for fam in ("imag", "real"):
        path = os.path.join(results, f"{fam}_cumulative.csv")
        if not os.path.exists(path):
            continue
        d = load(path)
        lines += [f"\n## {fam}\n",
                  "| statistic | F at largest X | F_oo (a) | F_oo (b) | theta (b) | Cohen-Lenstra |",
                  "|---|---|---|---|---|---|"]
        for s in STATS:
            if s not in d or not d[s][0]["cl"]:
                continue
            pts = [(float(r["X"]), float(r["cum_freq"])) for r in d[s]
                   if float(r["X"]) >= xmin and r["cum_freq"]]
            if len(pts) < 5:
                continue
            xs, ys = zip(*pts)
            f = fit(xs, ys)
            cl = float(d[s][0]["cl"])
            free = (f"{f['limfree']:.5f}", f"{f['theta']:.3f}") if "theta" in f else ("", "")
            lines.append(f"| {s} | {ys[-1]:.5f} | {f['lim16']:.5f} | {free[0]} | {free[1]} | {cl:.5f} |")
    out = os.path.join(results, "extrapolation.md")
    with open(out, "w") as fh:
        fh.write("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
