#!/usr/bin/env python3
"""
Convergence plots from the CSV files written by analyze.py.

usage: plots.py RESULTS_DIR       (reads imag_/real_cumulative.csv, writes RESULTS_DIR/plots/*.png)
"""
import csv
import math
import os
import sys
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.ticker import FuncFormatter, LogLocator, NullFormatter  # noqa: E402

sys.path.insert(0, os.path.dirname(__file__))
import cl_predictions as cl  # noqa: E402

# reference palette (light mode) -- categorical slots in fixed order
SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100"]
SURFACE, INK, INK2, MUTED = "#fcfcfb", "#0b0b0b", "#52514e", "#898781"
GRID, AXIS = "#e1e0d9", "#c3c2b7"
LW = 1.6          # ~2px lines
MS = 5.5          # >= 8px markers

plt.rcParams.update({
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE,
    "savefig.facecolor": SURFACE, "font.family": "DejaVu Sans",
    "font.size": 10, "text.color": INK, "axes.labelcolor": INK2,
    "axes.edgecolor": AXIS, "axes.linewidth": 0.8,
    "xtick.color": MUTED, "ytick.color": MUTED,
    "xtick.labelcolor": INK2, "ytick.labelcolor": INK2,
    "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.8,
    "grid.linestyle": "-", "axes.spines.top": False, "axes.spines.right": False,
    "axes.titlesize": 11, "axes.titleweight": "bold", "axes.titlecolor": INK,
    "axes.titlelocation": "left", "legend.frameon": False,
    "legend.labelcolor": INK2, "lines.solid_capstyle": "round",
    "lines.solid_joinstyle": "round",
})

FAMILY_TITLE = {"imag": "Imaginary quadratic fields (D < 0)",
                "real": "Real quadratic fields (D > 0)"}


def load(results):
    data = {}
    for fam in ("imag", "real"):
        path = os.path.join(results, f"{fam}_cumulative.csv")
        if not os.path.exists(path):
            continue
        d = defaultdict(list)
        with open(path) as f:
            for row in csv.DictReader(f):
                d[row["stat"]].append(row)
        data[fam] = d
    return data


def series(rows, key="cum_freq", xkey="X", min_n=1000, nkey="N"):
    xs, ys = [], []
    for r in rows:
        if r[key] in ("", None) or int(r[nkey]) < min_n:
            continue
        xs.append(float(r[xkey]))
        ys.append(float(r[key]))
    return xs, ys


def pow10_formatter(x, _pos):
    if x <= 0:
        return ""
    e = math.log10(x)
    if abs(e - round(e)) < 1e-9:
        return f"$10^{{{int(round(e))}}}$"
    return ""


def style_logx(ax):
    ax.set_xscale("log")
    ax.xaxis.set_major_locator(LogLocator(base=10))
    ax.xaxis.set_major_formatter(FuncFormatter(pow10_formatter))
    ax.xaxis.set_minor_formatter(NullFormatter())
    ax.grid(True, which="major")
    ax.grid(False, which="minor")


def end_label(ax, x, y, text, dy=0):
    ax.annotate(text, (x, y), xytext=(6, dy), textcoords="offset points",
                va="center", ha="left", fontsize=9, color=INK2)


def fit_slope(xs, ys, xmin):
    """Least-squares slope of log y against log x for x >= xmin."""
    pts = [(math.log(x), math.log(y)) for x, y in zip(xs, ys) if x >= xmin and y > 0]
    n = len(pts)
    mx = sum(a for a, _ in pts) / n
    my = sum(b for _, b in pts) / n
    return (sum((a - mx) * (b - my) for a, b in pts)
            / sum((a - mx) ** 2 for a, _ in pts))


def fmt_x(X):
    """10^k in mathtext when X is a power of 10, else with thousands commas."""
    k = math.log10(X)
    if abs(k - round(k)) < 1e-9:
        return f"$10^{{{int(round(k))}}}$"
    return f"{X:,}"


def save(fig, outdir, name):
    path = os.path.join(outdir, name)
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print("wrote", path)


# --------------------------------------------------------------------------

def plot_divisibility(data, outdir):
    """Relative deviation of P(p | h) from Cohen-Lenstra, p = 3, 5, 7."""
    fams = [f for f in ("imag", "real") if f in data]
    fig, axes = plt.subplots(1, len(fams), figsize=(5.6 * len(fams), 4.0), squeeze=False)
    for ax, fam in zip(axes[0], fams):
        for i, p in enumerate((3, 5, 7)):
            rows = data[fam][f"p{p}_div"]
            pred = float(rows[0]["cl"])
            xs, ys = series(rows)
            dev = [100 * (y / pred - 1) for y in ys]
            ax.plot(xs, dev, color=SERIES[i], lw=LW,
                    label=f"p = {p}: CL {pred:.4f}, data {ys[-1]:.4f} ({dev[-1]:+.2f}%)")
        ax.axhline(0, color=INK2, lw=1.0)
        ax.annotate("Cohen–Lenstra", (xs[0], 0), xytext=(0, 5), textcoords="offset points",
                    fontsize=9, color=INK2)
        style_logx(ax)
        ax.set_xlabel("X  (fields with |D| < X)")
        ax.set_ylabel("P(p | h): empirical / CL − 1  (%)")
        ax.set_title(FAMILY_TITLE[fam])
        ax.legend(loc="lower right", fontsize=9)
    fig.suptitle("Proportion of fields with p | h, relative to the Cohen–Lenstra prediction",
                 x=0.01, ha="left", fontsize=12, fontweight="bold", color=INK)
    fig.tight_layout()
    save(fig, outdir, "divisibility_convergence.png")


def plot_three_torsion(data, outdir):
    """Mean of |Cl[3]| vs Davenport-Heilbronn limit and the BST/TT second-order term."""
    fams = [f for f in ("imag", "real") if f in data]
    fig, axes = plt.subplots(1, len(fams), figsize=(5.6 * len(fams), 4.0), squeeze=False)
    for ax, fam in zip(axes[0], fams):
        rows = data[fam]["p3_sum_p^rank"]
        limit = float(rows[0]["cl"])
        xs, ys = series(rows, min_n=200)
        xs2, ys2 = series(rows, key="second_order", min_n=200)
        ax.plot(xs2, ys2, color=SERIES[1], lw=LW,
                label="main term + $X^{5/6}$ term (BST, Taniguchi–Thorne)")
        ax.plot(xs, ys, "o", color=SERIES[0], ms=MS, mec=SURFACE, mew=1.2,
                label="data (all fields with |D| < X)")
        ax.axhline(limit, color=INK2, lw=1.0)
        ax.annotate(f"Davenport–Heilbronn / CL limit = {limit:g}" if fam == "imag"
                    else "Davenport–Heilbronn / CL limit = 4/3",
                    (xs[0], limit), xytext=(0, 5), textcoords="offset points",
                    fontsize=9, color=INK2)
        end_label(ax, xs[-1], ys[-1], f"{ys[-1]:.4f}")
        style_logx(ax)
        ax.set_xlabel("X")
        ax.set_ylabel("mean of |Cl(D)[3]| over |D| < X")
        ax.set_title(FAMILY_TITLE[fam])
        ax.legend(loc="lower right", fontsize=8.5)
    fig.suptitle("3-torsion: the gap to Cohen–Lenstra is the proven $X^{-1/6}$ secondary term",
                 x=0.01, ha="left", fontsize=12, fontweight="bold", color=INK)
    fig.tight_layout()
    save(fig, outdir, "three_torsion_mean.png")


def plot_deviation_scaling(data, outdir):
    """|emp/CL - 1| for E|Cl[p]| (windowed) on log-log axes with an X^{-1/6} guide."""
    fams = [f for f in ("imag", "real") if f in data]
    fig, axes = plt.subplots(1, len(fams), figsize=(5.6 * len(fams), 4.0), squeeze=False)
    for ax, fam in zip(axes[0], fams):
        last = None
        for i, p in enumerate((3, 5, 7)):
            rows = data[fam][f"p{p}_sum_p^rank"]
            pred = float(rows[0]["cl"])
            pts = [(float(r["X"]), abs(float(r["cum_freq"]) / pred - 1)) for r in rows
                   if int(r["N"]) >= 5000 and float(r["cum_freq"]) != pred]
            xs, ys = zip(*pts)
            slope = fit_slope(xs, ys, xs[-1] / 100)
            ax.plot(xs, ys, color=SERIES[i], lw=LW,
                    label=f"p = {p}  (log-log slope over last 2 decades: {slope:.2f})")
            last = (xs, ys) if p == 3 else last
        xs, ys = last
        x0, y0 = xs[len(xs) // 2], ys[len(ys) // 2] * 0.62
        guide = [y0 * (x / x0) ** (-1 / 6) for x in xs]
        ax.plot(xs, guide, color=MUTED, lw=1.0)
        ax.annotate("reference slope −1/6", (xs[-1], guide[-1]), xytext=(-4, -12),
                    textcoords="offset points", fontsize=9, color=INK2, ha="right")
        style_logx(ax)
        ax.set_yscale("log")
        ax.set_xlabel("X")
        ax.set_ylabel("|mean of |Cl[p]| / CL − 1|   (|D| < X)")
        ax.set_title(FAMILY_TITLE[fam])
        ax.legend(loc="lower left", fontsize=8.5)
    fig.suptitle("How fast the p-torsion averages approach Cohen–Lenstra",
                 x=0.01, ha="left", fontsize=12, fontweight="bold", color=INK)
    fig.tight_layout()
    save(fig, outdir, "deviation_scaling.png")


def plot_sylow_types(data, outdir, p=3):
    """Dot plot: frequency of each p-Sylow type at the largest X vs CL."""
    fams = [f for f in ("imag", "real") if f in data]
    fig, axes = plt.subplots(1, len(fams), figsize=(5.6 * len(fams), 4.2), squeeze=False)
    for ax, fam in zip(axes[0], fams):
        types = []
        for stat, rows in data[fam].items():
            if stat.startswith(f"p{p}_type_") and not stat.endswith("larger"):
                r = rows[-1]
                key = stat[len(f"p{p}_type_"):]
                emp, pred = float(r["cum_freq"]), float(r["cl"])
                types.append((pred, key, emp))
        types.sort(reverse=True)
        names = []
        for j, (pred, key, emp) in enumerate(types):
            lam = [] if key == "0" else [int(e) for e in key.split(",")]
            names.append("trivial" if not lam else
                         " × ".join(f"Z/{p ** e}" for e in lam))
            ax.plot([pred, emp], [j, j], color=GRID, lw=2, zorder=1)
            ax.plot(pred, j, "o", color=SERIES[1], ms=MS, mec=SURFACE, mew=1.2, zorder=2,
                    label="Cohen–Lenstra" if j == 0 else None)
            if emp > 0:
                ax.plot(emp, j, "o", color=SERIES[0], ms=MS, mec=SURFACE, mew=1.2, zorder=3,
                        label="data" if j == 0 else None)
        ax.set_yticks(range(len(names)))
        ax.set_yticklabels(names)
        ax.invert_yaxis()
        ax.set_xscale("log")
        ax.grid(True, axis="x", which="major")
        ax.grid(False, axis="y")
        N = int(types and data[fam][f"p{p}_type_0"][-1]["N"])
        X = int(data[fam][f"p{p}_type_0"][-1]["X"])
        ax.set_xlabel(f"proportion of fields  ({N:,} fields, |D| < {fmt_x(X)})")
        ax.set_title(FAMILY_TITLE[fam])
        ax.legend(loc="lower right", fontsize=9)
    fig.suptitle(f"{p}-Sylow subgroup of the class group: data vs Cohen–Lenstra",
                 x=0.01, ha="left", fontsize=12, fontweight="bold", color=INK)
    fig.tight_layout()
    save(fig, outdir, f"sylow{p}_types.png")


def plot_four_rank(data, outdir):
    """Imaginary fields: 4-rank distribution vs Gerth limit and the finite-t mixture."""
    if "imag" not in data:
        return
    d = data["imag"]
    rs = [r for r in range(4) if f"two_rank4_{r}_vs_mixture" in d]
    emp = [float(d[f"two_rank4_{r}_vs_mixture"][-1]["cum_freq"]) for r in rs]
    lim = [float(d[f"two_rank4_{r}_vs_mixture"][-1]["cl"]) for r in rs]
    mix = [float(d[f"two_rank4_{r}_vs_mixture"][-1]["gerth_mixture"]) for r in rs]
    X = int(d["two_rank4_0_vs_mixture"][-1]["X"])
    fig, ax = plt.subplots(figsize=(6.4, 4.0))
    w = 0.26
    for k, (vals, lab) in enumerate(((emp, "data"),
                                     (mix, "Gerth, finite t = 2-rank (mixture over data)"),
                                     (lim, "Gerth / Cohen–Lenstra limit (t → ∞)"))):
        xs = [r + (k - 1) * w for r in rs]
        ax.bar(xs, vals, width=w * 0.9, color=SERIES[k], label=lab, zorder=2)
        for x, v in zip(xs, vals):
            ax.annotate(f"{v:.3f}" if v >= 0.0005 else "<0.001", (x, v), xytext=(0, 3),
                        textcoords="offset points", ha="center", fontsize=7.5, color=INK2)
    ax.set_xticks(rs)
    ax.set_xticklabels([f"4-rank {r}" for r in rs])
    ax.grid(True, axis="y")
    ax.grid(False, axis="x")
    ax.set_ylabel("proportion of fields")
    ax.set_title(f"2-part, imaginary fields with |D| < {fmt_x(X)}: 4-rank of Cl(D)")
    ax.set_ylim(0, max(emp + mix + lim) * 1.38)
    ax.legend(loc="upper right", fontsize=8.5)
    fig.tight_layout()
    save(fig, outdir, "four_rank.png")


def plot_real_h1(data, outdir):
    """Real fields: odd part of h trivial (all fields) and h = 1 (prime discriminants)."""
    if "real" not in data:
        return
    d = data["real"]
    pred = cl.prob_odd_part_trivial_real()
    fig, ax = plt.subplots(figsize=(6.4, 4.0))
    for i, (stat, lab) in enumerate((("prime_h_eq_1", "h = 1, prime discriminants D = p"),
                                     ("odd_part_trivial", "odd part of h trivial, all D > 0"))):
        xs, ys = series(d[stat], min_n=200)
        ax.plot(xs, [100 * y for y in ys], color=SERIES[i], lw=LW, label=lab)
        end_label(ax, xs[-1], 100 * ys[-1], f"{100 * ys[-1]:.2f}%")
    ax.axhline(100 * pred, color=INK2, lw=1.0)
    ax.annotate(f"Cohen–Lenstra {100 * pred:.3f}%", (xs[0], 100 * pred), xytext=(0, 5),
                textcoords="offset points", fontsize=9, color=INK2)
    style_logx(ax)
    ax.set_xlabel("X  (fields with D < X)")
    ax.set_ylabel("% of fields")
    ax.set_title("Real quadratic fields: class number one")
    ax.legend(loc="lower right", fontsize=9)
    fig.tight_layout()
    save(fig, outdir, "real_class_number_one.png")


def main():
    results = sys.argv[1] if len(sys.argv) > 1 else "results"
    outdir = os.path.join(results, "plots")
    os.makedirs(outdir, exist_ok=True)
    data = load(results)
    plot_divisibility(data, outdir)
    plot_three_torsion(data, outdir)
    plot_deviation_scaling(data, outdir)
    plot_sylow_types(data, outdir, 3)
    plot_sylow_types(data, outdir, 5)
    plot_four_rank(data, outdir)
    plot_real_h1(data, outdir)


if __name__ == "__main__":
    main()
