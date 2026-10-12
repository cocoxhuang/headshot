#!/usr/bin/env python3
"""
Build a self-contained HTML report (interactive charts) from the CSV files
written by analyze.py.

usage: make_report.py RESULTS_DIR OUT.html
"""
import csv
import json
import os
import sys
from collections import defaultdict

sys.path.insert(0, os.path.dirname(__file__))
import cl_predictions as cl  # noqa: E402


def load(results, fam):
    path = os.path.join(results, f"{fam}_cumulative.csv")
    if not os.path.exists(path):
        return None
    d = defaultdict(list)
    with open(path) as f:
        for row in csv.DictReader(f):
            d[row["stat"]].append(row)
    return d


def fnum(x):
    return None if x in ("", None) else float(x)


INTERVAL_GROUPS = {
    "div": ("P(p | h)", [(f"p = {p}", f"p{p}_div", "N") for p in (3, 5, 7, 11, 13)]),
    "rank3": ("3-rank", [(f"rank {r}", f"p3_rank_{r}", "N") for r in ("0", "1", "2", "3+")]),
    "syl3": ("3-Sylow subgroup", [(n, f"p3_type_{k}", "N") for n, k in
                                 (("trivial", "0"), ("Z/3", "1"), ("Z/9", "2"), ("Z/27", "3"),
                                  ("Z/3 × Z/3", "1,1"), ("Z/9 × Z/3", "2,1"))]),
    "rank4": ("4-rank", [(f"4-rank {r}", f"two_rank4_{r}", "N") for r in range(4)]),
    "h1": ("class number one", [("odd part trivial", "odd_part_trivial", "N"),
                                ("h = 1, prime D", "prime_h_eq_1", "prime_N")]),
}


def interval_payload(results, fam):
    from analyze import prediction
    path = os.path.join(results, f"{fam}_intervals.csv")
    if not os.path.exists(path):
        return None
    u = 0 if fam == "imag" else 1
    vals, ns = {}, {}
    with open(path) as f:
        for r in csv.DictReader(f):
            key = (int(r["lo"]), int(r["hi"]))
            vals[key + (r["stat"],)] = float(r["freq"])
            if r["stat"] == "p3_div":
                ns[key] = int(float(r["denominator"]))
    ivs = sorted(ns)
    out = {"intervals": [[lo, hi, ns[(lo, hi)]] for lo, hi in ivs], "groups": {}}
    for gk, (title, cols) in INTERVAL_GROUPS.items():
        if gk == "rank4" and u == 1 or gk == "h1" and u == 0:
            continue
        out["groups"][gk] = {
            "title": title, "cols": [c[0] for c in cols],
            "rows": [[vals.get((lo, hi, c[1])) for c in cols] for lo, hi in ivs],
            "cl": [prediction(c[1], u) for c in cols]}
    return out


def family_payload(d, fam):
    u = 0 if fam == "imag" else 1
    rows_N = d["p3_div"]
    xs = [int(r["X"]) for r in rows_N if int(r["N"]) >= 200]
    keep = set(xs)

    def ser(stat, key="cum_freq"):
        return [fnum(r[key]) for r in d[stat] if int(r["X"]) in keep]

    out = {"X": xs, "N": [int(r["N"]) for r in rows_N if int(r["X"]) in keep],
           "series": {}}
    for st in ("p3_div", "p5_div", "p7_div", "p11_div", "p3_sum_p^rank", "p5_sum_p^rank",
               "p7_sum_p^rank", "p3_rank_2", "odd_part_cyclic", "prime_h_eq_1",
               "odd_part_trivial", "p3_div_p2", "p3_sum_p^2rank"):
        if st not in d:
            continue
        clv = fnum(d[st][0]["cl"])
        entry = {"cum": ser(st), "win": ser(st, "win_freq"), "cl": clv}
        if st == "p3_sum_p^rank":
            entry["second"] = ser(st, "second_order")
            entry["second_win"] = ser(st, "second_order_win")
        if st.startswith("prime_"):
            entry["n"] = [int(r["N"]) for r in d[st] if int(r["X"]) in keep]
        out["series"][st] = entry
    # p-Sylow types at the largest X
    for p in (3, 5):
        types = []
        for st, rows in d.items():
            pref = f"p{p}_type_"
            if st.startswith(pref) and not st.endswith("larger"):
                key = st[len(pref):]
                lam = [] if key == "0" else [int(e) for e in key.split(",")]
                name = "trivial" if not lam else " × ".join(f"Z/{p ** e}" for e in lam)
                types.append({"name": name, "cl": fnum(rows[-1]["cl"]),
                              "emp": fnum(rows[-1]["cum_freq"]), "order": p ** sum(lam)})
        types.sort(key=lambda t: -t["cl"])
        out[f"sylow{p}"] = types
    # p-rank distributions at largest X
    out["ranks"] = {}
    for p in (3, 5, 7):
        out["ranks"][p] = [{"r": r, "emp": fnum(d[f"p{p}_rank_{r}"][-1]["cum_freq"]),
                            "cl": fnum(d[f"p{p}_rank_{r}"][-1]["cl"])}
                           for r in ("0", "1", "2", "3+")]
    # 4-ranks (imaginary: Gerth limit and finite-t mixture)
    if u == 0 and "two_rank4_0_vs_mixture" in d:
        out["rank4"] = [{"r": r, "emp": fnum(d[f"two_rank4_{r}_vs_mixture"][-1]["cum_freq"]),
                         "mix": fnum(d[f"two_rank4_{r}_vs_mixture"][-1]["gerth_mixture"]),
                         "lim": fnum(d[f"two_rank4_{r}_vs_mixture"][-1]["cl"])}
                        for r in range(4)]
        out["rank2"] = [{"t": t, "frac": fnum(d[f"two_rank2_{t}"][-1]["cum_freq"])}
                        for t in range(8) if f"two_rank2_{t}" in d]
    # independence
    out["pairs"] = [{"pair": st.split("_")[0], "emp": fnum(d[st][-1]["cum_freq"]),
                     "marg": fnum(d[st][-1]["marginals_product"]),
                     "cl": fnum(d[st][-1]["cl"])}
                    for st in ("p3p5_both_div", "p3p7_both_div", "p5p7_both_div") if st in d]
    # local conditions: behaviour of p in K
    out["local"] = []
    for p in (3, 5):
        for kind in ("split", "inert", "ram"):
            st = f"p{p}_sum_p^rank@{p}{kind}"
            if st in d:
                r = d[st][-1]
                out["local"].append({"p": p, "kind": kind, "N": int(r["N"]),
                                     "mean": fnum(r["cum_freq"]), "cl": fnum(r["cl"]),
                                     "div": fnum(d[f"p{p}_div@{p}{kind}"][-1]["cum_freq"]),
                                     "cl_div": fnum(d[f"p{p}_div@{p}{kind}"][-1]["cl"])})
    out["Nmax"] = out["N"][-1]
    out["Xmax"] = out["X"][-1]
    return out


TEMPLATE = r"""<title>Cohen–Lenstra Convergence</title>
<link rel="preconnect" href="https://fonts.googleapis.com">
<link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=STIX+Two+Text:ital,wght@0,400;0,600;1,400&family=IBM+Plex+Sans:wght@400;500;600&family=IBM+Plex+Mono:wght@400;500&display=swap">
<style>
/* Layout: a single reading column (~70ch) with full-width figures; figures are the content, prose frames them. */
:root {
  --bg: #f6f7f9; --surface: #fdfdfc; --ink: #15171c; --ink-2: #4c5160; --muted: #858a96;
  --rule: #dfe2e8; --grid: #e6e8ec; --axis: #c3c6ce; --accent: #2b4f9e;
  --s1: #2a78d6; --s2: #eb6834; --s3: #1baf7a; --s4: #eda100;
  --chip: #eceff4;
  --f-display: "STIX Two Text", "Times New Roman", Georgia, serif;
  --f-body: "IBM Plex Sans", system-ui, -apple-system, "Segoe UI", sans-serif;
  --f-mono: "IBM Plex Mono", ui-monospace, "SFMono-Regular", Menlo, monospace;
}
@media (prefers-color-scheme: dark) {
  :root:not([data-theme="light"]) {
    --bg: #121418; --surface: #1a1c21; --ink: #f1f2f5; --ink-2: #bfc3cd; --muted: #8a8f9b;
    --rule: #2c2f36; --grid: #2a2d33; --axis: #3a3e46; --accent: #8fb0f2;
    --s1: #3987e5; --s2: #d95926; --s3: #199e70; --s4: #c98500; --chip: #23262d;
    color-scheme: dark;
  }
}
:root[data-theme="dark"] {
  --bg: #121418; --surface: #1a1c21; --ink: #f1f2f5; --ink-2: #bfc3cd; --muted: #8a8f9b;
  --rule: #2c2f36; --grid: #2a2d33; --axis: #3a3e46; --accent: #8fb0f2;
  --s1: #3987e5; --s2: #d95926; --s3: #199e70; --s4: #c98500; --chip: #23262d;
  color-scheme: dark;
}
* { box-sizing: border-box; }
body { background: var(--bg); color: var(--ink); font: 15px/1.6 var(--f-body); }
.page { max-width: 980px; margin: 0 auto; padding-inline: 20px; padding-block: 40px 72px; display: grid; gap: 44px; }
.prose { max-width: 68ch; }
h1, h2, h3 { font-family: var(--f-display); font-weight: 600; text-wrap: balance; margin: 0; line-height: 1.2; }
h1 { font-size: clamp(30px, 5vw, 44px); letter-spacing: -0.01em; }
h2 { font-size: 26px; }
h3 { font-size: 18px; }
p { margin: 0; }
.stack { display: grid; gap: 14px; }
.eyebrow { font: 500 12px/1.4 var(--f-body); letter-spacing: 0.08em; text-transform: uppercase; color: var(--ink-2); }
.lede { font-size: 17px; color: var(--ink-2); max-width: 62ch; }
.chips { display: flex; flex-wrap: wrap; gap: 8px; }
.chip { background: var(--chip); border-radius: 6px; padding: 6px 10px; font-size: 13px; color: var(--ink-2); }
.chip b { color: var(--ink); font-weight: 600; font-variant-numeric: tabular-nums; }
.math { font-family: var(--f-display); font-style: italic; }
.formula { font-family: var(--f-display); font-size: 17px; background: var(--surface); border: 1px solid var(--rule); border-radius: 8px; padding: 14px 16px; overflow-x: auto; }
.figure { background: var(--surface); border: 1px solid var(--rule); border-radius: 10px; padding: 18px 18px 12px; display: grid; gap: 10px; min-width: 0; }
.figure-head { display: flex; flex-wrap: wrap; justify-content: space-between; align-items: baseline; gap: 10px; }
.caption { font-size: 13px; color: var(--ink-2); max-width: 80ch; }
.controls { display: flex; flex-wrap: wrap; gap: 6px; }
.seg { display: inline-flex; border: 1px solid var(--rule); border-radius: 7px; overflow: hidden; }
.seg button { font: 500 12.5px var(--f-body); color: var(--ink-2); background: transparent; border: 0; padding: 6px 11px; cursor: pointer; }
.seg button[aria-pressed="true"] { background: var(--chip); color: var(--ink); }
.seg button:focus-visible { outline: 2px solid var(--accent); outline-offset: -2px; }
.chart { position: relative; width: 100%; min-height: 260px; }
.chart svg { display: block; width: 100%; height: auto; overflow: visible; }
.chart text { font-family: var(--f-body); }
.legend { display: flex; flex-wrap: wrap; gap: 6px 16px; font-size: 12.5px; color: var(--ink-2); }
.legend span { display: inline-flex; align-items: center; gap: 6px; }
.key-line { width: 18px; height: 2px; border-radius: 1px; display: inline-block; }
.key-dot { width: 9px; height: 9px; border-radius: 50%; display: inline-block; }
.tip { position: absolute; pointer-events: none; background: var(--surface); border: 1px solid var(--rule); border-radius: 7px; padding: 8px 10px; font-size: 12px; color: var(--ink); box-shadow: 0 4px 14px rgba(0,0,0,0.12); min-width: 160px; z-index: 2; }
.tip .row { display: flex; justify-content: space-between; gap: 14px; font-variant-numeric: tabular-nums; }
.tip .row span:first-child { display: inline-flex; align-items: center; gap: 6px; color: var(--ink-2); }
.grid2 { display: grid; grid-template-columns: repeat(auto-fit, minmax(min(100%, 420px), 1fr)); gap: 18px; }
.tablewrap { overflow-x: auto; }
table { border-collapse: collapse; width: 100%; font-size: 13.5px; }
th, td { padding: 6px 10px; text-align: right; border-bottom: 1px solid var(--rule); white-space: nowrap; font-variant-numeric: tabular-nums; }
th:first-child, td:first-child { text-align: left; }
th { font-weight: 500; color: var(--ink-2); font-size: 12.5px; }
td.neg { color: var(--ink); }
.dev { font-family: var(--f-mono); font-size: 12.5px; }
ol.consider { margin: 0; padding-left: 1.4em; display: grid; gap: 14px; max-width: 74ch; }
ol.consider li::marker { font-family: var(--f-display); font-weight: 600; color: var(--accent); }
ol.consider b { font-weight: 600; }
code { font-family: var(--f-mono); font-size: 0.9em; background: var(--chip); padding: 1px 5px; border-radius: 4px; }
pre { font-family: var(--f-mono); font-size: 12.5px; background: var(--surface); border: 1px solid var(--rule); border-radius: 8px; padding: 12px 14px; overflow-x: auto; margin: 0; }
.foot { font-size: 12.5px; color: var(--muted); }
a { color: var(--accent); }
@media (max-width: 520px) { .page { padding-inline: 16px; gap: 34px; } .figure { padding: 14px 12px 10px; } }
</style>

<main class="page">
  <header class="stack">
    <div class="eyebrow">Class groups of quadratic fields · census to 10<sup>9</sup></div>
    <h1>Cohen–Lenstra Convergence</h1>
    <p class="lede" id="lede">Every imaginary and every real quadratic field up to the bounds below, compared with the Cohen–Lenstra predictions. The data does approach the predictions, but slowly. For p = 3 the remaining gap matches the proven secondary term.</p>
    <div class="chips" id="chips"></div>
  </header>

  <section class="stack">
    <h2>How far is each prime from its prediction?</h2>
    <div class="figure">
      <div class="figure-head">
        <h3 id="div-title">Proportion of fields with p | h, relative to Cohen–Lenstra</h3>
        <div class="controls">
          <div class="seg" role="group" aria-label="Family">
            <button id="div-imag" aria-pressed="true">Imaginary</button><button id="div-real" aria-pressed="false">Real</button>
          </div>
          <div class="seg" role="group" aria-label="Statistic">
            <button id="div-cum" aria-pressed="true">|D| &lt; X</button><button id="div-win" aria-pressed="false">Window at X</button>
          </div>
        </div>
      </div>
      <div class="legend" id="div-legend"></div>
      <div class="chart" id="div-chart"></div>
      <p class="caption">Each line is (observed frequency) / (CL probability) − 1. Zero means agreement. “|D| &lt; X” counts all fields below X. “Window at X” uses only the fields in [X/1.33, X), which shows the behaviour at |D| ≈ X. Sampling noise is about 10<sup>−4</sup> at the right edge, so the gaps are systematic, not random.</p>
    </div>
  </section>

  <section class="stack">
    <h2>p = 3: the gap is a theorem</h2>
    <p class="prose">Davenport and Heilbronn proved that the mean of |Cl(D)[3]| tends to 2 (imaginary) and 4/3 (real), as Cohen–Lenstra predict. Bhargava–Shankar–Tsimerman and Taniguchi–Thorne found the next term:</p>
    <div class="formula">Σ<sub>0 &lt; ±D &lt; X</sub> |Cl(D)[3]| = M<sub>±</sub>X + K<sub>±</sub>X<sup>5/6</sup> + o(X<sup>5/6</sup>), &nbsp; K<sub>±</sub> = c<sub>±</sub> · <span class="math">8ζ(1/3)</span> / (5Γ(2/3)<sup>3</sup>) · Π<sub>p</sub>(1 − (p<sup>1/3</sup>+1)/(p(p+1)))</div>
    <p class="prose">with M<sub>−</sub> = 6/π², M<sub>+</sub> = 4/π², c<sub>−</sub> = √3, c<sub>+</sub> = 1. Since ζ(1/3) &lt; 0, the correction is negative: K<sub>−</sub> = −0.4148 and K<sub>+</sub> = −0.2395. The mean therefore sits about 1.36·X<sup>−1/6</sup> below 2 for imaginary fields. A factor of 10 in X closes only a third of the gap.</p>
    <div class="figure">
      <div class="figure-head">
        <h3>Mean of |Cl(D)[3]| over |D| &lt; X</h3>
        <div class="controls"><div class="seg" role="group" aria-label="Family">
          <button id="t3-imag" aria-pressed="true">Imaginary</button><button id="t3-real" aria-pressed="false">Real</button>
        </div></div>
      </div>
      <div class="legend" id="t3-legend"></div>
      <div class="chart" id="t3-chart"></div>
      <p class="caption" id="t3-caption"></p>
    </div>
  </section>

  <section class="stack">
    <h2>How fast, prime by prime</h2>
    <div class="figure">
      <div class="figure-head">
        <h3>|mean of |Cl[p]| / CL − 1|, log-log</h3>
        <div class="controls"><div class="seg" role="group" aria-label="Family">
          <button id="sc-imag" aria-pressed="true">Imaginary</button><button id="sc-real" aria-pressed="false">Real</button>
        </div></div>
      </div>
      <div class="legend" id="sc-legend"></div>
      <div class="chart" id="sc-chart"></div>
      <p class="caption">A straight line of slope −θ means the gap shrinks like X<sup>−θ</sup>. For p = 3 the slope is the proven −1/6. For p = 5 and 7 there is no theorem, and the observed decay is faster.</p>
    </div>
  </section>

  <section class="stack">
    <h2>Which groups occur</h2>
    <div class="grid2">
      <div class="figure">
        <h3>3-Sylow subgroup, imaginary fields</h3>
        <div class="tablewrap"><table id="sylow3-imag"></table></div>
        <p class="caption">Groups with large order or rank lag furthest behind. Class numbers near |D| ≈ X are only about √X/π, so big groups are still rare at this range.</p>
      </div>
      <div class="figure">
        <h3>p-rank of Cl(D), imaginary fields</h3>
        <div class="tablewrap"><table id="ranks-imag"></table></div>
        <p class="caption">Predicted P(rank = r) = p<sup>−r²</sup> η<sub>∞</sub>(p)/η<sub>r</sub>(p)². Rank ≥ 2 needs h ≥ p², and it converges last.</p>
      </div>
    </div>
  </section>

  <section class="stack">
    <h2>Distribution by interval of |D|</h2>
    <p class="prose">Each row uses only the fields whose |D| lies in that decade, so it shows the distribution at that size rather than the running average.</p>
    <div class="figure">
      <div class="figure-head">
        <div class="controls">
          <div class="seg" role="group" aria-label="Family">
            <button id="iv-imag" aria-pressed="true">Imaginary</button><button id="iv-real" aria-pressed="false">Real</button>
          </div>
          <div class="seg" role="group" aria-label="Statistic" id="iv-groups"></div>
        </div>
      </div>
      <div class="tablewrap"><table id="iv-table"></table></div>
      <p class="caption">The bottom row is the Cohen–Lenstra prediction (Gerth's for 4-ranks). Every column moves toward it decade by decade. Per-window counts (8 windows per decade) are in <code>results/*_windows.csv</code>.</p>
    </div>
  </section>

  <section class="stack">
    <h2>Local conditions change the speed, not the limit</h2>
    <p class="prose">Restrict to fields where p splits, stays inert or ramifies in Q(√D). Cohen–Lenstra predict the same limit in each case, and for 3-torsion this is a theorem (Bhargava–Varma). The secondary terms differ, so the subfamilies sit at different distances from the limit.</p>
    <div class="figure">
      <div class="figure-head">
        <h3>Mean of |Cl(D)[p]| by behaviour of p</h3>
        <div class="controls"><div class="seg" role="group" aria-label="Family">
          <button id="lc-imag" aria-pressed="true">Imaginary</button><button id="lc-real" aria-pressed="false">Real</button>
        </div></div>
      </div>
      <div class="tablewrap"><table id="local-table"></table></div>
    </div>
  </section>

  <section class="stack">
    <h2>p = 2: compare at finite rank</h2>
    <p class="prose">Genus theory fixes the 2-rank of Cl(D) at ω(D) − 1, so Cohen–Lenstra does not apply to the 2-part itself. Gerth's version applies to the squares 2·Cl(D); it is proven for 4-ranks (Fouvry–Klüners) and for the full 2-part of imaginary fields (Smith). Because ω(D) grows like log log X, most fields here have 2-rank 1 to 3, far from the t → ∞ limit. The fair comparison is Gerth's law at the observed 2-rank t (corank of a random t×t matrix over F<sub>2</sub>), averaged over the data.</p>
    <div class="figure">
      <h3 id="r4-title">4-rank of Cl(D), imaginary fields</h3>
      <div class="legend" id="r4-legend"></div>
      <div class="chart" id="r4-chart"></div>
      <p class="caption">The finite-t prediction closes most of the gap to the data. The rest comes from slow convergence at fixed t: among fields with 2-rank 2 and |D| &lt; 10<sup>8</sup>, 44.9% have 4-rank 0, against the predicted 37.5%.</p>
    </div>
  </section>

  <section class="stack">
    <h2>Real fields: class number one</h2>
    <div class="figure">
      <h3>Share of real quadratic fields with trivial odd class group</h3>
      <div class="legend" id="h1-legend"></div>
      <div class="chart" id="h1-chart"></div>
      <p class="caption" id="h1-caption">For prime discriminants D = p the class number is odd, and Cohen–Lenstra predict h(p) = 1 with probability Π<sub>q odd</sub> Π<sub>k≥2</sub>(1 − q<sup>−k</sup>) = 75.446%.</p>
    </div>
  </section>

  <section class="stack">
    <h2>What to keep in mind</h2>
    <ol class="consider">
      <li><b>Fields, signs and ordering.</b> The heuristics are about maximal orders (fundamental discriminants), ordered by |D|. Imaginary fields have unit rank u = 0 and real fields u = 1, which changes every prediction: P(3 | h) is 43.99% versus 15.98%.</li>
      <li><b>p = 2 needs Gerth, at finite rank.</b> Compare 2·Cl(D) with the prediction at the observed 2-rank, not the limit. For real fields the theorem concerns the narrow class group, and these tables hold the ordinary one.</li>
      <li><b>Gaps are systematic.</b> With 3·10<sup>7</sup> to 3·10<sup>8</sup> fields, sampling error is about 10<sup>−4</sup>. Percent-level gaps are secondary terms, and for p = 3 they are proven (X<sup>−1/6</sup>).</li>
      <li><b>Small groups converge first.</b> h(D) ≈ √|D|·L(1,χ)/π, so p<sup>k</sup> | h needs h ≥ p<sup>k</sup>. Real fields are much worse: the regulator absorbs most of h·R. Trust a statistic only when its groups are much smaller than typical class numbers.</li>
      <li><b>Cumulative versus windowed.</b> Frequencies over |D| &lt; X mix all scales below X. Windows near X converge faster but are noisier. Both are in <code>results/*_cumulative.csv</code>.</li>
      <li><b>Moments versus probabilities.</b> The theorems are about moments, E[#Sur(Cl<sub>p</sub>, A)] = |A|<sup>−u</sup>. Higher moments are driven by rare large groups and converge more slowly: E|Cl[3]|² is 5.49 at 10<sup>8</sup>, against 6.</li>
      <li><b>Independence and local conditions.</b> Test independence of different primes against the product of the <i>observed</i> marginals. Main terms should not depend on local conditions (proven for 3-torsion by Bhargava–Varma), but secondary terms do, so subfamilies converge at different rates.</li>
      <li><b>Hypotheses.</b> Imaginary data is unconditional and was cross-checked against PARI. Real data relies on GRH via PARI's <code>quadclassunit</code>.</li>
    </ol>
  </section>

  <section class="stack">
    <h2>The data</h2>
    <p class="prose">Files live in <code>cohen-lenstra/</code> on branch <code>claude/class-group-cohen-lenstra-ekqekv</code> of <code>cocoxhuang/headshot</code>. Each line is <code>D c1,c2,…</code>, giving the invariant factors of Cl(Q(√D)) with c<sub>i+1</sub> | c<sub>i</sub>; the trivial group is <code>1</code>.</p>
    <pre>$ xzcat data/imag/imag_0000000000_0010000000.txt.xz | grep -m1 '^-3299 '
-3299 9,3</pre>
    <div class="tablewrap"><table id="data-table"></table></div>
    <p class="foot" id="foot"></p>
  </section>
</main>

<script>
const DATA = __DATA__;
const css = n => getComputedStyle(document.documentElement).getPropertyValue(n).trim();
const SUP = {"-":"⁻","0":"⁰","1":"¹","2":"²","3":"³","4":"⁴","5":"⁵","6":"⁶","7":"⁷","8":"⁸","9":"⁹"};
const pow10 = x => { const e = Math.log10(x); return Math.abs(e - Math.round(e)) < 1e-9 ? "10" + String(Math.round(e)).split("").map(c => SUP[c]).join("") : null; };
const fmtX = x => pow10(x) || x.toLocaleString("en-US");
const fmtInt = n => n.toLocaleString("en-US");
const pct = v => (v >= 0 ? "+" : "−") + Math.abs(100 * v).toFixed(2) + "%";
const fix = (v, d) => v == null ? "–" : v.toFixed(d);
const NS = "http://www.w3.org/2000/svg";
function el(tag, attrs, parent) { const e = document.createElementNS(NS, tag); for (const k in attrs) e.setAttribute(k, attrs[k]); if (parent) parent.appendChild(e); return e; }

/* Generic chart: log-x, linear or log y, lines / dots, reference lines, hover readout. */
function lineChart(host, cfg) {
  const draw = () => {
    host.innerHTML = "";
    const W = Math.max(300, host.clientWidth), H = Math.round(Math.min(380, Math.max(240, W * 0.5)));
    const m = { l: 58, r: cfg.rightPad || 18, t: 12, b: 40 };
    const svg = el("svg", { viewBox: `0 0 ${W} ${H}`, role: "img", "aria-label": cfg.aria || "" }, host);
    const xs = cfg.x, lx = v => Math.log10(v);
    const x0 = lx(xs[0]), x1 = lx(xs[xs.length - 1]);
    const vals = [];
    cfg.series.forEach(s => s.y.forEach(v => v != null && isFinite(v) && (!cfg.ylog || v > 0) && vals.push(v)));
    (cfg.refs || []).forEach(r => vals.push(r.y));
    let y0 = cfg.ymin != null ? cfg.ymin : Math.min(...vals), y1 = cfg.ymax != null ? cfg.ymax : Math.max(...vals);
    const ty = v => cfg.ylog ? Math.log10(v) : v;
    let a = ty(y0), b = ty(y1); const pad = (b - a) * 0.08 || 0.1; a -= pad; b += pad;
    const X = v => m.l + (lx(v) - x0) / (x1 - x0) * (W - m.l - m.r);
    const Y = v => H - m.b - (ty(v) - a) / (b - a) * (H - m.t - m.b);
    const grid = css("--grid"), axis = css("--axis"), ink2 = css("--ink-2"), muted = css("--muted");
    // y ticks
    let ticks = [];
    if (cfg.ylog) { for (let e = Math.floor(a); e <= Math.ceil(b); e++) { const v = Math.pow(10, e); if (ty(v) >= a && ty(v) <= b) ticks.push(v); } }
    else { const span = b - a, raw = span / 5, mag = Math.pow(10, Math.floor(Math.log10(raw))); const step = [1, 2, 2.5, 5, 10].map(k => k * mag).find(s => span / s <= 6); for (let v = Math.ceil(a / step) * step; v <= b; v += step) ticks.push(+v.toPrecision(12)); }
    ticks.forEach(v => {
      el("line", { x1: m.l, x2: W - m.r, y1: Y(v), y2: Y(v), stroke: grid, "stroke-width": 1 }, svg);
      const t = el("text", { x: m.l - 8, y: Y(v) + 4, "text-anchor": "end", "font-size": 11.5, fill: ink2 }, svg);
      t.textContent = cfg.yfmt ? cfg.yfmt(v) : v;
    });
    // x ticks: powers of ten
    for (let e = Math.ceil(x0); e <= Math.floor(x1); e++) {
      const v = Math.pow(10, e);
      el("line", { x1: X(v), x2: X(v), y1: m.t, y2: H - m.b, stroke: grid, "stroke-width": 1 }, svg);
      const t = el("text", { x: X(v), y: H - m.b + 18, "text-anchor": "middle", "font-size": 11.5, fill: ink2 }, svg);
      t.textContent = pow10(v);
    }
    el("line", { x1: m.l, x2: W - m.r, y1: H - m.b, y2: H - m.b, stroke: axis, "stroke-width": 1 }, svg);
    const xl = el("text", { x: (m.l + W - m.r) / 2, y: H - 4, "text-anchor": "middle", "font-size": 12, fill: muted }, svg);
    xl.textContent = cfg.xlabel || "X";
    // reference lines
    (cfg.refs || []).forEach(r => {
      el("line", { x1: m.l, x2: W - m.r, y1: Y(r.y), y2: Y(r.y), stroke: ink2, "stroke-width": 1.2 }, svg);
      if (r.label) { const t = el("text", { x: m.l + 6, y: Y(r.y) - 6, "font-size": 11.5, fill: ink2 }, svg); t.textContent = r.label; }
    });
    // series
    cfg.series.forEach(s => {
      const color = css(s.color);
      const pts = xs.map((x, i) => [x, s.y[i]]).filter(p => p[1] != null && isFinite(p[1]) && (!cfg.ylog || p[1] > 0));
      if (s.dots) pts.forEach(p => el("circle", { cx: X(p[0]), cy: Y(p[1]), r: 4, fill: color, stroke: css("--surface"), "stroke-width": 2 }, svg));
      else el("path", { d: pts.map((p, i) => (i ? "L" : "M") + X(p[0]).toFixed(1) + " " + Y(p[1]).toFixed(1)).join(""), fill: "none", stroke: color, "stroke-width": 2, "stroke-linejoin": "round", "stroke-linecap": "round" }, svg);
    });
    // hover layer
    const cross = el("line", { y1: m.t, y2: H - m.b, stroke: muted, "stroke-width": 1, visibility: "hidden" }, svg);
    const marks = cfg.series.map(s => el("circle", { r: 4.5, fill: css(s.color), stroke: css("--surface"), "stroke-width": 2, visibility: "hidden" }, svg));
    const tip = document.createElement("div"); tip.className = "tip"; tip.hidden = true; host.appendChild(tip);
    const hit = el("rect", { x: m.l, y: m.t, width: W - m.l - m.r, height: H - m.t - m.b, fill: "transparent", tabindex: 0 }, svg);
    let cur = xs.length - 1;
    const show = i => {
      cur = i; const x = xs[i];
      cross.setAttribute("x1", X(x)); cross.setAttribute("x2", X(x)); cross.setAttribute("visibility", "visible");
      let rows = `<div class="row" style="margin-bottom:4px"><span><b>X = ${fmtX(x)}</b></span><span>${cfg.nlabel ? cfg.nlabel(i) : ""}</span></div>`;
      cfg.series.forEach((s, k) => {
        const v = s.y[i];
        if (v == null || !isFinite(v) || (cfg.ylog && v <= 0)) { marks[k].setAttribute("visibility", "hidden"); return; }
        marks[k].setAttribute("cx", X(x)); marks[k].setAttribute("cy", Y(v)); marks[k].setAttribute("visibility", "visible");
        rows += `<div class="row"><span><i class="key-dot" style="background:${css(s.color)}"></i>${s.name}</span><span>${(s.tfmt || cfg.tfmt)(v)}</span></div>`;
      });
      tip.innerHTML = rows; tip.hidden = false;
      const px = X(x) / W * host.clientWidth, w = tip.offsetWidth;
      tip.style.left = Math.min(Math.max(0, px + 12 + w > host.clientWidth ? px - w - 12 : px + 12), host.clientWidth - w) + "px";
      tip.style.top = "8px";
    };
    const hide = () => { tip.hidden = true; cross.setAttribute("visibility", "hidden"); marks.forEach(mk => mk.setAttribute("visibility", "hidden")); };
    hit.addEventListener("pointermove", ev => {
      const r = svg.getBoundingClientRect(), px = (ev.clientX - r.left) / r.width * W;
      const lv = x0 + (px - m.l) / (W - m.l - m.r) * (x1 - x0);
      let best = 0, bd = Infinity; xs.forEach((x, i) => { const d = Math.abs(lx(x) - lv); if (d < bd) { bd = d; best = i; } });
      show(best);
    });
    hit.addEventListener("pointerleave", hide);
    hit.addEventListener("focus", () => show(cur));
    hit.addEventListener("blur", hide);
    hit.addEventListener("keydown", ev => { if (ev.key === "ArrowLeft") show(Math.max(0, cur - 1)); if (ev.key === "ArrowRight") show(Math.min(xs.length - 1, cur + 1)); });
  };
  draw();
  if (window.ResizeObserver) { let w = host.clientWidth; new ResizeObserver(() => { if (Math.abs(host.clientWidth - w) > 4) { w = host.clientWidth; draw(); } }).observe(host); }
  return draw;
}

function legend(host, items) {
  host.innerHTML = items.map(it => `<span><i class="${it.dot ? "key-dot" : "key-line"}" style="background:${css(it.color)}"></i>${it.name}</span>`).join("");
}

function segment(ids, onPick) {
  ids.forEach(id => document.getElementById(id).addEventListener("click", () => {
    ids.forEach(j => document.getElementById(j).setAttribute("aria-pressed", String(j === id)));
    onPick(id);
  }));
}

/* header chips */
(function () {
  const I = DATA.imag, R = DATA.real;
  const chips = [];
  if (I) chips.push(`<span class="chip"><b>${fmtInt(I.Nmax)}</b> imaginary fields, |D| &lt; ${fmtX(I.Xmax)} · unconditional</span>`);
  if (R) chips.push(`<span class="chip"><b>${fmtInt(R.Nmax)}</b> real fields, D &lt; ${fmtX(R.Xmax)} · GRH</span>`);
  if (I && DATA.imag_raw_max) chips.push(`<span class="chip">raw class groups stored to |D| &lt; <b>${fmtX(DATA.imag_raw_max)}</b></span>`);
  document.getElementById("chips").innerHTML = chips.join("");
  if (I && R) document.getElementById("lede").innerHTML = `Every imaginary quadratic field with |D| &lt; ${fmtX(I.Xmax)} and every real quadratic field with D &lt; ${fmtX(R.Xmax)}, compared with the Cohen–Lenstra predictions. The data does approach the predictions, but slowly. For p = 3 the remaining gap matches the proven secondary term.`;
})();

/* figure 1: divisibility */
(function () {
  let fam = "imag", mode = "cum";
  const host = document.getElementById("div-chart");
  const render = () => {
    const F = DATA[fam]; if (!F) return;
    const ps = [[3, "--s1"], [5, "--s2"], [7, "--s3"]];
    const series = ps.map(([p, c]) => { const s = F.series["p" + p + "_div"]; return { name: "p = " + p, color: c, y: s[mode].map(v => v == null ? null : v / s.cl - 1) }; });
    legend(document.getElementById("div-legend"), ps.map(([p, c]) => ({ name: `p = ${p} (CL ${F.series["p" + p + "_div"].cl.toFixed(4)})`, color: c })));
    const start = F.X.findIndex(x => x >= (fam === "imag" ? 3000 : 10000));
    const sl = a => a.slice(start);
    lineChart(host, { x: sl(F.X), series: series.map(s => ({ ...s, y: sl(s.y) })), refs: [{ y: 0, label: "Cohen–Lenstra" }],
      yfmt: v => { const r = Math.round(v * 100); return (r > 0 ? "+" : r < 0 ? "−" : "") + Math.abs(r) + "%"; }, tfmt: pct,
      nlabel: i => fmtInt(sl(F.N)[i]) + " fields", xlabel: mode === "cum" ? "X (fields with |D| < X)" : "X (fields with X/1.33 ≤ |D| < X)",
      aria: "Relative deviation of P(p divides h) from Cohen–Lenstra for p = 3, 5, 7" });
  };
  segment(["div-imag", "div-real"], id => { fam = id === "div-imag" ? "imag" : "real"; render(); });
  segment(["div-cum", "div-win"], id => { mode = id === "div-cum" ? "cum" : "win"; render(); });
  render();
})();

/* figure 2: 3-torsion mean */
(function () {
  let fam = "imag";
  const host = document.getElementById("t3-chart");
  const render = () => {
    const F = DATA[fam]; if (!F) return;
    const s = F.series["p3_sum_p^rank"];
    const start = F.X.findIndex(x => x >= 1000), sl = a => a.slice(start);
    legend(document.getElementById("t3-legend"), [{ name: "data", color: "--s1", dot: true }, { name: "main term + X⁵ᐟ⁶ term", color: "--s2" }]);
    lineChart(host, { x: sl(F.X), series: [{ name: "main + second term", color: "--s2", y: sl(s.second) }, { name: "data", color: "--s1", y: sl(s.cum), dots: true }],
      refs: [{ y: s.cl, label: fam === "imag" ? "limit 2" : "limit 4/3" }], yfmt: v => v.toFixed(2), tfmt: v => v.toFixed(5),
      nlabel: i => fmtInt(sl(F.N)[i]) + " fields", aria: "Mean of the 3-torsion of class groups against the two-term asymptotic" });
    const n = s.cum.length - 1;
    document.getElementById("t3-caption").textContent = `At X = ${fmtX(F.X[n])}: observed mean ${s.cum[n].toFixed(6)}, two-term prediction ${s.second[n].toFixed(6)}, limit ${fam === "imag" ? "2" : "4/3"}. The two-term formula follows the data at every scale.`;
  };
  segment(["t3-imag", "t3-real"], id => { fam = id === "t3-imag" ? "imag" : "real"; render(); });
  render();
})();

/* figure 3: deviation scaling */
(function () {
  let fam = "imag";
  const host = document.getElementById("sc-chart");
  const slope = (xs, ys) => { const p = xs.map((x, i) => [Math.log(x), Math.log(ys[i])]).filter(q => isFinite(q[1])); const n = p.length; const mx = p.reduce((s, q) => s + q[0], 0) / n, my = p.reduce((s, q) => s + q[1], 0) / n; return p.reduce((s, q) => s + (q[0] - mx) * (q[1] - my), 0) / p.reduce((s, q) => s + (q[0] - mx) ** 2, 0); };
  const render = () => {
    const F = DATA[fam]; if (!F) return;
    const ps = [[3, "--s1"], [5, "--s2"], [7, "--s3"]];
    const start = F.X.findIndex((x, i) => F.N[i] >= 5000), sl = a => a.slice(start);
    const xs = sl(F.X);
    const series = ps.map(([p, c]) => { const s = F.series["p" + p + "_sum_p^rank"]; return { p, name: "p = " + p, color: c, y: sl(s.cum).map(v => Math.abs(v / s.cl - 1)) }; });
    legend(document.getElementById("sc-legend"), series.map(s => {
      const k = xs.findIndex(x => x >= xs[xs.length - 1] / 100);
      return { name: `${s.name} (slope over last 2 decades ${slope(xs.slice(k), s.y.slice(k)).toFixed(2)})`, color: s.color };
    }));
    lineChart(host, { x: xs, series, ylog: true, yfmt: v => pow10(v) || v, tfmt: v => (100 * v).toFixed(3) + "%", nlabel: i => fmtInt(sl(F.N)[i]) + " fields",
      aria: "Relative gap of the mean p-torsion to Cohen–Lenstra on log-log axes" });
  };
  segment(["sc-imag", "sc-real"], id => { fam = id === "sc-imag" ? "imag" : "real"; render(); });
  render();
})();

/* tables */
(function () {
  const I = DATA.imag; if (!I) return;
  const t = document.getElementById("sylow3-imag");
  t.innerHTML = `<thead><tr><th>Cl[3<sup>∞</sup>]</th><th>CL</th><th>data, |D| &lt; ${fmtX(I.Xmax)}</th><th>data / CL − 1</th></tr></thead><tbody>` +
    I.sylow3.filter(r => r.cl > 2e-5).map(r => `<tr><td>${r.name}</td><td>${r.cl.toFixed(5)}</td><td>${r.emp.toFixed(5)}</td><td class="dev">${pct(r.emp / r.cl - 1)}</td></tr>`).join("") + "</tbody>";
  const rt = document.getElementById("ranks-imag");
  let html = `<thead><tr><th>p</th><th>rank</th><th>CL</th><th>data</th><th>data / CL − 1</th></tr></thead><tbody>`;
  [3, 5, 7].forEach(p => I.ranks[p].slice(0, 3).forEach(r => { html += `<tr><td>${r.r === "0" ? p : ""}</td><td>${r.r}</td><td>${r.cl.toFixed(5)}</td><td>${r.emp.toFixed(5)}</td><td class="dev">${pct(r.emp / r.cl - 1)}</td></tr>`; }));
  rt.innerHTML = html + "</tbody>";
})();

/* interval tables */
(function () {
  let fam = "imag", group = "div";
  const holder = document.getElementById("iv-groups");
  const labels = { div: "p | h", rank3: "3-rank", syl3: "3-Sylow", rank4: "4-rank", h1: "h = 1" };
  const render = () => {
    const F = DATA[fam], t = document.getElementById("iv-table");
    if (!F || !F.iv) { t.innerHTML = ""; holder.innerHTML = ""; return; }
    const keys = Object.keys(F.iv.groups);
    if (!keys.includes(group)) group = keys[0];
    holder.innerHTML = keys.map(k => `<button data-g="${k}" aria-pressed="${k === group}">${labels[k]}</button>`).join("");
    holder.querySelectorAll("button").forEach(b => b.addEventListener("click", () => { group = b.dataset.g; render(); }));
    const G = F.iv.groups[group];
    const iv = F.iv.intervals;
    const lab = (lo, hi) => `[${fmtX(lo)}, ${fmtX(hi)})`;
    t.innerHTML = `<thead><tr><th>|D| in</th><th>fields</th>${G.cols.map(c => `<th>${c}</th>`).join("")}</tr></thead><tbody>` +
      iv.map((r, i) => `<tr><td>${lab(r[0], r[1])}</td><td>${fmtInt(r[2])}</td>${G.rows[i].map(v => `<td>${v == null ? "–" : v.toFixed(5)}</td>`).join("")}</tr>`).join("") +
      `<tr><td><b>Cohen–Lenstra</b></td><td></td>${G.cl.map(v => `<td><b>${v == null ? "–" : v.toFixed(5)}</b></td>`).join("")}</tr></tbody>`;
  };
  segment(["iv-imag", "iv-real"], id => { fam = id === "iv-imag" ? "imag" : "real"; render(); });
  render();
})();

/* local conditions table */
(function () {
  const render = fam => {
    const F = DATA[fam], t = document.getElementById("local-table");
    if (!F || !F.local) { t.innerHTML = ""; return; }
    const name = { split: "splits", inert: "inert", ram: "ramifies" };
    t.innerHTML = `<thead><tr><th>p in Q(√D)</th><th>fields</th><th>mean |Cl[p]|</th><th>CL</th><th>gap</th><th>P(p | h)</th><th>gap</th></tr></thead><tbody>` +
      F.local.map(r => `<tr><td>p = ${r.p} ${name[r.kind]}</td><td>${fmtInt(r.N)}</td><td>${r.mean.toFixed(5)}</td><td>${r.cl.toFixed(4)}</td><td class="dev">${pct(r.mean / r.cl - 1)}</td><td>${r.div.toFixed(5)}</td><td class="dev">${pct(r.div / r.cl_div - 1)}</td></tr>`).join("") + "</tbody>";
  };
  segment(["lc-imag", "lc-real"], id => render(id === "lc-imag" ? "imag" : "real"));
  render("imag");
})();

/* figure 4: 4-rank bars */
(function () {
  const I = DATA.imag; if (!I || !I.rank4) return;
  const host = document.getElementById("r4-chart");
  document.getElementById("r4-title").textContent = `4-rank of Cl(D), imaginary fields with |D| < ${fmtX(I.Xmax)}`;
  const groups = [["data", "emp", "--s1"], ["Gerth at the observed 2-ranks", "mix", "--s2"], ["Gerth limit (t → ∞)", "lim", "--s3"]];
  legend(document.getElementById("r4-legend"), groups.map(g => ({ name: g[0], color: g[2], dot: true })));
  const draw = () => {
    host.innerHTML = "";
    const W = Math.max(300, host.clientWidth), H = 260, m = { l: 46, r: 10, t: 16, b: 34 };
    const svg = el("svg", { viewBox: `0 0 ${W} ${H}`, role: "img", "aria-label": "4-rank distribution: data, finite-rank Gerth prediction and limit" }, host);
    const ymax = 0.65, Y = v => H - m.b - v / ymax * (H - m.t - m.b);
    [0, 0.2, 0.4, 0.6].forEach(v => { el("line", { x1: m.l, x2: W - m.r, y1: Y(v), y2: Y(v), stroke: css("--grid") }, svg); const t = el("text", { x: m.l - 8, y: Y(v) + 4, "text-anchor": "end", "font-size": 11.5, fill: css("--ink-2") }, svg); t.textContent = v.toFixed(1); });
    const band = (W - m.l - m.r) / I.rank4.length, bw = Math.min(24, band / 4.2);
    I.rank4.forEach((r, i) => {
      const cx = m.l + band * (i + 0.5);
      groups.forEach((g, k) => {
        const v = r[g[1]], x = cx + (k - 1) * (bw + 4) - bw / 2, y = Y(v), h = Y(0) - y;
        const rr = Math.min(4, h / 2);
        el("path", { d: `M${x} ${Y(0)}V${y + rr}Q${x} ${y} ${x + rr} ${y}H${x + bw - rr}Q${x + bw} ${y} ${x + bw} ${y + rr}V${Y(0)}Z`, fill: css(g[2]) }, svg);
        const t = el("text", { x: x + bw / 2, y: y - 5, "text-anchor": "middle", "font-size": 10.5, fill: css("--ink-2") }, svg); t.textContent = v < 0.0005 ? "<.001" : v.toFixed(3).replace(/^0/, "");
      });
      const t = el("text", { x: cx, y: H - m.b + 18, "text-anchor": "middle", "font-size": 12, fill: css("--ink-2") }, svg); t.textContent = "4-rank " + r.r;
    });
    el("line", { x1: m.l, x2: W - m.r, y1: Y(0), y2: Y(0), stroke: css("--axis") }, svg);
  };
  draw();
  if (window.ResizeObserver) { let w = host.clientWidth; new ResizeObserver(() => { if (Math.abs(host.clientWidth - w) > 4) { w = host.clientWidth; draw(); } }).observe(host); }
})();

/* figure 5: real h = 1 */
(function () {
  const R = DATA.real; if (!R) return;
  const host = document.getElementById("h1-chart");
  const a = R.series["prime_h_eq_1"], b = R.series["odd_part_trivial"];
  const start = R.X.findIndex(x => x >= 1000), sl = v => v.slice(start);
  legend(document.getElementById("h1-legend"), [{ name: "h = 1 among prime discriminants D = p", color: "--s1" }, { name: "odd part of h trivial, all D > 0", color: "--s2" }]);
  lineChart(host, { x: sl(R.X), series: [{ name: "h(p) = 1", color: "--s1", y: sl(a.cum) }, { name: "odd part trivial", color: "--s2", y: sl(b.cum) }],
    refs: [{ y: a.cl, label: "Cohen–Lenstra 75.446%" }], yfmt: v => (100 * v).toFixed(0) + "%", tfmt: v => (100 * v).toFixed(3) + "%",
    nlabel: i => fmtInt(sl(R.N)[i]) + " fields", xlabel: "X (fields with D < X)", aria: "Share of real quadratic fields with class number one" });
  const n = a.cum.length - 1, v = a.cum[n], np = a.n[n];
  document.getElementById("h1-caption").textContent += ` Among the ${fmtInt(np)} prime discriminants below ${fmtX(R.X[n])}, ${(100 * v).toFixed(3)}% have h = 1: ${v > a.cl ? "still above the prediction, because small class groups dominate at this range, and drifting down toward it" : "below the prediction"}.`;
})();

/* data table and footer */
(function () {
  const t = document.getElementById("data-table");
  t.innerHTML = `<thead><tr><th>Path</th><th style="text-align:left">Contents</th></tr></thead><tbody>
    <tr><td><code>data/imag/</code></td><td style="text-align:left">all imaginary fields, |D| &lt; ${fmtX(DATA.imag_raw_max || 1e8)}, one xz file per 10⁷ discriminants</td></tr>
    <tr><td><code>data/real/</code></td><td style="text-align:left">all real fields, D &lt; ${DATA.real ? fmtX(DATA.real.Xmax) : "10⁸"}, one xz file per 10⁷ discriminants</td></tr>
    <tr><td><code>results/</code></td><td style="text-align:left">summaries, per-window counts, fits, plots</td></tr>
    <tr><td><code>src/</code></td><td style="text-align:left">computation, verification and analysis code</td></tr></tbody>`;
  document.getElementById("foot").textContent = "Imaginary class groups: reduced-form counting plus explicit p-Sylow construction (no GRH), identical to PARI on all 18,237,779 fields with |D| < 6·10⁷. Real class groups: PARI/GP 2.15 quadclassunit (GRH). Generated " + DATA.generated + ".";
})();
</script>
"""


def main():
    results, out = sys.argv[1], sys.argv[2]
    payload = {"generated": __import__("datetime").date.today().isoformat(), "imag_raw_max": 10 ** 8}
    for fam in ("imag", "real"):
        d = load(results, fam)
        payload[fam] = family_payload(d, fam) if d else None
        if payload[fam]:
            payload[fam]["iv"] = interval_payload(results, fam)
    html = TEMPLATE.replace("__DATA__", json.dumps(payload, separators=(",", ":")))
    with open(out, "w") as f:
        f.write(html)
    print("wrote", out, len(html), "bytes")


if __name__ == "__main__":
    main()
