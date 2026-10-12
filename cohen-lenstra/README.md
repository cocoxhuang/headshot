# Class groups of quadratic fields vs. Cohen–Lenstra

Complete tables of class groups of quadratic fields, plus the tools to compare
them with the Cohen–Lenstra heuristics and to watch (how slowly) they converge.

| Path | Contents |
|---|---|
| `data/imag/` | Cl(D) for **all 30,396,385 imaginary quadratic fields with \|D\| < 10⁸** (unconditional) |
| `data/real/` | Cl(D) for **all 3,234,172,229 real quadratic fields with D < 1.064·10¹⁰** (under GRH) |
| `results/` | statistics vs. Cohen–Lenstra, imaginary fields up to **\|D\| < 10⁹ (303,963,510 fields)**: `*_summary.md` (at the largest X), `*_intervals.md` (distributions decade by decade), `*_cumulative.csv`, `*_windows.csv` (8 windows per decade), `extrapolation.md`, `plots/`, `report.html` |
| `src/` | everything needed to recompute, check or extend all of the above |

**Compared with the LMFDB.** The LMFDB's complete list of quadratic fields stops
at |D| < 2·10⁶, so the real-field table here is far larger. For imaginary
fields the LMFDB also hosts the unconditional Mosunov–Jacobson class group
tables for all |D| < 2⁴⁰ ≈ 1.1·10¹² (about 2.2 TB,
<https://www.lmfdb.org/NumberField/QuadraticImaginaryClassGroups>), far beyond
what is stored here. Use those to push imaginary statistics further. The raw
imaginary class groups for 10⁸ ≤ |D| < 10⁹ (about 1 GB compressed) are not
stored here; `src/run_imag_fast.sh` regenerates them in under an hour on 4 cores.

## Data format

One line per field, in increasing |D|:

```
D c1,c2,...,ck
```

`D` is the fundamental discriminant and `c1,...,ck` are the invariant factors of
the class group of Q(√D): Cl(D) ≅ Z/c1 × ⋯ × Z/ck with c(i+1) | c(i), so
h(D) = c1⋯ck. The trivial group is written `1`. Every fundamental discriminant in
the range occurs exactly once (D = 1 is excluded). Files cover 10⁷ consecutive
values of |D| and are xz-compressed:

```console
$ xzcat data/imag/imag_0000000000_0010000000.txt.xz | head -4
-3 1
-4 1
-7 1
-8 1
$ xzcat data/imag/imag_0000000000_0010000000.txt.xz | grep -m1 '^-3299 '
-3299 9,3
```

```python
import lzma
with lzma.open("data/real/real_0000000000_0010000000.txt.xz", "rt") as f:
    for line in f:
        D, cyc = line.split()
        cyc = [int(c) for c in cyc.split(",")]      # [1] means trivial
```

## How the data was computed, and how far to trust it

* **Imaginary fields: unconditional.** `src/clgrp_imag.c` computes h(D) for every
  |D| < X at once by counting reduced forms (a, b, c), |b| ≤ a ≤ c (about
  (π/18)·X^{3/2} additions in total, fewer after skipping forms whose
  discriminant cannot be fundamental). For each odd prime p with p² | h it then
  builds the p-Sylow subgroup explicitly. It composes powers of small prime forms
  until the subgroup has order p^{v_p(h)}, which is known exactly, and reads the
  isomorphism type off the element orders. For p = 2, genus theory gives the
  2-rank ω(D) − 1, so only the subgroup of squares has to be built. No
  hypothesis is used anywhere. Cost: about 15 µs per field around 10⁸ and
  35 µs around 10⁹ on one core.
* **Real fields: GRH.** PARI/GP 2.15 `quadclassunit` (Buchmann's subexponential
  algorithm) via `src/clgrp_pari.py`. Its output is certified under GRH.
  Regulators are not stored.
* **Cross-checks.** For |D| < 6·10⁷ (18,237,779 fields), PARI and the
  unconditional program were run independently and agree on every class group,
  byte for byte. Separately, `src/countforms.c` (a second, simpler form counter)
  reproduces every imaginary class number. Random samples near 10⁹ agree with
  PARI. For real fields, 2,000 randomly chosen class numbers agree with PARI's
  `qfbclassno`, which uses a different algorithm (forms and baby-step
  giant-step).

## Results at a glance

Cumulative frequencies over all fields with |D| < X, against the Cohen–Lenstra
(CL) prediction. Full tables: `results/imag_summary.md`, `results/real_summary.md`.

**Imaginary quadratic fields** (u = 0)

| statistic | X = 10⁴ | 10⁵ | 10⁶ | 10⁷ | 10⁸ | 10⁹ | CL |
|---|---|---|---|---|---|---|---|
| P(3 \| h) | 0.3421 | 0.3784 | 0.4002 | 0.4137 | 0.4224 | 0.4281 | 0.4399 |
| P(5 \| h) | 0.2028 | 0.2213 | 0.2295 | 0.2347 | 0.2368 | 0.2381 | 0.2397 |
| P(7 \| h) | 0.1351 | 0.1482 | 0.1569 | 0.1605 | 0.1620 | 0.1626 | 0.1632 |
| mean of \|Cl[3]\| | 1.6980 | 1.7990 | 1.8633 | 1.9072 | 1.9366 | 1.9568 | 2 |
| P(3-rank = 2) | 0.0023 | 0.0070 | 0.0105 | 0.0133 | 0.0152 | 0.0166 | 0.0197 |
| P(odd part of Cl cyclic) | 0.9977 | 0.9916 | 0.9878 | 0.9845 | 0.9824 | 0.9808 | 0.9776 |
| P(4-rank = 1) | 0.3687 | 0.4105 | 0.4342 | 0.4495 | 0.4602 | 0.4681 | 0.5776 (Gerth) |

**Real quadratic fields** (u = 1)

| statistic | X = 10⁴ | 10⁵ | 10⁶ | 10⁷ | 10⁸ | 10⁹ | 10¹⁰ | 1.064·10¹⁰ | CL |
|---|---|---|---|---|---|---|---|---|---|
| P(3 \| h) | 0.0815 | 0.1087 | 0.1253 | 0.1370 | 0.1445 | 0.1494 | 0.1528 | 0.1529 | 0.1598 |
| P(5 \| h) | 0.0230 | 0.0353 | 0.0413 | 0.0454 | 0.0474 | 0.0484 | 0.0490 | 0.0490 | 0.0496 |
| P(7 \| h) | 0.0095 | 0.0142 | 0.0188 | 0.0213 | 0.0226 | 0.0232 | 0.0235 | 0.0235 | 0.0237 |
| mean of \|Cl[3]\| | 1.1630 | 1.2185 | 1.2537 | 1.2796 | 1.2968 | 1.3084 | 1.3164 | 1.3165 | 1.3333 |
| P(odd part of h trivial) | 0.8847 | 0.8353 | 0.8045 | 0.7842 | 0.7727 | 0.7659 | 0.7618 | 0.7617 | 0.7545 |
| P(h = 1), prime discriminants D = p | 0.8459 | 0.7956 | 0.7840 | 0.7717 | 0.7655 | 0.7618 | 0.7594 | 0.7593 | 0.7545 |

![p | h](results/plots/divisibility_convergence.png)

![3-torsion](results/plots/three_torsion_mean.png)

More plots: `results/plots/` (deviation scaling, Sylow-type distributions,
4-ranks, class number one for real fields).

## What to consider

**1. Compare like with like.**
- *Fields, not orders.* The heuristics are about class groups of maximal orders,
  i.e. fundamental discriminants (all of them are here). Non-maximal orders behave
  differently.
- *Keep the two signs apart.* For odd p the predicted law of the p-part is
  P(Cl_p ≅ G) ∝ 1/(|G|^u·|Aut G|), with unit rank u = 0 (imaginary) or 1 (real).
  For example P(3 | h) = 43.99% for imaginary fields but 15.98% for real ones.
- *Fix the ordering.* "Converges" means: over the fields with |D| < X, as X → ∞.
  Other orderings, such as real fields ordered by regulator, are different questions.

**2. p = 2 is not Cohen–Lenstra; use Gerth's version, at finite rank.**
Genus theory forces Cl(D)[2] ≅ (Z/2)^{ω(D)−1} (imaginary fields; for real
fields this holds for the narrow class group). So the 2-part is determined by
the number of prime factors, not random. Gerth's extension predicts that the
group of squares 2·Cl(D) follows the CL law with p = 2. This is proved for the
4-rank (Fouvry–Klüners) and for the whole 2^∞-part of imaginary fields (Smith).
In practice the limit is invisible: ω(D) grows like log log X, so the 2-rank is
mostly 1–3 for |D| < 10⁸. Even for fixed 2-rank t the 4-rank approaches its
finite-t law (the corank of a random t×t matrix over F₂) only slowly. So compare
with the *finite-t* prediction averaged over the observed 2-ranks (the "Gerth
mixture" rows in `imag_summary.md`, and `plots/four_rank.png`). That explains
most of the gap (4-rank 0 at |D| < 10⁹: data 0.488, finite-t mixture 0.444,
limit 0.289).
For real fields the theorem is about the *narrow* class group. These tables
contain the ordinary class group, so do not compare their 2-parts with Gerth
directly (the difference is tied to the negative Pell equation).

**3. Convergence is slow, and for p = 3 we know exactly why.**
With ≈3·10⁷ fields the sampling error of a frequency is ≈10⁻⁴, so the
percent-level gaps in the table are systematic, not noise. For p = 3 they are a
theorem. Davenport–Heilbronn give the main term, and Bhargava–Shankar–Tsimerman
and Taniguchi–Thorne give the second term:

  Σ_{0 < ±D < X} |Cl(D)[3]| = M± X + K± X^{5/6} + o(X^{5/6}),

with M₋ = 6/π², M₊ = 4/π², and
K± = c± · 8ζ(1/3)/(5Γ(2/3)³) · ∏_p (1 − (p^{1/3}+1)/(p(p+1))), c₋ = √3, c₊ = 1.
Numerically K₋ = −0.41478 and K₊ = −0.23947, since ζ(1/3) < 0. So the mean of
|Cl[3]| is about 2 − 1.36·X^{−1/6} for imaginary fields and 4/3 − 0.79·X^{−1/6}
for real ones. Each factor of 10 in X shrinks the gap by only a third. The data
follows main term + second term to within 2·10⁻⁵ at every scale: 1.936642
observed vs. 1.936660 predicted at X = 10⁸, and 1.956836 vs. 1.956849 at 10⁹
(plot above; decade by decade in `results/imag_intervals.md`). So the deficit is the predicted secondary
term, not evidence against Cohen–Lenstra. It is worth computing the second term
alongside any p = 3 statistic.
For p ≥ 5 there is no theorem. Empirically the relative gaps of the p-torsion
means shrink faster than X^{−1/6}: log-log slopes over the last two decades are
in `plots/deviation_scaling.png`. Fits of F(X) = F∞ + a·X^{−θ} extrapolate to
limits close to CL (`results/extrapolation.md`). Treat such extrapolations with
care: the same fits put real h(p) = 1 at 0.751–0.755, slightly *below* 75.446%,
although the data is still above it at 1.064·10¹⁰.

**4. Finite-size effects: class groups are small.**
For imaginary fields h(D) ≈ √|D|·L(1,χ_D)/π, so anything that needs a large
group is suppressed at accessible |D|. For example p^k | h needs h ≥ p^k, and
p-rank 2 needs h ≥ p². At |D| < 10⁹, P(19² | h) is still 12% below CL and
P(3-rank ≥ 3) is 56% below CL, while P(7 | h) is within 0.4%. Rule of thumb:
trust a statistic only when the groups it is about are much smaller than the
typical class number in your range. For real fields this is much worse:
h·R ≈ √D·L(1,χ_D)/2 and the regulator R is usually of size about √D, so class
numbers are usually tiny.

**5. Cumulative vs. windowed statistics.**
Cumulative frequencies (|D| < X) average over all scales below X. Windowed
frequencies (fields in [X/1.33, X), 8 windows per decade) show the behaviour at
|D| ≈ X and converge faster, at the price of more noise. Both are in
`*_cumulative.csv` (`cum_freq`, `win_freq`, with `win_N` fields per window).
Raw counts per window are in `*_windows.csv`, so any other binning can be
recomputed without the raw data.

**6. Moments vs. probabilities.**
The theorems for p = 3 (Davenport–Heilbronn) and for 4-ranks (Fouvry–Klüners)
are about *moments*: CL is equivalent to E[#Sur(Cl_p, A)] = |A|^{−u} for every finite
abelian p-group A. For example E|Cl[p]| = 1 + p^{−u} and, for imaginary fields,
E|Cl[p]|² = p + 3. Higher moments are dominated by rare large groups and
converge more slowly: the observed E|Cl[3]|² at 10⁹ is 5.65, against 6. The most
frequent events (p ∤ h, p-rank 1) converge fastest.

**7. Independence across primes, and local conditions.**
CL predicts that the p-parts for different p are independent. In the data,
P(3 | h and 5 | h) = 0.10179 at 10⁹ vs. 0.10194 for the product of the *observed*
marginals: compare with the observed marginals, not the CL values. The odd part
is also expected not to care about local conditions (D mod 8, how p splits,
...). Bhargava–Varma proved this for the mean 3-torsion. But secondary terms do
depend on local conditions (Taniguchi–Thorne find biases in arithmetic
progressions), so subfamilies converge at different speeds. For example, at
|D| < 10⁹ the mean of |Cl[3]| over imaginary fields is 1.9583 when 3 splits,
1.9583 when 3 is inert and 1.9524 when 3 ramifies (limit 2 in each case). The
`@3split`, `@3inert`, `@3ram` (and `@5…`) rows of the summaries have these splits.

**8. Real quadratic fields.**
For prime discriminants D = p ≡ 1 (mod 4) the class number is odd, and CL
predicts h(p) = 1 with probability ∏_{q odd} ∏_{k≥2}(1 − q^{−k}) = 75.446%.
Among the 241,401,087 prime discriminants p < 1.064·10¹⁰, 75.93% have h(p) = 1. That
is still above the prediction, but drifting down toward it interval by interval
(0.788, 0.782, 0.770, 0.765 for the decades 10⁴–10⁸, then 0.761 for [10⁸, 10⁹), then 0.759 for [10⁹, 10¹⁰), then 0.758 for [10¹⁰, 1.064·10¹⁰);
`results/real_intervals.md`, `plots/real_class_number_one.png`). te Riele and
Williams (2003) pushed this test much further.

**9. Larger computations exist.**
Jacobson–Ramachandran–Williams tabulated imaginary class groups to |D| < 10¹¹
(under ERH), and Mosunov–Jacobson did so unconditionally to |D| < 2⁴⁰. Their
tables are a good benchmark for the trends seen here, and for how far the
extrapolations in `results/extrapolation.md` can be trusted.

**10. Beyond quadratic fields.**
Cohen–Martinet extend the heuristics to other Galois groups. Malle observed that
they must be modified when the base field contains p-th roots of unity (the
p = 2 phenomenon above is the simplest case). Over function fields, Cohen–Lenstra
is a theorem in a large-q limit (Ellenberg–Venkatesh–Westerland).

## Reproducing or extending

Requirements: PARI via `pip install cypari` (or cypari2), numpy, scipy,
matplotlib, gcc with OpenMP.

```bash
cd cohen-lenstra
# imaginary fields, unconditional (about 1 h for |D| < 1e9 on 4 cores)
src/run_imag_fast.sh work/imag 0 1000000000 10000000
# real fields with PARI (GRH): about 1 h for D < 1e8 on 4 cores; resumable
src/run_all.sh work 100000000 1000000 4 "1"
# continue the real fields beyond the stored 1.7e8 (about 1 h per 1e8 at first)
src/run_all.sh work 1000000000 1000000 4 "1" 170000000
# independent check of imaginary class numbers
gcc -O3 -march=native -fopenmp -o countforms src/countforms.c -lm
./countforms 100000000 | python3 src/verify_classnumbers.py data/imag/*.xz
# statistics (count in parallel on disjoint files, then merge), fits, plots
python3 src/analyze.py count imag work/c1.csv <files...>
python3 src/analyze.py report imag results work/c*.csv
python3 src/extrapolate.py results 1e5
python3 src/plots.py results
python3 src/cl_predictions.py          # table of CL predictions
```

`src/cl_predictions.py` has the exact CL probabilities (p-rank laws, Sylow-type
laws via |Aut G|, moments, cyclicity of the odd part, ...) for u = 0, 1 and any p.

**Going further.** The form count grows like X^{3/2}, so |D| < 10¹⁰ costs about
30 times the counting work of 10⁹. That is roughly a day on 4 cores, or an hour
on a ~100-core cluster: split the range with `run_imag_fast.sh` and merge the
`analyze.py count` outputs. Past that, compute class numbers with theta-series/FFT
methods (Mosunov–Jacobson). Real fields cost more per field, because PARI must
also find the regulator.

## References

- H. Cohen, H. W. Lenstra Jr., *Heuristics on class groups of number fields*, Number Theory Noordwijkerhout 1983, LNM 1068 (1984).
- H. Davenport, H. Heilbronn, *On the density of discriminants of cubic fields II*, Proc. R. Soc. A 322 (1971).
- M. Bhargava, A. Shankar, J. Tsimerman, *On the Davenport–Heilbronn theorems and second order terms*, Invent. Math. 193 (2013).
- T. Taniguchi, F. Thorne, *Secondary terms in counting functions for cubic fields*, Duke Math. J. 162 (2013).
- F. Gerth III, *Extension of conjectures of Cohen and Lenstra*, Expo. Math. 5 (1987).
- É. Fouvry, J. Klüners, *On the 4-rank of class groups of quadratic number fields*, Invent. Math. 167 (2007).
- A. Smith, *2^∞-Selmer groups, 2^∞-class groups, and Goldfeld's conjecture*, arXiv:1702.02325.
- M. Bhargava, I. Varma, *The mean number of 3-torsion elements in the class groups and ideal groups of quadratic orders*, Proc. LMS 112 (2016).
- M. J. Jacobson Jr., S. Ramachandran, H. C. Williams, *Numerical results on class groups of imaginary quadratic fields*, ANTS VII, LNCS 4076 (2006).
- A. S. Mosunov, M. J. Jacobson Jr., *Unconditional class group tabulation of imaginary quadratic fields to \|Δ\| < 2⁴⁰*, Math. Comp. 85 (2016).
- H. J. J. te Riele, H. C. Williams, *New computations concerning the Cohen–Lenstra heuristics*, Experiment. Math. 12 (2003).
- G. Malle, *Cohen–Lenstra heuristic and roots of unity*, J. Number Theory 128 (2008).
- J. Ellenberg, A. Venkatesh, C. Westerland, *Homological stability for Hurwitz spaces and the Cohen–Lenstra conjecture over function fields*, Ann. of Math. 183 (2016).
- C. Hillar, D. Rhea, *Automorphisms of finite abelian groups*, Amer. Math. Monthly 114 (2007).
