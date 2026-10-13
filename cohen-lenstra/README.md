# Class groups of quadratic fields vs. Cohen–Lenstra

Class groups of every quadratic field up to a large bound, and the tools to
compare them with the Cohen–Lenstra (CL) heuristics and to watch how slowly
they converge.


| Path       | Contents                                             |
| ---------- | ---------------------------------------------------- |
| `results/` | statistics vs. CL, with imaginary fields up to       |
| `src/`     | code to recompute, check and extend all of the above |


The LMFDB's complete list of quadratic fields stops at |D| < 2·10⁶. For  
imaginary fields it also hosts the Mosunov–Jacobson class groups for all  
|D| < 2⁴⁰, far beyond what is stored here.

## How it was computed

- **Imaginary fields, unconditional.** `src/clgrp_imag.c` gets every h(D) at
once by counting reduced forms, then builds each p-Sylow subgroup with p² | h
explicitly (for p = 2, genus theory leaves only the squares to build).
- **Real fields, under GRH.** PARI 2.15 `quadclassunit` via `src/clgrp_pari.py`.
- **Checks.** For |D| < 6·10⁷ the two methods agree on every class group. A  
second form counter (`src/countforms.c`) reproduces every imaginary class  
number. Real class groups agree with PARI's `qfbclassno` and `bnfinit` on  
random samples, including above 2³² and 10¹⁰.



## Reproducing or extending

Needs cypari (or cypari2), numpy, scipy, matplotlib and gcc with OpenMP.

```bash
cd cohen-lenstra
src/run_imag_fast.sh work/imag 0 1000000000 10000000  # imaginary |D| < 1e9 (~1 h on 4 cores)
src/run_all.sh work XMAX 1000000 4 "1" XMIN           # real XMIN <= D < XMAX; resumable
```

