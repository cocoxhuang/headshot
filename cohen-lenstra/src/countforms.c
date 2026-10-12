/*
 * countforms.c -- unconditional class numbers h(D) of all imaginary quadratic
 * fields with fundamental discriminant D, 3 <= |D| < X, by counting reduced
 * forms.
 *
 * For a fundamental discriminant D < 0 every form of discriminant D is
 * primitive, and h(D) is the number of reduced forms (a,b,c), b^2 - 4ac = D:
 *
 *     |b| <= a <= c,   and  b >= 0  whenever |b| = a or a = c.
 *
 * We sweep over all (a, b, c) at once, block by block in d = -D, adding 1 to
 * cnt[d] for each reduced form (a, +-b, c).  The total work is about
 * (pi/18) X^{3/2} increments (1.7e11 for X = 1e8), and needs no hypothesis.
 *
 * Output (stdout), one line per fundamental D, in increasing |D|:   D h
 *
 * Build:  gcc -O3 -march=native -fopenmp -o countforms countforms.c -lm
 * Usage:  countforms X [blocksize]
 */
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef unsigned long ulong;

static ulong *primes;
static long nprimes;

static void
init_primes(ulong lim)
{
  unsigned char *comp = calloc(lim + 1, 1);
  ulong p, q;
  primes = malloc(sizeof(ulong) * (lim / 2 + 10));
  for (p = 2; p <= lim; p++)
  {
    if (comp[p]) continue;
    primes[nprimes++] = p;
    for (q = p * p; q <= lim; q += p) comp[q] = 1;
  }
  free(comp);
}

/* f[n - lo] = 1 iff n squarefree, lo <= n < hi */
static void
squarefree(unsigned char *f, ulong lo, ulong hi)
{
  long i;
  memset(f, 1, hi - lo);
  for (i = 0; i < nprimes; i++)
  {
    ulong p2 = primes[i] * primes[i], m;
    if (p2 >= hi) break;
    for (m = ((lo + p2 - 1) / p2) * p2; m < hi; m += p2) f[m - lo] = 0;
  }
}

/* is -d a fundamental discriminant? sf1 covers [lo,hi), sf4 covers [lo/4, ..) */
static int
negfund(ulong d, ulong lo, const unsigned char *sf1, const unsigned char *sf4)
{
  if (d & 1) return (d & 3) == 3 && sf1[d - lo];
  if (d & 3) return 0;
  {
    ulong n = d >> 2, r = n & 3;
    return (r == 1 || r == 2) && sf4[n - lo / 4];
  }
}

/* count reduced forms of discriminant -d for d0 <= d < d1 into cnt[d - d0] */
static void
count_block(uint32_t *cnt, ulong d0, ulong d1)
{
  ulong a, b, amax = (ulong)sqrt((double)(d1 - 1) / 3.0) + 1;
  memset(cnt, 0, (d1 - d0) * sizeof(uint32_t));
  for (a = 1; a <= amax; a++)
  {
    ulong a4 = 4 * a;
    if (3 * a * a > d1 - 1) break;          /* d = 4ac - b^2 >= 3a^2 */
    for (b = 0; b <= a; b++)
    {
      ulong b2 = b * b, cmin, cmax, d;
      uint32_t w = (b > 0 && b < a) ? 2 : 1; /* (a, +-b, c) */
      cmin = (d0 + b2 + a4 - 1) / a4;       /* 4ac - b^2 >= d0 */
      if (cmin < a) cmin = a;
      cmax = (d1 - 1 + b2) / a4;            /* 4ac - b^2 <= d1 - 1 */
      if (cmin > cmax) continue;
      d = a4 * cmin - b2 - d0;
      if (cmin == a && w == 2)
      { /* a = c: only b >= 0 is reduced */
        cnt[d] += 1; d += a4; cmin++;
      }
      for (; cmin <= cmax; cmin++, d += a4) cnt[d] += w;
    }
  }
}

int
main(int argc, char **argv)
{
  ulong X, B, nblocks, k;
  if (argc < 2) { fprintf(stderr, "usage: %s X [blocksize]\n", argv[0]); return 1; }
  X = strtoul(argv[1], NULL, 10);
  B = argc > 2 ? strtoul(argv[2], NULL, 10) : (1UL << 18);
  init_primes((ulong)sqrt((double)X) + 2);
  nblocks = (X + B - 1) / B;

#pragma omp parallel
  {
    uint32_t *cnt = malloc(B * sizeof(uint32_t));
    unsigned char *sf1 = malloc(B), *sf4 = malloc(B / 4 + 2);
    char *out = malloc(B * 24);
#pragma omp for ordered schedule(dynamic, 1)
    for (k = 0; k < nblocks; k++)
    {
      ulong d0 = k * B, d1 = d0 + B, d;
      size_t len = 0;
      if (d1 > X) d1 = X;
      if (d0 < 3) d0 = 3;
      count_block(cnt, d0, d1);
      squarefree(sf1, d0, d1);
      squarefree(sf4, d0 / 4, (d1 - 1) / 4 + 1);
      for (d = d0; d < d1; d++)
        if (negfund(d, d0, sf1, sf4))
          len += sprintf(out + len, "-%lu %u\n", d, cnt[d - d0]);
#pragma omp ordered
      fwrite(out, 1, len, stdout);
    }
    free(cnt); free(sf1); free(sf4); free(out);
  }
  return 0;
}
