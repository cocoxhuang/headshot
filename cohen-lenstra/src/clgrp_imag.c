/*
 * clgrp_imag.c -- class groups of all imaginary quadratic fields with
 * fundamental discriminant D, lo <= |D| < hi, WITHOUT assuming GRH.
 *
 * 1. h(D) is obtained by counting reduced forms (a,b,c), |b| <= a <= c,
 *    b >= 0 if |b| = a or a = c, block by block in |D| (as in countforms.c).
 * 2. For every prime p with p^2 | h, the p-Sylow subgroup S_p is built
 *    explicitly: the (h/p^e)-th powers of prime forms of small norm lie in S_p,
 *    and we add them (as explicit cosets) until the subgroup has order p^e =
 *    |S_p|, which is known exactly from step 1.  The isomorphism type of S_p is
 *    read off from the numbers of elements of order dividing p^j.  If p || h,
 *    S_p = Z/p.  For p = 2, genus theory gives the 2-rank r = omega(D) - 1, so
 *    only the subgroup of squares 2*S_2, of order 2^(v_2(h) - r), has to be
 *    built.  No hypothesis is used: h is exact and the stopping rules are exact.
 *
 * Output (stdout), one line per field, in increasing |D|, in the same format
 * as clgrp_pari.py:   D c1,c2,...,ck   (invariant factors, "1" if trivial).
 * With -s, a one-line summary of the elapsed work is printed to stderr.
 *
 * Build:  gcc -O3 -march=native -fopenmp -o clgrp_imag clgrp_imag.c -lm
 * Usage:  clgrp_imag lo hi [blocksize]
 */
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef int64_t i64;
typedef uint64_t u64;
typedef __int128 i128;

/* ------------------------------------------------------------------------ */
/* small primes and square roots modulo them                                 */

#define NQ 600                  /* prime forms of norm < ~4400 as generators */
static int nsmall;
static i64 small_primes[NQ];
static int *sqrt_tab[NQ];       /* sqrt_tab[i][x] = r with r^2 = x mod q, or -1 */

static u64 *primes;             /* primes up to sqrt(hi) for the sieves */
static long nprimes;

static void
init_tables(u64 lim)
{
  unsigned char *comp = calloc(lim + 10000, 1);
  u64 p, q;
  long i;
  primes = malloc(sizeof(u64) * (lim / 2 + 10000));
  for (p = 2; p < lim + 10000; p++)
  {
    if (comp[p]) continue;
    if (p <= lim) primes[nprimes++] = p;
    if (nsmall < NQ) small_primes[nsmall++] = (i64)p;
    for (q = p * p; q < lim + 10000; q += p) comp[q] = 1;
  }
  free(comp);
  for (i = 0; i < nsmall; i++)
  {
    i64 m = small_primes[i], r;
    sqrt_tab[i] = malloc(sizeof(int) * m);
    for (r = 0; r < m; r++) sqrt_tab[i][r] = -1;
    for (r = 0; r < m; r++) sqrt_tab[i][(r * r) % m] = (int)r;
  }
}

/* ------------------------------------------------------------------------ */
/* positive definite binary quadratic forms of discriminant D < 0            */

typedef struct { i64 a, b, c; } form;

static inline i64
floordiv(i64 x, i64 y)          /* y > 0 */
{
  i64 q = x / y;
  if ((x % y) && (x < 0)) q--;
  return q;
}

static inline i64
xgcd(i64 a, i64 b, i64 *u, i64 *v)   /* u a + v b = g = gcd(a,b) >= 0 */
{
  i64 u0 = 1, v0 = 0, u1 = 0, v1 = 1;
  while (b)
  {
    i64 q = a / b, t;
    t = a - q * b; a = b; b = t;
    t = u0 - q * u1; u0 = u1; u1 = t;
    t = v0 - q * v1; v0 = v1; v1 = t;
  }
  if (a < 0) { a = -a; u0 = -u0; v0 = -v0; }
  *u = u0; *v = v0;
  return a;
}

/* Cohen, Algorithm 5.4.2 */
static inline void
reduce(form *f)
{
  i64 a = f->a, b = f->b, c = f->c;
  for (;;)
  {
    if (b <= -a || b > a)
    {
      i64 q = floordiv(b, 2 * a), r = b - 2 * a * q;
      if (r > a) { r -= 2 * a; q++; }
      c -= (b + r) / 2 * q;
      b = r;
    }
    if (a > c) { i64 t = a; a = c; c = t; b = -b; continue; }
    if (a == c && b < 0) b = -b;
    break;
  }
  f->a = a; f->b = b; f->c = c;
}

/* Cohen, Algorithm 5.4.7: composition of primitive forms, then reduction */
static void
compose(form *r, const form *f1, const form *f2, i64 D)
{
  const form *g1 = f1, *g2 = f2;
  i64 s, n, y1, d, y2, x2, d1, v1, v2, rr, u, v;
  if (g1->a > g2->a) { g1 = f2; g2 = f1; }
  s = (g1->b + g2->b) / 2;
  n = g2->b - s;
  if (g2->a % g1->a == 0) { y1 = 0; d = g1->a; }
  else { d = xgcd(g2->a, g1->a, &u, &v); y1 = u; }
  if (s % d == 0) { y2 = -1; x2 = 0; d1 = d; }
  else { d1 = xgcd(s, d, &x2, &y2); y2 = -y2; }
  v1 = g1->a / d1;
  v2 = g2->a / d1;
  {
    i128 t = ((i128)y1 * y2 % v1) * n - (i128)x2 * g2->c;
    t %= v1;
    if (t < 0) t += v1;
    rr = (i64)t;
  }
  {
    /* normalise b mod 2a before computing c, to keep everything in range */
    form t;
    i64 a3 = v1 * v2, b3 = g2->b + 2 * v2 * rr, q = floordiv(b3, 2 * a3);
    b3 -= 2 * a3 * q;
    if (b3 > a3) b3 -= 2 * a3;
    t.a = a3; t.b = b3;
    t.c = (i64)(((i128)b3 * b3 - D) / (4 * (i128)a3));
#ifdef CHECK
    if (((i128)b3 * b3 - D) % (4 * (i128)a3)) { fprintf(stderr, "bad composition\n"); exit(3); }
#endif
    reduce(&t);
    *r = t;
  }
}

static void
power(form *r, const form *f, u64 m, i64 D)
{
  form x = *f, acc;
  int started = 0;
  while (m)
  {
    if (m & 1)
    {
      if (!started) { acc = x; started = 1; }
      else compose(&acc, &acc, &x, D);
    }
    m >>= 1;
    if (m) compose(&x, &x, &x, D);
  }
  if (!started) { acc.a = 1; acc.b = D & 1; acc.c = ((D & 1) - D) / 4; }
  *r = acc;
}

/* prime form of norm q = small_primes[i], or 0 if q is inert */
static int
prime_form(form *f, int i, i64 D)
{
  i64 q = small_primes[i], b;
  if (q == 2)
  {
    i64 m8 = ((D % 8) + 8) % 8;
    if (m8 == 1) b = 1;
    else if (m8 == 0) b = 0;      /* D = 8 mod 16 */
    else if (m8 == 4) b = 2;      /* D = 12 mod 16 */
    else return 0;                /* D = 5 mod 8: inert */
  }
  else
  {
    i64 x = ((D % q) + q) % q;
    int r = sqrt_tab[i][x];
    if (r < 0) return 0;
    b = r;
    if ((b & 1) != (D & 1)) b = q - b;   /* b = D mod 2 */
  }
  f->a = q; f->b = b; f->c = (b * b - D) / (4 * q);
  reduce(f);
  return 1;
}

/* ------------------------------------------------------------------------ */
/* explicit subgroups: hash set of reduced forms (keyed by (a,b))             */

typedef struct {
  u64 *keys;      /* 0 = empty */
  int *idx;
  form *elt;
  long n, cap, mask;
} fset;

static inline u64
fkey(const form *f) { return ((u64)f->a << 32) ^ (u64)(uint32_t)f->b ^ (1ULL << 63); }

static void
fset_init(fset *S, long maxn)
{
  long cap = 64;
  while (cap < 2 * maxn + 2) cap <<= 1;
  S->keys = calloc(cap, sizeof(u64));
  S->idx = malloc(cap * sizeof(int));
  S->elt = malloc((maxn + 1) * sizeof(form));
  S->n = 0; S->cap = cap; S->mask = cap - 1;
}

static void
fset_free(fset *S) { free(S->keys); free(S->idx); free(S->elt); }

static inline long
fset_find(const fset *S, const form *f)
{
  u64 k = fkey(f), h = (k * 0x9E3779B97F4A7C15ULL) >> 20;
  for (;; h++)
  {
    u64 kk = S->keys[h & S->mask];
    if (!kk) return -1;
    if (kk == k) return S->idx[h & S->mask];
  }
}

static inline void
fset_add(fset *S, const form *f)
{
  u64 k = fkey(f), h = (k * 0x9E3779B97F4A7C15ULL) >> 20;
  for (;; h++)
    if (!S->keys[h & S->mask])
    {
      S->keys[h & S->mask] = k;
      S->idx[h & S->mask] = (int)S->n;
      S->elt[S->n++] = *f;
      return;
    }
}

static inline int
is_one(const form *f) { return f->a == 1; }

/* The subgroup T = {x^m} has order p^2 (m = h/p^2 gives the p-Sylow S_p).
 * T = Z/p^2 iff some element has order p^2; otherwise two elements of order
 * p with g' not in <g> prove T = (Z/p)^2.  O(p) work. */
static int
sylow_p2(int *lambda, i64 D, u64 m, u64 p)
{
  fset C;                         /* the cyclic subgroup <g>, once found */
  int i, have = 0, res = -1;
  fset_init(&C, (long)p);
  for (i = 0; i < nsmall; i++)
  {
    form f, g, gp;
    if (!prime_form(&f, i, D)) continue;
    power(&g, &f, m, D);
    if (is_one(&g)) continue;
    power(&gp, &g, p, D);
    if (!is_one(&gp)) { lambda[0] = 2; res = 1; break; }   /* order p^2 */
    if (!have)
    {
      form x = g;
      form id = { 1, D & 1, ((D & 1) - D) / 4 };
      u64 t;
      fset_add(&C, &id);
      for (t = 1; t < p; t++) { fset_add(&C, &x); compose(&x, &x, &g, D); }
      have = 1;
    }
    else if (fset_find(&C, &g) < 0) { lambda[0] = 1; lambda[1] = 1; res = 2; break; }
  }
  fset_free(&C);
  return res;
}

/* Partition (exponents, decreasing) of the p-group T = {x^m : x in Cl},
 * known to have order p^e (m = h/p^e gives the p-Sylow subgroup S_p).
 * Returns its length, or -1 on failure. */
static int
sylow(int *lambda, i64 D, u64 m, u64 p, int e)
{
  u64 order = 1;
  int i, j, len = 0;
  fset S;
  for (i = 0; i < e; i++) order *= p;
  fset_init(&S, (long)order);
  {
    form id = { 1, D & 1, ((D & 1) - D) / 4 };
    fset_add(&S, &id);
  }
  for (i = 0; i < nsmall && (u64)S.n < order; i++)
  {
    form f, g, x;
    long base, t;
    if (!prime_form(&f, i, D)) continue;
    power(&g, &f, m, D);
    if (fset_find(&S, &g) >= 0) continue;
    base = S.n;
    x = g;
    while (fset_find(&S, &x) < 0)          /* add the coset x * H */
    {
      if ((u64)(S.n + base) > order) { fset_free(&S); return -1; }
      for (t = 0; t < base; t++)
      {
        form y;
        compose(&y, &x, &S.elt[t], D);
        fset_add(&S, &y);
      }
      compose(&x, &x, &g, D);
    }
  }
  if ((u64)S.n != order) { fset_free(&S); return -1; }
  {
    /* ordexp[t] = log_p of the order of element t, via x -> x^p */
    long n = S.n, t;
    int *par = malloc(n * sizeof(int)), *oe = malloc(n * sizeof(int));
    long cnt[64] = { 0 };
    for (t = 0; t < n; t++)
    {
      form y;
      power(&y, &S.elt[t], p, D);
      par[t] = (int)fset_find(&S, &y);
      oe[t] = -1;
      if (par[t] < 0) { free(par); free(oe); fset_free(&S); return -1; }
    }
    oe[0] = 0;                              /* element 0 is the identity */
    for (t = 0; t < n; t++)
    {
      long u = t, depth = 0;
      while (oe[u] < 0) { u = par[u]; depth++; }
      {
        int base_e = oe[u];
        long w = t;
        while (oe[w] < 0) { oe[w] = base_e + (int)depth; depth--; w = par[w]; }
      }
    }
    for (t = 0; t < n; t++) cnt[oe[t]]++;
    /* N_j = #{x : x^{p^j} = 1};  r_j = log_p(N_j / N_{j-1}) = #{i : lambda_i >= j} */
    {
      long N = cnt[0];
      int r[64], k;
      for (j = 1; j <= e; j++)
      {
        long Nj = N + cnt[j];
        u64 ratio = (u64)(Nj / N);
        k = 0;
        while (ratio > 1) { ratio /= p; k++; }
        r[j] = k;
        N = Nj;
      }
      for (i = 0; ; i++)
      {
        int lam = 0;
        for (j = 1; j <= e; j++) if (r[j] > i) lam++;
        if (!lam) break;
        lambda[len++] = lam;
      }
    }
    free(par); free(oe);
  }
  fset_free(&S);
  return len;
}

/* invariant factors of Cl(D), given h(D) and the 2-rank r2 = omega(D) - 1
 * (genus theory); returns number of factors (0 if trivial), or -1 on failure */
static int
class_group(u64 *cyc, i64 D, u64 h, int r2)
{
  u64 rem = h, p;
  int k = 0, i;
  for (i = 0; i < 64; i++) cyc[i] = 1;
  for (p = 2; rem > 1; p++)
  {
    int e = 0, lambda[64], len;
    if (p * p > rem) p = rem;               /* remaining cofactor is prime */
    while (rem % p == 0) { rem /= p; e++; }
    if (p == 2 && e < r2) return -1;         /* contradicts genus theory */
    if (!e) continue;
    if (p == 2)
    { /* S_2 has rank r2, so 2*S_2 has order 2^(e - r2); find its type mu and
         lift: lambda_i = mu_i + 1 for i < len(mu), then 1s up to r2 parts */
      int mu[64], lm = 0, f = e - r2;
      u64 m2 = h >> (e - 1);                /* x -> x^(2 h/2^e) maps Cl onto 2*S_2 */
      if (f == 1) { mu[0] = 1; lm = 1; }
      else if (f == 2) lm = sylow_p2(mu, D, m2, 2);
      else if (f > 2) lm = sylow(mu, D, m2, 2, f);
      if (lm < 0 || lm > r2) return -1;
      len = r2;
      for (i = 0; i < len; i++) lambda[i] = 1 + (i < lm ? mu[i] : 0);
    }
    else if (e == 1) { lambda[0] = 1; len = 1; }
    else if (e == 2) { if ((len = sylow_p2(lambda, D, h / (p * p), p)) < 0) return -1; }
    else
    {
      u64 m = h;
      for (i = 0; i < e; i++) m /= p;
      if ((len = sylow(lambda, D, m, p, e)) < 0) return -1;
    }
    for (i = 0; i < len; i++)
    {
      int t;
      for (t = 0; t < lambda[i]; t++) cyc[i] *= p;
    }
    if (len > k) k = len;
  }
  return k;
}

/* ------------------------------------------------------------------------ */
/* counting reduced forms in a block of |D|                                  */

static void
count_block(uint32_t *cnt, u64 d0, u64 d1)
{
  u64 a;
  memset(cnt, 0, (d1 - d0) * sizeof(uint32_t));
  for (a = 1; 3 * a * a <= d1 - 1; a++)
  {
    u64 a4 = 4 * a, b;
    /* lo(b) = ceil((d0 + b^2) / 4a) and hi(b) = floor((d1 - 1 + b^2) / 4a),
       updated incrementally in b (b^2 grows by 2b+1 < 4a) */
    u64 qlo = (d0 + a4 - 1) / a4, rlo = qlo * a4 - d0;          /* rlo in [0,4a) */
    u64 qhi = (d1 - 1) / a4, rhi = (d1 - 1) - qhi * a4;        /* rhi in [0,4a) */
    for (b = 0; b <= a; b++)
    {
      u64 cmin, cmax, d;
      uint32_t w;
      if (b)
      {
        u64 delta = 2 * b - 1;            /* b^2 - (b-1)^2 */
        if (rlo >= delta) rlo -= delta; else { rlo += a4 - delta; qlo++; }
        rhi += delta; if (rhi >= a4) { rhi -= a4; qhi++; }
      }
      cmin = qlo < a ? a : qlo;
      cmax = qhi;
      if (cmin > cmax) continue;
      w = (b > 0 && b < a) ? 2 : 1;
      if (b & 1)
      { /* d = 3 mod 4: every c */
        d = a4 * cmin - b * b - d0;
        if (cmin == a && w == 2) { cnt[d] += 1; d += a4; cmin++; }
        for (; cmin <= cmax; cmin++, d += a4) cnt[d] += w;
      }
      else
      { /* b = 2b', d = 4n with n = ac - b'^2; fundamental needs n = 1, 2 mod 4.
           Counts for other d are left incomplete (they are never read). */
        u64 bp2 = (b >> 1) * (b >> 1), am = a & 3, c, step, r[2];
        int nr, j;
        if (am == 0) continue;               /* n = -b'^2 = 0, 3 mod 4 */
        if (am & 1)
        { /* c = a^-1 (n + b'^2) mod 4, and a^-1 = a mod 4 */
          r[0] = (am * (1 + bp2)) & 3; r[1] = (am * (2 + bp2)) & 3;
          nr = 2; step = 4;
        }
        else { r[0] = 1; nr = 1; step = 2; } /* a = 2 mod 4: n = 1, 2 iff c odd */
        if (cmin == a && w == 2)
        {
          for (j = 0; j < nr; j++)
            if ((a & (step - 1)) == r[j]) cnt[a4 * a - b * b - d0] += 1;
          cmin++;
        }
        for (j = 0; j < nr; j++)
          for (c = cmin + ((r[j] - cmin) & (step - 1)), d = a4 * c - b * b - d0;
               c <= cmax; c += step, d += a4 * step)
            cnt[d] += w;
      }
    }
  }
}

static void
squarefree(unsigned char *f, u64 lo, u64 hi)
{
  long i;
  memset(f, 1, hi - lo);
  for (i = 0; i < nprimes; i++)
  {
    u64 p2 = primes[i] * primes[i], m;
    if (p2 >= hi) break;
    for (m = ((lo + p2 - 1) / p2) * p2; m < hi; m += p2) f[m - lo] = 0;
  }
}

/* om[i] = number of distinct prime factors of lo + i, lo <= n < hi */
static void
omega_sieve(unsigned char *om, u64 *rem, u64 lo, u64 hi)
{
  u64 n, len = hi - lo;
  long i;
  for (n = 0; n < len; n++) { rem[n] = lo + n; om[n] = 0; }
  for (i = 0; i < nprimes; i++)
  {
    u64 p = primes[i], m;
    if (p * p >= hi) break;
    for (m = ((lo + p - 1) / p) * p; m < hi; m += p)
    {
      om[m - lo]++;
      do rem[m - lo] /= p; while (rem[m - lo] % p == 0);
    }
  }
  for (n = 0; n < len; n++) if (rem[n] > 1) om[n]++;
}

static int
negfund(u64 d, u64 lo, const unsigned char *sf1, const unsigned char *sf4)
{
  if (d & 1) return (d & 3) == 3 && sf1[d - lo];
  if (d & 3) return 0;
  {
    u64 n = d >> 2, r = n & 3;
    return (r == 1 || r == 2) && sf4[n - lo / 4];
  }
}

int
main(int argc, char **argv)
{
  u64 lo, hi, B, nblocks, k;
  long failures = 0;
  if (argc < 3) { fprintf(stderr, "usage: %s lo hi [blocksize]\n", argv[0]); return 1; }
  lo = strtoull(argv[1], NULL, 10);
  hi = strtoull(argv[2], NULL, 10);
  B = argc > 3 ? strtoull(argv[3], NULL, 10) : (1UL << 18);
  if (lo < 3) lo = 3;
  if (hi <= lo) return 0;
  init_tables((u64)sqrt((double)hi) + 2);
  nblocks = (hi - lo + B - 1) / B;

#pragma omp parallel reduction(+:failures)
  {
    uint32_t *cnt = malloc(B * sizeof(uint32_t));
    unsigned char *sf1 = malloc(B), *sf4 = malloc(B / 4 + 2), *om = malloc(B);
    u64 *rem = malloc(B * sizeof(u64));
    char *out = malloc(B * 48);
#pragma omp for ordered schedule(dynamic, 1)
    for (k = 0; k < nblocks; k++)
    {
      u64 d0 = lo + k * B, d1 = d0 + B, d;
      size_t len = 0;
      if (d1 > hi) d1 = hi;
      count_block(cnt, d0, d1);
      squarefree(sf1, d0, d1);
      squarefree(sf4, d0 / 4, (d1 - 1) / 4 + 1);
      omega_sieve(om, rem, d0, d1);
      for (d = d0; d < d1; d++)
      {
        u64 cyc[64];
        int nc, i;
        if (!negfund(d, d0, sf1, sf4)) continue;
        nc = class_group(cyc, -(i64)d, cnt[d - d0], om[d - d0] - 1);
        len += sprintf(out + len, "-%lu ", (unsigned long)d);
        if (nc < 0)
        {
          failures++;
          fprintf(stderr, "FAILED at D = -%lu (h = %u)\n", (unsigned long)d, cnt[d - d0]);
          len += sprintf(out + len, "FAILED\n");
          continue;
        }
        if (nc == 0) len += sprintf(out + len, "1");
        for (i = 0; i < nc; i++)
          len += sprintf(out + len, i ? ",%lu" : "%lu", (unsigned long)cyc[i]);
        out[len++] = '\n';
      }
#pragma omp ordered
      fwrite(out, 1, len, stdout);
    }
    free(cnt); free(sf1); free(sf4); free(om); free(rem); free(out);
  }
  if (failures) fprintf(stderr, "%ld failures\n", failures);
  return failures ? 2 : 0;
}
