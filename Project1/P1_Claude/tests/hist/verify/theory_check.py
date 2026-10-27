"""Verifier: independent theory for the hist family.

Populations c_b are computed analytically (no enumeration of the generator):
  RANDU/RANDM seed 777: orbit = {y in [0,2^31): y = 1 or 3 mod 8}, x = y/2^31
  seed 65536: orbit = {2^16 u : u in [0,2^15), u = 1 or 3 mod 8}
  MINSTD: all s in [1, M-1], x = s/M, bin = floor(s B / M) (exact up to fp ties)
  RAND: k in [0, 2^22), weight 512 (top 510), x = k/2^22
Then E[D] and the empty fraction from the finite-population model.
"""
import numpy as np
from scipy.special import gammaln
from scipy import stats


def count13(n):
    """#{0 <= y < n : y mod 8 in {1,3}} for integer array n."""
    n = np.asarray(n, dtype=object)
    q = n // 8
    r = n % 8
    return np.array([2 * qq + (1 if rr > 1 else 0) + (1 if rr > 3 else 0) for qq, rr in zip(q, r)], dtype=np.int64)


def count13_fast(n):
    n = np.asarray(n, dtype=np.int64)
    q = n // 8
    r = n % 8
    return 2 * q + (r > 1) + (r > 3)


def pop_randu(B, scale=2**31, unit=1, b0=0, b1=None):
    # bin b holds y with b*scale/B <= y*unit < (b+1)*scale/B ; y counted in units
    b1 = B if b1 is None else b1
    b = np.arange(b0, b1 + 1, dtype=np.int64)
    # smallest integer u with u*unit*B >= b*scale  -> ceil(b*scale/(B*unit))
    lo = -((-b * scale) // (B * unit))
    F = count13_fast(lo)
    return np.diff(F)


def pop_minstd(B, M=2**31 - 1, b0=0, b1=None):
    b1 = B if b1 is None else b1
    b = np.arange(b0, b1 + 1, dtype=np.int64)
    # s in [1, M-1] with floor(s B / M) = b  <=> s >= ceil(b M / B)
    lo = -((-b * M) // B)
    lo = np.clip(lo, 1, M)
    return np.diff(lo)


def pop_rand(B, b0=0, b1=None):
    T = 2**22
    b1 = B if b1 is None else b1
    b = np.arange(b0, b1 + 1, dtype=np.int64)
    lo = -((-b * T) // B)          # first k in bin b
    lo = np.clip(lo, 0, T)
    nk = np.diff(lo)
    c = 512 * nk
    # top value k = T-1 has weight 510
    last = np.searchsorted(lo, T - 1, side="right") - 1
    if 0 <= last < len(c) and lo[last] <= T - 1 < lo[last + 1]:
        c[last] -= 2
    return c


def coc(popfun, B, chunk=10**7, **kw):
    """count-of-counts of a population function evaluated in chunks of bins"""
    from collections import Counter
    tot = Counter()
    for b0 in range(0, B, chunk):
        b1 = min(B, b0 + chunk)
        c = popfun(B, b0=b0, b1=b1, **kw)
        u, n = np.unique(c, return_counts=True)
        for a, m in zip(u, n):
            tot[int(a)] += int(m)
    k = np.array(sorted(tot)); m = np.array([tot[i] for i in k], float)
    return k, m


def pred(c, P, N, B, m=None):
    c = np.asarray(c, float)
    if m is None:
        m = np.ones_like(c)
    q, r = divmod(int(N), int(P))
    lam = N / B
    v = r * (c / P) * (1 - c / P) * (P - r) / (P - 1.0)
    a = q * c + r * c / P - lam
    S = np.sum(m * (v + a**2))
    D = S / (lam * (B - 1))
    sd = np.sqrt(np.sum(m * (2 * v**2 + 4 * v * a**2))) / (lam * (B - 1))
    if q >= 1:
        pe = np.sum(m[c == 0]) / B
    else:
        cc = c.astype(np.int64)
        lp = gammaln(P - cc + 1) - gammaln(P - cc - r + 1) - gammaln(P + 1) + gammaln(P - r + 1)
        pe = np.sum(m * np.exp(lp)) / B
    return D, sd, pe


if __name__ == "__main__":
    PU = 2**29
    PM = 2**31 - 2
    print("RANDU orbit population check: total", pop_randu(10).sum(), "= 2^29", 2**29)
    for B, N in [(10**5, 10**7), (10**6, 10**8), (10**7, 10**9), (10**8, 10**9)]:
        c, m = coc(pop_randu, B)
        D, sd, pe = pred(c, PU, N, B, m)
        print(f"RANDU B={B:.0e} N={N:.0e}: c range {c.min()}..{c.max()}  Dpred={D:.6f} sd={sd:.6f} empty={pe*B:.3f}"
              f"  z_pred={(D-1)*np.sqrt((B-1)/2):.3f}")
    for B, N in [(10**5, 10**7), (10**6, 10**8), (10**7, 10**9), (10**8, 10**9)]:
        c, m = coc(pop_minstd, B)
        D, sd, pe = pred(c, PM, N, B, m)
        print(f"MINSTD B={B:.0e} N={N:.0e}: c {c.min()}..{c.max()} sum {np.sum(c*m):.0f} Dpred={D:.6f} sd={sd:.6f} "
              f"empty={pe*B:.2f} z_pred={(D-1)*np.sqrt((B-1)/2):.3f}")
    f = 1e9 / PM
    print("  naive 1e8 (1-f)^21.47 =", 1e8 * (1 - f) ** (PM / 1e8))
    for B, N in [(10**5, 10**7), (10**6, 10**8), (10**7, 10**9), (10**8, 10**9)]:
        c, m = coc(pop_rand, B)
        D, sd, pe = pred(c, PM, N, B, m)
        print(f"RAND B={B:.0e} N={N:.0e}: c {c.min()}..{c.max()} sum {np.sum(c*m):.0f} Dpred={D:.6f} sd={sd:.6f} "
              f"empty={pe*B:.1f} z_pred={(D-1)*np.sqrt((B-1)/2):.3f}")
    # RAND filled bins at B=1e7: binomial(512, N/P)
    p = 1e9 / PM
    print("RAND filled: Bin(512,%.6f) mean %.3f sd %.4f" % (p, 512 * p, np.sqrt(512 * p * (1 - p))))
    # seed 65536: orbit 2^16 * {u in [0,2^15): u = 1,3 mod 8}; x = 2^16 u / 2^31 = u / 2^15
    P6 = 8192
    for d in range(1, 9):
        B = 10**d
        N = 100 * B if d < 8 else 10**9
        c, m = coc(pop_randu, B, scale=2**15, unit=1)
        D, sd, pe = pred(c, P6, N, B, m)
        print(f"seed65536 B=1e{d} N={N:.0e}: Dpred={D:.6f} sd={sd:.6f} empty={pe*B:.0f} q,r={divmod(N,P6)}")
    # exact tail expectations for the 4-sigma criterion at B=1e7, N=1e9
    B, N = 10**7, 10**9
    lo = B * stats.binom.cdf(60, N, 1 / B)
    hi = B * stats.binom.sf(139, N, 1 / B)
    print("E4LO, E4HI at B=1e7:", lo, hi, "sum", lo + hi)
    print("Poisson empties at lambda=10, B=1e8:", 1e8 * np.exp(-10), np.sqrt(1e8 * np.exp(-10) * (1 - np.exp(-10))))
    # scaling law
    for P in [PU, PM]:
        Nst = (3 * np.sqrt(200) * P) ** (2 / 3)
        print("N* for P=%d: %.3g  (12.16 P^(2/3) = %.3g)" % (P, Nst, 12.16 * P ** (2 / 3)))
