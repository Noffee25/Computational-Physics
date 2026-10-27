"""Exact finite-lattice internal energy and specific heat of the 2D Ising
model on an n x m torus (Kaufman 1949; Ferdinand and Fisher 1969).

Conventions (J = 1, k_B = 1, beta J = K):
  H = - sum_<ij> s_i s_j ,  Z = sum exp(K sum_<ij> s_i s_j)
  e = <sum_<ij> s_i s_j> / N  = (1/N) d lnZ / dK      (positive, about 1.45 at K_c)
  C = K^2 (<E^2> - <E>^2) / N = (K^2/N) d^2 lnZ / dK^2

Validation:
  (1) brute-force enumeration of all 2^(nm) states for small tori (4x4, 3x5, 4x5)
  (2) float64 row-to-row transfer matrix for 6x6, 8x8, 10x10
"""
import itertools, sys
import numpy as np
import mpmath as mp

mp.mp.dps = 60


def lnZ_kaufman(K, n, m):
    """ln Z for an n (rows) x m (columns) periodic lattice, Kaufman form.
    gamma_l: cosh gamma_l = cosh2K coth2K - cos(pi l / n), gamma_0 = 2K + ln tanh K (signed)."""
    K = mp.mpf(K)
    c = mp.cosh(2 * K) * mp.coth(2 * K)

    def gam(l):
        if l % (2 * n) == 0:
            return 2 * K + mp.log(mp.tanh(K))
        return mp.acosh(c - mp.cos(mp.pi * l / n))

    # work with logs of the four products, factoring out the largest
    Z1 = mp.mpf(1); Z2 = mp.mpf(1); Z3 = mp.mpf(1); Z4 = mp.mpf(1)
    for r in range(n):
        g_odd = gam(2 * r + 1)
        g_even = gam(2 * r)
        Z1 *= 2 * mp.cosh(m * g_odd / 2)
        Z2 *= 2 * mp.sinh(m * g_odd / 2)
        Z3 *= 2 * mp.cosh(m * g_even / 2)
        Z4 *= 2 * mp.sinh(m * g_even / 2)
    N = n * m
    return mp.log(mp.mpf(1) / 2) + (mp.mpf(N) / 2) * mp.log(2 * mp.sinh(2 * K)) + mp.log(Z1 + Z2 + Z3 + Z4)


def e_C_kaufman(K, n, m):
    N = n * m
    f = lambda k: lnZ_kaufman(k, n, m)
    d1 = mp.diff(f, mp.mpf(K), 1)
    d2 = mp.diff(f, mp.mpf(K), 2)
    return d1 / N, mp.mpf(K) ** 2 * d2 / N


def dos_bruteforce(n, m):
    """Exact density of states g(B) with B = sum_<ij> s_i s_j, by enumeration."""
    N = n * m
    idx = np.arange(2 ** N, dtype=np.int64)
    spins = 1 - 2 * ((idx[:, None] >> np.arange(N)) & 1)       # (2^N, N)
    S = spins.reshape(-1, n, m)
    B = (S * np.roll(S, 1, axis=1)).sum(axis=(1, 2)) + (S * np.roll(S, 1, axis=2)).sum(axis=(1, 2))
    vals, counts = np.unique(B, return_counts=True)
    return vals, counts


def e_C_from_dos(K, vals, counts, N):
    K = mp.mpf(K)
    w = [mp.mpf(int(c)) * mp.exp(K * int(b)) for b, c in zip(vals, counts)]
    Z = mp.fsum(w)
    m1 = mp.fsum(wi * int(b) for wi, b in zip(w, vals)) / Z
    m2 = mp.fsum(wi * int(b) ** 2 for wi, b in zip(w, vals)) / Z
    return m1 / N, K ** 2 * (m2 - m1 ** 2) / N


def lnZ_transfer(K, L):
    """float64 transfer matrix, L x L torus (rows of L spins, L rows)."""
    nconf = 2 ** L
    idx = np.arange(nconf)
    s = 1 - 2 * ((idx[:, None] >> np.arange(L)) & 1)          # (nconf, L)
    intra = (s * np.roll(s, 1, axis=1)).sum(axis=1)             # horizontal bonds in a row
    inter = s @ s.T                                             # vertical bonds between rows
    # symmetric T = exp(K intra/2) exp(K inter) exp(K intra/2)
    shift = 2.0 * L * K
    T = np.exp(K * (inter + 0.5 * intra[:, None] + 0.5 * intra[None, :]) - shift / L)
    lam = np.linalg.eigvalsh(T)
    lmax = np.max(np.abs(lam))
    return L * np.log(lmax) + np.log(np.sum((lam / lmax) ** L)) + shift


def e_C_transfer(K, L, h=2e-3):
    N = L * L
    f = lambda k: lnZ_transfer(k, L)
    # 5-point stencils, Richardson with h and h/2
    def d1(hh):
        return (-f(K + 2 * hh) + 8 * f(K + hh) - 8 * f(K - hh) + f(K - 2 * hh)) / (12 * hh)
    def d2(hh):
        return (-f(K + 2 * hh) + 16 * f(K + hh) - 30 * f(K) + 16 * f(K - hh) - f(K - 2 * hh)) / (12 * hh * hh)
    D1 = (16 * d1(h / 2) - d1(h)) / 15
    D2 = (16 * d2(h / 2) - d2(h)) / 15
    return D1 / N, K * K * D2 / N


if __name__ == "__main__":
    Kc = mp.log(1 + mp.sqrt(2)) / 2
    out = []
    print("K_c =", mp.nstr(Kc, 20), "  p_add = 1 - exp(-2K_c) =", mp.nstr(1 - mp.exp(-2 * Kc), 20))
    print("\n(1) brute force vs Kaufman")
    maxdev = 0.0
    for (n, m) in [(4, 4), (3, 5), (5, 3), (4, 5), (2, 6)]:
        vals, counts = dos_bruteforce(n, m)
        for K in [mp.mpf("0.2"), Kc, mp.mpf("0.7")]:
            eb, Cb = e_C_from_dos(K, vals, counts, n * m)
            ek, Ck = e_C_kaufman(K, n, m)
            dev = max(abs(eb - ek), abs(Cb - Ck))
            maxdev = max(maxdev, float(dev))
            print(f"  {n}x{m} K={mp.nstr(K,8):>10}: e_brute={mp.nstr(eb,15)} e_exact={mp.nstr(ek,15)}  C_brute={mp.nstr(Cb,15)} C_exact={mp.nstr(Ck,15)}  |dev|={mp.nstr(dev,3)}")
    print("  max |deviation| brute force vs Kaufman:", maxdev)
    out.append(("bruteforce_maxdev", maxdev))
    print("\n(2) float64 transfer matrix vs Kaufman at K_c")
    for L in [6, 8, 10]:
        et, Ct = e_C_transfer(float(Kc), L)
        ek, Ck = e_C_kaufman(Kc, L, L)
        print(f"  {L}x{L}: e_TM={et:.12f} e_exact={mp.nstr(ek,13)}  C_TM={Ct:.10f} C_exact={mp.nstr(Ck,11)}  de={float(et-ek):.2e} dC={float(Ct-Ck):.2e}")
        out.append((f"TM{L}_de", float(et - ek))); out.append((f"TM{L}_dC", float(Ct - Ck)))
    print("\n(3) Exact values at K_c (L x L torus)")
    with open("/Users/noah/Downloads/ComputationalPhysicsPHYS25000/Project1/P1_Claude/tests/ising/exact_Kc.txt", "w") as fh:
        fh.write("# L  e_exact  C_exact   (K = K_c, e = <sum s_i s_j>/N, C = K^2 var(E)/N)\n")
        for L in [4, 8, 16, 24, 32, 48, 64]:
            ek, Ck = e_C_kaufman(Kc, L, L)
            print(f"  L={L:3d}: e = {mp.nstr(ek, 16)}   C = {mp.nstr(Ck, 16)}")
            fh.write(f"{L} {mp.nstr(ek, 20)} {mp.nstr(Ck, 20)}\n")
    print("  infinite lattice e(K_c) = sqrt(2) =", mp.nstr(mp.sqrt(2), 16))
