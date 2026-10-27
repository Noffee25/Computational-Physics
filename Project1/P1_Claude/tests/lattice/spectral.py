"""Exact spectral test for multiplicative / mixed linear congruential generators.

All arithmetic is exact (Python integers and fractions.Fraction):
  * dual lattice  L*_t(a, M) = { s in Z^t : s_1 + a s_2 + ... + a^{t-1} s_t = 0 (mod M) }
  * LLL reduction (delta = 99/100) with exact rational Gram-Schmidt
  * Fincke-Pohst enumeration of every lattice vector with |s|^2 <= R^2, where the
    interval bounds are found with float sqrt and then re-checked exactly, so the
    enumeration can neither miss nor invent a vector.

nu_t(a, M) = min |s| over nonzero s in L*_t.  Distance between adjacent covering
hyperplanes (unit cube) = 1/nu_t when the point set is a full coset of the lattice.
"""
from fractions import Fraction
from math import isqrt, sqrt, floor, ceil, pi, gamma
import itertools

HERMITE_POW = {1: Fraction(1), 2: Fraction(4, 3), 3: Fraction(2), 4: Fraction(4), 5: Fraction(8),
               6: Fraction(64, 3), 7: Fraction(64), 8: Fraction(256)}  # gamma_t^t (exact)


def hermite(t):
    return float(HERMITE_POW[t]) ** (1.0 / t)


def dual_basis(a, M, t):
    B = []
    row = [0] * t
    row[0] = M
    B.append(row)
    p = 1
    for j in range(1, t):
        p = (p * a) % M
        row = [0] * t
        row[0] = (-p) % M
        if row[0] > M // 2:
            row[0] -= M
        row[j] = 1
        B.append(row)
    return B


def dot(u, v):
    return sum(x * y for x, y in zip(u, v))


def gram_schmidt(B):
    n = len(B)
    Bs = []
    mu = [[Fraction(0)] * n for _ in range(n)]
    Bn = []
    for i in range(n):
        v = [Fraction(x) for x in B[i]]
        for j in range(i):
            mu[i][j] = Fraction(dot(B[i], Bs[j])) / Bn[j] if Bn[j] != 0 else Fraction(0)
            v = [vi - mu[i][j] * bj for vi, bj in zip(v, Bs[j])]
        Bs.append(v)
        Bn.append(sum(x * x for x in v))
    return Bs, mu, Bn


def lll(B, delta=Fraction(99, 100)):
    B = [list(r) for r in B]
    n = len(B)
    Bs, mu, Bn = gram_schmidt(B)
    k = 1
    while k < n:
        for j in range(k - 1, -1, -1):
            q = round(mu[k][j])
            if q != 0:
                B[k] = [x - q * y for x, y in zip(B[k], B[j])]
                Bs, mu, Bn = gram_schmidt(B)
        if Bn[k] >= (delta - mu[k][k - 1] ** 2) * Bn[k - 1]:
            k += 1
        else:
            B[k], B[k - 1] = B[k - 1], B[k]
            Bs, mu, Bn = gram_schmidt(B)
            k = max(k - 1, 1)
    return B


def enumerate_ball(B, R2, maxcount=10 ** 6):
    """All nonzero lattice vectors v (up to sign) with |v|^2 <= R2 (exact).
    B must be a basis (preferably LLL-reduced). Returns list of integer vectors."""
    n = len(B)
    Bs, mu, Bn = gram_schmidt(B)
    R2 = Fraction(R2)
    out = []
    x = [0] * n

    def rec(i, rem):
        # center c_i = - sum_{j>i} x_j mu[j][i]
        c = -sum((x[j] * mu[j][i] for j in range(i + 1, n)), Fraction(0))
        # (x_i - c)^2 Bn[i] <= rem
        w = sqrt(float(rem / Bn[i])) if rem > 0 else 0.0
        lo = floor(float(c) - w) - 1
        hi = ceil(float(c) + w) + 1
        for xi in range(lo, hi + 1):
            d = (xi - c) ** 2 * Bn[i]
            if d <= rem:
                x[i] = xi
                if i == 0:
                    if any(x):
                        v = [sum(x[k] * B[k][col] for k in range(n)) for col in range(len(B[0]))]
                        out.append(v)
                        if len(out) > maxcount:
                            raise RuntimeError("enumeration exceeded maxcount")
                else:
                    rec(i - 1, rem - d)
        x[i] = 0

    rec(n - 1, R2)
    # keep one representative of each +/- pair
    seen = set()
    res = []
    for v in out:
        key = tuple(v)
        neg = tuple(-c for c in v)
        if neg in seen or key in seen:
            continue
        seen.add(key)
        res.append(v)
    return res


def canon(v):
    """Sign convention: first nonzero component positive."""
    for c in v:
        if c != 0:
            return list(v) if c > 0 else [-x for x in v]
    return list(v)


def shortest(a, M, t):
    """Exact nu_t^2 and all shortest dual vectors (up to sign)."""
    if t == 1:
        return M * M, [[M]]
    B = lll(dual_basis(a, M, t))
    R2 = min(dot(b, b) for b in B)
    vs = enumerate_ball(B, R2)
    n2 = min(dot(v, v) for v in vs)
    best = sorted({tuple(canon(v)) for v in vs if dot(v, v) == n2})
    return n2, [list(b) for b in best]


def in_dual(s, a, M):
    p, tot = 1, 0
    for c in s:
        tot += c * p
        p = (p * a) % M
    return tot % M == 0


def sdotv(s, a, M):
    p, tot = 1, 0
    for c in s:
        tot += c * p
        p = (p * a) % M
    return tot % M


def mu_knuth(nu2, t, M):
    return pi ** (t / 2) * float(nu2) ** (t / 2) / (gamma(t / 2 + 1) * M)


def merit(nu2, t, M):
    """Normalized figure of merit S_t = nu_t / (gamma_t^{1/2} M^{1/t}) in (0, 1]."""
    return sqrt(nu2) / (sqrt(hermite(t)) * M ** (1.0 / t))


# ---------------------------------------------------------------------------
# Generators with m = 2^e, a = 3 (mod 8), odd seed x0: reachable states are
# y = x0 or 3 x0 (mod 8) (2^{e-2} of them).  For an integer normal s, put
# g = s.v mod 2^e; the values s.X mod 2^e over all tuples form the two cosets
# x0 g + H, 3 x0 g + H with H = 2^{min(k+3,e)} Z, k = v2(g).  h(s) below is the
# largest cyclic gap of that set divided by 2^e; the largest empty slab between
# adjacent planes orthogonal to s has width h(s)/|s| in the unit cube.
# ---------------------------------------------------------------------------
def v2(g, e):
    if g == 0:
        return e
    k = 0
    while g % 2 == 0:
        g //= 2
        k += 1
    return k


def h_gap(s, a, e, x0=1):
    m = 1 << e
    g = sdotv(s, a, m)
    k = v2(g, e)
    d = 1 << min(k + 3, e)
    delta = (2 * x0 * g) % d
    gap = d if delta == 0 else max(delta, d - delta)
    return Fraction(gap, m), k


def coset_spectral(a, e, t, x0=1):
    """Exact max distance between adjacent covering hyperplanes for the 2^{e-2}
    tuples of a multiplicative LCG mod 2^e with a = 3 (mod 8), odd seed x0.
    Returns dict with nu_t for moduli 2^e, 2^{e-1}, 2^{e-2}, the uniform-family
    answer 1/nu(2^{e-1}) and the true max gap over all families."""
    assert a % 8 == 3 and x0 % 2 == 1
    res = {}
    for shift in (0, 1, 2, 3):
        n2, vs = shortest(a, 1 << (e - shift), t)
        res['nu2_m%d' % shift] = n2
        res['vec_m%d' % shift] = vs
    nu30 = sqrt(res['nu2_m1'])
    best = 1.0 / nu30
    best_vec = res['vec_m1'][0]
    best_h = Fraction(1)
    # vectors with k <= e-2 can only beat 'best' if |s| < h_k * nu30
    cand = []
    for k in range(0, e - 1):
        hk = min(Fraction(6 * (1 << k), 1 << e), Fraction(3, 4)) if k <= e - 3 else Fraction(1, 2)
        R = float(hk) * nu30
        if R < 1:
            continue
        B = lll(dual_basis(a, 1 << k, t)) if k > 0 else [[1 if i == j else 0 for j in range(t)] for i in range(t)]
        R2 = int(R * R) + 1
        vs = enumerate_ball(B, R2, maxcount=2 * 10 ** 5)
        for v in vs:
            h, kk = h_gap(v, a, e, x0)
            r = float(h) / sqrt(dot(v, v))
            cand.append((r, v, h, kk))
            if r > best + 1e-15:
                best, best_vec, best_h = r, canon(v), h
    res['dmax'] = best
    res['dmax_vec'] = best_vec
    res['dmax_h'] = best_h
    res['uniform_d'] = 1.0 / nu30
    return res


def brute_force_nu2(a, M, t):
    R = int(sqrt(hermite(t)) * M ** (1.0 / t)) + 1
    best = None
    for s in itertools.product(range(-R, R + 1), repeat=t):
        if not any(s):
            continue
        if in_dual(s, a, M):
            n2 = dot(s, s)
            if best is None or n2 < best:
                best = n2
    return best


if __name__ == "__main__":
    import random
    random.seed(1)
    # self-validation against brute force on small moduli
    nbad = 0
    ntest = 0
    for M in (1009, 1024, 2039, 4096):
        for t in (2, 3, 4, 5):
            for _ in range(3):
                a = random.randrange(2, M)
                if M % 2 == 0:
                    a |= 1
                n2, _v = shortest(a, M, t)
                bf = brute_force_nu2(a, M, t)
                ntest += 1
                if n2 != bf:
                    nbad += 1
                    print("MISMATCH", M, a, t, n2, bf)
    print("brute-force validation: %d cases, %d mismatches" % (ntest, nbad))
