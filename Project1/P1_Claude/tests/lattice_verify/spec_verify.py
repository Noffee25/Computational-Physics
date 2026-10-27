"""Independent spectral test (verifier).

LLL: sympy DomainMatrix.lll (third-party implementation, not the lattice agent's code).
Shortest vectors: box enumeration over LLL coordinates, |c_i| <= R * ||column i of B^{-1}||,
which provably contains every lattice vector with |v| <= R. Exact integer recheck.

Coset gap: d_t = max_k h(k) / nu_t(a, 2^k), with h(k) = 6*2^k/2^31 (k<=28), 1/2 (k=29), 1 (k>=30).
(For s with 2-adic valuation k of s.v, the occupied residues of s.X mod 2^31 have largest
cyclic gap h(k); h is nondecreasing in k, so the max over s reduces to the max over k of
h(k)/shortest vector of the modulus-2^k dual lattice.)
"""
import itertools, json, math, sys
import numpy as np
from sympy import ZZ
from sympy.polys.matrices import DomainMatrix

GAM_POW = {2: 4/3, 3: 2., 4: 4., 5: 8., 6: 64/3, 7: 64., 8: 256.}


def dual_basis(a, M, t):
    B = [[M] + [0] * (t - 1)]
    p = 1
    for j in range(1, t):
        p = (p * a) % M
        r = [0] * t
        r[0] = (-p) % M
        r[j] = 1
        B.append(r)
    return B


def lll(B):
    t = len(B)
    D = DomainMatrix([[ZZ(int(x)) for x in row] for row in B], (t, t), ZZ)
    R = D.lll()
    return [[int(x) for x in row] for row in R.to_Matrix().tolist()]


def shortest(a, M, t):
    B = lll(dual_basis(a, M, t))
    Bf = np.array(B, dtype=float)
    R2 = min(sum(x * x for x in b) for b in B)
    R = math.sqrt(R2) * (1 + 1e-9)
    Binv = np.linalg.inv(Bf)
    bounds = [int(math.floor(R * np.linalg.norm(Binv[:, i]) + 1e-6)) for i in range(t)]
    ranges = [np.arange(-b, b + 1) for b in bounds]
    ncomb = np.prod([len(r) for r in ranges])
    best = None
    vecs = []
    # enumerate in chunks over the first coordinate
    grids = np.meshgrid(*ranges[1:], indexing='ij') if t > 1 else []
    rest = np.stack([g.ravel() for g in grids], axis=1) if t > 1 else np.zeros((1, 0))
    Bi = np.array(B, dtype=object)
    for c0 in ranges[0]:
        C = np.concatenate([np.full((rest.shape[0], 1), c0), rest], axis=1).astype(np.int64)
        V = C @ np.array(B, dtype=np.int64)  # components small enough? check below
        n2 = (V.astype(float) ** 2).sum(axis=1)
        sel = (n2 <= R2 * (1 + 1e-9)) & (n2 > 0)
        for row in C[sel]:
            v = [sum(int(row[i]) * B[i][j] for i in range(t)) for j in range(t)]
            nn = sum(x * x for x in v)
            if nn == 0:
                continue
            vecs.append((nn, v))
    n2min = min(nn for nn, v in vecs)
    sv = set()
    for nn, v in vecs:
        if nn == n2min:
            for x in v:
                if x != 0:
                    sgn = 1 if x > 0 else -1
                    break
            sv.add(tuple(sgn * x for x in v))
    # exact membership check
    for v in sv:
        tot, p = 0, 1
        for c in v:
            tot += c * p
            p = p * a % M
        assert tot % M == 0, (a, M, v)
    return n2min, sorted(sv), int(ncomb)


def merit(n2, M, t):
    nu = math.sqrt(n2)
    S = nu / (math.sqrt(GAM_POW[t] ** (1.0 / t)) * M ** (1.0 / t))
    mu = math.pi ** (t / 2) * nu ** t / (math.gamma(t / 2 + 1) * M)
    return S, mu


def hk(k):
    if k >= 30:
        return 1.0
    if k == 29:
        return 0.5
    return 6.0 * 2 ** k / 2 ** 31


if __name__ == '__main__':
    out = {}
    gens = [('RANDU', 65539, 2 ** 30), ('RANDM', 6539, 2 ** 30), ('MINSTD', 16807, 2 ** 31 - 1),
            ('MINSTD2', 48271, 2 ** 31 - 1), ('LCG69069', 69069, 2 ** 32)]
    for name, a, M in gens:
        for t in range(2, 9):
            n2, sv, nc = shortest(a, M, t)
            S, mu = merit(n2, M, t)
            print(f"{name:9s} t={t} nu2={n2:12d} S={S:.4f} mu={mu:.4g} 1/nu={1/math.sqrt(n2):.4e} sv={sv[:2]} box={nc}")
            out[f"{name}_{t}"] = dict(n2=n2, S=S, mu=mu, sv=sv)
    # naive / reduced moduli for RANDU / RANDM
    for name, a in [('RANDU', 65539), ('RANDM', 6539)]:
        for t in range(2, 9):
            r = {}
            for e in (31, 30, 29):
                n2, sv, _ = shortest(a, 2 ** e, t)
                r[e] = n2
            # true coset gap
            best = (0, None, None)
            for k in range(1, 31):
                n2, sv, _ = shortest(a, 2 ** k, t)
                g = hk(k) / math.sqrt(n2)
                if g > best[0]:
                    best = (g, k, sv[0])
            nu30 = math.sqrt(r[30])
            St = (1 / best[0]) / (math.sqrt(GAM_POW[t] ** (1.0 / t)) * (2 ** 30) ** (1.0 / t))
            print(f"{name} t={t} nu2(2^31)={r[31]} nu2(2^30)={r[30]} nu2(2^29)={r[29]} "
                  f"1/nu30={1/nu30:.4e} 1/nu29={1/math.sqrt(r[29]):.4e} true d={best[0]:.4e} k={best[1]} s={best[2]} "
                  f"ratio={best[0]*nu30:.3f} S_true={St:.4f}")
            out[f"{name}_{t}_coset"] = dict(n31=r[31], n30=r[30], n29=r[29], d=best[0], k=best[1], s=best[2])
    json.dump(out, open('/Users/noah/Downloads/ComputationalPhysicsPHYS25000/Project1/P1_Claude/tests/lattice_verify/spec_verify.json', 'w'), indent=1)
