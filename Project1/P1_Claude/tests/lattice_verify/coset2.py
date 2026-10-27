"""Corrected coset reduction: h(k) is monotone except h(e-3)=3/4 > h(e-2)=1/2.
d = max( max_{k != e-3} h(k)/nu(L_k),  0.75 / min{|s| : s in L_{e-3}, v2(s.v) == e-3} )."""
import math, itertools, sys
import numpy as np
sys.path.insert(0, '/Users/noah/Downloads/ComputationalPhysicsPHYS25000/Project1/P1_Claude/tests/lattice_verify')
from spec_verify import shortest, dual_basis, lll, GAM_POW


def h_e(k, e):
    if k >= e - 1: return 1.0
    if k == e - 2: return 0.5
    return 6.0 * 2 ** k / 2 ** e


def enum_upto(a, M, t, R2):
    B = lll(dual_basis(a, M, t))
    Bf = np.array(B, dtype=float)
    R = math.sqrt(R2) * (1 + 1e-9)
    Binv = np.linalg.inv(Bf)
    bounds = [int(math.floor(R * np.linalg.norm(Binv[:, i]) + 1e-6)) for i in range(t)]
    grids = np.meshgrid(*[np.arange(-b, b + 1) for b in bounds], indexing='ij')
    C = np.stack([g.ravel() for g in grids], axis=1).astype(np.int64)
    out = []
    for i0 in range(0, len(C), 2_000_000):
        Cc = C[i0:i0 + 2_000_000]
        V = Cc @ np.array(B, dtype=np.int64)
        n2 = (V.astype(float) ** 2).sum(1)
        for v in V[(n2 <= R2 * (1 + 1e-9)) & (n2 > 0)]:
            out.append([int(x) for x in v])
    return out


def v2(g, e):
    g %= 2 ** e
    return e if g == 0 else (g & -g).bit_length() - 1


def coset_gap(a, e, t):
    M = 2 ** e
    best, arg = 0.0, None
    for k in range(1, e):
        if k == e - 3:
            continue
        n2, sv, _ = shortest(a, 2 ** k, t)
        val = h_e(k, e) / math.sqrt(n2)
        if val > best:
            best, arg = val, (k, sv[0])
    R2 = (0.75 / best) ** 2
    cand = enum_upto(a, 2 ** (e - 3), t, R2)
    for s in cand:
        g = sum(c * pow(a, j, M) for j, c in enumerate(s))
        if v2(g, e) == e - 3:
            val = 0.75 / math.sqrt(sum(c * c for c in s))
            if val > best + 1e-15:
                best, arg = val, (e - 3, tuple(s))
    return best, arg


if __name__ == '__main__':
    if sys.argv[1:] == ['small']:
        from coset_bruteforce import brute
        import coset_bruteforce as cb
        cb.formula = lambda a, e, t: coset_gap(a, e, t)[0]
        mism = n = 0
        for e in (9, 10, 11, 12):
            for a in sorted({2 ** (e // 2) + 3, 6539 % 2 ** e, 65539 % 2 ** e, 2 ** e - 5, 11, 27, 19, 35 % 2**e, 43}):
                if a % 8 != 3: continue
                for t in (2, 3, 4):
                    if t == 4 and e > 10: continue
                    for x0 in (1, 7, 777 % 2 ** e):
                        P, d0, b = brute(a, e, t, x0)
                        ok = abs(d0 - b) < 1e-12
                        n += 1; mism += (not ok)
                        if not ok: print('MISMATCH', e, a, t, x0, d0, b)
        print('small-modulus cases', n, 'mismatches', mism)
    else:
        for name, a in [('RANDU', 65539), ('RANDM', 6539)]:
            for t in range(2, 9):
                d, arg = coset_gap(a, 31, t)
                n30, sv30, _ = shortest(a, 2 ** 30, t)
                St = (1 / d) / (math.sqrt(GAM_POW[t] ** (1.0 / t)) * (2 ** 30) ** (1.0 / t))
                print(f"{name} t={t} true d={d:.4e} k={arg[0]} s={arg[1]} ratio={d*math.sqrt(n30):.4f} S_true={St:.4f}")
