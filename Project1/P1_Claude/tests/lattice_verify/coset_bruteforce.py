"""Brute-force check of the coset gap formula d_t = max_k h_e(k)/nu_t(a,2^k)
for small moduli 2^e, a = 3 mod 8, odd seeds: compute the empirical largest cyclic gap
of s.X mod 2^e over all tuples of the full period, for every integer s in a box large
enough to contain every s that could beat the formula value."""
import math, itertools, sys
import numpy as np
sys.path.insert(0, '/Users/noah/Downloads/ComputationalPhysicsPHYS25000/Project1/P1_Claude/tests/lattice_verify')
from spec_verify import shortest


def h_e(k, e):
    if k >= e - 1:
        return 1.0
    if k == e - 2:
        return 0.5
    return 6.0 * 2 ** k / 2 ** e


def formula(a, e, t):
    best = 0
    for k in range(1, e):
        n2, sv, _ = shortest(a, 2 ** k, t)
        best = max(best, h_e(k, e) / math.sqrt(n2))
    return best


def brute(a, e, t, x0):
    M = 2 ** e
    xs = [x0]
    while True:
        y = a * xs[-1] % M
        if y == x0:
            break
        xs.append(y)
    P = len(xs)
    X = np.array([[xs[(i + j) % P] for j in range(t)] for i in range(P)], dtype=np.int64)
    d0 = formula(a, e, t)
    B = int(math.ceil(1.0 / d0)) + 1
    best = 0.0
    rng = range(-B, B + 1)
    for s in itertools.product(rng, repeat=t):
        if not any(s):
            continue
        n = math.sqrt(sum(c * c for c in s))
        if 1.0 / n <= best:
            continue
        r = np.unique((X @ np.array(s, dtype=np.int64)) % M)
        gaps = np.diff(np.concatenate([r, [r[0] + M]]))
        g = gaps.max() / M / n
        best = max(best, g)
    return P, d0, best


if __name__ == '__main__':
    mism = 0
    ncase = 0
    for e in (10, 11, 12):
        for a in (2 ** (e // 2) + 3, 6539 % 2 ** e, 65539 % 2 ** e, 2 ** e - 5, 11, 27):
            if a % 8 != 3:
                continue
            for t in (2, 3):
                for x0 in (1, 7, 777 % 2 ** e):
                    P, d0, b = brute(a, e, t, x0)
                    ok = abs(d0 - b) < 1e-12
                    ncase += 1
                    mism += (not ok)
                    print(f"e={e} a={a} t={t} x0={x0} P={P} formula={d0:.6e} brute={b:.6e} {'OK' if ok else 'MISMATCH'}")
    print('cases', ncase, 'mismatches', mism)
