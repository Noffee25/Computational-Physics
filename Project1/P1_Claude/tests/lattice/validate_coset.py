"""Brute-force check of coset_spectral (a = 3 mod 8, m = 2^e, odd seed) on small e."""
import sys; sys.path.insert(0, "/Users/noah/Downloads/ComputationalPhysicsPHYS25000/Project1/P1_Claude/tests/lattice")
import itertools, numpy as np
from math import sqrt
from spectral import coset_spectral, shortest
import random
random.seed(3)
nbad = 0; n = 0
for e in (9, 10, 11, 12):
    m = 1 << e
    for trial in range(3):
        a = random.randrange(0, m // 8) * 8 + 3
        x0 = random.randrange(0, m // 2) * 2 + 1
        # full period of states
        xs = [x0]
        while True:
            nx = (a * xs[-1]) % m
            if nx == x0: break
            xs.append(nx)
        xs = np.array(xs, dtype=np.int64)
        assert len(xs) == m // 4
        for t in (2, 3, 4):
            T = np.stack([np.roll(xs, -j) for j in range(t)], axis=1)  # cyclic tuples
            res = coset_spectral(a, e, t, x0 % 8 if False else x0)
            nu30 = sqrt(res['nu2_m1'])
            R = int(nu30) + 1
            best = 0.0; bestv = None
            for s in itertools.product(range(-R, R + 1), repeat=t):
                if not any(s): continue
                ns = sqrt(sum(c * c for c in s))
                if ns > nu30 + 1e-9: continue
                vals = np.unique((T @ np.array(s, dtype=np.int64)) % m)
                gaps = np.diff(np.concatenate([vals, [vals[0] + m]]))
                g = gaps.max() / m / ns
                if g > best + 1e-15:
                    best, bestv = g, s
            n += 1
            ok = abs(best - res['dmax']) < 1e-12
            if not ok:
                nbad += 1
            print(e, a, x0, t, 'brute %.6f %s' % (best, bestv), 'exact %.6f %s h=%s' % (res['dmax'], res['dmax_vec'], res['dmax_h']),
                  'uniform %.6f' % res['uniform_d'], 'OK' if ok else 'MISMATCH')
print('cases', n, 'mismatches', nbad)
