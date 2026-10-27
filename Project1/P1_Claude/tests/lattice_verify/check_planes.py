import sys, math
import numpy as np
from scipy import stats
sys.path.insert(0, '/Users/noah/Downloads/ComputationalPhysicsPHYS25000/Project1/P1_Claude/tests/lattice_verify')
from lcg_np import randu, randm, minstd

M = 1 << 31
print('RANDU seed 1:', randu(1, 3))
print('MINSTD z_10000 seed 1:', minstd(1, 10000)[-1])

N = 10 ** 7
X = randu(777, 3 * N + 2)
# overlapping
S_ov = 9 * X[:-2] - 6 * X[1:-1] + X[2:]
print('overlapping divisible by 2^31:', np.all(S_ov % M == 0), 'count', len(S_ov))
T = X[:3 * N].reshape(N, 3)
S = 9 * T[:, 0] - 6 * T[:, 1] + T[:, 2]
print('non-overlapping divisible:', np.all(S % M == 0))
j = S // M
cnt = np.bincount(j + 5, minlength=15)
w = np.array([1, 3, 5, 7, 9, 11, 12, 12, 12, 11, 9, 7, 5, 3, 1]) / 108
E = N * w
chi2 = ((cnt - E) ** 2 / E).sum()
print('planes j=-5..9 counts', cnt.tolist())
print('chi2 non-overlap', chi2, 'p', stats.chi2.sf(chi2, 14), 'max|z|', np.max(np.abs((cnt - E) / np.sqrt(E * (1 - w)))),
      'argmax j', np.argmax(np.abs((cnt - E) / np.sqrt(E * (1 - w)))) - 5)
jo = S_ov[:N] // M
cnto = np.bincount(jo + 5, minlength=15)
chi2o = ((cnto - E) ** 2 / E).sum()
print('chi2 overlapping (first 1e7)', chi2o, stats.chi2.sf(chi2o, 14))

# exact plane shares by Monte Carlo of the density f_S at integers: f_S(j) via exact piecewise formula
# f_S = density of 9U1 - 6U2 + U3; compute by numerical convolution check at integers
from fractions import Fraction
def F3(x):
    # CDF of 9U1-6U2+U3 by inclusion-exclusion: P(sum c_i U_i <= x), c=(9,-6,1)
    # transform: -6U2 = -6 + 6U2' -> S = -6 + 9U1 + 6U2' + U3 with all positive coefficients
    c = [9, 6, 1]
    y = x + 6
    tot = Fraction(0)
    import itertools
    for sub in itertools.product([0, 1], repeat=3):
        z = y - sum(ci for ci, si in zip(c, sub) if si)
        if z > 0:
            tot += (-1) ** sum(sub) * z ** 3
    return tot / (6 * 9 * 6 * 1)
def f3(x):
    c = [9, 6, 1]; y = x + 6
    tot = Fraction(0)
    import itertools
    for sub in itertools.product([0, 1], repeat=3):
        z = y - sum(ci for ci, si in zip(c, sub) if si)
        if z > 0:
            tot += (-1) ** sum(sub) * 3 * z ** 2
    return tot / (6 * 54)
pj = [f3(Fraction(k)) for k in range(-5, 10)]
print('exact f_S(j)*108:', [p * 108 for p in pj], 'sum', sum(pj))

# single precision residuals
x = (X[:3 * N] / M).astype(np.float32).astype(np.float64).reshape(N, 3)
Sf = 9 * x[:, 0] - 6 * x[:, 1] + x[:, 2]
r = np.abs(Sf - np.rint(Sf))
print('float32: max resid', r.max(), 'frac exactly 0', np.mean(r == 0), 'median', np.median(r), 'all <0.005', np.all(r < 0.005))
xd = (X[:3 * N] / M).reshape(N, 3)
Sd = 9 * xd[:, 0] - 6 * xd[:, 1] + xd[:, 2]
print('float64 resid max', np.abs(Sd - np.rint(Sd)).max())
# MT control
rs = np.random.RandomState(777)
u = rs.random_sample(3 * N).reshape(N, 3)
Sm = 9 * u[:, 0] - 6 * u[:, 1] + u[:, 2]
fr = np.mean(np.abs(Sm - np.rint(Sm)) < 0.005)
print('MT19937(RandomState 777) near-integer fraction', fr, 'z', (fr - 0.01) / math.sqrt(0.01 * 0.99 / N))
