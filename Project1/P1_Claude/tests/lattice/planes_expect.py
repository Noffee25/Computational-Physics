"""Expected population of each RANDU plane 9x-6y+z = j from the geometry."""
import numpy as np, itertools, json
from fractions import Fraction
from scipy import stats
c = (9, 6, 1)   # S + 6 = 9U1 + 6(1-U2) + U3
def fT(x):
    tot = Fraction(0)
    for eps in itertools.product((0, 1), repeat=3):
        sh = sum(e * ci for e, ci in zip(eps, c))
        d = x - sh
        if d > 0:
            tot += (-1) ** sum(eps) * Fraction(d) ** 2
    return tot / (2 * 9 * 6 * 1)
p = {j: fT(j + 6) for j in range(-6, 11)}
print('sum p_j =', sum(p.values()))
for j in range(-5, 10):
    print(j, p[j], float(p[j]))
d = np.loadtxt('triplet_planes.txt')
J = d[:, 0].astype(int); ov = d[:, 1]; no = d[:, 2]
mask = (J >= -5) & (J <= 9)
pe = np.array([float(p[j]) for j in J[mask]])
out = {}
for name, cnt in (('overlapping', ov[mask]), ('nonoverlapping', no[mask])):
    N = cnt.sum(); E = N * pe
    chi2 = ((cnt - E) ** 2 / E).sum()
    pval = stats.chi2.sf(chi2, 14)
    z = (cnt - E) / np.sqrt(E * (1 - pe))
    print(name, 'N=%d chi2=%.2f dof=14 p=%.4f maxabs z=%.2f' % (N, chi2, pval, np.abs(z).max()))
    out[name] = dict(N=int(N), chi2=chi2, p=pval, z=z.tolist())
# full period (2^29) exact populations from lat_verify
full = {-5:4970421,-4:14912775,-3:24855137,-2:34797190,-1:44739543,0:54681904,1:59652325,2:59652323,3:59652324,4:54681901,5:44739550,6:34797187,7:24855133,8:14912780,9:4970419}
N = sum(full.values()); print('full period N =', N, 2**29)
for j in range(-5, 10):
    E = N * p[j]
    print(j, full[j], float(E), '%.2f' % (full[j] - float(E)), 'rel %.2e' % ((full[j] - float(E)) / float(E)))
out['p'] = {j: float(p[j]) for j in range(-5, 10)}
out['p_exact'] = {j: str(p[j]) for j in range(-5, 10)}
out['full'] = full
json.dump(out, open('planes_expect.json', 'w'), indent=1)
