# Calibrate the birthday-spacings null distribution (m=512, n=2^24, cyclic
# spacings, J = #spacings - #distinct spacing values) with two numpy bit
# generators unrelated to anything in rnglib (Philox4x64 and PCG64).
import numpy as np, json
from scipy import stats
m, n = 512, 2**24
def simulate(bitgen, R, chunk=8192):
    rng = np.random.Generator(bitgen)
    counts = np.zeros(8, np.int64)
    for start in range(0, R, chunk):
        d = np.sort(rng.integers(0, n, size=(chunk, m)), axis=1)
        sp = np.empty_like(d)
        sp[:, 1:] = np.diff(d, axis=1)
        sp[:, 0] = d[:, 0] + n - d[:, -1]
        sp.sort(axis=1)
        J = (np.diff(sp, axis=1) == 0).sum(axis=1)
        counts += np.bincount(np.minimum(J, 7), minlength=8)
    return counts
R = 2**19
res = {}
for name, bg in [('Philox', np.random.Philox(20251023)), ('PCG64', np.random.PCG64(777))]:
    c = simulate(bg, R)
    res[name] = c.tolist()
    print(name, c, 'mean J %.5f' % ((c*np.arange(8)).sum()/c.sum()))
c1, c2 = np.array(res['Philox']), np.array(res['PCG64'])
tab = np.vstack([c1, c2])
chi2, p, dof, _ = stats.chi2_contingency(tab)
print('homogeneity Philox vs PCG64: chi2 %.2f dof %d p %.3f' % (chi2, dof, p))
pois = stats.poisson.pmf(np.arange(7), 2.0); pois = np.append(pois, 1 - pois.sum())
tot = c1 + c2
e = tot.sum()*pois
print('pooled vs Poisson(2): chi2 %.1f p %.2e' % (((tot-e)**2/e).sum(), stats.chi2.sf(((tot-e)**2/e).sum(), 7)))
probs = tot / tot.sum()
json.dump({'probs': probs.tolist(), 'counts': tot.tolist(), 'R': int(tot.sum()),
           'mean_J': float((tot*np.arange(8)).sum()/tot.sum()), 'poisson2': pois.tolist()},
          open('out/bday_null.json', 'w'), indent=1)
