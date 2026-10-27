import numpy as np, json
from scipy import stats
D = "/Users/noah/Downloads/ComputationalPhysicsPHYS25000/Project1/P1_Claude/tests/lattice/"
c = np.loadtxt(D + 'ball_blocks.txt')
names = ['RANDU', 'RANDM', 'MINSTD', 'MINSTD2', 'LCG69069', 'MT19937', 'RANNUM']
B0 = 10000
res = {}
for j, g in enumerate(names):
    x = c[:, j]
    p = x.sum() / (len(x) * B0)
    res[g] = {}
    line = '%-9s p=%.6f ' % (g, p)
    for f in (1, 10, 100, 1000):
        nb = len(x) // f
        y = x[:nb * f].reshape(nb, f).sum(1)
        B = B0 * f
        D_ = y.var(ddof=1) / (B * p * (1 - p))
        chi = (nb - 1) * D_
        pv = stats.chi2.sf(chi, nb - 1)
        res[g][B] = dict(nblocks=nb, dispersion=D_, p_upper=pv, p_lower=stats.chi2.cdf(chi, nb - 1))
        line += '| B=%.0e nb=%d D=%.3f p>=%.2g ' % (B, nb, D_, pv)
    # z of the first 1e7 block relative to pi/6
    first = x[:1000].sum(); N = 1e7; q = np.pi / 6
    res[g]['z_first_1e7'] = (first / N - q) / np.sqrt(q * (1 - q) / N)
    zb = (x[:17000].reshape(17, 1000).sum(1) / N - q) / np.sqrt(q * (1 - q) / N)
    res[g]['z_blocks_1e7'] = zb.tolist()
    line += '| z(first 1e7)=%.2f' % res[g]['z_first_1e7']
    print(line)
for g in names:
    print(g, ' '.join('%.1f' % z for z in res[g]['z_blocks_1e7']))
json.dump(res, open(D + 'ball_blocks.json', 'w'), indent=1, default=float)
