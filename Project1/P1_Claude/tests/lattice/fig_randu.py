import sys, json
sys.path.insert(0, "/Users/noah/Downloads/ComputationalPhysicsPHYS25000/Project1/P1_Claude")
import numpy as np
import plotstyle as ps
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
ps.apply()
D = "/Users/noah/Downloads/ComputationalPhysicsPHYS25000/Project1/P1_Claude/tests/lattice/"
smp = np.loadtxt(D + 'triplet_sample.txt')
ru, mt = smp[:, :3], smp[:, 3:]
nshow = 12000
fig = plt.figure(figsize=(10, 7.2))
gs = GridSpec(2, 6, figure=fig, height_ratios=[1.05, 1.0], hspace=0.22, wspace=0.9)
# viewing direction d inside the planes (d . (9,-6,1) = 0) with 20 degrees elevation
e1 = np.array([2, 3, 0]) / np.sqrt(13); e2 = np.array([-3, 2, 39]) / np.sqrt(1534)
beta = np.arcsin(np.sin(np.radians(20)) / e2[2]); dv = np.cos(beta) * e1 + np.sin(beta) * e2
assert abs(dv @ np.array([9, -6, 1])) < 1e-12
azim_edge = np.degrees(np.arctan2(dv[1], dv[0])); elev_edge = np.degrees(np.arcsin(dv[2]))
print('edge-on view azim %.3f elev %.3f' % (azim_edge, elev_edge))
views = [(ru, 'RANDU', 22, -58, 'RANDU, generic view'),
         (ru, 'RANDU', elev_edge, azim_edge, 'RANDU, edge-on'),
         (mt, 'MT19937', elev_edge, azim_edge, 'MT19937, same view')]
for i, (P, g, el, az, title) in enumerate(views):
    ax = fig.add_subplot(gs[0, 2 * i:2 * i + 2], projection='3d', proj_type='ortho')
    ax.scatter(P[:nshow, 0], P[:nshow, 1], P[:nshow, 2], s=0.4, color=ps.COLORS[g], alpha=0.6, linewidths=0, rasterized=True)
    ax.view_init(elev=el, azim=az)
    ax.set_box_aspect((1, 1, 1))
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.set_zlim(0, 1)
    ax.set_xlabel('$x_k$', labelpad=-8); ax.set_ylabel('$x_{k+1}$', labelpad=-8); ax.set_zlabel('$x_{k+2}$', labelpad=-8)
    ax.tick_params(labelsize=7, pad=-3)
    ax.set_title(title, pad=0)
    for a in (ax.xaxis, ax.yaxis, ax.zaxis):
        a.pane.set_alpha(0.0)
# ---- histogram of S = 9x-6y+z
h = np.loadtxt(D + 'triplet_hist.txt')
xc = h[:, 0]; N = 1e7
ax = fig.add_subplot(gs[1, 0:3])
fr_mt = h[:, 2] / N
ax.plot(xc, fr_mt, color=ps.COLORS['MT19937'], lw=0.8, label='MT19937')
# theory density of 9U-6V+W times bin width
from itertools import product
def fS(x):
    x = np.asarray(x) + 6.0; c = (9, 6, 1); tot = np.zeros_like(x)
    for eps in product((0, 1), repeat=3):
        sh = sum(e * ci for e, ci in zip(eps, c)); d = np.clip(x - sh, 0, None)
        tot += (-1) ** sum(eps) * d ** 2
    return tot / 108.0
xx = np.linspace(-6, 10, 2000)
ax.plot(xx, 0.01 * fS(xx), color=ps.EXPECT, ls='--', lw=1.0, label='uniform theory, $0.01 f_S$')
sel = h[:, 1] > 0
ax.vlines(xc[sel], 0, h[sel, 1] / N, color=ps.COLORS['RANDU'], lw=2.0, label='RANDU')
pj = json.load(open(D + 'planes_expect.json'))['p']
J = np.array([int(k) for k in pj]); P = np.array([pj[k] for k in pj])
ax.plot(J, P, 'o', mfc='none', mec=ps.EXPECT, ms=6, label='expected plane share $f_S(j)$')
ax.set_yscale('log'); ax.set_ylim(3e-5, 4.0); ax.set_xlim(-6.5, 10.5)
ax.set_xlabel('$S = 9x_k - 6x_{k+1} + x_{k+2}$'); ax.set_ylabel('fraction of triplets per 0.01 bin')
ax.set_title('Distribution of $9x_k-6x_{k+1}+x_{k+2}$ ($10^7$ triplets)')
ax.legend(loc='upper center', fontsize=8, ncol=2)
# ---- residuals of plane populations
d = np.loadtxt(D + 'triplet_planes.txt'); m = (d[:, 0] >= -5) & (d[:, 0] <= 9)
j = d[m, 0]; cnt = d[m, 2]; Nn = cnt.sum(); E = Nn * P
z = (cnt - E) / np.sqrt(E * (1 - P))
chi2 = ((cnt - E) ** 2 / E).sum()
from scipy import stats
ax = fig.add_subplot(gs[1, 3:6])
ax.bar(j, z, color=ps.COLORS['RANDU'], width=0.7)
ax.axhline(2, color=ps.EXPECT, ls='--', lw=1); ax.axhline(-2, color=ps.EXPECT, ls='--', lw=1)
ax.set_xlabel('plane index $j$'); ax.set_ylabel('(observed $-$ expected) / $\\sigma$')
ax.set_xticks(range(-5, 10))
ax.set_ylim(-3.5, 3.5)
ax.set_title('Plane populations vs geometry')
ax.text(0.03, 0.95, '$\\chi^2 = %.1f$, 14 dof, $p = %.2f$\nexpected shares $(1,3,5,7,9,11,12,12,12,\\dots)/108$' % (chi2, stats.chi2.sf(chi2, 14)),
        transform=ax.transAxes, va='top', fontsize=8.5)
print('chi2', chi2, stats.chi2.sf(chi2, 14), 'max|z|', abs(z).max())
p = ps.save(fig, 'lattice_randu_planes'); print(p)
