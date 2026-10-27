import sys, json
sys.path.insert(0, "/Users/noah/Downloads/ComputationalPhysicsPHYS25000/Project1/P1_Claude")
import numpy as np
import plotstyle as ps
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, TwoSlopeNorm
from scipy import stats
ps.apply()
D = "/Users/noah/Downloads/ComputationalPhysicsPHYS25000/Project1/P1_Claude/tests/lattice/"
out = {}
# ------------------------------------------------ diffraction
d = np.loadtxt(D + 'structure_factor.txt')
n = d[:, :3].astype(int)
names = ['RANDU', 'RANDM', 'MINSTD', 'MT19937', 'RANNUM']
u = (n[:, 2] > 0) | ((n[:, 2] == 0) & ((n[:, 1] > 0) | ((n[:, 1] == 0) & (n[:, 0] > 0))))
NP = 1e5
fig = plt.figure(figsize=(10, 3.9))
gs = fig.add_gridspec(1, 5, width_ratios=[1, 1, 0.05, 0.32, 1.3], wspace=0.25)
for i, g in enumerate(['RANDU', 'MT19937']):
    ax = fig.add_subplot(gs[0, i])
    k = names.index(g)
    sl = n[:, 2] == 1
    grid = np.zeros((21, 21))
    for (a, b, c), S in zip(n[sl], d[sl, 3 + k]):
        grid[b + 10, a + 10] = S
    im = ax.imshow(grid, origin='lower', extent=(-10.5, 10.5, -10.5, 10.5), cmap='cividis',
                   norm=LogNorm(vmin=1e-2, vmax=1e5))
    ax.set_xlabel('$n_1$'); ax.set_ylabel('$n_2$'); ax.grid(False)
    ax.set_title('%s gas, $S(n_1,n_2,1)$' % g, fontsize=10)
    if g == 'RANDU':
        ax.annotate('$S(9,-6,1)=N$', xy=(9, -6), xytext=(-2, -9.3), color='white', fontsize=8.5,
                    arrowprops=dict(arrowstyle='->', color='white', lw=0.8))
cax = fig.add_subplot(gs[0, 2])
cb = fig.colorbar(im, cax=cax)
cb.set_label('structure factor $S(\\mathbf{n})$')
ax = fig.add_subplot(gs[0, 4])
ss = np.logspace(-3, 5.3, 400)
for g in ['RANDU', 'RANDM', 'MINSTD', 'MT19937', 'RANNUM']:
    S = np.sort(d[u, 3 + names.index(g)])
    surv = 1 - np.arange(len(S)) / len(S)
    ax.step(S, surv, where='post', color=ps.COLORS[g], lw=1.2, label=g)
    out['S_max_' + g] = float(S[-1]); out['S_second_' + g] = float(S[-2])
    # KS test against Exp(1) excluding Bragg peaks (S > 1000)
    Sb = S[S < 1000]
    out['KS_p_' + g] = float(stats.kstest(Sb, 'expon').pvalue)
ax.plot(ss, np.exp(-ss), color=ps.EXPECT, ls='--', lw=1, label='ideal gas, $e^{-s}$')
ax.set_xscale('log'); ax.set_yscale('log'); ax.set_ylim(1.5e-4, 1.5); ax.set_xlim(1e-3, 3e5)
ax.set_xlabel('$s$'); ax.set_ylabel('fraction of wave vectors with $S > s$')
ax.set_title('All 4630 wave vectors, $|n_i| \\leq 10$', fontsize=10)
ax.legend(fontsize=7.5, loc='lower left')
ax.annotate('RANDU Bragg peak\n$S = 10^5$', xy=(1e5, 1 / 4630), xytext=(1.2e2, 1.2e-3), fontsize=8,
            arrowprops=dict(arrowstyle='->', color=ps.INK2, lw=0.8))
p = ps.save(fig, 'lattice_diffraction'); print(p)
# ------------------------------------------------ isotropy
h = np.loadtxt(D + 'isotropy_costheta.txt')
x = h[:, 0]; NB = len(x); ND = 2e6
meth = ['RANDU, cube rejection', 'MT19937, cube rejection', 'RANDU, $(\\cos\\theta,\\phi)$ from pairs', 'RANDM, cube rejection']
for j in range(4):
    c = h[:, 1 + j]; E = ND / NB
    chi2 = ((c - E) ** 2 / E).sum()
    out['iso_chi2_%d' % j] = float(chi2); out['iso_p_%d' % j] = float(stats.chi2.sf(chi2, NB - 1))
    # same without the central bins (one bin either side of 0)
    m = np.abs(x) > 0.006
    chi2b = ((c[m] - E) ** 2 / E).sum()
    out['iso_chi2_nocenter_%d' % j] = float(chi2b); out['iso_p_nocenter_%d' % j] = float(stats.chi2.sf(chi2b, m.sum() - 1))
sky = np.loadtxt(D + 'isotropy_sky.txt')
for j in range(4):
    c = sky[:, 2 + j]; E = ND / len(c)
    chi2 = ((c - E) ** 2 / E).sum()
    out['sky_chi2_%d' % j] = float(chi2); out['sky_p_%d' % j] = float(stats.chi2.sf(chi2, len(c) - 1))
ex = np.loadtxt(D + 'isotropy_exact.txt'); out['frac_exact'] = float(ex[0, 1] / ex[0, 2])
mband = (np.abs(x) > 0.005) & (np.abs(x) < 0.18)
out['band_counts'] = [float(h[mband, 1 + j].sum()) for j in range(4)]
out['band_expected'] = float(ND * mband.sum() / NB)
fig = plt.figure(figsize=(10, 7.6))
gs = fig.add_gridspec(2, 3, width_ratios=[1, 1, 0.045], hspace=0.38, wspace=0.28)
# (a) pi from the ball volume
ax = fig.add_subplot(gs[0, 0])
rows = [l.split() for l in open(D + 'pi_running.txt') if not l.startswith('#')]
for g, off in (('RANDU', 1.0), ('RANDM', 0.9), ('MINSTD', 0.95), ('MT19937', 1.05), ('RANNUM', 1.1)):
    a = np.array([[float(r[1]), float(r[2]), float(r[3])] for r in rows if r[0] == g])
    a = a[a[:, 0] >= 1e4]
    ax.errorbar(a[:, 0] * off, a[:, 1] - np.pi, yerr=a[:, 2], fmt='o', ms=2.8, lw=0.9, color=ps.COLORS[g], label=g)
    out['pi_final_' + g] = [float(a[-1, 1]), float(a[-1, 2]), float((a[-1, 1] - np.pi) / a[-1, 2])]
dk = 2 * np.arange(-5, 6) / np.sqrt(118); W = (1 - dk ** 2).sum()
pi_lat = 6 * np.pi * W * (2 / np.sqrt(118)) / 8
out['pi_lattice_theory'] = float(pi_lat)
ax.axhline(pi_lat - np.pi, color=ps.COLORS['RANDU'], ls=':', lw=1.0, label='RANDU 15-plane quadrature')
nn = np.logspace(3.7, 8.6, 50)
ax.plot(nn, 2 * 3.0 / np.sqrt(nn) * np.sqrt(np.pi / 6 * (1 - np.pi / 6)) / 0.5, color=ps.EXPECT, ls='--', lw=0.9, label='$\\pm 2\\sigma$')
ax.plot(nn, -2 * 3.0 / np.sqrt(nn) * np.sqrt(np.pi / 6 * (1 - np.pi / 6)) / 0.5, color=ps.EXPECT, ls='--', lw=0.9)
ax.axhline(0, color=ps.INK2, lw=0.6)
ax.set_xscale('log'); ax.set_ylim(-0.035, 0.035); ax.set_xlim(5e3, 5e8)
ax.set_xlabel('number of triplets $N$'); ax.set_ylabel('$\\hat\\pi - \\pi$')
ax.set_title('$\\pi$ from the volume of the unit ball', fontsize=10)
ax.legend(fontsize=7, loc='lower right', ncol=2, columnspacing=0.8)
# (b) cos theta density
ax = fig.add_subplot(gs[0, 1:3])
dens = lambda c: c / ND / (2.0 / NB)
ax.step(x, dens(h[:, 2]), where='mid', color=ps.COLORS['MT19937'], lw=0.9, label='MT19937')
ax.step(x, dens(h[:, 1]), where='mid', color=ps.COLORS['RANDU'], lw=1.1, label='RANDU')
cc = np.linspace(-1, 1, 4001); th = np.zeros_like(cc)
for d in dk[dk > 0]:
    th += np.where(np.abs(cc) >= d, 2 * d ** 2 / (W * np.abs(cc) ** 3 + 1e-300), 0)
ax.plot(cc, th, color=ps.EXPECT, ls='--', lw=0.9, label='15-plane theory')
ax.axhline(0.5, color=ps.INK2, ls=':', lw=0.9, label='isotropic, 1/2')
ax.set_ylim(0, 1.6); ax.set_xlim(-1, 1)
ax.annotate('%.1f%% of all directions\nexactly in one plane' % (100 * out['frac_exact']), xy=(0.0, 1.6), xytext=(0.3, 1.42), fontsize=8,
            arrowprops=dict(arrowstyle='->', color=ps.INK2, lw=0.8))
ax.annotate('empty band\n$0<|\\cos\\theta|<2/\\sqrt{118}$', xy=(0.09, 0.03), xytext=(0.3, 1.08), fontsize=8,
            arrowprops=dict(arrowstyle='->', color=ps.INK2, lw=0.8))
ax.set_xlabel('$\\cos\\theta$ relative to $(9,-6,1)/\\sqrt{118}$'); ax.set_ylabel('probability density')
ax.set_title('Isotropic emission by cube rejection, $2\\times10^6$ directions', fontsize=10)
ax.legend(fontsize=7.5, loc='upper left')
for k, (j, g) in enumerate([(0, 'RANDU'), (1, 'MT19937')]):
    ax = fig.add_subplot(gs[1, k])
    c = sky[:, 2 + j].reshape(72, 144); E = ND / c.size
    im = ax.imshow(c / E, origin='lower', extent=(-180, 180, -1, 1), aspect='auto', cmap='RdBu_r',
                   norm=TwoSlopeNorm(vmin=0, vcenter=1, vmax=2))
    ax.set_xlabel('azimuth $\\phi$ (deg)'); ax.grid(False)
    if k == 0: ax.set_ylabel('$\\cos\\theta$')
    ax.set_title('%s sky map (equal area)' % g, fontsize=10)
cax = fig.add_subplot(gs[1, 2])
cb = fig.colorbar(im, cax=cax)
cb.set_label('counts / isotropic expectation')
p = ps.save(fig, 'lattice_isotropy'); print(p)
json.dump(out, open(D + 'phys_stats.json', 'w'), indent=1)
