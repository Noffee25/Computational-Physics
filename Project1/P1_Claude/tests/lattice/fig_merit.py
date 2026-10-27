import sys, json
sys.path.insert(0, "/Users/noah/Downloads/ComputationalPhysicsPHYS25000/Project1/P1_Claude")
import numpy as np
import plotstyle as ps
import matplotlib.pyplot as plt
from math import factorial
ps.apply()
D = "/Users/noah/Downloads/ComputationalPhysicsPHYS25000/Project1/P1_Claude/tests/lattice/"
EXTRA = {"MINSTD2": "#9c6b00", "LCG69069": "#8c5a2b"}
col = dict(ps.COLORS); col.update(EXTRA)
rows = json.load(open(D + 'summary.json'))
gens = ['RANDU', 'RANDM', 'MINSTD', 'MINSTD2', 'LCG69069']
lab = {'RANDU': 'RANDU ($a=65539$)', 'RANDM': 'RANDM ($a=6539$)', 'MINSTD': 'MINSTD ($a=16807$)',
       'MINSTD2': 'MINSTD2 ($a=48271$)', 'LCG69069': 'LCG69069 ($a=69069$)'}
mk = {'RANDU': 'o', 'RANDM': 's', 'MINSTD': 'D', 'MINSTD2': 'v', 'LCG69069': '^'}
fig, axs = plt.subplots(1, 2, figsize=ps.WIDE)
ax = axs[0]
for g in gens:
    r = [x for x in rows if x['gen'] == g]
    t = [x['t'] for x in r]; S = [x['S'] for x in r]
    ax.plot(t, S, marker=mk[g], ms=5, color=col[g], label=lab[g])
    if g in ('RANDU', 'RANDM'):
        te = [x['t'] for x in r if abs(x['S_eff'] - x['S']) > 1e-9]
        Se = [x['S_eff'] for x in r if abs(x['S_eff'] - x['S']) > 1e-9]
        ax.plot(te, Se, ls='none', marker=mk[g], ms=7, mfc='white', mec=col[g], mew=1.3)
ax.plot([], [], ls='none', marker='o', mfc='white', mec=ps.INK2, label='true largest gap (coset correction)')
ax.set_yscale('log'); ax.set_ylim(5e-3, 1.3)
ax.set_xlabel('dimension $t$'); ax.set_ylabel('$S_t = \\nu_t / (\\gamma_t^{1/2} m^{1/t})$')
ax.set_title('Normalized spectral figure of merit')
ax.legend(fontsize=7.5, loc='lower right')
ax = axs[1]
for g in gens:
    r = [x for x in rows if x['gen'] == g]
    ax.plot([x['t'] for x in r], [x['planes'] for x in r], marker=mk[g], ms=5, color=col[g], label=lab[g])
tt = np.arange(2, 9)
ax.plot(tt, [(factorial(t) * 2 ** 31) ** (1 / t) for t in tt], color=ps.EXPECT, ls='--', lw=1.1,
        label="Marsaglia bound $(t!\\,m)^{1/t}$, $m=2^{31}$")
ax.set_yscale('log')
ax.set_xlabel('dimension $t$'); ax.set_ylabel('occupied hyperplanes (shortest dual vector)')
ax.set_title('Hyperplanes covering all $t$-tuples (full period)')
ax.legend(fontsize=7.5, loc='upper right')
p = ps.save(fig, 'lattice_spectral_merit'); print(p)
