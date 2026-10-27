import sys
sys.path.insert(0, "/Users/noah/Downloads/ComputationalPhysicsPHYS25000/Project1/P1_Claude")
import numpy as np
import plotstyle as ps
import matplotlib.pyplot as plt
ps.apply()
D = "/Users/noah/Downloads/ComputationalPhysicsPHYS25000/Project1/P1_Claude/tests/lattice/"
EXTRA = {"MINSTD2": "#9c6b00", "LCG69069": "#8c5a2b"}
col = dict(ps.COLORS); col.update(EXTRA)
rows = [l.split() for l in open(D + 'mc_running.txt') if not l.startswith('#')]
data = {}
for g, k, n, m, e in rows:
    data.setdefault((g, int(k)), []).append((int(n), float(m), float(e)))
K = {1: ('RANDU', '$\\cos 2\\pi(9x_1-6x_2+x_3)$'),
     2: ('RANDM', '$\\cos 2\\pi(6539x_1-x_2)$'),
     3: ('RANDM', '$-\\cos 2\\pi(25x_1+115x_2+678x_3)$'),
     5: ('MINSTD', '$\\cos 2\\pi(16807x_1-x_2)$'),
     7: ('MINSTD2', '$\\cos 2\\pi(3x_1-4x_2+4x_3+4x_4+2x_5-x_6+2x_7+4x_8)$'),
     8: ('LCG69069', '$\\cos 2\\pi(3x_1-11x_2+4x_3-8x_4+4x_5-4x_6-\\varphi)$')}
fig, axs = plt.subplots(2, 3, figsize=(10, 6.6), sharex=True, sharey=True)
for ax, (k, (owner, lab)) in zip(axs.flat, K.items()):
    for g, off in ((owner, 1.0), ('MT19937', 0.93), ('RANNUM', 1.07)):
        a = np.array(data[(g, k)])
        ax.errorbar(a[:, 0] * off, a[:, 1], yerr=a[:, 2], fmt='o', ms=2.8, lw=0.9, capsize=0,
                    color=col[g], label=g)
    nn = np.logspace(1, 7, 100)
    ax.plot(nn, 2 / np.sqrt(2 * nn), color=ps.EXPECT, ls='--', lw=0.9, label='$\\pm 2/\\sqrt{2N}$')
    ax.plot(nn, -2 / np.sqrt(2 * nn), color=ps.EXPECT, ls='--', lw=0.9)
    ax.axhline(0, color=ps.INK2, lw=0.6)
    ax.set_xscale('log')
    ax.set_title('%s killer, $t=%d$' % (owner, {1:3,2:2,3:3,5:2,7:8,8:6}[k]), fontsize=10)
    ax.text(0.03, 0.06, lab, transform=ax.transAxes, fontsize=7.2, color=ps.INK2)
    ax.legend(loc='center right', fontsize=7)
    ax.set_ylim(-0.75, 1.2); ax.set_xlim(5, 2e7)
for ax in axs[1]:
    ax.set_xlabel('number of tuples $N$')
for ax in axs[:, 0]:
    ax.set_ylabel('MC estimate (exact value 0)')
p = ps.save(fig, 'lattice_mc_killer'); print(p)
