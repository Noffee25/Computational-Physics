import sys
sys.path.insert(0, "/Users/noah/Downloads/ComputationalPhysicsPHYS25000/Project1/P1_Claude")
import numpy as np
import plotstyle as ps
import matplotlib.pyplot as plt
ps.apply()
D = "/Users/noah/Downloads/ComputationalPhysicsPHYS25000/Project1/P1_Claude/tests/lattice/"
EXTRA = {"MINSTD2": "#9c6b00", "LCG69069": "#8c5a2b"}
col = dict(ps.COLORS); col.update(EXTRA)
panels = [('RANDU_t3', 'RANDU', 't=3', '(9,-6,1)', '15 planes'),
          ('RANDM_t6', 'RANDM', 't=6', '(3,6,-8,-14,9,-4)', '43 planes'),
          ('MINSTD_t8', 'MINSTD', 't=8', '(5,4,8,1,7,-2,0,1)', '27 planes'),
          ('MINSTD2_t8', 'MINSTD2', 't=8', '(3,-4,4,4,2,-1,2,4)', '23 planes'),
          ('LCG69069_t6', 'LCG69069', 't=6', '(3,-11,4,-8,4,-4)', '34 planes'),
          ('MT_t6', 'MT19937', 't=6', 'RANDM normal (3,6,-8,-14,9,-4)', 'no planes')]
label = {'MINSTD2': 'MINSTD2 ($a=48271$)', 'LCG69069': 'LCG69069 ($a=69069$, $c=1$)'}
fig, axs = plt.subplots(2, 3, figsize=(10, 6.4))
for ax, (f, g, t, sv, npl) in zip(axs.flat, panels):
    d = np.loadtxt(D + 'edgeon_' + f + '.txt')
    ax.scatter(d[:, 0], d[:, 1], s=0.5, color=col[g], alpha=0.7, linewidths=0, rasterized=True)
    ax.set_title('%s, %s: %s' % (label.get(g, g), t, npl), fontsize=10)
    ax.text(0.02, 0.98, 's = ' + sv, transform=ax.transAxes, va='top', fontsize=7.5, color=ps.INK2)
    ax.set_xlabel('$\\mathbf{s}\\cdot\\mathbf{x}/|\\mathbf{s}|$ (distance along normal)')
    ax.set_ylabel('orthogonal coordinate')
    ax.grid(False)
p = ps.save(fig, 'lattice_edgeon_gallery'); print(p)
