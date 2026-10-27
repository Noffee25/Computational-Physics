import sys
sys.path.insert(0, "/Users/noah/Downloads/ComputationalPhysicsPHYS25000/Project1/P1_Claude")
import numpy as np
import plotstyle as ps
import matplotlib.pyplot as plt
ps.apply()
D = "/Users/noah/Downloads/ComputationalPhysicsPHYS25000/Project1/P1_Claude/tests/lattice/"
fig, axs = plt.subplots(2, 2, figsize=(10, 7.4))
lim = 3 / 6539
for ax, g in zip(axs[0], ['RANDM', 'MT19937']):
    d = np.loadtxt(D + 'pairs_origin_%s.txt' % g)
    ax.scatter(d[:, 0] * 1e4, d[:, 1], s=0.6, color=ps.COLORS[g], alpha=0.6, linewidths=0, rasterized=True)
    ax.set_xlim(0, lim * 1e4); ax.set_ylim(0, 1)
    ax.set_xlabel('$x_k$  ($\\times 10^{-4}$)'); ax.set_ylabel('$x_{k+1}$')
    ax.set_title('%s pairs with $x_k < 3/6539$ (%d of $10^8$)' % (g, len(d)), fontsize=10)
    ax.grid(False)
    if g == 'RANDM':
        fig.canvas.draw()
        p0 = ax.transData.transform((0, 0)); p1 = ax.transData.transform((1e4 / 6539, 1))
        ang = np.degrees(np.arctan2(p1[1] - p0[1], p1[0] - p0[0]))
        for j in range(3):
            lab = '$x_{k+1}=6539\\,x_k$' if j == 0 else '$x_{k+1}=6539\\,x_k-%d$' % j
            ax.text((j + 0.5) / 6539 * 1e4 + 0.24, 0.5, lab, rotation=ang, fontsize=8,
                    color=ps.INK2, ha='center', va='center')
sn = np.sqrt(25 ** 2 + 115 ** 2 + 678 ** 2)
for ax, g in zip(axs[1], ['RANDM', 'MT19937']):
    d = np.loadtxt(D + 'randm3_slab_%s.txt' % g)
    ax.scatter((d[:, 0] - 0.6) * 1e3, d[:, 1], s=0.5, color=ps.COLORS[g], alpha=0.6, linewidths=0, rasterized=True)
    ax.set_xlabel('$(\\mathbf{s}\\cdot\\mathbf{x}/|\\mathbf{s}| - 0.6)\\times 10^{3}$,  $\\mathbf{s}=(25,115,678)$')
    ax.set_ylabel('orthogonal coordinate')
    ax.set_title('%s triplets in a slab of width $8/|\\mathbf{s}|$ (%d of $3\\times10^6$)' % (g, len(d)), fontsize=10)
    ax.set_xlim(0, 8 / sn * 1e3)
    ax.grid(False)
p = ps.save(fig, 'lattice_randm_structure'); print(p)
