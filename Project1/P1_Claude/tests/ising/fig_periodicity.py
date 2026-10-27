"""Exact periodicity of Metropolis chains driven by RANDU / RANDM: batch b
(16384 sweeps) versus batch b+128 (one generator period 2^29 draws later)."""
import sys, os
import numpy as np
sys.path.insert(0, "/Users/noah/Downloads/ComputationalPhysicsPHYS25000/Project1/P1_Claude")
sys.path.insert(0, "/Users/noah/Downloads/ComputationalPhysicsPHYS25000/Project1/P1_Claude/tests/ising")
import plotstyle as ps
import matplotlib.pyplot as plt
from analyze import read, EXACT
ps.apply()
HERE = "/Users/noah/Downloads/ComputationalPhysicsPHYS25000/Project1/P1_Claude/tests/ising"
e0 = EXACT[16][0]
fig, axes = plt.subplots(1, 3, figsize=(10, 3.7), sharex=True, sharey=True)
for ax, (fn, g) in zip(axes, [("B_m_mt", "MT19937"), ("B_m_randu", "RANDU"), ("B_m_randm", "RANDM")]):
    hdr, rows, _, _ = read(os.path.join(HERE, "runs", fn + ".dat"))
    eb = rows[:, 2].astype(float) / rows[:, 1].astype(float) / 256 - e0
    P = 128
    x, y = eb[:-P], eb[P:]
    same = np.array([all(rows[b, 2:] == rows[b + P, 2:]) for b in range(len(rows) - P)])
    ax.plot(x * 1e3, y * 1e3, "o", ms=3.2, color=ps.COLORS[g], alpha=0.7, mew=0)
    lim = 20
    ax.plot([-lim, lim], [-lim, lim], ls="--", color=ps.EXPECT, lw=1.0)
    ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim)
    ax.set_title(g)
    r = np.corrcoef(x, y)[0, 1]
    ax.text(0.04, 0.96, f"identical pairs: {same.sum()} of {len(same)}\ncorrelation {r:.3f}", transform=ax.transAxes,
            va="top", fontsize=8.5, color=ps.INK2, bbox=dict(fc="white", ec="none", alpha=0.85))
    ax.set_xlabel(r"batch $b$: $(\bar e_b - e_{\rm exact}) \times 10^3$")
    ax.set_aspect("equal")
axes[0].set_ylabel(r"batch $b+128$: $(\bar e - e_{\rm exact}) \times 10^3$")
print(ps.save(fig, "ising_metropolis_periodicity"))
