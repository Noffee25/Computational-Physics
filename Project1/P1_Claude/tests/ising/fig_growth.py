"""Figure: running estimate of e and C versus number of Wolff clusters,
R250 versus MT19937, with the +-2 sigma band expected for an unbiased run."""
import sys, os
import numpy as np
sys.path.insert(0, "/Users/noah/Downloads/ComputationalPhysicsPHYS25000/Project1/P1_Claude")
sys.path.insert(0, "/Users/noah/Downloads/ComputationalPhysicsPHYS25000/Project1/P1_Claude/tests/ising")
import plotstyle as ps
import matplotlib.pyplot as plt
from analyze import read, analyze, EXACT
ps.apply()
HERE = "/Users/noah/Downloads/ComputationalPhysicsPHYS25000/Project1/P1_Claude/tests/ising"
e0, C0 = EXACT[16]
N = 256; K = 0.44068679350977151

def running(fn):
    hdr, rows, _, _ = read(fn)
    n = rows[:, 1].astype(np.float64); sB = rows[:, 2].astype(np.float64); sB2 = rows[:, 3].astype(np.float64)
    cn = np.cumsum(n); cB = np.cumsum(sB); cB2 = np.cumsum(sB2)
    e = cB / cn / N
    C = K * K * (cB2 / cn - (cB / cn) ** 2) / N
    return cn, e, C

fig, axes = plt.subplots(1, 2, figsize=ps.WIDE)
r_mt = analyze(os.path.join(HERE, "runs/A_w_mt.dat"))
# expected 2 sigma envelope from the MT run's integrated autocorrelation (sigma ~ n^-1/2)
cn = None
for name, gen in [("A_w_r250", "R250"), ("A_w_mt", "MT19937")]:
    cn, e, C = running(os.path.join(HERE, "runs", name + ".dat"))
    axes[0].plot(cn, e - e0, color=ps.COLORS[gen], label=gen)
    axes[1].plot(cn, C - C0, color=ps.COLORS[gen], label=gen)
nn = np.logspace(np.log10(cn[0]), np.log10(cn[-1]), 200)
sig_e = r_mt["se"] * np.sqrt(r_mt["nmeas"] / nn)
sig_C = r_mt["sC"] * np.sqrt(r_mt["nmeas"] / nn)
for ax, sg in [(axes[0], sig_e), (axes[1], sig_C)]:
    ax.fill_between(nn, -2 * sg, 2 * sg, color=ps.GRID, alpha=0.8, lw=0, label=r"$\pm 2\sigma$ (unbiased)")
    ax.axhline(0, color=ps.EXPECT, ls="--", lw=1.0, label="exact")
    ax.set_xscale("log")
    ax.set_xlabel("number of Wolff cluster updates")
axes[0].set_ylabel(r"$e_{\rm MC} - e_{\rm exact}$")
axes[1].set_ylabel(r"$C_{\rm MC} - C_{\rm exact}$")
axes[0].set_title("Energy per spin, running estimate")
axes[1].set_title("Specific heat, running estimate")
axes[0].set_ylim(-0.006, 0.006); axes[1].set_ylim(-0.1, 0.06)
axes[0].legend(loc="lower right"); axes[1].legend(loc="upper right")
print(ps.save(fig, "ising_hidden_error_growth"))
# where does R250 leave the 2 sigma band for good?
for name in ["A_w_r250"]:
    cn, e, C = running(os.path.join(HERE, "runs", name + ".dat"))
    se = r_mt["se"] * np.sqrt(r_mt["nmeas"] / cn); sC = r_mt["sC"] * np.sqrt(r_mt["nmeas"] / cn)
    out_e = np.abs(e - e0) > 2 * se; out_C = np.abs(C - C0) > 2 * sC
    def last_in(o):
        idx = np.where(~o)[0]
        return cn[idx[-1]] if len(idx) else 0
    print("R250: last batch inside 2 sigma: e at", last_in(out_e), " C at", last_in(out_C))
    for nc in [1e4, 1e5, 1e6, 1e7]:
        i = np.searchsorted(cn, nc)
        i = min(i, len(cn) - 1)
        print(f"  n={cn[i]:.0e}: z_e={(e[i]-e0)/se[i]:+.1f}  z_C={(C[i]-C0)/sC[i]:+.1f}")
