"""Lag scan: bias of the Wolff specific heat and energy for two-tap XOR
shift registers x_n = x_{n-p} XOR x_{n-q} with increasing register length q."""
import sys, os
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
sys.path.insert(0, "/Users/noah/Downloads/ComputationalPhysicsPHYS25000/Project1/P1_Claude")
import plotstyle as ps
import matplotlib.pyplot as plt
ps.apply()
HERE = "/Users/noah/Downloads/ComputationalPhysicsPHYS25000/Project1/P1_Claude/tests/ising"
df = pd.read_csv(os.path.join(HERE, "results.csv")).set_index("name")
runs = [("C_q89", 38, 89), ("C_q127", 63, 127), ("A_w_r250", 103, 250), ("C_q521", 168, 521),
        ("C_q607", 273, 607), ("C_q1279", 418, 1279), ("C_q2281", 1029, 2281), ("C_q9689", 4187, 9689)]
q = np.array([r[2] for r in runs], float)
dC = np.array([df.loc[r[0], "dC"] for r in runs]); sC = np.array([df.loc[r[0], "sC"] for r in runs])
de = np.array([df.loc[r[0], "de"] for r in runs]); se = np.array([df.loc[r[0], "se"] for r in runs])
dpc = df.loc["A_w_mt", "draws_per_update"]
f = lambda x, A, q0: A * np.exp(-x / q0)
sel = q >= 250
pC, cC = curve_fit(f, q[sel], dC[sel], p0=[-0.1, 300], sigma=sC[sel], absolute_sigma=True)
chiC = (((dC[sel] - f(q[sel], *pC)) / sC[sel]) ** 2).sum()
pe, ce = curve_fit(f, q[sel], de[sel], p0=[0.005, 300], sigma=se[sel], absolute_sigma=True)
chie = (((de[sel] - f(q[sel], *pe)) / se[sel]) ** 2).sum()
print(f"C fit (q>=250): A={pC[0]:.4f}+-{np.sqrt(cC[0,0]):.4f} q0={pC[1]:.1f}+-{np.sqrt(cC[1,1]):.1f} chi2={chiC:.2f}/{sel.sum()-2}")
print(f"e fit (q>=250): A={pe[0]:.5f}+-{np.sqrt(ce[0,0]):.5f} q0={pe[1]:.1f}+-{np.sqrt(ce[1,1]):.1f} chi2={chie:.2f}/{sel.sum()-2}")
print("draws per cluster", dpc, " q0/dpc (C) =", pC[1] / dpc)
for (nm, p, qq), a, b, c, d in zip(runs, de, se, dC, sC):
    print(f"  GFSR({p},{qq}): de={a:+.2e}+-{b:.1e} ({a/b:+.1f} sigma)  dC={c:+.4f}+-{d:.4f} ({c/d:+.1f} sigma)  pred_C={f(qq,*pC):+.4f}")
fig, axes = plt.subplots(1, 2, figsize=ps.WIDE)
xx = np.logspace(np.log10(80), np.log10(11000), 300)
for ax, y, s, pp, lab, nm in [(axes[0], de, se, pe, r"$e_{\rm MC} - e_{\rm exact}$", "Energy per spin"),
                              (axes[1], dC, sC, pC, r"$C_{\rm MC} - C_{\rm exact}$", "Specific heat")]:
    cols = [ps.COLORS["R250"] if qq == 250 else ps.SINGLE for qq in q]
    for xi, yi, si, ci in zip(q, y, s, cols):
        ax.errorbar(xi, yi, yerr=si, fmt="o", ms=7, color=ci, capsize=3, lw=1.4)
    ax.plot(xx, f(xx, *pp), color=ps.EXPECT, ls="--", lw=1.0,
            label=f"fit $A\\,e^{{-q/q_0}}$ ($q\\geq 250$), $q_0$ = {pp[1]:.0f} draws")
    ax.axhline(0, color=ps.MUTED, lw=0.8)
    ax.set_xscale("log")
    ax.set_xlabel(r"register length $q$ (lag of $x_{n-q}$)")
    ax.set_ylabel(lab); ax.set_title(nm)
    # mark the largest possible draw count of one cluster on 16x16 (2N+1)
    ax.axvline(2 * 256 + 1, color=ps.MUTED, ls=":", lw=0.9)
    ax.plot([], [], "o", color=ps.SINGLE, label=r"two-tap XOR GFSR$(p,q)$")
    ax.plot([], [], "o", color=ps.COLORS["R250"], label="R250 = GFSR(103, 250)")
axes[0].legend(loc="upper right", fontsize=8.5)
axes[1].legend(loc="lower right", fontsize=8.5)
axes[1].text(2 * 256 + 40, -0.045, "largest possible\ndraw count of\none cluster (2N+1)", fontsize=8, color=ps.INK2, va="center")
print(ps.save(fig, "ising_lag_scan"))
