"""Remedies and traps: Wolff bias for R250 variants (L=16, K_c), measured
versus the prediction that the bias is linear in the three-point cumulant
kappa_3 of the acceptance indicators at lags (103, 250)."""
import sys, os
import numpy as np
import pandas as pd
sys.path.insert(0, "/Users/noah/Downloads/ComputationalPhysicsPHYS25000/Project1/P1_Claude")
import plotstyle as ps
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
ps.apply()
HERE = "/Users/noah/Downloads/ComputationalPhysicsPHYS25000/Project1/P1_Claude/tests/ising"
df = pd.read_csv(os.path.join(HERE, "results.csv")).set_index("name")

def comb(names, col, scol):
    v = np.array([df.loc[n, col] for n in names]); s = np.array([df.loc[n, scol] for n in names])
    w = 1 / s ** 2
    return (v * w).sum() / w.sum(), 1 / np.sqrt(w.sum())

K3_R250 = 0.0710828
rows = [  # label, runs, kappa3 for prediction (None = no prediction), color
    ("MT19937", ["A_w_mt"], 0.0, ps.COLORS["MT19937"]),
    ("R250 = GFSR(103, 250)", ["A_w_r250"], None, ps.COLORS["R250"]),
    ("R250 decimated by 2", ["C_dec2"], 0.0710588, ps.SINGLE),
    ("R250 decimated by 4", ["C_dec4"], 0.0710828, ps.SINGLE),
    ("R250 XOR R250 (two seeds)", ["C_xor2"], 0.0710850, ps.SINGLE),
    ("time-reversed taps (147, 250)", ["C_rev"], None, ps.SINGLE),
    ("R250 + R250 mod $2^{32}$", ["C_add2"], 0.0263870, ps.SINGLE),
    ("additive LFG(103, 250)", ["C_alfg", "C_alfg_s2"], -0.0147216, ps.SINGLE),
    ("R250 decimated by 3", ["C_dec3"], 0.0, ps.SINGLE),
    ("R250 decimated by 5", ["C_dec5"], 0.0, ps.SINGLE),
    ("R250 XOR MT19937", ["C_gmt"], 0.0, ps.SINGLE),
    ("Ziff four-tap (9689)", ["C_ziff", "C_ziff_s2"], 0.0, ps.SINGLE),
]
dC_R, sC_R = df.loc["A_w_r250", "dC"], df.loc["A_w_r250", "sC"]
de_R, se_R = df.loc["A_w_r250", "de"], df.loc["A_w_r250", "se"]
fig, axes = plt.subplots(1, 2, figsize=(10, 5.6), sharey=True, gridspec_kw=dict(width_ratios=[1, 1]))
y = np.arange(len(rows))[::-1]
table = []
for (lab, names, k3, col), yi in zip(rows, y):
    for ax, c, sc, ref, fac in [(axes[0], "de", "se", de_R, 1e3), (axes[1], "dC", "sC", dC_R, 1.0)]:
        v, s = comb(names, c, sc)
        ax.errorbar(v * fac, yi, xerr=s * fac, fmt="o", ms=7, color=col, capsize=3, lw=1.4)
        if k3 is not None:
            ax.plot(fac * ref * k3 / K3_R250, yi, marker="D", ms=7, mfc="none", mec=ps.EXPECT, mew=1.3, ls="")
    vC, sC = comb(names, "dC", "sC"); ve, se = comb(names, "de", "se")
    predC = dC_R * k3 / K3_R250 if k3 is not None else np.nan
    table.append((lab, ve, se, ve / se, vC, sC, vC / sC, predC, (vC - predC) / np.hypot(sC, abs(predC) * sC_R / abs(dC_R)) if k3 else np.nan))
for ax in axes:
    ax.axvline(0, color=ps.EXPECT, ls="--", lw=1.0)
    ax.axhspan(y[4] - 0.5, y[1] + 0.5, color=ps.GRID, alpha=0.35, lw=0)
axes[0].set_yticks(y); axes[0].set_yticklabels([r[0] for r in rows])
axes[0].set_xlabel(r"$(e_{\rm MC} - e_{\rm exact}) \times 10^{3}$"); axes[1].set_xlabel(r"$C_{\rm MC} - C_{\rm exact}$")
axes[0].set_title("Energy per spin"); axes[1].set_title("Specific heat")
h = [Line2D([], [], marker="o", ls="", color=ps.SINGLE, ms=7, label=r"measured ($\pm 1\sigma$)"),
     Line2D([], [], marker="D", ls="", mfc="none", mec=ps.EXPECT, mew=1.3, ms=7, label=r"prediction $\propto \kappa_3(103,250)$"),
     Line2D([], [], ls="--", color=ps.EXPECT, label="exact")]
axes[1].legend(handles=h, loc="lower left", fontsize=8.5)
print(ps.save(fig, "ising_remedies"))
for t in table:
    print(f"{t[0]:32s} de={t[1]:+.2e}+-{t[2]:.1e} ({t[3]:+.1f}s)  dC={t[4]:+.4f}+-{t[5]:.4f} ({t[6]:+.1f}s)  predC={t[7]:+.4f}  (meas-pred)/sig={t[8]:+.1f}")
