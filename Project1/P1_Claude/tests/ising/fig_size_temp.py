"""R250 Wolff bias of the specific heat versus lattice size L (at K_c) and versus coupling K (L=16)."""
import sys, os
import numpy as np
import pandas as pd
sys.path.insert(0, "/Users/noah/Downloads/ComputationalPhysicsPHYS25000/Project1/P1_Claude")
import plotstyle as ps
import matplotlib.pyplot as plt
ps.apply()
HERE = "/Users/noah/Downloads/ComputationalPhysicsPHYS25000/Project1/P1_Claude/tests/ising"
df = pd.read_csv(os.path.join(HERE, "results.csv")).set_index("name")
def get(names, col="dC", scol="sC", rel=True):
    v = np.array([df.loc[n, col] for n in names]); s = np.array([df.loc[n, scol] for n in names])
    ref = np.array([df.loc[n, "C"] - df.loc[n, "dC"] for n in names]) if col == "dC" else np.array([df.loc[n, "e"] - df.loc[n, "de"] for n in names])
    w = 1 / s ** 2
    m = (v * w).sum() / w.sum(); sm = 1 / np.sqrt(w.sum())
    return (100 * m / ref[0], 100 * sm / ref[0]) if rel else (m, sm)
Ls = [8, 16, 32, 64]
size = {"R250": [["D_L8_r250"], ["A_w_r250"], ["D_L32_r250"], ["D_L64_r250"]],
        "MT19937": [["D_L8_mt"], ["A_w_mt"], ["D_L32_mt", "D_L32_mt_s2"], ["D_L64_mt"]]}
Ks = [0.40, 0.44068679350977151, 0.50]
temp = {"R250": [["H_K040_r250"], ["A_w_r250"], ["H_K050_r250"]], "MT19937": [["H_K040_mt"], ["A_w_mt"], ["H_K050_mt"]]}
fig, axes = plt.subplots(1, 2, figsize=ps.WIDE)
for g, off in [("R250", 1.0), ("MT19937", 1.06)]:
    vals = [get(n) for n in size[g]]
    axes[0].errorbar(np.array(Ls) * off, [v[0] for v in vals], yerr=[v[1] for v in vals], fmt="o-", ms=7,
                     color=ps.COLORS[g], capsize=3, label=g)
    for L, n, v in zip(Ls, size[g], vals):
        print(g, "L=", L, f"dC/C = {v[0]:+.2f} +- {v[1]:.2f} %", " draws/cluster", df.loc[n[0], "draws_per_update"])
    vals = [get(n) for n in temp[g]]
    axes[1].errorbar(np.array(Ks) + (0.002 if g == "MT19937" else 0), [v[0] for v in vals], yerr=[v[1] for v in vals],
                     fmt="o-", ms=7, color=ps.COLORS[g], capsize=3, label=g)
    for K, n, v in zip(Ks, temp[g], vals):
        ev = get(n, "de", "se")
        print(g, "K=", K, f"dC/C = {v[0]:+.2f} +- {v[1]:.2f} %   de/e = {ev[0]:+.3f} +- {ev[1]:.3f} %", " <c>", df.loc[n[0], "csize"], "draws", df.loc[n[0], "draws_per_update"])
for ax in axes:
    ax.axhline(0, color=ps.EXPECT, ls="--", lw=1.0)
    ax.set_ylabel(r"$(C_{\rm MC} - C_{\rm exact}) / C_{\rm exact}$  (%)")
axes[0].legend(loc="lower left"); axes[1].legend(loc="lower right")
axes[0].set_xscale("log", base=2); axes[0].set_xticks(Ls); axes[0].set_xticklabels([str(L) for L in Ls])
axes[0].set_xlabel("lattice size $L$  (coupling $K_c$)")
axes[0].set_title("Wolff specific heat bias versus size")
axes[1].axvline(Ks[1], color=ps.MUTED, ls=":", lw=0.9)
axes[1].text(Ks[1] + 0.002, 0.6, "$K_c$", color=ps.INK2, fontsize=9)
axes[1].set_xlabel("coupling $K = J/k_BT$  ($L = 16$)")
axes[1].set_title("Wolff specific heat bias versus coupling")
print(ps.save(fig, "ising_size_coupling"))
