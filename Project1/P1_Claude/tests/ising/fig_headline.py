"""Headline figure: e_MC - e_exact and C_MC - C_exact per generator and algorithm (L=16, K_c).
Metropolis runs with RANDU/RANDM use only independent data (one RNG period per seed),
because the chain becomes exactly periodic after one period (see fig_periodicity.py)."""
import sys, os, glob
import numpy as np
import pandas as pd
sys.path.insert(0, "/Users/noah/Downloads/ComputationalPhysicsPHYS25000/Project1/P1_Claude")
sys.path.insert(0, "/Users/noah/Downloads/ComputationalPhysicsPHYS25000/Project1/P1_Claude/tests/ising")
import plotstyle as ps
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from period_check import stats_first
from analyze import EXACT
ps.apply()
HERE = "/Users/noah/Downloads/ComputationalPhysicsPHYS25000/Project1/P1_Claude/tests/ising"
os.chdir(HERE)
df = pd.read_csv("results.csv").set_index("name")
e0, C0 = EXACT[16]

def wmean(vals):
    v = np.array(vals)
    w1 = 1 / v[:, 1] ** 2; w2 = 1 / v[:, 3] ** 2
    return ((v[:, 0] * w1).sum() / w1.sum(), 1 / np.sqrt(w1.sum()), (v[:, 2] * w2).sum() / w2.sum(), 1 / np.sqrt(w2.sum()))

def from_df(n):
    r = df.loc[n]; return (r["de"], r["se"], r["dC"], r["sC"])

def universes(base, pattern, nb=128):
    out = []
    if base:
        e, se, _, C, sC, _ = stats_first(base, nb); out.append((e - e0, se, C - C0, sC))
    for f in sorted(glob.glob(pattern)):
        e, se, _, C, sC, _ = stats_first(f, nb); out.append((e - e0, se, C - C0, sC))
    return out

gens = ["MT19937", "RANNUM", "R250", "MINSTD", "RAND", "RANDU", "RANDM"]
wolff = {"MT19937": from_df("A_w_mt"), "RANNUM": from_df("A_w_rannum"), "R250": from_df("A_w_r250"),
         "MINSTD": from_df("A_w_minstd"), "RAND": from_df("A_w_rand"), "RANDU": from_df("A_w_randu"),
         "RANDM": from_df("A_w_randm")}
metro = {"MT19937": from_df("B_m_mt"), "RANNUM": from_df("B_m_rannum"), "R250": from_df("B_m_r250"),
         "MINSTD": wmean([from_df("B_m_minstd")] + universes(None, "runs/M_g4_s*.dat")),
         "RAND": from_df("B_m_rand"),
         "RANDU": wmean(universes("runs/B_m_randu.dat", "runs/M_g1_s*.dat")),
         "RANDM": wmean(universes("runs/B_m_randm.dat", "runs/M_g2_s*.dat"))}
rows = []
for g in gens:
    for alg, d in [("Wolff", wolff), ("Metropolis", metro)]:
        de, se, dC, sC = d[g]
        rows.append((g, alg, de, se, de / se, dC, sC, dC / sC))
tab = pd.DataFrame(rows, columns=["gen", "alg", "de", "se", "zE", "dC", "sC", "zC"])
tab.to_csv("headline_table.csv", index=False)
print(tab.to_string(float_format=lambda x: f"{x:.3g}"))

fig, axes = plt.subplots(1, 2, figsize=(10, 4.8), sharey=True)
ypos = np.arange(len(gens))[::-1]
for ax, col, scol, fac, lab in [(axes[0], "de", "se", 1e3, r"$(e_{\rm MC} - e_{\rm exact}) \times 10^{3}$"),
                                (axes[1], "dC", "sC", 1.0, r"$C_{\rm MC} - C_{\rm exact}$")]:
    for g, y in zip(gens, ypos):
        for alg, off, mk, filled in [("Wolff", +0.17, "o", True), ("Metropolis", -0.17, "s", False)]:
            r = tab[(tab.gen == g) & (tab.alg == alg)].iloc[0]
            ax.errorbar(r[col] * fac, y + off, xerr=r[scol] * fac, fmt=mk, ms=7, color=ps.COLORS[g],
                        mfc=ps.COLORS[g] if filled else "white", mew=1.6, capsize=3, lw=1.4)
            z = r[col] / r[scol]
            if abs(z) > 4:
                ax.annotate(f"{z:+.0f}$\\sigma$", (r[col] * fac, y + off), xytext=(0, 7), textcoords="offset points",
                            ha="center", fontsize=8.5, color=ps.INK2)
    ax.axvline(0, color=ps.EXPECT, ls="--", lw=1.0)
    ax.set_xlabel(lab)
axes[0].set_yticks(ypos); axes[0].set_yticklabels(gens)
axes[0].set_title("Energy per spin"); axes[1].set_title("Specific heat")
h = [Line2D([], [], marker="o", ls="", color=ps.INK2, ms=7, label="Wolff cluster"),
     Line2D([], [], marker="s", ls="", color=ps.INK2, mfc="white", mew=1.6, ms=7, label="Metropolis sweep"),
     Line2D([], [], ls="--", color=ps.EXPECT, label="exact (Kaufman)")]
axes[1].legend(handles=h, loc="lower left")
print(ps.save(fig, "ising_headline"))
