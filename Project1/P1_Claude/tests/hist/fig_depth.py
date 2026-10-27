"""Figure: the depth scan (D, chi-square z and empty-bin fraction vs depth)."""
import sys
sys.path.insert(0, "/Users/noah/Downloads/ComputationalPhysicsPHYS25000/Project1/P1_Claude/tests/hist")
from hist_common import *
import matplotlib.pyplot as plt
import plotstyle as ps
ps.apply()

df = load_sum("depth_*_sum.out")
d = df[(df.prec == 8) & (df.seed == 777)]

fig, axs = plt.subplots(1, 3, figsize=(13, 4.4))
dep = np.arange(1, 9)
B = 10.0**dep
N = np.where(dep <= 7, 100 * B, 10 * B)
sd = null_sd(B, N)

ax = axs[0]
ax.fill_between(dep, np.maximum(1 - 3 * sd, 0.02), 1 + 3 * sd, color=ps.GRID, alpha=0.9, lw=0,
                label="ideal, $\\pm 3\\sigma$")
for g in ps.ORDER:
    if g == "MIDSQ":
        continue
    s = d[d.gen == g].sort_values("depth")
    ax.plot(s.depth, s.D, "o-", color=ps.COLORS[g], ms=4, label=g)
ax.text(1.1, 600, "MIDSQ: $D = N$ (off scale)", color=ps.COLORS["MIDSQ"], fontsize=8.5, va="center")
ax.set_yscale("log")
ax.set_xlabel("depth $\\log_{10} B$")
ax.set_ylabel("dispersion index $D = s^2/\\bar n$")
ax.set_title("Dispersion index")
ax.set_ylim(3e-2, 1.5e3)

ax = axs[1]
ax.axhspan(-3, 3, color=ps.GRID, alpha=0.9, lw=0)
for g in ps.ORDER:
    if g == "MIDSQ":
        continue
    s = d[d.gen == g].sort_values("depth")
    ax.plot(s.depth, s.z, "o-", color=ps.COLORS[g], ms=4, label=g)
ax.text(1.1, 3e6, "MIDSQ: $z > 10^3$ (off scale)", color=ps.COLORS["MIDSQ"], fontsize=8.5, va="center")
ax.set_yscale("symlog", linthresh=3, linscale=0.6)
ax.set_ylim(-3e4, 2e7)
ax.set_xlabel("depth $\\log_{10} B$")
ax.set_ylabel("chi-square z-score")
ax.set_title("Significance of the deviation")
ax.set_yticks([-1e4, -1e2, -3, 0, 3, 1e2, 1e4, 1e6])

ax = axs[2]
ax.plot(dep, np.exp(-N / B), "--", color=ps.EXPECT, lw=1.2, label="ideal $e^{-\\lambda}$")
for g in ps.ORDER:
    s = d[d.gen == g].sort_values("depth")
    ax.plot(s.depth, s.nempty / s.B, "o-", color=ps.COLORS[g], ms=4, label=g)
ax.set_xlabel("depth $\\log_{10} B$")
ax.set_ylabel("fraction of empty bins")
ax.set_title("Empty bins")
ax.set_ylim(-0.03, 1.03)

for ax in axs:
    ax.set_xticks(dep)
    ax.axvline(7.5, color=ps.MUTED, lw=0.6, ls=":")
axs[2].text(7.6, 0.3, "$\\lambda=10$", ha="left", fontsize=8, color=ps.INK2)
axs[2].text(7.4, 0.3, "$\\lambda=100$", ha="right", fontsize=8, color=ps.INK2)
h, l = axs[0].get_legend_handles_labels()
h2, l2 = axs[2].get_legend_handles_labels()
h = h[:1] + [h2[1 + ps.ORDER.index(g)] for g in ps.ORDER]
l = l[:1] + ps.ORDER
fig.legend(h, l, loc="lower center", ncol=10, bbox_to_anchor=(0.5, -0.06), frameon=False)
fig.tight_layout(rect=(0, 0.04, 1, 1))
print(ps.save(fig, "hist_depth_scan"))

# table for the report
tab = d.pivot_table(index="gen", columns="depth", values="D")
print(tab.loc[ps.ORDER].to_string(float_format=lambda v: f"{v:.4g}"))
tz = d.pivot_table(index="gen", columns="depth", values="z")
print(tz.loc[ps.ORDER].to_string(float_format=lambda v: f"{v:.3g}"))
te = d.pivot_table(index="gen", columns="depth", values="nempty")
print(te.loc[ps.ORDER].to_string())
