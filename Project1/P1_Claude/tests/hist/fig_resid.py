"""Figure: residuals (measured - predicted D) / predicted SD for every series and depth."""
import sys
sys.path.insert(0, "/Users/noah/Downloads/ComputationalPhysicsPHYS25000/Project1/P1_Claude/tests/hist")
from hist_common import *
import matplotlib.pyplot as plt
import plotstyle as ps
ps.apply()

r = pd.read_csv(os.path.join(H, "theory_residuals.csv"))
order = ["RANDU_d", "RANDM_d", "RAND_d", "MINSTD_d", "R250_d", "RANNUM_d", "MT19937_d", "LOGISTIC_d",
         "RANDU_d_65536", "RANDM_d_65536", "RANDU_s", "MINSTD_s", "R250_s", "RANNUM_s", "MT19937_s"]
label = {s: s.replace("_d_65536", " (seed 65536)").replace("_d", "").replace("_s", " (REAL*4)") for s in order}
gen = {s: s.split("_")[0] for s in order}
fig, ax = plt.subplots(figsize=(10, 5.6))
ax.axvspan(-2, 2, color=ps.GRID, lw=0, label="$\\pm2\\sigma$")
for i, s in enumerate(order):
    t = r[r.series == s]
    y = len(order) - 1 - i
    det = t.sd_pred < 1e-3 * np.maximum(np.abs(t.Dpred - 1), 1e-12)
    ax.scatter(t.resid[~det], np.full((~det).sum(), y), s=26, color=ps.COLORS[gen[s]], zorder=3)
    ax.scatter(t.resid[det], np.full(det.sum(), y), s=26, facecolor="none", edgecolor=ps.COLORS[gen[s]], zorder=3)
    for _, row in t.iterrows():
        if abs(row.resid) > 2.5:
            ax.annotate(f"depth {row.depth}", (row.resid, y), xytext=(6, 4), textcoords="offset points", fontsize=8,
                        color=ps.INK2)
ax.set_yticks(range(len(order)))
ax.set_yticklabels([label[s] for s in order][::-1])
ax.set_xlabel("(measured $D$ $-$ predicted $D$) / predicted standard deviation")
ax.set_title("Theory versus measurement, all 120 depth-scan points")
ax.grid(axis="y", visible=False)
ax.set_xlim(-5, 5)
ax.scatter([], [], s=26, color=ps.MUTED, label="stochastic prediction")
ax.scatter([], [], s=26, facecolor="none", edgecolor=ps.MUTED, label="nearly deterministic histogram")
ax.legend(loc="lower right", fontsize=8)
fig.tight_layout()
print(ps.save(fig, "hist_residuals"))
ok = r[np.isfinite(r.resid)]
print("points", len(ok), "within 2 sigma:", (ok.resid.abs() < 2).sum(), "max", ok.resid.abs().max())
ex = ok[ok.series != "LOGISTIC_d"]
print("excluding LOGISTIC: max", ex.resid.abs().max(), "rms", np.sqrt((ex.resid**2).mean()),
      "chi2", (ex.resid**2).sum(), "n", len(ex))
