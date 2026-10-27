"""Figure: count-of-counts at B = 1e7, N = 1e9 for all nine generators."""
import sys
sys.path.insert(0, "/Users/noah/Downloads/ComputationalPhysicsPHYS25000/Project1/P1_Claude/tests/hist")
from hist_common import *
import matplotlib.pyplot as plt
import plotstyle as ps
ps.apply()

df = load_sum("depth_*_sum.out")
cc = load_coc("depth_*_coc.out")
coc, psum = load_period()
B, N = 10**7, 10**9
lam = N / B

model = {
    "RANDU": ("periodic", "RANDU_d", P_RANDU),
    "RANDM": ("periodic", "RANDU_d", P_RANDU),
    "RAND": ("periodic", "RAND_d", P_MINSTD),
    "MINSTD": ("periodic", "MINSTD_d", P_MINSTD),
    "LOGISTIC": ("logistic",),
}
xr = {"RAND": (-5, 300), "LOGISTIC": (-10, 1450), "MIDSQ": (-0.5, 3.5)}

fig, axs = plt.subplots(3, 3, figsize=(12, 9.6))
for ax, g in zip(axs.flat, ps.ORDER):
    t = cc[cc.tag == f"{g}_d_e7"]
    t = t[t.k >= 0].sort_values("k")
    over = cc[(cc.tag == f"{g}_d_e7") & (cc.k < 0)].nb.sum()
    lo, hi = xr.get(g, (35, 175))
    kmax = int(max(hi, t.k.max() if len(t) else 0)) + 1
    k = np.arange(kmax + 1)
    frac = np.zeros(kmax + 1)
    kk = t.k.values.astype(int)
    ok = kk <= kmax
    frac[kk[ok]] = t.nb.values[ok] / B
    if g == "LOGISTIC":
        nzk = frac > 0
        ax.plot(k[nzk], frac[nzk], "o", ms=2.2, color=ps.COLORS[g], mew=0, label=f"{g} (measured)")
    else:
        ax.bar(k, frac, width=1.0, color=ps.COLORS[g], lw=0, label=f"{g} (measured)")
    pk = stats.poisson.pmf(k, lam)
    ax.plot(k, pk, "--", color=ps.EXPECT, lw=1.2, label="Poisson(100)")
    if g in model:
        m = model[g]
        if m[0] == "periodic":
            tt = coc[(coc.tag == m[1]) & (coc.depth == 7)]
            kp, pp = coc_periodic(tt.c.values, tt.m.values, m[2], N, B, kmax)
        else:
            kp, pp = coc_logistic(N, B, kmax)
        ax.plot(kp, pp, "-", color=ps.INK, lw=0.9 if g != "LOGISTIC" else 0.6, label="finite-period model")
    ax.set_yscale("log")
    ax.set_ylim(1e-7 / 2, 1.5)
    ax.set_xlim(lo, hi)
    r = df[df.tag == f"{g}_d_e7"].iloc[0]
    ax.set_title(f"{g}:  $D$ = {r.D:.3g}")
    if g == "MIDSQ":
        ax.text(0.97, 0.45, "one bin holds\n999,999,999 hits\n(state 0 from draw 2)", transform=ax.transAxes,
                ha="right", va="center", fontsize=8.5, color=ps.INK2)
        ax.set_xticks([0, 1, 2, 3])
    ax.legend(loc="upper right" if g not in ("MIDSQ",) else "upper center", fontsize=7.5)
for ax in axs[-1]:
    ax.set_xlabel("hits in a bin, $k$")
for ax in axs[:, 0]:
    ax.set_ylabel("fraction of bins")
fig.tight_layout()
print(ps.save(fig, "hist_count_of_counts"))

# numbers for the text
for g in ps.ORDER:
    t = cc[(cc.tag == f"{g}_d_e7") & (cc.k >= 0)]
    r = df[df.tag == f"{g}_d_e7"].iloc[0]
    k = t.k.values.astype(float); w = t.nb.values.astype(float)
    nz = k > 0
    if nz.sum() > 0:
        mfill = np.sum(k[nz] * w[nz]) / np.sum(w[nz])
        sfill = np.sqrt(np.sum(w[nz] * (k[nz] - mfill) ** 2) / (np.sum(w[nz]) - 1))
    else:
        mfill = sfill = np.nan
    print(f"{g:9s} D={r.D:.6g} empty={r.nempty} min={r.minc} max={r.maxc} n4lo={r.n4lo} n4hi={r.n4hi} "
          f"e4lo={r.e4lo:.1f} e4hi={r.e4hi:.1f} filled-mean={mfill:.2f} filled-sd={sfill:.3f}")
