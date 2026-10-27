"""Figure: the student's configuration (B=1e5, N=1e7) repeated over 200 seeds,
and the logistic map's inflated chi-square spread at B=100, N=1e4."""
import sys
sys.path.insert(0, "/Users/noah/Downloads/ComputationalPhysicsPHYS25000/Project1/P1_Claude/tests/hist")
from hist_common import *
import matplotlib.pyplot as plt
import plotstyle as ps
ps.apply()

s5 = pd.read_csv(os.path.join(H, "seeds_e5_sum.out"), sep=r"\s+", names=SUMCOLS)
sl = pd.read_csv(os.path.join(H, "seeds_logis_sum.out"), sep=r"\s+", names=SUMCOLS)
coc, psum = load_period()
B, N = 10**5, 10**7
sd0 = null_sd(B, N)
pred = {}
for g, tag, P in [("RANDU", "RANDU_d", P_RANDU), ("RANDM", "RANDU_d", P_RANDU),
                  ("MINSTD", "MINSTD_d", P_MINSTD), ("RAND", "RAND_d", P_MINSTD)]:
    t = coc[(coc.tag == tag) & (coc.depth == 5)]
    if len(t):
        D, sd, pe = pred_periodic(t.c.values, t.m.values, P, N, B)
        pred[g] = ((D - 1) / sd0, sd / sd0)
for g in ["R250", "RANNUM", "MT19937"]:
    D, sd, pe = pred_lattice(32 if g == "R250" else 53, N, B)
    pred[g] = ((D - 1) / sd0, sd / sd0)

gens = ["RANDU", "RANDM", "RAND", "MINSTD", "R250", "RANNUM", "MT19937"]
fig, axs = plt.subplots(1, 2, figsize=(12, 4.6), gridspec_kw=dict(width_ratios=[1.25, 1]))
ax = axs[0]
ax.axvspan(-1.96, 1.96, color=ps.GRID, lw=0, label="ideal 95% range")
rng = np.random.default_rng(1)
stats_out = []
for i, g in enumerate(gens):
    z = s5[s5.gen == g].z.values
    y = len(gens) - 1 - i
    ax.plot(z, y + rng.uniform(-0.28, 0.28, len(z)), "o", ms=2.6, alpha=0.6, color=ps.COLORS[g], mew=0)
    if g in pred:
        ax.plot([pred[g][0]] * 2, [y - 0.4, y + 0.4], "-", color=ps.INK, lw=1.6,
                label="predicted mean" if i == 0 else None)
    ax.plot([z.mean()] * 2, [y - 0.4, y + 0.4], ":", color=ps.INK2, lw=1.4,
            label="observed mean" if i == 0 else None)
    lo = np.mean(z < -3)
    tt = stats.ttest_1samp(z, pred[g][0])
    ks = stats.kstest(z, stats.norm(pred[g][0], pred[g][1]).cdf)
    stats_out.append((g, len(z), z.mean(), z.std(ddof=1), pred[g][0], pred[g][1], np.mean(z < -1.96),
                      lo, tt.pvalue, ks.pvalue,
                      stats.ttest_1samp(z, 0.0).pvalue))
ax.set_yticks(range(len(gens)))
ax.set_yticklabels(gens[::-1])
ax.set_xlabel("chi-square z-score  ($B=10^5$, $N=10^7$)")
ax.set_title("Student's histogram test, 200 seeds")
ax.legend(loc="lower left", fontsize=8)
ax.grid(axis="y", visible=False)

ax = axs[1]
bins = np.linspace(-5, 9, 57)
xx = np.linspace(-5, 9, 400)
ax.plot(xx, stats.norm.pdf(xx), "--", color=ps.EXPECT, lw=1.3, label="ideal N(0,1)")
for g in ["MT19937", "LOGISTIC"]:
    z = sl[sl.gen == g].z.values
    ax.hist(z, bins=bins, density=True, histtype="step", lw=1.6, color=ps.COLORS[g],
            label=f"{g}: sd(z) = {z.std(ddof=1):.2f}")
ax.set_yscale("log")
ax.set_ylim(5e-4, 0.8)
ax.set_xlabel("chi-square z-score  ($B=100$, $N=10^4$)")
ax.set_ylabel("density")
ax.set_title("Serial correlation widens the null, 1000 seeds")
ax.legend(loc="upper right", fontsize=8)
fig.tight_layout()
print(ps.save(fig, "hist_seed_study"))

print("gen n mean sd pred_mean pred_sd frac(z<-1.96) frac(z<-3) p(mean=pred) KS_p(model) p(mean=0)")
for r in stats_out:
    print(" ".join(f"{v:.4g}" if isinstance(v, float) else str(v) for v in r))
for g in ["LOGISTIC", "MT19937", "RANNUM"]:
    s = sl[sl.gen == g]
    z = s.z.values
    ks = stats.kstest(s.chi2.values, stats.chi2(99).cdf)
    lev = stats.levene(z, sl[sl.gen == "MT19937"].z.values)
    print(g, "meanD=%.4f sdD=%.4f nullsd=%.4f sdz=%.3f frac|z|>3=%.4f KS(chi2_99) p=%.3g  F-test var ratio=%.3f" %
          (s.D.mean(), s.D.std(ddof=1), null_sd(100, 1e4), z.std(ddof=1), np.mean(np.abs(z) > 3), ks.pvalue,
           z.var(ddof=1)))
