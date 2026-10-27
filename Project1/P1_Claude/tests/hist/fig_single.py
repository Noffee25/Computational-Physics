"""Figure: anatomy of the REAL*4 storage artifact at B = 1e7, N = 1e9."""
import sys
sys.path.insert(0, "/Users/noah/Downloads/ComputationalPhysicsPHYS25000/Project1/P1_Claude/tests/hist")
from hist_common import *
import matplotlib.pyplot as plt
import plotstyle as ps
ps.apply()

B, N = 10**7, 10**9
lam = N / B
prof = pd.concat([pd.read_csv(f, sep=r"\s+", names=["tag", "ir", "xlo", "xhi", "mean", "D"])
                  for f in sorted(glob.glob(os.path.join(H, "depth_*_pro.out")))], ignore_index=True)
cc = load_coc("depth_*_coc.out")
coc, psum = load_period()
fws, fww = load_fw()
df = load_sum("depth_*_sum.out")


def octave_D(x):
    """Continuous-uniform input rounded to REAL*4: within the octave
    [2^-k-1, 2^-k) every bin holds floor(g) or ceil(g) floats, g = h/u."""
    k = np.floor(-np.log2(x))
    u = 2.0 ** (-k - 24)
    g = (1.0 / B) / u
    fr = g - np.floor(g)
    return 1 + lam * fr * (1 - fr) / g**2


fig = plt.figure(figsize=(12, 8.2))
gs = fig.add_gridspec(2, 3, height_ratios=[1.0, 1.0])
ax = fig.add_subplot(gs[0, :])
xx = np.logspace(-3, 0, 3000)[:-1]
for ser, g, mk in [("MT19937_s", "MT19937", "o"), ("MINSTD_s", "MINSTD", "s"), ("RANDU_s", "RANDU", "^")]:
    p = prof[prof.tag == f"{ser}_e7"]
    xm = 0.5 * (p.xlo + p.xhi)
    ax.plot(xm, p.D, mk, ms=2.6, color=ps.COLORS[g], label=f"{g}, REAL*4", alpha=0.85, mew=0)
ax.plot(xx, octave_D(xx), "--", color=ps.EXPECT, lw=1.2, label="rounding model (continuous input)")
for e in [0.5, 0.25, 0.125, 0.0625]:
    ax.axvline(e, color=ps.MUTED, lw=0.6, ls=":")
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlim(1.5e-3, 1)
ax.set_xlabel("$x$ (centre of a block of $10^4$ consecutive bins)")
ax.set_ylabel("dispersion index in block")
ax.set_title("Where the REAL*4 artifact lives ($B=10^7$, $N=10^9$)")
ax.legend(loc="upper left", fontsize=8, ncol=2)
ax.set_xticks([0.002, 0.005, 0.01, 0.02, 0.0625, 0.125, 0.25, 0.5, 1])
ax.set_xticklabels(["0.002", "0.005", "0.01", "0.02", "1/16", "1/8", "1/4", "1/2", "1"])

kmax = 220
k = np.arange(kmax + 1)
for j, (ser, g) in enumerate([("MT19937_s", "MT19937"), ("MINSTD_s", "MINSTD"), ("RANDU_s", "RANDU")]):
    ax = fig.add_subplot(gs[1, j])
    t = cc[(cc.tag == f"{ser}_e7") & (cc.k >= 0)]
    frac = np.zeros(kmax + 1)
    kk = t.k.values.astype(int); ok = kk <= kmax
    frac[kk[ok]] = t.nb.values[ok] / B
    ax.bar(k, frac, width=1.0, color=ps.COLORS[g], lw=0, label="measured")
    ax.plot(k, stats.poisson.pmf(k, lam), "--", color=ps.EXPECT, lw=1.1, label="Poisson(100)")
    if ser == "MT19937_s":
        t7 = fww[fww.depth == 7]
        kp, pp = coc_weights(t7.w.values, t7.nb.values, N, B, kmax)
        lab = "rounding model"
    else:
        P = P_RANDU if g == "RANDU" else P_MINSTD
        tt = coc[(coc.tag == ser) & (coc.depth == 7)]
        kp, pp = coc_periodic(tt.c.values, tt.m.values, P, N, B, kmax)
        lab = "full-period model"
    ax.plot(kp, pp, "-", color=ps.INK, lw=0.9, label=lab)
    r = df[df.tag == f"{ser}_e7"].iloc[0]
    ax.set_title(f"{g} in REAL*4:  $D$ = {r.D:.3g}")
    ax.set_yscale("log")
    ax.set_ylim(1e-7 / 2, 0.1)
    ax.set_xlim(20, kmax)
    ax.set_xlabel("hits in a bin, $k$")
    if j == 0:
        ax.set_ylabel("fraction of bins")
    ax.legend(loc="upper right", fontsize=7.5)
fig.tight_layout()
print(ps.save(fig, "hist_single_precision"))

# numbers: octave D values and measured block D averages per octave
for kk_ in range(0, 6):
    lo, hi = 2.0 ** (-kk_ - 1), 2.0 ** (-kk_)
    row = [f"[{lo:.4g},{hi:.4g}) model D={octave_D(np.array([0.75 * hi]))[0]:.4f}"]
    for ser in ["MT19937_s", "MINSTD_s", "RANDU_s", "RANNUM_s", "R250_s"]:
        p = prof[prof.tag == f"{ser}_e7"]
        sel = (p.xlo >= lo) & (p.xhi <= hi)
        row.append(f"{ser}={p[sel].D.mean():.4f}")
    print("  ".join(row))
