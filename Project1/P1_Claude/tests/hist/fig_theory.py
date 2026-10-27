"""Figure: measured D versus the finite-population / rounding theory."""
import sys
sys.path.insert(0, "/Users/noah/Downloads/ComputationalPhysicsPHYS25000/Project1/P1_Claude/tests/hist")
from hist_common import *
import matplotlib.pyplot as plt
import plotstyle as ps
from matplotlib.ticker import FuncFormatter, NullFormatter
ps.apply()
PLAIN = FuncFormatter(lambda v, pos: f"{v:g}")

df = load_sum("depth_*_sum.out")
coc, psum = load_period()
fws, fww = load_fw()

dep = np.arange(1, 9)
Bs = 10**dep
Ns = np.where(dep <= 7, 100 * Bs, 10 * Bs)


def predict(model, d):
    B = int(10**d)
    N = int(100 * B if d <= 7 else 10 * B)
    kind = model[0]
    if kind == "periodic":
        t = coc[(coc.tag == model[1]) & (coc.depth == d)]
        return pred_periodic(t.c.values, t.m.values, model[2], N, B)
    if kind == "weights":
        t = fww[fww.depth == d]
        D, sd, pe = pred_weights(t.w.values, t.nb.values, N, B)
        # use the exact sum of squares (the w histogram is binned at 1e-6)
        s2 = fws[fws.depth == d].sdel2.values[0]
        D = 1 + (N / B - 1 / B) * s2 / (B - 1)
        return D, sd, pe
    if kind == "lattice":
        return pred_lattice(model[1], N, B)
    if kind == "logistic":
        return pred_logistic(N, B)
    if kind == "midsq":
        return float(N), 0.0, 1 - 1 / B
    raise ValueError(kind)


PU = ("periodic", "RANDU_d", P_RANDU)
P6 = ("periodic", "RANDU_d_65536", P_65536)
panels = [
    ("RANDU, RANDM (seed 777)", [("RANDU_d", "RANDU", PU), ("RANDM_d", "RANDM", PU)]),
    ("MINSTD", [("MINSTD_d", "MINSTD", ("periodic", "MINSTD_d", P_MINSTD))]),
    ("RAND", [("RAND_d", "RAND", ("periodic", "RAND_d", P_MINSTD))]),
    ("RANDU, RANDM (seed 65536)", [("RANDU_d_65536", "RANDU", P6), ("RANDM_d_65536", "RANDM", P6)]),
    ("LOGISTIC", [("LOGISTIC_d", "LOGISTIC", ("logistic",))]),
    ("R250, RANNUM, MT19937", [("R250_d", "R250", ("lattice", 32)), ("RANNUM_d", "RANNUM", ("lattice", 53)),
                               ("MT19937_d", "MT19937", ("lattice", 53))]),
    ("RANDU stored in REAL*4", [("RANDU_s", "RANDU", ("periodic", "RANDU_s", P_RANDU))]),
    ("MINSTD stored in REAL*4", [("MINSTD_s", "MINSTD", ("periodic", "MINSTD_s", P_MINSTD))]),
    ("R250, RANNUM, MT19937 in REAL*4", [("R250_s", "R250", ("weights",)), ("RANNUM_s", "RANNUM", ("weights",)),
                                         ("MT19937_s", "MT19937", ("weights",))]),
]

rows = []
fig, axs = plt.subplots(3, 3, figsize=(12, 10.5), sharex=True)
for ax, (title, series) in zip(axs.flat, panels):
    pr = np.array([predict(series[0][2], d) for d in dep])
    sdn = null_sd(Bs, Ns)
    ax.fill_between(dep, np.maximum(1 - 3 * sdn, 1e-3), 1 + 3 * sdn, color=ps.GRID, lw=0,
                    label="ideal $\\pm3\\sigma$")
    ax.plot(dep, pr[:, 0], "--", color=ps.EXPECT, lw=1.3, label="prediction", zorder=3)
    mk = ["o", "s", "^"]
    for j, (ser, g, model) in enumerate(series):
        prj = np.array([predict(model, d) for d in dep])
        s = df[df.series == ser].sort_values("depth")
        lab = g if len(series) > 1 else "measured"
        ax.plot(s.depth, s.D, mk[j], color=ps.COLORS[g], ms=6 - j, mfc="none" if j else ps.COLORS[g],
                mew=1.3, label=lab, zorder=4)
        for _, r in s.iterrows():
            Dp, sdp, pe = prj[r.depth - 1]
            sdu = max(sdp, 1e-300)
            rows.append(dict(series=ser, depth=r.depth, B=r.B, N=r.N, D=r.D, Dpred=Dp, sd_pred=sdp,
                             resid=(r.D - Dp) / sdu if sdp > 0 else np.nan,
                             empty=r.nempty / r.B, empty_pred=pe,
                             z_null=r.z))
    ax.set_yscale("log")
    lo_, hi_ = ax.get_ylim()
    ax.set_ylim(max(lo_, 0.5 * np.nanmin(pr[:, 0]) if np.nanmin(pr[:, 0]) < 0.5 else 0.45), hi_)
    lo_, hi_ = ax.get_ylim()
    if hi_ / lo_ < 30:
        cand = [0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 1, 1.5, 2, 3, 5]
        ax.set_yticks([c for c in cand if lo_ <= c <= hi_])
    ax.yaxis.set_major_formatter(PLAIN)
    ax.yaxis.set_minor_formatter(NullFormatter())
    ax.set_title(title)
    ax.legend(loc="best", fontsize=8)
for ax in axs[-1]:
    ax.set_xlabel("depth $\\log_{10} B$")
    ax.set_xticks(dep)
for ax in axs[:, 0]:
    ax.set_ylabel("dispersion index $D$")
fig.tight_layout()
print(ps.save(fig, "hist_theory"))

res = pd.DataFrame(rows)
res.to_csv(os.path.join(H, "theory_residuals.csv"), index=False)
pd.set_option("display.width", 200)
print(res.to_string(float_format=lambda v: f"{v:.6g}"))
