"""Report figures for the 'spacing' family. Run with /usr/local/bin/python3 -I."""
import sys, os
sys.path.insert(0, "/Users/noah/Downloads/ComputationalPhysicsPHYS25000/Project1/P1_Claude")
sys.path.insert(0, "/Users/noah/Downloads/ComputationalPhysicsPHYS25000/Project1/P1_Claude/tests/spacing")
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import plotstyle as ps
from ms_stats import *
from analyze_scan import load
from blocks import CONFIGS

ps.apply()
D = "/Users/noah/Downloads/ComputationalPhysicsPHYS25000/Project1/P1_Claude/tests/spacing/"
B = pd.read_csv(D + "blocks.csv")
LABEL = {t: l for t, l, c in CONFIGS}
CKEY = {t: c for t, l, c in CONFIGS}
EULER = 0.5772156649


def exact_density_T(n, Tgrid):
    """density of T = (n+1)M - ln(n+1) from the exact CDF (finite differences)."""
    x = (Tgrid + np.log(n + 1)) / (n + 1)
    lF, lS = maxsp_cdf_sf(np.full_like(x, n), x)
    F = np.exp(lF)
    return np.gradient(F, Tgrid)


def fmt_p(lp):
    if not np.isfinite(lp) or lp < -300:
        return r"$< 10^{-300}$"
    if lp > -2:
        return "= %.2f" % (10 ** lp)
    return r"$= 10^{%.0f}$" % lp


# ---------------------------------------------------------------- Fig 1
def fig_dist():
    cases = [("MT", 1000), ("RANNUM", 1000000), ("LOGISTIC", 1000), ("RANDU65536", 1000)]
    fig, axes = plt.subplots(2, 2, figsize=(10, 6.6), sharex=True)
    for ax, (tag, N) in zip(axes.flat, cases):
        df = load(tag)
        b = df[(df.N == N) & (df.flag == 0)]
        n = b.n.values.astype(float)
        T = gumbel_T(n, b.M.values)
        row = B[(B.tag == tag) & (B.N == N)].iloc[0]
        bins = np.linspace(-3, 9, 61)
        col = ps.COLORS[CKEY[tag]]
        ax.hist(np.clip(T, -3, 9), bins=bins, density=True, histtype="step", color=col, lw=1.6,
                label="%s, sd(T) = %.2f" % (LABEL[tag], T.std(ddof=1)))
        tg = np.linspace(-3, 9, 600)
        ax.plot(tg, exact_density_T(float(N), tg), "--", color=ps.EXPECT, lw=1.3,
                label="exact null law")
        ax.set_title(r"%s, $N=10^{%d}$, $R=%d$" % (LABEL[tag], round(np.log10(N)), len(b)))
        txt = "KS p %s\nAD p %s" % (fmt_p(row.max_log10p_ks), fmt_p(row.max_log10p_ad))
        ax.text(0.97, 0.55, txt, transform=ax.transAxes, ha="right", va="top", fontsize=9, color=ps.INK2)
        ax.legend(loc="upper right")
        ax.set_ylim(0, 0.52)
    for ax in axes[1]:
        ax.set_xlabel(r"$T = (n+1)M_n - \ln(n+1)$")
    for ax in axes[:, 0]:
        ax.set_ylabel("density")
    ps.save(fig, "spacing_maxgap_dist")


# ---------------------------------------------------------------- Fig 2
def fig_scan():
    tags = ["MT", "RANNUM", "R250", "MINSTD", "RANDU", "RANDM", "RAND", "MTsingle",
            "RANDUsingle", "LOGISTIC", "RANDU65536", "MIDSQ12345678"]
    period = {"MINSTD": 2 ** 31 - 2, "RANDU": 2 ** 29, "RANDM": 2 ** 29, "RAND": 2 ** 31 - 2,
              "RANDUsingle": 2 ** 29, "LOGISTIC": 14324259 + 5638349, "RANDU65536": 8192,
              "MIDSQ12345678": 10352 + 100}
    grid = {"RAND": 2 ** 22, "MTsingle": 2 ** 24, "RANDUsingle": 2 ** 24}
    fig, axes = plt.subplots(3, 4, figsize=(11, 8.2), sharex=True, sharey=True)
    lo, hi = stats_gumbel_band()
    rng = np.random.default_rng(3)
    for ax, tag in zip(axes.flat, tags):
        df = load(tag)
        df = df[df.flag == 0]
        col = ps.COLORS[CKEY[tag]]
        Ns, means = [], []
        for N, b in df.groupby("N"):
            T = gumbel_T(b.n.values.astype(float), b.M.values)
            if len(T) > 150:
                T_show = rng.choice(T, 150, replace=False)
            else:
                T_show = T
            xj = N * np.exp(rng.normal(0, 0.06, len(T_show)))
            ax.scatter(xj, T_show, s=4, color=col, alpha=0.35, lw=0)
            Ns.append(N); means.append(T.mean())
        ax.plot(Ns, means, "-o", color=col, ms=3.5, lw=1.3, label="block mean")
        xs = np.array([5e2, 4e9])
        ax.fill_between(xs, lo, hi, color=ps.GRID, alpha=0.7, lw=0, label="null 99.8% band")
        ax.axhline(EULER, color=ps.EXPECT, ls="--", lw=1.0)
        if tag in period:
            ax.axvline(period[tag], color=ps.INK2, ls=":", lw=1.2, label="period")
        if tag in grid:
            ax.axvline(grid[tag], color=ps.INK2, ls="-.", lw=1.0, label="grid size $G$")
        ax.set_xscale("log")
        ax.set_yscale("symlog", linthresh=8, linscale=1.6)
        ax.set_title(LABEL[tag])
        ax.set_xlim(5e2, 4e9)
        ax.set_ylim(-20, 2e4)
        ax.set_yticks([-10, -5, 0, 5, 10, 100, 1000, 10000])
        ax.set_yticklabels(["-10", "-5", "0", "5", "10", "100", "1000", r"$10^4$"])
        ax.yaxis.set_minor_locator(matplotlib.ticker.NullLocator())
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    h = [Line2D([], [], marker="o", ls="", color=ps.INK2, alpha=0.5, ms=3, label="replicate (subsample)"),
         Line2D([], [], marker="o", color=ps.INK2, ms=4, label="block mean"),
         Patch(color=ps.GRID, label="Gumbel 0.1% to 99.9%"),
         Line2D([], [], ls="--", color=ps.EXPECT, label=r"Gumbel mean $\gamma$"),
         Line2D([], [], ls=":", color=ps.INK2, label="period or $\\mu+\\lambda$"),
         Line2D([], [], ls="-.", color=ps.INK2, label="number of grid values $G$")]
    fig.legend(handles=h, loc="lower center", ncol=6, bbox_to_anchor=(0.5, -0.035), fontsize=9)
    for ax in axes[-1]:
        ax.set_xlabel("draws per replicate $N$")
    for ax in axes[:, 0]:
        ax.set_ylabel(r"$T=(n+1)M_n-\ln(n+1)$")
    ps.save(fig, "spacing_maxgap_scan")


def stats_gumbel_band():
    # 0.1% and 99.9% quantiles of the standard Gumbel law
    q = np.array([0.001, 0.999])
    v = -np.log(-np.log(q))
    return v[0], v[1]


# ---------------------------------------------------------------- Fig 3
def fig_breakmap():
    tags = ["MT", "RANNUM", "R250", "MINSTD", "RANDU", "RANDM", "RAND", "MTsingle", "RANDUsingle",
            "LOGISTIC", "RANDU65536", "RANDM65536", "MIDSQ12345678", "MIDSQ777"]
    Nlist = sorted(B.N.unique())
    fig, axes = plt.subplots(1, 2, figsize=(12.5, 6.2), gridspec_kw=dict(width_ratios=[13, 8]))
    cmap = plt.get_cmap("magma_r")
    for ax, (stat, title, Nmax) in zip(axes, [("max_log10p", "Maximum spacing", 3e9),
                                               ("min_log10p", "Minimum spacing", 1.1e7)]):
        cols = [N for N in Nlist if N <= Nmax]
        M = np.full((len(tags), len(cols)), np.nan)
        for i, t in enumerate(tags):
            for k, N in enumerate(cols):
                r = B[(B.tag == t) & (B.N == N)]
                if len(r):
                    v = r[stat].iloc[0]
                    M[i, k] = -max(v, -300.0) if np.isfinite(v) else 300.0
        disp = np.clip(M, 0, 20)
        im = ax.imshow(disp, cmap=cmap, vmin=0, vmax=20, aspect="auto")
        for i in range(len(tags)):
            for k in range(len(cols)):
                if np.isnan(M[i, k]):
                    ax.text(k, i, "", ha="center", va="center")
                    continue
                v = M[i, k]
                s = ">300" if v >= 300 else ("%.0f" % v if v >= 2 else "")
                ax.text(k, i, s, ha="center", va="center", fontsize=7,
                        color="white" if v > 11 else ps.INK)
        ax.set_xticks(range(len(cols)))
        ax.set_xticklabels([n_label(N) for N in cols], rotation=60, fontsize=8)
        ax.set_yticks(range(len(tags)))
        ax.set_yticklabels([LABEL[t] for t in tags] if ax is axes[0] else [], fontsize=9)
        ax.set_title(title)
        ax.set_xlabel("draws per replicate $N$")
        ax.grid(False)
        for sp in ax.spines.values():
            sp.set_visible(False)
    cb = fig.colorbar(im, ax=axes, shrink=0.8, pad=0.02)
    cb.set_label(r"$-\log_{10} p$ (block, clipped at 20)")
    os.makedirs(ps.FIGDIR, exist_ok=True)
    fig.savefig(os.path.join(ps.FIGDIR, "spacing_breakmap.png"), bbox_inches="tight")
    plt.close(fig)


def n_label(N):
    special = {536870912: r"$2^{29}$", 2147483646: r"$2^{31}\!-\!2$"}
    if N in special:
        return special[N]
    e = np.log10(N)
    if abs(e - round(e)) < 1e-9:
        return r"$10^{%d}$" % round(e)
    m = N / 10 ** np.floor(e)
    return r"$%g\cdot10^{%d}$" % (m, np.floor(e))


# ---------------------------------------------------------------- Fig 4
def fig_minsp():
    fig, axes = plt.subplots(1, 2, figsize=ps.WIDE)
    ax = axes[0]
    for tag in ["MT", "RANNUM", "LOGISTIC", "RANDU", "MINSTD"]:
        g = B[(B.tag == tag) & (B.j == 0)]
        if tag == "LOGISTIC":
            g = g[g.N <= 3e5]
        ax.plot(g.N, g.mean_nnm, "o-", color=ps.COLORS[CKEY[tag]], ms=4, label=LABEL[tag])
    Ng = np.logspace(3, 6.1, 50)
    ax.axhline(1.0, color=ps.EXPECT, ls="--", lw=1.1, label="iid: 1")
    ax.axhline(2.0, color=ps.EXPECT, ls=":", lw=1.3, label="tent-map clusters: 2")
    ax.plot(Ng, Ng ** 2 * 2.0 ** -30, "-.", color=ps.COLORS["RANDU"], lw=1.0, label=r"floor $n^2\,2^{-30}$ (RANDU)")
    ax.plot(Ng, Ng ** 2 / (2.0 ** 31 - 1), "-.", color=ps.COLORS["MINSTD"], lw=1.0, label=r"floor $n^2/(2^{31}-1)$ (MINSTD)")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("draws per replicate $N$")
    ax.set_ylabel(r"mean of $n(n+1)\,m_n$")
    ax.set_title("Minimum spacing: too regular, too clustered")
    ax.legend(fontsize=8, loc="upper left")
    ax.set_ylim(0.5, 3e3)

    ax = axes[1]
    theo = {"RAND": 2.0 ** 22, "MTsingle": 0.75 * 2 ** 25, "RANDUsingle": 2.638e7, "R250": 2.0 ** 32}
    for tag in ["RAND", "MTsingle", "RANDUsingle", "R250", "MT"]:
        g = B[(B.tag == tag) & (B.j == 0)]
        col = ps.COLORS[CKEY[tag]]
        mk = "s" if "single" in tag else "o"
        ax.plot(g.N, g.frac_m0, mk, color=col, ms=4.5, label=LABEL[tag],
                mfc="white" if "single" in tag else col)
        if tag in theo:
            Ng = np.logspace(2.8, 6.2, 200)
            ax.plot(Ng, -np.expm1(-Ng * (Ng - 1) / (2 * theo[tag])), "--", color=col, lw=1.0)
    ax.set_xscale("log")
    ax.set_xlabel("draws per replicate $N$")
    ax.set_ylabel("fraction of replicates with a tie ($m_n=0$)")
    ax.set_title("Ties reveal the output grid")
    ax.plot([], [], "--", color=ps.EXPECT, lw=1.0, label=r"$1-e^{-N^2/2G}$, $G$ from the grid")
    ax.legend(fontsize=8, loc="center right")
    ps.save(fig, "spacing_minspacing")


# ---------------------------------------------------------------- Fig 5
def fig_dups():
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.6), gridspec_kw=dict(width_ratios=[1.25, 1]))
    ax = axes[0]
    Ng = np.logspace(3, 9.1, 300)
    curves = [("RAND", "dups/RAND.txt", 2.0 ** 22), ("MTsingle", "dups/MTsingle.txt", 0.75 * 2 ** 25),
              ("RANDU", "dups/RANDU.txt", None)]
    for tag, f, G in curves:
        a = np.loadtxt(D + f, comments="#")
        col = ps.COLORS[CKEY[tag]]
        ok = a[:, 1] > 0
        ax.plot(a[ok, 0], a[ok, 1], "-", color=col, lw=1.8, label=LABEL[tag] + " (one stream)")
    for tag, G in [("R250", 2.0 ** 32)]:
        g = B[(B.tag == tag) & (B.j == 0) & (B.dup_mean > 0)]
        ax.plot(g.N, g.dup_mean, "o", color=ps.COLORS[CKEY[tag]], ms=4.5, label=LABEL[tag] + " (replicate mean)")
    for G, lab in [(2.0 ** 22, r"$G=2^{22}$"), (0.75 * 2 ** 25, r"$G=0.75\cdot 2^{25}$"),
                   (2.0 ** 32, r"$G=2^{32}$"), (2.0 ** 53, r"$G=2^{53}$")]:
        ax.plot(Ng, expected_dups(Ng, G), "--", color=ps.EXPECT, lw=0.9)
        if G == 2.0 ** 22:
            ax.text(np.sqrt(2 * G * 300) / 1.4, 300, lab, fontsize=8, color=ps.INK2, ha="right", va="center")
        elif G == 0.75 * 2 ** 25:
            ax.text(np.sqrt(2 * G * 1e4) * 1.4, 1e4, lab, fontsize=8, color=ps.INK2, ha="left", va="center")
        else:
            ax.text(np.sqrt(2 * G * 3) * 1.4, 3, lab, fontsize=8, color=ps.INK2, ha="left", va="center")
    ax.plot(Ng[Ng > 2 ** 29], Ng[Ng > 2 ** 29] - 2 ** 29, ":", color=ps.COLORS["RANDU"], lw=1.4,
            label=r"$N-2^{29}$")
    ax.axvline(2 ** 29, color=ps.COLORS["RANDU"], lw=0.8, ls=":")
    ax.text(2 ** 29 * 0.85, 3e3, "RANDU period\n(no repeats before)", ha="right", fontsize=8, color=ps.COLORS["RANDU"])
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlim(1e3, 1.3e9); ax.set_ylim(0.5, 2e9)
    ax.set_xlabel("number of draws $N$")
    ax.set_ylabel(r"repeated values $D(N)$")
    ax.set_title("Birthday collisions")
    ax.legend(fontsize=7.8, loc="upper left")

    # G_eff dot plot
    ax = axes[1]
    rows = geff_rows()
    y = np.arange(len(rows))[::-1]
    for yi, r in zip(y, rows):
        col = ps.COLORS[r["ckey"]]
        if r["kind"] == "measured":
            ax.errorbar(r["log2"], yi, xerr=[[r["log2"] - r["lo"]], [r["hi"] - r["log2"]]], fmt="o",
                        color=col, ms=5, capsize=3, zorder=4)
        else:
            ax.annotate("", xy=(r["log2"] + 3.2, yi), xytext=(r["log2"], yi),
                        arrowprops=dict(arrowstyle="->", color=col, lw=1.6))
            ax.plot(r["log2"], yi, "|", color=col, ms=12, mew=2)
        if r.get("nominal") is not None:
            ax.plot(r["nominal"], yi, "D", mfc="none", mec=ps.INK, ms=9, mew=1.0, zorder=3)
    ax.plot([], [], "o", color=ps.INK2, label="measured (95% CI)")
    ax.plot([], [], "D", mfc="none", mec=ps.INK, label="theory")
    ax.annotate("", xy=(0, 0), xytext=(0, 0))
    ax.plot([], [], marker=r"$\rightarrow$", ls="", color=ps.INK2, ms=12, label="lower bound (0 repeats)")
    ax.set_yticks(y)
    ax.set_yticklabels([r["label"] for r in rows], fontsize=9)
    ax.set_xlabel(r"$\log_2 G_{\rm eff}$")
    ax.set_title("Effective number of distinct values")
    ax.set_xlim(0, 56)
    ax.legend(fontsize=7.8, loc="lower left")
    ps.save(fig, "spacing_duplicates")


def geff_rows():
    from scipy import stats as st

    def garwood(k):
        lo = 0 if k == 0 else st.chi2.ppf(0.025, 2 * k) / 2
        hi = st.chi2.ppf(0.975, 2 * k + 2) / 2
        return lo, hi
    out = []
    spec = [("MIDSQ12345678", 1000, None), ("RANDU65536", 1e5, 13.0), ("RAND", 1e4, 22.0),
            ("MTsingle", 1e5, np.log2(0.75 * 2 ** 25)),
            ("RANDUsingle", 1e5, np.log2(2.638e7)), ("R250", 1e6, 32.0)]
    for tag, N, nom in spec:
        r = B[(B.tag == tag) & (B.N == N)].iloc[0]
        G = invert_geff(r.n_mean, r.dup_mean) * 2 ** r.j
        lo_k, hi_k = garwood(r.dup_sum)
        # scale CI approximately with D (G ~ 1/D)
        Glo = G * r.dup_sum / hi_k
        Ghi = G * r.dup_sum / lo_k if lo_k > 0 else np.inf
        out.append(dict(label=LABEL[tag], ckey=CKEY[tag], kind="measured", log2=np.log2(G),
                        lo=np.log2(Glo), hi=np.log2(Ghi), nominal=nom))
    out[0]["nominal"] = np.log2(100.0)
    # 53-bit generators: all blocks pooled
    for tag in ["MT", "RANNUM"]:
        g = B[B.tag == tag]
        P = g.pairs_equiv.sum(); k = g.dup_sum.sum()
        lo_k, hi_k = garwood(k)
        out.append(dict(label=LABEL[tag], ckey=CKEY[tag], kind="measured", log2=np.log2(P / k),
                        lo=np.log2(P / hi_k), hi=np.log2(P / lo_k), nominal=53.0))
    # permutation generators: zero repeats below the period, lower bound
    for tag, nom in [("RANDU", 29.0), ("MINSTD", np.log2(2 ** 31 - 2.0))]:
        g = B[(B.tag == tag) & (B.dup_sum == 0)]
        P = g.pairs_equiv.max()
        out.append(dict(label=LABEL[tag], ckey=CKEY[tag], kind="bound", log2=np.log2(P / 2.996),
                        nominal=nom))
    return out


# ---------------------------------------------------------------- Fig 6
def fig_lcg_period():
    d = pd.read_csv(D + "period_lcg.txt", sep=r"\s+", comment="#", header=None,
                    names=["gen", "seed", "tz", "lam", "mu", "sec"])
    fig, ax = plt.subplots(figsize=(7.2, 4.6))
    k = np.arange(0, 25)
    ax.plot(k, 29 - k, "--", color=ps.EXPECT, lw=1.2, label=r"theory $\log_2\lambda = 29-k$")
    for gen, mk, off in [("RANDU", "o", -0.12), ("RANDM", "s", 0.12)]:
        g = d[(d.gen == gen) & (d.seed > 0)]
        ax.plot(g.tz + off, np.log2(g.lam), mk, color=ps.COLORS[gen], ms=5.5, mfc="white" if gen == "RANDM" else ps.COLORS[gen],
                ls="", label=gen + r" (55 nonzero seeds)")
    g = d[(d.gen == "MINSTD")]
    ax.plot([0, 0], np.log2(g.lam), "^", color=ps.COLORS["MINSTD"], ms=7, ls="", label=r"MINSTD (seeds 1, 777): $2^{31}-2$")
    ax.annotate("seed 65536 = $2^{16}$:\n$\\lambda = 8192$", xy=(16, 13), xytext=(17.5, 20),
                arrowprops=dict(arrowstyle="->", color=ps.INK2), fontsize=9, color=ps.INK2)
    ax.annotate("odd seeds 3, 77, 777, 7777,\n12345, 77777, 314159:\nall $\\lambda = 2^{29}$", xy=(0.3, 29.2), xytext=(6.5, 28.3),
                arrowprops=dict(arrowstyle="->", color=ps.INK2), fontsize=9, color=ps.INK2)
    ax.text(24.5, 31, "seed 0: fixed point, $\\lambda=1$", ha="right", fontsize=9, color=ps.INK2)
    ax.set_xlabel(r"trailing zero bits $k$ of the seed (seed $=2^k\times$odd)")
    ax.set_ylabel(r"$\log_2$ period $\lambda$")
    ax.set_title("Exact periods by Brent's algorithm")
    ax.set_xlim(-1, 25); ax.set_ylim(0, 32)
    ax.legend(fontsize=8.5, loc="lower left")
    ps.save(fig, "spacing_lcg_period")


# ---------------------------------------------------------------- Fig 7
def fig_maps():
    c = pd.read_csv(D + "midsq_cycles.txt", sep=r"\s+", comment="#", header=None,
                    names=["id", "len", "min", "basin", "mean_tail", "max_tail"]).sort_values("basin", ascending=False)
    rh = pd.read_csv(D + "midsq_rho_hist.txt", sep=r"\s+", comment="#", header=None, names=["rho", "count"])
    L = pd.read_csv(D + "logis_survey.txt", sep=r"\s+", comment="#", header=None,
                    names=["idx", "bits", "x0", "mu", "lam", "cmin", "sec"])
    fig, axes = plt.subplots(2, 2, figsize=(10.5, 7.4))
    col_m = ps.COLORS["MIDSQ"]; col_l = ps.COLORS["LOGISTIC"]
    ax = axes[0, 0]
    top = c.head(8)
    lab = ["%d (len %d)" % (m, l) if l > 1 else "%d (fixed)" % m for m, l in zip(top["min"], top["len"])]
    ax.barh(range(len(top))[::-1], top.basin / 1e8, color=col_m)
    ax.set_yticks(range(len(top))[::-1]); ax.set_yticklabels(lab, fontsize=8.5)
    for i, v in zip(range(len(top))[::-1], top.basin / 1e8):
        ax.text(v + 0.005, i, "%.1f%%" % (100 * v) if v > 0.001 else "%.2f%%" % (100 * v), va="center", fontsize=8)
    ax.set_xlabel("fraction of all $10^8$ seeds")
    ax.set_title("Middle-square: where every seed ends (17 cycles)")
    ax.set_xlim(0, 0.62)
    ax.grid(axis="y", visible=False)

    ax = axes[0, 1]
    w = rh["count"] / 1e8
    bins = np.arange(0, 17500, 250)
    ax.hist(rh.rho, bins=bins, weights=w, color=col_m, histtype="stepfilled", alpha=0.85,
            label="all $10^8$ seeds")
    ax.axvline(np.sqrt(np.pi * 1e8 / 8) * 2, color=ps.EXPECT, ls="--", lw=1.2,
               label=r"random map: $E[\mu+\lambda]=\sqrt{\pi N/2}$")
    ax.axvline(10452, color=ps.INK2, ls=":", lw=1.3, label="seed 12345678: 10452")
    ax.set_xlabel(r"distinct values before repeating, $\mu+\lambda$")
    ax.set_ylabel("fraction of seeds per bin")
    ax.set_title("Middle-square: at most 16,903 values")
    ax.legend(fontsize=8)

    ax = axes[1, 0]
    g = L.groupby("lam").agg(n=("idx", "size")).reset_index().sort_values("n", ascending=False)
    lab = ["fixed point 0" if l == 1 else "cycle %s" % format(int(l), ",") for l in g.lam]
    frac = g.n / len(L)
    err = np.sqrt(frac * (1 - frac) / len(L))
    ax.barh(range(len(g))[::-1], frac, xerr=1.96 * err, color=col_l, capsize=2, ecolor=ps.INK2)
    ax.set_yticks(range(len(g))[::-1]); ax.set_yticklabels(lab, fontsize=8.5)
    for i, v in zip(range(len(g))[::-1], frac):
        ax.text(v + 0.04, i, "%.1f%%" % (100 * v), va="center", fontsize=8)
    ax.set_xlabel("fraction of %d random starting points" % len(L))
    ax.set_title("Logistic map in double precision: final attractor")
    ax.set_xlim(0, 0.85)
    ax.grid(axis="y", visible=False)

    ax = axes[1, 1]
    ax.hist(L.mu / 1e6, bins=np.linspace(0, 100, 51), color=col_l, histtype="stepfilled", alpha=0.85,
            label="transient $\\mu$, all starts")
    ax.axvline(14.324259, color=ps.INK2, ls=":", lw=1.3, label=r"seed 777: $\mu=14.3\times10^6$")
    ax.set_xlabel(r"transient length $\mu$ (millions of draws)")
    ax.set_ylabel("number of starts")
    ax.set_title("Logistic map: draws before the orbit is periodic")
    ax.legend(fontsize=8)
    ps.save(fig, "spacing_maps")




# ---------------------------------------------------------------- Fig 8
def fig_birthday():
    S = pd.read_csv(D + "bday_summary.csv")
    fig, axes = plt.subplots(1, 3, figsize=(12, 4.3), sharey=True)
    common = [("MT", "MT19937"), ("RANNUM", "RANNUM"), ("R250", "R250"), ("MINSTD", "MINSTD"),
              ("RANDU", "RANDU"), ("RANDM", "RANDM"), ("LOGISTIC", "LOGISTIC")]
    extra1 = [("RAND", "RAND"), ("MTsingle", "MT19937 (REAL*4)")]
    for ax, t in zip(axes, [1, 2, 3]):
        lst = common + (extra1 if t == 1 else [])
        for tag, lab in lst:
            g = S[(S.tag == tag) & (S.t == t)].sort_values("n")
            ck = "MT19937" if tag in ("MT", "MTsingle") else tag
            ls = "--" if tag == "MTsingle" else "-"
            mk = "s" if tag == "MTsingle" else "o"
            ax.plot(g.n, g.Ymean / g.lam, ls, marker=mk, ms=3.5, color=ps.COLORS[ck], label=lab,
                    mfc="white" if tag == "MTsingle" else ps.COLORS[ck])
        ax.axhline(1.0, color=ps.EXPECT, ls="--", lw=1.0)
        if t > 1:
            for rho, ck in [(2.0 ** 31 - 2, "MINSTD"), (2.0 ** 29, "RANDU")]:
                ax.axvline((8 * rho) ** (1 / 3), color=ps.COLORS[ck], ls=":", lw=1.2)
        ax.set_xscale("log", base=2)
        ax.set_yscale("log")
        ax.set_xlabel("birthdays per sample $n$")
        ax.set_title({1: "t = 1 (one number per point)", 2: "t = 2 (pairs)", 3: "t = 3 (triplets)"}[t])
    axes[0].set_ylabel(r"mean collisions $\bar Y\,/\,\lambda$, $\lambda=n^3/4k=2$")
    axes[0].legend(fontsize=7.8, loc="upper left", ncol=1)
    axes[1].text(0.03, 0.97, r"dotted: $n=(8\rho)^{1/3}$, i.e. $k=\rho$" "\n(MINSTD orange, RANDU blue)",
                 transform=axes[1].transAxes, va="top", fontsize=8, color=ps.INK2)
    axes[0].set_ylim(0.5, 3000)
    ps.save(fig, "spacing_birthday")


if __name__ == "__main__":
    which = sys.argv[1:] or ["dist", "scan", "breakmap", "minsp", "dups", "lcg", "maps", "birthday"]
    for w in which:
        {"dist": fig_dist, "scan": fig_scan, "breakmap": fig_breakmap, "minsp": fig_minsp,
         "dups": fig_dups, "lcg": fig_lcg_period, "maps": fig_maps, "birthday": fig_birthday}[w]()
        print("done", w)
