# Report figures for the 'digits' family (decimal digits and bit-level structure)
import sys, json
import numpy as np
sys.path.insert(0, "/Users/noah/Downloads/ComputationalPhysicsPHYS25000/Project1/P1_Claude")
sys.path.insert(0, "/Users/noah/Downloads/ComputationalPhysicsPHYS25000/Project1/P1_Claude/tests/digits")
import plotstyle as ps
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, ListedColormap
from matplotlib.patches import Patch
from scipy import stats
from dig_analysis import load as dload, chi2 as dchi2, signif
from bit_analysis import load as bload, stats as bstats

ps.apply()
D = "/Users/noah/Downloads/ComputationalPhysicsPHYS25000/Project1/P1_Claude/tests/digits/out/"
NAMES = ['RANDU', 'RANDM', 'RAND', 'MINSTD', 'R250', 'RANNUM', 'MT19937', 'MIDSQ', 'LOGISTIC']
GID = {n: i + 1 for i, n in enumerate(NAMES)}
C = ps.COLORS
CRIT_HI = stats.chi2.isf(1e-4, 9)
CRIT_LO = stats.chi2.ppf(1e-4, 9)

# diverging map: blue arm, neutral gray midpoint, red arm (validated ramps)
DIV = LinearSegmentedColormap.from_list(
    "div", ["#104281", "#2a78d6", "#86b6ef", "#f0efec", "#f2a3a2", "#e34948", "#8f1f1f"])
BIN = ListedColormap(["#fcfcfb", ps.INK])


def sig_table(fn, kind, pmax):
    tab = dload(fn)
    s = np.full(pmax, np.nan)
    x2 = np.full(pmax, np.nan)
    for p in range(1, pmax + 1):
        if p in tab[kind]:
            x2[p - 1] = dchi2(tab[kind][p])
            s[p - 1] = signif(x2[p - 1])[0]
    return s, x2


# ---------------------------------------------------------------- Figure 1
def fig_heatmap():
    rows = NAMES + ['LOGISTIC*']
    panels = [("Exact digits of the double value", 0, 'INT', 24),
              ("Naive digits, double arithmetic", 0, 'NAI', 18),
              ("Exact digits of the REAL*4 value", 1, 'STO', 16),
              ("Naive digits, REAL*4 arithmetic", 1, 'NAI', 12)]
    fig, axes = plt.subplots(2, 2, figsize=(10, 7.8),
                             gridspec_kw=dict(width_ratios=[24, 18]))
    cap = 15.0
    for ax, (title, prec, kind, pmax) in zip(axes.flatten(), panels):
        M = np.zeros((len(rows), pmax))
        for r, name in enumerate(rows):
            if name == 'LOGISTIC*':
                fn = D + f"dig_g9_p{prec}_N1e7.txt"
            else:
                fn = D + f"dig_g{GID[name]}_p{prec}.txt"
            s, _ = sig_table(fn, kind, pmax)
            M[r] = s
        im = ax.imshow(np.clip(M, -cap, cap), cmap=DIV, vmin=-cap, vmax=cap,
                       aspect='auto', interpolation='nearest')
        ax.set_title(title)
        ax.set_xticks(range(0, pmax, 2))
        ax.set_xticklabels(range(1, pmax + 1, 2))
        ax.set_yticks(range(len(rows)))
        if pmax in (24, 16):
            ax.set_yticklabels([r if r != 'LOGISTIC*' else 'LOGISTIC, $10^7$' for r in rows])
        else:
            ax.set_yticklabels([])
        ax.set_xlabel("decimal position $p$")
        ax.grid(False)
        for k in range(len(rows) + 1):
            ax.axhline(k - 0.5, color="#fcfcfb", lw=0.6)
        for k in range(pmax + 1):
            ax.axvline(k - 0.5, color="#fcfcfb", lw=0.6)
    fig.subplots_adjust(right=0.88, hspace=0.32, wspace=0.06)
    cb = fig.colorbar(im, ax=axes, shrink=0.75, pad=0.02)
    cb.set_label(r"signed significance $\pm\log_{10}(1/p_{\rm tail})$, capped at 15")
    cb.ax.text(1.0, -0.04, "too uniform", transform=cb.ax.transAxes, ha='center', va='top', fontsize=8, color=ps.INK2)
    cb.ax.text(1.0, 1.03, "non-uniform", transform=cb.ax.transAxes, ha='center', va='bottom', fontsize=8, color=ps.INK2)
    import os
    os.makedirs(ps.FIGDIR, exist_ok=True)
    fig.savefig(ps.FIGDIR + "/digits_heatmap.png", bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------- Figure 2
def fig_depth():
    pr = json.load(open(D + "predict.json"))
    P = np.arange(1, 61)
    fig, axes = plt.subplots(2, 3, figsize=(10, 6.6), sharey=True)
    axes = axes.flatten()

    def band(ax, pmax):
        ax.axhspan(CRIT_LO, CRIT_HI, color="#f0efec", zorder=0)
        ax.set_yscale('log')
        ax.set_ylim(1e-2, 3e9)
        ax.set_xlim(0, pmax + 1)
        ax.set_xlabel("decimal position $p$")

    def meas(ax, name, kind, pmax, prec=0, fn=None, marker='o', label=None):
        fn = fn or D + f"dig_g{GID[name]}_p{prec}.txt"
        _, x2 = sig_table(fn, kind, pmax)
        ax.plot(np.arange(1, pmax + 1), x2, marker=marker, ms=4.5, lw=1.0, color=C[name],
                label=label or name, mfc=C[name] if marker == 'o' else 'none')

    # RAND
    ax = axes[0]; band(ax, 30)
    meas(ax, 'RAND', 'STO', 30, label='RAND measured')
    lam = np.array(pr['uniform']['22'][:30])
    ax.plot(P[:30], 9 + lam, ls='--', color=ps.EXPECT, lw=1.2, label=r'$9+\lambda_p$, $B=22$')
    ax.set_title("RAND ($2^{-22}$ grid)"); ax.set_ylabel(r"$\chi^2$ (9 dof)")
    ax.legend(loc='upper left')
    # RANDU / RANDM
    ax = axes[1]; band(ax, 30)
    meas(ax, 'RANDU', 'STO', 30, label='RANDU')
    meas(ax, 'RANDM', 'STO', 30, marker='s', label='RANDM')
    lr = np.array([v if v is not None else np.nan for v in pr['randu'][:30]])
    m = P[:30] >= 8
    ax.plot(P[:30][m], lr[m], ls='--', color=ps.EXPECT, lw=1.2, label=r'full-period $\lambda_p$')
    ax.set_title("RANDU and RANDM ($k/2^{31}$)")
    ax.legend(loc='upper left')
    # R250
    ax = axes[2]; band(ax, 30)
    meas(ax, 'R250', 'STO', 30, label='R250 measured')
    lam = np.array(pr['uniform']['32'][:30])
    ax.plot(P[:30], 9 + lam, ls='--', color=ps.EXPECT, lw=1.2, label=r'$9+\lambda_p$, $B=32$')
    ax.set_title("R250 ($2^{-32}$ grid)")
    ax.legend(loc='upper left')
    # MT / RANNUM
    ax = axes[3]; band(ax, 60)
    meas(ax, 'MT19937', 'STO', 60, label='MT19937')
    meas(ax, 'RANNUM', 'STO', 60, marker='s', label='RANNUM')
    lam = np.array(pr['uniform']['53'][:60])
    ax.plot(P, 9 + lam, ls='--', color=ps.EXPECT, lw=1.2, label=r'$9+\lambda_p$, $B=53$')
    ax.set_title("MT19937 and RANNUM ($2^{-53}$ grid)"); ax.set_ylabel(r"$\chi^2$ (9 dof)")
    ax.legend(loc='upper left')
    # MINSTD intended vs stored
    ax = axes[4]; band(ax, 60)
    _, x2i = sig_table(D + "dig_g4_p0.txt", 'INT', 60)
    _, x2s = sig_table(D + "dig_g4_p0.txt", 'STO', 60)
    ax.plot(P, x2s, marker='o', ms=4.5, lw=1.0, color=C['MINSTD'], label='stored double')
    ax.plot(P, x2i, marker='s', ms=4.5, lw=1.0, color=C['MINSTD'], mfc='none', label=r'exact $k/(2^{31}-1)$')
    ax.plot(P, 9 + lam, ls='--', color=ps.EXPECT, lw=1.2, label=r'$9+\lambda_p$, $B=53$')
    ax.set_title("MINSTD: rational value vs double")
    ax.legend(loc='upper left')
    # single precision
    ax = axes[5]; band(ax, 20)
    for name, mk in (('MT19937', 'o'), ('MINSTD', 's'), ('RANDU', '^')):
        _, x2 = sig_table(D + f"dig_g{GID[name]}_p1.txt", 'STO', 20)
        ax.plot(np.arange(1, 21), x2, marker=mk, ms=4.5, lw=1.0, color=C[name], label=name,
                mfc=C[name] if mk == 'o' else 'none')
    ax.plot(np.arange(1, 21), 9 + np.array(pr['single'][:20]), ls='--', color=ps.EXPECT, lw=1.2,
            label=r'$9+\lambda_p$, REAL*4')
    ax.set_title("Values stored in REAL*4")
    ax.legend(loc='upper left')
    ps.save(fig, "digits_depth")


# ---------------------------------------------------------------- Figure 3
def fig_bitplane():
    show = ['RANDU', 'RANDM', 'RAND', 'MINSTD', 'R250', 'MT19937', 'LOGISTIC', 'ZOOM']
    nrow = 128
    fig, axes = plt.subplots(2, 4, figsize=(10, 7.6))
    for ax, name in zip(axes.flatten(), show):
        if name == 'ZOOM':
            m = np.loadtxt(D + "dump_g1_p0.txt", dtype=np.int64)[:48]
            jj = np.arange(17, 32)
            bits = ((m[:, None] >> (53 - jj)[None, :]) & 1)
            ax.imshow(bits, cmap=BIN, vmin=0, vmax=1, aspect='auto', interpolation='nearest',
                      extent=(16.5, 31.5, 47.5, -0.5))
            ax.set_title("RANDU, bits 17 to 31 (zoom)")
            ax.set_xticks([17, 20, 25, 28, 31])
            ax.set_xlabel("bit position $j$")
            ax.grid(False)
            continue
        m = np.loadtxt(D + f"dump_g{GID[name]}_p0.txt", dtype=np.int64)[:nrow]
        bits = ((m[:, None] >> (53 - np.arange(1, 54))[None, :]) & 1)
        ax.imshow(bits, cmap=BIN, vmin=0, vmax=1, aspect='auto', interpolation='nearest',
                  extent=(0.5, 53.5, nrow - 0.5, -0.5))
        ax.set_title(name)
        ax.grid(False)
        ax.set_xticks([1, 10, 20, 31, 40, 53])
        if name in ('RANDU', 'R250'):
            ax.set_ylabel("output index $n$")
        if name in ('R250', 'MT19937', 'LOGISTIC'):
            ax.set_xlabel("bit position $j$")
    ps.save(fig, "digits_bitplane")


# ---------------------------------------------------------------- Figure 4
def fig_bitstats():
    rows = NAMES
    fig, axes = plt.subplots(2, 2, figsize=(10, 7.0), sharex=True, sharey=True)
    cap = 10.0
    for ci, prec in enumerate((0, 1)):
        F = np.zeros((9, 53)); R = np.zeros((9, 53)); K = np.zeros((9, 53), bool)
        for r, name in enumerate(rows):
            d = bload(GID[name], prec)
            f, zf, const, cor, zc = bstats(d)
            # signed significance of the frequency test and of the lag-1 correlation
            pf = 2 * stats.norm.sf(np.abs(zf))
            F[r] = np.sign(zf) * np.minimum(-np.log10(np.maximum(pf, 1e-300)), 300)
            zz = np.nan_to_num(zc[:, 0])
            pc = 2 * stats.norm.sf(np.abs(zz))
            R[r] = np.sign(zz) * np.minimum(-np.log10(np.maximum(pc, 1e-300)), 300)
            K[r] = const
        for ri, (M, ttl) in enumerate(((F, "frequency of ones"), (R, "lag-1 autocorrelation"))):
            ax = axes[ri, ci]
            Mc = np.clip(M, -cap, cap)
            Mc[K] = np.nan
            cmap = DIV.copy(); cmap.set_bad(ps.MUTED)
            im = ax.imshow(Mc, cmap=cmap, vmin=-cap, vmax=cap, aspect='auto', interpolation='nearest',
                           extent=(0.5, 53.5, 8.5, -0.5))
            ax.set_title(f"{ttl}, {'double' if prec == 0 else 'REAL*4'}")
            ax.set_yticks(range(9)); ax.set_yticklabels(rows)
            ax.grid(False)
            for k in range(10):
                ax.axhline(k - 0.5, color="#fcfcfb", lw=0.8)
            if ri == 1:
                ax.set_xlabel("bit position $j$")
            ax.set_xticks([1, 10, 20, 24, 31, 40, 53])
    fig.subplots_adjust(right=0.87)
    cb = fig.colorbar(im, ax=axes, shrink=0.8, pad=0.02)
    cb.set_label(r"signed $\log_{10}(1/p)$, capped at 10")
    fig.legend(handles=[Patch(color=ps.MUTED, label="constant bit (stuck at 0 or 1)")],
               loc='lower center', bbox_to_anchor=(0.45, -0.02), ncol=1)
    fig.savefig(ps.FIGDIR + "/digits_bitstats.png", bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------- Figure 5
def parse_period(g):
    bit = {}; dig = {}; lag = {}
    for line in open(D + f"period_g{g}.txt"):
        t = line.split()
        if t[0] == 'BITPERIOD':
            bit[int(t[1])] = int(t[2])
        elif t[0] == 'DIGPERIOD':
            dig[int(t[1])] = (int(t[2]), int(t[3]), int(t[4]))
        elif t[0] == 'RETURN':
            lag.setdefault('ret', {})[int(t[1])] = int(t[2])
    return bit, dig, lag


def fig_lcgbits():
    pred = json.load(open(D + "bitcor_predict.json"))
    fig, axes = plt.subplots(2, 2, figsize=(10, 7.0))
    # (a) bit periods
    ax = axes[0, 0]
    for g, name, mk in ((1, 'RANDU', 'o'), (2, 'RANDM', 's')):
        bit, dig, ret = parse_period(g)
        # bits from brute force (state bit i = output bit 31 - i), high bits from residue return times
        J = []; Pm = []
        for i, per in bit.items():
            J.append(31 - i); Pm.append(per)
        for i in range(23, 31):        # state bit i: period = return time of residue mod 2^(i+1)
            J.append(31 - i); Pm.append(ret['ret'][i + 1])
        o = np.argsort(J)
        ax.plot(np.array(J)[o], np.log2(np.array(Pm)[o]), marker=mk, ms=5, lw=0, color=C[name],
                mfc=C[name] if mk == 'o' else 'none', label=name + " measured")
    jj = np.arange(1, 29)
    ax.plot(jj, 30 - jj, ls='--', color=ps.EXPECT, lw=1.2, label=r"$2^{30-j}$ ($j\leq 28$)")
    ax.set_xlabel("bit position $j$ of $x=k/2^{31}$"); ax.set_ylabel(r"$\log_2$(period of bit $j$)")
    ax.set_title("Period of each bit")
    ax.legend()
    # (b) digit periods
    ax = axes[0, 1]
    for g, name, mk in ((1, 'RANDU', 'o'), (2, 'RANDM', 's')):
        bit, dig, ret = parse_period(g)
        pp = sorted(p for p in dig if dig[p][0] > 0)
        ax.plot(pp, [np.log2(dig[p][0]) for p in pp], marker=mk, ms=5, lw=0, color=C[name],
                mfc=C[name] if mk == 'o' else 'none', label=name + " measured")
    pp = np.arange(9, 29)
    ax.plot(pp, 30 - pp, ls='--', color=ps.EXPECT, lw=1.2, label=r"$2^{30-p}$ ($p\leq 28$)")
    ax.set_xlabel("decimal position $p$"); ax.set_ylabel(r"$\log_2$(period of digit $p$)")
    ax.set_title("Period of each decimal digit")
    ax.legend()
    # (c) lag-1 correlation vs prediction
    ax = axes[1, 0]
    for g, name, mk in ((1, 'RANDU', 'o'), (2, 'RANDM', 's')):
        d = bload(g, 0); f, zf, const, cor, zc = bstats(d)
        j = np.arange(1, 31)
        ax.plot(j, cor[:30, 0], marker=mk, ms=5, lw=0, color=C[name], mfc=C[name] if mk == 'o' else 'none',
                label=name + " measured")
        pj = np.array([int(k) for k in pred[str(g)].keys()])
        pv = np.array([np.nan if v is None else v for v in pred[str(g)].values()])
        ax.plot(pj, pv, ls='--', color=C[name], lw=1.0, alpha=0.8, label=name + " full-period theory")
    ax.axhline(1 / 3, color=ps.EXPECT, lw=0.8, ls=':')
    ax.text(1.5, 1 / 3 + 0.04, r"$1/3$", color=ps.INK2, fontsize=9)
    ax.set_xlabel("bit position $j$"); ax.set_ylabel(r"lag-1 autocorrelation of $b_j$")
    ax.set_title("Lag-1 autocorrelation of each bit")
    ax.set_ylim(-1.1, 0.75)
    ax.legend(loc='lower left', fontsize=8)
    # (d) linear complexity
    ax = axes[1, 1]
    for name, g, js in (('RANDU', 1, [1, 14, 16, 18, 20, 24, 28]), ('RANDM', 2, [1, 16, 20, 24])):
        lc = []
        for j in js:
            for line in open(D + f"bm_g{g}_j{j}.txt"):
                if line.startswith('FINAL'):
                    lc.append(int(line.split()[1]))
        ax.plot(js, np.log2(lc), marker='o' if name == 'RANDU' else 's', ms=6, lw=0, color=C[name],
                mfc=C[name] if name == 'RANDU' else 'none', label=name + " (Berlekamp-Massey)")
    jj = np.arange(1, 29)
    ax.plot(jj, np.log2(2.0 ** (29 - jj) + 1), ls='--', color=ps.EXPECT, lw=1.2, label=r"$2^{29-j}+1$")
    ax.axhline(np.log2(25000), color=ps.MUTED, lw=0.8, ls=':')
    ax.text(18.5, np.log2(25000) + 0.4, r"random: $L/2$ for $L=5\times10^4$", fontsize=8, color=ps.INK2)
    ax.set_xlabel("bit position $j$"); ax.set_ylabel(r"$\log_2$(linear complexity)")
    ax.set_title("Linear complexity of each bit stream")
    ax.set_ylim(0, 18)
    ax.legend(loc='lower left')
    ps.save(fig, "digits_lcgbits")


# ---------------------------------------------------------------- Figure 6
def fig_xor():
    fig, axes = plt.subplots(2, 2, figsize=(10, 7.2))
    # (a) R250 recurrence per bit
    ax = axes[0, 0]
    for name, mk in (('R250', 'o'), ('MT19937', 's'), ('MINSTD', '^')):
        d = bload(GID[name], 0)
        fr = d['v250'][:32] / d['n250']
        ax.plot(np.arange(1, 33), fr, marker=mk, ms=5, lw=0, color=C[name],
                mfc=C[name] if mk == 'o' else 'none', label=name)
    ax.axhline(0.5, color=ps.EXPECT, ls='--', lw=1.0, label='random: 1/2')
    ax.set_ylim(-0.05, 0.65)
    ax.set_xlabel("bit position $j$")
    ax.set_ylabel("fraction of $n$ violating the identity")
    ax.set_title(r"R250: $b_j(n)\oplus b_j(n-103)\oplus b_j(n-250)=0$")
    ax.legend(loc='center right')
    # (b) logistic tent-shift identity
    ax = axes[0, 1]
    for name, prec, mk, lab in (('LOGISTIC', 0, 'o', 'LOGISTIC, double'), ('LOGISTIC', 1, 's', 'LOGISTIC, REAL*4'),
                                ('MT19937', 0, '^', 'MT19937, double')):
        d = bload(GID[name], prec)
        fr = d['vtent'][:52] / (d['N'] - 1)
        j = np.arange(1, 53)
        nz = fr > 0
        col = C[name]
        ax.plot(j[nz], fr[nz], marker=mk, ms=4.5, lw=0.8, color=col, mfc=col if mk == 'o' else 'none', label=lab)
        if (~nz).any() and prec == 0 and name == 'LOGISTIC':
            ax.plot(j[~nz], np.full((~nz).sum(), 1.5e-8), marker='v', ms=4.5, lw=0, color=col)
    ax.axhline(1.0 / 5e7, color=ps.MUTED, ls=':', lw=0.8)
    ax.text(6.5, 4.5e-8, "zero violations (markers at bottom)", fontsize=8, color=ps.INK2)
    ax.set_yscale('log'); ax.set_ylim(1e-8, 1.5)
    ax.set_xlabel("bit position $j$")
    ax.set_ylabel("fraction of $n$ violating the identity")
    ax.set_title(r"Logistic: $b_j(n)\oplus b_{j+1}(n-1)\oplus b_1(n-1)=0$")
    ax.legend(loc='lower right', fontsize=8)
    # (c) MINSTD copy identity
    ax = axes[1, 0]
    J = np.arange(32, 54)
    for name, mk in (('MINSTD', 'o'), ('MT19937', 's')):
        d = bload(GID[name], 0)
        fr = d['vcopy'][31:53] / d['N']
        ax.plot(J, np.where(fr > 0, fr, np.nan), marker=mk, ms=5, lw=0.8, color=C[name],
                mfc=C[name] if mk == 'o' else 'none', label=name)
    ax.plot(J, (1 / 3) * 2.0 ** -(53 - J), ls='--', color=ps.EXPECT, lw=1.2, label=r"$\frac{1}{3}2^{-(53-j)}$ (rounding carry)")
    ax.set_yscale('log')
    ax.set_xlabel("bit position $j$")
    ax.set_ylabel(r"fraction with $b_j\neq b_{j-31}$")
    ax.set_title(r"MINSTD double: bits 32 to 53 copy bits 1 to 22")
    ax.legend(loc='lower right')
    # (d) lag-pair scan
    ax = axes[1, 1]
    for name, mk in (('MT19937', 'o'), ('RANDU', 's'), ('MINSTD', '^')):
        z = []; c = []
        for line in open(D + f"lag_g{GID[name]}.txt"):
            t = line.split()
            if t[0] == 'HIST':
                z.append(float(t[1])); c.append(int(t[2]))
        z = np.array(z); c = np.array(c, float); m = np.abs(z) < 10
        ax.plot(z[m], c[m], marker=mk, ms=4, lw=0.8, color=C[name], mfc=C[name] if mk == 'o' else 'none',
                label=name)
    zz = np.linspace(-6, 6, 241)
    ax.plot(zz, 179700 * 0.25 * stats.norm.pdf(zz), ls='--', color=ps.EXPECT, lw=1.2, label=r"$N(0,1)$")
    ax.set_yscale('log'); ax.set_ylim(0.5, 2e6); ax.set_xlim(-6.5, 6.5)
    ax.set_xlabel(r"$z$ of the three-term XOR count, all $1\leq s<t\leq600$")
    ax.set_ylabel("number of lag pairs (bin 0.25)")
    ax.set_title("Blind scan of 179,700 lag pairs")
    ax.text(-6.2, 1.2e6, "R250 (off scale): only (103, 250) and (206, 500)\nare hit, both with $z=-800$", fontsize=8,
            color=ps.INK2, va='top')
    ax.legend(loc='upper right', fontsize=8)
    ps.save(fig, "digits_xor")


# ---------------------------------------------------------------- Figure 7
def fig_lincomp():
    fig, ax = plt.subplots(figsize=(7.5, 5.0))
    series = [('RANNUM', 6, 1, '-', 'RANNUM bit 1'), ('R250', 5, 1, '-', 'R250 bit 1'),
              ('MT19937', 7, 1, '-', 'MT19937 bit 1'), ('RANDU', 1, 16, '-', 'RANDU bit 16'),
              ('RANDU', 1, 24, '--', 'RANDU bit 24')]
    for name, g, j, ls, lab in series:
        a = np.loadtxt(D + f"bm_g{g}_j{j}.txt", comments=('#', 'F'))
        ax.plot(a[:, 0], np.maximum(a[:, 1], 0.8), color=C[name], ls=ls, lw=1.6, label=lab)
        ax.text(5.2e4, a[-1, 1] * (1.15 if name == 'RANNUM' else 1.0), f"{int(a[-1, 1])}",
                va='bottom' if name == 'RANNUM' else 'center', fontsize=8, color=ps.INK2)
    n = np.logspace(1, np.log10(5e4), 100)
    ax.plot(n, n / 2, ls=':', color=ps.EXPECT, lw=1.2, label=r"$n/2$ (random)")
    ax.set_xscale('log'); ax.set_yscale('log')
    ax.set_xlim(10, 9e4); ax.set_ylim(1, 5e4)
    ax.set_xlabel("number of bits $n$")
    ax.set_ylabel("linear complexity $L(n)$")
    ax.set_title("Berlekamp-Massey linear complexity profiles")
    ax.legend(loc='upper left')
    ps.save(fig, "digits_lincomp")


if __name__ == "__main__":
    which = sys.argv[1:] or ['heatmap', 'depth', 'bitplane', 'bitstats', 'lcgbits', 'xor', 'lincomp']
    for w in which:
        globals()['fig_' + w]()
        print('made', w)
