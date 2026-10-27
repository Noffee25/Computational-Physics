"""figures.py : report figures of the battery family (battery_*.png)."""
import sys
sys.path.insert(0, "/Users/noah/Downloads/ComputationalPhysicsPHYS25000/Project1/P1_Claude")
sys.path.insert(0, "/Users/noah/Downloads/ComputationalPhysicsPHYS25000/Project1/P1_Claude/tests/battery")
import json, glob
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import plotstyle as ps
from readc3 import readc3
from summary import NAMES, load_lag, PERIOD, BASE, S

ps.apply()
C = ps.COLORS
DIV = LinearSegmentedColormap.from_list('div', [(0, '#123f78'), (0.25, '#8fb3e3'), (0.5, '#f0efec'), (0.75, '#efa3a2'), (1, '#8f2321')])


def fig_lagcorr():
    ded = json.load(open(BASE + 'out/dedekind_lags.json'))
    panels = [(1, 'out/lag/lag_1.txt', 'RANDU, N = 2^29 (full period)'),
              (2, 'out/lag/lag_2.txt', 'RANDM, N = 2^29 (full period)'),
              (4, 'out/lag/lag_4.txt', 'MINSTD, N = 2^29'),
              (5, 'out/lag/lag_5.txt', 'R250, N = 2^29'),
              (7, 'out/lag/lag_7.txt', 'MT19937, N = 2^29'),
              (9, 'out/lagN/lag_9_23.txt', 'LOGISTIC, N = 2^23 (before its cycle)')]
    fig, axes = plt.subplots(2, 3, figsize=(10, 5.6), sharex=True, sharey=True)
    for ax, (g, fn, title) in zip(axes.ravel(), panels):
        N, d = load_lag(BASE + fn)
        k = d[:, 0]; z = d[:, 1] / np.sqrt(N - k)
        name = NAMES[g]
        ax.axhspan(-3, 3, color=ps.GRID, alpha=0.45, lw=0)
        ax.axhline(3, color=ps.EXPECT, ls='--', lw=0.9); ax.axhline(-3, color=ps.EXPECT, ls='--', lw=0.9)
        ax.vlines(k, 0, z, color=C[name], lw=0.8)
        ax.plot(k, z, 'o', ms=2.2, color=C[name], label='measured')
        if name in ded:
            zp = np.array(ded[name]['pred']) * np.sqrt(N - k)
            ax.plot(k, zp, 'o', ms=4.5, mfc='none', mec=ps.EXPECT, mew=0.7, label='Dedekind-sum theory')
            ax.legend(loc='lower right', fontsize=8)
        Q = (z**2).sum()
        ax.text(0.02, 0.96, r'$Q_{300}$ = %.1f' % Q, transform=ax.transAxes, va='top', fontsize=9, color=ps.INK)
        ax.set_title(title)
        ax.set_ylim(-5.2, 5.2)
    for ax in axes[1]: ax.set_xlabel('lag k')
    for ax in axes[:, 0]: ax.set_ylabel(r'$z_k = \sqrt{N-k}\; r_k$')
    return ps.save(fig, 'battery_lagcorr')


def cpz(fn):
    d = np.loadtxt(fn)
    n = d[:, 0]; z = d[:, 1] / np.sqrt(n - 1)
    return n, z


def fig_convergence():
    ex = S['exact']
    fig, axes = plt.subplots(1, 3, figsize=(10, 3.8), sharey=True)
    specs = [('RANDM', 2, ['777', '1', '314159', '65539'], ex['RANDM'][1]['exact']),
             ('RANDU', 1, ['777', '1', '314159', '65539'], ex['RANDU'][1]['exact']),
             ('MT19937', 7, ['777', '1', '314159'], 0.0)]
    for ax, (name, g, seeds, rex) in zip(axes, specs):
        Ngrid = np.logspace(np.log10(3e4), np.log10(2**29), 200)
        ax.fill_between(Ngrid, -3, 3, color=ps.GRID, alpha=0.45, lw=0)
        for i, s in enumerate(seeds):
            fn = BASE + (f'out/lag/lag_{g}.txt.cp' if s == '777' else f'out/lag/lagconv_{g}_{s}.txt.cp')
            n, z = cpz(fn)
            ax.plot(n, z, color=C[name], lw=1.2, alpha=1.0 if i == 0 else 0.55,
                    label=f'{name}, {len(seeds)} seeds' if i == 0 else None)
        if rex != 0.0:
            ax.plot(Ngrid, rex * np.sqrt(Ngrid), color=ps.EXPECT, ls='--', lw=1.2,
                    label=r'exact $r_1\sqrt{N}$ (full period)')
            ax.axvline(2**29, color=ps.MUTED, lw=0.8, ls=':')
            ax.text(2**29, -4.4, ' P = 2^29', color=ps.INK2, fontsize=8, ha='right')
        else:
            ax.plot(Ngrid, 0 * Ngrid, color=ps.EXPECT, ls='--', lw=1.2, label='ideal: 0')
        ax.axhline(3, color=ps.EXPECT, lw=0.6); ax.axhline(-3, color=ps.EXPECT, lw=0.6)
        ax.set_xscale('log'); ax.set_xlabel('number of draws N')
        ax.set_title(f'{name} lag-1 correlation')
        ax.legend(loc='upper left', fontsize=8)
        ax.set_ylim(-4.8, 6)
    axes[0].set_ylabel(r'$z_1 = \sqrt{N}\, r_1$')
    return ps.save(fig, 'battery_lag1_convergence')


def fig_c3maps():
    cols = [(1, 'RANDU'), (2, 'RANDM'), (5, 'R250'), (7, 'MT19937')]
    fig, axes = plt.subplots(2, 4, figsize=(12, 6.4))
    th = {nm: np.load(BASE + f'out/c3theory_{nm}.npz') for nm in ('RANDU', 'RANDM')}
    for c, (g, name) in enumerate(cols):
        K, N, T1, T2, Z1, Z2 = readc3(BASE + f'out/c3/c3_{g}.bin')
        jj = np.arange(1, K + 1); cnt = N - np.maximum.outer(jj, jj)
        for r, (Z, lab, key) in enumerate(((Z1, r'$\langle P_1 P_1 P_1\rangle$ (odd)', 'T1'),
                                           (Z2, r'$\langle P_1 P_1 P_2\rangle$ (even)', 'T2'))):
            ax = axes[r, c]
            M = np.where(np.triu(np.ones((K, K), bool), 1), Z, np.nan)   # row j < col k
            im = ax.imshow(np.clip(M, -6, 6), origin='lower', cmap=DIV, vmin=-6, vmax=6,
                           extent=[0.5, K + 0.5, 0.5, K + 0.5], interpolation='nearest')
            ax.grid(False)
            ax.set_title(f'{name}: {lab}', fontsize=10)
            # theory-predicted cells (|z| > 3 at this N)
            marks = []
            if name in th:
                zp = th[name][key] * np.sqrt(cnt)
                iu = np.triu_indices(K, 1)
                sel = np.abs(zp[iu]) > 3
                marks = [(iu[0][t] + 1, iu[1][t] + 1, zp[iu][t], Z[iu][t]) for t in np.where(sel)[0]]
            if name == 'R250' and key == 'T1':
                marks = [(103, 250, None, Z[102, 249])]
            marks.sort(key=lambda t: -abs(t[3]))
            slots = [(12, 228), (12, 186), (12, 144)]
            for (j, k, zpred, zmeas), (tx, ty) in zip(marks, slots):
                ax.plot(k, j, 'o', ms=11, mfc='none', mec=ps.INK, mew=1.0)
                txt = f'({j},{k}): z = {zmeas:.1f}'
                ax.annotate(txt, (k, j), xytext=(tx, ty), fontsize=7.5, color=ps.INK,
                            arrowprops=dict(arrowstyle='-', color=ps.INK2, lw=0.6, shrinkB=6))
            if c == 0: ax.set_ylabel('smaller lag j')
            if r == 1: ax.set_xlabel('larger lag k')
    cb = fig.colorbar(im, ax=axes, shrink=0.75, pad=0.015)
    cb.set_label('z (clipped at 6)')
    fig.savefig(ps.FIGDIR + '/battery_c3maps.png', bbox_inches='tight')
    plt.close(fig)
    return ps.FIGDIR + '/battery_c3maps.png'


def tent_theory():
    u = (np.arange(4_000_000) + 0.5) / 4_000_000
    T = lambda v: 1 - np.abs(1 - 2 * v)
    P = {1: lambda v: np.sqrt(3) * (2 * v - 1), 2: lambda v: np.sqrt(5) * (6 * v * v - 6 * v + 1)}
    x0, x1, x2 = u, T(u), T(T(u))
    return {f'{p}{q}{r}': float(np.mean(P[p](x0) * P[q](x1) * P[r](x2))) for p in (1, 2) for q in (1, 2) for r in (1, 2)}


def fig_triples():
    tri = S['triples']; th = S['triples_theory']
    combos = ['111', '112', '121', '211', '122', '212', '221', '222']
    panels = [('RANDU', 'RANDU_1_2', 'RANDU, lags (1,2)', th['RANDU_1_1_2']),
              ('R250', '5_103_250', 'R250, lags (103,250)', th['R250_103_250']),
              ('LOGISTIC', '9_1_2', 'LOGISTIC, lags (1,2)', tent_theory()),
              ('MT19937', '7_1_2', 'MT19937, lags (1,2)', None)]
    keymap = {'RANDU_1_2': '1_1_2'}
    fig, axes = plt.subplots(1, 4, figsize=(12, 3.6), sharey=True)
    for ax, (name, key, title, theo) in zip(axes, panels):
        rec = tri[keymap.get(key, key)]
        N = rec['N']
        v = np.array([abs(rec['vals'][c]['mean']) for c in combos])
        x = np.arange(8)
        ax.bar(x, np.maximum(v, 1e-12), color=C[name], width=0.62, label='measured |T|')
        floor = 1 / np.sqrt(N)
        ax.axhline(floor, color=ps.EXPECT, ls='--', lw=1.0, label=r'null s.d. $1/\sqrt{N}$')
        if theo is not None:
            tx = [i for i, c in enumerate(combos) if c in theo and abs(theo[c]) > 1e-12]
            ax.plot(tx, [abs(theo[combos[i]]) for i in tx], 'D', ms=5, mfc='none', mec=ps.INK, mew=1.0,
                    label='theory')
        ax.set_yscale('log'); ax.set_ylim(1e-10, 2)
        ax.set_xticks(x); ax.set_xticklabels(combos, fontsize=8)
        ax.set_xlabel('kernel  p q r')
        ax.set_title(title)
        if name == 'RANDU':
            ax.legend(loc='upper right', fontsize=7.5)
    axes[0].set_ylabel(r'$|\langle P_p(x_{n-k})P_q(x_{n-j})P_r(x_n)\rangle|$')
    return ps.save(fig, 'battery_triples')


def fig_cells():
    cells = S['cells']
    fig, axes = plt.subplots(2, 2, figsize=(10, 7.2))
    gens = ['RANDU', 'RANDM', 'MINSTD', 'R250', 'MT19937']
    for ax, suf, title in ((axes[0, 0], 'd2', 'd = 2 (pairs)'), (axes[0, 1], 'd3', 'd = 3 (triples)'),
                           (axes[1, 0], 'd4', 'd = 4 (quadruples)')):
        for name in gens:
            rows = cells[f'{name}_{suf}']['rows']
            k = [r['k'] for r in rows]; z = [r['z'] for r in rows]
            ax.plot(k, z, '-o', ms=3.5, color=C[name], label=name)
        # without-replacement prediction for RANDU (d=2) or RANDM (d=4)
        tname = {'d2': 'RANDU', 'd3': None, 'd4': 'RANDM'}[suf]
        if tname:
            rows = cells[f'{tname}_{suf}']['rows']
            ax.plot([r['k'] for r in rows], [r['zdef'] for r in rows], ls='--', color=ps.EXPECT, lw=1.3,
                    label=r'$-f\sqrt{(C-1)/2}$ (%s)' % tname)
        ax.axhspan(-3, 3, color=ps.GRID, alpha=0.5, lw=0)
        ax.set_xscale('log'); ax.set_yscale('symlog', linthresh=3, linscale=0.6)
        tk = [-1000, -100, -10, -3, 0, 3, 10, 100, 1000, 10000]
        ax.set_yticks(tk); ax.set_yticklabels([f'{t:g}' for t in tk])
        ax.set_xlabel('cells per axis k'); ax.set_ylabel(r'$z = (X^2-\nu)/\sqrt{2\nu}$')
        ax.set_title(title + r', $10^7$ tuples')
        if suf == 'd4':
            ax.legend(fontsize=7.5, loc='upper left', bbox_to_anchor=(0.0, 0.72), ncol=2)
        else:
            ax.legend(fontsize=7.5, loc='lower left', ncol=2 if suf == 'd3' else 1)
    ax = axes[1, 1]
    rows = cells['RANDU_d3']['rows']
    k = np.array([r['k'] for r in rows])
    occ = 1 - np.array([r['empty'] for r in rows]); occ_exp = 1 - np.array([r['empty_exp'] for r in rows])
    kk = np.logspace(np.log10(4), np.log10(256), 100)
    ax.plot(kk, np.minimum(1, 15 / kk), ls='--', color=ps.EXPECT, label='plane theory min(1, 15/k)')
    ax.plot(k, occ_exp, ls=':', color=ps.MUTED, lw=1.4, label=r'ideal: $1-e^{-\lambda}$')
    ax.plot(k, occ, 'o-', color=C['RANDU'], ms=4, label='RANDU measured')
    rows = cells['MT19937_d3']['rows']
    ax.plot(k, 1 - np.array([r['empty'] for r in rows]), 's', color=C['MT19937'], ms=4, mfc='none', label='MT19937 measured')
    ax.set_xscale('log'); ax.set_yscale('log')
    ax.set_xlabel('cells per axis k'); ax.set_ylabel('fraction of occupied cells')
    ax.set_title('Occupied cells, d = 3')
    ax.legend(fontsize=7.5, loc='lower left')
    return ps.save(fig, 'battery_cells')


def fig_toogood():
    sl = S['secondlevel']; K = 2**20; Lb = 2**23
    fig, axes = plt.subplots(1, 3, figsize=(12, 3.9))
    ax = axes[0]
    for name in ('RANDU', 'RANDM', 'MINSTD', 'MT19937', 'R250'):
        cum = np.array(sl[name]['cum']); b = np.arange(1, len(cum) + 1)
        P = PERIOD.get([g for g, n in NAMES.items() if n == name][0], None)
        x = b * Lb / 2**29
        ax.plot(x, cum, 'o-' if name in ('RANDU', 'MINSTD') else '-', ms=2.5, color=C[name], label=name)
    xx = np.linspace(0, 1, 50)
    ax.plot(xx, 1 - xx, ls='--', color=ps.EXPECT, label=r'$1-N/P$ (P = 2^29)')
    ax.plot(xx, 1 - xx * 2**29 / (2**31 - 2), ls=':', color=ps.EXPECT, label=r'$1-N/P$ (P = 2^31-2)')
    ax.set_xlabel(r'prefix length N / 2^29'); ax.set_ylabel(r'$X^2/\nu$, cumulative histogram')
    ax.set_title(r'1D histogram, $2^{20}$ bins'); ax.legend(fontsize=7.5, loc='lower left')
    ax.set_ylim(-0.05, 1.12)

    ax = axes[1]
    for name in ('RANDU', 'RANDM', 'MINSTD', 'MT19937'):
        rows = np.array(S['Q_vs_N'][name]); g = [gg for gg, n in NAMES.items() if n == name][0]
        P = PERIOD.get(g, None)
        x = rows[:, 0] / (P if P else 2**29)
        ax.plot(rows[:, 0], rows[:, 1] / 300, 'o-', ms=3.5, color=C[name], label=name)
    ax.axhline(1, color=ps.EXPECT, ls='--', lw=1.0, label='ideal mean')
    ax.axhspan(1 - 3 * np.sqrt(2 / 300), 1 + 3 * np.sqrt(2 / 300), color=ps.GRID, alpha=0.5, lw=0)
    ax.axvline(2**29, color=ps.MUTED, ls=':', lw=0.8)
    ax.text(2**29, 1.25, 'P(RANDU) ', ha='right', fontsize=8, color=ps.INK2)
    ax.axvline(2**30, color=ps.MUTED, ls=':', lw=0.8)
    ax.text(2**30, 1.25, ' P(MINSTD)/2', ha='left', fontsize=8, color=ps.INK2)
    ax.set_xscale('log'); ax.set_xlabel('number of draws N'); ax.set_ylabel(r'$Q_{300}/300$')
    ax.set_title('Lag correlations, lags 1 to 300'); ax.legend(fontsize=7.5, loc='lower left')
    ax.set_ylim(-0.05, 1.35)

    ax = axes[2]
    bins = np.linspace(0, 1, 21)
    for name, ls in (('RANDU', '-'), ('MINSTD', '-'), ('MT19937', '-')):
        p = np.array(sl[name]['pvals'])
        h, _ = np.histogram(p, bins)
        ax.stairs(h, bins, color=C[name], lw=1.8, label=f"{name} (KS p = {sl[name]['KS_p']:.1e})" if sl[name]['KS_p'] > 0 else f'{name} (KS p < 1e-300)')
    ax.axhline(64 / 20, color=ps.EXPECT, ls='--', lw=1.0, label='uniform')
    ax.set_xlabel('first-level p-value'); ax.set_ylabel('count (64 blocks)')
    ax.set_title(r'Second level: 64 blocks of $2^{23}$')
    ax.legend(fontsize=7.5, loc='upper left')
    return ps.save(fig, 'battery_toogood')


def fig_grid():
    kn = S['knuth']; cells = S['cells']; lag = S['lag']; c3 = S['c3']; sl = S['secondlevel']
    order = ['RANDU', 'RANDM', 'RAND', 'MINSTD', 'R250', 'RANNUM', 'MT19937', 'MIDSQ', 'LOGISTIC']
    cols = ['lag-1 r', r'$Q_{300}$', '3-pt odd', '3-pt even', 'cells d=2', 'cells d=3', 'cells d=4',
            'cells (0,103,250)', '2nd level', 'runs', 'gap', 'perm', 'coupon', 'max-of-8', 'birthday',
            'runs', 'gap', 'perm', 'coupon', 'max-of-8', 'birthday']
    V = np.zeros((len(order), len(cols))); LOW = np.zeros_like(V, bool)
    def cellp(name, suf, k):
        r = [x for x in cells[f'{name}_{suf}']['rows'] if x['k'] == k][0]
        return r['logp'], r['tail'] == 'low'
    for i, name in enumerate(order):
        L = lag[name]
        vals = [(L['logp_z1'], False), (L['logp_Q'], L['tail_Q'] == 'low')]
        if name in c3:
            vals += [(c3[name]['T1']['logp_max'], False), (c3[name]['T2']['logp_max'], False)]
        else:
            vals += [(np.nan, False), (np.nan, False)]
        vals += [cellp(name, 'd2', 1000), cellp(name, 'd3', 100), cellp(name, 'd4', 32), cellp(name, 'lag', 64)]
        s = sl[name]
        vals += [(np.log10(max(s['KS_p'], 1e-300)), s['mean_z'] < 0)]
        for v in (0, 1):
            for t in ('RUNS', 'GAP', 'PERM', 'COUP', 'MAXT', 'BDAY'):
                rr = kn[f'{name}_v{v}'][t]
                vals += [(rr['log10p2'], rr['tail'] == 'low')]
        for j, (lp, low) in enumerate(vals):
            V[i, j] = -lp if np.isfinite(lp) else (300 if not np.isnan(lp) else np.nan)
            LOW[i, j] = low
    Vc = np.clip(V, 0, 12)
    seq = LinearSegmentedColormap.from_list('seq', ['#f4f3ef', '#b9cfee', '#2a78d6', '#123f78'])
    fig, ax = plt.subplots(figsize=(12, 4.6))
    im = ax.imshow(Vc, cmap=seq, vmin=0, vmax=12, aspect='auto')
    ax.grid(False)
    for i in range(V.shape[0]):
        for j in range(V.shape[1]):
            if np.isnan(V[i, j]):
                ax.text(j, i, 'n/a', ha='center', va='center', fontsize=6.5, color=ps.MUTED)
            elif V[i, j] >= 3:
                s = '>12' if V[i, j] > 12 else f'{V[i, j]:.0f}'
                if LOW[i, j]: s += '\nlow'
                ax.text(j, i, s, ha='center', va='center', fontsize=6.8,
                        color='white' if V[i, j] > 6 else ps.INK)
    ax.set_yticks(range(len(order))); ax.set_yticklabels(order)
    ax.set_xticks(range(len(cols))); ax.set_xticklabels(cols, rotation=55, ha='right', fontsize=8.5)
    for xline in (8.5, 14.5):
        ax.axvline(xline, color='white', lw=3)
    ax.text(4, -0.72, 'correlation and serial tests (this family)', ha='center', fontsize=9, color=ps.INK2)
    ax.text(11.5, -0.72, 'Knuth battery, all bits', ha='center', fontsize=9, color=ps.INK2)
    ax.text(17.5, -0.72, r'Knuth battery, frac$(2^8x)$', ha='center', fontsize=9, color=ps.INK2)
    cb = fig.colorbar(im, ax=ax, shrink=0.85, pad=0.01)
    cb.set_label(r'$-\log_{10} p$ (two-sided, clipped at 12)')
    ax.set_title('Summary of all tests: numbers shown where p < 0.001 ("low" = too uniform)', pad=24)
    return ps.save(fig, 'battery_summary_grid')


if __name__ == '__main__':
    for f in (fig_lagcorr, fig_convergence, fig_c3maps, fig_triples, fig_cells, fig_toogood, fig_grid):
        print(f())
