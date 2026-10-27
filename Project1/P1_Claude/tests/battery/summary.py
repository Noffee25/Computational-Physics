"""summary.py : turn the raw Fortran outputs of the battery family into the
numbers quoted in the report (written to out/summary.json and printed)."""
import sys
sys.path.insert(0, "/Users/noah/Downloads/ComputationalPhysicsPHYS25000/Project1/P1_Claude/tests/battery")
import json, glob, re
import numpy as np
from scipy import stats
import theory
from knuth_table import analyse as kn_analyse, TESTS as KN_TESTS
from readc3 import readc3

BASE = "/Users/noah/Downloads/ComputationalPhysicsPHYS25000/Project1/P1_Claude/tests/battery/"
NAMES = {1: 'RANDU', 2: 'RANDM', 3: 'RAND', 4: 'MINSTD', 5: 'R250', 6: 'RANNUM',
         7: 'MT19937', 8: 'MIDSQ', 9: 'LOGISTIC'}
PERIOD = {1: 2**29, 2: 2**29, 3: 2**31 - 2, 4: 2**31 - 2}
MULT = {1: 65539, 2: 6539, 4: 16807}
LN10 = np.log(10)


def two_sided_from_z(z):
    return (stats.norm.logsf(abs(z)) + np.log(2)) / LN10


def chi2_two_sided(X, df):
    """log10 of the two-sided p-value and the tail ('high' or 'low')."""
    up = stats.chi2.logsf(X, df); lo = stats.chi2.logcdf(X, df)
    if up < lo:
        return (min(0.0, up + np.log(2))) / LN10, 'high'
    return (min(0.0, lo + np.log(2))) / LN10, 'low'


def load_lag(fn):
    hdr = {}
    for l in open(fn):
        if l.startswith('# N '): hdr['N'] = int(l.split()[2])
    d = np.loadtxt(fn, comments='#')
    return hdr['N'], d


S = {}

# ------------------------------------------------------------------ lags
lag = {}
for g in range(1, 10):
    N, d = load_lag(BASE + f'out/lag/lag_{g}.txt')
    k = d[:, 0]; cnt = N - k
    z = d[:, 1:] / np.sqrt(cnt)[:, None]
    r = d[:, 1:] / cnt[:, None]
    Q = float((z[:, 0]**2).sum())
    lq, tail = chi2_two_sided(Q, 300)
    imax = int(np.argmax(np.abs(z[:, 0])))
    lag[NAMES[g]] = dict(N=N, r1=float(r[0, 0]), z1=float(z[0, 0]), Q300=Q, logp_Q=lq, tail_Q=tail,
                         maxabsz=float(abs(z[imax, 0])), kmax=int(k[imax]),
                         r22_1=float(r[0, 3]), r12_1=float(r[0, 1]), r21_1=float(r[0, 2]),
                         logp_z1=float(two_sided_from_z(z[0, 0])))
# logistic before its floating point cycle
N, d = load_lag(BASE + 'out/lagN/lag_9_23.txt')
z = d[:, 1] / np.sqrt(N - d[:, 0]); r = d[:, 1:] / (N - d[:, 0])[:, None]
lag['LOGISTIC_pre'] = dict(N=N, Q300=float((z**2).sum()), r1=float(r[0, 0]), z1=float(z[0]),
                           r22_1=float(r[0, 3]), r12_1=float(r[0, 1]), r21_1=float(r[0, 2]))
# Box-Pierce versus N
qN = {}
for g in (1, 2, 4, 7):
    rows = []
    for lg in range(24, 31):
        fn = BASE + (f'out/lagN/lag_{g}_{lg}.txt' if lg != 29 else f'out/lag/lag_{g}.txt')
        try:
            N, d = load_lag(fn)
        except Exception:
            continue
        z = d[:, 1] / np.sqrt(N - d[:, 0])
        rows.append((N, float((z**2).sum())))
    qN[NAMES[g]] = rows
S['lag'] = lag; S['Q_vs_N'] = qN

# exact full period values and Dedekind theory
exact = {}
for g, fn in ((1, 'out/exact_randu.txt'), (2, 'out/exact_randm.txt')):
    rows = {}
    for l in open(BASE + fn):
        m = re.match(r'lag\s+(\d+)\s+a\^k mod m =\s+(\d+)\s+r_k =\s+(\S+)', l)
        if m:
            kk, b, rk = int(m.group(1)), int(m.group(2)), float(m.group(3))
            rows[kk] = dict(b=b, exact=rk, dedekind=theory.serial_dedekind(b, 2**31),
                            z_at_P=rk * np.sqrt(2**29))
    exact[NAMES[g]] = rows
exact['MINSTD'] = {1: dict(b=16807, exact=None, dedekind=theory.serial_dedekind(16807, 2**31 - 1),
                           z_at_P=theory.serial_dedekind(16807, 2**31 - 1) * np.sqrt(2**31 - 2))}
# Dedekind prediction for all 300 lags, compared with the full period run
ded = {}
for g in (1, 2):
    a = MULT[g]; m = 2**31
    pred = np.array([theory.serial_dedekind(pow(a, kk, m), m) for kk in range(1, 301)])
    N, d = load_lag(BASE + f'out/lag/lag_{g}.txt')
    meas = d[:, 1] / (N - d[:, 0])
    ded[NAMES[g]] = dict(pred=pred.tolist(), meas=meas.tolist(),
                         maxdiff=float(np.max(np.abs(pred - meas))),
                         maxabs=float(np.max(np.abs(meas))))
S['exact'] = exact; S['dedekind_all'] = {k: dict(maxdiff=v['maxdiff'], maxabs=v['maxabs']) for k, v in ded.items()}
json.dump(ded, open(BASE + 'out/dedekind_lags.json', 'w'))

# ------------------------------------------------------------------ C3 maps
c3 = {}
for g in (1, 2, 3, 4, 5, 6, 7, 9):
    K, N, T1, T2, Z1, Z2 = readc3(BASE + f'out/c3/c3_{g}.bin')
    iu = np.triu_indices(K, 1)
    M = len(iu[0])
    out = {}
    for name, Z, T in (('T1', Z1, T1), ('T2', Z2, T2)):
        z = Z[iu]
        o = np.argsort(-np.abs(z))[:3]
        zmax = float(np.abs(z).max())
        # Sidak-corrected two-sided p of the maximum
        lp1 = (stats.norm.logsf(zmax) + np.log(2)) / LN10
        p1 = 10**lp1
        lp = np.log10(-np.expm1(M * np.log1p(-p1))) if p1 > 1e-300 and M * p1 < 1e3 else min(0, lp1 + np.log10(M))
        # null check: fraction of |z|>3 and KS of z against N(0,1) excluding top 5
        zz = np.sort(z)[3:-3]
        out[name] = dict(zmax=zmax, logp_max=float(min(0.0, lp)),
                         top=[(int(iu[0][t] + 1), int(iu[1][t] + 1), float(z[t]), float(T[iu][t])) for t in o],
                         sd=float(zz.std()), frac_gt3=float(np.mean(np.abs(z) > 3)))
    c3[NAMES[g]] = dict(N=int(N), K=int(K), pairs=M, **out)
S['c3'] = c3

# ------------------------------------------------------------------ triples
tri = {}
pred = {
    'RANDU_1_2': (9, -6, 1), 'RANDU_2_4': (81, -18, 1)}
for fn in sorted(glob.glob(BASE + 'out/tri/tri_*.txt')):
    key = fn.split('/')[-1][4:-4]
    lines = open(fn).read().split('\n')
    Nn = int(re.search(r'N\s+(\d+)', lines[0]).group(1))
    vals = {}
    for l in lines[2:]:
        if l.strip():
            p, q, r, mu, m2 = l.split()
            vals[p + q + r] = dict(mean=float(mu), sd=float(np.sqrt(max(float(m2) - float(mu)**2, 0))),
                                   z=float(mu) * np.sqrt(Nn))
    tri[key] = dict(N=Nn, vals=vals)
# theory
th = {}
for nm, c in (('1_1_2', (9, -6, 1)), ('1_2_4', (81, -18, 1))):
    th['RANDU_' + nm] = {f'{p}{q}{r}': theory.triple_fourier(*c, (p, q, r))
                         for p in (1, 2) for q in (1, 2) for r in (1, 2)}
th['R250_103_250'] = {'111': -12**1.5 / 56, '222': 0.6845106053570786}
S['triples'] = tri; S['triples_theory'] = th

# ------------------------------------------------------------------ cells
cells = {}
for g in range(1, 10):
    for suf in ('d2', 'd3', 'd4', 'lag'):
        fn = BASE + f'out/cells/c_{g}_{suf}.txt'
        L = open(fn).read().split('\n')
        Nt = int(re.search(r'Nt\s+(\d+)', L[0]).group(1))
        d = int(re.search(r' d\s+(\d+)', L[0]).group(1))
        rows = []
        for l in L[2:]:
            if l.strip():
                k, C, s2, ne, mx = map(int, l.split())
                lam = Nt / C
                X2 = s2 / lam - Nt
                z = (X2 - (C - 1)) / np.sqrt(2 * (C - 1))
                lp, tail = chi2_two_sided(X2, C - 1)
                # predicted 'without replacement' deficit for the LCGs
                zdef = None
                if g in PERIOD:
                    P = PERIOD[g]
                    Ppop = P // np.gcd(P, d) if suf != 'lag' else P
                    zdef = -(Nt / Ppop) * (C - 1) / np.sqrt(2 * (C - 1))
                rows.append(dict(k=k, C=C, X2=X2, z=float(z), logp=float(lp), tail=tail,
                                 empty=ne / C, empty_exp=float(np.exp(-lam)), maxc=mx, zdef=zdef))
        cells[f'{NAMES[g]}_{suf}'] = dict(Nt=Nt, d=d, rows=rows)
S['cells'] = cells

# ------------------------------------------------------------------ Knuth battery
kn = {}
for v in (0, 1):
    for g in range(1, 10):
        r = kn_analyse(BASE + f'out/kn/kn_{g}_v{v}.txt')
        row = {}
        for t in KN_TESTS:
            row[t] = dict(stat=float(r[t][0]), log10p=float(r[t][1]))
            if t != 'LAG1':
                # two-sided version for the summary grid
                df = {'RUNS': 5, 'GAP': 40, 'PERM': 23, 'COUP': 32, 'MAXT': 99, 'BDAY': 7}[t]
                lp, tail = chi2_two_sided(r[t][0], df) if np.isfinite(r[t][0]) else (-np.inf, 'high')
                row[t]['log10p2'] = float(lp); row[t]['tail'] = tail
            else:
                row[t]['log10p2'] = float(r[t][1]); row[t]['tail'] = 'high'
        kn[f'{NAMES[g]}_v{v}'] = row
S['knuth'] = kn
# pooled birthday spacings over five seeds
from knuth_table import parse as kn_parse, PR
pool = {}
for g in (1, 5, 7):
    files = [BASE + f'out/kn/kn_{g}_v0.txt'] + sorted(glob.glob(BASE + f'out/val2/kn_{g}_*.txt'))
    tot = np.zeros(8)
    for f in files: tot += kn_parse(f)['BDAY']
    e = tot.sum() * PR['BDAY']; X = float(((tot - e)**2 / e).sum())
    pool[NAMES[g]] = dict(R=int(tot.sum()), X2=X, log10p=float(stats.chi2.logsf(X, 7) / LN10),
                          meanJ=float((tot * np.arange(8)).sum() / tot.sum()),
                          nullJ=float((PR['BDAY'] * np.arange(8)).sum()))
S['bday_pooled'] = pool
S['bday_null'] = json.load(open(BASE + 'out/bday_null.json'))
# validation: distribution of the 8+8 MT/RANNUM battery p-values
vals = []
for f in sorted(glob.glob(BASE + 'out/val/kn_*.txt')):
    r = kn_analyse(f)
    vals += [10**r[t][1] for t in KN_TESTS]
vals = np.array(vals)
S['knuth_validation'] = dict(n=len(vals), ks_p=float(stats.kstest(vals, 'uniform').pvalue),
                             min_p=float(vals.min()), frac_below_0p01=float(np.mean(vals < 0.01)))

# ------------------------------------------------------------------ second level
K = 2**20; Lb = 2**23
sl = {}
for g in range(1, 10):
    d = np.loadtxt(BASE + f'out/sl/sl_{g}.txt', comments='#')
    zf = (d[:, 1] - (K - 1)) / np.sqrt(2 * (K - 1))
    zc = (d[:, 2] - 1023) / np.sqrt(2 * 1023)
    pf = stats.chi2.sf(d[:, 1], K - 1); pc = stats.chi2.sf(d[:, 2], 1023)
    ksf = stats.kstest(pf, 'uniform'); ksc = stats.kstest(pc, 'uniform')
    pred = -(K - 1) * (Lb / PERIOD[g]) / np.sqrt(2 * (K - 1)) if g in PERIOD else 0.0
    sl[NAMES[g]] = dict(mean_z=float(zf.mean()), sd_z=float(zf.std(ddof=1)), pred_z=float(pred),
                        KS_D=float(ksf.statistic), KS_log10p=float(np.log10(max(ksf.pvalue, 1e-300))),
                        KS_p=float(ksf.pvalue),
                        coarse_mean_z=float(zc.mean()), coarse_KS_D=float(ksc.statistic), coarse_KS_p=float(ksc.pvalue),
                        cum=(d[:, 3] / (K - 1)).tolist(), pvals=pf.tolist(), min_p=float(pf.min()), max_p=float(pf.max()))
S['secondlevel'] = sl

# deep bins of the student's histogram
deep = {}
for g in (1, 2, 7):
    L = open(BASE + f'out/deep/deep_{g}.txt').read().split('\n')
    m = re.search(r'mean sd min max\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)', L[1])
    deep[NAMES[g]] = dict(mean=float(m.group(1)), sd=float(m.group(2)), min=int(m.group(3)), max=int(m.group(4)))
P = 2**29; Nd = 10**9; Kb = 10**7
f = Nd / P - 1; EA = P / Kb; pA = EA - np.floor(EA)
varA = pA * (1 - pA)
deep['prediction'] = dict(f=f, EA=EA, sd=float(np.sqrt(varA * (1 + f)**2 + EA * f * (1 - f))), poisson_sd=10.0)
S['deepbins'] = deep

json.dump(S, open(BASE + 'out/summary.json', 'w'), indent=1, default=float)

# ------------------------------------------------------------------ print key numbers
if __name__ == '__main__':
    print('LAGS'); [print(' ', k, {kk: (round(vv, 7) if isinstance(vv, float) else vv) for kk, vv in v.items()}) for k, v in lag.items()]
    print('EXACT', json.dumps(exact, indent=0)[:1500])
    print('DEDEKIND all lags', S['dedekind_all'])
    print('C3'); [print(' ', k, v['T1']['zmax'], v['T1']['top'][0], '|', v['T2']['zmax'], v['T2']['top'][:2], 'sd', v['T1']['sd'], v['T2']['sd']) for k, v in c3.items()]
    print('BDAY pooled', pool)
    print('KNUTH validation', S['knuth_validation'])
    print('SECOND LEVEL'); [print(' ', k, {kk: vv for kk, vv in v.items() if kk not in ('cum', 'pvals')}) for k, v in sl.items()]
    print('DEEP', deep)
