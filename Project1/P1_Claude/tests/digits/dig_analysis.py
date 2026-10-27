# chi-square analysis of the digit count tables written by digscan
import numpy as np, glob, re, os, json
from scipy import stats
D = "/Users/noah/Downloads/ComputationalPhysicsPHYS25000/Project1/P1_Claude/tests/digits/out/"
NAMES = ['RANDU', 'RANDM', 'RAND', 'MINSTD', 'R250', 'RANNUM', 'MT19937', 'MIDSQ', 'LOGISTIC']

def load(fn):
    tab = {}
    for line in open(fn):
        t = line.split()
        if t[0] in ('STO', 'INT', 'NAI'):
            tab.setdefault(t[0], {})[int(t[1])] = np.array([int(v) for v in t[2:]], dtype=float)
        elif t[0] in ('MISNAI', 'MISINT'):
            tab.setdefault(t[0], {})[int(t[1])] = int(t[2])
        elif t[0] == '#':
            tab['N'] = int(t[t.index('N') + 1])
    return tab

def chi2(c):
    n = c.sum(); e = n / 10.0
    return ((c - e) ** 2 / e).sum()

def signif(x2, dof=9):
    """signed significance: +(-log10 p_right) if chi2 above median, -(-log10 p_left) below"""
    pr = stats.chi2.sf(x2, dof); pl = stats.chi2.cdf(x2, dof)
    if pr < pl:
        return -np.log10(max(pr, 1e-300)), pr, pl
    return np.log10(max(pl, 1e-300)), pr, pl

def table(fn, kind, pmax):
    tab = load(fn)
    out = []
    for p in range(1, pmax + 1):
        if p not in tab[kind]:
            out.append((p, np.nan, np.nan, np.nan, np.nan)); continue
        x2 = chi2(tab[kind][p]); s, pr, pl = signif(x2)
        out.append((p, x2, s, pr, pl))
    return out, tab

if __name__ == "__main__":
    for prec in (0, 1):
        for kind, pm in (('INT', 25), ('STO', 25), ('NAI', 18 if prec == 0 else 12)):
            if prec == 1 and kind == 'INT':
                continue
            print("\n=== prec", prec, kind, " (signed significance, chi2)")
            for g in range(1, 10):
                fn = D + f"dig_g{g}_p{prec}.txt"
                out, tab = table(fn, kind, pm)
                print(f"{NAMES[g-1]:9s}", " ".join(f"{o[2]:6.1f}" for o in out))
