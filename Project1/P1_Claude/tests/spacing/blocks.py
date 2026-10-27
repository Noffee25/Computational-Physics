"""Per-block test results for the spacing scan (used by figures and tables).

For every (configuration, N) block of spacing_scan output:
  * per-replicate exact PIT of the maximum spacing M (exact law, Poisson
    tail approximation below t = -2.5, which is conservative);
  * per-replicate PIT of the minimum spacing m, made resolution-aware:
    the lower-tail probability is evaluated at m + delta with delta the
    spacing of doubles in the window (2^(j-53)), so that an exact tie is
    scored as "m below the double-precision resolution";
  * block p-values: Fisher combination of two-sided per-replicate p,
    plus KS and Anderson-Darling on the PIT values when R >= 40;
    p_block = min(1, k * min(p_1..p_k)) (Bonferroni over the k tests used).
Writes blocks.csv.
"""
import sys, os
sys.path.insert(0, "/Users/noah/Downloads/ComputationalPhysicsPHYS25000/Project1/P1_Claude/tests/spacing")
import numpy as np
import pandas as pd
from scipy import stats
from ms_stats import *
from analyze_scan import load

D = "/Users/noah/Downloads/ComputationalPhysicsPHYS25000/Project1/P1_Claude/tests/spacing"

CONFIGS = [  # tag, label, generator color key, period (draws) or None, grid size G or None
    ("MT", "MT19937", "MT19937"),
    ("RANNUM", "RANNUM", "RANNUM"),
    ("R250", "R250", "R250"),
    ("MINSTD", "MINSTD", "MINSTD"),
    ("RANDU", "RANDU", "RANDU"),
    ("RANDM", "RANDM", "RANDM"),
    ("RAND", "RAND", "RAND"),
    ("MTsingle", "MT19937 (REAL*4)", "MT19937"),
    ("RANDUsingle", "RANDU (REAL*4)", "RANDU"),
    ("LOGISTIC", "LOGISTIC", "LOGISTIC"),
    ("RANDU65536", "RANDU seed 65536", "RANDU"),
    ("RANDM65536", "RANDM seed 65536", "RANDM"),
    ("MIDSQ12345678", "MIDSQ seed 12345678", "MIDSQ"),
    ("MIDSQ777", "MIDSQ seed 777", "MIDSQ"),
]


def fisher_log10(lp):
    lp = np.minimum(lp, 0.0)
    X = -2.0 * np.sum(lp)
    if not np.isfinite(X):
        return -np.inf
    return stats.chi2.logsf(X, 2 * len(lp)) / LN10


def combine(pf, pks, pad):
    ps = [p for p in (pf, pks, pad) if p is not None and not np.isnan(p)]
    k = len(ps)
    m = min(ps)
    return min(0.0, m + np.log10(k))


def replicate_logs(b):
    n = b.n.values.astype(float)
    lF, lS = maxsp_cdf_sf(n, b.M.values)
    delta = 2.0 ** (b.j.values.astype(float) - 53.0)
    lFm, _ = minsp_cdf_sf(n, b.m.values + delta)
    _, lSm = minsp_cdf_sf(n, b.m.values)
    return n, lF, lS, lFm, lSm


def block_table():
    rows = []
    for tag, label, ckey in CONFIGS:
        df = load(tag)
        df = df[df.flag == 0]
        for N, b in df.groupby("N", sort=True):
            n, lF, lS, lFm, lSm = replicate_logs(b)
            R = len(b)
            row = dict(tag=tag, label=label, N=N, R=R, j=int(b.j.iloc[0]), n_mean=n.mean())
            T = gumbel_T(n, b.M.values)
            row["T_mean"] = T.mean()
            row["T_sd"] = T.std(ddof=1) if R > 1 else np.nan
            for name, a, s in (("max", lF, lS), ("min", lFm, lSm)):
                lp = np.log(2.0) + np.minimum(a, s)
                pf = fisher_log10(lp)
                pks = pad = None
                if R >= 40:
                    u = pit_from_logs(a, s)
                    Dks, p1 = ks_test(u)
                    A2, p2 = ad_test(u)
                    pks = np.log10(max(p1, 1e-320))
                    pad = np.log10(max(p2, 1e-320))
                    row[name + "_KS_D"] = Dks
                    row[name + "_AD"] = A2
                row[name + "_log10p_fisher"] = pf
                row[name + "_log10p_ks"] = pks if pks is not None else np.nan
                row[name + "_log10p_ad"] = pad if pad is not None else np.nan
                row[name + "_log10p"] = combine(pf, pks, pad)
                Z, lpz = stouffer(a, s)
                row[name + "_Z"] = Z
            m = b.m.values
            row["frac_m0"] = float(np.mean(m == 0))
            row["mean_nnm"] = float(np.mean(n * (n + 1) * m))
            row["dup_sum"] = int(b.dup.sum())
            row["dup_mean"] = float(b.dup.mean())
            row["pairs_equiv"] = float(np.sum(n * (n - 1) / 2.0) * 2.0 ** b.j.iloc[0])
            rows.append(row)
    return pd.DataFrame(rows)


if __name__ == "__main__":
    B = block_table()
    B.to_csv(os.path.join(D, "blocks.csv"), index=False)
    pd.set_option("display.width", 250); pd.set_option("display.max_rows", 500)
    print(B[["label", "N", "R", "T_mean", "T_sd", "max_Z", "max_log10p_fisher", "max_log10p_ks", "max_log10p_ad",
             "max_log10p", "min_log10p", "frac_m0", "mean_nnm", "dup_mean"]].to_string(float_format=lambda v: "%.3g" % v))
