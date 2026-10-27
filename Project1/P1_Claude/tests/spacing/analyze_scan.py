"""Summarize spacing_scan output: max-spacing and min-spacing goodness of fit
per (configuration, N) block, and duplicate counts.

Writes scan_summary.csv in tests/spacing.
"""
import sys, os, glob
sys.path.insert(0, "/Users/noah/Downloads/ComputationalPhysicsPHYS25000/Project1/P1_Claude/tests/spacing")
import numpy as np
import pandas as pd
from ms_stats import *

D = "/Users/noah/Downloads/ComputationalPhysicsPHYS25000/Project1/P1_Claude/tests/spacing"


def load(tag):
    f = os.path.join(D, "scan", tag + ".txt")
    a = np.loadtxt(f, comments="#")
    df = pd.DataFrame(a, columns=["N", "r", "j", "n", "M", "m", "dup", "flag"])
    for c in ["N", "r", "j", "n", "dup", "flag"]:
        df[c] = df[c].astype(np.int64)
    return df


def block_stats(b):
    out = {}
    R = len(b)
    out["R"] = R
    good = b[b.flag == 0]
    out["overflow"] = int((b.flag == 2).sum())
    n = good.n.values.astype(float)
    M = good.M.values
    m = good.m.values
    out["n_mean"] = n.mean() if len(n) else np.nan
    # max spacing
    lF, lS = maxsp_cdf_sf(n, M)
    okk = np.isfinite(lF) & (n >= 1)
    lF, lS = lF[okk], lS[okk]
    T = gumbel_T(n[okk], M[okk])
    out["T_mean"] = T.mean() if len(T) else np.nan
    out["T_sd"] = T.std(ddof=1) if len(T) > 1 else np.nan
    if len(T) >= 1:
        u = pit_from_logs(lF, lS)
        Z, lp = stouffer(lF, lS)
        out["max_Z"] = Z
        out["max_log10p_Z"] = lp
        if len(T) >= 5:
            Dks, pks = ks_test(u)
            A2, pad = ad_test(u)
        else:
            Dks = pks = A2 = pad = np.nan
        out["max_KS_D"] = Dks; out["max_KS_p"] = pks
        out["max_AD"] = A2; out["max_AD_p"] = pad
        # single-replicate extreme: two-sided log10 p of the most extreme replicate
        l2 = (np.log(2) + np.minimum(lF, lS)) / LN10
        out["max_min_log10p_single"] = l2.min()
        out["max_median_log10p_single"] = np.median(l2)
    # min spacing
    lFm, lSm = minsp_cdf_sf(n, m)
    um = pit_from_logs(lFm, lSm)
    out["frac_m0"] = float(np.mean(m == 0)) if len(m) else np.nan
    out["mean_nnm"] = float(np.mean(n * (n + 1) * m)) if len(m) else np.nan
    if len(m) >= 5:
        Dks, pks = ks_test(um)
        out["min_KS_D"] = Dks; out["min_KS_p"] = pks
        A2, pad = ad_test(um)
        out["min_AD"] = A2; out["min_AD_p"] = pad
    # duplicates
    out["dup_sum"] = int(good.dup.sum())
    out["dup_mean"] = good.dup.mean()
    out["n2_sum"] = float(np.sum(n * (n - 1) / 2.0))
    return out


def summarize(tags):
    rows = []
    for tag in tags:
        df = load(tag)
        for N, b in df.groupby("N", sort=True):
            s = block_stats(b)
            s["tag"] = tag
            s["N"] = N
            s["j"] = int(b.j.iloc[0])
            rows.append(s)
    return pd.DataFrame(rows)


if __name__ == "__main__":
    tags = sorted(os.path.basename(f)[:-4] for f in glob.glob(os.path.join(D, "scan", "*.txt")))
    S = summarize(tags)
    S.to_csv(os.path.join(D, "scan_summary.csv"), index=False)
    pd.set_option("display.width", 250)
    pd.set_option("display.max_rows", 500)
    cols = ["tag", "N", "R", "n_mean", "T_mean", "T_sd", "max_Z", "max_log10p_Z", "max_KS_p", "max_AD_p",
            "max_min_log10p_single", "frac_m0", "mean_nnm", "min_KS_D", "min_KS_p", "dup_mean"]
    print(S[cols].to_string(float_format=lambda v: "%.4g" % v))
