"""Summarize birthday.f90 output: per (generator, t, b, n) block the mean
collision count Y, and a two-sided Poisson p-value for sum(Y) against
Poisson(R*lambda), lambda = n^3/(4k). Writes bday_summary.csv."""
import glob, os
import numpy as np
import pandas as pd
from scipy import stats

D = "/Users/noah/Downloads/ComputationalPhysicsPHYS25000/Project1/P1_Claude/tests/spacing/"
rows = []
for f in sorted(glob.glob(D + "bday/*.txt")):
    tag = os.path.basename(f)[:-4]
    if tag.startswith("check"):
        continue
    d = pd.read_csv(f, sep=r"\s+", comment="#", header=None, names=["t", "b", "n", "lam", "r", "Y"])
    for (t, b, n), g in d.groupby(["t", "b", "n"]):
        lam = g.lam.iloc[0]; R = len(g); Ys = g.Y.sum(); mu = R * lam
        lp_hi = stats.poisson.logsf(Ys - 1, mu) / np.log(10)
        lp_lo = stats.poisson.logcdf(Ys, mu) / np.log(10)
        rows.append(dict(tag=tag, t=t, b=b, n=n, R=R, lam=lam, Ymean=g.Y.mean(), Yvar=g.Y.var(),
                         log10p=min(0.0, np.log10(2) + min(lp_hi, lp_lo))))
pd.DataFrame(rows).to_csv(D + "bday_summary.csv", index=False)
print("wrote", len(rows), "blocks")
