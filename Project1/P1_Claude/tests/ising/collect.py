"""Analyze every run in runs/ and write results.csv (one row per run)."""
import glob, os, sys
import numpy as np
import pandas as pd
HERE = "/Users/noah/Downloads/ComputationalPhysicsPHYS25000/Project1/P1_Claude/tests/ising"
sys.path.insert(0, HERE)
from analyze import analyze, read

KC = 0.44068679350977151
EXK = {}
for line in open(os.path.join(HERE, "exact_K.txt")):
    if line.startswith("#"):
        continue
    a = line.split()
    EXK[round(float(a[1]), 6)] = (float(a[2]), float(a[3]))

rows = []
for fn in sorted(glob.glob(os.path.join(HERE, "runs", "[A-H]_*.dat"))):
    hdr, _, _, _ = read(fn)
    if "cpu" not in hdr:
        continue
    K = float(hdr["K"])
    if abs(K - KC) > 1e-9:
        e_ref, C_ref = EXK[round(K, 6)]
        r = analyze(fn, e_ref, C_ref)
    else:
        r = analyze(fn)
    # binning check: error of e from batches merged by 1, 4, 10
    for m in (4, 10):
        rm = analyze(fn, r["e"] - r["de"], r["C"] - r["dC"], merge=m)
        r[f"se_m{m}"] = rm["se"]; r[f"sC_m{m}"] = rm["sC"]
    r["name"] = os.path.basename(fn)[:-4]
    r["seed"] = int(hdr.get("seed", 777))
    rows.append(r)
df = pd.DataFrame(rows)
df.to_csv(os.path.join(HERE, "results.csv"), index=False)
pd.set_option("display.width", 250); pd.set_option("display.max_rows", 200)
cols = ["name", "nmeas", "e", "se", "se_m10", "de", "zE", "C", "sC", "sC_m10", "dC", "zC", "tau_acf", "csize", "draws_per_update", "cpu"]
print(df[cols].to_string(index=False, float_format=lambda x: f"{x:.6g}"))
