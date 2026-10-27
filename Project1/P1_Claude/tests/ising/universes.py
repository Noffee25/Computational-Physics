"""Metropolis (typewriter) with RANDU / RANDM / MINSTD: one RNG period of
2^21 sweeps per seed ('universe'); combine universes honestly."""
import sys, os, glob
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from period_check import stats_first, period_report
from analyze import EXACT
e0, C0 = EXACT[16]
out = {}
for g, name, base in [(1, "RANDU", "runs/B_m_randu.dat"), (2, "RANDM", "runs/B_m_randm.dat"), (4, "MINSTD", None)]:
    files = sorted(glob.glob(f"runs/M_g{g}_s*.dat"))
    res = []
    if base:
        res.append(stats_first(base, 128))
    for f in files:
        res.append(stats_first(f, 128))
    res = np.array(res)  # e, se, ze, C, sC, zC
    w = 1 / res[:, 1] ** 2; e = (res[:, 0] * w).sum() / w.sum(); se = 1 / np.sqrt(w.sum())
    wc = 1 / res[:, 4] ** 2; C = (res[:, 3] * wc).sum() / wc.sum(); sC = 1 / np.sqrt(wc.sum())
    chi_e = (((res[:, 0] - e) / res[:, 1]) ** 2).sum(); chi_C = (((res[:, 3] - C) / res[:, 4]) ** 2).sum()
    print(f"{name}: {len(res)} universes, per-universe z_e = {np.round(res[:,2],2)}, z_C = {np.round(res[:,5],2)}")
    print(f"   combined e = {e:.6f} +- {se:.6f}  (z = {(e-e0)/se:+.2f}), C = {C:.5f} +- {sC:.5f} (z = {(C-C0)/sC:+.2f}); "
          f"chi2 between universes e: {chi_e:.1f}/{len(res)-1}, C: {chi_C:.1f}/{len(res)-1}")
    out[name] = (e - e0, se, C - C0, sC, len(res))
    # periodicity of the extra-seed runs (each is exactly one period -> cannot test); report ok
np.save("scratch/universes.npy", out, allow_pickle=True)
