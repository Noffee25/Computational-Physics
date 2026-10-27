"""Exact-periodicity test of a Markov chain driven by a periodic generator.
If batch b and b+P have identical integer sums for all b, the chain has
become exactly periodic (common-random-number coupling of X_t and X_{t+P})."""
import sys, os
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from analyze import read, EXACT

def rows_of(fn):
    hdr, rows, _, _ = read(fn)
    return hdr, rows

def period_report(fn, P):
    hdr, rows = rows_of(fn)
    nb = len(rows)
    same = np.array([all(rows[b, 2:] == rows[b + P, 2:]) for b in range(nb - P)])
    return same.sum(), len(same), (np.argmax(same) if same.any() else None)

def stats_first(fn, nbat):
    hdr, rows = rows_of(fn)
    N = int(hdr["N"]); K = float(hdr["K"])
    r = rows[:nbat]
    n = r[:, 1].astype(float); sB = r[:, 2].astype(float); sB2 = r[:, 3].astype(float)
    eb = sB / n / N
    e = sB.sum() / n.sum() / N; se = eb.std(ddof=1) / np.sqrt(nbat)
    mB = sB.sum() / n.sum(); mB2 = sB2.sum() / n.sum()
    C = K * K * (mB2 - mB * mB) / N
    Cj = []
    for i in range(nbat):
        nn = n.sum() - n[i]; b1 = (sB.sum() - sB[i]) / nn; b2 = (sB2.sum() - sB2[i]) / nn
        Cj.append(K * K * (b2 - b1 * b1) / N)
    Cj = np.array(Cj); sC = np.sqrt((nbat - 1) / nbat * ((Cj - Cj.mean()) ** 2).sum())
    e0, C0 = EXACT[16]
    return e, se, (e - e0) / se, C, sC, (C - C0) / sC

if __name__ == "__main__":
    for fn in ["runs/B_m_randu.dat", "runs/B_m_randm.dat"]:
        for P in (64, 128):
            print(fn, "P =", P, "identical (b, b+P) pairs, total, first:", period_report(fn, P))
        e, se, ze, C, sC, zC = stats_first(fn, 128)
        print(f"   first period only (128 batches = 2^21 sweeps): e={e:.6f}({se*1e6:.0f}) z={ze:+.2f}   C={C:.5f}({sC*1e5:.0f}) z={zC:+.2f}")
        hdr, rows = rows_of(fn)
        e, se, ze, C, sC, zC = stats_first(fn, len(rows))
        print(f"   all {len(rows)} batches (naive):              e={e:.6f}({se*1e6:.0f}) z={ze:+.2f}   C={C:.5f}({sC*1e5:.0f}) z={zC:+.2f}")
    # Wolff + RANDU: any duplicated batches?
    for fn in ["runs/A_w_randu.dat", "runs/A_w_randm.dat"]:
        hdr, rows = rows_of(fn)
        keys = [tuple(r[2:]) for r in rows]
        print(fn, "duplicate batch rows:", len(keys) - len(set(keys)))
