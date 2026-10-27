"""Analysis of ising.f90 output files.
e = <B>/N with B = sum_<ij> s_i s_j ;  C = K^2 (<B^2> - <B>^2)/N.
Errors: batch means for e (batches of nper consecutive measurements), and
delete-one jackknife over batches for C. tau_int from the energy
autocorrelation sums (windowed, Sokal's automatic window c=6) and from
the batch-means / naive variance ratio.
"""
import sys, re, os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
EXACT = {}
for line in open(os.path.join(HERE, "exact_Kc.txt")):
    if line.startswith("#"):
        continue
    a = line.split()
    EXACT[int(a[0])] = (float(a[1]), float(a[2]))


def read(fn):
    hdr = {}
    rows = []
    acf = {}
    acf_n = None
    with open(fn) as fh:
        for line in fh:
            if line.startswith("# acf_nterms"):
                acf_n = int(line.split()[2])
            elif line.startswith("# acf "):
                a = line.split()
                acf[int(a[2])] = int(a[3])
            elif line.startswith("# ndraw"):
                a = line.split()
                hdr["ndraw"] = int(a[2]); hdr["cpu"] = float(a[4])
            elif line.startswith("#"):
                for k, v in re.findall(r"(\w+)=\s*([^\s]+)", line):
                    hdr[k] = v
            else:
                rows.append([int(x) for x in line.split()])
    rows = np.array(rows, dtype=object)
    return hdr, rows, acf, acf_n


def analyze(fn, e_ref=None, C_ref=None, merge=1):
    hdr, rows, acf, acf_n = read(fn)
    N = int(hdr["N"]); K = float(hdr["K"]); L = int(hdr["L"])
    n = rows[:, 1].astype(np.int64)
    sB = rows[:, 2].astype(np.float64); sB2 = rows[:, 3].astype(np.float64)
    sC = rows[:, 4].astype(np.float64)
    # merge consecutive batches
    nb = len(n) // merge
    n = n[:nb * merge].reshape(nb, merge).sum(1)
    sB = sB[:nb * merge].reshape(nb, merge).sum(1)
    sB2 = sB2[:nb * merge].reshape(nb, merge).sum(1)
    sC = sC[:nb * merge].reshape(nb, merge).sum(1)
    ntot = n.sum()
    mB = sB.sum() / ntot; mB2 = sB2.sum() / ntot
    e = mB / N
    eb = sB / n / N
    se = eb.std(ddof=1) / np.sqrt(nb)
    C = K * K * (mB2 - mB * mB) / N
    # jackknife for C
    Cj = np.empty(nb)
    for i in range(nb):
        nn = ntot - n[i]
        b1 = (sB.sum() - sB[i]) / nn; b2 = (sB2.sum() - sB2[i]) / nn
        Cj[i] = K * K * (b2 - b1 * b1) / N
    sCj = np.sqrt((nb - 1) / nb * ((Cj - Cj.mean()) ** 2).sum())
    var_naive = (mB2 - mB * mB) / N**2
    tau_bm = 0.5 * se**2 * ntot / var_naive
    # tau_int from acf sums
    tau_acf = np.nan
    if acf_n:
        mu = mB
        c0 = acf[0] / acf_n - mu * mu
        rho = np.array([(acf[k] / acf_n - mu * mu) / c0 for k in sorted(acf)])
        tau = 0.5
        for W in range(1, len(rho)):
            tau = 0.5 + rho[1:W + 1].sum()
            if W >= 6 * tau:
                break
        tau_acf = tau
    if e_ref is None:
        e_ref, C_ref = EXACT.get(L, (np.nan, np.nan))
    res = dict(file=os.path.basename(fn), alg=hdr.get("alg"), order=hdr.get("order"), src=int(hdr["src"]),
               dec=int(hdr.get("dec", 1)), lags=hdr.get("lags"), L=L, K=K, nmeas=int(ntot), nbatch=nb,
               e=e, se=se, C=C, sC=sCj, de=e - e_ref, zE=(e - e_ref) / se, dC=C - C_ref, zC=(C - C_ref) / sCj,
               tau_bm=tau_bm, tau_acf=tau_acf, csize=sC.sum() / ntot, ndraw=hdr.get("ndraw"),
               draws_per_update=hdr.get("ndraw", 0) / ntot if "ndraw" in hdr else np.nan, cpu=hdr.get("cpu"))
    return res


def fmt(r):
    return (f"{r['file']:28s} {r['alg']:6s} src={r['src']:2d} L={r['L']:3d} n={r['nmeas']:.2e} "
            f"e={r['e']:.6f}({r['se']*1e6:.0f}) de={r['de']:+.2e} z={r['zE']:+7.1f} | C={r['C']:.5f}({r['sC']*1e5:.0f}) "
            f"dC={r['dC']:+.2e} z={r['zC']:+7.1f} | tau_bm={r['tau_bm']:.2f} tau_acf={r['tau_acf']:.2f} "
            f"<c>={r['csize']:.1f} draws/upd={r['draws_per_update']:.1f}")


if __name__ == "__main__":
    for fn in sys.argv[1:]:
        print(fmt(analyze(fn)))
