"""Loaders and theory for the 'hist' family (deep histogram tests).

Run scripts that import this with /usr/local/bin/python3 -I.
"""
import os
import sys
import glob
import numpy as np
import pandas as pd
from scipy import stats
from scipy.special import gammaln

ROOT = "/Users/noah/Downloads/ComputationalPhysicsPHYS25000/Project1/P1_Claude"
H = os.path.join(ROOT, "tests", "hist")
sys.path.insert(0, ROOT)

SUMCOLS = ["tag", "gen", "prec", "seed", "B", "N", "nempty", "avg", "var", "D",
           "chi2", "z", "n4lo", "n4hi", "e4lo", "e4hi", "minc", "maxc", "secs"]

# periods (number of draws in one full cycle of the state)
P_RANDU = 2**29
P_MINSTD = 2**31 - 2
P_65536 = 8192
# logistic map in double precision from seed 777 (measured by logis_cycle.f90)
LOG_MU = 14324259
LOG_L = 5638349


def load_sum(pattern):
    fs = sorted(glob.glob(os.path.join(H, pattern)))
    df = pd.concat([pd.read_csv(f, sep=r"\s+", names=SUMCOLS) for f in fs], ignore_index=True)
    df["series"] = df["tag"].str.replace(r"_e\d$", "", regex=True)
    df["depth"] = np.round(np.log10(df["B"])).astype(int)
    df["lam"] = df["N"] / df["B"]
    return df


def load_coc(pattern):
    fs = sorted(glob.glob(os.path.join(H, pattern)))
    return pd.concat([pd.read_csv(f, sep=r"\s+", names=["tag", "k", "nb"]) for f in fs],
                     ignore_index=True)


def load_period():
    coc = pd.read_csv(os.path.join(H, "per_coc.out"), sep=r"\s+", names=["tag", "depth", "c", "m"])
    sm = pd.read_csv(os.path.join(H, "per_sum.out"), sep=r"\s+",
                     names=["tag", "gen", "prec", "seed", "B", "P", "tot", "cmin", "cmax",
                            "nret", "okrep", "secs"])
    return coc, sm


def null_sd(B, N):
    """Exact multinomial SD of D = chi2/(B-1)."""
    return np.sqrt(2.0 * (1.0 - 1.0 / N) / (B - 1.0))


# ---------------------------------------------------------------- periodic
def pred_periodic(c, m, P, N, B):
    """Expected D, its SD, expected empty fraction for a generator whose
    one-period histogram has m[i] bins holding c[i] reachable values.
    N = qP + r draws: n_b = q c_b + Hypergeometric(P, c_b, r)."""
    c = np.asarray(c, float)
    m = np.asarray(m, float)
    q, r = divmod(int(N), int(P))
    lam = N / B
    mu = c * N / P
    v = r * (c / P) * (1 - c / P) * (P - r) / (P - 1.0)
    a = mu - lam
    S = np.sum(m * (v + a**2))
    V = np.sum(m * (2 * v**2 + 4 * v * a**2))
    D = S / (lam * (B - 1))
    sd = np.sqrt(V) / (lam * (B - 1))
    if q >= 1:
        pe = np.sum(m[c == 0]) / B
    else:
        # P(no hit) = C(P-c, r)/C(P, r)
        lp = np.where(c <= P - r,
                      gammaln(P - c + 1) - gammaln(P - c - r + 1) - gammaln(P + 1) + gammaln(P - r + 1),
                      -np.inf)
        pe = np.sum(m * np.exp(lp)) / B
    return D, sd, pe


def coc_periodic(c, m, P, N, B, kmax):
    """Predicted count-of-counts pmf (fraction of bins with k hits)."""
    q, r = divmod(int(N), int(P))
    f = r / P
    k = np.arange(kmax + 1)
    out = np.zeros(kmax + 1)
    for ci, mi in zip(c, m):
        ci = int(ci)
        if ci == 0:
            out[0] += mi
            continue
        out += mi * stats.binom.pmf(k - q * ci, ci, f)
    return k, out / B


# ------------------------------------------------- rounded continuous (REAL*4)
def load_fw():
    sm = pd.read_csv(os.path.join(H, "fw_sum.out"), sep=r"\s+",
                     names=["depth", "B", "sdel2", "nzero", "sumw"])
    w = pd.read_csv(os.path.join(H, "fw_w.out"), sep=r"\s+", names=["depth", "w", "nb"])
    return sm, w


def pred_weights(w, nb, N, B):
    """Bins with relative weight w (mean 1): n_b ~ Poisson(lam w)."""
    w = np.asarray(w, float)
    nb = np.asarray(nb, float)
    lam = N / B
    mm = lam * w
    a = mm - lam
    S = np.sum(nb * (mm * (1 - w / B) + a**2))
    V = np.sum(nb * (mm + 2 * mm**2 + 4 * a**2 * mm + 4 * a * mm))
    D = S / (lam * (B - 1))
    sd = np.sqrt(V) / (lam * (B - 1))
    pe = np.sum(nb * np.exp(-mm)) / B
    return D, sd, pe


def coc_weights(w, nb, N, B, kmax):
    lam = N / B
    k = np.arange(kmax + 1)
    out = np.zeros(kmax + 1)
    for wi, ni in zip(w, nb):
        out += ni * stats.poisson.pmf(k, lam * wi)
    return k, out / B


# ------------------------------------------------------ finite word lattice
def pred_lattice(bits, N, B):
    """Outputs k/2^bits, all values equally likely, drawn with replacement
    (period >> N).  Bins hold floor(g) or floor(g)+1 lattice points."""
    T = 2**bits
    g0 = T // B
    n1 = T % B                     # bins with g0+1 points
    n0 = B - n1
    c = np.array([g0, g0 + 1], float)
    m = np.array([n0, n1], float)
    w = c * B / T
    return pred_weights(w, m, N, B)


# ------------------------------------------------------------- logistic map
def pred_logistic(N, B, mu=LOG_MU, L=LOG_L):
    """Tail of mu states visited once, then a cycle of L states visited
    (N-mu)/L times.  Tail and cycle points fall in bins as independent
    Poisson samples (densities mu/B and L/B)."""
    if N <= mu + L:
        lam = N / B
        return 1.0, null_sd(B, N), np.exp(-lam)
    q, r = divmod(N - mu, L)
    f = r / L
    D = (mu + L * f * (1 - f) + (N - mu) ** 2 / L) / N
    pe = np.exp(-(mu + L) / B)
    # SD: Monte Carlo free estimate from the Poisson mixture
    lam = N / B
    k, p = coc_logistic(N, B, kmax=int(lam * 40 + 200), mu=mu, L=L)
    m2 = np.sum(p * (k - lam) ** 2)
    m4 = np.sum(p * (k - lam) ** 4)
    sd = np.sqrt(B * (m4 - m2**2)) / (lam * (B - 1))
    return D, sd, pe


def coc_logistic(N, B, kmax, mu=LOG_MU, L=LOG_L):
    q, r = divmod(N - mu, L)
    f = r / L
    k = np.arange(kmax + 1)
    out = np.zeros(kmax + 1)
    lt, lc = mu / B, L / B
    tmax = int(lt + 12 * np.sqrt(lt) + 15)
    cmax = int(lc + 12 * np.sqrt(lc) + 15)
    for c in range(cmax + 1):
        pc = stats.poisson.pmf(c, lc)
        if c == 0:
            base = np.zeros(kmax + 1); base[0] = 1.0
        else:
            base = stats.binom.pmf(k - q * c, c, f)
        for t in range(tmax + 1):
            pt = stats.poisson.pmf(t, lt)
            if pt * pc < 1e-18:
                continue
            sh = np.zeros(kmax + 1)
            if t <= kmax:
                sh[t:] = base[:kmax + 1 - t]
            out += pc * pt * sh
    return k, out
