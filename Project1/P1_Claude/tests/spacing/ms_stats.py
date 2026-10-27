"""Exact null distributions and goodness-of-fit tools for the spacing tests.

n points iid U(0,1) plus the endpoints 0 and 1 give n+1 spacings
(D_0..D_n) ~ uniform on the simplex (Dirichlet(1,...,1)).

Maximum spacing M (Whitworth 1897; Fisher 1929; Darling 1953):
    P(M <= x) = sum_{k=0}^{n+1} (-1)^k C(n+1,k) (1-kx)_+^n
    (n+1) M - ln(n+1) -> standard Gumbel.
Minimum spacing m:
    P(m > x) = (1-(n+1)x)_+^n,  (n+1) n m -> Exp(1).
"""
import numpy as np
from scipy import special, stats

LN10 = np.log(10.0)


def maxsp_cdf_sf(n, x, kmax=120):
    """Return (log F, log S) with F = P(M<=x), S = 1-F, for arrays n, x.

    Alternating inclusion-exclusion sum when it is numerically safe
    (t = (n+1)x - ln(n+1) >= -2.5, terms bounded by ~e^{12}); otherwise
    the Poisson (extreme-value) approximation F = exp(-lambda),
    lambda = (n+1)(1-x)^n, which is the accurate regime for the far lower
    tail. log values are returned so that astronomically small tail
    probabilities keep their magnitude.
    """
    n = np.asarray(n, dtype=float)
    x = np.asarray(x, dtype=float)
    n, x = np.broadcast_arrays(n, x)
    logF = np.full(n.shape, np.nan)
    logS = np.full(n.shape, np.nan)
    ok = (n >= 1) & (x > 0)
    # x >= 1: impossible except degenerate
    top = ok & (x >= 1.0)
    logF[top] = 0.0
    logS[top] = -np.inf
    ok = ok & (x < 1.0)
    t = (n + 1) * x - np.log(n + 1)
    pois = ok & (t < -2.5)
    lam = np.exp(np.log(n[pois] + 1) + n[pois] * np.log1p(-x[pois]))
    logF[pois] = -lam
    logS[pois] = np.log(-np.expm1(-lam))
    exact = ok & ~pois
    ne = n[exact]
    xe = x[exact]
    S = np.zeros(ne.shape)
    logC = np.zeros(ne.shape)          # log C(n+1, k), built up exactly term by term
    for k in range(1, kmax + 1):
        valid = (k <= ne + 1) & (k * xe < 1.0)
        if not valid.any():
            break
        logC = logC + np.log(np.maximum(ne + 2 - k, 1.0)) - np.log(k)
        lt = logC + ne * np.log1p(-np.minimum(k * xe, 1 - 1e-16))
        term = np.where(valid, np.exp(lt), 0.0)
        S += ((-1) ** (k + 1)) * term
        if np.all(term < 1e-30 * np.maximum(np.abs(S), 1e-300)):
            break
    S = np.clip(S, 0.0, 1.0)
    F = 1.0 - S
    # use the leading term for S when it is tiny (relative accuracy)
    lead = np.exp(np.log(ne + 1) + ne * np.log1p(-xe))
    S_small = lead < 1e-8
    with np.errstate(divide="ignore"):
        lS = np.where(S_small, np.log(np.maximum(lead, 1e-320)) , np.log(np.maximum(S, 1e-320)))
        lF = np.log(np.maximum(F, 1e-320))
    # for very small S use log of leading term exactly (no underflow)
    lS = np.where(S_small, np.log(ne + 1) + ne * np.log1p(-xe), lS)
    logF[exact] = lF
    logS[exact] = lS
    return logF, logS


def minsp_cdf_sf(n, x):
    """(log F, log S) for the minimum spacing m: S = P(m > x) = (1-(n+1)x)^n."""
    n = np.asarray(n, dtype=float)
    x = np.asarray(x, dtype=float)
    z = (n + 1) * x
    with np.errstate(divide="ignore", invalid="ignore"):
        logS = np.where(z < 1, n * np.log1p(-np.minimum(z, 1 - 1e-16)), -np.inf)
        logF = np.log(-np.expm1(logS))
    logF = np.where(x <= 0, -np.inf, logF)
    logS = np.where(x <= 0, 0.0, logS)
    return logF, logS


def pit_from_logs(logF, logS):
    """PIT value u = F, computed accurately in both tails."""
    return np.where(logF < np.log(0.5), np.exp(logF), -np.expm1(logS))


def ks_test(u):
    u = np.sort(np.asarray(u))
    r = stats.kstest(u, "uniform")
    return r.statistic, r.pvalue


def ad_inf_cdf(z):
    """Marsaglia & Marsaglia (2004) asymptotic Anderson-Darling CDF."""
    z = float(z)
    if z <= 0:
        return 0.0
    if z < 2:
        return (np.exp(-1.2337141 / z) / np.sqrt(z) *
                (2.00012 + (.247105 - (.0649821 - (.0347962 - (.0116720 - .00168691 * z) * z) * z) * z) * z))
    return np.exp(-np.exp(1.0776 - (2.30695 - (.43424 - (.082433 - (.008056 - .0003146 * z) * z) * z) * z) * z))


def ad_test(u):
    """Anderson-Darling A^2 for U(0,1) and its asymptotic p-value."""
    u = np.sort(np.clip(np.asarray(u, dtype=float), 1e-300, 1 - 1e-16))
    R = len(u)
    i = np.arange(1, R + 1)
    A2 = -R - np.mean((2 * i - 1) * (np.log(u) + np.log1p(-u[::-1])))
    z = A2
    if z < 2:
        p = 1 - ad_inf_cdf(z)
    else:
        # upper tail directly: 1 - exp(-exp(g)) = -expm1(-exp(g))
        g = 1.0776 - (2.30695 - (.43424 - (.082433 - (.008056 - .0003146 * z) * z) * z) * z) * z
        p = -np.expm1(-np.exp(g))
    return A2, max(p, 0.0)


def norm_score_from_logs(logF, logS):
    """Phi^{-1}(F) computed from log F or log S without underflow."""
    out = np.empty(np.shape(logF))
    lo = logF < np.log(0.5)
    out[lo] = special.ndtri_exp(logF[lo])
    out[~lo] = -special.ndtri_exp(logS[~lo])
    return out


def stouffer(logF, logS):
    """Stouffer Z = sum Phi^{-1}(F_i)/sqrt(R); returns Z and log10 two-sided p."""
    z = norm_score_from_logs(np.asarray(logF), np.asarray(logS))
    Z = z.sum() / np.sqrt(len(z))
    lp = (np.log(2.0) + special.log_ndtr(-abs(Z))) / LN10
    return Z, lp


def gumbel_T(n, M):
    """Standardized statistic T = (n+1) M - ln(n+1) (-> Gumbel)."""
    n = np.asarray(n, dtype=float)
    return (n + 1) * np.asarray(M) - np.log(n + 1)


def expected_dups(N, G):
    """E[N - #distinct] for N draws uniform on G values."""
    N = np.asarray(N, dtype=float)
    return N - G * (-np.expm1(N * np.log1p(-1.0 / G)))


def invert_geff(N, D):
    """Solve expected_dups(N, G) = D for G (D>0)."""
    from scipy.optimize import brentq
    if D <= 0:
        return np.inf
    f = lambda lg: expected_dups(N, np.exp(lg)) - D
    lo, hi = np.log(1.0 + 1e-9), np.log(1e30)
    lo = np.log(max(N / 50.0, 1.0001))
    if f(lo) < 0:
        lo = np.log(1.0001)
    return float(np.exp(brentq(f, lo, hi)))
