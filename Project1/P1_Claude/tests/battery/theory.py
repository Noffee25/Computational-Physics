"""Closed-form predictions used by the battery family."""
import numpy as np
from fractions import Fraction as F
from math import gcd

def dedekind(h, k):
    """Dedekind sum s(h,k) by the reciprocity law (exact rational)."""
    s = F(0); sign = 1; h %= k
    while h != 0:
        s += sign * ((F(h, k) + F(k, h) + F(1, h * k)) / 12 - F(1, 4))
        h, k = k % h, h
        sign = -sign
    return s

def serial_dedekind(b, m):
    """Exact lag correlation of y = b x mod m over the full residue set 1..m-1."""
    return float(12 * m * dedekind(b, m) / ((m - 1) * (m - 2)))

def P1hat(n): return 1j * np.sqrt(3) / (np.pi * n)
def P2hat(n): return 3 * np.sqrt(5) / (np.pi ** 2 * n ** 2)
PH = {1: P1hat, 2: P2hat}

def triple_fourier(c_old, c_mid, c_new, pqr, hmax=20000):
    """<P_p(x_old) P_q(x_mid) P_r(x_new)> for points on the planes
    c_old x_old + c_mid x_mid + c_new x_new = integer (continuum limit; only the
    dual vectors h*(c_old,c_mid,c_new) are kept)."""
    p, q, r = pqr
    h = np.concatenate([np.arange(-hmax, 0), np.arange(1, hmax + 1)]).astype(float)
    # E[f(x) e^{2 pi i n x}] = conj(fhat(n)) = fhat(-n);  sum over h of prod fhat(-h c)
    val = PH[p](-h * c_old) * PH[q](-h * c_mid) * PH[r](-h * c_new)
    return float(np.real(val.sum()))

def randu_plane_cells(k):
    """Fraction of the k^3 cells that RANDU can ever occupy.  For a cell column
    (x,y) in [i/k,(i+1)/k) x [j/k,(j+1)/k), z = frac(6y - 9x) sweeps an interval
    of length (9+6)/k whose end points lie on the grid, i.e. exactly 15 of the
    k z-cells.  Hence the occupied fraction is min(1, 15/k)."""
    return min(1.0, 15.0 / k)

if __name__ == '__main__':
    print('RANDU 112', triple_fourier(9, -6, 1, (1, 1, 2)), np.sqrt(5) / 270)
    print('RANDU 121', triple_fourier(9, -6, 1, (1, 2, 1)), -np.sqrt(5) / 1620)
    print('RANDU 211', triple_fourier(9, -6, 1, (2, 1, 1)), np.sqrt(5) / 2430)
    print('RANDU 222', triple_fourier(9, -6, 1, (2, 2, 2)))
    print('RANDU 111', triple_fourier(9, -6, 1, (1, 1, 1)))
    print('RANDU (2,4) 112', triple_fourier(81, -18, 1, (1, 1, 2)), np.sqrt(5) / 7290)
    print('RANDU (1,4) 112', triple_fourier(27, -4, 1, (1, 1, 2)))
    for k in (8, 16, 32, 64, 100):
        print('k', k, 'crossed fraction', randu_plane_cells(k))
