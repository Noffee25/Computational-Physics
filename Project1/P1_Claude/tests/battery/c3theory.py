"""Predicted three-point Legendre correlations of RANDU / RANDM at every lag pair.

The orbit of an odd seed under a = 3 mod 8 is the coset C = {s : s = 1,3 mod 8}
of Z/2^31.  For x = s/2^31 on the orbit, x_mid = b1 x_old, x_new = b2 x_old
(b1 = a^(k-j), b2 = a^k mod 2^31).  The average of exp(2 pi i h.x) over the orbit
is exactly zero unless g = h_old + h_mid b1 + h_new b2 = 0 mod 2^28, in which case
it equals phi(u) = (w^u + w^(3u))/2, w = exp(2 pi i/8), u = (g mod 2^31)/2^28.
Hence  T_pqr = sum_h  P_p^(h_old) P_q^(h_mid) P_r^(h_new) phi(u(h)),
a sum over the 'coset dual lattice'.  We keep |h_mid|,|h_new| <= H and the
smallest representative of h_old, which captures all non-negligible terms."""
import numpy as np, sys
H = 48
MMAX = 400
M28 = 2**28; M31 = 2**31
def P1h(n): return 1j*np.sqrt(3)/(np.pi*n)
def P2h(n): return 3*np.sqrt(5)/(np.pi**2*n**2)
hm, hn = np.meshgrid(np.arange(-H, H+1), np.arange(-H, H+1), indexing='ij')
hm = hm.ravel().astype(np.int64); hn = hn.ravel().astype(np.int64)
keep = (hm != 0) & (hn != 0)
hm = hm[keep]; hn = hn[keep]
w = np.exp(2j*np.pi/8)
def predict(a, K):
    T1 = np.zeros((K+1, K+1)); T2 = np.zeros((K+1, K+1))
    apow = [pow(a, e, M31) for e in range(K+1)]
    for k in range(2, K+1):
        b2 = apow[k]
        for j in range(1, k):
            b1 = apow[k-j]
            r = (hm*b1 + hn*b2) % M28          # need h_old = -r mod 2^28
            ho = -r
            ho = np.where(ho < -M28//2, ho + M28, ho)
            sel = (ho != 0) & (np.abs(ho) <= 4000)
            if not sel.any(): continue
            h0 = ho[sel]; h1 = hm[sel]; h2 = hn[sel]
            # keep primitive vectors, one per +-pair, then sum all multiples m h
            g0 = np.gcd(np.gcd(np.abs(h0), np.abs(h1)), np.abs(h2))
            prim = (g0 == 1) & ((h2 > 0) | ((h2 == 0) & (h1 > 0)))
            h0 = h0[prim]; h1 = h1[prim]; h2 = h2[prim]
            g = (h0 + h1*b1 + h2*b2) % M31
            u = (g // M28) % 8
            mm = np.concatenate([np.arange(-MMAX, 0), np.arange(1, MMAX+1)])[:, None]
            um = (mm*u[None, :]) % 8
            phi = (w**um + w**(3*um))/2
            A0 = P1h(mm*h0[None, :]); A1 = P1h(mm*h1[None, :])
            t1 = np.real((A0*A1*P1h(mm*h2[None, :])*phi).sum())
            t2 = np.real((A0*A1*P2h(mm*h2[None, :])*phi).sum())
            T1[j, k] = T1[k, j] = t1; T2[j, k] = T2[k, j] = t2
    return T1[1:, 1:], T2[1:, 1:]
if __name__ == '__main__':
    K = int(sys.argv[1]) if len(sys.argv) > 1 else 256
    for name, a in (('RANDU', 65539), ('RANDM', 6539)):
        T1, T2 = predict(a, K)
        np.savez(f'out/c3theory_{name}.npz', T1=T1, T2=T2)
        for lab, T in (('T1', T1), ('T2', T2)):
            iu = np.triu_indices(K, 1); v = T[iu]; o = np.argsort(-np.abs(v))[:6]
            print(name, lab, [(iu[0][t]+1, iu[1][t]+1, '%.3e' % v[t]) for t in o])
