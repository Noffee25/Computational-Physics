"""Pair correlation of Metropolis acceptance indicators I = [u < w] at lag k
for the power-of-two LCGs, w = exp(-4K_c) and exp(-8K_c).
rho_k(w) = P(I_{n+k}=1 | I_n=1)/w - 1 ; null sd ~ sqrt((1-w)/(n w))."""
import numpy as np
M = 2 ** 31
w4 = np.exp(-4 * np.log(1 + np.sqrt(2)) / 2); w8 = w4 ** 2
def stream(a, seed, nblk, B=256):
    x = np.empty(B, dtype=np.uint64); s = seed
    for i in range(B):
        s = (a * s) % M; x[i] = s
    A = np.uint64(pow(a, B, M)); out = np.empty((nblk, B), dtype=np.uint64)
    out[0] = x
    for j in range(1, nblk):
        out[j] = (A * out[j - 1]) % np.uint64(M)
    return out.reshape(-1).astype(np.float64) / M
res = {}
for name, a in [("RANDM", 6539), ("RANDU", 65539)]:
    u = stream(a, 777, 200000)
    n = len(u)
    print(name, "n =", n)
    for w, wl in [(w4, "w4"), (w8, "w8")]:
        I = u < w
        line = []
        for k in [1, 2, 16, 255, 256, 257, 512]:
            a1 = I[:-k]; b1 = I[k:]
            p = (a1 & b1).sum() / a1.sum()
            rho = p / w - 1; sd = np.sqrt((1 - w) / (a1.sum() * w))
            line.append(f"k={k}: {rho:+.4f} ({rho/sd:+.1f}sd)")
        print(f"  {wl}={w:.4f}: " + "  ".join(line))
