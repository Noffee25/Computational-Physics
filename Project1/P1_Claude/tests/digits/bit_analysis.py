# per-bit statistics from bitscan output
import numpy as np
D = "/Users/noah/Downloads/ComputationalPhysicsPHYS25000/Project1/P1_Claude/tests/digits/out/"
NAMES = ['RANDU', 'RANDM', 'RAND', 'MINSTD', 'R250', 'RANNUM', 'MT19937', 'MIDSQ', 'LOGISTIC']

def load(g, prec):
    fn = D + f"bit_g{g}_p{prec}.txt"
    hdr = open(fn).readline().split()
    N = int(hdr[hdr.index('N') + 1]); nlag = [int(v) for v in hdr[hdr.index('nlag') + 1: hdr.index('nlag') + 4]]
    n250 = int(hdr[hdr.index('n250') + 1])
    a = np.loadtxt(fn, comments='#')
    j = a[:, 0].astype(int)
    order = np.argsort(j); a = a[order]
    return dict(N=N, nlag=nlag, n250=n250, j=a[:, 0].astype(int), ones=a[:, 1], a11=a[:, 2:5],
                v250=a[:, 5], vtent=a[:, 6], vcopy=a[:, 7])

def stats(d):
    N = d['N']; f = d['ones'] / N
    zf = (d['ones'] - N / 2) / np.sqrt(N / 4)
    const = (d['ones'] == 0) | (d['ones'] == N)
    cor = np.full((53, 3), np.nan)
    for il in range(3):
        nl = d['nlag'][il]
        with np.errstate(invalid='ignore', divide='ignore'):
            cor[:, il] = (d['a11'][:, il] / nl - f * f) / (f * (1 - f))
    cor[const, :] = np.nan
    zc = cor * np.sqrt(N)
    return f, zf, const, cor, zc

if __name__ == "__main__":
    for prec in (0, 1):
        print("\n==== prec", prec)
        for g in range(1, 10):
            d = load(g, prec); f, zf, const, cor, zc = stats(d)
            print(f"{NAMES[g-1]:9s} const bits:", list(d['j'][const]))
            print("   freq z  :", " ".join(f"{v:.0f}" if abs(v) < 1e9 else "x" for v in zf))
            print("   lag1 cor:", " ".join("c" if np.isnan(v) else f"{v:.3f}" for v in cor[:, 0]))
            print("   v250 frac:", " ".join(f"{v / d['n250']:.2f}" for v in d['v250']))
            print("   vtent frac:", " ".join(f"{v / (d['N'] - 1):.3f}" for v in d['vtent']))
            print("   vcopy frac (j=32..53):", " ".join(f"{v / d['N']:.3f}" for v in d['vcopy'][31:]))
