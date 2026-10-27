"""Independent numpy LCG implementations (from the definitions in rnglib.f90)."""
import numpy as np


def lcg_stream(a, m, seed, n, c=0, block=1 << 20):
    """Return integer states X_1..X_n of X_{k+1} = (a X_k + c) mod m, X_0 = seed (multiplicative if c=0)."""
    out = np.empty(n, dtype=np.int64)
    assert c == 0
    A = np.empty(block, dtype=np.int64)
    p = 1
    for j in range(block):
        p = p * a % m
        A[j] = p
    x = seed % m
    i = 0
    while i < n:
        k = min(block, n - i)
        if m == 1 << 31:
            out[i:i + k] = (x * A[:k]) & (m - 1)
        else:
            # x*A < 2^62 for m < 2^31
            out[i:i + k] = (x * A[:k]) % m
        x = int(out[i + k - 1])
        i += k
    return out


def randu(seed, n):
    return lcg_stream(65539, 1 << 31, seed, n)


def randm(seed, n):
    return lcg_stream(6539, 1 << 31, seed, n)


def minstd(seed, n):
    return lcg_stream(16807, (1 << 31) - 1, seed, n)
