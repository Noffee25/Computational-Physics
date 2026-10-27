import numpy as np
m = 2**31 - 1; a = 16807
def states(seed, n):
    out = np.empty(n, dtype=np.int64); x = seed
    # vectorised by jumping: simple loop in chunks
    for i in range(n):
        x = (a * x) % m; out[i] = x
    return out
def S_of(X, nvec, N):
    T = X[:3 * N].reshape(N, 3) / m
    ph = 2 * np.pi * (T @ np.array(nvec, float))
    return (np.cos(ph).sum() ** 2 + np.sin(ph).sum() ** 2) / N
X = states(777, 3 * 1000000)
for N in (100000, 300000, 1000000):
    print('seed 777 N=%d S(1,-10,2)=%.3f' % (N, S_of(X, (1, -10, 2), N)))
# second block of 1e5 particles from the same stream and other seeds
for blk in range(1, 6):
    T = X[3 * 100000 * blk: 3 * 100000 * (blk + 1)]
    print('block', blk, 'S=%.3f' % S_of(T, (1, -10, 2), 100000))
for seed in (1, 12345, 31415, 2718, 99991):
    Y = states(seed, 300000)
    print('seed', seed, 'S=%.3f' % S_of(Y, (1, -10, 2), 100000))
