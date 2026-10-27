import numpy as np
M = 2 ** 31; a = 6539
def seq(seed, n):
    L = 1 << 12
    first = np.empty(L, dtype=np.uint64); x = seed
    for i in range(L):
        x = (a * x) % M; first[i] = x
    AL = np.uint64(pow(a, L, M))
    nb = (n + L - 1) // L
    out = np.empty(nb * L, dtype=np.uint64); out[:L] = first
    for b in range(1, nb):
        out[b * L:(b + 1) * L] = (out[(b - 1) * L:b * L] * AL) % np.uint64(M)
    return out[:n].astype(np.float64) / M
q = np.pi / 6; N = 10 ** 7
x = seq(777, 3 * N + 3)
v = 2 * x - 1
for off in (0, 1, 2):
    w = v[off:off + 3 * N].reshape(N, 3)
    p = ((w ** 2).sum(1) <= 1).mean()
    print('offset %d: z = %.2f' % (off, (p - q) / np.sqrt(q * (1 - q) / N)))
# overlapping triplets
r2 = v[:-2] ** 2 + v[1:-1] ** 2 + v[2:] ** 2
p = (r2[:3 * N] <= 1).mean()
print('overlapping 3e7 triplets: z(naive, iid) = %.2f' % ((p - q) / np.sqrt(q * (1 - q) / (3 * N))))
# sanity: rnglib sequence check
print('first values', x[:3])
