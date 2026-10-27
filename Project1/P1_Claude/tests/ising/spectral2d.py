"""2D spectral test at lag k for multiplicative LCGs: shortest nonzero (h1,h2)
with h1 + a^k h2 = 0 mod m (Gauss-Lagrange reduction). Pairs (u_n, u_{n+k})
lie on parallel lines spaced 1/nu_2(k) apart."""
import math
def nu2(a, m):
    # basis of dual lattice: (m, 0) and (-a, 1)
    u = [m, 0]; v = [-(a % m), 1]
    def n2(x): return x[0] * x[0] + x[1] * x[1]
    if n2(u) < n2(v): u, v = v, u
    while True:
        # reduce u by v
        mu = round((u[0] * v[0] + u[1] * v[1]) / n2(v))
        u = [u[0] - mu * v[0], u[1] - mu * v[1]]
        if n2(u) >= n2(v):
            return math.sqrt(n2(v)), v
        u, v = v, u
gens = {"RANDU": (65539, 2**31), "RANDM": (6539, 2**31), "MINSTD": (16807, 2**31 - 1)}
for name, (a, m) in gens.items():
    res = []
    for k in range(1, 1025):
        nu, vec = nu2(pow(a, k, m), m)
        res.append((nu, k, vec))
    res.sort()
    print(name, "lags relevant to the 16x16 sweep:", {k: round(nu2(pow(a, k, m), m)[0], 1) for k in [1, 2, 15, 16, 17, 240, 255, 256, 257, 272]})
    print("   worst lags k<=1024:", [(k, round(nu, 1), vec) for nu, k, vec in res[:6]])
