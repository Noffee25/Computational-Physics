"""Search for short integer relations h0 u_n + sum_j h_j u_{n+k_j} = integer
(i.e. h0 + sum h_j a^{k_j} = 0 mod m) at the lags that a 16x16 typewriter
sweep makes physically adjacent (1, 16, N-16, N-1, N for N = 256 or 257)."""
import itertools, math
def search(a, m, lags, H):
    A = [pow(a, k, m) for k in lags]
    best = None
    rng = range(-H, H + 1)
    for hs in itertools.product(rng, repeat=len(lags)):
        if all(h == 0 for h in hs):
            continue
        r = (-sum(h * Ak for h, Ak in zip(hs, A))) % m
        for h0 in (r, r - m):
            if abs(h0) <= H:
                nrm = math.sqrt(h0 * h0 + sum(h * h for h in hs))
                if best is None or nrm < best[0]:
                    best = (nrm, (h0,) + hs)
    return best
gens = {"RANDU": (65539, 2**31), "RANDM": (6539, 2**31), "MINSTD": (16807, 2**31 - 1)}
lagsets = {"N=256 sweep": [1, 16, 240, 255, 256], "N=257 sweep": [1, 16, 241, 256, 257]}
for gname, (a, m) in gens.items():
    for lname, ls in lagsets.items():
        out = []
        for r in (1, 2):
            for comb in itertools.combinations(ls, r):
                H = 40 if r == 1 else 25
                b = search(a, m, list(comb), H)
                if b:
                    out.append((round(b[0], 1), comb, b[1]))
        out.sort()
        print(f"{gname:7s} {lname}: shortest relations (norm, lags, coeffs h0,h..): {out[:3]}")
