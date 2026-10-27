# exact full-period lag-1 autocorrelation of output bit j for RANDU / RANDM
# bit j of x = k/2^31 is the top bit of k mod 2^m, m = 32 - j; the residue
# orbit (seed = 1 mod 8, a = 3 mod 8) is {r : r = 1 or 3 mod 8}
import numpy as np, json, sys
sys.path.insert(0, "/Users/noah/Downloads/ComputationalPhysicsPHYS25000/Project1/P1_Claude/tests/digits")
from bit_analysis import load, stats
OUT = "/Users/noah/Downloads/ComputationalPhysicsPHYS25000/Project1/P1_Claude/tests/digits/out/"
res = {}
for g, a in ((1, 65539), (2, 6539)):
    pred = {}
    for j in range(5, 31):
        m = 32 - j
        r = np.arange(2 ** m, dtype=np.int64)
        r = r[(r % 8 == 1) | (r % 8 == 3)] if m >= 3 else r[r % 2 == 1]
        b0 = (r >> (m - 1)) & 1
        b1 = (((a % 2 ** m) * r) % 2 ** m >> (m - 1)) & 1
        f = b0.mean()
        pred[j] = float(((b0 * b1).mean() - f * f) / (f * (1 - f))) if 0 < f < 1 else None
    d = load(g, 0); f, zf, const, cor, zc = stats(d)
    print("gen", g, "a =", a)
    for j in range(5, 31):
        print(f"  j={j:2d} pred {pred[j] if pred[j] is not None else float('nan'):+.6f}  meas {cor[j-1,0]:+.6f}  1/a={1/a:.2e}")
    res[g] = pred
json.dump(res, open(OUT + "bitcor_predict.json", "w"))
