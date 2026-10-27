# Exact predictions of the digit-frequency chi-square from the output resolution
#  (a) uniform B-bit dyadic source k/2^B, k uniform on [0, 2^B)
#  (b) RANDU / RANDM orbit: k = 1 or 3 mod 8 (seed 777 = 1 mod 8), full-period average
#  (c) continuous uniform x rounded to nearest REAL*4 (24-bit significand)
import numpy as np, json
N = 1e8
OUT = "/Users/noah/Downloads/ComputationalPhysicsPHYS25000/Project1/P1_Claude/tests/digits/out/"

def pi_uniform(B, p):
    if p > B:
        c = np.zeros(10); c[0] = 1.0; return c
    q = 2 ** (B - p)
    # c_d = #{m in [0, 2^(B+1-p)) : floor(5 m / q) = d} = ceil((d+1) q/5) - ceil(d q/5)
    c = np.array([-(-(d + 1) * q // 5) + (d * q // 5 if (d * q) % 5 == 0 else -(-(d * q) // 5)) * -1 for d in range(10)], dtype=object)
    c = np.array([(-(-(d + 1) * q // 5)) - (-(-(d * q) // 5)) for d in range(10)], dtype=object)
    tot = sum(c)
    assert tot == 2 ** (B + 1 - p)
    return np.array([float(v) / float(tot) for v in c])

def pi_randu(p):
    # residues r mod 2^(32-p) with r = 1, 3 mod 8; digit = floor(r 5^p / 2^(31-p)) mod 10
    nb = 32 - p
    if nb < 3:
        return None
    r = np.arange(2 ** nb, dtype=np.int64)
    r = r[(r % 8 == 1) | (r % 8 == 3)]
    mod = 10 * 2 ** (31 - p)
    c5 = pow(5, p, mod)
    d = ((r * c5) % mod) >> (31 - p)
    cnt = np.bincount(d, minlength=10).astype(float)
    return cnt / cnt.sum()

def pi_single(p, kmax=34):
    # x uniform on [0,1) rounded to nearest float: binade k holds m/2^(24+k), m in [2^23, 2^24), weight 2^-(k+1)
    m0 = np.arange(2 ** 23, 2 ** 24, dtype=np.int64)
    tot = np.zeros(10)
    for k in range(kmax):
        B = 24 + k
        r = m0.copy(); mask = (1 << B) - 1
        for _ in range(p):
            r = r * 10
            d = r >> B
            r &= mask
        tot += 2.0 ** -(k + 1) * np.bincount(d, minlength=10) / 2 ** 23
    return tot / tot.sum()

def lam(pi, N=N):
    return 10 * N * ((pi - 0.1) ** 2).sum()

if __name__ == "__main__":
    res = {"N": N, "uniform": {}, "randu": {}, "single": {}}
    for B in (22, 24, 31, 32, 53):
        res["uniform"][B] = [lam(pi_uniform(B, p)) for p in range(1, 61)]
    res["randu"] = [lam(pi_randu(p)) if p >= 4 and p <= 29 else None for p in range(1, 61)]
    print("uniform B=22", " ".join(f"{v:.3g}" for v in res["uniform"][22][:16]))
    print("uniform B=32", " ".join(f"{v:.3g}" for v in res["uniform"][32][14:26]))
    print("uniform B=53", " ".join(f"{v:.3g}" for v in res["uniform"][53][34:50]))
    print("randu", " ".join(f"{v:.3g}" if v is not None else "-" for v in res["randu"][:30]))
    sing = []
    for p in range(1, 21):
        sing.append(lam(pi_single(p)))
    res["single"] = sing
    print("single", " ".join(f"{v:.3g}" for v in sing))
    json.dump(res, open(OUT + "predict.json", "w"))
