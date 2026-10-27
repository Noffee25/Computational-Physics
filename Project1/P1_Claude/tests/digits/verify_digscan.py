# Independent check of digscan exact-digit counts (MT19937 via numpy, RANDU, MINSTD)
import numpy as np
from fractions import Fraction
D = "/Users/noah/Downloads/ComputationalPhysicsPHYS25000/Project1/P1_Claude/tests/digits/"
def load(fn, kind):
    out = {}
    for line in open(D + fn):
        t = line.split()
        if t[0] == kind:
            out[int(t[1])] = np.array([int(v) for v in t[2:]])
    return out
N = 100000
P = [1, 5, 10, 15, 16, 20, 30, 45, 52, 53, 54]
def digits_dyadic(k, B, p):
    return (k * 10**p >> B) % 10
# MT: numpy RandomState(777).random_sample == rnglib mt_res53 with seed 777
x = np.random.RandomState(777).random_sample(N)
k = [int(v) for v in (x * 2.0**53).astype(np.int64)]
sto = load("test_7.txt", "STO")
for p in P:
    c = np.bincount([digits_dyadic(kk, 53, p) for kk in k], minlength=10)
    assert (c == sto[p]).all(), ("MT", p, c, sto[p])
print("MT19937 exact digits OK")
# RANDU
s = 777; ks = []
for i in range(N):
    s = (65539 * s) % 2**31; ks.append(s)
sto = load("test_1.txt", "STO")
for p in P:
    c = np.bincount([digits_dyadic(kk, 31, p) for kk in ks], minlength=10)
    assert (c == sto[p]).all(), ("RANDU", p)
print("RANDU exact digits OK")
# MINSTD intended k/(2^31-1) and stored double
M = 2**31 - 1; s = 777; ks = []
for i in range(N):
    s = (16807 * s) % M; ks.append(s)
it = load("test_4.txt", "INT"); sto = load("test_4.txt", "STO")
for p in P:
    c = np.bincount([(kk * 10**p // M) % 10 for kk in ks], minlength=10)
    assert (c == it[p]).all(), ("MINSTD int", p)
    xs = [Fraction(float(kk) / float(M)) for kk in ks]  # exact value of the double
    c2 = np.bincount([int(xx * 10**p) % 10 for xx in xs], minlength=10)
    assert (c2 == sto[p]).all(), ("MINSTD sto", p)
print("MINSTD intended and stored digits OK")
