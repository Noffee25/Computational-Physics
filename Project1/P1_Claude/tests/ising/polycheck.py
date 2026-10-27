"""GF(2) polynomial checks for the shift-register generators used.
(1) Irreducibility of x^q + sum x^{taps} + 1 for prime q (Rabin test):
    f irreducible  <=>  x^(2^q) = x mod f  and  gcd(x^2 - x, f) = 1.
    For Mersenne-prime q (2^q - 1 prime) irreducible implies primitive.
(2) Berlekamp-Massey on a single bit column of R250 decimated by d:
    the minimal recurrence of the decimated sequence (degree and number of taps).
"""
import numpy as np

def sq(a):
    return int("0".join(bin(a)[2:]), 2)

def red(a, q, taps):
    mask = (1 << q) - 1
    while a >> q:
        H = a >> q
        a &= mask
        a ^= H
        for t in taps:
            a ^= H << t
    return a

def gcd2(a, b):
    while b:
        # a mod b in GF(2)
        db = b.bit_length()
        while a.bit_length() >= db:
            a ^= b << (a.bit_length() - db)
        a, b = b, a
    return a

def irreducible_prime_degree(q, taps):
    f = (1 << q) | 1
    for t in taps:
        f |= 1 << t
    x = 2
    a = x
    for _ in range(q):
        a = red(sq(a), q, taps)
    if a != x:
        return False
    return gcd2(f, 0b110) == 1  # x^2 + x

def berlekamp_massey(s):
    s = list(s)
    n = len(s)
    C = [1] + [0] * n; B = [1] + [0] * n
    Lc, m, b = 0, 1, 1
    for i in range(n):
        d = s[i]
        for j in range(1, Lc + 1):
            d ^= C[j] & s[i - j]
        if d == 0:
            m += 1
        elif 2 * Lc <= i:
            T = C[:]
            for j in range(n - m + 1):
                C[j + m] ^= B[j]
            Lc = i + 1 - Lc; B = T; m = 1
        else:
            for j in range(n - m + 1):
                C[j + m] ^= B[j]
            m += 1
    return Lc, C[:Lc + 1]

if __name__ == "__main__":
    MERSENNE = {17, 31, 89, 127, 521, 607, 1279, 2203, 2281, 3217, 4253, 4423, 9689, 9941}
    cases = [(89, [38]), (127, [63]), (127, [30]), (127, [1]), (521, [168]), (521, [32]), (607, [273]),
             (1279, [418]), (2281, [1029]), (9689, [4187]), (9689, [471, 1586, 6988])]
    for q, taps in cases:
        ok = irreducible_prime_degree(q, taps)
        print(f"x^{q} + " + " + ".join(f"x^{t}" for t in sorted(taps, reverse=True)) + f" + 1 : irreducible={ok}, q Mersenne exponent={q in MERSENNE} -> primitive={ok and q in MERSENNE}")
    # Berlekamp-Massey on bit column of the R250 sequence x_n = x_{n-103} ^ x_{n-250}
    rng = np.random.default_rng(1)
    nbits = 3000 * 5 + 1000
    b = list(rng.integers(0, 2, 250))
    while len(b) < nbits * 1:
        b.append(b[-103] ^ b[-250])
    for d in [1, 2, 3, 4, 5, 7]:
        sub = b[300::d][:1200]
        Lc, C = berlekamp_massey(sub)
        taps = [j for j in range(1, Lc + 1) if C[j]]
        print(f"R250 bit column decimated by {d}: minimal recurrence degree {Lc}, number of feedback taps {len(taps)}, smallest taps {taps[:6]}")
