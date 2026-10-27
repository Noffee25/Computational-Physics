# rho structure (tail, cycle) of the 8-digit middle-square map used by rnglib
def f(s): return (s * s // 10000) % 100000000
for seed in [777, 12345678, 31415926, 77777777, 98765432]:
    seen = {}; s = seed % 100000000; n = 0
    while s not in seen:
        seen[s] = n; s = f(s); n += 1
    print(seed, "tail", seen[s], "cycle", n - seen[s], "cycle entry", s)
