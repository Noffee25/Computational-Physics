import sys
sys.path.insert(0, "/Users/noah/Downloads/ComputationalPhysicsPHYS25000/Project1/P1_Claude/tests/ising")
import mpmath as mp
from exact import e_C_kaufman
mp.mp.dps = 40
with open("exact_K.txt", "w") as fh:
    fh.write("# L K e C\n")
    for K in ["0.36", "0.40", "0.42", "0.44068679350977151", "0.46", "0.50", "0.55"]:
        e, C = e_C_kaufman(mp.mpf(K), 16, 16)
        fh.write(f"16 {K} {mp.nstr(e, 18)} {mp.nstr(C, 18)}\n")
        print(16, K, mp.nstr(e, 14), mp.nstr(C, 14))
