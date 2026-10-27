"""Exact spectral test, t = 2..8, for the course LCGs and two reference multipliers.
Writes spectral_results.json and prints a table."""
import sys, json, time
sys.path.insert(0, "/Users/noah/Downloads/ComputationalPhysicsPHYS25000/Project1/P1_Claude/tests/lattice")
from math import sqrt, factorial
from fractions import Fraction
from spectral import (shortest, coset_spectral, merit, mu_knuth, lll, dual_basis,
                      enumerate_ball, dot, h_gap, canon, sdotv, in_dual)

def plane_count(s, residues, m):
    """Number of hyperplanes s.x in (r + mZ)/m (r in residues) meeting the open cube."""
    neg = sum(c for c in s if c < 0)
    pos = sum(c for c in s if c > 0)
    n = 0
    for r in residues:
        f = Fraction(r, m)
        # integers j with neg < j + f < pos
        lo = neg - f; hi = pos - f
        jlo = int(lo) + 1 if lo == int(lo) or lo > 0 else int(lo)
        # do it simply
        import math
        jlo = math.floor(lo) + 1
        jhi = math.ceil(hi) - 1
        n += max(0, jhi - jlo + 1)
    return n

def residues_coset(s, a, e, x0=1):
    m = 1 << e
    g = sdotv(s, a, m)
    vals = set()
    from spectral import v2
    k = v2(g, e)
    d = 1 << min(k + 3, e)
    base = [(x0 * g) % d, (3 * x0 * g) % d]
    if d > (1 << 12):  # only small residue sets are needed (k >= e-3)
        pass
    out = set()
    for b in base:
        r = b
        while r < m:
            out.add(r); r += d
            if len(out) > 64: break
    return sorted(out)

GENS = [
    # name, a, c, M (true modulus), kind
    ("RANDU", 65539, 0, 2 ** 31, "pow2_3mod8"),
    ("RANDM", 6539, 0, 2 ** 31, "pow2_3mod8"),
    ("MINSTD", 16807, 0, 2 ** 31 - 1, "prime"),
    ("MINSTD2", 48271, 0, 2 ** 31 - 1, "prime"),
    ("LCG69069", 69069, 1, 2 ** 32, "mixed"),
]
out = []
for name, a, c, M, kind in GENS:
    for t in range(2, 9):
        t0 = time.time()
        rec = dict(gen=name, a=a, c=c, M=M, t=t)
        if kind == "pow2_3mod8":
            e = M.bit_length() - 1
            r = coset_spectral(a, e, t)
            Meff = M // 2
            nu2 = r['nu2_m1']; vecs = r['vec_m1']
            rec.update(nu2_m31=r['nu2_m0'], vec_m31=r['vec_m0'][0], nu2_m29=r['nu2_m2'], vec_m29=r['vec_m2'][0],
                       nu2_m28=r['nu2_m3'],
                       dmax=r['dmax'], dmax_vec=r['dmax_vec'], dmax_h=str(r['dmax_h']), uniform_d=r['uniform_d'])
            # plane counts
            s = vecs[0]
            res = residues_coset(s, a, e)
            rec['residues_short'] = [x for x in res]
            rec['planes_short'] = plane_count(s, res, M)
            s2 = r['dmax_vec']
            res2 = residues_coset(s2, a, e)
            rec['residues_dmax'] = res2 if len(res2) <= 8 else len(res2)
            rec['planes_dmax'] = plane_count(s2, res2, M)
        else:
            Meff = M
            nu2, vecs = shortest(a, M, t)
            rec['dmax'] = 1 / sqrt(nu2); rec['dmax_vec'] = vecs[0]
            s = vecs[0]
            rec['planes_short'] = plane_count(s, [0], M)
        rec.update(Meff=Meff, nu2=nu2, nu=sqrt(nu2), vecs=vecs, d=1 / sqrt(nu2),
                   S=merit(nu2, t, Meff), mu=mu_knuth(nu2, t, Meff),
                   marsaglia=(factorial(t) * Meff) ** (1.0 / t), L1=sum(abs(x) for x in vecs[0]))
        rec['sec'] = time.time() - t0
        out.append(rec)
        extra = ''
        if kind == "pow2_3mod8":
            extra = ' | nu2(2^31)=%d nu2(2^29)=%d  dmax=%.3e vec=%s h=%s planes=%d' % (
                rec['nu2_m31'], rec['nu2_m29'], rec['dmax'], rec['dmax_vec'], rec['dmax_h'], rec['planes_dmax'])
        print('%-8s t=%d nu2=%-12d nu=%10.2f 1/nu=%.3e S=%.4f mu=%.3g planes=%d vec=%s%s (%.1fs)' % (
            name, t, nu2, sqrt(nu2), 1 / sqrt(nu2), rec['S'], rec['mu'], rec['planes_short'], vecs, extra, rec['sec']))
        sys.stdout.flush()
json.dump(out, open('/Users/noah/Downloads/ComputationalPhysicsPHYS25000/Project1/P1_Claude/tests/lattice/spectral_results.json', 'w'), indent=1, default=str)
