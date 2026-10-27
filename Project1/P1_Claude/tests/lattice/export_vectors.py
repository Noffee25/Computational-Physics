"""Write vectors.txt for lat_verify.f90: one line per normal vector:
gen t mode G expected_h_num expected_h_den s_1..s_t   (padded to 8)"""
import sys, json
sys.path.insert(0, "/Users/noah/Downloads/ComputationalPhysicsPHYS25000/Project1/P1_Claude/tests/lattice")
from fractions import Fraction
from spectral import h_gap, sdotv, dot, canon
R = json.load(open('spectral_results.json'))
code = {'RANDU': 1, 'RANDM': 2, 'MINSTD': 3, 'MINSTD2': 4, 'LCG69069': 5}
lines = []
for r in R:
    g = code[r['gen']]; t = r['t']
    vecs = [('short', r['vecs'][0])]
    if r['gen'] in ('RANDU', 'RANDM'):
        if r['dmax_vec'] != r['vecs'][0]:
            vecs.append(('dmax', r['dmax_vec']))
        if canon(r['vec_m31']) != canon(r['vecs'][0]):
            vecs.append(('naive31', r['vec_m31']))
    for tag, s in vecs:
        if r['gen'] in ('RANDU', 'RANDM'):
            h, k = h_gap(s, r['a'], 31, 777)
            mode, G = (1, 0) if k >= 28 else (2, 2 ** k)
        else:
            h = Fraction(1); mode, G = 1, 0
        lines.append('%d %d %d %d %d %d %s   ! %s %s' % (g, t, mode, G, h.numerator, h.denominator,
                     ' '.join(str(c) for c in (s + [0] * (8 - t))), r['gen'], tag))
open('vectors.txt', 'w').write('\n'.join(lines) + '\n')
print('\n'.join(lines))
