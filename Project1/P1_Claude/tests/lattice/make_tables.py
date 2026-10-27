import json, numpy as np
D = "/Users/noah/Downloads/ComputationalPhysicsPHYS25000/Project1/P1_Claude/tests/lattice/"
rows = json.load(open(D + 'summary.json'))
R = {(r['gen'], r['t']): r for r in rows}
def vec(s): return '(' + ','.join(str(c) for c in s) + ')'
def off(r):
    res = r['residues']
    if not res: return '-'
    M = 2 ** 31 if r['gen'] in ('RANDU', 'RANDM') else (2 ** 31 - 1 if 'MINSTD' in r['gen'] else 2 ** 32)
    v = res[0] / M
    return '0' if v == 0 else ('1/2' if abs(v - 0.5) < 1e-12 else '%.4f' % v)
out = []
out.append('### Table A. RANDU versus RANDM: shortest dual vectors (lattice coset, modulus 2^30)\n')
out.append('| $t$ | RANDU $\\mathbf{s}$ | planes | $1/\\nu_t$ | $S_t$ | RANDM $\\mathbf{s}$ | offset $c$ | planes | $1/\\nu_t$ | true $d_t$ | $S_t$ |')
out.append('|---|---|---|---|---|---|---|---|---|---|---|')
for t in range(2, 9):
    u = R[('RANDU', t)]; m = R[('RANDM', t)]
    out.append('| %d | %s | %d | %.3e | %.4f | %s | %s | %d | %.3e | %.3e | %.3f |' % (
        t, vec(u['s']), u['planes'], u['d'], u['S'], vec(m['s']), off(m), m['planes'], m['d'], m['dmax'], m['S']))
out.append('\nPlanes: number of distinct hyperplanes $\\mathbf{s}\\cdot\\mathbf{x} \\in c + \\mathbb{Z}$ actually occupied over the full period ($2^{29}$ tuples). Offset $c$ for RANDU is 0 in every $t$ except $t=2$ where it is 1/2.\n')
out.append('### Table B. Reference LCGs (exact spectral test, full-period plane counts)\n')
out.append('| $t$ | MINSTD $\\nu_t^2$ | planes | $S_t$ | $\\mu_t$ | MINSTD2 $\\nu_t^2$ | planes | $S_t$ | $\\mu_t$ | LCG69069 $\\nu_t^2$ | planes | $S_t$ | $\\mu_t$ |')
out.append('|---|---|---|---|---|---|---|---|---|---|---|---|---|')
for t in range(2, 9):
    cells = []
    for g in ('MINSTD', 'MINSTD2', 'LCG69069'):
        r = R[(g, t)]
        cells.append('%d | %d | %.3f | %.3g' % (r['nu2'], r['planes'], r['S'], r['mu']))
    out.append('| %d | %s |' % (t, ' | '.join(cells)))
out.append('\n### Table C. Coset correction for $m=2^{31}$, $a \\equiv 3 \\pmod 8$ (rows where conventions disagree)\n')
out.append('| gen | $t$ | $\\nu_t^2(2^{31})$ | $\\nu_t^2(2^{30})$ | $\\nu_t^2(2^{29})$ | $1/\\nu_t(2^{30})$ | $1/\\nu_t(2^{29})$ | true max gap $d_t$ | from $\\mathbf{s}$ | $h(\\mathbf{s})$ |')
out.append('|---|---|---|---|---|---|---|---|---|---|')
for g in ('RANDU', 'RANDM'):
    for t in range(2, 9):
        r = R[(g, t)]
        if r['nu2_m31'] != r['nu2'] or r['nu2_m29'] != r['nu2'] or abs(r['dmax'] - r['d']) > 1e-15:
            out.append('| %s | %d | %d | %d | %d | %.4e | %.4e | %.4e | %s | %s |' % (g, t, r['nu2_m31'], r['nu2'], r['nu2_m29'],
                       r['d'], 1 / np.sqrt(r['nu2_m29']), r['dmax'], vec(r['dmax_vec']), r['dmax_h']))
print('\n'.join(out))
open(D + 'tables.md', 'w').write('\n'.join(out) + '\n')
