import json, re
from math import sqrt
D = "/Users/noah/Downloads/ComputationalPhysicsPHYS25000/Project1/P1_Claude/tests/lattice/"
R = json.load(open(D + 'spectral_results.json'))
ver = {}
for g in ('RANDU', 'RANDM', 'MINSTD', 'MINSTD2', 'LCG69069'):
    for line in open(D + 'verify_%s.txt' % g):
        if line.startswith('#') or 'populations' in line:
            continue
        parts = line.split('|')
        head = parts[0].split()
        t = int(head[0]); h = head[2]; s = [int(x) for x in head[3:3 + t]]
        res = [int(x) for x in parts[1].split()[1:]]
        pl = parts[2].split(); planes = int(pl[0])
        gp = parts[3].split(); emp, pred, ratio = float(gp[0]), float(gp[1]), float(gp[2])
        ver[(g, tuple(s))] = dict(t=t, h=h, residues=res, planes=planes, emp=emp, pred=pred, ratio=ratio)
rows = []
for r in R:
    g = r['gen']; t = r['t']; s = r['vecs'][0]
    v = ver[(g, tuple(s))]
    row = dict(gen=g, t=t, nu2=r['nu2'], nu=r['nu'], s=s, d=r['d'], S=r['S'], mu=r['mu'], planes=v['planes'],
               residues=v['residues'], emp=v['emp'], ratio=v['ratio'], Meff=r['Meff'], marsaglia=r['marsaglia'])
    if g in ('RANDU', 'RANDM'):
        vd = ver[(g, tuple(r['dmax_vec']))]
        row.update(dmax=r['dmax'], dmax_vec=r['dmax_vec'], dmax_h=r['dmax_h'], dmax_planes=vd['planes'],
                   dmax_emp=vd['emp'], dmax_ratio=vd['ratio'], nu2_m31=r['nu2_m31'], nu2_m29=r['nu2_m29'])
        row['S_eff'] = r['S'] * (r['d'] / r['dmax'])
    rows.append(row)
json.dump(rows, open(D + 'summary.json', 'w'), indent=1)
# markdown table
print('| generator | t | $\\nu_t^2$ | shortest dual vector $s$ | $1/\\nu_t$ | $S_t$ | planes (occupied) | true max gap $d_t$ |')
print('|---|---|---|---|---|---|---|---|')
for row in rows:
    dm = row.get('dmax', row['d'])
    extra = '' if abs(dm - row['d']) < 1e-15 else ' via %s (h=%s)' % (tuple(row['dmax_vec']), row['dmax_h'])
    print('| %s | %d | %d | %s | %.3e | %.3f | %d | %.3e%s |' % (row['gen'], row['t'], row['nu2'], tuple(row['s']), row['d'], row['S'],
          row['planes'], dm, extra))
print()
for row in rows:
    if row['gen'] in ('RANDU', 'RANDM'):
        print(row['gen'], row['t'], 'nu2(2^31)=%d nu2(2^30)=%d nu2(2^29)=%d' % (row['nu2_m31'], row['nu2'], row['nu2_m29']),
              'd=1/nu30=%.4e dmax=%.4e ratio %.3f Seff=%.3f residues=%s' % (row['d'], row['dmax'], row['dmax'] / row['d'], row['S_eff'], row['residues']))
