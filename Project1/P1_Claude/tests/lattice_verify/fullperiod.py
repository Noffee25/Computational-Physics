"""Exact full-period counts for RANDU and RANDM (seed class x0 = 777 = 1 mod 8):
all 2^29 overlapping triplets are (x, ax, a^2 x) mod 2^31 for x = 1 or 3 mod 8."""
import numpy as np, math, sys
M = 1 << 31
R62 = np.uint64(1 << 62)
for name, a in [('RANDU', 65539), ('RANDM', 6539)]:
    planes = np.zeros(15, dtype=np.int64)
    nin = nin_f = ncen = 0
    ntot = 0
    # residue classes for the RANDM two-offset normal (58,-375,-69)
    res = set()
    CH = 1 << 23
    for cls in (1, 3):
        for i0 in range(0, 1 << 28, CH):
            x = (np.arange(i0, i0 + CH, dtype=np.int64) << 3) + cls
            y = (a * x) & (M - 1)
            z = (a * y) & (M - 1)
            ntot += len(x)
            if name == 'RANDU':
                S = 9 * x - 6 * y + z
                assert np.all(S % M == 0)
                planes += np.bincount(S // M + 5, minlength=15)
            else:
                g = (58 * x - 375 * y - 69 * z) % M
                res |= set(np.unique(g).tolist())
            v = [np.abs(2 * w - M).astype(np.uint64) for w in (x, y, z)]
            r2 = v[0] * v[0] + v[1] * v[1] + v[2] * v[2]
            inb = r2 <= R62
            nin += int(inb.sum())
            vf = [2.0 * w / M - 1.0 for w in (x, y, z)]
            r2f = vf[0] ** 2 + vf[1] ** 2 + vf[2] ** 2
            nin_f += int((r2f <= 1.0).sum())
            if name == 'RANDU':
                ncen += int((inb & (9 * x - 6 * y + z == 2 * M)).sum())
    print(name, 'tuples', ntot, 'in ball exact', nin, 'float', nin_f, 'pi_full =', 6 * nin / ntot,
          'z(iid, N=2^29)', (6 * nin / ntot - math.pi) / (6 * math.sqrt((math.pi / 6) * (1 - math.pi / 6) / ntot)))
    if name == 'RANDU':
        w = np.array([1, 3, 5, 7, 9, 11, 12, 12, 12, 11, 9, 7, 5, 3, 1]) / 108
        print(' planes', planes.tolist())
        print(' dev from 2^29 p_j', np.round(planes - ntot * w, 1).tolist())
        print(' central plane in ball', ncen, 'fraction', ncen / nin)
        d = np.array([2 * k / math.sqrt(118) for k in range(-5, 6)])
        W = np.sum(1 - d ** 2)
        print(' W', W, '1/W', 1 / W, ' plane-sum pi', 6 * math.pi * (2 / math.sqrt(118)) * W / 8)
    else:
        print(' residues of (58,-375,-69).X mod 2^31 (units of 2^28):', sorted(r / 2 ** 28 for r in res))
