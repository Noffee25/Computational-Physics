import numpy as np, glob
from scipy import stats
def parse(fn):
    L=open(fn).read().split('\n'); d={}
    for i,l in enumerate(L):
        if l.startswith('BDAY'): d['BDAY']=np.array(L[i+1].split(),int)
        if l.startswith('RUNS'): d['RUNS']=np.array(L[i+1].split(),int)
    return d
lam=2.0
p=stats.poisson.pmf(np.arange(7),lam); p=np.append(p,1-p.sum())
tot=np.zeros(8)
for fn in sorted(glob.glob('out/val/kn_*.txt')):
    b=parse(fn)['BDAY']; tot+=b
    e=b.sum()*p; c=((b-e)**2/e).sum()
    print(fn, b, 'chi2 %.1f p %.3f'%(c, stats.chi2.sf(c,7)), 'mean J %.4f'%((b*np.arange(8)).sum()/b.sum()))
e=tot.sum()*p
print('pooled',tot, ((tot-e)**2/e).sum(), stats.chi2.sf(((tot-e)**2/e).sum(),7))
print('expected',np.round(e))
