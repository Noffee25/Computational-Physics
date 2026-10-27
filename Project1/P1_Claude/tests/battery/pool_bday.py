import sys; sys.path.insert(0, "/Users/noah/Downloads/ComputationalPhysicsPHYS25000/Project1/P1_Claude/tests/battery")
import numpy as np, glob
from scipy import stats
from knuth_table import parse, PR
for g,files in [('R250',['out/kn/kn_5_v0.txt']+sorted(glob.glob('out/val2/kn_5_*.txt'))),('MT',['out/kn/kn_7_v0.txt']+sorted(glob.glob('out/val2/kn_7_*.txt'))),('RANDU',['out/kn/kn_1_v0.txt']+sorted(glob.glob('out/val2/kn_1_*.txt')))]:
    tot=np.zeros(8)
    for f in files: tot+=parse(f)['BDAY']
    e=tot.sum()*PR['BDAY']; X=((tot-e)**2/e).sum()
    meanJ=(tot*np.arange(8)).sum()/tot.sum()
    print(g,'pooled R=%d'%tot.sum(),'obs',tot.astype(int),'\n   exp',np.round(e).astype(int),' X2=%.1f p=%.2e meanJ=%.4f (null %.4f)'%(X,stats.chi2.sf(X,7),meanJ,(PR['BDAY']*np.arange(8)).sum()))
