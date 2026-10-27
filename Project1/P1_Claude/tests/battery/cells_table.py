import numpy as np, glob, re
from scipy import stats
names={1:'RANDU',2:'RANDM',3:'RAND',4:'MINSTD',5:'R250',6:'RANNUM',7:'MT19937',8:'MIDSQ',9:'LOGISTIC'}
def load(fn):
    L=open(fn).read().split('\n')
    Nt=int(re.search(r'Nt\s+(\d+)',L[0]).group(1))
    rows=[]
    for l in L[2:]:
        if l.strip():
            k,C,s2,ne,mx=map(int,l.split()); rows.append((k,C,s2,ne,mx))
    return Nt,rows
for suf in ['d2','d3','d4','lag']:
    print('=====',suf)
    for g in range(1,10):
        Nt,rows=load(f'out/cells/c_{g}_{suf}.txt')
        out=[]
        for k,C,s2,ne,mx in rows:
            lam=Nt/C; X2=s2/lam-Nt; z=(X2-(C-1))/np.sqrt(2*(C-1))
            emp=ne/C; empexp=np.exp(-lam)
            out.append('k=%d z=%.1f emp=%.3f(%.3f)'%(k,z,emp,empexp))
        print('%-9s'%names[g],' | '.join(out))
